"""
furnace_model.py
================
Physics model of the Verzinnungsanlage (tin-coating line) tin bath.

Provides the individual thermal sub-models that make up the energy balance:

    C_B * dT/dt = P_heating(t-τ)
                  - P_losses(T)
                  - P_band(T, strip_params)
                  - P_airknife(state_%)
                  - P_water_cooling

All functions are stateless (except TransportDelay). The simulation runner
in simulation.py wires them together and steps through time.
"""

import numpy as np
import math
from collections import deque
import parameters as p


# ---------------------------------------------------------------------------
# 1-D lookup with linear interpolation (clamped at boundaries)
# ---------------------------------------------------------------------------
def lookup1d(x_bp, y_tbl, x: float) -> float:
    """Linear interpolation on a breakpoint table, clamped at boundaries."""
    return float(np.interp(x, x_bp, y_tbl))


# ---------------------------------------------------------------------------
# Air-knife operating state  (mirrors Simulink MATLAB function)
# ---------------------------------------------------------------------------
def airknife_percent(druck1: float, druck2: float) -> float:
    """Convert two nozzle pressures [mbar] → operating state [%]."""
    mean = (druck1 + druck2) / 2.0
    if mean < p.airknife_pressure_min:
        return 0.0
    return 100.0 * min(mean, float(p.airknife_pressure_max)) / p.airknife_pressure_max


# ---------------------------------------------------------------------------
# Strip heat extraction  (Wärmestrom durch Band)
# ---------------------------------------------------------------------------
def strip_mass_flow(speed_m_min: float, width_mm: float,
                    thickness_mm: float, rho_strip: float) -> float:
    """Strip mass flow rate (Massenstrom) [kg/s] = ρ · v · width · thickness."""
    v = speed_m_min * (1.0 / 60.0)   # m/min → m/s
    w = width_mm    * 1e-3           # mm → m
    t = thickness_mm * 1e-3          # mm → m
    return rho_strip * v * w * t     # kg/s


def strip_heat_flow(T_bath: float, T_strip_in: float,
                    speed_m_min: float, width_mm: float,
                    thickness_mm: float, cp_strip: float,
                    rho_strip: float) -> float:
    """Heat extracted by the steel strip passing through the tin bath [W]."""
    mdot = strip_mass_flow(speed_m_min, width_mm, thickness_mm, rho_strip)
    dT   = max(0.0, T_bath - T_strip_in)
    return mdot * cp_strip * dT          # W


# ---------------------------------------------------------------------------
# Heating correction  (mirrors korrHeizleistung Simulink MATLAB function)
# ---------------------------------------------------------------------------
def korr_heizleistung(fraction: float, strip_running: bool,
                      strip_length_m: float) -> float:
    """Limit heater to 20 % during standstill or when strip is too short."""
    if (not strip_running) or (strip_length_m < p.strip_length_min):
        return min(fraction, p.heating_limit_standstill)
    return fraction


# ---------------------------------------------------------------------------
# Setpoint selection  (mirrors select_setpoint Simulink MATLAB function)
# ---------------------------------------------------------------------------
def select_setpoint(strip_running: bool,
                    T_soll_betrieb: float,
                    T_soll_stillstand: float) -> float:
    """Return the active temperature setpoint based on strip running state."""
    return T_soll_betrieb if strip_running else T_soll_stillstand


# ---------------------------------------------------------------------------
# Transport-delay buffer  (mirrors Simulink TransportDelay block)
# ---------------------------------------------------------------------------
class TransportDelay:
    """
    FIFO transport delay implemented as a circular deque.

    All values entering the buffer at time t are emitted unchanged after
    exactly `delay_s` seconds — matching the Simulink TransportDelay block.
    """

    def __init__(self, delay_s: float, dt: float, init: float = 0.0):
        steps = max(1, int(round(delay_s / dt)))
        self._buf = deque([init] * steps, maxlen=steps)

    def step(self, value: float) -> float:
        """Push new value in; return the delayed (oldest) value."""
        out = self._buf[0]
        self._buf.append(value)
        return out

    def fill(self, value: float):
        """Re-initialise the entire buffer with a constant value."""
        for k in range(len(self._buf)):
            self._buf[k] = value


# ---------------------------------------------------------------------------
# PT1 first-order lag  (Verzögerungsglied 1. Ordnung)
# ---------------------------------------------------------------------------
class PT1:
    """
    First-order lag element:  τ·dy/dt = u − y   ⇔   Y(s)/U(s) = 1/(1 + τ·s)

    Unlike a pure transport (dead-time) delay, which holds its initial output
    and then jumps to the buffered value after τ seconds, a PT1 lets the output
    approach the input *smoothly* (exponentially), with no abrupt step.

    Discretised with the exact zero-order-hold solution (unconditionally
    stable for any dt):
        α      = exp(-dt / τ)
        y[k]   = α · y[k-1] + (1 − α) · u[k]

    The time constant τ can be changed between steps (e.g. Betrieb vs.
    Stillstand) via the `tau` attribute.
    """

    def __init__(self, tau_s: float, dt: float, init: float = 0.0):
        self.tau = tau_s
        self.dt  = dt
        self.y   = init

    def step(self, u: float) -> float:
        if self.tau <= 0.0:
            self.y = u                      # τ=0 → no lag, pass through
        else:
            alpha  = math.exp(-self.dt / self.tau)
            self.y = alpha * self.y + (1.0 - alpha) * u
        return self.y

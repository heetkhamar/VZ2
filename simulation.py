"""
simulation.py
=============
Physics-based thermal simulation of the Verzinnungsanlage (tin-coating line).

The model captures the energy balance of the tin bath (Zinnbad):

    C_B * dT_bath/dt = P_heating(t - τ_heat)
                       - P_losses(T_bath)          [radiation + convection]
                       - P_band(T_bath, band_params, t - τ_band)
                       - P_airknife(state_%)
                       - P_water_cooling(t - τ_water)

All parameters are defined in parameters.py and imported here.

Usage
-----
    python simulation.py

To change simulation time, heating power, or strip parameters edit the
section marked "SIMULATION INPUTS" at the bottom of this file, or import
run_simulation() from another script.
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import deque

import parameters as p


# ---------------------------------------------------------------------------
# Helper: 1-D lookup table with linear interpolation + extrapolation clamp
# ---------------------------------------------------------------------------
def lookup1d(x_bp: list, y_tbl: list, x: float) -> float:
    return float(np.interp(x, x_bp, y_tbl))


# ---------------------------------------------------------------------------
# Air-knife operating state calculation
# ---------------------------------------------------------------------------
def airknife_percent(druck1: float, druck2: float) -> float:
    """Convert two nozzle pressures [mbar] to operating state [%]."""
    mean = (druck1 + druck2) / 2.0
    if mean < p.airknife_pressure_min:
        mean = 0.0
    elif mean > p.airknife_pressure_max:
        mean = float(p.airknife_pressure_max)
    if mean == 0.0:
        return 0.0
    return 100.0 * mean / p.airknife_pressure_max


# ---------------------------------------------------------------------------
# Strip heat flow (Wärmestrom durch Band)
# ---------------------------------------------------------------------------
def strip_heat_flow(T_bath: float,
                    T_strip_in: float,
                    strip_speed_m_min: float,
                    strip_width_mm: float,
                    strip_thickness_mm: float,
                    cp_strip: float,
                    rho_strip: float) -> float:
    """
    Return the heat extracted from the bath by the moving steel strip [W].

    mass_flow = rho * v [m/s] * w [m] * t [m]   [kg/s]
    Q_dot     = mass_flow * cp * (T_bath - T_strip_in)
    """
    v_m_s  = strip_speed_m_min * (1.0 / 60.0)      # m/min → m/s
    w_m    = strip_width_mm    * 1e-3               # mm → m
    t_m    = strip_thickness_mm * 1e-3              # mm → m

    mass_flow = rho_strip * v_m_s * w_m * t_m      # kg/s

    delta_T = T_bath - T_strip_in
    if delta_T < 0.0:
        delta_T = 0.0                               # strip can only cool the bath

    return mass_flow * cp_strip * delta_T           # W


# ---------------------------------------------------------------------------
# Transport-delay buffer
# ---------------------------------------------------------------------------
class TransportDelay:
    """Fixed-time transport delay implemented as a circular buffer."""

    def __init__(self, delay_s: float, dt: float, initial_value: float = 0.0):
        steps = max(1, int(round(delay_s / dt)))
        self._buf = deque([initial_value] * steps, maxlen=steps)

    def step(self, value: float) -> float:
        """Push new value, return delayed value."""
        delayed = self._buf[0]
        self._buf.append(value)
        return delayed


# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------
def run_simulation(
        t_start: float       = None,
        t_end: float         = None,
        dt: float            = None,
        heating_fraction: float = None,
        strip_running: bool  = True,
        strip_speed:   float = None,
        strip_width:   float = None,
        strip_thickness: float = None,
        T_strip_in: float    = None,
        airknife_druck1: float = None,
        airknife_druck2: float = None,
        water_cooling_power: float = None,
        T_bath_init: float   = None,
):
    """
    Run the tin-bath thermal simulation.

    Parameters
    ----------
    t_start, t_end, dt          : simulation window [s]
    heating_fraction            : inductive heater demand [0..1]
    strip_running               : whether steel strip is passing through the bath
    strip_speed                 : strip speed [m/min]
    strip_width                 : strip width [mm]
    strip_thickness             : strip thickness [mm]
    T_strip_in                  : strip entry temperature [°C]
    airknife_druck1/2           : air-knife pressures OS / US [mbar]
    water_cooling_power         : constant water cooling power [W]
    T_bath_init                 : initial bath temperature [°C]

    Returns
    -------
    dict with keys:
        't'           – time array [s]
        'T_bath'      – bath temperature [°C]
        'P_heating'   – inductive heating power applied [W]
        'P_losses'    – radiation + convection losses [W]
        'P_band'      – strip heat extraction [W]
        'P_airknife'  – air-knife cooling [W]
        'P_water'     – water cooling [W]
        'P_net'       – net power into bath [W]  (= dT/dt * C_B)
    """
    # --- resolve parameters (fall back to parameters.py defaults) -----------
    t0   = t_start          if t_start          is not None else p.t_start
    t1   = t_end            if t_end            is not None else p.t_end
    step = dt               if dt               is not None else p.dt
    h_fr = heating_fraction if heating_fraction is not None else p.heating_fraction_default
    v    = strip_speed      if strip_speed      is not None else p.strip_speed
    bw   = strip_width      if strip_width      is not None else p.strip_width
    bt   = strip_thickness  if strip_thickness  is not None else p.strip_thickness
    T_in = T_strip_in       if T_strip_in       is not None else p.T_strip_in
    dk1  = airknife_druck1  if airknife_druck1  is not None else p.airknife_druck1_default
    dk2  = airknife_druck2  if airknife_druck2  is not None else p.airknife_druck2_default
    P_w  = water_cooling_power if water_cooling_power is not None else p.water_cooling_power_default
    T0   = T_bath_init      if T_bath_init      is not None else p.T_bath_init

    # --- clamp heater fraction -----------------------------------------------
    h_fr = max(0.0, min(1.0, h_fr))
    P_heat_demand = h_fr * p.P_heating_max      # W

    # --- transport delays ----------------------------------------------------
    delay_heat  = TransportDelay(p.delay_heating,    step, P_heat_demand)
    delay_band  = TransportDelay(p.delay_band,       step, 0.0)
    delay_water = TransportDelay(p.delay_water_cool, step, P_w)

    # --- air-knife cooling (constant during open-loop run) -------------------
    ak_state  = airknife_percent(dk1, dk2)
    P_ak      = lookup1d(p.airknife_state_bp, p.airknife_power_tbl, ak_state)

    # --- allocate result arrays ----------------------------------------------
    n_steps = int(round((t1 - t0) / step)) + 1
    t_arr   = np.linspace(t0, t0 + (n_steps - 1) * step, n_steps)

    T_arr       = np.zeros(n_steps)
    P_heat_arr  = np.zeros(n_steps)
    P_loss_arr  = np.zeros(n_steps)
    P_band_arr  = np.zeros(n_steps)
    P_ak_arr    = np.zeros(n_steps)
    P_water_arr = np.zeros(n_steps)
    P_net_arr   = np.zeros(n_steps)

    # --- initial state -------------------------------------------------------
    T = T0

    for i in range(n_steps):
        # Delayed heating power reaching the bath
        P_h = delay_heat.step(P_heat_demand)

        # Radiation + convection losses
        P_loss = lookup1d(p.loss_temp_bp, p.loss_power_tbl, T)

        # Strip heat extraction (with transport delay)
        if strip_running:
            Q_band_raw = strip_heat_flow(T, T_in, v, bw, bt, p.cp_strip, p.rho_strip)
        else:
            Q_band_raw = 0.0
        P_band = delay_band.step(Q_band_raw)

        # Delayed water cooling
        P_wc = delay_water.step(P_w)

        # Net power balance
        P_net = P_h - P_loss - P_band - P_ak - P_wc

        # Forward Euler integration: C_B * dT/dt = P_net
        dT = (P_net / p.C_B) * step
        T  = T + dT

        # Store
        T_arr[i]       = T
        P_heat_arr[i]  = P_h
        P_loss_arr[i]  = P_loss
        P_band_arr[i]  = P_band
        P_ak_arr[i]    = P_ak
        P_water_arr[i] = P_wc
        P_net_arr[i]   = P_net

    return {
        't':         t_arr,
        'T_bath':    T_arr,
        'P_heating': P_heat_arr,
        'P_losses':  P_loss_arr,
        'P_band':    P_band_arr,
        'P_airknife':P_ak_arr,
        'P_water':   P_water_arr,
        'P_net':     P_net_arr,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_results(results: dict):
    t_min = results['t'] / 60.0   # convert to minutes for readability

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # --- Temperature plot ---
    ax1 = axes[0]
    ax1.plot(t_min, results['T_bath'], color='tab:red', linewidth=1.8,
             label='Zinnbad-Temperatur')
    ax1.set_ylabel('Temperatur [°C]')
    ax1.set_title('Verzinnungsanlage – Thermische Simulation')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.4)
    ax1.set_ylim(bottom=200)

    # --- Power plot ---
    ax2 = axes[1]
    ax2.plot(t_min, results['P_heating'] / 1000, label='Heizleistung (ind.)',
             color='tab:orange', linewidth=1.5)
    ax2.plot(t_min, results['P_losses']  / 1000, label='Strahlungs-/Konvektionsverluste',
             color='tab:blue',   linewidth=1.2, linestyle='--')
    ax2.plot(t_min, results['P_band']    / 1000, label='Wärmestrom Band',
             color='tab:green',  linewidth=1.2, linestyle='-.')
    ax2.plot(t_min, results['P_airknife']/ 1000, label='Air-knife Kühlung',
             color='tab:purple', linewidth=1.2, linestyle=':')
    ax2.plot(t_min, results['P_water']   / 1000, label='Wasserkühlung',
             color='tab:cyan',   linewidth=1.2, linestyle=':')
    ax2.plot(t_min, results['P_net']     / 1000, label='Netto-Leistung (dT/dt)',
             color='black',      linewidth=1.0, linestyle='-', alpha=0.5)
    ax2.set_ylabel('Leistung [kW]')
    ax2.set_xlabel('Zeit [min]')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig('simulation_results.png', dpi=150)
    print("Plot saved to simulation_results.png")
    plt.show()


# ---------------------------------------------------------------------------
# SIMULATION INPUTS  –  edit these to configure an open-loop run
# ---------------------------------------------------------------------------
if __name__ == '__main__':

    # ---- Time ---------------------------------------------------------------
    SIM_DURATION_HOURS = 1          # hours
    DT_SECONDS         = 0.1        # integration step [s]

    # ---- Inductive heating --------------------------------------------------
    # Fraction of installed max power (P_heating_max = 7000 W)
    # 1.0 = full power, 0.0 = off
    HEATING_FRACTION = 0.75

    # ---- Steel strip --------------------------------------------------------
    STRIP_RUNNING   = True
    STRIP_SPEED     = 70.0    # m/min
    STRIP_WIDTH     = 350.0   # mm
    STRIP_THICKNESS = 0.340     # mm
    T_STRIP_ENTRY   = 30.0   # °C – strip temperature entering the bath

    # ---- Air knife (Luftmesser) ----------------------------------------------
    # Set both pressures to 0 to disable air-knife cooling
    AIRKNIFE_DRUCK1 = 300.0   # mbar  (OS – Oberseite)
    AIRKNIFE_DRUCK2 = 300.0   # mbar  (US – Unterseite)

    # ---- Water cooling -------------------------------------------------------
    WATER_COOLING = 10000.0       # W  (0 = off)

    # ---- Initial bath temperature --------------------------------------------
    T_BATH_START = 270.0      # °C

    # -------------------------------------------------------------------------
    print("Starting Verzinnungsanlage simulation …")
    print(f"  Duration       : {SIM_DURATION_HOURS} h")
    print(f"  Heating        : {HEATING_FRACTION * 100:.0f} % of {p.P_heating_max:.0f} W max"
          f" = {HEATING_FRACTION * p.P_heating_max:.0f} W")
    print(f"  Strip running  : {STRIP_RUNNING}")
    if STRIP_RUNNING:
        print(f"    Speed        : {STRIP_SPEED} m/min")
        print(f"    Width        : {STRIP_WIDTH} mm")
        print(f"    Thickness    : {STRIP_THICKNESS} mm")
        print(f"    Entry temp   : {T_STRIP_ENTRY} °C")
    ak_pct = airknife_percent(AIRKNIFE_DRUCK1, AIRKNIFE_DRUCK2)
    print(f"  Air knife      : {ak_pct:.1f} % operating state")
    print(f"  Water cooling  : {WATER_COOLING:.0f} W")
    print(f"  T_bath(0)      : {T_BATH_START} °C")
    print()

    results = run_simulation(
        t_start          = 0.0,
        t_end            = SIM_DURATION_HOURS * 3600.0,
        dt               = DT_SECONDS,
        heating_fraction = HEATING_FRACTION,
        strip_running    = STRIP_RUNNING,
        strip_speed      = STRIP_SPEED,
        strip_width      = STRIP_WIDTH,
        strip_thickness  = STRIP_THICKNESS,
        T_strip_in       = T_STRIP_ENTRY,
        airknife_druck1  = AIRKNIFE_DRUCK1,
        airknife_druck2  = AIRKNIFE_DRUCK2,
        water_cooling_power = WATER_COOLING,
        T_bath_init      = T_BATH_START,
    )

    T_final = results['T_bath'][-1]
    print(f"Simulation complete.")
    print(f"  Final bath temperature : {T_final:.2f} °C")
    print(f"  Mean heating power     : {results['P_heating'].mean():.1f} W")
    print(f"  Mean losses            : {results['P_losses'].mean():.1f} W")
    print(f"  Mean band heat flow    : {results['P_band'].mean():.1f} W")
    print()

    plot_results(results)

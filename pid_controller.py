"""
pid_controller.py
=================
Gain-scheduled PID controller for the Verzinnungsanlage tin-bath temperature.

Gains (Kp, Ti, Td) are looked up as a function of strip speed using the
1-D breakpoint tables extracted from the Simulink model (system_1349.xml)
and stored in parameters.py.

Control structure (simple parallel PID form):
    u = Kp * e + Ki * ∫e dt - Kd * (dT_meas/dt)   (clamped to [u_min, u_max])

Anti-windup  : NONE – the integrator accumulates freely, even when the
               output is saturated (a plain, simple controller).
Derivative   : on measurement only — avoids setpoint-kick on step changes.
Output clamp : [pid_output_min, pid_output_max] from parameters.py

Imported by simulation.py; gain tables live in parameters.py.
"""

import numpy as np
import parameters as p
from furnace_model import lookup1d


class PIDController:
    """
    Simple discrete PID controller in parallel form (NO anti-windup), with two
    gain sources:

      * SCHEDULED mode (default) – Kp is looked up from the Massenstrom table
        (pid_Kp_table over pid_massflow_bp); Ki, Kd are the constants pid_Ki/pid_Kd.
      * MANUAL mode (manual=True) – fixed ZEITBEREICH (time-domain) parameters
        Kp, Ti, Td, bypassing the table. Converted internally to the parallel
        gains Ki = Kp/Ti and Kd = Kp·Td.

    Parameters
    ----------
    dt     : simulation time step [s]
    u_min  : lower output clamp (default 0 – heater can't go negative)
    u_max  : upper output clamp (default 100 – output is in PERCENT, matching
             the Simulink output scale; the caller converts %/100 → fraction)
    manual : if True, use the fixed Zeitbereich params instead of the table
    Kp     : proportional gain (manual mode). If None → p.pid_manual_Kp
    Ti     : integration time [s] (manual mode). If None → p.pid_manual_Ti
    Td     : derivative time  [s] (manual mode). If None → p.pid_manual_Td
             (these are ignored unless manual=True)
    """

    def __init__(self, dt: float,
                 u_min: float = p.pid_output_min,
                 u_max: float = p.pid_output_max,
                 manual: bool = False,
                 Kp: float = None, Ti: float = None, Td: float = None):
        self.dt    = dt
        self.u_min = u_min
        self.u_max = u_max

        # Manual-mode configuration (Zeitbereich: Kp, Ti, Td)
        self.manual    = manual
        self.Kp_manual = Kp if Kp is not None else p.pid_manual_Kp
        self.Ti_manual = Ti if Ti is not None else p.pid_manual_Ti
        self.Td_manual = Td if Td is not None else p.pid_manual_Td

        self._integral  = 0.0
        self._prev_meas = None

    def reset(self, initial_output: float = 0.0,
              measurement: float = None):
        """
        Pre-load the integrator so the first output equals initial_output.
        Call once before the time loop to avoid a start-up transient.
        """
        self._integral  = initial_output
        self._prev_meas = measurement

    def step(self, setpoint: float, measurement: float,
             mass_flow: float) -> float:
        """
        Compute one control step.

        Parameters
        ----------
        setpoint    : temperature setpoint [°C]
        measurement : current bath temperature [°C]
        mass_flow   : strip mass flow (Massenstrom) for Kp scheduling [kg/s]

        Returns
        -------
        u : control output in [u_min, u_max]
        """
        # Parallel-form gains – source depends on mode:
        if self.manual:
            # MANUAL – Zeitbereich params (Kp, Ti, Td) converted to parallel gains:
            #   Ki = Kp / Ti      Kd = Kp · Td
            Kp = self.Kp_manual
            Ti = self.Ti_manual
            Td = self.Td_manual
            Ki = Kp / Ti if Ti > 0.0 else 0.0   # Ti<=0 → integral off (safe guard)
            Kd = Kp * Td
            self.Kp, self.Ti, self.Td = Kp, Ti, Td
        else:
            # SCHEDULED – Kp from the Massenstrom table; Ki, Kd constant.
            Kp = lookup1d(p.pid_massflow_bp, p.pid_Kp_table, mass_flow)
            Ki = p.pid_Ki
            Kd = p.pid_Kd
            # Derived time-domain (Zeitbereich) values, for reference/logging:
            #   Ti = Kp / Ki   (integration time [s])
            #   Td = Kd / Kp   (derivative time  [s])
            self.Kp = Kp
            self.Ti = Kp / Ki if Ki != 0.0 else float('inf')
            self.Td = Kd / Kp if Kp != 0.0 else 0.0

        e = setpoint - measurement

        # Derivative on measurement (avoids setpoint-kick). On the first call
        # there is no previous sample, so the derivative is 0.
        if self._prev_meas is None:
            self._prev_meas = measurement
        d_meas = (measurement - self._prev_meas) / self.dt

        # Simple controller: integrate every step (NO anti-windup).
        # The integral accumulates freely even when the output is saturated.
        self._integral += Ki * e * self.dt

        # P + I + D output, then clamp to the actuator limits [u_min, u_max].
        # Derivative on measurement → term is −Kd·d(measurement)/dt.
        u = Kp * e + self._integral - Kd * d_meas
        u = max(self.u_min, min(self.u_max, u))

        self._prev_meas = measurement
        return u

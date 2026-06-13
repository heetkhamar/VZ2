"""
simulation.py
=============
Physics-based thermal simulation of the Verzinnungsanlage (tin-coating line).

Energy balance of the tin bath (Zinnbad):

    C_B * dT_bath/dt = P_heating(t - τ)
                       - P_losses(T_bath)
                       - P_band(T_bath, strip_params)
                       - P_airknife(state_%)
                       - P_water_cooling

Two modes
---------
1. Manual / open-loop  → run_simulation()
   Heating driven by a time-varying schedule array (defined in parameters.py).

2. CSV + PID           → run_simulation_with_data()
   All process variables come from a CSV file (see data_loader.py).
   A gain-scheduled PID controller drives the inductive heater.
   Simulated temperature is plotted against the actual measured temperature
   from the CSV for direct comparison.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import deque

import parameters as p
import data_loader as dl


# ---------------------------------------------------------------------------
# Helper: 1-D lookup with linear interpolation (clamped at boundaries)
# ---------------------------------------------------------------------------
def lookup1d(x_bp, y_tbl, x: float) -> float:
    return float(np.interp(x, x_bp, y_tbl))


# ---------------------------------------------------------------------------
# Air-knife operating state  (from Simulink MATLAB function)
# ---------------------------------------------------------------------------
def airknife_percent(druck1: float, druck2: float) -> float:
    """Convert two nozzle pressures [mbar] → operating state [%]."""
    mean = (druck1 + druck2) / 2.0
    if mean < p.airknife_pressure_min:
        return 0.0
    mean = min(mean, float(p.airknife_pressure_max))
    return 100.0 * mean / p.airknife_pressure_max


# ---------------------------------------------------------------------------
# Strip heat flow  (Wärmestrom durch Band)
# ---------------------------------------------------------------------------
def strip_heat_flow(T_bath, T_strip_in, speed_m_min, width_mm,
                    thickness_mm, cp_strip, rho_strip) -> float:
    v   = speed_m_min * (1.0 / 60.0)   # m/min → m/s
    w   = width_mm    * 1e-3            # mm → m
    t   = thickness_mm * 1e-3           # mm → m
    mdot = rho_strip * v * w * t        # kg/s
    dT   = max(0.0, T_bath - T_strip_in)
    return mdot * cp_strip * dT         # W


# ---------------------------------------------------------------------------
# Heating correction  (korrHeizleistung – Simulink MATLAB function)
# Limits heater to 20 % when strip is not running or strip is too short.
# ---------------------------------------------------------------------------
def korr_heizleistung(fraction: float, strip_running: bool,
                      strip_length_m: float) -> float:
    if (not strip_running) or (strip_length_m < p.strip_length_min):
        return min(fraction, p.heating_limit_standstill)
    return fraction


# ---------------------------------------------------------------------------
# Setpoint selection  (Simulink MATLAB function)
# ---------------------------------------------------------------------------
def select_setpoint(strip_running: bool, T_soll_betrieb: float,
                    T_soll_stillstand: float) -> float:
    return T_soll_betrieb if strip_running else T_soll_stillstand


# ---------------------------------------------------------------------------
# Transport-delay buffer (circular deque)
# ---------------------------------------------------------------------------
class TransportDelay:
    def __init__(self, delay_s: float, dt: float, init: float = 0.0):
        steps = max(1, int(round(delay_s / dt)))
        self._buf = deque([init] * steps, maxlen=steps)

    def step(self, value: float) -> float:
        out = self._buf[0]
        self._buf.append(value)
        return out

    def fill(self, value: float):
        """Re-initialise entire buffer with a constant value."""
        for k in range(len(self._buf)):
            self._buf[k] = value


# ---------------------------------------------------------------------------
# PID controller – gain-scheduled, anti-windup, derivative on measurement
# ---------------------------------------------------------------------------
class PIDController:
    """
    Discrete PID in parallel form.

    Gains are scheduled by strip speed from the Simulink lookup tables
    (system_1349.xml):
        Kp = f(v),  Ti = f(v),  Td = f(v)
        Ki = Kp / Ti,  Kd = Kp * Td

    Anti-windup via integral clamping (back-calculation).
    Derivative is taken on the measurement to avoid set-point kick.
    """

    def __init__(self, dt: float, u_min: float = 0.0, u_max: float = 1.0):
        self.dt    = dt
        self.u_min = u_min
        self.u_max = u_max

        self._integral      = 0.0
        self._prev_meas     = None

    def reset(self, initial_output: float = 0.0, measurement: float = None):
        """Pre-load integral so the first output equals initial_output."""
        self._integral  = initial_output
        self._prev_meas = measurement

    def step(self, setpoint: float, measurement: float,
             strip_speed_m_min: float) -> float:
        """
        Compute one PID step.

        Returns control output clamped to [u_min, u_max].
        Gains are looked up from speed-scheduled tables in parameters.py.
        """
        v_m_s = strip_speed_m_min / 60.0

        Kp = lookup1d(p.pid_speed_bp, p.pid_Kp_table, v_m_s)
        Ti = lookup1d(p.pid_speed_bp, p.pid_Ti_table, v_m_s)
        Td = lookup1d(p.pid_speed_bp, p.pid_Td_table, v_m_s)
        Ki = Kp / max(Ti, 1.0)
        Kd = Kp * Td

        if self._prev_meas is None:
            self._prev_meas = measurement

        e       = setpoint - measurement
        d_meas  = (measurement - self._prev_meas) / self.dt   # derivative on meas

        # Unsaturated output
        u_raw = Kp * e + self._integral - Kd * d_meas

        # Clamp
        u = max(self.u_min, min(self.u_max, u_raw))

        # Integrate with anti-windup: only accumulate while not saturated
        if u_raw == u:
            self._integral += Ki * e * self.dt
        else:
            # Back-calculation anti-windup
            self._integral += Ki * e * self.dt - (u_raw - u)

        self._prev_meas = measurement
        return u


# ---------------------------------------------------------------------------
# MODE 1 – Manual / open-loop  (heating_schedule array, no CSV)
# ---------------------------------------------------------------------------
def run_simulation(
        t_start=None, t_end=None, dt=None,
        heating_schedule=None,
        strip_running: bool = True,
        strip_speed=None, strip_width=None, strip_thickness=None,
        T_strip_in=None,
        airknife_druck1=None, airknife_druck2=None,
        water_cooling_power=None,
        T_bath_init=None,
):
    """
    Open-loop simulation driven by a heating-fraction schedule array.
    All other inputs are constant scalars taken from parameters.py.
    """
    t0   = t_start or p.t_start
    t1   = t_end   or p.t_end
    step = dt      or p.dt
    v    = strip_speed     or p.strip_speed
    bw   = strip_width     or p.strip_width
    bt   = strip_thickness or p.strip_thickness
    T_in = T_strip_in      if T_strip_in      is not None else p.T_strip_in
    dk1  = airknife_druck1 if airknife_druck1 is not None else p.airknife_druck1_default
    dk2  = airknife_druck2 if airknife_druck2 is not None else p.airknife_druck2_default
    P_w  = water_cooling_power if water_cooling_power is not None else p.water_cooling_power_default
    T0   = T_bath_init     if T_bath_init     is not None else p.T_bath_init

    n_steps = int(round((t1 - t0) / step)) + 1

    if heating_schedule is not None:
        sched = np.clip(np.asarray(heating_schedule, float), 0, 1)
        if len(sched) != n_steps:
            raise ValueError(f"heating_schedule length {len(sched)} ≠ n_steps {n_steps}")
    else:
        sched = (p.heating_schedule if len(p.heating_schedule) == n_steps
                 else p.generate_heating_schedule(n_steps))

    t_arr       = np.linspace(t0, t0 + (n_steps - 1) * step, n_steps)
    T_arr       = np.zeros(n_steps)
    frac_arr    = np.zeros(n_steps)
    P_heat_arr  = np.zeros(n_steps)
    P_loss_arr  = np.zeros(n_steps)
    P_band_arr  = np.zeros(n_steps)
    P_ak_arr    = np.zeros(n_steps)
    P_water_arr = np.zeros(n_steps)
    P_net_arr   = np.zeros(n_steps)

    ak_state = airknife_percent(dk1, dk2)
    P_ak     = lookup1d(p.airknife_state_bp, p.airknife_power_tbl, ak_state)

    delay_heat  = TransportDelay(p.delay_heating,    step, sched[0] * p.P_heating_max)
    delay_band  = TransportDelay(p.delay_band,       step, 0.0)
    delay_water = TransportDelay(p.delay_water_cool, step, P_w)

    T = T0

    for i in range(n_steps):
        P_demand = sched[i] * p.P_heating_max
        P_h      = delay_heat.step(P_demand)
        P_loss   = lookup1d(p.loss_temp_bp, p.loss_power_tbl, T)
        Q_raw    = strip_heat_flow(T, T_in, v, bw, bt, p.cp_strip, p.rho_strip) if strip_running else 0.0
        P_band   = delay_band.step(Q_raw)
        P_wc     = delay_water.step(P_w)
        P_net    = P_h - P_loss - P_band - P_ak - P_wc

        T += (P_net / p.C_B) * step

        T_arr[i]      = T;       frac_arr[i]    = sched[i]
        P_heat_arr[i] = P_h;    P_loss_arr[i]  = P_loss
        P_band_arr[i] = P_band; P_ak_arr[i]    = P_ak
        P_water_arr[i]= P_wc;   P_net_arr[i]   = P_net

    return dict(t=t_arr, T_bath=T_arr, heating_frac=frac_arr,
                P_heating=P_heat_arr, P_losses=P_loss_arr, P_band=P_band_arr,
                P_airknife=P_ak_arr, P_water=P_water_arr, P_net=P_net_arr)


# ---------------------------------------------------------------------------
# MODE 2 – CSV + PID  (full model, mirrors Simulink with-data system)
# ---------------------------------------------------------------------------
def run_simulation_with_data(
        csv_path: str,
        dt: float = 1.0,
        T_bath_init: float = None,
        mode: str = 'pid',           # 'pid' or 'open_loop'
):
    """
    Full simulation driven by ibA CSV data with a gain-scheduled PID controller.

    Parameters
    ----------
    csv_path    : path to the ibA CSV file (see data_loader.py for format)
    dt          : simulation time step [s]
    T_bath_init : initial bath temperature [°C]; if None, uses CSV first value
    mode        : 'pid'       – PID controller drives the heater (default)
                  'open_loop' – oven power taken directly from CSV Ofenleistung

    Returns
    -------
    dict with all simulated time-series plus 'actual_temp' and 'setpoint_temp'
    from the CSV (for comparison plotting).
    """
    # --- Load and interpolate CSV -------------------------------------------
    df     = dl.load_csv(csv_path)
    t_end  = df['t_s'].iloc[-1]
    n_steps = int(round(t_end / dt)) + 1
    t_grid  = np.linspace(0.0, (n_steps - 1) * dt, n_steps)
    data    = dl.interpolate_to_grid(df, t_grid)

    T0 = T_bath_init if T_bath_init is not None else float(data['actual_temp'][0])

    # --- Allocate result arrays ---------------------------------------------
    T_arr       = np.zeros(n_steps)
    frac_arr    = np.zeros(n_steps)
    P_heat_arr  = np.zeros(n_steps)
    P_loss_arr  = np.zeros(n_steps)
    P_band_arr  = np.zeros(n_steps)
    P_ak_arr    = np.zeros(n_steps)
    P_water_arr = np.zeros(n_steps)
    P_net_arr   = np.zeros(n_steps)

    # --- Initialise state ---------------------------------------------------
    T = T0

    init_frac  = data['oven_power_pct'][0] / 100.0
    init_P     = init_frac * p.P_heating_max

    # Two delay buffers – Simulink uses 94 s (Betrieb) / 150 s (Stillstand)
    delay_betrieb    = TransportDelay(p.delay_heating,            dt, init_P)
    delay_stillstand = TransportDelay(p.delay_heating_stillstand, dt, init_P)
    delay_band       = TransportDelay(p.delay_band,               dt, 0.0)

    # PID (pre-loaded so first output ≈ initial oven power)
    pid = PIDController(dt, u_min=p.pid_output_min, u_max=p.pid_output_max)
    pid.reset(initial_output=init_frac, measurement=T0)

    prev_running = bool(data['strip_running'][0] > 0.5)

    # --- Time loop ----------------------------------------------------------
    for i in range(n_steps):
        # Current CSV inputs
        v        = data['strip_speed'][i]
        bw       = data['strip_width'][i]
        bt       = data['strip_thickness'][i]
        T_in     = data['temp_under_roller'][i]
        ak_os    = data['airknife_os_mbar'][i]
        ak_us    = data['airknife_us_mbar'][i]
        P_w      = data['water_cooling_W'][i]
        running  = data['strip_running'][i] > 0.5
        s_len    = data['strip_length_in_oven'][i]
        T_soll_b = data['setpoint_temp'][i]
        T_soll_s = data['setpoint_standstill'][i]

        # Setpoint selection (Simulink MATLAB function)
        T_soll = select_setpoint(running, T_soll_b, T_soll_s)

        # Heater demand ---------------------------------------------------
        if mode == 'pid':
            raw_frac = pid.step(T_soll, T, strip_speed_m_min=v)
        else:
            raw_frac = data['oven_power_pct'][i] / 100.0

        # Heating correction (korrHeizleistung MATLAB function)
        frac = korr_heizleistung(raw_frac, running, s_len)

        P_demand = frac * p.P_heating_max

        # Transport delay: Betrieb 94 s, Stillstand 150 s
        # Both buffers are advanced every step; only the active one is used.
        P_h_betrieb    = delay_betrieb.step(P_demand)
        P_h_stillstand = delay_stillstand.step(P_demand)
        P_h = P_h_betrieb if running else P_h_stillstand

        # Strip heat extraction
        if running:
            Q_raw  = strip_heat_flow(T, T_in, v, bw, bt, p.cp_strip, p.rho_strip)
        else:
            Q_raw  = 0.0
        P_band = delay_band.step(Q_raw)

        # Air-knife cooling
        ak_pct = airknife_percent(ak_os, ak_us)
        P_ak   = lookup1d(p.airknife_state_bp, p.airknife_power_tbl, ak_pct)

        # Radiation + convection losses (lookup table)
        P_loss = lookup1d(p.loss_temp_bp, p.loss_power_tbl, T)

        # Net power → integrate bath temperature
        P_net = P_h - P_loss - P_band - P_ak - P_w
        T    += (P_net / p.C_B) * dt

        prev_running = running

        # Store
        T_arr[i]      = T;       frac_arr[i]    = frac
        P_heat_arr[i] = P_h;    P_loss_arr[i]  = P_loss
        P_band_arr[i] = P_band; P_ak_arr[i]    = P_ak
        P_water_arr[i]= P_w;    P_net_arr[i]   = P_net

    return dict(
        t            = t_grid,
        T_bath       = T_arr,
        heating_frac = frac_arr,
        P_heating      = P_heat_arr,
        P_losses       = P_loss_arr,
        P_band         = P_band_arr,
        P_airknife     = P_ak_arr,
        P_water        = P_water_arr,
        P_net          = P_net_arr,
        # From CSV (for comparison)
        actual_temp    = data['actual_temp'],
        setpoint_temp  = np.array([
            select_setpoint(data['strip_running'][i] > 0.5,
                            data['setpoint_temp'][i],
                            data['setpoint_standstill'][i])
            for i in range(n_steps)
        ]),
        strip_running  = data['strip_running'],
        oven_power_pct = data['oven_power_pct'],
    )


# ---------------------------------------------------------------------------
# Plotting – open-loop results
# ---------------------------------------------------------------------------
def plot_results(results: dict, save_path: str = 'simulation_results.png'):
    t_min = results['t'] / 60.0

    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

    ax1 = axes[0]
    ax1.plot(t_min, results['T_bath'], color='tab:red', lw=1.8,
             label='Zinnbad-Temperatur')
    ax1.set_ylabel('Temperatur [°C]')
    ax1.set_title('Verzinnungsanlage – Thermische Simulation (open-loop)')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.4)

    ax2 = axes[1]
    ax2.plot(t_min, results['P_heating'] / 1000, color='tab:orange', lw=1.5,
             label='Heizleistung (ind.)')
    ax2.plot(t_min, results['P_losses']  / 1000, color='tab:blue',   lw=1.2,
             ls='--', label='Strahlung/Konvektion')
    ax2.plot(t_min, results['P_band']    / 1000, color='tab:green',  lw=1.2,
             ls='-.', label='Wärmestrom Band')
    ax2.plot(t_min, results['P_airknife']/ 1000, color='tab:purple', lw=1.2,
             ls=':', label='Air-knife Kühlung')
    ax2.plot(t_min, results['P_water']   / 1000, color='tab:cyan',   lw=1.2,
             ls=':', label='Wasserkühlung')
    ax2.plot(t_min, results['P_net']     / 1000, color='black',      lw=1.0,
             alpha=0.5, label='Netto-Leistung')
    ax2.set_ylabel('Leistung [kW]')
    ax2.set_xlabel('Zeit [min]')
    ax2.grid(True, alpha=0.4)

    ax2r = ax2.twinx()
    ax2r.step(t_min, results['heating_frac'] * 100, color='tab:red',
              lw=1.0, alpha=0.6, where='post', label='Heizer-Sollwert [%]')
    ax2r.set_ylabel('Heizer-Sollwert [%]', color='tab:red')
    ax2r.tick_params(axis='y', labelcolor='tab:red')
    ax2r.set_ylim(0, 110)

    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2r.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, loc='upper center',
               bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=8, framealpha=0.9)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(save_path, dpi=150)
    plt.close()
    print(f"Plot saved to {save_path}")


# ---------------------------------------------------------------------------
# Plotting – CSV + PID comparison
# ---------------------------------------------------------------------------
def plot_results_comparison(results: dict,
                             save_path: str = 'simulation_results.png'):
    """
    Three-panel comparison plot:
      1. Temperatures: simulated vs actual, setpoint, vessel wall temps
      2. Heating power: PID output vs actual CSV oven power + all loss terms
      3. Strip operating state
    """
    t_min = results['t'] / 60.0

    fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True,
                              gridspec_kw={'height_ratios': [3, 3, 1]})

    # --- Panel 1: Temperatures ----------------------------------------------
    ax1 = axes[0]
    ax1.plot(t_min, results['actual_temp'],   color='tab:green', lw=1.5,
             ls='-',  label='Isttemperatur (Messung)')
    ax1.plot(t_min, results['T_bath'],       color='tab:red',   lw=1.8,
             ls='-',  label='Zinnbad-Temperatur (Simulation)')
    ax1.plot(t_min, results['setpoint_temp'],color='black',     lw=1.0,
             ls='--', alpha=0.7, label='Solltemperatur')
    ax1.set_ylabel('Temperatur [°C]')
    ax1.set_title('Verzinnungsanlage – Simulation vs. Messung (PID-Regelung)')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.4)

    # RMSE annotation
    rmse = float(np.sqrt(np.mean((results['T_bath'] - results['actual_temp'])**2)))
    ax1.text(0.01, 0.04, f'RMSE = {rmse:.2f} °C', transform=ax1.transAxes,
             fontsize=8, color='tab:red',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # --- Panel 2: Heating power ---------------------------------------------
    ax2 = axes[1]
    ax2.plot(t_min, results['P_heating']    / 1000, color='tab:red',    lw=1.5,
             label='Heizleistung PID (sim)')
    ax2.plot(t_min, results['oven_power_pct'] * p.P_heating_max / 100_000,
             color='tab:green', lw=1.2, ls='--', label='Ofenleistung ist (CSV) [skaliert]')
    ax2.plot(t_min, results['P_losses']     / 1000, color='tab:blue',   lw=1.0,
             ls='--', label='Strahlung/Konvektion')
    ax2.plot(t_min, results['P_band']       / 1000, color='tab:purple', lw=1.0,
             ls='-.', label='Wärmestrom Band')
    ax2.plot(t_min, results['P_airknife']   / 1000, color='tab:cyan',   lw=1.0,
             ls=':', label='Air-knife')
    ax2.plot(t_min, results['P_net']        / 1000, color='black',      lw=0.8,
             alpha=0.5, label='Netto-Leistung')
    ax2.set_ylabel('Leistung [kW]')
    ax2.grid(True, alpha=0.4)

    ax2r = ax2.twinx()
    ax2r.step(t_min, results['heating_frac'] * 100, color='tab:red',
              lw=1.0, alpha=0.5, where='post', label='PID-Ausgang [%]')
    ax2r.plot(t_min, results['oven_power_pct'], color='tab:green',
              lw=1.0, alpha=0.5, ls='--', label='Ofenleistung ist [%]')
    ax2r.set_ylabel('Heizer [%]', color='gray')
    ax2r.tick_params(axis='y', labelcolor='gray')
    ax2r.set_ylim(0, 110)

    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2r.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, loc='upper center',
               bbox_to_anchor=(0.5, -0.05), ncol=4, fontsize=7, framealpha=0.9)

    # --- Panel 3: Strip running state ----------------------------------------
    ax3 = axes[2]
    ax3.fill_between(t_min, results['strip_running'], step='post',
                     color='tab:green', alpha=0.5, label='Band läuft')
    ax3.set_yticks([0, 1])
    ax3.set_yticklabels(['Stillstand', 'Betrieb'])
    ax3.set_ylabel('Zustand')
    ax3.set_xlabel('Zeit [min]')
    ax3.set_ylim(-0.1, 1.3)
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()
    print(f"Plot saved to {save_path}")


# ---------------------------------------------------------------------------
# SIMULATION INPUTS
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import os

    CSV_FILE    = 'sample_data.csv'
    MODE        = 'pid'       # 'pid' or 'open_loop'
    DT_SECONDS  = 1.0
    T_BATH_INIT = None        # None → use first CSV value

    if not os.path.exists(CSV_FILE):
        print(f"'{CSV_FILE}' not found – generating sample data first …")
        import create_sample_csv
        create_sample_csv.main()

    print(f"Running simulation in '{MODE}' mode from '{CSV_FILE}' …")
    results = run_simulation_with_data(
        csv_path    = CSV_FILE,
        dt          = DT_SECONDS,
        T_bath_init = T_BATH_INIT,
        mode        = MODE,
    )

    T_final = results['T_bath'][-1]
    rmse    = float(np.sqrt(np.mean((results['T_bath'] - results['actual_temp'])**2)))
    print(f"Simulation complete.")
    print(f"  Final bath temperature : {T_final:.2f} °C")
    print(f"  RMSE vs measured       : {rmse:.2f} °C")
    print(f"  Mean heating power     : {results['P_heating'].mean()/1000:.1f} kW")
    print()

    plot_results_comparison(results, save_path='simulation_results.png')

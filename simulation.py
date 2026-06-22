"""
simulation.py
=============
Simulation runner for the Verzinnungsanlage (tin-coating line).

Physics model  → furnace_model.py   (thermal sub-models, TransportDelay)
PID controller → pid_controller.py  (gain-scheduled PID, anti-windup)
Parameters     → parameters.py      (all physical constants and tuning tables)
CSV loading    → data_loader.py     (ibA export parser and interpolator)

Two operating modes
-------------------
1. Manual / open-loop  → run_simulation()
   Heating driven by a time-varying schedule array (defined in parameters.py).

2. CSV + PID / open-loop → run_simulation_with_data()
   Process variables come from a real ibA CSV file.
   In 'pid' mode   : a gain-scheduled PID drives the heater.
   In 'open_loop'  : the measured Ofenleistung replays the actual heater profile.
   Simulated temperature is plotted against the measured bath temperature.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import parameters as p
import data_loader as dl
from furnace_model import (
    lookup1d,
    airknife_percent,
    strip_heat_flow,
    strip_mass_flow,
    korr_heizleistung,
    select_setpoint,
    TransportDelay,
    PT1,
)
from pid_controller import PIDController


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
        # Fluxes are evaluated at the current state T (= temperature at t_arr[i]).
        # T_arr[0] holds the initial condition T0; the integration below advances
        # T from t_arr[i] to t_arr[i+1].
        P_demand = sched[i] * p.P_heating_max
        P_h      = delay_heat.step(P_demand)
        P_loss   = lookup1d(p.loss_temp_bp, p.loss_power_tbl, T)
        Q_raw    = strip_heat_flow(T, T_in, v, bw, bt, p.cp_strip, p.rho_strip) if strip_running else 0.0
        P_band   = delay_band.step(Q_raw)
        P_wc     = delay_water.step(P_w)
        P_net    = P_h - P_loss - P_band - P_ak - P_wc

        # Store state and fluxes at time t_arr[i] (before advancing)
        T_arr[i]       = T;       frac_arr[i]   = sched[i]
        P_heat_arr[i]  = P_h;    P_loss_arr[i] = P_loss
        P_band_arr[i]  = P_band; P_ak_arr[i]   = P_ak
        P_water_arr[i] = P_wc;   P_net_arr[i]  = P_net

        # Forward-Euler step: ΔT = P_net/C_B * dt, added to the previous T
        T += (P_net / p.C_B) * step

    return dict(t=t_arr, T_bath=T_arr, heating_frac=frac_arr,
                P_heating=P_heat_arr, P_losses=P_loss_arr, P_band=P_band_arr,
                P_airknife=P_ak_arr, P_water=P_water_arr, P_net=P_net_arr)


# ---------------------------------------------------------------------------
# MODE 2 – CSV-driven simulation (PID or open-loop replay)
# ---------------------------------------------------------------------------
def run_simulation_with_data(
        csv_path: str,
        dt: float = 1.0,
        T_bath_init: float = None,
        mode: str = 'pid',
        pid_manual: bool = False,
        pid_Kp: float = None,
        pid_Ti: float = None,
        pid_Td: float = None,
):
    """
    Simulation driven by a real ibA CSV export.

    Parameters
    ----------
    csv_path    : path to the ibA CSV file (see data_loader.py for format)
    dt          : simulation time step [s]
    T_bath_init : initial bath temperature [°C]; None → use first CSV value
    mode        : 'pid'       – PID controller drives the heater
                  'open_loop' – Ofenleistung from CSV drives the heater directly
    pid_manual  : if True, use fixed PID parameters instead of the Massenstrom table
                  (only relevant when mode='pid')
    pid_Kp,     : manual ZEITBEREICH parameters (Kp, Ti [s], Td [s]). If None
    pid_Ti,       while pid_manual=True, fall back to
    pid_Td        p.pid_manual_Kp / _Ti / _Td.

    Notes on open_loop mode
    -----------------------
    The CSV column 'Ofenleistung' is the total electrical power drawn by the
    induction heater [kW]. Only a fraction of that reaches the tin bath
    (p.furnace_efficiency). The effective bath heating power is:

        P_h_eff = oven_power_W * p.furnace_efficiency
    """
    df      = dl.load_csv(csv_path)
    t_end   = df['t_s'].iloc[-1]
    n_steps = int(round(t_end / dt)) + 1
    t_grid  = np.linspace(0.0, (n_steps - 1) * dt, n_steps)
    data    = dl.interpolate_to_grid(df, t_grid)

    T0 = T_bath_init if T_bath_init is not None else float(data['actual_temp'][0])

    # Effective bath power at t=0 (pre-loads delay buffers and PID)
    init_P_eff = data['oven_power_W'][0] * p.furnace_efficiency
    init_frac  = min(1.0, init_P_eff / p.P_heating_max)

    # Allocate result arrays
    T_arr       = np.zeros(n_steps)
    frac_arr    = np.zeros(n_steps)
    P_heat_arr  = np.zeros(n_steps)
    P_loss_arr  = np.zeros(n_steps)
    P_band_arr  = np.zeros(n_steps)
    P_ak_arr    = np.zeros(n_steps)
    P_water_arr = np.zeros(n_steps)
    P_net_arr   = np.zeros(n_steps)
    mdot_arr    = np.zeros(n_steps)

    # PID – pre-loaded so its first output starts from the initial heater fraction.
    # manual=True uses fixed gains (pid_Kp/Ki/Kd); otherwise the scheduled table.
    pid = PIDController(dt, manual=pid_manual,
                        Kp=pid_Kp, Ti=pid_Ti, Td=pid_Td)
    # PID output is in percent [0..100]; preload integrator with the initial
    # heater command in percent (init_frac is a 0..1 fraction).
    pid.reset(initial_output=init_frac * 100.0, measurement=T0)

    # Heater → bath thermal lag as a PT1 (first-order lag) instead of a pure
    # dead-time delay: the heating power approaches the demand smoothly rather
    # than jumping after a fixed delay. A single PT1 whose time constant
    # switches between Betrieb (tau_heating) and Stillstand (tau_heating_stillstand).
    pt1_heat = PT1(p.tau_heating, dt, 0.0)
    # Strip heat extraction keeps its pure transport delay (Bandleistung).
    delay_band = TransportDelay(p.delay_band, dt, 0.0)

    T = T0

    for i in range(n_steps):
        # Read CSV inputs for this time step
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

        T_soll = select_setpoint(running, T_soll_b, T_soll_s)

        # Strip mass flow (Massenstrom) – used for PID Kp scheduling
        mdot = strip_mass_flow(v, bw, bt, p.rho_strip) if running else 0.0

        # Heater demand -------------------------------------------------------
        if mode == 'pid':
            pct      = pid.step(T_soll, T, mass_flow=mdot)  # PID output in % (0..100)
            raw_frac = pct / 100.0                          # → heater fraction (0..1)
            frac     = korr_heizleistung(raw_frac, running, s_len)
            P_demand = frac * p.P_heating_max
        else:
            # Open-loop: use measured Ofenleistung directly in watts.
            # furnace_efficiency scales electrical → effective bath power.
            # No P_heating_max cap so the full measured profile is replayed.
            P_demand_raw = data['oven_power_W'][i] * p.furnace_efficiency
            # Apply standstill power limit (physical PLC constraint)
            limit_W  = p.heating_limit_standstill * p.P_heating_max
            P_demand = min(P_demand_raw, limit_W) if (not running or s_len < p.strip_length_min) else P_demand_raw
            frac     = P_demand / p.P_heating_max  # stored for reference; may exceed 1.0

        # PT1 first-order lag – time constant depends on mode (Betrieb/Stillstand)
        pt1_heat.tau = p.tau_heating if running else p.tau_heating_stillstand
        P_h = pt1_heat.step(P_demand)

        # Strip heat extraction
        Q_raw  = strip_heat_flow(T, T_in, v, bw, bt, p.cp_strip, p.rho_strip) if running else 0.0
        P_band = delay_band.step(Q_raw)

        # Air-knife cooling
        ak_pct = airknife_percent(ak_os, ak_us)
        P_ak   = lookup1d(p.airknife_state_bp, p.airknife_power_tbl, ak_pct)

        # Radiation + convection losses (lookup table, at current state T)
        P_loss = lookup1d(p.loss_temp_bp, p.loss_power_tbl, T)
        P_constant = p.loss_constant  # constant offset (added to the lookup table output)

        # Net power at current state T (= temperature at t_grid[i])
        P_net = P_h - P_loss - P_band - P_ak - P_w - P_constant

        # Store state and fluxes at time t_grid[i] (before advancing).
        # T_arr[0] holds the initial condition T0.
        T_arr[i]       = T;       frac_arr[i]   = frac
        P_heat_arr[i]  = P_h;    P_loss_arr[i] = P_loss
        P_band_arr[i]  = P_band; P_ak_arr[i]   = P_ak
        P_water_arr[i] = P_w;    P_net_arr[i]  = P_net
        mdot_arr[i]    = mdot

        # Forward-Euler step: ΔT = P_net/C_B * dt, added to the previous T
        T += (P_net / p.C_B) * dt

    return dict(
        t            = t_grid,
        T_bath       = T_arr,
        heating_frac = frac_arr,
        P_heating    = P_heat_arr,
        P_losses     = P_loss_arr,
        P_band       = P_band_arr,
        P_airknife   = P_ak_arr,
        P_water      = P_water_arr,
        P_net        = P_net_arr,
        mdot         = mdot_arr,   # Massenstrom [kg/s] per step
        # From CSV (for comparison)
        actual_temp   = data['actual_temp'],
        setpoint_temp = np.array([
            select_setpoint(data['strip_running'][i] > 0.5,
                            data['setpoint_temp'][i],
                            data['setpoint_standstill'][i])
            for i in range(n_steps)
        ]),
        strip_running = data['strip_running'],
        oven_power_W  = data['oven_power_W'],   # measured electrical power [W]
    )


# ---------------------------------------------------------------------------
# Plot – open-loop results
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
    ax2.plot(t_min, results['P_losses']   / 1000, color='tab:blue',   lw=1.2,
             ls='--', label='Strahlung/Konvektion')
    ax2.plot(t_min, results['P_band']     / 1000, color='tab:green',  lw=1.2,
             ls='-.', label='Wärmestrom Band')
    ax2.plot(t_min, results['P_airknife'] / 1000, color='tab:purple', lw=1.2,
             ls=':', label='Air-knife Kühlung')
    ax2.plot(t_min, results['P_water']    / 1000, color='tab:cyan',   lw=1.2,
             ls=':', label='Wasserkühlung')
    ax2.plot(t_min, results['P_net']      / 1000, color='black',      lw=1.0,
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
# Plot – CSV-driven comparison (simulated vs. measured)
# ---------------------------------------------------------------------------
def plot_results_comparison(results: dict,
                             save_path: str = 'simulation_results.png',
                             mode: str = 'pid'):
    """
    Four-panel comparison plot:
      1. Temperatures      – simulated vs. measured + setpoint
      2. Heating power      – effective bath power (sim) vs. scaled CSV Ofenleistung
      3. Temperature error  – (simulated − measured) at every time point
      4. Strip running state
    """
    t_min = results['t'] / 60.0

    fig, axes = plt.subplots(4, 1, figsize=(14, 13), sharex=True,
                              gridspec_kw={'height_ratios': [3, 3, 1.6, 1]})

    # --- Panel 1: Temperatures ----------------------------------------------
    ax1 = axes[0]
    ax1.plot(t_min, results['actual_temp'],   color='tab:green', lw=1.5,
             label='Isttemperatur (Messung)')
    ax1.plot(t_min, results['T_bath'],        color='tab:red',   lw=1.8,
             label='Zinnbad-Temperatur (Simulation)')
    ax1.plot(t_min, results['setpoint_temp'], color='black',     lw=1.0,
             ls='--', alpha=0.7, label='Solltemperatur')
    ax1.set_ylabel('Temperatur [°C]')
    mode_label = 'PID-Regelung' if mode == 'pid' else 'Open-Loop Replay'
    ax1.set_title(f'Verzinnungsanlage – Simulation vs. Messung ({mode_label})')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.4)

    rmse = float(np.sqrt(np.mean((results['T_bath'] - results['actual_temp'])**2)))
    ax1.text(0.01, 0.04, f'RMSE = {rmse:.2f} °C', transform=ax1.transAxes,
             fontsize=8, color='tab:red',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # --- Panel 2: Power [kW effective] --------------------------------------
    ax2 = axes[1]
    ofen_eff_kW = results['oven_power_W'] * p.furnace_efficiency / 1000.0
    ax2.plot(t_min, results['P_heating']  / 1000, color='tab:red',    lw=1.5,
             label='Heizleistung sim [kW eff]')
    ax2.plot(t_min, ofen_eff_kW,                  color='tab:green',  lw=1.2,
             ls='--', label=f'Ofenleistung × η={p.furnace_efficiency} [kW eff]')
    ax2.plot(t_min, results['P_losses']   / 1000, color='tab:blue',   lw=1.0,
             ls='--', label='Strahlung/Konvektion')
    ax2.plot(t_min, results['P_band']     / 1000, color='tab:purple', lw=1.0,
             ls='-.', label='Wärmestrom Band')
    ax2.plot(t_min, results['P_airknife'] / 1000, color='tab:cyan',   lw=1.0,
             ls=':', label='Air-knife')
    ax2.plot(t_min, results['P_net']      / 1000, color='black',      lw=0.8,
             alpha=0.5, label='Netto-Leistung')
    ax2.set_ylabel('Leistung [kW eff.]')
    ax2.grid(True, alpha=0.4)

    ax2r = ax2.twinx()
    ax2r.step(t_min, results['heating_frac'] * 100, color='tab:red',
              lw=1.0, alpha=0.5, where='post', label='Heizer-Ausgang [%]')
    ax2r.set_ylabel('Heizer-Ausgang [%]', color='gray')
    ax2r.tick_params(axis='y', labelcolor='gray')
    ax2r.set_ylim(0, 110)

    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2r.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, loc='upper center',
               bbox_to_anchor=(0.5, -0.05), ncol=4, fontsize=7, framealpha=0.9)

    # --- Panel 3: Differences (sim−measured) and control error (Regeldiff.) --
    ax3 = axes[2]
    dT  = results['T_bath'] - results['actual_temp']            # Sim − Messung
    reg = results['setpoint_temp'] - results['T_bath']          # Regeldifferenz e = Soll − Sim
    ax3.axhline(0.0, color='black', lw=0.8, alpha=0.6)
    ax3.plot(t_min, dT,  color='tab:red',  lw=1.2, label='ΔT = Sim − Messung')
    ax3.fill_between(t_min, dT, 0.0, color='tab:red', alpha=0.2)
    ax3.plot(t_min, reg, color='tab:blue', lw=1.2, label='Regeldifferenz = Soll − Sim')
    ax3.set_ylabel('Δ [°C]')
    mae = float(np.mean(np.abs(dT)))
    ax3.text(0.01, 0.06,
             f'ΔT: mean={dT.mean():+.2f}  MAE={mae:.2f}  min={dT.min():+.2f}  max={dT.max():+.2f}  |  '
             f'Regeldiff.: mean={reg.mean():+.2f}  max|·|={np.max(np.abs(reg)):.2f} °C',
             transform=ax3.transAxes, fontsize=7.5, color='black',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    ax3.legend(loc='upper right', fontsize=8, ncol=2)
    ax3.grid(True, alpha=0.4)

    # --- Panel 4: Strip running state ----------------------------------------
    ax4 = axes[3]
    ax4.fill_between(t_min, results['strip_running'], step='post',
                     color='tab:green', alpha=0.5, label='Band läuft')
    ax4.set_yticks([0, 1])
    ax4.set_yticklabels(['Stillstand', 'Betrieb'])
    ax4.set_ylabel('Zustand')
    ax4.set_xlabel('Zeit [min]')
    ax4.set_ylim(-0.1, 1.3)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()
    print(f"Plot saved to {save_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import os

    CSV_FILE    = 'CUSN6_0,4-Ofen1.txt'
    MODE        = 'pid' # 'pid' or 'open_loop'
    DT_SECONDS  = 1
    T_BATH_INIT = p.T_bath_init   # initial bath temperature from parameters.py
                                  # (set to None to start from the first CSV value)

    # --- PID parameter source (only used when MODE == 'pid') ----------------
    # PID_MANUAL = False → Kp scheduled from the Massenstrom table, Ki/Kd constant
    # PID_MANUAL = False  → use the fixed ZEITBEREICH params below (bypasses table)
    PID_MANUAL  = False     # manual params ON – set False to use the scheduled table
    PID_KP      = 5.0     # proportional gain          (used only if PID_MANUAL)
    PID_TI      = 220.0    # integration time [s]       (used only if PID_MANUAL)
    PID_TD      = 40.0     # derivative time  [s]       (used only if PID_MANUAL)

    if not os.path.exists(CSV_FILE):
        raise FileNotFoundError(
            f"'{CSV_FILE}' not found. Place the ibA export in the working directory."
        )

    src = 'manual gains' if (MODE == 'pid' and PID_MANUAL) else \
          'scheduled gains' if MODE == 'pid' else 'CSV power'
    print(f"Running simulation in '{MODE}' mode ({src}) from '{CSV_FILE}' …")
    results = run_simulation_with_data(
        csv_path    = CSV_FILE,
        dt          = DT_SECONDS,
        T_bath_init = T_BATH_INIT,
        mode        = MODE,
        pid_manual  = PID_MANUAL,
        pid_Kp      = PID_KP,
        pid_Ti      = PID_TI,
        pid_Td      = PID_TD,
    )

    T_final  = results['T_bath'][-1]
    mdot_avg = float(results['mdot'].mean())
    rmse     = float(np.sqrt(np.mean((results['T_bath'] - results['actual_temp'])**2)))
    print(f"Simulation complete.")
    print(f"  Final bath temperature : {T_final:.2f} °C")
    print(f"  Avg. Massenstrom       : {mdot_avg:.3f} kg/s")
    print(f"  RMSE vs measured       : {rmse:.2f} °C")
    print(f"  Mean heating power     : {results['P_heating'].mean()/1000:.1f} kW  (effective bath)")
    print(f"  Mean Ofenleistung      : {results['oven_power_W'].mean()/1000:.1f} kW  (electrical)")
    print()

    plot_results_comparison(results, save_path='simulation_results.png', mode=MODE)

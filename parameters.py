"""
parameters.py
=============
All physical and operational parameters for the Verzinnungsanlage (tin-coating line)
thermal simulation.

Unit conventions
----------------
Temperature  : °C
Power        : W
Mass         : kg
Length       : m
Time         : s
"""

import numpy as np

# ---------------------------------------------------------------------------
# Tin bath (Zinnbad) – pure liquid tin above 231.93 °C
# ---------------------------------------------------------------------------
rho_Sn = 6960          # kg/m³  – density at ~270 °C
cp_Sn  = 230           # J/(kg·K) – specific heat capacity (liquid)

# Bath geometry
l_B = 1.44             # m  – length
b_B = 0.89             # m  – width
h_B = 0.86             # m  – height
V_B = l_B * b_B * h_B # m³ – volume

m_B = V_B * rho_Sn    # kg – tin mass
C_B = m_B * cp_Sn     # J/K – effective thermal capacity of the bath

# Alternatively the mass can be set directly from a measured fill level.
# The Simulink model uses 7671 kg (system_1659) / 7700 kg (system_1904).
# Uncomment to override the geometry-derived value:
# m_B = 7671
# C_B = m_B * cp_Sn

# ---------------------------------------------------------------------------
# Initial condition
# ---------------------------------------------------------------------------
T_bath_init = 280.0    # °C – initial tin-bath temperature

# ---------------------------------------------------------------------------
# Ambient (Störgröße – disturbance)
# ---------------------------------------------------------------------------
T_ambient = 35.0       # °C – ambient temperature

# ---------------------------------------------------------------------------
# Steel strip (Band) parameters  –  used when strip is running
# ---------------------------------------------------------------------------
rho_strip = 8800.0     # kg/m³  – mean density (as used in Simulink constant)
cp_strip  = 500.0      # J/(kg·K) – mean specific heat of the strip material
                       # (temperature-averaged; exact value depends on alloy)

# Default strip dimensions & speed (can be overridden at runtime)
strip_width     = 350.0   # mm
strip_thickness = 0.5     # mm
strip_speed     = 60.0    # m/min

# Temperature of the strip entering the bath (below the roller, "unter Rolle").
# The strip is pre-heated in the annealing furnace and enters the tin bath
# close to bath temperature.  Typical delta = 5–15 °C below the bath setpoint.
T_strip_in = T_ambient     # °C – strip temperature just before entering the bath

# ---------------------------------------------------------------------------
# Transport delays (Transportverzögerungen)
# ---------------------------------------------------------------------------
# Values match the Simulink with-data subsystem (system_1904):
delay_heating     = 200.0   # s – inductive heater → bath (Heizleistung Betrieb)
delay_band        = 12.0   # s – strip-heat calculation → bath (Bandleistung)
delay_water_cool  = delay_heating  # s – water cooling → bath
delay_airknife    = 0.0    # s – air-knife cooling → bath (instantaneous)

# ---------------------------------------------------------------------------
# Heat losses: radiation + convection  (Wärmeverluste Strahlung + Konvektion)
# Lookup table – function of bath temperature [°C] → power loss [W]
# Source: Simulink Lookup_n-D block in the main thermal subsystem
# ---------------------------------------------------------------------------
loss_temp_bp    = [  0,  250,  260,  270,  280,  290,  300,  310,  320]  # °C
loss_power_tbl  = [  0, 20500, 22000, 23500, 25000, 27500, 30000, 33000, 37000]  # W
loss_constant = 95000.0  # W – constant offset (added to the lookup table output)

# ---------------------------------------------------------------------------
# Air-knife cooling  (Kühlleistung Air knife)
# Lookup table – operating state [%] → cooling power [W]
# ---------------------------------------------------------------------------
airknife_state_bp   = [  0,  50,   70,   80,   90,  100]   # %
airknife_power_tbl  = [  0, 25000, 26500, 26700, 26900, 26800]  # W

# Air-knife operating pressures (mbar) used to compute operating state [%]
# State (%) = 100 * mean_pressure / 600  (clamped 0–100)
airknife_pressure_min   = 20    # mbar – below this → treated as 0
airknife_pressure_max   = 700   # mbar – saturation

# Default operating pressures for each air-knife nozzle when running open-loop
airknife_druck1_default = 590.0  # mbar (OS – Oberseite)
airknife_druck2_default = 590.0  # mbar (US – Unterseite)

# ---------------------------------------------------------------------------
# Water cooling  (Wasserkühlleistung)
# ---------------------------------------------------------------------------
water_cooling_power_default = 10000.0  # W – off by default; set > 0 to enable

# ---------------------------------------------------------------------------
# Inductive heating  (Heizleistung)
# The heater output is specified as a fraction [0 … 1] of the installed
# maximum power.
#
# Sizing note: at steady state the heater must cover all losses:
#   - Radiation + convection at 270 °C : ~23 500 W
#   - Air-knife cooling at 50 % state  : ~25 000 W
#   - Strip heat flow at typical conditions: ~20 000–40 000 W
#   → total ≈ 70–90 kW  → installed max ~100 kW
# ---------------------------------------------------------------------------
P_heating_max = 700_000.0  # W – effective bath heating power at 100 % output
                           # (≠ total electrical power; see furnace_efficiency below)

# Inductive heating fraction [0..1] for open-loop (manual) operation.
# 0.75 → 90 kW, which roughly balances losses at ~270 °C with strip running.
heating_fraction_default = 0.75

# Furnace thermal efficiency
# Fraction of total electrical power (Ofenleistung) that reaches the tin bath.
# The remainder heats the furnace chamber, structural steel, and surrounding air.
# Calibrated from ibA data: at steady state Ofenleistung ≈ 200 kW maintains
# ~270 °C with ~50 kW total bath losses  →  η ≈ 50/200 = 0.25
# Used in open_loop mode:  P_h_eff = oven_power_W * furnace_efficiency
furnace_efficiency = 0.9

# ---------------------------------------------------------------------------
# Transport delays – mode-dependent  (from Simulink system_1904.xml)
# ---------------------------------------------------------------------------
delay_heating_stillstand = 150.0  # s – longer delay when strip is not running

# ---------------------------------------------------------------------------
# Heater → bath thermal lag as a PT1 (first-order lag) element
# Replaces the pure transport (dead-time) delay on the heating path: the
# heating power approaches the demand smoothly (τ·dy/dt = u − y) instead of
# jumping after a fixed dead time. τ is mode-dependent (Betrieb / Stillstand).
# ---------------------------------------------------------------------------
tau_heating            = 150.0    # s – PT1 time constant, Betrieb (strip running)
tau_heating_stillstand = 150.0   # s – PT1 time constant, Stillstand

# Minimum strip length in oven before full heating is allowed
strip_length_min        = 50.0    # m  (korrHeizleistung MATLAB function)
heating_limit_standstill = 0.20   # max heater fraction during standstill / short strip

# ---------------------------------------------------------------------------
# PID controller – PARALLEL form
#
# Control law (parallel / ideal-parallel form):
#     u = Kp·e + Ki·∫e dt + Kd·de/dt
#
# Gain handling:
#   Kp : PROPORTIONAL gain – scheduled by the strip MASS FLOW (Massenstrom),
#        looked up from pid_Kp_table over pid_massflow_bp
#   Ki : INTEGRAL gain (parallel)   – CONSTANT
#   Kd : DERIVATIVE gain (parallel) – CONSTANT
#
# The time-domain (Zeitbereich) representations are DERIVED from the above:
#     Ti = Kp / Ki      (integration time  [s])
#     Td = Kd / Kp      (derivative time   [s])
# Because Kp is mass-flow-scheduled, Ti and Td vary with the Massenstrom even
# though Ki and Kd are constant.
# ---------------------------------------------------------------------------
# Massenstrom breakpoints: 0.0 … 4.5 kg/s in 0.1 kg/s steps (46 points)
pid_massflow_bp = [round(0.1 * i, 1) for i in range(46)]

# P-Anteil – proportional gain scheduled by Massenstrom [kg/s]
# Linear ramp Kp = 5.00 + 0.667·Massenstrom (5.00 at 0 kg/s → 8.00 at 4.5 kg/s)
pid_Kp_table = [
    5.00, 5.07, 5.13, 5.20, 5.27, 5.33, 5.40, 5.47, 5.53, 5.60,
    5.67, 5.73, 5.80, 5.87, 5.93, 6.00, 6.07, 6.13, 6.20, 6.27,
    6.33, 6.40, 6.47, 6.53, 6.60, 6.67, 6.73, 6.80, 6.87, 6.93,
    7.00, 7.07, 7.13, 7.20, 7.27, 7.33, 7.40, 7.47, 7.53, 7.60,
    7.67, 7.73, 7.80, 7.87, 7.93, 8.00,
]

# Constant parallel-form integral and derivative gains
pid_Ki = 0.02    # integral gain   (parallel form), constant
pid_Kd = 300.0   # derivative gain (parallel form), constant

# ---------------------------------------------------------------------------
# Manual PID parameters – ZEITBEREICH (time-domain / standard form)
# Used when the controller runs in MANUAL mode (PIDController(manual=True)),
# bypassing the Massenstrom-scheduled table.
#
#   u = Kp·( e + (1/Ti)·∫e dt + Td·de/dt )
#
# The controller converts these to its internal parallel-form gains:
#     Ki = Kp / Ti      Kd = Kp · Td
# ---------------------------------------------------------------------------
pid_manual_Kp = 6.0      # proportional gain
pid_manual_Ti = 300.0    # integration time [s]  (Nachstellzeit)
pid_manual_Td = 50.0     # derivative time  [s]  (Vorhaltzeit)

# PID output saturation – PERCENT [0..100], matching the Simulink output scale.
# The controller output is a heater command in %, converted to a 0..1 fraction
# (frac = output/100) before scaling by P_heating_max in simulation.py.
# This gives the table/manual gains (Kp ≈ 5–8) a sensible proportional band
# (e.g. Kp=6 → full output at ~16 °C error) instead of a 0.16 °C band.
pid_output_min = 0.0
pid_output_max = 100.0

# ---------------------------------------------------------------------------
# Simulation time settings
# ---------------------------------------------------------------------------
t_start = 0.0          # s
t_end   = 3600.0       # s  (default: 1 hour)
dt      = 1.0          # s  – integration time step

# ---------------------------------------------------------------------------
# Heating schedule  (Heizleistungs-Zeitplan)
# ---------------------------------------------------------------------------
# Instead of a constant fraction, the inductive heater follows a time-varying
# schedule defined as an array of fractions [0..1], one value per time step.
#
# generate_heating_schedule() builds the array for a given simulation window.
# The schedule is piecewise-constant: power changes every `step_duration_s`
# seconds and is held flat within each segment (realistic for PLC setpoints).
#
# Parameters you can tune:
heating_schedule_min           = 0.40   # minimum heater fraction (40 % of P_max)
heating_schedule_max           = 1.00   # maximum heater fraction (100 % of P_max)
heating_schedule_step_duration = 300    # s – how often the power level changes (5 min)
heating_schedule_seed          = 42     # random seed for reproducibility


def generate_heating_schedule(n_steps: int,
                               step_duration_s: int  = None,
                               frac_min: float       = None,
                               frac_max: float       = None,
                               seed: int             = None) -> np.ndarray:
    """
    Build a piecewise-constant random heating-fraction schedule.

    Returns an array of length `n_steps` with values in [frac_min, frac_max].
    Each segment lasts `step_duration_s` time steps before changing.
    """
    seg_dur  = step_duration_s if step_duration_s is not None else heating_schedule_step_duration
    lo       = frac_min        if frac_min        is not None else heating_schedule_min
    hi       = frac_max        if frac_max        is not None else heating_schedule_max
    rng      = np.random.default_rng(seed if seed is not None else heating_schedule_seed)

    n_segments  = int(np.ceil(n_steps / seg_dur))
    seg_values  = rng.uniform(lo, hi, size=n_segments)

    # Repeat each segment value for `seg_dur` steps, then trim to n_steps
    schedule = np.repeat(seg_values, seg_dur)[:n_steps]
    return schedule


# Pre-built schedule for the default simulation window (t_end, dt defined above)
# Re-call generate_heating_schedule() if you change t_end or dt.
_default_n_steps = int(round((t_end - t_start) / dt)) + 1
heating_schedule = generate_heating_schedule(_default_n_steps)

# ---------------------------------------------------------------------------
# Vessel wall (Wanne) – 2-node thermal model
# ---------------------------------------------------------------------------
# The vessel wall is modelled as two thermal nodes:
#
#   T_bath ─[G_tin_inner]─ T_wanne_innen ─[G_wall]─ T_wanne_außen ─[G_outer]─ T_ambient
#
# Node 1 (innen): inner surface of the steel vessel wall
# Node 2 (außen): outer surface of the insulation layer
#
# Heat flows are computed from bath → inner wall → outer insulation → ambient.
# These nodes are *observational*: they do not feed back into the bath energy
# balance (which already accounts for total losses via the lookup table).

# --- Vessel geometry ---------------------------------------------------------
# Contact area between liquid tin and vessel walls (bottom + 4 sides)
A_wanne = 2 * (l_B * h_B + b_B * h_B) + l_B * b_B   # m²

# --- Steel vessel wall -------------------------------------------------------
d_steel      = 0.020    # m   – wall thickness
lambda_steel = 50.0     # W/(m·K) – thermal conductivity
rho_steel    = 7800.0   # kg/m³
cp_steel     = 500.0    # J/(kg·K)

m_wanne_innen = A_wanne * d_steel * rho_steel   # kg – steel wall mass
C_wanne_innen = m_wanne_innen * cp_steel         # J/K – thermal capacity of inner wall

# --- Insulation layer (außen) ------------------------------------------------
d_ins      = 0.050      # m   – insulation thickness (50 mm mineral wool)
lambda_ins = 0.10       # W/(m·K) – thermal conductivity of insulation
rho_ins    = 200.0      # kg/m³
cp_ins     = 840.0      # J/(kg·K)

m_wanne_aussen = A_wanne * d_ins * rho_ins   # kg
C_wanne_aussen = m_wanne_aussen * cp_ins      # J/K

# --- Heat transfer coefficients ----------------------------------------------
h_tin_to_wall  = 1000.0  # W/(m²·K) – liquid tin → inner steel surface (metal convection)
h_outer_to_air =   10.0  # W/(m²·K) – outer insulation surface → ambient air

# --- Derived thermal conductances [W/K] --------------------------------------
# From tin bath to inner wall (convection)
G_tin_inner = h_tin_to_wall * A_wanne

# Through steel wall + insulation (series conduction)
G_wall = A_wanne / (d_steel / lambda_steel + d_ins / lambda_ins)

# From outer insulation surface to ambient (convection)
G_outer = h_outer_to_air * A_wanne

# --- Initial temperatures ----------------------------------------------------
T_wanne_innen_init = T_bath_init   # °C – starts in thermal equilibrium with bath
T_wanne_aussen_init = T_ambient    # °C – starts at ambient

# ---------------------------------------------------------------------------
# Alloy number → specific heat (cp) lookup  [example values for common alloys]
# In the Simulink model, cp at 20 °C and 400 °C are fetched from a
# reference table indexed by Bandnummer (strip/alloy number).
# The mean cp used in the enthalpy calculation is:
#   cp_m = cp_20 + (T_ref – T_ambient) / (400 – 20) * (cp_400 – cp_20)
# Values below are representative; refine with your material data.
# ---------------------------------------------------------------------------
# Alloy index: 0 = generic steel, 1 = DC01, 2 = DX54D, ...
alloy_cp_20  = {0: 481, 1: 490, 2: 495}  # J/(kg·K) at 20 °C
alloy_cp_400 = {0: 620, 1: 630, 2: 625}  # J/(kg·K) at 400 °C

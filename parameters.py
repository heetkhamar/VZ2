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
T_bath_init = 264.0    # °C – initial tin-bath temperature

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
delay_heating     = 120.0  # s – inductive heater → bath (thermal lag in the system)
delay_band        = 5.0   # s – strip-heat calculation → bath integration
delay_water_cool  = delay_heating  # s – water cooling → bath
delay_airknife    = 0.0   # s – air-knife cooling → bath (instantaneous)

# ---------------------------------------------------------------------------
# Heat losses: radiation + convection  (Wärmeverluste Strahlung + Konvektion)
# Lookup table – function of bath temperature [°C] → power loss [W]
# Source: Simulink Lookup_n-D block in the main thermal subsystem
# ---------------------------------------------------------------------------
loss_temp_bp    = [  0,  250,  260,  270,  280,  290,  300,  310,  320]  # °C
loss_power_tbl  = [  0, 20500, 22000, 23500, 25000, 27500, 30000, 33000, 37000]  # W

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
water_cooling_power_default = 0.0  # W – off by default; set > 0 to enable

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
P_heating_max = 120_000.0  # W – maximum installed inductive heating power (120 kW)

# Inductive heating fraction [0..1] for open-loop (manual) operation.
# 0.75 → 90 kW, which roughly balances losses at ~270 °C with strip running.
heating_fraction_default = 0.75

# ---------------------------------------------------------------------------
# Transport delays – mode-dependent  (from Simulink system_1904.xml)
# ---------------------------------------------------------------------------
delay_heating_stillstand = 150.0  # s – longer delay when strip is not running

# Minimum strip length in oven before full heating is allowed
strip_length_min        = 50.0    # m  (korrHeizleistung MATLAB function)
heating_limit_standstill = 0.20   # max heater fraction during standstill / short strip

# ---------------------------------------------------------------------------
# PID controller – gain-scheduled by strip speed
# Source: Simulink 1-D Lookup Tables in system_1349.xml
#
# Breakpoints: strip speed [m/s] at 0.3 m/s spacing, 0 → 4.5 m/s (16 pts)
# Kp  : proportional gain          (P-Anteil)
# Ti  : integration time [s]       (I-Anteil Zeitdarstellung)
# Td  : derivative time [s]        (D-Anteil Zeitdarstellung)
# ---------------------------------------------------------------------------
pid_speed_bp = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1,
                2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.2, 4.5]  # m/s

pid_Kp_table = [5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4,
                6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0]

pid_Ti_table = [250, 260, 270, 280, 290, 300, 310, 320,
                330, 340, 350, 360, 370, 380, 390, 400]   # s

pid_Td_table = [60.0, 57.7, 55.6, 53.6, 51.7, 50.0, 48.4, 46.9,
                45.5, 44.1, 42.9, 41.7, 40.5, 39.5, 38.5, 37.5]  # s

# PID output saturation (heater fraction [0..1])
pid_output_min = 0.0
pid_output_max = 1.0

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

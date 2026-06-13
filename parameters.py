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
delay_water_cool  = dealy_heating  # s – water cooling → bath
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
P_heating_max = 700_000.0  # W – maximum installed inductive heating power

# Inductive heating fraction [0..1] for open-loop (manual) operation.
# 0.75 → 75 kW, which roughly balances losses at ~270 °C with strip running.
heating_fraction_default = 0.75

# ---------------------------------------------------------------------------
# Simulation time settings
# ---------------------------------------------------------------------------
t_start = 0.0          # s
t_end   = 3600.0       # s  (default: 1 hour)
dt      = 1.0          # s  – integration time step

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

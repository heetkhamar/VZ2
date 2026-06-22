"""
data_loader.py
==============
Load and preprocess ibA CSV exports for the Verzinnungsanlage simulation.

Default format (CUSN6_0,4-Ofen1.txt and similar ibA exports):
  - Separator      : semicolon (;)
  - Decimal        : dot (.)
  - Encoding       : Latin-1 (Windows-1252 compatible)
  - Time column    : 'time'  – already integer elapsed seconds from 0
  - Power columns  : Ofenleistung [kW], Kühlleistung [kW]  → converted to W

Adapt CSV_COLUMNS below to match your ibA export's column headers.
"""

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Column name mapping  (internal name → CSV column header)
# Edit right-hand values to match your ibA export.
# ---------------------------------------------------------------------------
CSV_COLUMNS = {
    # Time – already integer seconds in the ibA export
    't_s':                      'time',

    # Strip geometry
    'strip_thickness':          'Banddicke',                                    # mm
    'strip_width':              'Bandbreite',                                   # mm

    # Heating – Ofenleistung is the total electrical power [kW]; converted to W below
    'oven_power_kW':            'Ofenleistung',                                 # kW → W

    # Temperatures
    'actual_temp':              'Wannentemperatur innen',                       # °C (bath, inner sensor)
    'T_wanne_aussen':           'Wannentemperatur außen',                       # °C (bath, outer sensor)
    'setpoint_temp':            'Laufprogramm Online Zinnbad / Anlagenstartwert korrigiert',  # °C
    'setpoint_standstill':      'Laufprogramm Online Zinnbad / Anlagenstartwert Rezept',     # °C
    'temp_under_roller':        'Temperatur unter Zinnbadrolle',                # °C (strip under roller)
    'temp_under_roller_setpoint': 'Solltemperatur unter Zinnbadrolle',          # °C (strip setpoint)

    # Air knife
    'airknife_os_mbar':         'Abblasung OS Luftdruck',                       # mbar
    'airknife_us_mbar':         'Abblasung US Luftdruck',                       # mbar

    # Water cooling circuit 1 – Kühlleistung [kW] converted to W below
    'water_cooling_kW':         'Verzinnungsstation: Kühlleistung',             # kW → W
    'coolant_flow_1_lpm':       'Verzinnungsstation: Durchflussmenge',          # L/min
    'coolant_return_temp_1':    'Verzinnungsstation: Rücklauftemperatur',       # °C

    # Water cooling circuit 2 (FU – frequency-converter side)
    'coolant_flow_2_lpm':       'Verzinnungsstation FU: Durchflussmenge',       # L/min
    'coolant_return_temp_2':    'Verzinnungsstation FU: Rücklauftemperatur',    # °C

    # Strip process
    'strip_speed':              'Bandgeschwindigkeit',                          # m/min
    'strip_running':            'Bandtransport',                                # 0 / 1
    'strip_length_in_oven':     'Bandlänge auf Abhaspel',                       # m (uncoiler proxy)
    'traverse_down':            'Hubtraverse unten',                            # 0 / 1

    # Oven / furnace active flags
    'oven_1_active':            'MON Verzinnungseinheit Verzinnungsstation Ofen 1 aktiv',  # 0 / 1
    'oven_2_active':            'MON Verzinnungseinheit Verzinnungsstation Ofen 2 aktiv',  # 0 / 1

    # Material
    'alloy_nr':                 'Legierung',                                    # text / int
    # ChannelInfoFieldText is text metadata – intentionally not mapped
}

_INV_MAP = {v: k for k, v in CSV_COLUMNS.items()}

# Columns that must be present after loading; all others fall back to defaults.
REQUIRED_COLUMNS = {'actual_temp', 'setpoint_temp'}

COLUMN_DEFAULTS = {
    'strip_speed':                  60.0,
    'strip_width':                  360.0,
    'strip_thickness':              0.4,
    'oven_power_W':                 0.0,
    'setpoint_standstill':          265.0,
    'airknife_os_mbar':             310.0,
    'airknife_us_mbar':             310.0,
    'water_cooling_W':              0.0,
    'temp_under_roller':            265.0,
    'temp_under_roller_setpoint':   265.0,
    'T_wanne_aussen':               35.0,
    'coolant_flow_1_lpm':           0.0,
    'coolant_return_temp_1':        25.0,
    'coolant_flow_2_lpm':           0.0,
    'coolant_return_temp_2':        25.0,
    'traverse_down':                0.0,
    'oven_1_active':                0.0,
    'oven_2_active':                0.0,
    'alloy_nr':                     0.0,
    'strip_running':                1.0,
    'strip_length_in_oven':         100.0,
}


def load_csv(filepath: str,
             sep: str = ';',
             decimal: str = '.',
             encoding: str = 'latin-1') -> pd.DataFrame:
    """
    Load an ibA CSV export and return a DataFrame with internal column names
    and an elapsed-time column 't_s' [s].

    Parameters
    ----------
    filepath : path to the CSV file
    sep      : column separator (';' for ibA exports)
    decimal  : decimal separator ('.' for ibA exports, ',' for sample_data.csv)
    encoding : file encoding ('latin-1' handles German umlauts in ibA exports)
    """
    df = pd.read_csv(filepath, sep=sep, decimal=decimal, encoding=encoding)

    # Rename CSV column headers → internal names
    df = df.rename(columns=_INV_MAP)

    # Time column handling
    if 't_s' in df.columns:
        # ibA export: 'time' column already contains elapsed seconds
        df['t_s'] = pd.to_numeric(df['t_s'], errors='coerce').fillna(0.0)
    elif 'timestamp' in df.columns:
        # Datetime string format (sample_data.csv)
        df['timestamp'] = pd.to_datetime(df['timestamp'], dayfirst=True)
        df['t_s'] = (df['timestamp'] - df['timestamp'].iloc[0]).dt.total_seconds()
    else:
        df['t_s'] = np.arange(len(df), dtype=float)

    # Unit conversions: kW → W
    if 'oven_power_kW' in df.columns:
        df['oven_power_W'] = (
            pd.to_numeric(df['oven_power_kW'], errors='coerce').fillna(0.0) * 1000.0
        )
        df = df.drop(columns=['oven_power_kW'])

    if 'water_cooling_kW' in df.columns:
        df['water_cooling_W'] = (
            pd.to_numeric(df['water_cooling_kW'], errors='coerce').fillna(0.0) * 1000.0
        )
        df = df.drop(columns=['water_cooling_kW'])

    # Fill missing optional columns with defaults and coerce to numeric
    for col, default in COLUMN_DEFAULTS.items():
        if col not in df.columns:
            df[col] = default
        else:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(default)

    # Validate required columns
    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            raise ValueError(
                f"Required column '{col}' not found in {filepath}.\n"
                f"Available columns: {list(df.columns)}"
            )

    # Ensure strip_running is binary
    df['strip_running'] = df['strip_running'].clip(0, 1).round()

    df = df.sort_values('t_s').reset_index(drop=True)
    return df


def interpolate_to_grid(df: pd.DataFrame, t_grid: np.ndarray) -> dict:
    """
    Linearly interpolate all signal columns onto the simulation time grid.

    Returns a dict {internal_name: np.ndarray of length len(t_grid)}.
    Binary signals (strip_running) are rounded back to 0/1 after interpolation.
    """
    t_src = df['t_s'].values
    signals = {}

    skip = {'timestamp', 't_s'}
    for col in df.columns:
        if col in skip or not pd.api.types.is_numeric_dtype(df[col]):
            continue
        signals[col] = np.interp(t_grid, t_src, df[col].values,
                                  left=float(df[col].iloc[0]),
                                  right=float(df[col].iloc[-1]))

    if 'strip_running' in signals:
        signals['strip_running'] = np.round(signals['strip_running']).astype(float)

    return signals

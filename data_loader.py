"""
data_loader.py
==============
Load and preprocess ibA CSV exports for the Verzinnungsanlage simulation.

Expected CSV format (semicolon-separated, German decimal comma):
    Zeitstempel;Bandgeschwindigkeit_m_min;Bandbreite_mm;Banddicke_mm;
    Ofenleistung_Prozent;Solltemperatur_Betrieb_C;Solltemperatur_Stillstand_C;
    IstTemperatur_C;LuftmesserOS_mbar;LuftmesserUS_mbar;Wasserkuehlung_W;
    TemperaturUnterRolle_C;LegierungsNr;ZustandBandlauf;BandlaengeImOfen_m

Adapt CSV_COLUMNS below if your ibA export uses different column headers.
"""

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Column name mapping: internal name → CSV header in your ibA export
# Edit the right-hand values to match your actual file headers.
# ---------------------------------------------------------------------------
CSV_COLUMNS = {
    'timestamp':              'Zeitstempel',
    'strip_speed':            'Bandgeschwindigkeit_m_min',    # m/min
    'strip_width':            'Bandbreite_mm',                # mm
    'strip_thickness':        'Banddicke_mm',                 # mm
    'oven_power_pct':         'Ofenleistung_Prozent',         # 0–100 %
    'setpoint_temp':          'Solltemperatur_Betrieb_C',     # °C
    'setpoint_standstill':    'Solltemperatur_Stillstand_C',  # °C
    'actual_temp':            'IstTemperatur_C',              # °C  (measured)
    'airknife_os_mbar':       'LuftmesserOS_mbar',            # mbar
    'airknife_us_mbar':       'LuftmesserUS_mbar',            # mbar
    'water_cooling_W':        'Wasserkuehlung_W',             # W
    'temp_under_roller':      'TemperaturUnterRolle_C',       # °C
    'alloy_nr':               'LegierungsNr',                 # integer
    'strip_running':          'ZustandBandlauf',              # 0 or 1
    'strip_length_in_oven':   'BandlaengeImOfen_m',           # m
}

_INV_MAP = {v: k for k, v in CSV_COLUMNS.items()}

# Columns that must be present; others are optional (defaults used if missing)
REQUIRED_COLUMNS = {'actual_temp', 'setpoint_temp'}

COLUMN_DEFAULTS = {
    'strip_speed':          60.0,
    'strip_width':          900.0,
    'strip_thickness':      0.5,
    'oven_power_pct':       70.0,
    'setpoint_standstill':  265.0,
    'airknife_os_mbar':     300.0,
    'airknife_us_mbar':     300.0,
    'water_cooling_W':      0.0,
    'temp_under_roller':    35.0,
    'alloy_nr':             0.0,
    'strip_running':        1.0,
    'strip_length_in_oven': 100.0,
}


def load_csv(filepath: str,
             sep: str = ';',
             decimal: str = ',') -> pd.DataFrame:
    """
    Load an ibA CSV export and return a DataFrame with internal column names
    and an elapsed-time column 't_s' [s].

    Parameters
    ----------
    filepath : path to the CSV file
    sep      : column separator (';' is standard for German ibA exports)
    decimal  : decimal separator (',' for German format)
    """
    df = pd.read_csv(filepath, sep=sep, decimal=decimal)

    # Rename CSV columns → internal names
    df = df.rename(columns=_INV_MAP)

    # Parse timestamp → elapsed seconds
    if 'timestamp' in df.columns:
        df['timestamp'] = pd.to_datetime(df['timestamp'], dayfirst=True)
        df['t_s'] = (df['timestamp'] - df['timestamp'].iloc[0]).dt.total_seconds()
    elif 't_s' not in df.columns:
        df['t_s'] = np.arange(len(df), dtype=float)

    # Fill missing optional columns with defaults
    for col, default in COLUMN_DEFAULTS.items():
        if col not in df.columns:
            df[col] = default

    # Validate required columns
    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            raise ValueError(
                f"Required column '{col}' (CSV header: '{CSV_COLUMNS[col]}') "
                f"not found in {filepath}.\n"
                f"Available columns: {list(df.columns)}"
            )

    # Ensure binary
    df['strip_running'] = df['strip_running'].clip(0, 1).round()

    df = df.sort_values('t_s').reset_index(drop=True)
    return df


def interpolate_to_grid(df: pd.DataFrame, t_grid: np.ndarray) -> dict:
    """
    Linearly interpolate all signal columns onto the simulation time grid.

    Returns a dict {internal_name: np.ndarray of length len(t_grid)}.
    Binary signals (strip_running) are rounded after interpolation.
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

    # Re-binarise after interpolation
    if 'strip_running' in signals:
        signals['strip_running'] = np.round(signals['strip_running']).astype(float)

    return signals

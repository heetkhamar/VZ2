"""
create_sample_csv.py
====================
Generates a realistic sample ibA CSV file (sample_data.csv) for testing the
Verzinnungsanlage simulation without real ibA data.

The scenario simulates 2 hours of operation:
  0–15 min   : warmup – standstill, bath heating from 264 °C to standstill setpoint
  15–90 min  : production – strip running at varying speed and width
  90–120 min : cooldown – strip stops, bath maintained at standstill setpoint

Run once to create the file:
    python create_sample_csv.py
"""

import numpy as np
import pandas as pd
import parameters as p

SEED    = 7
RNG     = np.random.default_rng(SEED)
DT_CSV  = 30        # s – ibA logging interval
T_TOTAL = 7200      # s – 2 hours


SAMPLE_P_MAX = 120_000.0   # W – realistic installed power for sample data


def _run_physics(t, strip_speed, strip_width, strip_thickness,
                 T_strip_in, ak_os, ak_us, oven_pct, water_W,
                 strip_running, dt=DT_CSV):
    """
    Forward simulation to produce realistic 'actual' temperatures for the CSV.
    Uses SAMPLE_P_MAX (not the global P_heating_max) so temperatures stay in the
    realistic 260–290 °C operating range.
    """
    from collections import deque

    T = p.T_bath_init
    delay_steps = max(1, int(round(p.delay_heating / dt)))
    heat_buf = deque([oven_pct[0] / 100 * SAMPLE_P_MAX] * delay_steps, maxlen=delay_steps)

    temps = []
    for i in range(len(t)):
        v   = strip_speed[i] / 60.0
        bw  = strip_width[i] * 1e-3
        bt  = strip_thickness[i] * 1e-3
        mf  = p.rho_strip * v * bw * bt
        dT  = max(0.0, T - T_strip_in[i])
        Q_band = mf * p.cp_strip * dT if strip_running[i] > 0.5 else 0.0

        ak_mean = (ak_os[i] + ak_us[i]) / 2.0
        ak_pct  = 0.0 if ak_mean < 20 else min(100.0, 100 * ak_mean / 700)
        Q_ak    = float(np.interp(ak_pct, p.airknife_state_bp, p.airknife_power_tbl))
        Q_loss  = float(np.interp(T,      p.loss_temp_bp,      p.loss_power_tbl))

        P_demand = oven_pct[i] / 100.0 * SAMPLE_P_MAX
        P_h = heat_buf[0]; heat_buf.append(P_demand)

        P_net = P_h - Q_loss - Q_band - Q_ak - water_W[i]
        T    += (P_net / p.C_B) * dt
        temps.append(T)

    noise = RNG.normal(0, 0.3, len(temps))
    return np.array(temps) + noise


def main():
    n = int(T_TOTAL / DT_CSV) + 1
    t = np.linspace(0, T_TOTAL, n)

    # --- Operating scenario --------------------------------------------------
    strip_running = np.zeros(n)
    strip_running[(t >= 900) & (t < 5400)] = 1.0

    strip_speed = np.zeros(n)
    strip_speed[(t >= 900)  & (t < 1800)] = 40.0
    strip_speed[(t >= 1800) & (t < 3600)] = 60.0
    strip_speed[(t >= 3600) & (t < 4500)] = 80.0
    strip_speed[(t >= 4500) & (t < 5400)] = 60.0

    strip_width     = np.full(n, 900.0)
    strip_width[(t >= 3600) & (t < 4500)] = 1100.0

    strip_thickness = np.full(n, 0.50)
    strip_thickness[(t >= 3600) & (t < 4500)] = 0.35

    T_strip_in = np.full(n, p.T_ambient)   # strip pre-heated in annealing furnace
    T_strip_in[strip_running > 0.5] = 255.0

    ak_os = np.full(n, 0.0)    # air knife off during standstill
    ak_us = np.full(n, 0.0)
    ak_os[strip_running > 0.5] = RNG.uniform(280, 380, n)[strip_running > 0.5]
    ak_us[strip_running > 0.5] = RNG.uniform(280, 380, n)[strip_running > 0.5]

    water_W = np.zeros(n)

    T_soll_betrieb    = np.full(n, 270.0)
    T_soll_standstill = np.full(n, 265.0)

    alloy_nr   = np.zeros(n)
    strip_len  = np.zeros(n)
    strip_len[strip_running > 0.5] = 120.0

    # Oven power: realistic schedule (PLC output in %)
    # During standstill: lower power to hold 265 °C
    # During production: higher, varying with strip load
    oven_pct = np.full(n, 45.0)
    oven_pct[strip_running > 0.5] = 72.0
    # Add slow drift + random step changes
    oven_pct += np.cumsum(RNG.normal(0, 0.08, n))
    oven_pct  = np.clip(oven_pct, 30, 95)

    # Run physics to get "actual" measured temperature
    T_actual = _run_physics(t, strip_speed, strip_width, strip_thickness,
                             T_strip_in, ak_os, ak_us, oven_pct, water_W, strip_running)

    # --- Assemble DataFrame --------------------------------------------------
    timestamps = pd.date_range('2024-01-15 06:00:00', periods=n, freq=f'{DT_CSV}s')

    df = pd.DataFrame({
        'Zeitstempel':                    timestamps.strftime('%d.%m.%Y %H:%M:%S'),
        'Bandgeschwindigkeit_m_min':      np.round(strip_speed, 1),
        'Bandbreite_mm':                  np.round(strip_width, 0).astype(int),
        'Banddicke_mm':                   np.round(strip_thickness, 2),
        'Ofenleistung_Prozent':           np.round(oven_pct, 1),
        'Solltemperatur_Betrieb_C':       np.round(T_soll_betrieb, 1),
        'Solltemperatur_Stillstand_C':    np.round(T_soll_standstill, 1),
        'IstTemperatur_C':                np.round(T_actual, 2),
        'LuftmesserOS_mbar':             np.round(ak_os, 0).astype(int),
        'LuftmesserUS_mbar':             np.round(ak_us, 0).astype(int),
        'Wasserkuehlung_W':               np.round(water_W, 0).astype(int),
        'TemperaturUnterRolle_C':         np.round(T_strip_in, 1),
        'LegierungsNr':                   alloy_nr.astype(int),
        'ZustandBandlauf':                strip_running.astype(int),
        'BandlaengeImOfen_m':             np.round(strip_len, 1),
    })

    out = 'sample_data.csv'
    df.to_csv(out, sep=';', decimal=',', index=False)
    print(f"Created {out}  ({len(df)} rows, {DT_CSV}s interval, {T_TOTAL/3600:.0f} h total)")
    print(f"  Bath temp range in data: {T_actual.min():.1f} – {T_actual.max():.1f} °C")


if __name__ == '__main__':
    main()

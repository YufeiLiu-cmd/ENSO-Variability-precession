#!/usr/bin/env python3
"""02 Thermocline feedback α: Z20 → Niño3 SST (peak lag>0); compute + plot."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.paths import (
    CI_LEVEL,
    EXP_IDS,
    NINO3_MAT,
    RESULT_FIG,
    RESULT_MAT,
    SST_EQ_NC,
    Z20_NC,
    ensure_result_dirs,
)
from lib.stats import align_lag, monthly_anomaly, ols_slope_with_intercept
from lib.style import COLORS, load_nino3, plot_dual_nino3

# Identical to analyze_thermocline_feedback_regions.m
SST_LON = (210.0, 270.0)   # Niño3
Z20_LON = (140.0, 260.0)   # equatorial thermocline box
LAGS = np.arange(-6, 7)
PEAK_LAGS = "positive"     # peak among lag > 0


def thermocline_regress(T, H, lag, ci_level=CI_LEVEL):
    """T'(t+lag) = α · H'(t) when lag>=0 (H leads); matches MATLAB."""
    # MATLAB: if lag>=0: T_use=T(1+lag:end), H_use=H(1:end-lag)
    # → response T shifted forward = H leads by lag
    if lag >= 0:
        H_use, T_use = align_lag(H, T, lag)
    else:
        # lag<0: T_use=T(1:end+lag), H_use=H(1-lag:end) → T leads
        T_use, H_use = align_lag(T, H, -lag)
    return ols_slope_with_intercept(H_use, T_use, ci_level)


def main():
    ensure_result_dirs()
    if not Z20_NC.exists() or not SST_EQ_NC.exists():
        raise FileNotFoundError("Need z20_eq_lon_time_exp.nc and sst_eq_lon_time_exp.nc in data/prep")

    with Dataset(Z20_NC, "r") as nc:
        lon_z = np.asarray(nc.variables["lon"][:], float)
        exp_ids = np.asarray(nc.variables["exp_id"][:], int)
        time_range = np.asarray(nc.variables["time"][:], float)
        z20_nc = np.asarray(nc.variables["z20"][:], float)  # lon, time, exp
        target_t = float(nc.getncattr("target_temperature")) if "target_temperature" in nc.ncattrs() else 20.0

    with Dataset(SST_EQ_NC, "r") as nc:
        lon_s = np.asarray(nc.variables["lon"][:], float)
        exp_s = np.asarray(nc.variables["exp_id"][:], int)
        sst_nc = np.asarray(nc.variables["sst"][:], float)  # lon, time, exp

    if not np.array_equal(exp_ids, exp_s):
        raise RuntimeError("z20 and sst exp_id mismatch")
    # Reorder to EXP_IDS if needed
    order = [int(np.where(exp_ids == e)[0][0]) for e in EXP_IDS]
    z20_nc = z20_nc[:, :, order]
    sst_nc = sst_nc[:, :, order]
    exp_ids = np.array(EXP_IDS)

    n_exp = len(EXP_IDS)
    n_time = z20_nc.shape[1]
    n_lags = LAGS.size

    lon_z20_idx = np.where((lon_z >= Z20_LON[0]) & (lon_z <= Z20_LON[1]))[0]
    lon_sst_idx = np.where((lon_s >= SST_LON[0]) & (lon_s <= SST_LON[1]))[0]
    print(
        f"Thermocline feedback\n"
        f"  SST: 5S–5N pre-averaged, {SST_LON[0]:.0f}–{SST_LON[1]:.0f}°E ({lon_sst_idx.size} lon)\n"
        f"  Z20: 5S–5N pre-averaged, {Z20_LON[0]:.0f}–{Z20_LON[1]:.0f}°E ({lon_z20_idx.size} lon)"
    )

    sst_ts = np.full((n_time, n_exp), np.nan)
    z20_ts = np.full((n_time, n_exp), np.nan)
    sst_anom = np.full((n_time, n_exp), np.nan)
    z20_anom = np.full((n_time, n_exp), np.nan)
    alpha_lag = np.full((n_exp, n_lags), np.nan)
    alpha_ci_lag = np.full((n_exp, 2, n_lags), np.nan)
    corr_lag = np.full((n_exp, n_lags), np.nan)
    n_samples = np.full((n_exp, n_lags), np.nan)
    alpha_peak = np.full(n_exp, np.nan)
    alpha_ci = np.full((n_exp, 2), np.nan)
    lag_peak = np.full(n_exp, np.nan)

    for k in range(n_exp):
        ts_z = np.nanmean(z20_nc[lon_z20_idx, :, k], axis=0)
        ts_s = np.nanmean(sst_nc[lon_sst_idx, :, k], axis=0)
        z20_ts[:, k] = ts_z
        sst_ts[:, k] = ts_s
        z20_anom[:, k] = monthly_anomaly(ts_z)
        sst_anom[:, k] = monthly_anomaly(ts_s)

        for il, lag in enumerate(LAGS):
            a, ci, r, n = thermocline_regress(sst_anom[:, k], z20_anom[:, k], int(lag))
            alpha_lag[k, il] = a
            alpha_ci_lag[k, :, il] = ci
            corr_lag[k, il] = r
            n_samples[k, il] = n

        if PEAK_LAGS == "positive":
            peak_mask = LAGS > 0
        else:
            peak_mask = np.ones_like(LAGS, dtype=bool)
        peak_idx = np.where(peak_mask)[0]
        i_max = int(np.nanargmax(alpha_lag[k, peak_idx]))
        j = peak_idx[i_max]
        alpha_peak[k] = alpha_lag[k, j]
        lag_peak[k] = LAGS[j]
        alpha_ci[k, :] = alpha_ci_lag[k, :, j]
        print(f"EXP {exp_ids[k]:03d} | alpha_peak={alpha_peak[k]:+.4f} K/m (lag={int(lag_peak[k]):+d})")

    mat_path = RESULT_MAT / "thermocline_feedback_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": exp_ids.reshape(-1, 1).astype(float),
            "timeRange": time_range.reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "method": "SST_on_Z20_leadlag_eq_latmean",
            "peakLags": PEAK_LAGS,
            "sstLat": np.array([-5.0, 5.0]),
            "sstLon": np.array(SST_LON),
            "z20Lat": np.array([-5.0, 5.0]),
            "z20Lon": np.array(Z20_LON),
            "targetT": target_t,
            "lags": LAGS.reshape(1, -1),
            "lonSstIdx": lon_sst_idx.reshape(-1, 1) + 1,
            "lonZ20Idx": lon_z20_idx.reshape(-1, 1) + 1,
            "lonS": lon_s.reshape(-1, 1),
            "lonZ": lon_z.reshape(-1, 1),
            "sst_ts": sst_ts,
            "z20_ts": z20_ts,
            "sst_anom": sst_anom,
            "z20_anom": z20_anom,
            "alpha_lag": alpha_lag,
            "alpha_ci_lag": alpha_ci_lag,
            "corr_lag": corr_lag,
            "nSamples": n_samples,
            "alpha_peak": alpha_peak.reshape(-1, 1),
            "alpha_ci": alpha_ci,
            "lag_peak": lag_peak.reshape(-1, 1),
        },
    )
    print(f"Saved {mat_path}")

    err_lo = alpha_peak - alpha_ci[:, 0]
    err_hi = alpha_ci[:, 1] - alpha_peak

    enso, n_lo, n_hi = load_nino3(NINO3_MAT)
    plot_dual_nino3(
        EXP_IDS, enso, n_lo, n_hi, alpha_peak, err_lo, err_hi,
        color_right=COLORS["alpha"],
        ylabel_right=r"Thermocline feedback $\alpha$ (K m$^{-1}$)",
        label_right=r"Thermocline feedback ($\alpha$)",
        stem="03_nino3_thermocline_feedback_dual_axis",
        out_dir=RESULT_FIG,
        mean_fmt_right=".3e",
    )
    print(f"mean α = {np.nanmean(alpha_peak):.6e} K/m")
    print("Thermocline feedback done.")


if __name__ == "__main__":
    main()

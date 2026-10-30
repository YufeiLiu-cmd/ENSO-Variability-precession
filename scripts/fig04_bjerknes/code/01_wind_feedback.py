#!/usr/bin/env python3
"""01 Wind feedback μ: SST(Niño3) → τx(West Pac); compute + plot."""
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
    N_TIME,
    NINO3_MAT,
    RESULT_FIG,
    RESULT_MAT,
    TSW_DIR,
    USTRW_DIR,
    ensure_result_dirs,
)
from lib.stats import align_lag, area_weights, monthly_anomaly, ols_slope_through_origin, regional_mean_llt
from lib.style import COLORS, load_nino3, plot_dual_nino3

# Method A regions (0–360°E) — identical to analyze_wind_feedback.m
LAT_EQ = (-5.0, 5.0)
SST_LON = (210.0, 270.0)   # Niño3: 150°W–90°W
TX_LON = (140.0, 180.0)    # West Pacific wind-response box
LAGS = (0, 1)


def wind_feedback_regress(T, TX, lag, ci_level=CI_LEVEL):
    """τx'(t+lag) = μ · T'(t)  (MATLAB windFeedbackRegress)."""
    T_use, TX_use = align_lag(T, TX, lag)
    return ols_slope_through_origin(T_use, TX_use, ci_level)


def main():
    ensure_result_dirs()
    n_exp = len(EXP_IDS)
    n_lags = len(LAGS)

    first_tsw = TSW_DIR / f"{EXP_IDS[0]:03d}tsw.nc"
    first_tx = USTRW_DIR / f"{EXP_IDS[0]:03d}.nc"
    with Dataset(first_tsw, "r") as nc:
        lon = np.asarray(nc.variables["lon"][:], float)
        lat = np.asarray(nc.variables["lat"][:], float)

    lat_sst_idx = np.where((lat >= LAT_EQ[0]) & (lat <= LAT_EQ[1]))[0]
    lon_sst_idx = np.where((lon >= SST_LON[0]) & (lon <= SST_LON[1]))[0]
    lat_tx_idx = lat_sst_idx
    lon_tx_idx = np.where((lon >= TX_LON[0]) & (lon <= TX_LON[1]))[0]
    if lat_sst_idx.size == 0 or lon_sst_idx.size == 0 or lon_tx_idx.size == 0:
        raise RuntimeError("Empty index for SST or τx region")

    w_sst = area_weights(lat[lat_sst_idx], lon_sst_idx.size)
    w_tx = area_weights(lat[lat_tx_idx], lon_tx_idx.size)

    print(
        f"Wind feedback regions\n"
        f"  SST Niño3: {LAT_EQ[0]}–{LAT_EQ[1]}°N, {SST_LON[0]}–{SST_LON[1]}°E "
        f"| lon {lon_sst_idx.size}, lat {lat_sst_idx.size}\n"
        f"  τx WestPac: {LAT_EQ[0]}–{LAT_EQ[1]}°N, {TX_LON[0]}–{TX_LON[1]}°E "
        f"| lon {lon_tx_idx.size}, lat {lat_tx_idx.size}\n"
        f"  time 1:{N_TIME} | lags {LAGS}"
    )

    tsw_ts = np.full((N_TIME, n_exp), np.nan)
    ustrw_ts = np.full((N_TIME, n_exp), np.nan)
    tsw_anom = np.full((N_TIME, n_exp), np.nan)
    ustrw_anom = np.full((N_TIME, n_exp), np.nan)
    mu = np.full((n_exp, n_lags), np.nan)
    mu_ci = np.full((n_exp, 2, n_lags), np.nan)
    corr_rt = np.full((n_exp, n_lags), np.nan)
    n_samples = np.full((n_exp, n_lags), np.nan)

    for k, eid in enumerate(EXP_IDS):
        f_tsw = TSW_DIR / f"{eid:03d}tsw.nc"
        f_tx = USTRW_DIR / f"{eid:03d}.nc"
        print(f"EXP {eid:03d} | reading tsw & ustrw ...", flush=True)
        with Dataset(f_tsw, "r") as nc:
            raw_t = np.asarray(nc.variables["tsw"][:N_TIME], dtype=float)  # time,lat,lon
        with Dataset(f_tx, "r") as nc:
            raw_u = np.asarray(nc.variables["ustrw"][:N_TIME], dtype=float)

        ts_t = regional_mean_llt(raw_t, lat_sst_idx, lon_sst_idx, w_sst)
        ts_u = regional_mean_llt(raw_u, lat_tx_idx, lon_tx_idx, w_tx)
        tsw_ts[:, k] = ts_t
        ustrw_ts[:, k] = ts_u
        tsw_anom[:, k] = monthly_anomaly(ts_t)
        ustrw_anom[:, k] = monthly_anomaly(ts_u)

        for il, lag in enumerate(LAGS):
            m, ci, r, n = wind_feedback_regress(tsw_anom[:, k], ustrw_anom[:, k], lag)
            mu[k, il] = m
            mu_ci[k, :, il] = ci
            corr_rt[k, il] = r
            n_samples[k, il] = n
        print(f"  μ(lag0)={mu[k,0]:+.4e}  μ(lag1)={mu[k,1]:+.4e}  r(lag0)={corr_rt[k,0]:+.3f}")

    mat_path = RESULT_MAT / "wind_feedback_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": np.array(EXP_IDS, float).reshape(-1, 1),
            "timeRange": np.arange(1, N_TIME + 1).reshape(1, -1),
            "ciLevel": CI_LEVEL,
            "method": "A",
            "latEq": np.array(LAT_EQ),
            "sstLon": np.array(SST_LON),
            "txLon": np.array(TX_LON),
            "latSstIdx": lat_sst_idx.reshape(-1, 1) + 1,  # 1-based like MATLAB
            "lonSstIdx": lon_sst_idx.reshape(-1, 1) + 1,
            "latTxIdx": lat_tx_idx.reshape(-1, 1) + 1,
            "lonTxIdx": lon_tx_idx.reshape(-1, 1) + 1,
            "lon": lon.reshape(-1, 1),
            "lat": lat.reshape(-1, 1),
            "lags": np.array(LAGS).reshape(1, -1),
            "lagLabels": np.array(
                ["lag 0: SST(t) -> tau_x(t)", "lag 1: SST(t) -> tau_x(t+1)"], dtype=object
            ),
            "tsw_ts": tsw_ts,
            "ustrw_ts": ustrw_ts,
            "tsw_anom": tsw_anom,
            "ustrw_anom": ustrw_anom,
            "mu": mu,
            "mu_ci": mu_ci,
            "corr_rt": corr_rt,
            "nSamples": n_samples,
        },
    )
    print(f"Saved {mat_path}")

    # Primary diagnostic: lag=0 (as in dual-axis paper figure)
    lag0 = 0
    mu_plot = mu[:, lag0]
    ci_plot = mu_ci[:, :, lag0]
    err_lo = mu_plot - ci_plot[:, 0]
    err_hi = ci_plot[:, 1] - mu_plot

    enso, n_lo, n_hi = load_nino3(NINO3_MAT)
    plot_dual_nino3(
        EXP_IDS, enso, n_lo, n_hi, mu_plot, err_lo, err_hi,
        color_right=COLORS["mu"],
        ylabel_right=r"Wind feedback $\mu$ (Pa K$^{-1}$)",
        label_right=r"Wind feedback ($\mu$, lag 0)",
        stem="01_nino3_wind_feedback_dual_axis",
        out_dir=RESULT_FIG,
        mean_fmt_right=".3e",
        legend_ncol=2,
    )
    print(f"mean μ(lag0) = {np.nanmean(mu_plot):.6e} Pa/K")
    print("Wind feedback done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""03 Ocean response β: West-Pac τx → Z20; compute + plot (Bjerknes middle link)."""
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
    USTRW_DIR,
    Z20_NC,
    ensure_result_dirs,
)
from lib.stats import (
    align_lag,
    area_weights,
    monthly_anomaly,
    ols_slope_with_intercept,
    regional_mean_llt,
)
from lib.style import COLORS, load_nino3, plot_dual_nino3

# τx: West Pac (same as μ); Z20: eastern Pacific Niño3 longitudes
# (Deppenmeier et al. 2016; Lübbecke & McPhaden — eastern thermocline response)
LAT_EQ = (-5.0, 5.0)
TX_LON = (140.0, 180.0)
Z20_LON = (210.0, 270.0)  # Niño3 150°W–90°W — eastern Pacific Z20
LAGS = np.arange(0, 7)          # wind leads thermocline
PEAK_LAGS = "positive"          # peak among lag >= 0 (all scanned)


def ocean_response_regress(TX, H, lag, ci_level=CI_LEVEL):
    """H'(t+lag) = β · τx'(t)  — wind leads Z20."""
    TX_use, H_use = align_lag(TX, H, lag)
    return ols_slope_with_intercept(TX_use, H_use, ci_level)


def main():
    ensure_result_dirs()
    # Prefer reusing τx anomalies from wind-feedback mat if present
    wind_mat = RESULT_MAT / "wind_feedback_results.mat"

    n_exp = len(EXP_IDS)
    n_lags = LAGS.size

    if wind_mat.exists():
        W = sio.loadmat(wind_mat)
        ustrw_anom = np.asarray(W["ustrw_anom"], float)
        print(f"Loaded τx anomalies from {wind_mat.name}")
    else:
        print("Computing τx anomalies from ustrw nc ...")
        first = USTRW_DIR / f"{EXP_IDS[0]:03d}.nc"
        with Dataset(first, "r") as nc:
            lon = np.asarray(nc.variables["lon"][:], float)
            lat = np.asarray(nc.variables["lat"][:], float)
        lat_idx = np.where((lat >= LAT_EQ[0]) & (lat <= LAT_EQ[1]))[0]
        lon_idx = np.where((lon >= TX_LON[0]) & (lon <= TX_LON[1]))[0]
        w = area_weights(lat[lat_idx], lon_idx.size)
        ustrw_anom = np.full((N_TIME, n_exp), np.nan)
        for k, eid in enumerate(EXP_IDS):
            with Dataset(USTRW_DIR / f"{eid:03d}.nc", "r") as nc:
                raw = np.asarray(nc.variables["ustrw"][:N_TIME], float)
            ts = regional_mean_llt(raw, lat_idx, lon_idx, w)
            ustrw_anom[:, k] = monthly_anomaly(ts)

    # Always recompute Z20 over Niño3 longitudes from prep nc
    # (do NOT reuse thermocline mat: that Z20 is 140–260°E)
    print(f"Computing Z20 anomalies from prep nc over {Z20_LON[0]:.0f}–{Z20_LON[1]:.0f}°E ...")
    with Dataset(Z20_NC, "r") as nc:
        lon_z = np.asarray(nc.variables["lon"][:], float)
        exp_ids = np.asarray(nc.variables["exp_id"][:], int)
        z20_nc = np.asarray(nc.variables["z20"][:], float)
    order = [int(np.where(exp_ids == e)[0][0]) for e in EXP_IDS]
    z20_nc = z20_nc[:, :, order]
    lon_idx = np.where((lon_z >= Z20_LON[0]) & (lon_z <= Z20_LON[1]))[0]
    if lon_idx.size == 0:
        raise RuntimeError(f"No Z20 lon points in {Z20_LON}")
    n_time_z = z20_nc.shape[1]
    z20_anom = np.full((n_time_z, n_exp), np.nan)
    for k in range(n_exp):
        sub = z20_nc[lon_idx, :, k]
        # skip all-NaN lon columns (land / missing)
        valid_lon = np.isfinite(sub).any(axis=1)
        if not np.any(valid_lon):
            raise RuntimeError(f"EXP {EXP_IDS[k]:03d}: no valid Z20 in Niño3 box")
        sub_v = sub[valid_lon, :]
        cnt = np.sum(np.isfinite(sub_v), axis=0).astype(float)
        s = np.where(np.isfinite(sub_v), sub_v, 0.0).sum(axis=0)
        ts = np.divide(s, cnt, out=np.full(cnt.shape, np.nan), where=cnt > 0)
        z20_anom[:, k] = monthly_anomaly(ts)
    print(f"  Z20 Niño3 box: {lon_idx.size} lon points "
          f"({int(np.isfinite(z20_nc[lon_idx, :, 0]).any(axis=1).sum())} with data)")

    # Align lengths
    n_time = min(ustrw_anom.shape[0], z20_anom.shape[0])
    ustrw_anom = ustrw_anom[:n_time]
    z20_anom = z20_anom[:n_time]

    beta_lag = np.full((n_exp, n_lags), np.nan)
    beta_ci_lag = np.full((n_exp, 2, n_lags), np.nan)
    corr_lag = np.full((n_exp, n_lags), np.nan)
    n_samples = np.full((n_exp, n_lags), np.nan)
    beta_peak = np.full(n_exp, np.nan)
    beta_ci = np.full((n_exp, 2), np.nan)
    lag_peak = np.full(n_exp, np.nan)

    print(
        f"Ocean response β\n"
        f"  τx: 5S–5N, {TX_LON[0]:.0f}–{TX_LON[1]:.0f}°E\n"
        f"  Z20: 5S–5N, {Z20_LON[0]:.0f}–{Z20_LON[1]:.0f}°E\n"
        f"  lags {LAGS.tolist()} (wind leads); peak |β| among lags"
    )

    for k, eid in enumerate(EXP_IDS):
        for il, lag in enumerate(LAGS):
            b, ci, r, n = ocean_response_regress(ustrw_anom[:, k], z20_anom[:, k], int(lag))
            beta_lag[k, il] = b
            beta_ci_lag[k, :, il] = ci
            corr_lag[k, il] = r
            n_samples[k, il] = n
        # Peak: maximum β (westerly → deeper Z20 expected positive)
        i_max = int(np.nanargmax(beta_lag[k, :]))
        beta_peak[k] = beta_lag[k, i_max]
        lag_peak[k] = LAGS[i_max]
        beta_ci[k, :] = beta_ci_lag[k, :, i_max]
        print(f"EXP {eid:03d} | beta_peak={beta_peak[k]:+.4e} m/Pa (lag={int(lag_peak[k]):+d})")

    mat_path = RESULT_MAT / "ocean_response_beta_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": np.array(EXP_IDS, float).reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "method": "Z20_on_tau_x_leadlag_Bjerknes_middle_link",
            "latEq": np.array(LAT_EQ),
            "txLon": np.array(TX_LON),
            "z20Lon": np.array(Z20_LON),
            "lags": LAGS.reshape(1, -1),
            "ustrw_anom": ustrw_anom,
            "z20_anom": z20_anom,
            "beta_lag": beta_lag,
            "beta_ci_lag": beta_ci_lag,
            "corr_lag": corr_lag,
            "nSamples": n_samples,
            "beta_peak": beta_peak.reshape(-1, 1),
            "beta_ci": beta_ci,
            "lag_peak": lag_peak.reshape(-1, 1),
            "refs": (
                "Deppenmeier et al. 2016 Clim Dyn; "
                "Luebbecke & McPhaden 2013/2017; "
                "Bjerknes 1969 three-element loop"
            ),
        },
    )
    print(f"Saved {mat_path}")

    err_lo = beta_peak - beta_ci[:, 0]
    err_hi = beta_ci[:, 1] - beta_peak

    enso, n_lo, n_hi = load_nino3(NINO3_MAT)
    plot_dual_nino3(
        EXP_IDS, enso, n_lo, n_hi, beta_peak, err_lo, err_hi,
        color_right=COLORS["beta"],
        ylabel_right=r"Ocean response $\beta$ (m Pa$^{-1}$)",
        label_right=r"Ocean response ($\beta$)",
        stem="02_nino3_ocean_response_beta_dual_axis",
        out_dir=RESULT_FIG,
        mean_fmt_right=".3e",
        legend_ncol=2,
    )
    print(f"mean β = {np.nanmean(beta_peak):.6e} m/Pa")
    print("Ocean response β done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""05 Upper-ocean stratification (Cai et al. 2018 style): mean(T0-100m)-T120m.

Region: 5S–5N, 140E–150W. Ocean lon -180..180 converted to 0..360.
Means and 95% CIs both from monthly Scheme 1 (no annual substitute;
no EXP120 special case).
"""
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
    NINO3_DIR,
    RESULT_FIG,
    RESULT_MAT,
    THETAO_DIR,
    ensure_result_dirs,
)
from lib.preprocess import lon_to_360, scheme1_mean_and_ci
from lib.style import COLORS, dual_axis_errorbar, save_figure, single_errorbar

LAT_MIN, LAT_MAX = -5.0, 5.0
LON_MIN, LON_MAX = 140.0, 210.0
FILL = -9999.0
DEPTH_INTERP = np.arange(0, 341, 1)  # 0..340 m
N_TIME_DEFAULT = 1200


def _interp_many(z: np.ndarray, T: np.ndarray, zf: np.ndarray) -> np.ndarray:
    """Linear interpolate profiles T(n, nz) onto zf."""
    n = T.shape[0]
    out = np.full((n, zf.size), np.nan)
    for i in range(n):
        prof = T[i]
        valid = np.isfinite(prof)
        if valid.sum() < 2:
            continue
        out[i] = np.interp(zf, z[valid], prof[valid])
    return out


def monthly_stratification_fast(fpath: Path, n_time: int = N_TIME_DEFAULT) -> np.ndarray:
    with Dataset(fpath, "r") as ds:
        lon = np.asarray(ds.variables["lon"][:], dtype=float)
        lat = np.asarray(ds.variables["lat"][:], dtype=float)
        depth = np.asarray(ds.variables["depth"][:], dtype=float)
        lon360 = lon_to_360(lon)
        lat_idx = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
        lon_idx = np.where((lon360 >= LON_MIN) & (lon360 <= LON_MAX))[0]
        lat0, lat1 = int(lat_idx.min()), int(lat_idx.max()) + 1
        block = np.asarray(ds.variables["thetao"][0:n_time, :, lat0:lat1, :], dtype=float)
        lat_sub = lat[lat0:lat1]
        lat_mask = (lat_sub >= LAT_MIN) & (lat_sub <= LAT_MAX)
        block = block[:, :, lat_mask, :][:, :, :, lon_idx]
        block[block <= FILL] = np.nan

    T = np.nanmean(block, axis=2)  # time, depth, lon
    nt, nz, nlon = T.shape
    Tp = np.transpose(T, (0, 2, 1)).reshape(nt * nlon, nz)
    Ti = _interp_many(depth, Tp, DEPTH_INTERP).reshape(nt, nlon, DEPTH_INTERP.size)
    t_upper = np.nanmean(Ti[:, :, DEPTH_INTERP <= 100], axis=(1, 2))
    i120 = int(np.where(DEPTH_INTERP == 120)[0][0])
    t120 = np.nanmean(Ti[:, :, i120], axis=1)
    return (t_upper - t120).astype(float)


def main():
    ensure_result_dirs()
    n_exp = len(EXP_IDS)
    ts_all = np.full((1200, n_exp), np.nan)
    means = np.full(n_exp, np.nan)
    half_lo = np.full(n_exp, np.nan)
    half_hi = np.full(n_exp, np.nan)
    half = np.full(n_exp, np.nan)
    n_samp = np.full(n_exp, np.nan)
    n_time_used = np.full(n_exp, N_TIME_DEFAULT, dtype=int)

    for k, eid in enumerate(EXP_IDS):
        fpath = THETAO_DIR / f"thetao_{eid:03d}_300m.nc"
        if not fpath.exists():
            raise FileNotFoundError(fpath)
        print(f"[STRAT] {fpath.name} ...", flush=True)
        ts = monthly_stratification_fast(fpath, n_time=N_TIME_DEFAULT)
        ts_all[:, k] = ts
        mu, lo, hi, diag = scheme1_mean_and_ci(ts, CI_LEVEL)
        means[k] = mu
        half_lo[k] = mu - lo
        half_hi[k] = hi - mu
        half[k] = diag["half_width"]
        n_samp[k] = diag["n"]
        n_finite = int(np.isfinite(ts).sum())
        n_time_used[k] = n_finite
        print(
            f"  EXP {eid:03d} mean={mu:.4f} C  CI half={half[k]:.4f}  "
            f"finite_months={n_finite}",
            flush=True,
        )

    ci = np.column_stack([means - half_lo, means + half_hi])
    sio.savemat(
        RESULT_MAT / "stratification_precession_results.mat",
        {
            "expIDs": np.array(EXP_IDS, float).reshape(-1, 1),
            "stratification_month": ts_all.T,
            "stratification": means.reshape(-1, 1),
            "CI95": half.reshape(-1, 1),
            "strat_ci": ci,
            "ci_half_lo": half_lo.reshape(-1, 1),
            "ci_half_hi": half_hi.reshape(-1, 1),
            "nSamples": n_samp.reshape(-1, 1),
            "nTimeUsed": n_time_used.reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "latMin": LAT_MIN,
            "latMax": LAT_MAX,
            "lonMin": LON_MIN,
            "lonMax": LON_MAX,
            "depth_interp": DEPTH_INTERP.reshape(-1, 1),
            "definition": "mean(T0-100m) - T120m",
            "preprocess": "monthly Scheme1 mean (3mo runmean); 95% CI from annual means of smoothed (n~100)",
            "lon_note": "thetao lon -180..180 -> 0..360 before region select",
            "thetao_source": str(THETAO_DIR.resolve()),
        },
    )

    fig, _ = single_errorbar(
        EXP_IDS, means, half_lo, half_hi, color=COLORS["strat"],
        ylabel=r"Upper-ocean stratification ($^\circ$C)", label="Stratification",
        mean_val=float(np.nanmean(means)),
    )
    save_figure(fig, "05_stratification_vs_precession", RESULT_FIG)

    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    enso = nino["enso_amp"].ravel()
    nci = nino["chi2_ci_std"]
    fig, _, _ = dual_axis_errorbar(
        EXP_IDS, enso, enso - nci[:, 0], nci[:, 1] - enso, means, half_lo, half_hi,
        color_left=COLORS["nino3"], color_right=COLORS["strat"],
        ylabel_left="Nino3 amplitude (K)", ylabel_right=r"Stratification ($^\circ$C)",
        label_left="Nino3 amplitude", label_right="Stratification",
        mean_left=float(np.nanmean(enso)), mean_right=float(np.nanmean(means)),
    )
    save_figure(fig, "05_nino3_stratification_dual_axis", RESULT_FIG)
    print("Stratification done.")


if __name__ == "__main__":
    main()

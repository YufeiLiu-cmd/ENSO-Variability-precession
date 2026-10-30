#!/usr/bin/env python3
"""02 USTRW area-mean zonal wind stress: compute + plot (Scheme 1)."""
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
    USTRW_DIR,
    ensure_result_dirs,
)
from lib.preprocess import scheme1_mean_and_ci
from lib.style import COLORS, dual_axis_errorbar, save_figure, single_errorbar


LAT_RANGE = (-5.0, 5.0)
LON_RANGE = (140.0, 210.0)  # 140E–150W in 0–360 (atmosphere)


def area_mean_ustrw_fast(fpath: Path) -> np.ndarray:
    with Dataset(fpath, "r") as ds:
        lon = np.asarray(ds.variables["lon"][:], dtype=float)
        lat = np.asarray(ds.variables["lat"][:], dtype=float)
        lat_idx = np.where((lat >= LAT_RANGE[0]) & (lat <= LAT_RANGE[1]))[0]
        lon_idx = np.where((lon >= LON_RANGE[0]) & (lon <= LON_RANGE[1]))[0]
        la0, la1 = int(lat_idx.min()), int(lat_idx.max()) + 1
        lo0, lo1 = int(lon_idx.min()), int(lon_idx.max()) + 1
        u = np.asarray(ds.variables["ustrw"][0:1200, la0:la1, lo0:lo1], dtype=float)
        lat_sub = lat[la0:la1]
        lon_sub = lon[lo0:lo1]
        lat_mask = (lat_sub >= LAT_RANGE[0]) & (lat_sub <= LAT_RANGE[1])
        lon_mask = (lon_sub >= LON_RANGE[0]) & (lon_sub <= LON_RANGE[1])
        u = u[:, lat_mask][:, :, lon_mask]
        wlat = np.cos(np.deg2rad(lat_sub[lat_mask]))
        w2d = np.ones((wlat.size, int(lon_mask.sum())), dtype=float) * wlat[:, None]
        w2d /= w2d.sum()
        ts = np.nansum(u * w2d[None, :, :], axis=(1, 2))
    return ts.astype(float)


def main():
    ensure_result_dirs()
    n_exp = len(EXP_IDS)
    ts_all = np.full((1200, n_exp), np.nan)
    means = np.full(n_exp, np.nan)
    ci = np.full((n_exp, 2), np.nan)
    n_samp = np.full(n_exp, np.nan)

    for k, eid in enumerate(EXP_IDS):
        fpath = USTRW_DIR / f"{eid:03d}.nc"
        print(f"Reading {fpath.name} ...", flush=True)
        ts = area_mean_ustrw_fast(fpath)
        ts_all[:, k] = ts
        mu, lo, hi, diag = scheme1_mean_and_ci(ts, CI_LEVEL)
        means[k] = mu
        ci[k] = (lo, hi)
        n_samp[k] = diag["n"]
        print(f"  EXP {eid:03d} mean={mu:.4e} CI=[{lo:.4e},{hi:.4e}] n={diag['n']}", flush=True)

    exp = np.array(EXP_IDS, float)
    mat_path = RESULT_MAT / "ustrw_precession_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": exp.reshape(-1, 1),
            "ustrw_ts": ts_all,
            "ustrw_mean": means.reshape(-1, 1),
            "ustrw_ci": ci,
            "nSamples": n_samp.reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "latRange": np.array(LAT_RANGE),
            "lonRange": np.array(LON_RANGE),
            "preprocess": "Scheme1 mean (3mo runmean); 95% CI from annual means of smoothed (n~100)",
        },
    )
    print(f"Saved {mat_path}")

    err_lo = means - ci[:, 0]
    err_hi = ci[:, 1] - means
    fig, _ = single_errorbar(
        EXP_IDS,
        means,
        err_lo,
        err_hi,
        color=COLORS["ustrw"],
        ylabel=r"Area-mean zonal wind stress (Pa)",
        label=r"Area-mean $\tau_x$",
        mean_val=float(np.nanmean(means)),
    )
    save_figure(fig, "02_ustrw_vs_precession", RESULT_FIG)

    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    enso = nino["enso_amp"].ravel()
    nci = nino["chi2_ci_std"]
    fig, _, _ = dual_axis_errorbar(
        EXP_IDS,
        enso,
        enso - nci[:, 0],
        nci[:, 1] - enso,
        means,
        err_lo,
        err_hi,
        color_left=COLORS["nino3"],
        color_right=COLORS["ustrw"],
        ylabel_left="Nino3 amplitude (K)",
        ylabel_right="Area-mean zonal wind stress (Pa)",
        label_left="Nino3 amplitude",
        label_right=r"Area-mean $\tau_x$",
        mean_left=float(np.nanmean(enso)),
        mean_right=float(np.nanmean(means)),
    )
    save_figure(fig, "02_nino3_ustrw_dual_axis", RESULT_FIG)
    print("USTRW done.")


if __name__ == "__main__":
    main()

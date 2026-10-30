#!/usr/bin/env python3
"""04 MLD: 0.5 K threshold over 1200 months, Scheme-1 mean + t-CI, plot.

MLD crossing is evaluated on piecewise-linear native depth profiles
(equivalent to 1 m fine-grid search for monotonic segments).
Ocean lon (-180..180) is converted to 0..360 before selecting 140E–150W.

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

LAT_MIN, LAT_MAX = -2.5, 2.5
LON_MIN, LON_MAX = 140.0, 210.0
DT = 0.5
FILL = -9999.0
N_TIME_DEFAULT = 1200


def calc_mld_vectorized(z: np.ndarray, T: np.ndarray, dT: float) -> np.ndarray:
    """T shape (n, nz); z shape (nz,). Returns (n,)."""
    z = np.asarray(z, dtype=float)
    T = np.asarray(T, dtype=float)
    n, nz = T.shape
    out = np.full(n, np.nan)

    tref = T[:, 0]
    good = np.isfinite(tref)
    crit = tref - dT
    Twork = np.where(np.isfinite(T), T, np.inf)

    below = Twork < crit[:, None]
    any_below = below.any(axis=1) & good
    first = np.argmax(below, axis=1)

    finite = np.isfinite(T)
    rev = finite[:, ::-1]
    last_from_end = np.argmax(rev, axis=1)
    last = nz - 1 - last_from_end
    no_cross = good & ~any_below
    out[no_cross] = z[last[no_cross]]

    at0 = any_below & (first == 0)
    out[at0] = 0.0

    mid = any_below & (first > 0)
    if np.any(mid):
        i2 = first[mid]
        i1 = i2 - 1
        rows = np.where(mid)[0]
        z1 = z[i1]
        z2 = z[i2]
        t1 = T[rows, i1]
        t2 = T[rows, i2]
        c = crit[rows]
        denom = t2 - t1
        mld = np.where(
            np.abs(denom) < 1e-12,
            z2,
            z1 + (c - t1) * (z2 - z1) / denom,
        )
        out[rows] = mld
    return out


def monthly_area_mean_mld(fpath: Path, n_time: int = N_TIME_DEFAULT) -> np.ndarray:
    with Dataset(fpath, "r") as ds:
        lon = np.asarray(ds.variables["lon"][:], dtype=float)
        lat = np.asarray(ds.variables["lat"][:], dtype=float)
        depth = np.asarray(ds.variables["depth"][:], dtype=float)
        lon360 = lon_to_360(lon)
        lat_idx = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
        lon_idx = np.where((lon360 >= LON_MIN) & (lon360 <= LON_MAX))[0]
        if lat_idx.size == 0 or lon_idx.size == 0:
            raise RuntimeError(f"Empty MLD region in {fpath}")

        lat0, lat1 = int(lat_idx.min()), int(lat_idx.max()) + 1
        block = np.asarray(ds.variables["thetao"][0:n_time, :, lat0:lat1, :], dtype=float)
        lat_sub = lat[lat0:lat1]
        lat_mask = (lat_sub >= LAT_MIN) & (lat_sub <= LAT_MAX)
        block = block[:, :, lat_mask, :]
        block = block[:, :, :, lon_idx]
        block[block <= FILL] = np.nan

    nt, nz, nlat, nlon = block.shape
    T = np.transpose(block, (0, 2, 3, 1)).reshape(nt * nlat * nlon, nz)
    mld = calc_mld_vectorized(depth, T, DT).reshape(nt, nlat, nlon)
    return np.nanmean(mld, axis=(1, 2)).astype(float)


def main():
    ensure_result_dirs()
    n_exp = len(EXP_IDS)
    ts_all = np.full((1200, n_exp), np.nan)
    means = np.full(n_exp, np.nan)
    half_lo = np.full(n_exp, np.nan)
    half_hi = np.full(n_exp, np.nan)
    n_samp = np.full(n_exp, np.nan)
    n_time_used = np.full(n_exp, N_TIME_DEFAULT, dtype=int)

    for k, eid in enumerate(EXP_IDS):
        fpath = THETAO_DIR / f"thetao_{eid:03d}_300m.nc"
        if not fpath.exists():
            raise FileNotFoundError(fpath)
        print(f"[MLD] {fpath.name} ...", flush=True)
        ts = monthly_area_mean_mld(fpath, n_time=N_TIME_DEFAULT)
        ts_all[:, k] = ts
        mu, lo, hi, diag = scheme1_mean_and_ci(ts, CI_LEVEL)
        means[k] = mu
        half_lo[k] = mu - lo
        half_hi[k] = hi - mu
        n_samp[k] = diag["n"]
        n_finite = int(np.isfinite(ts).sum())
        n_time_used[k] = n_finite
        print(
            f"  EXP {eid:03d} mean={mu:.3f} m  CI half=[{half_lo[k]:.3f},{half_hi[k]:.3f}]  "
            f"finite_months={n_finite}",
            flush=True,
        )

    ci = np.column_stack([means - half_lo, means + half_hi])
    sio.savemat(
        RESULT_MAT / "mld_precession_results.mat",
        {
            "expIDs": np.array(EXP_IDS, float).reshape(-1, 1),
            "mld_ts": ts_all,
            "mld_mean": means.reshape(-1, 1),
            "mld_ci_mean": ci,
            "ci_half_lo": half_lo.reshape(-1, 1),
            "ci_half_hi": half_hi.reshape(-1, 1),
            "nSamples": n_samp.reshape(-1, 1),
            "nTimeUsed": n_time_used.reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "latMin": LAT_MIN,
            "latMax": LAT_MAX,
            "lonMin": LON_MIN,
            "lonMax": LON_MAX,
            "dT_threshold": DT,
            "preprocess": "monthly Scheme1 mean (3mo runmean); 95% CI from annual means of smoothed (n~100)",
            "lon_note": "thetao lon -180..180 -> 0..360 before region select",
            "thetao_source": str(THETAO_DIR.resolve()),
        },
    )

    fig, _ = single_errorbar(
        EXP_IDS, means, half_lo, half_hi, color=COLORS["mld"],
        ylabel="Mixed layer depth (m)", label="Mixed layer depth",
        mean_val=float(np.nanmean(means)),
    )
    save_figure(fig, "04_mld_vs_precession", RESULT_FIG)

    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    enso = nino["enso_amp"].ravel()
    nci = nino["chi2_ci_std"]
    fig, _, _ = dual_axis_errorbar(
        EXP_IDS, enso, enso - nci[:, 0], nci[:, 1] - enso, means, half_lo, half_hi,
        color_left=COLORS["nino3"], color_right=COLORS["mld"],
        ylabel_left="Nino3 amplitude (K)", ylabel_right="Mixed layer depth (m)",
        label_left="Nino3 amplitude", label_right="Mixed layer depth",
        mean_left=float(np.nanmean(enso)), mean_right=float(np.nanmean(means)),
    )
    save_figure(fig, "04_nino3_mld_dual_axis", RESULT_FIG)
    print("MLD done.")


if __name__ == "__main__":
    main()

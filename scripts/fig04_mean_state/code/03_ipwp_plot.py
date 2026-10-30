#!/usr/bin/env python3
"""03 IPWP mean position: plot only (precomputed barycenter + error24)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.paths import EXP_IDS, IPWP_DIR, NINO3_DIR, RESULT_FIG, RESULT_MAT, ensure_result_dirs
from lib.style import COLORS, dual_axis_errorbar, save_figure, single_errorbar


def main():
    ensure_result_dirs()
    ipwp = sio.loadmat(IPWP_DIR / "exp24aprczhixin.mat")
    bary = ipwp["barycenter"].ravel().astype(float)
    err = ipwp["error24"].ravel().astype(float)

    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    enso = nino["enso_amp"].ravel()
    nci = nino["chi2_ci_std"]

    exp = np.array(EXP_IDS, float)
    sio.savemat(
        RESULT_MAT / "ipwp_precession_results.mat",
        {
            "expIDs": exp.reshape(-1, 1),
            "barycenter": bary.reshape(-1, 1),
            "error24": err.reshape(-1, 1),
            "note": "Plotted from exp24aprczhixin.mat; no recomputation",
        },
    )

    fig, _ = single_errorbar(
        EXP_IDS,
        bary,
        err,
        err,
        color=COLORS["ipwp"],
        ylabel=r"IPWP mean position ($^\circ$E)",
        label="IPWP mean position",
        mean_val=float(np.nanmean(bary)),
    )
    save_figure(fig, "03_ipwp_vs_precession", RESULT_FIG)

    fig, _, _ = dual_axis_errorbar(
        EXP_IDS,
        enso,
        enso - nci[:, 0],
        nci[:, 1] - enso,
        bary,
        err,
        err,
        color_left=COLORS["nino3"],
        color_right=COLORS["ipwp"],
        ylabel_left="Nino3 amplitude (K)",
        ylabel_right=r"IPWP mean position ($^\circ$E)",
        label_left="Nino3 amplitude",
        label_right="IPWP mean position",
        mean_left=float(np.nanmean(enso)),
        mean_right=float(np.nanmean(bary)),
    )
    save_figure(fig, "03_nino3_ipwp_dual_axis", RESULT_FIG)
    print("IPWP done.")


if __name__ == "__main__":
    main()

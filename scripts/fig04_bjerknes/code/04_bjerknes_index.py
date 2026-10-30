#!/usr/bin/env python3
"""04 Simplified Bjerknes loop gain BJ = μ · β · α (+ 95% CI); compute + plot."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.paths import CI_LEVEL, EXP_IDS, NINO3_MAT, RESULT_FIG, RESULT_MAT, RESULT_TAB, ensure_result_dirs
from lib.stats import product_ci_delta
from lib.style import COLORS, load_nino3, plot_dual_nino3


def main():
    ensure_result_dirs()
    wind = sio.loadmat(RESULT_MAT / "wind_feedback_results.mat")
    beta = sio.loadmat(RESULT_MAT / "ocean_response_beta_results.mat")
    alpha = sio.loadmat(RESULT_MAT / "thermocline_feedback_results.mat")

    # μ: lag=0 (paper dual-axis convention)
    mu = wind["mu"][:, 0].astype(float)
    mu_ci = wind["mu_ci"][:, :, 0].astype(float)

    b = beta["beta_peak"].ravel().astype(float)
    b_ci = beta["beta_ci"].astype(float)

    a = alpha["alpha_peak"].ravel().astype(float)
    a_ci = alpha["alpha_ci"].astype(float)

    n = len(EXP_IDS)
    bj = np.full(n, np.nan)
    bj_ci = np.full((n, 2), np.nan)
    for k in range(n):
        val, ci = product_ci_delta(mu[k], mu_ci[k], b[k], b_ci[k], a[k], a_ci[k])
        bj[k] = val
        bj_ci[k] = ci
        print(
            f"EXP {EXP_IDS[k]:03d} | μ={mu[k]:+.3e}  β={b[k]:+.3e}  α={a[k]:+.3e}  "
            f"BJ={bj[k]:+.4f}  CI=[{bj_ci[k,0]:+.4f},{bj_ci[k,1]:+.4f}]"
        )

    mat_path = RESULT_MAT / "bjerknes_loop_gain_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": np.array(EXP_IDS, float).reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "method": "BJ_loop_gain_mu_x_beta_x_alpha",
            "definition": "BJ = mu(lag0) * beta_peak * alpha_peak (dimensionless)",
            "note": (
                "Simplified three-element Bjerknes loop gain "
                "(Bjerknes 1969; Deppenmeier et al. 2016; "
                "not the full Jin et al. 2006 BJ stability index)"
            ),
            "mu": mu.reshape(-1, 1),
            "mu_ci": mu_ci,
            "beta": b.reshape(-1, 1),
            "beta_ci": b_ci,
            "alpha": a.reshape(-1, 1),
            "alpha_ci": a_ci,
            "BJ": bj.reshape(-1, 1),
            "BJ_ci": bj_ci,
            "BJ_err_method": "relative delta-method from 95% CI half-widths",
        },
    )
    print(f"Saved {mat_path}")

    err_lo = bj - bj_ci[:, 0]
    err_hi = bj_ci[:, 1] - bj

    enso, n_lo, n_hi = load_nino3(NINO3_MAT)
    plot_dual_nino3(
        EXP_IDS, enso, n_lo, n_hi, bj, err_lo, err_hi,
        color_right=COLORS["bj"],
        ylabel_right=r"Bjerknes loop gain $BJ$ (–)",
        label_right=r"$BJ=\mu\beta\alpha$",
        stem="bj_nino3_bjerknes_dual_axis",
        out_dir=RESULT_FIG,
        mean_fmt_right=".3g",
    )

    # Parameter summary table (CSV + MD)
    tab_csv = RESULT_TAB / "feedback_parameters_summary.csv"
    with tab_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "exp_id", "nino3_amp_K", "mu_Pa_per_K", "mu_ci_lo", "mu_ci_hi",
            "beta_m_per_Pa", "beta_ci_lo", "beta_ci_hi", "beta_lag",
            "alpha_K_per_m", "alpha_ci_lo", "alpha_ci_hi", "alpha_lag",
            "BJ", "BJ_ci_lo", "BJ_ci_hi",
        ])
        beta_lag = beta["lag_peak"].ravel()
        alpha_lag = alpha["lag_peak"].ravel()
        for k in range(n):
            w.writerow([
                EXP_IDS[k], f"{enso[k]:.6f}",
                f"{mu[k]:.8e}", f"{mu_ci[k,0]:.8e}", f"{mu_ci[k,1]:.8e}",
                f"{b[k]:.8e}", f"{b_ci[k,0]:.8e}", f"{b_ci[k,1]:.8e}", int(beta_lag[k]),
                f"{a[k]:.8e}", f"{a_ci[k,0]:.8e}", f"{a_ci[k,1]:.8e}", int(alpha_lag[k]),
                f"{bj[k]:.8e}", f"{bj_ci[k,0]:.8e}", f"{bj_ci[k,1]:.8e}",
            ])
    print(f"Saved {tab_csv}")

    # Method selection table
    method_md = RESULT_TAB / "parameter_methods_table.md"
    method_md.write_text(
        """# Parameter selection and calculation table

| Parameter | Symbol | Physical meaning | Region / box | Formula | Lag choice | Units | 95% CI |
|-----------|--------|------------------|--------------|---------|------------|-------|--------|
| Nino3 amplitude | A | ENSO std. amplitude | 5°S–5°N, 210–270°E | χ² CI on std of Niño3 SSTA | — | K | χ² on variance |
| Wind feedback | μ | Atmospheric Bjerknes link: SST → τₓ | SST: Niño3; τₓ: 5°S–5°N, 140–180°E | τₓ' = μ T' (cov/var OLS) | **lag 0** (contemporaneous; lag1 also stored) | Pa K⁻¹ | Student-t on OLS slope |
| Ocean response | β | Wind → thermocline (Kelvin-wave link) | τₓ: 140–180°E; Z₂₀: **Niño3 210–270°E** (5°S–5°N) | Z₂₀' = β τₓ' | peak β among lag 0…6 (wind leads) | m Pa⁻¹ | Student-t on OLS slope |
| Thermocline feedback | α | Thermocline → SST | Z₂₀: 140–260°E; SST: Niño3 | T' = α Z₂₀' | peak α among lag > 0 (Z₂₀ leads) | K m⁻¹ | Student-t on OLS slope |
| Bjerknes loop gain | BJ | Three-element positive-loop gain | product of above | **BJ = μ · β · α** | uses μ(lag0), β_peak, α_peak | dimensionless | relative delta-method from component CIs |

**Literature**
- Three-element loop: Bjerknes (1969); Deppenmeier et al. (2016, Clim. Dyn.); Lübbecke & McPhaden (2013, 2017).
- Atmospheric μ diagnostics: Lloyd et al. (2009); Bellenger et al. (2014).
- Full Jin BJ stability index (Jin et al. 2006) needs mean currents / upwelling / thermal damping — **not** computed here; BJ here is the simplified loop-gain diagnostic only.

**Anomaly definition:** monthly anomaly after removing 12-month climatology (same as original MATLAB).
""",
        encoding="utf-8",
    )
    print(f"Saved {method_md}")
    print(f"mean BJ = {np.nanmean(bj):.4g}")
    print("Bjerknes loop gain done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""09 Two-column combined figure: mean-state drivers | feedbacks.

Layout (2 columns × 6 rows):
  Left  (a–f): Nino3, Walker, τx, IPWP, MLD, stratification
  Right (g–l): WWB, WWV, μ, β, α, BJ   (no duplicate Nino3)

Y-axes on the outer side of each column for a clean journal look.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.paths import EXP_IDS, NINO3_DIR, RESULT_FIG, RESULT_MAT, ROOT, ensure_result_dirs
from lib.style import COLORS

# Sibling project mats
SIBLING = ROOT.parent / "fig04_bjerknes"
WWB_MAT = RESULT_MAT / "wwb_precession_results.mat"
WWV_MAT = RESULT_MAT / "wwv_precession_results.mat"
WIND_MAT = SIBLING / "result" / "mat" / "wind_feedback_results.mat"
BETA_MAT = SIBLING / "result" / "mat" / "ocean_response_beta_results.mat"
ALPHA_MAT = SIBLING / "result" / "mat" / "thermocline_feedback_results.mat"
BJ_MAT = SIBLING / "result" / "mat" / "bjerknes_loop_gain_results.mat"

OUT_STEM = "09_combined_meanstate_feedbacks_2col"
# ~ journal 2-col width, tall enough for 6 rows
FIGSIZE = (26.0 / 2.54, 20.0 / 2.54)

# Extra colors used only on the right column
COLORS_FB = {
    "wwb": np.array([0.55, 0.20, 0.35]),
    "wwv": np.array([0.20, 0.35, 0.55]),
    "mu": np.array([0.12, 0.35, 0.58]),
    "beta": np.array([0.10, 0.45, 0.55]),
    "alpha": np.array([0.10, 0.45, 0.38]),
    "bj": np.array([0.55, 0.25, 0.55]),
}


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
            "font.size": 9,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
        }
    )


def _phase_ticklabels(every: int = 2):
    labs = []
    for i, v in enumerate(EXP_IDS):
        labs.append(f"{v:03d}\u00b0" if i % every == 0 else "")
    return labs


def _panel(
    ax,
    s,
    *,
    row: int,
    n_rows: int,
    y_on_left: bool,
    show_xtop: bool,
    show_xbot: bool,
):
    c = tuple(s["color"])
    y = np.asarray(s["y"], float)
    elo = np.asarray(s["elo"], float)
    ehi = np.asarray(s["ehi"], float)
    has_err = s.get("has_err", True)

    if has_err and (np.any(elo > 0) or np.any(ehi > 0)):
        ax.errorbar(
            EXP_IDS, y, yerr=np.vstack([elo, ehi]),
            fmt="-o", color=c, lw=1.35, ms=3.6, mfc=c, mec="w", mew=0.35,
            capsize=2.2, elinewidth=0.9, zorder=3,
        )
        y0 = np.nanmin(y - elo)
        y1 = np.nanmax(y + ehi)
    else:
        ax.plot(
            EXP_IDS, y, "-o", color=c, lw=1.35, ms=3.6, mfc=c, mec="w", mew=0.35, zorder=3,
        )
        y0, y1 = np.nanmin(y), np.nanmax(y)

    ax.axhline(np.nanmean(y), ls="--", color=c, lw=0.95, alpha=0.8)
    pad = max((y1 - y0) * 0.18, 1e-6)
    ymin, ymax = y0 - pad, y1 + pad
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-5, 350)
    ax.set_xticks(EXP_IDS)
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    ax.grid(axis="y", visible=False)
    ax.set_facecolor("none")

    ax.spines["top"].set_visible(row == 0)
    ax.spines["bottom"].set_visible(row == n_rows - 1)
    ax.spines["top"].set_linewidth(1.05)
    ax.spines["bottom"].set_linewidth(1.05)

    if y_on_left:
        ax.yaxis.set_label_position("left")
        ax.yaxis.tick_left()
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color(c)
        ax.spines["left"].set_linewidth(1.35)
        ax.spines["left"].set_bounds(ymin, ymax)
    else:
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(c)
        ax.spines["right"].set_linewidth(1.35)
        ax.spines["right"].set_bounds(ymin, ymax)

    ax.tick_params(axis="y", colors=c, labelsize=7.5, direction="in", length=3)
    ax.set_ylabel(s["ylabel"], color=c, fontsize=8.5, fontweight="bold", labelpad=2)

    if show_xtop:
        ax.tick_params(
            axis="x", bottom=False, top=True, labelbottom=False, labeltop=True,
            labelsize=6.5, direction="in", pad=1.5,
        )
        ax.set_xticklabels(_phase_ticklabels(2), fontsize=6.5)
    elif show_xbot:
        ax.tick_params(
            axis="x", bottom=True, top=False, labelbottom=True, labeltop=False,
            labelsize=6.5, direction="in", pad=1.5,
        )
        ax.set_xticklabels(_phase_ticklabels(2), fontsize=6.5)
        ax.set_xlabel(r"Precession phase ($^\circ$)", fontsize=8.5, fontweight="bold", labelpad=3)
    else:
        ax.tick_params(axis="x", bottom=False, top=False, labelbottom=False)
        ax.set_xticklabels([])

    ax.text(
        0.015, 0.90, s["tag"], transform=ax.transAxes,
        fontsize=9.5, fontweight="bold", va="top", ha="left",
    )


def load_left_column():
    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    walk = sio.loadmat(RESULT_MAT / "walker_edges_precession_results.mat")
    ustrw = sio.loadmat(RESULT_MAT / "ustrw_precession_results.mat")
    ipwp = sio.loadmat(RESULT_MAT / "ipwp_precession_results.mat")
    mld = sio.loadmat(RESULT_MAT / "mld_precession_results.mat")
    strat = sio.loadmat(RESULT_MAT / "stratification_precession_results.mat")

    enso = nino["enso_amp"].ravel().astype(float)
    nci = nino["chi2_ci_std"]
    center = walk["walk_center_lon"].ravel().astype(float)
    u_mean = ustrw["ustrw_mean"].ravel().astype(float)
    ipwp_y = ipwp["barycenter"].ravel().astype(float)
    ipwp_e = ipwp["error24"].ravel().astype(float)
    m_mean = mld["mld_mean"].ravel().astype(float)
    s_mean = strat["stratification"].ravel().astype(float)

    z = np.zeros(len(EXP_IDS))
    series = [
        dict(y=enso, elo=enso - nci[:, 0], ehi=nci[:, 1] - enso,
             color=COLORS["nino3"], ylabel="Nino3 amp. (K)", tag="(a)", has_err=True),
        dict(y=center, elo=z, ehi=z,
             color=COLORS["walker"], ylabel=r"Walker center ($^\circ$E)", tag="(b)", has_err=False),
        dict(y=u_mean, elo=u_mean - ustrw["ustrw_ci"][:, 0], ehi=ustrw["ustrw_ci"][:, 1] - u_mean,
             color=COLORS["ustrw"], ylabel=r"$\tau_x$ (Pa)", tag="(c)", has_err=True),
        dict(y=ipwp_y, elo=ipwp_e, ehi=ipwp_e,
             color=COLORS["ipwp"], ylabel=r"IPWP pos. ($^\circ$E)", tag="(d)", has_err=True),
        dict(y=m_mean, elo=m_mean - mld["mld_ci_mean"][:, 0], ehi=mld["mld_ci_mean"][:, 1] - m_mean,
             color=COLORS["mld"], ylabel="MLD (m)", tag="(e)", has_err=True),
        dict(y=s_mean, elo=s_mean - strat["strat_ci"][:, 0], ehi=strat["strat_ci"][:, 1] - s_mean,
             color=COLORS["strat"], ylabel=r"Stratification ($^\circ$C)", tag="(f)", has_err=True),
    ]
    return series


def load_right_column():
    for p in (WWB_MAT, WWV_MAT, WIND_MAT, BETA_MAT, ALPHA_MAT, BJ_MAT):
        if not Path(p).is_file():
            raise FileNotFoundError(p)

    wwb = sio.loadmat(WWB_MAT)
    wwv = sio.loadmat(WWV_MAT)
    wind = sio.loadmat(WIND_MAT)
    beta = sio.loadmat(BETA_MAT)
    alpha = sio.loadmat(ALPHA_MAT)
    bj = sio.loadmat(BJ_MAT)

    inten = wwb["wwb_inten"].ravel().astype(float)
    inten_ci = wwb["wwb_inten_ci"]
    wwv_y = wwv["wwv_mean"].ravel().astype(float)
    wwv_ci = wwv["wwv_ci"]
    mu = wind["mu"][:, 0].astype(float)
    mu_ci = wind["mu_ci"][:, :, 0]
    b = beta["beta_peak"].ravel().astype(float)
    b_ci = beta["beta_ci"]
    a = alpha["alpha_peak"].ravel().astype(float)
    a_ci = alpha["alpha_ci"]
    bj_y = bj["BJ"].ravel().astype(float)
    bj_ci = bj["BJ_ci"]

    return [
        dict(y=inten, elo=inten - inten_ci[:, 0], ehi=inten_ci[:, 1] - inten,
             color=COLORS_FB["wwb"], ylabel="WWB intensity (Pa)", tag="(g)", has_err=True),
        dict(y=wwv_y, elo=wwv_y - wwv_ci[:, 0], ehi=wwv_ci[:, 1] - wwv_y,
             color=COLORS_FB["wwv"], ylabel="WWV (m)", tag="(h)", has_err=True),
        dict(y=mu, elo=mu - mu_ci[:, 0], ehi=mu_ci[:, 1] - mu,
             color=COLORS_FB["mu"], ylabel=r"$\mu$ (Pa K$^{-1}$)", tag="(i)", has_err=True),
        dict(y=b, elo=b - b_ci[:, 0], ehi=b_ci[:, 1] - b,
             color=COLORS_FB["beta"], ylabel=r"$\beta$ (m Pa$^{-1}$)", tag="(j)", has_err=True),
        dict(y=a, elo=a - a_ci[:, 0], ehi=a_ci[:, 1] - a,
             color=COLORS_FB["alpha"], ylabel=r"$\alpha$ (K m$^{-1}$)", tag="(k)", has_err=True),
        dict(y=bj_y, elo=bj_y - bj_ci[:, 0], ehi=bj_ci[:, 1] - bj_y,
             color=COLORS_FB["bj"], ylabel=r"$BJ=\mu\beta\alpha$", tag="(l)", has_err=True),
    ]


def plot_two_column(left_series, right_series):
    apply_style()
    n_rows = 6
    assert len(left_series) == n_rows and len(right_series) == n_rows

    # margins: room for outer y-labels; slightly tighter center gap
    top, bottom = 0.905, 0.07
    gap = 0.010
    usable = top - bottom - gap * (n_rows - 1)
    h = usable / n_rows
    # left / right columns a bit closer
    x0, x1, x2, x3 = 0.09, 0.495, 0.525, 0.93

    fig = plt.figure(figsize=FIGSIZE)

    # Column titles (slightly lower than before)
    fig.text(
        (x0 + x1) / 2, 0.948, "Mean state",
        ha="center", va="center", fontsize=11, fontweight="bold", color="0.25",
    )
    fig.text(
        (x2 + x3) / 2, 0.948, "Coupled feedbacks",
        ha="center", va="center", fontsize=11, fontweight="bold", color="0.25",
    )

    for i in range(n_rows):
        y0 = top - (i + 1) * h - i * gap
        ax_l = fig.add_axes([x0, y0, x1 - x0, h])
        ax_r = fig.add_axes([x2, y0, x3 - x2, h])

        _panel(
            ax_l, left_series[i],
            row=i, n_rows=n_rows, y_on_left=True,
            show_xtop=(i == 0), show_xbot=(i == n_rows - 1),
        )
        _panel(
            ax_r, right_series[i],
            row=i, n_rows=n_rows, y_on_left=False,
            show_xtop=(i == 0), show_xbot=(i == n_rows - 1),
        )

    out = RESULT_FIG / f"{OUT_STEM}.png"
    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    print("saved", out)
    return out


def main():
    ensure_result_dirs()
    left = load_left_column()
    right = load_right_column()
    plot_two_column(left, right)
    print("Two-column combined figure done.")


if __name__ == "__main__":
    main()

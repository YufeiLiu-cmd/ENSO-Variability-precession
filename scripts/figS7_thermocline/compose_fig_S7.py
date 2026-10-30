#!/usr/bin/env python3
"""
Supplementary Figure S7 — panels (a)(c)(e) only.
Style strictly follows original A6thermocline/thermocline.ipynb (thermocline_supply.png):
  figsize=(10,10), Arial, fontweight normal, axes.linewidth=1.5,
  depth 0–220 m, hspace=0.18, colorbar positions as in A6.
Anomaly = 90° − 270° (thetao105 − thetao285).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
DATA = Path(
    __import__("os").environ.get(
        "ENSO_DATA_ROOT",
        str(HERE.parents[2] / "data"),
    )
) / "figS7"
RGB_B2R = DATA / "cmp_b2r.rgb"

LON0, LON1 = 140.0, 210.0
Z_UPPER, Z_REF = 100.0, 120.0


def load_anomaly_cmap():
    """A6 panel (c) used cmaps.div1_blue_orange; fall back to local cmp_b2r."""
    raw = np.loadtxt(RGB_B2R)
    if raw.max() > 1.5:
        raw = raw / 255.0
    x_old = np.linspace(0, 1, len(raw))
    x_new = np.linspace(0, 1, 256)
    cols = np.column_stack([np.interp(x_new, x_old, raw[:, i]) for i in range(3)])
    return LinearSegmentedColormap.from_list("cmp_b2r", cols)


def cai_stratification(T, lon, depth) -> float:
    mask = (lon >= LON0) & (lon <= LON1)
    di = np.arange(0.0, 341.0, 1.0)
    vals = []
    for i in np.where(mask)[0]:
        p = T[i].astype(float)
        valid = ~np.isnan(p)
        if valid.sum() < 2:
            continue
        pi = np.interp(di, depth[valid], p[valid], left=np.nan, right=np.nan)
        vals.append(np.nanmean(pi[di <= Z_UPPER]) - pi[int(Z_REF)])
    return float(np.nanmean(vals))


def compose() -> None:
    # —— exact A6 rcParams sequence ——
    fig = plt.figure(figsize=(10, 10))
    plt.rcParams["figure.figsize"] = (16, 4)
    plt.rcParams["font.size"] = 11
    plt.rcParams["axes.linewidth"] = 1.5
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.weight"] = "normal"
    plt.rcParams.update({"font.size": 16})

    cmap = load_anomaly_cmap()

    tm = sio.loadmat(DATA / "thermocline_diff.mat")
    lon = tm["lon_need"].ravel().astype(float)
    depth = tm["depth_need"].ravel().astype(float)
    X, Y = np.meshgrid(lon, -depth)

    # 90−270
    specs = [
        ("(a) JJA Anomaly (90° − 270°)", (-tm["diff_JJA"]).T, tm["h_JJA"],
         tm["JJA"][:, :, 0].T, tm["JJA"][:, :, 1].T,
         tm["JJA"][:, :, 0], tm["JJA"][:, :, 1], "JJA"),
        ("(b) DJF Anomaly (90° − 270°)", (-tm["diff_NDJ"]).T, tm["h_NDJ"],
         tm["NDJ"][:, :, 0].T, tm["NDJ"][:, :, 1].T,
         tm["NDJ"][:, :, 0], tm["NDJ"][:, :, 1], "DJF"),
        ("(c) Climatological Anomaly (90° − 270°)", (-tm["diff_annual"]).T, tm["h_annual"],
         tm["annual"][:, :, 0].T, tm["annual"][:, :, 1].T,
         tm["annual"][:, :, 0], tm["annual"][:, :, 1], "Annual"),
    ]

    rows = []
    for lab, diff, h, T90, T270, T90_ld, T270_ld, tag in specs:
        s90 = cai_stratification(T90_ld, lon, depth)
        s270 = cai_stratification(T270_ld, lon, depth)
        rows.append((tag, s90, s270, s90 - s270))
        print(f"{tag}: Strat90={s90:.3f} Strat270={s270:.3f} Δ={s90-s270:.3f}")

    np.savetxt(
        HERE / "cai_stratification_summary.csv",
        rows,
        fmt="%s",
        delimiter=",",
        header="season,strat_90,strat_270,delta_90_minus_270",
        comments="",
    )
    dStrat = {r[0]: r[3] for r in rows}

    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1])
    plt.subplots_adjust(hspace=0.18)

    xtick = np.arange(40, 300, 20)
    xtlab = [
        "40°E", "60°E", "80°E", "100°E", "120°E", "140°E", "160°E",
        "180°", "160°W", "140°W", "120°W", "100°W", "80°W",
    ]
    ytick = np.arange(-220, 1, 20)
    ytl = ["220", "", "180", "", "140", "", "100", "", "60", "", "20", "0"]

    # A6 anomaly levels
    levels = np.arange(-4, 4, 0.1)
    cfs = []

    for i, (lab, diff, h, T90, T270, T90_ld, T270_ld, tag) in enumerate(specs):
        ax = fig.add_subplot(gs[i, 0])

        # filled anomaly — same as A6 panel (c)
        cf = ax.contourf(X, Y, diff, levels=levels, cmap=cmap, extend="both")
        cfs.append(cf)

        # 20°C isotherms — A6 panel (c) style
        h90 = ax.contour(
            X, Y, T90, levels=[20], colors="b", linestyles="dashed", linewidths=2
        )
        ax.clabel(h90, inline=True, fontsize=11, fmt="20 C°", colors="b")
        h270 = ax.contour(
            X, Y, T270, levels=[20], colors="r", linewidths=2
        )
        ax.clabel(h270, inline=True, fontsize=11, fmt="20 C°", colors="r")

        # significance stippling — A6 style
        levh, lath = np.where(np.asarray(h).T == 1)
        ax.scatter(
            lon[lath],
            -depth[levh],
            color="lightgray",
            s=0.1,
            marker="p",
            alpha=0.8,
        )

        # Cai domain (added content; line weights match A6 contour lw=2)
        ax.add_patch(
            Rectangle(
                (LON0, -Z_UPPER),
                LON1 - LON0,
                Z_UPPER - 1.0,
                fill=False,
                edgecolor="g",
                linewidth=2,
                linestyle="-",
            )
        )
        ax.plot([LON0, LON1], [-Z_REF, -Z_REF], "m-", lw=2, alpha=0.8)
        ax.text(
            165,
            -88,
            f"ΔStrat = {dStrat[tag]:.3f} °C",
            fontsize=11,
            fontweight="normal",
            color="k",
            bbox=dict(facecolor="white", alpha=0.6, edgecolor="none"),
        )

        ax.set_xlim([40, 280])
        ax.set_ylim([-220, 0])
        ax.set_xticks(xtick)
        ax.set_yticks(ytick)
        ax.set_yticklabels(ytl, fontsize=11, fontweight="normal")
        ax.set_ylabel("Depth [m]", fontsize=11, color="k", fontweight="normal")
        ax.text(40, 5, lab, color="black", fontsize=11, weight="normal")

        if i < 2:
            ax.set_xticklabels(xtlab, rotation=0, fontsize=3, fontweight="normal")
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xticklabels(xtlab, rotation=0, fontsize=11, fontweight="normal")
            ax.set_xlabel("Longitude", fontsize=11, color="k", fontweight="normal")

        for sp in ax.spines.values():
            sp.set_linewidth(1.5)

    # Colorbar — A6 anomaly colorbar geometry (bottom right)
    # With three anomaly panels of identical scale, use one bar spanning lower+mid
    # plus match A6 dual-bar positions: tall upper + short lower, both anomaly.
    cax1 = fig.add_axes([0.92, 0.41, 0.015, 0.45])
    cb1 = fig.colorbar(cfs[0], cax=cax1, orientation="vertical")
    cb1.set_label("Temperature Anomaly [C°]", fontsize=11, fontweight="normal")
    cb1.set_ticks(np.arange(-4, 4.01, 2))
    cb1.set_ticklabels([int(t) for t in np.arange(-4, 4.01, 2)], fontsize=11)
    for t in cb1.ax.get_yticklabels():
        t.set_fontweight("normal")
        t.set_fontsize(11)

    cax2 = fig.add_axes([0.92, 0.12, 0.015, 0.2])
    cb2 = fig.colorbar(cfs[2], cax=cax2, orientation="vertical")
    cb2.set_label("Temperature Anomaly [C°]", fontsize=11, fontweight="normal")
    cb2.set_ticks(np.arange(-4, 4.01, 2))
    cb2.set_ticklabels([int(t) for t in np.arange(-4, 4.01, 2)], fontsize=11)
    for t in cb2.ax.get_yticklabels():
        t.set_fontweight("normal")
        t.set_fontsize(11)

    # Default dpi=100 → 1000×1000 like thermocline_supply.png
    out_png = HERE / "Fig_S7_thermocline_stratification.png"
    out_pdf = HERE / "Fig_S7_thermocline_stratification.pdf"
    fig.savefig(out_png, dpi=100, facecolor="white")
    fig.savefig(out_pdf, dpi=100, facecolor="white")
    # hi-res copy for manuscript insertion
    fig.savefig(HERE / "Fig_S7_thermocline_stratification_300dpi.png", dpi=300, facecolor="white")
    plt.rcParams.update({"font.size": 11})
    plt.close(fig)
    print(f"Saved: {out_png} (1000×1000, A6 match)")
    print(f"Saved: {out_pdf}")
    print(f"Saved: {HERE / 'Fig_S7_thermocline_stratification_300dpi.png'}")


if __name__ == "__main__":
    compose()

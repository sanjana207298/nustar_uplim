"""
visualization.py
----------------
Step-by-step diagnostic plots for NuSTAR upper limit analysis.

All functions produce self-contained figures and either display them
(show=True) or save to a specified output directory.
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm, PowerNorm
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization import ImageNormalize, ZScaleInterval, LogStretch


# ── Colour scheme ──────────────────────────────────────────────────────────
BG_DARK    = "#0f0f1a"
BG_PANEL   = "#1a1a2e"
SRC_COLOR  = "#ff6ec7"   # magenta-pink
BKG_COLOR  = "#00e5ff"   # cyan
ACCENT     = "#ffd700"   # gold
TEXT_COLOR = "#e0e0e0"


def _setup_style():
    plt.rcParams.update({
        "figure.facecolor":  BG_PANEL,
        "axes.facecolor":    BG_DARK,
        "text.color":        TEXT_COLOR,
        "axes.labelcolor":   TEXT_COLOR,
        "xtick.color":       TEXT_COLOR,
        "ytick.color":       TEXT_COLOR,
        "axes.edgecolor":    "gray",
        "grid.color":        "gray",
        "grid.alpha":        0.3,
        "font.family":       "monospace",
        "axes.titlepad":     10,
    })


def _save_or_show(fig, out_dir, filename, show):
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, filename)
        fig.savefig(path, dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"  Saved: {path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────
# Step 1 — Raw sky image + source position
# ──────────────────────────────────────────────────────────────────────────

def plot_sky_image(image, wcs, src_ra, src_dec,
                   src_r_pix, src_cx, src_cy,
                   obsid="", module="", band="full",
                   out_dir=None, show=True):
    """
    Plot the NuSTAR sky image with the target source circle overlaid.
    """
    _setup_style()
    fig = plt.figure(figsize=(9, 8), facecolor=BG_PANEL)

    try:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        _use_wcs = True
    except Exception:
        ax = fig.add_subplot(1, 1, 1)
        _use_wcs = False

    ax.set_facecolor(BG_DARK)

    # Log-scale normalisation
    img_log = np.log1p(image.copy())
    vmin, vmax = np.nanpercentile(img_log[img_log > 0], [5, 99.5])

    ax.imshow(img_log, origin="lower", cmap="afmhot",
              vmin=vmin, vmax=vmax, interpolation="nearest")

    # Source circle
    circ = mpatches.Circle(
        (src_cx, src_cy), src_r_pix,
        transform=ax.get_transform("pixel") if _use_wcs else ax.transData,
        edgecolor=SRC_COLOR, facecolor=SRC_COLOR, linewidth=2,
        linestyle="-", alpha=0.2, zorder=4,
    )
    outline = mpatches.Circle(
        (src_cx, src_cy), src_r_pix,
        transform=ax.get_transform("pixel") if _use_wcs else ax.transData,
        edgecolor=SRC_COLOR, facecolor="none", linewidth=2,
        linestyle="-", zorder=5,
    )
    ax.add_patch(circ)
    ax.add_patch(outline)
    ax.plot(src_cx, src_cy,
            transform=ax.get_transform("pixel") if _use_wcs else ax.transData,
            marker="+", color=SRC_COLOR, markersize=14, markeredgewidth=2, zorder=6)

    ax.set_title(
        f"Sky Image  ·  ObsID {obsid}  FPM{module}  ·  Band: {band.upper()}\n"
        f"Target RA={src_ra:.5f}°  Dec={src_dec:.5f}°  |  r_src={src_r_pix:.1f} px",
        color=TEXT_COLOR, fontsize=10,
    )

    if _use_wcs:
        ax.coords[0].set_axislabel("RA")
        ax.coords[1].set_axislabel("Dec")
    else:
        ax.set_xlabel("X (pixels)")
        ax.set_ylabel("Y (pixels)")

    # Legend
    src_p = mpatches.Patch(edgecolor=SRC_COLOR, facecolor=SRC_COLOR,
                            alpha=0.4, linewidth=2, label="Source region")
    ax.legend(handles=[src_p], fontsize=9, facecolor="#222",
              edgecolor="gray", labelcolor="white", loc="upper right")

    fig.tight_layout()
    _save_or_show(fig, out_dir, f"step1_sky_image_{module}.png", show)


# ──────────────────────────────────────────────────────────────────────────
# Step 2 — Source + background regions
# ──────────────────────────────────────────────────────────────────────────

def plot_regions(image, wcs,
                 src_cx, src_cy, src_r_pix,
                 bkg_cx, bkg_cy, bkg_r_pix,
                 bkg_type="circle", bkg_inner_pix=None,
                 obsid="", module="",
                 out_dir=None, show=True):
    """
    Plot sky image with both source and background regions overlaid.
    """
    _setup_style()
    fig, ax = plt.subplots(figsize=(9, 8), facecolor=BG_PANEL)
    ax.set_facecolor(BG_DARK)
    fig.patch.set_facecolor(BG_PANEL)

    img_log = np.log1p(image.copy())
    vmin, vmax = np.nanpercentile(img_log[img_log > 0], [5, 99.5])
    ax.imshow(img_log, origin="lower", cmap="afmhot",
              vmin=vmin, vmax=vmax, interpolation="nearest")

    # Source circle
    for kw in [dict(facecolor=SRC_COLOR, alpha=0.18), dict(facecolor="none")]:
        ax.add_patch(mpatches.Circle(
            (src_cx, src_cy), src_r_pix,
            edgecolor=SRC_COLOR, linewidth=2, linestyle="-", zorder=5, **kw))
    ax.plot(src_cx, src_cy, "+", color=SRC_COLOR,
            markersize=14, markeredgewidth=2, zorder=6)

    # Background region
    if bkg_type == "annulus" and bkg_inner_pix:
        for kw in [dict(facecolor=BKG_COLOR, alpha=0.12), dict(facecolor="none")]:
            ax.add_patch(mpatches.Annulus(
                (bkg_cx, bkg_cy), bkg_r_pix, bkg_r_pix - bkg_inner_pix,
                edgecolor=BKG_COLOR, linewidth=2, linestyle="--", zorder=4, **kw))
    else:
        for kw in [dict(facecolor=BKG_COLOR, alpha=0.12), dict(facecolor="none")]:
            ax.add_patch(mpatches.Circle(
                (bkg_cx, bkg_cy), bkg_r_pix,
                edgecolor=BKG_COLOR, linewidth=2, linestyle="--", zorder=4, **kw))
        ax.plot(bkg_cx, bkg_cy, "x", color=BKG_COLOR,
                markersize=10, markeredgewidth=2, zorder=5)

    ax.set_title(
        f"Source & Background Regions  ·  ObsID {obsid}  FPM{module}",
        color=TEXT_COLOR, fontsize=11,
    )
    ax.set_xlabel("X (pixels)", color=TEXT_COLOR)
    ax.set_ylabel("Y (pixels)", color=TEXT_COLOR)

    src_p = mpatches.Patch(edgecolor=SRC_COLOR, facecolor=SRC_COLOR,
                            alpha=0.4, linewidth=2, label="Source")
    bkg_p = mpatches.Patch(edgecolor=BKG_COLOR, facecolor=BKG_COLOR,
                            alpha=0.4, linewidth=2, label="Background")
    ax.legend(handles=[src_p, bkg_p], fontsize=9, facecolor="#222",
              edgecolor="gray", labelcolor="white", loc="upper right")

    fig.tight_layout()
    _save_or_show(fig, out_dir, f"step2_regions_{module}.png", show)


# ──────────────────────────────────────────────────────────────────────────
# Step 3 — Event scatter in source + background
# ──────────────────────────────────────────────────────────────────────────

def plot_event_scatter(ev_x, ev_y, src_mask, bkg_mask,
                       src_cx, src_cy, src_r_pix,
                       bkg_cx, bkg_cy, bkg_r_pix, bkg_type,
                       bkg_inner_pix=None,
                       obsid="", module="", band="full",
                       out_dir=None, show=True):
    """
    Scatter plot of individual events colour-coded by region membership.
    """
    _setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor(BG_PANEL)
    fig.suptitle(
        f"Event Distribution  ·  ObsID {obsid}  FPM{module}  ·  Band: {band.upper()}",
        color=TEXT_COLOR, fontsize=12,
    )

    for ax, (mask, cx, cy, r, label, color) in zip(
        axes,
        [
            (src_mask, src_cx, src_cy, src_r_pix, "Source Region", SRC_COLOR),
            (bkg_mask, bkg_cx, bkg_cy, bkg_r_pix, "Background Region", BKG_COLOR),
        ],
    ):
        ax.set_facecolor(BG_DARK)
        x_sel = ev_x[mask]
        y_sel = ev_y[mask]

        # All events as faint background
        ax.scatter(ev_x, ev_y, s=0.3, c="white", alpha=0.05, rasterized=True)
        # Region events
        ax.scatter(x_sel, y_sel, s=1.5, c=color, alpha=0.7, rasterized=True,
                   label=f"{label} ({int(mask.sum())} cts)")

        # Region outline
        ax.add_patch(mpatches.Circle(
            (cx, cy), r, edgecolor=color, facecolor="none",
            linewidth=2, linestyle="--", zorder=5))
        ax.plot(cx, cy, "+", color=color, markersize=12, markeredgewidth=2, zorder=6)

        ax.set_xlim(cx - r * 3, cx + r * 3)
        ax.set_ylim(cy - r * 3, cy + r * 3)
        ax.set_title(label, color=TEXT_COLOR, fontsize=10)
        ax.set_xlabel("X (pixels)", color=TEXT_COLOR)
        ax.set_ylabel("Y (pixels)", color=TEXT_COLOR)
        ax.legend(fontsize=8, facecolor="#222", edgecolor="gray", labelcolor="white")

    fig.tight_layout()
    _save_or_show(fig, out_dir, f"step3_event_scatter_{module}.png", show)


# ──────────────────────────────────────────────────────────────────────────
# Step 4 — Upper limit summary plot
# ──────────────────────────────────────────────────────────────────────────

def plot_upper_limit_summary(stats_results_A, stats_results_B=None,
                              obsid="", band="full",
                              out_dir=None, show=True):
    """
    Summary bar chart of upper limits for FPMA (and optionally FPMB).
    """
    _setup_style()

    modules = {"FPM-A": stats_results_A}
    if stats_results_B is not None:
        modules["FPM-B"] = stats_results_B

    n_modules = len(modules)
    fig, axes = plt.subplots(1, n_modules, figsize=(7 * n_modules, 6),
                              squeeze=False)
    fig.patch.set_facecolor(BG_PANEL)
    fig.suptitle(
        f"Upper Limit Summary  ·  ObsID {obsid}  ·  Band: {band.upper()}",
        color=TEXT_COLOR, fontsize=12,
    )

    for col, (mod_name, res) in enumerate(modules.items()):
        ax = axes[0][col]
        ax.set_facecolor(BG_DARK)

        ul = res["upper_limits"]
        labels  = list(ul.keys())
        cr_vals = [ul[k]["count_rate_ul"] for k in labels]
        # flux_ul removed — use WebPIMMS for flux conversion

        colors = ["#ffa07a", "#ff6347", "#dc143c"][:len(labels)]

        bars = ax.bar(labels, cr_vals, color=colors, edgecolor="white",
                      linewidth=0.8, alpha=0.85)
        for bar, v in zip(bars, cr_vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.03,
                    f"{v:.2e}", ha="center", va="bottom",
                    color="white", fontsize=9)

        # Text box with diagnostics
        info = (
            f"Src counts: {res['src_counts']}\n"
            f"Bkg counts: {res['bkg_counts']}\n"
            f"α (area ratio): {res['alpha']:.4f}\n"
            f"Expected bkg: {res['expected_bkg']:.2f}\n"
            f"Net counts: {res['net_counts']:.1f}\n"
            f"Li-Ma signif.: {res['lima_sig']:.2f}σ\n"
            f"Exposure: {res['exposure_s']:.0f} s"
        )
        ax.text(0.97, 0.97, info, transform=ax.transAxes,
                fontsize=8, va="top", ha="right",
                color=TEXT_COLOR,
                bbox=dict(boxstyle="round,pad=0.4", fc="#111", ec="gray", alpha=0.9),
                fontfamily="monospace")

        ax.set_title(mod_name, color=TEXT_COLOR, fontsize=11)
        ax.set_ylabel("Count-rate upper limit [cts/s]", color=TEXT_COLOR)
        ax.set_xlabel("Confidence level", color=TEXT_COLOR)
        ax.grid(axis="y", alpha=0.2, color="gray")

    fig.tight_layout()
    _save_or_show(fig, out_dir, f"step4_upper_limits_{band}.png", show)


# ──────────────────────────────────────────────────────────────────────────
# Step 5 — Radial profile (source vs background)
# ──────────────────────────────────────────────────────────────────────────

def plot_radial_profile(ev_x, ev_y, cx, cy, r_max,
                        obsid="", module="", band="full",
                        out_dir=None, show=True):
    """
    Annular radial profile of counts around the source position.
    """
    _setup_style()
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(BG_PANEL)
    ax.set_facecolor(BG_DARK)

    radii = np.sqrt((ev_x - cx)**2 + (ev_y - cy)**2)
    bins  = np.linspace(0, r_max * 2, 30)
    areas = np.pi * (bins[1:]**2 - bins[:-1]**2)
    counts, _ = np.histogram(radii, bins=bins)
    surf_density = counts / areas

    bin_centres = (bins[:-1] + bins[1:]) / 2
    colors = np.where(bin_centres <= r_max, SRC_COLOR, BKG_COLOR)

    ax.step(bin_centres, surf_density, where="mid", color=ACCENT,
            linewidth=1.5, label="Surface density")
    ax.scatter(bin_centres, surf_density, c=colors, s=20, zorder=5)

    ax.axvline(r_max, color=SRC_COLOR, linestyle="--", linewidth=1.5,
               label=f"Source radius ({r_max:.0f} px)")
    ax.set_xlabel("Radius from source position (pixels)", color=TEXT_COLOR)
    ax.set_ylabel("Count surface density [cts/px²]", color=TEXT_COLOR)
    ax.set_title(
        f"Radial Profile  ·  ObsID {obsid}  FPM{module}  ·  Band: {band.upper()}",
        color=TEXT_COLOR,
    )
    ax.legend(fontsize=9, facecolor="#222", edgecolor="gray", labelcolor="white")
    ax.grid(alpha=0.2)

    fig.tight_layout()
    _save_or_show(fig, out_dir, f"step5_radial_profile_{module}.png", show)


# ──────────────────────────────────────────────────────────────────────────
# Final — Combined report figure
# ──────────────────────────────────────────────────────────────────────────

def plot_combined_report(image_A, wcs_A,
                          src_cx_A, src_cy_A, src_r_A,
                          bkg_cx_A, bkg_cy_A, bkg_r_A,
                          stats_A,
                          image_B=None, wcs_B=None,
                          src_cx_B=None, src_cy_B=None, src_r_B=None,
                          bkg_cx_B=None, bkg_cy_B=None, bkg_r_B=None,
                          stats_B=None,
                          obsid="", src_ra=0.0, src_dec=0.0, band="full",
                          out_dir=None, show=True):
    """
    Single combined report figure: images + regions + UL bar chart.
    """
    _setup_style()
    has_B = image_B is not None and stats_B is not None

    n_cols = 3 if has_B else 2
    fig = plt.figure(figsize=(7 * n_cols, 12), facecolor=BG_PANEL)
    gs  = gridspec.GridSpec(2, n_cols, figure=fig,
                            hspace=0.35, wspace=0.3)

    # ── FPMA image ──
    ax_A = fig.add_subplot(gs[0, 0])
    _plot_image_with_regions(ax_A, image_A,
                             src_cx_A, src_cy_A, src_r_A,
                             bkg_cx_A, bkg_cy_A, bkg_r_A,
                             title=f"FPM-A  ·  {band.upper()}")

    # ── FPMB image ──
    if has_B:
        ax_B = fig.add_subplot(gs[0, 1])
        _plot_image_with_regions(ax_B, image_B,
                                 src_cx_B, src_cy_B, src_r_B,
                                 bkg_cx_B, bkg_cy_B, bkg_r_B,
                                 title=f"FPM-B  ·  {band.upper()}")

    # ── UL bar chart ──
    ax_ul = fig.add_subplot(gs[0, -1])
    ax_ul.set_facecolor(BG_DARK)
    _plot_ul_bars(ax_ul, stats_A, stats_B if has_B else None, band)

    # ── Text summary ──
    ax_txt = fig.add_subplot(gs[1, :])
    ax_txt.set_facecolor(BG_DARK)
    ax_txt.axis("off")
    summary = _build_text_summary(stats_A, stats_B if has_B else None,
                                  obsid, src_ra, src_dec, band)
    ax_txt.text(0.02, 0.95, summary,
                transform=ax_txt.transAxes,
                fontsize=9, va="top", ha="left",
                color=TEXT_COLOR, fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.5", fc="#0a0a15", ec="gray"))

    fig.suptitle(
        f"NuSTAR Non-Detection Upper Limit Report\n"
        f"ObsID: {obsid}   RA={src_ra:.5f}°  Dec={src_dec:.5f}°",
        color=ACCENT, fontsize=13, fontweight="bold", y=0.98,
    )

    _save_or_show(fig, out_dir, "nustar_uplim_report.png", show)
    return fig


def _plot_image_with_regions(ax, image, src_cx, src_cy, src_r,
                              bkg_cx, bkg_cy, bkg_r, title=""):
    ax.set_facecolor(BG_DARK)
    img_log = np.log1p(image)
    vmin, vmax = np.nanpercentile(img_log[img_log > 0], [5, 99.5])
    ax.imshow(img_log, origin="lower", cmap="afmhot",
              vmin=vmin, vmax=vmax, interpolation="nearest")

    for (cx, cy, r, color) in [(src_cx, src_cy, src_r, SRC_COLOR),
                                (bkg_cx, bkg_cy, bkg_r, BKG_COLOR)]:
        ax.add_patch(mpatches.Circle((cx, cy), r, edgecolor=color,
                                     facecolor=color, linewidth=2,
                                     alpha=0.15))
        ax.add_patch(mpatches.Circle((cx, cy), r, edgecolor=color,
                                     facecolor="none", linewidth=2,
                                     linestyle="--"))
    ax.plot(src_cx, src_cy, "+", color=SRC_COLOR,
            markersize=12, markeredgewidth=2, zorder=6)
    ax.set_title(title, color=TEXT_COLOR, fontsize=10)


def _plot_ul_bars(ax, stats_A, stats_B, band):
    ax.set_facecolor(BG_DARK)
    ul_A = stats_A["upper_limits"]
    labels = list(ul_A.keys())

    x = np.arange(len(labels))
    w = 0.35
    cr_A = [ul_A[k]["count_rate_ul"] for k in labels]

    bars_A = ax.bar(x - w/2 if stats_B else x, cr_A,
                    width=w if stats_B else 0.6,
                    color="#ff6347", edgecolor="white", alpha=0.85,
                    label="FPM-A")
    if stats_B:
        ul_B = stats_B["upper_limits"]
        cr_B = [ul_B[k]["count_rate_ul"] for k in labels]
        ax.bar(x + w/2, cr_B, width=w, color="#4169e1",
               edgecolor="white", alpha=0.85, label="FPM-B")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=TEXT_COLOR)
    ax.set_ylabel("UL count rate [cts/s]", color=TEXT_COLOR)
    ax.set_title("Upper Limits", color=TEXT_COLOR)
    ax.legend(fontsize=8, facecolor="#222", edgecolor="gray", labelcolor="white")
    ax.grid(axis="y", alpha=0.2)


def _build_text_summary(stats_A, stats_B, obsid, src_ra, src_dec, band):
    lines = [
        f"  NuSTAR Upper Limit Analysis — Summary",
        f"  {'─'*55}",
        f"  ObsID  : {obsid}",
        f"  Source : RA={src_ra:.6f}°  Dec={src_dec:.6f}°",
        f"  Band   : {band.upper()}",
        f"  {'─'*55}",
    ]
    for label, stats in [("FPM-A", stats_A), ("FPM-B", stats_B)]:
        if stats is None:
            continue
        ul = stats["upper_limits"]
        lines += [
            f"\n  {label}:",
            f"    Exposure        : {stats['exposure_s']:.0f} s",
            f"    Source counts   : {stats['src_counts']}",
            f"    Bkg counts      : {stats['bkg_counts']}  (α={stats['alpha']:.4f})",
            f"    Expected bkg    : {stats['expected_bkg']:.2f}",
            f"    Net counts      : {stats['net_counts']:.1f}",
            f"    Li-Ma signif.   : {stats['lima_sig']:.2f}σ",
        ]
        for sig, d in ul.items():
            lines.append(
                f"    UL ({sig:6s})  : {d['counts_ul']:.1f} cts  |  "
                f"{d['count_rate_ul']:.3e} cts/s  |  "
                f"→ WebPIMMS for flux"
            )
    lines += [
        f"\n  {'─'*55}",
        f"  → Convert count rates to flux at: https://cxc.harvard.edu/toolkit/pimms.jsp",
        f"  → Mission: NUSTAR | choose your spectral model + Galactic NH",
    ]
    return "\n".join(lines)

"""
interactive.py
--------------
Interactive matplotlib-based tool for placing the background region.

Usage
-----
Called when the user requests interactive background placement or when
the source is near a chip edge / corner where an auto-annulus would fail.

The user sees the sky image with the source circle already drawn, and
can click to place a background circle anywhere on the detector.
Multiple clicks cycle through candidates; pressing Enter confirms.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Button, Slider, TextBox
from astropy.wcs import WCS


class InteractiveRegionSelector:
    """
    Interactive tool for background region placement.

    Parameters
    ----------
    image        : 2-D array — sky image
    wcs          : astropy.wcs.WCS
    src_cx, src_cy : float — source centre in image pixels
    src_r          : float — source radius in image pixels
    bkg_r          : float — background radius in image pixels (initial guess)
    obsid          : str
    module         : str ('A' or 'B')
    """

    INSTRUCTIONS = (
        "LEFT CLICK  → place background circle\n"
        "RIGHT CLICK → remove last placement\n"
        "SCROLL      → resize background circle\n"
        "ENTER / close window → confirm selection"
    )

    def __init__(self, image, wcs, src_cx, src_cy, src_r,
                 bkg_r, obsid="", module=""):
        self.image   = image
        self.wcs     = wcs
        self.src_cx  = src_cx
        self.src_cy  = src_cy
        self.src_r   = src_r
        self.bkg_r   = bkg_r
        self.obsid   = obsid
        self.module  = module

        self.bkg_cx  = None
        self.bkg_cy  = None
        self.confirmed = False

        self._build_figure()

    # ------------------------------------------------------------------
    def _build_figure(self):
        """Build the matplotlib figure."""
        self.fig, self.ax = plt.subplots(figsize=(10, 9))
        self.fig.canvas.manager.set_window_title(
            f"NuSTAR UL — Background Region Selector  |  ObsID {self.obsid}  FPM{self.module}"
        )
        self.fig.patch.set_facecolor("#1a1a2e")
        self.ax.set_facecolor("#0f0f1a")

        # Image: log-scale with clipping
        img = np.log1p(self.image)
        vmin, vmax = np.nanpercentile(img, [2, 99])
        self.im = self.ax.imshow(
            img, origin="lower", cmap="afmhot",
            vmin=vmin, vmax=vmax,
            interpolation="nearest",
        )

        # Source circle (fixed, magenta)
        self._src_circle = mpatches.Circle(
            (self.src_cx, self.src_cy), self.src_r,
            edgecolor="magenta", facecolor="none",
            linewidth=2, linestyle="-", label="Source",
            zorder=5,
        )
        self.ax.add_patch(self._src_circle)
        self.ax.plot(self.src_cx, self.src_cy, "+", color="magenta",
                     markersize=12, markeredgewidth=2, zorder=6)

        # Background circle (will be drawn on click)
        self._bkg_patch  = None
        self._bkg_marker = None

        # Instruction text
        self.ax.set_title(
            self.INSTRUCTIONS,
            fontsize=9, color="white", loc="left",
            fontfamily="monospace", pad=10,
        )

        # Status box
        self._status_text = self.ax.text(
            0.98, 0.02,
            "No background region placed yet.",
            transform=self.ax.transAxes,
            color="yellow", fontsize=9, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="#111111", ec="gray", alpha=0.8),
        )

        # Colourbar
        cbar = self.fig.colorbar(self.im, ax=self.ax, pad=0.01, shrink=0.8)
        cbar.set_label("log(counts+1)", color="white", fontsize=9)
        cbar.ax.yaxis.set_tick_params(color="white")
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

        self.ax.tick_params(colors="white")
        for spine in self.ax.spines.values():
            spine.set_edgecolor("gray")

        # Legend
        src_patch = mpatches.Patch(edgecolor="magenta", facecolor="none",
                                   linewidth=2, label="Source")
        bkg_patch = mpatches.Patch(edgecolor="cyan",    facecolor="none",
                                   linewidth=2, label="Background (to place)")
        self.ax.legend(handles=[src_patch, bkg_patch],
                       loc="upper right", fontsize=8,
                       facecolor="#222", edgecolor="gray",
                       labelcolor="white")

        # Tight layout
        self.fig.tight_layout()

        # Connect events
        self._cid_click  = self.fig.canvas.mpl_connect(
            "button_press_event", self._on_click)
        self._cid_scroll = self.fig.canvas.mpl_connect(
            "scroll_event", self._on_scroll)
        self._cid_key    = self.fig.canvas.mpl_connect(
            "key_press_event", self._on_key)
        self._cid_close  = self.fig.canvas.mpl_connect(
            "close_event", self._on_close)

    # ------------------------------------------------------------------
    def _draw_bkg_circle(self):
        """Redraw background circle at current position."""
        if self._bkg_patch is not None:
            self._bkg_patch.remove()
        if self._bkg_marker is not None:
            self._bkg_marker.remove()
            self._bkg_marker = None

        self._bkg_patch = mpatches.Circle(
            (self.bkg_cx, self.bkg_cy), self.bkg_r,
            edgecolor="cyan", facecolor="cyan",
            linewidth=2, linestyle="--", alpha=0.15,
            zorder=4,
        )
        # Outline
        outline = mpatches.Circle(
            (self.bkg_cx, self.bkg_cy), self.bkg_r,
            edgecolor="cyan", facecolor="none",
            linewidth=2, linestyle="--",
            zorder=5,
        )
        self.ax.add_patch(self._bkg_patch)
        self.ax.add_patch(outline)
        self._bkg_patch = outline   # track the outline for removal

        m, = self.ax.plot(self.bkg_cx, self.bkg_cy, "x",
                          color="cyan", markersize=10, markeredgewidth=2, zorder=6)
        self._bkg_marker = m

        # Update WCS label if available
        try:
            ra_bkg, dec_bkg = self.wcs.all_pix2world([[self.bkg_cx, self.bkg_cy]], 0)[0]
            status = (
                f"Background:  ({self.bkg_cx:.1f}, {self.bkg_cy:.1f}) pix\n"
                f"             RA={ra_bkg:.5f}°  Dec={dec_bkg:.5f}°\n"
                f"             r = {self.bkg_r:.1f} pix  |  Scroll to resize"
            )
        except Exception:
            status = f"Background: ({self.bkg_cx:.1f}, {self.bkg_cy:.1f}) pix, r={self.bkg_r:.1f}"

        self._status_text.set_text(status)
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------------
    def _on_click(self, event):
        if event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        if event.button == 1:   # left click → place
            self.bkg_cx = event.xdata
            self.bkg_cy = event.ydata
            self._draw_bkg_circle()
        elif event.button == 3: # right click → clear
            self.bkg_cx = None
            self.bkg_cy = None
            if self._bkg_patch is not None:
                self._bkg_patch.remove()
                self._bkg_patch = None
            if self._bkg_marker is not None:
                self._bkg_marker.remove()
                self._bkg_marker = None
            self._status_text.set_text("Background cleared. Left-click to place.")
            self.fig.canvas.draw_idle()

    def _on_scroll(self, event):
        if event.inaxes != self.ax:
            return
        factor = 1.1 if event.button == "up" else 0.9
        self.bkg_r = max(5, self.bkg_r * factor)
        if self.bkg_cx is not None:
            self._draw_bkg_circle()

    def _on_key(self, event):
        if event.key in ("enter", "return"):
            self._confirm()

    def _on_close(self, event):
        self.confirmed = self.bkg_cx is not None

    def _confirm(self):
        if self.bkg_cx is None:
            print("\n  [WARNING] No background region placed yet!")
            return
        self.confirmed = True
        plt.close(self.fig)

    # ------------------------------------------------------------------
    def run(self):
        """
        Show the interactive window and block until user closes it.

        Returns
        -------
        bkg_cx, bkg_cy, bkg_r : floats — background centre (pixels) and radius
        confirmed              : bool
        """
        print("\n" + "─"*60)
        print("  INTERACTIVE BACKGROUND REGION SELECTOR")
        print("─"*60)
        print(f"  Source at pixel ({self.src_cx:.1f}, {self.src_cy:.1f}), r={self.src_r:.1f} px")
        print(f"  Default background radius: {self.bkg_r:.1f} px")
        print()
        print("  Instructions:")
        for line in self.INSTRUCTIONS.split("\n"):
            print(f"    {line}")
        print("─"*60 + "\n")

        plt.show(block=True)

        if not self.confirmed or self.bkg_cx is None:
            print("\n  No background region confirmed — will use auto-annulus fallback.")
            return None, None, self.bkg_r, False

        print(f"\n  ✓ Background confirmed at pixel ({self.bkg_cx:.1f}, {self.bkg_cy:.1f}), "
              f"r = {self.bkg_r:.1f} px")
        return self.bkg_cx, self.bkg_cy, self.bkg_r, True


# ---------------------------------------------------------------------------
# Convenience function for external use
# ---------------------------------------------------------------------------

def select_background_interactively(image, wcs, src_cx, src_cy,
                                    src_r_pix, bkg_r_pix,
                                    obsid="", module=""):
    """
    Launch interactive background selector and return RA/Dec of chosen region.

    Returns
    -------
    bkg_ra, bkg_dec, bkg_r_pix, confirmed
    If not confirmed, returns (None, None, bkg_r_pix, False).
    """
    sel = InteractiveRegionSelector(
        image, wcs, src_cx, src_cy, src_r_pix,
        bkg_r_pix, obsid=obsid, module=module,
    )
    bkg_cx, bkg_cy, bkg_r, confirmed = sel.run()

    if not confirmed:
        return None, None, bkg_r_pix, False

    # Convert pixel → RA/Dec
    bkg_ra, bkg_dec = wcs.all_pix2world([[bkg_cx, bkg_cy]], 0)[0]
    return float(bkg_ra), float(bkg_dec), float(bkg_r), True

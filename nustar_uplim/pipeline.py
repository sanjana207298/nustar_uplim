"""
pipeline.py
-----------
End-to-end analysis pipeline for NuSTAR non-detection upper limits.

Orchestrates:
  1. File discovery and loading
  2. WCS setup
  3. Source position conversion
  4. Interactive/auto background region selection
  5. Count extraction (source + background)
  6. Poisson upper limit computation  → count-rate upper limits
  7. WebPIMMS instructions printed for flux conversion
  8. Step-by-step visualisation
  9. Text + JSON report output
"""

import os
import json
import numpy as np
from datetime import datetime

from .io          import (find_event_cl_dir, load_events, load_sky_image,
                           get_events_wcs, summary_table)
from .coords      import parse_radec, arcsec_to_pixels, radec_to_pixels
from .extraction  import extract_counts, ENERGY_BANDS, parse_band, load_exposure_map
from .statistics  import compute_upper_limits, format_pimms_instructions
from .visualization import (plot_sky_image, plot_regions, plot_event_scatter,
                             plot_upper_limit_summary, plot_radial_profile,
                             plot_combined_report)
from .interactive import select_background_interactively


MODULES = ("A", "B")


# ──────────────────────────────────────────────────────────────────────────

class NuSTARUpperLimitPipeline:
    """
    Main pipeline class.

    Parameters
    ----------
    base_path      : str  — parent directory containing the obsid folder
    obsid          : str  — NuSTAR observation ID
    ra             : str or float  — source RA (any format)
    dec            : str or float  — source Dec (any format)
    src_radius_as  : float — source region radius in arcseconds (default 30")
    bkg_radius_as  : float — background region radius in arcseconds (default 90")
    energy_band    : str  — 'full', 'soft', 'hard', 'ultrahard' (default 'soft')
    bkg_mode       : str  — 'interactive', 'annulus', or 'auto'
    modules        : list of str, e.g. ['A', 'B']
    confidence_levels : tuple — default (0.9545, 0.9973) for 2σ, 3σ
    out_dir        : str or None — save plots/reports here
    show_plots     : bool — display plots interactively
    stat_method    : str  — 'background_inclusive' | 'frequentist' | 'bayesian'

    Note: Flux conversion is intentionally omitted. Use the output count-rate
    upper limits with WebPIMMS (https://cxc.harvard.edu/toolkit/pimms.jsp)
    and your preferred spectral model.
    """

    def __init__(self,
                 base_path,
                 obsid,
                 ra,
                 dec,
                 src_radius_as=20.0,
                 bkg_radius_as=100.0,
                 energy_band="soft",  # 3-10 keV default
                 bkg_mode="auto",
                 modules=("A", "B"),
                 confidence_levels=(0.9545, 0.9973),
                 out_dir=None,
                 show_plots=True,
                 stat_method="background_inclusive"):

        self.base_path     = str(base_path)
        self.obsid         = str(obsid).strip()
        self.ra_input      = ra
        self.dec_input     = dec
        self.src_radius_as = float(src_radius_as)
        self.bkg_radius_as = float(bkg_radius_as)
        self.energy_band   = energy_band
        self.bkg_mode      = bkg_mode
        self.modules       = list(modules)
        self.confidence_levels = confidence_levels
        self.out_dir       = out_dir
        self.show_plots    = show_plots
        self.stat_method   = stat_method

        self.results       = {}

    # ── Run ─────────────────────────────────────────────────────────────

    def run(self):
        """Execute the full analysis pipeline."""
        print("\n" + "═"*60)
        print("  NuSTAR Non-Detection Upper Limit Pipeline")
        print("  nustar_uplim  v1.0")
        print("═"*60)

        self._parse_coordinates()

        self.event_cl_dir = find_event_cl_dir(self.base_path, self.obsid)
        summary_table(self.event_cl_dir, self.obsid)

        per_module = {}
        images, wcss = {}, {}

        for mod in self.modules:
            print(f"\n{'─'*60}")
            print(f"  Processing FPM-{mod}")
            print(f"{'─'*60}")
            try:
                result = self._process_module(mod)
                per_module[mod] = result
                images[mod] = result["image"]
                wcss[mod]   = result["wcs"]
            except FileNotFoundError as e:
                print(f"  [SKIP] FPM-{mod}: {e}")

        if not per_module:
            raise RuntimeError("No modules could be processed. Check file paths.")

        self._plot_combined(per_module, images, wcss)
        self._write_report(per_module)

        # Print WebPIMMS guidance for all modules + confidence levels
        self._print_pimms_block(per_module)

        self.results = per_module
        return per_module

    # ── Coordinate parsing ───────────────────────────────────────────────

    def _parse_coordinates(self):
        print(f"\n  Parsing coordinates: RA='{self.ra_input}', Dec='{self.dec_input}'")
        self.src_ra, self.src_dec = parse_radec(self.ra_input, self.dec_input)
        print(f"  → RA = {self.src_ra:.6f}°,  Dec = {self.src_dec:.6f}°")

    # ── Per-module processing ────────────────────────────────────────────

    def _process_module(self, mod):
        events, ev_hdr, exposure = load_events(self.event_cl_dir, self.obsid, mod)
        image, wcs, img_hdr      = load_sky_image(self.event_cl_dir, self.obsid, mod)

        try:
            ev_wcs, _, _ = get_events_wcs(ev_hdr)
        except Exception:
            ev_wcs = wcs

        src_cx, src_cy = radec_to_pixels(self.src_ra, self.src_dec, wcs)
        src_r_pix      = arcsec_to_pixels(self.src_radius_as, wcs)
        bkg_r_pix      = arcsec_to_pixels(self.bkg_radius_as, wcs)

        print(f"  Source pixel:   ({src_cx:.1f}, {src_cy:.1f})")
        print(f"  Source radius:  {src_r_pix:.1f} px ({self.src_radius_as:.1f}\")")
        print(f"  Bkg radius:     {bkg_r_pix:.1f} px ({self.bkg_radius_as:.1f}\")")

        # Step 1: sky image
        print(f"\n  [Step 1] Plotting sky image...")
        plot_sky_image(
            image, wcs, self.src_ra, self.src_dec,
            src_r_pix, src_cx, src_cy,
            obsid=self.obsid, module=mod, band=self.energy_band,
            out_dir=self.out_dir, show=self.show_plots,
        )

        # Background region
        bkg_ra, bkg_dec, bkg_r_pix, bkg_type, bkg_inner = \
            self._select_background(mod, image, wcs, src_cx, src_cy,
                                    src_r_pix, bkg_r_pix)

        # Step 2: regions overlay
        print(f"\n  [Step 2] Plotting regions...")
        bkg_cx, bkg_cy = radec_to_pixels(bkg_ra, bkg_dec, wcs)
        plot_regions(
            image, wcs,
            src_cx, src_cy, src_r_pix,
            bkg_cx, bkg_cy, bkg_r_pix,
            bkg_type=bkg_type, bkg_inner_pix=bkg_inner,
            obsid=self.obsid, module=mod,
            out_dir=self.out_dir, show=self.show_plots,
        )

        # Try to load optional exposure map (vignetting correction)
        exp_map, exp_wcs_obj = load_exposure_map(
            self.event_cl_dir, self.obsid, mod)
        if exp_map is None:
            print("  Exposure map not found — skipping vignetting correction.")
            print("  (Run nuexpomap in HEASoft to generate it if needed.)")

        # Extract counts
        band_info = parse_band(self.energy_band)
        print(f"\n  [Step 3] Extracting counts — {band_info[4]}  "
              f"(PI {band_info[0]}–{band_info[1]})...")
        extr = extract_counts(
            events, ev_wcs,
            self.src_ra, self.src_dec, src_r_pix,
            bkg_ra, bkg_dec, bkg_r_pix,
            bkg_type=bkg_type, bkg_inner_pix=bkg_inner,
            energy_band=self.energy_band,
            exp_map=exp_map, exp_wcs=exp_wcs_obj,
        )
        print(f"  → Source counts:   {extr['src_counts']}")
        print(f"  → Background cts:  {extr['bkg_counts']}  (α={extr['alpha']:.4f})")
        print(f"  → Expected bkg:    {extr['expected_bkg']:.2f}")
        print(f"  → Net counts:      {extr['net_counts']:.2f}")

        # Event scatter
        print(f"\n  [Step 3b] Plotting event scatter...")
        plot_event_scatter(
            extr["ev_x"], extr["ev_y"],
            extr["src_mask"], extr["bkg_mask"],
            src_cx, src_cy, src_r_pix,
            bkg_cx, bkg_cy, bkg_r_pix, bkg_type,
            bkg_inner_pix=bkg_inner,
            obsid=self.obsid, module=mod, band=self.energy_band,
            out_dir=self.out_dir, show=self.show_plots,
        )

        # Radial profile
        print(f"\n  [Step 4] Radial profile...")
        plot_radial_profile(
            extr["ev_x"], extr["ev_y"], src_cx, src_cy,
            r_max=src_r_pix * 3,
            obsid=self.obsid, module=mod, band=self.energy_band,
            out_dir=self.out_dir, show=self.show_plots,
        )

        # Compute upper limits (count rates only)
        print(f"\n  [Step 5] Computing Poisson upper limits...")
        stats = compute_upper_limits(
            extr, exposure,
            band=self.energy_band,
            confidence_levels=self.confidence_levels,
            method=self.stat_method,
        )

        # UL plot
        plot_upper_limit_summary(
            stats, obsid=self.obsid, band=self.energy_band,
            out_dir=self.out_dir, show=self.show_plots,
        )

        self._print_module_summary(mod, stats)

        return {
            "events":    events,
            "image":     image,
            "wcs":       wcs,
            "exposure":  exposure,
            "extraction": extr,
            "stats":     stats,
            "src_cx": src_cx, "src_cy": src_cy, "src_r": src_r_pix,
            "bkg_cx": bkg_cx, "bkg_cy": bkg_cy, "bkg_r": bkg_r_pix,
            "bkg_type": bkg_type, "bkg_inner": bkg_inner,
        }

    # ── Background selection ─────────────────────────────────────────────

    def _select_background(self, mod, image, wcs, src_cx, src_cy,
                           src_r_pix, bkg_r_pix):
        mode = self.bkg_mode.lower()

        if mode == "annulus":
            print("  Background mode: ANNULUS around source")
            bkg_inner = src_r_pix
            bkg_outer = bkg_r_pix if bkg_r_pix > bkg_inner else bkg_inner * 3
            return (self.src_ra, self.src_dec,
                    bkg_outer, "annulus", bkg_inner)

        if mode == "interactive":
            return self._interactive_background(
                mod, image, wcs, src_cx, src_cy, src_r_pix, bkg_r_pix)

        # auto
        h, w = image.shape
        margin = bkg_r_pix * 2
        near_edge = (src_cx < margin or src_cx > w - margin or
                     src_cy < margin or src_cy > h - margin)

        if near_edge:
            print("  Source is near image edge — launching interactive selector.")
            return self._interactive_background(
                mod, image, wcs, src_cx, src_cy, src_r_pix, bkg_r_pix)
        else:
            print("  Background mode: AUTO-ANNULUS around source")
            bkg_inner = src_r_pix * 1.2
            bkg_outer = max(bkg_r_pix, bkg_inner * 2.5)
            return (self.src_ra, self.src_dec,
                    bkg_outer, "annulus", bkg_inner)

    def _interactive_background(self, mod, image, wcs,
                                 src_cx, src_cy, src_r_pix, bkg_r_pix):
        bkg_ra, bkg_dec, bkg_r, confirmed = select_background_interactively(
            image, wcs, src_cx, src_cy, src_r_pix, bkg_r_pix,
            obsid=self.obsid, module=mod,
        )
        if not confirmed:
            print("  Falling back to annulus background.")
            bkg_inner = src_r_pix * 1.2
            bkg_outer = max(bkg_r_pix, bkg_inner * 2.5)
            return (self.src_ra, self.src_dec,
                    bkg_outer, "annulus", bkg_inner)
        return bkg_ra, bkg_dec, bkg_r, "circle", None

    # ── Combined report figure ───────────────────────────────────────────

    def _plot_combined(self, per_module, images, wcss):
        print(f"\n  Generating combined report figure...")
        mods = list(per_module.keys())
        r_A  = per_module[mods[0]]
        r_B  = per_module[mods[1]] if len(mods) > 1 else None

        plot_combined_report(
            images[mods[0]], wcss[mods[0]],
            r_A["src_cx"], r_A["src_cy"], r_A["src_r"],
            r_A["bkg_cx"], r_A["bkg_cy"], r_A["bkg_r"],
            r_A["stats"],
            image_B=images[mods[1]] if r_B else None,
            wcs_B=wcss[mods[1]] if r_B else None,
            src_cx_B=r_B["src_cx"] if r_B else None,
            src_cy_B=r_B["src_cy"] if r_B else None,
            src_r_B=r_B["src_r"]   if r_B else None,
            bkg_cx_B=r_B["bkg_cx"] if r_B else None,
            bkg_cy_B=r_B["bkg_cy"] if r_B else None,
            bkg_r_B=r_B["bkg_r"]   if r_B else None,
            stats_B=r_B["stats"]   if r_B else None,
            obsid=self.obsid,
            src_ra=self.src_ra, src_dec=self.src_dec,
            band=self.energy_band,
            out_dir=self.out_dir, show=self.show_plots,
        )

    # ── Printing ─────────────────────────────────────────────────────────

    def _print_module_summary(self, mod, stats):
        print(f"\n  ┌─ FPM-{mod} Results ─────────────────────────────────┐")
        print(f"  │  Exposure:        {stats['exposure_s']:.0f} s")
        print(f"  │  Source counts:   {stats['src_counts']}")
        print(f"  │  Background:      {stats['bkg_counts']} cts  (α={stats['alpha']:.4f})")
        print(f"  │  Net counts:      {stats['net_counts']:.1f}")
        print(f"  │  Li-Ma signif:    {stats['lima_sig']:.2f}σ")
        for sig, d in stats["upper_limits"].items():
            print(f"  │  UL ({sig:6s}):    "
                  f"{d['counts_ul']:.2f} cts  |  "
                  f"{d['count_rate_ul']:.4e} cts/s  ← use in WebPIMMS")
        print(f"  └─────────────────────────────────────────────────────┘")

    def _print_pimms_block(self, per_module):
        """Print WebPIMMS instructions for each module and confidence level."""
        print("\n" + "═"*60)
        print("  FLUX CONVERSION via WebPIMMS")
        print("  https://cxc.harvard.edu/toolkit/pimms.jsp")
        print("═"*60)
        print("  Take the count-rate upper limit below and enter it into")
        print("  WebPIMMS with your spectral model (power law, BB, APEC...)")
        print("  and Galactic NH to get the flux upper limit.\n")

        for mod, res in per_module.items():
            ul = res["stats"]["upper_limits"]
            # Print instructions for the most commonly used level (3σ)
            preferred = list(ul.keys())[-1]   # last = highest sigma
            cr = ul[preferred]["count_rate_ul"]
            print(format_pimms_instructions(
                cr, self.energy_band,
                obsid=self.obsid, module=mod,
                sigma_label=preferred,
            ))

            # Also print all levels in a compact table
            print(f"  FPM-{mod} — all confidence levels:")
            print(f"  {'Level':<10} {'Counts UL':>12} {'Count Rate UL':>16}")
            print(f"  {'─'*40}")
            for sig, d in ul.items():
                print(f"  {sig:<10} {d['counts_ul']:>12.2f} {d['count_rate_ul']:>16.4e} cts/s")
            print()

    # ── Report files ──────────────────────────────────────────────────────

    def _write_report(self, per_module):
        if self.out_dir is None:
            return

        os.makedirs(self.out_dir, exist_ok=True)

        # JSON
        report_data = {
            "obsid":        self.obsid,
            "src_ra_deg":   self.src_ra,
            "src_dec_deg":  self.src_dec,
            "energy_band":  self.energy_band,
            "stat_method":  self.stat_method,
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "pimms_url":    "https://cxc.harvard.edu/toolkit/pimms.jsp",
            "note":         (
                "Use count_rate_ul values in WebPIMMS with Mission=NuSTAR, "
                "your spectral model and Galactic NH to obtain flux upper limits."
            ),
            "modules": {},
        }
        for mod, res in per_module.items():
            s  = res["stats"]
            ul = s["upper_limits"]
            report_data["modules"][f"FPM{mod}"] = {
                "exposure_s":   s["exposure_s"],
                "src_counts":   s["src_counts"],
                "bkg_counts":   s["bkg_counts"],
                "alpha":        s["alpha"],
                "expected_bkg": s["expected_bkg"],
                "net_counts":   s["net_counts"],
                "lima_sigma":   s["lima_sig"],
                "upper_limits": {
                    k: {
                        "confidence":    v["confidence"],
                        "counts_ul":     v["counts_ul"],
                        "count_rate_ul_cts_s": v["count_rate_ul"],
                        "webpimms_note": (
                            f"Enter {v['count_rate_ul']:.4e} cts/s as input "
                            f"count rate in WebPIMMS with Mission=NUSTAR FPM{mod}."
                        ),
                    }
                    for k, v in ul.items()
                },
            }

        json_path = os.path.join(self.out_dir, "nustar_uplim_results.json")
        with open(json_path, "w") as f:
            json.dump(report_data, f, indent=2)
        print(f"\n  Results saved → {json_path}")

        # Plain text
        txt_path = os.path.join(self.out_dir, "nustar_uplim_results.txt")
        with open(txt_path, "w") as f:
            f.write("NuSTAR Non-Detection Upper Limit Analysis\n")
            f.write("=" * 60 + "\n")
            f.write(f"ObsID:         {self.obsid}\n")
            f.write(f"Source RA:     {self.src_ra:.6f} deg\n")
            f.write(f"Source Dec:    {self.src_dec:.6f} deg\n")
            f.write(f"Energy band:   {self.energy_band.upper()}\n")
            f.write(f"Stat method:   {self.stat_method}\n")
            f.write(f"Generated:     {datetime.utcnow().isoformat()}Z\n")
            f.write("=" * 60 + "\n")
            f.write("\nFLUX CONVERSION: Use count_rate_ul in WebPIMMS\n")
            f.write("  https://cxc.harvard.edu/toolkit/pimms.jsp\n")
            f.write("  Mission: NUSTAR, choose your spectral model + NH\n\n")
            for mod, res in per_module.items():
                s  = res["stats"]
                ul = s["upper_limits"]
                f.write(f"FPM-{mod}\n")
                f.write(f"  Exposure:    {s['exposure_s']:.0f} s\n")
                f.write(f"  Src counts:  {s['src_counts']}\n")
                f.write(f"  Bkg counts:  {s['bkg_counts']}  (alpha={s['alpha']:.4f})\n")
                f.write(f"  Net counts:  {s['net_counts']:.1f}\n")
                f.write(f"  Li-Ma sig:   {s['lima_sig']:.2f} sigma\n")
                for sig, d in ul.items():
                    f.write(
                        f"  UL({sig:6s}): {d['counts_ul']:.2f} cts  "
                        f"{d['count_rate_ul']:.4e} cts/s\n"
                    )
                f.write("\n")
        print(f"  Report saved  → {txt_path}")

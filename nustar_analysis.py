#!/usr/bin/env python3
"""
nustar_analysis.py
------------------
Command-line entry point for the NuSTAR upper limit pipeline.

Usage examples
--------------
# Basic (interactive background, show all plots)
python nustar_analysis.py \\
    --path /data/nustar \\
    --obsid 90601606002 \\
    --ra "05 00 13.72" \\
    --dec "-03 20 51.20"

# Specify source and background radii (arcseconds)
python nustar_analysis.py \\
    --path /data/nustar \\
    --obsid 90601606002 \\
    --ra "05h00m13.72s" \\
    --dec "-03d20m51.2s" \\
    --src-radius 30 \\
    --bkg-radius 90 \\
    --band soft \\
    --bkg-mode annulus \\
    --out-dir ./results

# Decimal degrees RA/Dec, Bayesian upper limit, known distance
python nustar_analysis.py \\
    --path . --obsid 90601606002 \\
    --ra 75.057167 --dec -3.347556 \\
    --src-radius 25 --bkg-radius 80 \\
    --distance-mpc 10.0 \\
    --stat-method bayesian \\
    --no-show-plots --out-dir ./results

# Interactive mode (asks all questions)
python nustar_analysis.py --interactive
"""

import argparse
import sys
import os


# ── Banner ─────────────────────────────────────────────────────────────────

BANNER = r"""
  ╔══════════════════════════════════════════════════════════╗
  ║   nustar_uplim  ·  NuSTAR Non-Detection Upper Limits    ║
  ║   v1.0.0   github.com/your-org/nustar_uplim             ║
  ╚══════════════════════════════════════════════════════════╝
"""


# ── Argument parser ─────────────────────────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(
        prog="nustar_analysis",
        description="Compute NuSTAR count-rate and flux upper limits for "
                    "non-detected sources.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Required (or interactive)
    p.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Prompt for all parameters interactively (overrides other flags).",
    )
    p.add_argument("--path",   type=str, default=None,
                   help="Parent directory containing the ObsID folder.")
    p.add_argument("--obsid",  type=str, default=None,
                   help="NuSTAR ObsID (e.g. 90601606002).")
    p.add_argument("--ra",     type=str, default=None,
                   help="Source RA — any format: '05 00 13.72', '05:00:13.72', "
                        "'05h00m13.72s', or decimal degrees.")
    p.add_argument("--dec",    type=str, default=None,
                   help="Source Dec — any format: '-03 20 51.20', '-03:20:51.20', "
                        "'-03d20m51.2s', or decimal degrees.")

    # Region parameters
    p.add_argument("--src-radius", type=float, default=30.0,
                   help="Source region radius in arcseconds (default: 30).")
    p.add_argument("--bkg-radius", type=float, default=90.0,
                   help="Background region radius in arcseconds (default: 90).")
    p.add_argument("--bkg-mode", type=str, default="auto",
                   choices=["auto", "interactive", "annulus"],
                   help="Background region selection: 'auto' (annulus unless near "
                        "edge), 'interactive' (click on image), 'annulus' "
                        "(always use annulus). Default: auto.")

    # Energy band
    p.add_argument("--band", type=str, default="soft",
                   choices=["full", "soft", "hard", "ultrahard"],
                   help="Energy band: full (3–79 keV), soft (3–10 keV), "
                        "hard (10–30 keV), ultrahard (30–79 keV). Default: soft.")

    # Modules
    p.add_argument("--modules", type=str, nargs="+", default=["A", "B"],
                   choices=["A", "B"],
                   help="Focal-plane modules to process (default: A B).")

    # Statistical options
    p.add_argument("--stat-method", type=str, default="background_inclusive",
                   choices=["background_inclusive", "frequentist", "bayesian"],
                   help="Upper limit method (default: background_inclusive).")
    p.add_argument("--confidence", type=float, nargs="+",
                   default=[0.9545, 0.9973],
                   help="Confidence levels for upper limits (default: 0.9545 0.9973 "
                        "corresponding to 2σ and 3σ).")

    # Output
    p.add_argument("--out-dir", type=str, default=None,
                   help="Directory to save plots and report files.")
    p.add_argument("--no-show-plots", action="store_true",
                   help="Do not display interactive plots (useful for batch runs).")

    return p


# ── Interactive prompts ─────────────────────────────────────────────────────

def prompt(msg, default=None, cast=str, choices=None):
    """Simple input prompt with default and optional choice validation."""
    while True:
        if default is not None:
            full_msg = f"  {msg} [{default}]: "
        else:
            full_msg = f"  {msg}: "
        val = input(full_msg).strip()
        if val == "" and default is not None:
            val = default
        if val == "":
            print("  (required — please enter a value)")
            continue
        if choices and val not in choices:
            print(f"  Must be one of: {choices}")
            continue
        try:
            return cast(val)
        except (ValueError, TypeError):
            print(f"  Could not parse as {cast.__name__}, please retry.")


def interactive_setup():
    """Gather all parameters interactively."""
    print("\n  Fill in the analysis parameters (press Enter to accept defaults):\n")

    path       = prompt("Path to NuSTAR data directory")
    obsid      = prompt("ObsID (e.g. 90601606002)")
    ra         = prompt("Source RA  (e.g. '05 00 13.72' or '75.057167')")
    dec        = prompt("Source Dec (e.g. '-03 20 51.20' or '-3.348')")
    src_r      = prompt("Source region radius (arcsec)", default="30", cast=float)
    bkg_r      = prompt("Background region radius (arcsec)", default="90", cast=float)
    bkg_mode   = prompt("Background mode (auto/interactive/annulus)",
                        default="auto", choices=["auto", "interactive", "annulus"])
    band       = prompt("Energy band (full/soft/hard/ultrahard)", default="soft",
                        choices=["full", "soft", "hard", "ultrahard"])
    stat_meth  = prompt("Upper limit method (background_inclusive/frequentist/bayesian)",
                        default="background_inclusive",
                        choices=["background_inclusive", "frequentist", "bayesian"])

    out_dir    = prompt("Output directory (leave blank for no saving)", default="")
    out_dir    = out_dir if out_dir.strip() != "" else None

    show_str   = prompt("Show plots interactively? (yes/no)", default="yes",
                        choices=["yes", "no"])
    show_plots = show_str.lower() == "yes"

    return dict(
        base_path=path, obsid=obsid, ra=ra, dec=dec,
        src_radius_as=src_r, bkg_radius_as=bkg_r,
        bkg_mode=bkg_mode, energy_band=band,
        stat_method=stat_meth, out_dir=out_dir, show_plots=show_plots,
        confidence_levels=(0.9545, 0.9973),
        modules=["A", "B"],
    )


# ── Main ────────────────────────────────────────────────────────────────────

def main(argv=None):
    print(BANNER)

    parser = build_parser()
    args   = parser.parse_args(argv)

    # --- Interactive mode ---
    if args.interactive:
        kwargs = interactive_setup()
    else:
        # Validate required args
        missing = [f for f, v in [("--path", args.path),
                                   ("--obsid", args.obsid),
                                   ("--ra", args.ra),
                                   ("--dec", args.dec)]
                   if v is None]
        if missing:
            print(f"\n  ERROR: Missing required arguments: {', '.join(missing)}")
            print("  Use --interactive for guided input, or --help for usage.\n")
            parser.print_help()
            sys.exit(1)

        kwargs = dict(
            base_path=args.path,
            obsid=args.obsid,
            ra=args.ra,
            dec=args.dec,
            src_radius_as=args.src_radius,
            bkg_radius_as=args.bkg_radius,
            bkg_mode=args.bkg_mode,
            energy_band=args.band,
            modules=args.modules,
            confidence_levels=tuple(args.confidence),
            out_dir=args.out_dir,
            show_plots=not args.no_show_plots,
            stat_method=args.stat_method,
        )

    # --- Run pipeline ---
    from nustar_uplim import NuSTARUpperLimitPipeline

    pipeline = NuSTARUpperLimitPipeline(**kwargs)

    try:
        results = pipeline.run()
    except Exception as e:
        print(f"\n  FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(2)

    # --- Final summary ---
    print("\n" + "═"*60)
    print("  Analysis complete!")
    if kwargs.get("out_dir"):
        print(f"  All outputs saved to: {kwargs['out_dir']}")
    print("═"*60 + "\n")

    return results


if __name__ == "__main__":
    main()

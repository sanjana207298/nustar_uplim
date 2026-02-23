#!/usr/bin/env python3
"""
nustar_analysis.py
------------------
CLI entry point for the NuSTAR upper limit pipeline.

Usage examples
--------------
# Named band
python nustar_analysis.py \\
    --path /data/nustar --obsid 90601606002 \\
    --ra "05 00 13.72" --dec "-03 20 51.20" \\
    --band soft

# Custom keV range
python nustar_analysis.py \\
    --path /data/nustar --obsid 90601606002 \\
    --ra "05 00 13.72" --dec "-03 20 51.20" \\
    --band 8-30

# Non-default radii, annulus background, save output
python nustar_analysis.py \\
    --path /data/nustar --obsid 90601606002 \\
    --ra "05 00 13.72" --dec "-03 20 51.20" \\
    --src-radius 25 --bkg-radius 120 \\
    --band 3-10 --bkg-mode annulus \\
    --out-dir ./results

# Interactive mode (asks everything)
python nustar_analysis.py --interactive
"""

import argparse
import sys
import os

BANNER = r"""
  ╔══════════════════════════════════════════════════════════╗
  ║   nustar_uplim  ·  NuSTAR Non-Detection Upper Limits    ║
  ║   v1.0.0                                                ║
  ╚══════════════════════════════════════════════════════════╝
"""

BAND_HELP = (
    "Energy band. Named: full (3–79 keV), soft (3–10 keV), "
    "hard (10–30 keV), ultrahard (30–79 keV). "
    "Custom keV range: e.g. '8-30', '15:50'. Default: soft (3–10 keV)."
)


def build_parser():
    p = argparse.ArgumentParser(
        prog="nustar_analysis",
        description="Compute NuSTAR count-rate upper limits for non-detected sources.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    p.add_argument("--interactive", "-i", action="store_true",
                   help="Prompt for all parameters interactively.")
    p.add_argument("--path",  type=str, default=None,
                   help="Parent directory containing the ObsID folder.")
    p.add_argument("--obsid", type=str, default=None,
                   help="NuSTAR ObsID (e.g. 90601606002).")
    p.add_argument("--ra",    type=str, default=None,
                   help="Source RA — '05 00 13.72', '05:00:13.72', "
                        "'05h00m13.72s', or decimal degrees.")
    p.add_argument("--dec",   type=str, default=None,
                   help="Source Dec — '-03 20 51.20', '-03:20:51', "
                        "'-03d20m51s', or decimal degrees.")

    p.add_argument("--src-radius", type=float, default=20.0,
                   help="Source region radius in arcseconds (default: 20).")
    p.add_argument("--bkg-radius", type=float, default=100.0,
                   help="Background region radius in arcseconds (default: 100).")
    p.add_argument("--bkg-mode", type=str, default="auto",
                   choices=["auto", "interactive", "annulus"],
                   help="Background mode: auto, interactive, annulus (default: auto).")

    p.add_argument("--band", type=str, default="soft",
                   help=BAND_HELP)

    p.add_argument("--modules", type=str, nargs="+", default=["A", "B"],
                   choices=["A", "B"],
                   help="FPM modules to process (default: A B).")

    p.add_argument("--stat-method", type=str, default="background_inclusive",
                   choices=["background_inclusive", "frequentist", "bayesian"],
                   help="Upper limit method (default: background_inclusive).")
    p.add_argument("--confidence", type=float, nargs="+",
                   default=[0.9545, 0.9973],
                   help="Confidence levels (default: 0.9545 0.9973 = 2σ 3σ).")

    p.add_argument("--out-dir", type=str, default=None,
                   help="Directory to save plots and report files.")
    p.add_argument("--no-show-plots", action="store_true",
                   help="Suppress interactive plot windows.")

    return p


def prompt(msg, default=None, cast=str, choices=None):
    while True:
        suffix = f" [{default}]" if default is not None else ""
        val = input(f"  {msg}{suffix}: ").strip()
        if val == "" and default is not None:
            val = default
        if val == "":
            print("  (required)")
            continue
        if choices and val not in choices:
            print(f"  Must be one of: {choices}")
            continue
        try:
            return cast(val)
        except (ValueError, TypeError):
            print(f"  Could not parse as {cast.__name__}, please retry.")


def interactive_setup():
    print("\n  Fill in parameters (Enter to accept defaults):\n")
    path     = prompt("Path to NuSTAR data directory")
    obsid    = prompt("ObsID (e.g. 90601606002)")
    ra       = prompt("Source RA  (e.g. '05 00 13.72' or '75.057167')")
    dec      = prompt("Source Dec (e.g. '-03 20 51.20' or '-3.348')")
    src_r    = prompt("Source region radius (arcsec)", default="20", cast=float)
    bkg_r    = prompt("Background region radius (arcsec)", default="100", cast=float)
    bkg_mode = prompt("Background mode (auto/interactive/annulus)",
                      default="auto", choices=["auto", "interactive", "annulus"])
    print(f"\n  Band options:")
    print(f"    Named : full (3–79 keV), soft (3–10 keV), hard (10–30 keV), ultrahard (30–79 keV)")
    print(f"    Custom: enter keV range e.g. '8-30' or '15:50'")
    band     = prompt("Energy band", default="soft")
    stat_m   = prompt("Upper limit method (background_inclusive/frequentist/bayesian)",
                      default="background_inclusive",
                      choices=["background_inclusive", "frequentist", "bayesian"])
    out_dir  = prompt("Output directory (blank = no saving)", default="")
    out_dir  = out_dir if out_dir.strip() else None
    show_str = prompt("Show plots interactively? (yes/no)", default="yes",
                      choices=["yes", "no"])

    return dict(
        base_path=path, obsid=obsid, ra=ra, dec=dec,
        src_radius_as=src_r, bkg_radius_as=bkg_r,
        bkg_mode=bkg_mode, energy_band=band,
        stat_method=stat_m, out_dir=out_dir,
        show_plots=(show_str.lower() == "yes"),
        confidence_levels=(0.9545, 0.9973),
        modules=["A", "B"],
    )


def main(argv=None):
    print(BANNER)
    parser = build_parser()
    args   = parser.parse_args(argv)

    if args.interactive:
        kwargs = interactive_setup()
    else:
        missing = [f for f, v in [("--path", args.path), ("--obsid", args.obsid),
                                   ("--ra", args.ra), ("--dec", args.dec)]
                   if v is None]
        if missing:
            print(f"\n  ERROR: Missing required arguments: {', '.join(missing)}")
            print("  Use --interactive for guided input, or --help for usage.\n")
            parser.print_help()
            sys.exit(1)

        kwargs = dict(
            base_path=args.path, obsid=args.obsid, ra=args.ra, dec=args.dec,
            src_radius_as=args.src_radius, bkg_radius_as=args.bkg_radius,
            bkg_mode=args.bkg_mode, energy_band=args.band,
            modules=args.modules, confidence_levels=tuple(args.confidence),
            out_dir=args.out_dir, show_plots=not args.no_show_plots,
            stat_method=args.stat_method,
        )

    # Validate band before running
    from nustar_uplim.extraction import parse_band
    try:
        parsed = parse_band(kwargs["energy_band"])
        print(f"\n  Energy band: {parsed[4]}  (PI {parsed[0]}–{parsed[1]})")
    except ValueError as e:
        print(f"\n  ERROR: {e}")
        sys.exit(1)

    from nustar_uplim import NuSTARUpperLimitPipeline
    pipeline = NuSTARUpperLimitPipeline(**kwargs)

    try:
        results = pipeline.run()
    except Exception as e:
        print(f"\n  FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(2)

    print("\n" + "═"*60)
    print("  Analysis complete!")
    if kwargs.get("out_dir"):
        print(f"  Outputs saved to: {kwargs['out_dir']}")
    print("═"*60 + "\n")
    return results


if __name__ == "__main__":
    main()

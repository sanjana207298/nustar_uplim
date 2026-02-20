"""
nustar_uplim
============
Python package for NuSTAR non-detection upper limit analysis.

Outputs count-rate upper limits which you then convert to flux using
WebPIMMS (https://cxc.harvard.edu/toolkit/pimms.jsp) with your chosen
spectral model and column density.

Quick start
-----------
>>> from nustar_uplim import NuSTARUpperLimitPipeline
>>> pipe = NuSTARUpperLimitPipeline(
...     base_path  = "/data/nustar",
...     obsid      = "90601606002",
...     ra         = "05 00 13.72",
...     dec        = "-03 20 51.20",
...     src_radius_as = 30.0,
...     bkg_radius_as = 90.0,
...     energy_band   = "soft",
...     bkg_mode      = "auto",
...     out_dir       = "./results",
... )
>>> results = pipe.run()
"""

__version__ = "1.0.0"
__author__  = "nustar_uplim contributors"
__license__ = "MIT"

from .pipeline     import NuSTARUpperLimitPipeline
from .coords       import parse_radec
from .statistics   import (compute_upper_limits, lima_significance,
                            background_inclusive_upper_limit,
                            bayesian_upper_limit,
                            frequentist_upper_limit,
                            format_pimms_instructions,
                            ENERGY_BANDS, PIMMS_BAND_INFO, PIMMS_URL)
from .extraction   import extract_counts, filter_events_energy
from .io           import (find_event_cl_dir, load_events,
                            load_sky_image, get_events_wcs)
from .visualization import (plot_sky_image, plot_regions, plot_event_scatter,
                             plot_upper_limit_summary, plot_combined_report)

__all__ = [
    "NuSTARUpperLimitPipeline",
    "parse_radec",
    "compute_upper_limits",
    "lima_significance",
    "background_inclusive_upper_limit",
    "bayesian_upper_limit",
    "frequentist_upper_limit",
    "format_pimms_instructions",
    "extract_counts",
    "filter_events_energy",
    "find_event_cl_dir",
    "load_events",
    "load_sky_image",
    "ENERGY_BANDS",
    "PIMMS_BAND_INFO",
    "PIMMS_URL",
]

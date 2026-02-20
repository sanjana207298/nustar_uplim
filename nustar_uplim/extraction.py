"""
extraction.py
-------------
Extract counts from source and background regions in NuSTAR event files.

Regions are defined as circles in sky (RA/Dec) coordinates. Events are
filtered by their X,Y pixel positions using the WCS transformation.

Energy filtering is also supported (3–79 keV full NuSTAR band, or sub-bands).
"""

import numpy as np
from astropy.wcs import WCS


# NuSTAR standard energy bands (PI channels; 1 PI = 0.04 keV, offset 1.6 keV)
# PI = (E_keV - 1.6) / 0.04
ENERGY_BANDS = {
    "full":    (  35, 1935),   # 3–79 keV
    "soft":    (  35,  210),   # 3–10 keV
    "hard":    ( 210,  710),   # 10–30 keV
    "ultrahard":(710, 1935),   # 30–79 keV
}


def energy_to_pi(e_kev):
    """Convert energy in keV to NuSTAR PI channel."""
    return int((e_kev - 1.6) / 0.04)


def filter_events_energy(events, band="full"):
    """
    Return boolean mask for events within a given energy band.

    Parameters
    ----------
    events : numpy recarray — events table from FITS
    band   : str key from ENERGY_BANDS, or tuple (pi_min, pi_max)

    Returns
    -------
    mask : boolean array, length = len(events)
    """
    if isinstance(band, str):
        pi_min, pi_max = ENERGY_BANDS[band]
    else:
        pi_min, pi_max = band

    pi = events["PI"]
    return (pi >= pi_min) & (pi < pi_max)


def circle_mask(x, y, cx, cy, r):
    """
    Boolean mask for events inside a circle.

    Parameters
    ----------
    x, y   : arrays — event coordinates (pixels)
    cx, cy : float  — circle centre (pixels)
    r      : float  — circle radius (pixels)

    Returns
    -------
    mask : boolean array
    """
    return (x - cx)**2 + (y - cy)**2 <= r**2


def annulus_mask(x, y, cx, cy, r_in, r_out):
    """Boolean mask for an annular region."""
    d2 = (x - cx)**2 + (y - cy)**2
    return (d2 >= r_in**2) & (d2 <= r_out**2)


def region_area_pixels(mask_fn, *args, **kwargs):
    """
    Approximate area of a region in square pixels by evaluating mask_fn
    on a dense grid.  Used for background scaling.

    Parameters
    ----------
    mask_fn : callable — one of circle_mask, annulus_mask, etc.
    *args   : passed to mask_fn (centre, radii)
    kwargs  : grid_size=2000 — size of evaluation grid
    """
    grid_size = kwargs.get("grid_size", 2000)
    # Determine bounding box from args (assume circle or annulus)
    # args for circle:  x_arr, y_arr, cx, cy, r
    # We just use a uniform sampling approach
    cx, cy = args[2], args[3]
    r_max  = args[-1]  # outermost radius
    xs = np.linspace(cx - r_max - 1, cx + r_max + 1, grid_size)
    ys = np.linspace(cy - r_max - 1, cy + r_max + 1, grid_size)
    xg, yg = np.meshgrid(xs, ys)
    mask = mask_fn(xg, yg, *args[2:])
    pixel_size = (xs[1] - xs[0]) * (ys[1] - ys[0])
    return float(np.sum(mask)) * pixel_size


def extract_counts(events, events_wcs,
                   src_ra, src_dec, src_radius_pix,
                   bkg_ra, bkg_dec, bkg_radius_pix,
                   bkg_type="circle",          # 'circle' or 'annulus'
                   bkg_inner_pix=None,         # for annulus only
                   energy_band="full"):
    """
    Extract source and background counts from an events table.

    Parameters
    ----------
    events            : numpy recarray
    events_wcs        : astropy.wcs.WCS  (from sky image or events header)
    src_ra/dec        : float — source position in degrees
    src_radius_pix    : float — source circle radius in pixels
    bkg_ra/dec        : float — background centre in degrees
    bkg_radius_pix    : float — background circle outer radius (or annulus outer)
    bkg_type          : 'circle' or 'annulus'
    bkg_inner_pix     : float — annulus inner radius (pixels); required if annulus
    energy_band       : str or tuple

    Returns
    -------
    result : dict with keys:
        src_counts, bkg_counts, alpha, net_counts,
        src_cx, src_cy, bkg_cx, bkg_cy,
        src_area, bkg_area,
        src_events_mask, bkg_events_mask,
        n_events_total
    """
    # Energy filter
    e_mask = filter_events_energy(events, energy_band)
    ev = events[e_mask]

    x = ev["X"].astype(float)
    y = ev["Y"].astype(float)

    # Convert RA/Dec → pixel coords using WCS
    src_pix = events_wcs.all_world2pix([[src_ra, src_dec]], 0)[0]
    src_cx, src_cy = src_pix[0], src_pix[1]

    bkg_pix = events_wcs.all_world2pix([[bkg_ra, bkg_dec]], 0)[0]
    bkg_cx, bkg_cy = bkg_pix[0], bkg_pix[1]

    # Source mask
    src_mask = circle_mask(x, y, src_cx, src_cy, src_radius_pix)
    src_counts = int(np.sum(src_mask))

    # Background mask
    if bkg_type == "annulus":
        if bkg_inner_pix is None:
            bkg_inner_pix = src_radius_pix
        bkg_mask = annulus_mask(x, y, bkg_cx, bkg_cy, bkg_inner_pix, bkg_radius_pix)
        bkg_area  = np.pi * (bkg_radius_pix**2 - bkg_inner_pix**2)
    else:  # circle
        bkg_mask = circle_mask(x, y, bkg_cx, bkg_cy, bkg_radius_pix)
        bkg_area  = np.pi * bkg_radius_pix**2

    bkg_counts = int(np.sum(bkg_mask))
    src_area   = np.pi * src_radius_pix**2

    # Scaling factor: source area / background area
    alpha = src_area / bkg_area if bkg_area > 0 else 1.0

    # Background-scaled counts
    expected_bkg   = alpha * bkg_counts
    net_counts     = src_counts - expected_bkg

    return {
        "src_counts":      src_counts,
        "bkg_counts":      bkg_counts,
        "alpha":           alpha,
        "expected_bkg":    expected_bkg,
        "net_counts":      net_counts,
        "src_cx":          src_cx,
        "src_cy":          src_cy,
        "bkg_cx":          bkg_cx,
        "bkg_cy":          bkg_cy,
        "src_radius_pix":  src_radius_pix,
        "bkg_radius_pix":  bkg_radius_pix,
        "bkg_inner_pix":   bkg_inner_pix,
        "bkg_type":        bkg_type,
        "src_area":        src_area,
        "bkg_area":        bkg_area,
        "src_mask":        src_mask,
        "bkg_mask":        bkg_mask,
        "n_events_total":  len(ev),
        "energy_band":     energy_band,
        "ev_x":            x,
        "ev_y":            y,
    }

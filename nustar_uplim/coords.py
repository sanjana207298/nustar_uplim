"""
coords.py
---------
Utilities for parsing RA/Dec in multiple formats and converting to degrees.
Supports: '12 34 56', '12:34:56', '12h34m56s', decimal degrees, etc.
"""

from astropy.coordinates import SkyCoord
import astropy.units as u
import re


def parse_radec(ra_input, dec_input):
    """
    Robustly parse RA and Dec from many common string or float formats.

    Parameters
    ----------
    ra_input  : str or float — RA in any common format
    dec_input : str or float — Dec in any common format

    Returns
    -------
    (ra_deg, dec_deg) : tuple of float, both in decimal degrees
    """
    # If already numeric, handle directly
    if isinstance(ra_input, (int, float)) and isinstance(dec_input, (int, float)):
        return float(ra_input), float(dec_input)

    ra_str  = str(ra_input).strip()
    dec_str = str(dec_input).strip()

    # Attempt 1: Sexagesimal / unit-string RA (hours) + Dec (degrees)
    try:
        coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        return coord.ra.deg, coord.dec.deg
    except Exception:
        pass

    # Attempt 2: Both as decimal degrees
    try:
        coord = SkyCoord(ra_str, dec_str, unit=(u.deg, u.deg))
        return coord.ra.deg, coord.dec.deg
    except Exception:
        pass

    # Attempt 3: ICRS string like "05h00m13.72s -03d20m51.20s"
    try:
        coord = SkyCoord(f"{ra_str} {dec_str}", unit=(u.hourangle, u.deg))
        return coord.ra.deg, coord.dec.deg
    except Exception:
        pass

    # Attempt 4: Full ICRS string in one shot (astropy auto-detect)
    try:
        coord = SkyCoord(f"{ra_str} {dec_str}")
        return coord.ra.deg, coord.dec.deg
    except Exception:
        pass

    raise ValueError(
        f"Could not parse RA='{ra_input}', Dec='{dec_input}'.\n"
        "Accepted formats: '05 00 13.72', '05:00:13.72', '05h00m13.72s',\n"
        "                  '+75.123' (decimal degrees), '75.123d', etc."
    )


def arcsec_to_pixels(arcsec, wcs):
    """
    Convert an angular radius in arcseconds to pixels using a WCS object.
    Assumes the pixel scale is approximately uniform (typical for NuSTAR).

    Parameters
    ----------
    arcsec : float — radius in arcseconds
    wcs    : astropy.wcs.WCS

    Returns
    -------
    radius_pix : float
    """
    import numpy as np
    from astropy.wcs.utils import proj_plane_pixel_scales
    pixel_scales = proj_plane_pixel_scales(wcs)          # degrees / pixel
    pixel_scale_arcsec = np.mean(pixel_scales) * 3600.0  # arcsec / pixel
    return arcsec / pixel_scale_arcsec


def radec_to_pixels(ra_deg, dec_deg, wcs):
    """
    Convert RA/Dec (degrees) to pixel coordinates using a WCS object.

    Returns
    -------
    (x_pix, y_pix) : float pixel coordinates (0-indexed)
    """
    import numpy as np
    x, y = wcs.all_world2pix([[ra_deg, dec_deg]], 0)[0]
    return float(x), float(y)

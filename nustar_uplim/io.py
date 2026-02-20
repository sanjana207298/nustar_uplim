"""
io.py
-----
File discovery and loading for NuSTAR cleaned event files and sky images.

NuSTAR directory structure (post-pipeline):
  {path}/{obsid}/
    event_cl/
      nu{obsid}A01_cl.evt.gz   ← FPMA cleaned events (main science mode)
      nu{obsid}B01_cl.evt.gz   ← FPMB cleaned events
      nu{obsid}A01_sk.img.gz   ← FPMA sky image (histogram of X,Y)
      nu{obsid}B01_sk.img.gz   ← FPMB sky image
      nu{obsid}A01_gti.fits.gz ← FPMA good time intervals
      nu{obsid}B01_gti.fits.gz ← FPMB good time intervals
"""

import os
import glob
import gzip
import shutil
import tempfile

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np


# NuSTAR focal-plane modules
MODULES = ("A", "B")
# Default science mode (01 = standard science)
DEFAULT_MODE = "01"


def find_event_cl_dir(base_path, obsid):
    """
    Locate the event_cl directory for a given observation.

    Parameters
    ----------
    base_path : str — parent directory containing the obsid folder
    obsid     : str — NuSTAR observation ID, e.g. '90601606002'

    Returns
    -------
    event_cl_dir : str — absolute path to event_cl/
    """
    obsid = str(obsid).strip()
    candidates = [
        os.path.join(base_path, obsid, "event_cl"),
        os.path.join(base_path, "event_cl"),           # user already in obsid dir
        os.path.join(base_path, f"{obsid}/event_cl"),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return os.path.abspath(c)
    raise FileNotFoundError(
        f"Cannot find event_cl directory for ObsID {obsid} under '{base_path}'.\n"
        "Expected: {base_path}/{obsid}/event_cl/"
    )


def _find_file(event_cl_dir, pattern):
    """Glob for a single file matching pattern; raise if not found."""
    matches = glob.glob(os.path.join(event_cl_dir, pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching '{pattern}' in {event_cl_dir}")
    return sorted(matches)[0]   # take first / lowest-numbered


def load_events(event_cl_dir, obsid, module, mode=DEFAULT_MODE):
    """
    Load the cleaned event list for one focal-plane module.

    Parameters
    ----------
    event_cl_dir : str
    obsid        : str
    module       : 'A' or 'B'
    mode         : str, default '01'

    Returns
    -------
    events       : numpy recarray  — EVENTS binary table
    header       : fits.Header
    exposure     : float — livetime in seconds
    """
    obsid = str(obsid).strip()
    fname = f"nu{obsid}{module}{mode}_cl.evt.gz"
    path = _find_file(event_cl_dir, fname)

    print(f"  Loading events: {os.path.basename(path)}")
    with fits.open(path, memmap=False) as hdul:
        events  = hdul[1].data
        header  = hdul[1].header
        # Primary header sometimes has LIVETIME
        try:
            exposure = hdul[1].header["LIVETIME"]
        except KeyError:
            try:
                exposure = hdul[0].header["LIVETIME"]
            except KeyError:
                exposure = hdul[1].header.get("EXPOSURE", 1.0)

    print(f"    → {len(events)} events, exposure = {exposure:.1f} s")
    return events, header, float(exposure)


def load_sky_image(event_cl_dir, obsid, module, mode=DEFAULT_MODE):
    """
    Load the sky image (X-Y histogram) for one module.

    Returns
    -------
    image_data : 2-D numpy array
    wcs        : astropy.wcs.WCS (2-D)
    header     : fits.Header
    """
    obsid = str(obsid).strip()
    fname = f"nu{obsid}{module}{mode}_sk.img.gz"
    path = _find_file(event_cl_dir, fname)

    print(f"  Loading sky image: {os.path.basename(path)}")
    with fits.open(path, memmap=False) as hdul:
        # Find the image extension (first IMAGE with NAXIS=2)
        for ext in hdul:
            if ext.data is not None and ext.data.ndim == 2:
                image_data = ext.data.astype(float)
                header = ext.header
                wcs = WCS(header, naxis=2)
                print(f"    → image shape: {image_data.shape}")
                return image_data, wcs, header

    raise ValueError(f"No 2-D image found in {path}")


def load_gti(event_cl_dir, obsid, module, mode=DEFAULT_MODE):
    """
    Load Good Time Intervals.

    Returns
    -------
    gti_start, gti_stop : numpy arrays (seconds)
    """
    obsid = str(obsid).strip()
    fname = f"nu{obsid}{module}{mode}_gti.fits.gz"
    path = _find_file(event_cl_dir, fname)

    with fits.open(path, memmap=False) as hdul:
        gti  = hdul[1].data
        return gti["START"], gti["STOP"]


def get_events_wcs(events_header):
    """
    Build a WCS from the binary table WCS keywords that NuSTAR embeds
    in the events FITS header (TCTYP*, TCRPX*, TCRVL*, TCDLT*).

    NuSTAR events use columns X (col 1-indexed ~11) and Y (col ~12).
    This function finds those columns and builds a 2-D WCS.
    """
    # Find X and Y column indices
    x_col = y_col = None
    for i in range(1, events_header.get("TFIELDS", 50) + 1):
        ttype = events_header.get(f"TTYPE{i}", "").strip().upper()
        if ttype == "X":
            x_col = i
        elif ttype == "Y":
            y_col = i

    if x_col is None or y_col is None:
        raise ValueError("Could not find X and Y columns in events header.")

    # Build minimal header for 2-D WCS
    from astropy.io.fits import Header as FitsHeader
    wcs_hdr = FitsHeader()
    wcs_hdr["NAXIS"]  = 2
    wcs_hdr["NAXIS1"] = 1000
    wcs_hdr["NAXIS2"] = 1000

    for ax_i, col in enumerate([x_col, y_col], start=1):
        wcs_hdr[f"CTYPE{ax_i}"] = events_header.get(f"TCTYP{col}", "RA---TAN" if ax_i == 1 else "DEC--TAN")
        wcs_hdr[f"CRPIX{ax_i}"] = events_header.get(f"TCRPX{col}", 500.0)
        wcs_hdr[f"CRVAL{ax_i}"] = events_header.get(f"TCRVL{col}", 0.0)
        wcs_hdr[f"CDELT{ax_i}"] = events_header.get(f"TCDLT{col}", -1.38889e-4)  # ~0.5 arcsec
        wcs_hdr[f"CUNIT{ax_i}"] = events_header.get(f"TCUNI{col}", "deg")

    wcs = WCS(wcs_hdr, naxis=2)
    return wcs, x_col, y_col


def summary_table(event_cl_dir, obsid):
    """Print a summary of available files."""
    print(f"\n{'='*60}")
    print(f"  NuSTAR ObsID: {obsid}")
    print(f"  event_cl directory: {event_cl_dir}")
    print(f"{'='*60}")
    files = sorted(os.listdir(event_cl_dir))
    for f in files:
        size = os.path.getsize(os.path.join(event_cl_dir, f))
        print(f"  {f:<50s}  {size/1024:>8.1f} KB")
    print(f"{'='*60}\n")

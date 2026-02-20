"""
statistics.py
-------------
Poisson statistics for NuSTAR non-detection upper limits.

Methods implemented
-------------------
1. Frequentist (Gehrels 1986) confidence intervals via astropy
2. Bayesian Poisson upper limit (uniform prior) — Kraft, Burrows & Nousek (1991)
3. Background-inclusive upper limit (signal only, marginalised over background)
4. Detection significance (Li & Ma 1983, their Eq. 17 — used in gamma-ray astronomy
   but widely adopted for X-ray non-detections)

References
----------
Gehrels 1986, ApJ 303, 336
Kraft, Burrows & Nousek 1991, ApJ 374, 344 (KBN)
Li & Ma 1983, ApJ 272, 317
"""

import numpy as np
from scipy import stats, optimize
from astropy.stats import poisson_conf_interval


# ---------------------------------------------------------------------------
# 1. Simple frequentist upper limit on *total* counts (ignores background)
# ---------------------------------------------------------------------------

def frequentist_upper_limit(n_obs, confidence=0.9973):
    """
    Frequentist (Gehrels) upper limit on the *observed* count rate.

    Uses the exact Poisson confidence interval:
        UL such that P(N <= n_obs | UL) = 1 - confidence

    Parameters
    ----------
    n_obs      : int — observed counts
    confidence : float — confidence level (default: 0.9973 ≈ 3σ)

    Returns
    -------
    ul : float — upper limit on total counts
    """
    ul = poisson_conf_interval(
        n_obs,
        interval="frequentist-confidence",
        sigma=_confidence_to_sigma(confidence),
    )[1]
    return float(ul)


def _confidence_to_sigma(confidence):
    """Convert a confidence level to equivalent Gaussian σ."""
    return float(stats.norm.ppf((1 + confidence) / 2))


# ---------------------------------------------------------------------------
# 2. Background-inclusive upper limit (signal, given background estimate B)
# ---------------------------------------------------------------------------

def background_inclusive_upper_limit(n_obs, B, confidence=0.9973):
    """
    Upper limit on the *signal* counts given an expected background B.

    Finds s such that:
        sum_{k=0}^{n_obs} Poisson(k | s + B) = 1 - confidence

    Uses the method of Kraft, Burrows & Nousek (1991).
    For s+B >= 0.

    Parameters
    ----------
    n_obs      : int   — observed counts in source region
    B          : float — expected background counts in source region (α × N_bkg)
    confidence : float — default 0.9973 (3σ)

    Returns
    -------
    s_ul : float — upper limit on net signal counts
    """
    if B < 0:
        B = 0.0

    # Lower bound: s must be >= 0
    # We solve: sum_{k=0}^{n_obs} e^{-(s+B)} (s+B)^k / k!  =  1-confidence
    alpha = 1.0 - confidence

    def equation(s):
        if s < 0:
            return 1.0
        mu = s + B
        # CDF of Poisson at n_obs: P(X <= n_obs | mu)
        cdf = stats.poisson.cdf(n_obs, mu)
        return cdf - alpha

    # If background alone already predicts < n_obs, signal upper limit is small
    # We search in [0, large_value]
    try:
        # Try bracketing
        s_max = max(n_obs * 10 + 50, 100)
        if equation(0) < 0:
            # Upper limit is 0 (B >> n_obs)
            return 0.0
        if equation(s_max) > 0:
            s_max *= 10
        s_ul = optimize.brentq(equation, 0, s_max, xtol=1e-3)
    except ValueError:
        # Fallback: use Gehrels minus B
        s_ul = max(0.0, frequentist_upper_limit(n_obs, confidence) - B)

    return max(0.0, float(s_ul))


# ---------------------------------------------------------------------------
# 3. Bayesian upper limit (uniform prior on signal, exact marginalisation)
# ---------------------------------------------------------------------------

def bayesian_upper_limit(n_obs, B, confidence=0.9973):
    """
    Bayesian upper limit on signal counts with uniform prior (KBN 1991).

    Computes posterior p(s | n_obs, B) ∝ Poisson(n_obs | s+B) and integrates
    to find the credible upper limit.

    Parameters
    ----------
    n_obs      : int
    B          : float
    confidence : float

    Returns
    -------
    s_ul : float
    """
    # Posterior mode is at s = max(0, n_obs - B)
    # Numerically integrate the posterior
    s_mode = max(0.0, n_obs - B)
    s_max  = s_mode + 5 * np.sqrt(max(n_obs, 1)) + 100

    # Grid integration
    n_grid = 10000
    s_vals = np.linspace(0, s_max, n_grid)
    # unnormalised posterior
    log_post = np.array([stats.poisson.logpmf(n_obs, max(s + B, 1e-300))
                         for s in s_vals])
    log_post -= log_post.max()    # numerical stability
    posterior = np.exp(log_post)
    # Normalise
    norm = np.trapz(posterior, s_vals)
    if norm == 0:
        return float(s_max)
    posterior /= norm

    # Integrate CDF and find UL
    cdf = np.cumsum(posterior) * (s_vals[1] - s_vals[0])
    cdf /= cdf[-1]
    idx = np.searchsorted(cdf, confidence)
    if idx >= len(s_vals):
        return float(s_max)
    return float(s_vals[idx])


# ---------------------------------------------------------------------------
# 4. Li & Ma significance
# ---------------------------------------------------------------------------

def lima_significance(n_on, n_off, alpha):
    """
    Li & Ma (1983) Eq. 17 significance of a detection.

    Parameters
    ----------
    n_on   : int   — counts in source region
    n_off  : int   — counts in background region
    alpha  : float — ratio of source to background region areas

    Returns
    -------
    S : float — significance in Gaussian sigma (positive = excess, negative = deficit)
    """
    if n_on == 0 and n_off == 0:
        return 0.0
    if n_off == 0:
        return np.sqrt(2 * n_on * np.log(1.0 / alpha))

    # Eq. 17
    total  = n_on + n_off
    term1  = n_on  * np.log((1 + alpha) / alpha * n_on  / total) if n_on  > 0 else 0.0
    term2  = n_off * np.log((1 + alpha)         * n_off / total) if n_off > 0 else 0.0

    s2 = 2 * (term1 + term2)
    if s2 < 0:
        return -np.sqrt(-s2)
    S = np.sqrt(s2)
    # Sign: positive if excess
    if n_on < alpha * n_off:
        S = -S
    return float(S)


# ---------------------------------------------------------------------------
# 5. WebPIMMS guidance
# ---------------------------------------------------------------------------

# NuSTAR energy band definitions used in WebPIMMS
PIMMS_BAND_INFO = {
    "full":      ("3.0", "79.0"),
    "soft":      ("3.0", "10.0"),
    "hard":      ("10.0", "30.0"),
    "ultrahard": ("30.0", "79.0"),
}

PIMMS_URL = "https://cxc.harvard.edu/toolkit/pimms.jsp"


def format_pimms_instructions(count_rate_ul, band, obsid="", module="", sigma_label="3sigma"):
    """
    Return a formatted string with WebPIMMS instructions for converting
    a NuSTAR count-rate upper limit to a flux upper limit.

    Parameters
    ----------
    count_rate_ul : float   — upper limit on count rate [cts/s]
    band          : str     — energy band key
    obsid         : str
    module        : str
    sigma_label   : str

    Returns
    -------
    instructions : str
    """
    e_lo, e_hi = PIMMS_BAND_INFO.get(band, ("3.0", "79.0"))

    lines = [
        "",
        "  ┌─ WebPIMMS Instructions (" + sigma_label + " UL) ─────────────────────────────┐",
        f"  │  URL: {PIMMS_URL}",
        f"  │",
        f"  │  ObsID {obsid}  FPM{module}  |  Band: {band.upper()} ({e_lo}–{e_hi} keV)",
        f"  │",
        f"  │  STEP 1 — Set INPUT to:",
        f"  │    • Input type:    Count Rate",
        f"  │    • Mission:       NUSTAR",
        f"  │    • Detector:      FPM{module}",
        f"  │    • Energy range:  {e_lo} to {e_hi} keV",
        f"  │    • Count Rate:    {count_rate_ul:.4e} cts/s   ← paste this",
        f"  │",
        f"  │  STEP 2 — Set MODEL (choose based on your source):",
        f"  │    • Power Law    → set Photon Index (e.g. 2.0)",
        f"  │    • Black Body   → set kT",
        f"  │    • Bremsstrahlung → set T",
        f"  │    • APEC/MEKAL  → set T, abundance",
        f"  │    • NH (Galactic): use e.g. https://www.swift.ac.uk/analysis/nhtot/",
        f"  │",
        f"  │  STEP 3 — Set OUTPUT to:",
        f"  │    • Output type:  Flux",
        f"  │    • Energy range: {e_lo} to {e_hi} keV  (or unabsorbed full band)",
        f"  │",
        f"  │  STEP 4 — Click 'Submit' → read off the flux UL",
        f"  └──────────────────────────────────────────────────────────────────┘",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 6. Complete upper limit calculation for a single module
# ---------------------------------------------------------------------------

def compute_upper_limits(extraction_result, exposure_s, band="full",
                         confidence_levels=(0.9545, 0.9973),
                         method="background_inclusive"):
    """
    Compute Poisson upper limits on the signal count rate.

    Flux conversion is intentionally omitted — use WebPIMMS with the
    output count rate and your chosen spectral model.

    Parameters
    ----------
    extraction_result : dict from extraction.extract_counts()
    exposure_s        : float
    band              : str
    confidence_levels : iterable of floats — default (0.9545, 0.9973) = 2σ, 3σ
    method            : 'background_inclusive' | 'frequentist' | 'bayesian'

    Returns
    -------
    results : dict
    """
    ns    = extraction_result["src_counts"]
    nb    = extraction_result["bkg_counts"]
    alpha = extraction_result["alpha"]
    B     = alpha * nb   # expected background in source region
    net   = ns - B
    sig   = lima_significance(ns, nb, alpha)

    ul_func = {
        "background_inclusive": background_inclusive_upper_limit,
        "frequentist":          lambda n, b, c: max(0.0, frequentist_upper_limit(n, c) - b),
        "bayesian":             bayesian_upper_limit,
    }[method]

    ul_dict = {}
    for cl in confidence_levels:
        s_ul        = ul_func(ns, B, cl)
        sigma_label = f"{_confidence_to_sigma(cl):.0f}sigma"
        count_rate_ul = s_ul / exposure_s if exposure_s > 0 else 0.0
        ul_dict[sigma_label] = {
            "confidence":     cl,
            "counts_ul":      s_ul,
            "count_rate_ul":  count_rate_ul,   # cts/s — plug into WebPIMMS
        }

    return {
        "src_counts":    ns,
        "bkg_counts":    nb,
        "expected_bkg":  B,
        "net_counts":    net,
        "alpha":         alpha,
        "lima_sig":      sig,
        "exposure_s":    exposure_s,
        "energy_band":   band,
        "upper_limits":  ul_dict,
        "method":        method,
        "pimms_url":     PIMMS_URL,
    }

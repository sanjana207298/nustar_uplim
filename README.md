# nustar_uplim

**Python pipeline for NuSTAR non-detection upper limits**

> When your source is too faint to appear in the standard NuSTAR pipeline products,
> `nustar_uplim` extracts counts from cleaned event files, computes rigorous Poisson
> upper limits, and converts them to physical flux and luminosity limits — with
> step-by-step diagnostic plots and an interactive background selector.

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Usage](#usage)
4. [The Science](#the-science)
5. [Statistical Methods](#statistical-methods)
6. [Energy Conversion Factors](#energy-conversion-factors)
7. [Output Files](#output-files)
8. [NuSTAR File Structure](#nustar-file-structure)
9. [Limitations and Caveats](#limitations-and-caveats)
10. [References](#references)

---

## Installation

```bash
git clone https://github.com/your-org/nustar_uplim.git
cd nustar_uplim
pip install -r requirements.txt
# Or install as a package:
pip install -e .
```

You will need Python ≥ 3.8 and a working `matplotlib` backend for interactive plots
(install `PyQt5` or `PyQt6` if plots don't appear).

---

## Quick Start

```bash
# Minimal invocation
python nustar_analysis.py \
    --path /data/nustar \
    --obsid 90601606002 \
    --ra "05 00 13.72" \
    --dec "-03 20 51.20"

# Interactive guided mode
python nustar_analysis.py --interactive

# Save results, no interactive window
python nustar_analysis.py \
    --path /data/nustar \
    --obsid 90601606002 \
    --ra "05 00 13.72" \
    --dec "-03 20 51.20" \
    --src-radius 30 \
    --bkg-radius 90 \
    --band soft \
    --bkg-mode annulus \
    --out-dir ./results_sn2024xyz \
    --no-show-plots
```

Or use it as a Python library:

```python
from nustar_uplim import NuSTARUpperLimitPipeline

pipe = NuSTARUpperLimitPipeline(
    base_path     = "/data/nustar",
    obsid         = "90601606002",
    ra            = "05 00 13.72",        # any format
    dec           = "-03 20 51.20",
    src_radius_as = 30.0,                 # arcseconds
    bkg_radius_as = 90.0,
    energy_band   = "soft",               # 3–10 keV
    bkg_mode      = "auto",               # interactive if near edge
    out_dir       = "./results",
    show_plots    = True,
)
results = pipe.run()

# Access numerical results
fpma = results["A"]["stats"]
print(f"3σ flux UL (FPMA): {fpma['upper_limits']['3sigma']['flux_ul']:.2e} erg/cm²/s")
```

### Supported RA/Dec formats

Any of the following work for both `--ra` and `--dec`:

| Format | Example RA | Example Dec |
|---|---|---|
| Sexagesimal (space) | `"05 00 13.72"` | `"-03 20 51.20"` |
| Sexagesimal (colon) | `"05:00:13.72"` | `"-03:20:51.20"` |
| Unit string | `"05h00m13.72s"` | `"-03d20m51.2s"` |
| Decimal degrees | `"75.057167"` | `"-3.347556"` |
| Float | `75.057167` | `-3.347556` |

---

## Usage

### Command-line arguments

| Argument | Default | Description |
|---|---|---|
| `--path` | required | Parent directory of the ObsID folder |
| `--obsid` | required | NuSTAR ObsID string |
| `--ra` | required | Source RA (any format) |
| `--dec` | required | Source Dec (any format) |
| `--src-radius` | 20 | Source circle radius (arcseconds). NuSTAR PSF HPD ~58" so 20" captures the brightest core |
| `--bkg-radius` | 100 | Background region radius (arcseconds) |
| `--bkg-mode` | auto | `auto`, `interactive`, or `annulus` |
| `--band` | soft | Named: `full` (3–79 keV), `soft` (3–10 keV), `hard` (10–30 keV), `ultrahard` (30–79 keV). Custom keV range: `8-30`, `15:50` |
| `--modules` | A B | Which FPM modules to process |
| `--stat-method` | background_inclusive | Statistical method (see below) |
| `--confidence` | 0.9545 0.9973 | Confidence levels (2σ and 3σ) |
| `--out-dir` | — | Save plots and report here |
| `--no-show-plots` | — | Suppress interactive display |
| `--interactive` | — | Prompt for all parameters |

### Background modes

- **`annulus`**: Background is extracted from an annular ring centred on the source
  (inner radius = source radius; outer radius = `--bkg-radius`). Best for isolated
  sources away from chip edges.

- **`interactive`**: A matplotlib window opens showing the sky image with the source
  circle. You click anywhere on the image to place a background circle. Scroll to
  resize it. Press Enter or close the window to confirm.
  Use this when the source is near an edge or there are nearby contaminating sources.

- **`auto`** (default): Uses `annulus` if the source is well away from the chip edge;
  automatically falls back to `interactive` if the source is near an edge.

---

## The Science

### Why non-detections require special treatment

The standard NuSTAR pipeline (`nupipeline` + `nuproducts`) is designed for detected
sources. For undetected sources, it produces no spectrum or meaningful count rate.
However, a *non-detection* is scientifically informative — it places an **upper limit**
on the source flux, which constrains physical models.

This tool bypasses the standard pipeline and works directly with the cleaned event
lists (`_cl.evt.gz`), which are the product of `nupipeline` and contain photon arrival
times, positions, and energies that have already been:

- Filtered for good time intervals (GTIs)
- Corrected for detector patterns (hot pixels, bad columns)
- Transformed into sky (RA/Dec) coordinates

### What the pipeline does

1. **Load cleaned events** from `nu{obsid}{A|B}01_cl.evt.gz`
2. **Define source region**: a circle of chosen radius (arcsec) centred on the
   optical/radio position of the target
3. **Define background region**: annulus around source or user-placed circle,
   avoiding chip gaps, other sources, and readout trails
4. **Extract counts**:
   - Count events inside the source circle: $N_s$
   - Count events inside the background region: $N_b$
   - Compute the area scaling factor: $\alpha = A_s / A_b$
   - Background-scaled expected counts in source region: $B = \alpha N_b$
   - Net (background-subtracted) source counts: $S = N_s - B$
5. **Compute upper limit** on the true signal counts using Poisson statistics
6. **Convert** counts → count rate → flux using an energy conversion factor (ECF)

### Energy bands

| Band | PI range | Energy |
|---|---|---|
| `soft` | 35–210 | 3–10 keV |
| `hard` | 210–710 | 10–30 keV |
| `full` | 35–1935 | 3–79 keV |
| `ultrahard` | 710–1935 | 30–79 keV |
| `8-30` (custom) | 160–710 | 8–30 keV |

NuSTAR channels: 1 PI channel = 0.04 keV; offset = 1.6 keV.
Formula: $\text{PI} = (E_{\text{keV}} - 1.6) / 0.04$

---

## Statistical Methods

### Detection significance: Li & Ma (1983)

The significance of a detection (or non-detection) is computed using Li & Ma (1983),
Eq. 17, which is the standard in high-energy astrophysics:

$$S = \sqrt{2} \left[ N_s \ln\!\left(\frac{1+\alpha}{\alpha} \cdot \frac{N_s}{N_s+N_b}\right) + N_b \ln\!\left((1+\alpha) \cdot \frac{N_b}{N_s+N_b}\right) \right]^{1/2}$$

where $\alpha = A_s / A_b$ is the ratio of source to background region areas.

A positive value indicates a source excess; negative indicates a deficit.
For non-detections you will typically see $|S| < 3$.

### Upper limit methods

Three methods are available via `--stat-method`:

#### 1. `background_inclusive` (default, recommended)

This is the most rigorous method for low-count X-ray astronomy. It finds the
signal count upper limit $s_\text{UL}$ such that the probability of observing
$N_s$ or fewer counts, given signal $s_\text{UL}$ and background $B = \alpha N_b$,
equals $1 - \text{CL}$:

$$\sum_{k=0}^{N_s} \frac{e^{-(s_\text{UL}+B)} \, (s_\text{UL}+B)^k}{k!} = 1 - \text{CL}$$

This correctly accounts for the background uncertainty by treating the expected
background as a known nuisance. Solved numerically via `scipy.optimize.brentq`.

**Method**: Kraft, Burrows & Nousek (1991, ApJ 374, 344)

#### 2. `frequentist`

Simple frequentist upper limit on the *total* count rate using the Gehrels (1986)
approximation, then subtracting the expected background:

$$s_\text{UL} = \text{UL}_\text{Gehrels}(N_s) - B$$

Uses `astropy.stats.poisson_conf_interval` with `interval='frequentist-confidence'`.

Fastest method but slightly conservative when $B$ is large relative to $N_s$.

**Method**: Gehrels (1986, ApJ 303, 336)

#### 3. `bayesian`

Bayesian upper limit using a flat (uniform) prior on the signal and the Poisson
likelihood. The posterior $p(s \,|\, N_s, B) \propto \text{Pois}(N_s; s+B)$ is
integrated numerically, and $s_\text{UL}$ is the credible-interval upper bound.

$$\int_0^{s_\text{UL}} p(s \,|\, N_s, B) \, ds = \text{CL}$$

**Method**: Kraft, Burrows & Nousek (1991) Bayesian version

### Confidence levels

The default outputs are at **2σ** ($\text{CL} = 0.9545$) and **3σ** ($\text{CL} = 0.9973$).
Custom levels can be specified with `--confidence`.

Note: In X-ray non-detection analyses, the **3σ upper limit** is the standard
quantity reported in the literature (e.g. "We report a 3σ upper limit of X erg/s").

---

## Flux Conversion with WebPIMMS

`nustar_uplim` deliberately **does not convert count rates to flux** — because the
conversion depends entirely on the spectral model you assume, and that's a scientific
choice. Instead, the pipeline prints step-by-step WebPIMMS instructions with the
exact count rate to paste in.

**URL**: https://cxc.harvard.edu/toolkit/pimms.jsp

### How to use

1. Run `nustar_uplim` — it prints something like:
   ```
   WebPIMMS Instructions (3sigma UL)
     Input type:  Count Rate
     Mission:     NUSTAR
     Detector:    FPMA
     Energy:      3.0 to 10.0 keV
     Count Rate:  1.7600e-04 cts/s   ← paste this
   ```
2. Go to WebPIMMS and set:
   - **Input**: Count Rate → paste the value above
   - **Mission**: NUSTAR, select the correct FPM (A or B)
   - **Energy range**: match the band you used
   - **Model**: choose your spectral model:
     - Power Law → set photon index Γ (e.g. 2.0)
     - Black Body → set kT
     - APEC/MEKAL → set temperature and abundance
   - **Galactic NH**: use e.g. [the NH calculator](https://www.swift.ac.uk/analysis/nhtot/) for your sky position
   - **Output**: Flux in your chosen band
3. Click **Submit** — read off the flux upper limit

This way you are in full control of the spectral assumptions, and WebPIMMS uses the
official, up-to-date NuSTAR calibration files maintained by the CXC.

---

## Output Files

When `--out-dir` is specified, the following files are written:

| File | Description |
|---|---|
| `step1_sky_image_A.png` | Sky image (log scale) with source region (FPMA) |
| `step1_sky_image_B.png` | Same for FPMB |
| `step2_regions_A.png` | Source + background regions overlaid on image |
| `step3_event_scatter_A.png` | Individual events in source and background regions |
| `step4_upper_limits_<band>.png` | Bar chart of upper limits at each confidence level |
| `step5_radial_profile_A.png` | Count surface density vs radius from source |
| `nustar_uplim_report.png` | Combined single-page report figure |
| `nustar_uplim_results.json` | Machine-readable results (JSON) |
| `nustar_uplim_results.txt` | Human-readable plain-text summary |

### JSON output structure

```json
{
  "obsid": "90601606002",
  "src_ra_deg": 75.057167,
  "src_dec_deg": -3.347556,
  "energy_band": "soft",
  "modules": {
    "FPMA": {
      "exposure_s": 45231,
      "src_counts": 3,
      "bkg_counts": 48,
      "alpha": 0.1111,
      "expected_bkg": 5.33,
      "net_counts": -2.33,
      "lima_sigma": -0.91,
      "upper_limits": {
        "2sigma": {
          "confidence": 0.9545,
          "counts_ul": 4.82,
          "count_rate_ul": 1.07e-4,
          "flux_ul_erg_cm2_s": 3.74e-16
        },
        "3sigma": {
          "confidence": 0.9973,
          "counts_ul": 7.95,
          "count_rate_ul": 1.76e-4,
          "flux_ul_erg_cm2_s": 6.16e-16
        }
      }
    }
  }
}
```

---

## NuSTAR File Structure

After downloading from the HEASARC archive, an observation has this structure:

```
{obsid}/
├── auxil/           ← auxiliary files (orbit, etc.)
├── event_cl/        ← CLEANED EVENTS (input to this pipeline)
│   ├── nu{obsid}A01_cl.evt.gz    ← FPMA mode-01 cleaned events ✓
│   ├── nu{obsid}B01_cl.evt.gz    ← FPMB mode-01 cleaned events ✓
│   ├── nu{obsid}A01_sk.img.gz    ← FPMA sky image (used for WCS) ✓
│   ├── nu{obsid}B01_sk.img.gz    ← FPMB sky image ✓
│   ├── nu{obsid}A01_gti.fits.gz  ← Good time intervals
│   └── ...
├── event_uf/        ← unfiltered events (not needed here)
└── hk/              ← housekeeping files
```

This pipeline uses:
- `*A01_cl.evt.gz` / `*B01_cl.evt.gz` — cleaned events with X,Y,PI,TIME columns
- `*A01_sk.img.gz` / `*B01_sk.img.gz` — sky images for WCS coordinate mapping

---

## Limitations and Caveats

1. **Vignetting** *(handled if exposure map present)*: NuSTAR's effective area
   decreases with off-axis angle — roughly a 10% reduction at 2', 30–40% at 6'.
   The pipeline will automatically apply a vignetting correction if the exposure
   map file (`nu{obsid}{mod}01_ex.img.gz`) exists in the `event_cl` directory.
   This file is **not** produced on download — you need to run `nuexpomap` from
   HEASoft first:
   ```bash
   nuexpomap infile=nu90601606002A01_cl.evt.gz \
             attfile=nu90601606002A.attorb.gz \
             instrprobmapfile=CALDB ...
   ```
   If the file is absent the pipeline prints a warning and proceeds without
   the correction — acceptable for sources within ~2–3' of the optical axis.

2. **PSF effects**: The NuSTAR PSF has a half-power diameter (HPD) of ~58" and
   a FWHM of ~18" at 10 keV, with significant power-law wings. A 20" source
   region (the default) captures ~50% of the HPD; 30" captures ~60%. This means
   the true source counts are higher than measured — divide the inferred count
   rate by the enclosed energy fraction (EEF) for your chosen radius. EEF tables
   for NuSTAR are in the CALDB (`nustar/fpm/bcf/psf/`).

3. **Background spatial uniformity**: The simple area-scaled background assumes
   the background is spatially flat within the region. The NuSTAR background has
   gentle gradients, is higher near chip edges and readout nodes, and contains
   stray light from the optics bench. Use the interactive background mode to
   place the background region carefully — same chip, similar off-axis angle,
   avoiding other sources and chip edges.

4. **Mode 01 only**: This pipeline reads mode-01 (standard science mode) events
   by default. Modes 02–06 (low-efficiency readout, calibration source, etc.) are
   present in `event_cl` but not merged. For very faint sources you can manually
   merge the cleaned event files from all modes before running.

5. **No ARF/RMF convolution**: This tool computes count-rate upper limits, not
   spectra. For a full spectral analysis or to derive ECFs precisely, use
   Sherpa or pyXSPEC with the observation-specific ARF and RMF files.

---

## References

Gehrels, N. 1986, *ApJ*, 303, 336
— Confidence limits for small numbers of events in astrophysical data
[Frequentist Poisson intervals]

Kraft, R. P., Burrows, D. N., & Nousek, J. A. 1991, *ApJ*, 374, 344
— Determination of confidence intervals for the characteristic count rate
of a Poisson process
[Background-inclusive and Bayesian upper limits]

Li, T.-P. & Ma, Y.-Q. 1983, *ApJ*, 272, 317
— Analysis methods for results in gamma-ray astronomy
[Li-Ma detection significance, Eq. 17]

Wik, D. R. et al. 2014, *ApJ*, 792, 48
— NuSTAR Observations of the Bullet Cluster
[NuSTAR energy conversion factors and calibration]

Harrison, F. A. et al. 2013, *ApJ*, 770, 103
— The Nuclear Spectroscopic Telescope Array (NuSTAR) High-energy X-ray Mission
[NuSTAR instrument description, PSF, effective area]

---

## License

MIT License. See `LICENSE` for details.

## Contributing

Issues and pull requests welcome at https://github.com/your-org/nustar_uplim.

For questions about the science or statistics, please open a GitHub Discussion.

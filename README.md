# apass2johnson

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19153359.svg)](https://doi.org/10.5281/zenodo.19153359)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![JAAVSO 2026](https://img.shields.io/badge/paper-JAAVSO%202026-success)](https://doi.org/10.5281/zenodo.19153359)

**Empirical CCD photometric transformations from APASS DR9 (Sloan _g, r, i_) to
the Johnson–Kron–Cousins _B, V, R<sub>C</sub>, I<sub>C</sub>_ system.**

A pure-Python pipeline that downloads APASS DR9 photometry from VizieR, fetches
per-star reddening _E(B−V)_ from the NASA/IPAC IRSA SFD service, applies the
new empirical transformations published in **Montigiani (2026), JAAVSO**, and
returns _V, B, R<sub>C</sub>, I<sub>C</sub>_ magnitudes with full
quadrature-propagated uncertainties — ready to use as comparison stars in
MaxIm DL, AstroImageJ, APT, Muniwin, LesvePhotometry, etc.

---

## Table of contents

- [Why this exists](#why-this-exists)
- [Validation results](#validation-results)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input format](#input-format)
- [Output format](#output-format)
- [The transformation formulae](#the-transformation-formulae)
- [Uncertainty model](#uncertainty-model)
- [Limitations](#limitations)
- [Citing this work](#citing-this-work)
- [Changelog](#changelog)
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## Why this exists

The most widely used _gri → BVR<sub>C</sub>I<sub>C</sub>_ transformations in
the astronomical community are those by [Lupton (2005)](https://www.sdss.org/dr12/algorithms/sdssubvritransform/)
and [Jordi et al. (2006)](https://ui.adsabs.harvard.edu/abs/2006A%26A...460..339J),
both calibrated on samples of about 100–300 stars and exhibiting statistically
significant colour-dependent residuals.

This work re-derives those transformations using:

- a calibration sample of **N = 41,311 APASS DR9 stars** (100–400× the
  previous samples), selected from Landolt-field regions at galactic
  latitudes _|b|_ > 15°;
- **Orthogonal Distance Regression** (ODR), which accounts for measurement
  errors in both predictor and response variables (essential here because
  APASS colours and Landolt magnitudes have comparable uncertainties);
- **AIC-based model selection** between candidate polynomial orders, to
  protect against overfitting;
- independent **validation on N = 287 standards** from four catalogues:
  Landolt (1992, 2013), Clem & Landolt (2013), Stetson/Pancino (2022);
- **bootstrap uncertainty estimates** (N = 1000) on all coefficients.

## Validation results

Improvements over Lupton (2005), measured as the percentage reduction of the
RMS scatter against the 287-star validation sample (all values in mag):

| Band              | RMS Lupton | RMS Montigiani 2026 | 95% CI           | Improvement |
|-------------------|:----------:|:-------------------:|:----------------:|:-----------:|
| **V**             | 0.032      | **0.030**           | 0.026 – 0.033    | **+7%**     |
| **B**             | 0.083      | **0.053**           | 0.045 – 0.061    | **+36%**    |
| **R<sub>C</sub>** | 0.058      | **0.043**           | 0.037 – 0.048    | **+26%**    |
| **I<sub>C</sub>** | 0.105      | **0.073**           | 0.059 – 0.086    | **+30%** *  |

\* _I<sub>C</sub>_ is **provisional**: scatter is still ~0.07 mag and the
validation sample contains only N = 20 M-type stars.

Residual chromatic trends against _(g−r)<sub>0</sub>_ are statistically not
significant for V (p = 0.483), R<sub>C</sub> (p = 0.273), and I<sub>C</sub>
(p = 0.767). For B a residual chromatic slope persists, but at
+0.025 mag/mag is five times smaller than Lupton's (+0.138 mag/mag) and has
no practical impact.

---

## Installation

Requirements: **Python ≥ 3.8** and the following packages:

```bash
pip install numpy pandas astropy astroquery
```

Then clone or download the repository:

```bash
git clone https://github.com/Goldrake71/apass2johnson.git
cd apass2johnson
```

That's it — `apass2johnson.py` is a self-contained script, no further setup
needed.

## Quick start

```bash
python apass2johnson.py
```

You will be prompted for the input CSV file. A minimal `targets.csv` needs
only a name (optional) and equatorial coordinates:

```csv
Name,RA,DEC
SS_Cyg,21 42 42.800,+43 35 09.90
RR_Lyr,19 25 27.913,+42 47 03.69
```

The pipeline will then:

1. **download** APASS DR9 _g, r, i_ photometry for each star from VizieR
   (CDS Strasbourg);
2. **fetch** per-star _E(B−V)_ from the NASA/IPAC IRSA SFD service
   ([Schlafly & Finkbeiner 2011](https://ui.adsabs.harvard.edu/abs/2011ApJ...737..103S)),
   with an analytic fallback if the network is unreachable;
3. **compute** _V, B, R<sub>C</sub>, I<sub>C</sub>_ with full error
   propagation;
4. **save** the results to `johnson_results.csv` (UTF-8 BOM, ready for
   Excel).

A typical output line for a single star looks like:

```
  ► SS_Cyg  (RA=21:42:42.80  Dec=+43:35:09.9)
    APASS DR9: querying VizieR... OK  g=12.020  r=11.812  i=11.715
    E(B-V): attempt 1/3 astroquery IrsaDust... OK  E(B-V) = 0.0974
    ──────────────────────────────────────────────────────────────
    RESULTS  —  SS_Cyg
    g=12.020±0.030  r=11.812±0.029  i=11.715±0.034  E(B-V)=0.0974
    (g-r)=0.208
    ──────────────────────────────────────────────────────────────
    Band    Magnitude   ± σ_tot   Source
    V        11.8861   0.0364    Montigiani (2026)  N=287
    B        12.3122   0.0633    Montigiani (2026)  N=287
    Rc       11.6738   0.0584    Montigiani (2026)  N=287
    Ic       11.4196   0.0922    Montigiani (2026)  N=287  [PROVISIONAL]
    ──────────────────────────────────────────────────────────────
```

## Input format

The pipeline reads a CSV file with comma or semicolon as delimiter
(auto-detected). Lines starting with `#` are treated as comments.

**Minimal format:**

```csv
Name,RA,DEC
SS_Cyg,21 42 42.800,+43 35 09.90
```

**Complete format** (skip the APASS download if you already have photometry):

```csv
Name,RA,DEC,g,eg,r,er,i,ei,EBV
SS_Cyg,21 42 42.800,+43 35 09.90,12.020,0.030,11.812,0.029,11.715,0.034,0.097
```

Coordinates accept space or colon separators (`21 42 42.800` ≡ `21:42:42.800`).
The `Name` column is optional; if missing, stars are labelled `Star_2`, `Star_3`...
based on their line number.

## Output format

The output file `johnson_results.csv` contains the following columns:

| Column        | Description                                                  |
|---------------|--------------------------------------------------------------|
| Name, RA, DEC | identification                                               |
| V, sV, B, sB, Rc, sRc, Ic, sIc | JKC magnitudes and σ_tot           |
| g-r, r-i      | dereddened colour indices used by the formulae               |
| g, eg, r, er, i, ei | APASS DR9 input photometry                             |
| E_BV          | _E(B−V)_ adopted                                             |
| ebv_method    | `SFD/astroquery`, `SFD/urllib`, `analitica|b|=...`, `CSV input`, `manual` |
| in_range      | `YES` / `NO` — whether _(g−r)_ falls in the validity range   |
| ebv_warn      | warning for high extinction                                  |
| warning       | combined warning string (extrapolation, high extinction, ...) |
| source_VB, source_Rc, source_Ic | provenance of each formula                |

Compatible with MaxIm DL, AstroImageJ, APT, Muniwin, LesvePhotometry.

## The transformation formulae

Two equivalent versions are provided, depending on whether the input
APASS magnitudes have been corrected for interstellar reddening or not.

### Version A — dereddened input magnitudes

For _g<sub>0</sub>, r<sub>0</sub>, i<sub>0</sub>_ already corrected with
the SFD maps and Cardelli et al. (1989) extinction coefficients:

```
V₀  = g₀ − 0.01137 − 0.55868·(g−r)₀ + 0.00563·(g−r)₀²
B₀  = g₀ + 0.16841 + 0.44387·(g−r)₀ + 0.02810·(g−r)₀²
Rc,₀ = r₀ − 0.1462  − 0.0763·(g−r)₀  − 0.0426·(g−r)₀²
Ic,₀ = i₀ − 0.3235  + 0.4591·(r−i)₀  − 0.2421·(r−i)₀²  − 0.3115·(g−r)₀
```

Recommended for fields with _E(B−V)_ > 0.20 mag.

### Version B — observed magnitudes (CCD operational use)

Same coefficients, applied directly to APASS DR9 observed _g, r, i_ — no
prior reddening correction required. The second-order bias introduced by
this choice is bounded analytically (Sec. 3.3 of the paper):

| _E(B−V)_ | \|ΔV<sub>2</sub>\| | \|ΔB<sub>2</sub>\| | \|ΔR<sub>C,2</sub>\| | \|ΔI<sub>C,2</sub>\| |
|:---------:|:------------------:|:------------------:|:--------------------:|:--------------------:|
| 0.10 mag  | < 0.0001 mag       | < 0.0003 mag       | < 0.0004 mag         | < 0.001 mag          |
| 0.20 mag  | < 0.0003 mag       | < 0.001 mag        | < 0.002 mag          | < 0.004 mag          |
| 0.50 mag  | < 0.002 mag        | < 0.007 mag        | < 0.010 mag          | < 0.025 mag          |

For typical high-_b_ fields the bias is negligible; for _E(B−V)_ > 0.20 mag,
**Version A with per-star SFD corrections** is strongly recommended.

**Validity:** −0.28 ≤ (g−r) ≤ +1.10 (spectral types OB to late K). Do not
apply outside this range, or to Be, Am, Ap stars, or eclipsing binaries.

## Uncertainty model

The total per-band uncertainty σ<sub>tot</sub> is computed in quadrature
as:

```
σ_tot² = σ_phot² + σ_calib² + σ_ext²
```

where:

- **σ_phot** — APASS photometric uncertainties propagated through the
  analytic partial derivatives of the transformation;
- **σ_calib** — the validation RMS of the formula on the 287-star sample
  (0.030, 0.053, 0.043, 0.073 mag for V, B, R<sub>C</sub>, I<sub>C</sub>);
- **σ_ext** = \|R<sub>band</sub><sup>eff</sup> − R<sub>band</sub>\| ·
  _E(B−V)_, accounting for the residual extinction-related discrepancy
  between Version B output and the nominal JKC system.

The effective extinction coefficients are:

| Band              | _R<sup>eff</sup>_ | _R_ nominal     |
|-------------------|:-----------------:|:---------------:|
| V                 | +3.115            | 3.100           |
| B                 | +4.057            | 4.061           |
| R<sub>C</sub>     | +2.628            | 2.320           |
| I<sub>C</sub>     | +2.061            | 1.490           |

For V and B these match the nominal JKC values to within a few thousandths.
For R<sub>C</sub> and I<sub>C</sub> they differ substantially because those
formulae mix two SDSS bands with different extinction properties — see
Sec. 3.3 of the paper for the full algebraic derivation.

## Limitations

1. **Version B reddening bias.** Use Version A with per-star SFD for fields
   with _E(B−V)_ > 0.20 mag.
2. **Spectral coverage.** The validation sample contains N = 20 M-type
   stars (7%); reliability at the reddest colours is reduced. Do not
   extrapolate beyond −0.28 ≤ (g−r) ≤ +1.10.
3. **I<sub>C</sub> formula is provisional.** Scatter (~0.073 mag) remains
   larger than V/B/R<sub>C</sub>. For I<sub>C</sub> photometry requiring
   σ < 0.05 mag, local calibration with Landolt standards is recommended.
4. **Validation sample size.** N = 287 is limited by availability of
   well-calibrated JKC standards with APASS DR9 counterparts in
   uncrowded fields.
5. **APASS DR9 systematics.** Spatially-correlated zero-point errors at
   the ~0.02 mag level (Henden et al. 2016) constitute an accuracy floor
   that the transformations cannot remove.

## Citing this work

If you use this software or the transformations in published work, please
cite:

> Montigiani, N. 2026, "Empirical Transformation Formulae from APASS DR9
> to the Johnson–Kron–Cousins BVR<sub>C</sub>I<sub>C</sub> Photometric
> System", JAAVSO, in press.
> DOI: [10.5281/zenodo.19153359](https://doi.org/10.5281/zenodo.19153359)

BibTeX:

```bibtex
@article{Montigiani2026,
  author       = {Montigiani, Nico},
  title        = {Empirical Transformation Formulae from APASS DR9 to
                  the Johnson--Kron--Cousins $BVR_{C}I_{C}$
                  Photometric System},
  journal      = {Journal of the AAVSO},
  year         = {2026},
  note         = {in press},
  doi          = {10.5281/zenodo.19153359},
  url          = {https://github.com/Goldrake71/apass2johnson}
}
```

## Changelog

### v5.0.1 — Bug fix & paper alignment

- **Bug fix (I<sub>C</sub> uncertainty):** the residual extinction term in
  the I<sub>C</sub> error propagation was missing the contribution of the
  (g−r) covariate. The corrected expression is
  `|R_Ic − R_i − A1IC·(R_r−R_i) − A3IC·(R_g−R_r)|`, in accordance with
  Sec. 3.3 (eq. 13) of the paper. The σ<sub>Ic</sub> values produced
  by previous releases were overestimated by ~50%.
- All validation statistics, percentages, and sample sizes in code
  docstrings, banner output, and printed results are now aligned with
  the as-accepted JAAVSO manuscript:
  RMS = 0.030/0.053/0.043/0.073 mag for V/B/R<sub>C</sub>/I<sub>C</sub>,
  improvements +7/+36/+26/+30% over Lupton (2005), validated on N = 287
  independent standards (calibration on N = 41,311 APASS DR9 stars).
- Removed obsolete "v4.0 confirmed" / "v5.0 updated" version labels for
  individual bands.

### v5.0 — Initial public release

- ODR-based polynomial transformations from APASS DR9 _gri_ to JKC
  _BVR<sub>C</sub>I<sub>C</sub>_.
- AIC model selection, bootstrap uncertainty estimates.
- Automatic APASS retrieval via VizieR.
- Automatic per-star SFD reddening retrieval via NASA/IPAC IRSA.
- Compatible output for MaxIm DL, AstroImageJ, APT, Muniwin,
  LesvePhotometry.

## Acknowledgements

This software relies on data and infrastructure from the following
projects, gratefully acknowledged:

- **APASS** — AAVSO Photometric All-Sky Survey
  ([Henden et al. 2016](https://ui.adsabs.harvard.edu/abs/2016yCat.2336....0H))
- **VizieR** — Centre de Données astronomiques de Strasbourg
  ([Ochsenbein et al. 2000](https://ui.adsabs.harvard.edu/abs/2000A%26AS..143...23O))
- **NASA/IPAC IRSA** — Dust Extinction Calculator
  ([Schlafly & Finkbeiner 2011](https://ui.adsabs.harvard.edu/abs/2011ApJ...737..103S))
- The Landolt, Clem & Landolt, and Stetson/Pancino photometric standard
  catalogues
- **astropy** ([Astropy Collaboration 2018](https://www.astropy.org/)) and
  **astroquery** ([Ginsburg et al. 2019](https://astroquery.readthedocs.io/))

### AI assistance disclosure

In accordance with the disclosure included in the JAAVSO paper:
portions of the statistical analysis pipeline, the present software, and
parts of the manuscript were developed with the assistance of **Claude
(Anthropic)**. AI assistance was used for: (i) implementation and
testing of the ODR fitting and bootstrap resampling routines;
(ii) error-propagation algebra and the second-order reddening-bias
derivation; (iii) generation of diagnostic figures and validation
scripts; (iv) manuscript drafting and language editing. All scientific
results, numerical values, physical interpretations, and conclusions
were verified independently by the author, who takes full responsibility
for the accuracy and integrity of this work.

## License

Released under the [MIT License](LICENSE). You are free to use, copy,
modify, and distribute this software, with attribution.

---

**Contact:** Nico Montigiani — Osservatorio Astronomico Margherita Hack,
Firenze (MPC code A57). Bug reports, suggestions, and contributions are
welcome via the
[GitHub issue tracker](https://github.com/Goldrake71/apass2johnson/issues).

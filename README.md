# apass2johnson

**Empirical photometric transformation from APASS DR9 *gri* to Johnson-Kron-Cousins *BVR*<sub>C</sub>*I*<sub>C</sub>**

[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## Overview

`apass2johnson` converts [APASS DR9](https://www.aavso.org/apass) Sloan *gri* magnitudes
into the Johnson-Kron-Cousins *BVR*<sub>C</sub>*I*<sub>C</sub> photometric system
using the empirical polynomial formulae of **N. Montigiani (2026)**.

For each input target the software:

1. Retrieves APASS DR9 *gri* photometry automatically via VizieR/CDS
2. Fetches per-star *E*(*B*−*V*) from SFD dust maps via NASA/IPAC IRSA
3. Applies the transformation formulae with full error propagation
4. Outputs results to screen and to `johnson_results.csv`

---

## Transformation formulae (v5.0)

**Version B** — observed APASS magnitudes as input (recommended for CCD use):

```
V   = g  − 0.01137 − 0.55868·(g−r) + 0.00563·(g−r)²
B   = g  + 0.16841 + 0.44387·(g−r) + 0.02810·(g−r)²
Rc  = r  − 0.1462  − 0.0763·(g−r)  − 0.0426·(g−r)²
Ic  = i  − 0.3235  + 0.4591·(r−i)  − 0.2421·(r−i)²  − 0.3115·(g−r)
```

Valid for −0.28 ≤ (*g*−*r*) ≤ +1.10 (spectral types OB to late K).

**Version A** — reddening-corrected *g*₀, *r*₀, *i*₀ as input:
same coefficients applied to extinction-corrected magnitudes.
Recommended for *E*(*B*−*V*) > 0.20 mag.

---

## Validation results

| Band | RMS (this work) | 95% CI | RMS Lupton (2005) | Improvement |
|------|----------------|--------|-------------------|-------------|
| *V*  | 0.032 ± 0.002  | [0.029, 0.036] | 0.035 | +9%  |
| *B*  | 0.058 ± 0.004  | [0.051, 0.066] | 0.096 | **+40%** |
| *R*<sub>C</sub> | 0.051 ± 0.004 | [0.043, 0.059] | 0.060 | +16% |
| *I*<sub>C</sub> | 0.083 ± 0.007 | [0.071, 0.097] | 0.087 | +4% ⚠ provisional |

Bootstrap 95% CI, N=1000 resamplings. Validation on N=281 independent stars:
Landolt (1992, 2013), Clem & Landolt (2013), Stetson/Pancino (2022).

> **Note:** The *I*<sub>C</sub> transformation is provisional — larger scatter
> and residual systematics. For σ < 0.05 mag in *I*<sub>C</sub>, use local
> Landolt calibration.

---

## Installation

```bash
git clone https://github.com/montigiani/apass2johnson.git
cd apass2johnson
pip install -r requirements.txt
```

**Requirements:** Python ≥ 3.8, numpy, pandas, astropy, astroquery.
Internet access required for automatic APASS/IRSA retrieval.

**Windows note:** The software auto-detects terminal colour support.
It works correctly on both Windows Terminal/VS Code (with colours)
and classic cmd.exe (plain text output, no escape codes).

---

## Usage

### Command-line

```bash
python apass2johnson.py
```

The software will prompt for:
1. Path to your input CSV file
2. Whether to download APASS photometry automatically (if not in CSV)
3. How to obtain E(B−V) (automatic from IRSA or manual entry)
4. Output filename (default: `johnson_results.csv`)

### Input file format

Minimum required columns: `RA` and `DEC` (sexagesimal):

```csv
Name,RA,DEC
SS_Cyg,21 42 42.800,+43 35 09.90
RR_Lyr,19 25 27.913,+42 47 03.69
Mira,02 19 20.793,-02 58 39.50
```

Optional columns (retrieved from APASS DR9 if absent):
`g`, `eg`, `r`, `er`, `i`, `ei`, `EBV`

Coordinates accept spaces or colons (`21:42:42.800` or `21 42 42.800`).
Delimiter: comma or semicolon (auto-detected).
Lines starting with `#` are treated as comments.

See `example_input.csv` for a complete example.

### Output

Results printed to screen and saved to `johnson_results.csv`:

```
RESULTS  —  SS_Cyg
g=11.823±0.021  r=11.501±0.018  i=11.374±0.024  E(B-V)=0.0520
(g-r)=0.322  (r-i)=0.127
────────────────────────────────────────────────────────────────
Band        Magnitude   ± σ_tot   Source
────────────────────────────────────────────────────────────────
V             11.6723    0.0439   N. Montigiani (2026) v5.0  N=281
B             12.0214    0.0800   N. Montigiani (2026) v5.0  N=281
Rc            11.3781    0.0791   N. Montigiani (2026) v4.0  N=281
Ic            11.1526    0.0975   N. Montigiani (2026) v4.0  N=281
────────────────────────────────────────────────────────────────
```

---

## Warnings

- ⚠ **E(B−V) ≥ 0.20 mag**: Version A recommended (large extinction uncertainty)
- 🔴 **E(B−V) ≥ 0.50 mag**: star NOT recommended as photometric reference
- ⚠ **(g−r) outside [−0.28, +1.10]**: result is extrapolated beyond valid range
- ⚠ **I**<sub>C</sub>: provisional — use local Landolt calibration for σ < 0.05 mag
- ⚠ **Do not apply** to Be, Am, Ap, or eclipsing binary stars

---

## Citation

If you use this software, please cite:

**N. Montigiani (2026)**, *Empirical Transformation Formulae from APASS DR9
to the Johnson-Kron-Cousins BVR*<sub>C</sub>*I*<sub>C</sub> *Photometric System*,
PASP/AJ, submitted.

```bibtex
@article{Montigiani2026,
  author  = {N. Montigiani},
  title   = {Empirical Transformation Formulae from {APASS} {DR9}
             to the {Johnson}-{Kron}-{Cousins} {BVR$_C$I$_C$}
             Photometric System},
  journal = {PASP},
  year    = {2026},
  note    = {submitted}
}

@software{montigiani_software_2026,
  author    = {N. Montigiani},
  title     = {apass2johnson v5.0},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX},
  url       = {https://github.com/montigiani/apass2johnson}
}
```

---

## References

- Cardelli, J.A., Clayton, G.C. & Mathis, J.S. (1989), ApJ, 345, 245
- Clem, J.L. & Landolt, A.U. (2013), AJ, 146, 88
- Henden, A.A. et al. (2016), yCat 2336 (APASS DR9)
- Landolt, A.U. (1992), AJ, 104, 340
- Landolt, A.U. (2013), AJ, 146, 131
- Lupton, R. (2005), https://www.sdss.org/dr12/algorithms/sdssubvritransform/
- Pancino, E. et al. (2022), A&A, 664, A109
- Schlafly, E.F. & Finkbeiner, D.P. (2011), ApJ, 737, 103

---

## License

[MIT License](LICENSE) — © 2026 N. Montigiani

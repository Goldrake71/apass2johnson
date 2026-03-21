# Changelog

All notable changes to `apass2johnson` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/).

---

## [5.0] — 2026

### Fixed (v5.0 final)
- **KeyError 'Nome'** in `save_results()`: the summary table was still
  referencing the Italian key `'Nome'` instead of the translated `'Name'`.
  This caused a crash after all computations were complete. Fixed.
- **ANSI colour codes on Windows cmd.exe**: codes like `[96m`, `[92m`, `[0m`
  were printed literally on classic Windows Command Prompt instead of being
  interpreted as colours. Added `_supports_color()` auto-detection: colours
  are enabled only when the terminal supports VT100 (Linux/Mac always;
  Windows Terminal / VS Code; disabled on cmd.exe and redirected output).
- **Italian text residuals**: several strings in the banner, usage notes,
  and inline comments were still in Italian. All translated to English.
- **Incorrect RMS/percentage values** in banner: Ic showed `+7%` (should
  be `+4%`); all four bands now show the definitive validation values from
  the published article (N=281 stars, bootstrap 95% CI).

### Changed (v5.0)
- **V and B coefficients updated** from v4.0:
  - V: c₀=−0.01137, c₁=−0.55868, c₂=+0.00563  (RMS=0.032 mag, +9% vs Lupton)
  - B: c₀=+0.16841, c₁=+0.44387, c₂=+0.02810  (RMS=0.058 mag, +40% vs Lupton)
- Validation extended to N=281 stars from 4 independent catalogues
  (Landolt 1992, 2013; Clem & Landolt 2013; Stetson/Pancino 2022)
- Bootstrap 95% CI reported for all RMS values
- Output file renamed from `risultati_johnson.csv` to `johnson_results.csv`
- All messages, prompts, comments, and output fully translated to English
- Reference updated to N. Montigiani (2026)

### Unchanged (v5.0)
- Rc and Ic coefficients: confirmed from v4.0
- Error propagation framework (photometric + calibration + extinction)
- APASS DR9 / IRSA retrieval pipeline

---

## [4.0] — 2024

### Added
- **Rc band** (v4.0): c₀=−0.1462, c₁=−0.0763, c₂=−0.0426  (RMS=0.051 mag)
- **Ic band — 4-parameter model** (v4.0):
  c₀=−0.3235, c₁=+0.4591, c₂=−0.2421, c₃=−0.3115  (RMS=0.083 mag, provisional)
- Automatic E(B−V) retrieval from SFD via NASA/IPAC IRSA
- Version A (dereddened) and Version B (observed magnitudes) modes
- Regression upgraded from OLS to ODR (Orthogonal Distance Regression)
- Model selection via AIC (Akaike Information Criterion)

---

## [2.0] — 2023

### Added
- Automatic APASS DR9 retrieval via VizieR/CDS
- CSV input/output format compatible with MaxIm DL, APT, Muniwin
- V and B initial coefficients (later updated in v5.0)
- Calibration on 41,311 APASS stars cross-matched with Landolt fields

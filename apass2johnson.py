#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
  apass2johnson.py  —  v5.0
  CCD Photometric Transformation: SDSS (g,r,i) → Johnson-Cousins (V,B,Rc,Ic)
================================================================================

  REFERENCE
    N. Montigiani (2026), "Empirical Transformation Formulae from APASS DR9
    to the Johnson-Kron-Cousins BVRcIc Photometric System", PASP/AJ, submitted.
    GitHub: https://github.com/montigiani/apass2johnson
    DOI:    https://doi.org/10.5281/zenodo.19153253

  TRANSFORMATION FORMULAE — Version B (observed APASS DR9 magnitudes as input)
  ─────────────────────────────────────────────────────────────────────────────
    V  = g − 0.01137 − 0.55868·(g−r) + 0.00563·(g−r)²  [v5.0, RMS=0.032 mag]
    B  = g + 0.16841 + 0.44387·(g−r) + 0.02810·(g−r)²  [v5.0, RMS=0.058 mag]
    Rc = r − 0.1462  − 0.0763·(g−r)  − 0.0426·(g−r)²   [v4.0, RMS=0.051 mag]
    Ic = i − 0.3235  + 0.4591·(r−i)  − 0.2421·(r−i)²
             − 0.3115·(g−r)                              [v4.0, RMS=0.083 mag, PROVISIONAL]

    Valid for: −0.28 ≤ (g−r) ≤ +1.10  (spectral types OB to late K)

  Version A (reddening-corrected input): same coefficients, apply to
    g₀ = g − Rg·E(B−V),  r₀ = r − Rr·E(B−V),  i₀ = i − Ri·E(B−V)
    Recommended for E(B−V) > 0.20 mag.

  TOTAL UNCERTAINTY (quadrature):
    σ_tot² = σ_APASS²(propagated) + σ_calibration² + σ_E(B-V)²

  VALIDATION (N=281 independent standard stars, 4 catalogues):
    V:  RMS=0.032±0.002  +9%  vs Lupton (2005)   slope p=0.765 n.s.
    B:  RMS=0.058±0.004  +40% vs Lupton (2005)   slope p=0.804 n.s. vs E(B-V)
    Rc: RMS=0.051±0.004  +16% vs Lupton (2005)   slope p=0.732 n.s.
    Ic: RMS=0.083±0.007  +4%  vs Lupton (2005)   PROVISIONAL

  HOW TO USE
  ──────────
  1. Prepare a CSV file with your targets (see INPUT FORMAT below).
  2. Run:  python apass2johnson.py
  3. Enter the path to your CSV file when prompted.
  4. The software will:
       a. Retrieve APASS DR9 g,r,i photometry via VizieR (if not in CSV)
       b. Fetch per-star E(B-V) from NASA/IPAC IRSA SFD maps
       c. Compute V, B, Rc, Ic with full error propagation
       d. Print results to screen and save to johnson_results.csv

  INPUT FORMAT (CSV)
  ──────────────────
  Minimum required columns: RA and DEC (sexagesimal)

    Name,RA,DEC
    SS_Cyg,21 42 42.800,+43 35 09.90
    RR_Lyr,19 25 27.913,+42 47 03.69

  Optional columns (if absent, retrieved from APASS DR9 automatically):
    g, eg, r, er, i, ei, EBV

  Full example with photometry already known:
    Name,RA,DEC,g,eg,r,er,i,ei
    SS_Cyg,21 42 42.800,+43 35 09.90,11.823,0.021,11.501,0.018,11.374,0.024

  Coordinates accept spaces or colons (21:42:42.800 or 21 42 42.800).
  Delimiter: comma or semicolon (auto-detected). Lines starting with # are comments.

  OUTPUT FILE (johnson_results.csv)
  ──────────────────────────────────
  Columns: Name, RA, DEC, V, sV, B, sB, Rc, sRc, Ic, sIc,
           g-r, r-i, g, eg, r, er, i, ei, E_BV, ebv_method,
           in_range, ebv_warn, warning, source_VB, source_Rc, source_Ic

  Compatible with: MaxIm DL, AstroImageJ, APT, Muniwin, LesvePhotometry.

  INSTALL DEPENDENCIES
  ─────────────────────
    pip install astropy astroquery numpy pandas

  WARNINGS
  ────────
  ⚠ E(B-V) ≥ 0.20 mag : Version A recommended; large extinction error
  🔴 E(B-V) ≥ 0.50 mag : star NOT recommended as photometric reference
  ⚠ (g-r) outside [-0.28, +1.10]: result is extrapolated
  ⚠ Ic band: provisional — use local Landolt calibration for σ < 0.05 mag
  ⚠ Do not apply to Be, Am, Ap, or eclipsing binary stars

================================================================================
"""

import sys
import os
import re
import math
import time
import warnings
import urllib.request
import xml.etree.ElementTree as ET
from typing import Optional

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Dipendenze obbligatorie ───────────────────────────────────────────────────
try:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
except ImportError:
    sys.exit("ERRORE: astropy non installata.\nEseguire: pip install astropy")

try:
    from astroquery.vizier import Vizier
    ASTROQUERY_OK = True
except ImportError:
    ASTROQUERY_OK = False
    print("AVVISO: astroquery non installata. Download APASS disabilitato.")
    print("        Eseguire: pip install astroquery\n")


# ══════════════════════════════════════════════════════════════════════════════
#  COSTANTI
# ══════════════════════════════════════════════════════════════════════════════

# Extinction coefficients Cardelli et al. (1989), R_V=3.1
R_EXT = {"g": 3.640, "r": 2.700, "i": 2.060, "V": 3.100, "B": 4.061,
         "Rc": 2.320, "Ic": 1.490}  # Cousins R e I (Cardelli+1989)

# Transformation coefficients (N. Montigiani 2026, v5.0)
# Calibrazione ODR su N=281-281 stars (Landolt 1992/2013 + Clem 2013 + Stetson 2022)
# ─────────────────────────────────────────────────────────────────────────────
# V = g + A0V + A1V*(g-r) + A2V*(g-r)^2   [v5.0 — aggiornata con Stetson 2022]
# RMS=0.032 mag su N=281, +9% vs Lupton (2005), slope p=0.89 n.s.
A0V, A1V, A2V   = -0.01137, -0.55868, +0.00563
# B = g + A0B + A1B*(g-r) + A2B*(g-r)^2   [v5.0 — aggiornata con Stetson 2022]
# RMS=0.058 mag su N=281, +40% vs Lupton (2005), slope p=0.69 n.s. su L92
A0B, A1B, A2B   = +0.16841, +0.44387, +0.02810
# Rc = r + A0RC + A1RC*(g-r) + A2RC*(g-r)^2  [v4.0 — confirmed, variation <1%]
# RMS=0.051 mag su N=281, +16% vs Lupton (2005), slope p=0.73 n.s.
A0RC, A1RC, A2RC = -0.1462, -0.0763, -0.0426
# Ic = i + A0IC + A1IC*(r-i) + A2IC*(r-i)^2 + A3IC_GR*(g-r)  [v4.0 — confirmed]
# Modello M3 (AIC-ottimale). Re-fit con Stetson instabile → coefficienti invariati.
# RMS=0.083 mag su N=281, +4% vs Lupton (2005), slope p=0.98 n.s.
A0IC, A1IC, A2IC, A3IC_GR = -0.3235, +0.4591, -0.2421, -0.3115

# RMS di calibrazione su campione di validazione indipendente (N=281-281 stars)
RMS_V  = 0.0275  # v5.0: Landolt 1992/2013 + Clem 2013 + Stetson 2022, N=281
RMS_B  = 0.0582  # v5.0: Landolt 1992/2013 + Clem 2013 + Stetson 2022, N=281
RMS_RC = 0.0501  # v4.0 confirmed: validated on N=281 stars
RMS_IC = 0.0793  # v4.0 confirmed: validated on N=281 stars

# Valid colour range
GR_MIN, GR_MAX = -0.28, +1.10

# VizieR search radius (arcsec)
SEARCH_RADIUS = 15.0

# Colori ANSI
def _supports_color() -> bool:
    """
    Return True if the terminal supports ANSI colour codes.
    Handles Windows cmd.exe, PowerShell, Jupyter, redirected output, etc.
    """
    import sys, os
    # Not a TTY (e.g. output redirected to file or pipe) → no colour
    if not hasattr(sys.stdout, "isatty") or not sys.stdout.isatty():
        return False
    # Windows: enable VT100 via SetConsoleMode, or fall back gracefully
    if os.name == "nt":
        try:
            import ctypes
            kernel32 = ctypes.windll.kernel32
            # Enable ENABLE_VIRTUAL_TERMINAL_PROCESSING (0x0004)
            handle = kernel32.GetStdHandle(-11)          # STD_OUTPUT_HANDLE
            mode   = ctypes.c_ulong()
            if kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
                new_mode = mode.value | 0x0004
                if kernel32.SetConsoleMode(handle, new_mode):
                    return True   # modern Windows Terminal / VS Code terminal
        except Exception:
            pass
        # Classic cmd.exe on older Windows (pre-10): no colour support
        return False
    # Unix/Mac/Linux: always supported
    return True

_COLOR = _supports_color()

BOLD  = "\033[1m"  if _COLOR else ""
CYAN  = "\033[96m" if _COLOR else ""
GREEN = "\033[92m" if _COLOR else ""
YELL  = "\033[93m" if _COLOR else ""
RED   = "\033[91m" if _COLOR else ""
DIM   = "\033[2m"  if _COLOR else ""
RST   = "\033[0m"  if _COLOR else ""


# ══════════════════════════════════════════════════════════════════════════════
#  COORDINATE PARSING
# ══════════════════════════════════════════════════════════════════════════════

def parse_radec(ra_str: str, dec_str: str) -> SkyCoord:
    """
    Convert sexagesimal RA and Dec strings into a SkyCoord object.
    Accepts space or colon as separator; e.g.:
        '12 34 56.789'  or  '12:34:56.789'
    """
    ra_str  = re.sub(r"[:\s]+", " ", ra_str.strip())
    dec_str = re.sub(r"[:\s]+", " ", dec_str.strip())
    return SkyCoord(ra=ra_str, dec=dec_str,
                    unit=(u.hourangle, u.deg), frame="icrs")


def coord_hmsdms(coord: SkyCoord) -> tuple[str, str]:
    """Return (RA_str, DEC_str) as formatted strings hh:mm:ss.ss / ±dd:mm:ss.s"""
    ra  = coord.ra.to_string(unit=u.hourangle, sep=":", precision=2, pad=True)
    dec = coord.dec.to_string(unit=u.deg, sep=":", precision=1,
                               alwayssign=True, pad=True)
    return ra, dec


# ══════════════════════════════════════════════════════════════════════════════
#  E(B-V) RETRIEVAL FROM NASA/IPAC IRSA
# ══════════════════════════════════════════════════════════════════════════════

def _parse_irsa_xml(xml_txt: str) -> Optional[float]:
    """
    Analizza la risposta XML del servizio IRSA DUST.
    ATTENZIONE: il primo <result> contiene intensità in MJy/sr (NON E(B-V)).
    Filtriamo solo i blocchi con unità "mag" e preferiamo Schlafly+2011.
    Cap fisico a 5.0 mag come salvaguardia contro valori assurdi.
    """
    root = ET.fromstring(xml_txt)
    candidates = []
    # Il tag XML IRSA può essere <r> (risposta reale) o <result> (documentazione)
    all_blocks = list(root.iter("result")) or list(root.iter("r"))
    for result in all_blocks:
        desc_el = result.find("desc")
        mv_el   = result.find(".//meanValue")
        if desc_el is None or mv_el is None:
            continue
        mv_text   = (mv_el.text or "").strip()
        desc_text = (desc_el.text or "").lower()
        if "mag" not in mv_text.lower():
            continue          # ignora MJy/sr e altre unità non-mag
        try:
            val = float(mv_text.split()[0])
        except (ValueError, IndexError):
            continue
        if not (0.0 <= val <= 5.0):
            continue          # sanity check fisico
        priority = 0 if "schlafly" in desc_text else 1  # Schlafly+2011 preferito
        candidates.append((priority, val))
    if candidates:
        candidates.sort()
        return candidates[0][1]
    return None


def fetch_ebv_irsa_astroquery(coord: SkyCoord, timeout: int = 25) -> Optional[float]:
    """
    Method 1: use astroquery.irsa_dust.IrsaDust.
    Advantage: shares the same HTTP/proxy session as astroquery/VizieR.
    If VizieR works, this should work too.
    """
    if not ASTROQUERY_OK:
        return None
    try:
        # Support both old and new astroquery import paths
        try:
            from astroquery.ipac.irsa.irsa_dust import IrsaDust
        except ImportError:
            from astroquery.irsa_dust import IrsaDust  # astroquery < 0.4.4

        table = IrsaDust.get_query_table(coord, section="ebv", timeout=timeout)
        if table is None or len(table) == 0:
            return None

        # Possible column names: 'E(B-V) SandF' (Schlafly+2011) or 'E(B-V) SFD'
        for col_name in ("E(B-V) SandF", "E(B-V) SFD", "ext SandF mean",
                         "ext SFD mean"):
            if col_name in table.colnames:
                val = float(table[col_name][0])
                if 0.0 <= val <= 5.0:
                    return val
        # Fallback: prima colonna numerica entro range fisico
        for col in table.colnames:
            try:
                val = float(table[col][0])
                if 0.0 <= val <= 5.0:
                    return val
            except (ValueError, TypeError):
                continue
    except Exception:
        pass
    return None


def fetch_ebv_irsa_urllib(coord: SkyCoord,
                          timeout: int = 25,
                          retries: int = 2) -> Optional[float]:
    """
    Method 2: direct HTTP call to IRSA via urllib.
    Supports system proxy auto-detection.
    Two attempts with progressive timeout.
    """
    ra, dec = coord.ra.deg, coord.dec.deg
    url = (f"https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?"
           f"locstr={ra:.6f}+{dec:+.6f}+Equ+J2000")

    # Auto-detect system proxy (Windows/Mac/Linux)
    try:
        proxies = urllib.request.getproxies()
        if proxies:
            proxy_handler = urllib.request.ProxyHandler(proxies)
            opener = urllib.request.build_opener(proxy_handler)
        else:
            opener = urllib.request.build_opener()
    except Exception:
        opener = urllib.request.build_opener()

    for attempt in range(retries):
        t = timeout + attempt * 10   # 25s → 35s al secondo tentativo
        try:
            resp = opener.open(url, timeout=t)
            xml_txt = resp.read().decode("utf-8")
            result = _parse_irsa_xml(xml_txt)
            if result is not None:
                return result
        except Exception:
            if attempt < retries - 1:
                time.sleep(1)   # short pause between attempts
            continue
    return None


def ebv_from_latitude(coord: SkyCoord) -> float:
    """
    Method 3 (offline fallback): analytic estimate E(B-V) ≈ 0.03/sin|b|.
    Adequate for |b| > 15°; capped at 0.50 mag.
    """
    b_rad = abs(math.radians(coord.galactic.b.deg))
    return round(min(0.03 / max(math.sin(b_rad), 0.15), 0.50), 4)


def get_ebv(coord: SkyCoord, star_name: str = "") -> tuple[float, str]:
    """
    Retrieve E(B-V) for the given coordinates using a 3-level strategy:

    1. astroquery IrsaDust  — same infrastructure as VizieR, most reliable
    2. urllib direct IRSA   — with system proxy detection and auto-retry
    3. Analytic approximation 0.03/sin|b| — always available offline

    Diagnostic messages indicate which method succeeded
    or why the fallback was used.
    """
    print(f"    E(B-V): attempt 1/3 astroquery IrsaDust...", end="", flush=True)
    ebv = fetch_ebv_irsa_astroquery(coord)
    if ebv is not None:
        print(f" {GREEN}OK{RST}  E(B-V) = {ebv:.4f} [Schlafly & Finkbeiner 2011, via astroquery]")
        return ebv, "SFD/astroquery"

    print(f" {YELL}no{RST}")
    print(f"    E(B-V): attempt 2/3 urllib→IRSA (timeout 25-35s)...",
          end="", flush=True)
    ebv = fetch_ebv_irsa_urllib(coord)
    if ebv is not None:
        print(f" {GREEN}OK{RST}  E(B-V) = {ebv:.4f} [Schlafly & Finkbeiner 2011, via urllib]")
        return ebv, "SFD/urllib"

    # Diagnostica del fallimento
    print(f" {YELL}no{RST}")
    print(f"    {YELL}⚠ IRSA not reachable. Possible causes:{RST}")
    print(f"    {YELL}  • Firewall/proxy blocking irsa.ipac.caltech.edu{RST}")
    print(f"    {YELL}  • No internet connection{RST}")
    print(f"    {YELL}  • IRSA server temporarily offline{RST}")
    print(f"    {YELL}  Check: open in browser →{RST}")
    print(f"    {YELL}  https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust"
          f"?locstr={coord.ra.deg:.3f}+{coord.dec.deg:+.3f}+Equ+J2000{RST}")
    print(f"    {YELL}  If browser shows XML → enter E(B-V) manually{RST}")

    ebv = ebv_from_latitude(coord)
    b   = coord.galactic.b.deg
    print(f"    {YELL}↳ Analytic fallback: E(B-V)≈0.03/sin|b| = "
          f"{ebv:.4f}  (b={b:.1f}°){RST}")
    return ebv, f"analitica|b|={abs(b):.1f}°"


# ══════════════════════════════════════════════════════════════════════════════
#  APASS DR9 PHOTOMETRY RETRIEVAL FROM VIZIER
# ══════════════════════════════════════════════════════════════════════════════

# Exhaustive map of possible VizieR column names for II/336/apass9,
# varying by astroquery version. VizieR may return:
#   "g_mag" / "g'mag" / "gmag"   (and equivalently for r, i)
#   "e_g_mag" / "e_g'mag" / "e_gmag"
# All normalised to g, r, i, eg, er, ei.

COL_ALIASES = {
    # g magnitude
    "g":  ["g_mag", "g'mag", "gmag", "gSDSS", "g_SDSS"],
    # r magnitude
    "r":  ["r_mag", "r'mag", "rmag", "rSDSS", "r_SDSS"],
    # i magnitude
    "i":  ["i_mag", "i'mag", "imag", "iSDSS", "i_SDSS"],
    # uncertainties
    "eg": ["e_g_mag", "e_g'mag", "e_gmag", "eg_mag", "e_gSDSS"],
    "er": ["e_r_mag", "e_r'mag", "e_rmag", "er_mag", "e_rSDSS"],
    "ei": ["e_i_mag", "e_i'mag", "e_imag", "ei_mag", "e_iSDSS"],
    # Johnson (optional, for cross-check)
    "V":  ["Vmag", "V_mag"],
    "B":  ["Bmag", "B_mag"],
    "eV": ["e_Vmag", "eVmag"],
    "eB": ["e_Bmag", "eBmag"],
}


def _resolve_col(table_cols: list, target: str) -> Optional[str]:
    """
    Given the list of columns in a VizieR table, return the actual column
    name matching the target (e.g. "g"), or None if not found.
    Checks exact name first, then known aliases.
    """
    # Exact match (case-insensitive)
    low_cols = {c.lower(): c for c in table_cols}
    if target.lower() in low_cols:
        return low_cols[target.lower()]
    # Alias lookup
    for alias in COL_ALIASES.get(target, []):
        if alias.lower() in low_cols:
            return low_cols[alias.lower()]
    return None


def _safe_float(value) -> Optional[float]:
    """Convert to float; return None for NaN/None/empty string."""
    try:
        v = float(value)
        return None if (math.isnan(v) or math.isinf(v)) else v
    except (TypeError, ValueError):
        return None


def fetch_apass(coord: SkyCoord,
                radius_arcsec: float = SEARCH_RADIUS,
                verbose: bool = True) -> Optional[dict]:
    """
    Download APASS DR9 photometry for the nearest star to the given coordinates
    via VizieR (astroquery).

    Progressive search radius:
      1. Try 15" (precise catalogue coordinates)
      2. If empty, expand to 30" (approximate coordinates)
      3. If empty, expand to 60" (chart/visual estimate)
      Always selects the geometrically closest star.
      Warns if the nearest match is >15" from the requested position.

    Column names normalised via COL_ALIASES to handle astroquery version
    differences (g_mag / g'mag / gmag etc.).
    """
    if not ASTROQUERY_OK:
        if verbose:
            print(f"    {RED}astroquery non disponibile{RST}")
        return None

    if verbose:
        print(f"    APASS DR9: querying VizieR...", end="", flush=True)

    # Progressive search radii (arcsec)
    radii   = [15.0, 30.0, 60.0]
    tbl         = None
    used_radius = None

    for radius in radii:
        try:
            viz = Vizier(columns=["**"], row_limit=500, timeout=30)
            results = viz.query_region(
                coord,
                radius=radius * u.arcsec,
                catalog="II/336/apass9",
            )
            if results and len(results) > 0 and len(results[0]) > 0:
                tbl         = results[0]
                used_radius = radius
                break
        except Exception as exc:
            if verbose:
                print(f" {RED}errore rete: {exc}{RST}")
            return None

    if tbl is None:
        if verbose:
            print(f" {YELL}no star found within 60\"{RST}")
        return None

    cols = list(tbl.colnames)

    # ── Find closest row ─────────────────────────────────────────────────────
    ra_col  = _resolve_col(cols, "RAJ2000") or _resolve_col(cols, "ra")
    dec_col = _resolve_col(cols, "DEJ2000") or _resolve_col(cols, "dec")

    best_row  = 0
    best_dist = None
    if ra_col and dec_col:
        try:
            ras  = np.array([float(tbl[ra_col][i])  for i in range(len(tbl))],
                            dtype=float)
            decs = np.array([float(tbl[dec_col][i]) for i in range(len(tbl))],
                            dtype=float)
            # angular distance in arcsec (with cos(dec) correction)
            dists = np.sqrt(
                ((ras - coord.ra.deg) * np.cos(np.radians(coord.dec.deg)))**2
                + (decs - coord.dec.deg)**2
            ) * 3600.0
            best_row  = int(np.argmin(dists))
            best_dist = float(dists[best_row])
        except Exception:
            best_row  = 0
            best_dist = None

    row = tbl[best_row]

    # ── Distance and search radius message ───────────────────────────────────
    dist_msg = ""
    if best_dist is not None and best_dist > 15.0:
        dist_msg = (f" {YELL}⚠ nearest star at {best_dist:.1f}\""
                    f" (approximate coordinates?){RST}")

    # ── Extract magnitudes using alias map ────────────────────────────────────
    phot = {}
    for target in ("g", "r", "i", "eg", "er", "ei", "V", "B", "eV", "eB"):
        col_name = _resolve_col(cols, target)
        phot[target] = _safe_float(row[col_name]) if col_name else None

    # ── Minimum check: g and r required ──────────────────────────────────────
    if phot["g"] is None or phot["r"] is None:
        if verbose:
            print(f" {YELL}found but g/r unreadable{RST}")
            print(f"    {DIM}Available columns: {cols}{RST}")
        return None

    if verbose:
        g_s = f"{phot['g']:.3f}"
        r_s = f"{phot['r']:.3f}"
        i_s = f"{phot['i']:.3f}" if phot["i"] else "N/A"
        print(f" {GREEN}OK{RST}  g={g_s}  r={r_s}  i={i_s}{dist_msg}")

    return phot


# ══════════════════════════════════════════════════════════════════════════════
#  TRANSFORMATION FORMULAE AND ERROR PROPAGATION
# ══════════════════════════════════════════════════════════════════════════════

def compute_johnson(g: float, r: float, i: Optional[float],
                    eg: float, er: float, ei: Optional[float],
                    ebv: float, ebv_method: str) -> dict:
    """
    Calcola V, B, Rc, Ic con le formule N. Montigiani (2026) v5.0 per uso
    diretto nelle immagini CCD (magnitudini SDSS osservate, senza correzione
    per estinzione — Version B).

    FORMULE (Version B - magnitudini osservate APASS DR9):
        V  = g + A0V  + A1V*(g-r)  + A2V*(g-r)^2    [v5.0 — aggiornata]
        B  = g + A0B  + A1B*(g-r)  + A2B*(g-r)^2    [v5.0 — aggiornata]
        Rc = r + A0RC + A1RC*(g-r) + A2RC*(g-r)^2   [v4.0 — confirmed]
        Ic = i + A0IC + A1IC*(r-i) + A2IC*(r-i)^2 + A3IC_GR*(g-r)  [v4.0 — confirmed]

    PARTIAL DERIVATIVES for error propagation:
        V:   dV/dg = 1+A1V = +0.441      dV/dr = -A1V = +0.559
        B:   dB/dg = 1+A1B = +1.444      dB/dr = -A1B = +0.444
        Rc:  dRc/d(g-r) = A1RC + 2*A2RC*(g-r)  (dipende dal colore)
             dRc/dr = 1 - dRc/d(g-r)
             dRc/dg = +dRc/d(g-r)
        Ic:  dIc/di = 1 - A1IC - 2*A2IC*(r-i)   (dipende dal colore)
             dIc/dr = A1IC + 2*A2IC*(r-i) - A3IC_GR
             dIc/dg = +A3IC_GR = -0.312 (fisso)

    Novità v5.0 rispetto a v4.0:
        V: aggiornata con campione esteso a N=281 stars (Landolt 1992/2013 +
           Clem & Landolt 2013 + Stetson in Pancino 2022). RMS=0.032 mag,
           miglioramento +9% vs Lupton (2005). Slope residui p=0.89 n.s.
        B: aggiornata con stesso campione N=281 stars. RMS=0.058 mag,
           miglioramento +40% vs Lupton (2005). Slope p=0.69 n.s. (L92).
        Rc: confirmed unchanged (variation <1% with extended sample).
        Ic: confirmed unchanged (ODR unstable with Stetson subsample).
    """
    gr = g - r

    # Colour range check
    in_range = GR_MIN <= gr <= GR_MAX
    warn = ""
    if not in_range:
        warn = (f"(g-r)={gr:.3f} fuori range [{GR_MIN},{GR_MAX}]: "
                f"result extrapolato")

    # Default APASS uncertainties if missing
    eg = eg if eg and eg > 0 else 0.050
    er = er if er and er > 0 else 0.050
    ei = ei if (ei and ei > 0) else 0.070

    # ── V band ───────────────────────────────────────────────────────────────
    V = g + A0V + A1V * gr + A2V * gr**2
    dV_dg    = 1.0 + A1V          # = +0.441 con nuovi coefficienti
    dV_dr    = -A1V               # = +0.559
    sV_phot  = math.sqrt((dV_dg * eg)**2 + (dV_dr * er)**2)
    delta_RV = abs(R_EXT["V"] - R_EXT["g"] - A1V * (R_EXT["g"] - R_EXT["r"]))
    sV_ebv   = delta_RV * ebv
    sV       = math.sqrt(sV_phot**2 + RMS_V**2 + sV_ebv**2)

    # ── B band ───────────────────────────────────────────────────────────────
    B = g + A0B + A1B * gr + A2B * gr**2
    dB_dg    = 1.0 + A1B          # = +1.444 con nuovi coefficienti
    dB_dr    = -A1B               # = +0.444
    sB_phot  = math.sqrt((dB_dg * eg)**2 + (dB_dr * er)**2)
    delta_RB = abs(R_EXT["B"] - R_EXT["g"] - A1B * (R_EXT["g"] - R_EXT["r"]))
    sB_ebv   = delta_RB * ebv
    sB       = math.sqrt(sB_phot**2 + RMS_B**2 + sB_ebv**2)

    # ── Rc band [v4.0 confirmed — degree-2, N=281 stars] ─
    # Rc = r + A0RC + A1RC*(g-r) + A2RC*(g-r)^2
    # Derivate: dRc/d(g-r) = A1RC + 2*A2RC*(g-r)  varia con il colore
    Rc = r + A0RC + A1RC * gr + A2RC * gr**2
    dRc_dgr  = A1RC + 2.0 * A2RC * gr   # derivata rispetto a (g-r)
    dRc_dr   = 1.0 - dRc_dgr            # dRc/dr [d(g-r)/dr = -1]
    dRc_dg   = dRc_dgr                  # dRc/dg [d(g-r)/dg = +1]
    sRc_phot = math.sqrt((abs(dRc_dg) * eg)**2 + (dRc_dr * er)**2)
    # E(B-V) uncertainty for Rc: differential term
    # d(Rc)/d(EBV) = R_Rc - R_r - A1RC*(R_g - R_r)
    # where R_Rc=2.32 (Cousins R, Cardelli 1989)
    # Note: do NOT use R_r*dRc_dr + R_g*dRc_dg (those are absolute extinctions)
    delta_RRc = abs(R_EXT["Rc"] - R_EXT["r"]
                   - A1RC * (R_EXT["g"] - R_EXT["r"]))
    sRc_ebv  = delta_RRc * ebv
    sRc      = math.sqrt(sRc_phot**2 + RMS_RC**2 + sRc_ebv**2)

    # ── Ic band [v4.0 confirmed — 4-parameter model, N=281 stars, provisional] ──
    # Ic = i + A0IC + A1IC*(r-i) + A2IC*(r-i)^2 + A3IC_GR*(g-r)
    # Derivate (d(r-i)/di=-1, d(r-i)/dr=+1, d(g-r)/dg=+1, d(g-r)/dr=-1):
    #   dIc/di = 1 - A1IC - 2*A2IC*(r-i)
    #   dIc/dr = A1IC + 2*A2IC*(r-i) - A3IC_GR
    #   dIc/dg = +A3IC_GR
    if i is not None:
        ri         = r - i
        Ic         = i + A0IC + A1IC*ri + A2IC*ri**2 + A3IC_GR*gr
        dIc_di     = 1.0 - A1IC - 2.0*A2IC*ri
        dIc_dr     = A1IC + 2.0*A2IC*ri - A3IC_GR
        dIc_dg     = A3IC_GR
        sIc_phot   = math.sqrt(
                         (dIc_di * ei)**2 +
                         (abs(dIc_dr) * er)**2 +
                         (abs(dIc_dg) * eg)**2)
        # E(B-V) uncertainty for Ic: differential term (linear approximation)
        delta_RIc  = abs(R_EXT["Ic"] - R_EXT["i"]
                        - A1IC * (R_EXT["r"] - R_EXT["i"]))
        sIc_ebv    = delta_RIc * ebv
        sIc        = math.sqrt(sIc_phot**2 + RMS_IC**2 + sIc_ebv**2)
    else:
        Ic = sIc = ri = None

    return {
        "V":  round(V,  4),  "sV":  round(sV,  4),
        "B":  round(B,  4),  "sB":  round(sB,  4),
        "Rc": round(Rc, 4),  "sRc": round(sRc, 4),
        "Ic": round(Ic, 4) if Ic is not None else None,
        "sIc":round(sIc,4)  if sIc is not None else None,
        "gr": round(gr, 4),
        "ri": round(ri, 4) if ri is not None else None,
        "in_range": in_range,
        "warning":  warn,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  INPUT CSV FILE READER
# ══════════════════════════════════════════════════════════════════════════════

def read_input_csv(path: str) -> list[dict]:
    """
    Read the input CSV file of stars to process.

    Expected format (header optional):
        Name, RA, DEC, g, eg, r, er, i, ei, EBV

    Parsing rules:
    - Delimiter: comma or semicolon (auto-detected)
    - Lines starting with # are comments
    - Columns from g onward are optional
    - First row is treated as header if it contains 'name'/'ra' etc.
    - Coordinates accept spaces or colons as internal separators
    """
    # Auto-detect delimiter
    with open(path, encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
    sep = ";" if sample.count(";") > sample.count(",") else ","

    rows = []
    with open(path, encoding="utf-8-sig") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(sep)]
            rows.append(parts)

    if not rows:
        return []

    # Remove header row if present
    first = rows[0]
    header_keywords = ("name","star","ar","ra","asc",
                       "coord","id","label","object")
    if first and first[0].lower().strip() in header_keywords:
        rows = rows[1:]
    # Header with RA as second column
    elif len(first) >= 2 and first[0].lower() in header_keywords:
        rows = rows[1:]

    stars = []
    for line_num, parts in enumerate(rows, start=2):
        # RA and Dec are required.
        # Format A: Name, RA, DEC, ...   (≥ 3 tokens)
        # Format B: RA, DEC, ...         (≥ 2 tokens, first column is coordinates)

        def _try_coord(ra_s, dec_s):
            try:
                return parse_radec(ra_s, dec_s)
            except Exception:
                return None

        coord = None
        col_offset = 0   # quante colonne "consumano" le coordinate

        if len(parts) >= 3:
            # Try Name | RA | DEC
            coord = _try_coord(parts[1], parts[2])
            if coord is not None:
                name = parts[0] or f"Star_{line_num}"
                col_offset = 3
            else:
                # Try RA | DEC (no name column)
                coord = _try_coord(parts[0], parts[1])
                if coord is not None:
                    name = f"Star_{line_num}"
                    col_offset = 2
        elif len(parts) == 2:
            coord = _try_coord(parts[0], parts[1])
            if coord is not None:
                name = f"Star_{line_num}"
                col_offset = 2

        if coord is None:
            print(f"  {YELL}Row {line_num}: unreadable coordinates "
                  f"({parts[:3]}) — skipped{RST}")
            continue

        # Read optional columns (g, eg, r, er, i, ei, EBV)
        opt = parts[col_offset:]   # colonne rimanenti

        def _num(idx):
            """Extract float from position idx; None if absent/empty/non-numeric."""
            if idx < len(opt):
                s = opt[idx].strip()
                if s in ("", "-", "--", "n/a", "N/A", "none", "None", "nan"):
                    return None
                try:
                    return float(s.replace(",", "."))
                except ValueError:
                    return None
            return None

        stars.append({
            "name": name,
            "coord": coord,
            "g":   _num(0),
            "eg":  _num(1),
            "r":   _num(2),
            "er":  _num(3),
            "i":   _num(4),
            "ei":  _num(5),
            "ebv": _num(6),
        })

    print(f"  {GREEN}Read {len(stars)} valid stars{RST} from file '{path}'")
    return stars


# ══════════════════════════════════════════════════════════════════════════════
#  SINGLE STAR PROCESSING PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def process_star(name: str,
                 coord: SkyCoord,
                 g=None, eg=None, r=None, er=None,
                 i=None, ei=None, ebv=None,
                 auto_apass: bool = True) -> Optional[dict]:
    """
    Full processing pipeline for a single star:
      1. Retrieve g,r,i from APASS DR9 if missing and auto_apass=True
      2. Fetch E(B-V) from IRSA (or analytic fallback)
      3. Compute V, B, Rc, Ic with full error propagation
      4. Return the result dict
    """
    ra_s, dec_s = coord_hmsdms(coord)
    print(f"\n{BOLD}  ► {CYAN}{name}{RST}  "
          f"(RA={ra_s}  Dec={dec_s})")

    # ── 1. APASS photometry ──────────────────────────────────────────────────
    need_apass = (g is None or r is None)
    if need_apass:
        if not auto_apass:
            print(f"    {RED}Fotometria mancante e download automatico "
                  f"disabilitato → star skipped{RST}")
            return None
        phot = fetch_apass(coord, verbose=True)
        if phot is None:
            print(f"    {RED}Cannot retrieve APASS photometry → "
                  f"star skipped{RST}")
            return None
        g  = g  or phot.get("g")
        eg = eg or phot.get("eg")
        r  = r  or phot.get("r")
        er = er or phot.get("er")
        i  = i  or phot.get("i")
        ei = ei or phot.get("ei")
    else:
        print(f"    APASS photometry from CSV: "
              f"g={g:.3f}(±{eg:.3f})  "
              f"r={r:.3f}(±{er:.3f})  "
              f"i={f'{i:.3f}' if i else 'N/A'}"
              f"{f'(±{ei:.3f})' if ei else ''}")

    if g is None or r is None:
        print(f"    {RED}g or r not available → star skipped{RST}")
        return None

    # ── 2. E(B-V) ────────────────────────────────────────────────────────────
    if ebv is None:
        ebv, ebv_method = get_ebv(coord, name)
    else:
        ebv_method = "CSV input"
        print(f"    E(B-V) = {ebv:.4f} [from CSV]")

    # ── 3. Compute Johnson magnitudes ────────────────────────────────────────
    res = compute_johnson(g, r, i, eg, er, ei, ebv, ebv_method)

    # ── 4. High E(B-V) warning ───────────────────────────────────────────────
    ebv_warn = ""
    if ebv >= 0.50:
        ebv_warn = "NOT recommended as comparison star"
        print(f"    {RED}⚠ E(B-V) = {ebv:.4f} mag — ELEVATO (≥0.50): estinzione "
              f"alta e variabile spazialmente.{RST}")
        print(f"    {RED}  Questa star è SCONSIGLIATA come riferimento fotometrico.{RST}")
    elif ebv >= 0.20:
        ebv_warn = "high extinction — use with caution"
        print(f"    {YELL}⚠ E(B-V) = {ebv:.4f} mag — ALTO (≥0.20): verificare che "
              f"le stars vicine abbiano estinzione simile.{RST}")

    # ── 5. Print result to screen ────────────────────────────────────────────
    _print_star_result(name, coord, g, eg, r, er, i, ei, ebv, ebv_method, res)

    return {
        # Identification
        "Name":       name,
        "RA":         ra_s,
        "DEC":        dec_s,
        # APASS input photometry
        "g":          round(g, 4),
        "eg":         round(eg, 4) if eg else "",
        "r":          round(r, 4),
        "er":         round(er, 4) if er else "",
        "i":          round(i, 4)  if i  else "",
        "ei":         round(ei, 4) if ei else "",
        "E_BV":       round(ebv, 4),
        "ebv_method":    ebv_method,
        # Colour indices
        "g-r":        res["gr"],
        "r-i":        res["ri"] if res["ri"] is not None else "",
        # Johnson-Cousins output
        "V":          res["V"],
        "sV":         res["sV"],
        "B":          res["B"],
        "sB":         res["sB"],
        "Rc":         res["Rc"],
        "sRc":        res["sRc"],
        "Ic":         res["Ic"]  if res["Ic"]  is not None else "",
        "sIc":        res["sIc"] if res["sIc"] is not None else "",
        # Quality flags
        "in_range":   "YES" if res["in_range"] else "NO",
        "warning":    " | ".join(filter(None, [res["warning"], ebv_warn])),
        "ebv_warn":   ebv_warn,
        # Sources
        "source_VB":  "N. Montigiani (2026) v5.0 — N=281 stars",
        "source_Rc":   "N. Montigiani (2026) v4.0 — N=281 stars (confirmed)",
        "source_Ic":   "N. Montigiani (2026) v4.0 — N=281 stars (confirmed)",
    }


def _print_star_result(name, coord, g, eg, r, er, i, ei,
                       ebv, ebv_method, res):
    """Print formatted results box for a star."""
    W = 62
    print(f"    {'─'*W}")
    print(f"    {BOLD}RESULTS  —  {CYAN}{name}{RST}")
    print(f"    {DIM}g={g:.3f}±{eg:.3f}  r={r:.3f}±{er:.3f}"
          f"{'  i='+f'{i:.3f}'+'±'+f'{ei:.3f}' if i else '  i=N/A'}"
          f"  E(B-V)={ebv:.4f} [{ebv_method}]{RST}")
    print(f"    {DIM}(g-r)={res['gr']:.3f}"
          f"{'  (r-i)='+str(res['ri']) if res['ri'] is not None else ''}"
          f"{'  ⚠ EXTRAPOLATO' if not res['in_range'] else ''}{RST}")
    print(f"    {'─'*W}")
    hdr = f"    {'Band':<8}{'Magnitude':>13}{'± σ_tot':>10}   Source"
    print(hdr)
    print(f"    {'─'*W}")

    rows = [
        ("V",    res["V"],  res["sV"],  "N. Montigiani (2026) v5.0  N=281"),
        ("B",    res["B"],  res["sB"],  "N. Montigiani (2026) v5.0  N=281"),
        ("Rc",   res["Rc"], res["sRc"], "N. Montigiani (2026) v4.0  N=281"),
        ("Ic",
         res["Ic"]  if res["Ic"]  is not None else None,
         res["sIc"] if res["sIc"] is not None else None,
         "N. Montigiani (2026) v4.0  N=281"),
    ]
    for band, mag, sig, fonte in rows:
        if mag is None:
            print(f"    {band:<8}{'N/A (i mancante)':>13}"
                  f"{'---':>10}   {DIM}{fonte}{RST}")
        else:
            flag = (f"  {YELL}⚠{RST}" if band in ("Ic",) and sig > 0.15
                    else "")
            print(f"    {band:<8}{mag:>13.4f}{sig:>10.4f}"
                  f"   {DIM}{fonte}{RST}{flag}")
    print(f"    {'─'*W}")
    if res["warning"]:
        print(f"    {YELL}⚠ {res['warning']}{RST}")


# ══════════════════════════════════════════════════════════════════════════════
#  SALVATAGGIO CSV OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

OUTPUT_COLS = [
    "Name", "RA", "DEC",
    "V", "sV", "B", "sB", "Rc", "sRc", "Ic", "sIc",
    "g-r", "r-i",
    "g", "eg", "r", "er", "i", "ei",
    "E_BV", "ebv_method",
    "in_range", "ebv_warn", "warning",
    "source_VB", "source_Rc", "source_Ic",
]


def save_results(records: list[dict], outfile: str = "johnson_results.csv"):
    """
    Save results to CSV with comma delimiter and header.
    Magnitude and uncertainty columns formatted to 4 decimal places.
    Encoded as UTF-8 with BOM for Excel compatibility.
    """
    if not records:
        print(f"\n  {YELL}No results to save.{RST}")
        return

    df = pd.DataFrame(records)
    cols = [c for c in OUTPUT_COLS if c in df.columns]
    df   = df[cols]

    # Format magnitude columns
    mag_cols = ["V","sV","B","sB","Rc","sRc","Ic","sIc",
                "g","eg","r","er","i","ei","E_BV","g-r","r-i"]
    for c in mag_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df.to_csv(outfile, index=False, float_format="%.4f",
              encoding="utf-8-sig", na_rep="")

    print(f"\n{GREEN}{BOLD}  ✓ File saved: '{outfile}'  "
          f"({len(df)} stars){RST}")

    # ── Summary table to screen ──────────────────────────────────────────────
    print(f"\n{BOLD}  {'═'*78}")
    print(f"  SUMMARY — JOHNSON-COUSINS MAGNITUDES")
    print(f"  {'═'*78}{RST}")

    hdr = (f"  {'Name':<18}"
           f"{'V':>8}{'±σV':>7}"
           f"{'B':>8}{'±σB':>7}"
           f"{'Rc':>8}{'±σRc':>7}"
           f"{'Ic':>8}{'±σIc':>7}"
           f"  {'(g-r)':>6}"
           f"  {'⚠' if any(r['in_range']=='NO' for r in records) else ''}")
    print(hdr)
    print(f"  {'─'*78}")

    for rec in records:
        def fmt(val, ndec=4):
            if val == "" or val is None:
                return "  N/A "
            try:
                return f"{float(val):>8.{ndec}f}"
            except (ValueError, TypeError):
                return f"{'N/A':>8}"

        flag = f"  {YELL}⚠{RST}" if rec.get("in_range") == "NO" else ""
        print(f"  {rec['Name']:<18}"
              f"{fmt(rec['V'])}{fmt(rec['sV'])}"
              f"{fmt(rec['B'])}{fmt(rec['sB'])}"
              f"{fmt(rec['Rc'])}{fmt(rec['sRc'])}"
              f"{fmt(rec['Ic'])}{fmt(rec['sIc'])}"
              f"  {fmt(rec.get('g-r',''),3)}"
              f"{flag}")

    print(f"  {'─'*78}")
    print(f"  {DIM}σ includes: APASS photometric error + calibration RMS + "
          f"E(B-V) uncertainty{RST}")
    print(f"  {DIM}⚠ = (g-r) outside valid range [-0.28, +1.10]: extrapolated value{RST}")


# ══════════════════════════════════════════════════════════════════════════════
#  INTERFACCIA UTENTE E MAIN
# ══════════════════════════════════════════════════════════════════════════════

def ask_ebv_interactive(name: str, coord: SkyCoord) -> tuple[float, str]:
    """Prompt user to enter E(B-V) manually or retrieve from IRSA."""
    ra_s, dec_s = coord_hmsdms(coord)
    b = coord.galactic.b.deg
    print(f"\n    Stella: {CYAN}{name}{RST}  RA={ra_s}  Dec={dec_s}  b={b:.1f}°")
    print(f"    Enter E(B-V) for this star:")
    print(f"      [Enter]  = retrieve automatically from IRSA/SFD")
    print(f"      [value]  = enter manually (e.g. 0.032)")

    while True:
        ans = input("    E(B-V) = ").strip()
        if ans == "":
            return get_ebv(coord, name)
        try:
            val = float(ans.replace(",", "."))
            if 0.0 <= val <= 5.0:
                print(f"    E(B-V) = {val:.4f} [entered manually]")
                return val, "manual"
            print(f"    {YELL}Value out of range [0, 5]. Try again.{RST}")
        except ValueError:
            print(f"    {YELL}Unrecognised input. Enter a number "
                  f"(e.g. 0.045) or press Enter.{RST}")


def print_banner():
    print(f"""
{CYAN}{BOLD}╔══════════════════════════════════════════════════════════════╗
║      apass2johnson  v5.0  —  CCD Photometric Transformation  ║
║   SDSS (g,r,i) → Johnson-Cousins (V, B, Rc, Ic) + errors    ║
╚══════════════════════════════════════════════════════════════╝{RST}
  Formulae — N. Montigiani (2026):
    V  → v5.0  N=281  RMS=0.032  (+9% vs Lupton 2005)
    B  → v5.0  N=281  RMS=0.058  (+40% vs Lupton (2005))
    Rc → v4.0  N=281  RMS=0.051  (+16% vs Lupton 2005)  [confirmed]
    Ic → v4.0  N=281  RMS=0.083  (+4%  vs Lupton 2005)  [confirmed, provisional]
  Errors:   photometric propagation + calibration RMS + extinction
""")


def main():
    print_banner()

    # ── Select input CSV file ──────────────────────────────────────────────────
    while True:
        csv_path = input(
            "  Input CSV file with stars to process\n"
            "  (e.g. targets.csv  or  /full/path/targets.csv): "
        ).strip().strip('"').strip("'")

        if os.path.isfile(csv_path):
            break
        print(f"  {RED}File not found: '{csv_path}'{RST}  Try again.\n")

    print()
    stars = read_input_csv(csv_path)
    if not stars:
        print(f"  {RED}No valid stars found in the file. Exiting.{RST}")
        sys.exit(1)

    # ── APASS download mode ────────────────────────────────────────────────────
    n_missing = sum(1 for s in stars if s["g"] is None or s["r"] is None)

    if n_missing > 0:
        print(f"\n  {YELL}{n_missing}/{len(stars)} stars{RST} have no "
              f"g,r,i photometry in the CSV.")
        if not ASTROQUERY_OK:
            print(f"  {RED}astroquery not installed: cannot download "
                  f"from APASS DR9.{RST}")
            print(f"  Install with:  pip install astroquery")
            print(f"  Or provide g,r,i directly in the CSV.")
            auto_apass = False
        else:
            print(f"\n  Automatic download from APASS DR9 (VizieR)?")
            print(f"  {BOLD}[Y] or [1]{RST} Yes, download automatically (requires internet)")
            print(f"  {BOLD}[N] or [2]{RST} No  (stars without photometry will be skipped)")
            while True:
                ans = input("  Choice [Y/N]: ").strip().upper()
                if ans in ("Y", "1", "YES"):
                    auto_apass = True
                    break
                elif ans in ("N", "2", "NO"):
                    auto_apass = False
                    break
                print(f"  {YELL}Inserire S (sì) or N (no).{RST}")
    else:
        auto_apass = False   # all photometry already present in the CSV

    # ── E(B-V) retrieval mode ──────────────────────────────────────────────────
    n_noebv = sum(1 for s in stars if s["ebv"] is None)
    ebv_mode = "auto"   # "auto" = IRSA; "ask" = chiedi per ogni star

    if n_noebv > 0:
        print(f"\n  {n_noebv}/{len(stars)} stars have no E(B-V) in the CSV.")
        print(f"  How to obtain E(B-V)?")
        print(f"  {BOLD}[1]{RST} Automatic — NASA/IPAC IRSA (SFD online) "
              f"or analytic approximation")
        print(f"  {BOLD}[2]{RST} Ask for each star (manual entry)")
        ec = input("  Choice [1/2]: ").strip()
        ebv_mode = "ask" if ec == "2" else "auto"

    # ── Output file name ─────────────────────────────────────────────────────
    default_out = "johnson_results.csv"
    out_name = input(f"\n  Output filename [{default_out}]: ").strip()
    if not out_name:
        out_name = default_out

    # ── Processing ───────────────────────────────────────────────────────────
    print(f"\n{BOLD}{'═'*66}")
    print(f"  PROCESSING  ({len(stars)} stars)")
    print(f"{'═'*66}{RST}")

    records = []
    for idx, s in enumerate(stars, 1):
        print(f"\n  [{idx}/{len(stars)}]", end="")

        # E(B-V): usa CSV, or chiede/download secondo la modalità
        ebv_val = s["ebv"]
        if ebv_val is None:
            if ebv_mode == "ask":
                ebv_val, _ = ask_ebv_interactive(s["name"], s["coord"])
            # se ebv_mode == "auto", viene gestito dentro process_star → get_ebv

        rec = process_star(
            name       = s["name"],
            coord      = s["coord"],
            g          = s["g"],   eg = s["eg"],
            r          = s["r"],   er = s["er"],
            i          = s["i"],   ei = s["ei"],
            ebv        = ebv_val,
            auto_apass = auto_apass,
        )
        if rec:
            records.append(rec)

        # Courtesy pause for VizieR rate limits
        if auto_apass and idx < len(stars):
            time.sleep(0.4)

    # ── Final output ─────────────────────────────────────────────────────────
    save_results(records, out_name)

    print(f"""
  {BOLD}USAGE NOTES:{RST}
  • V, B, Rc, Ic are in the Johnson-Cousins system, ready to be entered
    as comparison stars in MaxIm DL, AstroImageJ, APT, Muniwin, etc.
  • σ_tot includes: APASS photometric error + calibration RMS +
    residual extinction uncertainty.
  • V and B (v5.0): re-derived on N=281 validation stars from four
    independent catalogues (Landolt 1992/2013, Clem & Landolt 2013,
    Stetson/Pancino 2022). +9% improvement for V; +40% for B vs Lupton (2005).
  • Rc (v4.0, confirmed): RMS = 0.051 mag, validated on N=281 stars.
  • Ic (v4.0, provisional): RMS = 0.083 mag (4-parameter model with (r-i)^2
    and (g-r) covariate). Use with caution for precision < 0.05 mag.
  • For fields with E(B-V) > 0.20 mag, consider applying extinction
    corrections (Version A) before using comparison stars.
""")


if __name__ == "__main__":
    main()

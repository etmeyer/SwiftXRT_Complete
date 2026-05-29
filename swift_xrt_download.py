#!/usr/bin/env python3
"""
Swift XRT Archive Data Downloader
==================================
Downloads observation data from the Swift XRT archive based on a source name
(resolved to coordinates via SIMBAD), explicit RA/Dec position, a single
observation ID, or a file containing a list of observation IDs.

Uses the HEASARC Browse interface to query the Swift Master Catalog
and downloads data products from the HEASARC archive.

Features:
    - Name resolution via SIMBAD, NED, CDS Sesame, and astroquery (chained)
    - Direct obsid input: single ID or batch from a file
    - Test mode: download only the first N observations (--test N)
    - Resume support: skips files that already exist with the correct size
    - Overwrite mode: force re-download of existing files (--overwrite)
    - Date-window filtering of catalog results (--start-date / --end-date)
    - Progress tracker with ETA for multi-observation downloads

Requirements:
    pip install requests
    pip install astroquery astropy   # optional, fallback for name resolution

Usage:
    # By source name
    python swift_xrt_download.py --name "Crab Nebula" --radius 12

    # By coordinates (decimal degrees)
    python swift_xrt_download.py --ra 83.633 --dec 22.014 --radius 12

    # Single observation ID
    python swift_xrt_download.py --obsid 00035393001

    # Batch download from a file of obsids (one per line)
    python swift_xrt_download.py --obsid-file my_obsids.txt

    # Test mode: grab only the first 3 observations
    python swift_xrt_download.py --name "GRS 1915+105" --test 3

    # Limit number of observations and choose output directory
    python swift_xrt_download.py --name "GRS 1915+105" --max-obs 5 --outdir ./data

    # Force re-download even if local files exist
    python swift_xrt_download.py --name "Cyg X-1" --test 2 --overwrite

    # Only PC mode, cleaned event files
    python swift_xrt_download.py --name "Mrk 421" --mode pc --clean-only

    # Only WT and PC modes (skip Image mode)
    python swift_xrt_download.py --obsid 00035393001 --mode wt pc

    # List observations without downloading
    python swift_xrt_download.py --name "Cyg X-1" --list-only

    # List only observations in a date window (end is exclusive)
    python swift_xrt_download.py --name "3C 273" --list-only \
        --start-date 2008-08-04 --end-date 2011-07-06
"""

import argparse
import math
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import date, timedelta
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Name resolution
# ---------------------------------------------------------------------------

def resolve_name_simbad(name: str) -> tuple[float, float]:
    """Resolve via the SIMBAD TAP service (robust, handles most names)."""
    url = "https://simbad.u-strasbg.fr/simbad/sim-id"
    params = {
        "Ident": name,
        "output.format": "votable",
        "output.params": "main_id,ra(d),dec(d)",
    }
    print(f"[info] Resolving '{name}' via SIMBAD …")
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    # Parse the VOTable response
    root = ET.fromstring(resp.text)
    ns_candidates = [
        "http://www.ivoa.net/xml/VOTable/v1.3",
        "http://www.ivoa.net/xml/VOTable/v1.2",
        "http://www.ivoa.net/xml/VOTable/v1.1",
        "",
    ]
    for ns_uri in ns_candidates:
        ns = {"v": ns_uri} if ns_uri else {}
        prefix = "v:" if ns_uri else ""
        fields = root.findall(f".//{prefix}FIELD", ns)
        if fields:
            break
    else:
        raise RuntimeError(f"SIMBAD returned no parseable VOTable for '{name}'")

    col_names = [f.attrib.get("name", f.attrib.get("ID", "")).lower()
                 for f in fields]
    tds = root.findall(f".//{prefix}TR/{prefix}TD", ns)
    if not tds or len(tds) < len(col_names):
        raise RuntimeError(f"SIMBAD returned no results for '{name}'")

    values = {c: td.text.strip() if td.text else ""
              for c, td in zip(col_names, tds)}

    # SIMBAD uses column names like "ra_d" / "dec_d" or "ra(d)" / "dec(d)"
    ra_str = values.get("ra_d") or values.get("ra(d)", "")
    dec_str = values.get("dec_d") or values.get("dec(d)", "")
    if not ra_str or not dec_str:
        raise RuntimeError(f"SIMBAD result for '{name}' missing coordinates "
                           f"(columns: {list(values.keys())})")

    ra, dec = float(ra_str), float(dec_str)
    print(f"[info] SIMBAD resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
    return ra, dec


def resolve_name_ned(name: str) -> tuple[float, float]:
    """Resolve via the NASA/IPAC Extragalactic Database (NED)."""
    url = "https://ned.ipac.caltech.edu/srs/ObjectLookup"
    params = {"name": name, "of": "json"}
    print(f"[info] Resolving '{name}' via NED …")
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    # NED returns a "ResultCode" of 3 for success, with position in "Preferred"
    result_code = data.get("ResultCode", 0)
    if result_code != 3:
        raise RuntimeError(f"NED could not resolve '{name}' "
                           f"(ResultCode={result_code})")

    pos = data.get("Preferred", {}).get("Position", {})
    ra = pos.get("RA")
    dec = pos.get("Dec")
    if ra is None or dec is None:
        raise RuntimeError(f"NED result for '{name}' missing coordinates")

    ra, dec = float(ra), float(dec)
    print(f"[info] NED resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
    return ra, dec


def resolve_name_sesame(name: str) -> tuple[float, float]:
    """Fallback resolver using CDS Sesame (queries SIMBAD, NED, VizieR)."""
    url = "https://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-ox/SNV"
    params = {"obj": name}
    print(f"[info] Resolving '{name}' via CDS Sesame …")
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)
    ns = {"s": "http://vizier.u-strasbg.fr/xml/sesame_xml.xsd"}
    resolver = root.find(".//s:Resolver", ns)
    if resolver is None:
        raise RuntimeError(f"Sesame returned no result for '{name}'")

    jradeg = resolver.find("s:jradeg", ns)
    jdedeg = resolver.find("s:jdedeg", ns)
    if jradeg is None or jdedeg is None:
        raise RuntimeError(f"Sesame result for '{name}' missing coordinates")

    ra, dec = float(jradeg.text), float(jdedeg.text)
    print(f"[info] Sesame resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
    return ra, dec


def resolve_name_astroquery(name: str) -> tuple[float, float]:
    """Last-resort resolver using astroquery (if installed)."""
    try:
        from astroquery.simbad import Simbad
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        print(f"[info] Resolving '{name}' via astroquery/SIMBAD …")
        result = Simbad.query_object(name)
        if result is None:
            raise RuntimeError(f"Could not resolve source name '{name}'")

        # astroquery ≥0.4.8 uses lowercase "ra"/"dec" (degrees);
        # older versions use uppercase "RA"/"DEC" (sexagesimal).
        colnames = [c.lower() for c in result.colnames]
        if "ra" in colnames and "dec" in colnames:
            ra_col = result.colnames[colnames.index("ra")]
            dec_col = result.colnames[colnames.index("dec")]
            ra_val, dec_val = result[ra_col][0], result[dec_col][0]
            # New astroquery returns degrees directly; old returns strings
            try:
                ra, dec = float(ra_val), float(dec_val)
            except (TypeError, ValueError):
                coord = SkyCoord(ra_val, dec_val, unit=(u.hourangle, u.deg))
                ra, dec = coord.ra.deg, coord.dec.deg
        else:
            raise RuntimeError(
                f"Unexpected columns in astroquery result: {result.colnames}")

        print(f"[info] astroquery resolved to RA={ra:.5f}°, Dec={dec:.5f}°")
        return ra, dec
    except ImportError:
        raise RuntimeError("astroquery is not installed")


# Ordered chain of resolvers – tried in sequence until one succeeds
_RESOLVERS = [
    ("SIMBAD", resolve_name_simbad),
    ("NED", resolve_name_ned),
    ("Sesame", resolve_name_sesame),
    ("astroquery", resolve_name_astroquery),
]


def resolve_name(name: str) -> tuple[float, float]:
    """Try each resolver in turn; return the first successful result."""
    errors: list[str] = []
    for label, func in _RESOLVERS:
        try:
            return func(name)
        except Exception as exc:
            print(f"[warn] {label} resolution failed: {exc}")
            errors.append(f"{label}: {exc}")
    raise RuntimeError(
        f"All resolvers failed for '{name}':\n  " + "\n  ".join(errors)
    )


# ---------------------------------------------------------------------------
# HEASARC query
# ---------------------------------------------------------------------------

HEASARC_CONESEARCH_URL = "https://heasarc.gsfc.nasa.gov/xamin/vo/cone"
HEASARC_TAP_URL = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
HEASARC_BROWSE_URL = "https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3query.pl"


def query_swift_master(ra: float, dec: float, radius_arcmin: float = 12.0) -> list[dict]:
    """
    Query the HEASARC swiftmastr catalog for Swift observations near (ra, dec).

    Tries multiple strategies in order:
      1. VO Cone Search (simple, standard, positional by design)
      2. W3Browse BatchDisplay text output
      3. HEASARC TAP with coordinate box

    Parameters
    ----------
    ra, dec : float
        Position in decimal degrees (J2000).
    radius_arcmin : float
        Search cone radius in arcminutes (default 12′, roughly the XRT FOV).

    Returns
    -------
    list of dict
        Each dict contains observation metadata fields.
    """
    strategies = [
        ("VO Cone Search",  lambda: _query_via_conesearch(ra, dec, radius_arcmin)),
        ("W3Browse (text)", lambda: _query_via_w3browse_text(ra, dec, radius_arcmin)),
        ("TAP (box)",       lambda: _query_via_tap_box(ra, dec, radius_arcmin)),
    ]
    for label, func in strategies:
        try:
            result = func()
            if result:
                return result
            print(f"[warn] {label} returned 0 results, trying next strategy …")
        except Exception as exc:
            print(f"[warn] {label} failed ({exc}), trying next strategy …")

    # All strategies exhausted
    return []


# --- Strategy 1: VO Cone Search (standard IVOA protocol) -------------------

def _query_via_conesearch(ra: float, dec: float,
                           radius_arcmin: float) -> list[dict]:
    """Query using the IVOA Simple Cone Search protocol.

    This is the most straightforward positional query — RA, Dec, and
    search radius are the only parameters.  Returns VOTable XML.
    """
    radius_deg = radius_arcmin / 60.0
    params = {
        "table": "swiftmastr",
        "RA": str(ra),
        "DEC": str(dec),
        "SR": str(radius_deg),
    }
    print(f"[info] Querying HEASARC VO Cone Search (swiftmastr) within "
          f"{radius_arcmin}′ of RA={ra:.5f}, Dec={dec:.5f} …")
    resp = requests.get(HEASARC_CONESEARCH_URL, params=params, timeout=90)
    resp.raise_for_status()
    return _parse_votable(resp.text)


# --- Strategy 2: W3Browse text output (pipe / plus-delimited) ---------------

def _query_via_w3browse_text(ra: float, dec: float,
                              radius_arcmin: float) -> list[dict]:
    """Query using W3Browse with pipe-delimited text output."""
    # Format coordinates explicitly as decimal degrees for W3Browse.
    # The sign on Dec must be explicit; append 'd' to each value.
    dec_str = f"+{dec}" if dec >= 0 else f"{dec}"
    params = {
        "tablehead": "name=heasarc_swiftmastr&description=Swift Master Catalog",
        "Action": "Query",
        "Coordinates": f"{ra}d {dec_str}d",
        "Equinox": "2000",
        "Radius": str(radius_arcmin),
        "Radius_unit": "arcmin",
        "NR": "CheckCaches/GRB/SIMBAD/NED",
        "ResultMax": "10000",
        "displaymode": "BatchDisplay",
        "Fields": "All",
        "vession": "img",
        "gifsize": "0",
    }
    print(f"[info] Querying HEASARC W3Browse (text) within {radius_arcmin}′ of "
          f"RA={ra:.5f}, Dec={dec:.5f} …")
    resp = requests.get(HEASARC_BROWSE_URL, params=params, timeout=60)
    resp.raise_for_status()
    return _parse_batch_text(resp.text)


def _parse_batch_text(text: str) -> list[dict]:
    """Parse the pipe-delimited BatchDisplay output from W3Browse.

    Handles both separator styles::

        col1|col2|col3          (pipes only)
        ---|---|---

        |col1|col2|col3|        (pipes with leading/trailing)
        +----+----+----+        (plus-delimited separator)
    """
    lines = text.strip().splitlines()

    # Find the separator row.  It looks like one of:
    #   ---|---|---        (dashes and pipes)
    #   +------+------+   (dashes and plus signs)
    header_idx = None
    data_start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        # A separator line is composed entirely of '-', '|', '+', and spaces
        if (stripped
                and all(c in "-|+ " for c in stripped)
                and "---" in stripped
                and i > 0):
            header_idx = i - 1
            data_start = i + 1
            break

    if header_idx is None or header_idx < 0:
        preview = "\n".join(lines[:10])
        raise RuntimeError(
            f"Could not find table header in W3Browse output:\n{preview}")

    col_names = [c.strip() for c in lines[header_idx].split("|") if c.strip()]

    observations = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.startswith("<") or line.startswith("Search"):
            continue
        # Stop at any trailing footer
        if line.startswith("BatchEnd") or line.startswith("***"):
            break
        values = [v.strip() for v in line.split("|")]
        # Strip empty tokens from leading/trailing pipes
        if values and values[0] == "":
            values = values[1:]
        if values and values[-1] == "":
            values = values[:-1]
        if len(values) < len(col_names):
            continue
        row = {c: values[j] for j, c in enumerate(col_names)}
        observations.append(row)

    print(f"[info] Found {len(observations)} observation(s).")
    return observations


# --- Strategy 3: TAP with coordinate box ------------------------------------

def _query_via_tap_box(ra: float, dec: float,
                        radius_arcmin: float) -> list[dict]:
    """Query HEASARC TAP using a simple coordinate-range WHERE clause."""
    radius_deg = radius_arcmin / 60.0
    cos_dec = math.cos(math.radians(dec))
    ra_range = radius_deg / max(cos_dec, 0.01)

    adql = (
        f"SELECT * FROM swiftmastr "
        f"WHERE ra BETWEEN {ra - ra_range} AND {ra + ra_range} "
        f"AND dec BETWEEN {dec - radius_deg} AND {dec + radius_deg}"
    )
    params = {
        "request": "doQuery",
        "version": "1.0",
        "lang": "ADQL",
        "format": "votable",
        "query": adql,
    }
    print(f"[info] Querying HEASARC TAP (box) (swiftmastr) within "
          f"{radius_arcmin}′ of RA={ra:.5f}, Dec={dec:.5f} …")
    resp = requests.get(HEASARC_TAP_URL, params=params, timeout=60)
    resp.raise_for_status()
    return _parse_votable(resp.text)


def _parse_votable(xml_text: str) -> list[dict]:
    """Minimal VOTable parser – extracts TABLEDATA rows."""
    # Guard against non-XML responses (HTML error pages, plain text, etc.)
    stripped = xml_text.strip()
    if not stripped.startswith("<?xml") and not stripped.startswith("<VOTABLE") \
       and not stripped.startswith("<vo:VOTABLE"):
        # Show a useful snippet for debugging
        preview = stripped[:300].replace("\n", " ")
        raise RuntimeError(
            f"Expected VOTable XML but got unexpected response: {preview!r}…"
        )

    root = ET.fromstring(xml_text)

    # Check for an INFO element with an error
    for ns_uri in ["http://www.ivoa.net/xml/VOTable/v1.3",
                    "http://www.ivoa.net/xml/VOTable/v1.2",
                    "http://www.ivoa.net/xml/VOTable/v1.1", ""]:
        ns = {"v": ns_uri} if ns_uri else {}
        prefix = "v:" if ns_uri else ""
        for info in root.findall(f".//{prefix}INFO", ns):
            if info.attrib.get("name") == "QUERY_STATUS" \
               and info.attrib.get("value") == "ERROR":
                msg = info.text.strip() if info.text else "unknown error"
                raise RuntimeError(f"HEASARC query error: {msg}")

    # VOTable namespace
    ns_candidates = [
        "http://www.ivoa.net/xml/VOTable/v1.3",
        "http://www.ivoa.net/xml/VOTable/v1.2",
        "http://www.ivoa.net/xml/VOTable/v1.1",
        "",
    ]
    for ns_uri in ns_candidates:
        ns = {"v": ns_uri} if ns_uri else {}
        prefix = "v:" if ns_uri else ""
        fields = root.findall(f".//{prefix}FIELD", ns)
        if fields:
            break
    else:
        print("[warn] No FIELD elements found in VOTable – "
              "the query may have returned no results.")
        return []

    col_names = [f.attrib.get("name", f"col{i}") for i, f in enumerate(fields)]
    rows_el = root.findall(f".//{prefix}TR", ns)
    observations = []
    for tr in rows_el:
        tds = tr.findall(f"{prefix}TD", ns)
        row = {}
        for cname, td in zip(col_names, tds):
            row[cname] = td.text.strip() if td.text else ""
        observations.append(row)

    print(f"[info] Found {len(observations)} observation(s).")
    return observations


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

HEASARC_DATA_BASE = "https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs"
UKSSDC_DATA_BASE = "https://www.swift.ac.uk/archive/reproc"

# Cache: maps obsid → verified base URL to avoid repeated probing
_url_cache: dict[str, str] = {}
# Which URL strategy worked last — used to order subsequent attempts
_preferred_strategy: list[str] = ["heasarc_date", "ukssdc", "heasarc_flat"]


def resolve_data_url(obsid: str, obs_meta: dict | None = None) -> str | None:
    """Discover the correct archive base URL for a given observation.

    Tries multiple URL patterns against HEASARC and UKSSDC until one
    responds successfully.  Caches results and learns which strategy
    works to speed up subsequent lookups.

    Parameters
    ----------
    obsid : str
        Swift observation ID (e.g. '00035393001').
    obs_meta : dict or None
        Catalog row for this observation.  If available, date fields
        are used to construct HEASARC date-based URLs.

    Returns
    -------
    str or None
        The base URL (ending with '/') for the observation directory,
        or None if no working URL could be found.
    """
    global _preferred_strategy
    obsid = obsid.strip().zfill(11)

    if obsid in _url_cache:
        return _url_cache[obsid]

    # Build candidate URLs keyed by strategy name
    candidates: dict[str, list[str]] = {
        "heasarc_date": [],
        "ukssdc": [f"{UKSSDC_DATA_BASE}/{obsid}/"],
        "heasarc_flat": [f"{HEASARC_DATA_BASE}/{obsid}/"],
    }

    # Extract dates from catalog metadata for HEASARC date-based paths
    if obs_meta:
        seen_months = set()
        for field in ("start_time", "archive_date", "processing_date"):
            val = obs_meta.get(field, "").strip()
            if val and len(val) >= 7:
                try:
                    year_month = val[:4] + "_" + val[5:7]
                    if year_month not in seen_months:
                        seen_months.add(year_month)
                        candidates["heasarc_date"].append(
                            f"{HEASARC_DATA_BASE}/{year_month}/{obsid}/")
                except (IndexError, ValueError):
                    pass

    # Try strategies in preferred order (learned from previous successes)
    for strategy in _preferred_strategy:
        for url in candidates.get(strategy, []):
            if _probe_url(url):
                _url_cache[obsid] = url
                # Move this strategy to front for future calls
                if _preferred_strategy[0] != strategy:
                    _preferred_strategy.remove(strategy)
                    _preferred_strategy.insert(0, strategy)
                return url

    # Nothing worked — return None
    return None


def _probe_url(url: str) -> bool:
    """Return True if *url* exists (HTTP 200 on HEAD or GET)."""
    try:
        resp = requests.head(url, timeout=10, allow_redirects=True)
        if resp.status_code == 200:
            return True
        # Some servers don't support HEAD on directories — try GET
        if resp.status_code == 405:
            resp = requests.get(url, timeout=10, stream=True)
            resp.close()
            return resp.status_code == 200
        return False
    except requests.RequestException:
        return False


def get_remote_size(url: str) -> int | None:
    """Return Content-Length from a HEAD request, or None if unavailable."""
    try:
        resp = requests.head(url, timeout=15, allow_redirects=True)
        resp.raise_for_status()
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl else None
    except (requests.RequestException, ValueError):
        return None


def download_file(url: str, dest: Path, overwrite: bool = False,
                  chunk_size: int = 1 << 16) -> tuple[Path, bool]:
    """Stream-download *url* to *dest*, showing progress.

    Returns (dest, was_downloaded).  If the file already exists with the
    correct size and *overwrite* is False, the download is skipped.
    """
    # --- Resume / skip logic ---
    if dest.exists() and not overwrite:
        local_size = dest.stat().st_size
        remote_size = get_remote_size(url)
        if remote_size is not None and local_size == remote_size:
            print(f"  {dest.name}: already exists ({local_size/1e6:.2f} MB), skipping.")
            return dest, False
        elif remote_size is None:
            # Can't verify size – keep existing file to be safe
            print(f"  {dest.name}: exists (size unverifiable), skipping. Use --overwrite to force.")
            return dest, False
        else:
            print(f"  {dest.name}: local size ({local_size}) ≠ remote ({remote_size}), re-downloading.")

    # --- Download ---
    resp = requests.get(url, stream=True, timeout=120)
    resp.raise_for_status()
    total = int(resp.headers.get("Content-Length", 0))
    downloaded = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=chunk_size):
            fh.write(chunk)
            downloaded += len(chunk)
            if total:
                pct = downloaded / total * 100
                print(f"\r  {dest.name}: {downloaded/1e6:.1f}/{total/1e6:.1f} MB "
                      f"({pct:.0f}%)", end="", flush=True)
    print()
    return dest, True


def list_remote_dir(url: str) -> list[str]:
    """Scrape a simple Apache/nginx directory listing for links."""
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return re.findall(r'href="([^"?]+)"', resp.text)


def build_file_filter(modes: list[str] | None = None,
                      clean_only: bool = False):
    """Return a predicate that accepts/rejects filenames based on filters.

    Parameters
    ----------
    modes : list of str or None
        XRT mode codes to keep, e.g. ["pc", "wt"].  Matched against the
        standard Swift XRT filename convention where the mode appears as
        ``xpc`` (Photon Counting), ``xwt`` (Windowed Timing), or ``xim``
        (Image) right after the obsid.  If None, all modes pass.
    clean_only : bool
        If True, reject unfiltered event files (those ending in ``_uf.evt``
        or ``_uf.evt.gz``).

    Returns
    -------
    callable(str) -> bool
        True if the file should be downloaded, False to skip it.

    Note
    ----
    An auxil-payload filter (keeping only the files xrtpipeline reads) was
    considered, but xrtpipeline's ``prefilter`` step requires the orbit file
    ``SWIFT_TLE_ARCHIVE.txt`` and the spacecraft-housekeeping ``sen.hk`` —
    both of which the bug log had assumed were dispensable.  Pruning them
    breaks every OBSID, so the full auxil/ directory is downloaded as-is.
    """
    # Pre-compile patterns
    # Swift event filenames look like: sw{obsid}x{mode}{…}.evt[.gz]
    # mode portion: pc, wt, im  (after the leading 'x')
    mode_pattern = None
    if modes:
        codes = "|".join(m.lower() for m in modes)
        # Match the mode code that appears after 'x' in the XRT filename
        mode_pattern = re.compile(rf"x(?:{codes})", re.IGNORECASE)

    uf_pattern = re.compile(r"_uf\.evt(\.gz)?$", re.IGNORECASE) if clean_only else None

    def _accept(fname: str) -> bool:
        # Mode filter only applies to event files
        is_event = fname.endswith((".evt", ".evt.gz"))
        if mode_pattern and is_event:
            if not mode_pattern.search(fname):
                return False
        if uf_pattern and is_event:
            if uf_pattern.search(fname):
                return False
        return True

    return _accept


def download_xrt_products(obsid: str, outdir: Path,
                          products: list[str] | None = None,
                          overwrite: bool = False,
                          file_filter=None,
                          obs_meta: dict | None = None) -> dict:
    """
    Download XRT data products for a single observation.

    Parameters
    ----------
    obsid : str
        Swift observation ID.
    outdir : Path
        Local directory to save files.
    products : list of str or None
        Subdirectories to download from, e.g. ["xrt/event", "xrt/products"].
        If None, downloads the xrt/event directory (cleaned event files).
    overwrite : bool
        If True, re-download files even when they already exist with the
        correct size.
    file_filter : callable or None
        A predicate ``f(filename) -> bool``.  Files for which it returns
        False are skipped.  If None, all files are downloaded.
    obs_meta : dict or None
        Catalog metadata for this observation (used for URL discovery).

    Returns
    -------
    dict with keys: downloaded (int), skipped (int), failed (int), filtered (int)
    """
    if products is None:
        # NOTE: auxil/ lives at the OBSID top level (<obsid>/auxil/), NOT under
        # xrt/ — xrtpipeline needs the attitude file sw<obsid>sat.fits.gz, the
        # orbit file SWIFT_TLE_ARCHIVE.txt, and sw<obsid>sen.hk from here. A
        # previous "xrt/auxil" default 404'd and broke every OBSID.
        products = ["xrt/event", "xrt/hk", "auxil"]

    stats = {"downloaded": 0, "skipped": 0, "failed": 0, "filtered": 0}

    # --- Discover the correct archive URL for this obsid ---
    base = resolve_data_url(obsid, obs_meta=obs_meta)
    if base is None:
        print(f"[warn] Could not find {obsid} in any known archive location. Skipping.")
        stats["failed"] += 1
        return stats

    print(f"[info] Archive URL: {base}")
    obs_dir = outdir / obsid
    obs_dir.mkdir(parents=True, exist_ok=True)

    for prod in products:
        prod_url = base + prod + "/"
        print(f"[info] Listing {prod_url}")
        try:
            links = list_remote_dir(prod_url)
        except requests.HTTPError as exc:
            print(f"[warn] Could not list {prod_url} ({exc})")
            continue

        files = [l for l in links if not l.endswith("/") and l not in (".", "..")]
        if not files:
            print(f"  (no files found in {prod})")
            continue

        local_prod_dir = obs_dir / prod.replace("/", os.sep)
        local_prod_dir.mkdir(parents=True, exist_ok=True)

        for fname in files:
            # Apply mode / clean-only filter
            if file_filter and not file_filter(fname):
                print(f"  {fname}: filtered out, skipping.")
                stats["filtered"] += 1
                continue

            file_url = prod_url + fname
            dest = local_prod_dir / fname
            try:
                _, was_downloaded = download_file(file_url, dest, overwrite=overwrite)
                if was_downloaded:
                    stats["downloaded"] += 1
                else:
                    stats["skipped"] += 1
            except requests.HTTPError as exc:
                print(f"  [warn] Failed to download {fname}: {exc}")
                stats["failed"] += 1

    return stats


# ---------------------------------------------------------------------------
# Progress tracker with ETA
# ---------------------------------------------------------------------------

class ProgressTracker:
    """Tracks elapsed time per observation and provides ETA estimates."""

    def __init__(self, total: int):
        self.total = total
        self.completed = 0
        self.start_time = time.monotonic()
        self.obs_times: list[float] = []  # seconds per observation
        self._obs_start: float = 0.0

    def begin_obs(self):
        """Call at the start of each observation download."""
        self._obs_start = time.monotonic()

    def end_obs(self):
        """Call when an observation finishes. Updates stats."""
        elapsed = time.monotonic() - self._obs_start
        self.obs_times.append(elapsed)
        self.completed += 1

    def _fmt_duration(self, seconds: float) -> str:
        """Human-readable duration string."""
        if seconds < 60:
            return f"{seconds:.0f}s"
        elif seconds < 3600:
            m, s = divmod(int(seconds), 60)
            return f"{m}m {s:02d}s"
        else:
            h, rem = divmod(int(seconds), 3600)
            m, s = divmod(rem, 60)
            return f"{h}h {m:02d}m {s:02d}s"

    def summary_line(self) -> str:
        """Return a one-line progress + ETA string."""
        remaining = self.total - self.completed
        wall = time.monotonic() - self.start_time

        if self.obs_times:
            avg = sum(self.obs_times) / len(self.obs_times)
            eta = avg * remaining
            eta_str = self._fmt_duration(eta)
        else:
            eta_str = "calculating…"

        elapsed_str = self._fmt_duration(wall)
        return (f"[progress] {self.completed}/{self.total} observations done | "
                f"elapsed {elapsed_str} | "
                f"~{eta_str} remaining ({remaining} left)")

    def final_summary(self) -> str:
        """Summary line for when everything is finished."""
        wall = time.monotonic() - self.start_time
        return (f"[progress] {self.completed}/{self.total} observations completed in "
                f"{self._fmt_duration(wall)}")




def _get_obsid(obs: dict) -> str:
    """Extract the observation ID from a catalog row, trying common column names."""
    for key in ("obsid", "obs_id", "OBSID", "OBS_ID"):
        val = obs.get(key, "").strip()
        if val:
            return val
    return ""


def print_observations(observations: list[dict]):
    """Print a summary table of observations."""
    # Preferred display columns, with aliases (first match wins per group).
    col_groups = [
        ["obsid", "obs_id"],
        ["name", "target_name"],
        ["start_time"],
        ["xrt_exposure", "xrt_expo"],
        ["xrt_expo_pc"],
        ["xrt_expo_wt"],
    ]

    # Determine which columns are actually present
    available: list[str] = []
    sample_keys = set()
    for obs in observations:
        sample_keys.update(obs.keys())
    for group in col_groups:
        for alias in group:
            if alias in sample_keys:
                available.append(alias)
                break

    if not available:
        available = list(observations[0].keys())[:5]

    # Compute column widths
    col_widths = {c: max(len(c), 11) for c in available}
    for obs in observations:
        for c in available:
            col_widths[c] = max(col_widths[c], len(obs.get(c, "")))
    # Cap width to keep things readable
    for c in col_widths:
        col_widths[c] = min(col_widths[c], 30)

    header = " | ".join(f"{c:>{col_widths[c]}s}" for c in available)
    print(f"\n{'='*len(header)}")
    print(header)
    print(f"{'-'*len(header)}")
    for obs in observations:
        row = " | ".join(f"{obs.get(c,'')[:col_widths[c]]:>{col_widths[c]}s}"
                         for c in available)
        print(row)
    print(f"{'='*len(header)}\n")


# ---------------------------------------------------------------------------
# Obsid file parser
# ---------------------------------------------------------------------------

def parse_obsid_file(filepath: str) -> list[str]:
    """Read observation IDs from a text file (one per line).

    Blank lines and lines starting with '#' are ignored.  If a line contains
    whitespace-separated columns (e.g. "00035393001  Crab"), only the first
    token is taken as the obsid.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"Obsid file not found: {filepath}")

    obsids: list[str] = []
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            token = line.split()[0]
            # Basic sanity: obsids should be numeric (possibly with leading zeros)
            if not token.replace("-", "").isdigit():
                print(f"[warn] {filepath}:{lineno}: '{token}' doesn't look like "
                      f"an obsid, skipping.")
                continue
            obsids.append(token)

    if not obsids:
        raise ValueError(f"No valid observation IDs found in {filepath}")

    print(f"[info] Read {len(obsids)} observation ID(s) from {filepath}")
    return obsids


# ---------------------------------------------------------------------------
# Date-window filtering (client-side, applied to catalog query results)
# ---------------------------------------------------------------------------

# MJD 40587 == 1970-01-01 (the Unix epoch), used to turn an MJD into a date.
_MJD_UNIX_EPOCH = 40587
# Catalog columns that may carry an observation date, in preference order.
_DATE_FIELDS = ("start_time", "time", "start_date", "obs_start")


def _obs_to_date(obs: dict) -> date | None:
    """Best-effort UTC calendar date for one catalog row.

    swiftmastr returns ``start_time`` as an MJD float (e.g. '56725.7111');
    some other HEASARC tables return an ISO 'YYYY-MM-DD …' string.  Returns a
    ``datetime.date`` or None if no date field is parseable.
    """
    for field in _DATE_FIELDS:
        val = (obs.get(field) or "").strip()
        if not val:
            continue
        # ISO date string?
        m = re.match(r"(\d{4}-\d{2}-\d{2})", val)
        if m:
            try:
                return date.fromisoformat(m.group(1))
            except ValueError:
                pass
        # MJD float?
        try:
            mjd = float(val.split()[0])
        except ValueError:
            continue
        if 30000 <= mjd <= 90000:  # sane MJD range (~1958 .. ~2123)
            return date(1970, 1, 1) + timedelta(days=int(mjd - _MJD_UNIX_EPOCH))
    return None


def filter_by_date(observations: list[dict],
                   start_date: date | None,
                   end_date: date | None) -> tuple[list[dict], int]:
    """Keep observations with start date in [start_date, end_date).

    The end is exclusive (the conventional half-open time window).  Rows with
    no parseable date are dropped when a filter is active; the count of such
    rows is returned so the caller can warn about them.
    """
    if not start_date and not end_date:
        return observations, 0
    kept: list[dict] = []
    undated = 0
    for obs in observations:
        d = _obs_to_date(obs)
        if d is None:
            undated += 1
            continue
        if start_date and d < start_date:
            continue
        if end_date and d >= end_date:
            continue
        kept.append(obs)
    return kept, undated


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _isodate(s: str) -> date:
    """argparse type converter: parse YYYY-MM-DD, error clearly otherwise."""
    try:
        return date.fromisoformat(s)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"invalid date '{s}': expected YYYY-MM-DD")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Download Swift XRT data from the HEASARC archive.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--name", "-n", type=str,
                      help="Source name to resolve (e.g. 'Crab Nebula', "
                           "'GRS 1915+105'). Resolves to coordinates then "
                           "performs a 12' cone search; multiple targets in "
                           "the same field will all be returned. Filter on "
                           "observation date or target name yourself in the "
                           "listing output if you need only one source.")
    grp.add_argument("--ra", type=float,
                      help="Right Ascension in decimal degrees (J2000). Must also give --dec.")
    grp.add_argument("--obsid", type=str,
                      help="Download a single observation by its ID (e.g. 00035393001).")
    grp.add_argument("--obsid-file", type=str, metavar="FILE",
                      help="Path to a text file with one observation ID per line.")

    p.add_argument("--dec", type=float, default=None,
                   help="Declination in decimal degrees (J2000).")
    p.add_argument("--radius", "-r", type=float, default=12.0,
                   help="Search radius in arcminutes (default: 12).")
    p.add_argument("--max-obs", "-m", type=int, default=None,
                   help="Maximum number of observations to download.")
    p.add_argument("--test", "-t", type=int, default=None, metavar="N",
                   help="Test mode: download only the first N observations.")
    p.add_argument("--overwrite", action="store_true",
                   help="Force re-download of files even if they already exist "
                        "with the correct size.")
    p.add_argument("--outdir", "-o", type=str, default="./swift_xrt_data",
                   help="Output directory (default: ./swift_xrt_data).")
    p.add_argument("--list-only", "-l", action="store_true",
                   help="List matching observations without downloading.")
    p.add_argument("--start-date", type=_isodate, default=None,
                   metavar="YYYY-MM-DD",
                   help="Only keep observations on or after this date (UTC, "
                        "inclusive). Applies to --list-only and to the "
                        "--name / --ra / --dec catalog query.")
    p.add_argument("--end-date", type=_isodate, default=None,
                   metavar="YYYY-MM-DD",
                   help="Only keep observations strictly before this date "
                        "(UTC, exclusive end). Applies to --list-only and to "
                        "the --name / --ra / --dec catalog query.")
    p.add_argument("--products", type=str, nargs="+",
                   default=None,
                   help="XRT product subdirs to download, e.g. 'xrt/event xrt/products'. "
                        "Default: xrt/event xrt/hk auxil")
    p.add_argument("--mode", type=str, nargs="+", default=None,
                   metavar="MODE",
                   help="Only download event files for these XRT modes. "
                        "Accepted codes: pc (Photon Counting), wt (Windowed "
                        "Timing), im (Image). E.g. --mode pc wt. "
                        "Non-event files (housekeeping etc.) are unaffected.")
    p.add_argument("--clean-only", action="store_true",
                   help="Skip unfiltered (_uf) event files and download only "
                        "cleaned/screened event files.")
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    # --- Determine input mode and build obsid list ---
    # Mode A: direct obsid(s) – skip catalog query entirely
    # Mode B: name or coordinates – query HEASARC, then select obsids

    obsid_list: list[str] = []          # final list of obsids to download
    observations: list[dict] = []       # catalog rows (only for Mode B)

    # Date filtering needs the catalog metadata, which only the --name/--ra/--dec
    # path produces. For direct obsid input there is nothing to filter against.
    if (args.start_date or args.end_date) and (args.obsid or args.obsid_file):
        print("[warn] --start-date/--end-date are ignored for --obsid / "
              "--obsid-file (no catalog dates available); use --name / --ra "
              "/ --dec to filter by date.", file=sys.stderr)

    if args.obsid:
        # --- Single obsid ---
        obsid_list = [args.obsid.strip()]
        print(f"[info] Single obsid mode: {obsid_list[0]}")

    elif args.obsid_file:
        # --- Obsid file ---
        obsid_list = parse_obsid_file(args.obsid_file)

    else:
        # --- Name or RA/Dec → catalog query ---
        if args.name:
            ra, dec = resolve_name(args.name)
        else:
            if args.dec is None:
                parser.error("--dec is required when using --ra")
            ra, dec = args.ra, args.dec

        observations = query_swift_master(ra, dec, radius_arcmin=args.radius)
        if not observations:
            print("[info] No observations found. Try increasing --radius.")
            sys.exit(0)

        # Client-side date-window filter (MJD-aware) on the catalog results.
        if args.start_date or args.end_date:
            n_before = len(observations)
            observations, undated = filter_by_date(
                observations, args.start_date, args.end_date)
            print(f"Listed {n_before} obs ({len(observations)} after date filter)",
                  file=sys.stderr)
            if undated:
                print(f"[warn] {undated} observation(s) had no parseable date "
                      f"and were dropped by the date filter.", file=sys.stderr)
            if not observations:
                print("[info] No observations remain after date filtering.")
                sys.exit(0)

        print_observations(observations)

        if args.list_only:
            print(f"[info] {len(observations)} observation(s) listed. "
                  "Use without --list-only to download.")
            sys.exit(0)

        obsid_list = [_get_obsid(obs) for obs in observations]
        obsid_list = [o for o in obsid_list if o]

    if not obsid_list:
        print("[info] No observation IDs to download.")
        sys.exit(0)

    # --- Apply --test / --max-obs limits ---
    limit = args.test if args.test is not None else args.max_obs
    if limit is not None:
        obsid_list = obsid_list[:limit]

    if args.test is not None:
        print(f"[test] TEST MODE – downloading only the first {args.test} "
              f"observation(s).")

    # Build a lookup for pretty names (only available in catalog mode)
    obs_lookup: dict[str, dict] = {}
    for obs in observations:
        oid = _get_obsid(obs)
        if oid:
            obs_lookup[oid] = obs

    # --- Download with progress tracking ---
    outdir = Path(args.outdir)
    overwrite = args.overwrite
    file_filter = build_file_filter(modes=args.mode, clean_only=args.clean_only)
    tracker = ProgressTracker(total=len(obsid_list))

    total_stats = {"downloaded": 0, "skipped": 0, "failed": 0, "filtered": 0}

    print(f"[info] Will download {len(obsid_list)} observation(s) to {outdir.resolve()}")
    if args.mode:
        print(f"[info] Mode filter: keeping only {', '.join(m.upper() for m in args.mode)} "
              f"event files.")
    if args.clean_only:
        print("[info] Clean-only: unfiltered (_uf) event files will be skipped.")
    if overwrite:
        print("[info] Overwrite mode ON – existing files will be re-downloaded.")
    else:
        print("[info] Existing files with correct size will be skipped "
              "(use --overwrite to force).")
    print()

    try:
        for i, obsid in enumerate(obsid_list, 1):
            label = obs_lookup.get(obsid, {}).get("name", "")
            name_str = f" ({label})" if label else ""
            print(f"=== Observation {i}/{len(obsid_list)}: {obsid}{name_str} ===")

            tracker.begin_obs()
            obs_meta = obs_lookup.get(obsid)
            stats = download_xrt_products(obsid, outdir, products=args.products,
                                          overwrite=overwrite,
                                          file_filter=file_filter,
                                          obs_meta=obs_meta)
            tracker.end_obs()

            for k in total_stats:
                total_stats[k] += stats[k]

            print(tracker.summary_line())
            print()

            # Be polite to the server
            if i < len(obsid_list):
                time.sleep(1)

    except KeyboardInterrupt:
        print("\n\n[interrupted] Download interrupted by user.")
        print(f"[interrupted] Completed {tracker.completed}/{tracker.total} "
              f"observations before interruption.")
        print("[interrupted] Re-run the same command to resume – existing "
              "files will be skipped automatically.")

    # --- Final summary ---
    print()
    print(tracker.final_summary())
    filtered_msg = (f", {total_stats['filtered']} filtered out"
                    if total_stats['filtered'] else "")
    print(f"[summary] Files: {total_stats['downloaded']} downloaded, "
          f"{total_stats['skipped']} skipped, "
          f"{total_stats['failed']} failed{filtered_msg}")
    print(f"[done] Data saved to {outdir.resolve()}")


if __name__ == "__main__":
    main()

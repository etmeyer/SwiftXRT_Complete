#!/usr/bin/env python3
"""
swift_xrt_extract_spectra.py

Automated spectral extraction pipeline for Swift XRT data,
supporting both PC and WT modes.

PC mode: Uses pile-up radii from swift_xrt_king_profile.py and
    pc_master_table.txt. Source regions are annuli (if piled up)
    or circles. Background from a large annulus far from source.

WT mode: Uses extraction parameters from swift_wt_summary_viewer.py
    and wt_master_table.txt. Source region is a circle (default 20
    pix). Background is an annulus symmetric about the WT window
    center (~100 pix). BACKSCAL keywords are corrected for WT 1D
    geometry per the UK SSDC standard.

Usage:
    # Both modes
    python swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524
    # PC only
    python swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524 --mode pc
    # WT only
    python swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524 --mode wt

Environment:
    HEASoft must be initialized (HEADAS set, headas-init.sh sourced)
    CALDB must be initialized (CALDB set, caldbinit.sh sourced)

Required input files:
    pc_master_table.txt      - PC observation table (for --mode pc/both)
    wt_master_table.txt      - WT observation table (for --mode wt/both)
    *_pileup.txt             - PC pile-up analysis results
    *_wt_profile.txt         - WT profile analysis results (optional)
"""

import os
import sys
import re
import glob
import subprocess
import argparse
import shutil
import numpy as np

try:
    from astropy.io import fits
except ImportError:
    print("ERROR: astropy is required.")
    sys.exit(1)


# ---------------------------------------------------------------
# Use the current working directory as the base directory.
# All OBSID folders must be subdirectories of this location.
# We resolve it to an absolute path so that relative references
# work correctly even if we cd into OBSID subdirectories.
# ---------------------------------------------------------------
BASE_DIR = os.path.abspath(os.getcwd())


# ---------------------------------------------------------------
# Environment checks
# ---------------------------------------------------------------

def check_environment():
    """
    Verify that HEASoft and CALDB are properly configured.
    We check for the environment variables and the availability
    of key executables.
    """
    headas = os.environ.get('HEADAS')
    caldb = os.environ.get('CALDB')

    if not headas:
        print("ERROR: HEADAS environment variable not set.")
        print("  Source your HEASoft init script first.")
        sys.exit(1)

    if not caldb:
        print("ERROR: CALDB environment variable not set.")
        print("  Source your CALDB init script first.")
        sys.exit(1)

    # Check that key tools are available
    for tool in ['xselect', 'xrtexpomap', 'xrtmkarf', 'grppha']:
        if shutil.which(tool) is None:
            print(f"ERROR: '{tool}' not found in PATH.")
            print("  Make sure HEASoft is fully installed and "
                  "initialized.")
            sys.exit(1)

    print(f"HEADAS: {headas}")
    print(f"CALDB:  {caldb}")
    return caldb


# ---------------------------------------------------------------
# Parse the master table
# ---------------------------------------------------------------

def read_master_table(table_file):
    """
    Read pc_master_table.txt and return list of dicts for
    observations marked as include=yes.

    Expected format (whitespace-separated, comment field in quotes):
        OBSID  filename  include  badstripe  "comment"
    """
    entries = []
    with open(table_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Skip header line
            if line.startswith('OBSID'):
                continue

            # Parse: handle quoted comment field
            # Split on quotes first to isolate comment
            quote_match = re.search(r'"([^"]*)"', line)
            comment = quote_match.group(1) if quote_match else ""

            # Remove the quoted part for parsing other fields
            before_quote = line[:quote_match.start()].strip() \
                if quote_match else line
            parts = before_quote.split()

            if len(parts) < 4:
                print(f"  WARNING: skipping malformed line {line_num}")
                continue

            entry = {
                'obsid': parts[0],
                'filename': parts[1],
                'include': parts[2].lower(),
                'badstripe': parts[3].lower(),
                'comment': comment,
            }
            entries.append(entry)

    included = [e for e in entries if e['include'] == 'yes']
    excluded = [e for e in entries if e['include'] != 'yes']
    print(f"Master table: {len(entries)} entries, "
          f"{len(included)} included, {len(excluded)} excluded.")

    return included


# ---------------------------------------------------------------
# Read pile-up information
# ---------------------------------------------------------------

def read_pileup_info(pileup_file):
    """
    Parse the _pileup.txt file produced by swift_xrt_king_profile.py.
    Returns dict with centroid position, plate scale, radii.
    """
    info = {
        'xc': None, 'yc': None,
        'ra': None, 'dec': None,
        'plate_scale': None,
        'pileup_radius': None,
        'override_radius': None,
        'count_rate': None,
    }

    with open(pileup_file, 'r') as f:
        for line in f:
            line = line.strip()

            m = re.search(
                r'Centroid position \(pix\):\s*X=([\d.]+)\s+Y=([\d.]+)',
                line)
            if m:
                info['xc'] = float(m.group(1))
                info['yc'] = float(m.group(2))

            m = re.search(
                r'Source position \(input\):\s*RA=([\d.]+)\s+'
                r'Dec=([\d.-]+)', line)
            if m:
                info['ra'] = float(m.group(1))
                info['dec'] = float(m.group(2))

            m = re.search(
                r'Plate scale:\s*([\d.]+)\s*arcsec/pixel', line)
            if m:
                info['plate_scale'] = float(m.group(1))

            m = re.search(
                r'Count rate:\s*([\d.]+)\s*ct/s', line)
            if m:
                info['count_rate'] = float(m.group(1))

            m = re.match(
                r'pileup_radius_arcsec\s*=\s*([\d.]+)', line)
            if m:
                info['pileup_radius'] = float(m.group(1))

            m = re.match(
                r'override_radius_arcsec\s*=\s*([\d.]+)', line)
            if m:
                info['override_radius'] = float(m.group(1))

    return info


# ---------------------------------------------------------------
# Determine inner exclusion radius
# ---------------------------------------------------------------

def get_inner_radius(pileup_info, default_inner=0.0):
    """
    Determine the inner radius for the source extraction annulus.

    Priority:
      1. Override radius (user-specified via pileup_overrides.txt)
      2. Auto-detected pile-up radius
      3. Zero (no pile-up exclusion)

    If the count rate is below ~0.5 ct/s, pile-up is negligible
    for PC mode and we use 0 regardless.
    """
    # If count rate is low enough, no pile-up correction needed
    rate = pileup_info.get('count_rate', 0)
    if rate is not None and rate < 0.5:
        return 0.0

    # Override takes priority
    if pileup_info.get('override_radius') is not None:
        return pileup_info['override_radius']

    # Then auto-detected
    if pileup_info.get('pileup_radius') is not None:
        return pileup_info['pileup_radius']

    return default_inner


# ---------------------------------------------------------------
# Determine optimal outer radius
# ---------------------------------------------------------------

def estimate_optimal_outer_radius(evt_file, xc, yc, plate_scale,
                                   r_inner_arcsec, r_max=80.0,
                                   rbin=2.0):
    """
    Estimate the optimal outer radius by finding where the
    cumulative signal-to-noise ratio peaks.

    As you increase the outer radius, you collect more source
    photons (signal), but the background area also grows. The
    S/N is approximately:

        S/N(r) = N_src(r) / sqrt(N_src(r) + N_bkg(r) * (A_src/A_bkg))

    In practice, for a point source on a low background, the S/N
    rises quickly as you capture more PSF flux, then flattens as
    you're mostly adding background. We look for the radius where
    adding another bin increases S/N by less than 1%.

    Returns:
        optimal outer radius in arcsec
    """
    with fits.open(evt_file) as hdul:
        evt_ext = None
        for i, ext in enumerate(hdul):
            if ext.name == 'EVENTS':
                evt_ext = i
                break
        if evt_ext is None:
            evt_ext = 1
        data = hdul[evt_ext].data
        x_events = data['X'].astype(float)
        y_events = data['Y'].astype(float)

    dist = np.sqrt((x_events - xc)**2 + (y_events - yc)**2)
    dist_arcsec = dist * plate_scale

    # Compute cumulative counts in expanding annuli starting
    # from r_inner
    edges = np.arange(r_inner_arcsec, r_max + rbin, rbin)
    if len(edges) < 3:
        return r_max

    cum_counts = np.zeros(len(edges) - 1)
    for i in range(len(edges) - 1):
        cum_counts[i] = np.sum(dist_arcsec < edges[i + 1])
        # Subtract the inner excluded region
        cum_counts[i] -= np.sum(dist_arcsec < r_inner_arcsec)

    # Estimate background from the outer portion (60-80 arcsec)
    bkg_mask = (dist_arcsec >= 60) & (dist_arcsec < 80)
    bkg_area = np.pi * (80**2 - 60**2)  # arcsec^2
    bkg_counts = np.sum(bkg_mask)
    bkg_rate_per_arcsec2 = bkg_counts / bkg_area if bkg_area > 0 else 0

    # Compute S/N for each outer radius
    r_outer = edges[1:]
    src_area = np.pi * (r_outer**2 - r_inner_arcsec**2)
    expected_bkg = bkg_rate_per_arcsec2 * src_area
    snr = np.where(cum_counts > 0,
                   cum_counts / np.sqrt(cum_counts + expected_bkg),
                   0)

    # Find where S/N flattens: marginal gain < 1%
    if np.max(snr) == 0:
        return 47.0  # safe default

    best_idx = 0
    for i in range(1, len(snr)):
        if snr[i] > snr[best_idx]:
            best_idx = i
        # Check if we've passed the peak or marginal gain < 1%
        if i > 2 and snr[i] < snr[i - 1]:
            best_idx = i - 1
            break
        if i > 2 and snr[best_idx] > 0:
            marginal = (snr[i] - snr[i - 1]) / snr[best_idx]
            if marginal < 0.01:
                best_idx = i
                break

    optimal_r = r_outer[best_idx]

    # Clamp to reasonable range
    optimal_r = max(optimal_r, 20.0)
    optimal_r = min(optimal_r, 70.0)

    return round(optimal_r, 1)


# ---------------------------------------------------------------
# Find a source-free background region
# ---------------------------------------------------------------

def find_background_region(evt_file, xc, yc, plate_scale,
                           r_outer_src, bkg_inner=100.0,
                           bkg_outer=160.0):
    """
    Determine a background region that avoids the source.

    Strategy: Use a large annulus centered on the source, far
    enough away that the PSF contribution is negligible. The
    default range of 100-160 arcsec is well outside the XRT PSF
    (which is effectively zero beyond ~60-70 arcsec) but still
    on the detector.

    We verify that the region actually contains events (i.e.,
    is on the detector). If not, we shrink the outer radius
    until it does.

    For observations where the source is near the edge, a
    circular offset region might be needed instead. We fall
    back to this if the annulus doesn't work.

    Returns:
        dict with 'type' ('annulus' or 'circle'),
        'params' (region-specific), and 'counts' (background cts)
    """
    with fits.open(evt_file) as hdul:
        evt_ext = None
        for i, ext in enumerate(hdul):
            if ext.name == 'EVENTS':
                evt_ext = i
                break
        if evt_ext is None:
            evt_ext = 1
        data = hdul[evt_ext].data
        x_events = data['X'].astype(float)
        y_events = data['Y'].astype(float)

    dist = np.sqrt((x_events - xc)**2 + (y_events - yc)**2)
    dist_arcsec = dist * plate_scale

    # Try annulus, shrinking outer radius if needed
    for bkg_out in [bkg_outer, 140.0, 120.0, 100.0]:
        if bkg_out <= bkg_inner:
            continue
        mask = (dist_arcsec >= bkg_inner) & (dist_arcsec < bkg_out)
        bkg_counts = np.sum(mask)
        if bkg_counts >= 5:  # need some counts to be usable
            r_inner_pix = bkg_inner / plate_scale
            r_outer_pix = bkg_out / plate_scale
            return {
                'type': 'annulus',
                'r_inner_arcsec': bkg_inner,
                'r_outer_arcsec': bkg_out,
                'r_inner_pix': r_inner_pix,
                'r_outer_pix': r_outer_pix,
                'counts': bkg_counts,
                'xc': xc, 'yc': yc,
            }

    # Fallback: offset circle to the side with lowest source
    # contamination. Place it at 120 arcsec from source with
    # radius 60 arcsec.
    offset_r = 120.0 / plate_scale
    bkg_r = 60.0 / plate_scale
    # Try four cardinal directions
    for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        bkg_xc = xc + dx * offset_r
        bkg_yc = yc + dy * offset_r
        d = np.sqrt((x_events - bkg_xc)**2 + (y_events - bkg_yc)**2)
        mask = d * plate_scale < 60.0
        if np.sum(mask) >= 5:
            return {
                'type': 'circle',
                'xc': bkg_xc, 'yc': bkg_yc,
                'r_pix': bkg_r,
                'r_arcsec': 60.0,
                'counts': int(np.sum(mask)),
            }

    # Last resort: use whatever annulus we can
    print("    WARNING: could not find good background region.")
    return {
        'type': 'annulus',
        'r_inner_arcsec': bkg_inner,
        'r_outer_arcsec': bkg_outer,
        'r_inner_pix': bkg_inner / plate_scale,
        'r_outer_pix': bkg_outer / plate_scale,
        'counts': 0,
        'xc': xc, 'yc': yc,
    }


# ---------------------------------------------------------------
# Write DS9 region files
# ---------------------------------------------------------------

def write_source_region(filepath, xc, yc, r_inner_pix, r_outer_pix):
    """
    Write a DS9-format source region file.
    If r_inner_pix > 0, write an annulus (for pile-up exclusion).
    Otherwise write a simple circle.
    """
    with open(filepath, 'w') as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=green dashlist=8 3 width=1 "
                "font=\"helvetica 10 normal roman\" select=1 "
                "highlite=1 dash=0 fixed=0 edit=1 move=1 "
                "delete=1 include=1 source=1\n")
        f.write("physical\n")
        if r_inner_pix > 0:
            f.write(f"annulus({xc:.4f},{yc:.4f},"
                    f"{r_inner_pix:.4f},{r_outer_pix:.4f})\n")
        else:
            f.write(f"circle({xc:.4f},{yc:.4f},"
                    f"{r_outer_pix:.4f})\n")


def write_background_region(filepath, bkg_info):
    """
    Write a DS9-format background region file.
    """
    with open(filepath, 'w') as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=cyan dashlist=8 3 width=1 "
                "font=\"helvetica 10 normal roman\" select=1 "
                "highlite=1 dash=0 fixed=0 edit=1 move=1 "
                "delete=1 include=1 source=1\n")
        f.write("physical\n")
        if bkg_info['type'] == 'annulus':
            f.write(f"annulus({bkg_info['xc']:.4f},"
                    f"{bkg_info['yc']:.4f},"
                    f"{bkg_info['r_inner_pix']:.4f},"
                    f"{bkg_info['r_outer_pix']:.4f})\n")
        else:
            f.write(f"circle({bkg_info['xc']:.4f},"
                    f"{bkg_info['yc']:.4f},"
                    f"{bkg_info['r_pix']:.4f})\n")


# ---------------------------------------------------------------
# Run xselect for spectrum extraction
# ---------------------------------------------------------------

def extract_spectrum_xselect(evt_file, region_file, output_pha,
                              session_name="xsel_tmp"):
    """
    Use xselect to extract a spectrum from an event file
    within a given region.

    xselect is driven via stdin commands piped to the process.
    """
    evt_dir = os.path.dirname(os.path.abspath(evt_file))
    evt_name = os.path.basename(evt_file)
    output_pha_abs = os.path.abspath(output_pha)
    region_abs = os.path.abspath(region_file)

    # xselect commands
    # Note: we use 'no' for saved session, set the data directory
    # and file, filter by region, extract spectrum, and save.
    commands = f"""{session_name}
no
read event
{evt_dir}
{evt_name}
yes
filter region {region_abs}
extract spectrum
save spectrum {output_pha_abs}
exit
no
"""

    # Remove any existing output file (xselect won't overwrite)
    if os.path.exists(output_pha_abs):
        os.remove(output_pha_abs)

    print(f"    [CMD] xselect (piped commands):")
    for line in commands.strip().split('\n'):
        print(f"           {line}")

    result = subprocess.run(
        ['xselect'],
        input=commands, capture_output=True, text=True,
        timeout=120)

    if not os.path.exists(output_pha_abs):
        print(f"    ERROR: xselect did not produce {output_pha}")
        print(f"    STDOUT (last 1000 chars):")
        for line in result.stdout[-1000:].split('\n'):
            print(f"      {line}")
        print(f"    STDERR (last 500 chars):")
        for line in result.stderr[-500:].split('\n'):
            print(f"      {line}")
        return False

    return True


# ---------------------------------------------------------------
# Run xrtexpomap
# ---------------------------------------------------------------

def find_auxiliary_files(obsid_path, obsid):
    """
    Locate the attitude file and housekeeping file needed by
    xrtexpomap and xrtmkarf. These are typically in the
    original data directory structure under auxil/ and xrt/hk/.

    We search recursively since directory structures vary.
    """
    # Attitude file: *pat.fits* or *sat.fits*
    att_files = glob.glob(os.path.join(
        obsid_path, '**', f'sw{obsid}*pat.fits*'), recursive=True)
    if not att_files:
        att_files = glob.glob(os.path.join(
            obsid_path, '**', f'sw{obsid}*sat.fits*'), recursive=True)
    att_file = att_files[0] if att_files else None

    # Housekeeping file: *xhd.hk*
    hk_files = glob.glob(os.path.join(
        obsid_path, '**', f'sw{obsid}*xhd.hk*'), recursive=True)
    hk_file = hk_files[0] if hk_files else None

    return att_file, hk_file


def run_xrtexpomap(evt_file, att_file, hk_file, out_dir):
    """
    Generate an exposure map using xrtexpomap.
    Returns the path to the exposure map, or None on failure.

    The exposure map accounts for bad pixels/columns,
    vignetting, and CCD window to produce a map of effective
    exposure time per pixel.
    """
    # xrtexpomap names the output based on the input filename
    evt_base = os.path.basename(evt_file).replace('.gz', '')
    expo_name = evt_base.replace('_cl.evt', '_ex.img')
    expo_path = os.path.join(out_dir, expo_name)

    # Skip if already exists
    if os.path.exists(expo_path):
        print(f"    Exposure map exists: {expo_name}")
        return expo_path

    cmd = [
        'xrtexpomap',
        f'infile={os.path.abspath(evt_file)}',
        f'attfile={os.path.abspath(att_file)}',
        f'hdfile={os.path.abspath(hk_file)}',
        f'outdir={os.path.abspath(out_dir)}',
        'clobber=yes',
    ]

    print(f"    [CMD] {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=300)

    if os.path.exists(expo_path):
        return expo_path
    else:
        print(f"    WARNING: xrtexpomap failed.")
        print(f"    STDOUT (last 500 chars):")
        for line in result.stdout[-500:].split('\n'):
            print(f"      {line}")
        print(f"    STDERR (last 500 chars):")
        for line in result.stderr[-500:].split('\n'):
            print(f"      {line}")
        return None


# ---------------------------------------------------------------
# Run xrtmkarf
# ---------------------------------------------------------------

def run_xrtmkarf(src_pha, expo_file, out_arf, src_region):
    """
    Generate an ARF using xrtmkarf.

    We set rmffile=CALDB so it automatically selects the correct
    RMF for the observation date/mode, and enable PSF correction
    for point sources (psfflag=yes).

    xrtmkarf also reports which RMF it used, which we capture
    from its output for later use in grppha.

    Returns:
        (arf_path, rmf_path) on success, (None, None) on failure.
    """
    src_pha_abs = os.path.abspath(src_pha)
    expo_abs = os.path.abspath(expo_file)
    out_arf_abs = os.path.abspath(out_arf)

    cmd = [
        'xrtmkarf',
        f'phafile={src_pha_abs}',
        f'expofile={expo_abs}',
        f'outfile={out_arf_abs}',
        'srcx=-1',               # auto-detect from WMAP in PHA
        'srcy=-1',
        'psfflag=yes',           # PSF correction for point source
        'rmffile=CALDB',         # auto-select from CALDB
        'extended=no',
        'clobber=yes',
    ]

    print(f"    [CMD] {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=300)

    # Parse output to find which RMF was used
    rmf_path = None
    for line in result.stdout.split('\n'):
        if 'Processing' in line and 'rmf' in line.lower():
            # Extract the RMF filename from the CALDB path
            m = re.search(r"'([^']*rmf[^']*)'", line)
            if m:
                rmf_full = m.group(1)
                # Convert remote URL or absolute path to $CALDB
                # relative path for portability
                if 'caldb/data/' in rmf_full:
                    idx = rmf_full.index('caldb/data/')
                    rmf_path = '$CALDB/data/' + \
                        rmf_full[idx + len('caldb/data/'):]
                elif rmf_full.startswith('/'):
                    # Local absolute path — try to make $CALDB relative
                    caldb = os.environ.get('CALDB', '')
                    if caldb and rmf_full.startswith(caldb):
                        rmf_path = '$CALDB' + \
                            rmf_full[len(caldb):]
                    else:
                        rmf_path = rmf_full
                else:
                    rmf_path = rmf_full

    if not os.path.exists(out_arf_abs):
        print(f"    WARNING: xrtmkarf failed.")
        print(f"    STDOUT (last 1000 chars):")
        for line in result.stdout[-1000:].split('\n'):
            print(f"      {line}")
        print(f"    STDERR (last 500 chars):")
        for line in result.stderr[-500:].split('\n'):
            print(f"      {line}")
        return None, None

    if rmf_path is None:
        print("    WARNING: could not determine RMF from xrtmkarf "
              "output. Full stdout:")
        for line in result.stdout.split('\n'):
            print(f"      {line}")

    print(f"    RMF: {rmf_path}")
    return out_arf_abs, rmf_path


# ---------------------------------------------------------------
# Run grppha for grouping and linking response files
# ---------------------------------------------------------------

def run_grppha(src_pha, out_pha, bkg_pha, arf_file, rmf_path,
               min_counts=20):
    """
    Use grppha to:
      1. Associate the background PHA, ARF, and RMF with the
         source spectrum (using $CALDB-relative path for RMF).
      2. Mark channels 0-29 as bad (below ~0.3 keV, unreliable).
      3. Group to a minimum number of counts per bin.

    The grouping threshold controls the statistics:
      - min_counts=20: suitable for chi-squared fitting
      - min_counts=1:  use with Cash statistics (cstat) in XSPEC

    Returns True on success.
    """
    src_abs = os.path.abspath(src_pha)
    out_abs = os.path.abspath(out_pha)
    bkg_abs = os.path.abspath(bkg_pha)
    arf_abs = os.path.abspath(arf_file)

    # Build grppha commands
    # grppha prompts: input file, output file, then GRPPHA[] commands
    # Note: RMF path uses $CALDB for portability
    commands = (
        f"{src_abs}\n"
        f"{out_abs}\n"
        f"bad 0-29\n"
        f"chkey backfile {bkg_abs}\n"
        f"chkey ancrfile {arf_abs}\n"
        f"chkey respfile {rmf_path}\n"
        f"group min {min_counts}\n"
        f"exit\n"
    )

    if os.path.exists(out_abs):
        os.remove(out_abs)

    print(f"    [CMD] grppha (piped commands):")
    for line in commands.strip().split('\n'):
        print(f"           {line}")

    result = subprocess.run(
        ['grppha'],
        input=commands, capture_output=True, text=True,
        timeout=60)

    if os.path.exists(out_abs):
        return True
    else:
        print(f"    WARNING: grppha failed.")
        print(f"    STDOUT (last 1000 chars):")
        for line in result.stdout[-1000:].split('\n'):
            print(f"      {line}")
        print(f"    STDERR (last 500 chars):")
        for line in result.stderr[-500:].split('\n'):
            print(f"      {line}")
        return False


# ---------------------------------------------------------------
# Get spectrum summary information
# ---------------------------------------------------------------

def get_spectrum_info(pha_file):
    """
    Read key information from a PHA spectrum file:
    total counts, exposure, channel range, grouping info.
    """
    info = {}
    with fits.open(pha_file) as hdul:
        header = hdul[1].header
        data = hdul[1].data

        info['exposure'] = header.get('EXPOSURE', 0)
        info['backscal'] = header.get('BACKSCAL', 0)
        info['total_counts'] = int(np.sum(data['COUNTS']))
        info['n_channels'] = len(data['COUNTS'])

        # Count grouped bins if GROUPING column exists
        if 'GROUPING' in data.columns.names:
            # In OGIP standard, GROUPING=1 marks start of a new bin
            info['n_grouped_bins'] = int(np.sum(data['GROUPING'] == 1))
        else:
            info['n_grouped_bins'] = None

    return info


# ---------------------------------------------------------------
# Read WT profile information
# ---------------------------------------------------------------

def read_wt_profile_info(profile_file):
    """
    Parse the _wt_profile.txt file produced by
    swift_wt_summary_viewer.py.

    Returns dict with source position, extraction radii,
    and BACKSCAL values.
    """
    info = {
        'xc': None, 'yc': None,
        'src_radius': 20,
        'bkg_inner': 80, 'bkg_outer': 120,
        'backscal_src': None, 'backscal_bkg': None,
        'plate_scale': 2.36,
    }

    with open(profile_file, 'r') as f:
        for line in f:
            line = line.strip()
            m = re.match(r'source_x\s*=\s*([\d.]+)', line)
            if m:
                info['xc'] = float(m.group(1))
            m = re.match(r'source_y\s*=\s*([\d.]+)', line)
            if m:
                info['yc'] = float(m.group(1))
            m = re.match(r'source_radius_pix\s*=\s*(\d+)', line)
            if m:
                info['src_radius'] = int(m.group(1))
            m = re.match(r'bkg_inner_pix\s*=\s*(\d+)', line)
            if m:
                info['bkg_inner'] = int(m.group(1))
            m = re.match(r'bkg_outer_pix\s*=\s*(\d+)', line)
            if m:
                info['bkg_outer'] = int(m.group(1))
            m = re.match(r'backscal_src\s*=\s*(\d+)', line)
            if m:
                info['backscal_src'] = int(m.group(1))
            m = re.match(r'backscal_bkg\s*=\s*(\d+)', line)
            if m:
                info['backscal_bkg'] = int(m.group(1))
            m = re.search(r'Plate scale:\s*([\d.]+)', line)
            if m:
                info['plate_scale'] = float(m.group(1))

    return info


# ---------------------------------------------------------------
# Fix WT BACKSCAL keywords
# ---------------------------------------------------------------

def fix_wt_backscal(pha_file, backscal_value):
    """
    Correct the BACKSCAL keyword in a WT-mode PHA file.

    XSELECT sets BACKSCAL based on the 2D area of the extraction
    region, but for WT mode it should reflect the 1D extent in
    the DETX direction. See:
    https://www.swift.ac.uk/analysis/xrt/backscal.php

    For a circular source region of radius r:
        BACKSCAL = 2 * r  (diameter in pixels)
    For an annular background region (r_inner to r_outer):
        BACKSCAL = r_outer - r_inner - 1
    """
    with fits.open(pha_file, mode='update') as hdul:
        old_val = hdul[1].header.get('BACKSCAL', 'N/A')
        hdul[1].header['BACKSCAL'] = backscal_value
        hdul[1].header.add_comment(
            f'BACKSCAL corrected for WT 1D geometry '
            f'(was {old_val})')


# ---------------------------------------------------------------
# Write WT region files
# ---------------------------------------------------------------

def write_wt_source_region(filepath, xc, yc, radius):
    """Write a circular source region for WT mode."""
    with open(filepath, 'w') as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=green dashlist=8 3 width=1 "
                "font=\"helvetica 10 normal roman\" select=1 "
                "highlite=1 dash=0 fixed=0 edit=1 move=1 "
                "delete=1 include=1 source=1\n")
        f.write("physical\n")
        f.write(f"circle({xc:.4f},{yc:.4f},{radius:.4f})\n")


def write_wt_background_region(filepath, xc, yc,
                                bkg_inner, bkg_outer):
    """Write an annular background region for WT mode."""
    with open(filepath, 'w') as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=cyan dashlist=8 3 width=1 "
                "font=\"helvetica 10 normal roman\" select=1 "
                "highlite=1 dash=0 fixed=0 edit=1 move=1 "
                "delete=1 include=1 source=1\n")
        f.write("physical\n")
        f.write(f"annulus({xc:.4f},{yc:.4f},"
                f"{bkg_inner:.4f},{bkg_outer:.4f})\n")


# ---------------------------------------------------------------
# Process one PC-mode observation
# ---------------------------------------------------------------

def process_pc_observation(entry, ra_src, dec_src, r_outer,
                           min_counts, bkg_inner, bkg_outer):
    """
    Full extraction pipeline for a single PC-mode observation.
    """
    obsid = entry['obsid']
    stem = entry['filename']
    obsid_path = os.path.join(BASE_DIR, obsid)

    print(f"\n{'='*65}")
    print(f"  OBSID: {obsid}    File: {stem}")
    print(f"{'='*65}")

    # ---- Locate event file ----
    evt_file = None
    for ext in ['_cl.evt', '_cl.evt.gz']:
        candidate = os.path.join(obsid_path, stem + ext)
        if os.path.exists(candidate):
            evt_file = candidate
            break
    # Search recursively if not in top level
    if evt_file is None:
        candidates = glob.glob(os.path.join(
            obsid_path, '**', stem + '_cl.evt*'), recursive=True)
        if candidates:
            evt_file = candidates[0]

    if evt_file is None:
        print(f"  ERROR: event file not found for {stem}")
        return None

    # ---- Read observation date from header ----
    with fits.open(evt_file) as hdul:
        evt_header = hdul[1].header
        date_obs = evt_header.get('DATE-OBS', 'N/A')
    print(f"  Date: {date_obs}")

    # ---- Read pile-up information ----
    pileup_file = os.path.join(obsid_path, f'{stem}_pileup.txt')
    if not os.path.exists(pileup_file):
        print(f"  ERROR: {stem}_pileup.txt not found. "
              f"Run swift_xrt_king_profile.py first.")
        return None

    pileup_info = read_pileup_info(pileup_file)
    if pileup_info['xc'] is None:
        print(f"  ERROR: no centroid in pileup file.")
        return None

    xc = pileup_info['xc']
    yc = pileup_info['yc']
    plate_scale = pileup_info.get('plate_scale', 2.36)

    # ---- Determine extraction radii ----
    r_inner = get_inner_radius(pileup_info)
    r_inner_pix = r_inner / plate_scale

    if r_outer == 'auto':
        r_out = estimate_optimal_outer_radius(
            evt_file, xc, yc, plate_scale, r_inner)
        print(f"  Auto outer radius: {r_out:.1f}\"")
    else:
        r_out = float(r_outer)

    r_outer_pix = r_out / plate_scale

    if r_inner > 0:
        print(f"  Source region: annulus {r_inner:.1f}\" - "
              f"{r_out:.1f}\" ({r_inner_pix:.1f} - "
              f"{r_outer_pix:.1f} pix)")
    else:
        print(f"  Source region: circle {r_out:.1f}\" "
              f"({r_outer_pix:.1f} pix)")

    # ---- Write source region file ----
    src_reg = os.path.join(obsid_path, f'{stem}_src.reg')
    write_source_region(src_reg, xc, yc, r_inner_pix, r_outer_pix)
    print(f"  Wrote: {os.path.basename(src_reg)}")

    # ---- Find background region ----
    bkg_info = find_background_region(
        evt_file, xc, yc, plate_scale, r_out,
        bkg_inner=bkg_inner, bkg_outer=bkg_outer)

    bkg_reg = os.path.join(obsid_path, f'{stem}_bkg.reg')
    write_background_region(bkg_reg, bkg_info)
    if bkg_info['type'] == 'annulus':
        print(f"  Background: annulus "
              f"{bkg_info['r_inner_arcsec']:.0f}\" - "
              f"{bkg_info['r_outer_arcsec']:.0f}\" "
              f"({bkg_info['counts']} cts)")
    else:
        print(f"  Background: offset circle "
              f"r={bkg_info['r_arcsec']:.0f}\" "
              f"({bkg_info['counts']} cts)")

    # ---- Extract source spectrum ----
    src_pha = os.path.join(obsid_path, f'{stem}_src.pha')
    print(f"  Extracting source spectrum...")
    if not extract_spectrum_xselect(evt_file, src_reg, src_pha,
                                     session_name="src"):
        print(f"  ERROR: source spectrum extraction failed.")
        return None

    # ---- Extract background spectrum ----
    bkg_pha = os.path.join(obsid_path, f'{stem}_bkg.pha')
    print(f"  Extracting background spectrum...")
    if not extract_spectrum_xselect(evt_file, bkg_reg, bkg_pha,
                                     session_name="bkg"):
        print(f"  ERROR: background spectrum extraction failed.")
        return None

    # ---- Generate exposure map ----
    att_file, hk_file = find_auxiliary_files(obsid_path, obsid)
    expo_file = None
    if att_file and hk_file:
        print(f"  Generating exposure map...")
        expo_file = run_xrtexpomap(evt_file, att_file, hk_file,
                                    obsid_path)
    else:
        # Check if pipeline already made one
        expo_candidates = glob.glob(os.path.join(
            obsid_path, '**', '*xpc*_ex.img*'), recursive=True)
        if expo_candidates:
            expo_file = expo_candidates[0]
            print(f"  Using existing exposure map: "
                  f"{os.path.basename(expo_file)}")
        else:
            print(f"  WARNING: cannot generate exposure map "
                  f"(missing att/hk files).")

    # ---- Generate ARF ----
    arf_file = os.path.join(obsid_path, f'{stem}.arf')
    rmf_path = None
    if expo_file:
        print(f"  Generating ARF...")
        arf_file, rmf_path = run_xrtmkarf(
            src_pha, expo_file, arf_file, src_reg)
    else:
        print(f"  WARNING: skipping ARF (no exposure map).")
        arf_file = None

    if arf_file is None or rmf_path is None:
        print(f"  ERROR: ARF/RMF generation failed.")
        return None

    # ---- Group spectrum ----
    grp_pha = os.path.join(obsid_path, f'{stem}_grp.pha')
    print(f"  Grouping spectrum (min {min_counts} counts/bin)...")
    if not run_grppha(src_pha, grp_pha, bkg_pha, arf_file,
                       rmf_path, min_counts):
        print(f"  ERROR: grppha failed.")
        return None

    # ---- Summary ----
    src_info = get_spectrum_info(src_pha)
    bkg_spec_info = get_spectrum_info(bkg_pha)
    grp_info = get_spectrum_info(grp_pha)

    summary = {
        'obsid': obsid,
        'stem': stem,
        'r_inner': r_inner,
        'r_outer': r_out,
        'src_counts': src_info['total_counts'],
        'bkg_counts': bkg_spec_info['total_counts'],
        'exposure': src_info['exposure'],
        'n_channels': src_info['n_channels'],
        'n_grouped_bins': grp_info.get('n_grouped_bins'),
        'backscal_src': src_info['backscal'],
        'backscal_bkg': bkg_spec_info['backscal'],
        'rmf': rmf_path,
        'count_rate': pileup_info.get('count_rate', 0),
        'date_obs': date_obs,
    }

    # ---- Write extraction log ----
    log_file = os.path.join(obsid_path, f'{stem}_extraction.log')
    with open(log_file, 'w') as f:
        f.write(f"# Spectral extraction log for {stem}\n")
        f.write(f"# Generated by swift_xrt_extract_spectra.py\n\n")
        f.write(f"Event file     : {os.path.basename(evt_file)}\n")
        f.write(f"Source RA/Dec  : {ra_src:.6f} / {dec_src:.6f}\n")
        f.write(f"Centroid (pix) : ({xc:.2f}, {yc:.2f})\n")
        f.write(f"Plate scale    : {plate_scale:.4f} arcsec/pix\n\n")
        f.write(f"Source region  : ")
        if r_inner > 0:
            f.write(f"annulus {r_inner:.1f}\" - {r_out:.1f}\"\n")
        else:
            f.write(f"circle {r_out:.1f}\"\n")
        f.write(f"Bkg region     : {bkg_info['type']} "
                f"({bkg_info['counts']} cts)\n\n")
        f.write(f"Source counts  : {summary['src_counts']}\n")
        f.write(f"Bkg counts     : {summary['bkg_counts']}\n")
        f.write(f"Exposure       : {summary['exposure']:.1f} s\n")
        f.write(f"Count rate     : "
                f"{summary.get('count_rate', 0):.3f} ct/s\n\n")
        f.write(f"BACKSCAL src   : {summary['backscal_src']:.6e}\n")
        f.write(f"BACKSCAL bkg   : {summary['backscal_bkg']:.6e}\n")
        f.write(f"Bkg scaling    : "
                f"{summary['backscal_src']/summary['backscal_bkg']:.4f}"
                f"\n\n") if summary['backscal_bkg'] > 0 else None
        f.write(f"Channels       : {summary['n_channels']}\n")
        f.write(f"Grouped bins   : {summary['n_grouped_bins']}\n")
        f.write(f"Min cts/bin    : {min_counts}\n\n")
        f.write(f"ARF            : {os.path.basename(arf_file)}\n")
        f.write(f"RMF            : {rmf_path}\n\n")
        f.write(f"Output files:\n")
        f.write(f"  {stem}_src.pha   (source spectrum)\n")
        f.write(f"  {stem}_bkg.pha   (background spectrum)\n")
        f.write(f"  {stem}.arf       (ancillary response)\n")
        f.write(f"  {stem}_grp.pha   (grouped, ready for XSPEC)\n")

    print(f"\n  --- Summary ---")
    print(f"  Source counts : {summary['src_counts']}")
    print(f"  Bkg counts    : {summary['bkg_counts']}")
    print(f"  Exposure      : {summary['exposure']:.1f} s")
    print(f"  Grouped bins  : {summary['n_grouped_bins']}")
    print(f"  RMF           : {rmf_path}")
    print(f"  Output        : {stem}_grp.pha (ready for XSPEC)")

    summary['mode'] = 'PC'
    return summary


# ---------------------------------------------------------------
# Process one WT-mode observation
# ---------------------------------------------------------------

def process_wt_observation(entry, min_counts,
                           wt_srcrad, wt_bkginner, wt_bkgouter):
    """
    Full extraction pipeline for a single WT-mode observation.

    Key differences from PC mode:
      - Source region is a circle (no pile-up exclusion needed
        below ~150 ct/s)
      - Background is an annulus symmetric about ~100 pixels
        (the WT window half-width)
      - BACKSCAL must be corrected for WT 1D geometry after
        extraction
      - Exposure map pattern uses *xwt* instead of *xpc*
    """
    obsid = entry['obsid']
    stem = entry['filename']
    obsid_path = os.path.join(BASE_DIR, obsid)

    print(f"\n{'='*65}")
    print(f"  OBSID: {obsid}    File: {stem}  [WT]")
    print(f"{'='*65}")

    # ---- Locate event file ----
    evt_file = None
    for ext in ['_cl.evt', '_cl.evt.gz']:
        candidate = os.path.join(obsid_path, stem + ext)
        if os.path.exists(candidate):
            evt_file = candidate
            break
    if evt_file is None:
        candidates = glob.glob(os.path.join(
            obsid_path, '**', stem + '_cl.evt*'), recursive=True)
        if candidates:
            evt_file = candidates[0]

    if evt_file is None:
        print(f"  ERROR: event file not found for {stem}")
        return None

    # ---- Read observation date from header ----
    with fits.open(evt_file) as hdul:
        evt_header = hdul[1].header
        date_obs = evt_header.get('DATE-OBS', 'N/A')
    print(f"  Date: {date_obs}")

    # ---- Read WT profile info (from viewer script output) ----
    profile_file = os.path.join(obsid_path,
                                 f'{stem}_wt_profile.txt')
    if os.path.exists(profile_file):
        wt_info = read_wt_profile_info(profile_file)
        xc = wt_info['xc']
        yc = wt_info['yc']
        src_radius = wt_info['src_radius']
        bkg_inner = wt_info['bkg_inner']
        bkg_outer = wt_info['bkg_outer']
        plate_scale = wt_info['plate_scale']
        print(f"  Using profile info from {stem}_wt_profile.txt")
    else:
        # Fall back to command-line defaults and find source
        # position from event file centroid
        print(f"  NOTE: no _wt_profile.txt found, using defaults.")
        src_radius = wt_srcrad
        bkg_inner = wt_bkginner
        bkg_outer = wt_bkgouter
        plate_scale = 2.36

        # Simple centroid from events
        with fits.open(evt_file) as hdul:
            data = hdul[1].data
            xc = float(np.median(data['X']))
            yc = float(np.median(data['Y']))

    if xc is None or yc is None:
        print(f"  ERROR: no source position available.")
        return None

    print(f"  Source: ({xc:.1f}, {yc:.1f})  "
          f"r={src_radius} pix")
    print(f"  Background: annulus {bkg_inner}-{bkg_outer} pix")

    # ---- Compute BACKSCAL values ----
    # For WT mode, BACKSCAL must reflect 1D extent in DETX
    backscal_src = 2 * src_radius
    backscal_bkg = bkg_outer - bkg_inner - 1
    print(f"  BACKSCAL: src={backscal_src}, bkg={backscal_bkg}")

    # ---- Write region files ----
    src_reg = os.path.join(obsid_path, f'{stem}_src.reg')
    write_wt_source_region(src_reg, xc, yc, src_radius)
    print(f"  Wrote: {os.path.basename(src_reg)}")

    bkg_reg = os.path.join(obsid_path, f'{stem}_bkg.reg')
    write_wt_background_region(bkg_reg, xc, yc,
                                bkg_inner, bkg_outer)
    print(f"  Wrote: {os.path.basename(bkg_reg)}")

    # ---- Extract source spectrum ----
    src_pha = os.path.join(obsid_path, f'{stem}_src.pha')
    print(f"  Extracting source spectrum...")
    if not extract_spectrum_xselect(evt_file, src_reg, src_pha,
                                     session_name="src"):
        print(f"  ERROR: source spectrum extraction failed.")
        return None

    # ---- Extract background spectrum ----
    bkg_pha = os.path.join(obsid_path, f'{stem}_bkg.pha')
    print(f"  Extracting background spectrum...")
    if not extract_spectrum_xselect(evt_file, bkg_reg, bkg_pha,
                                     session_name="bkg"):
        print(f"  ERROR: background spectrum extraction failed.")
        return None

    # ---- Correct BACKSCAL for WT 1D geometry ----
    # XSELECT writes BACKSCAL based on 2D region area, which is
    # wrong for WT mode where the data are effectively 1D.
    print(f"  Correcting BACKSCAL keywords...")
    fix_wt_backscal(src_pha, backscal_src)
    fix_wt_backscal(bkg_pha, backscal_bkg)

    # ---- Generate exposure map ----
    att_file, hk_file = find_auxiliary_files(obsid_path, obsid)
    expo_file = None
    if att_file and hk_file:
        print(f"  Generating exposure map...")
        expo_file = run_xrtexpomap(evt_file, att_file, hk_file,
                                    obsid_path)
    else:
        # Check for existing WT exposure map
        expo_candidates = glob.glob(os.path.join(
            obsid_path, '**', '*xwt*_ex.img*'), recursive=True)
        if expo_candidates:
            expo_file = expo_candidates[0]
            print(f"  Using existing exposure map: "
                  f"{os.path.basename(expo_file)}")
        else:
            print(f"  WARNING: cannot generate exposure map.")

    # ---- Generate ARF ----
    arf_file = os.path.join(obsid_path, f'{stem}.arf')
    rmf_path = None
    if expo_file:
        print(f"  Generating ARF...")
        arf_file, rmf_path = run_xrtmkarf(
            src_pha, expo_file, arf_file, src_reg)
    else:
        print(f"  WARNING: skipping ARF (no exposure map).")
        arf_file = None

    if arf_file is None or rmf_path is None:
        print(f"  ERROR: ARF/RMF generation failed.")
        return None

    # ---- Group spectrum ----
    grp_pha = os.path.join(obsid_path, f'{stem}_grp.pha')
    print(f"  Grouping spectrum (min {min_counts} counts/bin)...")
    if not run_grppha(src_pha, grp_pha, bkg_pha, arf_file,
                       rmf_path, min_counts):
        print(f"  ERROR: grppha failed.")
        return None

    # ---- Correct BACKSCAL in grouped file too ----
    fix_wt_backscal(grp_pha, backscal_src)

    # ---- Summary ----
    src_info = get_spectrum_info(src_pha)
    bkg_spec_info = get_spectrum_info(bkg_pha)
    grp_info = get_spectrum_info(grp_pha)

    summary = {
        'obsid': obsid,
        'stem': stem,
        'mode': 'WT',
        'r_inner': 0.0,
        'r_outer': float(src_radius),
        'src_counts': src_info['total_counts'],
        'bkg_counts': bkg_spec_info['total_counts'],
        'exposure': src_info['exposure'],
        'n_channels': src_info['n_channels'],
        'n_grouped_bins': grp_info.get('n_grouped_bins'),
        'backscal_src': backscal_src,
        'backscal_bkg': backscal_bkg,
        'rmf': rmf_path,
        'count_rate': src_info['total_counts'] / src_info['exposure']
            if src_info['exposure'] > 0 else 0,
        'date_obs': date_obs,
    }

    # ---- Write extraction log ----
    log_file = os.path.join(obsid_path, f'{stem}_extraction.log')
    with open(log_file, 'w') as f:
        f.write(f"# Spectral extraction log for {stem} [WT]\n")
        f.write(f"# Generated by swift_xrt_extract_spectra.py\n\n")
        f.write(f"Event file     : {os.path.basename(evt_file)}\n")
        f.write(f"Mode           : WT\n")
        f.write(f"Source (pix)   : ({xc:.2f}, {yc:.2f})\n")
        f.write(f"Plate scale    : {plate_scale:.4f} arcsec/pix\n\n")
        f.write(f"Source region  : circle r={src_radius} pix\n")
        f.write(f"Bkg region     : annulus "
                f"{bkg_inner}-{bkg_outer} pix\n\n")
        f.write(f"Source counts  : {summary['src_counts']}\n")
        f.write(f"Bkg counts     : {summary['bkg_counts']}\n")
        f.write(f"Exposure       : {summary['exposure']:.1f} s\n")
        f.write(f"Count rate     : "
                f"{summary['count_rate']:.3f} ct/s\n\n")
        f.write(f"BACKSCAL src   : {backscal_src} "
                f"(WT 1D: 2*r_src)\n")
        f.write(f"BACKSCAL bkg   : {backscal_bkg} "
                f"(WT 1D: r_out-r_in-1)\n\n")
        f.write(f"Channels       : {summary['n_channels']}\n")
        f.write(f"Grouped bins   : {summary['n_grouped_bins']}\n")
        f.write(f"Min cts/bin    : {min_counts}\n\n")
        f.write(f"ARF            : {os.path.basename(arf_file)}\n")
        f.write(f"RMF            : {rmf_path}\n\n")
        f.write(f"Output files:\n")
        f.write(f"  {stem}_src.pha   (source spectrum)\n")
        f.write(f"  {stem}_bkg.pha   (background spectrum)\n")
        f.write(f"  {stem}.arf       (ancillary response)\n")
        f.write(f"  {stem}_grp.pha   (grouped, ready for XSPEC)\n")

    print(f"\n  --- Summary ---")
    print(f"  Source counts : {summary['src_counts']}")
    print(f"  Bkg counts    : {summary['bkg_counts']}")
    print(f"  Exposure      : {summary['exposure']:.1f} s")
    print(f"  Grouped bins  : {summary['n_grouped_bins']}")
    print(f"  RMF           : {rmf_path}")
    print(f"  Output        : {stem}_grp.pha (ready for XSPEC)")

    return summary


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Extract Swift XRT spectra (PC and/or WT) '
                    'for XSPEC/Sherpa.')
    parser.add_argument('--ra', type=float, required=True,
                        help='Source RA in degrees')
    parser.add_argument('--dec', type=float, required=True,
                        help='Source Dec in degrees')
    parser.add_argument('--mode', type=str, default='both',
                        choices=['pc', 'wt', 'both'],
                        help='Which mode(s) to extract: '
                             'pc, wt, or both (default: both)')

    # PC-specific options
    parser.add_argument('--rout', type=str, default='47',
                        help='[PC] Outer extraction radius in '
                             'arcsec, or "auto" (default: 47)')
    parser.add_argument('--bkg-inner', type=float, default=100.0,
                        help='[PC] Background annulus inner radius '
                             'in arcsec (default: 100)')
    parser.add_argument('--bkg-outer', type=float, default=160.0,
                        help='[PC] Background annulus outer radius '
                             'in arcsec (default: 160)')

    # WT-specific options
    parser.add_argument('--wt-srcrad', type=int, default=20,
                        help='[WT] Source extraction radius in '
                             'pixels (default: 20)')
    parser.add_argument('--wt-bkginner', type=int, default=80,
                        help='[WT] Background annulus inner radius '
                             'in pixels (default: 80)')
    parser.add_argument('--wt-bkgouter', type=int, default=120,
                        help='[WT] Background annulus outer radius '
                             'in pixels (default: 120)')

    # Common options
    parser.add_argument('--mincounts', type=int, default=20,
                        help='Minimum counts per grouped bin '
                             '(default: 20 for chi-sq; 1 for cstat)')
    parser.add_argument('--pctable', type=str,
                        default='pc_master_table.txt',
                        help='PC master table filename')
    parser.add_argument('--wttable', type=str,
                        default='wt_master_table.txt',
                        help='WT master table filename')
    args = parser.parse_args()

    # --- Environment checks ---
    caldb = check_environment()

    # --- Read master tables ---
    pc_entries = []
    wt_entries = []

    if args.mode in ('pc', 'both'):
        pc_path = os.path.join(BASE_DIR, args.pctable)
        if os.path.exists(pc_path):
            pc_entries = read_master_table(pc_path)
            print(f"PC table: {len(pc_entries)} observations "
                  f"(from {args.pctable})")
        elif args.mode == 'pc':
            print(f"ERROR: {args.pctable} not found.")
            sys.exit(1)
        else:
            print(f"NOTE: {args.pctable} not found, skipping PC.")

    if args.mode in ('wt', 'both'):
        wt_path = os.path.join(BASE_DIR, args.wttable)
        if os.path.exists(wt_path):
            wt_entries = read_master_table(wt_path)
            print(f"WT table: {len(wt_entries)} observations "
                  f"(from {args.wttable})")
        elif args.mode == 'wt':
            print(f"ERROR: {args.wttable} not found.")
            sys.exit(1)
        else:
            print(f"NOTE: {args.wttable} not found, skipping WT.")

    if not pc_entries and not wt_entries:
        print("No observations to process.")
        sys.exit(0)

    print(f"\nSource: RA={args.ra:.6f}, Dec={args.dec:.6f}")
    if pc_entries:
        print(f"PC: outer radius={args.rout}, "
              f"bkg={args.bkg_inner}-{args.bkg_outer}\"")
    if wt_entries:
        print(f"WT: src r={args.wt_srcrad} pix, "
              f"bkg={args.wt_bkginner}-{args.wt_bkgouter} pix")
    print(f"Grouping: min {args.mincounts} counts/bin\n")

    # --- Process PC observations ---
    summaries = []

    if pc_entries:
        print(f"\n{'#'*65}")
        print(f"  PROCESSING PC MODE ({len(pc_entries)} observations)")
        print(f"{'#'*65}")

        for entry in pc_entries:
            result = process_pc_observation(
                entry, args.ra, args.dec, args.rout,
                args.mincounts, args.bkg_inner, args.bkg_outer)
            if result:
                summaries.append(result)

    # --- Process WT observations ---
    if wt_entries:
        print(f"\n{'#'*65}")
        print(f"  PROCESSING WT MODE ({len(wt_entries)} observations)")
        print(f"{'#'*65}")

        for entry in wt_entries:
            result = process_wt_observation(
                entry, args.mincounts,
                args.wt_srcrad, args.wt_bkginner,
                args.wt_bkgouter)
            if result:
                summaries.append(result)

    # --- Grand summary table ---
    if summaries:
        print(f"\n\n{'='*130}")
        print(f"  EXTRACTION SUMMARY")
        print(f"{'='*130}")
        hdr = (f"  {'OBSID':<14} {'File':<30} {'Mode':>4} "
               f"{'Date-Obs':<12} "
               f"{'Rin':>5} {'Rout':>5} "
               f"{'SrcCts':>7} {'BkgCts':>7} {'Exp(s)':>8} "
               f"{'Bins':>5}  {'RMF'}")
        print(hdr)
        print(f"  {'-'*128}")

        total_src = 0
        total_bkg = 0
        total_exp = 0.0

        for s in summaries:
            mode = s.get('mode', 'PC')
            date_short = s.get('date_obs', 'N/A')[:10]
            rmf_short = os.path.basename(
                s.get('rmf', 'N/A').replace('$CALDB/', ''))
            print(f"  {s['obsid']:<14} {s['stem']:<30} {mode:>4} "
                  f"{date_short:<12} "
                  f"{s['r_inner']:>5.1f} {s['r_outer']:>5.1f} "
                  f"{s['src_counts']:>7} {s['bkg_counts']:>7} "
                  f"{s['exposure']:>8.1f} "
                  f"{s['n_grouped_bins'] or 'N/A':>5}  "
                  f"{rmf_short}")
            total_src += s['src_counts']
            total_bkg += s['bkg_counts']
            total_exp += s['exposure']

        print(f"  {'-'*128}")
        print(f"  {'TOTAL':<14} {'':30} {'':>4} "
              f"{'':12} "
              f"{'':>5} {'':>5} "
              f"{total_src:>7} {total_bkg:>7} "
              f"{total_exp:>8.1f}")
        print(f"{'='*130}")

        n_pc = sum(1 for s in summaries if s.get('mode') == 'PC')
        n_wt = sum(1 for s in summaries if s.get('mode') == 'WT')
        print(f"\n  {len(summaries)} spectra extracted "
              f"({n_pc} PC, {n_wt} WT).")
        print(f"  Total source counts: {total_src}")
        print(f"  Total exposure: {total_exp:.1f} s "
              f"({total_exp/1000:.2f} ks)")
    else:
        print("\nNo spectra were successfully extracted.")


if __name__ == '__main__':
    main()

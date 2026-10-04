#!/usr/bin/env python3
"""
swift_xrt_summary.py

Crawl subdirectories (named by Swift OBSID) in the current working directory,
find all cleaned level-2 XRT event files (*_cl.evt), and produce a summary
table showing the target, pointing offset, mode sequence, start/end times,
and exposure for each.

Usage:
    python swift_xrt_summary.py
    python swift_xrt_summary.py --compact --ra 187.2779 --dec 2.0524

    Run from the directory that contains OBSID subdirectories. The script
    searches recursively for *_cl.evt files within each OBSID directory.
    Pointing offsets are measured from --ra/--dec if given, otherwise from
    each observation's own target position (RA_OBJ/DEC_OBJ).

Requirements:
    astropy (for FITS header reading)
"""

import os
import sys
import glob
import math
import re
import argparse
from collections import defaultdict

try:
    from astropy.io import fits
except ImportError:
    print("ERROR: astropy is required. Install with: pip install astropy")
    sys.exit(1)


def separation_arcmin(ra1, dec1, ra2, dec2):
    """Angular separation in arcmin between two positions in degrees."""
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1, dec1, ra2, dec2))
    # Haversine form, accurate at the arcminute scales of interest here
    h = (math.sin((dec2 - dec1) / 2) ** 2
         + math.cos(dec1) * math.cos(dec2) * math.sin((ra2 - ra1) / 2) ** 2)
    return math.degrees(2 * math.asin(math.sqrt(h))) * 60.0


def identify_mode(filename):
    """
    Identify the XRT mode and observation type from the event filename.

    Swift XRT cleaned event files follow a naming convention like:
        sw<OBSID>x{pc,wt}w{N}{st,sl,po,...}_cl.evt

    Returns a human-readable string like 'WT_SETTLING', 'PC_POINTED', etc.
    """
    basename = os.path.basename(filename).lower()

    # Determine mode: PC or WT
    if 'xpc' in basename:
        mode = 'PC'
    elif 'xwt' in basename:
        mode = 'WT'
    else:
        mode = 'UNKNOWN'

    # Determine observation type from the two-letter code before _cl
    # Common codes: st = settling, sl = slew, po = pointed
    # Extract the two-letter type code
    match = re.search(r'x(?:pc|wt)w\d([a-z]{2})_cl\.evt', basename)
    if match:
        type_code = match.group(1)
        type_map = {
            'st': 'SETTLING',
            'sl': 'SLEW',
            'po': 'POINTED',
        }
        obs_type = type_map.get(type_code, type_code.upper())
    else:
        obs_type = 'UNKNOWN'

    return f"{mode}_{obs_type}"


def get_event_info(filepath):
    """
    Read key header information from a cleaned event file,
    including GTI (Good Time Interval) analysis to detect
    multiple orbit segments.

    Returns a dict with mode, start time, end time, exposure,
    event count, and GTI information.
    """
    info = {
        'filepath': filepath,
        'filename': os.path.basename(filepath),
        'mode': identify_mode(filepath),
    }

    try:
        with fits.open(filepath) as hdul:
            # The EVENTS extension is typically extension 1
            header = hdul[1].header

            info['date_obs'] = header.get('DATE-OBS', 'N/A')
            info['date_end'] = header.get('DATE-END', 'N/A')
            info['exposure'] = header.get('EXPOSURE', 0.0)
            info['datamode'] = header.get('DATAMODE', 'N/A')
            info['obs_id'] = header.get('OBS_ID', 'N/A')

            # Target name and positions, for the pointing offset
            info['object'] = str(header.get('OBJECT', 'N/A')).strip()
            info['ra_pnt'] = header.get('RA_PNT')
            info['dec_pnt'] = header.get('DEC_PNT')
            info['ra_obj'] = header.get('RA_OBJ')
            info['dec_obj'] = header.get('DEC_OBJ')

            # Get number of events (rows in the EVENTS table)
            info['nevents'] = header.get('NAXIS2', 0)

            # Calculate approximate count rate
            if info['exposure'] > 0:
                info['count_rate'] = info['nevents'] / info['exposure']
            else:
                info['count_rate'] = 0.0

            # --- GTI analysis ---
            # The GTI extension lists time intervals when data
            # was collected. Multiple GTIs within a single event
            # file indicate multiple orbit segments separated by
            # Earth occultation gaps (~60 min Swift orbit period).
            info['n_gti'] = 0
            info['gti_durations'] = []
            info['gti_gaps'] = []
            info['total_elapsed'] = 0.0

            # Find the GTI extension
            gti_ext = None
            for i, ext in enumerate(hdul):
                if ext.name == 'GTI':
                    gti_ext = i
                    break

            if gti_ext is not None:
                gti_data = hdul[gti_ext].data
                if gti_data is not None and len(gti_data) > 0:
                    starts = gti_data['START']
                    stops = gti_data['STOP']
                    info['n_gti'] = len(starts)

                    # Duration of each GTI segment
                    durations = stops - starts
                    info['gti_durations'] = durations.tolist()

                    # Gaps between consecutive GTI segments
                    # (these are the Earth occultation periods)
                    if len(starts) > 1:
                        gaps = starts[1:] - stops[:-1]
                        info['gti_gaps'] = gaps.tolist()

                    # Total elapsed time from first start to
                    # last stop (includes gaps)
                    info['total_elapsed'] = \
                        float(stops[-1] - starts[0])

    except Exception as e:
        print(f"  WARNING: Could not read {filepath}: {e}")
        return None

    return info


def format_exposure(seconds):
    """Format exposure time in a readable way."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds:.1f}s ({seconds/60:.1f}m)"
    else:
        return f"{seconds:.1f}s ({seconds/3600:.2f}h)"


def pointing_summary(file_infos, ref=None):
    """
    Target name and pointing offset (arcmin) for one OBSID.

    The offset is the angle between where Swift pointed (RA_PNT/DEC_PNT)
    and ``ref`` (the --ra/--dec position) or, without one, the
    observation's own target (RA_OBJ/DEC_OBJ). A large offset puts the
    source far off-axis; in WT mode, more than ~4' can put it outside
    the 8'-wide window, as for 3C 273 OBSID 00035017041 (5.1').
    Returns (target, offset), with offset None if a position is missing.
    """
    pointed = [fi for fi in file_infos if 'POINTED' in fi['mode']]
    info = (pointed or file_infos)[0]
    if ref is None:
        ref = (info.get('ra_obj'), info.get('dec_obj'))
    pnt = (info.get('ra_pnt'), info.get('dec_pnt'))
    if None in ref or None in pnt:
        return info.get('object', 'N/A'), None
    return info.get('object', 'N/A'), separation_arcmin(*pnt, *ref)


def print_obsid_table(obsid, file_infos, ref=None, ref_label='target'):
    """Print a formatted table for a single OBSID."""
    # Sort by start time
    file_infos.sort(key=lambda x: x['date_obs'])

    # Header
    target, offset = pointing_summary(file_infos, ref)
    off_str = f"{offset:.1f}'" if offset is not None else 'N/A'
    print(f"\n{'='*100}")
    print(f"  OBSID: {obsid}    Target: {target}    "
          f"Pointing offset: {off_str} from {ref_label}")
    print(f"{'='*100}")

    # Column headers
    print(f"  {'#':<4} {'Mode':<14} {'DATAMODE':<12} {'Start (UTC)':<22} "
          f"{'End (UTC)':<22} {'Exposure':>12} {'Events':>8} {'ct/s':>8}")
    print(f"  {'-'*4} {'-'*14} {'-'*12} {'-'*22} {'-'*22} "
          f"{'-'*12} {'-'*8} {'-'*8}")

    total_exposure = 0.0
    total_events = 0

    for i, info in enumerate(file_infos, 1):
        exp_str = format_exposure(info['exposure'])
        rate_str = f"{info['count_rate']:.2f}" if info['count_rate'] > 0 else "N/A"

        print(f"  {i:<4} {info['mode']:<14} {info['datamode']:<12} "
              f"{info['date_obs']:<22} {info['date_end']:<22} "
              f"{exp_str:>12} {info['nevents']:>8} {rate_str:>8}")

        total_exposure += info['exposure']
        total_events += info['nevents']

    print(f"  {'-'*4} {'-'*14} {'-'*12} {'-'*22} {'-'*22} "
          f"{'-'*12} {'-'*8} {'-'*8}")
    print(f"  {'':4} {'TOTAL':<14} {'':12} {'':22} {'':22} "
          f"{format_exposure(total_exposure):>12} {total_events:>8}")

    # Show the sequence summary
    sequence = " -> ".join(info['mode'] for info in file_infos)
    print(f"\n  Sequence: {sequence}")

    # Pile-up warning for PC mode
    for info in file_infos:
        if 'PC' in info['mode'] and info['count_rate'] > 0.5:
            print(f"\n  *** WARNING: {info['mode']} has {info['count_rate']:.2f} ct/s "
                  f"-- possible pile-up (threshold ~0.5 ct/s) ***")
        if 'WT' in info['mode'] and info['count_rate'] > 150:
            print(f"\n  *** WARNING: {info['mode']} has {info['count_rate']:.2f} ct/s "
                  f"-- possible pile-up (threshold ~150 ct/s) ***")

    # GTI / orbit analysis
    # Show orbit segment details for each event file
    for info in file_infos:
        n_gti = info.get('n_gti', 0)
        if n_gti > 1:
            durations = info.get('gti_durations', [])
            gaps = info.get('gti_gaps', [])
            elapsed = info.get('total_elapsed', 0)
            duty = (info['exposure'] / elapsed * 100) \
                if elapsed > 0 else 0

            print(f"\n  {info['mode']} — {n_gti} orbit segments "
                  f"(elapsed: {elapsed:.0f}s, "
                  f"duty cycle: {duty:.0f}%):")
            for j, dur in enumerate(durations):
                gap_str = ""
                if j < len(gaps):
                    gap_min = gaps[j] / 60.0
                    gap_str = f"  [gap: {gap_min:.1f} min]"
                print(f"    Segment {j+1}: {dur:.1f}s "
                      f"({dur/60:.1f} min){gap_str}")
        elif n_gti == 1:
            print(f"\n  {info['mode']} — single orbit segment "
                  f"({info['exposure']:.1f}s)")


def get_compact_row(obsid, file_infos, ref=None):
    """
    Summarize a single OBSID into a compact row dict.
    """
    file_infos.sort(key=lambda x: x['date_obs'])
    target, offset = pointing_summary(file_infos, ref)

    total_exposure = sum(fi['exposure'] for fi in file_infos)
    date_start = file_infos[0]['date_obs']

    # Identify slew at start and end
    slew_start = ''
    slew_end = ''
    if 'SLEW' in file_infos[0]['mode']:
        slew_start = file_infos[0]['mode'].split('_')[0]  # 'PC' or 'WT'
    if len(file_infos) > 1 and 'SLEW' in file_infos[-1]['mode']:
        slew_end = file_infos[-1]['mode'].split('_')[0]

    # Count WT and PC pointed observations (excluding slew/settling)
    n_wt_pointed = sum(1 for fi in file_infos
                       if fi['mode'] == 'WT_POINTED')
    n_pc_pointed = sum(1 for fi in file_infos
                       if fi['mode'] == 'PC_POINTED')

    # Total WT pointed exposure
    wt_pointed_exp = sum(fi['exposure'] for fi in file_infos
                         if fi['mode'] == 'WT_POINTED')

    # Total PC pointed exposure
    pc_exp = sum(fi['exposure'] for fi in file_infos
                 if fi['mode'] == 'PC_POINTED')

    # Build sequence code:
    # 1 = WT_SLEW, 2 = PC_SLEW
    # 3 = WT_SETTLING, 4 = PC_SETTLING
    # 5 = WT_POINTED, 6 = PC_POINTED
    def _seq_code(mode_str):
        parts = mode_str.split('_')
        prefix = parts[0]  # WT or PC
        obs_type = parts[1] if len(parts) > 1 else 'UNKNOWN'
        codes = {
            ('WT', 'SLEW'): '1', ('PC', 'SLEW'): '2',
            ('WT', 'SETTLING'): '3', ('PC', 'SETTLING'): '4',
            ('WT', 'POINTED'): '5', ('PC', 'POINTED'): '6',
        }
        return codes.get((prefix, obs_type), '?')

    sequence = ''
    for fi in file_infos:
        sequence += _seq_code(fi['mode'])

    # Total count rate across all files
    total_events = sum(fi['nevents'] for fi in file_infos)
    total_cts_per_s = total_events / total_exposure if total_exposure > 0 else 0.0

    # Max orbit segments across pointed event files
    # (indicates multi-orbit observations)
    pointed_gtis = [fi.get('n_gti', 1) for fi in file_infos
                    if 'POINTED' in fi['mode']]
    max_orbits = max(pointed_gtis) if pointed_gtis else 1

    return {
        'obsid': obsid,
        'date_start': date_start,
        'total_exp_ks': total_exposure / 1000.0,
        'total_cts_per_s': total_cts_per_s,
        'slew_start': slew_start,
        'slew_end': slew_end,
        'n_wt_pointed': n_wt_pointed,
        'n_pc': n_pc_pointed,
        'wt_pointed_exp_ks': wt_pointed_exp / 1000.0,
        'pc_exp_ks': pc_exp / 1000.0,
        'sequence': sequence,
        'max_orbits': max_orbits,
        'offset': offset,
        'target': target,
    }


def print_compact_table(rows, ref_label='target', no_files=()):
    """Print the compact one-row-per-OBSID table."""
    # Column headers. Target goes last: names like
    # SDSSJ122933.69+015810.4 are too long for a fixed-width column.
    off_hdr = "Off(')"
    tw = max(len('Target'), max(len(r['target']) for r in rows))
    hdr = (f"  {'OBSID':<14} {'Date/Time':<22} {'Total(ks)':>10} {'ct/s':>8} "
           f"{'Slew_i':>7} {'Slew_f':>7} {'N_WT':>5} {'N_PC':>5} "
           f"{'WT_exp(ks)':>11} {'PC_exp(ks)':>11} {'Orb':>4} {'Seq':>6} "
           f"{off_hdr:>6}  {'Target':<{tw}}")
    sep = (f"  {'-'*14} {'-'*22} {'-'*10} {'-'*8} "
           f"{'-'*7} {'-'*7} {'-'*5} {'-'*5} {'-'*11} {'-'*11} {'-'*4} {'-'*6} "
           f"{'-'*6}  {'-'*tw}")
    width = len(hdr)

    print(f"\n{'='*width}")
    print(f"  COMPACT SUMMARY")
    print(f"  Sequence codes: 1=WT_SLEW  2=PC_SLEW  3=WT_SETTLING  4=PC_SETTLING  5=WT_POINTED  6=PC_POINTED")
    print(f"  Off(') = pointing offset in arcmin from {ref_label}")
    print(f"{'='*width}")
    print(hdr)
    print(sep)

    total_exp = 0.0
    total_wt_exp = 0.0
    total_pc_exp = 0.0

    for r in rows:
        off = f"{r['offset']:.1f}" if r['offset'] is not None else 'N/A'
        print(f"  {r['obsid']:<14} {r['date_start']:<22} {r['total_exp_ks']:>10.3f} "
              f"{r['total_cts_per_s']:>8.2f} "
              f"{r['slew_start']:>7} {r['slew_end']:>7} {r['n_wt_pointed']:>5} "
              f"{r['n_pc']:>5} {r['wt_pointed_exp_ks']:>11.3f} "
              f"{r['pc_exp_ks']:>11.3f} {r['max_orbits']:>4} {r['sequence']:>6} "
              f"{off:>6}  {r['target']}")
        total_exp += r['total_exp_ks']
        total_wt_exp += r['wt_pointed_exp_ks']
        total_pc_exp += r['pc_exp_ks']

    print(sep)
    print(f"  {'TOTAL':<14} {'':22} {total_exp:>10.3f} "
          f"{'':>8} "
          f"{'':>7} {'':>7} {'':>5} {'':>5} {total_wt_exp:>11.3f} "
          f"{total_pc_exp:>11.3f}")
    print(f"{'='*width}")
    # An OBSID whose Step 3 failed has no cleaned files; say so rather
    # than leave it silently missing from the table.
    if no_files:
        print(f"  Not in the table (no cleaned event files): "
              f"{', '.join(no_files)}")


def main():
    parser = argparse.ArgumentParser(
        description='Summarize Swift XRT cleaned event files across OBSID directories.')
    parser.add_argument('--compact', action='store_true',
                        help='Output a single-row-per-OBSID summary table')
    parser.add_argument('--ra', type=float, default=None,
                        help='Source RA in decimal degrees; pointing '
                             'offsets are measured from --ra/--dec '
                             '(default: from each observation\'s target)')
    parser.add_argument('--dec', type=float, default=None,
                        help='Source Dec in decimal degrees')
    args = parser.parse_args()
    if (args.ra is None) != (args.dec is None):
        parser.error('give --ra and --dec together')
    ref = (args.ra, args.dec) if args.ra is not None else None
    ref_label = '--ra/--dec' if ref else 'target'

    script_dir = os.getcwd()

    # Find all subdirectories that look like Swift OBSIDs (8-11 digit numbers)
    obsid_pattern = re.compile(r'^\d{8,11}$')
    obsid_dirs = sorted([
        d for d in os.listdir(script_dir)
        if os.path.isdir(os.path.join(script_dir, d)) and obsid_pattern.match(d)
    ])

    if not obsid_dirs:
        print("No OBSID directories found in the current directory.")
        print("Expected directories named with 8-11 digit Swift OBSIDs (e.g., 00031659107)")
        sys.exit(1)

    print(f"\nFound {len(obsid_dirs)} OBSID director{'y' if len(obsid_dirs)==1 else 'ies'}.")

    grand_total_exposure = 0.0
    grand_total_events = 0
    obsids_processed = 0
    compact_rows = []
    no_files = []

    for obsid in obsid_dirs:
        obsid_path = os.path.join(script_dir, obsid)

        # Search recursively for cleaned event files
        cl_files = glob.glob(os.path.join(obsid_path, '**', '*_cl.evt'), recursive=True)
        # Also check for .evt.gz
        cl_files += glob.glob(os.path.join(obsid_path, '**', '*_cl.evt.gz'), recursive=True)

        if not cl_files:
            no_files.append(obsid)
            if not args.compact:
                print(f"\n  OBSID {obsid}: No cleaned event files found.")
            continue

        # Read info from each file
        file_infos = []
        for f in cl_files:
            info = get_event_info(f)
            if info is not None:
                file_infos.append(info)

        if file_infos:
            if args.compact:
                compact_rows.append(get_compact_row(obsid, file_infos, ref))
            else:
                print_obsid_table(obsid, file_infos, ref, ref_label)

            obsids_processed += 1
            grand_total_exposure += sum(fi['exposure'] for fi in file_infos)
            grand_total_events += sum(fi['nevents'] for fi in file_infos)
        else:
            no_files.append(obsid)

    if args.compact:
        if compact_rows:
            print_compact_table(compact_rows, ref_label, no_files)
        else:
            print("\nNo cleaned event files found in any OBSID directory.")
    else:
        # Grand summary
        if obsids_processed > 1:
            print(f"\n{'='*100}")
            print(f"  GRAND TOTAL: {obsids_processed} OBSIDs, "
                  f"{format_exposure(grand_total_exposure)} total exposure, "
                  f"{grand_total_events} total events")
            print(f"{'='*100}")
        elif obsids_processed == 0:
            print("\nNo cleaned event files found in any OBSID directory.")


if __name__ == '__main__':
    main()

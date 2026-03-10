#!/usr/bin/env python3
"""
make_wt_master_table.py

Generate a master table for WT-mode pointed observations,
analogous to pc_master_table.txt. Each row lists one event
file with an include flag that the user can edit to deselect
observations.

By default, observations with exposure < 20 seconds are
marked include=no. All others are marked include=yes.

Usage:
    python make_wt_master_table.py
    python make_wt_master_table.py --expmin 50
    python make_wt_master_table.py --output my_wt_table.txt

Output format (whitespace-separated, comment in quotes):
    OBSID  filename  include  exp(s)  ct/s  n_gti  "comment"
"""

import os
import sys
import re
import glob
import argparse

try:
    from astropy.io import fits
except ImportError:
    print("ERROR: astropy is required.")
    sys.exit(1)


BASE_DIR = os.path.abspath(os.getcwd())


def find_wt_po_files():
    """Find all WT pointed cleaned event files."""
    obsid_pattern = re.compile(r'^\d{8,11}$')
    obsid_dirs = sorted([
        d for d in os.listdir(BASE_DIR)
        if os.path.isdir(os.path.join(BASE_DIR, d))
        and obsid_pattern.match(d)
    ])

    results = []
    for obsid in obsid_dirs:
        obsid_path = os.path.join(BASE_DIR, obsid)
        wt_files = glob.glob(os.path.join(
            obsid_path, '**', '*xwt*po*_cl.evt'), recursive=True)
        wt_files += glob.glob(os.path.join(
            obsid_path, '**', '*xwt*po*_cl.evt.gz'), recursive=True)
        wt_files = sorted(set(wt_files))

        for wf in wt_files:
            stem = os.path.basename(wf).replace(
                '_cl.evt.gz', '').replace('_cl.evt', '')
            results.append({
                'obsid': obsid,
                'filepath': wf,
                'stem': stem,
            })
    return results


def get_file_info(filepath):
    """Read exposure, event count, and GTI count from event file."""
    info = {}
    try:
        with fits.open(filepath) as hdul:
            header = hdul[1].header
            info['exposure'] = header.get('EXPOSURE', 0.0)
            info['nevents'] = header.get('NAXIS2', 0)
            info['count_rate'] = info['nevents'] / info['exposure'] \
                if info['exposure'] > 0 else 0

            # Count GTI segments
            n_gti = 0
            for ext in hdul:
                if ext.name == 'GTI':
                    if ext.data is not None:
                        n_gti = len(ext.data)
                    break
            info['n_gti'] = n_gti
    except Exception as e:
        print(f"  WARNING: could not read {filepath}: {e}")
        info = {'exposure': 0, 'nevents': 0,
                'count_rate': 0, 'n_gti': 0}
    return info


def main():
    parser = argparse.ArgumentParser(
        description='Generate WT master table for pointed '
                    'observations.')
    parser.add_argument('--expmin', type=float, default=20.0,
                        help='Minimum exposure (s) for default '
                             'include=yes (default: 20)')
    parser.add_argument('--output', type=str,
                        default='wt_master_table.txt',
                        help='Output filename '
                             '(default: wt_master_table.txt)')
    args = parser.parse_args()

    entries = find_wt_po_files()
    if not entries:
        print("No WT pointed event files found.")
        sys.exit(1)

    print(f"Found {len(entries)} WT pointed event file(s).")
    print(f"Exposure threshold for include=yes: {args.expmin}s")

    # Read metadata for each file
    n_yes = 0
    n_no = 0
    lines = []

    for entry in entries:
        info = get_file_info(entry['filepath'])
        exposure = info['exposure']
        ct_rate = info['count_rate']
        n_gti = info['n_gti']

        if exposure >= args.expmin:
            include = 'yes'
            n_yes += 1
        else:
            include = 'no'
            n_no += 1

        lines.append({
            'obsid': entry['obsid'],
            'stem': entry['stem'],
            'include': include,
            'exposure': exposure,
            'ct_rate': ct_rate,
            'n_gti': n_gti,
        })

    # Write table
    output_path = os.path.join(BASE_DIR, args.output)
    with open(output_path, 'w') as f:
        # Header
        f.write(f"{'OBSID':<14} {'filename':<35} "
                f"{'include':<10} {'exp(s)':>8} {'ct/s':>8} "
                f"{'n_gti':>6} {'comment'}\n")

        for row in lines:
            f.write(f"{row['obsid']:<14} {row['stem']:<35} "
                    f"{row['include']:<10} "
                    f"{row['exposure']:>8.1f} "
                    f"{row['ct_rate']:>8.2f} "
                    f"{row['n_gti']:>6} "
                    f'""\n')

    print(f"\nWrote {args.output}:")
    print(f"  {n_yes} observations include=yes "
          f"(exposure >= {args.expmin}s)")
    print(f"  {n_no} observations include=no "
          f"(exposure < {args.expmin}s)")
    print(f"\nEdit the 'include' column to deselect observations.")


if __name__ == '__main__':
    main()

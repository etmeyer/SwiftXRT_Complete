#!/usr/bin/env python3
"""
swift_wt_summary_viewer.py

Summary table and viewer for Swift XRT WT-mode pointed cleaned
event files. For each observation:
  1. Reports exposure, count rate, date
  2. Extracts the 1D DETX profile (the physically meaningful
     cross-strip direction, immune to roll angle changes)
  3. Estimates source extraction radius and background annulus
  4. Produces sky-coordinate images with region overlays
  5. Outputs a multi-page PDF and per-observation text files

Extraction regions follow the UK Swift Science Data Centre
standard (https://www.swift.ac.uk/analysis/xrt/backscal.php):
  - Source: circular region (default 20 pix radius)
  - Background: annulus centered on source, symmetric about
    100 pixels from center (default 80-120 pix)
  - BACKSCAL: must be manually set to reflect 1D extent
    (source = 2*r_src, bkg = r_outer - r_inner - 1)

Usage:
    python swift_wt_summary_viewer.py --ra 187.2779 --dec 2.0524
    python swift_wt_summary_viewer.py --ra 187.2779 --dec 2.0524 \\
        --srcrad 20 --bkginner 80 --bkgouter 120
    python swift_wt_summary_viewer.py --ra 187.2779 --dec 2.0524 \\
        --compact --nmax 10

Requirements:
    astropy, matplotlib, numpy, scipy
"""

import os
import sys
import glob
import re
import argparse
import numpy as np

try:
    from astropy.io import fits
    from astropy.time import Time
except ImportError:
    print("ERROR: astropy is required.")
    sys.exit(1)

try:
    from scipy.ndimage import gaussian_filter1d
except ImportError:
    print("ERROR: scipy is required.")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Circle
    from matplotlib.colors import PowerNorm
except ImportError:
    print("ERROR: matplotlib is required.")
    sys.exit(1)


BASE_DIR = os.path.abspath(os.getcwd())


# ---------------------------------------------------------------
# Find WT pointed cleaned event files
# ---------------------------------------------------------------

def find_wt_files():
    """
    Find all WT-mode pointed cleaned event files.
    Only 'po' (pointed) files are used:
      - 'sl' = slew (excluded)
      - 'st' = settling (excluded)
      - 'po' = pointed, stable observation (what we want)
    """
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
            results.append({
                'obsid': obsid,
                'filepath': wf,
                'filename': os.path.basename(wf),
                'stem': os.path.basename(wf).replace(
                    '_cl.evt.gz', '').replace('_cl.evt', ''),
            })
    return results


# ---------------------------------------------------------------
# Read event file metadata and coordinates
# ---------------------------------------------------------------

def get_wt_metadata(filepath):
    """
    Extract metadata and event coordinates from a WT event file.
    Reads both sky (X, Y) and detector (DETX) coordinates.
    """
    info = {}
    with fits.open(filepath) as hdul:
        evt_ext = None
        for i, ext in enumerate(hdul):
            if ext.name == 'EVENTS':
                evt_ext = i
                break
        if evt_ext is None:
            evt_ext = 1

        header = hdul[evt_ext].header
        data = hdul[evt_ext].data

        info['exposure'] = header.get('EXPOSURE', 0)
        info['n_events'] = header.get('NAXIS2', 0)
        info['datamode'] = header.get('DATAMODE', 'N/A')
        info['date_obs'] = header.get('DATE-OBS', None)
        info['date_end'] = header.get('DATE-END', None)

        # Count rate
        info['count_rate'] = info['n_events'] / info['exposure'] \
            if info['exposure'] > 0 else 0

        # Midpoint time
        if info['date_obs'] and info['date_end']:
            try:
                t0 = Time(info['date_obs'], format='isot', scale='utc')
                t1 = Time(info['date_end'], format='isot', scale='utc')
                t_mid = t0 + (t1 - t0) / 2.0
                info['date_mid'] = t_mid.isot
                info['mjd_mid'] = t_mid.mjd
            except Exception:
                info['date_mid'] = info['date_obs']
                info['mjd_mid'] = None
        else:
            info['date_mid'] = info.get('date_obs', 'N/A')
            info['mjd_mid'] = None

        # WCS for sky X column
        tfields = header.get('TFIELDS', 0)
        x_col = None
        for ci in range(1, tfields + 1):
            ttype = header.get(f'TTYPE{ci}', '').strip().upper()
            if ttype == 'X':
                x_col = ci
                break

        if x_col:
            info['tcrvl_x'] = header.get(f'TCRVL{x_col}', 0)
            info['tcrpx_x'] = header.get(f'TCRPX{x_col}', 500)
            info['tcdlt_x'] = header.get(f'TCDLT{x_col}',
                                          -0.0006548089)
            info['plate_scale'] = abs(info['tcdlt_x']) * 3600
        else:
            info['plate_scale'] = 2.36

        # Read all coordinate columns
        info['x_events'] = data['X'].astype(float)
        info['y_events'] = data['Y'].astype(float)
        info['time_events'] = data['TIME'].astype(float)

        # DETX for the 1D cross-strip profile
        if 'DETX' in data.columns.names:
            info['detx_events'] = data['DETX'].astype(float)
        else:
            info['detx_events'] = None

        # --- GTI (Good Time Intervals) ---
        # Multiple GTIs indicate multiple orbit segments
        info['gti_starts'] = []
        info['gti_stops'] = []
        info['n_gti'] = 0
        gti_ext = None
        for i, ext in enumerate(hdul):
            if ext.name == 'GTI':
                gti_ext = i
                break
        if gti_ext is not None:
            gti_data = hdul[gti_ext].data
            if gti_data is not None and len(gti_data) > 0:
                info['gti_starts'] = gti_data['START'].tolist()
                info['gti_stops'] = gti_data['STOP'].tolist()
                info['n_gti'] = len(info['gti_starts'])

    return info


# ---------------------------------------------------------------
# Find source position in sky coordinates
# ---------------------------------------------------------------

def find_sky_position(x_events, y_events, ra_src, dec_src, info):
    """
    Convert RA/Dec to sky pixel coordinates and refine by
    finding the event centroid nearby.
    """
    # WCS conversion
    x_guess = info['tcrpx_x'] + \
        (ra_src - info['tcrvl_x']) / info['tcdlt_x']
    # For Y, find the column
    # Approximate: use centroid of events near x_guess
    near_x = np.abs(x_events - x_guess) < 15
    if np.sum(near_x) > 5:
        y_guess = np.median(y_events[near_x])
    else:
        y_guess = np.median(y_events)

    # Refine centroid iteratively
    xc, yc = x_guess, y_guess
    for _ in range(3):
        dist = np.sqrt((x_events - xc)**2 + (y_events - yc)**2)
        near = dist < 15
        if np.sum(near) < 5:
            break
        xc = np.mean(x_events[near])
        yc = np.mean(y_events[near])

    return xc, yc


# ---------------------------------------------------------------
# Estimate BACKSCAL values
# ---------------------------------------------------------------

def compute_backscal(src_radius, bkg_inner, bkg_outer):
    """
    Compute the correct BACKSCAL values for WT mode.

    In WT mode, BACKSCAL should reflect the 1D extent of the
    extraction region in the DETX direction, NOT the 2D area.
    (See https://www.swift.ac.uk/analysis/xrt/backscal.php)

    Source BACKSCAL = 2 * r_src (diameter in pixels)
    Bkg BACKSCAL    = r_outer - r_inner - 1
      (minus 1 for the end-of-window bad pixel; approximate)

    Returns:
        backscal_src, backscal_bkg
    """
    backscal_src = 2 * src_radius
    backscal_bkg = bkg_outer - bkg_inner - 1

    return backscal_src, backscal_bkg


# ---------------------------------------------------------------
# Sky coordinate image with regions
# ---------------------------------------------------------------

def make_combined_plot(stem, x_events, y_events, detx_events,
                       xc, yc, src_radius, bkg_inner, bkg_outer,
                       plate_scale, exposure, count_rate,
                       n_events, date_mid, n_gti):
    """
    Combined page: sky image (left) + DETX cross-strip profile (right).
    """
    fig = plt.figure(figsize=(14, 7))
    fig.patch.set_facecolor('black')

    # --- Left: sky coordinate image with regions ---
    ax1 = fig.add_subplot(121)

    margin = max(bkg_outer + 30, 150)
    x_lo = int(xc - margin)
    x_hi = int(xc + margin)
    y_lo = int(yc - margin)
    y_hi = int(yc + margin)

    mask = ((x_events >= x_lo) & (x_events < x_hi) &
            (y_events >= y_lo) & (y_events < y_hi))

    if np.sum(mask) > 0:
        image, _, _ = np.histogram2d(
            y_events[mask], x_events[mask],
            bins=[np.arange(y_lo, y_hi + 1),
                  np.arange(x_lo, x_hi + 1)])
        vmax = max(np.max(image), 1)
        ax1.imshow(image, origin='lower', cmap='viridis',
                   extent=[x_lo, x_hi, y_lo, y_hi],
                   interpolation='nearest',
                   norm=PowerNorm(gamma=0.5, vmin=0, vmax=vmax),
                   aspect='equal')

    # Source circle (magenta)
    src_c = Circle((xc, yc), src_radius, fill=False,
                    edgecolor='magenta', linewidth=1.5,
                    label=f'Source (r={src_radius})')
    ax1.add_patch(src_c)
    ax1.plot(xc, yc, 'x', color='magenta', markersize=10,
             markeredgewidth=2, zorder=5)

    # Background annulus (cyan)
    ax1.add_patch(Circle((xc, yc), bkg_inner, fill=False,
                          edgecolor='cyan', linewidth=1.2,
                          linestyle='--'))
    ax1.add_patch(Circle((xc, yc), bkg_outer, fill=False,
                          edgecolor='cyan', linewidth=1.2,
                          linestyle='--',
                          label=f'Bkg ({bkg_inner}-{bkg_outer})'))

    ax1.legend(loc='upper left', fontsize=8, facecolor='black',
               edgecolor='white', labelcolor='white', framealpha=0.7)
    ax1.set_xlabel('Sky X (pixels)', color='white')
    ax1.set_ylabel('Sky Y (pixels)', color='white')
    ax1.set_facecolor('black')
    ax1.tick_params(colors='white')
    for spine in ax1.spines.values():
        spine.set_edgecolor('white')

    # --- Right: 1D DETX cross-strip profile ---
    ax2 = fig.add_subplot(122)

    if detx_events is not None and len(detx_events) > 0:
        # Find peak
        edges_d = np.arange(detx_events.min() - 1,
                             detx_events.max() + 2, 1)
        hist_d, _ = np.histogram(detx_events, bins=edges_d)
        centers_d = 0.5 * (edges_d[:-1] + edges_d[1:])
        smoothed = gaussian_filter1d(hist_d.astype(float), sigma=1.5)
        detx_center = centers_d[np.argmax(smoothed)]

        # Profile centered on peak
        hw = 100
        x_lo_d = detx_center - hw
        x_hi_d = detx_center + hw
        edges_p = np.arange(x_lo_d, x_hi_d + 1, 1)
        counts_p, _ = np.histogram(detx_events, bins=edges_p)
        x_rel = 0.5 * (edges_p[:-1] + edges_p[1:]) - detx_center

        ax2.step(x_rel, counts_p, where='mid', color='black',
                 linewidth=0.8, label='DETX profile')

        # Source region
        ax2.axvspan(-src_radius, src_radius, alpha=0.15,
                     color='green',
                     label=f'Source ({src_radius*plate_scale:.0f}")')
        # Background regions
        ax2.axvspan(-bkg_outer, -bkg_inner, alpha=0.15, color='cyan',
                     label=f'Bkg ({bkg_inner}-{bkg_outer})')
        ax2.axvspan(bkg_inner, bkg_outer, alpha=0.15, color='cyan')

        # Background level
        wing = max(int(0.15 * len(counts_p)), 3)
        bkg_level = 0.5 * (np.mean(counts_p[:wing]) +
                            np.mean(counts_p[-wing:]))
        ax2.axhline(bkg_level, color='red', linewidth=0.8,
                     linestyle='--', alpha=0.7,
                     label=f'Bkg: {bkg_level:.1f}')

        ax2.set_yscale('log')
        ax2.set_ylim(bottom=max(0.5, bkg_level * 0.1))
        ax2.legend(loc='upper right', fontsize=7)
        ax2.grid(True, alpha=0.3, which='both')
    else:
        ax2.text(0.5, 0.5, 'No DETX data', transform=ax2.transAxes,
                 ha='center', va='center', fontsize=14)

    ax2.set_xlabel('DETX offset from peak (pixels)')
    ax2.set_ylabel('Counts per pixel')

    # Title
    from math import floor, log10
    if count_rate > 0:
        rr = round(count_rate,
                    -int(floor(log10(abs(count_rate)))) + 2)
    else:
        rr = 0
    orbits_str = f'  ({n_gti} orbit{"s" if n_gti > 1 else ""})' \
        if n_gti > 0 else ''
    fig.suptitle(f'{stem}{orbits_str}', fontsize=12,
                 color='white', y=0.99)
    ax1.set_title(
        f'Exp: {exposure:.1f}s  Events: {n_events}  '
        f'Date: {date_mid[:10] if date_mid else "N/A"}',
        fontsize=9, pad=8, loc='left', color='white')
    ax1.text(0.99, 1.03, f'Rate: {rr} ct/s',
             transform=ax1.transAxes, fontsize=10, color='red',
             ha='right', va='bottom', fontweight='bold')

    fig.subplots_adjust(wspace=0.3)
    return fig


# ---------------------------------------------------------------
# Per-orbit sky images in 2x2 grid
# ---------------------------------------------------------------

def split_events_by_gti(info):
    """
    Split event arrays into per-orbit subsets using GTI intervals.

    Returns list of dicts, each with x, y, n_events, duration,
    gti_index.
    """
    orbits = []
    times = info['time_events']

    for i, (tstart, tstop) in enumerate(
            zip(info['gti_starts'], info['gti_stops'])):
        mask = (times >= tstart) & (times <= tstop)
        if np.sum(mask) < 1:
            continue
        orbits.append({
            'x': info['x_events'][mask],
            'y': info['y_events'][mask],
            'n_events': int(np.sum(mask)),
            'duration': tstop - tstart,
            'gti_index': i + 1,
        })

    return orbits


def make_orbit_grid(stem, orbits, xc, yc, src_radius,
                     bkg_inner, bkg_outer):
    """
    Create 2x2 grid pages of per-orbit sky images.
    Returns a list of figures (one per page of 4 orbits).
    """
    figs = []
    n_orbits = len(orbits)
    n_pages = (n_orbits + 3) // 4  # ceil division

    for page in range(n_pages):
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        fig.patch.set_facecolor('black')
        fig.suptitle(f'{stem}  —  per-orbit sky images '
                     f'(page {page+1}/{n_pages})',
                     fontsize=11, color='white', y=0.98)

        for slot in range(4):
            ax = axes[slot // 2][slot % 2]
            orbit_idx = page * 4 + slot

            if orbit_idx >= n_orbits:
                ax.set_visible(False)
                continue

            orb = orbits[orbit_idx]
            x_ev = orb['x']
            y_ev = orb['y']

            # Image bounds centered on source
            margin = max(bkg_outer + 30, 150)
            xl = int(xc - margin)
            xh = int(xc + margin)
            yl = int(yc - margin)
            yh = int(yc + margin)

            m = ((x_ev >= xl) & (x_ev < xh) &
                 (y_ev >= yl) & (y_ev < yh))

            if np.sum(m) > 0:
                img, _, _ = np.histogram2d(
                    y_ev[m], x_ev[m],
                    bins=[np.arange(yl, yh + 1),
                          np.arange(xl, xh + 1)])
                vmax = max(np.max(img), 1)
                ax.imshow(img, origin='lower', cmap='viridis',
                          extent=[xl, xh, yl, yh],
                          interpolation='nearest',
                          norm=PowerNorm(gamma=0.5, vmin=0,
                                         vmax=vmax),
                          aspect='equal')

            # Regions
            ax.add_patch(Circle((xc, yc), src_radius, fill=False,
                                 edgecolor='magenta', linewidth=1.2))
            ax.plot(xc, yc, 'x', color='magenta', markersize=8,
                    markeredgewidth=1.5, zorder=5)
            ax.add_patch(Circle((xc, yc), bkg_inner, fill=False,
                                 edgecolor='cyan', linewidth=0.8,
                                 linestyle='--'))
            ax.add_patch(Circle((xc, yc), bkg_outer, fill=False,
                                 edgecolor='cyan', linewidth=0.8,
                                 linestyle='--'))

            ax.set_title(
                f'Orbit {orb["gti_index"]}: '
                f'{orb["n_events"]} evt, '
                f'{orb["duration"]:.0f}s',
                fontsize=9, color='white')
            ax.set_facecolor('black')
            ax.tick_params(colors='white', labelsize=7)
            for spine in ax.spines.values():
                spine.set_edgecolor('white')

        fig.subplots_adjust(hspace=0.15, wspace=0.15)
        figs.append(fig)

    return figs


# ---------------------------------------------------------------
# Write extraction info text file
# ---------------------------------------------------------------

def write_wt_info(obsid_path, stem, info, xc, yc,
                  src_radius, bkg_inner, bkg_outer,
                  backscal_src, backscal_bkg, plate_scale):
    """Write WT extraction parameters to a text file."""
    output_file = os.path.join(obsid_path, f'{stem}_wt_profile.txt')
    with open(output_file, 'w') as f:
        f.write(f"# WT profile analysis for {stem}\n")
        f.write(f"# Regions in sky coordinates (pixels)\n")
        f.write(f"# Date: {info.get('date_mid', 'N/A')}\n")
        f.write(f"# Exposure: {info['exposure']:.1f} s\n")
        f.write(f"# Events: {info['n_events']}\n")
        f.write(f"# Count rate: {info['count_rate']:.3f} ct/s\n")
        f.write(f"# Plate scale: {plate_scale:.4f} "
                f"arcsec/pixel\n\n")
        f.write(f"source_x = {xc:.2f}\n")
        f.write(f"source_y = {yc:.2f}\n")
        f.write(f"source_radius_pix = {src_radius}\n")
        f.write(f"source_radius_arcsec = "
                f"{src_radius * plate_scale:.1f}\n")
        f.write(f"bkg_inner_pix = {bkg_inner}\n")
        f.write(f"bkg_outer_pix = {bkg_outer}\n\n")
        f.write(f"# BACKSCAL values for WT mode (1D extent):\n")
        f.write(f"# These must be manually set in grppha if\n")
        f.write(f"# XSELECT's auto-values are incorrect.\n")
        f.write(f"backscal_src = {backscal_src}\n")
        f.write(f"backscal_bkg = {backscal_bkg}\n")
    return output_file


# ---------------------------------------------------------------
# Compact summary table
# ---------------------------------------------------------------

def print_compact_table(all_info):
    """Print a one-row-per-observation summary table."""
    src_label = 'Rsrc'
    hdr = (f"  {'OBSID':<14} {'filename':<32} {'Date':<12} "
           f"{'Exp(s)':>8} {'ct/s':>8} {'Events':>7} "
           f"{src_label:>6} {'BkgIn':>6} {'BkgOut':>6} "
           f"{'BkSrc':>6} {'BkBkg':>6}")
    sep = f"  {'-' * (len(hdr) - 2)}"

    print(f"\n{'=' * len(hdr)}")
    print(f"  WT POINTED OBSERVATIONS SUMMARY")
    print(f"{'=' * len(hdr)}")
    print(hdr)
    print(sep)

    total_exp = 0
    total_events = 0

    for info in all_info:
        date_short = info['date_mid'][:10] \
            if info.get('date_mid') else 'N/A'
        print(f"  {info['obsid']:<14} {info['stem']:<32} "
              f"{date_short:<12} "
              f"{info['exposure']:>8.1f} "
              f"{info['count_rate']:>8.2f} "
              f"{info['n_events']:>7} "
              f"{info.get('src_radius', ''):>6} "
              f"{info.get('bkg_inner', ''):>6} "
              f"{info.get('bkg_outer', ''):>6} "
              f"{info.get('backscal_src', ''):>6} "
              f"{info.get('backscal_bkg', ''):>6}")
        total_exp += info['exposure']
        total_events += info['n_events']

    print(sep)
    print(f"  {'TOTAL':<14} {'':32} {'':12} "
          f"{total_exp:>8.1f} {'':>8} {total_events:>7}")
    print(f"{'=' * len(hdr)}")


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Summary and viewer for Swift XRT WT-mode '
                    'pointed event files.')
    parser.add_argument('--ra', type=float, required=True,
                        help='Source RA in degrees')
    parser.add_argument('--dec', type=float, required=True,
                        help='Source Dec in degrees')
    parser.add_argument('--srcrad', type=int, default=20,
                        help='Source extraction radius in pixels '
                             '(default: 20)')
    parser.add_argument('--bkginner', type=int, default=80,
                        help='Background annulus inner radius in '
                             'pixels (default: 80)')
    parser.add_argument('--bkgouter', type=int, default=120,
                        help='Background annulus outer radius in '
                             'pixels (default: 120)')
    parser.add_argument('--compact', action='store_true',
                        help='Print compact summary only, no plots')
    parser.add_argument('--expgt', type=float, default=20.0,
                        help='Minimum exposure time in seconds '
                             '(default: 20)')
    parser.add_argument('--pdf', type=str,
                        default='wt_profiles.pdf',
                        help='Output PDF filename')
    parser.add_argument('--nmax', type=int, default=None,
                        help='Process only first N observations')
    args = parser.parse_args()

    # Validate background annulus is symmetric about ~100 pix
    bkg_mid = (args.bkginner + args.bkgouter) / 2.0
    if abs(bkg_mid - 100) > 20:
        print(f"NOTE: Background annulus midpoint is {bkg_mid:.0f} "
              f"pix (ideal is ~100 for WT window half-width).")

    # Compute BACKSCAL values
    backscal_src, backscal_bkg = compute_backscal(
        args.srcrad, args.bkginner, args.bkgouter)

    # Find WT pointed files
    wt_entries = find_wt_files()
    if not wt_entries:
        print("No WT pointed event files found.")
        sys.exit(1)

    print(f"Found {len(wt_entries)} WT pointed event file(s).")
    print(f"Source position: RA={args.ra:.6f}, Dec={args.dec:.6f}")
    print(f"Source region: circle r={args.srcrad} pix")
    print(f"Background: annulus {args.bkginner}-{args.bkgouter} pix")
    print(f"BACKSCAL: src={backscal_src}, bkg={backscal_bkg}")
    print(f"Min exposure: {args.expgt} s")

    if args.nmax is not None:
        wt_entries = wt_entries[:args.nmax]
        print(f"  (limited to first {args.nmax})")

    all_info = []
    all_figs = []
    n_total = len(wt_entries)

    for idx, entry in enumerate(wt_entries, 1):
        obsid = entry['obsid']
        stem = entry['stem']
        filepath = entry['filepath']
        obsid_path = os.path.join(BASE_DIR, obsid)

        print(f"\n  [{idx}/{n_total}] {obsid} / {stem}")

        try:
            info = get_wt_metadata(filepath)
        except Exception as e:
            print(f"    ERROR reading file: {e}")
            continue

        info['obsid'] = obsid
        info['stem'] = stem

        print(f"    Exposure: {info['exposure']:.1f}s  "
              f"Events: {info['n_events']}  "
              f"Rate: {info['count_rate']:.2f} ct/s  "
              f"Date: {info.get('date_mid', 'N/A')}")

        if info['n_events'] < 10:
            print(f"    SKIPPED: too few events.")
            continue

        if info['exposure'] < args.expgt:
            print(f"    SKIPPED: exposure {info['exposure']:.1f}s "
                  f"< {args.expgt}s minimum.")
            continue

        plate_scale = info['plate_scale']

        # Find source position in sky coordinates
        xc, yc = find_sky_position(
            info['x_events'], info['y_events'],
            args.ra, args.dec, info)
        print(f"    Source (sky): X={xc:.1f}, Y={yc:.1f}")

        n_gti = info.get('n_gti', 1)
        if n_gti > 1:
            print(f"    Multi-orbit: {n_gti} GTI segments")

        # Store extraction parameters
        info['src_radius'] = args.srcrad
        info['bkg_inner'] = args.bkginner
        info['bkg_outer'] = args.bkgouter
        info['backscal_src'] = backscal_src
        info['backscal_bkg'] = backscal_bkg

        # Count events in source and background regions
        dist = np.sqrt((info['x_events'] - xc)**2 +
                        (info['y_events'] - yc)**2)
        src_cts = np.sum(dist <= args.srcrad)
        bkg_cts = np.sum((dist >= args.bkginner) &
                          (dist <= args.bkgouter))
        print(f"    Source counts (r={args.srcrad}): {src_cts}")
        print(f"    Bkg counts ({args.bkginner}-{args.bkgouter}): "
              f"{bkg_cts}")

        # Write info file
        txt_file = write_wt_info(
            obsid_path, stem, info, xc, yc,
            args.srcrad, args.bkginner, args.bkgouter,
            backscal_src, backscal_bkg, plate_scale)
        print(f"    Saved: {os.path.basename(txt_file)}")

        all_info.append(info)

        # Generate plots
        if not args.compact:
            # Combined sky image + DETX profile
            fig_main = make_combined_plot(
                stem, info['x_events'], info['y_events'],
                info['detx_events'], xc, yc,
                args.srcrad, args.bkginner, args.bkgouter,
                plate_scale, info['exposure'],
                info['count_rate'], info['n_events'],
                info.get('date_mid'), n_gti)

            png_main = os.path.join(obsid_path,
                                     f'{stem}_wt_combined.png')
            fig_main.savefig(png_main, dpi=150,
                             bbox_inches='tight', facecolor='black')
            print(f"    Saved: {os.path.basename(png_main)}")
            all_figs.append(fig_main)

            # Per-orbit sky images for multi-orbit observations
            if n_gti > 1:
                orbits = split_events_by_gti(info)
                if len(orbits) > 1:
                    orbit_figs = make_orbit_grid(
                        stem, orbits, xc, yc,
                        args.srcrad, args.bkginner,
                        args.bkgouter)
                    for ofig in orbit_figs:
                        all_figs.append(ofig)
                    print(f"    Created {len(orbit_figs)} "
                          f"orbit grid page(s) "
                          f"({len(orbits)} orbits)")

    # Print compact table
    if all_info:
        print_compact_table(all_info)

    # Collate into PDF
    if all_figs and not args.compact:
        pdf_path = os.path.join(BASE_DIR, args.pdf)
        print(f"\nCollating {len(all_figs)} figures into {pdf_path}")
        with PdfPages(pdf_path) as pdf:
            for fig in all_figs:
                pdf.savefig(fig, dpi=150, bbox_inches='tight',
                            facecolor=fig.get_facecolor())
                plt.close(fig)
        print(f"PDF saved: {pdf_path}")

    print(f"\nDone. {len(all_info)} WT observations processed.")


if __name__ == '__main__':
    main()

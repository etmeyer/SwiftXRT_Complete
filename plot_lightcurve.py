#!/usr/bin/env python3
"""
plot_lightcurve.py

Remake the νFν / Γ light curve plot from an existing
fit_results.txt file. Allows customization of axis limits
and vertical reference lines without re-running the fits.

Usage:
    python plot_lightcurve.py
    python plot_lightcurve.py --ylim_flux 1e-12 1e-10
    python plot_lightcurve.py --ylim_gamma 1.0 2.5
    python plot_lightcurve.py --xlim 2005 2024
    python plot_lightcurve.py --times 2007.75 2015.5 2023.3
    python plot_lightcurve.py --times 2007.75 2015.5 --tlabels "Epoch A" "Epoch B"
"""

import os
import sys
import re
import argparse
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print("ERROR: matplotlib is required.")
    sys.exit(1)


def read_fit_results(filepath):
    """
    Parse fit_results.txt into a list of dicts.
    Handles both the header row and data rows, skipping
    comment lines and separator lines.
    """
    results = []
    col_names = None

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Skip separator lines
            if set(line.replace(' ', '')) == {'-'}:
                continue
            # First non-comment, non-separator line is the header
            if col_names is None:
                col_names = line.split()
                continue
            # Data rows
            parts = line.split()
            if len(parts) < len(col_names):
                continue
            row = {}
            for i, name in enumerate(col_names):
                val = parts[i]
                row[name] = val
            results.append(row)

    return results


def parse_float(val):
    """Convert string to float, returning None for N/A or (frozen)."""
    if val is None or val == 'N/A' or val.startswith('('):
        return None
    try:
        return float(val)
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(
        description='Plot νFν and Γ light curve from fit results.')
    parser.add_argument('--input', type=str, default='fit_results.txt',
                        help='Input results table '
                             '(default: fit_results.txt)')
    parser.add_argument('--output', type=str,
                        default='flux_lightcurve.pdf',
                        help='Output PDF (default: flux_lightcurve.pdf)')
    parser.add_argument('--xlim', type=float, nargs=2, default=None,
                        metavar=('XMIN', 'XMAX'),
                        help='X-axis limits in decimal years')
    parser.add_argument('--ylim_flux', type=float, nargs=2, default=None,
                        metavar=('YMIN', 'YMAX'),
                        help='Upper panel (νFν) y-axis limits')
    parser.add_argument('--ylim_gamma', type=float, nargs=2, default=None,
                        metavar=('YMIN', 'YMAX'),
                        help='Lower panel (Γ) y-axis limits')
    parser.add_argument('--times', type=float, nargs='+', default=None,
                        metavar='YEAR',
                        help='Vertical reference lines at these '
                             'decimal years')
    parser.add_argument('--tlabels', type=str, nargs='+', default=None,
                        metavar='LABEL',
                        help='Labels for --times lines (must match '
                             'number of --times values)')
    parser.add_argument('--tcolor', type=str, default='green',
                        help='Color for reference lines '
                             '(default: green)')
    parser.add_argument('--title', type=str, default=None,
                        help='Plot title (default: auto)')
    parser.add_argument('--figsize', type=float, nargs=2,
                        default=[14, 7], metavar=('W', 'H'),
                        help='Figure size in inches (default: 14 7)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Output DPI (default: 150)')
    args = parser.parse_args()

    # Validate --tlabels
    if args.tlabels and args.times:
        if len(args.tlabels) != len(args.times):
            print("ERROR: --tlabels must have same number of "
                  "entries as --times.")
            sys.exit(1)

    # Read results
    if not os.path.exists(args.input):
        print(f"ERROR: {args.input} not found.")
        sys.exit(1)

    rows = read_fit_results(args.input)
    if not rows:
        print("ERROR: no data rows found in results file.")
        sys.exit(1)

    print(f"Read {len(rows)} entries from {args.input}")

    # Extract arrays
    valid = []
    for r in rows:
        mjd = parse_float(r.get('MJD'))
        f1 = parse_float(r.get('flux_1keV'))
        gam = parse_float(r.get('gamma'))
        if mjd is not None and f1 is not None and gam is not None:
            valid.append({
                'mjd': mjd,
                'flux_1keV': f1,
                'flux_1keV_err': parse_float(r.get('f1_err')),
                'gamma': gam,
                'gamma_err': parse_float(r.get('gamma_err')),
                'gamma_frozen': r.get('gamma_err', '') == '(frozen)',
                'mode': r.get('mode', 'PC'),
            })

    if not valid:
        print("ERROR: no valid data points.")
        sys.exit(1)

    print(f"  {len(valid)} valid data points")

    mjd = np.array([v['mjd'] for v in valid])
    flux = np.array([v['flux_1keV'] for v in valid])
    flux_err = np.array([
        v['flux_1keV_err'] if v['flux_1keV_err'] else 0
        for v in valid])
    gamma = np.array([v['gamma'] for v in valid])
    gamma_err = np.array([
        v['gamma_err'] if v['gamma_err'] and not v['gamma_frozen']
        else 0 for v in valid])
    frozen = np.array([v['gamma_frozen'] for v in valid])
    modes = np.array([v['mode'] for v in valid])

    # Sort by time
    order = np.argsort(mjd)
    mjd = mjd[order]; flux = flux[order]; flux_err = flux_err[order]
    gamma = gamma[order]; gamma_err = gamma_err[order]
    frozen = frozen[order]; modes = modes[order]

    # Convert to decimal years and νFν
    decimal_year = 2000.0 + (mjd - 51544.0) / 365.25
    nu_1keV = 2.418e17
    vfv = flux * nu_1keV
    vfv_err = flux_err * nu_1keV

    # Mode masks
    pc = modes == 'PC'
    wt = modes == 'WT'

    # --- Create plot ---
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=args.figsize, sharex=True,
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.05})

    # Upper: νFν
    if np.any(pc):
        ax1.errorbar(decimal_year[pc], vfv[pc], yerr=vfv_err[pc],
                     fmt='o', ms=5, color='navy',
                     ecolor='cornflowerblue',
                     elinewidth=1, capsize=3, label='PC')
    if np.any(wt):
        ax1.errorbar(decimal_year[wt], vfv[wt], yerr=vfv_err[wt],
                     fmt='D', ms=5, color='darkorange',
                     ecolor='sandybrown',
                     elinewidth=1, capsize=3, label='WT')

    ax1.set_ylabel(r'$\nu F_{\nu}$(1 keV)  [erg cm$^{-2}$ s$^{-1}$]',
                    fontsize=11)
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, which='both')
    title = args.title or 'Swift XRT Light Curve'
    ax1.set_title(title, fontsize=13)
    ax1.legend(loc='best', fontsize=9)

    if args.ylim_flux:
        ax1.set_ylim(args.ylim_flux[0], args.ylim_flux[1])

    # Lower: Γ
    # PC
    m = pc & ~frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=gamma_err[m],
                     fmt='o', ms=5, color='darkred',
                     ecolor='salmon', elinewidth=1, capsize=3,
                     label=r'PC free $\Gamma$')
    m = pc & frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=0,
                     fmt='o', ms=5, color='silver',
                     elinewidth=1, capsize=3,
                     label=r'PC frozen $\Gamma$')
    # WT
    m = wt & ~frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=gamma_err[m],
                     fmt='D', ms=5, color='darkorange',
                     ecolor='sandybrown', elinewidth=1, capsize=3,
                     label=r'WT free $\Gamma$')
    m = wt & frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=0,
                     fmt='D', ms=5, color='wheat',
                     elinewidth=1, capsize=3,
                     label=r'WT frozen $\Gamma$')

    ax2.set_ylabel(r'$\Gamma$', fontsize=12)
    ax2.set_xlabel('Year', fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=8, ncol=2)

    if args.ylim_gamma:
        ax2.set_ylim(args.ylim_gamma[0], args.ylim_gamma[1])

    if args.xlim:
        ax1.set_xlim(args.xlim[0], args.xlim[1])

    # Set x-axis major ticks to integer years
    from matplotlib.ticker import MultipleLocator, AutoMinorLocator
    ax2.xaxis.set_major_locator(MultipleLocator(1))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))

    # Vertical reference lines
    if args.times:
        for i, t in enumerate(args.times):
            ax1.axvline(t, color=args.tcolor, linewidth=1.2,
                        linestyle='--', alpha=0.7, zorder=0)
            ax2.axvline(t, color=args.tcolor, linewidth=1.2,
                        linestyle='--', alpha=0.7, zorder=0)
            if args.tlabels:
                # Place label above the upper panel so it doesn't
                # get clipped by the axes
                ax1.annotate(args.tlabels[i],
                             xy=(t, 1.0),
                             xycoords=('data', 'axes fraction'),
                             xytext=(3, 4), textcoords='offset points',
                             fontsize=8, color=args.tcolor,
                             ha='left', va='bottom',
                             rotation=90, clip_on=False)

    fig.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)

    n_pc = np.sum(pc)
    n_wt = np.sum(wt)
    print(f"Saved: {args.output} ({n_pc} PC, {n_wt} WT points)")


if __name__ == '__main__':
    main()

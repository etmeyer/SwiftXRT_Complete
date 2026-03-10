#!/usr/bin/env python3
"""
swift_xrt_king_profile.py

For each OBSID directory, find PC-mode cleaned level-2 event files,
localize a known point source, fit a King profile to the radial surface
brightness profile in a specified annular region, and produce diagnostic
plots showing the full radial profile with the King fit overlaid.

Usage:
    python swift_xrt_king_profile.py --ra 187.2779 --dec 2.0524 \
        --rmin 20 --rmax 60 --rbin 2 --maxplot 80

Required arguments:
    --ra        Source RA in degrees
    --dec       Source Dec in degrees

Optional arguments:
    --rmin      Inner radius of fitting annulus in arcsec (default: 20)
    --rmax      Outer radius of fitting annulus in arcsec (default: 60)
    --rbin      Radial bin width in arcsec (default: 2)
    --maxplot   Maximum radius to plot in arcsec (default: rmax + 20)
    --centroid  Centroid search radius in arcsec (default: 15)
    --rc        King profile core radius, fixed (default: 5.8)
    --beta      King profile beta slope, fixed (default: 1.55)
    --sigma     Sigma threshold, 2 consecutive bins (default: 3.0)
    --sigma2    Sigma threshold, single bin (default: 4.0)
    --pdf       Output PDF filename (default: king_profiles.pdf)

Optional override file:
    If a file named 'pileup_overrides.txt' exists in the working
    directory, it will be read. Format: one entry per line with the
    event file stem in column 1 and the override radius in arcsec
    in column 2. These are plotted as purple dashed lines.

Requirements:
    astropy, scipy, matplotlib, numpy
"""

import os
import sys
import glob
import re
import argparse
import numpy as np

try:
    from astropy.io import fits
    from astropy.wcs import WCS
except ImportError:
    print("ERROR: astropy is required. Install with: pip install astropy")
    sys.exit(1)

try:
    from scipy.optimize import curve_fit
except ImportError:
    print("ERROR: scipy is required. Install with: pip install scipy")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("ERROR: matplotlib is required. Install with: pip install matplotlib")
    sys.exit(1)


# ---------------------------------------------------------------
# Window mode identification
# ---------------------------------------------------------------

def get_window_info(filename):
    """
    Extract CCD window configuration from XRT event filename.

    Swift XRT PC mode uses different window sizes that determine the
    frame readout time, which directly affects the pile-up threshold.

    Returns:
        dict with 'window' (e.g. 'w3'), 'size' (e.g. '300x300'),
        'frame_time' (approx seconds), and 'label' (display string)
    """
    basename = os.path.basename(filename).lower()
    match = re.search(r'xpc(w\d)', basename)
    if not match:
        return {'window': '?', 'size': '?', 'frame_time': None,
                'label': 'Window: unknown'}

    window = match.group(1)
    # Approximate window sizes and frame times for PC mode
    window_info = {
        'w1': {'size': '100x100', 'frame_time': 0.28},
        'w2': {'size': '200x200', 'frame_time': 0.89},
        'w3': {'size': '300x300', 'frame_time': 1.77},
        'w4': {'size': '480x480', 'frame_time': 2.51},
    }
    info = window_info.get(window, {'size': '?', 'frame_time': None})
    ft = info['frame_time']
    label = f"Window: {window} ({info['size']}px"
    if ft is not None:
        label += f", {ft:.2f}s/frame"
    label += ")"

    return {'window': window, 'size': info['size'],
            'frame_time': ft, 'label': label}


# ---------------------------------------------------------------
# King profile model
# ---------------------------------------------------------------

def king_profile(r, S0, rc, beta, bkg):
    """
    King (beta) profile for X-ray surface brightness.

    S(r) = S0 * (1 + (r/rc)^2)^(-beta) + bkg

    Parameters:
        r     : radius (arcsec)
        S0    : central surface brightness (cts/s/arcmin^2)
        rc    : core radius (arcsec)
        beta  : slope parameter
        bkg   : constant background level (cts/s/arcmin^2)
    """
    return S0 * (1.0 + (r / rc) ** 2) ** (-beta) + bkg


# ---------------------------------------------------------------
# Source centroiding
# ---------------------------------------------------------------

def refine_centroid(x_events, y_events, x0, y0, search_radius_pix,
                    iterations=3):
    """
    Iteratively refine source centroid using mean of events within
    a shrinking aperture.

    Parameters:
        x_events, y_events : event pixel coordinates
        x0, y0             : initial guess (pixels)
        search_radius_pix  : initial search radius (pixels)
        iterations         : number of refinement steps

    Returns:
        (x_cen, y_cen) : refined centroid in pixels
    """
    xc, yc = x0, y0
    for i in range(iterations):
        radius = search_radius_pix / (1.0 + i * 0.3)
        dist = np.sqrt((x_events - xc) ** 2 + (y_events - yc) ** 2)
        mask = dist < radius
        if np.sum(mask) < 5:
            break
        xc = np.mean(x_events[mask])
        yc = np.mean(y_events[mask])
    return xc, yc


# ---------------------------------------------------------------
# Radial profile extraction
# ---------------------------------------------------------------

def extract_radial_profile(x_events, y_events, xc, yc, plate_scale,
                           rbin_arcsec, rmax_arcsec, exposure):
    """
    Compute azimuthally-averaged surface brightness profile.

    Returns:
        r_mid     : bin centers (arcsec)
        sb        : surface brightness (cts/s/arcmin^2)
        sb_err    : 1-sigma Poisson error on sb
        counts    : raw counts per bin
    """
    # Radial distance in arcsec
    dist_pix = np.sqrt((x_events - xc) ** 2 + (y_events - yc) ** 2)
    dist_arcsec = dist_pix * plate_scale

    # Define bins
    bin_edges = np.arange(0, rmax_arcsec + rbin_arcsec, rbin_arcsec)
    n_bins = len(bin_edges) - 1

    r_mid = np.zeros(n_bins)
    sb = np.zeros(n_bins)
    sb_err = np.zeros(n_bins)
    counts = np.zeros(n_bins, dtype=int)

    for i in range(n_bins):
        r_in = bin_edges[i]
        r_out = bin_edges[i + 1]
        r_mid[i] = 0.5 * (r_in + r_out)

        # Annular area in arcmin^2 (divide arcsec^2 by 3600)
        area_arcsec2 = np.pi * (r_out ** 2 - r_in ** 2)
        area_arcmin2 = area_arcsec2 / 3600.0

        # Count events in this annular bin
        mask = (dist_arcsec >= r_in) & (dist_arcsec < r_out)
        cts = np.sum(mask)
        counts[i] = cts

        # Surface brightness and Poisson error in cts/s/arcmin^2
        if exposure > 0 and area_arcmin2 > 0:
            sb[i] = cts / (exposure * area_arcmin2)
            sb_err[i] = np.sqrt(max(cts, 1)) / (exposure * area_arcmin2)
        else:
            sb[i] = 0
            sb_err[i] = 0

    return r_mid, sb, sb_err, counts


# ---------------------------------------------------------------
# Fitting routine
# ---------------------------------------------------------------

def fit_king_profile(r_mid, sb, sb_err, rmin_fit, rmax_fit, rc_fixed,
                     beta_fixed):
    """
    Fit King profile to data within the specified radial range,
    with rc and beta fixed to known Swift PSF values. Only S0
    (normalization) and bkg (background) are free parameters.

    Returns:
        popt  : best-fit parameters (S0, rc, beta, bkg)
        pcov  : covariance matrix (2x2 for S0 and bkg)
        mask  : boolean mask of bins used in fit
    """
    mask = (r_mid >= rmin_fit) & (r_mid <= rmax_fit) & (sb_err > 0)

    if np.sum(mask) < 3:
        print("    WARNING: fewer than 3 bins in fitting range.")
        return None, None, mask

    r_fit = r_mid[mask]
    sb_fit = sb[mask]
    err_fit = sb_err[mask]

    # Wrapper with rc and beta fixed
    def king_fixed(r, S0, bkg):
        return king_profile(r, S0, rc_fixed, beta_fixed, bkg)

    # Initial guesses
    S0_guess = np.max(sb_fit) * 2
    bkg_guess = np.min(sb_fit)

    p0 = [S0_guess, bkg_guess]
    bounds = ([0, 0], [np.inf, np.inf])

    try:
        popt_free, pcov_free = curve_fit(king_fixed, r_fit, sb_fit,
                                          p0=p0, sigma=err_fit,
                                          absolute_sigma=True,
                                          bounds=bounds, maxfev=10000)
        # Pack into full 4-parameter format for consistency
        popt = np.array([popt_free[0], rc_fixed, beta_fixed,
                         popt_free[1]])
        # Build a 4x4 covariance with zeros for fixed params
        pcov = np.zeros((4, 4))
        pcov[0, 0] = pcov_free[0, 0]  # S0 variance
        pcov[0, 3] = pcov_free[0, 1]  # S0-bkg covariance
        pcov[3, 0] = pcov_free[1, 0]
        pcov[3, 3] = pcov_free[1, 1]  # bkg variance

        return popt, pcov, mask
    except (RuntimeError, ValueError) as e:
        print(f"    WARNING: King profile fit failed: {e}")
        return None, None, mask


# ---------------------------------------------------------------
# Pile-up radius detection
# ---------------------------------------------------------------

def detect_pileup_radius(r_mid, sb, sb_err, popt, sigma1, sigma2,
                         rbin):
    """
    Detect the radius at which pile-up affects the radial profile
    using a two-tier sigma threshold.

    Starting from the bin closest to 30 arcsec, walk inward.
    Tier 1 (sigma1): 2+ consecutive bins exceeding sigma1 triggers.
    Tier 2 (sigma2): a single bin exceeding sigma2 triggers.
    Whichever is encountered first (outermost) sets the radius.

    Parameters:
        r_mid   : bin centers (arcsec)
        sb      : surface brightness
        sb_err  : errors
        popt    : King profile best-fit parameters
        sigma1  : lower threshold requiring 2 consecutive (e.g. 3.0)
        sigma2  : higher threshold requiring only 1 (e.g. 4.0)
        rbin    : bin width in arcsec

    Returns:
        pileup_radius    : recommended inner exclusion radius (arcsec),
                           or None if no pile-up detected
        flag_level       : integer array (0=normal, 1=sigma1, 2=sigma2)
        residuals_sigma  : residuals in sigma units for all bins
    """
    flag_level = np.zeros(len(r_mid), dtype=int)

    if popt is None:
        return None, flag_level, np.zeros(len(r_mid))

    # Compute residuals in sigma for all bins with valid data
    model = king_profile(r_mid, *popt)
    residuals_sigma = np.zeros(len(r_mid))
    valid = sb_err > 0
    residuals_sigma[valid] = (sb[valid] - model[valid]) / sb_err[valid]

    # Assign flag levels
    abs_resid = np.abs(residuals_sigma)
    flag_level[abs_resid >= sigma1] = 1
    flag_level[abs_resid >= sigma2] = 2

    # Find the bin closest to 30 arcsec as starting point
    start_idx = np.argmin(np.abs(r_mid - 30.0))

    # Walk inward, checking both tiers
    # Track: tier1 needs 2 consecutive, tier2 needs 1
    pileup_radius = None
    consecutive_t1 = 0

    for i in range(start_idx, -1, -1):
        if not valid[i]:
            consecutive_t1 = 0
            continue

        # Check tier 2 first (single bin, higher threshold)
        if flag_level[i] >= 2:
            candidate = r_mid[i] + rbin / 2.0
            if pileup_radius is None or candidate > pileup_radius:
                pileup_radius = candidate
            # Also reset consecutive counter since this bin counts
            # for tier 1 as well
            consecutive_t1 = 0
            break  # outermost trigger found, stop

        # Check tier 1 (consecutive bins, lower threshold)
        if flag_level[i] >= 1:
            consecutive_t1 += 1
            if consecutive_t1 >= 2:
                outer_bin_idx = i + consecutive_t1 - 1
                candidate = r_mid[outer_bin_idx] + rbin / 2.0
                if pileup_radius is None or candidate > pileup_radius:
                    pileup_radius = candidate
                break  # outermost trigger found, stop
        else:
            consecutive_t1 = 0

    return pileup_radius, flag_level, residuals_sigma


# ---------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------

def make_profile_plot(obsid, evt_stem, r_mid, sb, sb_err, popt, pcov,
                      mask_fit, rmin_fit, rmax_fit, rmax_plot, exposure,
                      n_events, window_info, sigma1, sigma2,
                      pileup_radius, override_radius, flag_level,
                      output_png):
    """
    Create the radial profile plot with King fit overlay.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 10),
                                    gridspec_kw={'height_ratios': [3, 1],
                                                 'hspace': 0.08},
                                    sharex=True)

    # Mask for data to plot (non-zero surface brightness)
    plot_mask = (r_mid <= rmax_plot) & (sb > 0)

    # --- Main panel: data + fit ---
    ax1.errorbar(r_mid[plot_mask], sb[plot_mask], yerr=sb_err[plot_mask],
                 fmt='o', ms=4, color='black', ecolor='gray',
                 elinewidth=0.8, capsize=2, label='Data', zorder=2)

    # Shade the fitting region
    ax1.axvspan(rmin_fit, rmax_fit, alpha=0.08, color='blue',
                label=f'Fit region ({rmin_fit:.0f}-{rmax_fit:.0f}")')

    # Overplot the King profile
    if popt is not None:
        r_model = np.linspace(0.1, rmax_plot, 500)
        sb_model = king_profile(r_model, *popt)
        ax1.plot(r_model, sb_model, '-', color='red', linewidth=1.8,
                 label='King profile fit', zorder=3)

        # Parameter text
        perr = np.sqrt(np.diag(pcov)) if pcov is not None else [0, 0, 0, 0]
        param_text = (
            f"$S_0$ = {popt[0]:.2e} $\\pm$ {perr[0]:.1e}\n"
            f"$r_c$ = {popt[1]:.2f}\" (fixed)\n"
            f"$\\beta$ = {popt[2]:.3f} (fixed)\n"
            f"bkg = {popt[3]:.2e} $\\pm$ {perr[3]:.1e}"
        )
        ax1.text(0.97, 0.97, param_text, transform=ax1.transAxes,
                 fontsize=9, verticalalignment='top',
                 horizontalalignment='right',
                 bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat',
                           alpha=0.8),
                 fontfamily='monospace')

    ax1.set_yscale('log')
    ax1.set_ylabel('Surface Brightness (cts/s/arcmin$^2$)')
    count_rate = n_events / exposure if exposure > 0 else 0
    # Round to 3 significant figures
    if count_rate > 0:
        from math import floor, log10
        rounded_rate = round(count_rate,
                             -int(floor(log10(abs(count_rate)))) + 2)
    else:
        rounded_rate = 0
    # Title: filename on top, details below with rate in red overlay
    fig.suptitle(evt_stem, fontsize=12, y=0.98)
    ax1.set_title(f'Exposure: {exposure:.1f}s    Events: {n_events}'
                  f'    {window_info["label"]}',
                  fontsize=10, pad=12, loc='left')
    ax1.text(0.99, 1.025, f'Rate: {rounded_rate} ct/s',
             transform=ax1.transAxes, fontsize=11, color='red',
             ha='right', va='bottom', fontweight='bold')
    ax1.legend(loc='right', bbox_to_anchor=(1.0, 0.62), fontsize=9)
    ax1.set_xlim(0, rmax_plot)
    ax1.grid(True, alpha=0.3, which='both')

    # --- Residual panel ---
    if popt is not None:
        residual_mask = plot_mask & (sb_err > 0)
        r_res = r_mid[residual_mask]
        sb_res = sb[residual_mask]
        err_res = sb_err[residual_mask]
        model_at_data = king_profile(r_res, *popt)
        residuals = (sb_res - model_at_data) / err_res
        flags_res = flag_level[residual_mask]

        # Four categories: normal in-fit, normal outside-fit,
        # tier1 (orange), tier2 (red)
        in_fit = (r_res >= rmin_fit) & (r_res <= rmax_fit)
        is_t1 = flags_res == 1  # between sigma1 and sigma2
        is_t2 = flags_res == 2  # above sigma2
        is_normal = flags_res == 0

        # Plot normal points first
        norm_fit = in_fit & is_normal
        norm_out = ~in_fit & is_normal
        if np.any(norm_fit):
            ax2.errorbar(r_res[norm_fit], residuals[norm_fit],
                         yerr=1, fmt='o', ms=4, color='blue',
                         ecolor='lightblue', elinewidth=0.8, capsize=2,
                         label='In fit region')
        if np.any(norm_out):
            ax2.errorbar(r_res[norm_out], residuals[norm_out],
                         yerr=1, fmt='s', ms=4, color='gray',
                         ecolor='lightgray', elinewidth=0.8, capsize=2,
                         label='Outside fit region')
        # Tier 1: orange diamonds (sigma1 <= |resid| < sigma2)
        if np.any(is_t1):
            ax2.errorbar(r_res[is_t1], residuals[is_t1],
                         yerr=1, fmt='D', ms=5, color='darkorange',
                         ecolor='moccasin', elinewidth=0.8, capsize=2,
                         label=f'>{sigma1:.1f}$\\sigma$', zorder=4)
        # Tier 2: red diamonds (|resid| >= sigma2)
        if np.any(is_t2):
            ax2.errorbar(r_res[is_t2], residuals[is_t2],
                         yerr=1, fmt='D', ms=5, color='red',
                         ecolor='lightsalmon', elinewidth=0.8, capsize=2,
                         label=f'>{sigma2:.1f}$\\sigma$', zorder=5)

        ax2.axhline(0, color='red', linewidth=1)
        ax2.axvspan(rmin_fit, rmax_fit, alpha=0.08, color='blue')
        ax2.set_ylabel('Residual ($\\sigma$)')
        ax2.legend(loc='lower right', fontsize=8)

        # Draw pile-up radius line in both panels
        if pileup_radius is not None:
            for ax in (ax1, ax2):
                ax.axvline(pileup_radius, color='red', linestyle='--',
                           linewidth=1.5, alpha=0.8)
            ax1.text(pileup_radius + 0.5, 0.5,
                     f'Pile-up: {pileup_radius:.1f}"',
                     transform=ax1.get_xaxis_transform(),
                     fontsize=9, color='red', va='center',
                     fontweight='bold')

        # Draw override radius in purple if present
        if override_radius is not None:
            for ax in (ax1, ax2):
                ax.axvline(override_radius, color='purple',
                           linestyle='--', linewidth=1.5, alpha=0.8)
            ax1.text(override_radius + 0.5, 0.35,
                     f'Override: {override_radius:.1f}"',
                     transform=ax1.get_xaxis_transform(),
                     fontsize=9, color='purple', va='center',
                     fontweight='bold')
    else:
        ax2.text(0.5, 0.5, 'Fit failed', transform=ax2.transAxes,
                 ha='center', va='center', fontsize=14, color='red')

    ax2.set_xlabel('Radius (arcsec)')
    ax2.grid(True, alpha=0.3)

    fig.savefig(output_png, dpi=150, bbox_inches='tight')

    return fig


# ---------------------------------------------------------------
# Process a single event file
# ---------------------------------------------------------------

def process_event_file(filepath, ra_src, dec_src, rmin_fit, rmax_fit,
                       rbin, rmax_plot, centroid_radius, rc_fixed,
                       beta_fixed, sigma1, sigma2, obsid, output_dir,
                       overrides):
    """
    Process one PC mode event file: centroid, extract profile, fit, plot.

    Returns:
        (output_png, fig) if successful, (None, None) otherwise.
    """
    filename = os.path.basename(filepath)
    print(f"  Processing: {filename}")

    try:
        with fits.open(filepath) as hdul:
            # Find the EVENTS extension
            evt_ext = None
            for i, ext in enumerate(hdul):
                if ext.name == 'EVENTS':
                    evt_ext = i
                    break
            if evt_ext is None:
                evt_ext = 1  # fallback

            header = hdul[evt_ext].header
            data = hdul[evt_ext].data

            if data is None or len(data) == 0:
                print("    WARNING: No events in file.")
                return None, None

            exposure = header.get('EXPOSURE', 0)
            if exposure <= 0:
                print("    WARNING: Zero or missing exposure time.")
                return None, None

            # Get event coordinates
            x_events = data['X'].astype(float)
            y_events = data['Y'].astype(float)

            # Find column indices for X and Y to get their WCS keywords
            # FITS table WCS uses TCRVL{n}, TCRPX{n}, TCDLT{n}
            x_col = None
            y_col = None
            tfields = header.get('TFIELDS', 0)
            for ci in range(1, tfields + 1):
                ttype = header.get(f'TTYPE{ci}', '').strip().upper()
                if ttype == 'X':
                    x_col = ci
                elif ttype == 'Y':
                    y_col = ci

            # Convert RA/Dec to pixel using table column WCS
            if x_col and y_col:
                tcrvl_x = header.get(f'TCRVL{x_col}', 0)
                tcrpx_x = header.get(f'TCRPX{x_col}', 500)
                tcdlt_x = header.get(f'TCDLT{x_col}', -0.0006548089)

                tcrvl_y = header.get(f'TCRVL{y_col}', 0)
                tcrpx_y = header.get(f'TCRPX{y_col}', 500)
                tcdlt_y = header.get(f'TCDLT{y_col}', 0.0006548089)

                # Simple linear WCS (valid for small XRT FOV)
                x0 = tcrpx_x + (ra_src - tcrvl_x) / tcdlt_x
                y0 = tcrpx_y + (dec_src - tcrvl_y) / tcdlt_y

                plate_scale = abs(tcdlt_x) * 3600  # deg/pixel -> arcsec/pixel
            else:
                # Fallback: try standard image WCS via astropy
                try:
                    wcs = WCS(header, naxis=2)
                    x0, y0 = wcs.all_world2pix(ra_src, dec_src, 1)
                    cdelt = abs(header.get('CDELT1', 0.0006548089))
                    plate_scale = cdelt * 3600
                except Exception:
                    print("    WARNING: Could not determine WCS.")
                    return None, None

            print(f"    Initial position (pix): ({x0:.2f}, {y0:.2f})")
            print(f"    Plate scale: {plate_scale:.3f} arcsec/pixel")

            # Refine centroid
            search_pix = centroid_radius / plate_scale
            xc, yc = refine_centroid(x_events, y_events, x0, y0,
                                      search_pix)
            shift = np.sqrt((xc - x0)**2 + (yc - y0)**2) * plate_scale
            print(f"    Refined centroid (pix): ({xc:.2f}, {yc:.2f}), "
                  f"shift: {shift:.2f} arcsec")

            # Extract radial profile out to max plotting radius
            r_mid, sb, sb_err, counts = extract_radial_profile(
                x_events, y_events, xc, yc, plate_scale,
                rbin, rmax_plot, exposure)

            total_cts = np.sum(counts)
            print(f"    Events within {rmax_plot}\": {total_cts}")

            # Fit King profile
            print(f"    Fitting King profile in {rmin_fit}-{rmax_fit}\" ...")
            popt, pcov, mask_fit = fit_king_profile(
                r_mid, sb, sb_err, rmin_fit, rmax_fit, rc_fixed,
                beta_fixed)

            if popt is not None:
                print(f"    Best fit: S0={popt[0]:.3e}  "
                      f"rc={popt[1]:.2f}\"  beta={popt[2]:.3f}  "
                      f"bkg={popt[3]:.3e}")

            # Detect pile-up radius
            pileup_radius, flag_level, resid_sigma = \
                detect_pileup_radius(
                    r_mid, sb, sb_err, popt, sigma1, sigma2, rbin)

            if pileup_radius is not None:
                print(f"    Pile-up detected: recommended inner radius "
                      f"= {pileup_radius:.1f}\"")
            else:
                print(f"    No pile-up detected at "
                      f"{sigma1:.1f}/{sigma2:.1f} sigma thresholds.")

            # Plot — use event filename stem for unique naming
            evt_stem = os.path.basename(filepath).replace('_cl.evt.gz', '').replace('_cl.evt', '')
            window_info = get_window_info(filepath)
            print(f"    {window_info['label']}")

            # Check for user override of pile-up radius
            override_radius = overrides.get(evt_stem, None)
            if override_radius is not None:
                print(f"    Override pile-up radius: {override_radius:.1f}\""
                      f" (matched '{evt_stem}')")
            elif overrides:
                print(f"    No override for '{evt_stem}'")

            output_png = os.path.join(output_dir,
                                       f'{evt_stem}_king_profile.png')
            fig = make_profile_plot(obsid, evt_stem, r_mid, sb, sb_err,
                              popt, pcov,
                              mask_fit, rmin_fit, rmax_fit, rmax_plot,
                              exposure, len(x_events), window_info,
                              sigma1, sigma2, pileup_radius,
                              override_radius, flag_level, output_png)
            print(f"    Saved: {output_png}")

            # Write pile-up radius to text file
            output_txt = os.path.join(output_dir,
                                       f'{evt_stem}_pileup.txt')
            with open(output_txt, 'w') as ftxt:
                ftxt.write(f"# Pile-up analysis for {evt_stem}\n")
                ftxt.write(f"# Source position (input): "
                           f"RA={ra_src:.6f} Dec={dec_src:.6f}\n")
                ftxt.write(f"# Centroid position (pix): "
                           f"X={xc:.2f} Y={yc:.2f}\n")
                ftxt.write(f"# Plate scale: "
                           f"{plate_scale:.4f} arcsec/pixel\n")
                ftxt.write(f"# Sigma thresholds: "
                           f"{sigma1:.1f} (2 consec), "
                           f"{sigma2:.1f} (single)\n")
                ftxt.write(f"# King PSF: rc={popt[1]:.2f}\" "
                           f"beta={popt[2]:.3f} (fixed)\n"
                           if popt is not None else
                           f"# King PSF fit failed\n")
                ftxt.write(f"# Window: {window_info['label']}\n")
                ftxt.write(f"# Count rate: "
                           f"{len(x_events)/exposure:.3f} ct/s\n")
                if pileup_radius is not None:
                    ftxt.write(f"pileup_radius_arcsec = "
                               f"{pileup_radius:.1f}\n")
                else:
                    ftxt.write(f"pileup_radius_arcsec = none\n")
                if override_radius is not None:
                    ftxt.write(f"override_radius_arcsec = "
                               f"{override_radius:.1f}\n")
            print(f"    Saved: {output_txt}")

            return output_png, fig

    except Exception as e:
        print(f"    ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None, None


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fit King profiles to Swift XRT PC-mode radial '
                    'profiles across OBSID directories.')
    parser.add_argument('--ra', type=float, required=True,
                        help='Source RA in degrees')
    parser.add_argument('--dec', type=float, required=True,
                        help='Source Dec in degrees')
    parser.add_argument('--rmin', type=float, default=20.0,
                        help='Inner radius of fitting annulus (arcsec, '
                             'default: 20)')
    parser.add_argument('--rmax', type=float, default=60.0,
                        help='Outer radius of fitting annulus (arcsec, '
                             'default: 60)')
    parser.add_argument('--rbin', type=float, default=2.0,
                        help='Radial bin width (arcsec, default: 2)')
    parser.add_argument('--maxplot', type=float, default=None,
                        help='Maximum radius to plot (arcsec, '
                             'default: rmax + 20)')
    parser.add_argument('--centroid', type=float, default=15.0,
                        help='Centroid search radius (arcsec, default: 15)')
    parser.add_argument('--rc', type=float, default=5.8,
                        help='King profile core radius in arcsec, fixed '
                             '(default: 5.8)')
    parser.add_argument('--beta', type=float, default=1.55,
                        help='King profile beta slope, fixed '
                             '(default: 1.55)')
    parser.add_argument('--sigma', type=float, default=3.0,
                        help='Sigma threshold for pile-up flagging, '
                             '2 consecutive bins required '
                             '(default: 3.0)')
    parser.add_argument('--sigma2', type=float, default=4.0,
                        help='Higher sigma threshold for pile-up, '
                             'single bin sufficient '
                             '(default: 4.0)')
    parser.add_argument('--pdf', type=str, default='king_profiles.pdf',
                        help='Output multi-page PDF filename '
                             '(default: king_profiles.pdf)')
    args = parser.parse_args()

    if args.maxplot is None:
        args.maxplot = args.rmax + 20.0

    script_dir = os.getcwd()

    # Find OBSID directories
    obsid_pattern = re.compile(r'^\d{8,11}$')
    obsid_dirs = sorted([
        d for d in os.listdir(script_dir)
        if os.path.isdir(os.path.join(script_dir, d))
        and obsid_pattern.match(d)
    ])

    if not obsid_dirs:
        print("No OBSID directories found.")
        sys.exit(1)

    # Load pile-up radius overrides if file exists
    overrides = {}
    override_file = os.path.join(script_dir, 'pileup_overrides.txt')
    if os.path.exists(override_file):
        print(f"Loading pile-up overrides from {override_file}")
        with open(override_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        overrides[parts[0]] = float(parts[1])
                        print(f"  Override: '{parts[0]}' -> "
                              f"{parts[1]} arcsec")
                    except ValueError:
                        print(f"  WARNING: could not parse: {line}")
        print(f"  Loaded {len(overrides)} override(s).")
    
    print(f"\nSource position: RA={args.ra:.6f}, Dec={args.dec:.6f}")
    print(f"Fitting annulus: {args.rmin}-{args.rmax} arcsec")
    print(f"King PSF: rc={args.rc}\" (fixed), beta={args.beta} (fixed)")
    print(f"Pile-up sigma threshold: {args.sigma} (2 consec), "
          f"{args.sigma2} (single)")
    print(f"Radial bin width: {args.rbin} arcsec")
    print(f"Plot range: 0-{args.maxplot} arcsec")
    print(f"Found {len(obsid_dirs)} OBSID directories.\n")

    all_pngs = []
    all_figs = []

    for obsid in obsid_dirs:
        obsid_path = os.path.join(script_dir, obsid)
        print(f"\n{'='*60}")
        print(f"OBSID: {obsid}")
        print(f"{'='*60}")

        # Find PC-mode cleaned event files
        # Match patterns: xpc*_cl.evt and xpc*_cl.evt.gz
        pc_files = glob.glob(os.path.join(obsid_path, '**', '*xpc*_cl.evt'),
                             recursive=True)
        pc_files += glob.glob(os.path.join(obsid_path, '**',
                                           '*xpc*_cl.evt.gz'),
                              recursive=True)

        # Exclude slew files — focus on pointed/settled data
        pc_files = [f for f in pc_files
                    if 'sl_cl' not in os.path.basename(f).lower()]

        if not pc_files:
            print("  No PC-mode pointed/settled event files found.")
            continue

        for pc_file in pc_files:
            result = process_event_file(
                pc_file, args.ra, args.dec,
                args.rmin, args.rmax, args.rbin, args.maxplot,
                args.centroid, args.rc, args.beta,
                args.sigma, args.sigma2, obsid, obsid_path,
                overrides)
            png, fig = result
            if png is not None:
                all_pngs.append(png)
                all_figs.append(fig)

    # Collate into multi-page PDF
    if all_figs:
        pdf_path = os.path.join(script_dir, args.pdf)
        print(f"\n{'='*60}")
        print(f"Collating {len(all_figs)} figures into {pdf_path}")
        print(f"{'='*60}")

        with PdfPages(pdf_path) as pdf:
            for fig in all_figs:
                pdf.savefig(fig, dpi=150, bbox_inches='tight')
                plt.close(fig)

        print(f"PDF saved: {pdf_path}")
    else:
        print("\nNo profiles were successfully generated.")


if __name__ == '__main__':
    main()

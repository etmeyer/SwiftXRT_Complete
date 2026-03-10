#!/usr/bin/env python3
"""
swift_xrt_fit_spectra.py

Sequentially fit grouped Swift XRT spectra (PC and/or WT mode)
using Sherpa.

Two model options (--model):
  absorbed: xstbabs * xsztbabs * xspowerlaw (requires --redshift)
  simple:   xstbabs * xspowerlaw (Galactic absorption only)

Reads pc_master_table.txt and/or wt_master_table.txt to
determine which observations to fit. Both modes are processed
in a single run with results combined in the output table
and light curve plot.

For low-count spectra (<200 cts), gamma is frozen to a default value.
Spectra with <40 counts are skipped entirely.

Usage:
    python swift_xrt_fit_spectra.py --nh 0.0179 --redshift 0.158
    python swift_xrt_fit_spectra.py --nh 0.0179 --model simple
    python swift_xrt_fit_spectra.py --nh 0.0179 --redshift 0.158 \\
        --modes pc          # PC only
    python swift_xrt_fit_spectra.py --nh 0.0179 --redshift 0.158 \\
        --modes wt --nmax 5 # WT only, first 5

Output:
    {stem}_sherpa_fit.log   per-observation fit log
    fit_results.txt         summary table of all fits
    flux_lightcurve.pdf     νFν and gamma vs time plot
"""

import os
import sys
import re
import glob
import argparse
import warnings
import numpy as np

try:
    from astropy.io import fits
    from astropy.time import Time
except ImportError:
    print("ERROR: astropy is required.")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("ERROR: matplotlib is required.")
    sys.exit(1)

# Sherpa imports — fail early if not installed
try:
    from sherpa.astro import ui as shp
    from sherpa.utils.err import EstErr, FitErr
except ImportError:
    print("ERROR: sherpa is required.")
    print("  Install via CIAO or: pip install sherpa")
    sys.exit(1)


BASE_DIR = os.path.abspath(os.getcwd())


# ---------------------------------------------------------------
# Read master tables
# ---------------------------------------------------------------

def read_master_table(table_file, mode='PC'):
    """
    Read a master table (PC or WT format), return included entries.

    Both formats have OBSID and filename as the first two columns,
    include as the third, and a quoted comment at the end. The
    middle columns differ but are not needed for fitting.

    Each entry gets a 'mode' field set to the supplied mode tag.
    """
    entries = []
    with open(table_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or \
                    line.startswith('OBSID'):
                continue
            quote_match = re.search(r'"([^"]*)"', line)
            comment = quote_match.group(1) if quote_match else ""
            before_quote = line[:quote_match.start()].strip() \
                if quote_match else line
            parts = before_quote.split()
            if len(parts) < 3:
                continue
            entries.append({
                'obsid': parts[0],
                'filename': parts[1],
                'include': parts[2].lower(),
                'mode': mode,
                'comment': comment,
            })
    return [e for e in entries if e['include'] == 'yes']


# ---------------------------------------------------------------
# Observation metadata
# ---------------------------------------------------------------

def get_obs_metadata(grp_pha):
    """
    Extract counts, exposure, dates from grouped spectrum.
    Computes midpoint time and MJD.
    """
    info = {}
    with fits.open(grp_pha) as hdul:
        header = hdul[1].header
        data = hdul[1].data
        info['total_counts'] = int(np.sum(data['COUNTS']))
        info['exposure'] = header.get('EXPOSURE', 0)
        info['date_obs'] = header.get('DATE-OBS', None)
        info['date_end'] = header.get('DATE-END', None)

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
    elif info['date_obs']:
        try:
            t = Time(info['date_obs'], format='isot', scale='utc')
            info['date_mid'] = t.isot
            info['mjd_mid'] = t.mjd
        except Exception:
            info['date_mid'] = info['date_obs']
            info['mjd_mid'] = None
    else:
        info['date_mid'] = 'N/A'
        info['mjd_mid'] = None

    info['count_rate'] = info['total_counts'] / info['exposure'] \
        if info['exposure'] > 0 else 0

    return info


# ---------------------------------------------------------------
# Resolve $CALDB in PHA header paths
# ---------------------------------------------------------------

def _resolve_caldb(path_str, caldb_override=None):
    """
    Resolve $CALDB in a file path string to the actual directory.

    Priority:
      1. caldb_override (from --caldb command line argument)
      2. CALDB environment variable

    This is needed because running in a CIAO environment sets
    $CALDB to the Chandra CALDB, not the HEASoft CALDB where
    Swift XRT calibration files live.
    """
    if '$CALDB' in path_str:
        caldb = caldb_override or os.environ.get('CALDB', '')
        if caldb:
            return path_str.replace('$CALDB', caldb)
    return path_str


def _get_header_paths(pha_file, caldb_override=None):
    """
    Read RESPFILE, ANCRFILE, BACKFILE from PHA header and
    resolve any $CALDB references.

    Returns dict with resolved paths (or None if not set).
    """
    paths = {'rmf': None, 'arf': None, 'bkg': None}
    with fits.open(pha_file) as hdul:
        header = hdul[1].header
        resp = header.get('RESPFILE', 'none').strip()
        ancr = header.get('ANCRFILE', 'none').strip()
        back = header.get('BACKFILE', 'none').strip()

    pha_dir = os.path.dirname(os.path.abspath(pha_file))

    for key, val, label in [('rmf', resp, 'RESPFILE'),
                             ('arf', ancr, 'ANCRFILE'),
                             ('bkg', back, 'BACKFILE')]:
        if val.lower() == 'none' or val == '':
            continue
        resolved = _resolve_caldb(val, caldb_override)
        # If not absolute, make relative to PHA directory
        if not os.path.isabs(resolved):
            resolved = os.path.join(pha_dir, resolved)
        if os.path.exists(resolved):
            paths[key] = resolved
        else:
            print(f"    WARNING: {label} not found: {resolved}")

    return paths


# ---------------------------------------------------------------
# Fit one spectrum with Sherpa
# ---------------------------------------------------------------

def fit_spectrum(grp_pha, nh_gal, redshift, gamma_value,
                 freeze_gamma, emin, emax, caldb_override=None,
                 bkg_mode='subtract', model_type='absorbed'):
    """
    Fit a single grouped PHA spectrum using Sherpa.

    Model: xstbabs.gal * xsztbabs.intr * xspowerlaw.pl

    Parameters set:
      gal.nH       = nh_gal    (frozen)
      intr.nH      = 0.01      (free, initial guess)
      intr.Redshift = redshift  (frozen)
      pl.PhoIndex  = gamma_value (free or frozen)
      pl.norm      = 1e-3      (free)

    Returns dict with fit results, or None on failure.
    """
    results = {
        'gamma': None, 'gamma_err': None,
        'gamma_lo': None, 'gamma_hi': None,
        'nh_int': None, 'nh_int_lo': None, 'nh_int_hi': None,
        'norm': None, 'norm_lo': None, 'norm_hi': None,
        'chi2': None, 'dof': None, 'reduced_chi2': None,
        'flux_band': None, 'flux_band_err': None,
        'flux_1keV': None, 'flux_1keV_err': None,
        'gamma_frozen': freeze_gamma,
    }

    try:
        # Clean sherpa state from any previous fit
        shp.clean()

        # Resolve $CALDB in any file paths referenced by the PHA
        # header, since Sherpa cannot expand shell variables.
        resolved = _get_header_paths(grp_pha, caldb_override)

        # Load the grouped spectrum.
        # Sherpa will warn about $CALDB paths it can't find —
        # suppress all warnings/output since we manually load
        # responses below with resolved paths.
        import logging
        sherpa_logger = logging.getLogger('sherpa')
        prev_level = sherpa_logger.level
        sherpa_logger.setLevel(logging.ERROR)
        _saved_stderr = sys.stderr
        _saved_stdout = sys.stdout
        try:
            devnull = open(os.devnull, 'w')
            sys.stderr = devnull
            sys.stdout = devnull
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                shp.load_pha(1, grp_pha)
        finally:
            sys.stderr = _saved_stderr
            sys.stdout = _saved_stdout
            devnull.close()
            sherpa_logger.setLevel(prev_level)

        # Manually load RMF/ARF/BKG with resolved paths.
        # This overwrites whatever Sherpa auto-loaded (or failed
        # to load) from the PHA header.
        # Suppress ENERG_LO warnings — these are harmless; Sherpa
        # replaces zero-energy boundaries with a tiny value.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if resolved['rmf'] is not None:
                print(f"    RMF: {os.path.basename(resolved['rmf'])}")
                shp.load_rmf(1, resolved['rmf'])
            if resolved['arf'] is not None:
                print(f"    ARF: {os.path.basename(resolved['arf'])}")
                shp.load_arf(1, resolved['arf'])
            if resolved['bkg'] is not None and bkg_mode == 'subtract':
                shp.load_bkg(1, resolved['bkg'])
                # Subtract the area-scaled background from the source.
                # BACKSCAL keywords handle the area scaling.
                shp.subtract(1)
            elif bkg_mode == 'none':
                pass

        # Set analysis to energy units (keV) so that ignore/notice
        # commands accept energy values rather than channel integers.
        # This requires a valid RMF to define the energy grid.
        shp.set_analysis(1, "energy")

        # Ignore bad channels first, then restrict to our range
        shp.ignore_bad(1)
        shp.notice(emin, emax)

        # Check we have enough noticed channels
        n_noticed = shp.get_data(1).get_dep(True).size
        if n_noticed < 3:
            print(f"    WARNING: only {n_noticed} noticed bins, "
                  f"skipping.")
            return None

        # Define the model.
        # Suppress the tbvabs version banner that the XSPEC model
        # library prints to stdout/stderr on first use.
        _so, _se = sys.stdout, sys.stderr
        devnull = open(os.devnull, 'w')
        try:
            sys.stdout = devnull
            sys.stderr = devnull
            if model_type == 'absorbed':
                # Full model: Galactic + intrinsic absorption
                #   xstbabs:    Galactic ISM absorption (frozen)
                #   xsztbabs:   redshifted intrinsic absorption (free)
                #   xspowerlaw: power law continuum
                shp.set_source(1,
                    "xstbabs.gal * xsztbabs.intr * xspowerlaw.pl")
            else:
                # Simple model: Galactic absorption only
                #   xstbabs:    Galactic ISM absorption (frozen)
                #   xspowerlaw: power law continuum
                shp.set_source(1,
                    "xstbabs.gal * xspowerlaw.pl")
        finally:
            sys.stdout, sys.stderr = _so, _se
            devnull.close()

        # --- Galactic absorption (always frozen) ---
        gal = shp.get_model_component("gal")
        gal.nH = nh_gal
        shp.freeze(gal.nH)

        # --- Intrinsic absorption (absorbed model only) ---
        if model_type == 'absorbed':
            intr = shp.get_model_component("intr")
            intr.nH = 0.01
            intr.nH.min = 0.0
            intr.nH.max = 100.0
            intr.Redshift = redshift
            shp.freeze(intr.Redshift)
            shp.thaw(intr.nH)

        # --- Power law ---
        pl = shp.get_model_component("pl")
        pl.PhoIndex = gamma_value
        pl.PhoIndex.min = 0.5
        pl.PhoIndex.max = 5.0
        pl.norm = 1.0e-3
        pl.norm.min = 1.0e-10

        if freeze_gamma:
            shp.freeze(pl.PhoIndex)
        else:
            shp.thaw(pl.PhoIndex)

        # --- Statistic and method ---
        # For subtracted data, use chi2datavar which uses the
        # observed variance (appropriate when background has been
        # subtracted and bins have enough counts from grouping).
        # For unsubtracted data, use chi2gehrels which applies
        # the Gehrels (1986) Poisson approximation.
        if bkg_mode == 'subtract':
            shp.set_stat("chi2datavar")
        else:
            shp.set_stat("chi2gehrels")
        shp.set_method("levmar")

        # --- Fit ---
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shp.fit(1)

        # Extract fit statistic
        fr = shp.get_fit_results()
        results['chi2'] = fr.statval
        results['dof'] = fr.dof
        if fr.dof > 0:
            results['reduced_chi2'] = fr.statval / fr.dof

        # Extract best-fit parameters
        results['gamma'] = float(pl.PhoIndex.val)
        results['norm'] = float(pl.norm.val)
        if model_type == 'absorbed':
            results['nh_int'] = float(intr.nH.val)
        else:
            results['nh_int'] = None

        # --- Confidence intervals ---
        # Use conf() for 90% confidence (delta chi2 = 2.706),
        # which matches XSPEC's default "error" behavior.
        # The "hard minimum/maximum hit" warnings are normal —
        # they mean conf() explored to a parameter boundary while
        # mapping the error surface. The fit values are unaffected.
        try:
            shp.set_conf_opt("sigma", 1.6)  # ~90% for 1 param
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if model_type == 'absorbed':
                    if not freeze_gamma:
                        shp.conf(1, pl.PhoIndex, intr.nH, pl.norm)
                    else:
                        shp.conf(1, intr.nH, pl.norm)
                else:
                    if not freeze_gamma:
                        shp.conf(1, pl.PhoIndex, pl.norm)
                    else:
                        shp.conf(1, pl.norm)

            cr = shp.get_conf_results()

            for i, pname in enumerate(cr.parnames):
                pmin = cr.parmins[i]
                pmax = cr.parmaxes[i]

                if 'PhoIndex' in pname:
                    if pmin is not None and pmax is not None:
                        results['gamma_lo'] = \
                            results['gamma'] + pmin
                        results['gamma_hi'] = \
                            results['gamma'] + pmax
                        results['gamma_err'] = (pmax - pmin) / 2.0
                elif 'nH' in pname and 'gal' not in pname:
                    if pmin is not None and pmax is not None:
                        results['nh_int_lo'] = \
                            results['nh_int'] + pmin
                        results['nh_int_hi'] = \
                            results['nh_int'] + pmax
                elif 'norm' in pname:
                    if pmin is not None and pmax is not None:
                        results['norm_lo'] = \
                            results['norm'] + pmin
                        results['norm_hi'] = \
                            results['norm'] + pmax

        except (EstErr, FitErr) as e:
            print(f"    WARNING: confidence failed: {e}")

        if freeze_gamma:
            results['gamma_err'] = 0.0

        # --- Band-integrated flux ---
        # calc_energy_flux returns absorbed flux in erg/cm²/s
        try:
            results['flux_band'] = float(
                shp.calc_energy_flux(id=1, lo=emin, hi=emax))
        except Exception as e:
            print(f"    WARNING: flux calculation failed: {e}")

        # Use sample_energy_flux for flux error via MC sampling.
        # This properly accounts for parameter correlations.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                samples = shp.sample_energy_flux(
                    id=1, lo=emin, hi=emax, num=500,
                    correlated=True)
            flux_dist = samples[:, 0]
            results['flux_band'] = float(np.median(flux_dist))
            results['flux_band_err'] = float(np.std(flux_dist))
        except Exception as e_sample:
            # Log why MC sampling failed — common causes include
            # covariance matrix issues or parameter space problems
            print(f"    NOTE: MC flux sampling failed ({e_sample}), "
                  f"using covariance fallback.")
            # Fallback: estimate flux error from the covariance
            # matrix via covar(). This is less rigorous than MC
            # sampling but gives a reasonable error estimate.
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    shp.covar(1)
                cov = shp.get_covar_results()
                # Fractional error from norm (dominant for powerlaw)
                for i, pname in enumerate(cov.parnames):
                    if 'norm' in pname:
                        sigma = cov.parmaxes[i]
                        if sigma is not None and \
                                results['norm'] is not None and \
                                results['norm'] > 0 and \
                                results['flux_band'] is not None:
                            frac_err = abs(sigma) / results['norm']
                            results['flux_band_err'] = \
                                results['flux_band'] * frac_err
                        break
            except Exception:
                pass

        # --- 1 keV flux density ---
        # The XSPEC powerlaw norm is photons/keV/cm²/s at 1 keV.
        # Convert to F_ν (erg/cm²/s/Hz):
        #   1 keV = 2.418e17 Hz
        #   1 keV = 1.602e-9 erg
        #   F_ν = norm × (1.602e-9 erg) / (2.418e17 Hz)
        #       = norm × 6.626e-27 erg/cm²/s/Hz
        keV_to_Hz = 2.418e17
        keV_to_erg = 1.602e-9

        if results['norm'] is not None:
            results['flux_1keV'] = \
                results['norm'] * keV_to_erg / keV_to_Hz

            if results['norm_lo'] is not None and \
                    results['norm_hi'] is not None:
                norm_err = (results['norm_hi'] -
                            results['norm_lo']) / 2.0
                results['flux_1keV_err'] = \
                    abs(norm_err) * keV_to_erg / keV_to_Hz

    except Exception as e:
        print(f"    ERROR during fit: {e}")
        import traceback
        traceback.print_exc()
        return None

    return results


# ---------------------------------------------------------------
# Process one observation
# ---------------------------------------------------------------

def process_one(entry, nh_gal, redshift, defgamma,
                min_counts_fit, min_counts_gamma, emin, emax,
                caldb_override=None, bkg_mode='subtract',
                model_type='absorbed'):
    """Fit one grouped spectrum. Returns results dict or None."""
    obsid = entry['obsid']
    stem = entry['filename']
    obsid_path = os.path.join(BASE_DIR, obsid)
    grp_pha = os.path.join(obsid_path, f'{stem}_grp.pha')

    if not os.path.exists(grp_pha):
        print(f"    Grouped spectrum not found: {stem}_grp.pha")
        return None

    meta = get_obs_metadata(grp_pha)
    counts = meta['total_counts']
    print(f"    Counts: {counts}  Exposure: "
          f"{meta['exposure']:.1f}s  Date: {meta['date_mid']}")

    if counts < min_counts_fit:
        print(f"    SKIPPED: {counts} counts < {min_counts_fit}")
        return None

    if counts < min_counts_gamma:
        freeze_gamma = True
        gamma_value = defgamma
        print(f"    Low counts ({counts} < {min_counts_gamma}): "
              f"freezing gamma={gamma_value}")
    else:
        freeze_gamma = False
        gamma_value = 2.0
        print(f"    Free gamma fit (initial guess: {gamma_value})")

    # Run fit
    fit = fit_spectrum(grp_pha, nh_gal, redshift,
                       gamma_value, freeze_gamma, emin, emax,
                       caldb_override, bkg_mode, model_type)

    if fit is None:
        print(f"    FIT FAILED.")
        return None

    # Attach metadata
    fit['obsid'] = obsid
    fit['stem'] = stem
    fit['mode'] = entry.get('mode', 'PC')
    fit['count_rate'] = meta['count_rate']
    fit['exposure'] = meta['exposure']
    fit['total_counts'] = counts
    fit['date_mid'] = meta['date_mid']
    fit['mjd_mid'] = meta['mjd_mid']

    # Write per-observation log
    log_file = os.path.join(obsid_path, f'{stem}_sherpa_fit.log')
    model_str = 'xstbabs * xsztbabs * xspowerlaw' \
        if model_type == 'absorbed' else 'xstbabs * xspowerlaw'
    with open(log_file, 'w') as f:
        f.write(f"# Sherpa fit log for {stem}\n")
        f.write(f"# Mode: {fit['mode']}\n")
        f.write(f"# Model: {model_str}\n\n")
        f.write(f"Spectrum     : {stem}_grp.pha\n")
        f.write(f"Counts       : {counts}\n")
        f.write(f"Exposure     : {meta['exposure']:.1f} s\n")
        f.write(f"Count rate   : {meta['count_rate']:.3f} ct/s\n")
        f.write(f"Date (mid)   : {meta['date_mid']}\n")
        f.write(f"Energy range : {emin}-{emax} keV\n\n")
        f.write(f"nH (Galactic): {nh_gal} (frozen)\n")
        if model_type == 'absorbed':
            f.write(f"nH (intrinsic): {fit['nh_int']}")
            if fit['nh_int_lo'] is not None:
                f.write(f"  ({fit['nh_int_lo']:.4f} - "
                        f"{fit['nh_int_hi']:.4f})")
            f.write(f"\n")
            f.write(f"Redshift     : {redshift} (frozen)\n")
        f.write(f"Gamma        : {fit['gamma']:.4f}")
        if fit['gamma_err'] and fit['gamma_err'] > 0:
            f.write(f"  +/- {fit['gamma_err']:.4f}")
        elif freeze_gamma:
            f.write(f"  (frozen)")
        f.write(f"\n")
        f.write(f"Norm         : {fit['norm']:.4e}")
        if fit['norm_lo'] is not None:
            f.write(f"  ({fit['norm_lo']:.4e} - "
                    f"{fit['norm_hi']:.4e})")
        f.write(f"\n\n")
        if fit['chi2'] is not None:
            f.write(f"Chi2/dof     : {fit['chi2']:.2f}/"
                    f"{fit['dof']} = "
                    f"{fit['reduced_chi2']:.3f}\n")
        if fit['flux_band'] is not None:
            f.write(f"Flux ({emin}-{emax} keV): "
                    f"{fit['flux_band']:.4e} erg/cm2/s")
            if fit['flux_band_err'] is not None:
                f.write(f"  +/- {fit['flux_band_err']:.4e}")
            f.write(f"\n")
        if fit['flux_1keV'] is not None:
            f.write(f"F_nu(1keV)   : {fit['flux_1keV']:.4e} "
                    f"erg/cm2/s/Hz")
            if fit['flux_1keV_err'] is not None:
                f.write(f"  +/- {fit['flux_1keV_err']:.4e}")
            f.write(f"\n")

    # Print summary
    if fit['gamma'] is not None:
        g_str = f"{fit['gamma']:.3f}"
        if fit['gamma_err'] and fit['gamma_err'] > 0:
            g_str += f" +/- {fit['gamma_err']:.3f}"
        elif freeze_gamma:
            g_str += " (frozen)"
        print(f"    Gamma: {g_str}")
    if fit['flux_band'] is not None:
        fb_str = f"{fit['flux_band']:.3e}"
        if fit['flux_band_err']:
            fb_str += f" +/- {fit['flux_band_err']:.3e}"
        print(f"    Flux ({emin}-{emax} keV): {fb_str} erg/cm²/s")
    if fit['flux_1keV'] is not None:
        f1_str = f"{fit['flux_1keV']:.3e}"
        if fit['flux_1keV_err']:
            f1_str += f" +/- {fit['flux_1keV_err']:.3e}"
        print(f"    F_ν(1keV): {f1_str} erg/cm²/s/Hz")
    if fit['reduced_chi2'] is not None:
        print(f"    χ²/dof: {fit['chi2']:.1f}/{fit['dof']} = "
              f"{fit['reduced_chi2']:.2f}")

    return fit


# ---------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------

def write_summary_table(results, output_file, model_type='absorbed'):
    """Write results table to file and terminal."""
    header = (
        f"{'OBSID':<14} {'filename':<30} {'mode':>4} {'ctrate':>8} "
        f"{'exp(s)':>8} {'DateObs':<22} {'MJD':>12} "
        f"{'flux_band':>12} {'fband_err':>12} "
        f"{'flux_1keV':>12} {'f1_err':>12} "
        f"{'gamma':>7} {'gamma_err':>10}"
    )
    sep = '-' * len(header)
    lines = [header, sep]

    for r in results:
        ctrate = f"{r['count_rate']:.3f}" \
            if r.get('count_rate') else "N/A"
        exp = f"{r['exposure']:.1f}" \
            if r.get('exposure') else "N/A"
        date = r.get('date_mid', 'N/A') or 'N/A'
        mjd = f"{r['mjd_mid']:.5f}" \
            if r.get('mjd_mid') else "N/A"
        fb = f"{r['flux_band']:.3e}" \
            if r.get('flux_band') else "N/A"
        fbe = f"{r['flux_band_err']:.3e}" \
            if r.get('flux_band_err') else "N/A"
        f1 = f"{r['flux_1keV']:.3e}" \
            if r.get('flux_1keV') else "N/A"
        f1e = f"{r['flux_1keV_err']:.3e}" \
            if r.get('flux_1keV_err') else "N/A"
        gam = f"{r['gamma']:.3f}" \
            if r.get('gamma') is not None else "N/A"
        ge = f"{r['gamma_err']:.3f}" \
            if r.get('gamma_err') is not None else "N/A"
        if r.get('gamma_frozen'):
            ge = "(frozen)"
        mode = r.get('mode', 'PC')

        lines.append(
            f"{r['obsid']:<14} {r['stem']:<30} {mode:>4} {ctrate:>8} "
            f"{exp:>8} {date:<22} {mjd:>12} "
            f"{fb:>12} {fbe:>12} "
            f"{f1:>12} {f1e:>12} "
            f"{gam:>7} {ge:>10}"
        )

    lines.append(sep)
    table = '\n'.join(lines)

    print(f"\n{'='*len(header)}")
    print("  FIT RESULTS SUMMARY")
    print(f"{'='*len(header)}")
    print(table)

    with open(output_file, 'w') as f:
        model_str = 'xstbabs * xsztbabs * xspowerlaw' \
            if model_type == 'absorbed' else 'xstbabs * xspowerlaw'
        f.write("# Swift XRT spectral fit results (Sherpa)\n")
        f.write(f"# Model: {model_str}\n")
        f.write("# flux_band in erg/cm2/s, "
                "flux_1keV in erg/cm2/s/Hz\n\n")
        f.write(table + '\n')

    print(f"\nSaved: {output_file}")


# ---------------------------------------------------------------
# Light curve plot
# ---------------------------------------------------------------

def make_lightcurve_plot(results, output_pdf):
    """Two-panel landscape plot: νFν at 1 keV and gamma vs time,
    with PC and WT mode data in distinct colors."""
    valid = [r for r in results
             if r.get('mjd_mid') is not None
             and r.get('flux_1keV') is not None
             and r.get('gamma') is not None]

    if not valid:
        print("No valid results for light curve plot.")
        return

    mjd = np.array([r['mjd_mid'] for r in valid])
    flux = np.array([r['flux_1keV'] for r in valid])
    flux_err = np.array([
        r['flux_1keV_err'] if r.get('flux_1keV_err') else 0
        for r in valid])
    gamma = np.array([r['gamma'] for r in valid])
    gamma_err = np.array([
        r['gamma_err'] if r.get('gamma_err') and r['gamma_err'] > 0
        else 0 for r in valid])
    frozen = np.array([r.get('gamma_frozen', False) for r in valid])
    modes = np.array([r.get('mode', 'PC') for r in valid])

    order = np.argsort(mjd)
    mjd = mjd[order]; flux = flux[order]; flux_err = flux_err[order]
    gamma = gamma[order]; gamma_err = gamma_err[order]
    frozen = frozen[order]; modes = modes[order]

    # Convert MJD to decimal years
    decimal_year = 2000.0 + (mjd - 51544.0) / 365.25

    # Convert Fν to νFν
    nu_1keV = 2.418e17  # Hz
    vfv = flux * nu_1keV
    vfv_err = flux_err * nu_1keV

    # Mode masks
    pc_mask = modes == 'PC'
    wt_mask = modes == 'WT'

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True,
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.05})

    # --- Upper panel: νFν at 1 keV ---
    # PC mode: navy circles
    if np.any(pc_mask):
        ax1.errorbar(decimal_year[pc_mask], vfv[pc_mask],
                     yerr=vfv_err[pc_mask],
                     fmt='o', ms=5, color='navy',
                     ecolor='cornflowerblue',
                     elinewidth=1, capsize=3, label='PC')
    # WT mode: dark orange diamonds
    if np.any(wt_mask):
        ax1.errorbar(decimal_year[wt_mask], vfv[wt_mask],
                     yerr=vfv_err[wt_mask],
                     fmt='D', ms=5, color='darkorange',
                     ecolor='sandybrown',
                     elinewidth=1, capsize=3, label='WT')

    ax1.set_ylabel(r'$\nu F_{\nu}$(1 keV)  [erg cm$^{-2}$ s$^{-1}$]',
                    fontsize=11)
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_title('Swift XRT Light Curve', fontsize=13)
    ax1.legend(loc='best', fontsize=9)

    # --- Lower panel: photon index ---
    # PC free gamma: dark red circles
    m = pc_mask & ~frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=gamma_err[m],
                     fmt='o', ms=5, color='darkred',
                     ecolor='salmon', elinewidth=1, capsize=3,
                     label=r'PC free $\Gamma$')
    # PC frozen gamma: gray circles
    m = pc_mask & frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=0,
                     fmt='o', ms=5, color='silver',
                     elinewidth=1, capsize=3,
                     label=r'PC frozen $\Gamma$')
    # WT free gamma: dark orange diamonds
    m = wt_mask & ~frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=gamma_err[m],
                     fmt='D', ms=5, color='darkorange',
                     ecolor='sandybrown', elinewidth=1, capsize=3,
                     label=r'WT free $\Gamma$')
    # WT frozen gamma: light orange diamonds
    m = wt_mask & frozen
    if np.any(m):
        ax2.errorbar(decimal_year[m], gamma[m], yerr=0,
                     fmt='D', ms=5, color='wheat',
                     elinewidth=1, capsize=3,
                     label=r'WT frozen $\Gamma$')

    ax2.set_ylabel(r'$\Gamma$', fontsize=12)
    ax2.set_xlabel('Year', fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=8, ncol=2)

    # Set x-axis major ticks to integer years
    from matplotlib.ticker import MultipleLocator, AutoMinorLocator
    ax2.xaxis.set_major_locator(MultipleLocator(1))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))

    fig.savefig(output_pdf, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_pdf}")


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fit Swift XRT spectra (PC and/or WT) with Sherpa.')
    parser.add_argument('--nh', type=float, required=True,
                        help='Galactic nH in units of 10^22 cm^-2')
    parser.add_argument('--redshift', type=float, default=None,
                        help='Source redshift (required for '
                             '"absorbed" model)')
    parser.add_argument('--model', type=str, default='absorbed',
                        choices=['absorbed', 'simple'],
                        help='Spectral model: "absorbed" = '
                             'tbabs*ztbabs*powerlaw (Galactic + '
                             'intrinsic absorption), "simple" = '
                             'tbabs*powerlaw (Galactic only) '
                             '(default: absorbed)')
    parser.add_argument('--defgamma', type=float, default=2.0,
                        help='Default photon index for low-count '
                             'spectra (default: 2.0)')
    parser.add_argument('--mincounts', type=int, default=40,
                        help='Minimum counts to attempt fitting '
                             '(default: 40)')
    parser.add_argument('--mingamma', type=int, default=200,
                        help='Minimum counts for free gamma '
                             '(default: 200)')
    parser.add_argument('--bkg', type=str, default='subtract',
                        choices=['subtract', 'none'],
                        help='Background handling: "subtract" for '
                             'standard area-scaled subtraction, '
                             '"none" to ignore background entirely '
                             '(default: subtract)')
    parser.add_argument('--caldb', type=str, default=None,
                        help='Path to HEASoft CALDB (if $CALDB '
                             'points to CIAO CALDB instead). '
                             'e.g. /opt/CALDB')
    parser.add_argument('--emin', type=float, default=0.3,
                        help='Lower energy bound in keV '
                             '(default: 0.3)')
    parser.add_argument('--emax', type=float, default=10.0,
                        help='Upper energy bound in keV '
                             '(default: 10.0)')
    parser.add_argument('--modes', type=str, default='both',
                        choices=['pc', 'wt', 'both'],
                        help='Which mode(s) to fit: pc, wt, or both '
                             '(default: both)')
    parser.add_argument('--pctable', type=str,
                        default='pc_master_table.txt',
                        help='PC master table filename')
    parser.add_argument('--wttable', type=str,
                        default='wt_master_table.txt',
                        help='WT master table filename')
    parser.add_argument('--nmax', type=int, default=None,
                        help='Only process the first N observations '
                             '(for quick testing)')
    parser.add_argument('--output', type=str,
                        default='fit_results.txt',
                        help='Output results table filename')
    parser.add_argument('--plot', type=str,
                        default='flux_lightcurve.pdf',
                        help='Output light curve plot filename')
    args = parser.parse_args()

    # Read master tables based on --modes
    entries = []

    if args.modes in ('pc', 'both'):
        pc_path = os.path.join(BASE_DIR, args.pctable)
        if os.path.exists(pc_path):
            pc_entries = read_master_table(pc_path, mode='PC')
            entries.extend(pc_entries)
            print(f"PC table: {len(pc_entries)} observations "
                  f"(from {args.pctable})")
        elif args.modes == 'pc':
            print(f"ERROR: {args.pctable} not found.")
            sys.exit(1)
        else:
            print(f"NOTE: {args.pctable} not found, skipping PC.")

    if args.modes in ('wt', 'both'):
        wt_path = os.path.join(BASE_DIR, args.wttable)
        if os.path.exists(wt_path):
            wt_entries = read_master_table(wt_path, mode='WT')
            entries.extend(wt_entries)
            print(f"WT table: {len(wt_entries)} observations "
                  f"(from {args.wttable})")
        elif args.modes == 'wt':
            print(f"ERROR: {args.wttable} not found.")
            sys.exit(1)
        else:
            print(f"NOTE: {args.wttable} not found, skipping WT.")

    if not entries:
        print("No observations to fit.")
        sys.exit(0)

    # Sort all entries by OBSID for chronological processing
    entries.sort(key=lambda e: e['obsid'])

    # Validate: absorbed model requires redshift
    if args.model == 'absorbed' and args.redshift is None:
        print("ERROR: --redshift is required for the 'absorbed' model.")
        sys.exit(1)

    model_str = 'tbabs * ztbabs * powerlaw' \
        if args.model == 'absorbed' else 'tbabs * powerlaw'
    print(f"\nModel: {model_str}")
    print(f"Galactic nH: {args.nh} x 10^22 cm^-2")
    if args.model == 'absorbed':
        print(f"Redshift: {args.redshift}")
    print(f"Default gamma: {args.defgamma}")
    print(f"Energy range: {args.emin}-{args.emax} keV")
    print(f"Min counts to fit: {args.mincounts}")
    print(f"Min counts for free gamma: {args.mingamma}")
    if args.caldb:
        print(f"HEASoft CALDB: {args.caldb}")
    print(f"Background: {args.bkg}")
    print(f"Total observations: {len(entries)}")
    if args.nmax is not None:
        entries = entries[:args.nmax]
        print(f"  (limited to first {args.nmax})")

    # Suppress sherpa chatter during batch processing
    import logging
    logging.getLogger('sherpa').setLevel(logging.WARNING)

    # Fit each observation
    n_total = len(entries)
    all_results = []
    for i, entry in enumerate(entries, 1):
        print(f"\n  [{i}/{n_total}] [{entry['mode']}] "
              f"{entry['obsid']} / {entry['filename']}")
        result = process_one(
            entry, args.nh, args.redshift, args.defgamma,
            args.mincounts, args.mingamma, args.emin, args.emax,
            args.caldb, args.bkg, args.model)
        if result is not None:
            all_results.append(result)

    if not all_results:
        print("\nNo successful fits.")
        sys.exit(0)

    # Summary table
    output_path = os.path.join(BASE_DIR, args.output)
    write_summary_table(all_results, output_path, args.model)

    # Light curve plot
    plot_path = os.path.join(BASE_DIR, args.plot)
    make_lightcurve_plot(all_results, plot_path)

    n_pc = sum(1 for r in all_results if r.get('mode') == 'PC')
    n_wt = sum(1 for r in all_results if r.get('mode') == 'WT')
    print(f"\nDone. {len(all_results)} spectra fitted successfully "
          f"({n_pc} PC, {n_wt} WT).")


if __name__ == '__main__':
    main()

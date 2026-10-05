#!/usr/bin/env python3
"""
parallel_fit.py

Run swift_xrt_fit_spectra.py in parallel by splitting the
master tables into chunks. Each worker fits a subset of
observations independently.

As in a single run, spectra with at least --mingamma counts are
fitted first, with gamma free; the rest are then fitted with gamma
frozen at the median of all those fits (or --defgamma).

After all workers complete, results are merged into a single
fit_results.txt and the light curve plot is generated.

Usage:
    python parallel_fit.py --nh 0.0179 --model simple --caldb /opt/CALDB --nproc 16
    python parallel_fit.py --nh 0.0179 --redshift 0.158 --caldb /opt/CALDB --nproc 32
    python parallel_fit.py --nh 0.0179 --model simple --caldb /opt/CALDB --dryrun
"""

import os
import sys
import re
import argparse
import subprocess
import shutil
import tempfile
import importlib.util
import statistics
from concurrent.futures import ProcessPoolExecutor, as_completed

from swift_xrt_env import require_fit_caldb, sherpa_missing_message


BASE_DIR = os.path.abspath(os.getcwd())

# Directory where this script (and sibling scripts) live.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def read_master_table_raw(table_file):
    """Read master table, returning header and included entries."""
    header = None
    entries = []
    with open(table_file, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            if stripped.startswith('OBSID'):
                header = line
                continue
            quote_match = re.search(r'"([^"]*)"', stripped)
            before_quote = stripped[:quote_match.start()].strip() \
                if quote_match else stripped
            parts = before_quote.split()
            if len(parts) >= 3 and parts[2].lower() == 'yes':
                entries.append(line)
    return header, entries


def write_mini_table(header, entries, filepath):
    """Write a subset of entries to a file."""
    with open(filepath, 'w') as f:
        if header:
            f.write(header)
        for entry in entries:
            f.write(entry)


def run_fit_chunk(chunk_id, mini_tables, fit_args, base_dir, env):
    """
    Run the fitting script on a chunk.

    Each chunk gets its own temp directory with symlinks to
    OBSID dirs (where the _grp.pha files live). Output goes
    to a chunk-specific results file.

    env is the parent process's os.environ, passed explicitly
    to ensure CIAO/HEASoft environment variables propagate
    correctly through ProcessPoolExecutor workers.
    """
    tmp_dir = tempfile.mkdtemp(
        prefix=f'fit_chunk{chunk_id:02d}_', dir=base_dir)

    try:
        # Symlink OBSID directories
        obsid_pattern = re.compile(r'^\d{8,11}$')
        for d in os.listdir(base_dir):
            full = os.path.join(base_dir, d)
            if os.path.isdir(full) and obsid_pattern.match(d):
                link = os.path.join(tmp_dir, d)
                if not os.path.exists(link):
                    os.symlink(full, link)

        # Copy mini-tables into temp dir
        for mode, mpath in mini_tables.items():
            shutil.copy2(mpath, os.path.join(
                tmp_dir, os.path.basename(mpath)))

        # Build command
        output_name = f'.fit_results_chunk{chunk_id:02d}.txt'
        plot_name = f'.fit_plot_chunk{chunk_id:02d}.pdf'

        # A leftover from an interrupted earlier run must not be
        # merged as if this run had produced it.
        stale = os.path.join(base_dir, output_name)
        if os.path.exists(stale):
            os.remove(stale)

        cmd = [sys.executable,
               os.path.join(SCRIPT_DIR, 'swift_xrt_fit_spectra.py')]
        cmd.extend(fit_args)
        cmd.extend(['--output', output_name, '--plot', plot_name])

        # Override table paths
        for mode, mpath in mini_tables.items():
            if mode == 'pc':
                cmd.extend(['--pctable', os.path.basename(mpath)])
            elif mode == 'wt':
                cmd.extend(['--wttable', os.path.basename(mpath)])

        result = subprocess.run(
            cmd, cwd=tmp_dir,
            capture_output=True, text=True, timeout=7200,
            env=env)

        # Copy results file back to base dir
        results_path = os.path.join(tmp_dir, output_name)
        if os.path.exists(results_path):
            shutil.copy2(results_path,
                          os.path.join(base_dir, output_name))

        return {
            'chunk_id': chunk_id,
            'returncode': result.returncode,
            'results_file': output_name,
            'stdout_tail': result.stdout[-500:] if result.stdout else '',
            'stderr_tail': result.stderr[-500:] if result.stderr else '',
        }

    except Exception as e:
        return {
            'chunk_id': chunk_id,
            'returncode': -1,
            'results_file': None,
            'stdout_tail': '',
            'stderr_tail': str(e),
        }
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def read_results_file(filepath):
    """
    Split a fit_results table into its '#' comment lines, column
    header line, separator line and data lines.
    """
    comments, colhdr, sep, rows = [], None, None, []
    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                comments.append(line)
            elif set(stripped.replace(' ', '')) == {'-'}:
                sep = sep or line
            elif stripped.startswith('OBSID'):
                colhdr = colhdr or line
            else:
                rows.append(line)
    return comments, colhdr, sep, rows


def free_gammas(colhdr, rows):
    """Photon indices of the rows fitted with gamma free."""
    names = colhdr.split()
    ig, ie = names.index('gamma'), names.index('gamma_err')
    out = []
    for row in rows:
        parts = row.split()
        if len(parts) > ie and parts[ie] != '(frozen)':
            try:
                out.append(float(parts[ig]))
            except ValueError:
                pass
    return out


def merge_results(tables, output_file):
    """
    Write the chunk tables as one fit_results.txt: every distinct
    comment line, one column header, and the data rows sorted by
    OBSID and mode. Returns the number of rows.
    """
    comments, colhdr, sep, rows = [], None, None, []
    for c, h, s_, r in tables:
        comments += [line for line in c if line not in comments]
        colhdr, sep = colhdr or h, sep or s_
        rows += r
    # The units line last, as in a single run's table
    comments.sort(key=lambda line: line.startswith('# flux_band'))
    # By OBSID, then mode (an OBSID can have both), as a single run
    rows.sort(key=lambda l: (l.split() + ['', '', ''])[0:3:2])
    with open(os.path.join(BASE_DIR, output_file), 'w') as f:
        f.writelines(comments)
        f.write('\n')
        if colhdr:
            f.write(colhdr)
        if sep:
            f.write(sep)
        f.writelines(rows)
        if sep:
            f.write(sep)
    return len(rows)


def spectrum_counts(line):
    """Total counts in a master-table line's grouped spectrum, or None."""
    from astropy.io import fits
    import numpy as np
    parts = line.split()
    grp = os.path.join(BASE_DIR, parts[0], f'{parts[1]}_grp.pha')
    try:
        with fits.open(grp) as hdul:
            return int(np.sum(hdul[1].data['COUNTS']))
    except (OSError, KeyError, IndexError, TypeError):
        return None


def run_phase(entries, passthrough, nproc, first_chunk, parent_env):
    """
    Fit entries (mode, line, header) in up to nproc chunks. Returns
    (worker results, next free chunk id).
    """
    n_chunks = min(nproc, len(entries))
    chunk_size = (len(entries) + n_chunks - 1) // n_chunks
    jobs = []
    cid = first_chunk
    for i in range(n_chunks):
        chunk = entries[i*chunk_size:(i+1)*chunk_size]
        if not chunk:
            continue
        # Separate PC and WT entries
        pc_lines = [(hdr, line) for mode, line, hdr in chunk
                    if mode == 'pc']
        wt_lines = [(hdr, line) for mode, line, hdr in chunk
                    if mode == 'wt']
        mini_tables = {}
        if pc_lines:
            mpath = os.path.join(BASE_DIR, f'.pc_fit_chunk_{cid:02d}.txt')
            write_mini_table(pc_lines[0][0],
                             [l for _, l in pc_lines], mpath)
            mini_tables['pc'] = mpath
        if wt_lines:
            mpath = os.path.join(BASE_DIR, f'.wt_fit_chunk_{cid:02d}.txt')
            write_mini_table(wt_lines[0][0],
                             [l for _, l in wt_lines], mpath)
            mini_tables['wt'] = mpath
        chunk_args = list(passthrough)
        if mini_tables.keys() == {'pc'}:
            chunk_args.extend(['--modes', 'pc'])
        elif mini_tables.keys() == {'wt'}:
            chunk_args.extend(['--modes', 'wt'])
        else:
            chunk_args.extend(['--modes', 'both'])
        jobs.append((cid, mini_tables, chunk_args))
        cid += 1

    print(f"Launching {len(jobs)} fitting workers...", flush=True)
    results = []
    with ProcessPoolExecutor(max_workers=nproc) as executor:
        futures = {}
        for jid, mtabs, cargs in jobs:
            future = executor.submit(
                run_fit_chunk, jid, mtabs, cargs, BASE_DIR, parent_env)
            futures[future] = jid
        print(f"  All {len(futures)} workers submitted. "
              f"Waiting for results...\n", flush=True)
        for future in as_completed(futures):
            jid = futures[future]
            result = future.result()
            results.append(result)
            status = 'OK' if result['returncode'] == 0 else 'FAIL'
            print(f"  Chunk {jid:02d}: {status} "
                  f"[{len(results)}/{len(futures)} done]", flush=True)
            if result['returncode'] != 0:
                # Per-OBSID failures are reported on stdout,
                # crashes on stderr.
                for tail in (result['stdout_tail'],
                             result['stderr_tail']):
                    for line in tail.strip().splitlines()[-3:]:
                        print(f"    {line}", flush=True)

    for _, mtabs, _ in jobs:
        for mpath in mtabs.values():
            if os.path.exists(mpath):
                os.remove(mpath)
    return results, cid


def collect_tables(results):
    """Read and remove the chunk tables a phase wrote."""
    tables = []
    for r in results:
        if not r['results_file']:
            continue
        path = os.path.join(BASE_DIR, r['results_file'])
        if os.path.exists(path):
            tables.append(read_results_file(path))
            os.remove(path)
    return tables


def main():
    parser = argparse.ArgumentParser(
        description='Run spectral fitting in parallel.')
    parser.add_argument('--nproc', type=int, default=8,
                        help='Number of parallel workers (default: 8)')
    parser.add_argument('--dryrun', action='store_true',
                        help='Show splitting without running')

    # Pass-through to fit script
    parser.add_argument('--nh', type=float, required=True)
    parser.add_argument('--redshift', type=float, default=None)
    parser.add_argument('--model', type=str, default=None,
                        choices=['absorbed', 'simple'])
    parser.add_argument('--defgamma', type=float, default=None)
    parser.add_argument('--stat', type=str, default='wstat',
                        choices=['wstat', 'chi2'])
    parser.add_argument('--mincounts', type=int, default=40)
    parser.add_argument('--mingamma', type=int, default=200)
    parser.add_argument('--emin', type=float, default=0.3)
    parser.add_argument('--emax', type=float, default=10.0)
    parser.add_argument('--bkg', type=str, default='use',
                        choices=['use', 'none', 'subtract'])
    parser.add_argument('--caldb', type=str, default=None)
    parser.add_argument('--modes', type=str, default='both',
                        choices=['pc', 'wt', 'both'])
    parser.add_argument('--pctable', type=str,
                        default='pc_master_table.txt')
    parser.add_argument('--wttable', type=str,
                        default='wt_master_table.txt')
    parser.add_argument('--output', type=str,
                        default='fit_results.txt')
    parser.add_argument('--plot', type=str,
                        default='flux_lightcurve.pdf')
    args = parser.parse_args()

    # Same default as swift_xrt_fit_spectra.py: the intrinsic
    # absorber only when there is a redshift to put it at.
    if args.model is None:
        args.model = 'absorbed' if args.redshift is not None \
            else 'simple'
    if args.model == 'absorbed' and args.redshift is None:
        print("ERROR: --redshift is required for the 'absorbed' model.")
        sys.exit(1)

    # Workers run with this same Python and environment, so check
    # Sherpa and the CALDB once here instead of failing in every chunk.
    if not args.dryrun:
        if importlib.util.find_spec('sherpa') is None:
            print(sherpa_missing_message(), file=sys.stderr)
            sys.exit(1)
        require_fit_caldb(args.caldb)

    # Build pass-through arguments
    passthrough = ['--nh', str(args.nh), '--model', args.model,
                   '--stat', args.stat,
                   '--mincounts', str(args.mincounts),
                   '--mingamma', str(args.mingamma),
                   '--emin', str(args.emin),
                   '--emax', str(args.emax),
                   '--bkg', args.bkg]
    if args.redshift is not None:
        passthrough.extend(['--redshift', str(args.redshift)])
    if args.caldb:
        passthrough.extend(['--caldb', args.caldb])

    # Read all entries
    all_entries = []  # (mode, line, header)

    if args.modes in ('pc', 'both'):
        pc_path = os.path.join(BASE_DIR, args.pctable)
        if os.path.exists(pc_path):
            header_pc, entries_pc = read_master_table_raw(pc_path)
            for e in entries_pc:
                all_entries.append(('pc', e, header_pc))
            print(f"PC: {len(entries_pc)} observations")

    if args.modes in ('wt', 'both'):
        wt_path = os.path.join(BASE_DIR, args.wttable)
        if os.path.exists(wt_path):
            header_wt, entries_wt = read_master_table_raw(wt_path)
            for e in entries_wt:
                all_entries.append(('wt', e, header_wt))
            print(f"WT: {len(entries_wt)} observations")

    if not all_entries:
        print("No observations to fit.")
        sys.exit(0)

    # Two phases, as in swift_xrt_fit_spectra.py: gamma free for
    # spectra with at least --mingamma counts, then the rest with gamma
    # frozen at the median of all those fits (or --defgamma). Splitting
    # one run across chunks would otherwise give each chunk its own
    # median.
    low = []
    high = []
    for entry in all_entries:
        c = spectrum_counts(entry[1])
        if c is not None and args.mincounts <= c < args.mingamma:
            low.append(entry)
        else:
            high.append(entry)
    print(f"Total: {len(all_entries)} observations: {len(high)} with "
          f"gamma free, then {len(low)} under {args.mingamma} counts "
          f"with gamma frozen")

    if args.dryrun:
        print("\n[DRY RUN] No fitting performed.")
        return

    # Capture the current environment so it propagates correctly
    # to subprocesses spawned by workers. CIAO and HEASoft set
    # many environment variables (ASCDS_CALIB, CALDB, HEADAS,
    # LD_LIBRARY_PATH, etc.) that must be present for Sherpa
    # and XSPEC model libraries to function.
    parent_env = os.environ.copy()

    results, tables = [], []
    next_chunk = 0
    if high:
        print()
        r, next_chunk = run_phase(high, passthrough, args.nproc,
                                  next_chunk, parent_env)
        results += r
        tables += collect_tables(r)

    if low:
        if args.defgamma is not None:
            gamma, note = args.defgamma, '--defgamma'
        else:
            gammas = []
            for _, colhdr, _, rows in tables:
                if colhdr:
                    gammas += free_gammas(colhdr, rows)
            if gammas:
                gamma = statistics.median(gammas)
                note = f'median of {len(gammas)} free fits'
            else:
                gamma, note = 2.0, \
                    'no free fits to take a median from; default'
        print(f"\n{len(low)} spectra under {args.mingamma} counts: "
              f"gamma frozen at {gamma:.3f} ({note})")
        # '=' keeps a note such as '--defgamma' from reading as an option
        r, next_chunk = run_phase(
            low, passthrough + ['--defgamma', f'{gamma:.4f}',
                                f'--defgamma-note={note}'],
            args.nproc, next_chunk, parent_env)
        results += r
        tables += collect_tables(r)

    # Merge results. A chunk with one failed OBSID exits non-zero but
    # still wrote valid results for the rest, so merge every table
    # that exists.
    if tables:
        n_merged = merge_results(tables, args.output)
        print(f"\nMerged {n_merged} results into {args.output}")

        # Generate the combined plot using plot_lightcurve.py
        plot_script = os.path.join(SCRIPT_DIR, 'plot_lightcurve.py')
        if os.path.exists(plot_script):
            print(f"Generating light curve plot...")
            subprocess.run(
                [sys.executable, plot_script,
                 '--input', args.output,
                 '--output', args.plot],
                cwd=BASE_DIR)
        else:
            print(f"NOTE: plot_lightcurve.py not found, "
                  f"skipping plot. Run it manually.")
    else:
        print("\nNo results to merge.")

    # Clean up any leftover chunk plots
    for f in os.listdir(BASE_DIR):
        if f.startswith('.fit_plot_chunk') and f.endswith('.pdf'):
            os.remove(os.path.join(BASE_DIR, f))

    n_ok = sum(1 for r in results if r['returncode'] == 0)
    n_fail = sum(1 for r in results if r['returncode'] != 0)
    print(f"\nDone. {n_ok} chunks succeeded, {n_fail} failed.")
    if n_fail:
        print("  Failed chunks' output is above; their successful fits "
              "are still in the merged table.")
        sys.exit(1)


if __name__ == '__main__':
    main()

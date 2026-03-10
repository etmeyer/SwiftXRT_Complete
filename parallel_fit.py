#!/usr/bin/env python3
"""
parallel_fit.py

Run swift_xrt_fit_spectra.py in parallel by splitting the
master tables into chunks. Each worker fits a subset of
observations independently.

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
from concurrent.futures import ProcessPoolExecutor, as_completed


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


def merge_results(result_files, output_file):
    """
    Merge chunk result files into a single fit_results.txt.
    Preserves the header from the first file and concatenates
    data rows, sorted by OBSID.
    """
    header_lines = []
    data_lines = []
    header_done = False

    for rf in sorted(result_files):
        filepath = os.path.join(BASE_DIR, rf)
        if not os.path.exists(filepath):
            continue
        with open(filepath, 'r') as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    if not header_done:
                        header_lines.append(line)
                    continue
                if stripped.startswith('#'):
                    if not header_done:
                        header_lines.append(line)
                    continue
                # Separator lines
                if set(stripped.replace(' ', '')) == {'-'}:
                    if not header_done:
                        header_lines.append(line)
                    continue
                # Column header
                if stripped.startswith('OBSID'):
                    if not header_done:
                        header_lines.append(line)
                        # Next separator
                        header_done = True
                    continue
                # Data line
                data_lines.append(line)

        # Clean up chunk file
        os.remove(filepath)

    # Sort data lines by OBSID (first column)
    data_lines.sort(key=lambda l: l.split()[0] if l.split() else '')

    # Write merged file
    with open(os.path.join(BASE_DIR, output_file), 'w') as f:
        for hl in header_lines:
            f.write(hl)
        for dl in data_lines:
            f.write(dl)
        # Final separator
        if header_lines:
            for hl in header_lines:
                if set(hl.strip().replace(' ', '')) == {'-'}:
                    f.write(hl)
                    break

    return len(data_lines)


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
    parser.add_argument('--model', type=str, default='simple',
                        choices=['absorbed', 'simple'])
    parser.add_argument('--defgamma', type=float, default=2.0)
    parser.add_argument('--mincounts', type=int, default=40)
    parser.add_argument('--mingamma', type=int, default=200)
    parser.add_argument('--emin', type=float, default=0.3)
    parser.add_argument('--emax', type=float, default=10.0)
    parser.add_argument('--bkg', type=str, default='subtract',
                        choices=['subtract', 'none'])
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

    # Build pass-through arguments
    passthrough = ['--nh', str(args.nh), '--model', args.model,
                   '--defgamma', str(args.defgamma),
                   '--mincounts', str(args.mincounts),
                   '--mingamma', str(args.mingamma),
                   '--emin', str(args.emin),
                   '--emax', str(args.emax),
                   '--bkg', args.bkg]
    if args.redshift is not None:
        passthrough.extend(['--redshift', str(args.redshift)])
    if args.caldb:
        passthrough.extend(['--caldb', args.caldb])

    # Read all entries and split
    all_entries = []  # (mode, line)

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

    total = len(all_entries)
    n_chunks = min(args.nproc, total)
    chunk_size = (total + n_chunks - 1) // n_chunks

    print(f"Total: {total} observations → {n_chunks} chunks "
          f"of ~{chunk_size}")

    if args.dryrun:
        print("\n[DRY RUN] No fitting performed.")
        return

    # Split into chunks, each chunk may have both PC and WT entries
    jobs = []
    for i in range(n_chunks):
        chunk = all_entries[i*chunk_size:(i+1)*chunk_size]
        if not chunk:
            continue

        # Separate PC and WT entries
        pc_lines = [(hdr, line) for mode, line, hdr in chunk
                     if mode == 'pc']
        wt_lines = [(hdr, line) for mode, line, hdr in chunk
                     if mode == 'wt']

        mini_tables = {}
        if pc_lines:
            mpath = os.path.join(BASE_DIR,
                                  f'.pc_fit_chunk_{i:02d}.txt')
            write_mini_table(pc_lines[0][0],
                              [l for _, l in pc_lines], mpath)
            mini_tables['pc'] = mpath
        if wt_lines:
            mpath = os.path.join(BASE_DIR,
                                  f'.wt_fit_chunk_{i:02d}.txt')
            write_mini_table(wt_lines[0][0],
                              [l for _, l in wt_lines], mpath)
            mini_tables['wt'] = mpath

        # Set mode for this chunk
        chunk_args = list(passthrough)
        if mini_tables.keys() == {'pc'}:
            chunk_args.extend(['--modes', 'pc'])
        elif mini_tables.keys() == {'wt'}:
            chunk_args.extend(['--modes', 'wt'])
        else:
            chunk_args.extend(['--modes', 'both'])

        jobs.append((i, mini_tables, chunk_args))

    print(f"\nLaunching {len(jobs)} fitting workers...", flush=True)

    # Capture the current environment so it propagates correctly
    # to subprocesses spawned by workers. CIAO and HEASoft set
    # many environment variables (ASCDS_CALIB, CALDB, HEADAS,
    # LD_LIBRARY_PATH, etc.) that must be present for Sherpa
    # and XSPEC model libraries to function.
    parent_env = os.environ.copy()

    results = []
    with ProcessPoolExecutor(max_workers=args.nproc) as executor:
        futures = {}
        for cid, mtabs, cargs in jobs:
            future = executor.submit(
                run_fit_chunk, cid, mtabs, cargs, BASE_DIR,
                parent_env)
            futures[future] = cid

        print(f"  All {len(futures)} workers submitted. "
              f"Waiting for results...\n", flush=True)

        for future in as_completed(futures):
            cid = futures[future]
            result = future.result()
            results.append(result)
            status = 'OK' if result['returncode'] == 0 else 'FAIL'
            print(f"  Chunk {cid:02d}: {status} "
                  f"[{len(results)}/{len(futures)} done]",
                  flush=True)
            if result['returncode'] != 0 and result['stderr_tail']:
                print(f"    {result['stderr_tail'][:200]}",
                      flush=True)

    # Clean up mini-tables
    for _, mtabs, _ in jobs:
        for mpath in mtabs.values():
            if os.path.exists(mpath):
                os.remove(mpath)

    # Merge results
    result_files = [r['results_file'] for r in results
                    if r['results_file'] and r['returncode'] == 0]

    if result_files:
        n_merged = merge_results(result_files, args.output)
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


if __name__ == '__main__':
    main()

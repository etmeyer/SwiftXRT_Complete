#!/usr/bin/env python3
"""
parallel_extract.py

Run swift_xrt_extract_spectra.py in parallel by splitting
the master tables into chunks and processing each chunk in
a separate subprocess. Each worker operates in a temporary
working directory to avoid xselect session file conflicts.

Usage:
    python parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 16
    python parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 32 --mode pc
    python parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 8 --dryrun

This script:
  1. Reads the master table(s)
  2. Splits entries into --nproc chunks
  3. Writes a temporary mini-table for each chunk
  4. Launches parallel instances of the extraction script
  5. Collects and reports results

Requirements:
    The extraction script and all its dependencies (HEASoft, CALDB)
    must be available in the current environment.
"""

import os
import sys
import re
import argparse
import subprocess
import tempfile
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed


BASE_DIR = os.path.abspath(os.getcwd())

# Directory where this script (and sibling scripts) live.
# This allows running from the data directory while the
# scripts are installed elsewhere (e.g., on PATH).
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def read_master_table_raw(table_file):
    """Read a master table, returning header line and included entries."""
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
            # Check include column
            quote_match = re.search(r'"([^"]*)"', stripped)
            before_quote = stripped[:quote_match.start()].strip() \
                if quote_match else stripped
            parts = before_quote.split()
            if len(parts) >= 3 and parts[2].lower() == 'yes':
                entries.append(line)
    return header, entries


def write_mini_table(header, entries, filepath):
    """Write a subset of master table entries to a file."""
    with open(filepath, 'w') as f:
        if header:
            f.write(header)
        for entry in entries:
            f.write(entry)


def run_chunk(chunk_id, mini_table, script_args, base_dir, env):
    """
    Run the extraction script on a chunk in a temporary
    working directory.

    We create a temp dir, symlink all OBSID directories into it,
    copy the mini-table, and run from there. This avoids xselect
    session file conflicts between parallel workers.

    env is the parent process's os.environ, passed explicitly
    to ensure HEASoft/CALDB environment variables propagate
    correctly through ProcessPoolExecutor workers.
    """
    # Create a temporary working directory
    tmp_dir = tempfile.mkdtemp(
        prefix=f'xrt_chunk{chunk_id:02d}_', dir=base_dir)

    try:
        # Symlink all OBSID directories into the temp dir
        # (so the script can find event files without copying)
        obsid_pattern = re.compile(r'^\d{8,11}$')
        for d in os.listdir(base_dir):
            full = os.path.join(base_dir, d)
            if os.path.isdir(full) and obsid_pattern.match(d):
                link = os.path.join(tmp_dir, d)
                if not os.path.exists(link):
                    os.symlink(full, link)

        # Also symlink pileup_overrides.txt if it exists
        overrides = os.path.join(base_dir, 'pileup_overrides.txt')
        if os.path.exists(overrides):
            os.symlink(overrides,
                        os.path.join(tmp_dir, 'pileup_overrides.txt'))

        # Copy the mini-table
        table_basename = os.path.basename(mini_table)
        shutil.copy2(mini_table, os.path.join(tmp_dir, table_basename))

        # Build the command
        cmd = [sys.executable,
               os.path.join(SCRIPT_DIR, 'swift_xrt_extract_spectra.py')]
        cmd.extend(script_args)

        # Override the table argument
        if '--pctable' in cmd or '--wttable' in cmd:
            # Already specified — the mini table name matches
            pass
        else:
            # Determine which table type this is from the filename
            if 'pc_chunk' in table_basename:
                cmd.extend(['--pctable', table_basename,
                            '--mode', 'pc'])
            elif 'wt_chunk' in table_basename:
                cmd.extend(['--wttable', table_basename,
                            '--mode', 'wt'])

        # Run from the temp directory
        result = subprocess.run(
            cmd, cwd=tmp_dir,
            capture_output=True, text=True, timeout=7200,
            env=env)

        # Since we symlinked OBSID dirs, outputs are written
        # directly into the real directories.

        return {
            'chunk_id': chunk_id,
            'returncode': result.returncode,
            'stdout_tail': result.stdout[-500:] if result.stdout else '',
            'stderr_tail': result.stderr[-500:] if result.stderr else '',
        }

    except Exception as e:
        return {
            'chunk_id': chunk_id,
            'returncode': -1,
            'stdout_tail': '',
            'stderr_tail': str(e),
        }
    finally:
        # Clean up temp directory (symlinks only, no real data)
        shutil.rmtree(tmp_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description='Run spectral extraction in parallel.')
    parser.add_argument('--nproc', type=int, default=8,
                        help='Number of parallel workers (default: 8)')
    parser.add_argument('--mode', type=str, default='both',
                        choices=['pc', 'wt', 'both'],
                        help='Which mode(s) to extract')
    parser.add_argument('--pctable', type=str,
                        default='pc_master_table.txt')
    parser.add_argument('--wttable', type=str,
                        default='wt_master_table.txt')
    parser.add_argument('--dryrun', action='store_true',
                        help='Show chunk splitting without running')

    # All remaining arguments are passed through to the
    # extraction script
    parser.add_argument('--ra', type=float, required=True)
    parser.add_argument('--dec', type=float, required=True)
    parser.add_argument('--rout', type=str, default='47')
    parser.add_argument('--mincounts', type=int, default=20)
    parser.add_argument('--bkg-inner', type=float, default=100.0)
    parser.add_argument('--bkg-outer', type=float, default=160.0)
    parser.add_argument('--wt-srcrad', type=int, default=20)
    parser.add_argument('--wt-bkginner', type=int, default=80)
    parser.add_argument('--wt-bkgouter', type=int, default=120)

    args = parser.parse_args()

    # Build pass-through arguments for the extraction script
    passthrough = [
        '--ra', str(args.ra),
        '--dec', str(args.dec),
        '--rout', args.rout,
        '--mincounts', str(args.mincounts),
        '--bkg-inner', str(args.bkg_inner),
        '--bkg-outer', str(args.bkg_outer),
        '--wt-srcrad', str(args.wt_srcrad),
        '--wt-bkginner', str(args.wt_bkginner),
        '--wt-bkgouter', str(args.wt_bkgouter),
    ]

    # Read master tables and split into chunks
    jobs = []  # list of (chunk_id, mini_table_path, mode)

    chunk_id = 0

    if args.mode in ('pc', 'both'):
        pc_path = os.path.join(BASE_DIR, args.pctable)
        if os.path.exists(pc_path):
            header, entries = read_master_table_raw(pc_path)
            n_pc = len(entries)
            # Split into chunks
            n_chunks = min(args.nproc, n_pc)
            if n_chunks > 0:
                chunk_size = (n_pc + n_chunks - 1) // n_chunks
                for i in range(n_chunks):
                    chunk_entries = entries[i*chunk_size:
                                            (i+1)*chunk_size]
                    if not chunk_entries:
                        continue
                    mini_path = os.path.join(
                        BASE_DIR,
                        f'.pc_chunk_{chunk_id:02d}.txt')
                    write_mini_table(header, chunk_entries, mini_path)
                    jobs.append((chunk_id, mini_path, 'pc'))
                    chunk_id += 1

                print(f"PC: {n_pc} observations → {n_chunks} chunks")
        else:
            if args.mode == 'pc':
                print(f"ERROR: {args.pctable} not found.")
                sys.exit(1)

    if args.mode in ('wt', 'both'):
        wt_path = os.path.join(BASE_DIR, args.wttable)
        if os.path.exists(wt_path):
            header, entries = read_master_table_raw(wt_path)
            n_wt = len(entries)
            n_chunks = min(args.nproc, n_wt)
            if n_chunks > 0:
                chunk_size = (n_wt + n_chunks - 1) // n_chunks
                for i in range(n_chunks):
                    chunk_entries = entries[i*chunk_size:
                                            (i+1)*chunk_size]
                    if not chunk_entries:
                        continue
                    mini_path = os.path.join(
                        BASE_DIR,
                        f'.wt_chunk_{chunk_id:02d}.txt')
                    write_mini_table(header, chunk_entries, mini_path)
                    jobs.append((chunk_id, mini_path, 'wt'))
                    chunk_id += 1

                print(f"WT: {n_wt} observations → {n_chunks} chunks")
        else:
            if args.mode == 'wt':
                print(f"ERROR: {args.wttable} not found.")
                sys.exit(1)

    if not jobs:
        print("No observations to process.")
        sys.exit(0)

    print(f"Total: {len(jobs)} chunks across {args.nproc} workers")

    if args.dryrun:
        print("\n[DRY RUN] Chunk contents:")
        for cid, mpath, mode in jobs:
            with open(mpath, 'r') as f:
                n_lines = sum(1 for l in f
                              if l.strip() and
                              not l.strip().startswith('OBSID'))
            print(f"  Chunk {cid:02d} [{mode.upper()}]: "
                  f"{n_lines} observations")
            os.remove(mpath)
        print("\nNo extraction performed (--dryrun).")
        return

    # Run in parallel
    print(f"\nLaunching {len(jobs)} extraction workers...", flush=True)

    # Capture the current environment so it propagates correctly
    # to subprocesses. HEASoft tools (xselect, xrtmkarf, etc.)
    # need HEADAS, CALDB, LD_LIBRARY_PATH, etc.
    parent_env = os.environ.copy()

    results = []
    with ProcessPoolExecutor(max_workers=args.nproc) as executor:
        futures = {}
        for cid, mpath, mode in jobs:
            # Build mode-specific args
            chunk_args = list(passthrough)
            if mode == 'pc':
                chunk_args.extend([
                    '--pctable', os.path.basename(mpath),
                    '--mode', 'pc'])
            else:
                chunk_args.extend([
                    '--wttable', os.path.basename(mpath),
                    '--mode', 'wt'])

            future = executor.submit(
                run_chunk, cid, mpath, chunk_args, BASE_DIR,
                parent_env)
            futures[future] = (cid, mode)

        print(f"  All {len(futures)} workers submitted. "
              f"Waiting for results...\n", flush=True)

        for future in as_completed(futures):
            cid, mode = futures[future]
            result = future.result()
            results.append(result)
            status = 'OK' if result['returncode'] == 0 else 'FAIL'
            print(f"  Chunk {cid:02d} [{mode.upper()}]: {status} "
                  f"[{len(results)}/{len(futures)} done]",
                  flush=True)
            if result['returncode'] != 0 and result['stderr_tail']:
                print(f"    {result['stderr_tail'][:200]}")

    # Clean up mini-tables
    for _, mpath, _ in jobs:
        if os.path.exists(mpath):
            os.remove(mpath)

    # Summary
    n_ok = sum(1 for r in results if r['returncode'] == 0)
    n_fail = sum(1 for r in results if r['returncode'] != 0)
    print(f"\nDone. {n_ok} chunks succeeded, {n_fail} failed.")

    if n_fail > 0:
        print("Failed chunks:")
        for r in results:
            if r['returncode'] != 0:
                print(f"  Chunk {r['chunk_id']:02d}: "
                      f"exit code {r['returncode']}")


if __name__ == '__main__':
    main()

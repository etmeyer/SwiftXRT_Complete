#!/usr/bin/env python3
"""
xrt_pipeline.py

Run xrtpipeline on Swift XRT observations to produce cleaned
level-2 event files, exposure maps, and auxiliary products.

Supports single-OBSID, sequential batch, and parallel batch modes.

Usage:
    # Single OBSID
    python xrt_pipeline.py --indir /path/to/00035017001 \
        --outdir /path/to/output/00035017001 \
        --ra 187.2779 --dec 2.0524

    # Batch: all OBSIDs under input directory
    python xrt_pipeline.py --batch \
        --indir /path/to/XRT_input \
        --outdir /path/to/XRT_output \
        --ra 187.2779 --dec 2.0524

    # Parallel batch
    python xrt_pipeline.py --batch --nproc 16 \
        --indir /path/to/XRT_input \
        --outdir /path/to/XRT_output \
        --ra 187.2779 --dec 2.0524

Requirements:
    HEASoft (xrtpipeline must be in PATH)
    CALDB initialized

Author: Eileen T. Meyer
"""

import os
import re
import sys
import argparse
import subprocess
import time
from pathlib import Path
from typing import Optional, Union, Dict
from concurrent.futures import ProcessPoolExecutor, as_completed


# ---------------------------------------------------------------
# ObsID detection
# ---------------------------------------------------------------

def get_obs_id(data_path: Path) -> Optional[str]:
    """
    Extract the ObsID from the observation directory.

    Strategy (in order of priority):
      1. The directory name itself IS the ObsID (8-11 digit number)
      2. Infer from files matching sw<obsid>* inside the tree
    """
    dir_name = data_path.name

    # 1. Directory name is the ObsID (8-11 digits, zero-padded)
    if re.fullmatch(r'\d{8,11}', dir_name):
        return dir_name

    # 2. Infer from files inside the directory
    for pattern in ('**/sw*.img.gz', '**/sw*_cl.evt*',
                    '**/sw*.evt*'):
        matches = list(data_path.glob(pattern))
        if matches:
            stem = matches[0].name
            m = re.match(r'sw(\d{8,11})', stem)
            if m:
                return m.group(1)

    return None


# ---------------------------------------------------------------
# Product verification
# ---------------------------------------------------------------

def verify_level2_products(output_path: Path, obs_id: str) -> dict:
    """
    Check the output directory for key Level-2 products.
    Uses wildcards for window modes since the XRT auto-selects
    w1-w4 based on count rate.
    """
    expected = {
        'cleaned_pc_evt':  list(output_path.glob(
            f'**/*{obs_id}xpc*po_cl.evt*')),
        'cleaned_wt_evt':  list(output_path.glob(
            f'**/*{obs_id}xwt*po_cl.evt*')),
        'exposure_map_pc': list(output_path.glob(
            f'**/*{obs_id}xpc*_ex.img*')),
        'attitude_file':   list(output_path.glob(
            f'**/*{obs_id}*pat.fits*')),
        'hk_file':         list(output_path.glob(
            f'**/*{obs_id}*xhd.hk*')),
    }
    return {k: ([f.name for f in v] if v else None)
            for k, v in expected.items()}


# ---------------------------------------------------------------
# Run xrtpipeline on a single OBSID
# ---------------------------------------------------------------

def run_pipeline(
    data_path: Union[Path, str],
    output_path: Union[Path, str],
    srcra: float,
    srcdec: float,
    createexpomap: str = 'yes',
    cleanup: str = 'no',
    clobber: str = 'yes',
    exprpcgrade: str = '0-12',
    exprwtgrade: str = '0-2',
    exprpdgrade: str = '0-2',
    logdir: Optional[Union[Path, str]] = None,
    quiet: bool = False,
    env: Optional[dict] = None,
) -> dict:
    """
    Run xrtpipeline on a single Swift XRT ObsID directory.

    Parameters
    ----------
    data_path     : Path to the raw ObsID input directory
    output_path   : Path to the output directory
    srcra/srcdec  : Source coordinates in decimal degrees
    createexpomap : Create exposure map? ('yes'/'no')
    cleanup       : Remove intermediate files? ('yes'/'no')
    clobber       : Overwrite existing output? ('yes'/'no')
    exprpcgrade   : PC-mode grade selection
    exprwtgrade   : WT-mode grade selection
    exprpdgrade   : PD-mode grade selection
    logdir        : Directory for log files (default: output_path)
    quiet         : Suppress terminal output (for parallel mode)

    Returns
    -------
    dict with keys: obs_id, success, products, elapsed_s, log_file
    """
    data_path = Path(data_path).resolve()
    output_path = Path(output_path).resolve()

    result = {
        'obs_id': None,
        'success': False,
        'products': {},
        'elapsed_s': 0,
        'log_file': None,
    }

    # Validate input
    if not data_path.exists():
        if not quiet:
            print(f"[ERROR] Input directory not found: {data_path}")
        return result

    # Extract ObsID
    obs_id = get_obs_id(data_path)
    if obs_id is None:
        if not quiet:
            print(f"[ERROR] Could not determine ObsID from "
                  f"{data_path}")
        return result

    result['obs_id'] = obs_id
    stem_inputs = f"sw{obs_id}"

    if not quiet:
        print(f"\n{'='*60}")
        print(f"  ObsID: {obs_id}")
        print(f"  Input:  {data_path}")
        print(f"  Output: {output_path}")
        print(f"{'='*60}")

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Log file
    log_path = Path(logdir) if logdir else output_path
    log_file = log_path / f'xrtpipeline_{obs_id}.log'
    result['log_file'] = str(log_file)

    # Build xrtpipeline command
    xrt_params = dict(
        indir=str(data_path),
        outdir=str(output_path),
        steminputs=stem_inputs,
        srcra=srcra,
        srcdec=srcdec,
        createexpomap=createexpomap,
        cleanup=cleanup,
        clobber=clobber,
        exprpcgrade=exprpcgrade,
        exprwtgrade=exprwtgrade,
        exprpdgrade=exprpdgrade,
    )

    cmd = ["xrtpipeline"] + [f"{k}={v}" for k, v in xrt_params.items()]
    cmd_str = " ".join(cmd)

    if not quiet:
        print(f"  [CMD] {cmd_str}")

    # Run xrtpipeline
    # CRITICAL: Use cwd= instead of os.chdir() — os.chdir is
    # process-wide and NOT safe for parallel execution. The cwd=
    # parameter sets the working directory only for the subprocess.
    t0 = time.time()

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(data_path),  # safe per-subprocess working dir
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout per OBSID
            env=env,  # None = inherit parent env; explicit dict for parallel
        )

        returncode = proc.returncode
        result['success'] = (returncode == 0)
        result['elapsed_s'] = time.time() - t0

        # Write log file
        with open(log_file, 'w') as f:
            f.write(f"# xrtpipeline log for ObsID {obs_id}\n")
            f.write(f"# Command: {cmd_str}\n")
            f.write(f"# Return code: {returncode}\n")
            f.write(f"# Elapsed: {result['elapsed_s']:.1f} s\n\n")
            if proc.stdout:
                f.write("=== STDOUT ===\n")
                f.write(proc.stdout)
            if proc.stderr:
                f.write("\n=== STDERR ===\n")
                f.write(proc.stderr)

        if not quiet:
            if returncode == 0:
                print(f"  [SUCCESS] {obs_id} "
                      f"({result['elapsed_s']:.0f}s)")
            else:
                print(f"  [FAILED]  {obs_id} exit code {returncode}")
                # Print last few lines of stderr for debugging
                if proc.stderr:
                    lines = proc.stderr.strip().split('\n')
                    for line in lines[-5:]:
                        print(f"    {line}")

    except subprocess.TimeoutExpired:
        result['elapsed_s'] = time.time() - t0
        if not quiet:
            print(f"  [TIMEOUT] {obs_id} exceeded 1 hour limit")
        with open(log_file, 'w') as f:
            f.write(f"# xrtpipeline TIMEOUT for ObsID {obs_id}\n")
            f.write(f"# Command: {cmd_str}\n")
        return result

    except Exception as exc:
        result['elapsed_s'] = time.time() - t0
        if not quiet:
            print(f"  [ERROR] {obs_id}: {exc}")
        return result

    # Verify products
    products = verify_level2_products(output_path, obs_id)
    result['products'] = products

    if not quiet:
        for product_type, filenames in products.items():
            if filenames:
                for fn in filenames:
                    print(f"  [FOUND]   {product_type:<20s} {fn}")
            else:
                print(f"  [MISSING] {product_type}")

    return result


# ---------------------------------------------------------------
# Worker function for parallel execution
# ---------------------------------------------------------------

def _parallel_worker(args_tuple):
    """
    Wrapper for ProcessPoolExecutor. Unpacks arguments and
    calls run_pipeline with quiet=True.
    """
    obs_dir, out_dir, srcra, srcdec, kwargs, env = args_tuple
    return run_pipeline(
        data_path=obs_dir,
        output_path=out_dir,
        srcra=srcra,
        srcdec=srcdec,
        quiet=True,
        env=env,
        **kwargs,
    )


# ---------------------------------------------------------------
# Batch runner (sequential or parallel)
# ---------------------------------------------------------------

def batch_run_pipeline(
    root_input_dir: Union[Path, str],
    root_output_dir: Union[Path, str],
    srcra: float,
    srcdec: float,
    nproc: int = 1,
    **pipeline_kwargs,
) -> Dict[str, dict]:
    """
    Run xrtpipeline for every ObsID subdirectory found under
    root_input_dir.

    Parameters
    ----------
    root_input_dir  : Parent directory containing ObsID subdirs
    root_output_dir : Parent output directory
    srcra/srcdec    : Source coordinates
    nproc           : Number of parallel workers (1 = sequential)
    **pipeline_kwargs : Extra args forwarded to run_pipeline()

    Returns
    -------
    dict mapping ObsID string -> result dict
    """
    root_input_dir = Path(root_input_dir).resolve()
    root_output_dir = Path(root_output_dir).resolve()

    # Collect ObsID subdirectories (8-11 digit names)
    obs_dirs = sorted([
        d for d in root_input_dir.iterdir()
        if d.is_dir() and re.fullmatch(r'\d{8,11}', d.name)
    ])

    if not obs_dirs:
        print(f"[WARNING] No ObsID directories found under "
              f"{root_input_dir}")
        return {}

    n_obs = len(obs_dirs)
    print(f"\nFound {n_obs} ObsID "
          f"director{'y' if n_obs == 1 else 'ies'} to process.")

    if nproc > 1:
        print(f"Running in parallel with {nproc} workers.",
              flush=True)
    else:
        print(f"Running sequentially.\n")

    # Build job list
    # Capture the environment explicitly for parallel workers
    parent_env = os.environ.copy()

    jobs = []
    for obs_dir in obs_dirs:
        obs_id = obs_dir.name
        out_dir = root_output_dir / obs_id
        jobs.append((obs_dir, out_dir, srcra, srcdec,
                      pipeline_kwargs, parent_env))

    results = {}
    t0_batch = time.time()

    if nproc <= 1:
        # Sequential mode
        for i, (obs_dir, out_dir, ra, dec, kwargs, _env) in \
                enumerate(jobs, 1):
            obs_id = obs_dir.name
            print(f"[{i}/{n_obs}] {obs_id}")
            r = run_pipeline(
                data_path=obs_dir, output_path=out_dir,
                srcra=ra, srcdec=dec, **kwargs)
            results[obs_id] = r
    else:
        # Parallel mode
        with ProcessPoolExecutor(max_workers=nproc) as executor:
            future_map = {}
            for job in jobs:
                obs_id = job[0].name
                future = executor.submit(_parallel_worker, job)
                future_map[future] = obs_id

            print(f"  All {len(future_map)} workers submitted. "
                  f"Waiting for results...\n", flush=True)

            completed = 0
            for future in as_completed(future_map):
                obs_id = future_map[future]
                completed += 1
                try:
                    r = future.result()
                    results[obs_id] = r
                    status = 'OK' if r['success'] else 'FAIL'
                    elapsed = f"{r['elapsed_s']:.0f}s"
                    print(f"  [{completed}/{n_obs}] {obs_id}: "
                          f"{status} ({elapsed})", flush=True)
                except Exception as exc:
                    results[obs_id] = {
                        'obs_id': obs_id, 'success': False,
                        'products': {}, 'elapsed_s': 0,
                        'log_file': None,
                    }
                    print(f"  [{completed}/{n_obs}] {obs_id}: "
                          f"EXCEPTION: {exc}", flush=True)

    batch_elapsed = time.time() - t0_batch

    # Summary
    passed = [k for k, v in results.items() if v['success']]
    failed = [k for k, v in results.items() if not v['success']]

    print(f"\n{'='*60}")
    print(f"  BATCH SUMMARY  ({n_obs} observations, "
          f"{batch_elapsed:.0f}s total)")
    print(f"{'='*60}")
    print(f"  SUCCESS : {len(passed)}")
    for obs in sorted(passed):
        t = results[obs]['elapsed_s']
        print(f"            {obs}  ({t:.0f}s)")
    if failed:
        print(f"  FAILED  : {len(failed)}")
        for obs in sorted(failed):
            print(f"            {obs}")
            log = results[obs].get('log_file')
            if log:
                print(f"              log: {log}")
    print(f"{'='*60}\n")

    return results


# ---------------------------------------------------------------
# Main (CLI)
# ---------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Run xrtpipeline on Swift XRT observations.')
    parser.add_argument('--indir', type=str, required=True,
                        help='Input directory (single OBSID or '
                             'parent with --batch)')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Output directory (single OBSID or '
                             'parent with --batch)')
    parser.add_argument('--ra', type=float, required=True,
                        help='Source RA in decimal degrees')
    parser.add_argument('--dec', type=float, required=True,
                        help='Source Dec in decimal degrees')
    parser.add_argument('--batch', action='store_true',
                        help='Process all OBSID subdirs under '
                             '--indir')
    parser.add_argument('--nproc', type=int, default=1,
                        help='Number of parallel workers '
                             '(default: 1 = sequential)')
    parser.add_argument('--createexpomap', type=str, default='yes',
                        choices=['yes', 'no'],
                        help='Create exposure maps (default: yes)')
    parser.add_argument('--cleanup', type=str, default='no',
                        choices=['yes', 'no'],
                        help='Remove intermediate files '
                             '(default: no)')
    parser.add_argument('--clobber', type=str, default='yes',
                        choices=['yes', 'no'],
                        help='Overwrite existing output '
                             '(default: yes)')
    args = parser.parse_args()

    kwargs = dict(
        createexpomap=args.createexpomap,
        cleanup=args.cleanup,
        clobber=args.clobber,
    )

    if args.batch:
        batch_run_pipeline(
            root_input_dir=args.indir,
            root_output_dir=args.outdir,
            srcra=args.ra,
            srcdec=args.dec,
            nproc=args.nproc,
            **kwargs,
        )
    else:
        result = run_pipeline(
            data_path=args.indir,
            output_path=args.outdir,
            srcra=args.ra,
            srcdec=args.dec,
            **kwargs,
        )
        sys.exit(0 if result['success'] else 1)


if __name__ == '__main__':
    main()

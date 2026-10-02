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
import pty
import signal
import shutil
import argparse
import tempfile
import subprocess
import time
from pathlib import Path
from typing import Optional, Union, Dict, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

from swift_xrt_env import require_heasoft_shell


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

def verify_level2_products(output_path: Path, input_path: Path,
                           obs_id: str) -> dict:
    """
    Check the OUTPUT directory for the cleaned Level-2 products that
    downstream tools actually consume: the cleaned event file
    (``*_cl.evt``) and exposure map (``*_ex.img``) for each observing
    mode. Window modes (w1-w4) are wild-carded since the XRT
    auto-selects them based on count rate.

    The check is *mode-aware*: only the modes actually present in the
    INPUT observation (PC and/or WT) are verified, so a PC-only or
    WT-only observation no longer reports the absent mode as
    ``[MISSING]``.

    NB: the attitude file (``sw<OBSID>sat.fits.gz``) and housekeeping
    files (``xrt/hk/*.hk``) are pipeline *inputs* that live in the
    input tree, not products that xrtpipeline writes to ``output_path``.
    The previous version globbed for them in ``output_path`` and so
    logged a bogus ``[MISSING] attitude_file`` / ``[MISSING] hk_file``
    on every successful run. They are dropped here: xrtpipeline aborts
    with a non-zero exit if either input is absent, so a post-run check
    for them was both misplaced and redundant.
    """
    # Detect which modes the observation actually contains from the
    # raw input event files (e.g. sw<obsid>xpcw3po*.evt*,
    # sw<obsid>xwtw2po*.evt*).
    modes = []
    if list(input_path.glob(f'**/sw{obs_id}xpc*.evt*')):
        modes.append('pc')
    if list(input_path.glob(f'**/sw{obs_id}xwt*.evt*')):
        modes.append('wt')

    expected = {}
    for m in modes:
        expected[f'cleaned_{m}_evt'] = list(output_path.glob(
            f'**/*{obs_id}x{m}*po_cl.evt*'))
        expected[f'exposure_map_{m}'] = list(output_path.glob(
            f'**/*{obs_id}x{m}*po_ex.img*'))

    return {k: ([f.name for f in v] if v else None)
            for k, v in expected.items()}


# ---------------------------------------------------------------
# Subprocess execution (headless-safe, with timeout + clean kill)
# ---------------------------------------------------------------

def _run_xrtpipeline(
    cmd, cwd: str, env: dict, timeout: int,
) -> Tuple[Optional[int], str, str, bool]:
    """
    Run ``cmd`` with a private pseudo-terminal as its controlling
    terminal and a hard wall-clock timeout.

    Why a pty (controlling-terminal fix)
    ------------------------------------
    Several HEASoft tasks (xrthkproc, xrtfilter, ... anything that goes
    through headas_stdio.c) call ``open("/dev/tty")`` at start-up to
    redirect prompts. With no controlling terminal -- cron, nohup,
    systemd, a detached SSH session -- that open fails with ENXIO
    ("No such device or address") and the task aborts:

        ERROR: No such device or address
        Task xrthkproc 0.0 terminating with status 6
        Unable to redirect prompts to the /dev/tty (headas_stdio.c:152)

    ``mode=h`` does NOT fix this: the open happens regardless of mode,
    and xrtpipeline's own (Perl) argument parser rejects ``mode=h``
    outright. The robust fix is to give the child a real controlling
    terminal. We allocate a pty, and in the child (post-fork, pre-exec)
    start a new session and open the pty slave as session leader so it
    becomes the controlling terminal -- exactly what ``script -qec``
    did as the Task-1.1 band-aid, but native and per-process.

    The pty is used ONLY as the controlling terminal. Real stdout and
    stderr are captured on ordinary pipes so the caller can log them
    and surface them on failure (no stream swallowing). stdin is
    /dev/null so any tool that ignores the pty and reads stdin gets EOF
    rather than blocking.

    Timeout / clean kill
    --------------------
    ``os.setsid()`` makes the child a process-group leader, so on
    timeout we ``killpg(SIGKILL)`` the whole group and every descendant
    it spawned (xrtproducts, xselect, ...) dies too -- not just the
    xrtpipeline driver.

    Returns ``(returncode, stdout, stderr, timed_out)``.
    """
    master_fd, slave_fd = pty.openpty()
    slave_name = os.ttyname(slave_fd)

    def _preexec():
        # New session (drops any inherited controlling tty), then open
        # the pty slave as the session leader so it becomes this
        # process's controlling terminal and open("/dev/tty") succeeds.
        os.setsid()
        fd = os.open(slave_name, os.O_RDWR)
        os.close(fd)

    proc = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        preexec_fn=_preexec,
    )
    # The parent does not use the slave end. Keep the master open for
    # the child's lifetime -- closing it would tear down the pty and
    # break any further open("/dev/tty") in the child.
    os.close(slave_fd)

    timed_out = False
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        stdout, stderr = proc.communicate()
    finally:
        os.close(master_fd)

    return proc.returncode, stdout, stderr, timed_out


def _write_run_log(log_file, obs_id, cmd_str, returncode,
                   elapsed_s, stdout, stderr, reason=None):
    """Write the full xrtpipeline stdout+stderr verbatim to log_file."""
    with open(log_file, 'w') as f:
        f.write(f"# xrtpipeline log for ObsID {obs_id}\n")
        f.write(f"# Command: {cmd_str}\n")
        if reason:
            f.write(f"# Result: FAILED ({reason})\n")
        f.write(f"# Return code: {returncode}\n")
        f.write(f"# Elapsed: {elapsed_s:.1f} s\n\n")
        if stdout:
            f.write("=== STDOUT ===\n")
            f.write(stdout)
        if stderr:
            f.write("\n=== STDERR ===\n")
            f.write(stderr)


# ---------------------------------------------------------------
# Run xrtpipeline on a single OBSID
# ---------------------------------------------------------------

def run_pipeline(
    data_path: Union[Path, str],
    output_path: Union[Path, str],
    srcra: float,
    srcdec: float,
    createexpomap: str = 'yes',
    extractproducts: str = 'no',
    cleanup: str = 'no',
    clobber: str = 'yes',
    exprpcgrade: str = '0-12',
    exprwtgrade: str = '0-2',
    exprpdgrade: str = '0-2',
    timeout: int = 600,
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
    extractproducts : Run xrtproducts to extract source spectra/light
                    curves? ('yes'/'no'). Default 'no': downstream
                    (parallel_extract.py) regenerates all source
                    products from the cleaned event files, so the
                    per-OBSID xrtproducts step is wasted work -- and it
                    is where the xselect hang on certain OBSIDs occurs.
    cleanup       : Remove intermediate files? ('yes'/'no')
    clobber       : Overwrite existing output? ('yes'/'no')
    exprpcgrade   : PC-mode grade selection
    exprwtgrade   : WT-mode grade selection
    exprpdgrade   : PD-mode grade selection
    timeout       : Per-OBSID wall-clock limit in seconds. On timeout
                    the whole xrtpipeline process group is killed and
                    the OBSID is recorded as FAILED so the batch can
                    continue instead of one stuck worker hanging it.
    logdir        : Directory for log files (default: output_path)
    quiet         : Suppress terminal output (for parallel mode)

    Returns
    -------
    dict with keys: obs_id, success, products, elapsed_s, log_file,
    reason, stderr_tail
    """
    data_path = Path(data_path).resolve()
    output_path = Path(output_path).resolve()

    result = {
        'obs_id': None,
        'success': False,
        'products': {},
        'elapsed_s': 0,
        'log_file': None,
        'reason': None,       # why it failed (timeout / exit code / ...)
        'stderr_tail': None,  # tail of the real diagnostic on failure
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
        extractproducts=extractproducts,
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

    # Per-OBSID private, writable PFILES directory.
    # WHY: HEASoft FTOOLs read/write their parameter (.par) files under
    # the first entry of $PFILES. The shared default ($HOME/pfiles) is
    # easily polluted -- on amorgos, CIAO's bundled prefilter learned a
    # non-existent CIAO leapsec path into $HOME/pfiles/prefilter.par, so
    # xrtpipeline's prefilter step aborts with
    #     couldn't get parameter 'leapname' [file not found ...]
    #     (PIL_BAD_FILE_ACCESS)
    # even though $HEADAS/refdata/leapsec.fits is perfectly readable.
    # A fresh private pfiles dir prepended to $HEADAS/syspfiles forces
    # every tool to copy the pristine system defaults, dodging the
    # pollution -- and makes parallel workers race-free (no shared .par).
    run_env = (env or os.environ).copy()
    headas = run_env.get('HEADAS', '')
    pfiles_dir = tempfile.mkdtemp(prefix=f'xrtpipe_pf_{obs_id}_')
    run_env['PFILES'] = (f"{pfiles_dir};{headas}/syspfiles"
                         if headas else pfiles_dir)

    # Run xrtpipeline.
    # CRITICAL: pass cwd= rather than os.chdir() -- os.chdir is
    # process-wide and unsafe for parallel execution; cwd= sets the
    # working directory for this subprocess only.
    t0 = time.time()
    returncode, stdout, stderr, timed_out = None, '', '', False
    run_error = None
    try:
        returncode, stdout, stderr, timed_out = _run_xrtpipeline(
            cmd, cwd=str(data_path), env=run_env, timeout=timeout)
    except Exception as exc:
        run_error = exc
    finally:
        shutil.rmtree(pfiles_dir, ignore_errors=True)

    result['elapsed_s'] = time.time() - t0

    # Classify the outcome (loud, explicit -- never a silent skip).
    if run_error is not None:
        result['reason'] = f"wrapper error: {run_error}"
        stderr = stderr or str(run_error)
    elif timed_out:
        result['reason'] = f"timeout after {timeout} seconds"
    elif returncode != 0:
        result['reason'] = f"exit code {returncode}"
    result['success'] = result['reason'] is None

    # Always write the full log (stdout + stderr verbatim, never
    # swallowed -- swallowing is what hid the headless leapname error).
    _write_run_log(log_file, obs_id, cmd_str, returncode,
                   result['elapsed_s'], stdout, stderr,
                   reason=result['reason'])

    if not result['success']:
        # Surface the real diagnostic instead of burying it in the log.
        diag = stderr.strip() or stdout.strip()
        result['stderr_tail'] = "\n".join(diag.splitlines()[-15:])
        if not quiet:
            print(f"  [FAILED]  {obs_id} ({result['reason']})")
            if result['stderr_tail']:
                print("    --- xrtpipeline output (tail) ---",
                      file=sys.stderr)
                for line in result['stderr_tail'].splitlines():
                    print(f"    {line}", file=sys.stderr)
            print(f"    (full log: {log_file})")
        return result

    if not quiet:
        print(f"  [SUCCESS] {obs_id} ({result['elapsed_s']:.0f}s)")

    # Verify products (only meaningful on a successful run).
    products = verify_level2_products(output_path, data_path, obs_id)
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
                    elapsed = f"{r['elapsed_s']:.0f}s"
                    if r['success']:
                        print(f"  [{completed}/{n_obs}] {obs_id}: "
                              f"OK ({elapsed})", flush=True)
                    else:
                        print(f"  [{completed}/{n_obs}] {obs_id}: "
                              f"FAILED ({elapsed}) -- "
                              f"{r.get('reason') or 'unknown'}",
                              flush=True)
                except Exception as exc:
                    results[obs_id] = {
                        'obs_id': obs_id, 'success': False,
                        'products': {}, 'elapsed_s': 0,
                        'log_file': None,
                        'reason': f"worker exception: {exc}",
                        'stderr_tail': None,
                    }
                    print(f"  [{completed}/{n_obs}] {obs_id}: "
                          f"FAILED -- worker exception: {exc}",
                          flush=True)

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
            reason = results[obs].get('reason') or 'unknown'
            print(f"    FAILED: {obs} ({reason})")
            tail = results[obs].get('stderr_tail')
            if tail:
                # Last couple of lines of the real diagnostic, so the
                # cause is visible in the summary, not just the log.
                for line in tail.splitlines()[-2:]:
                    print(f"            | {line}")
            log = results[obs].get('log_file')
            if log:
                print(f"            log: {log}")
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
    parser.add_argument('--extractproducts', type=str, default='no',
                        choices=['yes', 'no'],
                        help='Run xrtproducts to extract source '
                             'spectra/light curves (default: no). '
                             'Downstream regenerates these from the '
                             'cleaned event files, and this step is '
                             'where the xselect hang occurs.')
    parser.add_argument('--cleanup', type=str, default='no',
                        choices=['yes', 'no'],
                        help='Remove intermediate files '
                             '(default: no)')
    parser.add_argument('--clobber', type=str, default='yes',
                        choices=['yes', 'no'],
                        help='Overwrite existing output '
                             '(default: yes)')
    parser.add_argument('--timeout', type=int, default=600,
                        help='Per-OBSID wall-clock limit in seconds '
                             '(default: 600). On timeout the stuck '
                             'xrtpipeline process group is killed, the '
                             'OBSID is marked FAILED, and the batch '
                             'continues.')
    args = parser.parse_args()

    # xrtpipeline inherits our HEADAS/CALDB; refuse to start in a
    # terminal where CIAO has replaced them (see swift_xrt_env.py).
    require_heasoft_shell(['xrtpipeline'])

    kwargs = dict(
        createexpomap=args.createexpomap,
        extractproducts=args.extractproducts,
        cleanup=args.cleanup,
        clobber=args.clobber,
        timeout=args.timeout,
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

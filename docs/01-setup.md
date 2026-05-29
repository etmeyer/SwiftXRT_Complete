# Step 1 — Setup

Install the pipeline scripts, make the two external analysis environments
(HEASoft and CIAO) available, and verify the whole thing with
`swift_xrt_doctor.py`. This is the foundation every later step assumes: the
download and reduction steps shell out to HEASoft FTOOLS, and the fit step
runs inside CIAO's bundled Python (the only supported home for Sherpa). Get a
clean doctor run here and the rest of the workflow has what it needs.

## What runs

### Requirements

The pipeline shells out to **two** separate analysis environments — HEASoft
and CIAO — so both need to be available before you start. `swift_xrt_doctor.py`
(see [How it works](#how-it-works)) verifies all of this in one shot.

**Operating environment:**
- Python 3.9+
- [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) (tested with 6.36) — provides `xrtpipeline`, `xrtmkarf`, `grppha`, `xselect`, `ftlist`, used by the download, reduction, and extraction steps.
- [CIAO](https://cxc.cfa.harvard.edu/ciao/) (tested with 4.16) — used **only** for [Sherpa](https://sherpa.readthedocs.io/), the fitting engine in [Step 8](08-fit-and-plot.md). CIAO's bundled Python is the only supported home for Sherpa.

**CALDB:**
- [HEASoft Swift CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/) — the calibration tree containing the Swift XRT response files at `data/swift/xrt/cpf/rmf`.
- Sourcing CIAO repoints `$CALDB` at *its* Chandra calibration tree. If you run the fit step from a CIAO shell, you'll need `--caldb /path/to/heasoft/caldb` so the Swift RMFs can be found; see the [`$CALDB` and CIAO](#gotchas) gotcha below.

**Python packages:**
- `astropy`
- `numpy`
- `scipy`
- `matplotlib`
- `requests` (for data download only)
- `astroquery` (optional, fallback name resolver)

`scipy` must be importable from **whichever Python ends up running the
King-profile step ([Step 5](05-pc-inspection.md)) and the fit step
([Step 8](08-fit-and-plot.md))**. Because the fit step needs Sherpa, and Sherpa
lives only in CIAO's bundled Python, the recommended route is to let that one
Python run both — which means installing scipy into the CIAO environment once
with `conda install -p /opt/ciao/ciao-4.16 scipy` (see below). A single CIAO
shell can then run the entire pipeline.

> Sherpa can also be obtained via `pip install sherpa` into a standalone
> environment, but that is not the tested/recommended route here — CIAO is. If
> you go the pip route, you are responsible for making `scipy` and the other
> packages importable in that same environment.

**Data:** Swift XRT observations downloaded from the HEASARC archive (see
[Step 2 — Download data](02-download.md)), organized as OBSID subdirectories.
The download and pipeline scripts handle this automatically.

### Install the pipeline

Clone or download the repository, place the scripts somewhere permanent, and
make them executable:

```bash
# Example: install to /opt/swift-xrt-pipeline
sudo cp -r swift-xrt-pipeline /opt/swift-xrt-pipeline
sudo chmod +x /opt/swift-xrt-pipeline/*.py
```

### Shell setup

On a shared machine, put the environment setup in `/etc/bash.bashrc.local` so
every user picks it up, rather than in a single user's `~/.bashrc`. The block
below defines `setup_swiftxrt` (which manages only the pipeline's own `PATH`)
and sources HEASoft and the HEASoft Swift CALDB:

```bash
# /etc/bash.bashrc.local — Swift XRT pipeline site-wide setup
setup_swiftxrt() {
    local pipedir="/opt/swift-xrt-pipeline"
    case ":$PATH:" in *":$pipedir:"*) ;; *) export PATH="$pipedir:$PATH" ;; esac
    command -v swift_xrt_summary.py >/dev/null && \
        echo "Swift XRT Pipeline ready ($(ls $pipedir/*.py 2>/dev/null | wc -l) scripts in $pipedir)"
}

# HEASoft (adjust path / version to your install)
export HEADAS=/opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39
source $HEADAS/headas-init.sh

# Swift CALDB (HEASoft side, NOT CIAO's)
source /opt/CALDB/software/tools/caldbinit.sh

# Pipeline on PATH
setup_swiftxrt

# CIAO is intentionally NOT sourced here — users source it when they need
# Sherpa (the fit step), e.g. by running `ciao` (alias to /opt/ciao/.../bin/ciao.sh).
```

`setup_swiftxrt` only puts the pipeline scripts on `PATH` — it does **not**
touch HEASoft, CIAO, or CALDB. Those are the separate `source` lines above,
sourced by the user or admin.

### One-time scipy into CIAO

The recommended Python for the pipeline is CIAO's bundled one, because that is
where Sherpa lives. CIAO does not ship `scipy`, which the King-profile step
([Step 5](05-pc-inspection.md)) and a few other utilities need, so install it
into the CIAO environment once:

```bash
# Make CIAO's bundled python ship scipy too (one-time)
conda install -p /opt/ciao/ciao-4.16 scipy
# (substitute the prefix path for your CIAO install location)
```

With this done, a single CIAO shell has HEASoft (via the shell setup), Sherpa
(via CIAO), and scipy — enough to run the entire pipeline end to end without
juggling two shells.

## How it works

### Verify with the doctor

`swift_xrt_doctor.py` runs a green/yellow/red checklist over everything the
pipeline assumes about the environment and exits non-zero if anything fails:

```bash
swift_xrt_doctor.py             # full output, colored if on a terminal
swift_xrt_doctor.py --quiet     # only failures printed
swift_xrt_doctor.py --no-color  # plain output for logs
# Exit code 0 only if every check passes — usable in CI / cron preambles.
```

It checks that:

- the pipeline scripts are on `PATH` (`setup_swiftxrt`),
- HEASoft is loaded (`$HEADAS` set) and the FTOOLS resolve,
- the CALDB environment variables are set (and warns if `$CALDB` points at CIAO's),
- the Swift XRT response files are present under `$CALDB`,
- Sherpa is importable,
- the required Python packages are importable,
- there is free disk at `/opt`.

After completing the install, shell setup, and scipy step, a healthy run from a
CIAO shell with the HEASoft CALDB in front looks like this (the lone `[WARN]`
is the optional `astroquery` resolver):

```
[OK] Pipeline on PATH: swift_xrt_summary.py -> /opt/swift-xrt-pipeline/swift_xrt_summary.py
[OK] HEASoft loaded ($HEADAS set, all FTOOLS on PATH)
       xrtpipeline  /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/xrtpipeline
       xrtmkarf     /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/xrtmkarf
       grppha       /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/grppha
       xselect      /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/xselect
       ftlist       /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/ftlist
[OK] HEASoft version 6.36
[OK] CALDB configured
       $CALDB=/opt/CALDB
       $CALDBCONFIG=/opt/ciao/ciao-4.16/CALDB/software/tools/caldb.config
[OK] Swift XRT response files present under $CALDB
[OK] Python 3.11.6 (/opt/ciao/ciao-4.16/binexe/python3.11)
[OK] astropy      7.2.0
[OK] numpy        1.26.2
[OK] scipy        1.17.1
[OK] matplotlib   3.8.2
[OK] requests     2.31.0
[WARN] astroquery not installed (optional; download script falls back to SIMBAD/NED/Sesame)
[OK] sherpa       4.16.0
[OK] Disk free at /opt: 26.5 GB

14 checks: 13 ok, 1 warn, 0 fail
```

The single most useful diagnostic is the **yellow** CALDB clash: you're in a
CIAO shell and `$CALDB` is still pointing at CIAO's Chandra tree instead of
HEASoft's Swift CALDB:

```
[OK] Pipeline on PATH: swift_xrt_summary.py -> /opt/swift-xrt-pipeline/swift_xrt_summary.py
[OK] HEASoft loaded ($HEADAS set, all FTOOLS on PATH)
       xrtpipeline  /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/xrtpipeline
       xrtmkarf     /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/xrtmkarf
       grppha       /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/grppha
       xselect      /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/xselect
       ftlist       /opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39/bin/ftlist
[OK] HEASoft version 6.36
[WARN] $CALDB points inside a CIAO install (/opt/ciao/ciao-4.16/CALDB)
       This is the Chandra CALDB, not HEASoft's Swift CALDB.
       Pass --caldb /opt/CALDB to parallel_fit.py / swift_xrt_fit_spectra.py.
[WARN] Swift XRT RMFs found at /opt/CALDB but not under $CALDB
       Point $CALDB at /opt/CALDB or pass --caldb /opt/CALDB.
[OK] Python 3.11.6 (/opt/ciao/ciao-4.16/binexe/python3.11)
[OK] astropy      7.2.0
[OK] numpy        1.26.2
[OK] scipy        1.17.1
[OK] matplotlib   3.8.2
[OK] requests     2.31.0
[WARN] astroquery not installed (optional; download script falls back to SIMBAD/NED/Sesame)
[OK] sherpa       4.16.0
[OK] Disk free at /opt: 26.5 GB

14 checks: 11 ok, 3 warn, 0 fail
```

In that case, pass `--caldb /opt/CALDB` to the fit step (or re-source the
HEASoft CALDB ahead of CIAO's). If you see any `[FAIL]` lines, the doctor
prints exactly what to run to fix each one.

## Inputs and outputs

**Inputs:** none — this step takes nothing from a prior stage. It consumes only
the external installs you point it at: the HEASoft tree (`$HEADAS`), the
HEASoft Swift CALDB (`/opt/CALDB`), and the CIAO install (`/opt/ciao/...`).

**Outputs:** a working environment, not files. Concretely:
- the pipeline scripts on `PATH` (via `setup_swiftxrt`),
- `$HEADAS`, `$CALDB`, `$CALDBCONFIG` exported,
- scipy importable from CIAO's Python,
- a green `swift_xrt_doctor.py` run (exit code 0).

Everything downstream — [Step 2](02-download.md) onward — assumes this state.

## Common variants

```bash
# Full colored checklist (default)
swift_xrt_doctor.py

# Only show what's broken
swift_xrt_doctor.py --quiet

# Plain text for logs / CI
swift_xrt_doctor.py --no-color

# Single-user install (no root): put the /etc/bash.bashrc.local block in
# ~/.bashrc or ~/.bash_profile instead. On a multi-user machine, a ~/.bashrc
# install sets up the environment for only the one user who did it.
```

## Gotchas

- **`$CALDB` and CIAO coexistence.** Sourcing CIAO repoints `$CALDB` at its own
  Chandra calibration tree. If you run the fit step from a CIAO shell, the Swift
  XRT RMFs won't be found under `$CALDB`. Either re-source the HEASoft CALDB
  ahead of CIAO's, or pass `--caldb /opt/CALDB` to `swift_xrt_fit_spectra.py` /
  `parallel_fit.py`. This is the most common yellow-flag the doctor reports, and
  it is cross-referenced from [Step 2 — Download](02-download.md#gotchas) and the
  fit step.
- **`setup_swiftxrt` only touches `PATH`.** It does not source HEASoft, CIAO, or
  CALDB — those are separate `source` lines in `/etc/bash.bashrc.local`. If the
  doctor reports HEASoft missing even though `setup_swiftxrt` ran, you skipped
  (or mis-ordered) the `source $HEADAS/headas-init.sh` line.
- **Use the conda `-p` (prefix) form, not `-n` (name).** Install scipy into CIAO
  with `conda install -p /opt/ciao/ciao-4.16 scipy`, pointing at the install
  prefix. The named-environment form (`-n`) targets a different conda env and
  won't land scipy where CIAO's Python actually looks.
- **Sherpa lives only in CIAO's Python.** The fit step has no supported home
  outside CIAO's bundled interpreter. That's why scipy is installed *into* CIAO
  rather than into a separate venv — so one Python runs both the King-profile
  step and the fit step.
- **Site-wide vs. single-user.** On a shared box, install into
  `/etc/bash.bashrc.local` so every user inherits the environment. A `~/.bashrc`
  install only sets things up for the one user who did it.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_

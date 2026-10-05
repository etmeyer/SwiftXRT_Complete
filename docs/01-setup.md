# Step 1 — Setup

Summary: Here you will find instructions on installing the pipeline scripts,
making the two external analysis environments (HEASoft and CIAO) available, and
verifying them with `swift_xrt_doctor.py`. The pipeline runs in **two
terminals**: a HEASoft terminal for download, reduction, inspection and
extraction, and a CIAO terminal for the fit, because Sherpa (the fitting
engine) lives only in CIAO's Python, and CIAO breaks the HEASoft tools. Read
[Two terminals](#two-terminals) before you start.

## What runs

### Requirements

The pipeline shells out to **two** separate analysis environments — HEASoft
and CIAO — so both need to be installed before you start. `swift_xrt_doctor.py`
(see [How it works](#how-it-works)) checks each terminal.

**Operating environment:**
- Python 3.9+
- [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) (tested with 6.36) — provides `xrtpipeline`, `xrtmkarf`, `grppha`, `xselect`, `ftlist`, used by the reduction and extraction steps.
- [CIAO](https://cxc.cfa.harvard.edu/ciao/) (tested with 4.16) — used **only** for [Sherpa](https://sherpa.readthedocs.io/), the fitting engine in [Step 8](08-fit-and-plot.md). CIAO's bundled Python is the only supported home for Sherpa.

**CALDB:**
- [HEASoft Swift CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/) — the calibration tree containing the Swift XRT response files at `data/swift/xrt/cpf/rmf`.
- Setting up CIAO repoints `$CALDB` at *its* Chandra calibration tree, which has no Swift files. The fit step therefore takes `--caldb /path/to/heasoft/caldb`; see [Gotchas](#gotchas).

**Python packages:**
- `astropy`
- `numpy`
- `scipy`
- `matplotlib`
- `requests` (for data download only)
- `astroquery` (optional, fallback name resolver)

Each terminal runs the scripts with whatever `python3` comes first on its
`PATH`, so these packages must be importable there. On amorgos the HEASoft
terminal's Python (the `heasoft` conda env) has everything except Sherpa, and CIAO's Python has everything including Sherpa.

> Sherpa can also be obtained via `pip install sherpa` into a standalone
> environment, but that is not the tested route here — CIAO is.

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

On a shared machine, define the setup commands in `/etc/bash.bashrc.local` so
every user gets them. On amorgos that file defines three commands, which you
**run yourself** in each new terminal — nothing is set up automatically at
login:

```bash
# /etc/bash.bashrc.local (excerpt, amorgos)

# Puts the pipeline scripts on PATH. Touches nothing else.
setup_swiftxrt() {
    local pipedir="/opt/swift-xrt-pipeline"
    case ":$PATH:" in *":$pipedir:"*) ;; *) export PATH="$pipedir:$PATH" ;; esac
    command -v swift_xrt_summary.py >/dev/null && \
        echo "Swift XRT Pipeline ready ($(ls $pipedir/*.py 2>/dev/null | wc -l) scripts in $pipedir)"
}

# HEASoft + the HEASoft (Swift) CALDB.
heainit() {
    conda activate heasoft
    export HEADAS=/opt/heasoft/heasoft-6.36/x86_64-pc-linux-gnu-libc2.39
    . $HEADAS/headas-init.sh
    source /opt/CALDB/software/tools/caldbinit.sh
}

# CIAO (Sherpa).
alias ciao='source /opt/ciao/ciao-4.16/bin/ciao.sh'
```

These only exist in interactive terminals (`/etc/bash.bashrc` skips the file
otherwise), so a shell script that calls `heainit` or `ciao` must set up the
environment itself, e.g. by sourcing `headas-init.sh` / `ciao.sh` directly.

### Two terminals

| Terminal | Set up with | Runs |
| -------- | ----------- | ---- |
| **HEASoft** | `setup_swiftxrt; heainit` | Steps 2–7: download, xrtpipeline, survey, inspection, extraction |
| **CIAO** | `setup_swiftxrt; ciao` | Step 8: fit, with `--caldb /opt/CALDB` |

```bash
# Terminal 1 — HEASoft. Never run `ciao` in this terminal.
setup_swiftxrt
heainit
swift_xrt_download.py ...      # Step 2
xrt_pipeline.py ...            # Steps 3-7
...

# Terminal 2 — CIAO, for the fit only.
setup_swiftxrt
ciao
swift_xrt_fit_spectra.py ... --caldb /opt/CALDB     # Step 8
```

Why two:

- **CIAO breaks the HEASoft tools.** Once you run `ciao`, `python3` is CIAO's
  wrapper script, which resets `$HEADAS` to CIAO's own `spectral` directory
  and `$CALDB` to the Chandra CALDB *inside every pipeline script*. Every
  HEASoft tool the script launches then fails, e.g.
  `Failed to open /opt/ciao/ciao-4.16/spectral/bin//xselect.mdb`. Running
  `heainit` again does not undo it; open a new terminal.
- **Without CIAO there is no Sherpa.** In the HEASoft terminal the fit step
  stops with `ERROR: sherpa is not importable from this Python`.

Steps 2 and 4–6 are plain Python and also work in the CIAO terminal.

The scripts check this themselves: the HEASoft steps refuse to start in a CIAO
terminal, and the fit refuses to start without Sherpa or without Swift
responses in its CALDB, each saying which terminal to use.

### One-time scipy into CIAO

CIAO does not ship `scipy`, which the inspection steps (4–6) need. It is only
required if you want to run those steps from the CIAO terminal, but it is
cheap; on amorgos it is already done:

```bash
conda install -p /opt/ciao/ciao-4.16 scipy
# (substitute the prefix path for your CIAO install location)
```

## How it works

### Verify with the doctor

`swift_xrt_doctor.py` runs a checklist over everything the pipeline assumes
about the current terminal, then lists which steps that terminal can run. Run
it once in each terminal:

```bash
swift_xrt_doctor.py             # full output, colored if on a terminal
swift_xrt_doctor.py --quiet     # only failures printed
swift_xrt_doctor.py --no-color  # plain output for logs
# Exit code 0 only if no check FAILs -- usable in CI / cron preambles.
```

It checks that:

- the pipeline scripts are on `PATH` (`setup_swiftxrt`),
- HEASoft is loaded (`$HEADAS` set) and the FTOOLS resolve — or reports that
  CIAO is set up and this is a fit-only terminal,
- the CALDB environment variables are set (and warns if `$CALDB` points at CIAO's),
- the Swift XRT response files are present under `$CALDB`,
- the Python packages are importable (and Sherpa, for the fit),
- there is free disk at `/opt`.

A healthy **HEASoft terminal** (`setup_swiftxrt; heainit`) on amorgos. The
Sherpa warning is expected: Sherpa belongs to the CIAO terminal.

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
       $CALDBCONFIG=/opt/CALDB/software/tools/caldb.config
[OK] Swift XRT response files present under $CALDB
[OK] Python 3.12.13 (/opt/anaconda3/envs/heasoft/bin/python3)
[OK] astropy      7.2.0
[OK] numpy        2.4.2
[OK] scipy        1.15.2
[OK] matplotlib   3.10.8
[OK] requests     2.34.2
[WARN] astroquery not installed (optional; download script falls back to SIMBAD/NED/Sesame)
[WARN] sherpa not importable -- not needed until the fit step (Step 8), which runs in a CIAO terminal: setup_swiftxrt; ciao
[OK] Disk free at /opt: 14.8 GB

This terminal can run:
  Step 2     download                 yes
  Steps 3, 7 xrtpipeline, extraction  yes
  Steps 4-6  survey, inspection       yes
  Step 8     fit                      no -- needs CIAO: setup_swiftxrt; ciao

14 checks: 12 ok, 2 warn, 0 fail
```

A healthy **CIAO terminal** (`setup_swiftxrt; ciao`). The CALDB warnings are
why the fit takes `--caldb /opt/CALDB`:

```
[OK] Pipeline on PATH: swift_xrt_summary.py -> /opt/swift-xrt-pipeline/swift_xrt_summary.py
[WARN] CIAO is set up in this terminal: use it for the fit step (Step 8) only
       CIAO's python replaced $HEADAS with /opt/ciao/ciao-4.16/spectral,
       so xrtpipeline and extraction fail here. Run those in
       a separate terminal: setup_swiftxrt; heainit (no ciao).
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
[OK] Disk free at /opt: 14.8 GB

This terminal can run:
  Step 2     download                 yes
  Steps 3, 7 xrtpipeline, extraction  no -- needs HEASoft without CIAO: setup_swiftxrt; heainit
  Steps 4-6  survey, inspection       yes
  Step 8     fit                      yes (CALDB: --caldb /opt/CALDB)

13 checks: 9 ok, 4 warn, 0 fail
```

If you see any `[FAIL]` lines, the doctor prints what to run to fix each one.

## Inputs and outputs

**Inputs:** none — this step takes nothing from a prior stage. It consumes only
the external installs you point it at: the HEASoft tree (`$HEADAS`), the
HEASoft Swift CALDB (`/opt/CALDB`), and the CIAO install (`/opt/ciao/...`).

**Outputs:** two working terminals, not files:
- a HEASoft terminal: pipeline on `PATH`, `$HEADAS`, `$CALDB` (Swift),
  `$CALDBCONFIG` exported, no CIAO;
- a CIAO terminal: pipeline on `PATH`, CIAO's Python with Sherpa;
- a doctor run with no `[FAIL]` in each.

Everything downstream — [Step 2](02-download.md) onward — assumes this state.

## Common variants

```bash
# Full colored checklist (default)
swift_xrt_doctor.py

# Only show what's broken
swift_xrt_doctor.py --quiet

# Plain text for logs / CI
swift_xrt_doctor.py --no-color

# Single-user install (no root): put the /etc/bash.bashrc.local definitions
# in ~/.bashrc instead. On a multi-user machine, a ~/.bashrc install sets up
# the environment for only the one user who did it.
```

## Gotchas

- **Never run `ciao` in the HEASoft terminal.** CIAO's Python resets `$HEADAS`
  and `$CALDB` inside every pipeline script, so `xrt_pipeline.py` and the
  extraction step fail (they now refuse to start and tell you so). Running
  `heainit` afterwards does not undo it — CIAO's `python3` stays first on
  `PATH` and `$CALDB` stays CIAO's. Open a new terminal.
- **The fit needs `--caldb`.** In the CIAO terminal `$CALDB` is the Chandra
  CALDB with no Swift files. Pass `--caldb /opt/CALDB` to
  `swift_xrt_fit_spectra.py` / `parallel_fit.py`; without it they stop with an
  error naming the problem.
- **`setup_swiftxrt` only touches `PATH`.** It does not set up HEASoft, CIAO,
  or CALDB — that is `heainit` / `ciao`.
- **Use the conda `-p` (prefix) form, not `-n` (name), for CIAO.** CIAO is
  registered by path, so `conda install -n ciao-4.16 ...` silently does
  nothing; use `conda install -p /opt/ciao/ciao-4.16 ...`.
- **Keep the path to `XRT_output` under about 150 characters.** Beyond
  that, xselect cannot name its working files and Step 7 fails
  ([Step 7, gotcha 7](07-extract.md#gotchas)).
- **Site-wide vs. single-user.** On a shared box, define the commands in
  `/etc/bash.bashrc.local` so every user inherits them. A `~/.bashrc` install
  only sets things up for the one user who did it.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_

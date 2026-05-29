# Swift XRT Spectral Analysis Pipeline

A Python-based pipeline for automated spectral extraction, fitting, and light curve generation from Swift X-Ray Telescope (XRT) data in both Photon Counting (PC) and Windowed Timing (WT) modes.

Developed for analysis of point sources (e.g., blazars, AGN) observed across multiple epochs. The pipeline handles pile-up assessment, optimal extraction region determination, proper BACKSCAL correction for WT mode, and batch spectral fitting with Sherpa.

**Author:** Eileen T. Meyer ([@etmeyer](https://github.com/etmeyer))  
**License:** MIT

---

## Requirements

The pipeline shells out to **two** separate analysis environments — HEASoft and
CIAO — so both need to be available before you start. `swift_xrt_doctor.py`
(see Installation, Step 4) verifies all of this in one shot.

**Operating environment:**
- Python 3.9+
- [HEASoft](https://heasarc.gsfc.nasa.gov/lheasoft/) (tested with 6.36) — provides `xrtpipeline`, `xrtmkarf`, `grppha`, `xselect`, `ftlist`, used by the download, reduction, and extraction steps.
- [CIAO](https://cxc.cfa.harvard.edu/ciao/) (tested with 4.16) — used **only** for [Sherpa](https://sherpa.readthedocs.io/), the fitting engine in Step 6. CIAO's bundled Python is the only supported home for Sherpa.

**CALDB:**
- [HEASoft Swift CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/) — the calibration tree containing the Swift XRT response files at `data/swift/xrt/cpf/rmf`.
- Sourcing CIAO repoints `$CALDB` at *its* Chandra calibration tree. If you run the fit step (Step 6) from a CIAO shell, you'll need `--caldb /path/to/heasoft/caldb` so the Swift RMFs can be found; see Installation Step 2 and the [CALDB and CIAO coexistence](#caldb-and-ciao-coexistence) note below.

**Python packages:**
- `astropy`
- `numpy`
- `scipy`
- `matplotlib`
- `requests` (for data download only)
- `astroquery` (optional, fallback name resolver)

`scipy` must be importable from **whichever Python ends up running the King-profile step (Step 2a) and the fit step (Step 6)**. Because the fit step needs Sherpa, and Sherpa lives only in CIAO's bundled Python, the recommended route is to let that one Python run both — which means installing scipy into the CIAO environment once with `conda install -n ciao-4.16 scipy` (see Installation, Step 3). A single CIAO shell can then run the entire pipeline.

> Sherpa can also be obtained via `pip install sherpa` into a standalone environment, but that is not the tested/recommended route here — CIAO is. If you go the pip route, you are responsible for making `scipy` and the other packages importable in that same environment.

**Data:** Swift XRT observations downloaded from the HEASARC archive (see Step 0 below), organized as OBSID subdirectories. The download and pipeline scripts handle this automatically.

---

### Step 1 — Place the pipeline

Clone or download the repository, place the scripts somewhere permanent, and make them executable:

```bash
# Example: install to /opt/swift-xrt-pipeline
sudo cp -r swift-xrt-pipeline /opt/swift-xrt-pipeline
sudo chmod +x /opt/swift-xrt-pipeline/*.py
```

### Step 2 — Site-wide shell setup (multi-user boxes)

On a shared machine, put the environment setup in `/etc/bash.bashrc.local` so every user picks it up, rather than in a single user's `~/.bashrc`. The block below defines `setup_swiftxrt` (which manages only the pipeline's own `PATH`) and sources HEASoft and the HEASoft Swift CALDB:

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

`setup_swiftxrt` only puts the pipeline scripts on `PATH` — it does **not** touch HEASoft, CIAO, or CALDB. Those are the separate `source` lines above, sourced by the user or admin.

> **Single-user install:** If you lack root, or you're the only user, put the same block in your `~/.bashrc` (or `~/.bash_profile`) instead. On a multi-user machine, a `~/.bashrc` install sets up the environment for only the one user who did it.

If you run the fit step from a CIAO shell, see the [CALDB and CIAO coexistence](#caldb-and-ciao-coexistence) note below for `--caldb`.

### Step 3 — One-time scipy into CIAO

The recommended Python for the pipeline is CIAO's bundled one, because that is where Sherpa lives. CIAO does not ship `scipy`, which the King-profile step (Step 2a) and a few other utilities need, so install it into the CIAO environment once:

```bash
# Make CIAO's bundled python ship scipy too (one-time)
conda install -n ciao-4.16 scipy
# (substitute the conda env name for your CIAO version)
```

With this done, a single CIAO shell has HEASoft (via Step 2), Sherpa (via CIAO), and scipy — enough to run the entire pipeline end to end without juggling two shells.

### Step 4 — Verify with the doctor

`swift_xrt_doctor.py` runs a green/yellow/red checklist over everything the pipeline assumes about the environment and exits non-zero if anything fails:

```bash
swift_xrt_doctor.py             # full output, colored if on a terminal
swift_xrt_doctor.py --quiet     # only failures printed
swift_xrt_doctor.py --no-color  # plain output for logs
# Exit code 0 only if every check passes — usable in CI / cron preambles.
```

It checks that:

- the pipeline scripts are on `PATH` (Step 1 / `setup_swiftxrt`),
- HEASoft is loaded (`$HEADAS` set) and the FTOOLS resolve,
- the CALDB environment variables are set (and warns if `$CALDB` points at CIAO's),
- the Swift XRT response files are present under `$CALDB`,
- Sherpa is importable,
- the required Python packages are importable,
- there is free disk at `/opt`.

After completing Steps 1–3, a healthy run from a CIAO shell with the HEASoft CALDB in front looks like this (the lone `[WARN]` is the optional `astroquery` resolver):

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
[OK] Python 3.11.6 (/opt/ciao/ciao-4.16/binexe/python3.11)
[OK] astropy      7.2.0
[OK] numpy        1.26.2
[OK] scipy        1.11.4
[OK] matplotlib   3.8.2
[OK] requests     2.31.0
[WARN] astroquery not installed (optional; download script falls back to SIMBAD/NED/Sesame)
[OK] sherpa       4.16.0
[OK] Disk free at /opt: 26.9 GB

14 checks: 13 ok, 1 warn, 0 fail
```

The single most useful diagnostic is the **yellow** CALDB clash: you're in a CIAO shell and `$CALDB` is still pointing at CIAO's Chandra tree instead of HEASoft's Swift CALDB:

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
[OK] scipy        1.11.4
[OK] matplotlib   3.8.2
[OK] requests     2.31.0
[WARN] astroquery not installed (optional; download script falls back to SIMBAD/NED/Sesame)
[OK] sherpa       4.16.0
[OK] Disk free at /opt: 26.9 GB

14 checks: 11 ok, 3 warn, 0 fail
```

In that case, pass `--caldb /opt/CALDB` to the fit step (or re-source the HEASoft CALDB ahead of CIAO's). If you see any `[FAIL]` lines, the doctor prints exactly what to run to fix each one.

---

## Workflow Overview

The pipeline is designed to be run step-by-step, with visual inspection at each stage. The recommended workflow is:

```
0. Download & reduce data
   a. Download from archive     →  swift_xrt_download.py
   b. Run xrtpipeline           →  xrt_pipeline.py
1. Survey observations          →  swift_xrt_summary.py
2. PC-mode inspection
   a. Pile-up / PSF analysis    →  swift_xrt_king_profile.py
   b. Source images              →  swift_pc_source_viewer.py
   c. Create PC master table    →  (shell command)
3. WT-mode inspection
   a. Profile viewer            →  swift_wt_summary_viewer.py
   b. Create WT master table    →  make_wt_master_table.py
4. Edit master tables           →  (manual review)
5. Extract spectra              →  swift_xrt_extract_spectra.py
6. Fit spectra & plot           →  swift_xrt_fit_spectra.py
7. Customize plot (optional)    →  plot_lightcurve.py
```

### Step 0: Download and reduce data

**0a.** Download raw observations from the HEASARC archive:

```bash
# By source name (resolves coordinates automatically)
swift_xrt_download.py --name "3C 273" --outdir XRT_input

# By coordinates
swift_xrt_download.py --ra 187.2779 --dec 2.0524 --radius 12 --outdir XRT_input

# Preview what's available without downloading
swift_xrt_download.py --name "3C 273" --list-only

# Download specific observations from a file
swift_xrt_download.py --obsid-file my_obsids.txt --outdir XRT_input
```

**0b.** Run the Swift XRT pipeline to produce cleaned level-2 event files:

```bash
# Process all OBSIDs (sequential)
xrt_pipeline.py --batch --indir XRT_input --outdir XRT_output \
    --ra 187.2779 --dec 2.0524

# Process in parallel (16 workers)
xrt_pipeline.py --batch --nproc 16 --indir XRT_input --outdir XRT_output \
    --ra 187.2779 --dec 2.0524
```

All subsequent steps are run from within the output directory:

```bash
cd XRT_output
```

### Step 1: Survey all observations

Get an overview of what data exists across all OBSIDs, including mode sequences, exposures, count rates, and orbit structure:

```bash
swift_xrt_summary.py              # detailed per-OBSID tables
swift_xrt_summary.py --compact    # one row per OBSID
```

This reveals which observations have PC data, WT data, or both, and flags potential pile-up.

### Step 2: PC-mode inspection

**2a.** Fit King profiles to assess pile-up and determine extraction radii:

```bash
swift_xrt_king_profile.py --ra 187.2779 --dec 2.0524
```

This produces per-observation diagnostic plots and `_pileup.txt` files with centroid positions and pile-up radii. For observations where the automated pile-up radius needs adjustment, create a `pileup_overrides.txt` file (see script details below).

**2b.** Generate zoomed source images for visual inspection:

```bash
swift_pc_source_viewer.py
```

Inspect the PDF for bad columns through the source, anomalous PSF shapes, or other issues.

**2c.** Create the PC master table:

```bash
find . -name '*xpc*po*_cl.evt' | sort | \
  awk -F'/' '{obsid=$2; file=$NF; gsub(/^\.\//, "", obsid); \
  sub(/_cl\.evt$/, "", file); \
  printf "%-14s %-35s %-10s %-10s \"%s\"\n", obsid, file, "yes", "no", ""}' | \
  (printf "%-14s %-35s %-10s %-10s %s\n" "OBSID" "filename" "include" "badstripe" "comment"; cat) \
  > pc_master_table.txt
```

### Step 3: WT-mode inspection

**3a.** Run the WT viewer to inspect strip profiles and extraction regions:

```bash
swift_wt_summary_viewer.py --ra 187.2779 --dec 2.0524
```

**3b.** Generate the WT master table:

```bash
make_wt_master_table.py
```

Observations under 20 seconds exposure are automatically set to `include=no`.

### Step 4: Edit master tables

Open `pc_master_table.txt` and `wt_master_table.txt` in a text editor. Set `include` to `no` for any observations you want to exclude (bad columns through source, anomalous data, etc.). Add notes in the comment field.

### Step 5: Extract spectra

```bash
# Initialize HEASoft and CALDB first, then:
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524

# Or one mode at a time:
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524 --mode pc
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524 --mode wt

# For large datasets, run in parallel:
parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 16
```

### Step 6: Fit spectra and generate light curve

```bash
# If running within CIAO (Sherpa), specify the HEASoft CALDB path:
swift_xrt_fit_spectra.py --nh 0.0179 --model simple \
    --caldb /path/to/heasoft/caldb

# With intrinsic absorption:
swift_xrt_fit_spectra.py --nh 0.0179 --redshift 0.158 \
    --caldb /path/to/heasoft/caldb

# For large datasets, run in parallel:
parallel_fit.py --nh 0.0179 --model simple --caldb /opt/CALDB --nproc 16
```

Output: `fit_results.txt` (table) and `flux_lightcurve.pdf` (νFν and Γ vs. time, with PC and WT points color-coded).

### Step 7 (optional): Customize the light curve plot

```bash
# Remake with custom axis limits
plot_lightcurve.py --ylim_flux 5e-12 5e-11 --ylim_gamma 1.2 2.2

# Zoom to a specific time range with reference lines
plot_lightcurve.py --xlim 2005 2024 \
    --times 2007.75 2015.5 2023.3 \
    --tlabels "Voltage change" "Flare A" "Flare B"
```

---

## Script Reference

### `swift_xrt_download.py`

Download Swift XRT observation data from the HEASARC archive. Supports source name resolution, coordinate-based queries, and direct OBSID specification. Automatically handles archive URL discovery (date-based HEASARC paths and UKSSDC mirror), resume of interrupted downloads, and XRT mode filtering.

```
Usage:
    swift_xrt_download.py --name "3C 273" --outdir XRT_input
    swift_xrt_download.py --ra 187.2779 --dec 2.0524 --outdir XRT_input
    swift_xrt_download.py --obsid 00035017001 --outdir XRT_input
    swift_xrt_download.py --obsid-file obsids.txt --outdir XRT_input

Input modes (mutually exclusive):
    --name        Resolve source name to coordinates, query catalog
    --ra/--dec    Query catalog at given position
    --obsid       Download a single observation by ID
    --obsid-file  Read observation IDs from file (one per line)

Key options:
    --radius      Search radius in arcmin (default: 12)
    --outdir      Output directory (default: current directory)
    --list-only   Show matching observations without downloading
    --max-obs     Maximum number of observations to download
    --test N      Download only N observations (for testing)
    --mode        Filter XRT modes: pc, wt, im (default: all)
    --clean-only  Skip unfiltered event files (*_uf.evt)
    --overwrite   Re-download existing files
    --products    Product types to download (default: xrt)

OBSID file format:
    # Comments and blank lines are ignored
    00035017001
    00035017002    # inline comments OK
    00050900011

Dependencies: requests (required), astroquery (optional)
```

Name resolution tries SIMBAD, NED, CDS Sesame, and astroquery in sequence. Archive URL discovery learns the optimal strategy (HEASARC date-based vs. UKSSDC flat) after the first successful download and reuses it for subsequent observations.

### `xrt_pipeline.py`

Run the HEASoft `xrtpipeline` task on raw Swift XRT data to produce cleaned level-2 event files, exposure maps, and auxiliary products. Supports single-OBSID, sequential batch, and parallel batch modes.

```
Usage:
    # Single OBSID
    xrt_pipeline.py --indir /path/to/00035017001 \
        --outdir /path/to/output/00035017001 --ra 187.2779 --dec 2.0524

    # Batch: all OBSIDs under input directory
    xrt_pipeline.py --batch --indir XRT_input --outdir XRT_output \
        --ra 187.2779 --dec 2.0524

    # Parallel batch
    xrt_pipeline.py --batch --nproc 16 --indir XRT_input \
        --outdir XRT_output --ra 187.2779 --dec 2.0524

Key options:
    --batch         Process all OBSID subdirs under --indir
    --nproc         Number of parallel workers (default: 1)
    --createexpomap Create exposure maps (default: yes)
    --cleanup       Remove intermediate files (default: no)
    --clobber       Overwrite existing output (default: yes)

Output (per OBSID):
    Cleaned event files for all modes (PC, WT) and observation
    types (slew, settling, pointed), exposure maps, attitude
    files, housekeeping files.

Requires: HEASoft (xrtpipeline), CALDB
```

Each `xrtpipeline` call processes everything within an OBSID — all modes and observation types — producing separate cleaned event files for each. Per-OBSID log files are written to the output directory. In parallel mode, the environment (HEASoft, CALDB paths) is explicitly propagated to worker processes.

### `swift_xrt_summary.py`

Crawl OBSID directories and report on all cleaned event files. Shows mode sequences (WT settling → WT pointed → PC pointed), exposures, count rates, GTI/orbit structure, and pile-up warnings.

```
Usage:
    swift_xrt_summary.py              # detailed tables
    swift_xrt_summary.py --compact    # one row per OBSID

Output (terminal only):
    Per-OBSID tables with mode, exposure, events, count rate
    GTI orbit analysis (segment durations, gap lengths, duty cycle)
    Pile-up warnings (PC >0.5 ct/s, WT >150 ct/s)
    Compact table with columns:
        OBSID, Date/Time, Total(ks), ct/s, Slew_i, Slew_f,
        N_WT, N_PC, WT_exp(ks), PC_exp(ks), Orb, Seq

Sequence codes:
    1=WT_SLEW  2=PC_SLEW  3=WT_SETTLING  4=PC_SETTLING
    5=WT_POINTED  6=PC_POINTED
```

### `swift_xrt_king_profile.py`

Fit King profiles to PC-mode radial surface brightness profiles to assess pile-up and determine source extraction regions.

```
Usage:
    swift_xrt_king_profile.py --ra <RA> --dec <DEC> [options]

Required:
    --ra        Source RA in degrees
    --dec       Source Dec in degrees

Key options:
    --rmin      Inner fit annulus radius in arcsec (default: 20)
    --rmax      Outer fit annulus radius in arcsec (default: 60)
    --rbin      Radial bin width in arcsec (default: 2)
    --sigma     Pile-up detection threshold, 2 consecutive bins (default: 3.0)
    --sigma2    Pile-up detection threshold, single bin (default: 4.0)
    --rc        King core radius, fixed (default: 5.8")
    --beta      King beta slope, fixed (default: 1.55)
    --pdf       Output PDF filename (default: king_profiles.pdf)

Output (per OBSID):
    {stem}_king_profile.png   - diagnostic plot
    {stem}_pileup.txt         - centroid, plate scale, pile-up radius

Override file (optional):
    pileup_overrides.txt      - manual pile-up radius overrides
    Format: <stem> <radius_arcsec>
    Example: sw00031659107xpcw3po 6.0
```

The King model `S(r) = S0 * (1 + (r/rc)²)^(-β) + bkg` is fit to the outer wings only (rmin–rmax), then extrapolated inward. Where the data fall below the model indicates the pile-up boundary. The core radius (rc=5.8") and slope (β=1.55) are fixed to the Swift XRT calibration values; only S0 and background are free.

### `swift_pc_source_viewer.py`

Generate zoomed viridis images of the source from PC-mode event files. Overlays source circles, pile-up radii, and optionally `xrtcentroid` positions.

```
Usage:
    swift_pc_source_viewer.py [options]

Requires:
    _pileup.txt files from swift_xrt_king_profile.py

Key options:
    --dimension   Image size: '100px' or '60arcsec' (default: 100px)
    --radius      Overlay circle radius in arcsec (default: 8)
    --sosta       Also plot xrtcentroid positions from
                  source_extraction_OBSID.txt files
    --pdf         Output PDF filename (default: source_images.pdf)

Output (per OBSID):
    {stem}_source.png         - zoomed source image
    Collated PDF with all images
```

Useful for identifying bad columns through the source, anomalous PSF shapes, nearby contaminating sources, or other issues that should lead to excluding an observation.

### `swift_wt_summary_viewer.py`

Summary table and visual diagnostic viewer for WT-mode pointed observations. Shows sky-coordinate images with extraction regions and 1D DETX cross-strip profiles. For multi-orbit observations, produces per-orbit sky image grids.

```
Usage:
    swift_wt_summary_viewer.py --ra <RA> --dec <DEC> [options]

Required:
    --ra        Source RA in degrees
    --dec       Source Dec in degrees

Key options:
    --srcrad      Source circle radius in pixels (default: 20)
    --bkginner    Background annulus inner radius in pixels (default: 80)
    --bkgouter    Background annulus outer radius in pixels (default: 120)
    --expgt       Minimum exposure in seconds (default: 20)
    --compact     Print summary table only, no plots
    --nmax        Process only first N observations
    --pdf         Output PDF filename (default: wt_profiles.pdf)

Output (per OBSID):
    {stem}_wt_combined.png    - sky image + DETX profile
    {stem}_wt_profile.txt     - source position, extraction parameters,
                                BACKSCAL values
    Collated PDF (combined pages + per-orbit grids for multi-orbit obs)
```

The background annulus should be symmetric about 100 pixels (the WT window half-width). The default 80–120 pixel annulus gives 40 pixels of 1D background regardless of source position in the window. See the [UK SSDC BACKSCAL guide](https://www.swift.ac.uk/analysis/xrt/backscal.php) for details.

### `make_wt_master_table.py`

Generate `wt_master_table.txt` with include/exclude flags for WT pointed observations.

```
Usage:
    make_wt_master_table.py [options]

Key options:
    --expmin    Minimum exposure for include=yes (default: 20 seconds)
    --output    Output filename (default: wt_master_table.txt)

Output format:
    OBSID  filename  include  exp(s)  ct/s  n_gti  "comment"
```

Observations below the exposure threshold are set to `include=no` by default. Edit the file to exclude additional observations based on your visual inspection.

### `swift_xrt_extract_spectra.py`

Automated spectral extraction for both PC and WT modes. Calls HEASoft FTOOLS (`xselect`, `xrtexpomap`, `xrtmkarf`, `grppha`) to produce grouped spectra ready for fitting.

```
Usage:
    swift_xrt_extract_spectra.py --ra <RA> --dec <DEC> [options]

Required:
    --ra        Source RA in degrees
    --dec       Source Dec in degrees

Mode selection:
    --mode      pc, wt, or both (default: both)

PC options:
    --rout        Outer radius in arcsec or "auto" (default: 47)
    --bkg-inner   Background inner radius in arcsec (default: 100)
    --bkg-outer   Background outer radius in arcsec (default: 160)

WT options:
    --wt-srcrad     Source radius in pixels (default: 20)
    --wt-bkginner   Background inner radius in pixels (default: 80)
    --wt-bkgouter   Background outer radius in pixels (default: 120)

Common options:
    --mincounts   Minimum counts per grouped bin (default: 20)
    --pctable     PC master table (default: pc_master_table.txt)
    --wttable     WT master table (default: wt_master_table.txt)

Required input files:
    pc_master_table.txt       (for PC mode)
    wt_master_table.txt       (for WT mode)
    {stem}_pileup.txt         (for PC: from king_profile script)
    {stem}_wt_profile.txt     (for WT: from wt_summary_viewer, optional)

Output (per observation):
    {stem}_src.reg            DS9 source region
    {stem}_bkg.reg            DS9 background region
    {stem}_src.pha            Source spectrum
    {stem}_bkg.pha            Background spectrum
    {stem}.arf                Ancillary response file
    {stem}_grp.pha            Grouped spectrum (ready for fitting)
    {stem}_extraction.log     Extraction details
```

**PC mode** uses pile-up radii from `_pileup.txt` to set annular source regions when needed. The outer radius can be fixed or auto-optimized based on signal-to-noise.

**WT mode** uses circular source regions and annular backgrounds. After extraction, BACKSCAL keywords are corrected for WT 1D geometry (source BACKSCAL = 2×r, background BACKSCAL = r_outer − r_inner − 1), following the [UK SSDC standard](https://www.swift.ac.uk/analysis/xrt/backscal.php).

The summary table includes observation dates and RMF filenames, useful for verifying that the CALDB selects the correct response (e.g., `s0` for pre-Sept 2007, `s6` for post-Sept 2007 substrate voltage change).

### `swift_xrt_fit_spectra.py`

Batch spectral fitting using Sherpa. Fits each grouped spectrum independently with a powerlaw model and produces a combined results table and light curve plot.

**Note:** If running within a CIAO environment, `$CALDB` points to the Chandra CALDB, not the HEASoft CALDB. Use `--caldb` to specify the HEASoft CALDB path so that Swift RMFs can be found.

```
Usage:
    swift_xrt_fit_spectra.py --nh <nH> [options]

Required:
    --nh          Galactic nH in units of 10^22 cm^-2

Model selection:
    --model       absorbed (tbabs*ztbabs*powerlaw, requires --redshift)
                  simple (tbabs*powerlaw, default)
    --redshift    Source redshift (required for absorbed model)

Fitting options:
    --defgamma    Frozen gamma for low-count spectra (default: 2.0)
    --mincounts   Minimum counts to fit at all (default: 40)
    --mingamma    Minimum counts for free gamma (default: 200)
    --emin        Lower energy bound in keV (default: 0.3)
    --emax        Upper energy bound in keV (default: 10.0)
    --bkg         subtract or none (default: subtract)

Mode and table selection:
    --modes       pc, wt, or both (default: both)
    --pctable     PC master table (default: pc_master_table.txt)
    --wttable     WT master table (default: wt_master_table.txt)

Other:
    --caldb       Path to HEASoft CALDB (if $CALDB is CIAO's)
    --nmax        Process only first N observations (for testing)
    --output      Results table filename (default: fit_results.txt)
    --plot        Light curve PDF filename (default: flux_lightcurve.pdf)

Output:
    fit_results.txt           Summary table (OBSID, mode, fluxes, gamma, ...)
    flux_lightcurve.pdf       νFν(1 keV) and Γ vs. time (decimal years)
    {stem}_sherpa_fit.log     Per-observation fit log (in OBSID directory)

Fitting logic:
    <40 counts:    skipped entirely
    40–200 counts: gamma frozen to --defgamma, nH and norm free
    >200 counts:   all parameters free

Light curve:
    Upper panel: νFν at 1 keV (erg/cm²/s) on log scale
    Lower panel: photon index Γ
    PC mode: navy circles | WT mode: orange diamonds
    Frozen gamma shown as distinct symbols
```

### `plot_lightcurve.py`

Standalone script to remake the light curve plot from an existing `fit_results.txt` without re-running the fits. Provides options for axis limits, vertical reference lines, and custom titles.

```
Usage:
    plot_lightcurve.py [options]

Key options:
    --input         Input results file (default: fit_results.txt)
    --output        Output PDF (default: flux_lightcurve.pdf)
    --xlim          X-axis limits in decimal years (e.g., --xlim 2005 2024)
    --ylim_flux     Upper panel y-axis limits (e.g., --ylim_flux 1e-12 1e-10)
    --ylim_gamma    Lower panel y-axis limits (e.g., --ylim_gamma 1.0 2.5)
    --times         Vertical reference lines at decimal years
                    (e.g., --times 2007.75 2015.5 2023.3)
    --tlabels       Labels for reference lines (must match --times count)
    --tcolor        Reference line color (default: green)
    --title         Custom plot title
    --figsize       Figure dimensions in inches (default: 14 7)
    --dpi           Output resolution (default: 150)

Examples:
    # Default plot
    plot_lightcurve.py

    # Customized
    plot_lightcurve.py --xlim 2005 2024 \
        --ylim_flux 5e-12 5e-11 --ylim_gamma 1.2 2.2 \
        --times 2007.75 2015.5 --tlabels "Voltage change" "Flare" \
        --title "3C 273 X-ray Light Curve" --output lc_custom.pdf
```

### `parallel_extract.py`

Run spectral extraction in parallel by splitting master tables into chunks, each processed by a separate worker. Each worker runs in its own temporary directory with symlinks to the OBSID data, avoiding xselect session file conflicts.

```
Usage:
    parallel_extract.py --ra <RA> --dec <DEC> --nproc <N> [extraction options]

Key options:
    --nproc       Number of parallel workers (default: 8)
    --mode        pc, wt, or both (default: both)
    --dryrun      Show chunk splitting without running

All other arguments (--rout, --mincounts, --bkg-inner, etc.) are
passed through to swift_xrt_extract_spectra.py.

Example:
    parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 16
    parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 32 --mode pc --dryrun
```

### `parallel_fit.py`

Run spectral fitting in parallel, then merge results into a single `fit_results.txt` and generate the combined light curve plot.

```
Usage:
    parallel_fit.py --nh <nH> --nproc <N> [fitting options]

Key options:
    --nproc       Number of parallel workers (default: 8)
    --dryrun      Show splitting without running

All other arguments (--model, --redshift, --caldb, etc.) are
passed through to swift_xrt_fit_spectra.py.

Example:
    parallel_fit.py --nh 0.0179 --model simple --caldb /opt/CALDB --nproc 16
    parallel_fit.py --nh 0.0179 --redshift 0.158 --caldb /opt/CALDB --nproc 32
```

After all chunks complete, results are merged and sorted by OBSID, and `plot_lightcurve.py` is called automatically to generate the combined plot.

---

## Notes

### Parallelization

For datasets with 100+ observations, the extraction and fitting steps can take hours when run sequentially. The `parallel_extract.py` and `parallel_fit.py` wrappers split the master tables into chunks and run workers in parallel. Each worker operates in an isolated temporary directory with symlinks to the OBSID data directories, avoiding file conflicts (particularly the xselect session files that would collide if multiple instances ran in the same directory).

A reasonable starting point is `--nproc` equal to half the number of CPU cores, since each xselect/xrtmkarf process is itself somewhat I/O bound. For a 48-core machine, `--nproc 16` to `--nproc 24` is a good range. Use `--dryrun` first to verify the chunk splitting.

### CALDB and CIAO coexistence

CIAO sets `$CALDB` to its own Chandra calibration directory. If you run the fitting script from a CIAO environment, use `--caldb /path/to/heasoft/caldb` to point to the HEASoft CALDB containing Swift XRT response files.

The extraction script requires HEASoft tools and should be run from a HEASoft environment (not CIAO).

### WT mode BACKSCAL

In WT mode, the CCD reads out as a 1D strip. XSELECT sets BACKSCAL based on the 2D area of the extraction region, which is incorrect. The extraction script automatically corrects this to reflect the 1D extent in the DETX direction. See the [UK SSDC documentation](https://www.swift.ac.uk/analysis/xrt/backscal.php) for details.

### File naming conventions

Swift XRT event files follow the pattern `sw[OBSID]x{pc,wt}w{N}{sl,st,po}_cl.evt`:
- `xpc` / `xwt`: Photon Counting / Windowed Timing mode
- `w1`–`w4`: CCD window size (w1 fastest readout, w4 slowest)
- `sl`: slew, `st`: settling, `po`: pointed (only `po` files are used for science)
- `_cl`: cleaned level-2 data

### Pile-up thresholds

- **PC mode:** ~0.5 ct/s (depends on window mode)
- **WT mode:** ~150 ct/s

The King profile script automatically detects and measures pile-up for PC mode. For WT mode at typical blazar count rates (<100 ct/s), pile-up is not a concern.

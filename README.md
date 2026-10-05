# Swift XRT Spectral Analysis Pipeline

A Python-based pipeline for automated spectral extraction, fitting, and light curve generation from Swift X-Ray Telescope (XRT) data in both Photon Counting (PC) and Windowed Timing (WT) modes.

Developed for analysis of point sources (e.g., blazars, AGN) observed across multiple epochs. The pipeline handles pile-up assessment, optimal extraction region determination, proper BACKSCAL correction for WT mode, and batch spectral fitting with Sherpa.

**Author:** Eileen T. Meyer ([@etmeyer](https://github.com/etmeyer))  
**License:** MIT

> **Documentation:** see [docs/](docs/index.md) for the full step-by-step walkthrough.

---

## Quick start

Download and reduce one epoch of 3C273 data, end to end. The pipeline uses
**two terminals** — HEASoft for Steps 2–7, CIAO for the Step 8 fit — because
CIAO breaks the HEASoft tools and Sherpa only exists in CIAO. See
[docs/01-setup.md](docs/01-setup.md#two-terminals).

```bash
# Terminal 1 (HEASoft) -- never run `ciao` here
setup_swiftxrt; heainit
swift_xrt_download.py --name "3C 273" \
    --start-date 2008-08-04 --end-date 2011-07-06 \
    --outdir XRT_input                       # Step 2
# Steps 3-7: see docs/index.md

# Terminal 2 (CIAO) -- Step 8 only
setup_swiftxrt; ciao
```

---

## Setup

See [docs/01-setup.md](docs/01-setup.md) for installation, environment
setup, the two-terminal rule, and the `swift_xrt_doctor.py` verification
script.

---

## Workflow Overview

The pipeline is designed to be run step-by-step, with visual inspection at
each stage. Step numbers match [docs/](docs/index.md).

```
Step  What                         Script                        Terminal
1     Setup                        swift_xrt_doctor.py           both
2     Download                     swift_xrt_download.py         either
3     Run xrtpipeline              xrt_pipeline.py               HEASoft
4     Survey observations          swift_xrt_summary.py          either
5     PC-mode inspection
        a. Pile-up / PSF analysis  swift_xrt_king_profile.py     either
        b. Source images           swift_pc_source_viewer.py     either
        c. PC master table         (shell command)               either
6     WT-mode inspection
        a. Profile viewer          swift_wt_summary_viewer.py    either
        b. WT master table         make_wt_master_table.py       either
        c. Edit master tables      (manual review)
7     Extract spectra              swift_xrt_extract_spectra.py  HEASoft
8     Fit spectra & plot           swift_xrt_fit_spectra.py      CIAO (--caldb)
        (optional) Customize plot  plot_lightcurve.py            either
```

### Step 2 — Download → see [docs/02-download.md](docs/02-download.md)

### Step 3 — Run xrtpipeline → see [docs/03-xrtpipeline.md](docs/03-xrtpipeline.md)

In the HEASoft terminal:

```bash
xrt_pipeline.py --batch --nproc 16 --indir XRT_input --outdir XRT_output \
    --ra 187.2779 --dec 2.0524
```

All subsequent steps are run from within the output directory:

```bash
cd XRT_output
```

### Step 4 — Survey all observations → see [docs/04-survey.md](docs/04-survey.md)

```bash
swift_xrt_summary.py --compact --ra 187.2779 --dec 2.0524   # one row per OBSID
swift_xrt_summary.py --ra 187.2779 --dec 2.0524             # detailed per-OBSID tables
```

### Step 5 — PC-mode inspection → see [docs/05-pc-inspection.md](docs/05-pc-inspection.md)

```bash
swift_xrt_king_profile.py --ra 187.2779 --dec 2.0524   # 5a: pile-up radii
swift_pc_source_viewer.py                              # 5b: source images

# 5c: the PC master table
find . -name '*xpc*po*_cl.evt' | sort | \
  awk -F'/' '{obsid=$2; file=$NF; gsub(/^\.\//, "", obsid); \
  sub(/_cl\.evt$/, "", file); \
  printf "%-14s %-35s %-10s %-10s \"%s\"\n", obsid, file, "yes", "no", ""}' | \
  (printf "%-14s %-35s %-10s %-10s %s\n" "OBSID" "filename" "include" "badstripe" "comment"; cat) \
  > pc_master_table.txt
```

To change a pile-up radius, add the file stem and radius to `pileup_overrides.txt` and re-run `swift_xrt_king_profile.py`; extraction reads the radius from `_pileup.txt`.

### Step 6 — WT-mode inspection → see [docs/06-wt-inspection.md](docs/06-wt-inspection.md)

```bash
swift_wt_summary_viewer.py --ra 187.2779 --dec 2.0524   # 6a: WT plots, regions, detection
make_wt_master_table.py                                 # 6b: the WT master table
```

**6c.** Edit `pc_master_table.txt` and `wt_master_table.txt`: set `include` to `no` for any observation you want to exclude, with the reason in the comment. Re-running `make_wt_master_table.py` overwrites your edits; use `--output` to write a new file instead.

### Step 7 — Extract spectra → see [docs/07-extract.md](docs/07-extract.md)

In the HEASoft terminal (`setup_swiftxrt; heainit`, **no** `ciao`); the
scripts refuse to start in a CIAO terminal.

```bash
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524     # one at a time
parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 16   # or in parallel
```

Both exit non-zero and list the observations that failed, if any.

### Step 8 — Fit spectra and plot the light curve → see [docs/08-fit-and-plot.md](docs/08-fit-and-plot.md)

In a separate CIAO terminal (`setup_swiftxrt; ciao`) — Sherpa only exists
there. CIAO points `$CALDB` at the Chandra CALDB, so always pass the HEASoft
CALDB with `--caldb`:

```bash
swift_xrt_fit_spectra.py --nh 0.0179 --caldb /opt/CALDB             # one at a time
parallel_fit.py --nh 0.0179 --caldb /opt/CALDB --nproc 16           # or in parallel

# With absorption in the source, at its redshift:
parallel_fit.py --nh 0.0179 --redshift 0.158 --caldb /opt/CALDB --nproc 16
```

Output: `fit_results.txt` (table) and `flux_lightcurve.pdf` (νFν and Γ vs. time, with PC and WT points color-coded).

**Optional:** customize the light curve plot:

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

Run the HEASoft `xrtpipeline` task on raw Swift XRT data to produce cleaned level-2 event files and exposure maps. Supports single-OBSID, sequential batch, and parallel batch modes. Details: [docs/03-xrtpipeline.md](docs/03-xrtpipeline.md).

```
Usage:
    # Single OBSID
    xrt_pipeline.py --indir XRT_input/00035017001 \
        --outdir XRT_output/00035017001 --ra 187.2779 --dec 2.0524

    # Batch: all OBSIDs under input directory
    xrt_pipeline.py --batch --indir XRT_input --outdir XRT_output \
        --ra 187.2779 --dec 2.0524

    # Parallel batch
    xrt_pipeline.py --batch --nproc 16 --indir XRT_input \
        --outdir XRT_output --ra 187.2779 --dec 2.0524

Key options:
    --batch            Process all OBSID subdirs under --indir
    --nproc            Number of parallel workers (default: 1)
    --timeout          Per-OBSID time limit in seconds (default: 600)
    --createexpomap    Create exposure maps (default: yes)
    --extractproducts  Also run xrtproducts (default: no; unused downstream)
    --cleanup          Remove intermediate files (default: no)
    --clobber          Overwrite existing output (default: yes)

Output (per OBSID, in <outdir>/<OBSID>/):
    Cleaned event files (*_cl.evt) for all modes (PC, WT) and
    segments (pointed, settling, slew), exposure maps (*_ex.img),
    intermediate files, and xrtpipeline_<OBSID>.log.

Exit status: 0 only if every OBSID succeeded.

Requires: HEASoft (xrtpipeline), CALDB; run in a HEASoft terminal
```

Each `xrtpipeline` call processes everything within an OBSID — all modes and observation types — producing separate cleaned event files for each. Each run gets a private parameter-file folder and its own pseudo-terminal, so parallel batches don't collide and also run under `nohup` or `cron`. An OBSID counts as OK only if `xrtpipeline` exits 0 and leaves a cleaned event file and exposure map for each mode with pointed data.

### `swift_xrt_summary.py`

Crawl OBSID directories and report on all cleaned event files. Shows each OBSID's target and pointing offset, mode sequences (WT settling → WT pointed → PC pointed), exposures, count rates, GTI/orbit structure, and pile-up warnings. Details: [docs/04-survey.md](docs/04-survey.md).

```
Usage:
    swift_xrt_summary.py --ra 187.2779 --dec 2.0524              # detailed tables
    swift_xrt_summary.py --compact --ra 187.2779 --dec 2.0524    # one row per OBSID

Options:
    --compact     One row per OBSID
    --ra/--dec    Source position; pointing offsets are measured from it
                  (default: from each observation's own target)

Output (terminal only):
    Per-OBSID tables with target, pointing offset, mode, exposure,
    events, count rate
    GTI orbit analysis (segment durations, gap lengths, duty cycle)
    Pile-up warnings (PC >0.5 ct/s, WT >150 ct/s)
    Compact table with columns:
        OBSID, Date/Time, Total(ks), ct/s, Slew_i, Slew_f,
        N_WT, N_PC, WT_exp(ks), PC_exp(ks), Orb, Seq, Off('), Target
    OBSIDs with no cleaned event files are listed under the compact table

Sequence codes:
    1=WT_SLEW  2=PC_SLEW  3=WT_SETTLING  4=PC_SETTLING
    5=WT_POINTED  6=PC_POINTED
```

### `swift_xrt_king_profile.py`

Fit King profiles to PC-mode radial surface brightness profiles to assess pile-up and determine source extraction regions. Details: [docs/05-pc-inspection.md](docs/05-pc-inspection.md).

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
    --rc        King core radius, fixed (default: 5.8")
    --beta      King beta slope, fixed (default: 1.55)
    --sbthresh  Pile-up threshold, counts/frame/arcmin^2 (default: 4.5)
    --sigma     Diagnostic ring flags, 2 consecutive bins (default: 3.0)
    --sigma2    Diagnostic ring flags, single bin (default: 4.0)
    --pdf       Output PDF filename (default: king_profiles.pdf)

Output (per OBSID):
    {stem}_king_profile.png   - diagnostic plot
    {stem}_pileup.txt         - centroid, plate scale, pile-up radius

Override file (optional):
    pileup_overrides.txt      - manual pile-up radius overrides
    Format: <stem> <radius_arcsec>
    Example: sw00031659107xpcw3po 6.0
    Re-run this script after editing it: extraction reads the
    radius from _pileup.txt, not from the override file.
```

The radial profile divides each ring by its exposed area, from the xrtpipeline exposure map (`{stem}_ex.img`), so bad columns crossing the source don't distort it; without the map it falls back to the geometric area, with a warning. The King model `S(r) = S0 * (1 + (r/rc)²)^(-β) + bkg` is fit to the outer wings only (rmin–rmax), where pile-up doesn't reach. The core radius (rc=5.8") and slope (β=1.55) are fixed; only S0 and background are free. The pile-up radius is where this profile, with S0 at the top of its 1σ range, falls below `--sbthresh` counts per frame per arcmin² (pile-up is per frame, so this scales with the window's frame time). If the radius reaches into the fitting annulus, the wings are refit further out.

The threshold was calibrated on spectra. For five 3C 273 PC observations, spectra extracted with inner radii of 0–28″ show the 1 keV flux levelling off once the piled-up core is excluded, at 10–16″. At 4.5, with the exposure-corrected profile, the radius lands 1.5–5.4″ beyond that in all five, and the fluxes agree with the levelled-off values. Rings that differ from the model by `--sigma`/`--sigma2` are still flagged in the plot and `_pileup.txt`, but only as a diagnostic: they mark where the core deficit becomes significant (10–12″ for these), a few arcsec short of where the flux levels off. On the uncorrected profile, bad columns made the flags set 20–24″ where 10–14″ suffices, and 10″ for an off-axis source whose flux was then 17% low. Details and the test: [docs/05-pc-inspection.md](docs/05-pc-inspection.md#how-the-radius-is-chosen).

### `swift_pc_source_viewer.py`

Generate zoomed viridis images of the source from PC-mode event files. Overlays source circles, pile-up radii, and optionally `xrtcentroid` positions. Details: [docs/05-pc-inspection.md](docs/05-pc-inspection.md).

```
Usage:
    swift_pc_source_viewer.py [options]

Requires:
    _pileup.txt files from swift_xrt_king_profile.py

Key options:
    --dimension   Image size: '100px' or '60arcsec' (default: 100px)
    --radius      Overlay circle radius in arcsec (default: 8)
    --sosta       Also plot xrtcentroid positions from
                  source_extraction_OBSID.txt files (not made by
                  this pipeline)
    --pdf         Output PDF filename (default: source_images.pdf)

Output (per OBSID):
    {stem}_source.png         - zoomed source image
    Collated PDF with all images
```

Useful for identifying bad columns through the source, anomalous PSF shapes, nearby contaminating sources, or other issues that should lead to excluding an observation.

### `swift_wt_summary_viewer.py`

Summary table and visual diagnostic viewer for WT-mode pointed observations. Shows sky-coordinate images with extraction regions and 1D DETX cross-strip profiles. For multi-orbit observations, produces per-orbit sky image grids. Details: [docs/06-wt-inspection.md](docs/06-wt-inspection.md).

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
    --detsigma    Minimum source significance above background to
                  count as detected (default: 3)
    --compact     Print summary table only, no plots
    --nmax        Process only first N observations
    --pdf         Output PDF filename (default: wt_profiles.pdf)

Output (per OBSID):
    {stem}_wt_combined.png    - sky image + DETX profile
    {stem}_wt_profile.txt     - source position, extraction parameters,
                                BACKSCAL values, source detection
                                (counts, sigma, source_detected)
    Collated PDF (combined pages + per-orbit grids for multi-orbit obs)
```

The background annulus should be symmetric about 100 pixels (the WT window half-width). The default 80–120 pixel annulus gives 40 pixels of 1D background regardless of source position in the window. See the [UK SSDC BACKSCAL guide](https://www.swift.ac.uk/analysis/xrt/backscal.php) for details.

### `make_wt_master_table.py`

Generate `wt_master_table.txt` with include/exclude flags for WT pointed observations. Details: [docs/06-wt-inspection.md](docs/06-wt-inspection.md).

```
Usage:
    make_wt_master_table.py [options]

Key options:
    --expmin    Minimum exposure for include=yes (default: 20 seconds)
    --output    Output filename (default: wt_master_table.txt)

Output format:
    OBSID  filename  include  exp(s)  ct/s  n_gti  "comment"
```

Observations below the exposure threshold, or whose `_wt_profile.txt` says the source was not detected, are set to `include=no` by default, with the reason in the comment. Run `swift_wt_summary_viewer.py` first so the detection results exist. Edit the file to exclude additional observations based on your visual inspection; the script rewrites the file from scratch, so re-run it with `--output` to avoid losing those edits.

### `swift_xrt_extract_spectra.py`

Automated spectral extraction for both PC and WT modes. Calls HEASoft FTOOLS (`xselect`, `xrtmkarf`, `grppha`), with the exposure maps xrtpipeline made in Step 3, to produce grouped spectra ready for fitting. See [docs/07-extract.md](docs/07-extract.md).

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

WT options (used only where there is no _wt_profile.txt):
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
    {stem}_ex.img             (exposure map, from xrtpipeline in Step 3)

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

**WT mode** uses the source position and radii from `_wt_profile.txt` (circular source region, annular background). After extraction, BACKSCAL keywords are corrected for WT 1D geometry (source BACKSCAL = 2×r, background BACKSCAL = r_outer − r_inner − 1), following the [UK SSDC standard](https://www.swift.ac.uk/analysis/xrt/backscal.php).

The summary table includes observation dates and RMF filenames, useful for verifying that the CALDB selects the correct response (e.g., `s0` for pre-Sept 2007, `s6` for post-Sept 2007 substrate voltage change).

### `swift_xrt_fit_spectra.py`

Batch spectral fitting using Sherpa. Fits each grouped spectrum independently with an absorbed power law (W-stat, `wilm` abundances, 1σ errors) and produces a combined results table and light curve plot. See [docs/08-fit-and-plot.md](docs/08-fit-and-plot.md).

**Note:** Run this from the CIAO terminal (Sherpa lives only in CIAO's Python). There, `$CALDB` points to the Chandra CALDB, not the HEASoft CALDB, so pass `--caldb` (e.g. `--caldb /opt/CALDB`) so that Swift RMFs can be found; the script stops with an explanation if neither has them.

```
Usage:
    swift_xrt_fit_spectra.py --nh <nH> [options]

Required:
    --nh          Galactic nH in units of 10^22 cm^-2

Model selection:
    --model       absorbed (tbabs*ztbabs*powerlaw, requires --redshift)
                  simple (tbabs*powerlaw)
                  (default: absorbed with --redshift, else simple)
    --redshift    Source redshift (required for absorbed model)
    --abund       Abundance table for tbabs (default: wilm)

Fitting options:
    --stat        wstat (default) or chi2 (the method before Oct 2026)
    --defgamma    Gamma for spectra under --mingamma counts
                  (default: median of the free fits in this run)
    --mincounts   Minimum counts to fit at all (default: 40)
    --mingamma    Minimum counts for free gamma (default: 200)
    --emin        Lower energy bound in keV (default: 0.3)
    --emax        Upper energy bound in keV (default: 10.0)
    --bkg         use or none (default: use)

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
    ≥200 counts:   fitted first, all parameters free
    40–199 counts: fitted next, gamma frozen at the median of the
                   free fits (or --defgamma)
    All errors are 1 sigma.

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

Full per-chunk logs: parallel_extract_logs/chunkNN_<mode>.log
(the console shows only a short tail for failed chunks).
At the end it checks that every included observation got a new
_grp.pha, and exits non-zero if not or if any chunk failed.

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
    parallel_fit.py --nh 0.0179 --caldb /opt/CALDB --nproc 16
    parallel_fit.py --nh 0.0179 --redshift 0.158 --caldb /opt/CALDB --nproc 32
```

It fits in two phases, like a single run: the spectra with at least `--mingamma` counts first, then the rest with Γ frozen at the median of all those fits. After all chunks complete, results are merged and sorted by OBSID, and `plot_lightcurve.py` is called automatically to generate the combined plot.

---

## Notes

### Parallelization

For datasets with 100+ observations, the extraction and fitting steps can take hours when run sequentially. The `parallel_extract.py` and `parallel_fit.py` wrappers split the master tables into chunks and run workers in parallel. Each worker operates in an isolated temporary directory with symlinks to the OBSID data directories, avoiding file conflicts (particularly the xselect session files that would collide if multiple instances ran in the same directory).

A reasonable starting point is `--nproc` equal to half the number of CPU cores, since each xselect/xrtmkarf process is itself somewhat I/O bound. For a 48-core machine, `--nproc 16` to `--nproc 24` is a good range. Use `--dryrun` first to verify the chunk splitting.

### HEASoft and CIAO: two terminals

HEASoft and CIAO cannot share a terminal for this pipeline. Once CIAO is set
up, `python3` is CIAO's wrapper script, which resets `$HEADAS` and `$CALDB` to
CIAO's own trees inside every pipeline script, so `xrtpipeline`, `xselect`,
`xrtmkarf` and `grppha` fail. The fit, in turn, needs Sherpa, which only CIAO's
Python has. Run Steps 2–7 in a HEASoft terminal and the Step 8 fit in a CIAO
terminal with `--caldb /path/to/heasoft/caldb`. The scripts check this and say
which terminal to use; `swift_xrt_doctor.py` lists which steps the current
terminal can run. Details: [docs/01-setup.md](docs/01-setup.md#two-terminals).

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

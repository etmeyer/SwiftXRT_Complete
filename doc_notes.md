# Doc-session notes (out-of-scope items for future `.doc` sessions)

Logged during **1.2.1.doc** (Step 1 README rewrite). None of these were fixed in
that PR — they belong to later doc sessions or to the fix side.

- **Sample doctor blocks reflect the post-Step-3 state.** The green and yellow
  output blocks in the rewritten Installation section show `scipy` as `[OK]`
  (the state *after* `conda install -n ciao-4.16 scipy`). The literal captures in
  `test_runs/3c273_epoch1/doctor_green.txt` and `doctor_yellow.txt` still show
  `[FAIL] scipy NOT IMPORTABLE` (and fail counts of 1), because scipy is **not
  yet installed** into the `ciao-4.16` env on amorgos. The scipy version in the
  README (`1.11.4`) is illustrative — regenerate the captures with
  `run_doctor_scenarios.sh` after the scipy install and reconcile the exact
  version string if it differs. (So `doctor_green.txt` is, today, not actually
  all-green despite its name.)

- **Workflow Step 5** ("Initialize HEASoft and CALDB first, then:") could
  forward-reference the new Installation Steps 2–3 once those step descriptions
  are revised (1.2.7.doc).

- **Requirements** historically listed `xrtexpomap` as a HEASoft tool the
  pipeline uses; the extraction Script Reference confirms it (`xselect`,
  `xrtexpomap`, `xrtmkarf`, `grppha`). `swift_xrt_doctor.py` does not probe
  `xrtexpomap` specifically — not wrong, just noting the checklist isn't
  exhaustive of every FTOOL.

Logged during **1.2.2.fix** (Step 2/3 `swift_xrt_download.py`). Out of scope for
that branch (`fix/step2-download-script`); flagged for the relevant owners.

- ~~**`xrt_pipeline.py` wrapper fails where bare `xrtpipeline` succeeds (env).**~~
  FIXED in PR #5; the by-hand workaround is gotcha 8 in
  `docs/03-xrtpipeline.md`. Running the wrapper headless on amorgos aborts in the `prefilter`/`xrtfilter`
  step with `couldn't get parameter 'leapname' [file not found (or has wrong
  access type)]` (PIL_BAD_FILE_ACCESS) — even though `pget prefilter leapname`
  resolves to `$HEADAS/refdata/leapsec.fits` and that file is readable. Running
  `xrtpipeline` **directly** with the same input tree, CALDB, and a fresh
  `PFILES=/tmp/...;$HEADAS/syspfiles` completes cleanly (exit 0, produces
  `*_cl.evt` + `*_ex.img` for PC and WT). So the leapname failure is in how the
  wrapper sets up the environment for the spawned FTOOLS, not in the data or the
  download. Belongs to the Step-4 (`xrt_pipeline.py`) session. NB: this is
  separate from the known `[MISSING] attitude_file/hk_file`-on-success post-check
  bug already in the 1.1 bug log (Step 4 entry).

- ~~**`xrt_pipeline.py` suppresses xrtpipeline's own stdout/stderr.**~~ FIXED in
  PR #5; 1.2.3.doc made the printed lines skip xrtpipeline's closing
  banner so they show the error itself. On failure
  the wrapper prints only `[FAILED] … exit code N` + `[MISSING] …` lines; the
  real xrtpipeline error is buried in `<outdir>/xrtpipeline_<obsid>.log`. Made
  diagnosing the leapname issue slower than it needed to be. Step-4 session.

- **`swift_xrt_download.py resolve_data_url()` date-based HEASARC strategy is
  dead with the live catalog.** `start_time` from `swiftmastr` is an **MJD
  float** (e.g. `56725.7111`), but `resolve_data_url` builds the HEASARC
  date path via `val[:4] + "_" + val[5:7]` (assumes `YYYY-MM-DD`), yielding a
  bogus `5672_5.` month. Harmless today because the UKSSDC mirror strategy
  succeeds first, but the `heasarc_date` branch never produces a valid URL. Not
  touched in 1.2.2.fix (URL discovery was out of scope); worth a real fix or
  removal in a later download-script session.

Logged during **1.2.3.fix** (Step 3 `xrt_pipeline.py` wrapper hardening). The
two `xrt_pipeline.py` items above (headless leapname/PIL abort, swallowed
stdout/stderr) were the targets of this session and are now FIXED on
`fix/step3-xrtpipeline-wrapper` — they can be struck from this file by the
paired `1.2.3.doc` session. Diagnosis for the record: the "leapname" error was
**not** matplotlib/Python-Imaging; it was HEASoft's Parameter Interface Library
reading a CIAO-polluted `$HOME/pfiles/prefilter.par` (CIAO and HEASoft share
`$HOME/pfiles`, and CIAO's prefilter had learned a non-existent
`/opt/ciao/.../leapsec.fits` path). Fixed with a fresh per-run private PFILES.
The headless `/dev/tty` abort is a separate cause (controlling-terminal
requirement), fixed with a pty — `mode=h` does **not** help and is rejected by
xrtpipeline's parser. New out-of-scope finding:

- **`xrt_pipeline.py` runs xrtpipeline with `cwd=<input OBSID dir>`**, so the
  cleaning step's xselect drops session files (`xselect.log`,
  `xsel_timefile.asc`, `xselect<pid>*.flt`) into the *input* tree rather than a
  scratch dir. Harmless today (per-OBSID dirs are unique, downstream is
  unaffected, and the data tree is gitignored), but it mutates the read-only-ish
  archive on every run and isn't reproducible-clean. Left as-is in 1.2.3.fix
  (out of the four-issue scope; changing cwd needs its own validation). A future
  hardening could run each OBSID from a private scratch cwd. Still open;
  documented as gotcha 7 in `docs/03-xrtpipeline.md` (1.2.3.doc).

Logged during **two-terminal env fix** (`fix/two-shell-env`, prompted by a
user report of "step 5" errors blamed on the environment not finding Sherpa).
Diagnosis: no single terminal could run the documented pipeline. In a CIAO
terminal, CIAO's `python3` wrapper resets `$HEADAS`/`$CALDB` inside every
script, so extraction (and xrtpipeline) fail; without CIAO there is no Sherpa.
`docs/01-setup.md`'s "a single CIAO shell can run the entire pipeline" was
wrong. Fixed on that branch: docs rewritten around two terminals; scripts
refuse to start in the wrong terminal with an explanation; doctor lists what
the current terminal can run; README renumbered to the docs/ step numbers (the
"Workflow Step 5" forward-reference item above is resolved by this); Bug A
(dead ARF/BACKFILE paths → ~100× fluxes) fixed; extraction/fit exit non-zero
on failure. Still open:

- ~~Bug B race~~ — FIXED on `fix/bugb-parallel-extract`. Cause: all
  parallel workers shared `$HOME/pfiles`; an xselect occasionally read
  `extractor.par` while another rewrote it (`Can't stat user parameter file
  .../extractor.par`, `Error in extractor`) and that OBSID got no spectrum.
  Each extraction run now uses a private PFILES dir (as xrt_pipeline.py
  does). Reproduced 4/7 runs dropping an OBSID at `--nproc 16` on the
  88-spectrum epoch-1 set; 0/8 after the fix. Full chunk logs now kept in
  `parallel_extract_logs/`.
- ~~`requests` missing from the `heasoft` conda env~~ — installed
  2026-10-02 (`requests` + `charset-normalizer` only, with openssl/certifi/
  ca-certificates pinned to their existing versions); the "download before
  `heainit`" workaround was removed from the docs.
- ~~`swift_xrt_download.py --obsid X --list-only` downloads anyway~~ —
  FIXED on `fix/download-list-only`: the `--list-only` exit lived only in the
  catalog (`--name`/`--ra`) branch; `--obsid`/`--obsid-file` now list the
  OBSIDs and exit. The docs flowchart had the same bypass and is corrected.
- **xselect silently loses output for long paths.** At ~140 characters it
  printed "Wrote spectrum to ..." but no file appeared (and it crashed in
  `xsl_exit`); ~100 characters work. Documented as a gotcha only.
- ~~OBSID 041-style 0-count WT spectra~~ — FIXED on
  `fix/wt-source-detection`. 041 was not a centroiding bug: it was pointed
  5.1' off target with 3C 273 outside the WT window (0 events within 20 px
  of the target; flat DETX profile), so the data genuinely contain no
  source. `swift_wt_summary_viewer.py` now tests the source against the
  1D-scaled background (041: -8.4 sigma; faintest real epoch-1 obs: 39.5
  sigma), records it in `_wt_profile.txt`, and `make_wt_master_table.py`
  sets include=no with the reason in the comment column.
  Side note: `find_sky_position` derives the target X from RA only and Y
  from the median of nearby events (never from Dec). It works for on-axis
  data but is crude; a proper TAN-WCS seed would be more robust.
- ~~King-profile low-count guard~~ — FIXED on `fix/king-lowcount-guard`.
  043 (122 s) and 080 (75 s) got 2"/4" because with ~200-300 events the
  residuals can't reach 3 sigma. The four long PC exposures of the same
  ~2.5 ct/s source put the pile-up edge at ~1 ct/s/arcmin^2 of the wing-fit
  King profile (w3; 1.75 counts/frame/arcmin^2); applying that to 043/080
  gives 17.8"/17.2". When S0 error > 10% the script now uses that radius.
  Evidence the old radii were wrong: 080's PC flux (4" exclusion) was 2.4x
  below its own WT segment; 045's well-excised PC agrees with its WT to 3%.
  With proper excision 043/080 have 57/39 counts left, too few to fit, so
  they drop out of the light curve (080 keeps its WT point).
  The fitter now counts "< 3 noticed bins" as a skip, not a failure.
- 073's 24" radius (flagged in May as possibly over-estimated) is a
  profile measurement and is untouched; the PSF-threshold radius for it is
  20.4", consistent within the method's spread.

Logged during **1.2.3.doc** (Step 3 page, `docs/step3-xrtpipeline`). Writing
the page turned up four `xrt_pipeline.py` problems, fixed on that branch:
`--batch` always exited 0; the failure lines showed xrtpipeline's closing
banner instead of the error; batch mode never checked products (an exit 0
without a cleaned file would have dropped the OBSID silently), and the
single-OBSID check expected WT products from WT settling-only data (073,
074); `--createexpomap no` always failed because xrtpipeline's `useexpomap`
stayed `yes`. Other findings:

- **037 no longer hangs, even with `--extractproducts yes`** (31 s on
  2026-10-04). Probably PR #5's stdin from `/dev/null`: in May the hung
  xselect was waiting on the terminal. The timeout stays as a backstop.
- **On a timeout the printed lines say little.** xrtpipeline's output is
  block-buffered into the pipe, so a killed run has usually flushed only its
  opening banner. The reason line still says "timeout after N seconds".
- **Extraction's own exposure-map path never runs** in the `XRT_output`
  layout: `find_auxiliary_files` looks for `sw<OBSID>*pat|sat.fits*` and
  `*xhd.hk*` in the output folder, but the attitude file stays in
  `XRT_input/<OBSID>/auxil/` and xrtpipeline writes `xhdtc.hk`. So Step 7
  always uses xrtpipeline's `*_ex.img`, and fails without it. For 1.2.7.doc.
- **The King profile and PC viewer read every `*xpc*_cl.evt`**, not just
  pointed ones (the WT tools and both master tables take `po` only).
  Harmless for 3C 273, whose PC data are all pointed, but a PC settling or
  slew file would be profiled. For 1.2.5.doc.

Logged during **1.2.4.doc** (Step 4 page, `docs/step4-survey`).

- **To do: a documentation website.** Build a GitHub Pages site (sidebar,
  search; e.g. MkDocs) from `docs/` so the pages read like a manual. Asked
  for on 2026-10-04; planned after the remaining doc pages. Until then the
  docs read best on GitHub (the repo is public), or in a terminal with
  `python3 -m rich.markdown -p docs/<page>.md`.
- ~~`docs/02-download.md` said the survey shows each pointing's off-axis
  angle; it didn't~~ — FIXED on `docs/step4-survey`: the survey now shows
  each OBSID's target and pointing offset (from `--ra`/`--dec`, else the
  observation's own target). On epoch 1 this flags OBSID 041 at 5.1′ (all
  others within 3.6′) already at Step 4.
- ~~`--compact` silently drops OBSIDs with no cleaned event files~~ — FIXED
  on the same branch; they are listed under the table.
- **Steps 4–7 also read compressed files.** The survey, King profile, PC and
  WT viewers and WT master table also glob `*_cl.evt.gz`, and extraction
  globs `*_ex.img*`. Run in `XRT_input` instead of `XRT_output`, they would
  quietly use the archive's own cleaned files (the PC master-table one-liner
  finds nothing there). The Step 4 page warns about it; later pages should
  say "run inside `XRT_output`" as well.

Logged during **1.2.5.doc** (Step 5 page, `docs/step5-pc-inspection`).
Findings:

- ~~The pile-up radius counts excesses too~~ — REPLACED on
  `docs/step5-pc-inspection` (the author asked for a smarter automatic
  radius; overrides had been set by hand). The residual flags (|residual|,
  3σ/4σ) set 018 and 073 at 20″/24″ from an excess at 13–23″, and the
  off-axis SDSS J122933 pointing at 10″. Ground truth came from spectra:
  for 018, 045, 073, 074 and 00091742013, inner radii 0–28″ were extracted
  and fit with the pipeline (overrides + Steps 7–8), and the 1 keV flux
  levels off at 10–16″ (2σ test against 20″, using
  σ_diff² = σ_outer² − σ_inner² for nested annuli). The flags missed in
  both directions, and a deficit-only test (with either PSF shape) gave
  8–10″ for 018/073, where their flux was still 9–18% low. The radius is
  now the PSF surface-brightness threshold (S0 + 1σ) at 4.0
  counts/frame/arcmin² for every observation; PR #11's 1.75 had been
  calibrated to the inflated flag radii. Fluxes at the new radii match the
  levelled-off values within 0.5σ, with 21–34% smaller errors for 018, 073
  and 074. The calibration covers one source at 2.3–3.6 ct/s.
  Test data stay on amorgos in `test_runs/pu/` (gitignored): `R00`–`R28`
  extractions and fits, `scan.json`, `NEW/` (new method), `analysis/`.
- **The CALDB PC PSF (`swxpsf20010101v006`, 2020) is not the script's
  King shape.** v006 is a King with rc = 1.581 px (3.7″) and β = 1.305 plus
  a 7.5% Gaussian with σ = 3.149 px; the script's fixed rc = 5.8″,
  β = 1.55 matches v004/v005. Neither describes the 12–20″ rings
  consistently: data/model there runs from +34% (018) to −10% (045) with
  the script's shape and +16% to −20% with v006, varying by observation.
  Likely causes: bad columns crossing the source (the profile isn't
  exposure-corrected), pile-up reaching ~16″ in the brightest (045), maybe
  the 3C 273 jet at 13–22″ (a wedge test was inconclusive, confounded by
  exposure). That per-observation structure, not pile-up, is what the
  flags picked up; v006 in the profile gave the same radii.
- **The fitter freezes Γ at `--defgamma` 2.0 below `--mingamma` 200
  counts.** For 3C 273 (Γ ≈ 1.5) that biases the 1 keV flux of short PC
  snapshots such as 043 (82 counts at its new 14.2″). For 1.2.8.
- **Overrides reach extraction only through `_pileup.txt`.** Extraction
  never reads `pileup_overrides.txt`, though `parallel_extract.py` symlinks
  it into each chunk dir as if it did. Documented ("re-run 5a"); for
  1.2.7 to decide whether extraction should read the file directly.
- **Extraction ignores the pile-up radius, overrides included, below 0.5
  ct/s** of whole-field rate (`get_inner_radius`). Documented.
- **`badstripe` is parsed by the extractor and never used.** Documented as
  a note-only column.
- **`--sosta` reads `source_extraction_<OBSID>.txt` files that nothing in
  the pipeline writes.** Documented; a candidate for removal.
- **First images in `docs/img/`**: Step 5 profile and source plots, as
  256-colour PNGs (25–32 KB each) that keep the plot colours exact.

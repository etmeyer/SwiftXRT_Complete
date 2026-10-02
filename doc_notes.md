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

- **`xrt_pipeline.py` wrapper fails where bare `xrtpipeline` succeeds (env).**
  Running the wrapper headless on amorgos aborts in the `prefilter`/`xrtfilter`
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

- **`xrt_pipeline.py` suppresses xrtpipeline's own stdout/stderr.** On failure
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
  hardening could run each OBSID from a private scratch cwd.

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

- **Bug B race itself is not fixed**, only surfaced: `parallel_extract.py` now
  checks that every included OBSID got a fresh `_grp.pha` and exits non-zero
  listing the misses, but the underlying intermittent drop is undiagnosed.
- ~~`requests` missing from the `heasoft` conda env~~ — installed
  2026-10-02 (`requests` + `charset-normalizer` only, with openssl/certifi/
  ca-certificates pinned to their existing versions); the "download before
  `heainit`" workaround was removed from the docs.
- **`swift_xrt_download.py --obsid X --list-only` downloads anyway.**
  `--list-only` is ignored with `--obsid` (26 files landed in
  `./swift_xrt_data`). Docs say `--obsid` skips the catalog query, but
  `--list-only` should still mean "don't download".
- **xselect silently loses output for long paths.** At ~140 characters it
  printed "Wrote spectrum to ..." but no file appeared (and it crashed in
  `xsl_exit`); ~100 characters work. Documented as a gotcha only.
- OBSID 041-style 0-count WT spectra and the King-profile low-count guard
  (from the May bug log) are untouched.

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

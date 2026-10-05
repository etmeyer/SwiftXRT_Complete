# Swift XRT Spectral Analysis Pipeline — Documentation

A walk-through of the pipeline, step by step. Each page covers one stage and
includes mechanism, inputs/outputs, gotchas, and space for your own notes as
you work through it.

For installation and environment setup, start at [Step 1](01-setup.md).
For a one-source quick-start example, see the [README](../README.md).

## Workflow

| Step | Doc | Scripts | Terminal | What it does |
| ---- | --- | ------- | -------- | ------------ |
| 1 | [Setup](01-setup.md) | (env + `swift_xrt_doctor.py`) | both | Install the pipeline and verify each terminal |
| 2 | [Download data](02-download.md) | `swift_xrt_download.py` | either | List, filter, and download Swift XRT observations |
| 3 | [Run xrtpipeline](03-xrtpipeline.md) | `xrt_pipeline.py` | HEASoft | Produce cleaned level-2 event files |
| 4 | [Survey observations](04-survey.md) | `swift_xrt_summary.py` | either | Per-OBSID targets, pointing offsets, modes, exposures, orbits |
| 5 | [PC-mode inspection](05-pc-inspection.md) | `swift_xrt_king_profile.py`, `swift_pc_source_viewer.py`, PC master table | either | Pile-up + source images + PC master selection table |
| 6 | [WT-mode inspection](06-wt-inspection.md) | `swift_wt_summary_viewer.py`, `make_wt_master_table.py` | either | WT 1D-strip profiles + WT master selection table |
| 7 | [Extract spectra](07-extract.md) | `parallel_extract.py` (or `swift_xrt_extract_spectra.py`) | HEASoft | Source and background spectra, ARFs, RMFs, grouping |
| 8 | [Fit and plot](08-fit-and-plot.md) | `parallel_fit.py` (or `swift_xrt_fit_spectra.py`), `plot_lightcurve.py` | CIAO, with `--caldb` | Spectral fits + νFν light curve |

**Two terminals.** Steps 2–7 run in a HEASoft terminal
(`setup_swiftxrt; heainit`; never `ciao`), and Step 8 in a CIAO terminal
(`setup_swiftxrt; ciao`). Why, and what goes wrong otherwise:
[Step 1 — Two terminals](01-setup.md#two-terminals).

## How to use these docs

Read top-to-bottom in order — each page assumes the previous step's outputs.
Each page ends with a "Notes" section where you can drop your own
observations as you work through; use it freely.

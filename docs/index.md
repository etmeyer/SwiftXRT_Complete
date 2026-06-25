# Swift XRT Spectral Analysis Pipeline — Documentation

A walk-through of the pipeline, step by step. Each page covers one stage and
includes mechanism, inputs/outputs, gotchas, and space for your own notes as
you work through it.

For installation and environment setup, start at [Step 1](01-setup.md).
For a one-source quick-start example, see the [README](../README.md).

## Workflow

| Step | Doc | Scripts | What it does |
| ---- | --- | ------- | ------------ |
| 1 | [Setup](01-setup.md) | (env + `swift_xrt_doctor.py`) | Install the pipeline and verify the environment |
| 2 | [Download data](02-download.md) | `swift_xrt_download.py` | List, filter, and download Swift XRT observations |
| 3 | [Run xrtpipeline](03-xrtpipeline.md) | `xrt_pipeline.py` | Produce cleaned level-2 event files |
| 4 | [Survey observations](04-survey.md) | `swift_xrt_summary.py` | Per-OBSID summary of modes, exposures, orbits |
| 5 | [PC-mode inspection](05-pc-inspection.md) | `swift_xrt_king_profile.py`, `swift_pc_source_viewer.py`, PC master table | Pile-up + source images + PC master selection table |
| 6 | [WT-mode inspection](06-wt-inspection.md) | `swift_wt_summary_viewer.py`, `make_wt_master_table.py` | WT 1D-strip profiles + WT master selection table |
| 7 | [Extract spectra](07-extract.md) | `parallel_extract.py` (or `swift_xrt_extract_spectra.py`) | Per-OBSID source/background extraction + grouping |
| 8 | [Fit and plot](08-fit-and-plot.md) | `parallel_fit.py` (or `swift_xrt_fit_spectra.py`), `plot_lightcurve.py` | Spectral fits + νFν light curve |

## How to use these docs

Read top-to-bottom in order — each page assumes the previous step's outputs.
Each page ends with a "Notes" section where you can drop your own
observations as you work through; use it freely.

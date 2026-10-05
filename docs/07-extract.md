# Step 7 — Extract spectra

> Run in the **HEASoft terminal** (`setup_swiftxrt; heainit`, never
> `ciao`); see [Step 1 — Two terminals](01-setup.md#two-terminals).

For every observation with `include` = `yes` in the master tables from
[Step 5](05-pc-inspection.md) and [Step 6](06-wt-inspection.md),
`swift_xrt_extract_spectra.py` makes a spectrum ready to fit:

- a source and a background spectrum, cut from the cleaned event file;
- the ARF, the effective area for that source region;
- the RMF, the detector's energy response, chosen from CALDB;
- a grouped copy of the source spectrum whose header points to the other
  three. [Step 8](08-fit-and-plot.md) fits this file.

`parallel_extract.py` does the same, several observations at a time. Run
either from inside `XRT_output`. Each observation takes a few seconds;
`parallel_extract.py` extracted the 87 spectra of 3C 273's first epoch in
about 20 s on 16 workers.

## What runs

```bash
cd XRT_output
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524      # one at a time
parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 16    # or in parallel
```

For each observation the script:

1. **Sets the regions** and writes them as ds9 region files, `_src.reg`
   and `_bkg.reg`.
   - **PC**: around the source position in `_pileup.txt` (Step 5a). The
     source region is a 47″ circle, or an annulus from the pile-up radius
     to 47″. The background is an annulus 100–160″ from the source.
   - **WT**: the position and radii in `_wt_profile.txt` (Step 6a): a
     20-pixel (47″) source circle and an 80–120 pixel background annulus.
2. **Extracts the spectra** with xselect: `_src.pha` and `_bkg.pha`. For
   WT it then replaces BACKSCAL, the region size used to scale the
   background, with the strip lengths Step 6 explains: 40 and 39 pixels.
3. **Makes the ARF** with xrtmkarf, from the exposure map xrtpipeline made
   for the same event file in [Step 3](03-xrtpipeline.md).
4. **Picks the RMF.** xrtmkarf chooses it from CALDB for the observation's
   mode and date.
5. **Groups the spectrum** with grppha into `_grp.pha`: channels below
   0.3 keV are marked bad, the background, ARF and RMF are written into its
   header, and channels are combined into bins of at least 20 counts.
6. **Writes `_extraction.log`**: regions, counts, BACKSCAL, bins and
   response files.

It ends with a summary table. If any observation failed it lists them and
exits with status 1.

| Flag | Purpose |
| ---- | ------- |
| `--ra` / `--dec` | Required, though the positions come from Steps 5 and 6; the PC log records them |
| `--mode` | `pc`, `wt` or `both` (default) |
| `--rout` | PC source outer radius in arcsec (default 47), or `auto` |
| `--bkg-inner` / `--bkg-outer` | PC background annulus in arcsec (default 100 and 160) |
| `--wt-srcrad`, `--wt-bkginner`, `--wt-bkgouter` | WT radii in pixels (default 20, 80, 120), used only where there is no `_wt_profile.txt` |
| `--mincounts` | Minimum counts per bin (default 20) |
| `--pctable` / `--wttable` | Master tables to read (default `pc_master_table.txt` and `wt_master_table.txt`) |

`parallel_extract.py` takes the same flags, plus:

| Flag | Purpose |
| ---- | ------- |
| `--nproc` | Number of observations extracted at once (default 8) |
| `--dryrun` | Show how the observations would be split up, without extracting |

## How it works

### PC regions

From the log of OBSID 073:

```
Source region  : annulus 15.5" - 47.0"
Bkg region     : annulus (546 cts)

Source counts  : 4236
Bkg counts     : 546
Exposure       : 6303.1 s
Count rate     : 2.574 ct/s

BACKSCAL src   : 1.117000e-03
BACKSCAL bkg   : 8.828000e-03
Bkg scaling    : 0.1265
```

The outer radius, 47″ (20 pixels), holds about 90% of a point source's
counts. The inner radius is the pile-up radius from Step 5, so this annulus
holds only about a fifth of them; the ARF accounts for the rest (below).
Step 7 uses a full circle when `_pileup.txt` has no pile-up radius, or when
its count rate is under 0.5 ct/s
([Step 5, gotcha 3](05-pc-inspection.md#gotchas)).

BACKSCAL records each region's area, and the background is scaled by their
ratio, here 0.1265. If the annulus holds fewer than 5 counts (for
example, when the source is near the edge of the CCD), the script shrinks its
outer radius to 140″, then 120″, and failing that uses a 60″ circle 120″
from the source.

For a source as bright as 3C 273 the annulus is not empty sky. In 073 its
surface brightness is about 8 times that 300–400″ from the source, so most
of its 546 counts come from the source's PSF wings. Scaled to the source region they are 69
counts, 1.6% of the source's, so the subtraction lowers piled-up PC fluxes
by about 1.5%. With a full 47″ circle the effect is a few tenths of a
percent.

`--rout auto` sets the outer radius where the signal-to-noise ratio stops
growing by 1% per 2″ step, between 20″ and 70″. For the five PC
observations here it chose 30–44″.

### WT regions

Step 7 takes the source position, the radii and the plate scale from
`_wt_profile.txt`, and computes BACKSCAL from the radii as Step 6 does:
2 × 20 = 40 for the source and 120 − 80 − 1 = 39 for the background. An
observation without a `_wt_profile.txt` gets the `--wt-*` radii, and the
median position of all the events in its file as the source position, a
rough guess.

### The ARF

xrtmkarf builds the ARF (effective area against energy) for the source
region from:

- the CALDB effective area of the telescope and CCD at the source's position
  on the detector, which falls off away from the centre of the field;
- the exposure map, which lowers it by the share of the region's exposure
  lost to bad columns, bad pixels and the edges of the readout window;
- a PSF correction (`psfflag=yes`), which scales it by the fraction of a
  point source's counts that fall inside the region.

So the ARF of 073's 15.5–47″ annulus is 27 cm² at 1.5 keV, against 146 cm²
for the 47″ WT circle of OBSID 029, and both fits give the flux of the
whole source.

The exposure map has to be the one xrtpipeline made for the same event file,
`<stem>_ex.img`. An OBSID often has several, for the WT settling segment,
the pointed segment and each PC window. Before October 2026 Step 7 used the
first one it found, which for 34 of the 82 WT observations of epoch 1 was
the wrong one: re-extracted with the right map, their fluxes changed by
−24% to +13%. The log records the map used:

```
ARF            : sw00035017073xpcw3po.arf
Exposure map   : sw00035017073xpcw3po_ex.img
RMF            : $CALDB/data/swift/xrt/cpf/rmf/swxpc0to12s6_20090101v014.rmf
```

### The RMF

xrtmkarf picks the RMF from CALDB for the mode, the event grades and the
date. Its name says which one: `swxpc0to12s6_20090101v014.rmf` is for PC
mode, grades 0–12, substrate voltage 6 V, valid from 1 January 2009,
version 14. WT spectra get `swxwt0to2s6_…`. In the test set the 2009–2010
observations got the 2009 RMF, 080 (2011) the 2011 one and 00091742013
(2014) the 2013 one. Observations from before the 2007 change of the CCD's
substrate voltage get `s0` RMFs.

The spectrum records the RMF as a `$CALDB/...` path. Step 8 runs in the
CIAO terminal, where `$CALDB` is Chandra's, so it needs
`--caldb /opt/CALDB` to find it.

### Grouping

grppha makes `_grp.pha` with these commands:

```
bad 0-29
chkey backfile sw00035017073xpcw3po_bkg.pha
chkey ancrfile sw00035017073xpcw3po.arf
chkey respfile $CALDB/data/swift/xrt/cpf/rmf/swxpc0to12s6_20090101v014.rmf
group min 20
```

Channels 0–29 (below 0.3 keV) are marked bad. Channels are then combined,
from the bottom up, into bins of at least 20 counts; the last channels,
which never reach 20, are marked bad too. The background and ARF are
recorded by file name only, so they have to stay in the same folder as
`_grp.pha`.

### The summary table

For the test set used on these pages (5 PC and 3 WT observations):

```
  OBSID          File                           Mode Date-Obs      Rin" Rout"  SrcCts  BkgCts   Exp(s)  Bins  RMF
  --------------------------------------------------------------------------------------------------------------------------------
  00035017018    sw00035017018xpcw3po             PC 2009-01-12    16.1  47.0    1773     232   2544.7    75  swxpc0to12s6_20090101v014.rmf
  00035017043    sw00035017043xpcw3po             PC 2009-05-17    13.7  47.0      84      14    122.4     4  swxpc0to12s6_20090101v014.rmf
  00035017073    sw00035017073xpcw3po             PC 2010-12-09    15.5  47.0    4236     546   6303.1   161  swxpc0to12s6_20090101v014.rmf
  00035017080    sw00035017080xpcw3po             PC 2011-02-19    13.8  47.0      50       8     74.9     2  swxpc0to12s6_20110101v014.rmf
  00091742013    sw00091742013xpcw3po             PC 2014-02-05    15.2  47.0     761      76   1098.8    36  swxpc0to12s6_20130101v014.rmf
  00035017029    sw00035017029xwtw2po             WT 2009-02-09     0.0  47.1    5602      79   1342.4   207  swxwt0to2s6_20090101v015.rmf
  00035017037    sw00035017037xwtw2po             WT 2009-03-01     0.0  47.1    3710      57    960.6   142  swxwt0to2s6_20090101v015.rmf
  00035017080    sw00035017080xwtw2po             WT 2011-02-19     0.0  47.1    3878      77    972.4   145  swxwt0to2s6_20110101v015.rmf
```

`Rin"` and `Rout"` are the source region radii in arcsec (`Rin"` 0 means a
full circle). `SrcCts` and `BkgCts` are the counts in each region, before
scaling. `Bins` counts the usable bins. The short PC snapshots have few: 043
(84 counts) has 4 and 080 (50 counts) has 2. Step 8 skips spectra with
fewer than 3 bins and fits the photon index only above 200 counts.

### Running in parallel

`parallel_extract.py` splits each master table into at most `--nproc`
chunks and runs `swift_xrt_extract_spectra.py` on each, `--nproc` at a time.
Each chunk runs in its own folder, `XRT_output/xrt_chunkNN_*`, which holds
links to the OBSID folders, so the results land in the real OBSID folders.
The chunk folders are removed afterwards. Each chunk's full output goes to
`parallel_extract_logs/chunkNN_pc.log` or `_wt.log`.

At the end it checks that every included observation has a `_grp.pha`
written during this run. For the epoch-1 data:

```
PC: 6 observations → 6 chunks
WT: 81 observations → 14 chunks
Total: 20 chunks across 16 workers
...
Done. 20 chunks succeeded, 0 failed.
All 87 grouped spectra written.
```

When an observation fails (here a table row for data that don't exist), it
shows the end of that chunk's output, lists what is missing, and exits with
status 1:

```
  Chunk 03 [WT]: FAIL [1/4 done]
    No spectra were successfully extracted.
      FAILED (1): sw00035017099xwtw2po
    full log: parallel_extract_logs/chunk03_wt.log
  ...
Done. 3 chunks succeeded, 1 failed.

ERROR: 1 of 4 observations have no new _grp.pha:
    sw00035017099xwtw2po
  Re-run them (e.g. swift_xrt_extract_spectra.py with a master table listing only these) and check the output.
Failed chunks:
  Chunk 03: exit code 1
```

## Inputs and outputs

**Inputs:**

- the cleaned event files and exposure maps in `XRT_output` (Step 3);
- `_pileup.txt` for each PC file (Step 5a); a PC observation without one
  fails;
- `_wt_profile.txt` for each WT file (Step 6a);
- `pc_master_table.txt` and `wt_master_table.txt` (Steps 5c and 6b, edited
  in 6c).

**Outputs**, in each OBSID folder:

```
XRT_output/
├── parallel_extract_logs/                 per-chunk logs (parallel_extract.py)
└── 00035017073/
    ├── sw00035017073xpcw3po_src.reg       source region (ds9)
    ├── sw00035017073xpcw3po_bkg.reg       background region (ds9)
    ├── sw00035017073xpcw3po_src.pha       source spectrum
    ├── sw00035017073xpcw3po_bkg.pha       background spectrum
    ├── sw00035017073xpcw3po.arf           ARF
    ├── sw00035017073xpcw3po_grp.pha       grouped spectrum   ← read by Step 8
    └── sw00035017073xpcw3po_extraction.log
```

A run of `swift_xrt_extract_spectra.py` also leaves xselect's `xselect.log`
in `XRT_output`.

## Common variants

```bash
cd XRT_output

# Everything, one observation at a time
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524

# Everything, 16 at a time
parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 16

# One mode only
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524 --mode wt

# Redo a few observations: a table with the header and just their rows
head -1 pc_master_table.txt > redo_pc.txt
grep -E '00035017018|00035017073' pc_master_table.txt >> redo_pc.txt
swift_xrt_extract_spectra.py --ra 187.2779 --dec 2.0524 --mode pc --pctable redo_pc.txt

# See how parallel_extract.py would split the work
parallel_extract.py --ra 187.2779 --dec 2.0524 --nproc 16 --dryrun
```

Each run overwrites the outputs of every observation it extracts.

## Gotchas

1. **If you extracted before October 2026, re-run this step, then Step 8.**
   Many WT ARFs were built from the wrong exposure map (34 of 82 in 3C 273's
   first epoch, whose fluxes changed by −24% to +13% when re-extracted),
   and Steps 5 and 6 have since changed the PC pile-up radii and, far from
   the celestial equator, the source positions.

2. **Change the WT regions in Step 6, not here.** The `--wt-*` flags only
   apply to observations without a `_wt_profile.txt`. To change a WT region,
   re-run the WT viewer with `--srcrad`, `--bkginner` or `--bkgouter`, then
   this step.

3. **Pile-up settings come from `_pileup.txt`.** To change a radius, edit
   `pileup_overrides.txt` and re-run Step 5a
   ([Step 5, gotcha 1](05-pc-inspection.md#gotchas)). Below 0.5 ct/s Step 7
   ignores both and extracts a full circle.

4. **Each spectrum needs its Step 3 exposure map.** Without `<stem>_ex.img`
   (for example after `xrt_pipeline.py --createexpomap no`) there is no ARF
   and the observation fails.

5. **Nothing keeps other sources out of the background region.** Another
   X-ray source in the PC background annulus, or on the WT strip, is
   subtracted along with the background. For a crowded field, open the
   event file with `_bkg.reg` in ds9; for WT, look for a second peak in the
   Step 6 profile.

6. **Keep `--mincounts` at 20 unless you also change how Step 8 fits.**
   Step 8 fits these bins as they are, and with fewer counts per bin its
   fit becomes biased.

7. **Keep the path to `XRT_output` under about 150 characters**
   (`pwd | wc -c`). Longer paths leave xselect unable to name its working
   files, and every extraction fails: `parallel_extract.py` failed at 180
   characters, `swift_xrt_extract_spectra.py` at 200.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_

# Step 8 — Fit and plot

> Run in the **CIAO terminal** (`setup_swiftxrt; ciao`) and pass
> `--caldb /opt/CALDB`; see [Step 1 — Two terminals](01-setup.md#two-terminals).

`swift_xrt_fit_spectra.py` fits each grouped spectrum from
[Step 7](07-extract.md) with an absorbed power law in Sherpa. It writes the
fluxes and photon indices to `fit_results.txt` and plots them against time
in `flux_lightcurve.pdf`. `parallel_fit.py` does the same with several
workers, and `plot_lightcurve.py` redraws the plot from the table.

Run them from inside `XRT_output`. `parallel_fit.py` fitted the 87 spectra
of 3C 273's first epoch in about 30 s on 16 workers.

## What runs

```bash
cd XRT_output
swift_xrt_fit_spectra.py --nh 0.0179 --caldb /opt/CALDB             # one at a time
parallel_fit.py --nh 0.0179 --caldb /opt/CALDB --nproc 16           # or in parallel
```

For each observation with `include` = `yes` in the master tables, the
script:

1. **Loads** `_grp.pha` with the background spectrum, ARF and RMF its header
   names, finding the RMF in `--caldb`.
2. **Keeps 0.3–10 keV**, minus the bins Step 7 marked bad.
3. **Fits** a power law absorbed by the Galaxy: `tbabs × powerlaw`, with
   the Galactic column frozen at `--nh` and the photon index Γ and
   normalization free. With `--redshift` the model also includes absorption
   in the source itself (`ztbabs`, column free).
4. **Measures** the 0.3–10 keV flux and the flux density at 1 keV, with
   1σ errors on these and on Γ.
5. **Writes** `_sherpa_fit.log` in the OBSID folder.

The fits come in two passes. First, every spectrum with at least 200 counts
(`--mingamma`) is fitted with Γ free. Then the spectra with 40–199 counts
are fitted with Γ frozen at the median of those free fits. Spectra under 40
counts (`--mincounts`), or with fewer than 3 usable bins, are skipped.

At the end it writes the table and the plot. If any spectrum failed it
lists them and exits with status 1.

| Flag | Purpose |
| ---- | ------- |
| `--nh` | Galactic column density in 10²² cm⁻² (required) |
| `--caldb` | Where the Swift CALDB is; always `/opt/CALDB` in the CIAO terminal |
| `--redshift` | Source redshift; adds absorption in the source |
| `--model` | `absorbed` (`tbabs × ztbabs × powerlaw`) or `simple` (`tbabs × powerlaw`); default `absorbed` with `--redshift`, else `simple` |
| `--stat` | Fit statistic: `wstat` (default) or `chi2`, the method before October 2026 |
| `--abund` | Abundance table for `tbabs` (default `wilm`) |
| `--mingamma` | Spectra with fewer counts get a frozen Γ (default 200) |
| `--defgamma` | The Γ to freeze them at (default: the median of the free fits) |
| `--mincounts` | Spectra with fewer counts are skipped (default 40) |
| `--bkg` | `use` the background spectrum (default) or `none` |
| `--emin` / `--emax` | Energy range in keV (default 0.3 and 10) |
| `--modes` | `pc`, `wt` or `both` (default) |
| `--pctable` / `--wttable` | Master tables to read |
| `--output` / `--plot` | Output names (default `fit_results.txt` and `flux_lightcurve.pdf`) |
| `--nmax` | Fit only the first N observations (`swift_xrt_fit_spectra.py` only) |

`parallel_fit.py` takes the same flags except `--nmax`, plus `--nproc`
(number of workers, default 8) and `--dryrun`.

## How it works

### The model

`tbabs` is the absorption by the gas between us and the source. Its column,
`--nh`, comes from a survey of Galactic hydrogen (for example HEASARC's
nH tool); for 3C 273 it is 1.79 × 10²⁰ cm⁻², so `--nh 0.0179`. It is held
fixed, so only the power law is fitted.

`tbabs` is built for the element abundances of Wilms, Allen & McCray (2000),
`wilm`, and the fits use those. Sherpa's default is `angr` (Anders &
Grevesse 1989), which has more of the heavier elements per hydrogen atom and
so absorbs more for the same column. For 3C 273 the difference is small
(Γ 0.012 lower, 1 keV flux 1.5% lower with `wilm`). Refitting its spectra
with the column set to 3 × 10²¹ cm⁻², it was 0.2 in Γ and 20% in the 1 keV
flux.

### The fit statistic

The fit uses W-stat: the Poisson likelihood of the counts in the source
spectrum and in the background spectrum together, with the background level
in each bin worked out from both. (It is what XSPEC's `cstat` becomes when a
background spectrum is loaded.)

Fits before October 2026 used χ² on background-subtracted spectra instead.
χ² weights each bin by its own counts, so bins that came out low get more
weight, and the fit ends up about one count per bin too low. Step 7 groups
every spectrum to about 20–35 counts per bin however bright the source, so
that is a 3–4% bias at every count level, which is largest compared with
the error bars for the brightest spectra (Humphrey, Liu & Buote 2009). On
simulated spectra built from the fits of OBSIDs 029 (WT) and 073 (PC),
400 per row, whose true flux is known:

| Counts per spectrum | χ², 20-count bins | W-stat, same bins |
| ------------------- | ----------------- | ----------------- |
| ~100 | −26% (WT), −36% (PC), with Γ frozen at 2.0 | +5%, +3%, Γ free (scatter ±18–28%) |
| ~250 | −3%, −9% | +1%, +1% |
| ~1000 | −4%, −4% | +0.4%, +0.6% |
| ~5000 | −4%, −4% | +0.1%, 0.0% |
| ~20,000 | −3%, −3% | 0.0%, 0.0% |

W-stat does need a few counts in each bin. On 1-count bins it was 19% low
for the ~100-count WT spectra, because the WT background region is as small
as the source region and most of its bins are empty. Step 7's 20-count bins
avoid that.

### Spectra with few counts

Below 200 counts Γ is poorly measured, so it is frozen and only the
normalization is fitted. Until October 2026 it was frozen at 2.0. For
3C 273, whose Γ is 1.4–1.8, that made the 0.3–10 keV flux of ~100-count
spectra 26–36% too low. Now the value is the median Γ of the spectra fitted
with Γ free in the same run, or `--defgamma` if you give one. In the first
epoch the median of the 85 free fits is 1.565, and OBSID 043's 84-count PC
spectrum went from 8.8 × 10⁻¹¹ to 1.24 × 10⁻¹⁰ erg cm⁻² s⁻¹.

The 1 keV flux density barely depends on the frozen value (043's changed by
1%), because most of the counts are near 1 keV. That makes it the steadier
measure for a light curve with short observations in it.

### Fluxes and errors

- **`flux_band`**: the 0.3–10 keV energy flux of the fitted model, as
  observed, with the Galactic absorption (erg cm⁻² s⁻¹).
- **`flux_1keV`**: the power law's flux density at 1 keV, without the
  absorption (erg cm⁻² s⁻¹ Hz⁻¹). The plot shows νF_ν = ν × `flux_1keV`.
- **Errors** are all 1σ. Γ and `flux_1keV` come from Sherpa's `conf`, which
  lets the other free parameters adjust. The `flux_band` error is the spread
  of the flux over 500 random draws of the fitted parameters.

Fits before October 2026 gave 90% errors on Γ and `flux_1keV` but 1σ on
`flux_band`.

### Reading the output

The script prints a header with the settings, then each fit:

```
Model: tbabs * powerlaw
Galactic nH: 0.0179 x 10^22 cm^-2
Statistic: W-stat (source and background as Poisson data)
Abundances: wilm
...
  [4/8] [PC] 00035017073 / sw00035017073xpcw3po
    Counts: 4236  Exposure: 6303.1s  Date: 2010-12-09T20:27:20.839
    Free gamma fit (initial guess: 2.0)
    RMF: swxpc0to12s6_20090101v014.rmf
    ARF: sw00035017073xpcw3po.arf
    Gamma: 1.481 +/- 0.022
    Flux (0.3-10.0 keV): 1.718e-10 +/- 3.815e-12 erg/cm²/s
    F_ν(1keV): 1.363e-28 +/- 2.489e-30 erg/cm²/s/Hz
    wstat/dof: 175.3/159 = 1.10
...
  2 spectra under 200 counts: gamma frozen at 1.571 (median of 6 free fits)

  [7/8] [PC] 00035017043 / sw00035017043xpcw3po
    Counts: 84  Exposure: 122.4s  Date: 2009-05-17T05:35:00.337
    Low counts (84 < 200): freezing gamma=1.571 (median of 6 free fits)
    ...
    Gamma: 1.571 (frozen)
```

`wstat/dof` is the statistic at the best fit over the degrees of freedom.
With 20-count bins it is near 1 for a good fit (median 1.11 for the first
epoch, range 0.85–1.58). It is a rough guide, not a probability.

`fit_results.txt`, for the test set used on these pages:

```
# Swift XRT spectral fit results (Sherpa)
# Model: xstbabs * xspowerlaw
# Statistic: W-stat (source and background as Poisson data)
# Abundances: wilm; errors are 1 sigma
# Under 200 counts gamma is frozen at 1.571 (median of 6 free fits)
# flux_band in erg/cm2/s, flux_1keV in erg/cm2/s/Hz

OBSID          filename                       mode   ctrate   exp(s) DateObs                         MJD    flux_band    fband_err    flux_1keV       f1_err   gamma  gamma_err stat/dof
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
00035017018    sw00035017018xpcw3po             PC    0.697   2544.7 2009-01-12T12:44:29.779  54843.53090    1.777e-10    5.787e-12    1.573e-28    4.296e-30   1.586      0.035    1.351
00035017029    sw00035017029xwtw2po             WT    4.173   1342.4 2009-02-09T16:52:51.078  54871.70337    1.778e-10    3.238e-12    1.609e-28    2.499e-30   1.607      0.020    1.130
00035017037    sw00035017037xwtw2po             WT    3.862    960.6 2009-03-01T23:15:32.352  54891.96912    1.668e-10    3.833e-12    1.434e-28    2.800e-30   1.556      0.025    1.019
00035017043    sw00035017043xpcw3po             PC    0.686    122.4 2009-05-17T05:35:00.337  54968.23264    1.235e-10    1.345e-11    1.077e-28    1.222e-29   1.571   (frozen)    1.268
00035017073    sw00035017073xpcw3po             PC    0.672   6303.1 2010-12-09T20:27:20.839  55539.85232    1.718e-10    3.815e-12    1.363e-28    2.489e-30   1.481      0.022    1.102
00035017080    sw00035017080xwtw2po             WT    3.988    972.4 2011-02-19T00:34:54.775  55611.02425    1.699e-10    3.891e-12    1.522e-28    2.880e-30   1.598      0.025    1.071
00091742013    sw00091742013xpcw3po             PC    0.693   1098.8 2014-02-05T15:00:38.015  56693.62544    1.731e-10    8.804e-12    1.372e-28    5.819e-30   1.480      0.052    1.372
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

`ctrate` is the counts in the source region per second. For a piled-up PC
source that region leaves out the core, so `ctrate` is only part of the
source's rate. `DateObs` and `MJD` are the middle of the observation. The
080 PC spectrum is missing: it has 2 usable bins and was skipped.

### The light curve

![Light curve of 3C 273's first epoch, 2009 to mid-2011: the upper panel shows νFν at 1 keV between about 2.6e-11 and 6e-11 erg/cm²/s, varying by a factor of two over months, with WT points as orange diamonds and PC points as navy circles; the lower panel shows Γ between 1.4 and 1.8](img/step8_lightcurve_e1.png)

The upper panel is νF_ν at 1 keV and the lower panel Γ, against time, with
PC as circles and WT as diamonds. A frozen Γ is drawn in grey (PC) or pale
orange (WT) and has no error bar. Above, the PC point in May 2009 with the
large error bar is OBSID 043.

### Running in parallel

`parallel_fit.py` splits the observations into up to `--nproc` chunks and
runs `swift_xrt_fit_spectra.py` on each, in its own folder
(`XRT_output/fit_chunkNN_*`, removed afterwards). It does the two passes
across all chunks: first the spectra with at least 200 counts, then it
takes the median Γ of all those fits and fits the rest with it.

```
Total: 87 observations: 85 with gamma free, then 2 under 200 counts with gamma frozen

Launching 15 fitting workers...
...
2 spectra under 200 counts: gamma frozen at 1.565 (median of 85 free fits)
Launching 2 fitting workers...
...
Merged 86 results into fit_results.txt
Generating light curve plot...

Done. 17 chunks succeeded, 0 failed.
```

It merges the chunks' tables into one `fit_results.txt`, sorted by OBSID,
and draws the plot with `plot_lightcurve.py`. A chunk that fails is
reported with the end of its output, its successful fits are still merged,
and the script exits with status 1.

## Inputs and outputs

**Inputs:** the grouped spectra with their backgrounds and ARFs (Step 7),
the master tables (Steps 5c, 6b, 6c), and the Swift CALDB.

**Outputs:**

```
XRT_output/
├── fit_results.txt                       the table above
├── flux_lightcurve.pdf                   the light curve
└── 00035017073/
    └── sw00035017073xpcw3po_sherpa_fit.log
```

## Common variants

```bash
cd XRT_output

# With absorption in the source, at its redshift
parallel_fit.py --nh 0.0179 --redshift 0.158 --caldb /opt/CALDB --nproc 16

# Freeze Γ for the short spectra at a value of your own
parallel_fit.py --nh 0.0179 --caldb /opt/CALDB --nproc 16 --defgamma 1.6

# The method before October 2026, for comparison (errors are 1σ either way)
swift_xrt_fit_spectra.py --nh 0.0179 --caldb /opt/CALDB --stat chi2 \
    --abund angr --defgamma 2.0 --output fit_results_old.txt --plot lc_old.pdf

# Redraw the plot: axis limits, reference lines, title
plot_lightcurve.py --xlim 2009 2011.6 --ylim_flux 2e-11 7e-11 \
    --times 2010.0 --tlabels "Flare" --title "3C 273" --output lc.pdf
```

`plot_lightcurve.py` also takes `--input`, `--ylim_gamma`, `--tcolor`,
`--figsize` and `--dpi`; an `--output` ending in `.png` makes a PNG.

## Gotchas

1. **CIAO terminal, with `--caldb /opt/CALDB`.** Sherpa exists only in
   CIAO, and CIAO's `$CALDB` is Chandra's. Without `--caldb` the scripts stop
   and say so.

2. **If you fitted before October 2026, re-run this step** (after Steps 5–7
   if they changed). The fit statistic, abundances, the frozen Γ and the
   errors have all changed. For 3C 273's first epoch the band fluxes rose by
   3–7% (median 4.6%), the 1 keV fluxes by 1–5% (median 2.3%), Γ fell by
   0.018 (median), and the Γ and 1 keV errors are 1.6 times smaller (1σ
   instead of 90%).

3. **`flux_band` includes the Galactic absorption; `flux_1keV` does not.**
   For 3C 273 the absorption is small, but for a source behind a large
   column the two differ a lot below 2 keV.

4. **The frozen Γ comes from the spectra in the same run.** With
   `--modes pc` it is the median of the PC fits only; with few bright
   spectra it is a noisy number, and with none it falls back to 2.0 (the
   output says so). Use `--defgamma` in those cases.

5. **Each run overwrites `fit_results.txt` and `flux_lightcurve.pdf`.** Use
   `--output` and `--plot` to keep a second version.

6. **A `wstat/dof` well above about 1.5 is worth a look** at the spectrum and
   its residuals in Sherpa or XSPEC: the source may not be a simple power law
   there, or the extraction may have gone wrong.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_

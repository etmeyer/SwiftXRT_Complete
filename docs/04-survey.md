# Step 4 — Survey observations

`swift_xrt_summary.py` reads the headers of every cleaned event file from
[Step 3](03-xrtpipeline.md) and prints what you have, OBSID by OBSID: the
target, how far Swift pointed from your source, the modes and segments in time
order, exposures, count rates, and orbit structure. It changes nothing on
disk. Use it to see what you are working with before inspecting the PC and WT
data in Steps 5 and 6.

Run it from inside `XRT_output`, in either terminal. It is plain Python and
takes about a second for 86 OBSIDs.

## What runs

The script finds every OBSID folder (a name of 8–11 digits) in the current
folder, reads each cleaned event file in it (`*_cl.evt`, all segments), and
prints one of two reports:

- **Detailed** (default): a table per OBSID with one row per event file, the
  mode sequence, pile-up warnings, and each file's orbit segments.
- **Compact** (`--compact`): one row per OBSID. This is the one to read first.

| Flag | Purpose |
| ---- | ------- |
| `--compact` | One row per OBSID instead of a table per OBSID |
| `--ra` / `--dec` | Your source position, J2000 decimal degrees. Pointing offsets are measured from it; without it, from each observation's own target |

## How it works

Everything comes from the headers and good-time tables of the cleaned event
files, so the survey is quick and never reads the events themselves:

- **Mode and segment** come from the file name: `pc` or `wt`, and `po`
  (pointed), `st` (settling) or `sl` (slew).
- **Exposure, start and end** are the `EXPOSURE`, `DATE-OBS` and `DATE-END`
  keywords. The number of events is the length of the events table.
- **Count rate** is events divided by exposure for the whole file. That
  includes background and every source in the field, so it is only a rough
  brightness.
- **Pile-up warnings** appear when a file's count rate passes about 0.5 ct/s
  in PC mode or 150 ct/s in WT mode. Step 5 measures PC pile-up properly.
- **Orbit segments** are the good-time intervals (GTIs). Swift orbits every
  ~96 minutes and loses the target behind the Earth for part of each orbit,
  so a long observation comes in pieces separated by gaps of an hour or more.
- **Target** is the `OBJECT` keyword: what Swift was observing.
- **Pointing offset** is the angle between where Swift pointed (`RA_PNT`,
  `DEC_PNT`) and your `--ra`/`--dec`, or, without them, the observation's own
  target (`RA_OBJ`, `DEC_OBJ`).

The compact report for the five test OBSIDs from Step 3 plus one pointing of
the nearby AGN SDSS J122933+015810 (`00035017030` failed in Step 3 on
purpose):

```
$ swift_xrt_summary.py --compact --ra 187.2779 --dec 2.0524

Found 6 OBSID directories.

======================================================================================================================================================
  COMPACT SUMMARY
  Sequence codes: 1=WT_SLEW  2=PC_SLEW  3=WT_SETTLING  4=PC_SETTLING  5=WT_POINTED  6=PC_POINTED
  Off(') = pointing offset in arcmin from --ra/--dec
======================================================================================================================================================
  OBSID          Date/Time               Total(ks)     ct/s  Slew_i  Slew_f  N_WT  N_PC  WT_exp(ks)  PC_exp(ks)  Orb    Seq Off(')  Target
  -------------- ---------------------- ---------- -------- ------- ------- ----- ----- ----------- ----------- ---- ------ ------  ------------------
  00035017018    2009-01-12T12:22:59.0765      2.556     2.87                     1     1       0.012       2.545    1     56    0.6  3C273
  00035017029    2009-02-09T16:41:38.7730      1.342     4.57                     1     0       1.342       0.000    1      5    3.1  3C273
  00035017037    2009-03-01T22:24:02.9188      0.961     4.24                     1     0       0.961       0.000    2      5    1.6  3C273
  00035017073    2010-12-09T14:33:47.3664      7.317     2.27      WT             0     1       0.000       6.303    5    136    3.0  3C273
  00091742013    2014-02-05T12:30:47.1663      1.350     2.11      WT             0     1       0.000       1.099    3    136    7.3  SDSSJ122933+015810
  -------------- ---------------------- ---------- -------- ------- ------- ----- ----- ----------- ----------- ---- ------ ------  ------------------
  TOTAL                                     13.526                                            2.315       9.947
======================================================================================================================================================
  Not in the table (no cleaned event files): 00035017030
```

| Column | Meaning |
| ------ | ------- |
| `Date/Time` | Start of the first segment (UTC) |
| `Total(ks)` | Exposure summed over all segments, in ks |
| `ct/s` | All events over the total exposure; a rough brightness |
| `Slew_i`, `Slew_f` | Mode of a slew segment at the start or end, if any |
| `N_WT`, `N_PC` | Number of pointed WT and PC event files |
| `WT_exp(ks)`, `PC_exp(ks)` | Pointed exposure in each mode |
| `Orb` | Most orbit segments (GTIs) in any pointed file |
| `Seq` | Modes in time order, using the codes in the legend; `136` is WT slew, WT settling, then PC pointed |
| `Off(')` | Pointing offset in arcmin |
| `Target` | What Swift was observing (`OBJECT`) |

The detailed report for one OBSID, here a 6 ks PC observation in five pieces:

```
====================================================================================================
  OBSID: 00035017073    Target: 3C273    Pointing offset: 3.0' from --ra/--dec
====================================================================================================
  #    Mode           DATAMODE     Start (UTC)            End (UTC)                  Exposure   Events     ct/s
  ---- -------------- ------------ ---------------------- ---------------------- ------------ -------- --------
  1    WT_SLEW        WINDOWED     2010-12-09T14:33:47.3664 2010-12-09T22:38:27.3849 957.1s (16.0m)      207     0.22
  2    WT_SETTLING    WINDOWED     2010-12-09T14:36:27.3631 2010-12-09T22:38:37.9846        56.3s      212     3.77
  3    PC_POINTED     PHOTON       2010-12-09T17:49:39.7863 2010-12-09T23:05:01.8926 6303.1s (1.75h)    16225     2.57
  ---- -------------- ------------ ---------------------- ---------------------- ------------ -------- --------
       TOTAL                                                                     7316.5s (2.03h)    16644

  Sequence: WT_SLEW -> WT_SETTLING -> PC_POINTED

  *** WARNING: PC_POINTED has 2.57 ct/s -- possible pile-up (threshold ~0.5 ct/s) ***
  ...
  PC_POINTED — 5 orbit segments (elapsed: 18922s, duty cycle: 33%):
    Segment 1: 1080.6s (18.0 min)  [gap: 0.0 min]
    Segment 2: 501.5s (8.4 min)  [gap: 69.6 min]
    Segment 3: 1582.1s (26.4 min)  [gap: 70.6 min]
    Segment 4: 1582.1s (26.4 min)  [gap: 69.6 min]
    Segment 5: 1582.1s (26.4 min)
```

The duty cycle is exposure over elapsed time. A `gap: 0.0 min` means one
orbit's data were split into two intervals, so `Orb` can count a few more
pieces than there were orbits.

**What to look for:**

- **Which OBSIDs have PC data** (`N_PC` 1) and which have WT data (`N_WT` 1).
  They go to [Step 5](05-pc-inspection.md) and [Step 6](06-wt-inspection.md).
  Some observations open with a few seconds of WT before switching to PC,
  like 018's 12 s; Step 6 excludes WT exposures under 20 s automatically.
- **Pile-up warnings on PC data.** Expected for a bright source like 3C 273;
  Step 5 measures how much of the PSF core to exclude.
- **A different `Target`.** The Step 2 cone search also returns pointings of
  other sources near yours, here SDSS J122933+015810. Your source is then
  off-axis by `Off(')`; you decide in Steps 5 and 6 whether to keep them.
- **Large offsets in WT mode.** The WT window is a strip about 8′ wide, so a
  pointing more than ~4′ from your source can leave it outside the window. In
  epoch 1, OBSID 041 is at 5.1′ (all others within 3.6′) and contains no
  source; Step 6 detects and excludes it.
- **OBSIDs not in the table.** Their xrtpipeline run failed; see
  [Step 3](03-xrtpipeline.md), gotcha 3.

## Inputs and outputs

**Inputs:** the `XRT_output` folder from Step 3 (run the script inside it),
and optionally your source position (`--ra`, `--dec`).

**Outputs:** none on disk; the report goes to the terminal. Redirect it to a
file to keep it.

## Common variants

```bash
cd XRT_output

# One row per OBSID: start here
swift_xrt_summary.py --compact --ra 187.2779 --dec 2.0524

# The wide table without line wrapping (arrow keys scroll sideways)
swift_xrt_summary.py --compact --ra 187.2779 --dec 2.0524 | less -S

# The detailed report, one screen at a time
swift_xrt_summary.py --ra 187.2779 --dec 2.0524 | less

# Keep a copy
swift_xrt_summary.py --compact --ra 187.2779 --dec 2.0524 > survey_compact.txt
```

## Gotchas

1. **Run it inside `XRT_output`.** The script reads compressed
   `*_cl.evt.gz` files too, so in `XRT_input` it summarizes the archive's own
   cleaned files instead of yours. The table looks normal and even lists
   OBSIDs whose Step 3 run failed. In a folder without OBSID folders it stops
   with `No OBSID directories found`.

2. **Pass `--ra`/`--dec`.** Without them, each offset is measured from that
   observation's own target, which says nothing about your source in a
   pointing of something else. The SDSS pointing above shows 2.0′ without
   them (from SDSS J122933) and 7.3′ with them (from 3C 273).

3. **Count rates are for the whole field.** In the SDSS pointing above, the
   2.50 ct/s PC rate and its pile-up warning come mostly from 3C 273, 7.3′
   away, not from the target. Treat the rates and warnings as hints.

4. **The compact table is about 150 characters wide.** Widen the terminal or
   use `less -S`, or the rows wrap and the columns become hard to read.

## Notes

<!-- Eileen: drop observations here as you walk through. Format suggestion:
     - 2026-MM-DD — observation / gotcha / "I ran this on X and Y happened"
-->

_(no notes yet)_

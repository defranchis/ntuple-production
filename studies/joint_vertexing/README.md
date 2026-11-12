# Joint PV + SV vertexer for ALEPH (research prototype)

An offline Python package (`jv`) that finds and fits the primary vertex (PV) and the secondary vertices (SVs) of an
event together: every track is shared softly between the PV, the SV candidates and an outlier term in one annealed fit.

**Status:** research prototype, not part of the stage1 pipeline. It runs on files harvested from stage1 outputs. The
configurations were tuned on simulated events (MC) only; `configs/` holds the frozen PV + SV fitter ("frozen final",
`final_*.json`) and the configurations with the V0 class (`v0_*.json`, section "V0 class").

## Algorithm (frozen final, stage by stage)

Abbreviations: IP = impact parameter; HF = heavy flavour (b or c hadron); V0 = neutral two-prong decay (K0s, Lambda,
photon conversion); VDET / ITC = ALEPH vertex detector / inner tracking chamber; svn = the SV collection of stage1.

1. **S0 inputs:** stage1 baseline tracks with finite parameters, minus the daughters of tight stage1 V0s; each track gets
   an upper bound on its production radius from its innermost VDET/ITC hits (hit-radius rule, applied in the later
   stages).
2. **S1 reference PV:** jet-projected PV with the beam-spot constraint, as in the ALEPH routine QFNDIP (each track
   constrains the PV mainly perpendicular to its jet, where lifetime offsets do not pull); the beam-constrained fit with
   worst-track removal (as stage1's PV) is the fallback.
3. **Lifetime sign:** 3D IP significance of every track to the reference PV, signed by its jet axis.
4. **S2 + S3 candidates (seedless):** maxima of a vertex-density function of all S0 tracks (as in ZVTOP, the SLD
   topological vertex finder); each maximum takes the tracks closer to it than to the PV that fly forward; then a fit
   with worst-track removal, mass, flight significance, V0 veto, acceptance window (step 7) and a logistic candidate
   score; overlapping candidates merged.
5. **S4 joint fit:** annealed multi-vertex fit of PV + candidates (temperatures 8, 2, 1) with Student-t soft weights; beam
   spot and a lifetime-sign prior on the PV term; SV compatibility includes the SV covariance; SV tracks must point along
   the flight direction; close vertices merged; SVs with < 2 tracks, outside the window or below the final-SV score dropped.
6. **S5 iteration (one round):** unassigned tracks and PV tracks with large positive significance seed two-track fits,
   grown as in the CMS inclusive vertex finder (IVF: a track joins if closer to the candidate than to the PV), scored,
   then S4 again.
7. **S6 acceptance window** (used in steps 4-6): distance to the PV 0.01-3 cm, pointing, position error along the momentum,
   <= 8 tracks, flight significance > 3 for two-track SVs.
8. **Cascade split:** an SV with 4-7 tracks becomes two vertices (B -> D) when the best bipartition lowers the chi2 by
   > 9, both halves are separated, displaced from the PV, aligned along the flight direction and pass the final-SV score.
9. **Output:** PV and SVs (position, covariance; SV mass and track count), per track its PV weight and hard assignment
   (weight > 0.5).

Other algorithms of the study are configurations of the same code: stage1 as stored (`nominal`), the first
end-to-end version (`hybrid_v0`), rows modelled on LCFIPlus (the vertex finder of the International Linear Collider
flavour tagging), a multi-vertex fit (MVF) and IVF, and the truth-primary oracle PV (MC truth selects the primary
tracks); see `jv/hybrid.py` (`DEFAULT`, `ROWS`) and `jv/run.py`.

## Inputs

Per file stem (e.g. `mc_ZM4212.41`), the same stem in two directories:
- the stage1 output of this repository (`<stem>.root`; PV, V0s, svn SVs, jets, track membership bits),
- the EDM4hep file stage1 was run on (`<stem>.root`; EDM4hep = the event data model of the Key4hep software): track
  states with their covariance (stage1 does not store it), hit counts and MC truth.

`jv.harvest` joins them by (run, event) on stage1's event list and writes one `<stem>.npz` per file (layout in its
docstring). The samples are defined in `jv/io.py` (`SAMPLES`): DEV = `mc_ZM4212.41-50` (all tuning), HELDOUT =
`mc_ZM4212.51-60` (evaluated once after the freeze), DATA = `data_Y15224.10-19`, FAST = the first 20000 events of
`mc_ZM4212.41`, TRAIN = `mc_ZM4212.50` (fit of the V0 tables), DEV9 = DEV without TRAIN, REG2 = `mc_ZM4212.41-42`
(regression runs).

## Environment

LCG (CERN software stack on /cvmfs) view LCG_107 (`x86_64-el9-gcc13-opt`: Python 3 with numpy, numba, uproot, awkward,
matplotlib). `run_py.sh` sets it up (override with `JV_LCG_SETUP=<setup.sh>`), single-threaded numerical libraries, a
local numba cache (`NUMBA_CACHE_DIR`, default `/tmp/jv_numba_$USER`), runs the module from this folder (relative paths
are relative to it), prints one status line and exits with the module's return code. `jv.harvest` reports a failed file as
a per-stem `FAILED` line and always exits 0: check its output for `FAILED`.
Directories are passed by environment variable (no defaults); give absolute paths, because `run_py.sh` changes into this
folder:

| variable | content |
|---|---|
| `JV_STAGE1` | stage1 outputs (`<stem>.root`), read by `jv.harvest` |
| `JV_RAW` | EDM4hep inputs of stage1 (`<stem>.root`), read by `jv.harvest` |
| `JV_INPUTS` | harvested `<stem>.npz`, written by `jv.harvest`, read by the other modules |
| `JV_RUNS` | run outputs `<tag>/<sample>/<algo>/*.npz`, written by `jv.run`, read by `jv.evaluate` and `jv.evaluate_v0` |

## Commands

```
bash run_py.sh 3600 jv.harvest 16 mc_ZM4212.41 mc_ZM4212.42 ...           # one npz per stem
bash run_py.sh 3600 jv.run 32 DEV mytag final_LR=configs/final_LR.json final_S1=configs/final_S1.json nominal
bash run_py.sh 3600 jv.evaluate 32 mytag DEV $JV_RUNS/mytag/rows_DEV.json --ref nominal --nboot 100
bash run_py.sh 600 jv.figs $JV_RUNS/mytag/rows_DEV.json $JV_RUNS/mytag/figs_DEV
bash reproduce.sh 32 DEV mytag                                                # run + evaluate the frozen final rows
bash reproduce.sh 32 DEV mytag2 v0                                            # the same + the nominal V0 configuration
bash run_py.sh 1800 jv.evaluate_v0 8 DEV $JV_RUNS/mytag/v0eval_DEV.json       # stage1 V0 lists vs the true-V0 table
```

`reproduce.sh` refuses an existing `$JV_RUNS/<tag>/<sample>` and writes `$JV_RUNS/<tag>/rows_<sample>.json`, with `v0`
also `$JV_RUNS/<tag>/v0_<sample>.json` (one row per S1 working point plus the reference row `stage1_tight`, see its header). `jv.evaluate_v0` and the `v0` mode
need inputs harvested with the true-V0 arrays (`t_v0t`).

The harvest writes a true-V0 table (`v0t_*`, MC; layout in `jv/harvest.py`) and follows GALEPH continuations
(restarted particles): a track linked to a continuation of a V0 daughter counts as that daughter (truth class 3, V0
daughter). Truth-based `jv.evaluate` numbers therefore differ from those on harvests of commit 0a61092, with which the
frozen final was evaluated. The fit reads no truth arrays; only the truth-seeded supply and the oracle rows do.

`jv.run` arguments: workers, sample (`DEV`, `FAST:20` = 20 events per unit for a smoke test), a new run tag, then
algorithms (`nominal`, `oracle`, a `jv.hybrid.ROWS` name, or `name=config.json` with overrides of `jv.hybrid.DEFAULT`;
a name starting with `oracle` runs the oracle PV with the config's `bs_sig`). Results per event do not depend on the
worker count; the output file boundaries do, so use a new tag for a different worker count. `JV_UPCACHE=1` reuses
stages S0-S3 between configurations that differ only in S4 keys (variant scans; the CPU times are then not valid).
`jv.evaluate` compares every algorithm under `<tag>/<sample>` with `--ref` on the same events (paired bootstrap;
metric definitions in its docstring). `jv.figs` as shipped does not plot the final rows: it draws its standard
figures only for the comparison rows (`nominal`, `hybrid_v0`, oracle, LCFIPlus-, MVF-, IVF-like). To plot `final_LR` /
`final_S1`, add them to `COL` and `LAB` (and `ORDER`) in `jv/figs.py` and call its functions with `algos=`.

## Frozen final configuration (`configs/final_LR.json`, `configs/final_S1.json`)

Overrides of `jv.hybrid.DEFAULT`; the two files differ only in the beam-spot constraint: LR = luminous region
(125 um x 16 um x 0.7 cm), S1 = the widths of the stage1 PV fit (200 um x 100 um x 2 cm; not to be confused with stage S1).

| key | value | meaning |
|---|---|---|
| `bs_sig` | LR: 0.0125, 0.0016, 0.7; S1: 0.02, 0.01, 2.0 | beam-spot constraint widths x, y, z (cm) on the PV fits of stages S1 and S4 |
| `pv_ref` | jetproj_hard | reference PV of stage S1: jet-projected PV (QFNDIP rule) |
| `supply` | vd | S2 + S3 candidates from vertex-density maxima (no pair seeds) |
| `vd_n_climb` | 4 | maximum fixed-point iterations per climb towards a density maximum |
| `sl_score_min` | -1.0 | cut on the candidate score (log-odds) |
| `sl_score_seeded` | true | the candidate score also applies to pair-seeded candidates (S5) |
| `sl_final_min` | 0.55 | cut on the final-SV score (log-odds) |
| `temps` | 8, 2, 1 | S4 annealing temperatures |
| `weight_form` | student | S4 soft weights Student-t instead of Gaussian |
| `student_nu` | 4.0 | degrees of freedom of the Student-t weights |
| `chi2_mode_sv` | vcov | track-SV chi2 with the SV covariance added to the track's |
| `chi2_cut_sv` | 16.0 | outlier constant (chi2) of the SV weights |
| `chi2_cut_neg` | 7.0 | PV outlier constant for tracks with negative lifetime sign (others: 5) |
| `pv_Tmin` | 4.0 | the PV is refitted with its weights evaluated at temperature max(T, 4) |
| `prior_s` | 2.0 | lifetime-sign prior: PV weight x 1 / (1 + (s / prior_s)^2), s = positive signed IP significance |
| `sv_forward` | 0.8 | minimum cos(track direction at the SV, flight direction) for SV tracks |
| `sv_min_vdet` | 0 | minimum VDET hits of SV tracks (0 = none) |
| `acc_dis_lo` | 0.01 | S6: minimum PV-SV distance (cm) |
| `acc_sigl` | 0.3 | S6: maximum SV position error along its momentum (cm) |
| `acc_cos` | 0.5 | S6: minimum cos(flight direction, SV momentum) |
| `acc_fsig2` | 3.0 | S6: minimum flight significance of two-track SVs |
| `split_casc` | true | cascade split of final SVs on |
| `split_order` | 0.98 | split: minimum cos between the near-to-far half direction and the PV-to-near direction |
| `split_pick` | valid | split: best chi2 among the bipartitions that pass the geometric cuts |

The two logistic scores use `jv/score_coefs.json` (candidates) and `jv/score_coefs_final.json` (final SVs), fitted on
MC truth of the FAST sample.

## V0 class (`configs/v0_S1.json`, `configs/v0_LR.json`)

With `v0_class` (`jv/v0.py`; default off = frozen final) K0S, Lambda and photon conversions (`jv/conv.py`) are two-track
vertex hypotheses inside the S4 joint fit: each candidate pair competes with the PV, the SVs and the outlier term for
its two tracks, with odds from a prior, mass, pointing, inner-hit and heavy-flavour (HF) term (model in the `jv/v0.py`
docstring). `jv.run` then also writes a V0 list per event (`v0_*`, `t_v0`; `jv/run.py`). Tables: `jv/v0_params.json`,
fitted on `mc_ZM4212.50`.

The nominal configuration is the frozen final plus:

| key | value | meaning |
|---|---|---|
| `pt_const` | stage1 | momentum constant 0.0029979 x 1.5 T as in stage1 (`legacy` = 0.0045 of the frozen final) |
| `remove_tight_v0`, `v0_veto` | false | required: otherwise stage1's tight V0 daughters are removed and V0-like pairs vetoed before the V0 class sees them |
| `v0_class`, `v0_conv` | true | V0 class, with photon conversions |
| `v0_prior`, `v0_prior_nm` | S1, S1_nomass (LR: LR, LR_nomass) | prior tables of the weight and of the weight without the mass term (`v0_w0`) |
| `v0_hf` | S1 (LR: LR) | HF term: displaced tracks at the V0 vertex and in the daughters' jets, decay-length density; the table holds the PV-pointing factor 1.75 |
| `v0_twin_sig`, `v0_twin_ext` | 3.0, 1.5 | twin rule: an SV within 3 sigma of the V0 vertex is the pair's own hypothesis unless its weights on other tracks sum to more than 1.5 |
| `v0_conv_rmin` | 4.5 | lower bound on the conversion radius from the beam axis (cm) |
| `v0_conv_params`, `v0_hit_table`, `v0_conv_hscale` | default, default, 1.3364 / 1.3615 | conversion pointing and inner-hit terms, track-term chi2 scales |
| `v0_qt_min` | 0.015, 0.015 | minimum Armenteros qT of K0S, Lambda candidates (GeV) |

Variations, each the nominal with one change: `v0_twinfrac_*` (twin rule without a parameter, `v0_twin_ext` frac, HF
table `*_twinfrac`), `v0_svpoint_*` (SV-pointing odds: `v0_pt_point`, `v0_pt_params`), `v0_decayangle_*` (HF term with
the |cos theta*| feature, table `*_decayangle`).

The stored weight `v0_w` is not calibrated as a probability. Working points of the nominal configuration, i.e. the
thresholds on `v0_w` at which its list has the purity of stage1's tight V0 list (conversions: of the tight list of
`jv/conv.py`) [measured, simulation, DEV; DEV includes `mc_ZM4212.50`, the file the V0 tables were fitted on]:

| beam spot | K0S | Lambda | conversions |
|---|---|---|---|
| S1 | 0.398 | 0.329 | 0.381 |
| LR | 0.399 | 0.329 | 0.380 |

Configuration errors raise when the configuration is built (`jv.hybrid.config`): a V0 configuration without the joint
fit or with a missing table (ValueError / KeyError), and an `oracle*` row name together with `v0_class` (`jv.run`).

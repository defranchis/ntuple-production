# Aleph

Version of the FCCAnalyses code that supports command line arguments, to be able to process both data and MC with the same script. Nightlies version of the key4hep stack is required for this. 

## Setup

Folow the steps described in the general [README](../README.md) at the top level to setup the code, then just `cd src`. 
 
<!-- ```bash
git clone https://github.com/Apranikstar/Aleph.git
cd Aleph
git submodule update --init --recursive
cd FCCAnalyses
fccanalysis build -j 8
cd ../src
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh #or compile and source the FCCAnalyses module.
```

Need nightlies because updated FCCAnalyses version after this commit is needed: https://github.com/HEP-FCC/FCCAnalyses/pull/474 -->

## Stage1: Produce ntuples

Note: Change the data fraction based on your needs.

### Run on MC:
```bash
fccanalysis run stage1.py -- --tag <version_tag>  --MCflavour <flavour_index>
```

Output files will be in: 
`/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/<year>/<mc_type>/stage1/<version_tag>/<flavour_name>.root`; set `ALEPH_OUT_DIR=<dir>` to write to `<dir>/wp1_stage1/<version_tag>/` instead (e.g. when `/eos/experiment` is read-only).

Fraction of events to process can be set via `--fraction <val>`, default is to process all events. 

`--year <year>` and `--MCtype <type>` are also supported command line arguments, currently we only have `1994` and `zqq` here. 

### Run on data:
```bash
fccanalysis run stage1.py -- --tag <version_tag> --doData 
```

Output files will be in the working directory (`--batch`: `/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/<year>/stage1/<version_tag>/`); set `ALEPH_OUT_DIR=<dir>` to write to `<dir>/wp2_data/<version_tag>/`.

`--year` and `--fraction` is also supported as an argument here. 

### Run on batch:
```
fccanalysis submit stage1.py -- --tag VXX-XX --MCflavour X --batch --chunks X
```

### Reconstruction modules: defaults and opt-outs

A flag-less `stage1.py` runs everything: the standalone primary-vertex fitter ([`analyzer_pvnew.h`](analyzer_pvnew.h)), the two-tier V0 module ([`analyzer_v0new.h`](analyzer_v0new.h)), the V0-first secondary-vertex module ([`analyzer_svnew.h`](analyzer_svnew.h)), the φ→K⁺K⁻ finder ([`analyzer_phikk.h`](analyzer_phikk.h)) and the D*→D⁰π finder ([`analyzer_dstar.h`](analyzer_dstar.h)), writing the `pv_*`, `v0n_*`/`v0njet_*`, `svn_*`/`svm_*`, `phikk_*` and `dstar_*`/`d0_*` branches. On MC the truth-matching branches (`truev0_*`, `v0c_*`, `v0n_class`, `truephi_*`, `truedstar_*`, `trued0_*`, ...) are added automatically; they are skipped under `--doData`.

The legacy code paths remain available as opt-outs:

| flag | meaning |
| --- | --- |
| `--oldPV` | legacy PV chain exactly as implemented in `FCCAnalyses`: `get_PrimaryTracks` + `VertexFitter_Tk`, and the origin-referenced track pre-selection instead of the beamspot-referenced one. No `pv_*` flag branches. Note that the beamspot constraint of the `get_PrimaryTracks` selection fit is passed in 10 µm units while its track parameters are read in cm, so that constraint is off by a factor 1000 and is effectively absent; the final `VertexFitter_Tk` fit is unaffected. |
| `--oldSV` | drop the standalone SV module: no `svn_*`/`svm_*` branches. |
| `--oldV0` | drop the two-tier V0 module: no `v0n_*`/`v0njet_*` and no V0 truth branches, and on MC no φ→KK or D* truth blocks either (`truephi_*`, `phikk_class`, `truedstar_*`, `trued0_*`, `dstar_class`, `d0_class`), since those share the V0 truth-matching input. Implies `--oldSV`, since the SV finder consumes the tight-V0 track veto. The φ→KK and D* finders themselves still run, with their Ks/Λ track veto empty (that veto list is a V0-module product). |
| `--noPhiKK` | skip the φ→K⁺K⁻ finder: no `phikk_*`/`truephi_*` branches. |
| `--noDstar` | skip the D*→D⁰π finder: no `dstar_*`/`d0_*`/`truedstar_*`/`trued0_*` branches. |

The selection is not configurable from the command line: every tuned value is a named `constexpr` in the header that uses it, and that declaration is its only definition. The paragraphs below describe the selections in words; the numbers quoted are those declarations.

### The primary-vertex selection

Tracks enter the fit through a pre-selection window on the impact parameters, `|D0| < 0.75 cm` and `|Z0| < 2 cm` (`PVN_D0_MAX`, `PVN_Z0_MAX`), referenced to the run beamspot — or to the origin under `--oldPV`, which leaves it off-centre by the beamspot offset in data. Track/vertex compatibility is then judged at `chi2max = 5` (`PVN_CHI2_MAX`); lowering it claims fewer tracks as primary and so leaves more to the secondary finders.

Both PV fits — the selection fit that prunes the track list and the final position fit — share one beamspot constraint, of Gaussian widths 200 um transverse in x, 100 um in y and 2 cm along the beam (`PVN_BS_SIGMA_X/Y/Z`, declared in cm). All of these live in [`analyzer_pvnew.h`](analyzer_pvnew.h), which is loaded for either PV chain because the legacy chain reads the same numbers, rescaled to its own unit convention.

### The secondary-vertex selection

The SV values live in [`analyzer_svnew.h`](analyzer_svnew.h) as the `SVN_*` constants, read directly by `findSVs`, so `stage1.py` passes none of them. A seed or a growth step is accepted at normalised vertex `chi2 < 10` (`SVN_CHI2`) with a per-track chi2 contribution capped at 5 (`SVN_TRK_CHI2`; <=0 disables the cap). A candidate is kept if it sits between 0.03 and 3 cm from the PV (`SVN_DIS_LO`, `SVN_DIS_HI`), its longitudinal vertex sigma stays under 0.10 cm (`SVN_SIGL_MAX`), and it points back to the PV with `cosPointing > 0.7` (`SVN_COS_POINT`). Growth stops at 8 tracks per candidate (`SVN_MAX_TRK`). Seeds are claimed best-chi2 first.

Two SV collections are written from the same event: `svn_*` runs on the tracks left after the tight-claimed V0 tracks are masked out (`SVN_MASK_MODE`), and `svm_*` is the unmasked control twin (`SVN_MASK_NONE`), for studying the V0/SV interplay.

### Other options

`--excludeRuns RUN [RUN ...]` (data only) vetoes the listed run numbers before any selection; `eventsProcessed` still counts the raw input. `--nthreads N` overrides the built-in RDataFrame thread count for the run. `--procfile NAME` restricts a run to one input file of the process (`QQB/NAME` on MC, the year directory on data), which is how one condor job maps to one file.

Environment overrides for local (non-`--batch`) runs: `ALEPH_RECLUS_DIR=<dir>` reads the input files from `<dir>` instead of `/eos/experiment` (a re-clustered copy lifts the RDataFrame thread cap set by the few TTree clusters of the raw files); `ALEPH_OUT_DIR=<dir>` writes the output under `<dir>` (see above).

### The two-tier V0 module

Standalone V0 (Ks/Λ) reconstruction in [`analyzer_v0new.h`](analyzer_v0new.h). It is *V0-first* by design: its track claims are meant to be consumed downstream (secondary-vertex finding runs on the unclaimed tracks), so the module optimises the correctness of each claim, not just the candidate list.

**Candidate building.** All opposite-charge secondary track pairs are vertexed with a single consistent fit (`VertexFitter_Tk`); every downstream quantity (momenta, invariant masses under both hypotheses, Armenteros–Podolanski (AP) variables, pointing) is derived from the refitted momenta at that vertex — there is no second fit.

**Two selection tiers, evaluated per hypothesis (Ks and Λ):**

- **Tight** — the adopted physics selection: mass window; momentum-tiered pointing cut (separate Ks and Λ ladders, `ksPointThr` / `lamPointThr`); a qT veto against photon conversions (Λ only); and a resolution-scaled AP-band cut around the exact kinematic locus — the band half-width follows the measured σ_ell(p) of each species (`ksBandThr`, `lamBandThrTight`; the Λ band is floored at low p and capped at the nominal ramp edge), plus common fit-quality (χ²) and displacement requirements.
- **Loose** — the ML-training tier: same windows/χ²/displacement, but flat pointing, a wider Ks AP band, a Λ AP band equal to a fixed fraction (0.8) of the ramp half-width floored at the tight band, and a relaxed Λ qT veto. It is a strict superset of tight and is what gets stored, so *any* tighter selection (including the historical ones) can be re-derived offline from any production.

**Hypothesis arbitration.** A pair passing both hypotheses is booked as the one whose invariant mass is closer to its window centre (normalised by the window half-width).

**Exclusive claiming, tight first.** Candidates claim their tracks exclusively in quality order: all tight candidates claim before any loose one, and within a tier the best-χ² candidate claims first. A track is claimed once; later candidates using it are dropped. This preserves the tight-only output exactly regardless of the loose tier, and defines the track set left to the SV finder.

**Stored flags and ML inputs.** `v0n_tight` re-derives the adopted tight package offline from the same single-source helpers (booking a candidate ≠ selecting it — variant productions still carry the adopted-package flag). Any tighter selection is re-derivable offline from the stored loose tier. `v0n_bandSig` and `v0n_massSig` store the AP-band and mass cut variables as signed pulls in resolution units for training; all other cut variables (cosPointing, pointSig, qT, χ², displacement, p, invM) are stored raw.

**Output branches.** Two groups, both written by default and both dropped by `--oldV0`:

- `n_v0n_event`, `v0n_pdg`, `v0n_invM`, `v0n_alpha`, `v0n_qt`, `v0n_chi2`, `v0n_dxyz`, `v0n_p`, `v0n_px/py/pz`, `v0n_cosPointing`, `v0n_pointSig`, `v0n_tight`, `v0n_bandSig`, `v0n_massSig`, `v0n_vx/vy/vz`, the vertex-fit covariance `v0n_cov_*`, the daughter joins `v0n_trk1_origIdx`/`v0n_trk2_origIdx` and their `v0n_trk{1,2}_dEdx_{pads,wires}_{value,error}` — event-order candidate quantities, independent of jet assignment.
- `n_v0njet_jets`, `n_v0njet_ks`, `n_v0njet_lambda` and `v0njet_*` — the new candidates pushed through the same per-jet assignment and jet-relative getters as the existing `v0_*` block, so old vs new is an apples-to-apples comparison at jet level.

Two conditional additions on top: the V0→nearest-SV pointing features `v0n_svnCosPoint`, `v0n_svnPointSig`, `v0n_svnIdx` need the SV module (dropped by `--oldSV`), and the truth-matching branches (`v0n_class`, `v0n_trueidx`, `v0n_pairmult`, `v0n_trackshared`, `v0n_trk1`, `v0n_trk2`, `truev0_foundnew_*`, `truev0_*`, `v0c_*`) are MC only — they are skipped under `--doData`.

**Further utilities in the headers.** [`analyzer_truth.h`](analyzer_truth.h) carries the MC truth-matching utilities (true-V0 finding, track↔MC index recovery, candidate truth classification) used to derive the tunings above; it is loaded unconditionally, because its truth-FREE candidate accessors (`candChi2`, `candDxyz`, `candP`, `candPcomp`, `candCosPointing`, `candVtxPos`) are also used on data.

### The φ→K⁺K⁻ module

Standalone φ(1020)→K⁺K⁻ reconstruction in [`analyzer_phikk.h`](analyzer_phikk.h), on by default and skipped with `--noPhiKK` (`phikk_*` branches; truth branches added on MC). It is an extension of the standalone V0 machinery, intended to deliver a **kinematically tagged kaon sample for dE/dx calibration** — so **no dE/dx quantity enters any selection**; the daughters' dE/dx measurements are stored, never cut on.

**Candidate building.** Track pairs are formed from the *full* baseline-selected track list — primary and secondary tracks alike, with no masking by the PV split, minus the tight-claimed V0 daughters (see below) — and vertexed with the same single consistent `VertexFitter_Tk` call as the V0 module (momenta rescaled once by the cm-as-mm factor 10). A pre-fit K⁺K⁻ mass window on the perigee momenta removes the bulk of the pair combinatorics before any fit. Same-charge pairs are reconstructed too and flagged (`phikk_same_sign`), as the data-driven combinatorial control. There is **no exclusive claiming**: a track may appear in several candidates.

**V0 veto.** Tracks already claimed as daughters of a tight Ks/Λ candidate are removed from the pairing pool, the same claim list the D∗ finder uses. It is a product of the V0 module, so under `--oldV0` the veto is simply empty and the finder still runs. Candidates still claim nothing themselves: within the surviving pool a track may appear in several φ candidates.

**Promptness is not required.** The φ has cτ ≈ 46 fm, but φ from b/c decays (B→φX, D_s→φπ) are genuinely displaced and must be kept, so no displacement window and no pointing cut are applied. |vtx − PV| and its 3D significance are stored (`phikk_dpv`, `phikk_dpvSig`) and never cut on. A *storage fiducial* `DPV_FID` (≈ 3× the truth fit-error tail) bounds the stored sample: near-collinear K⁺K⁻ pairs leave the vertex unconstrained along the flight direction, so a genuine φ can be reconstructed out to ≈ 17 cm.

**Armenteros–Podolanski for equal masses.** With equal daughter masses the AP locus is the ellipse (α/α_max)² + (q_T/p\*)² = 1 centred at α = 0, with p\* = √(m_φ²/4 − m_K²) = 126.9 MeV and α_max = p\*/(β E\*), E\* = m_φ/2. `phikk_bandEll` stores the left-hand side (1 on the exact locus). Note that for equal masses (α, q_T) together with the pair momentum determine the invariant mass exactly, so the AP band is a momentum-dependent reparametrisation of the mass window rather than independent information; it is stored, not cut on.

**Selection.** Not configurable from the command line: every value is a named `constexpr` in `analyzer_phikk.h`, namespace `AlephPhiKK`, read directly by `findPhiKK`, so `stage1.py` passes none of them. The stored sample is deliberately loose, so that any working point is re-derivable offline: K⁺K⁻ mass window 0.98–1.10 GeV (`M_LO`/`M_HI`, with a 20 MeV pre-fit margin `PRE_MARGIN`), vertex χ² < 25 (`CHI2_CUT`, loose sanity, ndf = 1), storage fiducial |vtx − PV| < 50 cm (`DPV_FID`), and a daughter momentum floor of 0.3 GeV (`P_MIN_DEF`, applied to the *perigee* momentum while the stored `phikk_trk*_p` is at-vertex — they differ by ~1% in the tail). No displacement, pointing, AP-band or track-quality cut is applied, and same-charge pairs are kept as the control sample.

**Working-point flags.** Two per-candidate integer flags are stored alongside the loose sample, evaluated on the stored quantities themselves (post-fit mass, at-vertex daughter momenta): `phikk_wp` = |m − m_φ| < `WP_DM` and both daughters p > `WP_PDAU` and |vtx − PV| < `WP_DPV`; `phikk_tight` = the same with |m − m_φ| < `TIGHT_DM` and both daughters σ(d0) < `TIGHT_SIGD0` (constants in `analyzer_phikk.h`, namespace `AlephPhiKK`). Both are **charge-blind** on purpose: the signal sample is `phikk_wp && !phikk_same_sign` and the same-charge control sample under *identical* cuts is `phikk_wp && phikk_same_sign` (likewise for `phikk_tight`), which is what makes the sideband subtraction well defined. They are labels, not cuts — no candidate is dropped by them.

**Stored per-daughter block.** Original-track index, **charge**, momentum, cosθ, d0/z0/σ(d0), nVDET/nITC hits, track χ²/ndf, whether the track was in the fitted primary set, and the pads/wires dE/dx value and error — everything the offline quality and purity scans need. On MC: `phikk_class` (both daughters are kaons from the same true φ→K⁺K⁻), the daughter MC PDG and its LUND mother PDG, plus a `truephi_*` block with the true φ→K⁺K⁻ list, its origin class (5 = from a b hadron, 4 = from c, 0 = other) and how many of its kaons were tracked. Truth provenance is decoded from `generatorStatus` = 10000·KS + LUND mother line, since the MCParticles daughter relations are empty in these files.

### The D*→D⁰π module

Standalone D*⁺→D⁰(K⁻π⁺)π⁺_slow reconstruction in [`analyzer_dstar.h`](analyzer_dstar.h), on by default and skipped with `--noDstar` (`dstar_*` and `d0_*` branches; truth branches added on MC). Like the φ→KK module it is an extension of the standalone V0 machinery whose purpose is a **kinematically tagged kaon sample for dE/dx calibration** — so **no dE/dx quantity enters any selection**; the daughters' dE/dx measurements are stored, never cut on. The mass difference Δm = m(Kππ_s) − m(Kπ) is the handle: its resolution is set by the slow pion alone, so a narrow Δm window together with the D⁰ mass window isolates a sample in which the track given the kaon mass really is a kaon.

**Two output collections.** The stand-alone **D⁰→Kπ** list (`d0_*`, one entry per track pair *and mass assignment*) needs no slow pion at all; lacking the Δm handle its purity has to come from displacement, pointing and the helicity angle, all of which are stored. The **D\*** list (`dstar_*`, one entry per D⁰ candidate × third track) repeats the D⁰ quantities and carries `dstar_d0idx`, the index of its parent entry in the D⁰ list. `n_d0_event` / `n_dstar_event` count them.

**Candidate building.** Track pairs are formed from the *full* baseline-selected track list — primary and secondary tracks alike, with no masking by the PV split — and vertexed with the same single consistent `VertexFitter_Tk` call as the V0 module (momenta rescaled once by the cm-as-mm factor 10). **Both mass assignments** of every opposite-charge pair are separate candidates: the kaon hypothesis is what defines the tag, so (a=K, b=π) and (a=π, b=K) are different objects, not a symmetry to be resolved. A pre-fit Kπ mass window on the perigee momenta removes the bulk of the pair combinatorics before any fit, and the two-track fits are memoised per pair, so a pair reached by both mass assignments and by several slow pions is fitted exactly once; `n_d0fits_event` counts the fits actually performed. There is **no exclusive claiming**: a track may appear in several candidates. Per candidate, `dstar_nsec`/`d0_nsec` record how many legs sit in the secondary pool, and each daughter carries `..._pool` (0 primary set, 1 secondary set, 2 neither) beside `..._isprim`.

**No three-track fit.** The D⁰ flies ≈0.6 mm while the slow pion comes from the D\* decay point, i.e. from the PV region, so a common three-track vertex would be wrong. The D\* is built from the D⁰ momenta **at the fitted D⁰ vertex** plus the slow pion's **perigee** momentum, and Δm is required below `DM_MAX`.

**Right-sign vs wrong-sign.** The slow pion of a true D*⁺ carries the charge of the pion from the D⁰ (both opposite to the kaon). Both slow-pion charges are kept and flagged (`dstar_rs` = 1 for right-sign); the wrong-sign combination has no signal and measures the combinatorial background under the Δm peak from the data themselves.

**Promptness is not required.** D\* from b decays give a D⁰ that does not point at the PV, so there is no displacement window and no pointing cut: `dpv` and `dpvSig` (3D significance vs the PV) are **stored, never cut on**, and only the wide storage fiducial `DPV_FID` (10 cm) bounds |vtx − PV|; `cosPoint` is stored too and enters the `dstar_tight` label alone, which rejects random pairs without requiring promptness. `cosThetaStar`, the cosine of the kaon direction in the D⁰ rest frame w.r.t. the D⁰ lab flight direction, is stored too: flat for a true two-body decay, peaked at |cos| = 1 for combinatorics.

**Selection.** Not configurable from the command line: every value is a named `constexpr` in `analyzer_dstar.h`, namespace `AlephDstar`, read directly by `findDstar`, so `stage1.py` passes none of them. What is *stored* is deliberately loose, so any working point is re-derivable offline: Kπ mass window 1.70–2.03 GeV (`M_LO`/`M_HI`, 50 MeV pre-fit margin `PRE_MARGIN`), D⁰ vertex χ² < 25 (`CHI2_CUT`, loose sanity), storage fiducial |vtx − PV| < 10 cm (`DPV_FID`), Δm ceiling 0.20 GeV (`DM_MAX`), and perigee momentum floors of 0.3 GeV for the K and π (`P_MIN`) and 0.1 GeV for the slow pion (`PS_MIN`). No track-quality prefilter is applied.

**V0 veto.** Tracks already claimed as daughters of a tight Ks/Λ candidate are removed from the pool before any pairing. The claim list is a product of the V0 module, so under `--oldV0` the veto is simply empty and the finder still runs.

**Working-point flags.** Labels, not cuts — no candidate is dropped by them, and they are evaluated on the stored post-fit quantities. `d0_loose` = |m(Kπ) − m_D⁰| < 60 MeV (`D0_LOOSE_DM`); `d0_tight` = |m(Kπ) − m_D⁰| < 30 MeV (`D0_TIGHT_DM`) plus |vtx−PV| significance > 3 (`D0_TIGHT_DPVSIG`), `cosPoint` > 0.99 (`D0_TIGHT_COSPOINT`), |cosθ\*| < 0.8 (`D0_TIGHT_COSSTAR`) and the shared daughter-momentum and χ² requirements. `dstar_loose` = |m(Kπ) − m_D⁰| < 50 MeV (`DS_LOOSE_DM`) and |Δm − 145.426 MeV| < 3.0 MeV (`DS_LOOSE_DDM`), with no requirement on the slow pion or the pointing. `dstar_tight` tightens those to 25 MeV (`DS_TIGHT_DM`) and 1.5 MeV (`DS_TIGHT_DDM`) and adds p(K) > 1 GeV (`TIGHT_PK`), p(π) > 1 GeV (`TIGHT_PPI`), D⁰ vertex χ² < 10 (`TIGHT_CHI2`), p(π_s) > 0.3 GeV (`DS_TIGHT_PS`), `cosPoint` > 0.95 (`DS_TIGHT_COSPOINT`), and a **primary-pattern veto**: a candidate whose K and π are both in the fitted primary set while its slow pion is not is rejected, that pattern being combinatorial rather than a real D\* topology.

**Stored per-daughter block** (`trkK_`, `trkPi_`, and `trkPis_` for the slow pion). Original-track index, **charge**, momentum, cosθ, d0/z0/σ(d0), nVDET/nITC hits, track χ²/ndf, whether the track was in the fitted primary set, and the pads/wires dE/dx value and error. On MC: `dstar_class` (1 = all three tracks link to the K, π and π_s of one true D\* with the correct K/π assignment, 2 = the same D\* with K↔π swapped, 3 = a true D⁰ from a D\* but the wrong slow pion, 4 = a true D⁰→Kπ whose mother is not a D\*, 0 = otherwise), `d0_class` (1 = correct assignment, 2 = swapped), the `*_trueidx` back-pointers, the daughter MC PDG and its LUND mother PDG, and the `truedstar_*` / `trued0_*` lists of the generated decays with their origin class (5 = from a b hadron, 4 = from c, 0 = other), the D⁰ flight distance, how many daughters were tracked, the per-daughter track pool (`*_pool`: −1 unlinked, 0 primary set, 1 secondary set, 2 neither) and the `found_loose`/`found_tight` efficiency flags. Truth provenance is decoded from `generatorStatus` = 10000·KS + LUND mother line, since the MCParticles daughter relations are empty in these files. Track→MC matching uses every link of a track (the relation carries several unweighted partners for ~10% of linked tracks); for a candidate of class ≥1 the stored daughter PDG is that of the matched link, otherwise the first link. The φ→KK truth follows the same convention.

### Index maps, vertex flags and other stored branches

`prim2origIdx` and `sec2origIdx` map the primary and secondary track collections back to the original `Tracks` index space, which is the frame every finder's `*_origIdx` branch uses; `recopart_tracks_index`, `recopart_tracks_begin` and `recopart_tracks_end` expose the `ReconstructedParticle`→`Track` relation itself, so the same join can be redone offline. The primary-vertex flags are `pv_converged` (the position fit converged), `pv_split_converged` (every pruning pass converged, not just the final fit), `pv_trivial` (too few tracks entered the fit, so the vertex carries no event information even when both converged) and `pv_good`, the single named predicate `pv_converged && pv_split_converged && !pv_trivial` that all downstream users should test. `chi2_o_ndf_tracks_all` holds χ²/ndf; it previously held ndf. The legacy `sv_*` block is empty for an event whose PV fit did not converge.

### Joining the candidate lists to the tracks and to the pfcand block

Three store-only branch groups make the reconstructed objects addressable in the **original `Tracks`** index space, which is what every finder's `*_origIdx` branch already uses. `pfcand_trackIdx` gives that index for each jet constituent (−1 for a constituent with no track), so the `pfcand_*` block and any candidate leg can be matched directly instead of through the `RecoParticle`→`Track` relation. Every candidate leg additionally carries `<leg>_isChargedHad`, a tri-state particle-flow label: 1 = the leg's linked reconstructed particle is a PF charged hadron, 0 = it is an electron, muon or another type, −1 = the track has no linked reconstructed particle (mostly very soft tracks, below the PF momentum reach). The label uses the same particle-flow type codes as `pfcand_isChargedHad`.

`trk_member` and `trk_nCand` are per-track arrays over the whole `Tracks` collection. `trk_member` is a bitmask recording every set a track belongs to: bit 0 fitted primary-vertex set, bit 1 daughter of any stored V0 candidate, bit 2 daughter of a *tight* V0, bit 3 leg of any stored φ→KK candidate, bit 4 leg of a φ candidate passing `phikk_wp`, bit 5 leg of any stored D⁰ candidate, bit 6 leg of any stored D\* candidate (slow pion included), bit 7 leg of a D\* passing `dstar_tight`, bit 8 track of a secondary vertex (`svn`), bit 9 baseline-selected track. `trk_nCand` counts how many stored candidates (V0 + φ + D⁰ + D\*) use the track, which is the multiplicity the offline 1/n de-duplication weight needs — the finders claim no tracks exclusively, so one track can serve many candidates. Both are built in a single pass over the finished candidate lists, so they add no reconstruction work.

### STAGE 2:

Don't touch stage2.py!
Open up stage2_all.py and change the desired input and output directories. 
Set the number of cpus.
Now you can decide if you want to divide each flavor into multiple files, then you can change the argument `n_final_files`.
Run it with nightlies.






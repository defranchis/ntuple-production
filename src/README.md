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
`/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/<year>/<mc_type>/stage1/<version_tag>/<flavour_name>.root`. 

Fraction of events to process can be set via `--fraction <val>`, default is to process all events. 

`--year <year>` and `--MCtype <type>` are also supported command line arguments, currently we only have `1994` and `zqq` here. 

### Run on data:
```bash
fccanalysis run stage1.py -- --tag <version_tag> --doData 
```

Output files will be in: `/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/<year>/stage1/<version_tag>/`

`--year` and `--fraction` is also supported as an argument here. 

`--noDedxGate` (data and MC) accepts every linked dE/dx measurement as valid, i.e. switches off the failed-leg omega sentinel gate. Use it with converters that no longer copy the track omega into a failed leg.

`pfcand_d0` and `pfcand_z0` are the **raw** perigee impact parameters of the track linked to each jet constituent, taken from the stored track state: referenced to the coordinate origin, in cm, in the LCIO sign convention (the ALEPH→LCIO flip of `D0` and `omega` is applied; the `pfcand_*` covariance branches are read from the same flipped collection). They are *not* recomputed at the primary vertex — `pfcand_dxy` and `pfcand_dz` are the PV-referenced ones. A constituent with no track (a neutral) carries the guard value −9 in both, again as for the covariance branches; the PV-referenced `pfcand_dxy/dz/phi0/C/ct` also read −9 for a charge-0 constituent that carries a track link.

The −9 sentinel of `pfcand_d0`/`pfcand_z0` lies inside the physical range of those variables, so the safe test for "no track" is the track link itself (or exact equality with −9), not a range cut.

`pfcand_C` is now half the fitted curvature, |omega|/2 in 1/cm, signed by the charge. It was previously Bz·c/(2 pT) from the energy-flow pT in 1/mm, so both the scale (a factor 10) and the source differ — the latter matters for V0 daughters, whose energy-flow momentum is quoted at the V0 vertex.

`pfcand_dEdx_pads_type` and `pfcand_dEdx_wires_type` are no longer a validity mask: an accepted leg can carry any type and a rejected one reads −9 in value, error and type, so test the value branch rather than `type == 0`.

### Run on batch:
```
fccanalysis submit stage1.py -- --tag VXX-XX --MCflavour X --batch --chunks X
```

### STAGE 2:

Don't touch stage2.py!
Open up stage2_all.py and change the desired input and output directories. 
Set the number of cpus.
Now you can decide if you want to divide each flavor into multiple files, then you can change the argument `n_final_files`.
Run it with nightlies.






## Charged PF candidate dump

`stage1_chargedpf.py` writes one flat vector per quantity over the **charged particle-flow
candidates** of the event, i.e. the `RecoParticles` with `charge != 0` that own a track through
the `RecoParticles -> Tracks` relation. There is no jet clustering, no vertexing and no
SV/V0 reconstruction; all branches carry the prefix `cpf_` and are ordered identically, so
entry `i` of every branch belongs to the same candidate.

```bash
fccanalysis run stage1_chargedpf.py -- --tag <version_tag> --MCflavour <flavour_index>
fccanalysis run stage1_chargedpf.py -- --tag <version_tag> --doData
```

Defaults follow `stage1.py`, with `stage1_chargedpf` in place of `stage1` in the output path.
`--inputDir` and `--outdir` override the input and output directories, `--fraction` the
fraction of input to process, and `--nthreads` the number of threads (default 16).
`--year` and `--MCtype` behave as for `stage1.py`.

The PID flags, the dE/dx values and the track quantities are produced by the same helpers
`stage1.py` uses for its `pfcand_*` branches, so the two agree candidate by candidate.

Units are as stored in the input, nothing is rescaled: lengths in cm, momenta and masses
in GeV, angles in rad, magnetic field in T (see `aleph_units.h`).

### Event level

| branch | meaning |
|---|---|
| `run_number`, `event_number` | run and event number from the `EventHeader` |
| `n_cpf` | number of charged PF candidates in the event |
| `jetPID` | truth flavour of the event as defined in `stage1.py` (`-999` on data) |

### Kinematics and identification

| branch | meaning / unit |
|---|---|
| `cpf_px`, `cpf_py`, `cpf_pz`, `cpf_e` | energy-flow four-momentum [GeV] |
| `cpf_p`, `cpf_pt` | momentum and transverse momentum [GeV] |
| `cpf_theta`, `cpf_phi` | polar and azimuthal angle of the momentum [rad], `phi` in (-pi, pi] |
| `cpf_charge` | electric charge [e] |
| `cpf_mass` | mass assigned by the energy-flow algorithm [GeV] |
| `cpf_type` | `ReconstructedParticle::type` as stored (0 in the 1994 files) |
| `cpf_pidType` | raw `ParticleID::type` of the candidate |
| `cpf_isMu`, `cpf_isEl`, `cpf_isChargedHad` | 1/0 flags for `ParticleID::type` 2, 1 and 0 |

### dE/dx and PID p-values

| branch | meaning / unit |
|---|---|
| `cpf_dEdx_pads_value`, `cpf_dEdx_pads_error` | TPC pads dE/dx and its error [normalised dE/dx, ~1 for a MIP] |
| `cpf_dEdx_pads_type` | `dQdx.type` word of the pads measurement |
| `cpf_dEdx_wires_*` | the same three for the TPC wires |
| `cpf_PID_pval_pads_{ele,mu,pi,kaon,proton}` | signed p-value of the pads dE/dx under each mass hypothesis, from the Bethe-Bloch fits of `analyzer.h` |
| `cpf_PID_pval_wires_{ele,mu,pi,kaon,proton}` | the same for the wires |

A rejected or missing dE/dx measurement reads `-9` in value, error, type and in all five
p-values, so test the value branch rather than the type.

### Track parameters

All taken from the first track state of the candidate's own track, after the ALEPH -> LCIO
sign flip of `D0` and `omega` that `stage1.py` also applies.

| branch | meaning / unit |
|---|---|
| `cpf_d0`, `cpf_z0` | perigee impact parameters w.r.t. the coordinate origin [cm] |
| `cpf_phi0` | azimuth of the momentum at the perigee [rad, (0, 2pi) in the 1994 files] |
| `cpf_omega` | signed curvature [1/cm]; pT [GeV] = `kPtPerTeslaCm` (`aleph_units.h`) * Bz [T] / abs(omega) |
| `cpf_tanLambda` | tangent of the dip angle [dimensionless] |
| `cpf_d0d0`, `cpf_z0z0` | covariance diagonal of `d0` and `z0` [cm^2] |
| `cpf_phi0phi0`, `cpf_tanLtanL` | covariance diagonal of `phi0` [rad^2] and of `tanLambda` |
| `cpf_omegaomega` | covariance diagonal of `omega` [1/cm^2] |
| `cpf_trackChi2`, `cpf_trackNdof` | fit chi2 and number of degrees of freedom |
| `cpf_nTrackHits_VDET/ITC/TPC` | number of hits in the vertex detector, inner tracking chamber and time projection chamber |
| `cpf_trackIdx` | index of the track in the `Tracks` collection |

### MC truth

Resolved through the `trackMCLink` association (`_trackMCLink_from` = `Tracks` index,
`_trackMCLink_to` = `MCParticles` index), taking the first link of the candidate's track.
On `--doData` the branches exist with the same names and types and read `-999`
(`cpf_mc_n_links` reads 0).

| branch | meaning / unit |
|---|---|
| `cpf_mc_pdg` | PDG code of the linked MC particle |
| `cpf_mc_parent_pdg` | PDG code of its first parent |
| `cpf_mc_p` | momentum of the linked MC particle [GeV] |
| `cpf_mc_vertex_x/y/z` | its production vertex [cm] |
| `cpf_mc_n_links` | number of `trackMCLink` entries of that track |

The 1994 EDM4HEP files carry no `MCParticles` parent or daughter relations, so
`cpf_mc_parent_pdg` is `-999` for every candidate until the converter fills them.

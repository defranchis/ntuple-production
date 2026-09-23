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

`pfcand_mc_pdg`, `pfcand_mc_p` [GeV], `pfcand_mc_vtx_r` and `pfcand_mc_vtx_z` [cm, transverse radius and z of the production vertex] describe the MC particle of the first `trackMCLink` entry of the constituent's track; `pfcand_mc_nlinks` counts that track's links. A track links to every MC particle with a hit on it, so the first link is not always the only contributor. All five read −999 for neutral constituents, for tracks without a link, and on data.

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






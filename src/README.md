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

The runs listed in `run_veto.py` (vertexing-unusable, tracking-quality and off-plateau-field runs) are vetoed before any selection; `--excludeRuns RUN [RUN ...]` adds runs, `--noRunVeto` switches the list off. `eventsProcessed` still counts the raw input.

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






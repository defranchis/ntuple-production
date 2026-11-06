#!/bin/bash
# reproduce.sh <nproc> <sample> <tag>: run the frozen final at both beam-spot settings (final_LR, final_S1) and nominal
# stage1 on <sample> (e.g. DEV) under the new run tag <tag>, then evaluate both finals against nominal on the same
# events (100 paired bootstrap resamples) -> $JV_RUNS/<tag>/rows_<sample>.json. Needs JV_INPUTS and JV_RUNS (README.md).
H=$(dirname "$(readlink -f "$0")")
N=${1:?usage: reproduce.sh <nproc> <sample> <tag>}; S=${2:?usage: reproduce.sh <nproc> <sample> <tag>}; T=${3:?usage: reproduce.sh <nproc> <sample> <tag>}
: "${JV_INPUTS:?set JV_INPUTS}" "${JV_RUNS:?set JV_RUNS}"
if [ -e "$JV_RUNS/$T/$S" ]; then echo "$JV_RUNS/$T/$S exists: choose a new tag"; exit 2; fi
bash "$H/run_py.sh" 7200 jv.run "$N" "$S" "$T" final_LR=configs/final_LR.json final_S1=configs/final_S1.json nominal || exit 1
bash "$H/run_py.sh" 3600 jv.evaluate "$N" "$T" "$S" "$JV_RUNS/$T/rows_$S.json" --ref nominal --nboot 100

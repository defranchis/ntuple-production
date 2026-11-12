#!/bin/bash
# reproduce.sh <nproc> <sample> <tag> [v0]: run the frozen final at both beam-spot settings (final_LR, final_S1) and nominal
# stage1 on <sample> (e.g. DEV) under the new run tag <tag>, then evaluate both finals against nominal on the same
# events (100 paired bootstrap resamples) -> $JV_RUNS/<tag>/rows_<sample>.json. With 'v0' also the nominal V0
# configuration (v0_LR, v0_S1) in the same run and evaluation, and the V0 list of v0_S1 against stage1's tight V0 list
# -> $JV_RUNS/<tag>/v0_<sample>.json, with one row per S1 working point (thresholds quoted from the working-point table
# in README section "V0 class"): read K0S from row run:<tag>:v0_S1:0.398, Lambda from :0.329, conversions from :0.381
# (each row applies its threshold to every species). 'v0' needs inputs harvested with the true-V0 arrays (t_v0t). Needs
# JV_INPUTS and JV_RUNS (README.md); refuses an existing $JV_RUNS/<tag>/<sample>. Exit status: 2 = bad argument or
# existing tag, 1 = failed step.
H=$(dirname "$(readlink -f "$0")")
USAGE="usage: reproduce.sh <nproc> <sample> <tag> [v0]"
N=${1:?$USAGE}; S=${2:?$USAGE}; T=${3:?$USAGE}
case "$4" in ""|v0) ;; *) echo "$USAGE (unknown 4th argument '$4')"; exit 2;; esac
: "${JV_INPUTS:?set JV_INPUTS}" "${JV_RUNS:?set JV_RUNS}"
if [ -e "$JV_RUNS/$T/$S" ]; then echo "$JV_RUNS/$T/$S exists: choose a new tag"; exit 2; fi
V0ROWS=""
[ "$4" = v0 ] && V0ROWS="v0_LR=configs/v0_LR.json v0_S1=configs/v0_S1.json"
bash "$H/run_py.sh" 7200 jv.run "$N" "$S" "$T" final_LR=configs/final_LR.json final_S1=configs/final_S1.json $V0ROWS nominal || exit 1
bash "$H/run_py.sh" 3600 jv.evaluate "$N" "$T" "$S" "$JV_RUNS/$T/rows_$S.json" --ref nominal --nboot 100 || exit 1
if [ "$4" = v0 ]; then
    WP="0.398 0.329 0.381"          # S1 working points K0S, Lambda, conversions: working-point table in README section "V0 class"
    R="stage1_tight"; for w in $WP; do R="$R,run:$T:v0_S1:$w"; done
    bash "$H/run_py.sh" 3600 jv.evaluate_v0 "$N" "$S" "$JV_RUNS/$T/v0_$S.json" --rows "$R" --ref stage1_tight || exit 1
fi

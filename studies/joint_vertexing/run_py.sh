#!/bin/bash
# run_py.sh <timeout s> <python -m module + args...>: LCG view, single-threaded BLAS, numba cache on local disk,
# cwd = this folder (package jv; relative paths in the arguments are relative to it), one status line at the end
# (DONE / TIMEOUT / FAILED), exit code of the module.
# Environment: JV_LCG_SETUP (setup script of an LCG view with numpy, numba, uproot, awkward, matplotlib; default LCG_107),
# NUMBA_CACHE_DIR (default /tmp/jv_numba_$USER), and the JV_* directories the module needs (README.md).
source ${JV_LCG_SETUP:-/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh} > /dev/null 2>&1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-/tmp/jv_numba_$USER}
cd "$(dirname "$(readlink -f "$0")")" || exit 1
T=$1; shift
timeout $T python3 -m "$@"
rc=$?
if [ $rc = 0 ]; then echo "DONE $* $(date)"; elif [ $rc = 124 ]; then echo "TIMEOUT $* $(date)"; else echo "FAILED rc=$rc $* $(date)"; fi
exit $rc

#!/bin/bash
# Wrapper script for HTCondor.
# Usage: run_cc_proton.sh <replica_id>
# Called by cc_proton_xsec.sub via $(Process).

REPLICA=$1
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON=/nfs/pic.es/user/f/fmontagn/miniconda3/envs/nuDIS/bin/python

cd "$REPO_ROOT"

$PYTHON scripts/run_cc_proton_xsec.py \
    --replica     "$REPLICA"  \
    --x_min_bin   1e-3        \
    --x_max_bin   1           \
    --num_bins_x  10          \
    --Q2_min_bin  10          \
    --Q2_max_bin  1e4         \
    --num_bins_Q2 10          \
    --max_workers 4

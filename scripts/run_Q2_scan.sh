#!/bin/bash

# -----------------------------------------------------------------------------
# Run 4 processes x 3 Q2 bins = 12 jobs in parallel.
#
# Parallelism settings:
#   MAX_PARALLEL    = max concurrent jobs
#   NCORES_PER_JOB  = CPU cores passed to each VEGAS run
#
# Constraint: MAX_PARALLEL * NCORES_PER_JOB <= total cores on your machine.
# Default: 4 parallel jobs x 2 cores = 8 cores total.
#
# Logs: data/logs/<folder>.log  (one file per job)
# -----------------------------------------------------------------------------

MAX_PARALLEL=4
NCORES_PER_JOB=2

CARD="card/events_card.dat"

PROCESSES=("d_p" "ubar_p" "s_p" "cbar_p")

Q2_BINS=(
    "8 100"
    "100 1000"
    "1000 10000"
)

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nudis

mkdir -p data/logs

pids=()
labels=()

for proc in "${PROCESSES[@]}"; do
    for bin in "${Q2_BINS[@]}"; do
        q2min=$(echo $bin | awk '{print $1}')
        q2max=$(echo $bin | awk '{print $2}')
        folder="${proc}/Q2_${q2min}_${q2max}"
        logfile="data/logs/${proc}_Q2_${q2min}_${q2max}.log"

        echo "Launching: $folder  (log: $logfile)"

        python3 events.py \
            --card        "$CARD"          \
            --process     "$proc"          \
            --folder      "$folder"        \
            --Q2_min_bin  "$q2min"         \
            --Q2_max_bin  "$q2max"         \
            --ncores      "$NCORES_PER_JOB" \
            > "$logfile" 2>&1 &

        pids+=($!)
        labels+=("$folder")

        # When MAX_PARALLEL jobs are running, wait for all of them before launching more
        if (( ${#pids[@]} >= MAX_PARALLEL )); then
            for i in "${!pids[@]}"; do
                wait "${pids[$i]}" && echo "Done:   ${labels[$i]}" \
                                   || echo "FAILED: ${labels[$i]}"
            done
            pids=()
            labels=()
        fi
    done
done

# Wait for any remaining jobs
for i in "${!pids[@]}"; do
    wait "${pids[$i]}" && echo "Done:   ${labels[$i]}" \
                        || echo "FAILED: ${labels[$i]}"
done

echo "All 12 runs completed."

#!/bin/bash

# Master job submission script for configuration: T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600

set -e

echo "=== Submitting jobs for detector configuration: T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600 ==="

# Submit build jobs in a sequential chain (to avoid source conflicts)
echo "Submitting 100 build jobs in sequential chain..."
BUILD_JOB_IDS=()
PREV_BUILD_ID=""

for run in $(seq 1 100); do
    if [ -z "$PREV_BUILD_ID" ]; then
        # First build has no dependency
        BUILD_ID=$(sbatch /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_${run}_T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600.job | awk '{print $4}')
        echo "Submitted build $run: Job ID $BUILD_ID (no dependency)"
    else
        # Subsequent builds depend on previous build completion
        BUILD_ID=$(sbatch --dependency=afterok:$PREV_BUILD_ID /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_${run}_T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600.job | awk '{print $4}')
        [ $(($run % 10)) -eq 0 ] && echo "Submitted builds $(($run-9))-$run (last: $BUILD_ID)"
    fi
    BUILD_JOB_IDS+=($BUILD_ID)
    PREV_BUILD_ID=$BUILD_ID
    
    # Submit corresponding run job (depends only on its own build completing)
    RUN_ID=$(sbatch --dependency=afterok:$BUILD_ID /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/run_${run}_T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600.job | awk '{print $4}')
    
    # Store run job ID for later dependency
    if [ $run -eq 1 ]; then
        RUN_JOB_IDS=($RUN_ID)
    else
        RUN_JOB_IDS+=($RUN_ID)
    fi
done

echo "All build jobs submitted in sequential chain"
echo "All run jobs submitted in parallel (each depends on its own build)"

# Create dependency string for parsing jobs (wait for all run jobs)
RUN_DEPENDENCY=$(IFS=:; echo "${RUN_JOB_IDS[*]}")

# Submit parsing jobs with dependency on all run jobs
echo "Submitting parsing jobs (dependent on all run jobs)..."
PARSE_JOB_IDS=()
for group in $(seq 1 10); do
    JOB_ID=$(sbatch --dependency=afterok:$RUN_DEPENDENCY /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/parse_group_${group}_T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600.job | awk '{print $4}')
    PARSE_JOB_IDS+=($JOB_ID)
    echo "Submitted parsing group $group: Job ID $JOB_ID"
done

# Create dependency string for combination job (wait for all parsing jobs)
PARSE_DEPENDENCY=$(IFS=:; echo "${PARSE_JOB_IDS[*]}")

# Submit combination job with dependency on parsing jobs
echo "Submitting combination job (dependent on parsing jobs)..."
COMBINE_JOB_ID=$(sbatch --dependency=afterok:$PARSE_DEPENDENCY /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/combine_T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600.job | awk '{print $4}')
echo "Submitted combination job: Job ID $COMBINE_JOB_ID"

echo
echo "=== Job submission completed! ==="
echo "Build jobs: ${BUILD_JOB_IDS[@]}"
echo "Run jobs: ${RUN_JOB_IDS[@]}"
echo "Parsing jobs: ${PARSE_JOB_IDS[@]}"
echo "Combination job: $COMBINE_JOB_ID"
echo
echo "NOTE: Build jobs run sequentially (~10 min each) to avoid source conflicts"
echo "      Run jobs start in parallel as soon as their build completes (~2 days each for 100k events)"
echo "      Total wall time: ~17 minutes (builds) + ~2 days (parallel runs) for 10M total events"
echo
echo "Monitor progress with: squeue -u $USER"
echo "Results will be in: /users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-550_T4-560_P1-215_P2-230_F1-260_F2-270_E1-600"

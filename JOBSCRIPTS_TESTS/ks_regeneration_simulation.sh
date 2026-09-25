#!/bin/bash

# KS Regeneration Simulation Script
#
# Simulates KL→KS regeneration events by replacing the KL beam with a KS beam
# at a specified z-position between the dipole magnet exit and the flux monitor.
#
# The detector geometry is fixed to the dissertation final configuration:
#   T1=425, T2=435, T3=465, T4=475, P1=400, P2=415, F1=455, F2=465, E1=700
#
# Only the KS beam start position (KS_Z0) is varied, scanning the region
# from the dipole exit (z=200 cm) to the start of the detector array (z=400 cm).
#
# Usage:  ./ks_regeneration_simulation.sh <KS_Z0_cm>
# Where:  KS_Z0_cm is the beam start z-position in cm (must be in range 200–400)
#
# Example: ./ks_regeneration_simulation.sh 300
#   → Simulates KS beam starting at z=300 cm, detector fixed at dissertation config.

set -e

# ── Fixed detector configuration (dissertation final output) ─────────────────
T1=425; T2=435; T3=465; T4=475
P1=400; P2=415
F1=455; F2=465
E1=700
# ────────────────────────────────────────────────────────────────────────────

DIPOLE_END=200
DETECTOR_START=400  # P1 position — KS must decay before reaching this

# ── Argument parsing ─────────────────────────────────────────────────────────
if [ $# -ne 1 ]; then
    echo "Usage: $0 <KS_Z0_cm>"
    echo "  KS_Z0_cm : KS beam start z-position in cm (integer, 200–400)"
    echo ""
    echo "Example: $0 300"
    echo ""
    echo "Fixed detector config: T1=${T1} T2=${T2} T3=${T3} T4=${T4}"
    echo "                       P1=${P1} P2=${P2} F1=${F1} F2=${F2} E1=${E1}"
    exit 1
fi

KS_Z0=$1

# Validate KS_Z0 is numeric and in range
if ! [[ "$KS_Z0" =~ ^[0-9]+$ ]]; then
    echo "ERROR: KS_Z0_cm must be a non-negative integer (got: $KS_Z0)"
    exit 1
fi

if (( KS_Z0 < DIPOLE_END || KS_Z0 > DETECTOR_START )); then
    echo "ERROR: KS_Z0=${KS_Z0} cm is outside the valid range [${DIPOLE_END}, ${DETECTOR_START}] cm"
    echo "       The KS beam must start after the dipole (z>${DIPOLE_END} cm)"
    echo "       and before the first detector (z<=${DETECTOR_START} cm)"
    exit 1
fi

# ── Configuration strings ────────────────────────────────────────────────────
DETECTOR_CONFIG="T1-${T1}_T2-${T2}_T3-${T3}_T4-${T4}_P1-${P1}_P2-${P2}_F1-${F1}_F2-${F2}_E1-${E1}"
CONFIG_STR="KS_Z0-${KS_Z0}_${DETECTOR_CONFIG}"

echo "=== KS Regeneration Simulation Pipeline ==="
echo "KS beam start position : z = ${KS_Z0} cm"
echo "Detector configuration : ${DETECTOR_CONFIG}"
echo "Full config string     : ${CONFIG_STR}"
echo ""

# ── Path definitions ─────────────────────────────────────────────────────────
BASE_SCENARIO_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/BASE_SCENARIO_DIR"
SCENARIO_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/${CONFIG_STR}"
RESULTS_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/${CONFIG_STR}"
PARSING_DIR="/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/${CONFIG_STR}_DATA_READ"
LOG_BASE_DIR="/users/bp969/scratch/JOB_LOGS/${CONFIG_STR}"

# ── Create directory structure ───────────────────────────────────────────────
echo "=== Creating directory structure ==="
mkdir -p "${LOG_BASE_DIR}/build"
mkdir -p "${LOG_BASE_DIR}/simulation"
mkdir -p "${LOG_BASE_DIR}/parsing"
mkdir -p "${LOG_BASE_DIR}/combination"
mkdir -p "${RESULTS_DIR}"
mkdir -p "${PARSING_DIR}"
mkdir -p "/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs"
echo "Created log directories in ${LOG_BASE_DIR}"

# Copy base simulation directory if not already present
if [ ! -d "${SCENARIO_DIR}/src" ]; then
    echo "Copying simulation files to: ${SCENARIO_DIR}"
    mkdir -p "${SCENARIO_DIR}"
    cp -r "${BASE_SCENARIO_DIR}"/* "${SCENARIO_DIR}/"
    echo "Simulation files copied."
else
    echo "Configuration directory already exists: ${SCENARIO_DIR}"
fi

# Copy ROOT parsing scripts to parsing directory
echo "Copying ROOT parsing scripts to ${PARSING_DIR}"
cp "/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C"  "${PARSING_DIR}/"
cp "/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C"              "${PARSING_DIR}/"
cp "/users/bp969/scratch/VIKING_FOLDER/HISTOGRAM_MAKING/parse_diag_events.py"          "${PARSING_DIR}/"

# ── Step 1: Edit detector positions ─────────────────────────────────────────
echo ""
echo "=== Step 1: Editing detector construction file ==="

DETECTOR_FILE="${SCENARIO_DIR}/src/JLabKDetectorConstruction.cc"
BACKUP_FILE="${DETECTOR_FILE}.backup"

if [ ! -f "$BACKUP_FILE" ]; then
    cp "$DETECTOR_FILE" "$BACKUP_FILE"
    echo "Created backup: $BACKUP_FILE"
fi

echo "Updating detector positions in $DETECTOR_FILE"
sed -i "s|G4double trkrPosZ\[4\] = {[^}]*};|G4double trkrPosZ[4] = {${T1}.*cm, ${T2}.*cm, ${T3}.*cm, ${T4}.*cm};|" "$DETECTOR_FILE"
sed -i "s|G4double pizzaPosZ\[4\] = {[^}]*};|G4double pizzaPosZ[4] = {${P1}.*cm, ${P2}.*cm, 1100.*cm, 1200.*cm};|" "$DETECTOR_FILE"
sed -i "s|G4double friPosZ\[2\] = {[^}]*};|G4double friPosZ[2] = {${F1}.*cm, ${F2}.*cm};|" "$DETECTOR_FILE"
sed -i "s|G4double tofPosZ = [^;]*;|G4double tofPosZ = ${E1}.*cm;|" "$DETECTOR_FILE"
echo "Detector positions updated."

# ── Step 2: Edit primary generator to use KS at z=KS_Z0 ─────────────────────
echo ""
echo "=== Step 2: Editing primary generator action for KS beam ==="

PRIMGEN_FILE="${SCENARIO_DIR}/src/JLabKPrimaryGeneratorAction.cc"
PRIMGEN_BACKUP="${PRIMGEN_FILE}.backup"

if [ ! -f "$PRIMGEN_BACKUP" ]; then
    cp "$PRIMGEN_FILE" "$PRIMGEN_BACKUP"
    echo "Created backup: $PRIMGEN_BACKUP"
fi

# Switch particle from kaon0L to kaon0S
sed -i 's|FindParticle("kaon0L")|FindParticle("kaon0S")|g' "$PRIMGEN_FILE"

# Set the beam start z-position.  The original file has: G4double z0 = 0.0;
# We replace it with the requested KS start position (in cm Geant4 units).
sed -i "s|G4double z0 = [^;]*;|G4double z0 = ${KS_Z0} \* cm;|" "$PRIMGEN_FILE"

echo "Primary generator updated:"
echo "  Particle : kaon0S"
echo "  z0       : ${KS_Z0} cm"
grep 'FindParticle'     "$PRIMGEN_FILE" | head -1
grep 'G4double z0'      "$PRIMGEN_FILE" | head -1

# ── Step 3: Create batch macro ───────────────────────────────────────────────
echo ""
echo "=== Step 3: Creating batch macro ==="

cat > "${SCENARIO_DIR}/batch.mac" << 'EOF'
# Batch mode macro for headless execution
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

/run/numberOfThreads 4

/run/initialize

/run/beamOn 100000
EOF

echo "Created batch.mac"

# ── Step 4: Create build and run jobscripts ──────────────────────────────────
echo ""
echo "=== Step 4: Creating build and run jobscripts ==="

create_build_jobscript() {
    local run_num=$1
    local jobscript="/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_${run_num}_${CONFIG_STR}.job"

    cat > "$jobscript" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=KSBuild_${CONFIG_STR:0:30}_${run_num}
#SBATCH --partition=nodes
#SBATCH --time=0-04:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${LOG_BASE_DIR}/build/Build_${CONFIG_STR}_${run_num}_%j.log
#SBATCH --error=${LOG_BASE_DIR}/build/Build_${CONFIG_STR}_${run_num}_%j.err

set -e

module purge
module load Geant4/11.2.2-GCC-12.3.0
module load ROOT/6.30.06-foss-2023a

OUTFILE="${CONFIG_STR}_${run_num}.root"
RUNACTION_FILE="${SCENARIO_DIR}/src/JLabKRunAction.cc"
EXECUTABLE_DEST="${SCENARIO_DIR}/jlabk_${run_num}"

echo "=== KS Build ${run_num} starting at \$(date) ==="
echo "Output filename: \$OUTFILE"

cd ${SCENARIO_DIR}

if [ ${run_num} -eq 1 ]; then
    echo "First build: Cleaning CMake cache and reconfiguring..."
    rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake
    cmake .
fi

echo "Setting output filename to \$OUTFILE in source code..."
sed -i 's|G4String fileName = "[^"]*";|G4String fileName = "'\$OUTFILE'";|' "\$RUNACTION_FILE"

echo "Rebuilding simulation executable..."
make -j4

cp jlabk "\$EXECUTABLE_DEST"
echo "Created executable: \$EXECUTABLE_DEST"

echo "=== KS Build ${run_num} completed at \$(date) ==="
EOF

    echo "Created build jobscript ${run_num}: $jobscript"
}

create_run_jobscript() {
    local run_num=$1
    local jobscript="/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/run_${run_num}_${CONFIG_STR}.job"

    cat > "$jobscript" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=KSRun_${CONFIG_STR:0:30}_${run_num}
#SBATCH --partition=nodes
#SBATCH --time=0-04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${LOG_BASE_DIR}/simulation/Run_${CONFIG_STR}_${run_num}_%j.log
#SBATCH --error=${LOG_BASE_DIR}/simulation/Run_${CONFIG_STR}_${run_num}_%j.err

set -e

module purge
module load Geant4/11.2.2-GCC-12.3.0
module load ROOT/6.30.06-foss-2023a

SEED=\$((1000000 + ${run_num}))
OUTFILE="${CONFIG_STR}_${run_num}.root"
EXECUTABLE="${SCENARIO_DIR}/jlabk_${run_num}"

echo "=== KS Simulation ${run_num} starting at \$(date) ==="
echo "KS beam start: z = ${KS_Z0} cm"
echo "Seed: \$SEED  |  Output: \$OUTFILE"

cd ${SCENARIO_DIR}

if [ ! -f "\$EXECUTABLE" ]; then
    echo "ERROR: Executable \$EXECUTABLE not found!"
    exit 1
fi

\$EXECUTABLE batch.mac \$SEED
EXIT_CODE=\$?

if [ -f "\$OUTFILE" ]; then
    FILESIZE=\$(stat -c%s "\$OUTFILE" 2>/dev/null)
    if [ "\$FILESIZE" -gt 1000000 ]; then
        mv "\$OUTFILE" "${RESULTS_DIR}/\$OUTFILE"
        echo "Completed simulation ${run_num}: \$OUTFILE (\${FILESIZE} bytes)"
        [ \$EXIT_CODE -eq 139 ] && echo "NOTE: Exit code 139 (cleanup segfault) ignored — output OK"
        rm -f "\$EXECUTABLE"
    else
        echo "ERROR: Output file too small (\${FILESIZE} bytes) — simulation likely failed"
        rm -f "\$EXECUTABLE"
        exit 1
    fi
else
    echo "ERROR: Output file '\$OUTFILE' not created (exit code: \$EXIT_CODE)"
    rm -f "\$EXECUTABLE"
    exit 1
fi

echo "=== KS Simulation ${run_num} completed at \$(date) ==="
EOF

    echo "Created run jobscript ${run_num}: $jobscript"
}

for run in $(seq 1 100); do
    create_build_jobscript $run
    create_run_jobscript $run
done

# ── Step 5: Create parsing jobscripts ────────────────────────────────────────
echo ""
echo "=== Step 5: Creating parsing jobscripts ==="

create_parsing_jobscript() {
    local file_num=$1
    local jobscript="/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/parse_file_${file_num}_${CONFIG_STR}.job"

    cat > "$jobscript" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=KSParse_${CONFIG_STR:0:28}_F${file_num}
#SBATCH --partition=nodes
#SBATCH --time=0-04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${LOG_BASE_DIR}/parsing/Parse_${CONFIG_STR}_F${file_num}_%j.log
#SBATCH --error=${LOG_BASE_DIR}/parsing/Parse_${CONFIG_STR}_F${file_num}_%j.err

set -e

module purge
module load ROOT/6.30.06-foss-2023a

cd ${PARSING_DIR}

INPUT_FILE="${RESULTS_DIR}/${CONFIG_STR}_${file_num}.root"
echo "Waiting for: \$INPUT_FILE"
while [ ! -f "\$INPUT_FILE" ]; do
    sleep 30
done

ACCEPTANCE_OUTPUT="${RESULTS_DIR}/${CONFIG_STR}_${file_num}_acceptance.root"
VECTORS_OUTPUT="${RESULTS_DIR}/${CONFIG_STR}_${file_num}_vectors.root"

INPUT_FILENAME="${CONFIG_STR}_${file_num}.root"
cp "\$INPUT_FILE" "\$INPUT_FILENAME"

root -l -b -q "KLong_save_momentum_acceptance.C(\"\$INPUT_FILENAME\")"

DIAG_LOG="${RESULTS_DIR}/${CONFIG_STR}_${file_num}_diag_clean.txt"
root -l -b -q "KLong_save_vectors.C(\"\$INPUT_FILENAME\")" 2>&1 | \
    python3 "${PARSING_DIR}/parse_diag_events.py" --log "\$DIAG_LOG"
ROOT_EXIT=\${PIPESTATUS[0]}
if [ \$ROOT_EXIT -ne 0 ]; then
    echo "ERROR: KLong_save_vectors.C failed (exit code \$ROOT_EXIT)"
    exit \$ROOT_EXIT
fi

mv "${CONFIG_STR}_${file_num}_acceptance.root" "\$ACCEPTANCE_OUTPUT" 2>/dev/null || true
mv "${CONFIG_STR}_${file_num}_vectors.root"    "\$VECTORS_OUTPUT"    2>/dev/null || true
rm -f "\$INPUT_FILENAME"

[ -f "\$ACCEPTANCE_OUTPUT" ] || { echo "ERROR: \$ACCEPTANCE_OUTPUT not created!"; exit 1; }
[ -f "\$VECTORS_OUTPUT"    ] || { echo "ERROR: \$VECTORS_OUTPUT not created!";    exit 1; }

echo "Completed parsing for file ${file_num}"
EOF

    echo "Created parsing jobscript: $jobscript"
}

for file_num in $(seq 1 100); do
    create_parsing_jobscript $file_num
done

# ── Step 6: Create combination jobscript ─────────────────────────────────────
echo ""
echo "=== Step 6: Creating combination jobscript ==="

cat > "/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/combine_${CONFIG_STR}.job" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=KSCombine_${CONFIG_STR:0:28}
#SBATCH --partition=nodes
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${LOG_BASE_DIR}/combination/Combine_${CONFIG_STR}_%j.log
#SBATCH --error=${LOG_BASE_DIR}/combination/Combine_${CONFIG_STR}_%j.err

set -e

module purge
module load ROOT/6.30.06-foss-2023a

echo "Starting file combination for ${CONFIG_STR}"

cd ${PARSING_DIR}

echo "Checking for all acceptance and vector files..."
for i in \$(seq 1 100); do
    while [ ! -f "${RESULTS_DIR}/${CONFIG_STR}_\${i}_acceptance.root" ]; do
        sleep 30
    done
    while [ ! -f "${RESULTS_DIR}/${CONFIG_STR}_\${i}_vectors.root" ]; do
        sleep 30
    done
done

echo "Combining acceptance files..."
hadd -f "${RESULTS_DIR}/${CONFIG_STR}_combined_acceptance.root" \
    ${RESULTS_DIR}/${CONFIG_STR}_[0-9]*_acceptance.root

echo "Combining vector files..."
hadd -f "${RESULTS_DIR}/${CONFIG_STR}_combined_vectors.root" \
    ${RESULTS_DIR}/${CONFIG_STR}_[0-9]*_vectors.root

echo ""
echo "=== Combination complete ==="
echo "  Acceptance : ${RESULTS_DIR}/${CONFIG_STR}_combined_acceptance.root"
echo "  Vectors    : ${RESULTS_DIR}/${CONFIG_STR}_combined_vectors.root"

rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_*_${CONFIG_STR}.job   || true
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/run_*_${CONFIG_STR}.job     || true
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/parse_file_*_${CONFIG_STR}.job || true
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/combine_${CONFIG_STR}.job   || true
EOF

echo "Created combination jobscript."

# ── Step 7: Create submission script ─────────────────────────────────────────
echo ""
echo "=== Step 7: Creating submission script ==="

cat > "/users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh" << EOF
#!/bin/bash

# KS Regeneration job submission script
# Configuration: ${CONFIG_STR}

set -e

echo "=== Submitting KS Regeneration jobs ==="
echo "KS beam z0 = ${KS_Z0} cm  |  Detector: ${DETECTOR_CONFIG}"

echo "Submitting 100 build jobs in sequential chain..."
BUILD_JOB_IDS=()
PREV_BUILD_ID=""

for run in \$(seq 1 100); do
    if [ -z "\$PREV_BUILD_ID" ]; then
        BUILD_ID=\$(sbatch /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_\${run}_${CONFIG_STR}.job | awk '{print \$4}')
        echo "Submitted build \$run: Job ID \$BUILD_ID"
    else
        BUILD_ID=\$(sbatch --dependency=afterok:\$PREV_BUILD_ID \
            /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_\${run}_${CONFIG_STR}.job | awk '{print \$4}')
        [ \$((\$run % 10)) -eq 0 ] && echo "Submitted builds \$((\$run-9))-\$run (last: \$BUILD_ID)"
    fi
    BUILD_JOB_IDS+=(\$BUILD_ID)
    PREV_BUILD_ID=\$BUILD_ID

    RUN_ID=\$(sbatch --dependency=afterok:\$BUILD_ID \
        /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/run_\${run}_${CONFIG_STR}.job | awk '{print \$4}')
    RUN_JOB_IDS+=(\$RUN_ID)
done

echo "All build/run jobs submitted."

RUN_DEPENDENCY=\$(IFS=:; echo "\${RUN_JOB_IDS[*]}")

echo "Submitting 100 parsing jobs..."
PARSE_JOB_IDS=()
for file_num in \$(seq 1 100); do
    JOB_ID=\$(sbatch --dependency=afterok:\$RUN_DEPENDENCY \
        /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/parse_file_\${file_num}_${CONFIG_STR}.job | awk '{print \$4}')
    PARSE_JOB_IDS+=(\$JOB_ID)
done
echo "Parsing jobs submitted."

PARSE_DEPENDENCY=\$(IFS=:; echo "\${PARSE_JOB_IDS[*]}")

COMBINE_JOB_ID=\$(sbatch --dependency=afterok:\$PARSE_DEPENDENCY \
    /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/combine_${CONFIG_STR}.job | awk '{print \$4}')
echo "Combination job: \$COMBINE_JOB_ID"

echo ""
echo "=== KS Regeneration submission complete ==="
echo "KS beam z0   : ${KS_Z0} cm"
echo "Results in   : ${RESULTS_DIR}"
echo "Monitor with : squeue -u \$USER"
echo ""
echo "Final outputs:"
echo "  ${CONFIG_STR}_combined_acceptance.root"
echo "  ${CONFIG_STR}_combined_vectors.root"

rm -f "\$0"
EOF

chmod +x "/users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh"
echo "Created submission script: submit_${CONFIG_STR}.sh"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "=== KS Regeneration pipeline setup complete ==="
echo ""
echo "  KS beam start      : z = ${KS_Z0} cm"
echo "  Particle           : kaon0S  (KS)"
echo "  Detector config    : ${DETECTOR_CONFIG}"
echo "  Full config string : ${CONFIG_STR}"
echo ""
echo "  Scenario dir : ${SCENARIO_DIR}"
echo "  Results dir  : ${RESULTS_DIR}"
echo "  Parsing dir  : ${PARSING_DIR}"
echo ""
echo "To submit the pipeline to SLURM, run:"
echo "  /users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh"
echo ""
echo "To scan multiple KS start positions, use the sweep helper:"
echo "  ./ks_regeneration_sweep.sh"

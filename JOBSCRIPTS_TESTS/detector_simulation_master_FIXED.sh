#!/bin/bash

# Master Detector Simulation Script - FIXED VERSION
# Usage: ./detector_simulation_master_FIXED.sh T1 T2 T3 T4 P1 P2 F1 F2 E1
# Where:
# T1-T4: Tracker positions (cm from dipole end)
# P1-P2: Pizza detector positions (cm from dipole end)  
# F1-F2: FRI wall positions (cm from dipole end)
# E1: TOF wall (End detector) position (cm from dipole end)

set -e

# Check if correct number of arguments provided
if [ $# -ne 9 ]; then
    echo "Usage: $0 T1 T2 T3 T4 P1 P2 F1 F2 E1"
    echo "Example: $0 250 260 270 280 220 230 290 300 400"
    echo "All positions are in cm from dipole magnet origin (z=0)"
    echo "IMPORTANT: Dipole extends from z=0 to z=200cm"
    echo "           All detectors must be positioned at z > 200cm"
    exit 1
fi

# Input parameters
T1=$1  # Tracker 1 position
T2=$2  # Tracker 2 position  
T3=$3  # Tracker 3 position
T4=$4  # Tracker 4 position
P1=$5  # Pizza 1 position
P2=$6  # Pizza 2 position
F1=$7  # FRI wall 1 position
F2=$8  # FRI wall 2 position
E1=$9  # TOF wall (End detector) position

# Geometry validation
DIPOLE_END=200  # Dipole magnet ends at z=200cm
MIN_DETECTOR_POS=210  # Minimum safe detector position

echo "=== Geometry Validation ==="
echo "Dipole magnet extends from z=0 to z=${DIPOLE_END}cm"
echo "Minimum detector position: z>${MIN_DETECTOR_POS}cm"

# Check all positions are safe
ALL_POSITIONS=($T1 $T2 $T3 $T4 $P1 $P2 $F1 $F2 $E1)
POSITION_NAMES=("T1" "T2" "T3" "T4" "P1" "P2" "F1" "F2" "E1")

for i in ${!ALL_POSITIONS[@]}; do
    pos=${ALL_POSITIONS[i]}
    name=${POSITION_NAMES[i]}
    if (( $(echo "$pos <= $MIN_DETECTOR_POS" | bc -l) )); then
        echo "ERROR: $name position ($pos cm) is too close to dipole magnet!"
        echo "       Must be > $MIN_DETECTOR_POS cm to avoid geometry overlaps"
        exit 1
    fi
done

echo "✓ All detector positions are geometrically valid"

# Create configuration string for filenames
CONFIG_STR="T1-${T1}_T2-${T2}_T3-${T3}_T4-${T4}_P1-${P1}_P2-${P2}_F1-${F1}_F2-${F2}_E1-${E1}"

echo "=== Master Detector Simulation Script - FIXED VERSION ==="
echo "Configuration: $CONFIG_STR"
echo "Tracker positions (cm): T1=$T1, T2=$T2, T3=$T3, T4=$T4"
echo "Pizza positions (cm): P1=$P1, P2=$P2"
echo "FRI positions (cm): F1=$F1, F2=$F2" 
echo "TOF position (cm): E1=$E1"

# Set paths
BASE_SCENARIO_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM"
SCENARIO_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/${CONFIG_STR}"
DETECTOR_FILE="${BASE_SCENARIO_DIR}/src/JLabKDetectorConstruction.cc"
BACKUP_FILE="${DETECTOR_FILE}.backup"
RESULTS_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/${CONFIG_STR}"
PARSING_DIR="/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/${CONFIG_STR}_DATA_READ"

# Create organized log directories and temp jobs directory
LOG_BASE_DIR="/users/bp969/scratch/JOB_LOGS/${CONFIG_STR}"
mkdir -p "${LOG_BASE_DIR}/build"
mkdir -p "${LOG_BASE_DIR}/simulation"  
mkdir -p "${LOG_BASE_DIR}/parsing"
mkdir -p "${LOG_BASE_DIR}/combination"
mkdir -p "${LOG_BASE_DIR}/histograms"
mkdir -p "/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs"
echo "📁 Created log directories in ${LOG_BASE_DIR}"

# Create simulation directory for this configuration
echo "Creating simulation directory: $SCENARIO_DIR"
mkdir -p "$SCENARIO_DIR"

# Copy entire SCENARIO_5_SIM directory to configuration-specific directory
if [ ! -d "$SCENARIO_DIR/src" ]; then
    echo "Copying simulation files to configuration directory..."
    cp -r "$BASE_SCENARIO_DIR"/* "$SCENARIO_DIR/"
    echo "Simulation files copied to $SCENARIO_DIR"
else
    echo "Configuration directory already exists with source files"
fi

# Create results directories
mkdir -p "$RESULTS_DIR"
mkdir -p "$PARSING_DIR"
mkdir -p "/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs"

# Update detector file path to use configuration-specific directory
DETECTOR_FILE="${SCENARIO_DIR}/src/JLabKDetectorConstruction.cc"
BACKUP_FILE="${DETECTOR_FILE}.backup"

# Cleanup behavior: set to "true" to remove intermediate files/logs after successful completion
# Set to "false" to keep all intermediate files for debugging
# Parse optional flags
CLEANUP_ON_SUCCESS=false
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --no-cleanup)
            CLEANUP_ON_SUCCESS=false
            shift
            ;;
        *)
            break
            ;;
    esac
done

echo "=== Step 1: Backing up and editing detector construction file ==="

# Create backup of original detector construction file (if needed for reference)
if [ ! -f "$BACKUP_FILE" ]; then
    cp "$DETECTOR_FILE" "$BACKUP_FILE"
    echo "Created backup: $BACKUP_FILE"
else
    echo "Backup already exists: $BACKUP_FILE"
fi

# Note: We do NOT restore from backup since the template now has correct positions
# and we want to preserve our edits for parallel pipeline execution

# FIXED function to edit detector positions in the construction file
edit_detector_positions() {
    local file="$1"
    
    echo "Editing detector positions in $file"
    
    # Show original positions before editing
    echo "=== Original positions in file ==="
    echo "Original tracker positions:"
    grep "G4double trkrPosZ\[4\]" "$file" | head -1 || echo "  Not found"
    echo "Original pizza positions:"
    grep "G4double pizzaPosZ\[4\]" "$file" | head -1 || echo "  Not found"
    echo "Original FRI positions:"
    grep "G4double friPosZ\[2\]" "$file" | head -1 || echo "  Not found"
    echo "Original TOF position:"
    grep "G4double tofPosZ =" "$file" | head -1 || echo "  Not found"
    
    # FIXED: Edit tracker positions (4 trackers) - corrected regex patterns
    sed -i "s|G4double trkrPosZ\[4\] = {[^}]*};|G4double trkrPosZ[4] = {${T1}.*cm, ${T2}.*cm, ${T3}.*cm, ${T4}.*cm};|" "$file"
    
    # FIXED: Edit pizza positions (2 pizzas) - only first two positions, keep downstream at 1100,1200
    sed -i "s|G4double pizzaPosZ\[4\] = {[^}]*};|G4double pizzaPosZ[4] = {${P1}.*cm, ${P2}.*cm, 1100.*cm, 1200.*cm};|" "$file"
    
    # FIXED: Edit FRI positions (2 FRI walls)
    sed -i "s|G4double friPosZ\[2\] = {[^}]*};|G4double friPosZ[2] = {${F1}.*cm, ${F2}.*cm};|" "$file"
    
    # FIXED: Edit TOF position (1 TOF wall) - corrected pattern
    sed -i "s|G4double tofPosZ = [^;]*;|G4double tofPosZ = ${E1}.*cm;|" "$file"
    
    echo "Detector positions updated successfully"
    
    # Verify changes were applied correctly
    echo "=== Verifying position updates ==="
    echo "Tracker positions:"
    grep "G4double trkrPosZ\[4\]" "$file" | head -1
    echo "Pizza positions:"
    grep "G4double pizzaPosZ\[4\]" "$file" | head -1
    echo "FRI positions:"
    grep "G4double friPosZ\[2\]" "$file" | head -1
    echo "TOF position:"
    grep "G4double tofPosZ =" "$file" | head -1
}

edit_detector_positions "$DETECTOR_FILE"

echo "=== Step 2: Creating batch macro and building Geant4 simulation ==="

# Create batch mode macro file for headless execution
cat > "${SCENARIO_DIR}/batch.mac" << 'EOF'
# Batch mode macro for headless execution
/control/verbose 0
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

# Set number of threads for parallel processing
/run/numberOfThreads 4

# Initialize the run
/run/initialize

# Run simulation with specified number of events (1M total = 100 jobs × 10,000 events)
/run/beamOn 10000
EOF

echo "Created batch.mac for headless simulation"

# Note: CMake configuration and build will be handled by individual build jobs
# where the Geant4 module is properly loaded
echo "Skipping initial build - will be done by build jobs with proper module environment"

echo "=== Step 3: Creating build and run jobscripts (100 builds sequential, 100 runs parallel) ==="

# Function to create build jobscript (edits source and rebuilds)
create_build_jobscript() {
    local run_num=$1
    
    local jobscript="/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_${run_num}_${CONFIG_STR}.job"
    
    cat > "$jobscript" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=Build_${CONFIG_STR}_${run_num}
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

# Load required modules for Geant4 simulation
module purge
module load Geant4/11.2.2-GCC-12.3.0
module load ROOT/6.30.06-foss-2023a

OUTFILE="${CONFIG_STR}_${run_num}.root"
RUNACTION_FILE="${SCENARIO_DIR}/src/JLabKRunAction.cc"
EXECUTABLE_DEST="${SCENARIO_DIR}/jlabk_${run_num}"

echo "=== Build ${run_num} starting at \$(date) ==="
echo "Output filename: \$OUTFILE"

cd ${SCENARIO_DIR}

# For first build, clean CMake cache and reconfigure
if [ ${run_num} -eq 1 ]; then
    echo "First build: Cleaning CMake cache and reconfiguring..."
    rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake
    cmake .
fi

# Edit output filename in source code using sed - FIXED pattern
echo "Setting output filename to \$OUTFILE in source code..."
sed -i 's|G4String fileName = "[^"]*";|G4String fileName = "'\$OUTFILE'";|' "\$RUNACTION_FILE"

# Rebuild the executable
echo "Rebuilding simulation executable..."
make -j4

# Copy executable to unique name so parallel runs don't interfere
cp jlabk "\$EXECUTABLE_DEST"
echo "Created executable: \$EXECUTABLE_DEST"

echo "=== Build ${run_num} completed at \$(date) ==="
EOF

    echo "Created build jobscript ${run_num}: $jobscript"
}

# Function to create run jobscript (runs the simulation)
create_run_jobscript() {
    local run_num=$1
    
    local jobscript="/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/run_${run_num}_${CONFIG_STR}.job"
    
    cat > "$jobscript" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=Run_${CONFIG_STR}_${run_num}
#SBATCH --partition=nodes
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${LOG_BASE_DIR}/simulation/Run_${CONFIG_STR}_${run_num}_%j.log
#SBATCH --error=${LOG_BASE_DIR}/simulation/Run_${CONFIG_STR}_${run_num}_%j.err

set -e

# Load required modules for Geant4 simulation
module purge
module load Geant4/11.2.2-GCC-12.3.0
module load ROOT/6.30.06-foss-2023a

SEED=\$((1000000 + ${run_num}))
OUTFILE="${CONFIG_STR}_${run_num}.root"
EXECUTABLE="${SCENARIO_DIR}/jlabk_${run_num}"

echo "=== Simulation ${run_num} starting at \$(date) ==="
echo "Seed: \$SEED"
echo "Output: \$OUTFILE"
echo "Executable: \$EXECUTABLE"

cd ${SCENARIO_DIR}

# Check executable exists
if [ ! -f "\$EXECUTABLE" ]; then
    echo "ERROR: Executable \$EXECUTABLE not found!"
    exit 1
fi

# Run simulation
echo "Running simulation..."
\$EXECUTABLE batch.mac \$SEED
EXIT_CODE=\$?

# Check output file was created with reasonable size
if [ -f "\$OUTFILE" ]; then
    FILESIZE=\$(stat -f%z "\$OUTFILE" 2>/dev/null || stat -c%s "\$OUTFILE" 2>/dev/null)
    if [ "\$FILESIZE" -gt 1000000 ]; then
        # File exists and is large enough (>1MB) - move to results
        mv "\$OUTFILE" "${RESULTS_DIR}/\$OUTFILE"
        echo "Completed simulation ${run_num}: \$OUTFILE (\${FILESIZE} bytes)"
        
        if [ \$EXIT_CODE -eq 139 ]; then
            echo "NOTE: Exit code 139 (cleanup segfault) ignored - output file created successfully"
        fi
        
        # Cleanup: remove unique executable after successful run
        rm -f "\$EXECUTABLE"
        echo "Cleaned up executable: \$EXECUTABLE"
    else
        echo "ERROR: Output file is too small (\${FILESIZE} bytes, need >1MB) - simulation likely failed"
        rm -f "\$EXECUTABLE"
        exit 1
    fi
else
    echo "ERROR: Output file '\$OUTFILE' not created (exit code: \$EXIT_CODE)"
    rm -f "\$EXECUTABLE"
    exit 1
fi

echo "=== Simulation ${run_num} completed at \$(date) ==="
EOF

    echo "Created run jobscript ${run_num}: $jobscript"
}

# Create 100 build and run jobscripts
for run in $(seq 1 100); do
    create_build_jobscript $run
    create_run_jobscript $run
done

echo "=== Step 4: Creating data parsing jobscripts ==="

# Copy parsing scripts to parsing directory
cp "/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C" "$PARSING_DIR/"
cp "/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C" "$PARSING_DIR/"

# Function to create parsing jobscript
create_parsing_jobscript() {
    local group_start=$1
    local group_end=$2
    local group_num=$3
    
    local jobscript="/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/parse_group_${group_num}_${CONFIG_STR}.job"
    
    cat > "$jobscript" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=Parse_${CONFIG_STR}_G${group_num}
#SBATCH --partition=nodes
#SBATCH --time=1-23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=4G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${LOG_BASE_DIR}/parsing/Parse_${CONFIG_STR}_G${group_num}_%j.log
#SBATCH --error=${LOG_BASE_DIR}/parsing/Parse_${CONFIG_STR}_G${group_num}_%j.err

set -e

echo "Starting parsing group ${group_num} (files ${group_start}-${group_end})"

module purge
module load ROOT/6.30.06-foss-2023a

cd ${PARSING_DIR}

# Wait for simulation files to be ready
echo "Checking for simulation files..."
for i in \$(seq ${group_start} ${group_end}); do
    INPUT_FILE="${RESULTS_DIR}/${CONFIG_STR}_\${i}.root"
    while [ ! -f "\$INPUT_FILE" ]; do
        echo "Waiting for \$INPUT_FILE..."
        sleep 30
    done
done

# Parallel parsing (restored - plotting script fix resolved the data access issue)
for i in \$(seq ${group_start} ${group_end}); do
    (
        INPUT_FILE="${RESULTS_DIR}/${CONFIG_STR}_\${i}.root"
        ACCEPTANCE_OUTPUT="${RESULTS_DIR}/${CONFIG_STR}_\${i}_acceptance.root"
        VECTORS_OUTPUT="${RESULTS_DIR}/${CONFIG_STR}_\${i}_vectors.root"
        
        echo "Processing file \$i: \$INPUT_FILE"
        
        # Run acceptance analysis (creates output in RESULTS_DIR)
        root -l -b -q "KLong_save_momentum_acceptance.C(\"\$INPUT_FILE\")"
        
        # Run vector analysis (creates output in RESULTS_DIR)
        root -l -b -q "KLong_save_vectors.C(\"\$INPUT_FILE\")"
        
        # Check that output files were created (they're already in RESULTS_DIR)
        if [ -f "\$ACCEPTANCE_OUTPUT" ]; then
            echo "Confirmed \$ACCEPTANCE_OUTPUT created successfully"
        else
            echo "ERROR: Expected output file \$ACCEPTANCE_OUTPUT not created!"
        fi
        
        if [ -f "\$VECTORS_OUTPUT" ]; then
            echo "Confirmed \$VECTORS_OUTPUT created successfully"
        else
            echo "ERROR: Expected output file \$VECTORS_OUTPUT not created!"
        fi
        
        echo "Completed parsing for file \$i"
        # Cleanup: remove raw simulation file after successful parsing (if enabled)
        if [ "$CLEANUP_ON_SUCCESS" = true ] && [ -f "\$INPUT_FILE" ]; then
            echo "Removing raw simulation file: \$INPUT_FILE"
            rm -f "\$INPUT_FILE"
        fi
    ) &
done

wait
echo "Completed parsing group ${group_num}"
EOF

    echo "Created parsing jobscript: $jobscript"
}

# Create 10 parsing jobscripts (10 files per script)
for group in $(seq 1 10); do
    start=$(( (group-1)*10 + 1 ))
    end=$(( group*10 ))
    create_parsing_jobscript $start $end $group
done

echo "=== Step 5: Creating combination jobscript ==="

cat > "/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/combine_${CONFIG_STR}.job" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=Combine_${CONFIG_STR}
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

# Wait for all acceptance files
echo "Checking for all acceptance files..."
for i in \$(seq 1 100); do
    ACCEPTANCE_FILE="${RESULTS_DIR}/${CONFIG_STR}_\${i}_acceptance.root"
    while [ ! -f "\$ACCEPTANCE_FILE" ]; do
        echo "Waiting for \$ACCEPTANCE_FILE..."
        sleep 30
    done
done

# Wait for all vector files
echo "Checking for all vector files..."
for i in \$(seq 1 100); do
    VECTORS_FILE="${RESULTS_DIR}/${CONFIG_STR}_\${i}_vectors.root"
    while [ ! -f "\$VECTORS_FILE" ]; do
        echo "Waiting for \$VECTORS_FILE..."
        sleep 30
    done
done

# Combine acceptance files
echo "Combining acceptance files..."
hadd -f "${RESULTS_DIR}/${CONFIG_STR}_combined_acceptance.root" ${RESULTS_DIR}/${CONFIG_STR}_[0-9]*_acceptance.root

# Combine vector files  
echo "Combining vector files..."
hadd -f "${RESULTS_DIR}/${CONFIG_STR}_combined_vectors.root" ${RESULTS_DIR}/${CONFIG_STR}_[0-9]*_vectors.root

if [ "$CLEANUP_ON_SUCCESS" = true ]; then
    # Cleanup per-run acceptance and vector files after successful combination
    echo "Cleaning up individual acceptance and vector files..."
    rm -f ${RESULTS_DIR}/${CONFIG_STR}_[0-9]*_acceptance.root
    rm -f ${RESULTS_DIR}/${CONFIG_STR}_[0-9]*_vectors.root

    # Cleanup build job logs and error files
    echo "Cleaning up build job logs..."
    find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Build_${CONFIG_STR}_*.log" -delete || true
    find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Build_${CONFIG_STR}_*.err" -delete || true

    # Cleanup run job logs and error files
    echo "Cleaning up run job logs..."
    find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Run_${CONFIG_STR}_*.log" -delete || true
    find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Run_${CONFIG_STR}_*.err" -delete || true

    # Cleanup parsing job logs
    echo "Cleaning up parsing job logs..."
    find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Parse_${CONFIG_STR}_G*.log" -delete || true
    find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Parse_${CONFIG_STR}_G*.err" -delete || true
else
    echo "CLEANUP_ON_SUCCESS=false: keeping individual run outputs and logs"
fi

echo "File combination completed successfully!"
echo "Combined acceptance file: ${RESULTS_DIR}/${CONFIG_STR}_combined_acceptance.root"
echo "Combined vectors file: ${RESULTS_DIR}/${CONFIG_STR}_combined_vectors.root"

# Cleanup combine job logs
echo "Cleaning up combine job logs..."
find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Combine_${CONFIG_STR}*.log" -delete || true
find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "Combine_${CONFIG_STR}*.err" -delete || true

# Remove temporary job scripts
echo "Removing temporary job scripts for ${CONFIG_STR}..."
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_*_${CONFIG_STR}.job || true
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/run_*_${CONFIG_STR}.job || true
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/parse_group_*_${CONFIG_STR}.job || true
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/combine_${CONFIG_STR}.job || true
EOF

echo "=== Step 6: Creating master job submission script ==="

cat > "/users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh" << EOF
#!/bin/bash

# Master job submission script for configuration: ${CONFIG_STR}

set -e

echo "=== Submitting jobs for detector configuration: ${CONFIG_STR} ==="

# Submit build jobs in a sequential chain (to avoid source conflicts)
echo "Submitting 100 build jobs in sequential chain..."
BUILD_JOB_IDS=()
PREV_BUILD_ID=""

for run in \$(seq 1 100); do
    if [ -z "\$PREV_BUILD_ID" ]; then
        # First build has no dependency
        BUILD_ID=\$(sbatch /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_\${run}_${CONFIG_STR}.job | awk '{print \$4}')
        echo "Submitted build \$run: Job ID \$BUILD_ID (no dependency)"
    else
        # Subsequent builds depend on previous build completion
        BUILD_ID=\$(sbatch --dependency=afterok:\$PREV_BUILD_ID /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/build_\${run}_${CONFIG_STR}.job | awk '{print \$4}')
        [ \$((\$run % 10)) -eq 0 ] && echo "Submitted builds \$((\$run-9))-\$run (last: \$BUILD_ID)"
    fi
    BUILD_JOB_IDS+=(\$BUILD_ID)
    PREV_BUILD_ID=\$BUILD_ID
    
    # Submit corresponding run job (depends only on its own build completing)
    RUN_ID=\$(sbatch --dependency=afterok:\$BUILD_ID /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/run_\${run}_${CONFIG_STR}.job | awk '{print \$4}')
    
    # Store run job ID for later dependency
    if [ \$run -eq 1 ]; then
        RUN_JOB_IDS=(\$RUN_ID)
    else
        RUN_JOB_IDS+=(\$RUN_ID)
    fi
done

echo "All build jobs submitted in sequential chain"
echo "All run jobs submitted in parallel (each depends on its own build)"

# Create dependency string for parsing jobs (wait for all run jobs)
RUN_DEPENDENCY=\$(IFS=:; echo "\${RUN_JOB_IDS[*]}")

# Submit parsing jobs with dependency on all run jobs
echo "Submitting parsing jobs (dependent on all run jobs)..."
PARSE_JOB_IDS=()
for group in \$(seq 1 10); do
    JOB_ID=\$(sbatch --dependency=afterok:\$RUN_DEPENDENCY /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/parse_group_\${group}_${CONFIG_STR}.job | awk '{print \$4}')
    PARSE_JOB_IDS+=(\$JOB_ID)
    echo "Submitted parsing group \$group: Job ID \$JOB_ID"
done

# Create dependency string for combination job (wait for all parsing jobs)
PARSE_DEPENDENCY=\$(IFS=:; echo "\${PARSE_JOB_IDS[*]}")

# Submit combination job with dependency on parsing jobs
echo "Submitting combination job (dependent on parsing jobs)..."
COMBINE_JOB_ID=\$(sbatch --dependency=afterok:\$PARSE_DEPENDENCY /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/combine_${CONFIG_STR}.job | awk '{print \$4}')
echo "Submitted combination job: Job ID \$COMBINE_JOB_ID"

echo
echo "=== Job submission completed! ==="
echo "Build jobs: \${BUILD_JOB_IDS[@]}"
echo "Run jobs: \${RUN_JOB_IDS[@]}"
echo "Parsing jobs: \${PARSE_JOB_IDS[@]}"
echo "Combination job: \$COMBINE_JOB_ID"
echo
echo "NOTE: Build jobs run sequentially (~10 min each) to avoid source conflicts"
echo "      Run jobs start in parallel as soon as their build completes (~2 days each for 100k events)"
echo "      Total wall time: ~17 minutes (builds) + ~2 days (parallel runs) for 10M total events"
echo
echo "Monitor progress with: squeue -u \$USER"
echo "Results will be in: ${RESULTS_DIR}"
EOF

chmod +x "/users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh"

echo
echo "=== FIXED Setup completed successfully! ==="
echo
echo "Configuration: $CONFIG_STR"
echo "Results directory: $RESULTS_DIR"
echo "Parsing directory: $PARSING_DIR"
echo
echo "FIXED ISSUES:"
echo "  ✓ Corrected SCENARIO_5_SIM template with proper detector positions"
echo "  ✓ Fixed sed patterns to use | delimiters instead of /"
echo "  ✓ Added verification output to confirm position updates"
echo "  ✓ Template now starts with correct baseline positions"
echo
echo "To start the FIXED simulation pipeline, run:"
echo "  /users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh"
echo
echo "This will:"
echo "  1. Submit 100 build jobs sequentially (each edits source, rebuilds, creates unique executable)"
echo "  2. Submit 100 run jobs in parallel (each starts as soon as its build completes)"
echo "  3. Submit 10 parsing jobs (dependent on all runs completing)"
echo "  4. Submit 1 combination job (dependent on parsing completion)"
echo
echo "NOTE: Builds run sequentially to avoid conflicts, but runs execute in parallel for speed"
echo "      Total wall time: ~17 minutes (sequential builds) + ~2 days (parallel runs)"
echo "      TOTAL EVENTS: 10 million (100 jobs × 100,000 events each)"
echo
echo "Final outputs will be:"
echo "  - ${CONFIG_STR}_combined_acceptance.root"
echo "  - ${CONFIG_STR}_combined_vectors.root"

echo
echo "NOTE: Detector construction file in ${SCENARIO_DIR} has been updated with your specified positions."
echo "      Each pipeline creates its own configuration directory for parallel execution."
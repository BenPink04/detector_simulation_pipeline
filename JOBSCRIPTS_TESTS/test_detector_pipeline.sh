#!/bin/bash

# Test script for detector simulation pipeline
# This creates a smaller test run (5 simulations instead of 100)

set -e

echo "=== Testing Detector Simulation Pipeline ==="
echo "This test runs 5 simulations instead of 100 for verification"
echo

# Test parameters (SAFE detector positions - after dipole magnet ends at z=200cm)
T1=250  # Tracker 1 - safe distance after dipole
T2=260  # Tracker 2
T3=270  # Tracker 3  
T4=280  # Tracker 4
P1=220  # Pizza 1 - after dipole end
P2=230  # Pizza 2
F1=290  # FRI 1 - after trackers
F2=300  # FRI 2
E1=400  # TOF wall - furthest downstream

CONFIG_STR="T1-${T1}_T2-${T2}_T3-${T3}_T4-${T4}_P1-${P1}_P2-${P2}_F1-${F1}_F2-${F2}_E1-${E1}"
echo "Test configuration: $CONFIG_STR"

# Set paths
SCENARIO_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM"
DETECTOR_FILE="${SCENARIO_DIR}/src/JLabKDetectorConstruction.cc"
BACKUP_FILE="${DETECTOR_FILE}.test_backup"
RESULTS_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/TEST_${CONFIG_STR}"
PARSING_DIR="/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/TEST_${CONFIG_STR}_DATA_READ"

# Create test directories
mkdir -p "$RESULTS_DIR"
mkdir -p "$PARSING_DIR"

echo "Test directories created:"
echo "  Results: $RESULTS_DIR"
echo "  Parsing: $PARSING_DIR"

# Test detector file editing
echo
echo "=== Testing detector file editing ==="

# Create backup
cp "$DETECTOR_FILE" "$BACKUP_FILE"

# Test the sed commands that will be used
echo "Testing tracker position edit..."
sed -i "s/G4double trkrPosZ\[4\] = {[^}]*};/G4double trkrPosZ[4] = {${T1}.*cm, ${T2}.*cm, ${T3}.*cm, ${T4}.*cm};/" "$DETECTOR_FILE"

echo "Testing pizza position edit..."
sed -i "s/G4double pizzaPosZ\[4\] = {[^}]*};/G4double pizzaPosZ[4] = {${P1}.*cm, ${P2}.*cm, 1100.*cm, 1200.*cm};/" "$DETECTOR_FILE"

echo "Testing FRI position edit..."
sed -i "s/G4double friPosZ\[2\] = {[^}]*};/G4double friPosZ[2] = {${F1}.*cm, ${F2}.*cm};/" "$DETECTOR_FILE"

echo "Testing TOF position edit..."
sed -i "s/G4double tofPosZ = [^;]*;/G4double tofPosZ = ${E1}.*cm;/" "$DETECTOR_FILE"

# Check what changed
echo "Checking edited positions:"
grep "trkrPosZ\[4\]" "$DETECTOR_FILE" || echo "Tracker positions not found"
grep "pizzaPosZ\[4\]" "$DETECTOR_FILE" || echo "Pizza positions not found"  
grep "friPosZ\[2\]" "$DETECTOR_FILE" || echo "FRI positions not found"
grep "tofPosZ.*=" "$DETECTOR_FILE" || echo "TOF position not found"

# Restore original file
cp "$BACKUP_FILE" "$DETECTOR_FILE"
echo "Original detector file restored"

# Create test parsing scripts
echo
echo "=== Creating test parsing scripts ==="
cp "/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C" "$PARSING_DIR/"
cp "/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C" "$PARSING_DIR/"
echo "Parsing scripts copied to test directory"

# Create a simple test jobscript (single simulation)
echo
echo "=== Creating test simulation jobscript ==="
cat > "/users/bp969/scratch/JOBSCRIPTS_TESTS/test_sim_${CONFIG_STR}.job" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=TEST_Sim_${CONFIG_STR}
#SBATCH --partition=nodes
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

set -e

# Load required modules for Geant4 simulation
module purge
module load Geant4/11.2.2-GCC-12.3.0
module load ROOT/6.30.06-foss-2023a

echo "TEST: Starting single simulation for configuration: ${CONFIG_STR}"
echo "Loaded modules: \$(module list 2>&1 | grep -E '(Geant4|ROOT)' || echo 'None')"
cd ${SCENARIO_DIR}

# Apply detector configuration
sed -i "s/G4double trkrPosZ\[4\] = {[^}]*};/G4double trkrPosZ[4] = {${T1}.*cm, ${T2}.*cm, ${T3}.*cm, ${T4}.*cm};/" src/JLabKDetectorConstruction.cc
sed -i "s/G4double pizzaPosZ\[4\] = {[^}]*};/G4double pizzaPosZ[4] = {${P1}.*cm, ${P2}.*cm, 1100.*cm, 1200.*cm};/" src/JLabKDetectorConstruction.cc
sed -i "s/G4double friPosZ\[2\] = {[^}]*};/G4double friPosZ[2] = {${F1}.*cm, ${F2}.*cm};/" src/JLabKDetectorConstruction.cc
sed -i "s/G4double tofPosZ = [^;]*;/G4double tofPosZ = ${E1}.*cm;/" src/JLabKDetectorConstruction.cc

# Build if needed
if [ ! -f "jlabk" ] || [ "src/JLabKDetectorConstruction.cc" -nt "jlabk" ]; then
    echo "Building Geant4 simulation..."
    make clean
    make -j4
fi

# Run simulation
SEED=1000001
OUTFILE="${CONFIG_STR}_1.root"

if [ -f "run.mac" ]; then
    echo "Running simulation with seed \$SEED"
    ./jlabk run.mac \$SEED
    
    if [ -f "out.root" ]; then
        mv out.root "${RESULTS_DIR}/\$OUTFILE"
        echo "SUCCESS: Test simulation completed -> \$OUTFILE"
    else
        echo "ERROR: No output file generated"
        exit 1
    fi
else
    echo "ERROR: run.mac not found"
    exit 1
fi

echo "TEST: Simulation completed successfully"
EOF

echo "Test simulation jobscript created"

# Create test parsing jobscript
echo "=== Creating test parsing jobscript ==="
cat > "/users/bp969/scratch/JOBSCRIPTS_TESTS/test_parse_${CONFIG_STR}.job" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=TEST_Parse_${CONFIG_STR}
#SBATCH --partition=nodes
#SBATCH --time=30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

set -e

echo "TEST: Starting parsing for configuration: ${CONFIG_STR}"

module purge
module load ROOT/6.30.06-foss-2023a

cd ${PARSING_DIR}

INPUT_FILE="${RESULTS_DIR}/${CONFIG_STR}_1.root"

# Wait for input file
while [ ! -f "\$INPUT_FILE" ]; do
    echo "Waiting for \$INPUT_FILE..."
    sleep 10
done

echo "Processing \$INPUT_FILE"

# Run analyses
root -l -b -q "KLong_save_momentum_acceptance.C(\\"\\$INPUT_FILE\\")" 
root -l -b -q "KLong_save_vectors.C(\\"\\$INPUT_FILE\\")"

echo "TEST: Parsing completed successfully"

# Check output files
ACCEPTANCE_FILE="${RESULTS_DIR}/${CONFIG_STR}_1_acceptance.root"
VECTORS_FILE="${RESULTS_DIR}/${CONFIG_STR}_1_vectors.root"

if [ -f "\$ACCEPTANCE_FILE" ]; then
    echo "SUCCESS: Acceptance file created -> \$ACCEPTANCE_FILE"
else
    echo "WARNING: Acceptance file not found"
fi

if [ -f "\$VECTORS_FILE" ]; then
    echo "SUCCESS: Vectors file created -> \$VECTORS_FILE"
else
    echo "WARNING: Vectors file not found"
fi

# Cleanup test logs and input file
echo "Cleaning up test logs and input files..."
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/TEST_Sim_*_${CONFIG_STR}-*.log 2>/dev/null || true
rm -f /users/bp969/scratch/JOBSCRIPTS_TESTS/TEST_Sim_*_${CONFIG_STR}-*.err 2>/dev/null || true
rm -f "${RESULTS_DIR}/${CONFIG_STR}_1.root" 2>/dev/null || true
EOF

echo "Test parsing jobscript created"

echo
echo "=== Test Setup Complete ==="
echo 
echo "To run the test:"
echo "1. Submit test simulation:"
echo "   sbatch /users/bp969/scratch/JOBSCRIPTS_TESTS/test_sim_${CONFIG_STR}.job"
echo
echo "2. After simulation completes, submit parsing:"
echo "   sbatch /users/bp969/scratch/JOBSCRIPTS_TESTS/test_parse_${CONFIG_STR}.job"
echo
echo "3. Check results in:"
echo "   $RESULTS_DIR"
echo
echo "If the test works, you can run the full pipeline with:"
echo "   ./detector_simulation_master.sh $T1 $T2 $T3 $T4 $P1 $P2 $F1 $F2 $E1"

# Clean up
rm "$BACKUP_FILE" 2>/dev/null || true

echo
echo "Test setup completed successfully!"
#!/bin/bash

# Geometry Validation and Fix Script for Detector Simulation
# This script helps debug and fix geometry issues

echo "=== Detector Geometry Validation Tool ==="
echo

# Define geometry constraints
DIPOLE_START=0
DIPOLE_END=200
MIN_SAFE_POS=210
GUIDE_TUBE_POS=850  # Guide tube is at z=850cm

echo "Geometry constraints:"
echo "  Dipole magnet: z = ${DIPOLE_START} to ${DIPOLE_END} cm"
echo "  Guide tube: z = ${GUIDE_TUBE_POS} cm" 
echo "  Minimum safe detector position: z > ${MIN_SAFE_POS} cm"
echo

# Function to validate a detector configuration
validate_config() {
    local T1=$1 T2=$2 T3=$3 T4=$4 P1=$5 P2=$6 F1=$7 F2=$8 E1=$9
    
    echo "Validating configuration: T1=$T1 T2=$T2 T3=$T3 T4=$T4 P1=$P1 P2=$P2 F1=$F1 F2=$F2 E1=$E1"
    
    local valid=true
    local positions=($T1 $T2 $T3 $T4 $P1 $P2 $F1 $F2 $E1)
    local names=("T1" "T2" "T3" "T4" "P1" "P2" "F1" "F2" "E1")
    
    for i in ${!positions[@]}; do
        pos=${positions[i]}
        name=${names[i]}
        
        if (( $(echo "$pos <= $MIN_SAFE_POS" | bc -l) )); then
            echo "  ❌ $name=$pos cm (TOO CLOSE to dipole, must be > $MIN_SAFE_POS cm)"
            valid=false
        elif (( $(echo "$pos > $GUIDE_TUBE_POS" | bc -l) )); then
            echo "  ⚠️  $name=$pos cm (beyond guide tube at $GUIDE_TUBE_POS cm)"
        else
            echo "  ✅ $name=$pos cm (safe)"
        fi
    done
    
    # Check for detector overlaps (assuming 5cm thickness per detector)
    local sorted_pos=($(printf '%s\n' "${positions[@]}" | sort -n))
    for i in $(seq 0 $((${#sorted_pos[@]}-2))); do
        pos1=${sorted_pos[i]}
        pos2=${sorted_pos[i+1]}
        if (( $(echo "$pos2 - $pos1 < 5" | bc -l) )); then
            echo "  ❌ Potential overlap: detectors at $pos1 cm and $pos2 cm are too close (< 5cm apart)"
            valid=false
        fi
    done
    
    if [ "$valid" = true ]; then
        echo "  ✅ Configuration is geometrically valid!"
        return 0
    else
        echo "  ❌ Configuration has geometry issues!"
        return 1
    fi
}

# Function to suggest safe configurations
suggest_configs() {
    echo "Suggested safe detector configurations:"
    echo
    
    echo "1. Compact setup (good resolution):"
    echo "   ./detector_simulation_master.sh 250 260 270 280 220 230 290 300 350"
    validate_config 250 260 270 280 220 230 290 300 350
    echo
    
    echo "2. Spread out setup (wide acceptance):"
    echo "   ./detector_simulation_master.sh 300 350 400 450 220 250 500 550 600"
    validate_config 300 350 400 450 220 250 500 550 600
    echo
    
    echo "3. Minimal separation (debugging):"
    echo "   ./detector_simulation_master.sh 215 225 235 245 220 230 255 265 275"
    validate_config 215 225 235 245 220 230 255 265 275
    echo
}

# Function to fix the current detector construction file
fix_detector_file() {
    local detector_file="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM/src/JLabKDetectorConstruction.cc"
    local backup_file="${detector_file}.backup"
    
    echo "=== Fixing detector construction file ==="
    
    if [ ! -f "$backup_file" ]; then
        echo "No backup file found. Cannot restore original."
        return 1
    fi
    
    echo "Restoring original detector construction file..."
    cp "$backup_file" "$detector_file"
    
    echo "Setting safe default positions..."
    # Use working positions from other scenarios
    sed -i "s/G4double trkrPosZ\[4\] = {[^}]*};/G4double trkrPosZ[4] = {250.*cm, 260.*cm, 270.*cm, 280.*cm};/" "$detector_file"
    sed -i "s/G4double pizzaPosZ\[4\] = {[^}]*};/G4double pizzaPosZ[4] = {220.*cm, 230.*cm, 1100.*cm, 1200.*cm};/" "$detector_file"
    sed -i "s/G4double friPosZ\[2\] = {[^}]*};/G4double friPosZ[2] = {290.*cm, 300.*cm};/" "$detector_file"
    sed -i "s/G4double tofPosZ = [^;]*;/G4double tofPosZ = 350.*cm;/" "$detector_file"
    
    echo "✅ Detector construction file fixed with safe positions"
    echo "   Trackers: 250, 260, 270, 280 cm"
    echo "   Pizzas: 220, 230 cm"
    echo "   FRI: 290, 300 cm" 
    echo "   TOF: 350 cm"
}

# Function to test build
test_build() {
    echo "=== Testing Geant4 build ==="
    cd /users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM
    
    # Load required modules
    echo "Loading Geant4 and ROOT modules..."
    module purge
    module load Geant4/11.2.2-GCC-12.3.0
    module load ROOT/6.30.06-foss-2023a
    
    if [ -f "Makefile" ]; then
        echo "Cleaning previous build..."
        make clean > /dev/null 2>&1
        
        echo "Building with fixed geometry..."
        if make -j4 > build.log 2>&1; then
            echo "✅ Build successful!"
            if [ -f "jlabk" ]; then
                echo "✅ Executable 'jlabk' created"
            fi
        else
            echo "❌ Build failed. Check build.log for details."
            tail -10 build.log
            return 1
        fi
    else
        echo "❌ No Makefile found in scenario directory"
        return 1
    fi
}

# Main menu
if [ $# -eq 9 ]; then
    echo "Validating provided configuration..."
    validate_config $1 $2 $3 $4 $5 $6 $7 $8 $9
elif [ $# -eq 1 ] && [ "$1" = "fix" ]; then
    fix_detector_file
    test_build
elif [ $# -eq 1 ] && [ "$1" = "suggest" ]; then
    suggest_configs
else
    echo "Usage:"
    echo "  $0 T1 T2 T3 T4 P1 P2 F1 F2 E1    # Validate configuration"
    echo "  $0 fix                              # Fix detector file and test build"
    echo "  $0 suggest                          # Show suggested configurations"
    echo
    echo "Examples:"
    echo "  $0 250 260 270 280 220 230 290 300 350"
    echo "  $0 fix"
    echo "  $0 suggest"
fi
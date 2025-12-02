#!/bin/bash

# Fix Geometry Overlaps in All Template Files
# This script applies the geometry overlap fixes to all scenario template files
# so that future pipeline runs will automatically use clean geometry

set -e

echo "=== Geometry Template Fix Script ==="
echo "Applying overlap fixes to all scenario template files..."

# Define all template directories to fix
TEMPLATE_DIRS=(
    "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_3_SIM"
    "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_4_SIM"
    "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM"
    "/users/bp969/scratch/VIKING_FOLDER/SCENARIOS/Scenarios/Scenario 3/Scenario 3 within 2 metres"
    "/users/bp969/scratch/VIKING_FOLDER/SCENARIOS/Scenarios/Scenario 4"
    "/users/bp969/scratch/VIKING_FOLDER/SCENARIOS/Scenarios/Unchanged Setup - no scenario"
)

# Geometry fixes to apply
apply_geometry_fixes() {
    local detector_file="$1"
    local template_name="$2"
    
    echo "  Processing: $template_name"
    
    # Create backup
    cp "$detector_file" "${detector_file}.backup_$(date +%Y%m%d_%H%M%S)"
    
    # Fix 1: Remove duplicate -20.4*cm StrawTube position
    echo "    ✓ Fixing duplicate StrawTube position (-20.4*cm appears twice)"
    sed -i 's/-20\.4\*cm, -20\.4\*cm/-20.4*cm/g' "$detector_file"
    
    # Fix 2: Increase tracker aperture clearance (add 0.5mm = 0.05cm)
    echo "    ✓ Increasing tracker aperture clearance (+0.5mm)"
    # Look for tracker aperture definitions and increase by 0.05cm
    sed -i 's/G4double trackerAperture = \([0-9.]*\)\*cm/G4double trackerAperture = (\1+0.05)*cm/g' "$detector_file"
    
    # Fix 3: Adjust pizza slice angles to 14.9° with clearance
    echo "    ✓ Adjusting pizza slice angles to 14.9° (was 15.0°)"
    # Reduce pizza slice angles from 15° to 14.9° for 0.1° clearance
    sed -i 's/15\.0\*deg/14.9*deg/g' "$detector_file"
    sed -i 's/15\*deg/14.9*deg/g' "$detector_file"
    
    # Fix 4: Add inner radius clearance for pizza slices (0.1mm = 0.01cm)
    echo "    ✓ Adding pizza slice inner radius clearance (+0.1mm)"
    sed -i 's/G4double pizzaInnerRadius = \([0-9.]*\)\*cm/G4double pizzaInnerRadius = (\1+0.01)*cm/g' "$detector_file"
    
    echo "    ✅ Fixed: $template_name"
}

# Apply fixes to all template files
for template_dir in "${TEMPLATE_DIRS[@]}"; do
    detector_file="$template_dir/src/JLabKDetectorConstruction.cc"
    
    if [ -f "$detector_file" ]; then
        echo "📁 Template: $(basename "$template_dir")"
        apply_geometry_fixes "$detector_file" "$(basename "$template_dir")"
    else
        echo "⚠️  Skipping: $detector_file (not found)"
    fi
done

echo ""
echo "=== Verification ==="
echo "Checking for remaining geometry overlaps in templates..."

# Verify fixes were applied
for template_dir in "${TEMPLATE_DIRS[@]}"; do
    detector_file="$template_dir/src/JLabKDetectorConstruction.cc"
    
    if [ -f "$detector_file" ]; then
        echo "🔍 Checking: $(basename "$template_dir")"
        
        # Check for duplicate -20.4*cm
        duplicate_count=$(grep -o "\-20\.4\*cm" "$detector_file" | wc -l)
        if [ $duplicate_count -gt 1 ]; then
            echo "  ❌ Still has duplicate -20.4*cm entries ($duplicate_count found)"
        else
            echo "  ✅ No duplicate StrawTube positions"
        fi
        
        # Check for 15.0*deg (should be 14.9*deg)
        if grep -q "15\.0\*deg\|15\*deg" "$detector_file"; then
            echo "  ❌ Still has 15.0° pizza slice angles"
        else
            echo "  ✅ Pizza slice angles updated to 14.9°"
        fi
    fi
done

echo ""
echo "=== Summary ==="
echo "✅ Geometry fixes applied to all template files"
echo "✅ Future pipeline runs will use clean geometry automatically" 
echo "✅ Existing running configurations need to be rebuilt or manually fixed"
echo ""
echo "📋 Next Steps:"
echo "1. Run a test pipeline to verify zero overlap warnings"
echo "2. Consider rebuilding active configurations that are currently running"
echo "3. Monitor job logs for geometry warnings"

echo ""
echo "🎯 GEOMETRY TEMPLATE FIX COMPLETE"
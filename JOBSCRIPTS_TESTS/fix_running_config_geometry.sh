#!/bin/bash

# Fix Geometry Overlaps in Existing Running Configurations
# This script fixes geometry overlaps in already-created configuration directories
# Usage: ./fix_running_config_geometry.sh [CONFIG_NAME]
# Example: ./fix_running_config_geometry.sh T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600

set -e

CONFIG_NAME="$1"

if [ -z "$CONFIG_NAME" ]; then
    echo "Usage: $0 CONFIG_NAME"
    echo "Example: $0 T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600"
    echo ""
    echo "Available configurations:"
    ls -1 /users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/ | grep "^T1-" | head -10
    exit 1
fi

CONFIG_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/$CONFIG_NAME"
DETECTOR_FILE="$CONFIG_DIR/src/JLabKDetectorConstruction.cc"

if [ ! -f "$DETECTOR_FILE" ]; then
    echo "❌ Configuration not found: $CONFIG_DIR"
    exit 1
fi

echo "=== Fixing Geometry Overlaps in Running Configuration ==="
echo "Configuration: $CONFIG_NAME"
echo "Detector file: $DETECTOR_FILE"

# Create backup
BACKUP_FILE="${DETECTOR_FILE}.backup_$(date +%Y%m%d_%H%M%S)"
cp "$DETECTOR_FILE" "$BACKUP_FILE"
echo "✅ Backup created: $BACKUP_FILE"

# Check current overlap issues
echo ""
echo "🔍 Analyzing current geometry issues..."

# Count duplicate -20.4*cm entries
duplicate_count=$(grep -o "\-20\.4\*cm" "$DETECTOR_FILE" | wc -l)
echo "  StrawTube duplicate positions: $duplicate_count occurrences of -20.4*cm"

# Check pizza slice angles  
if grep -q "15\.0\*deg\|15\*deg" "$DETECTOR_FILE"; then
    echo "  Pizza slice angles: 15.0° (needs reduction to 14.9°)"
else
    echo "  Pizza slice angles: Already optimized"
fi

# Apply geometry fixes
echo ""
echo "🔧 Applying geometry fixes..."

# Fix 1: Remove duplicate -20.4*cm StrawTube position
echo "  1. Fixing duplicate StrawTube positions..."
sed -i.tmp1 's/-20\.4\*cm, -20\.4\*cm/-20.4*cm/g' "$DETECTOR_FILE"

# Fix 2: Reduce pizza slice angles for clearance
echo "  2. Adjusting pizza slice angles to 14.9°..."
sed -i.tmp2 's/15\.0\*deg/14.9*deg/g' "$DETECTOR_FILE"
sed -i.tmp3 's/15\*deg/14.9*deg/g' "$DETECTOR_FILE"

# Clean up intermediate files
rm -f "${DETECTOR_FILE}.tmp1" "${DETECTOR_FILE}.tmp2" "${DETECTOR_FILE}.tmp3"

# Verify fixes
echo ""
echo "✅ Verification of applied fixes..."

# Check duplicate fix
new_duplicate_count=$(grep -o "\-20\.4\*cm" "$DETECTOR_FILE" | wc -l)
if [ $new_duplicate_count -eq 1 ]; then
    echo "  ✅ StrawTube duplicates fixed (1 occurrence of -20.4*cm)"
elif [ $new_duplicate_count -gt 1 ]; then
    echo "  ⚠️  StrawTube duplicates partially fixed ($new_duplicate_count occurrences remain)"
else
    echo "  ❌ No -20.4*cm entries found (may need manual review)"
fi

# Check pizza slice fix
if grep -q "15\.0\*deg\|15\*deg" "$DETECTOR_FILE"; then
    echo "  ⚠️  Some 15.0° pizza slice angles remain"
else
    echo "  ✅ Pizza slice angles updated to 14.9°"
fi

# Show diff summary
echo ""
echo "📊 Changes made:"
diff_lines=$(diff "$BACKUP_FILE" "$DETECTOR_FILE" | wc -l)
echo "  Total lines changed: $diff_lines"

# If configuration is currently running, warn about rebuild
if [ -f "$CONFIG_DIR/CMakeCache.txt" ]; then
    echo ""
    echo "⚠️  WARNING: Configuration appears to be built/running"
    echo "   You may need to rebuild the simulation to apply geometry changes:"
    echo "   cd $CONFIG_DIR"
    echo "   make clean && make"
fi

echo ""
echo "🎯 GEOMETRY FIX COMPLETE FOR $CONFIG_NAME"
echo "📁 Backup available at: $BACKUP_FILE"
echo "🔄 Consider running a test simulation to verify zero overlap warnings"
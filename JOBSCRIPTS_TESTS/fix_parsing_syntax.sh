#!/bin/bash

echo "Fixing syntax errors in all parsing macros..."

# Find all C files in DATA_PARSING directory and fix the dangling << operators
find VIKING_FOLDER/DATA_PARSING -name "*.C" -type f | while read file; do
    echo "Fixing $file"
    
    # Fix the dangling << " | True p: " << true_p_mag; lines
    sed -i '/\/\/ std::cout.*Selected Event/N;s/\/\/ std::cout.*Selected Event.*\n[[:space:]]*<< " | True p: " << true_p_mag;/\/\/ std::cout << "Selected Event #" << event_counter << " | True p: " << true_p_mag;/' "$file"
    
    # Also fix any similar issues with vector files that might have dangling operators
    sed -i '/\/\/ std::cout.*Event/N;s/\/\/ std::cout.*Event.*\n[[:space:]]*<< " | True p: " << true_p_mag/\/\/ std::cout << "Event #" << event_counter << " | True p: " << true_p_mag/' "$file"
    
    # Fix any other dangling << operators by commenting them out completely
    sed -i 's/^[[:space:]]*<< " | True p: " << true_p_mag;$/\/\/ << " | True p: " << true_p_mag;/' "$file"
    sed -i 's/^[[:space:]]*<< " | True p: " << true_p_mag$/\/\/ << " | True p: " << true_p_mag/' "$file"
done

echo "Fixed syntax errors in all parsing macros"
echo "Verifying fixes..."

# Check for remaining syntax issues
echo "Checking for remaining dangling operators..."
if find VIKING_FOLDER/DATA_PARSING -name "*.C" -exec grep -l "^[[:space:]]*<<.*true_p_mag" {} \; | grep -v "//"; then
    echo "WARNING: Some files still have dangling operators"
else
    echo "✅ All dangling operators have been fixed"
fi
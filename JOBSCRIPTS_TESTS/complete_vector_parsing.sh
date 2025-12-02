#!/bin/bash

# Complete Vector Parsing Script for 10M Event Pipeline
# This script generates all missing vector files for the completed simulation

set -e

CONFIG_STR="T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600"
RESULTS_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/${CONFIG_STR}"
PARSING_DIR="/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING"

echo "=== COMPLETING VECTOR PARSING FOR 10M EVENT PIPELINE ==="
echo "Configuration: ${CONFIG_STR}"
echo "Results directory: ${RESULTS_DIR}"
echo "Starting vector generation for 100 simulation files..."

# Load ROOT module
module load ROOT/6.30.06-foss-2023a

cd ${PARSING_DIR}

# Check how many vector files already exist
EXISTING_VECTORS=$(ls ${RESULTS_DIR}/*_vectors.root 2>/dev/null | wc -l || echo 0)
echo "Found ${EXISTING_VECTORS} existing vector files"

# Generate missing vector files
PROCESSED=0
FAILED=0
START_TIME=$(date +%s)

for i in {1..100}; do
    SIMULATION_FILE="${RESULTS_DIR}/${CONFIG_STR}_${i}.root"
    VECTOR_FILE="${RESULTS_DIR}/${CONFIG_STR}_${i}_vectors.root"
    
    # Check if simulation file exists
    if [ ! -f "${SIMULATION_FILE}" ]; then
        echo "❌ Missing simulation file: ${CONFIG_STR}_${i}.root"
        ((FAILED++))
        continue
    fi
    
    # Check if vector file already exists
    if [ -f "${VECTOR_FILE}" ]; then
        echo "✅ Vector file already exists: ${CONFIG_STR}_${i}_vectors.root"
        ((PROCESSED++))
        continue
    fi
    
    echo "Processing vectors for file ${i}/100: ${CONFIG_STR}_${i}.root"
    
    # Run vector processing
    if root -l -b -q "KLong_save_vectors.C(\"${SIMULATION_FILE}\")" > /tmp/vector_${i}.log 2>&1; then
        if [ -f "${VECTOR_FILE}" ]; then
            FILE_SIZE=$(stat -f%z "${VECTOR_FILE}" 2>/dev/null || stat -c%s "${VECTOR_FILE}" 2>/dev/null)
            echo "✅ Created vector file: ${CONFIG_STR}_${i}_vectors.root (${FILE_SIZE} bytes)"
            ((PROCESSED++))
        else
            echo "❌ Vector processing failed for file ${i}: No output file created"
            ((FAILED++))
        fi
    else
        echo "❌ Vector processing failed for file ${i}: ROOT execution error"
        ((FAILED++))
    fi
    
    # Progress update every 10 files
    if [ $((i % 10)) -eq 0 ]; then
        ELAPSED=$(($(date +%s) - START_TIME))
        echo "Progress: ${i}/100 files processed (${PROCESSED} successful, ${FAILED} failed) - ${ELAPSED}s elapsed"
    fi
done

# Final summary
TOTAL_VECTORS=$(ls ${RESULTS_DIR}/*_vectors.root 2>/dev/null | wc -l || echo 0)
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo "🎯 ==================================="
echo "🎯 VECTOR PARSING COMPLETED!"
echo "🎯 ==================================="
echo ""
echo "Total simulation files processed: 100"
echo "Successful vector files created: ${PROCESSED}"
echo "Failed vector generations: ${FAILED}" 
echo "Total vector files available: ${TOTAL_VECTORS}"
echo "Total processing time: ${TOTAL_TIME} seconds"
echo ""

if [ ${TOTAL_VECTORS} -eq 100 ]; then
    echo "✅ SUCCESS: All 100 vector files are ready for combination!"
    echo "Next step: Run the combination script"
else
    echo "⚠️  WARNING: Only ${TOTAL_VECTORS}/100 vector files available"
    echo "Check failed files before proceeding to combination"
fi

echo ""
echo "Vector files location: ${RESULTS_DIR}/*_vectors.root"
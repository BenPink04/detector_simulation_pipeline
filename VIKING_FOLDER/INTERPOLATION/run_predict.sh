#!/bin/bash
# =============================================================================
# run_predict.sh
#
# Predict resolution and acceptance curves for a given pizza position,
# using previously-simulated archives as training data.
#
# USAGE:
#   ./run_predict.sh <archives_parent_dir> <target_P1> [t3_mode] [output_dir]
#
# ARGUMENTS:
#   archives_parent_dir  — directory containing multiple archive sub-dirs
#                          (each with *combined_vectors.root / *combined_acceptance.root)
#   target_P1            — pizza position to predict for [cm]
#   t3_mode              — "upstream" | "downstream" | "auto"  (default: auto)
#   output_dir           — where to save the PNG  (default: <parent_dir>/plots)
#
# EXAMPLE:
#   ./run_predict.sh /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS 350.0
#   ./run_predict.sh /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS 300.0 downstream
#   ./run_predict.sh /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS 420.0 upstream /tmp/out
#
# OUTPUT:
#   <output_dir>/KLong_predicted_P1<N>_<mode>.png
#     Left panel:  FWHM(p) = A * p^n  with training curves (dashed) for context
#     Right panel: double-logistic acceptance with training curves (dashed)
# =============================================================================

set -e

PARENT_DIR="${1:-}"
TARGET_P1="${2:-}"
T3_MODE="${3:-auto}"
OUTPUT_DIR="${4:-}"

if [[ -z "$PARENT_DIR" || -z "$TARGET_P1" ]]; then
    echo "Usage: $0 <archives_parent_dir> <target_P1> [t3_mode] [output_dir]"
    echo "  t3_mode: upstream | downstream | auto (default: auto)"
    exit 1
fi

if [[ ! -d "$PARENT_DIR" ]]; then
    echo "ERROR: directory not found: $PARENT_DIR"
    exit 1
fi

# Validate t3_mode
if [[ "$T3_MODE" != "upstream" && "$T3_MODE" != "downstream" && \
      "$T3_MODE" != "no_T3"    && "$T3_MODE" != "auto" ]]; then
    echo "ERROR: t3_mode must be 'upstream', 'downstream', 'no_T3', or 'auto' (got: $T3_MODE)"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREDICT_MACRO="${SCRIPT_DIR}/KLong_predict_from_geometry.C"

if [[ ! -f "$PREDICT_MACRO" ]]; then
    echo "ERROR: macro not found: $PREDICT_MACRO"
    exit 1
fi

# Check for at least one training file
N_VEC=$(find "$PARENT_DIR" -maxdepth 2 -name '*combined_vectors.root' 2>/dev/null | wc -l)
N_ACC=$(find "$PARENT_DIR" -maxdepth 2 -name '*combined_acceptance.root' 2>/dev/null | wc -l)
echo "  Training data: $N_VEC vector files, $N_ACC acceptance files found under $PARENT_DIR"
if [[ "$N_VEC" -eq 0 && "$N_ACC" -eq 0 ]]; then
    echo "ERROR: no *combined_vectors.root or *combined_acceptance.root files found."
    exit 1
fi

echo "=== Predicting curves for P1=${TARGET_P1} cm (T3 mode: ${T3_MODE}) ==="
echo "    Training source: $PARENT_DIR"
echo "    Output dir:      ${OUTPUT_DIR:-<parent>/plots}"
echo ""

module purge 2>/dev/null || true
module load ROOT/6.30.06-foss-2023a 2>/dev/null || {
    echo "WARNING: could not load ROOT module; assuming ROOT is already in PATH"
}

# Clean any stale ACLiC artefacts for this macro
rm -f "${PREDICT_MACRO%.C}_C"* 2>/dev/null || true

ROOT_CMD="KLong_predict_from_geometry.C+(\"${PARENT_DIR}\", ${TARGET_P1}, \"${T3_MODE}\", \"${OUTPUT_DIR}\")"
root -l -b -q "${SCRIPT_DIR}/${ROOT_CMD}"

echo ""
echo "=== Done ==="

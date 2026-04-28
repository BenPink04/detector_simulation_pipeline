#!/usr/bin/env bash
# ==============================================================================
# run_pizza_fit.sh
#
# Fits parametric models to the momentum resolution and acceptance curves for
# all detector configurations found in a given archive directory.
#
# USAGE:
#   ./run_pizza_fit.sh <archive_dir>
#
# EXAMPLE:
#   ./run_pizza_fit.sh /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/Diss_Draft_2_20260415
#   ./run_pizza_fit.sh /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/Upstream_pizzas_20260416
#
# OUTPUT:
#   <archive_dir>/plots/KLong_fit_resolution.png
#   <archive_dir>/plots/KLong_fit_acceptance.png
#
# The two ROOT macros are expected to live alongside this script.
# ==============================================================================

set -euo pipefail

# --------------------------------------------------------------------------
# Environment: load ROOT if available via module system
# --------------------------------------------------------------------------
if command -v module &>/dev/null; then
    module purge
    module load ROOT/6.30.06-foss-2023a
fi

# --------------------------------------------------------------------------
# Argument check
# --------------------------------------------------------------------------
if [[ $# -lt 1 ]]; then
    echo ""
    echo "  ERROR: No archive directory supplied."
    echo ""
    echo "  Usage:  ./run_pizza_fit.sh <archive_dir>"
    echo ""
    echo "  Example:"
    echo "    ./run_pizza_fit.sh /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/Diss_Draft_2_20260415"
    echo ""
    exit 1
fi

ARCHIVE_DIR="$1"

# Strip trailing slash for cleanliness
ARCHIVE_DIR="${ARCHIVE_DIR%/}"

# --------------------------------------------------------------------------
# Validate directory
# --------------------------------------------------------------------------
if [[ ! -d "$ARCHIVE_DIR" ]]; then
    echo "ERROR: Directory not found: $ARCHIVE_DIR"
    exit 1
fi

# Quick check that some expected files exist
N_VEC=$(find "$ARCHIVE_DIR" -maxdepth 1 -name '*combined_vectors.root' 2>/dev/null | wc -l)
N_ACC=$(find "$ARCHIVE_DIR" -maxdepth 1 -name '*combined_acceptance.root' 2>/dev/null | wc -l)

if [[ $N_VEC -eq 0 && $N_ACC -eq 0 ]]; then
    echo "ERROR: No *combined_vectors.root or *combined_acceptance.root files found in:"
    echo "  $ARCHIVE_DIR"
    exit 1
fi

echo ""
echo "=============================================="
echo "  Pizza Position Fit"
echo "  Archive : $ARCHIVE_DIR"
echo "  Vectors files  : $N_VEC"
echo "  Acceptance files: $N_ACC"
echo "=============================================="

# --------------------------------------------------------------------------
# Locate macros (same directory as this script)
# --------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RES_MACRO="$SCRIPT_DIR/KLong_fit_resolution.C"
ACC_MACRO="$SCRIPT_DIR/KLong_fit_acceptance.C"

for macro in "$RES_MACRO" "$ACC_MACRO"; do
    if [[ ! -f "$macro" ]]; then
        echo "ERROR: Macro not found: $macro"
        exit 1
    fi
done

# --------------------------------------------------------------------------
# Create output plot directory
# --------------------------------------------------------------------------
mkdir -p "$ARCHIVE_DIR/plots"

# --------------------------------------------------------------------------
# Run resolution fit
# --------------------------------------------------------------------------
if [[ $N_VEC -gt 0 ]]; then
    echo ""
    echo "--- Running resolution fit ---"
    root -l -b -q "${RES_MACRO}+(\"${ARCHIVE_DIR}\")"
    echo ""
else
    echo "Skipping resolution fit (no combined_vectors.root files found)."
fi

# --------------------------------------------------------------------------
# Run acceptance fit
# --------------------------------------------------------------------------
if [[ $N_ACC -gt 0 ]]; then
    echo ""
    echo "--- Running acceptance fit ---"
    root -l -b -q "${ACC_MACRO}+(\"${ARCHIVE_DIR}\")"
    echo ""
else
    echo "Skipping acceptance fit (no combined_acceptance.root files found)."
fi

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "  Done.  Output plots:"
for f in "$ARCHIVE_DIR/plots/KLong_fit_resolution.png" \
          "$ARCHIVE_DIR/plots/KLong_fit_resolution.root" \
          "$ARCHIVE_DIR/plots/KLong_fit_acceptance.png" \
          "$ARCHIVE_DIR/plots/KLong_fit_acceptance.root"; do
    if [[ -f "$f" ]]; then
        echo "    [OK]  $f"
    fi
done
echo "=============================================="
echo ""

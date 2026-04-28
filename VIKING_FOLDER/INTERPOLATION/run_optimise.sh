#!/usr/bin/env bash
# =============================================================================
# run_optimise.sh
#
# Scans pizza position P1 over the interpolated range and identifies the
# "ideal" detector setup — maximising acceptance and minimising resolution —
# for upstream, downstream, or both T3 tracking modes.
#
# USAGE:
#   ./run_optimise.sh <archives_parent_dir> [t3_mode] [output_dir] [fom_mode]
#
# ARGUMENTS:
#   archives_parent_dir  — directory containing archive sub-dirs
#                          (each with *combined_vectors.root / *combined_acceptance.root)
#   t3_mode              — "upstream" | "downstream" | "both"  (default: both)
#   output_dir           — where to save outputs (default: <parent_dir>/plots)
#   fom_mode             — 1 = N/<FWHM>  (default)  |  2 = N/<FWHM^2>
#
# EXAMPLES:
#   ./run_optimise.sh /path/to/data
#   ./run_optimise.sh /path/to/data both /tmp/out 1
#   ./run_optimise.sh /path/to/data both /tmp/out 2
#
# OUTPUT:
#   <output_dir>/KLong_optimise_<mode>.png   — three-panel scan plot
#   <output_dir>/KLong_optimise_<mode>.root  — canvas saved as ROOT file
#   Console summary of optimal P1 per metric per mode
# =============================================================================

set -euo pipefail

PARENT_DIR="${1:-}"
T3_MODE="${2:-both}"
OUTPUT_DIR="${3:-}"
FOM_MODE="${4:-1}"

if [[ -z "$PARENT_DIR" ]]; then
    echo ""
    echo "  ERROR: No archive parent directory supplied."
    echo ""
    echo "  Usage: ./run_optimise.sh <archives_parent_dir> [t3_mode] [output_dir]"
    echo "  t3_mode options: upstream | downstream | both  (default: both)"
    echo ""
    echo "  Example:"
    echo "    ./run_optimise.sh /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS"
    echo ""
    exit 1
fi

if [[ ! -d "$PARENT_DIR" ]]; then
    echo "ERROR: Directory not found: $PARENT_DIR"
    exit 1
fi

if [[ "$T3_MODE" != "upstream" && "$T3_MODE" != "downstream" && "$T3_MODE" != "both" ]]; then
    echo "ERROR: t3_mode must be 'upstream', 'downstream', or 'both' (got: $T3_MODE)"
    exit 1
fi

if [[ "$FOM_MODE" != "1" && "$FOM_MODE" != "2" ]]; then
    echo "ERROR: fom_mode must be 1 or 2 (got: $FOM_MODE)"
    exit 1
fi

# Load ROOT
if command -v module &>/dev/null; then
    module purge
    module load ROOT/6.30.06-foss-2023a
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="${SCRIPT_DIR}/KLong_optimise_geometry.C"

if [[ ! -f "$MACRO" ]]; then
    echo "ERROR: Macro not found: $MACRO"
    exit 1
fi

N_VEC=$(find "$PARENT_DIR" -maxdepth 2 -name '*combined_vectors.root' 2>/dev/null | wc -l)
N_ACC=$(find "$PARENT_DIR" -maxdepth 2 -name '*combined_acceptance.root' 2>/dev/null | wc -l)

echo ""
echo "=============================================="
echo "  Geometry Optimisation Scan"
echo "  Source  : $PARENT_DIR"
echo "  Mode    : $T3_MODE"
echo "  Vectors : $N_VEC files"
echo "  Acceptance: $N_ACC files"
echo "=============================================="

if [[ $N_VEC -eq 0 && $N_ACC -eq 0 ]]; then
    echo "ERROR: No training files found under $PARENT_DIR"
    exit 1
fi

# Clean stale ACLiC artefacts
rm -f "${MACRO%.C}_C"* 2>/dev/null || true

root -l -b -q "${MACRO}+(\"${PARENT_DIR}\", \"${T3_MODE}\", \"${OUTPUT_DIR}\", ${FOM_MODE})"

echo ""
echo "=============================================="
OUT="${OUTPUT_DIR:-${PARENT_DIR}/plots}"
echo "  Done. Outputs written to: $OUT"
echo "=============================================="
echo ""

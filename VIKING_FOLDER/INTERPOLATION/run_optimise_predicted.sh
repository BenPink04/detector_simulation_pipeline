#!/usr/bin/env bash
# ============================================================================
# run_optimise_predicted.sh
#
# Wrapper for KLong_optimise_predicted.C
# Uses predicted TF1 curves (numerical integration) for metrics.
#
# Usage:
#   ./run_optimise_predicted.sh <archives_parent_dir> [upstream|downstream|both] [output_dir] [fom_mode]
#   fom_mode: 1 = N/<FWHM> (default), 2 = N/<FWHM^2>
# ============================================================================

set -euo pipefail

PARENT_DIR="${1:?Usage: $0 <archives_parent_dir> [upstream|downstream|both] [output_dir] [fom_mode]}"
MODE="${2:-both}"
OUT_DIR="${3:-}"
FOM_MODE="${4:-1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="${SCRIPT_DIR}/KLong_optimise_predicted.C"

# Load ROOT if available via modules
if command -v module &>/dev/null; then
    module load ROOT 2>/dev/null || true
fi

# Clean ACLiC artefacts so recompilation is forced
rm -f "${SCRIPT_DIR}/KLong_optimise_predicted_C"* 2>/dev/null || true

if [ -z "${OUT_DIR}" ]; then
    root -l -b -q "${MACRO}+(\"${PARENT_DIR}\",\"${MODE}\",\"\",${FOM_MODE})"
else
    root -l -b -q "${MACRO}+(\"${PARENT_DIR}\",\"${MODE}\",\"${OUT_DIR}\",${FOM_MODE})"
fi

#!/usr/bin/env bash
# generate_truez_files.sh
#
# Generates *_combined_vectors_truez.root from every existing
# *_combined_vectors.root in SIMULATION_RESULTS, without rerunning
# the pipeline.  Run this once after fixing the ROOTEOF heredoc bug.
#
# Usage:
#   cd HISTOGRAM_MAKING
#   ./generate_truez_files.sh
#
# The script reads reco_p_truez from the standard combined file and
# writes a second file where:
#   reco_p       = reco_p_truez  (truth-vertex reconstruction)
#   reco_p_poca  = reco_p        (original PoCA reconstruction, preserved)
#   reco_p_truet, true_p, reco_vertex_*, true_vertex_* — unchanged

set -euo pipefail

module load ROOT/6.30.06-foss-2023a 2>/dev/null || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="${SCRIPT_DIR}/make_truez.C"
SIM_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS"

OK=0
SKIP=0
FAIL=0

echo "============================================================"
echo "  generate_truez_files.sh"
echo "  Source: ${SIM_DIR}"
echo "============================================================"
echo ""

for CONFIG_DIR in "${SIM_DIR}"/*/; do
    [[ -d "${CONFIG_DIR}" ]] || continue
    CONFIG=$(basename "${CONFIG_DIR}")
    IN="${CONFIG_DIR}${CONFIG}_combined_vectors.root"
    OUT="${CONFIG_DIR}${CONFIG}_combined_vectors_truez.root"

    if [[ ! -f "${IN}" ]]; then
        echo "  [SKIP] ${CONFIG} — no combined_vectors.root"
        SKIP=$((SKIP+1))
        continue
    fi

    echo "  Processing: ${CONFIG}"

    root -l -b -q "${MACRO}(\"${IN}\",\"${OUT}\")" 2>&1 | grep -v "^$\|^Processing\|^Info\|ROOT v"

    if [[ -f "${OUT}" ]]; then
        echo "    -> $(du -h "${OUT}" | cut -f1)  ${CONFIG}_combined_vectors_truez.root"
        OK=$((OK+1))
    else
        echo "  [FAIL] No output produced for ${CONFIG}" >&2
        FAIL=$((FAIL+1))
    fi
    echo ""
done

echo "============================================================"
echo "  Done."
echo "  Generated : ${OK}"
echo "  Skipped   : ${SKIP}  (no combined_vectors.root)"
echo "  Failed    : ${FAIL}"
echo "============================================================"

#!/usr/bin/env bash
# KLong_run_all.sh
#
# Runs all non-by-pizza histogram and validation macros, substituting the
# archive folder name at runtime, and collects every output PNG into a
# single directory.
#
# Usage:
#   ./KLong_run_all.sh [archive_name] [output_dir]
#
# Arguments:
#   archive_name   Name of the folder inside ARCHIVED_RESULTS to read from.
#                  Default: TGRAPH_TEST_20260316
#   output_dir     Directory to write all output PNGs into.
#                  Default: <archive_name>_plots  (created if it doesn't exist)
#
# Example:
#   ./KLong_run_all.sh TGRAPH_TEST_20260317 my_plots
#
# Macros run (in order):
#   KLong_plot_compare_acceptance.C       -> KLong_acceptance_comparison.png
#   KLong_plot_compare_resolution.C       -> KLong_resolution_comparison.png
#   KLong_plot_compare_spreads.C          -> KLong_plot_compare_spreads.png
#   KLong_plot_gaussian_spreads.C         -> KLong_plot_gaussian_spreads.png
#   KLong_plot_gaussian_spreads_truncated.C -> KLong_plot_gaussian_spreads_truncated.png
#   KLong_plot_resolution_histbar.C       -> KLong_histbar_plot_*.png
#   KLong_validate_reconstruction.C       -> KLong_validate_ratio.png
#                                            KLong_validate_delta_p.png
#   KLong_validate_vertex.C               -> KLong_validate_vertex_residuals.png
#                                            KLong_validate_vertex_z_profile.png
#
# Note: KLong_plot_acceptance.C reads from SIMULATION_RESULTS (not the archive)
# and is intentionally excluded — run it separately if needed.

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
ARCHIVE="${1:-TGRAPH_TEST_20260316}"
OUTDIR="${2:-${ARCHIVE}_plots}"

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR_ABS="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"

# Verify the archive folder actually exists before spending time running macros
ARCHIVE_BASE="/users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS"
if [[ ! -d "${ARCHIVE_BASE}/${ARCHIVE}" ]]; then
    echo "ERROR: Archive folder not found: ${ARCHIVE_BASE}/${ARCHIVE}" >&2
    exit 1
fi

echo "============================================================"
echo "  KLong_run_all.sh"
echo "  Archive : ${ARCHIVE}"
echo "  Output  : ${OUTDIR_ABS}"
echo "============================================================"
echo ""

# Count how many ROOT files are in the archive
N_CONFIGS=$(find "${ARCHIVE_BASE}/${ARCHIVE}" -name '*combined_vectors.root' | wc -l)
echo "  Found ${N_CONFIGS} configuration(s) in archive"
echo ""

# ── Helper: sed-substitute archive name into a temp file, run it, move PNGs ──
#
# run_macro <source_macro.C> <png_glob_pattern> [<png_glob_pattern> ...]
#
# Each png_glob_pattern may contain shell wildcards (e.g. "KLong_histbar_plot_*.png").
# Patterns are evaluated in the SCRIPT_DIR after root finishes.

run_macro() {
    local macro="$1"; shift
    local pngs=("$@")     # remaining args = PNG glob patterns to collect

    local src="${SCRIPT_DIR}/${macro}"

    # Put the temp file in a tmpdir but keep the original filename so ROOT's
    # function-name-matches-filename auto-call mechanism works correctly.
    local tmpdir
    tmpdir="$(mktemp -d "${SCRIPT_DIR}/tmp_XXXXXX")"
    local tmp="${tmpdir}/${macro}"

    # Substitute the default archive name with the requested one
    sed "s|TGRAPH_TEST_20260316|${ARCHIVE}|g" "${src}" > "${tmp}"

    echo "──────────────────────────────────────────────────────────"
    echo "  Running: ${macro}"

    # Run from SCRIPT_DIR so relative SaveAs paths land there for easy collection
    (cd "${SCRIPT_DIR}" && root -l -b -q "${tmp}" 2>&1) || {
        echo "  WARNING: ${macro} exited with non-zero status" >&2
    }

    # Move produced PNGs to output directory
    local moved=0
    for pattern in "${pngs[@]}"; do
        for f in "${SCRIPT_DIR}"/${pattern}; do
            if [[ -f "${f}" ]]; then
                mv -f "${f}" "${OUTDIR_ABS}/"
                echo "  -> $(basename "${f}")"
                (( moved++ )) || true
            fi
        done
    done
    [[ ${moved} -eq 0 ]] && echo "  WARNING: no PNG output found for ${macro}" >&2

    rm -rf "${tmpdir}"
    echo ""
}

# ── Run each macro ─────────────────────────────────────────────────────────────

run_macro "KLong_plot_compare_acceptance.C" \
    "KLong_acceptance_comparison.png"

run_macro "KLong_plot_compare_resolution.C" \
    "KLong_resolution_comparison.png"

run_macro "KLong_plot_compare_spreads.C" \
    "KLong_plot_compare_spreads.png"

run_macro "KLong_plot_gaussian_spreads.C" \
    "KLong_plot_gaussian_spreads.png"

run_macro "KLong_plot_gaussian_spreads_truncated.C" \
    "KLong_plot_gaussian_spreads_truncated.png"

run_macro "KLong_plot_resolution_histbar.C" \
    "KLong_histbar_plot_*.png"

run_macro "KLong_validate_reconstruction.C" \
    "KLong_validate_ratio.png" \
    "KLong_validate_delta_p.png"

run_macro "KLong_validate_vertex.C" \
    "KLong_validate_vertex_residuals.png" \
    "KLong_validate_vertex_z_profile.png"

# ── Summary ────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "  Done. Output PNGs in: ${OUTDIR_ABS}"
echo ""
ls -1 "${OUTDIR_ABS}"/*.png 2>/dev/null | while read -r f; do echo "  $(basename "${f}")"; done
echo "============================================================"

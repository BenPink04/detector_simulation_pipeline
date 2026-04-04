#!/usr/bin/env bash
# KLong_run_all.sh
#
# Runs all non-by-pizza histogram and validation macros, substituting the
# archive folder name at runtime, and collects every output PNG into a
# single directory.
#
# Two passes are run automatically:
#   Pass 1 (standard): macros read *_combined_vectors.root       -> output_dir/standard/
#   Pass 2 (truez):    macros read *_combined_vectors_truez.root -> output_dir/truez/
# Acceptance macros (which read *_combined_acceptance.root) are run in both
# passes; their output is identical but placed in the respective subdirectory.
#
# Usage:
#   ./KLong_run_all.sh [archive_name] [output_dir]
#
# Arguments:
#   archive_name        Name of the folder inside ARCHIVED_RESULTS to read from.
#                       Default: TGRAPH_TEST_20260316
#   output_dir          Root directory for output PNGs.
#                       Default: <archive_name>_plots  (created if it doesn't exist)
#                       Subdirs standard/ and truez/ are created automatically.
#   anomaly_threshold   Maximum allowed |delta_p / true_p| before an event is
#                       filtered out of all histogram macros.
#                       Default: 1.0  (100% — removes only wildly unphysical events)
#                       Set lower (e.g. 0.5) to focus plots on the core population.
#
# Example:
#   ./KLong_run_all.sh TGRAPH_TEST_20260317 my_plots
#   ./KLong_run_all.sh TGRAPH_TEST_20260317 my_plots 0.5   # tighter 50% cut
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
ANOMALY_THRESHOLD="${3:-1.0}"  # fractional |delta_p/true_p| cut applied in all macros

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${OUTDIR}/standard" "${OUTDIR}/truez"
OUTDIR_ROOT="$(cd "${OUTDIR}" && pwd)"
OUTDIR_ABS="${OUTDIR_ROOT}/standard"   # overridden per-pass

# Verify the archive folder actually exists before spending time running macros
ARCHIVE_BASE="/users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS"
if [[ ! -d "${ARCHIVE_BASE}/${ARCHIVE}" ]]; then
    echo "ERROR: Archive folder not found: ${ARCHIVE_BASE}/${ARCHIVE}" >&2
    exit 1
fi

echo "============================================================"
echo "  KLong_run_all.sh"
echo "  Archive            : ${ARCHIVE}"
echo "  Output             : ${OUTDIR_ROOT}"
echo "    standard/ -> PoCA reconstruction"
echo "    truez/    -> truth-vertex reconstruction"
echo "  Anomaly threshold  : ${ANOMALY_THRESHOLD}  (|delta_p/true_p| cut)"
echo "============================================================"
echo ""

# Count ROOT files in archive
N_CONFIGS=$(find "${ARCHIVE_BASE}/${ARCHIVE}" -name '*combined_vectors.root' ! -name '*_truez*' | wc -l)
N_TRUEZ=$(find "${ARCHIVE_BASE}/${ARCHIVE}" -name '*combined_vectors_truez.root' | wc -l)
echo "  Found ${N_CONFIGS} standard configuration(s) in archive"
echo "  Found ${N_TRUEZ} truez configuration(s) in archive"
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

    # Substitute the default archive name, vector filename (truez pass),
    # and anomaly_threshold value so all macros use the same cut.
    if [[ "${VECTORS_MODE}" == "truez" ]]; then
        sed -e "s|TGRAPH_TEST_20260316|${ARCHIVE}|g" \
            -e "s|_combined_vectors\.root|_combined_vectors_truez.root|g" \
            -e "s|anomaly_threshold = [0-9.]*|anomaly_threshold = ${ANOMALY_THRESHOLD}|g" \
            "${src}" > "${tmp}"
    else
        sed -e "s|TGRAPH_TEST_20260316|${ARCHIVE}|g" \
            -e "s|anomaly_threshold = [0-9.]*|anomaly_threshold = ${ANOMALY_THRESHOLD}|g" \
            "${src}" > "${tmp}"
    fi

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

# ── Inner helper: runs all macros once (uses globals OUTDIR_ABS, VECTORS_MODE) ─
run_all_macros() {
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
}

# ── Pass 1: standard PoCA reconstruction ──────────────────────────────────────
VECTORS_MODE=standard
OUTDIR_ABS="${OUTDIR_ROOT}/standard"
echo "============================================================"
echo "  Pass 1/2 — Standard PoCA reconstruction"
echo "  Output  : ${OUTDIR_ABS}"
echo "============================================================"
echo ""
run_all_macros

# ── Pass 2: truth-vertex reconstruction ───────────────────────────────────────
VECTORS_MODE=truez
OUTDIR_ABS="${OUTDIR_ROOT}/truez"
echo "============================================================"
echo "  Pass 2/2 — Truth-vertex reconstruction (reco_p_truez)"
echo "  Output  : ${OUTDIR_ABS}"
echo "============================================================"
if [[ "${N_TRUEZ}" -eq 0 ]]; then
    echo "  WARNING: No *_combined_vectors_truez.root files found in archive." >&2
    echo "           Run a new simulation with the updated KLong_save_vectors.C first." >&2
    echo ""
else
    echo ""
    run_all_macros
fi

# ── Summary ────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "  Done."
echo ""
echo "  Standard plots : ${OUTDIR_ROOT}/standard/"
ls -1 "${OUTDIR_ROOT}/standard/"*.png 2>/dev/null | wc -l | xargs -I{} echo "    {} PNG(s)"
echo "  TrueZ plots    : ${OUTDIR_ROOT}/truez/"
ls -1 "${OUTDIR_ROOT}/truez/"*.png 2>/dev/null | wc -l | xargs -I{} echo "    {} PNG(s)"
echo "============================================================"

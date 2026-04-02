#!/usr/bin/env bash

# ============================================================
#   File_Archiver.sh
#   Copies combined_acceptance and combined_vectors root files
#   from SIMULATION_RESULTS into a timestamped archive folder.
#
#   Usage:
#     ./File_Archiver.sh <ARCHIVE_NAME>
#     ./File_Archiver.sh          # prompts for name interactively
#
#   Arguments:
#     ARCHIVE_NAME   Base name for the archive directory.
#                    A date suffix (_YYYYMMDD) is appended automatically.
#                    Example: TGRAPH_TEST_5  ->  TGRAPH_TEST_5_20260402
#
#   What it does:
#     For every subdirectory in SIMULATION_RESULTS/ it looks for:
#       <CONFIG>_combined_acceptance.root
#       <CONFIG>_combined_vectors.root
#     and copies both into ARCHIVED_RESULTS/<ARCHIVE_NAME>_<DATE>/.
#     Missing files are flagged as [MISSING] and counted in the summary.
#
#   Example:
#     cd HISTOGRAM_MAKING
#     ./File_Archiver.sh TGRAPH_TEST_5
#     # -> ARCHIVED_RESULTS/TGRAPH_TEST_5_20260402/
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIKING_DIR="$(dirname "$SCRIPT_DIR")"
SIM_DIR="$VIKING_DIR/SIMULATION_RESULTS"
ARCHIVE_BASE="$VIKING_DIR/ARCHIVED_RESULTS"

# ------------------------------------------------------------
#   Get archive name from user
# ------------------------------------------------------------
if [[ -n "$1" ]]; then
    BASE_NAME="$1"
else
    read -rp "Enter archive name: " BASE_NAME
fi

if [[ -z "$BASE_NAME" ]]; then
    echo "ERROR: No archive name provided. Exiting."
    exit 1
fi

DATE_SUFFIX=$(date +%Y%m%d)
ARCHIVE_NAME="${BASE_NAME}_${DATE_SUFFIX}"
ARCHIVE_DIR="$ARCHIVE_BASE/$ARCHIVE_NAME"

# ------------------------------------------------------------
#   Safety checks
# ------------------------------------------------------------
if [[ ! -d "$SIM_DIR" ]]; then
    echo "ERROR: Simulation results directory not found:"
    echo "       $SIM_DIR"
    exit 1
fi

if [[ -d "$ARCHIVE_DIR" ]]; then
    echo "WARNING: Archive directory already exists:"
    echo "         $ARCHIVE_DIR"
    read -rp "Overwrite / add to it? [y/N]: " CONFIRM
    [[ "$CONFIRM" =~ ^[Yy]$ ]] || { echo "Aborted."; exit 0; }
fi

mkdir -p "$ARCHIVE_DIR"

# ------------------------------------------------------------
#   Find and copy combined files
# ------------------------------------------------------------
echo ""
echo "============================================================"
echo "  File_Archiver.sh"
echo "  Archive : $ARCHIVE_NAME"
echo "  Source  : $SIM_DIR"
echo "  Dest    : $ARCHIVE_DIR"
echo "============================================================"
echo ""

CONFIGS_FOUND=0
CONFIGS_OK=0
CONFIGS_WARN=0

for CONFIG_DIR in "$SIM_DIR"/*/; do
    [[ -d "$CONFIG_DIR" ]] || continue

    CONFIG_NAME=$(basename "$CONFIG_DIR")
    CONFIGS_FOUND=$((CONFIGS_FOUND + 1))

    ACCEPTANCE_FILE="$CONFIG_DIR/${CONFIG_NAME}_combined_acceptance.root"
    VECTORS_FILE="$CONFIG_DIR/${CONFIG_NAME}_combined_vectors.root"

    echo "  Processing: $CONFIG_NAME"

    MISSING=0

    if [[ -f "$ACCEPTANCE_FILE" ]]; then
        cp "$ACCEPTANCE_FILE" "$ARCHIVE_DIR/"
        echo "    [OK] Copied combined_acceptance.root"
    else
        echo "    [MISSING] No combined_acceptance.root found"
        MISSING=$((MISSING + 1))
    fi

    if [[ -f "$VECTORS_FILE" ]]; then
        cp "$VECTORS_FILE" "$ARCHIVE_DIR/"
        echo "    [OK] Copied combined_vectors.root"
    else
        echo "    [MISSING] No combined_vectors.root found"
        MISSING=$((MISSING + 1))
    fi

    if [[ "$MISSING" -eq 0 ]]; then
        CONFIGS_OK=$((CONFIGS_OK + 1))
    else
        CONFIGS_WARN=$((CONFIGS_WARN + 1))
    fi

    echo ""
done

# ------------------------------------------------------------
#   Summary
# ------------------------------------------------------------
echo "============================================================"
echo "  Done."
echo ""
echo "  Configurations found   : $CONFIGS_FOUND"
echo "  Fully archived         : $CONFIGS_OK"
echo "  Missing files          : $CONFIGS_WARN"
echo ""
echo "  Output: $ARCHIVE_DIR"
echo "============================================================"

# Warn if archive is empty
if [[ -z "$(ls -A "$ARCHIVE_DIR")" ]]; then
    echo ""
    echo "  WARNING: Archive directory is empty. No combined files"
    echo "           were found. Has the full pipeline been run?"
    echo ""
    exit 1
fi

exit 0

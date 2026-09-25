#!/bin/bash

# KS Regeneration Position Sweep
#
# Runs ks_regeneration_simulation.sh for a set of KS beam start positions
# spanning the region from the dipole exit (z=200 cm) to the first detector (z=400 cm).
#
# Default: 5 positions at 200, 250, 300, 350, 400 cm.
# Override by passing positions as arguments:
#   ./ks_regeneration_sweep.sh 210 260 310 360

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SIM_SCRIPT="${SCRIPT_DIR}/ks_regeneration_simulation.sh"

if [ ! -f "$SIM_SCRIPT" ]; then
    echo "ERROR: ks_regeneration_simulation.sh not found at ${SIM_SCRIPT}"
    exit 1
fi

chmod +x "$SIM_SCRIPT"

# Default sweep positions (cm); override via command-line arguments
if [ $# -gt 0 ]; then
    POSITIONS=("$@")
else
    POSITIONS=(200 250 300 350 400)
fi

echo "=== KS Regeneration Position Sweep ==="
echo "Positions to sweep (cm): ${POSITIONS[*]}"
echo ""

for z0 in "${POSITIONS[@]}"; do
    echo "────────────────────────────────────────"
    echo "Setting up KS pipeline for z0 = ${z0} cm"
    echo "────────────────────────────────────────"
    "$SIM_SCRIPT" "$z0"
    echo ""
done

echo "=== Sweep setup complete ==="
echo ""
echo "A separate submission script was created for each z-position."
echo "Submit all at once with:"
echo ""
for z0 in "${POSITIONS[@]}"; do
    echo "  /users/bp969/scratch/JOBSCRIPTS_TESTS/submit_KS_Z0-${z0}_T1-425_T2-435_T3-465_T4-475_P1-400_P2-415_F1-455_F2-465_E1-700.sh"
done
echo ""
echo "Or submit all in one go:"
echo "  for f in /users/bp969/scratch/JOBSCRIPTS_TESTS/submit_KS_Z0-*.sh; do bash \"\$f\"; done"

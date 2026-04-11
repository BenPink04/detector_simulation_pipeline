#!/usr/bin/env bash
# reparse_all_configs.sh
#
# Re-submits parsing (and combination) jobs for all 4 detector configurations
# using the already-fixed KLong_save_vectors.C / KLong_save_momentum_acceptance.C
# scripts that live in each config's DATA_READ directory.
#
# This script does NOT re-run Geant4 simulations.  It overwrites only the
# per-file _vectors.root and _acceptance.root outputs, then re-combines them
# into _combined_vectors.root, _combined_acceptance.root and the truez
# diagnostic file.
#
# Usage:
#   bash reparse_all_configs.sh [--dry-run]
#
# With --dry-run the job scripts are written but not submitted.
#
# Prerequisites:
#   - Raw simulation ROOT files already in SIMULATION_RESULTS/<config>/
#   - Fixed KLong_save_vectors.C + KLong_save_momentum_acceptance.C already in
#     DATA_PARSING/<config>_DATA_READ/ (no pip_tof_z = z_tof remaining)
#   - parse_diag_events.py present in the same DATA_READ directory

set -euo pipefail

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

BASE="/users/bp969/scratch/VIKING_FOLDER"
SCRATCH="/users/bp969/scratch"
TEMP_DIR="${SCRATCH}/JOBSCRIPTS_TESTS/temp_jobs"
mkdir -p "$TEMP_DIR"

# ── The 4 configurations ──────────────────────────────────────────────────────
CONFIGS=(
    "T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700"
    "T1-320_T2-330_T3-680_T4-690_P1-295_P2-310_F1-340_F2-350_E1-700"
    "T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700"
    "T1-480_T2-490_T3-680_T4-690_P1-455_P2-470_F1-500_F2-510_E1-700"
)

# ── Helper: verify the fixed scripts are in place ─────────────────────────────
check_scripts() {
    local cfg="$1"
    local pdir="${BASE}/DATA_PARSING/${cfg}_DATA_READ"
    for f in KLong_save_vectors.C KLong_save_momentum_acceptance.C; do
        if [[ ! -f "${pdir}/${f}" ]]; then
            echo "ERROR: ${pdir}/${f} not found!" >&2
            return 1
        fi
    done
    # Check that the geometry-face bug is gone
    if grep -q "pip_tof_z = z_tof" "${pdir}/KLong_save_vectors.C" 2>/dev/null; then
        echo "ERROR: ${pdir}/KLong_save_vectors.C still has pip_tof_z = z_tof!" >&2
        return 1
    fi
    # Ensure parse_diag_events.py is present (copy from HISTOGRAM_MAKING if not)
    if [[ ! -f "${pdir}/parse_diag_events.py" ]]; then
        echo "Copying parse_diag_events.py to ${pdir}/"
        cp "${BASE}/HISTOGRAM_MAKING/parse_diag_events.py" "${pdir}/"
    fi
    return 0
}

# ── Create a parsing job for one file ────────────────────────────────────────
create_parse_job() {
    local cfg="$1"
    local fnum="$2"
    local pdir="${BASE}/DATA_PARSING/${cfg}_DATA_READ"
    local rdir="${BASE}/SIMULATION_RESULTS/${cfg}"
    local logdir="${SCRATCH}/JOB_LOGS/${cfg}/parsing"
    mkdir -p "$logdir"
    local jobscript="${TEMP_DIR}/reparse_${cfg}_F${fnum}.job"

    cat > "$jobscript" << JOBEOF
#!/usr/bin/env bash
#SBATCH --job-name=Reparse_${cfg}_F${fnum}
#SBATCH --partition=nodes
#SBATCH --time=0-04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${logdir}/Reparse_${cfg}_F${fnum}_%j.log
#SBATCH --error=${logdir}/Reparse_${cfg}_F${fnum}_%j.err

set -e

echo "=== Reparse job F${fnum} for ${cfg} ==="
module purge
module load ROOT/6.30.06-foss-2023a

cd "${pdir}"

INPUT_FILE="${rdir}/${cfg}_${fnum}.root"
echo "Waiting for: \$INPUT_FILE"
while [[ ! -f "\$INPUT_FILE" ]]; do
    echo "  ...sleeping 30s"
    sleep 30
done

INPUT_FILENAME="${cfg}_${fnum}.root"
cp "\$INPUT_FILE" "\$INPUT_FILENAME"

ACCEPTANCE_OUT="${rdir}/${cfg}_${fnum}_acceptance.root"
VECTORS_OUT="${rdir}/${cfg}_${fnum}_vectors.root"

echo "Running acceptance macro..."
root -l -b -q "KLong_save_momentum_acceptance.C(\"\${INPUT_FILENAME}\")"

echo "Running vectors macro..."
DIAG_LOG="${rdir}/${cfg}_${fnum}_diag_clean.txt"
root -l -b -q "KLong_save_vectors.C(\"\${INPUT_FILENAME}\")" 2>&1 | \
    python3 "${pdir}/parse_diag_events.py" --log "\$DIAG_LOG"
ROOT_EXIT=\${PIPESTATUS[0]}
if [[ \$ROOT_EXIT -ne 0 ]]; then
    echo "ERROR: KLong_save_vectors.C failed (exit \$ROOT_EXIT)"
    exit \$ROOT_EXIT
fi

mv -f "${cfg}_${fnum}_acceptance.root" "\$ACCEPTANCE_OUT" 2>/dev/null || true
mv -f "${cfg}_${fnum}_vectors.root"    "\$VECTORS_OUT"    2>/dev/null || true
rm -f "\$INPUT_FILENAME"

[[ -f "\$ACCEPTANCE_OUT" ]] || { echo "ERROR: missing \$ACCEPTANCE_OUT"; exit 1; }
[[ -f "\$VECTORS_OUT"    ]] || { echo "ERROR: missing \$VECTORS_OUT";    exit 1; }

echo "=== Done F${fnum} ==="
JOBEOF
    echo "$jobscript"
}

# ── Create the combination job for one config ─────────────────────────────────
create_combine_job() {
    local cfg="$1"
    local pdir="${BASE}/DATA_PARSING/${cfg}_DATA_READ"
    local rdir="${BASE}/SIMULATION_RESULTS/${cfg}"
    local logdir="${SCRATCH}/JOB_LOGS/${cfg}/combination"
    mkdir -p "$logdir"
    local jobscript="${TEMP_DIR}/recombine_${cfg}.job"

    cat > "$jobscript" << JOBEOF
#!/usr/bin/env bash
#SBATCH --job-name=Recombine_${cfg}
#SBATCH --partition=nodes
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=${logdir}/Recombine_${cfg}_%j.log
#SBATCH --error=${logdir}/Recombine_${cfg}_%j.err

set -e

echo "=== Combination job for ${cfg} ==="
module purge
module load ROOT/6.30.06-foss-2023a

cd "${pdir}"

# Verify all 100 per-file outputs exist
echo "Checking acceptance files..."
for i in \$(seq 1 100); do
    f="${rdir}/${cfg}_\${i}_acceptance.root"
    while [[ ! -f "\$f" ]]; do echo "waiting \$f..."; sleep 30; done
done
echo "Checking vectors files..."
for i in \$(seq 1 100); do
    f="${rdir}/${cfg}_\${i}_vectors.root"
    while [[ ! -f "\$f" ]]; do echo "waiting \$f..."; sleep 30; done
done

echo "Combining acceptance files..."
hadd -f "${rdir}/${cfg}_combined_acceptance.root" \
    ${rdir}/${cfg}_[0-9]*_acceptance.root

echo "Combining vectors files..."
hadd -f "${rdir}/${cfg}_combined_vectors.root" \
    ${rdir}/${cfg}_[0-9]*_vectors.root

echo "Creating truth-vertex diagnostic combined file..."
root -l -b -q << 'ROOTEOF'
{
    TFile *fin = TFile::Open("${rdir}/${cfg}_combined_vectors.root");
    if (!fin || fin->IsZombie()) { std::cerr << "ERROR open\n"; return; }
    TTree *tin = (TTree*)fin->Get("kaonVectors");
    if (!tin) { std::cerr << "ERROR: no kaonVectors tree\n"; return; }

    TFile *fout = new TFile(
        "${rdir}/${cfg}_combined_vectors_truez.root", "RECREATE");
    TTree *tout = new TTree("kaonVectors",
        "Kaon vectors - truth-vertex diagnostic (reco_p = reco_p_truez)");

    std::vector<double> *reco_p=nullptr, *true_p=nullptr;
    std::vector<double> *reco_p_truez=nullptr, *reco_p_truet=nullptr;
    std::vector<double> *reco_vx=nullptr, *reco_vy=nullptr, *reco_vz=nullptr;
    std::vector<double> *true_vx=nullptr, *true_vy=nullptr, *true_vz=nullptr;

    tin->SetBranchAddress("reco_p",        &reco_p);
    tin->SetBranchAddress("true_p",        &true_p);
    tin->SetBranchAddress("reco_p_truez",  &reco_p_truez);
    tin->SetBranchAddress("reco_p_truet",  &reco_p_truet);
    tin->SetBranchAddress("reco_vertex_x", &reco_vx);
    tin->SetBranchAddress("reco_vertex_y", &reco_vy);
    tin->SetBranchAddress("reco_vertex_z", &reco_vz);
    tin->SetBranchAddress("true_vertex_x", &true_vx);
    tin->SetBranchAddress("true_vertex_y", &true_vy);
    tin->SetBranchAddress("true_vertex_z", &true_vz);

    std::vector<double> o_rp, o_rp_poca, o_rp_truet, o_tp;
    std::vector<double> o_rvx, o_rvy, o_rvz, o_tvx, o_tvy, o_tvz;
    tout->Branch("reco_p",        &o_rp);
    tout->Branch("reco_p_poca",   &o_rp_poca);
    tout->Branch("reco_p_truet",  &o_rp_truet);
    tout->Branch("true_p",        &o_tp);
    tout->Branch("reco_vertex_x", &o_rvx);
    tout->Branch("reco_vertex_y", &o_rvy);
    tout->Branch("reco_vertex_z", &o_rvz);
    tout->Branch("true_vertex_x", &o_tvx);
    tout->Branch("true_vertex_y", &o_tvy);
    tout->Branch("true_vertex_z", &o_tvz);

    Long64_t nEntries = tin->GetEntries();
    for (Long64_t ev = 0; ev < nEntries; ++ev) {
        tin->GetEntry(ev);
        o_rp       = *reco_p_truez;
        o_rp_poca  = *reco_p;
        o_rp_truet = *reco_p_truet;
        o_tp       = *true_p;
        o_rvx = *reco_vx; o_rvy = *reco_vy; o_rvz = *reco_vz;
        o_tvx = *true_vx; o_tvy = *true_vy; o_tvz = *true_vz;
        tout->Fill();
    }
    fout->Write(); fout->Close(); fin->Close();
    std::cout << "truez file written (" << nEntries << " entries)\n";
}
ROOTEOF

echo "=== Combination complete for ${cfg} ==="
JOBEOF
    echo "$jobscript"
}

# ── Main loop — process each config ──────────────────────────────────────────
ALL_PARSE_IDS=()

for CFG in "${CONFIGS[@]}"; do
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Config: $CFG"

    if ! check_scripts "$CFG"; then
        echo "SKIPPING $CFG due to script check failure" >&2
        continue
    fi

    # Determine which file numbers have raw sim files
    RDIR="${BASE}/SIMULATION_RESULTS/${CFG}"
    PARSE_JOB_IDS=()

    for FNUM in $(seq 1 100); do
        RAW_FILE="${RDIR}/${CFG}_${FNUM}.root"
        if [[ ! -f "$RAW_FILE" ]]; then
            echo "  [SKIP F${FNUM}] raw sim file not found: $RAW_FILE"
            continue
        fi
        JOBSCRIPT=$(create_parse_job "$CFG" "$FNUM")
        if $DRY_RUN; then
            echo "  [DRY-RUN] would submit: $JOBSCRIPT"
        else
            JOB_ID=$(sbatch "$JOBSCRIPT" | awk '{print $4}')
            echo "  Submitted F${FNUM} → job $JOB_ID"
            PARSE_JOB_IDS+=("$JOB_ID")
            ALL_PARSE_IDS+=("$JOB_ID")
        fi
    done

    # Create & submit combination job, dependent on all parse jobs for this config
    COMBINE_SCRIPT=$(create_combine_job "$CFG")
    if $DRY_RUN; then
        echo "  [DRY-RUN] would submit combine: $COMBINE_SCRIPT"
    elif [[ ${#PARSE_JOB_IDS[@]} -gt 0 ]]; then
        DEP=$(IFS=:; echo "afterok:${PARSE_JOB_IDS[*]}")
        COMBINE_ID=$(sbatch --dependency="$DEP" "$COMBINE_SCRIPT" | awk '{print $4}')
        echo "  Submitted combination → job $COMBINE_ID (dep: $DEP)"
    else
        echo "  No parse jobs submitted for $CFG — skipping combine"
    fi
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Done.  Total parse jobs submitted: ${#ALL_PARSE_IDS[@]}"
$DRY_RUN && echo "(DRY-RUN mode — no jobs were actually submitted)"

#!/bin/bash

# Overnight Complete Pipeline Script
# This script runs the full simulation pipeline with automatic histogram generation
# Usage: ./overnight_pipeline.sh T1 T2 T3 T4 P1 P2 F1 F2 E1

set -e

# Check if correct number of arguments provided
if [ $# -ne 9 ]; then
    echo "Usage: $0 T1 T2 T3 T4 P1 P2 F1 F2 E1"
    echo "Example: $0 250 260 330 340 215 230 290 300 400"
    exit 1
fi

# Input parameters
T1=$1 T2=$2 T3=$3 T4=$4 P1=$5 P2=$6 F1=$7 F2=$8 E1=$9
CONFIG_STR="T1-${T1}_T2-${T2}_T3-${T3}_T4-${T4}_P1-${P1}_P2-${P2}_F1-${F1}_F2-${F2}_E1-${E1}"
RESULTS_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/${CONFIG_STR}"

echo "=== OVERNIGHT PIPELINE SETUP ==="
echo "Configuration: $CONFIG_STR"
echo "Setting up complete pipeline with automatic histogram generation..."

# Step 1: Run the main pipeline setup
echo "=== Step 1: Setting up main simulation pipeline ==="
/users/bp969/scratch/JOBSCRIPTS_TESTS/detector_simulation_master.sh $T1 $T2 $T3 $T4 $P1 $P2 $F1 $F2 $E1

# Step 2: Create final histogram generation job
echo "=== Step 2: Creating automatic histogram generation job ==="

cat > "/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/histograms_${CONFIG_STR}.job" << EOF
#!/usr/bin/env bash
#SBATCH --job-name=Histograms_${CONFIG_STR}
#SBATCH --partition=nodes
#SBATCH --time=0-00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --account=pet-hadron-2019
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bp969@york.ac.uk
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

set -e

module purge
module load ROOT/6.30.06-foss-2023a

echo "=== AUTOMATIC HISTOGRAM GENERATION ==="
echo "Configuration: ${CONFIG_STR}"
echo "Results directory: ${RESULTS_DIR}"

cd /users/bp969/scratch/VIKING_FOLDER/HISTOGRAM_MAKING

# Check for required combined files
ACCEPTANCE_FILE="${RESULTS_DIR}/${CONFIG_STR}_combined_acceptance.root"
VECTORS_FILE="${RESULTS_DIR}/${CONFIG_STR}_combined_vectors.root"

echo "Checking for combined data files..."
if [ ! -f "\$ACCEPTANCE_FILE" ]; then
    echo "❌ ERROR: Acceptance file not found: \$ACCEPTANCE_FILE"
    exit 1
fi

if [ ! -f "\$VECTORS_FILE" ]; then
    echo "❌ ERROR: Vectors file not found: \$VECTORS_FILE"  
    exit 1
fi

echo "✅ Both combined files found!"
echo "  Acceptance: \$ACCEPTANCE_FILE"
echo "  Vectors: \$VECTORS_FILE"

# Create configuration-specific acceptance plotting script
cat > plot_acceptance_${CONFIG_STR}.C << 'ACCEPTANCE_EOF'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <string>

void plot_acceptance_${CONFIG_STR}() {
    std::vector<std::string> filenames = {
        "${RESULTS_DIR}/${CONFIG_STR}_combined_acceptance.root"
    };

    double p_min = 0.0, p_max = 10.0;
    double bin_width = 0.5;
    int n_bins = int((p_max - p_min) / bin_width);
    std::vector<int> reco_in_bin(n_bins, 0);
    std::vector<int> total_in_bin(n_bins, 0);

    for (const auto& fname : filenames) {
        TFile *file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            std::cout << "Cannot open " << fname << std::endl;
            continue;
        }

        TTree *tree = (TTree*)file->Get("kaonEventInfo");
        if (!tree) {
            std::cout << "Tree kaonEventInfo not found!" << std::endl;
            file->Close();
            continue;
        }

        std::vector<double> *true_mom_vec = nullptr;
        std::vector<int> *reco_flag_vec = nullptr;
        int n_triple_pion_events = 0;
        tree->SetBranchAddress("true_mom_vec", &true_mom_vec);
        tree->SetBranchAddress("reco_flag_vec", &reco_flag_vec);
        tree->SetBranchAddress("n_triple_pion_events", &n_triple_pion_events);

        tree->GetEntry(0);
        std::cout << "Processing " << n_triple_pion_events << " kaon events..." << std::endl;

        for (int i = 0; i < n_triple_pion_events; ++i) {
            double p = (*true_mom_vec)[i];
            int bin = int((p - p_min) / bin_width);
            if (bin < 0 || bin >= n_bins) continue;
            total_in_bin[bin]++;
            if ((*reco_flag_vec)[i]) {
                reco_in_bin[bin]++;
            }
        }
        file->Close();
    }

    TH1D *h_acceptance = new TH1D("h_acceptance", "Kaon Reconstruction Acceptance (${CONFIG_STR});True Momentum [GeV/c];Acceptance", n_bins, p_min, p_max);

    int total_events = 0, total_reconstructed = 0;
    for (int bin = 0; bin < n_bins; ++bin) {
        int reconstructed_vertices = reco_in_bin[bin];
        int total_vertices = total_in_bin[bin];
        total_events += total_vertices;
        total_reconstructed += reconstructed_vertices;
        
        double acceptance = (total_vertices > 0) ? double(reconstructed_vertices) / double(total_vertices) : 0.0;
        h_acceptance->SetBinContent(bin + 1, acceptance);
    }

    TCanvas *c1 = new TCanvas("c1", "Kaon Acceptance ${CONFIG_STR}", 800, 600);
    h_acceptance->SetFillColor(kBlue-9);
    h_acceptance->SetLineColor(kBlue+2);
    h_acceptance->GetYaxis()->SetRangeUser(0, 1.1);
    h_acceptance->Draw("hist");
    c1->Update();

    gPad->SaveAs("${CONFIG_STR}_acceptance_plot.png");
    std::cout << "Acceptance plot saved as ${CONFIG_STR}_acceptance_plot.png" << std::endl;
    
    TFile *outfile = new TFile("${CONFIG_STR}_acceptance_histogram.root", "RECREATE");
    h_acceptance->Write();
    c1->Write();
    outfile->Close();
    
    std::cout << "=== ACCEPTANCE SUMMARY ===" << std::endl;
    std::cout << "Total kaon events: " << total_events << std::endl;
    std::cout << "Total reconstructed: " << total_reconstructed << std::endl;
    std::cout << "Overall acceptance: " << (total_events > 0 ? (double)total_reconstructed/total_events : 0) << std::endl;
}
ACCEPTANCE_EOF

# Create configuration-specific vectors plotting script  
cat > plot_vectors_${CONFIG_STR}.C << 'VECTORS_EOF'
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <cmath>

void plot_vectors_${CONFIG_STR}() {
    std::vector<std::string> filenames = {
        "${RESULTS_DIR}/${CONFIG_STR}_combined_vectors.root"
    };

    int nbins = 40;
    double min_p = 0.0, max_p = 10.0;
    TH2D *h2 = new TH2D("h2", "Resolution vs True Kaon Momentum (${CONFIG_STR});True Kaon Momentum [GeV/c];|Reco - True|/True", nbins, min_p, max_p, nbins, 0, 1);

    int total_kaons = 0;
    for (const auto& fname : filenames) {
        TFile *inFile = TFile::Open(fname.c_str());
        if (!inFile || inFile->IsZombie()) {
            std::cout << "Cannot open " << fname << std::endl;
            continue;
        }

        TTree *tree = (TTree*)inFile->Get("kaonVectors");
        if (!tree) { 
            std::cout << "Tree kaonVectors not found!" << std::endl;
            inFile->Close(); 
            continue; 
        }

        std::vector<double> *reco_p = nullptr;
        std::vector<double> *true_p = nullptr;
        tree->SetBranchAddress("reco_p", &reco_p);
        tree->SetBranchAddress("true_p", &true_p);
        
        if (tree->GetEntries() > 0) {
            tree->GetEntry(0);
            if (reco_p && true_p) {
                std::cout << "Found " << reco_p->size() << " kaon momentum pairs" << std::endl;
                total_kaons += reco_p->size();
                
                for (size_t i = 0; i < reco_p->size(); ++i) {
                    double true_momentum = (*true_p)[i];
                    double reco_momentum = (*reco_p)[i];
                    if (true_momentum > 0) {
                        double res = std::abs((reco_momentum - true_momentum) / true_momentum);
                        if (res <= 1.0) h2->Fill(true_momentum, res);
                    }
                }
            }
        }
        inFile->Close();
    }

    TH1D *h_mean = new TH1D("h_mean", "Mean Relative Kaon Momentum Resolution (${CONFIG_STR});True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)", nbins, min_p, max_p);
    
    for (int i = 1; i <= nbins; ++i) {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
        h_mean->SetBinContent(i, mean);
        delete proj;
    }

    TCanvas *c1 = new TCanvas("c1", "Kaon Resolution ${CONFIG_STR}", 800, 600);
    h_mean->SetFillColor(kRed-9);
    h_mean->SetBarWidth(0.9);
    h_mean->SetBarOffset(0.05);
    h_mean->Draw("BAR");
    c1->Update();

    gPad->SaveAs("${CONFIG_STR}_resolution_plot.png");
    std::cout << "Resolution plot saved as ${CONFIG_STR}_resolution_plot.png" << std::endl;
    
    TFile *outfile = new TFile("${CONFIG_STR}_resolution_histogram.root", "RECREATE");
    h_mean->Write();
    h2->Write();
    c1->Write();
    outfile->Close();
    
    std::cout << "=== RESOLUTION SUMMARY ===" << std::endl;
    std::cout << "Total kaons processed: " << total_kaons << std::endl;
}
VECTORS_EOF

# Generate histograms
echo "Generating acceptance histogram..."
root -l -b -q "plot_acceptance_${CONFIG_STR}.C()"

echo "Generating resolution histogram..."  
root -l -b -q "plot_vectors_${CONFIG_STR}.C()"

# Copy results to a summary directory
SUMMARY_DIR="/users/bp969/scratch/VIKING_FOLDER/PIPELINE_RESULTS/${CONFIG_STR}"
mkdir -p "\$SUMMARY_DIR"

echo "Copying results to summary directory: \$SUMMARY_DIR"
cp ${CONFIG_STR}_*.png "\$SUMMARY_DIR/" 2>/dev/null || echo "No PNG files to copy"
cp ${CONFIG_STR}_*.root "\$SUMMARY_DIR/" 2>/dev/null || echo "No ROOT files to copy"
cp "\$ACCEPTANCE_FILE" "\$SUMMARY_DIR/" 2>/dev/null || echo "Could not copy acceptance file"
cp "\$VECTORS_FILE" "\$SUMMARY_DIR/" 2>/dev/null || echo "Could not copy vectors file"

# Cleanup temp scripts
rm -f plot_acceptance_${CONFIG_STR}.C plot_vectors_${CONFIG_STR}.C

echo ""
echo "🎉 ==================================="
echo "🎉 PIPELINE COMPLETED SUCCESSFULLY!"
echo "🎉 ==================================="
echo ""
echo "Configuration: ${CONFIG_STR}"  
echo "Results directory: \$SUMMARY_DIR"
echo ""
echo "Generated files:"
ls -la "\$SUMMARY_DIR/"
echo ""
echo "✅ Acceptance histogram: ${CONFIG_STR}_acceptance_plot.png"
echo "✅ Resolution histogram: ${CONFIG_STR}_resolution_plot.png"
echo "✅ Combined data files copied to summary directory"
echo ""
echo "Your overnight run is complete!"
EOF

# Step 3: Update the submission script to include histogram generation
echo "=== Step 3: Updating submission script for overnight run ==="

# Modify the submission script to add histogram generation dependency
sed -i '/^COMBINE_JOB_ID=/a\
\
# Submit automatic histogram generation job (dependent on combination job)\
echo "Submitting automatic histogram generation job (dependent on combination job)..."\
HISTOGRAM_JOB_ID=$(sbatch --dependency=afterok:$COMBINE_JOB_ID /users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs/histograms_'${CONFIG_STR}'.job | awk '"'"'{print $4}'"'"')\
echo "Submitted histogram generation job: Job ID $HISTOGRAM_JOB_ID"' "/users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh"

# Update final status message
sed -i '/^echo "Results will be in: ${RESULTS_DIR}"/c\
echo "Results will be in: ${RESULTS_DIR}"\
echo "Final summary will be in: /users/bp969/scratch/VIKING_FOLDER/PIPELINE_RESULTS/'${CONFIG_STR}'"' "/users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh"

echo ""
echo "🌙 ==================================="
echo "🌙 OVERNIGHT PIPELINE READY!"
echo "🌙 ==================================="
echo ""
echo "Configuration: $CONFIG_STR"
echo ""
echo "To start your overnight run, execute:"
echo "  /users/bp969/scratch/JOBSCRIPTS_TESTS/submit_${CONFIG_STR}.sh"
echo ""
echo "This will:"
echo "  1. 🔨 Build 100 simulations (sequential, ~17 min total)"
echo "  2. 🚀 Run 100 simulations (parallel, ~2 days total)"  
echo "  3. 📊 Parse data from all runs (~2 hours)"
echo "  4. 🔗 Combine all data files (~30 min)"
echo "  5. 📈 Generate histograms automatically (~5 min)"
echo ""  
echo "Expected completion time: ~2 days (10 million events)"
echo "Final results location: /users/bp969/scratch/VIKING_FOLDER/PIPELINE_RESULTS/${CONFIG_STR}/"
echo ""
echo "💤 Perfect for 2-day production runs! You'll have:"
echo "  ✅ Acceptance histogram (PNG + ROOT)"
echo "  ✅ Resolution histogram (PNG + ROOT)" 
echo "  ✅ Combined data files"
echo "  ✅ Email notifications on completion/failure"
echo ""
echo "Start now with: ./submit_${CONFIG_STR}.sh"
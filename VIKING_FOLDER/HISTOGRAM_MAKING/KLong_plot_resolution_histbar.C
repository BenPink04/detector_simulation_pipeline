#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <cmath>

void KLong_plot_resolution_histbar() {
    std::vector<std::string> filenames = {
    "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700_combined_vectors.root"
    };

    // Extract configuration string from filename
    std::string config_string = "";
    if (!filenames.empty()) {
        std::string fname = filenames[0];
        
        // Find the last occurrence of the configuration pattern (after the last '/')
        size_t last_slash = fname.find_last_of('/');
        std::string basename = (last_slash != std::string::npos) ? fname.substr(last_slash + 1) : fname;
        
        size_t start = basename.find("T1-");
        size_t end = basename.find("_combined_vectors.root");
        if (start != std::string::npos && end != std::string::npos) {
            config_string = basename.substr(start, end - start);
        }
    }

    int nbins = 40;
    double min_p = 0.0;
    double max_p = 10.0;
    
    std::string h2_title = "Resolution vs True Kaon Momentum";
    if (!config_string.empty()) {
        h2_title += " (Config: " + config_string + ")";
    }
    h2_title += ";True Kaon Momentum [GeV/c];|Reco - True|/True";
    
    TH2D *h2 = new TH2D("h2", h2_title.c_str(), nbins, min_p, max_p, nbins, 0, 1);

    // Fill histogram with all data, skipping anomalous results
    double anomaly_threshold = 1.0; // Ignore resolution > 1 (100%)
    for (const auto& fname : filenames) {
        TFile *inFile = TFile::Open(fname.c_str());
        if (!inFile || inFile->IsZombie()) continue;

        TTree *tree = (TTree*)inFile->Get("kaonVectors");
        if (!tree) { inFile->Close(); continue; }

        std::vector<double> *reco_p = nullptr;
        std::vector<double> *true_p = nullptr;
        tree->SetBranchAddress("reco_p", &reco_p);
        tree->SetBranchAddress("true_p", &true_p);
        
        // Loop through ALL entries in the combined file (each entry contains vectors from one original file)
        Long64_t nEntries = tree->GetEntries();
        std::cout << "Processing " << nEntries << " entries from " << fname << std::endl;
        
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            
            // Process all events in this entry's vectors
            for (size_t i = 0; i < reco_p->size(); ++i) {
                double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
                if (res > anomaly_threshold) continue; // Skip anomalous results
                h2->Fill((*true_p)[i], res);
            }
        }
        inFile->Close();
    }

    // Make a 1D histogram of mean resolution per bin
    std::string h_mean_title = "Mean Relative Kaon Momentum Resolution";
    if (!config_string.empty()) {
        h_mean_title += " (Config: " + config_string + ")";
    }
    h_mean_title += ";True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)";
    
    TH1D *h_mean = new TH1D("h_mean", h_mean_title.c_str(), nbins, min_p, max_p);
    for (int i = 1; i <= nbins; ++i) {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
        h_mean->SetBinContent(i, mean);
        delete proj;
    }

    TCanvas *c1 = new TCanvas("c1", "Mean Kaon Relative Momentum Resolution", 800, 600);
    h_mean->SetFillColor(kBlue-9);
    h_mean->SetBarWidth(0.9);
    h_mean->SetBarOffset(0.05);
    h_mean->Draw("BAR");
    c1->Update();

    // Create filename with configuration string
    std::string output_filename = "KLong_histbar_plot";
    if (!config_string.empty()) {
        output_filename += "_" + config_string;
    }
    output_filename += ".png";
    
    // Save the plot as PNG files
    gPad->SaveAs(output_filename.c_str());
    
    std::cout << "Plot saved as " << output_filename << std::endl;
}
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <string>

void KLong_plot_acceptance() {
    std::vector<std::string> filenames = {
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700_combined_acceptance.root"
        // Add more filenames as needed
    };

    // Extract configuration string from filename
    std::string config_string = "";
    if (!filenames.empty()) {
        std::string fname = filenames[0];
        
        // Find the last occurrence of the configuration pattern (after the last '/')
        size_t last_slash = fname.find_last_of('/');
        std::string basename = (last_slash != std::string::npos) ? fname.substr(last_slash + 1) : fname;
        
        size_t start = basename.find("T1-");
        size_t end = basename.find("_combined_acceptance.root");
        if (start != std::string::npos && end != std::string::npos) {
            config_string = basename.substr(start, end - start);
        }
    }

    double p_min = 0.0, p_max = 10.0;
    double bin_width = 0.5;
    int n_bins = int((p_max - p_min) / bin_width);
    std::vector<int> reco_in_bin(n_bins, 0);
    std::vector<int> total_in_bin(n_bins, 0);

    for (const auto& fname : filenames) {
        TFile *file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            std::cout << "Cannot open " << fname << "\n";
            continue;
        }

        TTree *tree = (TTree*)file->Get("kaonEventInfo");
        if (!tree) {
            std::cout << "Tree kaonEventInfo not found in " << fname << "!\n";
            file->Close();
            continue;
        }

        std::vector<double> *true_mom_vec = nullptr;
        std::vector<int> *reco_flag_vec = nullptr;
        int n_triple_pion_events = 0;
        tree->SetBranchAddress("true_mom_vec", &true_mom_vec);
        tree->SetBranchAddress("reco_flag_vec", &reco_flag_vec);
        tree->SetBranchAddress("n_triple_pion_events", &n_triple_pion_events);

        // Loop through ALL entries in the combined file (each entry contains data from one original file)
        Long64_t nEntries = tree->GetEntries();
        std::cout << "Processing " << nEntries << " entries from " << fname << std::endl;
        
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            
            // Process all events in this entry
            for (int i = 0; i < n_triple_pion_events; ++i) {
                double p = (*true_mom_vec)[i];
                int bin = int((p - p_min) / bin_width);
                if (bin < 0 || bin >= n_bins) continue;
                total_in_bin[bin]++;
                if ((*reco_flag_vec)[i]) {
                    reco_in_bin[bin]++;
                }
            }
        }
        file->Close();
    }

    std::string title = "Reconstruction Acceptance";
    if (!config_string.empty()) {
        title += " (Config: " + config_string + ")";
    }
    title += ";True Momentum [GeV/c];Acceptance (Reconstructed/Total)";
    
    TH1D *h_proportion = new TH1D("h_proportion", title.c_str(), n_bins, p_min, p_max);

    for (int bin = 0; bin < n_bins; ++bin) {
        int reconstructed_vertices = reco_in_bin[bin];
        int total_vertices = total_in_bin[bin];
        
        double acceptance;
        if (total_vertices > 0) {
            acceptance = double(reconstructed_vertices) / double(total_vertices);
        } else {
            acceptance = 0.0; // No vertices in this bin
        }
        
        h_proportion->SetBinContent(bin + 1, acceptance);
        
        // Debug output to see what's happening in each bin
        if (total_vertices > 0) {
            std::cout << "Bin " << bin << " (p=" << p_min + bin*bin_width << "-" << p_min + (bin+1)*bin_width 
                      << " GeV/c): Total=" << total_vertices 
                      << ", Reconstructed=" << reconstructed_vertices 
                      << ", Acceptance=" << acceptance << std::endl;
        }
    }

    h_proportion->Draw("hist");

    // Create filename with configuration string
    std::string output_filename = "KLong_acceptance_plot";
    if (!config_string.empty()) {
        output_filename += "_" + config_string;
    }
    output_filename += ".png";
    
    // Save the plot as PNG files
    gPad->SaveAs(output_filename.c_str());

    std::cout << "Plot saved as " << output_filename << std::endl;
}
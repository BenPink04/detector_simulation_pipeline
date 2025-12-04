#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TText.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

void KLong_plot_compare_acceptance() {
    // Define multiple files to compare - using actual existing paths
    std::vector<std::string> filenames = {
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-260_F2-270_E1-600_combined_acceptance.root",
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700_combined_acceptance.root",
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600_combined_acceptance.root",
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-400_F2-410_E1-600/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-400_F2-410_E1-600_combined_acceptance.root"
    };

    // Extract configuration labels from filenames
    std::vector<std::string> config_labels;
    for (const auto& fname : filenames) {
        std::string label = "";
        
        // Find the last occurrence of the configuration pattern (after the last '/')
        size_t last_slash = fname.find_last_of('/');
        std::string basename = (last_slash != std::string::npos) ? fname.substr(last_slash + 1) : fname;
        
        size_t start = basename.find("T1-");
        size_t end = basename.find("_combined_acceptance.root");
        if (start != std::string::npos && end != std::string::npos) {
            label = basename.substr(start, end - start);
            // Keep full configuration label (all detector positions)
        } else {
            label = "Config_" + std::to_string(config_labels.size() + 1);
        }
        
        config_labels.push_back(label);
    }

    // Histogram parameters (matching your acceptance plot)
    double p_min = 0.0, p_max = 10.0;
    double bin_width = 0.5;
    int n_bins = int((p_max - p_min) / bin_width);

    // Colors for different configurations
    std::vector<int> colors = {kBlue, kRed, kGreen+2, kMagenta, kCyan+2, kOrange+2, kViolet, kGray+2};
    
    // Create canvas with divided layout
    TCanvas *c1 = new TCanvas("c1", "Kaon Reconstruction Acceptance Comparison", 1400, 700);
    c1->Divide(2, 1);  // Divide canvas: left for histogram, right for legend
    c1->cd(1);         // Switch to left pad for histogram
    gPad->SetPad(0.0, 0.0, 0.75, 1.0);  // Left pad takes 75% of width
    c1->cd(2);         // Switch to right pad for legend
    gPad->SetPad(0.75, 0.0, 1.0, 1.0);  // Right pad takes 25% of width
    c1->cd(1);         // Go back to left pad for drawing histogram
    
    // Vector to store histograms
    std::vector<TH1D*> h_acceptance_vec;
    
    // Process each file
    for (size_t file_idx = 0; file_idx < filenames.size(); ++file_idx) {
        const std::string& fname = filenames[file_idx];
        const std::string& config_label = config_labels[file_idx];
        
        std::cout << "Processing file " << (file_idx + 1) << "/" << filenames.size() 
                  << ": " << config_label << std::endl;
        
        // Initialize counters for this configuration
        std::vector<int> reco_in_bin(n_bins, 0);
        std::vector<int> total_in_bin(n_bins, 0);

        // Process the acceptance file
        TFile *file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            std::cout << "Warning: Cannot open file " << fname << std::endl;
            continue;
        }

        std::cout << "  Successfully opened file: " << fname << std::endl;
        file->ls(); // List contents for debugging

        TTree *tree = (TTree*)file->Get("kaonEventInfo");
        if (!tree) {
            std::cout << "Warning: Tree kaonEventInfo not found in " << fname << std::endl;
            file->Close();
            continue;
        }

        std::cout << "  Tree found with " << tree->GetEntries() << " entries" << std::endl;

        std::vector<double> *true_mom_vec = nullptr;
        std::vector<int> *reco_flag_vec = nullptr;
        int n_triple_pion_events = 0;
        tree->SetBranchAddress("true_mom_vec", &true_mom_vec);
        tree->SetBranchAddress("reco_flag_vec", &reco_flag_vec);
        tree->SetBranchAddress("n_triple_pion_events", &n_triple_pion_events);

        // Loop through ALL entries in the combined file
        Long64_t nEntries = tree->GetEntries();
        std::cout << "  Processing " << nEntries << " entries..." << std::endl;
        
        int total_events = 0;
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
                total_events++;
            }
        }
        
        std::cout << "  Total kaon events processed: " << total_events << std::endl;
        file->Close();

        // Create acceptance histogram for this configuration
        std::string h_name = "h_acceptance_" + std::to_string(file_idx);
        std::string h_title = "Reconstruction Acceptance Comparison";
        h_title += ";True Momentum [GeV/c];Acceptance (Reconstructed/Total)";
        
        TH1D *h_acceptance = new TH1D(h_name.c_str(), h_title.c_str(), n_bins, p_min, p_max);

        // Fill acceptance histogram
        for (int bin = 0; bin < n_bins; ++bin) {
            int reconstructed_vertices = reco_in_bin[bin];
            int total_vertices = total_in_bin[bin];
            
            double acceptance;
            if (total_vertices > 0) {
                acceptance = double(reconstructed_vertices) / double(total_vertices);
            } else {
                acceptance = 0.0; // No vertices in this bin
            }
            
            h_acceptance->SetBinContent(bin + 1, acceptance);
        }
        
        // Style the histogram for outlined bars (no fill)
        h_acceptance->SetFillStyle(0); // No fill - hollow bars
        h_acceptance->SetLineColor(colors[file_idx % colors.size()]);
        h_acceptance->SetLineWidth(3); // Thick outline
        h_acceptance->SetBarWidth(0.8);
        h_acceptance->SetBarOffset(0.1);
        
        h_acceptance_vec.push_back(h_acceptance);
    }
    
    if (h_acceptance_vec.empty()) {
        std::cout << "Error: No valid histograms created!" << std::endl;
        return;
    }
    
    // Draw histograms
    gStyle->SetOptStat(0); // Turn off statistics box
    
    // Set up the first histogram
    h_acceptance_vec[0]->SetTitle("Kaon Reconstruction Acceptance Comparison");
    h_acceptance_vec[0]->GetYaxis()->SetRangeUser(0, 0.010); // Acceptance range 0-1.5% to better show the data
    h_acceptance_vec[0]->Draw("BAR");
    
    // Draw additional histograms
    for (size_t i = 1; i < h_acceptance_vec.size(); ++i) {
        h_acceptance_vec[i]->Draw("BAR SAME");
    }
    
    // Switch to right pad for legend
    c1->cd(2);
    
    // Clear the pad and set up for manual drawing
    gPad->Clear();
    gPad->Range(0, 0, 1, 1);
    
    // Draw simple legend with smaller text positioned to the left
    double y_start = 0.9;
    double y_step = 0.8 / h_acceptance_vec.size();  // Distribute entries evenly
    
    for (size_t i = 0; i < h_acceptance_vec.size(); ++i) {
        double entry_y = y_start - i * y_step;
        
        // Draw the color line indicator
        TLine *colorLine = new TLine(0.01, entry_y, 0.08, entry_y);
        colorLine->SetLineColor(h_acceptance_vec[i]->GetLineColor());
        colorLine->SetLineWidth(3);
        colorLine->Draw();
        
        // Draw text with smaller size and positioned further left
        TText *text = new TText(0.10, entry_y, config_labels[i].c_str());
        text->SetTextSize(0.025);  // Smaller text size
        text->SetTextColor(kBlack);
        text->Draw();
    }
    
    // Switch back to left pad
    c1->cd(1);
    
    // Add grid for better readability
    c1->SetGrid(1, 1);
    c1->Update();
    
    // Save the plot - save entire canvas to include both pads
    std::string output_filename = "KLong_acceptance_comparison.png";
    c1->SaveAs(output_filename.c_str());
    
    std::cout << "\nComparison plot saved as: " << output_filename << std::endl;
    
    // Print summary statistics
    std::cout << "\nSummary Statistics:" << std::endl;
    std::cout << "Configuration\t\tMean Acceptance (0.5-2.0 GeV)\tMean Acceptance (2.0-5.0 GeV)" << std::endl;
    std::cout << "-------------\t\t-----------------------------\t-----------------------------" << std::endl;
    
    for (size_t i = 0; i < h_acceptance_vec.size(); ++i) {
        TH1D* h = h_acceptance_vec[i];
        
        // Calculate mean acceptance in different momentum ranges
        double sum_low = 0, count_low = 0;
        double sum_high = 0, count_high = 0;
        
        for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
            double momentum = h->GetBinCenter(bin);
            double acceptance = h->GetBinContent(bin);
            
            if (acceptance > 0) {
                if (momentum >= 0.5 && momentum <= 2.0) {
                    sum_low += acceptance;
                    count_low++;
                } else if (momentum >= 2.0 && momentum <= 5.0) {
                    sum_high += acceptance;
                    count_high++;
                }
            }
        }
        
        double mean_low = (count_low > 0) ? sum_low / count_low : 0;
        double mean_high = (count_high > 0) ? sum_high / count_high : 0;
        
        std::cout << config_labels[i] << "\t\t" 
                  << std::fixed << std::setprecision(4) << mean_low << "\t\t\t\t"
                  << mean_high << std::endl;
    }
}
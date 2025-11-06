#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <string>

// Function to create an acceptance histogram from an acceptance file
TH1D* createAcceptanceHist(const std::string& filename, const std::string& histname, int color) {
    double p_min = 0.0, p_max = 10.0;
    double bin_width = 0.5;
    int n_bins = int((p_max - p_min) / bin_width);
    std::vector<int> reco_in_bin(n_bins, 0);
    std::vector<int> total_in_bin(n_bins, 0);
    
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return nullptr;
    }
    
    TTree *tree = (TTree*)file->Get("kaonEventInfo");
    if (!tree) {
        std::cerr << "Error: Tree kaonEventInfo not found in " << filename << std::endl;
        file->Close();
        return nullptr;
    }
    
    std::vector<double> *true_mom_vec = nullptr;
    std::vector<int> *reco_flag_vec = nullptr;
    int n_triple_pion_events = 0;
    tree->SetBranchAddress("true_mom_vec", &true_mom_vec);
    tree->SetBranchAddress("reco_flag_vec", &reco_flag_vec);
    tree->SetBranchAddress("n_triple_pion_events", &n_triple_pion_events);
    
    Long64_t nentries = tree->GetEntries();
    std::cout << "Processing " << filename << " with " << nentries << " entries..." << std::endl;
    
    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
        tree->GetEntry(ientry);
        
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
    
    TH1D *h_proportion = new TH1D(histname.c_str(), 
                                   "Reconstruction Acceptance;True Momentum [GeV/c];Acceptance (Reconstructed/Total)", 
                                   n_bins, p_min, p_max);
    
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
    }
    
    // Set histogram style
    h_proportion->SetLineColor(color);
    h_proportion->SetLineWidth(2);
    h_proportion->SetStats(0);
    
    return h_proportion;
}

void KLong_plot_acceptance_comparison_test() {
    // Define the combined acceptance files to compare
    // Format: {filename, legend label, color}
    std::vector<std::tuple<std::string, std::string, int>> configurations = {
        {"/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/SCEN5_DATA_READ/out_seed_Scen5_1_acceptance.root", 
         "Config 1:", kBlue},
        {"/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/SCEN5_DATA_READ/out_seed_Scen5_2_acceptance.root", 
         "Config 2:", kRed}, 
        {"/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/SCEN5_DATA_READ/out_seed_Scen5_3_acceptance.root", 
         "Config 3:", kGreen+2}

        // Add more configurations here as needed
        // Example:
        // {"T1-300_T2-310_T3-380_T4-390_P1-215_P2-230_F1-290_F2-300_E1-400_combined_acceptance.root", 
        //  "Config 2: T1=300, T2=310, T3=380, T4=390", kRed},
    };
    
    TCanvas *c1 = new TCanvas("c1", "Kaon Acceptance Comparison", 1000, 700);
    TLegend *legend = new TLegend(0.15, 0.15, 0.55, 0.38);
    legend->SetBorderSize(1);
    legend->SetFillStyle(0);
    
    std::vector<TH1D*> histograms;
    bool first = true;
    double max_y = 0;
    
    // Create histograms for each configuration
    for (size_t i = 0; i < configurations.size(); ++i) {
        const auto& config = configurations[i];
        std::string filename = std::get<0>(config);
        std::string label = std::get<1>(config);
        int color = std::get<2>(config);
        
        TH1D* hist = createAcceptanceHist(filename, Form("h_acceptance_%zu", i), color);
        
        if (!hist) {
            std::cerr << "Warning: Skipping " << filename << " (could not create histogram)" << std::endl;
            continue;
        }
        
        histograms.push_back(hist);
        
        // Track maximum y value for proper scaling
        double this_max = hist->GetMaximum();
        if (this_max > max_y) max_y = this_max;
        
        // Draw the histogram
        if (first) {
            hist->Draw("HIST");
            first = false;
        } else {
            hist->Draw("HIST SAME");
        }
        
        legend->AddEntry(hist, label.c_str(), "l");
    }
    
    if (histograms.empty()) {
        std::cerr << "Error: No valid histograms created!" << std::endl;
        return;
    }
    
    // Set y-axis range with some padding (acceptance is 0-1)
    histograms[0]->SetMaximum(max_y * 1.15 > 1.0 ? 1.0 : max_y * 1.15);
    histograms[0]->SetMinimum(0);
    
    legend->Draw();
    c1->Update();
    
    // Save the canvas
    c1->SaveAs("kaon_acceptance_comparison.png");
    c1->SaveAs("kaon_acceptance_comparison.pdf");
    
    std::cout << "\nComparison plot saved as kaon_acceptance_comparison.png/pdf" << std::endl;
}

// Alternative version that takes filenames as arguments
void KLong_plot_acceptance_comparison(std::vector<std::string> filenames, 
                                       std::vector<std::string> labels,
                                       std::vector<int> colors = {}) {
    if (filenames.size() != labels.size()) {
        std::cerr << "Error: Number of filenames must match number of labels!" << std::endl;
        return;
    }
    
    // Default colors if not provided
    if (colors.empty()) {
        colors = {kBlue, kRed, kGreen+2, kMagenta+2, kOrange+1, kCyan+2, kViolet, kPink+1};
    }
    
    TCanvas *c1 = new TCanvas("c1", "Kaon Acceptance Comparison", 1000, 700);
    TLegend *legend = new TLegend(0.15, 0.15, 0.55, 0.38);
    legend->SetBorderSize(1);
    legend->SetFillStyle(0);
    
    std::vector<TH1D*> histograms;
    bool first = true;
    double max_y = 0;
    
    // Create histograms for each file
    for (size_t i = 0; i < filenames.size(); ++i) {
        int color = colors[i % colors.size()];
        
        TH1D* hist = createAcceptanceHist(filenames[i], Form("h_acceptance_%zu", i), color);
        
        if (!hist) {
            std::cerr << "Warning: Skipping " << filenames[i] << " (could not create histogram)" << std::endl;
            continue;
        }
        
        histograms.push_back(hist);
        
        // Track maximum y value for proper scaling
        double this_max = hist->GetMaximum();
        if (this_max > max_y) max_y = this_max;
        
        // Draw the histogram
        if (first) {
            hist->Draw("HIST");
            first = false;
        } else {
            hist->Draw("HIST SAME");
        }
        
        legend->AddEntry(hist, labels[i].c_str(), "l");
    }
    
    if (histograms.empty()) {
        std::cerr << "Error: No valid histograms created!" << std::endl;
        return;
    }
    
    // Set y-axis range with some padding (acceptance is 0-1)
    histograms[0]->SetMaximum(max_y * 1.15 > 1.0 ? 1.0 : max_y * 1.15);
    histograms[0]->SetMinimum(0);
    
    legend->Draw();
    c1->Update();
    
    // Save the canvas
    c1->SaveAs("kaon_acceptance_comparison.png");
    c1->SaveAs("kaon_acceptance_comparison.pdf");
    
    std::cout << "\nComparison plot saved as kaon_acceptance_comparison.png/pdf" << std::endl;
}

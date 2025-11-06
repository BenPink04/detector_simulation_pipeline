#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLegend.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>

// Function to create a mean resolution histogram from a vectors file
TH1D* createMeanResolutionHist(const std::string& filename, const std::string& histname, int color) {
    int nbins = 40;
    double min_p = 0.0;
    double max_p = 10.0;
    
    TH2D *h2 = new TH2D(Form("h2_%s", histname.c_str()), 
                        "Resolution vs True Kaon Momentum;True Kaon Momentum [GeV/c];|Reco - True|/True", 
                        nbins, min_p, max_p, nbins, 0, 1);
    
    // Fill histogram with data, skipping anomalous results
    double anomaly_threshold = 1.0; // Ignore resolution > 1 (100%)
    
    TFile *inFile = TFile::Open(filename.c_str());
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        delete h2;
        return nullptr;
    }
    
    TTree *tree = (TTree*)inFile->Get("kaonVectors");
    if (!tree) {
        std::cerr << "Error: Could not find kaonVectors tree in " << filename << std::endl;
        inFile->Close();
        delete h2;
        return nullptr;
    }
    
    std::vector<double> *reco_p = nullptr;
    std::vector<double> *true_p = nullptr;
    tree->SetBranchAddress("reco_p", &reco_p);
    tree->SetBranchAddress("true_p", &true_p);
    
    Long64_t nentries = tree->GetEntries();
    std::cout << "Processing " << filename << " with " << nentries << " entries..." << std::endl;
    
    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
        tree->GetEntry(ientry);
        for (size_t i = 0; i < reco_p->size(); ++i) {
            double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
            if (res > anomaly_threshold) continue; // Skip anomalous results
            h2->Fill((*true_p)[i], res);
        }
    }
    
    inFile->Close();
    
    // Make a 1D histogram of mean resolution per bin
    TH1D *h_mean = new TH1D(histname.c_str(), 
                            "Mean Relative Kaon Momentum Resolution;True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)", 
                            nbins, min_p, max_p);
    
    for (int i = 1; i <= nbins; ++i) {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
        h_mean->SetBinContent(i, mean);
        delete proj;
    }
    
    delete h2;
    
    // Set histogram style
    h_mean->SetLineColor(color);
    h_mean->SetLineWidth(2);
    h_mean->SetStats(0);
    
    return h_mean;
}

void KLong_plot_resolution_comparison() {
    // Define the combined vector files to compare
    // Format: {filename, legend label, color}
    std::vector<std::tuple<std::string, std::string, int>> configurations = {
        {"T1-250_T2-260_T3-330_T4-340_P1-215_P2-230_F1-290_F2-300_E1-400_combined_vectors.root", 
         "Config 1: T1=250, T2=260, T3=330, T4=340", kBlue},
        // Add more configurations here as needed
        // Example:
        // {"T1-300_T2-310_T3-380_T4-390_P1-215_P2-230_F1-290_F2-300_E1-400_combined_vectors.root", 
        //  "Config 2: T1=300, T2=310, T3=380, T4=390", kRed},
    };
    
    TCanvas *c1 = new TCanvas("c1", "Kaon Momentum Resolution Comparison", 1000, 700);
    TLegend *legend = new TLegend(0.15, 0.65, 0.55, 0.88);
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
        
        TH1D* hist = createMeanResolutionHist(filename, Form("h_mean_%zu", i), color);
        
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
    
    // Set y-axis range with some padding
    histograms[0]->SetMaximum(max_y * 1.15);
    histograms[0]->SetMinimum(0);
    
    legend->Draw();
    c1->Update();
    
    // Save the canvas
    c1->SaveAs("kaon_resolution_comparison.png");
    c1->SaveAs("kaon_resolution_comparison.pdf");
    
    std::cout << "\nComparison plot saved as kaon_resolution_comparison.png/pdf" << std::endl;
}

// Alternative version that takes filenames as arguments
void KLong_plot_resolution_comparison(std::vector<std::string> filenames, 
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
    
    TCanvas *c1 = new TCanvas("c1", "Kaon Momentum Resolution Comparison", 1000, 700);
    TLegend *legend = new TLegend(0.15, 0.65, 0.55, 0.88);
    legend->SetBorderSize(1);
    legend->SetFillStyle(0);
    
    std::vector<TH1D*> histograms;
    bool first = true;
    double max_y = 0;
    
    // Create histograms for each file
    for (size_t i = 0; i < filenames.size(); ++i) {
        int color = colors[i % colors.size()];
        
        TH1D* hist = createMeanResolutionHist(filenames[i], Form("h_mean_%zu", i), color);
        
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
    
    // Set y-axis range with some padding
    histograms[0]->SetMaximum(max_y * 1.15);
    histograms[0]->SetMinimum(0);
    
    legend->Draw();
    c1->Update();
    
    // Save the canvas
    c1->SaveAs("kaon_resolution_comparison.png");
    c1->SaveAs("kaon_resolution_comparison.pdf");
    
    std::cout << "\nComparison plot saved as kaon_resolution_comparison.png/pdf" << std::endl;
}

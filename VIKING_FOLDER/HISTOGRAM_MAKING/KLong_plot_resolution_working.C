#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TPaveText.h"
#include <vector>
#include <iostream>
#include <cmath>

void KLong_plot_resolution_working() {
    // Use the specified combined vectors file
    std::string filename = "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_combined_vectors.root";
    
    std::cout << "Processing combined vectors file: " << filename << std::endl;

    int nbins = 20;
    double min_p = 0.0;
    double max_p = 10.0;
    TH2D *h2 = new TH2D("h2", "Kaon Momentum Resolution;True Kaon Momentum [GeV/c];Relative Resolution |Reco - True|/True", 
                        nbins, min_p, max_p, 100, 0, 1);

    // Fill histogram with all data, skipping anomalous results
    double anomaly_threshold = 1.0; // Ignore resolution > 100%
    int total_entries = 0;
    
    TFile *inFile = TFile::Open(filename.c_str());
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        return;
    }

    TTree *tree = (TTree*)inFile->Get("kaonVectors");
    if (!tree) { 
        std::cerr << "Error: No kaonVectors tree in: " << filename << std::endl;
        inFile->Close(); 
        return; 
    }

    std::vector<double> *reco_p = nullptr;
    std::vector<double> *true_p = nullptr;
    tree->SetBranchAddress("reco_p", &reco_p);
    tree->SetBranchAddress("true_p", &true_p);
    
    Long64_t nentries = tree->GetEntries();
    std::cout << "Processing " << nentries << " entries from combined vectors file..." << std::endl;
    
    for (Long64_t entry = 0; entry < nentries; entry++) {
        tree->GetEntry(entry);
        
        if (!reco_p || !true_p) continue;
        
        for (size_t i = 0; i < reco_p->size() && i < true_p->size(); ++i) {
            if ((*true_p)[i] <= 0) continue; // Skip invalid true momentum
            
            double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
            if (res > anomaly_threshold) continue; // Skip anomalous results
            
            h2->Fill((*true_p)[i], res);
            total_entries++;
        }
        
        if (entry % 1000 == 0) {
            std::cout << "Processed " << entry << " entries..." << std::endl;
        }
    }
    
    inFile->Close();
    std::cout << "Processed combined file with " << total_entries << " valid kaon momentum pairs" << std::endl;

    // Make a 1D histogram of mean resolution per bin
    TH1D *h_mean = new TH1D("h_mean", "Mean Relative Kaon Momentum Resolution;True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)", 
                            nbins, min_p, max_p);
    
    for (int i = 1; i <= nbins; ++i) {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
        h_mean->SetBinContent(i, mean);
        delete proj;
    }

    // Create canvas and plot
    TCanvas *c1 = new TCanvas("c1", "Kaon Momentum Resolution - Combined Data", 800, 600);
    h_mean->SetFillColor(kBlue-9);
    h_mean->SetBarWidth(0.9);
    h_mean->SetBarOffset(0.05);
    h_mean->Draw("BAR");
    
    // Add statistics text
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.85, 0.85, "NDC");
    pt->AddText("Combined vectors file");
    pt->AddText(Form("Total momentum pairs: %d", total_entries));
    pt->AddText("T1=240, T2=250, T3=570, T4=580");
    pt->AddText("P1=215, P2=230, F1=260, F2=270");
    pt->Draw();
    
    c1->Update();

    // Save the plot
    c1->SaveAs("kaon_resolution_working.png");
    std::cout << "Plot saved as kaon_resolution_working.png" << std::endl;
}
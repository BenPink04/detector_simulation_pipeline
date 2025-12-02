#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <cmath>

void plot_vectors_T1_240_T2_250_T3_570_T4_580_P1_215_P2_230_F1_260_F2_270_E1_600() {
    std::vector<std::string> filenames = {
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_combined_vectors.root"
    };

    int nbins = 40;
    double min_p = 0.0, max_p = 10.0;
    TH2D *h2 = new TH2D("h2", "Resolution vs True Kaon Momentum (T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600);True Kaon Momentum [GeV/c];|Reco - True|/True", nbins, min_p, max_p, nbins, 0, 1);

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

    TH1D *h_mean = new TH1D("h_mean", "Mean Relative Kaon Momentum Resolution (T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600);True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)", nbins, min_p, max_p);
    
    for (int i = 1; i <= nbins; ++i) {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
        h_mean->SetBinContent(i, mean);
        delete proj;
    }

    TCanvas *c1 = new TCanvas("c1", "Kaon Resolution T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600", 800, 600);
    h_mean->SetFillColor(kRed-9);
    h_mean->SetBarWidth(0.9);
    h_mean->SetBarOffset(0.05);
    h_mean->Draw("BAR");
    c1->Update();

    gPad->SaveAs("T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_resolution_plot.png");
    std::cout << "Resolution plot saved as T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_resolution_plot.png" << std::endl;
    
    TFile *outfile = new TFile("T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_resolution_histogram.root", "RECREATE");
    h_mean->Write();
    h2->Write();
    c1->Write();
    outfile->Close();
    
    std::cout << "=== RESOLUTION SUMMARY ===" << std::endl;
    std::cout << "Total kaons processed: " << total_kaons << std::endl;
}

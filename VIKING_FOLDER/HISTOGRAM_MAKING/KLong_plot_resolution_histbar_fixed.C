#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <cmath>

void KLong_plot_resolution_histbar_fixed() {
    std::vector<std::string> filenames = {
        "/users/bp969/scratch/VIKING_FOLDER/PIPELINE_RESULTS/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600_combined_vectors.root"
        };

    int nbins = 40;
    double min_p = 0.0;
    double max_p = 10.0;
    TH2D *h2 = new TH2D("h2", "Resolution vs True Kaon Momentum;True Kaon Momentum [GeV/c];|Reco - True|/True", nbins, min_p, max_p, nbins, 0, 1);

    // Fill histogram with all data, skipping anomalous results
    double anomaly_threshold = 1.0; // Ignore resolution > 1 (100%)
    int total_kaons = 0;
    
    for (const auto& fname : filenames) {
        TFile *inFile = TFile::Open(fname.c_str());
        if (!inFile || inFile->IsZombie()) {
            std::cout << "Cannot open " << fname << std::endl;
            continue;
        }

        TTree *tree = (TTree*)inFile->Get("kaonVectors");
        if (!tree) { 
            std::cout << "Tree kaonVectors not found in " << fname << std::endl;
            inFile->Close(); 
            continue; 
        }

        std::cout << "Processing " << fname << " with " << tree->GetEntries() << " entries" << std::endl;

        // Updated branch addresses for fixed structure
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
                        if (res <= anomaly_threshold) { // Skip anomalous results
                            h2->Fill(true_momentum, res);
                        }
                    }
                }
            } else {
                std::cout << "Failed to load momentum vectors!" << std::endl;
            }
        } else {
            std::cout << "No entries in tree!" << std::endl;
        }
        
        inFile->Close();
    }

    std::cout << "Total kaons processed: " << total_kaons << std::endl;

    // Make a 1D histogram of mean resolution per bin
    TH1D *h_mean = new TH1D("h_mean", "Mean Relative Kaon Momentum Resolution;True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)", nbins, min_p, max_p);
    
    for (int i = 1; i <= nbins; ++i) {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
        h_mean->SetBinContent(i, mean);
        
        if (proj->GetEntries() > 0) {
            std::cout << "Bin " << i << " (p=" << min_p + (i-1)*(max_p-min_p)/nbins 
                      << "-" << min_p + i*(max_p-min_p)/nbins << " GeV/c): "
                      << proj->GetEntries() << " entries, mean resolution=" << mean << std::endl;
        }
        
        delete proj;
    }

    TCanvas *c1 = new TCanvas("c1", "Mean Kaon Relative Momentum Resolution", 800, 600);
    h_mean->SetFillColor(kBlue-9);
    h_mean->SetBarWidth(0.9);
    h_mean->SetBarOffset(0.05);
    h_mean->Draw("BAR");
    c1->Update();

    // Save the plot as PNG files
    gPad->SaveAs("KLong_resolution_histbar_fixed.png");
    std::cout << "Plot saved as KLong_resolution_histbar_fixed.png" << std::endl;
    
    // Also save histogram as ROOT file
    TFile *outfile = new TFile("KLong_resolution_histogram_fixed.root", "RECREATE");
    h_mean->Write();
    h2->Write();
    c1->Write();
    outfile->Close();
    std::cout << "Histograms saved as KLong_resolution_histogram_fixed.root" << std::endl;
}
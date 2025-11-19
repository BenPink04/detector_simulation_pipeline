#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>

void plot_10M_acceptance() {
    std::cout << "Creating acceptance histogram for 10M event dataset..." << std::endl;
    
    TFile *file = TFile::Open("/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600_combined_acceptance.root");
    
    TTree *tree = (TTree*)file->Get("kaonEventInfo");
    std::vector<double> *true_mom_vec = nullptr;
    std::vector<int> *reco_flag_vec = nullptr;
    int n_triple_pion_events = 0;
    
    tree->SetBranchAddress("true_mom_vec", &true_mom_vec);
    tree->SetBranchAddress("reco_flag_vec", &reco_flag_vec);
    tree->SetBranchAddress("n_triple_pion_events", &n_triple_pion_events);
    tree->GetEntry(0);
    
    std::cout << "Processing " << n_triple_pion_events << " kaon events..." << std::endl;

    TH1D *h_acceptance = new TH1D("h_acceptance", "Kaon Reconstruction Acceptance (10M Events);True Momentum [GeV/c];Acceptance", 20, 0.0, 10.0);

    std::vector<int> reco_in_bin(20, 0);
    std::vector<int> total_in_bin(20, 0);

    for (int i = 0; i < n_triple_pion_events; ++i) {
        double p = (*true_mom_vec)[i];
        int bin = int(p / 0.5);
        if (bin >= 0 && bin < 20) {
            total_in_bin[bin]++;
            if ((*reco_flag_vec)[i]) reco_in_bin[bin]++;
        }
    }

    for (int bin = 0; bin < 20; ++bin) {
        double acceptance = (total_in_bin[bin] > 0) ? double(reco_in_bin[bin]) / double(total_in_bin[bin]) : 0.0;
        h_acceptance->SetBinContent(bin + 1, acceptance);
    }

    TCanvas *c1 = new TCanvas("c1", "10M Event Acceptance", 800, 600);
    h_acceptance->SetFillColor(kBlue-9);
    h_acceptance->Draw("hist");
    c1->SaveAs("acceptance_10M_events.png");
    
    TFile *outfile = new TFile("acceptance_10M_histogram.root", "RECREATE");
    h_acceptance->Write();
    outfile->Close();
    
    std::cout << "Acceptance histogram created!" << std::endl;
    file->Close();
}

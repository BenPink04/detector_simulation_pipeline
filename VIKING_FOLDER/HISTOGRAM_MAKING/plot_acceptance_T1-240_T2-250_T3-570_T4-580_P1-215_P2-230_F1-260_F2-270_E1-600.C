#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <string>

void plot_acceptance_T1_240_T2_250_T3_570_T4_580_P1_215_P2_230_F1_260_F2_270_E1_600() {
    std::vector<std::string> filenames = {
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_combined_acceptance.root"
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

    TH1D *h_acceptance = new TH1D("h_acceptance", "Kaon Reconstruction Acceptance (T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600);True Momentum [GeV/c];Acceptance", n_bins, p_min, p_max);

    int total_events = 0, total_reconstructed = 0;
    for (int bin = 0; bin < n_bins; ++bin) {
        int reconstructed_vertices = reco_in_bin[bin];
        int total_vertices = total_in_bin[bin];
        total_events += total_vertices;
        total_reconstructed += reconstructed_vertices;
        
        double acceptance = (total_vertices > 0) ? double(reconstructed_vertices) / double(total_vertices) : 0.0;
        h_acceptance->SetBinContent(bin + 1, acceptance);
    }

    TCanvas *c1 = new TCanvas("c1", "Kaon Acceptance T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600", 800, 600);
    h_acceptance->SetFillColor(kBlue-9);
    h_acceptance->SetLineColor(kBlue+2);
    h_acceptance->GetYaxis()->SetRangeUser(0, 1.1);
    h_acceptance->Draw("hist");
    c1->Update();

    gPad->SaveAs("T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_acceptance_plot.png");
    std::cout << "Acceptance plot saved as T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_acceptance_plot.png" << std::endl;
    
    TFile *outfile = new TFile("T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_acceptance_histogram.root", "RECREATE");
    h_acceptance->Write();
    c1->Write();
    outfile->Close();
    
    std::cout << "=== ACCEPTANCE SUMMARY ===" << std::endl;
    std::cout << "Total kaon events: " << total_events << std::endl;
    std::cout << "Total reconstructed: " << total_reconstructed << std::endl;
    std::cout << "Overall acceptance: " << (total_events > 0 ? (double)total_reconstructed/total_events : 0) << std::endl;
}

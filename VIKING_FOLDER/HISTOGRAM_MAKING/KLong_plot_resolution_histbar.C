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
    "/users/bp969/scratch/VIKING_FOLDER/PIPELINE_RESULTS/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600_combined_vectors.root"
    };

    int nbins = 40;
    double min_p = 0.0;
    double max_p = 10.0;
    TH2D *h2 = new TH2D("h2", "Resolution vs True Kaon Momentum;True Kaon Momentum [GeV/c];|Reco - True|/True", nbins, min_p, max_p, nbins, 0, 1);

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
        tree->GetEntry(0);

        for (size_t i = 0; i < reco_p->size(); ++i) {
            double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
            if (res > anomaly_threshold) continue; // Skip anomalous results
            h2->Fill((*true_p)[i], res);
        }
        inFile->Close();
    }

    // Make a 1D histogram of mean resolution per bin
    TH1D *h_mean = new TH1D("h_mean", "Mean Relative Kaon Momentum Resolution;True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)", nbins, min_p, max_p);
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

      // Save the plot as PNG files
    gPad->SaveAs("KLong_histbar_plot.png");
}
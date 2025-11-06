#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>
#include <cmath>

void KLong_plot_resolution_histogram() {
    std::vector<std::string> filenames = {
         "out_seed1_vectors.root",
        "out_seed2_vectors.root",
        "out_seed3_vectors.root",
        "out_seed4_vectors.root",
        "out_seed5_vectors.root",
        "out_seed6_vectors.root",
        "out_seed7_vectors.root",
        "out_seed8_vectors.root",
        "out_seed9_vectors.root",
        "out_seed10_vectors.root"
    };

    // Histogram settings
    int nbins = 40;
    double min_p = 0.0;
    double max_p = 10.0;
    TH2D *h2 = new TH2D("h2", "Resolution vs True Kaon Momentum;True Kaon Momentum [GeV/c];|Reco - True|/True", nbins, min_p, max_p, nbins, 0, 1);

    // Fill histogram with all data
    for (const auto& fname : filenames) {
        TFile *inFile = TFile::Open(fname.c_str());
        if (!inFile || inFile->IsZombie()) {
            std::cout << "Cannot open root file: " << fname << std::endl;
            continue;
        }

        TTree *tree = (TTree*)inFile->Get("kaonVectors");
        if (!tree) {
            std::cout << "Tree kaonVectors not found in " << fname << std::endl;
            inFile->Close();
            continue;
        }

        std::vector<double> *reco_p = nullptr;
        std::vector<double> *true_p = nullptr;

        tree->SetBranchAddress("reco_p", &reco_p);
        tree->SetBranchAddress("true_p", &true_p);

        tree->GetEntry(0);

        for (size_t i = 0; i < reco_p->size(); ++i) {
            double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
            h2->Fill((*true_p)[i], res);
        }

        inFile->Close();
    }

    // Calculate average resolution per bin
    std::vector<double> bin_centers, bin_means;
    for (int i = 1; i <= nbins; ++i) {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
        double center = h2->GetXaxis()->GetBinCenter(i);
        bin_centers.push_back(center);
        bin_means.push_back(mean);
        delete proj;
    }

    // Plot average resolution vs true momentum
    TCanvas *c1 = new TCanvas("c1", "Average Kaon Relative Momentum Resolution", 800, 600);
    TGraph *g_avg = new TGraph(nbins);
    for (int i = 0; i < nbins; ++i)
        g_avg->SetPoint(i, bin_centers[i], bin_means[i]);
    g_avg->SetTitle("Average Relative Kaon Momentum Resolution;True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)");
    g_avg->SetMarkerStyle(20);
    g_avg->SetMarkerColor(kRed);
    g_avg->Draw("AP");

    // Optionally, also draw the 2D histogram
    TCanvas *c2 = new TCanvas("c2", "2D Histogram", 800, 600);
    h2->Draw("COLZ");

    c1->Update();
    c2->Update();
}
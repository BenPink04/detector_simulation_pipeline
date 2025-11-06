#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVector3.h"
#include <vector>
#include <iostream>
#include <cmath>

void KLong_plot_resolution() {
    // List of files to process
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

    std::vector<double> all_true_p, all_resolution;

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
        

        tree->GetEntry(0); // Only one entry, all vectors

        for (size_t i = 0; i < reco_p->size(); ++i) {
            double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
            all_true_p.push_back((*true_p)[i]);
            all_resolution.push_back(res);
        }

        inFile->Close();
    }

    TCanvas *c1 = new TCanvas("c1", "Kaon Relative Momentum Resolution", 800, 600);
    TGraph *g_res = new TGraph(all_resolution.size());
    for (size_t i = 0; i < all_resolution.size(); ++i)
        g_res->SetPoint(i, all_true_p[i], all_resolution[i]);
    g_res->SetTitle("Relative Kaon Momentum Resolution;True Kaon Momentum [GeV/c];|Reco - True|/True");
    g_res->SetMarkerStyle(20);
    g_res->SetMarkerColor(kBlue);
    g_res->Draw("AP");
    c1->Update();
}
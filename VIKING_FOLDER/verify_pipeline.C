void verify_pipeline() {
    // Verify the 10M event pipeline by plotting resolution from combined vectors
    std::string filename = "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600_combined_vectors.root";
    
    TFile *inFile = TFile::Open(filename.c_str());
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    TTree *tree = (TTree*)inFile->Get("vectors");
    if (!tree) {
        std::cerr << "Error: Tree 'vectors' not found in file" << std::endl;
        inFile->Close();
        return;
    }

    // Set up branch addresses
    std::vector<double> *reco_p = nullptr;
    std::vector<double> *true_p = nullptr;
    
    tree->SetBranchAddress("reco_p", &reco_p);
    tree->SetBranchAddress("true_p", &true_p);

    // Parameters for histogram
    int nbins = 20;
    double min_p = 0.0;
    double max_p = 10.0;
    double anomaly_threshold = 1.0; // Skip resolution > 100%

    // Create 2D histogram for resolution vs true momentum
    TH2D *h2 = new TH2D("h2", "Kaon Momentum Resolution;True Kaon Momentum [GeV/c];Relative Resolution |Reco - True|/True", 
                        nbins, min_p, max_p, 100, 0, 1);

    Long64_t nentries = tree->GetEntries();
    std::cout << "Processing " << nentries << " entries from combined vectors file..." << std::endl;

    for (Long64_t entry = 0; entry < nentries; ++entry) {
        tree->GetEntry(entry);
        
        if (!reco_p || !true_p) continue;
        
        for (size_t i = 0; i < reco_p->size() && i < true_p->size(); ++i) {
            if ((*true_p)[i] <= 0) continue; // Skip invalid true momentum
            
            double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
            if (res > anomaly_threshold) continue; // Skip anomalous results
            
            h2->Fill((*true_p)[i], res);
        }
    }

    std::cout << "Total entries in 2D histogram: " << h2->GetEntries() << std::endl;

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
    TCanvas *c1 = new TCanvas("c1", "10M Event Pipeline - Mean Kaon Relative Momentum Resolution", 800, 600);
    h_mean->SetFillColor(kBlue-9);
    h_mean->SetBarWidth(0.9);
    h_mean->SetBarOffset(0.05);
    h_mean->Draw("BAR");
    
    // Add some statistics text
    TPaveText *pt = new TPaveText(0.6, 0.7, 0.85, 0.85, "NDC");
    pt->AddText(Form("Total Events: %lld", nentries));
    pt->AddText(Form("Histogram Entries: %d", (int)h2->GetEntries()));
    pt->AddText("10M Event Pipeline");
    pt->Draw();
    
    c1->Update();

    // Save the plot
    c1->SaveAs("pipeline_verification_histbar.png");
    std::cout << "Verification plot saved as pipeline_verification_histbar.png" << std::endl;

    inFile->Close();
}
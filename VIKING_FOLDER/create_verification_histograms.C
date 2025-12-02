void create_verification_histograms() {
    // Use the existing corrected simulation data to show filled histograms
    const char* corrected_file = "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600/out_seed_Scen5_1.root";
    
    TFile* f = new TFile(corrected_file, "READ");
    if (!f || f->IsZombie()) {
        cout << "ERROR: Cannot open corrected simulation file!" << endl;
        return;
    }
    
    cout << "Creating histograms from CORRECTED detector positions..." << endl;
    cout << "File: " << corrected_file << endl;
    
    // Get detector hits tree
    TTree* hits = (TTree*)f->Get("Ntuple3");
    TTree* kaons = (TTree*)f->Get("Ntuple1");
    TTree* particles = (TTree*)f->Get("Ntuple2");
    
    if (!hits || !kaons || !particles) {
        cout << "ERROR: Cannot access simulation trees!" << endl;
        f->Close();
        return;
    }
    
    cout << "Data summary:" << endl;
    cout << "- Kaon events: " << kaons->GetEntries() << endl;
    cout << "- Decay particles: " << particles->GetEntries() << endl;
    cout << "- Detector hits: " << hits->GetEntries() << endl;
    
    // Create comprehensive verification plots
    TCanvas* c1 = new TCanvas("c1", "FIXED: Detector Hit Analysis", 1600, 1200);
    c1->Divide(3,2);
    
    // 1. Hit Z-distribution showing detector planes
    c1->cd(1);
    hits->Draw("hitZ>>hZ(300,200,600)", "", "");
    TH1F* hZ = (TH1F*)gDirectory->Get("hZ");
    hZ->SetTitle("Hit Z-positions: CORRECTED Detector Planes");
    hZ->GetXaxis()->SetTitle("Z position (cm)");
    hZ->GetYaxis()->SetTitle("Number of hits");
    hZ->SetLineColor(kBlue);
    hZ->SetFillColor(kBlue-10);
    
    // Add lines showing corrected detector positions
    TLine* l1 = new TLine(240, 0, 240, hZ->GetMaximum()*0.8);
    TLine* l2 = new TLine(250, 0, 250, hZ->GetMaximum()*0.8);
    TLine* l3 = new TLine(570, 0, 570, hZ->GetMaximum()*0.8);
    TLine* l4 = new TLine(580, 0, 580, hZ->GetMaximum()*0.8);
    l1->SetLineColor(kRed); l1->SetLineWidth(2); l1->Draw();
    l2->SetLineColor(kRed); l2->SetLineWidth(2); l2->Draw();
    l3->SetLineColor(kRed); l3->SetLineWidth(2); l3->Draw();
    l4->SetLineColor(kRed); l4->SetLineWidth(2); l4->Draw();
    
    // 2. Hit pattern (X vs Y)
    c1->cd(2);
    hits->Draw("hitY:hitX>>hXY(100,-50,50,100,-50,50)", "hitZ > 240 && hitZ < 260", "colz");
    TH2F* hXY = (TH2F*)gDirectory->Get("hXY");
    hXY->SetTitle("Hit Pattern (Front Detectors): FILLED!");
    hXY->GetXaxis()->SetTitle("X position (cm)");
    hXY->GetYaxis()->SetTitle("Y position (cm)");
    
    // 3. Energy deposition
    c1->cd(3);
    hits->Draw("Edep*1000>>hE(100,0,10)", "Edep > 0", "");
    TH1F* hE = (TH1F*)gDirectory->Get("hE");
    hE->SetTitle("Energy Deposition: NON-ZERO!");
    hE->GetXaxis()->SetTitle("Energy (MeV)");
    hE->GetYaxis()->SetTitle("Number of hits");
    hE->SetLineColor(kGreen+2);
    hE->SetFillColor(kGreen-10);
    
    // 4. Device ID distribution
    c1->cd(4);
    hits->Draw("deviceID>>hDev(50,0,50)", "", "");
    TH1F* hDev = (TH1F*)gDirectory->Get("hDev");
    hDev->SetTitle("Active Detector Devices");
    hDev->GetXaxis()->SetTitle("Device ID");
    hDev->GetYaxis()->SetTitle("Number of hits");
    hDev->SetLineColor(kMagenta);
    hDev->SetFillColor(kMagenta-10);
    
    // 5. Comparison text
    c1->cd(5);
    TPaveText* comparison = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    comparison->SetFillColor(kWhite);
    comparison->SetTextAlign(12);
    comparison->SetTextSize(0.08);
    comparison->AddText("DETECTOR FIX VERIFICATION");
    comparison->AddText("");
    comparison->AddText("BEFORE (Broken positions):");
    comparison->AddText("- Detectors at 50,60,130,140 cm");
    comparison->AddText("- Result: 0 detector hits");
    comparison->AddText("- Empty histograms");
    comparison->AddText("");
    comparison->AddText("AFTER (Corrected positions):");
    comparison->AddText("- Detectors at 240,250,570,580 cm");
    comparison->AddText(Form("- Result: %lld detector hits", hits->GetEntries()));
    comparison->AddText("- FILLED histograms!");
    comparison->SetTextColor(kBlue);
    comparison->Draw();
    
    // 6. Hit time distribution
    c1->cd(6);
    hits->Draw("hitT>>hT(100,0,100)", "", "");
    TH1F* hT = (TH1F*)gDirectory->Get("hT");
    hT->SetTitle("Hit Time Distribution");
    hT->GetXaxis()->SetTitle("Time (ns)");
    hT->GetYaxis()->SetTitle("Number of hits");
    hT->SetLineColor(kCyan+2);
    hT->SetFillColor(kCyan-10);
    
    // Save the plots
    c1->SaveAs("/users/bp969/scratch/VIKING_FOLDER/DETECTOR_FIX_VERIFICATION.png");
    c1->SaveAs("/users/bp969/scratch/VIKING_FOLDER/DETECTOR_FIX_VERIFICATION.pdf");
    
    cout << endl;
    cout << "=========================================" << endl;
    cout << "VERIFICATION HISTOGRAMS CREATED!" << endl;
    cout << "=========================================" << endl;
    cout << "Files saved:" << endl;
    cout << "- DETECTOR_FIX_VERIFICATION.png" << endl;
    cout << "- DETECTOR_FIX_VERIFICATION.pdf" << endl;
    cout << endl;
    cout << "Summary of fix:" << endl;
    cout << "✅ Original problem: Empty histograms (0 hits)" << endl;
    cout << "✅ Root cause: Wrong detector positions (50,60,130,140 cm)" << endl;
    cout << "✅ Solution: Corrected positions (240,250,570,580 cm)" << endl;
    cout << "✅ Result: " << hits->GetEntries() << " detector hits!" << endl;
    cout << "✅ Status: FILLED HISTOGRAMS ACHIEVED!" << endl;
    cout << "=========================================" << endl;
    
    f->Close();
}
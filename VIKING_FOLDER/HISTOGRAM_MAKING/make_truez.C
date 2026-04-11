// make_truez.C
// Standalone ROOT macro — run as:
//   root -l -b -q 'make_truez.C("path/to/_combined_vectors.root", "path/to/_combined_vectors_truez.root")'
//
// Reads reco_p_truez from the input and writes a new kaonVectors tree where
//   reco_p       = reco_p_truez   (truth-vertex reconstruction)
//   reco_p_poca  = reco_p         (original PoCA momentum, preserved)
// All other branches are copied unchanged.

void make_truez(const char* inPath, const char* outPath) {
    TFile *fin = TFile::Open(inPath);
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open " << inPath << "\n"; return; }

    TTree *tin = (TTree*)fin->Get("kaonVectors");
    if (!tin) {
        std::cerr << "ERROR: kaonVectors tree not found in " << inPath << "\n";
        fin->Close(); return; }

    if (!tin->GetBranch("reco_p_truez")) {
        std::cerr << "SKIP: reco_p_truez branch missing in " << inPath
                  << " — rerun pipeline with updated KLong_save_vectors.C first\n";
        fin->Close(); return; }

    TFile *fout = new TFile(outPath, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot create " << outPath << "\n";
        fin->Close(); return; }

    TTree *tout = new TTree("kaonVectors",
        "Kaon vectors - truth-vertex diagnostic (reco_p = reco_p_truez)");

    std::vector<double> *reco_p       = nullptr;
    std::vector<double> *true_p       = nullptr;
    std::vector<double> *reco_p_truez = nullptr;
    std::vector<double> *reco_p_truet = nullptr;
    std::vector<double> *rvx = nullptr, *rvy = nullptr, *rvz = nullptr;
    std::vector<double> *tvx = nullptr, *tvy = nullptr, *tvz = nullptr;

    tin->SetBranchAddress("reco_p",        &reco_p);
    tin->SetBranchAddress("true_p",        &true_p);
    tin->SetBranchAddress("reco_p_truez",  &reco_p_truez);
    tin->SetBranchAddress("reco_p_truet",  &reco_p_truet);
    tin->SetBranchAddress("reco_vertex_x", &rvx);
    tin->SetBranchAddress("reco_vertex_y", &rvy);
    tin->SetBranchAddress("reco_vertex_z", &rvz);
    tin->SetBranchAddress("true_vertex_x", &tvx);
    tin->SetBranchAddress("true_vertex_y", &tvy);
    tin->SetBranchAddress("true_vertex_z", &tvz);

    std::vector<double> out_reco_p, out_reco_p_poca, out_reco_p_truet;
    std::vector<double> out_true_p;
    std::vector<double> out_rvx, out_rvy, out_rvz;
    std::vector<double> out_tvx, out_tvy, out_tvz;

    tout->Branch("reco_p",        &out_reco_p);
    tout->Branch("reco_p_poca",   &out_reco_p_poca);
    tout->Branch("reco_p_truet",  &out_reco_p_truet);
    tout->Branch("true_p",        &out_true_p);
    tout->Branch("reco_vertex_x", &out_rvx);
    tout->Branch("reco_vertex_y", &out_rvy);
    tout->Branch("reco_vertex_z", &out_rvz);
    tout->Branch("true_vertex_x", &out_tvx);
    tout->Branch("true_vertex_y", &out_tvy);
    tout->Branch("true_vertex_z", &out_tvz);

    Long64_t n = tin->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tin->GetEntry(i);
        out_reco_p       = *reco_p_truez;
        out_reco_p_poca  = *reco_p;
        out_reco_p_truet = *reco_p_truet;
        out_true_p       = *true_p;
        out_rvx = *rvx; out_rvy = *rvy; out_rvz = *rvz;
        out_tvx = *tvx; out_tvy = *tvy; out_tvz = *tvz;
        tout->Fill();
    }

    fout->Write();
    fout->Close();
    fin->Close();
    std::cout << "Created " << outPath << " (" << n << " entries)\n";
}

void pdg_check() {
    TFile *f = TFile::Open("/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700_100.root");
    if(!f || f->IsZombie()){ printf("Cannot open file\n"); return; }
    TTree *t3 = (TTree*)f->Get("Ntuple3");
    if(!t3){ printf("No Ntuple3\n"); return; }
    Double_t ev3, dev3, pdg3, ht, hz;
    t3->SetBranchAddress("evtNb",       &ev3);
    t3->SetBranchAddress("deviceID",    &dev3);
    t3->SetBranchAddress("PDGEncoding", &pdg3);
    t3->SetBranchAddress("hitT",        &ht);
    t3->SetBranchAddress("hitZ",        &hz);
    int n_pz=0, n_tf=0, cpp=0, ctf=0;
    std::map<int,int> pdgcount;
    std::set<int> devset;
    for(Long64_t i=0; i<t3->GetEntries(); i++){
        t3->GetEntry(i);
        pdgcount[(int)pdg3]++;
        if((int)dev3>=1900 && (int)dev3<=2100) devset.insert((int)dev3);
        if((int)dev3>=1953 && (int)dev3<=2000){
            n_pz++;
            if(cpp<5){ printf("PIZZA evt=%d dev=%d pdg=%d z=%.2f t=%.6f\n",(int)ev3,(int)dev3,(int)pdg3,hz,ht); cpp++; }
        }
        if((int)dev3>=2053 && (int)dev3<=2070){
            n_tf++;
            if(ctf<5){ printf("TOF   evt=%d dev=%d pdg=%d z=%.2f t=%.6f\n",(int)ev3,(int)dev3,(int)pdg3,hz,ht); ctf++; }
        }
    }
    printf("n_pizza=%d  n_tof=%d\n",n_pz,n_tf);
    printf("Device IDs in range 1900-2100:");
    for(int d:devset) printf(" %d",d);
    printf("\n");
    printf("PDG code counts (top 20):\n");
    int shown=0;
    for(auto &kv:pdgcount){ if(shown++<20) printf("  pdg=%-8d  count=%d\n",kv.first,kv.second); }
    TTree *t2=(TTree*)f->Get("Ntuple2");
    if(!t2){ printf("No Ntuple2\n"); f->Close(); return; }
    Double_t ev2, pdg2, ppx, ppy, ppz;
    t2->SetBranchAddress("evtNb",                  &ev2);
    t2->SetBranchAddress("DKparticle_PDGEncoding", &pdg2);
    t2->SetBranchAddress("DKparticle_momX",        &ppx);
    t2->SetBranchAddress("DKparticle_momY",        &ppy);
    t2->SetBranchAddress("DKparticle_momZ",        &ppz);
    printf("Ntuple2 first 10 entries:\n");
    for(Long64_t i=0; i<10; i++){
        t2->GetEntry(i);
        printf("  evt=%d pdg=%d p=(%.4f,%.4f,%.4f)\n",(int)ev2,(int)pdg2,ppx,ppy,ppz);
    }
    f->Close();
}

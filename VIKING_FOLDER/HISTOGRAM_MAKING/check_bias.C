void check_bias() {
    const char* file = "/users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/ToF_FIX_1_20260404/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700_combined_vectors.root";
    TFile *f = TFile::Open(file);
    TTree *t = (TTree*)f->Get("kaonVectors");
    std::vector<double> *rp=nullptr, *tp=nullptr, *rptz=nullptr, *rptt=nullptr;
    t->SetBranchAddress("reco_p",       &rp);
    t->SetBranchAddress("true_p",       &tp);
    t->SetBranchAddress("reco_p_truez", &rptz);
    t->SetBranchAddress("reco_p_truet", &rptt);

    const int nb = 6;
    double edges[] = {0.5,1.5,2.5,3.5,5.0,7.0,11.0};
    double sr[nb]={}, stz[nb]={}, stt[nb]={};
    double sr_c[nb]={}, stz_c[nb]={};
    double n[nb]={}, nc[nb]={}, ntt[nb]={};

    for (Long64_t i=0; i<t->GetEntries(); i++) {
        t->GetEntry(i);
        for (size_t j=0; j<rp->size(); j++) {
            double r  = (*rp)[j];
            double tr = (*tp)[j];
            double tz = (*rptz)[j];
            double tt = (*rptt)[j];
            if (tr<=0 || tz<=0 || !std::isfinite(r)) continue;
            int bin=-1;
            for (int b=0; b<nb; b++) if (tr>=edges[b] && tr<edges[b+1]) { bin=b; break; }
            if (bin<0) continue;
            sr[bin]  += (r -tr)/tr;
            stz[bin] += (tz-tr)/tr;
            n[bin]++;
            if (tt>0) { stt[bin]+=(tt-tr)/tr; ntt[bin]++; }
            if (fabs((tz-tr)/tr) < 0.4) {
                sr_c[bin]  += (r -tr)/tr;
                stz_c[bin] += (tz-tr)/tr;
                nc[bin]++;
            }
        }
    }
    printf("\n%-14s %5s  %8s %8s %8s  |  %8s %8s  %5s\n",
           "p_true range","N","std","truez","truet","cut_std","cut_tz","N_cut");
    for (int b=0; b<nb; b++)
        printf("[%4.1f,%4.1f)  %5d  %8.4f %8.4f %8.4f  |  %8.4f %8.4f  %5d\n",
               edges[b],edges[b+1],(int)n[b],
               n[b]>0?sr[b]/n[b]:0., n[b]>0?stz[b]/n[b]:0.,
               ntt[b]>0?stt[b]/ntt[b]:0.,
               nc[b]>0?sr_c[b]/nc[b]:0., nc[b]>0?stz_c[b]/nc[b]:0., (int)nc[b]);
    f->Close();
}

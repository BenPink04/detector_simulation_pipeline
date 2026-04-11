void debug_reco() {
    const char* file = "/users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/ToF_FIX_1_20260404/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700_combined_vectors.root";
    TFile *f = TFile::Open(file);
    TTree *t = (TTree*)f->Get("kaonVectors");

    std::vector<double> *rp=nullptr, *tp=nullptr, *rptz=nullptr;
    std::vector<double> *rvx=nullptr, *rvy=nullptr, *rvz=nullptr;
    std::vector<double> *tvx=nullptr, *tvy=nullptr, *tvz=nullptr;
    t->SetBranchAddress("reco_p",       &rp);
    t->SetBranchAddress("true_p",       &tp);
    t->SetBranchAddress("reco_p_truez", &rptz);
    t->SetBranchAddress("reco_vertex_x",&rvx);
    t->SetBranchAddress("reco_vertex_y",&rvy);
    t->SetBranchAddress("reco_vertex_z",&rvz);
    t->SetBranchAddress("true_vertex_x",&tvx);
    t->SetBranchAddress("true_vertex_y",&tvy);
    t->SetBranchAddress("true_vertex_z",&tvz);

    int printed=0;
    const double c_cm_ns = 29.9792458;  // cm/ns
    const double mK = 0.497611;         // GeV/c^2

    for (Long64_t i=0; i<t->GetEntries() && printed<10; i++) {
        t->GetEntry(i);
        for (size_t j=0; j<rp->size() && printed<10; j++) {
            double r = (*rp)[j], tr = (*tp)[j], tz = (*rptz)[j];
            double rx=(*rvx)[j],ry=(*rvy)[j],rz=(*rvz)[j];
            double tx=(*tvx)[j],ty=(*tvy)[j],tzv=(*tvz)[j];
            if (tr <= 0.) continue;

            // reco flight = distance from origin to reco vertex (what standard uses)
            double flight_reco = sqrt(rx*rx+ry*ry+rz*rz);
            // true flight = distance from origin to true vertex
            double flight_true = sqrt(tx*tx+ty*ty+tzv*tzv);

            // What kaon_v truez must have used: flight_true / t_decay_truez
            // We can back out t_decay_truez from kaon_p_truez:
            // tz = gamma*mK*beta => tz/mK = gamma*beta = p/m => beta = tz/sqrt(tz^2+mK^2)
            // v_tz = beta * c
            double beta_tz_expect = tr / sqrt(tr*tr + mK*mK);  // from true_p
            double t_decay_true = (flight_true*1e-2) / (beta_tz_expect * 2.998e8);  // s

            // What truez actually produces: must have used t_decay such that
            // v_tz = flight_true / t_decay_tz => beta_tz = v_tz/c
            // tz = gamma*mK*beta_tz
            double beta_tz_reco = tz > 0 ? tz / sqrt(tz*tz + mK*mK) : -1.;
            double t_decay_tz_implied = (tz > 0) ? (flight_true*1e-2) / (beta_tz_reco * 2.998e8) : -1.;

            printf("Evt %d/%zu: true_p=%.3f tz=%.3f r=%.3f\n"
                   "  true_vtx=(%.1f,%.1f,%.1f) reco_vtx=(%.1f,%.1f,%.1f)\n"
                   "  flight_true=%.2f cm  flight_reco=%.2f cm\n"
                   "  t_decay_true=%.4f ns  t_decay_tz_implied=%.4f ns\n"
                   "  dp/p_tz=%.3f\n",
                   (int)i, j, tr, tz, r,
                   tx,ty,tzv, rx,ry,rz,
                   flight_true, flight_reco,
                   t_decay_true*1e9, t_decay_tz_implied*1e9,
                   (tz-tr)/tr);
            printed++;
        }
    }
    f->Close();
}

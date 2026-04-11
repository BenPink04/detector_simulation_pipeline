void inspect_ntuple() {
    const double mpi  = 0.13957;
    const double c_cms = 29.9792458;
    const char* raw = "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700_100.root";
    TFile *f = TFile::Open(raw);

    // Read Ntuple1: kaon truth
    TTree *nt1 = (TTree*)f->Get("Ntuple1");
    Int_t ev1; Double_t kpx,kpy,kpz;
    nt1->SetBranchAddress("evtNb",       &ev1);
    nt1->SetBranchAddress("kaonDK_momX", &kpx);
    nt1->SetBranchAddress("kaonDK_momY", &kpy);
    nt1->SetBranchAddress("kaonDK_momZ", &kpz);
    std::map<int,double> kaon_p;
    for (Long64_t i=0; i<nt1->GetEntries(); i++) {
        nt1->GetEntry(i);
        kaon_p[ev1] = sqrt(kpx*kpx+kpy*kpy+kpz*kpz);
    }

    // Read Ntuple2: true pion momenta
    TTree *nt2 = (TTree*)f->Get("Ntuple2");
    Int_t ev2, pdg2; Double_t ppx,ppy,ppz;
    nt2->SetBranchAddress("evtNb",               &ev2);
    nt2->SetBranchAddress("DKparticle_PDGEncoding",&pdg2);
    nt2->SetBranchAddress("DKparticle_momX",      &ppx);
    nt2->SetBranchAddress("DKparticle_momY",      &ppy);
    nt2->SetBranchAddress("DKparticle_momZ",      &ppz);
    std::map<int,double> vpip_true, vpim_true;  // evt -> v_z_true (cm/ns)
    for (Long64_t i=0; i<nt2->GetEntries(); i++) {
        nt2->GetEntry(i);
        double p = sqrt(ppx*ppx+ppy*ppy+ppz*ppz);
        double e = sqrt(p*p + mpi*mpi);
        double vz = (ppz/e) * c_cms;
        if (pdg2 == 211)  vpip_true[ev2] = vz;
        if (pdg2 == -211) vpim_true[ev2] = vz;
    }

    // Read Ntuple3: hits at pizza and tof
    TTree *nt3 = (TTree*)f->Get("Ntuple3");
    Int_t ev3, dev3, pdg3; Double_t hz, ht;
    nt3->SetBranchAddress("evtNb",       &ev3);
    nt3->SetBranchAddress("deviceID",    &dev3);
    nt3->SetBranchAddress("PDGEncoding", &pdg3);
    nt3->SetBranchAddress("hitZ",        &hz);
    nt3->SetBranchAddress("hitT",        &ht);

    // per-event: store min-time pizza and tof hits
    std::map<int,double> pip_t_pz, pip_z_pz, pip_t_tf, pip_z_tf;
    std::map<int,double> pim_t_pz, pim_z_pz, pim_t_tf, pim_z_tf;
    for (Long64_t i=0; i<nt3->GetEntries(); i++) {
        nt3->GetEntry(i);
        bool is_pz = (dev3>=1953 && dev3<=2000);
        bool is_tf = (dev3>=2053 && dev3<=2070);
        if (!is_pz && !is_tf) continue;
        if (pdg3 == 211 && is_pz) {
            if (!pip_t_pz.count(ev3) || ht<pip_t_pz[ev3]) { pip_t_pz[ev3]=ht; pip_z_pz[ev3]=hz; }
        }
        if (pdg3 == 211 && is_tf) {
            if (!pip_t_tf.count(ev3) || ht<pip_t_tf[ev3]) { pip_t_tf[ev3]=ht; pip_z_tf[ev3]=hz; }
        }
        if (pdg3 == -211 && is_pz) {
            if (!pim_t_pz.count(ev3) || ht<pim_t_pz[ev3]) { pim_t_pz[ev3]=ht; pim_z_pz[ev3]=hz; }
        }
        if (pdg3 == -211 && is_tf) {
            if (!pim_t_tf.count(ev3) || ht<pim_t_tf[ev3]) { pim_t_tf[ev3]=ht; pim_z_tf[ev3]=hz; }
        }
    }

    // Compare
    int printed=0; double sfp=0, sfm=0; int ng=0;
    for (auto& kv : pip_t_pz) {
        int evt = kv.first;
        if (!pip_t_tf.count(evt)||!pim_t_pz.count(evt)||!pim_t_tf.count(evt)) continue;
        if (!vpip_true.count(evt)||!vpim_true.count(evt)||!kaon_p.count(evt)) continue;
        double pip_dt = pip_t_tf[evt] - pip_t_pz[evt];  // ns
        double pip_dz = pip_z_tf[evt] - pip_z_pz[evt];  // cm
        double pim_dt = pim_t_tf[evt] - pim_t_pz[evt];
        double pim_dz = pim_z_tf[evt] - pim_z_pz[evt];
        if (pip_dt<=0||pim_dt<=0) continue;
        double pip_vz_reco = pip_dz / pip_dt;  // cm/ns, z-component proxy
        double pim_vz_reco = pim_dz / pim_dt;
        double pip_frac = (pip_vz_reco - vpip_true[evt]) / vpip_true[evt];
        double pim_frac = (pim_vz_reco - vpim_true[evt]) / vpim_true[evt];
        sfp += pip_frac; sfm += pim_frac; ng++;
        if (printed<6) {
            printf("Evt%4d K_p=%.2f:\n"
                   "  pip: dz/dt=%.3f  vz_true=%.3f  frac=%+.4f  piz_z=%.2f tof_z=%.2f\n"
                   "  pim: dz/dt=%.3f  vz_true=%.3f  frac=%+.4f\n",
                   evt, kaon_p[evt],
                   pip_vz_reco, vpip_true[evt], pip_frac, pip_z_pz[evt], pip_z_tf[evt],
                   pim_vz_reco, vpim_true[evt], pim_frac);
            printed++;
        }
    }
    printf("\nStats over %d events:\n  pip_frac mean=%+.4f\n  pim_frac mean=%+.4f\n",
           ng, sfp/ng, sfm/ng);
    f->Close();
}
//   Ntuple1: kaon truth (vx,vy,vz, momentum, decay time)
//   Ntuple2: daughter pion truth momenta
//   Ntuple3: detector hits (pizza time, tof time, positions)
// Then computes:
//   v_pi_truth  = |p_pi| / E_pi * c  (from Ntuple2)
//   pip_v_reco  = (path pizza->tof) / (t_tof - t_pizza)  (from Ntuple3, same as code)

void inspect_ntuple() {
    const double mpi  = 0.13957; // GeV/c^2
    const double c_cms = 29.9792458; // cm/ns
    const char* raw = "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700_100.root";
    TFile *f = TFile::Open(raw);

    // --- Read Ntuple1 (kaon truth) ---
    TTree *nt1 = (TTree*)f->Get("Ntuple1");
    Int_t evt1; Double_t kpx,kpy,kpz,kvx,kvy,kvz,kdt;
    nt1->SetBranchAddress("evtNb",       &evt1);
    nt1->SetBranchAddress("kaonDK_momX", &kpx);
    nt1->SetBranchAddress("kaonDK_momY", &kpy);
    nt1->SetBranchAddress("kaonDK_momZ", &kpz);
    nt1->SetBranchAddress("kaonDK_posX", &kvx);
    nt1->SetBranchAddress("kaonDK_posY", &kvy);
    nt1->SetBranchAddress("kaonDK_posZ", &kvz);
    nt1->SetBranchAddress("kaonDK_time", &kdt);  // s -> need to convert
    std::map<int,std::tuple<double,double,double,double,double,double,double>> kaon_truth;
    for (Long64_t i=0; i<nt1->GetEntries(); i++) {
        nt1->GetEntry(i);
        kaon_truth[evt1] = {kpx,kpy,kpz,kvx,kvy,kvz,kdt};
    }

    // --- Read Ntuple2 (pion truth momenta) ---
    TTree *nt2 = (TTree*)f->Get("Ntuple2");
    Int_t evt2, pdg2; Double_t ppx,ppy,ppz;
    nt2->SetBranchAddress("evtNb",               &evt2);
    nt2->SetBranchAddress("DKparticle_PDGEncoding",&pdg2);
    nt2->SetBranchAddress("DKparticle_momX",      &ppx);
    nt2->SetBranchAddress("DKparticle_momY",      &ppy);
    nt2->SetBranchAddress("DKparticle_momZ",      &ppz);
    std::map<int,std::pair<double,double>> pion_v_true;  // evt -> (v_pip_truth, v_pim_truth)
    std::map<int,std::tuple<double,double,double>> pip_mom_true, pim_mom_true;
    for (Long64_t i=0; i<nt2->GetEntries(); i++) {
        nt2->GetEntry(i);
        double p = sqrt(ppx*ppx+ppy*ppy+ppz*ppz);
        double e = sqrt(p*p + mpi*mpi);
        double v = (p/e) * c_cms;  // in cm/ns
        if (pdg2 == 211)  pip_mom_true[evt2] = {ppx,ppy,ppz};
        if (pdg2 == -211) pim_mom_true[evt2] = {ppx,ppy,ppz};
    }

    // --- Read Ntuple3 (hit data) — collect pizza and tof times per event ---
    TTree *nt3 = (TTree*)f->Get("Ntuple3");
    Int_t evt3, dev3, pdg3; Double_t hx,hy,hz,ht;
    nt3->SetBranchAddress("evtNb",       &evt3);
    nt3->SetBranchAddress("deviceID",    &dev3);
    nt3->SetBranchAddress("PDGEncoding", &pdg3);
    nt3->SetBranchAddress("hitX",        &hx);
    nt3->SetBranchAddress("hitY",        &hy);
    nt3->SetBranchAddress("hitZ",        &hz);
    nt3->SetBranchAddress("hitT",        &ht);  // in ns

    // pizza IDs: 1953-2000; tof IDs: 2053-2070
    const double z_tof = 700.;  // cm
    struct HitSummary {
        double pip_t_pizza=-1, pip_z_pizza=0, pip_t_tof=-1;
        double pim_t_pizza=-1, pim_z_pizza=0, pim_t_tof=-1;
        double pip_z_tof=z_tof, pim_z_tof=z_tof;
    };
    std::map<int,HitSummary> hits;
    for (Long64_t i=0; i<nt3->GetEntries(); i++) {
        nt3->GetEntry(i);
        bool is_pizza = (dev3>=1953 && dev3<=2000);
        bool is_tof   = (dev3>=2053 && dev3<=2070);
        if (!is_pizza && !is_tof) continue;
        auto& h = hits[evt3];
        if (pdg3 == 211 && is_pizza)  { if(h.pip_t_pizza<0||ht<h.pip_t_pizza){h.pip_t_pizza=ht; h.pip_z_pizza=hz;} }
        if (pdg3 == 211 && is_tof)    { if(h.pip_t_tof<0  ||ht<h.pip_t_tof)  {h.pip_t_tof=ht;  h.pip_z_tof=hz;}  }
        if (pdg3 == -211 && is_pizza) { if(h.pim_t_pizza<0||ht<h.pim_t_pizza){h.pim_t_pizza=ht; h.pim_z_pizza=hz;} }
        if (pdg3 == -211 && is_tof)   { if(h.pim_t_tof<0  ||ht<h.pim_t_tof)  {h.pim_t_tof=ht;  h.pim_z_tof=hz;}  }
    }

    // --- Compare ---
    int printed = 0;
    double sum_frac_pip=0, sum_frac_pim=0; int n_good=0;
    for (auto& [evt, hs] : hits) {
        if (hs.pip_t_pizza<0 || hs.pip_t_tof<0) continue;
        if (hs.pim_t_pizza<0 || hs.pim_t_tof<0) continue;
        auto km = kaon_truth.find(evt);
        if (km == kaon_truth.end()) continue;
        auto pip_m = pip_mom_true.find(evt);
        if (pip_m == pip_mom_true.end()) continue;
        auto pim_m = pim_mom_true.find(evt);
        if (pim_m == pim_mom_true.end()) continue;

        // True pion velocities (in cm/ns)
        auto [ppx_t,ppy_t,ppz_t] = pip_m->second;
        double pp_t = sqrt(ppx_t*ppx_t+ppy_t*ppy_t+ppz_t*ppz_t);
        double v_pip_true = pp_t / sqrt(pp_t*pp_t+mpi*mpi) * c_cms;
        auto [mmpx_t,mmpy_t,mmpz_t] = pim_m->second;
        double mp_t = sqrt(mmpx_t*mmpx_t+mmpy_t*mmpy_t+mmpz_t*mmpz_t);
        double v_pim_true = mp_t / sqrt(mp_t*mp_t+mpi*mpi) * c_cms;

        // Reconstructed pion velocities: path from pizza_z to z_tof along track fit
        // BUT here we just use the raw hit positions (no track fit)
        // to check if the basic delta_z / delta_t gives the right velocity
        double pip_dz = hs.pip_z_tof - hs.pip_z_pizza;  // cm, z-only distance
        double pip_dt = hs.pip_t_tof - hs.pip_t_pizza;   // ns
        double pip_v_zdz = (pip_dt>0) ? pip_dz / pip_dt : -1;  // cm/ns, z-component only

        double pim_dz = hs.pim_z_tof - hs.pim_z_pizza;
        double pim_dt = hs.pim_t_tof - hs.pim_t_pizza;
        double pim_v_zdz = (pim_dt>0) ? pim_dz / pim_dt : -1;

        // True velocity z-component
        double v_pip_z_true = ppz_t/sqrt(pp_t*pp_t+mpi*mpi) * c_cms;
        double v_pim_z_true = mmpz_t/sqrt(mp_t*mp_t+mpi*mpi) * c_cms;

        if (pip_dt > 0 && pim_dt > 0) {
            sum_frac_pip += (pip_v_zdz - v_pip_z_true) / v_pip_z_true;
            sum_frac_pim += (pim_v_zdz - v_pim_z_true) / v_pim_z_true;
            n_good++;
        }
        if (printed < 8) {
            auto [kpx2,kpy2,kpz2,kvx2,kvy2,kvz2,kdt2] = km->second;
            double km_p = sqrt(kpx2*kpx2+kpy2*kpy2+kpz2*kpz2);
            printf("Evt%4d: K_p=%.2f  piz_z_pizza=%.1f piz_z_tof=%.1f\n"
                   "  pip: v_true=%.4f  v_z_reco=%.4f  dt=%.4f ns  dz=%.1f cm  frac=(%.4f)\n"
                   "  pim: v_true=%.4f  v_z_reco=%.4f  dt=%.4f ns  dz=%.1f cm  frac=(%.4f)\n",
                   evt, km_p, hs.pip_z_pizza, hs.pip_z_tof,
                   v_pip_true, pip_v_zdz, pip_dt, pip_dz, (pip_v_zdz-v_pip_z_true)/v_pip_z_true,
                   v_pim_true, pim_v_zdz, pim_dt, pim_dz, (pim_v_zdz-v_pim_z_true)/v_pim_z_true);
            printed++;
        }
    }
    if (n_good>0) {
        printf("\nMean fractional bias on v_z (= dz/dt vs v_z_true):\n");
        printf("  pip: %.4f (N=%d)\n", sum_frac_pip/n_good, n_good);
        printf("  pim: %.4f (N=%d)\n", sum_frac_pim/n_good, n_good);
    }
    f->Close();
}

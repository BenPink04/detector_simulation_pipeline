void pion_v_check() {
    const double mpi   = 0.13957;
    const double c_cns = 29.9792458; // cm/ns
    const char* raw = "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700/T1-400_T2-410_T3-680_T4-690_P1-375_P2-390_F1-420_F2-430_E1-700_100.root";
    TFile *f = TFile::Open(raw);

    // kaon truth
    TTree *nt1 = (TTree*)f->Get("Ntuple1");
    Double_t ev1; Double_t kpx,kpy,kpz;
    nt1->SetBranchAddress("evtNb",&ev1);
    nt1->SetBranchAddress("kaonDK_momX",&kpx);
    nt1->SetBranchAddress("kaonDK_momY",&kpy);
    nt1->SetBranchAddress("kaonDK_momZ",&kpz);
    std::map<int,double> kp;
    for(Long64_t i=0;i<nt1->GetEntries();i++){nt1->GetEntry(i);kp[(int)ev1]=sqrt(kpx*kpx+kpy*kpy+kpz*kpz);}

    // pion truth momenta
    TTree *nt2 = (TTree*)f->Get("Ntuple2");
    Double_t ev2,pdg2; Double_t ppx,ppy,ppz;
    nt2->SetBranchAddress("evtNb",&ev2);
    nt2->SetBranchAddress("DKparticle_PDGEncoding",&pdg2);
    nt2->SetBranchAddress("DKparticle_momX",&ppx);
    nt2->SetBranchAddress("DKparticle_momY",&ppy);
    nt2->SetBranchAddress("DKparticle_momZ",&ppz);
    std::map<int,double> vpip,vpim, vpip_tot,vpim_tot;
    for(Long64_t i=0;i<nt2->GetEntries();i++){
        nt2->GetEntry(i);
        double p=sqrt(ppx*ppx+ppy*ppy+ppz*ppz),e=sqrt(p*p+mpi*mpi);
        double vz=(ppz/e)*c_cns;
        double vtot=(p/e)*c_cns;
        if((int)pdg2==211)  { vpip[(int)ev2]=vz;  vpip_tot[(int)ev2]=vtot; }
        if((int)pdg2==-211) { vpim[(int)ev2]=vz;  vpim_tot[(int)ev2]=vtot; }
    }

    // hit data (3D)
    TTree *nt3 = (TTree*)f->Get("Ntuple3");
    Double_t ev3,dev3,pdg3; Double_t hx,hy,hz,ht;
    nt3->SetBranchAddress("evtNb",&ev3);
    nt3->SetBranchAddress("deviceID",&dev3);
    nt3->SetBranchAddress("PDGEncoding",&pdg3);
    nt3->SetBranchAddress("hitX",&hx);
    nt3->SetBranchAddress("hitY",&hy);
    nt3->SetBranchAddress("hitZ",&hz);
    nt3->SetBranchAddress("hitT",&ht);

    // pip: pizza time/x/y/z and tof time/x/y/z
    std::map<int,double> ptpz,pxpz,pypz,pzzpz,pttf,pxtf,pytf,pztf;
    // pim
    std::map<int,double> mtpz,mxpz,mypz,mzzpz,mttf,mxtf,mytf,mztf;
    for(Long64_t i=0;i<nt3->GetEntries();i++){
        nt3->GetEntry(i);
        int idev=(int)dev3, ipdg=(int)pdg3, iev=(int)ev3;
        bool ispz=(idev>=1953&&idev<=2000),istf=(idev>=2053&&idev<=2070);
        if(!ispz&&!istf) continue;
        if(ipdg==211&&ispz){if(!ptpz.count(iev)||ht<ptpz[iev]){ptpz[iev]=ht;pxpz[iev]=hx;pypz[iev]=hy;pzzpz[iev]=hz;}}
        if(ipdg==211&&istf){if(!pttf.count(iev)||ht<pttf[iev]){pttf[iev]=ht;pxtf[iev]=hx;pytf[iev]=hy;pztf[iev]=hz;}}
        if(ipdg==-211&&ispz){if(!mtpz.count(iev)||ht<mtpz[iev]){mtpz[iev]=ht;mxpz[iev]=hx;mypz[iev]=hy;mzzpz[iev]=hz;}}
        if(ipdg==-211&&istf){if(!mttf.count(iev)||ht<mttf[iev]){mttf[iev]=ht;mxtf[iev]=hx;mytf[iev]=hy;mztf[iev]=hz;}}
    }

    int pr=0;
    double sfp_z=0,sfm_z=0, sfp_3d=0,sfm_3d=0;
    int ng=0;
    for(auto& kv:ptpz){
        int evt=kv.first;
        if(!pttf.count(evt)||!mtpz.count(evt)||!mttf.count(evt)) continue;
        if(!vpip.count(evt)||!vpim.count(evt)||!kp.count(evt)) continue;
        if(!vpip_tot.count(evt)||!vpim_tot.count(evt)) continue;
        double pdz=pztf[evt]-pzzpz[evt], pdt=pttf[evt]-ptpz[evt];
        double mdz=mztf[evt]-mzzpz[evt], mdt=mttf[evt]-mtpz[evt];
        if(pdt<=0||mdt<=0) continue;
        // z-only velocity
        double pvz=pdz/pdt, mvz=mdz/mdt;
        // 3D velocity from raw hits
        double pdx=pxtf[evt]-pxpz[evt], pdy=pytf[evt]-pypz[evt];
        double mdx=mxtf[evt]-mxpz[evt], mdy=mytf[evt]-mypz[evt];
        double pv3d = sqrt(pdx*pdx+pdy*pdy+pdz*pdz)/pdt;
        double mv3d = sqrt(mdx*mdx+mdy*mdy+mdz*mdz)/mdt;
        // true total speed (cm/ns)
        double vp_true = c_cns*sqrt(vpip[evt]*vpip[evt])/(c_cns); // vz only — get true total
        // recompute true total speed from true momentum
        // (vpip stored as vz; need v_total from ntuple2 — stored as (pz/E)*c so not total)
        // Use the stored vz as proxy and also show 3D raw / vz_true
        double pfz = (pvz - vpip[evt])/fabs(vpip[evt]);
        double mfz = (mvz - vpim[evt])/fabs(vpim[evt]);
        sfp_z+=pfz; sfm_z+=mfz;
        // 3D raw hit speed vs true total speed
        double pf3d = (pv3d - vpip_tot[evt])/vpip_tot[evt];
        double mf3d = (mv3d - vpim_tot[evt])/vpim_tot[evt];
        sfp_3d+=pf3d; sfm_3d+=mf3d;
        ng++;
        if(pr<6){
            printf("Evt%4d Kp=%.2f: pip z=%.1f->%.1f dt=%.4f vz=%.3f(tru_z=%.3f,tru_tot=%.3f) frac_3d=%+.4f\n"
                   "              pim z=%.1f->%.1f dt=%.4f vz=%.3f(tru_z=%.3f,tru_tot=%.3f) frac_3d=%+.4f\n",
                   evt,kp[evt],
                   pzzpz[evt],pztf[evt],pdt,pvz,vpip[evt],vpip_tot[evt],pf3d,
                   mzzpz[evt],mztf[evt],mdt,mvz,vpim[evt],vpim_tot[evt],mf3d);
            pr++;
        }
    }
    printf("\nN=%d  pip_vz_frac=%+.5f  pim_vz_frac=%+.5f\n",ng,sfp_z/ng,sfm_z/ng);
    printf("Mean 3D raw speed frac (v3d_raw/v_true-1): pip=%+.5f  pim=%+.5f\n",sfp_3d/ng,sfm_3d/ng);
    f->Close();
}

// ============================================================================
// KLong_find_simple_events.C — find events where pion tracks lie in a
//                              coordinate plane (easiest to verify by hand)
//
// PURPOSE:
//   Reconstructs every event in a file and evaluates how closely each pion
//   track aligns with the XZ or YZ plane by checking its azimuthal angle phi:
//
//     phi = atan2(dir.Y(), dir.X())
//     phi ~ 0/180  =>  track in XZ plane (y(z) ~ 0)
//     phi ~ +-90   =>  track in YZ plane (x(z) ~ 0)
//
//   Events where BOTH pions lie in the same coordinate plane are the simplest
//   to verify numerically: the 3D path length reduces to sqrt(dH^2 + dz^2)
//   where H is either X or Y.
//
//   Output is sorted by max(|phi deviation|) of both tracks, smallest first.
//   Events within the tolerance are marked with *.
//   Use the event number with KLong_visualize_event to see the 4-panel display.
//
// INPUT:  one raw simulation .root file (Ntuple1/2/3)
// OUTPUT: <stem>_simple_geometry.txt
//
// USAGE (from VISUALISER/ directory):
//   root -l -q 'KLong_find_simple_events.C("/path/to/file.root")'
//   root -l -q 'KLong_find_simple_events.C("/path/to/file.root", 5.0)'
//   root -l -q 'KLong_find_simple_events.C("/path/to/file.root", 10.0, "out.txt")'
//     phi_tol_deg (default 10): max allowed deviation from XZ or YZ plane
// ============================================================================

#include "vis_common.h"

void KLong_find_simple_events(
    const char* filename    = "Scenario3_Seed1.root",
    double      phi_tol_deg = 10.0,
    const char* output_txt  = "")
{
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) { std::cout << "Cannot open file\n"; return; }

    std::string fname(filename);
    auto extract = [](const std::string& s, const std::string& key) -> double {
        size_t pos = s.find(key + "-");
        if (pos == std::string::npos) return 0.;
        pos += key.size() + 1;
        size_t end = pos;
        while (end < s.size() && std::isdigit((unsigned char)s[end])) end++;
        return (end == pos) ? 0. : std::stod(s.substr(pos, end-pos));
    };
    double z_tof    = extract(fname, "E1");
    double z_pizza1 = extract(fname, "P1");
    double z_pizza2 = extract(fname, "P2");
    (void)z_pizza1; (void)z_pizza2;

    // Select pi+pi-pi0 events
    TTree* tree2 = (TTree*)file->Get("Ntuple2");
    if (!tree2) { file->Close(); return; }
    Double_t evtNb2, dkPDG;
    tree2->SetBranchAddress("evtNb", &evtNb2);
    tree2->SetBranchAddress("DKparticle_PDGEncoding", &dkPDG);
    std::map<int, std::vector<int>> evt_products;
    for (Long64_t i=0; i<tree2->GetEntries(); i++) {
        tree2->GetEntry(i);
        evt_products[(int)evtNb2].push_back((int)dkPDG);
    }
    std::vector<int> sel_events;
    for (auto& kv : evt_products) {
        bool hp=false, hm=false, h0=false;
        for (int c : kv.second) {
            if (c== 211) hp=true; if (c==-211) hm=true; if (c== 111) h0=true;
        }
        if (hp && hm && h0) sel_events.push_back(kv.first);
    }
    std::unordered_set<int> sel_set(sel_events.begin(), sel_events.end());
    std::cout << "Selected " << sel_events.size() << " pi+pi-pi0 events\n";

    // MC truth
    std::unordered_map<int,double> truth_px, truth_py, truth_pz;
    TTree* tree1 = (TTree*)file->Get("Ntuple1");
    if (tree1) {
        Double_t en1, kpx, kpy, kpz, kvx, kvy, kvz;
        tree1->SetBranchAddress("evtNb",       &en1);
        tree1->SetBranchAddress("kaonDK_momX", &kpx);
        tree1->SetBranchAddress("kaonDK_momY", &kpy);
        tree1->SetBranchAddress("kaonDK_momZ", &kpz);
        tree1->SetBranchAddress("kaonDK_posX", &kvx);
        tree1->SetBranchAddress("kaonDK_posY", &kvy);
        tree1->SetBranchAddress("kaonDK_posZ", &kvz);
        for (Long64_t i=0; i<tree1->GetEntries(); i++) {
            tree1->GetEntry(i);
            int en=(int)en1;
            truth_px[en]=kpx; truth_py[en]=kpy; truth_pz[en]=kpz;
        }
    }

    // Load hits
    TTree* tree3 = (TTree*)file->Get("Ntuple3");
    if (!tree3) { file->Close(); return; }
    Double_t evtNb3, hx, hy, hz, ht, pdg3, devID3;
    tree3->SetBranchAddress("evtNb",       &evtNb3);
    tree3->SetBranchAddress("hitX",        &hx);
    tree3->SetBranchAddress("hitY",        &hy);
    tree3->SetBranchAddress("hitZ",        &hz);
    tree3->SetBranchAddress("hitT",        &ht);
    tree3->SetBranchAddress("PDGEncoding", &pdg3);
    tree3->SetBranchAddress("deviceID",    &devID3);

    auto is_tracker = [](int id){ return id>=1    && id<=1952; };
    auto is_fri     = [](int id){ return id>=2001 && id<=2052; };
    auto is_pizza   = [](int id){ return id>=1953 && id<=2000; };
    auto is_tof     = [](int id){ return id>=2053 && id<=2070; };

    TRandom3 rng(42);
    const double smear_t = 0.0015;
    std::unordered_map<int, EventReco> ev_data;
    for (Long64_t i=0; i<tree3->GetEntries(); i++) {
        tree3->GetEntry(i);
        int en=(int)evtNb3;
        if (sel_set.find(en)==sel_set.end()) continue;
        auto& ev=ev_data[en];
        int id=(int)devID3;
        for (int pion=0; pion<2; pion++) {
            if ((int)pdg3 != (pion==0 ? 211 : -211)) continue;
            if (is_tracker(id)||is_fri(id))
                (pion==0 ? ev.hits_pip : ev.hits_pim).push_back({hx,hy,hz,ht,id});
            if (is_pizza(id)) {
                double st=rng.Gaus(ht,smear_t);
                if (pion==0) { if (!ev.has_pip_pizza||st<ev.pip_pizza_time) {
                    ev.has_pip_pizza=true; ev.pip_pizza_time=st;
                    ev.pip_pizza_x=hx; ev.pip_pizza_y=hy; ev.pip_pizza_z=hz; }
                } else       { if (!ev.has_pim_pizza||st<ev.pim_pizza_time) {
                    ev.has_pim_pizza=true; ev.pim_pizza_time=st;
                    ev.pim_pizza_x=hx; ev.pim_pizza_y=hy; ev.pim_pizza_z=hz; } }
            }
            if (is_tof(id)) {
                double st=rng.Gaus(ht,smear_t);
                double smx=rng.Gaus(tof_bar_xc(id),TOF_BAR_HALF_WIDTH);
                double smy=rng.Gaus(tof_bar_yc(id),TOF_Y_UNCERTAINTY);
                if (pion==0) { if (!ev.has_pip_tof||st<ev.pip_tof_time) {
                    ev.has_pip_tof=true; ev.pip_tof_time=st; ev.pip_tof_deviceID=id;
                    ev.pip_tof_x=smx; ev.pip_tof_y=smy; ev.pip_tof_z=z_tof; }
                } else       { if (!ev.has_pim_tof||st<ev.pim_tof_time) {
                    ev.has_pim_tof=true; ev.pim_tof_time=st; ev.pim_tof_deviceID=id;
                    ev.pim_tof_x=smx; ev.pim_tof_y=smy; ev.pim_tof_z=z_tof; } }
            }
        }
    }

    // Minimum angular distance to any of {0, 90, 180, 270} degrees (mod 180)
    auto phi_plane_dev = [](double phi_deg) -> double {
        double p = std::fmod(std::fabs(phi_deg), 180.);
        double d0  = std::min(p, 180.-p);
        double d90 = std::fabs(p - 90.);
        return std::min(d0, d90);
    };
    auto plane_name = [](double phi_deg) -> const char* {
        double p = std::fmod(std::fabs(phi_deg), 180.);
        double d0  = std::min(p, 180.-p);
        double d90 = std::fabs(p - 90.);
        return (d0 <= d90) ? "XZ" : "YZ";
    };

    struct SimpleResult {
        int    event_num;
        double p_true, p_reco, delta_p;
        double phi_pip, phi_pim, theta_pip, theta_pim;
        double max_dev;
        std::string plane;
        double vtx_z;
    };
    std::vector<SimpleResult> results;

    for (int en : sel_events) {
        auto ev_it = ev_data.find(en);
        if (ev_it==ev_data.end()) continue;
        const auto& ev = ev_it->second;
        if (!ev.has_pip_pizza||!ev.has_pip_tof||!ev.has_pim_pizza||!ev.has_pim_tof) continue;

        TrackFitV fp = fit_track_vis(ev.hits_pip, ev.pip_tof_deviceID, z_tof);
        TrackFitV fm = fit_track_vis(ev.hits_pim, ev.pim_tof_deviceID, z_tof);
        if (!fp.valid||!fm.valid) continue;

        TVector3 decay_vtx = poca_vis(fp.orig, fp.dir, fm.orig, fm.dir);
        if (fp.dir.Dot(fm.dir) > 0.9998) continue;

        TVector3 pip_pza = eval_at_z(fp, ev.pip_pizza_z);
        TVector3 pim_pza = eval_at_z(fm, ev.pim_pizza_z);
        TVector3 pip_tof(ev.pip_tof_x, ev.pip_tof_y, ev.pip_tof_z);
        TVector3 pim_tof(ev.pim_tof_x, ev.pim_tof_y, ev.pim_tof_z);
        double pipL =(pip_tof-pip_pza).Mag(), pimL =(pim_tof-pim_pza).Mag();
        double pipdt=ev.pip_tof_time-ev.pip_pizza_time;
        double pimdt=ev.pim_tof_time-ev.pim_pizza_time;
        if (pipdt<=0||pimdt<=0) continue;
        double pipv=(pipL*1e-2)/(pipdt*1e-9), pimv=(pimL*1e-2)/(pimdt*1e-9);
        double pipdp=(pip_pza-decay_vtx).Mag(), pimdp=(pim_pza-decay_vtx).Mag();
        double piptd=ev.pip_pizza_time*1e-9-(pipdp*1e-2)/pipv;
        double pimtd=ev.pim_pizza_time*1e-9-(pimdp*1e-2)/pimv;
        double kdt=0.5*(piptd+pimtd);
        double fl=decay_vtx.Mag();
        double kv=(fl*1e-2)/kdt;
        double bk=kv/2.99792458e8; if(bk>=1.)bk=0.9999;
        double gk=1./std::sqrt(1.-bk*bk);
        double kp=gk*0.497611*bk;
        if (kp<=0||kp>11.) continue;
        if (truth_px.find(en)==truth_px.end()) continue;
        double tp=std::sqrt(truth_px[en]*truth_px[en]+truth_py[en]*truth_py[en]+truth_pz[en]*truth_pz[en]);
        if (tp==0.) continue;

        double phi_p   = std::atan2(fp.dir.Y(), fp.dir.X()) * 180./M_PI;
        double phi_m   = std::atan2(fm.dir.Y(), fm.dir.X()) * 180./M_PI;
        double theta_p = std::acos(std::min(1.0, std::fabs(fp.dir.Z()))) * 180./M_PI;
        double theta_m = std::acos(std::min(1.0, std::fabs(fm.dir.Z()))) * 180./M_PI;
        double devp    = phi_plane_dev(phi_p);
        double devm    = phi_plane_dev(phi_m);
        double maxdev  = std::max(devp, devm);

        std::string plabel = "MIXED";
        if (std::string(plane_name(phi_p)) == plane_name(phi_m))
            plabel = plane_name(phi_p);

        results.push_back({en, tp, kp, kp-tp, phi_p, phi_m, theta_p, theta_m,
                           maxdev, plabel, decay_vtx.Z()});
    }

    std::sort(results.begin(), results.end(),
              [](const SimpleResult& a, const SimpleResult& b){ return a.max_dev < b.max_dev; });

    auto base_name = [](const std::string& p) -> std::string {
        size_t sl=p.find_last_of("/\\");
        std::string fn=(sl==std::string::npos)?p:p.substr(sl+1);
        size_t dot=fn.find_last_of('.'); return (dot==std::string::npos)?fn:fn.substr(0,dot);
    };
    std::string stem = base_name(fname);
    std::string out_path = (std::string(output_txt).empty())
                           ? (stem + "_simple_geometry.txt") : std::string(output_txt);

    std::ofstream ofs(out_path);
    ofs << std::fixed << std::setprecision(3);
    ofs << "# Simple-geometry events — sorted by max phi-plane deviation (smallest = most aligned)\n";
    ofs << "# File: " << filename << "   phi_tol = " << phi_tol_deg << " deg\n";
    ofs << "# Plane: XZ = phi~0/180 (tracks in XZ plane), YZ = phi~+-90 (tracks in YZ plane)\n";
    ofs << "# Pass event number to KLong_visualize_event.C for a diagram or print_event_detail breakdown.\n";
    ofs << "#\n";
    ofs << "# " << std::setw(8)  << "event"
                << std::setw(8)  << "p_true"
                << std::setw(8)  << "p_reco"
                << std::setw(8)  << "delta_p"
                << std::setw(8)  << "phi+"
                << std::setw(8)  << "phi-"
                << std::setw(8)  << "theta+"
                << std::setw(8)  << "theta-"
                << std::setw(9)  << "max_dev"
                << std::setw(7)  << "plane"
                << std::setw(9)  << "vtx_z"
                << "\n";
    ofs << "# " << std::string(95, '-') << "\n";

    int n_xz=0, n_yz=0, n_in_tol=0;
    for (const auto& r : results) {
        bool in_tol = (r.max_dev <= phi_tol_deg);
        if (in_tol) { n_in_tol++; if (r.plane=="XZ") n_xz++; if (r.plane=="YZ") n_yz++; }
        ofs << "  " << std::setw(8)  << r.event_num
                    << std::setw(8)  << r.p_true
                    << std::setw(8)  << r.p_reco
                    << std::setw(8)  << r.delta_p
                    << std::setw(8)  << r.phi_pip
                    << std::setw(8)  << r.phi_pim
                    << std::setw(8)  << r.theta_pip
                    << std::setw(8)  << r.theta_pim
                    << std::setw(9)  << r.max_dev
                    << std::setw(7)  << r.plane
                    << std::setw(9)  << r.vtx_z
                    << (in_tol ? "  *" : "")
                    << "\n";
    }
    ofs.close();

    std::cout << "Reconstructed:  " << results.size() << " events\n";
    std::cout << "Within tol (" << phi_tol_deg << " deg):  " << n_in_tol
              << "  (XZ=" << n_xz << "  YZ=" << n_yz << ")\n";
    std::cout << "Wrote: " << out_path << "\n";
    if (n_in_tol > 0) {
        std::cout << "Top candidate: event " << results[0].event_num
                  << "  max_dev=" << results[0].max_dev << " deg"
                  << "  phi+=" << results[0].phi_pip
                  << "  phi-=" << results[0].phi_pim
                  << "  plane=" << results[0].plane
                  << "  p_true=" << results[0].p_true << " GeV/c\n";
        std::cout << "  -> visualise: root -l -q 'KLong_visualize_event.C(\""
                  << filename << "\", " << results[0].event_num << ")'\n";
    }
    file->Close();
}

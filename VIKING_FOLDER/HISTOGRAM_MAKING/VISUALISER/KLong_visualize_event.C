// ============================================================================
// KLong_visualize_event.C — single-file event display for K_L -> pi+pi-pi0
//
// PURPOSE:
//   For a given raw simulation .root file (one seed):
//     1. Runs the full K_L -> pi+pi-pi0 reconstruction on every event
//     2. Writes a text file sorted by |delta_p| (largest bias first)
//     3. For a specified event number draws a 4-panel figure:
//          Panel 1 (top-left):  3D perspective view
//          Panel 2 (top-right): XZ projection ("side view")
//          Panel 3 (bottom-left): YZ projection ("top view")
//          Panel 4 (bottom-right): XY projection ("front view")
//        and prints a step-by-step kinematic breakdown to stdout.
//
// INPUT:  a raw simulation .root file with Ntuple1, Ntuple2, Ntuple3
//         e.g. T1-240_..._E1-700_1.root   (NOT a _combined_ or _vectors_ file)
// OUTPUT: <stem>_delta_p_sorted.txt     — all events sorted by |delta_p|
//         KLong_event<N>_view.png        — 4-panel figure for event N
//
// USAGE (from VISUALISER/ directory):
//   root -l -q 'KLong_visualize_event.C("/path/to/file.root")'
//   root -l -q 'KLong_visualize_event.C("/path/to/file.root", 42)'
//   root -l -q 'KLong_visualize_event.C("/path/to/file.root", 42, "out.txt")'
// ============================================================================

#include "vis_drawing.h"

void KLong_visualize_event(
    const char* filename     = "Scenario3_Seed1.root",
    int         target_event = -1,
    const char* output_txt   = "")
{
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cout << "Cannot open: " << filename << "\n"; return;
    }

    // === Parse geometry from filename ===
    std::string fname(filename);
    auto extract = [](const std::string& s, const std::string& key) -> double {
        size_t pos = s.find(key + "-");
        if (pos == std::string::npos) return 0.;
        pos += key.size() + 1;
        size_t end = pos;
        while (end < s.size() && std::isdigit((unsigned char)s[end])) end++;
        return (end == pos) ? 0. : std::stod(s.substr(pos, end-pos));
    };
    double z_t1     = extract(fname, "T1");
    double z_t2     = extract(fname, "T2");
    double z_t3     = extract(fname, "T3");
    double z_t4     = extract(fname, "T4");
    double z_pizza1 = extract(fname, "P1");
    double z_pizza2 = extract(fname, "P2");
    double z_fri1   = extract(fname, "F1");
    double z_fri2   = extract(fname, "F2");
    double z_tof    = extract(fname, "E1");

    auto base_name = [](const std::string& p) -> std::string {
        size_t sl = p.find_last_of("/\\");
        std::string fn = (sl == std::string::npos) ? p : p.substr(sl+1);
        size_t dot = fn.find_last_of('.');
        return (dot == std::string::npos) ? fn : fn.substr(0, dot);
    };
    std::string stem = base_name(fname);

    std::string txt_path = (std::string(output_txt).empty())
                           ? (stem + "_delta_p_sorted.txt")
                           : std::string(output_txt);

    std::cout << "Detector geometry parsed:\n"
              << "  T1=" << z_t1 << " T2=" << z_t2
              << " T3=" << z_t3 << " T4=" << z_t4 << " cm\n"
              << "  P1=" << z_pizza1 << " P2=" << z_pizza2 << " cm\n"
              << "  F1=" << z_fri1   << " F2=" << z_fri2   << " cm\n"
              << "  TOF=" << z_tof << " cm\n";

    // === Load event filter from Ntuple2 ===
    TTree* tree2 = (TTree*)file->Get("Ntuple2");
    if (!tree2) { std::cout << "Ntuple2 not found\n"; file->Close(); return; }
    Double_t evtNb2, dkPDG;
    tree2->SetBranchAddress("evtNb",                  &evtNb2);
    tree2->SetBranchAddress("DKparticle_PDGEncoding", &dkPDG);
    std::map<int, std::vector<int>> evt_products;
    for (Long64_t i = 0; i < tree2->GetEntries(); i++) {
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

    // === Load MC truth from Ntuple1 ===
    std::unordered_map<int,double> truth_px, truth_py, truth_pz;
    std::unordered_map<int, TVector3> truth_vtx;
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
        for (Long64_t i = 0; i < tree1->GetEntries(); i++) {
            tree1->GetEntry(i);
            int en = (int)en1;
            truth_px[en] = kpx; truth_py[en] = kpy; truth_pz[en] = kpz;
            truth_vtx[en] = TVector3(kvx, kvy, kvz);
        }
    }

    // === Load hits from Ntuple3 ===
    TTree* tree3 = (TTree*)file->Get("Ntuple3");
    if (!tree3) { std::cout << "Ntuple3 not found\n"; file->Close(); return; }
    Double_t evtNb3, hx, hy, hz, ht, pdg3, devID3;
    tree3->SetBranchAddress("evtNb",       &evtNb3);
    tree3->SetBranchAddress("hitX",        &hx);
    tree3->SetBranchAddress("hitY",        &hy);
    tree3->SetBranchAddress("hitZ",        &hz);
    tree3->SetBranchAddress("hitT",        &ht);
    tree3->SetBranchAddress("PDGEncoding", &pdg3);
    tree3->SetBranchAddress("deviceID",    &devID3);

    auto is_tracker = [](int id){ return id >= 1    && id <= 1952; };
    auto is_fri     = [](int id){ return id >= 2001 && id <= 2052; };
    auto is_pizza   = [](int id){ return id >= 1953 && id <= 2000; };
    auto is_tof     = [](int id){ return id >= 2053 && id <= 2070; };

    TRandom3 rng(42);
    const double smear_t_sigma = 0.0015; // ns

    std::unordered_map<int, EventReco> ev_data;
    for (Long64_t i = 0; i < tree3->GetEntries(); i++) {
        tree3->GetEntry(i);
        int en = (int)evtNb3;
        if (sel_set.find(en) == sel_set.end()) continue;
        auto& ev = ev_data[en];
        int id = (int)devID3;
        for (int pion = 0; pion < 2; pion++) {
            int pdg_want = (pion == 0) ? 211 : -211;
            if ((int)pdg3 != pdg_want) continue;
            if (is_tracker(id) || is_fri(id)) {
                auto& vec = (pion==0) ? ev.hits_pip : ev.hits_pim;
                vec.push_back({hx, hy, hz, ht, id});
            }
            if (is_pizza(id)) {
                double st = rng.Gaus(ht, smear_t_sigma);
                if (pion == 0) {
                    if (!ev.has_pip_pizza || st < ev.pip_pizza_time) {
                        ev.has_pip_pizza=true; ev.pip_pizza_time=st;
                        ev.pip_pizza_x=hx; ev.pip_pizza_y=hy; ev.pip_pizza_z=hz;
                    }
                } else {
                    if (!ev.has_pim_pizza || st < ev.pim_pizza_time) {
                        ev.has_pim_pizza=true; ev.pim_pizza_time=st;
                        ev.pim_pizza_x=hx; ev.pim_pizza_y=hy; ev.pim_pizza_z=hz;
                    }
                }
            }
            if (is_tof(id)) {
                double st  = rng.Gaus(ht, smear_t_sigma);
                double smx = rng.Gaus(tof_bar_xc(id), TOF_BAR_HALF_WIDTH);
                double smy = rng.Gaus(tof_bar_yc(id), TOF_Y_UNCERTAINTY);
                if (pion == 0) {
                    if (!ev.has_pip_tof || st < ev.pip_tof_time) {
                        ev.has_pip_tof=true; ev.pip_tof_time=st;
                        ev.pip_tof_deviceID=id;
                        ev.pip_tof_x=smx; ev.pip_tof_y=smy; ev.pip_tof_z=z_tof;
                    }
                } else {
                    if (!ev.has_pim_tof || st < ev.pim_tof_time) {
                        ev.has_pim_tof=true; ev.pim_tof_time=st;
                        ev.pim_tof_deviceID=id;
                        ev.pim_tof_x=smx; ev.pim_tof_y=smy; ev.pim_tof_z=z_tof;
                    }
                }
            }
        }
    }

    // === Reconstruction loop ===
    std::vector<VisEvent> reco_events;

    for (int en : sel_events) {
        auto ev_it = ev_data.find(en);
        if (ev_it == ev_data.end()) continue;
        const auto& ev = ev_it->second;
        if (!ev.has_pip_pizza || !ev.has_pip_tof ||
            !ev.has_pim_pizza || !ev.has_pim_tof) continue;

        TrackFitV fp = fit_track_vis(ev.hits_pip, ev.pip_tof_deviceID, z_tof);
        TrackFitV fm = fit_track_vis(ev.hits_pim, ev.pim_tof_deviceID, z_tof);
        if (!fp.valid || !fm.valid) continue;

        TVector3 decay_vtx = poca_vis(fp.orig, fp.dir, fm.orig, fm.dir);
        if (fp.dir.Dot(fm.dir) > 0.9998) continue;

        TVector3 pip_pizza_pos = eval_at_z(fp, ev.pip_pizza_z);
        TVector3 pim_pizza_pos = eval_at_z(fm, ev.pim_pizza_z);
        TVector3 pip_tof_pos(ev.pip_tof_x, ev.pip_tof_y, ev.pip_tof_z);
        TVector3 pim_tof_pos(ev.pim_tof_x, ev.pim_tof_y, ev.pim_tof_z);

        double pip_L  = (pip_tof_pos - pip_pizza_pos).Mag();
        double pim_L  = (pim_tof_pos - pim_pizza_pos).Mag();
        double pip_dt = ev.pip_tof_time - ev.pip_pizza_time;
        double pim_dt = ev.pim_tof_time - ev.pim_pizza_time;
        if (pip_dt <= 0 || pim_dt <= 0) continue;

        double pip_v = (pip_L * 1e-2) / (pip_dt * 1e-9);
        double pim_v = (pim_L * 1e-2) / (pim_dt * 1e-9);

        double pip_path = (pip_pizza_pos - decay_vtx).Mag();
        double pim_path = (pim_pizza_pos - decay_vtx).Mag();
        double pip_td = ev.pip_pizza_time * 1e-9 - (pip_path * 1e-2) / pip_v;
        double pim_td = ev.pim_pizza_time * 1e-9 - (pim_path * 1e-2) / pim_v;
        double kaon_decay_t = 0.5 * (pip_td + pim_td);

        double flight = decay_vtx.Mag();
        double kaon_v = (flight * 1e-2) / kaon_decay_t;
        double beta_K = kaon_v / 2.99792458e8;
        if (beta_K >= 1.0) beta_K = 0.9999;
        double gamma_K = 1. / std::sqrt(1. - beta_K * beta_K);
        double kaon_p  = gamma_K * 0.497611 * beta_K;
        if (kaon_p <= 0 || kaon_p > 11.) continue;

        auto tv_it = truth_vtx.find(en);
        if (tv_it == truth_vtx.end()) continue;
        double tp = std::sqrt(truth_px[en]*truth_px[en] +
                              truth_py[en]*truth_py[en] +
                              truth_pz[en]*truth_pz[en]);
        if (tp == 0.) continue;

        VisEvent ve;
        ve.event_num    = en;
        ve.reco_p       = kaon_p;
        ve.true_p       = tp;
        ve.delta_p      = kaon_p - tp;
        ve.reco_vtx     = decay_vtx;
        ve.true_vtx     = tv_it->second;
        ve.hits_pip     = ev.hits_pip;
        ve.hits_pim     = ev.hits_pim;
        ve.pip_pizza_raw = TVector3(ev.pip_pizza_x, ev.pip_pizza_y, ev.pip_pizza_z);
        ve.pim_pizza_raw = TVector3(ev.pim_pizza_x, ev.pim_pizza_y, ev.pim_pizza_z);
        ve.pip_tof_pos  = pip_tof_pos;
        ve.pim_tof_pos  = pim_tof_pos;
        ve.pip_pizza_pos = pip_pizza_pos;
        ve.pim_pizza_pos = pim_pizza_pos;
        ve.fit_pip      = fp;
        ve.fit_pim      = fm;
        ve.true_kaon_dir = TVector3(truth_px[en]/tp, truth_py[en]/tp, truth_pz[en]/tp);
        ve.pip_pizza_time    = ev.pip_pizza_time;
        ve.pim_pizza_time    = ev.pim_pizza_time;
        ve.pip_tof_time      = ev.pip_tof_time;
        ve.pim_tof_time      = ev.pim_tof_time;
        ve.pip_tof_devID     = ev.pip_tof_deviceID;
        ve.pim_tof_devID     = ev.pim_tof_deviceID;
        ve.pip_track_cm      = pip_L;
        ve.pim_track_cm      = pim_L;
        ve.pip_dt_ns         = pip_dt;
        ve.pim_dt_ns         = pim_dt;
        ve.pip_v             = pip_v;
        ve.pim_v             = pim_v;
        ve.pip_decay_path_cm = pip_path;
        ve.pim_decay_path_cm = pim_path;
        ve.pip_decay_time    = pip_td;
        ve.pim_decay_time    = pim_td;
        ve.kaon_decay_time   = kaon_decay_t;
        ve.kaon_flight_cm    = flight;
        ve.kaon_v            = kaon_v;
        ve.kaon_beta         = beta_K;
        ve.kaon_gamma        = gamma_K;
        ve.source_file       = "";
        reco_events.push_back(ve);
    }

    std::cout << "Reconstructed " << reco_events.size() << " events\n";

    // === Sort by |delta_p| descending and write text file ===
    std::vector<VisEvent*> sorted_refs;
    for (auto& ve : reco_events) sorted_refs.push_back(&ve);
    std::sort(sorted_refs.begin(), sorted_refs.end(),
              [](const VisEvent* a, const VisEvent* b){
                  return std::fabs(a->delta_p) > std::fabs(b->delta_p);
              });

    std::ofstream ofs(txt_path);
    ofs << std::fixed << std::setprecision(4);
    ofs << "# Events sorted by |delta_p| = |p_reco - p_true| (largest first)\n";
    ofs << "# File: " << filename << "\n";
    ofs << "# "
        << std::setw(8)  << "event"    << "  "
        << std::setw(10) << "p_true"   << "  "
        << std::setw(10) << "p_reco"   << "  "
        << std::setw(10) << "delta_p"  << "  "
        << std::setw(10) << "|delta_p|"<< "  "
        << std::setw(10) << "reco_vx"  << "  "
        << std::setw(10) << "reco_vy"  << "  "
        << std::setw(10) << "reco_vz"  << "  "
        << std::setw(10) << "true_vz"  << "\n";
    ofs << std::string(110, '-') << "\n";
    for (const VisEvent* ve : sorted_refs) {
        ofs << "  "
            << std::setw(8)  << ve->event_num          << "  "
            << std::setw(10) << ve->true_p              << "  "
            << std::setw(10) << ve->reco_p              << "  "
            << std::setw(10) << ve->delta_p             << "  "
            << std::setw(10) << std::fabs(ve->delta_p)  << "  "
            << std::setw(10) << ve->reco_vtx.X()        << "  "
            << std::setw(10) << ve->reco_vtx.Y()        << "  "
            << std::setw(10) << ve->reco_vtx.Z()        << "  "
            << std::setw(10) << ve->true_vtx.Z()        << "\n";
    }
    for (auto& ve : reco_events) ve.source_file = std::string(filename);
    ofs.close();
    std::cout << "Wrote sorted event list to: " << txt_path << "\n";

    // === Visualization ===
    if (target_event < 0) {
        std::cout << "No target event specified. Pass an event number as the "
                     "second argument to generate the 3D view.\n"
                     "Hint: check the top entries in " << txt_path << "\n";
        file->Close();
        return;
    }

    const VisEvent* draw_ev = nullptr;
    for (const auto& ve : reco_events) {
        if (ve.event_num == target_event) { draw_ev = &ve; break; }
    }
    if (!draw_ev) {
        std::cout << "Event " << target_event
                  << " not found in reconstructed events (may fail PIZZA/TOF/track cuts).\n";
        file->Close();
        return;
    }

    print_event_detail(*draw_ev);

    // === 4-panel canvas ===
    TCanvas* c = new TCanvas(Form("c_ev%d", target_event),
                             Form("Event %d  |  file: %s", target_event, stem.c_str()),
                             1400, 1000);
    c->SetFillColor(10);
    TPad* pad3d = new TPad("pad3d", "3D view", 0.00, 0.50, 0.60, 1.00);
    TPad* padXZ = new TPad("padXZ", "XZ",      0.60, 0.50, 1.00, 1.00);
    TPad* padYZ = new TPad("padYZ", "YZ",      0.00, 0.00, 0.50, 0.50);
    TPad* padXY = new TPad("padXY", "XY",      0.50, 0.00, 1.00, 0.50);
    for (TPad* p : {pad3d, padXZ, padYZ, padXY}) {
        c->cd(); p->Draw();
        p->SetLeftMargin(0.10); p->SetRightMargin(0.02);
        p->SetTopMargin(0.09);  p->SetBottomMargin(0.10);
    }
    draw_3d_panel (pad3d, *draw_ev, z_t1, z_t2, z_t3, z_t4, z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);
    draw_proj_panel(padXZ, *draw_ev, "XZ", z_t1, z_t2, z_t3, z_t4, z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);
    draw_proj_panel(padYZ, *draw_ev, "YZ", z_t1, z_t2, z_t3, z_t4, z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);
    draw_proj_panel(padXY, *draw_ev, "XY", z_t1, z_t2, z_t3, z_t4, z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);

    c->cd();
    TLatex tit;
    tit.SetNDC(true); tit.SetTextSize(0.022); tit.SetTextFont(42);
    tit.DrawLatex(0.01, 0.003,
                  Form("File: %s    Event: %d    p_{true}=%.4f GeV/c    "
                       "p_{reco}=%.4f GeV/c    #Deltap=%.4f GeV/c",
                       stem.c_str(), target_event,
                       draw_ev->true_p, draw_ev->reco_p, draw_ev->delta_p));

    std::string png_name = Form("KLong_event%d_view.png", target_event);
    c->SaveAs(png_name.c_str());
    std::cout << "Saved: " << png_name << "\n";

    file->Close();
}

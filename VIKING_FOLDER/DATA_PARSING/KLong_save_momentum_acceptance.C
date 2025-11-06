#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <vector>
#include <iostream>
#include <cmath>

struct HitInfo { double x, y, z, t; int deviceID; };

TVector3 closest_point_between_lines(const TVector3& p1, const TVector3& v1,
                                    const TVector3& p2, const TVector3& v2) {
    TVector3 w0 = p1 - p2;
    double a = v1.Dot(v1);
    double b = v1.Dot(v2);
    double c = v2.Dot(v2);
    double d = v1.Dot(w0);
    double e = v2.Dot(w0);
    double denom = a*c - b*b;
    double sc = (b*e - c*d) / denom;
    double tc = (a*e - b*d) / denom;
    TVector3 point1 = p1 + v1 * sc;
    TVector3 point2 = p2 + v2 * tc;
    return 0.5 * (point1 + point2); // midpoint is PoCA
}

void KLong_save_momentum_acceptance(const char* filename = "Scenario3_Seed1.root") {
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) { std::cout << "Cannot open root file\n"; return; }

    // --- Step 1: Find events with triple pion decay (π+, π-, π0) ---
    TTree *tree2 = (TTree*)file->Get("Ntuple2");
    if (!tree2) { std::cout << "Tree Ntuple2 not found!\n"; file->Close(); return; }

    Double_t evtNb2, DKparticle_PDGEncoding;
    tree2->SetBranchAddress("evtNb", &evtNb2);
    tree2->SetBranchAddress("DKparticle_PDGEncoding", &DKparticle_PDGEncoding);

    std::map<int, std::vector<int>> event_products;
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        event_products[(int)evtNb2].push_back((int)DKparticle_PDGEncoding);
    }

    std::vector<int> selected_events;
    for (const auto& kv : event_products) {
        std::vector<int> products = kv.second;
        bool has_pim = false, has_pi0 = false, has_pip = false;
        for (int pdg : products) {
            if (pdg == -211) has_pim = true;
            if (pdg == 111) has_pi0 = true;
            if (pdg == 211) has_pip = true;
        }
        if (has_pim && has_pi0 && has_pip)
            selected_events.push_back(kv.first);
    }

    std::cout << "Found " << selected_events.size() << " events with π+, π-, π0 as direct kaon decay products.\n";

    TTree *tree3 = (TTree*)file->Get("Ntuple3");
    if (!tree3) { std::cout << "Tree Ntuple3 not found!\n"; file->Close(); return; }

    Double_t evtNb, x, y, z, t, pdg, deviceID;
    tree3->SetBranchAddress("evtNb", &evtNb);
    tree3->SetBranchAddress("hitX", &x);
    tree3->SetBranchAddress("hitY", &y);
    tree3->SetBranchAddress("hitZ", &z);
    tree3->SetBranchAddress("hitT", &t);
    tree3->SetBranchAddress("PDGEncoding", &pdg);
    tree3->SetBranchAddress("deviceID", &deviceID);

    auto is_pizza = [](int id) { return (id >= 489 && id <= 536); };
    auto is_tof   = [](int id) { return (id >= 589 && id <= 606); };
    auto is_tracker = [](int id) { return (id >= 1 && id <= 488); };
    auto is_fri = [](int id) { return (id >= 537 && id <= 588); };

    TRandom3 randGen(0);
    double smear_sigma = 5;      // cm
    double smear_time_sigma = 0.0015; // ns

    std::vector<double> reco_p, true_p;
    std::vector<TVector3> reco_vertex, true_vertex;

    std::vector<double> all_true_p;      // True momentum for ALL selected events
    std::vector<int> all_reco_flags;     // Reconstruction flag for ALL selected events

    int event_counter = 0;
    for (int event_number : selected_events) {
        event_counter++; // 1-based index for selected events

        // --- True kaon decay momentum and vertex from MC truth (Ntuple1) ---
        TTree *tree1 = (TTree*)file->Get("Ntuple1");
        double true_px = 0, true_py = 0, true_pz = 0;
        double true_vx = 0, true_vy = 0, true_vz = 0;
        bool found_truth = false;
        if (tree1) {
            Double_t evtNb1, kaonDK_momX, kaonDK_momY, kaonDK_momZ;
            Double_t kaonDK_posX, kaonDK_posY, kaonDK_posZ;
            tree1->SetBranchAddress("evtNb", &evtNb1);
            tree1->SetBranchAddress("kaonDK_momX", &kaonDK_momX);
            tree1->SetBranchAddress("kaonDK_momY", &kaonDK_momY);
            tree1->SetBranchAddress("kaonDK_momZ", &kaonDK_momZ);
            tree1->SetBranchAddress("kaonDK_posX", &kaonDK_posX);
            tree1->SetBranchAddress("kaonDK_posY", &kaonDK_posY);
            tree1->SetBranchAddress("kaonDK_posZ", &kaonDK_posZ);

            Long64_t nEntries1 = tree1->GetEntries();
            for (Long64_t i = 0; i < nEntries1; ++i) {
                tree1->GetEntry(i);
                if ((int)evtNb1 == event_number) {
                    true_px = kaonDK_momX;
                    true_py = kaonDK_momY;
                    true_pz = kaonDK_momZ;
                    true_vx = kaonDK_posX;
                    true_vy = kaonDK_posY;
                    true_vz = kaonDK_posZ;
                    found_truth = true;
                    break;
                }
            }
        }
        double true_p_mag = std::sqrt(true_px*true_px + true_py*true_py + true_pz*true_pz);
        TVector3 true_vertex_vec(true_vx, true_vy, true_vz);

        // --- Try to reconstruct the kaon momentum ---
        double pip_pizza_time = -1, pip_tof_time = -1;
        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;
        double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;

        double pim_pizza_time = -1, pim_tof_time = -1;
        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;
        double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;

        std::vector<HitInfo> hits_pip, hits_pim;

        Long64_t nEntries3 = tree3->GetEntries();
        for (Long64_t i = 0; i < nEntries3; ++i) {
            tree3->GetEntry(i);
            if ((int)evtNb != event_number) continue;
            if ((int)pdg == 211) {
                if (is_tracker((int)deviceID) || is_fri((int)deviceID))
                    hits_pip.push_back({x, y, z, t, (int)deviceID});
                if (is_pizza((int)deviceID)) {
                    double smeared_t = randGen.Gaus(t, smear_time_sigma);
                    if (pip_pizza_time < 0 || smeared_t < pip_pizza_time) {
                        pip_pizza_time = smeared_t;
                        pip_pizza_x = x; pip_pizza_y = y; pip_pizza_z = z;
                    }
                }
                if (is_tof((int)deviceID)) {
                    double smeared_t = randGen.Gaus(t, smear_time_sigma);
                    if (pip_tof_time < 0 || smeared_t < pip_tof_time) {
                        pip_tof_time = smeared_t;
                        pip_tof_x = x; pip_tof_y = y; pip_tof_z = z;
                    }
                }
            }
            if ((int)pdg == -211) {
                if (is_tracker((int)deviceID) || is_fri((int)deviceID))
                    hits_pim.push_back({x, y, z, t, (int)deviceID});
                if (is_pizza((int)deviceID)) {
                    double smeared_t = randGen.Gaus(t, smear_time_sigma);
                    if (pim_pizza_time < 0 || smeared_t < pim_pizza_time) {
                        pim_pizza_time = smeared_t;
                        pim_pizza_x = x; pim_pizza_y = y; pim_pizza_z = z;
                    }
                }
                if (is_tof((int)deviceID)) {
                    double smeared_t = randGen.Gaus(t, smear_time_sigma);
                    if (pim_tof_time < 0 || smeared_t < pim_tof_time) {
                        pim_tof_time = smeared_t;
                        pim_tof_x = x; pim_tof_y = y; pim_tof_z = z;
                    }
                }
            }
        }

        bool reconstructable = false;
        double kaon_p = -1;
        TVector3 decay_vertex;

        if (hits_pip.size() >= 2 && hits_pim.size() >= 2 &&
            pip_pizza_time >= 0 && pip_tof_time >= 0 && pim_pizza_time >= 0 && pim_tof_time >= 0) {

            TVector3 pip_pizza_pos(
                randGen.Gaus(pip_pizza_x, smear_sigma),
                randGen.Gaus(pip_pizza_y, smear_sigma),
                pip_pizza_z
            );
            TVector3 pim_pizza_pos(
                randGen.Gaus(pim_pizza_x, smear_sigma),
                randGen.Gaus(pim_pizza_y, smear_sigma),
                pim_pizza_z
            );

            TVector3 pip_start(hits_pip.front().x, hits_pip.front().y, hits_pip.front().z);
            TVector3 pip_end(hits_pip.back().x, hits_pip.back().y, hits_pip.back().z);
            TVector3 pip_dir = (pip_end - pip_start).Unit();

            TVector3 pim_start(hits_pim.front().x, hits_pim.front().y, hits_pim.front().z);
            TVector3 pim_end(hits_pim.back().x, hits_pim.back().y, hits_pim.back().z);
            TVector3 pim_dir = (pim_end - pim_start).Unit();

            decay_vertex = closest_point_between_lines(pip_start, pip_dir, pim_start, pim_dir);

            TVector3 pip_tof_pos(
                randGen.Gaus(pip_tof_x, smear_sigma),
                randGen.Gaus(pip_tof_y, smear_sigma),
                pip_tof_z
            );
            TVector3 pim_tof_pos(
                randGen.Gaus(pim_tof_x, smear_sigma),
                randGen.Gaus(pim_tof_y, smear_sigma),
                pim_tof_z
            );

            double pip_track_length_cm = (pip_tof_pos - pip_pizza_pos).Mag(); // cm
            double pim_track_length_cm = (pim_tof_pos - pim_pizza_pos).Mag(); // cm

            double pip_dt_ns = pip_tof_time - pip_pizza_time;
            double pim_dt_ns = pim_tof_time - pim_pizza_time;

            double pip_v = (pip_dt_ns > 0) ? (pip_track_length_cm * 1e-2) / (pip_dt_ns * 1e-9) : 0; // m/s
            double pim_v = (pim_dt_ns > 0) ? (pim_track_length_cm * 1e-2) / (pim_dt_ns * 1e-9) : 0; // m/s

            double pip_path_cm = (pip_pizza_pos - decay_vertex).Mag();
            double pip_dt_s = pip_path_cm * 1e-2 / pip_v; // s
            double pip_decay_time = pip_pizza_time * 1e-9 - pip_dt_s; // s

            double pim_path_cm = (pim_pizza_pos - decay_vertex).Mag();
            double pim_dt_s = pim_path_cm * 1e-2 / pim_v; // s
            double pim_decay_time = pim_pizza_time * 1e-9 - pim_dt_s; // s

            double kaon_decay_time = 0.5 * (pip_decay_time + pim_decay_time); // s

            TVector3 kaon_prod(0,0,0);
            double kaon_flight_length_cm = (decay_vertex - kaon_prod).Mag();
            double kaon_flight_time_s = kaon_decay_time; // assuming kaon produced at t=0

            double kaon_velocity = (kaon_flight_length_cm * 1e-2) / kaon_flight_time_s; // m/s

            double m_K = 0.497611; // GeV/c^2
            double beta_K = kaon_velocity / 2.99792458e8;
            if (beta_K >= 1.0) beta_K = 0.9999;
            double gamma_K = 1.0 / std::sqrt(1 - beta_K*beta_K);
            kaon_p = gamma_K * m_K * beta_K; // GeV/c

            if (kaon_p > 0 && kaon_p <= 11) {
                reconstructable = true;
            }
        }

        // Print event info
        // std::cout << "Selected Event #" << event_counter
                  << " | True p: " << true_p_mag;
        if (reconstructable) {
            // std::cout << " | Reco p: " << kaon_p << std::endl;
            reco_p.push_back(kaon_p);
            true_p.push_back(true_p_mag);
            reco_vertex.push_back(decay_vertex);
            true_vertex.push_back(true_vertex_vec);
        } else {
            // std::cout << " | Reco p: not reconstructable" << std::endl;
        }

        // Always add true momentum
        all_true_p.push_back(true_p_mag);
        // Add reconstruction flag
        all_reco_flags.push_back(reconstructable ? 1 : 0);
    }

    // --- Save event info to file ---
    #include <string>
    std::string inFileName(filename);
    size_t lastdot = inFileName.find_last_of(".");
    std::string base = (lastdot == std::string::npos) ? inFileName : inFileName.substr(0, lastdot);
    std::string outFileName = base + "_acceptance.root";
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    TTree *outTree = new TTree("kaonEventInfo", "Kaon event info for acceptance study");

    std::vector<double> true_mom_vec = all_true_p;
    std::vector<int> reco_flag_vec = all_reco_flags;
    int n_triple_pion_events = selected_events.size();

    outTree->Branch("n_triple_pion_events", &n_triple_pion_events);
    outTree->Branch("true_mom_vec", &true_mom_vec);
    outTree->Branch("reco_flag_vec", &reco_flag_vec);

    outTree->Fill();
    outFile->Write();
    outFile->Close();
    file->Close();

    std::cout << "Saved kaon event info to acceptance_seed_.root\n";
}
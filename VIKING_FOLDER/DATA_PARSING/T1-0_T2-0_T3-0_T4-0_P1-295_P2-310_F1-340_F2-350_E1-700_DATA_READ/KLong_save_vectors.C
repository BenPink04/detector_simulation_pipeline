#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "ROOT/RDataFrame.hxx"
#include <chrono>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>

struct HitInfo { double x, y, z, t; int deviceID; };

struct TruthInfo {
    double px, py, pz;
    double vx, vy, vz;
};

struct EventReco {
    std::vector<HitInfo> hits_pip;
    std::vector<HitInfo> hits_pim;
    bool has_pip_pizza = false;
    double pip_pizza_time = 0;
    double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;
    bool has_pip_tof = false;
    double pip_tof_time = 0;
    double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;
    bool has_pim_pizza = false;
    double pim_pizza_time = 0;
    double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;
    bool has_pim_tof = false;
    double pim_tof_time = 0;
    double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;
};

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

void KLong_save_vectors(const char* filename = "Scenario3_Seed1.root") {
    ROOT::EnableImplicitMT(4);
    auto start_time = std::chrono::steady_clock::now();

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

    size_t total_events = selected_events.size();
    std::cout << "Found " << total_events << " events with π+, π-, π0 as direct kaon decay products.\n";
    std::cout << "Processing " << total_events << " selected events...\n";

    std::unordered_set<int> selected_event_set(selected_events.begin(), selected_events.end());


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
    std::vector<double> reco_vertex_x, reco_vertex_y, reco_vertex_z;
    std::vector<double> true_vertex_x, true_vertex_y, true_vertex_z;

    // --- Build truth lookup from Ntuple1 (once) ---
    std::unordered_map<int, TruthInfo> truth_map;
    TTree *tree1 = (TTree*)file->Get("Ntuple1");
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
            truth_map[(int)evtNb1] = {
                kaonDK_momX, kaonDK_momY, kaonDK_momZ,
                kaonDK_posX, kaonDK_posY, kaonDK_posZ
            };
        }
    }

    // --- Build hit lookup from Ntuple3 (once) ---
    std::unordered_map<int, EventReco> event_data;
    Long64_t nEntries3 = tree3->GetEntries();
    for (Long64_t i = 0; i < nEntries3; ++i) {
        tree3->GetEntry(i);
        int evt_id = (int)evtNb;
        if (selected_event_set.find(evt_id) == selected_event_set.end()) continue;

        auto &ev = event_data[evt_id];
        if ((int)pdg == 211) {
            if (is_tracker((int)deviceID) || is_fri((int)deviceID)) {
                ev.hits_pip.push_back({x, y, z, t, (int)deviceID});
            }
            if (is_pizza((int)deviceID)) {
                double smeared_t = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pip_pizza || smeared_t < ev.pip_pizza_time) {
                    ev.has_pip_pizza = true;
                    ev.pip_pizza_time = smeared_t;
                    ev.pip_pizza_x = x; ev.pip_pizza_y = y; ev.pip_pizza_z = z;
                }
            }
            if (is_tof((int)deviceID)) {
                double smeared_t = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pip_tof || smeared_t < ev.pip_tof_time) {
                    ev.has_pip_tof = true;
                    ev.pip_tof_time = smeared_t;
                    ev.pip_tof_x = x; ev.pip_tof_y = y; ev.pip_tof_z = z;
                }
            }
        }
        if ((int)pdg == -211) {
            if (is_tracker((int)deviceID) || is_fri((int)deviceID)) {
                ev.hits_pim.push_back({x, y, z, t, (int)deviceID});
            }
            if (is_pizza((int)deviceID)) {
                double smeared_t = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pim_pizza || smeared_t < ev.pim_pizza_time) {
                    ev.has_pim_pizza = true;
                    ev.pim_pizza_time = smeared_t;
                    ev.pim_pizza_x = x; ev.pim_pizza_y = y; ev.pim_pizza_z = z;
                }
            }
            if (is_tof((int)deviceID)) {
                double smeared_t = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pim_tof || smeared_t < ev.pim_tof_time) {
                    ev.has_pim_tof = true;
                    ev.pim_tof_time = smeared_t;
                    ev.pim_tof_x = x; ev.pim_tof_y = y; ev.pim_tof_z = z;
                }
            }
        }
    }

    int event_counter = 0;
    for (int event_number : selected_events) {
        event_counter++;
        if (event_counter == 1 || event_counter % 1000 == 0 || event_counter == (int)total_events) {
            auto now = std::chrono::steady_clock::now();
            double elapsed_s = std::chrono::duration<double>(now - start_time).count();
            double pct = (total_events > 0) ? (100.0 * event_counter / total_events) : 0.0;
            std::cout << "Progress " << event_counter << "/" << total_events
                      << " (" << pct << "%) - elapsed " << elapsed_s << " s" << std::endl;
        }
        double pip_pizza_time = -1, pip_tof_time = -1;
        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;
        double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;

        double pim_pizza_time = -1, pim_tof_time = -1;
        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;
        double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;

        std::vector<HitInfo> hits_pip, hits_pim;
        auto ev_it = event_data.find(event_number);
        if (ev_it != event_data.end()) {
            hits_pip = ev_it->second.hits_pip;
            hits_pim = ev_it->second.hits_pim;
            if (ev_it->second.has_pip_pizza) {
                pip_pizza_time = ev_it->second.pip_pizza_time;
                pip_pizza_x = ev_it->second.pip_pizza_x;
                pip_pizza_y = ev_it->second.pip_pizza_y;
                pip_pizza_z = ev_it->second.pip_pizza_z;
            }
            if (ev_it->second.has_pip_tof) {
                pip_tof_time = ev_it->second.pip_tof_time;
                pip_tof_x = ev_it->second.pip_tof_x;
                pip_tof_y = ev_it->second.pip_tof_y;
                pip_tof_z = ev_it->second.pip_tof_z;
            }
            if (ev_it->second.has_pim_pizza) {
                pim_pizza_time = ev_it->second.pim_pizza_time;
                pim_pizza_x = ev_it->second.pim_pizza_x;
                pim_pizza_y = ev_it->second.pim_pizza_y;
                pim_pizza_z = ev_it->second.pim_pizza_z;
            }
            if (ev_it->second.has_pim_tof) {
                pim_tof_time = ev_it->second.pim_tof_time;
                pim_tof_x = ev_it->second.pim_tof_x;
                pim_tof_y = ev_it->second.pim_tof_y;
                pim_tof_z = ev_it->second.pim_tof_z;
            }
        }

        if (hits_pip.size() < 2 || hits_pim.size() < 2) continue;
        if (pip_pizza_time < 0 || pip_tof_time < 0 || pim_pizza_time < 0 || pim_tof_time < 0)
            continue;

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

        TVector3 decay_vertex = closest_point_between_lines(pip_start, pip_dir, pim_start, pim_dir);

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
        double kaon_p = gamma_K * m_K * beta_K; // GeV/c

        // --- True kaon decay momentum and vertex from MC truth (lookup) ---
        double true_px = 0, true_py = 0, true_pz = 0;
        double true_vx = 0, true_vy = 0, true_vz = 0;
        bool found_truth = false;
        auto truth_it = truth_map.find(event_number);
        if (truth_it != truth_map.end()) {
            found_truth = true;
            true_px = truth_it->second.px;
            true_py = truth_it->second.py;
            true_pz = truth_it->second.pz;
            true_vx = truth_it->second.vx;
            true_vy = truth_it->second.vy;
            true_vz = truth_it->second.vz;
        }

        double true_p_mag = std::sqrt(true_px*true_px + true_py*true_py + true_pz*true_pz);
        TVector3 true_vertex_vec(true_vx, true_vy, true_vz);

        if (kaon_p > 11) continue;
        if (found_truth && true_p_mag != 0) {
            // Print as soon as the event is processed
            std::cout << "Event " << event_number
                      << " | Reco p: " << kaon_p
                      << " | True p: " << true_p_mag
                      << " | Reco vertex: (" << decay_vertex.X() << ", " << decay_vertex.Y() << ", " << decay_vertex.Z() << ")"
                      << " | True vertex: (" << true_vertex_vec.X() << ", " << true_vertex_vec.Y() << ", " << true_vertex_vec.Z() << ")"
                      << std::endl;

            reco_p.push_back(kaon_p);
            true_p.push_back(true_p_mag);
            reco_vertex_x.push_back(decay_vertex.X());
            reco_vertex_y.push_back(decay_vertex.Y());
            reco_vertex_z.push_back(decay_vertex.Z());
            true_vertex_x.push_back(true_vertex_vec.X());
            true_vertex_y.push_back(true_vertex_vec.Y());
            true_vertex_z.push_back(true_vertex_vec.Z());
        }
    }

    // --- Save vectors to file ---
    std::string inFileName(filename);
    size_t lastdot = inFileName.find_last_of(".");
    std::string base = (lastdot == std::string::npos) ? inFileName : inFileName.substr(0, lastdot);
    std::string outFileName = base + "_vectors.root";
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    TTree *outTree = new TTree("kaonVectors", "Kaon momentum and vertex vectors");

    std::vector<double> *p_reco = &reco_p;
    std::vector<double> *p_true = &true_p;
    std::vector<double> *v_reco_x = &reco_vertex_x;
    std::vector<double> *v_reco_y = &reco_vertex_y;
    std::vector<double> *v_reco_z = &reco_vertex_z;
    std::vector<double> *v_true_x = &true_vertex_x;
    std::vector<double> *v_true_y = &true_vertex_y;
    std::vector<double> *v_true_z = &true_vertex_z;

    outTree->Branch("reco_p", &p_reco);
    outTree->Branch("true_p", &p_true);
    outTree->Branch("reco_vertex_x", &v_reco_x);
    outTree->Branch("reco_vertex_y", &v_reco_y);
    outTree->Branch("reco_vertex_z", &v_reco_z);
    outTree->Branch("true_vertex_x", &v_true_x);
    outTree->Branch("true_vertex_y", &v_true_y);
    outTree->Branch("true_vertex_z", &v_true_z);

    outTree->Fill();
    outFile->Write();
    outFile->Close();
    file->Close();

    std::cout << "Saved kaon momentum and vertex vectors to " << outFileName << "\n";
}
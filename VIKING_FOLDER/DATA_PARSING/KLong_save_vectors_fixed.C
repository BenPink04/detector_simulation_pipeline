#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>

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

void KLong_save_vectors_fixed(const char* filename = "Scenario3_Seed1.root") {
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
            if (pdg == 211) has_pip = true;
            if (pdg == -211) has_pim = true;
            if (pdg == 111) has_pi0 = true;
        }
        if (has_pip && has_pim && has_pi0) {
            selected_events.push_back(kv.first);
        }
    }
    std::cout << "Found " << selected_events.size() << " events with π+, π-, π0 as direct kaon decay products.\n";

    // --- Data storage vectors (using separate coordinate vectors instead of TVector3) ---
    std::vector<double> reco_p, true_p;
    std::vector<double> reco_vertex_x, reco_vertex_y, reco_vertex_z;
    std::vector<double> true_vertex_x, true_vertex_y, true_vertex_z;

    // Detector definitions
    auto is_tracker = [](int id) { return id >= 1 && id <= 4; };
    auto is_pizza = [](int id) { return id >= 5 && id <= 8; };
    auto is_fri = [](int id) { return id >= 9 && id <= 10; };
    auto is_tof = [](int id) { return id == 11; };

    // Smearing parameters
    TRandom3 randGen(0);
    double smear_sigma = 0.01; // 100 μm position smearing
    double smear_time_sigma = 0.1; // 100 ps time smearing

    // --- Step 2: Process hits for selected events ---
    TTree *tree3 = (TTree*)file->Get("Ntuple3");
    if (!tree3) { std::cout << "Tree Ntuple3 not found!\n"; file->Close(); return; }

    Double_t evtNb, x, y, z, t, pdg, deviceID;
    tree3->SetBranchAddress("evtNb", &evtNb);
    tree3->SetBranchAddress("x", &x);
    tree3->SetBranchAddress("y", &y);
    tree3->SetBranchAddress("z", &z);
    tree3->SetBranchAddress("t", &t);
    tree3->SetBranchAddress("particle_PDGEncoding", &pdg);
    tree3->SetBranchAddress("deviceID", &deviceID);

    Long64_t nEntries3 = tree3->GetEntries();
    for (int event_number : selected_events) {
        std::vector<HitInfo> hits_pip, hits_pim;
        double pip_pizza_time = -1, pip_tof_time = -1, pim_pizza_time = -1, pim_tof_time = -1;
        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;
        double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;
        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;
        double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;

        // Collect hits for this event
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

        // Approximate momentum vectors using first two hits per particle
        if (hits_pip.size() >= 2) {
            TVector3 pip_point1(hits_pip[0].x, hits_pip[0].y, hits_pip[0].z);
            TVector3 pip_point2(hits_pip[1].x, hits_pip[1].y, hits_pip[1].z);
            TVector3 pip_direction = (pip_point2 - pip_point1).Unit();

            TVector3 pim_point1(hits_pim[0].x, hits_pim[0].y, hits_pim[0].z);
            TVector3 pim_point2(hits_pim[1].x, hits_pim[1].y, hits_pim[1].z);
            TVector3 pim_direction = (pim_point2 - pim_point1).Unit();

            TVector3 decay_vertex = closest_point_between_lines(pip_point1, pip_direction, pim_point1, pim_direction);

            // TOF-based momentum calculation
            double pip_flight_time = pip_tof_time - pip_pizza_time;
            double pim_flight_time = pim_tof_time - pim_pizza_time;

            if (pip_flight_time <= 0 || pim_flight_time <= 0) continue;

            double pip_flight_distance = (pip_tof_pos - pip_pizza_pos).Mag() / 100.0;
            double pim_flight_distance = (pim_tof_pos - pim_pizza_pos).Mag() / 100.0;

            double pip_velocity = pip_flight_distance / (pip_flight_time * 1e-9);
            double pim_velocity = pim_flight_distance / (pim_flight_time * 1e-9);

            double m_pi = 0.13957; // GeV/c^2
            double beta_pip = pip_velocity / 2.99792458e8;
            double beta_pim = pim_velocity / 2.99792458e8;
            
            if (beta_pip >= 1.0) beta_pip = 0.9999;
            if (beta_pim >= 1.0) beta_pim = 0.9999;
            
            double gamma_pip = 1.0 / std::sqrt(1 - beta_pip*beta_pip);
            double gamma_pim = 1.0 / std::sqrt(1 - beta_pim*beta_pim);
            
            TVector3 pip_momentum = gamma_pip * m_pi * beta_pip * pip_direction;
            TVector3 pim_momentum = gamma_pim * m_pi * beta_pim * pim_direction;

            TVector3 kaon_momentum = pip_momentum + pim_momentum;
            double kaon_p = kaon_momentum.Mag();

            // Get true kaon momentum and vertex from MC truth (Ntuple1)
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

            // Apply reasonable momentum cut and require truth match
            if (kaon_p > 11) continue;
            if (found_truth && true_p_mag != 0) {
                std::cout << "Event " << event_number
                          << " | Reco p: " << kaon_p
                          << " | True p: " << true_p_mag
                          << " | Reco vertex: (" << decay_vertex.X() << ", " << decay_vertex.Y() << ", " << decay_vertex.Z() << ")"
                          << " | True vertex: (" << true_vx << ", " << true_vy << ", " << true_vz << ")"
                          << std::endl;

                // Store data using separate coordinate vectors (avoids TVector3 dictionary issue)
                reco_p.push_back(kaon_p);
                true_p.push_back(true_p_mag);
                reco_vertex_x.push_back(decay_vertex.X());
                reco_vertex_y.push_back(decay_vertex.Y());
                reco_vertex_z.push_back(decay_vertex.Z());
                true_vertex_x.push_back(true_vx);
                true_vertex_y.push_back(true_vy);  
                true_vertex_z.push_back(true_vz);
            }
        }
    }

    // --- Save vectors to file ---
    std::string inFileName(filename);
    size_t lastdot = inFileName.find_last_of(".");
    std::string base = (lastdot == std::string::npos) ? inFileName : inFileName.substr(0, lastdot);
    std::string outFileName = base + "_vectors.root";
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    TTree *outTree = new TTree("kaonVectors", "Kaon momentum and vertex vectors");

    // Use pointers to vectors for ROOT tree branches
    std::vector<double> *p_reco = &reco_p;
    std::vector<double> *p_true = &true_p;
    std::vector<double> *v_reco_x = &reco_vertex_x;
    std::vector<double> *v_reco_y = &reco_vertex_y;
    std::vector<double> *v_reco_z = &reco_vertex_z;
    std::vector<double> *v_true_x = &true_vertex_x;
    std::vector<double> *v_true_y = &true_vertex_y;
    std::vector<double> *v_true_z = &true_vertex_z;

    // Create branches for momentum and vertex coordinates (no TVector3 dictionary needed)
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

    std::cout << "Saved " << reco_p.size() << " kaon momentum and vertex vectors to " << outFileName << std::endl;
}
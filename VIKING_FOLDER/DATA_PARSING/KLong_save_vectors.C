#include "TFile.h"#include "TFile.h"

#include "TTree.h"#include "TTree.h"

#include "TVector3.h"#include "TVector3.h"

#include <vector>#include <vector>

#include <iostream>#include <iostream>

#include <cmath>#include <cmath>

#include <string>#include <string>

#include <map>

struct HitInfo { double x, y, z, t; int deviceID; };

struct HitInfo { double x, y, z, t; int deviceID; };

TVector3 closest_point_between_lines(const TVector3& p1, const TVector3& v1,

TVector3 closest_point_between_lines(const TVector3& p1, const TVector3& v1,                                    const TVector3& p2, const TVector3& v2) {

                                    const TVector3& p2, const TVector3& v2) {    TVector3 w0 = p1 - p2;

    TVector3 w0 = p1 - p2;    double a = v1.Dot(v1);

    double a = v1.Dot(v1);    double b = v1.Dot(v2);

    double b = v1.Dot(v2);    double c = v2.Dot(v2);

    double c = v2.Dot(v2);    double d = v1.Dot(w0);

    double d = v1.Dot(w0);    double e = v2.Dot(w0);

    double e = v2.Dot(w0);    double denom = a*c - b*b;

    double denom = a*c - b*b;    double sc = (b*e - c*d) / denom;

    if (std::abs(denom) < 1e-6) return 0.5 * (p1 + p2); // Lines parallel, return midpoint    double tc = (a*e - b*d) / denom;

    double sc = (b*e - c*d) / denom;    TVector3 point1 = p1 + v1 * sc;

    double tc = (a*e - b*d) / denom;    TVector3 point2 = p2 + v2 * tc;

    TVector3 point1 = p1 + v1 * sc;    return 0.5 * (point1 + point2); // midpoint is PoCA

    TVector3 point2 = p2 + v2 * tc;}

    return 0.5 * (point1 + point2); // midpoint is PoCA

}void KLong_save_vectors(const char* filename = "Scenario3_Seed1.root") {

    std::cout << "Processing vectors for file: " << filename << std::endl;

void KLong_save_vectors(const char* filename = "Scenario3_Seed1.root") {    TFile *file = TFile::Open(filename);

    std::cout << "Processing vectors for file: " << filename << std::endl;    if (!file || file->IsZombie()) { std::cout << "Cannot open root file\n"; return; }

    TFile *file = TFile::Open(filename);

    if (!file || file->IsZombie()) {     // --- Step 1: Find events with triple pion decay (π+, π-, π0) ---

        std::cout << "ERROR: Cannot open root file: " << filename << std::endl;     TTree *tree2 = (TTree*)file->Get("Ntuple2");

        return;     if (!tree2) { std::cout << "Tree Ntuple2 not found!\n"; file->Close(); return; }

    }

    Double_t evtNb2, DKparticle_PDGEncoding;

    // --- Step 1: Find events with triple pion decay (π+, π-, π0) ---    tree2->SetBranchAddress("evtNb", &evtNb2);

    TTree *tree2 = (TTree*)file->Get("Ntuple2");    tree2->SetBranchAddress("DKparticle_PDGEncoding", &DKparticle_PDGEncoding);

    if (!tree2) { 

        std::cout << "ERROR: Tree Ntuple2 not found!" << std::endl;     std::map<int, std::vector<int>> event_products;

        file->Close();     Long64_t nEntries2 = tree2->GetEntries();

        return;     for (Long64_t i = 0; i < nEntries2; ++i) {

    }        tree2->GetEntry(i);

        event_products[(int)evtNb2].push_back((int)DKparticle_PDGEncoding);

    Double_t evtNb2, DKparticle_PDGEncoding;    }

    tree2->SetBranchAddress("evtNb", &evtNb2);

    tree2->SetBranchAddress("DKparticle_PDGEncoding", &DKparticle_PDGEncoding);    std::vector<int> selected_events;

    for (const auto& kv : event_products) {

    std::map<int, std::vector<int>> event_products;        std::vector<int> products = kv.second;

    Long64_t nEntries2 = tree2->GetEntries();        bool has_pim = false, has_pi0 = false, has_pip = false;

            for (int pdg : products) {

    for (Long64_t i = 0; i < nEntries2; ++i) {            if (pdg == -211) has_pim = true;

        tree2->GetEntry(i);            if (pdg == 111) has_pi0 = true;

        event_products[(int)evtNb2].push_back((int)DKparticle_PDGEncoding);            if (pdg == 211) has_pip = true;

    }        }

        if (has_pim && has_pi0 && has_pip)

    std::vector<int> selected_events;            selected_events.push_back(kv.first);

    for (const auto& kv : event_products) {    }

        std::vector<int> products = kv.second;

        bool has_pim = false, has_pi0 = false, has_pip = false;    std::cout << "Found " << selected_events.size() << " events with π+, π-, π0 as direct kaon decay products.\n";

        for (int pdg : products) {

            if (pdg == 211) has_pip = true;

            if (pdg == -211) has_pim = true;    TTree *tree3 = (TTree*)file->Get("Ntuple3");

            if (pdg == 111) has_pi0 = true;    if (!tree3) { std::cout << "Tree Ntuple3 not found!\n"; file->Close(); return; }

        }

        if (has_pip && has_pim && has_pi0) {    Double_t evtNb, x, y, z, t, pdg, deviceID;

            selected_events.push_back(kv.first);    tree3->SetBranchAddress("evtNb", &evtNb);

        }    tree3->SetBranchAddress("hitX", &x);

    }    tree3->SetBranchAddress("hitY", &y);

    std::cout << "Found " << selected_events.size() << " events with π+, π-, π0 as direct kaon decay products." << std::endl;    tree3->SetBranchAddress("hitZ", &z);

    tree3->SetBranchAddress("hitT", &t);

    // --- Data storage vectors (using separate coordinates to avoid TVector3 dictionary issues) ---    tree3->SetBranchAddress("PDGEncoding", &pdg);

    std::vector<double> reco_p, true_p;    tree3->SetBranchAddress("deviceID", &deviceID);

    std::vector<double> reco_vertex_x, reco_vertex_y, reco_vertex_z;

    std::vector<double> true_vertex_x, true_vertex_y, true_vertex_z;    auto is_pizza = [](int id) { return (id >= 489 && id <= 536); };

    auto is_tof   = [](int id) { return (id >= 589 && id <= 606); };

    // Detector definitions    auto is_tracker = [](int id) { return (id >= 1 && id <= 488); };

    auto is_tracker = [](int id) { return id >= 1 && id <= 4; };    auto is_fri = [](int id) { return (id >= 537 && id <= 588); };

    auto is_pizza = [](int id) { return id >= 5 && id <= 8; };

    auto is_fri = [](int id) { return id >= 9 && id <= 10; };    TRandom3 randGen(0);

    auto is_tof = [](int id) { return id == 11; };    double smear_sigma = 5;      // cm

    double smear_time_sigma = 0.0015; // ns

    // --- Step 2: Process hits for selected events ---

    TTree *tree3 = (TTree*)file->Get("Ntuple3");    std::vector<double> reco_p, true_p;

    if (!tree3) {     std::vector<TVector3> reco_vertex, true_vertex;

        std::cout << "ERROR: Tree Ntuple3 not found!" << std::endl; 

        file->Close();     for (int event_number : selected_events) {

        return;         double pip_pizza_time = -1, pip_tof_time = -1;

    }        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;

        double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;

    Double_t evtNb, x, y, z, t, pdg, deviceID;

    tree3->SetBranchAddress("evtNb", &evtNb);        double pim_pizza_time = -1, pim_tof_time = -1;

    tree3->SetBranchAddress("x", &x);        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;

    tree3->SetBranchAddress("y", &y);        double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;

    tree3->SetBranchAddress("z", &z);

    tree3->SetBranchAddress("t", &t);        std::vector<HitInfo> hits_pip, hits_pim;

    tree3->SetBranchAddress("particle_PDGEncoding", &pdg);

    tree3->SetBranchAddress("deviceID", &deviceID);        Long64_t nEntries3 = tree3->GetEntries();

        for (Long64_t i = 0; i < nEntries3; ++i) {

    // Smearing parameters            tree3->GetEntry(i);

    TRandom3 randGen(0);            if ((int)evtNb != event_number) continue;

    double smear_sigma = 0.01; // 100 μm position smearing            if ((int)pdg == 211) {

    double smear_time_sigma = 0.1; // 100 ps time smearing                if (is_tracker((int)deviceID) || is_fri((int)deviceID))

                    hits_pip.push_back({x, y, z, t, (int)deviceID});

    Long64_t nEntries3 = tree3->GetEntries();                if (is_pizza((int)deviceID)) {

                        double smeared_t = randGen.Gaus(t, smear_time_sigma);

    for (int event_number : selected_events) {                    if (pip_pizza_time < 0 || smeared_t < pip_pizza_time) {

        std::vector<HitInfo> hits_pip, hits_pim;                        pip_pizza_time = smeared_t;

        double pip_pizza_time = -1, pip_tof_time = -1, pim_pizza_time = -1, pim_tof_time = -1;                        pip_pizza_x = x; pip_pizza_y = y; pip_pizza_z = z;

        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;                    }

        double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;                }

        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;                if (is_tof((int)deviceID)) {

        double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;                    double smeared_t = randGen.Gaus(t, smear_time_sigma);

                    if (pip_tof_time < 0 || smeared_t < pip_tof_time) {

        // Collect hits for this event                        pip_tof_time = smeared_t;

        for (Long64_t i = 0; i < nEntries3; ++i) {                        pip_tof_x = x; pip_tof_y = y; pip_tof_z = z;

            tree3->GetEntry(i);                    }

            if ((int)evtNb != event_number) continue;                }

                        }

            if ((int)pdg == 211) {  // π+ hits            if ((int)pdg == -211) {

                if (is_tracker((int)deviceID) || is_fri((int)deviceID))                if (is_tracker((int)deviceID) || is_fri((int)deviceID))

                    hits_pip.push_back({x, y, z, t, (int)deviceID});                    hits_pim.push_back({x, y, z, t, (int)deviceID});

                if (is_pizza((int)deviceID)) {                if (is_pizza((int)deviceID)) {

                    double smeared_t = randGen.Gaus(t, smear_time_sigma);                    double smeared_t = randGen.Gaus(t, smear_time_sigma);

                    if (pip_pizza_time < 0 || smeared_t < pip_pizza_time) {                    if (pim_pizza_time < 0 || smeared_t < pim_pizza_time) {

                        pip_pizza_time = smeared_t;                        pim_pizza_time = smeared_t;

                        pip_pizza_x = x; pip_pizza_y = y; pip_pizza_z = z;                        pim_pizza_x = x; pim_pizza_y = y; pim_pizza_z = z;

                    }                    }

                }                }

                if (is_tof((int)deviceID)) {                if (is_tof((int)deviceID)) {

                    double smeared_t = randGen.Gaus(t, smear_time_sigma);                    double smeared_t = randGen.Gaus(t, smear_time_sigma);

                    if (pip_tof_time < 0 || smeared_t < pip_tof_time) {                    if (pim_tof_time < 0 || smeared_t < pim_tof_time) {

                        pip_tof_time = smeared_t;                        pim_tof_time = smeared_t;

                        pip_tof_x = x; pip_tof_y = y; pip_tof_z = z;                        pim_tof_x = x; pim_tof_y = y; pim_tof_z = z;

                    }                    }

                }                }

            }            }

            if ((int)pdg == -211) {  // π- hits        }

                if (is_tracker((int)deviceID) || is_fri((int)deviceID))

                    hits_pim.push_back({x, y, z, t, (int)deviceID});        if (hits_pip.size() < 2 || hits_pim.size() < 2) continue;

                if (is_pizza((int)deviceID)) {        if (pip_pizza_time < 0 || pip_tof_time < 0 || pim_pizza_time < 0 || pim_tof_time < 0)

                    double smeared_t = randGen.Gaus(t, smear_time_sigma);            continue;

                    if (pim_pizza_time < 0 || smeared_t < pim_pizza_time) {

                        pim_pizza_time = smeared_t;        TVector3 pip_pizza_pos(

                        pim_pizza_x = x; pim_pizza_y = y; pim_pizza_z = z;            randGen.Gaus(pip_pizza_x, smear_sigma),

                    }            randGen.Gaus(pip_pizza_y, smear_sigma),

                }            pip_pizza_z

                if (is_tof((int)deviceID)) {        );

                    double smeared_t = randGen.Gaus(t, smear_time_sigma);        TVector3 pim_pizza_pos(

                    if (pim_tof_time < 0 || smeared_t < pim_tof_time) {            randGen.Gaus(pim_pizza_x, smear_sigma),

                        pim_tof_time = smeared_t;            randGen.Gaus(pim_pizza_y, smear_sigma),

                        pim_tof_x = x; pim_tof_y = y; pim_tof_z = z;            pim_pizza_z

                    }        );

                }

            }        TVector3 pip_start(hits_pip.front().x, hits_pip.front().y, hits_pip.front().z);

        }        TVector3 pip_end(hits_pip.back().x, hits_pip.back().y, hits_pip.back().z);

        TVector3 pip_dir = (pip_end - pip_start).Unit();

        // Require minimum hits and timing info

        if (hits_pip.size() < 2 || hits_pim.size() < 2) continue;        TVector3 pim_start(hits_pim.front().x, hits_pim.front().y, hits_pim.front().z);

        if (pip_pizza_time < 0 || pip_tof_time < 0 || pim_pizza_time < 0 || pim_tof_time < 0) continue;        TVector3 pim_end(hits_pim.back().x, hits_pim.back().y, hits_pim.back().z);

        TVector3 pim_dir = (pim_end - pim_start).Unit();

        // Smear detector positions

        TVector3 pip_pizza_pos(        TVector3 decay_vertex = closest_point_between_lines(pip_start, pip_dir, pim_start, pim_dir);

            randGen.Gaus(pip_pizza_x, smear_sigma),

            randGen.Gaus(pip_pizza_y, smear_sigma),        TVector3 pip_tof_pos(

            pip_pizza_z            randGen.Gaus(pip_tof_x, smear_sigma),

        );            randGen.Gaus(pip_tof_y, smear_sigma),

        TVector3 pim_pizza_pos(            pip_tof_z

            randGen.Gaus(pim_pizza_x, smear_sigma),        );

            randGen.Gaus(pim_pizza_y, smear_sigma),        TVector3 pim_tof_pos(

            pim_pizza_z            randGen.Gaus(pim_tof_x, smear_sigma),

        );            randGen.Gaus(pim_tof_y, smear_sigma),

        TVector3 pip_tof_pos(            pim_tof_z

            randGen.Gaus(pip_tof_x, smear_sigma),        );

            randGen.Gaus(pip_tof_y, smear_sigma),

            pip_tof_z        double pip_track_length_cm = (pip_tof_pos - pip_pizza_pos).Mag(); // cm

        );        double pim_track_length_cm = (pim_tof_pos - pim_pizza_pos).Mag(); // cm

        TVector3 pim_tof_pos(

            randGen.Gaus(pim_tof_x, smear_sigma),        double pip_dt_ns = pip_tof_time - pip_pizza_time;

            randGen.Gaus(pim_tof_y, smear_sigma),        double pim_dt_ns = pim_tof_time - pim_pizza_time;

            pim_tof_z

        );        double pip_v = (pip_dt_ns > 0) ? (pip_track_length_cm * 1e-2) / (pip_dt_ns * 1e-9) : 0; // m/s

        double pim_v = (pim_dt_ns > 0) ? (pim_track_length_cm * 1e-2) / (pim_dt_ns * 1e-9) : 0; // m/s

        // Calculate track directions using first two hits

        TVector3 pip_point1(hits_pip[0].x, hits_pip[0].y, hits_pip[0].z);        double pip_path_cm = (pip_pizza_pos - decay_vertex).Mag();

        TVector3 pip_point2(hits_pip[1].x, hits_pip[1].y, hits_pip[1].z);        double pip_dt_s = pip_path_cm * 1e-2 / pip_v; // s

        TVector3 pip_direction = (pip_point2 - pip_point1).Unit();        double pip_decay_time = pip_pizza_time * 1e-9 - pip_dt_s; // s



        TVector3 pim_point1(hits_pim[0].x, hits_pim[0].y, hits_pim[0].z);        double pim_path_cm = (pim_pizza_pos - decay_vertex).Mag();

        TVector3 pim_point2(hits_pim[1].x, hits_pim[1].y, hits_pim[1].z);        double pim_dt_s = pim_path_cm * 1e-2 / pim_v; // s

        TVector3 pim_direction = (pim_point2 - pim_point1).Unit();        double pim_decay_time = pim_pizza_time * 1e-9 - pim_dt_s; // s



        // Find decay vertex        double kaon_decay_time = 0.5 * (pip_decay_time + pim_decay_time); // s

        TVector3 decay_vertex = closest_point_between_lines(pip_point1, pip_direction, pim_point1, pim_direction);

        TVector3 kaon_prod(0,0,0);

        // Calculate velocities from TOF        double kaon_flight_length_cm = (decay_vertex - kaon_prod).Mag();

        double pip_flight_time = pip_tof_time - pip_pizza_time;        double kaon_flight_time_s = kaon_decay_time; // assuming kaon produced at t=0

        double pim_flight_time = pim_tof_time - pim_pizza_time;

        double kaon_velocity = (kaon_flight_length_cm * 1e-2) / kaon_flight_time_s; // m/s

        if (pip_flight_time <= 0 || pim_flight_time <= 0) continue;

        double m_K = 0.497611; // GeV/c^2

        double pip_flight_distance = (pip_tof_pos - pip_pizza_pos).Mag() / 100.0; // convert to meters        double beta_K = kaon_velocity / 2.99792458e8;

        double pim_flight_distance = (pim_tof_pos - pim_pizza_pos).Mag() / 100.0; // convert to meters        if (beta_K >= 1.0) beta_K = 0.9999;

        double gamma_K = 1.0 / std::sqrt(1 - beta_K*beta_K);

        double pip_velocity = pip_flight_distance / (pip_flight_time * 1e-9); // m/s        double kaon_p = gamma_K * m_K * beta_K; // GeV/c

        double pim_velocity = pim_flight_distance / (pim_flight_time * 1e-9); // m/s

        // --- True kaon decay momentum and vertex from MC truth (Ntuple1) ---

        // Calculate momenta        TTree *tree1 = (TTree*)file->Get("Ntuple1");

        double m_pi = 0.13957; // GeV/c^2        double true_px = 0, true_py = 0, true_pz = 0;

        double c = 2.99792458e8; // m/s        double true_vx = 0, true_vy = 0, true_vz = 0;

                bool found_truth = false;

        double beta_pip = pip_velocity / c;        if (tree1) {

        double beta_pim = pim_velocity / c;            Double_t evtNb1, kaonDK_momX, kaonDK_momY, kaonDK_momZ;

                    Double_t kaonDK_posX, kaonDK_posY, kaonDK_posZ;

        if (beta_pip >= 1.0) beta_pip = 0.9999;            tree1->SetBranchAddress("evtNb", &evtNb1);

        if (beta_pim >= 1.0) beta_pim = 0.9999;            tree1->SetBranchAddress("kaonDK_momX", &kaonDK_momX);

                    tree1->SetBranchAddress("kaonDK_momY", &kaonDK_momY);

        double gamma_pip = 1.0 / std::sqrt(1 - beta_pip*beta_pip);            tree1->SetBranchAddress("kaonDK_momZ", &kaonDK_momZ);

        double gamma_pim = 1.0 / std::sqrt(1 - beta_pim*beta_pim);            tree1->SetBranchAddress("kaonDK_posX", &kaonDK_posX);

                    tree1->SetBranchAddress("kaonDK_posY", &kaonDK_posY);

        TVector3 pip_momentum = gamma_pip * m_pi * beta_pip * pip_direction;            tree1->SetBranchAddress("kaonDK_posZ", &kaonDK_posZ);

        TVector3 pim_momentum = gamma_pim * m_pi * beta_pim * pim_direction;

            Long64_t nEntries1 = tree1->GetEntries();

        // Reconstructed kaon momentum            for (Long64_t i = 0; i < nEntries1; ++i) {

        TVector3 kaon_momentum = pip_momentum + pim_momentum;                tree1->GetEntry(i);

        double kaon_p = kaon_momentum.Mag();                if ((int)evtNb1 == event_number) {

                    true_px = kaonDK_momX;

        // Get true kaon momentum and vertex from MC truth (Ntuple1)                    true_py = kaonDK_momY;

        TTree *tree1 = (TTree*)file->Get("Ntuple1");                    true_pz = kaonDK_momZ;

        double true_px = 0, true_py = 0, true_pz = 0;                    true_vx = kaonDK_posX;

        double true_vx = 0, true_vy = 0, true_vz = 0;                    true_vy = kaonDK_posY;

        bool found_truth = false;                    true_vz = kaonDK_posZ;

                            found_truth = true;

        if (tree1) {                    break;

            Double_t evtNb1, kaonDK_momX, kaonDK_momY, kaonDK_momZ;                }

            Double_t kaonDK_posX, kaonDK_posY, kaonDK_posZ;            }

            tree1->SetBranchAddress("evtNb", &evtNb1);        }

            tree1->SetBranchAddress("kaonDK_momX", &kaonDK_momX);

            tree1->SetBranchAddress("kaonDK_momY", &kaonDK_momY);        double true_p_mag = std::sqrt(true_px*true_px + true_py*true_py + true_pz*true_pz);

            tree1->SetBranchAddress("kaonDK_momZ", &kaonDK_momZ);        TVector3 true_vertex_vec(true_vx, true_vy, true_vz);

            tree1->SetBranchAddress("kaonDK_posX", &kaonDK_posX);

            tree1->SetBranchAddress("kaonDK_posY", &kaonDK_posY);        if (kaon_p > 11) continue;

            tree1->SetBranchAddress("kaonDK_posZ", &kaonDK_posZ);        if (found_truth && true_p_mag != 0) {

            // Print as soon as the event is processed

            Long64_t nEntries1 = tree1->GetEntries();            std::cout << "Event " << event_number

            for (Long64_t i = 0; i < nEntries1; ++i) {                      << " | Reco p: " << kaon_p

                tree1->GetEntry(i);                      << " | True p: " << true_p_mag

                if ((int)evtNb1 == event_number) {                      << " | Reco vertex: (" << decay_vertex.X() << ", " << decay_vertex.Y() << ", " << decay_vertex.Z() << ")"

                    true_px = kaonDK_momX;                      << " | True vertex: (" << true_vertex_vec.X() << ", " << true_vertex_vec.Y() << ", " << true_vertex_vec.Z() << ")"

                    true_py = kaonDK_momY;                      << std::endl;

                    true_pz = kaonDK_momZ;

                    true_vx = kaonDK_posX;            reco_p.push_back(kaon_p);

                    true_vy = kaonDK_posY;            true_p.push_back(true_p_mag);

                    true_vz = kaonDK_posZ;            reco_vertex.push_back(decay_vertex);

                    found_truth = true;            true_vertex.push_back(true_vertex_vec);

                    break;        }

                }    }

            }

        }    // --- Save vectors to file ---

    std::string inFileName(filename);

        double true_p_mag = std::sqrt(true_px*true_px + true_py*true_py + true_pz*true_pz);    size_t lastdot = inFileName.find_last_of(".");

    std::string base = (lastdot == std::string::npos) ? inFileName : inFileName.substr(0, lastdot);

        // Apply cuts and store data    std::string outFileName = base + "_vectors.root";

        if (kaon_p > 11) continue; // Momentum cut    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");

        if (found_truth && true_p_mag != 0) {    TTree *outTree = new TTree("kaonVectors", "Kaon momentum and vertex vectors");

            std::cout << "Event " << event_number

                      << " | Reco p: " << kaon_p    std::vector<double> *p_reco = &reco_p;

                      << " | True p: " << true_p_mag << std::endl;    std::vector<double> *p_true = &true_p;

    std::vector<TVector3> *v_reco = &reco_vertex;

            // Store data using separate coordinate vectors (avoids TVector3 dictionary issue)    std::vector<TVector3> *v_true = &true_vertex;

            reco_p.push_back(kaon_p);

            true_p.push_back(true_p_mag);    outTree->Branch("reco_p", &p_reco);

            reco_vertex_x.push_back(decay_vertex.X());    outTree->Branch("true_p", &p_true);

            reco_vertex_y.push_back(decay_vertex.Y());    outTree->Branch("reco_vertex", &v_reco);

            reco_vertex_z.push_back(decay_vertex.Z());    outTree->Branch("true_vertex", &v_true);

            true_vertex_x.push_back(true_vx);

            true_vertex_y.push_back(true_vy);    outTree->Fill();

            true_vertex_z.push_back(true_vz);    outFile->Write();

        }    outFile->Close();

    }    file->Close();



    // --- Save vectors to file ---    std::cout << "Saved kaon momentum and vertex vectors to kaon_vectors.root\n";

    std::string inFileName(filename);}
    size_t lastdot = inFileName.find_last_of(".");
    std::string base = (lastdot == std::string::npos) ? inFileName : inFileName.substr(0, lastdot);
    std::string outFileName = base + "_vectors.root";
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    TTree *outTree = new TTree("kaonVectors", "Kaon momentum and vertex vectors");

    // Use pointers to vectors for ROOT tree branches (fixed version)
    std::vector<double> *p_reco = &reco_p;
    std::vector<double> *p_true = &true_p;
    std::vector<double> *v_reco_x = &reco_vertex_x;
    std::vector<double> *v_reco_y = &reco_vertex_y;
    std::vector<double> *v_reco_z = &reco_vertex_z;
    std::vector<double> *v_true_x = &true_vertex_x;
    std::vector<double> *v_true_y = &true_vertex_y;
    std::vector<double> *v_true_z = &true_vertex_z;

    // Create branches (no TVector3 - avoids dictionary compilation issues)
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

    std::cout << "SUCCESS: Saved " << reco_p.size() << " kaon momentum and vertex vectors to " << outFileName << std::endl;
}
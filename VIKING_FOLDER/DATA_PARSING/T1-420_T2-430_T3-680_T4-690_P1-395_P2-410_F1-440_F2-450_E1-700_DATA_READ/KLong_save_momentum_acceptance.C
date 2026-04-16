// ============================================================================
// KLong_save_momentum_acceptance.C  —  3D TGRAPH TRACK RECONSTRUCTION
//
// PURPOSE:
//   Same reconstruction as KLong_save_vectors.C, but records ALL events that
//   have the correct decay channel (pi+ pi- pi0), not only reconstructed ones.
//   A per-event reco_flag (1=reconstructed, 0=not) is saved alongside the
//   true momentum, enabling acceptance/efficiency calculations.
//
//   Track fitting uses TGraphErrors x(z) and y(z) fitted with pol1 (straight
//   line), applied to the raw Geant4 hit positions.  See KLong_save_vectors.C
//   and TRACKING_OVERHAUL_README.md for full geometric and algorithmic details.
//
// INPUT:  a simulation .root file containing Ntuple1, Ntuple2, Ntuple3
//         Filename must contain detector position tags, e.g.:
//           T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700
//
// OUTPUT: <input_base>_acceptance.root  (one TTree named "kaonEventInfo")
//         Branches:
//           n_triple_pion_events  : total events with pi+pi-pi0 decay
//           true_mom_vec          : vector of true kaon momenta (GeV/c)
//           reco_flag_vec         : vector of ints (1=reco success, 0=failure)
//
// SEE ALSO:
//   KLong_save_momentum_acceptance_pixelation.C — DSSSD/pixel alternative
//   KLong_save_vectors.C                        — reco-only companion macro
//   TRACKING_OVERHAUL_README.md                 — geometry and methodology
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "ROOT/RDataFrame.hxx"
#include <chrono>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <set>

// ============================================================================
// GEOMETRY CONSTANTS  (taken directly from JLabKDetectorConstruction.cc)
// ============================================================================

const double STRAW_HALF_WIDTH  = 0.40;   // cm
const double STRAW_HALF_LENGTH = 45.0;   // cm

const double FRI_STRIP_X[26] = {
    -66., -60., -54., -48., -42., -36., -30., -24., -18., -12., -7.5,
     -4.5, -1.5,  1.5,  4.5,  7.5,  12.,  18.,  24.,  30.,  36.,  42.,
     48.,  54.,  60.,  66.
};

const double FRI_HALF_LENGTH[26] = {
    37.25, 42.75, 53.75, 59.25, 59.25, 64.75,
    70.25, 70.25, 70.25, 70.25,
    70.25, 70.25, 70.25, 70.25, 70.25, 70.25,
    70.25, 70.25, 70.25, 70.25,
    64.75, 59.25, 59.25, 53.75, 42.75, 37.25
};

const double TOF_BAR_HALF_WIDTH = 6.0;   // cm
const double TOF_Y_UNCERTAINTY  = 10.0;  // cm
const double SQRT2 = 1.41421356237;

// ============================================================================
// STRUCTS
// ============================================================================

struct HitInfo { double x, y, z, t; int deviceID; };

struct TruthInfo {
    double px, py, pz;
    double vx, vy, vz;
};

struct EventReco {
    std::vector<HitInfo> hits_pip;
    std::vector<HitInfo> hits_pim;

    bool   has_pip_pizza = false;
    double pip_pizza_time = 0;
    double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;

    bool   has_pip_tof = false;
    double pip_tof_time = 0;
    int    pip_tof_deviceID = -1;
    double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;

    bool   has_pim_pizza = false;
    double pim_pizza_time = 0;
    double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;

    bool   has_pim_tof = false;
    double pim_tof_time = 0;
    int    pim_tof_deviceID = -1;
    double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;
};

// Result of a 3D TGraph straight-line track fit.
struct TrackFit {
    TVector3 dir;
    TVector3 orig;
    bool     valid;
};

// ============================================================================
// GEOMETRY HELPER FUNCTIONS
// ============================================================================

double fri_half_width(int strip_i) {
    return (strip_i >= 10 && strip_i <= 15) ? 1.5 : 3.0;
}

double tof_bar_x_centre(int copyNumber) {
    int i = copyNumber - 2053;
    return -96.0 + i * 12.0;
}

double tof_bar_y_centre(int copyNumber) {
    int i = copyNumber - 2053;
    if (i == 8) return -72.0;
    if (i == 9) return +72.0;
    return 0.0;
}

// ============================================================================
// 3D TGRAPH TRACK FIT
//
// Builds TGraphErrors x(z) and y(z) from raw Geant4 hit positions and fits
// pol1 to each.  Returns the track direction and reference point at z=0.
// See KLong_save_vectors.C for detailed documentation of the algorithm.
//
// Detector assignments:
//   T1  (1-488)    : measures X, sigma = 0.04 cm (effective resolution better than half-width due to time info)
//   T2  (489-976)  : measures Y, sigma = 0.04 cm (effective resolution better than half-width due to time info)
//   T3  (977-1464) : stereo +45, both axes, sigma = 0.04 * sqrt(2) cm (propagated)
//   T4  (1465-1952): stereo -45, both axes, sigma = 0.04 * sqrt(2) cm (propagated)
//   FRI-W1 (2001-2026): measures X, sigma = fri_half_width(strip_i) * 0.5 cm (uniform distribution across strip width -> sigma = half-width * 0.5)
//   FRI-W2 (2027-2052): measures Y, sigma = fri_half_width(strip_i) * 0.5 cm (uniform distribution across strip width -> sigma = half-width * 0.5)
//   TOF (via tof_devID): bar centre X/Y added with geometry-derived sigmas
//
// All four sub-layers of each tracker station contribute hit points to the fit.
// Requires >= 2 points at >= 2 distinct z-planes in each projection.
// ============================================================================
TrackFit fit_3dgraph_track(const std::vector<HitInfo>& hits,
                           int    tof_devID,
                           double z_tof)
{
    const TrackFit FAIL = {{0,0,1}, {0,0,0}, false};

    std::vector<double> xz, xv, xe;
    std::vector<double> yz, yv, ye;

    for (const auto& h : hits) {
        int id = h.deviceID;

        if (id >= 1 && id <= 488) {
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(0.04); // cm  — effective resolution better than half-width due to time info
        }
        else if (id >= 489 && id <= 976) {
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(0.04); // cm  — effective resolution better than half-width due to time info
        }
        else if (id >= 977 && id <= 1464) {
            double sxy = 0.04 * SQRT2;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(sxy);
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(sxy);
        }
        else if (id >= 1465 && id <= 1952) {
            double sxy = 0.04 * SQRT2;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(sxy);
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(sxy);
        }
        else if (id >= 2001 && id <= 2026) {
            int strip_i = id - 2001;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(fri_half_width(strip_i) * 0.5); // cm  — uniform distribution across strip width -> sigma = half-width * 0.5
        }
        else if (id >= 2027 && id <= 2052) {
            int strip_i = id - 2027;
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(fri_half_width(strip_i) * 0.5); // cm  — uniform distribution across strip width -> sigma = half-width * 0.5
        }
        // TOF hits are intentionally excluded from the track fit.
        // The bar centre has large uncertainties (±6 cm in X, ±10 cm in Y)
        // which, at the ~460 cm lever arm from the trackers, bias the fitted
        // slope and degrade backward-extrapolated vertex resolution.
        // The TOF hit is still used for the velocity measurement.
    }

    // tof_devID and z_tof kept in signature for call-site compatibility.
    (void)tof_devID; (void)z_tof;

    if ((int)xz.size() < 2 || (int)yz.size() < 2) return FAIL;

    // Check for at least 2 distinct z-planes in each fit
    auto has_distinct_z = [](const std::vector<double>& zv) -> bool {
        for (size_t k = 1; k < zv.size(); ++k)
            if (std::fabs(zv[k] - zv[0]) > 1.0) return true;
        return false;
    };
    if (!has_distinct_z(xz) || !has_distinct_z(yz)) return FAIL;

    // Fit x(z) = p0 + p1*z
    TGraphErrors gx((int)xz.size(), xz.data(), xv.data(), nullptr, xe.data());
    gx.Fit("pol1", "Q");
    TF1* fx = gx.GetFunction("pol1");
    if (!fx) return FAIL;
    double bx = fx->GetParameter(0);
    double ax = fx->GetParameter(1);

    // Fit y(z) = p0 + p1*z
    TGraphErrors gy((int)yz.size(), yz.data(), yv.data(), nullptr, ye.data());
    gy.Fit("pol1", "Q");
    TF1* fy = gy.GetFunction("pol1");
    if (!fy) return FAIL;
    double by = fy->GetParameter(0);
    double ay = fy->GetParameter(1);

    if (!std::isfinite(ax) || !std::isfinite(ay) ||
        !std::isfinite(bx) || !std::isfinite(by)) return FAIL;

    TVector3 dir(ax, ay, 1.0);
    dir = dir.Unit();
    TVector3 orig(bx, by, 0.);
    return {dir, orig, true};
}

// ============================================================================
// POINT OF CLOSEST APPROACH (PoCA) between two straight 3D lines.
// ============================================================================
TVector3 closest_point_between_lines(const TVector3& p1, const TVector3& v1,
                                      const TVector3& p2, const TVector3& v2)
{
    TVector3 w0 = p1 - p2;
    double a = v1.Dot(v1), b = v1.Dot(v2), c = v2.Dot(v2);
    double d = v1.Dot(w0), e = v2.Dot(w0);
    double denom = a*c - b*b;
    if (std::abs(denom) < 1e-14) return 0.5 * (p1 + p2);
    double sc = (b*e - c*d) / denom;
    double tc = (a*e - b*d) / denom;
    return 0.5 * ((p1 + v1*sc) + (p2 + v2*tc));
}

// ===========================================================================
// MAIN FUNCTION
// [TO CHANGE] Update the default filename for a different scenario or seed.
// ===========================================================================
void KLong_save_momentum_acceptance(const char* filename = "Scenario3_Seed1.root") {

    ROOT::EnableImplicitMT(4);
    auto start_time = std::chrono::steady_clock::now();

    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) { std::cout << "Cannot open root file\n"; return; }

    // ========================================================================
    // PARSE DETECTOR Z POSITIONS FROM FILENAME
    // [TO CHANGE] If filename format changes, update extract_param below.
    // ========================================================================
    auto extract_param = [](const std::string& fname,
                             const std::string& key) -> double {
        size_t pos = fname.find(key + "-");
        if (pos == std::string::npos) return 0.0;
        pos += key.size() + 1;
        size_t end = pos;
        while (end < fname.size() && std::isdigit((unsigned char)fname[end])) end++;
        if (end == pos) return 0.0;
        return std::stod(fname.substr(pos, end - pos));
    };

    std::string fname_str(filename);
    double z_tof = extract_param(fname_str, "E1");

    std::cout << "Detector Z positions:\n"
              << "  T1=" << extract_param(fname_str,"T1")
              << "  T2=" << extract_param(fname_str,"T2")
              << "  T3=" << extract_param(fname_str,"T3")
              << "  T4=" << extract_param(fname_str,"T4") << " cm\n"
              << "  F1=" << extract_param(fname_str,"F1")
              << "  F2=" << extract_param(fname_str,"F2") << " cm\n"
              << "  E1=" << z_tof << " cm\n";

    // ========================================================================
    // STEP 1: SELECT EVENTS WITH K_L -> pi+ pi- pi0
    // [TO CHANGE] For a different decay channel, update PDG codes below.
    // ========================================================================
    TTree *tree2 = (TTree*)file->Get("Ntuple2");
    if (!tree2) { std::cout << "Tree Ntuple2 not found!\n"; file->Close(); return; }

    Double_t evtNb2, DKparticle_PDGEncoding;
    tree2->SetBranchAddress("evtNb",                  &evtNb2);
    tree2->SetBranchAddress("DKparticle_PDGEncoding", &DKparticle_PDGEncoding);

    std::map<int, std::vector<int>> event_products;
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        event_products[(int)evtNb2].push_back((int)DKparticle_PDGEncoding);
    }

    std::vector<int> selected_events;
    for (const auto& kv : event_products) {
        bool has_pip = false, has_pim = false, has_pi0 = false;
        for (int code : kv.second) {
            if (code ==  211) has_pip = true;
            if (code == -211) has_pim = true;
            if (code ==  111) has_pi0 = true;
        }
        if (has_pip && has_pim && has_pi0)
            selected_events.push_back(kv.first);
    }

    size_t total_events = selected_events.size();
    std::cout << "Found " << total_events
              << " events with pi+, pi-, pi0 as direct kaon decay products.\n"
              << "Processing " << total_events << " selected events...\n";

    std::unordered_set<int> selected_event_set(
        selected_events.begin(), selected_events.end());

    // ========================================================================
    // STEP 2: SET UP Ntuple3
    // ========================================================================
    TTree *tree3 = (TTree*)file->Get("Ntuple3");
    if (!tree3) { std::cout << "Tree Ntuple3 not found!\n"; file->Close(); return; }

    Double_t evtNb, x, y, z, t, pdg, devID;
    tree3->SetBranchAddress("evtNb",       &evtNb);
    tree3->SetBranchAddress("hitX",        &x);
    tree3->SetBranchAddress("hitY",        &y);
    tree3->SetBranchAddress("hitZ",        &z);
    tree3->SetBranchAddress("hitT",        &t);
    tree3->SetBranchAddress("PDGEncoding", &pdg);
    tree3->SetBranchAddress("deviceID",    &devID);

    // -----------------------------------------------------------------------
    // DEVICE ID CLASSIFICATION
    // [TO CHANGE] Update ranges if geometry changes.
    // -----------------------------------------------------------------------
    auto is_tracker  = [](int id) { return (id >=    1 && id <=  1952); };
    auto is_fri      = [](int id) { return (id >= 2001 && id <=  2052); };
    auto is_pizza    = [](int id) { return (id >= 1953 && id <=  2000); };
    auto is_tof      = [](int id) { return (id >= 2053 && id <=  2070); };

    // -----------------------------------------------------------------------
    // RESOLUTION PARAMETERS
    // [TO CHANGE] Tune smear_time_sigma for your timing resolution.
    // -----------------------------------------------------------------------
    TRandom3 randGen(0);
    double smear_sigma      = 3.0;    // cm  — PIZZA x/y position smear (~1-bar-width uncertainty, simulates real detector response)
    double smear_time_sigma = 0.0015; // ns  — 1.5 ps timing smear

    // Output storage
    std::vector<double> reco_p, true_p;
    std::vector<double> reco_vertex_x, reco_vertex_y, reco_vertex_z;
    std::vector<double> true_vertex_x, true_vertex_y, true_vertex_z;
    std::vector<double> all_true_p;
    std::vector<int>    all_reco_flags;

    // ========================================================================
    // STEP 3: BUILD MC TRUTH LOOKUP FROM Ntuple1
    // ========================================================================
    std::unordered_map<int, TruthInfo> truth_map;
    TTree *tree1 = (TTree*)file->Get("Ntuple1");
    if (tree1) {
        Double_t evtNb1, kpx, kpy, kpz, kvx, kvy, kvz;
        tree1->SetBranchAddress("evtNb",       &evtNb1);
        tree1->SetBranchAddress("kaonDK_momX", &kpx);
        tree1->SetBranchAddress("kaonDK_momY", &kpy);
        tree1->SetBranchAddress("kaonDK_momZ", &kpz);
        tree1->SetBranchAddress("kaonDK_posX", &kvx);
        tree1->SetBranchAddress("kaonDK_posY", &kvy);
        tree1->SetBranchAddress("kaonDK_posZ", &kvz);
        Long64_t nEntries1 = tree1->GetEntries();
        for (Long64_t i = 0; i < nEntries1; ++i) {
            tree1->GetEntry(i);
            truth_map[(int)evtNb1] = { kpx, kpy, kpz, kvx, kvy, kvz };
        }
    }

    // ========================================================================
    // STEP 4: BUILD HIT LOOKUP FROM Ntuple3
    // ========================================================================
    std::unordered_map<int, EventReco> event_data;
    Long64_t nEntries3 = tree3->GetEntries();
    for (Long64_t i = 0; i < nEntries3; ++i) {
        tree3->GetEntry(i);
        int evt_id = (int)evtNb;
        if (selected_event_set.find(evt_id) == selected_event_set.end()) continue;

        auto &ev = event_data[evt_id];
        int id = (int)devID;

        if ((int)pdg == 211) {
            if (is_tracker(id) || is_fri(id))
                ev.hits_pip.push_back({x, y, z, t, id});
            if (is_pizza(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pip_pizza || st < ev.pip_pizza_time) {
                    ev.has_pip_pizza = true; ev.pip_pizza_time = st;
                    ev.pip_pizza_x = randGen.Gaus(x, smear_sigma);
                    ev.pip_pizza_y = randGen.Gaus(y, smear_sigma);
                    ev.pip_pizza_z = z;
                }
            }
            if (is_tof(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pip_tof || st < ev.pip_tof_time) {
                    ev.has_pip_tof = true; ev.pip_tof_time = st;
                    ev.pip_tof_deviceID = id;
                    ev.pip_tof_x = randGen.Gaus(tof_bar_x_centre(id), TOF_BAR_HALF_WIDTH);
                    ev.pip_tof_y = randGen.Gaus(tof_bar_y_centre(id), TOF_Y_UNCERTAINTY);
                    ev.pip_tof_z = z; // raw Geant4 hit z (not geometry face z_tof)
                }
            }
        }
        if ((int)pdg == -211) {
            if (is_tracker(id) || is_fri(id))
                ev.hits_pim.push_back({x, y, z, t, id});
            if (is_pizza(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pim_pizza || st < ev.pim_pizza_time) {
                    ev.has_pim_pizza = true; ev.pim_pizza_time = st;
                    ev.pim_pizza_x = randGen.Gaus(x, smear_sigma);
                    ev.pim_pizza_y = randGen.Gaus(y, smear_sigma);
                    ev.pim_pizza_z = z;
                }
            }
            if (is_tof(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pim_tof || st < ev.pim_tof_time) {
                    ev.has_pim_tof = true; ev.pim_tof_time = st;
                    ev.pim_tof_deviceID = id;
                    ev.pim_tof_x = randGen.Gaus(tof_bar_x_centre(id), TOF_BAR_HALF_WIDTH);
                    ev.pim_tof_y = randGen.Gaus(tof_bar_y_centre(id), TOF_Y_UNCERTAINTY);
                    ev.pim_tof_z = z; // raw Geant4 hit z (not geometry face z_tof)
                }
            }
        }
    }

    // ========================================================================
    // STEP 5: RECONSTRUCT — EVENT LOOP (records ALL selected events)
    // ========================================================================
    int event_counter = 0;
    for (int event_number : selected_events) {
        event_counter++;
        if (event_counter == 1 || event_counter % 1000 == 0
                || event_counter == (int)total_events) {
            auto now = std::chrono::steady_clock::now();
            double elapsed_s =
                std::chrono::duration<double>(now - start_time).count();
            double pct = total_events > 0
                ? (100.0 * event_counter / total_events) : 0.0;
            std::cout << "Progress " << event_counter << "/" << total_events
                      << " (" << pct << "%) - elapsed " << elapsed_s << " s\n";
        }

        // --- MC truth ---
        double true_px = 0, true_py = 0, true_pz = 0;
        double true_vx = 0, true_vy = 0, true_vz = 0;
        bool found_truth = false;
        auto truth_it = truth_map.find(event_number);
        if (truth_it != truth_map.end()) {
            found_truth = true;
            true_px = truth_it->second.px; true_py = truth_it->second.py;
            true_pz = truth_it->second.pz;
            true_vx = truth_it->second.vx; true_vy = truth_it->second.vy;
            true_vz = truth_it->second.vz;
        }
        double true_p_mag = std::sqrt(
            true_px*true_px + true_py*true_py + true_pz*true_pz);
        TVector3 true_vertex_vec(true_vx, true_vy, true_vz);

        // Always record true momentum (even for unreconstructable events)
        all_true_p.push_back(true_p_mag);

        // --- Unpack hits ---
        double pip_pizza_time = -1, pim_pizza_time = -1;
        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;
        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;
        double pip_tof_time = -1, pim_tof_time = -1;
        double pip_tof_x = 0, pip_tof_y = 0, pip_tof_z = 0;
        double pim_tof_x = 0, pim_tof_y = 0, pim_tof_z = 0;
        int    pip_tof_devID = -1, pim_tof_devID = -1;
        std::vector<HitInfo> hits_pip, hits_pim;

        auto ev_it = event_data.find(event_number);
        if (ev_it != event_data.end()) {
            const auto& ev = ev_it->second;
            hits_pip = ev.hits_pip; hits_pim = ev.hits_pim;
            if (ev.has_pip_pizza) {
                pip_pizza_time = ev.pip_pizza_time;
                pip_pizza_x = ev.pip_pizza_x; pip_pizza_y = ev.pip_pizza_y;
                pip_pizza_z = ev.pip_pizza_z;
            }
            if (ev.has_pim_pizza) {
                pim_pizza_time = ev.pim_pizza_time;
                pim_pizza_x = ev.pim_pizza_x; pim_pizza_y = ev.pim_pizza_y;
                pim_pizza_z = ev.pim_pizza_z;
            }
            if (ev.has_pip_tof) {
                pip_tof_time  = ev.pip_tof_time;
                pip_tof_devID = ev.pip_tof_deviceID;
                pip_tof_x = ev.pip_tof_x; pip_tof_y = ev.pip_tof_y; pip_tof_z = ev.pip_tof_z;
            }
            if (ev.has_pim_tof) {
                pim_tof_time  = ev.pim_tof_time;
                pim_tof_devID = ev.pim_tof_deviceID;
                pim_tof_x = ev.pim_tof_x; pim_tof_y = ev.pim_tof_y; pim_tof_z = ev.pim_tof_z;
            }
        }

        bool reconstructable = false;
        double kaon_p = -1;
        TVector3 decay_vertex;

        // Attempt reconstruction — all failures fall through to reco_flag=0
        do {
            // PIZZA + TOF required for both pions
            if (pip_pizza_time < 0 || pip_tof_time < 0 ||
                pim_pizza_time < 0 || pim_tof_time < 0) break;

            // ----------------------------------------------------------------
            // FIT 3D TRACK LINES via TGraphErrors x(z) and y(z) for each pion.
            // See fit_3dgraph_track() above and KLong_save_vectors.C for details.
            // ----------------------------------------------------------------
            TrackFit fit_pip = fit_3dgraph_track(hits_pip, pip_tof_devID, z_tof);
            TrackFit fit_pim = fit_3dgraph_track(hits_pim, pim_tof_devID, z_tof);

            if (!fit_pip.valid || !fit_pim.valid) break;

            // PoCA -> decay vertex
            decay_vertex = closest_point_between_lines(
                fit_pip.orig, fit_pip.dir,
                fit_pim.orig, fit_pim.dir);

            // Near-parallel track cut: reject near-collinear pion pairs
            // (unstable PoCA denominator at small opening angles / high p_K)
            double cos_opening = fit_pip.dir.Dot(fit_pim.dir);
            if (cos_opening > 0.9998) break;   // cos(1.1 deg) ~ 0.9998

            // Track evaluation helper — evaluates fitted track at given z-plane
            auto eval_track_at_z = [](const TrackFit& f, double z) -> TVector3 {
                double t_param = (z - f.orig.Z()) / f.dir.Z();
                return f.orig + t_param * f.dir;
            };

            // PIZZA positions — track-extrapolated (replaces 5 cm Gaussian smear)
            TVector3 pip_pizza_pos = eval_track_at_z(fit_pip, pip_pizza_z);
            TVector3 pim_pizza_pos = eval_track_at_z(fit_pim, pim_pizza_z);

            // TOF positions — bar-ID smeared positions stored at hit-collection time.
            // x = Gaus(bar_centre_x, TOF_BAR_HALF_WIDTH), y = Gaus(bar_centre_y,
            // TOF_Y_UNCERTAINTY), z = bar face from geometry. Consistent with
            // the independently smeared time measurement.
            TVector3 pip_tof_pos(pip_tof_x, pip_tof_y, pip_tof_z);
            TVector3 pim_tof_pos(pim_tof_x, pim_tof_y, pim_tof_z);

            // Pion velocities — path = 3D distance between raw measured hit positions:
            //   pip_pizza_x/y/z = raw Geant4 hit on pizza (no fitting)
            //   pip_tof_x/y/z   = bar-ID-smeared TOF position (detector resolution)
            // Avoids two competing systematic biases:
            //  (1) Track-slope Mag(): Jensen's inequality (slope noise a_x^2/a_y^2 ≥ 0
            //      always) → systematic pip_track_cm↑ → pip_v↑ → p↓
            //  (2) z-projection |Δz|: underestimates by cos θ < 1 → pip_v↓ → p↑
            // Residual Jensen's bias from TOF bar-smearing ~σ_TOF²/(2L) ≈ 0.3%.
            TVector3 pip_pizza_raw(pip_pizza_x, pip_pizza_y, pip_pizza_z);
            TVector3 pim_pizza_raw(pim_pizza_x, pim_pizza_y, pim_pizza_z);
            // Jensen's correction: subtract known smearing variances from Mag2()
            const double sigma2_transverse =
                (TOF_BAR_HALF_WIDTH * TOF_BAR_HALF_WIDTH + smear_sigma * smear_sigma)
              + (TOF_Y_UNCERTAINTY  * TOF_Y_UNCERTAINTY  + smear_sigma * smear_sigma);
            const double dz_pip = pip_tof_z - pip_pizza_z;
            const double dz_pim = pim_tof_z - pim_pizza_z;
            double pip_track_cm = std::sqrt(std::max(
                (pip_tof_pos - pip_pizza_raw).Mag2() - sigma2_transverse,
                dz_pip * dz_pip));
            double pim_track_cm = std::sqrt(std::max(
                (pim_tof_pos - pim_pizza_raw).Mag2() - sigma2_transverse,
                dz_pim * dz_pim));
            double pip_dt_ns    = pip_tof_time - pip_pizza_time;
            double pim_dt_ns    = pim_tof_time - pim_pizza_time;
            double pip_v = (pip_dt_ns > 0)
                ? (pip_track_cm * 1e-2) / (pip_dt_ns * 1e-9) : 0.;
            double pim_v = (pim_dt_ns > 0)
                ? (pim_track_cm * 1e-2) / (pim_dt_ns * 1e-9) : 0.;

            // Kaon decay time
            double pip_path_cm = (pip_pizza_pos - decay_vertex).Mag();
            double pip_decay_t = pip_pizza_time * 1e-9 - (pip_path_cm * 1e-2 / pip_v);
            double pim_path_cm = (pim_pizza_pos - decay_vertex).Mag();
            double pim_decay_t = pim_pizza_time * 1e-9 - (pim_path_cm * 1e-2 / pim_v);
            double kaon_decay_time = 0.5 * (pip_decay_t + pim_decay_t);

            // Kaon momentum
            TVector3 kaon_prod(0., 0., 0.);
            double flight_cm = (decay_vertex - kaon_prod).Mag();
            double kaon_v    = (flight_cm * 1e-2) / kaon_decay_time;
            double m_K    = 0.497611;
            double beta_K = kaon_v / 2.99792458e8;
            if (beta_K >= 1.0) beta_K = 0.9999;
            double gamma_K = 1.0 / std::sqrt(1. - beta_K * beta_K);
            kaon_p = gamma_K * m_K * beta_K;

            if (kaon_p > 0 && kaon_p <= 11.) reconstructable = true;

        } while (false);  // Single-pass do-while used as a breakable block

        // Record per-event reconstruction flag
        all_reco_flags.push_back(reconstructable ? 1 : 0);

        if (reconstructable && found_truth && true_p_mag != 0.) {
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

    // ========================================================================
    // STEP 6: SAVE OUTPUT
    // ========================================================================
    std::string base = fname_str.substr(0, fname_str.find_last_of("."));
    std::string outFileName = base + "_acceptance.root";

    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    TTree *outTree = new TTree("kaonEventInfo", "Kaon event info for acceptance study");

    std::vector<double> true_mom_vec  = all_true_p;
    std::vector<int>    reco_flag_vec = all_reco_flags;
    int n_triple_pion_events = (int)selected_events.size();

    outTree->Branch("n_triple_pion_events", &n_triple_pion_events);
    outTree->Branch("true_mom_vec",         &true_mom_vec);
    outTree->Branch("reco_flag_vec",        &reco_flag_vec);

    outTree->Fill();
    outFile->Write();
    outFile->Close();
    file->Close();

    std::cout << "Saved kaon event info to " << outFileName << "\n";
}

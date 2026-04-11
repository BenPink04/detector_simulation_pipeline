// ============================================================================
// KLong_save_momentum_acceptance.C — MODULE-BASED RECONSTRUCTION (Option A)
//
// PURPOSE:
//   Same reconstruction as KLong_save_vectors.C, but records ALL events that
//   have the correct decay channel (pi+ pi- pi0), not just reconstructed ones.
//   A per-event reco_flag (1=reconstructed, 0=not) is saved alongside the
//   true momentum, enabling acceptance/efficiency calculations.
//
//   Uses geometry-aware module-based hit reconstruction (Option A overhaul).
//   See KLong_save_vectors.C and TRACKING_OVERHAUL_README.md for full details.
//
// INPUT:  a simulation .root file containing Ntuple1, Ntuple2, Ntuple3
//         Filename must contain detector position tags, e.g.:
//           T1-240_T2-250_T3-0_T4-0_P1-215_P2-230_F1-260_F2-270_E1-700
//
// OUTPUT: <input_base>_acceptance.root  (one TTree named "kaonEventInfo")
//         Branches:
//           n_triple_pion_events  : total number of events with pi+pi-pi0 decay
//           true_mom_vec          : vector of true kaon momenta (GeV/c)
//           reco_flag_vec         : vector of ints (1=reco success, 0=failure)
// ============================================================================

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
    double pip_tof_z    = 0;  // raw Geant4 hit z (cm)
    int    pip_tof_deviceID = -1;

    bool   has_pim_pizza = false;
    double pim_pizza_time = 0;
    double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;

    bool   has_pim_tof = false;
    double pim_tof_time = 0;
    double pim_tof_z    = 0;  // raw Geant4 hit z (cm)
    int    pim_tof_deviceID = -1;
};

struct ModuleMeasurement {
    double x_centre, y_centre, z_pos;
    double x_half_range, y_half_range;
    bool   valid;
};

struct TrackPoint {
    double x, y, z;   // centre of the 2D overlap area (cm)
    double sx, sy;    // half-widths of the overlap area (cm)
};

// 1D interval [centre - half_width, centre + half_width].
// valid=false means the interval is empty (disjoint constraints).
struct Interval1D {
    double centre;
    double half_width;
    bool   valid;
};

// Intersect a set of 1D constraint regions supplied as {centre, half_width} pairs.
// Returns the tightest interval consistent with ALL inputs.
// Returns {0,0,false} if any two regions are disjoint.
Interval1D intersect_constraints(
    const std::vector<std::pair<double,double>>& regions)
{
    if (regions.empty()) return {0., 9999., false};
    double lo = regions[0].first - regions[0].second;
    double hi = regions[0].first + regions[0].second;
    for (size_t k = 1; k < regions.size(); ++k) {
        double new_lo = regions[k].first - regions[k].second;
        double new_hi = regions[k].first + regions[k].second;
        lo = std::max(lo, new_lo);
        hi = std::min(hi, new_hi);
        if (lo > hi) return {0., 0., false};  // disjoint — reject
    }
    return {0.5*(lo + hi), 0.5*(hi - lo), true};
}

// ============================================================================
// GEOMETRY HELPER FUNCTIONS
// ============================================================================

// Globally unique IDs encode both station and straw (station * 488 + straw_cn):
//   Station 0 (IDs    1-488) : straw_cn = ID,        sub-layer = (straw_cn-1)/122
//   Station 1 (IDs  489-976) : straw_cn = ID - 488,  sub-layer = (straw_cn-1)/122
//   Station 2 (IDs 977-1464) : straw_cn = ID - 976,  sub-layer = (straw_cn-1)/122
//   Station 3 (IDs 1465-1952): straw_cn = ID - 1464, sub-layer = (straw_cn-1)/122
double tracker_local_x(int deviceID) {
    int straw_cn  = ((deviceID - 1) % 488) + 1;  // 1-488 within station box
    int sub_layer = (straw_cn - 1) / 122;        // 0, 1, 2, or 3
    int i         = (straw_cn - 1) % 122;        // straw index 0..121 within sub-layer
    double x = -48.4 + i * 0.8;                 // cm, base position
    // Sub-layers 1 and 3 are staggered +0.4 cm (half of 0.8 cm pitch)
    if (sub_layer == 1 || sub_layer == 3) x += 0.4;
    return x;
}

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

ModuleMeasurement get_module_measurement(int deviceID,
                                          const double z_tracker[4],
                                          const double z_fri[2],
                                          double z_tof)
{
    ModuleMeasurement m;
    m.x_centre = 0; m.y_centre = 0; m.z_pos = 0;
    m.x_half_range = 9999; m.y_half_range = 9999;
    m.valid = false;

    if (deviceID >= 1 && deviceID <= 1952)
    {
        int station = (deviceID - 1) / 488;
        double lx   = tracker_local_x(deviceID);
        m.z_pos = z_tracker[station];
        if (m.z_pos <= 0.) return m;
        m.valid = true;
        if (station == 0) {
            m.x_centre = lx; m.y_centre = 0.;
            m.x_half_range = STRAW_HALF_WIDTH; m.y_half_range = STRAW_HALF_LENGTH;
        }
        else if (station == 1) {
            m.x_centre = 0.; m.y_centre = lx;
            m.x_half_range = STRAW_HALF_LENGTH; m.y_half_range = STRAW_HALF_WIDTH;
        }
        else if (station == 2) {
            m.x_centre = lx / SQRT2; m.y_centre = lx / SQRT2;
            m.x_half_range = STRAW_HALF_LENGTH / SQRT2;
            m.y_half_range = STRAW_HALF_LENGTH / SQRT2;
        }
        else {
            m.x_centre =  lx / SQRT2; m.y_centre = -lx / SQRT2;
            m.x_half_range = STRAW_HALF_LENGTH / SQRT2;
            m.y_half_range = STRAW_HALF_LENGTH / SQRT2;
        }
    }
    else if (deviceID >= 2001 && deviceID <= 2052)
    {
        int fri_index = deviceID - 2001;
        int wall      = fri_index / 26;
        int strip_i   = fri_index % 26;
        m.z_pos = z_fri[wall];
        if (m.z_pos <= 0.) return m;
        m.valid = true;
        double lx = FRI_STRIP_X[strip_i];
        double hw = fri_half_width(strip_i);
        double hl = FRI_HALF_LENGTH[strip_i];
        if (wall == 0) {
            m.x_centre = lx; m.y_centre = 0.;
            m.x_half_range = hw; m.y_half_range = hl;
        }
        else {
            m.x_centre = 0.; m.y_centre = lx;
            m.x_half_range = hl; m.y_half_range = hw;
        }
    }
    else if (deviceID >= 2053 && deviceID <= 2070)
    {
        m.z_pos = z_tof;
        if (m.z_pos <= 0.) return m;
        m.valid = true;
        m.x_centre     = tof_bar_x_centre(deviceID);
        m.y_centre     = tof_bar_y_centre(deviceID);
        m.x_half_range = TOF_BAR_HALF_WIDTH;
        m.y_half_range = TOF_Y_UNCERTAINTY;
    }

    return m;
}

// ============================================================================
// TRACK FITTING FUNCTIONS
// ============================================================================

TVector3 fit_track_direction(const std::vector<TrackPoint>& pts)
{
    if (pts.size() < 2) return TVector3(0., 0., 1.);

    double sw=0, swz=0, swz2=0;
    double swx=0, swzx=0, swy=0, swzy=0;

    for (const auto& p : pts) {
        double wx = 1.0 / (p.sx * p.sx + 1e-6);
        double wy = 1.0 / (p.sy * p.sy + 1e-6);
        double w  = 0.5 * (wx + wy);
        sw   += w;  swz  += w*p.z;  swz2 += w*p.z*p.z;
        swx  += wx*p.x;  swzx += wx*p.z*p.x;
        swy  += wy*p.y;  swzy += wy*p.z*p.y;
    }

    double denom = sw * swz2 - swz * swz;
    if (std::abs(denom) < 1e-12) return TVector3(0., 0., 1.);

    double ax = (sw * swzx - swz * swx) / denom;
    double ay = (sw * swzy - swz * swy) / denom;
    return TVector3(ax, ay, 1.0).Unit();
}

TVector3 fit_track_origin(const std::vector<TrackPoint>& pts)
{
    if (pts.empty()) return TVector3(0., 0., 0.);
    double mx = 0, my = 0, mz = 0;
    for (const auto& p : pts) { mx += p.x; my += p.y; mz += p.z; }
    int n = (int)pts.size();
    return TVector3(mx/n, my/n, mz/n);
}

TVector3 closest_point_between_lines(const TVector3& p1, const TVector3& v1,
                                      const TVector3& p2, const TVector3& v2)
{
    TVector3 w0 = p1 - p2;
    double a = v1.Dot(v1), b = v1.Dot(v2), c = v2.Dot(v2);
    double d = v1.Dot(w0), e = v2.Dot(w0);
    double denom = a*c - b*b;
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
    double z_tracker[4] = {
        extract_param(fname_str, "T1"), extract_param(fname_str, "T2"),
        extract_param(fname_str, "T3"), extract_param(fname_str, "T4")
    };
    double z_fri[2] = {
        extract_param(fname_str, "F1"), extract_param(fname_str, "F2")
    };
    double z_tof = extract_param(fname_str, "E1");

    std::cout << "Detector Z positions:\n"
              << "  Tracker: T1=" << z_tracker[0] << " T2=" << z_tracker[1]
              << " T3=" << z_tracker[2] << " T4=" << z_tracker[3] << " cm\n"
              << "  FRI: F1=" << z_fri[0] << " F2=" << z_fri[1] << " cm\n"
              << "  TOF: E1=" << z_tof << " cm\n";

    // ========================================================================
    // STEP 1: SELECT EVENTS WITH K_L -> pi+ pi- pi0
    // [TO CHANGE] For a different decay channel change the PDG codes below.
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
    // smear_sigma applies to PIZZA only; tracker/FRI/TOF use geometry lookup.
    // [TO CHANGE] Tune smear_time_sigma for your timing resolution.
    // -----------------------------------------------------------------------
    TRandom3 randGen(0);
    double smear_sigma      = 5.0;    // cm  — PIZZA position smear only
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
                    ev.pip_pizza_x = x; ev.pip_pizza_y = y; ev.pip_pizza_z = z;
                }
            }
            if (is_tof(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pip_tof || st < ev.pip_tof_time) {
                    ev.has_pip_tof = true; ev.pip_tof_time = st;
                    ev.pip_tof_z = z; // raw Geant4 hit z (not geometry face z_tof)
                    ev.pip_tof_deviceID = id;
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
                    ev.pim_pizza_x = x; ev.pim_pizza_y = y; ev.pim_pizza_z = z;
                }
            }
            if (is_tof(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pim_tof || st < ev.pim_tof_time) {
                    ev.has_pim_tof = true; ev.pim_tof_time = st;
                    ev.pim_tof_z = z; // raw Geant4 hit z (not geometry face z_tof)
                    ev.pim_tof_deviceID = id;
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

        // Always record true momentum (even if not reconstructable)
        all_true_p.push_back(true_p_mag);

        // --- Unpack hits ---
        double pip_pizza_time = -1, pim_pizza_time = -1;
        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;
        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;
        double pip_tof_time = -1, pim_tof_time = -1;
        double pip_tof_z_raw = z_tof, pim_tof_z_raw = z_tof; // default to geometry face
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
                pip_tof_z_raw = ev.pip_tof_z;
                pip_tof_devID = ev.pip_tof_deviceID;
            }
            if (ev.has_pim_tof) {
                pim_tof_time  = ev.pim_tof_time;
                pim_tof_z_raw = ev.pim_tof_z;
                pim_tof_devID = ev.pim_tof_deviceID;
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
            // Track reconstruction: DSSSD-style pixel building.
            // Mirrors KLong_save_vectors.C — see that file and README for details.
            // Sub-layer 0 only; T1+T2 → DSSSD pixel; F1+F2 → FRI pixel;
            // T3+T4 → stereo pixel via u/v unfolding.
            // Arrival-angle correction applied using truth momentum direction.
            // ----------------------------------------------------------------
            double truth_dxdz = (std::abs(true_pz) > 1e-6) ? true_px / true_pz : 0.;
            double truth_dydz = (std::abs(true_pz) > 1e-6) ? true_py / true_pz : 0.;

            auto build_track_points = [&](const std::vector<HitInfo>& hits,
                                          double td_xdz, double td_ydz)
                -> std::vector<TrackPoint>
            {
                // Sub-layer 0 filter
                std::vector<HitInfo> filtered;
                for (const auto& h : hits) {
                    if (h.deviceID >= 1 && h.deviceID <= 1952) {
                        int straw_cn = ((h.deviceID - 1) % 488) + 1;
                        if ((straw_cn - 1) / 122 != 0) continue;
                    }
                    filtered.push_back(h);
                }

                std::map<int, std::vector<std::pair<double,double>>> x_regs, y_regs;
                std::map<int, double> x_zpos, y_zpos;
                std::vector<std::pair<double,double>> t3_u_regs, t4_v_regs;

                for (const auto& h : filtered) {
                    if (h.deviceID >= 977 && h.deviceID <= 1098 && z_tracker[2] > 0) {
                        t3_u_regs.push_back({tracker_local_x(h.deviceID), STRAW_HALF_WIDTH});
                        continue;
                    }
                    if (h.deviceID >= 1465 && h.deviceID <= 1586 && z_tracker[3] > 0) {
                        t4_v_regs.push_back({tracker_local_x(h.deviceID), STRAW_HALF_WIDTH});
                        continue;
                    }
                    ModuleMeasurement m = get_module_measurement(
                        h.deviceID, z_tracker, z_fri, z_tof);
                    if (!m.valid) continue;

                    int key;
                    if      (h.deviceID >= 1    && h.deviceID <= 488)  key = 0;
                    else if (h.deviceID >= 489  && h.deviceID <= 976)  key = 1;
                    else if (h.deviceID >= 2001 && h.deviceID <= 2026) key = 10;
                    else if (h.deviceID >= 2027 && h.deviceID <= 2052) key = 11;
                    else continue;

                    bool is_xmeas = (m.x_half_range < m.y_half_range);
                    if (is_xmeas) {
                        x_regs[key].push_back({m.x_centre, m.x_half_range});
                        x_zpos[key] = m.z_pos;
                    } else {
                        y_regs[key].push_back({m.y_centre, m.y_half_range});
                        y_zpos[key] = m.z_pos;
                    }
                }

                std::vector<TrackPoint> tps;
                bool has_x_c = false, has_y_c = false;

                auto add_dsssd_pixel = [&](int xk, int yk, double free_half) {
                    bool hx = x_regs.count(xk) && !x_regs[xk].empty();
                    bool hy = y_regs.count(yk) && !y_regs[yk].empty();
                    if (hx && hy) {
                        Interval1D xi = intersect_constraints(x_regs[xk]);
                        Interval1D yi = intersect_constraints(y_regs[yk]);
                        if (xi.valid && yi.valid) {
                            double zx = x_zpos[xk], zy = y_zpos[yk];
                            double zm = 0.5 * (zx + zy);
                            double xm = xi.centre + td_xdz * (zm - zx);
                            double ym = yi.centre + td_ydz * (zm - zy);
                            tps.push_back({xm, ym, zm, xi.half_width, yi.half_width});
                            has_x_c = true; has_y_c = true;
                        }
                    } else if (hx) {
                        Interval1D xi = intersect_constraints(x_regs[xk]);
                        if (xi.valid) {
                            tps.push_back({xi.centre, 0., x_zpos[xk],
                                           xi.half_width, free_half});
                            has_x_c = true;
                        }
                    } else if (hy) {
                        Interval1D yi = intersect_constraints(y_regs[yk]);
                        if (yi.valid) {
                            tps.push_back({0., yi.centre, y_zpos[yk],
                                           free_half, yi.half_width});
                            has_y_c = true;
                        }
                    }
                };

                add_dsssd_pixel(0, 1, STRAW_HALF_LENGTH);  // T1(X)+T2(Y)
                add_dsssd_pixel(10, 11, 70.25);             // F1(X)+F2(Y)

                // T3+T4 stereo: x=(u+v)/√2, y=(u-v)/√2
                if (!t3_u_regs.empty() && !t4_v_regs.empty()) {
                    Interval1D ui = intersect_constraints(t3_u_regs);
                    Interval1D vi = intersect_constraints(t4_v_regs);
                    if (ui.valid && vi.valid) {
                        double zm  = 0.5 * (z_tracker[2] + z_tracker[3]);
                        double xs  = (ui.centre + vi.centre) / SQRT2;
                        double ys  = (ui.centre - vi.centre) / SQRT2;
                        double sxy = std::hypot(ui.half_width, vi.half_width) / SQRT2;
                        tps.push_back({xs, ys, zm, sxy, sxy});
                        has_x_c = true; has_y_c = true;
                    }
                }

                if (!has_x_c || !has_y_c) return {};
                return tps;
            };

            std::vector<TrackPoint> track_points_pip = build_track_points(hits_pip, truth_dxdz, truth_dydz);
            std::vector<TrackPoint> track_points_pim = build_track_points(hits_pim, truth_dxdz, truth_dydz);

            if (track_points_pip.size() < 2 || track_points_pim.size() < 2) break;

            // Weighted line-of-best-fit track directions
            TVector3 pip_dir   = fit_track_direction(track_points_pip);
            TVector3 pip_start = fit_track_origin(track_points_pip);
            TVector3 pim_dir   = fit_track_direction(track_points_pim);
            TVector3 pim_start = fit_track_origin(track_points_pim);

            // PoCA -> decay vertex
            decay_vertex = closest_point_between_lines(
                pip_start, pip_dir, pim_start, pim_dir);

            // PIZZA positions (smeared)
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

            // TOF positions from bar geometry
            TVector3 pip_tof_pos(
                tof_bar_x_centre(pip_tof_devID),
                tof_bar_y_centre(pip_tof_devID) + randGen.Gaus(0., TOF_Y_UNCERTAINTY),
                z_tof
            );
            TVector3 pim_tof_pos(
                tof_bar_x_centre(pim_tof_devID),
                tof_bar_y_centre(pim_tof_devID) + randGen.Gaus(0., TOF_Y_UNCERTAINTY),
                z_tof
            );

            // Pion velocities — path uses track-fit extrapolation to z_tof,
            // not the smeared ToF hit position.  TOF_Y_UNCERTAINTY = 10 cm
            // always inflates the 3D path length in quadrature, systematically
            // overestimating pip_v and therefore underestimating kaon momentum.
            double t_pip_tof = (pip_dir.Z() != 0.) ? (pip_tof_z_raw - pip_start.Z()) / pip_dir.Z() : 0.;
            TVector3 pip_tof_fit = pip_start + t_pip_tof * pip_dir;
            double t_pim_tof = (pim_dir.Z() != 0.) ? (pim_tof_z_raw - pim_start.Z()) / pim_dir.Z() : 0.;
            TVector3 pim_tof_fit = pim_start + t_pim_tof * pim_dir;
            double pip_track_cm = (pip_tof_fit - pip_pizza_pos).Mag();
            double pim_track_cm = (pim_tof_fit - pim_pizza_pos).Mag();
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

        } while (false);  // Single-pass do-while used as a break-able block

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

    std::vector<double> true_mom_vec = all_true_p;
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

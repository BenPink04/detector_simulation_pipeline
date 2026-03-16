// ============================================================================
// KLong_save_vectors.C  —  3D TGRAPH TRACK RECONSTRUCTION
//
// PURPOSE:
//   Reconstructs K_L (K-long) momentum and decay vertex from simulated
//   detector hits, for events where the kaon decays via K_L -> pi+ pi- pi0.
//
//   Track directions are found by fitting straight lines to the raw Geant4
//   3D hit positions using TGraphErrors.  For each pion, two independent
//   TGraphErrors fits are built:
//     - x(z) : using X-measuring detectors (T1, FRI-Wall1, TOF bar centre)
//     - y(z) : using Y-measuring detectors (T2, FRI-Wall2, TOF bar centre)
//   Each fit is a pol1 (straight line).  The resulting slopes and intercepts
//   define the 3D track direction for each pion.
//
//   T3/T4 stereo hits contribute to both x(z) and y(z) fits with a larger
//   positional uncertainty (STRAW_HALF_WIDTH * sqrt(2)) to account for the
//   +/-45 degree rotation of the stereo straws.
//
//   All four sub-layers of each tracker station contribute to the fit.
//   The TOF bar hit is added to both fits to ensure at least two distinct
//   z-planes are available for the line fit.
//
//   For each reconstructable event the macro saves:
//     - reco_p      : reconstructed kaon momentum magnitude (GeV/c)
//     - true_p      : true (MC truth) kaon momentum magnitude (GeV/c)
//     - reco_vertex : reconstructed decay vertex (x, y, z) in cm
//     - true_vertex : true (MC truth) decay vertex (x, y, z) in cm
//
// INPUT:  a simulation .root file containing Ntuple1, Ntuple2, Ntuple3
//         Filename must contain detector position tags, e.g.:
//           T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700
//
// OUTPUT: <input_base>_vectors.root  (one TTree named "kaonVectors")
//
// RECONSTRUCTION STRATEGY (summary):
//   1. Select events with K_L -> pi+ pi- pi0 (from Ntuple2)
//   2. Load MC truth decay momentum/vertex (from Ntuple1)
//   3. Load detector hits for pi+ and pi- (from Ntuple3)
//   4. For each pion: build TGraphErrors x(z) and y(z) from raw Geant4 hit
//      positions, weighted by per-detector measurement precision
//   5. Fit pol1 to each graph -> track direction and reference origin
//   6. PoCA of the two track lines -> reconstructed decay vertex
//   7. Pion velocities from PIZZA (TOF start) and TOF wall (TOF stop)
//   8. Back-propagate to decay vertex -> kaon decay time -> kaon momentum
//
// SEE ALSO:
//   KLong_save_vectors_pixelation.C  — alternative DSSSD/pixel intersection approach
//   TRACKING_OVERHAUL_README.md      — full geometry and reconstruction details
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
#include <iomanip>

// ============================================================================
// GEOMETRY CONSTANTS  (taken directly from JLabKDetectorConstruction.cc)
// ============================================================================

// --- Tracker straws ---
const double STRAW_HALF_WIDTH  = 0.40;   // cm  — measurement precision transverse to tube
const double STRAW_HALF_LENGTH = 45.0;   // cm  — unconstrained axis along tube

// --- FRI strips ---
// 26 strips per wall.  X centre of each strip in wall's local frame.
const double FRI_STRIP_X[26] = {
    -66., -60., -54., -48., -42., -36., -30., -24., -18., -12., -7.5,
     -4.5, -1.5,  1.5,  4.5,  7.5,  12.,  18.,  24.,  30.,  36.,  42.,
     48.,  54.,  60.,  66.
};

// Half-length of each FRI strip along its long axis (varies near beam hole).
const double FRI_HALF_LENGTH[26] = {
    37.25, 42.75, 53.75, 59.25, 59.25, 64.75,
    70.25, 70.25, 70.25, 70.25,
    70.25, 70.25, 70.25, 70.25, 70.25, 70.25,
    70.25, 70.25, 70.25, 70.25,
    64.75, 59.25, 59.25, 53.75, 42.75, 37.25
};

// --- TOF wall ---
const double TOF_BAR_HALF_WIDTH = 6.0;   // cm  (X, bar half-width)
const double TOF_Y_UNCERTAINTY  = 10.0;  // cm  (Y, from light-path timing along bar)

const double SQRT2 = 1.41421356237;

// ============================================================================
// STRUCTS
// ============================================================================

// One raw detector hit.  x/y/z are the raw Geant4 positions in the lab frame.
struct HitInfo { double x, y, z, t; int deviceID; };

struct TruthInfo {
    double px, py, pz;   // GeV/c
    double vx, vy, vz;   // cm
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

    bool   has_pim_pizza = false;
    double pim_pizza_time = 0;
    double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;

    bool   has_pim_tof = false;
    double pim_tof_time = 0;
    int    pim_tof_deviceID = -1;
};

// Result of a 3D TGraph straight-line track fit.
struct TrackFit {
    TVector3 dir;    // unit direction vector  (dx/dz, dy/dz, 1).Unit()
    TVector3 orig;   // reference point on track, evaluated at z = 0:
                     //   orig = (bx, by, 0)  where x(z) = ax*z + bx
    bool     valid;  // false if fit failed or insufficient z-span
};

// ============================================================================
// GEOMETRY HELPER FUNCTIONS
// ============================================================================

// Measurement precision (half-width) of a FRI strip. strip_i is 0-based within a wall.
double fri_half_width(int strip_i) {
    return (strip_i >= 10 && strip_i <= 15) ? 1.5 : 3.0;
}

// X centre of a TOF bar in lab frame. copyNumber range 2053..2070.
double tof_bar_x_centre(int copyNumber) {
    int i = copyNumber - 2053;      // 0..17
    return -96.0 + i * 12.0;        // cm
}

// Nominal Y centre of a TOF bar (cm).
double tof_bar_y_centre(int copyNumber) {
    int i = copyNumber - 2053;
    if (i == 8) return -72.0;
    if (i == 9) return +72.0;
    return 0.0;
}

// ============================================================================
// 3D TGRAPH TRACK FIT
//
// Builds two TGraphErrors objects from the raw Geant4 hit positions:
//   gx : z-axis = hit z,  y-axis = hit x,  y-errors = detector x-precision
//   gy : z-axis = hit z,  y-axis = hit y,  y-errors = detector y-precision
//
// Each is fitted with pol1:  f(z) = p0 + p1*z
//   => track: x(z) = bx + ax*z,  y(z) = by + ay*z
//   => direction vector: (ax, ay, 1).Unit()
//   => reference point at z=0: (bx, by, 0)
//
// DETECTOR CONTRIBUTIONS:
//   T1  (IDs   1-488)   : measures X; contributes to gx only
//                         sigma_x = STRAW_HALF_WIDTH (0.4 cm)
//   T2  (IDs 489-976)   : measures Y; contributes to gy only
//                         sigma_y = STRAW_HALF_WIDTH (0.4 cm)
//   T3  (IDs 977-1464)  : stereo +45deg; contributes to both gx and gy
//                         sigma = STRAW_HALF_WIDTH * sqrt(2)  (propagated)
//   T4  (IDs 1465-1952) : stereo -45deg; same as T3
//   FRI-W1 (2001-2026)  : measures X; contributes to gx only
//                         sigma_x = fri_half_width(strip_i)  (1.5 or 3.0 cm)
//   FRI-W2 (2027-2052)  : measures Y; contributes to gy only
//                         sigma_y = fri_half_width(strip_i)
//   TOF    (2053-2070)  : added via tof_devID; bar centre X/Y from geometry
//                         sigma_x = TOF_BAR_HALF_WIDTH (6 cm)
//                         sigma_y = TOF_Y_UNCERTAINTY  (10 cm)
//
// All four sub-layers of each tracker station contribute hit points to the fit.
// The fit requires >= 2 points AND >= 2 distinct z-planes in each of gx and gy.
//
// Parameters:
//   hits     : raw tracker + FRI HitInfo objects (h.x/h.y/h.z from Geant4)
//   tof_devID: deviceID of the TOF bar hit (-1 if none)
//   z_tof    : z position of TOF wall parsed from filename (cm)
// ============================================================================
TrackFit fit_3dgraph_track(const std::vector<HitInfo>& hits,
                           int    tof_devID,
                           double z_tof)
{
    const TrackFit FAIL = {{0,0,1}, {0,0,0}, false};

    // Accumulate (z_pos, measurement, sigma) for x-fit and y-fit separately
    std::vector<double> xz, xv, xe;   // for x(z) = ax*z + bx
    std::vector<double> yz, yv, ye;   // for y(z) = ay*z + by

    for (const auto& h : hits) {
        int id = h.deviceID;

        if (id >= 1 && id <= 488) {
            // T1 (station 0): straws oriented along Y, measures X
            // Raw Geant4 h.x is the true particle x position in the straw.
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(STRAW_HALF_WIDTH);
        }
        else if (id >= 489 && id <= 976) {
            // T2 (station 1): straws oriented along X, measures Y
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(STRAW_HALF_WIDTH);
        }
        else if (id >= 977 && id <= 1464) {
            // T3 (station 2): stereo +45 deg.
            // Straw oriented along (1/sqrt2, 1/sqrt2, 0).
            // Measurement precision STRAW_HALF_WIDTH is along the perpendicular
            // (rotated) axis; projecting onto lab X and Y gives sigma * sqrt(2).
            double sxy = STRAW_HALF_WIDTH * SQRT2;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(sxy);
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(sxy);
        }
        else if (id >= 1465 && id <= 1952) {
            // T4 (station 3): stereo -45 deg, same error propagation as T3
            double sxy = STRAW_HALF_WIDTH * SQRT2;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(sxy);
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(sxy);
        }
        else if (id >= 2001 && id <= 2026) {
            // FRI Wall 1: strips oriented along Y, measures X
            int strip_i = id - 2001;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(fri_half_width(strip_i));
        }
        else if (id >= 2027 && id <= 2052) {
            // FRI Wall 2: strips oriented along X, measures Y
            int strip_i = id - 2027;
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(fri_half_width(strip_i));
        }
        // TOF hits are added below via tof_devID (bar positions from geometry)
    }

    // Add TOF bar hit: positions from geometry lookup (not raw Geant4 position)
    if (tof_devID >= 2053 && tof_devID <= 2070 && z_tof > 0) {
        xz.push_back(z_tof);
        xv.push_back(tof_bar_x_centre(tof_devID));
        xe.push_back(TOF_BAR_HALF_WIDTH);

        yz.push_back(z_tof);
        yv.push_back(tof_bar_y_centre(tof_devID));
        ye.push_back(TOF_Y_UNCERTAINTY);
    }

    // Require >= 2 fit points in each projection
    if ((int)xz.size() < 2 || (int)yz.size() < 2) return FAIL;

    // Require >= 2 distinct z-planes for a meaningful slope fit
    // (all points at the same z -> degenerate: slope undefined)
    auto has_distinct_z = [](const std::vector<double>& zv) -> bool {
        for (size_t k = 1; k < zv.size(); ++k)
            if (std::fabs(zv[k] - zv[0]) > 1.0) return true;
        return false;
    };
    if (!has_distinct_z(xz) || !has_distinct_z(yz)) return FAIL;

    // ----------------------------------------------------------------
    // Fit x(z) with pol1:  x = p0 + p1*z
    //   p0 = bx (intercept at z=0), p1 = ax (slope dx/dz)
    // ----------------------------------------------------------------
    TGraphErrors gx((int)xz.size(), xz.data(), xv.data(), nullptr, xe.data());
    gx.Fit("pol1", "Q");   // "Q" = quiet (no printout)
    TF1* fx = gx.GetFunction("pol1");
    if (!fx) return FAIL;
    double bx = fx->GetParameter(0);
    double ax = fx->GetParameter(1);

    // ----------------------------------------------------------------
    // Fit y(z) with pol1:  y = p0 + p1*z
    // ----------------------------------------------------------------
    TGraphErrors gy((int)yz.size(), yz.data(), yv.data(), nullptr, ye.data());
    gy.Fit("pol1", "Q");
    TF1* fy = gy.GetFunction("pol1");
    if (!fy) return FAIL;
    double by = fy->GetParameter(0);
    double ay = fy->GetParameter(1);

    // Sanity: reject non-finite parameters (fit diverged)
    if (!std::isfinite(ax) || !std::isfinite(ay) ||
        !std::isfinite(bx) || !std::isfinite(by)) return FAIL;

    // Direction vector: (ax, ay, 1).Unit()
    TVector3 dir(ax, ay, 1.0);
    dir = dir.Unit();
    // Reference point: track evaluated at z = 0
    TVector3 orig(bx, by, 0.);
    return {dir, orig, true};
}

// ============================================================================
// POINT OF CLOSEST APPROACH (PoCA) between two straight 3D lines.
// Line 1: p1 + s*v1   Line 2: p2 + t*v2
// Returns the midpoint of the two closest points -> approximated decay vertex.
// ============================================================================
TVector3 closest_point_between_lines(const TVector3& p1, const TVector3& v1,
                                      const TVector3& p2, const TVector3& v2)
{
    TVector3 w0 = p1 - p2;
    double a = v1.Dot(v1), b = v1.Dot(v2), c = v2.Dot(v2);
    double d = v1.Dot(w0), e = v2.Dot(w0);
    double denom = a*c - b*b;
    if (std::abs(denom) < 1e-14) return 0.5 * (p1 + p2); // parallel lines
    double sc = (b*e - c*d) / denom;
    double tc = (a*e - b*d) / denom;
    return 0.5 * ((p1 + v1*sc) + (p2 + v2*tc));
}

// ===========================================================================
// MAIN FUNCTION
// [TO CHANGE] Update the default filename for a different scenario or seed.
// ===========================================================================
void KLong_save_vectors(const char* filename = "Scenario3_Seed1.root") {

    // [TO CHANGE] Thread count for implicit MT.
    ROOT::EnableImplicitMT(4);
    auto start_time = std::chrono::steady_clock::now();

    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) { std::cout << "Cannot open root file\n"; return; }

    // ========================================================================
    // PARSE DETECTOR Z POSITIONS FROM FILENAME
    //
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

    // Suppress unused variable warning (z_tracker, z_fri used only as reference;
    // actual z values for the track fit come from raw Geant4 h.z)
    (void)z_tracker; (void)z_fri;

    std::cout << "Detector Z positions (from filename):\n"
              << "  T1=" << extract_param(fname_str,"T1")
              << "  T2=" << extract_param(fname_str,"T2")
              << "  T3=" << extract_param(fname_str,"T3")
              << "  T4=" << extract_param(fname_str,"T4") << " cm\n"
              << "  F1=" << extract_param(fname_str,"F1")
              << "  F2=" << extract_param(fname_str,"F2") << " cm\n"
              << "  E1=" << z_tof << " cm\n";

    // ========================================================================
    // STEP 1: SELECT EVENTS WITH K_L -> pi+ pi- pi0
    //
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
    // [TO CHANGE] Update ranges if detector geometry changes.
    // -----------------------------------------------------------------------
    auto is_tracker  = [](int id) { return (id >=    1 && id <=  1952); };
    auto is_fri      = [](int id) { return (id >= 2001 && id <=  2052); };
    auto is_pizza    = [](int id) { return (id >= 1953 && id <=  2000); };
    auto is_tof      = [](int id) { return (id >= 2053 && id <=  2070); };

    // -----------------------------------------------------------------------
    // RESOLUTION PARAMETERS
    // [TO CHANGE] Tune smear_time_sigma for your timing resolution.
    // [TO CHANGE] Set smear_sigma=0 for ideal PIZZA position.
    // -----------------------------------------------------------------------
    TRandom3 randGen(0);
    double smear_sigma      = 5.0;    // cm  — PIZZA position smear (TOF start)
    double smear_time_sigma = 0.0015; // ns  — 1.5 ps timing smear

    // [TO CHANGE] Set DEBUG_MODE = false to suppress per-event diagnostics.
    const bool DEBUG_MODE       = true;
    const int  DEBUG_MAX_EVENTS = 20;
    int        debug_event_count = 0;

    // Output storage
    std::vector<double> reco_p, true_p;
    std::vector<double> reco_vertex_x, reco_vertex_y, reco_vertex_z;
    std::vector<double> true_vertex_x, true_vertex_y, true_vertex_z;

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
    //
    // Tracker + FRI : store raw HitInfo (raw Geant4 x/y/z used in track fit).
    // PIZZA         : keep raw x/y/z + earliest smeared time (TOF start).
    // TOF           : keep deviceID + earliest smeared time (TOF stop).
    //                 Bar centre position is derived from deviceID in fit.
    // ========================================================================
    std::unordered_map<int, EventReco> event_data;
    Long64_t nEntries3 = tree3->GetEntries();
    for (Long64_t i = 0; i < nEntries3; ++i) {
        tree3->GetEntry(i);
        int evt_id = (int)evtNb;
        if (selected_event_set.find(evt_id) == selected_event_set.end()) continue;

        auto &ev = event_data[evt_id];
        int id = (int)devID;

        // --- pi+ (PDG 211) ---
        if ((int)pdg == 211) {
            if (is_tracker(id) || is_fri(id)) {
                ev.hits_pip.push_back({x, y, z, t, id});
            }
            if (is_pizza(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pip_pizza || st < ev.pip_pizza_time) {
                    ev.has_pip_pizza  = true;
                    ev.pip_pizza_time = st;
                    ev.pip_pizza_x = x; ev.pip_pizza_y = y; ev.pip_pizza_z = z;
                }
            }
            if (is_tof(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pip_tof || st < ev.pip_tof_time) {
                    ev.has_pip_tof      = true;
                    ev.pip_tof_time     = st;
                    ev.pip_tof_deviceID = id;
                }
            }
        }

        // --- pi- (PDG -211) ---
        if ((int)pdg == -211) {
            if (is_tracker(id) || is_fri(id)) {
                ev.hits_pim.push_back({x, y, z, t, id});
            }
            if (is_pizza(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pim_pizza || st < ev.pim_pizza_time) {
                    ev.has_pim_pizza  = true;
                    ev.pim_pizza_time = st;
                    ev.pim_pizza_x = x; ev.pim_pizza_y = y; ev.pim_pizza_z = z;
                }
            }
            if (is_tof(id)) {
                double st = randGen.Gaus(t, smear_time_sigma);
                if (!ev.has_pim_tof || st < ev.pim_tof_time) {
                    ev.has_pim_tof      = true;
                    ev.pim_tof_time     = st;
                    ev.pim_tof_deviceID = id;
                }
            }
        }
    }

    // ========================================================================
    // STEP 5: RECONSTRUCT KAON MOMENTUM — EVENT LOOP
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

        // ----------------------------------------------------------------
        // Unpack pre-built event data
        // ----------------------------------------------------------------
        double pip_pizza_time = -1, pim_pizza_time = -1;
        double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;
        double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;
        double pip_tof_time = -1, pim_tof_time = -1;
        int    pip_tof_devID = -1, pim_tof_devID = -1;
        std::vector<HitInfo> hits_pip, hits_pim;

        auto ev_it = event_data.find(event_number);
        if (ev_it != event_data.end()) {
            const auto& ev = ev_it->second;
            hits_pip = ev.hits_pip;
            hits_pim = ev.hits_pim;
            if (ev.has_pip_pizza) {
                pip_pizza_time = ev.pip_pizza_time;
                pip_pizza_x = ev.pip_pizza_x;
                pip_pizza_y = ev.pip_pizza_y;
                pip_pizza_z = ev.pip_pizza_z;
            }
            if (ev.has_pim_pizza) {
                pim_pizza_time = ev.pim_pizza_time;
                pim_pizza_x = ev.pim_pizza_x;
                pim_pizza_y = ev.pim_pizza_y;
                pim_pizza_z = ev.pim_pizza_z;
            }
            if (ev.has_pip_tof) {
                pip_tof_time  = ev.pip_tof_time;
                pip_tof_devID = ev.pip_tof_deviceID;
            }
            if (ev.has_pim_tof) {
                pim_tof_time  = ev.pim_tof_time;
                pim_tof_devID = ev.pim_tof_deviceID;
            }
        }

        // ----------------------------------------------------------------
        // Device ID to detector name (for debug output).
        // ----------------------------------------------------------------
        auto device_id_info = [](int id) -> std::string {
            if (id >=    1 && id <=  488) return "T1 (tracker st.0, measures X)";
            if (id >=  489 && id <=  976) return "T2 (tracker st.1, measures Y)";
            if (id >=  977 && id <= 1464) return "T3 (tracker st.2, stereo U +45deg)";
            if (id >= 1465 && id <= 1952) return "T4 (tracker st.3, stereo V -45deg)";
            if (id >= 1953 && id <= 1976) return "PIZZA-flat (TOF start)";
            if (id >= 1977 && id <= 2000) return "PIZZA-conical (TOF start)";
            if (id >= 2001 && id <= 2026) return "FRI-Wall1/F1 (measures X)";
            if (id >= 2027 && id <= 2052) return "FRI-Wall2/F2 (measures Y)";
            if (id >= 2053 && id <= 2070) return "TOF-wall bar (X=centre, Y~10cm)";
            return "UNKNOWN";
        };

        // ----------------------------------------------------------------
        // DEBUG: per-event hit summary (first DEBUG_MAX_EVENTS events only)
        // ----------------------------------------------------------------
        if (DEBUG_MODE && debug_event_count < DEBUG_MAX_EVENTS) {
            ++debug_event_count;
            std::cout << std::fixed << std::setprecision(1);
            std::cout << "\n=== DEBUG EVENT " << event_number
                      << " (" << debug_event_count << "/" << DEBUG_MAX_EVENTS
                      << ") ===  [pi+pi-pi0 decay]\n";
            std::cout << "  Reconstruction (3D TGraph method):\n"
                         "    - PIZZA + TOF hit required for each pion\n"
                         "    - TGraphErrors x(z) and y(z) fitted with pol1\n"
                         "    - Need >= 2 fit points at >= 2 distinct z-planes per axis\n"
                         "    - TOF bar hit always added to both fits\n"
                         "    - All tracker sub-layers used\n";

            for (int pion = 0; pion < 2; ++pion) {
                const std::vector<HitInfo>& phits = (pion==0) ? hits_pip : hits_pim;
                bool has_pz = (pion==0) ? (pip_pizza_time>=0) : (pim_pizza_time>=0);
                bool has_tf = (pion==0) ? (pip_tof_time  >=0) : (pim_tof_time  >=0);
                int  tf_dev = (pion==0) ? pip_tof_devID       : pim_tof_devID;
                const char* pn = (pion==0) ? "pi+" : "pi-";

                int cnt_t1=0, cnt_t2=0, cnt_t3=0, cnt_t4=0, cnt_f1=0, cnt_f2=0;
                std::set<int> uid_set;
                for (const auto& h : phits) {
                    int id = h.deviceID;
                    uid_set.insert(id);
                    if      (id <=  488) cnt_t1++;
                    else if (id <=  976) cnt_t2++;
                    else if (id <= 1464) cnt_t3++;
                    else if (id <= 1952) cnt_t4++;
                    else if (id <= 2026) cnt_f1++;
                    else if (id <= 2052) cnt_f2++;
                }
                std::cout << "  " << pn << ": "
                          << phits.size() << " tracker/FRI hits"
                          << "  PIZZA=" << (has_pz ? "YES" : "NO ")
                          << "  TOF=" << (has_tf ? "YES" : "NO") << "\n";
                std::cout << "    Hit counts:"
                          << "  T1=" << cnt_t1 << "  T2=" << cnt_t2
                          << "  T3=" << cnt_t3 << "  T4=" << cnt_t4
                          << "  F1=" << cnt_f1 << "  F2=" << cnt_f2 << "\n";
                if (!uid_set.empty()) {
                    std::cout << "    Fired device IDs (" << uid_set.size() << "):\n";
                    for (int uid : uid_set)
                        std::cout << "      ID=" << std::setw(4) << uid
                                  << "  ->  " << device_id_info(uid) << "\n";
                }
                // Preview TGraph fit
                TrackFit dbg_fit = fit_3dgraph_track(phits, tf_dev, z_tof);
                if (dbg_fit.valid)
                    std::cout << "    TGraph fit OK:"
                              << "  dir=(" << dbg_fit.dir.X() << ", "
                                           << dbg_fit.dir.Y() << ", "
                                           << dbg_fit.dir.Z() << ")"
                              << "  orig_xy=(" << dbg_fit.orig.X() << ", "
                                               << dbg_fit.orig.Y() << ")\n";
                else
                    std::cout << "    TGraph fit FAILED\n";
                if (!has_pz)
                    std::cout << "    [FAIL-TOF: no PIZZA hit for " << pn << "]\n";
                if (!has_tf)
                    std::cout << "    [FAIL-TOF: no TOF wall hit for " << pn << "]\n";
                if (!dbg_fit.valid) {
                    std::cout << "    [FAIL-TRACKING: need >= 2 x-measurements and"
                                 " >= 2 y-measurements at distinct z-planes]\n";
                }
            }
        }

        // ----------------------------------------------------------------
        // RECONSTRUCTION CUT: require PIZZA + TOF hits for both pions
        // ----------------------------------------------------------------
        if (pip_pizza_time < 0 || pip_tof_time < 0 ||
            pim_pizza_time < 0 || pim_tof_time < 0) continue;

        // ----------------------------------------------------------------
        // FIT 3D TRACK LINES: TGraphErrors x(z) and y(z), pol1 for each pion.
        // The TOF bar hit is included in each fit to extend the z-span.
        // ----------------------------------------------------------------
        TrackFit fit_pip = fit_3dgraph_track(hits_pip, pip_tof_devID, z_tof);
        TrackFit fit_pim = fit_3dgraph_track(hits_pim, pim_tof_devID, z_tof);

        if (!fit_pip.valid || !fit_pim.valid) continue;

        // ----------------------------------------------------------------
        // DECAY VERTEX: Point of Closest Approach between the two track lines
        // ----------------------------------------------------------------
        TVector3 decay_vertex = closest_point_between_lines(
            fit_pip.orig, fit_pip.dir,
            fit_pim.orig, fit_pim.dir);

        // ----------------------------------------------------------------
        // PIZZA hit positions — smeared (TOF start position).
        // [TO CHANGE] Set smear_sigma=0 for ideal reconstruction.
        // ----------------------------------------------------------------
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

        // ----------------------------------------------------------------
        // TOF wall hit positions — derived from bar geometry (not raw smear).
        // X: bar centre exact.  Y: nominal bar Y + Gaussian smear.
        // [TO CHANGE] Switch to Uniform for flat Y prior if preferred.
        // ----------------------------------------------------------------
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

        // ----------------------------------------------------------------
        // PION VELOCITY  v = path_length / delta_t
        // [TO CHANGE] Replace Mag() with arc length for curved trajectories.
        // ----------------------------------------------------------------
        double pip_track_cm = (pip_tof_pos - pip_pizza_pos).Mag();
        double pim_track_cm = (pim_tof_pos - pim_pizza_pos).Mag();
        double pip_dt_ns = pip_tof_time - pip_pizza_time;
        double pim_dt_ns = pim_tof_time - pim_pizza_time;
        double pip_v = (pip_dt_ns > 0)
            ? (pip_track_cm * 1e-2) / (pip_dt_ns * 1e-9) : 0.;
        double pim_v = (pim_dt_ns > 0)
            ? (pim_track_cm * 1e-2) / (pim_dt_ns * 1e-9) : 0.;

        // ----------------------------------------------------------------
        // KAON DECAY TIME — back-propagate each pion from PIZZA to vertex.
        // [TO CHANGE] Use weighted average if measurement qualities differ.
        // ----------------------------------------------------------------
        double pip_path_cm    = (pip_pizza_pos - decay_vertex).Mag();
        double pip_decay_time = pip_pizza_time * 1e-9
                                - (pip_path_cm * 1e-2 / pip_v);

        double pim_path_cm    = (pim_pizza_pos - decay_vertex).Mag();
        double pim_decay_time = pim_pizza_time * 1e-9
                                - (pim_path_cm * 1e-2 / pim_v);

        double kaon_decay_time = 0.5 * (pip_decay_time + pim_decay_time); // s

        // ----------------------------------------------------------------
        // KAON MOMENTUM  p = gamma * m_K * beta
        //
        // [TO CHANGE] Update kaon_prod if kaon not produced at origin.
        // [TO CHANGE] m_K = K_L mass. Use 0.493677 GeV/c^2 for K+/K-.
        // ----------------------------------------------------------------
        TVector3 kaon_prod(0., 0., 0.);
        double flight_cm = (decay_vertex - kaon_prod).Mag();
        double kaon_v    = (flight_cm * 1e-2) / kaon_decay_time;

        double m_K    = 0.497611; // GeV/c^2 — K_L mass
        double beta_K = kaon_v / 2.99792458e8;
        if (beta_K >= 1.0) beta_K = 0.9999;
        double gamma_K = 1.0 / std::sqrt(1. - beta_K * beta_K);
        double kaon_p  = gamma_K * m_K * beta_K;

        // ----------------------------------------------------------------
        // MC TRUTH LOOKUP
        // ----------------------------------------------------------------
        auto truth_it = truth_map.find(event_number);
        if (truth_it == truth_map.end()) continue;

        double true_px = truth_it->second.px;
        double true_py = truth_it->second.py;
        double true_pz = truth_it->second.pz;
        double true_p_mag = std::sqrt(
            true_px*true_px + true_py*true_py + true_pz*true_pz);
        if (true_p_mag == 0.) continue;

        TVector3 true_vertex_vec(
            truth_it->second.vx,
            truth_it->second.vy,
            truth_it->second.vz);

        // ----------------------------------------------------------------
        // RECONSTRUCTION CUT: reject unphysically large reco momenta.
        // [TO CHANGE] Adjust upper bound to match beam momentum range.
        // ----------------------------------------------------------------
        if (kaon_p > 11.) continue;

        std::cout << "Event " << event_number
                  << " | Reco p: " << kaon_p
                  << " | True p: " << true_p_mag
                  << " | Reco vertex: (" << decay_vertex.X()      << ", "
                                          << decay_vertex.Y()     << ", "
                                          << decay_vertex.Z()     << ")"
                  << " | True vertex: (" << true_vertex_vec.X()  << ", "
                                          << true_vertex_vec.Y() << ", "
                                          << true_vertex_vec.Z() << ")\n";

        reco_p.push_back(kaon_p);
        true_p.push_back(true_p_mag);
        reco_vertex_x.push_back(decay_vertex.X());
        reco_vertex_y.push_back(decay_vertex.Y());
        reco_vertex_z.push_back(decay_vertex.Z());
        true_vertex_x.push_back(true_vertex_vec.X());
        true_vertex_y.push_back(true_vertex_vec.Y());
        true_vertex_z.push_back(true_vertex_vec.Z());

    } // end event loop

    // ========================================================================
    // STEP 6: SAVE OUTPUT TO ROOT FILE
    //
    // [TO CHANGE] Add outTree->Branch() calls to save additional quantities.
    // [TO CHANGE] Rename "_vectors" suffix if saving a different set of info.
    // ========================================================================
    std::string base = fname_str.substr(0, fname_str.find_last_of("."));
    std::string outFileName = base + "_vectors.root";

    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    TTree *outTree = new TTree("kaonVectors", "Kaon momentum and vertex vectors");

    std::vector<double> *p_reco   = &reco_p;
    std::vector<double> *p_true   = &true_p;
    std::vector<double> *v_reco_x = &reco_vertex_x;
    std::vector<double> *v_reco_y = &reco_vertex_y;
    std::vector<double> *v_reco_z = &reco_vertex_z;
    std::vector<double> *v_true_x = &true_vertex_x;
    std::vector<double> *v_true_y = &true_vertex_y;
    std::vector<double> *v_true_z = &true_vertex_z;

    outTree->Branch("reco_p",        &p_reco);
    outTree->Branch("true_p",        &p_true);
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

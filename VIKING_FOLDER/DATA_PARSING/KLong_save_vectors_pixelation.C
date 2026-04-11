// ============================================================================
// KLong_save_vectors.C  —  MODULE-BASED RECONSTRUCTION (Option A overhaul)
//
// PURPOSE:
//   Reconstructs K_L (K-long) momentum and decay vertex from simulated
//   detector hits, for events where the kaon decays via K_L -> pi+ pi- pi0.
//
//   Hit positions are derived from the known geometry of the module that
//   fired (tracker straw / FRI strip / TOF bar) rather than from the raw
//   Geant4 hit coordinates.  The raw hit position in the ROOT ntuple is
//   intentionally ignored for tracking; only deviceID, hitT, and PDG are
//   used from Ntuple3.
//
//   For each reconstructable event it saves:
//     - reco_p      : reconstructed kaon momentum magnitude (GeV/c)
//     - true_p      : true (MC truth) kaon momentum magnitude (GeV/c)
//     - reco_vertex : reconstructed decay vertex (x, y, z) in cm
//     - true_vertex : true (MC truth) decay vertex (x, y, z) in cm
//
// INPUT:  a simulation .root file containing Ntuple1, Ntuple2, Ntuple3
//         The filename must contain the detector position information in the
//         standard format, e.g.:
//           T1-240_T2-250_T3-0_T4-0_P1-215_P2-230_F1-260_F2-270_E1-700
//         These values are parsed to set the Z positions of each sub-detector.
//
// OUTPUT: <input_base>_vectors.root  (one TTree named "kaonVectors")
//
// RECONSTRUCTION STRATEGY (summary):
//   1. Select events with K_L -> pi+ pi- pi0 decay products (from Ntuple2)
//   2. Load MC truth decay momentum/vertex for each event (from Ntuple1)
//   3. Load all detector hits for pi+ and pi- (from Ntuple3, using deviceID)
//   4. For each tracker/FRI hit, look up the module centre position and
//      which coordinate axis it constrains (X or Y)
//   5. Require at least one X-constraining AND one Y-constraining hit per pion
//      (two-detector coincidence) to form a 2D track measurement
//   6. Build TrackPoints (best-estimate 3D positions) from X/Y hit pairs at
//      the same detector plane (within Z_PAIR_TOLERANCE)
//   7. Fit a weighted straight line through all TrackPoints -> track direction
//   8. PoCA between pi+ and pi- track lines -> reconstructed decay vertex
//   9. Pion velocities from PIZZA (TOF start) and TOF wall (TOF stop)
//  10. Back-propagate to decay vertex -> kaon decay time -> kaon momentum
//
// GEOMETRY HANDLED (device ID ranges from JLabKDetectorConstruction.cc):
//   Tracker station 0  (IDs    1-488) : measures X  (phi=0 deg,   T1)
//   Tracker station 1  (IDs  489-976) : measures Y  (phi=-90 deg, T2)
//   Tracker station 2  (IDs 977-1464) : stereo U +45 deg (T3, add to fit later)
//   Tracker station 3  (IDs 1465-1952): stereo V -45 deg (T4, add to fit later)
//   PIZZA flat         (IDs 1953-1976): TOF start only, not used for tracking
//   PIZZA conical      (IDs 1977-2000): TOF start only, not used for tracking
//   FRI wall 1         (IDs 2001-2026): measures X  (phi=0 deg)
//   FRI wall 2         (IDs 2027-2052): measures Y  (phi=+90 deg)
//   TOF wall           (IDs 2053-2070): X from bar centre, Y +/- TOF_Y_UNCERTAINTY
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
#include <iomanip>

// ============================================================================
// GEOMETRY CONSTANTS  (taken directly from JLabKDetectorConstruction.cc)
// ============================================================================

// --- Tracker straws ---
// 122 straws per station; outer radius 0.40 cm; half-length 45 cm along tube.
const double STRAW_HALF_WIDTH  = 0.40;   // cm (measurement precision, transverse to tube)
const double STRAW_HALF_LENGTH = 45.0;   // cm (unconstrained axis, along tube)

// --- FRI strips ---
// 26 strips per wall.  Outer (wide) strips 6 cm; inner strips 3 cm wide.
// X centre of each strip in the wall's local frame (symmetric about 0).
const double FRI_STRIP_X[26] = {
    -66., -60., -54., -48., -42., -36., -30., -24., -18., -12., -7.5,
     -4.5, -1.5,  1.5,  4.5,  7.5,  12.,  18.,  24.,  30.,  36.,  42.,
     48.,  54.,  60.,  66.
};

// Half-length of each FRI strip along its long axis (varies due to beam hole cutout).
// Indexed 0..25 within a single wall.
const double FRI_HALF_LENGTH[26] = {
    37.25, 42.75, 53.75, 59.25, 59.25, 64.75,   // i = 0..5  (edge, short strips)
    70.25, 70.25, 70.25, 70.25,                  // i = 6..9  (full-length wide)
    70.25, 70.25, 70.25, 70.25, 70.25, 70.25,    // i = 10..15 (inner 3 cm strips)
    70.25, 70.25, 70.25, 70.25,                  // i = 16..19 (full-length wide)
    64.75, 59.25, 59.25, 53.75, 42.75, 37.25     // i = 20..25 (edge, short strips)
};

// --- TOF wall ---
// 18 bars, 12 cm pitch.  Bar half-width 6 cm.
// Y uncertainty modelled as Gaussian with sigma = TOF_Y_UNCERTAINTY (10 cm).
const double TOF_BAR_HALF_WIDTH = 6.0;   // cm (X, bar half-width)
const double TOF_Y_UNCERTAINTY  = 10.0;  // cm (Y, from light-path timing along bar)

const double SQRT2 = 1.41421356237;

// ============================================================================
// STRUCTS
// ============================================================================

// One raw detector hit stored during the Ntuple3 scan.
// x/y/z are retained but are NOT used for tracker/FRI/TOF reconstruction
// (geometry lookup replaces them). They are used only for PIZZA position.
struct HitInfo { double x, y, z, t; int deviceID; };

// MC-truth info for a kaon decay.
struct TruthInfo {
    double px, py, pz;   // GeV/c
    double vx, vy, vz;   // cm
};

// Per-event collected hits and timing.
struct EventReco {
    std::vector<HitInfo> hits_pip;  // tracker + FRI hits for pi+
    std::vector<HitInfo> hits_pim;  // tracker + FRI hits for pi-

    // PIZZA: TOF start — raw position smeared at reconstruction time
    bool   has_pip_pizza = false;
    double pip_pizza_time = 0;      // ns (smeared)
    double pip_pizza_x = 0, pip_pizza_y = 0, pip_pizza_z = 0;

    // TOF wall: TOF stop — only deviceID stored; position from geometry lookup
    bool   has_pip_tof = false;
    double pip_tof_time = 0;        // ns (smeared)
    double pip_tof_z    = 0;        // raw Geant4 hit z (cm)
    int    pip_tof_deviceID = -1;

    bool   has_pim_pizza = false;
    double pim_pizza_time = 0;
    double pim_pizza_x = 0, pim_pizza_y = 0, pim_pizza_z = 0;

    bool   has_pim_tof = false;
    double pim_tof_time = 0;
    double pim_tof_z    = 0;        // raw Geant4 hit z (cm)
    int    pim_tof_deviceID = -1;
};

// Result of a module geometry lookup.
struct ModuleMeasurement {
    double x_centre;      // best-estimate X in lab frame (cm)
    double y_centre;      // best-estimate Y in lab frame (cm)
    double z_pos;         // Z of the detector plane (cm)
    double x_half_range;  // +-uncertainty in X; large value = unconstrained axis
    double y_half_range;  // +-uncertainty in Y; large value = unconstrained axis
    bool   valid;         // false if deviceID unrecognised or detector disabled (z=0)
};

// A 2D+Z position built from an intra-module intersection of X and Y constraints.
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

// Local-frame X centre of a tracker straw (cm).
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

// Half-width (measurement axis) of a FRI strip (cm).
// strip_i is the 0-based index within a single wall (0..25).
double fri_half_width(int strip_i) {
    return (strip_i >= 10 && strip_i <= 15) ? 1.5 : 3.0;
}

// X centre of a TOF bar in lab frame (cm). copyNumber range: 2053..2070.
double tof_bar_x_centre(int copyNumber) {
    int i = copyNumber - 2053;         // 0..17
    return -96.0 + i * 12.0;           // cm
}

// Nominal Y centre of a TOF bar (cm).
// Bars i=8 (ID 2061) and i=9 (ID 2062) straddle the beam gap.
double tof_bar_y_centre(int copyNumber) {
    int i = copyNumber - 2053;
    if (i == 8) return -72.0;
    if (i == 9) return +72.0;
    return 0.0;
}

// ============================================================================
// MODULE MEASUREMENT LOOKUP
//
// For a given deviceID and detector Z positions (parsed from the filename),
// returns the best-estimate 2D position and per-axis uncertainties.
//
// z_tracker[4] = {T1, T2, T3, T4} in cm — 0 means that station is disabled.
// z_fri[2]     = {F1, F2}          in cm — 0 means that wall is disabled.
// z_tof        = E1                in cm — 0 means TOF wall is disabled.
// ============================================================================
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
        // --- TRACKER STRAW ---
        int station = (deviceID - 1) / 488;       // 0=X(T1), 1=Y(T2), 2=U(T3), 3=V(T4)
        double lx   = tracker_local_x(deviceID);  // local X before rotation

        m.z_pos = z_tracker[station];
        if (m.z_pos <= 0.) return m;   // station disabled
        m.valid = true;

        if (station == 0) {
            // phi=0 deg: tube runs along Y, measures X
            m.x_centre     = lx;
            m.y_centre     = 0.;
            m.x_half_range = STRAW_HALF_WIDTH;
            m.y_half_range = STRAW_HALF_LENGTH;
        }
        else if (station == 1) {
            // phi=-90 deg: box rotated, local X maps to lab Y
            m.x_centre     = 0.;
            m.y_centre     = lx;
            m.x_half_range = STRAW_HALF_LENGTH;
            m.y_half_range = STRAW_HALF_WIDTH;
        }
        else if (station == 2) {
            // phi=+45 deg: stereo U plane, constrains x+y = sqrt(2)*lx
            // [TO CHANGE] Use proper stereo unfolding when U/V stations are enabled.
            m.x_centre     = lx / SQRT2;
            m.y_centre     = lx / SQRT2;
            m.x_half_range = STRAW_HALF_LENGTH / SQRT2;
            m.y_half_range = STRAW_HALF_LENGTH / SQRT2;
        }
        else {
            // phi=-45 deg: stereo V plane
            m.x_centre     =  lx / SQRT2;
            m.y_centre     = -lx / SQRT2;
            m.x_half_range = STRAW_HALF_LENGTH / SQRT2;
            m.y_half_range = STRAW_HALF_LENGTH / SQRT2;
        }
    }
    else if (deviceID >= 2001 && deviceID <= 2052)
    {
        // --- FRI STRIP ---
        int fri_index = deviceID - 2001;         // 0..51 across both walls
        int wall      = fri_index / 26;          // 0=Wall1(measures X), 1=Wall2(measures Y)
        int strip_i   = fri_index % 26;          // 0..25 within one wall

        m.z_pos = z_fri[wall];
        if (m.z_pos <= 0.) return m;   // wall disabled
        m.valid = true;

        double lx = FRI_STRIP_X[strip_i];
        double hw = fri_half_width(strip_i);
        double hl = FRI_HALF_LENGTH[strip_i];

        if (wall == 0) {
            // Wall 1: strips run along Y, measures X
            m.x_centre     = lx;
            m.y_centre     = 0.;
            m.x_half_range = hw;
            m.y_half_range = hl;
        }
        else {
            // Wall 2: rotated 90 deg, strips run along X, measures Y
            m.x_centre     = 0.;
            m.y_centre     = lx;
            m.x_half_range = hl;
            m.y_half_range = hw;
        }
    }
    else if (deviceID >= 2053 && deviceID <= 2070)
    {
        // --- TOF BAR ---
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
// WEIGHTED LINEAR REGRESSION — TRACK DIRECTION FIT
//
// Fits x = ax*z + bx and y = ay*z + by independently through a set of
// TrackPoints, weighting each point by 1/sigma^2.
// Returns the unit direction vector (ax, ay, 1).Unit().
// Falls back to (0,0,1) if fewer than 2 points.
// ============================================================================
TVector3 fit_track_direction(const std::vector<TrackPoint>& pts)
{
    if (pts.size() < 2) return TVector3(0., 0., 1.);

    double sw=0, swz=0, swz2=0;
    double swx=0, swzx=0;
    double swy=0, swzy=0;

    for (const auto& p : pts)
    {
        double wx = 1.0 / (p.sx * p.sx + 1e-6);
        double wy = 1.0 / (p.sy * p.sy + 1e-6);
        double w  = 0.5 * (wx + wy);

        sw   += w;
        swz  += w  * p.z;
        swz2 += w  * p.z * p.z;
        swx  += wx * p.x;
        swzx += wx * p.z * p.x;
        swy  += wy * p.y;
        swzy += wy * p.z * p.y;
    }

    double denom = sw * swz2 - swz * swz;
    if (std::abs(denom) < 1e-12) return TVector3(0., 0., 1.);

    double ax = (sw * swzx - swz * swx) / denom;  // dx/dz
    double ay = (sw * swzy - swz * swy) / denom;  // dy/dz

    return TVector3(ax, ay, 1.0).Unit();
}

// Return the centroid of the TrackPoints as a reference point on the fitted line.
TVector3 fit_track_origin(const std::vector<TrackPoint>& pts)
{
    if (pts.empty()) return TVector3(0., 0., 0.);
    double mx = 0, my = 0, mz = 0;
    for (const auto& p : pts) { mx += p.x; my += p.y; mz += p.z; }
    int n = (int)pts.size();
    return TVector3(mx/n, my/n, mz/n);
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
    return 0.5 * (point1 + point2);
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
    // The filename must contain the standard position tags, e.g.:
    //   T1-240_T2-250_T3-0_T4-0_P1-215_P2-230_F1-260_F2-270_E1-700
    // A value of 0 means that detector plane is disabled.
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
        extract_param(fname_str, "T1"),
        extract_param(fname_str, "T2"),
        extract_param(fname_str, "T3"),
        extract_param(fname_str, "T4")
    };
    double z_fri[2] = {
        extract_param(fname_str, "F1"),
        extract_param(fname_str, "F2")
    };
    double z_tof = extract_param(fname_str, "E1");

    std::cout << "Detector Z positions parsed from filename:\n"
              << "  Tracker : T1=" << z_tracker[0] << "  T2=" << z_tracker[1]
              << "  T3=" << z_tracker[2] << "  T4=" << z_tracker[3] << " cm\n"
              << "  FRI     : F1=" << z_fri[0]  << "  F2=" << z_fri[1]  << " cm\n"
              << "  TOF     : E1=" << z_tof << " cm\n";

    // ========================================================================
    // STEP 1: SELECT EVENTS WITH K_L -> pi+ pi- pi0
    //
    // [TO CHANGE] For a different decay channel, change the PDG codes and
    //             has_pip/has_pim/has_pi0 booleans below.
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
            if (code ==  211) has_pip = true;  // pi+  [TO CHANGE for other channels]
            if (code == -211) has_pim = true;  // pi-
            if (code ==  111) has_pi0 = true;  // pi0
        }
        if (has_pip && has_pim && has_pi0)    // [TO CHANGE] condition for desired mode
            selected_events.push_back(kv.first);
    }

    size_t total_events = selected_events.size();
    std::cout << "Found " << total_events
              << " events with pi+, pi-, pi0 as direct kaon decay products.\n"
              << "Processing " << total_events << " selected events...\n";

    std::unordered_set<int> selected_event_set(
        selected_events.begin(), selected_events.end());

    // ========================================================================
    // STEP 2: SET UP Ntuple3  (raw detector hits)
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
    //
    // [TO CHANGE] If detector geometry changes, update these ID ranges to
    //             match the new copy number assignments in
    //             JLabKDetectorConstruction.cc.
    // -----------------------------------------------------------------------
    auto is_tracker    = [](int id) { return (id >=    1 && id <=  1952); };
    auto is_fri        = [](int id) { return (id >= 2001 && id <=  2052); };
    auto is_pizza      = [](int id) { return (id >= 1953 && id <=  2000); };
    auto is_tof        = [](int id) { return (id >= 2053 && id <=  2070); };

    // Sub-classify FRI walls and tracker stations for X vs Y measurement:
    auto is_fri_wall1  = [](int id) { return (id >= 2001 && id <= 2026); }; // measures X
    auto is_fri_wall2  = [](int id) { return (id >= 2027 && id <= 2052); }; // measures Y

    // Returns station index: 0=X(T1), 1=Y(T2), 2=U(T3), 3=V(T4)
    auto tracker_station = [](int id) -> int {
        if (id < 1 || id > 1952) return -1;
        return (id - 1) / 488;
    };

    // -----------------------------------------------------------------------
    // RESOLUTION PARAMETERS
    //
    // smear_sigma is used ONLY for PIZZA hit positions (TOF start).
    // Tracker, FRI and TOF wall positions come from geometry lookup, not smear.
    //
    // [TO CHANGE] Tune smear_time_sigma for the timing resolution of your study.
    // [TO CHANGE] Set smear_sigma=0 for ideal PIZZA position.
    // -----------------------------------------------------------------------
    TRandom3 randGen(0);
    double smear_sigma      = 5.0;    // cm  — PIZZA position smear
    double smear_time_sigma = 0.0015; // ns  — timing smear (1.5 ps) all detectors

    // [TO CHANGE] Set DEBUG_MODE = false to suppress per-event hit diagnostics.
    // [TO CHANGE] Increase DEBUG_MAX_EVENTS to see more event summaries.
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
    // Tracker + FRI hits  : store raw HitInfo (deviceID + t used; x/y/z ignored
    //                       in reconstruction — geometry lookup replaces them).
    // PIZZA hits          : keep raw x/y/z + earliest smeared time (TOF start).
    // TOF hits            : keep deviceID + earliest smeared time (TOF stop).
    //                       Bar centre position is looked up from deviceID later.
    //
    // [TO CHANGE] If you want a different hit selection strategy (e.g. use
    //             the hit closest to the extrapolated track rather than
    //             earliest in time), modify the is/has logic below.
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
                    ev.pip_tof_z        = z; // raw Geant4 hit z (not geometry face z_tof)
                    ev.pip_tof_deviceID = id;  // bar centre looked up in event loop
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
                    ev.pim_tof_z        = z; // raw Geant4 hit z (not geometry face z_tof)
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
        double pip_tof_z_raw = z_tof, pim_tof_z_raw = z_tof; // default to geometry face
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
                pip_tof_z_raw = ev.pip_tof_z;
                pip_tof_devID = ev.pip_tof_deviceID;
            }
            if (ev.has_pim_tof) {
                pim_tof_time  = ev.pim_tof_time;
                pim_tof_z_raw = ev.pim_tof_z;
                pim_tof_devID = ev.pim_tof_deviceID;
            }
        }

        // ----------------------------------------------------------------
        // Device ID to module name (for debug output).
        // ----------------------------------------------------------------
        auto device_id_info = [](int id) -> std::string {
            if (id >=    1 && id <=  488) return "T1 (tracker st.0, measures X)";
            if (id >=  489 && id <=  976) return "T2 (tracker st.1, measures Y)";
            if (id >=  977 && id <= 1464) return "T3 (tracker st.2, stereo U +45deg)";
            if (id >= 1465 && id <= 1952) return "T4 (tracker st.3, stereo V -45deg)";
            if (id >= 1953 && id <= 1976) return "PIZZA-flat (TOF start, not tracking)";
            if (id >= 1977 && id <= 2000) return "PIZZA-conical (TOF start, not tracking)";
            if (id >= 2001 && id <= 2026) return "FRI-Wall1/F1 (measures X)";
            if (id >= 2027 && id <= 2052) return "FRI-Wall2/F2 (measures Y)";
            if (id >= 2053 && id <= 2070) return "TOF-wall bar (X=centre, Y~10cm)";
            return "UNKNOWN";
        };

        // ----------------------------------------------------------------
        // Track reconstruction: build DSSSD-style TrackPoints from hits.
        //
        // SIMPLIFIED MODEL (per supervisor guidance, March 2026):
        //   - Tracker: use sub-layer 0 only (straws 1-122 per station).
        //     This models one layer of X-tubes (T1) and one of Y-tubes (T2),
        //     matching realistic electronics readout constraints.
        //   - DSSSD assumption: pair T1(X)+T2(Y) into one tight 2D pixel at
        //     their mid-Z.  Similarly pair F1(X)+F2(Y) and T3(U)+T4(V).
        //   - Arrival angle correction: the particle is not travelling purely
        //     along Z.  When projecting the X measurement at z=T1 and the Y
        //     measurement at z=T2 to the common mid-Z, we correct for the
        //     transverse slope using the MC-truth momentum direction (dx/dz,
        //     dy/dz).  This is an approved approximation for the simplified
        //     tracker model (supervisor email, March 2026).
        //   - T3+T4: proper stereo unfolding (u+v → x,y) replaces the
        //     approximated equal-half-range treatment.
        //
        // truth_dxdz, truth_dydz: true px/pz, py/pz for the kaon event
        //   (passed in; set to 0 if truth is unavailable).
        // ----------------------------------------------------------------
        auto build_track_points = [&](const std::vector<HitInfo>& hits,
                                      double truth_dxdz, double truth_dydz)
            -> std::vector<TrackPoint>
        {
            // --- A: Filter tracker to sub-layer 0 only ---
            // Sub-layer 0: straw_cn 1-122 within each 488-ID station block.
            std::vector<HitInfo> filtered;
            for (const auto& h : hits) {
                if (h.deviceID >= 1 && h.deviceID <= 1952) {
                    int straw_cn = ((h.deviceID - 1) % 488) + 1;
                    if ((straw_cn - 1) / 122 != 0) continue; // skip sub-layers 1,2,3
                }
                filtered.push_back(h);
            }

            // --- B: Separate T3/T4 stereo hits; accumulate X/Y for T1,T2,F1,F2 ---
            std::map<int, std::vector<std::pair<double,double>>> x_regs, y_regs;
            std::map<int, double> x_zpos, y_zpos;
            // T3 u-measurements (sub-layer-0 IDs 977-1098) and
            // T4 v-measurements (sub-layer-0 IDs 1465-1586)
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
                if      (h.deviceID >= 1    && h.deviceID <= 488)  key = 0;  // T1 X
                else if (h.deviceID >= 489  && h.deviceID <= 976)  key = 1;  // T2 Y
                else if (h.deviceID >= 2001 && h.deviceID <= 2026) key = 10; // F1 X
                else if (h.deviceID >= 2027 && h.deviceID <= 2052) key = 11; // F2 Y
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

            // --- C: Build DSSSD pixels for each station pair ---
            std::vector<TrackPoint> tps;
            bool has_x_c = false, has_y_c = false;

            // Helper: form one pixel from an X station (xk) and Y station (yk).
            // Applies arrival-angle correction to project both axes to mid-Z.
            auto add_dsssd_pixel = [&](int xk, int yk, double free_half) {
                bool hx = x_regs.count(xk) && !x_regs[xk].empty();
                bool hy = y_regs.count(yk) && !y_regs[yk].empty();
                if (hx && hy) {
                    Interval1D xi = intersect_constraints(x_regs[xk]);
                    Interval1D yi = intersect_constraints(y_regs[yk]);
                    if (xi.valid && yi.valid) {
                        double zx = x_zpos[xk], zy = y_zpos[yk];
                        double zm = 0.5 * (zx + zy);
                        // Project X (measured at zx) and Y (measured at zy) to mid-Z
                        double xm = xi.centre + truth_dxdz * (zm - zx);
                        double ym = yi.centre + truth_dydz * (zm - zy);
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

            // Tracker DSSSD: T1(X, key=0) + T2(Y, key=1)
            add_dsssd_pixel(0, 1, STRAW_HALF_LENGTH);
            // FRI DSSSD:     F1(X, key=10) + F2(Y, key=11)
            add_dsssd_pixel(10, 11, 70.25);

            // --- D: T3+T4 stereo pixel (proper u/v unfolding) ---
            // U = (x+y)/sqrt(2) [T3], V = (x-y)/sqrt(2) [T4]
            // => x = (u+v)/sqrt(2),  y = (u-v)/sqrt(2)
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

        // ----------------------------------------------------------------
        // Truth momentum slopes for arrival-angle correction in DSSSD pixel.
        // dx/dz = px/pz, dy/dz = py/pz from the kaon truth 4-vector.
        // ----------------------------------------------------------------
        double truth_dxdz = 0., truth_dydz = 0.;
        {
            auto slope_it = truth_map.find(event_number);
            if (slope_it != truth_map.end()) {
                const auto& tr = slope_it->second;
                if (std::abs(tr.pz) > 1e-6) {
                    truth_dxdz = tr.px / tr.pz;
                    truth_dydz = tr.py / tr.pz;
                }
            }
        }

        // ----------------------------------------------------------------
        // DEBUG: per-event hit summary (first DEBUG_MAX_EVENTS events only)
        // ----------------------------------------------------------------
        if (DEBUG_MODE && debug_event_count < DEBUG_MAX_EVENTS) {
            ++debug_event_count;
            std::cout << std::fixed << std::setprecision(1);
            std::cout << "\n=== DEBUG EVENT " << event_number
                      << " (" << debug_event_count << "/" << DEBUG_MAX_EVENTS << ") ==="
                      << "  [3-pion decay event]\n";
            std::cout << "  Reconstruction requires (DSSSD model):\n"
                         "    - PIZZA hit + TOF hit for each pion\n"
                         "    - >= 2 TrackPoints per pion (DSSSD pixels/FRI pixel)\n"
                         "    - at least 1 X-measuring module (T1/F1 sub-layer-0) AND\n"
                         "      at least 1 Y-measuring module (T2/F2 sub-layer-0) per pion\n"
                         "    - T1+T2 paired as DSSSD pixel; F1+F2 paired as FRI pixel;\n"
                         "      T3+T4 stereo unfolded (u+v->x, u-v->y) if both present\n";

            for (int pion = 0; pion < 2; ++pion) {
                const std::vector<HitInfo>& phits = (pion==0) ? hits_pip : hits_pim;
                bool has_pz = (pion==0) ? (pip_pizza_time>=0) : (pim_pizza_time>=0);
                bool has_tf = (pion==0) ? (pip_tof_time  >=0) : (pim_tof_time  >=0);
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
                std::cout << "    Module hit counts:"
                          << "  T1(X)=" << cnt_t1 << "  T2(Y)=" << cnt_t2
                          << "  T3(U)=" << cnt_t3 << "  T4(V)=" << cnt_t4
                          << "  F1(X)=" << cnt_f1 << "  F2(Y)=" << cnt_f2 << "\n";

                if (!uid_set.empty()) {
                    std::cout << "    Unique device IDs fired (" << uid_set.size() << "):\n";
                    for (int uid : uid_set)
                        std::cout << "      ID=" << std::setw(4) << uid
                                  << "  ->  " << device_id_info(uid) << "\n";
                }

                bool x_ok = (cnt_t1>0 || cnt_t3>0 || cnt_f1>0);
                bool y_ok = (cnt_t2>0 || cnt_t4>0 || cnt_f2>0);
                std::cout << "    Constraint status:"
                          << "  X=" << (x_ok?"OK    ":"MISSING")
                          << "  Y=" << (y_ok?"OK    ":"MISSING")
                          << "  PIZZA=" << (has_pz?"OK    ":"MISSING")
                          << "  TOF=" << (has_tf?"OK    ":"MISSING") << "\n";

                // Preview TrackPoints for this pion
                std::vector<TrackPoint> dbg_tps = build_track_points(phits, truth_dxdz, truth_dydz);
                std::cout << "    TrackPoints built: " << dbg_tps.size();
                for (const auto& tp : dbg_tps)
                    std::cout << "  [z=" << tp.z
                              << " x=" << tp.x << "+-" << tp.sx
                              << " y=" << tp.y << "+-" << tp.sy << "]";
                std::cout << "\n";

                if (dbg_tps.size() < 2)
                    std::cout << "    [FAIL-TRACKING: only " << dbg_tps.size()
                              << " TrackPoint(s) – need >= 2 distinct Z planes]\n";
                if (!x_ok)
                    std::cout << "    [FAIL-TRACKING: no X-measuring module (T1/T3/F1) hit]\n";
                if (!y_ok)
                    std::cout << "    [FAIL-TRACKING: no Y-measuring module (T2/T4/F2) hit]\n";
                if (!has_pz)
                    std::cout << "    [FAIL-TOF: no PIZZA hit for " << pn << "]\n";
                if (!has_tf)
                    std::cout << "    [FAIL-TOF: no TOF wall hit for " << pn << "]\n";
            }
        }

        // ----------------------------------------------------------------
        // RECONSTRUCTION CUT: require PIZZA + TOF hits for both pions
        // ----------------------------------------------------------------
        if (pip_pizza_time < 0 || pip_tof_time < 0 ||
            pim_pizza_time < 0 || pim_tof_time < 0) continue;

        std::vector<TrackPoint> track_points_pip = build_track_points(hits_pip, truth_dxdz, truth_dydz);
        std::vector<TrackPoint> track_points_pim = build_track_points(hits_pim, truth_dxdz, truth_dydz);

        // Need TrackPoints from at least 2 distinct detector planes.
        if (track_points_pip.size() < 2 || track_points_pim.size() < 2) continue;

        // ----------------------------------------------------------------
        // WEIGHTED LINE-OF-BEST-FIT for each pion track
        // ----------------------------------------------------------------
        TVector3 pip_dir   = fit_track_direction(track_points_pip);
        TVector3 pip_start = fit_track_origin(track_points_pip);

        TVector3 pim_dir   = fit_track_direction(track_points_pim);
        TVector3 pim_start = fit_track_origin(track_points_pim);

        // ----------------------------------------------------------------
        // DECAY VERTEX via PoCA of the two track lines
        // ----------------------------------------------------------------
        TVector3 decay_vertex = closest_point_between_lines(
            pip_start, pip_dir, pim_start, pim_dir);

        // ----------------------------------------------------------------
        // PIZZA hit positions — smeared as before (TOF start position).
        // [TO CHANGE] Remove Gaus() calls (or set smear_sigma=0) for ideal reco.
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
        // TOF wall hit positions — derived from bar geometry (no raw smear).
        //
        // X : bar centre from device ID (exact).
        // Y : bar nominal Y centre + Gaussian smear (sigma = TOF_Y_UNCERTAINTY)
        //     to model light-path timing uncertainty along the bar.
        //
        // [TO CHANGE] Switch to Uniform(-TOF_Y_UNCERTAINTY, +TOF_Y_UNCERTAINTY)
        //             for a flat prior on Y position.
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
        // Path uses track-fit extrapolation to z_tof (not the smeared ToF
        // hit position — TOF_Y_UNCERTAINTY = 10 cm would otherwise inflate
        // the 3D path length in quadrature on every event, causing a
        // systematic overestimate of pip_v and therefore a systematic
        // underestimate of kaon momentum via the back-propagation chain).
        // ----------------------------------------------------------------
        double t_pip_tof = (pip_dir.Z() != 0.) ? (pip_tof_z_raw - pip_start.Z()) / pip_dir.Z() : 0.;
        TVector3 pip_tof_fit = pip_start + t_pip_tof * pip_dir;
        double t_pim_tof = (pim_dir.Z() != 0.) ? (pim_tof_z_raw - pim_start.Z()) / pim_dir.Z() : 0.;
        TVector3 pim_tof_fit = pim_start + t_pim_tof * pim_dir;
        double pip_track_cm = (pip_tof_fit - pip_pizza_pos).Mag();
        double pim_track_cm = (pim_tof_fit - pim_pizza_pos).Mag();

        double pip_dt_ns = pip_tof_time - pip_pizza_time;
        double pim_dt_ns = pim_tof_time - pim_pizza_time;

        double pip_v = (pip_dt_ns > 0)
            ? (pip_track_cm * 1e-2) / (pip_dt_ns * 1e-9) : 0.; // m/s
        double pim_v = (pim_dt_ns > 0)
            ? (pim_track_cm * 1e-2) / (pim_dt_ns * 1e-9) : 0.;

        // ----------------------------------------------------------------
        // KAON DECAY TIME — back-propagate each pion from PIZZA to decay vertex.
        // [TO CHANGE] Use weighted average if measurement qualities differ.
        // ----------------------------------------------------------------
        double pip_path_cm    = (pip_pizza_pos - decay_vertex).Mag();
        double pip_decay_time = pip_pizza_time * 1e-9
                                - (pip_path_cm * 1e-2 / pip_v);  // s

        double pim_path_cm    = (pim_pizza_pos - decay_vertex).Mag();
        double pim_decay_time = pim_pizza_time * 1e-9
                                - (pim_path_cm * 1e-2 / pim_v);

        double kaon_decay_time = 0.5 * (pip_decay_time + pim_decay_time); // s

        // ----------------------------------------------------------------
        // KAON MOMENTUM  p = gamma * m_K * beta
        //
        // Assumes kaon produced at origin at t=0 and travels in a straight line.
        //
        // [TO CHANGE] Update kaon_prod if the production point is not the origin.
        // [TO CHANGE] m_K = K_L mass. Use 0.493677 GeV/c^2 for K+/K-.
        // ----------------------------------------------------------------
        TVector3 kaon_prod(0., 0., 0.); // [TO CHANGE] if not produced at origin

        double flight_cm = (decay_vertex - kaon_prod).Mag();
        double kaon_v    = (flight_cm * 1e-2) / kaon_decay_time; // m/s

        double m_K    = 0.497611; // GeV/c^2 — K_L mass
        double beta_K = kaon_v / 2.99792458e8;
        if (beta_K >= 1.0) beta_K = 0.9999;
        double gamma_K = 1.0 / std::sqrt(1. - beta_K * beta_K);
        double kaon_p  = gamma_K * m_K * beta_K; // GeV/c

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
        // RECONSTRUCTION CUT: momentum upper bound
        // Reject unphysically large reco momenta (reconstruction failures).
        // [TO CHANGE] Adjust to match your beam momentum range.
        // ----------------------------------------------------------------
        if (kaon_p > 11.) continue;

        std::cout << "Event " << event_number
                  << " | Reco p: " << kaon_p
                  << " | True p: " << true_p_mag
                  << " | Reco vertex: (" << decay_vertex.X()     << ", "
                                         << decay_vertex.Y()     << ", "
                                         << decay_vertex.Z()     << ")"
                  << " | True vertex: (" << true_vertex_vec.X() << ", "
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
    // [TO CHANGE] Add outTree->Branch() calls here to save additional quantities.
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

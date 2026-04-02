// ============================================================================
// vis_common.h — shared includes, constants, structs, and physics for the
//               KLong visualiser suite.
//
// Include this in every macro that does reconstruction but not drawing.
// vis_drawing.h already includes this, so don't double-include both manually.
//
// Geometry constants, detector IDs, and track-fit logic mirror
// KLong_save_vectors.C exactly.
// ============================================================================
#ifndef VIS_COMMON_H
#define VIS_COMMON_H

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TView.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TMarker.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArc.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TMultiGraph.h"
#include <fstream>
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
#include <sstream>
#include "TSystemDirectory.h"
#include "TList.h"
#include "TFileMerger.h"

// ============================================================================
// GEOMETRY CONSTANTS  (must match KLong_save_vectors.C)
// ============================================================================
static const double STRAW_HALF_WIDTH  = 0.40;
static const double STRAW_HALF_LENGTH = 45.0;
static const double FRI_STRIP_X[26] = {
    -66., -60., -54., -48., -42., -36., -30., -24., -18., -12., -7.5,
     -4.5, -1.5,  1.5,  4.5,  7.5,  12.,  18.,  24.,  30.,  36.,  42.,
     48.,  54.,  60.,  66.
};
static const double FRI_HALF_LENGTH[26] = {
    37.25, 42.75, 53.75, 59.25, 59.25, 64.75,
    70.25, 70.25, 70.25, 70.25,
    70.25, 70.25, 70.25, 70.25, 70.25, 70.25,
    70.25, 70.25, 70.25, 70.25,
    64.75, 59.25, 59.25, 53.75, 42.75, 37.25
};
static const double TOF_BAR_HALF_WIDTH = 6.0;
static const double TOF_Y_UNCERTAINTY  = 10.0;
static const double SQRT2_VIS = 1.41421356237;

// ============================================================================
// STRUCTS
// ============================================================================
struct HitInfoV { double x, y, z, t; int deviceID; };

struct TrackFitV {
    TVector3 dir, orig;
    bool valid;
};

struct EventReco {
    std::vector<HitInfoV> hits_pip, hits_pim;
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

// All data needed to display one event
struct VisEvent {
    int event_num;
    double reco_p, true_p, delta_p;
    // Vertices
    TVector3 reco_vtx, true_vtx;
    // Raw tracker/FRI hits
    std::vector<HitInfoV> hits_pip, hits_pim;
    // Raw PIZZA hit positions
    TVector3 pip_pizza_raw, pim_pizza_raw;
    // Bar-ID smeared TOF positions (used in reco)
    TVector3 pip_tof_pos, pim_tof_pos;
    // Track-extrapolated PIZZA positions (used in reco path)
    TVector3 pip_pizza_pos, pim_pizza_pos;
    // Track fits
    TrackFitV fit_pip, fit_pim;
    // MC truth kaon direction (unit vector)
    TVector3 true_kaon_dir;
    // Source file (seed) this event came from
    std::string source_file;

    // --- Kinematic intermediate quantities (for print_event_detail) ---
    double pip_pizza_time = 0, pim_pizza_time = 0;
    double pip_tof_time   = 0, pim_tof_time   = 0;
    int    pip_tof_devID  = -1, pim_tof_devID = -1;
    double pip_track_cm   = 0, pim_track_cm   = 0;  // PIZZA→TOF path length
    double pip_dt_ns      = 0, pim_dt_ns      = 0;  // TOF - PIZZA timing
    double pip_v          = 0, pim_v          = 0;  // pion velocity (m/s)
    double pip_decay_path_cm = 0, pim_decay_path_cm = 0; // PIZZA→vertex
    double pip_decay_time = 0, pim_decay_time = 0;  // kaon decay time from each
    double kaon_decay_time = 0;
    double kaon_flight_cm = 0;
    double kaon_v = 0, kaon_beta = 0, kaon_gamma = 0;
};

// ============================================================================
// GEOMETRY HELPERS
// ============================================================================
double fri_hw(int i) { return (i >= 10 && i <= 15) ? 1.5 : 3.0; }
double tof_bar_xc(int id) { return -96.0 + (id - 2053) * 12.0; }
double tof_bar_yc(int id) {
    int i = id - 2053;
    if (i == 8) return -72.0;
    if (i == 9) return +72.0;
    return 0.0;
}

// ============================================================================
// 3D TRACK FIT  (identical to KLong_save_vectors.C)
// ============================================================================
TrackFitV fit_track_vis(const std::vector<HitInfoV>& hits, int tof_devID, double z_tof)
{
    const TrackFitV FAIL = {{0,0,1},{0,0,0},false};
    std::vector<double> xz, xv, xe, yz, yv, ye;

    for (const auto& h : hits) {
        int id = h.deviceID;
        if (id >= 1 && id <= 488) {
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(STRAW_HALF_WIDTH);
        } else if (id >= 489 && id <= 976) {
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(STRAW_HALF_WIDTH);
        } else if (id >= 977 && id <= 1464) {
            double sxy = STRAW_HALF_WIDTH * SQRT2_VIS;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(sxy);
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(sxy);
        } else if (id >= 1465 && id <= 1952) {
            double sxy = STRAW_HALF_WIDTH * SQRT2_VIS;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(sxy);
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(sxy);
        } else if (id >= 2001 && id <= 2026) {
            int si = id - 2001;
            xz.push_back(h.z); xv.push_back(h.x); xe.push_back(fri_hw(si));
        } else if (id >= 2027 && id <= 2052) {
            int si = id - 2027;
            yz.push_back(h.z); yv.push_back(h.y); ye.push_back(fri_hw(si));
        }
        // TOF excluded from track fit (large bar uncertainties bias slope at high lever arm)
    }
    // tof_devID and z_tof kept in signature for call-site compatibility.
    (void)tof_devID; (void)z_tof;
    if ((int)xz.size() < 2 || (int)yz.size() < 2) return FAIL;
    auto distinct = [](const std::vector<double>& v) {
        for (size_t k=1; k<v.size(); ++k) if (std::fabs(v[k]-v[0]) > 1.0) return true;
        return false;
    };
    if (!distinct(xz) || !distinct(yz)) return FAIL;

    TGraphErrors gx((int)xz.size(), xz.data(), xv.data(), nullptr, xe.data());
    gx.Fit("pol1","Q"); TF1* fx = gx.GetFunction("pol1"); if (!fx) return FAIL;
    double bx = fx->GetParameter(0), ax = fx->GetParameter(1);

    TGraphErrors gy((int)yz.size(), yz.data(), yv.data(), nullptr, ye.data());
    gy.Fit("pol1","Q"); TF1* fy = gy.GetFunction("pol1"); if (!fy) return FAIL;
    double by = fy->GetParameter(0), ay = fy->GetParameter(1);

    if (!std::isfinite(ax)||!std::isfinite(ay)||!std::isfinite(bx)||!std::isfinite(by)) return FAIL;
    TVector3 dir(ax, ay, 1.0); dir = dir.Unit();
    return {dir, {bx, by, 0.}, true};
}

TVector3 poca_vis(const TVector3& p1, const TVector3& v1, const TVector3& p2, const TVector3& v2)
{
    TVector3 w0 = p1 - p2;
    double a = v1.Dot(v1), b = v1.Dot(v2), c = v2.Dot(v2);
    double d = v1.Dot(w0), e = v2.Dot(w0);
    double denom = a*c - b*b;
    if (std::abs(denom) < 1e-14) return 0.5*(p1+p2);
    double sc = (b*e - c*d)/denom, tc = (a*e - b*d)/denom;
    return 0.5*((p1 + v1*sc) + (p2 + v2*tc));
}

TVector3 eval_at_z(const TrackFitV& f, double z) {
    double tp = (z - f.orig.Z()) / f.dir.Z();
    return f.orig + tp * f.dir;
}

#endif // VIS_COMMON_H

// ============================================================================
// KLong_visualize_event.C — 3D reconstruction event display
//
// PURPOSE:
//   For a given simulation .root file:
//     1. Re-runs the full K_L -> pi+pi-pi0 reconstruction on every event
//     2. Writes a text file of all reconstructed events sorted by |delta_p|
//        (largest bias first) so outliers are easy to locate
//     3. For a specified event number draws a 4-panel figure:
//          Panel 1 (top-left):  3D perspective view — detector planes, hit
//                               positions, fitted track lines, vertices
//          Panel 2 (top-right): XZ projection ("side view")
//          Panel 3 (bottom-left): YZ projection ("top view")
//          Panel 4 (bottom-right): XY projection ("front view")
//
// INPUT:  a simulation .root file with Ntuple1, Ntuple2, Ntuple3
//         (e.g., T1-240_..._E1-700_1.root, NOT the combined_ file)
// OUTPUT: <stem>_delta_p_sorted.txt  — events sorted by |delta_p|
//         KLong_event<N>_view.png    — 4-panel figure for event N
//
// USAGE:
//   Generate text file only (no display):
//     root -l -q 'KLong_visualize_event.C("file.root")'
//
//   Generate text file + display event 42 + print step-by-step breakdown:
//     root -l -q 'KLong_visualize_event.C("file.root", 42)'
//
//   Custom output filename:
//     root -l -q 'KLong_visualize_event.C("file.root", 42, "outliers.txt")'
//
//   Run over all seed files in a config directory (combined sorted output):
//     root -l -q 'KLong_visualize_config("/path/to/config_dir")'         // text only
//     root -l -q 'KLong_visualize_config("/path/to/config_dir", 46480)'  // + display
//     root -l -q 'KLong_visualize_config("/path/to/config_dir", 46480, "out.txt")'
//   The directory must contain raw simulation files (*_N.root, NOT *_vectors/acceptance*).
//   Output text file has an extra "seed" column identifying which file each event came from.
//   The PNG is saved as KLong_event<N>_<stem>_view.png.
//
//   Find events where both pion tracks lie in a coordinate plane (verification):
//     root -l -q 'KLong_find_simple_events("file.root")'
//     root -l -q 'KLong_find_simple_events("file.root", 5.0, "simple.txt")'
//       phi_tol_deg (default 10): max allowed phi deviation from XZ or YZ plane
//       Output lists events sorted by how aligned they are; events marked * pass tol.
//       Event with phi ~ 0/180 => both tracks in XZ plane (y~0, 2D path check)
//       Event with phi ~ +-90  => both tracks in YZ plane (x~0, 2D path check)
//
// Geometry and reconstruction logic mirrors KLong_save_vectors.C exactly.
// ============================================================================

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
// 3D TRACK FIT (identical to KLong_save_vectors.C)
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

// ============================================================================
// 3D DRAWING HELPERS
// ============================================================================

// Draw a rectangular wireframe at fixed z
void draw_rect3d(double xmin, double xmax, double ymin, double ymax, double z,
                 int color, int width = 1, int style = 1)
{
    TPolyLine3D* frame = new TPolyLine3D(5);
    frame->SetPoint(0, xmin, ymin, z);
    frame->SetPoint(1, xmax, ymin, z);
    frame->SetPoint(2, xmax, ymax, z);
    frame->SetPoint(3, xmin, ymax, z);
    frame->SetPoint(4, xmin, ymin, z);
    frame->SetLineColor(color); frame->SetLineWidth(width); frame->SetLineStyle(style);
    frame->Draw();
}

// Draw a circular ring at fixed z (outer radius only, or inner+outer)
void draw_ring3d(double z, double r_outer, int color, int n_seg = 36)
{
    TPolyLine3D* circle = new TPolyLine3D(n_seg + 1);
    for (int i = 0; i <= n_seg; i++) {
        double phi = 2. * M_PI * i / n_seg;
        circle->SetPoint(i, r_outer * std::cos(phi), r_outer * std::sin(phi), z);
    }
    circle->SetLineColor(color); circle->SetLineWidth(2);
    circle->Draw();
}

// Draw a 3D straight line between two z-planes along a TrackFit
void draw_track3d(const TrackFitV& f, double z_start, double z_end,
                  int color, int style = 1, int width = 2)
{
    TVector3 p0 = eval_at_z(f, z_start);
    TVector3 p1 = eval_at_z(f, z_end);
    TPolyLine3D* line = new TPolyLine3D(2);
    line->SetPoint(0, p0.X(), p0.Y(), p0.Z());
    line->SetPoint(1, p1.X(), p1.Y(), p1.Z());
    line->SetLineColor(color); line->SetLineStyle(style); line->SetLineWidth(width);
    line->Draw();
}

// Draw individual TOF bars
void draw_tof_bars3d(double z_tof, int color)
{
    for (int id = 2053; id <= 2070; id++) {
        double xc = tof_bar_xc(id);
        double yc = tof_bar_yc(id);
        double bar_half_y = (id == 2061 || id == 2062) ? 2.0 : 60.0; // gap bars shorter
        draw_rect3d(xc - TOF_BAR_HALF_WIDTH, xc + TOF_BAR_HALF_WIDTH,
                    yc - bar_half_y,          yc + bar_half_y,
                    z_tof, color, 1);
    }
}

// ============================================================================
// 3D PANEL DRAWING
// ============================================================================
void draw_3d_panel(TPad* pad, const VisEvent& ev,
                   double z_t1, double z_t2, double z_t3, double z_t4,
                   double z_pizza1, double z_pizza2,
                   double z_fri1,   double z_fri2,
                   double z_tof)
{
    pad->cd();
    pad->SetFillColor(kBlack);  // dark background for 3D

    // Create 3D view
    double margin = 20.;
    double xext = 110., yext = 80.;
    double rmin[3] = { -xext, -yext, -30. };
    double rmax[3] = {  xext,  yext,  z_tof + margin };
    TView* view = TView::CreateView(1, rmin, rmax);
    Int_t irep;
    view->SetView(30., 75., 0., irep);  // phi=30, theta=75 (slightly from above)
    view->ShowAxis();

    // ---- Detector planes ----
    // Tracker stations (gray wireframes, ±STRAW_HALF_LENGTH = 45 cm extent)
    const double TRK = STRAW_HALF_LENGTH; // 45 cm
    int col_trk = kGray+1;
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t1, col_trk, 1);
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t2, col_trk, 1);
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t3, col_trk, 1);
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t4, col_trk, 1);

    // PIZZA detectors (orange rings, r ~ 20 cm)
    draw_ring3d(z_pizza1, 20., kOrange+1);
    draw_ring3d(z_pizza2, 20., kOrange+1);

    // FRI walls (cyan wireframes, x: -66 to +66, y: ±70 cm)
    draw_rect3d(-66., 66., -70.25, 70.25, z_fri1, kCyan+1, 1);
    draw_rect3d(-66., 66., -70.25, 70.25, z_fri2, kCyan+1, 1);

    // TOF wall bars (magenta)
    draw_tof_bars3d(z_tof, kMagenta+1);

    // ---- Kaon line: origin to true vertex ----
    {
        TVector3 kend = ev.true_vtx;
        TPolyLine3D* klin = new TPolyLine3D(2);
        klin->SetPoint(0, 0., 0., 0.);
        klin->SetPoint(1, kend.X(), kend.Y(), kend.Z());
        klin->SetLineColor(kYellow); klin->SetLineStyle(2); klin->SetLineWidth(2);
        klin->Draw();
    }

    // ---- pi+ track line (red, from z_pizza1-5 to z_tof+5) ----
    if (ev.fit_pip.valid)
        draw_track3d(ev.fit_pip, z_pizza1 - 5., z_tof + 5., kRed,   1, 2);
    // ---- pi- track line (blue) ----
    if (ev.fit_pim.valid)
        draw_track3d(ev.fit_pim, z_pizza1 - 5., z_tof + 5., kBlue+1, 1, 2);

    // ---- Tracker + FRI hits pi+ (red crosses) ----
    if (!ev.hits_pip.empty()) {
        TPolyMarker3D* m = new TPolyMarker3D((int)ev.hits_pip.size(), 2);
        for (int i = 0; i < (int)ev.hits_pip.size(); i++)
            m->SetPoint(i, ev.hits_pip[i].x, ev.hits_pip[i].y, ev.hits_pip[i].z);
        m->SetMarkerColor(kRed); m->SetMarkerStyle(2); m->SetMarkerSize(1.2); m->Draw();
    }
    // ---- Tracker + FRI hits pi- (blue crosses) ----
    if (!ev.hits_pim.empty()) {
        TPolyMarker3D* m = new TPolyMarker3D((int)ev.hits_pim.size(), 2);
        for (int i = 0; i < (int)ev.hits_pim.size(); i++)
            m->SetPoint(i, ev.hits_pim[i].x, ev.hits_pim[i].y, ev.hits_pim[i].z);
        m->SetMarkerColor(kBlue+1); m->SetMarkerStyle(2); m->SetMarkerSize(1.2); m->Draw();
    }

    // ---- PIZZA hit positions (star markers) ----
    {
        TPolyMarker3D* mpz = new TPolyMarker3D(1, 29);  // 29 = star
        mpz->SetPoint(0, ev.pip_pizza_raw.X(), ev.pip_pizza_raw.Y(), ev.pip_pizza_raw.Z());
        mpz->SetMarkerColor(kRed); mpz->SetMarkerSize(2.); mpz->Draw();
    }
    {
        TPolyMarker3D* mpz = new TPolyMarker3D(1, 29);
        mpz->SetPoint(0, ev.pim_pizza_raw.X(), ev.pim_pizza_raw.Y(), ev.pim_pizza_raw.Z());
        mpz->SetMarkerColor(kBlue+1); mpz->SetMarkerSize(2.); mpz->Draw();
    }

    // ---- TOF smeared positions (square markers) ----
    {
        TPolyMarker3D* mt = new TPolyMarker3D(1, 21);
        mt->SetPoint(0, ev.pip_tof_pos.X(), ev.pip_tof_pos.Y(), ev.pip_tof_pos.Z());
        mt->SetMarkerColor(kRed); mt->SetMarkerSize(2.); mt->Draw();
    }
    {
        TPolyMarker3D* mt = new TPolyMarker3D(1, 21);
        mt->SetPoint(0, ev.pim_tof_pos.X(), ev.pim_tof_pos.Y(), ev.pim_tof_pos.Z());
        mt->SetMarkerColor(kBlue+1); mt->SetMarkerSize(2.); mt->Draw();
    }

    // ---- Reconstructed decay vertex (large green filled circle) ----
    {
        TPolyMarker3D* mv = new TPolyMarker3D(1, 20);
        mv->SetPoint(0, ev.reco_vtx.X(), ev.reco_vtx.Y(), ev.reco_vtx.Z());
        mv->SetMarkerColor(kGreen+2); mv->SetMarkerSize(3.); mv->Draw();
    }
    // ---- True decay vertex (large open circle, black) ----
    {
        TPolyMarker3D* mv = new TPolyMarker3D(1, 24);
        mv->SetPoint(0, ev.true_vtx.X(), ev.true_vtx.Y(), ev.true_vtx.Z());
        mv->SetMarkerColor(kWhite); mv->SetMarkerSize(3.); mv->Draw();
    }

    // Legend as TLatex in NDC coordinates
    TLatex ltx;
    ltx.SetNDC(true);
    ltx.SetTextSize(0.033);
    ltx.SetTextColor(kWhite);
    ltx.DrawLatex(0.02, 0.97, Form("Event %d   p_{true}=%.3f GeV/c   p_{reco}=%.3f GeV/c   #Deltap=%.3f GeV/c",
                                   ev.event_num, ev.true_p, ev.reco_p, ev.delta_p));
    ltx.SetTextColor(kRed);         ltx.DrawLatex(0.02, 0.91, "#pi^{+} trk hits");
    ltx.SetTextColor(kBlue+1);      ltx.DrawLatex(0.15, 0.91, "#pi^{-} trk hits");
    ltx.SetTextColor(kGreen+2);     ltx.DrawLatex(0.22, 0.91, "Reco vtx");
    ltx.SetTextColor(kWhite);       ltx.DrawLatex(0.34, 0.91, "True vtx");
    ltx.SetTextColor(kYellow);      ltx.DrawLatex(0.44, 0.91, "Kaon path");
    ltx.SetTextColor(kGray+1);      ltx.DrawLatex(0.58, 0.91, "Trackers");
    ltx.SetTextColor(kCyan+1);      ltx.DrawLatex(0.68, 0.91, "FRI");
    ltx.SetTextColor(kMagenta+1);   ltx.DrawLatex(0.74, 0.91, "TOF");
    ltx.SetTextColor(kOrange+1);    ltx.DrawLatex(0.83, 0.91, "PIZZA");
}

// ============================================================================
// 2D PROJECTION PANEL
// ============================================================================
// proj: "XZ", "YZ", or "XY"
void draw_proj_panel(TPad* pad, const VisEvent& ev,
                     const char* proj,
                     double z_t1, double z_t2, double z_t3, double z_t4,
                     double z_pizza1, double z_pizza2,
                     double z_fri1,   double z_fri2,
                     double z_tof)
{
    pad->cd();
    pad->SetFillColor(10);  // white background for 2D
    pad->SetLeftMargin(0.12); pad->SetBottomMargin(0.12);

    bool is_xz = (std::string(proj) == "XZ");
    bool is_yz = (std::string(proj) == "YZ");
    bool is_xy = (std::string(proj) == "XY");

    // Axis ranges
    double hmin = 0, hmax = 0, vmin = 0, vmax = 0;
    const char *hlab = "", *vlab = "", *title = "";
    if (is_xz) { hmin=0; hmax=z_tof+20; vmin=-110; vmax=110; hlab="z (cm)"; vlab="x (cm)"; title="Side view (XZ)"; }
    if (is_yz) { hmin=0; hmax=z_tof+20; vmin=-90;  vmax=90;  hlab="z (cm)"; vlab="y (cm)"; title="Top view (YZ)"; }
    if (is_xy) { hmin=-110; hmax=110; vmin=-90; vmax=90; hlab="x (cm)"; vlab="y (cm)"; title="Front view (XY)"; }

    TH2F* frame = new TH2F(Form("hframe_%s_%d", proj, ev.event_num),
                           Form("%s;%s;%s", title, hlab, vlab),
                           1, hmin, hmax, 1, vmin, vmax);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.048);
    frame->GetYaxis()->SetLabelSize(0.048);
    frame->Draw();

    // Helper to extract (h,v) from a 3D point
    auto hv = [&](const TVector3& p) -> std::pair<double,double> {
        if (is_xz) return {p.Z(), p.X()};
        if (is_yz) return {p.Z(), p.Y()};
        return {p.X(), p.Y()};
    };

    // ---- Detector z-lines (only for XZ and YZ) ----
    if (!is_xy) {
        auto draw_det_line = [&](double z, int color, int style) {
            TLine* l = new TLine(z, vmin, z, vmax);
            l->SetLineColor(color); l->SetLineStyle(style); l->SetLineWidth(1);
            l->Draw();
        };
        draw_det_line(z_pizza1, kOrange+1, 2);
        draw_det_line(z_pizza2, kOrange+1, 2);
        draw_det_line(z_t1,  kGray+1, 1);
        draw_det_line(z_t2,  kGray+1, 1);
        draw_det_line(z_fri1, kCyan+1, 1);
        draw_det_line(z_fri2, kCyan+1, 1);
        draw_det_line(z_t3,  kGray+1, 1);
        draw_det_line(z_t4,  kGray+1, 1);
        draw_det_line(z_tof, kMagenta+1, 1);
        // Labels
        TLatex tl; tl.SetTextSize(0.038); tl.SetTextAngle(90); tl.SetTextAlign(12);
        tl.SetTextColor(kGray+2);
        tl.DrawLatex(z_t1,  vmax*0.6, "T1");
        tl.DrawLatex(z_t2,  vmax*0.6, "T2");
        tl.DrawLatex(z_fri1,vmax*0.6, "F1");
        tl.DrawLatex(z_fri2,vmax*0.6, "F2");
        tl.DrawLatex(z_t3,  vmax*0.6, "T3");
        tl.DrawLatex(z_t4,  vmax*0.6, "T4");
        tl.SetTextColor(kMagenta+1);
        tl.DrawLatex(z_tof, vmax*0.6, "TOF");
        tl.SetTextColor(kOrange+1);
        tl.DrawLatex(z_pizza1, vmax*0.6, "P1");
        tl.DrawLatex(z_pizza2, vmax*0.6, "P2");
    }
    // For XY: draw detector circles/boxes
    if (is_xy) {
        // PIZZA ring
        TEllipse* e = new TEllipse(0., 0., 20., 20.); e->SetFillStyle(0);
        e->SetLineColor(kOrange+1); e->SetLineWidth(2); e->Draw("same");
        // Tracker outline
        TBox* bt = new TBox(-45., -45., 45., 45.);
        bt->SetFillStyle(0); bt->SetLineColor(kGray+1); bt->Draw("l same");
        // FRI outline
        TBox* bf = new TBox(-66., -70.25, 66., 70.25);
        bf->SetFillStyle(0); bf->SetLineColor(kCyan+1); bf->Draw("l same");
    }

    // ---- Tracker + FRI hits ----
    auto make_hit_graph = [&](const std::vector<HitInfoV>& hits, int color, int marker) -> TGraph* {
        TGraph* g = new TGraph();
        g->SetMarkerColor(color); g->SetMarkerStyle(marker); g->SetMarkerSize(0.9);
        for (const auto& h : hits) {
            auto [ph, pv] = hv(TVector3(h.x, h.y, h.z));
            g->AddPoint(ph, pv);
        }
        return g;
    };
    TGraph* gp = make_hit_graph(ev.hits_pip, kRed,    2);  // 2 = cross
    TGraph* gm = make_hit_graph(ev.hits_pim, kBlue+1, 2);
    gp->Draw("P same");
    gm->Draw("P same");

    // ---- PIZZA hit positions ----
    {
        auto [h,v] = hv(ev.pip_pizza_raw);
        TMarker* m = new TMarker(h, v, 29); m->SetMarkerColor(kRed); m->SetMarkerSize(2.0); m->Draw();
    }
    {
        auto [h,v] = hv(ev.pim_pizza_raw);
        TMarker* m = new TMarker(h, v, 29); m->SetMarkerColor(kBlue+1); m->SetMarkerSize(2.0); m->Draw();
    }

    // ---- TOF smeared positions ----
    {
        auto [h,v] = hv(ev.pip_tof_pos);
        TMarker* m = new TMarker(h, v, 21); m->SetMarkerColor(kRed); m->SetMarkerSize(1.5); m->Draw();
    }
    {
        auto [h,v] = hv(ev.pim_tof_pos);
        TMarker* m = new TMarker(h, v, 21); m->SetMarkerColor(kBlue+1); m->SetMarkerSize(1.5); m->Draw();
    }

    // ---- Track projection lines ----
    if (ev.fit_pip.valid && !is_xy) {
        double z0 = z_pizza1 - 5., z1 = z_tof + 5.;
        auto [ph0,pv0] = hv(eval_at_z(ev.fit_pip, z0));
        auto [ph1,pv1] = hv(eval_at_z(ev.fit_pip, z1));
        TLine* l = new TLine(ph0, pv0, ph1, pv1);
        l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw();
    }
    if (ev.fit_pim.valid && !is_xy) {
        double z0 = z_pizza1 - 5., z1 = z_tof + 5.;
        auto [ph0,pv0] = hv(eval_at_z(ev.fit_pim, z0));
        auto [ph1,pv1] = hv(eval_at_z(ev.fit_pim, z1));
        TLine* l = new TLine(ph0, pv0, ph1, pv1);
        l->SetLineColor(kBlue+1); l->SetLineWidth(2); l->Draw();
    }
    // XY: draw track extrapolated to a few detector z-planes
    if (is_xy && ev.fit_pip.valid) {
        double dzvals[] = {z_t1, z_t2, z_fri1, z_fri2, z_t3, z_t4};
        TGraph* gline_p = new TGraph();
        TGraph* gline_m = new TGraph();
        gline_p->SetLineColor(kRed);    gline_p->SetLineWidth(2);
        gline_m->SetLineColor(kBlue+1); gline_m->SetLineWidth(2);
        for (double dz : dzvals) {
            auto [h, v] = hv(eval_at_z(ev.fit_pip, dz));
            gline_p->AddPoint(h, v);
            auto [h2, v2] = hv(eval_at_z(ev.fit_pim, dz));
            gline_m->AddPoint(h2, v2);
        }
        gline_p->Draw("L same");
        gline_m->Draw("L same");
    }

    // ---- Vertices ----
    {
        auto [h,v] = hv(ev.reco_vtx);
        TMarker* m = new TMarker(h, v, 20); m->SetMarkerColor(kGreen+2); m->SetMarkerSize(2.2); m->Draw();
    }
    {
        auto [h,v] = hv(ev.true_vtx);
        TMarker* m = new TMarker(h, v, 24); m->SetMarkerColor(kBlack); m->SetMarkerSize(2.2); m->Draw();
    }

    // Kaon origin-to-true-vertex line
    if (!is_xy) {
        auto [h0, v0] = hv(TVector3(0., 0., 0.));
        auto [h1, v1] = hv(ev.true_vtx);
        TLine* kl = new TLine(h0, v0, h1, v1);
        kl->SetLineColor(kOrange-3); kl->SetLineStyle(2); kl->SetLineWidth(1); kl->Draw();
    }

    // Panel title
    TLatex ltx;
    ltx.SetNDC(true); ltx.SetTextSize(0.052); ltx.SetTextFont(62);
    ltx.DrawLatex(0.15, 0.94, title);
}

// ============================================================================
// DETAILED EVENT BREAKDOWN
//
// Prints a step-by-step numerical breakdown to stdout: all detector hits with
// their fit contributions, track fit parameters, theta/phi of each track,
// PIZZA/TOF positions, path lengths, velocities, decay times, and the full
// kaon momentum derivation.  Use this to verify the reconstruction by hand.
// ============================================================================
void print_event_detail(const VisEvent& ev)
{
    // Device name and fit role helpers
    auto det_name  = [](int id) -> const char* {
        if (id >=    1 && id <=  488) return "T1 (X-meas)   ";
        if (id >=  489 && id <=  976) return "T2 (Y-meas)   ";
        if (id >=  977 && id <= 1464) return "T3 (stereo+45)";
        if (id >= 1465 && id <= 1952) return "T4 (stereo-45)";
        if (id >= 1953 && id <= 2000) return "PIZZA         ";
        if (id >= 2001 && id <= 2026) return "FRI-W1 (X)    ";
        if (id >= 2027 && id <= 2052) return "FRI-W2 (Y)    ";
        if (id >= 2053 && id <= 2070) return "TOF bar       ";
        return "UNKNOWN       ";
    };
    auto fit_axis = [](int id) -> const char* {
        if (id >=    1 && id <=  488) return "x(z) ";
        if (id >=  489 && id <=  976) return "y(z) ";
        if (id >=  977 && id <= 1952) return "x+y  ";
        if (id >= 2001 && id <= 2026) return "x(z) ";
        if (id >= 2027 && id <= 2052) return "y(z) ";
        if (id >= 2053 && id <= 2070) return "excl ";
        return "?    ";
    };
    auto sigma_str = [](int id) -> std::string {
        if (id >=    1 && id <=  488) return "0.400 cm (x)  ";
        if (id >=  489 && id <=  976) return "0.400 cm (y)  ";
        if (id >=  977 && id <= 1952) return "0.566 cm (x+y)";
        if (id >= 2001 && id <= 2026) {
            int si = id-2001; return (si>=10&&si<=15) ? "1.500 cm (x)  " : "3.000 cm (x)  ";
        }
        if (id >= 2027 && id <= 2052) {
            int si = id-2027; return (si>=10&&si<=15) ? "1.500 cm (y)  " : "3.000 cm (y)  ";
        }
        if (id >= 2053 && id <= 2070) return "excluded      ";
        return "?             ";
    };

    std::cout << std::fixed;
    const std::string bar(78, '=');
    std::cout << "\n" << bar << "\n";
    std::cout << "  DETAILED RECONSTRUCTION BREAKDOWN  —  Event " << ev.event_num << "\n";
    std::cout << bar << "\n";
    std::cout << std::setprecision(4);
    std::cout << "  p_true      = " << ev.true_p  << " GeV/c\n";
    std::cout << "  p_reco      = " << ev.reco_p  << " GeV/c\n";
    std::cout << "  delta_p     = " << ev.delta_p << " GeV/c\n";
    std::cout << "  Reco vertex : (x=" << ev.reco_vtx.X() << ", y="
              << ev.reco_vtx.Y() << ", z=" << ev.reco_vtx.Z() << ") cm\n";
    std::cout << "  True vertex : (x=" << ev.true_vtx.X() << ", y="
              << ev.true_vtx.Y() << ", z=" << ev.true_vtx.Z() << ") cm\n";
    std::cout << "  Vtx dz (reco-true) = " << ev.reco_vtx.Z()-ev.true_vtx.Z()
              << " cm  |dr| = " << (ev.reco_vtx-ev.true_vtx).Mag() << " cm\n";

    for (int pion = 0; pion < 2; pion++) {
        const std::vector<HitInfoV>& hits = pion==0 ? ev.hits_pip : ev.hits_pim;
        const char* pn = pion==0 ? "PI+" : "PI-";
        int cnt[7] = {};
        for (const auto& h : hits) {
            if      (h.deviceID >=    1 && h.deviceID <=  488) cnt[0]++;
            else if (h.deviceID >=  489 && h.deviceID <=  976) cnt[1]++;
            else if (h.deviceID >=  977 && h.deviceID <= 1464) cnt[2]++;
            else if (h.deviceID >= 1465 && h.deviceID <= 1952) cnt[3]++;
            else if (h.deviceID >= 2001 && h.deviceID <= 2026) cnt[5]++;
            else if (h.deviceID >= 2027 && h.deviceID <= 2052) cnt[6]++;
        }
        std::cout << "\n  ── " << pn << " TRACKER/FRI HITS  ("
                  << hits.size() << " total: T1=" << cnt[0] << " T2=" << cnt[1]
                  << " T3=" << cnt[2] << " T4=" << cnt[3]
                  << " F1=" << cnt[5] << " F2=" << cnt[6] << ") ──\n";
        std::cout << "  " << std::setw(5)  << "DevID"
                  << "  " << std::setw(14) << "Detector"
                  << "  " << std::setw(9)  << "x (cm)"
                  << "  " << std::setw(9)  << "y (cm)"
                  << "  " << std::setw(9)  << "z (cm)"
                  << "  " << std::setw(5)  << "Fit"
                  << "  Sigma\n";
        std::cout << "  " << std::string(70, '-') << "\n";
        std::vector<HitInfoV> sh = hits;
        std::sort(sh.begin(), sh.end(), [](const HitInfoV& a, const HitInfoV& b){ return a.z < b.z; });
        for (const auto& h : sh) {
            std::cout << "  " << std::setw(5)  << h.deviceID
                      << "  " << std::setw(14) << det_name(h.deviceID)
                      << "  " << std::setw(9)  << std::setprecision(3) << h.x
                      << "  " << std::setw(9)  << h.y
                      << "  " << std::setw(9)  << h.z
                      << "  " << std::setw(5)  << fit_axis(h.deviceID)
                      << "  " << sigma_str(h.deviceID) << "\n";
        }
        std::cout << std::setprecision(4);

        // PIZZA
        const TVector3& piz_raw = pion==0 ? ev.pip_pizza_raw  : ev.pim_pizza_raw;
        const TVector3& piz_pos = pion==0 ? ev.pip_pizza_pos  : ev.pim_pizza_pos;
        double piz_t  = pion==0 ? ev.pip_pizza_time : ev.pim_pizza_time;
        std::cout << "  PIZZA raw hit      : (" << piz_raw.X() << ", " << piz_raw.Y()
                  << ", " << piz_raw.Z() << ") cm   t_smeared = " << piz_t << " ns\n";
        std::cout << "  PIZZA track extrap : (" << piz_pos.X() << ", " << piz_pos.Y()
                  << ", " << piz_pos.Z() << ") cm   [used for path & decay time]\n";

        // TOF
        const TVector3& tof_pos = pion==0 ? ev.pip_tof_pos  : ev.pim_tof_pos;
        int   tof_dev = pion==0 ? ev.pip_tof_devID   : ev.pim_tof_devID;
        double tof_t  = pion==0 ? ev.pip_tof_time    : ev.pim_tof_time;
        std::cout << "  TOF bar devID      : " << tof_dev
                  << "   bar_centre_x = " << tof_bar_xc(tof_dev)
                  << " cm   bar_centre_y = " << tof_bar_yc(tof_dev) << " cm\n";
        std::cout << "  TOF smeared pos    : (" << tof_pos.X() << ", " << tof_pos.Y()
                  << ", " << tof_pos.Z() << ") cm   t_smeared = " << tof_t << " ns\n";

        // Track fit
        const TrackFitV& fit = pion==0 ? ev.fit_pip : ev.fit_pim;
        double ax = fit.dir.X() / fit.dir.Z();
        double ay = fit.dir.Y() / fit.dir.Z();
        double phi_deg   = std::atan2(fit.dir.Y(), fit.dir.X()) * 180. / M_PI;
        double theta_deg = std::acos(std::min(1.0, std::fabs(fit.dir.Z()))) * 180. / M_PI;
        std::cout << "  Track fit  x(z)    = " << fit.orig.X() << " + (" << ax << ") * z\n";
        std::cout << "  Track fit  y(z)    = " << fit.orig.Y() << " + (" << ay << ") * z\n";
        std::cout << "  Direction (unit)   : (" << fit.dir.X() << ", "
                  << fit.dir.Y() << ", " << fit.dir.Z() << ")\n";
        std::cout << "  Track theta        : " << theta_deg
                  << " deg  (angle from Z axis; 0=beam-parallel)\n";
        std::cout << "  Track phi          : " << phi_deg
                  << " deg  (azimuth; 0/180=XZ-plane, +-90=YZ-plane)\n";

        // Velocity chain
        double L   = pion==0 ? ev.pip_track_cm      : ev.pim_track_cm;
        double dt  = pion==0 ? ev.pip_dt_ns          : ev.pim_dt_ns;
        double v   = pion==0 ? ev.pip_v              : ev.pim_v;
        double dp  = pion==0 ? ev.pip_decay_path_cm  : ev.pim_decay_path_cm;
        double td  = pion==0 ? ev.pip_decay_time     : ev.pim_decay_time;
        double beta_pi = v / 2.99792458e8;
        double gamma_pi = (beta_pi < 1.) ? 1./std::sqrt(1.-beta_pi*beta_pi) : 999.;
        std::cout << "  Path PIZZA→TOF     : " << L  << " cm\n";
        std::cout << "  delta_t (TOF-PIZZA): " << dt << " ns\n";
        std::cout << "  Pion velocity      : " << std::setprecision(8) << v << " m/s"
                  << "  (beta = " << beta_pi << ")\n";
        std::cout << "  Pion gamma         : " << std::setprecision(4) << gamma_pi << "\n";
        std::cout << "  Pion p (derived)   : "
                  << gamma_pi * 0.13957 * beta_pi << " GeV/c\n";
        std::cout << "  Path PIZZA→vertex  : " << dp << " cm\n";
        std::cout << "  Decay time from " << pn << " : "
                  << std::setprecision(11) << td << " s\n";
    }

    std::cout << "\n  ── KAON MOMENTUM ──\n";
    std::cout << std::setprecision(4);
    std::cout << "  Avg kaon decay time : "
              << std::setprecision(11) << ev.kaon_decay_time << " s\n";
    std::cout << std::setprecision(4);
    std::cout << "  Kaon flight dist    : " << ev.kaon_flight_cm << " cm  (from origin)\n";
    std::cout << "  Kaon velocity       : "
              << std::setprecision(8) << ev.kaon_v << " m/s\n";
    std::cout << std::setprecision(6);
    std::cout << "  Kaon beta           : " << ev.kaon_beta  << "\n";
    std::cout << "  Kaon gamma          : " << ev.kaon_gamma << "\n";
    std::cout << std::setprecision(4);
    std::cout << "  p_reco = gamma*m_K*beta = " << ev.kaon_gamma
              << " * 0.497611 * " << ev.kaon_beta
              << " = " << ev.reco_p << " GeV/c\n";
    std::cout << std::string(78, '=') << "\n\n";
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================
void KLong_visualize_event(
    const char* filename = "Scenario3_Seed1.root",
    int   target_event   = -1,
    const char* output_txt = "")
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

    // Build base name for output files
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

    TRandom3 rng(42);  // fixed seed for reproducible smear in visualization
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
                double st = rng.Gaus(ht, smear_t_sigma);
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
        if (fp.dir.Dot(fm.dir) > 0.9998) continue;  // near-parallel cut

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

        TVector3 kaon_prod(0., 0., 0.);
        double flight = (decay_vtx - kaon_prod).Mag();
        double kaon_v = (flight * 1e-2) / kaon_decay_t;
        double m_K = 0.497611;
        double beta_K = kaon_v / 2.99792458e8;
        if (beta_K >= 1.0) beta_K = 0.9999;
        double gamma_K = 1. / std::sqrt(1. - beta_K * beta_K);
        double kaon_p  = gamma_K * m_K * beta_K;
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
        // Kaon true direction (for display)
        double tp_mag = tp; // already computed
        ve.true_kaon_dir = TVector3(truth_px[en]/tp_mag,
                                    truth_py[en]/tp_mag,
                                    truth_pz[en]/tp_mag);
        // Kinematic intermediates (for print_event_detail)
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
        ve.source_file       = ""; // filled by caller if multi-file
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
    // mark source file (single-file mode: just the filename)
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

    // Find the requested event
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

    // Print step-by-step numerical breakdown to stdout
    print_event_detail(*draw_ev);

    // === Build 4-panel canvas ===
    TCanvas* c = new TCanvas(Form("c_ev%d", target_event),
                             Form("Event %d  |  file: %s", target_event, stem.c_str()),
                             1400, 1000);
    c->SetFillColor(10);

    // Create 4 pads manually for flexible sizing
    // Left column (3D, 60%): top half
    TPad* pad3d = new TPad("pad3d", "3D view",   0.00, 0.50, 0.60, 1.00);
    TPad* padXZ = new TPad("padXZ", "XZ",        0.60, 0.50, 1.00, 1.00);
    TPad* padYZ = new TPad("padYZ", "YZ",        0.00, 0.00, 0.50, 0.50);
    TPad* padXY = new TPad("padXY", "XY",        0.50, 0.00, 1.00, 0.50);

    for (TPad* p : {pad3d, padXZ, padYZ, padXY}) {
        c->cd(); p->Draw();
        p->SetLeftMargin(0.10); p->SetRightMargin(0.02);
        p->SetTopMargin(0.09);  p->SetBottomMargin(0.10);
    }

    draw_3d_panel (pad3d, *draw_ev,
                   z_t1, z_t2, z_t3, z_t4,
                   z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);
    draw_proj_panel(padXZ, *draw_ev, "XZ",
                    z_t1, z_t2, z_t3, z_t4,
                    z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);
    draw_proj_panel(padYZ, *draw_ev, "YZ",
                    z_t1, z_t2, z_t3, z_t4,
                    z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);
    draw_proj_panel(padXY, *draw_ev, "XY",
                    z_t1, z_t2, z_t3, z_t4,
                    z_pizza1, z_pizza2, z_fri1, z_fri2, z_tof);

    // Main title bar
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

// ============================================================================
// MULTI-SEED CONFIG VISUALISER
//
// Merges all raw simulation seed files (*.root, excluding *_vectors*,
// *_acceptance*, *_combined*, *_raw_merged*) in a given directory into a
// single merged file using ROOT's TFileMerger (equivalent to hadd), then
// runs the standard KLong_visualize_event reconstruction on that one file.
//
// The merged file is written to <config_dir>/<stem>_raw_merged.root and is
// REUSED on subsequent calls — the merge step is skipped if the file already
// exists.  Delete it manually to force a re-merge.
//
// Geometry is parsed from the directory name (same tag format as filenames).
//
// USAGE:
//   root -l -q 'KLong_visualize_config("/path/to/dir")'           // text only
//   root -l -q 'KLong_visualize_config("/path/to/dir", 46480)'    // + display event
//   Custom output:
//   root -l -q 'KLong_visualize_config("/path/to/dir", 46480, "out.txt")'
// ============================================================================
void KLong_visualize_config(
    const char* config_dir  = ".",
    int         target_event = -1,
    const char* output_txt  = "")
{
    // --- Collect all raw seed files from the directory ---
    TSystemDirectory tdir("d", config_dir);
    TList* files = tdir.GetListOfFiles();
    if (!files) { std::cout << "Cannot list directory: " << config_dir << "\n"; return; }

    std::vector<std::string> seed_files;
    TIter next(files);
    TObject* obj;
    while ((obj = next())) {
        std::string name(obj->GetName());
        if (name.size() < 5 || name.substr(name.size()-5) != ".root") continue;
        if (name.find("_vectors")     != std::string::npos) continue;
        if (name.find("_acceptance")  != std::string::npos) continue;
        if (name.find("_combined")    != std::string::npos) continue;
        if (name.find("_raw_merged")  != std::string::npos) continue;
        seed_files.push_back(std::string(config_dir) + "/" + name);
    }
    std::sort(seed_files.begin(), seed_files.end());

    if (seed_files.empty()) {
        std::cout << "No raw seed files found in: " << config_dir << "\n";
        std::cout << "(Looking for *.root not matching *_vectors*, *_acceptance*, "
                     "*_combined*, *_raw_merged*)\n";
        return;
    }
    std::cout << "Found " << seed_files.size() << " seed files in " << config_dir << "\n";

    // Directory base name (for filenames and geometry parsing)
    auto dir_basename = [](const std::string& p) -> std::string {
        size_t sl = p.find_last_of("/\\");
        return (sl == std::string::npos) ? p : p.substr(sl+1);
    };
    std::string stem = dir_basename(std::string(config_dir));

    // Merged file lives alongside the seed files
    std::string merged_path = std::string(config_dir) + "/" + stem + "_raw_merged.root";

    // --- Merge seed files with TFileMerger (skip if already done) ---
    bool need_merge = true;
    {
        TFile* probe = TFile::Open(merged_path.c_str(), "READ");
        if (probe && !probe->IsZombie()) { need_merge = false; probe->Close(); delete probe; }
    }
    if (!need_merge) {
        std::cout << "Reusing cached merged file: " << merged_path << "\n";
        std::cout << "  (delete it to force a re-merge)\n";
    } else {
        std::cout << "Merging " << seed_files.size() << " seed files -> "
                  << merged_path << "\n";
        TFileMerger merger(/*isLocal=*/kFALSE);
        merger.OutputFile(merged_path.c_str(), "RECREATE");
        for (const auto& sf : seed_files) merger.AddFile(sf.c_str());
        if (!merger.Merge()) {
            std::cout << "ERROR: TFileMerger failed.\n";
            return;
        }
        std::cout << "Merge complete.\n";
    }

    // Output text file: default to <stem>_all_seeds_delta_p.txt
    std::string txt_path = (std::string(output_txt).empty())
                           ? (stem + "_all_seeds_delta_p.txt") : std::string(output_txt);

    // --- Hand off to the standard single-file visualiser ---
    std::cout << "Running reconstruction on merged file...\n";
    KLong_visualize_event(merged_path.c_str(), target_event, txt_path.c_str());
}
// FIND SIMPLE-GEOMETRY EVENTS
//
// Scans every reconstructed event and evaluates how close each pion track
// lies to a coordinate plane by checking its azimuthal angle phi:
//
//   phi = atan2(dir.Y(), dir.X())  [range: -180 to +180 degrees]
//
//   phi ~   0° or ±180°  =>  track lies in XZ plane  (y(z) ~ 0)
//   phi ~  ±90°          =>  track lies in YZ plane  (x(z) ~ 0)
//
// Events where BOTH pions lie in the same coordinate plane are simplest to
// verify by hand: the 3D path length reduces to sqrt(dH^2 + dz^2) where H
// is either X or Y.
//
// Output: text file sorted by max(|phi deviation|) of both tracks (best first)
//   Columns: event  p_true  p_reco  delta_p  phi+  phi-  theta+  theta-
//            max_dev  plane  vtx_z
//
// USAGE:
//   root -l -q 'KLong_find_simple_events("file.root")'
//   root -l -q 'KLong_find_simple_events("file.root", 10.0, "simple.txt")'
// ============================================================================
void KLong_find_simple_events(
    const char* filename    = "Scenario3_Seed1.root",
    double      phi_tol_deg = 10.0,
    const char* output_txt  = "")
{
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) { std::cout << "Cannot open file\n"; return; }

    // Parse geometry from filename
    std::string fname(filename);
    auto extract = [](const std::string& s, const std::string& key) -> double {
        size_t pos = s.find(key + "-");
        if (pos == std::string::npos) return 0.;
        pos += key.size() + 1;
        size_t end = pos;
        while (end < s.size() && std::isdigit((unsigned char)s[end])) end++;
        return (end == pos) ? 0. : std::stod(s.substr(pos, end-pos));
    };
    double z_tof    = extract(fname, "E1");
    double z_pizza1 = extract(fname, "P1");
    double z_pizza2 = extract(fname, "P2");
    (void)z_pizza1; (void)z_pizza2;

    // Select pi+pi-pi0 events
    TTree* tree2 = (TTree*)file->Get("Ntuple2");
    if (!tree2) { file->Close(); return; }
    Double_t evtNb2, dkPDG;
    tree2->SetBranchAddress("evtNb", &evtNb2);
    tree2->SetBranchAddress("DKparticle_PDGEncoding", &dkPDG);
    std::map<int, std::vector<int>> evt_products;
    for (Long64_t i=0; i<tree2->GetEntries(); i++) {
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

    // MC truth
    std::unordered_map<int,double> truth_px, truth_py, truth_pz;
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
        for (Long64_t i=0; i<tree1->GetEntries(); i++) {
            tree1->GetEntry(i);
            int en=(int)en1;
            truth_px[en]=kpx; truth_py[en]=kpy; truth_pz[en]=kpz;
        }
    }

    // Load hits
    TTree* tree3 = (TTree*)file->Get("Ntuple3");
    if (!tree3) { file->Close(); return; }
    Double_t evtNb3, hx, hy, hz, ht, pdg3, devID3;
    tree3->SetBranchAddress("evtNb",       &evtNb3);
    tree3->SetBranchAddress("hitX",        &hx);
    tree3->SetBranchAddress("hitY",        &hy);
    tree3->SetBranchAddress("hitZ",        &hz);
    tree3->SetBranchAddress("hitT",        &ht);
    tree3->SetBranchAddress("PDGEncoding", &pdg3);
    tree3->SetBranchAddress("deviceID",    &devID3);

    auto is_tracker = [](int id){ return id>=1    && id<=1952; };
    auto is_fri     = [](int id){ return id>=2001 && id<=2052; };
    auto is_pizza   = [](int id){ return id>=1953 && id<=2000; };
    auto is_tof     = [](int id){ return id>=2053 && id<=2070; };

    TRandom3 rng(42);
    const double smear_t = 0.0015;
    std::unordered_map<int, EventReco> ev_data;
    for (Long64_t i=0; i<tree3->GetEntries(); i++) {
        tree3->GetEntry(i);
        int en=(int)evtNb3;
        if (sel_set.find(en)==sel_set.end()) continue;
        auto& ev=ev_data[en];
        int id=(int)devID3;
        for (int pion=0; pion<2; pion++) {
            if ((int)pdg3 != (pion==0 ? 211 : -211)) continue;
            if (is_tracker(id)||is_fri(id)) {
                (pion==0 ? ev.hits_pip : ev.hits_pim).push_back({hx,hy,hz,ht,id});
            }
            if (is_pizza(id)) {
                double st=rng.Gaus(ht,smear_t);
                if (pion==0) { if (!ev.has_pip_pizza||st<ev.pip_pizza_time) {
                    ev.has_pip_pizza=true; ev.pip_pizza_time=st;
                    ev.pip_pizza_x=hx; ev.pip_pizza_y=hy; ev.pip_pizza_z=hz; }
                } else     { if (!ev.has_pim_pizza||st<ev.pim_pizza_time) {
                    ev.has_pim_pizza=true; ev.pim_pizza_time=st;
                    ev.pim_pizza_x=hx; ev.pim_pizza_y=hy; ev.pim_pizza_z=hz; } }
            }
            if (is_tof(id)) {
                double st=rng.Gaus(ht,smear_t);
                double smx=rng.Gaus(tof_bar_xc(id),TOF_BAR_HALF_WIDTH);
                double smy=rng.Gaus(tof_bar_yc(id),TOF_Y_UNCERTAINTY);
                if (pion==0) { if (!ev.has_pip_tof||st<ev.pip_tof_time) {
                    ev.has_pip_tof=true; ev.pip_tof_time=st; ev.pip_tof_deviceID=id;
                    ev.pip_tof_x=smx; ev.pip_tof_y=smy; ev.pip_tof_z=z_tof; }
                } else     { if (!ev.has_pim_tof||st<ev.pim_tof_time) {
                    ev.has_pim_tof=true; ev.pim_tof_time=st; ev.pim_tof_deviceID=id;
                    ev.pim_tof_x=smx; ev.pim_tof_y=smy; ev.pim_tof_z=z_tof; } }
            }
        }
    }

    // Minimum angular distance to any of {0, 90, 180, 270} degrees (mod 360)
    auto phi_plane_dev = [](double phi_deg) -> double {
        double p = std::fmod(std::fabs(phi_deg), 180.);  // fold to [0,180)
        // distance to nearest of 0° or 90°
        double d0 = std::min(p, 180.-p);      // distance to 0° or 180°
        double d90 = std::fabs(p - 90.);      // distance to 90° (equiv -90°)
        return std::min(d0, d90);
    };
    // Which plane: "XZ" if nearest to 0/180, "YZ" if nearest to +-90
    auto plane_name = [](double phi_deg) -> const char* {
        double p = std::fmod(std::fabs(phi_deg), 180.);
        double d0  = std::min(p, 180.-p);
        double d90 = std::fabs(p - 90.);
        return (d0 <= d90) ? "XZ" : "YZ";
    };

    // Result struct
    struct SimpleResult {
        int    event_num;
        double p_true, p_reco, delta_p;
        double phi_pip, phi_pim, theta_pip, theta_pim;
        double max_dev;   // max(phi_plane_dev(phi+), phi_plane_dev(phi-))
        std::string plane;
        double vtx_z;
    };
    std::vector<SimpleResult> results;

    for (int en : sel_events) {
        auto ev_it = ev_data.find(en);
        if (ev_it==ev_data.end()) continue;
        const auto& ev = ev_it->second;
        if (!ev.has_pip_pizza||!ev.has_pip_tof||!ev.has_pim_pizza||!ev.has_pim_tof) continue;

        TrackFitV fp = fit_track_vis(ev.hits_pip, ev.pip_tof_deviceID, z_tof);
        TrackFitV fm = fit_track_vis(ev.hits_pim, ev.pim_tof_deviceID, z_tof);
        if (!fp.valid||!fm.valid) continue;

        TVector3 decay_vtx = poca_vis(fp.orig, fp.dir, fm.orig, fm.dir);
        if (fp.dir.Dot(fm.dir) > 0.9998) continue;

        TVector3 pip_pza = eval_at_z(fp, ev.pip_pizza_z);
        TVector3 pim_pza = eval_at_z(fm, ev.pim_pizza_z);
        TVector3 pip_tof(ev.pip_tof_x, ev.pip_tof_y, ev.pip_tof_z);
        TVector3 pim_tof(ev.pim_tof_x, ev.pim_tof_y, ev.pim_tof_z);
        double pipL=(pip_tof-pip_pza).Mag(), pimL=(pim_tof-pim_pza).Mag();
        double pipdt=ev.pip_tof_time-ev.pip_pizza_time;
        double pimdt=ev.pim_tof_time-ev.pim_pizza_time;
        if (pipdt<=0||pimdt<=0) continue;
        double pipv=(pipL*1e-2)/(pipdt*1e-9), pimv=(pimL*1e-2)/(pimdt*1e-9);
        double pipdp=(pip_pza-decay_vtx).Mag(), pimdp=(pim_pza-decay_vtx).Mag();
        double piptd=ev.pip_pizza_time*1e-9-(pipdp*1e-2)/pipv;
        double pimtd=ev.pim_pizza_time*1e-9-(pimdp*1e-2)/pimv;
        double kdt=0.5*(piptd+pimtd);
        double fl=(decay_vtx).Mag();
        double kv=(fl*1e-2)/kdt;
        double bk=kv/2.99792458e8; if(bk>=1.)bk=0.9999;
        double gk=1./std::sqrt(1.-bk*bk);
        double kp=gk*0.497611*bk;
        if (kp<=0||kp>11.) continue;
        if (truth_px.find(en)==truth_px.end()) continue;
        double tp=std::sqrt(truth_px[en]*truth_px[en]+truth_py[en]*truth_py[en]+truth_pz[en]*truth_pz[en]);
        if (tp==0.) continue;

        double phi_p   = std::atan2(fp.dir.Y(), fp.dir.X()) * 180./M_PI;
        double phi_m   = std::atan2(fm.dir.Y(), fm.dir.X()) * 180./M_PI;
        double theta_p = std::acos(std::min(1.0, std::fabs(fp.dir.Z()))) * 180./M_PI;
        double theta_m = std::acos(std::min(1.0, std::fabs(fm.dir.Z()))) * 180./M_PI;
        double devp    = phi_plane_dev(phi_p);
        double devm    = phi_plane_dev(phi_m);
        double maxdev  = std::max(devp, devm);

        // Label plane: both pions must point toward same family of planes
        std::string plabel = "MIXED";
        if (std::string(plane_name(phi_p)) == plane_name(phi_m)) {
            plabel = plane_name(phi_p);
        }

        results.push_back({en, tp, kp, kp-tp, phi_p, phi_m, theta_p, theta_m,
                           maxdev, plabel, decay_vtx.Z()});
    }

    std::sort(results.begin(), results.end(),
              [](const SimpleResult& a, const SimpleResult& b){ return a.max_dev < b.max_dev; });

    // Determine output filename
    auto base_name = [](const std::string& p) -> std::string {
        size_t sl=p.find_last_of("/\\");
        std::string fn=(sl==std::string::npos)?p:p.substr(sl+1);
        size_t dot=fn.find_last_of('.'); return (dot==std::string::npos)?fn:fn.substr(0,dot);
    };
    std::string stem = base_name(fname);
    std::string out_path = (std::string(output_txt).empty())
                           ? (stem + "_simple_geometry.txt") : std::string(output_txt);

    std::ofstream ofs(out_path);
    ofs << std::fixed << std::setprecision(3);
    ofs << "# Simple-geometry events — sorted by max phi-plane deviation (smallest = most aligned)\n";
    ofs << "# File: " << filename << "   phi_tol = " << phi_tol_deg << " deg\n";
    ofs << "# Plane: XZ = phi~0/180 (tracks in XZ plane), YZ = phi~+-90 (tracks in YZ plane)\n";
    ofs << "# Pass event number to KLong_visualize_event to see a diagram;\n";
    ofs << "#   pass it to print_event_detail (via visualize) for the numeric breakdown.\n";
    ofs << "#\n";
    ofs << "# " << std::setw(8)  << "event"
                << std::setw(8)  << "p_true"
                << std::setw(8)  << "p_reco"
                << std::setw(8)  << "delta_p"
                << std::setw(8)  << "phi+"
                << std::setw(8)  << "phi-"
                << std::setw(8)  << "theta+"
                << std::setw(8)  << "theta-"
                << std::setw(9)  << "max_dev"
                << std::setw(7)  << "plane"
                << std::setw(9)  << "vtx_z"
                << "\n";
    ofs << "# " << std::string(95, '-') << "\n";

    int n_xz=0, n_yz=0, n_in_tol=0;
    for (const auto& r : results) {
        bool in_tol = (r.max_dev <= phi_tol_deg);
        if (in_tol) {
            n_in_tol++;
            if (r.plane=="XZ") n_xz++;
            if (r.plane=="YZ") n_yz++;
        }
        ofs << "  " << std::setw(8)  << r.event_num
                    << std::setw(8)  << r.p_true
                    << std::setw(8)  << r.p_reco
                    << std::setw(8)  << r.delta_p
                    << std::setw(8)  << r.phi_pip
                    << std::setw(8)  << r.phi_pim
                    << std::setw(8)  << r.theta_pip
                    << std::setw(8)  << r.theta_pim
                    << std::setw(9)  << r.max_dev
                    << std::setw(7)  << r.plane
                    << std::setw(9)  << r.vtx_z
                    << (in_tol ? "  *" : "")
                    << "\n";
    }
    ofs.close();

    std::cout << "Reconstructed:   " << results.size() << " events\n";
    std::cout << "Within tol (" << phi_tol_deg << " deg):  " << n_in_tol
              << "  (XZ=" << n_xz << "  YZ=" << n_yz << ")\n";
    std::cout << "Wrote: " << out_path << "\n";
    if (n_in_tol > 0) {
        std::cout << "Top candidate: event " << results[0].event_num
                  << "  max_dev=" << results[0].max_dev << " deg"
                  << "  phi+=" << results[0].phi_pip
                  << "  phi-=" << results[0].phi_pim
                  << "  plane=" << results[0].plane
                  << "  p_true=" << results[0].p_true << " GeV/c\n";
        std::cout << "  -> visualise with: "
                     "KLong_visualize_event(\"" << filename << "\", "
                  << results[0].event_num << ")\n";
    }
    file->Close();
}

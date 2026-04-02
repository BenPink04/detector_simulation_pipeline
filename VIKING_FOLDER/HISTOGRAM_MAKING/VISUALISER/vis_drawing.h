// ============================================================================
// vis_drawing.h — all drawing routines and the detailed event printout for
//                the KLong visualiser suite.
//
// Includes vis_common.h (do not include both manually in the same translation
// unit — just include this file and you get everything).
// ============================================================================
#ifndef VIS_DRAWING_H
#define VIS_DRAWING_H

#include "vis_common.h"

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

// Draw a circular ring at fixed z (outer radius only)
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
        double bar_half_y = (id == 2061 || id == 2062) ? 2.0 : 60.0;
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
    pad->SetFillColor(kBlack);

    double margin = 20.;
    double xext = 110., yext = 80.;
    double rmin[3] = { -xext, -yext, -30. };
    double rmax[3] = {  xext,  yext,  z_tof + margin };
    TView* view = TView::CreateView(1, rmin, rmax);
    Int_t irep;
    view->SetView(30., 75., 0., irep);
    view->ShowAxis();

    // ---- Detector planes ----
    const double TRK = STRAW_HALF_LENGTH;
    int col_trk = kGray+1;
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t1, col_trk, 1);
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t2, col_trk, 1);
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t3, col_trk, 1);
    draw_rect3d(-TRK, TRK, -TRK, TRK, z_t4, col_trk, 1);

    draw_ring3d(z_pizza1, 20., kOrange+1);
    draw_ring3d(z_pizza2, 20., kOrange+1);

    draw_rect3d(-66., 66., -70.25, 70.25, z_fri1, kCyan+1, 1);
    draw_rect3d(-66., 66., -70.25, 70.25, z_fri2, kCyan+1, 1);

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

    // ---- pi+ / pi- track lines ----
    if (ev.fit_pip.valid)
        draw_track3d(ev.fit_pip, z_pizza1 - 5., z_tof + 5., kRed,    1, 2);
    if (ev.fit_pim.valid)
        draw_track3d(ev.fit_pim, z_pizza1 - 5., z_tof + 5., kBlue+1, 1, 2);

    // ---- Tracker + FRI hits ----
    if (!ev.hits_pip.empty()) {
        TPolyMarker3D* m = new TPolyMarker3D((int)ev.hits_pip.size(), 2);
        for (int i = 0; i < (int)ev.hits_pip.size(); i++)
            m->SetPoint(i, ev.hits_pip[i].x, ev.hits_pip[i].y, ev.hits_pip[i].z);
        m->SetMarkerColor(kRed); m->SetMarkerStyle(2); m->SetMarkerSize(1.2); m->Draw();
    }
    if (!ev.hits_pim.empty()) {
        TPolyMarker3D* m = new TPolyMarker3D((int)ev.hits_pim.size(), 2);
        for (int i = 0; i < (int)ev.hits_pim.size(); i++)
            m->SetPoint(i, ev.hits_pim[i].x, ev.hits_pim[i].y, ev.hits_pim[i].z);
        m->SetMarkerColor(kBlue+1); m->SetMarkerStyle(2); m->SetMarkerSize(1.2); m->Draw();
    }

    // ---- PIZZA hit positions ----
    {
        TPolyMarker3D* mpz = new TPolyMarker3D(1, 29);
        mpz->SetPoint(0, ev.pip_pizza_raw.X(), ev.pip_pizza_raw.Y(), ev.pip_pizza_raw.Z());
        mpz->SetMarkerColor(kRed); mpz->SetMarkerSize(2.); mpz->Draw();
    }
    {
        TPolyMarker3D* mpz = new TPolyMarker3D(1, 29);
        mpz->SetPoint(0, ev.pim_pizza_raw.X(), ev.pim_pizza_raw.Y(), ev.pim_pizza_raw.Z());
        mpz->SetMarkerColor(kBlue+1); mpz->SetMarkerSize(2.); mpz->Draw();
    }

    // ---- TOF smeared positions ----
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

    // ---- Vertices ----
    {
        TPolyMarker3D* mv = new TPolyMarker3D(1, 20);
        mv->SetPoint(0, ev.reco_vtx.X(), ev.reco_vtx.Y(), ev.reco_vtx.Z());
        mv->SetMarkerColor(kGreen+2); mv->SetMarkerSize(3.); mv->Draw();
    }
    {
        TPolyMarker3D* mv = new TPolyMarker3D(1, 24);
        mv->SetPoint(0, ev.true_vtx.X(), ev.true_vtx.Y(), ev.true_vtx.Z());
        mv->SetMarkerColor(kWhite); mv->SetMarkerSize(3.); mv->Draw();
    }

    // Legend
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
    pad->SetFillColor(10);
    pad->SetLeftMargin(0.12); pad->SetBottomMargin(0.12);

    bool is_xz = (std::string(proj) == "XZ");
    bool is_yz = (std::string(proj) == "YZ");
    bool is_xy = (std::string(proj) == "XY");

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

    auto hv = [&](const TVector3& p) -> std::pair<double,double> {
        if (is_xz) return {p.Z(), p.X()};
        if (is_yz) return {p.Z(), p.Y()};
        return {p.X(), p.Y()};
    };

    // ---- Detector z-lines (XZ and YZ only) ----
    if (!is_xy) {
        auto draw_det_line = [&](double z, int color, int style) {
            TLine* l = new TLine(z, vmin, z, vmax);
            l->SetLineColor(color); l->SetLineStyle(style); l->SetLineWidth(1);
            l->Draw();
        };
        draw_det_line(z_pizza1, kOrange+1, 2);
        draw_det_line(z_pizza2, kOrange+1, 2);
        draw_det_line(z_t1,   kGray+1,    1);
        draw_det_line(z_t2,   kGray+1,    1);
        draw_det_line(z_fri1, kCyan+1,    1);
        draw_det_line(z_fri2, kCyan+1,    1);
        draw_det_line(z_t3,   kGray+1,    1);
        draw_det_line(z_t4,   kGray+1,    1);
        draw_det_line(z_tof,  kMagenta+1, 1);
        TLatex tl; tl.SetTextSize(0.038); tl.SetTextAngle(90); tl.SetTextAlign(12);
        tl.SetTextColor(kGray+2);
        tl.DrawLatex(z_t1,   vmax*0.6, "T1");
        tl.DrawLatex(z_t2,   vmax*0.6, "T2");
        tl.DrawLatex(z_fri1, vmax*0.6, "F1");
        tl.DrawLatex(z_fri2, vmax*0.6, "F2");
        tl.DrawLatex(z_t3,   vmax*0.6, "T3");
        tl.DrawLatex(z_t4,   vmax*0.6, "T4");
        tl.SetTextColor(kMagenta+1);
        tl.DrawLatex(z_tof, vmax*0.6, "TOF");
        tl.SetTextColor(kOrange+1);
        tl.DrawLatex(z_pizza1, vmax*0.6, "P1");
        tl.DrawLatex(z_pizza2, vmax*0.6, "P2");
    }
    if (is_xy) {
        TEllipse* e = new TEllipse(0., 0., 20., 20.); e->SetFillStyle(0);
        e->SetLineColor(kOrange+1); e->SetLineWidth(2); e->Draw("same");
        TBox* bt = new TBox(-45., -45., 45., 45.);
        bt->SetFillStyle(0); bt->SetLineColor(kGray+1); bt->Draw("l same");
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
    TGraph* gp = make_hit_graph(ev.hits_pip, kRed,    2);
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

    TLatex ltx;
    ltx.SetNDC(true); ltx.SetTextSize(0.052); ltx.SetTextFont(62);
    ltx.DrawLatex(0.15, 0.94, title);
}

// ============================================================================
// DETAILED EVENT BREAKDOWN
//
// Prints a step-by-step numerical breakdown to stdout: all detector hits,
// track fit parameters, theta/phi of each track, PIZZA/TOF positions, path
// lengths, velocities, decay times, and the full kaon momentum derivation.
// Call this to verify the reconstruction by hand.
// ============================================================================
void print_event_detail(const VisEvent& ev)
{
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

        const TVector3& piz_raw = pion==0 ? ev.pip_pizza_raw  : ev.pim_pizza_raw;
        const TVector3& piz_pos = pion==0 ? ev.pip_pizza_pos  : ev.pim_pizza_pos;
        double piz_t  = pion==0 ? ev.pip_pizza_time : ev.pim_pizza_time;
        std::cout << "  PIZZA raw hit      : (" << piz_raw.X() << ", " << piz_raw.Y()
                  << ", " << piz_raw.Z() << ") cm   t_smeared = " << piz_t << " ns\n";
        std::cout << "  PIZZA track extrap : (" << piz_pos.X() << ", " << piz_pos.Y()
                  << ", " << piz_pos.Z() << ") cm   [used for path & decay time]\n";

        const TVector3& tof_pos = pion==0 ? ev.pip_tof_pos  : ev.pim_tof_pos;
        int   tof_dev = pion==0 ? ev.pip_tof_devID   : ev.pim_tof_devID;
        double tof_t  = pion==0 ? ev.pip_tof_time    : ev.pim_tof_time;
        std::cout << "  TOF bar devID      : " << tof_dev
                  << "   bar_centre_x = " << tof_bar_xc(tof_dev)
                  << " cm   bar_centre_y = " << tof_bar_yc(tof_dev) << " cm\n";
        std::cout << "  TOF smeared pos    : (" << tof_pos.X() << ", " << tof_pos.Y()
                  << ", " << tof_pos.Z() << ") cm   t_smeared = " << tof_t << " ns\n";

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

#endif // VIS_DRAWING_H

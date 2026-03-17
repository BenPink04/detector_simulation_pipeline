// KLong_validate_vertex.C
//
// Compares the reconstructed decay vertex to the true (Geant4) decay vertex.
// Reads the kaonVectors TTree from the TGRAPH_TEST archive and produces:
//
//   Canvas 1 — KLong_validate_vertex_residuals.png
//     2×2 grid of residual distributions per configuration:
//       Δx = reco_x − true_x
//       Δy = reco_y − true_y
//       Δz = reco_z − true_z
//       |Δr| = sqrt(Δx²+Δy²+Δz²)   (3D distance)
//     Each panel shows all configurations overlaid with a Gaussian fit.
//
//   Canvas 2 — KLong_validate_vertex_z_profile.png
//     Mean reconstructed z vs true z (TProfile) per configuration.
//     Ideal response = diagonal. Deviation shows systematic z-position bias.
//
// Usage:
//   root -l KLong_validate_vertex.C

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cmath>

// ── Detector layout helpers ───────────────────────────────────────────────────

struct DetectorLayout {
    std::vector<double> trackers, pizzas, fris;
    double end_detector = 0;
};

DetectorLayout parseConfigString(const std::string& config) {
    DetectorLayout layout;
    std::istringstream ss(config);
    std::string token;
    while (std::getline(ss, token, '_')) {
        if (token.empty()) continue;
        size_t pos = token.find('-');
        if (pos == std::string::npos) continue;
        std::string dtype = token.substr(0, pos);
        double position   = std::stod(token.substr(pos + 1));
        if (position == 0) continue;
        if      (dtype[0] == 'T') layout.trackers.push_back(position);
        else if (dtype[0] == 'P') layout.pizzas.push_back(position);
        else if (dtype[0] == 'F') layout.fris.push_back(position);
        else if (dtype[0] == 'E') layout.end_detector = position;
    }
    return layout;
}

void drawDetectorLayout(const std::vector<std::string>& config_labels,
                        const std::vector<int>& line_colors) {
    double min_pos = 1e9, max_pos = -1e9;
    for (const auto& cfg : config_labels) {
        DetectorLayout l = parseConfigString(cfg);
        for (double p : l.trackers) { min_pos = std::min(min_pos,p); max_pos = std::max(max_pos,p); }
        for (double p : l.pizzas)   { min_pos = std::min(min_pos,p); max_pos = std::max(max_pos,p); }
        for (double p : l.fris)     { min_pos = std::min(min_pos,p); max_pos = std::max(max_pos,p); }
        min_pos = std::min(min_pos, l.end_detector);
        max_pos = std::max(max_pos, l.end_detector);
    }
    double range = max_pos - min_pos;
    min_pos -= range * 0.1;
    max_pos += range * 0.1;

    TText *title = new TText(0.5, 0.97, "Detector Layouts");
    title->SetTextAlign(23); title->SetTextSize(0.04); title->Draw();

    int n = (int)config_labels.size();
    double y_spacing = 0.80 / (n + 1);

    for (int i = 0; i < n; ++i) {
        DetectorLayout l = parseConfigString(config_labels[i]);
        double yc = 0.88 - (i + 1) * y_spacing;

        TLine *ind = new TLine(0.02, yc, 0.08, yc);
        ind->SetLineColor(line_colors[i]); ind->SetLineWidth(4); ind->Draw();

        const double xs = 0.80, xo = 0.15;
        auto xp = [&](double pos){ return xo + (pos - min_pos) / (max_pos - min_pos) * xs; };

        for (double p : l.trackers) {
            TBox *b = new TBox(xp(p)-0.008, yc-0.015, xp(p)+0.008, yc+0.015);
            b->SetFillColor(kBlue); b->SetLineColor(kBlue); b->Draw(); }
        for (double p : l.pizzas) {
            TBox *b = new TBox(xp(p)-0.008, yc-0.015, xp(p)+0.008, yc+0.015);
            b->SetFillColor(kRed); b->SetLineColor(kRed); b->Draw(); }
        for (double p : l.fris) {
            TBox *b = new TBox(xp(p)-0.008, yc-0.015, xp(p)+0.008, yc+0.015);
            b->SetFillColor(kMagenta); b->SetLineColor(kMagenta); b->Draw(); }
        TBox *be = new TBox(xp(l.end_detector)-0.008, yc-0.015,
                            xp(l.end_detector)+0.008, yc+0.015);
        be->SetFillColor(kGreen); be->SetLineColor(kGreen); be->Draw();
    }

    double ly = 0.06, lx = 0.03, lsp = 0.22;
    auto addLegItem = [&](double x, int col, const char *lbl) {
        TBox *b = new TBox(x, ly, x+0.025, ly+0.035); b->SetFillColor(col); b->Draw();
        TText *t = new TText(x+0.035, ly+0.018, lbl); t->SetTextSize(0.030); t->Draw();
    };
    addLegItem(lx,          kBlue,    "Tracker");
    addLegItem(lx+lsp,      kRed,     "Pizza");
    addLegItem(lx+2*lsp,    kMagenta, "FRI");
    addLegItem(lx+3*lsp,    kGreen,   "End");
}

// ── Main macro ────────────────────────────────────────────────────────────────

void KLong_validate_vertex() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // ── Discover files ─────────────────────────────────────────────────────
    std::string findCmd =
        "find /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/TGRAPH_TEST_20260316"
        " -name '*combined_vectors.root' 2>/dev/null | sort";
    FILE *pipe = popen(findCmd.c_str(), "r");
    if (!pipe) { std::cerr << "Error: popen failed\n"; return; }

    std::vector<std::string> root_files, config_labels;
    char buf[512];
    while (fgets(buf, sizeof(buf), pipe)) {
        std::string fp(buf);
        fp.erase(std::remove(fp.begin(), fp.end(), '\n'), fp.end());
        if (fp.empty()) continue;
        root_files.push_back(fp);
        size_t sl  = fp.find_last_of('/');
        std::string fn = fp.substr(sl + 1);
        size_t sfx = fn.find("_combined_vectors.root");
        config_labels.push_back((sfx != std::string::npos) ? fn.substr(0, sfx) : fn);
    }
    pclose(pipe);

    std::cout << "Found " << root_files.size() << " configuration(s)\n";
    if (root_files.empty()) { std::cerr << "No files found — check archive path\n"; return; }

    int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kCyan+1, kOrange+1, kViolet, kTeal-5};
    std::vector<int> line_colors;

    // ── Per-config storage ─────────────────────────────────────────────────
    struct ConfigResult {
        std::string label;
        int color;
        TH1D    *h_dx = nullptr, *h_dy = nullptr, *h_dz = nullptr, *h_dr = nullptr;
        TProfile *prof_z = nullptr;   // mean reco_z per true_z bin
        // Gaussian fit results
        double mu_dx=0, sig_dx=0, mu_dy=0, sig_dy=0, mu_dz=0, sig_dz=0;
        double mean_dr=0, sig_dr=0;
    };
    std::vector<ConfigResult> results;

    // ── Loop over configurations ───────────────────────────────────────────
    for (size_t idx = 0; idx < root_files.size(); ++idx) {
        ConfigResult res;
        res.label = config_labels[idx];
        res.color = colors[idx % 8];
        line_colors.push_back(res.color);

        TFile *file = TFile::Open(root_files[idx].c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening " << root_files[idx] << "\n"; continue; }

        TTree *tree = (TTree*)file->Get("kaonVectors");
        if (!tree) {
            std::cerr << "kaonVectors tree not found\n"; file->Close(); continue; }

        std::vector<double> *reco_x=nullptr, *reco_y=nullptr, *reco_z=nullptr;
        std::vector<double> *true_x=nullptr, *true_y=nullptr, *true_z=nullptr;
        tree->SetBranchAddress("reco_vertex_x", &reco_x);
        tree->SetBranchAddress("reco_vertex_y", &reco_y);
        tree->SetBranchAddress("reco_vertex_z", &reco_z);
        tree->SetBranchAddress("true_vertex_x", &true_x);
        tree->SetBranchAddress("true_vertex_y", &true_y);
        tree->SetBranchAddress("true_vertex_z", &true_z);

        // Histograms — wide ranges to capture outlier tails
        std::string n = Form("%zu", idx);
        res.h_dx   = new TH1D(("h_dx_"+n).c_str(),   "", 120, -30.0,  30.0);   // cm
        res.h_dy   = new TH1D(("h_dy_"+n).c_str(),   "", 120, -30.0,  30.0);   // cm
        res.h_dz   = new TH1D(("h_dz_"+n).c_str(),   "", 120, -300.0, 300.0);  // cm
        res.h_dr   = new TH1D(("h_dr_"+n).c_str(),   "", 100,   0.0,  300.0);  // cm
        // Profile: true_z  0 → 250 cm, 25 bins (10 cm each)
        res.prof_z = new TProfile(("prof_z_"+n).c_str(), "", 25, 0.0, 250.0);

        for (auto *h : {res.h_dx, res.h_dy, res.h_dz, res.h_dr})
            h->SetDirectory(nullptr);
        res.prof_z->SetDirectory(nullptr);

        Long64_t nEnt = tree->GetEntries();
        for (Long64_t e = 0; e < nEnt; ++e) {
            tree->GetEntry(e);
            for (size_t i = 0; i < reco_z->size(); ++i) {
                double rx = (*reco_x)[i], ry = (*reco_y)[i], rz = (*reco_z)[i];
                double tx = (*true_x)[i], ty = (*true_y)[i], tz = (*true_z)[i];

                if (!std::isfinite(rx)||!std::isfinite(ry)||!std::isfinite(rz)) continue;
                if (!std::isfinite(tx)||!std::isfinite(ty)||!std::isfinite(tz)) continue;

                double dx = rx - tx, dy = ry - ty, dz = rz - tz;
                double dr = std::sqrt(dx*dx + dy*dy + dz*dz);

                res.h_dx->Fill(dx);
                res.h_dy->Fill(dy);
                res.h_dz->Fill(dz);
                res.h_dr->Fill(dr);
                res.prof_z->Fill(tz, rz);
            }
        }
        file->Close();

        // Style
        for (auto *h : {res.h_dx, res.h_dy, res.h_dz, res.h_dr}) {
            h->SetLineColor(res.color); h->SetLineWidth(2); }
        res.prof_z->SetLineColor(res.color);
        res.prof_z->SetLineWidth(2);
        res.prof_z->SetMarkerColor(res.color);
        res.prof_z->SetMarkerStyle(20 + (int)idx);
        res.prof_z->SetMarkerSize(0.9);

        // Gaussian fits: fit core (±2 RMS) for each residual
        auto fitGaus = [](TH1D *h, double &mu, double &sig) {
            if (h->GetEntries() < 50) return;
            double m = h->GetMean(), r = h->GetRMS();
            TF1 *g = new TF1("gfit_tmp", "gaus", m - 2*r, m + 2*r);
            g->SetParameters(h->GetMaximum(), m, r);
            h->Fit(g, "RQN");
            mu = g->GetParameter(1); sig = g->GetParameter(2);
            delete g;
        };
        fitGaus(res.h_dx, res.mu_dx, res.sig_dx);
        fitGaus(res.h_dy, res.mu_dy, res.sig_dy);
        fitGaus(res.h_dz, res.mu_dz, res.sig_dz);
        // For |Δr| use histogram mean and RMS (not a Gaussian — 3D distance is strictly positive)
        if (res.h_dr->GetEntries() > 0) {
            res.mean_dr = res.h_dr->GetMean();
            res.sig_dr  = res.h_dr->GetRMS();
        }

        results.push_back(res);
    }

    if (results.empty()) { std::cerr << "No valid results\n"; return; }

    // ─────────────────────────────────────────────────────────────────────────
    // CANVAS 1 — Residual histograms (2×2 grid + detector layout column)
    // ─────────────────────────────────────────────────────────────────────────
    {
        TCanvas *c1 = new TCanvas("c_vtx", "Vertex Residuals", 1600, 900);

        // Outer pads: left 75% = 2×2 residuals, right 25% = detector layout
        TPad *p_plots  = new TPad("p_plots",  "", 0.00, 0.00, 0.75, 1.00);
        TPad *p_layout = new TPad("p_layout", "", 0.75, 0.00, 1.00, 1.00);
        p_plots->Draw(); p_layout->Draw();

        p_plots->Divide(2, 2, 0.005, 0.005);

        // Helper: draw one residual panel
        auto drawPanel = [&](int padN, TH1D* ConfigResult::*hPtr,
                              const char *title, const char *xtitle,
                              double xlo, double xhi,
                              double (ConfigResult::*mu_ptr),
                              double (ConfigResult::*sig_ptr)) {
            p_plots->cd(padN);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
            gPad->SetGrid(1, 0);

            double hmax = 0;
            for (const auto &r : results)
                hmax = std::max(hmax, (r.*hPtr)->GetMaximum());

            bool first = true;
            for (const auto &r : results) {
                TH1D *h = r.*hPtr;
                h->GetXaxis()->SetRangeUser(xlo, xhi);
                h->GetYaxis()->SetRangeUser(0, hmax * 1.25);
                h->SetTitle(Form("%s;%s;Entries", title, xtitle));
                h->GetXaxis()->SetTitleSize(0.05);
                h->GetYaxis()->SetTitleSize(0.05);
                h->Draw(first ? "HIST" : "HIST SAME");
                first = false;
            }

            // Reference line at 0
            TLine *zero = new TLine(0, 0, 0, hmax * 1.25);
            zero->SetLineColor(kBlack); zero->SetLineStyle(2); zero->SetLineWidth(2);
            zero->Draw("SAME");
        };

        // Δx
        drawPanel(1, &ConfigResult::h_dx,
                  "#Deltax = x_{reco} #minus x_{true}", "#Deltax [cm]",
                  -30, 30, &ConfigResult::mu_dx, &ConfigResult::sig_dx);
        // Δy
        drawPanel(2, &ConfigResult::h_dy,
                  "#Deltay = y_{reco} #minus y_{true}", "#Deltay [cm]",
                  -30, 30, &ConfigResult::mu_dy, &ConfigResult::sig_dy);
        // Δz
        drawPanel(3, &ConfigResult::h_dz,
                  "#Deltaz = z_{reco} #minus z_{true}", "#Deltaz [cm]",
                  -300, 300, &ConfigResult::mu_dz, &ConfigResult::sig_dz);
        // |Δr| — no zero line (starts at 0)
        {
            p_plots->cd(4);
            gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
            gPad->SetGrid(1, 0);

            double hmax = 0;
            for (const auto &r : results)
                hmax = std::max(hmax, r.h_dr->GetMaximum());

            bool first = true;
            for (const auto &r : results) {
                r.h_dr->GetXaxis()->SetRangeUser(0, 300);
                r.h_dr->GetYaxis()->SetRangeUser(0, hmax * 1.25);
                r.h_dr->SetTitle(
                    "|#Deltar| = #sqrt{#Deltax^{2}+#Deltay^{2}+#Deltaz^{2}};|#Deltar| [cm];Entries");
                r.h_dr->GetXaxis()->SetTitleSize(0.05);
                r.h_dr->GetYaxis()->SetTitleSize(0.05);
                r.h_dr->Draw(first ? "HIST" : "HIST SAME");
                first = false;
            }
        }

        // Detector layout panel
        p_layout->cd();
        drawDetectorLayout(config_labels, line_colors);

        c1->Update();
        c1->SaveAs("KLong_validate_vertex_residuals.png");
        std::cout << "\n=== Saved KLong_validate_vertex_residuals.png ===" << std::endl;
        delete c1;
    }

    // ─────────────────────────────────────────────────────────────────────────
    // CANVAS 2 — reco_z vs true_z profile
    // ─────────────────────────────────────────────────────────────────────────
    {
        TCanvas *c2 = new TCanvas("c_prof", "Reco vertex z vs True vertex z", 1400, 700);
        c2->Divide(2, 1);

        c2->cd(1);
        gPad->SetPad(0.00, 0.00, 0.75, 1.00);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        gPad->SetGrid(1, 1);

        // Find y range from profiles
        double ymin = 1e9, ymax = -1e9;
        for (const auto &r : results) {
            for (int b = 1; b <= r.prof_z->GetNbinsX(); ++b) {
                if (r.prof_z->GetBinEntries(b) < 5) continue;
                double v = r.prof_z->GetBinContent(b);
                if (!std::isfinite(v)) continue;
                ymin = std::min(ymin, v);
                ymax = std::max(ymax, v);
            }
        }
        double pad = (ymax - ymin) * 0.1;
        ymin = std::min(ymin - pad, -10.0);
        ymax = std::max(ymax + pad, 260.0);

        // Ideal diagonal
        TLine *diag = new TLine(0, 0, 250, 250);
        diag->SetLineColor(kBlack); diag->SetLineStyle(2); diag->SetLineWidth(2);

        bool first = true;
        for (const auto &r : results) {
            // Mask bins with < 5 entries to avoid noisy points
            for (int b = 1; b <= r.prof_z->GetNbinsX(); ++b)
                if (r.prof_z->GetBinEntries(b) < 5)
                    r.prof_z->SetBinContent(b, 0), r.prof_z->SetBinEntries(b, 1e-9);

            if (first) {
                r.prof_z->SetTitle(
                    "Reconstructed vs True decay vertex z;"
                    "True vertex z [cm];"
                    "Mean reconstructed vertex z [cm]");
                r.prof_z->GetXaxis()->SetRangeUser(0, 250);
                r.prof_z->GetYaxis()->SetRangeUser(ymin, ymax);
                r.prof_z->Draw("EP");
                diag->Draw("SAME");
                first = false;
            } else {
                r.prof_z->Draw("EP SAME");
            }
        }

        TLatex *ann = new TLatex(10, 235, "Dashed: ideal (z_{reco} = z_{true})");
        ann->SetTextSize(0.030); ann->SetTextColor(kBlack); ann->Draw();

        c2->cd(2);
        gPad->SetPad(0.75, 0.00, 1.00, 1.00);
        drawDetectorLayout(config_labels, line_colors);

        c2->Update();
        c2->SaveAs("KLong_validate_vertex_z_profile.png");
        std::cout << "\n=== Saved KLong_validate_vertex_z_profile.png ===" << std::endl;
        delete c2;
    }

    // ── Summary printout ──────────────────────────────────────────────────
    std::cout << "\n=== VERTEX RESIDUAL SUMMARY ===" << std::endl;
    std::cout << Form("%-60s  %8s %8s  %8s %8s  %8s %8s  %8s %8s\n",
        "Config", "mu_dx", "sig_dx", "mu_dy", "sig_dy", "mu_dz", "sig_dz", "mean_dr", "rms_dr");
    for (const auto &r : results) {
        std::cout << Form("%-60s  %+7.2f  %7.2f  %+7.2f  %7.2f  %+7.2f  %7.2f  %7.2f  %7.2f\n",
            r.label.c_str(),
            r.mu_dx, r.sig_dx, r.mu_dy, r.sig_dy, r.mu_dz, r.sig_dz, r.mean_dr, r.sig_dr);
    }
    std::cout << "(all values in cm)" << std::endl;
}

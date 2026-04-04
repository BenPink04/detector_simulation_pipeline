// KLong_validate_reconstruction.C
//
// Validation tests for the momentum reconstruction fixes documented in
// MOMENTUM_RECONSTRUCTION_ISSUES.md.
//
// Test 1: p_reco / p_true vs p_true (per momentum bin, per configuration)
//         Expected after fix: flat at ~1.0 with no slope.
//         Before fix: ratio falls below 1.0, worsening at high p.
//
// Test 2: Integrated Dp = p_reco - p_true distribution (all momenta, all configs)
//         Expected after fix: Gaussian centred at 0 for all configurations.
//         Before fix: mean is significantly negative (systematic underestimate).
//
// Usage:
//   root -l KLong_validate_reconstruction.C
//
// Output files (written to the current directory):
//   KLong_validate_ratio.png      -- p_reco/p_true vs p_true
//   KLong_validate_delta_p.png    -- integrated Dp distribution

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
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

// ── Detector layout helpers (shared with other macros) ────────────────────────

struct DetectorLayout {
    std::vector<double> trackers;
    std::vector<double> pizzas;
    std::vector<double> fris;
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
    for (const auto& config : config_labels) {
        DetectorLayout l = parseConfigString(config);
        for (double p : l.trackers) { min_pos = std::min(min_pos, p); max_pos = std::max(max_pos, p); }
        for (double p : l.pizzas)   { min_pos = std::min(min_pos, p); max_pos = std::max(max_pos, p); }
        for (double p : l.fris)     { min_pos = std::min(min_pos, p); max_pos = std::max(max_pos, p); }
        min_pos = std::min(min_pos, l.end_detector);
        max_pos = std::max(max_pos, l.end_detector);
    }
    double range = max_pos - min_pos;
    min_pos -= range * 0.1;
    max_pos += range * 0.1;

    TText *title = new TText(0.5, 0.98, "Detector Layouts");
    title->SetTextAlign(23); title->SetTextSize(0.04); title->Draw();

    int n = (int)config_labels.size();
    double y_spacing = 0.85 / (n + 1);

    for (int i = 0; i < n; ++i) {
        DetectorLayout l = parseConfigString(config_labels[i]);
        double yc = 0.90 - (i + 1) * y_spacing;

        TLine *ind = new TLine(0.02, yc, 0.08, yc);
        ind->SetLineColor(line_colors[i]); ind->SetLineWidth(4); ind->Draw();

        const double xs = 0.80, xo = 0.15;
        auto xpos = [&](double pos) { return xo + (pos - min_pos) / (max_pos - min_pos) * xs; };

        for (double p : l.trackers) {
            double x = xpos(p);
            TBox *b = new TBox(x-0.008, yc-0.015, x+0.008, yc+0.015);
            b->SetFillColor(kBlue); b->SetLineColor(kBlue); b->Draw();
        }
        for (double p : l.pizzas) {
            double x = xpos(p);
            TBox *b = new TBox(x-0.008, yc-0.015, x+0.008, yc+0.015);
            b->SetFillColor(kRed); b->SetLineColor(kRed); b->Draw();
        }
        for (double p : l.fris) {
            double x = xpos(p);
            TBox *b = new TBox(x-0.008, yc-0.015, x+0.008, yc+0.015);
            b->SetFillColor(kMagenta); b->SetLineColor(kMagenta); b->Draw();
        }
        double x = xpos(l.end_detector);
        TBox *b = new TBox(x-0.008, yc-0.015, x+0.008, yc+0.015);
        b->SetFillColor(kGreen); b->SetLineColor(kGreen); b->Draw();
    }

    // Detector-type legend
    double ly = 0.08, lx = 0.15, lsp = 0.18;
    auto addLegItem = [&](double x, int col, const char *label) {
        TBox *b = new TBox(x, ly, x+0.02, ly+0.03);
        b->SetFillColor(col); b->Draw();
        TText *t = new TText(x+0.03, ly+0.015, label);
        t->SetTextSize(0.03); t->Draw();
    };
    addLegItem(lx,          kBlue,    "Tracker");
    addLegItem(lx+lsp,      kRed,     "Pizza");
    addLegItem(lx+2*lsp,    kMagenta, "FRI");
    addLegItem(lx+3*lsp,    kGreen,   "End");
}

// ── Main macro ────────────────────────────────────────────────────────────────

void KLong_validate_reconstruction() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // ── Discover combined_vectors.root files ───────────────────────────────
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

    // ── Momentum binning for ratio plot ───────────────────────────────────
    //    18 bins, 0.5–10 GeV/c (matching gaussian_spreads_truncated)
    const int nBins = 18;
    double pBins[nBins+1] = {0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
                              2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 10.0};

    int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kCyan+1, kOrange+1, kViolet, kTeal-5};

    // Per-config storage
    struct ConfigResult {
        std::string label;
        int color;
        TGraph  *g_ratio   = nullptr;   // Test 1: mean(p_reco/p_true) per bin
        TGraph  *g_ratio_lo = nullptr;  // ±1σ band lower edge
        TGraph  *g_ratio_hi = nullptr;  // ±1σ band upper edge
        TH1D    *h_delta_p  = nullptr;  // Test 2: integrated Dp
        double   fit_mean   = 0;
        double   fit_sigma  = 0;
    };

    std::vector<ConfigResult> results;
    std::vector<int> line_colors;

    // ── Loop over configurations ───────────────────────────────────────────
    for (size_t idx = 0; idx < root_files.size(); ++idx) {
        ConfigResult res;
        res.label = config_labels[idx];
        res.color = colors[idx % 8];
        line_colors.push_back(res.color);

        TFile *file = TFile::Open(root_files[idx].c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening " << root_files[idx] << "\n";
            continue;
        }

        TTree *tree = (TTree*)file->Get("kaonVectors");
        if (!tree) {
            std::cerr << "Error: kaonVectors tree not found in " << root_files[idx] << "\n";
            file->Close(); continue;
        }

        std::vector<double> *reco_p = nullptr, *true_p = nullptr;
        tree->SetBranchAddress("reco_p", &reco_p);
        tree->SetBranchAddress("true_p", &true_p);

        // Per-bin accumulators for ratio
        std::vector<double> sum_ratio(nBins, 0), sum_ratio2(nBins, 0);
        std::vector<int>    count(nBins, 0);

        // Integrated Dp histogram: range −3 to +3 GeV/c, 120 bins
        std::string hname = Form("h_dp_%zu", idx);
        res.h_delta_p = new TH1D(hname.c_str(), "", 120, -3.0, 3.0);
        res.h_delta_p->SetDirectory(nullptr);

        const double anomaly_threshold = 1.0; // Filter |delta_p/true_p| > 100% as unphysical outliers
        Long64_t nEnt = tree->GetEntries();
        for (Long64_t e = 0; e < nEnt; ++e) {
            tree->GetEntry(e);
            for (size_t i = 0; i < reco_p->size(); ++i) {
                double pt  = (*true_p)[i];
                double pr  = (*reco_p)[i];

                // Skip unphysical / non-finite entries
                if (!std::isfinite(pt) || !std::isfinite(pr)) continue;
                if (pt <= 0 || pr <= 0) continue;
                if (std::abs(pr - pt) / pt > anomaly_threshold) continue;

                double dp  = pr - pt;
                double rat = pr / pt;

                // Fill ratio bins
                for (int b = 0; b < nBins; ++b) {
                    if (pt >= pBins[b] && pt < pBins[b+1]) {
                        sum_ratio[b]  += rat;
                        sum_ratio2[b] += rat * rat;
                        count[b]++;
                        break;
                    }
                }

                // Fill integrated Dp
                res.h_delta_p->Fill(dp);
            }
        }
        file->Close();

        // Build ratio TGraph (mean ± RMS per bin)
        std::vector<double> px, py, py_lo, py_hi;
        for (int b = 0; b < nBins; ++b) {
            if (count[b] < 20) continue;
            double mean  = sum_ratio[b] / count[b];
            double mean2 = sum_ratio2[b] / count[b];
            double rms   = std::sqrt(std::max(0.0, mean2 - mean*mean));
            double pcen  = 0.5 * (pBins[b] + pBins[b+1]);
            px.push_back(pcen);
            py.push_back(mean);
            py_lo.push_back(mean - rms);
            py_hi.push_back(mean + rms);
        }

        if (!px.empty()) {
            res.g_ratio    = new TGraph((int)px.size(), &px[0], &py[0]);
            res.g_ratio_lo = new TGraph((int)px.size(), &px[0], &py_lo[0]);
            res.g_ratio_hi = new TGraph((int)px.size(), &px[0], &py_hi[0]);

            for (auto *g : {res.g_ratio, res.g_ratio_lo, res.g_ratio_hi}) {
                g->SetLineColor(res.color);
                g->SetLineWidth(2);
                g->SetMarkerColor(res.color);
            }
            res.g_ratio->SetMarkerStyle(20 + (int)idx);
            res.g_ratio->SetMarkerSize(1.1);
            res.g_ratio_lo->SetLineStyle(2);
            res.g_ratio_hi->SetLineStyle(2);
        }

        // Fit Dp histogram
        if (res.h_delta_p->GetEntries() > 50) {
            double hmean = res.h_delta_p->GetMean();
            double hrms  = res.h_delta_p->GetRMS();
            TF1 *gfit = new TF1(Form("gfit_%zu", idx), "gaus",
                                hmean - 2*hrms, hmean + 2*hrms);
            gfit->SetParameters(res.h_delta_p->GetMaximum(), hmean, hrms);
            res.h_delta_p->Fit(gfit, "RQN");
            res.fit_mean  = gfit->GetParameter(1);
            res.fit_sigma = gfit->GetParameter(2);
            delete gfit;

            res.h_delta_p->SetLineColor(res.color);
            res.h_delta_p->SetLineWidth(2);
        }

        results.push_back(res);
    }

    if (results.empty()) { std::cerr << "No valid results\n"; return; }

    // ─────────────────────────────────────────────────────────────────────────
    // CANVAS 1: p_reco / p_true  vs  p_true
    // ─────────────────────────────────────────────────────────────────────────
    {
        TCanvas *c1 = new TCanvas("c_ratio", "p_reco / p_true vs p_true", 1400, 700);
        c1->Divide(2, 1);

        // Left: ratio plot
        c1->cd(1);
        gPad->SetPad(0.0, 0.0, 0.75, 1.0);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.12);
        gPad->SetGrid(1, 1);

        // Reference line at 1.0
        TLine *ref = new TLine(1.0, 1.0, 4.0, 1.0);
        ref->SetLineColor(kBlack); ref->SetLineStyle(2); ref->SetLineWidth(2);

        // Find y range (restrict x to 1–4 GeV/c)
        double ymin = 1e9, ymax = -1e9;
        for (const auto &r : results) {
            if (!r.g_ratio_lo || !r.g_ratio_hi) continue;
            for (int i = 0; i < r.g_ratio_lo->GetN(); ++i) {
                double x = r.g_ratio_lo->GetX()[i];
                if (x < 1.0 || x > 4.0) continue;
                ymin = std::min(ymin, r.g_ratio_lo->GetY()[i]);
                ymax = std::max(ymax, r.g_ratio_hi->GetY()[i]);
            }
        }
        double ypad = (ymax - ymin) * 0.15;
        ymin = std::min(ymin - ypad, 0.85);
        ymax = std::max(ymax + ypad, 1.15);

        bool first = true;
        for (const auto &r : results) {
            if (!r.g_ratio) continue;
            if (first) {
                r.g_ratio->SetTitle(
                    "Momentum Ratio;True momentum p_{true} [GeV/c];#LTp_{reco}/p_{true}#GT");
                r.g_ratio->GetXaxis()->SetLimits(1.0, 4.0);
                r.g_ratio->GetYaxis()->SetRangeUser(ymin, ymax);
                r.g_ratio->Draw("ALP");
                ref->Draw("SAME");
                first = false;
            } else {
                r.g_ratio->Draw("LP SAME");
            }
            r.g_ratio_lo->Draw("L SAME");
            r.g_ratio_hi->Draw("L SAME");
        }

        // // Legend
        // TLegend *leg = new TLegend(0.14, 0.14, 0.55, 0.14 + 0.05 * results.size());
        // leg->SetBorderSize(0);
        // leg->SetFillStyle(0);
        // leg->SetTextSize(0.028);
        // for (const auto &r : results) {
        //     if (!r.g_ratio) continue;
        //     // Shorten config label: just show T1-XXX
        //     std::string shortlabel = r.label.substr(0, r.label.find('_'));
        //     leg->AddEntry(r.g_ratio, shortlabel.c_str(), "lp");
        // }
        // leg->Draw();

        // Annotation: ideal line
        TLatex *ann = new TLatex(1.05, 1.01, "Ideal (ratio = 1.0)");
        ann->SetTextSize(0.027); ann->SetTextColor(kBlack); ann->Draw();

        // Right: detector layout schematic
        c1->cd(2);
        gPad->SetPad(0.75, 0.0, 1.0, 1.0);
        drawDetectorLayout(config_labels, line_colors);

        c1->Update();
        c1->SaveAs("KLong_validate_ratio.png");
        std::cout << "\n=== Saved KLong_validate_ratio.png ===" << std::endl;

        // Print per-bin summary
        std::cout << "\n=== RATIO SUMMARY (mean p_reco/p_true per bin) ===" << std::endl;
        for (const auto &r : results) {
            if (!r.g_ratio) continue;
            std::cout << r.label << "\n";
            for (int i = 0; i < r.g_ratio->GetN(); ++i) {
                std::cout << Form("  p_true = %4.2f GeV/c  ->  ratio = %.4f\n",
                                  r.g_ratio->GetX()[i], r.g_ratio->GetY()[i]);
            }
        }
        delete c1;
    }

    // ─────────────────────────────────────────────────────────────────────────
    // CANVAS 2: Integrated Dp = p_reco - p_true  distribution
    // ─────────────────────────────────────────────────────────────────────────
    {
        TCanvas *c2 = new TCanvas("c_dp", "#Deltap Distribution (all momenta)", 1400, 700);
        c2->Divide(2, 1);

        // Left: Dp histogram
        c2->cd(1);
        gPad->SetPad(0.0, 0.0, 0.75, 1.0);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.12);
        gPad->SetGrid(1, 0);

        // Find display x range: ±max(|fit_mean| + 3*fit_sigma) across configs
        double dp_max = 0;
        for (const auto &r : results)
            dp_max = std::max(dp_max, std::abs(r.fit_mean) + 3.5 * r.fit_sigma);
        dp_max = std::min(dp_max + 0.2, 3.0);

        double hmax = 0;
        for (const auto &r : results) {
            if (!r.h_delta_p) continue;
            r.h_delta_p->GetXaxis()->SetRangeUser(-dp_max, dp_max);
            hmax = std::max(hmax, r.h_delta_p->GetMaximum());
        }

        // Reference line at Dp = 0
        TLine *zero = new TLine(0.0, 0.0, 0.0, hmax * 1.25);
        zero->SetLineColor(kBlack); zero->SetLineStyle(2); zero->SetLineWidth(2);

        bool first = true;
        for (const auto &r : results) {
            if (!r.h_delta_p) continue;
            if (first) {
                r.h_delta_p->SetTitle(
                    "#Deltap = p_{reco} #minus p_{true} (all p_{true});#Deltap [GeV/c];Entries");
                r.h_delta_p->GetXaxis()->SetRangeUser(-dp_max, dp_max);
                r.h_delta_p->GetYaxis()->SetRangeUser(0, hmax * 1.25);
                r.h_delta_p->Draw("HIST");
                zero->Draw("SAME");
                first = false;
            } else {
                r.h_delta_p->Draw("HIST SAME");
            }
        }

        // // Legend with fit results
        // TLegend *leg = new TLegend(0.14, 0.60, 0.62, 0.60 + 0.065 * results.size());
        // leg->SetBorderSize(0);
        // leg->SetFillStyle(0);
        // leg->SetTextSize(0.027);
        // for (const auto &r : results) {
        //     if (!r.h_delta_p) continue;
        //     std::string shortlabel = r.label.substr(0, r.label.find('_'));
        //     std::string entry = Form("%s   #mu=%.3f  #sigma=%.3f GeV/c",
        //                               shortlabel.c_str(), r.fit_mean, r.fit_sigma);
        //     leg->AddEntry(r.h_delta_p, entry.c_str(), "l");
        // }
        // leg->Draw();

        // Right: detector layout schematic
        c2->cd(2);
        gPad->SetPad(0.75, 0.0, 1.0, 1.0);
        drawDetectorLayout(config_labels, line_colors);

        c2->Update();
        c2->SaveAs("KLong_validate_delta_p.png");
        std::cout << "\n=== Saved KLong_validate_delta_p.png ===" << std::endl;

        // Print fit summary
        std::cout << "\n=== DELTA_P FIT SUMMARY ===" << std::endl;
        for (const auto &r : results) {
            if (!r.h_delta_p) continue;
            std::cout << Form("%-70s  mu = %+.4f GeV/c   sigma = %.4f GeV/c\n",
                              r.label.c_str(), r.fit_mean, r.fit_sigma);
        }
        delete c2;
    }
}

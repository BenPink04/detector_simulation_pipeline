// ============================================================================
// KLong_interpolate_pizza.C
//
// PURPOSE:
//   Builds a continuous approximation of the KLong detector's momentum
//   resolution and acceptance as a function of pizza detector position.
//
//   For each simulated pizza configuration found in ARCHIVED_RESULTS the
//   macro:
//     1. Reads the combined_vectors.root file (tree: kaonVectors,
//        branches: reco_p, true_p) and builds a resolution-vs-momentum
//        profile.  The profile mean per bin is fitted with the standard
//        two-term tracking resolution formula:
//
//            sigma_p/p(p) = sqrt( a^2 + (b*p)^2 )
//
//        where 'a' captures the multiple-scattering dominated low-p regime
//        and 'b' the hit-resolution dominated high-p regime.
//
//     2. Reads the combined_acceptance.root file (tree: kaonEventInfo,
//        branches: true_mom_vec, reco_flag_vec) and builds an efficiency
//        histogram (accepted / total per momentum bin).  This is fitted
//        with a sigmoid:
//
//            acc(p) = norm / (1 + exp(-k*(p - p0)))
//
//     3. Writes all extracted fit parameters to:
//          - pizza_fit_params.csv   (human-readable, for Python / inspection)
//          - pizza_fit_params.root  (TTree + TSpline3 objects for ROOT queries)
//
//   After building the parameter table a demo "query" is run that, given a
//   requested pizza position P1, interpolates a_res, b_res and the acceptance
//   sigmoid parameters and draws the resulting resolution and acceptance
//   curves.
//
// USAGE (from ROOT prompt or .x):
//   root -l -q 'KLong_interpolate_pizza.C'
//
//   To customise which archive directories to include, edit the vector
//   'ARCHIVE_DIRS' near the top of the main function.
//
// OUTPUT:
//   INTERPOLATION/pizza_fit_params.csv
//   INTERPOLATION/pizza_fit_params.root   (TTree "fitParams", TSpline3 objects)
//   INTERPOLATION/query_demo_P1_<X>.png   (resolution + acceptance for queried P1)
//
// DEPENDENCIES:
//   Standard ROOT (TFile, TTree, TH1, TH2, TF1, TSpline, TCanvas, TGraphErrors)
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

// ============================================================================
// CONFIGURATION
// ============================================================================

// Momentum binning for profiles
static const int    N_P_BINS   = 30;
static const double P_MIN      = 0.5;   // GeV/c  — avoid very low stats near 0
static const double P_MAX      = 9.0;   // GeV/c

// Resolution quality cuts
static const double RES_MAX    = 0.8;   // reject |reco-true|/true > this
static const int    MIN_COUNTS = 10;    // minimum entries per bin for fit

// ============================================================================
// HELPERS
// ============================================================================

// Extract a named detector position (e.g. "P1", "P2") from a config string
// like  T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700
double extractPosition(const std::string& label, const std::string& detName) {
    // look for  _<detName>-<value>_  or  <detName>-<value> at start/end
    std::string srch = detName + "-";
    size_t pos = label.find(srch);
    if (pos == std::string::npos) return -1.0;
    pos += srch.size();
    // read digits (may be 3 digit, could also be 0 for disabled)
    size_t end = label.find_first_not_of("0123456789", pos);
    return std::stod(label.substr(pos, end - pos));
}

// Strip path to get the bare config label (without _combined_vectors.root etc.)
std::string baseLabel(const std::string& filepath) {
    size_t slash = filepath.find_last_of('/');
    std::string base = (slash != std::string::npos) ? filepath.substr(slash + 1) : filepath;
    for (const std::string& suf : {"_combined_vectors.root", "_combined_acceptance.root", ".root"}) {
        size_t p = base.find(suf);
        if (p != std::string::npos) { base = base.substr(0, p); break; }
    }
    return base;
}

// ============================================================================
// PER-CONFIG DATA
// ============================================================================

struct ConfigParams {
    std::string label;
    std::string archive;
    double P1, P2;

    // Resolution fit:  sigma_p/p = sqrt(a^2 + (b*p)^2)
    double a_res, a_res_err;
    double b_res, b_res_err;
    double res_fit_chi2ndf;
    bool   res_fit_ok;

    // Acceptance fit:  norm / (1 + exp(-k*(p-p0)))
    double acc_norm,  acc_norm_err;
    double acc_k,     acc_k_err;
    double acc_p0,    acc_p0_err;
    double acc_fit_chi2ndf;
    bool   acc_fit_ok;
};

// ============================================================================
// FIT RESOLUTION PROFILE
// ============================================================================

ConfigParams fitResolution(const std::string& vecFile, const std::string& accFile) {
    ConfigParams cp;
    cp.label      = baseLabel(vecFile);
    // extract archive name from path
    {
        const std::string marker = "ARCHIVED_RESULTS/";
        size_t pos = vecFile.find(marker);
        if (pos != std::string::npos) {
            pos += marker.size();
            size_t end = vecFile.find('/', pos);
            cp.archive = vecFile.substr(pos, end == std::string::npos ? std::string::npos : end - pos);
        } else {
            cp.archive = "UNKNOWN";
        }
    }
    cp.P1 = extractPosition(cp.label, "P1");
    cp.P2 = extractPosition(cp.label, "P2");
    cp.res_fit_ok = false;
    cp.acc_fit_ok = false;
    cp.a_res = cp.b_res = cp.a_res_err = cp.b_res_err = cp.res_fit_chi2ndf = 0;
    cp.acc_norm = cp.acc_k = cp.acc_p0 = cp.acc_norm_err = cp.acc_k_err = cp.acc_p0_err = cp.acc_fit_chi2ndf = 0;

    // ---- Resolution from vectors file ----
    TFile* fv = TFile::Open(vecFile.c_str(), "READ");
    if (!fv || fv->IsZombie()) {
        std::cerr << "  [WARN] Cannot open vectors file: " << vecFile << std::endl;
        if (fv) delete fv;
        goto DO_ACCEPTANCE;  // still try acceptance
    }
    {
        TTree* tv = (TTree*)fv->Get("kaonVectors");
        if (!tv) {
            std::cerr << "  [WARN] No kaonVectors tree in " << vecFile << std::endl;
            fv->Close(); delete fv;
            goto DO_ACCEPTANCE;
        }

        std::vector<double>* reco_p = nullptr;
        std::vector<double>* true_p = nullptr;
        tv->SetBranchAddress("reco_p", &reco_p);
        tv->SetBranchAddress("true_p", &true_p);

        // 2D histogram: true_p vs relative resolution
        TH2D h2("h2res", "", N_P_BINS, P_MIN, P_MAX, 80, 0, RES_MAX);

        Long64_t nEnt = tv->GetEntries();
        for (Long64_t ev = 0; ev < nEnt; ++ev) {
            tv->GetEntry(ev);
            if (!reco_p || !true_p) continue;
            for (size_t i = 0; i < reco_p->size(); ++i) {
                double tp = (*true_p)[i];
                if (tp <= 0) continue;
                double res = std::abs((*reco_p)[i] - tp) / tp;
                if (res < RES_MAX) h2.Fill(tp, res);
            }
        }
        fv->Close(); delete fv;

        // Build profile: mean resolution per bin
        // Use TGraphErrors for fitting
        std::vector<double> gx, gy, gex, gey;
        for (int ib = 1; ib <= N_P_BINS; ++ib) {
            TH1D* proj = h2.ProjectionY("_py_tmp", ib, ib);
            if (proj->GetEntries() >= MIN_COUNTS) {
                double bc  = h2.GetXaxis()->GetBinCenter(ib);
                double bw  = h2.GetXaxis()->GetBinWidth(ib) * 0.5;
                gx.push_back(bc);
                gy.push_back(proj->GetMean());
                gex.push_back(bw);
                gey.push_back(proj->GetMeanError());
            }
            delete proj;
        }

        if ((int)gx.size() < 4) {
            std::cerr << "  [WARN] Too few bins for resolution fit in " << cp.label << std::endl;
            goto DO_ACCEPTANCE;
        }

        TGraphErrors gr((int)gx.size(), gx.data(), gy.data(), gex.data(), gey.data());
        gr.SetName("gr_res_tmp");

        // Fit:  f(p) = sqrt(a^2 + (b*p)^2)
        // Written as sqrt([0]^2 + ([1]*x)^2)
        TF1 fRes("fRes", "sqrt([0]*[0] + ([1]*x)*([1]*x))", P_MIN, P_MAX);
        fRes.SetParNames("a", "b");
        fRes.SetParameters(0.05, 0.02);
        fRes.SetParLimits(0, 1e-4, 1.0);
        fRes.SetParLimits(1, 1e-5, 0.5);

        TFitResultPtr rfit = gr.Fit(&fRes, "QSRN");

        if ((int)rfit == 0 || rfit->IsValid()) {
            cp.a_res      = fRes.GetParameter(0);
            cp.a_res_err  = fRes.GetParError(0);
            cp.b_res      = fRes.GetParameter(1);
            cp.b_res_err  = fRes.GetParError(1);
            double ndf    = fRes.GetNDF();
            cp.res_fit_chi2ndf = (ndf > 0) ? fRes.GetChisquare() / ndf : -1;
            cp.res_fit_ok = (cp.a_res > 0 && cp.b_res > 0);
        } else {
            std::cerr << "  [WARN] Resolution fit failed for " << cp.label << std::endl;
        }
    }

    DO_ACCEPTANCE:
    // ---- Acceptance from acceptance file ----
    if (accFile.empty()) return cp;

    TFile* fa = TFile::Open(accFile.c_str(), "READ");
    if (!fa || fa->IsZombie()) {
        std::cerr << "  [WARN] Cannot open acceptance file: " << accFile << std::endl;
        if (fa) delete fa;
        return cp;
    }
    {
        TTree* ta = (TTree*)fa->Get("kaonEventInfo");
        if (!ta) {
            std::cerr << "  [WARN] No kaonEventInfo tree in " << accFile << std::endl;
            fa->Close(); delete fa;
            return cp;
        }

        std::vector<double>* true_mom = nullptr;
        std::vector<int>*    reco_flag = nullptr;
        ta->SetBranchAddress("true_mom_vec",   &true_mom);
        ta->SetBranchAddress("reco_flag_vec",  &reco_flag);

        TH1D hAll ("hAll",  "", N_P_BINS, P_MIN, P_MAX);
        TH1D hReco("hReco", "", N_P_BINS, P_MIN, P_MAX);

        Long64_t nEnt = ta->GetEntries();
        for (Long64_t ev = 0; ev < nEnt; ++ev) {
            ta->GetEntry(ev);
            if (!true_mom || !reco_flag) continue;
            for (size_t i = 0; i < true_mom->size(); ++i) {
                double tp = (*true_mom)[i];
                hAll.Fill(tp);
                if ((*reco_flag)[i] == 1) hReco.Fill(tp);
            }
        }
        fa->Close(); delete fa;

        // Build efficiency graph
        std::vector<double> gx, gy, gex, gey;
        for (int ib = 1; ib <= N_P_BINS; ++ib) {
            double nAll  = hAll.GetBinContent(ib);
            double nReco = hReco.GetBinContent(ib);
            if (nAll < MIN_COUNTS) continue;
            double eff = nReco / nAll;
            double err = std::sqrt(eff * (1 - eff) / nAll);
            gx.push_back(hAll.GetXaxis()->GetBinCenter(ib));
            gy.push_back(eff);
            gex.push_back(hAll.GetXaxis()->GetBinWidth(ib) * 0.5);
            gey.push_back(err);
        }

        if ((int)gx.size() < 4) {
            std::cerr << "  [WARN] Too few bins for acceptance fit in " << cp.label << std::endl;
            return cp;
        }

        TGraphErrors gacc((int)gx.size(), gx.data(), gy.data(), gex.data(), gey.data());
        gacc.SetName("gacc_tmp");

        // Sigmoid fit:  norm / (1 + exp(-k*(x-p0)))
        TF1 fAcc("fAcc", "[0] / (1 + exp(-[1]*(x-[2])))", P_MIN, P_MAX);
        fAcc.SetParNames("norm", "k", "p0");
        // Initial guesses: plateau near mean gy, k~1, midpoint near middle of range
        double gy_max = *std::max_element(gy.begin(), gy.end());
        fAcc.SetParameters(gy_max > 0 ? gy_max : 0.5, 2.0, 2.0);
        fAcc.SetParLimits(0, 0.0, 1.0);
        fAcc.SetParLimits(1, 0.01, 50.0);
        fAcc.SetParLimits(2, P_MIN, P_MAX);

        TFitResultPtr racc = gacc.Fit(&fAcc, "QSRN");

        if ((int)racc == 0 || racc->IsValid()) {
            cp.acc_norm      = fAcc.GetParameter(0);
            cp.acc_norm_err  = fAcc.GetParError(0);
            cp.acc_k         = fAcc.GetParameter(1);
            cp.acc_k_err     = fAcc.GetParError(1);
            cp.acc_p0        = fAcc.GetParameter(2);
            cp.acc_p0_err    = fAcc.GetParError(2);
            double ndf       = fAcc.GetNDF();
            cp.acc_fit_chi2ndf = (ndf > 0) ? fAcc.GetChisquare() / ndf : -1;
            cp.acc_fit_ok    = true;
        } else {
            std::cerr << "  [WARN] Acceptance fit failed for " << cp.label << std::endl;
        }
    }

    return cp;
}

// ============================================================================
// QUERY FUNCTION
// Interpolates fitted parameters at an arbitrary pizza position P1_query
// and draws the resulting resolution and acceptance curves.
// ============================================================================

void queryPizzaPosition(double P1_query,
                        TSpline3* sp_a,    TSpline3* sp_b,
                        TSpline3* sp_norm, TSpline3* sp_k, TSpline3* sp_p0,
                        const std::string& outDir) {
    double a    = sp_a->Eval(P1_query);
    double b    = sp_b->Eval(P1_query);
    double norm = sp_norm->Eval(P1_query);
    double k    = sp_k->Eval(P1_query);
    double p0   = sp_p0->Eval(P1_query);

    // Clamp to physically meaningful values
    if (a    < 0)    a    = 0;
    if (b    < 0)    b    = 0;
    if (norm < 0)    norm = 0;
    if (norm > 1.0)  norm = 1.0;
    if (k    < 0.01) k    = 0.01;

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Query: P1 = " << P1_query << " cm" << std::endl;
    std::cout << "  Resolution params:  a = " << a << ",  b = " << b << std::endl;
    std::cout << "  Acceptance params:  norm = " << norm
              << ",  k = " << k << ",  p0 = " << p0 << std::endl;
    std::cout << "  Resolution at 1 GeV/c : " << std::sqrt(a*a + (b*1.0)*(b*1.0)) << std::endl;
    std::cout << "  Resolution at 3 GeV/c : " << std::sqrt(a*a + (b*3.0)*(b*3.0)) << std::endl;
    std::cout << "  Resolution at 5 GeV/c : " << std::sqrt(a*a + (b*5.0)*(b*5.0)) << std::endl;
    std::cout << "  Acceptance plateau    : " << norm << std::endl;
    std::cout << "========================================\n" << std::endl;

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cQuery", "Pizza Interpolation Query", 1200, 500);
    c->Divide(2, 1);

    // -- Resolution plot --
    c->cd(1);
    gPad->SetGrid(1, 1);
    TF1* fRes = new TF1("fResQuery", "sqrt([0]*[0]+([1]*x)*([1]*x))", P_MIN, P_MAX);
    fRes->SetParameters(a, b);
    fRes->SetLineColor(kBlue+1);
    fRes->SetLineWidth(3);
    fRes->GetXaxis()->SetTitle("True Kaon Momentum [GeV/c]");
    fRes->GetYaxis()->SetTitle("#sigma_{p}/p");
    fRes->SetTitle(Form("Momentum Resolution  (P1 = %.0f cm)", P1_query));
    fRes->GetYaxis()->SetRangeUser(0, 0.3);
    fRes->Draw();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.DrawLatex(0.15, 0.82, Form("a = %.4f,  b = %.4f", a, b));

    // -- Acceptance plot --
    c->cd(2);
    gPad->SetGrid(1, 1);
    TF1* fAcc = new TF1("fAccQuery", "[0]/(1+exp(-[1]*(x-[2])))", P_MIN, P_MAX);
    fAcc->SetParameters(norm, k, p0);
    fAcc->SetLineColor(kRed+1);
    fAcc->SetLineWidth(3);
    fAcc->GetXaxis()->SetTitle("True Kaon Momentum [GeV/c]");
    fAcc->GetYaxis()->SetTitle("Acceptance");
    fAcc->SetTitle(Form("Acceptance  (P1 = %.0f cm)", P1_query));
    fAcc->GetYaxis()->SetRangeUser(0, 1.05);
    fAcc->Draw();
    latex.DrawLatex(0.15, 0.82, Form("norm=%.3f  k=%.3f  p_{0}=%.2f", norm, k, p0));

    std::string outFile = outDir + Form("/query_demo_P1_%.0f.png", P1_query);
    c->SaveAs(outFile.c_str());
    std::cout << "  Plot saved: " << outFile << std::endl;

    delete c;
}

// ============================================================================
// MAIN
// ============================================================================

void KLong_interpolate_pizza() {

    gStyle->SetOptStat(0);

    const std::string BASE_DIR   = "/users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS";
    const std::string OUTPUT_DIR = "/users/bp969/scratch/VIKING_FOLDER/INTERPOLATION";
    gSystem->mkdir(OUTPUT_DIR.c_str(), kTRUE);

    // ---- Select which archive directories to include ----
    // Comment out any you don't want.
    std::vector<std::string> ARCHIVE_DIRS = {
        "Diss_Draft_2_20260415",
        "Upstream_pizzas_20260416",
        "EARLY_TRACKERS_20260217",
        "KLong_Diss_Draft_Figs_20260413",
        // "BETA_VALIDATION_20260404",
        // "FINAL_TEST_20260413",
    };

    // -------------------------------------------------------------------------
    // Step 1: collect all combined_vectors.root files and match with
    //         combined_acceptance.root where available
    // -------------------------------------------------------------------------
    // Map from config label -> {vectors file, acceptance file}
    std::map<std::string, std::pair<std::string, std::string>> configFiles;

    for (const auto& archDir : ARCHIVE_DIRS) {
        std::string fullDir = BASE_DIR + "/" + archDir;
        std::string findCmd = "find " + fullDir + " -name '*combined_vectors.root' 2>/dev/null | sort";
        FILE* pipe = popen(findCmd.c_str(), "r");
        if (!pipe) continue;
        char buf[1024];
        while (fgets(buf, sizeof(buf), pipe)) {
            std::string vf(buf);
            vf.erase(std::remove(vf.begin(), vf.end(), '\n'), vf.end());
            if (vf.empty()) continue;
            std::string lbl = baseLabel(vf);
            // prefer most recent archive (map insertion only if not already present)
            if (configFiles.find(lbl) == configFiles.end()) {
                configFiles[lbl] = {vf, ""};
            }
        }
        pclose(pipe);

        // Now find acceptance files in same directory
        std::string findAcc = "find " + fullDir + " -name '*combined_acceptance.root' 2>/dev/null | sort";
        pipe = popen(findAcc.c_str(), "r");
        if (!pipe) continue;
        while (fgets(buf, sizeof(buf), pipe)) {
            std::string af(buf);
            af.erase(std::remove(af.begin(), af.end(), '\n'), af.end());
            if (af.empty()) continue;
            std::string lbl = baseLabel(af);
            if (configFiles.find(lbl) != configFiles.end() &&
                configFiles[lbl].second.empty()) {
                configFiles[lbl].second = af;
            }
        }
        pclose(pipe);
    }

    std::cout << "Found " << configFiles.size() << " unique configurations to process." << std::endl;

    // -------------------------------------------------------------------------
    // Step 2: extract fit parameters for every configuration
    // -------------------------------------------------------------------------
    std::vector<ConfigParams> allParams;
    allParams.reserve(configFiles.size());

    int idx = 0;
    for (const auto& entry : configFiles) {
        ++idx;
        const std::string& lbl = entry.first;
        const std::string& vf  = entry.second.first;
        const std::string& af  = entry.second.second;

        std::cout << "\n[" << idx << "/" << configFiles.size() << "] " << lbl << std::endl;
        if (!af.empty()) std::cout << "  acceptance: " << af << std::endl;

        ConfigParams cp = fitResolution(vf, af);
        allParams.push_back(cp);

        std::cout << "  P1=" << cp.P1 << "  P2=" << cp.P2;
        if (cp.res_fit_ok)
            std::cout << "  a=" << cp.a_res << "  b=" << cp.b_res
                      << "  chi2/ndf=" << cp.res_fit_chi2ndf;
        else
            std::cout << "  [res fit FAILED]";
        if (cp.acc_fit_ok)
            std::cout << "  acc_norm=" << cp.acc_norm
                      << "  k=" << cp.acc_k << "  p0=" << cp.acc_p0;
        else
            std::cout << "  [acc fit FAILED]";
        std::cout << std::endl;
    }

    // -------------------------------------------------------------------------
    // Step 3: write CSV
    // -------------------------------------------------------------------------
    std::string csvPath = OUTPUT_DIR + "/pizza_fit_params.csv";
    std::ofstream csv(csvPath);
    csv << "label,archive,P1,P2,"
        << "a_res,a_res_err,b_res,b_res_err,res_chi2ndf,res_fit_ok,"
        << "acc_norm,acc_norm_err,acc_k,acc_k_err,acc_p0,acc_p0_err,acc_chi2ndf,acc_fit_ok\n";
    for (const auto& cp : allParams) {
        csv << cp.label   << ","
            << cp.archive << ","
            << cp.P1      << ","
            << cp.P2      << ","
            << cp.a_res      << "," << cp.a_res_err << ","
            << cp.b_res      << "," << cp.b_res_err << ","
            << cp.res_fit_chi2ndf << "," << (int)cp.res_fit_ok << ","
            << cp.acc_norm   << "," << cp.acc_norm_err << ","
            << cp.acc_k      << "," << cp.acc_k_err << ","
            << cp.acc_p0     << "," << cp.acc_p0_err << ","
            << cp.acc_fit_chi2ndf << "," << (int)cp.acc_fit_ok << "\n";
    }
    csv.close();
    std::cout << "\nCSV written: " << csvPath << std::endl;

    // -------------------------------------------------------------------------
    // Step 4: build TSpline3 objects over P1 for each parameter
    //         (only use configs where fits succeeded)
    // -------------------------------------------------------------------------
    // Collect good points, sorted by P1
    struct SplinePoint { double P1, a_res, b_res, acc_norm, acc_k, acc_p0; };
    std::vector<SplinePoint> goodRes, goodAcc;

    for (const auto& cp : allParams) {
        if (cp.res_fit_ok) {
            // only one entry per P1 value — keep first encountered
            bool dup = false;
            for (const auto& sp : goodRes) if (std::abs(sp.P1 - cp.P1) < 0.5) { dup = true; break; }
            if (!dup) goodRes.push_back({cp.P1, cp.a_res, cp.b_res, 0, 0, 0});
        }
        if (cp.acc_fit_ok) {
            bool dup = false;
            for (const auto& sp : goodAcc) if (std::abs(sp.P1 - cp.P1) < 0.5) { dup = true; break; }
            if (!dup) goodAcc.push_back({cp.P1, 0, 0, cp.acc_norm, cp.acc_k, cp.acc_p0});
        }
    }

    // Sort by P1
    auto byP1 = [](const SplinePoint& a, const SplinePoint& b) { return a.P1 < b.P1; };
    std::sort(goodRes.begin(), goodRes.end(), byP1);
    std::sort(goodAcc.begin(), goodAcc.end(), byP1);

    std::cout << "\nBuilding splines from " << goodRes.size()
              << " resolution points and " << goodAcc.size() << " acceptance points." << std::endl;

    if (goodRes.size() < 2) {
        std::cerr << "ERROR: Not enough successful resolution fits to build a spline. "
                     "Check that the archive paths are correct and data files exist." << std::endl;
        return;
    }

    // Resolution splines
    std::vector<double> xp, ya, yb;
    for (const auto& sp : goodRes) { xp.push_back(sp.P1); ya.push_back(sp.a_res); yb.push_back(sp.b_res); }

    TSpline3* sp_a = new TSpline3("sp_a_res", xp.data(), ya.data(), (int)xp.size());
    TSpline3* sp_b = new TSpline3("sp_b_res", xp.data(), yb.data(), (int)xp.size());
    sp_a->SetTitle("Spline: a_res vs P1;P1 [cm];a");
    sp_b->SetTitle("Spline: b_res vs P1;P1 [cm];b");

    // Acceptance splines (may be limited)
    TSpline3* sp_norm = nullptr;
    TSpline3* sp_k    = nullptr;
    TSpline3* sp_p0   = nullptr;
    bool hasAccSpline = (goodAcc.size() >= 2);

    if (hasAccSpline) {
        std::vector<double> xpa, ynorm, yk, yp0;
        for (const auto& sp : goodAcc) {
            xpa.push_back(sp.P1); ynorm.push_back(sp.acc_norm);
            yk.push_back(sp.acc_k); yp0.push_back(sp.acc_p0);
        }
        sp_norm = new TSpline3("sp_acc_norm", xpa.data(), ynorm.data(), (int)xpa.size());
        sp_k    = new TSpline3("sp_acc_k",    xpa.data(), yk.data(),    (int)xpa.size());
        sp_p0   = new TSpline3("sp_acc_p0",   xpa.data(), yp0.data(),   (int)xpa.size());
        sp_norm->SetTitle("Spline: acc_norm vs P1;P1 [cm];norm");
        sp_k->SetTitle("Spline: acc_k vs P1;P1 [cm];k");
        sp_p0->SetTitle("Spline: acc_p0 vs P1;P1 [cm];p0");
    }

    // -------------------------------------------------------------------------
    // Step 5: write ROOT file with TTree (raw params) + TSpline3 objects
    // -------------------------------------------------------------------------
    std::string rootPath = OUTPUT_DIR + "/pizza_fit_params.root";
    TFile* fOut = TFile::Open(rootPath.c_str(), "RECREATE");

    // TTree of raw parameters
    TTree* tParams = new TTree("fitParams", "Pizza position fit parameters");
    char   t_label[512]; char t_archive[256];
    double t_P1, t_P2;
    double t_a_res, t_a_res_err, t_b_res, t_b_res_err, t_res_chi2ndf;
    int    t_res_ok;
    double t_acc_norm, t_acc_norm_err, t_acc_k, t_acc_k_err, t_acc_p0, t_acc_p0_err, t_acc_chi2ndf;
    int    t_acc_ok;

    tParams->Branch("label",        t_label,       "label/C");
    tParams->Branch("archive",      t_archive,      "archive/C");
    tParams->Branch("P1",           &t_P1,          "P1/D");
    tParams->Branch("P2",           &t_P2,          "P2/D");
    tParams->Branch("a_res",        &t_a_res,       "a_res/D");
    tParams->Branch("a_res_err",    &t_a_res_err,   "a_res_err/D");
    tParams->Branch("b_res",        &t_b_res,       "b_res/D");
    tParams->Branch("b_res_err",    &t_b_res_err,   "b_res_err/D");
    tParams->Branch("res_chi2ndf",  &t_res_chi2ndf, "res_chi2ndf/D");
    tParams->Branch("res_fit_ok",   &t_res_ok,      "res_fit_ok/I");
    tParams->Branch("acc_norm",     &t_acc_norm,    "acc_norm/D");
    tParams->Branch("acc_norm_err", &t_acc_norm_err,"acc_norm_err/D");
    tParams->Branch("acc_k",        &t_acc_k,       "acc_k/D");
    tParams->Branch("acc_k_err",    &t_acc_k_err,   "acc_k_err/D");
    tParams->Branch("acc_p0",       &t_acc_p0,      "acc_p0/D");
    tParams->Branch("acc_p0_err",   &t_acc_p0_err,  "acc_p0_err/D");
    tParams->Branch("acc_chi2ndf",  &t_acc_chi2ndf, "acc_chi2ndf/D");
    tParams->Branch("acc_fit_ok",   &t_acc_ok,      "acc_fit_ok/I");

    for (const auto& cp : allParams) {
        strncpy(t_label,   cp.label.c_str(),   511);
        strncpy(t_archive, cp.archive.c_str(), 255);
        t_P1 = cp.P1; t_P2 = cp.P2;
        t_a_res = cp.a_res; t_a_res_err = cp.a_res_err;
        t_b_res = cp.b_res; t_b_res_err = cp.b_res_err;
        t_res_chi2ndf = cp.res_fit_chi2ndf;
        t_res_ok = (int)cp.res_fit_ok;
        t_acc_norm = cp.acc_norm; t_acc_norm_err = cp.acc_norm_err;
        t_acc_k = cp.acc_k;       t_acc_k_err = cp.acc_k_err;
        t_acc_p0 = cp.acc_p0;     t_acc_p0_err = cp.acc_p0_err;
        t_acc_chi2ndf = cp.acc_fit_chi2ndf;
        t_acc_ok = (int)cp.acc_fit_ok;
        tParams->Fill();
    }
    tParams->Write();

    // Write splines
    sp_a->Write("sp_a_res");
    sp_b->Write("sp_b_res");
    if (hasAccSpline) {
        sp_norm->Write("sp_acc_norm");
        sp_k->Write("sp_acc_k");
        sp_p0->Write("sp_acc_p0");
    }

    fOut->Close();
    std::cout << "ROOT file written: " << rootPath << std::endl;

    // -------------------------------------------------------------------------
    // Step 6: diagnostic plot — spline curves over data points
    // -------------------------------------------------------------------------
    {
        TCanvas* cSplines = new TCanvas("cSplines", "Parameter Splines", 1400, 500);
        cSplines->Divide(2, 1);

        // -- a_res --
        cSplines->cd(1);
        gPad->SetGrid(1,1);
        std::vector<double> xp_a, ya_a, yb_a;
        for (const auto& sp : goodRes) { xp_a.push_back(sp.P1); ya_a.push_back(sp.a_res); yb_a.push_back(sp.b_res); }
        TGraph* gra = new TGraph((int)xp_a.size(), xp_a.data(), ya_a.data());
        gra->SetMarkerStyle(20); gra->SetMarkerColor(kBlue+1); gra->SetMarkerSize(1.2);
        gra->SetTitle("a_{res} vs Pizza Position;P1 [cm];a_{res}");
        gra->Draw("AP");
        sp_a->SetLineColor(kBlue+1); sp_a->SetLineWidth(2);
        sp_a->Draw("same");

        TGraph* grb = new TGraph((int)xp_a.size(), xp_a.data(), yb_a.data());
        grb->SetMarkerStyle(21); grb->SetMarkerColor(kRed+1); grb->SetMarkerSize(1.2);
        grb->SetTitle("b_{res} vs Pizza Position;P1 [cm];b_{res}");
        grb->Draw("P same");
        sp_b->SetLineColor(kRed+1); sp_b->SetLineWidth(2);
        sp_b->Draw("same");

        TLegend* leg1 = new TLegend(0.55, 0.70, 0.88, 0.88);
        leg1->AddEntry(gra, "a (MS term)", "P");
        leg1->AddEntry(sp_a, "Spline a", "L");
        leg1->AddEntry(grb, "b (hit-res term)", "P");
        leg1->AddEntry(sp_b, "Spline b", "L");
        leg1->Draw();

        // -- acceptance plateau (if available) --
        cSplines->cd(2);
        gPad->SetGrid(1,1);
        if (hasAccSpline) {
            std::vector<double> xpa2, ynorm2;
            for (const auto& sp : goodAcc) { xpa2.push_back(sp.P1); ynorm2.push_back(sp.acc_norm); }
            TGraph* gran = new TGraph((int)xpa2.size(), xpa2.data(), ynorm2.data());
            gran->SetMarkerStyle(20); gran->SetMarkerColor(kGreen+2); gran->SetMarkerSize(1.2);
            gran->SetTitle("Acceptance plateau vs Pizza Position;P1 [cm];Acceptance plateau");
            gran->GetYaxis()->SetRangeUser(0, 1.0);
            gran->Draw("AP");
            sp_norm->SetLineColor(kGreen+2); sp_norm->SetLineWidth(2);
            sp_norm->Draw("same");
        } else {
            TText* tt = new TText(0.3, 0.5, "No acceptance data available");
            tt->SetNDC(); tt->Draw();
        }

        cSplines->SaveAs((OUTPUT_DIR + "/pizza_splines_diagnostic.png").c_str());
        std::cout << "Diagnostic spline plot saved." << std::endl;
        delete cSplines;
    }

    // -------------------------------------------------------------------------
    // Step 7: demo query at a few pizza positions
    // -------------------------------------------------------------------------
    if (hasAccSpline) {
        double P1_min = goodRes.front().P1;
        double P1_max = goodRes.back().P1;
        double P1_mid = 0.5 * (P1_min + P1_max);

        for (double qP1 : {P1_min, P1_mid, P1_max}) {
            queryPizzaPosition(qP1, sp_a, sp_b, sp_norm, sp_k, sp_p0, OUTPUT_DIR);
        }
    } else {
        // resolution-only queries
        double P1_min = goodRes.front().P1;
        double P1_max = goodRes.back().P1;
        for (double qP1 : {P1_min, 0.5*(P1_min+P1_max), P1_max}) {
            double a = sp_a->Eval(qP1);
            double b = sp_b->Eval(qP1);
            std::cout << "\nQuery P1=" << qP1
                      << "  -> a=" << a << "  b=" << b
                      << "  sigma(3GeV)=" << std::sqrt(a*a+(b*3)*(b*3)) << std::endl;
        }
    }

    std::cout << "\n=== KLong_interpolate_pizza complete ===" << std::endl;
    std::cout << "Outputs in: " << OUTPUT_DIR << std::endl;
}

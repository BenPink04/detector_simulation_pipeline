// ============================================================================
// KLong_optimise_geometry.C
//
// PURPOSE:
//   Scan over pizza position P1 using the same spline interpolation
//   infrastructure as KLong_predict_from_geometry.C and identify the
//   "ideal" detector setup — i.e. the P1 that simultaneously maximises
//   acceptance and minimises momentum resolution.
//
// USAGE (via run_optimise.sh, or directly):
//   root -l -b -q 'KLong_optimise_geometry.C+("/path/to/archives_parent")'
//   root -l -b -q 'KLong_optimise_geometry.C+("/path/to/archives_parent","upstream")'
//   root -l -b -q 'KLong_optimise_geometry.C+("/path/to/archives_parent","both","/output/dir")'
//
// ARGUMENTS:
//   parentDir   — directory containing one or more archive sub-directories
//   t3Mode      — "upstream" | "downstream" | "both" (default "both")
//   outputDir   — where to write outputs (default: parentDir/plots)
//
// METRICS:
//   Resolution : mean FWHM integrated over 1–4 GeV/c  (minimise)
//                  = [a0*3 + (exp(a1+4*a2) - exp(a1+a2)) / a2] / 3
//   Acceptance : peak plateau N from the double-logistic fit  (maximise)
//   FOM        : peakN / meanFWHM  (maximise — composite optimum)
//
// OUTPUTS per mode:
//   KLong_optimise_<mode>.png   — three-panel scan plot
//   KLong_optimise_<mode>.root  — canvas + TGraph objects
//   Console summary of optimal P1 values for each metric
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TColor.h"
#include "TROOT.h"
#include "TFitResult.h"
#include "TPad.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>

// ============================================================================
// GEOMETRY HELPERS  (identical to KLong_predict_from_geometry.C)
// ============================================================================

static double opt_extractPos(const std::string& label, const std::string& det) {
    std::string srch = det + "-";
    size_t pos = label.find(srch);
    if (pos == std::string::npos) return -1.0;
    pos += srch.size();
    size_t end = label.find_first_not_of("0123456789", pos);
    return std::stod(label.substr(pos, end - pos));
}

static std::string opt_t3Mode(const std::string& label) {
    double t3 = opt_extractPos(label, "T3");
    if (t3 > 0. && t3 < 660.) return "upstream";
    if (t3 >= 660.)            return "downstream";
    return "no_T3";
}

// ============================================================================
// SPLINE INTERPOLANT  (identical to KLong_predict_from_geometry.C)
// ============================================================================

struct OptParamSpline {
    std::vector<double> xs, ys;
    TSpline3* sp = nullptr;

    void build(const std::string& name) {
        if (xs.size() < 2) { sp = nullptr; return; }
        std::vector<int> idx(xs.size());
        for (int i = 0; i < (int)idx.size(); i++) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b){ return xs[a] < xs[b]; });
        std::vector<double> sx, sy;
        for (int i : idx) { sx.push_back(xs[i]); sy.push_back(ys[i]); }
        xs = sx; ys = sy;
        sp = new TSpline3(name.c_str(), &xs[0], &ys[0], (int)xs.size());
    }

    double eval(double x) const {
        if (!sp || xs.empty()) return 0.;
        int n = (int)xs.size();
        if (x < xs[0]) {
            if (n < 2) return ys[0];
            double slope = (ys[1]-ys[0])/(xs[1]-xs[0]);
            return ys[0] + slope*(x - xs[0]);
        }
        if (x > xs[n-1]) {
            if (n < 2) return ys[n-1];
            double slope = (ys[n-1]-ys[n-2])/(xs[n-1]-xs[n-2]);
            return ys[n-1] + slope*(x - xs[n-1]);
        }
        return sp->Eval(x);
    }

    bool isExtrapolating(double x) const {
        if (xs.empty()) return true;
        return (x < xs[0] || x > xs.back());
    }
    double xMin() const { return xs.empty() ? 0. : xs[0]; }
    double xMax() const { return xs.empty() ? 0. : xs.back(); }
    bool   empty() const { return xs.empty(); }
};

// ============================================================================
// PER-CONFIG FITTING  (identical to KLong_predict_from_geometry.C)
// ============================================================================

struct OptResFitResult {
    double p1, t1;
    std::string mode;
    double a0, a1, a2, chi2ndf;
    double a0e, a1e, a2e;   // parameter errors (diagonal)
    double cov_a0a1, cov_a0a2, cov_a1a2;  // off-diagonal covariances
    bool ok;
};

struct OptAccFitResult {
    double p1, t1;
    std::string mode;
    double N, k, p0, m, p1param, chi2ndf;
    double Ne;              // error on plateau parameter N
    bool ok;
};

OptResFitResult opt_fitResolution(const std::string& vecFile,
                                   const std::string& label) {
    OptResFitResult r;
    r.p1 = opt_extractPos(label,"P1");
    r.t1 = opt_extractPos(label,"T1");
    r.mode = opt_t3Mode(label);
    r.ok = false;
    r.a0 = r.a1 = r.a2 = r.chi2ndf = 0.;
    r.a0e = r.a1e = r.a2e = 0.;
    r.cov_a0a1 = r.cov_a0a2 = r.cov_a1a2 = 0.;

    const int NB = 18;
    double pB[NB+1] = {0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,
                       3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,10.0};
    const double ACUT = 1.0;
    const int  MINENT = 30;
    const double FIT_PMIN = 0.0, FIT_PMAX = 5.0;

    TFile* tf = TFile::Open(vecFile.c_str(),"READ");
    if (!tf || tf->IsZombie()) { if(tf) delete tf; return r; }
    TTree* tree = (TTree*)tf->Get("kaonVectors");
    if (!tree) { tf->Close(); delete tf; return r; }

    std::vector<double>* reco_p = nullptr;
    std::vector<double>* true_p = nullptr;
    tree->SetBranchAddress("reco_p",&reco_p);
    tree->SetBranchAddress("true_p",&true_p);

    std::vector<TH1D*> hb(NB);
    for (int i = 0; i < NB; i++) {
        hb[i] = new TH1D(Form("opt_rh_%s_%d",label.c_str(),i),"",100,-3.,3.);
        hb[i]->SetDirectory(nullptr);
    }
    Long64_t nEnt = tree->GetEntries();
    for (Long64_t ev = 0; ev < nEnt; ev++) {
        tree->GetEntry(ev);
        if (!reco_p||!true_p) continue;
        for (int k = 0; k < (int)reco_p->size(); k++) {
            double tp = (*true_p)[k], dp = (*reco_p)[k]-tp;
            if (tp<=0.) continue;
            if (std::abs(dp)/tp > ACUT) continue;
            for (int b = 0; b < NB; b++)
                if (tp>=pB[b]&&tp<pB[b+1]) { hb[b]->Fill(dp); break; }
        }
    }
    tf->Close(); delete tf;

    std::vector<double> gx, gy, gex, gey;
    for (int i = 0; i < NB; i++) {
        if (hb[i]->GetEntries()<MINENT) { delete hb[i]; continue; }
        double mean = hb[i]->GetMean(), rms = hb[i]->GetRMS();
        TF1* gf = new TF1(Form("opt_rgf_%s_%d",label.c_str(),i),"gaus",
                          mean-2.*rms, mean+2.*rms);
        gf->SetParameters(hb[i]->GetMaximum(), mean, rms);
        hb[i]->Fit(gf,"RQN");
        double sig = std::abs(gf->GetParameter(2));
        double sige = gf->GetParError(2);
        double mid = 0.5*(pB[i]+pB[i+1]);
        if (mid >= FIT_PMIN && mid <= FIT_PMAX) {
            const double MIN_REL_ERR = 0.20;
            double fwhm   = 2.3548 * sig;
            double fwhm_e = 2.3548 * (sige > 0. ? sige : 1e-6);
            if (fwhm_e < MIN_REL_ERR * fwhm) fwhm_e = MIN_REL_ERR * fwhm;
            gx.push_back(mid); gy.push_back(fwhm);
            gex.push_back(0.5*(pB[i+1]-pB[i])); gey.push_back(fwhm_e);
        }
        hb[i]->GetListOfFunctions()->Remove(gf);
        delete gf; delete hb[i];
    }
    if ((int)gx.size() < 3) return r;

    int nb = (int)gx.size();
    std::vector<double> edges;
    for (int i = 0; i < nb; i++) edges.push_back(gx[i]-gex[i]);
    edges.push_back(gx[nb-1]+gex[nb-1]);
    TH1F* hd = new TH1F(Form("opt_hd_%s",label.c_str()),"",nb,edges.data());
    hd->SetDirectory(nullptr);
    for (int i = 0; i < nb; i++) { hd->SetBinContent(i+1,gy[i]); hd->SetBinError(i+1,gey[i]); }

    TF1* fm = new TF1(Form("opt_rfm_%s",label.c_str()),"pol0(0)+expo(1)",FIT_PMIN,FIT_PMAX);
    fm->SetParameter(0, 0.05);
    fm->SetParameter(1, -0.5);
    fm->SetParameter(2, -0.5);
    TFitResultPtr fr = hd->Fit(fm,"SRLN");
    hd->GetListOfFunctions()->Remove(fm);
    delete hd;

    if ((int)fr==0 || fr->IsValid()) {
        r.a0=fm->GetParameter(0); r.a1=fm->GetParameter(1); r.a2=fm->GetParameter(2);
        r.a0e=fm->GetParError(0); r.a1e=fm->GetParError(1); r.a2e=fm->GetParError(2);
        if (fr->IsValid()) {
            r.cov_a0a1 = fr->CovMatrix(0,1);
            r.cov_a0a2 = fr->CovMatrix(0,2);
            r.cov_a1a2 = fr->CovMatrix(1,2);
        }
        double ndf = fm->GetNDF();
        r.chi2ndf = (ndf>0) ? fm->GetChisquare()/ndf : -1.;
        r.ok = true;
    }
    delete fm;
    return r;
}

OptAccFitResult opt_fitAcceptance(const std::string& accFile,
                                   const std::string& label) {
    OptAccFitResult r;
    r.p1 = opt_extractPos(label,"P1");
    r.t1 = opt_extractPos(label,"T1");
    r.mode = opt_t3Mode(label);
    r.ok = false;
    r.N=r.k=r.p0=r.m=r.p1param=r.chi2ndf=0.;
    r.Ne=0.;

    const double P_MIN=0., P_MAX=10., BIN_W=0.25;
    const int NBINS = (int)((P_MAX-P_MIN)/BIN_W);

    TFile* tf = TFile::Open(accFile.c_str(),"READ");
    if (!tf||tf->IsZombie()) { if(tf)delete tf; return r; }
    TTree* tree = (TTree*)tf->Get("kaonEventInfo");
    if (!tree) { tf->Close(); delete tf; return r; }

    std::vector<double>* true_mom_vec = nullptr;
    std::vector<int>*    reco_flag_vec = nullptr;
    int n_triple = 0;
    tree->SetBranchAddress("true_mom_vec",&true_mom_vec);
    tree->SetBranchAddress("reco_flag_vec",&reco_flag_vec);
    tree->SetBranchAddress("n_triple_pion_events",&n_triple);

    std::vector<double> totalCnt(NBINS,0.), recoCnt(NBINS,0.);
    Long64_t nEnt = tree->GetEntries(); double totalTriple = 0.;
    for (Long64_t ev = 0; ev < nEnt; ev++) {
        tree->GetEntry(ev);
        totalTriple += n_triple;
        if (!true_mom_vec||!reco_flag_vec) continue;
        for (int k = 0; k < (int)true_mom_vec->size(); k++) {
            double p = (*true_mom_vec)[k];
            int bn = (int)((p-P_MIN)/BIN_W);
            if (bn<0||bn>=NBINS) continue;
            totalCnt[bn]++;
            if ((*reco_flag_vec)[k]==1) recoCnt[bn]++;
        }
    }
    tf->Close(); delete tf;
    if (totalTriple <= 0.) return r;

    std::vector<double> bx, by, bey;
    for (int i = 0; i < NBINS; i++) {
        if (totalCnt[i]<5.) continue;
        double eff = recoCnt[i]/totalTriple;
        double mid = P_MIN+(i+0.5)*BIN_W;
        bx.push_back(mid); by.push_back(eff);
        double err = std::sqrt(eff*(1.-std::min(eff,1.))/totalCnt[i]);
        bey.push_back(err>0.?err:1e-6);
    }
    if ((int)bx.size() < 6) return r;

    std::vector<double> bez(bx.size(),0.);
    TGraphErrors* gg = new TGraphErrors((int)bx.size(),&bx[0],&by[0],&bez[0],&bey[0]);
    TF1* fm = new TF1(Form("opt_afm_%s",label.c_str()),
        "[0]/((1.+exp(-[1]*(x-[2])))*(1.+exp([3]*(x-[4]))))",0.5,10.);
    fm->SetParNames("N","k","p0","m","p1");
    double ymax = *std::max_element(by.begin(),by.end());
    fm->SetParameters(ymax,4.,1.1,0.2,5.);
    fm->SetParLimits(0,0.,1.); fm->SetParLimits(1,0.5,20.);
    fm->SetParLimits(2,0.5,2.5); fm->SetParLimits(3,0.,2.);
    fm->SetParLimits(4,1.,9.);
    TFitResultPtr fr = gg->Fit(fm,"QSRN");
    delete gg;

    if ((int)fr==0||fr->IsValid()) {
        r.N=fm->GetParameter(0); r.k=fm->GetParameter(1);
        r.p0=fm->GetParameter(2); r.m=fm->GetParameter(3);
        r.p1param=fm->GetParameter(4);
        r.Ne=fm->GetParError(0);
        double ndf = fm->GetNDF();
        r.chi2ndf = (ndf>0) ? fm->GetChisquare()/ndf : -1.;
        r.ok = true;
    }
    delete fm;
    return r;
}

// ============================================================================
// HELPER: compute mean FWHM over 1–4 GeV/c given pol0+expo parameters
//   FWHM(p) = a0 + exp(a1 + a2*p)
//   mean = (1/3) * integral_1^4 FWHM(p) dp
//         = a0 + (exp(a1+4*a2) - exp(a1+a2)) / (3*a2)
// ============================================================================
static double meanFWHM_1to4(double a0, double a1, double a2) {
    if (std::abs(a2) < 1e-10)
        return a0 + std::exp(a1);        // degenerate: expo is constant
    return a0 + (std::exp(a1 + 4.*a2) - std::exp(a1 + 1.*a2)) / (3.*a2);
}

// ============================================================================
// HELPER: run one mode scan, fill result vectors, print summary
// ============================================================================
struct ScanResult {
    std::vector<double> p1, fwhm, acc, fom, fom_err;
    double optP1_fwhm, optP1_acc, optP1_fom;
    double optVal_fwhm, optVal_acc, optVal_fom;
    double p1Min, p1Max;  // interpolation range
    bool valid = false;
};

// Simple symmetric moving-average smoother. halfWin is in number of points.
static std::vector<double> smoothVector(const std::vector<double>& v, int halfWin) {
    int n = (int)v.size();
    std::vector<double> s(n);
    for (int i = 0; i < n; i++) {
        int lo = std::max(0, i - halfWin);
        int hi = std::min(n - 1, i + halfWin);
        double sum = 0.; int cnt = 0;
        for (int j = lo; j <= hi; j++) { sum += v[j]; cnt++; }
        s[i] = sum / cnt;
    }
    return s;
}

static ScanResult runScan(const std::string& mode,
                           const std::vector<OptResFitResult>& resResults,
                           const std::vector<OptAccFitResult>& accResults,
                           double scanMin, double scanMax, double scanStep,
                           int fomMode = 1) {
    ScanResult out;

    // ---- build splines ----
    OptParamSpline sp_a0, sp_a1, sp_a2, spN, sp_k, sp_p0, sp_m, sp_p1p;
    OptParamSpline sp_a0e, sp_a1e, sp_a2e, spNe;  // diagonal errors
    OptParamSpline sp_cov01, sp_cov02, sp_cov12;   // off-diagonal covariances
    {
        std::map<double,double> ma0, ma1, ma2, ma0e, ma1e, ma2e;
        std::map<double,double> mcov01, mcov02, mcov12;
        for (auto& r : resResults) {
            if (r.mode != mode) continue;
            ma0[r.p1]=r.a0; ma1[r.p1]=r.a1; ma2[r.p1]=r.a2;
            ma0e[r.p1]=r.a0e; ma1e[r.p1]=r.a1e; ma2e[r.p1]=r.a2e;
            mcov01[r.p1]=r.cov_a0a1; mcov02[r.p1]=r.cov_a0a2; mcov12[r.p1]=r.cov_a1a2;
        }
        for (auto& [p,v] : ma0)    { sp_a0.xs.push_back(p);   sp_a0.ys.push_back(v); }
        for (auto& [p,v] : ma1)    { sp_a1.xs.push_back(p);   sp_a1.ys.push_back(v); }
        for (auto& [p,v] : ma2)    { sp_a2.xs.push_back(p);   sp_a2.ys.push_back(v); }
        for (auto& [p,v] : ma0e)   { sp_a0e.xs.push_back(p);  sp_a0e.ys.push_back(std::abs(v)); }
        for (auto& [p,v] : ma1e)   { sp_a1e.xs.push_back(p);  sp_a1e.ys.push_back(std::abs(v)); }
        for (auto& [p,v] : ma2e)   { sp_a2e.xs.push_back(p);  sp_a2e.ys.push_back(std::abs(v)); }
        for (auto& [p,v] : mcov01) { sp_cov01.xs.push_back(p); sp_cov01.ys.push_back(v); }
        for (auto& [p,v] : mcov02) { sp_cov02.xs.push_back(p); sp_cov02.ys.push_back(v); }
        for (auto& [p,v] : mcov12) { sp_cov12.xs.push_back(p); sp_cov12.ys.push_back(v); }
    }
    {
        std::map<double,double> mN, mk, mp0, mm, mp1, mNe;
        for (auto& r : accResults) {
            if (r.mode != mode) continue;
            mN[r.p1]=r.N; mk[r.p1]=r.k; mp0[r.p1]=r.p0;
            mm[r.p1]=r.m; mp1[r.p1]=r.p1param;
            mNe[r.p1]=r.Ne;
        }
        for (auto& [p,v]:mN)  { spN.xs.push_back(p);   spN.ys.push_back(v); }
        for (auto& [p,v]:mk)  { sp_k.xs.push_back(p);  sp_k.ys.push_back(v); }
        for (auto& [p,v]:mp0) { sp_p0.xs.push_back(p); sp_p0.ys.push_back(v); }
        for (auto& [p,v]:mm)  { sp_m.xs.push_back(p);  sp_m.ys.push_back(v); }
        for (auto& [p,v]:mp1) { sp_p1p.xs.push_back(p);sp_p1p.ys.push_back(v); }
        for (auto& [p,v]:mNe) { spNe.xs.push_back(p);  spNe.ys.push_back(std::abs(v)); }
    }

    if (sp_a0.empty() || spN.empty()) {
        std::cerr << "  [" << mode << "] Not enough training data — skipping." << std::endl;
        return out;
    }

    sp_a0.build("opt_sp_a0_"+mode); sp_a1.build("opt_sp_a1_"+mode);
    sp_a2.build("opt_sp_a2_"+mode);
    sp_a0e.build("opt_sp_a0e_"+mode); sp_a1e.build("opt_sp_a1e_"+mode);
    sp_a2e.build("opt_sp_a2e_"+mode);
    sp_cov01.build("opt_sp_cov01_"+mode); sp_cov02.build("opt_sp_cov02_"+mode);
    sp_cov12.build("opt_sp_cov12_"+mode);
    spN.build("opt_spN_"+mode); spNe.build("opt_spNe_"+mode);
    sp_k.build("opt_sp_k_"+mode); sp_p0.build("opt_sp_p0_"+mode);
    sp_m.build("opt_sp_m_"+mode); sp_p1p.build("opt_sp_p1p_"+mode);

    // Use the interpolation range (no extrapolation in the scan)
    double lo = std::max(scanMin, std::max(sp_a0.xMin(), spN.xMin()));
    double hi = std::min(scanMax, std::min(sp_a0.xMax(), spN.xMax()));
    out.p1Min = lo; out.p1Max = hi;

    if (hi <= lo) {
        std::cerr << "  [" << mode << "] Scan range empty after clamping to training data." << std::endl;
        return out;
    }

    // ---- scan ----
    double bestFWHM_val = 1e9, bestAcc_val = -1e9, bestFOM_val = -1e9;
    double bestFWHM_p1 = lo, bestAcc_p1 = lo, bestFOM_p1 = lo;

    for (double p1 = lo; p1 <= hi + 1e-9; p1 += scanStep) {
        double a0 = sp_a0.eval(p1);
        double a1 = sp_a1.eval(p1);
        double a2 = sp_a2.eval(p1);
        double mf = meanFWHM_1to4(a0, a1, a2);
        // Clamp: negative FWHM is unphysical (can happen at extrapolation edges)
        if (mf < 0.) mf = 0.;

        double peakAcc = std::max(std::min(spN.eval(p1), 1.), 0.);
        double denom = (fomMode == 2) ? mf * mf : mf;
        double fom = (denom > 1e-10) ? peakAcc / denom : 0.;

        // ---- error propagation (including covariance cross-terms) ----
        double sa0v = sp_a0e.eval(p1), sa1v = sp_a1e.eval(p1), sa2v = sp_a2e.eval(p1);
        double cov01 = sp_cov01.eval(p1), cov02 = sp_cov02.eval(p1), cov12 = sp_cov12.eval(p1);
        double sigma_mf = 0.;
        if (mf > 0.) {
            double E4 = std::exp(a1 + 4.*a2), E1 = std::exp(a1 + a2);
            double d0 = 1.0;
            double d1 = (std::abs(a2) > 1e-10) ? (E4 - E1) / (3.*a2) : std::exp(a1);
            double d2 = (std::abs(a2) > 1e-10)
                ? (4.*E4 - E1) / (3.*a2) - (E4 - E1) / (3.*a2*a2) : 0.;
            double var = d0*d0*sa0v*sa0v + d1*d1*sa1v*sa1v + d2*d2*sa2v*sa2v
                       + 2.*d0*d1*cov01 + 2.*d0*d2*cov02 + 2.*d1*d2*cov12;
            sigma_mf = (var > 0.) ? std::sqrt(var) : 0.;
        }
        double sigma_acc = spNe.eval(p1);
        double sigma_fom = (mf > 1e-10 && peakAcc > 1e-10)
            ? fom * std::sqrt(std::pow(sigma_acc / peakAcc, 2)
                            + std::pow((fomMode==2 ? 2. : 1.) * sigma_mf / mf, 2))
            : 0.;
        // Clamp to 50% relative — spline edges overshoot at sparse training points
        sigma_fom = std::min(sigma_fom, 0.1 * fom);

        out.p1.push_back(p1);
        out.fwhm.push_back(mf);
        out.acc.push_back(peakAcc);
        out.fom.push_back(fom);
        out.fom_err.push_back(sigma_fom);

        if (mf < bestFWHM_val) { bestFWHM_val = mf; bestFWHM_p1 = p1; }
        if (peakAcc > bestAcc_val) { bestAcc_val = peakAcc; bestAcc_p1 = p1; }
        if (fom > bestFOM_val) { bestFOM_val = fom; bestFOM_p1 = p1; }
    }

    out.optP1_fwhm = bestFWHM_p1; out.optVal_fwhm = bestFWHM_val;
    out.optP1_acc  = bestAcc_p1;  out.optVal_acc  = bestAcc_val;
    out.optP1_fom  = bestFOM_p1;  out.optVal_fom  = bestFOM_val;

    // Smooth all output curves over a ±15 cm window (15 pts at 1 cm/step)
    // to suppress TSpline3 oscillations between training knots.
    const int HALF_WIN = 15;
    out.fwhm    = smoothVector(out.fwhm,    HALF_WIN);
    out.acc     = smoothVector(out.acc,     HALF_WIN);
    out.fom     = smoothVector(out.fom,     HALF_WIN);
    out.fom_err = smoothVector(out.fom_err, HALF_WIN);

    // Re-derive optimal positions from smoothed curves
    out.optP1_fwhm = out.p1[0]; out.optVal_fwhm = out.fwhm[0];
    out.optP1_acc  = out.p1[0]; out.optVal_acc  = out.acc[0];
    out.optP1_fom  = out.p1[0]; out.optVal_fom  = out.fom[0];
    for (int i = 1; i < (int)out.p1.size(); i++) {
        if (out.fwhm[i] < out.optVal_fwhm) { out.optVal_fwhm = out.fwhm[i]; out.optP1_fwhm = out.p1[i]; }
        if (out.acc[i]  > out.optVal_acc)  { out.optVal_acc  = out.acc[i];  out.optP1_acc  = out.p1[i]; }
        if (out.fom[i]  > out.optVal_fom)  { out.optVal_fom  = out.fom[i];  out.optP1_fom  = out.p1[i]; }
    }

    out.valid = true;

    std::cout << "\n  ┌─ Optimisation results [" << mode << "] "
              << "(P1 range: " << lo << "–" << hi << " cm, step=" << scanStep << " cm)" << std::endl;
    std::cout << "  │  Best RESOLUTION  → P1 = " << bestFWHM_p1
              << " cm   (mean FWHM 1–4 GeV/c = " << bestFWHM_val << " GeV/c)" << std::endl;
    std::cout << "  │  Best ACCEPTANCE  → P1 = " << bestAcc_p1
              << " cm   (peak N = " << bestAcc_val << ")" << std::endl;
    std::cout << "  └  Best FOM (N/FWHM) → P1 = " << bestFOM_p1
              << " cm   (FOM = " << bestFOM_val << ")" << std::endl;

    return out;
}

// ============================================================================
// HELPER: draw a single mode's scan onto three existing TPads
// ============================================================================
static void drawScanPanels(TPad* pFWHM, TPad* pAcc, TPad* pFOM,
                            const ScanResult& s, const std::string& mode,
                            int colFWHM, int colAcc, int colFOM,
                            bool isFirst, int fomMode = 1) {
    if (!s.valid || s.p1.empty()) return;
    int n = (int)s.p1.size();

    auto makeGraph = [&](const std::vector<double>& y) -> TGraph* {
        return new TGraph(n, s.p1.data(), y.data());
    };

    // --- FWHM panel ---
    pFWHM->cd();
    TGraph* gF = makeGraph(s.fwhm);
    gF->SetName(Form("gFWHM_%s",mode.c_str()));
    gF->SetTitle("Mean FWHM vs Pizza Position;Pizza position P1 [cm];Mean FWHM 1#font[122]{-}4 GeV/c [GeV/c]");
    gF->SetLineColor(colFWHM); gF->SetLineWidth(2); gF->SetLineStyle(1);
    gF->Draw(isFirst ? "AL" : "L same");
    if (isFirst) {
        gF->GetXaxis()->SetTitleSize(0.055); gF->GetYaxis()->SetTitleSize(0.052);
        gF->GetXaxis()->SetLabelSize(0.048); gF->GetYaxis()->SetLabelSize(0.048);
        gF->GetYaxis()->SetTitleOffset(1.9);
        gF->GetXaxis()->SetNdivisions(5);
        gF->GetYaxis()->SetNdivisions(5);
    }
    // Compute axis range with headroom
    double yMinF = *std::min_element(s.fwhm.begin(),s.fwhm.end()) * 0.90;
    double yMaxF = *std::max_element(s.fwhm.begin(),s.fwhm.end()) * 1.15;
    gF->GetHistogram()->SetMinimum(yMinF);
    gF->GetHistogram()->SetMaximum(yMaxF);
    gPad->Update();
    // Vertical line spanning full frame height
    TLine* lF = new TLine(s.optP1_fwhm, gPad->GetUymin(), s.optP1_fwhm, gPad->GetUymax());
    lF->SetLineColor(colFWHM); lF->SetLineStyle(2); lF->SetLineWidth(2);
    lF->Draw();

    // --- Acceptance panel ---
    pAcc->cd();
    TGraph* gA = makeGraph(s.acc);
    gA->SetName(Form("gAcc_%s",mode.c_str()));
    gA->SetTitle("Peak Acceptance vs Pizza Position;Pizza position P1 [cm];Peak acceptance N");
    gA->SetLineColor(colAcc); gA->SetLineWidth(2); gA->SetLineStyle(1);
    gA->Draw(isFirst ? "AL" : "L same");
    if (isFirst) {
        gA->GetXaxis()->SetTitleSize(0.055); gA->GetYaxis()->SetTitleSize(0.052);
        gA->GetXaxis()->SetLabelSize(0.048); gA->GetYaxis()->SetLabelSize(0.048);
        gA->GetYaxis()->SetTitleOffset(1.9);
        gA->GetXaxis()->SetNdivisions(5);
        gA->GetYaxis()->SetNdivisions(5);
    }
    double yMinA = *std::min_element(s.acc.begin(),s.acc.end()) * 0.90;
    double yMaxA = *std::max_element(s.acc.begin(),s.acc.end()) * 1.15;
    gA->GetHistogram()->SetMinimum(yMinA);
    gA->GetHistogram()->SetMaximum(yMaxA);
    gPad->Update();
    TLine* lA = new TLine(s.optP1_acc, gPad->GetUymin(), s.optP1_acc, gPad->GetUymax());
    lA->SetLineColor(colAcc); lA->SetLineStyle(2); lA->SetLineWidth(2);
    lA->Draw();

    // --- FOM panel ---
    pFOM->cd();
    int nf = (int)s.p1.size();
    // Compute axis range accounting for error bars
    double yMinFOM = 1e30, yMaxFOM = -1e30;
    for (int i = 0; i < nf; i++) {
        yMinFOM = std::min(yMinFOM, s.fom[i] - s.fom_err[i]);
        yMaxFOM = std::max(yMaxFOM, s.fom[i] + s.fom_err[i]);
    }
    yMinFOM = std::max(0., yMinFOM * 0.90);
    yMaxFOM *= 1.15;
    TGraphErrors* gFOMband = new TGraphErrors(nf, s.p1.data(), s.fom.data(),
                                               nullptr, s.fom_err.data());
    gFOMband->SetName(Form("gFOMband_%s",mode.c_str()));
    gFOMband->SetTitle(Form("Figure of Merit vs Pizza Position;Pizza position P1 [cm];Figure of merit  N / #LT FWHM%s #GT",
                             fomMode==2 ? "^{2}" : ""));
    gFOMband->SetFillColorAlpha(colFOM, 0.25);
    gFOMband->SetFillStyle(1001);
    gFOMband->SetLineWidth(0);
    gFOMband->Draw(isFirst ? "A E3" : "E3 same");
    if (isFirst) {
        gFOMband->GetXaxis()->SetTitleSize(0.055); gFOMband->GetYaxis()->SetTitleSize(0.050);
        gFOMband->GetXaxis()->SetLabelSize(0.048); gFOMband->GetYaxis()->SetLabelSize(0.048);
        gFOMband->GetYaxis()->SetTitleOffset(1.9);
        gFOMband->GetXaxis()->SetNdivisions(5);
        gFOMband->GetYaxis()->SetNdivisions(5);
    }
    gFOMband->GetHistogram()->SetMinimum(yMinFOM);
    gFOMband->GetHistogram()->SetMaximum(yMaxFOM);
    gPad->Update();
    TGraph* gFOM = new TGraph(nf, s.p1.data(), s.fom.data());
    gFOM->SetName(Form("gFOM_%s",mode.c_str()));
    gFOM->SetLineColor(colFOM); gFOM->SetLineWidth(2); gFOM->SetLineStyle(1);
    gFOM->Draw("L same");
    TLine* lFOM = new TLine(s.optP1_fom, gPad->GetUymin(), s.optP1_fom, gPad->GetUymax());
    lFOM->SetLineColor(colFOM); lFOM->SetLineStyle(2); lFOM->SetLineWidth(2);
    lFOM->Draw();
}

// Draw legend + optimal position labels in the bottom annotation strip
static void drawBottomStrip(TPad* pBot,
    bool hasUp, bool hasDn,
    int colFWHM_up, int colFWHM_dn,
    int colAcc_up,  int colAcc_dn,
    int colFOM_up,  int colFOM_dn,
    double optFWHM_up, double optAcc_up, double optFOM_up,
    double optFWHM_dn, double optAcc_dn, double optFOM_dn) {

    pBot->cd();
    pBot->SetFillStyle(0); pBot->SetBorderSize(0);

    // Column x-centres in NDC (matching the three equal plot pads above)
    const double cx[3]  = {0.167, 0.500, 0.833};
    const char*  hdr[3] = {"Best FWHM position", "Best Acc. position", "Best FOM position"};
    const int    cup[3] = {colFWHM_up, colAcc_up,  colFOM_up};
    const int    cdn[3] = {colFWHM_dn, colAcc_dn,  colFOM_dn};
    const double vup[3] = {optFWHM_up, optAcc_up,  optFOM_up};
    const double vdn[3] = {optFWHM_dn, optAcc_dn,  optFOM_dn};

    for (int i = 0; i < 3; i++) {
        TLatex* th = new TLatex(cx[i], 0.82, hdr[i]);
        th->SetNDC(true); th->SetTextAlign(21);
        th->SetTextSize(0.22); th->SetTextColor(kGray+2); th->Draw();

        double y = 0.58;
        if (hasUp) {
            TLatex* tu = new TLatex(cx[i], y, Form("Upstream: %.0f cm", vup[i]));
            tu->SetNDC(true); tu->SetTextAlign(21);
            tu->SetTextSize(0.26); tu->SetTextColor(cup[i]); tu->Draw();
            y -= 0.28;
        }
        if (hasDn) {
            TLatex* td = new TLatex(cx[i], y, Form("Downstream: %.0f cm", vdn[i]));
            td->SetNDC(true); td->SetTextAlign(21);
            td->SetTextSize(0.26); td->SetTextColor(cdn[i]); td->Draw();
        }
    }
}

// ============================================================================
// MAIN
// ============================================================================

void KLong_optimise_geometry(const char* parentDirC,
                               const char* t3ModeC   = "both",
                               const char* outputDirC = "",
                               int         fomMode   = 1) {

    std::string parentDir(parentDirC);
    while (!parentDir.empty() && parentDir.back()=='/') parentDir.pop_back();

    std::string t3Mode(t3ModeC);
    std::string outputDir(outputDirC);
    if (outputDir.empty()) outputDir = parentDir + "/plots";
    gSystem->mkdir(outputDir.c_str(), kTRUE);

    // Scan settings — stay within the interpolation range
    const double SCAN_MIN  = 100.;   // cm  (clipped to training range automatically)
    const double SCAN_MAX  = 600.;   // cm
    const double SCAN_STEP = 1.0;    // cm

    // -------------------------------------------------------------------------
    // 1. Discover all training files
    // -------------------------------------------------------------------------
    std::cout << "\n=== KLong Geometry Optimisation ===" << std::endl;
    std::cout << "    Source: " << parentDir << std::endl;
    std::cout << "    Mode:   " << t3Mode << std::endl;

    auto findFiles = [&](const std::string& suffix)
        -> std::vector<std::pair<std::string,std::string>> {
        std::string cmd = "find \""+parentDir+"\" -maxdepth 2 -name '*"+suffix
                         +"' 2>/dev/null | sort";
        FILE* pipe = popen(cmd.c_str(),"r");
        std::vector<std::pair<std::string,std::string>> out;
        if (!pipe) return out;
        char buf[2048];
        while (fgets(buf,sizeof(buf),pipe)) {
            std::string fp(buf);
            fp.erase(std::remove(fp.begin(),fp.end(),'\n'),fp.end());
            if (fp.empty()) continue;
            size_t sl = fp.find_last_of('/');
            std::string base = (sl!=std::string::npos)?fp.substr(sl+1):fp;
            size_t sfx = base.find("_"+suffix);
            if (sfx==std::string::npos) sfx=base.find(suffix);
            std::string lbl = (sfx!=std::string::npos)?base.substr(0,sfx):base;
            out.push_back({fp,lbl});
        }
        pclose(pipe);
        return out;
    };

    auto vecFiles = findFiles("combined_vectors.root");
    auto accFiles = findFiles("combined_acceptance.root");
    std::cout << "  Found " << vecFiles.size() << " vector files, "
              << accFiles.size() << " acceptance files." << std::endl;

    if (vecFiles.empty() && accFiles.empty()) {
        std::cerr << "\nERROR: No training files found under " << parentDir << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // 2. Fit all training configurations
    // -------------------------------------------------------------------------
    std::vector<OptResFitResult> resResults;
    std::vector<OptAccFitResult> accResults;
    std::map<std::string,bool> seenVec, seenAcc;

    std::cout << "\n--- Fitting training data ---" << std::endl;
    for (auto& [fp,lbl] : vecFiles) {
        if (seenVec[lbl]) continue;
        seenVec[lbl] = true;
        OptResFitResult r = opt_fitResolution(fp,lbl);
        if (r.ok) {
            resResults.push_back(r);
            std::cout << "  [RES OK] P1=" << r.p1 << "  mode=" << r.mode << std::endl;
        }
    }
    for (auto& [fp,lbl] : accFiles) {
        if (seenAcc[lbl]) continue;
        seenAcc[lbl] = true;
        OptAccFitResult r = opt_fitAcceptance(fp,lbl);
        if (r.ok) {
            accResults.push_back(r);
            std::cout << "  [ACC OK] P1=" << r.p1 << "  N=" << r.N
                      << "  mode=" << r.mode << std::endl;
        }
    }

    if (resResults.empty() || accResults.empty()) {
        std::cerr << "\nERROR: Need both resolution and acceptance fits. "
                  << "Check your archive directory." << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // 3. Determine which modes to run
    // -------------------------------------------------------------------------
    bool doUpstream   = (t3Mode=="upstream"   || t3Mode=="both" || t3Mode=="auto");
    bool doDownstream = (t3Mode=="downstream" || t3Mode=="both" || t3Mode=="auto");

    // Check which modes actually have data
    auto hasModeData = [&](const std::string& m) {
        for (auto& r : resResults) if (r.mode==m && r.ok) return true;
        return false;
    };
    if (!hasModeData("upstream"))   doUpstream   = false;
    if (!hasModeData("downstream")) doDownstream = false;

    if (!doUpstream && !doDownstream) {
        std::cerr << "\nERROR: No data for requested mode(s)." << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // 4. Run scans
    // -------------------------------------------------------------------------
    ScanResult scanUp, scanDown;
    if (doUpstream)   scanUp   = runScan("upstream",   resResults, accResults, SCAN_MIN, SCAN_MAX, SCAN_STEP, fomMode);
    if (doDownstream) scanDown = runScan("downstream", resResults, accResults, SCAN_MIN, SCAN_MAX, SCAN_STEP, fomMode);

    // -------------------------------------------------------------------------
    // 5. Plot — three-panel canvas, one curve per mode
    // -------------------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetFrameLineWidth(1);

    // Per-panel colours: upstream / downstream are distinct complementary colours
    // FWHM panel   : red      / cyan
    // Acceptance   : green    / magenta
    // FOM          : violet   / orange
    const int colFWHM_up = kRed+1;     const int colFWHM_dn = kCyan+1;
    const int colAcc_up  = kGreen+2;   const int colAcc_dn  = kMagenta+1;
    const int colFOM_up  = kViolet+1;  const int colFOM_dn  = kOrange+1;

    TCanvas* c = new TCanvas("cOpt","KLong Geometry Optimisation",1800,750);
    TPad* pFWHM = new TPad("pFWHM","",0.000,0.26,0.333,1.00); pFWHM->Draw();
    TPad* pAcc  = new TPad("pAcc", "",0.333,0.26,0.667,1.00); pAcc->Draw();
    TPad* pFOM  = new TPad("pFOM", "",0.667,0.26,1.000,1.00); pFOM->Draw();
    TPad* pBot  = new TPad("pBot", "",0.000,0.00,1.000,0.26); pBot->Draw();

    for (TPad* p : {pFWHM, pAcc, pFOM}) {
        p->SetGrid(1,1);
        p->SetLeftMargin(0.20); p->SetRightMargin(0.04);
        p->SetTopMargin(0.11);  p->SetBottomMargin(0.14);
    }
    pBot->SetMargin(0.01, 0.01, 0.01, 0.01);

    // Titles are carried by each TGraph's SetTitle() — no separate header needed

    bool firstDrawn = false;
    if (doUpstream && scanUp.valid) {
        drawScanPanels(pFWHM, pAcc, pFOM, scanUp,   "upstream",
                       colFWHM_up, colAcc_up, colFOM_up, !firstDrawn, fomMode);
        firstDrawn = true;
    }
    if (doDownstream && scanDown.valid) {
        drawScanPanels(pFWHM, pAcc, pFOM, scanDown, "downstream",
                       colFWHM_dn, colAcc_dn, colFOM_dn, !firstDrawn, fomMode);
        firstDrawn = true;
    }

    // Bottom strip: legend + optimal position labels
    {
        double up_fwhm = (doUpstream   && scanUp.valid)   ? scanUp.optP1_fwhm   : 0.;
        double up_acc  = (doUpstream   && scanUp.valid)   ? scanUp.optP1_acc    : 0.;
        double up_fom  = (doUpstream   && scanUp.valid)   ? scanUp.optP1_fom    : 0.;
        double dn_fwhm = (doDownstream && scanDown.valid) ? scanDown.optP1_fwhm : 0.;
        double dn_acc  = (doDownstream && scanDown.valid) ? scanDown.optP1_acc  : 0.;
        double dn_fom  = (doDownstream && scanDown.valid) ? scanDown.optP1_fom  : 0.;
        drawBottomStrip(pBot,
            doUpstream && scanUp.valid, doDownstream && scanDown.valid,
            colFWHM_up, colFWHM_dn, colAcc_up, colAcc_dn, colFOM_up, colFOM_dn,
            up_fwhm, up_acc, up_fom, dn_fwhm, dn_acc, dn_fom);
    }

    // ---- Save ----
    std::string suffix = (doUpstream && doDownstream) ? "both"
                       : (doUpstream ? "upstream" : "downstream");
    std::string outPNG  = outputDir + "/KLong_optimise_" + suffix + ".png";
    std::string outROOT = outputDir + "/KLong_optimise_" + suffix + ".root";

    c->SaveAs(outPNG.c_str());
    std::cout << "\nSaved: " << outPNG << std::endl;

    TFile* fOut = TFile::Open(outROOT.c_str(),"RECREATE");
    c->Write();
    fOut->Close(); delete fOut;
    std::cout << "Saved: " << outROOT << std::endl;

    // -------------------------------------------------------------------------
    // 6. Final summary
    // -------------------------------------------------------------------------
    std::cout << "\n╔══════════════════════════════════════════════════════════╗" << std::endl;
    std::cout <<   "║          IDEAL SETUP SUMMARY                             ║" << std::endl;
    std::cout <<   "╠══════════════════════════════════════════════════════════╣" << std::endl;
    if (doUpstream && scanUp.valid) {
        std::cout << "║  UPSTREAM" << std::endl;
        std::cout << "║    Best resolution : P1 = " << std::setw(6) << scanUp.optP1_fwhm
                  << " cm  (mean FWHM = " << scanUp.optVal_fwhm << " GeV/c)" << std::endl;
        std::cout << "║    Best acceptance : P1 = " << std::setw(6) << scanUp.optP1_acc
                  << " cm  (peak N = "   << scanUp.optVal_acc  << ")" << std::endl;
        std::cout << "║    Best FOM        : P1 = " << std::setw(6) << scanUp.optP1_fom
                  << " cm  (N/FWHM = "   << scanUp.optVal_fom  << ")" << std::endl;
    }
    if (doDownstream && scanDown.valid) {
        std::cout << "║  DOWNSTREAM" << std::endl;
        std::cout << "║    Best resolution : P1 = " << std::setw(6) << scanDown.optP1_fwhm
                  << " cm  (mean FWHM = " << scanDown.optVal_fwhm << " GeV/c)" << std::endl;
        std::cout << "║    Best acceptance : P1 = " << std::setw(6) << scanDown.optP1_acc
                  << " cm  (peak N = "   << scanDown.optVal_acc  << ")" << std::endl;
        std::cout << "║    Best FOM        : P1 = " << std::setw(6) << scanDown.optP1_fom
                  << " cm  (N/FWHM = "   << scanDown.optVal_fom  << ")" << std::endl;
    }
    std::cout <<   "╚══════════════════════════════════════════════════════════╝" << std::endl;
}

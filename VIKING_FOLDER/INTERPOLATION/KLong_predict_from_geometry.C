// ============================================================================
// KLong_predict_from_geometry.C
//
// PURPOSE:
//   Given a parent directory containing multiple simulation archives (each
//   with *combined_vectors.root and *combined_acceptance.root), build
//   TSpline3 interpolants of the resolution and acceptance fit parameters
//   as a function of pizza position P1.  Then, for a user-supplied target
//   P1 (and T3 tracking mode), predict the resolution and acceptance curves
//   WITHOUT requiring any new simulation data.
//
// USAGE (via run_predict.sh, or directly):
//   root -l -b -q 'KLong_predict_from_geometry.C+("/path/to/ARCHIVED_RESULTS", 350.0)'
//   root -l -b -q 'KLong_predict_from_geometry.C+("/path/to/ARCHIVED_RESULTS", 350.0, "upstream")'
//   root -l -b -q 'KLong_predict_from_geometry.C+("/path/to/ARCHIVED_RESULTS", 350.0, "downstream", "/output/dir")'
//
// ARGUMENTS:
//   parentDir  — directory containing one or more archive sub-directories
//   targetP1   — pizza position to predict curves for [cm]
//   t3Mode     — "upstream" | "downstream" | "auto" (default "auto")
//                  upstream  : T3 ≈ T1+80  (pizzas plus close-pair downstream tracker)
//                  downstream: T3 = 680    (pizzas with far fixed tracker)
//                  auto      : choose the mode with more training points
//   outputDir  — where to write the PNG; defaults to parentDir/plots
//
// MODELS:
//   Resolution: FWHM(p) = a0 + exp(a1 + a2*p)  (pol0+expo)
//   Acceptance: A(p) = N / [(1+exp(-k*(p-p0))) * (1+exp(m*(p-p1)))]
//
// INTERPOLATION:
//   For each fit parameter, a TSpline3 is constructed over the training P1
//   values.  If targetP1 lies outside the training range, a linear
//   extrapolation from the two nearest endpoints is used and a warning is
//   printed.
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TColor.h"
#include "TROOT.h"
#include "TFitResult.h"
#include "TMultiGraph.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <map>

// ============================================================================
// GEOMETRY HELPERS
// ============================================================================

static double pred_extractPos(const std::string& label, const std::string& det) {
    std::string srch = det + "-";
    size_t pos = label.find(srch);
    if (pos == std::string::npos) return -1.0;
    pos += srch.size();
    size_t end = label.find_first_not_of("0123456789", pos);
    return std::stod(label.substr(pos, end - pos));
}

// Returns a mode string that identifies the T3 configuration:
//   "upstream"    : T3 > 0 and T3 < 600  (close downstream T3/T4 pair)
//   "downstream"  : T3 >= 600             (far fixed T3/T4 pair)
//   "no_T3"       : T3 == 0 or missing    (no downstream tracker pair)
static std::string pred_t3Mode(const std::string& label) {
    double t3 = pred_extractPos(label, "T3");
    if (t3 > 0. && t3 < 660.) return "upstream";
    if (t3 >= 660.)           return "downstream";
    return "no_T3";
}

// ============================================================================
// SPLINE INTERPOLANT  (with linear extrapolation outside range)
// ============================================================================

struct ParamSpline {
    std::vector<double> xs, ys;  // sorted training points
    TSpline3* sp = nullptr;

    void build(const std::string& name) {
        if (xs.size() < 2) { sp = nullptr; return; }
        // Sort by x
        std::vector<int> idx(xs.size()); for (int i=0;i<(int)idx.size();i++) idx[i]=i;
        std::sort(idx.begin(),idx.end(),[&](int a,int b){return xs[a]<xs[b];});
        std::vector<double> sx, sy;
        for (int i : idx) { sx.push_back(xs[i]); sy.push_back(ys[i]); }
        xs = sx; ys = sy;
        sp = new TSpline3(name.c_str(),&xs[0],&ys[0],(int)xs.size());
    }

    double eval(double x) const {
        if (!sp || xs.empty()) return 0.;
        int n = (int)xs.size();
        if (x < xs[0]) {
            // Linear extrapolation from first two points
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
};

// ============================================================================
// PER-CONFIG FITTING  (inlined from KLong_fit_resolution / _acceptance)
// ============================================================================

struct ResFitResult {
    double p1, t1;
    std::string mode;
    double a0, a1, a2, chi2ndf;  // FWHM(p) = a0 + exp(a1 + a2*p)  (pol0+expo)
    bool ok;
};

struct AccFitResult {
    double p1, t1;
    std::string mode;
    double N, k, p0, m, p1param, chi2ndf;
    bool ok;
};

ResFitResult pred_fitResolution(const std::string& vecFile,
                                 const std::string& label) {
    ResFitResult r;
    r.p1 = pred_extractPos(label,"P1");
    r.t1 = pred_extractPos(label,"T1");
    r.mode = pred_t3Mode(label);
    r.ok = false;
    r.a0 = r.a1 = r.a2 = r.chi2ndf = 0.;

    const int NB = 18;
    double pB[NB+1] = {0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,
                       3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,10.0};
    const double ACUT = 1.0;
    const int MINENT = 30;
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
    for (int i=0;i<NB;i++) {
        hb[i] = new TH1D(Form("pred_rh_%s_%d",label.c_str(),i),"",100,-3.,3.);
        hb[i]->SetDirectory(nullptr);
    }

    Long64_t nEnt = tree->GetEntries();
    for (Long64_t ev=0;ev<nEnt;ev++) {
        tree->GetEntry(ev);
        if (!reco_p||!true_p) continue;
        for (int k=0;k<(int)reco_p->size();k++) {
            double tp=(*true_p)[k], dp=(*reco_p)[k]-tp;
            if (tp<=0.) continue;
            if (std::abs(dp)/tp > ACUT) continue;
            for (int b=0;b<NB;b++) if(tp>=pB[b]&&tp<pB[b+1]){hb[b]->Fill(dp);break;}
        }
    }
    tf->Close(); delete tf;

    std::vector<double> gx, gy, gex, gey;
    for (int i=0;i<NB;i++) {
        if (hb[i]->GetEntries()<MINENT){delete hb[i];continue;}
        double mean=hb[i]->GetMean(), rms=hb[i]->GetRMS();
        TF1* gf = new TF1(Form("pred_rgf_%s_%d",label.c_str(),i),"gaus",mean-2.*rms,mean+2.*rms);
        gf->SetParameters(hb[i]->GetMaximum(),mean,rms);
        hb[i]->Fit(gf,"RQN");
        double sig=std::abs(gf->GetParameter(2)), sige=gf->GetParError(2);
        double mid=0.5*(pB[i]+pB[i+1]);
        if (mid>=FIT_PMIN && mid<=FIT_PMAX) {
            const double MIN_REL_ERR = 0.20;
            double fwhm   = 2.3548 * sig;
            double fwhm_e = 2.3548 * (sige > 0. ? sige : 1e-6);
            if (fwhm_e < MIN_REL_ERR * fwhm) fwhm_e = MIN_REL_ERR * fwhm;
            gx.push_back(mid); gy.push_back(fwhm);
            gex.push_back(0.5*(pB[i+1]-pB[i]));
            gey.push_back(fwhm_e);
        }
        hb[i]->GetListOfFunctions()->Remove(gf);
        delete gf; delete hb[i];
    }
    if ((int)gx.size()<3) return r;

    // Build TH1F from FWHM data and fit it with pol0+expo "SRL"
    // (matches supervisor's suggestion: hResData->Fit("fo11","RL"))
    int nb = (int)gx.size();
    std::vector<double> edges;
    for (int i=0;i<nb;i++) edges.push_back(gx[i]-gex[i]);
    edges.push_back(gx[nb-1]+gex[nb-1]);
    TH1F* hd = new TH1F(Form("pred_hd_%s",label.c_str()),
        "",nb,edges.data());
    hd->SetDirectory(nullptr);
    for (int i=0;i<nb;i++) { hd->SetBinContent(i+1,gy[i]); hd->SetBinError(i+1,gey[i]); }

    TF1* fm = new TF1(Form("pred_rfm_%s",label.c_str()),"pol0(0)+expo(1)",FIT_PMIN,FIT_PMAX);
    fm->SetParameter(0, 0.05);
    fm->SetParameter(1, -0.5);
    fm->SetParameter(2, -0.5);
    TFitResultPtr fr = hd->Fit(fm,"SRLN");
    hd->GetListOfFunctions()->Remove(fm);
    delete hd;

    if ((int)fr==0 || fr->IsValid()) {
        r.a0=fm->GetParameter(0); r.a1=fm->GetParameter(1); r.a2=fm->GetParameter(2);
        double ndf=fm->GetNDF();
        r.chi2ndf=(ndf>0)?fm->GetChisquare()/ndf:-1.;
        r.ok=true;
    }
    delete fm;
    return r;
}

AccFitResult pred_fitAcceptance(const std::string& accFile,
                                 const std::string& label) {
    AccFitResult r;
    r.p1 = pred_extractPos(label,"P1");
    r.t1 = pred_extractPos(label,"T1");
    r.mode = pred_t3Mode(label);
    r.ok = false;
    r.N=r.k=r.p0=r.m=r.p1param=r.chi2ndf=0.;

    const double P_MIN=0., P_MAX=10., BIN_W=0.25;
    const int NBINS=(int)((P_MAX-P_MIN)/BIN_W);

    TFile* tf = TFile::Open(accFile.c_str(),"READ");
    if (!tf||tf->IsZombie()){if(tf)delete tf;return r;}
    TTree* tree=(TTree*)tf->Get("kaonEventInfo");
    if(!tree){tf->Close();delete tf;return r;}

    std::vector<double>* true_mom_vec=nullptr;
    std::vector<int>*    reco_flag_vec=nullptr;
    int n_triple=0;
    tree->SetBranchAddress("true_mom_vec",&true_mom_vec);
    tree->SetBranchAddress("reco_flag_vec",&reco_flag_vec);
    tree->SetBranchAddress("n_triple_pion_events",&n_triple);

    std::vector<double> totalCnt(NBINS,0.), recoCnt(NBINS,0.);
    Long64_t nEnt=tree->GetEntries(); double totalTriple=0.;
    for (Long64_t ev=0;ev<nEnt;ev++) {
        tree->GetEntry(ev);
        totalTriple+=n_triple;
        if(!true_mom_vec||!reco_flag_vec) continue;
        for(int k=0;k<(int)true_mom_vec->size();k++){
            double p=(*true_mom_vec)[k];
            int bn=(int)((p-P_MIN)/BIN_W);
            if(bn<0||bn>=NBINS) continue;
            totalCnt[bn]++;
            if((*reco_flag_vec)[k]==1) recoCnt[bn]++;
        }
    }
    tf->Close(); delete tf;
    if(totalTriple<=0.) return r;

    std::vector<double> bx,by,bey;
    for(int i=0;i<NBINS;i++){
        if(totalCnt[i]<5.) continue;
        double eff=recoCnt[i]/totalTriple;
        double mid=P_MIN+(i+0.5)*BIN_W;
        bx.push_back(mid); by.push_back(eff);
        double err=std::sqrt(eff*(1.-std::min(eff,1.))/totalCnt[i]);
        bey.push_back(err>0.?err:1e-6);
    }
    if((int)bx.size()<6) return r;

    std::vector<double> bez(bx.size(),0.);
    TGraphErrors* gg=new TGraphErrors((int)bx.size(),&bx[0],&by[0],&bez[0],&bey[0]);
    TF1* fm=new TF1(Form("pred_afm_%s",label.c_str()),
        "[0]/((1.+exp(-[1]*(x-[2])))*(1.+exp([3]*(x-[4]))))",0.5,10.);
    fm->SetParNames("N","k","p0","m","p1");
    double ymax=*std::max_element(by.begin(),by.end());
    fm->SetParameters(ymax,4.,1.1,0.2,5.);
    fm->SetParLimits(0,0.,1.); fm->SetParLimits(1,0.5,20.);
    fm->SetParLimits(2,0.5,2.5); fm->SetParLimits(3,0.,2.);
    fm->SetParLimits(4,1.,9.);
    TFitResultPtr fr=gg->Fit(fm,"QSRN");
    delete gg;

    if((int)fr==0||fr->IsValid()){
        r.N=fm->GetParameter(0); r.k=fm->GetParameter(1);
        r.p0=fm->GetParameter(2); r.m=fm->GetParameter(3);
        r.p1param=fm->GetParameter(4);
        double ndf=fm->GetNDF();
        r.chi2ndf=(ndf>0)?fm->GetChisquare()/ndf:-1.;
        r.ok=true;
    }
    delete fm;
    return r;
}

// ============================================================================
// MAIN
// ============================================================================

void KLong_predict_from_geometry(const char* parentDirC,
                                   double targetP1,
                                   const char* t3ModeC   = "auto",
                                   const char* outputDirC = "") {

    std::string parentDir(parentDirC);
    while (!parentDir.empty() && parentDir.back()=='/') parentDir.pop_back();

    std::string t3Mode(t3ModeC);
    std::string outputDir(outputDirC);
    if (outputDir.empty()) outputDir = parentDir + "/plots";
    gSystem->mkdir(outputDir.c_str(), kTRUE);

    // -------------------------------------------------------------------------
    // 1. Discover all training files
    // -------------------------------------------------------------------------
    std::cout << "\n=== Scanning training archives under: " << parentDir << " ===" << std::endl;

    auto findFiles = [&](const std::string& suffix) -> std::vector<std::pair<std::string,std::string>> {
        // Returns (filepath, label) pairs
        std::string cmd = "find \""+parentDir+"\" -maxdepth 2 -name '*"+suffix+"' 2>/dev/null | sort";
        FILE* pipe = popen(cmd.c_str(),"r");
        std::vector<std::pair<std::string,std::string>> out;
        if (!pipe) return out;
        char buf[2048];
        while (fgets(buf,sizeof(buf),pipe)) {
            std::string fp(buf);
            fp.erase(std::remove(fp.begin(),fp.end(),'\n'),fp.end());
            if (fp.empty()) continue;
            size_t sl=fp.find_last_of('/');
            std::string base=(sl!=std::string::npos)?fp.substr(sl+1):fp;
            size_t sfx=base.find("_"+suffix);
            if (sfx==std::string::npos) sfx=base.find(suffix);
            std::string lbl=(sfx!=std::string::npos)?base.substr(0,sfx):base;
            // Remove leading archive dir name from label if present
            out.push_back({fp,lbl});
        }
        pclose(pipe);
        return out;
    };

    auto vecFiles = findFiles("combined_vectors.root");
    auto accFiles = findFiles("combined_acceptance.root");
    std::cout << "  Found " << vecFiles.size() << " vector files, "
              << accFiles.size() << " acceptance files." << std::endl;

    // -------------------------------------------------------------------------
    // 2. Run fits on all training data
    // -------------------------------------------------------------------------
    std::vector<ResFitResult> resResults;
    std::vector<AccFitResult> accResults;

    // Deduplicate by label (different archives may have the same config)
    std::map<std::string,bool> seenVec, seenAcc;

    std::cout << "\n--- Resolution fits ---" << std::endl;
    for (auto& [fp,lbl] : vecFiles) {
        if (seenVec[lbl]) continue;
        seenVec[lbl]=true;
        std::cout << "  " << lbl << " ... " << std::flush;
        ResFitResult r = pred_fitResolution(fp,lbl);
        if (r.ok) {
            resResults.push_back(r);
            std::cout << "a0=" << r.a0 << "  a1=" << r.a1 << "  a2=" << r.a2
                      << "  mode=" << r.mode << std::endl;
        } else {
            std::cout << "FAILED" << std::endl;
        }
    }

    std::cout << "\n--- Acceptance fits ---" << std::endl;
    for (auto& [fp,lbl] : accFiles) {
        if (seenAcc[lbl]) continue;
        seenAcc[lbl]=true;
        std::cout << "  " << lbl << " ... " << std::flush;
        AccFitResult r = pred_fitAcceptance(fp,lbl);
        if (r.ok) {
            accResults.push_back(r);
            std::cout << "N=" << r.N << "  k=" << r.k << "  mode=" << r.mode << std::endl;
        } else {
            std::cout << "FAILED" << std::endl;
        }
    }

    if (resResults.empty() && accResults.empty()) {
        std::cerr << "\nERROR: No fits succeeded; cannot build interpolants." << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // 3. Determine T3 mode
    // -------------------------------------------------------------------------
    // Count available training points per mode
    auto countMode = [](const std::vector<ResFitResult>& v, const std::string& m){
        int c=0; for(auto& r:v) if(r.mode==m) c++; return c;
    };
    int nUp   = countMode(resResults,"upstream");
    int nDown = countMode(resResults,"downstream");
    int nNo   = countMode(resResults,"no_T3");

    std::string chosenMode = t3Mode;
    if (chosenMode=="auto") {
        // Pick whichever has the most training points
        if (nUp >= nDown && nUp >= nNo)       chosenMode = "upstream";
        else if (nDown >= nUp && nDown >= nNo) chosenMode = "downstream";
        else                                   chosenMode = "no_T3";
        std::cout << "\nauto mode: upstream=" << nUp
                  << "  downstream=" << nDown
                  << "  no_T3=" << nNo
                  << " → choosing '" << chosenMode << "'" << std::endl;
    }
    std::cout << "Using T3 mode: " << chosenMode << std::endl;

    // -------------------------------------------------------------------------
    // 4. Build interpolation splines
    // -------------------------------------------------------------------------
    ParamSpline sp_a0, sp_a1, sp_a2, spN_p, sp_k, sp_p0, sp_m, sp_p1p;

    // Deduplicate by P1: if two configs have the same P1 in this mode,
    // keep the last one seen (map overwrites).  This prevents TSpline3
    // from receiving duplicate x values which cause NaN evaluation.
    {
        std::map<double,double> ma0, ma1, ma2;
        for (auto& r : resResults) {
            if (r.mode != chosenMode) continue;
            ma0[r.p1] = r.a0;  ma1[r.p1] = r.a1;  ma2[r.p1] = r.a2;
        }
        for (auto& [p,v] : ma0) { sp_a0.xs.push_back(p); sp_a0.ys.push_back(v); }
        for (auto& [p,v] : ma1) { sp_a1.xs.push_back(p); sp_a1.ys.push_back(v); }
        for (auto& [p,v] : ma2) { sp_a2.xs.push_back(p); sp_a2.ys.push_back(v); }
    }
    {
        std::map<double,double> mN,mk,mp0,mm,mp1;
        for (auto& r : accResults) {
            if (r.mode != chosenMode) continue;
            mN[r.p1]=r.N; mk[r.p1]=r.k; mp0[r.p1]=r.p0;
            mm[r.p1]=r.m; mp1[r.p1]=r.p1param;
        }
        for (auto& [p,v]:mN)  { spN_p.xs.push_back(p); spN_p.ys.push_back(v); }
        for (auto& [p,v]:mk)  { sp_k.xs.push_back(p);  sp_k.ys.push_back(v); }
        for (auto& [p,v]:mp0) { sp_p0.xs.push_back(p); sp_p0.ys.push_back(v); }
        for (auto& [p,v]:mm)  { sp_m.xs.push_back(p);  sp_m.ys.push_back(v); }
        for (auto& [p,v]:mp1) { sp_p1p.xs.push_back(p);sp_p1p.ys.push_back(v); }
    }

    sp_a0.build("sp_a0");  sp_a1.build("sp_a1");  sp_a2.build("sp_a2");
    spN_p.build("spN");  sp_k.build("sp_k");
    sp_p0.build("sp_p0"); sp_m.build("sp_m"); sp_p1p.build("sp_p1");

    if (sp_a0.xs.empty() && spN_p.xs.empty()) {
        std::cerr << "\nERROR: No training data for mode '" << chosenMode
                  << "'. Try a different t3Mode." << std::endl;
        return;
    }

    // -------------------------------------------------------------------------
    // 5. Predict at targetP1
    // -------------------------------------------------------------------------
    bool extrapRes = sp_a0.isExtrapolating(targetP1);
    bool extrapAcc = spN_p.isExtrapolating(targetP1);
    if (extrapRes || extrapAcc) {
        std::cout << "\n  WARNING: target P1=" << targetP1
                  << " is outside the training range ["
                  << std::min(sp_a0.xMin(),spN_p.xMin()) << ", "
                  << std::max(sp_a0.xMax(),spN_p.xMax())
                  << "].  Using linear extrapolation — results may be unreliable." << std::endl;
    }

    double pred_a0 = sp_a0.eval(targetP1);
    double pred_a1 = sp_a1.eval(targetP1);
    double pred_a2 = sp_a2.eval(targetP1);
    // No clamping needed for pol0+expo parameters
    double predAcc_N  = spN_p.eval(targetP1);
    double pred_k     = sp_k.eval(targetP1);
    double pred_p0    = sp_p0.eval(targetP1);
    double pred_m     = sp_m.eval(targetP1);
    double pred_p1val = sp_p1p.eval(targetP1);

    predAcc_N= std::max(std::min(predAcc_N,1.), 1e-4);
    pred_k   = std::max(pred_k,  0.5);
    pred_p0  = std::max(pred_p0, 0.1);
    pred_m   = std::max(pred_m,  0.);
    pred_p1val = std::max(pred_p1val, 0.5);

    std::cout << "\n=== Predicted parameters for P1=" << targetP1
              << " cm  (mode=" << chosenMode << ") ===" << std::endl;
    std::cout << "  Resolution:  a0=" << pred_a0 << "  a1=" << pred_a1
              << "  a2=" << pred_a2 << std::endl;
    std::cout << "  Acceptance:  N=" << predAcc_N << "  k=" << pred_k
              << "  p0=" << pred_p0 << "  m=" << pred_m
              << "  p1=" << pred_p1val << std::endl;

    // -------------------------------------------------------------------------
    // 6. Draw predicted curves
    // -------------------------------------------------------------------------
    // Canvas: two plot pads (top 70%) + annotation/legend pad (bottom 30%)
    TCanvas* c = new TCanvas("cPred","KLong Predicted Curves",1600,960);

    TPad* pRes = new TPad("pRes","", 0.00, 0.30, 0.50, 0.94);
    TPad* pAcc = new TPad("pAcc","", 0.50, 0.30, 1.00, 0.94);
    TPad* pBot = new TPad("pBot","", 0.00, 0.00, 1.00, 0.30);

    for (TPad* p : {pRes, pAcc}) {
        p->SetGrid(1,1);
        p->SetLeftMargin(0.15); p->SetRightMargin(0.06);
        p->SetTopMargin(0.10);  p->SetBottomMargin(0.15);
        // Remove the pad frame border box
        p->SetFrameLineColor(0); p->SetFrameLineWidth(0);
    }
    pBot->SetLeftMargin(0.01); pBot->SetRightMargin(0.01);
    pBot->SetTopMargin(0.04);  pBot->SetBottomMargin(0.04);

    pRes->Draw(); pAcc->Draw(); pBot->Draw();

    // ---- Resolution plot ----
    pRes->cd();

    TF1* fRes = new TF1("fPredRes","pol0(0)+expo(1)",0.0,5.0);
    fRes->SetParameters(pred_a0, pred_a1, pred_a2);
    fRes->SetLineColor(kAzure+2); fRes->SetLineWidth(3);

    TH1D* hResFrame = new TH1D("hPredResFrame",
        ";True Momentum [GeV/c];FWHM(#Deltap) [GeV/c]",
        100, 0.0, 5.5);
    hResFrame->SetMinimum(0.);
    hResFrame->SetMaximum(std::min(fRes->GetMaximum(0.0,5.0)*1.3, 1.5));
    hResFrame->SetLineColor(0);
    hResFrame->SetStats(0);
    hResFrame->GetXaxis()->SetTitleSize(0.048); hResFrame->GetXaxis()->SetLabelSize(0.040);
    hResFrame->GetYaxis()->SetTitleSize(0.048); hResFrame->GetYaxis()->SetLabelSize(0.040);
    hResFrame->GetXaxis()->SetNdivisions(3);   // labels every ~2 GeV/c
    hResFrame->GetYaxis()->SetNdivisions(5);
    hResFrame->Draw("AXIS");

    for (auto& r : resResults) {
        if (r.mode!=chosenMode) continue;
        TF1* ft = new TF1(Form("fTrRes_%s_%.0f",chosenMode.c_str(),r.p1),
            "pol0(0)+expo(1)",0.0,5.0);
        ft->SetParameters(r.a0, r.a1, r.a2);
        ft->SetLineColor(kGray+1); ft->SetLineWidth(1); ft->SetLineStyle(2);
        ft->Draw("same");
    }
    fRes->Draw("same");

    // Panel title
    TLatex pltx; pltx.SetNDC(); pltx.SetTextAlign(22);
    pltx.SetTextSize(0.055); pltx.SetTextFont(62);
    pltx.DrawLatex(0.565, 0.955, "Momentum Resolution");

    // ---- Acceptance plot ----
    pAcc->cd();

    TF1* fAcc = new TF1("fPredAcc",
        "[0]/((1.+exp(-[1]*(x-[2])))*(1.+exp([3]*(x-[4]))))",
        0.25, 9.);
    fAcc->SetParameters(predAcc_N, pred_k, pred_p0, pred_m, pred_p1val);
    fAcc->SetLineColor(kOrange+2); fAcc->SetLineWidth(3);

    TH1D* hAccFrame = new TH1D("hPredAccFrame",
        ";True Momentum [GeV/c];Acceptance",
        100, 0.25, 9.);
    hAccFrame->SetMinimum(0.);
    hAccFrame->SetMaximum(predAcc_N * 1.3);
    hAccFrame->SetLineColor(0);
    hAccFrame->SetStats(0);
    hAccFrame->GetXaxis()->SetTitleSize(0.048); hAccFrame->GetXaxis()->SetLabelSize(0.040);
    hAccFrame->GetYaxis()->SetTitleSize(0.048); hAccFrame->GetYaxis()->SetLabelSize(0.040);
    hAccFrame->GetXaxis()->SetNdivisions(4);   // labels every ~2 GeV/c
    hAccFrame->GetYaxis()->SetNdivisions(5);
    hAccFrame->Draw("AXIS");

    for (auto& r : accResults) {
        if (r.mode!=chosenMode) continue;
        TF1* ft = new TF1(Form("fTrAcc_%s_%.0f",chosenMode.c_str(),r.p1),
            "[0]/((1.+exp(-[1]*(x-[2])))*(1.+exp([3]*(x-[4]))))",0.25,9.);
        ft->SetParameters(r.N,r.k,r.p0,r.m,r.p1param);
        ft->SetLineColor(kGray+1); ft->SetLineWidth(1); ft->SetLineStyle(2);
        ft->Draw("same");
    }
    fAcc->Draw("same");

    pltx.DrawLatex(0.565, 0.955, "Acceptance");

    // ---- Bottom annotation pad ----
    pBot->cd();

    TLatex ltx; ltx.SetNDC();

    // Resolution block (left half of bottom pad, x 0.02 – 0.49)
    ltx.SetTextFont(62); ltx.SetTextSize(0.13);
    ltx.DrawLatex(0.025, 0.82, "Resolution model:");
    ltx.SetTextFont(42); ltx.SetTextSize(0.115);
    ltx.DrawLatex(0.025, 0.62, "FWHM(p) = a_{0} + exp(a_{1} + a_{2} #upoint p)");
    ltx.SetTextSize(0.095);
    ltx.DrawLatex(0.025, 0.44, Form("a_{0} = %.5f", pred_a0));
    ltx.DrawLatex(0.025, 0.30, Form("a_{1} = %.5f", pred_a1));
    ltx.DrawLatex(0.025, 0.16, Form("a_{2} = %.5f", pred_a2));
    if (extrapRes) {
        ltx.SetTextColor(kRed+1); ltx.SetTextSize(0.09);
        ltx.DrawLatex(0.025, 0.04, "[extrapolated]");
        ltx.SetTextColor(kBlack);
    }

    // Line key for resolution (solid = predicted, dashed = training)
    // --- Tune these four values to move/resize the whole legend block ---
    double legX      = 0.30;   // left edge of the legend lines (NDC x)
    double legY      = 0.71;   // y of the top (predicted) row  (NDC y)
    double legLineW  = 0.055;  // length of each legend line segment
    double legRowSep = 0.25;   // vertical gap between rows (positive = downward)
    // ---------------------------------------------------------------------
    double legTextX  = legX + legLineW + 0.008;   // x where label text starts
    double legTextOY = -0.03;                      // label y offset below the line

    TLine lnRes; lnRes.SetNDC();
    lnRes.SetLineColor(kAzure+2); lnRes.SetLineWidth(3); lnRes.SetLineStyle(1);
    lnRes.DrawLineNDC(legX, legY, legX + legLineW, legY);
    ltx.SetTextFont(42); ltx.SetTextSize(0.075);
    ltx.DrawLatex(legTextX, legY + legTextOY, Form("Predicted (P1=%.0f cm)", targetP1));

    double legY2 = legY - legRowSep;
    TLine lnTrRes; lnTrRes.SetNDC();
    lnTrRes.SetLineColor(kGray+1); lnTrRes.SetLineWidth(1); lnTrRes.SetLineStyle(2);
    lnTrRes.DrawLineNDC(legX, legY2, legX + legLineW, legY2);
    ltx.DrawLatex(legTextX, legY2 + legTextOY, "Training configs");

    // Vertical separator
    TLine sep; sep.SetNDC();
    sep.SetLineColor(kGray+2); sep.SetLineWidth(1); sep.SetLineStyle(1);
    sep.DrawLineNDC(0.50, 0.05, 0.50, 0.95);

    // Acceptance block (right half, x 0.52 – 0.98)
    ltx.SetTextFont(62); ltx.SetTextSize(0.13);
    ltx.DrawLatex(0.525, 0.82, "Acceptance model:");
    ltx.SetTextFont(42); ltx.SetTextSize(0.100);
    ltx.DrawLatex(0.525, 0.62, "A(p) = N / [(1+e^{-k(p-p_{0})})(1+e^{m(p-p_{1})})]");
    ltx.SetTextSize(0.090);
    ltx.DrawLatex(0.525, 0.44, Form("N = %.5f    k = %.3f    p_{0} = %.3f",
                                     predAcc_N, pred_k, pred_p0));
    ltx.DrawLatex(0.525, 0.27, Form("m = %.3f    p_{1} = %.3f", pred_m, pred_p1val));
    if (extrapAcc) {
        ltx.SetTextColor(kRed+1); ltx.SetTextSize(0.09);
        ltx.DrawLatex(0.525, 0.11, "[extrapolated]");
        ltx.SetTextColor(kBlack);
    }

    // Line key for acceptance (mirrors resolution legend layout on right half)
    double accLegX     = 0.80;   // left edge of acceptance legend lines (NDC x)
    double accLegLineW = legLineW;
    double accLegTextX = accLegX + accLegLineW + 0.008;

    TLine lnAcc; lnAcc.SetNDC();
    lnAcc.SetLineColor(kOrange+2); lnAcc.SetLineWidth(3); lnAcc.SetLineStyle(1);
    lnAcc.DrawLineNDC(accLegX, legY, accLegX + accLegLineW, legY);
    ltx.SetTextSize(0.075);
    ltx.DrawLatex(accLegTextX, legY + legTextOY, Form("Predicted (P1=%.0f cm)", targetP1));

    TLine lnTrAcc; lnTrAcc.SetNDC();
    lnTrAcc.SetLineColor(kGray+1); lnTrAcc.SetLineWidth(1); lnTrAcc.SetLineStyle(2);
    lnTrAcc.DrawLineNDC(accLegX, legY2, accLegX + accLegLineW, legY2);
    ltx.DrawLatex(accLegTextX, legY2 + legTextOY, "Training configs");

    // ---- Canvas-level title ----
    c->cd(0);
    TLatex tTitle; tTitle.SetNDC(); tTitle.SetTextAlign(22);
    tTitle.SetTextFont(62); tTitle.SetTextSize(0.036);
    tTitle.DrawLatex(0.5, 0.978,
        Form("Predicted detector response: P1 = %.0f cm   |   T3 mode: %s   |   Training: %d res + %d acc configs",
             targetP1, chosenMode.c_str(), (int)sp_a0.xs.size(), (int)spN_p.xs.size()));

    // -------------------------------------------------------------------------
    // 7. Save
    // -------------------------------------------------------------------------
    std::string outFile = Form("%s/KLong_predicted_P1%.0f_%s.png",
                               outputDir.c_str(), targetP1, chosenMode.c_str());
    c->SaveAs(outFile.c_str());
    std::string outRoot = outFile.substr(0, outFile.size()-4) + ".root";
    c->SaveAs(outRoot.c_str());
    std::cout << "Saved: " << outFile << std::endl;
    std::cout << "Saved: " << outRoot << std::endl;
    delete c;
}

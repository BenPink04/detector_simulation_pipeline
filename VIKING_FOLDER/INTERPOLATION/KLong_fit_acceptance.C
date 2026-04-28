// ============================================================================
// KLong_fit_acceptance.C
//
// PURPOSE:
//   For every detector configuration found in a given archive directory,
//   compute the reconstruction acceptance as a function of true kaon momentum,
//   then fit a logistic (sigmoid) model:
//
//       A(p) = N / ( 1 + exp(-k * (p - p0)) )
//
//   The pizza position P1 appears as a numerical constant in the per-config
//   legend entry.
//
//   Output:
//     <archive_dir>/plots/KLong_fit_acceptance.png
//       Dashed = binned efficiency data, Solid = fitted sigmoid, same colour
//       per config.  General equation shown as TLatex.
//
// USAGE:
//   root -l -b -q 'KLong_fit_acceptance.C+("/path/to/archive")'
//   (called automatically by run_pizza_fit.sh)
//
// INPUT TREE:  kaonEventInfo
//   branches: true_mom_vec (vector<double>), reco_flag_vec (vector<int>),
//             n_triple_pion_events (int)
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TColor.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"
#include "TROOT.h"
#include "TFitResult.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>

// ============================================================================
// HELPERS  (all prefixed acc_ to avoid name clashes with resolution macro)
// ============================================================================

double acc_extractPos(const std::string& label, const std::string& det) {
    std::string srch = det + "-";
    size_t pos = label.find(srch);
    if (pos == std::string::npos) return -1.0;
    pos += srch.size();
    size_t end = label.find_first_not_of("0123456789", pos);
    return std::stod(label.substr(pos, end - pos));
}

struct AccDetLayout {
    std::vector<double> trackers, pizzas, fris;
    double end_detector;
    AccDetLayout() : end_detector(0.) {}
};

AccDetLayout acc_parseConfig(const std::string& config) {
    AccDetLayout layout;
    std::istringstream ss(config);
    std::string token;
    while (std::getline(ss, token, '_')) {
        if (token.empty()) continue;
        size_t pos = token.find('-');
        if (pos == std::string::npos) continue;
        std::string dt = token.substr(0, pos);
        double val = std::stod(token.substr(pos + 1));
        if (val == 0.) continue;
        if      (dt[0]=='T') layout.trackers.push_back(val);
        else if (dt[0]=='P') layout.pizzas.push_back(val);
        else if (dt[0]=='F') layout.fris.push_back(val);
        else if (dt[0]=='E') layout.end_detector = val;
    }
    return layout;
}

int acc_getColor(double t1, double t1Lo, double t1Hi) {
    // Two-colour gradient: teal (0,180,180) → orange (255,140,0)
    // Normalised dynamically over the T1 range present in this archive.
    const float r0=0,   g0=180, b0=180;  // teal
    const float r1=255, g1=140, b1=0;    // orange
    double range = (t1Hi > t1Lo) ? (t1Hi - t1Lo) : 1.;
    float f = (float)std::max(0., std::min(1., (t1-t1Lo)/range));
    return TColor::GetColor((r0+f*(r1-r0))/255.f,
                            (g0+f*(g1-g0))/255.f,
                            (b0+f*(b1-b0))/255.f);
}

void acc_drawLayout(const std::vector<std::string>& labels,
                    const std::vector<int>& colours) {
    double mn = 1e9, mx = -1e9;
    for (int i = 0; i < (int)labels.size(); ++i) {
        AccDetLayout l = acc_parseConfig(labels[i]);
        for (int j=0;j<(int)l.trackers.size();++j){mn=std::min(mn,l.trackers[j]);mx=std::max(mx,l.trackers[j]);}
        for (int j=0;j<(int)l.pizzas.size();++j)  {mn=std::min(mn,l.pizzas[j]);  mx=std::max(mx,l.pizzas[j]);}
        for (int j=0;j<(int)l.fris.size();++j)    {mn=std::min(mn,l.fris[j]);    mx=std::max(mx,l.fris[j]);}
        if (l.end_detector > 0.) { mn=std::min(mn,l.end_detector); mx=std::max(mx,l.end_detector); }
    }
    double rng = mx - mn;
    if (rng <= 0.) rng = 1.;
    mn -= rng*0.1; mx += rng*0.1;
    gPad->Range(0,0,1,1);

    TText* ttl = new TText(0.5,0.97,"Detector Layouts");
    ttl->SetTextAlign(23); ttl->SetTextSize(0.045); ttl->Draw();

    int n = (int)labels.size();
    double ysp = 0.80 / (n + 1);
    for (int i = 0; i < n; ++i) {
        AccDetLayout l = acc_parseConfig(labels[i]);
        double yc = 0.90 - (i+1)*ysp;
        TLine* ind = new TLine(0.02,yc,0.09,yc);
        ind->SetLineColor(colours[i]); ind->SetLineWidth(3); ind->Draw();
        double xs = 0.78, xo = 0.14;
        for (int j=0;j<(int)l.trackers.size();++j){
            double x = xo+(l.trackers[j]-mn)/(mx-mn)*xs;
            TBox* b=new TBox(x-0.007,yc-0.013,x+0.007,yc+0.013);
            b->SetFillColor(kBlue); b->SetLineColor(kBlue); b->Draw();
        }
        for (int j=0;j<(int)l.pizzas.size();++j){
            double x = xo+(l.pizzas[j]-mn)/(mx-mn)*xs;
            TBox* b=new TBox(x-0.007,yc-0.013,x+0.007,yc+0.013);
            b->SetFillColor(kRed); b->SetLineColor(kRed); b->Draw();
        }
        for (int j=0;j<(int)l.fris.size();++j){
            double x = xo+(l.fris[j]-mn)/(mx-mn)*xs;
            TBox* b=new TBox(x-0.007,yc-0.013,x+0.007,yc+0.013);
            b->SetFillColor(kMagenta); b->SetLineColor(kMagenta); b->Draw();
        }
        if (l.end_detector > 0.) {
            double x = xo+(l.end_detector-mn)/(mx-mn)*xs;
            TBox* b=new TBox(x-0.007,yc-0.013,x+0.007,yc+0.013);
            b->SetFillColor(kGreen+1); b->SetLineColor(kGreen+1); b->Draw();
        }
    }
    double ly=0.07, lx=0.12, lsp=0.21;
    TBox* b1=new TBox(lx,      ly,lx+0.025,      ly+0.03); b1->SetFillColor(kBlue);    b1->Draw();
    TText*t1=new TText(lx+0.03,      ly+0.015,"Tracker"); t1->SetTextSize(0.032); t1->Draw();
    TBox* b2=new TBox(lx+lsp,  ly,lx+lsp+0.025,  ly+0.03); b2->SetFillColor(kRed);     b2->Draw();
    TText*t2=new TText(lx+lsp+0.03,  ly+0.015,"Pizza");   t2->SetTextSize(0.032); t2->Draw();
    TBox* b3=new TBox(lx+2*lsp,ly,lx+2*lsp+0.025,ly+0.03); b3->SetFillColor(kMagenta); b3->Draw();
    TText*t3=new TText(lx+2*lsp+0.03,ly+0.015,"FRI");     t3->SetTextSize(0.032); t3->Draw();
    TBox* b4=new TBox(lx+3*lsp,ly,lx+3*lsp+0.025,ly+0.03); b4->SetFillColor(kGreen+1); b4->Draw();
    TText*t4=new TText(lx+3*lsp+0.03,ly+0.015,"ToF");     t4->SetTextSize(0.032); t4->Draw();
}

// ============================================================================
// MAIN
// ============================================================================

void KLong_fit_acceptance(const char* archiveDirC = "") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    std::string archiveDir(archiveDirC);
    if (archiveDir.empty()) {
        std::cerr << "Usage: root -l -b -q 'KLong_fit_acceptance.C+(\"/path/to/archive\")'" << std::endl;
        return;
    }
    while (!archiveDir.empty() && archiveDir.back()=='/') archiveDir.pop_back();
    std::string plotDir = archiveDir + "/plots";
    gSystem->mkdir(plotDir.c_str(), kTRUE);

    // Collect acceptance files
    std::string cmd = "find "+archiveDir+" -maxdepth 1 -name '*combined_acceptance.root' 2>/dev/null | sort";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) { std::cerr << "popen failed" << std::endl; return; }
    std::vector<std::string> files, labels;
    char buf[1024];
    while (fgets(buf, sizeof(buf), pipe)) {
        std::string fp(buf);
        fp.erase(std::remove(fp.begin(),fp.end(),'\n'),fp.end());
        if (fp.empty()) continue;
        files.push_back(fp);
        size_t sl = fp.find_last_of('/');
        std::string base = (sl!=std::string::npos) ? fp.substr(sl+1) : fp;
        size_t sfx = base.find("_combined_acceptance.root");
        labels.push_back(sfx!=std::string::npos ? base.substr(0,sfx) : base);
    }
    pclose(pipe);
    if (files.empty()) { std::cerr << "No *combined_acceptance.root in " << archiveDir << std::endl; return; }
    std::cout << "Found " << (int)files.size() << " config(s)." << std::endl;

    // Momentum bins matching KLong_plot_compare_acceptance.C
    const double P_MIN = 0.0, P_MAX = 10.0, BIN_W = 0.25;
    const int NB = (int)((P_MAX - P_MIN) / BIN_W);
    const int MINENT = 5;

    // Canvas
    // Three-pad canvas: main plot (top-left), layout (top-right), legend (bottom)
    TCanvas* c = new TCanvas("cAccfit","KLong Acceptance Fit",1600,960);
    TPad* pPlot   = new TPad("pPlot",  "", 0.00, 0.30, 0.72, 1.00);
    TPad* pLayout = new TPad("pLayout","", 0.72, 0.30, 1.00, 1.00);
    TPad* pLeg    = new TPad("pLeg",   "", 0.00, 0.00, 1.00, 0.30);
    pPlot->SetLeftMargin(0.13); pPlot->SetBottomMargin(0.13); pPlot->SetGrid(1,1);
    pLeg->SetLeftMargin(0.01);  pLeg->SetRightMargin(0.01);
    pLeg->SetTopMargin(0.02);   pLeg->SetBottomMargin(0.02);
    pPlot->Draw(); pLayout->Draw(); pLeg->Draw();
    pPlot->cd();

    // Frame
    TH1D* hFr = new TH1D("hAccFr",
        "Reconstruction Acceptance: Data (dashed) vs Model (solid)"
        ";True Kaon Momentum [GeV/c];Acceptance",
        (int)(P_MAX/BIN_W), P_MIN, P_MAX);
    hFr->Draw("AXIS");  // y-range set after loop once data range is known

    // Accumulated during loop; all drawn in bottom panel after
    std::vector<std::string> legEntries;
    std::vector<int>         legColors;

    std::vector<std::string> layoutLabels;
    std::vector<int>         layoutColours;
    double globalAccYmin = 1.0, globalAccYmax = 0.0;
    std::vector<TH1F*> hToSave;

    // Pre-compute T1 range so the gradient spans the actual configs in this archive
    double t1Lo = 1e9, t1Hi = -1e9;
    for (auto& lbl : labels) {
        double t1 = acc_extractPos(lbl,"T1");
        if (t1 > 0.) { t1Lo=std::min(t1Lo,t1); t1Hi=std::max(t1Hi,t1); }
    }
    if (t1Lo > t1Hi) { t1Lo = 240.; t1Hi = 480.; }  // fallback

    for (int idx = 0; idx < (int)files.size(); ++idx) {
        const std::string& fp  = files[idx];
        const std::string& lbl = labels[idx];
        double p1  = acc_extractPos(lbl,"P1");
        double t1  = acc_extractPos(lbl,"T1");
        int    col = acc_getColor(t1, t1Lo, t1Hi);

        std::cout << "\n[" << idx+1 << "/" << (int)files.size() << "] " << lbl << "  P1=" << p1 << std::endl;

        TFile* tf = TFile::Open(fp.c_str(),"READ");
        if (!tf || tf->IsZombie()) {
            std::cerr << "  WARN: cannot open " << fp << std::endl;
            if (tf) delete tf;
            continue;
        }
        TTree* tree = (TTree*)tf->Get("kaonEventInfo");
        if (!tree) {
            std::cerr << "  WARN: no kaonEventInfo tree" << std::endl;
            tf->Close(); delete tf; continue;
        }

        std::vector<double>* true_mom  = nullptr;
        std::vector<int>*    reco_flag = nullptr;
        int n_triple = 0;
        tree->SetBranchAddress("true_mom_vec",  &true_mom);
        tree->SetBranchAddress("reco_flag_vec", &reco_flag);
        tree->SetBranchAddress("n_triple_pion_events", &n_triple);

        // Count all generated and reconstructed per bin using integer arrays
        std::vector<long long> nAll(NB, 0LL), nReco(NB, 0LL);

        Long64_t nEnt = tree->GetEntries();
        std::cout << "  Entries: " << nEnt << std::endl;
        for (Long64_t ev = 0; ev < nEnt; ++ev) {
            tree->GetEntry(ev);
            if (!true_mom || !reco_flag) continue;
            for (int k = 0; k < n_triple && k < (int)true_mom->size(); ++k) {
                double tp = (*true_mom)[k];
                if (tp < P_MIN || tp >= P_MAX) continue;
                int b = (int)((tp - P_MIN) / BIN_W);
                if (b < 0 || b >= NB) continue;
                nAll[b]++;
                if (k < (int)reco_flag->size() && (*reco_flag)[k] > 0) nReco[b]++;
            }
        }
        tf->Close(); delete tf;

        // Build efficiency arrays
        std::vector<double> gx, gy, gex, gey;
        for (int i = 0; i < NB; ++i) {
            if (nAll[i] < MINENT) continue;
            double eff  = (double)nReco[i] / nAll[i];
            double effe = std::sqrt(eff*(1.-eff)/nAll[i]);
            if (effe <= 0.) effe = 1e-4;
            gx.push_back( P_MIN + (i + 0.5)*BIN_W );
            gy.push_back( eff );
            gex.push_back( 0.5*BIN_W );
            gey.push_back( effe );
        }

        if ((int)gx.size() < 3) {
            std::cerr << "  Too few bins for fit, skipping." << std::endl;
            continue;
        }

        // Track y range
        for (int i = 0; i < (int)gy.size(); ++i) {
            if (gy[i] < globalAccYmin) globalAccYmin = gy[i];
            if (gy[i] > globalAccYmax) globalAccYmax = gy[i];
        }

        // DATA: dashed TGraphErrors
        TGraphErrors* gd = new TGraphErrors((int)gx.size(),
            &gx[0], &gy[0], &gex[0], &gey[0]);
        gd->SetLineColor(col); gd->SetLineWidth(2); gd->SetLineStyle(2);
        gd->SetMarkerColor(col); gd->SetMarkerStyle(20); gd->SetMarkerSize(0.9);
        gd->Draw("LP same");

        // Fit MODEL: double-logistic — sharp rise * gentle fall
        //   A(p) = N / [(1+exp(-k*(p-p0))) * (1+exp(m*(p-p1)))]
        // k >> m gives sharp turn-on, gentle turn-off with flat peak between
        TGraphErrors* gfit = new TGraphErrors((int)gx.size(),
            &gx[0], &gy[0], &gex[0], &gey[0]);
        TF1* fm = new TF1(Form("afm_%d",idx),
            "[0]/((1.+exp(-[1]*(x-[2])))*(1.+exp([3]*(x-[4]))))", P_MIN, P_MAX);
        fm->SetParNames("N","k","p0","m","p1");
        fm->SetParameters(0.015, 10.0, 0.5, 0.8, 4.0);
        fm->SetParLimits(0, 1e-5, 1.0);
        fm->SetParLimits(1, 0.5,  100.0);   // sharp rise
        fm->SetParLimits(2, 0.0,  3.0);
        fm->SetParLimits(3, 0.01, 10.0);    // gentle fall
        fm->SetParLimits(4, 1.0,  P_MAX);
        TFitResultPtr rfit = gfit->Fit(fm,"QSRN");
        delete gfit;

        if ((int)rfit == 0 || rfit->IsValid()) {
            double N   = fm->GetParameter(0);
            double k   = fm->GetParameter(1);
            double p0  = fm->GetParameter(2);
            double m   = fm->GetParameter(3);
            double p1v = fm->GetParameter(4);
            double ndf = fm->GetNDF();
            double c2n = (ndf > 0) ? fm->GetChisquare()/ndf : -1.;
            std::cout << "  Fit OK: N=" << N << "  k=" << k << "  p0=" << p0
                      << "  m=" << m << "  p1=" << p1v << "  chi2/ndf=" << c2n << std::endl;
            fm->SetLineColor(col); fm->SetLineWidth(2); fm->SetLineStyle(1);
            fm->Draw("same");
            // Persist data and fit curve as TH1Fs for ROOT output
            {
                int nb = (int)gx.size();
                std::vector<double> edges;
                for (int i=0;i<nb;++i) edges.push_back(gx[i]-gex[i]);
                edges.push_back(gx[nb-1]+gex[nb-1]);
                TH1F* hd = new TH1F(Form("hAccData_P1%.0f",p1),
                    Form("Acceptance data P1=%.0f cm;True Momentum [GeV/c];Acceptance",p1),
                    nb, edges.data());
                hd->SetDirectory(nullptr);
                for (int i=0;i<nb;++i) { hd->SetBinContent(i+1,gy[i]); hd->SetBinError(i+1,gey[i]); }
                hToSave.push_back(hd);
                double xMax = gx.back()+gex.back();
                TH1F* hf = new TH1F(Form("hAccFit_P1%.0f",p1),
                    Form("Acceptance fit P1=%.0f cm;True Momentum [GeV/c];Acceptance",p1),
                    100, gx.front()-gex.front(), xMax);
                hf->SetDirectory(nullptr);
                for (int i=1;i<=100;++i) hf->SetBinContent(i, fm->Eval(hf->GetBinCenter(i)));
                hToSave.push_back(hf);
            }
            // Two-row legend entry: row 1 gets the coloured line marker
            legEntries.push_back(Form("P1=%.0f   N=%.4f   k=%.2f   p_{0}=%.2f", p1, N, k, p0));
            legColors.push_back(col);
            // Row 2: text only (sentinel colour -1 = no marker)
            legEntries.push_back(Form("      m=%.2f   p_{1}=%.2f   #chi^{2}/ndf=%.2f", m, p1v, c2n));
            legColors.push_back(-1);
        } else {
            std::cerr << "  Fit failed." << std::endl;
            legEntries.push_back(Form("P1=%.0f   [fit failed]", p1));
            legColors.push_back(col);
            legEntries.push_back(std::string(""));
            legColors.push_back(-1);
            delete fm;
        }
        layoutLabels.push_back(lbl);
        layoutColours.push_back(col);
    }

    // Adjust y-axis to data range
    if (globalAccYmax > globalAccYmin) {
        double apad = (globalAccYmax - globalAccYmin) * 0.15;
        hFr->GetYaxis()->SetRangeUser(
            std::max(0., globalAccYmin - apad),
            globalAccYmax + apad * 1.5);
    }
    pPlot->RedrawAxis();
    pPlot->Update();

    // Detector layout pad
    pLayout->cd();
    if (!layoutLabels.empty()) acc_drawLayout(layoutLabels, layoutColours);

    // ── Bottom legend panel ─────────────────────────────────────────────────
    pLeg->cd();
    pLeg->Range(0, 0, 1, 1);

    // Equation + key at the top of the pad
    TLatex eqLat;
    eqLat.SetNDC(); eqLat.SetTextSize(0.07); eqLat.SetTextAlign(12);
    eqLat.DrawLatex(0.02, 0.84,
        "Model:  A(p) = N / [(1 + e^{-k(p-p_{0})}) #upoint (1 + e^{m(p-p_{1})})]");

    // Dashed / solid key
    TLine aldash(0.50, 0.86, 0.56, 0.86);
    aldash.SetLineStyle(2); aldash.SetLineWidth(3); aldash.SetLineColor(kGray+2);
    aldash.Draw();
    TLatex altd; altd.SetNDC(); altd.SetTextSize(0.065); altd.SetTextAlign(12);
    altd.DrawLatex(0.565, 0.84, "Binned data");
    TLine alsol(0.70, 0.86, 0.76, 0.86);
    alsol.SetLineStyle(1); alsol.SetLineWidth(3); alsol.SetLineColor(kGray+2);
    alsol.Draw();
    TLatex alts; alts.SetNDC(); alts.SetTextSize(0.065); alts.SetTextAlign(12);
    alts.DrawLatex(0.765, 0.84, "Fitted model");

    // TLegend for all config entries (two rows per config)
    int nLE    = (int)legEntries.size();
    int nConfs = nLE / 2;  // each config contributes 2 rows
    int nCols  = (nConfs > 4) ? 3 : (nConfs > 2) ? 2 : 1;
    TLegend* legBot = new TLegend(0.01, 0.01, 0.99, 0.72);
    legBot->SetBorderSize(1); legBot->SetFillStyle(1001); legBot->SetFillColor(0);
    legBot->SetNColumns(nCols); legBot->SetTextSize(0.065);
    legBot->SetHeader("Config entries (P1 in cm)", "C");
    for (int i = 0; i < nLE; ++i) {
        if (legColors[i] == -1) {
            // Continuation row: text only, no marker
            legBot->AddEntry((TObject*)nullptr, legEntries[i].c_str(), "");
        } else {
            TLine* ld = new TLine();
            ld->SetLineColor(legColors[i]); ld->SetLineWidth(3);
            legBot->AddEntry(ld, legEntries[i].c_str(), "l");
        }
    }
    legBot->Draw();
    pLeg->Update();

    c->Update();
    std::string out = plotDir + "/KLong_fit_acceptance.png";
    c->SaveAs(out.c_str());
    std::string outRoot = plotDir + "/KLong_fit_acceptance.root";
    {
        TFile* fOut = TFile::Open(outRoot.c_str(),"RECREATE");
        if (fOut && !fOut->IsZombie()) {
            c->Write();
            for (TH1F* h : hToSave) h->Write();
            fOut->Close(); delete fOut;
        }
    }
    for (TH1F* h : hToSave) delete h;
    std::cout << "\nSaved: " << out << std::endl;
    std::cout << "Saved: " << outRoot << std::endl;
    delete c;
}

// ============================================================================
// KLong_fit_resolution.C
//
// PURPOSE:
//   For every detector configuration found in a given archive directory,
//   reconstruct the Gaussian FWHM of delta_p = (reco_p - true_p) as a
//   function of true momentum, then fit a two-term tracking-resolution
//   model:
//
//       FWHM(p) = 2.355 * sqrt( a^2 + (b * p)^2 )
//
//   where 'a' captures the multiple-scattering plateau and 'b' captures the
//   hit-position-resolution growth.  The pizza position P1 appears as a
//   numerical constant in the per-config legend entry.
//
//   Output:
//     <archive_dir>/plots/KLong_fit_resolution.png
//       Dashed = binned FWHM data, Solid = fitted model, same colour per
//       config.  General equation shown as TLatex.
//
// USAGE:
//   root -l -b -q 'KLong_fit_resolution.C+("/path/to/archive")'
//   (called automatically by run_pizza_fit.sh)
//
// INPUT TREE:  kaonVectors  (branches: reco_p, true_p  — vector<double>)
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
// HELPERS  (all prefixed res_ to avoid name clashes with acceptance macro)
// ============================================================================

double res_extractPos(const std::string& label, const std::string& det) {
    std::string srch = det + "-";
    size_t pos = label.find(srch);
    if (pos == std::string::npos) return -1.0;
    pos += srch.size();
    size_t end = label.find_first_not_of("0123456789", pos);
    return std::stod(label.substr(pos, end - pos));
}

struct ResDetLayout {
    std::vector<double> trackers, pizzas, fris;
    double end_detector;
    ResDetLayout() : end_detector(0.) {}
};

ResDetLayout res_parseConfig(const std::string& config) {
    ResDetLayout layout;
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

int res_getColor(double t1, double t1Lo, double t1Hi) {
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

void res_drawLayout(const std::vector<std::string>& labels,
                    const std::vector<int>& colours) {
    double mn = 1e9, mx = -1e9;
    for (int i = 0; i < (int)labels.size(); ++i) {
        ResDetLayout l = res_parseConfig(labels[i]);
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
        ResDetLayout l = res_parseConfig(labels[i]);
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

void KLong_fit_resolution(const char* archiveDirC = "") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    std::string archiveDir(archiveDirC);
    if (archiveDir.empty()) {
        std::cerr << "Usage: root -l -b -q 'KLong_fit_resolution.C+(\"/path/to/archive\")'" << std::endl;
        return;
    }
    while (!archiveDir.empty() && archiveDir.back()=='/') archiveDir.pop_back();
    std::string plotDir = archiveDir + "/plots";
    gSystem->mkdir(plotDir.c_str(), kTRUE);

    // Collect files
    std::string cmd = "find "+archiveDir+" -maxdepth 1 -name '*combined_vectors.root' 2>/dev/null | sort";
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
        size_t sfx = base.find("_combined_vectors.root");
        labels.push_back(sfx!=std::string::npos ? base.substr(0,sfx) : base);
    }
    pclose(pipe);
    if (files.empty()) { std::cerr << "No *combined_vectors.root in " << archiveDir << std::endl; return; }
    std::cout << "Found " << (int)files.size() << " config(s)." << std::endl;

    // Momentum bins (same as KLong_plot_gaussian_spreads macros)
    const int NB = 18;
    double pB[NB+1] = {0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,
                       3.0,3.5,4.0,4.5,5.0,5.5,6.0,7.0,10.0};
    const double PMIN = 0.5, PMAX = 10.0, ACUT = 1.0;
    const int    MINENT = 30;
    // Display window (frame x-axis and visible data points)
    const double DISP_PMIN = 1.0, DISP_PMAX = 4.0;
    // Fit window for pol0+expo TF1 (supervisor: use 0-5)
    const double FIT_PMIN = 0.0, FIT_PMAX = 5.0;

    // Three-pad canvas: main plot (top-left), layout (top-right), legend (bottom)
    TCanvas* c = new TCanvas("cResfit","KLong Resolution Fit",1600,960);
    TPad* pPlot   = new TPad("pPlot",  "", 0.00, 0.30, 0.72, 1.00);
    TPad* pLayout = new TPad("pLayout","", 0.72, 0.30, 1.00, 1.00);
    TPad* pLeg    = new TPad("pLeg",   "", 0.00, 0.00, 1.00, 0.30);
    pPlot->SetLeftMargin(0.13); pPlot->SetBottomMargin(0.13); pPlot->SetGrid(1,1);
    pLeg->SetLeftMargin(0.01);  pLeg->SetRightMargin(0.01);
    pLeg->SetTopMargin(0.02);   pLeg->SetBottomMargin(0.02);
    pPlot->Draw(); pLayout->Draw(); pLeg->Draw();
    pPlot->cd();

    // Frame histogram — x-axis spans the display range only
    TH1D* hFr = new TH1D("hResFr",
        "Momentum Resolution: Data (dashed) vs Model (solid)"
        ";True Kaon Momentum [GeV/c];FWHM(#Deltap) [GeV/c]",
        100, DISP_PMIN, DISP_PMAX);
    hFr->Draw("AXIS");  // y-range set after loop once data max is known

    // Accumulated during loop; all drawn in bottom panel after
    std::vector<std::string> legEntries;
    std::vector<int>         legColors;

    std::vector<std::string> layoutLabels;
    std::vector<int>         layoutColours;
    double globalYmax = 0.;
    std::vector<TH1F*> hToSave;  // TH1Fs written into the .root output

    // Pre-compute T1 range so the gradient spans the actual configs in this archive
    double t1Lo = 1e9, t1Hi = -1e9;
    for (auto& lbl : labels) {
        double t1 = res_extractPos(lbl,"T1");
        if (t1 > 0.) { t1Lo=std::min(t1Lo,t1); t1Hi=std::max(t1Hi,t1); }
    }
    if (t1Lo > t1Hi) { t1Lo = 240.; t1Hi = 480.; }  // fallback

    // Process each file
    for (int idx = 0; idx < (int)files.size(); ++idx) {
        const std::string& fp  = files[idx];
        const std::string& lbl = labels[idx];
        double p1  = res_extractPos(lbl,"P1");
        double t1  = res_extractPos(lbl,"T1");
        int    col = res_getColor(t1, t1Lo, t1Hi);

        std::cout << "\n[" << idx+1 << "/" << (int)files.size() << "] " << lbl << "  P1=" << p1 << std::endl;

        TFile* tf = TFile::Open(fp.c_str(),"READ");
        if (!tf || tf->IsZombie()) {
            std::cerr << "  WARN: cannot open " << fp << std::endl;
            if (tf) delete tf;
            continue;
        }
        TTree* tree = (TTree*)tf->Get("kaonVectors");
        if (!tree) {
            std::cerr << "  WARN: no kaonVectors tree" << std::endl;
            tf->Close(); delete tf; continue;
        }

        std::vector<double>* reco_p = nullptr;
        std::vector<double>* true_p = nullptr;
        tree->SetBranchAddress("reco_p", &reco_p);
        tree->SetBranchAddress("true_p", &true_p);

        // Per-bin histograms of delta_p.
        // SetDirectory(nullptr) detaches them from the TFile so they survive
        // tf->Close() and ROOT does not attempt to delete them when the file
        // is closed.
        std::vector<TH1D*> hb(NB);
        for (int i = 0; i < NB; ++i) {
            hb[i] = new TH1D(Form("rh_%d_%d",idx,i),"",100,-3.,3.);
            hb[i]->SetDirectory(nullptr);
        }

        Long64_t nEnt = tree->GetEntries();
        std::cout << "  Entries: " << nEnt << std::endl;
        for (Long64_t ev = 0; ev < nEnt; ++ev) {
            tree->GetEntry(ev);
            if (!reco_p || !true_p) continue;
            for (int k = 0; k < (int)reco_p->size(); ++k) {
                double tp = (*true_p)[k];
                double dp = (*reco_p)[k] - tp;
                if (tp <= 0.) continue;
                if (std::abs(dp)/tp > ACUT) continue;
                for (int b = 0; b < NB; ++b) {
                    if (tp >= pB[b] && tp < pB[b+1]) { hb[b]->Fill(dp); break; }
                }
            }
        }
        tf->Close(); delete tf;

        // Gaussian fit per bin → FWHM
        // Use heap TF1 and remove from histogram list before deleting, to
        // prevent ROOT attempting to delete a pointer it does not own.
        std::vector<double> gx, gy, gex, gey;
        for (int i = 0; i < NB; ++i) {
            if (hb[i]->GetEntries() < MINENT) { delete hb[i]; continue; }
            double mean = hb[i]->GetMean();
            double rms  = hb[i]->GetRMS();
            TF1* gf = new TF1(Form("rgf_%d_%d",idx,i),"gaus",mean-2.*rms,mean+2.*rms);
            gf->SetParameters(hb[i]->GetMaximum(), mean, rms);
            hb[i]->Fit(gf,"RQN");
            double sig  = std::abs(gf->GetParameter(2));
            double sige = gf->GetParError(2);
            // Enforce a minimum relative error of 5% so that high-statistics
            // high-p bins do not completely dominate the fit at the expense of
            // the low-p region.
            const double MIN_REL_ERR = 0.20;
            double fwhm    = 2.3548 * sig;
            double fwhm_e  = 2.3548 * (sige > 0. ? sige : 1e-6);
            if (fwhm_e < MIN_REL_ERR * fwhm) fwhm_e = MIN_REL_ERR * fwhm;
            double mid = 0.5*(pB[i]+pB[i+1]);
            gx.push_back( mid );
            gy.push_back( fwhm );
            gex.push_back( 0.5*(pB[i+1]-pB[i]) );
            gey.push_back( fwhm_e );
            // Only track y-max within the display range
            if (mid >= DISP_PMIN && mid <= DISP_PMAX && fwhm > globalYmax) globalYmax = fwhm;
            // Detach gf from histogram before deleting either object
            hb[i]->GetListOfFunctions()->Remove(gf);
            delete gf;
            delete hb[i];
        }

        if ((int)gx.size() < 3) {
            std::cerr << "  Too few bins for fit, skipping." << std::endl;
            continue;
        }

        // Build hResData TH1F (variable-bin, FWHM values + errors)
        // and fit it with pol0(0)+expo(1) using "RL" as per supervisor instruction.
        int nb = (int)gx.size();
        std::vector<double> edges;
        for (int i=0;i<nb;++i) edges.push_back(gx[i]-gex[i]);
        edges.push_back(gx[nb-1]+gex[nb-1]);
        TH1F* hd = new TH1F(Form("hResData_P1%.0f",p1),
            Form("FWHM data P1=%.0f cm;True Kaon Momentum [GeV/c];FWHM(#Deltap) [GeV/c]",p1),
            nb, edges.data());
        hd->SetDirectory(nullptr);
        for (int i=0;i<nb;++i) { hd->SetBinContent(i+1,gy[i]); hd->SetBinError(i+1,gey[i]); }
        hd->SetStats(0);

        // Fit MODEL:  FWHM(p) = a0 + exp(a1 + a2*p)  (pol0+expo range 0-5, "RL")
        TF1* fm = new TF1(Form("rfm_%d",idx), "pol0(0)+expo(1)", FIT_PMIN, FIT_PMAX);
        fm->SetParameter(0, 0.05);   // constant floor
        fm->SetParameter(1, -0.5);   // log-amplitude of expo
        fm->SetParameter(2, -0.5);   // expo decay slope
        TFitResultPtr rfit = hd->Fit(fm,"SRLN");  // S=save, R=range, L=likelihood, N=no draw
        hd->GetListOfFunctions()->Remove(fm);      // detach fm so we control its lifetime

        // Draw DATA as dashed TGraphErrors — only points within display range
        std::vector<double> dgx, dgy, dgex, dgey;
        for (int i=0;i<nb;++i) {
            if (gx[i] >= DISP_PMIN && gx[i] <= DISP_PMAX) {
                dgx.push_back(gx[i]); dgy.push_back(gy[i]);
                dgex.push_back(gex[i]); dgey.push_back(gey[i]);
            }
        }
        if (!dgx.empty()) {
            TGraphErrors* gd = new TGraphErrors((int)dgx.size(),
                &dgx[0], &dgy[0], &dgex[0], &dgey[0]);
            gd->SetLineColor(col); gd->SetLineWidth(2); gd->SetLineStyle(2);
            gd->SetMarkerColor(col); gd->SetMarkerStyle(20); gd->SetMarkerSize(0.9);
            gd->Draw("LP same");
        }

        if ((int)rfit == 0 || rfit->IsValid()) {
            double a0  = fm->GetParameter(0);
            double a1  = fm->GetParameter(1);
            double a2  = fm->GetParameter(2);
            double ndf = fm->GetNDF();
            double c2n = (ndf > 0) ? fm->GetChisquare()/ndf : -1.;
            std::cout << "  Fit OK: a0=" << a0 << "  a1=" << a1
                      << "  a2=" << a2 << "  chi2/ndf=" << c2n << std::endl;
            fm->SetLineColor(col); fm->SetLineWidth(2); fm->SetLineStyle(1);
            fm->Draw("same");
            // Persist histograms
            hToSave.push_back(hd);
            TH1F* hf = new TH1F(Form("hResFit_P1%.0f",p1),
                Form("FWHM pol0+expo fit P1=%.0f cm;True Kaon Momentum [GeV/c];FWHM(#Deltap) [GeV/c]",p1),
                100, FIT_PMIN, FIT_PMAX);
            hf->SetDirectory(nullptr);
            for (int i=1;i<=100;++i) hf->SetBinContent(i, fm->Eval(hf->GetBinCenter(i)));
            hToSave.push_back(hf);
            // Two-row legend entry: row 1 gets the coloured line marker
            legEntries.push_back(Form("P1=%.0f   a_{0}=%.3f   a_{1}=%.3f", p1, a0, a1));
            legColors.push_back(col);
            // Row 2: text only (sentinel colour -1 = no marker)
            legEntries.push_back(Form("      a_{2}=%.3f   #chi^{2}/ndf=%.2f", a2, c2n));
            legColors.push_back(-1);
        } else {
            std::cerr << "  Fit failed." << std::endl;
            delete hd;
            legEntries.push_back(Form("P1=%.0f   [fit failed]", p1));
            legColors.push_back(col);
            legEntries.push_back(std::string(""));
            legColors.push_back(-1);
            delete fm;
        }
        layoutLabels.push_back(lbl);
        layoutColours.push_back(col);
    }

    // Adjust y-axis — cap at 2 GeV/c regardless of data
    if (globalYmax > 0.) {
        hFr->GetYaxis()->SetRangeUser(0, std::min(globalYmax * 1.2, 1.5));
        hFr->Draw("AXIS same");
    }
    pPlot->RedrawAxis();
    pPlot->Update();

    // Detector layout pad
    pLayout->cd();
    if (!layoutLabels.empty()) res_drawLayout(layoutLabels, layoutColours);

    // ── Bottom legend panel ─────────────────────────────────────────────────
    pLeg->cd();
    pLeg->Range(0, 0, 1, 1);

    // Equation + key line as TLatex at the top of the pad
    TLatex eqLat;
    eqLat.SetNDC(); eqLat.SetTextSize(0.09); eqLat.SetTextAlign(12);
    eqLat.DrawLatex(0.02, 0.84,
        "Model:  FWHM(p) = a_{0} + exp(a_{1} + a_{2} #upoint p)");

    // Dashed / solid key, drawn with TLine (uses pad Range 0-1 coords)
    TLine rldash(0.56, 0.86, 0.62, 0.86);
    rldash.SetLineStyle(2); rldash.SetLineWidth(3); rldash.SetLineColor(kGray+2);
    rldash.Draw();
    TLatex rltd; rltd.SetNDC(); rltd.SetTextSize(0.08); rltd.SetTextAlign(12);
    rltd.DrawLatex(0.635, 0.84, "Binned data");
    TLine rlsol(0.76, 0.86, 0.82, 0.86);
    rlsol.SetLineStyle(1); rlsol.SetLineWidth(3); rlsol.SetLineColor(kGray+2);
    rlsol.Draw();
    TLatex rlts; rlts.SetNDC(); rlts.SetTextSize(0.08); rlts.SetTextAlign(12);
    rlts.DrawLatex(0.825, 0.84, "Fitted model");

    // TLegend for all config entries
    int nLE    = (int)legEntries.size();
    int nConfs = nLE / 2;  // each config contributes 2 rows
    int nCols  = (nConfs > 4) ? 3 : (nConfs > 2) ? 2 : 1;
    TLegend* legBot = new TLegend(0.01, 0.01, 0.99, 0.72);
    legBot->SetBorderSize(1); legBot->SetFillStyle(1001); legBot->SetFillColor(0);
    legBot->SetNColumns(nCols); legBot->SetTextSize(0.08);
    legBot->SetHeader("Config entries (P1 in cm)", "C");
    for (int i = 0; i < nLE; ++i) {
        if (legColors[i] == -1) {
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
    std::string out = plotDir + "/KLong_fit_resolution.png";
    c->SaveAs(out.c_str());
    std::string outRoot = plotDir + "/KLong_fit_resolution.root";
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

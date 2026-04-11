#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TText.h"
#include "TLine.h"
#include "TBox.h"
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>

// Structure to hold detector positions
struct DetectorLayout {
    std::vector<double> trackers;
    std::vector<double> pizzas;
    std::vector<double> fris;
    double end_detector;
};

// Parse configuration string to extract detector positions
DetectorLayout parseConfigString(const std::string& config) {
    DetectorLayout layout;
    std::istringstream ss(config);
    std::string token;
    
    while (std::getline(ss, token, '_')) {
        if (token.empty()) continue;
        size_t pos = token.find('-');
        if (pos == std::string::npos) continue;
        
        std::string detector_type = token.substr(0, pos);
        double position = std::stod(token.substr(pos + 1));
        
        if (position == 0) continue;
        
        if (detector_type[0] == 'T') {
            layout.trackers.push_back(position);
        } else if (detector_type[0] == 'P') {
            layout.pizzas.push_back(position);
        } else if (detector_type[0] == 'F') {
            layout.fris.push_back(position);
        } else if (detector_type[0] == 'E') {
            layout.end_detector = position;
        }
    }
    
    return layout;
}

// Draw detector layout schematic
void drawDetectorLayout(const std::vector<std::string>& config_labels,
                        double x_min, double x_max, double y_min, double y_max) {
    gPad->Range(x_min, y_min, x_max, y_max);
    
    // Find global min/max positions
    double min_pos = 1e9, max_pos = -1e9;
    for (const auto& config : config_labels) {
        DetectorLayout layout = parseConfigString(config);
        for (double pos : layout.trackers) {
            if (pos < min_pos) min_pos = pos;
            if (pos > max_pos) max_pos = pos;
        }
        for (double pos : layout.pizzas) {
            if (pos < min_pos) min_pos = pos;
            if (pos > max_pos) max_pos = pos;
        }
        for (double pos : layout.fris) {
            if (pos < min_pos) min_pos = pos;
            if (pos > max_pos) max_pos = pos;
        }
        if (layout.end_detector < min_pos) min_pos = layout.end_detector;
        if (layout.end_detector > max_pos) max_pos = layout.end_detector;
    }
    
    double range = max_pos - min_pos;
    min_pos -= range * 0.1;
    max_pos += range * 0.1;
    
    TText *title = new TText(0.5, 0.98, "Detector Layouts");
    title->SetTextAlign(23);
    title->SetTextSize(0.04);
    title->Draw();
    
    int n_configs = config_labels.size();
    double y_spacing = 0.85 / (n_configs + 1);
    
    for (size_t i = 0; i < config_labels.size(); ++i) {
        DetectorLayout layout = parseConfigString(config_labels[i]);
        double y_center = 0.90 - (i + 1) * y_spacing;
        
        TText *config_num = new TText(0.05, y_center, Form("%d", (int)i + 1));
        config_num->SetTextAlign(22);
        config_num->SetTextSize(0.045);
        config_num->SetTextFont(62);
        config_num->Draw();
        
        double x_scale = 0.80;
        double x_offset = 0.15;
        
        // Trackers (blue)
        for (double pos : layout.trackers) {
            double x = x_offset + (pos - min_pos) / (max_pos - min_pos) * x_scale;
            TBox *box = new TBox(x - 0.008, y_center - 0.015, x + 0.008, y_center + 0.015);
            box->SetFillColor(kBlue);
            box->SetLineColor(kBlue);
            box->Draw();
        }
        
        // Pizzas (red)
        for (double pos : layout.pizzas) {
            double x = x_offset + (pos - min_pos) / (max_pos - min_pos) * x_scale;
            TBox *box = new TBox(x - 0.008, y_center - 0.015, x + 0.008, y_center + 0.015);
            box->SetFillColor(kRed);
            box->SetLineColor(kRed);
            box->Draw();
        }
        
        // FRIs (magenta)
        for (double pos : layout.fris) {
            double x = x_offset + (pos - min_pos) / (max_pos - min_pos) * x_scale;
            TBox *box = new TBox(x - 0.008, y_center - 0.015, x + 0.008, y_center + 0.015);
            box->SetFillColor(kMagenta);
            box->SetLineColor(kMagenta);
            box->Draw();
        }
        
        // End detector (green)
        double x = x_offset + (layout.end_detector - min_pos) / (max_pos - min_pos) * x_scale;
        TBox *box = new TBox(x - 0.008, y_center - 0.015, x + 0.008, y_center + 0.015);
        box->SetFillColor(kGreen);
        box->SetLineColor(kGreen);
        box->Draw();
    }
    
    // Legend for detector types
    double legend_y = 0.08;
    double legend_x_start = 0.15;
    double legend_spacing = 0.18;
    
    TBox *t_box = new TBox(legend_x_start, legend_y, legend_x_start + 0.02, legend_y + 0.03);
    t_box->SetFillColor(kBlue);
    t_box->Draw();
    TText *t_text = new TText(legend_x_start + 0.03, legend_y + 0.015, "Tracker");
    t_text->SetTextSize(0.03);
    t_text->Draw();
    
    TBox *p_box = new TBox(legend_x_start + legend_spacing, legend_y, 
                           legend_x_start + legend_spacing + 0.02, legend_y + 0.03);
    p_box->SetFillColor(kRed);
    p_box->Draw();
    TText *p_text = new TText(legend_x_start + legend_spacing + 0.03, legend_y + 0.015, "Pizza");
    p_text->SetTextSize(0.03);
    p_text->Draw();
    
    TBox *f_box = new TBox(legend_x_start + 2*legend_spacing, legend_y, 
                           legend_x_start + 2*legend_spacing + 0.02, legend_y + 0.03);
    f_box->SetFillColor(kMagenta);
    f_box->Draw();
    TText *f_text = new TText(legend_x_start + 2*legend_spacing + 0.03, legend_y + 0.015, "FRI");
    f_text->SetTextSize(0.03);
    f_text->Draw();
    
    TBox *e_box = new TBox(legend_x_start + 3*legend_spacing, legend_y, 
                           legend_x_start + 3*legend_spacing + 0.02, legend_y + 0.03);
    e_box->SetFillColor(kGreen);
    e_box->Draw();
    TText *e_text = new TText(legend_x_start + 3*legend_spacing + 0.03, legend_y + 0.015, "End");
    e_text->SetTextSize(0.03);
    e_text->Draw();
}

void KLong_plot_compare_spreads() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    
    // Find combined vectors files
    std::string findCmd = "find /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/TGRAPH_TEST_20260316 -name '*combined_vectors.root' 2>/dev/null | sort";
    FILE* pipe = popen(findCmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Could not run find command" << std::endl;
        return;
    }
    
    std::vector<std::string> root_files;
    std::vector<std::string> config_labels;
    char buffer[512];
    
    while (fgets(buffer, sizeof(buffer), pipe)) {
        std::string filepath(buffer);
        filepath.erase(std::remove(filepath.begin(), filepath.end(), '\n'), filepath.end());
        if (!filepath.empty()) {
            root_files.push_back(filepath);
            
            // Extract config name from filename
            size_t lastSlash = filepath.find_last_of('/');
            std::string filename = filepath.substr(lastSlash + 1);
            // Remove _combined_vectors.root from the end
            size_t suffix_pos = filename.find("_combined_vectors.root");
            std::string config = (suffix_pos != std::string::npos) ? filename.substr(0, suffix_pos) : filename;
            config_labels.push_back(config);
        }
    }
    pclose(pipe);
    
    std::cout << "Found " << root_files.size() << " configuration(s)" << std::endl;
    
    if (root_files.empty()) {
        std::cerr << "No combined_vectors.root files found" << std::endl;
        return;
    }
    
    // Color scheme
    int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kCyan+1, kOrange+1, kViolet, kTeal-5};
    
    // Store histograms so they persist
    std::vector<TH2D*> histograms;
    std::vector<int> line_colors;
    
    // Process each configuration
    for (size_t idx = 0; idx < root_files.size() && idx < 7; ++idx) {
        std::string filepath = root_files[idx];
        std::string config = config_labels[idx];
        
        std::cout << "Processing: " << config << std::endl;
        
        TFile *file = TFile::Open(filepath.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            continue;
        }
        
        TTree *tree = (TTree*)file->Get("kaonVectors");
        if (!tree) {
            std::cerr << "Error: Could not find kaonVectors tree" << std::endl;
            file->Close();
            continue;
        }
        
        // Create scatter plot
        std::string h_name = Form("h_scatter_%d", (int)idx);
        TH2D *h_scatter = new TH2D(h_name.c_str(), 
                                   Form("Config %d;#Deltap (Reco - True) [GeV/c];True Momentum [GeV/c]", (int)idx+1),
                                   200, -2.0, 2.0, 200, 0, 10);
        
        // Variables
        std::vector<double> *reco_p = nullptr;
        std::vector<double> *true_p = nullptr;
        tree->SetBranchAddress("reco_p", &reco_p);
        tree->SetBranchAddress("true_p", &true_p);
        
        // Fill scatter plot
        const double anomaly_threshold = 1.0; // Filter |delta_p/true_p| > 100% as unphysical outliers
        Long64_t nEntries = tree->GetEntries();
        int points_filled = 0;
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            if (reco_p && true_p) {
                for (size_t i = 0; i < reco_p->size(); ++i) {
                    double p_true  = (*true_p)[i];
                    double delta_p = (*reco_p)[i] - p_true;
                    if (p_true <= 0.) continue;
                    if (std::abs(delta_p) / p_true > anomaly_threshold) continue;
                    h_scatter->Fill(delta_p, p_true);
                    points_filled++;
                }
            }
        }
        
        std::cout << "  Filled " << points_filled << " points into scatter plot" << std::endl;
        
        // Style histogram
        int color = colors[idx % 8];
        line_colors.push_back(color);
        h_scatter->SetMarkerColor(color);
        h_scatter->SetMarkerStyle(20);
        h_scatter->SetMarkerSize(0.8);
        h_scatter->SetTitle(Form("Config %d;#Deltap = p_{reco} - p_{true} [GeV/c];p_{true} [GeV/c]", (int)idx+1));
        h_scatter->GetXaxis()->SetTitleSize(0.05);
        h_scatter->GetYaxis()->SetTitleSize(0.05);
        h_scatter->GetXaxis()->SetLabelSize(0.04);
        h_scatter->GetYaxis()->SetLabelSize(0.04);
        h_scatter->SetDirectory(0); // Detach from file
        
        histograms.push_back(h_scatter);
        file->Close();
    }
    
    std::cout << "\nCreating canvas and drawing histograms..." << std::endl;
    
    // Create canvas with 4x2 grid for 7 scatter plots + 1 legend
    TCanvas *c1 = new TCanvas("c1", "Momentum Spread Comparison", 1800, 1200);
    c1->Divide(4, 2);
    
    // Draw all histograms
    for (size_t idx = 0; idx < histograms.size(); ++idx) {
        c1->cd(idx + 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        
        TH2D *h_scatter = histograms[idx];
        
        if (h_scatter->GetEntries() > 0) {
            h_scatter->Draw("COLZ"); // Draw as color map showing density
        } else {
            std::cout << "  WARNING: No entries in histogram " << idx << "!" << std::endl;
            h_scatter->Draw("AXIS");
        }
        
        // Add vertical line at x=0
        TLine *zero_line = new TLine(0, 0, 0, 10);
        zero_line->SetLineStyle(2);
        zero_line->SetLineColor(kBlack);
        zero_line->Draw("SAME");
        
        gPad->Modified();
        gPad->Update();
    }
    
    // Last panel: detector layout legend
    if (root_files.size() > 0) {
        c1->cd(8);
        drawDetectorLayout(config_labels, 0, 1, 0, 1);
    }
    
    c1->Update();
    c1->SaveAs("KLong_plot_compare_spreads.png");
    std::cout << "\n=== SUCCESS: Saved KLong_plot_compare_spreads.png ===" << std::endl;
    
    delete c1;
}

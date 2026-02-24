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
#include <map>

// Structure to hold detector positions
struct DetectorLayout {
    std::vector<double> trackers;
    std::vector<double> pizzas;
    std::vector<double> fris;
    double end_detector;
};

// Structure to hold pizza configuration
struct PizzaConfig {
    double p1;
    double p2;
    
    bool operator<(const PizzaConfig& other) const {
        if (p1 != other.p1) return p1 < other.p1;
        return p2 < other.p2;
    }
    
    bool operator==(const PizzaConfig& other) const {
        return p1 == other.p1 && p2 == other.p2;
    }
    
    std::string toString() const {
        std::ostringstream ss;
        ss << "P1-" << (int)p1 << "_P2-" << (int)p2;
        return ss.str();
    }
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

// Extract pizza configuration from config string
PizzaConfig extractPizzaConfig(const std::string& config) {
    PizzaConfig pc;
    pc.p1 = 0;
    pc.p2 = 0;
    
    std::istringstream ss(config);
    std::string token;
    
    while (std::getline(ss, token, '_')) {
        if (token.empty()) continue;
        
        size_t pos = token.find('-');
        if (pos == std::string::npos) continue;
        
        std::string detector_type = token.substr(0, pos);
        double position = std::stod(token.substr(pos + 1));
        
        if (detector_type == "P1") {
            pc.p1 = position;
        } else if (detector_type == "P2") {
            pc.p2 = position;
        }
    }
    
    return pc;
}

// Extract dataset type from filepath
std::string extractDatasetType(const std::string& filepath) {
    if (filepath.find("STANDARD_PIZZA") != std::string::npos) return "STANDARD_PIZZA";
    if (filepath.find("INVERTED_PIZZA") != std::string::npos) return "INVERTED_PIZZA";
    if (filepath.find("NO_TRACKERS") != std::string::npos) return "NO_TRACKERS";
    if (filepath.find("EARLY_TRACKERS") != std::string::npos) return "EARLY_TRACKERS";
    return "UNKNOWN";
}

// Get color for dataset type
int getDatasetColor(const std::string& dataset_type) {
    if (dataset_type == "STANDARD_PIZZA") return kBlue;
    if (dataset_type == "INVERTED_PIZZA") return kRed;
    if (dataset_type == "NO_TRACKERS") return kGreen+2;
    if (dataset_type == "EARLY_TRACKERS") return kMagenta;
    return kBlack;
}

// Get dataset type order for sorting
int getDatasetOrder(const std::string& dataset_type) {
    if (dataset_type == "STANDARD_PIZZA") return 0;
    if (dataset_type == "INVERTED_PIZZA") return 1;
    if (dataset_type == "NO_TRACKERS") return 2;
    if (dataset_type == "EARLY_TRACKERS") return 3;
    return 4;
}

// Draw detector layout schematic
void drawDetectorLayout(const std::vector<std::string>& config_labels,
                        const std::vector<int>& line_colors,
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
        
        TLine *indicator = new TLine(0.02, y_center, 0.08, y_center);
        indicator->SetLineColor(line_colors[i]);
        indicator->SetLineWidth(4);
        indicator->Draw();
        
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

void KLong_plot_compare_spreads_by_pizza() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    
    // Create output directory
    std::string output_dir = "spreads_by_pizza";
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    
    // Find combined vectors files
    std::string findCmd = "find /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS -name '*combined_vectors.root' 2>/dev/null | sort";
    FILE* pipe = popen(findCmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Could not run find command" << std::endl;
        return;
    }
    
    std::vector<std::string> root_files;
    char buffer[512];
    
    while (fgets(buffer, sizeof(buffer), pipe)) {
        std::string filepath(buffer);
        filepath.erase(std::remove(filepath.begin(), filepath.end(), '\n'), filepath.end());
        if (!filepath.empty()) {
            root_files.push_back(filepath);
        }
    }
    pclose(pipe);
    
    std::cout << "Found " << root_files.size() << " combined vectors files" << std::endl;
    
    if (root_files.empty()) {
        std::cerr << "No combined_vectors.root files found" << std::endl;
        return;
    }
    
    // Group files by pizza configuration
    std::map<PizzaConfig, std::vector<std::tuple<std::string, std::string, std::string>>> pizza_groups;
    
    for (const auto& fname : root_files) {
        // Extract config name from filename
        size_t lastSlash = fname.find_last_of('/');
        std::string filename = fname.substr(lastSlash + 1);
        // Remove _combined_vectors.root from the end
        size_t suffix_pos = filename.find("_combined_vectors.root");
        std::string config = (suffix_pos != std::string::npos) ? filename.substr(0, suffix_pos) : filename;
        
        // Extract pizza configuration
        PizzaConfig pc = extractPizzaConfig(config);
        
        // Extract dataset type
        std::string dataset_type = extractDatasetType(fname);
        
        // Group by pizza configuration
        pizza_groups[pc].push_back(std::make_tuple(fname, config, dataset_type));
    }
    
    // Sort each pizza group by dataset type for consistent ordering
    for (auto& pizza_entry : pizza_groups) {
        std::sort(pizza_entry.second.begin(), pizza_entry.second.end(),
            [](const std::tuple<std::string, std::string, std::string>& a,
               const std::tuple<std::string, std::string, std::string>& b) {
                return getDatasetOrder(std::get<2>(a)) < getDatasetOrder(std::get<2>(b));
            });
    }
    
    std::cout << "Found " << pizza_groups.size() << " unique pizza configurations" << std::endl;
    
    // Color scheme
    int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kCyan+1, kOrange+1, kViolet, kTeal-5};
    
    // Process each pizza configuration
    for (const auto& pizza_entry : pizza_groups) {
        PizzaConfig pc = pizza_entry.first;
        std::vector<std::tuple<std::string, std::string, std::string>> files_configs = pizza_entry.second;
        
        std::cout << "\n=== Processing pizza configuration: " << pc.toString() << " ===" << std::endl;
        
        // Store histograms so they persist
        std::vector<TH2D*> histograms;
        std::vector<int> line_colors;
        std::vector<std::string> config_labels;
        
        // Process each configuration
        for (size_t idx = 0; idx < files_configs.size() && idx < 7; ++idx) {
            std::string filepath = std::get<0>(files_configs[idx]);
            std::string config = std::get<1>(files_configs[idx]);
            std::string dataset_type = std::get<2>(files_configs[idx]);
            
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
            Long64_t nEntries = tree->GetEntries();
            int points_filled = 0;
            for (Long64_t entry = 0; entry < nEntries; ++entry) {
                tree->GetEntry(entry);
                if (reco_p && true_p) {
                    for (size_t i = 0; i < reco_p->size(); ++i) {
                        double delta_p = (*reco_p)[i] - (*true_p)[i];
                        h_scatter->Fill(delta_p, (*true_p)[i]);
                        points_filled++;
                    }
                }
            }
            
            std::cout << "  Filled " << points_filled << " points into scatter plot" << std::endl;
            
            // Style histogram
            int color = getDatasetColor(dataset_type);
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
            config_labels.push_back(config);
            file->Close();
        }
        
        if (histograms.empty()) {
            std::cout << "No valid histograms for " << pc.toString() << std::endl;
            continue;
        }
        
        std::cout << "\nCreating canvas and drawing histograms..." << std::endl;
        
        // Create canvas with 4x2 grid for up to 7 scatter plots + 1 legend
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
        if (histograms.size() > 0) {
            c1->cd(8);
            drawDetectorLayout(config_labels, line_colors, 0, 1, 0, 1);
        }
        
        c1->Update();
        std::string output_filename = output_dir + "/KLong_plot_compare_spreads_" + pc.toString() + ".png";
        c1->SaveAs(output_filename.c_str());
        std::cout << "\n=== SUCCESS: Saved " << output_filename << " ===" << std::endl;
        
        // Clean up
        for (auto h : histograms) {
            delete h;
        }
        delete c1;
    }
    
    std::cout << "\n=== All pizza configurations processed ===" << std::endl;
}

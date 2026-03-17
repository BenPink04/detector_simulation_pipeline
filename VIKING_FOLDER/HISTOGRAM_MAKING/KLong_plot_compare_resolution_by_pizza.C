#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TBox.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>

// Structure to hold detector positions
struct DetectorLayout {
    std::vector<double> trackers;  // T1-T4
    std::vector<double> pizzas;    // P1-P2
    std::vector<double> fris;      // F1-F2
    double end_detector;            // E1
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
    
    // Parse format: T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700
    std::istringstream ss(config);
    std::string token;
    
    while (std::getline(ss, token, '_')) {
        if (token.empty()) continue;
        
        // Find the hyphen separator
        size_t pos = token.find('-');
        if (pos == std::string::npos) continue;
        
        std::string detector_type = token.substr(0, pos);
        double position = std::stod(token.substr(pos + 1));
        
        // Skip disabled components (position = 0)
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
    // Set up coordinate system in user coordinates
    gPad->Range(x_min, y_min, x_max, y_max);
    
    // Find global min/max positions for scaling
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
    
    // Add padding
    double range = max_pos - min_pos;
    min_pos -= range * 0.1;
    max_pos += range * 0.1;
    
    // Title (in user coordinates)
    TText *title = new TText(0.5, 0.98, "Detector Layouts (Top View)");
    title->SetTextAlign(23);
    title->SetTextSize(0.04);
    title->Draw();
    
    // Draw each configuration
    int n_configs = config_labels.size();
    double y_spacing = 0.85 / (n_configs + 1);
    
    for (size_t i = 0; i < config_labels.size(); ++i) {
        DetectorLayout layout = parseConfigString(config_labels[i]);
        double y_center = 0.90 - (i + 1) * y_spacing;
        
        // Draw configuration line indicator
        TLine *indicator = new TLine(0.02, y_center, 0.08, y_center);
        indicator->SetLineColor(line_colors[i]);
        indicator->SetLineWidth(4);
        indicator->Draw();
        
        // Draw detector components
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
        
        // FRIs (pink/magenta)
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
    
    // Add legend for detector types at bottom
    double legend_y = 0.08;
    double legend_x_start = 0.15;
    double legend_spacing = 0.18;
    
    // Tracker legend
    TBox *t_box = new TBox(legend_x_start, legend_y, legend_x_start + 0.02, legend_y + 0.03);
    t_box->SetFillColor(kBlue);
    t_box->Draw();
    TText *t_text = new TText(legend_x_start + 0.03, legend_y + 0.015, "Tracker");
    t_text->SetTextSize(0.03);
    t_text->Draw();
    
    // Pizza legend
    TBox *p_box = new TBox(legend_x_start + legend_spacing, legend_y, legend_x_start + legend_spacing + 0.02, legend_y + 0.03);
    p_box->SetFillColor(kRed);
    p_box->Draw();
    TText *p_text = new TText(legend_x_start + legend_spacing + 0.03, legend_y + 0.015, "Pizza");
    p_text->SetTextSize(0.03);
    p_text->Draw();
    
    // FRI legend
    TBox *f_box = new TBox(legend_x_start + 2*legend_spacing, legend_y, legend_x_start + 2*legend_spacing + 0.02, legend_y + 0.03);
    f_box->SetFillColor(kMagenta);
    f_box->Draw();
    TText *f_text = new TText(legend_x_start + 2*legend_spacing + 0.03, legend_y + 0.015, "FRI");
    f_text->SetTextSize(0.03);
    f_text->Draw();
    
    // End detector legend
    TBox *e_box = new TBox(legend_x_start + 3*legend_spacing, legend_y, legend_x_start + 3*legend_spacing + 0.02, legend_y + 0.03);
    e_box->SetFillColor(kGreen);
    e_box->Draw();
    TText *e_text = new TText(legend_x_start + 3*legend_spacing + 0.03, legend_y + 0.015, "End");
    e_text->SetTextSize(0.03);
    e_text->Draw();
}

void KLong_plot_compare_resolution_by_pizza() {
    // Create output directory
    std::string output_dir = "resolution_by_pizza";
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    
    // Find all subdirectories in ARCHIVED_RESULTS and all combined vectors files
    std::string findCmd = "find /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/TGRAPH_TEST_20260316 -name '*combined_vectors.root' 2>/dev/null | sort";
    FILE* pipe = popen(findCmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Could not run find command" << std::endl;
        return;
    }
    
    std::vector<std::string> filenames;
    char buffer[512];
    
    while (fgets(buffer, sizeof(buffer), pipe)) {
        std::string filepath(buffer);
        filepath.erase(std::remove(filepath.begin(), filepath.end(), '\n'), filepath.end());
        if (!filepath.empty()) {
            filenames.push_back(filepath);
        }
    }
    pclose(pipe);
    
    std::cout << "Found " << filenames.size() << " combined vectors files across all directories" << std::endl;
    
    if (filenames.empty()) {
        std::cout << "No combined vectors files found" << std::endl;
        return;
    }

    // Extract configuration labels and group by pizza configuration
    std::map<PizzaConfig, std::vector<std::tuple<std::string, std::string, std::string>>> pizza_groups;
    
    for (const auto& fname : filenames) {
        // Extract config name from filename
        size_t last_slash = fname.find_last_of('/');
        std::string basename = (last_slash != std::string::npos) ? fname.substr(last_slash + 1) : fname;
        
        // Remove _combined_vectors.root from the end
        size_t suffix_pos = basename.find("_combined_vectors.root");
        std::string label = (suffix_pos != std::string::npos) ? basename.substr(0, suffix_pos) : basename;
        
        // Extract pizza configuration
        PizzaConfig pc = extractPizzaConfig(label);
        
        // Extract dataset type
        std::string dataset_type = extractDatasetType(fname);
        
        // Group by pizza configuration
        pizza_groups[pc].push_back(std::make_tuple(fname, label, dataset_type));
    }
    
    // Sort each pizza group by dataset type for consistent ordering
    for (auto& pizza_entry : pizza_groups) {
        std::sort(pizza_entry.second.begin(), pizza_entry.second.end(),
            [](const std::tuple<std::string, std::string, std::string>& a,
               const std::tuple<std::string, std::string, std::string>& b) {
                return getDatasetOrder(std::get<2>(a)) < getDatasetOrder(std::get<2>(b));
            });
    }

    std::cout << "\nFound " << pizza_groups.size() << " unique pizza configurations:" << std::endl;
    for (const auto& entry : pizza_groups) {
        std::cout << "  " << entry.first.toString() << ": " << entry.second.size() << " files" << std::endl;
    }

    // Histogram parameters
    int nbins = 40;
    double min_p = 0.0;
    double max_p = 10.0;
    double anomaly_threshold = 1.0; // Ignore resolution > 1 (100%)

    // Colors for different configurations
    std::vector<int> colors = {kBlue, kRed, kGreen+2, kMagenta, kCyan+2, kOrange+2, kViolet, kGray+2};
    
    // Process each pizza configuration
    for (const auto& pizza_entry : pizza_groups) {
        PizzaConfig pc = pizza_entry.first;
        std::vector<std::tuple<std::string, std::string, std::string>> files_labels = pizza_entry.second;
        
        std::cout << "\n=== Processing pizza configuration: " << pc.toString() << " ===" << std::endl;
        
        // Create canvas
        TCanvas *c1 = new TCanvas("c1", "Kaon Momentum Resolution Comparison", 1400, 700);
        c1->Divide(2, 1);  // Divide canvas: left for histogram, right for legend
        c1->cd(1);         // Switch to left pad for histogram
        gPad->SetPad(0.0, 0.0, 0.75, 1.0);  // Left pad takes 75% of width
        c1->cd(2);         // Switch to right pad for legend
        gPad->SetPad(0.75, 0.0, 1.0, 1.0);  // Right pad takes 25% of width
        c1->cd(1);         // Go back to left pad for drawing histogram
        
        // Vector to store histograms
        std::vector<TH1D*> h_mean_vec;
        std::vector<std::string> config_labels;
        
        // Process each file for this pizza configuration
        for (size_t file_idx = 0; file_idx < files_labels.size(); ++file_idx) {
            const std::string& fname = std::get<0>(files_labels[file_idx]);
            const std::string& config_label = std::get<1>(files_labels[file_idx]);
            const std::string& dataset_type = std::get<2>(files_labels[file_idx]);
            
            std::cout << "Processing file " << (file_idx + 1) << "/" << files_labels.size() 
                      << ": " << config_label << std::endl;
            
            // Create 2D histogram for this configuration
            std::string h2_name = "h2_" + std::to_string(file_idx);
            std::string h2_title = "Resolution vs True Kaon Momentum (" + config_label + ")";
            h2_title += ";True Kaon Momentum [GeV/c];|Reco - True|/True";
            
            TH2D *h2 = new TH2D(h2_name.c_str(), h2_title.c_str(), nbins, min_p, max_p, nbins, 0, 1);
            
            // Fill histogram with data from this file
            TFile *inFile = TFile::Open(fname.c_str());
            if (!inFile || inFile->IsZombie()) {
                std::cerr << "  Error: Could not open file: " << fname << std::endl;
                continue;
            }

            TTree *tree = (TTree*)inFile->Get("kaonVectors");
            if (!tree) {
                std::cerr << "  Error: Could not find kaonVectors tree" << std::endl;
                inFile->Close();
                continue;
            }

            std::vector<double> *reco_p = nullptr;
            std::vector<double> *true_p = nullptr;
            tree->SetBranchAddress("reco_p", &reco_p);
            tree->SetBranchAddress("true_p", &true_p);
            
            // Loop through ALL entries in the combined file
            Long64_t nEntries = tree->GetEntries();
            std::cout << "  Processing " << nEntries << " entries..." << std::endl;
            
            int total_events = 0;
            for (Long64_t entry = 0; entry < nEntries; ++entry) {
                tree->GetEntry(entry);
                if (reco_p && true_p) {
                    for (size_t i = 0; i < reco_p->size(); ++i) {
                        double rel_resolution = std::abs((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i];
                        if (rel_resolution < anomaly_threshold) {
                            h2->Fill((*true_p)[i], rel_resolution);
                            total_events++;
                        }
                    }
                }
            }
            
            std::cout << "  Total kaon events processed: " << total_events << std::endl;
            inFile->Close();

            // Create 1D histogram of mean resolution per bin
            std::string h_mean_name = "h_mean_" + std::to_string(file_idx);
            std::string h_mean_title = "Mean Relative Kaon Momentum Resolution";
            h_mean_title += ";True Kaon Momentum [GeV/c];Mean(|Reco - True|/True)";
            
            TH1D *h_mean = new TH1D(h_mean_name.c_str(), h_mean_title.c_str(), nbins, min_p, max_p);
            
            for (int i = 1; i <= nbins; ++i) {
                TH1D *proj = h2->ProjectionY("_py", i, i);
                if (proj->GetEntries() > 0) {
                    h_mean->SetBinContent(i, proj->GetMean());
                }
                delete proj;
            }
            
            // Style the histogram as a line graph
            int color = getDatasetColor(dataset_type);
            h_mean->SetLineColor(color);
            h_mean->SetLineWidth(3);
            h_mean->SetMarkerColor(color);
            h_mean->SetMarkerStyle(20 + file_idx); // Different marker style for each
            h_mean->SetMarkerSize(1.2);
            
            h_mean_vec.push_back(h_mean);
            config_labels.push_back(config_label);
            
            // Clean up 2D histogram
            delete h2;
        }
        
        if (h_mean_vec.empty()) {
            std::cout << "Error: No valid histograms created for " << pc.toString() << std::endl;
            delete c1;
            continue;
        }
        
        // Draw histograms as line graphs
        gStyle->SetOptStat(0); // Turn off statistics box
        
        // Set up the first histogram
        std::string plot_title = "Kaon Momentum Resolution Comparison (" + pc.toString() + ")";
        h_mean_vec[0]->SetTitle(plot_title.c_str());
        h_mean_vec[0]->GetYaxis()->SetRangeUser(0, 0.5); // Adjust range as needed
        h_mean_vec[0]->Draw("HIST L"); // Draw as histogram line
        
        // Draw additional histograms
        for (size_t i = 1; i < h_mean_vec.size(); ++i) {
            h_mean_vec[i]->Draw("HIST L SAME"); // Draw as histogram line on same axes
        }
        
        // Switch to right pad for detector layout visualization
        c1->cd(2);
        
        // Clear the pad
        gPad->Clear();
        
        // Collect line colors
        std::vector<int> line_colors;
        for (size_t i = 0; i < h_mean_vec.size(); ++i) {
            line_colors.push_back(h_mean_vec[i]->GetLineColor());
        }
        
        // Draw detector layouts
        drawDetectorLayout(config_labels, line_colors, 0, 1, 0, 1);
        
        // Switch back to left pad
        c1->cd(1);
        
        // Add grid for better readability
        gPad->SetGrid(1, 1);
        c1->Update();
        
        // Save the plot - save entire canvas to include both pads
        std::string output_filename = output_dir + "/KLong_resolution_comparison_" + pc.toString() + ".png";
        c1->SaveAs(output_filename.c_str());
        
        std::cout << "\nComparison plot saved as: " << output_filename << std::endl;
        
        // Clean up
        for (auto h : h_mean_vec) {
            delete h;
        }
        delete c1;
    }
    
    std::cout << "\n=== All pizza configurations processed ===" << std::endl;
}

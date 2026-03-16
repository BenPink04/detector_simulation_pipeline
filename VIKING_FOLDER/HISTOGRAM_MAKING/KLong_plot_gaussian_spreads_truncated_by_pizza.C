#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TText.h"
#include "TLine.h"
#include "TBox.h"
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cmath>
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

void KLong_plot_gaussian_spreads_truncated_by_pizza() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    // Create output directory
    std::string output_dir = "gaussian_spreads_truncated_by_pizza";
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    
    // Find combined vectors files
    std::string findCmd = "find /users/bp969/scratch/VIKING_FOLDER/ARCHIVED_RESULTS/STANDARD_PIZZA_NEW_4_20260313 -name '*combined_vectors.root' 2>/dev/null | sort";
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
    
    std::cout << "Found " << root_files.size() << " configuration(s)" << std::endl;
    
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
    
    // Momentum binning for Gaussian fits
    const int nBins = 18;
    double pBins[nBins+1] = {0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 10.0};
    
    // Process each pizza configuration
    for (const auto& pizza_entry : pizza_groups) {
        PizzaConfig pc = pizza_entry.first;
        std::vector<std::tuple<std::string, std::string, std::string>> files_configs = pizza_entry.second;
        
        std::cout << "\n=== Processing pizza configuration: " << pc.toString() << " ===" << std::endl;
        
        std::vector<TGraph*> graphs;
        std::vector<int> line_colors;
        std::vector<std::string> config_labels;
        
        // Process each configuration
        for (size_t idx = 0; idx < files_configs.size(); ++idx) {
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
            
            // Variables
            std::vector<double> *reco_p = nullptr;
            std::vector<double> *true_p = nullptr;
            tree->SetBranchAddress("reco_p", &reco_p);
            tree->SetBranchAddress("true_p", &true_p);
            
            // Create histograms for each momentum bin
            std::vector<TH1D*> h_bins;
            for (int i = 0; i < nBins; ++i) {
                std::string h_name = Form("h_bin_%d_%d", (int)idx, i);
                TH1D *h = new TH1D(h_name.c_str(), "", 100, -2.0, 2.0);
                h_bins.push_back(h);
            }
            
            // Fill histograms
            Long64_t nEntries = tree->GetEntries();
            for (Long64_t entry = 0; entry < nEntries; ++entry) {
                tree->GetEntry(entry);
                if (reco_p && true_p) {
                    for (size_t i = 0; i < reco_p->size(); ++i) {
                        double true_mom = (*true_p)[i];
                        double delta_p = (*reco_p)[i] - true_mom;
                        
                        // Find which bin this momentum falls into
                        for (int bin = 0; bin < nBins; ++bin) {
                            if (true_mom >= pBins[bin] && true_mom < pBins[bin+1]) {
                                h_bins[bin]->Fill(delta_p);
                                break;
                            }
                        }
                    }
                }
            }
            
            // Fit Gaussians and extract FWHM
            std::vector<double> p_centers, fwhm_values, fwhm_errors;
            
            for (int i = 0; i < nBins; ++i) {
                if (h_bins[i]->GetEntries() < 50) continue; // Skip bins with too few entries
                
                // Fit Gaussian
                double mean = h_bins[i]->GetMean();
                double rms = h_bins[i]->GetRMS();
                
                TF1 *gauss = new TF1(Form("gauss_%d_%d", (int)idx, i), "gaus", mean - 2*rms, mean + 2*rms);
                gauss->SetParameters(h_bins[i]->GetMaximum(), mean, rms);
                
                h_bins[i]->Fit(gauss, "RQN"); // R=range, Q=quiet, N=no draw
                
                double sigma = gauss->GetParameter(2);
                double sigma_err = gauss->GetParError(2);
                double fwhm = 2.355 * sigma; // FWHM = 2*sqrt(2*ln(2)) * sigma
                double fwhm_err = 2.355 * sigma_err;
                
                double p_center = (pBins[i] + pBins[i+1]) / 2.0;
                
                p_centers.push_back(p_center);
                fwhm_values.push_back(fwhm);
                fwhm_errors.push_back(fwhm_err);
                
                delete gauss;
                delete h_bins[i];
            }
            
            file->Close();
            
            // Create TGraph (without error bars)
            if (p_centers.size() > 0) {
                TGraph *graph = new TGraph(p_centers.size(), 
                                           &p_centers[0], &fwhm_values[0]);
                int color = getDatasetColor(dataset_type);
                line_colors.push_back(color);
                graph->SetLineColor(color);
                graph->SetLineWidth(3);
                graph->SetMarkerColor(color);
                graph->SetMarkerStyle(20 + idx);
                graph->SetMarkerSize(1.2);
                graphs.push_back(graph);
                config_labels.push_back(config);
            }
        }
        
        if (graphs.empty()) {
            std::cerr << "Error: No valid graphs created for " << pc.toString() << std::endl;
            continue;
        }
        
        // Create canvas
        TCanvas *c1 = new TCanvas("c1", "Gaussian FWHM vs Momentum", 1400, 700);
        c1->Divide(2, 1);
        
        // Left panel: FWHM plot
        c1->cd(1);
        gPad->SetPad(0.0, 0.0, 0.75, 1.0);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetGrid(1, 1);
        
        // Draw graphs - calculate max FWHM only for points in the 1-4 GeV/c range
        double max_fwhm = 0;
        for (auto g : graphs) {
            for (int i = 0; i < g->GetN(); ++i) {
                double x, y;
                g->GetPoint(i, x, y);
                if (x >= 1.0 && x <= 4.0 && y > max_fwhm) {
                    max_fwhm = y;
                }
            }
        }
        
        std::string plot_title = "Momentum Resolution (Gaussian FWHM) - " + pc.toString();
        graphs[0]->SetTitle((plot_title + ";True Momentum [GeV/c];FWHM(#Deltap) [GeV/c]").c_str());
        graphs[0]->GetXaxis()->SetLimits(1.0, 4.0);
        graphs[0]->GetYaxis()->SetRangeUser(0, max_fwhm * 1.2);
        graphs[0]->Draw("ALP");
        
        for (size_t i = 1; i < graphs.size(); ++i) {
            graphs[i]->Draw("LP SAME");
        }
        
        // Right panel: detector layout legend
        c1->cd(2);
        gPad->SetPad(0.75, 0.0, 1.0, 1.0);
        drawDetectorLayout(config_labels, line_colors, 0, 1, 0, 1);
        
        c1->cd(1);
        c1->Update();
        std::string output_filename = output_dir + "/KLong_plot_gaussian_spreads_truncated_" + pc.toString() + ".png";
        c1->SaveAs(output_filename.c_str());
        std::cout << "\n=== SUCCESS: Saved " << output_filename << " ===" << std::endl;
        
        // Clean up
        for (auto g : graphs) {
            delete g;
        }
        delete c1;
    }
    
    std::cout << "\n=== All pizza configurations processed ===" << std::endl;
}

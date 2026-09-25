#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TBox.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <sstream>
#include <map>

// Structure to hold detector positions
struct DetectorLayout {
    std::vector<double> trackers;  // T1-T4
    std::vector<double> pizzas;    // P1-P2
    std::vector<double> fris;      // F1-F2
    double end_detector;            // E1
    double z0;                      // Tungsten foil / KS start position
};

// Parse configuration string to extract detector positions
DetectorLayout parseConfigString(const std::string& config) {
    DetectorLayout layout;
    layout.end_detector = 0.0;
    layout.z0 = -1.0;
    
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
        } else if (detector_type == "Z0") {
            layout.z0 = position;
        }
    }
    
    return layout;
}

// Extract Z0 position from a config label string (e.g. "KS_Z0-300_T1-..." -> 300)
int extractZ0(const std::string& config) {
    size_t pos = config.find("Z0-");
    if (pos == std::string::npos) return 0;
    size_t end = config.find('_', pos + 3);
    std::string val = config.substr(pos + 3, end == std::string::npos ? std::string::npos : end - (pos + 3));
    return std::stoi(val);
}

// Assign color by interpolating in RGB along the Z0 scan.
int getConfigColor(int z0_val) {
    struct Anchor { int z0; float r, g, b; };
    const Anchor anchors[] = {
        {200,  63,  63, 255},
        {250,   0, 170, 255},
        {300,   0, 196, 196},
        {350, 255, 140,   0},
        {400, 220,  20,  60},
    };
    const int nAnchors = 5;
    if (z0_val <= anchors[0].z0)
        return TColor::GetColor(anchors[0].r/255.f, anchors[0].g/255.f, anchors[0].b/255.f);
    if (z0_val >= anchors[nAnchors-1].z0)
        return TColor::GetColor(anchors[nAnchors-1].r/255.f, anchors[nAnchors-1].g/255.f, anchors[nAnchors-1].b/255.f);
    for (int i = 0; i < nAnchors - 1; ++i) {
        if (z0_val >= anchors[i].z0 && z0_val <= anchors[i+1].z0) {
            float frac = float(z0_val - anchors[i].z0) / float(anchors[i+1].z0 - anchors[i].z0);
            float r = anchors[i].r + frac * (anchors[i+1].r - anchors[i].r);
            float g = anchors[i].g + frac * (anchors[i+1].g - anchors[i].g);
            float b = anchors[i].b + frac * (anchors[i+1].b - anchors[i].b);
            return TColor::GetColor(r/255.f, g/255.f, b/255.f);
        }
    }
    return kBlack;
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

        int z0_val = extractZ0(config_labels[i]);
        TText *z0_label = new TText(0.085, y_center, Form("Z0=%d", z0_val));
        z0_label->SetTextAlign(12);
        z0_label->SetTextSize(0.022);
        z0_label->Draw();
        
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
    TText *e_text = new TText(legend_x_start + 3*legend_spacing + 0.03, legend_y + 0.015, "ToF");
    e_text->SetTextSize(0.03);
    e_text->Draw();

}

void KLong_plot_compare_acceptance_truncated() {
    // Find combined acceptance files using system find command
    std::string findCmd = "find /users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS -name 'KS_Z0-*combined_acceptance.root' 2>/dev/null | sort";
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
    
    // Sort filenames by Z0 for consistent scan ordering
    std::sort(filenames.begin(), filenames.end(), [](const std::string& a, const std::string& b) {
        auto get_z0 = [](const std::string& s) {
            size_t p = s.find("Z0-");
            if (p == std::string::npos) return 0;
            size_t e = s.find('_', p + 3);
            return std::stoi(s.substr(p + 3, e == std::string::npos ? std::string::npos : e - (p + 3)));
        };
        return get_z0(a) < get_z0(b);
    });
    
    std::cout << "Found " << filenames.size() << " combined acceptance files:" << std::endl;
    for (const auto& f : filenames) {
        std::cout << "  " << f << std::endl;
    }
    
    if (filenames.empty()) {
        std::cout << "No combined acceptance files found" << std::endl;
        return;
    }

    // Extract configuration labels from filenames
    std::vector<std::string> config_labels;
    for (const auto& fname : filenames) {
        // Extract config name from filename
        size_t last_slash = fname.find_last_of('/');
        std::string basename = (last_slash != std::string::npos) ? fname.substr(last_slash + 1) : fname;
        
        // Remove _combined_acceptance.root from the end
        size_t suffix_pos = basename.find("_combined_acceptance.root");
        std::string label = (suffix_pos != std::string::npos) ? basename.substr(0, suffix_pos) : basename;
        
        config_labels.push_back(label);
    }

    // Histogram parameters (matching your acceptance plot)
    double p_min = 0.0, p_max = 10.0;
    double bin_width = 0.5;
    int n_bins = int((p_max - p_min) / bin_width);

    // Create canvas with divided layout
    TCanvas *c1 = new TCanvas("c1", "Kaon Reconstruction Acceptance Comparison", 1400, 700);
    c1->Divide(2, 1);  // Divide canvas: left for histogram, right for legend
    c1->cd(1);         // Switch to left pad for histogram
    gPad->SetPad(0.0, 0.0, 0.75, 1.0);  // Left pad takes 75% of width
    c1->cd(2);         // Switch to right pad for legend
    gPad->SetPad(0.75, 0.0, 1.0, 1.0);  // Right pad takes 25% of width
    c1->cd(1);         // Go back to left pad for drawing histogram
    
    // Vector to store histograms
    std::vector<TH1D*> h_acceptance_vec;
    
    // Process each file
    for (size_t file_idx = 0; file_idx < filenames.size(); ++file_idx) {
        const std::string& fname = filenames[file_idx];
        const std::string& config_label = config_labels[file_idx];
        
        std::cout << "Processing file " << (file_idx + 1) << "/" << filenames.size() 
                  << ": " << config_label << std::endl;
        
        // Initialize counters for this configuration
        std::vector<int> reco_in_bin(n_bins, 0);
        std::vector<int> total_in_bin(n_bins, 0);

        // Process the acceptance file
        TFile *file = TFile::Open(fname.c_str());
        if (!file || file->IsZombie()) {
            std::cout << "Warning: Cannot open file " << fname << std::endl;
            continue;
        }

        std::cout << "  Successfully opened file: " << fname << std::endl;
        file->ls(); // List contents for debugging

        TTree *tree = (TTree*)file->Get("kaonEventInfo");
        if (!tree) {
            std::cout << "Warning: Tree kaonEventInfo not found in " << fname << std::endl;
            file->Close();
            continue;
        }

        std::cout << "  Tree found with " << tree->GetEntries() << " entries" << std::endl;

        std::vector<double> *true_mom_vec = nullptr;
        std::vector<int> *reco_flag_vec = nullptr;
        int n_triple_pion_events = 0;
        tree->SetBranchAddress("true_mom_vec", &true_mom_vec);
        tree->SetBranchAddress("reco_flag_vec", &reco_flag_vec);
        tree->SetBranchAddress("n_triple_pion_events", &n_triple_pion_events);

        // Loop through ALL entries in the combined file
        Long64_t nEntries = tree->GetEntries();
        std::cout << "  Processing " << nEntries << " entries..." << std::endl;
        
        int total_events = 0;
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            
            // Process all events in this entry
            for (int i = 0; i < n_triple_pion_events; ++i) {
                double p = (*true_mom_vec)[i];
                int bin = int((p - p_min) / bin_width);
                if (bin < 0 || bin >= n_bins) continue;
                
                total_in_bin[bin]++;
                if ((*reco_flag_vec)[i]) {
                    reco_in_bin[bin]++;
                }
                total_events++;
            }
        }
        
        std::cout << "  Total kaon events processed: " << total_events << std::endl;
        file->Close();

        // Create acceptance histogram for this configuration
        std::string h_name = "h_acceptance_" + std::to_string(file_idx);
        std::string h_title = "Reconstruction Acceptance Comparison";
        h_title += ";True Momentum [GeV/c];Acceptance (Reconstructed/Total)";
        
        TH1D *h_acceptance = new TH1D(h_name.c_str(), h_title.c_str(), n_bins, p_min, p_max);

        // Fill acceptance histogram
        for (int bin = 0; bin < n_bins; ++bin) {
            int reconstructed_vertices = reco_in_bin[bin];
            int total_vertices = total_in_bin[bin];
            
            double acceptance;
            if (total_vertices > 0) {
                acceptance = double(reconstructed_vertices) / double(total_vertices);
            } else {
                acceptance = 0.0; // No vertices in this bin
            }
        
            
            h_acceptance->SetBinContent(bin + 1, acceptance);
        }
        
        // Style the histogram as a line graph
        int color = getConfigColor(extractZ0(config_label));
        h_acceptance->SetLineColor(color);
        h_acceptance->SetLineWidth(3);
        h_acceptance->SetMarkerColor(color);
        h_acceptance->SetMarkerStyle(20 + file_idx); // Different marker style for each
        h_acceptance->SetMarkerSize(1.2);
        
        h_acceptance_vec.push_back(h_acceptance);
    }
    
    if (h_acceptance_vec.empty()) {
        std::cout << "Error: No valid histograms created!" << std::endl;
        return;
    }
    
    // Draw histograms as line graphs
    gStyle->SetOptStat(0); // Turn off statistics box
    
    // Find maximum value across all histograms within 1-4 GeV/c to set y-axis range
    double max_acceptance = 0;
    for (size_t i = 0; i < h_acceptance_vec.size(); ++i) {
        TH1D *h = h_acceptance_vec[i];
        for (int b = h->FindBin(1.0); b <= h->FindBin(4.0); ++b) {
            double val = h->GetBinContent(b);
            if (val > max_acceptance) max_acceptance = val;
        }
    }
    
    // Set up the first histogram
    h_acceptance_vec[0]->SetTitle("Kaon Reconstruction Acceptance Comparison");
    h_acceptance_vec[0]->GetXaxis()->SetRangeUser(1.0, 4.0);
    h_acceptance_vec[0]->SetMinimum(0.0);
    double y_max = (max_acceptance > 0.0) ? (max_acceptance * 1.10) : 1.0;
    h_acceptance_vec[0]->SetMaximum(y_max);
    h_acceptance_vec[0]->Draw("HIST L"); // Draw as histogram line
    
    // Draw additional histograms
    for (size_t i = 1; i < h_acceptance_vec.size(); ++i) {
        h_acceptance_vec[i]->Draw("HIST L SAME"); // Draw as histogram line on same axes
    }
    
    // Switch to right pad for detector layout visualization
    c1->cd(2);
    
    // Clear the pad
    gPad->Clear();
    
    // Collect line colors
    std::vector<int> line_colors;
    for (size_t i = 0; i < h_acceptance_vec.size(); ++i) {
        line_colors.push_back(h_acceptance_vec[i]->GetLineColor());
    }
    
    // Draw detector layouts
    drawDetectorLayout(config_labels, line_colors, 0, 1, 0, 1);
    
    // Switch back to left pad
    c1->cd(1);
    
    // Add grid for better readability
    c1->SetGrid(1, 1);
    c1->Update();
    
    // Save the plot - save entire canvas to include both pads
    std::string output_filename = "KLong_acceptance_comparison_truncated.png";
    c1->SaveAs(output_filename.c_str());
    
    std::cout << "\nComparison plot saved as: " << output_filename << std::endl;
    
    // Print summary statistics
    std::cout << "\nSummary Statistics:" << std::endl;
    std::cout << "Configuration\t\tMean Acceptance (0.5-2.0 GeV)\tMean Acceptance (2.0-5.0 GeV)" << std::endl;
    std::cout << "-------------\t\t-----------------------------\t-----------------------------" << std::endl;
    
    for (size_t i = 0; i < h_acceptance_vec.size(); ++i) {
        TH1D* h = h_acceptance_vec[i];
        
        // Calculate mean acceptance in different momentum ranges
        double sum_low = 0, count_low = 0;
        double sum_high = 0, count_high = 0;
        
        for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
            double momentum = h->GetBinCenter(bin);
            double acceptance = h->GetBinContent(bin);
            
            if (acceptance > 0) {
                if (momentum >= 0.5 && momentum <= 2.0) {
                    sum_low += acceptance;
                    count_low++;
                } else if (momentum >= 2.0 && momentum <= 5.0) {
                    sum_high += acceptance;
                    count_high++;
                }
            }
        }
        
        double mean_low = (count_low > 0) ? sum_low / count_low : 0;
        double mean_high = (count_high > 0) ? sum_high / count_high : 0;
        
        std::cout << config_labels[i] << "\t\t" 
                  << std::fixed << std::setprecision(4) << mean_low << "\t\t\t\t"
                  << mean_high << std::endl;
    }
}
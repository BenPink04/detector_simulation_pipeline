#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>

void KLong_plot_compare_resolution() {
    // Define multiple files to compare - add your file paths here
    std::vector<std::string> filenames = {
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-260_F2-270_E1-600_combined_vectors.root",
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700/T1-240_T2-250_T3-680_T4-690_P1-215_P2-230_F1-260_F2-270_E1-700_combined_vectors.root",
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600/T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600_combined_vectors.root",
        "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-400_F2-410_E1-600/T1-240_T2-250_T3-580_T4-590_P1-215_P2-230_F1-400_F2-410_E1-600_combined_vectors.root"
        // Add more file paths here as needed
    };

    // Extract configuration labels from filenames
    std::vector<std::string> config_labels;
    for (const auto& fname : filenames) {
        std::string label = "";
        
        // Find the last occurrence of the configuration pattern (after the last '/')
        size_t last_slash = fname.find_last_of('/');
        std::string basename = (last_slash != std::string::npos) ? fname.substr(last_slash + 1) : fname;
        
        size_t start = basename.find("T1-");
        size_t end = basename.find("_combined_vectors.root");
        if (start != std::string::npos && end != std::string::npos) {
            label = basename.substr(start, end - start);
            // Keep full configuration label (all detector positions)
            // No truncation - includes T1-T4, P1-P2, F1-F2, E1 positions
        } else {
            label = "Config_" + std::to_string(config_labels.size() + 1);
        }
        
        config_labels.push_back(label);
    }

    // Histogram parameters
    int nbins = 40;
    double min_p = 0.0;
    double max_p = 10.0;
    double anomaly_threshold = 1.0; // Ignore resolution > 1 (100%)

    // Colors for different configurations
    std::vector<int> colors = {kBlue, kRed, kGreen+2, kMagenta, kCyan+2, kOrange+2, kViolet, kGray+2};
    
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
    
    // Process each file
    for (size_t file_idx = 0; file_idx < filenames.size(); ++file_idx) {
        const std::string& fname = filenames[file_idx];
        const std::string& config_label = config_labels[file_idx];
        
        std::cout << "Processing file " << (file_idx + 1) << "/" << filenames.size() 
                  << ": " << config_label << std::endl;
        
        // Create 2D histogram for this configuration
        std::string h2_name = "h2_" + std::to_string(file_idx);
        std::string h2_title = "Resolution vs True Kaon Momentum (" + config_label + ")";
        h2_title += ";True Kaon Momentum [GeV/c];|Reco - True|/True";
        
        TH2D *h2 = new TH2D(h2_name.c_str(), h2_title.c_str(), nbins, min_p, max_p, nbins, 0, 1);
        
        // Fill histogram with data from this file
        TFile *inFile = TFile::Open(fname.c_str());
        if (!inFile || inFile->IsZombie()) {
            std::cout << "Warning: Cannot open file " << fname << std::endl;
            continue;
        }

        TTree *tree = (TTree*)inFile->Get("kaonVectors");
        if (!tree) {
            std::cout << "Warning: Tree kaonVectors not found in " << fname << std::endl;
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
            
            // Process all events in this entry's vectors
            for (size_t i = 0; i < reco_p->size(); ++i) {
                double res = std::abs(((*reco_p)[i] - (*true_p)[i]) / (*true_p)[i]);
                if (res > anomaly_threshold) continue; // Skip anomalous results
                h2->Fill((*true_p)[i], res);
                total_events++;
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
            double mean = proj->GetEntries() > 0 ? proj->GetMean() : 0;
            h_mean->SetBinContent(i, mean);
            delete proj;
        }
        
        // Style the histogram for outlined bars (no fill)
        h_mean->SetFillStyle(0); // No fill - hollow bars
        h_mean->SetLineColor(colors[file_idx % colors.size()]);
        h_mean->SetLineWidth(3); // Thick outline
        h_mean->SetBarWidth(0.8);
        h_mean->SetBarOffset(0.1);
        
        h_mean_vec.push_back(h_mean);
        
        // Clean up 2D histogram
        delete h2;
    }
    
    if (h_mean_vec.empty()) {
        std::cout << "Error: No valid histograms created!" << std::endl;
        return;
    }
    
    // Draw histograms
    gStyle->SetOptStat(0); // Turn off statistics box
    
    // Set up the first histogram
    h_mean_vec[0]->SetTitle("Kaon Momentum Resolution Comparison");
    h_mean_vec[0]->GetYaxis()->SetRangeUser(0, 0.5); // Adjust range as needed
    h_mean_vec[0]->Draw("BAR");
    
    // Draw additional histograms
    for (size_t i = 1; i < h_mean_vec.size(); ++i) {
        h_mean_vec[i]->Draw("BAR SAME");
    }
    
    // Switch to right pad for legend
    c1->cd(2);
    
    // Function to wrap text into multiple lines after P1 and P2 (pizza parameters)
    auto wrapTextToLines = [](const std::string& text, size_t maxWidth) -> std::vector<std::string> {
        std::vector<std::string> lines;
        std::istringstream words(text);
        std::string word;
        std::string currentLine = "";
        bool foundP1 = false, foundP2 = false;
        
        while (words >> word) {
            if (currentLine.empty()) {
                currentLine = word;
            } else {
                std::string testLine = currentLine + " " + word;
                
                // Track if we've seen P1 and P2 parameters
                if (word == "P1") foundP1 = true;
                if (word == "P2") foundP2 = true;
                
                // Check if we should force a break after both pizza parameters (P1 and P2)
                bool shouldBreakAfterPizzas = false;
                if (foundP1 && foundP2 && (word == "F1" || word.find("F1") == 0)) {
                    shouldBreakAfterPizzas = true;
                }
                
                if (testLine.length() <= maxWidth && !shouldBreakAfterPizzas) {
                    currentLine = testLine;
                } else {
                    lines.push_back(currentLine);
                    currentLine = word;
                }
            }
        }
        if (!currentLine.empty()) {
            lines.push_back(currentLine);
        }
        
        return lines;
    };
    
    // Clear the pad and set up for manual drawing
    gPad->Clear();
    gPad->Range(0, 0, 1, 1);
    
    // Draw simple legend with smaller text positioned to the left
    double y_start = 0.9;
    double y_step = 0.8 / h_mean_vec.size();  // Distribute entries evenly
    
    for (size_t i = 0; i < h_mean_vec.size(); ++i) {
        double entry_y = y_start - i * y_step;
        
        // Draw the color line indicator
        TLine *colorLine = new TLine(0.01, entry_y, 0.08, entry_y);
        colorLine->SetLineColor(h_mean_vec[i]->GetLineColor());
        colorLine->SetLineWidth(3);
        colorLine->Draw();
        
        // Draw text with smaller size and positioned further left
        TText *text = new TText(0.10, entry_y, config_labels[i].c_str());
        text->SetTextSize(0.025);  // Smaller text size
        text->SetTextColor(kBlack);
        text->Draw();
    }
    
    // Switch back to left pad
    c1->cd(1);
    
    // Add grid for better readability
    c1->SetGrid(1, 1);
    c1->Update();
    
    // Save the plot - save entire canvas to include both pads
    std::string output_filename = "KLong_resolution_comparison.png";
    c1->SaveAs(output_filename.c_str());
    
    std::cout << "\nComparison plot saved as: " << output_filename << std::endl;
    
    // Print summary statistics
    std::cout << "\nSummary Statistics:" << std::endl;
    std::cout << "Configuration\t\tMean Resolution (0.5-2.0 GeV)\tMean Resolution (2.0-5.0 GeV)" << std::endl;
    std::cout << "-------------\t\t--------------------------\t--------------------------" << std::endl;
    
    for (size_t i = 0; i < h_mean_vec.size(); ++i) {
        TH1D* h = h_mean_vec[i];
        
        // Calculate mean resolution in different momentum ranges
        double sum_low = 0, count_low = 0;
        double sum_high = 0, count_high = 0;
        
        for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
            double momentum = h->GetBinCenter(bin);
            double resolution = h->GetBinContent(bin);
            
            if (resolution > 0) {
                if (momentum >= 0.5 && momentum <= 2.0) {
                    sum_low += resolution;
                    count_low++;
                } else if (momentum >= 2.0 && momentum <= 5.0) {
                    sum_high += resolution;
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
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>

void quick_test(const char* filename) {
    TFile* file = new TFile(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cout << "ERROR: Cannot open file " << filename << std::endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("KLongTree");
    if (!tree) {
        std::cout << "ERROR: Cannot find KLongTree in file" << std::endl;
        return;
    }
    
    // Get basic info
    Long64_t nentries = tree->GetEntries();
    std::cout << "Total entries in tree: " << nentries << std::endl;
    
    // Set up branches for hit counting
    std::vector<double> *trkrX = 0, *trkrY = 0, *trkrZ = 0;
    tree->SetBranchAddress("TrkrHitX", &trkrX);
    tree->SetBranchAddress("TrkrHitY", &trkrY);  
    tree->SetBranchAddress("TrkrHitZ", &trkrZ);
    
    int events_with_hits = 0;
    int total_hits = 0;
    
    // Quick analysis of first 1000 events
    int check_events = TMath::Min((Long64_t)1000, nentries);
    
    for (int i = 0; i < check_events; i++) {
        tree->GetEntry(i);
        
        if (trkrX && trkrX->size() > 0) {
            events_with_hits++;
            total_hits += trkrX->size();
        }
    }
    
    std::cout << "Events checked: " << check_events << std::endl;
    std::cout << "Events with detector hits: " << events_with_hits << std::endl;
    std::cout << "Total detector hits: " << total_hits << std::endl;
    std::cout << "Average hits per event with hits: " << (events_with_hits > 0 ? (double)total_hits/events_with_hits : 0) << std::endl;
    std::cout << "Hit rate: " << (double)events_with_hits/check_events * 100 << "%" << std::endl;
    
    file->Close();
}
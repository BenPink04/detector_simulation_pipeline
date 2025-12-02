#include "TFile.h"
#include "TTree.h"
#include <iostream>

void check_file_contents(const char* filename) {
    TFile* file = new TFile(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cout << "ERROR: Cannot open file " << filename << std::endl;
        return;
    }
    
    std::cout << "File: " << filename << std::endl;
    std::cout << "File contents:" << std::endl;
    file->ls();
    
    // Look for any TTree objects
    TList* keys = file->GetListOfKeys();
    for (int i = 0; i < keys->GetSize(); i++) {
        TObject* obj = file->Get(keys->At(i)->GetName());
        if (obj && obj->IsA() == TTree::Class()) {
            TTree* tree = (TTree*)obj;
            std::cout << "Found tree: " << tree->GetName() << " with " << tree->GetEntries() << " entries" << std::endl;
            
            // List branches
            TObjArray* branches = tree->GetListOfBranches();
            std::cout << "  Branches:" << std::endl;
            for (int j = 0; j < branches->GetEntries(); j++) {
                std::cout << "    " << branches->At(j)->GetName() << std::endl;
            }
        }
    }
    
    file->Close();
}
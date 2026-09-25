#define G4MULTITHREADED
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4PhysListFactory.hh"
#include "JLabKDetectorConstruction.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "JLabKActionInitialization.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"

#include <cstdio>

int main(int argc,char** argv)
{
    G4UIExecutive* ui = nullptr;
    bool useVis = false;

    // Interactive mode if no macro file is given
    if (argc == 1) {
        ui = new G4UIExecutive(argc, argv);
        useVis = true;
    } else {
        // If macro file is vis.mac, enable visualization
        std::string macroFile = argv[1];
        if (macroFile.find("vis.mac") != std::string::npos) {
            useVis = true;
        }
    }

    // Random engine setup
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    int seed = (argc > 2) ? atoi(argv[2]) : 8752391;
    CLHEP::HepRandom::setTheSeed(seed);
    G4cout << "Random engine seeded with " << seed << G4endl;

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    runManager->SetUserInitialization(new JLabKDetectorConstruction);

    G4int verbose;
    G4PhysListFactory factory;
    G4VModularPhysicsList* physList = factory.GetReferencePhysList("FTFP_BERT");
    physList->SetVerboseLevel(verbose = 1);
    physList->ReplacePhysics(new G4EmLivermorePolarizedPhysics);
    physList->ReplacePhysics(new G4DecayPhysics);
    physList->ReplacePhysics(new G4RadioactiveDecayPhysics);
    runManager->SetUserInitialization(physList);

    runManager->SetUserInitialization(new JLabKActionInitialization);

    // Only create visualization manager if needed
    G4VisManager* visManager = nullptr;
    if (useVis) {
        visManager = new G4VisExecutive();
        visManager->Initialize();
    }

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (ui == nullptr) {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    } else {
        // Interactive mode
        ui->SessionStart();
        delete ui;
    }

    if (visManager) delete visManager;
    delete runManager;
    return 0;
}

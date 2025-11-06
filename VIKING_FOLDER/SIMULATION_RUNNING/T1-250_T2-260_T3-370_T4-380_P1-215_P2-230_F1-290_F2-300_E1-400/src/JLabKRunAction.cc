#include "JLabKRunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"

JLabKRunAction::JLabKRunAction() : G4UserRunAction()
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  
  G4cout << "Using " << mgr->GetType() << " analysis manager" << G4endl;

  mgr->SetVerboseLevel(0);
  mgr->SetNtupleMerging(true);

  mgr->SetFirstNtupleId(1);

  mgr->CreateNtuple("Ntuple1", "JLabK");
  mgr->CreateNtupleDColumn("evtNb");
  mgr->CreateNtupleDColumn("kaonDK_momX");
  mgr->CreateNtupleDColumn("kaonDK_momY");
  mgr->CreateNtupleDColumn("kaonDK_momZ");
  mgr->CreateNtupleDColumn("kaonDK_posX");
  mgr->CreateNtupleDColumn("kaonDK_posY");
  mgr->CreateNtupleDColumn("kaonDK_posZ");
  mgr->CreateNtupleDColumn("kaonDK_time");
  mgr->FinishNtuple();

  mgr->CreateNtuple("Ntuple2", "JLabK");
  mgr->CreateNtupleDColumn("evtNb");  
  mgr->CreateNtupleDColumn("DKparticle_PDGEncoding");
  mgr->CreateNtupleDColumn("DKparticle_momX");
  mgr->CreateNtupleDColumn("DKparticle_momY");
  mgr->CreateNtupleDColumn("DKparticle_momZ");
  mgr->FinishNtuple();
  
  mgr->CreateNtuple("Ntuple3", "JLabK");
  mgr->CreateNtupleDColumn("evtNb");
  mgr->CreateNtupleDColumn("Edep");
  mgr->CreateNtupleDColumn("hitX");
  mgr->CreateNtupleDColumn("hitY");
  mgr->CreateNtupleDColumn("hitZ");
  mgr->CreateNtupleDColumn("hitT"); 
  mgr->CreateNtupleDColumn("deviceID");
  mgr->CreateNtupleDColumn("PDGEncoding");
  mgr->FinishNtuple();
}

JLabKRunAction::~JLabKRunAction()
{
  // G4AnalysisManager is a singleton managed by Geant4 - do not delete manually
}

void JLabKRunAction::BeginOfRunAction(const G4Run*)
{
  // Get analysis manager
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();

  // Open an output file
  G4String fileName = "T1-250_T2-260_T3-370_T4-380_P1-215_P2-230_F1-290_F2-300_E1-400_100.root";  
  mgr->OpenFile(fileName);
}

void JLabKRunAction::EndOfRunAction(const G4Run* run)
{
  // print histogram statistics
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();
  mgr->Write();
  mgr->CloseFile();
  
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 1) return;

  G4String runType;
  
  if (G4Threading::IsMasterThread()) runType = "Global Run";
  else runType = "Local Run-";

  G4cout
  << "\n----------------------End of " << runType << "------------------------"
  << "\n The run consists of " << nofEvents << " events."
  << "\n------------------------------------------------------------"
  << G4endl;
}

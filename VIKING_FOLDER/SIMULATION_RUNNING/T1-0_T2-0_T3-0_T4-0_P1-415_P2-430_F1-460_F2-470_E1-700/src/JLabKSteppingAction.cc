#include "JLabKSteppingAction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

JLabKSteppingAction::JLabKSteppingAction()
: fKaonDK(false)
, fKaonDK_mom()
, fKaonDK_pos()
, fDKparticleTrackID(0)
, fV_DKparticle_mom()
, fV_DKparticle_PDGEncoding()
, fV_Edep() 
, fV_hitPos()
, fV_deviceID()
, fV_PDGEncoding()
, fV_hitT()
, fKaonDK_time(-1) // Initialize decay time
{}

void JLabKSteppingAction::BeginOfEventAction()
{
  fKaonDK = false;
  fKaonDK_mom = G4ThreeVector(0., 0., 0.);
  fKaonDK_pos = G4ThreeVector(0., 0., 0.);

  fDKparticleTrackID = -1;

  // Removes all elements from the vectors at the beginning of the event
  fV_DKparticle_mom.clear();
  fV_DKparticle_PDGEncoding.clear();
  fV_Edep.clear();
  fV_hitPos.clear();
  fV_deviceID.clear();
  fV_PDGEncoding.clear();
  fV_hitT.clear();
}

void JLabKSteppingAction::EndOfEventAction()
{
  // Record the vector elements in the .root file at the end of the event 
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();

  G4int evtNb =
     G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  // Check if the kaon decayed during the event
  if (fKaonDK)
  {
    mgr->FillNtupleDColumn(1, 0, evtNb);
    mgr->FillNtupleDColumn(1, 1, fKaonDK_mom.x() / GeV);
    mgr->FillNtupleDColumn(1, 2, fKaonDK_mom.y() / GeV);
    mgr->FillNtupleDColumn(1, 3, fKaonDK_mom.z() / GeV);
    mgr->FillNtupleDColumn(1, 4, fKaonDK_pos.x());
    mgr->FillNtupleDColumn(1, 5, fKaonDK_pos.y());
    mgr->FillNtupleDColumn(1, 6, fKaonDK_pos.z());
    mgr->FillNtupleDColumn(1, 7, fKaonDK_time / ns);
    mgr->AddNtupleRow(1);

    G4int nbOfDKparticles = fV_DKparticle_mom.size();

    // Loop over the particle produced by the kaon decay
    for (G4int i = 0; i < nbOfDKparticles; i++)
    {
      G4ThreeVector DKparticle_mom = fV_DKparticle_mom.at(i);

      mgr->FillNtupleDColumn(2, 0, evtNb);
      mgr->FillNtupleDColumn(2, 1, fV_DKparticle_PDGEncoding.at(i));
      mgr->FillNtupleDColumn(2, 2, DKparticle_mom.x());
      mgr->FillNtupleDColumn(2, 3, DKparticle_mom.y());
      mgr->FillNtupleDColumn(2, 4, DKparticle_mom.z());
      mgr->AddNtupleRow(2);
    }
  }

  // Loop over the number of hits outside the guide tube during the events. This
  // is equal to the number of elements in the vector fV_Edep()
  G4int nbOfHits = fV_Edep.size();

  for (G4int i = 0; i < nbOfHits; i++)
  {
    G4ThreeVector hitPos = fV_hitPos.at(i);

    mgr->FillNtupleDColumn(3, 0, evtNb);
    mgr->FillNtupleDColumn(3, 1, fV_Edep.at(i));
    mgr->FillNtupleDColumn(3, 2, hitPos.x());
    mgr->FillNtupleDColumn(3, 3, hitPos.y());
    mgr->FillNtupleDColumn(3, 4, hitPos.z());
    mgr->FillNtupleDColumn(3, 5, fV_hitT.at(i)); // Assuming column 5 is for hitT
    mgr->FillNtupleDColumn(3, 6, fV_deviceID.at(i));
    mgr->FillNtupleDColumn(3, 7, fV_PDGEncoding.at(i));
    mgr->AddNtupleRow(3);
  }
}

void JLabKSteppingAction::UserSteppingAction(const G4Step* step)
{ 
  G4Track* track = step->GetTrack();

  G4int PDGEncoding = track->GetParticleDefinition()->GetPDGEncoding();

  // Store the momentum of emission of the particles produced by the kaon decay 
  // and check if this has not been already recorded
  if (track->GetParentID() == 1
      && track->GetCreatorProcess()->GetProcessName() == "Decay"
      && track->GetTrackID() != fDKparticleTrackID)
  {
    fDKparticleTrackID = track->GetTrackID();

    fV_DKparticle_mom.push_back(track->GetVertexMomentumDirection());
    fV_DKparticle_PDGEncoding.push_back(PDGEncoding);
  }

  G4double Edep = step->GetTotalEnergyDeposit();

  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  const G4VProcess* postProcessDefinedStep =
   postStepPoint->GetProcessDefinedStep();
  if (postProcessDefinedStep == nullptr) return;

  if (postStepPoint->GetMaterial() == nullptr) return;

  // Store the kaon momentum and position right when it decays
  if (track->GetParentID() == 0
      && postProcessDefinedStep->GetProcessName() == "Decay")
  {
    fKaonDK = true;
    fKaonDK_mom = step->GetPreStepPoint()->GetMomentum();
    fKaonDK_pos = postStepPoint->GetPosition() / cm;
    fKaonDK_time = postStepPoint->GetGlobalTime() / ns; // Store decay time in ns
  }

  if (postStepPoint->GetPhysicalVolume() == nullptr) return;
  G4String device = postStepPoint->GetPhysicalVolume()->GetName();

  // Explicitly filter out mother volumes
  if (device == "Tracker" || device == "Pizza1" || device == "Pizza2" ||
      device == "FriWall" || device == "TofWall" || device == "ExperimentalHall") {
    return;
  }

  // Only record device ID for specific elements
  if (device == "StrawTube" || device == "FriStrip" || device == "FriStrip1" ||
      device == "FriStrip2" || device == "FriStrip3" || device == "FriStrip4" ||
      device == "FriStripLong" || device == "FriStripInner" ||
      device == "FriStripInner1_cut" || device == "FriStripInner2_cut" ||
      device == "FriStripInner3_cut" || device == "FriStripInner4_cut" ||
      device == "PizzaSlice" || device == "TofBar") 
  {
      G4TouchableHandle theTouchable = postStepPoint->GetTouchableHandle();
      G4int depth = 0;
      if (device == "PizzaSlice") depth = 0; // was previously depth = 1 - changed as cannot understand reasoning
      G4int deviceID = theTouchable->GetCopyNumber(depth);
      fV_deviceID.push_back(deviceID);

      // Store the energy, position, and PDG encoding as before
      fV_Edep.push_back(Edep);
      fV_hitPos.push_back(postStepPoint->GetPosition() / cm);
      fV_PDGEncoding.push_back(PDGEncoding);
      fV_hitT.push_back(postStepPoint->GetGlobalTime() / ns); // Store time in ns
  }

  return;
}

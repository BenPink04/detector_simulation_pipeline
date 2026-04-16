// ============================================================================
// JLabKSteppingAction.cc
//
// PURPOSE:
//   Geant4 stepping action that runs once per simulation STEP.
//   It is responsible for extracting and buffering physics information
//   into per-event vectors, which are then written to ROOT ntuples at
//   the end of each event.
//
// THREE DATA STREAMS (corresponding to Ntuple1, Ntuple2, Ntuple3 in ROOT):
//
//   Ntuple1 (index 1) — Kaon decay info (one row per event where K decays):
//     evtNb, kaonDK_momX/Y/Z (GeV/c), kaonDK_posX/Y/Z (cm), kaonDK_time (ns)
//
//   Ntuple2 (index 2) — Decay products (one row per product per event):
//     evtNb, DKparticle_PDGEncoding, momX/Y/Z (direction, not momentum magnitude)
//
//   Ntuple3 (index 3) — Detector hits (one row per hit per event):
//     evtNb, Edep, hitX/Y/Z (cm), hitT (ns), deviceID, PDGEncoding
//
// HOW HITS GET INTO THE ROOT FILE:
//   1. UserSteppingAction() runs every step and fills per-event C++ vectors.
//   2. EndOfEventAction() drains those vectors into ROOT ntuple rows.
//   3. BeginOfEventAction() clears the vectors for the next event.
// ============================================================================

#include "JLabKSteppingAction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

// ----------------------------------------------------------------------------
// Constructor: initialise all member variables to safe defaults.
// The per-event vectors are empty by default; they are filled in
// UserSteppingAction() and cleared in BeginOfEventAction().
// ----------------------------------------------------------------------------
JLabKSteppingAction::JLabKSteppingAction()
: fKaonDK(false)               // flag: did the primary kaon decay this event?
, fKaonDK_mom()                // kaon 3-momentum at decay point (GeV/c)
, fKaonDK_pos()                // kaon decay position (cm)
, fDKparticleTrackID(0)        // track ID of the last recorded decay product
, fV_DKparticle_mom()          // vector of decay-product momentum DIRECTIONS
, fV_DKparticle_PDGEncoding()  // vector of decay-product PDG codes
, fV_Edep()                    // vector of energy deposits per hit (MeV)
, fV_hitPos()                  // vector of hit positions (cm)
, fV_deviceID()                // vector of detector element copy numbers
, fV_PDGEncoding()             // vector of PDG codes for each hit particle
, fV_hitT()                    // vector of global hit times (ns)
, fKaonDK_time(-1)             // global time of kaon decay (ns); -1 = not decayed
{}

// ----------------------------------------------------------------------------
// BeginOfEventAction — called once at the START of each event.
//
// Resets all per-event accumulators so that data from the previous event
// does not leak into the current one.
// [TO CHANGE] If you add new per-event member vectors, clear them here too.
// ----------------------------------------------------------------------------
void JLabKSteppingAction::BeginOfEventAction()
{
  fKaonDK = false;                              // reset decay flag
  fKaonDK_mom = G4ThreeVector(0., 0., 0.);      // reset kaon decay momentum
  fKaonDK_pos = G4ThreeVector(0., 0., 0.);      // reset kaon decay position

  fDKparticleTrackID = -1; // -1 means no decay product recorded yet

  // Clear all per-event hit buffers
  fV_DKparticle_mom.clear();
  fV_DKparticle_PDGEncoding.clear();
  fV_Edep.clear();
  fV_hitPos.clear();
  fV_deviceID.clear();
  fV_PDGEncoding.clear();
  fV_hitT.clear();
}

// ----------------------------------------------------------------------------
// EndOfEventAction — called once at the END of each event.
//
// This is where all accumulated step data is written to the ROOT file.
// The G4AnalysisManager ntuples must have been defined elsewhere (typically
// in the RunAction) with the column layout assumed below.
//
// NTUPLE COLUMN LAYOUT:
//
//   Ntuple 1 — "Kaon decay" (one row per event where the kaon decayed)
//     Col 0: evtNb
//     Col 1: kaonDK_momX  (GeV/c)
//     Col 2: kaonDK_momY  (GeV/c)
//     Col 3: kaonDK_momZ  (GeV/c)
//     Col 4: kaonDK_posX  (cm, already divided in UserSteppingAction)
//     Col 5: kaonDK_posY  (cm)
//     Col 6: kaonDK_posZ  (cm)
//     Col 7: kaonDK_time  (ns)
//
//   Ntuple 2 — "Decay products" (one row per product per event)
//     Col 0: evtNb
//     Col 1: DKparticle_PDGEncoding
//     Col 2: DKparticle_momX  (unit direction, NOT magnitude)
//     Col 3: DKparticle_momY
//     Col 4: DKparticle_momZ
//
//   Ntuple 3 — "Detector hits" (one row per qualifying hit per event)
//     Col 0: evtNb
//     Col 1: Edep       (MeV)
//     Col 2: hitX       (cm)
//     Col 3: hitY       (cm)
//     Col 4: hitZ       (cm)
//     Col 5: hitT       (ns)
//     Col 6: deviceID   (copy number of the detector element)
//     Col 7: PDGEncoding (particle type that made the hit)
//
// [TO CHANGE] If you add a new branch to any ntuple (in RunAction), you must
//             add a corresponding FillNtupleDColumn() call here with the
//             correct ntuple index and column index.
// ----------------------------------------------------------------------------
void JLabKSteppingAction::EndOfEventAction()
{
  G4AnalysisManager* mgr = G4AnalysisManager::Instance();

  // Get the current event number (used as a key in all ntuples)
  G4int evtNb =
     G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  // ---- Ntuple 1 + Ntuple 2: written only if the kaon decayed this event ----
  if (fKaonDK)
  {
    // --- Ntuple 1: one row summarising the kaon decay ---
    mgr->FillNtupleDColumn(1, 0, evtNb);
    mgr->FillNtupleDColumn(1, 1, fKaonDK_mom.x() / GeV); // convert Geant4 internal units -> GeV
    mgr->FillNtupleDColumn(1, 2, fKaonDK_mom.y() / GeV);
    mgr->FillNtupleDColumn(1, 3, fKaonDK_mom.z() / GeV);
    mgr->FillNtupleDColumn(1, 4, fKaonDK_pos.x()); // already in cm (divided in UserSteppingAction)
    mgr->FillNtupleDColumn(1, 5, fKaonDK_pos.y());
    mgr->FillNtupleDColumn(1, 6, fKaonDK_pos.z());
    mgr->FillNtupleDColumn(1, 7, fKaonDK_time / ns); // convert to ns
    mgr->AddNtupleRow(1); // commit this row to the ROOT file

    // --- Ntuple 2: one row per decay product ---
    G4int nbOfDKparticles = fV_DKparticle_mom.size();
    for (G4int i = 0; i < nbOfDKparticles; i++)
    {
      G4ThreeVector DKparticle_mom = fV_DKparticle_mom.at(i); // unit direction vector

      mgr->FillNtupleDColumn(2, 0, evtNb);
      mgr->FillNtupleDColumn(2, 1, fV_DKparticle_PDGEncoding.at(i));
      mgr->FillNtupleDColumn(2, 2, DKparticle_mom.x()); // NOTE: direction only, not |p|
      mgr->FillNtupleDColumn(2, 3, DKparticle_mom.y()); // [TO CHANGE] if you need |p|,
      mgr->FillNtupleDColumn(2, 4, DKparticle_mom.z()); //   store GetMomentum() instead
      mgr->AddNtupleRow(2);
    }
  }

  // ---- Ntuple 3: one row per detector hit (all events) --------------------
  // fV_Edep is the master vector; all others are filled in lockstep with it
  // so they always have the same length.
  G4int nbOfHits = fV_Edep.size();

  for (G4int i = 0; i < nbOfHits; i++)
  {
    G4ThreeVector hitPos = fV_hitPos.at(i);

    mgr->FillNtupleDColumn(3, 0, evtNb);
    mgr->FillNtupleDColumn(3, 1, fV_Edep.at(i));      // energy deposit (MeV)
    mgr->FillNtupleDColumn(3, 2, hitPos.x());          // hit position x (cm)
    mgr->FillNtupleDColumn(3, 3, hitPos.y());
    mgr->FillNtupleDColumn(3, 4, hitPos.z());
    mgr->FillNtupleDColumn(3, 5, fV_hitT.at(i));      // global hit time (ns)
    mgr->FillNtupleDColumn(3, 6, fV_deviceID.at(i));  // detector element ID (copy number)
    mgr->FillNtupleDColumn(3, 7, fV_PDGEncoding.at(i)); // particle type
    mgr->AddNtupleRow(3);
  }
}

// ----------------------------------------------------------------------------
// UserSteppingAction — called by Geant4 once for EVERY step of EVERY particle.
//
// A "step" is a particle moving from one point (preStepPoint) to the next
// (postStepPoint).  This function inspects each step and decides whether to
// buffer its data into the per-event vectors.
//
// THREE things are tracked here:
//   A) Decay products of the primary kaon  -> fV_DKparticle_*
//   B) Kaon decay point and momentum       -> fKaonDK_*, fKaonDK_time
//   C) Detector hits in sensitive volumes  -> fV_Edep, fV_hitPos, etc.
// ----------------------------------------------------------------------------
void JLabKSteppingAction::UserSteppingAction(const G4Step* step)
{ 
  G4Track* track = step->GetTrack();

  // PDG code of the particle making this step (e.g. 211 = pi+, -211 = pi-, 311 = K0)
  G4int PDGEncoding = track->GetParticleDefinition()->GetPDGEncoding();

  // ==========================================================================
  // A) RECORD KAON DECAY PRODUCTS (Ntuple 2)
  //
  // Condition: particle was created by a Decay process AND its parent is
  // track 1 (the primary kaon).  We also check fDKparticleTrackID to avoid
  // recording the same track multiple times (a track can have many steps).
  //
  // We store the VERTEX MOMENTUM DIRECTION (unit vector at birth) and the PDG
  // code.  The actual momentum magnitude is NOT stored here.
  //
  // [TO CHANGE] If you need the momentum magnitude of each decay product,
  //   replace GetVertexMomentumDirection() with GetVertexMomentum() and store
  //   the full 3-vector (or its magnitude) in a separate member variable.
  // [TO CHANGE] If the kaon is not always track ID 1 (e.g. in overlay events),
  //   you may need a more robust parent-matching strategy.
  // ==========================================================================
  if (track->GetParentID() == 1                                          // parent is primary kaon
      && track->GetCreatorProcess()->GetProcessName() == "Decay"         // created by decay
      && track->GetTrackID() != fDKparticleTrackID)                      // not already recorded
  {
    fDKparticleTrackID = track->GetTrackID(); // remember this track so we don't double-count

    fV_DKparticle_mom.push_back(track->GetVertexMomentumDirection()); // unit direction at birth
    fV_DKparticle_PDGEncoding.push_back(PDGEncoding);
  }

  // Energy deposited in this step (MeV in Geant4 internal units)
  G4double Edep = step->GetTotalEnergyDeposit();

  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  // Safety checks: skip steps with no defined process or material
  // (can happen at world boundaries or for optical photons)
  const G4VProcess* postProcessDefinedStep =
   postStepPoint->GetProcessDefinedStep();
  if (postProcessDefinedStep == nullptr) return;

  if (postStepPoint->GetMaterial() == nullptr) return;

  // ==========================================================================
  // B) RECORD KAON DECAY POINT (Ntuple 1)
  //
  // Condition: the PRIMARY particle (parentID == 0, i.e. the kaon) undergoes
  // a Decay process at the POST step point.
  // We capture the pre-step momentum (the kaon's momentum just before decay)
  // and the post-step position (where the decay happened).
  //
  // fKaonDK_time is stored in ns and read back in KLong_save_vectors.C as the
  // kaon decay time for the TOF-based momentum reconstruction.
  //
  // [TO CHANGE] If you want the kaon's momentum AT the decay point (vs just
  //   before), use postStepPoint->GetMomentum() instead of preStepPoint.
  //   For most cases (straight tracks) the difference is negligible.
  // ==========================================================================
  if (track->GetParentID() == 0                                          // primary particle
      && postProcessDefinedStep->GetProcessName() == "Decay")            // decayed this step
  {
    fKaonDK = true;
    fKaonDK_mom  = step->GetPreStepPoint()->GetMomentum();               // 3-momentum before decay (Geant4 units)
    fKaonDK_pos  = postStepPoint->GetPosition() / cm;                    // decay position in cm
    fKaonDK_time = postStepPoint->GetGlobalTime() / ns;                  // global time in ns
  }

  // Can't record a hit if there is no physical volume at the post-step point
  if (postStepPoint->GetPhysicalVolume() == nullptr) return;

  // Name of the physical volume the particle has just entered
  G4String device = postStepPoint->GetPhysicalVolume()->GetName();

  // ==========================================================================
  // C) RECORD DETECTOR HITS (Ntuple 3)
  //
  // We only record hits in specific SENSITIVE volumes (individual detector
  // elements).  Mother/envelope volumes are explicitly rejected first to avoid
  // double-counting when a step crosses a volume boundary.
  //
  // HOW DEVICE IDs ARE ASSIGNED:
  //   Each sensitive element has a unique copy number in the Geant4 geometry
  //   (set by G4PVPlacement's copyNo argument in the DetectorConstruction).
  //   GetCopyNumber(depth) retrieves the copy number at geometry depth `depth`.
  //   depth=0 means the volume itself; higher depths are parent volumes.
  //
  //   The globally unique deviceID ranges used in KLong_save_vectors.C are:
  //     StrawTube (tracker)  Station 0 (T1/X):  IDs 1-488
  //                          Station 1 (T2/Y):  IDs 489-976
  //                          Station 2 (T3/U):  IDs 977-1464
  //                          Station 3 (T4/V):  IDs 1465-1952
  //     PizzaSlice           PIZZA flat:         IDs 1953-1976
  //                          PIZZA conical:      IDs 1977-2000
  //     FriStrip*            FRI Wall 1 (X):     IDs 2001-2026
  //                          FRI Wall 2 (Y):     IDs 2027-2052
  //     TofBar               TOF wall:           IDs 2053-2070
  //
  // [TO CHANGE] If you add a new sensitive detector type, add its volume name
  //   to the allowlist below AND add the corresponding ID range to the
  //   is_* lambdas in KLong_save_vectors.C so hits are classified correctly.
  // [TO CHANGE] If the geometry changes and copy numbers shift, update the
  //   ID ranges in both this file's comments and in KLong_save_vectors.C.
  // ==========================================================================

  // --- Reject envelope/mother volumes (no sensitive material of their own) ---
  // [TO CHANGE] Add any new envelope volume names here to prevent double-counting.
  if (device == "Tracker" || device == "Pizza1" || device == "Pizza2" ||
      device == "FriWall" || device == "TofWall" || device == "ExperimentalHall") {
    return;
  }

  // --- Accept only known sensitive detector element volumes ---
  // [TO CHANGE] Add new sensitive volume names here if you extend the detector.
  if (device == "StrawTube" || device == "FriStrip" || device == "FriStrip1" ||
      device == "FriStrip2" || device == "FriStrip3" || device == "FriStrip4" ||
      device == "FriStripLong" || device == "FriStripInner" ||
      device == "FriStripInner1_cut" || device == "FriStripInner2_cut" ||
      device == "FriStripInner3_cut" || device == "FriStripInner4_cut" ||
      device == "PizzaSlice" || device == "TofBar") 
  {
      // Retrieve the copy number to identify WHICH detector element was hit.
      // depth=0 reads the copy number of the volume itself (the leaf node).
      // [TO CHANGE] If a volume's copy number is defined on a parent level
      //   (e.g. a slice inside a ring), increase depth accordingly.
      G4TouchableHandle theTouchable = postStepPoint->GetTouchableHandle();
      G4int deviceID;
      if (device == "StrawTube") {
          // Encode station identity into the deviceID so each physical straw has a
          // globally unique ID across all 4 tracker placements.
          //   GetCopyNumber(0) = straw copy number within station box (1-488)
          //   GetCopyNumber(1) = parent tracker box copy number (620-623)
          // Resulting deviceID ranges by station:
          //   Station 0 (T1/X): 1-488      Station 1 (T2/Y): 489-976
          //   Station 2 (T3/U): 977-1464   Station 3 (T4/V): 1465-1952
          G4int strawCN = theTouchable->GetCopyNumber(0);  // 1-488 within box
          G4int boxCN   = theTouchable->GetCopyNumber(1);  // 620-623
          G4int station = boxCN - 620;                     // 0, 1, 2, or 3
          deviceID = station * 488 + strawCN;              // 1-1952, globally unique
      } else {
          // PIZZA/FRI/TOF: original copy numbers are 489-606.
          // Shift up by 1464 to clear tracker ID space 1-1952:
          //   PIZZA flat    (was 489-512) -> 1953-1976
          //   PIZZA conical (was 513-536) -> 1977-2000
          //   FRI Wall 1    (was 537-562) -> 2001-2026
          //   FRI Wall 2    (was 563-588) -> 2027-2052
          //   TOF           (was 589-606) -> 2053-2070
          deviceID = theTouchable->GetCopyNumber(0) + 1464;
      }

      // Push all hit quantities into the per-event vectors.
      // All vectors are always pushed together so they stay the same length.
      fV_deviceID.push_back(deviceID);
      fV_Edep.push_back(Edep);                                      // MeV
      fV_hitPos.push_back(postStepPoint->GetPosition() / cm);       // cm
      fV_PDGEncoding.push_back(PDGEncoding);                        // particle type
      fV_hitT.push_back(postStepPoint->GetGlobalTime() / ns);       // ns
  }

  return;
}

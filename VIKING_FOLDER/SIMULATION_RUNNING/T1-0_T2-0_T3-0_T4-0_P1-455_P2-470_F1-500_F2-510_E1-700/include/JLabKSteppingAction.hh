#ifndef AnnihilationPhotonsSteppingAction_hh
#define AnnihilationPhotonsSteppingAction_hh

#include "JLabKVSteppingAction.hh"

#include "G4ThreeVector.hh"

#include <vector>

class JLabKSteppingAction: public JLabKVSteppingAction
{
public:

  JLabKSteppingAction();
  virtual void BeginOfEventAction();
  virtual void UserSteppingAction(const G4Step*);
  virtual void EndOfEventAction();

private:

  // These are used to remember quantities from call to call of
  // UserSteppingAction
  G4bool fKaonDK;
  G4ThreeVector fKaonDK_mom;
  G4ThreeVector fKaonDK_pos;
  G4int fDKparticleTrackID;
  std::vector<G4ThreeVector> fV_DKparticle_mom;
  std::vector<G4int> fV_DKparticle_PDGEncoding;

  std::vector<G4double> fV_Edep;
  std::vector<G4ThreeVector> fV_hitPos;
  std::vector<G4int> fV_deviceID;
  std::vector<G4int> fV_PDGEncoding;
  std::vector<G4double> fV_hitT;

  G4double fKaonDK_time; // Decay time in ns
};

#endif

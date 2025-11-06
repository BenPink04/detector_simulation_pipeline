#ifndef JLabKPrimaryGeneratorAction_hh
#define JLabKPrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class G4GeneralParticleSource;

class JLabKPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  JLabKPrimaryGeneratorAction();
  virtual ~JLabKPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun* fpParticleGun;
  G4int evtNb;
};

#endif

#include "JLabKPrimaryGeneratorAction.hh"

#include "Randomize.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"

JLabKPrimaryGeneratorAction::JLabKPrimaryGeneratorAction()
{
  fpParticleGun = new G4ParticleGun(1);
}

JLabKPrimaryGeneratorAction::~JLabKPrimaryGeneratorAction()
{
  delete fpParticleGun;
}

void JLabKPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("kaon0L");
  // G4ParticleDefinition* particle = particleTable->FindParticle("e-");

  fpParticleGun->SetParticleDefinition(particle);

  // Sample position uniformly in a 3 mm radius circle (x-y plane)
  G4double r_max = 3.0 * mm;
  G4double r = r_max * std::sqrt(G4UniformRand());
  G4double theta = 2.0 * M_PI * G4UniformRand();
  G4double x0 = r * std::cos(theta);
  G4double y0 = r * std::sin(theta);
  G4double z0 = 0.0;

  fpParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

  // Uniform sampling of the initial energy between 0.1 and 10 GeV 
  G4double Einit = (.1 + (10. - .1) * G4UniformRand()) * GeV;

  fpParticleGun->SetParticleEnergy(Einit);
  fpParticleGun->GeneratePrimaryVertex(anEvent);

  evtNb = anEvent->GetEventID();

  if ((evtNb + 1) % 1000 == 0)
   G4cout << (evtNb + 1) / 1e6 << " million kaons" << G4endl;
}
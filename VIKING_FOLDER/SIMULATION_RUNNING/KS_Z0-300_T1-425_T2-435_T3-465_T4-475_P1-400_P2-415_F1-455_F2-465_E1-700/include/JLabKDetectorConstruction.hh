#ifndef JLabKDetectorConstruction_hh
#define JLabKDetectorConstruction_hh 1

#include "G4VUserDetectorConstruction.hh"

#include "JLabKMagneticField.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class JLabKMagneticField;

class JLabKDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  JLabKDetectorConstruction();
  virtual ~JLabKDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  virtual void ConstructSDandField();

private:
  G4LogicalVolume* crystal_log;

  static G4ThreadLocal JLabKMagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fMagneticFieldMgr;

  std::vector<G4double> fV_hitT;
};

#endif


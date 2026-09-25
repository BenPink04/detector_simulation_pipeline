#ifndef JLabKMagneticField_h
#define JLabKMagneticField_h 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"

using namespace CLHEP;
using namespace std;

class JLabKMagneticField : public G4MagneticField
{  
  public:

    JLabKMagneticField();
    ~JLabKMagneticField();
  
    void  GetFieldValue(const G4double Point[3], G4double *Bfield) const;

  private:

    G4double fDipole_By;
    G4double fSolenoid_Bz;
  
    // Limits of the dipole magnetic field
    G4double dipoleB_minX, dipoleB_maxX, dipoleB_minY, dipoleB_maxY,
             dipoleB_minZ, dipoleB_maxZ;

    // Limits of the solenoid magnetic field
    G4double solenoidB_minZ, solenoidB_maxZ;
};

#endif

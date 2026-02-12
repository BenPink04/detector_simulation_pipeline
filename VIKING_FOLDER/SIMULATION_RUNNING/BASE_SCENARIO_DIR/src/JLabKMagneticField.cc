#include "JLabKMagneticField.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"

JLabKMagneticField::JLabKMagneticField() 
{
  // Set the magnetic fields in the dipole (y component) and the solenoid 
  // (z component)
  fDipole_By = 1.5*tesla;
  fSolenoid_Bz = 1.5*tesla;

  // Define the spatial extension of the magnetic fields (global coordinates)
  // Need to collect the position of the magnets
    
  // Dipole
  G4VPhysicalVolume* dipole_phys =
   G4PhysicalVolumeStore::GetInstance()->GetVolume("Dipole");

  if (dipole_phys)
  {
    // The magnetic field is assumed to be centered in the dipole
    // Position of the dipole center along the z-axis
    G4double dipoleB_posZ = dipole_phys->GetObjectTranslation().z();

    // Dimensions of the dipole magnetic field along the x-, y- and z-axes
    G4double dipoleB_lengthXZ = 6.*cm;
    G4double dipoleB_lengthY = 100.*cm;

    // Limits of the dipole magnetic field along the x-, y- and z-axes
    dipoleB_minX = -dipoleB_lengthXZ / 2.;
    dipoleB_maxX = dipoleB_lengthXZ / 2.;
    dipoleB_minY = -dipoleB_lengthY / 2.;
    dipoleB_maxY = dipoleB_lengthY / 2.;
    dipoleB_minZ = dipoleB_posZ - dipoleB_lengthXZ / 2.;
    dipoleB_maxZ = dipoleB_posZ + dipoleB_lengthXZ / 2.;
  }

  else 
  {
    G4ExceptionDescription msg1;
    msg1 << "Dipole physical volume not found." << G4endl;
    G4Exception("JLabKMagneticField::JLabKMagneticField()", "MyCode0002", 
                JustWarning, msg1);
  }

  // Solenoid
  G4VPhysicalVolume* solenoid_phys =
   G4PhysicalVolumeStore::GetInstance()->GetVolume("Solenoid");

  if (solenoid_phys)
  {
    // Position of the solenoid center along the z-axis
    G4double solenoidB_posZ = solenoid_phys->GetObjectTranslation().z();

    // Dimensions of the solenoid magnetic field along the x-, y- and z-axes
    G4double solenoidB_height = 200.*cm;

    // Limits of the dipole magnetic field along the x-, y- and z-axes
    solenoidB_minZ = solenoidB_posZ - solenoidB_height / 2.;
    solenoidB_maxZ = solenoidB_posZ + solenoidB_height / 2.;
  }

  else 
  {
    G4ExceptionDescription msg2;
    msg2 << "Dipole physical volume not found." << G4endl;
    G4Exception("JLabKMagneticField::JLabKMagneticField()", "MyCode0002", 
                JustWarning, msg2);
  }
}

JLabKMagneticField::~JLabKMagneticField()
{}

void JLabKMagneticField::GetFieldValue(const G4double point[3], 
                                       G4double *Bfield) const
{
  // Reset the x, y and z components of the magnetic field
  Bfield[0] = 0.;
  Bfield[1] = 0.;
  Bfield[2] = 0.;
  
  // Particle global coordinates
  G4double x = point[0];
  G4double y = point[1];
  G4double z = point[2];

  // Checks if the particle is inside the dipole
  if (x >= dipoleB_minX && x <= dipoleB_maxX 
      && y >= dipoleB_minY && y <= dipoleB_maxY
      && z >= dipoleB_minZ && z <= dipoleB_maxZ)
  {   
    Bfield[1] = fDipole_By;
  }

  // Radial distance
  G4double radDist = sqrt(x*x + y*y);

  // Checks if the particle is inside the solenoid
  if (radDist <= 45.*cm && z >= solenoidB_minZ && z <= solenoidB_maxZ)
  {   
    Bfield[2] = fSolenoid_Bz;
  }
}
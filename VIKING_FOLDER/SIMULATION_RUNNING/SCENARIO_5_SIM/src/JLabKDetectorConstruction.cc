#include "JLabKDetectorConstruction.hh"

#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"

G4ThreadLocal JLabKMagneticField*
 JLabKDetectorConstruction::fMagneticField = nullptr;

G4ThreadLocal G4FieldManager*
 JLabKDetectorConstruction::fMagneticFieldMgr = nullptr;

JLabKDetectorConstruction::JLabKDetectorConstruction()
{}

JLabKDetectorConstruction::~JLabKDetectorConstruction()
{
  if (fMagneticField)
  {
    delete fMagneticField;
    fMagneticField = nullptr;
  }
}

G4VPhysicalVolume* JLabKDetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;

  // Materials
  G4NistManager* nist = G4NistManager::Instance();
  //Including air, aluminium, vaccum and plastic materials
  G4Material* plastic = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* aluminium = nist->FindOrBuildMaterial("G4_Al");
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

  G4double z, a, fractionmass, density;
  G4String name, symbol;
  G4int ncomponents;


  // Defining and creating the 50/50 mixture of Argon and Carbon Dioxide for tracker tube
  a = 39.95*g/mole; //argon
  G4Element* elA = new G4Element(name = "Argon", symbol = "Ar", z = 18., a);

  a = 44.01*g/mole; //CO2
  G4Element* elC = new G4Element(name = "Carbon", symbol = "C", z = 6., a);

  density = 1.88*kg/m3; //density of CO2
  G4Material* TubeGas = new G4Material(name = "Tube_Gas", density, ncomponents = 2);
  TubeGas->AddElement(elA, fractionmass = 0.50);
  TubeGas->AddElement(elC, fractionmass = 0.50);

  //------------------------------------------------------------------------
  // Experimental hall (= Geant4 world)
  G4VSolid* expHall = new G4Box("ExperimentalHall", 100.*m, 100.*m, 100.*m);

  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall, air,
                                                     "ExperimentalHall");

  G4VPhysicalVolume* expHall_phys = new G4PVPlacement(G4Transform3D(), 
                                                      expHall_log,
                                                      "ExperimentalHall", 0,
                                                      false, 0);
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------  
  // Guide tube - Parameters defined but geometry NOT placed to avoid overlaps
  // Keep parameters for use by other components (tracker apertures, etc.)
  G4double gdTube_outerRad = 3.*cm;
  G4double gdTubeThk = .1*cm;
  G4double gdTube_inrRad = gdTube_outerRad - gdTubeThk;
  G4double gdTube_halfL = 1000.*cm;
  
  // Position along z-axis (defined but not used)
  G4double gdTubePosZ = 850.*cm;

  // GuideTube geometry creation and placement DISABLED to eliminate overlaps
  /*

  G4VSolid* gdTube = new G4Tubs("GuideTube", gdTube_inrRad, gdTube_outerRad, 
                                gdTube_halfL, 0., 360.*deg);

  G4LogicalVolume* gdTube_log = new G4LogicalVolume(gdTube, aluminium, 
                                                    "GuideTube");
  
  new G4PVPlacement(0, G4ThreeVector(0., 0., gdTubePosZ), gdTube_log, 
                    "GuideTube", expHall_log, false, 0, checkOverlaps);
// Line 135: THIS REFERENCES THE REMOVED GUIDETUBE
G4VSolid* apert = new G4Tubs("Aperture", 0., gdTube_outerRad,    // ← gdTube_outerRad undefined after commenting out GuideTube
                             trkrHalfZ + .1*mm, 0., 360.*deg);
  // Hollow shell of the guide tube where vacuum is created 
  G4VSolid* gdTubeHollow = new G4Tubs("GuideTubeHollow", 0., gdTube_inrRad,
                                      gdTube_halfL, 0., 360.*deg);

  G4LogicalVolume* gdTubeHollow_log = new G4LogicalVolume(gdTubeHollow, vacuum,
                                                          "GuideTubeHollow");

  new G4PVPlacement(0, G4ThreeVector(0., 0., gdTubePosZ), gdTubeHollow_log,
                    "GuideTubeHollow", expHall_log, false, 0, checkOverlaps);
  //------------------------------------------------------------------------
  */
  //------------------------------------------------------------------------  
  // Dipole = box with an extrusion 
  // Dipole's half length, thickness and extrusion half length
  // This should be modified to refprect the dimensions of the Hall D Pair Spectrometer magnet
  G4double dipHalfEdgeL = 100.*cm;
  G4double dipThk = 50.*cm;
  G4double dipExtrHalfXY = 7.62*cm; // 6 inch total width (3 inch half-width) - Can make edit to be 1.5 inch for 3 inch total width


  // Position along z-axis
  // The center of the dipole entrance is placed at the origin of the lab frame
  G4double dipPosZ = dipHalfEdgeL;

  G4VSolid* dipBox = new G4Box("DipoleBox", dipHalfEdgeL, dipHalfEdgeL, 
                                                          dipHalfEdgeL);

  // Volume to be extruded (need to extend the extrusion a bit along the z-axis 
  // so it appears in the visualisation)
  G4VSolid* dipExtr = new G4Box("DipoleExtrusion", dipExtrHalfXY, dipExtrHalfXY,
                                dipHalfEdgeL + .1*cm);

  // Dipole shape obtained by extruding the "DipoleExtrusion" from the 
  // "DipoleBox"  
  G4VSolid* dipole = new G4SubtractionSolid("Dipole", dipBox, dipExtr, 0,
                                            G4ThreeVector(0., 0., 0.));

  G4LogicalVolume* dipole_log = new G4LogicalVolume(dipole, aluminium, 
                                                    "Dipole");

  new G4PVPlacement(0, G4ThreeVector(0., 0., dipPosZ), dipole_log, "Dipole",
                    expHall_log, false, 0, checkOverlaps);
  //------------------------------------------------------------------------  


  //------------------------------------------------------------------------  
  // Trackers = box with a cylindrical extrusion for the guide tube (aperture) 
  // Half-length along x, y and z
  G4double trkrHalfXY = 50.*cm;
  G4double trkrHalfZ = 2.5*cm;

  G4VSolid* trkrBox = new G4Box("TrackerBox",
                                trkrHalfXY, trkrHalfXY, trkrHalfZ);

  G4LogicalVolume* trkrBox_log = new G4LogicalVolume(trkrBox, air, "TrackerBox");
  
  // Tracker aperture
  G4VSolid* apert = new G4Tubs("Aperture", 0., gdTube_outerRad,
                               trkrHalfZ + .1*mm, 0., 360.*deg);

  //Tracker straws - 4 layers of 122 straws, each layer offset by one straw radius
  // Outer radius, thickness, inner radius and half length
  G4double trkStraw_outerRad = 0.40*cm;
  //G4double trkStraw_outerRad = 4*cm;
  G4double trkStrawThk = .0026*cm;
  //G4double trkStrawThk = 2*cm;
  G4double trkStraw_inrRad = trkStraw_outerRad - trkStrawThk;
  G4double trkStraw_halfL = 45*cm;


  G4double strawPosX[122] = {-48.4*cm, -47.6*cm, -46.8*cm, -46.0*cm, -45.2*cm, -44.4*cm, -43.6*cm, -42.8*cm,
                             -42.0*cm, -41.2*cm, -40.4*cm, -39.6*cm, -38.8*cm, -38.0*cm, -37.2*cm, -36.4*cm,
			     -35.6*cm, -34.8*cm, -34.0*cm, -33.2*cm, -32.4*cm, -31.6*cm, -30.8*cm, -30.0*cm,
                             -29.2*cm, -28.4*cm, -27.6*cm, -26.8*cm, -26.0*cm, -25.2*cm, -24.4*cm, -23.6*cm,
			     -22.8*cm, -22.0*cm, -21.2*cm, -20.4*cm, -19.6*cm, -18.8*cm, -18.0*cm,
			     -17.2*cm, -16.4*cm, -15.6*cm, -14.8*cm, -14.0*cm, -13.2*cm, -12.4*cm, -11.6*cm,
			     -10.8*cm, -10.0*cm,  -9.2*cm,  -8.4*cm,  -7.6*cm,  -6.8*cm,  -6.0*cm,  -5.2*cm,
			      -4.4*cm,  -3.6*cm,  -2.8*cm,  -2.0*cm,  -1.2*cm,  -0.4*cm,   0.4*cm,   1.2*cm,
			       2.0*cm,   2.8*cm,   3.6*cm,   4.4*cm,   5.2*cm,   6.0*cm,   6.8*cm,   7.6*cm,
			       8.4*cm,   9.2*cm,  10.0*cm,  10.8*cm,  11.6*cm,  12.4*cm,  13.2*cm,  14.0*cm,
			      14.8*cm,  15.6*cm,  16.4*cm,  17.2*cm,  18.0*cm,  18.8*cm,  19.6*cm,  20.4*cm,
			      21.2*cm,  22.0*cm,  22.8*cm,  23.6*cm,  24.4*cm,  25.2*cm,  26.0*cm,  26.8*cm,
			      27.6*cm,  28.4*cm,  29.2*cm,  30.0*cm,  30.8*cm,  31.6*cm,  32.4*cm,  33.2*cm,
			      34.0*cm,  34.8*cm,  35.6*cm,  36.4*cm,  37.2*cm,  38.0*cm,  38.8*cm,  39.6*cm,
			      40.4*cm,  41.2*cm,  42.0*cm,  42.8*cm,  43.6*cm,  44.4*cm,  45.2*cm,  46.0*cm,
			      46.8*cm,  47.6*cm};

  G4double strawPosZ[4] = {-1.5*cm, -0.5*cm, 0.5*cm, 1.5*cm};
			     
  G4VSolid* trkStraw = new G4Tubs("StrawTube", trkStraw_inrRad, trkStraw_outerRad, 
                                trkStraw_halfL, 0., 360.*deg);

  G4LogicalVolume* trkStraw_log = new G4LogicalVolume(trkStraw, TubeGas, 
                                                    "StrawTube");

  //Place one layer of straws in tracker box
  for (G4int i = 0; i < 122; i++)
  {
    G4RotationMatrix* trkStrawRot = new G4RotationMatrix;
    trkStrawRot->rotateX(90.*deg);


    new G4PVPlacement(trkStrawRot, G4ThreeVector(strawPosX[i], 0., strawPosZ[0]), trkStraw_log,
		      "StrawTube", trkrBox_log, false, i+1, checkOverlaps);

    new G4PVPlacement(trkStrawRot, G4ThreeVector(strawPosX[i]+0.4*cm, 0., strawPosZ[1]), trkStraw_log,
		      "StrawTube", trkrBox_log, false, i+123, checkOverlaps);

    new G4PVPlacement(trkStrawRot, G4ThreeVector(strawPosX[i], 0., strawPosZ[2]), trkStraw_log,
		      "StrawTube", trkrBox_log, false, i+245, checkOverlaps);

    new G4PVPlacement(trkStrawRot, G4ThreeVector(strawPosX[i]+0.4*cm, 0., strawPosZ[3]), trkStraw_log,
		      "StrawTube", trkrBox_log, false, i+367, checkOverlaps);

 }

  
  // Tracker obtained by extruding the "Aperture" from the "TrackerBox"  
  //G4VSolid* trkr = new G4SubtractionSolid("Traker", trkrBox, apert, 0,
  //                                        G4ThreeVector(0., 0., 0.));

  G4LogicalVolume* trkr_log = new G4LogicalVolume(trkrBox, air, "Tracker");

  // Rotations around the z-axis and a positions of the trackers
  // Order: x, y, u, v (xy then uv orientations)
  G4double trkrPhi[4] = {0.*deg, -90.*deg, 45.*deg, -45.*deg};
  //G4double trkrPosZ[4] = {50.*cm, 60.*cm, 130.*cm, 140.*cm}; // OLD - too close to dipole
  G4double trkrPosZ[4] = {240.*cm, 250.*cm, 570.*cm, 580.*cm}; // CORRECTED - proper positions

  for (G4int i = 0; i < 4; i++)
  {
    G4RotationMatrix* trkrRot = new G4RotationMatrix;
    trkrRot->rotateZ(trkrPhi[i]);

    new G4PVPlacement(trkrRot, G4ThreeVector(0., 0., trkrPosZ[i]), trkrBox_log,
                      "Tracker", expHall_log, false, i+620, checkOverlaps); 
   //Setting all major components to have a copyno of i - might change later if causes conjestion to deviceno around 1-4
  }
  //------------------------------------------------------------------------  

  //------------------------------------------------------------------------    
  // Pizzas = mother volume for the pizza slices
  // Pizzas inner radius, outer radius and half thickness
  G4double pizza_outerRad = 37.2*cm;
  G4double pizza_halfThk = .15*cm;
  G4double pizza_halfVolThk = .95*cm;
  //G4double pizza_halfThk = .13*cm; //the 8 elements around the cross are slightly thinner

  //Flat layer
  G4VSolid* pizza1 = new G4Tubs("Pizza1", gdTube_outerRad, pizza_outerRad,
                               pizza_halfThk, 0., 360.*deg);

  G4LogicalVolume* pizza1_log = new G4LogicalVolume(pizza1, air, "Pizza1");


  //Flat layer (changed from conical to eliminate mother volume overlap issues)
  G4VSolid* pizza2 = new G4Tubs("Pizza2", gdTube_outerRad, pizza_outerRad,
                               pizza_halfThk, 0., 360.*deg);

  G4LogicalVolume* pizza2_log = new G4LogicalVolume(pizza2, air, "Pizza2");

  // Pizza slices to be placed inside the pizza
  // G4VSolid* pizzaSlice = new G4Tubs("PizzaSlice", gdTube_outerRad,
  //                                   pizza_outerRad, pizza_halfThk, 0., 57.*deg);
  G4VSolid* pizzaSlice = new G4Tubs("PizzaSlice", gdTube_outerRad,
                                    pizza_outerRad, pizza_halfThk, 0., 14.*deg);//<15 degrees for 24-slice

  G4LogicalVolume* pizzaSlice_log = new G4LogicalVolume(pizzaSlice, plastic, 
                                                        "PizzaSlice");

  for (G4int i = 0; i < 24; i++)   //Modify for number of slices in WASA detector (24)
  {
    // Angle to place the slices inside the pizza
    G4double pizzaSlicePhi = (30. + i * 15.)*deg;
    //15 degrees for 24-slice. Do we need to start at 30?

    G4RotationMatrix* pizzaSliceRot = new G4RotationMatrix;
    pizzaSliceRot->rotateZ(pizzaSlicePhi);
    //pizzaSliceRot->rotateY(-15.*deg);

    new G4PVPlacement(pizzaSliceRot, G4ThreeVector(0., 0., 0.), pizzaSlice_log,
                      "PizzaSlice", pizza1_log, false, 489 + i, checkOverlaps);
  }

  for (G4int i = 0; i < 24; i++)   //Modify for number of slices in WASA detector (24)
    {
      // Angle to place the slices inside the pizza
      G4double pizzaSlicePhi = (30. + i * 15.)*deg;
      //15 degrees for 24-slice. Do we need to start at 30?
      
      G4RotationMatrix* pizzaSliceRot = new G4RotationMatrix;
      pizzaSliceRot->rotateZ(pizzaSlicePhi);
      //pizzaSliceRot->rotateY(20.*deg);
      //Commented out Y rotation to make flat pizza layer
      new G4PVPlacement(pizzaSliceRot, G4ThreeVector(0., 0., 0.), pizzaSlice_log,
			"PizzaSlice", pizza2_log, false, 513 + i, checkOverlaps);
    }


  // Z position of the pizzas and rotation around the z-axis
  // Reduce number of layers
  G4double pizzaPhi[4] = {15.*deg, 0.*deg, 15.*deg, 0.*deg};
  //G4double pizzaPosZ[4] = {15.*cm, 30.*cm, 1100.*cm, 1200.*cm}; // OLD - too close to dipole
  G4double pizzaPosZ[4] = {215.*cm, 230.*cm, 1100.*cm, 1200.*cm}; // CORRECTED - proper positions

  
  //  for (G4int i = 0; i < 2; i++)
  // {

    G4RotationMatrix* rot2 = new G4RotationMatrix;
    rot2->rotateZ(pizzaPhi[0]);

    new G4PVPlacement(rot2, G4ThreeVector(0., 0., pizzaPosZ[1]), pizza1_log,
                      "Pizza1", expHall_log, false, 0, checkOverlaps);

     G4RotationMatrix* rot3 = new G4RotationMatrix;
     rot3->rotateZ(pizzaPhi[1]);

     new G4PVPlacement(rot3, G4ThreeVector(0., 0., pizzaPosZ[0]), pizza2_log,
                       "Pizza2", expHall_log, false, 0, checkOverlaps);

    
    //  }
  //------------------------------------------------------------------------  
  
  //------------------------------------------------------------------------  
  //Foward Range Hodoscope goes here
  // FriWall = mother volume for the FRI wall
  // Create volume as for the tracker; rectangle (with beam hole extruded?)
  // Half-length along x, y and z
  //Two layers? Rotated 90 degrees wrt each other?
  //G4double tofHalfX = 50.*cm;
  G4double FriHalfX = 70.*cm;
  G4double FriHalfY = 72.*cm;
  G4double FriHalfZ = 0.3*cm;

  G4VSolid* friWall = new G4Box("FriWall",
                                FriHalfX, FriHalfY, FriHalfZ);
 
  //Cylinder to be extracted for a beam hole in visualisation
  //Define beam hole cylinder
  G4double beamholeradius = 5.0*cm;  //CAN ADJUST
  G4double beamholehalflength = FriHalfZ + 5.0*cm;    // Larger than FriHalfZ to ensure it cuts right
  G4VSolid* beamhole = new G4Tubs(
        "BeamHole",            //name
        0.,                    //Inner radius - solid so 0
        beamholeradius,        //outer radius
        beamholehalflength,    //half length in z
        0.,                    //starting phi 
        360.* deg              // total phi angle
  );
  
  //Subtract beam hole from FriWall
  G4VSolid* friWallWithHole = new G4SubtractionSolid(
        "FriWallWithHole",      // New solid name
        friWall,                // Original box
        beamhole,               // Cylinder to subtract
        0,                      // Rotation (N/A so 0)
        G4ThreeVector(0., 0., 0.)  // Position - center 
  );

  //Replace original FriWall with new solid with subtracted center
  G4LogicalVolume* friWall1_log = new G4LogicalVolume(
        friWallWithHole,        // now using subtracted solid
        air,                    // material 
        "FriWall1"
  );

  //Duplicate the FRI wall 
  G4LogicalVolume*friWall2_log = new G4LogicalVolume(
        friWallWithHole,        // now using subtracted solid
        air,                    // material 
        "FriWall2"
  );
 
  // FRI strips to be placed inside the FRI wall
  //Wide outer strips (6 cm), shortest length
  G4double friStripHalfX = 3.*cm;
  G4double friStripHalfY = 37.25*cm;
  G4double friStripHalfZ = 0.25*cm;

  G4VSolid* friStrip = new G4Box("FriStrip",
			       friStripHalfX, friStripHalfY, friStripHalfZ);

  G4LogicalVolume* friStrip_log = new G4LogicalVolume(friStrip, plastic, 
                                                        "FriStrip");

  // FRI strips to be placed inside the FRI wall
  //Wide outer strips (6 cm), full length
  G4double friStripLongHalfY = 70.25*cm;

  G4VSolid* friStripLong = new G4Box("FriStripLong",
			       friStripHalfX, friStripLongHalfY, friStripHalfZ);

  G4LogicalVolume* friStripLong_log = new G4LogicalVolume(friStripLong, plastic,
                                                        "FriStripLong");
                                                        
  // FRI strips to be placed inside the FRI wall
  //Wide outer strips (6cm), Second shortest length
  G4double friStrip_1HalfY = 42.75*cm;
  
  G4VSolid* friStrip_1 = new G4Box("FriStrip1",
             friStripHalfX, friStrip_1HalfY, friStripHalfZ);

  G4LogicalVolume* friStrip_1_log = new G4LogicalVolume(friStrip_1, plastic, 
                                                        "FriStrip1");

  // FRI strips to be placed inside the FRI wall
  //Wide outer strips (6cm), Third shortest length
  G4double friStrip_2HalfY = 53.75*cm;

  G4VSolid* friStrip_2 = new G4Box("FriStrip2",
             friStripHalfX, friStrip_2HalfY, friStripHalfZ);

  G4LogicalVolume* friStrip_2_log = new G4LogicalVolume(friStrip_2, plastic, 
                                                        "FriStrip2");

  // FRI strips to be placed inside the FRI wall
  //Wide outer strips (6cm), Fourth shortest length
  G4double friStrip_3HalfY = 59.25*cm;

  G4VSolid* friStrip_3 = new G4Box("FriStrip3",
             friStripHalfX, friStrip_3HalfY, friStripHalfZ);

  G4LogicalVolume* friStrip_3_log = new G4LogicalVolume(friStrip_3, plastic, 
                                                        "FriStrip3");

  // FRI strips to be placed inside the FRI wall
  //Wide outer strips (6cm), Fifth shortest length
  G4double friStrip_4HalfY = 64.75*cm;

  G4VSolid* friStrip_4 = new G4Box("FriStrip4",
             friStripHalfX, friStrip_4HalfY, friStripHalfZ);

  G4LogicalVolume* friStrip_4_log = new G4LogicalVolume(friStrip_4, plastic, 
                                                        "FriStrip4");
  
  // FRI strips to be placed inside the FRI wall
  //Narrow inner strips (3 cm), No Cutout
  G4double friStripInnerHalfX = 1.5*cm;

  G4VSolid* friStripInner = new G4Box("FriStripInner",
			       friStripInnerHalfX, friStripLongHalfY, friStripHalfZ);

  G4LogicalVolume* friStripInner_log = new G4LogicalVolume(friStripInner, plastic, 
                                                        "FriStripInner");
  
  // FRI strips to be placed inside the FRI wall
  //Narrow inner strips (3 cm), Cutout of edge of cylinder

  G4VSolid* friStripInner_1 = new G4Box("FriStripInner1",
             friStripInnerHalfX, friStripLongHalfY, friStripHalfZ);

  G4LogicalVolume* friStripInner_1_log = new G4LogicalVolume(friStripInner_1, plastic, 
                                                        "FriStripInner1");

  // Subtract beam hole from fristripinner1 with offset
  G4VSolid* friStripInner_1_cut = new G4SubtractionSolid(
        "FriStripInner1_cut",      // New solid name
        friStripInner_1,                // Original box
        beamhole,                       // Cylinder to subtract
        0,                              // Rotation (N/A so 0)
        G4ThreeVector( 4.5*cm, 0., 0.)     // Position - offset by width and a half of strip
  );          
  
  G4LogicalVolume* friStripInner_1_cut_log = new G4LogicalVolume(
        friStripInner_1_cut,        // now using subtracted solid
        plastic,                      // material 
        "FriStripInner1_cut"
  );

  // FRI strips to be placed inside the FRI wall
  //Narrow inner strips (3 cm), Cutout of center of cylinder

  G4VSolid* friStripInner_2 = new G4Box("FriStripInner2",
             friStripInnerHalfX, friStripLongHalfY, friStripHalfZ);
  
  G4LogicalVolume* friStripInner_2_log = new G4LogicalVolume(friStripInner_2, plastic, 
                                                        "FriStripInner2");

  // Subtract beam hole from fristripinner2 with offset
  G4VSolid* friStripInner_2_cut = new G4SubtractionSolid(
        "FriStripInner2_cut",      // New solid name
        friStripInner_2,                // Original box
        beamhole,                       // Cylinder to subtract
        0,                              // Rotation (N/A so 0)
        G4ThreeVector( 1.5*cm, 0., 0.)     // Position - offset by half a width of strip
  );

  G4LogicalVolume* friStripInner_2_cut_log = new G4LogicalVolume(
        friStripInner_2_cut,        // now using subtracted solid
        plastic,                      // material 
        "FriStripInner2_cut"
  );

  // FRI strips to be placed inside the FRI wall
  //Narrow inner strips (3 cm), Cutout of center of cylinder other side
  
  G4VSolid* friStripInner_3 = new G4Box("FriStripInner3",
             friStripInnerHalfX, friStripLongHalfY, friStripHalfZ);

  G4LogicalVolume* friStripInner_3_log = new G4LogicalVolume(friStripInner_3, plastic, 
                                                        "FriStripInner3");

  // Subtract beam hole from fristripinner3 with offset
  G4VSolid* friStripInner_3_cut = new G4SubtractionSolid(
        "FriStripInner3_cut",      // New solid name
        friStripInner_3,                // Original box
        beamhole,                       // Cylinder to subtract   
        0,                              // Rotation (N/A so 0)
        G4ThreeVector( -1.5*cm, 0., 0.)     // Position - offset by half a width of strip
  );

  G4LogicalVolume* friStripInner_3_cut_log = new G4LogicalVolume(
        friStripInner_3_cut,        // now using subtracted solid
        plastic,                      // material 
        "FriStripInner3_cut"
  );

  // FRI strips to be placed inside the FRI wall
  //Narrow inner strips (3 cm), Cutout of edge of cylinder other side

  G4VSolid* friStripInner_4 = new G4Box("FriStripInner4",
             friStripInnerHalfX, friStripLongHalfY, friStripHalfZ);

  G4LogicalVolume* friStripInner_4_log = new G4LogicalVolume(friStripInner_4, plastic, 
                                                        "FriStripInner4");

  // Subtract beam hole from fristripinner4 with offset
  G4VSolid* friStripInner_4_cut = new G4SubtractionSolid(
        "FriStripInner4_cut",      // New solid name
        friStripInner_4,                // Original box
        beamhole,                       // Cylinder to subtract
        0,                              // Rotation (N/A so 0)
        G4ThreeVector( -4.5*cm, 0., 0.)     // Position - offset by a width and a half of strip
  );

  G4LogicalVolume* friStripInner_4_cut_log = new G4LogicalVolume(
        friStripInner_4_cut,        // now using subtracted solid
        plastic,                      // material 
        "FriStripInner4_cut"
  );

  //---------------------
  // Place the FRI wall in the experimental hall                        
  //G4double stripPosX[12] = {-55.*cm, -45.*cm, -35.*cm, -25.*cm, -15.*cm, -5.*cm, 5.*cm, 15.*cm, 25.*cm, 35.*cm, 45.*cm, 55.*cm};
  G4double stripPosX[52] = {-66.*cm, -60.*cm, -54.*cm, -48.*cm, -42.*cm, -36.*cm, -30.*cm, -24.*cm, -18.*cm, -12.*cm, -7.5*cm,
    -4.5*cm, -1.5*cm, 1.5*cm, 4.5*cm, 7.5*cm, 12.*cm, 18.*cm, 24.*cm, 30.*cm, 36.*cm, 42.*cm, 48.*cm, 54.*cm, 60.*cm, 66.*cm,
                            -66.*cm, -60.*cm, -54.*cm, -48.*cm, -42.*cm, -36.*cm, -30.*cm, -24.*cm, -18.*cm, -12.*cm, -7.5*cm,
    -4.5*cm, -1.5*cm, 1.5*cm, 4.5*cm, 7.5*cm, 12.*cm, 18.*cm, 24.*cm, 30.*cm, 36.*cm, 42.*cm, 48.*cm, 54.*cm, 60.*cm, 66.*cm};


  //G4double stripPosY[10] = {-54.*cm, -42.*cm, -30.*cm, -18.*cm, -6.*cm, 6.*cm, 18.*cm, 30.*cm, 42.*cm, 54.*cm};

  //TOF bar placement
  for (G4int i = 0; i < 26; i++)   //Modify for number of slices in WASA detector
  {
    if((i == 0) || (i == 25)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_log, "FriStrip",
			friWall1_log, false, i + 537, checkOverlaps);
    }
    else if((i == 1) || (i == 24)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_1_log, "FriStrip1",
      friWall1_log, false, i + 537, checkOverlaps);
    }
    else if((i == 2) || (i == 23)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_2_log, "FriStrip2",
      friWall1_log, false, i + 537, checkOverlaps);
    }
    else if((i == 3) || (i == 4) || (i == 21) || (i == 22)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_3_log, "FriStrip3",
      friWall1_log, false, i + 537, checkOverlaps);
    }
    else if((i == 5) || (i == 20)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_4_log, "FriStrip4",
      friWall1_log, false, i + 537, checkOverlaps);
    }
    else if(((i > 5) && (i < 10)) || ((i < 20) && (i > 15))){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripLong_log, "FriStripLong",
			friWall1_log, false, i + 537, checkOverlaps);
    }
    else if((i == 11)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_1_cut_log, "FriStripInner1_cut",
      friWall1_log, false, i + 537, checkOverlaps);  
    }
    else if((i == 12)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_2_cut_log, "FriStripInner2_cut",
      friWall1_log, false, i + 537, checkOverlaps);
    }
    else if((i == 13)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_3_cut_log, "FriStripInner3_cut",
      friWall1_log, false, i + 537, checkOverlaps);
    }
    else if((i == 14)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_4_cut_log, "FriStripInner4_cut",
      friWall1_log, false, i + 537, checkOverlaps);
    }
    else if(( i == 10) || (i==15)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_log, "FriStripInner",
			friWall1_log, false, i + 537, checkOverlaps);
    }
  }
  // Place the second FRI wall
  for (G4int i = 26; i < 52; i++)
  {
    if((i == 26) || (i == 51)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_log, "FriStrip",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 27) || (i == 50)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_1_log, "FriStrip1",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 28) || (i == 49)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_2_log, "FriStrip2",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 29) || (i == 30) || (i == 47) || (i == 48)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_3_log, "FriStrip3",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 31) || (i == 46)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStrip_4_log, "FriStrip4",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if(((i > 31) && (i < 36)) || ((i < 46) && (i > 41))){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripLong_log, "FriStripLong",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 37)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_1_cut_log, "FriStripInner1_cut",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 38)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_2_cut_log, "FriStripInner2_cut",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 39)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_3_cut_log, "FriStripInner3_cut",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if((i == 40)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_4_cut_log, "FriStripInner4_cut",
      friWall2_log, false, i + 537, checkOverlaps);
    }
    else if(( i == 36) || (i==41)){
      new G4PVPlacement(0, G4ThreeVector(stripPosX[i], 0., 0.), friStripInner_log, "FriStripInner",
      friWall2_log, false, i + 537, checkOverlaps);
    }
  }
  

  // Position of the FRI wall
  G4double friPosZ[2] = {260.*cm, 270.*cm};  //replaces downstream pizzas from conceptual design - Position of FRI 
  //G4double friPosZ[2] = {90.*cm, 100.*cm};  //OLD - too close to dipole
  G4double friPhi[2] = {0.*deg, 90.*deg}; // Rotation of 2 detectors

  // Define the 2 fri walls in an array
  G4LogicalVolume* friWall_log[2];
  friWall_log[0] = friWall1_log; // First wall
  friWall_log[1] = friWall2_log; // Second wall
  
    for (G4int i = 0; i < 2; i++)
  {
     G4RotationMatrix* rot3 = new G4RotationMatrix;
     rot3->rotateZ(friPhi[i]);

     new G4PVPlacement(rot3, G4ThreeVector(0., 0., friPosZ[i]), friWall_log[i],
                       "FriWall", expHall_log, false, i, checkOverlaps);
  }

  //------------------------------------------------------------------------


  //------------------------------------------------------------------------  
  //TOF wall goes here
  // TofWall = mother volume for the TOF wall
  // Create volume as for the tracker; rectangle (with beam hole extruded?)
  // Half-length along x, y and z
  //Two layers? Rotated 90 degrees wrt each other?
  //G4double tofHalfX = 50.*cm;
  G4double tofHalfX = 106.*cm;
  G4double tofHalfY = 131.*cm;
  G4double tofHalfZ = 1.1*cm;

  G4VSolid* tofWall = new G4Box("TofWall",
                                tofHalfX, tofHalfY, tofHalfZ);

  G4LogicalVolume* tofWall_log = new G4LogicalVolume(tofWall, air, "TofWall");

  // TOF bars to be placed inside the TOF wall
  //G4double tBarHalfX = 4.*cm;
  G4double tBarHalfX = 6.*cm;
  G4double tBarHalfY = 130.*cm;
  G4double tBarShortHalfY = 58.25*cm;
  G4double tBarShortUpperHalfY = 57.25*cm;
  G4double tBarHalfZ = 1.0*cm;

  G4VSolid* tofBar = new G4Box("TofBar",
			       tBarHalfX, tBarHalfY, tBarHalfZ);

  G4LogicalVolume* tofBar_log = new G4LogicalVolume(tofBar, plastic, 
                                                        "TofBar");

  G4VSolid* tofBarShort = new G4Box("TofBarShort",
			       tBarHalfX, tBarShortHalfY, tBarHalfZ);

  G4LogicalVolume* tofBarShort_log = new G4LogicalVolume(tofBarShort, plastic, 
                                                        "TofBarShort");

  G4VSolid* tofBarShortUpper = new G4Box("TofBarShortUpper",
			       tBarHalfX, tBarShortUpperHalfY, tBarHalfZ);

  G4LogicalVolume* tofBarShortUpper_log = new G4LogicalVolume(tofBarShortUpper, plastic, 
							      "TofBarShortUpper");

  //G4double barPosX[10] = {-45.*cm, -35.*cm, -25.*cm, -15.*cm, -5.*cm, 5.*cm, 15.*cm, 25.*cm, 35.*cm, 45.*cm};
  //G4double barPosX[12] = {-60.*cm, -48.*cm, -36.*cm, -24.*cm, -12.*cm, 0.*cm, 0.*cm, 12.*cm, 24.*cm, 36.*cm, 48.*cm, 60.*cm};
  G4double barPosX[18] = {-96.*cm, -84.*cm, -72.*cm, -60.*cm, -48.*cm, -36.*cm, -24.*cm, -12.*cm, 0.*cm, 0.*cm, 12.*cm, 24.*cm, 36.*cm, 48.*cm, 60.*cm, 72.*cm, 84.*cm, 96.*cm};
  
  
  //TOF bar placement
  for (G4int i = 0; i < 18; i++) //20/6
    //for (G4int i = 0; i < 12; i++)   //Modify for number of bars used from WASA detector
    {
  // if(i == 5){
      if(i == 8){//20/6
    //lower TOF bar
        new G4PVPlacement(0, G4ThreeVector(barPosX[i], -72.*cm, 0.), tofBarShort_log, "TofBarShort", 
                          tofWall_log, false, i+589, checkOverlaps);
      }
      //else if(i == 6){
    else if(i == 9){//20/6
    
        //upper TOF bar
        new G4PVPlacement(0, G4ThreeVector(barPosX[i], 72.*cm, 0.), tofBarShortUpper_log, "TofBarShortUpper",
        tofWall_log, false, i+589, checkOverlaps);
      }
      else{
        new G4PVPlacement(0, G4ThreeVector(barPosX[i], 0., 0.), tofBar_log, "TofBar",
                      tofWall_log, false, i+589, checkOverlaps);
      }
    }


  
  G4double tofPosZ = 600.*cm;  //replace downstream pizzas

  new G4PVPlacement(0, G4ThreeVector(0., 0., tofPosZ), tofWall_log,
		    "TofWall", expHall_log, false, 0, checkOverlaps);
  //------------------------------------------------------------------------
    
  //------------------------------------------------------------------------      
  // Solenoid - no longer part of KFM implementation
  G4double solenoid_inrRad = 45.*cm;
  G4double solenoid_outerRad = 55.*cm;
  G4double solenoid_halfThk = 100.*cm;

  G4double solenoidPosZ = 850.*cm;

  G4VSolid* solenoid = new G4Tubs("Solenoid", solenoid_inrRad, 
                                  solenoid_outerRad, solenoid_halfThk, 
                                  0., 360.*deg);

  G4LogicalVolume* solenoid_log = new G4LogicalVolume(solenoid, aluminium, 
                                                      "Solenoid");

  //uncomment to add solenoid
  //new G4PVPlacement(0, G4ThreeVector(0., 0., solenoidPosZ), solenoid_log, 
                    //"Solenoid", expHall_log, false, 11, checkOverlaps);
  //------------------------------------------------------------------------  

  
  // Visualisation attributes
  expHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  /*
  G4VisAttributes* gdTubeAtt = new G4VisAttributes(G4Colour(.5, .5, .5, .4));
  gdTubeAtt->SetDaughtersInvisible(false);
  gdTubeAtt->SetForceSolid(true);
  gdTube_log->SetVisAttributes(gdTubeAtt);
  */
  G4VisAttributes* dipoleAtt = new G4VisAttributes(G4Colour(.7, .7, .7));
  dipoleAtt->SetDaughtersInvisible(false);
  dipoleAtt->SetForceSolid(true);
  dipole_log->SetVisAttributes(dipoleAtt);

  G4VisAttributes* trkrAtt = new G4VisAttributes(G4Colour(0., 0., 1.));
  trkrAtt->SetDaughtersInvisible(false);
  trkrAtt->SetForceSolid(true);
  //trkr_log->SetVisAttributes(trkrAtt);
  trkStraw_log->SetVisAttributes(trkrAtt);

  pizza1_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  pizza2_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes* pizzaSliceAtt = new G4VisAttributes(G4Colour(1., 0., 0.));
  pizzaSliceAtt->SetDaughtersInvisible(false);
  pizzaSliceAtt->SetForceSolid(true);
  pizzaSlice_log->SetVisAttributes(pizzaSliceAtt);

  G4VisAttributes* friAtt = new G4VisAttributes(G4Colour(1., 0., 0.));
  friAtt->SetDaughtersInvisible(false);
  friAtt->SetForceSolid(true);
  friStrip_log->SetVisAttributes(friAtt);
  friStripLong_log->SetVisAttributes(friAtt);
  friStripInner_log->SetVisAttributes(friAtt);
  friStrip_1_log->SetVisAttributes(friAtt);
  friStrip_2_log->SetVisAttributes(friAtt);
  friStrip_3_log->SetVisAttributes(friAtt);
  friStrip_4_log->SetVisAttributes(friAtt);
  friStripInner_1_cut_log->SetVisAttributes(friAtt);
  friStripInner_2_cut_log->SetVisAttributes(friAtt);
  friStripInner_3_cut_log->SetVisAttributes(friAtt);
  friStripInner_4_cut_log->SetVisAttributes(friAtt);
  
  //Remove Solenoid
  G4VisAttributes* solenoidAtt = new G4VisAttributes(G4Colour(.7, .7, .7));
  solenoidAtt->SetDaughtersInvisible(false);
  solenoidAtt->SetForceSolid(true);
  solenoid_log->SetVisAttributes(solenoidAtt);

   G4VisAttributes* tofAtt = new G4VisAttributes(G4Colour(0., 1., 0.));
   tofAtt->SetDaughtersInvisible(false);
   tofAtt->SetForceSolid(true);
   //tofWall_log->SetVisAttributes(tofAtt);
   tofBar_log->SetVisAttributes(tofAtt);
   tofBarShort_log->SetVisAttributes(tofAtt);
   tofBarShortUpper_log->SetVisAttributes(tofAtt);

  return expHall_phys;
}

void JLabKDetectorConstruction::ConstructSDandField()
{
  if (!fMagneticField) fMagneticField = new JLabKMagneticField();
  fMagneticFieldMgr =
   G4TransportationManager::GetTransportationManager()->GetFieldManager();

  fMagneticFieldMgr->SetDetectorField(fMagneticField);
  fMagneticFieldMgr->CreateChordFinder(fMagneticField);
}

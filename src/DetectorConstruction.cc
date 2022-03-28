//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fNofLayers(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
	// Lead material defined using NIST Manager
	auto nistManager = G4NistManager::Instance();
	nistManager->FindOrBuildMaterial("G4_Si");

	// Vacuum
	G4double a, z, density;
	new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);

	// Print materials
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
	// Geometry parameters
	G4double annThick  = 0.3*mm;
	G4double annHole   = 6./2*mm;
	G4double annRadius = 23.9/2*mm;

	G4double cDetThick = 0.3*mm;
	G4double cDetRad   = sqrt(300./pi)*mm;

	auto worldSizeXY = 1.2*m;
	auto worldSizeZ  = 1.2*m; 

	//double annularPositionZ, circularPositionZ; 

	extern double annularPositionZ, circularPositionZ;

	// Get materials
	auto defaultMaterial = G4Material::GetMaterial("Galactic");
	auto detectorMaterial = G4Material::GetMaterial("G4_Si");

	if ( ! defaultMaterial || ! detectorMaterial )
	{
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined."; 
		G4Exception("DetectorConstruction::DefineVolumes()",
		"MyCode0001", FatalException, msg);
	}  

	//     
	// World
	//
	auto worldS  = new G4Box("World",worldSizeXY/2, worldSizeXY/2, worldSizeZ/2);
	auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");
	auto worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"World",0,false,0,fCheckOverlaps);

	//                               
	// Annular detector
	//  
	auto annularDetector
	= new G4Tubs("annDet",     // its name
	annHole,
	annRadius,
	annThick/2,
	0,
	360*deg); // its size

	auto annularDetectorLV
	= new G4LogicalVolume(
	annularDetector,     // its solid
	detectorMaterial,  // its material
	"annularDetectorLV");   // its name
		   
	new G4PVPlacement(
	0,                // no rotation
	G4ThreeVector(0,0,annularPositionZ*mm),  // at (0,0,0)
	annularDetectorLV,          // its logical volume                         
	"annularDetector",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	0,                // copy number
	fCheckOverlaps);  // checking overlaps 


	//                               
	// Circular detector
	//  
	auto circularDetector
	= new G4Tubs("Detectorimeter",     // its name
	0,
	cDetRad,
	cDetThick/2,
	0,
	360*deg); // its size

	auto circularDetectorLV
	= new G4LogicalVolume(
	circularDetector,     // its solid
	detectorMaterial,  // its material
	"circularDetectorLV");   // its name
		   
	new G4PVPlacement(
	0,                // no rotation
	G4ThreeVector(0,0,circularPositionZ*mm),  // at (0,0,0)
	circularDetectorLV,          // its logical volume                         
	"circularDetector",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	0,                // copy number
	fCheckOverlaps);  // checking overlaps 


	//                               
	// Foil - for visuals only
/*	//  
	auto foilTube
	= new G4Tubs("Detectorimeter",     // its name
	0,
	3.*mm,
	0.01*mm,
	0,
	360*deg); // its size

	auto foilLV
	= new G4LogicalVolume(
	foilTube,     // its solid
	defaultMaterial,  // its material
	"foilLV");   // its name
		   
	new G4PVPlacement(
	0,                // no rotation
	G4ThreeVector(),  // at (0,0,0)
	foilLV,          // its logical volume                         
	"foilTube",    // its name
	worldLV,          // its mother  volume
	false,            // no boolean operation
	0,                // copy number
	fCheckOverlaps);  // checking overlaps 
*/

	//                                        
	// Visualization attributes
	//
	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

	auto detectorVisAtt = new G4VisAttributes(G4Colour(180./255,125./255,0./255,0.7));
	auto foilVisAtt     = new G4VisAttributes(G4Colour(1.,1.,1.,0.7));
	//detectorVisAtt->SetVisibility(true);
	annularDetectorLV->SetVisAttributes(detectorVisAtt);
	circularDetectorLV->SetVisAttributes(detectorVisAtt);
//	foilLV->SetVisAttributes(foilVisAtt);

	//
	// Always return the physical World
	//
	return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // 
  // Sensitive detectors
  //
  auto annularSD = new DetectorSD("annularSD","annularHitsCollection",1);
  G4SDManager::GetSDMpointer()->AddNewDetector(annularSD);
  SetSensitiveDetector("annularDetectorLV",annularSD);

  auto circularSD = new DetectorSD("circularSD","circularHitsCollection",1);
  G4SDManager::GetSDMpointer()->AddNewDetector(circularSD);
  SetSensitiveDetector("circularDetectorLV",circularSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

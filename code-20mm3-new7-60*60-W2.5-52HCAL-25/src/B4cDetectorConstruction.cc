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
// $Id$
// 
/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class
/*
***************************************MAP FOR TILE NUMBERS************************************

	[00][10][20][30][40][50]
	[01][11][21][31][41][51]
	[02][12][22][32][42][52]
	[03][13][23][33][43][53]
	[04][14][24][34][44][54]
	[05][15][25][35][45][55]

*/
#include "B4cDetectorConstruction.hh"
#include "B4cCalorimeterSD.hh"
#include "EMSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>

const double offset=0.5;//in mm
const double dist=1.0;
const double tot_len=25.;
const double h=0.5*sqrt(3)*dist;

const int nx1=int((tot_len-2*offset)/(dist/2));//+1
const int ny1=int((tot_len-offset)/(2*h));//+1
const int ny2=int((tot_len-offset-h)/(2*h));//+1

const double x0=-((tot_len/2.0)-offset);
const double y01=((tot_len/2.0)-offset);
const double y02=((tot_len/2.0)-offset-h);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fMessenger(0),
   fMagField(0),
   fCheckOverlaps(false)
{
  // Define /B4/det commands using generic messenger class
  fMessenger 
    = new G4GenericMessenger(this, "/B4/det/", "Detector construction control");

  // Define /B4/det/setMagField command
  G4GenericMessenger::Command& setMagFieldCmd
    = fMessenger->DeclareMethod("setMagField", 
                                &B4cDetectorConstruction::SetMagField, 
                                "Define magnetic field value (in X direction");
  setMagFieldCmd.SetUnitCategory("Magnetic flux density");                                
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{ 
  delete fMagField;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_Fe",fromIsotopes);
  nistManager->FindOrBuildMaterial("G4_POLYSTYRENE",fromIsotopes);
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
    G4int nofLayers2 = 52;//36;
    G4double calorSizeXY  = 1000.*mm;
    G4double hadcell=10.*cm;
    G4double emcell=2.5*cm;//WScFi
    G4double absoThickness2 = 20. *mm;
    G4double gapThickness2 = 3. *mm;//used to be 6.0 mm
    G4double num_of_towers=6;//HCAL
    G4double num_of_towers_EM=24;//WScFi
    G4double layerThickness2 = absoThickness2+gapThickness2;
    G4double calorThickness2 =nofLayers2*layerThickness2;
    G4double worldSizeXY = 5.0* calorSizeXY;
    G4double calorEMZ=17.0*cm; 
    G4double worldSizeZ  = 5.0 *( calorThickness2+calorEMZ);
  // Get materials
    G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
    G4Material* absorberMaterial = G4Material::GetMaterial("G4_Pb");
    G4Material* gapMaterial = G4Material::GetMaterial("liquidArgon");
    G4Material* gapMaterial2 =G4Material::GetMaterial("G4_POLYSTYRENE");
    G4Material* absorberMaterial2= G4Material::GetMaterial("G4_Fe");
    
    
    G4double a=183.85*g/mole;
    G4Element* elW=new G4Element("Tungsten","W",74.,a);
    
    G4Material* EMCal_abs_mat=new G4Material("EMCal_fiber_mat",14.44*g/cm3,2);
    EMCal_abs_mat->AddElement(elW, 74.8*perCent);//96 74.8% W + 25.2
    EMCal_abs_mat->AddMaterial(gapMaterial2, 25.2*perCent);//4
    
  
    if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ||!gapMaterial2 ||!absorberMaterial2|| !EMCal_abs_mat) {
        G4cerr << "Cannot retrieve materials already defined. " << G4endl;
        G4cerr << "Exiting application " << G4endl;
        exit(1);
    }
    G4MaterialPropertiesTable* mtpt = new G4MaterialPropertiesTable();
    gapMaterial2->SetMaterialPropertiesTable(mtpt);
    gapMaterial2->GetIonisation()->SetBirksConstant(.2*mm/MeV);
  //     
  // World
  //
    G4VSolid* worldS= new G4Box("World",worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // name, its size
                         
    G4LogicalVolume* worldLV= new G4LogicalVolume(worldS, defaultMaterial, "World");         // solid, material, its name
                                   
    G4VPhysicalVolume* worldPV = new G4PVPlacement( 0, G4ThreeVector(), worldLV, "World", 0, false, 0, fCheckOverlaps);  //no rotation, at (0,0,0),its logical volume ,its name, its mother  volume, no boolean operation, copy number, checking overlaps
  
  //                               
  // Calorimeter
  //
    G4LogicalVolume* calorLV[6][6];
    
    char bu[200];
    char bu2[200];
    char bu3[200];

    for(int i=0;i<num_of_towers;i++){

        for(int j=0;j<num_of_towers;j++){
            sprintf(bu,"Calorimeter%d%d",i,j);
            G4VSolid* calorimeterS= new G4Box("Calorimeter", hadcell/2, hadcell/2, calorThickness2/2); // its name, its size
           
            calorLV[i][j] = new G4LogicalVolume(calorimeterS, defaultMaterial, bu);   //its solid, material, its name

                                   
            new G4PVPlacement(0, G4ThreeVector(-25.*cm+i*10.*cm, 45.*cm-j*10.*cm,calorEMZ/2+calorThickness2/2), calorLV[i][j],bu, worldLV, false, 0, fCheckOverlaps);  // no rotation, at (0,0,0), its logical volume,  its name, its mother  volume, no boolean operation, copy number,   checking overlaps
        }

    }


  //                                 
  // Layer
  //

    G4LogicalVolume* layerLV[6][6];

    for(int i=0;i<num_of_towers;i++){
	
        for(int j=0;j<num_of_towers;j++){
            sprintf(bu,"Layer%d%d",i,j);
            sprintf(bu2,"Layer2%d%d",i,j);
            G4VSolid* layerS = new G4Box("Layer", hadcell/2, hadcell/2, layerThickness2/2); //its name, its size
                         
            layerLV[i][j] = new G4LogicalVolume( layerS, defaultMaterial, bu);         // its solid, its material, its name
            new G4PVReplica( "Layer", layerLV[i][j], calorLV[i][j], kZAxis,nofLayers2,layerThickness2);  // name, logical volume, mother, axis of replication, number of replica, witdth of replica
        }
    }


  //                               
  // Absorber
  //
    G4LogicalVolume* absorberLV[6][6];

    for(int i=0;i<num_of_towers;i++){
        for(int j=0;j<num_of_towers;j++){
            sprintf(bu,"Abso%d%d",i,j);
            sprintf(bu2,"Abso2%d%d",i,j);

            sprintf(bu2,"Abso2%d%d",i,j);

            G4VSolid* absorberS= new G4Box("Abso", hadcell/2, hadcell/2, absoThickness2/2); // its name, its size
                         
            absorberLV[i][j] = new G4LogicalVolume(absorberS,absorberMaterial2,bu);          //its solid, material,  its name
                            

       
            new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness2/2), // its position
                 absorberLV[i][j],       // its logical volume                         
                 "Abso",           // its name
                 layerLV[i][j],          // its mother  volume
                 false,            // no boolean operation
                 i+j,                // copy number
                 fCheckOverlaps);  // checking overlaps 

	}
}

  //                               
  // Gap
  //
    char buffer [200];
    G4LogicalVolume* gapLV[6][6];
    for(int i=0;i<num_of_towers;i++){


        for(int j=0;j< num_of_towers;j++){
            sprintf(buffer,"gap%d%d",i,j);
            sprintf(bu2,"gap2%d%d",i,j);

            G4VSolid* gapS = new G4Box("Gap", hadcell/2, hadcell/2, gapThickness2/2); // its name, its size

                         
            gapLV[i][j] = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial2,      // its material
                 buffer);           // its name
                           
            new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness2/2), // its position
                 gapLV[i][j],            // its logical volume                         
                 buffer,            // its name
                 layerLV[i][j],          // its mother  volume
                 false,            // no boolean operation
                 i+j,                // copy number
                 fCheckOverlaps);  // checking overlaps 

        }
    }
    
    //EM

//EM LOGIC

G4LogicalVolume* calorEM[24][24];
G4int copynono=0;
for(int i=0;i<num_of_towers_EM;i++){
	for(int j=0;j<num_of_towers_EM;j++){
		sprintf(bu,"CalorimeterEM%02d%02d",i,j);
		  G4VSolid* calorimeterEM = new G4Box("CalorimeterEM_",     // its name
                 emcell/2, emcell/2, calorEMZ/2);
calorEM[i][j] = new G4LogicalVolume(calorimeterEM,defaultMaterial,bu);
        new G4PVPlacement(0,G4ThreeVector(-28.75*cm+i*2.5*cm , 48.75*cm-j*2.5*cm ,0.),calorEM[i][j],bu,worldLV,false,0,fCheckOverlaps);
	}
}

    
    //absorber


G4LogicalVolume* absorberEMLV[24][24];
for(int i=0;i<num_of_towers_EM;i++){
	for(int j=0;j<num_of_towers_EM;j++){
		sprintf(bu,"AbsoEM%02d%02d",i,j);
  G4VSolid* absorberEM
    = new G4Box("AbsoEM",            // its name
                 emcell/2, emcell/2, calorEMZ/2); // its size
absorberEMLV[i][j] = new G4LogicalVolume(absorberEM, EMCal_abs_mat,bu);
new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),absorberEMLV[i][j],"AbsoEM_p", calorEM[i][j],false,0,fCheckOverlaps);
	}

}

/*
    //Fibers
    G4LogicalVolume* gapEMLV[24][24];
   G4VSolid* gapEM = new G4Tubs("GapEM",             // its name
                                 0.0, 0.235*mm, calorEMZ/2,0.0,360.0 * deg); // its size//0.0 * deg, 360.0 * deg
    double step_x=(dist/2.0)*mm;
    double step_y=(2.0*h)*mm;
    for(int sx=0;sx<num_of_towers_EM;sx++){
   	for(int sy=0;sy<num_of_towers_EM;sy++){
        		sprintf(bu3,"gapEM%02d%02d",sx,sy);
        		gapEMLV[sx][sy] = new G4LogicalVolume(gapEM,gapMaterial2,bu3);

       		copynono=0;
       		for(int i=0;i<nx1;i++){
       			G4double pos_x=x0*mm+i*step_x;
        			G4double pos_y= 0.0 ;
        			int jend=(i%2==0) ? ny1 : ny2;
        				for(int j=0;j<jend;j++){   
           				sprintf(bu2,"gap2%02d%02d%03d%03d",sx, sy, i,j);
            				if(i%2==0) pos_y=y01*mm-j*step_y;
            				if(i%2!=0) pos_y=y02*mm-j*step_y;
            				new G4PVPlacement(0,G4ThreeVector(pos_x,pos_y, 0.),gapEMLV[sx][sy],bu2, absorberEMLV[sx][sy],false, copynono,fCheckOverlaps);
            				copynono++;
        				}
    				}

     			}
     		}
*/
  //
  // print parameters
  //
 /* G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << nofLayers << " layers of: [ "
         << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
         << " + "
         << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
  */
  
  // 
  // Sensitive detectors
  //

    B4cCalorimeterSD* absoSD[6][6];
    B4cCalorimeterSD* gapSD[6][6];


    char buffer2 [200];
    char buffer3 [200];
    char buffer4 [200];
    char buffer5 [200];
    char bu4[200];
    char bu5[200];




    for(int i=0;i<num_of_towers;i++){
        for(int j=0;j<num_of_towers;j++){
            sprintf(buffer2,"AbsorberHitsCollection%02d%02d",i,j);
            sprintf(buffer3,"AbsorberSD%02d%02d",i,j);
            absoSD[i][j] = new B4cCalorimeterSD(buffer3, buffer2, nofLayers2);
            G4SDManager::GetSDMpointer()->AddNewDetector(absoSD[i][j] );
            absorberLV[i][j]->SetSensitiveDetector(absoSD[i][j]);


            sprintf(buffer4,"GapHitsCollection%02d%02d",i,j);
            sprintf(buffer5,"GapSD%02d%02d",i,j);
   
            gapSD[i][j] = new B4cCalorimeterSD(buffer5, buffer4, nofLayers2);
            G4SDManager::GetSDMpointer()->AddNewDetector(gapSD[i][j] );
            gapLV[i][j]->SetSensitiveDetector(gapSD[i][j]);
 
        }
    }

   



//B4cCalorimeterSD* gap2EMSD[24][24];
B4cCalorimeterSD* absoEMSD[24][24];
copynono=0;
for(int i=0;i<num_of_towers_EM;i++){
	for(int j=0;j<num_of_towers_EM;j++){
/*
sprintf(buffer4,"gapEMHitsCollection%02d%02d",i,j);
sprintf(buffer5,"gapEMSD2%02d%02d",i,j);
copynono++;
  gap2EMSD[i][j] = new B4cCalorimeterSD(buffer5, buffer4, 1);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap2EMSD[i][j]);
  gapEMLV[i][j]->SetSensitiveDetector(gap2EMSD[i][j]);
*/

sprintf(buffer2,"AbsorberEMHitsCollection%02d%02d",i,j);
sprintf(buffer3,"AbsorberEMSD%02d%02d",i,j);
  absoEMSD[i][j] = new B4cCalorimeterSD(buffer3, buffer2,1);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoEMSD[i][j]);
  absorberEMLV[i][j]->SetSensitiveDetector(absoEMSD[i][j]);
	}

}
   



 
  //                                        
  // Visualization attributes
  //
    worldLV->SetVisAttributes (G4VisAttributes::Invisible);
    G4VisAttributes invis=G4VisAttributes::Invisible;
    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1,0,0));//red
    
    G4VisAttributes* abs2VisAtt= new G4VisAttributes(G4Colour(1,0,1));//magenta
    G4VisAttributes* absVisAtt= new G4VisAttributes(G4Colour(0,0,1));//blue
    
    G4VisAttributes* gap2VisAtt= new G4VisAttributes(G4Colour(1,1,0));//yellow
    G4VisAttributes* gapVisAtt= new G4VisAttributes(G4Colour(0,1,1));//cyan
    
    G4VisAttributes* hcal2= new G4VisAttributes(G4Colour(0.5,0.5,0.5));//gray
    
    G4bool draw_hcal=true;
    simpleBoxVisAtt->SetVisibility(draw_hcal);
    absVisAtt->SetVisibility(draw_hcal);
    hcal2->SetVisibility(draw_hcal);
    gap2VisAtt->SetVisibility(draw_hcal);
    gapVisAtt->SetVisibility(draw_hcal);
    abs2VisAtt->SetVisibility(draw_hcal);
    absVisAtt->SetVisibility(draw_hcal);
    for(int i=0;i<num_of_towers;i++){
        for(int j=0;j<num_of_towers;j++){
            
            calorLV[i][j]->SetVisAttributes(simpleBoxVisAtt);
            
            gapLV[i][j]->SetVisAttributes(gapVisAtt);
            absorberLV[i][j]->SetVisAttributes(absVisAtt);
                     
            layerLV[i][j]->SetVisAttributes(invis);
        }
    }
    
    //G4Colour myColourMag(1, 0 , 1, 0.5);
    G4VisAttributes* calorEMvis= new G4VisAttributes(G4Colour(1,0,1));//magenta calorimeter
 //   G4VisAttributes* gapEMvis= new G4VisAttributes(G4Colour(1,1,0));//yellow gap
    G4VisAttributes* absEMvis= new G4VisAttributes(G4Colour(0,1,1));//cyan absorber
    
    calorEMvis->SetVisibility(true);
    calorEMvis->SetForceWireframe (true);
 //   gapEMvis->SetVisibility(false);
 //   gapEMvis->SetForceSolid (true);
    absEMvis->SetVisibility(true);
    absEMvis->SetForceWireframe(true);

   for(int i=0;i<num_of_towers_EM;i++){
        for(int j=0;j<num_of_towers_EM;j++){
    calorEM[i][j]->SetVisAttributes(calorEMvis);
   // gapEMLV[i][j]->SetVisAttributes(gapEMvis);
    absorberEMLV[i][j]->SetVisAttributes(absEMvis);
        }
    }

  //
  // Always return the physical World
  //
    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::SetMagField(G4double fieldValue)
{
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Delete the existing magnetic field
  if ( fMagField )  delete fMagField; 

  if ( fieldValue != 0. ) {
    // create a new one if not null
    fMagField 
      = new G4UniformMagField(G4ThreeVector(fieldValue, 0., 0.));
      
    fieldManager->SetDetectorField(fMagField);
    fieldManager->CreateChordFinder(fMagField);
  } 
  else {
    fMagField = 0;
    fieldManager->SetDetectorField(fMagField);
  }
}

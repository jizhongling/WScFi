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
/// \file B4cCalorimeterSD.cc
/// \brief Implementation of the B4cCalorimeterSD class

#include <iostream>

#include "EMSD.hh"//

#include "B4cCalorHit.hh"
#include "EMSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"//
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"//
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"


#include "G4ios.hh"//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EMSD::EMSD(G4String name)
 : G4VSensitiveDetector(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EMSD::~EMSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EMSD::Initialize(G4HCofThisEvent* hce)
{
    EvisF1Tile = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool EMSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
    
    G4double stepLength = 0.;
    if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
        stepLength = step->GetStepLength();
    }
//G4cout<<"Edep in EM: "<<edep<<G4endl;
  if ( edep==0.) return false;

  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
  
    //G4VPhysicalVolume* volume =step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();//added this
    G4VPhysicalVolume* volume=touchable->GetVolume();
    
    if(strcmp(volume->GetName(),"EMGapPhysical")==0){
        //G4int F1LArGapId = physVol->GetCopyNo();
        //G4int F1TileId = EmModule->GetF1TileID(F1LArGapId);
        G4Material *mat = volume->GetLogicalVolume()->GetMaterial();
        G4double charge =step->GetTrack()->GetDefinition()->GetPDGCharge();
        G4double birk=mat->GetIonisation()->GetBirksConstant();
        if(birk*edep*stepLength*charge !=0){
            edep=edep/(1.+birk*edep/stepLength);
        }

//        G4cout<<"hits: "<<EvisF1Tile<<" "<<volume->GetCopyNo()<<G4endl;
        EvisF1Tile +=  edep;
    };
 
// G4cout<<"Edep in EM: "<<edep<<G4endl;     
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EMSD::EndOfEvent(G4HCofThisEvent*)
{

    //G4cout<<"SD: "<<EvisF1Tile/CLHEP::keV <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

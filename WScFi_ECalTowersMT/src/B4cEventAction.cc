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
/// \file B4cEventAction.cc
/// \brief Implementation of the B4cEventAction class

#include "B4cEventAction.hh"
#include "B4cCalorimeterSD.hh"
#include "B4cCalorHit.hh"
#include "B4Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"
#include "EMSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4UniformMagField.hh"

#include "Randomize.hh"
#include <iomanip>
const G4int nofLayers=52;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::B4cEventAction(EMSteppingAction* SA)
 : G4UserEventAction(),
   fMessenger(0),
   fPrintModulo(1),
  StepAction(SA)
{
  // // Define /B4/event commands using generic messenger class
  // fMessenger = new G4GenericMessenger(this, "/B4/event/", "Event control");

  // // Define /B4/event/setPrintModulo command
  // G4GenericMessenger::Command& setPrintModulo 
  //   = fMessenger->DeclareProperty("setPrintModulo", 
  //                                 fPrintModulo, 
  //                                "Print events modulo n");
  // setPrintModulo.SetRange("value>0");                                
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::~B4cEventAction()
{
  //delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHitsCollection* 
B4cEventAction::GetHitsCollection(const G4String& hcName,
                                  const G4Event* event) const
{
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(hcName);
  B4cCalorHitsCollection* hitsCollection 
    = static_cast<B4cCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4cerr << "Cannot access hitsCollection " << hcName << G4endl;
    exit(1);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength,double number1,double number2) const
{
  // print event statistics
  //G4cout
//	<<number1<<number2<<" "<<
//  G4cout   << "   Absorber: total energy: "
//	<< "Number of the tower: "<<number1<<number2 
  //   << std::setw(7) << G4BestUnit(absoEdep, "Energy")
   //  << "       total track length: " 
    // << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
    // << G4endl
    // << "        Gap: total energy: " 
     G4cout << gapEdep/MeV <<G4endl;
    // << "       total track length: " 
    // std::setw(7) << G4BestUnit(gapTrackLength, "Length")
      //G4endl;
}

void B4cEventAction::PrintTileStatistics(G4double absoEdep, G4double absoTrackLength, G4double gapEdep, G4double gapTrackLength, int tile_num) const
{

	G4cout<<gapEdep/MeV;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::BeginOfEventAction(const G4Event* event)
{  
  // G4FieldManager* fieldManager
  //   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  G4int eventID = event->GetEventID();
  if ( eventID % fPrintModulo == 0 )  { 
    G4cout << "---> Begin of event: " << eventID+1 << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
    StepAction->initialize(eventID);

    // const G4double p[4] = {0.,0.,0.,0.};
    // G4double b[3];
    // fieldManager->GetDetectorField()->GetFieldValue(p, b);
    // G4cout<<"Mag field is: "<<b[0]/tesla<<G4endl;

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections

    B4cCalorHitsCollection* absoHC[6][6];
    B4cCalorHitsCollection* gapHC[6][6];
//added another dimension
    B4cCalorHit* absoHit[6][6][nofLayers+1];
    B4cCalorHit* gapHit[6][6][nofLayers+1];
    char buffer[200];
    char buffer2[200];


    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            sprintf(buffer,"AbsorberHitsCollection%02d%02d",i,j);
            sprintf(buffer2,"GapHitsCollection%02d%02d",i,j);
            absoHC[i][j] = GetHitsCollection(buffer, event);
            gapHC[i][j] = GetHitsCollection(buffer2, event);

  // Get hit with total values
            absoHit[i][j][nofLayers] = (*absoHC[i][j])[absoHC[i][j]->entries()-1];
            gapHit[i][j][nofLayers] = (*gapHC[i][j])[absoHC[i][j]->entries()-1];
//now get the individula for each cell
            for(int k=0;k<nofLayers;k++){
                absoHit[i][j][k] = (*absoHC[i][j])[k];
                gapHit[i][j][k] = (*gapHC[i][j])[k];

            }


        	  }
    	}



//add EM hist
	B4cCalorHitsCollection* absoEMHC[24][24];
	B4cCalorHitsCollection* gapEMHC[24][24];
	B4cCalorHit* absoEMHit[24][24];
	B4cCalorHit* gapEMHit[24][24];

   	for(int i=0;i<24;i++){
       		for(int j=0;j<24;j++){


         	sprintf(buffer,"AbsorberEMHitsCollection%02d%02d",i,j);
         	sprintf(buffer2,"gapEMHitsCollection%02d%02d",i,j);
		absoEMHC[i][j]=GetHitsCollection(buffer,event);
		gapEMHC[i][j]=GetHitsCollection(buffer2,event);
		absoEMHit[i][j]=(*absoEMHC[i][j])[absoEMHC[i][j]->entries()-1];
		gapEMHit[i][j]=(*gapEMHC[i][j])[gapEMHC[i][j]->entries()-1];

		}
	}


  // Print per event (modulo n)
  //

    G4int eventID = event->GetEventID();
    if ( eventID % fPrintModulo == 0) {
  
        
    G4cout<<"EM: "<<G4endl;

	for(int i=0;i<24;i++){
		for(int j=0;j<24;j++){
        		PrintEventStatistics(absoEMHit[i][j]->GetEdep(),absoEMHit[i][j]->GetTrackLength(),gapEMHit[i][j]->GetEdep(),gapEMHit[i][j]->GetTrackLength(),i,j);
		}
	}
        
	G4cout<<"HAD: "<<G4endl;

	for(int i=0;i<6;i++){
		for(int j=0;j<6;j++){

    			PrintEventStatistics(
      	absoHit[i][j][nofLayers]->GetEdep(), absoHit[i][j][nofLayers]->GetTrackLength(),
      gapHit[i][j][nofLayers]->GetEdep(), gapHit[i][j][nofLayers]->GetTrackLength(),i,j);
			for(int k=0;k<nofLayers;k++){
				PrintTileStatistics(
        absoHit[i][j][k]->GetEdep(), absoHit[i][j][k]->GetTrackLength(),
      gapHit[i][j][k]->GetEdep(), gapHit[i][j][k]->GetTrackLength(),k);
				if(k<(nofLayers-1)) G4cout<<" ";
				if(k==(nofLayers-1)) G4cout<<G4endl;
				}
  			}  
		}
	
	}
 G4cout << "---> End of event: " << eventID+1 << G4endl;


//G4int eventID =event->GetEventID();
//PrintEventStatistics(absohit->GetEdep(),absohit->GetTrackLength(),gaphit->GetEdep(),gaphit->GetTrackLength());  
  // Fill histograms, ntuple
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 //G4double EmEdep  = StepAction->GetEdepEM("FCALEM");
  // fill histograms
  //analysisManager->FillH1(1, EmEdep);
 // analysisManager->FillH1(2, gapHit[0][0]->GetEdep());
 // analysisManager->FillH1(3, absoHit[0][0]->GetTrackLength());
 // analysisManager->FillH1(4, gapHit[0][0]->GetTrackLength());
  
  // fill ntuple
 // analysisManager->FillNtupleDColumn(0, absoHit[0][0]->GetEdep());
 // analysisManager->FillNtupleDColumn(1, gapHit[0][0]->GetEdep());
 // analysisManager->FillNtupleDColumn(2, absoHit[0][0]->GetTrackLength());
 // analysisManager->FillNtupleDColumn(3, gapHit[0][0]->GetTrackLength());
 // analysisManager->AddNtupleRow();  
}  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

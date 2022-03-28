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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "DetectorSD.hh"
#include "DetectorHit.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
 : G4UserEventAction(),
   fAnnularHCID(-1),
   fCircularHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
	auto hitsCollection = static_cast<DetectorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

	if ( ! hitsCollection )
	{
		G4ExceptionDescription msg;
		msg << "Cannot access hitsCollection ID " << hcID; 
		G4Exception("EventAction::GetHitsCollection()",
		"MyCode0003", FatalException, msg);
	}         

	return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{  
	// Get hits collections IDs (only once)
	if ( fAnnularHCID == -1 )
	{
		fAnnularHCID  = G4SDManager::GetSDMpointer()->GetCollectionID("annularHitsCollection");
		fCircularHCID = G4SDManager::GetSDMpointer()->GetCollectionID("circularHitsCollection");
	}

	// Get hits collections
	auto annularHC  = GetHitsCollection(fAnnularHCID, event);
	auto circularHC = GetHitsCollection(fCircularHCID, event);

	// Get hit with total values
	auto annularHit  = (*annularHC)[annularHC->entries()-1];
	auto circularHit = (*circularHC)[circularHC->entries()-1];

	// get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();
/*	analysisManager->FillH1(0, annularHit->GetEdep());
	analysisManager->FillH1(1, circularHit->GetEdep());
	analysisManager->FillH1(2, annularHit->GetTrackLength());
	analysisManager->FillH1(3, circularHit->GetTrackLength());
*/
	// fill ntuple
	analysisManager->FillNtupleDColumn(0, circularHit->GetEdep()/keV);
	analysisManager->FillNtupleDColumn(1, annularHit->GetEdep()/keV);
	analysisManager->FillNtupleDColumn(2, circularHit->GetTrackLength());
	analysisManager->FillNtupleDColumn(3, annularHit->GetTrackLength());
	analysisManager->FillNtupleDColumn(4, 1.);
	analysisManager->FillNtupleDColumn(5, 1.);
	analysisManager->AddNtupleRow();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

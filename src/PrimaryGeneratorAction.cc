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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <G4INCLRandom.hh>

#include <numeric>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(50.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	// Get particle definitions
	auto alphaParticle = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
	auto betaParticle  = G4ParticleTable::GetParticleTable()->FindParticle("e-");
	auto gammaParticle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
	
	// Generate position
	G4ThreeVector position = sourcePosition();
	//G4cout << "\n\n\n\n**********************************************\n\n";
	// Generate alpha decay path based on intensities from input file
	int decayPath = generateAlpha();
	//G4cout << "\n\n My decay path is = " << decayPath;
	if(optionAlpha == 1)
	{
		fParticleGun->SetParticleDefinition(alphaParticle);
		fParticleGun->SetParticleEnergy(alphaEnergy[decayPath]*keV);
		fParticleGun->SetParticlePosition(position);
		fParticleGun->SetParticleMomentumDirection(isotropicSource());
		fParticleGun->GeneratePrimaryVertex(anEvent);
		
		//G4cout << " Alpha energy: " << alphaEnergy[decayPath] << "\n";
	}
	
	// Generate gammas and/or conversion electrons + x rays
	if(optionGamma==1 || optionCE==1)
	{
		int numberOfGammas = gammaIntensity[decayPath].size();
		for(int i=0; i<numberOfGammas; i++)
		{
			// Check transition probability and generate either gamma or CE
			double gammaRandom = G4UniformRand();
			if(gammaRandom<gammaIntensity[decayPath][i])
			{
				// Generate either gamma [transitionNumber=0] or CE [transitionNumber>0]
				int transitionNumber = generateTransition(i);

				//G4cout << " Transition number: " << transitionNumber << "\n";
				
				// If not converted, generate gamma
				if(optionGamma==1 && transitionNumber==0)
				{
					fParticleGun->SetParticleDefinition(gammaParticle);
					fParticleGun->SetParticleEnergy(gammaEnergy[i]*keV);
					fParticleGun->SetParticleMomentumDirection(isotropicSource());
					fParticleGun->GeneratePrimaryVertex(anEvent);
					//G4cout << " Gamma energy: " << gammaEnergy[i] << "\n";
				}				
				// If converted, generate CE
				else if(optionCE==1 && transitionNumber>0)
				{
					fParticleGun->SetParticleDefinition(betaParticle);
					fParticleGun->SetParticleEnergy(ceEnergy[i][transitionNumber-1]*keV);
					fParticleGun->SetParticleMomentumDirection(isotropicSource());
					fParticleGun->GeneratePrimaryVertex(anEvent);
					//G4cout << " CE energy: " << ceEnergy[i][transitionNumber-1] << "\n";
					//testFile << "ceE:" << ceEnergy[i][transitionNumber-1] << "\t";
					
					// Generate x rays following K conversion
					if(optionX==1 && transitionNumber==1)
					{
						
						
						int kXRay = generateKXRay();
						G4cout << " K CONVERSION: " << kXRay << "\n";
						if(kXRay>-1)
						{
							fParticleGun->SetParticleDefinition(gammaParticle);
							fParticleGun->SetParticleEnergy(kEnergy[kXRay]*keV);
							fParticleGun->SetParticleMomentumDirection(isotropicSource());
							fParticleGun->GeneratePrimaryVertex(anEvent);
							//G4cout << " X ray energy: " << kEnergy[kXRay] << "\n";

							// If Ka x ray, generate L x ray from sub shell
							if(kXRay==0 || kXRay==1 || kXRay==2)
							{
								int lXRay = generateLXRay(kXRay);
								if(lXRay>0)
								{
									//G4cout << " Transition number: " << transitionNumber << "\n";
									//G4cout << " L X ray number: "  << kXRay << "   " << lXRay << "\n";
									fParticleGun->SetParticleDefinition(gammaParticle);
									fParticleGun->SetParticleEnergy(lEnergy[kXRay][lXRay]*keV);
									fParticleGun->SetParticleMomentumDirection(isotropicSource());
									fParticleGun->GeneratePrimaryVertex(anEvent);
									//G4cout << " L X ray energy: " << lEnergy[kXRay][lXRay] << "\n";
								}
							}
						}
					}
					
					// Generate x rays following L conversion from sub shell
					if(optionX==1 && transitionNumber>2 && transitionNumber<5)
					{
						int lXRay = generateLXRay(transitionNumber-2);
						//G4cout << " L CONVERSION: " << lXRay << "        " << transitionNumber-2 << "\n";
						if(lXRay>=0)
						{
							fParticleGun->SetParticleDefinition(gammaParticle);
							fParticleGun->SetParticleEnergy(lEnergy[transitionNumber-2][lXRay]*keV);
							fParticleGun->SetParticleMomentumDirection(isotropicSource());
							fParticleGun->GeneratePrimaryVertex(anEvent);
							//G4cout << " L X ray energy: " << lEnergy[transitionNumber-2][lXRay] << "\n";
						}
					}
					
				}
			}
		}
	}
/*
	G4cout << "\n\n**********************************************\n\n";
	G4cout << "\n\n L X RAY 1 2 :  " << lEnergy[0][1] << "\n";
	G4cout <<     " L X RAY 2 3 :  " << lEnergy[1][2] << "\n";
	G4cout <<     " L X RAY 3 4 :  " << lEnergy[2][3] << "\n\n\n";
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Gaussian source profile
G4ThreeVector sourcePosition()
{
	double foilWidth  = 3*mm;
	double sigma     = foilWidth/4.29193;
	double xPos = G4RandGauss::shoot(0,sigma);
	double yPos = G4RandGauss::shoot(0,sigma);
	double zPos = 0;
	//double distrDepth = (0.5/4.29193)*um;
	//double zPos = G4RandGauss::shoot(0,distrDepth);
	//double zPos = (-6.346)*mm;
	
	return G4ThreeVector(xPos,yPos,zPos);
}

// Isotropic source definition
G4ThreeVector isotropicSource()
{
	// Define isotropic source
	G4double pX, pY, pZ;
	G4double pTot;
	do
	{
		pX = (G4UniformRand()-0.5)/0.5;
		pY = (G4UniformRand()-0.5)/0.5;
		pZ = (G4UniformRand()-0.5)/0.5;
		pTot = (pX*pX) + (pY*pY) + (pZ*pZ);
	}while (pTot > 1 || pTot == 0.0);

	pTot = std::sqrt(pTot);
	pX /= pTot;
	pY /= pTot;
	pZ /= pTot;
	
	return G4ThreeVector(pX,pY,pZ);
}

// Returns the decay path for which alpha, and gamma transitions
G4int generateAlpha()
{
	double sumIntensities = std::accumulate(alphaIntensity.begin(), alphaIntensity.end(), 0.0);
	
	double alphaRandom = G4UniformRand();
	double prob = 1.0;
	int choice=0, decayPath;
	
	// Get decay path based on normalised sum of alpha intensities
	do
	{
		decayPath=choice;
		prob -= alphaIntensity[choice]/sumIntensities;
		choice++;
	}while(alphaRandom<prob);
	
	return decayPath;
}

// Returns the transition number for gamma ray (==0), or conversion electron (>0)
G4int generateTransition(int gamma)
{
	// Get conversion coefficient for gamma
	double totalConvCoeff = std::accumulate(ceIntensity[gamma].begin(), ceIntensity[gamma].end(), 0.0);

	// Default generate gamma
	int transition=0;
	
	// Check if transition is converted, if so choose CE for transition
	double randomNumber = G4UniformRand();
	double prob = 1.0;
	if(randomNumber < totalConvCoeff/(totalConvCoeff+1.))
	{
		randomNumber = G4UniformRand();
		do
		{
			transition++;
			prob -= ceIntensity[gamma][transition-1]/totalConvCoeff;
		}while(randomNumber<prob);
	}
			
	return transition;
}

// Return x ray energy after K conversion. Only use K x rays, if Ka1 or Ka2 produced then
// use generateLXRay function (see comment there).
G4int generateKXRay()
{
	// Get sum of x ray intensities
	double totalXRayIntensity = 0;
	for(int i=0; i<9; i++) totalXRayIntensity+=kIntensity[i];

	// Default to no x ray
	int transition=-1;
	
	// Check if transition is converted, if so choose CE for transition
	double randomNumber = G4UniformRand();
	double prob = 1.0;
	
	if(randomNumber<totalXRayIntensity/100)
	{
		randomNumber = G4UniformRand();
		do
		{
			transition++;
			prob -= kIntensity[transition]/totalXRayIntensity;
		}while(randomNumber<prob);
	}
	
	return transition;
}

// Generate x rays after L conversion, or following Ka1/Ka2 from K conversion
// Note that x ray intensities for direct L conversion from subshell are not
// consistent with those from K conversion tables (/100 K conversions), estimate
// ~15% uncertainty in L conversion intensity
G4int generateLXRay(int subShell)
{
	// Get sum of x ray intensities
	double totalXRayIntensity = 0;
	int arraySize = lIntensity[subShell].size();
	for(int i=0; i<arraySize; i++) totalXRayIntensity+=lIntensity[subShell][i];

	// Default to no x ray
	int transition=-1;
	
	// Check if transition is converted, if so choose CE for transition
	double randomNumber = G4UniformRand();
	double prob = 1.0;
	
	if(randomNumber<totalXRayIntensity/100)
	{
		randomNumber = G4UniformRand();
		do
		{
			transition++;
			prob -= lIntensity[subShell][transition]/totalXRayIntensity;
		}while(randomNumber<prob);
	}
	return transition;
}
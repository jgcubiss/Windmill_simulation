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
/// \file example.cc
/// \brief Main program of the  example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "myVariables.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " example [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
	// Read in variables from files
	readConfigFile();
	readGammaFile();
	readAlphaFile();
	readXRayFiles();
	
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
#ifdef G4MULTITHREADED
  G4int nThreads = 4;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Optionally: choose a different Random engine...
  //
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the default run manager
  //
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
#ifdef G4MULTITHREADED
  if ( nThreads > 0 ) { 
    runManager->SetNumberOfThreads(nThreads);
  }  
#endif

  // Set mandatory initialization classes
  //
  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  auto physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);
    
  auto actionInitialization = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);
  
  // Initialize visualization
  auto visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();
	
  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Read input from config file and output info to terminal
void readConfigFile()
{
	std::string dummy, dummy2;
	std::ifstream configFile;
	configFile.open("config.cfg");
	std::vector<double> var;
	G4cout << "\n\n";
	while(configFile >> dummy)
	{
		if(dummy=="#") getline(configFile,dummy);
		else var.push_back(std::stod(dummy));
	}

	Z = (int)var[0];
	annularPositionZ = var[1];
	circularPositionZ = var[2];
	optionAlpha = (int)var[3];
	optionGamma = (int)var[4];
	optionCE = (int)var[5];
	optionX = (int)var[6];

	var.clear();
	configFile.close();
	
	std::string alphas, gammas, ce, xrays;
	if(optionAlpha==1) alphas = "ON";
	if(optionAlpha==0) alphas = "OFF";
	if(optionGamma==1) gammas = "ON";
	if(optionGamma==0) gammas = "OFF";
	if(optionCE==1)        ce = "ON";
	if(optionCE==0)        ce = "OFF";
	if(optionX==1) xrays  = "ON";
	if(optionX==0) xrays  = "OFF";
	
	G4cout << "\n  ************************************************\n";
	G4cout << "  Z of daughter     = " << Z << "\n"
	       << "  Annular position  = " << annularPositionZ  << " mm\n"
		   << "  Circular position = " << circularPositionZ << " mm\n\n"
		   << "  Alphas   " << alphas << "\n"
		   << "  Gammas   " << gammas << "\n"
		   << "  CEs      " << ce     << "\n"
		   << "  X rays   " << xrays  << "\n"
		   << "  ************************************************\n\n";
}

// Read input from gammas file and store in matrices
void readGammaFile()
{
	std::ifstream  gammaFile;
	gammaFile.open("g_in.dat");
	double tmp;
	std::vector<double> var0, var1;
	int nColumns=0;
	while(gammaFile >> tmp)
	{
		if(nColumns==0) gammaEnergy.push_back(tmp);
		else if(nColumns%2!=0) var0.push_back(tmp);
		else var1.push_back(tmp);
		
		nColumns++;
		
		if(nColumns==17)
		{
			ceEnergy.push_back(var0);
			ceIntensity.push_back(var1);
			var0.clear();
			var1.clear();
			nColumns=0;
		}
	}
	gammaFile.close();
}

// Read input from alphas file
void readAlphaFile()
{
	std::ifstream  alphaFile;
	alphaFile.open("a_in.dat");
	double tmp;
	std::vector<double> var;
	nGammas = gammaEnergy.size();
	nAlphas = 0;
	int nColumns = 0;
	while(alphaFile >> tmp)
	{
		if(nColumns==0) alphaEnergy.push_back(tmp);
		else if(nColumns==1) alphaIntensity.push_back(tmp);
		else var.push_back(tmp);
		nColumns++;
		if(nColumns==nGammas+2)
		{
			nColumns = 0;
			gammaIntensity.push_back(var);
			var.clear();
			nAlphas++;
		}
	}
	alphaFile.close();
}

// Read input from x ray files
void readXRayFiles()
{
	// Read in x rays following K conversion for daughter nucleus
	std::ifstream  kXRayFile;
	std::string dummy;
	double tmpEnergy, tmpIntensity;
	kXRayFile.open("../include/K_conversion_xrays.dat");
	// Skip lines to daughter z
	for(int i=0; i<Z-5; i++) getline(kXRayFile,dummy);
	
	for(int i=0; i<24; i++)
	{
		kXRayFile >> tmpEnergy >> tmpIntensity;
		kEnergy.push_back(tmpEnergy);
		kIntensity.push_back(tmpIntensity);
	}
	kXRayFile.close();
	
	// Read in x rays following L1, L2 and L3 conversion for daughter nucleus
	int firstZ[3]   = {13, 13, 21};
	int nColumns[3] = {14, 10, 6};

	for(int i=0; i<3; i++)
	{
		std::ostringstream fileName;
		fileName << "../include/L" << i+1 << "_conversion_xrays.dat";
		std::ifstream  LXRayFile;
		LXRayFile.open(fileName.str());
		
		std::vector<double> tmpE, tmpI;
		
		// Skip lines to daughter Z
		for(int j=0; j<Z-firstZ[i]; j++) getline(LXRayFile,dummy);
				
		// Read in energies and intensities
		for(int j=0; j<nColumns[i]; j++)
		{
			LXRayFile >> tmpEnergy >> tmpIntensity;
			tmpE.push_back(tmpEnergy);
			tmpI.push_back(tmpIntensity);
		}
		
		// Push to matrix
		lEnergy.push_back(tmpE);
		lIntensity.push_back(tmpI);
	
		tmpE.clear();
		tmpI.clear();
		LXRayFile.close();
	}
}

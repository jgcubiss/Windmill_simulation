#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TError.h>

#include <Riostream.h>

double funcCB(double *x, double *par);

void broaden()
{
	// Parameters taken from 199At data
	Double_t parAlpha[2] = {1.75847, 1.5280};
	Double_t     parN[2] = {1.6587,  1.97958};	
	Double_t parSigma[2] = {11.285,  11.6651};	

	Double_t eShift[2] = {0.0, 0.0};

	ostringstream ipFile, opFile;
	ipFile << "Raw.root";
	opFile << "Broad.root";

	// Open input and output files
	TFile *iFile = TFile::Open(ipFile.str().c_str(),"UPDATE");
	TFile *oFile = TFile::Open(opFile.str().c_str(),"RECREATE");

	// Loads tree from in file
	TTree* iTree = (TTree*)iFile->Get("Raw");
	oFile->cd();
	iTree->CloneTree()->Write();
	iFile->Write();

	Int_t nEntries = iTree->GetEntries();
	
	// Define broadened E variables
	Double_t E_circ, L_circ, E_ann, L_ann;
	Int_t multc, multa, multt, decPath;

	// Define output variables
	Double_t E_c, E_c_b, L_c, E_a, E_a_b, L_a;
	Int_t multc_, multa_, multt_, decPath_;

	// Clone tree from in file
	iTree->SetBranchAddress("E_circ", &E_circ);
	iTree->SetBranchAddress("L_circ", &L_circ);
	iTree->SetBranchAddress("E_ann",  &E_ann);
	iTree->SetBranchAddress("L_ann",  &L_ann);
	iTree->SetBranchAddress("multC",  &multc);
	iTree->SetBranchAddress("multA",  &multa);
	iTree->SetBranchAddress("multT",  &multt);
	iTree->SetBranchAddress("path",   &decPath);

	// Define branches for out file
	TTree *oTree = new TTree("Broad","Broad");
	oTree->Branch("E_circ", &E_c_b,    "E_circ_b/D");
	oTree->Branch("L_circ", &L_c,      "L_circ/D");
	oTree->Branch("E_ann" , &E_a_b,    "E_ann_b/D");
	oTree->Branch("L_ann" , &L_a,      "L_ann/D");
	oTree->Branch("multC" , &multc_,   "multC/I");
	oTree->Branch("multA" , &multa_,   "multA/I");
	oTree->Branch("multT" , &multt_,   "multT/I");
	oTree->Branch("path"  , &decPath_, "path/I");

	TRandom3 r;

	TF1 *f[2];

	for(int i=0; i<2; i++)
	{
		ostringstream oss, fs;
		oss << "f[" << i << "]";
		fs << "ROOT::Math::crystalball_function(x," << parAlpha[i] << "," << parN[i] << "," << parSigma[i] << "," << 0 << ")";

		f[i] = new TF1(oss.str().c_str(),fs.str().c_str(),-10000,10000);
		f[i]->SetNpx(10000);
	}

	cout << "\n  -------------------Broadening simulated data---------------------\n\n";

	// Give each entry new param, add to new branch
	for(int i=0; i<nEntries; i++)
	{
		if(i % 1000 == 0) cout << "  Processed events " << i << "/"<< nEntries << "\r";
		iTree->GetEntry(i);
		E_c_b = 0;
		E_a_b = 0;

		if(E_ann>0)  E_a_b = f[0]->GetRandom()+E_ann+eShift[0];
		if(E_circ>0) E_c_b = f[1]->GetRandom()+E_circ+eShift[1];

		L_a      = L_ann;
		L_c      = L_circ;
		multc_   = multc;
		multa_   = multa;
		multt_   = multt;
		decPath_ = decPath;
		oTree->Fill();
	}
	cout << "  Processed events " << nEntries << "/"<< nEntries << "\r";
	cout << "\n\n  ----------------------Processed All Events-----------------------\n" << endl;

	iTree->Delete();
	oTree->Write();
	oFile->Write();
	gROOT->ProcessLine(".q");
}

double funcCB(double *x, double *par)
{
	double  xcur = x[0];
	double alpha = par[0];
	double     n = par[1];
	double sigma = par[2];
	double    mu = par[3];

	double f = ROOT::Math::crystalball_function(xcur,alpha,n,sigma,mu);

	return f;
}

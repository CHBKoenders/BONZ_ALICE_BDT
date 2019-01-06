#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"
#include "TList.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodBase.h"

void Select(float pTMin=0, float pTMax=90)
{	
	TFile *inputFile = TFile::Open("FullSimulation_LHC17pCENTwSDD_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");
	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	
	//set branches
	float pt_cand;
	float invM;
	bool sel;
	
	inputTree->SetBranchAddress("pt_cand", &pt_cand);	
	inputTree->SetBranchAddress("inv_mass", &invM);
	inputTree->SetBranchAddress("isselectedstd", &sel);

	//open file for histograms
	TFile *histFile = TFile::Open("ALICE_regular_selection_pT_"+pTRange+".root", "recreate");
	
	TH1D *selectedHist = new TH1D("selectedHist", "Invariant mass distribution of selected candidates", 500, 1, 2.5);
	TH1D *rejectedHist = new TH1D("rejectedHist", "Invariant mass distribution of rejected candidates", 500, 1, 2.5);
	
	int nEntries = inputTree->GetEntries();
	
	for(int i=1; i<nEntries; i++)
	{
		inputTree->GetEntry(i);
		if(not(pt_cand >= pTMin && pt_cand <= pTMax))
		{
			//outside of pTRange so do nothing
		}
		else if(sel)
		{
			selectedHist->Fill(invM);
		}
		else
		{
			rejectedHist->Fill(invM);
		}
	}
	
	histFile->Write();
	inputFile->Close();
}

void SelectTestingData(float pTMin=2, float pTMax = 90)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	
	TFile *fInput = new TFile(fName);
	
	//read Tree from file
	TTree *testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	//set branches
	float BDT_value, invM;
	bool isSig, isselectedstd;
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_value);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	testingOutputTree->SetBranchAddress("isselectedst", &isselectedstd);
	
	
}

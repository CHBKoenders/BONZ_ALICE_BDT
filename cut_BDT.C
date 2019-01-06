#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/TMVAGui.h"

#include <vector>
#include <iostream>

void Cut(float cutVal = 0, float pTMin = 0, float pTMax = 90, bool usepTBinTrained = 1)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString cutValue = Form("%.4f", cutVal);
	
	TString histFileName;
	TString histFileDirName = "BDT_results/";
	if(usepTBinTrained)
	{
		histFileName = "BDT_result_pT_"+pTRange+"_cut_at_"+cutValue+".root";
	}
	else
	{
		histFileName = "BDT_result_pT_"+pTRange+"_cut_at_"+cutValue+"_no_pT_train.root";
	}
	
	TFile *histFile = TFile::Open(histFileDirName+histFileName, "RECREATE");
	
	TH1D *sigHist = new TH1D("sigHist", "Classified Signal when cutting at "+cutValue, 100, 1.6,2.2);
	TH1D *bkgHist = new TH1D("bkgHist", "Classified Background when cutting at "+cutValue, 100, 1.6,2.2);
	
	//Open file
	TString inFileName;
	if(usepTBinTrained)
	{
		inFileName = "response_pT_"+pTRange+".root";
	}
	else
	{
		inFileName = "response_pT_"+pTRange+"_no_pT_train.root";
	}
	
	if(gSystem->AccessPathName(inFileName))
	{
		std::cout << "The file " << inFileName << " does not exist! Fixing this by running apply_BDT.C" << std::endl;
		gROOT->ProcessLine(".L apply_BDT.C");
		gROOT->ProcessLine(Form("Apply(%.2f, %.2f,%d)",pTMin, pTMax,usepTBinTrained));
	}
	
	TFile *inFile = TFile::Open(inFileName);
	
	TTree *inTree = (TTree*)inFile->Get("responseTree");
	
	float BDT_Response;
	float invM;
	
	inTree->SetBranchAddress("BDT_Response", &BDT_Response);
	inTree->SetBranchAddress("invM", &invM);
	
	int numEvents = inTree->GetEntries();
	
	for(int i=0; i<numEvents; i++)
	{
		inTree->GetEntry(i);
		
		if(BDT_Response >= cutVal)
		{
			//classify as signal
			sigHist->Fill(invM);
		
		}
		else
		{
			//classify as background
			bkgHist->Fill(invM);
		}
	}
	
	histFile->Write("", TObject::kOverwrite);
	
	histFile->Close();
	inFile->Close();
}

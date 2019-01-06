	#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodBase.h"
#include "TMVA/ResultsClassification.h"

#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "TVectorF.h"

#include <vector>
#include <iostream>

TTree *testingOutputTree = NULL;
float BDT_value;
float invM;
float pt_prong0, pt_prong1, pt_prong2;
float norm_dl_xy;
bool isSig;
bool isselectedstd;

double gSigEff, gSigEffErr;
double gPreCutEff;
double gCutValue;
double gCorrYield;
double gCorrYieldErr;

bool badCut;

float massDPlus = 1.86962; //In Gev

void GetPreCutEff(float pTMin, float pTMax)
{
	TCut preCut = "";
	
	TCut pTCut = Form("pt_cand >= %f && pt_cand <= %f",pTMin,pTMax);
	TCut sigCut = "cand_type!=0";
	
	//open file
	TFile *inputFile = TFile::Open("FullSimulation_LHC18a4a2CENT_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");

	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	int tot_sig = inputTree->GetEntries(sigCut&&pTCut);
	int pass_sig = inputTree->GetEntries(sigCut&&pTCut&&preCut);
	
	gPreCutEff = (double)pass_sig/tot_sig;
	
	//delete inputList;
	//TODO why must I delete inputTree if I close its file?
	//delete inputTree;
	//inputFile->Close();
	delete inputList;
	delete inputFile;
}

void GetSigEffBDT(float cutValue)
{		
	int numEvents = testingOutputTree -> GetEntries();
	
	int true_sig = 0;
	int false_sig = 0;
	int true_bkg = 0;
	int false_bkg = 0;
	int sigCount = 0;
	
	for (int iEvent = 0; iEvent < numEvents; iEvent++) 
	{
	 	testingOutputTree->GetEntry(iEvent);

	 	if(BDT_value >= cutValue)
	 	{
	 		sigCount++;
	 		//classified as signal
	 		true_sig += isSig; 	//+1 if it's actually signal
	 			 	
	 		false_sig += not isSig; //+1 if it's acutally background
	
	 	}
	 	else
	 	{
			//classified as background
	 		false_bkg += isSig; 	//+1 if it's actually signal
		 		
	 		true_bkg += not isSig; //+1 if it's acutally background
	 		
	 	}
	}
	
	if(sigCount < 100)
	{
		std::cout << "Less than a 100 events pass at this cut value!" << std::endl;
		badCut = 1;
		return;
	}
	
	badCut = 0;
	
	float sigEff = (float)true_sig/(true_sig + false_bkg);
	float bkgEff = (float)false_sig/(true_bkg + false_sig);
	float bkgRej = (float)true_bkg/(true_bkg + false_sig);
	
	gSigEff = sigEff;
	gSigEffErr = 1/TMath::Sqrt(numEvents);
}

void GetSigEffStd()
{	
	int numEvents = testingOutputTree -> GetEntries();
	
	int true_sig = 0;
	int false_sig = 0;
	int true_bkg = 0;
	int false_bkg = 0;
	
	for (int iEvent = 0; iEvent < numEvents; iEvent++) 
	{
	 	testingOutputTree->GetEntry(iEvent);
	 	
	 	if(isselectedstd)
	 	{
	 		//classified as signal
	 		true_sig += isSig; 	//+1 if it's actually signal
	 			 	
	 		false_sig += not isSig; //+1 if it's acutally background
	
	 	}
	 	else
	 	{
			//classified as background
	 		false_bkg += isSig; 	//+1 if it's actually signal
		 		
	 		true_bkg += not isSig; //+1 if it's acutally background
	 		
	 	}
	}
	
	float sigEff = (float)true_sig/(true_sig + false_bkg);
	float bkgEff = (float)false_sig/(true_bkg + false_sig);
	float bkgRej = (float)true_bkg/(true_bkg + false_sig);
	
	gSigEff = sigEff;
	gSigEffErr = 1/TMath::Sqrt(numEvents);

}

void CorrectBDTAtCutTesting(float pTMin = 2, float pTMax = 90, float cutValue = 0)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	
	//Load testing data
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	else
	{
		fInput = TFile::Open(fName);
	}
	
	//read Tree from file
	testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	//set global branches
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_value);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	testingOutputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	
	double sigEff, sigEffErr;
	double sigYield, sigYieldErr;
	double bkgYield, bkgYieldErr;
	double significance, significanceErr;
	
	TFile *fitFile;
	TString fileName;
	TString cutValueS = Form("%.4f", cutValue);

	//open BDT fit file
	fileName = "Fitting_results/BDT_fit_pT_"+pTRange+"_at_cut_"+cutValueS+"_testing.root";
	
	if(gSystem->AccessPathName(fileName))
	{
		std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
		
		return;
	}
	fitFile = TFile::Open(fileName);
	
	TVectorF *sigYieldV = (TVectorF*)fitFile->Get("sigYield");
	sigYield = (*sigYieldV)[0];
	sigYieldErr = (*sigYieldV)[1];
	TVectorF *bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
	bkgYield = (*bkgYieldV)[0];
	bkgYieldErr = (*bkgYieldV)[1];
	TVectorF *significanceV = (TVectorF*)fitFile->Get("significance");
	significance = (*significanceV)[0];
	significanceErr = (*significanceV)[1];
	
	TVectorF *cutValueV = (TVectorF*)fitFile->Get("cutValue");
	
	GetSigEffBDT((*cutValueV)[0]);
	GetPreCutEff(pTMin, pTMax);


	sigEff = gSigEff*gPreCutEff;
	sigEffErr = gSigEffErr;
	
	//std::cout << Form("Signal efficiency of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, sigEff, sigEffErr) << std::endl;
	
	gCorrYield = sigYield/sigEff;
	gCorrYieldErr = TMath::Abs(gCorrYield)*TMath::Sqrt((sigEffErr/sigEff)*(sigEffErr/sigEff)+(sigYieldErr/sigYield)*(sigYieldErr/sigYield));
	
	//std::cout << Form("Efficiency-corrected signal yield of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, corrYield, corrYieldErr) << std::endl;
	
	delete fitFile;
	delete fInput;
	testingOutputTree = NULL;
}

void CorrectBDTAtCut(float pTMin = 2, float pTMax = 90, float cutValue = 0)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	
	//Load testing data
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	else
	{
		fInput = TFile::Open(fName);
	}
	

	if(!testingOutputTree)
	{
	//read Tree from file
	testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	//set global branches
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_value);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	//testingOutputTree->SetBranchAddress("norm_dl_xy", &norm_dl_xy);
	//testingOutputTree->SetBranchAddress("pt_prong0", &pt_prong0);
	//testingOutputTree->SetBranchAddress("pt_prong1", &pt_prong1);
	//testingOutputTree->SetBranchAddress("pt_prong2", &pt_prong2);
	
	testingOutputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	}

	double sigEff, sigEffErr;
	double sigYield, sigYieldErr;
	double bkgYield, bkgYieldErr;
	double significance, significanceErr;
	
	TFile *fitFile;
	TString fileName;
	TString cutValueS = Form("%.4f", cutValue);

	//open BDT fit file
	fileName = "Fitting_results/BDT_fit_pT_"+pTRange+"_at_cut_"+cutValueS+".root";
	
	if(gSystem->AccessPathName(fileName))
	{
		std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
		
		return;
	}
	fitFile = TFile::Open(fileName);
	
	TVectorF *sigYieldV = (TVectorF*)fitFile->Get("sigYield");
	sigYield = (*sigYieldV)[0];
	sigYieldErr = (*sigYieldV)[1];
	TVectorF *bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
	bkgYield = (*bkgYieldV)[0];
	bkgYieldErr = (*bkgYieldV)[1];
	TVectorF *significanceV = (TVectorF*)fitFile->Get("significance");
	significance = (*significanceV)[0];
	significanceErr = (*significanceV)[1];
	
	TVectorF *cutValueV = (TVectorF*)fitFile->Get("cutValue");
	
	GetSigEffBDT((*cutValueV)[0]);
	GetPreCutEff(pTMin, pTMax);


	sigEff = gSigEff*gPreCutEff;
	sigEffErr = gSigEffErr;
	
	//std::cout << Form("Signal efficiency of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, sigEff, sigEffErr) << std::endl;
	
	gCorrYield = sigYield/sigEff;
	gCorrYieldErr = TMath::Abs(gCorrYield)*TMath::Sqrt((sigEffErr/sigEff)*(sigEffErr/sigEff)+(sigYieldErr/sigYield)*(sigYieldErr/sigYield));
	
	//std::cout << Form("Efficiency-corrected signal yield of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, corrYield, corrYieldErr) << std::endl;

	
	delete testingOutputTree;
	testingOutputTree = NULL;
	fitFile->Close();
	fInput->Close();
}

void CorrectTesting(float pTMin = 2, float pTMax = 90, bool BDT = 1)
{
	//TODO check wether precuteff is needed
	
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	
	//Load testing data
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	else
	{
		fInput = TFile::Open(fName);
	}
	
	//read Tree from file
	testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	//set global branches
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_value);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	testingOutputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	
	double sigEff, sigEffErr;
	double sigYield, sigYieldErr;
	double bkgYield, bkgYieldErr;
	double significance, significanceErr;
	
	TFile *fitFile;
	TString fileName;
	TString typeS;
	if(BDT)
	{
		typeS = "BDT";
		//open BDT fit file
		fileName = "Fitting_results/BDT_fit_pT_"+pTRange+"_testing.root";
		
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		fitFile = TFile::Open(fileName);
		
		TVectorF *sigYieldV = (TVectorF*)fitFile->Get("sigYield");
		sigYield = (*sigYieldV)[0];
		sigYieldErr = (*sigYieldV)[1];
		TVectorF *bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
		bkgYield = (*bkgYieldV)[0];
		bkgYieldErr = (*bkgYieldV)[1];
		TVectorF *significanceV = (TVectorF*)fitFile->Get("significance");
		significance = (*significanceV)[0];
		significanceErr = (*significanceV)[1];
		
		TVectorF *cutValueV = (TVectorF*)fitFile->Get("cutValue");
		gCutValue = (*cutValueV)[0];
		std::cout << "gCutValue is set at " << gCutValue << std::endl;
		GetSigEffBDT(gCutValue);
	}
	else
	{
		typeS = "std";
		//open Std fit file
		fileName = "Fitting_results/Std_fit_pT_"+pTRange+"_testing.root";
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		
		fitFile = TFile::Open(fileName);
		
		TVectorF *sigYieldV = (TVectorF*)fitFile->Get("sigYield");
		sigYield = (*sigYieldV)[0];
		sigYieldErr = (*sigYieldV)[1];
		TVectorF *bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
		bkgYield = (*bkgYieldV)[0];
		bkgYieldErr = (*bkgYieldV)[1];
		TVectorF *significanceV = (TVectorF*)fitFile->Get("significance");
		significance = (*significanceV)[0];
		significanceErr = (*significanceV)[1];
		
		GetSigEffStd();
	}

	sigEff = gSigEff;
	sigEffErr = gSigEffErr;
	
	//std::cout << Form("Signal efficiency of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, sigEff, sigEffErr) << std::endl;
	
	gCorrYield = sigYield/sigEff;
	gCorrYieldErr = TMath::Abs(gCorrYield)*TMath::Sqrt((sigEffErr/sigEff)*(sigEffErr/sigEff)+(sigYieldErr/sigYield)*(sigYieldErr/sigYield));
	
	//std::cout << Form("Efficiency-corrected signal yield of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, corrYield, corrYieldErr) << std::endl;
	
	delete fitFile;
	delete fInput;
	testingOutputTree = NULL;
}

int GetTrueYield(float pTMin, float pTMax)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";	
	
	//Load testing data
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return 0;
	
	}
	else
	{
		fInput = TFile::Open(fName);
	}
	
	//read Tree from file
	testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	//set global branches
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	
	int nEntries = testingOutputTree->GetEntries();
	int yield = 0;
	for(int i = 0; i < nEntries; i++)
	{
		testingOutputTree->GetEntry(i);
		
		yield += isSig;
		
	}
	
	
	delete fInput;
	testingOutputTree = NULL;
	
	return yield;
}

void Correct(float pTMin = 2, float pTMax = 90, bool BDT = 1)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	
	//Load testing data
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	else
	{
		fInput = TFile::Open(fName);
	}
	
	if(!testingOutputTree)
	{
	//read Tree from file
	testingOutputTree = (TTree*)fInput->Get("testOutputTree");

	
	//set global branches
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_value);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	//testingOutputTree->SetBranchAddress("norm_dl_xy", &norm_dl_xy);
	//testingOutputTree->SetBranchAddress("pt_prong0", &pt_prong0);
	//testingOutputTree->SetBranchAddress("pt_prong1", &pt_prong1);
	//testingOutputTree->SetBranchAddress("pt_prong2", &pt_prong2);
	testingOutputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	}
	
	double sigEff, sigEffErr;
	double sigYield, sigYieldErr;
	double bkgYield, bkgYieldErr;
	double significance, significanceErr;
	
	TFile *fitFile;
	TString fileName;
	TString typeS;
	if(BDT)
	{
		typeS = "BDT";
		//open BDT fit file
		fileName = "Fitting_results/BDT_fit_pT_"+pTRange+".root";
		
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		fitFile = TFile::Open(fileName);
		
		TVectorF *sigYieldV = (TVectorF*)fitFile->Get("sigYield");
		sigYield = (*sigYieldV)[0];
		sigYieldErr = (*sigYieldV)[1];
		TVectorF *bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
		bkgYield = (*bkgYieldV)[0];
		bkgYieldErr = (*bkgYieldV)[1];
		TVectorF *significanceV = (TVectorF*)fitFile->Get("significance");
		significance = (*significanceV)[0];
		significanceErr = (*significanceV)[1];
		
		TVectorF *cutValueV = (TVectorF*)fitFile->Get("cutValue");
		gCutValue = (*cutValueV)[0];
		std::cout << "gCutValue is set at " << gCutValue << std::endl;
		GetSigEffBDT(gCutValue);
	}
	else
	{
		typeS = "std";
		//open Std fit file
		fileName = "Fitting_results/Std_fit_pT_"+pTRange+".root";
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		
		fitFile = TFile::Open(fileName);
		
		TVectorF *sigYieldV = (TVectorF*)fitFile->Get("sigYield");
		sigYield = (*sigYieldV)[0];
		sigYieldErr = (*sigYieldV)[1];
		TVectorF *bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
		bkgYield = (*bkgYieldV)[0];
		bkgYieldErr = (*bkgYieldV)[1];
		TVectorF *significanceV = (TVectorF*)fitFile->Get("significance");
		significance = (*significanceV)[0];
		significanceErr = (*significanceV)[1];
		
		GetSigEffStd();
	}
	std::cout << "getting pre cut eff" << std::endl;
	GetPreCutEff(pTMin, pTMax);

	sigEff = gSigEff*gPreCutEff;
	sigEffErr = gSigEffErr;
	
	std::cout << Form("Signal efficiency of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, sigEff, sigEffErr) << std::endl;
	
	gCorrYield = sigYield/sigEff;
	gCorrYieldErr = TMath::Abs(gCorrYield)*TMath::Sqrt((sigEffErr/sigEff)*(sigEffErr/sigEff)+(sigYieldErr/sigYield)*(sigYieldErr/sigYield));
	
	std::cout << Form("Efficiency-corrected signal yield of "+typeS+" analysis in pTRange from %.2f to %.2f is %.4f +- %.4f",pTMin, pTMax, gCorrYield, gCorrYieldErr) << std::endl;
	
	delete testingOutputTree;
	testingOutputTree = NULL;
	delete fitFile;
	delete fInput;
}

void FillSigEffGraph(float horScale = 1, float vertScale = 0.5)
{
	gStyle->SetOptTitle(0);
	
	TGraphErrors *effGraphBDT = new TGraphErrors();
	TGraphErrors *effGraphStd = new TGraphErrors();
	
	int pTLowArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int numBins = 10;
	
	for(int i=0; i<numBins; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		//BDT
		Correct(pTLow, pTHigh, 1);
		effGraphBDT->SetPoint(i, (pTLow+pTHigh)/2., gSigEff);
		effGraphBDT->SetPointError(i, (pTHigh-pTLow)/2., gSigEffErr);
		
		//Std
		Correct(pTLow, pTHigh, 0);
		effGraphStd->SetPoint(i, (pTLow+pTHigh)/2., gSigEff);
		effGraphStd->SetPointError(i, (pTHigh-pTLow)/2., gSigEffErr);
	}
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	TCanvas *c1 = new TCanvas;
	c1->SetCanvasSize(horScale*horSize, vertScale*vertSize);
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(effGraphBDT);
	combinedGraph->Add(effGraphStd);
	
	combinedGraph->SetTitle("Signal efficiency in #it{p}_{T} intervals");
	
	effGraphBDT->SetTitle("Boosted Decision Tree");
	effGraphBDT->SetLineColor(2);
	effGraphBDT->SetMarkerStyle(8);
	
	effGraphStd->SetTitle("ALICE Standard");
	effGraphStd->SetLineColor(4);
	effGraphStd->SetMarkerStyle(22);
	
	combinedGraph->Draw("AP");
	
	combinedGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Signal Efficiency");
	combinedGraph->GetYaxis()->CenterTitle();
	
	gPad->BuildLegend();

}

void FillCorrYieldAtCutGraph(float pTMin, float pTMax, int nCuts, bool testing = 0)
{
	TGraphErrors *corrYieldGraphBDT = new TGraphErrors();
	TGraphErrors *corrYieldGraphStd = new TGraphErrors();
	
	//to get values at opt cut
	std::cout << "getting opt cut values" << std::endl;
	if(testing)
	{
		CorrectTesting(pTMin,pTMax,1);
	}
	else
	{
		Correct(pTMin, pTMax, 1);
	}
	
	std::cout << "opt cut value is: " << gCutValue << std::endl;
	
	double optCutValue = gCutValue;
	float optCutYield = gCorrYield;
	float optCutYieldErr = gCorrYieldErr;
	
	//to get Std value
	std::cout<< "getting std values" << std::endl;
	if(testing)
	{
		CorrectTesting(pTMin, pTMax, 0);
	}
	else
	{
		Correct(pTMin, pTMax, 0);
	}
	
	corrYieldGraphStd->SetPoint(0, -1, gCorrYield);
	corrYieldGraphStd->SetPointError(0, 0, gCorrYieldErr);
	
	corrYieldGraphStd->SetPoint(1, 1, gCorrYield);
	corrYieldGraphStd->SetPointError(1, 0, gCorrYieldErr);
	
	float cutStep = 2./(nCuts+1);
	
	float cutValue;
	
	int offset = 0;
	
	float progress = 0.;
	
	std::cout << "starting loop..." << std::endl;
	for(int i=0; i<nCuts; i++)
	{
	
		progress = (float)100*(i+1)/nCuts;
		std::cout << int(progress) << "% done\r";
		std::cout.flush();
	
		cutValue = -1 + (i+1)*cutStep;
		
		if(testing)
		{
			CorrectBDTAtCutTesting(pTMin, pTMax, cutValue);
		}
		else
		{
			CorrectBDTAtCut(pTMin, pTMax, cutValue);
		}
		
		if(!badCut)
		{
			corrYieldGraphBDT->SetPoint(i+offset, cutValue, gCorrYield);
			corrYieldGraphBDT->SetPointError(i+offset, 0, gCorrYieldErr);
		}
		
		if(cutValue <= optCutValue && (cutValue+cutStep) >= optCutValue)
		{
			std::cout << "insertion point reached" << std::endl;
			corrYieldGraphBDT->SetPoint(i+1, optCutValue, optCutYield);
			corrYieldGraphBDT->SetPointError(i+1, 0, optCutYieldErr);
			
			offset = 1;
		}
		
		//std::cout << Form("Corrected yield at cut value %.4f is %.2f", cutValue, gCorrYield) << std::endl;
	}

	TCanvas *c1= new TCanvas;
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(corrYieldGraphStd, "PL 3");
	combinedGraph->Add(corrYieldGraphBDT, "C");
	
	TString title;
	if(testing)
	{
		title = Form("Corrected signal yield at BDT cut values in #it{p}_{T} range from %.2f to %.2f GeV/#it{c} on MC testing data",pTMin, pTMax);
	}
	else
	{
		title = Form("Corrected signal yield at BDT cut values in #it{p}_{T} range from %.2f to %.2f GeV/#it{c}",pTMin, pTMax);
	}
	combinedGraph->SetTitle(title);
	
	corrYieldGraphBDT->SetTitle("BDT");
	corrYieldGraphBDT->SetLineColor(2);

	corrYieldGraphStd->SetTitle("Std");
	corrYieldGraphStd->SetLineColor(4);
	corrYieldGraphStd->SetFillStyle(3001);
	corrYieldGraphStd->SetFillColor(4);

	
	combinedGraph->Draw("AP");
	
	combinedGraph->GetXaxis()->SetTitle("BDT Cut Value");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Corrected Yield");
	combinedGraph->GetYaxis()->CenterTitle();
	
	c1->Update();
	
	TLine *cutLine = new TLine(optCutValue, c1->GetUymin(), optCutValue, c1->GetUymax());
	
	cutLine->SetLineColor(6);
	cutLine->Draw();
	
	TLegend *legend = new TLegend(0.3,0.15);
	legend->AddEntry(corrYieldGraphBDT, "BDT","LF");
	legend->AddEntry(corrYieldGraphStd, "Std","LF");
	legend->AddEntry(cutLine, Form("Optimal significance cut (%.4f)",optCutValue), "L");	
	
	legend->Draw();
}
void FillCorrYieldGraph(bool testing = 0, float horScale = 1, float vertScale = .5)
{
	gStyle->SetOptTitle(0);
	
	TGraphErrors *corrYieldGraphBDT = new TGraphErrors();
	TGraphErrors *corrYieldGraphStd = new TGraphErrors();
	TGraphErrors *trueYieldGraph = new TGraphErrors();
	
	int pTLowArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int numBins = 10;
	
	for(int i=0; i<numBins; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		//BDT
		if(testing)
		{
			CorrectTesting(pTLow, pTHigh, 1);
		}
		else
		{
			Correct(pTLow, pTHigh, 1);
		}
		
		corrYieldGraphBDT->SetPoint(i, (pTLow+pTHigh)/2., gCorrYield);
		corrYieldGraphBDT->SetPointError(i, (pTHigh-pTLow)/2., gCorrYieldErr);
		
		//Std
		if(testing)
		{
			CorrectTesting(pTLow, pTHigh, 0);
		}
		else
		{
			Correct(pTLow, pTHigh, 0);
		}
		
		corrYieldGraphStd->SetPoint(i, (pTLow+pTHigh)/2., gCorrYield);
		corrYieldGraphStd->SetPointError(i, (pTHigh-pTLow)/2., gCorrYieldErr);
		
		//true yield
		if(testing)
		{
			int trueYield = GetTrueYield(pTLow, pTHigh);
			
			std::cout << Form("The true yield in pT range from %.2d to %.2d is %d", pTLow, pTHigh, trueYield) << std::endl;
			
			trueYieldGraph->SetPoint(i, (pTLow+pTHigh)/2., trueYield);
			trueYieldGraph->SetPointError(i, (pTHigh-pTLow)/2.,0);
		}
	}
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	TCanvas *c = new TCanvas;
	c->SetCanvasSize(horScale*horSize, vertScale*vertSize);
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(corrYieldGraphBDT);
	combinedGraph->Add(corrYieldGraphStd);
	if(testing)
	{
		combinedGraph->Add(trueYieldGraph);
		trueYieldGraph->SetTitle("'True' Yield");
		trueYieldGraph->SetLineColor(3);
		trueYieldGraph->SetLineWidth(1);	
	}
	
	combinedGraph->SetTitle("Efficiency-corrected yield in #it{p}_{T} intervals");
	
	corrYieldGraphBDT->SetTitle("Boosted Decision Tree");
	corrYieldGraphBDT->SetLineColor(2);
	corrYieldGraphBDT->SetMarkerStyle(8);
	
	corrYieldGraphStd->SetTitle("ALICE Standard");
	corrYieldGraphStd->SetLineColor(4);
	corrYieldGraphStd->SetMarkerStyle(22);
	
	combinedGraph->Draw("AP");
	
	combinedGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Efficiency-Corrected Yield");
	combinedGraph->GetYaxis()->CenterTitle();
	
	gPad->BuildLegend();
	
	
}

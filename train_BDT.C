#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodBase.h"

TMVA::Factory *factory;

TString pTRange;
TString methodName;

float massDPlus = 1.86962; //In Gev

void SaveTestingOutput()
{	
	TMVA::IMethod * imethod = factory->GetMethod("dataset",methodName);
	TMVA::MethodBase * method = dynamic_cast<TMVA::MethodBase *>(imethod);

	TMVA::DataSet *dataset = method->Data();
	TMVA::DataSetInfo &dsi = method->DataInfo();
	TMVA::Results *results = dataset->GetResults(methodName, TMVA::Types::kTesting, TMVA::Types::kClassification);
	
	//create file for our Testing Output Tree
	TString outName = "testing_output_"+pTRange+".root";
	TFile *fTestingOutputTree = new TFile(outName, "recreate");
	fTestingOutputTree->cd();
	
	//create TTree
	TTree *testingOutputTree = new TTree("testOutputTree", "BDT output on testing data");
	
	//set branches
	float pt_cand, d_len, d_len_xy, norm_dl_xy, cos_p, cos_p_xy, imp_par_xy, pt_prong0, pt_prong1, pt_prong2, sig_vert;
	
	float BDT_value;
	float invM;
	bool isSig,isselectedstd;
	
	testingOutputTree->Branch("pt_cand", &pt_cand);
	testingOutputTree->Branch("d_len", &d_len);
	testingOutputTree->Branch("d_len_xy", &d_len_xy);
	testingOutputTree->Branch("norm_dl_xy", &norm_dl_xy);
	testingOutputTree->Branch("cos_p", &cos_p);
	testingOutputTree->Branch("cos_p_xy", &cos_p_xy);
	testingOutputTree->Branch("imp_par_xy", &imp_par_xy);
	testingOutputTree->Branch("pt_prong0", &pt_prong0);
	testingOutputTree->Branch("pt_prong1", &pt_prong1);
	testingOutputTree->Branch("pt_prong2", &pt_prong2);
	testingOutputTree->Branch("sig_vert", &sig_vert);
	
	testingOutputTree->Branch("invM", &invM);
	testingOutputTree->Branch("isselectedstd", &isselectedstd);
	
	testingOutputTree->Branch("isSig", &isSig);
	testingOutputTree->Branch("BDT_value", &BDT_value);

	
	//fill the tree
	int numEvents = dataset->GetNTestEvents();
	std::cout << "Number of test events is " << numEvents << std::endl;
	int nSig=0;
	int nBkg=0;
	int count=0;	
	for (int iEvent = 0; iEvent < numEvents; ++iEvent) 
	{
	
		float valueEvent = (*results)[iEvent][0];
		const TMVA::Event *ev = dataset->GetEvent(iEvent, TMVA::Types::kTesting);
	    	
	    	BDT_value = valueEvent;
		isSig = dsi.IsSignal(ev); //True if signal
		
		nSig+=isSig;
		nBkg+=!isSig;
		count++;
		
		//same order as when adding variables to the dataloader
		d_len = ev->GetValue(0);
		d_len_xy = ev->GetValue(1);
		norm_dl_xy = ev->GetValue(2);
		cos_p = ev->GetValue(3);
		cos_p_xy = ev->GetValue(4);
		imp_par_xy = ev->GetValue(5);
		pt_prong0 = ev->GetValue(6);
		pt_prong1 = ev->GetValue(7);
		pt_prong2 = ev->GetValue(8);
		sig_vert = ev->GetValue(9);
		pt_cand = ev->GetValue(10);
		
		//spectator vars
		invM = ev->GetSpectator(0);
		isselectedstd = ev->GetSpectator(1);
		
		testingOutputTree->Fill();
	}

	//save the TTree
	testingOutputTree->Write("", TObject::kOverwrite);
	
	//close the file
	fTestingOutputTree->Close();
	
	std::cout << "number of checked events: " << count << std::endl;
	std::cout << "number of signal testing samples: " << nSig << std::endl;
	std::cout << "number of background testing samples: " << nBkg << std::endl;
}

void Train(float PtMin=0, float PtMax=90)
{
	
	std::cout << Form("Training BDT for pT range from %.2f to %.2f", PtMin, PtMax) << std::endl;
	
	pTRange = Form("%.2f_%.2f",PtMin,PtMax);

	//pre selection cuts
	TCut preCut = ""; //Form("pt_prong0>=0.3&&pt_prong1>=0.3&&pt_prong2>=0.3 &&inv_mass<=1.98", massDPlus); // && cos_p<=0.98 && cos_p_xy<=0.93 && norm_dl_xy>=0&&norm_dl_xy<=90 && d_len>=0&&d_len<=0.4 && imp_par_xy>=0.03&&imp_par_xy<=0.14 && imp_par_xy>=0.03&&imp_par_xy<=0.14
	
	TMVA::Tools::Instance();
	
	//signal
	TFile *sigFile = TFile::Open("sig_from_MC_pT_"+pTRange+".root");
	TTree *sigTree = (TTree*)sigFile->Get("MC_sig");
	
	//background
	TFile *bkgSidebandFile = TFile::Open("bkg_from_sidebands_pT_"+pTRange+".root");
	TTree *bkgSidebandTree = (TTree*)bkgSidebandFile->Get("data_sidebands");
	
	//TFile *bkgMCFile = TFile::Open("bkg_from_MC_pT_"+pTRange+".root");
	//TTree *bkgMCTree = (TTree*)bkgMCFile->Get("MC_bkg");
	
	
	TString outFileName = "BDT_event_output_Pt_"+pTRange+".root";
	TFile *outputFile = TFile::Open(outFileName, "RECREATE");
	
	methodName = "BDT_Gradient_pT_"+pTRange;
	
	//TODO always check factory name	
	factory = new TMVA::Factory("Classification_bkg_from_sidebands", outputFile,"Color:DrawProgressBar:AnalysisType=Classification");

	TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

	
	//add the tree
	dataloader->AddSignalTree(sigTree, 1);	
	dataloader->AddBackgroundTree(bkgSidebandTree, 1);
	//dataloader->AddBackgroundTree(bkgMCTree, 1);

	//add variables
	
	dataloader->AddVariable("d_len");
	dataloader->AddVariable("d_len_xy");
	dataloader->AddVariable("norm_dl_xy");
	dataloader->AddVariable("cos_p");
	dataloader->AddVariable("cos_p_xy");
	dataloader->AddVariable("imp_par_xy");
	dataloader->AddVariable("pt_prong0");
	dataloader->AddVariable("pt_prong1");
	dataloader->AddVariable("pt_prong2");
	dataloader->AddVariable("sig_vert");
	dataloader->AddVariable("pt_cand");
	
	dataloader->AddSpectator("inv_mass");
	dataloader->AddSpectator("isselectedstd");
	
	//cut on pT
	TCut pTCutS = Form("pt_cand >= %f && pt_cand <= %f",PtMin,PtMax);
	TCut pTCutB = Form("pt_cand >= %f && pt_cand <= %f",PtMin,PtMax);
	
	//total cut
	TCut totCutS = pTCutS&&preCut;
	TCut totCutB = pTCutB&&preCut;

	//get values for amount of training samples
	int nS = sigTree->GetEntries(totCutS);
	int nB = bkgSidebandTree->GetEntries(totCutB);// + bkgMCTree->GetEntries(totCutB);
	
	TString nTrainSig = Form("%i",static_cast<int>(nS*0.90));
	TString nTrainBkg = Form("%i", static_cast<int>(nB*0.75));
	
	TString nTestSig = Form("%i",static_cast<int>(nS*0.10));
	TString nTestBkg = Form("%i", static_cast<int>(nB*0.25));
	
	TString prepareOptions = "nTrain_Signal="+nTrainSig+":nTrain_Background="+nTrainBkg+":nTest_Signal="+nTestSig+":nTest_Background="+nTestBkg+":NormMode=None:SplitMode=Random";
		
	//prepare
	dataloader->PrepareTrainingAndTestTree(totCutS, totCutB, prepareOptions);
	
	//book BDT
	factory->BookMethod(dataloader, TMVA::Types::kBDT, methodName,"NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=50:MaxDepth=2" );

	//Train, Test and Evaluate
	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();
		
	SaveTestingOutput();
	
	outputFile->Write("",TObject::kOverwrite);
	
	delete factory;
	delete dataloader;
	
	delete sigFile;
	delete bkgSidebandFile;
	delete outputFile;
}

void TrainAllPtBins()
{

	int pTLowArr[] = 	{2,3,4,5,6,6,7,8 ,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,8,10,12,12,16,24};
	
	for(int i=0; i<12; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		Train(pTLow, pTHigh);
	}	
	
	
}

#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/TMVAGui.h"

#include <vector>
#include <iostream>

void Apply(float pTMin = 0, float pTMax = 90, bool usepTBinTrained = 1)
{
	TH1D *responseHist;
	
	if(usepTBinTrained)
	{
		responseHist = new TH1D("responseHist", Form("Response of the non pT Bin trained  BDT in pT range from %.2f to %.2f", pTMin, pTMax), 100,-1,1);
	}
	else
	{
		responseHist = new TH1D("responseHist", Form("Response of the BDT in pT range from %.2f to %.2f", pTMin, pTMax), 100,-1,1);
	}
		
	TMVA::Tools::Instance();
	
	TMVA::Reader *reader = new TMVA::Reader();
	
	//directory stuff
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);	
	TString weightFileName;
	
	//TODO check this name before running
	weightFileName = "Classification_bkg_from_sidebands_BDT_Gradient_pT_"+pTRange+".weights.xml";

	TString dirName = "dataset/weights/";
	
	if(gSystem->AccessPathName(dirName+weightFileName))
	{
		std::cout << "the weight file "+weightFileName+" does not exist! Train the BDT first. Aborting..." << std::endl;
		return;
	}

	//variables
	float pt_cand, d_len, d_len_xy, norm_dl_xy, cos_p, cos_p_xy, imp_par_xy, pt_prong0, pt_prong1, pt_prong2, sig_vert;
	
	//spectator
	float invM;
	int  isselectedstd;
	
	reader->AddVariable("d_len", &d_len);	
	reader->AddVariable("d_len_xy", &d_len_xy);	
	reader->AddVariable("norm_dl_xy", &norm_dl_xy);
	reader->AddVariable("cos_p", &cos_p);
	reader->AddVariable("cos_p_xy", &cos_p_xy);
	reader->AddVariable("imp_par_xy", &imp_par_xy);
	reader->AddVariable("pt_prong0", &pt_prong0);
	reader->AddVariable("pt_prong1", &pt_prong1);
	reader->AddVariable("pt_prong2", &pt_prong2);
	reader->AddVariable("sig_vert", &sig_vert);
	reader->AddVariable("pt_cand", &pt_cand);
	reader->AddSpectator("inv_mass", &invM);
	reader->AddSpectator("isselectedstd", &isselectedstd);

	
	TString outFileName;
	if(usepTBinTrained)
	{
		outFileName = "response_pT_"+pTRange+".root";
	}
	else
	{
		outFileName = "response_pT_"+pTRange+"_no_pT_train.root";
	}
	
	TFile *treeFile = TFile::Open(outFileName, "RECREATE");
	
	//create the output TTree
	treeFile->cd();
	TTree *outputTree = new TTree("responseTree", Form("Tree with the response of the BDT in the pT range from %.2f to %.2f", pTMin, pTMax));
	
	float BDT_Response;

	outputTree->Branch("invM", &invM);
	outputTree->Branch("BDT_Response", &BDT_Response);
	
	outputTree->Branch("pt_cand", &pt_cand);
	outputTree->Branch("d_len", &d_len);
	outputTree->Branch("d_len_xy", &d_len_xy);
	outputTree->Branch("norm_dl_xy", &norm_dl_xy);
	outputTree->Branch("cos_p", &cos_p);
	outputTree->Branch("cos_p_xy", &cos_p_xy);
	outputTree->Branch("imp_par_xy", &imp_par_xy);
	outputTree->Branch("pt_prong0", &pt_prong0);
	outputTree->Branch("pt_prong1", &pt_prong1);
	outputTree->Branch("pt_prong2", &pt_prong2);
	outputTree->Branch("sig_vert", &sig_vert);
		
	//Book method	
	
	reader->BookMVA("BDT_Gradient",dirName+weightFileName);
	
	//Get the inputTree
	TFile *inputFile = TFile::Open("FullSimulation_LHC17pCENTwSDD_Dplus_AnalysisResults.root");
	inputFile->cd();
	gDirectory->cd("PWGHF_TreeCreator");


	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	//set branches
	
	inputTree->SetBranchAddress("inv_mass", &invM);
	
	inputTree->SetBranchAddress("pt_cand", &pt_cand);
	inputTree->SetBranchAddress("d_len", &d_len);
	inputTree->SetBranchAddress("d_len_xy", &d_len_xy);
	inputTree->SetBranchAddress("norm_dl_xy", &norm_dl_xy);
	inputTree->SetBranchAddress("cos_p", &cos_p);
	inputTree->SetBranchAddress("cos_p_xy", &cos_p_xy);
	inputTree->SetBranchAddress("imp_par_xy", &imp_par_xy);
	inputTree->SetBranchAddress("pt_prong0", &pt_prong0);
	inputTree->SetBranchAddress("pt_prong1", &pt_prong1);
	inputTree->SetBranchAddress("pt_prong2", &pt_prong2);
	inputTree->SetBranchAddress("sig_vert", &sig_vert);

	//loop through entries and save response in new TTree
	int nEntries = inputTree->GetEntries();
	int sigCount = 0;
	int bkgCount = 0;

	float progress =0;
	std::cout << "pushing data through BDT..." << std::endl;
	for(int i=0; i<nEntries; i++)
	{
		inputTree->GetEntry(i);
		
		//we don't bother looking at any pt_prong < 0.3
		if(pt_prong0>=0.3 && pt_prong1 >=0.3 && pt_prong2 >= 0.3)
			{
			BDT_Response = reader->EvaluateMVA("BDT_Gradient");
			
			//Only look at the response if pTCand is in range
			if(pt_cand >= pTMin && pt_cand <= pTMax)
			{
				outputTree->Fill();			
				//responseHist->Fill(BDT_Response);
			}
			
			progress = 100*(i+1)/nEntries;
			
			if(progress < 10)
			{
				std::cout << "  " << int(progress) << "% done\r";
				std::cout.flush();
			}
			else if(progress < 100)
			{
				std::cout << " " << int(progress) << "% done\r";
				std::cout.flush();
			}
			else
			{
				std::cout << int(progress) << "% done\r";
				std::cout.flush();
			}
		}	
	}
	std::cout << std::endl;

	inputFile->Close();
	
	treeFile->Write("",TObject::kOverwrite);
	
	//TCanvas *c1 = new TCanvas("c1","histogram canvas",1000,1000);
	
	//responseHist->Draw();
}

void ApplyMC(float pTMin, float pTMax)
{
	TH1D *responseHist = new TH1D("responseHist", Form("Response of the BDT in pT range from %.2f to %.2f", pTMin, pTMax), 100,-1,1);

TMVA::Tools::Instance();
	
	TMVA::Reader *reader = new TMVA::Reader();
	
	//directory stuff
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);	
	TString weightFileName;

	weightFileName = "LightEventDataClassification_BDT_Gradient_pT_"+pTRange+".weights.xml";

	TString dirName = "dataset/weights/";
	
	if(gSystem->AccessPathName(dirName+weightFileName))
	{
		std::cout << "the weight file does not exist! Train the BDT first. Aborting..." << std::endl;
		return;
	}

	//variables
	float pt_cand, d_len, d_len_xy, norm_dl_xy, cos_p, cos_p_xy, imp_par_xy, pt_prong0, pt_prong1, pt_prong2, sig_vert;
	
	//spectator
	float invM;
	int  isselectedstd;
	
	reader->AddVariable("d_len", &d_len);
	reader->AddVariable("d_len_xy", &d_len_xy);
	reader->AddVariable("norm_dl_xy", &norm_dl_xy);
	reader->AddVariable("cos_p", &cos_p);
	reader->AddVariable("cos_p_xy", &cos_p_xy);
	reader->AddVariable("imp_par_xy", &imp_par_xy);
	reader->AddVariable("pt_prong0", &pt_prong0);
	reader->AddVariable("pt_prong1", &pt_prong1);
	reader->AddVariable("pt_prong2", &pt_prong2);
	reader->AddSpectator("sig_vert", &sig_vert);
	reader->AddSpectator("pt_cand", &pt_cand);
	reader->AddSpectator("inv_mass", &invM);
	reader->AddSpectator("isselectedstd", &isselectedstd);

	
	TString outFileName = "response_pT_"+pTRange+"_MC.root";
	
	TFile *treeFile = TFile::Open(outFileName, "RECREATE");
	
	//create the output TTree
	treeFile->cd();
	TTree *outputTree = new TTree("responseTree", Form("Tree with the response of the BDT in the pT range from %.2f to %.2f", pTMin, pTMax));
	
	float BDT_Response;
	bool isSig;
	
	outputTree->Branch("invM", &invM);
	outputTree->Branch("BDT_Response", &BDT_Response);
	
	outputTree->Branch("pt_cand", &pt_cand);
	outputTree->Branch("d_len", &d_len);
	outputTree->Branch("d_len_xy", &d_len_xy);
	outputTree->Branch("norm_dl_xy", &norm_dl_xy);
	outputTree->Branch("cos_p", &cos_p);
	outputTree->Branch("cos_p_xy", &cos_p_xy);
	outputTree->Branch("imp_par_xy", &imp_par_xy);
	outputTree->Branch("pt_prong0", &pt_prong0);
	outputTree->Branch("pt_prong1", &pt_prong1);
	outputTree->Branch("pt_prong2", &pt_prong2);
	outputTree->Branch("sig_vert", &sig_vert);
	outputTree->Branch("isSig", &isSig);
		
	//Book method	
	
	reader->BookMVA("BDT_Gradient",dirName+weightFileName);
	
	//input tree
	TFile *inputFile = TFile::Open("FullSimulation_LHC18a4a2CENT_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");


	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	//set branches
	Char_t candType;
	
	inputTree->SetBranchAddress("inv_mass", &invM);
	
	inputTree->SetBranchAddress("pt_cand", &pt_cand);
	inputTree->SetBranchAddress("d_len", &d_len);
	inputTree->SetBranchAddress("d_len_xy", &d_len_xy);
	inputTree->SetBranchAddress("norm_dl_xy", &norm_dl_xy);
	inputTree->SetBranchAddress("cos_p", &cos_p);
	inputTree->SetBranchAddress("cos_p_xy", &cos_p_xy);
	inputTree->SetBranchAddress("imp_par_xy", &imp_par_xy);
	inputTree->SetBranchAddress("pt_prong0", &pt_prong0);
	inputTree->SetBranchAddress("pt_prong1", &pt_prong1);
	inputTree->SetBranchAddress("pt_prong2", &pt_prong2);
	inputTree->SetBranchAddress("sig_vert", &sig_vert);
	inputTree->SetBranchAddress("cand_type", &candType);

	//loop through entries and save response in new TTree
	int nEntries = inputTree->GetEntries();
	int sigCount = 0;
	int bkgCount = 0;

	float progress =0;
	std::cout << "pushing data through BDT..." << std::endl;
	for(int i=0; i<nEntries; i++)
	{
		inputTree->GetEntry(i);
		
			isSig = (candType != 0);
				
			BDT_Response = reader->EvaluateMVA("BDT_Gradient");
			
			//Only look at the response if pTCand is in range
			if(pt_cand >= pTMin && pt_cand <= pTMax)
			{
				outputTree->Fill();			
				//responseHist->Fill(BDT_Response);
			}
			
			progress = 100*(i+1)/nEntries;
			
			if(progress < 10)
			{
				std::cout << "  " << int(progress) << "% done\r";
				std::cout.flush();
			}
			else if(progress < 100)
			{
				std::cout << " " << int(progress) << "% done\r";
				std::cout.flush();
			}
			else
			{
				std::cout << int(progress) << "% done\r";
				std::cout.flush();
			}
			
	}
	
	std::cout << std::endl;

	inputFile->Close();
	
	treeFile->Write("",TObject::kOverwrite);
	
	//TCanvas *c1 = new TCanvas("c1","histogram canvas",1000,1000);

}

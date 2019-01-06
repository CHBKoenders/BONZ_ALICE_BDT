#include "TTree.h"
#include "TFile.h"

#include "TMath.h"

//sigma used is for the entire data set, nog a pT range -> next step
void CreateSidebandFile(float pTMin, float pTMax)
{
	std::cout << Form("Creating sideband training file for pT range from %.2f to %.2f", pTMin, pTMax) << std::endl;
	
	TString pTRange = Form("%.2f_%.2f",pTMin, pTMax);

	float sigma;
	float mean;
	
	//to get sigma and mean
	TString fSigmaName = Form("MC_true_signal_pT_%.2f_%.2f.root", pTMin,  pTMax);
	TString dirSigma = "MC_true_signal/";
	
	if(gSystem->AccessPathName(dirSigma+fSigmaName))
	{
		//file does not exist yet
		gROOT->ProcessLine(".L get_sigma.C");
		gROOT->ProcessLine(Form("GetSigma(%.2f,%.2f)",pTMin, pTMax));	
	}
	
	TFile *sigmaFile = TFile::Open(dirSigma+fSigmaName);
	TVectorF *vSigma = (TVectorF*)sigmaFile->Get("sigma");
	TVectorF *vMean = (TVectorF*)sigmaFile->Get("mean");
	sigma = TMath::Abs((*vSigma)[0]);
	mean = (*vMean)[0];
	
	
	//input
	TFile *inputFile = TFile::Open("FullSimulation_LHC17pCENTwSDD_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");


	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	//output
	TFile *outputFile = TFile::Open("bkg_from_sidebands_pT_"+pTRange+".root", "recreate");
	TTree *outputTree = new TTree("data_sidebands", "Sidebands from alice data for background training");
	
	//vars
	float pt_cand, d_len, d_len_xy, norm_dl_xy, cos_p, cos_p_xy, imp_par_xy, pt_prong0, pt_prong1, pt_prong2, sig_vert;
	float invM;
	bool isselectedstd;
	
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
	outputTree->Branch("inv_mass", &invM);
	outputTree->Branch("isselectedstd", &isselectedstd);
	
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
	inputTree->SetBranchAddress("inv_mass", &invM);
	inputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	
	//tree loop
	int nEntries = inputTree->GetEntries();

	for(int i = 0; i < nEntries; i++)
	{
		inputTree->GetEntry(i);
		
		if(invM<(mean-3*sigma) || invM>(mean+3*sigma))
		{
			if(pt_cand>=pTMin && pt_cand <=pTMax)
			{
				outputTree->Fill();
			}
		}
	}
	
	outputFile->Write("",TObject::kOverwrite);
	
	delete inputList;
	delete vSigma;
	delete vMean;
	delete inputFile;
	delete outputFile;
	delete sigmaFile;
}

void CreateSigFile(float pTMin, float pTMax)
{

	std::cout << Form("Creating Monte-Carlo signal training file for pT range from %.2f to %.2f", pTMin, pTMax) << std::endl;

	TString pTRange = Form("%.2f_%.2f",pTMin, pTMax);
	
	//input
	TFile *inputFile = TFile::Open("FullSimulation_LHC18a4a2CENT_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");


	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	//output
	TFile* outputFile = TFile::Open("sig_from_MC_pT_"+pTRange+".root", "recreate");
	TTree* outputTree = new TTree("MC_sig","signal from ALICE MC");
	
	//vars
	float pt_cand, d_len, d_len_xy, norm_dl_xy, cos_p, cos_p_xy, imp_par_xy, pt_prong0, pt_prong1, pt_prong2, sig_vert;
	float invM;
	bool isselectedstd;
	char cand_type;
	
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
	outputTree->Branch("inv_mass", &invM);
	outputTree->Branch("isselectedstd", &isselectedstd);
	
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
	inputTree->SetBranchAddress("inv_mass", &invM);
	inputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	inputTree->SetBranchAddress("cand_type", &cand_type);
	
	int nEntries = inputTree->GetEntries();
	
	for(int i = 0; i < nEntries; i++)
	{
		inputTree->GetEntry(i);
		
		if(cand_type!=0 && pt_cand >= pTMin && pt_cand <= pTMax)
		{
			outputTree->Fill();
		}
	
	}
	
	outputFile->Write("",TObject::kOverwrite);
	
	delete inputList;
	delete inputFile;
	delete outputFile;
}

void runSideBandAllpTBins()
{
	int pTLowArr[] = 	{2,3,4,5,6,6,7,8 ,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,8,10,12,12,16,24};
	
	for(int i=0; i<12; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		CreateSidebandFile(pTLow, pTHigh);
	}


}

void RunSigFileAllpTBins()
{
	int pTLowArr[] = 	{2,3,4,5,6,6,7,8 ,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,8,10,12,12,16,24};
	
	for(int i=0; i<12; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		CreateSigFile(pTLow, pTHigh);
	}
	

}

void CreateBkgFile(float pTMin, float pTMax)
{
	TString pTRange = Form("%.2f_%.2f",pTMin, pTMax);

	float sigma;
	float mean;
	
	//to get sigma and mean
	TString fSigmaName = Form("MC_true_signal_pT_%.2f_%.2f.root", pTMin,  pTMax);
	TString dirSigma = "MC_true_signal/";
	
	if(gSystem->AccessPathName(dirSigma+fSigmaName))
	{
		//file does not exist yet
		gROOT->ProcessLine(".L get_sigma.C");
		gROOT->ProcessLine(Form("GetSigma(%.2f,%.2f)",pTMin, pTMax));	
	}
	
	TFile *sigmaFile = TFile::Open(dirSigma+fSigmaName);
	TVectorF *vSigma = (TVectorF*)sigmaFile->Get("sigma");
	TVectorF *vMean = (TVectorF*)sigmaFile->Get("mean");
	sigma = TMath::Abs((*vSigma)[0]);
	mean = (*vMean)[0];
	
	sigmaFile->Close();
	
	//input
	TFile *inputFile = TFile::Open("FullSimulation_LHC18a4a2CENT_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");


	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	//output
	TFile* outputFile = TFile::Open("bkg_from_MC_pT_"+pTRange+".root", "recreate");
	TTree* outputTree = new TTree("MC_bkg","background from ALICE MC in 3 sigma around D+ invM mean");
	
	//vars
	float pt_cand, d_len, d_len_xy, norm_dl_xy, cos_p, cos_p_xy, imp_par_xy, pt_prong0, pt_prong1, pt_prong2, sig_vert;
	float invM;
	bool isselectedstd;
	char cand_type;
	
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
	outputTree->Branch("inv_mass", &invM);
	outputTree->Branch("isselectedstd", &isselectedstd);
	
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
	inputTree->SetBranchAddress("inv_mass", &invM);
	inputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	inputTree->SetBranchAddress("cand_type", &cand_type);
	
	int nEntries = inputTree->GetEntries();
	
	for(int i = 0; i < nEntries; i++)
	{
		inputTree->GetEntry(i);
		
		if(cand_type==0 && invM>=(mean-3*sigma) && invM<=(mean+3*sigma) && pt_cand>=pTMin && pt_cand <=pTMax)
		{
			outputTree->Fill();
		}
	
	}
	
	outputFile->Write("",TObject::kOverwrite);
	
	delete inputFile;
	delete outputFile;
}

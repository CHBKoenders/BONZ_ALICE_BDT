#include "TTree.h"
#include "TFile.h"

#include "TMath.h"

//sigma used is for the entire data set, nog a pT range -> next step
void createFile(float mean = 1.86991, float sigma = 8.72107e-03)
{
	sigma = TMath::Abs(sigma);
	
	//input
	TFile *inputFile = TFile::Open("FullSimulation_LHC17pCENTwSDD_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");


	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	//output
	TFile *outputFile = TFile::Open("bkg_from_sidebands.root", "recreate");
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

			outputTree->Fill();
		}
	}
	
	outputFile->Write("",TObject::kOverwrite);
	
	delete outputFile;
}

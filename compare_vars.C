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
#include "TH1D.h"
#include "TRatioPlot.h"

#include "TLegend.h"
#include "TStyle.h"

void Compare(float pTMin=0, float pTMax=90)
{	
	TFile *dataFile = TFile::Open("FullSimulation_LHC17pCENTwSDD_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");
	TList *dataList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *dataTree = (TTree*)dataList->First();
	
	TFile *MCFile = TFile::Open("FullSimulation_LHC18a4a2CENT_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");
	TList *MCList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *MCTree = (TTree*)MCList->First();
	
	//vars
	float d_len, d_len_xy, norm_dl_xy, cos_p, cos_p_xy, imp_par_xy, pt_prong0, pt_prong1, pt_prong2, sig_vert, pt_cand;
	
		
	std::string nameArr[] = {"sig_vert", "pT", "d_len", "d_len_xy", "norm_dl_xy", "cos_p", "cos_p_xy", "imp_par_xy", "pt_prong0", "pt_prong1", "pt_prong2"};
	int numVars = 11;
	
	std::map<std::string, float*>  vars;
	vars["sig_vert"] = &sig_vert;
	vars["pT"] = &pt_cand;	
	vars["d_len"] = &d_len;
	vars["d_len_xy"] = &d_len_xy;
	vars["norm_dl_xy"] = &norm_dl_xy;
	vars["cos_p"] = &cos_p;
	vars["cos_p_xy"] = &cos_p_xy;
	vars["imp_par_xy"] = &imp_par_xy;
	vars["pt_prong0"] = &pt_prong0;
	vars["pt_prong1"] = &pt_prong1;
	vars["pt_prong2"] = &pt_prong2;
	
	
	dataTree->SetBranchAddress("sig_vert", &sig_vert);
	dataTree->SetBranchAddress("pt_cand", &pt_cand);
	dataTree->SetBranchAddress("d_len", &d_len);
	dataTree->SetBranchAddress("d_len_xy", &d_len_xy);
	dataTree->SetBranchAddress("norm_dl_xy", &norm_dl_xy);
	dataTree->SetBranchAddress("cos_p", &cos_p);
	dataTree->SetBranchAddress("cos_p_xy", &cos_p_xy);
	dataTree->SetBranchAddress("imp_par_xy", &imp_par_xy);
	dataTree->SetBranchAddress("pt_prong0", &pt_prong0);
	dataTree->SetBranchAddress("pt_prong1", &pt_prong1);
	dataTree->SetBranchAddress("pt_prong2", &pt_prong2);

	
	MCTree->SetBranchAddress("sig_vert", &sig_vert);
	MCTree->SetBranchAddress("pt_cand", &pt_cand);
	MCTree->SetBranchAddress("d_len", &d_len);
	MCTree->SetBranchAddress("d_len_xy", &d_len_xy);
	MCTree->SetBranchAddress("norm_dl_xy", &norm_dl_xy);
	MCTree->SetBranchAddress("cos_p", &cos_p);
	MCTree->SetBranchAddress("cos_p_xy", &cos_p_xy);
	MCTree->SetBranchAddress("imp_par_xy", &imp_par_xy);
	MCTree->SetBranchAddress("pt_prong0", &pt_prong0);
	MCTree->SetBranchAddress("pt_prong1", &pt_prong1);
	MCTree->SetBranchAddress("pt_prong2", &pt_prong2);
	
	TFile *workFile = TFile::Open("workFile.root","RECREATE");
	
	//data loop
	int nDataEntries = dataTree->GetEntries();
	TH1D *sig_vertDataHist = new TH1D("sig_vertDataHist", "sig_vert", 200,0,0.06);
	sig_vertDataHist->SetLineColor(4);
	TH1D *pTDataHist = new TH1D("pTDataHist", "pT", 500,0,24);
	TH1D *d_lenDataHist = new TH1D("d_lenDataHist", "d_len", 200,0,5);
	TH1D *d_len_xyDataHist = new TH1D("d_len_xyDataHist", "d_len_xy", 200,0,3);
	TH1D *norm_dl_xyDataHist = new TH1D("norm_dl_xyDataHist", "norm_dl_xy", 200,0,600);
	TH1D *cos_pDataHist = new TH1D("cos_pDataHist", "cos_p", 200,0,1);
	TH1D *cos_p_xyDataHist = new TH1D("cos_p_xyDataHist", "cos_p_xy", 200,0,1);
	TH1D *imp_par_xyDataHist = new TH1D("imp_par_xyDataHist", "imp_par_xy", 200,0,2.5);
	TH1D *pt_prong0DataHist = new TH1D("pt_prong0DataHist", "pt_prong0",500,0,50);
	TH1D *pt_prong1DataHist = new TH1D("pt_prong1DataHist", "pt_prong1", 500,0,120);
	TH1D *pt_prong2DataHist = new TH1D("pt_prong2DataHist", "pt_prong2", 500,0,60);
	
	TString histName;
	TH1D* currHistData;
	float currVarVal;
	for(int i = 0; i<nDataEntries; i++)
	{
		dataTree->GetEntry(i);
		if(pt_cand >= pTMin && pt_cand <= pTMax)
		{
			for(int j = 0; j<numVars; j++)
			{
				histName = nameArr[j]+"DataHist";
				currHistData = (TH1D*)workFile->Get(histName);
				
				currVarVal = *(vars[nameArr[j]]);
				
				currHistData->Fill(currVarVal);
			}
		}
		
	}
	
	//MC loop
	int nMCEntries = MCTree->GetEntries();
	TH1D *sig_vertMCHist = new TH1D("sig_vertMCHist", "sig_vert", 200,0,0.06);
	TH1D *pTMCHist = new TH1D("pTMCHist", "pT", 500,0,24);
	TH1D *d_lenMCHist = new TH1D("d_lenMCHist", "d_len", 200,0,5);
	TH1D *d_len_xyMCHist = new TH1D("d_len_xyMCHist", "d_len_xy", 200,0,3);
	TH1D *norm_dl_xyMCHist = new TH1D("norm_dl_xyMCHist", "norm_dl_xy", 200,0,600);
	TH1D *cos_pMCHist = new TH1D("cos_pMCHist", "cos_p", 200,0,1);
	TH1D *cos_p_xyMCHist = new TH1D("cos_p_xyMCHist", "cos_p_xy", 200,0,1);
	TH1D *imp_par_xyMCHist = new TH1D("imp_par_xyMCHist", "imp_par", 200,0,2.5);
	TH1D *pt_prong0MCHist = new TH1D("pt_prong0MCHist", "pt_prong0", 500,0,50);
	TH1D *pt_prong1MCHist = new TH1D("pt_prong1MCHist", "pt_prong1", 500,0,120);
	TH1D *pt_prong2MCHist = new TH1D("pt_prong2MCHist", "pt_prong2", 500,0,60);
	
	TH1D* currHistMC;
	for(int i = 0; i<nMCEntries; i++)
	{
		MCTree->GetEntry(i);
		if(pt_cand >= pTMin && pt_cand <= pTMax)
		{
			for(int j = 0; j<numVars; j++)
			{
				histName = nameArr[j]+"MCHist";
				currHistMC = (TH1D*)workFile->Get(histName);
				
				currVarVal = *(vars[nameArr[j]]);
				
				currHistMC->Fill(currVarVal);
			}
		}
	}
	
	TCanvas *c1 = new TCanvas;
	
	TLegend *legend = new TLegend(3,15);
	
	legend->AddEntry(sig_vertDataHist,"Data","L");
	legend->AddEntry(sig_vertMCHist,"MC","L");
	
	c1->Divide(4,3);
	
	c1->cd(1);
	
	legend->Draw();
	
	TRatioPlot *currRatioPlot;
	double dataInt, MCInt;
	int firstBin, lastBin;
	double xLowData, xLowMC, xHighData,xHighMC;
	double xMin, xMax;
	for(int i = 0; i < numVars; i++)
	{
		histName = nameArr[i]+"DataHist";
		currHistData = (TH1D*)workFile->Get(histName);
		histName = nameArr[i]+"MCHist";
		currHistMC = (TH1D*)workFile->Get(histName);
		
		firstBin = currHistData->FindFirstBinAbove(0);
		xLowData = currHistData->GetBinLowEdge(firstBin);
		lastBin = currHistData->FindLastBinAbove(0);
		xHighData = currHistData->GetBinLowEdge(lastBin) + currHistData->GetBinWidth(lastBin);
		
		firstBin = currHistMC->FindFirstBinAbove(0);
		xLowMC = currHistMC->GetBinLowEdge(firstBin);
		lastBin = currHistMC->FindLastBinAbove(0);
		xHighMC = currHistMC->GetBinLowEdge(lastBin) + currHistMC->GetBinWidth(lastBin);
		
		xMin = TMath::Max(xLowData, xLowMC);
		xMax = TMath::Min(xHighData, xHighMC);
		
		currHistData->SetAxisRange(xMin, xMax);
		currHistMC->SetAxisRange(xMin, xMax);
		
		currHistData->SetLineColor(4);
		currHistMC->SetLineColor(2);
		
		currHistData->SetStats(0);
		currHistMC->SetStats(0);
			
		dataInt = currHistData->Integral();
		MCInt = currHistMC->Integral();
		
		currHistData->Scale(1/dataInt);
		currHistMC->Scale(1/MCInt);
		
		currRatioPlot = new TRatioPlot(currHistData, currHistMC);
		
		c1->cd(i+2);
	
		currRatioPlot->DrawClone();
	}
	
	
}

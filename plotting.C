#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"

#include "TPaveText.h"
#include "TStyle.h"

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

void PlotInvMassHist(float pTMin, float pTMax, bool BDT, float horScale=1, float vertScale=0.5)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	
	TFile *fitFile;
	TString fileName;
	TString typeS;
	TH1D *dataHist;
	TString histTitle;
	
	if(BDT)
	{
		typeS = "Boosted Decision Tree";
		//open BDT fit file
		fileName = "Fitting_results/BDT_fit_pT_"+pTRange+".root";
		
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		fitFile = TFile::Open(fileName);
		
		dataHist = (TH1D*)fitFile->Get("sigHist");
		histTitle = Form("Fitting result for BDT analysis in #it{p}_{T} range from %.2f to %.2f",pTMin, pTMax);
	}
	else
	{
		typeS = "ALICE Standard";
		//open Std fit file
		fileName = "Fitting_results/Std_fit_pT_"+pTRange+".root";
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		
		fitFile = TFile::Open(fileName);
		
		dataHist = (TH1D*)fitFile->Get("selectedHist");
		histTitle = Form("Fitting result for regular ALICE analysis in #it{p}_{T} range from %.2f to %.2f",pTMin, pTMax);
	}
	
	float binWidth = dataHist->GetBinWidth(dataHist->GetMaximumBin());
	//convert from GeV to MeV
	binWidth *= 1000;

	dataHist->SetTitle(histTitle);
	
	dataHist->GetXaxis()->SetTitle("Invariant Mass (GeV/#it{c}^{2})");
	dataHist->GetXaxis()->CenterTitle();
	dataHist->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}",binWidth));
	dataHist->GetYaxis()->CenterTitle();
	
	//so stat box does not get drawn
	dataHist->SetStats(0);
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	TCanvas *c = new TCanvas;
	c->SetCanvasSize(horScale*horSize, vertScale*vertSize);
	//No title
	gStyle->SetOptTitle(0);

	TPaveText *textBox = new TPaveText(.65,.65,.85,.75,"blNDC");
	
	textBox->AddText("Analysis: "+typeS);
	textBox->AddText(Form("%.2f #leq #it{p}_{T} #leq %.2f GeV/#it{c}", pTMin, pTMax));
	textBox->AddText(Form("Entries: %.0f", dataHist->GetEntries()));
	
	dataHist->Draw("E");
	dataHist->GetXaxis()->SetRangeUser(1.7, 2.06);
	
	textBox->SetBorderSize(1);
	textBox->SetFillStyle(0);
	textBox->SetTextFont(42);
	textBox->SetTextSize(.02);
	
	textBox->Draw();


}

void PlotFit(float pTMin, float pTMax, bool BDT, float horScale = 0.5, float vertScale = 0.5)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	
	TFile *fitFile;
	TString fileName;
	TString typeS;
	TH1D *dataHist;
	TString histTitle;
	
	if(BDT)
	{
		typeS = "Boosted Decision Tree";
		//open BDT fit file
		fileName = "Fitting_results/BDT_fit_pT_"+pTRange+".root";
		
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		fitFile = TFile::Open(fileName);
		
		dataHist = (TH1D*)fitFile->Get("sigHist");
		histTitle = Form("Fitting result for BDT analysis in #it{p}_{T} range from %.2f to %.2f",pTMin, pTMax);
	}
	else
	{
		typeS = "ALICE Standard";
		//open Std fit file
		fileName = "Fitting_results/Std_fit_pT_"+pTRange+".root";
		if(gSystem->AccessPathName(fileName))
		{
			std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
			
			return;
		}
		
		fitFile = TFile::Open(fileName);
		
		dataHist = (TH1D*)fitFile->Get("selectedHist");
		histTitle = Form("Fitting result for regular ALICE analysis in #it{p}_{T} range from %.2f to %.2f",pTMin, pTMax);
	}
	
	float binWidth = dataHist->GetBinWidth(dataHist->GetMaximumBin());
	//convert from GeV to MeV
	binWidth *= 1000;
	
	TVectorF *sigYieldV = (TVectorF*)fitFile->Get("sigYield");
	float sigYield = (*sigYieldV)[0];
	float sigYieldErr = (*sigYieldV)[1];
	TVectorF *bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
	float bkgYield = (*bkgYieldV)[0];
	float bkgYieldErr = (*bkgYieldV)[1];
	TVectorF *significanceV = (TVectorF*)fitFile->Get("significance");
	float significance = (*significanceV)[0];
	float significanceErr = (*significanceV)[1];
	
	//retrieve fit
	TF1 *fit = dataHist->GetFunction("signal");

	//now we retrieve the background function	
	TF1 *background = new TF1("background", "pol2(0)",1.673,2.067);
	background->SetParameters(fit->GetParameter(3),fit->GetParameter(4),fit->GetParameter(5));

	background->SetTitle("background");
	background->SetLineColor(1);
	background->SetLineStyle(2);

	fit->SetTitle("fit");
	
	dataHist->SetTitle(histTitle);
	
	dataHist->GetXaxis()->SetTitle("Invariant Mass (GeV/#it{c}^{2})");
	dataHist->GetXaxis()->CenterTitle();
	dataHist->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}",binWidth));
	dataHist->GetYaxis()->CenterTitle();
	
	//so stat box does not get drawn
	dataHist->SetStats(0);
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	TCanvas *c = new TCanvas;
	c->SetCanvasSize(horScale*horSize, vertScale*vertSize);
	//No title
	gStyle->SetOptTitle(0);

	TPaveText *textBox = new TPaveText(.65,.55,.85,.85,"blNDC");

	
	textBox->AddText("Analysis: "+typeS);
	textBox->AddText(Form("%.2f #leq #it{p}_{T} #leq %.2f GeV/#it{c}", pTMin, pTMax));
	textBox->AddText(Form("#mu = %.4f #pm %5f", fit->GetParameter("mean"), fit->GetParError(fit->GetParNumber("mean"))));
	textBox->AddText(Form("#sigma = %.4f #pm %5f", fit->GetParameter("sigma"), fit->GetParError(fit->GetParNumber("sigma"))));
	textBox->AddText(Form("S (3#sigma) = %.2f #pm %.2f",sigYield, sigYieldErr));
	textBox->AddText(Form("#frac{S}{#sqrt{S+B}} (3#sigma) = %.2f #pm %.2f", significance, significanceErr));
	
	dataHist->Draw("E");
	dataHist->GetXaxis()->SetRangeUser(1.7, 2.06);
	background->Draw("LSame");
	fit->Draw("LSame");
	
	TLegend *legend = new TLegend(0.15,0.75,0.4,0.85);
	legend->AddEntry(dataHist, "Data","L");
	legend->AddEntry(background, "Background fit", "L");
	legend->AddEntry(fit, "Total fit", "L");	
	
	legend->Draw();
	
	textBox->SetBorderSize(1);
	textBox->SetFillStyle(0);
	textBox->SetTextFont(42);
	textBox->SetTextSize(.02);
	
	textBox->Draw();

}

void PlotFitAllpTBins()
{
	int pTMinArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTMaxArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int nBins = 10;
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	TCanvas *cBDT = new TCanvas;
	cBDT->SetCanvasSize(horSize, vertSize);
	cBDT->Divide(2,5);
	gStyle->SetOptTitle(0);
	
	TCanvas *cStd = new TCanvas;
	cStd->SetCanvasSize(horSize, vertSize);
	cStd->Divide(2,5);
	gStyle->SetOptTitle(0);
	
	float pTMin;
	float pTMax;
	TFile *fitFile;
	TString pTRange;
	TString fileName;
	TString typeS;
	TString histTitle;
	TVectorF *sigYieldV;
	TVectorF *bkgYieldV;
	TVectorF *significanceV;
	float sigYield, sigYieldErr, bkgYield, bkgYieldErr, significance, significanceErr;
	for(int i=0; i<nBins; i++)
	{
		pTMin = pTMinArr[i];	
		pTMax = pTMaxArr[i];
		
		pTRange = Form("%.2f_%.2f", pTMin, pTMax);
		
		TString typeS;
		for(int j=0; j<2; j++)
		{
			TH1D *dataHist;
			
			if(!j)
			{
				typeS = "Boosted Decision Tree";
				cBDT->cd(i+1);
				
				//open BDT fit file
				fileName = "Fitting_results/BDT_fit_pT_"+pTRange+".root";
				
				if(gSystem->AccessPathName(fileName))
				{
					std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
					
					return;
				}
				fitFile = TFile::Open(fileName);
				
				dataHist = (TH1D*)fitFile->Get("sigHist");
				histTitle = Form("Fitting result for BDT analysis in #it{p}_{T} range from %.2f to %.2f GeV/#it{c}",pTMin, pTMax);
			}
			else
			{
				typeS = "ALICE Standard";
				cStd->cd(i+1);
				
				//open Std fit file
				fileName = "Fitting_results/Std_fit_pT_"+pTRange+".root";
				if(gSystem->AccessPathName(fileName))
				{
					std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
					
					return;
				}
				
				fitFile = TFile::Open(fileName);
				
				dataHist = (TH1D*)fitFile->Get("selectedHist");
				histTitle = Form("Fitting result for regular ALICE analysis in #it{p}_{T} range from %.2f to %.2f GeV/#it{c}",pTMin, pTMax);
			}
			
			float binWidth = dataHist->GetBinWidth(dataHist->GetMaximumBin());
			//convert from GeV to MeV
			binWidth *= 1000;
			
			sigYieldV = (TVectorF*)fitFile->Get("sigYield");
			sigYield = (*sigYieldV)[0];
			sigYieldErr = (*sigYieldV)[1];
			bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
			bkgYield = (*bkgYieldV)[0];
			bkgYieldErr = (*bkgYieldV)[1];
			significanceV = (TVectorF*)fitFile->Get("significance");
			significance = (*significanceV)[0];
			significanceErr = (*significanceV)[1];
			
			//retrieve fit
			TF1 *fit = dataHist->GetFunction("signal");

			//now we retrieve the background function	
			TF1 *background = new TF1("background", "pol2(0)",1.673,2.067);
			background->SetParameters(fit->GetParameter(3),fit->GetParameter(4),fit->GetParameter(5));

			background->SetTitle("background");
			background->SetLineColor(1);
			background->SetLineStyle(2);

			fit->SetTitle("fit");
			
			dataHist->SetTitle(histTitle);
			
			dataHist->GetXaxis()->SetTitle("Invariant Mass (GeV/#it{c}^{2})");
			dataHist->GetXaxis()->CenterTitle();
			dataHist->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}",binWidth));
			dataHist->GetYaxis()->CenterTitle();
			
			//so stat box does not get drawn
			dataHist->SetStats(0);
			
			TPaveText *textBox = new TPaveText(.65,.50,.85,.85,"blNDC");
			
			textBox->AddText("Analysis: "+typeS);
			textBox->AddText(Form("%.2f #leq #it{p}_{T} #leq %.2f GeV/#it{c}", pTMin, pTMax));
			textBox->AddText(Form("#mu = %.4f #pm %5f", fit->GetParameter("mean"), fit->GetParError(fit->GetParNumber("mean"))));
			textBox->AddText(Form("#sigma = %.4f #pm %5f", fit->GetParameter("sigma"), fit->GetParError(fit->GetParNumber("sigma"))));
			textBox->AddText(Form("S (3#sigma) = %.2f #pm %.2f",sigYield, sigYieldErr));
			textBox->AddText(Form("#frac{S}{#sqrt{S+B}} (3#sigma) = %.2f #pm %.2f", significance, significanceErr));
			
			dataHist->Draw("E");
			dataHist->GetXaxis()->SetRangeUser(1.7, 2.06);
			background->Draw("LSame");
			fit->Draw("LSame");
			
			TLegend *legend = new TLegend(0.15,0.75,0.4,0.85);
			legend->AddEntry(dataHist, "Data","L");
			legend->AddEntry(background, "Background fit", "L");
			legend->AddEntry(fit, "Total fit", "L");	
			
			legend->Draw();
			
			textBox->SetBorderSize(1);
			textBox->SetFillStyle(0);
			textBox->SetTextFont(42);
			//textBox->SetTextSizePixels(80);
			textBox->SetTextSize(0.025);
			
			textBox->Draw();
		}			
	}
	
	

}

void PlotRawYields(float horScale = 1, float vertScale = .25)
{
	int pTMinArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTMaxArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int nBins = 10;
	
	int horSize = 1920;
	int vertSize = TMath::Sqrt(2)*horSize/4;
	
	TCanvas *cSig = new TCanvas;
	
	TGraphErrors *rawSigYieldGraphBDT = new TGraphErrors;
	TGraphErrors *rawSigYieldGraphStd = new TGraphErrors;
	
	cSig->SetCanvasSize(horSize, vertSize);
	
	TCanvas *cBkg = new TCanvas;
	
	TGraphErrors *rawBkgYieldGraphBDT = new TGraphErrors;
	TGraphErrors *rawBkgYieldGraphStd = new TGraphErrors;
	
	cBkg->SetCanvasSize(horSize, vertSize);
	
	float pTMin;
	float pTMax;
	TFile *fitFile;
	TString pTRange;
	TString fileName;
	TString typeS;
	TString histTitle;
	TVectorF *sigYieldV;
	TVectorF *bkgYieldV;
	
	float sigYield, sigYieldErr, bkgYield, bkgYieldErr;
	
	//pT bin loop
	for(int i=0; i<nBins; i++)
	{
		pTMin = pTMinArr[i];	
		pTMax = pTMaxArr[i];
		
		pTRange = Form("%.2f_%.2f", pTMin, pTMax);
		
		//BDT & Std loop
		for(int j=0; j<2; j++)
		{
			TH1D *dataHist;
			
			if(!j)
			{				
				//open BDT fit file
				fileName = "Fitting_results/BDT_fit_pT_"+pTRange+".root";
				
				if(gSystem->AccessPathName(fileName))
				{
					std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
					
					return;
				}
				fitFile = TFile::Open(fileName);
				
				dataHist = (TH1D*)fitFile->Get("sigHist");
				histTitle = Form("Fitting result for BDT analysis in #it{p}_{T} range from %.2f to %.2f GeV/#it{c}",pTMin, pTMax);
			}
			else
			{				
				//open Std fit file
				fileName = "Fitting_results/Std_fit_pT_"+pTRange+".root";
				if(gSystem->AccessPathName(fileName))
				{
					std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
					
					return;
				}
				
				fitFile = TFile::Open(fileName);
				
				dataHist = (TH1D*)fitFile->Get("selectedHist");
				histTitle = Form("Fitting result for regular ALICE analysis in #it{p}_{T} range from %.2f to %.2f Gev/#it{c}",pTMin, pTMax);
			}
			
			
			sigYieldV = (TVectorF*)fitFile->Get("sigYield");
			sigYield = (*sigYieldV)[0];
			sigYieldErr = (*sigYieldV)[1];
			bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
			bkgYield = (*bkgYieldV)[0];
			bkgYieldErr = (*bkgYieldV)[1];
		
			
			if(!j)
			{
				rawSigYieldGraphBDT->SetPoint(i, (pTMin+pTMax)/2., sigYield);
				rawSigYieldGraphBDT->SetPointError(i,(pTMax-pTMin)/2.,sigYieldErr);
				
				rawBkgYieldGraphBDT->SetPoint(i, (pTMin+pTMax)/2., bkgYield);
				rawBkgYieldGraphBDT->SetPointError(i,(pTMax-pTMin)/2.,bkgYieldErr);
				

			}
			else
			{
				rawSigYieldGraphStd->SetPoint(i, (pTMin+pTMax)/2., sigYield);
				rawSigYieldGraphStd->SetPointError(i,(pTMax-pTMin)/2.,sigYieldErr);
				
				rawBkgYieldGraphStd->SetPoint(i, (pTMin+pTMax)/2., bkgYield);
				rawBkgYieldGraphStd->SetPointError(i,(pTMax-pTMin)/2.,bkgYieldErr);
			}
		}
	}
	
	//plotting
	
	cSig->cd();
	//no title
	gStyle->SetOptTitle(0);
	
	TMultiGraph *rawSigYieldGraph = new TMultiGraph();
	rawSigYieldGraph->Add(rawSigYieldGraphBDT);
	rawSigYieldGraph->Add(rawSigYieldGraphStd);
	
	rawSigYieldGraph->SetTitle("Raw signal yield in #it{p}_{T} intervals");
	
	rawSigYieldGraphBDT->SetTitle("Boosted Decision Tree");
	rawSigYieldGraphBDT->SetLineColor(2);
	rawSigYieldGraphBDT->SetMarkerStyle(8);
	
	rawSigYieldGraphStd->SetTitle("ALICE standard");
	rawSigYieldGraphStd->SetLineColor(4);
	rawSigYieldGraphStd->SetMarkerStyle(22);
	
	rawSigYieldGraph->Draw("AP");
	
	rawSigYieldGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	rawSigYieldGraph->GetXaxis()->CenterTitle();
	rawSigYieldGraph->GetYaxis()->SetTitle("Raw Signal Yield");
	rawSigYieldGraph->GetYaxis()->CenterTitle();	
	
	gPad->BuildLegend();
	
	//////////////////////////////////////////////////////////////////////////
	
	cBkg->cd();
	//no title
	gStyle->SetOptTitle(0);
	
	TMultiGraph *rawBkgYieldGraph = new TMultiGraph();
	rawBkgYieldGraph->Add(rawBkgYieldGraphBDT);
	rawBkgYieldGraph->Add(rawBkgYieldGraphStd);
	
	rawSigYieldGraph->SetTitle("Raw background yield in #it{p}_{T} intervals");
	
	rawBkgYieldGraphBDT->SetTitle("Boosted Decision Tree");
	rawBkgYieldGraphBDT->SetLineColor(2);
	rawBkgYieldGraphBDT->SetMarkerStyle(8);
	
	rawBkgYieldGraphStd->SetTitle("ALICE standard");
	rawBkgYieldGraphStd->SetLineColor(4);
	rawBkgYieldGraphStd->SetMarkerStyle(22);
	
	rawBkgYieldGraph->Draw("AP");
	
	rawBkgYieldGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	rawBkgYieldGraph->GetXaxis()->CenterTitle();
	rawBkgYieldGraph->GetYaxis()->SetTitle("Raw Background Yield");
	rawBkgYieldGraph->GetYaxis()->CenterTitle();	
	
	gPad->BuildLegend();	
	
}

void PlotStatGraphs()
{
	int pTMinArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTMaxArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int nBins = 10;
	
	int horSize = 1920;
	int vertSize = TMath::Sqrt(2)*horSize/4;
	
	
	TCanvas *cSig = new TCanvas;
	
	TGraphErrors *significanceGraphBDT = new TGraphErrors;
	TGraphErrors *significanceGraphStd = new TGraphErrors;
	
	cSig->SetCanvasSize(horSize, vertSize);
	
	
	TCanvas *cRatio = new TCanvas;
	
	TGraphErrors *ratioGraphBDT = new TGraphErrors;
	TGraphErrors *ratioGraphStd = new TGraphErrors;
	
	cRatio->SetCanvasSize(horSize, vertSize);
	
	
	TCanvas *cMean = new TCanvas;
	
	TGraphErrors *meanGraphBDT = new TGraphErrors;
	TGraphErrors *meanGraphStd = new TGraphErrors;
	TGraphErrors *meanGraphMC  = new TGraphErrors;
	
	cMean->SetCanvasSize(horSize, vertSize);
	
	
	TCanvas *cSigma = new TCanvas;
	
	TGraphErrors *sigmaGraphBDT = new TGraphErrors;
	TGraphErrors *sigmaGraphStd = new TGraphErrors;
	TGraphErrors *sigmaGraphMC  = new TGraphErrors;
	
	cSigma->SetCanvasSize(horSize, vertSize);
	
	
	
	
	float pTMin;
	float pTMax;
	TFile *fitFile;
	TFile *fMC;
	TString pTRange;
	TString fileName;
	TString typeS;
	TString histTitle;
	TVectorF *sigYieldV;
	TVectorF *bkgYieldV;
	TVectorF *significanceV;
	TVectorF *vMean;
	TVectorF *vSigma;
	float sigYield, sigYieldErr, bkgYield, bkgYieldErr, significance, significanceErr, ratio, ratioErr, mean, meanErr, sigma, sigmaErr;
	
	//pT bin loop
	for(int i=0; i<nBins; i++)
	{
		pTMin = pTMinArr[i];	
		pTMax = pTMaxArr[i];
		
		pTRange = Form("%.2f_%.2f", pTMin, pTMax);
		
		//BDT & Std loop
		for(int j=0; j<2; j++)
		{
			TH1D *dataHist;
			
			if(!j)
			{				
				//open BDT fit file
				fileName = "Fitting_results/BDT_fit_pT_"+pTRange+".root";
				
				if(gSystem->AccessPathName(fileName))
				{
					std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
					
					return;
				}
				fitFile = TFile::Open(fileName);
				
				dataHist = (TH1D*)fitFile->Get("sigHist");
				histTitle = Form("Fitting result for BDT analysis in #it{p}_{T} range from %.2f to %.2f GeV/#it{c}",pTMin, pTMax);
			}
			else
			{				
				//open Std fit file
				fileName = "Fitting_results/Std_fit_pT_"+pTRange+".root";
				if(gSystem->AccessPathName(fileName))
				{
					std::cout << "The file " << fileName << " does not exist! Run fit.C to fix this. Aborting operations" << std::endl;
					
					return;
				}
				
				fitFile = TFile::Open(fileName);
				
				dataHist = (TH1D*)fitFile->Get("selectedHist");
				histTitle = Form("Fitting result for regular ALICE analysis in #it{p}_{T} range from %.2f to %.2f Gev/#it{c}",pTMin, pTMax);
			}
			
			
			sigYieldV = (TVectorF*)fitFile->Get("sigYield");
			sigYield = (*sigYieldV)[0];
			sigYieldErr = (*sigYieldV)[1];
			bkgYieldV = (TVectorF*)fitFile->Get("bkgYield");
			bkgYield = (*bkgYieldV)[0];
			bkgYieldErr = (*bkgYieldV)[1];
			significanceV = (TVectorF*)fitFile->Get("significance");
			significance = (*significanceV)[0];
			significanceErr = (*significanceV)[1];
			
			ratio = sigYield/bkgYield;
			ratioErr = TMath::Abs(ratio)*TMath::Sqrt((sigYieldErr/sigYield)*(sigYieldErr/sigYield)+(bkgYieldErr/bkgYield)*(bkgYieldErr/bkgYield));
			
			//retrieve fit
			TF1 *fit = dataHist->GetFunction("signal");
			
			mean = fit->GetParameter("mean");
			meanErr = fit->GetParError(fit->GetParNumber("mean"));
			
			sigma = fit->GetParameter("sigma");
			sigmaErr = fit->GetParError(fit->GetParNumber("sigma"));
			
			if(!j)
			{
				significanceGraphBDT->SetPoint(i,(pTMin+pTMax)/2.,significance);
				significanceGraphBDT->SetPointError(i,(pTMax-pTMin)/2.,significanceErr);
				
				ratioGraphBDT->SetPoint(i,(pTMin+pTMax)/2.,ratio);
				ratioGraphBDT->SetPointError(i,(pTMax-pTMin)/2.,ratioErr);
				
				meanGraphBDT->SetPoint(i,(pTMin+pTMax)/2.,mean);
				meanGraphBDT->SetPointError(i,(pTMax-pTMin)/2.,meanErr);
				
				sigmaGraphBDT->SetPoint(i,(pTMin+pTMax)/2.,sigma);
				sigmaGraphBDT->SetPointError(i,(pTMax-pTMin)/2.,sigmaErr);
			}
			else
			{
				significanceGraphStd->SetPoint(i, (pTMin+pTMax)/2., significance);
				significanceGraphStd->SetPointError(i, (pTMax-pTMin)/2., significanceErr);
				
				ratioGraphStd->SetPoint(i, (pTMin+pTMax)/2., ratio);
				ratioGraphStd->SetPointError(i, (pTMax-pTMin)/2., ratioErr);
				
				meanGraphStd->SetPoint(i, (pTMin+pTMax)/2., mean);
				meanGraphStd->SetPointError(i, (pTMax-pTMin)/2., meanErr);
				
				sigmaGraphStd->SetPoint(i, (pTMin+pTMax)/2., sigma);
				sigmaGraphStd->SetPointError(i, (pTMax-pTMin)/2., sigmaErr);
			}
		}
		
		//MC values
		
		fileName = Form("MC_true_signal/MC_true_signal_pT_%.2f_%.2f.root", pTMin,  pTMax);
		
		if(gSystem->AccessPathName(fileName))
		{
			//file does not exist
			std::cout << "MC file not found" << std::endl;
			
			return;	
		}
		
		fMC = TFile::Open(fileName);
		vMean = (TVectorF*)fMC->Get("mean");
		vSigma = (TVectorF*)fMC->Get("sigma");
		
		mean = (*vMean)[0];
		meanErr = (*vMean)[1];
	
		sigma = (*vSigma)[0];
		sigmaErr = (*vSigma)[1];		
		
		meanGraphMC->SetPoint(i, (pTMin+pTMax)/2., mean);
		meanGraphMC->SetPointError(i, (pTMax-pTMin)/2., meanErr);
		
		sigmaGraphMC->SetPoint(i, (pTMin+pTMax)/2., sigma);
		sigmaGraphMC->SetPointError(i, (pTMax-pTMin)/2., sigmaErr);		
	}
	
	//plotting
	
	cSig->cd();
	//no title
	gStyle->SetOptTitle(0);
	
	TMultiGraph *significanceGraph = new TMultiGraph();
	significanceGraph->Add(significanceGraphBDT);
	significanceGraph->Add(significanceGraphStd);
	
	significanceGraph->SetTitle("Significance in #it{p}_{T} intervals");
	
	significanceGraphBDT->SetTitle("Boosted Decision Tree");
	significanceGraphBDT->SetLineColor(2);
	//significanceGraphBDT->SetFillColor(2);
	//significanceGraphBDT->SetFillStyle(3001);
	significanceGraphBDT->SetMarkerStyle(8);
	significanceGraphStd->SetTitle("ALICE standard");
	significanceGraphStd->SetLineColor(4);
	//significanceGraphStd->SetFillColor(4);
	//significanceGraphStd->SetFillStyle(3001);
	significanceGraphStd->SetMarkerStyle(22);
	
	significanceGraph->Draw("AP");
	
	significanceGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	significanceGraph->GetXaxis()->CenterTitle();
	significanceGraph->GetYaxis()->SetTitle("Significance (S/#sqrt{S+B})");
	significanceGraph->GetYaxis()->CenterTitle();	
	
	gPad->BuildLegend();
	
	/////////////////////////////////////////////////////////////////////
	
	cRatio->cd();
	//no title
	gStyle->SetOptTitle(0);
	
	TMultiGraph *ratioGraph = new TMultiGraph();
	ratioGraph->Add(ratioGraphBDT);
	ratioGraph->Add(ratioGraphStd);
	
	ratioGraph->SetTitle("Signal to Background Ratio in #it{p}_{T} intervals");
	
	ratioGraphBDT->SetTitle("Boosted Decision Tree");
	ratioGraphBDT->SetLineColor(2);
	//ratioGraphBDT->SetFillColor(2);
	//ratioGraphBDT->SetFillStyle(3001);
	ratioGraphBDT->SetMarkerStyle(8);
	ratioGraphStd->SetTitle("ALICE Standard");
	ratioGraphStd->SetLineColor(4);
	//ratioGraphStd->SetFillColor(4);
	//ratioGraphStd->SetFillStyle(3001);
	ratioGraphStd->SetMarkerStyle(22);
	
	ratioGraph->Draw("AP");
	
	ratioGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	ratioGraph->GetXaxis()->CenterTitle();
	ratioGraph->GetYaxis()->SetTitle("Ratio (S/B)");
	ratioGraph->GetYaxis()->CenterTitle();

	gPad->BuildLegend();
	
	//////////////////////////////////////////////////////////////////////
	
	cMean->cd();
	//no title
	gStyle->SetOptTitle(0);
	
	TMultiGraph *meanGraph = new TMultiGraph();
	meanGraph->Add(meanGraphMC,"PL 3");
	meanGraph->Add(meanGraphBDT);
	meanGraph->Add(meanGraphStd);
	
	meanGraph->SetTitle("Mean of Gaussian fit in #it{p}_{T} intervals");
	
	meanGraphBDT->SetTitle("Boosted Decision Tree");
	meanGraphBDT->SetLineColor(2);
	//ratioGraphBDT->SetFillColor(2);
	//ratioGraphBDT->SetFillStyle(3001);
	meanGraphBDT->SetMarkerStyle(8);
	meanGraphStd->SetTitle("ALICE Standard");
	meanGraphStd->SetLineColor(4);
	//ratioGraphStd->SetFillColor(4);
	//ratioGraphStd->SetFillStyle(3001);
	meanGraphStd->SetMarkerStyle(22);
	meanGraphMC->SetTitle("Monte-Carlo");
	meanGraphMC->SetLineColor(3);
	meanGraphMC->SetFillColor(3);
	meanGraphMC->SetFillStyle(3001);
	
	
	meanGraph->Draw("AP");
	
	meanGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	meanGraph->GetXaxis()->CenterTitle();
	meanGraph->GetYaxis()->SetTitle("#mu (GeV/#it{c}^{2})");
	meanGraph->GetYaxis()->CenterTitle();

	gPad->BuildLegend();
	
	////////////////////////////////////////////////////////////////////////
	
	cSigma->cd();
	//not title
	gStyle->SetOptTitle(0);
	
	TMultiGraph *sigmaGraph = new TMultiGraph();
	sigmaGraph->Add(sigmaGraphMC, "PL 3");
	sigmaGraph->Add(sigmaGraphBDT);
	sigmaGraph->Add(sigmaGraphStd);
	
	
	sigmaGraph->SetTitle("Sigma of Gaussian fit in #it{p}_{T} intervals");
	
	sigmaGraphBDT->SetTitle("Boosted Decision Tree");
	sigmaGraphBDT->SetLineColor(2);
	//ratioGraphBDT->SetFillColor(2);
	//ratioGraphBDT->SetFillStyle(3001);
	sigmaGraphBDT->SetMarkerStyle(8);
	sigmaGraphStd->SetTitle("ALICE Standard");
	sigmaGraphStd->SetLineColor(4);
	//ratioGraphStd->SetFillColor(4);
	//ratioGraphStd->SetFillStyle(3001);
	sigmaGraphStd->SetMarkerStyle(22);
	sigmaGraphMC->SetTitle("Monte-Carlo");
	sigmaGraphMC->SetLineColor(3);
	sigmaGraphMC->SetFillColor(3);
	sigmaGraphMC->SetFillStyle(3001);
	
	sigmaGraph->Draw("AP");
	
	sigmaGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	sigmaGraph->GetXaxis()->CenterTitle();
	sigmaGraph->GetYaxis()->SetTitle("#sigma (GeV/#it{c}^{2})");
	sigmaGraph->GetYaxis()->CenterTitle();

	gPad->BuildLegend();
}

void PlotCorrYieldAtCut(float pTMin, float pTMax, float horScale=1, float  vertScale=1)
{	
	TString fName = Form("eff_corrected_yield_at_BDT_cut/corr_yield_at_cut_pt_%.2f_%.2f.root", pTMin, pTMax);
	
	TFile *canvasFile = TFile::Open(fName);
	
	TCanvas *cPlot = (TCanvas*)canvasFile->Get("c1");
	
	cPlot->GetListOfPrimitives()->FindObject("title")->Clear();
	
	TLegend *legend = (TLegend*)cPlot->GetListOfPrimitives()->FindObject("TPave");
	
	((TLegendEntry*)legend->GetListOfPrimitives()->At(0))->SetLabel("Boosted Decision Tree");
	((TLegendEntry*)legend->GetListOfPrimitives()->At(1))->SetLabel("ALICE Standard");
	
	((TMultiGraph*)cPlot->GetListOfPrimitives()->FindObject(""))->GetYaxis()->SetTitle("Efficiency-Corrected Yield");
	((TMultiGraph*)cPlot->GetListOfPrimitives()->FindObject(""))->GetXaxis()->SetTitle("Boosted Decision Tree Cut Value");
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	TCanvas *c = new TCanvas;
	
	cPlot->DrawClonePad();
	
	TPaveText *textBox = new TPaveText(.45,.90,.55,.96,"blNDC");
	
	textBox->AddText(Form("%.2f #leq #it{p}_{T} #leq  %.2f GeV/#it{c}",pTMin, pTMax));
	
	textBox->SetBorderSize(0);
	textBox->SetFillStyle(0);
	textBox->SetTextFont(42);
	textBox->SetTextSize(.04);
	
	textBox->Draw();
	
	c->SetCanvasSize(horScale*horSize, vertScale*vertSize);	
}

void PlotCorrYieldAtCutAllpTBins()
{
	int pTMinArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTMaxArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int nBins = 10;
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	TCanvas *c = new TCanvas;
	c->Divide(2,5);
	
	float pTMin, pTMax;
	TString fName;
	TFile *canvasFile;
	TCanvas *cPlot;
	for(int i=0; i<nBins; i++)
	{
		pTMin = pTMinArr[i];	
		pTMax = pTMaxArr[i];
		
		fName = Form("eff_corrected_yield_at_BDT_cut/corr_yield_at_cut_pt_%.2f_%.2f.root", pTMin, pTMax);
		
		canvasFile = TFile::Open(fName);
	
		cPlot = (TCanvas*)canvasFile->Get("c1");
		
		cPlot->GetListOfPrimitives()->FindObject("title")->Clear();
		
		TLegend *legend = (TLegend*)cPlot->GetListOfPrimitives()->FindObject("TPave");
	
		((TLegendEntry*)legend->GetListOfPrimitives()->At(0))->SetLabel("Boosted Decision Tree");
		((TLegendEntry*)legend->GetListOfPrimitives()->At(1))->SetLabel("ALICE Standard");
		
		((TMultiGraph*)cPlot->GetListOfPrimitives()->FindObject(""))->GetYaxis()->SetTitle("Efficiency-Corrected Yield");
		((TMultiGraph*)cPlot->GetListOfPrimitives()->FindObject(""))->GetXaxis()->SetTitle("Boosted Decision Tree Cut Value");
		
		c->cd(i+1);
		
		cPlot->DrawClonePad();
		TPaveText *textBox = new TPaveText(.45,.90,.55,.96,"blNDC");
	
		textBox->AddText(Form("%.2f #leq #it{p}_{T} #leq  %.2f GeV/#it{c}",pTMin, pTMax));
	
		textBox->SetBorderSize(0);
		textBox->SetFillStyle(0);
		textBox->SetTextFont(42);
		textBox->SetTextSize(.04);
	
		textBox->Draw();		
	}
	
	c->SetCanvasSize(horSize, vertSize);
}

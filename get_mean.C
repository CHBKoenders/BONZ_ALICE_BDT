#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TROOT.h"

#include "TVectorF.h"

#include "TList.h"
#include "TFile.h"
#include "TMath.h"

//We will generate a data from a decaying exponential and a gaussian
//Then we will attempt to fit it like it's a real data

//first we need some global pointer variables
TH1D *trueSigDataMean;

//fits the data from the histogram, if the histogram is empty an error message will be displayed
void FittrueSigDataMean(bool silent = 0)
{
	if(!trueSigDataMean)
	{
		std::cout<<"Data has not been loaded yet! Terminating operation"<<std::endl;
		return;
	}
	
	//to find a good starting value for 'amplitude'
	//First we need to find the difference between the peak count and the background count at the peak
	//after that we scale that difference to get a good starting value
	float maxVal = trueSigDataMean->GetMaximum();
	//std::cout<< "maxVal is: " << maxVal << std::endl;
	float amp = maxVal*TMath::Sqrt(2.*TMath::Pi())*8e-3;
	//std::cout<< "amp is: " << amp << std::endl;
	
	//Fit	
	TF1 *fit_func = new TF1("signal", "[0]/([2]*TMath::Sqrt(2.*TMath::Pi()))*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",1.673,2.067);
	
	fit_func->SetParameters(amp,1.86,8e-3);
	fit_func->SetParNames("amplitude", "mean", "sigma");
	
	if(silent)
	{
		trueSigDataMean->Fit(fit_func,"QR0");
	}
	else
	{
		trueSigDataMean->Fit(fit_func,"R0");
	}	
}
	
void DrawFitMean()
{
	TCanvas *c = new TCanvas("c","Fit",1000,1000);
	
	//retrieve fit
	TF1 *fit = trueSigDataMean->GetFunction("signal");

	fit->SetTitle("fit");
	
	trueSigDataMean->SetTitle("Fitting result");	
	
	trueSigDataMean->Draw("E");
	trueSigDataMean->GetXaxis()->SetRangeUser(1.75, 2.2);
	fit->Draw("LSame");
	
	TLegend *legend = new TLegend(0.3,0.15);
	legend->AddEntry(trueSigDataMean, "Data","L");
	legend->AddEntry(fit, "fit", "L");	
	
	legend->Draw();
	
}

void ExtractTrueSignalMean(float pTMin, float pTMax)
{
	TString pTRange = Form("%.2f_%.2f",pTMin,pTMax);
	
	TFile *inputFile = TFile::Open("FullSimulation_LHC18a4a2CENT_Dplus_AnalysisResults.root");
	gDirectory->cd("PWGHF_TreeCreator");


	TList *inputList = (TList*)gDirectory->Get("coutputTreeHFTreeCreator");
	TTree *inputTree = (TTree*)inputList->First();
	
	TString outFileName = "MC_true_signal_pT_"+pTRange+".root";
	TString outDirName = "MC_true_signal/";
	TFile *outputFile = TFile::Open(outDirName+outFileName, "RECREATE");
	
	TH1F *trueSigHist;
	
	TCut sCut = "cand_type>0";
	TCut pTCut = Form("pt_cand >= %.2f && pt_cand <= %.2f",pTMin,pTMax);
	
	inputTree->Draw("inv_mass>>trueSigHist",sCut&&pTCut);
	
	trueSigHist = (TH1F*)gDirectory->Get("trueSigHist");
	
	trueSigHist->Write("", TObject::kOverwrite);
	
	outputFile->Close();
	inputFile->Close();
}


void GetMean(float pTMin = 0, float pTMax = 90, bool silent = 1)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	
	TString fName = "MC_true_signal_pT_"+pTRange+".root";
	TString dirName = "MC_true_signal/";
	
	if(gSystem->AccessPathName(dirName+fName))
	{
		//file does not exist
		ExtractTrueSignalMean(pTMin, pTMax);
	}

	TFile *dataFile = TFile::Open(dirName+fName, "UPDATE");
	trueSigDataMean = (TH1D*)dataFile->Get("trueSigHist");
	
	//fit
	FittrueSigDataMean(silent);

	//retrieve sigma and write it to the file in a vector
	TF1 *fit = trueSigDataMean->GetFunction("signal");
	TVectorF v(2);
	float mean = fit->GetParameter("mean");
	float meanErr = fit->GetParError(1);
	v[0] = mean;
	v[1] = meanErr;
	
	v.Write("mean", TObject::kOverwrite);
	
	dataFile->Close();
	
	std::cout << "mean is: " << mean << std::endl;
}


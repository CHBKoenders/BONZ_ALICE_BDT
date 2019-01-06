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
TH1D *trueSigData;

//fits the data from the histogram, if the histogram is empty an error message will be displayed
void FitTrueSigData(bool silent = 0)
{
	if(!trueSigData)
	{
		std::cout<<"Data has not been loaded yet! Terminating operation"<<std::endl;
		return;
	}
	
	//to find a good starting value for 'amplitude'
	//First we need to find the difference between the peak count and the background count at the peak
	//after that we scale that difference to get a good starting value
	float maxVal = trueSigData->GetMaximum();
	//std::cout<< "maxVal is: " << maxVal << std::endl;
	float amp = maxVal*TMath::Sqrt(2.*TMath::Pi())*8e-3;
	//std::cout<< "amp is: " << amp << std::endl;
	
	//Fit	
	TF1 *fit_func = new TF1("signal", "[0]/([2]*TMath::Sqrt(2.*TMath::Pi()))*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",1.673,2.067);
	
	fit_func->SetParameters(amp,1.86,8e-3);
	fit_func->SetParNames("amplitude", "mean", "sigma");
	
	if(silent)
	{
		trueSigData->Fit(fit_func,"QR0");
	}
	else
	{
		trueSigData->Fit(fit_func,"R0");
	}	
}
	
void DrawFit()
{
	TCanvas *c = new TCanvas("c","Fit",1000,1000);
	
	//retrieve fit
	TF1 *fit = trueSigData->GetFunction("signal");

	fit->SetTitle("fit");
	
	trueSigData->SetTitle("Fitting result");	
	
	trueSigData->Draw("E");
	trueSigData->GetXaxis()->SetRangeUser(1.75, 2.2);
	fit->Draw("LSame");
	
	TLegend *legend = new TLegend(0.3,0.15);
	legend->AddEntry(trueSigData, "Data","L");
	legend->AddEntry(fit, "fit", "L");	
	
	legend->Draw();
	
}

void ExtractTrueSignal(float pTMin, float pTMax)
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

//also gets mean now
void GetSigma(float pTMin = 0, float pTMax = 90, bool silent = 1)
{

	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	
	TString fName = "MC_true_signal_pT_"+pTRange+".root";
	TString dirName = "MC_true_signal/";
	
	if(gSystem->AccessPathName(dirName+fName))
	{
		//file does not exist
		ExtractTrueSignal(pTMin, pTMax);
	}

	TFile *dataFile = TFile::Open(dirName+fName, "UPDATE");
	trueSigData = (TH1D*)dataFile->Get("trueSigHist");
	
	//fit
	FitTrueSigData(silent);

	//retrieve sigma and write it to the file in a vector
	TF1 *fit = trueSigData->GetFunction("signal");
	TVectorF vs(2);
	TVectorF vm(2);
	float sigma = fit->GetParameter("sigma");
	float sigmaErr = fit->GetParError(2);
	float mean = fit->GetParameter("mean");
	float meanErr = fit->GetParError(1);
	vs[0] = sigma;
	vs[1] = sigmaErr;
	vm[0] = mean;
	vm[1] = meanErr;
	
	vs.Write("sigma", TObject::kOverwrite);
	vm.Write("mean", TObject::kOverwrite);
	
	dataFile->Close();
}


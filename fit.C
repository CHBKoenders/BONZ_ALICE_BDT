#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TROOT.h"

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLine.h"

#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TMath.h"
#include "TVectorF.h"


//first we need some global pointer variables
TH1D *data;
TGraphErrors *yieldGraph;

//other globals
double gpTMin;
double gpTMax;
double gCutValue;
bool gBDT;
bool gpTBinTrained;

double gMean;
double gMeanErr;
double gSigmaFit;
double gSigmaFitErr;

double gSignificance, gSignificanceErr;
double gRawSigYield, gRawSigYieldErr;
double gRawBkgYield, gRawBkgYieldErr;
double gRatio, gRatioErr;

const float massD0 = 1.865;
float gSigma;
float gSigmaErr;

//Fitting function for background only
double fit_bckgrnd(double *x, double *par)
{
	float xx = x[0];
	
	float xleft = massD0 - 5*gSigma;
	float xright = massD0 + 5*gSigma;	

	if (xleft<=xx && xx<=xright)
	{
		TF1::RejectPoint();
		return 0;
	}
	double f = par[0]+par[1]*xx+par[2]*xx*xx;

	return f;
}

//fits the data from the histogram, if the histogram is empty an error message will be displayed
void FitData(bool silent = 0)
{
	if(!data)
	{
		std::cout<<"Data has not been loaded yet! Terminating operation"<<std::endl;
		return;
	}

	//first we fit the background by leaving out the signal region
	TF1 *fit_func = new TF1("fit_background", fit_bckgrnd,1.673,2.067,3);
	
	fit_func->SetParameters(1,1,1);
	
	if(silent)
	{
		data->Fit(fit_func,"QR0");
	}
	else
	{
		data->Fit(fit_func,"R0");
	}
	//for checking
	//data->Draw("E");
	//fit_func->Draw("LSame");

	//now we get the parameters from the background fit
	double bkg_pars[3] = {fit_func->GetParameter(0),fit_func->GetParameter(1),fit_func->GetParameter(2)};

	//to find a good starting value for 'amplitude'
	//First we need to find the difference between the peak count and the background count at the peak
	//after that we scale that difference to get a good starting value
	float maxVal = data->GetMaximum();
	//std::cout<< "maxVal is: " << maxVal << std::endl;
	float regVal = bkg_pars[0] + bkg_pars[1]*massD0 + bkg_pars[2]*massD0*massD0;
	
	//std::cout<< "regVal is: " << regVal << std::endl;
	
	float amp = (maxVal-regVal)*TMath::Sqrt(2.*TMath::Pi())*8e-3;
	//std::cout<< "amp is: " << amp << std::endl;
	
	//Now we fit the entire signal	
	fit_func = new TF1("signal", "[0]/([2]*TMath::Sqrt(2.*TMath::Pi()))*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))+pol2(3)",1.673,2.067);
	
	
	
	fit_func->SetParameters(amp, massD0,gSigma,bkg_pars[0],bkg_pars[1],bkg_pars[2]);
	fit_func->SetParNames("amplitude", "mean", "sigma", "pol2 constant", "pol2 coefficient x","pol2 coefficient x^2");
	
	if(silent)
	{
		data->Fit(fit_func,"QR0");
	}
	else
	{
		data->Fit(fit_func,"R0");
	}
	
	gMean = fit_func->GetParameter(1);
	gMeanErr = fit_func->GetParError(1);
	gSigmaFit = fit_func->GetParameter(2);	
	gSigmaFitErr = fit_func->GetParError(2);
}

void ShowResponse(double pTMin, double pTMax, double cutValue)
{
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	
	TString fileName = "response_pT_"+pTRange+".root";
	
	if(gSystem->AccessPathName(fileName))
	{
		std::cout << fileName << " does not exist! Run Apply BDT to fix this." << std::endl;
		return;
	}
	
	TFile *responseFile = TFile::Open(fileName);
	
	TTree *responseTree = (TTree*)responseFile->Get("responseTree");
	
	float response;
	
	responseTree->SetBranchAddress("BDT_Response", &response);
	
	int nEntries = responseTree->GetEntries();
	
	TCanvas *c1 = new TCanvas;
	TH1D *responseHist = new TH1D("BDT response", Form("BDT response in pT range from %.2f to %.2f",pTMin, pTMax), 200, -1,1);
	
	for(int i=0; i<nEntries; i++)
	{
		responseTree->GetEntry(i);
		responseHist->Fill(response);
	}
	
	responseHist->Draw();
	
	c1->Update();
	
	TLine *cutLine = new TLine(cutValue, c1->GetUymin(), cutValue, c1->GetUymax());
	
	cutLine->SetLineColor(2);
	cutLine->Draw();
	
	TLegend *legend = new TLegend(0.3,0.15);
	legend->AddEntry(responseHist, "BDT response","LF");
	legend->AddEntry(cutLine, Form("Optimal cut (%.4f)",cutValue), "L");	
	
	legend->Draw();
}

void ShowFit()
{
	TString title;
	
	if(gBDT)
	{
		//BDT
		
		if(gpTBinTrained)
		{
			title = Form("Fitting result for BDT analysis in pT range from %.2f to %.2f",gpTMin, gpTMax);
		}
		else
		{
			title = Form("Fitting result for non pT bin trained BDT analysis in pT range from %.2f to %.2f",gpTMin, gpTMax);
		}
		
		
	}
	else
	{
		//Regular
		title = Form("Fitting result for regular ALICE analysis in pT range from %.2f to %.2f",gpTMin, gpTMax);
	}
	
	TCanvas *c = new TCanvas("c","Fit",1000,1000);
	
	//retrieve fit
	TF1 *fit = data->GetFunction("signal");

	//now we retrieve the background function	
	TF1 *background = new TF1("background", "pol2(0)",1.673,2.067);
	background->SetParameters(fit->GetParameter(3),fit->GetParameter(4),fit->GetParameter(5));

	background->SetTitle("background");
	background->SetLineColor(1);
	background->SetLineStyle(2);

	fit->SetTitle("fit");
	
	data->SetTitle(title);	
	
	data->Draw("E");
	data->GetXaxis()->SetRangeUser(1.7, 2.06);
	background->Draw("LSame");
	fit->Draw("LSame");
	
	TLegend *legend = new TLegend(0.3,0.15);
	legend->AddEntry(data, "Data","L");
	legend->AddEntry(background, "Background fit", "L");
	legend->AddEntry(fit, "Total fit", "L");	
	
	legend->Draw();
	
}

void GetYield(bool silent)
{
	
	//retrieve fit
	TF1 *fit = data->GetFunction("signal");

	double mean = fit->GetParameter(1);
	float sigma = fit->GetParameter(2);
	
	//std::cout<< "fitting sigma is: " << sigma << std::endl;
	
	double xLeft = mean-3*TMath::Abs(sigma);
	double xRight = mean+3*TMath::Abs(sigma);
	
	//std::cout << "Integration 3 sigma range is from " << xMin << " to " << xMax << std::endl;
	
	//now we retrieve the background function	
	TF1 *bkgFunc = new TF1("background", "pol2(0)",1.65,2.1);
	bkgFunc->SetParameters(fit->GetParameter(3),fit->GetParameter(4),fit->GetParameter(5));
	
	//and the signal function
	TF1 *sigFunc = new TF1("signal", "[0]/([2]*TMath::Sqrt(2.*TMath::Pi()))*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2]))", 1.65, 2.1);
	sigFunc->SetParameters(fit->GetParameter(0),fit->GetParameter(1),fit->GetParameter(2));
	
	//Background extraction
	int leftBand = data->FindBin(xLeft);
	int rightBand = data->FindBin(xRight);

	double leftBandInt = data->Integral(1,leftBand);
	double rightBandInt = data->Integral(rightBand,data->GetNbinsX());
	double bandInt = leftBandInt + rightBandInt;
	
	double sum2 = 0;
	for(int i = 1; i<=leftBand; i++)
	{
		sum2 += data->GetBinError(i)*data->GetBinError(i);
	}
	
	for(int i = rightBand; i<=data->GetNbinsX(); i++)
	{
		sum2 += data->GetBinError(i)*data->GetBinError(i);
	}
	
	double bandIntErr = TMath::Sqrt(sum2);
	
	gRawBkgYield = bkgFunc->Integral(xLeft, xRight)/(double)data->GetBinWidth(1);
	gRawBkgYieldErr = (bandIntErr/bandInt)*gRawBkgYield;
	
	//Signal extraction
	gRawSigYield = sigFunc->Integral(xLeft,xRight)/(double)data->GetBinWidth(1);
	double rawYield = fit->GetParameter(0)/data->GetBinWidth(1);
	double rawYieldErr = fit->GetParError(0)/data->GetBinWidth(1);
	gRawSigYieldErr = (rawYieldErr/rawYield)*gRawSigYield;
	
	//compute significance
	if(rawYield+gRawBkgYield <= 0.)
	{
		gSignificance = 1;
		gSignificanceErr = 0;
		
		return;
	}
	
	double sigErrSq = gRawSigYieldErr*gRawSigYieldErr;
	double bkgErrSq = gRawBkgYieldErr*gRawBkgYieldErr;
	double sigPlusBkg = gRawSigYield + gRawBkgYield;
	
	if(sigPlusBkg>0. && gRawSigYield>0.)
	{
		gSignificance = gRawSigYield/TMath::Sqrt(sigPlusBkg);
		gSignificanceErr = gSignificance*TMath::Sqrt((sigErrSq+bkgErrSq)/(4.*sigPlusBkg*sigPlusBkg) + (gRawBkgYield/sigPlusBkg)*sigErrSq/(gRawSigYield*gRawSigYield));
	}
	else
	{
		gSignificance = 0.;
		gSignificanceErr = 0.;
	}
	
	gRatio = gRawSigYield/gRawBkgYield;
	gRatioErr = TMath::Abs(gRatio)*TMath::Sqrt((gRawSigYieldErr/gRawSigYield)*(gRawSigYieldErr/gRawSigYield)+(gRawBkgYieldErr/gRawBkgYield)*(gRawBkgYieldErr/gRawBkgYield));
	
	//report the values
	if(not silent)
	{
		std::cout << "\nExtracted signal		: " << Form("%.1f +- %.2f", gRawSigYield, gRawSigYieldErr) << std::endl;
		std::cout << "Extracted background	: " << Form("%.1f +- %.2f", gRawBkgYield, gRawBkgYieldErr) << std::endl;
		
		std::cout << "Signal to background ratio: " << Form("%.2f", gRawSigYield/gRawBkgYield) << std::endl;
		
		std::cout << "Significance		: " << Form("%.1f +- %.2f", gSignificance, gSignificanceErr) << std::endl;
	}
	
	if(yieldGraph)
	{
		int nPoints = yieldGraph->GetN();
		yieldGraph->SetPoint(nPoints,gCutValue, gRawSigYield);
		yieldGraph->SetPointError(nPoints, 0, gRawSigYieldErr);
	}
}

void SaveFitTesting(bool atCut = 0)
{

	//create file
	TString pTRange = Form("%.2f_%.2f",gpTMin, gpTMax);
	TString fName;
	
	if(gBDT)
	{
		if(atCut)
		{
			TString cutValueS = Form("%.4f", gCutValue);
			fName = "Fitting_results/BDT_fit_pT_"+pTRange+"_at_cut_"+cutValueS+"_testing.root";
		}
		else
		{
			fName = "Fitting_results/BDT_fit_pT_"+pTRange+"_testing.root";
		}
	}
	
	else
	{
		fName = "Fitting_results/Std_fit_pT_"+pTRange+"_testing.root";	
	}
	TFile *fitSave = TFile::Open(fName,"RECREATE");
	
	//fill file
	//with: Data, total fit function, all yield results, BDT cutValue	
	data->SetDirectory(fitSave);
	data->Write();
	
	TVectorF sigYield(2); 
	sigYield[0] = gRawSigYield;
	sigYield[1] = gRawSigYieldErr;
	sigYield.Write("sigYield", TObject::kOverwrite);
	TVectorF bkgYield(2);
	bkgYield[0] = gRawBkgYield;
	bkgYield[1] = gRawBkgYieldErr;
	bkgYield.Write("bkgYield", TObject::kOverwrite);
	TVectorF significance(2);
	significance[0] = gSignificance;
	significance[1] = gSignificanceErr;
	significance.Write("significance", TObject::kOverwrite);
	TVectorF cutValue(1);
	cutValue[0] = gCutValue;
	cutValue.Write("cutValue", TObject::kOverwrite);
	
	delete fitSave;

}

void SaveFit(bool atCut = 0)
{

	//create file
	TString pTRange = Form("%.2f_%.2f",gpTMin, gpTMax);
	TString fName;
	
	if(gBDT)
	{
		if(atCut)
		{
			TString cutValueS = Form("%.4f", gCutValue);
			fName = "Fitting_results/BDT_fit_pT_"+pTRange+"_at_cut_"+cutValueS+".root";
		}
		else
		{
			fName = "Fitting_results/BDT_fit_pT_"+pTRange+".root";
		}
	}
	
	else
	{
		fName = "Fitting_results/Std_fit_pT_"+pTRange+".root";	
	}
	TFile *fitSave = TFile::Open(fName,"RECREATE");
	
	//fill file
	//with: Data, total fit function, all yield results, BDT cutValue	
	data->SetDirectory(fitSave);
	data->Write();
	
	TVectorF sigYield(2); 
	sigYield[0] = gRawSigYield;
	sigYield[1] = gRawSigYieldErr;
	sigYield.Write("sigYield", TObject::kOverwrite);
	TVectorF bkgYield(2);
	bkgYield[0] = gRawBkgYield;
	bkgYield[1] = gRawBkgYieldErr;
	bkgYield.Write("bkgYield", TObject::kOverwrite);
	TVectorF significance(2);
	significance[0] = gSignificance;
	significance[1] = gSignificanceErr;
	significance.Write("significance", TObject::kOverwrite);
	TVectorF cutValue(1);
	cutValue[0] = gCutValue;
	cutValue.Write("cutValue", TObject::kOverwrite);
	
	delete fitSave;

}


float Fit(float pTMin = 0, float pTMax = 90, bool BDT = 0, float cutValue = 0, bool silent = 0, bool pTBinTrained = 1)
{
	//Init
	gpTMin = pTMin;
	gpTMax = pTMax;
	gBDT = BDT;
	gpTBinTrained = pTBinTrained;
	gCutValue = cutValue;
	
	//to get sigma and mean
	TString fSigmaName = Form("MC_true_signal_pT_%.2f_%.2f.root", gpTMin,  gpTMax);
	TString dirSigma = "MC_true_signal/";
	
	if(gSystem->AccessPathName(dirSigma+fSigmaName))
	{
		//file does not exist yet
		gROOT->ProcessLine(".L get_sigma.C");
		gROOT->ProcessLine(Form("GetSigma(%.2f,%.2f)",gpTMin, gpTMax));	
	}
	
	TFile *sigmaFile = TFile::Open(dirSigma+fSigmaName);
	TVectorF *vSigma = (TVectorF*)sigmaFile->Get("sigma");
	gSigma = (*vSigma)[0];
	gSigmaErr = (*vSigma)[1];
	
	sigmaFile->Close();
	
	std::cout << "MC sigma is: " << gSigma << " (check)" << std::endl;
	
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString cutValueS = Form("%.4f", cutValue);
	
	TString fName;
	TString fDirName = "BDT_results/";
	TFile *dataFile;
	if(BDT)
	{
		if(pTBinTrained)
		{
			fName = "BDT_result_pT_"+pTRange+"_cut_at_"+cutValueS+".root";
		}
		else
		{
			fName = "BDT_result_pT_"+pTRange+"_cut_at_"+cutValueS+"_no_pT_train.root";
		}
			
		if(gSystem->AccessPathName(fDirName+fName))
		{
			//file does not exist
			
			std::cout << "The file " << fName << " does not exist! Running cut_BDT.C to fix this..." << std::endl;
			
			gROOT->ProcessLine(".L cut_BDT.C");
			gROOT->ProcessLine(Form("Cut(%.4f, %.2f, %.2f, %d)", cutValue, pTMin , pTMax, pTBinTrained));
		}
		
		dataFile = TFile::Open(fDirName+fName);
		data = (TH1D*)dataFile->Get("sigHist");
		
		if(!silent)
		{
		
			TCanvas *c1 = new TCanvas("c1", Form("invariant mass of signal and background classes in pT range from %.2f to %.2f",pTMin,pTMax), 1500, 1000);
	
			c1->Divide(1,2);
			
			c1->cd(2);
			((TH1D*)dataFile->Get("bkgHist"))->Draw();
			c1->cd(1);
			data->Draw();
		}
	
	}
	else
	{
		fName = "ALICE_regular_selection_pT_"+pTRange+".root";
		
		if(gSystem->AccessPathName(fName))
		{
			//file does not exist
			
			std::cout << "The file " << fName << " does not exist! Running regular_select.C to fix this..." << std::endl;
			
			gROOT->ProcessLine(".L regular_select.C");
			gROOT->ProcessLine(Form("Select(%.2f, %.2f)", pTMin , pTMax));
		}

		dataFile = TFile::Open(fName);
		data = (TH1D*)dataFile->Get("selectedHist");
	}
	
	if(data->GetEntries() == 0)
	{
		std::cout << "Nothing passes at this cut value! Aborting fit procedure..." << std::endl;
		return 0;
	}
	
	if(BDT)
	{
	
	}
	else
	{
		if(pTMin >= 6 && pTMin <= 10)
		{
			data->Rebin(2);
		}
		else if(pTMin == 12)
		{
			data->Rebin(4);
		}
		else if(pTMin >= 16)
		{
			data->Rebin(3);
		}
		
	}
	
	FitData(silent);
	GetYield(silent);
	
	//delete dataFile;
	
	return gSignificance;
}

float FitTestingData(float pTMin = 2, float pTMax = 90, float cutValue = 0, bool BDT = 0, bool silent = 1)
{
	//Init
	gpTMin = pTMin;
	gpTMax = pTMax;
	gBDT = BDT;
	gCutValue = cutValue;
	gpTBinTrained = 1;
	
	//to get sigma
	TString fSigmaName = Form("MC_true_signal_pT_%.2f_%.2f.root", gpTMin,  gpTMax);
	TString dirSigma = "MC_true_signal/";
	
	if(gSystem->AccessPathName(dirSigma+fSigmaName))
	{
		//file does not exist yet
		gROOT->ProcessLine(".L get_sigma.C");
		gROOT->ProcessLine(Form("GetSigma(%.2f,%.2f)",gpTMin, gpTMax));	
	}
	
	TFile *sigmaFile = TFile::Open(dirSigma+fSigmaName);
	TVectorF *vSigma = (TVectorF*)sigmaFile->Get("sigma");
	gSigma = (*vSigma)[0];
	
	sigmaFile->Close();
	
	//vars
	float BDT_value, invM;
	bool isSig, isselectedstd;
	
	//load in the testing data TTree
	
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return 0;
	
	}
	else
	{
		fInput = new TFile(fName);
	}
	
	//read Tree from file
	TTree *testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	//set branches
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_value);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	testingOutputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	
	//loop through tree and create de invM signal histogram
	TH1D *sigHist = new TH1D("sigHist", Form("Classified Signal when cutting at %.4f",cutValue), 100, 1.6,1.98);
	
	int nEntries = testingOutputTree->GetEntries();
	for(int i=0; i<nEntries; i++)
	{
		testingOutputTree->GetEntry(i);

		if(BDT && BDT_value >= cutValue)
		{
			sigHist->Fill(invM);
		}
		else if(!BDT && isselectedstd)
		{
			sigHist->Fill(invM);
		}
	}
	
	data = (TH1D*)sigHist->Clone("data");
	
		if(data->GetEntries() < 100)
	{
		std::cout << "Less than 100 entries at this cut value! Aborting fit procedure..." << std::endl;
		return 0;
	}
	
	FitData(silent);
	
	GetYield(silent);
	
	return gSignificance;
}

float ComputeSignificance(float cutValue)
{
	std::cout << Form("Computing Significance at cut value %.4f", cutValue) << std::endl;
	float significance = Fit(gpTMin, gpTMax, 1, cutValue,1, gpTBinTrained);
	
	//close file that was opened last
	((TFile*)gROOT->GetListOfFiles()->Last())->Close();
	
	return significance;

}

float ComputeSignificanceTesting(float cutValue)
{
	std::cout << Form("Computing Significance at cut value %.4f", cutValue) << std::endl;
	float significance = FitTestingData(gpTMin, gpTMax, cutValue,1);
	
	//close file that was opened last
	((TFile*)gROOT->GetListOfFiles()->Last())->Close();
	
	return significance;

}

void PlotSignificance(float pTMin,  float pTMax)
{
	gpTMin = pTMin;
	gpTMax = pTMax;
	gpTBinTrained = 1;
	
	TF1 *significanceFunc = new TF1("significance", "ComputeSignificance(x)", -1,1);
	
	TCanvas *c = new TCanvas;
	
	significanceFunc->SetTitle(Form("%.2f #leq #it{p}_{T} #leq %.2f",pTMin, pTMax));
	
	significanceFunc->GetYaxis()->SetTitle("Significance");
	significanceFunc->GetYaxis()->CenterTitle();
	
	significanceFunc->GetXaxis()->SetTitle("BDT cut value");
	significanceFunc->GetXaxis()->CenterTitle();
	
	float optCutVal = significanceFunc->GetMaximumX(-1, 1);
	
	significanceFunc->Draw();
	
	int horSize = 2000;
	int vertSize = TMath::Sqrt(2)*horSize;
	
	c->SetCanvasSize(.5*horSize, .25*vertSize);
	
	TLine *cutLine = new TLine(optCutVal, c->GetUymin(), optCutVal, c->GetUymax());
	
	cutLine->SetLineColor(8);
	cutLine->Draw();
}

double FindOptCut(float pTMin, float pTMax ,bool pTBinTrained = 1,float cutMin = -1, float cutMax = 1)
{
	gpTMin = pTMin;
	gpTMax = pTMax;
	gpTBinTrained = pTBinTrained;
	
	TF1 *significanceFunc = new TF1("significance", "ComputeSignificance(x)", -1,1);
	
	float maxSignificance = significanceFunc->GetMaximum(cutMin, cutMax);
	float optCutVal = significanceFunc->GetMaximumX(cutMin, cutMax);
	
	std::cout << "The maximum significance is " << Form("%4.4f",maxSignificance) << " when cutting at " << Form("%4.4f", optCutVal) << std::endl;
	
	Fit( pTMin, pTMax,1,optCutVal, 0, pTBinTrained);
	
	return optCutVal;
}

void FindOptCutTesting(float pTMin, float pTMax, float cutMin = -1, float cutMax = 1)
{
	gpTMin = pTMin;
	gpTMax = pTMax;
	gpTBinTrained = 1;
	
	TF1 *significanceFunc = new TF1("significance", "ComputeSignificanceTesting(x)", -1,1);
	
	float maxSignificance = significanceFunc->GetMaximum(cutMin, cutMax);
	float optCutVal = significanceFunc->GetMaximumX(cutMin, cutMax);
	
	std::cout << "The maximum significance is " << Form("%4.4f",maxSignificance) << " when cutting at " << Form("%4.4f", optCutVal) << std::endl;
	
	FitTestingData( pTMin, pTMax,optCutVal, 1,0);
}

void FillYieldGraph(float pTMin, float pTMax, int  nPoints,float minCut = -1, float maxCut = 1)
{
	yieldGraph = new TGraphErrors();
	TGraph *trueYieldGraph = new TGraph();
	
	float cutRange = maxCut - minCut;
	float stepSize = cutRange/(nPoints);
	
	float cutValue;
	
	//get testing data
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	else
	{
		fInput = new TFile(fName);
	}
	
	//read Tree from file
	TTree *testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	float BDT_response, invM;
	bool isSig;
	
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_response);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	
	
	float fitReturn;
	
	for(int i=0; i<nPoints; i++)
	{
		gCutValue = minCut + i*stepSize;
		fitReturn = FitTestingData(pTMin, pTMax, gCutValue,1,1);
		
		//std::cout << fitReturn << std::endl;
		
		if(fitReturn == 0)
		{
			//do nothing
		}
		else
		{
			//retrieve fit
			TF1 *fit = data->GetFunction("signal");

			double mean = fit->GetParameter(1);
			float sigma = fit->GetParameter(2);
			
			double minMass = mean-3*TMath::Abs(sigma);
			double maxMass = mean+3*TMath::Abs(sigma);
			
			//get the true values
			int nEntries = testingOutputTree->GetEntries();
			int nSig = 0;
			int nBkg = 0;
			for(int i=0; i<nEntries; i++)
			{
				testingOutputTree->GetEntry(i);
				
				if(BDT_response >= gCutValue && invM>=minMass && invM<=maxMass)
				{
					nSig += isSig;
					nBkg += !isSig;
				}
			
			}
			
			//std::cout<< "nSig is " << nSig << std::endl;
			
			trueYieldGraph->SetPoint(i,gCutValue,nSig);
		}
	}
	
	TCanvas *c1 = new TCanvas("c1", Form("extracted yield vs true yield in pT range from %.2f to %.2f",pTMin, pTMax),1024,720);
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(yieldGraph);
	combinedGraph->Add(trueYieldGraph);
	
	combinedGraph->SetTitle(Form("extracted yield vs true yield in pT range from %.2f to %.2f",pTMin, pTMax));
	
	yieldGraph->SetTitle("Extracted Signal yield");
	yieldGraph->SetLineColor(2);
	yieldGraph->SetMarkerStyle(7);
	yieldGraph->SetMarkerColor(2);
	yieldGraph->SetFillColor(2);
	yieldGraph->SetFillStyle(3001);
	trueYieldGraph->SetTitle("True signal yield");
	trueYieldGraph->SetLineColor(4);
	trueYieldGraph->SetMarkerColor(4);
	trueYieldGraph->SetMarkerStyle(7);
	
	combinedGraph->Draw("APL 3");
	
	combinedGraph->GetXaxis()->SetTitle("BDT cut value");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Extracted Yield Count");
	combinedGraph->GetYaxis()->CenterTitle();
	
	
	gPad->BuildLegend();
}

void FillSignificanceGraph()
{
	TGraphErrors *significanceGraphBDT = new TGraphErrors;
	TGraphErrors *significanceGraphStd = new TGraphErrors;
	
	int pTLowArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int numBins = 10;
	
	for(int i=0; i<numBins; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		double optCut = FindOptCut(pTLow, pTHigh);		
		Fit(pTLow,pTHigh,1,optCut,1);
		
		significanceGraphBDT->SetPoint(i,(pTLow+pTHigh)/2.,gSignificance);
		significanceGraphBDT->SetPointError(i,(pTHigh-pTLow)/2.,gSignificanceErr);
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
		
		Fit(pTLow, pTHigh);
		
		significanceGraphStd->SetPoint(i, (pTLow+pTHigh)/2., gSignificance);
		significanceGraphStd->SetPointError(i, (pTHigh-pTLow)/2., gSignificanceErr);
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
	}
	
	TCanvas *c1 = new TCanvas("c1", "Significance in pT intervals",1024,720);
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(significanceGraphBDT);
	combinedGraph->Add(significanceGraphStd);
	
	combinedGraph->SetTitle("Significance in #it{p}_{T} intervals");
	
	significanceGraphBDT->SetTitle("BDT Significance");
	significanceGraphBDT->SetLineColor(2);
	//significanceGraphBDT->SetFillColor(2);
	//significanceGraphBDT->SetFillStyle(3001);
	significanceGraphBDT->SetMarkerStyle(8);
	significanceGraphStd->SetTitle("Std Significance");
	significanceGraphStd->SetLineColor(4);
	//significanceGraphStd->SetFillColor(4);
	//significanceGraphStd->SetFillStyle(3001);
	significanceGraphStd->SetMarkerStyle(22);
	
	combinedGraph->Draw("AP");
	
	combinedGraph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c}");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Significance");
	combinedGraph->GetYaxis()->CenterTitle();	
	
	gPad->BuildLegend();

}

void FillRatioGraph()
{
	TGraphErrors *ratioGraphBDT = new TGraphErrors;
	TGraphErrors *ratioGraphStd = new TGraphErrors;
	
	int pTLowArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int numBins = 10;
	
	for(int i=0; i<numBins; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		double optCut = FindOptCut(pTLow, pTHigh);
				
		Fit(pTLow,pTHigh,1,optCut,1);
		
		ratioGraphBDT->SetPoint(i,(pTLow+pTHigh)/2.,gRatio);
		ratioGraphBDT->SetPointError(i,(pTHigh-pTLow)/2.,gRatioErr);

		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
		
		Fit(pTLow, pTHigh);
		
		ratioGraphStd->SetPoint(i, (pTLow+pTHigh)/2., gRatio);
		ratioGraphStd->SetPointError(i, (pTHigh-pTLow)/2., gRatioErr);
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
	}
	
	TCanvas *c1 = new TCanvas("c1", "Signal to Background Ratio in #it{p}_{T} intervals",1024,720);
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(ratioGraphBDT);
	combinedGraph->Add(ratioGraphStd);
	
	combinedGraph->SetTitle("Signal to Background Ratio in #it{p}_{T} intervals");
	
	ratioGraphBDT->SetTitle("BDT Ratio");
	ratioGraphBDT->SetLineColor(2);
	//ratioGraphBDT->SetFillColor(2);
	//ratioGraphBDT->SetFillStyle(3001);
	ratioGraphBDT->SetMarkerStyle(8);
	ratioGraphStd->SetTitle("Std Ratio");
	ratioGraphStd->SetLineColor(4);
	//ratioGraphStd->SetFillColor(4);
	//ratioGraphStd->SetFillStyle(3001);
	ratioGraphStd->SetMarkerStyle(22);
	
	combinedGraph->Draw("AP");
	
	combinedGraph->GetXaxis()->SetTitle("#{p}_{T} (GeV/#it{c})");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Ratio");
	combinedGraph->GetYaxis()->CenterTitle();

	gPad->BuildLegend();

}

void FillMeanGraph()
{
	TGraphErrors *meanGraphBDT = new TGraphErrors;
	TGraphErrors *meanGraphStd = new TGraphErrors;
	TGraphErrors *meanGraphMC  = new TGraphErrors;
	
	int pTLowArr[] = 	{2,3,4,5,6,7,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,10,12,16,24};
	
	int numBins = 10;
	
	for(int i=0; i<numBins; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		double optCut = FindOptCut(pTLow, pTHigh);
				
		Fit(pTLow,pTHigh,1,optCut,1);
		
		meanGraphBDT->SetPoint(i,(pTLow+pTHigh)/2.,gMean);
		meanGraphBDT->SetPointError(i,(pTHigh-pTLow)/2.,gMeanErr);

		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
		
		Fit(pTLow, pTHigh);
		
		meanGraphStd->SetPoint(i, (pTLow+pTHigh)/2., gMean);
		meanGraphStd->SetPointError(i, (pTHigh-pTLow)/2., gMeanErr);
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
		
		//retrieve MC mean
		TString fMeanName = Form("MC_true_signal_pT_%d.00_%d.00.root", pTLow,  pTHigh);
		TString dirMean = "MC_true_signal/";
		
		if(gSystem->AccessPathName(dirMean+fMeanName))
		{
			//file does not exist yet
			std::cout << "calculating MC mean..." << std::endl;
			gROOT->ProcessLine(".L get_mean.C");
			gROOT->ProcessLine(Form("GetMean(%d,%d)",pTLow, pTHigh));	
		}
		
		TFile *meanFile = TFile::Open(dirMean+fMeanName);
		TVectorF *vMean = (TVectorF*)meanFile->Get("mean");
		
		if(!vMean)
		{
			delete meanFile;
			delete vMean;
			//null
			std::cout << "calculating MC mean..." << std::endl;
			gROOT->ProcessLine(".L get_mean.C");
			gROOT->ProcessLine(Form("GetMean(%d,%d)",pTLow, pTHigh));	
		}
		
		meanFile = TFile::Open(dirMean+fMeanName);
		vMean = (TVectorF*)meanFile->Get("mean");
		
		float meanMC = (*vMean)[0];
		float meanMCErr = (*vMean)[1];
		
		meanGraphMC->SetPoint(i, (pTLow+pTHigh)/2., meanMC);
		meanGraphMC->SetPointError(i, (pTHigh-pTLow)/2., meanMCErr);
		
		meanFile->Close();
	}
	
	TCanvas *c1 = new TCanvas("c1", "Mean of Gaussian fit in pT intervals",1024,720);
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(meanGraphMC,"PL 3");
	combinedGraph->Add(meanGraphBDT);
	combinedGraph->Add(meanGraphStd);
	
	combinedGraph->SetTitle("Mean of Gaussian fit in pT intervals");
	
	meanGraphBDT->SetTitle("BDT Mean");
	meanGraphBDT->SetLineColor(2);
	//ratioGraphBDT->SetFillColor(2);
	//ratioGraphBDT->SetFillStyle(3001);
	meanGraphBDT->SetMarkerStyle(8);
	meanGraphStd->SetTitle("Std Mean");
	meanGraphStd->SetLineColor(4);
	//ratioGraphStd->SetFillColor(4);
	//ratioGraphStd->SetFillStyle(3001);
	meanGraphStd->SetMarkerStyle(22);
	meanGraphMC->SetTitle("MC Mean");
	meanGraphMC->SetLineColor(3);
	meanGraphMC->SetFillColor(3);
	meanGraphMC->SetFillStyle(3001);
	
	
	combinedGraph->Draw("AP");
	
	combinedGraph->GetXaxis()->SetTitle("pT");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Mean");
	combinedGraph->GetYaxis()->CenterTitle();

	gPad->BuildLegend();
}

void FillSigmaGraph()
{
	TGraphErrors *sigmaGraphBDT = new TGraphErrors;
	TGraphErrors *sigmaGraphStd = new TGraphErrors;
	TGraphErrors *sigmaGraphMC  = new TGraphErrors;
	
	int pTLowArr[] = 	{2,3,4,5,6,6,7,8 ,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,8,10,12,12,16,24};
	
	for(int i=0; i<12; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		double optCut = FindOptCut(pTLow, pTHigh);
				
		Fit(pTLow,pTHigh,1,optCut,1);
		
		sigmaGraphBDT->SetPoint(i,(pTLow+pTHigh)/2.,gSigmaFit);
		sigmaGraphBDT->SetPointError(i,(pTHigh-pTLow)/2.,gSigmaFitErr);

		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
		
		Fit(pTLow, pTHigh);
		
		sigmaGraphStd->SetPoint(i, (pTLow+pTHigh)/2., gSigmaFit);
		sigmaGraphStd->SetPointError(i, (pTHigh-pTLow)/2., gSigmaFitErr);
		
		sigmaGraphMC->SetPoint(i, (pTLow+pTHigh)/2., gSigma);
		sigmaGraphMC->SetPointError(i, (pTHigh-pTLow)/2., gSigmaErr);
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
	}
	
	TCanvas *c1 = new TCanvas("c1", "Sigma of Gaussian fit in pT intervals",1024,720);
	
	TMultiGraph *combinedGraph = new TMultiGraph();
	combinedGraph->Add(sigmaGraphMC, "PL 3");
	combinedGraph->Add(sigmaGraphBDT);
	combinedGraph->Add(sigmaGraphStd);
	
	
	combinedGraph->SetTitle("Sigma of Gaussian fit in pT intervals");
	
	sigmaGraphBDT->SetTitle("BDT Sigma");
	sigmaGraphBDT->SetLineColor(2);
	//ratioGraphBDT->SetFillColor(2);
	//ratioGraphBDT->SetFillStyle(3001);
	sigmaGraphBDT->SetMarkerStyle(8);
	sigmaGraphStd->SetTitle("Std Sigma");
	sigmaGraphStd->SetLineColor(4);
	//ratioGraphStd->SetFillColor(4);
	//ratioGraphStd->SetFillStyle(3001);
	sigmaGraphStd->SetMarkerStyle(22);
	sigmaGraphMC->SetTitle("MC Sigma");
	sigmaGraphMC->SetLineColor(3);
	sigmaGraphMC->SetFillColor(3);
	sigmaGraphMC->SetFillStyle(3001);
	
	combinedGraph->Draw("AP");
	
	combinedGraph->GetXaxis()->SetTitle("pT");
	combinedGraph->GetXaxis()->CenterTitle();
	combinedGraph->GetYaxis()->SetTitle("Sigma");
	combinedGraph->GetYaxis()->CenterTitle();

	gPad->BuildLegend();
}

void RunBDTCutRange(float pTMin, float pTMax, int nCuts, bool testing = 0)
{
	float cutStep = 2./(nCuts+1);
	
	float cutValue;
	for(int i=1; i<=nCuts; i++)
	{
		cutValue = -1 + i*cutStep;
		
		if(testing)
		{
			FitTestingData(pTMin, pTMax, cutValue, 1, 1);
			
			SaveFitTesting(1);
		}
		else
		{
			Fit(pTMin, pTMax, 1, cutValue,1);
			
			SaveFit(1);
		}
		
		//close file that was opened last -> Gives Segmentation Violation at i=2 when cut_BDT is used...
		//((TFile*)gROOT->GetListOfFiles()->Last())->;
	}
}

void RunAllpTBinsTesting()
{
	int pTLowArr[] = 	{2,3,4,5,6,6,7,8 ,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,8,10,12,12,16,24};
	
	for(int i=0; i<12; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		//Std
		FitTestingData(pTLow, pTHigh);
		SaveFitTesting();
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
		
		//BDT
		FindOptCutTesting(pTLow, pTHigh);
		SaveFitTesting();
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
	}
}

void RunAllpTBins()
{
	int pTLowArr[] = 	{2,3,4,5,6,6,7,8 ,8 ,10,12,16};
	int pTHighArr[] = 	{3,4,5,6,7,8,8,10,12,12,16,24};
	
	for(int i=0; i<12; i++)
	{
		int pTLow = pTLowArr[i];
		int pTHigh = pTHighArr[i];
		
		//Std
		Fit(pTLow, pTHigh);
		SaveFit();
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
		
		//BDT
		FindOptCut(pTLow, pTHigh);
		SaveFit();
		
		//close file that was opened last
		((TFile*)gROOT->GetListOfFiles()->Last())->Close();
	}
}

void checkStd(float pTMin, float pTMax)
{
	FitTestingData(pTMin, pTMax);
	std::cout << "When fitting the testing data we get: \n" << std::endl;
	GetYield(0);
	
	//get testing data
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	else
	{
		fInput = new TFile(fName);
	}
	
	//read Tree from file
	TTree *testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	float BDT_response, invM;
	bool isSig,isselectedstd;
	
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_response);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	testingOutputTree->SetBranchAddress("isselectedstd", &isselectedstd);
	
	//retrieve fit
	TF1 *fit = data->GetFunction("signal");
	double mean = fit->GetParameter(1);
	float sigma = fit->GetParameter(2);
	
	double minMass = mean-3*TMath::Abs(sigma);
	double maxMass = mean+3*TMath::Abs(sigma);
	
	//get the true values
	int nEntries = testingOutputTree->GetEntries();
	int nSig = 0;
	int nBkg = 0;
	TH1D *trueSigHist = new TH1D("trueSigHist", "True signal from standard selection", 100, 1.6,2.2);
	for(int i=0; i<nEntries; i++)
	{
		testingOutputTree->GetEntry(i);
		
		if(isselectedstd && invM>=minMass && invM<=maxMass)
		{
			nSig += isSig;
			nBkg += !isSig;
			
			trueSigHist->Fill(invM);
		}
	}
		
	std::cout <<  "\n\nThe true values are:\n" << std::endl;
	std::cout << "Extracted signal		: " << Form("%d", nSig) << std::endl;
	std::cout << "Extracted background	: " << Form("%d", nBkg) << std::endl;
	std::cout << "Signal to background ratio: " << Form("%.2f", (float)nSig/nBkg) << std::endl;		
	std::cout << "S/Srt(S+B)		: " << Form("%.1f", nSig/TMath::Sqrt(nSig+nBkg)) << std::endl;	
	
	trueSigHist->Draw();
}

void checkBDT(float pTMin, float pTMax, float cutValue)
{
	FitTestingData(pTMin,pTMax,cutValue,1,1);
	std::cout << "When fitting the testing data we get: \n" << std::endl;
	GetYield(0);
	
	//retrieve fit
	TF1 *fit = data->GetFunction("signal");

	double mean = fit->GetParameter(1);
	float sigma = fit->GetParameter(2);
	
	double minMass = mean-3*TMath::Abs(sigma);
	double maxMass = mean+3*TMath::Abs(sigma);
	
	//get testing data
	TString pTRange = Form("%.2f_%.2f", pTMin, pTMax);
	TString fName = "testing_output_"+pTRange+".root";
	TFile* fInput;
	if(gSystem->AccessPathName(fName))	
	{
		//file does not exist
		std::cout << "the file " << fName << " does not exist! Aborting operation..." << std::endl;
		return;
	
	}
	else
	{
		fInput = new TFile(fName);
	}
	
	//read Tree from file
	TTree *testingOutputTree = (TTree*)fInput->Get("testOutputTree");
	
	float BDT_response, invM;
	bool isSig;
	
	testingOutputTree->SetBranchAddress("BDT_value", &BDT_response);
	testingOutputTree->SetBranchAddress("invM", &invM);
	testingOutputTree->SetBranchAddress("isSig", &isSig);
	
	int nEntries = testingOutputTree->GetEntries();
	int nSig = 0;
	int nBkg = 0;
	TH1D *trueSigHist = new TH1D("trueSigHist", "True signal from BDT selection", 100, 1.6,2.2);
	for(int i=0; i<nEntries; i++)
	{
		testingOutputTree->GetEntry(i);
		
		if(BDT_response >= cutValue && invM>=minMass && invM<=maxMass)
		{
			nSig += isSig;
			nBkg += !isSig;
			
			trueSigHist->Fill(invM);
		}
	
	}
	
	std::cout << Form("\n\nThe true values are:\n") << std::endl;
	std::cout << "Extracted signal		: " << Form("%d", nSig) << std::endl;
	std::cout << "Extracted background	: " << Form("%d", nBkg) << std::endl;
	std::cout << "Signal to background ratio: " << Form("%.2f", (float)nSig/nBkg) << std::endl;		
	std::cout << "S/Srt(S+B)		: " << Form("%.1f", nSig/TMath::Sqrt(nSig+nBkg)) << std::endl;
	
	trueSigHist->Draw();
}

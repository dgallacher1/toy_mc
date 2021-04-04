#include "InverseCDF.h"

ClassImp(InverseCDF)

InverseCDF::InverseCDF()
{
  Init();
}

InverseCDF::~InverseCDF()
{
  Clear();
}


TGraph* InverseCDF::GetInverseHisto(TH1D* hPDF)
{
  hPDF->Scale(1.0/hPDF->Integral()); // Renormalize incase histogram isn't scaled
  TH1D *hCDF = (TH1D*) hPDF->GetCumulative();//Convert to CDF
  hCDF->SetBinContent(1,0.0); // Ensure histogram starts at 0
  TGraph *gCDF = new TGraph(hCDF);// Make TGraph from histogram
  Int_t num_points = gCDF->GetN();
  Double_t *x = gCDF->GetX();
  Double_t *y = gCDF->GetY();

  TGraph *gInvCDF = new TGraph(num_points,y,x); // Invert axes

  //Cleanup
  delete gCDF;
  delete hCDF;

  return gInvCDF;
}

TGraph* InverseCDF::GetInverse(vector<Double_t> x_vals, vector<Double_t> y_vals)
{
  //Generate PDF first
  TH1D *hPDF = GetPDF(x_vals,y_vals);
  //Use GetInverseHisto
  TGraph *gInvCDF = GetInverseHisto(hPDF);
  return gInvCDF;
}


/// Reset class
void InverseCDF::Clear(Option_t *option)
{
  
}

///Initialize parameters
void InverseCDF::Init()
{

}

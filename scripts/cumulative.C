


#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"
#include <math.h>       /* isnan, sqrt */


void cumulative(){
  fStyle->DEAPStyle();
  gStyle->SetNdivisions(508,"xyz");
  gSystem->Load("/home/davidg/pyrene_toymc/src/libtoymc.so");

  ///Get ToyMC input spectra
  ToyMC *tmc = new ToyMC(1234);
  TFile *DataFile = new TFile("../dat/spectras.root","READ");
  tmc->LoadPDFs(DataFile);
  tmc->LoadFunctions();


  int nbins = 2000;
  int high = 10000.;

  TF1 *fLArPS = tmc->GetLArPS();
  TH1D *hLArPS = new TH1D("hLArPS","",nbins,0,high);
  hLArPS->Add(fLArPS,1.0);
  hLArPS->Scale(1.0/hLArPS->Integral());
  hLArPS->SetLineColor(kBlue);

  TF1 *fTPBPS = tmc->GetTPBPS();
  TH1D *hTPBPS = new TH1D("hTPBPS","",nbins,0,high);
  hTPBPS->Add(fTPBPS,1.0);
  hTPBPS->Scale(1.0/hTPBPS->Integral());
  hTPBPS->SetLineColor(kMagenta);


  TH1D* hLArTPB = new TH1D("hLArTPB","LAr convolved with TPB PS",nbins,0,high);
  for(int i=1;i<nbins;i++){
    double val = hTPBPS->GetBinContent(i)*hLArPS->GetBinContent(i);
    hLArTPB->SetBinContent(i,val);
  }
  hLArTPB->SetLineColor(kBlue);

  TF1 *fPyrenePS = tmc->GetPyrenePS();
  TH1D *hPyrenePS = new TH1D("hPyrenePS","NR LAr Pulseshape Comparison;time[ns];Intensity(A.U)",nbins,0,high);
  hPyrenePS->Add(fPyrenePS,1.0);
  hPyrenePS->Scale(1.0/hPyrenePS->Integral());
  hPyrenePS->SetLineColor(kGreen);

  TH1D *hPyreneCumu = (TH1D*)hPyrenePS->GetCumulative();
  hPyreneCumu->Scale(1.0/hPyreneCumu->GetBinContent(hPyreneCumu->GetMaximumBin()));
  hPyreneCumu->SetLineColor(kGreen);

  TH1D *hLArTPBCumu = (TH1D*)hLArTPB->GetCumulative();
  hLArTPBCumu->Scale(1.0/hLArTPBCumu->GetBinContent(hLArTPBCumu->GetMaximumBin()));
  hLArTPBCumu->SetLineColor(kBlue);

  TH1D *hSubLArCumu = (TH1D*)hPyreneCumu->Clone();
  hSubLArCumu->Add(hLArTPBCumu,-1);
  hSubLArCumu->SetLineColor(kBlack);
  hSubLArCumu->SetTitle("Cumulative difference between Pyrene and NR LAr Pulseshapes;time[ns];Difference(AU)");


  TCanvas *c1 = new TCanvas("c1","NR PS Comparison",1200,1000);
  c1->Divide(3,1);
  c1->cd(1);
  hPyrenePS->Draw();
  hLArTPB->Draw("same");
  gPad->SetLogy();

  c1->cd(2);
  hPyreneCumu->Draw();
  hPyreneCumu->GetXaxis()->SetRangeUser(0,1000);
  hLArTPBCumu->Draw("same");
  //gPad->SetLogx();
  //gPad->SetLogy();


  c1->cd(3);
  hSubLArCumu->Draw();


  ///Repeat for ER parameters
  Double_t paramlar[] = {0.225,2.2,0.045,75.5,0.73,1445.0};
  tmc->SetLArParameters(paramlar);
  TF1 *fLArERPS = tmc->GetLArPS();

  TH1D *hLArERPS = new TH1D("hLArERPS","",nbins,0,high);
  hLArERPS->Add(fLArERPS,1.0);
  hLArERPS->Scale(1.0/hLArERPS->Integral());
  hLArERPS->SetLineColor(kBlue);


  TH1D* hLArERTPB = new TH1D("hLArERTPB","ER LAr convolved with TPB PS",nbins,0,high);
  for(int i=1;i<nbins;i++){
    double val = hTPBPS->GetBinContent(i)*hLArERPS->GetBinContent(i);
    hLArERTPB->SetBinContent(i,val);
  }
  hLArERTPB->SetLineColor(kBlue);

  TH1D *hLArERTPBCumu = (TH1D*)hLArERTPB->GetCumulative();
  hLArERTPBCumu->Scale(1.0/hLArERTPBCumu->GetBinContent(hLArERTPBCumu->GetMaximumBin()));
  hLArERTPBCumu->SetLineColor(kBlue);

  TH1D *hSubLArCumuER = (TH1D*)hPyreneCumu->Clone();
  hSubLArCumuER->Add(hLArERTPBCumu,-1);
  hSubLArCumuER->SetLineColor(kBlack);
  hSubLArCumuER->SetTitle("Cumulative difference between Pyrene and ER LAr Pulseshapes;time[ns];Difference(AU)");



  TCanvas *c2 = new TCanvas("c2","ER PS Comparison",1200,1000);
  c2->Divide(3,1);
  c2->cd(1);
  hLArERTPB->Draw();
  hPyrenePS->Draw("same");
  gPad->SetLogy();

  c2->cd(2);
  hPyreneCumu->Draw();
  hPyreneCumu->GetXaxis()->SetRangeUser(0,5000);
  hLArERTPBCumu->Draw("same");
  //gPad->SetLogx();
  //gPad->SetLogy();


  c2->cd(3);
  hSubLArCumuER->Draw();



}

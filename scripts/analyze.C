//Analyze some outputs

#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"
#include <math.h>       /* isnan, sqrt */

TH1D *CalculatefPyrene(TTree *t,double winLow, double winHigh);


void analyze(string filename="output_0.root"){


  string path = "../output/";
  path+=filename;
  cout << "Reading.. " << path <<endl;
  TFile *fileinNA = new TFile(path.c_str());
  string pathNR = "../output/output_1.root";
  TFile *fileinNR = new TFile(pathNR.c_str());
  string pathER = "../output/output_2.root";
  TFile *fileinER = new TFile(pathER.c_str());


  //Read in tree
  TTree *tree = (TTree*) fileinNA->Get("T");
  TTree *treeNR = (TTree*) fileinNR->Get("T");
  TTree *treeER = (TTree*) fileinER->Get("T");


  int nbins = 250;
  int high = 12000.;

  TCanvas *c1 = new TCanvas("c1","Canvas",1200,800);
  c1->Divide(2,1);
  c1->cd(1);
  tree->Draw(Form("times>>hp(%i,-100,%i)",nbins,high));
  TH1D *hFraction = (TH1D*)gPad->GetPrimitive("hp");
  hFraction->SetTitle("PMT Hit time distribution;Time[ns];Count");
  hFraction->Scale(1.0/hFraction->Integral());

  //Read in LAr NR's
  treeNR->Draw(Form("times>>hp1(%i,-100,%i)",nbins,high));
  TH1D *hFractionNR = (TH1D*)gPad->GetPrimitive("hp1");
  hFractionNR->SetTitle("PMT Hit time distribution;Time[ns];Count");
  hFractionNR->Scale(1.0/hFractionNR->Integral());
  hFractionNR->SetLineColor(kBlue);

  //Read in LAr NR's
  treeER->Draw(Form("times>>hp2(%i,-100,%i)",nbins,high));
  TH1D *hFractionER = (TH1D*)gPad->GetPrimitive("hp2");
  hFractionER->SetTitle("PMT Hit time distribution;Time[ns];Count");
  hFractionER->Scale(1.0/hFractionER->Integral());
  hFractionER->SetLineColor(kRed);

  c1->cd(1);
  hFractionNR->Draw();
  hFraction->Draw("same");
  hFractionER->Draw("same");

  TLegend *leg1 = new TLegend(0.3,0.7,0.9,0.9);
  leg1->AddEntry(hFraction,"Fast MC Pyrene coated Neck simulation","l");
  leg1->AddEntry(hFractionNR,"Fast MC NR LAr bulk simulation","l");
  leg1->AddEntry(hFractionER,"Fast MC ER LAr bulk simulation","l");
  leg1->SetFillStyle(3001);
  leg1->Draw("same");
  gPad->SetLogy();


  ///Get ToyMC input spectra
  ToyMC *tmc = new ToyMC(1234);
  TFile *DataFile = new TFile("../dat/spectras.root","READ");
  tmc->LoadPDFs(DataFile);
  tmc->LoadFunctions();

  TF1 *fLArPS = tmc->GetLArPS();
  TH1D *hLArPS = new TH1D("hLArPS","",nbins,0,high);
  hLArPS->Add(fLArPS,1.0);
  hLArPS->Scale(1.0/hLArPS->GetBinContent(hLArPS->GetMaximumBin()));
  hLArPS->SetLineColor(kBlue);

  TF1 *fTPBPS = tmc->GetTPBPS();
  TH1D *hTPBPS = new TH1D("hTPBPS","",nbins,0,high);
  hTPBPS->Add(fTPBPS,1.0);
  hTPBPS->Scale(1.0/hTPBPS->GetBinContent(hTPBPS->GetMaximumBin()));
  hTPBPS->SetLineColor(kMagenta);
  TF1 *fPyrenePS = tmc->GetPyrenePS();
  TH1D *hPyrenePS = new TH1D("hPyrenePS","",nbins,0,high);
  hPyrenePS->Add(fPyrenePS,1.0);
  hPyrenePS->Scale(1.0/hPyrenePS->GetBinContent(hPyrenePS->GetMaximumBin()));
  hPyrenePS->SetLineColor(kGreen);

  TF1 *fPMTPS = tmc->GetPMTPS();
  TH1D *hPMTPS = new TH1D("hPMTPS","",nbins,0,high);
  hPMTPS->Add(fPMTPS,1.0);
  hPMTPS->Scale(1.0/hPMTPS->GetBinContent(hPMTPS->GetMaximumBin()));
  hPMTPS->SetLineColor(kRed);

  c1->cd(2);
  gPad->SetLogy();
  hLArPS->Draw();
  hLArPS->SetTitle("Pulseshape comparison;Time[ns];Intensity(AU)");
  hPyrenePS->Draw("same");
  hTPBPS->Draw("same");
  hPMTPS->Draw("same");

  TLegend *leg = new TLegend(0.5,0.6,0.9,0.9);
  leg->AddEntry(hLArPS,"NR LAr Pulseshape","lp");
  leg->AddEntry(hTPBPS,"TPB WLS Pulseshape","lp");
  leg->AddEntry(hPyrenePS,"Pyrene 2 exp Pulseshape","lp");
  leg->AddEntry(hPMTPS,"PMT AP Pulseshape","lp");
  leg->SetFillStyle(3001);
  leg->Draw("same");

  c1->Print(Form("../plots/Pulseshapes_%d_%d.pdf",0,high));


  ////////////////
  //Read trees
  ////////////////

  Double_t windowLowNR = -20.0; //NR rejection
  Double_t windowHighNR = 40.0;
  Double_t windowLowER = 60.0; // ER rejection
  Double_t windowHighER = 400.0;

  TH1D *hFpyreneNANR = CalculatefPyrene(tree,windowLowNR,windowHighNR);
  TH1D *hFpyreneNR = CalculatefPyrene(treeNR,windowLowNR,windowHighNR);

  TH1D *hFpyreneNAER = CalculatefPyrene(tree,windowLowER,windowHighER);
  TH1D *hFpyreneER = CalculatefPyrene(treeER,windowLowER,windowHighER);



  //Draw PDFs
  TCanvas *c3 = new TCanvas("c2","fPyrene",1200,800);
  c3->Divide(2,1);
  c3->cd(1);
  hFpyreneNANR->Scale(1.0/hFpyreneNANR->Integral());
  hFpyreneNR->SetLineColor(kBlue);
  hFpyreneNR->Scale(1.0/hFpyreneNR->Integral());
  hFpyreneNR->Draw();
  hFpyreneNANR->SetTitle("FpyreneNR Optimization;FpyreneNR;Intensity(AU)");
  hFpyreneNANR->Draw("same");

  //Find overlap
  TH1D *hOv = (TH1D*) hFpyreneNANR->Clone();
  for(int i =1; i < hFpyreneNANR->GetNbinsX() ; i++){
    Double_t binContNA = hFpyreneNANR->GetBinContent(i);
    Double_t binContNR = hFpyreneNR->GetBinContent(i);
    Double_t min = (binContNA < binContNR) ? binContNA : binContNR;
    hOv->SetBinContent(i,min);
  }

  c3->cd(1);
  hOv->SetLineColor(kMagenta);
  hOv->SetFillStyle(3001);
  hOv->SetFillColorAlpha(kMagenta,0.5);
  hOv->Draw("same");
  cout <<"Overlap = " << hOv->Integral() << endl;


  TLegend *leg31 = new TLegend(0.5,0.7,0.9,0.9);
  leg31->AddEntry(hFpyreneNANR,"Neck-alpha with PyrenePS","lp");
  leg31->AddEntry(hFpyreneNR,"Ar40 NR","lp");
  leg31->AddEntry(hOv,"Overlap","lpf");
  leg31->SetFillStyle(3001);
  leg31->Draw("same");

  c3->cd(2);
  hFpyreneNAER->Scale(1.0/hFpyreneNAER->Integral());
  hFpyreneER->SetLineColor(kRed);
  hFpyreneER->Scale(1.0/hFpyreneER->Integral());
  hFpyreneER->Draw();
  hFpyreneNAER->SetTitle("FpyreneER Optimization;FpyreneER;Intensity(AU)");
  hFpyreneNAER->Draw("same");

  TLegend *leg32 = new TLegend(0.5,0.7,0.9,0.9);
  leg32->AddEntry(hFpyreneNAER,"Neck-alpha with PyrenePS","lp");
  leg32->AddEntry(hFpyreneER,"Ar39 Betas","lp");
  leg32->SetFillStyle(3001);
  leg32->Draw("same");

  c3->Print("../plots/Fpyrene.pdf");

  //Draw nice plots
  fStyle->DEAPStyle();
  //->drawDEAP();
  TCanvas *c4 = new TCanvas("c4","Paper Plot",1200,800);
  c4->Divide(2,1);
  c4->cd(1);
  hFraction->SetLineColor(fStyle->Color(0));
  hFraction->SetTitle("Fast MC Pulseshape Comparison;Photon Arrival Time [ns];Intensity [AU]");
  hFraction->Draw();
  hFraction->GetYaxis()->SetTitleOffset(2);
  hFraction->GetYaxis()->SetRangeUser(1e-5,1);
  hFractionNR->SetLineColor(fStyle->Color(1));
  hFractionNR->Draw("same");
  gPad->SetLogy();
  TLegend *leg4 = new TLegend(0.3,0.6,0.9,0.9);
  leg4->AddEntry(hFraction,"Pyrene+PS Inlet Alphas Pulseshape","lp");
  leg4->AddEntry(hFractionNR,"Bulk LAr NR-like Pulseshape","lp");
  leg4->SetFillStyle(3001);
  leg4->Draw("same");

  c4->cd(2);
  hFpyreneNANR->SetLineColor(fStyle->Color(0));
  hFpyreneNR->SetTitle("Fast MC Optimized PSD Comparison;fPPS; Intensity [AU]");
  hFpyreneNR->SetLineColor(fStyle->Color(1));
  hFpyreneNR->Draw();
  hFpyreneNR->GetYaxis()->SetTitleOffset(2);
  hFpyreneNANR->Draw("same");
  TLegend *leg5 = new TLegend(0.15,0.6,0.6,0.9);
  leg5->AddEntry(hFpyreneNANR,"Pyrene+PS Inlet Alphas","lp");
  leg5->AddEntry(hFpyreneNR,"Bulk LAr NRs","lp");
  leg5->SetFillStyle(3001);
  leg5->Draw("same");

  c4->Print("../plots/nice_plot.pdf");

  TFile *filePlots = new TFile("../plots/analysis_plots.root","RECREATE");
  filePlots->Write();
  filePlots->Close();
}

//Find fPyrene
TH1D *CalculatefPyrene(TTree *t, double winLow, double winHigh){

  Double_t windowEnd = 10000.0;
  TH1D *hFpyrene = new TH1D("",Form("Pyrene PSD [%3.2f,%3.2f ns];fPyrene;Count",winLow,winHigh),50,0,1);

  vector<Double_t> *times = 0;
  Int_t numHits;
  Double_t edep;
  Double_t shadowFraction;
  t->SetBranchAddress("times",&times);
  t->SetBranchAddress("numHits",&numHits);
  t->SetBranchAddress("edep",&edep);
  t->SetBranchAddress("shadowFraction",&shadowFraction);

  for(int iTrial = 0; iTrial< t->GetEntries();iTrial++){
    t->GetEntry(iTrial);
    Double_t fPyrene = 0.0;
    Double_t integral = 0.0;
    if(!times)continue;
    const int num = (*times).size();
    if(num==0)continue;
    for(int i=0;i<numHits-1;i++){
      if((*times)[i] > winLow && (*times)[i] < winHigh) fPyrene++;
      if((*times)[i] > winLow && (*times)[i] < windowEnd) integral++;
    }
    fPyrene = fPyrene/integral;
    hFpyrene->Fill(fPyrene);
  }
  return hFpyrene;
}

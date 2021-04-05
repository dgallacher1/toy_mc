//Analyze some outputs

#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"


void analyze(){

  TFile *filein = new TFile("../output/output.root");

  //Read in tree
  TTree *tree = (TTree*) filein->Get("T");

  TCanvas *c1 = new TCanvas("c1","Canvas",1200,800);
  c1->Divide(2,1);
  c1->cd(1);
  tree->Draw("times>>hp(100,0,15000)");
  TH1D *hFraction = (TH1D*)gPad->GetPrimitive("hp");
  hFraction->SetTitle("PMT Hit time distribution;Time[ns];Count");
  tree->SetLineColor(kRed);
  tree->Draw("times","shadowFraction>0.8","same");
  gPad->SetLogy();

  TH1D *hCumuFrac = hFraction->GetCumulative();
  c1->cd(2);
  hCumuFrac->Scale(1.0/hCumuFrac->Integral());
  hCumuFrac->SetTitle("Cumulative Timing distribution;Time[ns];Integral(AU)");
  hCumuFrac->Draw();

  vector<Double_t> *times=0;
  Int_t numHits;
  Double_t edep;
  Double_t shadowFraction;

  tree->SetBranchAddress("times",&times);
  tree->SetBranchAddress("numHits",&numHits);
  tree->SetBranchAddress("edep",&edep);
  tree->SetBranchAddress("shadowFraction",&shadowFraction);

  TH1D *hFpyrene = new TH1D("hFpryene","Pyrene PSD [0,100 ns];fPyrene;Count",50,0,1);

  Double_t windowLow = 0.0;
  Double_t windowHigh = 100.0;
  Double_t windowEnd = 10000.0;
  for(int iTrial = 0;iTrial<tree->GetEntries();iTrial++){
    tree->GetEntry(iTrial);
    Double_t fPyrene = 0.0;
    TH1D hTemp("","hTemp",1000,0,windowEnd);
    const int num = (*times).size();
    vector<Double_t> weights(num,1);
    hTemp.FillN(num,&(*times)[0],&weights[0]);
    fPyrene = hTemp.Integral(hTemp.FindBin(windowLow),hTemp.FindBin(windowHigh));
    Double_t integral = hTemp.Integral(hTemp.FindBin(windowLow),hTemp.FindBin(windowEnd));
    fPyrene/=integral;
    hFpyrene->Fill(fPyrene);
  }


  //Scan over windows
  TH2D *hFpWindow = new TH2D("hFpScan","fPyrene Scan;Window Length[ns];Window Start[ns]",50,10,1000,50,0,1000);


  TCanvas *c2 = new TCanvas("c2","fPyrene",1200,800);
  hFpyrene->Draw();


}

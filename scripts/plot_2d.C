/*
Make 2D plots showing energy vs PSD parameter for Fast MC of Neck-alphas with Pyrene coating and bulk NR's and ER's
Author: David Gallacher
Date: April 2021
*/

#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"
#include <math.h>       /* isnan, sqrt */


TH1D *CalculatefPyrenevsEnergy(TTree *t,double winLow, double winHigh);

void plot_2d(string filename="output_0.root"){
    fStyle->DEAPStyle();

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

    Double_t windowLowNR = -20.0; //NR rejection
    Double_t windowHighNR = 40.0;
    Double_t windowLowER = 60.0; // ER rejection
    Double_t windowHighER = 400.0;

    //NR discrimination
    TH2D *hFpENANR = CalculatefPyrenevsEnergy(tree,windowLowNR,windowHighNR);
    TH2D *hFpENR = CalculatefPyrenevsEnergy(treeNR,windowLowNR,windowHighNR);

    TCanvas *c1 = new TCanvas("c1","Canvas",1200,800);
    hFpENR->SetMarkerStyle(32);
    hFpENR->SetMarkerColor(fStyle->Color(1));
    hFpENANR->SetMarkerColor(fStyle->Color(0));
    hFpENANR->SetMarkerStyle(4);

    hFpENR->Draw();
    hFpENANR->Draw("same");


}

//Find fPyrene vs num Hits
TH2D *CalculatefPyrenevsEnergy(TTree *t, double winLow, double winHigh){

  Double_t windowEnd = 10000.0;
  TH2D *hFpyrene = new TH2D("",Form("PE vs Pyrene PSD [%3.2f,%3.2f ns];PE;fPyreneNR",winLow,winHigh),100,0,5000,50,0,1);

  vector<Double_t> *times = 0;
  Int_t numHits;
  Int_t numPh;
  Int_t numAP;
  Double_t edep;
  Double_t shadowFraction;
  t->SetBranchAddress("times",&times);
  t->SetBranchAddress("numHits",&numHits);
  t->SetBranchAddress("numPhotons",&numPh);
  t->SetBranchAddress("numAP",&numAP);
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
    Int_t PE = numHits;
    hFpyrene->Fill(PE,fPyrene);
    //cout << "Num photons = "<<numPh << " PDE = " << (numHits-numAP)/double(numPh)<<endl;
  }
  return hFpyrene;
}

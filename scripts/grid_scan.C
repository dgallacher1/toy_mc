//Grid scan for maximum separation between pyrene coated neck-alphas and LAr events
//Author: David Gallacher

#include <vector>


//Need to load shared library file to access class methods
gSystem->Load("/home/davidg/pyrene_toymc/src/libtoymc.so");

TH1D *CalculatefPyrene(TTree *t,double winLow, double winHigh);

void grid_scan(string filename="output_0.root"){

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
  //TTree *treeNR = (TTree*) fileinER->Get("T");

  Double_t windowLow[] = {0.0,25.0,50.0,100.0,150.0,200.0,250.0,300.0};
  Double_t windowLength[] = {25.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0};
  int nbins = 10;

  //FindOverlap(tree,treeNR,0.0,25.0)

  TCanvas *c1 = new TCanvas("c1","Canvas",1200,800);
  TH2D *hScan = new TH2D("hScan","fPyreneNR window scan;Window start [ns];Window Length [ns];Overlap (AU)",nbins,0.0,200.0,nbins,25.0,300.0);

  for(int iBinX=1;iBinX<nbins+1;iBinX++){
    Double_t winStart = hScan->GetXaxis()->GetBinLowEdge(iBinX);//windowLow[iBinX];
    for(int iBinY=1;iBinY<nbins+1;iBinY++){
      Double_t winStop = winStart+ hScan->GetYaxis()->GetBinLowEdge(iBinY);//windowLow[iBinX];
      //cout << "win start = "<< winStart <<endl;
      //cout << "win stop = "<< winStop <<endl;
      Double_t ov = FindOverlap(tree,treeNR,winStart,winStop);
      //cout << "ov = "<<ov <<endl;
      hScan->SetBinContent(iBinX,iBinY,ov);
    }//Stop times
  }//Start times

  fStyle->DEAPStyle();
  fStyle->chooseDSPalette(1);
  hScan->Draw("colz");
  gPad->SetLogz();
}

//Calculate overlap of h1 with h2, assuming same binning
Double_t FindOverlap(TTree *T1, TTree *T2,Double_t winLow,Double_t winHigh){

    Double_t windowEnd = 10000.0;
    TH1D hFpyrene("",Form("Pyrene PSD [%3.2f,%3.2f ns];fPyrene;Count",winLow,winHigh),50,0,1);
    vector<Double_t> *times = 0;
    Int_t numHits;
    Double_t edep;
    Double_t shadowFraction;
    T1->SetBranchAddress("times",&times);
    T1->SetBranchAddress("numHits",&numHits);
    T1->SetBranchAddress("edep",&edep);
    T1->SetBranchAddress("shadowFraction",&shadowFraction);

    int ntrials = 5000;

    for(int iTrial = 0;iTrial<ntrials;iTrial++){
      T1->GetEntry(iTrial);
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
      hFpyrene.Fill(fPyrene);
    }

    TH1D hNR("",Form("Pyrene PSD [%3.2f,%3.2f ns];fPyrene;Count",winLow,winHigh),50,0,1);
    times = 0;
    numHits = 0;
    edep = 0.0;
    shadowFraction = 0.0;
    T2->SetBranchAddress("times",&times);
    T2->SetBranchAddress("numHits",&numHits);
    T2->SetBranchAddress("edep",&edep);
    T2->SetBranchAddress("shadowFraction",&shadowFraction);

    for(int iTrial = 0;iTrial<ntrials;iTrial++){
      T2->GetEntry(iTrial);
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
      hNR.Fill(fPyrene);
    }

    //Scale to unit integral
    hFpyrene.Scale(1.0/hFpyrene.Integral());
    hNR.Scale(1.0/hNR.Integral());

    TH1D *hOv = (TH1D*) hFpyrene.Clone();
    for(int i =1;i<hFpyrene.GetNbinsX();i++){
      Double_t binContNA = hFpyrene.GetBinContent(i);
      Double_t binContNR = hNR.GetBinContent(i);
      Double_t min = (binContNA < binContNR) ? binContNA : binContNR;
      hOv->SetBinContent(i,min);
    }
    return hOv->Integral();


    //////////
    //Double_t integral = hFpyrene.Integral();

    // Double_t mean1 = hFpyrene.GetMean();
    // Double_t mean2 = hNR.GetMean();
    //
    // Int_t binStart = 0;
    // //Integral from left to right, so need to find which is where
    // Bool_t swap = mean1>mean2;
    // cout << "Swap = "<< swap <<endl;
    // if(swap) binStart = hFpyrene.FindFirstBinAbove(0.0);
    // else binStart = hNR.FindFirstBinAbove(0.0);
    // Double_t overlap = 0.0;
    // if(swap) overlap = hFpyrene.Integral(binStart,hNR.FindLastBinAbove(0.0));
    // else overlap = hFpyrene.Integral(binStart,hFpyrene.GetNbinsX());
    // if(overlap == 0)
    // {
    //   cout << "ov =0 "<<endl;
    //   return 0.0;
    // }
    // //What percentage of h1 overlaps with h2
    // return overlap/integral;
}

//Eg for use:
//root -l draw_plot_paper.C'("updated_ap_1_type_0.root")'


#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"
#include <math.h>       /* isnan, sqrt */

TH2D *CalculatefPyrenevsEnergy(TTree *t,double winLow, double winHigh);
gSystem->Load("libtoymc");

void leakage_vs_pe(string filename="output_0.root")
{
    //DEAP Style settings
    fStyle->DEAPStyle();
    //File in
    string path = "../output/";
    path+=filename;
    cout << "Reading.. " << path <<endl;
    TFile *fileinNA = new TFile(path.c_str());
    string pathNR = path.substr(0,path.find(".root")-1);
    pathNR+="1.root";
    TFile *fileinNR = new TFile(pathNR.c_str());

    //Read in tree
    TTree *tree = (TTree*) fileinNA->Get("T");
    TTree *treeNR = (TTree*) fileinNR->Get("T");

    Double_t windowLowNR = -20.0; //NR rejection
    Double_t windowHighNR = 40.0;

    //NR discrimination
    TH2D *hFpENANR = CalculatefPyrenevsEnergy(tree,windowLowNR,windowHighNR);
    TH2D *hFpENR = CalculatefPyrenevsEnergy(treeNR,windowLowNR,windowHighNR);

    //hFpENANR->Draw("colz");

    TH1D *hLeak = new TH1D("hLeak","90% UL Leakage Fraction;PE;LF",100,0,1);
    Int_t nEvts = tree->GetEntries();
    int binstart = 1;
    //Find upper limit
    if (!gROOT->GetClass("TFeldmanCousins")) gSystem->Load("libPhysics");
    TFeldmanCousins f;

    TH1D *hFpyreneNR = (TH1D*) hFpENR->ProjectionY("nr_py",1,hFpENR->GetNbinsY());
    Int_t bin_start = hFpyreneNR->FindFirstBinAbove(0);
    cout << "bin start = "<<bin_start<<endl;


    //Calculate leakage
    //TH1D *hFpyreneNR = NULL;
    TH1D *hFPPS = NULL;
    for(int iBin =1; iBin < hFpENANR->GetNbinsX();iBin++){
      // Get the projections, increasing one bin at a time, last leakage is then the total leakage
      //hFPPS = (TH1D*) hFpENANR->ProjectionY(Form("_na_%i",iBin),1,iBin+1);
      //printf("NR Distribution Start = %f \n",nr_start);
      //Integral of PPS distribution from start of NR distribution
      Double_t pps_leakage_n = 0;
      //pps_leakage_n = hFPPS->Integral(hFPPS->FindBin(bin_start),hFPPS->GetNbinsX());

      Int_t entries = 0;
      for(int iY=1;iY<hFpENANR->GetNbinsY();iY++){
        //Assume 100% acceptance for LAr NRs
        //Find first bin above 0 for NR distribution in this bin
        for(int iX = 1;iX < iBin;iX++){
          if(iY > bin_start) pps_leakage_n += hFpENANR->GetBinContent(iX,iY);
        }
        entries += hFpENANR->GetBinContent(iBin,iY);
      }

      //printf("Leakage of PPS into NR region = %f \n",pps_leakage_n);
      //90% Upper limit on number of events with 0 background
      Double_t upper_lim = f.CalculateUpperLimit(pps_leakage_n,0);
      if(upper_lim<pps_leakage_n) upper_lim = pps_leakage_n + sqrt(pps_leakage_n); //sqrt(N) errors for large N
      //printf("Upper limit PPS into NR region = %f \n",upper_lim);
      //printf("num entries = %i\n",entries);
      //Double_t leak = (entries==0) ? 0.0 : upper_lim/float(entries);//Leakage fraction
      Double_t leak = upper_lim/float(nEvts);
      //printf("90% UL on LF for PPS into NR region = %f \n",leak);
      hLeak->SetBinContent(iBin,leak);
    }

    //Formatting
    hLeak->SetTitle("90% UL Leakage fraction for PPS Inlet Alphas; Fraction of total PE; Leakage Fraction");
    hLeak->SetLineColor(fStyle->Color(0));
    hLeak->Draw();
    gPad->SetLogy();
    gPad->Print("../plots/leakage_vs_pe.pdf");


}

//Find fPyrene vs num Hits
TH2D *CalculatefPyrenevsEnergy(TTree *t, double winLow, double winHigh){

  Double_t windowEnd = 10000.0;
  Double_t rescale =1.0; //40000/2000
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
    Int_t PE = numHits*rescale;
    hFpyrene->Fill(PE,fPyrene);
    //cout << "Num photons = "<<numPh << " PDE = " << (numHits-numAP)/double(numPh)<<endl;
  }
  return hFpyrene;
}

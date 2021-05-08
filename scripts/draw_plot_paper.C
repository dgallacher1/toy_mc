//Nice plots
//Eg for use:
//root -l draw_plot_paper.C'("updated_ap_1_type_0.root")'


#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"
#include <math.h>       /* isnan, sqrt */

TH1D *CalculatefPyrene(TTree *t,double winLow, double winHigh,double PE_cut);
gSystem->Load("libtoymc");

void draw_plot_paper(string filename="output_0.root")
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

    //Pulseshape parameters
    int nbins = 250;
    int high = 10000.;

    TCanvas *c1 = new TCanvas("c1","Plots",1200,800);
    tree->Draw(Form("times>>hp(%i,-100,%i)",nbins,high));
    TH1D *hFraction = (TH1D*)gPad->GetPrimitive("hp");
    hFraction->SetTitle("Hit time distribution;Time[ns];Intensity (A.U)");
    hFraction->Scale(1.0/hFraction->Integral());

    //Read in LAr NR's
    treeNR->Draw(Form("times>>hp1(%i,-100,%i)",nbins,high));
    TH1D *hFractionNR = (TH1D*)gPad->GetPrimitive("hp1");
    hFractionNR->SetTitle("Hit time distribution;Time[ns];Intensity (A.U)");
    hFractionNR->Scale(1.0/hFractionNR->Integral());
    hFractionNR->SetLineColor(fStyle->Color(0));

    /////////
    ///Formatting
    /////////

    //Match colors from diagram in paper
    int o_index = 9999;
    TColor *c_orange = new TColor(o_index,0.902,0.404,0.173);
    hFraction->SetLineColor(o_index);
    int g_index = 9998;
    TColor *c_green = new TColor(g_index,33./255,122./255,11./255);
    hFractionNR->SetLineColor(g_index);

    gPad->SetLogy();
    hFractionNR->Draw();
    hFraction->SetTitle("Fast MC Pulseshape Comparison;Photon Arrival Time [ns];Intensity [AU]");
    hFractionNR->GetYaxis()->SetTitleOffset(1.5);
    hFraction->Draw("same");
    TLegend *leg = new TLegend(0.5,0.6,0.9,0.9);
    leg->AddEntry(hFraction,"PPS Inlet Alphas Pulseshape","lp");
    leg->AddEntry(hFractionNR,"Bulk LAr NR-like Pulseshape","lp");
    leg->SetFillStyle(3001);
    leg->Draw("same");
    c1->Print("../plots/Pulseshape_draw_plot_paper.pdf");

    TCanvas *c2 = new TCanvas("c2","Plots",1200,800);

    //Get Fpyrene
    Double_t windowLowNR = -20.0; //NR rejection
    Double_t windowHighNR = 40.0;
    Double_t windowLowER = 60.0; // ER rejection
    Double_t windowHighER = 400.0;
    Double_t PE_cut_lowest_quartile = 1400; //from Plot_2d.C, With 2000 Ph/MeV LY

    TH1D *hFpyreneNANR = CalculatefPyrene(tree,windowLowNR,windowHighNR,PE_cut_lowest_quartile);
    TH1D *hFpyreneNR = CalculatefPyrene(treeNR,windowLowNR,windowHighNR,1e5);//No PE Cut on NR events,(Full acceptance)


    hFpyreneNANR->Scale(1.0/hFpyreneNANR->Integral());
    hFpyreneNR->Scale(1.0/hFpyreneNR->Integral());
    hFpyreneNR->Draw();
    hFpyreneNANR->Draw("same");

    //Find overlap
    TH1D *hOv = (TH1D*) hFpyreneNANR->Clone();
    for(int i =1; i < hFpyreneNANR->GetNbinsX() ; i++){
      Double_t binContNA = hFpyreneNANR->GetBinContent(i);
      Double_t binContNR = hFpyreneNR->GetBinContent(i);
      Double_t min = (binContNA < binContNR) ? binContNA : binContNR;
      hOv->SetBinContent(i,min);
    }
    //hOv->Draw("same");
    //gPad->SetLogy();

    //Formatting
    hOv->SetLineColor(kMagenta);
    hOv->SetFillStyle(3001);
    hOv->SetFillColorAlpha(kMagenta,0.5);
    cout <<"Overlap = " << hOv->Integral() << endl;
    hFpyreneNR->SetTitle("Fast MC Optimized PSD Comparison; FPyrene; Intensity [AU]");
    hFpyreneNANR->SetLineColor(o_index);
    hFpyreneNR->SetLineColor(g_index);
    hFpyreneNR->GetYaxis()->SetTitleOffset(1.5);
    TLegend *legFp = new TLegend(0.15,0.7,0.4,0.9);
    legFp->AddEntry(hFpyreneNANR,"PPS Inlet Alphas","lp");
    legFp->AddEntry(hFpyreneNR,"Bulk LAr NRs","lp");
    //legFp->AddEntry(hOv,"Overlap","lp");
    legFp->SetFillStyle(3001);
    legFp->Draw("same");

    c2->Print("../plots/Fpyrene_draw_plot_paper.pdf");

}

//Find FPyrene
TH1D *CalculatefPyrene(TTree *t, double winLow, double winHigh,double PE_cut){

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
    if(numHits > PE_cut)continue; //introduce a PE Cut
    for(int i=0;i<numHits-1;i++){
      if((*times)[i] > winLow && (*times)[i] < winHigh) fPyrene++;
      if((*times)[i] > winLow && (*times)[i] < windowEnd) integral++;
    }
    fPyrene = fPyrene/integral;
    hFpyrene->Fill(fPyrene);
  }
  return hFpyrene;
}

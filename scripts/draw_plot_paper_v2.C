//Nice plots
//Eg for use:
//root -l draw_plot_paper.C'("updated_ap_1_type_0.root")'


#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"
#include <math.h>       /* isnan, sqrt */

TH1D *CalculatefPyrene(TTree *t,double winLow, double winHigh,double PE_cut);
gSystem->Load("libtoymc");

void draw_plot_paper_v2(string filename="default_PLQY_ap_1_type_0.root")
{
    //DEAP Style settings
    fStyle->DEAPStyle();

    //File in
    string path = "../output/";
    path+=filename;
    cout << "Reading.. " << path <<endl;
    //TFile *fileinNADef = new TFile("../output/update_plqy_ap_1_type_0.root");
    TFile *fileinNADef = new TFile("../output/pulse_shape_ap_1_type_0.root");

    TFile *fileinNALow = new TFile("../output/update_plqy_low_ap_1_type_0.root");
    TFile *fileinNAHigh = new TFile("../output/high_PLQY_ap_1_type_0.root");

    string pathNR = path.substr(0,path.find(".root")-1);
    pathNR+="1.root";
    TFile *fileinNR = new TFile(pathNR.c_str());

    //Read in tree
    TTree *treeDef = (TTree*) fileinNADef->Get("T");
    TTree *treeLow = (TTree*) fileinNALow->Get("T");
    TTree *treeHigh = (TTree*) fileinNAHigh->Get("T");

    TTree *treeNR = (TTree*) fileinNR->Get("T");

    //Pulseshape parameters
    int nbins = 250;
    int high = 6000.;

    //Match colors from diagram in paper
    // int o_index = 9999;
    // TColor *c_orange = new TColor(o_index,0.902,0.404,0.173);
    // int g_index = 9998;
    // TColor *c_green = new TColor(g_index,33./255,122./255,11./255);

    /*
    TCanvas *c1 = new TCanvas("c1","Plots",1200,800);
    treeDef->Draw(Form("times>>hp(%i,-100,%i)",nbins,high));
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

    hFraction->SetLineColor(kRed);
    hFractionNR->SetLineColor(kBlue);
    hFractionNR->SetLineStyle(3);

    gPad->SetLogy();
    gPad->SetLogx();
    hFractionNR->Draw();
//    hFractionNR->SetTitle("Toy MC Pulseshape Comparison;Photon Arrival Time [ns];Intensity [AU]");
    hFractionNR->SetTitle(";Photon Arrival Time [ns];Intensity [AU]");
    hFractionNR->GetYaxis()->SetTitleSize(0.08);
    hFractionNR->GetXaxis()->SetTitleSize(0.08);
    hFractionNR->GetYaxis()->SetTitleOffset(1.5);
    hFraction->Draw("same");
    TLegend *leg = new TLegend(0.4,0.6,0.9,0.9);
    leg->AddEntry(hFraction,"PPS Inlet Alphas Pulseshape","lp");
    leg->AddEntry(hFractionNR,"Bulk LAr NR-like Pulseshape","lp");
    leg->SetFillStyle(3001);
    leg->Draw("same");
    c1->Print("../plots/Pulseshape_draw_plot_paper.pdf");

    */
    TCanvas *c2 = new TCanvas("c2","Plots",1200,800);

    //Get Fpyrene
    Double_t windowLowNR = -20.0; //NR rejection
    Double_t windowHighNR = 40.0;
    Double_t windowLowER = 60.0; // ER rejection
    Double_t windowHighER = 400.0;
    //Double_t PE_cut_lowest_quartile = 8845; //from Plot_2d.C, With 20000 Ph/MeV LY
    Double_t PE_cut_lowest_quartile = 15166; //from Plot_2d.C, With 20000 Ph/MeV LY
    Double_t PE_cut_lowest_quartile_low = 12362; //from Plot_2d.C, With 20000 Ph/MeV LY
    Double_t PE_cut_lowest_quartile_high = 19100; //from Plot_2d.C, With 20000 Ph/MeV LY

    TH1D *hFpyreneNANR = CalculatefPyrene(treeDef,windowLowNR,windowHighNR,PE_cut_lowest_quartile);
    //TH1D *hFpyreneNANRLow = CalculatefPyrene(treeLow,windowLowNR,windowHighNR,PE_cut_lowest_quartile_low);
    //TH1D *hFpyreneNANRHigh = CalculatefPyrene(treeHigh,windowLowNR,windowHighNR,PE_cut_lowest_quartile_high);
    TH1D *hFpyreneNANRLow = new TH1D("","",50,0,1);
    TH1D *hFpyreneNANRHigh = new TH1D("","",50,0,1))
    TH1D *hFpyreneNR = CalculatefPyrene(treeNR,windowLowNR,windowHighNR,1e5);//No PE Cut on NR events,(Full acceptance)

    TH1D *hFPPS = (TH1D*)hFpyreneNANR->Clone();
    TH1D *hFPPSLow = (TH1D*)hFpyreneNANRLow->Clone();
    TH1D *hFPPSHigh = (TH1D*)hFpyreneNANRHigh->Clone();


    hFpyreneNANR->Scale(1.0/hFpyreneNANR->Integral());
    hFpyreneNANRLow->Scale(1.0/hFpyreneNANRLow->Integral());
    hFpyreneNANRHigh->Scale(1.0/hFpyreneNANRHigh->Integral());
    hFpyreneNR->Scale(1.0/hFpyreneNR->Integral());
    hFpyreneNR->Draw();
    hFpyreneNR->GetYaxis()->SetRangeUser(0,0.4);
    hFpyreneNANR->Draw("same");
    hFpyreneNANRLow->Draw("same");
    //hFpyreneNANRHigh->Draw("same");

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
    hFpyreneNR->SetTitle(" ; Fpyrene; Intensity [AU]");
    hFpyreneNANR->SetLineColor(kRed);
    hFpyreneNANRLow->SetLineColor(kBlack);
    hFpyreneNANRHigh->SetLineColor(kGreen);
    hFpyreneNR->SetLineStyle(2);
    hFpyreneNR->SetLineColor(kBlue);
    hFpyreneNANRLow->SetLineStyle(5);
    hFpyreneNR->GetYaxis()->SetTitleOffset(1.1);
    hFpyreneNR->GetYaxis()->SetTitleSize(32);
    hFpyreneNR->GetXaxis()->SetTitleSize(32);
    hFpyreneNR->GetYaxis()->SetLabelSize(28);
    hFpyreneNR->GetXaxis()->SetLabelSize(28);
    TLegend *legFp = new TLegend(0.55,0.65,0.95,0.9);
    legFp->AddEntry(hFpyreneNANR,"PPS Inlet Alphas: PLQY = 0.46 ","lp");
    legFp->AddEntry(hFpyreneNANRLow,"PPS Inlet Alphas: PLQY = 0.36 ","lp");
    //legFp->AddEntry(hFpyreneNANRHigh,"PPS Inlet Alphas - High PLQY","lp");
    legFp->AddEntry(hFpyreneNR,"Bulk LAr NRs","lp");
    //legFp->AddEntry(hOv,"Overlap","lp");
    legFp->SetFillStyle(3001);
    legFp->Draw("same");
    c2->Print("../plots/Fpyrene_draw_plot_paper.pdf");

    TFile *fout = new TFile("fpyrene_hist.root","RECREATE");
    fout->cd();
    hFpyreneNR->Write();
    hFpyreneNANR->Write();
    hFpyreneNANRLow->Write();
    hFpyreneNANRHigh->Write();
    fout->Close();
    //Calculate leakage

    //Assume 90% acceptance for LAr NRs
    TH1D *hLArCumu = (TH1D*) hFpyreneNR->GetCumulative();
    Double_t thresh = hLArCumu->GetBinCenter(hLArCumu->FindFirstBinAbove(0.1));
    cout << "90 % acceptance thresh = " <<thresh<<endl;
    //hLArCumu->Draw("same");

    //Double_t nr_start = hFpyreneNR->GetBinLowEdge(hFpyreneNR->FindFirstBinAbove(0));
    Double_t nr_start = thresh;
    printf("NR Distribution Start = %f \n",nr_start);

    Int_t nEvts = treeDef->GetEntries();
    //Integral of PPS distribution from start of NR distribution
    Double_t pps_leakage_n = hFPPS->Integral(hFPPS->FindBin(nr_start),hFPPS->GetNbinsX());
    printf("Leakage of PPS into NR region = %f \n",pps_leakage_n);

    //Find upper limit
    if (!gROOT->GetClass("TFeldmanCousins")) gSystem->Load("libPhysics");
    TFeldmanCousins f;
    Double_t upper_lim = f.CalculateUpperLimit(pps_leakage_n,0); //90% Upper limit on number of events with 0 background
    printf("Upper limit PPS into NR region = %f \n",upper_lim);
    Double_t leak = upper_lim/float(nEvts);//Leakage fraction
    printf("90% UL on LF for PPS into NR region = %f \n",leak);

    Int_t nEvtsLow = treeLow->GetEntries();
    //Integral of PPS distribution from start of NR distribution
    Double_t pps_leakage_n_low = hFPPSLow->Integral(hFPPSLow->FindBin(nr_start),hFPPSLow->GetNbinsX());
    printf("Leakage of low PLQY PPS into NR region = %f \n",pps_leakage_n_low);

    Double_t upper_lim_low = f.CalculateUpperLimit(pps_leakage_n_low,0); //90% Upper limit on number of events with 0 background
    printf("Upper limit low PLQY PPS into NR region = %f \n",upper_lim_low);
    Double_t leak_low= upper_lim_low/float(nEvtsLow);//Leakage fraction
    printf("90% UL on LF for low PLQY PPS into NR region = %f \n",leak_low);


    Int_t nEvtsHigh = treeHigh->GetEntries();
    //Integral of PPS distribution from start of NR distribution
    Double_t pps_leakage_n_high = hFPPSHigh->Integral(hFPPSHigh->FindBin(nr_start),hFPPSHigh->GetNbinsX());
    printf("Leakage of low PLQY PPS into NR region = %f \n",pps_leakage_n_high);

    Double_t upper_lim_high = f.CalculateUpperLimit(pps_leakage_n_high,0); //90% Upper limit on number of events with 0 background
    printf("Upper limit high PLQY PPS into NR region = %f \n",upper_lim_high);
    Double_t leak_high= upper_lim_high/float(nEvtsHigh);//Leakage fraction
    printf("90% UL on LF for high PLQY PPS into NR region = %f \n",leak_high);


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

  for(int iTrial = 0; iTrial<t->GetEntries();iTrial++){
    t->GetEntry(iTrial);
    Double_t fPyrene = 0.0;
    Double_t integral = 0.0;
    if(!times)continue;
    const int num = (*times).size();
    if(num==0)continue;
    if(numHits > PE_cut) continue; //introduce a PE Cut
    for(int i=0;i<numHits-1;i++){
      if((*times)[i] > winLow && (*times)[i] < winHigh) fPyrene++;
      if((*times)[i] > winLow && (*times)[i] < windowEnd) integral++;
    }
    fPyrene = fPyrene/integral;
    hFpyrene->Fill(fPyrene);
  }
  return hFpyrene;
}

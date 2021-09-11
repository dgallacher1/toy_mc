//Nice plots
//Eg for use:
//root -l draw_plot_paper.C'("updated_ap_1_type_0.root")'


#include "../src/ToyMC.h"
#include "../src/InverseCDF.h"
#include <math.h>       /* isnan, sqrt */

gSystem->Load("libtoymc");

void draw_pulseshape(string filename="default_PLQY_ap_1_type_0.root")
{
    //DEAP Style settings
    fStyle->DEAPStyle();

    //File in
    string path = "../output/";
    path+=filename;
    cout << "Reading.. " << path <<endl;
    TFile *fileinNADef = new TFile("../output/pulse_shape_ap_1_type_0.root");

    string pathNR = path.substr(0,path.find(".root")-1);
    pathNR+="1.root";
    TFile *fileinNR = new TFile(pathNR.c_str());

    //Read in tree
    TTree *treeDef = (TTree*) fileinNADef->Get("T");

    TTree *treeNR = (TTree*) fileinNR->Get("T");

    //Pulseshape parameters
    int nbins = 250;
    int high = 6000.;

    //Match colors from diagram in paper
    // int o_index = 9999;
    // TColor *c_orange = new TColor(o_index,0.902,0.404,0.173);
    // int g_index = 9998;
    // TColor *c_green = new TColor(g_index,33./255,122./255,11./255);


    TCanvas *c1 = new TCanvas("c1","Plots",1200,800);
    treeDef->Draw(Form("times>>hp(%i,-100,%i)",nbins,high));
    TH1D *hFraction = (TH1D*)gPad->GetPrimitive("hp");
    hFraction->SetTitle("Hit time distribution;Time[ns];Intensity (A.U)");
    hFraction->Scale(1.0/hFraction->Integral());

    //Read in LAr NR's
    treeNR->Draw(Form("times>>hp1(%i,-100,%i)",nbins,high));
    TH1D *hFractionNR = (TH1D*)gPad->GetPrimitive("hp1");
    hFractionNR->SetTitle(" ;Time[ns];Intensity (A.U)");
    hFractionNR->Scale(1.0/hFractionNR->Integral());
    hFractionNR->SetLineColor(fStyle->Color(0));

    /////////
    ///Formatting
    /////////

    hFraction->SetLineColor(kRed);
    hFractionNR->SetLineColor(kBlue);
    hFractionNR->SetLineStyle(2);

    gPad->SetLogy();
    //gPad->SetLogx();
    hFractionNR->Draw();
//    hFractionNR->SetTitle("Toy MC Pulseshape Comparison;Photon Arrival Time [ns];Intensity [AU]");
    hFractionNR->GetYaxis()->SetTitle("Intensity [AU]");
    hFractionNR->GetXaxis()->SetTitle("Photon Arrival Time [ns]");
    hFractionNR->GetYaxis()->SetTitleSize(36);
    hFractionNR->GetYaxis()->SetRangeUser(1e-5,1);
    hFractionNR->GetXaxis()->SetTitleSize(36);
    hFractionNR->GetYaxis()->SetTitleOffset(1.1);
    hFraction->Draw("same");
    TLegend *leg = new TLegend(0.4,0.6,0.9,0.9);
    leg->AddEntry(hFraction,"PPS Inlet Alphas pulse shape","lp");
    leg->AddEntry(hFractionNR,"Bulk LAr NR-like pulse shape","lp");
    leg->SetFillStyle(3001);
    leg->Draw("same");
    c1->Print("../plots/Pulseshape_draw_plot_paper_update.pdf");
    TFile *fout = new TFile("pulseshape_plot.root","RECREATE");
    fout->cd();
    hFractionNR->Write();
    hFraction->Write();
    fout->Close();

}

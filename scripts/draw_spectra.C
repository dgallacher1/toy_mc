void draw_spectra(){

  TFile *filein=new TFile("../dat/spectras.root","READ");

  TH1D *hAcry = (TH1D*)filein->Get("hAcrylicAttn");
  TH1D *hTPB = (TH1D*)filein->Get("hTPB");
  TH1D *hPMT = (TH1D*)filein->Get("hPMTEff");

  fStyle->DEAPStyle();

  hAcry->SetLineColor(fStyle->Color(0));
  hTPB->SetLineColor(fStyle->Color(1));
  hPMT->SetLineColor(fStyle->Color(2));

  hAcry->SetLineWidth(3);
  hTPB->SetLineWidth(3);
  hPMT->SetLineWidth(3);


  hTPB->Rebin(2);
  hTPB->Smooth();
  hTPB->Scale(1.0/hTPB->GetMaximum());

  TCanvas *c1 = new TCanvas("c1","Spectras",1200,800);
  gPad->SetGrid();


  hAcry->SetTitle("Wavelength Dependence for critical components;Wavelength [nm];Intensity/Efficiency (AU)");

  hAcry->Draw();
  hPMT->Draw("same");
  hTPB->Draw("same");

  TLegend *leg = new TLegend(0.6,0.5,0.95,0.8);
  leg->AddEntry(hAcry,"Acrylic Attenuation Probability (1 m path length)","lpf");
  leg->AddEntry(hPMT,"Hamamatsu R5912 PMT efficiency","lpf");
  leg->AddEntry(hTPB,"TPB Wavelength Emission Spectrum","lpf");
  leg->Draw("same");
  leg->SetFillStyle(3001);

  c1->Print("../plots/spectras.pdf");


  TH1D *hPyr = (TH1D*)filein->Get("hPyrene");
  hPyr->Draw();
  hPyr->Scale(1.0/hPyr->Integral(hPyr->FindBin(360.0),hPyr->FindBin(600.0)));
  ofstream fileout;
  fileout.open("../dat/pyrene.txt");
  //fileout << "QFlight_option: "energy",";
  fileout << "WLSCOMPONENT_option: wavelength,"<<endl;
  fileout << "WLSCOMPONENT_value1:[";

  for(int i=hPyr->FindBin(360.0);i < hPyr->FindBin(600.0);i++){
    fileout << Form(" %3.2f,",hPyr->GetBinCenter(i));
  }
  fileout << "],\n";
  fileout << "WLSCOMPONENT_value2:[";

  for(int i=hPyr->FindBin(360.0);i < hPyr->FindBin(600.0);i++){
    fileout << Form(" %2.5f,",hPyr->GetBinContent(i));
  }
  fileout << "],\n";
  fileout.close();




}

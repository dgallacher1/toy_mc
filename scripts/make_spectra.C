//Quick ROOT Macro to write wavelength dependent spectra and efficiency curves
//This requires RAT to run
//Author: David Gallacher

void make_spectra(){

  TFile *fileout = new TFile("../dat/spectras.root","RECREATE");

  RAT::DB *rdb = RAT::DB::Get();

  RAT::DBLinkPtr lpyrene = rdb->GetLink("OPTICS","pyrene_PS");
  vector< double > wlpy = lpyrene->GetDArray("WLSCOMPONENT_value1");
  vector< double > amppy = lpyrene->GetDArray("WLSCOMPONENT_value2");

  TH1D *hPyrene = new TH1D("hPyrene","hPyrene",wlpy.size()-1,&wlpy[0]);
  hPyrene->SetContent(&amppy[0]);

  RAT::DBLinkPtr llar = rdb->GetLink("OPTICS","liquid_Ar");
  vector< double > wllar = llar->GetDArray("SCINTILLATION_value1");
  vector< double > amplar = llar->GetDArray("SCINTILLATION_value2");

  TH1D *hLAr = new TH1D("hLAr","hLAr",wllar.size()-1,&wllar[0]);
  hLAr->SetContent(&amplar[0]);

  RAT::DBLinkPtr ltpb = rdb->GetLink("OPTICS","tpb");
  vector< double > wltpb = ltpb->GetDArray("WLSCOMPONENT_value1");
  vector< double > amptpb = ltpb->GetDArray("WLSCOMPONENT_value2");

  TH1D *hTPB = new TH1D("hTPB","hTPB",wltpb.size()-1,&wltpb[0]);
  hTPB->SetContent(&amptpb[0]);

  RAT::DBLinkPtr lpmt = rdb->GetLink("OPTICS","photocathode_R5912_HQE");
  vector< double > wlpmt = lpmt->GetDArray("EFFICIENCY_value1");
  vector< double > amppmt = lpmt->GetDArray("EFFICIENCY_value2");

  TH1D *hPMT = new TH1D("hPMT","hPMT",wlpmt.size()-1,&wlpmt[0]);
  hPMT->SetContent(&amppmt[0]);

  //Here we read the absorption lengths from the OPTICS table
  RAT::DBLinkPtr lLG = rdb->GetLink("OPTICS","acrylic_rpt");
  vector< double > WLabs = lLG->GetDArray("ABSLENGTH_value1");
  vector< double > Labs = lLG->GetDArray("ABSLENGTH_value2");

  TH1D *hAcrylicAttn = new TH1D("hAcrylicAttn","hAcrylicAttn",WLabs.size()-1,&WLabs[0]);
  hAcrylicAttn->SetContent(&Labs[0]);
  TGraph *gAcrylicAttn = new TGraph(hAcrylicAttn);

  Double_t distance = 100.0; // 100 cm travel distance

  //Apply Beer-Lambert law to acrylic absorption assuming 'distance' travelled by light
  for(int i=1;i<WLabs.size();i++){
    double prob = exp(-distance/gAcrylicAttn->Eval(hAcrylicAttn->GetBinCenter(i)));
    hAcrylicAttn->SetBinContent(i,prob);
  }

  fileout->cd();
  hPyrene->Write("hPyrene");
  hLAr->Write("hLAr");
  hTPB->Write("hTPB");
  hPMT->Write("hPMTEff");
  hAcrylicAttn->Write("hAcrylicAttn");

}

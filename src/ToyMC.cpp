#include "ToyMC.h"

ClassImp(ToyMC)

ToyMC::ToyMC(Int_t _seed)
{
  Init();
  seed = _seed;
}

ToyMC::ToyMC())
{
  Init();
}

ToyMC::~ToyMC()
{
  Clear();
}

//Run a single experiment
void ToyMC::RunToy(TTree *result_tree)
{
  //Ensure we're using TRandom3
  delete gRandom;
  gRandom = new TRandom3;
  //Set the seed
  gRandom->SetSeed(seed);

  //Prep work
  InverseCDF *invcdf= new InverseCDF();
  LoadFunctions();

  //Set TPB Time binning for inverse cdf
  vector<Double_t> bins;
  bins.push_back(-2e-12); bins.push_back(2e-12); bins.push_back(4e-12); bins.push_back(0.1);
  for (int k=4; bins[k-1]+0.1*pow(k-3.,1.7)<windowEnd; k++) bins.push_back(bins[k-1] + 0.1*pow(k-3.,1.5));
  bins.push_back(windowEnd);
  TH1D *HistTPB= new TH1D("h_tpb_wls","Hist",bins.size()-1,&bins[0]);
  HistTPB->Add(TPBPulseShape);
  gTPBTime = invcdf->GetInverseHisto(HistTPB);


  //Make result tree
  vector<Double_t> times;
  Double_t edep = 0.0;
  Double_t shadowFraction = 0.0;
  Int_t numPhotons = 0;
  Int_t numPyrenePhotons = 0;
  Int_t numAP = 0;
  Int_t trialNum = 0;
  Int_t experimentNum = 0;

  result_tree = new TTree("T","Result Tree");
  result_tree->Branch("times",&times);
  result_tree->Branch("edep",&edep,"edep/D");
  result_tree->Branch("numPhotons",&numPhotons,"numPhotons/I");
  result_tree->Branch("numPyrenePhotons",&numPyrenePhotons,"numPyrenePhotons/I");
  result_tree->Branch("numAP",&numAP,"numAP/I");
  result_tree->Branch("shadowFraction",&shadowFraction,"shadowFraction/D");
  result_tree->Branch("trialNum",&trialNum,"trialNum/I");
  result_tree->Branch("experimentNum",&numAP,"experimentNum/I");

  TStopwatch sw;
  //Loop over all experiments and trials and run a trial and fill tree.
  for(Int_t iExperiment=0;iExperiment<numExperiments;iExperiment++){
    for(Int_t iTrial=0;iTrial<numTrials;iTrial++){
      //Pass parameters
      trialNum = iTrial;
      experimentNum = iExperiment;
      //Run the trial
      DoTrial(&edep,&shadowFraction,&numPhotons,&numAP,&times);
      //Fill ttree with results
      tree->Fill();
    }//End loop over trials
  }//End loop over experiments
  sw.Stop();
  cout <<"Toy completed: "<< "Ran "<< numExperiments << " for " <<numTrials << "trials each. "<<endl;
  sw.Print();
}

//Here's where we do the toy simulation
void ToyMC::DoTrial(Double_t &edep, Double_t &shadowFraction, Int_t &numPhotons,Int_t &numPyrenePhotons, Double_t &numAP,vector<Double_t> &times)
{
  //Sample from our Landau distribution
  edep = fRandEnergy->GetRandom();
  //Can shadow half to all light
  shadowFraction = gRandom->Uniform(0.5,1.0);
  //Mean photon production based on this edep
  Double_t meanNumPhotons = edep*lightYield;
  numPhotons = int(gRandom->Gaus(meanNumPhotons,sqrt(meanNumPhotons)));
  //Photons that arent "shadowed" are LAr photons
  Int_t numLArPhotons = numPhotons*(1.0-shadowFraction);
  //Shadowed photons are pyrene WLSd photons
  numPyrenePhotons= numPhotons-numLArPhotons;
  //Loop over LAr photons
  for(Int_t iLAr = 0;iLAr<numLArPhotons;iLAr++)
  {
    Double_t t = 0.0; //T0 for each photon is 0
    t+=fLArPS->GetRandom();//Get LAr scint time offset
    Double_t wl = gLArSpec->Eval(gRandom->Rndm());//Evaluate the inverse cdf for a random wavelength to start with

    //Add delay for TPB timing and swap to TPB wavelength, assumes 100% WLSE for TPB
    t+=gTPBTime->Eval(gRandom->Rndm());
    wl=gTPBSpec->Eval(gRandom->Rndm());

    //Acrylic Attenuation is wavelength dependent (gAcrylicAttn is survival prob) so if we attenuate the photon, skip
    if(gAcrylicAttn->Eval(wl) < gRandom->Rndm())continue;
    //PMT Efficiency
    if(gPMTEff->Eval(wl) < gRandom->Rndm())continue;
    times.push_back(t);//hit pmt sucessfully

    //After-pulsing
    if(APProb < gRandom->Rndm())continue;
    t+=fPMT->GetRandom();
    times.push_back(t);
    numAP++;
  }

  //Loop over pyrene photons
  for(Int_t iPy = 0;iPy<numPyrenePhotons;iPy++)
  {
    Double_t t = 0.0; //T0 for each photon is 0
    t+=fLArPS->GetRandom();//Get LAr scint time offset
    //Evaluate the inverse cdf for a random wavelength to start with (Can be modified to WL dependent probabilities in future)
    Double_t wl = gLArSpec->Eval(gRandom->Rndm());

    //WLS the photon
    if(pyrene_WLSE < gRandom->Rndm()) continue;
    t+=fPyrenePS->GetRandom(); //Add delay time
    wl=gPyreneSpec->Eval(gRandom->Rndm());//Sample a random wavelength

    //Sample from the number of reflections as a binomial with probability of 0.5
    Int_t numReflections = gRandom->Binomial(2*meanReflections,0.5);
    Double_t reflProb = pow(reflectionProb,numReflections);
    //If we dont survive the successive reflections then kill this photon
    if(reflProb < gRandom->Rndm())continue;

    //Acrylic Attenuation is wavelength dependent (gAcrylicAttn is survival prob) so if we attenuate the photon, skip
    if(gAcrylicAttn->Eval(wl) < gRandom->Rndm())continue;
    //PMT Efficiency
    if(gPMTEff->Eval(wl) < gRandom->Rndm())continue;
    times.push_back(t); //hit pmt sucessfully

    //After-pulsing
    if(APProb < gRandom->Rndm())continue;
    t+=fPMT->GetRandom();
    times.push_back(t);
    numAP++;
  }

}

void ToyMC::LoadPDFs(TFile *SpectrumFile)
{

  InverseCDF *invcdf= new InverseCDF();
  //Load PDF Spectrums
  TH1D *hPyreneSpec = (TH1D*) SpectrumFile->Get("hPyrene");
  TH1D *hTPBSpec = (TH1D*) SpectrumFile->Get("hTPB");
  TH1D *hLArSpec = (TH1D*) SpectrumFile->Get("hLAr");
  TH1D *hAcrylicAttn = (TH1D*) SpectrumFile->Get("hAcrylic");
  TH1D *hPMTEff= (TH1D*) SpectrumFile->Get("hPMTEff");

  //Get Inverse cdf versions of spectrums
  gPyreneSpec = invcdf->GetInverseHisto(hPyreneSpec);
  gTPBSpec = invcdf->GetInverseHisto(hTPBSpec);
  gLArSpec = invcdf->GetInverseHisto(hLArSpec);

  //PDFs of acrylic survival probability and pmt efficiency
  gAcrylicAttn = new TGraph(hAcrylicAttn);
  gPMTEff = new TGraph(hPMTEff);

}

void ToyMC::LoadFunctions()
{
  //Default TPB 128nm Photon Values
  Double_t paramtpb[]= {4.6,12000.,2.e7,0.06,2518.,0.86,8.3,0.078,130.};
  fTPBPS = new TF1("fTPBPS",this,&ToyMC::TPBPulseShape,-2e-12,windowEnd,9,"ToyMC","TPBPulseshape");
  SetTPBParameters(paramtpb);

  //Default NR LAr parameters
  Double_t paramlar[]= {0.87,2.2,0.06,75.5,0.112,1445.0};
  fLArPS = new TF1("fLArPS",this,&ToyMC::LArPulseShape,0.001,windowEnd,6,"ToyMC","LArPulseShape");
  SetLArParameters(paramlar);

  //Default Pyrene_PS Pulse-shape, 2 exponential
  Double_t parampy[]= {0.8,250,0.2,105};
  fPyrenePS = new TF1("fPyrenePS","[0]*exp(-x/[1])+[2]*exp(-x/[3])",0.001,windowEnd);
  SetPyreneParameters(parampy);

  //Random Energy deposition following a landau distribution with tail to low energies
  fRandEnergy = new TF1("fland", "TMath::Landau(-x,[0],[1])",0,6);
  fRandEnergy->SetParameter(0,-meanEnergy);
  fRandEnergy->SetParameter(1,0.2); //A guess-FIXME, gives a reasonable tail

  //PMT AP Pulse-shape, from liquid Ar pulseshape paper (D3600)
  Double_t parampmt[]= {0.67,1660,680,0.33,6300,1350};
  fPMT = new TF1("fPMT","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])",0.001,windowEnd);
  SetPMTParameters(parampmt);

}

void ToyMC::SetPMTParameters(const Double_t *params)
{
  if(fPMT)fPMT->SetParameters(params);
}

void ToyMC::SetPyreneParameters(const Double_t *params)
{
  if(fPyrenePS)fPyrenePS->SetParameters(params);
}

void ToyMC::SetTPBParameters(const Double_t *params)
{
  if(fTPBPS)fTPBPS->SetParameters(params);
}

void ToyMC::SetLArParameters(const Double_t *params)
{
  if(fLArPS)fLArPS->SetParameters(params);
}

//Copied from DEAP-3600 uses model from https://arxiv.org/abs/2001.09855
Double_t ToyMC::LArPulseShape(Double_t *t,Double_t *par)
{
  Double_t tt = t[0];
  Double_t ampLAr = par[0]/par[1]*exp(-tt/par[1]);
  ampLAr+= (par[2])/((1.0+tt/par[3])*(1.0+tt/par[3]))*1.0/par[3];
  ampLAr+= par[4]/par[5]*exp(-tt/par[5]);
  return ampLAr;
}

Double_t ToyMC::TPBPulseShape(Double_t *t,Double_t *par)
{
  //par[n] = [A, t_a, tau_T, Nd, tau_d, Np, tau_S, Nm, tau_m]//Parameter names for reference from Stanford and Westerdale (1804.06895)
  Double_t tt = t[0];
  Double_t y = ROOT::Math::expint(-par[1]/par[2]);
  Double_t f = ROOT::Math::expint(-(tt+par[1])/par[2]) -y;
  Double_t f1 = exp(-2.*tt/par[2]);
  Double_t f2 = pow(1. + par[0]*f, 2.);
  Double_t f3 = (1. + tt/par[1]);
  Double_t fd= par[3]/par[4] * f1/(f2*f3);
  Double_t fp= par[5]/par[6] * exp(-tt/par[6]);
  Double_t fm= par[7]/par[8] * exp(-tt/par[8]);
  return fp + fm + fd;
}


/// Reset class
void ToyMC::Clear(Option_t *option)
{
  numTrials=0;
  numExperiments=0;
  reflectionProb = 0.0;
  windowEnd = 0.0; //ns
  pyrene_WLSE = 0.0; //Relative to TPB
  lightYield = 0.0; //Photons/MeV
  meanEnergy = 0.0;
  APProb = 0.0;

  delete gPyreneSpec;
  delete gLArSpec;
  delete gTPBSpec;
  delete gPMTEff;
  delete gAcrylicAttn;
  delete gTPBTime;

  delete fPyrenePS;
  delete fLArPS;
  delete fPMT;
  delete fTPBPS;
  delete fRandEnergy;
}

///Initialize parameters
void ToyMC::Init()
{
  seed = 1234;//Default seed
  numTrials=10;
  numExperiments=1;
  reflectionProb = 0.9;
  windowEnd = 10000.0; //ns
  pyrene_WLSE = 0.40; //Relative to TPB
  lightYield = 40000.0; //Photons/MeV
  meanReflections = 1; //Mean of a binomial for number of reflections
  meanEnergy = 5.4;//MeV
  APProb = 0.075;//

}

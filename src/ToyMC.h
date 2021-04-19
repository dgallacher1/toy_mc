#ifndef ROOT_ToyMC
#define ROOT_ToyMC

///////////////////////////////////
/// Toy MC program for Slow WLS coating in neck
/// Used to find optimum pulse-shape parameters
/// Author: David Gallacher
/// Date: April 3rd
///
///////////////////////////////////

#include "TObject.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TMath.h"
#include "TStopwatch.h"
#include <Math/SpecFuncMathMore.h>
#include "TROOT.h"

#include "Util.h"
#include "InverseCDF.h"

#include <vector>
#include <iostream>

using namespace std;

class ToyMC : public TObject{
private:
  // Private data members
  Int_t       seed;
  Int_t       numTrials;
  Int_t       numExperiments;
  Double_t    windowEnd;          // End of timing window in ns
  Double_t    reflectionProb;     // Probability to survive reflection
  Double_t    pyrene_WLSE;        // Probability of WLSing pyrene
  Double_t    lightYield;         // LAr light yield
  Int_t       meanReflections;    // Mean number of reflections Pyrene photons undergo
  Double_t    meanEnergy;         // Mean energy for the starting alpha
  Double_t    APProb;             // Combined Afterpulse probability
  Double_t    quenching;          // Quenching factor

  //Distributions loaded
  TGraph*       gPyreneSpec;      // Pyrene emission spectrum
  TGraph*       gTPBSpec;         // TPB emission spectrum
  TGraph*       gLArSpec;         // LAr emission specturm
  TGraph*       gAcrylicAttn;     // Acrylic survival probability vs wavelength assuming 1m of travel distance
  TGraph*       gPMTEff;          // PMT Efficiency vs wavelength
  TGraph*       gTPBTime;         // TPB Time in inverseCDF config, for speeding up sims, faster than TF1->GetRandom()
  TGraph*       gAr39Energy;      // Ar39 Beta energy spectrum from RAT

  //Timing functions
  TF1*        fPMT;               // PMT AP timing distribution, considering 2 main AP functions
  TF1*        fPyrenePS;          // Pyrene pulseshape, 2 exp model
  TF1*        fTPBPS;             // TPB Pulseshape, from rat-deap
  TF1*        fLArPS;             // LAr Pulseshape model from https://arxiv.org/abs/2001.09855
  TF1*        fGeo;               // Geometric broadening of DEAP-3600 ER events from https://arxiv.org/abs/2001.09855

  //Edep
  TF1*        fRandEnergy;        // Landau distribution model for energy deposition


public:
   ToyMC(Int_t _seed);
   ToyMC();
   ~ToyMC();
   void       Init();
   void       Clear(Option_t *option = "");


   ///Run the toyMC
   TTree*     RunToy(Int_t type);
   //Each trial is ran here
   void       DoNeckAlphaTrial(Double_t &edep, Double_t &shadowFraction, Int_t &numPhotons,Int_t &numPyrenePhotons, Int_t &numAP,vector<Double_t> &times,Int_t &numHits);
   void       DoLArTrial(Int_t type, Double_t &edep, Double_t &shadowFraction, Int_t &numPhotons,Int_t &numPyrenePhotons, Int_t &numAP,vector<Double_t> &times,Int_t &numHits);


   ///Set and Get the number of trials for this toy
   Int_t      GetNumTrials(){return numTrials;}
   void       SetNumTrials(Int_t _numTrials){numTrials = _numTrials;}
   ///Set and get the number of experiments for this toy
   Int_t      GetNumExperiments(){return numExperiments;}
   void       SetNumExperiments(Int_t _numExperiments){numExperiments = _numExperiments;}

   ///Load the appropriate PDFs for running
   void       LoadPDFs(TFile *SpectrumFile);

   ///Initialize timing functions
   void       LoadFunctions();
   ///Set parameters for PMT Timing Function
   void       SetPMTParameters(const Double_t *params);
   ///Set parameters for Pyrene Timing Function
   void       SetPyreneParameters(const Double_t *params);
   ///Set parameters for TPB Timing Function
   void       SetTPBParameters(const Double_t *params);
   ///Set parameters for LAr Timing function
   void       SetLArParameters(const Double_t *params);

   //Set Other parameters

   //Set and Get mean energy for edep
   Double_t   GetMeanEnergy(){return meanEnergy;}
   void       SetMeanEnergy(Double_t _meanEnergy){meanEnergy = _meanEnergy;}

   //Set and Get the quenching factors
   Double_t   GetQuenchingFactor(){return quenching;}
   void       SetQuenchingFactor(Double_t _quenching){quenching = _quenching;}

   //Set and get light yield for LAr
   Double_t   GetLightYield(){return lightYield;}
   void       SetLightYield(Double_t _lightYield){lightYield = _lightYield;}

   //Set and get Mean number of reflections for pyrene photons
   Int_t      GetMeanReflections(){return meanReflections;}
   void       SetMeanReflections(Int_t _meanReflections){meanReflections = _meanReflections;}

   //Set and get the reflection probability (assuming uniform in WL)
   Double_t   GetReflectionProb(){return reflectionProb;}
   void       SetReflectionProb(Double_t _reflectionProb){reflectionProb = _reflectionProb;}

   //Set and get the pyrene WLS efficiency (Assuming uniform in WL)
   Double_t   GetPyreneWLSE(){return pyrene_WLSE;}
   void       SetPyreneWLSE(Double_t _pyrene_WLSE){pyrene_WLSE = _pyrene_WLSE;}

   //Set and get the size of the window
   Double_t   GetWindow(){return windowEnd;}
   void       SetWindow(Double_t _windowEnd){ windowEnd = _windowEnd;}

   //TF1 Functions for pulseshapes
   Double_t   LArPulseShape(Double_t *t,Double_t *par);
   Double_t   TPBPulseShape(Double_t *t,Double_t *par);

   //Get Pulse-shape functions
   TF1*       GetPMTPS(){return fPMT;}
   TF1*       GetTPBPS(){return fTPBPS;}
   TF1*       GetLArPS(){return fLArPS;}
   TF1*       GetPyrenePS(){return fPyrenePS;}
   TF1*       GetGeoPS(){return fGeo;}



   ClassDef(ToyMC,1)
};

#endif

#ifndef ROOT_InverseCDF
#define ROOT_InverseCDF

///////////////////////////////////
/// Helper class to read p.d.fs into inverse CDF's for random sampling
/// Author: David Gallacher
/// Date: April 2021
///
///////////////////////////////////

#include "TObject.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TGraph.h"

#include "Util.h"

#include <vector>
#include <iostream>

using namespace std;

class InverseCDF : public TObject{
private:
  // Private data members


public:
   InverseCDF();
   ~InverseCDF();
   void       Init();
   void       Clear(Option_t *option = "");

   ///Generate the Inverse CDF from a histogram containing the probability density function
   TGraph*    GetInverseHisto(TH1D *hPDF);
   ///Generate the Inverse CDF from vectors containing the x and y data of a distribution
   TGraph*    GetInverse(vector<Double_t> x_vals, vector<Double_t> y_vals);
   ///Generate a histogram of a probability density function
   TH1D*      GetPDF(vector<Double_t> x_vals, vector<Double_t> y_vals);

   ClassDef(InverseCDF,1)
};

#endif

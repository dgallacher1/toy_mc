///////////////////////////////////
///Main program for ToyMC
///Takes in command line arguments for seed,
///
///
///////////////////////////////////

#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"

#include "Util.h"
#include "ToyMC.h"
#include "InverseCDF.h"

#include <iostream>   // std::cout
#include <string>     // std::string, std::stoi
#include <vector>
#include <sstream>

//Prototypes
vector<string> ParseCommandLineArguments(int argc, char** argv);
int main(int argc, char** argv);

using namespace std;
int main(int argc, char **argv)
{
  //Initial Parameters
  int seed = 1234;
  int type = 0;
  int numTrials = 10;
  int numExperiments = 1;
  int ap_mode = 1;//Default is afterpulsing enabled
  string foutname = "output";

  // Parse the command line arguments.
  vector< string > argsVec = ParseCommandLineArguments(argc, argv);
  if(argc<5){
    cerr << "Usage: " << argv[0] << " seed type(0==N.A, 1== NR, 2==ER) numTrials numExperiments fileoutName(optional) afterpulsing on/off(1 == on, 0 == off [optional])" << endl;
    return 1;
  }

  seed = atoi(argsVec[1].c_str());
  type = atoi(argsVec[2].c_str());
  numTrials = atoi(argsVec[3].c_str());
  numExperiments = atoi(argsVec[4].c_str());
  if(argc > 5) foutname = argsVec[5];
  if(argc==7) ap_mode = atoi(argsVec[6].c_str());

  TFile *DataFile = new TFile("../dat/spectras.root","READ");

  Double_t quench = 0.71;
  if(type==1) quench = 0.25;
  if(type==2) quench = 1.0;

  ToyMC *the_toy = new ToyMC(seed);
  the_toy->SetNumTrials(numTrials);
  the_toy->SetNumExperiments(numExperiments);
  the_toy->LoadPDFs(DataFile);
  the_toy->SetQuenchingFactor(quench);
  the_toy->SetLightYield(20000);//Artifically reduce light yield in order to speed up simulations
  the_toy->SetAfterPulsing(ap_mode);//Turn afterpulsing on/off (1 == on, 0 == off)
  //the_toy->SetPyreneWLSE(0.59);//TUM+Queens increase
  //the_toy->SetPyreneWLSE(0.46);//TUM
  //the_toy->SetPyreneWLSE(0.35);//TUM + reduction from Marcin

  //To change parameters, call it after LoadFunctions()

  string filename = "../output/"+foutname+Form("_ap_%i_type_%i.root",ap_mode,type);
  TFile *fileout = new TFile(filename.c_str(),"RECREATE");

  //Run the toy MC
  TTree *result_tree = the_toy->RunToy(type);

  //Save to file
  //fileout->cd();
  result_tree->Write("T",TObject::kOverwrite);
  fileout->Close();

  return 0;
}


//Function to parse the command line arguments
vector<string> ParseCommandLineArguments(int argc, char** argv)
{
  vector< string > vecStr;
  stringstream tmpStream;
  string tmpStr = "";
  for(Int_t iArg = 0; iArg < argc; iArg++){
    tmpStr = ""; tmpStream.clear(); tmpStream.str("");
	tmpStream << argv[iArg];
	tmpStream >> tmpStr;
	vecStr.push_back(tmpStr);
  }

  return vecStr;
}

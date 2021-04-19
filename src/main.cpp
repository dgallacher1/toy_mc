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
  string foutname = "output";

  // Parse the command line arguments.
  vector< string > argsVec = ParseCommandLineArguments(argc, argv);
  if(argc<5){
    cerr << "Usage: " << argv[0] << " seed type(0==N.A, 1== NR, 2==ER) numTrials numExperiments fileoutName(optional)" << endl;
    return 1;
  }

  seed = atoi(argsVec[1].c_str());
  type = atoi(argsVec[2].c_str());
  numTrials = atoi(argsVec[3].c_str());
  numExperiments = atoi(argsVec[4].c_str());
  if(argc==6) foutname = argsVec[5];

  TFile *DataFile = new TFile("../dat/spectras.root","READ");

  Double_t quench = 0.72;
  if(type==1) quench = 0.8;
  if(type==2) quench = 1.0;

  ToyMC *the_toy = new ToyMC(seed);
  the_toy->SetNumTrials(numTrials);
  the_toy->SetNumExperiments(numExperiments);
  the_toy->LoadPDFs(DataFile);
  the_toy->SetQuenchingFactor(quench);
  the_toy->SetLightYield(2000);//Artifically reduce light yield in order to speed up simulations
  //To change parameters, call it after LoadFunctions()

  string filename = "../output/"+foutname+"_"+Form("%i.root",type);
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

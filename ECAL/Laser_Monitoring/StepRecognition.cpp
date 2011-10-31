//Test program for the class StepRecognition.h
//Compile with: g++ `root-config --libs --glibs --cflags` StepRecognition.cpp -o TestStep
//Execute with: ./TestStep


#include "StepRecognition.h"

#include "TFile.h"
#include "TString.h"
#include "TProfile.h"


int main ()
{
  
  //Open the ROOT file
  TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_Mar.root");
  //TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_Apr.root");
  //TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_May.root");
  //TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_Jun.root");
  //TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_Jul.root");
  //TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_Aug.root");
  //TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_Sep.root");
  //TFile* inputFile = TFile::Open("LaserEEAnalysis_fed_605_Oct.root");
  
  
  //Create the output file
  TString *outputFile = new TString("outputFile.txt");
  
  //Open the TProfile
  TProfile *LaserAmplitude = (TProfile*) inputFile -> Get("p_matacqampli_las_fed605_harness6");
  
  //Call the identification function
  Step_ID(LaserAmplitude, *outputFile);
  
  
  return 1;
  
}

#include "Zutils.h"
#include "setTDRStyle.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TTree.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <strstream>


int main(int argc, char **argv){

   //set the style
  setTDRStyle();
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 
  gStyle->SetOptFit(1);

 /// Acquisition from cfg file
 
 if(argc != 2){
 std::cerr << " >>>>> analysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
 return 1;
 }

 parseConfigFile (argv[1]) ;
 
 std::string treeName  = gConfigParser -> readStringOption("Input::treeName");
 std::cout<<" Input Tree Name = "<<treeName<<std::endl;
 std::string inputDataFile =  gConfigParser -> readStringOption("Input::inputDataFile");
 std::cout<<" Input Data File = "<<inputDataFile<<std::endl;
 std::string inputMCFile =  gConfigParser -> readStringOption("Input::inputMCFile");
 std::cout<<" Input MC File = "<<inputMCFile<<std::endl;
 std::string WeightforMC =  gConfigParser -> readStringOption("Input::WeightforMC");
 std::cout<<" Weights for MC = "<<WeightforMC<<std::endl;
 std::string outputFile =  gConfigParser -> readStringOption("Output::outputFile");
 std::cout<<" Output Data File = "<<outputFile<<std::endl;

 /// Input data infos

 TChain* treeDATA = new TChain(treeName.c_str());
 TChain* treeMC = new TChain(treeName.c_str());

 FillChain(*treeDATA,inputDataFile.c_str());
 FillChain(*treeMC,inputMCFile.c_str());

 std::cout << " MC: " << std::setw(8) << treeMC->GetEntries() << " entries" << std::endl; 
 std::cout << " DATA: " << std::setw(8) << treeDATA->GetEntries() << " entries" << std::endl;
  
 if (treeDATA->GetEntries() == 0 || treeMC->GetEntries() == 0 ){
    std::cout << ">>>recalibZ::Error: at least one file is empty" << std::endl; 
    return -1;
  }

 std::vector<std::string> FitCategories;
 FitCategories = gConfigParser -> readStringListOption("Input::FitCategories");

 std::cout << " >>>>> Input::FitCategories size = " << FitCategories.size() << std::endl;  
 std::cout << " >>>>> >>>>>  "; 
 for (int iCat = 0; iCat < FitCategories.size(); iCat++){
  std::cout << " " << FitCategories.at(iCat) << ", ";
 }
 std::cout << std::endl; 
 
 //--- weights for MC
 TFile weightsFile (WeightforMC.c_str(),"READ"); 
 TH1F* hweights = (TH1F*)weightsFile.Get("hweights");
 float w[100];
 for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
 weightsFile.Close();

 /// Output infos
 TFile* outputTFile = new TFile(outputFile.c_str(),"RECREATE");

 /// option Infos
 int nbinZ  = gConfigParser -> readIntOption("Option::nbinZ");
 std::cout<<" nbinZ = "<<nbinZ<<std::endl;
 double mZ_Max =  gConfigParser -> readDoubleOption("Option::mZMax");
 std::cout<<" mZ_Max = "<<mZ_Max<<std::endl;
 double mZ_Min =  gConfigParser -> readDoubleOption("Option::mZMin");
 std::cout<<" mZ_Min = "<<mZ_Min<<std::endl;
 double scaleEB = gConfigParser -> readDoubleOption("Option::scaleEB");
 std::cout<<" scaleEB = "<<scaleEB<<std::endl;
 double scaleEE = gConfigParser -> readDoubleOption("Option::scaleEE");
 std::cout<<" scaleEE = "<<scaleEE<<std::endl;
 int nPoints  = gConfigParser -> readIntOption("Option::nPoints");
 std::cout<<" nPoints = "<<nPoints<<std::endl;


 ///**** Book histos

 std::map<std::string,TH1F*> ZmassDATA;
 std::map<std::string,TH1F*> ZmassDATA_regression;
 std::map<std::string,TH1F*> ZmassMC;
 std::map<std::string,TH1F*> ZmassMC_regression;


 for(unsigned int i = 0; i < FitCategories.size(); ++i){
   std::string category = FitCategories.at(i);
   std::string histoName1 = "h_ZmassDATA_"+category;
   ZmassDATA[category] = new TH1F(histoName1.c_str(),"",nbinZ,mZ_Min,mZ_Max);
   ZmassDATA[category] -> Sumw2();

   std::string histoName2 = "h_ZmassDATA_regression_"+category;
   ZmassDATA_regression[category] = new TH1F(histoName2.c_str(),"",nbinZ,mZ_Min,mZ_Max);
   ZmassDATA_regression[category] -> Sumw2();

   std::string histoName3 = "h_ZmassMC_"+category;
   ZmassMC[category] = new TH1F(histoName3.c_str(),"",nbinZ,mZ_Min,mZ_Max);
   ZmassMC[category] -> Sumw2();

   std::string histoName4 = "h_ZmassMC_regression_"+category;
   ZmassMC_regression[category] = new TH1F(histoName4.c_str(),"",nbinZ,mZ_Min,mZ_Max);
   ZmassMC_regression[category] -> Sumw2();

  }

 /// Set branch addresses
 int isZ;
 float ele1ele2_scM,ele1ele2_scM_regression;
 int ele1_isEB,ele2_isEB;
 float ele1_scEta,ele2_scEta,ele1_seedE,ele2_seedE,ele1_scE,ele2_scE,ele1_es,ele2_es,ele1_scERaw,ele2_scERaw, ele1_scE_regression,
       ele2_scE_regression;
 int ele1_seedIeta,ele1_seedIphi,ele2_seedIeta,ele2_seedIphi,ele1_seedIx,ele2_seedIx,ele1_seedIy,ele2_seedIy,ele1_seedZside,ele2_seedZside;
 int PUit_NumInteractions;
 
 treeDATA->SetBranchAddress("isZ", &isZ);
 treeDATA->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM);
 treeDATA->SetBranchAddress("ele1ele2_scM_regression", &ele1ele2_scM_regression);
 treeDATA->SetBranchAddress("ele1_isEB",  &ele1_isEB);
 treeDATA->SetBranchAddress("ele2_isEB",  &ele2_isEB);

 treeMC->SetBranchAddress("isZ", &isZ);
 treeMC->SetBranchAddress("ele1ele2_scM", &ele1ele2_scM);
 treeMC->SetBranchAddress("ele1_isEB",  &ele1_isEB);
 treeMC->SetBranchAddress("ele2_isEB",  &ele2_isEB);
 treeMC->SetBranchAddress("PUit_NumInteractions", &PUit_NumInteractions);

 treeDATA->SetBranchAddress("ele1_seedIeta", &ele1_seedIeta);
 treeDATA->SetBranchAddress("ele1_seedIphi", &ele1_seedIphi);
 treeDATA->SetBranchAddress("ele2_seedIeta", &ele2_seedIeta);
 treeDATA->SetBranchAddress("ele2_seedIphi",  &ele2_seedIphi);
 treeDATA->SetBranchAddress("ele1_seedIx",  &ele1_seedIx);
 treeDATA->SetBranchAddress("ele2_seedIx",  &ele2_seedIx);
 treeDATA->SetBranchAddress("ele1_seedIy",  &ele1_seedIy);
 treeDATA->SetBranchAddress("ele2_seedIy",  &ele2_seedIy);
 treeDATA->SetBranchAddress("ele1_seedZside",  &ele1_seedZside);
 treeDATA->SetBranchAddress("ele2_seedZside",  &ele2_seedZside);

 treeMC->SetBranchAddress("ele1_seedIeta", &ele1_seedIeta);
 treeMC->SetBranchAddress("ele1_seedIphi", &ele1_seedIphi);
 treeMC->SetBranchAddress("ele2_seedIeta", &ele2_seedIeta);
 treeMC->SetBranchAddress("ele2_seedIphi",  &ele2_seedIphi);
 treeMC->SetBranchAddress("ele1_seedIx",  &ele1_seedIx);
 treeMC->SetBranchAddress("ele2_seedIx",  &ele2_seedIx);
 treeMC->SetBranchAddress("ele1_seedIy",  &ele1_seedIy);
 treeMC->SetBranchAddress("ele2_seedIy",  &ele2_seedIy);
 treeMC->SetBranchAddress("ele1_seedZside",  &ele1_seedZside);
 treeMC->SetBranchAddress("ele2_seedZside",  &ele2_seedZside);


 treeDATA->SetBranchAddress("ele1_scE", &ele1_scE);
 treeDATA->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
 treeDATA->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
 treeDATA->SetBranchAddress("ele2_scE", &ele2_scE);
 treeDATA->SetBranchAddress("ele1_scERaw", &ele1_scERaw);
 treeDATA->SetBranchAddress("ele1_es", &ele1_es);
 treeDATA->SetBranchAddress("ele2_scERaw", &ele2_scERaw);
 treeDATA->SetBranchAddress("ele2_es", &ele2_es);

 treeMC->SetBranchAddress("ele1_scE", &ele1_scE);
 treeMC->SetBranchAddress("ele1_scE_regression", &ele1_scE_regression);
 treeMC->SetBranchAddress("ele2_scE_regression", &ele2_scE_regression);
 treeMC->SetBranchAddress("ele1_scERaw", &ele1_scERaw);
 treeMC->SetBranchAddress("ele1_es", &ele1_es);
 treeMC->SetBranchAddress("ele2_scE", &ele2_scE);
 treeMC->SetBranchAddress("ele2_scERaw", &ele2_scERaw);
 treeMC->SetBranchAddress("ele2_es", &ele2_es);
   
 //*** Loop on MC **//
 std::cout <<" Fill with MC Events "<<std::endl;
 int nEntriesMC = treeMC -> GetEntries();
 for(int iEntry = 0; iEntry < nEntriesMC; ++iEntry){

  if( (iEntry%100000 == 0) ) std::cout << " reading saved entry " << iEntry << std::endl;
  treeMC -> GetEntry(iEntry);
  double weight = w[PUit_NumInteractions];
   //only the Z
   if (isZ != 1) continue;
   if( (ele1_seedZside== 0) && (ele2_seedZside == 0) ){
	  ZmassMC["EB-EB"] ->  Fill( ele1ele2_scM * sqrt(scaleEB*scaleEB),weight );
	  ZmassMC_regression["EB-EB"] ->  Fill( ele1ele2_scM*sqrt((ele1_scE_regression/ele1_scE)*(ele2_scE_regression/ele2_scE)) * sqrt(scaleEB*scaleEB),weight);
	}    
   else if( fabs(ele1_seedZside== 1) && fabs(ele2_seedZside== 1) ){
	  ZmassMC["EE-EE"] ->  Fill( ele1ele2_scM * sqrt(scaleEE*scaleEE),weight );
	  ZmassMC_regression["EE-EE"] ->  Fill( ele1ele2_scM*sqrt((ele1_scE_regression/ele1_scE)*(ele2_scE_regression/ele2_scE)) * sqrt(scaleEE*scaleEE),weight );

	}    
   else{
	  ZmassMC["EB-EE"] ->  Fill( ele1ele2_scM * sqrt(scaleEB*scaleEE),weight );
	  ZmassMC_regression["EB-EE"] ->  Fill( ele1ele2_scM* sqrt(scaleEB*scaleEE)* sqrt((ele1_scE_regression/ele1_scE)*(ele2_scE_regression/ele2_scE)),weight );
       }

 }

 //*** Loop over Data **//
 int nEntriesDATA = treeDATA -> GetEntries();
 std::cout <<" Fill with DATA Events "<<std::endl;
 for(int iEntry = 0; iEntry < nEntriesDATA; ++iEntry){

  if( (iEntry%100000 == 0) ) std::cout << " reading saved entry " << iEntry << std::endl;
  treeDATA -> GetEntry(iEntry);
   //only the Z
   if (isZ != 1) continue;
      
   if( (ele1_seedZside== 0) && (ele2_seedZside== 0) ){
	  ZmassDATA["EB-EB"] ->  Fill( ele1ele2_scM * sqrt(scaleEB*scaleEB) );
	  ZmassDATA_regression["EB-EB"] ->  Fill( ele1ele2_scM_regression * sqrt(scaleEB*scaleEB));
	}    
   else if( fabs(ele1_seedZside== 1) && fabs(ele2_seedZside== 1) ){
	  ZmassDATA["EE-EE"] ->  Fill( ele1ele2_scM * sqrt(scaleEE*scaleEE) );
	  ZmassDATA_regression["EE-EE"] ->  Fill( ele1ele2_scM_regression * sqrt(scaleEE*scaleEE) );

	}    
   else{
	  ZmassDATA["EB-EE"] ->  Fill( ele1ele2_scM * sqrt(scaleEB*scaleEE) );
	  ZmassDATA_regression["EB-EE"] ->  Fill( ele1ele2_scM_regression * sqrt(scaleEB*scaleEE) );
       }

 }

 /// Z Lineshape Tool

  BinnedFitZPeak("EB-EB", 1, ZmassDATA["EB-EB"], ZmassMC["EB-EB"], nPoints, mZ_Min, mZ_Max);
//   BinnedFitZPeak("EE-EE", 2, ZmassDATA["EE-EE"], ZmassMC["EE-EE"], nPoints, mZ_Min, mZ_Max);
//   BinnedFitZPeak("EB-EE", 1, ZmassDATA["EB-EE"], ZmassMC["EB-EE"], nPoints, mZ_Min, mZ_Max);
//   BinnedFitZPeak("EB-EB", 1, ZmassDATA_regression["EB-EB"], ZmassMC_regression["EB-EB"], nPoints, mZ_Min, mZ_Max);
//   BinnedFitZPeak("EE-EE", 2, ZmassDATA_regression["EE-EE"], ZmassMC_regression["EE-EE"], nPoints, mZ_Min, mZ_Max);
//   BinnedFitZPeak("EB-EE", 1, ZmassDATA_regression["EB-EE"], ZmassMC_regression["EB-EE"], nPoints, mZ_Min, mZ_Max);

 /// Output save 
 std::cout << "recalibZ::Saving and Closing" << std::endl;

 ZmassDATA["EB-EB"]->Write();
 ZmassDATA["EB-EE"]->Write();
 ZmassDATA["EE-EE"]->Write();
 ZmassMC["EB-EB"]->Write();
 ZmassMC["EB-EE"]->Write();
 ZmassMC["EE-EE"]->Write();

 ZmassDATA_regression["EB-EB"]->Write();
 ZmassDATA_regression["EB-EE"]->Write();
 ZmassDATA_regression["EE-EE"]->Write();
 ZmassMC_regression["EB-EB"]->Write();
 ZmassMC_regression["EB-EE"]->Write();
 ZmassMC_regression["EE-EE"]->Write();

 outputTFile -> Close();

 return 0;
}

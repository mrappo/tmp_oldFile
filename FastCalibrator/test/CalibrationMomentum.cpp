// per compilare: g++ -Wall -o CalibrationMomentum `root-config --cflags --glibs` CalibrationMomentum.cpp

#include "../CommonTools/TEndcapRings.h"
#include "../CommonTools/histoFunc.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "treeReader.h"




#define xtalWidth 0.01745329
#define PI        3.1415926536 



using namespace std;

 
bool IsEtaGap(float eta){
  float feta = fabs(eta);
  if( fabs(feta - 0 )<3) return true;
  if( fabs(feta - 25)<3) return true;
  if( fabs(feta - 45)<3) return true;
  if( fabs(feta - 65)<3) return true;
  if( fabs(feta - 85)<3) return true;
  return false;
}

int templIndexEB(float eta){
    float feta = fabs(eta);
    if (feta <= 25)               {return 0;}
    if (feta>  25 && feta <=  45) {return 0;}
    if (feta>  45 && feta <=  65) {return 0;}
    if (feta>  65 && feta <=  85) {return 0;}
    return -1;
}

int templIndexEE(float eta){
    float feta = fabs(eta);
    if (feta>  85 && feta <=  98) {return 0;}
    if (feta>  98 && feta <= 100) {return 0;}
    if (feta> 100 && feta <= 118) {return 0;}
    if (feta> 118 )               {return 0;}
    
    return -1;
}
   


//**************  MAIN PROGRAM **************************************************************
int main(int argc, char** argv){

  /// Acquisition from cfg file
 if(argc != 2){
 std::cerr << ">>>>> analysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
 return 1;
 }

 parseConfigFile (argv[1]) ;
 
 std::string TreeName = gConfigParser -> readStringOption("Input::TreeName");
 std::string infileDATA = gConfigParser -> readStringOption("Input::infileDATA");
 std::string infileMC = gConfigParser -> readStringOption("Input::infileMC");
 std::string WeightforMC =  gConfigParser -> readStringOption("Input::WeightforMC");

 int  nPhiBinsEB = gConfigParser -> readIntOption("Input::nPhiBinsEB");
 int  nPhiBinsEE = gConfigParser -> readIntOption("Input::nPhiBinsEE");
 int  nEtaBinsEB = gConfigParser -> readIntOption("Input::nEtaBinsEB");
 int  nEtaBinsEE = gConfigParser -> readIntOption("Input::nEtaBinsEE");
 int  nPhiBinsTempEB = gConfigParser -> readIntOption("Input::nPhiBinsTempEB");
 int  nPhiBinsTempEE = gConfigParser -> readIntOption("Input::nPhiBinsTempEE");
 int  rebinEB = gConfigParser -> readIntOption("Input::rebinEB");
 int  rebinEE = gConfigParser -> readIntOption("Input::rebinEE");


 std::string outputFile = gConfigParser -> readStringOption("Output::outputFile");

 cout <<" Basic Configuration " <<endl;
 cout <<" Tree Name = "<<TreeName<<endl;
 cout <<" infileDATA = "<<infileDATA<<endl;
 cout <<" infileMC = "<<infileMC<<endl;
 cout <<" WeightforMC = "<<WeightforMC<<endl;
 cout <<" nPhiBinsEB = "<<nPhiBinsEB<<endl;
 cout <<" nPhiBinsEE = "<<nPhiBinsEE<<endl;
 cout <<" nEtaBinsEB = "<<nEtaBinsEB<<endl;
 cout <<" nEtaBinsEE = "<<nEtaBinsEE<<endl;
 cout <<" nPhiBinsTempEB = "<<nPhiBinsTempEB<<endl;
 cout <<" nPhiBinsTempEE = "<<nPhiBinsTempEE<<endl;
 cout <<" rebinEB = "<<rebinEB<<endl;
 cout <<" rebinEE = "<<rebinEE<<endl;



 cout << "Making calibration plots for Momentum scale studies "<< endl;

     
  //---- variables for selection
  float r9min = 0.00 ;
  float r9max = 9999 ;  
  float etaMax  = 2.5;
  float eta2Max = 2.5;

  bool usePUweights = true;

  //--- weights for MC
  TFile weightsFile (WeightforMC.c_str(),"READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");

  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();


  //----- NTUPLES--------------------
  TChain *ntu_DA = new TChain(TreeName.c_str());
  TChain *ntu_MC = new TChain(TreeName.c_str());

  if(!FillChain(*ntu_DA, infileDATA.c_str())) return 1;

  if(!FillChain(*ntu_MC, infileMC.c_str())) return 1;

  std::cout << "     DATA: " << ntu_DA->GetEntries() << " entries in Data sample" << std::endl;
  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;

  // observables  
  int isW;
  float PV_z;
  float EoP, scEta, scPhi, mZ;
  float scEta2,scEne2,scPhi2;
  float scE3x3, scE5x5, scEne, scERaw,scEt;  
  float charge, scLocalEta, scLocalPhi,crackCorr,localCorr; 
  float pTK,pTK2; 
  float scEneReg,scEneReg2;
  int ele1_iphi,ele1_ieta,ele1_ix,ele1_iy,ele1_iz; 
  int ele2_iphi,ele2_ieta,ele2_ix,ele2_iy,ele2_iz;
  int npu;
  
  // Set branch addresses for Data  
  ntu_DA->SetBranchAddress("isW", &isW);
  ntu_DA->SetBranchAddress("PV_z", &PV_z);
  ntu_DA->SetBranchAddress("ele1ele2_scM",     &mZ);
  ntu_DA->SetBranchAddress("ele1_scEta", &scEta);
  ntu_DA->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_DA->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_DA->SetBranchAddress("ele2_scPhi", &scPhi2);
  ntu_DA->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_DA->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_DA->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_DA->SetBranchAddress("ele1_scERaw", &scERaw);
  ntu_DA->SetBranchAddress("ele1_scE", &scEne);
  ntu_DA->SetBranchAddress("ele2_scE", &scEne2);
  ntu_DA->SetBranchAddress("ele1_scE_regression", &scEneReg);
  ntu_DA->SetBranchAddress("ele2_scE_regression", &scEneReg2);
  ntu_DA->SetBranchAddress("ele1_charge", &charge);
  ntu_DA->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_DA->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_DA->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_DA->SetBranchAddress("ele1_scLocalContCorr",&localCorr); 
  ntu_DA->SetBranchAddress("ele1_scEt",&scEt); 
  ntu_DA->SetBranchAddress("ele1_tkP",&pTK); 
  ntu_DA->SetBranchAddress("ele2_tkP",&pTK2); 
  ntu_DA->SetBranchAddress("ele1_seedIphi",&ele1_iphi); 
  ntu_DA->SetBranchAddress("ele1_seedIeta",&ele1_ieta); 
  ntu_DA->SetBranchAddress("ele1_seedIx",&ele1_ix); 
  ntu_DA->SetBranchAddress("ele1_seedIy",&ele1_iy); 
  ntu_DA->SetBranchAddress("ele1_seedZside",&ele1_iz); 
  ntu_DA->SetBranchAddress("ele2_seedIphi",&ele2_iphi); 
  ntu_DA->SetBranchAddress("ele2_seedIeta",&ele2_ieta); 
  ntu_DA->SetBranchAddress("ele2_seedIx",&ele2_ix); 
  ntu_DA->SetBranchAddress("ele2_seedIy",&ele2_iy); 
  ntu_DA->SetBranchAddress("ele2_seedZside",&ele2_iz); 

  // Set branch addresses for MC  
  ntu_MC->SetBranchAddress("PUit_NumInteractions", &npu);
  ntu_MC->SetBranchAddress("isW", &isW);
  ntu_MC->SetBranchAddress("PV_z", &PV_z);
  ntu_MC->SetBranchAddress("ele1ele2_scM",     &mZ);
  ntu_MC->SetBranchAddress("ele1_scEta", &scEta);
  ntu_MC->SetBranchAddress("ele2_scEta", &scEta2);
  ntu_MC->SetBranchAddress("ele1_scPhi", &scPhi);
  ntu_MC->SetBranchAddress("ele1_EOverP", &EoP);
  ntu_MC->SetBranchAddress("ele1_e3x3", &scE3x3);
  ntu_MC->SetBranchAddress("ele1_e5x5", &scE5x5);
  ntu_MC->SetBranchAddress("ele1_scE", &scEne);
  ntu_MC->SetBranchAddress("ele2_scE", &scEne2);
  ntu_MC->SetBranchAddress("ele1_scERaw", &scERaw);
  ntu_MC->SetBranchAddress("ele1_charge", &charge);
  ntu_MC->SetBranchAddress("ele1_scLocalPhi",&scLocalPhi); 
  ntu_MC->SetBranchAddress("ele1_scLocalEta",&scLocalEta); 
  ntu_MC->SetBranchAddress("ele1_scCrackCorr",&crackCorr); 
  ntu_MC->SetBranchAddress("ele1_scEt",&scEt); 
  ntu_MC->SetBranchAddress("ele1_tkP",&pTK); 
  ntu_MC->SetBranchAddress("ele2_tkP",&pTK2); 
  ntu_MC->SetBranchAddress("ele1_seedIphi",&ele1_iphi); 
  ntu_MC->SetBranchAddress("ele1_seedIeta",&ele1_ieta); 
  ntu_MC->SetBranchAddress("ele1_seedIx",&ele1_ix); 
  ntu_MC->SetBranchAddress("ele1_seedIy",&ele1_iy); 
  ntu_MC->SetBranchAddress("ele1_seedZside",&ele1_iz); 
    

  // histogram definition in EB and fit functions

  TH1F*** h_EoP_EB = new TH1F**[nPhiBinsEB];   
  TH1F*** h_EoC_EB = new TH1F**[nPhiBinsEB]; 
  TH1F*** h_Phi_EB = new TH1F**[nPhiBinsEB]; // used to map iEta (as defined for Barrel and Endcap geom) into eta 
  TF1*** f_EoP_EB = new TF1**[nPhiBinsEB];
  TF1*** f_EoC_EB = new TF1**[nPhiBinsEB];

  std::vector< std::pair<int,int> > refIdEB;
  std::pair<int,int> temp ; temp.first=0; temp.second=0;
  refIdEB.assign(nPhiBinsEB,temp) ;

  // histogram definition in EE and fit functions

  TH1F*** h_EoP_EE = new TH1F**[nPhiBinsEE];   
  TH1F*** h_EoC_EE = new TH1F**[nPhiBinsEE];  
  TH1F*** h_Phi_EE = new TH1F**[nPhiBinsEE]; // used to map iEta (as defined for Barrel and Endcap geom) into eta 
  TF1*** f_EoP_EE = new TF1**[nPhiBinsEE];
  TF1*** f_EoC_EE = new TF1**[nPhiBinsEE];

  std::vector< std::pair<int,int> > refIdEE;
  refIdEE.assign(nPhiBinsEE,temp) ;

  // Initializate histos in EB
  for(Int_t i = 0; i < nPhiBinsEB; ++i)
  {
    h_EoP_EB[i] = new TH1F* [nEtaBinsEB];
    h_EoC_EB[i] = new TH1F* [nEtaBinsEB];
    h_Phi_EB[i] = new TH1F* [nEtaBinsEB];

    for(Int_t j=0; j< nEtaBinsEB ; j++)
    {
     TString histoName;
     histoName= Form("EoP_%d_%d_EB", i,j);
     h_EoP_EB[i][j] = new TH1F (histoName, histoName, 2000, 0., 2.);
     h_EoP_EB[i][j]->SetFillColor(kRed+2);
     h_EoP_EB[i][j]->SetLineColor(kRed+2);
     h_EoP_EB[i][j]->SetFillStyle(3004);

    histoName=Form("EoC_%d_%d_EB", i,j);
    h_EoC_EB[i][j] = new TH1F(histoName, histoName, 2000, 0., 2.);
    h_EoC_EB[i][j]->SetFillColor(kGreen+2);
    h_EoC_EB[i][j]->SetLineColor(kGreen+2);
    h_EoC_EB[i][j]->SetFillStyle(3004);
   
    histoName=Form("Phi_%d_%d_EB", i,j);   
    h_Phi_EB[i][j] = new TH1F(histoName, histoName, 360, 0., 360.); 
   }
  }

  // Initializate histos in EE
  for(Int_t i = 0; i < nPhiBinsEE; ++i)
  {
    h_EoP_EE[i] = new TH1F* [nEtaBinsEE];
    h_EoC_EE[i] = new TH1F* [nEtaBinsEE];
    h_Phi_EE[i] = new TH1F* [nEtaBinsEE];

    for(Int_t j=0; j< nEtaBinsEE ; j++)
    {
 
    TString histoName;
    histoName=Form("EoP_%d_%d_EE", i,j);
    h_EoP_EE[i][j] = new TH1F(histoName, histoName, 2000, 0., 2.);
    h_EoP_EE[i][j]->SetFillColor(kRed+2);
    h_EoP_EE[i][j]->SetLineColor(kRed+2);
    h_EoP_EE[i][j]->SetFillStyle(3004);

    histoName=Form("EoC_%d_%d_EE", i,j);
    h_EoC_EE[i][j] = new TH1F(histoName, histoName, 2000, 0., 2.);
    h_EoC_EE[i][j]->SetFillColor(kGreen+2);
    h_EoC_EE[i][j]->SetLineColor(kGreen+2);
    h_EoC_EE[i][j]->SetFillStyle(3004);
   
    histoName=Form("Phi_%d_%d_EE", i,j);   
    h_Phi_EE[i][j] = new TH1F(histoName, histoName, 360, 0., 360.); 
  }
 }

  // Template in EE and EB

  TH1F*** h_template_EB = new TH1F**[nPhiBinsTempEB];
  TH1F*** h_template_EE = new TH1F**[nPhiBinsTempEE];
  
 
  for (Int_t imod = 0; imod<nPhiBinsTempEB; imod++){
    h_template_EB[imod] = new TH1F*[nEtaBinsEB];
    for(Int_t j = 0; j<nEtaBinsEB; j++){
    TString histoName;
    histoName=Form("template_%d_%d_EB", imod,j);
    h_template_EB[imod][j] = new TH1F(histoName, "", 2000, 0., 2.);
   }
  }

  for (Int_t imod = 0; imod<nPhiBinsTempEE; imod++){
   h_template_EE[imod] = new TH1F*[nEtaBinsEE];
   for(Int_t j = 0; j<nEtaBinsEE; j++){
    TString histoName;
    histoName=Form("template_%d_%d_EE", imod,j);
    h_template_EE[imod][j] = new TH1F(histoName, "", 2000, 0., 2.);  
  }
 }
  
  TH1F** h_phi_data_EB = new TH1F*[nEtaBinsEB];
  TH1F** h_phi_mc_EB   = new TH1F*[nEtaBinsEB];
  TH1F** h_phi_data_EE = new TH1F*[nEtaBinsEE];
  TH1F** h_phi_mc_EE   = new TH1F*[nEtaBinsEE];

  for (Int_t index = 0; index < nEtaBinsEB ; index++){
    TString name;
    name=Form("h_phi_data_EB_%d",index);
    h_phi_data_EB[index]= new TH1F(name,"h_phi_data",100,-TMath::Pi(),TMath::Pi());
    name=Form("h_phi_mc_EB_%d",index);
    h_phi_mc_EB[index]= new TH1F(name,"h_phi_mc",100,-TMath::Pi(),TMath::Pi());
   }

  for(Int_t index = 0; index < nEtaBinsEE ; index++){
    TString name;
    name=Form("h_phi_data_EE_%d",index);
    h_phi_data_EE[index] = new TH1F(name,"h_phi_data",100,-TMath::Pi(),TMath::Pi());
    name=Form("h_phi_mc_EE_%d",index);
    h_phi_mc_EE[index] =  new TH1F(name,"h_phi_mc",100,-TMath::Pi(),TMath::Pi());
  }

  
  TH1F* h_et_data  = new TH1F("h_et_data","h_et_data",100,0,100);
  TH1F* h_et_mc    = new TH1F("h_et_mc","h_et_mc",100,0,100);


  // Initialize endcap geometry
  TEndcapRings *eRings = new TEndcapRings(); 

  //**************************** loop on MC, make refernce and fit dist
  float ww = 1;
  float var = 0;
  std::cout << "Loop in MC events " << endl; 
  for(int entry = 0; entry < ntu_MC->GetEntries(); ++entry) {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;
    //    if (entry>1000) break;

    ntu_MC->GetEntry(entry);
    
    if (isW==1) continue;
    if( fabs(scEta) > etaMax ) continue;
    if( fabs(scEta2) > eta2Max ) continue;
   
    float R9 = scE3x3/scEne;
    if ( R9 < r9min || R9 > r9max ) continue; 

    //--- PU weights
    if (usePUweights) ww = w[npu];
 
    //--- set ieta for the Endcaps
    int iphi,ieta;

    var = (mZ * sqrt(pTK/scEne))/91.19;    /// use the momentum for ele1

    if(ele1_iz==0){
       int modPhi = int (ele1_iphi/(360./nPhiBinsTempEB));
       int modEta  = templIndexEB(ele1_ieta);
 
       if(modPhi == nPhiBinsTempEB) h_template_EB[0][modEta] ->  Fill(var,ww);
       else h_template_EB[modPhi][modEta]->Fill(var,ww);

       // fill MC histos in eta bins
       int PhibinEB = int (ele1_iphi/(360./nPhiBinsEB)); 
 
       if(PhibinEB==nPhiBinsEB){ h_EoP_EB[0][modEta]-> Fill(var,ww);  // This is MC 
                                 h_phi_mc_EB[modEta]->Fill(scPhi,ww);}
       else{ h_EoP_EB[PhibinEB][modEta] -> Fill(var,ww);  // This is MC
             std::pair<int,int> temp;
             temp.first=modPhi; temp.second=modEta;
             refIdEB.at(PhibinEB)=temp; 
             h_phi_mc_EB[modEta]->Fill(scPhi,ww);}     
    }

    else{ 
          iphi = eRings->GetEndcapIphi(ele1_ix,ele1_iy,ele1_iz);
          ieta = eRings->GetEndcapIeta(ele1_ix,ele1_iy,ele1_iz);
          int modPhi = int (iphi/(360./nPhiBinsTempEE));
          int modEta =  templIndexEE(ieta);
          if(modPhi==nPhiBinsTempEE) h_template_EE[0][modEta] ->  Fill(var,ww);
          else h_template_EE[modPhi][modEta] ->  Fill(var,ww);

          // fill MC histos in eta bins
          int PhibinEE = int (iphi/(360./nPhiBinsEE)); 
          if(PhibinEE==nPhiBinsEE){ h_EoP_EE[0][modEta] -> Fill(var,ww);  // This is MC
                                    h_phi_mc_EE[modEta]->Fill(scPhi,ww);}
          else{
               h_EoP_EE[PhibinEE][modEta] -> Fill(var,ww);  // This is MC
               std::pair<int,int> temp;
               temp.first=modPhi; temp.second=modEta;
               refIdEE.at(PhibinEE)=temp; 
               h_phi_mc_EE[modEta]->Fill(scPhi,ww);}

        }

    h_et_mc ->Fill(scEt,ww);
  } 

  // loop on events in Data
  std::cout << "Loop in Data events " << endl; 

  for(int entry = 0; entry < ntu_DA->GetEntries(); ++entry) {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;

    ntu_DA->GetEntry(entry);

    if (isW==1) continue;
    if ( fabs(scEta) > etaMax) continue;
    if ( fabs(scEta2) > eta2Max) continue;

    float R9 = scE3x3/scEne;
    if ( R9 < r9min || R9 > r9max ) continue; 

    //--- set ieta for the Endcaps
    int iphi,ieta;
   
    var  = (mZ * sqrt(pTK/scEne))/91.19 ;    /// use the momentum for ele1

    if(ele1_iz==0){

       int PhibinEB = int (ele1_iphi/(360./nPhiBinsEB)); 
       int modEta =  templIndexEB(ele1_ieta);

       if(PhibinEB==nPhiBinsEB){ h_EoC_EB[0][modEta]-> Fill(var);  // This is DATA
                                 h_phi_data_EB[modEta]->Fill(scPhi);
                          }
       else{h_EoC_EB[PhibinEB][modEta] -> Fill(var);  // This is DATA
            h_Phi_EB[PhibinEB][modEta]->Fill((double) ele1_iphi);
                         
            h_phi_data_EB[modEta]->Fill(scPhi);}
      
     }
    else{  iphi = eRings->GetEndcapIphi(ele1_ix,ele1_iy,ele1_iz); 
           ieta = eRings->GetEndcapIeta(ele1_ix,ele1_iy,ele1_iz);
           int PhibinEE = int (iphi/(360./nPhiBinsEE));
           int modEta =  templIndexEE(ieta);

           if(PhibinEE==nPhiBinsEE){
                              h_EoC_EE[0][modEta] -> Fill(var);  // This is DATA
                              h_phi_data_EE[modEta] ->Fill(scPhi);}
           else{
                  h_EoC_EE[PhibinEE][modEta] -> Fill(var);  // This is DATA
                  h_Phi_EE[PhibinEE][modEta] -> Fill((double) iphi); 
                  h_phi_data_EE[modEta] ->Fill(scPhi);}
           
        }
     
    //use also the other electron
    var = (mZ * sqrt(pTK2/scEne2))/91.19 ;    /// use the momentum for ele2

    if(ele2_iz==0){

       int PhibinEB = int (ele2_iphi/(360./nPhiBinsEB)); 
       int modEta =  templIndexEB(ele2_ieta);
  
       if(PhibinEB==nPhiBinsEB){ h_EoC_EB[0][modEta]-> Fill(var);  // This is DATA
                                 h_phi_data_EB[modEta]->Fill(scPhi);
                          }
       else{h_EoC_EB[PhibinEB][modEta] -> Fill(var);  // This is DATA
            h_Phi_EB[PhibinEB][modEta]->Fill((double) ele2_iphi);
            h_phi_data_EB[modEta]->Fill(scPhi);}
      
     }
    else{     
           iphi = eRings->GetEndcapIphi(ele2_ix,ele2_iy,ele2_iz); 
           ieta = eRings->GetEndcapIeta(ele2_ix,ele2_iy,ele2_iz); 

           int PhibinEE = int (iphi/(360./nPhiBinsEE)); 
           int modEta =  templIndexEE(ieta);
           if(PhibinEE==nPhiBinsEE){
                              h_EoC_EE[0][modEta] -> Fill(var);  // This is DATA
                              h_phi_data_EE[modEta] ->Fill(scPhi);}
           else{
                  h_EoC_EE[PhibinEE][modEta] -> Fill(var);  // This is DATA
                  h_Phi_EE[PhibinEE][modEta] -> Fill((double) iphi); 
                  h_phi_data_EE[modEta] ->Fill(scPhi);}
           
        }

   
    h_et_data ->Fill(scEt);

  }

  std::cout << "End loop: Analyze events " << endl; 
    
  // draw results   
  TGraphErrors** g_EoP_EB   = new TGraphErrors*[nEtaBinsEB];
  TGraphErrors** g_EoC_EB   = new TGraphErrors*[nEtaBinsEB];
  TGraphErrors** g_Rat_EB   = new TGraphErrors*[nEtaBinsEB];
  TGraphErrors** g_EoP_EE   = new TGraphErrors*[nEtaBinsEB];
  TGraphErrors** g_EoC_EE   = new TGraphErrors*[nEtaBinsEB];
  TGraphErrors** g_Rat_EE   = new TGraphErrors*[nEtaBinsEB];
 
  /// Template binned Functions  
  histoFunc ***templateHistoFuncEB= new histoFunc** [nPhiBinsTempEB]; 
  histoFunc ***templateHistoFuncEE= new histoFunc** [nPhiBinsTempEE]; 

  for (int mod=0; mod<nPhiBinsTempEB;mod++) {
    templateHistoFuncEB[mod]= new histoFunc* [nEtaBinsEB];
    for (int j=0; j<nEtaBinsEB; j++){
     h_template_EB[mod][j] -> Rebin(rebinEB);
     templateHistoFuncEB[mod][j] = new histoFunc(h_template_EB[mod][j]);
   }
 }

 for (int mod=0; mod<nPhiBinsTempEE;mod++) {
    templateHistoFuncEE[mod]= new histoFunc* [nEtaBinsEE];
    for (int j=0; j<nEtaBinsEE; j++){
    h_template_EE[mod][j] -> Rebin(rebinEE);
    templateHistoFuncEE[mod][j] = new histoFunc(h_template_EE[mod][j]);
  }
}
  // Template Fit in EB

 for(int j=0; j<nEtaBinsEB; j++){
   
   g_EoP_EB[j]= new TGraphErrors();
   g_EoC_EB[j]= new TGraphErrors();
   g_Rat_EB[j]= new TGraphErrors();
 }

 
 for(int i = 0; i < nPhiBinsEB; ++i){
   
   f_EoP_EB[i] = new TF1*[nEtaBinsEB];
   f_EoC_EB[i] = new TF1*[nEtaBinsEB];

    for(int j=0; j<nEtaBinsEB; j++){
 
    h_EoP_EB[i][j] -> Rebin(rebinEB);    
    h_EoC_EB[i][j] -> Rebin(rebinEB);    
   
    std::pair<int,int> mod = refIdEB.at(i); 
    
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
    
    char funcName[50];
    sprintf(funcName,"f_EoP_%d_%d_Ref_%d_%d_EB",i,j,mod.first,mod.second);
    f_EoP_EB[i][j] = new TF1(funcName, templateHistoFuncEB[mod.first][mod.second], 0.6, 1.3, 3, "histoFunc");
    f_EoP_EB[i][j] -> SetParName(0,"Norm"); 
    f_EoP_EB[i][j] -> SetParName(1,"Scale factor"); 
    f_EoP_EB[i][j] -> SetLineWidth(1); 
    f_EoP_EB[i][j] -> SetLineColor(kRed+2); 
    f_EoP_EB[i][j] -> SetNpx(10000);
    h_EoP_EB[i][j] -> Sumw2();
   
    // uncorrected    
    double xNorm = h_EoP_EB[i][j]->Integral()/h_template_EB[mod.first][mod.second]->Integral() *
                   h_EoP_EB[i][j]->GetBinWidth(1)/h_template_EB[mod.first][mod.second]->GetBinWidth(1); 

    f_EoP_EB[i][j] -> FixParameter(0, xNorm);
    f_EoP_EB[i][j] -> SetParameter(1, 1.05 );
    f_EoP_EB[i][j] -> FixParameter(2, 0.);
    
    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<10;trial++) {
      rp = h_EoP_EB[i][j] -> Fit(funcName, "QRLW+");
      fStatus = rp;
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
     }
    
    float flPhi = h_Phi_EB[i][j]->GetMean(); 

    if(i==0) g_EoP_EB[j] -> SetPoint(i, 0. , pow(f_EoP_EB[i][j]->GetParameter(1),2));
    else  g_EoP_EB[j] -> SetPoint(i, int(flPhi) , pow(f_EoP_EB[i][j]->GetParameter(1),2));

    g_EoP_EB[j] -> SetPointError(i, 0., 2*f_EoP_EB[i][j]->GetParError(1));
        
    //cout  << " ***** " <<  1./f_EoP[i]->GetParameter(1) << " " << f_EoP[i]->GetParError(1) << " " << fStatus <<  endl; 
  
    //ratio preparation
    float rat = f_EoP_EB[i][j]->GetParameter(1);
    float era = f_EoP_EB[i][j]->GetParError(1); 
    
    xNorm = h_EoC_EB[i][j]->Integral()/h_template_EB[mod.first][mod.second]->Integral() *
            h_EoC_EB[i][j]->GetBinWidth(1)/h_template_EB[mod.first][mod.second]->GetBinWidth(1); 

    sprintf(funcName,"f_EoC_%d_%d_Ref_%d_%d_EB",i,j,mod.first,mod.second);
    
    f_EoC_EB[i][j] = new TF1(funcName, templateHistoFuncEB[mod.first][mod.second], 0.6, 1.3, 3, "histoFunc");
    f_EoC_EB[i][j] -> SetParName(0,"Norm"); 
    f_EoC_EB[i][j] -> SetParName(1,"Scale factor"); 
    f_EoC_EB[i][j] -> SetLineWidth(1); 
    f_EoC_EB[i][j] -> SetLineColor(kGreen+2); 
    f_EoC_EB[i][j] -> SetNpx(10000);

    f_EoC_EB[i][j] -> FixParameter(0, xNorm);
    f_EoC_EB[i][j] -> SetParameter(1, 1.05 );
    f_EoC_EB[i][j] -> FixParameter(2, 0.);
    
    std::cout << "***** Re-Fitting ";
    for (int trial=0;trial<10;trial++) {
      rp = h_EoC_EB[i][j] -> Fit(funcName, "QRL+");
      fStatus = rp;
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
    }
    
    if(i==0) g_EoC_EB[j] -> SetPoint(i, 0., pow(f_EoC_EB[i][j]->GetParameter(1),2));
    else g_EoC_EB[j] -> SetPoint(i, int(flPhi), pow(f_EoC_EB[i][j]->GetParameter(1),2));
 
    g_EoC_EB[j] -> SetPointError(i, 0., 2*f_EoC_EB[i][j]->GetParError(1));
    //    cout << " ********** " <<  1./f_EoC[i]->GetParameter(1) << " " << f_EoC[i]->GetParError(1) << endl; 

    //ratio finalization
    rat /= f_EoC_EB[i][j]->GetParameter(1);
    era = rat*sqrt(era*era+f_EoC_EB[i][j]->GetParError(1)*f_EoC_EB[i][j]->GetParError(1)); 
    
    if(i==0) g_Rat_EB[j] -> SetPoint(i, 0. , rat);
    else  g_Rat_EB[j] -> SetPoint(i, int(flPhi) , rat);

    g_Rat_EB[j] -> SetPointError(i,  0. , era); 
    g_Rat_EB[j]->SetLineColor(kBlue+2); 
  }
 }

 // Template Fit in EE
 for(int j=0; j<nEtaBinsEE; j++){
   
   g_EoP_EE[j]= new TGraphErrors();
   g_EoC_EE[j]= new TGraphErrors();
   g_Rat_EE[j]= new TGraphErrors();
 }

 for(int i = 0; i < nPhiBinsEE; ++i){
   
   f_EoP_EE[i]=new TF1*[nEtaBinsEE];
   f_EoC_EE[i]=new TF1*[nEtaBinsEE];

    for(int j=0; j<nEtaBinsEE; j++){
    h_EoP_EE[i][j] -> Rebin(rebinEE);    
    h_EoC_EE[i][j] -> Rebin(rebinEE);    
    
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
    std::pair<int,int> mod = refIdEE.at(i); 
    
    char funcName[50];
    sprintf(funcName,"f_EoP_%d_%d_Ref_%d_%d_EE",i,j,mod.first,mod.second);
    f_EoP_EE[i][j] = new TF1(funcName, templateHistoFuncEE[mod.first][mod.second], 0.2, 1.3, 3, "histoFunc");
    f_EoP_EE[i][j] -> SetParName(0,"Norm"); 
    f_EoP_EE[i][j] -> SetParName(1,"Scale factor"); 
    f_EoP_EE[i][j] -> SetLineWidth(1); 
    f_EoP_EE[i][j] -> SetLineColor(kRed+2); 
    f_EoP_EE[i][j] -> SetNpx(10000);
    f_EoP_EE[i][j] -> SetNpx(10000);
    h_EoP_EE[i][j] -> Sumw2();
    
    // uncorrected    
    double xNorm = h_EoP_EE[i][j]->Integral()/h_template_EE[mod.first][mod.second]->Integral() *
                   h_EoP_EE[i][j]->GetBinWidth(1)/h_template_EE[mod.first][mod.second]->GetBinWidth(1); 
    f_EoP_EE[i][j] -> FixParameter(0, xNorm);
    f_EoP_EE[i][j] -> SetParameter(1, 1.05 );
    f_EoP_EE[i][j] -> FixParameter(2, 0.);
    
    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<10;trial++) {
      rp = h_EoP_EE[i][j] -> Fit(funcName, "QRLW+");
      fStatus = rp;
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
    }
    
    float flPhi = h_Phi_EE[i][j]->GetMean(); 
    
    if(i==0) g_EoP_EE[j] -> SetPoint(i, 0. , pow(f_EoP_EE[i][j]->GetParameter(1),2));
    else g_EoP_EE[j] -> SetPoint(i, int (flPhi) , pow(f_EoP_EE[i][j]->GetParameter(1),2));

    g_EoP_EE[j] -> SetPointError(i, 0., 2*f_EoP_EE[i][j]->GetParError(1));
    //cout  << " ***** " <<  1./f_EoP[i]->GetParameter(1) << " " << f_EoP[i]->GetParError(1) << " " << fStatus <<  endl; 

    //ratio preparation
    float rat = f_EoP_EE[i][j]->GetParameter(1);
    float era = f_EoP_EE[i][j]->GetParError(1); 
     
    // corrected    
    xNorm = h_EoC_EE[i][j]->Integral()/h_template_EE[mod.first][mod.second]->Integral() *
            h_EoC_EE[i][j]->GetBinWidth(1)/h_template_EE[mod.first][mod.second]->GetBinWidth(1); 
    
    sprintf(funcName,"f_EoC_%d_%d_Ref_%d_%d_EE",i,j,mod.first,mod.second);
    f_EoC_EE[i][j] = new TF1(funcName, templateHistoFuncEE[mod.first][mod.second], 0.2, 1.3, 3, "histoFunc");
    f_EoC_EE[i][j] -> SetParName(0,"Norm"); 
    f_EoC_EE[i][j] -> SetParName(1,"Scale factor"); 
    f_EoC_EE[i][j] -> SetLineWidth(1); 
    f_EoC_EE[i][j] -> SetLineColor(kGreen+2); 
    f_EoC_EE[i][j] -> SetNpx(10000);

    f_EoC_EE[i][j] -> FixParameter(0, xNorm);
    f_EoC_EE[i][j] -> SetParameter(1, 1.05 );
    f_EoC_EE[i][j] -> FixParameter(2, 0.);
    
    std::cout << "***** Re-Fitting ";
    for (int trial=0;trial<10;trial++) {
      rp = h_EoC_EE[i][j] -> Fit(funcName, "QRL+");
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
    }

    if(i==0) g_EoC_EE[j] -> SetPoint(i, 0., pow(f_EoC_EE[i][j]->GetParameter(1),2));
    else g_EoC_EE[j] -> SetPoint(i, int(flPhi), pow(f_EoC_EE[i][j]->GetParameter(1),2));
    g_EoC_EE[j] -> SetPointError(i, 0., 2*f_EoC_EE[i][j]->GetParError(1));
    //    cout << " ********** " <<  1./f_EoC[i]->GetParameter(1) << " " << f_EoC[i]->GetParError(1) << endl; 

    //ratio finalization
    rat /= f_EoC_EE[i][j]->GetParameter(1);
    era = rat*sqrt(era*era+f_EoC_EE[i][j]->GetParError(1)*f_EoC_EE[i][j]->GetParError(1)); 
    
    if(i==0) g_Rat_EE[j] -> SetPoint(i, 0. , rat);
    else  g_Rat_EE[j] -> SetPoint(i, int(flPhi) , rat);
    g_Rat_EE[j] -> SetPointError(i,  0. , era);
 
    g_Rat_EE[j]->SetLineColor(kBlue+2); 
  }
 }
  // Output 

  TFile* o = new TFile(outputFile.c_str(),"RECREATE");
  o -> cd();
  
  for(int i=0; i<nEtaBinsEB ; i++){
  TString Name;
  Name = Form("g_EoP_EB_%d",i);
  if(g_EoP_EB[i]->GetN()!=0) g_EoP_EB[i] -> Write(Name);
  Name = Form("g_EoC_EB_%d",i);
  if(g_EoC_EB[i]->GetN()!=0) g_EoC_EB[i] -> Write(Name);
  Name = Form("g_Rat_EB_%d",i);
  if(g_Rat_EB[i]->GetN()!=0) g_Rat_EB[i] -> Write(Name);
  }

  for(int i=0; i<nEtaBinsEE ; i++){
  TString Name;
  Name = Form("g_EoP_EE_%d",i);
  if(g_EoP_EE[i]->GetN()!=0) g_EoP_EE[i] -> Write(Name);
  Name = Form("g_EoC_EE_%d",i);
  if(g_EoC_EE[i]->GetN()!=0) g_EoC_EE[i] -> Write(Name);
  Name = Form("g_Rat_EE_%d",i);
  if(g_Rat_EE[i]->GetN()!=0) g_Rat_EE[i] -> Write(Name);
  }

  for (int imod = 0; imod<nPhiBinsTempEB; imod++){
    for (int imod2 = 0; imod2<nEtaBinsEB; imod2++){
    if(h_template_EB[imod][imod2]->GetEntries()!=0) h_template_EB[imod][imod2]->Write();
   }
  }

  for (int imod = 0; imod<nPhiBinsTempEE; imod++){
    for(int imod2 = 0; imod2<nEtaBinsEE; imod2++){
     if(h_template_EE[imod][imod2]->GetEntries()!=0) h_template_EE[imod][imod2]->Write();
   }
  }
 
  for (int imod = 0; imod<nPhiBinsEB; imod++){
   for(int imod2 = 0; imod2<nEtaBinsEB; imod2++){
    if(h_EoP_EB[imod][imod2]->GetEntries()!=0) h_EoP_EB[imod][imod2]->Write();
   }
  }

  for (int imod = 0; imod<nPhiBinsEB; imod++){
   for(int imod2 = 0; imod2<nEtaBinsEB; imod2++){
    if(h_EoC_EB[imod][imod2]->GetEntries()!=0) h_EoC_EB[imod][imod2]->Write();
   }
  }

  for (int imod = 0; imod<nPhiBinsEE; imod++){
   for(int imod2 = 0; imod2<nEtaBinsEE; imod2++){
    if(h_EoP_EE[imod][imod2]->GetEntries()!=0) h_EoP_EE[imod][imod2]->Write();
   }
  }

 
  for (int imod = 0; imod<nPhiBinsEE; imod++){
   for(int imod2 = 0; imod2<nEtaBinsEE; imod2++){
    if(h_EoC_EE[imod][imod2]->GetEntries()!=0) h_EoC_EE[imod][imod2]->Write();
   }
  }

  for(int imod =0; imod< nEtaBinsEB; imod++){
   if(h_phi_mc_EB[imod]->GetEntries()!=0) h_phi_mc_EB[imod]->Write();
   if(h_phi_data_EB[imod]->GetEntries()!=0) h_phi_data_EB[imod]->Write();
  }

  for(int imod =0; imod< nEtaBinsEE; imod++){
   if(h_phi_mc_EE[imod]->GetEntries()!=0) h_phi_mc_EE[imod]->Write(); 
    if(h_phi_data_EE[imod]->GetEntries()!=0) h_phi_data_EE[imod]->Write();
  }


h_et_mc->Write();
h_et_data->Write();

o -> Close();

return 0;
}

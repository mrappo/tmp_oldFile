// per compilare: g++ -Wall -o CalibrationMomentum `root-config --cflags --glibs` CalibrationMomentum.cpp

#include "../CommonTools/TEndcapRings.h"
#include "../CommonTools/histoFunc.h"
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

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

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

int templIndex(float eta){
    float feta = fabs(eta);
    if (feta <= 25)               {return 0;}
    if (feta>  25 && feta <=  45) {return 1;}
    if (feta>  45 && feta <=  65) {return 2;}
    if (feta>  65 && feta <=  85) {return 3;}
    if (feta>  85 && feta <=  98) {return 4;}
    if (feta>  98 && feta <= 108) {return 5;}
    if (feta> 108 && feta <= 118) {return 6;}
    if (feta> 118 )               {return 7;}

    return -1;
}

//**************  MAIN PROGRAM **************************************************************
int main(int argc, char** argv)
{
  //---- output file to save graphs
  char outfilename[100];
  sprintf(outfilename,"output/MomentumCalibration_vtx_min5_max100.root");
   
  //---- variables for selection
  float r9min = 0.00 ;
  float r9max = 9999 ;  
  float etaMax  = 2.5;
  float eta2Max = 2.5;
  float minPVz  = 5;
  float maxPVz  = 100;

  bool usePUweights = true;

  //--- weights for MC
  TFile weightsFile("CommonTools/weights/PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","READ"); 
  TH1F* hweights = (TH1F*)weightsFile.Get("hweights");

  float w[100];
  for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
    w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
  }
  weightsFile.Close();


  //----- NTUPLES--------------------
  TChain *ntu_DA = new TChain("ntu");
  TChain *ntu_MC = new TChain("ntu");

  // Data 
  ntu_DA->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-May10ReReco-v1_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v5_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("/data2/calibrator/NTUPLES/Run2011A/WZAnalysis/WZAnalysis_DoubleElectron_Run2011A-ZElectron-PromptSkim-v6_42XReReco_FT_R_42_V21B.root");
  ntu_DA->Add("/data2/calibrator/NTUPLES/Run2011B/WZAnalysis/WZAnalysis_DoubleElectron_Run2011B-ZElectron-PromptSkim-v1_42XReReco_FT_R_42_V21B.root");
 
  // --- MC Fall 2011
  ntu_MC->Add("/data2/calibrator/NTUPLES/Fall11/WZAnalysis/WZAnalysis_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1.root");


  std::cout << "     DATA: " << ntu_DA->GetEntries() << " entries in Data sample" << std::endl;
  std::cout << "     MC  : " << ntu_MC->GetEntries() << " entries in  MC  sample" << std::endl;

  // observables  
  int isW;
  float PV_z;
  float EoP, scEta, scPhi, mZ;
  float scEta2,scEne2,scPhi2;
  float scE3x3, scE5x5, scEne, scERaw, zVtx, scEt;  
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
  
  // Define number of phi bins
  Int_t iPhiBinsEB = 180; 
  Int_t nBinsEB = iPhiBinsEB/1 ;
  Int_t iPhiBinsEE = 40; 
  Int_t nBinsEE = iPhiBinsEE/1 ;
 
  std::cout << "nBinsEB = " << nBinsEB << std::endl;
  std::cout << "nBinsEE = " << nBinsEE << std::endl;

  // histogram definition
  TH1F** h_EoP_EB = new TH1F*[nBinsEB];   
  TH1F** h_EoC_EB = new TH1F*[nBinsEB]; 
  TH1F** h_Phi_EB = new TH1F*[nBinsEB]; // used to map iEta (as defined for Barrel and Endcap geom) into eta 
  TF1** f_EoP_EB = new TF1*[nBinsEB];
  TF1** f_EoC_EB = new TF1*[nBinsEB];
  std::vector<int> refIdEB(nBinsEB);

  TH1F** h_EoP_EE = new TH1F*[nBinsEE];   
  TH1F** h_EoC_EE = new TH1F*[nBinsEE];  
  TH1F** h_Phi_EE = new TH1F*[nBinsEE]; // used to map iEta (as defined for Barrel and Endcap geom) into eta 
  TF1** f_EoP_EE = new TF1*[nBinsEE];
  TF1** f_EoC_EE = new TF1*[nBinsEE];
  std::vector<int> refIdEE(nBinsEE);


  for(Int_t i = 0; i < nBinsEB; ++i)
  {
    TString histoName;
    histoName=Form("EoP_%d_EB", i);
    h_EoP_EB[i] = new TH1F(histoName, histoName, 900, 0., 2.);
    h_EoP_EB[i] -> SetFillColor(kRed+2);
    h_EoP_EB[i] -> SetLineColor(kRed+2);
    h_EoP_EB[i] -> SetFillStyle(3004);

    histoName=Form("EoC_%d_EB", i);
    h_EoC_EB[i] = new TH1F(histoName, histoName, 900, 0., 2.);
    h_EoC_EB[i] -> SetFillColor(kGreen+2);
    h_EoC_EB[i] -> SetLineColor(kGreen+2);
    h_EoC_EB[i] -> SetFillStyle(3004);
   
    histoName=Form("Phi_%d_EB", i);   
    h_Phi_EB[i] = new TH1F(histoName, histoName, 900, 0., 2.); 
  }

  for(Int_t i = 0; i < nBinsEE; ++i)
  {
    TString histoName;
    histoName=Form("EoP_%d_EE", i);
    h_EoP_EE[i] = new TH1F(histoName, histoName, 900, 0., 2.);
    h_EoP_EE[i] -> SetFillColor(kRed+2);
    h_EoP_EE[i] -> SetLineColor(kRed+2);
    h_EoP_EE[i] -> SetFillStyle(3004);

    histoName=Form("EoC_%d_EE", i);
    h_EoC_EE[i] = new TH1F(histoName, histoName, 900, 0., 2.);
    h_EoC_EE[i] -> SetFillColor(kGreen+2);
    h_EoC_EE[i] -> SetLineColor(kGreen+2);
    h_EoC_EE[i] -> SetFillStyle(3004);
   
    histoName=Form("Phi_%d_EE", i);   
    h_Phi_EE[i] = new TH1F(histoName, histoName, 1000, 0., 360.); 
  }

  Int_t nBinsTempEB = 180 ;
  Int_t nBinsTempEE = 40 ;
 
  TH1F* h_template_EB[nBinsTempEB];
  TH1F* h_template_EE[nBinsTempEE];

  for (Int_t imod = 0; imod<nBinsTempEB; imod++){
    TString histoName;
    histoName=Form("template_%d_EB", imod);
    h_template_EB[imod] = new TH1F(histoName, "", 900, 0., 2.);  
  }

  for (Int_t imod = 0; imod<nBinsTempEE; imod++){
    TString histoName;
    histoName=Form("template_%d_EE", imod);
    h_template_EE[imod] = new TH1F(histoName, "", 900, 0., 2.);  
  }

  
  TH1F* h_phi_data_EB = new TH1F("h_phi_data_EB","h_phi_data",100,-TMath::Pi(),TMath::Pi());
  TH1F* h_phi_mc_EB   = new TH1F("h_phi_mc_EB","h_phi_mc",100,-TMath::Pi(),TMath::Pi());
  TH1F* h_phi_data_EE = new TH1F("h_phi_data_EE","h_phi_data",100,-TMath::Pi(),TMath::Pi());
  TH1F* h_phi_mc_EE   = new TH1F("h_phi_mc_EE","h_phi_mc",100,-TMath::Pi(),TMath::Pi());
 
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
    if( fabs(scEta2) > eta2Max ) continue;
    float R9 = scE3x3/scEne;
    if ( R9 < r9min || R9 > r9max ) continue; 

    // PVz cut
    if (PV_z < minPVz || PV_z > maxPVz) continue;

    //--- PU weights
    if (usePUweights) ww = w[npu];
 
    //--- set ieta for the Endcaps
    int iphi;
    // //--- cut phi cracks
    // if (abs(ieta)<86) {
    //   float phi = (scPhi+PI)/xtalWidth;
    //   float modphi = (int)phi%20;
    //   if (fabs(modphi-10)<3.) continue;
    // }

    var = mZ * sqrt(pTK/scEne)/91.19;    /// use the momentum for ele1
    if(ele1_iz==0){
       int mod = int (ele1_iphi/(360./nBinsTempEB)); 
       if(mod == nBinsTempEB) h_template_EB[0] ->  Fill(var,ww);
       else h_template_EB[mod]->Fill(var,ww);
       // fill MC histos in eta bins
       int binEB = int (ele1_iphi/(360./nBinsEB)); 
 
       if(binEB==nBinsEB){ h_EoP_EB[0]-> Fill(var,ww);  // This is MC 
                           h_phi_mc_EB->Fill(scPhi*360,ww);}
       else{h_EoP_EB[binEB] -> Fill(var,ww);  // This is MC
            refIdEB.at(binEB) = mod; 
            h_phi_mc_EB->Fill(scPhi,ww);}

     
    }
    else{ iphi = eRings->GetEndcapIphi(ele1_ix,ele1_iy,ele1_iz);
          int mod = int (iphi/(360./nBinsTempEE)); 
          if(mod==nBinsTempEE) h_template_EE[0] ->  Fill(var,ww);
          else h_template_EE[mod] ->  Fill(var,ww);

          // fill MC histos in eta bins
          int binEE = int (iphi/(360./nBinsEE)); 
          if(binEE==nBinsEE){ h_EoP_EE[0] -> Fill(var,ww);  // This is MC
                              h_phi_mc_EE->Fill(scPhi,ww);}
          else{
               h_EoP_EE[binEE] -> Fill(var,ww);  // This is MC
               refIdEE.at(binEE) = mod; 
               h_phi_mc_EE->Fill(scPhi,ww);}

        }

    h_et_mc ->Fill(scEt,ww);
  } 
  // loop on events in Data
  std::cout << "Loop in Data events " << endl; 
  for(int entry = 0; entry < ntu_DA->GetEntries(); ++entry) {
    if( entry%10000 == 0 ) std::cout << "reading saved entry " << entry << "\r" << std::flush;
    //    if (entry>1000) break;

    ntu_DA->GetEntry(entry);

    if (isW==1) continue;
    if ( fabs(scEta2) > eta2Max) continue;
    float R9 = scE3x3/scEne;
    if ( R9 < r9min || R9 > r9max ) continue; 

    if (PV_z < minPVz || PV_z > maxPVz) continue;

    //--- set ieta for the Endcaps
    int iphi;
   
    // //--- cut phi cracks
    // if (abs(ieta)<86) {
    //   float phi = (scPhi+PI)/xtalWidth;
    //   float modphi = (int)phi%20;
    //   if (fabs(modphi-10)<3.) continue;
    // }
    var  = mZ * sqrt(pTK/scEne) / 91.19;    /// use the momentum for ele1
    //var = mZ/91.19;    /// use the SC energy for ele1
    if(ele1_iz==0){
       int binEB = int (ele1_iphi/(360./nBinsEB)); 
 
       if(binEB==nBinsEB){ h_EoC_EB[0]-> Fill(var);  // This is DATA
                           h_Phi_EB[0]->Fill((double) ele1_iphi); 
                           h_phi_data_EB->Fill(scPhi);
                          }
       else{h_EoC_EB[binEB] -> Fill(var);  // This is DATA
            h_Phi_EB[binEB]->Fill((double) ele1_iphi);
            h_phi_data_EB->Fill(scPhi);}
      
     }
    else{  iphi = eRings->GetEndcapIphi(ele1_ix,ele1_iy,ele1_iz); 
           int binEE = int (iphi/(360./nBinsEE)); 
           if(binEE==nBinsEE){
                              h_EoC_EE[0] -> Fill(var);  // This is DATA
                              h_Phi_EE[0] -> Fill((double) iphi); 
                              h_phi_data_EE ->Fill(scPhi);}
           else{
                  h_EoC_EE[binEE] -> Fill(var);  // This is DATA
                  h_Phi_EE[binEE] -> Fill((double) iphi); 
                  h_phi_data_EE ->Fill(scPhi);}
           
        }

    //use also the other electron
    var = mZ * sqrt(pTK2/scEne2) / 91.19;    /// use the momentum for ele2
    //var = mZ / 91.19;    /// use the SC energy for ele2

    // cut phi cracks in EB
    //    if (abs(ieta2)<86) {
    //     float phi = (scPhi2+3.1415926536)/0.01745329;
    //     float modphi = (int)phi%20;
    //     if (fabs(modphi-10)<3.) continue;
    //    }

    if(ele2_iz==0){

       int binEB = int (ele2_iphi/(360./nBinsEB)); 
           
       if(binEB==nBinsEB){ h_EoC_EB[0]-> Fill(var);  // This is DATA
                           h_Phi_EB[0]->Fill((double) ele2_iphi); 
                           h_phi_data_EB->Fill(scPhi);
                          }
       else{h_EoC_EB[binEB] -> Fill(var);  // This is DATA
            h_Phi_EB[binEB]->Fill((double) ele2_iphi);
            h_phi_data_EB->Fill(scPhi);}
      
     }
    else{     
           iphi = eRings->GetEndcapIphi(ele2_ix,ele2_iy,ele2_iz); 
           int binEE = int (iphi/(360./nBinsEE)); 
           if(binEE==nBinsEE){
                              h_EoC_EE[0] -> Fill(var);  // This is DATA
                              h_Phi_EE[0] -> Fill((double) iphi); 
                              h_phi_data_EE ->Fill(scPhi);}
           else{
                  h_EoC_EE[binEE] -> Fill(var);  // This is DATA
                  h_Phi_EE[binEE] -> Fill((double) iphi); 
                  h_phi_data_EE ->Fill(scPhi);}
           
        }

   
    h_et_data ->Fill(scEt);

  }

  std::cout << "End loop: Analyze events " << endl; 
    
  // draw results   
  TGraphErrors* g_EoP_EB   = new TGraphErrors();
  TGraphErrors* g_EoC_EB   = new TGraphErrors();
  TGraphErrors* g_Rat_EB   = new TGraphErrors();
  TGraphErrors* g_EoP_EE   = new TGraphErrors();
  TGraphErrors* g_EoC_EE   = new TGraphErrors();
  TGraphErrors* g_Rat_EE   = new TGraphErrors();
 
  int rebin = 10;
  
  histoFunc *templateHistoFuncEB[nBinsTempEB]; 
  histoFunc *templateHistoFuncEE[nBinsTempEE]; 

  for (int mod=0; mod<nBinsTempEB;mod++) {
    h_template_EB[mod] -> Rebin(rebin);
    templateHistoFuncEB[mod] = new histoFunc(h_template_EB[mod]);
  }

 for (int mod=0; mod<nBinsTempEE;mod++) {
    h_template_EE[mod] -> Rebin(rebin);
    templateHistoFuncEE[mod] = new histoFunc(h_template_EE[mod]);
  }

  // Template Fit in EB
  for(int i = 0; i < nBinsEB; ++i)
  {
    h_EoP_EB[i] -> Rebin(rebin);    
    h_EoC_EB[i] -> Rebin(rebin);    
    
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
    int mod = refIdEB.at(i); 
    
    char funcName[50];
    sprintf(funcName,"f_EoP_%d_Ref_%d_EB",i,mod);
    f_EoP_EB[i] = new TF1(funcName, templateHistoFuncEB[mod], 0.6, 1.3, 3, "histoFunc");
    f_EoP_EB[i] -> SetParName(0,"Norm"); 
    f_EoP_EB[i] -> SetParName(1,"Scale factor"); 
    f_EoP_EB[i] -> SetLineWidth(1); 
    f_EoP_EB[i] -> SetLineColor(kRed+2); 
    f_EoP_EB[i] -> SetNpx(10000);
    h_EoP_EB[i] -> Sumw2();
    // uncorrected    
    double xNorm = h_EoP_EB[i]->Integral()/h_template_EB[mod]->Integral() *
                   h_EoP_EB[i]->GetBinWidth(1)/h_template_EB[mod]->GetBinWidth(1); 

    f_EoP_EB[i] -> FixParameter(0, xNorm);
    f_EoP_EB[i] -> SetParameter(1, 1.05 );
    f_EoP_EB[i] -> FixParameter(2, 0.);
    
    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<10;trial++) {
      rp = h_EoP_EB[i] -> Fit(funcName, "QRLW+");
      fStatus = rp;
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
     }
    
    //    int bin = (ieta+iEtaBins/2) * (nBins/iEtaBins);  //questa da invertire   
    float flPhi = i * (iPhiBinsEB/nBinsEB) - iPhiBinsEB/2 + 0.5;
    flPhi = h_Phi_EB[i]->GetMean(); 
        
    // g_EoP -> SetPoint(i, flEta , 1./f_EoP[i]->GetParameter(1));
    // g_EoP -> SetPointError(i, 0., f_EoP[i]->GetParError(1));
    g_EoP_EB -> SetPoint(i, flPhi , pow(f_EoP_EB[i]->GetParameter(1),2));
    g_EoP_EB -> SetPointError(i, 0., 2*f_EoP_EB[i]->GetParError(1));
    //cout  << " ***** " <<  1./f_EoP[i]->GetParameter(1) << " " << f_EoP[i]->GetParError(1) << " " << fStatus <<  endl; 

    //ratio preparation
    float rat = f_EoP_EB[i]->GetParameter(1);
    float era = f_EoP_EB[i]->GetParError(1); 

    // corrected    
    xNorm = h_EoC_EB[i]->Integral()/h_template_EB[mod]->Integral() *
            h_EoC_EB[i]->GetBinWidth(1)/h_template_EB[mod]->GetBinWidth(1); 

    sprintf(funcName,"f_EoC_%d_Ref_%d_EB",i,mod);
    
    f_EoC_EB[i] = new TF1(funcName, templateHistoFuncEB[mod], 0.6, 1.3, 3, "histoFunc");
    f_EoC_EB[i] -> SetParName(0,"Norm"); 
    f_EoC_EB[i] -> SetParName(1,"Scale factor"); 
    f_EoC_EB[i] -> SetLineWidth(1); 
    f_EoC_EB[i] -> SetLineColor(kGreen+2); 
    f_EoC_EB[i] -> SetNpx(10000);

    f_EoC_EB[i] -> FixParameter(0, xNorm);
    f_EoC_EB[i] -> SetParameter(1, 1.05 );
    f_EoC_EB[i] -> FixParameter(2, 0.);
    
    std::cout << "***** Re-Fitting ";
    for (int trial=0;trial<10;trial++) {
      rp = h_EoC_EB[i] -> Fit(funcName, "QRL+");
      fStatus = rp;
      cout<<" Trial " <<trial<<endl;
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
    }

    // g_EoC -> SetPoint(i, flEta, 1./f_EoC[i]->GetParameter(1));
    // g_EoC -> SetPointError(i, 0., f_EoC[i]->GetParError(1));
    g_EoC_EB -> SetPoint(i, flPhi, pow(f_EoC_EB[i]->GetParameter(1),2));
    g_EoC_EB -> SetPointError(i, 0., 2*f_EoC_EB[i]->GetParError(1));
    //    cout << " ********** " <<  1./f_EoC[i]->GetParameter(1) << " " << f_EoC[i]->GetParError(1) << endl; 

    //ratio finalization
    rat /= f_EoC_EB[i]->GetParameter(1);
    //    rat = 1+2.*(rat-1);
    era = rat*sqrt(era*era+f_EoC_EB[i]->GetParError(1)*f_EoC_EB[i]->GetParError(1)); 
    
    g_Rat_EB -> SetPoint(i, flPhi , rat); 
    g_Rat_EB -> SetPointError(i,  0. , era); 
    g_Rat_EB->SetLineColor(kBlue+2); 
  }

 // Template Fit in EE

  for(int i = 0; i < nBinsEE; ++i)
  {
    h_EoP_EE[i] -> Rebin(rebin);    
    h_EoC_EE[i] -> Rebin(rebin);    
    
    // define the fitting function
    // N.B. [0] * ( [1] * f( [1]*(x-[2]) ) )
    int mod = refIdEE.at(i); 
    
    char funcName[50];
    sprintf(funcName,"f_EoP_%d_Ref_%d_EE",i,mod);
    f_EoP_EE[i] = new TF1(funcName, templateHistoFuncEE[mod], 0.6, 1.3, 3, "histoFunc");
    f_EoP_EE[i] -> SetParName(0,"Norm"); 
    f_EoP_EE[i] -> SetParName(1,"Scale factor"); 
    f_EoP_EE[i] -> SetLineWidth(1); 
    f_EoP_EE[i] -> SetLineColor(kRed+2); 
    f_EoP_EE[i] -> SetNpx(10000);
    f_EoP_EE[i] -> SetNpx(10000);
    h_EoP_EE[i] -> Sumw2();
    
    // uncorrected    
    double xNorm = h_EoP_EE[i]->Integral()/h_template_EE[mod]->Integral() *
                   h_EoP_EE[i]->GetBinWidth(1)/h_template_EE[mod]->GetBinWidth(1); 

    f_EoP_EE[i] -> FixParameter(0, xNorm);
    f_EoP_EE[i] -> SetParameter(1, 1.05 );
    f_EoP_EE[i] -> FixParameter(2, 0.);
    
    TFitResultPtr rp;
    int fStatus; 
    for (int trial=0;trial<10;trial++) {
      rp = h_EoP_EE[i] -> Fit(funcName, "QRLW+");
      fStatus = rp;
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
    }
    
    //    int bin = (ieta+iEtaBins/2) * (nBins/iEtaBins);  //questa da invertire   
    float flPhi = i * (iPhiBinsEE/nBinsEE) - iPhiBinsEE/2 + 0.5;
    flPhi = h_Phi_EE[i]->GetMean(); 
        
    // g_EoP -> SetPoint(i, flEta , 1./f_EoP[i]->GetParameter(1));
    // g_EoP -> SetPointError(i, 0., f_EoP[i]->GetParError(1));
    g_EoP_EE -> SetPoint(i, flPhi , pow(f_EoP_EE[i]->GetParameter(1),2));
    g_EoP_EE -> SetPointError(i, 0., 2*f_EoP_EE[i]->GetParError(1));
    //cout  << " ***** " <<  1./f_EoP[i]->GetParameter(1) << " " << f_EoP[i]->GetParError(1) << " " << fStatus <<  endl; 

    //ratio preparation
    float rat = f_EoP_EE[i]->GetParameter(1);
    float era = f_EoP_EE[i]->GetParError(1); 

    // corrected    
    xNorm = h_EoC_EE[i]->Integral()/h_template_EE[mod]->Integral() *
            h_EoC_EE[i]->GetBinWidth(1)/h_template_EE[mod]->GetBinWidth(1); 

    sprintf(funcName,"f_EoC_%d_Ref_%d_EE",i,mod);
    f_EoC_EE[i] = new TF1(funcName, templateHistoFuncEE[mod], 0.6, 1.3, 3, "histoFunc");
    f_EoC_EE[i] -> SetParName(0,"Norm"); 
    f_EoC_EE[i] -> SetParName(1,"Scale factor"); 
    f_EoC_EE[i] -> SetLineWidth(1); 
    f_EoC_EE[i] -> SetLineColor(kGreen+2); 
    f_EoC_EE[i] -> SetNpx(10000);

    f_EoC_EE[i] -> FixParameter(0, xNorm);
    f_EoC_EE[i] -> SetParameter(1, 1.05 );
    f_EoC_EE[i] -> FixParameter(2, 0.);
    
    std::cout << "***** Re-Fitting ";
    for (int trial=0;trial<10;trial++) {
      rp = h_EoC_EE[i] -> Fit(funcName, "QRL+");
      cout<<" Trial " <<trial<<endl;
      if (fStatus !=4) break;
      else if(trial==9) cout <<" No good Fit "<<endl;
  }

    // g_EoC -> SetPoint(i, flEta, 1./f_EoC[i]->GetParameter(1));
    // g_EoC -> SetPointError(i, 0., f_EoC[i]->GetParError(1));
    g_EoC_EE -> SetPoint(i, flPhi, pow(f_EoC_EE[i]->GetParameter(1),2));
    g_EoC_EE -> SetPointError(i, 0., 2*f_EoC_EE[i]->GetParError(1));
    //    cout << " ********** " <<  1./f_EoC[i]->GetParameter(1) << " " << f_EoC[i]->GetParError(1) << endl; 

    //ratio finalization
    rat /= f_EoC_EE[i]->GetParameter(1);
    //    rat = 1+2.*(rat-1);
    era = rat*sqrt(era*era+f_EoC_EE[i]->GetParError(1)*f_EoC_EE[i]->GetParError(1)); 
    
    g_Rat_EE -> SetPoint(i, flPhi , rat); 
    g_Rat_EE -> SetPointError(i,  0. , era); 
    g_Rat_EE->SetLineColor(kBlue+2); 
  }
  // Output 
  TFile* o = new TFile(outfilename,"RECREATE");
  o -> cd();
  
  g_EoP_EB -> Write("g_EoP_EB");
  g_EoC_EB -> Write("g_EoC_EB");
  g_Rat_EB -> Write("g_Rat_EB");
  g_EoP_EE -> Write("g_EoP_EE");
  g_EoC_EE -> Write("g_EoC_EE");
  g_Rat_EE -> Write("g_Rat_EE");
  
  for (int imod = 0; imod<nBinsTempEB; imod++){
    h_template_EB[imod]->Write();
  }
  

  for (int imod = 0; imod<nBinsTempEE; imod++){
    h_template_EE[imod]->Write();
  }
 
 for (int imod = 0; imod<nBinsEB; imod++){
     h_EoC_EB[imod]->Write();
  }
 
  for (int imod = 0; imod<nBinsEE; imod++){
    h_EoC_EE[imod]->Write();
  }
 
  h_phi_mc_EB->Write();
  h_phi_mc_EE->Write();
  h_phi_data_EB->Write();
  h_phi_data_EE->Write();

  h_et_mc->Write();
  h_et_data->Write();

  o -> Close();
 
return 0;
}
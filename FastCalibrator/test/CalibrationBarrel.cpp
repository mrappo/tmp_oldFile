/// Stand-alone program for Normalize IC EB by the mean on a eta ring + skipping xtal near dead channels and TT
/// in the normalization procedure
/// Folded Plots for Spread IC, Statistical Precision and spread
/// Correct IC near cracks and for momentum scale and produce txt IC values

#include <vector>
#include <utility>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "CalibrationUtils.h"

using namespace std;

int templIndexEB(float eta){
    float feta = fabs(eta);
    if (feta <= 25)               {return 0;}
    if (feta>  25 && feta <=  45) {return 0;}
    if (feta>  45 && feta <=  65) {return 0;}
    if (feta>  65 && feta <=  85) {return 0;}
    return -1;
}



/// Main function
int main(int argc, char **argv){


 // Set style options
 gROOT->Reset();
 gROOT->SetStyle("Plain");

 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);
 gStyle->SetOptTitle(1); 
 gStyle->SetOptStat(0); 
 gStyle->SetFitFormat("6.3g"); 
 gStyle->SetPalette(1); 
  
 gStyle->SetTextFont(42);
 gStyle->SetTextSize(0.05);
 gStyle->SetTitleFont(42,"xyz");
 gStyle->SetTitleSize(0.05);
 gStyle->SetLabelFont(42,"xyz");
 gStyle->SetLabelSize(0.05);
 gStyle->SetTitleXOffset(0.8);
 gStyle->SetTitleYOffset(1.1);
 gROOT->ForceStyle();

 /// Acquisition from cfg file
 
 if(argc != 2){
 std::cerr << ">>>>> analysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
 return 1;
 }

 parseConfigFile (argv[1]) ;
 
 std::string infile1  = gConfigParser -> readStringOption("Input::Inputfile1");
 std::string inputMomentumScale =  gConfigParser -> readStringOption("Input::inputMomentumScale");
 int evalStat = gConfigParser -> readIntOption("Input::evalStat");
 int nEtaBinsEB = gConfigParser -> readIntOption("Input::nEtaBinsEB");
 bool is2012Calib =  gConfigParser -> readBoolOption("Input::is2012Calib");

 std::string infile2, infile3;
 if(evalStat){
 infile2 = gConfigParser -> readStringOption("Input::Inputfile2");
 infile3 = gConfigParser -> readStringOption("Input::Inputfile3");
 }

 bool isMC = gConfigParser -> readBoolOption("Input::isMC");
 
 if ( infile1.empty()) {
    cout << " No input file specified !" << endl;
    return 1;
 }

 if ( evalStat && (infile2.empty() || infile3.empty() )){
    cout << " No input files to evaluate statistical precision specified !" << endl;
    return 1;
 }

 cout << "Making calibration plots for: " << infile1 << endl;
 
 std::string outputTxt = gConfigParser -> readStringOption("Output::outputTxt");
 std::string fileType = gConfigParser -> readStringOption("Output::fileType");
 std::string dirName = gConfigParser -> readStringOption("Output::dirName");
 bool printPlots = gConfigParser -> readBoolOption("Output::printPlots");

 TApplication* theApp = new TApplication("Application",&argc, argv);
 
 /// map for dead TT centre
 std::vector< std::pair<int,int> > TT_centre ;
 if(!is2012Calib) InitializeDeadTT_EB(TT_centre);
 if(is2012Calib)  InitializeDeadTT_EB2012(TT_centre);

 ///  Input file with full statistic
 
 TFile *f = new TFile(infile1.c_str());
 TH2F *h_scale_EB = (TH2F*)f->Get("h_scale_EB");
 TH2F *hcmap = (TH2F*) h_scale_EB->Clone("hcmap");
 hcmap -> Reset("ICEMS");
 hcmap -> ResetStats();

 NormalizeIC_EB(h_scale_EB, hcmap,TT_centre);
  
 ///-----------------------------------------------------------------
 ///--- Build the precision vs ieta plot starting from the TH2F of IC
 ///-----------------------------------------------------------------
  
 int ringGroupSize = 1;
 int nEtaRing = 85;
 TH1F **hspreadEtaFold= new TH1F*[nEtaRing];

 BookSpreadHistos_EB(hcmap,hspreadEtaFold,ringGroupSize,nEtaRing);

 /// Total spred Graph through Gaussian fit

 TGraphErrors *sigma_vs_EtaFold = new TGraphErrors();
 sigma_vs_EtaFold->SetMarkerStyle(20);
 sigma_vs_EtaFold->SetMarkerSize(1);
 sigma_vs_EtaFold->SetMarkerColor(kBlue+2);

 TF1 *fgaus = new TF1("fgaus","gaus",-10,10);
 TF1 *fgaus2 = new TF1("fgaus2","gaus",-100,100);

 int np = 0;
  
  for (int i = 1; i < nEtaRing+1; i++){
    float etaring = hcmap->GetYaxis()->GetBinCenter((ringGroupSize*i + ringGroupSize*(i-1))/2 + 1);
    float e     = 0.5*ringGroupSize;
    fgaus->SetParameter(1,1);
    fgaus->SetParameter(2,hspreadEtaFold[i-1]->GetRMS());
    fgaus->SetRange(1-5*hspreadEtaFold[i-1]->GetRMS(),1+5*hspreadEtaFold[i-1]->GetRMS());
    hspreadEtaFold[i-1]->Fit("fgaus","QRNME");

    sigma_vs_EtaFold-> SetPoint(np,fabs(etaring),fgaus->GetParameter(2));
    sigma_vs_EtaFold-> SetPointError(np,e,fgaus->GetParError(2));
    np++;
  }

  /// evaluate statistical precision and residual term
  
  TGraphErrors *statprecision_vs_EtaFold = new TGraphErrors();
  statprecision_vs_EtaFold->SetMarkerStyle(20);
  statprecision_vs_EtaFold->SetMarkerSize(1);
  statprecision_vs_EtaFold->SetMarkerColor(kRed+2);

  TGraphErrors *residual_vs_EtaFold = new TGraphErrors();
  residual_vs_EtaFold->SetMarkerStyle(20);
  residual_vs_EtaFold->SetMarkerSize(1);
  residual_vs_EtaFold->SetMarkerColor(kGreen+2);


  if (evalStat){

  /// Acqusition map from split-stat
  TFile *f2 = new TFile(infile2.c_str());
  TH2F *h_scale_EB_2 = (TH2F*)f2->Get("h_scale_EB");
  TH2F *hcmap2 = (TH2F*) h_scale_EB->Clone("hcmap2");
  hcmap2 -> Reset("ICEMS");
  hcmap2 -> ResetStats();

  NormalizeIC_EB(h_scale_EB_2, hcmap2,TT_centre);
  
  TFile *f3 = new TFile(infile3.c_str());
  TH2F *h_scale_EB_3 = (TH2F*)f3->Get("h_scale_EB");
  TH2F *hcmap3 = (TH2F*) h_scale_EB_3->Clone("hcmap3");
  hcmap3 -> Reset("ICEMS");
  hcmap3 -> ResetStats();
  
  NormalizeIC_EB(h_scale_EB_3, hcmap3,TT_centre);

 
  /// spread histos for statistical precision folding EB+ on EB-
  TH1F **hstatprecisionEtaFold= new TH1F*[nEtaRing];

  BookSpreadStatHistos_EB(hcmap2,hcmap3,hstatprecisionEtaFold,ringGroupSize,nEtaRing);

  np = 0;

  /// statical spread
  for (int i = 1; i < nEtaRing+1; i++){
      float etaring = hcmap2->GetYaxis()->GetBinCenter((ringGroupSize*i + ringGroupSize*(i-1))/2 + 1);
      float e     = 0.5*ringGroupSize;
      fgaus->SetParameter(1,1);
      fgaus->SetParameter(2,hstatprecisionEtaFold[i-1]->GetRMS());
      fgaus->SetRange(-5*hstatprecisionEtaFold[i-1]->GetRMS(),5*hstatprecisionEtaFold[i-1]->GetRMS());
      hstatprecisionEtaFold[i-1]->Fit("fgaus","QRNME");

      statprecision_vs_EtaFold-> SetPoint(np,fabs(etaring),fgaus->GetParameter(2));
      statprecision_vs_EtaFold-> SetPointError(np,e,fgaus->GetParError(2));
      np++;
   }
  
  ResidualSpread (statprecision_vs_EtaFold, sigma_vs_EtaFold, residual_vs_EtaFold);
 }
  /////////////////////////////////////////////////////////////
  /// Apply corrections for Momentum scale vs phi
  /////////////////////////////////////////////////////////////
  
  TFile *f4 = new TFile(inputMomentumScale.c_str());
 
  TGraphErrors** g_EoC_EB = new TGraphErrors* [nEtaBinsEB];

  for (int i =0 ; i< nEtaBinsEB ; i++){
  TString Name = Form ("g_EoC_EB_%d",i);
  g_EoC_EB[i] = (TGraphErrors*)f4->Get(Name);
  }

  TH2F *hcmapMomentumCorrected = (TH2F*) h_scale_EB->Clone("hcmap");
  
  hcmapMomentumCorrected -> Reset("ICEMS");
  hcmapMomentumCorrected -> ResetStats();
  
  TGraphErrors *phiProjection = new TGraphErrors();
  phiProjection->SetMarkerStyle(20);
  phiProjection->SetMarkerSize(1);
  phiProjection->SetMarkerColor(kBlack);

  TGraphErrors *phiCorrection = new TGraphErrors();
  phiCorrection->SetMarkerStyle(20);
  phiCorrection->SetMarkerSize(1);
  phiCorrection->SetMarkerColor(kGreen+2);


  /// For draw the projection value vs phi before and after correction
  PhiProfile(phiProjection,g_EoC_EB,hcmap);
  
  /// Correction of the map for momentum systematic
  for(int iPhi =1; iPhi<hcmap->GetNbinsX()+1; iPhi++){
   for(int iEta =1; iEta<hcmap->GetNbinsY()+1; iEta++){
     if(hcmap->GetBinContent(iPhi,iEta)==0) continue;
     double xPhi=0, yValue=0;
     int modEta = templIndexEB(iEta-85);
     if(modEta==-1) continue;
     g_EoC_EB[modEta]->GetPoint(int((iPhi-1)/(360./g_EoC_EB[0]->GetN())),xPhi,yValue);
     hcmapMomentumCorrected->SetBinContent(iPhi,iEta,hcmap->GetBinContent(iPhi,iEta)*yValue);
   }
  }

  /// Projection after momentum correction
  PhiProfile(phiCorrection,g_EoC_EB,hcmapMomentumCorrected);
   
  //////////////////////////////////////////////////////////////////
  ///Plot Folded %20 Phi for mean IC Normalized value and spread 
  //////////////////////////////////////////////////////////////////

  TGraphErrors *ic_vs_PhiFold_crack_EBp = new TGraphErrors();
  ic_vs_PhiFold_crack_EBp->SetMarkerStyle(20);
  ic_vs_PhiFold_crack_EBp->SetMarkerSize(1);
  ic_vs_PhiFold_crack_EBp->SetMarkerColor(kRed);
    
  TGraphErrors *ic_vs_PhiFold_crack_EBm = new TGraphErrors();
  ic_vs_PhiFold_crack_EBm->SetMarkerStyle(20);
  ic_vs_PhiFold_crack_EBm->SetMarkerSize(1);
  ic_vs_PhiFold_crack_EBm->SetMarkerColor(kBlue);

  TGraphErrors *spread_ic_vs_PhiFold_crack_EBp = new TGraphErrors();
  spread_ic_vs_PhiFold_crack_EBp->SetMarkerStyle(20);
  spread_ic_vs_PhiFold_crack_EBp->SetMarkerSize(1);
  spread_ic_vs_PhiFold_crack_EBp->SetMarkerColor(kRed);
    
  TGraphErrors *spread_ic_vs_PhiFold_crack_EBm = new TGraphErrors();
  spread_ic_vs_PhiFold_crack_EBm->SetMarkerStyle(20);
  spread_ic_vs_PhiFold_crack_EBm->SetMarkerSize(1);
  spread_ic_vs_PhiFold_crack_EBm->SetMarkerColor(kBlue);

  TH1F* hspreadPhiFold_crack_EBp[20];
  TH1F* hspreadPhiFold_crack_EBm[20];
  
  int nStep =0;
  char hname[100];
  int nbins =500;
    
  for(int jbin = 1; jbin < hcmapMomentumCorrected-> GetNbinsX()+1; jbin++){
   if (jbin <= 20) {
        nStep++;
        sprintf(hname,"hspread_iphiFolded_crack_EBp%02d",nStep);
        hspreadPhiFold_crack_EBp[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
        sprintf(hname,"hspread_iphiFolded_crack_EBm%02d",nStep);
        hspreadPhiFold_crack_EBm[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
      }

   for(int ibin = 1; ibin < hcmapMomentumCorrected-> GetNbinsY()+1; ibin++){
     float ic = hcmapMomentumCorrected->GetBinContent(jbin,ibin);
     
     bool isGood = CheckxtalIC_EB(hcmapMomentumCorrected,jbin,ibin);
     bool isGoodTT = CheckxtalTT_EB(jbin,ibin,TT_centre);

    if (ic>0 && ic<2 && isGood && isGoodTT ) {
      if(jbin <= 20) {
          if (ibin <= nEtaRing) hspreadPhiFold_crack_EBm[jbin-1]->Fill(ic); //from 1 to 85 included
          if (ibin >= nEtaRing+2) hspreadPhiFold_crack_EBp[jbin-1]->Fill(ic); //from 86 to 170 included
       }
       else{
          int kbin = (jbin-1)%20 ;
          if (ibin <= nEtaRing) hspreadPhiFold_crack_EBm[kbin]->Fill(ic);
          if (ibin >= nEtaRing+2) hspreadPhiFold_crack_EBp[kbin]->Fill(ic);
       }
      }
    }
  }

  np=0;
  /// Fit in each iphi%20 region (if you want you can use the mean of the distribution)
  for(int i=1; i<=20; i++){
      fgaus2->SetParameter(1,hspreadPhiFold_crack_EBp[i-1]->GetMean());
      fgaus2->SetParameter(2,hspreadPhiFold_crack_EBp[i-1]->GetRMS());

      fgaus2->SetRange(hspreadPhiFold_crack_EBp[i-1]->GetMean()-5*hspreadPhiFold_crack_EBp[i-1]->GetRMS(),
                       hspreadPhiFold_crack_EBp[i-1]->GetMean()+5*hspreadPhiFold_crack_EBp[i-1]->GetRMS());
      hspreadPhiFold_crack_EBp[i-1]->Fit("fgaus2","QRNME");

      ic_vs_PhiFold_crack_EBp-> SetPoint(np, i, fgaus2->GetParameter(1));
      ic_vs_PhiFold_crack_EBp-> SetPointError(np,0.5,fgaus2->GetParError(1));

      spread_ic_vs_PhiFold_crack_EBp -> SetPoint(np, i, fgaus2->GetParameter(2));
      spread_ic_vs_PhiFold_crack_EBp-> SetPointError(np,0.5,fgaus2->GetParError(2));

      fgaus2->SetParameter(1,hspreadPhiFold_crack_EBm[i-1]->GetMean());
      fgaus2->SetParameter(2,hspreadPhiFold_crack_EBm[i-1]->GetRMS());

      fgaus2->SetRange(hspreadPhiFold_crack_EBm[i-1]->GetMean()-5*hspreadPhiFold_crack_EBm[i-1]->GetRMS(),
                       hspreadPhiFold_crack_EBm[i-1]->GetMean()+5*hspreadPhiFold_crack_EBm[i-1]->GetRMS());
      hspreadPhiFold_crack_EBm[i-1]->Fit("fgaus2","QRNME");

      ic_vs_PhiFold_crack_EBm-> SetPoint(np, i, fgaus2->GetParameter(1));
      ic_vs_PhiFold_crack_EBm-> SetPointError(np,0.5,fgaus2->GetParError(1));
      
      spread_ic_vs_PhiFold_crack_EBm-> SetPoint(np, i, fgaus2->GetParameter(2));
      spread_ic_vs_PhiFold_crack_EBm-> SetPointError(np,0.5,fgaus2->GetParError(2));
      np++;
   }

  /// IC Correction for Crack structure 
  TF1 *pol0_EBp = new TF1("pol0_EBp","pol0",4,16);
  TF1 *pol0_EBm = new TF1("pol0_EBm","pol0",4,16);
  pol0_EBp->SetLineColor(kWhite);
  pol0_EBm->SetLineColor(kWhite);

  ic_vs_PhiFold_crack_EBp->Fit("pol0_EBp","QRNME");
  ic_vs_PhiFold_crack_EBm->Fit("pol0_EBm","QRNME");
  
  TH2F *hcmap_crackcorrected = (TH2F*) h_scale_EB->Clone("hcmap_crackcorrected");
  hcmap_crackcorrected->Reset("ICMES");
  hcmap_crackcorrected->ResetStats();

  /// crack corrected map
  for(int ibin =1 ; ibin <hcmapMomentumCorrected->GetNbinsX()+1 ; ibin++){
    for(int jbin =1; jbin <hcmapMomentumCorrected->GetNbinsY()+1 ; jbin++){
      float ic = hcmapMomentumCorrected->GetBinContent(ibin,jbin);
      int iPhi ;
      if(ibin<=20) iPhi=ibin-1;
      else iPhi = (ibin-1)%20;
      if (jbin <= nEtaRing){
      double ix,ic_crack;
      ic_vs_PhiFold_crack_EBm->GetPoint(iPhi,ix,ic_crack);
      if(iPhi>=16)
      hcmap_crackcorrected->SetBinContent(ibin,jbin,ic*pol0_EBm->GetParameter(0)/ic_crack);
      else hcmap_crackcorrected->SetBinContent(ibin,jbin,ic);
      }
    
      if (jbin >= nEtaRing+2){
      double ix,ic_crack;
      ic_vs_PhiFold_crack_EBp->GetPoint(iPhi,ix,ic_crack);
      if(iPhi<=4)
      hcmap_crackcorrected->SetBinContent(ibin,jbin,ic*pol0_EBp->GetParameter(0)/ic_crack);
      else hcmap_crackcorrected->SetBinContent(ibin,jbin,ic);
      }
    }
  }
  
  
  TGraphErrors *ic_vs_PhiFold_corrected_EBp = new TGraphErrors();
  ic_vs_PhiFold_corrected_EBp->SetMarkerStyle(20);
  ic_vs_PhiFold_corrected_EBp->SetMarkerSize(1);
  ic_vs_PhiFold_corrected_EBp->SetMarkerColor(kRed);
    
  TGraphErrors *ic_vs_PhiFold_corrected_EBm = new TGraphErrors();
  ic_vs_PhiFold_corrected_EBm->SetMarkerStyle(20);
  ic_vs_PhiFold_corrected_EBm->SetMarkerSize(1);
  ic_vs_PhiFold_corrected_EBm->SetMarkerColor(kBlue);
   

  TH1F* hspreadPhiFold_corrected_EBp[20];
  TH1F* hspreadPhiFold_corrected_EBm[20];
  nStep =0;
    
  for(int jbin = 1; jbin < hcmap_crackcorrected-> GetNbinsX()+1; jbin++){
     if (jbin <= 20) {
        nStep++;
        sprintf(hname,"hspreadPhiFold_corrected_EBp%02d",nStep);
        hspreadPhiFold_corrected_EBp[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
        sprintf(hname,"hspreadPhiFold_corrected_EBm%02d",nStep);
        hspreadPhiFold_corrected_EBm[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
       }
    for(int ibin = 1; ibin < hcmap_crackcorrected-> GetNbinsY()+1; ibin++){
     float ic = hcmap_crackcorrected->GetBinContent(jbin,ibin);
     bool isGood = CheckxtalIC_EB(hcmap_crackcorrected,jbin,ibin);
     bool isGoodTT = CheckxtalTT_EB(jbin,ibin,TT_centre);

     if (ic>0 && ic<2 && isGood && isGoodTT ) {
       if(jbin <= 20){
          if (ibin <= nEtaRing) hspreadPhiFold_corrected_EBm[jbin-1]->Fill(ic); //from 1 to 85 included
          if (ibin >= nEtaRing+2) hspreadPhiFold_corrected_EBp[jbin-1]->Fill(ic); //from 86 to 170 included
       }
       else{
          int kbin = (jbin-1)%20 ;
          if (ibin <= nEtaRing) hspreadPhiFold_corrected_EBm[kbin]->Fill(ic);
          if (ibin >= nEtaRing+2) hspreadPhiFold_corrected_EBp[kbin]->Fill(ic);
       }
      }
    }
  }

   np=0;
   for(int i=1; i<=20; i++) {
      fgaus2->SetParameter(1,hspreadPhiFold_corrected_EBp[i-1]->GetMean());
      fgaus2->SetParameter(2,hspreadPhiFold_corrected_EBp[i-1]->GetRMS());

      fgaus2->SetRange(hspreadPhiFold_corrected_EBp[i-1]->GetMean()-5*hspreadPhiFold_corrected_EBp[i-1]->GetRMS(),
                       hspreadPhiFold_corrected_EBp[i-1]->GetMean()+5*hspreadPhiFold_corrected_EBp[i-1]->GetRMS());
      hspreadPhiFold_corrected_EBp[i-1]->Fit("fgaus2","QRNME");

      ic_vs_PhiFold_corrected_EBp-> SetPoint(np, i, fgaus2->GetParameter(1));
      ic_vs_PhiFold_corrected_EBp-> SetPointError(np,0.5,fgaus2->GetParError(1));
  
      fgaus2->SetParameter(1,hspreadPhiFold_corrected_EBm[i-1]->GetMean());
      fgaus2->SetParameter(2,hspreadPhiFold_corrected_EBm[i-1]->GetRMS());

      fgaus2->SetRange(hspreadPhiFold_corrected_EBm[i-1]->GetMean()-5*hspreadPhiFold_corrected_EBm[i-1]->GetRMS(),
                       hspreadPhiFold_corrected_EBm[i-1]->GetMean()+5*hspreadPhiFold_corrected_EBm[i-1]->GetRMS());
      hspreadPhiFold_corrected_EBm[i-1]->Fit("fgaus2","QRNME");

      ic_vs_PhiFold_corrected_EBm-> SetPoint(np, i, fgaus2->GetParameter(1));
      ic_vs_PhiFold_corrected_EBm-> SetPointError(np,0.5,fgaus2->GetParError(1));
      
      np++;
   }

   /// Phi projection after crack correction
   TGraphErrors *PhiProjection_CrackCorrection = new TGraphErrors();
   PhiProjection_CrackCorrection->SetMarkerStyle(20);
   PhiProjection_CrackCorrection->SetMarkerSize(1);
   PhiProjection_CrackCorrection->SetMarkerColor(kRed);
   /// renormalize after correction
   TH2F *hcmapFinal = (TH2F*) hcmap_crackcorrected->Clone("hcmapFinal");
   hcmapFinal -> Reset("ICEMS");
   hcmapFinal -> ResetStats();

   NormalizeIC_EB(hcmap_crackcorrected, hcmapFinal,TT_centre,false);
   
   PhiProfile(PhiProjection_CrackCorrection,g_EoC_EB,hcmapFinal);

   /// Distribution of phi profile

   TH1F *Profile1 = new TH1F("Profile1","Profile1",100,0.97,1.03);
   TH1F *Profile2 = new TH1F("Profile2","Profile2",100,0.97,1.03);
   TH1F *Profile3 = new TH1F("Profile3","Profile3",100,0.97,1.03);

 
   for(int i=0; i<phiProjection->GetN() ; i++){
      double x=0,y=0;
      phiProjection->GetPoint(i,x,y);
      Profile1->Fill(y);
      phiCorrection->GetPoint(i,x,y);
      Profile2->Fill(y);
      PhiProjection_CrackCorrection->GetPoint(i,x,y);
      Profile3->Fill(y);
   }
 
  /// Spread after correction folding EB+ over EB-:
  TH1F **hspreadEtaFold2= new TH1F* [nEtaRing];
  nStep=0;

  BookSpreadHistos_EB(hcmapFinal, hspreadEtaFold2, ringGroupSize,nEtaRing);

  /// Total spred Graph through Gaussian fit
  TGraphErrors *sigma_vs_EtaFold_corrected = new TGraphErrors();
  sigma_vs_EtaFold_corrected->SetMarkerStyle(20);
  sigma_vs_EtaFold_corrected->SetMarkerSize(1);
  sigma_vs_EtaFold_corrected->SetMarkerColor(kBlue+2);
  np=0;

  for (int i = 1; i < 86; i++){
    float etaring = hcmap_crackcorrected->GetYaxis()->GetBinCenter((ringGroupSize*i + ringGroupSize*(i-1))/2 + 1);
    float e     = 0.5*ringGroupSize;
    fgaus->SetParameter(1,1);
    fgaus->SetParameter(2,hspreadEtaFold2[i-1]->GetRMS());
    fgaus->SetRange(1-5*hspreadEtaFold2[i-1]->GetRMS(),1+5*hspreadEtaFold2[i-1]->GetRMS());
    hspreadEtaFold2[i-1]->Fit("fgaus","QRNME");
    sigma_vs_EtaFold_corrected-> SetPoint(np,fabs(etaring),fgaus->GetParameter(2));
    sigma_vs_EtaFold_corrected-> SetPointError(np,e,fgaus->GetParError(2));
    np++;
  }
  
  TGraphErrors* residual_vs_EtaFold_Corrected = new TGraphErrors();
  residual_vs_EtaFold_Corrected->SetMarkerStyle(20);
  residual_vs_EtaFold_Corrected->SetMarkerSize(1);
  residual_vs_EtaFold_Corrected->SetMarkerColor(kGreen+2);

  if(evalStat) ResidualSpread (statprecision_vs_EtaFold, sigma_vs_EtaFold_corrected,residual_vs_EtaFold_Corrected);

  ///------------------------------------------------------------------------
  ///-----------------------------------------------------------------
  ///--- Draw plots
  ///-----------------------------------------------------------------
 
  TCanvas *c[30];

  c[0] = new TCanvas("hspreadEB","hspreadEB");
  c[0]->SetLeftMargin(0.1); 
  c[0]->SetRightMargin(0.13); 
  c[0]->SetGridx();
  
  h_scale_EB->GetXaxis()->SetNdivisions(1020);
  h_scale_EB->GetXaxis() -> SetLabelSize(0.03);
  h_scale_EB->GetXaxis() ->SetTitle("i#phi");
  h_scale_EB->GetYaxis() ->SetTitle("i#eta");
  h_scale_EB->GetZaxis() ->SetRangeUser(0.9,1.1);
  h_scale_EB->Draw("COLZ");
   
  
  c[1] = new TCanvas("hcmap","hcmap normalized");
  c[1]->SetGridx();
  c[1]->SetGridy();
  hcmap->GetXaxis()->SetNdivisions(1020);
  hcmap->GetXaxis() -> SetLabelSize(0.03);
  hcmap->GetXaxis() ->SetTitle("i#phi");
  hcmap->GetYaxis() ->SetTitle("i#eta");
  hcmap->GetZaxis() ->SetRangeUser(0.9,1.1);
  hcmap->Draw("colz");
 
  c[2] = new TCanvas("csigmaFold","csigmaFold");
  c[2]->SetGridx();
  c[2]->SetGridy();
  sigma_vs_EtaFold->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.10);
  sigma_vs_EtaFold->GetHistogram()->GetXaxis()-> SetRangeUser(0,85);
  sigma_vs_EtaFold->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
  sigma_vs_EtaFold->GetHistogram()->GetXaxis()-> SetTitle("|i#eta|");
  sigma_vs_EtaFold->Draw("ap");
  if (evalStat){
    statprecision_vs_EtaFold->Draw("psame");
    sigma_vs_EtaFold->Draw("psame");
    TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
    leg->SetFillColor(0);
    leg->AddEntry(statprecision_vs_EtaFold,"statistical precision", "LP");
    leg->AddEntry(sigma_vs_EtaFold,"spread", "LP");
    leg->Draw("same");
  }


  c[3] = new TCanvas("cresidualFold","cresidualFold");
  c[3]->SetGridx();
  c[3]->SetGridy();
  residual_vs_EtaFold->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.05);
  residual_vs_EtaFold->GetHistogram()->GetXaxis()-> SetRangeUser(0,85);
  residual_vs_EtaFold->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
  residual_vs_EtaFold->GetHistogram()->GetXaxis()-> SetTitle("|i#eta|");
  residual_vs_EtaFold->Draw("ap");

  c[4] = new TCanvas("ICPhiProjection","ICPhiProjection");
  c[4]->SetGridx();
  c[4]->SetGridy();
  phiProjection->GetYaxis()->SetRangeUser(0.97,1.03);
  phiProjection->GetXaxis()-> SetRangeUser(0,360);
  phiProjection->GetYaxis()-> SetTitle("#bar{IC}");
  phiProjection->GetXaxis()-> SetTitle("i#phi");
  phiProjection->Draw("ap");

  c[5] = new TCanvas("ICPhiProjectionCorrected","ICPhiProjectionCorrected");
  c[5]->SetGridx();
  c[5]->SetGridy();
  phiCorrection->GetYaxis()->SetRangeUser(0.97,1.03);
  phiCorrection->GetXaxis()-> SetRangeUser(0,360);
  phiCorrection->GetYaxis()-> SetTitle("#bar{IC}");
  phiCorrection->GetXaxis()-> SetTitle("i#phi");
  phiCorrection->Draw("ap");
 
  c[6] = new TCanvas("hcmapcorrected","hcmapcorrected");
  c[6]->SetGridx();
  c[6]->SetGridy();
  hcmapMomentumCorrected->GetXaxis()->SetNdivisions(1020);
  hcmapMomentumCorrected->GetXaxis() -> SetLabelSize(0.03);
  hcmapMomentumCorrected->GetXaxis() ->SetTitle("i#phi");
  hcmapMomentumCorrected->GetYaxis() ->SetTitle("i#eta");
  hcmapMomentumCorrected->GetZaxis() ->SetRangeUser(0.9,1.1);
  hcmapMomentumCorrected->Draw("colz");
 

  c[7] = new TCanvas("cphimeanfold_crack_EB+","cphimeanfold_crack_EB+");
  c[7]->SetGridx();
  c[7]->SetGridy();
  
  TLegend * legg1 = new TLegend(0.75,0.75,0.89, 0.89);
  legg1->AddEntry(ic_vs_PhiFold_crack_EBp,"EB+","LP");
  legg1->SetFillColor(0);
 
  ic_vs_PhiFold_crack_EBp->GetHistogram()->SetTitle(" Mean IC EB+");
  ic_vs_PhiFold_crack_EBp->GetHistogram()->GetYaxis()-> SetRangeUser(0.98,1.02);
  ic_vs_PhiFold_crack_EBp->GetHistogram()->GetXaxis()-> SetRangeUser(0.5,20.5);
  ic_vs_PhiFold_crack_EBp->GetHistogram()->GetYaxis()-> SetTitle("<IC>");
  ic_vs_PhiFold_crack_EBp->GetHistogram()->GetXaxis()-> SetTitle("i#phi%20");
  ic_vs_PhiFold_crack_EBp->GetHistogram()->SetTitle("EB+");
  ic_vs_PhiFold_crack_EBp->Draw("ap");
  legg1->Draw("same");

  
  c[8] = new TCanvas("cphimeanfold_crack_EB-","cphimeanfold_crackEB-");
  c[8]->cd();
  c[8]->SetGridx();
  c[8]->SetGridy();

  TLegend * legg2 = new TLegend(0.75,0.75,0.89, 0.89);
  legg2->AddEntry(ic_vs_PhiFold_crack_EBm,"EB-","LP");
  legg2->SetFillColor(0);
  
  ic_vs_PhiFold_crack_EBm->GetHistogram()->SetTitle(" Mean IC EB-");
  ic_vs_PhiFold_crack_EBm->GetHistogram()->GetYaxis()-> SetRangeUser(0.98,1.02);
  ic_vs_PhiFold_crack_EBm->GetHistogram()->GetXaxis()-> SetRangeUser(0.5,20.5);
  ic_vs_PhiFold_crack_EBm->GetHistogram()->SetTitle("EB-");
  ic_vs_PhiFold_crack_EBm->GetHistogram()->GetYaxis()-> SetTitle("<IC>");
  ic_vs_PhiFold_crack_EBm->GetHistogram()->GetXaxis()-> SetTitle("i#phi%20");
  ic_vs_PhiFold_crack_EBm->Draw("ap");
  legg2->Draw("same");

 /* TGraphErrors *ic_vs_PhiFold_crack_EBm_reflect = new TGraphErrors();
  ic_vs_PhiFold_crack_EBm_reflect->SetMarkerStyle(20);
  ic_vs_PhiFold_crack_EBm_reflect->SetMarkerSize(1);
  ic_vs_PhiFold_crack_EBm_reflect->SetMarkerColor(kBlue);
  
  /// Reflect EB+ over EB- iphi%20 plot
  int k=ic_vs_PhiFold_crack_EBm->GetN()-1;
  for(int iPoint=0; iPoint<ic_vs_PhiFold_crack_EBp->GetN(); iPoint++)
  { 
     double ix, iy;
     ic_vs_PhiFold_crack_EBm->GetPoint(k,ix,iy);
     double ex = ic_vs_PhiFold_crack_EBm ->GetErrorX(k);
     double ey = ic_vs_PhiFold_crack_EBm ->GetErrorY(k);
     ic_vs_PhiFold_crack_EBm_reflect->SetPoint(iPoint,iPoint+1,iy);
     ic_vs_PhiFold_crack_EBm_reflect->SetPointError(iPoint,ex,ey);
     
     k--;
  }
   c[12] = new TCanvas("cphimeanfold_crack_EB_ref","cphimeanfold_crackEB_ref");
   c[12]->cd();
   c[12]->SetGridx();
   c[12]->SetGridy();
   TLegend * legg = new TLegend(0.75,0.75,0.89, 0.89);
  
   legg->Clear();
   legg->AddEntry(ic_vs_PhiFold_crack_EBp,"EB+","LP");
   legg->AddEntry(ic_vs_PhiFold_crack_EBm_reflect,"EB-","LP");

   ic_vs_PhiFold_crack_EBp->Draw("apsame");
   ic_vs_PhiFold_crack_EBm_reflect->Draw("psame");
   legg->SetFillColor(0);
   legg->Draw("same");
   */

   c[9] = new TCanvas("hcmap_crackcorrected","hcmap_crackcorrected");
   c[9]->SetLeftMargin(0.1); 
   c[9]->SetRightMargin(0.13); 
   c[9]->SetGridx();
  
   hcmap_crackcorrected->GetXaxis()->SetNdivisions(1020);
   hcmap_crackcorrected->GetXaxis() -> SetLabelSize(0.03);
   hcmap_crackcorrected->GetXaxis() ->SetTitle("i#phi");
   hcmap_crackcorrected->GetYaxis() ->SetTitle("i#eta");
   hcmap_crackcorrected->GetZaxis() ->SetRangeUser(0.9,1.1);
   hcmap_crackcorrected->Draw("COLZ");

   c[10] = new TCanvas("hcmapFinal","hcmapFinal");
   c[10]->SetLeftMargin(0.1); 
   c[10]->SetRightMargin(0.13); 
   c[10]->SetGridx();
  
   hcmapFinal->GetXaxis()->SetNdivisions(1020);
   hcmapFinal->GetXaxis() -> SetLabelSize(0.03);
   hcmapFinal->GetXaxis() ->SetTitle("i#phi");
   hcmapFinal->GetYaxis() ->SetTitle("i#eta");
   hcmapFinal->GetZaxis() ->SetRangeUser(0.9,1.1);
   hcmapFinal->Draw("COLZ");
   
   
   c[11] = new TCanvas("cphimeanfold_corrected_EB+","cphimeanfold_crack_EB+");
   c[11]->SetGridx();
   c[11]->SetGridy();
   ic_vs_PhiFold_corrected_EBp->GetHistogram()->SetTitle(" Mean IC EB+");
   ic_vs_PhiFold_corrected_EBp->GetHistogram()->GetYaxis()-> SetRangeUser(0.95,1.05);
   ic_vs_PhiFold_corrected_EBp->GetHistogram()->GetXaxis()-> SetRangeUser(0.5,20.5);
   ic_vs_PhiFold_corrected_EBp->GetHistogram()->GetYaxis()-> SetTitle("mean IC");
   ic_vs_PhiFold_corrected_EBp->GetHistogram()->GetXaxis()-> SetTitle("i#phi%20");
   ic_vs_PhiFold_corrected_EBp->Draw("ap");

   c[12] = new TCanvas("cphimeanfold_corrected_EB-","cphimeanfold_crack_EB-");
   c[12]->SetGridx();
   c[12]->SetGridy();
   ic_vs_PhiFold_corrected_EBm->GetHistogram()->SetTitle(" Mean IC EB-");
   ic_vs_PhiFold_corrected_EBm->GetHistogram()->GetYaxis()-> SetRangeUser(0.95,1.05);
   ic_vs_PhiFold_corrected_EBm->GetHistogram()->GetXaxis()-> SetRangeUser(0.5,20.5);
   ic_vs_PhiFold_corrected_EBm->GetHistogram()->GetYaxis()-> SetTitle("mean IC");
   ic_vs_PhiFold_corrected_EBm->GetHistogram()->GetXaxis()-> SetTitle("i#phi%20");
   ic_vs_PhiFold_corrected_EBm->Draw("ap");

   c[13] = new TCanvas("PhiProjection_crackcorrected","PhiProjection_crackcorrected");
   c[13]->SetGridx();
   c[13]->SetGridy();
   PhiProjection_CrackCorrection->GetHistogram()->GetYaxis()-> SetRangeUser(0.97,1.03);
   PhiProjection_CrackCorrection->GetHistogram()->GetXaxis()-> SetRangeUser(0,360);
   PhiProjection_CrackCorrection->GetHistogram()->GetYaxis()-> SetTitle("Mean IC");
   PhiProjection_CrackCorrection->GetHistogram()->GetXaxis()-> SetTitle("i#phi");
   PhiProjection_CrackCorrection->Draw("ap");

   c[14] = new TCanvas("PhiProjection_same","PhiProjection_same");
   c[14]->SetGridx();
   c[14]->SetGridy();
   PhiProjection_CrackCorrection->GetXaxis()->SetRangeUser(100,200);
   PhiProjection_CrackCorrection->Draw("ap");
   phiProjection->Draw("psame");
 
   TLegend * legg3 = new TLegend(0.75,0.75,0.89, 0.89);
   legg3->AddEntry(phiProjection,"Original IC","LP");
   legg3->AddEntry(PhiProjection_CrackCorrection,"IC Crack Corrected","LP");
   legg3->SetFillColor(0);
   legg3->Draw("same");

   c[15] = new TCanvas("Profile1","Profile1");
   c[15]->SetGridx();
   c[15]->SetGridy();
   Profile1->GetXaxis()->SetTitle("#bar{IC}");
   Profile1->SetLineColor(kBlack);
   Profile1->SetMarkerSize(0.8);
   Profile1->SetLineWidth(2);
   Profile2->SetLineColor(kGreen+2);
   Profile2->SetMarkerSize(0.8);
   Profile2->SetLineWidth(2);

   Profile1->Draw();
   Profile2->Draw("same");

   TLegend * legg4 = new TLegend(0.75,0.75,0.89, 0.89);
   legg4->AddEntry(Profile1,"Original IC","LP");
   legg4->AddEntry(Profile2,"IC Momentum Corrected","LP");
   legg4->SetFillColor(0);
   legg4->Draw("same");


   c[16] = new TCanvas("Profile2","Profile2");
   c[16]->SetGridx();
   c[16]->SetGridy();
   Profile1->GetXaxis()->SetTitle("#bar{IC}");
   Profile1->SetLineColor(kBlack);
   Profile1->SetMarkerSize(0.8);
   Profile3->SetLineColor(kRed);
   Profile3->SetMarkerSize(0.8);
   Profile3->SetLineWidth(2);
    
   fgaus->SetParameter(1,1);
   fgaus->SetParameter(2,Profile1->GetRMS());
   fgaus->SetRange(1-5*Profile1->GetRMS(),1+5*Profile1->GetRMS());
   fgaus->SetLineColor(kBlack);
   Profile1->Fit("fgaus","QRME");
   cout<<" Mean Values : Uncorrected = "<<fgaus->GetParameter(1)<<" RMS = "<<fgaus->GetParameter(2)<<endl;
   
   fgaus->SetParameter(1,1);
   fgaus->SetParameter(2,Profile3->GetRMS());
   fgaus->SetRange(1-5*Profile3->GetRMS(),1+5*Profile3->GetRMS());
   fgaus->SetLineColor(kRed);
   Profile3->Fit("fgaus","QRME");
   cout<<" Mean Values : Corrected Crack = "<<Profile3->GetMean()<<" RMS "<<Profile3->GetRMS()<<endl;
   
   Profile1->Draw();
   Profile3->Draw("same");

   TLegend * legg5 = new TLegend(0.75,0.75,0.89, 0.89);
   legg5->AddEntry(Profile1,"Original IC","LP");
   legg5->AddEntry(Profile3,"IC Crack Corrected","LP");
   legg5->SetFillColor(0);
   legg5->Draw("same");
   
   

   c[17] = new TCanvas("csigmaFoldCorrected","csigmaFoldCorrected");
   c[17]->SetGridx();
   c[17]->SetGridy();
   sigma_vs_EtaFold_corrected->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.07);
   sigma_vs_EtaFold_corrected->GetHistogram()->GetXaxis()-> SetRangeUser(0,85);
   sigma_vs_EtaFold_corrected->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
   sigma_vs_EtaFold_corrected->GetHistogram()->GetXaxis()-> SetTitle("|i#eta|");
   sigma_vs_EtaFold_corrected->Draw("ap");
   if (evalStat){
    statprecision_vs_EtaFold->Draw("psame");
    sigma_vs_EtaFold->Draw("psame");
    TLegend * leg2 = new TLegend(0.6,0.7,0.89, 0.89);
    leg2->SetFillColor(0);
    leg2->AddEntry(statprecision_vs_EtaFold,"statistical precision", "LP");
    leg2->AddEntry(sigma_vs_EtaFold,"spread", "LP");
    leg2->Draw("same");
   }
   
   if(evalStat)
   {
    c[18] = new TCanvas("cResidualFoldCorrected","cResidualFoldCorrected");
    c[18]->SetGridx();
    c[18]->SetGridy();
    residual_vs_EtaFold_Corrected->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.04);
    residual_vs_EtaFold_Corrected->GetHistogram()->GetXaxis()-> SetRangeUser(0,85);
    residual_vs_EtaFold_Corrected->GetHistogram()->GetYaxis()-> SetTitle("residual");
    residual_vs_EtaFold_Corrected->GetHistogram()->GetXaxis()-> SetTitle("|i#eta|");
    residual_vs_EtaFold_Corrected->Draw("ap");
   }


   /// Dump IC in a txt file ---> IC from isolated electrons

   std::ofstream outTxt (outputTxt.c_str(),std::ios::out);
   outTxt << "---------------------------------------------------------------" << std::endl;
   outTxt << std::fixed << std::setprecision(0) << std::setw(10) << "iEta"
	  << std::fixed << std::setprecision(0) << std::setw(10) << "iPhi"
	  << std::fixed << std::setprecision(0) << std::setw(10) << "iZ"
	  << std::fixed << std::setprecision(6) << std::setw(15) << "IC"
	  << std::fixed << std::setprecision(6) << std::setw(15) << "error"
	  << std::endl;
   outTxt << "---------------------------------------------------------------" << std::endl;
   for (int iEta = 1; iEta < hcmapFinal->GetNbinsY()+1 ; iEta ++)
   {
     if (hcmapFinal->GetYaxis()->GetBinLowEdge(iEta) == 0) continue; //skip ieta=0

     double x,statPrec;
     if (hcmapFinal->GetYaxis()->GetBinLowEdge(iEta) < 0)
       statprecision_vs_EtaFold->GetPoint(int(fabs(hcmapFinal->GetYaxis()->GetBinLowEdge(iEta+85))),x,statPrec);  //mirroring of the folded precision
     else
       statprecision_vs_EtaFold->GetPoint(int(fabs(hcmapFinal->GetYaxis()->GetBinLowEdge(iEta-85))),x,statPrec);  //mirroring of the folded precision


     for (int iPhi =1 ; iPhi < hcmapFinal -> GetNbinsX()+1; iPhi++)
       {

	 outTxt << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinal->GetYaxis()->GetBinLowEdge(iEta)
		<< std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinal->GetXaxis()->GetBinLowEdge(iPhi) 
		<< std::fixed << std::setprecision(0) << std::setw(10) << "0"; //iz for the barrel

	   if(hcmapFinal->GetBinContent(iPhi,iEta) == 0.)
	     {
	       outTxt << std::fixed << std::setprecision(6) << std::setw(15) << "-1."
		      << std::fixed << std::setprecision(6) << std::setw(15) << "999."
		      << std::endl;
	     }
	   else
	     {
	       outTxt << std::fixed << std::setprecision(6) << std::setw(15) << hcmapFinal->GetBinContent(iPhi,iEta) 
		      << std::fixed << std::setprecision(6) << std::setw(15) << statPrec
		      << std::endl;
	     }
	 
       }
   }

theApp->Run();

return 0; 

}


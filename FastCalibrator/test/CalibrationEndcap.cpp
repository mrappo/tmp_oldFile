#include <vector>
#include <utility>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include "TFile.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TEndcapRings.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "TApplication.h"
#include "CalibrationUtils.h"


/// Stand-alone program to produce ECAL single electron calibration plots for EE
/// Input Files : MC or Data splistat and no splitstat, Momentum scale vs phi plot
/// Output Files :  StatPrec_MC_R9_EE.root --> stat precision MC usefull for CompareCalibMCTruth_EE.C only MC
//                  
using namespace std;

int templIndexEE(float eta){
    float feta = fabs(eta);
    if(eta<0){
    if (feta>  85 && feta <=  98) {return 0;}
    if (feta>  98 && feta <= 100) {return 0;}
    if (feta> 100 && feta <= 118) {return 0;}
    if (feta> 118 )               {return 0;}
    }
    else{
    if (feta>  85 && feta <=  98) {return 1;}
    if (feta>  98 && feta <= 100) {return 1;}
    if (feta> 100 && feta <= 118) {return 1;}
    if (feta> 118 )               {return 1;}
    }
    return -1;
}


/// Run ./bin/CalibrationEndacp.cpp cfg/calibrationEE_cfg.cfg

int main (int argc, char **argv)
{

  /// Set style options
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
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
  int evalStat = gConfigParser -> readIntOption("Input::evalStat");
  std::string infile2,infile3; 
  if(evalStat){
   infile2= gConfigParser -> readStringOption("Input::Inputfile2");
   infile3 = gConfigParser -> readStringOption("Input::Inputfile3");
  }

  std::string inputMomentumScale =  gConfigParser -> readStringOption("Input::inputMomentumScale");
  std::string SystematicToAdd =  gConfigParser -> readStringOption("Input::SystematicToAdd"); 
  bool isMC = gConfigParser -> readBoolOption("Input::isMC");
  bool is2012Calib = gConfigParser -> readBoolOption("Input::is2012Calib");

  int nEtaBinsEE = gConfigParser -> readIntOption("Input::nEtaBinsEE");

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
  std::string fileName  = gConfigParser -> readStringOption("Output::fileName");
  
  TApplication* theApp = new TApplication("Application",&argc, argv);
  
  /// imput file with full statistic normlized to the mean in a ring

  TFile *f = new TFile(infile1.c_str());
  TH2F **hcScale = new TH2F*[2];
  TH2F **hcmap = new TH2F*[2];

  std::vector< std::pair<int,int> > TT_centre_EEP; 
  std::vector< std::pair<int,int> > TT_centre_EEM;
 
  if(!is2012Calib) InitializeDeadTTEEP(TT_centre_EEP);
  if(!is2012Calib) InitializeDeadTTEEM(TT_centre_EEM);

  if(is2012Calib) InitializeDeadTTEEP2012(TT_centre_EEP);
  if(is2012Calib) InitializeDeadTTEEM2012(TT_centre_EEM);


  hcScale[0] = (TH2F*)f->Get("h_scale_EEM");
  hcScale[1] = (TH2F*)f->Get("h_scale_EEP");
 
  hcmap[0] = (TH2F*) hcScale[0]->Clone("hcmapEEM");
  hcmap[1] = (TH2F*) hcScale[1]->Clone("hcmapEEP");
  hcmap[0]->Reset();
  hcmap[1]->Reset();

  TEndcapRings *eRings = new TEndcapRings(); 
  NormalizeIC_EE(hcScale, hcmap, TT_centre_EEP,TT_centre_EEM, eRings);
  
  ///--------------------------------------------------------------------------------
  ///--- Build the precision vs ring plot starting from the TH2F of IC folded and not
  ///--------------------------------------------------------------------------------

  TH1F ***hspread = new TH1F**[2];
  for( int i =0 ; i<2 ; i ++)   hspread[i] = new TH1F*[40];
 
  TH1F **hspreadAll = new TH1F*[40];
 
  /// ring geometry for the endcap
  
  BookSpreadHistos_EE(hcmap,hspread,hspreadAll,eRings);
  /// Graph Error for spread EE+ and EE-

  TGraphErrors *sigma_vs_ring[3];
  sigma_vs_ring[0] = new TGraphErrors();
  sigma_vs_ring[0]->SetMarkerStyle(20);
  sigma_vs_ring[0]->SetMarkerSize(1);
  sigma_vs_ring[0]->SetMarkerColor(kBlue+2);
  
  sigma_vs_ring[1] = new TGraphErrors();
  sigma_vs_ring[1]->SetMarkerStyle(20);
  sigma_vs_ring[1]->SetMarkerSize(1);
  sigma_vs_ring[1]->SetMarkerColor(kBlue+2);

  sigma_vs_ring[2] = new TGraphErrors();
  sigma_vs_ring[2]->SetMarkerStyle(20);
  sigma_vs_ring[2]->SetMarkerSize(1);
  sigma_vs_ring[2]->SetMarkerColor(kBlue+2);
 
  /// Graph for scale vs ring EE+, EE- and folded

  TGraphErrors *scale_vs_ring[3];
  scale_vs_ring[0] = new TGraphErrors();
  scale_vs_ring[0]->SetMarkerStyle(20);
  scale_vs_ring[0]->SetMarkerSize(1);
  scale_vs_ring[0]->SetMarkerColor(kBlue+2);

  scale_vs_ring[1] = new TGraphErrors();
  scale_vs_ring[1]->SetMarkerStyle(20);
  scale_vs_ring[1]->SetMarkerSize(1);
  scale_vs_ring[1]->SetMarkerColor(kBlue+2);

  scale_vs_ring[2] = new TGraphErrors();
  scale_vs_ring[2]->SetMarkerStyle(20);
  scale_vs_ring[2]->SetMarkerSize(1);
  scale_vs_ring[2]->SetMarkerColor(kBlue+2);
    
  TF1 *fgaus = new TF1("fgaus","gaus",-10,10);
  int np[3] = {0};
  
  
  /// Gaussian fit for EE+ and EE- ---> not folded
  for (int k = 0; k < 2 ; k++){
    for (int iring = 0; iring < 40; iring++){
      if (hspread[k][iring]-> GetEntries() == 0) continue;
      float e     = 0.5*hcmap[k]-> GetYaxis()->GetBinWidth(1);
      fgaus->SetParameter(1,1);
      fgaus->SetParameter(2,hspread[k][iring]->GetRMS());
      fgaus->SetRange(1-5*hspread[k][iring]->GetRMS(),1+5*hspread[k][iring]->GetRMS());
      hspread[k][iring]->Fit("fgaus","QRN");
      sigma_vs_ring[k]-> SetPoint(np[k],iring,fgaus->GetParameter(2));
      sigma_vs_ring[k]-> SetPointError(np[k], e ,fgaus->GetParError(2));
      scale_vs_ring[k]-> SetPoint(np[k],iring,fgaus->GetParameter(1));
      scale_vs_ring[k]-> SetPointError(np[k],e,fgaus->GetParError(1));
      np[k]++;    
    }
  }

  /// Folded spread distribution
  for (int iring = 0; iring < 40; iring++){
    if (hspreadAll[iring]-> GetEntries() == 0) continue;
    float e     = 0.5*hcmap[0]-> GetYaxis()->GetBinWidth(1);
      fgaus->SetParameter(1,1);
      fgaus->SetParameter(2,hspreadAll[iring]->GetRMS());
      fgaus->SetRange(1-5*hspreadAll[iring]->GetRMS(),1+5*hspreadAll[iring]->GetRMS());
      hspreadAll[iring]->Fit("fgaus","QRN");
      sigma_vs_ring[2]-> SetPoint(np[2],iring,fgaus->GetParameter(2));
      sigma_vs_ring[2]-> SetPointError(np[2], e ,fgaus->GetParError(2));
      scale_vs_ring[2]-> SetPoint(np[2],iring,fgaus->GetParameter(1));
      scale_vs_ring[2]-> SetPointError(np[2],e,fgaus->GetParError(1));
      np[2]++;    
    }
 ///----------------- Statistical Precision  and Residual --------------------

  TGraphErrors *statprecision_vs_ring[3];
  statprecision_vs_ring[0]  = new TGraphErrors();
  statprecision_vs_ring[0]->SetMarkerStyle(20);
  statprecision_vs_ring[0]->SetMarkerSize(1);
  statprecision_vs_ring[0]->SetMarkerColor(kRed+2);

  statprecision_vs_ring[1]  = new TGraphErrors();
  statprecision_vs_ring[1]->SetMarkerStyle(20);
  statprecision_vs_ring[1]->SetMarkerSize(1);
  statprecision_vs_ring[1]->SetMarkerColor(kRed+2);

  statprecision_vs_ring[2]  = new TGraphErrors();
  statprecision_vs_ring[2]->SetMarkerStyle(20);
  statprecision_vs_ring[2]->SetMarkerSize(1);
  statprecision_vs_ring[2]->SetMarkerColor(kRed+2);
    
  TGraphErrors *residual_vs_ring[3];
  residual_vs_ring[0] = new TGraphErrors();
  residual_vs_ring[0]->SetMarkerStyle(20);
  residual_vs_ring[0]->SetMarkerSize(1);
  residual_vs_ring[0]->SetMarkerColor(kGreen+2);
  
  residual_vs_ring[1] = new TGraphErrors();
  residual_vs_ring[1]->SetMarkerStyle(20);
  residual_vs_ring[1]->SetMarkerSize(1);
  residual_vs_ring[1]->SetMarkerColor(kGreen+2);
  
  residual_vs_ring[2] = new TGraphErrors();
  residual_vs_ring[2]->SetMarkerStyle(20);
  residual_vs_ring[2]->SetMarkerSize(1);
  residual_vs_ring[2]->SetMarkerColor(kGreen+2);

  if (evalStat){
    
    /// acquisition file for statistical precision

   TFile *f2 = new TFile(infile2.c_str());
   TFile *f3 = new TFile(infile3.c_str());
   TH2F **hcScale2 = new TH2F*[2];
   TH2F **hcmap2 = new TH2F*[2];
   hcScale2[0] = (TH2F*)f2->Get("h_scale_EEM"); 
   hcScale2[1] = (TH2F*)f2->Get("h_scale_EEP");
 
   TH2F **hcScale3 = new TH2F*[2];
   TH2F **hcmap3 = new TH2F*[2];
   hcScale3[0] = (TH2F*)f3->Get("h_scale_EEM"); 
   hcScale3[1] = (TH2F*)f3->Get("h_scale_EEP");

   hcmap2[0] = (TH2F*) hcScale2[0]->Clone("hcmapEEM2");
   hcmap2[1] = (TH2F*) hcScale2[1]->Clone("hcmapEEP2");
   hcmap2[0]->Reset();
   hcmap2[1]->Reset();

   hcmap3[0] = (TH2F*) hcScale3[0]->Clone("hcmapEEM3");
   hcmap3[1] = (TH2F*) hcScale3[1]->Clone("hcmapEEP3");
   hcmap3[0]->Reset();
   hcmap3[1]->Reset();

   NormalizeIC_EE(hcScale2, hcmap2, TT_centre_EEP,TT_centre_EEM, eRings);
 
   NormalizeIC_EE(hcScale3, hcmap3, TT_centre_EEP,TT_centre_EEM, eRings);
 

   TH1F ***hstatprecision = new TH1F**[2];
   for(int i =0; i<2; i++) hstatprecision[i]=new TH1F*[40];

   TH1F **hstatprecisionAll= new TH1F*[40];
   
   BookSpreadStatHistos_EE(hcmap2,hcmap3,hstatprecision,hstatprecisionAll,eRings);
   
   /// Gaussian fit of the even/odd distribution (rms of the distribution can be also used)
   int n[3] = {0};
   for (int k = 0; k < 2; k++){
      for (int iring = 0; iring < 40 ; iring++){
	if ( hstatprecision[k][iring]->GetEntries() == 0) continue;
	float e     = 0.5*hcmap2[k]-> GetYaxis()->GetBinWidth(1);
	fgaus->SetParameter(1,1);
	fgaus->SetParameter(2,hstatprecision[k][iring]->GetRMS());
	fgaus->SetRange(-5*hstatprecision[k][iring]->GetRMS(),5*hstatprecision[k][iring]->GetRMS());
	TString name = Form("ff%d_%d",iring,k);

        hstatprecision[k][iring]->Fit("fgaus","QRN");
        statprecision_vs_ring[k]-> SetPoint(n[k],iring,fgaus->GetParameter(2));
	statprecision_vs_ring[k]-> SetPointError(n[k],e,fgaus->GetParError(2));
	n[k]++;
      }
    }
    
    for (int iring = 0; iring < 40 ; iring++){
	if ( hstatprecisionAll[iring]->GetEntries() == 0) continue;
	float e     = 0.5*hcmap2[0]-> GetYaxis()->GetBinWidth(1);
	fgaus->SetParameter(1,1);
	fgaus->SetParameter(2,hstatprecisionAll[iring]->GetRMS());
	fgaus->SetRange(-5*hstatprecisionAll[iring]->GetRMS(),5*hstatprecisionAll[iring]->GetRMS());
	TString name = Form("ffAll%d",iring);
        hstatprecisionAll[iring]->Fit("fgaus","QRN");
      
	statprecision_vs_ring[2]-> SetPoint(n[2],iring,fgaus->GetParameter(2));
	statprecision_vs_ring[2]-> SetPointError(n[2],e,fgaus->GetParError(2));
	n[2]++;
      }
     
    /// Residual spread plot
    for (int k = 0; k < 3 ; k++) ResidualSpread (statprecision_vs_ring[k], sigma_vs_ring[k], residual_vs_ring[k]);
      
  }
 /// Momentum scale correction
 
 TFile* input = new TFile(inputMomentumScale.c_str());
 
 TGraphErrors** g_EoC_EE = new TGraphErrors* [nEtaBinsEE];

 for (int i =0 ; i< nEtaBinsEE ; i++){
  TString Name = Form ("g_EoC_EE_%d",i);
  g_EoC_EE[i] = (TGraphErrors*)input->Get(Name);
 }

 TH2F* mapConversionEEp = (TH2F*) input->Get("mapConversionEEp");
 TH2F* mapConversionEEm = (TH2F*) input->Get("mapConversionEEm");

 TGraphErrors* PhiProjectionEEp = new TGraphErrors();
 TGraphErrors* PhiProjectionEEm = new TGraphErrors();
 
 PhiProjectionEEp->SetMarkerStyle(20);
 PhiProjectionEEp->SetMarkerSize(1);
 PhiProjectionEEp->SetMarkerColor(kRed);

 PhiProjectionEEm->SetMarkerStyle(20);
 PhiProjectionEEm->SetMarkerSize(1);
 PhiProjectionEEm->SetMarkerColor(kBlue);
 
 TGraphErrors* PhiProjectionEEp_Corrected = new TGraphErrors();
 TGraphErrors* PhiProjectionEEm_Corrected = new TGraphErrors();

 PhiProjectionEEp_Corrected->SetMarkerStyle(20);
 PhiProjectionEEp_Corrected->SetMarkerSize(1);
 PhiProjectionEEp_Corrected->SetMarkerColor(kRed);

 PhiProjectionEEm_Corrected->SetMarkerStyle(20);
 PhiProjectionEEm_Corrected->SetMarkerSize(1);
 PhiProjectionEEm_Corrected->SetMarkerColor(kBlue);

 PhiProfileEE(PhiProjectionEEm, g_EoC_EE, hcmap[0],eRings,-1);
 
 PhiProfileEE(PhiProjectionEEp, g_EoC_EE, hcmap[1],eRings,1);

 /// Correction EE+ and EE-

 TH2F* mapMomentumCorrected[2];
 mapMomentumCorrected[0] = (TH2F*) hcmap[0]->Clone("mapMomentumCorrected_EEM");
 mapMomentumCorrected[1] = (TH2F*) hcmap[1]->Clone("mapMomentumCorrected_EEP");
 mapMomentumCorrected[0]->Reset();
 mapMomentumCorrected[1]->Reset();

 for(int ix=1; ix<hcmap[0]->GetNbinsX()+1;ix++){
   for(int iy=1; iy<hcmap[0]->GetNbinsY()+1;iy++){
    if(hcmap[0]->GetBinContent(ix,iy)==0) continue;
     int iPhi = int(eRings->GetEndcapIphi(ix,iy,-1));
     int modEta = templIndexEE(int(eRings->GetEndcapIeta(ix,iy,-1)));
     double xphi,yphi;
     if(modEta == -1) continue;
     g_EoC_EE[modEta]->GetPoint(int(iPhi/(360./PhiProjectionEEm->GetN())),xphi,yphi);
     yphi = 1.;
     mapMomentumCorrected[0]->SetBinContent(ix,iy,hcmap[0]->GetBinContent(ix,iy)*yphi);
   }
  }

  for(int ix=1; ix<hcmap[1]->GetNbinsX()+1;ix++){
   for(int iy=1; iy<hcmap[1]->GetNbinsY()+1;iy++){
    if(hcmap[1]->GetBinContent(ix,iy)==0) continue;
     int iPhi = int(eRings->GetEndcapIphi(ix,iy,1));
     int modEta = templIndexEE(int(eRings->GetEndcapIeta(ix,iy,-1)));
     if(modEta == -1) continue;
     double xphi,yphi; 
     g_EoC_EE[modEta]->GetPoint(int(iPhi/(360./PhiProjectionEEp->GetN())),xphi,yphi);
     yphi = 1.;
     mapMomentumCorrected[1]->SetBinContent(ix,iy,hcmap[1]->GetBinContent(ix,iy)*yphi);
   }
  }
 /// New Normalization after momentum scale correction

 TH2F** hcmapFinalEE= new TH2F*[2];

 hcmapFinalEE[0] = (TH2F*) mapMomentumCorrected[1]->Clone("hcmapFinalEEp");
 hcmapFinalEE[1] = (TH2F*) mapMomentumCorrected[0]->Clone("hcmapFinalEEm");
  
 hcmapFinalEE[0] -> Reset("ICEMS");
 hcmapFinalEE[1] -> Reset("ICEMS");
 hcmapFinalEE[0] -> ResetStats();
 hcmapFinalEE[1] -> ResetStats();

 NormalizeIC_EE(mapMomentumCorrected, hcmapFinalEE, TT_centre_EEP,TT_centre_EEM, eRings,false);


 /// EE+ and EE- projection after correction

 PhiProfileEE(PhiProjectionEEm_Corrected, g_EoC_EE, hcmapFinalEE[0],eRings,-1);
 
 PhiProfileEE(PhiProjectionEEp_Corrected, g_EoC_EE, hcmapFinalEE[1],eRings,1);

 /// Projection Histos :
 TH1F* ProfileEEp = new TH1F ("ProfileEEp","ProfileEEp",60,0.9,1.1);
 TH1F* ProfileEEm = new TH1F ("ProfileEEm","ProfileEEm",60,0.9,1.1);

 TH1F* ProfileEEpCorrected = new TH1F ("ProfileEEpCorrected","ProfileEEpCorrected",60,0.9,1.1);
 TH1F* ProfileEEmCorrected = new TH1F ("ProfileEEmCorrected","ProfileEEmCorrected",60,0.9,1.1);

for( int i=0; i<PhiProjectionEEm->GetN(); i++){
   double x,y;
   PhiProjectionEEm->GetPoint(i,x,y);
   ProfileEEm->Fill(y);
 }

 for( int i=0; i<PhiProjectionEEp->GetN(); i++){
   double x,y;
   PhiProjectionEEp->GetPoint(i,x,y);
   ProfileEEp->Fill(y);
 }

 fgaus->SetParameter(1,1);
 fgaus->SetParameter(2,ProfileEEm->GetRMS());
 fgaus->SetRange(ProfileEEm->GetMean()-5*ProfileEEm->GetRMS(),ProfileEEm->GetMean()+5*ProfileEEm->GetRMS());
 fgaus->SetLineColor(kBlack);
 ProfileEEm->SetLineWidth(2);
 ProfileEEm->Fit("fgaus","QRME");

 cout<<" EEm Uncorrected: Mean = "<<fgaus->GetParameter(1)<<" Sigma "<<fgaus->GetParameter(2)<<endl;

 fgaus->SetParameter(1,1);
 fgaus->SetParameter(2,ProfileEEp->GetRMS());
 fgaus->SetRange(ProfileEEp->GetMean()-5*ProfileEEp->GetRMS(),ProfileEEp->GetMean()+5*ProfileEEm->GetRMS());
 ProfileEEp->SetLineWidth(kBlack);
 ProfileEEp->SetLineWidth(2);
 ProfileEEp->Fit("fgaus","QRME");

 cout<<" EEp Uncorrected: Mean = "<<fgaus->GetParameter(1)<<" Sigma "<<fgaus->GetParameter(2)<<endl;

 for( int i=0; i<PhiProjectionEEm_Corrected->GetN(); i++){
   double x,y;
   PhiProjectionEEm_Corrected->GetPoint(i,x,y);
   ProfileEEmCorrected->Fill(y);
 }

 for( int i=0; i<PhiProjectionEEp_Corrected->GetN(); i++){
   double x,y;
   PhiProjectionEEp_Corrected->GetPoint(i,x,y);
   ProfileEEpCorrected->Fill(y);
 }

 fgaus->SetParameter(1,1);
 fgaus->SetParameter(2,ProfileEEmCorrected->GetRMS());
 fgaus->SetRange(ProfileEEmCorrected->GetMean()-5*ProfileEEmCorrected->GetRMS(),ProfileEEmCorrected->GetMean()+5*ProfileEEmCorrected->GetRMS());
 fgaus->SetLineColor(kRed);
 ProfileEEmCorrected->Fit("fgaus","QRME");

 cout<<" EEm Corrected: Mean = "<<fgaus->GetParameter(1)<<" Sigma "<<fgaus->GetParameter(2)<<endl;

 
 fgaus->SetParameter(1,1);
 fgaus->SetParameter(2,ProfileEEpCorrected->GetRMS());
 fgaus->SetRange(ProfileEEpCorrected->GetMean()-5*ProfileEEpCorrected->GetRMS(),ProfileEEpCorrected->GetMean()+5*ProfileEEpCorrected->GetRMS());
 fgaus->SetLineColor(kBlue);
 ProfileEEpCorrected->Fit("fgaus","QRME");

 cout<<" EEp Corrected: Mean = "<<fgaus->GetParameter(1)<<" Sigma "<<fgaus->GetParameter(2)<<endl;

 /// reEvaluate precision of IC

 TH1F ***hspreadCorrected= new TH1F**[2];
 for(int i =0; i<2 ; i++) hspreadCorrected[i]= new TH1F*[50];

 TH1F** hspreadAllCorrected= new TH1F*[40];

 BookSpreadHistos_EE(hcmapFinalEE,hspreadCorrected,hspreadAllCorrected,eRings);

  /// Graph Error for spread EE+ and EE-

  TGraphErrors *sigma_vs_ringCorrected[3];
  sigma_vs_ringCorrected[0] = new TGraphErrors();
  sigma_vs_ringCorrected[0]->SetMarkerStyle(20);
  sigma_vs_ringCorrected[0]->SetMarkerSize(1);
  sigma_vs_ringCorrected[0]->SetMarkerColor(kBlue+2);
  
  sigma_vs_ringCorrected[1] = new TGraphErrors();
  sigma_vs_ringCorrected[1]->SetMarkerStyle(20);
  sigma_vs_ringCorrected[1]->SetMarkerSize(1);
  sigma_vs_ringCorrected[1]->SetMarkerColor(kBlue+2);

  sigma_vs_ringCorrected[2] = new TGraphErrors();
  sigma_vs_ringCorrected[2]->SetMarkerStyle(20);
  sigma_vs_ringCorrected[2]->SetMarkerSize(1);
  sigma_vs_ringCorrected[2]->SetMarkerColor(kBlue+2);
    
  int np2[3] = {0};

  /// Gaussian fit for EE+ and EE- ---> not folded
  for (int k = 0; k < 2 ; k++){
    for (int iring = 0; iring < 40; iring++){
      if (hspreadCorrected[k][iring]-> GetEntries() == 0) continue;
      float e     = 0.5*mapMomentumCorrected[k]-> GetYaxis()->GetBinWidth(1);
      fgaus->SetParameter(1,1);
      fgaus->SetParameter(2,hspreadCorrected[k][iring]->GetRMS());
      fgaus->SetRange(1-5*hspreadCorrected[k][iring]->GetRMS(),1+5*hspreadCorrected[k][iring]->GetRMS());
      hspreadCorrected[k][iring]->Fit("fgaus","QRN");
      sigma_vs_ringCorrected[k]-> SetPoint(np2[k],iring,fgaus->GetParameter(2));
      sigma_vs_ringCorrected[k]-> SetPointError(np2[k], e ,fgaus->GetParError(2));
      np2[k]++;    
    }
  }

  /// Folded spread distribution
  for (int iring = 0; iring < 40; iring++){
    if (hspreadAllCorrected[iring]-> GetEntries() == 0) continue;
    float e     = 0.5*mapMomentumCorrected[0]-> GetYaxis()->GetBinWidth(1);
      fgaus->SetParameter(1,1);
      fgaus->SetParameter(2,hspreadAllCorrected[iring]->GetRMS());
      fgaus->SetRange(1-5*hspreadAllCorrected[iring]->GetRMS(),1+5*hspreadAllCorrected[iring]->GetRMS());
      hspreadAllCorrected[iring]->Fit("fgaus","QRN");
      sigma_vs_ringCorrected[2]-> SetPoint(np2[2],iring,fgaus->GetParameter(2));
      sigma_vs_ringCorrected[2]-> SetPointError(np2[2], e ,fgaus->GetParError(2));
      np2[2]++;    
    }

  /// New residuals :

  TGraphErrors *residual_vs_ringCorrected[3];
  residual_vs_ringCorrected[0] = new TGraphErrors();
  residual_vs_ringCorrected[0]->SetMarkerStyle(20);
  residual_vs_ringCorrected[0]->SetMarkerSize(1);
  residual_vs_ringCorrected[0]->SetMarkerColor(kGreen+2);
  
  residual_vs_ringCorrected[1] = new TGraphErrors();
  residual_vs_ringCorrected[1]->SetMarkerStyle(20);
  residual_vs_ringCorrected[1]->SetMarkerSize(1);
  residual_vs_ringCorrected[1]->SetMarkerColor(kGreen+2);
  
  residual_vs_ringCorrected[2] = new TGraphErrors();
  residual_vs_ringCorrected[2]->SetMarkerStyle(20);
  residual_vs_ringCorrected[2]->SetMarkerSize(1);
  residual_vs_ringCorrected[2]->SetMarkerColor(kGreen+2);

  if(evalStat){  for (int k = 0; k < 3 ; k++) ResidualSpread (statprecision_vs_ring[k], sigma_vs_ringCorrected[k], residual_vs_ringCorrected[k]);
  }
 
  
  ///-----------------------------------------------------------------
  ///--- Draw plots
  ///-----------------------------------------------------------------
  TFile* outFile = new TFile(fileName.c_str(),"RECREATE");
  outFile -> cd();
  
  TCanvas *cEEP[15];
  TCanvas *cEEM[15];
  TCanvas *cAll[15];
  
  
  /// --- plot 0 : map of coefficients 
  //cEEM[0] = new TCanvas("cmapEEM","cmapEEM");
  //cEEM[0] -> cd();
  //cEEM[0]->SetLeftMargin(0.1); 
  //cEEM[0]->SetRightMargin(0.13); 
  //cEEM[0]->SetGridx();
  //cEEM[0]->SetGridy();
  
  hcmap[0]->GetXaxis() -> SetLabelSize(0.03);
  hcmap[0]->GetXaxis() ->SetTitle("ix");
  hcmap[0]->GetYaxis() ->SetTitle("iy");
  hcmap[0]->GetZaxis() ->SetRangeUser(0.8,1.2);
  //hcmap[0]->Draw("colz");
  hcmap[0]->Write();
  
  
  //cEEP[0] = new TCanvas("cmapEEP","cmapEEP");
  //cEEP[0] -> cd();
  //cEEP[0]->SetLeftMargin(0.1); 
  //cEEP[0]->SetRightMargin(0.13); 
  //cEEP[0]->SetGridx();
  //cEEP[0]->SetGridy();
  
  //hcmap[1]->GetXaxis()->SetNdivisions(1020);
  hcmap[1]->GetXaxis() -> SetLabelSize(0.03);
  hcmap[1]->GetXaxis() ->SetTitle("ix");
  hcmap[1]->GetYaxis() ->SetTitle("iy");
  hcmap[1]->GetZaxis() ->SetRangeUser(0.8,1.2);
  //hcmap[1]->Draw("colz");
  hcmap[1]->Write();
  
  
  /// --- plot 1 : ring precision vs ieta
  //cEEP[1] = new TCanvas("csigmaEEP","csigmaEEP");
  //cEEP[1]->SetGridx();
  //cEEP[1]->SetGridy();
  
  sigma_vs_ring[1]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
  sigma_vs_ring[1]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
  sigma_vs_ring[1]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
  sigma_vs_ring[1]->GetHistogram()->GetXaxis()-> SetTitle("ring");
  //sigma_vs_ring[1]->Draw("ap");
  if (evalStat){
    //statprecision_vs_ring[1]->Draw("psame");
    //sigma_vs_ring[1]->Draw("psame");
    statprecision_vs_ring[1]->Write("statprecision_vs_ring_EEp");
    sigma_vs_ring[1]->Write("sigma_vs_ring_EEp");
    TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
    leg->SetFillColor(0);
    leg->AddEntry(statprecision_vs_ring[1],"statistical precision", "LP");
    leg->AddEntry(sigma_vs_ring[1],"spread", "LP");
    //leg->Draw("same");
  }
  
  
  //cEEM[1] = new TCanvas("csigmaEEM","csigmaEEM");
  //cEEM[1]->SetGridx();
  //cEEM[1]->SetGridy();
  
  sigma_vs_ring[0]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
  sigma_vs_ring[0]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
  sigma_vs_ring[0]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
  sigma_vs_ring[0]->GetHistogram()->GetXaxis()-> SetTitle("ring");
  //sigma_vs_ring[0]->Draw("ap");
  if (evalStat){
    //statprecision_vs_ring[0]->Draw("psame");
    //sigma_vs_ring[0]->Draw("psame");
    statprecision_vs_ring[0]->Write("statprecision_vs_ring_EEm");
    sigma_vs_ring[0]->Write("sigma_vs_ring_EEm");
    TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
    leg->SetFillColor(0);
    leg->AddEntry(statprecision_vs_ring[0],"statistical precision", "LP");
    leg->AddEntry(sigma_vs_ring[0],"spread", "LP");
    //leg->Draw("same");
  }
  
  
  /// --- plot 2 : statistical precision vs ieta
  if (evalStat){
    //cEEP[2] = new TCanvas("cresidualP","cresidualP");
    //cEEP[2]->SetGridx();
    //cEEP[2]->SetGridy();
    residual_vs_ring[1]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    residual_vs_ring[1]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    residual_vs_ring[1]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
    residual_vs_ring[1]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    //residual_vs_ring[1]->Draw("ap");
    residual_vs_ring[1]->Write("residual_vs_ring_EEp");
    
    //cEEM[2] = new TCanvas("cresidualM","cresidualM");
    //cEEM[2]->SetGridx();
    //cEEM[2]->SetGridy();
    residual_vs_ring[0]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    residual_vs_ring[0]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    residual_vs_ring[0]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
    residual_vs_ring[0]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    //residual_vs_ring[0]->Draw("ap");
    residual_vs_ring[0]->Write("residual_vs_ring_EEm");
    
     
    //cAll[0] = new TCanvas("csigmaFolded","csigmaFolded");
    //cAll[0]->SetGridx();
    //cAll[0]->SetGridy();
    
    sigma_vs_ring[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
    sigma_vs_ring[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    sigma_vs_ring[2]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
    sigma_vs_ring[2]->GetHistogram()->GetXaxis()-> SetTitle("ring");
    //sigma_vs_ring[2]->Draw("ap");
    if (evalStat){
      //statprecision_vs_ring[2]->Draw("psame");
      //sigma_vs_ring[2]->Draw("psame");
      statprecision_vs_ring[2]->Write("statprecision_vs_ring");
      sigma_vs_ring[2]->Write("sigma_vs_ring");
      TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
      leg->SetFillColor(0);
      leg->AddEntry(statprecision_vs_ring[2],"statistical precision", "LP");
      leg->AddEntry(sigma_vs_ring[2],"spread", "LP");
      //leg->Draw("same");
    }
    
    
    //cAll[1] = new TCanvas("cresidualFolded","cresidualFolded");
    //cAll[1]->SetGridx();
    //cAll[1]->SetGridy();
    
    residual_vs_ring[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    residual_vs_ring[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    residual_vs_ring[2]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
    residual_vs_ring[2]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    //residual_vs_ring[2]->Draw("ap");
    residual_vs_ring[2]->Write("residual_vs_ring");
    
    
    /// save precision for MC comparison
    if(isMC==true)
    {
     TFile * output = new TFile ("output/StatPrec_MC_R9_EE.root","RECREATE");
     output->cd();
     statprecision_vs_ring[0]->SetName("gr_stat_prec_EEP");
     statprecision_vs_ring[1]->SetName("gr_stat_prec_EEM");
     statprecision_vs_ring[2]->SetName("gr_stat_prec");

     statprecision_vs_ring[0]->Write();
     statprecision_vs_ring[1]->Write();
     statprecision_vs_ring[2]->Write();
     output->Close();
    }
   }
   
  outFile -> cd();
   /// Orignal projection plot
   //cEEM[3] = new TCanvas("PhiProjectionEEM","PhiProjectionEEM");
   //cEEM[3]->SetGridx();
   //cEEM[3]->SetGridy();
   
   PhiProjectionEEm->GetHistogram()->GetYaxis()-> SetRangeUser(0.93,1.07);
   PhiProjectionEEm->GetHistogram()->GetXaxis()-> SetRangeUser(0,200);
   PhiProjectionEEm->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEm->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   //PhiProjectionEEm->Draw("apl");
   PhiProjectionEEm->Write("PhiProjectionEEm");
   
   
   //cEEP[3] = new TCanvas("PhiProjectionEEP","PhiProjectionEEP");
   //cEEP[3]->SetGridx();
   //cEEP[3]->SetGridy();
   
   PhiProjectionEEp->GetHistogram()->GetYaxis()-> SetRangeUser(0.93,1.07);
   PhiProjectionEEp->GetHistogram()->GetXaxis()-> SetRangeUser(0,200);
   PhiProjectionEEp->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEp->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   //PhiProjectionEEp->Draw("apl");
   PhiProjectionEEp->Write("PhiProjectionEEp");
   
   
   //cAll[2]= new TCanvas("PhiProjectionSame","PhiProjectionSame");
   //cAll[2]->SetGridx();
   //cAll[2]->SetGridy();
   
   //PhiProjectionEEp->Draw("apl");
   //PhiProjectionEEm->Draw("plsame");
   TLegend * leg2 = new TLegend(0.6,0.7,0.89, 0.89);
   leg2->SetFillColor(0);
   leg2->AddEntry(PhiProjectionEEm,"EE- Projection", "LP");
   leg2->AddEntry(PhiProjectionEEp,"EE+ Projection", "LP");
   //leg2->Draw("same");
   
   
   /// Map corrected
   //cEEM[4] = new TCanvas("cmapCorrectedEEM","cmapCorrectedEEM");
   //cEEM[4] -> cd();
   //cEEM[4]->SetLeftMargin(0.1); 
   //cEEM[4]->SetRightMargin(0.13); 
   //cEEM[4]->SetGridx();
   //cEEM[4]->SetGridy();
   hcmapFinalEE[0]->GetXaxis() -> SetLabelSize(0.03);
   hcmapFinalEE[0]->GetXaxis() ->SetTitle("ix");
   hcmapFinalEE[0]->GetYaxis() ->SetTitle("iy");
   hcmapFinalEE[0]->GetZaxis() ->SetRangeUser(0.8,1.2);
   //hcmapFinalEE[0]->Draw("colz");
   hcmapFinalEE[0]->Write();
   
   
   //cEEP[4] = new TCanvas("cmapCorrectedEEP","cmapCorrectedEEP");
   //cEEP[4] -> cd();
   //cEEP[4]->SetLeftMargin(0.1); 
   //cEEP[4]->SetRightMargin(0.13); 
   //cEEP[4]->SetGridx();
   //cEEP[4]->SetGridy();
   hcmapFinalEE[1]->GetXaxis() -> SetLabelSize(0.03);
   hcmapFinalEE[1]->GetXaxis() ->SetTitle("ix");
   hcmapFinalEE[1]->GetYaxis() ->SetTitle("iy");
   hcmapFinalEE[1]->GetZaxis() ->SetRangeUser(0.8,1.2);
   //hcmapFinalEE[1]->Draw("colz");
   hcmapFinalEE[1]->Write();
   
   
   /// Projection after correction
   //cEEM[5] = new TCanvas("PhiProjectionEEM_Corrected","PhiProjectionEEM_Corrected");
   //cEEM[5]->SetGridx();
   //cEEM[5]->SetGridy();
   PhiProjectionEEm_Corrected->GetHistogram()->GetYaxis()-> SetRangeUser(0.93,1.07);
   PhiProjectionEEm_Corrected->GetHistogram()->GetXaxis()-> SetRangeUser(100,200);
   PhiProjectionEEm_Corrected->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEm_Corrected->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   //PhiProjectionEEm_Corrected->Draw("apl");
   PhiProjectionEEm_Corrected->Write("PhiProjectionEEm_Corrected");
   
   
   //cEEP[5] = new TCanvas("PhiProjectionEEP_Corrected","PhiProjectionEEP_Corrected");
   //cEEP[5]->SetGridx();
   //cEEP[5]->SetGridy();
   PhiProjectionEEp_Corrected->GetHistogram()->GetYaxis()-> SetRangeUser(0.93,1.07);
   PhiProjectionEEp_Corrected->GetHistogram()->GetXaxis()-> SetRangeUser(0,200);
   PhiProjectionEEp_Corrected->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEp_Corrected->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   //PhiProjectionEEp_Corrected->Draw("apl");
   PhiProjectionEEp_Corrected->Write("PhiProjectionEEp_Corrected");
   
   
   //cAll[3] = new TCanvas("PhiProjectionSame_Corrected","PhiProjectionSame_Corrected");
   //cAll[3]->SetGridx();
   //cAll[3]->SetGridy();
   //PhiProjectionEEp_Corrected->Draw("apl");
   //PhiProjectionEEm_Corrected->Draw("plsame");
   TLegend * leg3 = new TLegend(0.6,0.7,0.89, 0.89);
   leg3->SetFillColor(0);
   leg3->AddEntry(PhiProjectionEEm,"EE- Projection Corrected", "LP");
   leg3->AddEntry(PhiProjectionEEp,"EE+ Projection Corrected", "LP");
   //leg3->Draw("same");
   
   
   //cAll[4] = new TCanvas("ProfileEEm","ProfileEEm");
   //cAll[4]->SetGridx();
   //cAll[4]->SetGridy();
   ProfileEEm->GetXaxis()->SetTitle("#bar{IC}");
   ProfileEEm->SetLineColor(kBlack);
   ProfileEEm->SetMarkerSize(0.8);
   ProfileEEmCorrected->SetLineColor(kRed);
   ProfileEEmCorrected->SetMarkerSize(0.8);
   ProfileEEmCorrected->SetLineWidth(2);
   //ProfileEEmCorrected->Draw();
   //ProfileEEm->Draw("same");
   ProfileEEmCorrected->Write();
   ProfileEEm->Write();
   
   TLegend * leg4 = new TLegend(0.6,0.7,0.89, 0.89);
   leg4->SetFillColor(0);
   leg4->AddEntry(ProfileEEm,"EE- Projection ", "LP");
   leg4->AddEntry(ProfileEEmCorrected,"EE- Projection Corrected ", "LP");
   //leg4->Draw("same");
   
   
   //cAll[5] = new TCanvas("ProfileEEp","ProfileEEp");
   //cAll[5]->SetGridx();
   //cAll[5]->SetGridy();
   ProfileEEp->GetXaxis()->SetTitle("#bar{IC}");
   ProfileEEp->SetLineColor(kBlack);
   ProfileEEp->SetMarkerSize(0.8);
   ProfileEEpCorrected->SetLineColor(kBlue);
   ProfileEEpCorrected->SetMarkerSize(0.8);
   ProfileEEpCorrected->SetLineWidth(2);
   //ProfileEEpCorrected->Draw();
   //ProfileEEp->Draw("same");
   ProfileEEpCorrected->Write();
   ProfileEEp->Write();
   
   TLegend * leg5 = new TLegend(0.6,0.7,0.89, 0.89);
   leg5->SetFillColor(0);
   leg5->AddEntry(ProfileEEp,"EE+ Projection ", "LP");
   leg5->AddEntry(ProfileEEpCorrected,"EE+ Projection Corrected ", "LP");
   //leg5->Draw("same");
   
   
   /// new Precision and residual:
   //cAll[6] = new TCanvas("csigmaFoldedCorrected","csigmaFoldedCorrected");
   //cAll[6]->SetGridx();
   //cAll[6]->SetGridy();
   sigma_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
   sigma_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
   sigma_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
   sigma_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetTitle("ring");
   //sigma_vs_ringCorrected[2]->Draw("ap");
   if (evalStat){
     //statprecision_vs_ring[2]->Draw("psame");
     //sigma_vs_ringCorrected[2]->Draw("psame");
     statprecision_vs_ring[2]->Write();
     sigma_vs_ringCorrected[2]->Write();
     TLegend * leg6 = new TLegend(0.6,0.7,0.89, 0.89);
     leg6->SetFillColor(0);
     leg6->AddEntry(statprecision_vs_ring[2],"statistical precision", "LP");
     leg6->AddEntry(sigma_vs_ringCorrected[2],"spread corrected", "LP");
     //leg6->Draw("same");
   }
   
   
   //cAll[7] = new TCanvas("cresidualFoldedCorrected","cresidualFoldedCorrected");
   //cAll[7]->SetGridx();
   //cAll[7]->SetGridy();
   residual_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
   residual_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
   residual_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
   residual_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
   //residual_vs_ringCorrected[2]->Draw("ap");
   residual_vs_ringCorrected[2]->Write("residual_vs_ringCorrected");
 
   /// IC constant set:

   TFile *exisistingEE = new TFile ("output/existingEE.root","READ");
   TH2F* exmap = (TH2F*) exisistingEE->Get("endcap");

   TH2F* warning_Map_EEP = (TH2F*) exmap->Clone("warning_Map_EEP");
   warning_Map_EEP->Reset("ICMES");
   warning_Map_EEP->Reset();

   TH2F* warning_Map_EEM = (TH2F*) exmap->Clone("warning_Map_EEM");
   warning_Map_EEM->Reset("ICMES");
   warning_Map_EEM->Reset();

   /// Acquisition residual systematic from scalib MC
  
   TFile *f4= new TFile(SystematicToAdd.c_str(),"READ");
   TGraphErrors* systematicEEM = (TGraphErrors*) f4->Get("residual_EEM");
   TGraphErrors* systematicEEP = (TGraphErrors*) f4->Get("residual_EEP");
   TGraphErrors* systematicAll = (TGraphErrors*) f4->Get("residual_All");


   TGraphErrors *tempEEP = new TGraphErrors();
   TGraphErrors *tempEEM = new TGraphErrors();
   TGraphErrors *tempAll = new TGraphErrors();


   tempEEP->SetMarkerStyle(20);
   tempEEP->SetMarkerSize(1);
   tempEEP->SetMarkerColor(kBlue);

   tempEEM->SetMarkerStyle(20);
   tempEEM->SetMarkerSize(1);
   tempEEM->SetMarkerColor(kBlue);
  
   tempAll->SetMarkerStyle(20);
   tempAll->SetMarkerSize(1);
   tempAll->SetMarkerColor(kBlue);
   
  
   for(int i=0 ; i<statprecision_vs_ring[0]->GetN(); i++){
    double k,h,j,m;
    statprecision_vs_ring[0]->GetPoint(i,k,h);
    systematicEEM->GetPoint(i,j,m);
    if(m!=0){ tempEEM->SetPoint(i,k,sqrt(h*h+m*m));
             tempEEM->SetPointError(i,sqrt(statprecision_vs_ring[0]->GetErrorX(i)*statprecision_vs_ring[0]->GetErrorX(i)+
                     systematicEEM->GetErrorX(i)*systematicEEM->GetErrorX(i)),sqrt(statprecision_vs_ring[0]->GetErrorY(i)*statprecision_vs_ring[0]->GetErrorY(i)+systematicEEM->GetErrorY(i)*systematicEEM->GetErrorY(i)));}
               
    else{ tempEEM->SetPoint(i,k,h);
          tempEEM->SetPointError(i,statprecision_vs_ring[0]->GetErrorX(i),statprecision_vs_ring[0]->GetErrorY(i));
         }
   }
  
    for(int i=0 ; i<statprecision_vs_ring[1]->GetN(); i++){
    double k,h,j,m;
    statprecision_vs_ring[1]->GetPoint(i,k,h);
    systematicEEP->GetPoint(i,j,m);
    if(m!=0){ tempEEP->SetPoint(i,k,sqrt(h*h+m*m));
             tempEEP->SetPointError(i,sqrt(statprecision_vs_ring[1]->GetErrorX(i)*statprecision_vs_ring[1]->GetErrorX(i)+
                     systematicEEP->GetErrorX(i)*systematicEEP->GetErrorX(i)),sqrt(statprecision_vs_ring[1]->GetErrorY(i)*statprecision_vs_ring[1]->GetErrorY(i)+systematicEEP->GetErrorY(i)*systematicEEP->GetErrorY(i)));}

    else{ tempEEP->SetPoint(i,k,h);
          tempEEP->SetPointError(i,statprecision_vs_ring[1]->GetErrorX(i),statprecision_vs_ring[1]->GetErrorY(i));
         }

   }

    for(int i=0 ; i<statprecision_vs_ring[2]->GetN(); i++){
    double k,h,j,m;
    statprecision_vs_ring[2]->GetPoint(i,k,h);
    systematicAll->GetPoint(i,j,m);
    if(m!=0){ tempAll->SetPoint(i,k,sqrt(h*h+m*m));
             tempAll->SetPointError(i,sqrt(statprecision_vs_ring[2]->GetErrorX(i)*statprecision_vs_ring[2]->GetErrorX(i)+
                     systematicAll->GetErrorX(i)*systematicAll->GetErrorX(i)),sqrt(statprecision_vs_ring[2]->GetErrorY(i)*statprecision_vs_ring[2]->GetErrorY(i)+systematicAll->GetErrorY(i)*systematicAll->GetErrorY(i)));}

    else{ tempAll->SetPoint(i,k,h);
          tempAll->SetPointError(i,statprecision_vs_ring[2]->GetErrorX(i),statprecision_vs_ring[2]->GetErrorY(i));}
   
   }
   
   
    outFile -> cd();
    
    //cEEM[7] = new TCanvas("UncertaintyEEM","UncertaintyEEM");
    //cEEM[7]->SetGridx();
    //cEEM[7]->SetGridy();
    
    tempEEM->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    tempEEM->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    tempEEM->GetHistogram()->GetYaxis()-> SetTitle("Total Error");
    tempEEM->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    //tempEEM->Draw("ap");
    tempEEM->Write("tempEEM");
    //statprecision_vs_ring[0]->Draw("psame");
    TLegend * leg7 = new TLegend(0.6,0.7,0.89, 0.89);
    leg7->SetFillColor(0);
    leg7->AddEntry(statprecision_vs_ring[0],"statistical precision EEM", "LP");
    leg7->AddEntry(tempEEM,"Statistical + Systematic MC ", "LP");
    //leg7->Draw("same");
    
    
    //cEEP[7] = new TCanvas("UncertaintyEEP","UncertaintyEEP");
    //cEEP[7]->SetGridx();
    //cEEP[7]->SetGridy();
    
    tempEEP->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    tempEEP->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    tempEEP->GetHistogram()->GetYaxis()-> SetTitle("Total Error");
    tempEEP->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    //tempEEP->Draw("ap");
    tempEEP->Write("tempEEP");
    //statprecision_vs_ring[1]->Draw("psame");
    TLegend * leg8 = new TLegend(0.6,0.7,0.89, 0.89);
    leg8->SetFillColor(0);
    leg8->AddEntry(statprecision_vs_ring[1],"statistical precision EEP", "LP");
    leg8->AddEntry(tempEEP,"Statistical + Systematic MC ", "LP");
    //leg8->Draw("same");
    
        
    //cAll[8] = new TCanvas("UncertaintyAll","UncertaintyAll");
    //cAll[8]->SetGridx();
    //cAll[8]->SetGridy();
    
    tempAll->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    tempAll->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    tempAll->GetHistogram()->GetYaxis()-> SetTitle("Total Error");
    tempAll->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    //tempAll->Draw("ap");
    tempAll->Write("tempAll");
    //statprecision_vs_ring[2]->Draw("psame");
    TLegend * leg9 = new TLegend(0.6,0.7,0.89, 0.89);
    leg9->SetFillColor(0);
    leg9->AddEntry(statprecision_vs_ring[2],"statistical precision folded", "LP");
    leg9->AddEntry(tempEEP,"Statistical + Systematic MC ", "LP");
    //leg9->Draw("same");
  

  if(isMC == false)
  {
   std::ofstream outTxt (outputTxt.c_str(),std::ios::out);

   outTxt << "---------------------------------------------------------------" << std::endl;
   outTxt << std::fixed << std::setprecision(0) << std::setw(10) << "iX"
          << std::fixed << std::setprecision(0) << std::setw(10) << "iY"
          << std::fixed << std::setprecision(0) << std::setw(10) << "iZ"
          << std::fixed << std::setprecision(6) << std::setw(15) << "IC"
          << std::fixed << std::setprecision(6) << std::setw(15) << "error"
          << std::endl;
   outTxt << "---------------------------------------------------------------" << std::endl;

    for (int ix = 1; ix < hcmapFinalEE[0]->GetNbinsX()+1 ; ix ++)
    {
      for (int iy = 1; iy < hcmapFinalEE[0] -> GetNbinsY()+1; iy++)
	{
          if( exmap->GetBinContent(ix,iy) !=1) continue;

	  double X,statPrec,Y,sysPrec;
	  statprecision_vs_ring[0]->GetPoint(int(eRings->GetEndcapRing(ix,iy,-1)),X,statPrec);
          systematicAll->GetPoint(int(eRings->GetEndcapRing(ix,iy,-1)),Y,sysPrec);

          if( (hcmapFinalEE[0]->GetBinContent(ix,iy)>0.4 && hcmapFinalEE[0]->GetBinContent(ix,iy)<2.)|| hcmapFinalEE[0]->GetBinContent(ix,iy)!=0. )
	    
	    outTxt << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[0]->GetXaxis()->GetBinLowEdge(ix)
		   << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[0]->GetYaxis()->GetBinLowEdge(iy)
		   << std::fixed << std::setprecision(0) << std::setw(10) << "-1"
		   << std::fixed << std::setprecision(6) << std::setw(15) << hcmapFinalEE[0]->GetBinContent(ix,iy)
		   << std::fixed << std::setprecision(6) << std::setw(15) << sqrt(statPrec*statPrec+sysPrec*sysPrec)
		   << std::endl;

          else{

            outTxt << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[0]->GetXaxis()->GetBinLowEdge(ix)
                   << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[0]->GetYaxis()->GetBinLowEdge(iy)
                   << std::fixed << std::setprecision(0) << std::setw(10) << "-1"
                   << std::fixed << std::setprecision(6) << std::setw(15) << "-1."
                   << std::fixed << std::setprecision(6) << std::setw(15) << "999."
                   << std::endl;

	    warning_Map_EEM->Fill(ix,iy);
	  }
	  
	}
    }
    
    for (int ix = 1; ix < hcmapFinalEE[1]->GetNbinsX()+1 ; ix ++)
      {
	for (int iy = 1; iy < hcmapFinalEE[1] -> GetNbinsY()+1; iy++)
	  {
	    if( exmap->GetBinContent(ix,iy) !=1) continue;

	  double X,statPrec,Y,sysPrec;
	  statprecision_vs_ring[1]->GetPoint(int(eRings->GetEndcapRing(ix,iy,1)),X,statPrec);
          systematicAll->GetPoint(int(eRings->GetEndcapRing(ix,iy,1)),Y,sysPrec);

	if((hcmapFinalEE[1]->GetBinContent(ix,iy)>0.4 && hcmapFinalEE[1]->GetBinContent(ix,iy)<2.)|| hcmapFinalEE[1]->GetBinContent(ix,iy)!=0.)

	      outTxt << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[1]->GetXaxis()->GetBinLowEdge(ix)
		     << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[1]->GetYaxis()->GetBinLowEdge(iy)
		     << std::fixed << std::setprecision(0) << std::setw(10) << "1"
		     << std::fixed << std::setprecision(6) << std::setw(15) << hcmapFinalEE[1]->GetBinContent(ix,iy)
		     << std::fixed << std::setprecision(6) << std::setw(15) << sqrt(statPrec*statPrec+sysPrec*sysPrec)
		     << std::endl;

	    else{

              outTxt << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[1]->GetXaxis()->GetBinLowEdge(ix)
                     << std::fixed << std::setprecision(0) << std::setw(10) << hcmapFinalEE[1]->GetYaxis()->GetBinLowEdge(iy)
                     << std::fixed << std::setprecision(0) << std::setw(10) << "1"
                     << std::fixed << std::setprecision(6) << std::setw(15) << "-1."
                     << std::fixed << std::setprecision(6) << std::setw(15) << "999."
                     << std::endl;


	      warning_Map_EEP->Fill(ix,iy);
	    }
	  }
      }
  }
  
  
  //cEEP[8] = new TCanvas("Warning_EEP","Warning_EEP");
  //cEEP[8]->SetGridx();
  //cEEP[8]->SetGridy();
  warning_Map_EEP->GetXaxis() -> SetLabelSize(0.03);
  //warning_Map_EEP->Draw("COLZ");
  warning_Map_EEP->GetXaxis() ->SetTitle("ix");
  warning_Map_EEP->GetYaxis() ->SetTitle("iy");
  warning_Map_EEP->GetZaxis() ->SetRangeUser(0.85,1.15);

  //cEEM[8] = new TCanvas("Warning_EEM","Warning_EEM");
  //cEEM[8]->SetGridx();
  //cEEM[8]->SetGridy();
  warning_Map_EEM->GetXaxis() -> SetLabelSize(0.03);
  //warning_Map_EEM->Draw("COLZ");
  warning_Map_EEM->GetXaxis() ->SetTitle("ix");
  warning_Map_EEM->GetYaxis() ->SetTitle("iy");
  warning_Map_EEM->GetZaxis() ->SetRangeUser(0.85,1.15);
  
  outFile -> Close();
  
  theApp->Run();
  
  return 0;
}


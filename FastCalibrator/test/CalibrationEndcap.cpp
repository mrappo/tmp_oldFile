#include <vector>
#include <utility>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "../CommonTools/TEndcapRings.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "TApplication.h"

/// Stand-alone program to produce ECAL single electron calibration plots for EE
/// Input Files : MC or Data splistat and no splitstat, Momentum scale vs phi plot
/// Output Files :  StatPrec_MC_R9_EE.root --> stat precision MC usefull for CompareCalibMCTruth_EE.C only MC
//                  
using namespace std;

/// Run ./bin/CalibrationEndacp.cpp cfg/calibrationEE_cfg.cfg

int main (int argc, char **argv)
{

  /// by xtal

  int nbins = 200;

  /// Set style options
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(0); 
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
  std::string infile2 = gConfigParser -> readStringOption("Input::Inputfile2");
  std::string infile3 = gConfigParser -> readStringOption("Input::Inputfile3");
  std::string inputMomentumScale =  gConfigParser -> readStringOption("Input::inputMomentumScale");
  int evalStat = gConfigParser -> readIntOption("Input::evalStat");
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
 
  std::string fileType = gConfigParser -> readStringOption("Output::fileType");
  std::string dirName = gConfigParser -> readStringOption("Output::dirName");
  bool printPlots = gConfigParser -> readBoolOption("Output::printPlots");

  TApplication* theApp = new TApplication("Application",&argc, argv);
  
  /// imput file with full statistic normlized to the mean in a ring

  TFile *f = new TFile(infile1.c_str());
  TH2F *hcmap[2];
  hcmap[0] = (TH2F*)f->Get("h_scale_map_EEM");
  hcmap[1] = (TH2F*)f->Get("h_scale_map_EEP");
    

  /// ring geometry for the endcap
  TH2F *hrings[2];
  hrings[0] = (TH2F*)hcmap[0]->Clone("hringsEEM");
  hrings[1] = (TH2F*)hcmap[0]->Clone("hringsEEP");
  hrings[0] ->Reset("ICMES");
  hrings[1] ->Reset("ICMES");
  hrings[0] ->ResetStats();
  hrings[1] ->ResetStats();

  FILE *fRing;
  fRing = fopen("macros/eerings.dat","r");
  int x,y,z,ir;
  while(fscanf(fRing,"(%d,%d,%d) %d \n",&x,&y,&z,&ir) !=EOF ) {
    if(z>0) hrings[1]->Fill(x,y,ir); 
    if(z<0) hrings[0]->Fill(x,y,ir);
  }

  ///--------------------------------------------------------------------------------
  ///--- Build the precision vs ring plot starting from the TH2F of IC folded and not
  ///--------------------------------------------------------------------------------

  char hname[100];
  TH1F *hspread[2][50];
  TH1F* hspreadAll [40];


  for (int k = 0; k < 2; k++){         
    for (int iring = 0; iring < 40 ; iring++){
      if (k==0){
       sprintf(hname,"hspreadAll_ring%02d",iring);
       hspreadAll[iring] = new TH1F(hname, hname, nbins,0.,2.);
       sprintf(hname,"hspreadEEM_ring%02d",iring);
       hspread[k][iring] = new TH1F(hname, hname, nbins,0.,2.);
      }
      else{ sprintf(hname,"hspreadEEP_ring%02d",iring);
            hspread[k][iring] = new TH1F(hname, hname, nbins,0.,2.); }
   }
  }
  
 /// spread all distribution, spread for EE+ and EE- and comparison with the MC truth
  for (int k = 0; k < 2 ; k++){
    for (int ix = 1; ix < 101; ix++){
      for (int iy = 1; iy < 101; iy++){
	int iz = k;
	if (k==0) iz = -1;
	int mybin = hcmap[k] -> FindBin(ix,iy);
	int ring  = int(hrings[1]-> GetBinContent(mybin));
	float ic = hcmap[k]->GetBinContent(mybin);
 	if (ic>0){
	  hspread[k][ring]->Fill(ic);
          hspreadAll[ring]->Fill(ic);}
      }
    }
  }
  
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
      hspread[k][iring]->Fit("fgaus","QR");
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
      hspreadAll[iring]->Fit("fgaus","QR");
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
   TH2F *hcmap2[2];
   hcmap2[0] = (TH2F*)f2->Get("h_scale_map_EEM"); 
   hcmap2[1] = (TH2F*)f2->Get("h_scale_map_EEP");

   TH2F *hcmap3[2];
   hcmap3[0] = (TH2F*)f3->Get("h_scale_map_EEM"); 
   hcmap3[1] = (TH2F*)f3->Get("h_scale_map_EEP");

   TH1F *hstatprecision[2][40];
   TH1F *hstatprecisionAll[40];

    /// stat precision histos for each EE ring
   for (int k = 0; k < 2; k++){
     for (int iring = 0; iring < 40 ; iring ++){
       if (k==0)
	 { sprintf(hname,"hstatprecisionAll_ring%02d",iring);
           hstatprecisionAll[iring] = new TH1F(hname, hname, nbins,-1.3,1.3);
           sprintf(hname,"hstatprecisionEEM_ring%02d",iring);
           hstatprecision[k][iring] = new TH1F(hname, hname, nbins,-1.3,1.3);}
	else {sprintf(hname,"hstatprecisionEEP_ring%02d",iring);
              hstatprecision[k][iring] = new TH1F(hname, hname, nbins,-1.3,1.3);}
       }
     }
//     
    for (int k = 0; k < 2 ; k++){
      for (int ix = 1; ix < 102; ix++){
	for (int iy = 1; iy < 102; iy++){
	  int iz = k;
	  if (k==0) iz = -1;
	  int mybin = hcmap2[k] -> FindBin(ix,iy);
	  int ring  = int(hrings[1]-> GetBinContent(mybin));
	  float ic1 = hcmap2[k]->GetBinContent(mybin);
	  float ic2 = hcmap3[k]->GetBinContent(mybin);
          if (ic1>0 && ic2 >0){
	    hstatprecision[k][ring]->Fill((ic1-ic2)/(ic1+ic2)); /// sigma (diff/sum) gives the stat. precision on teh entire sample
            hstatprecisionAll[ring]->Fill((ic1-ic2)/(ic1+ic2));
	  }
	}
      }
    }
 
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

        hstatprecision[k][iring]->Fit("fgaus","QR");
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
        hstatprecisionAll[iring]->Fit("fgaus","QR");
      
	statprecision_vs_ring[2]-> SetPoint(n[2],iring,fgaus->GetParameter(2));
	statprecision_vs_ring[2]-> SetPointError(n[2],e,fgaus->GetParError(2));
	n[2]++;
      }
     
    
    TH1F *hresidual[3];
    hresidual[0] = new TH1F("hresidualEEM","hresidualEEM",1000,0,1);
    hresidual[1] = new TH1F("hresidualEEP","hresidualEEP",1000,0,1);
    hresidual[2] = new TH1F("hresidualAll","hresidualAll",1000,0,1);

    TH1F *hstat[3];
    hstat[0] = new TH1F("hstatEEM","hstatEEM",1000,0,0.5);
    hstat[1] = new TH1F("hstatEEP","hstatEEP",1000,0,0.5);
    hstat[2] = new TH1F("hstatAll","hstatAll",1000,0,0.5);
   
    TH1F *hspre[3];
    hspre[0] = new TH1F("hspreEEM","hspreEEM",1000,0,0.5);
    hspre[1] = new TH1F("hspreEEP","hspreEEP",1000,0,0.5);
    hspre[2] = new TH1F("hspreAll","hspreAll",1000,0,0.5);

    /// Residual spread plot

    for (int k = 0; k < 3 ; k++){
      for (int i= 0; i < statprecision_vs_ring[k]-> GetN(); i++){
	double spread, espread;
	double stat, estat;
	double residual, eresidual;
	double xdummy,ex;
	sigma_vs_ring[k]-> GetPoint(i, xdummy, spread );
	espread = sigma_vs_ring[k]-> GetErrorY(i);
	statprecision_vs_ring[k]-> GetPoint(i, xdummy, stat );
	estat = statprecision_vs_ring[k]-> GetErrorY(i);
	ex = statprecision_vs_ring[k]-> GetErrorX(i);
	if (spread > stat ){
	  residual  = sqrt( spread*spread - stat*stat );
	  eresidual = sqrt( pow(spread*espread,2) + pow(stat*estat,2))/residual;
	}
	else {
	  residual = 0;
	  eresidual = 0;
	}
	residual_vs_ring[k]->SetPoint(i,xdummy, residual);
	residual_vs_ring[k]->SetPointError(i,ex,eresidual);
      }
    }
 
 }

 /// Momentum scale correction
 
 TFile* input = new TFile(inputMomentumScale.c_str());
 
 TGraphErrors* g_EoC_EE = (TGraphErrors*) input->Get("g_EoC_EE_0");
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


 TEndcapRings *eRings = new TEndcapRings();
 std::vector<double> vectSum;
 std::vector<double> vectCounter;
 
 vectCounter.assign(g_EoC_EE->GetN(),0.);
 vectSum.assign(g_EoC_EE->GetN(),0.);

 /// EE+ and EE- projection 

 for(int ix=1; ix<hcmap[0]->GetNbinsX()+1;ix++){
   for(int iy=1; iy<hcmap[0]->GetNbinsY()+1;iy++){
    if(hcmap[0]->GetBinContent(ix,iy)==0) continue;
      int iPhi = int(eRings->GetEndcapIphi(ix,iy,-1)/(360./g_EoC_EE->GetN()));
      vectSum.at(iPhi)=vectSum.at(iPhi)+hcmap[0]->GetBinContent(ix,iy);
      vectCounter.at(iPhi)=vectCounter.at(iPhi)+1;
  }
 }


 for(unsigned int i=0; i<vectCounter.size();i++)
  PhiProjectionEEm->SetPoint(i,int(i*(360./g_EoC_EE->GetN())),vectSum.at(i)/vectCounter.at(i));

 for(unsigned int i=0; i<vectSum.size(); i++){
  vectSum.at(i)=0; vectCounter.at(i)=0;
 }

 for(int ix=1; ix<hcmap[1]->GetNbinsX()+1;ix++){
   for(int iy=1; iy<hcmap[1]->GetNbinsY()+1;iy++){
    if(hcmap[1]->GetBinContent(ix,iy)==0) continue;
     int iPhi = int(eRings->GetEndcapIphi(ix,iy,1)/(360./g_EoC_EE->GetN()));
     vectSum.at(iPhi)=vectSum.at(iPhi)+hcmap[1]->GetBinContent(ix,iy);
     vectCounter.at(iPhi)=vectCounter.at(iPhi)+1;
  }
 }

 
 for(unsigned int i=0; i<vectCounter.size();i++)
  PhiProjectionEEp->SetPoint(i,int(i*(360./g_EoC_EE->GetN())),vectSum.at(i)/vectCounter.at(i));
 
 for(unsigned int i=0; i<vectSum.size(); i++){
  vectSum.at(i)=0; vectCounter.at(i)=0;
 }

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
    double xphi,yphi;
    g_EoC_EE->GetPoint(int(iPhi/(360./PhiProjectionEEm->GetN())),xphi,yphi);
    mapMomentumCorrected[0]->SetBinContent(ix,iy,hcmap[0]->GetBinContent(ix,iy)*yphi);
  }
 }

 for(int ix=1; ix<hcmap[1]->GetNbinsX()+1;ix++){
  for(int iy=1; iy<hcmap[1]->GetNbinsY()+1;iy++){
   if(hcmap[1]->GetBinContent(ix,iy)==0) continue;
    int iPhi = int(eRings->GetEndcapIphi(ix,iy,1));
    double xphi,yphi; 
    g_EoC_EE->GetPoint(int(iPhi/(360./PhiProjectionEEp->GetN())),xphi,yphi);
    mapMomentumCorrected[1]->SetBinContent(ix,iy,hcmap[1]->GetBinContent(ix,iy)*yphi);
  }
 }
  
 
 /// EE+ and EE- projection after correction

 for(int ix=1; ix<mapMomentumCorrected[0]->GetNbinsX()+1; ix++){
   for(int iy=1; iy<mapMomentumCorrected[0]->GetNbinsY()+1; iy++){
     if(mapMomentumCorrected[0]->GetBinContent(ix,iy)==0) continue;
     int iPhi = int(eRings->GetEndcapIphi(ix,iy,-1)/(360./g_EoC_EE->GetN()));
     vectSum.at(iPhi)=vectSum.at(iPhi)+mapMomentumCorrected[0]->GetBinContent(ix,iy);
     vectCounter.at(iPhi)=vectCounter.at(iPhi)+1;
  }
 }


 for(unsigned int i=0; i<vectCounter.size();i++)
  PhiProjectionEEm_Corrected->SetPoint(i,int(i*(360./g_EoC_EE->GetN())),vectSum.at(i)/vectCounter.at(i));
 
 for(unsigned int i=0; i<vectSum.size(); i++){
  vectSum.at(i)=0; vectCounter.at(i)=0;
 }

 for(int ix=1; ix<mapMomentumCorrected[1]->GetNbinsX()+1;ix++){
   for(int iy=1; iy<mapMomentumCorrected[1]->GetNbinsY()+1;iy++){
    if(mapMomentumCorrected[1]->GetBinContent(ix,iy)==0) continue;
     int iPhi = int(eRings->GetEndcapIphi(ix,iy,1)/(360./g_EoC_EE->GetN()));
     vectSum.at(iPhi)=vectSum.at(iPhi)+mapMomentumCorrected[1]->GetBinContent(ix,iy);
     vectCounter.at(iPhi)=vectCounter.at(iPhi)+1;
  }
 }
 
 for(unsigned int i=0; i<vectCounter.size();i++)
  PhiProjectionEEp_Corrected->SetPoint(i,int(i*(360./g_EoC_EE->GetN())),vectSum.at(i)/vectCounter.at(i));
 
 for(unsigned int i=0; i<vectSum.size(); i++){
  vectSum.at(i)=0; vectCounter.at(i)=0;
 }
    
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
 ProfileEEm->Fit("fgaus","QR");

 cout<<" EEm Uncorrected: Mean = "<<fgaus->GetParameter(1)<<" Sigma "<<fgaus->GetParameter(2)<<endl;

 fgaus->SetParameter(1,1);
 fgaus->SetParameter(2,ProfileEEp->GetRMS());
 fgaus->SetRange(ProfileEEp->GetMean()-5*ProfileEEp->GetRMS(),ProfileEEp->GetMean()+5*ProfileEEm->GetRMS());
 ProfileEEp->SetLineWidth(kBlack);
 ProfileEEp->SetLineWidth(2);
 ProfileEEp->Fit("fgaus","QR");

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
 ProfileEEmCorrected->Fit("fgaus","QR");

 cout<<" EEm Corrected: Mean = "<<fgaus->GetParameter(1)<<" Sigma "<<fgaus->GetParameter(2)<<endl;

 
 fgaus->SetParameter(1,1);
 fgaus->SetParameter(2,ProfileEEpCorrected->GetRMS());
 fgaus->SetRange(ProfileEEpCorrected->GetMean()-5*ProfileEEpCorrected->GetRMS(),ProfileEEpCorrected->GetMean()+5*ProfileEEpCorrected->GetRMS());
 fgaus->SetLineColor(kBlue);
 ProfileEEpCorrected->Fit("fgaus","QR");

 cout<<" EEp Corrected: Mean = "<<fgaus->GetParameter(1)<<" Sigma "<<fgaus->GetParameter(2)<<endl;

 /// reEvaluate precision of IC

 TH1F *hspreadCorrected[2][50];
 TH1F* hspreadAllCorrected [40];


  for (int k = 0; k < 2; k++){         
    for (int iring = 0; iring < 40 ; iring++){
      if (k==0){
       sprintf(hname,"hspreadAllCorrected_ring%02d",iring);
       hspreadAllCorrected[iring] = new TH1F(hname, hname, nbins,0.,2.);
       sprintf(hname,"hspreadEEMCorrected_ring%02d",iring);
       hspreadCorrected[k][iring] = new TH1F(hname, hname, nbins,0.,2.);
      }
      else{ sprintf(hname,"hspreadEEPCorrected_ring%02d",iring);
            hspreadCorrected[k][iring] = new TH1F(hname, hname, nbins,0.,2.); }
   }
  }
  
 /// spread all distribution, spread for EE+ and EE- 
  for (int k = 0; k < 2 ; k++){
    for (int ix = 1; ix < 101; ix++){
      for (int iy = 1; iy < 101; iy++){
	int iz = k;
	if (k==0) iz = -1;
	int mybin = mapMomentumCorrected[k] -> FindBin(ix,iy);
	int ring  = int(hrings[1]-> GetBinContent(mybin));
	float ic = mapMomentumCorrected[k]->GetBinContent(mybin);
 	if (ic>0){
	  hspreadCorrected[k][ring]->Fill(ic);
          hspreadAllCorrected[ring]->Fill(ic);}
      }
    }
  }
  
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
      hspreadCorrected[k][iring]->Fit("fgaus","QR");
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
      hspreadAllCorrected[iring]->Fit("fgaus","QR");
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

  if(evalStat)
  {
     for (int k = 0; k < 3 ; k++){
      for (int i= 0; i < statprecision_vs_ring[k]-> GetN(); i++){
	double spread, espread;
	double stat, estat;
	double residual, eresidual;
	double xdummy,ex;
	sigma_vs_ringCorrected[k]-> GetPoint(i, xdummy, spread );
	espread = sigma_vs_ringCorrected[k]-> GetErrorY(i);
	statprecision_vs_ring[k]-> GetPoint(i, xdummy, stat );
	estat = statprecision_vs_ring[k]-> GetErrorY(i);
	ex = statprecision_vs_ring[k]-> GetErrorX(i);
	if (spread > stat ){
	  residual  = sqrt( spread*spread - stat*stat );
	  eresidual = sqrt( pow(spread*espread,2) + pow(stat*estat,2))/residual;
	}
	else {
	  residual = 0;
	  eresidual = 0;
	}
	residual_vs_ringCorrected[k]->SetPoint(i,xdummy, residual);
	residual_vs_ringCorrected[k]->SetPointError(i,ex,eresidual);
      }
    }
  }

  
  ///-----------------------------------------------------------------
  ///--- Draw plots
  ///-----------------------------------------------------------------
  TCanvas *cEEP[15];
  TCanvas *cEEM[15];
  TCanvas *cAll[15];

  /// --- plot 0 : map of coefficients 
  cEEM[0] = new TCanvas("cmapEEM","cmapEEM");
  cEEM[0] -> cd();
  cEEM[0]->SetLeftMargin(0.1); 
  cEEM[0]->SetRightMargin(0.13); 
  cEEM[0]->SetGridx();
  cEEM[0]->SetGridy();
  hcmap[0]->GetXaxis() -> SetLabelSize(0.03);
  hcmap[0]->GetXaxis() ->SetTitle("ix");
  hcmap[0]->GetYaxis() ->SetTitle("iy");
  hcmap[0]->GetZaxis() ->SetRangeUser(0.8,1.2);
  hcmap[0]->Draw("colz");
 
  cEEP[0] = new TCanvas("cmapEEP","cmapEEP");
  cEEP[0] -> cd();
  cEEP[0]->SetLeftMargin(0.1); 
  cEEP[0]->SetRightMargin(0.13); 
  cEEP[0]->SetGridx();
  cEEP[0]->SetGridy();
  //  hcmap[1]->GetXaxis()->SetNdivisions(1020);
  hcmap[1]->GetXaxis() -> SetLabelSize(0.03);
  hcmap[1]->GetXaxis() ->SetTitle("ix");
  hcmap[1]->GetYaxis() ->SetTitle("iy");
  hcmap[1]->GetZaxis() ->SetRangeUser(0.8,1.2);
  hcmap[1]->Draw("colz");
 
 
  /// --- plot 1 : ring precision vs ieta
  cEEP[1] = new TCanvas("csigmaEEP","csigmaEEP");
  cEEP[1]->SetGridx();
  cEEP[1]->SetGridy();
  sigma_vs_ring[1]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
  sigma_vs_ring[1]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
  sigma_vs_ring[1]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
  sigma_vs_ring[1]->GetHistogram()->GetXaxis()-> SetTitle("ring");
  sigma_vs_ring[1]->Draw("ap");
  if (evalStat){
    statprecision_vs_ring[1]->Draw("psame");
    sigma_vs_ring[1]->Draw("psame");
    TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
    leg->SetFillColor(0);
    leg->AddEntry(statprecision_vs_ring[1],"statistical precision", "LP");
    leg->AddEntry(sigma_vs_ring[1],"spread", "LP");
    leg->Draw("same");
  }

  cEEM[1] = new TCanvas("csigmaEEM","csigmaEEM");
  cEEM[1]->SetGridx();
  cEEM[1]->SetGridy();
  sigma_vs_ring[0]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
  sigma_vs_ring[0]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
  sigma_vs_ring[0]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
  sigma_vs_ring[0]->GetHistogram()->GetXaxis()-> SetTitle("ring");
  sigma_vs_ring[0]->Draw("ap");
  if (evalStat){
    statprecision_vs_ring[0]->Draw("psame");
    sigma_vs_ring[0]->Draw("psame");
    TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
    leg->SetFillColor(0);
    leg->AddEntry(statprecision_vs_ring[0],"statistical precision", "LP");
    leg->AddEntry(sigma_vs_ring[0],"spread", "LP");
    leg->Draw("same");
  }


  /// --- plot 2 : statistical precision vs ieta

  if (evalStat){
    cEEP[2] = new TCanvas("cstatP","cstatP");
    cEEP[2]->SetGridx();
    cEEP[2]->SetGridy();
    statprecision_vs_ring[1]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    statprecision_vs_ring[1]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    statprecision_vs_ring[1]->GetHistogram()->GetYaxis()-> SetTitle("#sigma((c_{P}-c_{D})/(c_{P}+c_{D}))");
    statprecision_vs_ring[1]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    statprecision_vs_ring[1]->Draw("ap");
    
    cEEP[3] = new TCanvas("cresidualP","cresidualP");
    cEEP[3]->SetGridx();
    cEEP[3]->SetGridy();
    residual_vs_ring[1]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    residual_vs_ring[1]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    residual_vs_ring[1]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
    residual_vs_ring[1]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    residual_vs_ring[1]->Draw("ap");
 
    cEEM[2] = new TCanvas("cstatM","cstatM");
    cEEM[2]->SetGridx();
    cEEM[2]->SetGridy();
    statprecision_vs_ring[0]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    statprecision_vs_ring[0]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    statprecision_vs_ring[0]->GetHistogram()->GetYaxis()-> SetTitle("#sigma((c_{P}-c_{D})/(c_{P}+c_{D}))");
    statprecision_vs_ring[0]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    statprecision_vs_ring[0]->Draw("ap");
    
    cEEM[3] = new TCanvas("cresidualM","cresidualM");
    cEEM[3]->SetGridx();
    cEEM[3]->SetGridy();
    residual_vs_ring[0]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    residual_vs_ring[0]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    residual_vs_ring[0]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
    residual_vs_ring[0]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    residual_vs_ring[0]->Draw("ap");
 
    cAll[0] = new TCanvas("csigmaFolded","csigmaFolded");
    cAll[0]->SetGridx();
    cAll[0]->SetGridy();
    sigma_vs_ring[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
    sigma_vs_ring[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    sigma_vs_ring[2]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
    sigma_vs_ring[2]->GetHistogram()->GetXaxis()-> SetTitle("ring");
    sigma_vs_ring[2]->Draw("ap");
    if (evalStat){
    statprecision_vs_ring[2]->Draw("psame");
    sigma_vs_ring[2]->Draw("psame");
    TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
    leg->SetFillColor(0);
    leg->AddEntry(statprecision_vs_ring[2],"statistical precision", "LP");
    leg->AddEntry(sigma_vs_ring[2],"spread", "LP");
    leg->Draw("same");
    }

    cAll[1] = new TCanvas("cresidualFolded","cresidualFolded");
    cAll[1]->SetGridx();
    cAll[1]->SetGridy();
    residual_vs_ring[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    residual_vs_ring[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    residual_vs_ring[2]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
    residual_vs_ring[2]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    residual_vs_ring[2]->Draw("ap");
 
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
    }
   }
  
   /// Orignal projection plot
   cEEM[4] = new TCanvas("PhiProjectionEEM","PhiProjectionEEM");
   cEEM[4]->SetGridx();
   cEEM[4]->SetGridy();
   PhiProjectionEEm->GetHistogram()->GetYaxis()-> SetRangeUser(0.90,1.1);
   PhiProjectionEEm->GetHistogram()->GetXaxis()-> SetRangeUser(0,360);
   PhiProjectionEEm->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEm->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   PhiProjectionEEm->Draw("apl");
  
   cEEP[4] = new TCanvas("PhiProjectionEEP","PhiProjectionEEP");
   cEEP[4]->SetGridx();
   cEEP[4]->SetGridy();
   PhiProjectionEEp->GetHistogram()->GetYaxis()-> SetRangeUser(0.90,1.1);
   PhiProjectionEEp->GetHistogram()->GetXaxis()-> SetRangeUser(0,360);
   PhiProjectionEEp->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEp->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   PhiProjectionEEp->Draw("apl");

   cAll[2]= new TCanvas("PhiProjectionSame","PhiProjectionSame");
   cAll[2]->SetGridx();
   cAll[2]->SetGridy();
   PhiProjectionEEp->Draw("apl");
   PhiProjectionEEm->Draw("plsame");
   TLegend * leg2 = new TLegend(0.6,0.7,0.89, 0.89);
   leg2->SetFillColor(0);
   leg2->AddEntry(PhiProjectionEEm,"EE- Projection", "LP");
   leg2->AddEntry(PhiProjectionEEp,"EE+ Projection", "LP");
   leg2->Draw("same");
   

   /// Map corrected
   cEEM[5] = new TCanvas("cmapCorrectedEEM","cmapCorrectedEEM");
   cEEM[5] -> cd();
   cEEM[5]->SetLeftMargin(0.1); 
   cEEM[5]->SetRightMargin(0.13); 
   cEEM[5]->SetGridx();
   cEEM[5]->SetGridy();
   mapMomentumCorrected[0]->GetXaxis() -> SetLabelSize(0.03);
   mapMomentumCorrected[0]->GetXaxis() ->SetTitle("ix");
   mapMomentumCorrected[0]->GetYaxis() ->SetTitle("iy");
   mapMomentumCorrected[0]->GetZaxis() ->SetRangeUser(0.8,1.2);
   mapMomentumCorrected[0]->Draw("colz");

   cEEP[5] = new TCanvas("cmapCorrectedEEP","cmapCorrectedEEP");
   cEEP[5] -> cd();
   cEEP[5]->SetLeftMargin(0.1); 
   cEEP[5]->SetRightMargin(0.13); 
   cEEP[5]->SetGridx();
   cEEP[5]->SetGridy();
   mapMomentumCorrected[1]->GetXaxis() -> SetLabelSize(0.03);
   mapMomentumCorrected[1]->GetXaxis() ->SetTitle("ix");
   mapMomentumCorrected[1]->GetYaxis() ->SetTitle("iy");
   mapMomentumCorrected[1]->GetZaxis() ->SetRangeUser(0.8,1.2);
   mapMomentumCorrected[1]->Draw("colz");

   /// Projection after correction
   cEEM[6] = new TCanvas("PhiProjectionEEM_Corrected","PhiProjectionEEM_Corrected");
   cEEM[6]->SetGridx();
   cEEM[6]->SetGridy();
   PhiProjectionEEm_Corrected->GetHistogram()->GetYaxis()-> SetRangeUser(0.9,1.1);
   PhiProjectionEEm_Corrected->GetHistogram()->GetXaxis()-> SetRangeUser(0,360);
   PhiProjectionEEm_Corrected->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEm_Corrected->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   PhiProjectionEEm_Corrected->Draw("apl");
  
   cEEP[6] = new TCanvas("PhiProjectionEEP_Corrected","PhiProjectionEEP_Corrected");
   cEEP[6]->SetGridx();
   cEEP[6]->SetGridy();
   PhiProjectionEEp_Corrected->GetHistogram()->GetYaxis()-> SetRangeUser(0.9,1.1);
   PhiProjectionEEp_Corrected->GetHistogram()->GetXaxis()-> SetRangeUser(0,360);
   PhiProjectionEEp_Corrected->GetHistogram()->GetYaxis()-> SetTitle("#bar{IC}");
   PhiProjectionEEp_Corrected->GetHistogram()->GetXaxis()-> SetTitle("#phi");
   PhiProjectionEEp_Corrected->Draw("apl");

   cAll[3] = new TCanvas("PhiProjectionSame_Corrected","PhiProjectionSame_Corrected");
   cAll[3]->SetGridx();
   cAll[3]->SetGridy();
   PhiProjectionEEp_Corrected->Draw("apl");
   PhiProjectionEEm_Corrected->Draw("plsame");
   TLegend * leg3 = new TLegend(0.6,0.7,0.89, 0.89);
   leg3->SetFillColor(0);
   leg3->AddEntry(PhiProjectionEEm,"EE- Projection Corrected", "LP");
   leg3->AddEntry(PhiProjectionEEp,"EE+ Projection Corrected", "LP");
   leg3->Draw("same");

 
   cAll[4] = new TCanvas("ProfileEEm","ProfileEEm");
   cAll[4]->SetGridx();
   cAll[4]->SetGridy();
   ProfileEEm->GetXaxis()->SetTitle("#bar{IC}");
   ProfileEEm->SetLineColor(kBlack);
   ProfileEEm->SetMarkerSize(0.8);
   ProfileEEmCorrected->SetLineColor(kRed);
   ProfileEEmCorrected->SetMarkerSize(0.8);
   ProfileEEmCorrected->SetLineWidth(2);
   ProfileEEmCorrected->Draw();
   ProfileEEm->Draw("same");
   
   TLegend * leg4 = new TLegend(0.6,0.7,0.89, 0.89);
   leg4->SetFillColor(0);
   leg4->AddEntry(ProfileEEm,"EE- Projection ", "LP");
   leg4->AddEntry(ProfileEEmCorrected,"EE- Projection Corrected ", "LP");
   leg4->Draw("same");
  
   cAll[5] = new TCanvas("ProfileEEp","ProfileEEp");
   cAll[5]->SetGridx();
   cAll[5]->SetGridy();
   ProfileEEp->GetXaxis()->SetTitle("#bar{IC}");
   ProfileEEp->SetLineColor(kBlack);
   ProfileEEp->SetMarkerSize(0.8);
   ProfileEEpCorrected->SetLineColor(kBlue);
   ProfileEEpCorrected->SetMarkerSize(0.8);
   ProfileEEpCorrected->SetLineWidth(2);
   ProfileEEpCorrected->Draw();
   ProfileEEp->Draw("same");
  
   TLegend * leg5 = new TLegend(0.6,0.7,0.89, 0.89);
   leg5->SetFillColor(0);
   leg5->AddEntry(ProfileEEp,"EE+ Projection ", "LP");
   leg5->AddEntry(ProfileEEpCorrected,"EE+ Projection Corrected ", "LP");
   leg5->Draw("same");

   /// new Precision and residual:

    cAll[6] = new TCanvas("csigmaFoldedCorrected","csigmaFoldedCorrected");
    cAll[6]->SetGridx();
    cAll[6]->SetGridy();
    sigma_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.20);
    sigma_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    sigma_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
    sigma_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetTitle("ring");
    sigma_vs_ringCorrected[2]->Draw("ap");
    if (evalStat){
    statprecision_vs_ring[2]->Draw("psame");
    sigma_vs_ringCorrected[2]->Draw("psame");
    TLegend * leg6 = new TLegend(0.6,0.7,0.89, 0.89);
    leg6->SetFillColor(0);
    leg6->AddEntry(statprecision_vs_ring[2],"statistical precision", "LP");
    leg6->AddEntry(sigma_vs_ringCorrected[2],"spread corrected", "LP");
    leg6->Draw("same");
    }

    cAll[7] = new TCanvas("cresidualFoldedCorrected","cresidualFoldedCorrected");
    cAll[7]->SetGridx();
    cAll[7]->SetGridy();
    residual_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    residual_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    residual_vs_ringCorrected[2]->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
    residual_vs_ringCorrected[2]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    residual_vs_ringCorrected[2]->Draw("ap");
 
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
  
   TFile *f4= new TFile("output/ResidualForSystematic_EE_MC.root","READ");
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
    if(m>0) tempEEM->SetPoint(i,k,sqrt(h*h+m*m));
    else tempEEM->SetPoint(i,k,h);
   }
  
    for(int i=0 ; i<statprecision_vs_ring[1]->GetN(); i++){
    double k,h,j,m;
    statprecision_vs_ring[1]->GetPoint(i,k,h);
    systematicEEP->GetPoint(i,j,m);
    if(m>0) tempEEP->SetPoint(i,k,sqrt(h*h+m*m));
    else tempEEP->SetPoint(i,k,h);
   }

    for(int i=0 ; i<statprecision_vs_ring[2]->GetN(); i++){
    double k,h,j,m;
    statprecision_vs_ring[2]->GetPoint(i,k,h);
    systematicAll->GetPoint(i,j,m);
    if(m>0) tempAll->SetPoint(i,k,sqrt(h*h+m*m));
    else tempAll->SetPoint(i,k,h);
   }

    cAll[8] = new TCanvas("UncertaintyAll","UncertaintyAll");
    cAll[8]->SetGridx();
    cAll[8]->SetGridy();
    tempAll->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.10);
    tempAll->GetHistogram()->GetXaxis()-> SetRangeUser(0,40);
    tempAll->GetHistogram()->GetYaxis()-> SetTitle("Total Error");
    tempAll->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
    tempAll->Draw("ap");
    statprecision_vs_ring[2]->Draw("psame");


  if(isMC == false)
  {
   std::ofstream outTxt ("Calibration_Coefficient_EE.txt",std::ios::out);

   outTxt << "---------------------------------------------------------------" << std::endl;
   outTxt << std::fixed << std::setprecision(0) << std::setw(10) << "iX"
          << std::fixed << std::setprecision(0) << std::setw(10) << "iY"
          << std::fixed << std::setprecision(0) << std::setw(10) << "iZ"
          << std::fixed << std::setprecision(6) << std::setw(15) << "IC"
          << std::fixed << std::setprecision(6) << std::setw(15) << "error"
          << std::endl;
   outTxt << "---------------------------------------------------------------" << std::endl;

    for (int ix = 1; ix < mapMomentumCorrected[0]->GetNbinsX()+1 ; ix ++)
    {
      for (int iy = 1; iy < mapMomentumCorrected[0] -> GetNbinsY()+1; iy++)
	{
          if( exmap->GetBinContent(ix,iy) !=1) continue;

	  double X,statPrec,Y,sysPrec;
	  statprecision_vs_ring[0]->GetPoint(int(hrings[0]->GetBinContent(ix,iy)),X,statPrec);
          systematicAll->GetPoint(int(hrings[0]->GetBinContent(ix,iy)),Y,sysPrec);

          if( (mapMomentumCorrected[0]->GetBinContent(ix,iy)>0.4 && mapMomentumCorrected[0]->GetBinContent(ix,iy)<2.)|| mapMomentumCorrected[0]->GetBinContent(ix,iy)==0 )
	    
	    outTxt << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[0]->GetXaxis()->GetBinLowEdge(ix)
		   << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[0]->GetYaxis()->GetBinLowEdge(iy)
		   << std::fixed << std::setprecision(0) << std::setw(10) << "-1"
		   << std::fixed << std::setprecision(6) << std::setw(15) << mapMomentumCorrected[0]->GetBinContent(ix,iy)
		   << std::fixed << std::setprecision(6) << std::setw(15) << sqrt(statPrec*statPrec+sysPrec*sysPrec)
		   << std::endl;

          else{

            outTxt << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[0]->GetXaxis()->GetBinLowEdge(ix)
                   << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[0]->GetYaxis()->GetBinLowEdge(iy)
                   << std::fixed << std::setprecision(0) << std::setw(10) << "-1"
                   << std::fixed << std::setprecision(6) << std::setw(15) << "0"
                   << std::fixed << std::setprecision(6) << std::setw(15) << sqrt(statPrec*statPrec+sysPrec*sysPrec)
                   << std::endl;

	    warning_Map_EEM->Fill(ix,iy);
	  }
	  
	}
    }
    
    for (int ix = 1; ix < mapMomentumCorrected[1]->GetNbinsX()+1 ; ix ++)
      {
	for (int iy = 1; iy < mapMomentumCorrected[1] -> GetNbinsY()+1; iy++)
	  {
	    if( exmap->GetBinContent(ix,iy) !=1) continue;

	  double X,statPrec,Y,sysPrec;
	  statprecision_vs_ring[1]->GetPoint(int(hrings[1]->GetBinContent(ix,iy)),X,statPrec);
          systematicAll->GetPoint(int(hrings[1]->GetBinContent(ix,iy)),Y,sysPrec);

	    if((mapMomentumCorrected[1]->GetBinContent(ix,iy)>0.4 && mapMomentumCorrected[1]->GetBinContent(ix,iy)<2.)|| mapMomentumCorrected[1]->GetBinContent(ix,iy)==0)

	      outTxt << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[1]->GetXaxis()->GetBinLowEdge(ix)
		     << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[1]->GetYaxis()->GetBinLowEdge(iy)
		     << std::fixed << std::setprecision(0) << std::setw(10) << "1"
		     << std::fixed << std::setprecision(6) << std::setw(15) << mapMomentumCorrected[1]->GetBinContent(ix,iy)
		     << std::fixed << std::setprecision(6) << std::setw(15) << sqrt(statPrec*statPrec+sysPrec*sysPrec)
		     << std::endl;

	    else{

              outTxt << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[1]->GetXaxis()->GetBinLowEdge(ix)
                     << std::fixed << std::setprecision(0) << std::setw(10) << mapMomentumCorrected[1]->GetYaxis()->GetBinLowEdge(iy)
                     << std::fixed << std::setprecision(0) << std::setw(10) << "1"
                     << std::fixed << std::setprecision(6) << std::setw(15) << "0"
                     << std::fixed << std::setprecision(6) << std::setw(15) << sqrt(statPrec*statPrec+sysPrec*sysPrec)
                     << std::endl;


	      warning_Map_EEP->Fill(ix,iy);
	    }
	  }
      }
  }
  
  
  cEEP[8] = new TCanvas("Warning_EEP","Warning_EEP");
  cEEP[8]->SetGridx();
  cEEP[8]->SetGridy();
  warning_Map_EEP->GetXaxis() -> SetLabelSize(0.03);
  warning_Map_EEP->Draw("COLZ");
  warning_Map_EEP->GetXaxis() ->SetTitle("ix");
  warning_Map_EEP->GetYaxis() ->SetTitle("iy");
  warning_Map_EEP->GetZaxis() ->SetRangeUser(0.85,1.15);

  cEEM[8] = new TCanvas("Warning_EEM","Warning_EEM");
  cEEM[8]->SetGridx();
  cEEM[8]->SetGridy();
  warning_Map_EEM->GetXaxis() -> SetLabelSize(0.03);
  warning_Map_EEM->Draw("COLZ");
  warning_Map_EEM->GetXaxis() ->SetTitle("ix");
  warning_Map_EEM->GetYaxis() ->SetTitle("iy");
  warning_Map_EEM->GetZaxis() ->SetRangeUser(0.85,1.15);

  theApp->Run();
  
  return 0;
}


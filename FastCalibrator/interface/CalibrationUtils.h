#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TEndcapRings.h"


//############## ECAL BARREL ####################

/// Check if the crystal is near to a dead one
bool CheckxtalIC_EB (TH2F* h_scale_EB,int iPhi, int iEta );

/// Check if the crystal is near to a dead TT
bool CheckxtalTT_EB (int iPhi, int iEta, const std::vector<std::pair<int,int> >& TT_centre );

/// Initialize TT dead map
void InitializeDeadTT_EB(std::vector<std::pair<int,int> >& TT_centre);

void InitializeDeadTT_EB2012(std::vector<std::pair<int,int> >& TT_centre);


/// Normalize IC EB
void NormalizeIC_EB(TH2F* h_scale_EB, TH2F* hcmap,const std::vector< std::pair<int,int> > & TT_centre, bool skip = true);

/// Book spread Histos
void BookSpreadHistos_EB(TH2F* hcmap, TH1F **hspreadEtaFold, const int & ringGroupSize,const int & nEtaRing);

/// Book spread stat Histos
void BookSpreadStatHistos_EB(TH2F* hcmap2,TH2F* hcmap3,TH1F **hstatprecisionEtaFold,const int & ringGroupSize,const int & nEtaRing);

/// Phi Projection EB
void PhiProfile(TGraphErrors *phiProjection, TGraphErrors **MomentumScale, TH2F* hcmap);


/// Residual Spread 
void ResidualSpread (TGraphErrors *statprecision, TGraphErrors *Spread, TGraphErrors *Residual);

//################# ECAL ENDCAPS #####################

/// check if the xtal is near to a dead one
bool CheckxtalIC_EE(TH2F* h_scale_EE,int ix, int iy, int ir);

/// check if the xtal is neat to a dead TT
bool CheckxtalTT_EE(int ix, int iy, int ir,const std::vector<std::pair<int,int> >& TT_centre );

/// Map dead TT EE+
void InitializeDeadTTEEP(std::vector<std::pair<int,int> >& TT_centre);

void InitializeDeadTTEEP2012(std::vector<std::pair<int,int> >& TT_centre);

/// Map dead TT EE-
void InitializeDeadTTEEM(std::vector<std::pair<int,int> >& TT_centre);

void InitializeDeadTTEEM2012(std::vector<std::pair<int,int> >& TT_centre);

/// Normalize in function of ring
void NormalizeIC_EE(TH2F** hcmap, TH2F** hcmap2, const std::vector< std::pair<int,int> > & TT_centre_EEP,const  std::vector< std::pair<int,int> > & TT_centre_EEM, TEndcapRings *eRings, bool skip = true);

/// Book spread  Histos
void BookSpreadHistos_EE(TH2F** hcmap, TH1F ***hspread, TH1F **hspreadAll,  TEndcapRings *eRings);

/// Book spread stat Histos
void BookSpreadStatHistos_EE(TH2F** hcmap2,TH2F** hcmap3, TH1F ***hstatprecision, TH1F **hstatprecisionAll,  TEndcapRings *eRings);

/// Phi Projection EB
void PhiProfileEE(TGraphErrors *phiProjection, TGraphErrors **MomentumScale, TH2F* hcmap,TEndcapRings *eRings, const int & iz);

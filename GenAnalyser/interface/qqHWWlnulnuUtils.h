/*
qqHWWlnulnuUtils
*/
#ifndef qqHWWlnulnuUtils_h
#define qqHWWlnulnuUtils_h

#include <vector>
#include <cmath>

#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include "treeReader.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

#include "Math/GenVector/VectorUtil.h"

#include <algorithm>


///==== GetTrendInfo ====
///==== Transform TH1 with "trace" information to TH1 ====
TH1F* GetTrendInfo(TH1F* hTrend, double min = -1.5, double max = 1.5);
 
///==== Pull Plot: drawing utility ====
void PullPlot(TCanvas* canvas, TH1* hDATA, TH1* hMC);
void PullPlot(TCanvas* canvas, TH1* hDATA, THStack* hsMC);  
TH1F* PullPlot(TH1F* hDATA, TH1F* hMC);
std::pair<TGraphErrors*, TGraphErrors*> grPullPlot(TH1F* hDATA, TH1F* hMC);


///==== Draw Stack ====
void DrawStack(THStack* hs, int error = 0, double syst = 0);
void DrawStackError(THStack* hs, double syst = 0);
void DrawStackNormalized(THStack* hs);

///==== Add systrematic error ====
void AddError(THStack* hs, double syst = 0);


///==== MC Decay Channel of V V ====
///==== lepton - lepton ==== < 1 ,   > 
///==== emu - emu ========== < 1 , 1 >
///==== emu - tau ========== < 1 , 2 >
///==== tau - emu ========== < 1 , 3 >
///==== tau - tau ========== < 1 , 4 >
///==== lepton - quark ===== < 2 ,   >
///==== emu - quark ======== < 2 , 1 >
///==== quark - emu ======== < 2 , 2 >
///==== tau - quark ======== < 2 , 3 >
///==== quark - tau ======== < 2 , 4 >
///==== quark - quark ====== < 3 , -1 >

std::pair<int,int> GetMCDecayChannel(const std::vector<float>& pdgId);
std::pair<int,int> GetMCDecayChannel(float pdgId11,float pdgId12, float pdgId21, float pdgId22);




///==== get GoodCombination ====
///==== return the value of the Max MVA and the corresponding number in the combinations ====
std::pair<double,int> GetCombinationMVA(treeReader& reader, std::vector<std::vector<int> >& combinations,TMVA::Reader* TMVAreader, TString& methodName, Float_t* input_variables, std::vector<int>* whitelistJet = NULL);


///==== return the value of the Min Chi2 and the corresponding number in the combinations ====
std::pair<double,int> GetCombinationChi2(treeReader& reader, std::vector<std::vector<int> >& combinations,std::vector<int>* whitelistJet);
 


///==== get Selected Lepton ====
int SelectLeptonTTBar(std::vector<ROOT::Math::XYZTVector>& leptons,const std::string& method,const double& ptMin,const double& etaMax,int& numLeptons,const std::vector<int>* blacklist);

///==== get Selected Muon ====
int SelectMuonTTBar(std::vector<ROOT::Math::XYZTVector>& leptons,
		    std::vector<float>& tkIso,
		    std::vector<float>& emIso,
		    std::vector<float>& hadIso,
		    const std::string& method,
		    const double& ptMin,
		    const double& etaMax,
		    const double& isoMax,
		    int& numLeptons,
		    const std::vector<int>* blacklist
		    );

///==== get Selected Electron ====
int SelectElectronTTBar(std::vector<ROOT::Math::XYZTVector>& leptons,
		    std::vector<float>& tkIso,
		    std::vector<float>& emIso,
		    std::vector<float>& hadIso,
		    std::vector<float>& eleId,
		    std::vector<float>& eled0,
		    const std::string& method,
		    const double& ptMin,
		    const double& etaMax,
		    const double& isoMax,
		    const double& d0Max,
		    int& numLeptons,
		    const std::vector<int>* blacklist
		    );

///==== get num Jets ====
int getNumJets(const std::vector<int>& whitelist);


///==== CopyTree ====
void CopyTree();



///==== get Jet Energy Correction ====
double getJEC(const ROOT::Math::XYZTVector& Jet);


///==== get best Combination Jet ID with MVA ====
std::pair<double,int> GetCombination_Jets_ID_MVA(std::vector<ROOT::Math::XYZTVector> &jets, std::vector<std::vector<int> >& combinations,TMVA::Reader* TMVAreader, TString& methodName, Float_t* input_variables, std::vector<int>* whitelistJet);
// std::pair<double,int> GetCombination_Jets_ID_MVA(treeReader& reader, std::vector<std::vector<int> >& combinations,TMVA::Reader* TMVAreader, TString& methodName, Float_t* input_variables, std::vector<int>* whitelistJet);
 


///==== read list of systematics file ====
void ReadFileSystematics(std::string CutSystematicFile, std::vector< std::pair< int, std::pair<std::string, double> > >& listSystematics);
void ReadFileSystematicsWithRegion(std::string CutSystematicFile, std::vector< std::pair< std::string, std::string> >& listSystematics);


///==== std::string replace ====
// std::string::size_type repl(std::string& s, const std::string& from, const std::string& to);
void repl(std::string& s, const std::string& from, const std::string& to);

///==== modify list of cuts to include systematics ====
void ModifyCut(std::vector <std::string> & vCut, const std::vector< std::pair< int, std::pair<std::string, double> > >& listSystematics);
void ModifyCutWithRegion(std::vector <std::string> & vCut, const std::vector< std::pair< std::string, std::string> >& listSystematics);
 


///==== CJV with JES estimation ====
int getCJVJES(std::vector<ROOT::Math::XYZTVector>& jets,
	   int q1,
	   int q2,
	   const double& EtMin,
	   const std::vector<int>* blacklist = 0,
	   int kind = 1);   //---- kind = +1 or -1
 


/** Electron isolation / ID */
bool IsEleIsolatedID_VBF( treeReader& reader,const std::vector<double>& BarrelSelections, const std::vector<double>& EndCapSelections, int iEle);
bool IsEleIsolatedIDPUCorrected_VBF( treeReader& reader,const std::vector<double>& BarrelSelections, const std::vector<double>& EndCapSelections, int iEle);
bool IsEle_VBFMerged( treeReader& reader, int iEle);

/** Muon isolation  / ID */
bool IsMuIsolatedID_VBF( treeReader& reader,const std::vector<double>& Selections, int iMu);
bool IsMuIsolatedIDPUCorrected_VBF( treeReader& reader,const std::vector<double>& Selections, int iMu);
bool IsMu_VBFMerged( treeReader& reader, int iMu);
bool IsMu_Soft( treeReader& reader, int iMu);

/** Jet ID */
bool IsJetID( treeReader& reader, int iJet);

#endif


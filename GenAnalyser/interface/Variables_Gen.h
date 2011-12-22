#include "treeReader.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "hFactory.h"
#include "h2Factory.h"
#include "stdHisto.h"

#include "TH1F.h"
#include "TProfile.h"
#include "TObject.h"
#include "TTree.h"

#include "qqHWWlnulnuUtils.h"
#include <vector>
#include "Math/GenVector/VectorUtil.h"

struct Variables_Gen
{
 ///____ tree definition
 TFile* m_outputRootFile;
 TTree* m_reducedTree;
 TTree* m_efficiencyTree;
 
 
 ///____ input parameters
 int runId; 
 int lumiId; 
 int eventId; 
 int numEntriesBefore;
 double preselection_efficiency;
 double XSection;
 double XSectionErrorUp;
 double XSectionErrorDown;
 
 ///____ lepton variables
 
 
 double l1_pX;
 double l1_pY;
 double l1_pZ;
 double l1_pT;
 double l1_E;
 double l1_Eta;
 double l1_Phi;
 
 double l1_charge;
 int    l1_flavour;
 
 
 double l2_pX;
 double l2_pY;
 double l2_pZ;
 double l2_pT;
 double l2_E;
 double l2_Eta;
 double l2_Phi;
 
 double l2_charge;
 int    l2_flavour;
 

 double M_ll;
 double DEta_ll;
 double DPhi_ll;
 
 double l1_Z;
 double l2_Z;
 double Z_ll;
 
 int Nleptons_pT5;
 int Nleptons_pT10;
 int Nleptons_pT15;
 int Nleptons_pT20;
 int Nleptons_pT25;
 int Nleptons_pT30;
 
 
 ///____ Met variables
 double met_X;
 double met_Y;
 double met;
 double pmet; 
 double chmet; 
 double chmet_X; 
 double chmet_Y; 
 double pchmet; 
 double minMet;

 ///____ jetdeau [jd] variables: met and jet corrections

 double jd_dphi1;
 double jd_dphi2;

 double jd_Delta_q1_P;
 double jd_Delta_q2_P;  

 double jd_between;
 
 double jd_met;
 double jd_met_X;
 double jd_met_Y;
 double jd_pmet;
 
 double jd_q1_pX;
 double jd_q1_pY;
 double jd_q1_pZ;
 double jd_q1_pT;
 double jd_q1_E;
 double jd_q1_Eta;
 double jd_q1_Phi;
 
 double jd_q2_pX;
 double jd_q2_pY;
 double jd_q2_pZ;
 double jd_q2_pT;
 double jd_q2_E;
 double jd_q2_Eta;
 double jd_q2_Phi;
 
 double jd_M_qq;
 

 /// _______ transverse mass
 double mT;

 ///____ dphi_Jet_ll
 
 double DPhiJet_ll;
 double maxDPhiJet_ll;
 
 double DPhiSingleJet_ll;
 double DPhiDoubleJet_ll;
  
 ///____ jet variables
 
 double q1_pX;
 double q1_pY;
 double q1_pZ;
 double q1_pT;
 double q1_E;
 double q1_Eta;
 double q1_Phi;
 
 double q2_pX;
 double q2_pY;
 double q2_pZ;
 double q2_pT;
 double q2_E;
 double q2_Eta;
 double q2_Phi;
 
 double q3_pX ;
 double q3_pY ;
 double q3_pZ ;
 double q3_pT ;
 double q3_E ;
 double q3_Eta;
 double q3_Phi;

 
 double M_qq;
 double DEta_qq;
 double DPhi_qq;
  
 int JV_20;
 int JV_30;
 int JV_40;
 int CJV_20;
 int CJV_30;
 int CJV_40;
 
 int q_Z_01_20;
 int q_Z_03_20;
 int q_Z_05_20;
 int q_Z_07_20;
 int q_Z_09_20;
 int q_Z_10_20;
 int q_Z_12_20;
 int q_Z_14_20;
 
 int q_Z_01_30;
 int q_Z_03_30;
 int q_Z_05_30;
 int q_Z_07_30;
 int q_Z_09_30;
 int q_Z_10_30;
 int q_Z_12_30;
 int q_Z_14_30;
 
 int q_Z_01_40;
 int q_Z_03_40;
 int q_Z_05_40;
 int q_Z_07_40;
 int q_Z_09_40;
 int q_Z_10_40;
 int q_Z_12_40;
 int q_Z_14_40;
 
 ///===== Higgs info
 
 double H_pT;
 double H_E;
 double H_pX;
 double H_pY;
 double H_pZ;
 double H_eta;
 double H_phi;
 double H_mt;
 
 ///====== Higgs WW events
 
 double mc_Q1pT;
 double mc_Q1pX;
 double mc_Q1pY;
 double mc_Q1pZ;
 double mc_Q1eta;
 double mc_Q1phi;
 double mc_Q1E;
 double mc_Q1pdgTag;
 double mc_Q1charge;
 
 
 double mc_Q2pT;
 double mc_Q2pX;
 double mc_Q2pY;
 double mc_Q2pZ;
 double mc_Q2eta;
 double mc_Q2phi;
 double mc_Q2E;
 double mc_Q2pdgTag;
 double mc_Q2charge;
 
 double mc_V1pT;
 double mc_V1pX;
 double mc_V1pY;
 double mc_V1pZ;
 double mc_V1eta;
 double mc_V1phi;
 double mc_V1E;
 double mc_V1pdgTag;
 double mc_V1charge;
 
 double mc_V2pT;
 double mc_V2pX;
 double mc_V2pY;
 double mc_V2pZ;
 double mc_V2eta;
 double mc_V2phi;
 double mc_V2E;
 double mc_V2pdgTag;
 double mc_V2charge;
 
 double mc_l1_V1pT;
 double mc_l1_V1pX;
 double mc_l1_V1pY;
 double mc_l1_V1pZ;
 double mc_l1_V1eta;
 double mc_l1_V1phi;
 double mc_l1_V1E;
 double mc_l1_V1pdgTag;
 double mc_l1_V1charge;
 
 double mc_l2_V1pT;
 double mc_l2_V1pX;
 double mc_l2_V1pY;
 double mc_l2_V1pZ;
 double mc_l2_V1eta;
 double mc_l2_V1phi;
 double mc_l2_V1E;
 double mc_l2_V1pdgTag;
 double mc_l2_V1charge;
 
 double mc_l1_V2pT;
 double mc_l1_V2pX;
 double mc_l1_V2pY;
 double mc_l1_V2pZ;
 double mc_l1_V2eta;
 double mc_l1_V2phi;
 double mc_l1_V2E;
 double mc_l1_V2pdgTag;
 double mc_l1_V2charge;
 
 double mc_l2_V2pT;
 double mc_l2_V2pX;
 double mc_l2_V2pY;
 double mc_l2_V2pZ;
 double mc_l2_V2eta;
 double mc_l2_V2phi;
 double mc_l2_V2E;
 double mc_l2_V2pdgTag;
 double mc_l2_V2charge;
 
 
 ///=== Other MC info
 int numPUMCoot;
 int numPUMCit;
 int numPUMC;
 
 
   
};


void InitializeTree(Variables_Gen&, const std::string& );

void FillTree(Variables_Gen& vars);
void FillEfficiencyTree(Variables_Gen& vars);

void SetEventVariables(Variables_Gen& vars, treeReader& reader);

void SetLeptonsVariables(Variables_Gen& vars, treeReader& reader,const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2);

void SetMetVariables(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2);

void SetQJetVariables(Variables_Gen& vars, treeReader& reader, const int& q1, const int& q2, const std::vector<int>& blacklistJet_forCJV);

void SetMCVariables(Variables_Gen& vars, treeReader& reader);

void SetMTVariable(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2);

void SetDPhiJetll(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2, const int& q1, const int& q2);
void SetDPhiJetll(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2, const std::vector<int> & blacklistJet);

void SetJDVariables(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2, const int& q1, const int& q2);

void FindAddJet (treeReader& reader,const int& q1,const int& q2, const std::vector<int>* blacklistJet_forCJV,std::vector<ROOT::Math::XYZTVector> & Jet_Candidate, const double& EtMin);
void SetThirdJetVariables(Variables_Gen& vars, std::vector<ROOT::Math::XYZTVector> & Jet_Candidate);

void SaveTree(Variables_Gen& vars);

struct pT_sort:
public std::binary_function< ROOT::Math::XYZTVector, ROOT::Math::XYZTVector, bool >
{
  bool operator() (ROOT::Math::XYZTVector x, ROOT::Math::XYZTVector y)
  {
    return x.pt() < y.pt();
  }
};
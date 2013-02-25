#ifndef VBFMuonClass_h
#define VBFMuonClass_h

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

#include "TLorentzVector.h"

#include "treeReader.h"

#include <string>



class VBFMuonClass {

 public :


  VBFMuonClass() ;
  VBFMuonClass( TTree *inputTree = 0) ;
  VBFMuonClass( TFile *inputFile = 0, std::string inputTreeName = "WJet") ;
  ~VBFMuonClass();

  void SetBranchAddressAndStatus(TTree * inputTree);

  TTree* GetTree();

  void SetTree(TTree* inputTree =0);

  void SetNewBranches(TTree* inputTree =0);  

  void InitializateVariables();

  TTree* fTree ;
  treeReader * fReader ;
  // new Branch for kinematic fit result

  float fit_mu_px ,  fit_mu_py ,  fit_mu_pz ,  fit_mu_e ;
  float fit_nv_px ,  fit_nv_py ,  fit_nv_pz ,  fit_nv_e ;
  float fit_subjet1_px ,  fit_subjet1_py ,  fit_subjet1_pz ,  fit_subjet1_e ;
  float fit_subjet2_px ,  fit_subjet2_py ,  fit_subjet2_pz ,  fit_subjet2_e ;
  float fit_lvj_m  ,  fit_lv_m , fit_j_m , fit_subjet1_m, fit_subjet2_m, fit_lvj_pt, fit_lvj_phi, fit_lvj_eta, fit_lvj_e, fit_chi2  ;
  int   fit_NDF   , fit_status ;

  float boosted_lvj_m , boosted_j_m, boosted_subjet1_m, boosted_subjet2_m, boosted_lv_m , boosted_lvj_pt, boosted_lvj_phi, boosted_lvj_eta, boosted_lvj_e ;
  float boostedW_lvj_m , boostedW_j_m, boostedW_subjet1_m, boostedW_subjet2_m, boostedW_lv_m , boostedW_lvj_pt, boostedW_lvj_phi, boostedW_lvj_eta, boostedW_lvj_e ;

  float boosted_wjj_ang_ha, boosted_wjj_ang_hb, boosted_wjj_ang_hs, boosted_wjj_ang_phi, boosted_wjj_ang_phia, boosted_wjj_ang_phib;

  
  // vbf variables for couple with high pt
  
  float vbf_maxpt_jj_e ,   vbf_maxpt_jj_pt ,   vbf_maxpt_jj_eta ,  vbf_maxpt_jj_phi , vbf_maxpt_jj_m ;
  float vbf_maxpt_j1_e ,   vbf_maxpt_j1_pt ,   vbf_maxpt_j1_eta ,  vbf_maxpt_j1_phi , vbf_maxpt_j1_m ;
  float vbf_maxpt_j2_e ,   vbf_maxpt_j2_pt ,   vbf_maxpt_j2_eta ,  vbf_maxpt_j2_phi , vbf_maxpt_j2_m ;
  float vbf_maxpt_jj_deta ,vbf_maxpt_jj_dphi; 

  float vbf_maxpt_j1_QGLikelihood,  vbf_maxpt_j2_QGLikelihood;

  float vbf_maxpt_j1_isPileUpLoose  , vbf_maxpt_j2_isPileUpLoose  ;
  float vbf_maxpt_j1_isPileUpMedium , vbf_maxpt_j2_isPileUpMedium ;
  float vbf_maxpt_j1_isPileUpTight  , vbf_maxpt_j2_isPileUpTight  ;

  int vbf_maxpt_jj_type,   vbf_maxpt_n_excj,   vbf_maxpt_n_exfj;

  float vbf_maxpt_j1_bDiscriminatorSSVHE , vbf_maxpt_j1_bDiscriminatorTCHE, vbf_maxpt_j1_bDiscriminatorCSV;
  float vbf_maxpt_j1_bDiscriminatorSSVHP, vbf_maxpt_j1_bDiscriminatorTCHP ;
  float vbf_maxpt_j2_bDiscriminatorSSVHE, vbf_maxpt_j2_bDiscriminatorTCHE, vbf_maxpt_j2_bDiscriminatorCSV;
  float vbf_maxpt_j2_bDiscriminatorSSVHP, vbf_maxpt_j2_bDiscriminatorTCHP ;

  float vbf_maxpt_j1_ChargedHadronEnergy,  vbf_maxpt_j1_ChargedHadronEnergyFrac, vbf_maxpt_j1_NeutralHadronEnergy;
  float vbf_maxpt_j1_NeutralHadronEnergyFrac, vbf_maxpt_j1_ChargedEmEnergy, vbf_maxpt_j1_ChargedEmEnergyFrac, vbf_maxpt_j1_ChargedMuEnergy;
  float vbf_maxpt_j1_ChargedMuEnergyFrac, vbf_maxpt_j1_NeutralEmEnergy,vbf_maxpt_j1_NeutralEmEnergyFrac, vbf_maxpt_j1_ChargedMultiplicity;
  float vbf_maxpt_j1_NeutralMultiplicity, vbf_maxpt_j1_MuonMultiplicity,vbf_maxpt_j1_PhotonEnergy, vbf_maxpt_j1_PhotonEnergyFraction;
  float vbf_maxpt_j1_ElectronEnergy, vbf_maxpt_j1_ElectronEnergyFraction,vbf_maxpt_j1_MuonEnergy, vbf_maxpt_j1_MuonEnergyFraction;
  float vbf_maxpt_j1_HFHadronEnergy, vbf_maxpt_j1_HFHadronEnergyFraction,vbf_maxpt_j1_HFEMEnergy, vbf_maxpt_j1_HFEMEnergyFraction;
  float vbf_maxpt_j1_ChargedHadronMultiplicity, vbf_maxpt_j1_NeutralHadronMultiplicity, vbf_maxpt_j1_PhotonMultiplicity;
  float vbf_maxpt_j1_ElectronMultiplicity,vbf_maxpt_j1_HFHadronMultiplicity;

  float vbf_maxpt_j2_ChargedHadronEnergy,  vbf_maxpt_j2_ChargedHadronEnergyFrac, vbf_maxpt_j2_NeutralHadronEnergy;
  float vbf_maxpt_j2_NeutralHadronEnergyFrac, vbf_maxpt_j2_ChargedEmEnergy, vbf_maxpt_j2_ChargedEmEnergyFrac, vbf_maxpt_j2_ChargedMuEnergy;
  float vbf_maxpt_j2_ChargedMuEnergyFrac, vbf_maxpt_j2_NeutralEmEnergy,vbf_maxpt_j2_NeutralEmEnergyFrac, vbf_maxpt_j2_ChargedMultiplicity;
  float vbf_maxpt_j2_NeutralMultiplicity, vbf_maxpt_j2_MuonMultiplicity,vbf_maxpt_j2_PhotonEnergy, vbf_maxpt_j2_PhotonEnergyFraction;
  float vbf_maxpt_j2_ElectronEnergy, vbf_maxpt_j2_ElectronEnergyFraction,vbf_maxpt_j2_MuonEnergy, vbf_maxpt_j2_MuonEnergyFraction;
  float vbf_maxpt_j2_HFHadronEnergy, vbf_maxpt_j2_HFHadronEnergyFraction,vbf_maxpt_j2_HFEMEnergy, vbf_maxpt_j2_HFEMEnergyFraction;
  float vbf_maxpt_j2_ChargedHadronMultiplicity, vbf_maxpt_j2_NeutralHadronMultiplicity, vbf_maxpt_j2_PhotonMultiplicity;
  float vbf_maxpt_j2_ElectronMultiplicity,vbf_maxpt_j2_HFHadronMultiplicity;

  // vbf variables for couple with high Deta  

  float vbf_maxDeta_jj_e ,   vbf_maxDeta_jj_pt ,   vbf_maxDeta_jj_eta ,  vbf_maxDeta_jj_phi , vbf_maxDeta_jj_m ;
  float vbf_maxDeta_j1_e ,   vbf_maxDeta_j1_pt ,   vbf_maxDeta_j1_eta ,  vbf_maxDeta_j1_phi , vbf_maxDeta_j1_m ;
  float vbf_maxDeta_j2_e ,   vbf_maxDeta_j2_pt ,   vbf_maxDeta_j2_eta ,  vbf_maxDeta_j2_phi , vbf_maxDeta_j2_m ;
  float vbf_maxDeta_jj_deta ,vbf_maxDeta_jj_dphi; 

  float vbf_maxDeta_j1_QGLikelihood, vbf_maxDeta_j2_QGLikelihood;

  float vbf_maxDeta_j1_isPileUpLoose  , vbf_maxDeta_j2_isPileUpLoose  ;
  float vbf_maxDeta_j1_isPileUpMedium , vbf_maxDeta_j2_isPileUpMedium ;
  float vbf_maxDeta_j1_isPileUpTight  , vbf_maxDeta_j2_isPileUpTight  ;

  int vbf_maxDeta_jj_type,   vbf_maxDeta_n_excj,   vbf_maxDeta_n_exfj;
 
  float vbf_maxDeta_j1_bDiscriminatorSSVHE, vbf_maxDeta_j1_bDiscriminatorTCHE, vbf_maxDeta_j1_bDiscriminatorCSV;
  float vbf_maxDeta_j1_bDiscriminatorSSVHP, vbf_maxDeta_j1_bDiscriminatorTCHP ;
  float vbf_maxDeta_j2_bDiscriminatorSSVHE, vbf_maxDeta_j2_bDiscriminatorTCHE, vbf_maxDeta_j2_bDiscriminatorCSV;
  float vbf_maxDeta_j2_bDiscriminatorSSVHP, vbf_maxDeta_j2_bDiscriminatorTCHP ;

  float vbf_maxDeta_j1_ChargedHadronEnergy,  vbf_maxDeta_j1_ChargedHadronEnergyFrac, vbf_maxDeta_j1_NeutralHadronEnergy;
  float vbf_maxDeta_j1_NeutralHadronEnergyFrac, vbf_maxDeta_j1_ChargedEmEnergy, vbf_maxDeta_j1_ChargedEmEnergyFrac, vbf_maxDeta_j1_ChargedMuEnergy;
  float vbf_maxDeta_j1_ChargedMuEnergyFrac, vbf_maxDeta_j1_NeutralEmEnergy,vbf_maxDeta_j1_NeutralEmEnergyFrac, vbf_maxDeta_j1_ChargedMultiplicity;
  float vbf_maxDeta_j1_NeutralMultiplicity, vbf_maxDeta_j1_MuonMultiplicity,vbf_maxDeta_j1_PhotonEnergy, vbf_maxDeta_j1_PhotonEnergyFraction;
  float vbf_maxDeta_j1_ElectronEnergy, vbf_maxDeta_j1_ElectronEnergyFraction,vbf_maxDeta_j1_MuonEnergy, vbf_maxDeta_j1_MuonEnergyFraction;
  float vbf_maxDeta_j1_HFHadronEnergy, vbf_maxDeta_j1_HFHadronEnergyFraction,vbf_maxDeta_j1_HFEMEnergy, vbf_maxDeta_j1_HFEMEnergyFraction;
  float vbf_maxDeta_j1_ChargedHadronMultiplicity, vbf_maxDeta_j1_NeutralHadronMultiplicity, vbf_maxDeta_j1_PhotonMultiplicity;
  float vbf_maxDeta_j1_ElectronMultiplicity,vbf_maxDeta_j1_HFHadronMultiplicity;

  float vbf_maxDeta_j2_ChargedHadronEnergy,  vbf_maxDeta_j2_ChargedHadronEnergyFrac, vbf_maxDeta_j2_NeutralHadronEnergy;
  float vbf_maxDeta_j2_NeutralHadronEnergyFrac, vbf_maxDeta_j2_ChargedEmEnergy, vbf_maxDeta_j2_ChargedEmEnergyFrac, vbf_maxDeta_j2_ChargedMuEnergy;
  float vbf_maxDeta_j2_ChargedMuEnergyFrac, vbf_maxDeta_j2_NeutralEmEnergy,vbf_maxDeta_j2_NeutralEmEnergyFrac, vbf_maxDeta_j2_ChargedMultiplicity;
  float vbf_maxDeta_j2_NeutralMultiplicity, vbf_maxDeta_j2_MuonMultiplicity,vbf_maxDeta_j2_PhotonEnergy, vbf_maxDeta_j2_PhotonEnergyFraction;
  float vbf_maxDeta_j2_ElectronEnergy, vbf_maxDeta_j2_ElectronEnergyFraction,vbf_maxDeta_j2_MuonEnergy, vbf_maxDeta_j2_MuonEnergyFraction;
  float vbf_maxDeta_j2_HFHadronEnergy, vbf_maxDeta_j2_HFHadronEnergyFraction,vbf_maxDeta_j2_HFEMEnergy, vbf_maxDeta_j2_HFEMEnergyFraction;
  float vbf_maxDeta_j2_ChargedHadronMultiplicity, vbf_maxDeta_j2_NeutralHadronMultiplicity, vbf_maxDeta_j2_PhotonMultiplicity;
  float vbf_maxDeta_j2_ElectronMultiplicity,vbf_maxDeta_j2_HFHadronMultiplicity;


  // vbf variables for couple with high Mjj  

  float vbf_maxMjj_jj_e ,   vbf_maxMjj_jj_pt ,   vbf_maxMjj_jj_eta ,  vbf_maxMjj_jj_phi , vbf_maxMjj_jj_m ;
  float vbf_maxMjj_j1_e ,   vbf_maxMjj_j1_pt ,   vbf_maxMjj_j1_eta ,  vbf_maxMjj_j1_phi , vbf_maxMjj_j1_m ;
  float vbf_maxMjj_j2_e ,   vbf_maxMjj_j2_pt ,   vbf_maxMjj_j2_eta ,  vbf_maxMjj_j2_phi , vbf_maxMjj_j2_m ;
  float vbf_maxMjj_jj_deta ,vbf_maxMjj_jj_dphi; 

  float vbf_maxMjj_j1_QGLikelihood, vbf_maxMjj_j2_QGLikelihood;

  float vbf_maxMjj_j1_isPileUpLoose  , vbf_maxMjj_j2_isPileUpLoose  ;
  float vbf_maxMjj_j1_isPileUpMedium , vbf_maxMjj_j2_isPileUpMedium ;
  float vbf_maxMjj_j1_isPileUpTight  , vbf_maxMjj_j2_isPileUpTight ;

  int vbf_maxMjj_jj_type,   vbf_maxMjj_n_excj,   vbf_maxMjj_n_exfj;

  float vbf_maxMjj_j1_bDiscriminatorSSVHE, vbf_maxMjj_j1_bDiscriminatorTCHE, vbf_maxMjj_j1_bDiscriminatorCSV;
  float vbf_maxMjj_j1_bDiscriminatorSSVHP, vbf_maxMjj_j1_bDiscriminatorTCHP ;
  float vbf_maxMjj_j2_bDiscriminatorSSVHE, vbf_maxMjj_j2_bDiscriminatorTCHE, vbf_maxMjj_j2_bDiscriminatorCSV;
  float vbf_maxMjj_j2_bDiscriminatorSSVHP, vbf_maxMjj_j2_bDiscriminatorTCHP ;

  float vbf_maxMjj_j1_ChargedHadronEnergy,  vbf_maxMjj_j1_ChargedHadronEnergyFrac, vbf_maxMjj_j1_NeutralHadronEnergy;
  float vbf_maxMjj_j1_NeutralHadronEnergyFrac, vbf_maxMjj_j1_ChargedEmEnergy, vbf_maxMjj_j1_ChargedEmEnergyFrac, vbf_maxMjj_j1_ChargedMuEnergy;
  float vbf_maxMjj_j1_ChargedMuEnergyFrac, vbf_maxMjj_j1_NeutralEmEnergy,vbf_maxMjj_j1_NeutralEmEnergyFrac, vbf_maxMjj_j1_ChargedMultiplicity;
  float vbf_maxMjj_j1_NeutralMultiplicity, vbf_maxMjj_j1_MuonMultiplicity,vbf_maxMjj_j1_PhotonEnergy, vbf_maxMjj_j1_PhotonEnergyFraction;
  float vbf_maxMjj_j1_ElectronEnergy, vbf_maxMjj_j1_ElectronEnergyFraction,vbf_maxMjj_j1_MuonEnergy, vbf_maxMjj_j1_MuonEnergyFraction;
  float vbf_maxMjj_j1_HFHadronEnergy, vbf_maxMjj_j1_HFHadronEnergyFraction,vbf_maxMjj_j1_HFEMEnergy, vbf_maxMjj_j1_HFEMEnergyFraction;
  float vbf_maxMjj_j1_ChargedHadronMultiplicity, vbf_maxMjj_j1_NeutralHadronMultiplicity, vbf_maxMjj_j1_PhotonMultiplicity;
  float vbf_maxMjj_j1_ElectronMultiplicity,vbf_maxMjj_j1_HFHadronMultiplicity;

  float vbf_maxMjj_j2_ChargedHadronEnergy,  vbf_maxMjj_j2_ChargedHadronEnergyFrac, vbf_maxMjj_j2_NeutralHadronEnergy;
  float vbf_maxMjj_j2_NeutralHadronEnergyFrac, vbf_maxMjj_j2_ChargedEmEnergy, vbf_maxMjj_j2_ChargedEmEnergyFrac, vbf_maxMjj_j2_ChargedMuEnergy;
  float vbf_maxMjj_j2_ChargedMuEnergyFrac, vbf_maxMjj_j2_NeutralEmEnergy,vbf_maxMjj_j2_NeutralEmEnergyFrac, vbf_maxMjj_j2_ChargedMultiplicity;
  float vbf_maxMjj_j2_NeutralMultiplicity, vbf_maxMjj_j2_MuonMultiplicity,vbf_maxMjj_j2_PhotonEnergy, vbf_maxMjj_j2_PhotonEnergyFraction;
  float vbf_maxMjj_j2_ElectronEnergy, vbf_maxMjj_j2_ElectronEnergyFraction,vbf_maxMjj_j2_MuonEnergy, vbf_maxMjj_j2_MuonEnergyFraction;
  float vbf_maxMjj_j2_HFHadronEnergy, vbf_maxMjj_j2_HFHadronEnergyFraction,vbf_maxMjj_j2_HFEMEnergy, vbf_maxMjj_j2_HFEMEnergyFraction;
  float vbf_maxMjj_j2_ChargedHadronMultiplicity, vbf_maxMjj_j2_NeutralHadronMultiplicity, vbf_maxMjj_j2_PhotonMultiplicity;
  float vbf_maxMjj_j2_ElectronMultiplicity,vbf_maxMjj_j2_HFHadronMultiplicity;
 
};

#endif

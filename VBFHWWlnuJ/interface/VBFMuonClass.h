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

  void SetBranchAddressAndStatus(TTree * inputTree = 0 );

  TTree* GetTree();

  void SetTree(TTree* inputTree =0);

  void SetReader(TTree* inputTree = 0);

  void SetNewBranches(TTree* inputTree =0);  

  void InitializateVariables();


  TTree* fTree ;
  treeReader * fReader ;
  
  //variables for the hadronic W
  
  int WHadposition;  
  int numberJetBin;

  float Hadronic_W_Jet_mass_uncorr;   
  float Hadronic_W_Jet_mass_tr_uncorr;  
  float Hadronic_W_Jet_mass_ft_uncorr;   
  float Hadronic_W_Jet_mass_pr_uncorr;   
  float Hadronic_W_Jet_massdrop_pr_uncorr;   
  float Hadronic_W_Jet_tau2tau1;  
  float Hadronic_W_Jet_tau1;  
  float Hadronic_W_Jet_tau2;   
  float Hadronic_W_Jet_tau3;   
  float Hadronic_W_Jet_tau4;   
  float Hadronic_W_Jet_pt;   
  float Hadronic_W_Jet_eta;   
  float Hadronic_W_Jet_phi;   
  float Hadronic_W_Jet_e;   
  float Hadronic_W_Jet_pt_tr_uncorr;   
  float Hadronic_W_Jet_pt_tr;   
  float Hadronic_W_Jet_eta_tr;   
  float Hadronic_W_Jet_phi_tr;   
  float Hadronic_W_Jet_e_tr;  
  float Hadronic_W_Jet_pt_ft_uncorr;   
  float Hadronic_W_Jet_pt_ft;   
  float Hadronic_W_Jet_eta_ft;   
  float Hadronic_W_Jet_phi_ft;   
  float Hadronic_W_Jet_e_ft;   
  float Hadronic_W_Jet_pt_pr_uncorr;   
  float Hadronic_W_Jet_pt_pr;   
  float Hadronic_W_Jet_eta_pr;   
  float Hadronic_W_Jet_phi_pr;   
  float Hadronic_W_Jet_e_pr;   
  float Hadronic_W_Jet_prsubjet1_px;  
  float Hadronic_W_Jet_prsubjet1_py;  
  float Hadronic_W_Jet_prsubjet1_pz;   
  float Hadronic_W_Jet_prsubjet1_e;   
  float Hadronic_W_Jet_prsubjet2_px;   
  float Hadronic_W_Jet_prsubjet2_py;   
  float Hadronic_W_Jet_prsubjet2_pz;   
  float Hadronic_W_Jet_prsubjet2_e;   
  float Hadronic_W_Jet_mass;  
  float Hadronic_W_Jet_mass_tr;   
  float Hadronic_W_Jet_mass_ft;  
  float Hadronic_W_Jet_mass_pr;   
  float Hadronic_W_Jet_massdrop;   
  float Hadronic_W_Jet_area;   
  float Hadronic_W_Jet_area_tr;   
  float Hadronic_W_Jet_area_ft;   
  float Hadronic_W_Jet_area_pr;   
  float Hadronic_W_Jet_jetconsituents;   
  float Hadronic_W_Jet_jetcharge;   
  float Hadronic_W_Jet_rcores;   
  float Hadronic_W_Jet_ptcores;   
  float Hadronic_W_Jet_planarflow;   
  float Hadronic_W_Jet_qjetmass;   
  float Hadronic_W_Jet_qjetmassdrop;   
  float Hadronic_W_Jet_deltaR_ljet;   
  float Hadronic_W_Jet_deltaphi_METjet;   
  float Hadronic_W_Jet_deltaphi_Vca8jet;   
  float Hadronic_W_Jet_rcores01;  
  float Hadronic_W_Jet_rcores02;  
  float Hadronic_W_Jet_rcores03;   
  float Hadronic_W_Jet_rcores04;   
  float Hadronic_W_Jet_rcores05;   
  float Hadronic_W_Jet_rcores06;   
  float Hadronic_W_Jet_rcores07;  
  float Hadronic_W_Jet_rcores08;   
  float Hadronic_W_Jet_rcores09;   
  float Hadronic_W_Jet_rcores10;  
  float Hadronic_W_Jet_rcores11;   
  float Hadronic_W_Jet_ptcores01;  
  float Hadronic_W_Jet_ptcores02;   
  float Hadronic_W_Jet_ptcores03;  
  float Hadronic_W_Jet_ptcores04;   
  float Hadronic_W_Jet_ptcores05;   
  float Hadronic_W_Jet_ptcores06;   
  float Hadronic_W_Jet_ptcores07;   
  float Hadronic_W_Jet_ptcores08;   
  float Hadronic_W_Jet_ptcores09;   
  float Hadronic_W_Jet_ptcores10;   
  float Hadronic_W_Jet_ptcores11;  
  float Hadronic_W_Jet_planarflow01;  
  float Hadronic_W_Jet_planarflow02;   
  float Hadronic_W_Jet_planarflow03;   
  float Hadronic_W_Jet_planarflow04;   
  float Hadronic_W_Jet_planarflow05;  
  float Hadronic_W_Jet_planarflow06;  
  float Hadronic_W_Jet_planarflow07;   
  float Hadronic_W_Jet_planarflow08;   
  float Hadronic_W_Jet_planarflow09;  
  float Hadronic_W_Jet_planarflow10;   
  float Hadronic_W_Jet_planarflow11;   
  float Hadronic_W_Jet_mass_sensi_tr;   
  float Hadronic_W_Jet_mass_sensi_ft;  
  float Hadronic_W_Jet_mass_sensi_pr;   
  float Hadronic_W_Jet_qjetmassvolatility;  
  float Hadronic_W_Jet_prsubjet1ptoverjetpt;   
  float Hadronic_W_Jet_prsubjet2ptoverjetpt; 
  float Hadronic_W_Jet_prsubjet1subjet2_deltaR;

  // new branch for pz of the neutrino
  float W_mass_type0_met , W_pz_type0_met ,  W_nu1_pz_type0_met, W_nu2_pz_type0_met ;
  float W_mass_type2_met , W_pz_type2_met ,  W_nu1_pz_type2_met, W_nu2_pz_type2_met ;

  float W_mass_type0, W_pz_type0,  W_nu1_pz_type0, W_nu2_pz_type0 ;
  float W_mass_type2, W_pz_type2,  W_nu1_pz_type2, W_nu2_pz_type2 ;
  
  // new Branch for kinematic fit result

  float fit_mu_px_type0 ,  fit_mu_py_type0 ,  fit_mu_pz_type0 ,  fit_mu_e_type0 ;
  float fit_nv_px_type0 ,  fit_nv_py_type0 ,  fit_nv_pz_type0 ,  fit_nv_e_type0 ;
  float fit_subjet1_px_type0 ,  fit_subjet1_py_type0 ,  fit_subjet1_pz_type0 ,  fit_subjet1_e_type0 ;
  float fit_subjet2_px_type0 ,  fit_subjet2_py_type0 ,  fit_subjet2_pz_type0 ,  fit_subjet2_e_type0 ;
  float fit_lvj_m_type0  ,  fit_lv_m_type0 , fit_j_m_type0 , fit_subjet1_m_type0, fit_subjet2_m_type0, fit_lvj_pt_type0, fit_lvj_phi_type0, fit_lvj_eta_type0, fit_lvj_e_type0, fit_chi2_type0  ;
  int   fit_NDF_type0   , fit_status_type0 ;

  float fit_mu_px_type2 ,  fit_mu_py_type2 ,  fit_mu_pz_type2 ,  fit_mu_e_type2 ;
  float fit_nv_px_type2 ,  fit_nv_py_type2 ,  fit_nv_pz_type2 ,  fit_nv_e_type2 ;
  float fit_subjet1_px_type2 ,  fit_subjet1_py_type2 ,  fit_subjet1_pz_type2 ,  fit_subjet1_e_type2 ;
  float fit_subjet2_px_type2 ,  fit_subjet2_py_type2 ,  fit_subjet2_pz_type2 ,  fit_subjet2_e_type2 ;
  float fit_lvj_m_type2  ,  fit_lv_m_type2 , fit_j_m_type2 , fit_subjet1_m_type2, fit_subjet2_m_type2, fit_lvj_pt_type2, fit_lvj_phi_type2, fit_lvj_eta_type2, fit_lvj_e_type2, fit_chi2_type2  ;
  int   fit_NDF_type2   , fit_status_type2 ;


  float fit_mu_px_type0_met ,  fit_mu_py_type0_met ,  fit_mu_pz_type0_met ,  fit_mu_e_type0_met ;
  float fit_nv_px_type0_met ,  fit_nv_py_type0_met ,  fit_nv_pz_type0_met ,  fit_nv_e_type0_met ;
  float fit_subjet1_px_type0_met ,  fit_subjet1_py_type0_met ,  fit_subjet1_pz_type0_met ,  fit_subjet1_e_type0_met ;
  float fit_subjet2_px_type0_met ,  fit_subjet2_py_type0_met ,  fit_subjet2_pz_type0_met ,  fit_subjet2_e_type0_met ;
  float fit_lvj_m_type0_met  ,  fit_lv_m_type0_met , fit_j_m_type0_met , fit_subjet1_m_type0_met, fit_subjet2_m_type0_met, fit_lvj_pt_type0_met, fit_lvj_phi_type0_met, fit_lvj_eta_type0_met, fit_lvj_e_type0_met, fit_chi2_type0_met  ;
  int   fit_NDF_type0_met   , fit_status_type0_met ;

  float fit_mu_px_type2_met ,  fit_mu_py_type2_met ,  fit_mu_pz_type2_met ,  fit_mu_e_type2_met ;
  float fit_nv_px_type2_met ,  fit_nv_py_type2_met ,  fit_nv_pz_type2_met ,  fit_nv_e_type2_met ;
  float fit_subjet1_px_type2_met ,  fit_subjet1_py_type2_met ,  fit_subjet1_pz_type2_met ,  fit_subjet1_e_type2_met ;
  float fit_subjet2_px_type2_met ,  fit_subjet2_py_type2_met ,  fit_subjet2_pz_type2_met ,  fit_subjet2_e_type2_met ;
  float fit_lvj_m_type2_met  ,  fit_lv_m_type2_met , fit_j_m_type2_met , fit_subjet1_m_type2_met, fit_subjet2_m_type2_met, fit_lvj_pt_type2_met, fit_lvj_phi_type2_met, fit_lvj_eta_type2_met, fit_lvj_e_type2_met, fit_chi2_type2_met  ;
  int   fit_NDF_type2_met   , fit_status_type2_met ;


  //////////////////////

  float boosted_lvj_m_type0 , boosted_j_m_type0, boosted_subjet1_m_type0, boosted_subjet2_m_type0, boosted_lv_m_type0 , boosted_lvj_pt_type0, boosted_lvj_phi_type0, boosted_lvj_eta_type0, boosted_lvj_e_type0 ;
  float boostedW_lvj_m_type0 , boostedW_j_m_type0, boostedW_subjet1_m_type0, boostedW_subjet2_m_type0, boostedW_lv_m_type0 , boostedW_lvj_pt_type0, boostedW_lvj_phi_type0, boostedW_lvj_eta_type0, boostedW_lvj_e_type0 ;

  float boosted_wjj_ang_ha_type0, boosted_wjj_ang_hb_type0, boosted_wjj_ang_hs_type0, boosted_wjj_ang_phi_type0, boosted_wjj_ang_phia_type0, boosted_wjj_ang_phib_type0;

  float boosted_lvj_m_type0_met , boosted_j_m_type0_met, boosted_subjet1_m_type0_met, boosted_subjet2_m_type0_met, boosted_lv_m_type0_met , boosted_lvj_pt_type0_met, boosted_lvj_phi_type0_met, boosted_lvj_eta_type0_met, boosted_lvj_e_type0_met ;
  float boostedW_lvj_m_type0_met , boostedW_j_m_type0_met, boostedW_subjet1_m_type0_met, boostedW_subjet2_m_type0_met, boostedW_lv_m_type0_met , boostedW_lvj_pt_type0_met, boostedW_lvj_phi_type0_met, boostedW_lvj_eta_type0_met, boostedW_lvj_e_type0_met ;

  float boosted_wjj_ang_ha_type0_met, boosted_wjj_ang_hb_type0_met, boosted_wjj_ang_hs_type0_met, boosted_wjj_ang_phi_type0_met, boosted_wjj_ang_phia_type0_met, boosted_wjj_ang_phib_type0_met;

  float boosted_lvj_m_type2 , boosted_j_m_type2, boosted_subjet1_m_type2, boosted_subjet2_m_type2, boosted_lv_m_type2 , boosted_lvj_pt_type2, boosted_lvj_phi_type2, boosted_lvj_eta_type2, boosted_lvj_e_type2 ;
  float boostedW_lvj_m_type2 , boostedW_j_m_type2, boostedW_subjet1_m_type2, boostedW_subjet2_m_type2, boostedW_lv_m_type2 , boostedW_lvj_pt_type2, boostedW_lvj_phi_type2, boostedW_lvj_eta_type2, boostedW_lvj_e_type2 ;

  float boosted_wjj_ang_ha_type2, boosted_wjj_ang_hb_type2, boosted_wjj_ang_hs_type2, boosted_wjj_ang_phi_type2, boosted_wjj_ang_phia_type2, boosted_wjj_ang_phib_type2;

  float boosted_lvj_m_type2_met , boosted_j_m_type2_met, boosted_subjet1_m_type2_met, boosted_subjet2_m_type2_met, boosted_lv_m_type2_met , boosted_lvj_pt_type2_met, boosted_lvj_phi_type2_met, boosted_lvj_eta_type2_met, boosted_lvj_e_type2_met ;
  float boostedW_lvj_m_type2_met , boostedW_j_m_type2_met, boostedW_subjet1_m_type2_met, boostedW_subjet2_m_type2_met, boostedW_lv_m_type2_met , boostedW_lvj_pt_type2_met, boostedW_lvj_phi_type2_met, boostedW_lvj_eta_type2_met, boostedW_lvj_e_type2_met ;

  float boosted_wjj_ang_ha_type2_met, boosted_wjj_ang_hb_type2_met, boosted_wjj_ang_hs_type2_met, boosted_wjj_ang_phi_type2_met, boosted_wjj_ang_phia_type2_met, boosted_wjj_ang_phib_type2_met;

  
  // vbf variables for couple with high pt
  
  float vbf_maxpt_jj_e ,   vbf_maxpt_jj_pt ,   vbf_maxpt_jj_eta ,  vbf_maxpt_jj_phi , vbf_maxpt_jj_m ;
  float vbf_maxpt_j1_e ,   vbf_maxpt_j1_pt ,   vbf_maxpt_j1_eta ,  vbf_maxpt_j1_phi , vbf_maxpt_j1_m ;
  float vbf_maxpt_j2_e ,   vbf_maxpt_j2_pt ,   vbf_maxpt_j2_eta ,  vbf_maxpt_j2_phi , vbf_maxpt_j2_m ;
  float vbf_maxpt_jj_deta ,vbf_maxpt_jj_dphi; 

  float vbf_maxpt_j1_QGLikelihood,  vbf_maxpt_j2_QGLikelihood;

  bool vbf_maxpt_j1_isPileUpLoose  , vbf_maxpt_j2_isPileUpLoose  ;
  bool vbf_maxpt_j1_isPileUpMedium , vbf_maxpt_j2_isPileUpMedium ;
  bool vbf_maxpt_j1_isPileUpTight  , vbf_maxpt_j2_isPileUpTight  ;

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

  bool vbf_maxDeta_j1_isPileUpLoose  , vbf_maxDeta_j2_isPileUpLoose  ;
  bool vbf_maxDeta_j1_isPileUpMedium , vbf_maxDeta_j2_isPileUpMedium ;
  bool vbf_maxDeta_j1_isPileUpTight  , vbf_maxDeta_j2_isPileUpTight  ;

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

  bool vbf_maxMjj_j1_isPileUpLoose  , vbf_maxMjj_j2_isPileUpLoose  ;
  bool vbf_maxMjj_j1_isPileUpMedium , vbf_maxMjj_j2_isPileUpMedium ;
  bool vbf_maxMjj_j1_isPileUpTight  , vbf_maxMjj_j2_isPileUpTight ;

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

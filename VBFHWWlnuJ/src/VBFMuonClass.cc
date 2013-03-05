#include "VBFMuonClass.h"


VBFMuonClass::VBFMuonClass(){}

VBFMuonClass::VBFMuonClass(TTree* inputTree){

  if(inputTree == 0) {
                       TFile* f = new TFile("/data2/rgerosa/RD_Tree_v1/RD_mu_STopS_Tbar_CMSSW532.root");
                       fTree = (TTree*) f -> Get("WJet");
  }
  else fTree = inputTree ;
}


VBFMuonClass::VBFMuonClass(TFile* inputFile, std::string inputTreeName){

  if(inputFile == 0) {
                       TFile* f = new TFile("/data2/rgerosa/RD_Tree_v1/RD_mu_STopS_Tbar_CMSSW532.root");
                       fTree = (TTree*) f -> Get("WJet");
  }
  else fTree = (TTree*) inputFile -> Get(inputTreeName.c_str());

}


VBFMuonClass::~VBFMuonClass(){

  delete fTree ; 
}

TTree* VBFMuonClass::GetTree(){

  return fTree ;
}


void VBFMuonClass::SetReader(TTree* inputTree){

  if(inputTree ==0) return ;
 
  fReader = new treeReader((TTree*)(fTree), true);

  SetBranchAddressAndStatus(fTree);

}

void VBFMuonClass::SetTree(TTree* inputTree){

  if(inputTree==0){  TFile* f = new TFile("/data2/rgerosa/RD_Tree_v1/RD_mu_STopS_Tbar_CMSSW532.root");
                     fTree = (TTree*) f -> Get("WJet");
  }
  
  else fTree = inputTree ;

  fReader = new treeReader((TTree*)(fTree), false);

   SetBranchAddressAndStatus(fTree);

}

void VBFMuonClass::SetBranchAddressAndStatus ( TTree* inputTree){


  fTree->SetBranchStatus("MassV*",0);
  fTree->SetBranchStatus("Mass*j*",0);
  fTree->SetBranchStatus("cosTheta*",0);
  fTree->SetBranchStatus("cosJackson*",0);
  fTree->SetBranchStatus("colorCorr*", 0);
  fTree->SetBranchStatus("JetGen*",0);
  fTree->SetBranchStatus("numGenJets*",0);
  fTree->SetBranchStatus("Photon*",0);
  fTree->SetBranchStatus("*Photon*",0);
  fTree->SetBranchStatus("W_Hb*",0);
  fTree->SetBranchStatus("W_Lepton*",0);
  fTree->SetBranchStatus("W_Met*",0);
  fTree->SetBranchStatus("W_tLepton*",0);
  fTree->SetBranchStatus("W_tMet*",0);

  fTree->SetBranchStatus("fit_*",0);
  fTree->SetBranchStatus("ang_*",0);
  fTree->SetBranchStatus("masslvjj",0);
  fTree->SetBranchStatus("ptlvjj",0);
  fTree->SetBranchStatus("ylvjj",0);
  fTree->SetBranchStatus("philvjj",0);
  fTree->SetBranchStatus("mva2j*",0);
  fTree->SetBranchStatus("mva3j*",0);
  fTree->SetBranchStatus("mvavbf*",0);
  fTree->SetBranchStatus("GroomedJet_number*",0);
  fTree->SetBranchStatus("GroomedJet*ak7",0);
  fTree->SetBranchStatus("boostedW_*",0);
  fTree->SetBranchStatus("vbf_*",0);
  fTree->SetBranchStatus("ttH_*",0);
}




void VBFMuonClass::SetNewBranches ( TTree* inputTree){

   inputTree->Branch("fit_mu_px",&fit_mu_px,"fit_mu_px/F");
   inputTree->Branch("fit_mu_py",&fit_mu_py,"fit_mu_py/F");
   inputTree->Branch("fit_mu_pz",&fit_mu_pz,"fit_mu_pz/F"); 
   inputTree->Branch("fit_mu_e",&fit_mu_e,"fit_mu_e/F");
   inputTree->Branch("fit_nv_px",&fit_mu_px,"fit_nv_px/F");
   inputTree->Branch("fit_nv_py",&fit_nv_py,"fit_nv_py/F");
   inputTree->Branch("fit_nv_pz",&fit_nv_pz,"fit_nv_pz/F");
   inputTree->Branch("fit_nv_e",&fit_nv_px,"fit_nv_e/F");
   inputTree->Branch("fit_subjet1_px",&fit_subjet1_px,"fit_subjet1_px/F");
   inputTree->Branch("fit_subjet1_py",&fit_subjet1_py,"fit_subjet1_py/F");
   inputTree->Branch("fit_subjet1_pz",&fit_subjet1_pz,"fit_subjet1_pz/F");
   inputTree->Branch("fit_subjet1_e",&fit_subjet1_e,"fit_subjet1_e/F");
   inputTree->Branch("fit_subjet2_px",&fit_subjet2_px,"fit_subjet2_px/F");
   inputTree->Branch("fit_subjet2_py",&fit_subjet2_py,"fit_subjet2_py/F");
   inputTree->Branch("fit_subjet2_pz",&fit_subjet2_pz,"fit_subjet2_pz/F");
   inputTree->Branch("fit_subjet2_e",&fit_subjet2_e,"fit_subjet2_e/F");
   inputTree->Branch("fit_subjet1_m",&fit_subjet1_m,"fit_subjet1_m/F");
   inputTree->Branch("fit_subjet2_m",&fit_subjet2_m,"fit_subjet2_m/F");
   inputTree->Branch("fit_lvj_m",&fit_lvj_m,"fit_lvj_m/F");
   inputTree->Branch("fit_lv_m",&fit_lv_m,"fit_lv_m/F");
   inputTree->Branch("fit_j_m",&fit_j_m,"fit_j_m/F");
   inputTree->Branch("fit_lvj_pt",&fit_lvj_pt,"fit_lvj_pt/F");
   inputTree->Branch("fit_lvj_eta", &fit_lvj_eta,"fit_lvj_eta/F");
   inputTree->Branch("fit_lvj_phi",&fit_lvj_phi,"fit_lvj_phi/F");
   inputTree->Branch("fit_lvj_e",&fit_lvj_e,"fit_lvj_e/F");
   inputTree->Branch("fit_chi2",&fit_chi2,"fit_chi2/F");
   inputTree->Branch("fit_NDF",&fit_NDF,"fit_NDF/I");
   inputTree->Branch("fit_status",&fit_status,"fit_status/I");


   inputTree->Branch("boosted_lvj_e",&boosted_lvj_e,"boosted_lvj_e/F");
   inputTree->Branch("boosted_lvj_pt",&boosted_lvj_pt,"boosted_lvj_pt/F");
   inputTree->Branch("boosted_lvj_eta",&boosted_lvj_eta,"boosted_lvj_eta/F");
   inputTree->Branch("boosted_lvj_phi",&boosted_lvj_phi,"boosted_lvj_phi/F");
   inputTree->Branch("boosted_lvj_m",&boosted_lvj_m,"boosted_lvj_m/F");
   inputTree->Branch("boosted_j_m",&boosted_j_m,"boosted_j_m/F");
   inputTree->Branch("boosted_subjet1_m",&boosted_subjet1_m,"boosted_subjet1_m/F");
   inputTree->Branch("boosted_subjet2_m",&boosted_subjet2_m,"boosted_subjet2_m/F");
   inputTree->Branch("boosted_lv_m",&boosted_lv_m,"boosted_lv_m/F");

   inputTree->Branch("boostedW_lvj_e",&boostedW_lvj_e,"boostedW_lvj_e/F");
   inputTree->Branch("boostedW_lvj_pt",&boostedW_lvj_pt,"boostedW_lvj_pt/F");
   inputTree->Branch("boostedW_lvj_eta",&boostedW_lvj_eta,"boostedW_lvj_eta/F");
   inputTree->Branch("boostedW_lvj_phi",&boostedW_lvj_phi,"boostedW_lvj_phi/F");
   inputTree->Branch("boostedW_lvj_m",&boostedW_lvj_m,"boostedW_lvj_m/F");
   inputTree->Branch("boostedW_j_m",&boostedW_j_m,"boostedW_j_m/F");
   inputTree->Branch("boostedW_subjet1_m",&boostedW_subjet1_m,"boostedW_subjet1_m/F");
   inputTree->Branch("boostedW_subjet2_m",&boostedW_subjet2_m,"boostedW_subjet2_m/F");
   inputTree->Branch("boostedW_lv_m",&boostedW_lv_m,"boostedW_lv_m/F");

   inputTree->Branch("boosted_wjj_ang_ha",&boosted_wjj_ang_ha,"boosted_wjj_ang_ha/F");
   inputTree->Branch("boosted_wjj_ang_hb",&boosted_wjj_ang_ha,"boosted_wjj_ang_ha/F");
   inputTree->Branch("boosted_wjj_ang_hs",&boosted_wjj_ang_hs,"boosted_wjj_ang_hs/F");
   inputTree->Branch("boosted_wjj_ang_phi",&boosted_wjj_ang_phi,"boosted_wjj_ang_phi/F");
   inputTree->Branch("boosted_wjj_ang_phia",&boosted_wjj_ang_phia,"boosted_wjj_ang_phia/F");
   inputTree->Branch("boosted_wjj_ang_phib",&boosted_wjj_ang_phib,"boosted_wjj_ang_phib/F");

   // max pt
 
   inputTree->Branch("vbf_maxpt_jj_e",&vbf_maxpt_jj_e,"vbf_maxpt_jj_e/F");
   inputTree->Branch("vbf_maxpt_jj_pt",&vbf_maxpt_jj_pt,"vbf_maxpt_jj_pt/F");
   inputTree->Branch("vbf_maxpt_jj_eta",&vbf_maxpt_jj_eta,"vbf_maxpt_jj_eta/F");
   inputTree->Branch("vbf_maxpt_jj_phi",&vbf_maxpt_jj_phi,"vbf_maxpt_jj_phi/F");
   inputTree->Branch("vbf_maxpt_jj_m",&vbf_maxpt_jj_m,"vbf_maxpt_jj_m/F");
  
   inputTree->Branch("vbf_maxpt_j1_e",&vbf_maxpt_j1_e,"vbf_maxpt_j1_e/F");
   inputTree->Branch("vbf_maxpt_j1_pt",&vbf_maxpt_j1_pt,"vbf_maxpt_j1_pt/F");
   inputTree->Branch("vbf_maxpt_j1_eta",&vbf_maxpt_j1_eta,"vbf_maxpt_j1_eta/F");
   inputTree->Branch("vbf_maxpt_j1_phi",&vbf_maxpt_j1_phi,"vbf_maxpt_j1_phi/F");
   inputTree->Branch("vbf_maxpt_j1_m",&vbf_maxpt_j1_m,"vbf_maxpt_j1_m/F");
  
   inputTree->Branch("vbf_maxpt_j2_e",&vbf_maxpt_j2_e,"vbf_maxpt_j2_e/F");
   inputTree->Branch("vbf_maxpt_j2_pt",&vbf_maxpt_j2_pt,"vbf_maxpt_j2_pt/F");
   inputTree->Branch("vbf_maxpt_j2_eta",&vbf_maxpt_j2_eta,"vbf_maxpt_j2_eta/F");
   inputTree->Branch("vbf_maxpt_j2_phi",&vbf_maxpt_j2_phi,"vbf_maxpt_j2_phi/F");
   inputTree->Branch("vbf_maxpt_j2_m",&vbf_maxpt_j2_m,"vbf_maxpt_j2_m/F");
  
   inputTree->Branch("vbf_maxpt_jj_deta",&vbf_maxpt_jj_deta,"vbf_maxpt_jj_deta/F");
   inputTree->Branch("vbf_maxpt_jj_dphi",&vbf_maxpt_jj_dphi,"vbf_maxpt_jj_dphi/F");
   
   inputTree->Branch("vbf_maxpt_j1_QGLikelihood",&vbf_maxpt_j1_QGLikelihood,"vbf_maxpt_j1_QGLikelihood/F");
   inputTree->Branch("vbf_maxpt_j2_QGLikelihood",&vbf_maxpt_j2_QGLikelihood,"vbf_maxpt_j2_QGLikelihood/F");

   
   inputTree->Branch("vbf_maxpt_j1_isPileUpLoose",&vbf_maxpt_j1_isPileUpLoose,"vbf_maxpt_j1_isPileUpLoose/O");
   inputTree->Branch("vbf_maxpt_j1_isPileUpMedium",&vbf_maxpt_j1_isPileUpMedium,"vbf_maxpt_j1_isPileUpMedium/O");
   inputTree->Branch("vbf_maxpt_j1_isPileUpTight",&vbf_maxpt_j1_isPileUpTight,"vbf_maxpt_j1_isPileUpTight/O");

   inputTree->Branch("vbf_maxpt_j2_isPileUpLoose",&vbf_maxpt_j2_isPileUpLoose,"vbf_maxpt_j2_isPileUpLoose/O");
   inputTree->Branch("vbf_maxpt_j2_isPileUpMedium",&vbf_maxpt_j2_isPileUpMedium,"vbf_maxpt_j2_isPileUpMedium/O");
   inputTree->Branch("vbf_maxpt_j2_isPileUpTight",&vbf_maxpt_j2_isPileUpTight,"vbf_maxpt_j2_isPileUpTight/O");

   inputTree->Branch("vbf_maxpt_jj_type",&vbf_maxpt_jj_type,"vbf_maxpt_jj_type/I");
   inputTree->Branch("vbf_maxpt_n_excj",&vbf_maxpt_n_excj,"vbf_maxpt_n_excj/I");
   inputTree->Branch("vbf_maxpt_n_exfj",&vbf_maxpt_n_exfj,"vbf_maxpt_n_exfj/I");

   inputTree->Branch("vbf_maxpt_j1_bDiscriminatorSSVHE",&vbf_maxpt_j1_bDiscriminatorSSVHE,"vbf_maxpt_j1_bDiscriminatorSSVHE/F");
   inputTree->Branch("vbf_maxpt_j1_bDiscriminatorTCHE",&vbf_maxpt_j1_bDiscriminatorTCHE,"vbf_maxpt_j1_bDiscriminatorTCHE/F");  
   inputTree->Branch("vbf_maxpt_j1_bDiscriminatorCSV",&vbf_maxpt_j1_bDiscriminatorCSV,"vbf_maxpt_j1_bDiscriminatorCSV/F");
   inputTree->Branch("vbf_maxpt_j1_bDiscriminatorSSVHP",&vbf_maxpt_j1_bDiscriminatorSSVHP,"vbf_maxpt_j1_bDiscriminatorSSVHP/F");
   inputTree->Branch("vbf_maxpt_j1_bDiscriminatorTCHP",&vbf_maxpt_j1_bDiscriminatorTCHP,"vbf_maxpt_j1_bDiscriminatorTCHP/F");

   inputTree->Branch("vbf_maxpt_j2_bDiscriminatorSSVHE",&vbf_maxpt_j2_bDiscriminatorSSVHE,"vbf_maxpt_j2_bDiscriminatorSSVHE/F");
   inputTree->Branch("vbf_maxpt_j2_bDiscriminatorTCHE",&vbf_maxpt_j2_bDiscriminatorTCHE,"vbf_maxpt_j2_bDiscriminatorTCHE/F");  
   inputTree->Branch("vbf_maxpt_j2_bDiscriminatorCSV",&vbf_maxpt_j2_bDiscriminatorCSV,"vbf_maxpt_j2_bDiscriminatorCSV/F");
   inputTree->Branch("vbf_maxpt_j2_bDiscriminatorSSVHP",&vbf_maxpt_j2_bDiscriminatorSSVHP,"vbf_maxpt_j2_bDiscriminatorSSVHP/F");
   inputTree->Branch("vbf_maxpt_j2_bDiscriminatorTCHP",&vbf_maxpt_j2_bDiscriminatorTCHP,"vbf_maxpt_j2_bDiscriminatorTCHP/F");

   inputTree->Branch("vbf_maxpt_j1_ChargedHadronEnergy",&vbf_maxpt_j1_ChargedHadronEnergy,"vbf_maxpt_j1_ChargedHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_ChargedHadronEnergyFrac",&vbf_maxpt_j1_ChargedHadronEnergyFrac,"vbf_maxpt_j1_ChargedHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j1_NeutralHadronEnergy",&vbf_maxpt_j1_NeutralHadronEnergy,"vbf_maxpt_j1_NeutralHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_NeutralHadronEnergyFrac",&vbf_maxpt_j1_NeutralHadronEnergyFrac,"vbf_maxpt_j1_NeutralHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j1_ChargedEmEnergy",&vbf_maxpt_j1_ChargedEmEnergy,"vbf_maxpt_j1_ChargedEmEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_ChargedEmEnergyFrac",&vbf_maxpt_j1_ChargedEmEnergyFrac,"vbf_maxpt_j1_ChargedEmEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j1_ChargedMuEnergy",&vbf_maxpt_j1_ChargedMuEnergy,"vbf_maxpt_j1_ChargedMuEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_ChargedMuEnergyFrac",&vbf_maxpt_j1_ChargedMuEnergyFrac,"vbf_maxpt_j1_ChargedMuEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j1_NeutralEmEnergy",&vbf_maxpt_j1_NeutralEmEnergy,"vbf_maxpt_j1_NeutralEmEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_NeutralEmEnergyFrac",&vbf_maxpt_j1_NeutralEmEnergyFrac,"vbf_maxpt_j1_NeutralEmEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j1_ChargedMultiplicity",&vbf_maxpt_j1_ChargedMultiplicity,"vbf_maxpt_j1_ChargedMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j1_NeutralMultiplicity",&vbf_maxpt_j1_NeutralMultiplicity,"vbf_maxpt_j1_NeutralMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j1_MuonMultiplicity",&vbf_maxpt_j1_MuonMultiplicity,"vbf_maxpt_j1_MuonMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j1_PhotonEnergy",&vbf_maxpt_j1_PhotonEnergy,"vbf_maxpt_j1_PhotonEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_PhotonEnergyFraction",&vbf_maxpt_j1_PhotonEnergyFraction,"vbf_maxpt_j1_PhotonEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j1_ElectronEnergy",&vbf_maxpt_j1_ElectronEnergy,"vbf_maxpt_j1_ElectronEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_ElectronEnergyFraction",&vbf_maxpt_j1_ElectronEnergyFraction,"vbf_maxpt_j1_ElectronEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j1_HFHadronEnergy",&vbf_maxpt_j1_HFHadronEnergy,"vbf_maxpt_j1_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_HFHadronEnergyFraction",&vbf_maxpt_j1_HFHadronEnergyFraction,"vbf_maxpt_j1_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j1_HFEMEnergy",&vbf_maxpt_j1_HFEMEnergy,"vbf_maxpt_j1_HFEMEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_HFHadronEnergy",&vbf_maxpt_j1_HFHadronEnergy,"vbf_maxpt_j1_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j1_HFHadronEnergyFraction",&vbf_maxpt_j1_HFHadronEnergyFraction,"vbf_maxpt_j1_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j1_ChargedHadronMultiplicity",&vbf_maxpt_j1_ChargedHadronMultiplicity,"vbf_maxpt_j1_ChargedHadronMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j1_NeutralHadronMultiplicity",&vbf_maxpt_j1_NeutralHadronMultiplicity,"vbf_maxpt_j1_NeutralHadronMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j1_PhotonMultiplicity",&vbf_maxpt_j1_PhotonMultiplicity,"vbf_maxpt_j1_PhotonMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j1_ElectronMultiplicity",&vbf_maxpt_j1_ElectronMultiplicity,"vbf_maxpt_j1_ElectronMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j1_HFHadronMultiplicity",&vbf_maxpt_j1_HFHadronMultiplicity,"vbf_maxpt_j1_HFHadronMultiplicity/F");

   inputTree->Branch("vbf_maxpt_j2_ChargedHadronEnergy",&vbf_maxpt_j2_ChargedHadronEnergy,"vbf_maxpt_j2_ChargedHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_ChargedHadronEnergyFrac",&vbf_maxpt_j2_ChargedHadronEnergyFrac,"vbf_maxpt_j2_ChargedHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j2_NeutralHadronEnergy",&vbf_maxpt_j2_NeutralHadronEnergy,"vbf_maxpt_j2_NeutralHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_NeutralHadronEnergyFrac",&vbf_maxpt_j2_NeutralHadronEnergyFrac,"vbf_maxpt_j2_NeutralHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j2_ChargedEmEnergy",&vbf_maxpt_j2_ChargedEmEnergy,"vbf_maxpt_j2_ChargedEmEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_ChargedEmEnergyFrac",&vbf_maxpt_j2_ChargedEmEnergyFrac,"vbf_maxpt_j2_ChargedEmEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j2_ChargedMuEnergy",&vbf_maxpt_j2_ChargedMuEnergy,"vbf_maxpt_j2_ChargedMuEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_ChargedMuEnergyFrac",&vbf_maxpt_j2_ChargedMuEnergyFrac,"vbf_maxpt_j2_ChargedMuEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j2_NeutralEmEnergy",&vbf_maxpt_j2_NeutralEmEnergy,"vbf_maxpt_j2_NeutralEmEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_NeutralEmEnergyFrac",&vbf_maxpt_j2_NeutralEmEnergyFrac,"vbf_maxpt_j2_NeutralEmEnergyFrac/F");
   inputTree->Branch("vbf_maxpt_j2_ChargedMultiplicity",&vbf_maxpt_j2_ChargedMultiplicity,"vbf_maxpt_j2_ChargedMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j2_NeutralMultiplicity",&vbf_maxpt_j2_NeutralMultiplicity,"vbf_maxpt_j2_NeutralMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j2_MuonMultiplicity",&vbf_maxpt_j2_MuonMultiplicity,"vbf_maxpt_j2_MuonMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j2_PhotonEnergy",&vbf_maxpt_j2_PhotonEnergy,"vbf_maxpt_j2_PhotonEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_PhotonEnergyFraction",&vbf_maxpt_j2_PhotonEnergyFraction,"vbf_maxpt_j2_PhotonEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j2_ElectronEnergy",&vbf_maxpt_j2_ElectronEnergy,"vbf_maxpt_j2_ElectronEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_ElectronEnergyFraction",&vbf_maxpt_j2_ElectronEnergyFraction,"vbf_maxpt_j2_ElectronEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j2_HFHadronEnergy",&vbf_maxpt_j2_HFHadronEnergy,"vbf_maxpt_j2_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_HFHadronEnergyFraction",&vbf_maxpt_j2_HFHadronEnergyFraction,"vbf_maxpt_j2_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j2_HFEMEnergy",&vbf_maxpt_j2_HFEMEnergy,"vbf_maxpt_j2_HFEMEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_HFHadronEnergy",&vbf_maxpt_j2_HFHadronEnergy,"vbf_maxpt_j2_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxpt_j2_HFHadronEnergyFraction",&vbf_maxpt_j2_HFHadronEnergyFraction,"vbf_maxpt_j2_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxpt_j2_ChargedHadronMultiplicity",&vbf_maxpt_j2_ChargedHadronMultiplicity,"vbf_maxpt_j2_ChargedHadronMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j2_NeutralHadronMultiplicity",&vbf_maxpt_j2_NeutralHadronMultiplicity,"vbf_maxpt_j2_NeutralHadronMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j2_PhotonMultiplicity",&vbf_maxpt_j2_PhotonMultiplicity,"vbf_maxpt_j2_PhotonMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j2_ElectronMultiplicity",&vbf_maxpt_j2_ElectronMultiplicity,"vbf_maxpt_j2_ElectronMultiplicity/F");
   inputTree->Branch("vbf_maxpt_j2_HFHadronMultiplicity",&vbf_maxpt_j2_HFHadronMultiplicity,"vbf_maxpt_j2_HFHadronMultiplicity/F");

   // max Deta pair

   inputTree->Branch("vbf_maxDeta_jj_e",&vbf_maxDeta_jj_e,"vbf_maxDeta_jj_e/F");
   inputTree->Branch("vbf_maxDeta_jj_pt",&vbf_maxDeta_jj_pt,"vbf_maxDeta_jj_pt/F");
   inputTree->Branch("vbf_maxDeta_jj_eta",&vbf_maxDeta_jj_eta,"vbf_maxDeta_jj_eta/F");
   inputTree->Branch("vbf_maxDeta_jj_phi",&vbf_maxDeta_jj_phi,"vbf_maxDeta_jj_phi/F");
   inputTree->Branch("vbf_maxDeta_jj_m",&vbf_maxDeta_jj_m,"vbf_maxDeta_jj_m/F");
  
   inputTree->Branch("vbf_maxDeta_j1_e",&vbf_maxDeta_j1_e,"vbf_maxDeta_j1_e/F");
   inputTree->Branch("vbf_maxDeta_j1_pt",&vbf_maxDeta_j1_pt,"vbf_maxDeta_j1_pt/F");
   inputTree->Branch("vbf_maxDeta_j1_eta",&vbf_maxDeta_j1_eta,"vbf_maxDeta_j1_eta/F");
   inputTree->Branch("vbf_maxDeta_j1_phi",&vbf_maxDeta_j1_phi,"vbf_maxDeta_j1_phi/F");
   inputTree->Branch("vbf_maxDeta_j1_m",&vbf_maxDeta_j1_m,"vbf_maxDeta_j1_m/F");
  
   inputTree->Branch("vbf_maxDeta_j2_e",&vbf_maxDeta_j2_e,"vbf_maxDeta_j2_e/F");
   inputTree->Branch("vbf_maxDeta_j2_pt",&vbf_maxDeta_j2_pt,"vbf_maxDeta_j2_pt/F");
   inputTree->Branch("vbf_maxDeta_j2_eta",&vbf_maxDeta_j2_eta,"vbf_maxDeta_j2_eta/F");
   inputTree->Branch("vbf_maxDeta_j2_phi",&vbf_maxDeta_j2_phi,"vbf_maxDeta_j2_phi/F");
   inputTree->Branch("vbf_maxDeta_j2_m",&vbf_maxDeta_j2_m,"vbf_maxDeta_j2_m/F");
  
   inputTree->Branch("vbf_maxDeta_jj_deta",&vbf_maxDeta_jj_deta,"vbf_maxDeta_jj_deta/F");
   inputTree->Branch("vbf_maxDeta_jj_dphi",&vbf_maxDeta_jj_dphi,"vbf_maxDeta_jj_dphi/F");
   
   inputTree->Branch("vbf_maxDeta_j1_QGLikelihood",&vbf_maxDeta_j1_QGLikelihood,"vbf_maxDeta_j1_QGLikelihood/F");
   inputTree->Branch("vbf_maxDeta_j2_QGLikelihood",&vbf_maxDeta_j2_QGLikelihood,"vbf_maxDeta_j2_QGLikelihood/F");

   
   inputTree->Branch("vbf_maxDeta_j1_isPileUpLoose",&vbf_maxDeta_j1_isPileUpLoose,"vbf_maxDeta_j1_isPileUpLoose/O");
   inputTree->Branch("vbf_maxDeta_j1_isPileUpMedium",&vbf_maxDeta_j1_isPileUpMedium,"vbf_maxDeta_j1_isPileUpMedium/O");
   inputTree->Branch("vbf_maxDeta_j1_isPileUpTight",&vbf_maxDeta_j1_isPileUpTight,"vbf_maxDeta_j1_isPileUpTight/O");

   inputTree->Branch("vbf_maxDeta_j2_isPileUpLoose",&vbf_maxDeta_j2_isPileUpLoose,"vbf_maxDeta_j2_isPileUpLoose/O");
   inputTree->Branch("vbf_maxDeta_j2_isPileUpMedium",&vbf_maxDeta_j2_isPileUpMedium,"vbf_maxDeta_j2_isPileUpMedium/O");
   inputTree->Branch("vbf_maxDeta_j2_isPileUpTight",&vbf_maxDeta_j2_isPileUpTight,"vbf_maxDeta_j2_isPileUpTight/O");

   inputTree->Branch("vbf_maxDeta_jj_type",&vbf_maxDeta_jj_type,"vbf_maxDeta_jj_type/I");
   inputTree->Branch("vbf_maxDeta_n_excj",&vbf_maxDeta_n_excj,"vbf_maxDeta_n_excj/I");
   inputTree->Branch("vbf_maxDeta_n_exfj",&vbf_maxDeta_n_exfj,"vbf_maxDeta_n_exfj/I");

   inputTree->Branch("vbf_maxDeta_j1_bDiscriminatorSSVHE",&vbf_maxDeta_j1_bDiscriminatorSSVHE,"vbf_maxDeta_j1_bDiscriminatorSSVHE/F");
   inputTree->Branch("vbf_maxDeta_j1_bDiscriminatorTCHE",&vbf_maxDeta_j1_bDiscriminatorTCHE,"vbf_maxDeta_j1_bDiscriminatorTCHE/F");  
   inputTree->Branch("vbf_maxDeta_j1_bDiscriminatorCSV",&vbf_maxDeta_j1_bDiscriminatorCSV,"vbf_maxDeta_j1_bDiscriminatorCSV/F");
   inputTree->Branch("vbf_maxDeta_j1_bDiscriminatorSSVHP",&vbf_maxDeta_j1_bDiscriminatorSSVHP,"vbf_maxDeta_j1_bDiscriminatorSSVHP/F");
   inputTree->Branch("vbf_maxDeta_j1_bDiscriminatorTCHP",&vbf_maxDeta_j1_bDiscriminatorTCHP,"vbf_maxDeta_j1_bDiscriminatorTCHP/F");

   inputTree->Branch("vbf_maxDeta_j2_bDiscriminatorSSVHE",&vbf_maxDeta_j2_bDiscriminatorSSVHE,"vbf_maxDeta_j2_bDiscriminatorSSVHE/F");
   inputTree->Branch("vbf_maxDeta_j2_bDiscriminatorTCHE",&vbf_maxDeta_j2_bDiscriminatorTCHE,"vbf_maxDeta_j2_bDiscriminatorTCHE/F");  
   inputTree->Branch("vbf_maxDeta_j2_bDiscriminatorCSV",&vbf_maxDeta_j2_bDiscriminatorCSV,"vbf_maxDeta_j2_bDiscriminatorCSV/F");
   inputTree->Branch("vbf_maxDeta_j2_bDiscriminatorSSVHP",&vbf_maxDeta_j2_bDiscriminatorSSVHP,"vbf_maxDeta_j2_bDiscriminatorSSVHP/F");
   inputTree->Branch("vbf_maxDeta_j2_bDiscriminatorTCHP",&vbf_maxDeta_j2_bDiscriminatorTCHP,"vbf_maxDeta_j2_bDiscriminatorTCHP/F");

   inputTree->Branch("vbf_maxDeta_j1_ChargedHadronEnergy",&vbf_maxDeta_j1_ChargedHadronEnergy,"vbf_maxDeta_j1_ChargedHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_ChargedHadronEnergyFrac",&vbf_maxDeta_j1_ChargedHadronEnergyFrac,"vbf_maxDeta_j1_ChargedHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j1_NeutralHadronEnergy",&vbf_maxDeta_j1_NeutralHadronEnergy,"vbf_maxDeta_j1_NeutralHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_NeutralHadronEnergyFrac",&vbf_maxDeta_j1_NeutralHadronEnergyFrac,"vbf_maxDeta_j1_NeutralHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j1_ChargedEmEnergy",&vbf_maxDeta_j1_ChargedEmEnergy,"vbf_maxDeta_j1_ChargedEmEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_ChargedEmEnergyFrac",&vbf_maxDeta_j1_ChargedEmEnergyFrac,"vbf_maxDeta_j1_ChargedEmEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j1_ChargedMuEnergy",&vbf_maxDeta_j1_ChargedMuEnergy,"vbf_maxDeta_j1_ChargedMuEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_ChargedMuEnergyFrac",&vbf_maxDeta_j1_ChargedMuEnergyFrac,"vbf_maxDeta_j1_ChargedMuEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j1_NeutralEmEnergy",&vbf_maxDeta_j1_NeutralEmEnergy,"vbf_maxDeta_j1_NeutralEmEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_NeutralEmEnergyFrac",&vbf_maxDeta_j1_NeutralEmEnergyFrac,"vbf_maxDeta_j1_NeutralEmEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j1_ChargedMultiplicity",&vbf_maxDeta_j1_ChargedMultiplicity,"vbf_maxDeta_j1_ChargedMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j1_NeutralMultiplicity",&vbf_maxDeta_j1_NeutralMultiplicity,"vbf_maxDeta_j1_NeutralMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j1_MuonMultiplicity",&vbf_maxDeta_j1_MuonMultiplicity,"vbf_maxDeta_j1_MuonMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j1_PhotonEnergy",&vbf_maxDeta_j1_PhotonEnergy,"vbf_maxDeta_j1_PhotonEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_PhotonEnergyFraction",&vbf_maxDeta_j1_PhotonEnergyFraction,"vbf_maxDeta_j1_PhotonEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j1_ElectronEnergy",&vbf_maxDeta_j1_ElectronEnergy,"vbf_maxDeta_j1_ElectronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_ElectronEnergyFraction",&vbf_maxDeta_j1_ElectronEnergyFraction,"vbf_maxDeta_j1_ElectronEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j1_HFHadronEnergy",&vbf_maxDeta_j1_HFHadronEnergy,"vbf_maxDeta_j1_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_HFHadronEnergyFraction",&vbf_maxDeta_j1_HFHadronEnergyFraction,"vbf_maxDeta_j1_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j1_HFEMEnergy",&vbf_maxDeta_j1_HFEMEnergy,"vbf_maxDeta_j1_HFEMEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_HFHadronEnergy",&vbf_maxDeta_j1_HFHadronEnergy,"vbf_maxDeta_j1_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j1_HFHadronEnergyFraction",&vbf_maxDeta_j1_HFHadronEnergyFraction,"vbf_maxDeta_j1_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j1_ChargedHadronMultiplicity",&vbf_maxDeta_j1_ChargedHadronMultiplicity,"vbf_maxDeta_j1_ChargedHadronMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j1_NeutralHadronMultiplicity",&vbf_maxDeta_j1_NeutralHadronMultiplicity,"vbf_maxDeta_j1_NeutralHadronMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j1_PhotonMultiplicity",&vbf_maxDeta_j1_PhotonMultiplicity,"vbf_maxDeta_j1_PhotonMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j1_ElectronMultiplicity",&vbf_maxDeta_j1_ElectronMultiplicity,"vbf_maxDeta_j1_ElectronMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j1_HFHadronMultiplicity",&vbf_maxDeta_j1_HFHadronMultiplicity,"vbf_maxDeta_j1_HFHadronMultiplicity/F");

   inputTree->Branch("vbf_maxDeta_j2_ChargedHadronEnergy",&vbf_maxDeta_j2_ChargedHadronEnergy,"vbf_maxDeta_j2_ChargedHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_ChargedHadronEnergyFrac",&vbf_maxDeta_j2_ChargedHadronEnergyFrac,"vbf_maxDeta_j2_ChargedHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j2_NeutralHadronEnergy",&vbf_maxDeta_j2_NeutralHadronEnergy,"vbf_maxDeta_j2_NeutralHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_NeutralHadronEnergyFrac",&vbf_maxDeta_j2_NeutralHadronEnergyFrac,"vbf_maxDeta_j2_NeutralHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j2_ChargedEmEnergy",&vbf_maxDeta_j2_ChargedEmEnergy,"vbf_maxDeta_j2_ChargedEmEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_ChargedEmEnergyFrac",&vbf_maxDeta_j2_ChargedEmEnergyFrac,"vbf_maxDeta_j2_ChargedEmEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j2_ChargedMuEnergy",&vbf_maxDeta_j2_ChargedMuEnergy,"vbf_maxDeta_j2_ChargedMuEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_ChargedMuEnergyFrac",&vbf_maxDeta_j2_ChargedMuEnergyFrac,"vbf_maxDeta_j2_ChargedMuEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j2_NeutralEmEnergy",&vbf_maxDeta_j2_NeutralEmEnergy,"vbf_maxDeta_j2_NeutralEmEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_NeutralEmEnergyFrac",&vbf_maxDeta_j2_NeutralEmEnergyFrac,"vbf_maxDeta_j2_NeutralEmEnergyFrac/F");
   inputTree->Branch("vbf_maxDeta_j2_ChargedMultiplicity",&vbf_maxDeta_j2_ChargedMultiplicity,"vbf_maxDeta_j2_ChargedMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j2_NeutralMultiplicity",&vbf_maxDeta_j2_NeutralMultiplicity,"vbf_maxDeta_j2_NeutralMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j2_MuonMultiplicity",&vbf_maxDeta_j2_MuonMultiplicity,"vbf_maxDeta_j2_MuonMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j2_PhotonEnergy",&vbf_maxDeta_j2_PhotonEnergy,"vbf_maxDeta_j2_PhotonEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_PhotonEnergyFraction",&vbf_maxDeta_j2_PhotonEnergyFraction,"vbf_maxDeta_j2_PhotonEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j2_ElectronEnergy",&vbf_maxDeta_j2_ElectronEnergy,"vbf_maxDeta_j2_ElectronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_ElectronEnergyFraction",&vbf_maxDeta_j2_ElectronEnergyFraction,"vbf_maxDeta_j2_ElectronEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j2_HFHadronEnergy",&vbf_maxDeta_j2_HFHadronEnergy,"vbf_maxDeta_j2_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_HFHadronEnergyFraction",&vbf_maxDeta_j2_HFHadronEnergyFraction,"vbf_maxDeta_j2_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j2_HFEMEnergy",&vbf_maxDeta_j2_HFEMEnergy,"vbf_maxDeta_j2_HFEMEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_HFHadronEnergy",&vbf_maxDeta_j2_HFHadronEnergy,"vbf_maxDeta_j2_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxDeta_j2_HFHadronEnergyFraction",&vbf_maxDeta_j2_HFHadronEnergyFraction,"vbf_maxDeta_j2_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxDeta_j2_ChargedHadronMultiplicity",&vbf_maxDeta_j2_ChargedHadronMultiplicity,"vbf_maxDeta_j2_ChargedHadronMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j2_NeutralHadronMultiplicity",&vbf_maxDeta_j2_NeutralHadronMultiplicity,"vbf_maxDeta_j2_NeutralHadronMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j2_PhotonMultiplicity",&vbf_maxDeta_j2_PhotonMultiplicity,"vbf_maxDeta_j2_PhotonMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j2_ElectronMultiplicity",&vbf_maxDeta_j2_ElectronMultiplicity,"vbf_maxDeta_j2_ElectronMultiplicity/F");
   inputTree->Branch("vbf_maxDeta_j2_HFHadronMultiplicity",&vbf_maxDeta_j2_HFHadronMultiplicity,"vbf_maxDeta_j2_HFHadronMultiplicity/F");

   // max maxMjj pair

   inputTree->Branch("vbf_maxMjj_jj_e",&vbf_maxMjj_jj_e,"vbf_maxMjj_jj_e/F");
   inputTree->Branch("vbf_maxMjj_jj_pt",&vbf_maxMjj_jj_pt,"vbf_maxMjj_jj_pt/F");
   inputTree->Branch("vbf_maxMjj_jj_eta",&vbf_maxMjj_jj_eta,"vbf_maxMjj_jj_eta/F");
   inputTree->Branch("vbf_maxMjj_jj_phi",&vbf_maxMjj_jj_phi,"vbf_maxMjj_jj_phi/F");
   inputTree->Branch("vbf_maxMjj_jj_m",&vbf_maxMjj_jj_m,"vbf_maxMjj_jj_m/F");
  
   inputTree->Branch("vbf_maxMjj_j1_e",&vbf_maxMjj_j1_e,"vbf_maxMjj_j1_e/F");
   inputTree->Branch("vbf_maxMjj_j1_pt",&vbf_maxMjj_j1_pt,"vbf_maxMjj_j1_pt/F");
   inputTree->Branch("vbf_maxMjj_j1_eta",&vbf_maxMjj_j1_eta,"vbf_maxMjj_j1_eta/F");
   inputTree->Branch("vbf_maxMjj_j1_phi",&vbf_maxMjj_j1_phi,"vbf_maxMjj_j1_phi/F");
   inputTree->Branch("vbf_maxMjj_j1_m",&vbf_maxMjj_j1_m,"vbf_maxMjj_j1_m/F");
  
   inputTree->Branch("vbf_maxMjj_j2_e",&vbf_maxMjj_j2_e,"vbf_maxMjj_j2_e/F");
   inputTree->Branch("vbf_maxMjj_j2_pt",&vbf_maxMjj_j2_pt,"vbf_maxMjj_j2_pt/F");
   inputTree->Branch("vbf_maxMjj_j2_eta",&vbf_maxMjj_j2_eta,"vbf_maxMjj_j2_eta/F");
   inputTree->Branch("vbf_maxMjj_j2_phi",&vbf_maxMjj_j2_phi,"vbf_maxMjj_j2_phi/F");
   inputTree->Branch("vbf_maxMjj_j2_m",&vbf_maxMjj_j2_m,"vbf_maxMjj_j2_m/F");
  
   inputTree->Branch("vbf_maxMjj_jj_deta",&vbf_maxMjj_jj_deta,"vbf_maxMjj_jj_deta/F");
   inputTree->Branch("vbf_maxMjj_jj_dphi",&vbf_maxMjj_jj_dphi,"vbf_maxMjj_jj_dphi/F");
   
   inputTree->Branch("vbf_maxMjj_j1_QGLikelihood",&vbf_maxMjj_j1_QGLikelihood,"vbf_maxMjj_j1_QGLikelihood/F");
   inputTree->Branch("vbf_maxMjj_j2_QGLikelihood",&vbf_maxMjj_j2_QGLikelihood,"vbf_maxMjj_j2_QGLikelihood/F");

   
   inputTree->Branch("vbf_maxMjj_j1_isPileUpLoose",&vbf_maxMjj_j1_isPileUpLoose,"vbf_maxMjj_j1_isPileUpLoose/O");
   inputTree->Branch("vbf_maxMjj_j1_isPileUpMedium",&vbf_maxMjj_j1_isPileUpMedium,"vbf_maxMjj_j1_isPileUpMedium/O");
   inputTree->Branch("vbf_maxMjj_j1_isPileUpTight",&vbf_maxMjj_j1_isPileUpTight,"vbf_maxMjj_j1_isPileUpTight/O");

   inputTree->Branch("vbf_maxMjj_j2_isPileUpLoose",&vbf_maxMjj_j2_isPileUpLoose,"vbf_maxMjj_j2_isPileUpLoose/O");
   inputTree->Branch("vbf_maxMjj_j2_isPileUpMedium",&vbf_maxMjj_j2_isPileUpMedium,"vbf_maxMjj_j2_isPileUpMedium/O");
   inputTree->Branch("vbf_maxMjj_j2_isPileUpTight",&vbf_maxMjj_j2_isPileUpTight,"vbf_maxMjj_j2_isPileUpTight/O");

   inputTree->Branch("vbf_maxMjj_jj_type",&vbf_maxMjj_jj_type,"vbf_maxMjj_jj_type/I");
   inputTree->Branch("vbf_maxMjj_n_excj",&vbf_maxMjj_n_excj,"vbf_maxMjj_n_excj/I");
   inputTree->Branch("vbf_maxMjj_n_exfj",&vbf_maxMjj_n_exfj,"vbf_maxMjj_n_exfj/I");

   inputTree->Branch("vbf_maxMjj_j1_bDiscriminatorSSVHE",&vbf_maxMjj_j1_bDiscriminatorSSVHE,"vbf_maxMjj_j1_bDiscriminatorSSVHE/F");
   inputTree->Branch("vbf_maxMjj_j1_bDiscriminatorTCHE",&vbf_maxMjj_j1_bDiscriminatorTCHE,"vbf_maxMjj_j1_bDiscriminatorTCHE/F");  
   inputTree->Branch("vbf_maxMjj_j1_bDiscriminatorCSV",&vbf_maxMjj_j1_bDiscriminatorCSV,"vbf_maxMjj_j1_bDiscriminatorCSV/F");
   inputTree->Branch("vbf_maxMjj_j1_bDiscriminatorSSVHP",&vbf_maxMjj_j1_bDiscriminatorSSVHP,"vbf_maxMjj_j1_bDiscriminatorSSVHP/F");
   inputTree->Branch("vbf_maxMjj_j1_bDiscriminatorTCHP",&vbf_maxMjj_j1_bDiscriminatorTCHP,"vbf_maxMjj_j1_bDiscriminatorTCHP/F");

   inputTree->Branch("vbf_maxMjj_j2_bDiscriminatorSSVHE",&vbf_maxMjj_j2_bDiscriminatorSSVHE,"vbf_maxMjj_j2_bDiscriminatorSSVHE/F");
   inputTree->Branch("vbf_maxMjj_j2_bDiscriminatorTCHE",&vbf_maxMjj_j2_bDiscriminatorTCHE,"vbf_maxMjj_j2_bDiscriminatorTCHE/F");  
   inputTree->Branch("vbf_maxMjj_j2_bDiscriminatorCSV",&vbf_maxMjj_j2_bDiscriminatorCSV,"vbf_maxMjj_j2_bDiscriminatorCSV/F");
   inputTree->Branch("vbf_maxMjj_j2_bDiscriminatorSSVHP",&vbf_maxMjj_j2_bDiscriminatorSSVHP,"vbf_maxMjj_j2_bDiscriminatorSSVHP/F");
   inputTree->Branch("vbf_maxMjj_j2_bDiscriminatorTCHP",&vbf_maxMjj_j2_bDiscriminatorTCHP,"vbf_maxMjj_j2_bDiscriminatorTCHP/F");

   inputTree->Branch("vbf_maxMjj_j1_ChargedHadronEnergy",&vbf_maxMjj_j1_ChargedHadronEnergy,"vbf_maxMjj_j1_ChargedHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_ChargedHadronEnergyFrac",&vbf_maxMjj_j1_ChargedHadronEnergyFrac,"vbf_maxMjj_j1_ChargedHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j1_NeutralHadronEnergy",&vbf_maxMjj_j1_NeutralHadronEnergy,"vbf_maxMjj_j1_NeutralHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_NeutralHadronEnergyFrac",&vbf_maxMjj_j1_NeutralHadronEnergyFrac,"vbf_maxMjj_j1_NeutralHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j1_ChargedEmEnergy",&vbf_maxMjj_j1_ChargedEmEnergy,"vbf_maxMjj_j1_ChargedEmEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_ChargedEmEnergyFrac",&vbf_maxMjj_j1_ChargedEmEnergyFrac,"vbf_maxMjj_j1_ChargedEmEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j1_ChargedMuEnergy",&vbf_maxMjj_j1_ChargedMuEnergy,"vbf_maxMjj_j1_ChargedMuEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_ChargedMuEnergyFrac",&vbf_maxMjj_j1_ChargedMuEnergyFrac,"vbf_maxMjj_j1_ChargedMuEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j1_NeutralEmEnergy",&vbf_maxMjj_j1_NeutralEmEnergy,"vbf_maxMjj_j1_NeutralEmEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_NeutralEmEnergyFrac",&vbf_maxMjj_j1_NeutralEmEnergyFrac,"vbf_maxMjj_j1_NeutralEmEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j1_ChargedMultiplicity",&vbf_maxMjj_j1_ChargedMultiplicity,"vbf_maxMjj_j1_ChargedMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j1_NeutralMultiplicity",&vbf_maxMjj_j1_NeutralMultiplicity,"vbf_maxMjj_j1_NeutralMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j1_MuonMultiplicity",&vbf_maxMjj_j1_MuonMultiplicity,"vbf_maxMjj_j1_MuonMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j1_PhotonEnergy",&vbf_maxMjj_j1_PhotonEnergy,"vbf_maxMjj_j1_PhotonEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_PhotonEnergyFraction",&vbf_maxMjj_j1_PhotonEnergyFraction,"vbf_maxMjj_j1_PhotonEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j1_ElectronEnergy",&vbf_maxMjj_j1_ElectronEnergy,"vbf_maxMjj_j1_ElectronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_ElectronEnergyFraction",&vbf_maxMjj_j1_ElectronEnergyFraction,"vbf_maxMjj_j1_ElectronEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j1_HFHadronEnergy",&vbf_maxMjj_j1_HFHadronEnergy,"vbf_maxMjj_j1_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_HFHadronEnergyFraction",&vbf_maxMjj_j1_HFHadronEnergyFraction,"vbf_maxMjj_j1_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j1_HFEMEnergy",&vbf_maxMjj_j1_HFEMEnergy,"vbf_maxMjj_j1_HFEMEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_HFHadronEnergy",&vbf_maxMjj_j1_HFHadronEnergy,"vbf_maxMjj_j1_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j1_HFHadronEnergyFraction",&vbf_maxMjj_j1_HFHadronEnergyFraction,"vbf_maxMjj_j1_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j1_ChargedHadronMultiplicity",&vbf_maxMjj_j1_ChargedHadronMultiplicity,"vbf_maxMjj_j1_ChargedHadronMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j1_NeutralHadronMultiplicity",&vbf_maxMjj_j1_NeutralHadronMultiplicity,"vbf_maxMjj_j1_NeutralHadronMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j1_PhotonMultiplicity",&vbf_maxMjj_j1_PhotonMultiplicity,"vbf_maxMjj_j1_PhotonMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j1_ElectronMultiplicity",&vbf_maxMjj_j1_ElectronMultiplicity,"vbf_maxMjj_j1_ElectronMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j1_HFHadronMultiplicity",&vbf_maxMjj_j1_HFHadronMultiplicity,"vbf_maxMjj_j1_HFHadronMultiplicity/F");

   inputTree->Branch("vbf_maxMjj_j2_ChargedHadronEnergy",&vbf_maxMjj_j2_ChargedHadronEnergy,"vbf_maxMjj_j2_ChargedHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_ChargedHadronEnergyFrac",&vbf_maxMjj_j2_ChargedHadronEnergyFrac,"vbf_maxMjj_j2_ChargedHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j2_NeutralHadronEnergy",&vbf_maxMjj_j2_NeutralHadronEnergy,"vbf_maxMjj_j2_NeutralHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_NeutralHadronEnergyFrac",&vbf_maxMjj_j2_NeutralHadronEnergyFrac,"vbf_maxMjj_j2_NeutralHadronEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j2_ChargedEmEnergy",&vbf_maxMjj_j2_ChargedEmEnergy,"vbf_maxMjj_j2_ChargedEmEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_ChargedEmEnergyFrac",&vbf_maxMjj_j2_ChargedEmEnergyFrac,"vbf_maxMjj_j2_ChargedEmEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j2_ChargedMuEnergy",&vbf_maxMjj_j2_ChargedMuEnergy,"vbf_maxMjj_j2_ChargedMuEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_ChargedMuEnergyFrac",&vbf_maxMjj_j2_ChargedMuEnergyFrac,"vbf_maxMjj_j2_ChargedMuEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j2_NeutralEmEnergy",&vbf_maxMjj_j2_NeutralEmEnergy,"vbf_maxMjj_j2_NeutralEmEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_NeutralEmEnergyFrac",&vbf_maxMjj_j2_NeutralEmEnergyFrac,"vbf_maxMjj_j2_NeutralEmEnergyFrac/F");
   inputTree->Branch("vbf_maxMjj_j2_ChargedMultiplicity",&vbf_maxMjj_j2_ChargedMultiplicity,"vbf_maxMjj_j2_ChargedMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j2_NeutralMultiplicity",&vbf_maxMjj_j2_NeutralMultiplicity,"vbf_maxMjj_j2_NeutralMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j2_MuonMultiplicity",&vbf_maxMjj_j2_MuonMultiplicity,"vbf_maxMjj_j2_MuonMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j2_PhotonEnergy",&vbf_maxMjj_j2_PhotonEnergy,"vbf_maxMjj_j2_PhotonEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_PhotonEnergyFraction",&vbf_maxMjj_j2_PhotonEnergyFraction,"vbf_maxMjj_j2_PhotonEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j2_ElectronEnergy",&vbf_maxMjj_j2_ElectronEnergy,"vbf_maxMjj_j2_ElectronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_ElectronEnergyFraction",&vbf_maxMjj_j2_ElectronEnergyFraction,"vbf_maxMjj_j2_ElectronEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j2_HFHadronEnergy",&vbf_maxMjj_j2_HFHadronEnergy,"vbf_maxMjj_j2_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_HFHadronEnergyFraction",&vbf_maxMjj_j2_HFHadronEnergyFraction,"vbf_maxMjj_j2_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j2_HFEMEnergy",&vbf_maxMjj_j2_HFEMEnergy,"vbf_maxMjj_j2_HFEMEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_HFHadronEnergy",&vbf_maxMjj_j2_HFHadronEnergy,"vbf_maxMjj_j2_HFHadronEnergy/F");
   inputTree->Branch("vbf_maxMjj_j2_HFHadronEnergyFraction",&vbf_maxMjj_j2_HFHadronEnergyFraction,"vbf_maxMjj_j2_HFHadronEnergyFraction/F");
   inputTree->Branch("vbf_maxMjj_j2_ChargedHadronMultiplicity",&vbf_maxMjj_j2_ChargedHadronMultiplicity,"vbf_maxMjj_j2_ChargedHadronMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j2_NeutralHadronMultiplicity",&vbf_maxMjj_j2_NeutralHadronMultiplicity,"vbf_maxMjj_j2_NeutralHadronMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j2_PhotonMultiplicity",&vbf_maxMjj_j2_PhotonMultiplicity,"vbf_maxMjj_j2_PhotonMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j2_ElectronMultiplicity",&vbf_maxMjj_j2_ElectronMultiplicity,"vbf_maxMjj_j2_ElectronMultiplicity/F");
   inputTree->Branch("vbf_maxMjj_j2_HFHadronMultiplicity",&vbf_maxMjj_j2_HFHadronMultiplicity,"vbf_maxMjj_j2_HFHadronMultiplicity/F");

}

void VBFMuonClass::InitializateVariables(){

 fit_mu_px=0 ,  fit_mu_py=0 ,  fit_mu_pz=0 ,  fit_mu_e=0 ;
 fit_nv_px=0 ,  fit_nv_py=0 ,  fit_nv_pz=0 ,  fit_nv_e=0 ;
 fit_subjet1_px=0 ,  fit_subjet1_py=0 ,  fit_subjet1_pz=0 ,  fit_subjet1_e =0 ;
 fit_subjet2_px=0 ,  fit_subjet2_py=0 ,  fit_subjet2_pz=0 ,  fit_subjet2_e =0 ;

 fit_lvj_m=0  , fit_lv_m =-999 , fit_j_m=-999, fit_chi2=-999 , fit_lvj_pt = -999, fit_lvj_eta = -999, fit_lvj_phi = -999 ;
 fit_lvj_phi=-999 ; fit_subjet1_m =-999; fit_subjet2_m=-999; 
 fit_NDF=-999, fit_status=-999 ;


 boosted_lvj_m =-999 , boosted_j_m=-999, boosted_subjet1_m=-999, boosted_subjet2_m=-999, boosted_lv_m=-999 , boosted_lvj_pt=-999;
 boosted_lvj_phi=-999, boosted_lvj_eta=-999, boosted_lvj_e=-999 ;

 boostedW_lvj_m =-999 , boostedW_j_m=-999, boostedW_subjet1_m=-999, boostedW_subjet2_m=-999, boostedW_lv_m=-999 , boostedW_lvj_pt=-999;
 boostedW_lvj_phi=-999, boostedW_lvj_eta=-999, boostedW_lvj_e=-999 ;

 boosted_wjj_ang_ha=-999, boosted_wjj_ang_hb=-999, boosted_wjj_ang_hs=-999, boosted_wjj_ang_phi=-999, boosted_wjj_ang_phia=-999, boosted_wjj_ang_phib=-999;

 vbf_maxpt_jj_e=-999,   vbf_maxpt_jj_pt=-999 ,   vbf_maxpt_jj_eta=-999 ,  vbf_maxpt_jj_phi=-999 , vbf_maxpt_jj_m=-999 ;
 vbf_maxpt_j1_e=-999 ,   vbf_maxpt_j1_pt=-999 ,   vbf_maxpt_j1_eta=-999 ,  vbf_maxpt_j1_phi=-999 , vbf_maxpt_j1_m=-999 ;

 vbf_maxpt_jj_e=-999 ,   vbf_maxpt_jj_pt=-999 ,   vbf_maxpt_jj_eta=-999 ,  vbf_maxpt_jj_phi=-999 , vbf_maxpt_jj_m=-999 ;
 vbf_maxpt_j1_e=-999 ,   vbf_maxpt_j1_pt=-999 ,   vbf_maxpt_j1_eta=-999 ,  vbf_maxpt_j1_phi=-999 , vbf_maxpt_j1_m=-999 ;
 vbf_maxpt_j2_e=-999 ,   vbf_maxpt_j2_pt=-999 ,   vbf_maxpt_j2_eta=-999 ,  vbf_maxpt_j2_phi=-999 , vbf_maxpt_j2_m=-999 ;
 vbf_maxpt_jj_deta=-999 ,vbf_maxpt_jj_dphi=-999;

 vbf_maxpt_j1_QGLikelihood=-999,  vbf_maxpt_j2_QGLikelihood=-999;

 vbf_maxpt_j1_isPileUpLoose=false  , vbf_maxpt_j2_isPileUpLoose=false  ;
 vbf_maxpt_j1_isPileUpMedium=false , vbf_maxpt_j2_isPileUpMedium=false ;
 vbf_maxpt_j1_isPileUpTight=false  , vbf_maxpt_j2_isPileUpTight=false  ;

 vbf_maxpt_jj_type=-999,   vbf_maxpt_n_excj=-999,   vbf_maxpt_n_exfj=-999;

 vbf_maxpt_j1_bDiscriminatorSSVHE=-999, vbf_maxpt_j1_bDiscriminatorTCHE=-999, vbf_maxpt_j1_bDiscriminatorCSV=-999;
 vbf_maxpt_j1_bDiscriminatorSSVHP=-999, vbf_maxpt_j1_bDiscriminatorTCHP=-999;
 vbf_maxpt_j2_bDiscriminatorSSVHE=-999, vbf_maxpt_j2_bDiscriminatorTCHE=-999, vbf_maxpt_j2_bDiscriminatorCSV=-999;
 vbf_maxpt_j2_bDiscriminatorSSVHP=-999, vbf_maxpt_j2_bDiscriminatorTCHP=-999;

 vbf_maxpt_j1_ChargedHadronEnergy=-999,  vbf_maxpt_j1_ChargedHadronEnergyFrac=-999, vbf_maxpt_j1_NeutralHadronEnergy=-999;
 vbf_maxpt_j1_NeutralHadronEnergyFrac=-999, vbf_maxpt_j1_ChargedEmEnergy=-999, vbf_maxpt_j1_ChargedEmEnergyFrac=-999, vbf_maxpt_j1_ChargedMuEnergy=-999;
 vbf_maxpt_j1_ChargedMuEnergyFrac=-999, vbf_maxpt_j1_NeutralEmEnergy=-999,vbf_maxpt_j1_NeutralEmEnergyFrac=-999, vbf_maxpt_j1_ChargedMultiplicity=-999;
 vbf_maxpt_j1_NeutralMultiplicity=-999, vbf_maxpt_j1_MuonMultiplicity=-999,vbf_maxpt_j1_PhotonEnergy=-999, vbf_maxpt_j1_PhotonEnergyFraction=-999;
 vbf_maxpt_j1_ElectronEnergy=-999, vbf_maxpt_j1_ElectronEnergyFraction=-999,vbf_maxpt_j1_MuonEnergy=-999, vbf_maxpt_j1_MuonEnergyFraction=-999;
 vbf_maxpt_j1_HFHadronEnergy=-999, vbf_maxpt_j1_HFHadronEnergyFraction=-999,vbf_maxpt_j1_HFEMEnergy=-999, vbf_maxpt_j1_HFEMEnergyFraction=-999;
 vbf_maxpt_j1_ChargedHadronMultiplicity=-999, vbf_maxpt_j1_NeutralHadronMultiplicity=-999, vbf_maxpt_j1_PhotonMultiplicity=-999;
 vbf_maxpt_j1_ElectronMultiplicity=-999,vbf_maxpt_j1_HFHadronMultiplicity=-999;

 vbf_maxpt_j2_ChargedHadronEnergy=-999,  vbf_maxpt_j2_ChargedHadronEnergyFrac=-999, vbf_maxpt_j2_NeutralHadronEnergy=-999;
 vbf_maxpt_j2_NeutralHadronEnergyFrac=-999, vbf_maxpt_j2_ChargedEmEnergy=-999, vbf_maxpt_j2_ChargedEmEnergyFrac=-999, vbf_maxpt_j2_ChargedMuEnergy=-999;
 vbf_maxpt_j2_ChargedMuEnergyFrac=-999, vbf_maxpt_j2_NeutralEmEnergy=-999,vbf_maxpt_j2_NeutralEmEnergyFrac=-999, vbf_maxpt_j2_ChargedMultiplicity=-999;
 vbf_maxpt_j2_NeutralMultiplicity=-999, vbf_maxpt_j2_MuonMultiplicity=-999,vbf_maxpt_j2_PhotonEnergy=-999, vbf_maxpt_j2_PhotonEnergyFraction=-999;
 vbf_maxpt_j2_ElectronEnergy=-999, vbf_maxpt_j2_ElectronEnergyFraction=-999,vbf_maxpt_j2_MuonEnergy=-999, vbf_maxpt_j2_MuonEnergyFraction=-999;
 vbf_maxpt_j2_HFHadronEnergy=-999, vbf_maxpt_j2_HFHadronEnergyFraction=-999,vbf_maxpt_j2_HFEMEnergy=-999, vbf_maxpt_j2_HFEMEnergyFraction=-999;
 vbf_maxpt_j2_ChargedHadronMultiplicity=-999, vbf_maxpt_j2_NeutralHadronMultiplicity=-999, vbf_maxpt_j2_PhotonMultiplicity=-999;
 vbf_maxpt_j2_ElectronMultiplicity=-999,vbf_maxpt_j2_HFHadronMultiplicity=-999;


 vbf_maxDeta_jj_e=-999,   vbf_maxDeta_jj_pt=-999 ,   vbf_maxDeta_jj_eta=-999 ,  vbf_maxDeta_jj_phi=-999 , vbf_maxDeta_jj_m=-999 ;
 vbf_maxDeta_j1_e=-999 ,   vbf_maxDeta_j1_pt=-999 ,   vbf_maxDeta_j1_eta=-999 ,  vbf_maxDeta_j1_phi=-999 , vbf_maxDeta_j1_m=-999 ;

 vbf_maxDeta_jj_e=-999 ,   vbf_maxDeta_jj_pt=-999 ,   vbf_maxDeta_jj_eta=-999 ,  vbf_maxDeta_jj_phi=-999 , vbf_maxDeta_jj_m=-999 ;
 vbf_maxDeta_j1_e=-999 ,   vbf_maxDeta_j1_pt=-999 ,   vbf_maxDeta_j1_eta=-999 ,  vbf_maxDeta_j1_phi=-999 , vbf_maxDeta_j1_m=-999 ;
 vbf_maxDeta_j2_e=-999 ,   vbf_maxDeta_j2_pt=-999 ,   vbf_maxDeta_j2_eta=-999 ,  vbf_maxDeta_j2_phi=-999 , vbf_maxDeta_j2_m=-999 ;
 vbf_maxDeta_jj_deta=-999 ,vbf_maxDeta_jj_dphi=-999;

 vbf_maxDeta_j1_QGLikelihood=-999,  vbf_maxDeta_j2_QGLikelihood=-999;

 vbf_maxDeta_j1_isPileUpLoose=false  , vbf_maxDeta_j2_isPileUpLoose=false  ;
 vbf_maxDeta_j1_isPileUpMedium=false , vbf_maxDeta_j2_isPileUpMedium=false ;
 vbf_maxDeta_j1_isPileUpTight=false  , vbf_maxDeta_j2_isPileUpTight=false  ;

 vbf_maxDeta_jj_type=-999,   vbf_maxDeta_n_excj=-999,   vbf_maxDeta_n_exfj=-999;

 vbf_maxDeta_j1_bDiscriminatorSSVHE=-999, vbf_maxDeta_j1_bDiscriminatorTCHE=-999, vbf_maxDeta_j1_bDiscriminatorCSV=-999;
 vbf_maxDeta_j1_bDiscriminatorSSVHP=-999, vbf_maxDeta_j1_bDiscriminatorTCHP=-999;
 vbf_maxDeta_j2_bDiscriminatorSSVHE=-999, vbf_maxDeta_j2_bDiscriminatorTCHE=-999, vbf_maxDeta_j2_bDiscriminatorCSV=-999;
 vbf_maxDeta_j2_bDiscriminatorSSVHP=-999, vbf_maxDeta_j2_bDiscriminatorTCHP=-999;

 vbf_maxDeta_j1_ChargedHadronEnergy=-999,  vbf_maxDeta_j1_ChargedHadronEnergyFrac=-999, vbf_maxDeta_j1_NeutralHadronEnergy=-999;
 vbf_maxDeta_j1_NeutralHadronEnergyFrac=-999, vbf_maxDeta_j1_ChargedEmEnergy=-999, vbf_maxDeta_j1_ChargedEmEnergyFrac=-999, vbf_maxDeta_j1_ChargedMuEnergy=-999;
 vbf_maxDeta_j1_ChargedMuEnergyFrac=-999, vbf_maxDeta_j1_NeutralEmEnergy=-999,vbf_maxDeta_j1_NeutralEmEnergyFrac=-999, vbf_maxDeta_j1_ChargedMultiplicity=-999;
 vbf_maxDeta_j1_NeutralMultiplicity=-999, vbf_maxDeta_j1_MuonMultiplicity=-999,vbf_maxDeta_j1_PhotonEnergy=-999, vbf_maxDeta_j1_PhotonEnergyFraction=-999;
 vbf_maxDeta_j1_ElectronEnergy=-999, vbf_maxDeta_j1_ElectronEnergyFraction=-999,vbf_maxDeta_j1_MuonEnergy=-999, vbf_maxDeta_j1_MuonEnergyFraction=-999;
 vbf_maxDeta_j1_HFHadronEnergy=-999, vbf_maxDeta_j1_HFHadronEnergyFraction=-999,vbf_maxDeta_j1_HFEMEnergy=-999, vbf_maxDeta_j1_HFEMEnergyFraction=-999;
 vbf_maxDeta_j1_ChargedHadronMultiplicity=-999, vbf_maxDeta_j1_NeutralHadronMultiplicity=-999, vbf_maxDeta_j1_PhotonMultiplicity=-999;
 vbf_maxDeta_j1_ElectronMultiplicity=-999,vbf_maxDeta_j1_HFHadronMultiplicity=-999;

 vbf_maxDeta_j2_ChargedHadronEnergy=-999,  vbf_maxDeta_j2_ChargedHadronEnergyFrac=-999, vbf_maxDeta_j2_NeutralHadronEnergy=-999;
 vbf_maxDeta_j2_NeutralHadronEnergyFrac=-999, vbf_maxDeta_j2_ChargedEmEnergy=-999, vbf_maxDeta_j2_ChargedEmEnergyFrac=-999, vbf_maxDeta_j2_ChargedMuEnergy=-999;
 vbf_maxDeta_j2_ChargedMuEnergyFrac=-999, vbf_maxDeta_j2_NeutralEmEnergy=-999,vbf_maxDeta_j2_NeutralEmEnergyFrac=-999, vbf_maxDeta_j2_ChargedMultiplicity=-999;
 vbf_maxDeta_j2_NeutralMultiplicity=-999, vbf_maxDeta_j2_MuonMultiplicity=-999,vbf_maxDeta_j2_PhotonEnergy=-999, vbf_maxDeta_j2_PhotonEnergyFraction=-999;
 vbf_maxDeta_j2_ElectronEnergy=-999, vbf_maxDeta_j2_ElectronEnergyFraction=-999,vbf_maxDeta_j2_MuonEnergy=-999, vbf_maxDeta_j2_MuonEnergyFraction=-999;
 vbf_maxDeta_j2_HFHadronEnergy=-999, vbf_maxDeta_j2_HFHadronEnergyFraction=-999,vbf_maxDeta_j2_HFEMEnergy=-999, vbf_maxDeta_j2_HFEMEnergyFraction=-999;
 vbf_maxDeta_j2_ChargedHadronMultiplicity=-999, vbf_maxDeta_j2_NeutralHadronMultiplicity=-999, vbf_maxDeta_j2_PhotonMultiplicity=-999;
 vbf_maxDeta_j2_ElectronMultiplicity=-999,vbf_maxDeta_j2_HFHadronMultiplicity=-999;


 vbf_maxMjj_jj_e=-999,   vbf_maxMjj_jj_pt=-999 ,   vbf_maxMjj_jj_eta=-999 ,  vbf_maxMjj_jj_phi=-999 , vbf_maxMjj_jj_m=-999 ;
 vbf_maxMjj_j1_e=-999 ,   vbf_maxMjj_j1_pt=-999 ,   vbf_maxMjj_j1_eta=-999 ,  vbf_maxMjj_j1_phi=-999 , vbf_maxMjj_j1_m=-999 ;

 vbf_maxMjj_jj_e=-999 ,   vbf_maxMjj_jj_pt=-999 ,   vbf_maxMjj_jj_eta=-999 ,  vbf_maxMjj_jj_phi=-999 , vbf_maxMjj_jj_m=-999 ;
 vbf_maxMjj_j1_e=-999 ,   vbf_maxMjj_j1_pt=-999 ,   vbf_maxMjj_j1_eta=-999 ,  vbf_maxMjj_j1_phi=-999 , vbf_maxMjj_j1_m=-999 ;
 vbf_maxMjj_j2_e=-999 ,   vbf_maxMjj_j2_pt=-999 ,   vbf_maxMjj_j2_eta=-999 ,  vbf_maxMjj_j2_phi=-999 , vbf_maxMjj_j2_m=-999 ;
 vbf_maxMjj_jj_deta=-999 ,vbf_maxMjj_jj_dphi=-999;

 vbf_maxMjj_j1_QGLikelihood=-999,  vbf_maxMjj_j2_QGLikelihood=-999;

 vbf_maxMjj_j1_isPileUpLoose =false , vbf_maxMjj_j2_isPileUpLoose =false  ;
 vbf_maxMjj_j1_isPileUpMedium=false , vbf_maxMjj_j2_isPileUpMedium=false ;
 vbf_maxMjj_j1_isPileUpTight =false  , vbf_maxMjj_j2_isPileUpTight=false  ;

 vbf_maxMjj_jj_type=-999,   vbf_maxMjj_n_excj=-999,   vbf_maxMjj_n_exfj=-999;


 vbf_maxMjj_j1_bDiscriminatorSSVHE=-999, vbf_maxMjj_j1_bDiscriminatorTCHE=-999, vbf_maxMjj_j1_bDiscriminatorCSV=-999;
 vbf_maxMjj_j1_bDiscriminatorSSVHP=-999, vbf_maxMjj_j1_bDiscriminatorTCHP=-999;
 vbf_maxMjj_j2_bDiscriminatorSSVHE=-999, vbf_maxMjj_j2_bDiscriminatorTCHE=-999, vbf_maxMjj_j2_bDiscriminatorCSV=-999;
 vbf_maxMjj_j2_bDiscriminatorSSVHP=-999, vbf_maxMjj_j2_bDiscriminatorTCHP=-999;

 vbf_maxMjj_j1_ChargedHadronEnergy=-999,  vbf_maxMjj_j1_ChargedHadronEnergyFrac=-999, vbf_maxMjj_j1_NeutralHadronEnergy=-999;
 vbf_maxMjj_j1_NeutralHadronEnergyFrac=-999, vbf_maxMjj_j1_ChargedEmEnergy=-999, vbf_maxMjj_j1_ChargedEmEnergyFrac=-999, vbf_maxMjj_j1_ChargedMuEnergy=-999;
 vbf_maxMjj_j1_ChargedMuEnergyFrac=-999, vbf_maxMjj_j1_NeutralEmEnergy=-999,vbf_maxMjj_j1_NeutralEmEnergyFrac=-999, vbf_maxMjj_j1_ChargedMultiplicity=-999;
 vbf_maxMjj_j1_NeutralMultiplicity=-999, vbf_maxMjj_j1_MuonMultiplicity=-999,vbf_maxMjj_j1_PhotonEnergy=-999, vbf_maxMjj_j1_PhotonEnergyFraction=-999;
 vbf_maxMjj_j1_ElectronEnergy=-999, vbf_maxMjj_j1_ElectronEnergyFraction=-999,vbf_maxMjj_j1_MuonEnergy=-999, vbf_maxMjj_j1_MuonEnergyFraction=-999;
 vbf_maxMjj_j1_HFHadronEnergy=-999, vbf_maxMjj_j1_HFHadronEnergyFraction=-999,vbf_maxMjj_j1_HFEMEnergy=-999, vbf_maxMjj_j1_HFEMEnergyFraction=-999;
 vbf_maxMjj_j1_ChargedHadronMultiplicity=-999, vbf_maxMjj_j1_NeutralHadronMultiplicity=-999, vbf_maxMjj_j1_PhotonMultiplicity=-999;
 vbf_maxMjj_j1_ElectronMultiplicity=-999,vbf_maxMjj_j1_HFHadronMultiplicity=-999;

 vbf_maxMjj_j2_ChargedHadronEnergy=-999,  vbf_maxMjj_j2_ChargedHadronEnergyFrac=-999, vbf_maxMjj_j2_NeutralHadronEnergy=-999;
 vbf_maxMjj_j2_NeutralHadronEnergyFrac=-999, vbf_maxMjj_j2_ChargedEmEnergy=-999, vbf_maxMjj_j2_ChargedEmEnergyFrac=-999, vbf_maxMjj_j2_ChargedMuEnergy=-999;
 vbf_maxMjj_j2_ChargedMuEnergyFrac=-999, vbf_maxMjj_j2_NeutralEmEnergy=-999,vbf_maxMjj_j2_NeutralEmEnergyFrac=-999, vbf_maxMjj_j2_ChargedMultiplicity=-999;
 vbf_maxMjj_j2_NeutralMultiplicity=-999, vbf_maxMjj_j2_MuonMultiplicity=-999,vbf_maxMjj_j2_PhotonEnergy=-999, vbf_maxMjj_j2_PhotonEnergyFraction=-999;
 vbf_maxMjj_j2_ElectronEnergy=-999, vbf_maxMjj_j2_ElectronEnergyFraction=-999,vbf_maxMjj_j2_MuonEnergy=-999, vbf_maxMjj_j2_MuonEnergyFraction=-999;
 vbf_maxMjj_j2_HFHadronEnergy=-999, vbf_maxMjj_j2_HFHadronEnergyFraction=-999,vbf_maxMjj_j2_HFEMEnergy=-999, vbf_maxMjj_j2_HFEMEnergyFraction=-999;
 vbf_maxMjj_j2_ChargedHadronMultiplicity=-999, vbf_maxMjj_j2_NeutralHadronMultiplicity=-999, vbf_maxMjj_j2_PhotonMultiplicity=-999;
 vbf_maxMjj_j2_ElectronMultiplicity=-999,vbf_maxMjj_j2_HFHadronMultiplicity=-999;



}

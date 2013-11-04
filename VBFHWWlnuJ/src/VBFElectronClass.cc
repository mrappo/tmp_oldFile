#include "VBFElectronClass.h"


VBFElectronClass::VBFElectronClass(){}

VBFElectronClass::VBFElectronClass(TTree* inputTree){

  if(inputTree == 0) {
                       TFile* f = new TFile("/data2/rgerosa/RD_Tree_v1/RD_el_STopS_Tbar_CMSSW532.root");
                       fTree = (TTree*) f -> Get("WJet");
  }
  else fTree = inputTree ;

}


VBFElectronClass::VBFElectronClass(TFile* inputFile, std::string inputTreeName){

  if(inputFile == 0) {
                       TFile* f = new TFile("/data2/rgerosa/RD_Tree_v1/RD_el_STopS_Tbar_CMSSW532.root");
                       fTree = (TTree*) f -> Get("WJet");
  }
  else fTree = (TTree*) inputFile -> Get(inputTreeName.c_str());

}


VBFElectronClass::~VBFElectronClass(){

  delete fTree ; 
}

TTree* VBFElectronClass::GetTree(){

  return fTree ;
}

void VBFElectronClass::SetReader( TTree* inputTree){

  if(inputTree == 0) return ; 

  fReader = new treeReader((TTree*)(fTree), false);
  
  SetBranchAddressAndStatus(fTree);

}

void VBFElectronClass::SetTree(TTree* inputTree){

  if(inputTree==0){  TFile* f = new TFile("/data2/rgerosa/RD_Tree_v1/RD_el_STopS_Tbar_CMSSW532.root");
                     fTree = (TTree*) f -> Get("WJet");
  }
  
  else fTree = inputTree ;

  fReader = new treeReader((TTree*)(fTree), false);

  SetBranchAddressAndStatus(fTree);

}

void VBFElectronClass::SetBranchAddressAndStatus ( TTree* inputTree){

  fTree->SetBranchStatus("MassV*",0);
  fTree->SetBranchStatus("Mass*j*",0);
  fTree->SetBranchStatus("cosTheta*",0);
  fTree->SetBranchStatus("cosJackson*",0);
  fTree->SetBranchStatus("colorCorr*", 0);
  fTree->SetBranchStatus("Photon*",0);
  fTree->SetBranchStatus("*Photon*",0);
  fTree->SetBranchStatus("W_Hb*",0);
  fTree->SetBranchStatus("W_Lepton*",0);
  fTree->SetBranchStatus("W_Met*",0);
  fTree->SetBranchStatus("W_tLepton*",0);
  fTree->SetBranchStatus("W_tParton*",0);
  fTree->SetBranchStatus("W_tMet*",0);
  fTree->SetBranchStatus("W_tbb*",0);

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
  fTree->SetBranchStatus("GroomedJet_AK7*",0);
  fTree->SetBranchStatus("GenGroomedJet_AK7*",0);
  fTree->SetBranchStatus("GroomedJet_AK5*",0);
  fTree->SetBranchStatus("GenGroomedJet_AK5*",0);
  fTree->SetBranchStatus("vbf*",0);
  fTree->SetBranchStatus("boostedW*",0);
  fTree->SetBranchStatus("ttH*",0);

  fTree->SetBranchStatus("Hadronic_W_*",0);
  fTree->SetBranchStatus("WHadposition",0);


}


void VBFElectronClass::SetNewBranches ( TTree* inputTree){


    inputTree->Branch("WHadposition",&WHadposition,"WHadposition/I"); 
    inputTree->Branch("numberJetBin","std::vector<int>",&numberJetBin); 

    inputTree->Branch("Hadronic_W_Jet_mass_uncorr",&Hadronic_W_Jet_mass_uncorr,"Hadronic_W_Jet_mass_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_mass_tr_uncorr",&Hadronic_W_Jet_mass_tr_uncorr,"Hadronic_W_Jet_mass_tr_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_mass_ft_uncorr",&Hadronic_W_Jet_mass_ft_uncorr,"Hadronic_W_Jet_mass_ft_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_mass_pr_uncorr",&Hadronic_W_Jet_mass_pr_uncorr,"Hadronic_W_Jet_mass_pr_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_massdrop_pr_uncorr",&Hadronic_W_Jet_massdrop_pr_uncorr,"Hadronic_W_Jet_massdrop_pr_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_tau2tau1",&Hadronic_W_Jet_tau2tau1,"Hadronic_W_Jet_tau2tau1/F");
    inputTree->Branch("Hadronic_W_Jet_tau1",&Hadronic_W_Jet_tau1,"Hadronic_W_Jet_tau1/F");
    inputTree->Branch("Hadronic_W_Jet_tau2",&Hadronic_W_Jet_tau2,"Hadronic_W_Jet_tau2/F");
    inputTree->Branch("Hadronic_W_Jet_tau3",&Hadronic_W_Jet_tau3,"Hadronic_W_Jet_tau3/F");
    inputTree->Branch("Hadronic_W_Jet_tau4",&Hadronic_W_Jet_tau4,"Hadronic_W_Jet_tau4/F");
    inputTree->Branch("Hadronic_W_Jet_pt",&Hadronic_W_Jet_pt,"Hadronic_W_Jet_pt/F");
    inputTree->Branch("Hadronic_W_Jet_eta",&Hadronic_W_Jet_eta,"Hadronic_W_Jet_eta/F");
    inputTree->Branch("Hadronic_W_Jet_phi",&Hadronic_W_Jet_phi,"Hadronic_W_Jet_phi/F");
    inputTree->Branch("Hadronic_W_Jet_e",&Hadronic_W_Jet_e,"Hadronic_W_Jet_e/F");
    inputTree->Branch("Hadronic_W_Jet_pt_tr_uncorr",&Hadronic_W_Jet_pt_tr_uncorr,"Hadronic_W_Jet_pt_tr_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_pt_tr",&Hadronic_W_Jet_pt_tr,"Hadronic_W_Jet_pt_tr/F");
    inputTree->Branch("Hadronic_W_Jet_eta_tr",&Hadronic_W_Jet_eta_tr,"Hadronic_W_Jet_eta_tr/F");
    inputTree->Branch("Hadronic_W_Jet_phi_tr",&Hadronic_W_Jet_phi_tr,"Hadronic_W_Jet_phi_tr/F");
    inputTree->Branch("Hadronic_W_Jet_e_tr",&Hadronic_W_Jet_e_tr,"Hadronic_W_Jet_e_tr/F");
    inputTree->Branch("Hadronic_W_Jet_pt_ft_uncorr",&Hadronic_W_Jet_pt_ft_uncorr,"Hadronic_W_Jet_pt_ft_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_pt_ft",&Hadronic_W_Jet_pt_ft,"Hadronic_W_Jet_pt_ft/F");
    inputTree->Branch("Hadronic_W_Jet_eta_ft",&Hadronic_W_Jet_eta_ft,"Hadronic_W_Jet_eta_ft/F");
    inputTree->Branch("Hadronic_W_Jet_phi_ft",&Hadronic_W_Jet_phi_ft,"Hadronic_W_Jet_phi_ft/F");
    inputTree->Branch("Hadronic_W_Jet_e_ft",&Hadronic_W_Jet_e_ft,"Hadronic_W_Jet_e_ft/F");
    inputTree->Branch("Hadronic_W_Jet_pt_pr_uncorr",&Hadronic_W_Jet_pt_pr_uncorr,"Hadronic_W_Jet_pt_pr_uncorr/F");
    inputTree->Branch("Hadronic_W_Jet_pt_pr",&Hadronic_W_Jet_pt_pr,"Hadronic_W_Jet_pt_pr/F");
    inputTree->Branch("Hadronic_W_Jet_eta_pr",&Hadronic_W_Jet_eta_pr,"Hadronic_W_Jet_eta_pr/F");
    inputTree->Branch("Hadronic_W_Jet_phi_pr",&Hadronic_W_Jet_phi_pr,"Hadronic_W_Jet_phi_pr/F");
    inputTree->Branch("Hadronic_W_Jet_e_pr",&Hadronic_W_Jet_e_pr,"Hadronic_W_Jet_e_pr/F");
    inputTree->Branch("Hadronic_W_Jet_prsubjet1_px",&Hadronic_W_Jet_prsubjet1_px,"Hadronic_W_Jet_prsubjet1_px/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet1_py",&Hadronic_W_Jet_prsubjet1_py,"Hadronic_W_Jet_prsubjet1_py/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet1_pz",&Hadronic_W_Jet_prsubjet1_pz,"Hadronic_W_Jet_prsubjet1_pz/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet1_e",&Hadronic_W_Jet_prsubjet1_e,"Hadronic_W_Jet_prsubjet1_e/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet2_px",&Hadronic_W_Jet_prsubjet2_px,"Hadronic_W_Jet_prsubjet2_px/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet2_py",&Hadronic_W_Jet_prsubjet2_py,"Hadronic_W_Jet_prsubjet2_py/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet2_pz",&Hadronic_W_Jet_prsubjet2_pz,"Hadronic_W_Jet_prsubjet2_pz/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet2_e",&Hadronic_W_Jet_prsubjet2_e,"Hadronic_W_Jet_prsubjet2_e/F");  
    inputTree->Branch("Hadronic_W_Jet_mass",&Hadronic_W_Jet_mass_pr,"Hadronic_W_Jet_mass/F"); 
    inputTree->Branch("Hadronic_W_Jet_mass_tr",&Hadronic_W_Jet_mass_pr,"Hadronic_W_Jet_mass_tr/F"); 
    inputTree->Branch("Hadronic_W_Jet_mass_ft",&Hadronic_W_Jet_mass_pr,"Hadronic_W_Jet_mass_ft/F");     
    inputTree->Branch("Hadronic_W_Jet_mass_pr",&Hadronic_W_Jet_mass_pr,"Hadronic_W_Jet_mass_pr/F");  
    inputTree->Branch("Hadronic_W_Jet_massdrop",&Hadronic_W_Jet_massdrop,"Hadronic_W_Jet_massdrop/F");
    inputTree->Branch("Hadronic_W_Jet_area",&Hadronic_W_Jet_area_pr,"Hadronic_W_Jet_area/F"); 
    inputTree->Branch("Hadronic_W_Jet_area_tr",&Hadronic_W_Jet_area_pr,"Hadronic_W_Jet_area_tr/F"); 
    inputTree->Branch("Hadronic_W_Jet_area_ft",&Hadronic_W_Jet_area_pr,"Hadronic_W_Jet_area_ft/F"); 
    inputTree->Branch("Hadronic_W_Jet_area_pr",&Hadronic_W_Jet_area_pr,"Hadronic_W_Jet_area_pr/F"); 
    inputTree->Branch("Hadronic_W_Jet_jetconsituents",&Hadronic_W_Jet_jetconsituents,"Hadronic_W_Jet_jetconsituents/F"); 
    inputTree->Branch("Hadronic_W_Jet_jetcharge",&Hadronic_W_Jet_jetcharge,"Hadronic_W_Jet_jetcharge/F");
    inputTree->Branch("Hadronic_W_Jet_rcores",&Hadronic_W_Jet_rcores,"Hadronic_W_Jet_rcores/F");  
    inputTree->Branch("Hadronic_W_Jet_ptcores",&Hadronic_W_Jet_ptcores,"Hadronic_W_Jet_ptcores/F");  
    inputTree->Branch("Hadronic_W_Jet_planarflow",&Hadronic_W_Jet_planarflow,"Hadronic_W_Jet_planarflow/F");  
    inputTree->Branch("Hadronic_W_Jet_qjetmass",&Hadronic_W_Jet_qjetmass,"Hadronic_W_Jet_qjetmass/F");
    inputTree->Branch("Hadronic_W_Jet_qjetmassdrop",&Hadronic_W_Jet_qjetmassdrop,"Hadronic_W_Jet_qjetmassdrop/F");
    inputTree->Branch("Hadronic_W_Jet_deltaR_ljet",&Hadronic_W_Jet_deltaR_ljet,"Hadronic_W_Jet_deltaR_ljet/F");
    inputTree->Branch("Hadronic_W_Jet_deltaphi_METjet",&Hadronic_W_Jet_deltaphi_METjet,"Hadronic_W_Jet_deltaphi_METjet/F");
    inputTree->Branch("Hadronic_W_Jet_deltaphi_Vca8jet",&Hadronic_W_Jet_deltaphi_Vca8jet,"Hadronic_W_Jet_deltaphi_Vca8jet/F");
    inputTree->Branch("Hadronic_W_Jet_rcores01",&Hadronic_W_Jet_rcores01,"Hadronic_W_Jet_rcores01/F");
    inputTree->Branch("Hadronic_W_Jet_rcores02",&Hadronic_W_Jet_rcores02,"Hadronic_W_Jet_rcores02/F");
    inputTree->Branch("Hadronic_W_Jet_rcores03",&Hadronic_W_Jet_rcores03,"Hadronic_W_Jet_rcores03/F");
    inputTree->Branch("Hadronic_W_Jet_rcores04",&Hadronic_W_Jet_rcores04,"Hadronic_W_Jet_rcores04/F");
    inputTree->Branch("Hadronic_W_Jet_rcores05",&Hadronic_W_Jet_rcores05,"Hadronic_W_Jet_rcores05/F");
    inputTree->Branch("Hadronic_W_Jet_rcores06",&Hadronic_W_Jet_rcores06,"Hadronic_W_Jet_rcores06/F");
    inputTree->Branch("Hadronic_W_Jet_rcores07",&Hadronic_W_Jet_rcores07,"Hadronic_W_Jet_rcores07/F"); 
    inputTree->Branch("Hadronic_W_Jet_rcores08",&Hadronic_W_Jet_rcores08,"Hadronic_W_Jet_rcores08/F");
    inputTree->Branch("Hadronic_W_Jet_rcores09",&Hadronic_W_Jet_rcores09,"Hadronic_W_Jet_rcores09/F");
    inputTree->Branch("Hadronic_W_Jet_rcores10",&Hadronic_W_Jet_rcores10,"Hadronic_W_Jet_rcores10/F");
    inputTree->Branch("Hadronic_W_Jet_rcores11",&Hadronic_W_Jet_rcores11,"Hadronic_W_Jet_rcores11/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores01",&Hadronic_W_Jet_ptcores01,"Hadronic_W_Jet_ptcores01/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores02",&Hadronic_W_Jet_ptcores02,"Hadronic_W_Jet_ptcores02/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores03",&Hadronic_W_Jet_ptcores03,"Hadronic_W_Jet_ptcores03/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores04",&Hadronic_W_Jet_ptcores04,"Hadronic_W_Jet_ptcores04/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores05",&Hadronic_W_Jet_ptcores05,"Hadronic_W_Jet_ptcores05/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores06",&Hadronic_W_Jet_ptcores06,"Hadronic_W_Jet_ptcores06/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores07",&Hadronic_W_Jet_ptcores07,"Hadronic_W_Jet_ptcores07/F"); 
    inputTree->Branch("Hadronic_W_Jet_ptcores08",&Hadronic_W_Jet_ptcores08,"Hadronic_W_Jet_ptcores08/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores09",&Hadronic_W_Jet_ptcores09,"Hadronic_W_Jet_ptcores09/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores10",&Hadronic_W_Jet_ptcores10,"Hadronic_W_Jet_ptcores10/F");
    inputTree->Branch("Hadronic_W_Jet_ptcores11",&Hadronic_W_Jet_ptcores11,"Hadronic_W_Jet_ptcores11/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow01",&Hadronic_W_Jet_planarflow01,"Hadronic_W_Jet_planarflow01/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow02",&Hadronic_W_Jet_planarflow02,"Hadronic_W_Jet_planarflow02/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow03",&Hadronic_W_Jet_planarflow03,"Hadronic_W_Jet_planarflow03/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow04",&Hadronic_W_Jet_planarflow04,"Hadronic_W_Jet_planarflow04/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow05",&Hadronic_W_Jet_planarflow05,"Hadronic_W_Jet_planarflow05/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow06",&Hadronic_W_Jet_planarflow06,"Hadronic_W_Jet_planarflow06/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow07",&Hadronic_W_Jet_planarflow07,"Hadronic_W_Jet_planarflow07/F"); 
    inputTree->Branch("Hadronic_W_Jet_planarflow08",&Hadronic_W_Jet_planarflow08,"Hadronic_W_Jet_planarflow08/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow09",&Hadronic_W_Jet_planarflow09,"Hadronic_W_Jet_planarflow09/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow10",&Hadronic_W_Jet_planarflow10,"Hadronic_W_Jet_planarflow10/F");
    inputTree->Branch("Hadronic_W_Jet_planarflow11",&Hadronic_W_Jet_planarflow11,"Hadronic_W_Jet_planarflow11/F");
    inputTree->Branch("Hadronic_W_Jet_mass_sensi_tr",&Hadronic_W_Jet_mass_sensi_tr,"Hadronic_W_Jet_mass_sensi_tr/F"); 
    inputTree->Branch("Hadronic_W_Jet_mass_sensi_ft",&Hadronic_W_Jet_mass_sensi_ft,"Hadronic_W_Jet_mass_sensi_ft/F");     
    inputTree->Branch("Hadronic_W_Jet_mass_sensi_pr",&Hadronic_W_Jet_mass_sensi_pr,"Hadronic_W_Jet_mass_sensi_pr/F"); 
    inputTree->Branch("Hadronic_W_Jet_qjetmassvolatility",&Hadronic_W_Jet_qjetmassvolatility,"Hadronic_W_Jet_qjetmassvolatility/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet1ptoverjetpt",&Hadronic_W_Jet_prsubjet1ptoverjetpt,"Hadronic_W_Jet_prsubjet1ptoverjetpt/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet2ptoverjetpt",&Hadronic_W_Jet_prsubjet2ptoverjetpt,"Hadronic_W_Jet_prsubjet2ptoverjetpt/F"); 
    inputTree->Branch("Hadronic_W_Jet_prsubjet1subjet2_deltaR",&Hadronic_W_Jet_prsubjet1subjet2_deltaR,"Hadronic_W_Jet_prsubjet1subjet2_deltaR/F"); 


    /// kinematic fit branches                                                                                                                                                             
    inputTree->Branch("fit_el_px_type0_met",&fit_el_px_type0_met,"fit_el_px_type0_met/F");
    inputTree->Branch("fit_el_py_type0_met",&fit_el_py_type0_met,"fit_el_py_type0_met/F");
    inputTree->Branch("fit_el_pz_type0_met",&fit_el_pz_type0_met,"fit_el_pz_type0_met/F");
    inputTree->Branch("fit_el_e_type0_met",&fit_el_e_type0_met,"fit_el_e_type0_met/F");
    inputTree->Branch("fit_nv_px_type0_met",&fit_el_px_type0_met,"fit_nv_px_type0_met/F");
    inputTree->Branch("fit_nv_py_type0_met",&fit_nv_py_type0_met,"fit_nv_py_type0_met/F");
    inputTree->Branch("fit_nv_pz_type0_met",&fit_nv_pz_type0_met,"fit_nv_pz_type0_met/F");
    inputTree->Branch("fit_nv_e_type0_met",&fit_nv_px_type0_met,"fit_nv_e_type0_met/F");
    inputTree->Branch("fit_subjet1_px_type0_met",&fit_subjet1_px_type0_met,"fit_subjet1_px_type0_met/F");
    inputTree->Branch("fit_subjet1_py_type0_met",&fit_subjet1_py_type0_met,"fit_subjet1_py_type0_met/F");
    inputTree->Branch("fit_subjet1_pz_type0_met",&fit_subjet1_pz_type0_met,"fit_subjet1_pz_type0_met/F");
    inputTree->Branch("fit_subjet1_e_type0_met",&fit_subjet1_e_type0_met,"fit_subjet1_e_type0_met/F");
    inputTree->Branch("fit_subjet2_px_type0_met",&fit_subjet2_px_type0_met,"fit_subjet2_px_type0_met/F");
    inputTree->Branch("fit_subjet2_py_type0_met",&fit_subjet2_py_type0_met,"fit_subjet2_py_type0_met/F");
    inputTree->Branch("fit_subjet2_pz_type0_met",&fit_subjet2_pz_type0_met,"fit_subjet2_pz_type0_met/F");
    inputTree->Branch("fit_subjet2_e_type0_met",&fit_subjet2_e_type0_met,"fit_subjet2_e_type0_met/F");
    inputTree->Branch("fit_subjet1_m_type0_met",&fit_subjet1_m_type0_met,"fit_subjet1_m_type0_met/F");
    inputTree->Branch("fit_subjet2_m_type0_met",&fit_subjet2_m_type0_met,"fit_subjet2_m_type0_met/F");
    inputTree->Branch("fit_lvj_m_type0_met",&fit_lvj_m_type0_met,"fit_lvj_m_type0_met/F");
    inputTree->Branch("fit_lv_m_type0_met",&fit_lv_m_type0_met,"fit_lv_m_type0_met/F");
    inputTree->Branch("fit_j_m_type0_met",&fit_j_m_type0_met,"fit_j_m_type0_met/F");
    inputTree->Branch("fit_lvj_pt_type0_met",&fit_lvj_pt_type0_met,"fit_lvj_pt_type0_met/F");
    inputTree->Branch("fit_lvj_eta_type0_met", &fit_lvj_eta_type0_met,"fit_lvj_eta_type0_met/F");
    inputTree->Branch("fit_lvj_phi_type0_met",&fit_lvj_phi_type0_met,"fit_lvj_phi_type0_met/F");
    inputTree->Branch("fit_lvj_e_type0_met",&fit_lvj_e_type0_met,"fit_lvj_e_type0_met/F");
    inputTree->Branch("fit_chi2_type0_met",&fit_chi2_type0_met,"fit_chi2_type0_met/F");
    inputTree->Branch("fit_NDF_type0_met",&fit_NDF_type0_met,"fit_NDF_type0_met/I");
    inputTree->Branch("fit_status_type0_met",&fit_status_type0_met,"fit_status_type0_met/I");

    inputTree->Branch("fit_el_px_type0",&fit_el_px_type0,"fit_el_px_type0/F");
    inputTree->Branch("fit_el_py_type0",&fit_el_py_type0,"fit_el_py_type0/F");
    inputTree->Branch("fit_el_pz_type0",&fit_el_pz_type0,"fit_el_pz_type0/F"); 
    inputTree->Branch("fit_el_e_type0",&fit_el_e_type0,"fit_el_e_type0/F");
    inputTree->Branch("fit_nv_px_type0",&fit_el_px_type0,"fit_nv_px_type0/F");
    inputTree->Branch("fit_nv_py_type0",&fit_nv_py_type0,"fit_nv_py_type0/F");
    inputTree->Branch("fit_nv_pz_type0",&fit_nv_pz_type0,"fit_nv_pz_type0/F");
    inputTree->Branch("fit_nv_e_type0",&fit_nv_px_type0,"fit_nv_e_type0/F");
    inputTree->Branch("fit_subjet1_px_type0",&fit_subjet1_px_type0,"fit_subjet1_px_type0/F");
    inputTree->Branch("fit_subjet1_py_type0",&fit_subjet1_py_type0,"fit_subjet1_py_type0/F");
    inputTree->Branch("fit_subjet1_pz_type0",&fit_subjet1_pz_type0,"fit_subjet1_pz_type0/F");
    inputTree->Branch("fit_subjet1_e_type0",&fit_subjet1_e_type0,"fit_subjet1_e_type0/F");
    inputTree->Branch("fit_subjet2_px_type0",&fit_subjet2_px_type0,"fit_subjet2_px_type0/F");
    inputTree->Branch("fit_subjet2_py_type0",&fit_subjet2_py_type0,"fit_subjet2_py_type0/F");
    inputTree->Branch("fit_subjet2_pz_type0",&fit_subjet2_pz_type0,"fit_subjet2_pz_type0/F");
    inputTree->Branch("fit_subjet2_e_type0",&fit_subjet2_e_type0,"fit_subjet2_e_type0/F");
    inputTree->Branch("fit_subjet1_m_type0",&fit_subjet1_m_type0,"fit_subjet1_m_type0/F");
    inputTree->Branch("fit_subjet2_m_type0",&fit_subjet2_m_type0,"fit_subjet2_m_type0/F");
    inputTree->Branch("fit_lvj_m_type0",&fit_lvj_m_type0,"fit_lvj_m_type0/F");
    inputTree->Branch("fit_lv_m_type0",&fit_lv_m_type0,"fit_lv_m_type0/F");
    inputTree->Branch("fit_j_m_type0",&fit_j_m_type0,"fit_j_m_type0/F");
    inputTree->Branch("fit_lvj_pt_type0",&fit_lvj_pt_type0,"fit_lvj_pt_type0/F");
    inputTree->Branch("fit_lvj_eta_type0", &fit_lvj_eta_type0,"fit_lvj_eta_type0/F");
    inputTree->Branch("fit_lvj_phi_type0",&fit_lvj_phi_type0,"fit_lvj_phi_type0/F");
    inputTree->Branch("fit_lvj_e_type0",&fit_lvj_e_type0,"fit_lvj_e_type0/F");
    inputTree->Branch("fit_chi2_type0",&fit_chi2_type0,"fit_chi2_type0/F");
    inputTree->Branch("fit_NDF_type0",&fit_NDF_type0,"fit_NDF_type0/I");
    inputTree->Branch("fit_status_type0",&fit_status_type0,"fit_status_type0/I");


    inputTree->Branch("fit_el_px_type2",&fit_el_px_type2,"fit_el_px_type2/F");
    inputTree->Branch("fit_el_py_type2",&fit_el_py_type2,"fit_el_py_type2/F");
    inputTree->Branch("fit_el_pz_type2",&fit_el_pz_type2,"fit_el_pz_type2/F"); 
    inputTree->Branch("fit_el_e_type2",&fit_el_e_type2,"fit_el_e_type2/F");
    inputTree->Branch("fit_nv_px_type2",&fit_el_px_type2,"fit_nv_px_type2/F");
    inputTree->Branch("fit_nv_py_type2",&fit_nv_py_type2,"fit_nv_py_type2/F");
    inputTree->Branch("fit_nv_pz_type2",&fit_nv_pz_type2,"fit_nv_pz_type2/F");
    inputTree->Branch("fit_nv_e_type2",&fit_nv_px_type2,"fit_nv_e_type2/F");
    inputTree->Branch("fit_subjet1_px_type2",&fit_subjet1_px_type2,"fit_subjet1_px_type2/F");
    inputTree->Branch("fit_subjet1_py_type2",&fit_subjet1_py_type2,"fit_subjet1_py_type2/F");
    inputTree->Branch("fit_subjet1_pz_type2",&fit_subjet1_pz_type2,"fit_subjet1_pz_type2/F");
    inputTree->Branch("fit_subjet1_e_type2",&fit_subjet1_e_type2,"fit_subjet1_e_type2/F");
    inputTree->Branch("fit_subjet2_px_type2",&fit_subjet2_px_type2,"fit_subjet2_px_type2/F");
    inputTree->Branch("fit_subjet2_py_type2",&fit_subjet2_py_type2,"fit_subjet2_py_type2/F");
    inputTree->Branch("fit_subjet2_pz_type2",&fit_subjet2_pz_type2,"fit_subjet2_pz_type2/F");
    inputTree->Branch("fit_subjet2_e_type2",&fit_subjet2_e_type2,"fit_subjet2_e_type2/F");
    inputTree->Branch("fit_subjet1_m_type2",&fit_subjet1_m_type2,"fit_subjet1_m_type2/F");
    inputTree->Branch("fit_subjet2_m_type2",&fit_subjet2_m_type2,"fit_subjet2_m_type2/F");
    inputTree->Branch("fit_lvj_m_type2",&fit_lvj_m_type2,"fit_lvj_m_type2/F");
    inputTree->Branch("fit_lv_m_type2",&fit_lv_m_type2,"fit_lv_m_type2/F");
    inputTree->Branch("fit_j_m_type2",&fit_j_m_type2,"fit_j_m_type2/F");
    inputTree->Branch("fit_lvj_pt_type2",&fit_lvj_pt_type2,"fit_lvj_pt_type2/F");
    inputTree->Branch("fit_lvj_eta_type2", &fit_lvj_eta_type2,"fit_lvj_eta_type2/F");
    inputTree->Branch("fit_lvj_phi_type2",&fit_lvj_phi_type2,"fit_lvj_phi_type2/F");
    inputTree->Branch("fit_lvj_e_type2",&fit_lvj_e_type2,"fit_lvj_e_type2/F");
    inputTree->Branch("fit_chi2_type2",&fit_chi2_type2,"fit_chi2_type2/F");
    inputTree->Branch("fit_NDF_type2",&fit_NDF_type2,"fit_NDF_type2/I");
    inputTree->Branch("fit_status_type2",&fit_status_type2,"fit_status_type2/I");

    inputTree->Branch("fit_el_px_type2_met",&fit_el_px_type2_met,"fit_el_px_type2_met/F");
    inputTree->Branch("fit_el_py_type2_met",&fit_el_py_type2_met,"fit_el_py_type2_met/F");
    inputTree->Branch("fit_el_pz_type2_met",&fit_el_pz_type2_met,"fit_el_pz_type2_met/F"); 
    inputTree->Branch("fit_el_e_type2_met",&fit_el_e_type2_met,"fit_el_e_type2_met/F");
    inputTree->Branch("fit_nv_px_type2_met",&fit_el_px_type2_met,"fit_nv_px_type2_met/F");
    inputTree->Branch("fit_nv_py_type2_met",&fit_nv_py_type2_met,"fit_nv_py_type2_met/F");
    inputTree->Branch("fit_nv_pz_type2_met",&fit_nv_pz_type2_met,"fit_nv_pz_type2_met/F");
    inputTree->Branch("fit_nv_e_type2_met",&fit_nv_px_type2_met,"fit_nv_e_type2_met/F");
    inputTree->Branch("fit_subjet1_px_type2_met",&fit_subjet1_px_type2_met,"fit_subjet1_px_type2_met/F");
    inputTree->Branch("fit_subjet1_py_type2_met",&fit_subjet1_py_type2_met,"fit_subjet1_py_type2_met/F");
    inputTree->Branch("fit_subjet1_pz_type2_met",&fit_subjet1_pz_type2_met,"fit_subjet1_pz_type2_met/F");
    inputTree->Branch("fit_subjet1_e_type2_met",&fit_subjet1_e_type2_met,"fit_subjet1_e_type2_met/F");
    inputTree->Branch("fit_subjet2_px_type2_met",&fit_subjet2_px_type2_met,"fit_subjet2_px_type2_met/F");
    inputTree->Branch("fit_subjet2_py_type2_met",&fit_subjet2_py_type2_met,"fit_subjet2_py_type2_met/F");
    inputTree->Branch("fit_subjet2_pz_type2_met",&fit_subjet2_pz_type2_met,"fit_subjet2_pz_type2_met/F");
    inputTree->Branch("fit_subjet2_e_type2_met",&fit_subjet2_e_type2_met,"fit_subjet2_e_type2_met/F");
    inputTree->Branch("fit_subjet1_m_type2_met",&fit_subjet1_m_type2_met,"fit_subjet1_m_type2_met/F");
    inputTree->Branch("fit_subjet2_m_type2_met",&fit_subjet2_m_type2_met,"fit_subjet2_m_type2_met/F");
    inputTree->Branch("fit_lvj_m_type2_met",&fit_lvj_m_type2_met,"fit_lvj_m_type2_met/F");
    inputTree->Branch("fit_lv_m_type2_met",&fit_lv_m_type2_met,"fit_lv_m_type2_met/F");
    inputTree->Branch("fit_j_m_type2_met",&fit_j_m_type2_met,"fit_j_m_type2_met/F");
    inputTree->Branch("fit_lvj_pt_type2_met",&fit_lvj_pt_type2_met,"fit_lvj_pt_type2_met/F");
    inputTree->Branch("fit_lvj_eta_type2_met", &fit_lvj_eta_type2_met,"fit_lvj_eta_type2_met/F");
    inputTree->Branch("fit_lvj_phi_type2_met",&fit_lvj_phi_type2_met,"fit_lvj_phi_type2_met/F");
    inputTree->Branch("fit_lvj_e_type2_met",&fit_lvj_e_type2_met,"fit_lvj_e_type2_met/F");
    inputTree->Branch("fit_chi2_type2_met",&fit_chi2_type2_met,"fit_chi2_type2_met/F");
    inputTree->Branch("fit_NDF_type2_met",&fit_NDF_type2_met,"fit_NDF_type2_met/I");
    inputTree->Branch("fit_status_type2_met",&fit_status_type2_met,"fit_status_type2_met/I");

    // new branch for pz of the neutrino                                                                                                                                                     
    inputTree->Branch("W_mass_type0_met",&W_mass_type0_met,"W_mass_type0_met/F");
    inputTree->Branch("W_mass_type2_met",&W_mass_type2_met,"W_mass_type2_met/F");

    inputTree->Branch("W_mass_type0",&W_mass_type0,"W_mass_type0/F");
    inputTree->Branch("W_mass_type2",&W_mass_type2,"W_mass_type2/F");

    inputTree->Branch("W_pz_type0_met",&W_pz_type0_met,"W_pz_type0_met/F");
    inputTree->Branch("W_pz_type2_met",&W_pz_type2_met,"W_pz_type2_met/F");

    inputTree->Branch("W_pz_type0",&W_pz_type0,"W_pz_type0/F");
    inputTree->Branch("W_pz_type2",&W_pz_type2,"W_pz_type2/F");

    inputTree->Branch("W_nu1_pz_type0_met",&W_nu1_pz_type0_met,"W_nu1_pz_type0_met/F");
    inputTree->Branch("W_nu1_pz_type2_met",&W_nu1_pz_type2_met,"W_nu1_pz_type2_met/F");

    inputTree->Branch("W_nu1_pz_type0",&W_nu1_pz_type0,"W_nu1_pz_type0/F");
    inputTree->Branch("W_nu1_pz_type2",&W_nu1_pz_type2,"W_nu1_pz_type2/F");

    inputTree->Branch("W_nu2_pz_type0_met",&W_nu2_pz_type0_met,"W_nu2_pz_type0_met/F");
    inputTree->Branch("W_nu2_pz_type2_met",&W_nu2_pz_type2_met,"W_nu2_pz_type2_met/F");

    inputTree->Branch("W_nu2_pz_type0",&W_nu2_pz_type0,"W_nu2_pz_type0/F");
    inputTree->Branch("W_nu2_pz_type2",&W_nu2_pz_type2,"W_nu2_pz_type2/F");
 
    //////////////////////////////
    inputTree->Branch("boosted_lvj_e_type0",&boosted_lvj_e_type0,"boosted_lvj_e_type0/F");
    inputTree->Branch("boosted_lvj_pt_type0",&boosted_lvj_pt_type0,"boosted_lvj_pt_type0/F");
    inputTree->Branch("boosted_lvj_eta_type0",&boosted_lvj_eta_type0,"boosted_lvj_eta_type0/F");
    inputTree->Branch("boosted_lvj_phi_type0",&boosted_lvj_phi_type0,"boosted_lvj_phi_type0/F");
    inputTree->Branch("boosted_lvj_m_type0",&boosted_lvj_m_type0,"boosted_lvj_m_type0/F");
    inputTree->Branch("boosted_j_m_type0",&boosted_j_m_type0,"boosted_j_m_type0/F");
    inputTree->Branch("boosted_subjet1_m_type0",&boosted_subjet1_m_type0,"boosted_subjet1_m_type0/F");
    inputTree->Branch("boosted_subjet2_m_type0",&boosted_subjet2_m_type0,"boosted_subjet2_m_type0/F");
    inputTree->Branch("boosted_lv_m_type0",&boosted_lv_m_type0,"boosted_lv_m_type0/F");

    inputTree->Branch("boostedW_lvj_e_type0",&boostedW_lvj_e_type0,"boostedW_lvj_e_type0/F");
    inputTree->Branch("boostedW_lvj_pt_type0",&boostedW_lvj_pt_type0,"boostedW_lvj_pt_type0/F");
    inputTree->Branch("boostedW_lvj_eta_type0",&boostedW_lvj_eta_type0,"boostedW_lvj_eta_type0/F");
    inputTree->Branch("boostedW_lvj_phi_type0",&boostedW_lvj_phi_type0,"boostedW_lvj_phi_type0/F");
    inputTree->Branch("boostedW_lvj_m_type0",&boostedW_lvj_m_type0,"boostedW_lvj_m_type0/F");
    inputTree->Branch("boostedW_j_m_type0",&boostedW_j_m_type0,"boostedW_j_m_type0/F");
    inputTree->Branch("boostedW_subjet1_m_type0",&boostedW_subjet1_m_type0,"boostedW_subjet1_m_type0/F");
    inputTree->Branch("boostedW_subjet2_m_type0",&boostedW_subjet2_m_type0,"boostedW_subjet2_m_type0/F");
    inputTree->Branch("boostedW_lv_m_type0",&boostedW_lv_m_type0,"boostedW_lv_m_type0/F");

    inputTree->Branch("boosted_wjj_ang_ha_type0",&boosted_wjj_ang_ha_type0,"boosted_wjj_ang_ha_type0/F");
    inputTree->Branch("boosted_wjj_ang_hb_type0",&boosted_wjj_ang_ha_type0,"boosted_wjj_ang_ha_type0/F");
    inputTree->Branch("boosted_wjj_ang_hs_type0",&boosted_wjj_ang_hs_type0,"boosted_wjj_ang_hs_type0/F");
    inputTree->Branch("boosted_wjj_ang_phi_type0",&boosted_wjj_ang_phi_type0,"boosted_wjj_ang_phi_type0/F");
    inputTree->Branch("boosted_wjj_ang_phia_type0",&boosted_wjj_ang_phia_type0,"boosted_wjj_ang_phia_type0/F");
    inputTree->Branch("boosted_wjj_ang_phib_type0",&boosted_wjj_ang_phib_type0,"boosted_wjj_ang_phib_type0/F");

    inputTree->Branch("boosted_lvj_e_type2",&boosted_lvj_e_type2,"boosted_lvj_e_type2/F");
    inputTree->Branch("boosted_lvj_pt_type2",&boosted_lvj_pt_type2,"boosted_lvj_pt_type2/F");
    inputTree->Branch("boosted_lvj_eta_type2",&boosted_lvj_eta_type2,"boosted_lvj_eta_type2/F");
    inputTree->Branch("boosted_lvj_phi_type2",&boosted_lvj_phi_type2,"boosted_lvj_phi_type2/F");
    inputTree->Branch("boosted_lvj_m_type2",&boosted_lvj_m_type2,"boosted_lvj_m_type2/F");
    inputTree->Branch("boosted_j_m_type2",&boosted_j_m_type2,"boosted_j_m_type2/F");
    inputTree->Branch("boosted_subjet1_m_type2",&boosted_subjet1_m_type2,"boosted_subjet1_m_type2/F");
    inputTree->Branch("boosted_subjet2_m_type2",&boosted_subjet2_m_type2,"boosted_subjet2_m_type2/F");
    inputTree->Branch("boosted_lv_m_type2",&boosted_lv_m_type2,"boosted_lv_m_type2/F");

    inputTree->Branch("boostedW_lvj_e_type2",&boostedW_lvj_e_type2,"boostedW_lvj_e_type2/F");
    inputTree->Branch("boostedW_lvj_pt_type2",&boostedW_lvj_pt_type2,"boostedW_lvj_pt_type2/F");
    inputTree->Branch("boostedW_lvj_eta_type2",&boostedW_lvj_eta_type2,"boostedW_lvj_eta_type2/F");
    inputTree->Branch("boostedW_lvj_phi_type2",&boostedW_lvj_phi_type2,"boostedW_lvj_phi_type2/F");
    inputTree->Branch("boostedW_lvj_m_type2",&boostedW_lvj_m_type2,"boostedW_lvj_m_type2/F");
    inputTree->Branch("boostedW_j_m_type2",&boostedW_j_m_type2,"boostedW_j_m_type2/F");
    inputTree->Branch("boostedW_subjet1_m_type2",&boostedW_subjet1_m_type2,"boostedW_subjet1_m_type2/F");
    inputTree->Branch("boostedW_subjet2_m_type2",&boostedW_subjet2_m_type2,"boostedW_subjet2_m_type2/F");
    inputTree->Branch("boostedW_lv_m_type2",&boostedW_lv_m_type2,"boostedW_lv_m_type2/F");


    inputTree->Branch("boosted_wjj_ang_ha_type2",&boosted_wjj_ang_ha_type2,"boosted_wjj_ang_ha_type2/F");
    inputTree->Branch("boosted_wjj_ang_hb_type2",&boosted_wjj_ang_ha_type2,"boosted_wjj_ang_ha_type2/F");
    inputTree->Branch("boosted_wjj_ang_hs_type2",&boosted_wjj_ang_hs_type2,"boosted_wjj_ang_hs_type2/F");
    inputTree->Branch("boosted_wjj_ang_phi_type2",&boosted_wjj_ang_phi_type2,"boosted_wjj_ang_phi_type2/F");
    inputTree->Branch("boosted_wjj_ang_phia_type2",&boosted_wjj_ang_phia_type2,"boosted_wjj_ang_phia_type2/F");
    inputTree->Branch("boosted_wjj_ang_phib_type2",&boosted_wjj_ang_phib_type2,"boosted_wjj_ang_phib_type2/F");

    inputTree->Branch("boosted_lvj_e_type0_met",&boosted_lvj_e_type0_met,"boosted_lvj_e_type0_met/F");
    inputTree->Branch("boosted_lvj_pt_type0_met",&boosted_lvj_pt_type0_met,"boosted_lvj_pt_type0_met/F");
    inputTree->Branch("boosted_lvj_eta_type0_met",&boosted_lvj_eta_type0_met,"boosted_lvj_eta_type0_met/F");
    inputTree->Branch("boosted_lvj_phi_type0_met",&boosted_lvj_phi_type0_met,"boosted_lvj_phi_type0_met/F");
    inputTree->Branch("boosted_lvj_m_type0_met",&boosted_lvj_m_type0_met,"boosted_lvj_m_type0_met/F");
    inputTree->Branch("boosted_j_m_type0_met",&boosted_j_m_type0_met,"boosted_j_m_type0_met/F");
    inputTree->Branch("boosted_subjet1_m_type0_met",&boosted_subjet1_m_type0_met,"boosted_subjet1_m_type0_met/F");
    inputTree->Branch("boosted_subjet2_m_type0_met",&boosted_subjet2_m_type0_met,"boosted_subjet2_m_type0_met/F");
    inputTree->Branch("boosted_lv_m_type0_met",&boosted_lv_m_type0_met,"boosted_lv_m_type0_met/F");

    inputTree->Branch("boostedW_lvj_e_type0_met",&boostedW_lvj_e_type0_met,"boostedW_lvj_e_type0_met/F");
    inputTree->Branch("boostedW_lvj_pt_type0_met",&boostedW_lvj_pt_type0_met,"boostedW_lvj_pt_type0_met/F");
    inputTree->Branch("boostedW_lvj_eta_type0_met",&boostedW_lvj_eta_type0_met,"boostedW_lvj_eta_type0_met/F");
    inputTree->Branch("boostedW_lvj_phi_type0_met",&boostedW_lvj_phi_type0_met,"boostedW_lvj_phi_type0_met/F");
    inputTree->Branch("boostedW_lvj_m_type0_met",&boostedW_lvj_m_type0_met,"boostedW_lvj_m_type0_met/F");
    inputTree->Branch("boostedW_j_m_type0_met",&boostedW_j_m_type0_met,"boostedW_j_m_type0_met/F");
    inputTree->Branch("boostedW_subjet1_m_type0_met",&boostedW_subjet1_m_type0_met,"boostedW_subjet1_m_type0_met/F");
    inputTree->Branch("boostedW_subjet2_m_type0_met",&boostedW_subjet2_m_type0_met,"boostedW_subjet2_m_type0_met/F");
    inputTree->Branch("boostedW_lv_m_type0_met",&boostedW_lv_m_type0_met,"boostedW_lv_m_type0_met/F");

    inputTree->Branch("boosted_wjj_ang_ha_type0_met",&boosted_wjj_ang_ha_type0_met,"boosted_wjj_ang_ha_type0_met/F");
    inputTree->Branch("boosted_wjj_ang_hb_type0_met",&boosted_wjj_ang_ha_type0_met,"boosted_wjj_ang_ha_type0_met/F");
    inputTree->Branch("boosted_wjj_ang_hs_type0_met",&boosted_wjj_ang_hs_type0_met,"boosted_wjj_ang_hs_type0_met/F");
    inputTree->Branch("boosted_wjj_ang_phi_type0_met",&boosted_wjj_ang_phi_type0_met,"boosted_wjj_ang_phi_type0_met/F");
    inputTree->Branch("boosted_wjj_ang_phia_type0_met",&boosted_wjj_ang_phia_type0_met,"boosted_wjj_ang_phia_type0_met/F");
    inputTree->Branch("boosted_wjj_ang_phib_type0_met",&boosted_wjj_ang_phib_type0_met,"boosted_wjj_ang_phib_type0_met/F");

    inputTree->Branch("boosted_lvj_e_type2_met",&boosted_lvj_e_type2_met,"boosted_lvj_e_type2_met/F");
    inputTree->Branch("boosted_lvj_pt_type2_met",&boosted_lvj_pt_type2_met,"boosted_lvj_pt_type2_met/F");
    inputTree->Branch("boosted_lvj_eta_type2_met",&boosted_lvj_eta_type2_met,"boosted_lvj_eta_type2_met/F");
    inputTree->Branch("boosted_lvj_phi_type2_met",&boosted_lvj_phi_type2_met,"boosted_lvj_phi_type2_met/F");
    inputTree->Branch("boosted_lvj_m_type2_met",&boosted_lvj_m_type2_met,"boosted_lvj_m_type2_met/F");
    inputTree->Branch("boosted_j_m_type2_met",&boosted_j_m_type2_met,"boosted_j_m_type2_met/F");
    inputTree->Branch("boosted_subjet1_m_type2_met",&boosted_subjet1_m_type2_met,"boosted_subjet1_m_type2_met/F");
    inputTree->Branch("boosted_subjet2_m_type2_met",&boosted_subjet2_m_type2_met,"boosted_subjet2_m_type2_met/F");
    inputTree->Branch("boosted_lv_m_type2_met",&boosted_lv_m_type2_met,"boosted_lv_m_type2_met/F");

    inputTree->Branch("boostedW_lvj_e_type2_met",&boostedW_lvj_e_type2_met,"boostedW_lvj_e_type2_met/F");
    inputTree->Branch("boostedW_lvj_pt_type2_met",&boostedW_lvj_pt_type2_met,"boostedW_lvj_pt_type2_met/F");
    inputTree->Branch("boostedW_lvj_eta_type2_met",&boostedW_lvj_eta_type2_met,"boostedW_lvj_eta_type2_met/F");
    inputTree->Branch("boostedW_lvj_phi_type2_met",&boostedW_lvj_phi_type2_met,"boostedW_lvj_phi_type2_met/F");
    inputTree->Branch("boostedW_lvj_m_type2_met",&boostedW_lvj_m_type2_met,"boostedW_lvj_m_type2_met/F");
    inputTree->Branch("boostedW_j_m_type2_met",&boostedW_j_m_type2_met,"boostedW_j_m_type2_met/F");
    inputTree->Branch("boostedW_subjet1_m_type2_met",&boostedW_subjet1_m_type2_met,"boostedW_subjet1_m_type2_met/F");
    inputTree->Branch("boostedW_subjet2_m_type2_met",&boostedW_subjet2_m_type2_met,"boostedW_subjet2_m_type2_met/F");
    inputTree->Branch("boostedW_lv_m_type2_met",&boostedW_lv_m_type2_met,"boostedW_lv_m_type2_met/F");

    inputTree->Branch("boosted_wjj_ang_ha_type2_met",&boosted_wjj_ang_ha_type2_met,"boosted_wjj_ang_ha_type2_met/F");
    inputTree->Branch("boosted_wjj_ang_hb_type2_met",&boosted_wjj_ang_ha_type2_met,"boosted_wjj_ang_ha_type2_met/F");
    inputTree->Branch("boosted_wjj_ang_hs_type2_met",&boosted_wjj_ang_hs_type2_met,"boosted_wjj_ang_hs_type2_met/F");
    inputTree->Branch("boosted_wjj_ang_phi_type2_met",&boosted_wjj_ang_phi_type2_met,"boosted_wjj_ang_phi_type2_met/F");
    inputTree->Branch("boosted_wjj_ang_phia_type2_met",&boosted_wjj_ang_phia_type2_met,"boosted_wjj_ang_phia_type2_met/F");
    inputTree->Branch("boosted_wjj_ang_phib_type2_met",&boosted_wjj_ang_phib_type2_met,"boosted_wjj_ang_phib_type2_met/F");

    ///////// max pt
 
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

void VBFElectronClass::InitializateVariables(){



    WHadposition=-999; 
    numberJetBin.clear();
 
    Hadronic_W_Jet_mass_uncorr = -999;   
    Hadronic_W_Jet_mass_tr_uncorr = -999;  
    Hadronic_W_Jet_mass_ft_uncorr = -999;   
    Hadronic_W_Jet_mass_pr_uncorr = -999;   
    Hadronic_W_Jet_massdrop_pr_uncorr = -999;   
    Hadronic_W_Jet_tau2tau1 = -999;  
    Hadronic_W_Jet_tau1 = -999;  
    Hadronic_W_Jet_tau2 = -999;   
    Hadronic_W_Jet_tau3 = -999;   
    Hadronic_W_Jet_tau4 = -999;   
    Hadronic_W_Jet_pt = -999;   
    Hadronic_W_Jet_eta = -999;   
    Hadronic_W_Jet_phi = -999;   
    Hadronic_W_Jet_e = -999;   
    Hadronic_W_Jet_pt_tr_uncorr = -999;   
    Hadronic_W_Jet_pt_tr = -999;   
    Hadronic_W_Jet_eta_tr = -999;   
    Hadronic_W_Jet_phi_tr = -999;   
    Hadronic_W_Jet_e_tr = -999;  
    Hadronic_W_Jet_pt_ft_uncorr = -999;   
    Hadronic_W_Jet_pt_ft = -999;   
    Hadronic_W_Jet_eta_ft = -999;   
    Hadronic_W_Jet_phi_ft = -999;   
    Hadronic_W_Jet_e_ft = -999;   
    Hadronic_W_Jet_pt_pr_uncorr = -999;   
    Hadronic_W_Jet_pt_pr = -999;   
    Hadronic_W_Jet_eta_pr = -999;   
    Hadronic_W_Jet_phi_pr = -999;   
    Hadronic_W_Jet_e_pr = -999;   
    Hadronic_W_Jet_prsubjet1_px = -999;  
    Hadronic_W_Jet_prsubjet1_py = -999;  
    Hadronic_W_Jet_prsubjet1_pz = -999;   
    Hadronic_W_Jet_prsubjet1_e = -999;   
    Hadronic_W_Jet_prsubjet2_px = -999;   
    Hadronic_W_Jet_prsubjet2_py = -999;   
    Hadronic_W_Jet_prsubjet2_pz = -999;   
    Hadronic_W_Jet_prsubjet2_e = -999;   
    Hadronic_W_Jet_mass = -999;  
    Hadronic_W_Jet_mass_tr = -999;   
    Hadronic_W_Jet_mass_ft = -999;  
    Hadronic_W_Jet_mass_pr = -999;   
    Hadronic_W_Jet_massdrop = -999;   
    Hadronic_W_Jet_area = -999;   
    Hadronic_W_Jet_area_tr = -999;   
    Hadronic_W_Jet_area_ft = -999;   
    Hadronic_W_Jet_area_pr = -999;   
    Hadronic_W_Jet_jetconsituents = -999;   
    Hadronic_W_Jet_jetcharge = -999;   
    Hadronic_W_Jet_rcores = -999;   
    Hadronic_W_Jet_ptcores = -999;   
    Hadronic_W_Jet_planarflow = -999;   
    Hadronic_W_Jet_qjetmass = -999;   
    Hadronic_W_Jet_qjetmassdrop = -999;   
    Hadronic_W_Jet_deltaR_ljet = -999;   
    Hadronic_W_Jet_deltaphi_METjet = -999;   
    Hadronic_W_Jet_deltaphi_Vca8jet = -999;   
    Hadronic_W_Jet_rcores01 = -999;  
    Hadronic_W_Jet_rcores02 = -999;  
    Hadronic_W_Jet_rcores03 = -999;   
    Hadronic_W_Jet_rcores04 = -999;   
    Hadronic_W_Jet_rcores05 = -999;   
    Hadronic_W_Jet_rcores06 = -999;   
    Hadronic_W_Jet_rcores07 = -999;  
    Hadronic_W_Jet_rcores08 = -999;   
    Hadronic_W_Jet_rcores09 = -999;   
    Hadronic_W_Jet_rcores10 = -999;  
    Hadronic_W_Jet_rcores11 = -999;   
    Hadronic_W_Jet_ptcores01 = -999;  
    Hadronic_W_Jet_ptcores02 = -999;   
    Hadronic_W_Jet_ptcores03 = -999;  
    Hadronic_W_Jet_ptcores04 = -999;   
    Hadronic_W_Jet_ptcores05 = -999;   
    Hadronic_W_Jet_ptcores06 = -999;   
    Hadronic_W_Jet_ptcores07 = -999;   
    Hadronic_W_Jet_ptcores08 = -999;   
    Hadronic_W_Jet_ptcores09 = -999;   
    Hadronic_W_Jet_ptcores10 = -999;   
    Hadronic_W_Jet_ptcores11 = -999;  
    Hadronic_W_Jet_planarflow01 = -999;  
    Hadronic_W_Jet_planarflow02 = -999;   
    Hadronic_W_Jet_planarflow03 = -999;   
    Hadronic_W_Jet_planarflow04 = -999;   
    Hadronic_W_Jet_planarflow05 = -999;  
    Hadronic_W_Jet_planarflow06 = -999;  
    Hadronic_W_Jet_planarflow07 = -999;   
    Hadronic_W_Jet_planarflow08 = -999;   
    Hadronic_W_Jet_planarflow09 = -999;  
    Hadronic_W_Jet_planarflow10 = -999;   
    Hadronic_W_Jet_planarflow11 = -999;   
    Hadronic_W_Jet_mass_sensi_tr = -999;   
    Hadronic_W_Jet_mass_sensi_ft = -999;  
    Hadronic_W_Jet_mass_sensi_pr = -999;   
    Hadronic_W_Jet_qjetmassvolatility = -999;  
    Hadronic_W_Jet_prsubjet1ptoverjetpt = -999;   
    Hadronic_W_Jet_prsubjet2ptoverjetpt = -999; 
    Hadronic_W_Jet_prsubjet1subjet2_deltaR = -999;


    // new branch for pz of the neutrino                                                                                                                                                     
    W_mass_type0_met = 0. , W_pz_type0_met = -9999.,  W_nu1_pz_type0_met = -9999., W_nu2_pz_type0_met = -9999.;
    W_mass_type2_met = 0. , W_pz_type2_met = -9999.,  W_nu1_pz_type2_met = -9999., W_nu2_pz_type2_met = -9999.;

    W_mass_type0 = 0. , W_pz_type0 = -9999.,  W_nu1_pz_type0 = -9999., W_nu2_pz_type0 = -9999.;
    W_mass_type2 = 0. , W_pz_type2 = -9999.,  W_nu1_pz_type2 = -9999., W_nu2_pz_type2 = -9999.;

    fit_el_px_type0=-999 ,  fit_el_py_type0=-999 ,  fit_el_pz_type0=-999 ,  fit_el_e_type0=-999 ;
    fit_nv_px_type0=-999 ,  fit_nv_py_type0=-999 ,  fit_nv_pz_type0=-999 ,  fit_nv_e_type0=-999 ;
    fit_subjet1_px_type0=-999 ,  fit_subjet1_py_type0=-999 ,  fit_subjet1_pz_type0=-999 ,  fit_subjet1_e_type0 =-999 ;
    fit_subjet2_px_type0=-999 ,  fit_subjet2_py_type0=-999 ,  fit_subjet2_pz_type0=-999 ,  fit_subjet2_e_type0 =-999 ;
    fit_lvj_m_type0= -999  , fit_lv_m_type0 =-999 , fit_j_m_type0=-999, fit_chi2_type0=-999 , fit_lvj_pt_type0 = -999, fit_lvj_eta_type0 = -999, fit_lvj_phi_type0 = -999 ;
    fit_lvj_phi_type0=-999 ; fit_subjet1_m_type0 =-999; fit_subjet2_m_type0=-999; 
    fit_NDF_type0=-999, fit_status_type0=-999 ;

    fit_el_px_type2=-999 ,  fit_el_py_type2=-999 ,  fit_el_pz_type2=-999 ,  fit_el_e_type2=-999 ;
    fit_nv_px_type2=-999 ,  fit_nv_py_type2=-999 ,  fit_nv_pz_type2=-999 ,  fit_nv_e_type2=-999 ;
    fit_subjet1_px_type2=-999 ,  fit_subjet1_py_type2=-999 ,  fit_subjet1_pz_type2=-999 ,  fit_subjet1_e_type2 =-999 ;
    fit_subjet2_px_type2=-999 ,  fit_subjet2_py_type2=-999 ,  fit_subjet2_pz_type2=-999 ,  fit_subjet2_e_type2 =-999 ;
    fit_lvj_m_type2= -999  , fit_lv_m_type2 =-999 , fit_j_m_type2=-999, fit_chi2_type2=-999 , fit_lvj_pt_type2 = -999, fit_lvj_eta_type2 = -999, fit_lvj_phi_type2 = -999 ;
    fit_lvj_phi_type2=-999 ; fit_subjet1_m_type2 =-999; fit_subjet2_m_type2=-999; 
    fit_NDF_type2=-999, fit_status_type2=-999 ;

    fit_el_px_type0_met=-999 ,  fit_el_py_type0_met=-999 ,  fit_el_pz_type0_met=-999 ,  fit_el_e_type0_met=-999 ;
    fit_nv_px_type0_met=-999 ,  fit_nv_py_type0_met=-999 ,  fit_nv_pz_type0_met=-999 ,  fit_nv_e_type0_met=-999 ;
    fit_subjet1_px_type0_met=-999 ,  fit_subjet1_py_type0_met=-999 ,  fit_subjet1_pz_type0_met=-999 ,  fit_subjet1_e_type0_met =-999 ;
    fit_subjet2_px_type0_met=-999 ,  fit_subjet2_py_type0_met=-999 ,  fit_subjet2_pz_type0_met=-999 ,  fit_subjet2_e_type0_met =-999 ;
    fit_lvj_m_type0_met= -999  , fit_lv_m_type0_met =-999 , fit_j_m_type0_met=-999, fit_chi2_type0_met=-999 , fit_lvj_pt_type0_met = -999, fit_lvj_eta_type0_met = -999, fit_lvj_phi_type0_met = -999 ;
    fit_lvj_phi_type0_met=-999 ; fit_subjet1_m_type0_met =-999; fit_subjet2_m_type0_met=-999; 
    fit_NDF_type0_met=-999, fit_status_type0_met=-999 ;

    fit_el_px_type2_met=-999 ,  fit_el_py_type2_met=-999 ,  fit_el_pz_type2_met=-999 ,  fit_el_e_type2_met=-999 ;
    fit_nv_px_type2_met=-999 ,  fit_nv_py_type2_met=-999 ,  fit_nv_pz_type2_met=-999 ,  fit_nv_e_type2_met=-999 ;
    fit_subjet1_px_type2_met=-999 ,  fit_subjet1_py_type2_met=-999 ,  fit_subjet1_pz_type2_met=-999 ,  fit_subjet1_e_type2_met =-999 ;
    fit_subjet2_px_type2_met=-999 ,  fit_subjet2_py_type2_met=-999 ,  fit_subjet2_pz_type2_met=-999 ,  fit_subjet2_e_type2_met =-999 ;
    fit_lvj_m_type2_met= -999  , fit_lv_m_type2_met =-999 , fit_j_m_type2_met=-999, fit_chi2_type2_met=-999 , fit_lvj_pt_type2_met = -999, fit_lvj_eta_type2_met = -999, fit_lvj_phi_type2_met = -999 ;
    fit_lvj_phi_type2_met=-999 ; fit_subjet1_m_type2_met =-999; fit_subjet2_m_type2_met=-999; 
    fit_NDF_type2_met=-999, fit_status_type2_met=-999 ;


    /////////////////////

    /////////////////////

    boosted_lvj_m_type0 =-999 , boosted_j_m_type0=-999, boosted_subjet1_m_type0=-999, boosted_subjet2_m_type0=-999, boosted_lv_m_type0=-999 , boosted_lvj_pt_type0=-999;
    boosted_lvj_phi_type0=-999, boosted_lvj_eta_type0=-999, boosted_lvj_e_type0=-999 ;

    boostedW_lvj_m_type0 =-999 , boostedW_j_m_type0=-999, boostedW_subjet1_m_type0=-999, boostedW_subjet2_m_type0=-999, boostedW_lv_m_type0=-999 , boostedW_lvj_pt_type0=-999;
    boostedW_lvj_phi_type0=-999, boostedW_lvj_eta_type0=-999, boostedW_lvj_e_type0=-999 ;

    boosted_wjj_ang_ha_type0=-999, boosted_wjj_ang_hb_type0=-999, boosted_wjj_ang_hs_type0=-999, boosted_wjj_ang_phi_type0=-999, boosted_wjj_ang_phia_type0=-999, boosted_wjj_ang_phib_type0=-999;

    boosted_lvj_m_type0_met =-999 , boosted_j_m_type0_met=-999, boosted_subjet1_m_type0_met=-999, boosted_subjet2_m_type0_met=-999, boosted_lv_m_type0_met=-999 , boosted_lvj_pt_type0_met=-999;
    boosted_lvj_phi_type0_met=-999, boosted_lvj_eta_type0_met=-999, boosted_lvj_e_type0_met=-999 ;

    boostedW_lvj_m_type0_met =-999 , boostedW_j_m_type0_met=-999, boostedW_subjet1_m_type0_met=-999, boostedW_subjet2_m_type0_met=-999, boostedW_lv_m_type0_met=-999 , boostedW_lvj_pt_type0_met=-999;
    boostedW_lvj_phi_type0_met=-999, boostedW_lvj_eta_type0_met=-999, boostedW_lvj_e_type0_met=-999 ;

    boosted_wjj_ang_ha_type0_met=-999, boosted_wjj_ang_hb_type0_met=-999, boosted_wjj_ang_hs_type0_met=-999, boosted_wjj_ang_phi_type0_met=-999, boosted_wjj_ang_phia_type0_met=-999, boosted_wjj_ang_phib_type0_met=-999;

    boosted_lvj_m_type2 =-999 , boosted_j_m_type2=-999, boosted_subjet1_m_type2=-999, boosted_subjet2_m_type2=-999, boosted_lv_m_type2=-999 , boosted_lvj_pt_type2=-999;
    boosted_lvj_phi_type2=-999, boosted_lvj_eta_type2=-999, boosted_lvj_e_type2=-999 ;

    boostedW_lvj_m_type2 =-999 , boostedW_j_m_type2=-999, boostedW_subjet1_m_type2=-999, boostedW_subjet2_m_type2=-999, boostedW_lv_m_type2=-999 , boostedW_lvj_pt_type2=-999;
    boostedW_lvj_phi_type2=-999, boostedW_lvj_eta_type2=-999, boostedW_lvj_e_type2=-999 ;

    boosted_wjj_ang_ha_type2=-999, boosted_wjj_ang_hb_type2=-999, boosted_wjj_ang_hs_type2=-999, boosted_wjj_ang_phi_type2=-999, boosted_wjj_ang_phia_type2=-999, boosted_wjj_ang_phib_type2=-999;

    boosted_lvj_m_type2_met =-999 , boosted_j_m_type2_met=-999, boosted_subjet1_m_type2_met=-999, boosted_subjet2_m_type2_met=-999, boosted_lv_m_type2_met=-999 , boosted_lvj_pt_type2_met=-999;
    boosted_lvj_phi_type2_met=-999, boosted_lvj_eta_type2_met=-999, boosted_lvj_e_type2_met=-999 ;

    boostedW_lvj_m_type2_met =-999 , boostedW_j_m_type2_met=-999, boostedW_subjet1_m_type2_met=-999, boostedW_subjet2_m_type2_met=-999, boostedW_lv_m_type2_met=-999 , boostedW_lvj_pt_type2_met=-999;
    boostedW_lvj_phi_type2_met=-999, boostedW_lvj_eta_type2_met=-999, boostedW_lvj_e_type2_met=-999 ;

    boosted_wjj_ang_ha_type2_met=-999, boosted_wjj_ang_hb_type2_met=-999, boosted_wjj_ang_hs_type2_met=-999, boosted_wjj_ang_phi_type2_met=-999, boosted_wjj_ang_phia_type2_met=-999, boosted_wjj_ang_phib_type2_met=-999;

    //////////////

 vbf_maxpt_jj_e=-999,   vbf_maxpt_jj_pt=-999 ,   vbf_maxpt_jj_eta=-999 ,  vbf_maxpt_jj_phi=-999 , vbf_maxpt_jj_m=-999 ;
 vbf_maxpt_j1_e=-999 ,   vbf_maxpt_j1_pt=-999 ,   vbf_maxpt_j1_eta=-999 ,  vbf_maxpt_j1_phi=-999 , vbf_maxpt_j1_m=-999 ;

 vbf_maxpt_jj_e=-999 ,   vbf_maxpt_jj_pt=-999 ,   vbf_maxpt_jj_eta=-999 ,  vbf_maxpt_jj_phi=-999 , vbf_maxpt_jj_m=-999 ;
 vbf_maxpt_j1_e=-999 ,   vbf_maxpt_j1_pt=-999 ,   vbf_maxpt_j1_eta=-999 ,  vbf_maxpt_j1_phi=-999 , vbf_maxpt_j1_m=-999 ;
 vbf_maxpt_j2_e=-999 ,   vbf_maxpt_j2_pt=-999 ,   vbf_maxpt_j2_eta=-999 ,  vbf_maxpt_j2_phi=-999 , vbf_maxpt_j2_m=-999 ;
 vbf_maxpt_jj_deta=-999 ,vbf_maxpt_jj_dphi=-999;

 vbf_maxpt_j1_QGLikelihood=-999,  vbf_maxpt_j2_QGLikelihood=-999;

 vbf_maxpt_j1_isPileUpLoose =false  , vbf_maxpt_j2_isPileUpLoose =false  ;
 vbf_maxpt_j1_isPileUpMedium=false  , vbf_maxpt_j2_isPileUpMedium=false  ;
 vbf_maxpt_j1_isPileUpTight =false  , vbf_maxpt_j2_isPileUpTight =false  ;

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

 vbf_maxDeta_j1_isPileUpLoose =false  , vbf_maxDeta_j2_isPileUpLoose =false  ;
 vbf_maxDeta_j1_isPileUpMedium=false  , vbf_maxDeta_j2_isPileUpMedium=false  ;
 vbf_maxDeta_j1_isPileUpTight =false  , vbf_maxDeta_j2_isPileUpTight =false  ;

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

 vbf_maxMjj_j1_isPileUpLoose =false  , vbf_maxMjj_j2_isPileUpLoose =false  ;
 vbf_maxMjj_j1_isPileUpMedium=false  , vbf_maxMjj_j2_isPileUpMedium=false ;
 vbf_maxMjj_j1_isPileUpTight =false  , vbf_maxMjj_j2_isPileUpTight =false  ;

 vbf_maxMjj_jj_type=-999,   vbf_maxMjj_n_excj=-999,   vbf_maxMjj_n_exfj=-999 ;


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

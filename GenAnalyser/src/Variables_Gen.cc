#include "Variables_Gen.h"


void InitializeTree(Variables_Gen& vars, const std::string& outputRootFileName)
{
 
 vars.m_outputRootFile = new TFile(outputRootFileName.c_str(), "RECREATE");  
 
 ///--------------------------
 ///---- Efficiency tree ----
 ///--------------------------
 
 vars.m_efficiencyTree = new TTree("outTreeSelections", "outTreeSelections");
 vars.m_efficiencyTree -> SetDirectory(vars.m_outputRootFile);
 
 
 vars.m_efficiencyTree -> Branch("XSection", &vars.XSection, "XSection/D");
 vars.m_efficiencyTree -> Branch("XSectionErrorUp", &vars.XSectionErrorUp, "XSectionErrorUp/D");
 vars.m_efficiencyTree -> Branch("XSectionErrorDown", &vars.XSectionErrorDown, "XSectionErrorDown/D");
 vars.m_efficiencyTree -> Branch("numEntriesBefore", &vars.numEntriesBefore, "numEntriesBefore/I");
 vars.m_efficiencyTree -> Branch("preselection_efficiency", &vars.preselection_efficiency, "preselection_efficiency/D");
 
 
 
 ///-------------------------
 ///---- Reduced tree ----
 ///-------------------------
 
 vars.m_reducedTree = new TTree("outTreeJetLep", "outTreeJetLep");
 vars.m_reducedTree -> SetDirectory(vars.m_outputRootFile);
 
 vars.m_reducedTree -> Branch("runId",        &vars.runId,               "runId/I");
 vars.m_reducedTree -> Branch("lumiId",       &vars.lumiId,             "lumiId/I");
 vars.m_reducedTree -> Branch("eventId",      &vars.eventId,           "eventId/I");
 
 ///~~~~ lepton variables
 vars.m_reducedTree -> Branch("l1_pX",  &vars.l1_pX,   "l1_pX/D");
 vars.m_reducedTree -> Branch("l1_pY",  &vars.l1_pY,   "l1_pY/D");
 vars.m_reducedTree -> Branch("l1_pZ",  &vars.l1_pZ,   "l1_pZ/D");
 vars.m_reducedTree -> Branch("l1_pT",  &vars.l1_pT,   "l1_pT/D");
 vars.m_reducedTree -> Branch("l1_E",  &vars.l1_E,   "l1_E/D");
 vars.m_reducedTree -> Branch("l1_Eta",  &vars.l1_Eta,   "l1_Eta/D");
 vars.m_reducedTree -> Branch("l1_Phi",  &vars.l1_Phi,   "l1_Phi/D"); 
 
 vars.m_reducedTree -> Branch("l1_charge",  &vars.l1_charge,   "l1_charge/D");
 vars.m_reducedTree -> Branch("l1_flavour", &vars.l1_flavour, "l1_flavour/I");
 
 vars.m_reducedTree -> Branch("l2_pX",  &vars.l2_pX,   "l2_pX/D");
 vars.m_reducedTree -> Branch("l2_pY",  &vars.l2_pY,   "l2_pY/D");
 vars.m_reducedTree -> Branch("l2_pZ",  &vars.l2_pZ,   "l2_pZ/D");
 vars.m_reducedTree -> Branch("l2_pT",  &vars.l2_pT,   "l2_pT/D");
 vars.m_reducedTree -> Branch("l2_E",  &vars.l2_E,   "l2_E/D");
 vars.m_reducedTree -> Branch("l2_Eta",  &vars.l2_Eta,   "l2_Eta/D");
 vars.m_reducedTree -> Branch("l2_Phi",  &vars.l2_Phi,   "l2_Phi/D"); 
 
 vars.m_reducedTree -> Branch("l2_charge",  &vars.l2_charge,   "l2_charge/D");
 vars.m_reducedTree -> Branch("l2_flavour", &vars.l2_flavour, "l2_flavour/I");

 vars.m_reducedTree -> Branch("M_ll", &vars.M_ll, "M_ll/D");
 vars.m_reducedTree -> Branch("DEta_ll", &vars.DEta_ll, "DEta_ll/D");
 vars.m_reducedTree -> Branch("DPhi_ll", &vars.DPhi_ll, "DPhi_ll/D");
 
 vars.m_reducedTree -> Branch("l1_Z", &vars.l1_Z, "l1_Z/D");
 vars.m_reducedTree -> Branch("l2_Z", &vars.l2_Z, "l2_Z/D");
 vars.m_reducedTree -> Branch("Z_ll", &vars.Z_ll, "Z_ll/D");
 
 
 vars.m_reducedTree -> Branch("Nleptons_pT5", &vars.Nleptons_pT5, "Nleptons_pT5/I");
 vars.m_reducedTree -> Branch("Nleptons_pT10", &vars.Nleptons_pT10, "Nleptons_pT10/I");
 vars.m_reducedTree -> Branch("Nleptons_pT15", &vars.Nleptons_pT15, "Nleptons_pT15/I");
 vars.m_reducedTree -> Branch("Nleptons_pT20", &vars.Nleptons_pT20, "Nleptons_pT20/I");
 vars.m_reducedTree -> Branch("Nleptons_pT25", &vars.Nleptons_pT25, "Nleptons_pT25/I");
 vars.m_reducedTree -> Branch("Nleptons_pT30", &vars.Nleptons_pT30, "Nleptons_pT30/I");
 
 ///~~~~ met variables
 vars.m_reducedTree-> Branch("met_X",    &vars.met_X,       "met_X/D");
 vars.m_reducedTree-> Branch("met_Y",    &vars.met_Y,       "met_Y/D");
 vars.m_reducedTree-> Branch("met",        &vars.met,           "met/D");
 vars.m_reducedTree-> Branch("pmet",      &vars.pmet,         "pmet/D");
 
 ///~~ additional variables 
 vars.m_reducedTree -> Branch("mT", &vars.mT, "mT/D");
 vars.m_reducedTree -> Branch("maxDPhiJet_ll", &vars.maxDPhiJet_ll, "maxDPhiJet_ll/D");
 vars.m_reducedTree -> Branch("DPhiSingleJet_ll",    &vars.DPhiSingleJet_ll,    "DPhiSingleJet_ll/D");
 vars.m_reducedTree -> Branch("DPhiDoubleJet_ll", &vars.DPhiDoubleJet_ll, "DPhiDoubleJet_ll/D");


 ///~~~~ jet variables
 
 vars.m_reducedTree -> Branch("q1_pX", &vars.q1_pX, "q1_pX/D"); 
 vars.m_reducedTree -> Branch("q1_pY", &vars.q1_pY, "q1_pY/D"); 
 vars.m_reducedTree -> Branch("q1_pZ", &vars.q1_pZ, "q1_pZ/D"); 
 vars.m_reducedTree -> Branch("q1_pT", &vars.q1_pT, "q1_pT/D"); 
 vars.m_reducedTree -> Branch("q1_E", &vars.q1_E, "q1_E/D"); 
 vars.m_reducedTree -> Branch("q1_Eta", &vars.q1_Eta, "q1_Eta/D"); 
 vars.m_reducedTree -> Branch("q1_Phi", &vars.q1_Phi, "q1_Phi/D"); 
 vars.m_reducedTree -> Branch("q2_pX", &vars.q2_pX, "q2_pX/D"); 
 vars.m_reducedTree -> Branch("q2_pY", &vars.q2_pY, "q2_pY/D"); 
 vars.m_reducedTree -> Branch("q2_pZ", &vars.q2_pZ, "q2_pZ/D"); 
 vars.m_reducedTree -> Branch("q2_pT", &vars.q2_pT, "q2_pT/D"); 
 vars.m_reducedTree -> Branch("q2_E", &vars.q2_E, "q2_E/D"); 
 vars.m_reducedTree -> Branch("q2_Eta", &vars.q2_Eta, "q2_Eta/D"); 
 vars.m_reducedTree -> Branch("q2_Phi", &vars.q2_Phi, "q2_Phi/D"); 
 vars.m_reducedTree -> Branch("q3_pX", &vars.q3_pX, "q3_pX/D"); 
 vars.m_reducedTree -> Branch("q3_pY", &vars.q3_pY, "q3_pY/D"); 
 vars.m_reducedTree -> Branch("q3_pZ", &vars.q3_pZ, "q3_pZ/D"); 
 vars.m_reducedTree -> Branch("q3_pT", &vars.q3_pT, "q3_pT/D"); 
 vars.m_reducedTree -> Branch("q3_E", &vars.q3_E, "q3_E/D"); 
 vars.m_reducedTree -> Branch("q3_Eta", &vars.q3_Eta, "q3_Eta/D"); 
 vars.m_reducedTree -> Branch("q3_Phi", &vars.q3_Phi, "q3_Phi/D"); 
 
 
 vars.m_reducedTree -> Branch("M_qq", &vars.M_qq, "M_qq/D"); 
 vars.m_reducedTree -> Branch("DEta_qq", &vars.DEta_qq, "DEta_qq/D"); 
 vars.m_reducedTree -> Branch("DPhi_qq", &vars.DPhi_qq, "DPhi_qq/D"); 
 
 vars.m_reducedTree -> Branch("JV_20", &vars.JV_20, "JV_20/I"); 
 vars.m_reducedTree -> Branch("JV_30", &vars.JV_30, "JV_30/I"); 
 vars.m_reducedTree -> Branch("JV_40", &vars.JV_40, "JV_40/I"); 
 vars.m_reducedTree -> Branch("CJV_20", &vars.CJV_20, "CJV_20/I"); 
 vars.m_reducedTree -> Branch("CJV_30", &vars.CJV_30, "CJV_30/I"); 
 vars.m_reducedTree -> Branch("CJV_40", &vars.CJV_40, "CJV_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_01_20", &vars.q_Z_01_20, "q_Z_01_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_03_20", &vars.q_Z_03_20, "q_Z_03_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_05_20", &vars.q_Z_05_20, "q_Z_05_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_07_20", &vars.q_Z_07_20, "q_Z_07_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_09_20", &vars.q_Z_09_20, "q_Z_09_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_10_20", &vars.q_Z_10_20, "q_Z_10_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_12_20", &vars.q_Z_12_20, "q_Z_12_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_14_20", &vars.q_Z_14_20, "q_Z_14_20/I"); 
 vars.m_reducedTree -> Branch("q_Z_01_30", &vars.q_Z_01_30, "q_Z_01_30/I"); 
 vars.m_reducedTree -> Branch("q_Z_03_30", &vars.q_Z_03_30, "q_Z_03_30/I"); 
 vars.m_reducedTree -> Branch("q_Z_05_30", &vars.q_Z_05_30, "q_Z_05_30/I"); 
 vars.m_reducedTree -> Branch("q_Z_07_30", &vars.q_Z_07_30, "q_Z_07_30/I"); 
 vars.m_reducedTree -> Branch("q_Z_09_30", &vars.q_Z_09_30, "q_Z_09_30/I"); 
 vars.m_reducedTree -> Branch("q_Z_10_30", &vars.q_Z_10_30, "q_Z_10_30/I"); 
 vars.m_reducedTree -> Branch("q_Z_12_30", &vars.q_Z_12_30, "q_Z_12_30/I"); 
 vars.m_reducedTree -> Branch("q_Z_14_30", &vars.q_Z_14_30, "q_Z_14_30/I");  
 vars.m_reducedTree -> Branch("q_Z_01_40", &vars.q_Z_01_40, "q_Z_01_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_03_40", &vars.q_Z_03_40, "q_Z_03_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_05_40", &vars.q_Z_05_40, "q_Z_05_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_07_40", &vars.q_Z_07_40, "q_Z_07_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_09_40", &vars.q_Z_09_40, "q_Z_09_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_10_40", &vars.q_Z_10_40, "q_Z_10_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_12_40", &vars.q_Z_12_40, "q_Z_12_40/I"); 
 vars.m_reducedTree -> Branch("q_Z_14_40", &vars.q_Z_14_40, "q_Z_14_40/I"); 
 
 ///==== MC Higgs and WW variables ====
 
 vars.m_reducedTree -> Branch("H_pT",&vars.H_pT,"H_pT/D");
 vars.m_reducedTree -> Branch("H_pX",&vars.H_pX,"H_pX/D");
 vars.m_reducedTree -> Branch("H_pY",&vars.H_pY,"H_pY/D");
 vars.m_reducedTree -> Branch("H_pZ",&vars.H_pZ,"H_pZ/D");
 vars.m_reducedTree -> Branch("H_E",&vars.H_E,"H_E/D");
 vars.m_reducedTree -> Branch("H_eta",&vars.H_eta,"H_eta/D");
 vars.m_reducedTree -> Branch("H_phi",&vars.H_phi,"H_phi/D");
 vars.m_reducedTree -> Branch("H_mt",&vars.H_mt,"H_mt/D");
 
 vars.m_reducedTree -> Branch("mc_Q1pT",&vars.mc_Q1pT,"mc_Q1pT/D");
 vars.m_reducedTree -> Branch("mc_Q1pX",&vars.mc_Q1pX,"mc_Q1pX/D");
 vars.m_reducedTree -> Branch("mc_Q1pY",&vars.mc_Q1pY,"mc_Q1pY/D");
 vars.m_reducedTree -> Branch("mc_Q1pZ",&vars.mc_Q1pZ,"mc_Q1pZ/D");
 vars.m_reducedTree -> Branch("mc_Q1E",&vars.mc_Q1E,"mc_Q1E/D");
 vars.m_reducedTree -> Branch("mc_Q1eta",&vars.mc_Q1eta,"mc_Q1eta/D");
 vars.m_reducedTree -> Branch("mc_Q1pdgTag",&vars.mc_Q1pdgTag,"mc_Q1pdgTag/D");
 vars.m_reducedTree -> Branch("mc_Q1charge",&vars.mc_Q1charge,"mc_Q1charge/D");
 
 vars.m_reducedTree -> Branch("mc_Q2pT",&vars.mc_Q2pT,"mc_Q2pT/D");
 vars.m_reducedTree -> Branch("mc_Q2pX",&vars.mc_Q2pX,"mc_Q2pX/D");
 vars.m_reducedTree -> Branch("mc_Q2pY",&vars.mc_Q2pY,"mc_Q2pY/D");
 vars.m_reducedTree -> Branch("mc_Q2pZ",&vars.mc_Q2pZ,"mc_Q2pZ/D");
 vars.m_reducedTree -> Branch("mc_Q2E",&vars.mc_Q2E,"mc_Q2E/D");
 vars.m_reducedTree -> Branch("mc_Q2eta",&vars.mc_Q2eta,"mc_Q2eta/D");
 vars.m_reducedTree -> Branch("mc_Q2pdgTag",&vars.mc_Q2pdgTag,"mc_Q2pdgTag/D");
 vars.m_reducedTree -> Branch("mc_Q2charge",&vars.mc_Q2charge,"mc_Q2charge/D");
 
 vars.m_reducedTree -> Branch("mc_V1pT",&vars.mc_V1pT,"mc_V1pT/D");
 vars.m_reducedTree -> Branch("mc_V1pX",&vars.mc_V1pX,"mc_V1pX/D");
 vars.m_reducedTree -> Branch("mc_V1pY",&vars.mc_V1pY,"mc_V1pY/D");
 vars.m_reducedTree -> Branch("mc_V1pZ",&vars.mc_V1pZ,"mc_V1pZ/D");
 vars.m_reducedTree -> Branch("mc_V1E",&vars.mc_V1E,"mc_V1E/D");
 vars.m_reducedTree -> Branch("mc_V1eta",&vars.mc_V1eta,"mc_V1eta/D");
 vars.m_reducedTree -> Branch("mc_V1pdgTag",&vars.mc_V1pdgTag,"mc_V1pdgTag/D");
 vars.m_reducedTree -> Branch("mc_V1charge",&vars.mc_V1charge,"mc_V1charge/D");
 
 vars.m_reducedTree -> Branch("mc_V2pT",&vars.mc_V2pT,"mc_V2pT/D");
 vars.m_reducedTree -> Branch("mc_V2pX",&vars.mc_V2pX,"mc_V2pX/D");
 vars.m_reducedTree -> Branch("mc_V2pY",&vars.mc_V2pY,"mc_v2pY/D");
 vars.m_reducedTree -> Branch("mc_V2pZ",&vars.mc_V2pZ,"mc_V2pZ/D");
 vars.m_reducedTree -> Branch("mc_V2E",&vars.mc_V2E,"mc_V2E/D");
 vars.m_reducedTree -> Branch("mc_V2eta",&vars.mc_V2eta,"mc_V2eta/D");
 vars.m_reducedTree -> Branch("mc_V2pdgTag",&vars.mc_V2pdgTag,"mc_V2pdgTag/D");
 vars.m_reducedTree -> Branch("mc_V2charge",&vars.mc_V2charge,"mc_V2charge/D");
 
 vars.m_reducedTree -> Branch("mc_l1_V1pT",&vars.mc_l1_V1pT,"mc_l1_V1pT/D");
 vars.m_reducedTree -> Branch("mc_l1_V1pX",&vars.mc_l1_V1pX,"mc_l1_V1pX/D");
 vars.m_reducedTree -> Branch("mc_l1_V1pY",&vars.mc_l1_V1pY,"mc_l1_v1pY/D");
 vars.m_reducedTree -> Branch("mc_l1_V1pZ",&vars.mc_l1_V1pZ,"mc_l1_V1pZ/D");
 vars.m_reducedTree -> Branch("mc_l1_V1E",&vars.mc_l1_V1E,"mc_l1_V1E/D");
 vars.m_reducedTree -> Branch("mc_l1_V1eta",&vars.mc_l1_V1eta,"mc_l1_V1eta/D");
 vars.m_reducedTree -> Branch("mc_l1_V1pdgTag",&vars.mc_l1_V1pdgTag,"mc_l1_V1pdgTag/D");
 vars.m_reducedTree -> Branch("mc_l1_V1charge",&vars.mc_l1_V1charge,"mc_l1_V1charge/D");
 
 vars.m_reducedTree -> Branch("mc_l2_V1pT",&vars.mc_l2_V1pT,"mc_l2_V1pT/D");
 vars.m_reducedTree -> Branch("mc_l2_V1pX",&vars.mc_l2_V1pX,"mc_l2_V1pX/D");
 vars.m_reducedTree -> Branch("mc_l2_V1pY",&vars.mc_l2_V1pY,"mc_l2_v1pY/D");
 vars.m_reducedTree -> Branch("mc_l2_V1pZ",&vars.mc_l2_V1pZ,"mc_l2_V1pZ/D");
 vars.m_reducedTree -> Branch("mc_l2_V1E",&vars.mc_l2_V1E,"mc_l2_V1E/D");
 vars.m_reducedTree -> Branch("mc_l2_V1eta",&vars.mc_l2_V1eta,"mc_l2-V1eta/D");
 vars.m_reducedTree -> Branch("mc_l2_V1pdgTag",&vars.mc_l2_V1pdgTag,"mc_l2_V1pdgTag/D");
 vars.m_reducedTree -> Branch("mc_l2_V1charge",&vars.mc_l2_V1charge,"mc_l2_V1charge/D");
 
 vars.m_reducedTree -> Branch("mc_l2_V2pT",&vars.mc_l2_V2pT,"mc_l2_V2pT/D");
 vars.m_reducedTree -> Branch("mc_l2_V2pX",&vars.mc_l2_V2pX,"mc_l2_V2pX/D");
 vars.m_reducedTree -> Branch("mc_l2_V2pY",&vars.mc_l2_V2pY,"mc_l2_v2pY/D");
 vars.m_reducedTree -> Branch("mc_l2_V2pZ",&vars.mc_l2_V2pZ,"mc_l2_V2pZ/D");
 vars.m_reducedTree -> Branch("mc_l2_V2E",&vars.mc_l2_V2E,"mc_l2_V2E/D");
 vars.m_reducedTree -> Branch("mc_l2_V2eta",&vars.mc_l2_V2eta,"mc_l2_V2eta/D");
 vars.m_reducedTree -> Branch("mc_l2_V2pdgTag",&vars.mc_l2_V2pdgTag,"mc_l2_V2pdgTag/D");
 vars.m_reducedTree -> Branch("mc_l2_V2charge",&vars.mc_l2_V2charge,"mc_l2_V2charge/D");
 
 vars.m_reducedTree -> Branch("mc_l1_V2pT",&vars.mc_l1_V2pT,"mc_l1_V2pT/D");
 vars.m_reducedTree -> Branch("mc_l1_V2pX",&vars.mc_l1_V2pX,"mc_l1_V2pX/D");
 vars.m_reducedTree -> Branch("mc_l1_V2pY",&vars.mc_l1_V2pY,"mc_l1_v2pY/D");
 vars.m_reducedTree -> Branch("mc_l1_V2pZ",&vars.mc_l1_V2pZ,"mc_l1_V2pZ/D");
 vars.m_reducedTree -> Branch("mc_l1_V2E",&vars.mc_l1_V2E,"mc_l1_V2E/D");
 vars.m_reducedTree -> Branch("mc_l1_V2eta",&vars.mc_l1_V2eta,"mc_l1_V2eta/D");
 vars.m_reducedTree -> Branch("mc_l1_V2pdgTag",&vars.mc_l1_V2pdgTag,"mc_l1_V2pdgTag/D");
 vars.m_reducedTree -> Branch("mc_l1_V2charge",&vars.mc_l1_V2charge,"mc_l1_V2charge/D");
 
}


void FillTree(Variables_Gen& vars){
 vars.m_reducedTree -> Fill();
}

void FillEfficiencyTree(Variables_Gen& vars){
 vars.m_efficiencyTree -> Fill();
}

void SaveTree(Variables_Gen& vars)
{
 // save tree
 vars.m_outputRootFile -> cd();
 vars.m_reducedTree -> Write();
 vars.m_efficiencyTree -> Write();
 vars.m_outputRootFile -> Close();
}


void SetLeptonsVariables(Variables_Gen& vars, treeReader& reader,const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2)
{

if (FlavourLep1 == 13) {
  vars.l1_pX = reader.Get4V("mcMu")->at(iLep1).X();
  vars.l1_pY = reader.Get4V("mcMu")->at(iLep1).Y();
  vars.l1_pZ = reader.Get4V("mcMu")->at(iLep1).Z();
  vars.l1_pT = reader.Get4V("mcMu")->at(iLep1).pt(); 
  vars.l1_E  = reader.Get4V("mcMu")->at(iLep1).E();
  vars.l1_Eta = reader.Get4V("mcMu")->at(iLep1).Eta();
  vars.l1_Phi = reader.Get4V("mcMu")->at(iLep1).Phi(); 
  
  vars.l1_charge = reader.GetFloat("mcMu_charge")->at(iLep1);
  vars.l1_flavour = 13;
 }
 
 if (FlavourLep2 == 13) {
  vars.l2_pX = reader.Get4V("mcMu")->at(iLep2).X();
  vars.l2_pY = reader.Get4V("mcMu")->at(iLep2).Y();
  vars.l2_pZ = reader.Get4V("mcMu")->at(iLep2).Z();
  vars.l2_pT = reader.Get4V("mcMu")->at(iLep2).pt(); 
  vars.l2_E = reader.Get4V("mcMu")->at(iLep2).E();
  vars.l2_Eta = reader.Get4V("mcMu")->at(iLep2).Eta();
  vars.l2_Phi = reader.Get4V("mcMu")->at(iLep2).Phi();
  
  vars.l2_charge = reader.GetFloat("mcMu_charge")->at(iLep2);
  vars.l2_flavour = 13;
  
 }
 
 if (FlavourLep1 == 11) {
  vars.l1_pX = reader.Get4V("mcEle")->at(iLep1).X();
  vars.l1_pY = reader.Get4V("mcEle")->at(iLep1).Y();
  vars.l1_pZ = reader.Get4V("mcEle")->at(iLep1).Z();
  vars.l1_pT = reader.Get4V("mcEle")->at(iLep1).pt(); 
  vars.l1_E = reader.Get4V("mcEle")->at(iLep1).E();
  vars.l1_Eta = reader.Get4V("mcEle")->at(iLep1).Eta();
  vars.l1_Phi = reader.Get4V("mcEle")->at(iLep1).Phi();
  
  vars.l1_charge = reader.GetFloat("mcEle_charge")->at(iLep1);
  vars.l1_flavour = 11;
 }
 if (FlavourLep2 == 11) {
  
  vars.l2_pX = reader.Get4V("mcEle")->at(iLep2).X();
  vars.l2_pY = reader.Get4V("mcEle")->at(iLep2).Y();
  vars.l2_pZ = reader.Get4V("mcEle")->at(iLep2).Z();
  vars.l2_pT = reader.Get4V("mcEle")->at(iLep2).pt(); 
  vars.l2_E = reader.Get4V("mcEle")->at(iLep2).E();
  vars.l2_Eta = reader.Get4V("mcEle")->at(iLep2).Eta();
  vars.l2_Phi = reader.Get4V("mcEle")->at(iLep2).Phi();
  vars.l2_charge = reader.GetFloat("mcEle_charge")->at(iLep2);
  vars.l2_flavour = 11;
 
 }
 
 double etaMean = (vars.q1_Eta + vars.q2_Eta) / 2.;
 double dEta = fabs(vars.q1_Eta - vars.q2_Eta);

 ///==== Zepp for lepton ====
 ///==== and ====
 ///==== Zepp for lepton system ====
 
 if (FlavourLep1 == 11 && FlavourLep2 == 11) {
  vars.M_ll = (reader.Get4V("mcEle")->at(iLep1) + reader.Get4V("mcEle")->at(iLep2)).mass();
  vars.DEta_ll = deltaEta(reader.Get4V("mcEle")->at(iLep1).Eta() , reader.Get4V("mcEle")->at(iLep2).Eta());
  vars.DPhi_ll = deltaPhi(reader.Get4V("mcEle")->at(iLep1).Phi() , reader.Get4V("mcEle")->at(iLep2).Phi());  
  vars.l1_Z = (reader.Get4V("mcEle")->at(iLep1).Eta() - etaMean)/dEta;
  vars.l2_Z = (reader.Get4V("mcEle")->at(iLep2).Eta() - etaMean)/dEta;
  vars.Z_ll = ((reader.Get4V("mcEle")->at(iLep1) + reader.Get4V("mcEle")->at(iLep2)).Eta() - etaMean)/dEta;
 }
 if (FlavourLep1 == 11 && FlavourLep2 == 13) {
  vars.M_ll = (reader.Get4V("mcEle")->at(iLep1) + reader.Get4V("mcMu")->at(iLep2)).mass();
  vars.DEta_ll = deltaEta(reader.Get4V("mcEle")->at(iLep1).Eta() , reader.Get4V("mcMu")->at(iLep2).Eta());
  vars.DPhi_ll = deltaPhi(reader.Get4V("mcEle")->at(iLep1).Phi() , reader.Get4V("mcMu")->at(iLep2).Phi());  
  vars.l1_Z = (reader.Get4V("mcEle")->at(iLep1).Eta() - etaMean)/dEta;
  vars.l2_Z = (reader.Get4V("mcMu")->at(iLep2).Eta() - etaMean)/dEta;
  vars.Z_ll = ((reader.Get4V("mcEle")->at(iLep1) + reader.Get4V("mcMu")->at(iLep2)).Eta() - etaMean)/dEta;
 }
 if (FlavourLep1 == 13 && FlavourLep2 == 11) {
  vars.M_ll = (reader.Get4V("mcMu")->at(iLep1) + reader.Get4V("mcEle")->at(iLep2)).mass();
  vars.DEta_ll = deltaEta(reader.Get4V("mcMu")->at(iLep1).Eta() , reader.Get4V("mcEle")->at(iLep2).Eta());
  vars.DPhi_ll = deltaPhi(reader.Get4V("mcMu")->at(iLep1).Phi() , reader.Get4V("mcEle")->at(iLep2).Phi());  
  vars.l1_Z = (reader.Get4V("mcMu")->at(iLep1).Eta() - etaMean)/dEta;
  vars.l2_Z = (reader.Get4V("mcEle")->at(iLep2).Eta() - etaMean)/dEta;
  vars.Z_ll = ((reader.Get4V("mcMu")->at(iLep1) + reader.Get4V("mcEle")->at(iLep2)).Eta() - etaMean)/dEta;
 }
 if (FlavourLep1 == 13 && FlavourLep2 == 13) {
  vars.M_ll = (reader.Get4V("mcMu")->at(iLep1) + reader.Get4V("mcMu")->at(iLep2)).mass();
  vars.DEta_ll = deltaEta(reader.Get4V("mcMu")->at(iLep1).Eta() , reader.Get4V("mcMu")->at(iLep2).Eta());
  vars.DPhi_ll = deltaPhi(reader.Get4V("mcMu")->at(iLep1).Phi() , reader.Get4V("mcMu")->at(iLep2).Phi());  
  vars.l1_Z = (reader.Get4V("mcMu")->at(iLep1).Eta() - etaMean)/dEta;
  vars.l2_Z = (reader.Get4V("mcMu")->at(iLep2).Eta() - etaMean)/dEta;
  vars.Z_ll = ((reader.Get4V("mcMu")->at(iLep1) + reader.Get4V("mcMu")->at(iLep2)).Eta() - etaMean)/dEta;
 }
}



void SetMetVariables(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2)
{
 double Lep1_phi = 100. ;
 if (FlavourLep1 == 11 ) Lep1_phi = reader.Get4V ("mcEle")->at (iLep1).Phi () ;
 else Lep1_phi = reader.Get4V ("mcMu")->at (iLep1).Phi () ;
 double Lep2_phi = 100. ;
 if (FlavourLep2 == 11 ) Lep2_phi = reader.Get4V ("mcEle")->at (iLep2).Phi () ;
 else Lep2_phi = reader.Get4V ("mcMu")->at (iLep2).Phi () ;

 double MET_phi = reader.Get4V ("mcMet")->at (0).Phi () ;
 double deltaPhiMin = std::min (deltaPhi (MET_phi, Lep1_phi), deltaPhi (MET_phi, Lep2_phi)) ;

 if (deltaPhiMin < 1.57079632679) vars.pmet = reader.Get4V("mcMet")->at(0).Et() * sin (deltaPhiMin) ;
 else vars.pmet = reader.Get4V("mcMet")->at(0).Et() ;

 vars.met_X = reader.Get4V("mcMet")->at(0).X();
 vars.met_Y = reader.Get4V("mcMet")->at(0).Y();
 vars.met = reader.Get4V("mcMet")->at(0).Et();
 
 
}

void SetEventVariables(Variables_Gen& vars, treeReader& reader)
{
 if (reader.GetInt("lumiId")->size() != 0) vars.lumiId = reader.GetInt("lumiId")->at(0);
 else vars.lumiId = -1;
 if (reader.GetInt("eventId")->size() != 0) vars.eventId = reader.GetInt("eventId")->at(0);
 else vars.eventId = -1;
 if (reader.GetInt("runId")->size() != 0) vars.runId = reader.GetInt("runId")->at(0);
 else vars.runId = -1;
}


void SetMCVariables(Variables_Gen& vars, treeReader& reader)
{
    
 if (reader.GetInt("mc_PUit_NumInteractions")->size() != 0) {
  vars.numPUMCit = reader.GetInt("mc_PUit_NumInteractions")->at(0);
 }
 else vars.numPUMCit = -1;
 
 if (reader.GetInt("mc_PUoot_NumInteractions")->size() != 0) {
 
  vars.numPUMCoot = 0;
  
  for (int iter = 0; iter<reader.GetInt("mc_PUoot_NumInteractions")->size(); iter++) {
     vars.numPUMCoot += reader.GetInt("mc_PUoot_NumInteractions")->at(iter);
  } 
 }
 else vars.numPUMCoot = -1;
 
 if ( vars.numPUMCit != -1 && vars.numPUMCoot != -1) {
  vars.numPUMC = (vars.numPUMCit + vars.numPUMCoot + 0.5) / 3;
 }
 else  if ( vars.numPUMCit != -1) {
   vars.numPUMC = vars.numPUMCit;
 } else   if ( vars.numPUMCoot != -1) {
   vars.numPUMC = (vars.numPUMCoot + 0.5) / 2;
 } else {
    vars.numPUMC = -1;
 }

 
 if (reader.Get4V("mcH")->size() != 0) {
    
  vars.H_pT = reader.Get4V("mcH")->at(0).pt();
  vars.H_E =  reader.Get4V("mcH")->at(0).E();
  vars.H_pX = reader.Get4V("mcH")->at(0).X();
  vars.H_pY = reader.Get4V("mcH")->at(0).Y();
  vars.H_pZ = reader.Get4V("mcH")->at(0).Z();
  vars.H_eta = reader.Get4V("mcH")->at(0).Eta();
  vars.H_phi = reader.Get4V("mcH")->at(0).Phi();
  vars.H_mt = reader.Get4V("mcH")->at(0).mt();
 }
 else {
  
  vars.H_pT = -1;
  vars.H_E = -1;
  vars.H_pX = -1;
  vars.H_pY = -1;
  vars.H_pZ = -1;
  vars.H_eta = -1;
  vars.H_phi = -1;
  vars.H_mt = -1;
  
 }
 
 if (reader.Get4V("mcQ1_tag")->size() != 0) {
  
  vars.mc_Q1pT = reader.Get4V("mcQ1_tag")->at(0).pt();
  vars.mc_Q1E = reader.Get4V("mcQ1_tag")->at(0).E();
  vars.mc_Q1pX = reader.Get4V("mcQ1_tag")->at(0).X();
  vars.mc_Q1pY = reader.Get4V("mcQ1_tag")->at(0).Y();
  vars.mc_Q1pZ = reader.Get4V("mcQ1_tag")->at(0).Z();
  vars.mc_Q1eta = reader.Get4V("mcQ1_tag")->at(0).Eta();
  vars.mc_Q1phi = reader.Get4V("mcQ1_tag")->at(0).Phi();
  vars.mc_Q1charge = reader.GetFloat("mcQ1_tag_charge")->at(0);
  vars.mc_Q1pdgTag = reader.GetFloat("mcQ1_tag_pdgId")->at(0);
  
  
 }
 else {
  
  vars.mc_Q1pT = -1;
  vars.mc_Q1E = -1;
  vars.mc_Q1pX = -1;
  vars.mc_Q1pY = -1;
  vars.mc_Q1pZ = -1;
  vars.mc_Q1eta = -1;
  vars.mc_Q1phi = -1;
  vars.mc_Q1pdgTag = -1;
  vars.mc_Q1charge =-1;
  
 }
 
 if (reader.Get4V("mcQ2_tag")->size() != 0) {
    
  vars.mc_Q2pT = reader.Get4V("mcQ2_tag")->at(0).Pt();
  vars.mc_Q2E = reader.Get4V("mcQ2_tag")->at(0).E();
  vars.mc_Q2pX = reader.Get4V("mcQ2_tag")->at(0).X();
  vars.mc_Q2pY = reader.Get4V("mcQ2_tag")->at(0).Y();
  vars.mc_Q2pZ = reader.Get4V("mcQ2_tag")->at(0).Z();
  vars.mc_Q2eta = reader.Get4V("mcQ2_tag")->at(0).Eta();
  vars.mc_Q2phi = reader.Get4V("mcQ2_tag")->at(0).Phi();
  vars.mc_Q2pdgTag = reader.GetFloat("mcQ2_tag_pdgId")->at(0);
  vars.mc_Q2charge = reader.GetFloat("mcQ2_tag_charge")->at(0);
 }
 else {
  
  vars.mc_Q2pT = -1;
  vars.mc_Q2E = -1;
  vars.mc_Q2pX = -1;
  vars.mc_Q2pY = -1;
  vars.mc_Q2pZ = -1;
  vars.mc_Q2eta = -1;
  vars.mc_Q2phi = -1;
  vars.mc_Q2pdgTag = -1;
  vars.mc_Q2charge =-1;
  
 }

if (reader.Get4V("mcV1")->size() != 0) {
    
  vars.mc_V1pT = reader.Get4V("mcV1")->at(0).Pt();
  vars.mc_V1E = reader.Get4V("mcV1")->at(0).E();
  vars.mc_V1pX = reader.Get4V("mcV1")->at(0).X();
  vars.mc_V1pY = reader.Get4V("mcV1")->at(0).Y();
  vars.mc_V1pZ = reader.Get4V("mcV1")->at(0).Z();
  vars.mc_V1eta = reader.Get4V("mcV1")->at(0).Eta();
  vars.mc_V1phi = reader.Get4V("mcV1")->at(0).Phi();
  vars.mc_V1pdgTag = reader.GetFloat("mcV1_pdgId")->at(0);
  vars.mc_V1charge = reader.GetFloat("mcV1_charge")->at(0);
  
 }
 else {
  
  vars.mc_V1pT = -1;
  vars.mc_V1E = -1;
  vars.mc_V1pX = -1;
  vars.mc_V1pY = -1;
  vars.mc_V1pZ = -1;
  vars.mc_V1eta = -1;
  vars.mc_V1phi = -1;
  vars.mc_V1pdgTag = -1;
  vars.mc_V1charge =-1;
  
 }

if (reader.Get4V("mcV2")->size() != 0) {
    
  vars.mc_V2pT = reader.Get4V("mcV2")->at(0).Pt();
  vars.mc_V2E = reader.Get4V("mcV2")->at(0).E();
  vars.mc_V2pX = reader.Get4V("mcV2")->at(0).X();
  vars.mc_V2pY = reader.Get4V("mcV2")->at(0).Y();
  vars.mc_V2pZ = reader.Get4V("mcV2")->at(0).Z();
  vars.mc_V2eta = reader.Get4V("mcV2")->at(0).Eta();
  vars.mc_V2phi = reader.Get4V("mcV2")->at(0).Phi();
  vars.mc_V2pdgTag = reader.GetFloat("mcV2_pdgId")->at(0);
  vars.mc_V2charge = reader.GetFloat("mcV2_charge")->at(0);
  
 }
 else {
  
  vars.mc_V2pT = -1;
  vars.mc_V2E = -1;
  vars.mc_V2pX = -1;
  vars.mc_V2pY = -1;
  vars.mc_V2pZ = -1;
  vars.mc_V2eta = -1;
  vars.mc_V2phi = -1;
  vars.mc_V2pdgTag = -1;
  vars.mc_V2charge =-1;
  
 }

if (reader.Get4V("mcF1_fromV1")->size() != 0) {
    
  vars.mc_l1_V1pT = reader.Get4V("mcF1_fromV1")->at(0).Pt();
  vars.mc_l1_V1E = reader.Get4V("mcF1_fromV1")->at(0).E();
  vars.mc_l1_V1pX = reader.Get4V("mcF1_fromV1")->at(0).X();
  vars.mc_l1_V1pY = reader.Get4V("mcF1_fromV1")->at(0).Y();
  vars.mc_l1_V1pZ = reader.Get4V("mcF1_fromV1")->at(0).Z();
  vars.mc_l1_V1eta = reader.Get4V("mcF1_fromV1")->at(0).Eta();
  vars.mc_l1_V1phi = reader.Get4V("mcF1_fromV1")->at(0).Phi();
  vars.mc_l1_V1pdgTag = reader.GetFloat("mcF1_fromV1_charge")->at(0);
  vars.mc_l1_V1charge = reader.GetFloat("mcF1_fromV1_pdgId")->at(0);
  
 }
 else {
  
  vars.mc_l1_V1pT = -1;
  vars.mc_l1_V1E = -1;
  vars.mc_l1_V1pX = -1;
  vars.mc_l1_V1pY = -1;
  vars.mc_l1_V1pZ = -1;
  vars.mc_l1_V1eta = -1;
  vars.mc_l1_V1phi = -1;
  vars.mc_l1_V1pdgTag = -1;
  vars.mc_l1_V1charge =-1;
  
 }

if (reader.Get4V("mcF2_fromV1")->size() != 0) {
    
  vars.mc_l2_V1pT = reader.Get4V("mcF2_fromV1")->at(0).Pt();
  vars.mc_l2_V1E = reader.Get4V("mcF2_fromV1")->at(0).E();
  vars.mc_l2_V1pX = reader.Get4V("mcF2_fromV1")->at(0).X();
  vars.mc_l2_V1pY = reader.Get4V("mcF2_fromV1")->at(0).Y();
  vars.mc_l2_V1pZ = reader.Get4V("mcF2_fromV1")->at(0).Z();
  vars.mc_l2_V1eta = reader.Get4V("mcF2_fromV1")->at(0).Eta();
  vars.mc_l2_V1phi = reader.Get4V("mcF2_fromV1")->at(0).Phi();
  vars.mc_l2_V1pdgTag = reader.GetFloat("mcF2_fromV1_charge")->at(0);
  vars.mc_l2_V1charge = reader.GetFloat("mcF2_fromV1_pdgId")->at(0);
  
 }
 else {
  
  vars.mc_l2_V1pT = -1;
  vars.mc_l2_V1E = -1;
  vars.mc_l2_V1pX = -1;
  vars.mc_l2_V1pY = -1;
  vars.mc_l2_V1pZ = -1;
  vars.mc_l2_V1eta = -1;
  vars.mc_l2_V1phi = -1;
  vars.mc_l2_V1pdgTag = -1;
  vars.mc_l2_V1charge =-1;
  
 }

if (reader.Get4V("mcF1_fromV2")->size() != 0) {
    
  vars.mc_l1_V2pT = reader.Get4V("mcF1_fromV2")->at(0).Pt();
  vars.mc_l1_V2E = reader.Get4V("mcF1_fromV2")->at(0).E();
  vars.mc_l1_V2pX = reader.Get4V("mcF1_fromV2")->at(0).X();
  vars.mc_l1_V2pY = reader.Get4V("mcF1_fromV2")->at(0).Y();
  vars.mc_l1_V2pZ = reader.Get4V("mcF1_fromV2")->at(0).Z();
  vars.mc_l1_V2eta = reader.Get4V("mcF1_fromV2")->at(0).Eta();
  vars.mc_l1_V2phi = reader.Get4V("mcF1_fromV2")->at(0).Phi();
  vars.mc_l1_V2pdgTag = reader.GetFloat("mcF1_fromV2_charge")->at(0);
  vars.mc_l1_V2charge = reader.GetFloat("mcF1_fromV2_pdgId")->at(0);
  
 }
 else {
  
  vars.mc_l1_V2pT = -1;
  vars.mc_l1_V2E = -1;
  vars.mc_l1_V2pX = -1;
  vars.mc_l1_V2pY = -1;
  vars.mc_l1_V2pZ = -1;
  vars.mc_l1_V2eta = -1;
  vars.mc_l1_V2phi = -1;
  vars.mc_l1_V2pdgTag = -1;
  vars.mc_l1_V2charge =-1;
  
 }
 
 
if (reader.Get4V("mcF2_fromV2")->size() != 0) {
    
  vars.mc_l2_V2pT = reader.Get4V("mcF2_fromV2")->at(0).Pt();
  vars.mc_l2_V2E = reader.Get4V("mcF2_fromV2")->at(0).E();
  vars.mc_l2_V2pX = reader.Get4V("mcF2_fromV2")->at(0).X();
  vars.mc_l2_V2pY = reader.Get4V("mcF2_fromV2")->at(0).Y();
  vars.mc_l2_V2pZ = reader.Get4V("mcF2_fromV2")->at(0).Z();
  vars.mc_l2_V2eta = reader.Get4V("mcF2_fromV2")->at(0).Eta();
  vars.mc_l2_V2phi = reader.Get4V("mcF2_fromV2")->at(0).Phi();
  vars.mc_l2_V2pdgTag = reader.GetFloat("mcF2_fromV2_charge")->at(0);
  vars.mc_l2_V2charge = reader.GetFloat("mcF2_fromV2_pdgId")->at(0);
  
 }
 else {
  
  vars.mc_l2_V2pT = -1;
  vars.mc_l2_V2E = -1;
  vars.mc_l2_V2pX = -1;
  vars.mc_l2_V2pY = -1;
  vars.mc_l2_V2pZ = -1;
  vars.mc_l2_V2eta = -1;
  vars.mc_l2_V2phi = -1;
  vars.mc_l2_V2pdgTag = -1;
  vars.mc_l2_V2charge =-1;
  
 }

 
}




void SetQJetVariables(Variables_Gen& vars, treeReader& reader, const int& q1, const int& q2, const std::vector<int>& blacklistJet_forCJV)
{

if(q1!=-1 && q2!=-1)
{  
 vars.q1_pX = reader.Get4V("mcJet")->at(q1).X();
 vars.q1_pY = reader.Get4V("mcJet")->at(q1).Y();
 vars.q1_pZ = reader.Get4V("mcJet")->at(q1).Z();
 vars.q1_pT = reader.Get4V("mcJet")->at(q1).pt();
 vars.q1_E = reader.Get4V("mcJet")->at(q1).E();
 vars.q1_Eta = reader.Get4V("mcJet")->at(q1).Eta();
 vars.q1_Phi = reader.Get4V("mcJet")->at(q1).Phi();
 
 vars.q2_pX = reader.Get4V("mcJet")->at(q2).X();
 vars.q2_pY = reader.Get4V("mcJet")->at(q2).Y();
 vars.q2_pZ = reader.Get4V("mcJet")->at(q2).Z();
 vars.q2_pT = reader.Get4V("mcJet")->at(q2).pt();
 vars.q2_E = reader.Get4V("mcJet")->at(q2).E();
 vars.q2_Eta = reader.Get4V("mcJet")->at(q2).Eta();
 vars.q2_Phi = reader.Get4V("mcJet")->at(q2).Phi();

 vars.M_qq = (reader.Get4V("mcJet")->at(q1) + reader.Get4V("mcJet")->at(q2)).mass();
 vars.DEta_qq = deltaEta(reader.Get4V("mcJet")->at(q1).Eta() , reader.Get4V("mcJet")->at(q2).Eta());
 vars.DPhi_qq = deltaPhi(reader.Get4V("mcJet")->at(q1).Phi() , reader.Get4V("mcJet")->at(q2).Phi());
 
 vars.CJV_20 = getCJV(*(reader.Get4V("mcJet")),q1,q2,20.,&blacklistJet_forCJV);
 vars.CJV_30 = getCJV(*(reader.Get4V("mcJet")),q1,q2,30.,&blacklistJet_forCJV);
 vars.CJV_40 = getCJV(*(reader.Get4V("mcJet")),q1,q2,40.,&blacklistJet_forCJV);
 
 vars.JV_20 = getJV(*(reader.Get4V("mcJet")),20.,&blacklistJet_forCJV);
 vars.JV_30 = getJV(*(reader.Get4V("mcJet")),30.,&blacklistJet_forCJV);
 vars.JV_40 = getJV(*(reader.Get4V("mcJet")),40.,&blacklistJet_forCJV);
 
 
 vars.q_Z_01_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,0.1,&blacklistJet_forCJV);
 vars.q_Z_03_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,0.3,&blacklistJet_forCJV);
 vars.q_Z_05_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,0.5,&blacklistJet_forCJV);
 vars.q_Z_07_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,0.7,&blacklistJet_forCJV);
 vars.q_Z_09_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,0.9,&blacklistJet_forCJV);
 vars.q_Z_10_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,1.0,&blacklistJet_forCJV);
 vars.q_Z_12_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,1.2,&blacklistJet_forCJV);
 vars.q_Z_14_20 = getZepp(*(reader.Get4V("mcJet")),q1,q2,20.,1.4,&blacklistJet_forCJV);
 
 vars.q_Z_01_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,0.1,&blacklistJet_forCJV);
 vars.q_Z_03_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,0.3,&blacklistJet_forCJV);
 vars.q_Z_05_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,0.5,&blacklistJet_forCJV);
 vars.q_Z_07_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,0.7,&blacklistJet_forCJV);
 vars.q_Z_09_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,0.9,&blacklistJet_forCJV);
 vars.q_Z_10_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,1.0,&blacklistJet_forCJV);
 vars.q_Z_12_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,1.2,&blacklistJet_forCJV);
 vars.q_Z_14_30 = getZepp(*(reader.Get4V("mcJet")),q1,q2,30.,1.4,&blacklistJet_forCJV);
 
 vars.q_Z_01_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,0.1,&blacklistJet_forCJV);
 vars.q_Z_03_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,0.3,&blacklistJet_forCJV);
 vars.q_Z_05_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,0.5,&blacklistJet_forCJV);
 vars.q_Z_07_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,0.7,&blacklistJet_forCJV);
 vars.q_Z_09_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,0.9,&blacklistJet_forCJV);
 vars.q_Z_10_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,1.0,&blacklistJet_forCJV);
 vars.q_Z_12_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,1.2,&blacklistJet_forCJV);
 vars.q_Z_14_40 = getZepp(*(reader.Get4V("mcJet")),q1,q2,40.,1.4,&blacklistJet_forCJV);
 
}

if(q1==-1 && q2==-1)
{
  vars.JV_20 = getJV(*(reader.Get4V("jets")),20.,&blacklistJet_forCJV);
  vars.JV_30 = getJV(*(reader.Get4V("jets")),30.,&blacklistJet_forCJV);
  vars.JV_40 = getJV(*(reader.Get4V("jets")),40.,&blacklistJet_forCJV);
}
  
}




///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SetMTVariable(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2){

 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  totalP4; 
 
  if(FlavourLep1==11 && FlavourLep2==11)
     totalP4=reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

  if(FlavourLep1==11 && FlavourLep2==13)
    totalP4=reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

  if(FlavourLep1==13 && FlavourLep2==11)
    totalP4=reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

   if(FlavourLep1==13 && FlavourLep2==13)
     totalP4=reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

    double DeltaPhipfMet=fabs(deltaPhi(totalP4.Phi(),reader.Get4V("mcMet")->at(0).Phi()));

    vars.mT=sqrt(2*totalP4.pt()*reader.Get4V("mcMet")->at(0).pt()*(1-cos(DeltaPhipfMet)));
    
}

///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void SetDPhiJetll(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2, const int& q1, const int& q2){

 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  totalP4; 

  if(FlavourLep1==11 && FlavourLep2==11)
     totalP4 = reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

  if(FlavourLep1==11 && FlavourLep2==13)
    totalP4 = reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

  if(FlavourLep1==13 && FlavourLep2==11)
    totalP4 = reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

   if(FlavourLep1==13 && FlavourLep2==13)
     totalP4 = reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

   
    vars.DPhiSingleJet_ll = fabs(deltaPhi(reader.Get4V("mcJet")->at(q1).Phi(),totalP4.Phi()));
    vars.DPhiDoubleJet_ll = fabs(deltaPhi((reader.Get4V("mcJet")->at(q1) + reader.Get4V("mcJet")->at(q2)).Phi(),totalP4.Phi()));


    vars.DPhiJet_ll = std::min (vars.DPhiSingleJet_ll, vars.DPhiDoubleJet_ll);
    vars.maxDPhiJet_ll = std::max (vars.DPhiSingleJet_ll, vars.DPhiDoubleJet_ll);
   
  
}

void SetDPhiJetll(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2, const std::vector<int>& blacklistJet ){

 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  totalP4; 

  if(FlavourLep1==11 && FlavourLep2==11)
     totalP4 = reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

  if(FlavourLep1==11 && FlavourLep2==13)
    totalP4 = reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

  if(FlavourLep1==13 && FlavourLep2==11)
    totalP4 = reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

   if(FlavourLep1==13 && FlavourLep2==13)
     totalP4 = reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

  
     double Jet_max_pt=-1;
     int Jet_max_pt_pos=0;

   for(int iJet=0; iJet<reader.Get4V("mcJet")->size(); iJet++)
    {
      bool skipJet = false;
      for(unsigned int kk = 0; kk < blacklistJet.size(); ++kk) {
   
	if(blacklistJet.at(kk) == static_cast<int>(iJet)) skipJet = true;}
	
      if(skipJet) continue;
	
      if(reader.Get4V("mcJet")->at(iJet).pt()<=15.0)continue;
   
      if(reader.Get4V("mcJet")->at(iJet).pt() > Jet_max_pt)
      {Jet_max_pt=reader.Get4V("mcJet")->at(iJet).pt();
       Jet_max_pt_pos=iJet;
        }
     }
     
  
  vars.DPhiJet_ll=fabs(deltaPhi(reader.Get4V("mcJet")->at(Jet_max_pt_pos).Phi(),totalP4.Phi()));
  
}



///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void SetJDVariables(Variables_Gen& vars, treeReader& reader, const int& iLep1, const int& iLep2, const int& FlavourLep1, const int& FlavourLep2, const int& q1, const int& q2){

 //-- Z direction
 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  totalZ; 
 
  if(FlavourLep1==11 && FlavourLep2==11)
     totalZ = reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

  if(FlavourLep1==11 && FlavourLep2==13)
    totalZ = reader.Get4V ("mcEle")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

  if(FlavourLep1==13 && FlavourLep2==11)
    totalZ = reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcEle")->at (iLep2);

   if(FlavourLep1==13 && FlavourLep2==13)
     totalZ = reader.Get4V ("mcMu")->at (iLep1)+reader.Get4V ("mcMu")->at (iLep2);

  //-- inversion of Z direction -> jet pair direction
   totalZ = - totalZ;

 //-- check that -Z is between jets -> angle -Z, jet1/2 < 90°
   vars.jd_dphi1 = fabs(deltaPhi(reader.Get4V("mcJet")->at(q1).Phi(),totalZ.Phi()));
   vars.jd_dphi2 = fabs(deltaPhi(reader.Get4V("mcJet")->at(q2).Phi(),totalZ.Phi()));
   
 //-- correction transverse direction of jets wrt -Z direction
  
   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  temp_q1 = reader.Get4V("mcJet")->at(q1); 
   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  temp_q2 = reader.Get4V("mcJet")->at(q2); 

   ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag>    V_q1 = temp_q1.Vect();
   ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag>    V_q2 = temp_q2.Vect();   
   ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag>    V_Z  = totalZ.Vect();
  
   V_q1.SetZ (0);
   V_q2.SetZ (0);
   V_Z.SetZ (0);
      
   ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> V_q1q2 = V_q1.Cross(V_q2);
   ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> V_q1Z = V_q1.Cross(V_Z);
   ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> V_Zq2 = V_Z.Cross(V_q2);

  if  ((V_q1q2.Z() * V_q1Z.Z()) > 0 && (V_q1Z.Z() * V_Zq2.Z()) > 0) {
   vars.jd_between = 1;
  }
  else {
   vars.jd_between = 0;
  }

//  if  (0) {
//  vars.jd_Delta_q1_P  = 1;
//  vars.jd_Delta_q2_P  = 1;
//
//   vars.jd_met = vars.met;
//   vars.jd_met_X = vars.met_X;
//   vars.jd_met_Y = vars.met_Y;
//
//   vars.jd_q1_pX = vars.q1_pX;
//   vars.jd_q1_pY = vars.q1_pY;
//   vars.jd_q1_pZ = vars.q1_pZ;
//   vars.jd_q1_pT = vars.q1_pT;
//   vars.jd_q1_E   = vars.q1_E;
//   vars.jd_q1_Eta = vars.q1_Eta;
//   vars.jd_q1_Phi = vars.q1_Phi;
// 
//   vars.jd_q2_pX = vars.q2_pX;
//   vars.jd_q2_pY = vars.q2_pY;
//   vars.jd_q2_pZ = vars.q2_pZ;
//   vars.jd_q2_pT = vars.q2_pT;
//   vars.jd_q2_E   = vars.q2_E;
//   vars.jd_q2_Eta = vars.q2_Eta;
//   vars.jd_q2_Phi = vars.q2_Phi;
   
//  }
//  else { 
 
  ///-- alpha * q1 + beta * q2 = -Z   
    
   double delta =  temp_q1.X() * temp_q2.Y()  -  temp_q1.Y() * temp_q2.X() ;
   double alpha =  (totalZ.X() * temp_q2.Y()  -  totalZ.Y() * temp_q2.X()) / delta;
   double beta   =  (totalZ.Y() * temp_q1.X()  -  totalZ.X() * temp_q1.Y()) / delta;

   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > new_q1 = alpha * temp_q1;
   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > new_q2 =   beta * temp_q2;   

   vars.jd_Delta_q1_P = alpha ;
   vars.jd_Delta_q2_P = beta ;

   vars.jd_q1_pX = new_q1.X() ;
   vars.jd_q1_pY = new_q1.Y();
   vars.jd_q1_pZ = new_q1.Z();
   vars.jd_q1_pT = new_q1.Rho();
   vars.jd_q1_E   = sqrt(vars.q1_E * vars.q1_E - temp_q1.Pt() * temp_q1.Pt() + new_q1.Rho() * new_q1.Rho());
   vars.jd_q1_Eta = vars.q1_Eta;
   vars.jd_q1_Phi = vars.q1_Phi;
 
   vars.jd_q2_pX = new_q2.X() ;
   vars.jd_q2_pY = new_q2.Y();
   vars.jd_q2_pZ = new_q2.Z();
   vars.jd_q2_pT = new_q2.Rho();
   vars.jd_q2_E   = sqrt(vars.q2_E * vars.q2_E - temp_q2.Pt() * temp_q2.Pt() + new_q2.Rho() * new_q2.Rho());
   vars.jd_q2_Eta = vars.q2_Eta;
   vars.jd_q2_Phi = vars.q2_Phi;

   vars.jd_met_X = vars.met_X + vars.q1_pX - vars.jd_q1_pX + vars.q2_pX - vars.jd_q2_pX;
   vars.jd_met_Y = vars.met_Y + vars.q1_pY - vars.jd_q1_pY + vars.q2_pY - vars.jd_q2_pY;

   vars.jd_met = sqrt ( vars.jd_met_X * vars.jd_met_X + vars.jd_met_Y * vars.jd_met_Y);

    vars.jd_M_qq = ( new_q1 + new_q2 ).mass();


   //---- project-jd-met ---- 
  double Lep1_phi = 100. ;
  if (FlavourLep1 == 11 ) Lep1_phi = reader.Get4V ("mcEle")->at (iLep1).Phi () ;
  else Lep1_phi = reader.Get4V ("mcMu")->at (iLep1).Phi () ;
  double Lep2_phi = 100. ;
  if (FlavourLep2 == 11 ) Lep2_phi = reader.Get4V ("mcEle")->at (iLep2).Phi () ;
  else Lep2_phi = reader.Get4V ("mcMu")->at (iLep2).Phi () ;

  double MET_phi = atan2(vars.jd_met_Y,vars.jd_met_X) ;
  double deltaPhiMin = std::min (deltaPhi (MET_phi, Lep1_phi), deltaPhi (MET_phi, Lep2_phi)) ;

 if (deltaPhiMin < 1.57079632679) vars.jd_pmet = vars.jd_met * sin (deltaPhiMin) ;
 else vars.jd_pmet = vars.jd_met;

//  }
}
 
///====== Find third jet of the event

void FindAddJet (treeReader& reader,const int& q1,const int& q2, const std::vector<int>* blacklistJet_forCJV,
		 std::vector<ROOT::Math::XYZTVector> & Jet_Candidate, const double& EtMin){
  
  for(int iJet=0; iJet<reader.Get4V("mcJet")->size() ; iJet ++){
      if (iJet==q1 || iJet==q2) continue;
      
      if (reader.Get4V("mcJet")->at(iJet).Et() < EtMin) continue;
       
      bool skipJet=false;
 
     if(blacklistJet_forCJV){
      for(unsigned int kk = 0; kk < blacklistJet_forCJV->size(); ++kk){
        if(blacklistJet_forCJV->at(kk) == static_cast<int>(iJet)) skipJet = true;
       } 
     }
    
    if(skipJet) continue;
 
    Jet_Candidate.push_back(reader.Get4V("mcJet")->at(iJet));
  }
  
  if(Jet_Candidate.size()!=0) sort(Jet_Candidate.begin(),Jet_Candidate.end(),pT_sort());

}
      
///======================= Set Third Jet variables


void SetThirdJetVariables(Variables_Gen& vars, std::vector<ROOT::Math::XYZTVector> & Jet_Candidate)
{

 if(Jet_Candidate.size()!=0){
   vars.q3_pX = Jet_Candidate.at(0).X();
   vars.q3_pY = Jet_Candidate.at(0).Y();
   vars.q3_pZ = Jet_Candidate.at(0).Z();
   vars.q3_pT = Jet_Candidate.at(0).pt();
   vars.q3_E = Jet_Candidate.at(0).E();
   vars.q3_Eta = Jet_Candidate.at(0).Eta();
   vars.q3_Phi = Jet_Candidate.at(0).Phi();
 }
  
}


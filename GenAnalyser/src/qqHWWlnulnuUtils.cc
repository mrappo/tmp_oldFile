/*
qqHWWlnulnuUtils
*/

#include "qqHWWlnulnuUtils.h"



///==== MC Decay Channel of VV ====
std::pair<int,int> GetMCDecayChannel(const std::vector<float>& pdgId){
 
 std::pair<int,int> result;
 result.first  = -1;
 result.second = -1;
 
 if (pdgId.size() == 4){
  ///=================   Vf
  ///== pdgId.at(0) ==> f11 
  ///== pdgId.at(1) ==> f12 (V = 1 & f = 2)
  ///== pdgId.at(2) ==> f21
  ///== pdgId.at(3) ==> f22
  
  int f11 = abs(pdgId.at(0));
  int f12 = abs(pdgId.at(1));
  int f21 = abs(pdgId.at(2));
  int f22 = abs(pdgId.at(3));
  
  
  ///==== lepton - lepton ====
  if (f11 >= 11 && f11 <= 16 && f21 >= 11 && f21 <= 16) {
   result.first  = 1;
   ///==== emu - emu ====
   if (f11 >= 11 && f11 <= 14 && f21 >= 11 && f21 <= 14) {
    result.second  = 1;
   }
   ///==== emu - tau ====
   if (f11 >= 11 && f11 <= 14 && f21 >=15 && f21 <= 16) {
    result.second  = 2;
   }
   ///==== tau - emu ====
   if (f21 >= 11 && f21 <= 14 && f11 >= 15 && f11 <= 16) {
    result.second  = 3;    
   }
   ///==== tau - tau ====
   if (f21 >= 15 && f21 <= 16 && f11 >= 15 && f11 <= 16) {
    result.second  = 4;
   }
  }
  
  ///==== lepton - quark ====
  if ((f11 >= 11 && f11 <= 16 && f21 >= 1 && f21 <= 4) || (f21 >= 11 && f21 <= 16 && f11 >= 1 && f11 <= 4)) {
   result.first  = 2;
   ///==== emu - quark ====
   if (f11 >= 11 && f11 <= 14 && f21 >= 1 && f21 <= 4){
    result.second  = 1;
   }
   ///==== quark - emu ====
   if (f21 >= 11 && f21 <= 14 && f11 >= 1 && f11 <= 4) {
    result.second  = 2;
   }
   ///==== tau - quark ====
   if (f21 >= 1 && f21 <= 4 && f11 >= 15 && f11 <= 16){
    result.second  = 3;
   }
   ///==== quark -tau ====
   if (f11 >= 1 && f11 <= 4 && f21 >=15 && f21 <= 16) {
    result.second  = 4;
   }   
  }
  
  ///==== quark - quark ====
  if (f11 >= 1 && f11 <= 4 && f21 >= 1 && f21 <= 4) {
   result.first  = 3;
   result.second  = -1;   
  }
 }
 
 return result;
}




std::pair<int,int> GetMCDecayChannel(float pdgId11,float pdgId12, float pdgId21, float pdgId22){
 std::vector<float> vectPdgId;
 vectPdgId.push_back(pdgId11);
 vectPdgId.push_back(pdgId12);
 vectPdgId.push_back(pdgId21);
 vectPdgId.push_back(pdgId22);
 return GetMCDecayChannel(vectPdgId);
}










///==== get GoodCombination MVA ====

std::pair<double,int> GetCombinationMVA(treeReader& reader, std::vector<std::vector<int> >& combinations,TMVA::Reader* TMVAreader, TString& methodName, Float_t* input_variables, std::vector<int>* whitelistJet){
 
 Float_t pT_RECO_q1;
 Float_t pT_RECO_q2;
 Float_t pT_RECO_b1;
 Float_t pT_RECO_b2;
 
 Float_t DR_RECO_qq;
 Float_t DR_RECO_qb;
 Float_t DR_RECO_bb;
 
 Float_t jets_trackCountingHighEffBJetTags_RECO_q1;
 Float_t jets_trackCountingHighEffBJetTags_RECO_b1;
 Float_t jets_trackCountingHighPurBJetTags_RECO_q1;
 Float_t jets_trackCountingHighPurBJetTags_RECO_b1;
 Float_t jets_simpleSecondaryVertexBJetTags_RECO_q1;
 Float_t jets_simpleSecondaryVertexBJetTags_RECO_b1;
 Float_t jets_combinedSecondaryVertexBJetTags_RECO_q1;
 Float_t jets_combinedSecondaryVertexBJetTags_RECO_b1;
 Float_t jets_combinedSecondaryVertexMVABJetTags_RECO_q1;
 Float_t jets_combinedSecondaryVertexMVABJetTags_RECO_b1;
 
 Float_t jets_trackCountingHighEffBJetTags_RECO_q2;
 Float_t jets_trackCountingHighEffBJetTags_RECO_b2;
 Float_t jets_trackCountingHighPurBJetTags_RECO_q2;
 Float_t jets_trackCountingHighPurBJetTags_RECO_b2;
 Float_t jets_simpleSecondaryVertexBJetTags_RECO_q2;
 Float_t jets_simpleSecondaryVertexBJetTags_RECO_b2;
 Float_t jets_combinedSecondaryVertexBJetTags_RECO_q2;
 Float_t jets_combinedSecondaryVertexBJetTags_RECO_b2;
 Float_t jets_combinedSecondaryVertexMVABJetTags_RECO_q2;
 Float_t jets_combinedSecondaryVertexMVABJetTags_RECO_b2; 
 
 
//  std::vector<ROOT::Math::XYZTVector>* muons = reader.Get4V("muons");
//  std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");
 
 std::vector<ROOT::Math::XYZTVector>* jets = reader.Get4V("jets");
 std::vector<float>* jets_trackCountingHighEffBJetTags = reader.GetFloat("jets_trackCountingHighEffBJetTags");
 std::vector<float>* jets_trackCountingHighEffBJetTagsDR = reader.GetFloat("jets_trackCountingHighEffBJetTagsDR");
 std::vector<float>* jets_trackCountingHighPurBJetTags = reader.GetFloat("jets_trackCountingHighPurBJetTags");
 std::vector<float>* jets_simpleSecondaryVertexBJetTags = reader.GetFloat("jets_simpleSecondaryVertexBJetTags");
 std::vector<float>* jets_combinedSecondaryVertexBJetTags = reader.GetFloat("jets_combinedSecondaryVertexBJetTags");
 std::vector<float>* jets_combinedSecondaryVertexMVABJetTags = reader.GetFloat("jets_combinedSecondaryVertexMVABJetTags");
 
 double TMVA_value = -100000000;
 int MVACombination = -1;
 
 int nComb = combinations.size();
 
 for (int iComb = 0; iComb<nComb; iComb++){
  int q1 = combinations.at(iComb).at(0);
  int q2 = combinations.at(iComb).at(1);
  int b1 = combinations.at(iComb).at(2);
  int b2 = combinations.at(iComb).at(3);    
  
  if(whitelistJet != NULL)
   if( (whitelistJet -> at(q1)) != 1 || (whitelistJet -> at(q2)) != 1 || (whitelistJet -> at(b1)) != 1 || (whitelistJet -> at(b2)) != 1)
    continue;
     
  pT_RECO_q1 = jets->at(q1).Pt();
  pT_RECO_q2 = jets->at(q2).Pt();
  pT_RECO_b1 = jets->at(b1).Pt();
  pT_RECO_b2 = jets->at(b2).Pt();
  DR_RECO_qq = ROOT::Math::VectorUtil::DeltaR(jets->at(q1),jets->at(q2));
  DR_RECO_bb = ROOT::Math::VectorUtil::DeltaR(jets->at(b1),jets->at(b2));
  DR_RECO_qb = std::min(ROOT::Math::VectorUtil::DeltaR(jets->at(q1),jets->at(b2)),ROOT::Math::VectorUtil::DeltaR(jets->at(b1),jets->at(q2)));
  
  if (jets->at(q1).Pt() > jets->at(q2).Pt()) {
   jets_trackCountingHighEffBJetTags_RECO_q1 = jets_trackCountingHighEffBJetTags->at(q1);
   jets_trackCountingHighEffBJetTags_RECO_q2 = jets_trackCountingHighEffBJetTags->at(q2);  
   jets_trackCountingHighPurBJetTags_RECO_q1 = jets_trackCountingHighPurBJetTags->at(q1);
   jets_trackCountingHighPurBJetTags_RECO_q2 = jets_trackCountingHighPurBJetTags->at(q2);
   jets_simpleSecondaryVertexBJetTags_RECO_q1 = jets_simpleSecondaryVertexBJetTags->at(q1);
   jets_simpleSecondaryVertexBJetTags_RECO_q2 = jets_simpleSecondaryVertexBJetTags->at(q2);
   jets_combinedSecondaryVertexBJetTags_RECO_q1 = jets_combinedSecondaryVertexBJetTags->at(q1);
   jets_combinedSecondaryVertexBJetTags_RECO_q2 = jets_combinedSecondaryVertexBJetTags->at(q2);
   jets_combinedSecondaryVertexMVABJetTags_RECO_q1 = jets_combinedSecondaryVertexMVABJetTags->at(q1);
   jets_combinedSecondaryVertexMVABJetTags_RECO_q2 = jets_combinedSecondaryVertexMVABJetTags->at(q2);
  }
  else  {
   jets_trackCountingHighEffBJetTags_RECO_q1 = jets_trackCountingHighEffBJetTags->at(q2);
   jets_trackCountingHighEffBJetTags_RECO_q2 = jets_trackCountingHighEffBJetTags->at(q1);  
   jets_trackCountingHighPurBJetTags_RECO_q1 = jets_trackCountingHighPurBJetTags->at(q2);
   jets_trackCountingHighPurBJetTags_RECO_q2 = jets_trackCountingHighPurBJetTags->at(q1);
   jets_simpleSecondaryVertexBJetTags_RECO_q1 = jets_simpleSecondaryVertexBJetTags->at(q2);
   jets_simpleSecondaryVertexBJetTags_RECO_q2 = jets_simpleSecondaryVertexBJetTags->at(q1);
   jets_combinedSecondaryVertexBJetTags_RECO_q1 = jets_combinedSecondaryVertexBJetTags->at(q2);
   jets_combinedSecondaryVertexBJetTags_RECO_q2 = jets_combinedSecondaryVertexBJetTags->at(q1);
   jets_combinedSecondaryVertexMVABJetTags_RECO_q1 = jets_combinedSecondaryVertexMVABJetTags->at(q2);
   jets_combinedSecondaryVertexMVABJetTags_RECO_q2 = jets_combinedSecondaryVertexMVABJetTags->at(q1);
  }
  
  if (jets->at(b1).Pt() > jets->at(b2).Pt()) {
   jets_trackCountingHighEffBJetTags_RECO_b1 = jets_trackCountingHighEffBJetTags->at(b1);
   jets_trackCountingHighEffBJetTags_RECO_b2 = jets_trackCountingHighEffBJetTags->at(b2);  
   jets_trackCountingHighPurBJetTags_RECO_b1 = jets_trackCountingHighPurBJetTags->at(b1);
   jets_trackCountingHighPurBJetTags_RECO_b2 = jets_trackCountingHighPurBJetTags->at(b2);
   jets_simpleSecondaryVertexBJetTags_RECO_b1 = jets_simpleSecondaryVertexBJetTags->at(b1);
   jets_simpleSecondaryVertexBJetTags_RECO_b2 = jets_simpleSecondaryVertexBJetTags->at(b2);
   jets_combinedSecondaryVertexBJetTags_RECO_b1 = jets_combinedSecondaryVertexBJetTags->at(b1);
   jets_combinedSecondaryVertexBJetTags_RECO_b2 = jets_combinedSecondaryVertexBJetTags->at(b2);
   jets_combinedSecondaryVertexMVABJetTags_RECO_b1 = jets_combinedSecondaryVertexMVABJetTags->at(b1);
   jets_combinedSecondaryVertexMVABJetTags_RECO_b2 = jets_combinedSecondaryVertexMVABJetTags->at(b2);
  }
  else  {
   jets_trackCountingHighEffBJetTags_RECO_b1 = jets_trackCountingHighEffBJetTags->at(b2);
   jets_trackCountingHighEffBJetTags_RECO_b2 = jets_trackCountingHighEffBJetTags->at(b1);  
   jets_trackCountingHighPurBJetTags_RECO_b1 = jets_trackCountingHighPurBJetTags->at(b2);
   jets_trackCountingHighPurBJetTags_RECO_b2 = jets_trackCountingHighPurBJetTags->at(b1);
   jets_simpleSecondaryVertexBJetTags_RECO_b1 = jets_simpleSecondaryVertexBJetTags->at(b2);
   jets_simpleSecondaryVertexBJetTags_RECO_b2 = jets_simpleSecondaryVertexBJetTags->at(b1);
   jets_combinedSecondaryVertexBJetTags_RECO_b1 = jets_combinedSecondaryVertexBJetTags->at(b2);
   jets_combinedSecondaryVertexBJetTags_RECO_b2 = jets_combinedSecondaryVertexBJetTags->at(b1);
   jets_combinedSecondaryVertexMVABJetTags_RECO_b1 = jets_combinedSecondaryVertexMVABJetTags->at(b2);
   jets_combinedSecondaryVertexMVABJetTags_RECO_b2 = jets_combinedSecondaryVertexMVABJetTags->at(b1);
  }
  
  input_variables[0] = pT_RECO_q1;
  input_variables[1] = pT_RECO_q2;
  input_variables[2] = pT_RECO_b1;
  input_variables[3] = pT_RECO_b2;
  input_variables[4] = DR_RECO_qq;
  input_variables[5] = DR_RECO_qb;
  input_variables[6] = DR_RECO_bb;
  input_variables[7] = jets_trackCountingHighEffBJetTags_RECO_q1;
  input_variables[8] = jets_trackCountingHighEffBJetTags_RECO_b1;
  input_variables[9] = jets_trackCountingHighEffBJetTags_RECO_q2;
  input_variables[11] = jets_trackCountingHighEffBJetTags_RECO_b2;
  
  double MVA_temp = TMVAreader->EvaluateMVA(methodName); 
  if (MVA_temp > TMVA_value) {
   TMVA_value = MVA_temp;
   MVACombination = iComb;
  }
 }
 
 
 
 std::pair<double,int> result;
 result.first = TMVA_value;
 result.second = MVACombination;
 
 return result;
}


///==== return the value of the Min Chi2 and the corresponding number in the combinations ====
std::pair<double,int> GetCombinationChi2(treeReader& reader, std::vector<std::vector<int> >& combinations,std::vector<int>* whitelistJet,const ROOT::Math::XYZTVector& lepton, const ROOT::Math::XYZTVector& MET){
 
 double Chi2 = -1;
 double Chi2Combination = -1;
 
 std::vector<ROOT::Math::XYZTVector>* jets = reader.Get4V("jets");
 std::vector<float>* jets_trackCountingHighEffBJetTags = reader.GetFloat("jets_trackCountingHighEffBJetTags");
 std::vector<float>* jets_trackCountingHighEffBJetTagsDR = reader.GetFloat("jets_trackCountingHighEffBJetTagsDR");
 std::vector<float>* jets_trackCountingHighPurBJetTags = reader.GetFloat("jets_trackCountingHighPurBJetTags");
 std::vector<float>* jets_simpleSecondaryVertexBJetTags = reader.GetFloat("jets_simpleSecondaryVertexBJetTags");
 std::vector<float>* jets_combinedSecondaryVertexBJetTags = reader.GetFloat("jets_combinedSecondaryVertexBJetTags");
 std::vector<float>* jets_combinedSecondaryVertexMVABJetTags = reader.GetFloat("jets_combinedSecondaryVertexMVABJetTags");
 
 
 int nComb = combinations.size();
 
 for (int iComb = 0; iComb<nComb; iComb++){
  int q1 = combinations.at(iComb).at(0);
  int q2 = combinations.at(iComb).at(1);
  int b1 = combinations.at(iComb).at(2); //---> in this analysis this is the bjet that couple with qq
  int b2 = combinations.at(iComb).at(3);    
  
  if(whitelistJet != NULL)
   if( (whitelistJet -> at(q1)) != 1 || (whitelistJet -> at(q2)) != 1 || (whitelistJet -> at(b1)) != 1 || (whitelistJet -> at(b2)) != 1)
    continue;
  
  double sigma1 = 1;
  double sigma2 = 1;
  double sigma3 = 1;
  
  double massW = 80.398;
  double massTop = 170; //--------- ??????????????????????????????????????????????????
  ROOT::Math::XYZTVector neutrino; // = GetNeutrino(MET,lepton,massW);
  
  double Chi2_temp = 0;
  Chi2_temp += ((jets->at(q1) + jets->at(q2)).M() - massW) * ((jets->at(q1) + jets->at(q2)).M() - massW)  / sigma1 / sigma1; 
  
  Chi2_temp += ((jets->at(q1) + jets->at(q2) + jets->at(b1)).M() - massTop) * ((jets->at(q1) + jets->at(q2)).M() - massTop)  / sigma1 / sigma1; 
  
  Chi2_temp += ((jets->at(b2) + lepton + neutrino).M() - massTop) * ((jets->at(q1) + jets->at(q2)).M() - massTop)  / sigma1 / sigma1; 
    
  if (Chi2_temp < Chi2 || Chi2 == -1) {
   Chi2 = Chi2_temp;
   Chi2Combination = iComb;
  }
 }
 
 std::pair<double,int> result;
 result.first = Chi2;
 result.second = Chi2Combination;
 return result;
}


///==== get Selected Lepton: electron or muon ====
int SelectLeptonTTBar(std::vector<ROOT::Math::XYZTVector>& leptons,const std::string& method,const double& ptMin,const double& etaMax,int& numLeptons,const std::vector<int>* blacklist){
 //---- initialize variable with result
 int it = -1;
 numLeptons = 0;
 // initialize the selection variable
 double maxPt = -999999.;
 double tempPt = 0.;
 //---- loop over leptons
 for(unsigned int i = 0; i < leptons.size(); ++i){
  if( sqrt(leptons.at(i).Perp2()) < ptMin ) continue;
  if( fabs(leptons.at(i).eta()) > etaMax ) continue;
  bool skipLep = false;
  if(blacklist){
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk){
    if(blacklist -> at(kk) == static_cast<int>(i)) skipLep = true;
   }
  }
  if(skipLep) continue;
  numLeptons ++; ///---- number of leptons that survive selections ----
  //---- select leptons with different techniques
  if(method == "maxPt")
  {
   tempPt = sqrt(leptons.at(i).perp2());
   if(tempPt > maxPt)
   {
    maxPt = tempPt;
    it = i;
   }
  }
  //-------------------------------------
 }//---- loop over leptons 
 if(method == "maxPt")
  return it;
 else return -1;
}


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
		    ){
 //---- initialize variable with result
 int it = -1;
 numLeptons = 0;
 // initialize the selection variable
 double maxPt = -999999.;
 double tempPt = 0.;
 //---- loop over leptons
 for(unsigned int i = 0; i < leptons.size(); ++i){
  if( sqrt(leptons.at(i).Perp2()) < ptMin ) continue;
  if( fabs(leptons.at(i).eta()) > etaMax ) continue;
  if( (tkIso.at(i) + emIso.at(i) + hadIso.at(i)) / sqrt(leptons.at(i).Perp2()) > isoMax ) continue; 
//   if( mud0.at(i) > d0Max ) continue; 
  bool skipLep = false;
  if(blacklist){
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk){
    if(blacklist -> at(kk) == static_cast<int>(i)) skipLep = true;
   }
  }
  if(skipLep) continue;
  numLeptons ++; ///---- number of leptons that survive selections ----
  //---- select leptons with different techniques
  if(method == "maxPt")
  {
   tempPt = sqrt(leptons.at(i).perp2());
   if(tempPt > maxPt)
   {
    maxPt = tempPt;
    it = i;
   }
  }
  //-------------------------------------
 }//---- loop over leptons 
 if(method == "maxPt")
  return it;
 else return -1;
}




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
		    ){
 //---- initialize variable with result
 int it = -1;
 numLeptons = 0;
 // initialize the selection variable
 double maxPt = -999999.;
 double tempPt = 0.;
 //---- loop over leptons
 for(unsigned int i = 0; i < leptons.size(); ++i){
  if( sqrt(leptons.at(i).Perp2()) < ptMin ) continue;
  if( fabs(leptons.at(i).eta()) > etaMax ) continue;
  
  ///---- H-like isolation ----
  if (((tkIso.at(i) + emIso.at(i) + hadIso.at(i)) < 8.) && (sqrt(leptons.at(i).Perp2()) < (10. + 15./8.*(tkIso.at(i) + emIso.at(i) + hadIso.at(i))))) continue;
  if ((tkIso.at(i) + emIso.at(i) + hadIso.at(i)) > 8.) continue;
  
//   if( (tkIso.at(i) + emIso.at(i) + hadIso.at(i)) / sqrt(leptons.at(i).Perp2()) > isoMax ) continue; 
  
  if( eleId.at(i) == 0 ) continue;
  if( eled0.at(i) > d0Max ) continue; 
  bool skipLep = false;
  if(blacklist){
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk){
    if(blacklist -> at(kk) == static_cast<int>(i)) skipLep = true;
   }
  }
  if(skipLep) continue;
  numLeptons ++; ///---- number of leptons that survive selections ----
  //---- select leptons with different techniques
  if(method == "maxPt")
  {
   tempPt = sqrt(leptons.at(i).perp2());
   if(tempPt > maxPt)
   {
    maxPt = tempPt;
    it = i;
   }
  }
  //-------------------------------------
 }//---- loop over leptons 
 if(method == "maxPt")
  return it;
 else return -1;
}



///==== get num Jets ====
int getNumJets(const std::vector<int>& whitelist){
 int result = 0;
 for (int i=0; i<whitelist.size(); i++) {
  result += whitelist.at(i);
 }
 return result;
}



///==== CopyTree ====
void CopyTree(){
  
}








///==== get Jet Energy Correction ====
double getJEC(const ROOT::Math::XYZTVector& Jet){
 double pT = Jet.Pt();
 double result = exp(-3.32898e-01  -2.34368e-02 * pT) + 7.40976e-01;
//  expo(3) is a substitute for exp([3]+[4]*x)
 return result;
}














///=================================== VBF =======================================

///==== get best Combination Jet ID with MVA ====

std::pair<double,int> GetCombination_Jets_ID_MVA(std::vector<ROOT::Math::XYZTVector> &jets, std::vector<std::vector<int> >& combinations,TMVA::Reader* TMVAreader, TString& methodName, Float_t* input_variables, std::vector<int>* whitelistJet){
 
 Float_t pT_RECO_q1;
 Float_t pT_RECO_q2;
 Float_t eta_RECO_q1;
 Float_t eta_RECO_q2;
 Float_t eta_RECO_q1_eta_RECO_q2;
 Float_t Deta_RECO_q12;
 Float_t Mjj;
 
//  std::vector<ROOT::Math::XYZTVector>* jets = reader.Get4V("jets");
  
 double TMVA_value = -100000000;
 int MVACombination = -1;
 
 int nComb = combinations.size();
 
 for (int iComb = 0; iComb<nComb; iComb++){
  int q1 = combinations.at(iComb).at(0);
  int q2 = combinations.at(iComb).at(1);
  
  if(whitelistJet != NULL)
   if((whitelistJet -> at(q1)) != 1 || (whitelistJet -> at(q2)) != 1) continue;
   
   if (jets.at(q1).Pt() < jets.at(q2).Pt()) {
    int tempq = q1;
    q1 = q2;
    q2 = tempq;
   }
   
   pT_RECO_q1 = jets.at(q1).Pt();
  pT_RECO_q2 = jets.at(q2).Pt();
  
  eta_RECO_q1 = jets.at(q1).Eta();
  eta_RECO_q2 = jets.at(q2).Eta();
  
  eta_RECO_q1_eta_RECO_q2 = eta_RECO_q1 * eta_RECO_q2;
  Deta_RECO_q12 = fabs(eta_RECO_q1-eta_RECO_q2);
  
  Mjj = (jets.at(q1) + jets.at(q2)).M();
  
  
  input_variables[0] = pT_RECO_q1;
  input_variables[1] = pT_RECO_q2;
  input_variables[2] = eta_RECO_q1;
  input_variables[3] = eta_RECO_q2;
  input_variables[4] = eta_RECO_q1_eta_RECO_q2;
  input_variables[5] = Deta_RECO_q12;
  input_variables[6] = Mjj;
  
  double MVA_temp = TMVAreader->EvaluateMVA(methodName); 
  if (MVA_temp > TMVA_value) {
   TMVA_value = MVA_temp;
   MVACombination = iComb;
  }
 }
 
 std::pair<double,int> result;
 result.first = TMVA_value;
 result.second = MVACombination;
 
 return result;
}





///==== Plot tools ====
///
///    _ \   |         |   
///   |   |  |   _ \   __| 
///   ___/   |  (   |  |   
///  _|     _| \___/  \__| 
///                        



///==== GetTrendInfo ====
///==== Transform TH1 with "trace" information to TH1 ====

TH1F* GetTrendInfo(TH1F* hTrend, double min, double max){
 int nbin = hTrend->GetNbinsX();
//  double maxX = hTrend->GetYaxis()->GetXmax();
//  double minX = hTrend->GetYaxis()->GetXmin();
//  std::cout << "minX = " << minX << " ::: " << maxX << std::endl;
//  std::cout << "min  = " << min  << " ::: " << max  << std::endl;
//  maxX = std::min(max,maxX);
//  minX = std::max(min,minX);
//  std::cout << "minX = " << minX << " ::: " << maxX << std::endl;
 double maxX = max;
 double minX = min;
 
 std::string name = hTrend->GetName(); 
 std::string nameNew;
 nameNew = name + "_Info";
 
 TH1F* hTrendInfo = new TH1F (nameNew.c_str(),nameNew.c_str(),nbin,minX,maxX);
 for (int iBin = 0; iBin<nbin; iBin++){
  hTrendInfo->Fill(hTrend->GetBinContent(iBin+1));
 }
 
 hTrendInfo->SetLineColor(hTrend->GetLineColor());
 hTrendInfo->SetLineWidth(2);
 hTrendInfo->SetMarkerColor(hTrend->GetMarkerColor());
 hTrendInfo->SetMarkerStyle(20);
 hTrendInfo->SetMarkerSize(1);
 hTrendInfo->SetFillStyle(1001);
 hTrendInfo->GetXaxis()->SetTitle(hTrend->GetYaxis()->GetTitle()); 
 return hTrendInfo;
}



///==== Pull Plot ====
std::pair<TGraphErrors*, TGraphErrors*> grPullPlot(TH1F* hDATA, TH1F* hMC){
 int nbin = hDATA->GetNbinsX();
 double max = hDATA->GetXaxis()->GetXmax();
 double min = hDATA->GetXaxis()->GetXmin();
 std::string name1 = hDATA->GetName();
 std::string name2 = hMC->GetName(); 
 std::string nameNew = name1 + name2; 
 
 TGraphErrors* grPool = new TGraphErrors();
 TGraphErrors* grPoolMC = new TGraphErrors();
 int point = 0; 
 for (int iBin = 0; iBin<nbin; iBin++){
  double X = hDATA->GetBinCenter(iBin+1);
  double DATA = hDATA->GetBinContent(iBin+1);
  double MC = hMC->GetBinContent(iBin+1);
  double errMC = hMC->GetBinError(iBin+1);
  
  grPool->SetPoint      (point, X, (MC ? DATA/MC : 0));
  grPool->SetPointError (point, 0, (MC ? sqrt(DATA)/MC : 0));
  
  grPoolMC->SetPoint      (point, X, 1.);
  grPoolMC->SetPointError (point, 0, (MC ? errMC / MC : 0));
  point++;
 }
 grPool->SetLineColor(kRed);
 grPool->SetLineWidth(2);
 grPool->SetMarkerColor(kRed);
 grPool->SetMarkerStyle(20);
 grPool->SetMarkerSize(1);
 grPool->GetXaxis()->SetTitle(hMC->GetXaxis()->GetTitle()); 
 grPool->GetYaxis()->SetTitle("DATA/MC"); 
 grPool->GetYaxis()->SetRangeUser(0.,3.); 
 
 grPoolMC->SetLineWidth(0);
 grPoolMC->SetMarkerSize(0);
 grPoolMC->SetFillColor(kGreen);
 grPoolMC->GetXaxis()->SetTitle(hMC->GetXaxis()->GetTitle()); 
 grPoolMC->GetYaxis()->SetTitle("DATA/MC"); 
 grPoolMC->GetYaxis()->SetRangeUser(0.,3.); 
 
 return std::pair<TGraphErrors*, TGraphErrors*> (grPool, grPoolMC); 
}


//==== with TH1 ====

TH1F* PullPlot(TH1F* hDATA, TH1F* hMC){
 int nbin = hDATA->GetNbinsX();
 double max = hDATA->GetXaxis()->GetXmax();
 double min = hDATA->GetXaxis()->GetXmin();
 std::string name1 = hDATA->GetName();
 std::string name2 = hMC->GetName(); 
 std::string nameNew = name1 + name2; 
 TH1F* hPool = new TH1F (nameNew.c_str(),nameNew.c_str(),nbin,min,max);
 for (int iBin = 0; iBin<nbin; iBin++){
  double A = hDATA->GetBinContent(iBin+1);
  double B = hMC->GetBinContent(iBin+1);
  if (A+B != 0) {
   hPool->SetBinContent(iBin,(A-B)/(A+B)*2.);
//    hPool->SetBinError(iBin,4. * A / (A+B) / (A+B) * sqrt(A+B));
   hPool->SetBinError(iBin,4. * B / (A+B) / (A+B) * sqrt(A));
  }
 }
 hPool->SetLineColor(kRed);
 hPool->SetLineWidth(2);
 hPool->SetMarkerColor(kRed);
 hPool->SetMarkerStyle(20);
 hPool->SetMarkerSize(1);
 hPool->GetXaxis()->SetTitle(hMC->GetXaxis()->GetTitle()); 
 hPool->GetYaxis()->SetTitle("2(DATA-MC)/(DATA+MC)"); 
 hPool->GetYaxis()->SetRangeUser(-2.,2.); 
 return hPool;
}


///==== Pull Plot: drawing utility ====
void PullPlot(TCanvas* canvas, TH1* hDATA, TH1* hMC){
 canvas->Divide(1,2);
 int nbin = hDATA->GetNbinsX();
 double max = hDATA->GetXaxis()->GetXmax();
 double min = hDATA->GetXaxis()->GetXmin();
 std::string name1 = hDATA->GetName();
 std::string name2 = hMC->GetName();
 
 std::string nameNew = name1 + name2;
 
 TH1F* hPool = new TH1F (nameNew.c_str(),nameNew.c_str(),nbin,min,max);
 for (int iBin = 0; iBin<nbin; iBin++){
  double A = hDATA->GetBinContent(iBin);
  double B = hMC->GetBinContent(iBin);
  if (A+B != 0) {
   hPool->SetBinContent(iBin,(A-B)/(A+B)*2.);
   hPool->SetBinError(iBin,4. * B / (A+B) / (A+B) * sqrt(A));
  }
 }
 
 double maxY_DATA = hDATA->GetMaximum();
 double maxY_MC = hMC->GetMaximum();
 
 double minY_DATA = hDATA->GetMinimum();
 double minY_MC = hMC->GetMinimum();
 
 canvas->cd(1); 
 if (maxY_MC > maxY_DATA) {
  hMC->Draw();
  hDATA->Draw("EsameP");
 }
 else {
  hDATA->Draw("EP");
  hMC->Draw("same");  
 }
 hDATA->GetYaxis()->SetRangeUser(std::min(minY_MC,minY_DATA),std::max(maxY_MC,maxY_DATA) * 1.1);
 hMC  ->GetYaxis()->SetRangeUser(std::min(minY_MC,minY_DATA),std::max(maxY_MC,maxY_DATA) * 1.1);
 
 gPad->SetGrid();
 
 canvas->cd(2); 
 hPool->SetLineColor(kRed);
 hPool->SetLineWidth(2);
 hPool->SetMarkerColor(kRed);
 hPool->SetMarkerStyle(20);
 hPool->SetMarkerSize(1);
 hPool->Draw("EP");
 hPool->GetXaxis()->SetTitle(hMC->GetXaxis()->GetTitle()); 
 hPool->GetYaxis()->SetTitle("2(DATA-MC)/(DATA+MC)"); 
 hPool->Draw("EP");
 gPad->SetGrid();
}


void PullPlot(TCanvas* canvas, TH1* hDATA, THStack* hsMC){
 canvas->Divide(1,2);
 int nbin = hDATA->GetNbinsX();
 double max = hDATA->GetXaxis()->GetXmax();
 double min = hDATA->GetXaxis()->GetXmin();
 std::string name1 = hDATA->GetName();
 std::string name2 = hsMC->GetName();
 
 std::string nameNew = name1 + name2;
 
 TH1F* hPool = new TH1F (nameNew.c_str(),nameNew.c_str(),nbin,min,max);
 for (int iBin = 0; iBin<nbin; iBin++){
  double A = hDATA->GetBinContent(iBin);
  double B = ((TH1*)(hsMC->GetStack()->Last()))->GetBinContent(iBin);
  if (A+B != 0) {
   hPool->SetBinContent(iBin,(A-B)/(A+B)*2.);
   hPool->SetBinError(iBin,4. * B / (A+B) / (A+B) * sqrt(A));
  }
 }
 
 double maxY_DATA = hDATA->GetMaximum();
 double maxY_MC = ((TH1*)(hsMC->GetStack()->Last()))->GetMaximum();
 
 double minY_DATA = hDATA->GetMinimum();
 double minY_MC = ((TH1*)(hsMC->GetStack()->Last()))->GetMinimum();
 
 canvas->cd(1); 
 hDATA->Draw("EP");
 DrawStack(hsMC);
 hDATA->Draw("EsameP");
 gPad->SetGrid();
 
 canvas->cd(2); 
 hPool->SetLineColor(kRed);
 hPool->SetLineWidth(2);
 hPool->SetMarkerColor(kRed);
 hPool->SetMarkerStyle(20);
 hPool->SetMarkerSize(1);
 hPool->Draw("EP");
 hPool->GetXaxis()->SetTitle(hDATA->GetXaxis()->GetTitle()); 
 hPool->GetYaxis()->SetTitle("2(DATA-MC)/(DATA+MC)"); 
 hPool->Draw("EP");
 gPad->SetGrid();
}

///==== Draw Stack ====
///---- if error == 1 than consider systematic error
void DrawStack(THStack* hs, int error, double syst){ 
 if (error == 1) {
   DrawStackError(hs, syst);
 }
 else {
  TObjArray* histos = hs->GetStack () ;
  Int_t number = histos->GetEntries();
  TH1F* last = (TH1F*) histos->At (number-1) ;
  last->Draw ("hist") ;
  for (int i = number-2 ; i >=0 ; --i) 
  {
   TH1F * histo = (TH1F*) histos->At (i) ;
   Style_t origStyle = histo->GetFillStyle ();
   Color_t origColor = histo->GetFillColor ();
   TH1F* dummy = (TH1F*) histo->Clone () ;
   dummy->SetFillStyle (1001) ; 
   dummy->SetFillColor (10) ;        
   dummy->Draw ("same hist") ;
   histo->Draw ("same hist") ;
  }
 }
}


///==== Draw Stack ====
void DrawStackError(THStack* hs, double syst){ 
 TObjArray* histos = hs->GetStack () ;
 Int_t number = histos->GetEntries();
 TH1F* last = (TH1F*) histos->At (number-1) ;
 last->DrawClone ("hist") ;
 for (int i = number-2 ; i >=0 ; --i) 
 {
  TH1F * histo = (TH1F*) histos->At (i) ;
  Style_t origStyle = histo->GetFillStyle ();
  Color_t origColor = histo->GetFillColor ();
  TH1F* dummy = (TH1F*) histo->Clone () ;
  dummy->SetFillStyle (1001) ; 
  dummy->SetFillColor (10) ;        
  dummy->Draw ("same hist") ;
  histo->Draw ("same hist") ;
 }
 Style_t origStyleLast = last->GetFillStyle ();
 Color_t origColorLast = last->GetFillColor ();
 last->SetFillStyle(3335);
 last->SetFillColor(kBlack);
 last->SetMarkerSize(0);
 
 std::vector <double> vErr ;
 for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
  double additionalError = last->GetBinContent(iBin+1) * syst;
  vErr.push_back(last->GetBinError(iBin+1));
  last->SetBinError(iBin+1,sqrt(additionalError*additionalError + last->GetBinError(iBin+1) * last->GetBinError(iBin+1)) );
 }
 last->DrawClone ("sameE2") ;
 //---- restore hist ----
 last->SetFillStyle(origStyleLast);
 last->SetFillColor(origColorLast);
 for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
  last->SetBinError(iBin+1, vErr.at(iBin));
 }
}

///==== Add systrematic error ====
void AddError(THStack* hs, double syst){ 
 TObjArray* histos = hs->GetStack () ;
 Int_t number = histos->GetEntries();
 TH1F* last = (TH1F*) histos->At (number-1) ;
  for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
   double additionalError = last->GetBinContent(iBin+1) * syst;
   last->SetBinError(iBin+1,sqrt(additionalError*additionalError + last->GetBinError(iBin+1) * last->GetBinError(iBin+1)) );
 }
}




///==== Draw Stack Normalized ====
///---- if error == 1 than consider systematic error
void DrawStackNormalized(THStack* hs){ 
 TObjArray* histos = hs->GetStack () ;
 Int_t number = histos->GetEntries();
 TH1F* last = (TH1F*) histos->At (number-1) ;
 last->DrawNormalized ("L") ;
 for (int i = number-2 ; i >=0 ; --i){
  TH1F * histo = (TH1F*) histos->At (i) ;
  histo->DrawNormalized ("same L") ;
 } 
}





///==== read list of systematics file ====
void ReadFileSystematics(std::string CutSystematicFile, std::vector< std::pair< int, std::pair<std::string, double> > >& listSystematics){
 std::ifstream inFile(CutSystematicFile.c_str());
 std::string buffer; 
 std::string tempVar;
 double tempSyst;
 int tempKind;
 while(!inFile.eof()){
  getline(inFile,buffer);
  if (buffer != ""){ ///---> save from empty line at the end!
   if (buffer.at(0) != '#'){     
    std::stringstream line( buffer );      
    line >> tempVar; 
    line >> tempSyst;
    line >> tempKind;
    std::cout << " Syst = " << tempVar << " ~ " << tempSyst << " ~ " << tempKind << std::endl;
    std::pair<std::string, double> tempPair(tempVar,tempSyst);
    std::pair< int, std::pair<std::string, double> > tempPairPair(tempKind, tempPair);
    listSystematics.push_back (tempPairPair);
   } 
  }
 }
}

void ReadFileSystematicsWithRegion(std::string CutSystematicFile, std::vector< std::pair< std::string, std::string> >& listSystematics){
 std::ifstream inFile(CutSystematicFile.c_str());
 std::string buffer; 
 std::string tempVarOld;
 std::string tempVarNew;
 while(!inFile.eof()){
  getline(inFile,buffer);
  if (buffer != ""){ ///---> save from empty line at the end!
   if (buffer.at(0) != '#'){     
    std::stringstream line( buffer );      
    line >> tempVarOld; 
    
    getline(inFile,buffer);
    std::stringstream lineNew( buffer );      
    tempVarNew = buffer; 
    std::cout << " Syst = " << tempVarOld << " ~> " << tempVarNew << std::endl;
    std::pair< std::string, std::string> tempPair(tempVarOld, tempVarNew);
    listSystematics.push_back (tempPair);
   } 
  }
 }
}


///==== std::string replace ====
void repl(std::string& s, const std::string& from, const std::string& to){
 std::string::size_type cnt(std::string::npos);
 if(from != to && !from.empty())
 {
  std::string::size_type pos1(0);
  std::string::size_type pos2(0);
  const std::string::size_type from_len(from.size());
  const std::string::size_type to_len(to.size());
  cnt = 0;
  
  while((pos1 = s.find(from, pos2)) != std::string::npos)
  {
   s.replace(pos1, from_len, to);
   pos2 = pos1 + to_len;
   ++cnt;
  }
 }
}

///==== modify list of cuts to include systematics ====
void ModifyCut(std::vector <std::string> & vCut, const std::vector< std::pair< int, std::pair<std::string, double> > >& listSystematics){
 std::cout << std::endl;
 for (int iSyst = 0; iSyst < listSystematics.size(); iSyst++) {
  std::string var = listSystematics.at(iSyst).second.first;
  double syst = listSystematics.at(iSyst).second.second;
  int kind = listSystematics.at(iSyst).first;
  char systChar[128];
  sprintf(systChar,"%0.2f",syst); 
  std::string newvar = "(";
  if (kind == 1){ //---- dilation
   newvar += var;
   newvar += "*(1+(";
   newvar += systChar;
   newvar += ")))";
  }
  else { //---- shift
   newvar += var;
   newvar += "+(";
   newvar += syst;
   newvar += "))";
  }
  std::cout << " newvar = " << newvar << std::endl;
  std::cout << "    var = " <<    var << std::endl;
  std::cout << "   syst = " <<   syst << std::endl;
  for (int iCut = 0; iCut < vCut.size(); iCut++) {
   repl(vCut.at(iCut),var,newvar);
  }
 } 
 for (int iCut = 0; iCut < vCut.size(); iCut++) {
  std::cout << " newCut = " << vCut.at(iCut) << std::endl;
 }
 std::cout << std::endl;
}



///==== modify list of cuts to include systematics ====
void ModifyCutWithRegion(std::vector <std::string> & vCut, const std::vector< std::pair< std::string, std::string> >& listSystematics){
 std::cout << std::endl;
 for (int iSyst = 0; iSyst < listSystematics.size(); iSyst++) {
  std::string variable = listSystematics.at(iSyst).first;
  std::string newvariable = listSystematics.at(iSyst).second;
  std::cout << "    variable = " <<    variable << std::endl;
  std::cout << " newvariable = " << newvariable << std::endl;
  for (int iCut = 0; iCut < vCut.size(); iCut++) {
   repl(vCut.at(iCut),variable,newvariable);
  }
 } 
 for (int iCut = 0; iCut < vCut.size(); iCut++) {
  std::cout << " newCut = " << vCut.at(iCut) << std::endl;
 }
 std::cout << std::endl;
}






///==== CJV with JES estimation ====
#include "../JESUnc.cxx"

int getCJVJES(std::vector<ROOT::Math::XYZTVector>& jets,
	   int q1,
	   int q2,
	   const double& EtMin,
	   const std::vector<int>* blacklist,
	   int kind){
 
 int CJV = 0;
 double etaMin = jets.at(q1).Eta();
 double etaMax = jets.at(q2).Eta();
 
 if (etaMax < etaMin) std::swap(etaMin,etaMax);
 
 for(unsigned int i = 0; i < jets.size(); ++i)
 {
  if (i==q1 || i==q2) continue;
  
  double unc = JESUnc(jets.at(i).eta(),jets.at(i).Pt());
  
  if ((jets.at(i).Et() * (1.+ (kind) * unc)) < EtMin) continue;
  
  bool skipJet = false;
  if(blacklist){
   for(unsigned int kk = 0; kk < blacklist -> size(); ++kk) {
    if(blacklist -> at(kk) == static_cast<int>(i)) skipJet = true;
   }
  }
  
  if (skipJet) continue;
  
  if(jets.at(i).Eta() > etaMax || jets.at(i).Eta() < etaMin) continue;
  
  CJV++;
 } 
 
 return CJV;
}



/** Electron isolation / ID */

bool IsEleIsolatedID_VBF( treeReader& reader,const std::vector<double>& BarrelSelections, const std::vector<double>& EndcapSelections, int iEle){
 
 bool skipEle = false;
 
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && (reader.GetFloat("electrons_tkIsoR03")->at(iEle) + reader.GetFloat("electrons_emIsoR03")->at(iEle) + reader.GetFloat("electrons_hadIsoR03_depth1")->at(iEle) + reader.GetFloat("electrons_hadIsoR03_depth2")->at(iEle))/reader.Get4V("electrons")->at(iEle).pt() > BarrelSelections.at(0)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) &&  fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) > BarrelSelections.at(1)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) &&  fabs(reader.GetFloat("electrons_deltaPhiIn")->at(iEle)) > BarrelSelections.at(2)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) &&  fabs(reader.GetFloat("electrons_deltaEtaIn")->at(iEle)) > BarrelSelections.at(3)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetInt("electrons_mishits")->at(iEle) > BarrelSelections.at(4)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dist")->at(iEle) > BarrelSelections.at(5)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dcot")->at(iEle) > BarrelSelections.at(6)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && (fabs(reader.GetFloat("electrons_dist")->at(iEle)) < BarrelSelections.at(5) && fabs(reader.GetFloat("electrons_dcot")->at(iEle)) < BarrelSelections.at(6)) ) skipEle = true;     
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dz_PV")->at(iEle) > BarrelSelections.at(7)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dxy_PV")->at(iEle) > BarrelSelections.at(8)) skipEle = true;    

 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && (reader.GetFloat("electrons_tkIsoR03")->at(iEle) + reader.GetFloat("electrons_emIsoR03")->at(iEle) + reader.GetFloat("electrons_hadIsoR03_depth1")->at(iEle) + reader.GetFloat("electrons_hadIsoR03_depth2")->at(iEle))/reader.Get4V("electrons")->at(iEle).pt() > EndcapSelections.at(0)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) &&  fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) > EndcapSelections.at(1)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) &&  fabs(reader.GetFloat("electrons_deltaPhiIn")->at(iEle)) > EndcapSelections.at(2)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) &&  fabs(reader.GetFloat("electrons_deltaEtaIn")->at(iEle)) > EndcapSelections.at(3)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetInt("electrons_mishits")->at(iEle) > EndcapSelections.at(4)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dist")->at(iEle) > EndcapSelections.at(5)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dcot")->at(iEle) > EndcapSelections.at(6)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && (fabs(reader.GetFloat("electrons_dist")->at(iEle)) < EndcapSelections.at(5) && fabs(reader.GetFloat("electrons_dcot")->at(iEle)) < EndcapSelections.at(6)) ) skipEle = true;     
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dz_PV")->at(iEle) > EndcapSelections.at(7)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dxy_PV")->at(iEle) > EndcapSelections.at(8)) skipEle = true;    
 
 return (!skipEle);
 
}

/** Electron isolation / ID 
 * with PU correction 
 */

bool IsEleIsolatedIDPUCorrected_VBF( treeReader& reader,const std::vector<double>& BarrelSelections, const std::vector<double>& EndcapSelections, int iEle){
 
 bool skipEle = false;
 
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && (reader.GetFloat("electrons_tkIsoR03")->at(iEle) + std::max((float) 0., reader.GetFloat("electrons_emIsoR03")->at(iEle)-1) + reader.GetFloat("electrons_hadIsoR03_depth1")->at(iEle) + reader.GetFloat("electrons_hadIsoR03_depth2")->at(iEle) - reader.GetFloat("rho_isolation")->at(0) * 0.3 * 0.3 * reader.GetFloat("PV_normalizedChi2")->size())/reader.Get4V("electrons")->at(iEle).pt() > BarrelSelections.at(0)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) &&  fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) > BarrelSelections.at(1)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) &&  fabs(reader.GetFloat("electrons_deltaPhiIn")->at(iEle)) > BarrelSelections.at(2)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) &&  fabs(reader.GetFloat("electrons_deltaEtaIn")->at(iEle)) > BarrelSelections.at(3)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetInt("electrons_mishits")->at(iEle) > BarrelSelections.at(4)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dist")->at(iEle) > BarrelSelections.at(5)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dcot")->at(iEle) > BarrelSelections.at(6)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && (fabs(reader.GetFloat("electrons_dist")->at(iEle)) < BarrelSelections.at(5) && fabs(reader.GetFloat("electrons_dcot")->at(iEle)) < BarrelSelections.at(6)) ) skipEle = true;     
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dz_PV")->at(iEle) > BarrelSelections.at(7)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) < 1.5) && reader.GetFloat("electrons_dxy_PV")->at(iEle) > BarrelSelections.at(8)) skipEle = true;    

 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && (reader.GetFloat("electrons_tkIsoR03")->at(iEle) + reader.GetFloat("electrons_emIsoR03")->at(iEle) + reader.GetFloat("electrons_hadIsoR03_depth1")->at(iEle) + reader.GetFloat("electrons_hadIsoR03_depth2")->at(iEle) - reader.GetFloat("rho_isolation")->at(0) * 0.3 * 0.3 * reader.GetFloat("PV_normalizedChi2")->size())/reader.Get4V("electrons")->at(iEle).pt() > EndcapSelections.at(0)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) &&  fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) > EndcapSelections.at(1)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) &&  fabs(reader.GetFloat("electrons_deltaPhiIn")->at(iEle)) > EndcapSelections.at(2)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) &&  fabs(reader.GetFloat("electrons_deltaEtaIn")->at(iEle)) > EndcapSelections.at(3)) skipEle = true;
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetInt("electrons_mishits")->at(iEle) > EndcapSelections.at(4)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dist")->at(iEle) > EndcapSelections.at(5)) skipEle = true;    
//  if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dcot")->at(iEle) > EndcapSelections.at(6)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && (fabs(reader.GetFloat("electrons_dist")->at(iEle)) < EndcapSelections.at(5) && fabs(reader.GetFloat("electrons_dcot")->at(iEle)) < EndcapSelections.at(6)) ) skipEle = true;     
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dz_PV")->at(iEle) > EndcapSelections.at(7)) skipEle = true;    
 if ((fabs(reader.Get4V("electrons")->at(iEle).Eta()) > 1.5) && reader.GetFloat("electrons_dxy_PV")->at(iEle) > EndcapSelections.at(8)) skipEle = true;    

 return (!skipEle);
 
}


/** Electron definition in merged analysis */
bool IsEle_VBFMerged( treeReader& reader, int iEle){
 bool skipEle = false;
 
  if ( ! (  reader.GetInt("electrons_numberOfLostHits")-> at(iEle)  == 0  &&   reader.GetFloat("electrons_passVtxConvert")  -> at(iEle)  != 0   )   )  skipEle = true;
  
//    ELE_MERGE_CONV =  " ( gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0 && userFloat('convValueMapProd:passVtxConvert') != 0 ) "

  if ( ! (
         ( reader.GetInt("electrons_isEB")-> at(iEle)  == 1   
           && reader.Get4V("electrons")->at(iEle).pt() < 20 
           && fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) < 0.01
           && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) > -0.03   && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) < 0.03
           && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) > -0.004 && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) < 0.004
           && reader.GetFloat("electrons_hOverE")->at(iEle) < 0.025       && (reader.GetFloat("electrons_fbrem")->at(iEle) > 0.15  ||   ( fabs(reader.GetFloat("electrons_etaSC")->at(iEle)) < 1. &&    reader.GetFloat("electrons_eSCOverP")->at(iEle) > 0.95))
          )   ||
          ( reader.GetInt("electrons_isEB")-> at(iEle)  == 0   
           && reader.Get4V("electrons")->at(iEle).pt() < 20 
           && fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) < 0.03
           && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) > -0.02   && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) < 0.02
           && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) > -0.005 && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) < 0.005
           && reader.GetFloat("electrons_hOverE")->at(iEle) < 0.1 && reader.GetFloat("electrons_fbrem")->at(iEle) > 0.15 
          )  ||
          ( reader.GetInt("electrons_isEB")-> at(iEle)  == 1 
           && reader.Get4V("electrons")->at(iEle).pt() > 20 
           && fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) < 0.01
           && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) > -0.06   && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) < 0.06
           && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) > -0.004 && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) < 0.004
           && reader.GetFloat("electrons_hOverE")->at(iEle) < 0.04
          )  ||
          ( reader.GetInt("electrons_isEB")-> at(iEle)  == 0 
           && reader.Get4V("electrons")->at(iEle).pt() > 20 
           && fabs(reader.GetFloat("electrons_sigmaIetaIeta")->at(iEle)) < 0.03
           && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) > -0.03   && reader.GetFloat("electrons_deltaPhiIn")->at(iEle) < 0.03
           && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) > -0.007 && reader.GetFloat("electrons_deltaEtaIn")->at(iEle) < 0.007
           && reader.GetFloat("electrons_hOverE")->at(iEle) < 0.1
          )
       )
      )    skipEle = true;
      
//          ELE_MERGE_ID   =  ("(( isEB && pt < 20 && sigmaIetaIeta < 0.01 &&" +           
//                   " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
//                   " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
//                   " hadronicOverEm < 0.025 && (fbrem > 0.15 || ( abs(superCluster.eta) < 1. && eSuperClusterOverP > 0.95 )) ) ||" +
//                   " ( (!isEB) && pt < 20 && sigmaIetaIeta < 0.03 &&" +           
//                   " deltaPhiSuperClusterTrackAtVtx > -0.02 && deltaPhiSuperClusterTrackAtVtx < 0.02 &&" +
//                   " deltaEtaSuperClusterTrackAtVtx > -0.005 && deltaEtaSuperClusterTrackAtVtx < 0.005 &&" +
//                   " hadronicOverEm < 0.1 && fbrem > 0.15) ||" +
//                   " ( isEB && pt > 20 && sigmaIetaIeta < 0.01 &&" +           
//                   " deltaPhiSuperClusterTrackAtVtx > -0.06 && deltaPhiSuperClusterTrackAtVtx < 0.06 &&" +
//                   " deltaEtaSuperClusterTrackAtVtx > -0.004 && deltaEtaSuperClusterTrackAtVtx < 0.004 &&" +
//                   " hadronicOverEm < 0.04) ||" +
//                   " ( (!isEB) && pt > 20 && sigmaIetaIeta < 0.03  &&  " +
//                   " deltaPhiSuperClusterTrackAtVtx > -0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 &&" +
//                   " deltaEtaSuperClusterTrackAtVtx > -0.007 && deltaEtaSuperClusterTrackAtVtx < 0.007 && " + 
//                   " hadronicOverEm < 0.1 ) ) ")
//
          
  if ( !    (  fabs(reader.GetFloat("electrons_tip")-> at(iEle)) < 0.02    &&    fabs(reader.GetFloat("electrons_dzPV")-> at(iEle)) < 0.1   )  )  skipEle = true;

//     ELE_MERGE_IP   =   "( abs(userFloat('tip')) < 0.02 && abs(userFloat('dzPV')) < 0.1 )"
        
        
  if ( !    (  
                 ( reader.GetInt("electrons_isEB")-> at(iEle)  == 1    &&   reader.GetFloat("electrons_PFIso")-> at(iEle) /  reader.Get4V("electrons")->at(iEle).pt()  < 0.13  )    ||
                 ( reader.GetInt("electrons_isEB")-> at(iEle)  == 0    &&   reader.GetFloat("electrons_PFIso")-> at(iEle) /  reader.Get4V("electrons")->at(iEle).pt()  < 0.09  )
              )
      )  skipEle = true;
        
//    ELE_MERGE_ISO  =   "( (isEB && " + SMURF_ISO + " < 0.13) || ( !isEB && " + SMURF_ISO + " < 0.09 ) )"
//    SMURF_ISO = ("( ( userFloat('eleSmurfPF') ) / pt )")
        
        
/*        
  double lik =  reader.GetFloat("egammaIDLikelihood")-> at(iEle) ;
  bool ELE_ID_LH_80_2011 = (  
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 1  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)  == 0 && lik > 0.193 ) ||
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 1  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)   >  0 && lik > 1.345 ) ||
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 0  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)  == 0 && lik > 0.810 ) ||
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 0  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)   >  0 && lik > 3.021 )
         );
         
  bool ELE_ID_LH_90_2011 = (  
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 1  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)  == 0 && lik > -1.497 ) ||
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 1  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)   >  0 && lik > -1.521 ) ||
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 0  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)  == 0 && lik > -2.571 ) ||
         (reader.GetInt("electrons_isEB")-> at(iEle)  == 0  &&  reader.GetInt("electrons_numberOfBrems")-> at(iEle)   >  0 && lik > -0.657 )
         );
        
  if ( !    (  
                 ( reader.Get4V("electrons")->at(iEle).pt()   < 20 &&   ELE_ID_LH_80_2011  )    ||
                 ( reader.Get4V("electrons")->at(iEle).pt() >= 20 &&   ELE_ID_LH_90_2011  ) 
              )
      )  skipEle = true;

*/
           
///ELE_MERGE_ID2   =  ("( (pt < 20 && " + ELE_ID_LH_80_2011 +") || (pt >= 20 && "+ ELE_ID_LH_90_2011 + ") )")
//           
//           ELE_ID_LH_80_2011=("( (  isEB  && numberOfBrems == 0 && electronID('egammaIDLikelihood') >  0.193 ) ||" +           
//                   "  (  isEB  && numberOfBrems  > 0 && electronID('egammaIDLikelihood') >  1.345 ) ||" + 
//                   "  ((!isEB) && numberOfBrems == 0 && electronID('egammaIDLikelihood') >  0.810 ) ||" + 
//                   "  ((!isEB) && numberOfBrems  > 0 && electronID('egammaIDLikelihood') >  3.021 ) )" )
//
//           ELE_ID_LH_90_2011=("( (  isEB  && numberOfBrems == 0 && electronID('egammaIDLikelihood') > -1.497 ) ||" +           
//                   "  (  isEB  && numberOfBrems  > 0 && electronID('egammaIDLikelihood') > -1.521 ) ||" + 
//                   "  ((!isEB) && numberOfBrems == 0 && electronID('egammaIDLikelihood') > -2.571 ) ||" + 
//                   "  ((!isEB) && numberOfBrems  > 0 && electronID('egammaIDLikelihood') > -0.657 ) )" )
//
           
 return (!skipEle);
}




/** Muon definition in merged analysis */

bool IsMu_VBFMerged( treeReader& reader, int iMu){
 bool skipMu = false;
 
  if (!(
       (   (
             reader.GetInt("muons_global")->at(iMu)   && 
             reader.GetFloat("muons_normalizedChi2")->at(iMu)<10  &&   
             reader.GetInt("muons_numberOfValidMuonHits")->at(iMu)>0 && 
             reader.GetInt("muons_numberOfMatches")->at(iMu)>1
             ) ||
           (
           reader.GetInt("muons_tracker")->at(iMu) && reader.GetFloat("muons_TMLastStationTight")->at(iMu) > 0.5  
           ) 
        )
      && reader.GetInt("muons_innerTrack_found")->at(iMu)>10 
      && reader.GetInt("muons_numberOfValidPixelHits")->at(iMu)>0       //------------------------> to be added!
      && fabs(reader.GetFloat("muons_trackPtErrorOverPt")->at(iMu)) < 0.10
      )   
     ) skipMu=true;

//    MUON_ID_CUT=("(( (isGlobalMuon() && "
//                 "    globalTrack.normalizedChi2 <10 &&" +
//                 "    globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + 
//                 "    numberOfMatches > 1 ) || " + 
//                 "   (isTrackerMuon() && muonID('TMLastStationTight')) ) && " + 
//                 " innerTrack.found >10 &&" +
//                 " innerTrack.hitPattern().numberOfValidPixelHits > 0 && " + 
//                 " abs(track.ptError / pt) < 0.10 )")
//    


 float SMURF_ISO = reader.GetFloat("muons_PFIso")->at(iMu) / (reader.Get4V("muons")->at(iMu).pt());
 
 if(!(
    (fabs(reader.Get4V("muons")->at(iMu).eta())<1.479     && reader.Get4V("muons")->at(iMu).pt()>20 && SMURF_ISO<0.13) ||
    (fabs(reader.Get4V("muons")->at(iMu).eta())>=1.479   &&  reader.Get4V("muons")->at(iMu).pt()>20 && SMURF_ISO<0.09) ||
    (fabs(reader.Get4V("muons")->at(iMu).eta())<1.479     &&  reader.Get4V("muons")->at(iMu).pt()<=20 && SMURF_ISO<0.06) ||
    (fabs(reader.Get4V("muons")->at(iMu).eta())>=1.479   &&  reader.Get4V("muons")->at(iMu).pt()<=20 && SMURF_ISO<0.05)
    )
   )  skipMu=true;

//    SMURF_ISO = ("( userFloat('muSmurfPF') )/ pt")
//    MUON_MERGE_ISO  =   ("( (abs(eta) < 1.479 && pt >  20 && " + SMURF_ISO + " < 0.13) || ( abs(eta) >= 1.479 && pt >  20 && " + SMURF_ISO + " < 0.09 ) || " + 
//                         "  (abs(eta) < 1.479 && pt <= 20 && " + SMURF_ISO + " < 0.06) || ( abs(eta) >= 1.479 && pt <= 20 && " + SMURF_ISO + " < 0.05 ) )  ")
//    

   if(!(
       (
        (reader.Get4V("muons")->at(iMu).pt()>=20 && fabs(reader.GetFloat("muons_tip")->at(iMu))<0.02)||
        (reader.Get4V("muons")->at(iMu).pt()<20 && fabs(reader.GetFloat("muons_tip")->at(iMu))<0.01)
       ) &&
      fabs(reader.GetFloat("muons_dzPV")->at(iMu))<0.1)
   ) skipMu=true;

//    MUON_MERGE_IP  = ("( ( (pt >= 20 && abs(userFloat('tip')) < 0.02) || (pt < 20 && abs(userFloat('tip')) < 0.01) ) && " +
//                      "  abs(userFloat('dzPV'))  < 0.1 )" )
//    

 return (!skipMu); 
}


//  ------------------------------------------------------------

/** Muon isolation / ID */

bool IsMuIsolatedID_VBF( treeReader& reader,const std::vector<double>& Selections, int iMu){
 
 bool skipMu = false;
 
 if ( (reader.GetFloat("muons_tkIsoR03")->at(iMu) + reader.GetFloat("muons_emIsoR03")->at(iMu) + reader.GetFloat("muons_hadIsoR03")->at(iMu)) / reader.Get4V("muons")->at(iMu).pt() > Selections.at(0) ) skipMu = true;
 if (reader.GetFloat("muons_normalizedChi2")->at(iMu) >= Selections.at(1) )           skipMu = true;
 if (reader.GetInt("muons_numberOfValidTrackerHits")->at(iMu) <=  Selections.at(2)  )     skipMu = true;
 if (reader.GetInt("muons_numberOfValidMuonHits")->at(iMu) <=  Selections.at(3)  )     skipMu = true;
 if (reader.GetInt("muons_tracker")->at(iMu) != Selections.at(4)  )     skipMu = true;
 if (reader.GetInt("muons_standalone")->at(iMu) != Selections.at(5)  )     skipMu = true;
 if (reader.GetInt("muons_global")->at(iMu) != Selections.at(6)  )     skipMu = true;
 if (reader.GetFloat("muons_dz_PV")->at(iMu)  > Selections.at(7)) skipMu = true;    
 if (reader.GetFloat("muons_dxy_PV")->at(iMu) > Selections.at(8)) skipMu = true;    
 
 return (!skipMu);
}



bool IsMuIsolatedIDPUCorrected_VBF( treeReader& reader,const std::vector<double>& Selections, int iMu){
 
 bool skipMu = false;
 
 if ( (reader.GetFloat("muons_tkIsoR03")->at(iMu) + reader.GetFloat("muons_emIsoR03")->at(iMu) + reader.GetFloat("muons_hadIsoR03")->at(iMu)  - reader.GetFloat("rho_isolation")->at(0) * 0.3 * 0.3 * reader.GetFloat("PV_normalizedChi2")->size()) / reader.Get4V("muons")->at(iMu).pt() > Selections.at(0) ) skipMu = true;
 if (reader.GetFloat("muons_normalizedChi2")->at(iMu) >= Selections.at(1) )           skipMu = true;
 if (reader.GetInt("muons_numberOfValidTrackerHits")->at(iMu) <=  Selections.at(2)  )     skipMu = true;
 if (reader.GetInt("muons_numberOfValidMuonHits")->at(iMu) <=  Selections.at(3)  )     skipMu = true;
 if (reader.GetInt("muons_tracker")->at(iMu) != Selections.at(4)  )     skipMu = true;
 if (reader.GetInt("muons_standalone")->at(iMu) != Selections.at(5)  )     skipMu = true;
 if (reader.GetInt("muons_global")->at(iMu) != Selections.at(6)  )     skipMu = true;
 if (reader.GetFloat("muons_dz_PV")->at(iMu)  > Selections.at(7)) skipMu = true;    
 if (reader.GetFloat("muons_dxy_PV")->at(iMu) > Selections.at(8)) skipMu = true;    
 
 return (!skipMu);
}



//  ------------------------------------------------------------

/** Soft Muon identification */

bool IsMu_Soft( treeReader& reader, int iMu){
 //==== soft mu veto
 //==== http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Mangano/WWAnalysis/Filters/src/SoftMuonVeto.cc
 //====
  bool softMu=false;
  float SMURF_ISO = reader.GetFloat("muons_PFIso")->at(iMu) / (reader.Get4V("muons")->at(iMu).pt());
  if (
     reader.Get4V("muons")->at(iMu).pt()>3 && 
     reader.GetInt("muons_tracker")->at(iMu) &&

     reader.GetFloat("muons_TMLastStationAngTight")->at(iMu) > 0.5 &&  //=========> angTIght to be added!
     
     reader.GetInt("muons_innerTrack_found")->at(iMu)>10 &&
     fabs(reader.GetFloat("muons_dz_PV")->at(iMu))<0.1&&
     fabs(reader.GetFloat("muons_tip")->at(iMu))<0.2 
     ){
     //if(reader.Get4V("muons")->at(iMu).pt()>20)
     //{ if(SMURF_ISO>0.1)
     //     softMu=true;
    //}
   //else {
    softMu=true;
   //}
  } 
 return(softMu);
}




//  ------------------------------------------------------------
/** Jet ID */


bool IsJetID( treeReader& reader, int iJet){
 //==== jet ID
 //==== http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Mangano/WWAnalysis/DataFormats/src/SkimEvent.cc
 //====

  bool skipJet =false;
  unsigned int multiplicity = reader.GetInt("jets_neutralMultiplicity")->at(iJet)+reader.GetInt("jets_chargedMultiplicity")->at(iJet);
  if(
     reader.GetFloat("jets_neutralEmEnergyFraction")->at(iJet)>=0.99 ||
     reader.GetFloat("jets_neutralHadronEnergyFraction")->at(iJet)>=0.99 ||
     multiplicity ==0
     )   skipJet=true;

  if(fabs(reader.Get4V("jets")->at(iJet).Eta()) <2.4) { 
     if( 
         reader.GetFloat("jets_chargedEmEnergyFraction")->at(iJet)>=0.99 ||
         reader.GetFloat("jets_chargedHadronEnergyFraction")->at(iJet) ==0 ||
         reader.GetInt("jets_chargedMultiplicity")->at(iJet)==0
        ) skipJet=true;
   }	
 
  return (!skipJet);
}














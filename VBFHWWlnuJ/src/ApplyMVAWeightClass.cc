#include "ApplyMVAWeightClass.h"

// constructor from files
ApplyMVAWeightClass::ApplyMVAWeightClass(const std::vector<TFile*> & SampleFileList, const std::string & TreeName,
                                         const std::string & inputFilePath , const std::string & Label ){

  (*this).SetInputTree (SampleFileList, TreeName) ;

  (*this).SetReaderTree();

  (*this).SetLabel (Label);

  (*this).SetInputFilePath (inputFilePath) ;
  
  reader_ = new TMVA::Reader( TreeName_+"_"+Label_);


}

// constructor from trees
ApplyMVAWeightClass::ApplyMVAWeightClass(const std::vector<TTree*> & SampleTreeList, const std::string & TreeName,
		                         const std::string & inputFilePath , const std::string & Label ){


  (*this).SetTreeName (TreeName) ;

  (*this).SetInputTree (SampleTreeList) ;

  (*this).SetReaderTree (SampleTreeList) ;

  (*this).SetLabel (Label);

  (*this).SetInputFilePath (inputFilePath) ;

  reader_ = new TMVA::Reader( TreeName_+"_"+Label_);

}

// Deconstructor                                                                                                                                                         
ApplyMVAWeightClass::~ApplyMVAWeightClass(){

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { 
   if(SampleTreeList_.at(iTree)!=0)  delete SampleTreeList_.at(iTree) ; 
  }

  if(reader_!=0) reader_->Delete() ;

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { if(treeReader_.at(iTree)!=0) delete treeReader_.at(iTree) ;}

}

// Set Input variables and spectator 
void ApplyMVAWeightClass::AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables){

  (*this).SetTrainingVariables(mapTrainingVariables);
  (*this).SetSpectatorVariables(mapSpectatorVariables);

  setTrainingVariables_  = new std::vector<Float_t> ( mapTrainingVariables_.size()) ;
  setSpectatorVariables_ = new std::vector<Float_t> ( mapSpectatorVariables_.size()) ;

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ){
    reader_->AddVariable(mapTrainingVariables_.at(iVar),&setTrainingVariables_->at(iVar));
  }

  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ){
    reader_->AddSpectator(mapSpectatorVariables_.at(iVar),&setSpectatorVariables_->at(iVar));
  }

  return;
}

// Prepare the tree -> set the preselection cut
void ApplyMVAWeightClass::AddPrepareReader (const std::string & LeptonType, const std::string & preselectionCutType,
					    std::vector<double> * JetPtBinOfTraining, const int & pTBin ){

  if(JetPtBinOfTraining!=NULL || !JetPtBinOfTraining ){
    pTJetMin_ = JetPtBinOfTraining->at(pTBin) ;
    pTJetMax_ = JetPtBinOfTraining->at(pTBin+1) ;
  }
  else{
    pTJetMin_ = 0;
    pTJetMax_ = 2000 ;
  }

  preselectionCutType_ = preselectionCutType ;
  LeptonType_ = LeptonType ;

  return ;
}


// Book the MVA method
void ApplyMVAWeightClass::BookMVAWeight (const std::string & methodName, const std::string & weightFile, const std::string & nameBranch) {

  methodName_=methodName ;
  weightFile_=weightFile ; 
  reader_->BookMVA(methodName_.c_str(),(inputFilePath_+"/"+weightFile_).c_str()) ;
  nameBranch_ = nameBranch ;  

  return ;

}

// Fill the MVA weight
void ApplyMVAWeightClass::FillMVAWeight (const std::string & LeptonType, const std::string & preselectionCutType) {

  LeptonType_ = LeptonType ;
  preselectionCutType_ = preselectionCutType ;

  std::cout<<std::endl;

  // loop on the sample list and create the branch
  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++){

      newBranch_  = SampleTreeList_.at(iTree)->Branch(nameBranch_.c_str(),&weight_,(nameBranch_+"/F").c_str());

      for( int iEntry = 0; iEntry < SampleTreeList_.at(iTree)->GetEntries(); iEntry++){
	// loop on the entries   
       if(iEntry%10000 == 0 ) std::cout<<" Entry = "<<iEntry<<std::endl;
      
       SampleTreeList_.at(iTree)->GetEntry(iEntry);  weight_ = 0 ;
       // fill the value of the training and spectator variable for each event
       for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ) 
	setTrainingVariables_->at(iVar) = treeReader_.at(iTree)->getFloat(mapTrainingVariables_.at(iVar).c_str())[0] ; 

      for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ) 
	setSpectatorVariables_->at(iVar) = treeReader_.at(iTree)->getFloat(mapSpectatorVariables_.at(iVar).c_str())[0] ;

      // preselection type definition --> only the basic cut SB+SR since we want to skip only the events that are never used in the analysis
      bool isGoodEvent = false ;
      if (preselectionCutType_ == "basicPreselectionCutEXO" && (LeptonType_ == "Mu" || LeptonType_ == "mu" || LeptonType_ == "Muon" || LeptonType_ == "muon") ) {         

	isGoodEvent = treeReader_.at(iTree)->getInt("issignal")[0] && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 40 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 50 && ( treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && 
                      treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      else if (preselectionCutType_ == "basicPreselectionCutEXO" && (LeptonType_ == "El" || LeptonType_ == "el" || LeptonType_ == "Electron" || LeptonType_ == "electron") ){ 
      
	isGoodEvent = treeReader_.at(iTree)->getFloat("issignal")[0]  && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 80 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 90 && (treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && 
                      treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      else if (preselectionCutType_ == "basicPreselectionCutEXO" && (LeptonType_ == "MuEl" || LeptonType_ == "muel" || LeptonType_ == "MuonEle" || LeptonType_ == "muonele") ){ 
      
	isGoodEvent = treeReader_.at(iTree)->getFloat("issignal")[0]  && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 80 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 90 && (treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && 
                      treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      
      else if (preselectionCutType_ == "basicPreselectionCutHiggs" && (LeptonType_ == "Mu" || LeptonType_ == "mu" || LeptonType_ == "Muon" || LeptonType_ == "muon") ) {         
	isGoodEvent = treeReader_.at(iTree)->getInt("issignal")[0] && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 50 && 
 	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 30 && treeReader_.at(iTree)->getFloat("numberJetBin")[0]<2 && 
                      ( treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      else if (preselectionCutType_ == "basicPreselectionCutHiggs" && (LeptonType_ == "El" || LeptonType_ == "el" || LeptonType_ == "Electron" || LeptonType_ == "electron") ){       
	isGoodEvent = treeReader_.at(iTree)->getFloat("issignal")[0]  && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 70 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 35 && treeReader_.at(iTree)->getFloat("numberJetBin")[0]<2 &&
                      (treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      else if (preselectionCutType_ == "basicPreselectionCutHiggs" && (LeptonType_ == "MuEl" || LeptonType_ == "muel" || LeptonType_ == "MuonEle" || LeptonType_ == "muonele") ){       
	isGoodEvent = treeReader_.at(iTree)->getFloat("issignal")[0]  && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 70 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 35 && treeReader_.at(iTree)->getFloat("numberJetBin")[0]<2 &&
                      (treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }

      else if (preselectionCutType_ == "basicVBFPreselectionCutHiggs" && (LeptonType_ == "Mu" || LeptonType_ == "mu" || LeptonType_ == "Muon" || LeptonType_ == "muon") ) {         
	isGoodEvent = treeReader_.at(iTree)->getInt("issignal")[0] && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 50 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 30 && treeReader_.at(iTree)->getFloat("numberJetBin")[0]>=2 && 
                      ( treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      else if (preselectionCutType_ == "basicVBFPreselectionCutHiggs" && (LeptonType_ == "El" || LeptonType_ == "el" || LeptonType_ == "Electron" || LeptonType_ == "electron") ){ 
      
	isGoodEvent = treeReader_.at(iTree)->getFloat("issignal")[0]  && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 70 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 35 && treeReader_.at(iTree)->getFloat("numberJetBin")[0]>=2 &&
                      (treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      else if (preselectionCutType_ == "basicVBFPreselectionCutHiggs" && (LeptonType_ == "MuEl" || LeptonType_ == "muel" || LeptonType_ == "MuonEle" || LeptonType_ == "muonele") ){       
	isGoodEvent = treeReader_.at(iTree)->getFloat("issignal")[0]  && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 70 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 35 && treeReader_.at(iTree)->getFloat("numberJetBin")[0]>=2 && 
                      (treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }

      // if is a good event -> fill the branch with the value of the output      
      if(isGoodEvent) { weight_ = reader_->EvaluateMVA(methodName_.c_str()); 
                        newBranch_->Fill() ;
      }
      else { weight_ = -100. ; 
            newBranch_->Fill() ; 
      } 
    }       
    SampleTreeList_.at(iTree)->Write("", TObject::kOverwrite);
  }

  return ;

}


// Set input tree
void ApplyMVAWeightClass::SetInputTree (const std::vector<TFile*> & SampleFileList,  const std::string & TreeName){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;

  for(size_t iFile = 0 ; iFile < SampleFileList.size() ; iFile ++){

    if(SampleFileList.at(iFile)!=0)  SampleTreeList_.push_back((TTree*) SampleFileList.at(iFile)->Get(TreeName_.c_str()));

  }

  return ;


}

void ApplyMVAWeightClass::SetInputTree (const std::vector<TTree*> & SampleTreeList){

  if(SampleTreeList.size()!=0) SampleTreeList_ = SampleTreeList ;

  return ;

}

// Set Training variables name
void ApplyMVAWeightClass::SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables){

  if(mapTrainingVariables.size()!=0) mapTrainingVariables_=mapTrainingVariables;

  return ;

}

void ApplyMVAWeightClass::SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables){

  if(mapSpectatorVariables.size()!=0) mapSpectatorVariables_=mapSpectatorVariables;
 
  return ;

}

// Set input path
void ApplyMVAWeightClass::SetInputFilePath ( const std::string & InputFilePath){

  if(!InputFilePath.empty()) inputFilePath_ = InputFilePath ;

  return;

}

void ApplyMVAWeightClass::SetTreeName ( const std::string & TreeName ){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;

  return ;

}

void ApplyMVAWeightClass::SetLabel ( const std::string & Label ){

  Label_ = Label ;

  return ;

}

void ApplyMVAWeightClass::SetReaderTree(){

  for(size_t iTree = 0; iTree<SampleTreeList_.size() ; iTree++)
    treeReader_.push_back(new treeReader((TTree*) SampleTreeList_.at(iTree), false) );

  return ;

}



void ApplyMVAWeightClass::SetReaderTree(const std::vector<TTree*> & SampleTreeList){

  for(size_t iTree = 0; iTree<SampleTreeList.size() ; iTree++)
    treeReader_.push_back(new treeReader((TTree*) SampleTreeList.at(iTree), false));

  return ;

}

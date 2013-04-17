#include "ApplyMVAWeightClass.h"

ApplyMVAWeightClass::ApplyMVAWeightClass(const std::vector<TFile*> & SampleFileList, const std::string & TreeName,
                                         const std::string & inputFilePath , const std::string & Label ){

  SetInputTree (SampleFileList, TreeName) ;

  SetReaderTree () ;

  SetLabel ( Label );

  SetInputFilePath ( inputFilePath ) ;
  
  reader_ = new TMVA::Reader( TreeName_+"_"+Label_);


}

ApplyMVAWeightClass::ApplyMVAWeightClass(const std::vector<TTree*> & SampleTreeList, const std::string & TreeName,
		                         const std::string & inputFilePath , const std::string & Label ){


  SetTreeName (TreeName) ;

  SetInputTree (SampleTreeList) ;

  SetReaderTree (SampleTreeList) ;

  SetLabel ( Label );

  SetInputFilePath ( inputFilePath ) ;

  reader_ = new TMVA::Reader( TreeName_+"_"+Label_);

}


// Deconstructor                                                                                                                                                                                 
ApplyMVAWeightClass::~ApplyMVAWeightClass(){

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { if(SampleTreeList_.at(iTree)!=0)  delete SampleTreeList_.at(iTree) ; }

  if(reader_!=0) reader_->Delete() ;

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++) { if(treeReader_.at(iTree)!=0) delete treeReader_.at(iTree) ;}

}

// Set Input variables and  spectator 

void ApplyMVAWeightClass::AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables){

  SetTrainingVariables(mapTrainingVariables);
  SetSpectatorVariables(mapSpectatorVariables);

  setTrainingVariables_  = new std::vector<Float_t> ( mapTrainingVariables_.size()) ;
  setSpectatorVariables_ = new std::vector<Float_t> ( mapSpectatorVariables_.size()) ;

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ){
    reader_->AddVariable(mapTrainingVariables_.at(iVar),&setTrainingVariables_->at(iVar));
  }

  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ){
    reader_->AddSpectator(mapSpectatorVariables_.at(iVar),&setSpectatorVariables_->at(iVar));
  }

}

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

}


void ApplyMVAWeightClass::BookMVAWeight (const std::string & methodName, const std::string & weightFile, const std::string & nameBranch) {

  methodName_=methodName ;
  weightFile_=weightFile ; 
  reader_->BookMVA(methodName_.c_str(),(inputFilePath_+"/"+weightFile_).c_str()) ;
  nameBranch_ = nameBranch ;  
}

void ApplyMVAWeightClass::FillMVAWeight (const std::string & LeptonType, const std::string & preselectionCutType) {

  LeptonType_ = LeptonType ;
  preselectionCutType_ = preselectionCutType ;
  std::cout<<std::endl;

  for(size_t iTree = 0; iTree < SampleTreeList_.size() ; iTree++){

      newBranch_  = SampleTreeList_.at(iTree)->Branch(nameBranch_.c_str(),&weight_,(nameBranch_+"/F").c_str());
      for( int iEntry = 0; iEntry < SampleTreeList_.at(iTree)->GetEntries(); iEntry++){


      if(iEntry%10000 == 0 ) std::cout<<" Entry = "<<iEntry<<std::endl;
      
      SampleTreeList_.at(iTree)->GetEntry(iEntry);  weight_ = 0 ;

      for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ) 
	setTrainingVariables_->at(iVar) = treeReader_.at(iTree)->getFloat(mapTrainingVariables_.at(iVar).c_str())[0] ; 

      for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ) 
	setSpectatorVariables_->at(iVar) = treeReader_.at(iTree)->getFloat(mapSpectatorVariables_.at(iVar).c_str())[0] ;

      bool isGoodEvent = false ;

      if (preselectionCutType_ == "basicPreselectionCut" && (LeptonType_ == "Mu" || LeptonType_ == "mu" || LeptonType_ == "Muon" || LeptonType_ == "muon") ) {         

	isGoodEvent = treeReader_.at(iTree)->getInt("issignal")[0] && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 40 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 50 && ( treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && 
                      treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      }
      else if (preselectionCutType_ == "basicPreselectionCut" && (LeptonType_ == "El" || LeptonType_ == "el" || LeptonType_ == "Electron" || LeptonType_ == "electron") ) 
      
	isGoodEvent = treeReader_.at(iTree)->getFloat("issignal")[0]  && treeReader_.at(iTree)->getFloat("v_pt")[0] > 200 && treeReader_.at(iTree)->getFloat("pfMET")[0] > 80 && 
	              treeReader_.at(iTree)->getFloat("l_pt")[0] > 90 && (treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] > pTJetMin_ && 
                      treeReader_.at(iTree)->getFloat("ungroomed_jet_pt")[0] < pTJetMax_ )  ;
      
      if(isGoodEvent) { weight_ = reader_->EvaluateMVA(methodName_.c_str()); newBranch_->Fill() ;}
      else { weight_ = -100. ; newBranch_->Fill() ; } 
    } 
      
    SampleTreeList_.at(iTree)->Write("", TObject::kOverwrite);
  }
}



void ApplyMVAWeightClass::SetInputTree (const std::vector<TFile*> & SampleFileList,  const std::string & TreeName){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;

  for(size_t iFile = 0 ; iFile < SampleFileList.size() ; iFile ++){

    if(SampleFileList.at(iFile)!=0)  SampleTreeList_.push_back((TTree*) SampleFileList.at(iFile)->Get(TreeName_.c_str()));

  }

}

void ApplyMVAWeightClass::SetInputTree (const std::vector<TTree*> & SampleTreeList){

  if(SampleTreeList.size()!=0) SampleTreeList_ = SampleTreeList ;

}


void ApplyMVAWeightClass::SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables){

  if(mapTrainingVariables.size()!=0) mapTrainingVariables_=mapTrainingVariables;

}

void ApplyMVAWeightClass::SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables){

  if(mapSpectatorVariables.size()!=0) mapSpectatorVariables_=mapSpectatorVariables;

}

void ApplyMVAWeightClass::SetInputFilePath ( const std::string & InputFilePath){

  if(!InputFilePath.empty()) inputFilePath_ = InputFilePath ;

}

void ApplyMVAWeightClass::SetTreeName ( const std::string & TreeName ){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;


}

void ApplyMVAWeightClass::SetLabel ( const std::string & Label ){

  Label_ = Label ;

}

void ApplyMVAWeightClass::SetReaderTree(){

  for(size_t iTree = 0; iTree<SampleTreeList_.size() ; iTree++)
    treeReader_.push_back(new treeReader((TTree*) SampleTreeList_.at(iTree), false) );
}



void ApplyMVAWeightClass::SetReaderTree(const std::vector<TTree*> & SampleTreeList){

  for(size_t iTree = 0; iTree<SampleTreeList.size() ; iTree++)
    treeReader_.push_back(new treeReader((TTree*) SampleTreeList.at(iTree), false));

}

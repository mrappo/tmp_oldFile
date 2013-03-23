#include "TrainingMVAClass.h"

// Constructor

TrainingMVAClass::TrainingMVAClass(const std::vector<TFile*> & signalFileList, const std::vector<TFile*> & backgroundFileList, const std::string & TreeName,
				   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label ){

   SetSignalTree (signalFileList,TreeName) ;

   SetBackgroundTree (backgroundFileList,TreeName) ;

   SetLabel ( Label );

   SetOutputFile ( outputFilePath , outputFileName ) ;

   factory_ = new TMVA::Factory (TreeName_+"_"+Label_,outputFile_, Form("!V:!Silent:%sColor:DrawProgressBar:AnalysisType=Classification", gROOT->IsBatch()?"!":""));

}

TrainingMVAClass::TrainingMVAClass(const std::vector<TTree*> & signalTreeList, const std::vector<TTree*> & backgroundTreeList,  const std::string & TreeName,
                                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label ){
   

   SetTreeName (TreeName) ;

   SetSignalTree (signalTreeList) ;

   SetBackgroundTree (backgroundTreeList) ;

   SetLabel ( Label );

   SetOutputFile ( outputFilePath , outputFileName ) ;

   factory_ = new TMVA::Factory (TreeName_+"_"+Label_,outputFile_, Form("!V:!Silent:%sColor:DrawProgressBar:AnalysisType=Classification", gROOT->IsBatch()?"!":""));


}

// Deconstructor

TrainingMVAClass::~TrainingMVAClass(){

  for(size_t iTree = 0; iTree < signalTreeList_.size() ; iTree++) { if(signalTreeList_.at(iTree)!=0)  delete signalTreeList_.at(iTree) ; }

  for(size_t iTree = 0; iTree < backgroundTreeList_.size() ; iTree++) {if(backgroundTreeList_.at(iTree)!=0)  delete backgroundTreeList_.at(iTree) ;}

  if(outputFile_!=0) outputFile_->Delete() ;

  if(preselectionCut_!=0) preselectionCut_->Delete() ;
 
  if(factory_!=0) factory_->Delete() ;

  if(treeReader_!=0) delete treeReader_;

}

// AddTrainingVariables in the MVA

void TrainingMVAClass::AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables){

  SetTrainingVariables(mapTrainingVariables);
  SetSpectatorVariables(mapSpectatorVariables);

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ ){
    std::cout<<" train " <<mapTrainingVariables_.at(iVar)<<std::endl;
    factory_->AddVariable(mapTrainingVariables_.at(iVar)+" := "+mapTrainingVariables_.at(iVar),'F');
  }

  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ ){
    std::cout<<" train " <<mapSpectatorVariables_.at(iVar)<<std::endl;
    factory_->AddSpectator(mapSpectatorVariables_.at(iVar),'F');
  }    
}


// Book MVA Training Variables 

void TrainingMVAClass::BookMVATrees (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight){

  SetGlobalSampleWeight(signalGlobalWeight,backgroundGlobalWeight);
 
  if(signalGlobalWeight.size() == signalTreeList_.size()){

    for(size_t iTree = 0; iTree<signalTreeList_.size(); iTree ++) 
      factory_->AddSignalTree (signalTreeList_.at(iTree),signalGlobalWeight.at(iTree)) ;
  }
  else{
        
    for(size_t iTree = 0; iTree<signalTreeList_.size(); iTree ++) 
      factory_->AddSignalTree (signalTreeList_.at(iTree),1.0) ;
  }

  if(backgroundGlobalWeight.size() == backgroundTreeList_.size()){

    for(size_t iTree = 0; iTree<backgroundTreeList_.size(); iTree ++) 
      factory_->AddBackgroundTree (backgroundTreeList_.at(iTree),backgroundGlobalWeight.at(iTree)) ;
  }
  else{
        
    for(size_t iTree = 0; iTree<backgroundTreeList_.size(); iTree ++) 
      factory_->AddBackgroundTree (backgroundTreeList_.at(iTree),1.0) ;
  }

}


void TrainingMVAClass::AddPrepareTraining ( const std::string & LeptonType, const std::string & preselectionCutType, const std::string & weightString,
                                            std::vector<double> * JetPtBinOfTraining, const int & pTBin,
                                            const int & nTraining, const int & nTesting, const std::string & splitMode, const std::string & NormMode){



  if(JetPtBinOfTraining!=NULL || !JetPtBinOfTraining ){

                                  pTJetMin_ = JetPtBinOfTraining->at(pTBin) ;
                                  pTJetMax_ = JetPtBinOfTraining->at(pTBin+1) ;
  }
  else{
                                  pTJetMin_ = 0;
                                  pTJetMax_ = 2000 ;

  }

  preselectionCut_ = new TCut (GetPreselectionCut(LeptonType,preselectionCutType).Data()) ;

  TString Option = Form("nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:SplitMode=%s:NormMode=%s:!V",
                         nTraining,nTesting,nTraining,nTesting,splitMode.c_str(),NormMode.c_str());

  SetEventWeight (weightString);

  factory_->PrepareTrainingAndTestTree( *(preselectionCut_),*(preselectionCut_), Option.Data() );

}


void TrainingMVAClass::BookandTrainRectangularCuts (const std::string & FitMethod ){

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  if(FitMethod!=""){ TString Option = Form("!H:!V:FitMethod=%s:EffSel", FitMethod.c_str());
                     TString Name = Form("Cuts%s",FitMethod.c_str());
                     factory_->BookMethod( TMVA::Types::kCuts, Name.Data(),Option.Data());
  }

  else{
        factory_->BookMethod( TMVA::Types::kCuts, "CutsMC","!H:!V:FitMethod=MC:EffSel:" );
        factory_->BookMethod( TMVA::Types::kCuts, "CutsGA","!H:!V:FitMethod=GA:EffSel:Steps=40:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");
        factory_->BookMethod( TMVA::Types::kCuts, "CutsSA","!H:!V:FitMethod=SA:EffSel:KernelTemp=IncAdaptive:Eps=1e-10:UseDefaultScale" );

  }

  factory_->OptimizeAllMethods();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout<< "==> TMVAClassification is done!" << std::endl;

}


void TrainingMVAClass::BookandTrainLikelihood ( const std::string & LikelihoodType ){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  if( LikelihoodType == "LikelihoodKDE") 
    factory_->BookMethod(TMVA::Types::kLikelihood, "LikelihoodKDE","!H:!V:VarTransform=I,D,P:IgnoreNegWeightsInTraining:!TransformOutput:"
                                                                   "PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None");

  else if( LikelihoodType == "PDERS")  
      factory_->BookMethod(TMVA::Types::kPDERS, LikelihoodType.c_str(),
                           "!H:!V:VarTransform=I,D,P,G:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:DeltaFrac=4:GaussSigma=0.3:NormTree=T");

  else if( LikelihoodType == "PDEFoam")  
      factory_->BookMethod(TMVA::Types::kPDEFoam, LikelihoodType.c_str(),"!H:!V::VarTransform=I,D,P:IgnoreNegWeightsInTraining:SigBgSeparate=F:TailCut=0.001"
                                                                         ":VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");

  else if( LikelihoodType == "PDEFoamBoost")  
      factory_->BookMethod(TMVA::Types::kPDEFoam, LikelihoodType.c_str(),
                           "!H:!V::VarTransform=I,D,P:IgnoreNegWeightsInTraining:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4"
                           ":UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=300:nBin=20:Nmin=300:Kernel=None:Compress=T");

  else factory_->BookMethod( TMVA::Types::kLikelihood, LikelihoodType.c_str(),"!H:!V:VarTransform=I,D,P:!TransformOutput:IgnoreNegWeightsInTraining:PDFInterpol=Spline2"
                                                                              ":NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");

  factory_->OptimizeAllMethods();                                                                                                                                                          

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainFisherDiscriminant(){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  factory_->BookMethod( TMVA::Types::kFisher, "Fisher",
                        "!H:!V:VarTransform=I,D,P:CreateMVAPdfs:IgnoreNegWeightsInTraining:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10:Fisher" );

  factory_->OptimizeAllMethods();                                                                                                                                                          

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainLinearDiscriminant(){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Training Testing and Evaluating   
  outputFile_->cd();

  factory_->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=I,D,P:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  factory_->OptimizeAllMethods();
  
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainMLP(const int & nCycles, const std::string & HiddenLayers, const std::string & NeuronType,
				       const std::string & TrainingMethod, const int & TestRate, const int & ConvergenceTests){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  TString Option = Form ("!H:!V:VarTransform=I,D,P,D:NCycles=%d:HiddenLayers=%s:NeuronType=%s:CreateMVAPdfs:"
                         "TrainingMethod=%s:TestRate=%d:ConvergenceTests=%d:!UseRegulator",nCycles,HiddenLayers.c_str(),NeuronType.c_str(),TrainingMethod.c_str(),TestRate,ConvergenceTests);

  factory_->BookMethod( TMVA::Types::kMLP, "MLP", Option.Data());

  factory_->OptimizeAllMethods();                                                                                                                                                          
  
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}


void TrainingMVAClass::BookandTrainBDT ( const int & NTrees, const std::string & BoostType, const float & AdaBoostBeta,
				         const std::string & PruneMethod, const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  TString Option = Form ("!H:!V:VarTransform=I,D,P,D:CreateMVAPdfs:IgnoreNegWeightsInTraining:NTrees=%d:BoostType=%s:AdaBoostBeta=%f:PruneMethod=%s:"
                         "PruneStrength=%d:MaxDepth=%d:SeparationType=%s:",NTrees,BoostType.c_str(),AdaBoostBeta,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDT", Option.Data());

  //  factory_->OptimizeAllMethods();                                                                                                                                                            
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainBDTG ( const int & NTrees, const float & GradBaggingFraction, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  TString Option = Form ("!H:!V:VarTransform=I,D,P,D:CreateMVAPdfs:IgnoreNegWeights:NTrees=%d:BoostType=Grad:UseBaggedGrad:GradBaggingFraction=%f:NNodesMax=5:"
                         "PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDTG", Option.Data());

  // factory_->OptimizeAllMethods();                                                                                                                                                           

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}


void TrainingMVAClass::BookandTrainBDTF ( const int & NTrees, const std::string & BoostType, const float & AdaBoostBeta, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  TString Option = Form ("!H:!V:VarTransform=I,D,P,D:CreateMVAPdfs:UseFisherCuts:IgnoreNegWeights:NTrees=%d:BoostType=%s:AdaBoostBeta=%f"
                         "PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s",NTrees,BoostType.c_str(),AdaBoostBeta,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDTMitFisher", Option.Data());

  //  factory_->OptimizeAllMethods();
                                                                                                                                                             
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  factory_->DeleteAllMethods();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
}


//Set Methods

void TrainingMVAClass::SetSignalTree (const std::vector<TFile*> & signalFileList, const std::string & TreeName){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;

   for(size_t iFile = 0 ; iFile < signalFileList.size() ; iFile ++){
    
     if(signalFileList.at(iFile)!=0)  signalTreeList_.push_back((TTree*) signalFileList.at(iFile)->Get(TreeName_.c_str()));

     }
 }


void TrainingMVAClass::SetSignalTree (const std::vector<TTree*> & signalTreeList){
  
  if(signalTreeList.size()!=0) signalTreeList_ = signalTreeList ; 

}


void TrainingMVAClass::SetBackgroundTree (const std::vector<TFile*> & backgroundFileList, const std::string & TreeName){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;
     
  for(size_t iFile = 0 ; iFile < backgroundFileList.size() ; iFile ++){

     if(backgroundFileList.at(iFile)!=0) backgroundTreeList_.push_back((TTree*) backgroundFileList.at(iFile)->Get(TreeName_.c_str()));

     }
}

void TrainingMVAClass::SetBackgroundTree (const std::vector<TTree*> & backgroundTreeList){

   if(backgroundTreeList.size()!=0) backgroundTreeList_ = backgroundTreeList ;
}

void TrainingMVAClass::SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables){

   if(mapTrainingVariables.size()!=0) mapTrainingVariables_=mapTrainingVariables;
}

void TrainingMVAClass::SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables){

   if(mapSpectatorVariables.size()!=0) mapSpectatorVariables_=mapSpectatorVariables;
}


void TrainingMVAClass::SetLabel (const std::string & Label ){

  Label_ = Label ;
}

void TrainingMVAClass::SetTreeName (const std::string & TreeName ){

  if(TreeName!="") TreeName_ = TreeName ;
  else TreeName_ = "WJet" ;

}

void TrainingMVAClass::SetOutputFile ( const std::string & outputFilePath , const std::string & outputFileName ){
  
  if( !outputFilePath.empty() && !outputFileName.empty()) { 
 
   outputFilePath_=outputFilePath; 
   outputFileName_=outputFileName;

   outputFileNameComplete_ = outputFilePath_+"/"+outputFileName_+"_"+Label_+".root" ;

   outputFile_ = new TFile((outputFilePath_+"/"+outputFileName_+"_"+Label_+".root").c_str(),"RECREATE");
   
   outputFile_->cd();
 }

}

void TrainingMVAClass::SetGlobalSampleWeight (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight){

  signalGlobalWeight_ =  signalGlobalWeight ;
  backgroundGlobalWeight_ = backgroundGlobalWeight ;

}


void TrainingMVAClass::SetEventWeight (const std::string & weightString){

  factory_->SetWeightExpression(weightString);

}

// Take Preselection Selection 

TString TrainingMVAClass::GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType){


  if(preselectionCutType == "basicPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);
	     
  else if(preselectionCutType == "basicPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 90 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=100 && jet_mass_pr <= 130 ) )"
		" && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 90 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=100 && jet_mass_pr <= 130 ) )"
                " && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form(" issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=60 && jet_mass_pr <= 105 )"
                " && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=60 && jet_mass_pr <= 105 ) && nbjets_csvm_veto == 0 && "
                "( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);
    
  else if( preselectionCutType == "basicVBFPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("ggdboostedWevt==1 && GroomedJet_CA8_pt[0]>200 && W_pt>200 && event_met_pfmet>50 &&  W_muon_pt > 30");
      
  else if( preselectionCutType == "basicVBFPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("ggdboostedWevt==1 && GroomedJet_CA8_pt[0]>200 && W_pt>200 && event_met_pfmet>70 && W_electron_pt > 35");


  else if(preselectionCutType == "basicVBFPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("ggdboostedWevt && GroomedJet_CA8_pt[0]>200 && W_pt>200 && event_met_pfmet>50 && GroomedJet_CA8_deltaphi_METca8jet > 2.0 && GroomedJet_CA8_deltaR_lca8jet > 1.57 &&"
                " W_muon_pt > 30");

  else if(preselectionCutType == "basicVBFPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("ggdboostedWevt && GroomedJet_CA8_pt[0]>200 && W_pt>200 && event_met_pfmet>70 && GroomedJet_CA8_deltaphi_METca8jet > 2.0 && GroomedJet_CA8_deltaR_lca8jet > 1.57 &&"  
                " W_electron_pt > 35");

    else return Form("v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

}



// print Training results 

void  TrainingMVAClass::PrintTrainingResults (){

  std::string command = " if [ ! -e plots ] ; then mkdir plots ; fi";
  system(command.c_str());


  std::cout << "******************************************************* "<<std::endl;  
  std::cout << "==> Print Output Plots For: " << outputFile_->GetName() << std::endl;
  std::cout << "******************************************************* "<<std::endl;  

  std::string ROOTStyle =  getenv ("ROOTStyle");

  gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());

  std::string PWD = getenv("PWD");
 
  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/variables.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  std::cout<<command<<std::endl;
  system(command.c_str());

  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/correlationscatters.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  system(command.c_str());

  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/correlations.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  system(command.c_str());

  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/mvas.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  system(command.c_str());

  command = "root -l -q -b "+PWD+"/macros/TMVAMacro/efficiencies.C\\(\\\""+PWD+"/"+outputFileNameComplete_+"\\\"\\)";
  system(command.c_str());

  command = " if [ ! -e "+outputFilePath_+"/trainingPlots ] ; then mkdir "+outputFilePath_+"/trainingPlots ; fi";
  system(command.c_str());

  command = " if [ ! -e "+outputFilePath_+"/trainingPlots/"+outputFileName_+"_"+Label_+" ] ; then mkdir "+outputFilePath_+"/trainingPlots/"+
                         outputFileName_+"_"+Label_+" ; fi";
  system(command.c_str());
  
  command = "mv ./plots/* "+outputFilePath_+"/trainingPlots/"+outputFileName_+"_"+Label_+"/" ;
  system(command.c_str());

  command = "mv weights/ "+outputFilePath_;
  system(command.c_str());

  command = "rm -rf ./plots" ;
  system(command.c_str());

  std::cout << "==> Wrote image files: " << outputFilePath_+"/trainingPlots/"+outputFileName_+"_"+Label_ << std::endl;
  std::cout << "==> TMVA Plots are done!" << std::endl;
   
}


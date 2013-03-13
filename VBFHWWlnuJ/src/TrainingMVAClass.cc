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

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ )
    factory_->AddVariable(mapTrainingVariables_.at(iVar)+" := "+mapTrainingVariables_.at(iVar),'F');
   
  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ )
    factory_->AddSpectator(mapSpectatorVariables_.at(iVar),'F');
    
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
                                            const int & nTraining, const int & nTesting, const std::string & splitMode, const std::string & NormMode){

  preselectionCut_ = new TCut (GetPreselectionCut(LeptonType,preselectionCutType).c_str()) ;

  TString Option = Form("nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:SplitMode=%s:NormMode=%s:!V",
                         nTraining,nTesting,nTraining,nTesting,splitMode.c_str(),NormMode.c_str());

  SetEventWeight (weightString);

  factory_->PrepareTrainingAndTestTree( *(preselectionCut_),*(preselectionCut_), Option.Data() );

}


void TrainingMVAClass::BookandTrainRectangularCuts (const std::string & FitMethod ){

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure
  outputFileWeightName_["Cuts"+FitMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_Cuts"+FitMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["Cuts"+FitMethod+"_"+Label_];

  // Training Testing and Evaluating 
  outputFile_->cd();

  if(FitMethod!=""){ TString Option = Form("!H:!V:CreateMVAPdfs:FitMethod=%s:EffSel", FitMethod.c_str());
                     TString Name = Form("Cuts%s",FitMethod.c_str());
                     factory_->BookMethod( TMVA::Types::kCuts, Name.Data(),Option.Data());
  }

  else{
        factory_->BookMethod( TMVA::Types::kCuts, "CutsMC","!H:!V:CreateMVAPdfs:FitMethod=MC:EffSel:" );
        factory_->BookMethod( TMVA::Types::kCuts, "CutsGA","!H:!V:CreateMVAPdfs:FitMethod=GA:EffSel:Steps=40:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");
        factory_->BookMethod( TMVA::Types::kCuts, "CutsSA","!H:!V:CreateMVAPdfs:FitMethod=SA:EffSel:KernelTemp=IncAdaptive:Eps=1e-10:UseDefaultScale" );

  }

  factory_->OptimizeAllMethods();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout<< "==> TMVAClassification is done!" << std::endl;

}


void TrainingMVAClass::BookandTrainLikelihood ( const std::string & LikelihoodType ){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure
  outputFileWeightName_[LikelihoodType+"_"+Label_] = outputFilePath_+"/TMVAWeight_"+LikelihoodType+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_[LikelihoodType+"_"+Label_];

  // Training Testing and Evaluating 
  outputFile_->cd();

  if( LikelihoodType == "LikelihoodKDE") 
    factory_->BookMethod(TMVA::Types::kLikelihood, "LikelihoodKDE","!H:!V:!VarTransform=I,D,P,G:TransformOutput::CreateMVAPdfs:IgnoreNegWeightsInTraining:"
                                                                   "PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=5");

  else if( LikelihoodType == "PDERS")  
      factory_->BookMethod(TMVA::Types::kPDERS, LikelihoodType.c_str(),
                           "!H:!V:VarTransforms=I;D;P;G:CreateMVAPdfs:IgnoreNegWeightsInTraining:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:DeltaFrac=4:GaussSigma=0.3:NormTree=T");

  else if( LikelihoodType == "PDEFoam")  
      factory_->BookMethod(TMVA::Types::kPDEFoam, LikelihoodType.c_str(),"!H:!V::VarTransform=I,D,P,G:CreateMVAPdfs:IgnoreNegWeightsInTraining:SigBgSeparate=F:TailCut=0.001"
                                                                         ":VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");

  else if( LikelihoodType == "PDEFoamBoost")  
      factory_->BookMethod(TMVA::Types::kPDEFoam, LikelihoodType.c_str(),
                           "!H:!V::VarTransforms=I;D;P;G:CreateMVAPdfs:IgnoreNegWeightsInTraining:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4"
                           ":UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T");

  else factory_->BookMethod( TMVA::Types::kLikelihood, LikelihoodType.c_str(),"!H:!V:VarTransform=I,D,P,G:CreateMVAPdfs:IgnoreNegWeightsInTraining:!TransformOutput:PDFInterpol=Spline2"
                                                                              ":NAvEvtPeDrBin=50");

  factory_->OptimizeAllMethods();                                                                                                                                                          

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainFisherDiscriminant(){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure                                                              

  outputFileWeightName_["Fisher"+Label_] = outputFilePath_+"/TMVAWeight_Fisher_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["Fisher"+Label_];

  // Training Testing and Evaluating  
  outputFile_->cd();

  factory_->BookMethod( TMVA::Types::kFisher, "Fisher",
                        "!H:!V:VarTransform=I,D,P,G,D:CreateMVAPdfs:IgnoreNegWeightsInTraining:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10:Fisher" );

  factory_->OptimizeAllMethods();                                                                                                                                                          

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainLinearDiscriminant(){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure

  outputFileWeightName_["LD"+Label_] = outputFilePath_+"/TMVAWeight_LD_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["LD"+Label_];

  // Training Testing and Evaluating   
  outputFile_->cd();

  factory_->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=I,D,P,G,D:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  factory_->OptimizeAllMethods();
  
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainMLP(const int & nCycles, const std::string & HiddenLayers, const std::string & NeuronType,
				       const std::string & TrainingMethod, const int & TestRate, const int & ConvergenceTests){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                                

  outputFileWeightName_["MLP_"+NeuronType+"_"+TrainingMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_MLP_"+NeuronType+"_"+TrainingMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["MLP_"+NeuronType+"_"+TrainingMethod+"_"+Label_];

  // Training Testing and Evaluating                                                  
  outputFile_->cd();

  TString Option = Form ("!H:!V:VarTransform=I,D,P,G:CreateMVAPdfs:IgnoreNegWeightsInTraining:NCycles=%d:HiddenLayers=%s:NeuronType=%s:"
                         "TrainingMethod=%s:TestRate=%d:ConvergenceTests=%d:!UseRegulator",nCycles,HiddenLayers.c_str(),NeuronType.c_str(),TrainingMethod.c_str(),TestRate,ConvergenceTests);

  factory_->BookMethod( TMVA::Types::kMLP, "MLP", Option.Data());

  factory_->OptimizeAllMethods();                                                                                                                                                          
  
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}


void TrainingMVAClass::BookandTrainBDT ( const int & NTrees, const std::string & BoostType, const float & AdaBoostBeta,
				         const std::string & PruneMethod, const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                                 

  outputFileWeightName_["BDT_"+BoostType+"_"+PruneMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_BDT_"+BoostType+"_"+PruneMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["BDT_"+BoostType+"_"+PruneMethod+"_"+Label_];

  // Training Testing and Evaluating                                                                                                                                           
  outputFile_->cd();

  TString Option = Form ("!H:!V:VarTransform=I,D,P,G:CreateMVAPdfs:IgnoreNegWeightsInTraining:NTrees=%d:BoostType=%s:AdaBoostBeta=%f:PruneMethod=%s:"
                         "PruneStrength=%d:MaxDepth=%d:SeparationType=%s",NTrees,BoostType.c_str(),AdaBoostBeta,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDT", Option.Data());

  factory_->OptimizeAllMethods();                                                                                                                                                            

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}

void TrainingMVAClass::BookandTrainBDTG ( const int & NTrees, const float & GradBaggingFraction, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                                

  outputFileWeightName_["BDTG_"+PruneMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_BDTG_"+PruneMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["BDTG_"+PruneMethod+"_"+Label_];

  // Training Testing and Evaluating 
  outputFile_->cd();

  TString Option = Form ("!H:!V:VarTransform=I,D,P,G:CreateMVAPdfs:IgnoreNegWeightsInTraining:NTrees=%d:BoostType=Grad:UseBaggedGrad:GradBaggingFraction=%f:"
                         "PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDTG", Option.Data());

  factory_->OptimizeAllMethods();                                                                                                                                                           

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

  std::cout << "==> Wrote root file: " << outputFile_->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

}


void TrainingMVAClass::BookandTrainBDTF ( const int & NTrees, const float & GradBaggingFraction, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure                                                                                                                                 
  outputFileWeightName_["BDTF_"+PruneMethod+"_"+Label_] = outputFilePath_+"/TMVAWeight_BDTF_"+PruneMethod+"_"+Label_;
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_["BDTF_"+PruneMethod+"_"+Label_];

  // Training Testing and Evaluating 
  outputFile_->cd();

  TString Option = Form ("!H:!V:VarTransform=I,D,P,G:CreateMVAPdfs:IgnoreNegWeightsInTraining:UseFisherCuts:NTrees=%d:BoostType=Grad:UseBaggedGrad:GradBaggingFraction=%f:"
                         "PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDTF", Option.Data());

  factory_->OptimizeAllMethods();                                                                                                                                                             
  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();

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

std::string TrainingMVAClass::GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType){

  if( preselectionCutType == "basicPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return "issignal && v_pt > 250 && pfMET > 50 && l_pt > 30 && ungroomed_jet_pt > 250" ;

  else if(preselectionCutType == "basicPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return "issignal && v_pt > 250 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 250" ;

  else if(preselectionCutType == "basicSRPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return "issignal && v_pt > 250 && pfMET > 50 && l_pt > 30 && ungroomed_jet_pt > 250 && ( jet_mass_pr >=65 && jet_mass_pr <= 100 )";

  else if(preselectionCutType == "basicSRPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return "issignal && v_pt > 250 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 250 && ( jet_mass_pr >=65 && jet_mass_pr <= 100 )";

  else return "v_pt > 250 && pfMET > 50 && l_pt > 30 && ungroomed_jet_pt > 250" ;

}



// print Training results 

void  TrainingMVAClass::PrintTrainingResults (){

  std::string command = " if [ ! -e plots ] ; then mkdir plots ; fi";
  system(command.c_str());


  std::cout << "******************************************************* "<<std::endl;  
  std::cout << "==> Print Output Plots For: " << outputFile_->GetName() << std::endl;
  std::cout << "******************************************************* "<<std::endl;  


  command = "root -l -b -q ../macros/TMVAMacros/variables.C("+outputFileNameComplete_+")";
  gROOT->ProcessLine(command.c_str());

  command = "root -l -b -q ../macros/TMVAMacros/correlationscatter.C("+outputFileNameComplete_+")";
  gROOT->ProcessLine(command.c_str());

  command = "root -l -b -q ../macros/TMVAMacros/correlations.C("+outputFileNameComplete_+")";
  gROOT->ProcessLine(command.c_str());

  command = "root -l -b -q ../macros/TMVAMacros/mvas.C("+outputFileNameComplete_+")";
  gROOT->ProcessLine(command.c_str());

  command = "root -l -b -q ../macros/TMVAMacros/mvaeffs.C("+outputFileNameComplete_+")";
  gROOT->ProcessLine(command.c_str());

  command = "root -l -b -q ../macros/TMVAMacros/efficiencies.C("+outputFileNameComplete_+")";
  gROOT->ProcessLine(command.c_str());

  command = " if [ ! -e "+outputFilePath_+"/trainingPlots ] ; then mkdir"+outputFilePath_+"/trainingPlots ; fi";
  system(command.c_str());

  command = " if [ ! -e "+outputFilePath_+"/trainingPlots/"+outputFileNameComplete_+" ] ; then mkdir"+outputFilePath_+"/trainingPlots/"+
                         outputFileNameComplete_+" ; fi";
  system(command.c_str());

  command = "mv ./plots/ "+outputFilePath_+"/trainingPlots/"+outputFileNameComplete_+"/" ;
  system(command.c_str());

  std::cout << "==> Wrote image files: " << outputFilePath_+"/trainingPlots/"+outputFileNameComplete_ << std::endl;
  std::cout << "==> TMVA Plots are done!" << std::endl;
   
}


void  TrainingMVAClass::ReadWeightFromXML ( const std::string & LeptonType, const std::string & preselectionCutType){


  std::cout << "***************************************************************** "<<std::endl;  
  std::cout << "==> Start Read MVA Ouput related to : " << outputFile_->GetName() << std::endl;
  std::cout << "***************************************************************** "<<std::endl;  
   

  reader_ = new TMVA::Reader( "!Color:!Silent" );

  setTrainingVariables_  = new std::vector<Float_t> ( mapTrainingVariables_.size()) ;
  setSpectatorVariables_ = new std::vector<Float_t> ( mapSpectatorVariables_.size()) ;

  // Set Cut to apply on the events read from initial tree

  int CutType = 0;

  if   (preselectionCutType == "basicPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )             CutType = 0;
      
  else if(preselectionCutType == "basicPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )   CutType = 1;
	
  else if(preselectionCutType == "basicSRPreselectionCut" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )         CutType = 2;
	
  else if(preselectionCutType == "basicSRPreselectionCut" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") ) CutType = 3 ;
	
  else CutType = 4;
 

  // Set Input and Spectator variables

  for( size_t iVar = 0; iVar< mapTrainingVariables_.size(); iVar++)
    reader_->AddVariable (mapTrainingVariables_.at(iVar)+" := "+mapTrainingVariables_.at(iVar), &setTrainingVariables_->at(iVar));

  for( size_t iVar = 0; iVar< mapSpectatorVariables_.size(); iVar++)
    reader_->AddSpectator (mapSpectatorVariables_.at(iVar)+" := "+mapSpectatorVariables_.at(iVar), &setSpectatorVariables_->at(iVar));

  
  // Fill MVA Value for Signal Tree

  for(size_t iSignal = 0 ; iSignal < signalTreeList_.size(); iSignal++){

   treeReader_ = new treeReader ((TTree*) signalTreeList_.at(iSignal), false) ;

   std::map<std::string,std::string>::const_iterator itMap = outputFileWeightName_.begin();

   int iMethod = 0;
  
   std::vector<float> weight (outputFileWeightName_.size(),0.);

   // Loop on all the training method used
   for( ; itMap != outputFileWeightName_.end() ; ++itMap){

    reader_->BookMVA(itMap->first,itMap->second+".xml");

    std::cout<<" Read for signal weight file : "<<itMap->second+".xml"<<std::endl;

    // Set new Branches for MVA output
    newBranches_.at(iMethod) =  signalTreeList_.at(iSignal)->Branch((itMap->first).c_str(),&weight.at(iMethod),(itMap->first+"/F").c_str());    

    for( int iEntry = 0; iEntry < signalTreeList_.at(iSignal)->GetEntries(); iEntry++){

      if (iEntry % 10000 == 0) std::cout << "reading event " << iEntry << std::endl ;

      signalTreeList_.at(iSignal)->GetEntry(iEntry);
    
      if( CutType == 0 ) {
	if( treeReader_->getFloat("issignal")[0] == 0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 50 || 
            treeReader_->getFloat("l_pt")[0] < 30 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill() ; }
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else if( CutType == 1 ) {
	if( treeReader_->getFloat("issignal")[0] == 0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 70 || 
            treeReader_->getFloat("l_pt")[0] < 35 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill(); }
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else if( CutType == 2 ) {
	if( treeReader_->getFloat("issignal")[0] == 0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 50 || 
            treeReader_->getFloat("l_pt")[0] < 30 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 || 
            (treeReader_->getFloat("jet_mass_pr")[0] <=65 || treeReader_->getFloat("jet_mass_pr")[0] >=100) ){ weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill() ;}
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else if( CutType == 3 ) {
      if( treeReader_->getFloat("issignal")[0] ==0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 70 || 
          treeReader_->getFloat("l_pt")[0] < 35 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 || 
          (treeReader_->getFloat("jet_mass_pr")[0] <=65 || treeReader_->getFloat("jet_mass_pr")[0] >=100) ){ weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill() ;}
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else {
	if( treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 50 || 
            treeReader_->getFloat("l_pt")[0] < 30 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill(); }
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      

    } // end event Loop
    
    iMethod++; 

   } // end method Loop

    delete treeReader_ ;

    signalTreeList_.at(iSignal)->Write("", TObject::kOverwrite);
    
  } // end signal tree loop

  // Fill MVA Value for Background Tree

  for(size_t iBack = 0 ; iBack < backgroundTreeList_.size(); iBack++){

   treeReader_ = new treeReader ((TTree*) backgroundTreeList_.at(iBack), false) ;

   std::map<std::string,std::string>::const_iterator itMap = outputFileWeightName_.begin();

   int iMethod = 0;
  
   std::vector<float> weight (outputFileWeightName_.size(),0.);

   // Loop on all the training method used
   for( ; itMap != outputFileWeightName_.end() ; ++itMap){

    reader_->BookMVA(itMap->first,itMap->second+".xml");

    std::cout<<" Read for background weight file : "<<itMap->second+".xml"<<std::endl;

    // Set new Branches for MVA output
    newBranches_.at(iMethod) =  backgroundTreeList_.at(iBack)->Branch((itMap->first).c_str(),&weight.at(iMethod),(itMap->first+"/F").c_str());    

    for( int iEntry = 0; iEntry < backgroundTreeList_.at(iBack)->GetEntries(); iEntry++){

      if (iEntry % 10000 == 0) std::cout << "reading event " << iEntry << std::endl ;

      backgroundTreeList_.at(iBack)->GetEntry(iEntry);
    
      if( CutType == 0 ) {
	if( treeReader_->getFloat("issignal")[0] == 0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 50 || 
            treeReader_->getFloat("l_pt")[0] < 30 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill() ; }
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else if( CutType == 1 ) {
	if( treeReader_->getFloat("issignal")[0] == 0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 70 || 
            treeReader_->getFloat("l_pt")[0] < 35 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill(); }
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else if( CutType == 2 ) {
	if( treeReader_->getFloat("issignal")[0] == 0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 50 || 
            treeReader_->getFloat("l_pt")[0] < 30 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 || 
            (treeReader_->getFloat("jet_mass_pr")[0] <=65 || treeReader_->getFloat("jet_mass_pr")[0] >=100) ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill() ;}
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else if( CutType == 3 ) {
      if( treeReader_->getFloat("issignal")[0] ==0 || treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 70 || 
          treeReader_->getFloat("l_pt")[0] < 35 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 || 
          (treeReader_->getFloat("jet_mass_pr")[0] <=65 || treeReader_->getFloat("jet_mass_pr")[0] >=100) ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill() ;}
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      else {
	if( treeReader_->getFloat("v_pt")[0] < 250 || treeReader_->getFloat("pfMET")[0] < 50 || 
            treeReader_->getFloat("l_pt")[0] < 30 || treeReader_->getFloat("ungroomed_jet_pt")[0] < 250 ) { weight.at(iMethod) = -100 ; newBranches_.at(iMethod)->Fill(); }
        else{ weight.at(iMethod) = reader_->EvaluateMVA(itMap->first); newBranches_.at(iMethod)->Fill() ;}
      }
      

    } // end event Loop
    
    iMethod++; 

   } // end method Loop

    delete treeReader_ ;

    backgroundTreeList_.at(iBack)->Write("", TObject::kOverwrite);
    
  } // end signal tree loop

}

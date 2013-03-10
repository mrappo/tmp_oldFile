#include "TrainingMVAClass.h"

// Constructor

TrainingMVAClass::TrainingMVAClass(const std::vector<TFile*> & signalFileList, const std::vector<TFile*> & backgroundFileList, const std::string & TreeName,
				   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label ){

   SetSignalTree (signalFileList,TreeName) ;

   SetBackgroundTree (backgroundFileList,TreeName) ;

   SetLabel ( Label );

   SetOutputFile ( outputFilePath , outputFileName ) ;

   factory_ = new TMVA::Factory (TreeName_+"_"+Label_,outputFile_, Form("!V:!Silent:%sColor:DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification", gROOT->IsBatch()?"!":""));

}

TrainingMVAClass::TrainingMVAClass(const std::vector<TTree*> & signalTreeList, const std::vector<TTree*> & backgroundTreeList,  const std::string & TreeName,
                                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label ){
   

   SetTreeName (TreeName) ;

   SetSignalTree (signalTreeList) ;

   SetBackgroundTree (backgroundTreeList) ;

   SetLabel ( Label );

   SetOutputFile ( outputFilePath , outputFileName ) ;

   factory_ = new TMVA::Factory (TreeName_+"_"+Label_,outputFile_, Form("!V:!Silent:%sColor:DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification", gROOT->IsBatch()?"!":""));

}

// Deconstructor

TrainingMVAClass::~TrainingMVAClass(){

  for(size_t iTree = 0; iTree < signalTreeList_.size() ; iTree++) { if(signalTreeList_.at(iTree)!=0)  delete signalTreeList_.at(iTree) ; }

  for(size_t iTree = 0; iTree < backgroundTreeList_.size() ; iTree++) {if(backgroundTreeList_.at(iTree)!=0)  delete backgroundTreeList_.at(iTree) ;}

  if(outputFile_!=0) outputFile_->Delete() ;

  if(preselectionCut_!=0) preselectionCut_->Delete() ;
 
  if(factory_!=0) factory_->Delete() ;

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


void TrainingMVAClass::AddPrepareTraining ( const std::string & cutString, const std::string & weightString, 
                                            const int & nTraining, const int & nTesting, const std::string & splitMode, const std::string & NormMode){

  preselectionCut_ = new TCut (cutString.c_str()) ;

  TString Option = Form("nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:SplitMode=%s:NormMode=%s:!V",
                         nTraining,nTesting,nTraining,nTesting,splitMode.c_str(),NormMode.c_str());

  SetEventWeight (weightString);

  factory_->PrepareTrainingAndTestTree( *(preselectionCut_),*(preselectionCut_), Option.Data() );

}


void TrainingMVAClass::BookandTrainRectangularCuts (const std::string & FitMethod ){

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  // Set Name of the Weight file for TMVA evaluating procedure
  outputFileWeightName_.push_back("TMVAWeight_Cuts"+FitMethod+"_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

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
  outputFileWeightName_.push_back("TMVAWeight_"+LikelihoodType+"_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

  // Training Testing and Evaluating 
  outputFile_->cd();

  if( LikelihoodType == "LikelihoodKDE") 
      factory_->BookMethod(TMVA::Types::kLikelihood, "LikelihoodKDE","!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=5");

  else if( LikelihoodType == "PDERS")  
      factory_->BookMethod(TMVA::Types::kPDERS, LikelihoodType.c_str(),
                           "!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:DeltaFrac=4:GaussSigma=0.3:NormTree=T");

  else factory_->BookMethod( TMVA::Types::kLikelihood, LikelihoodType.c_str(),"!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:!TransformOutput:PDFInterpol=Spline2:NAvEvtPeDrBin=50");

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

  outputFileWeightName_.push_back("TMVAWeight_Fisher_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

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

  outputFileWeightName_.push_back("TMVAWeight_LD_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

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

  outputFileWeightName_.push_back("TMVAWeight_MLP_"+NeuronType+"_"+TrainingMethod+"_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

  // Training Testing and Evaluating                                                  
  outputFile_->cd();

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:NCycles=%d:HiddenLayers=%s:NeuronType=%s:TrainingMethod=%s:TestRate=%d:ConvergenceTests=%d:!UseRegulator",
                         nCycles,HiddenLayers.c_str(),NeuronType.c_str(),TrainingMethod.c_str(),TestRate,ConvergenceTests);

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

  outputFileWeightName_.push_back("TMVAWeight_BDT_"+BoostType+"_"+PruneMethod+"_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

  // Training Testing and Evaluating                                                                                                                                           
  outputFile_->cd();

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:NTrees=%d:BoostType=%s:AdaBoostBeta=%f:PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s",
                         NTrees,BoostType.c_str(),AdaBoostBeta,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

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

  outputFileWeightName_.push_back("TMVAWeight_BDTG_"+PruneMethod+"_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

  // Training Testing and Evaluating 
  outputFile_->cd();

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:NTrees=%d:BoostType=Grad:UseBaggedGrad:GradBaggingFraction=%f:PruneMethod=%s:PruneStrength=%d:MaxDepth=%d"
                         ":SeparationType=%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

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
  outputFileWeightName_.push_back("TMVAWeight_BDTF_"+PruneMethod+"_"+Label_);
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outputFileWeightName_.back();

  // Training Testing and Evaluating 
  outputFile_->cd();

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:UseFisherCuts:NTrees=%d:BoostType=Grad:UseBaggedGrad:GradBaggingFraction=%f:PruneMethod=%s:PruneStrength=%d"
                         ":MaxDepth=%d:SeparationType=%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

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

  command = "mv ./plots/* "+outputFilePath_+"/trainingPlots/"+outputFileNameComplete_+"/" ;
  system(command.c_str());

  std::cout << "==> Wrote image files: " << outputFilePath_+"/trainingPlots/"+outputFileNameComplete_ << std::endl;
  std::cout << "==> TMVA Plots are done!" << std::endl;
   
}

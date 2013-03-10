#include "TrainingMVAClass.h"

// Constructor

TrainingMVAClass::TrainingMVAClass(const std::vector<TFile*> & signalFileList, const std::vector<TFile*> & backgroundFileList, const std::string & TreeName,
				   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label ){

   SetSignalTree (signalFileList,TreeName) ;

   SetBackgroundTree (backgroundFileList,TreeName) ;

   SetOutputFile ( outputFilePath , outputFileName ) ;

   SetLabel ( Label );
 
   factory_ = new TMVA::Factory ((TreeName_+Label).c_str(),outputFile_,"!V:!Silent:Color:DrawProgressBar:Transformations=I;N;D;P;G,D:AnalysisType=Classification");
    
}

TrainingMVAClass::TrainingMVAClass(const std::vector<TTree*> & signalTreeList, const std::vector<TTree*> & backgroundTreeList,  const std::string & TreeName,
                                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label ){
   

   SetTreeName (TreeName) ;

   SetSignalTree (signalTreeList) ;

   SetBackgroundTree (backgroundTreeList) ;

   SetOutputFile ( outputFilePath , outputFileName ) ;

   SetLabel ( Label );
 
   factory_ = new TMVA::Factory ((TreeName_+Label).c_str(),outputFile_,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

}

// Deconstructor

TrainingMVAClass::~TrainingMVAClass(){

  for(size_t iTree = 0; iTree < signalTreeList_.size() ; iTree++) { if(signalTreeList_.at(iTree)!=0)  delete signalTreeList_.at(iTree) ; }

  for(size_t iTree = 0; iTree < backgroundTreeList_.size() ; iTree++) {if(backgroundTreeList_.at(iTree)!=0)  delete backgroundTreeList_.at(iTree) ;}

  if(outputFile_!=0) outputFile_->Delete() ;

  if(preselectionCut_!=0) preselectionCut_->Delete() ;
 
  if(factory_!=0) factory_->Delete() ;

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

// AddTrainingVariables in the MVA

void TrainingMVAClass::AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables, const std::string & weightString ){

  SetTrainingVariables(mapTrainingVariables);
  SetSpectatorVariables(mapSpectatorVariables);

  for( size_t iVar = 0 ; iVar < mapTrainingVariables_.size() ; iVar ++ )
    factory_->AddVariable((mapTrainingVariables_.at(iVar)+" := "+mapTrainingVariables_.at(iVar)).c_str(), 'F');

  for( size_t iVar = 0 ; iVar < mapSpectatorVariables_.size() ; iVar ++ )
    factory_->AddSpectator(mapSpectatorVariables_.at(iVar).c_str());

  SetEventWeight(weightString);

}

void TrainingMVAClass::AddPrepareTraining ( const std::string & cutString, const int & nTraining, const int & nTesting, const std::string & splitMode, const std::string & NormMode ){

  preselectionCut_ = new TCut (cutString.c_str()) ;

  TString Option = Form("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:SplitMode=%s:NormMode=%s:!V",
                         nTraining,nTesting,nTraining,nTesting,splitMode.c_str(),NormMode.c_str());

  factory_->PrepareTrainingAndTestTree( *(preselectionCut_), Option.Data() );

}


void TrainingMVAClass::BookandTrainRectangularCuts (const std::string & FitMethod ){

  if(FitMethod!=""){ TString Option = Form("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:FitMethod=%s:EffSel", FitMethod.c_str());
                     TString Name = Form("Cuts%s",FitMethod.c_str());
                     factory_->BookMethod( TMVA::Types::kCuts, Name.Data(),Option.Data());
  }

  else{
        factory_->BookMethod( TMVA::Types::kCuts, "CutsMC","!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:FitMethod=MC:EffSel:" );
        factory_->BookMethod( TMVA::Types::kCuts, "CutsGA","!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:FitMethod=GA:EffSel:Steps=40:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");
        factory_->BookMethod( TMVA::Types::kCuts, "CutsSA","!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:FitMethod=SA:EffSel:KernelTemp=IncAdaptive:Eps=1e-10:UseDefaultScale" );

  }

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+"Cuts"+FitMethod+".root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+"Cuts"+FitMethod);

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

}


void TrainingMVAClass::BookandTrainLikelihood ( const std::string & LikelihoodType ){


  if( LikelihoodType == "LikelihoodKDE") 
      factory_->BookMethod(TMVA::Types::kLikelihood, "LikelihoodKDE","!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=5");

  else if( LikelihoodType == "PDERS")  
      factory_->BookMethod(TMVA::Types::kPDERS, LikelihoodType.c_str(),
                           "!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:DeltaFrac=4:GaussSigma=0.3:NormTree=T");

  else factory_->BookMethod( TMVA::Types::kLikelihood, LikelihoodType.c_str(),"!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:!TransformOutput:PDFInterpol=Spline2:NAvEvtPeDrBin=50");


  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());
  
  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+LikelihoodType+".root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+LikelihoodType);

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

}

void TrainingMVAClass::BookandTrainFisherDiscriminant(){


  factory_->BookMethod( TMVA::Types::kFisher, "Fisher", "!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10:Fisher" );

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());


  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+"Fisher.root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+"Fisher");

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

}

void TrainingMVAClass::BookandTrainLinearDiscriminant(){


  factory_->BookMethod( TMVA::Types::kLD, "LD", "!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+"LD.root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+"LD");

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

}

void TrainingMVAClass::BookandTrainMLP(const int & nCycles, const std::string & HiddenLayers, const std::string & NeuronType,
				       const std::string & TrainingMethod, const int & TestRate, const int & ConvergenceTests){

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:NCycles=%d:HiddenLayers=%s:NeuronType=%s:TrainingMethod=%s:TestRate=%d:ConvergenceTests=%d:!UseRegulator",
                         nCycles,HiddenLayers.c_str(),NeuronType.c_str(),TrainingMethod.c_str(),TestRate,ConvergenceTests);

  factory_->BookMethod( TMVA::Types::kMLP, "MLP", Option.Data());

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());


  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+"MLP.root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+"MLP");

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

}


void TrainingMVAClass::BookandTrainBDT ( const int & NTrees, const std::string & BoostType, const float & AdaBoostBeta,
				         const std::string & PruneMethod, const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:NTrees=%d:BoostType=%s:AdaBoostBeta=%f:PruneMethod=%s:PruneStrength=%d:MaxDepth=%d:SeparationType=%s",
                         NTrees,BoostType.c_str(),AdaBoostBeta,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDT", Option.Data());

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());


  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+"BDT.root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+"BDT");

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

}

void TrainingMVAClass::BookandTrainBDTG ( const int & NTrees, const float & GradBaggingFraction, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:NTrees=%d:BoostType=Grad:UseBaggedGrad:GradBaggingFraction=%f:PruneMethod=%s:PruneStrength=%d:MaxDepth=%d"
                         ":SeparationType=%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDTG", Option.Data());

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+"BDTG.root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+"BDTG");

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

}


void TrainingMVAClass::BookandTrainBDTF ( const int & NTrees, const float & GradBaggingFraction, const std::string & PruneMethod,
                       			  const int & PruneStrength, const int & MaxDepth, const std::string & SeparationType){

  TString Option = Form ("!H:!V:CreateMVAPdfs:IgnoreNegWeightsInTraining:UseFisherCuts:NTrees=%d:BoostType=Grad:UseBaggedGrad:GradBaggingFraction=%f:PruneMethod=%s:PruneStrength=%d"
                         ":MaxDepth=%d:SeparationType=%s",NTrees,GradBaggingFraction,PruneMethod.c_str(),PruneStrength,MaxDepth,SeparationType.c_str());

  factory_->BookMethod( TMVA::Types::kBDT, "BDTF", Option.Data());

  std::string command = " if [ ! -e "+outputFilePath_+" ] ; then mkdir "+outputFilePath_+" ; fi";
  system(command.c_str());

  outputFile_ = new TFile ((outputFilePath_+"/"+outputFileName_+"BDTF.root").c_str(),"RECREATE");

  outputFileNameComplete_.push_back(outputFilePath_+"/"+outputFileName_+"BDTF");

  outputFile_->cd();

  factory_->TrainAllMethods();

  factory_->TestAllMethods();

  factory_->EvaluateAllMethods();

  outputFile_->Close();
  outputFile_->Delete();

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
   
   }
}

void TrainingMVAClass::SetGlobalSampleWeight (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight){

  signalGlobalWeight_ =  signalGlobalWeight ;
  backgroundGlobalWeight_ = backgroundGlobalWeight ;

}


void TrainingMVAClass::SetEventWeight (const std::string & weightString){

  factory_->SetWeightExpression(weightString.c_str());

}

// print Training results 

void  TrainingMVAClass::PrintTrainingResults (){

  std::string command = " if [ ! -e plots ] ; then mkdir plots ; fi";
  system(command.c_str());


  for(size_t iRootFile = 0 ; iRootFile < outputFileNameComplete_.size() ; iRootFile++){

    command = "root -l -b -q ../macros/TMVAMacros/variables.C("+outputFileNameComplete_.at(iRootFile)+".root)";
    gROOT->ProcessLine(command.c_str());

    command = "root -l -b -q ../macros/TMVAMacros/correlationscatter.C("+outputFileNameComplete_.at(iRootFile)+".root)";
    gROOT->ProcessLine(command.c_str());

    command = "root -l -b -q ../macros/TMVAMacros/correlations.C("+outputFileNameComplete_.at(iRootFile)+".root)";
    gROOT->ProcessLine(command.c_str());

    command = "root -l -b -q ../macros/TMVAMacros/mvas.C("+outputFileNameComplete_.at(iRootFile)+".root)";
    gROOT->ProcessLine(command.c_str());

    command = "root -l -b -q ../macros/TMVAMacros/mvaeffs.C("+outputFileNameComplete_.at(iRootFile)+".root)";
    gROOT->ProcessLine(command.c_str());

    command = "root -l -b -q ../macros/TMVAMacros/efficiencies.C("+outputFileNameComplete_.at(iRootFile)+".root)";
    gROOT->ProcessLine(command.c_str());

    command = " if [ ! -e "+outputFilePath_+"/trainingPlots ] ; then mkdir"+outputFilePath_+"/trainingPlots ; fi";
    system(command.c_str());

    command = " if [ ! -e "+outputFilePath_+"/trainingPlots/"+outputFileNameComplete_.at(iRootFile)+" ] ; then mkdir"+outputFilePath_+"/trainingPlots/"+
                            outputFileNameComplete_.at(iRootFile)+" ; fi";
    system(command.c_str());

    command = "mv ./plots/"+outputFileNameComplete_.at(iRootFile)+".root* "+outputFilePath_+"/trainingPlots/"+outputFileNameComplete_.at(iRootFile)+"/" ;
    system(command.c_str());

  }


  
}

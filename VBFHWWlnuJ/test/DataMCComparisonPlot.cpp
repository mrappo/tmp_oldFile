#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TLatex.h"

#include "ntpleUtils.h"
#include "ConfigParser.h"
#include "ReadInputFile.h"
#include "DataMCPlotTool.h"

void LatexCMS (double lumi) ;

// To lunch the sequence for the Higgs : 

// 1) Ele : ls cfg/DataMCComparison_InputCfgFile/ | grep _el_ | grep HWW | grep cfg | grep -v ~ | awk '{print "./bin/DataMCComparisonPlot.exe cfg/DataMCComparison_InputCfgFile/"$1}' | /bin/sh
// 2) Mu  : ls cfg/DataMCComparison_InputCfgFile/ | grep _mu_ | grep HWW | grep cfg | grep -v ~ | awk '{print "./bin/DataMCComparisonPlot.exe cfg/DataMCComparison_InputCfgFile/"$1}' | /bin/sh


// To lunch the sequence for the RSG Graviton : 

// 1) Ele : ls cfg/DataMCComparison_InputCfgFile/ | grep _el_ | grep RSG | grep cfg | grep -v ~ | awk '{print "./bin/DataMCComparisonPlot.exe cfg/DataMCComparison_InputCfgFile/"$1}' | /bin/sh
// 2) Mu  : ls cfg/DataMCComparison_InputCfgFile/ | grep _mu_ | grep RSG | grep cfg | grep -v ~ | awk '{print "./bin/DataMCComparisonPlot.exe cfg/DataMCComparison_InputCfgFile/"$1}' | /bin/sh


int main (int argc, char **argv){

  if(argc!=2){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit "<<std::endl; return -1; }

  // Load TTree Lybrary                                                                                                                                                                      

  gSystem->Load("libTree.so");

  // Set Root style from global enviroment path                                                                                                                                               
 
  std::string ROOTStyle =  getenv ("ROOTStyle");

  gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());

  // parse config file parameter
  
  parseConfigFile(argv[1]);

  std::string InputDirectory     = gConfigParser -> readStringOption("Input::InputDirectory");
  std::string InputSampleList    = gConfigParser -> readStringOption("Input::InputSampleList");
  std::string InputVariableList  = gConfigParser -> readStringOption("Input::InputVariableList");

  std::string InputVariableListBlinded ;
  try{  InputVariableListBlinded = gConfigParser -> readStringOption("Input::InputVariableListBlinded");}
  catch(const char* exceptionString){ InputVariableListBlinded = "NULL" ;
                                      std::cerr<<" InputVariableListBlinded Set by default to --> NULL "<<std::endl;
  }                                     

  std::string InputCutList       = gConfigParser -> readStringOption("Input::InputCutList");

  std::cout<<"      "<<std::endl;
  std::cout<<" InputDirectory: "<<InputDirectory<<std::endl;
  std::cout<<"      "<<std::endl;
  std::cout<<" InputSampleList: "<<InputSampleList<<std::endl;
  std::cout<<"      "<<std::endl;
  std::cout<<" InputVariableList: "<<InputVariableList<<std::endl;
  std::cout<<"      "<<std::endl;
  std::cout<<" InputVariableListBlinded: "<<InputVariableListBlinded<<std::endl;
  std::cout<<"      "<<std::endl;
  std::cout<<" InputCutList: "<<InputCutList<<std::endl;
  std::cout<<"      "<<std::endl;


  std::string TreeName ;
  try{ TreeName  = gConfigParser -> readStringOption("Input::TreeName");}
  catch(char const* exceptionString){ TreeName = "WJet"; 
                                      std::cerr<<" TreeName Set by default to --> WJet "<<std::endl;
  }

  std::cout<<" TreeName: "<<TreeName<<std::endl;
  std::cout<<"      "<<std::endl;
  
  std::string SignalggHName ; 
  try{  SignalggHName = gConfigParser -> readStringOption("Input::SignalggHName");}
  catch(char const* exceptionString) { SignalggHName = "NULL";
                                       std::cerr<<" No Signal ggH Name --> NULL  "<<std::endl;
  }

  std::cout<<" SignalggHName: "<<SignalggHName<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string SignalqqHName ;
  try{  SignalqqHName = gConfigParser -> readStringOption("Input::SignalqqHName");}
  catch(char const* exceptionString) { SignalqqHName = "NULL";
                                       std::cerr<<" No Signal qqH Name --> NULL "<<std::endl;
  }

  std::cout<<" SignalqqHName: "<<SignalqqHName<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string SignalRSGPythiaName ;
  try{  SignalRSGPythiaName = gConfigParser -> readStringOption("Input::SignalPythiaName");}
  catch(char const* exceptionString){ SignalRSGPythiaName = "NULL";
                                      std::cerr<<" No Signal RSG Pythia Name --> NULL "<<std::endl;
  }

  std::cout<<" SignalRSGPythiaName: "<<SignalRSGPythiaName<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string SignalRSGHerwigName ;
  try{  SignalRSGHerwigName = gConfigParser -> readStringOption("Input::SignalHerwigName");}
  catch(char const* exceptionString) { SignalRSGHerwigName = "NULL";
                                       std::cerr<<" No Signal RSG Herwig Name --> NULL "<<std::endl;
  }

  std::cout<<" SignalRSGHerwigName: "<<SignalRSGHerwigName<<std::endl;
  std::cout<<"      "<<std::endl;


  bool   WithoutData  = gConfigParser -> readBoolOption("Input::WithoutData");


  std::string BackgroundWeight   = gConfigParser -> readStringOption("Option::BackgroundWeight");

  std::cout<<" BackgroundWeight: "<<BackgroundWeight<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string SignalggHWeight;
  try{ SignalggHWeight  = gConfigParser -> readStringOption("Option::SignalggHWeight");}
  catch(char const* exceptionString) { SignalggHWeight="1";
                                       std::cerr<<" Weight ggH set to --> 1 Default "<<std::endl;
  }

  std::cout<<" SignalggHWeight: "<<SignalggHWeight<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string SignalqqHWeight;
  try{ SignalqqHWeight  = gConfigParser -> readStringOption("Option::SignalqqHWeight");}
  catch(char const* exceptionString) { SignalqqHWeight="1";
                                       std::cerr<<" Weight qqH set to --> 1 Default "<<std::endl;
  }

  std::cout<<" SignalqqHWeight: "<<SignalqqHWeight<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string SignalRSGPythiaWeight;
  try{ SignalRSGPythiaWeight  = gConfigParser -> readStringOption("Option::SignalPythiaWeight");}
  catch(char const* exceptionString) { SignalRSGPythiaWeight="1";
                                       std::cerr<<" Weight Pyhtia set to --> 1 Default "<<std::endl;
  }

  std::cout<<" SignalRSGPythiaWeight: "<<SignalRSGPythiaWeight<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string SignalRSGHerwigWeight;
  try{ SignalRSGHerwigWeight  = gConfigParser -> readStringOption("Option::SignalHerwigWeight");}
  catch(char const* exceptionString) { SignalRSGHerwigWeight="1";
                                       std::cerr<<" Weight Herwig set to --> 1 Default "<<std::endl;
  }

  std::cout<<" SignalRSGHerwigWeight: "<<SignalRSGHerwigWeight<<std::endl;
  std::cout<<"      "<<std::endl;

 
  double Lumi         = gConfigParser -> readDoubleOption("Option::Lumi");

  std::cout<<" Lumi: "<<Lumi<<std::endl;
  std::cout<<"      "<<std::endl;


  double SignalScaleFactor;
  try{ SignalScaleFactor  = gConfigParser -> readDoubleOption("Option::SignalScaleFactor");}
  catch(char const* exceptionString) { SignalScaleFactor=100;
                                       std::cerr<<" Signal Scale Factor  --> 100 Default "<<std::endl;
  }

  std::cout<<" Signal Scale Factor : "<<SignalScaleFactor<<std::endl;
  std::cout<<"      "<<std::endl;


  std::string OutputRootDirectory   = gConfigParser -> readStringOption("Output::OutputRootDirectory");

  std::cout<<" OutputRootDirectory: "<<OutputRootDirectory<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string OutputRootFile        = gConfigParser -> readStringOption("Output::OutputRootFile");

  std::cout<<" OutputRootFile: "<<OutputRootFile<<std::endl;
  std::cout<<"      "<<std::endl;

  
  std::string command = "if [ ! -e "+OutputRootDirectory+" ] ; then mkdir "+OutputRootDirectory+" ; fi";
  std::cout<<" command = "<<command<<std::endl;
  std::cout<<"           "<<std::endl;

  system(command.c_str());

  
  std::string OutputPlotDirectory = OutputRootDirectory+"/"+OutputRootFile ;
  OutputPlotDirectory.replace(OutputPlotDirectory.end()-5,OutputPlotDirectory.end(),"_plot");
  std::cout<<" OutputPlotDirectory = "<<OutputPlotDirectory<<std::endl;
  std::cout<<"           "<<std::endl;
  
  
  command = "if [ ! -e "+OutputPlotDirectory+" ] ; then mkdir "+OutputPlotDirectory+" ; fi";
  std::cout<<" command = "<<command<<std::endl;
  std::cout<<"           "<<std::endl;

  system(command.c_str());

  
  command = "if [ ! -f "+OutputPlotDirectory+" ] ; then rm "+OutputPlotDirectory+"/* ; fi";
  std::cout<<" command = "<<command<<std::endl;
  std::cout<<"           "<<std::endl;
  
  system(command.c_str());
  

  TFile *outputFile = new TFile((OutputRootDirectory+"/"+OutputRootFile).c_str(),"RECREATE");


  // read sample input file list to Plot

  std::vector <std::string> NameSample;
  std::vector <std::string> NameReducedSample;
  std::vector <int> ColorSample;
  std::vector <double> SampleCrossSection;
  std::vector <int> NumEntriesBefore;


  if(ReadInputSampleFile(InputSampleList,NameSample,NameReducedSample,ColorSample,SampleCrossSection,NumEntriesBefore) <= 0){ 
    std::cerr<<" Empty Input Sample File or not Exisisting --> Exit "<<std::endl; return -1;}


  // read input variables to Plot

  std::vector <std::string> Variables;
  std::vector <double> VariablesMinValue;
  std::vector <double> VariablesMaxValue;
  std::vector <int> VariablesNbin;
  std::vector <std::string> VariablesTitle;

 
  if(ReadInputVariableFile(InputVariableList,Variables,VariablesNbin,VariablesMinValue,VariablesMaxValue,VariablesTitle) <= 0){ 
    std::cerr<<" Empty Variable List File or not Exisisting --> Exit "<<std::endl; return -1;}

  std::vector <double> VariablesBlindedMinValue(Variables.size(),-999.);
  std::vector <double> VariablesBlindedMaxValue(Variables.size(),-999.);

 
  if(InputVariableListBlinded == "NULL" || ReadInputVariableBlindedFile(InputVariableListBlinded,Variables,VariablesNbin,VariablesMinValue,VariablesMaxValue,
                                                                         VariablesBlindedMinValue, VariablesBlindedMaxValue, VariablesTitle) <=0)
    std::cerr<<" Empty Variable Blinded List File or not Exisisting --> Exit "<<std::endl;

 

  std::vector <std::string> CutList;

  if(ReadInputCutFile(InputCutList,CutList) <= 0){ 
    std::cerr<<" Empty Cut List File or not Exisisting --> Exit "<<std::endl; return -1;}

  // take input File and related tree from sampleList and do the plots

  std::vector <TTree*> TreeVect;
  std::vector <TFile*> FileVect;

  TH1F* histos[CutList.size()][Variables.size()][NameSample.size()];
  TString hname ;
   
  for (size_t iCut=0; iCut<CutList.size(); iCut++){
    
    std::cout<<std::endl;
    std::cout<<" Cut String "<<CutList.at(iCut)<<std::endl;
    std::cout<<std::endl;

    for (size_t iVar=0; iVar<Variables.size(); iVar++){

     std::cout<<std::endl;
     std::cout<<" Variable "<<Variables.at(iVar)<<std::endl;
     std::cout<<std::endl;

     for (size_t iSample=0; iSample<NameSample.size(); iSample++){

       TString NameFile = Form("%s/%s.root",InputDirectory.c_str(),NameSample.at(iSample).c_str());
       std::cout<<" Input File : "<< NameFile.Data()<<std::endl;

       FileVect.push_back ( new TFile (NameFile.Data(),"READ") );  
       TreeVect.push_back( (TTree*) FileVect.at(iSample)->Get(TreeName.c_str()));
        
       hname.Form ("%s_%s_%d",NameSample.at(iSample).c_str(),Variables.at(iVar).c_str(),int(iCut));
       histos[iCut][iVar][iSample] = new TH1F (hname.Data(),"",VariablesNbin.at(iVar),VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar));
       histos[iCut][iVar][iSample]->Sumw2();
       
       if( NameReducedSample.at(iSample) == "DATA" && !WithoutData) TreeVect.at(iSample)-> Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(), (CutList.at(iCut)).c_str() ,"goff");
	 
       else if(NameReducedSample.at(iSample) == SignalggHName && SignalggHName!="NULL")
                TreeVect.at(iSample)->Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(),("("+SignalggHWeight+")*( "+CutList.at(iCut)+")").c_str() ,"goff");
		
       else if(NameReducedSample.at(iSample) == SignalqqHName && SignalqqHName!="NULL") 
	 TreeVect.at(iSample)->Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(),("("+SignalqqHWeight+")*( "+CutList.at(iCut)+")").c_str() ,"goff");

       else if(NameReducedSample.at(iSample) == SignalRSGPythiaName && SignalRSGPythiaName!="NULL") 
	 TreeVect.at(iSample)->Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(),("("+SignalRSGPythiaWeight+")*( "+CutList.at(iCut)+")").c_str() ,"goff");

       else if(NameReducedSample.at(iSample) == SignalRSGHerwigName && SignalRSGHerwigName!="NULL") 
	 TreeVect.at(iSample)->Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(),("("+SignalRSGHerwigWeight+")*( "+CutList.at(iCut)+")").c_str() ,"goff");

       else  TreeVect.at(iSample)->Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(),("("+BackgroundWeight+") * ("+CutList.at(iCut)+")").c_str() ,"goff");

       histos[iCut][iVar][iSample]->SetFillColor(ColorSample.at(iSample));
       histos[iCut][iVar][iSample]->SetLineColor(ColorSample.at(iSample));
       histos[iCut][iVar][iSample]->GetXaxis()->SetTitle(Variables.at(iVar).c_str());
       histos[iCut][iVar][iSample]->GetXaxis()->SetTitleSize(0.02);
       histos[iCut][iVar][iSample]->GetYaxis()->SetTitle("Entries");
       histos[iCut][iVar][iSample]->GetYaxis()->SetTitleSize(0.02);

    } 
   }
  }
  
  // Normalization to the lumi of MC samples 
  double norm;

  for (size_t iCut=0; iCut<CutList.size(); iCut++){

    for (size_t iVar=0; iVar<Variables.size(); iVar++){
 
     for (size_t iSample=0; iSample<NameSample.size(); iSample++){
      
       if(NameReducedSample.at(iSample) == "DATA")  continue;
       norm =  Lumi*SampleCrossSection.at(iSample) / NumEntriesBefore.at(iSample);
       histos[iCut][iVar][iSample]->Scale(1.*norm);
     }	
   }
  }

  outputFile->cd();

  TCanvas* c[CutList.size()][Variables.size()];
  TCanvas* cLog[CutList.size()][Variables.size()];

  TLegend* leg[CutList.size()][Variables.size()];

  TH1F*    histo_top[CutList.size()][Variables.size()];
  TH1F*    histo_diboson[CutList.size()][Variables.size()];
  TH1F*    histo_WJets[CutList.size()][Variables.size()];
  TH1F*    histoSum[CutList.size()][Variables.size()];
  TH1F*    RatioDataMC[CutList.size()][Variables.size()];

  THStack* hs[CutList.size()][Variables.size()];

  int iSampleData = 0;
  int iSampleggH = -1;
  int iSampleqqH = -1;
  int iSampleRSGPythia = -1;
  int iSampleRSGHerwig = -1;
  
  for (size_t iCut=0; iCut<CutList.size(); iCut++){

     for (size_t iVar=0; iVar<Variables.size(); iVar++){

          TString CanvasName = Form("%s_%zu",Variables.at(iVar).c_str(),iCut);         
	  c[iCut][iVar] = new TCanvas (CanvasName.Data() ,"" ) ;
  
          TString CanvasNameLog = Form("%s_%zu_Log",Variables.at(iVar).c_str(),iCut);
          cLog[iCut][iVar] = new TCanvas (CanvasNameLog.Data() ,"" ) ;

	  leg[iCut][iVar] = new TLegend (0.81, 0.6, 0.99, 0.90) ;
	  leg[iCut][iVar]->SetFillColor(0);
 
          TString histoName = Form("%s_sTop_%d",Variables.at(iVar).c_str(),int(iCut));

	  histo_top[iCut][iVar] = new TH1F ( histoName.Data(),"",VariablesNbin.at(iVar),
                                             VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar)) ;

          histoName = Form("%s_diboson_%d",Variables.at(iVar).c_str(),int(iCut));

	  histo_diboson[iCut][iVar] = new TH1F ( histoName.Data(),"",VariablesNbin.at(iVar),
                                             VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar)) ;

          histoName = Form("%s_WJets_%d",Variables.at(iVar).c_str(),int(iCut));

	  histo_WJets[iCut][iVar] = new TH1F ( histoName.Data(),"",VariablesNbin.at(iVar),
                                             VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar)) ;

          histoName = Form("%s_stack_%d",Variables.at(iVar).c_str(),int(iCut));

	  hs[iCut][iVar] = new THStack (histoName.Data(),"") ;

          histoName = Form("%s_sum_%d",Variables.at(iVar).c_str(),int(iCut));

	  histoSum[iCut][iVar] = new TH1F ( histoName.Data(),"",
					    VariablesNbin.at(iVar),VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar)) ;


          TPad* upperPad ; TPad* lowerPad ; TPad* upperPadLog ; TPad* lowerPadLog ;  
	   
          if(!WithoutData){

              upperPad = new TPad("upperPad", "upperPad", .005, .250, .995, .910);
	      lowerPad = new TPad("lowerPad", "lowerPad", .005, .025, .995, .250);

              upperPadLog = new TPad("upperPadLog", "upperPadLog", .005, .250, .995, .910);
	      lowerPadLog = new TPad("lowerPadLog", "lowerPadLog", .005, .025, .995, .250);

          }
          else{
	        upperPad = new TPad("upperPad", "upperPad", .005, .050, .995, .910); 

                upperPadLog = new TPad("upperPadLog", "upperPadLog", .005, .050, .995, .910);

	  }

	  c[iCut][iVar] ->cd();
 	  if(!WithoutData) lowerPad->Draw();
	  upperPad->Draw();     

	  cLog[iCut][iVar] ->cd();
 	  if(!WithoutData) lowerPadLog->Draw();
	  upperPadLog->Draw();     
	  
	  for (size_t iSample = 0; iSample<NameSample.size(); iSample++){
	  
	    if( NameReducedSample.at(iSample) == "DATA" && !WithoutData){

	      histos[iCut][iVar][iSample]->SetLineColor(kBlack);
	      histos[iCut][iVar][iSample]->SetLineStyle(1);
              histos[iCut][iVar][iSample]->GetXaxis()->SetTitle((VariablesTitle.at(iVar)).c_str());
              histos[iCut][iVar][iSample]->GetXaxis()->SetTitleSize(0.04);
	      histos[iCut][iVar][iSample]->GetYaxis()->SetTitle("Entries");
              histos[iCut][iVar][iSample]->GetYaxis()->SetTitleSize(0.04);

	      iSampleData = iSample;                                                                       
	      leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSample], (NameReducedSample.at(iSample)).c_str(), "ple" ); 
	    }
	    
	    else if ( NameReducedSample.at(iSample)==SignalggHName && SignalggHName!="NULL"){
	      iSampleggH = iSample;
	    }

	    else if ( NameReducedSample.at(iSample)==SignalqqHName && SignalqqHName!="NULL"){
	      iSampleqqH = iSample;
	    }
            
            else if ( NameReducedSample.at(iSample)==SignalRSGPythiaName && SignalRSGPythiaName!="NULL"){
	      iSampleRSGPythia = iSample;
	    }

	     else if ( NameReducedSample.at(iSample)==SignalRSGHerwigName && SignalRSGHerwigName!="NULL"){
	      iSampleRSGHerwig = iSample;
	    }
           
	    
	    else if (( NameReducedSample.at(iSample)=="STop") || ( NameReducedSample.at(iSample)=="tt_bar") )
	    {  
		histo_top[iCut][iVar]->SetFillColor(ColorSample.at(iSample));
		histo_top[iCut][iVar]->SetLineColor(ColorSample.at(iSample));
		histo_top[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
		histoSum[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
	    }
	    else if ( NameReducedSample.at(iSample)=="W+Jets" )
	    {  
		histo_WJets[iCut][iVar]->SetFillColor(ColorSample.at(iSample));
		histo_WJets[iCut][iVar]->SetLineColor(ColorSample.at(iSample));
		histo_WJets[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
		histoSum[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
	    } 
	    else if (( NameReducedSample.at(iSample)=="WW") || ( NameReducedSample.at(iSample)=="WZ") || ( NameReducedSample.at(iSample)=="ZZ") )
	    {  
		histo_diboson[iCut][iVar]->SetFillColor(ColorSample.at(iSample));
		histo_diboson[iCut][iVar]->SetLineColor(ColorSample.at(iSample));
		histo_diboson[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
		histoSum[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
	    }
	    else
	    {
		hs[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
		leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSample], (NameReducedSample.at(iSample)).c_str(), "fl" );
		histoSum[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
	    }

	  }
	  
	  leg[iCut][iVar]->AddEntry( histo_top[iCut][iVar], "Top", "fl" );
	  leg[iCut][iVar]->AddEntry( histo_diboson[iCut][iVar], "diBoson", "fl" );
	  leg[iCut][iVar]->AddEntry( histo_WJets[iCut][iVar], "W+Jets", "fl" );

	  hs[iCut][iVar]->Add(histo_top[iCut][iVar]);
	  hs[iCut][iVar]->Add(histo_WJets[iCut][iVar]);
	  hs[iCut][iVar]->Add(histo_diboson[iCut][iVar]);
     
          upperPad->cd();


	  DrawStackError(hs[iCut][iVar],0,Variables.at(iVar));
                         
	  if(!WithoutData){ 

           TObjArray* histoList = hs[iCut][iVar] -> GetStack();
           TH1F* histo          = (TH1F*) histoList->At(histoList -> GetEntries()-1);

	   upperPad->RangeAxis(histos[iCut][iVar][iSampleData]->GetXaxis()->GetXmin(),0.,histos[iCut][iVar][iSampleData]->GetXaxis()->GetXmax(),
                               fabs(std::max(histos[iCut][iVar][iSampleData]->GetYaxis()->GetXmax(),histo->GetYaxis()->GetXmax())));
                                                               

            if((VariablesBlindedMinValue.at(iVar) != -999. && VariablesBlindedMaxValue.at(iVar) != -999.) && VariablesBlindedMinValue.at(iVar) != VariablesBlindedMaxValue.at(iVar)){
	     
                  for(int iBin = histos[iCut][iVar][iSampleData]->FindBin(VariablesBlindedMinValue.at(iVar)) ; 
                      iBin< histos[iCut][iVar][iSampleData]->FindBin(VariablesBlindedMaxValue.at(iVar)) ; iBin++) histos[iCut][iVar][iSampleData]->SetBinContent(iBin,-1.);

                  histos[iCut][iVar][iSampleData]->Draw("E same");
                
            }
            else histos[iCut][iVar][iSampleData]->Draw("E same");
        
          }
	  
	  if(SignalggHName!="NULL" && iSampleggH!=-1){ 
              
	            TString Name = Form("%s*%d",NameReducedSample.at(iSampleggH).c_str(),int(SignalScaleFactor));
                    leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSampleggH], Name.Data(), "l" );
                                     
                    histos[iCut][iVar][iSampleggH]->SetLineWidth(2);
                    histos[iCut][iVar][iSampleggH]->Scale(SignalScaleFactor*1.);
                    histos[iCut][iVar][iSampleggH]->SetFillStyle(0);

                    histos[iCut][iVar][iSampleggH]->Draw("hist same");
          }
       
          if(SignalqqHName!="NULL" && iSampleqqH!=-1){ 

          	    TString Name = Form("%s*%d",NameReducedSample.at(iSampleqqH).c_str(),int(SignalScaleFactor)); 
               	    leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSampleqqH], Name.Data(), "l" );
	            histos[iCut][iVar][iSampleqqH]->SetFillStyle(0);
                    histos[iCut][iVar][iSampleqqH]->SetLineWidth(2);
          
	            histos[iCut][iVar][iSampleqqH]->Scale(SignalScaleFactor*1.);
	            histos[iCut][iVar][iSampleqqH]->Draw("hist same");

          }

	  if(SignalRSGPythiaName!="NULL" && iSampleRSGPythia!=-1){ 
 
 	            TString Name = Form("%s*%d",NameReducedSample.at(iSampleRSGPythia).c_str(),int(SignalScaleFactor)); 
               	    leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSampleRSGPythia], Name.Data(), "l" );
                    histos[iCut][iVar][iSampleRSGPythia]->SetFillStyle(0);
                                     
                    histos[iCut][iVar][iSampleRSGPythia]->SetLineWidth(2);
                    histos[iCut][iVar][iSampleRSGPythia]->Scale(SignalScaleFactor*1.);

                    histos[iCut][iVar][iSampleRSGPythia]->Draw("hist same");
          }
       
	  if(SignalRSGHerwigName!="NULL" && iSampleRSGHerwig!=-1){ 

	            TString Name = Form("%s*%d",NameReducedSample.at(iSampleRSGHerwig).c_str(),int(SignalScaleFactor)); 
               	    leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSampleRSGHerwig], Name.Data(), "l" );
                                     
                    histos[iCut][iVar][iSampleRSGHerwig]->SetFillStyle(0);
                    histos[iCut][iVar][iSampleRSGHerwig]->SetLineWidth(2);
                    histos[iCut][iVar][iSampleRSGHerwig]->Scale(SignalScaleFactor*1.);

                    histos[iCut][iVar][iSampleRSGHerwig]->Draw("hist same");
          }
       
	  leg[iCut][iVar]->Draw("same");

          LatexCMS(Lumi);

	  if(!WithoutData) {
         
                            lowerPad->cd();
	                    lowerPad->SetGridx();
	                    lowerPad->SetGridy();
	  
                            RatioDataMC[iCut][iVar] = (TH1F*) histos[iCut][iVar][iSampleData]->Clone(("RatioDataMC-"+Variables.at(iVar)+"-"+CutList.at(iCut)).c_str()) ;
                            RatioDataMC[iCut][iVar]->Divide(histoSum[iCut][iVar]);
                            RatioDataMC[iCut][iVar]->SetMinimum(0.4);
                            RatioDataMC[iCut][iVar]->SetMaximum(1.8);
			
                            RatioDataMC[iCut][iVar]->GetXaxis()->SetLabelSize(0.09);

                            RatioDataMC[iCut][iVar]->GetYaxis()->SetLabelSize(0.09);
                            RatioDataMC[iCut][iVar]->GetYaxis()->CenterTitle();
			    RatioDataMC[iCut][iVar]->GetYaxis()->SetTitle("Data/MC");
			    RatioDataMC[iCut][iVar]->GetYaxis()->SetTitleSize(0.12);
			    RatioDataMC[iCut][iVar]->GetYaxis()->SetTitleOffset(0.4);

	                    RatioDataMC[iCut][iVar]->Draw("PE");
         
                            lowerPadLog->cd();
	                    lowerPadLog->SetGridx();
	                    lowerPadLog->SetGridy();

	                    RatioDataMC[iCut][iVar]->Draw("PE");

	  }
	  
	  c[iCut][iVar]->Write();
          
	  std::string canvasname = CanvasName.Data() ;
	  std::replace(canvasname.begin(),canvasname.end(),'[','_');
	  std::replace(canvasname.begin(),canvasname.end(),']','_');
	            
	  c[iCut][iVar]->Print( (OutputPlotDirectory+"/"+canvasname+".pdf").c_str(),"pdf");
	  c[iCut][iVar]->Print( (OutputPlotDirectory+"/"+canvasname+".png").c_str(),"png");
	  //	  c[iCut][iVar]->Print( (OutputPlotDirectory+"/"+canvasname+".eps").c_str(),"eps");

	  c[iCut][iVar]->Close();
	  
	  if(!WithoutData){  lowerPadLog->cd();
                             lowerPadLog->SetGridx();
	                     lowerPadLog->SetGridy();

                             RatioDataMC[iCut][iVar]->Draw("PE");
          }

	  upperPadLog->cd();
	  upperPadLog->SetLogy();

	  DrawStackError(hs[iCut][iVar],0,Variables.at(iVar));

	  if(!WithoutData){

           TObjArray* histoList = hs[iCut][iVar] -> GetStack();
           TH1F* histo          = (TH1F*) histoList->At(histoList -> GetEntries()-1);

	   upperPadLog->RangeAxis(histos[iCut][iVar][iSampleData]->GetXaxis()->GetXmin(),0.,histos[iCut][iVar][iSampleData]->GetXaxis()->GetXmax(),
                                  fabs(std::max(histos[iCut][iVar][iSampleData]->GetYaxis()->GetXmax(),histo->GetYaxis()->GetXmax())));
                              

           if((VariablesBlindedMinValue.at(iVar) != -999. && VariablesBlindedMaxValue.at(iVar) != -999.) && VariablesBlindedMinValue.at(iVar) != VariablesBlindedMaxValue.at(iVar)){
	     
                  for(int iBin = histos[iCut][iVar][iSampleData]->FindBin(VariablesBlindedMinValue.at(iVar)) ; 
                      iBin< histos[iCut][iVar][iSampleData]->FindBin(VariablesBlindedMaxValue.at(iVar)) ; iBin++) histos[iCut][iVar][iSampleData]->SetBinContent(iBin,-1.);

                  histos[iCut][iVar][iSampleData]->Draw("E same");
                
            }
            else histos[iCut][iVar][iSampleData]->Draw("E same");
          }

	  if(SignalggHName!="NULL" && iSampleggH!=-1) histos[iCut][iVar][iSampleggH]->Draw("hist same");
	  if(SignalqqHName!="NULL" && iSampleqqH!=-1) histos[iCut][iVar][iSampleqqH]->Draw("hist same");
	  if(SignalRSGPythiaName!="NULL" && iSampleRSGPythia!=-1) histos[iCut][iVar][iSampleRSGPythia]->Draw("hist same");
	  if(SignalRSGHerwigName!="NULL" && iSampleRSGHerwig!=-1) histos[iCut][iVar][iSampleRSGHerwig]->Draw("hist same");

	  leg[iCut][iVar]->Draw("same");
          LatexCMS(Lumi);
       
	  cLog[iCut][iVar]->Write();
          std::string canvasnameLog = CanvasNameLog.Data() ;
	  std::replace(canvasnameLog.begin(),canvasnameLog.end(),'[','_');
	  std::replace(canvasnameLog.begin(),canvasnameLog.end(),']','_');
          
	  cLog[iCut][iVar]->Print( (OutputPlotDirectory+"/"+canvasnameLog+".pdf").c_str(),"pdf");
	  cLog[iCut][iVar]->Print( (OutputPlotDirectory+"/"+canvasnameLog+".png").c_str(),"png");
	  //	  cLog[iCut][iVar]->Print( (OutputPlotDirectory+"/"+canvasnameLog+".eps").c_str(),"eps");

          cLog[iCut][iVar]->Close();
     }
    }
  
 outputFile->Close();

 return 0 ;
}



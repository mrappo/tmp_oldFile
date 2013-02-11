#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>
#include <map>

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


void LatexCMS (double lumi) ;

// To lunch the sequence : 

// 1) Electron : ls cfg/DataMCComparison_InputCfgFile/ | grep _el_ | grep cfg | grep -v ~ | awk '{print "./bin/DataMCComparisonPlot.exe cfg/DataMCComparison_InputCfgFile/"$1}' | /bin/sh
// 1) Muon     : ls cfg/DataMCComparison_InputCfgFile/ | grep _mu_ | grep cfg | grep .v ~ | awk '{print "./bin/DataMCComparisonPlot.exe cfg/DataMCComparison_InputCfgFile/"$1}' | /bin/sh

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
  std::string InputCutList       = gConfigParser -> readStringOption("Input::InputCutList");
  std::string TreeName           = gConfigParser -> readStringOption("Input::TreeName");

  std::string SignalqqHName      = gConfigParser -> readStringOption("Input::SignalqqHName");
  std::string SignalggHName      = gConfigParser -> readStringOption("Input::SignalggHName");

  std::string BackgroundWeight     = gConfigParser -> readStringOption("Input::BackgroundWeight");
  std::string SignalggHWeight      = gConfigParser -> readStringOption("Input::SignalggHWeight");

  std::string OutputRootDirectory     = gConfigParser -> readStringOption("Output::OutputRootDirectory");
  std::string OutputRootFile          = gConfigParser -> readStringOption("Output::OutputRootFile");

  double Lumi       = gConfigParser -> readDoubleOption("Input::Lumi");

  std::string command = "if [ ! -e "+OutputRootDirectory+" ] ; then mkdir "+OutputRootDirectory+" ; fi";
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
        
       hname.Form ("%s_%s_%s",NameSample.at(iSample).c_str(),Variables.at(iVar).c_str(),CutList.at(iCut).c_str() );
       histos[iCut][iVar][iSample] = new TH1F (hname.Data(),"",VariablesNbin.at(iVar),VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar));
       histos[iCut][iVar][iSample]->Sumw2();
       
       if( NameReducedSample.at(iSample) == "DATA") TreeVect.at(iSample)-> Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(), (CutList.at(iCut)).c_str() ,"goff");

       else if(NameReducedSample.at(iSample) != SignalggHName) TreeVect.at(iSample)->Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(),
                                                                                           ("("+BackgroundWeight+") * ("+CutList.at(iCut)+")").c_str() ,"goff");

       else TreeVect.at(iSample)->Draw((Variables.at(iVar)+" >> "+hname.Data()).c_str(),("("+SignalggHWeight+")*( "+CutList.at(iCut)+")").c_str() ,"goff");
       
       histos[iCut][iVar][iSample]->SetFillColor(ColorSample.at(iSample));
       histos[iCut][iVar][iSample]->SetLineColor(ColorSample.at(iSample));
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
  TLegend* leg[CutList.size()][Variables.size()];
  TH1F*    histo_top[CutList.size()][Variables.size()];
  TH1F*    histo_diboson[CutList.size()][Variables.size()];
  THStack* hs[CutList.size()][Variables.size()];
  TH1F*    histoSum[CutList.size()][Variables.size()];
  TH1F*    RatioDataMC[CutList.size()][Variables.size()];

  int iSampleData = 0;
  int iSampleggH = 0;
  int iSamplevbf = 0;

  for (size_t iCut=0; iCut<CutList.size(); iCut++){

     for (size_t iVar=0; iVar<Variables.size(); iVar++){

          TString CanvasName = Form("%s_%zu",Variables.at(iVar).c_str(),iCut);         
	  c[iCut][iVar] = new TCanvas (CanvasName.Data() ,"" ) ;
	   
	  TPad* upperPad = new TPad("upperPad", "upperPad", .005, .300, .995, .910);
	  TPad* lowerPad = new TPad("lowerPad", "lowerPad", .005, .005, .995, .295);

	  lowerPad->Draw();
	  upperPad->Draw();        

	  upperPad->cd();
 
	  leg[iCut][iVar] = new TLegend (0.81, 0.6, 0.99, 0.90) ;
	  leg[iCut][iVar]->SetFillColor(0);

	  histo_top[iCut][iVar] = new TH1F ( (Variables.at(iVar)+"sTop"+CutList.at(iCut)).c_str(),"",VariablesNbin.at(iVar),
                                             VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar)) ;

	  histo_diboson[iCut][iVar] = new TH1F ( (Variables.at(iVar)+"diboson"+CutList.at(iCut)).c_str(),"",VariablesNbin.at(iVar),
                                             VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar)) ;

	  hs[iCut][iVar] = new THStack ((Variables.at(iVar)+CutList.at(iCut)).c_str(),"") ;
	  histoSum[iCut][iVar] = new TH1F ((Variables.at(iVar)+"sum"+CutList.at(iCut)).c_str(),"",
					   VariablesNbin.at(iVar),VariablesMinValue.at(iVar),VariablesMaxValue.at(iVar)) ;

	  
	  for (size_t iSample = 0; iSample<NameSample.size(); iSample++){
	  
	    if( NameReducedSample.at(iSample) == "DATA"){

	      histos[iCut][iVar][iSample]->SetLineColor(kBlack);
	      histos[iCut][iVar][iSample]->SetLineStyle(1);
	      histos[iCut][iVar][iSample]->Draw("E");
              histos[iCut][iVar][iSample]->GetXaxis()->SetTitle((VariablesTitle.at(iVar)).c_str());
              histos[iCut][iVar][iSample]->GetXaxis()->SetTitleSize(0.06);
	      gPad->Modified();
	      iSampleData = iSample;                                                                       
	      leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSample], (NameReducedSample.at(iSample)).c_str(), "ple" ); 
	    }
	    
	    else if ( NameReducedSample.at(iSample)==SignalggHName  || NameReducedSample.at(iSample)==SignalqqHName ){

	      if(NameReducedSample.at(iSample)==SignalggHName) iSampleggH = iSample;
	      else iSamplevbf = iSample; 
	    }
	    
	    else if (( NameReducedSample.at(iSample)=="STop") || ( NameReducedSample.at(iSample)=="tt_bar") )
	      {  
		histo_top[iCut][iVar]->SetFillColor(ColorSample.at(iSample));
		histo_top[iCut][iVar]->SetLineColor(ColorSample.at(iSample));
		histo_top[iCut][iVar]->Add(histos[iCut][iVar][iSample]);
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
	  hs[iCut][iVar]->Add(histo_top[iCut][iVar]);
	  hs[iCut][iVar]->Add(histo_diboson[iCut][iVar]);
 
	  hs[iCut][iVar]->Draw("hist");
	  histos[iCut][iVar][iSampleData]->Draw("E same");
	  leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSampleggH], (NameReducedSample.at(iSampleggH)+"*10").c_str(), "l" );
	  leg[iCut][iVar]->AddEntry( histos[iCut][iVar][iSamplevbf], (NameReducedSample.at(iSamplevbf)+"*10").c_str(), "l" );

	  histos[iCut][iVar][iSampleggH]->SetLineWidth(2);
	  histos[iCut][iVar][iSamplevbf]->SetLineWidth(2);
          

	  histos[iCut][iVar][iSampleggH]->Scale(10*1.);
	  histos[iCut][iVar][iSamplevbf]->Scale(10*1.);
	  histos[iCut][iVar][iSampleggH]->SetFillStyle(1);
	  histos[iCut][iVar][iSamplevbf]->SetFillStyle(1);
	  histos[iCut][iVar][iSampleggH]->Draw("hist same");
	  histos[iCut][iVar][iSamplevbf]->Draw("hist same");

	  leg[iCut][iVar]->Draw("same");

          LatexCMS(Lumi);


	  lowerPad->cd();
	  lowerPad->SetGridx();
	  lowerPad->SetGridy();
	  
	  RatioDataMC[iCut][iVar] = (TH1F*) histos[iCut][iVar][iSampleData]->Clone(("RatioDataMC-"+Variables.at(iVar)+"-"+CutList.at(iCut)).c_str()) ;
          RatioDataMC[iCut][iVar]->Divide(histoSum[iCut][iVar]);
          RatioDataMC[iCut][iVar]->SetMinimum(0.5);
          RatioDataMC[iCut][iVar]->SetMaximum(1.5);
	  RatioDataMC[iCut][iVar]->Draw("PE");
	  RatioDataMC[iCut][iVar]->GetXaxis()->SetLabelSize(0.06);
	  RatioDataMC[iCut][iVar]->GetYaxis()->SetLabelSize(0.06);

	  c[iCut][iVar]->Write();
	  c[iCut][iVar]->Close();
	}
    }
 
 outputFile->Close();

 return 0 ;

}


void LatexCMS (double lumi){

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);

  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.80,0.962,"#sqrt{s} = 8 TeV");
  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.66,0.962,Form("#int #font[12]{L} dt = %.1f pb^{-1}", (float)lumi));

  latex.SetTextAlign(11); // align left
  //  latex.DrawLatex(0.15,0.93,"CMS,  #sqrt{s} = 7 TeV");//preliminary 2011");
  latex.DrawLatex(0.15,0.962,"CMS preliminary");

}


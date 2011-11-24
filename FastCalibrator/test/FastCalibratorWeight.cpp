
#include "FastCalibratorWeight.h"
#include <iostream>
#include <fstream>
#include <stdlib.h> 

#include "ConfigParser.h"
#include "ntpleUtils.h"

int main (int argc, char ** argv) 
{
    
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> FastCalibrator::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
    
  // Parse the config file
  parseConfigFile (argv[1]) ;

  std::string inputFile       = gConfigParser -> readStringOption("Input::inputFile");
  std::string inputTree       = gConfigParser -> readStringOption("Input::inputTree");
 
  bool isMiscalib = gConfigParser -> readBoolOption("Input::isMiscalib");
  bool isSaveEPDistribution = gConfigParser -> readBoolOption("Input::isSaveEPDistribution");

  std::string outputFile      = gConfigParser -> readStringOption("Output::outputFile");
  std::string outputFileDistributionEP       = gConfigParser -> readStringOption("Output::outputFileDistributionEP");

  int numberOfEvents       = gConfigParser -> readIntOption("Options::numberOfEvents");
  int useZ                 = gConfigParser -> readIntOption("Options::useZ");
  int useW                 = gConfigParser -> readIntOption("Options::useW");
  int splitStat            = gConfigParser -> readIntOption("Options::splitStat");
  int nLoops               = gConfigParser -> readIntOption("Options::nLoops");
  
  TChain * albero = new TChain (inputTree.c_str());
  albero -> Add(inputFile.c_str());
  
  //Use the whole sample statistics if numberOfEvents < 0
  if ( numberOfEvents < 0 ) numberOfEvents = albero->GetEntries(); 
  

  // run in normal mode: full statistics
  if ( splitStat == 0 ) {
   
    TString name ;
    TString outputTxtFile ;
    if(isMiscalib == true && useZ == 1 ) name = Form ("%s_Z_miscalib.root",outputFile.c_str());
    else{
           if(isMiscalib == false && useZ == 1)name = Form ("%s_Z.root",outputFile.c_str());
           else{
                 if(isMiscalib == true && useZ == 0) name = Form ("%s_miscalib.root",outputFile.c_str());
                 else  name = Form ("%s.root",outputFile.c_str());
               }
         }
     
    TFile *f1 = new TFile(name,"RECREATE");
    outputTxtFile = name - ".root" + ".txt";
    if(isSaveEPDistribution == true)
    {
     TFile *f2 = new TFile(outputFileDistributionEP.c_str(), "UPDATE");
     FastCalibratorWeight analyzer(albero,f2);
     analyzer.bookHistos(nLoops);
     analyzer.Loop(numberOfEvents, useZ, useW, splitStat, nLoops, isMiscalib,isSaveEPDistribution);
     analyzer.printOnTxt(outputTxtFile);
    }
    else
    {
     FastCalibratorWeight analyzer(albero);
     analyzer.bookHistos(nLoops);
     analyzer.Loop(numberOfEvents, useZ, useW, splitStat, nLoops, isMiscalib,isSaveEPDistribution);
     analyzer.saveHistos(f1);
     analyzer.printOnTxt(outputTxtFile);
    }
   
  }

  // run in even-odd mode: half statistics
  else if ( splitStat == 1 ) {
    
    // Prepare the outputs
    std::string evenFile = "Even_" + outputFile;
    std::string oddFile = "Odd_" + outputFile;
    TString name;
    TString name2;
    if(isMiscalib == true && useZ == 1 )
    { name = Form ("%s_Z_miscalib.root",evenFile.c_str());
      name2 = Form ("%s_Z_miscalib.root",oddFile.c_str());
    }

    if(isMiscalib == false && useZ == 1)
    {
      name = Form ("%s_Z.root",evenFile.c_str());
      name2 = Form ("%s_Z.root",oddFile.c_str());
    }
       
    if(isMiscalib == true && useZ == 0)
    { name = Form ("%s_miscalib.root",evenFile.c_str());
      name2 = Form ("%s_miscalib.root",oddFile.c_str());
       }
    
     if(isMiscalib == false && useZ == 0)
    {
      name = Form ("%s.root",evenFile.c_str());
      name2 = Form ("%s.root",oddFile.c_str());
    }

    TFile *f1 = new TFile(name,"RECREATE");
    TFile *f2 = new TFile(name2,"RECREATE");
     
    // Run on odd
    FastCalibratorWeight analyzer_even(albero);
    analyzer_even.bookHistos(nLoops);
    analyzer_even.Loop(numberOfEvents, useZ, useW, splitStat, nLoops,isMiscalib,isSaveEPDistribution);
    analyzer_even.saveHistos(f1);
  
    // Run on even
    FastCalibratorWeight analyzer_odd(albero);
    analyzer_odd.bookHistos(nLoops);
    analyzer_odd.Loop(numberOfEvents, useZ, useW, splitStat*(-1), nLoops,isMiscalib,isSaveEPDistribution);
    analyzer_odd.saveHistos(f2);
    
  }

  delete albero;
  return 0;
}

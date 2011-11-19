
#include "FastCalibratorEE.h"
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
    std::cerr << ">>>>> FastCalibratorEE::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
    
  // Parse the config file
  parseConfigFile (argv[1]) ;

  std::string inputFile    = gConfigParser -> readStringOption("Input::inputFile");
  std::string inputTree    = gConfigParser -> readStringOption("Input::inputTree");
  std::string outputFile   = gConfigParser -> readStringOption("Output::outputFile");
  std::string outputTxtFile   = gConfigParser -> readStringOption("Output::outputTxtFile");
   
  int numberOfEvents       = gConfigParser -> readIntOption("Options::numberOfEvents");
  int useZ                 = gConfigParser -> readIntOption("Options::useZ");
  int useW                 = gConfigParser -> readIntOption("Options::useW");
  int splitStat            = gConfigParser -> readIntOption("Options::splitStat");
  int nLoops               = gConfigParser -> readIntOption("Options::nLoops");
  
  TChain * albero = new TChain (inputTree.c_str()) ;
  albero -> Add(inputFile.c_str());
  
  //Use the whole sample statistics if numberOfEvents < 0
  if ( numberOfEvents < 0 ) numberOfEvents = albero->GetEntries (); 
  

  // run in normal mode: full statistics
  if ( splitStat == 0 ) {
    TFile *f1 = new TFile(outputFile.c_str(),"RECREATE");  
    FastCalibratorEE analyzer(albero);
    analyzer.bookHistos(nLoops);
    analyzer.Loop(numberOfEvents, useZ, useW, splitStat, nLoops);
    analyzer.saveHistos(f1);
  }
  // run in even-odd mode: halv statistics
  else if ( splitStat == 1 ) {
    
    // Prepare the outputs
    std::string evenFile = "Even_" + outputFile;
    std::string oddFile = "Odd_" + outputFile;
    TFile *f1 = new TFile(evenFile.c_str(),"RECREATE");  
    TFile *f2 = new TFile(oddFile.c_str(),"RECREATE");  
    
    // Run on odd
    FastCalibratorEE analyzer_even(albero);
    analyzer_even.bookHistos(nLoops);
    analyzer_even.Loop(numberOfEvents, useZ, useW, splitStat, nLoops);
    analyzer_even.saveHistos(f1);
  
    // Run on even
    FastCalibratorEE analyzer_odd(albero);
    analyzer_odd.bookHistos(nLoops);
    analyzer_odd.Loop(numberOfEvents, useZ, useW, splitStat*(-1), nLoops);
    analyzer_odd.saveHistos(f2);
    
  }

  delete albero;
  return 0;
}

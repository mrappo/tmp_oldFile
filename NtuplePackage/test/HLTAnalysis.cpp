///==== Compare two set of HLT triggers ====

#include "treeReader.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

int main(int argc, char** argv)	// chiede in ingresso il file di configurazione .cfg
{ 
 //Check if all nedeed arguments to parse are there
 if(argc != 2)
 {
  std::cerr << ">>>>> HLTAnalysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
  return 1;
 }
 
 // Parse the config file --> legge le info del file .cfg
 parseConfigFile (argv[1]) ;
 std::string treeName = gConfigParser -> readStringOption("Input::treeName");
 
 //---- Name of a txt file with the list of files to analyze. 
 //---- If this name is not given, a list of files (even a single file) is expected, if even this is not provided the program stops!
 //---- if both are provided, both lists are added!
 
 std::string inputFileList;
 std::vector<std::string> inputFileVector;
 
 bool InputExpected = true;
 bool InputExpectedList = true;
 bool InputExpectedVector = true;
 
 try {
  inputFileList = gConfigParser -> readStringOption("Input::inputFileList");
  InputExpected = true;
  InputExpectedList = true;
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
  InputExpected = false;
  InputExpectedList = false;
 }

 try {
  inputFileVector = gConfigParser -> readStringListOption("Input::inputFileVector");
  InputExpected = true;
  InputExpectedVector = true;
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
  InputExpected = false;
  InputExpectedVector  = false;
 }

 if (!InputExpected) {
  std::cerr << "*** no input file list" << std::endl;
  return 1;
 }
 
 
 TChain* chain = new TChain (treeName.c_str ());
 
 if (InputExpectedList) if (!FillChain (*chain, inputFileList.c_str ())) return 1 ;
 if (InputExpectedVector) if (!FillVectorChain (*chain, inputFileVector)) return 1 ;
 
 treeReader reader ( (TTree*) (chain));
 
 ///==== number of events ====
 int debug = gConfigParser -> readIntOption("Input::debug"); //---- if 1 then debug, else no
 int entryMAX = gConfigParser -> readIntOption("Input::entryMAX");
 int entryMIN = gConfigParser -> readIntOption("Input::entryMIN");
 int entryMOD = gConfigParser -> readIntOption("Input::entryMOD");
 
 if (entryMAX == -1) entryMAX = reader.GetEntries();
 else if (reader.GetEntries() < entryMAX) entryMAX = reader.GetEntries();
 
 std::cout << ">>>>> input::entryMIN  " << entryMIN  << std::endl;  
 std::cout << ">>>>> input::entryMAX  " << entryMAX  << std::endl;   
 std::cout << ">>>>> input::entryMOD  " << entryMOD  << std::endl;  
 
 ///==== HLT triggers paths ====
 std::vector<std::string> HLTVector_A;
 HLTVector_A = gConfigParser -> readStringListOption("Options::HLTVectorA");
 
 std::vector<std::string> HLTVector_B;
 HLTVector_B = gConfigParser -> readStringListOption("Options::HLTVectorB");

 ///==== Output ====
 std::string OutFileName    = gConfigParser -> readStringOption("Output::outFileName");
 std::cout << ">>>>> Output::outFileName  " << OutFileName  << std::endl;  
 
 TFile outFile(OutFileName.c_str(),"RECREATE");
 outFile.cd();


 TTree TreeHLT("TreeHLT","TreeHLT");

 int totAll_ = 0; 
 int totA_ = 0;
 int totB_ = 0;
 int totAandB_ = 0;
 int totAorB_ = 0;
 
 int All_; 
 int A_;
 int B_;
 int AandB_;
 int AorB_;
 
 std::vector<std::string> TotalHLT;
 
 
 TreeHLT.Branch("TotalHLT","std::vector<std::string>",&TotalHLT);
 TreeHLT.Branch("HLTVectorA","std::vector<std::string>",&HLTVector_A);
 TreeHLT.Branch("HLTVectorB","std::vector<std::string>",&HLTVector_B);
 TreeHLT.Branch("A",&A_,"A/I");
 TreeHLT.Branch("B",&B_,"B/I");
 TreeHLT.Branch("AandB",&AandB_,"AandB/I");
 TreeHLT.Branch("AorB",&AorB_,"AorB/I");
 TreeHLT.Branch("All",&All_,"All/I");
 
 for(int iEvent = entryMIN ; iEvent < entryMAX ; ++iEvent) {
  reader.GetEntry(iEvent);
  if((iEvent%entryMOD) == 0) std::cout << ">>>>> analysis::GetEntry " << iEvent  << ":" << reader.GetEntries() << " (" << entryMAX << ")" << std::endl;   
  
  All_ = 1;
  A_ = 0;
  B_ = 0;
  AandB_ = 0;
  AorB_ = 0;

  
  for (int iHLT = 0; iHLT < reader.GetString("HLT_Names")->size(); iHLT++){
   for (std::vector<std::string>::const_iterator iHLTA = HLTVector_A.begin(); iHLTA<HLTVector_A.end(); iHLTA++){
    if (reader.GetString("HLT_Names")->at(iHLT) == *iHLTA && reader.GetFloat("HLT_Accept")->at(iHLT) == 1) A_ = 1;
   }
   for (std::vector<std::string>::const_iterator iHLTB = HLTVector_B.begin(); iHLTB<HLTVector_B.end(); iHLTB++){
    if (reader.GetString("HLT_Names")->at(iHLT) == *iHLTB  && reader.GetFloat("HLT_Accept")->at(iHLT) == 1) B_ = 1;
   }
  }
  if (A_ == 1 && B_ == 1) { AandB_ = 1; totAandB_ ++;}
  if (A_ == 1 || B_ == 1) { AorB_ = 1;  totAorB_ ++; }
  if (A_ == 1) { totA_ ++;}
  if (B_ == 1) { totB_ ++;}
  totAll_ ++;
  
  TotalHLT = *(reader.GetString("HLT_Names"));
  TreeHLT.Fill();
 }   

 outFile.Write();
 
 
 std::cout << " *** Dump HLT names *** " << std::endl;
 for(int iEvent = entryMIN ; iEvent < entryMIN+1 ; ++iEvent) {
  reader.GetEntry(iEvent); 
  for (int iHLT = 0; iHLT < reader.GetString("HLT_Names")->size(); iHLT++){
   std::cout << " [" << iHLT << "] = " << reader.GetString("HLT_Names")->at(iHLT) << std::endl;
  }
 }
 
 std::cout << " *** Results HLT  *** " << std::endl;
 std::cout << " efficiency A = " << totA_ << " / " << totAll_ << " = " << static_cast<double>(totA_) / totAll_ << std::endl;
 std::cout << " efficiency B = " << totB_ << " / " << totAll_ << " = " << static_cast<double>(totB_) / totAll_ << std::endl;
 
 return 0;
 
 }

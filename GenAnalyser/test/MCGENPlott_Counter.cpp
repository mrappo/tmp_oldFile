///=========== Lunch in automatic
///== ls test/Latinos/HiggsMasses_GenLevel/*.cfg |   tr  "_" " " | tr "/" " " | awk '{print "./bin/MCGENPlott_Counter.exe test/Latinos/HiggsMasses_GenLevel/"$5" \n"}' | /bin/sh


///=== Command to generate Final table
///=== paste output/Final_Result.txt output/Cross_Section.txt >CMS.txt
///=== cat CMS.txt | awk '{print "|  "  $4*$10 "  |"}' > tmp.txt
///=== paste CMS.txt tmp.txt


#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TPie.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "qqHWWlnulnuUtils.h"
#include "Variables_Gen.h"

#include "../test/TDRStyle.cc"
#include "../test/Read.cc"
// #include "../test/DrawTools.h"

#include "PUclass.h"

#include <iomanip>



#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif


///
///
///      ___|                       |                     \  |   ___|         /      __ \     \  __ __|   \    
///     |       _ \   |   |  __ \   __|   _ \   __|      |\/ |  |            /       |   |   _ \    |    _ \   
///     |      (   |  |   |  |   |  |     __/  |         |   |  |           /        |   |  ___ \   |   ___ \  
///    \____| \___/  \__,_| _|  _| \__| \___| _|        _|  _| \____|     _/        ____/ _/    _\ _| _/    _\ 
///   
///   

// std::vector<double> PUWeight;

int GetNumList(std::vector<int> &list){
 int result = 0;
 for (int it = 0; it<list.size(); it++) result += list.at(it);
 return result;
}


int main(int argc, char** argv)
{ 
 TDRStyle();
 
 gStyle->SetPadTopMargin(0.11);
 gStyle->SetPadLeftMargin(0.07);
 gStyle->SetPadRightMargin(0.23);
 gStyle->cd(); 
 
 
 std::cout << " " << std::endl;
 std::cout << " " << std::endl;
 std::cout << "      ___|                       |                     \\  |   ___|        " << std::endl;
 std::cout << "     |       _ \\   |   |  __ \\   __|   _ \\   __|      |\\/ |  |         " << std::endl;
 std::cout << "     |      (   |  |   |  |   |  |     __/  |         |   |  |             " << std::endl;
 std::cout << "    \\____| \\___/  \\__,_| _|  _| \\__| \\___| _|        _|  _| \\____|   " << std::endl; 
 std::cout << " " << std::endl;
 std::cout << " " << std::endl; 

 char normal[] = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };
 char black[] = { 0x1b, '[', '0', ';', '3', '0', 'm', 0 };
 char red[] = { 0x1b, '[', '0', ';', '3', '1', 'm', 0 };
 char green[] = { 0x1b, '[', '0', ';', '3', '2', 'm', 0 };
 char yellow[] = { 0x1b, '[', '0', ';', '3', '3', 'm', 0 };
 char blue[] = { 0x1b, '[', '0', ';', '3', '4', 'm', 0 };
 char purple[] = { 0x1b, '[', '0', ';', '3', '5', 'm', 0 };
 char cyan[] = { 0x1b, '[', '0', ';', '3', '6', 'm', 0 };
 char Lgray[] = { 0x1b, '[', '0', ';', '3', '7', 'm', 0 };
 char Dgray[] = { 0x1b, '[', '0', ';', '3', '8', 'm', 0 };
 char Bred[] = { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
 
 ///==== for bold colors, just change the 0 after the [ to a 1
 
 EColor vColor[1000] = {
  (EColor) (kRed+1),
  (EColor) (kRed+3),
  (EColor) (kGray+1),
  (EColor) (kAzure-2),
  (EColor) (kAzure-9),
  (EColor) (kYellow),
  (EColor) (kGreen+2),
//   
  kGreen,
  //kMagenta,(EColor) (kMagenta+1),(EColor) (kMagenta+2),
  kTeal,//(EColor) (kTeal+1),
  kRed,
  kGray,
  kOrange,(EColor) (kOrange+1),
  kBlue,//(EColor)(kBlue+1),(EColor) (kBlue+2),
  (EColor) (kPink+2),//(EColor) (kPink+1),(EColor) (kPink+2),
  kViolet,
  kYellow,
  kGray,(EColor) (kGray+1),(EColor) (kViolet),(EColor) (kYellow),(EColor) (kGray)
 };
 
 
 
 
 ///==== Check if all nedeed arguments to parse are there                                                                                                                               
 if(argc != 2)
 {
  std::cerr << ">>>>> analysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
  return 1;
 }


 // Parse the config file                                                                                                                                                          
 parseConfigFile (argv[1]) ;
 
 std::string treeName  = gConfigParser -> readStringOption("Input::treeName");
 std::string treeNameSelections = gConfigParser -> readStringOption("Input::treeNameSelections");
 std::string fileSamples = gConfigParser -> readStringOption("Input::fileSamples");
 std::string inputDirectory = gConfigParser -> readStringOption("Input::inputDirectory");
 
 std::string inputBeginningFile = "out_NtupleProducer_"; 
 try {
  inputBeginningFile = gConfigParser -> readStringOption("Input::inputBeginningFile");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
 }
 std::cout << ">>>>> Input::inputBeginningFile  " << inputBeginningFile  << std::endl;  
 
 double LUMI = gConfigParser -> readDoubleOption("Input::Lumi");
  

 ///=== acuisition from sample file ---> parsing through the cfg
 
 char nameFileIn[1000];
 char *nameSample[1000];
 char *nameHumanReadable[1000];
 char* xsectionName[1000];
 

 sprintf(nameFileIn,"%s",fileSamples.c_str());

 int numberOfSamples = ReadFile(nameFileIn, nameSample, nameHumanReadable, xsectionName);

 
 ///==== Cut File read from cfg
 
 std::string CutFile = gConfigParser -> readStringOption("Selections::CutFile");
 std::vector<std::string> vCut;
 
 std::cout << " nCuts   = " << ReadFileCut(CutFile, vCut) << std::endl;
 
 
 ///==== output file ====

 std::string OutFileName    = gConfigParser -> readStringOption("Output::outFileName");
 std::cout << ">>>>> Output::outFileName  " << OutFileName  << std::endl;  
 
 TFile outFile(OutFileName.c_str(),"RECREATE");
 outFile.cd();
 
 
 
 ///==== pT Higgs reweight (begin) ====
 std::string nameptHWeight; 
 try {
  nameptHWeight = gConfigParser -> readStringOption("Input::nameptHWeight");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
 }
 std::cout << ">>>>> input::nameptHWeight  " << nameptHWeight  << std::endl;  
 if (nameptHWeight != ""){
  TString toLoad;
  toLoad = Form(".L %s+",nameptHWeight.c_str());
  gROOT->ProcessLine(toLoad.Data()); ///==== -----> Open root and execute the command toLoad.Data() on the root command line 
 }
 
 std::string nameptHWeightSample; 
 try {
  nameptHWeightSample = gConfigParser -> readStringOption("Input::nameptHWeightSample");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
 }
 std::cout << ">>>>> input::nameptHWeightSample  " << nameptHWeightSample  << std::endl;  
 
 
 ///==== pT Higgs reweight (end) ====
 
 
 ///==== debug flag ====
 
 bool  debug = false; 
 try {
  debug = gConfigParser -> readBoolOption("Input::debug");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
 }
 std::cout << ">>>>> input::debug  " << debug  << std::endl;
 
 
 
 ///==== list of variables to plot ====
 std::vector<double> vMin;
 std::vector<double> vMax;
 std::vector<int> vNBin;
 std::vector<std::string> vVarName;
 std::vector<std::string> vVarNameHR;
 std::string VarFile = gConfigParser -> readStringOption("Plot::VarFile");
 
 int numVar = ReadFileVar(VarFile,vMin,vMax,vNBin,vVarName,vVarNameHR);
  
 
 
 ///==== program ====
 
 
 double start, end;
 start = clock();
 
 int Entries[100];
 
 // Tree in order to take the infos from input file
 TTree *treeEffVect[100];
 TTree *treeJetLepVect[100];
 double Normalization[1000];
 double xsection[1000];
 double xsectionErrorUp[1000];
 double xsectionErrorDown[1000];
 
 
 ///=== save information of weight coefficient for MC related to the cross section of the process
 
 for (int iSample=0; iSample<numberOfSamples; iSample++){
  xsection[iSample] = atof(xsectionName[iSample]);
 }

for (int iSample=0; iSample<numberOfSamples; iSample++){
  
  char nameFile[20000];
  sprintf(nameFile,"%s/%s%s.root",inputDirectory.c_str(),inputBeginningFile.c_str(),nameSample[iSample]);  
  if (debug) std::cout << " nameFile = " << nameFile << std::endl;
  
  ///==== Open input file for each sample
  
  TFile* f = new TFile(nameFile, "READ");
  
  treeEffVect[iSample] = (TTree*) f->Get(treeNameSelections.c_str());
  if (treeEffVect[iSample] != 0) {
   char nameTreeEff[100];
   sprintf(nameTreeEff,"treeEff_%d",iSample); 
   treeEffVect[iSample]->SetName(nameTreeEff);      
  }
  
   treeJetLepVect[iSample] = (TTree*) f->Get(treeName.c_str());
  if (treeJetLepVect[iSample] != 0) {
   
   char nameTreeJetLep[100];
   sprintf(nameTreeJetLep,"treeJetLep_%d",iSample); 
   treeJetLepVect[iSample]->SetName(nameTreeJetLep);
  }
 }
 
 ///===== save info of human readble name ====
 
 std::vector<int> join_samples;
 std::vector<std::string> name_samples;
 
 for (int iSample=0; iSample<numberOfSamples; iSample++){
  name_samples.push_back(nameHumanReadable[iSample]);
  join_samples.push_back(-1);
 }
 
 std::vector<std::string> reduced_name_samples;
 std::vector<int>         reduced_name_samples_flag;
 for (int iSample = (numberOfSamples-1); iSample>= 0; iSample--){
  bool flag_name = false;
  for (unsigned int iName=0; iName<reduced_name_samples.size(); iName++){
   if (reduced_name_samples.at(iName) == name_samples.at(iSample)) flag_name = true;
  }
  if (flag_name == false) {
   reduced_name_samples.push_back(name_samples.at(iSample));
   reduced_name_samples_flag.push_back(-1);
  }

 }  ///=== Avoid double counting samples
 
 std::cout << " numberOfSamples = " << numberOfSamples << std::endl;
 
 
 for (int iSample = (numberOfSamples-1); iSample>= 0; iSample--){
  double XSection;
  double XSectionErrorUp;
  double XSectionErrorDown;
  int numEntriesBefore;
  double preselection_efficiency;
  if (treeEffVect[iSample] != 0) {   
   treeEffVect[iSample]->SetBranchAddress("XSection",&XSection); // indentify the branch of the tree with a plain variable
   treeEffVect[iSample]->SetBranchAddress("XSectionErrorUp",&XSectionErrorUp); 
   treeEffVect[iSample]->SetBranchAddress("XSectionErrorDown",&XSectionErrorDown); 
   treeEffVect[iSample]->SetBranchAddress("numEntriesBefore",&numEntriesBefore);
   treeEffVect[iSample]->SetBranchAddress("preselection_efficiency",&preselection_efficiency);  
   treeEffVect[iSample]->GetEntry(0);
  }
  if(XSection == 1)
  {
    std::cout << " Xsection = " << XSection << " ~~~> " << xsection[iSample] << std::endl;
    XSection = xsection[iSample];
    
  }
  
  if (numEntriesBefore != 0) {
   Normalization[iSample] = LUMI * XSection * preselection_efficiency / numEntriesBefore;
  }
  else {
   Normalization[iSample] = 0; 
  }    
  
  
 }

 
 
 TLegend* leg = new TLegend(0.8,0.25,0.98,0.78);
 bool LegendBuilt = false;

 TString lumiName = Form("#splitline{L = %.1f pb^{-1}}{#splitline{#sqrt{s} = 7}{CMS preliminary}}", LUMI);

 TLatex *latex = new TLatex(0.80, 0.90, lumiName); 
 latex->SetTextAlign(12);
 latex->SetNDC();
 latex->SetTextFont(42);
 latex->SetTextSize(0.03);
  
 if (debug) std::cout << " Cut size = " << vCut.size() << " ~~ " << std::endl;
 std::cout.precision (5) ;
 std::cout.unsetf(std::ios::scientific);

 ///========================================================================================================
 ///==================================== Plot of the different Variables ===================================
 ///========================================================================================================
 
 TH1F* histo_temp_Plot[100][103][43];
 TH1F* histo_temp_Normalized_Plot[100][103][43];
 
/* 
 ///==== cicle on selections ====
 for (unsigned int iCut = 0; iCut<vCut.size(); iCut++){
  TString Cut = Form ("%s",vCut.at(iCut).c_str());
  if (debug) std::cout << " Cut[" << iCut << ":" << vCut.size() << "] = " << Cut.Data() << " ~~ " << std::endl;
  ///==== cicle on variables to plot ====
  for (unsigned int iVar = 0; iVar<vVarName.size(); iVar++){
   if (debug) std::cout << " Var[" << iVar << ":" << vVarName.size() << "] = " << vVarName.at(iVar) << " ~~ " << std::endl;
   ///==== initialize ====
   for (unsigned int iName=0; iName<reduced_name_samples.size(); iName++){
    reduced_name_samples_flag.at(iName) = -1;
   }
   
   ///==== cicle on samples ====
   for (int iSample = (numberOfSamples-1); iSample>= 0; iSample--){
    if (debug) std::cout << " Sample[" << iSample << ":" << numberOfSamples << "] = " << nameSample[iSample] << " ~~ " << std::endl;
    TString name_histo_temp = Form("%s_%d_%d_temp_Plot",nameSample[iSample], iCut, iVar);
    histo_temp_Plot[iSample][iCut][iVar] = new TH1F(name_histo_temp,name_histo_temp,vNBin.at(iVar),vMin.at(iVar), vMax.at(iVar));
    char toDraw[1000];
    sprintf(toDraw,"%s >> %s",vVarName.at(iVar).c_str(),name_histo_temp.Data());      

    histo_temp_Plot[iSample][iCut][iVar] -> Sumw2(); //---- così mette l'errore giusto!
    
    TString CutExtended;
   
    if (nameptHWeight != "" && name_samples.at(iSample) == nameptHWeightSample){
      CutExtended = Form ("(%s) * ptHWeight(ptH)",Cut.Data());    
      }
    else {
     CutExtended = Form ("(%s)",Cut.Data());    
    }
    
    treeJetLepVect[iSample]->Draw(toDraw,CutExtended,"");
    
    if (Normalization[iSample]>0) { 
    
    TString Name_Normalized = Form("%s_Normalized",name_histo_temp.Data());
    histo_temp_Normalized_Plot[iSample][iCut][iVar]= (TH1F*) histo_temp_Plot[iSample][iCut][iVar]->Clone(Name_Normalized);
    histo_temp_Normalized_Plot[iSample][iCut][iVar] -> Scale(Normalization[iSample]); 
   }
    
   
//     std::cout << "Processing: " << blue << (((double) numberOfSamples - iSample)/numberOfSamples) << "% \r"  << normal << std::flush;
   } ///==== end cicle on samples ====
//    std::cout << "###";
  std::cout << "Processing: " << blue << (((double) iCut)/vCut.size())*100. << "% "  << normal <<  " -- " <<  blue << (((double) iVar)/vVarName.size())*100. << "% \r"  << normal << std::flush;   
  } ///==== end cicle on variables to plot ====
  //   std::cout << "***";
 } ///==== end cicle on selections ====
 std::cout << std::endl; 
 
 
 */
 
 ///========================================================================================================
 ///==================================== Efficiency Calculation ============================================
 ///========================================================================================================
 
 TH1F* histo_temp_Normalized[100][20];
 TH1F* histo_temp[100][20];
 TH1F* histo[100][20];
 TH1F* histo_Normalized[100][20];
 double XSection_reduced[100][20];
 double XSectionErrorUp_reduced[100][20];
 double XSectionErrorDown_reduced[100][20];
 
 long int numEntriesBefore_reduced[100][20];
 
 ///==== cicle on selections ====
 for (unsigned int iCut =0; iCut<vCut.size(); iCut++){
  TString Cut = Form ("%s",vCut.at(iCut).c_str());
  
  if (debug) std::cout << " Cut[" << iCut << ":" << vCut.size() << "] = " << Cut.Data() << " ~~ " << std::endl;
    for (unsigned int iName=0; iName<reduced_name_samples.size(); iName++){
   reduced_name_samples_flag.at(iName) = -1;
  }

  ///==== cicle on samples ====
  for (int iSample = (numberOfSamples-1); iSample>= 0; iSample--){
  
    if (debug) std::cout << " Sample[" << iSample << ":" << numberOfSamples << "] = " << nameSample[iSample] << " ~~ " << std::endl;
   
    TString name_histo_temp = Form("%s_%d_temp",nameSample[iSample], iCut);

    histo_temp[iSample][iCut] = new TH1F(name_histo_temp,name_histo_temp,100,-10,10000000000);
    char toDraw[1000];
    sprintf(toDraw,"1>>%s",name_histo_temp.Data());      
    
    histo_temp[iSample][iCut] -> Sumw2(); //---- così mette l'errore giusto!
   
    TString CutExtended;
    if (nameptHWeight != "" && name_samples.at(iSample) == nameptHWeightSample){  ///=== pt Higgs reweight option 
        
        CutExtended = Form ("(%s) * ptHWeight(H_pT)",Cut.Data());

    }
   else {
         CutExtended = Form ("%s",Cut.Data());  
	
   }
    
   treeJetLepVect[iSample]->Draw(toDraw,CutExtended); /// selections applied through Draw method of TTree
   
   if (Normalization[iSample]>0) { 
    
    TString Name_Normalized = Form("%s_Normalized",name_histo_temp.Data());
    histo_temp_Normalized[iSample][iCut]= (TH1F*) histo_temp[iSample][iCut]->Clone(Name_Normalized);
    histo_temp_Normalized[iSample][iCut] -> Scale(Normalization[iSample]); 
   }
   
   
  double XSection;
  double XSectionErrorUp;
  double XSectionErrorDown;
  
  int numEntriesBefore;
  double preselection_efficiency;
  
  if (treeEffVect[iSample] != 0) {   
   treeEffVect[iSample]->SetBranchAddress("XSection",&XSection); // indentify the branch of the tree with a plain variable
   treeEffVect[iSample]->SetBranchAddress("numEntriesBefore",&numEntriesBefore);
   treeEffVect[iSample]->SetBranchAddress("preselection_efficiency",&preselection_efficiency);  
   treeEffVect[iSample]->SetBranchAddress("XSectionErrorUp",&XSectionErrorUp); 
   treeEffVect[iSample]->SetBranchAddress("XSectionErrorDown",&XSectionErrorDown); 
  
   treeEffVect[iSample]->GetEntry(0);  
  }
  
  if(XSection == 1)
   {
     XSection= xsection[iSample];  
   }
  
    
   for (unsigned int iName=0; iName<reduced_name_samples.size(); iName++){
  
    if (name_samples.at(iSample) == reduced_name_samples.at(iName)){ 
      if (reduced_name_samples_flag.at(iName) == -1){
  
      TString name_histoTot_temp = Form("%s_%d_Tot_temp",reduced_name_samples.at(iName).c_str(),iCut);
      TString name_HR_histoTot_temp = Form("cut %d",iCut);
      histo[iName][iCut] = new TH1F(name_histoTot_temp,name_HR_histoTot_temp,100,-10,10000000000);
      histo[iName][iCut] -> Sumw2(); //---- così mette l'errore giusto!
      
      name_histoTot_temp = Form("%s_%d_Tot_temp_Normalized",reduced_name_samples.at(iName).c_str(),iCut);
      name_HR_histoTot_temp = Form("cut %d",iCut);
      histo_Normalized[iName][iCut] = new TH1F(name_histoTot_temp,name_HR_histoTot_temp,100,-10,10000000000);
      histo_Normalized[iName][iCut] -> Sumw2(); //---- così mette l'errore giusto!
      reduced_name_samples_flag.at(iName) = 1;
      
      XSection_reduced[iName][iCut] = XSection ; 
      XSectionErrorUp_reduced[iName][iCut] = XSectionErrorUp;
      XSectionErrorDown_reduced[iName][iCut] = XSectionErrorDown;
      
     }
     histo[iName][iCut] -> Add(histo_temp[iSample][iCut]);
     numEntriesBefore_reduced[iName][iCut] = numEntriesBefore_reduced[iName][iCut]+numEntriesBefore;
     histo_Normalized[iName][iCut] -> Add(histo_temp_Normalized[iSample][iCut]);
    
    }
   }
   std::cout <<"Processing: " << blue << (((double) iCut)/vCut.size())*100. << "% "  << normal <<  " -- " <<  red << (((double) numberOfSamples - iSample)/(numberOfSamples))*100. << "% \r"  << normal << std::flush;   
  } ///==== end cicle on samples ====
  } ///==== end cicle on selections ====
 
 if (debug) std::cout << " >>> Reprocessing ... " << std::endl;
 
 std::cout << std::endl;
 
 ///========== Efficiency at each step of the analysis
 double numEvents[100][20];
 double numEvents_Normalized[100][20];
 
 TH1F* hTrend[100];
 TH1F* hTrend_Normalized[100];

 
 std::cout<<" ############################################################################# "<<std::endl;
 std::cout<<"                        Monte Carlo Events                                     "<<std::endl;
 std::cout<<" ############################################################################# "<<std::endl;

 
 for (unsigned int iName=0; iName<reduced_name_samples.size(); iName++){
 
  TString nameTHTrend = Form("%s_Trend",reduced_name_samples.at(iName).c_str());
  hTrend[iName] = new TH1F (nameTHTrend,nameTHTrend,vCut.size()+1,0,vCut.size()+1);
  hTrend[iName]->GetXaxis()->SetTitle("Selections");
  
  hTrend[iName]->SetMarkerColor(vColor[iName]);
  hTrend[iName]->SetLineColor(vColor[iName]);
  hTrend[iName]->SetFillColor(vColor[iName]);
  hTrend[iName]->SetLineWidth(2);
  hTrend[iName]->SetFillStyle(3001);
  int numEntriesBefore;
  
  for (unsigned int iCut = 0; iCut<vCut.size(); iCut++){
  
   double error = 0;
   numEvents[iName][iCut] = histo[iName][iCut]->IntegralAndError(0,histo[iName][iCut]->GetNbinsX()+1,error);
   hTrend[iName]->SetBinContent(iCut+1,numEvents[iName][iCut]);
   hTrend[iName]->SetBinError(iCut+1,error);

   TString nameBin = Form("%d",iCut);
   hTrend[iName]->GetXaxis()->SetBinLabel(iCut+1,nameBin);
   
   if(treeEffVect[iName]!=0)
   {
    treeEffVect[iName]->SetBranchAddress("numEntriesBefore",&numEntriesBefore);
    treeEffVect[iName]->GetEntry(0);
   }
 
   
   std::cout << ">>>  numEvents[" << iName << "," << reduced_name_samples.at(iName) << "][" << iCut << "] = " << numEvents[iName][iCut] <<
   " , " << histo[iName][iCut]->GetEntries() << " , " << histo[iName][iCut]->GetEffectiveEntries()<<" , "
   <<" error: "<<sqrt(numEvents[iName][iCut])<<std::endl;
   std::cout<<">>> Total Efficiency["<<iName << "," <<reduced_name_samples.at(iName) 
   << "][" << iCut << "] = " << numEvents[iName][iCut]/(numEntriesBefore_reduced[iName][iCut]) << " , " 
   << histo[iName][iCut]->GetEntries()/(numEntriesBefore_reduced[iName][iCut]) << " , "
   << histo[iName][iCut]->GetEffectiveEntries()/(numEntriesBefore_reduced[iName][iCut]) 
   << " error: " << sqrt(numEvents[iName][iCut]*numEntriesBefore_reduced[iName][iCut]*numEntriesBefore_reduced[iName][iCut]+
   numEvents[iName][iCut]*numEvents[iName][iCut]*numEntriesBefore_reduced[iName][iCut])*
   (1/(numEntriesBefore_reduced[iName][iCut]*numEntriesBefore_reduced[iName][iCut])) << std::endl;
   std::cout<< ">>>  Cross Section [" << iName << "," << reduced_name_samples.at(iName) << "][" << iCut << "] = " <<XSection_reduced[iName][iCut] <<
               " Error Up " << XSectionErrorUp_reduced[iName][iCut] << " Error Down " <<XSectionErrorDown_reduced[iName][iCut]<<std::endl;
    std::cout <<  " Efficiency * Cross Section = " << XSection_reduced[iName][iCut]*numEvents[iName][iCut]/(numEntriesBefore_reduced[iName][iCut])<<std::endl;
       //    if(iCut!=0) std::cout<<  ">>> CutEfficiency[" << iName << "," <<  reduced_name_samples.at(iName) 
//     << "][" << iCut << "] = " << numEvents[iName][iCut]/(numEvents[iName][iCut-1]) << " , " 
//     << histo[iName][iCut]->GetEntries()/(histo[iName][iCut-1]->GetEntries()) << " , "
//     << histo[iName][iCut]->GetEffectiveEntries()/(histo[iName][iCut-1]->GetEffectiveEntries()) 
//     << " error: " << sqrt(numEvents[iName][iCut])/(numEvents[iName][iCut-1]) << std::endl; 
       
       
    }
  
 }
 
 if(LUMI!=1)
 {  
 
  std::cout<<" ############################################################################# "<<std::endl;
  std::cout<<"                        Monte Carlo Events Scaled to Lumi                      "<<std::endl;
  std::cout<<" ############################################################################# "<<std::endl;
   
 
 
  for (unsigned int iName=0; iName<reduced_name_samples.size(); iName++){
 
   TString nameTHTrend = Form("%s_Trend_Normalized",reduced_name_samples.at(iName).c_str());
   hTrend_Normalized[iName] = new TH1F (nameTHTrend,nameTHTrend,vCut.size()+1,0,vCut.size()+1);
   hTrend_Normalized[iName]->GetXaxis()->SetTitle("Selections");

   hTrend_Normalized[iName]->SetMarkerColor(vColor[iName]);
   hTrend_Normalized[iName]->SetLineColor(vColor[iName]);
   hTrend_Normalized[iName]->SetFillColor(vColor[iName]);
   hTrend_Normalized[iName]->SetLineWidth(2);
   hTrend_Normalized[iName]->SetFillStyle(3001);
   int numEntriesBefore;
  
   for (unsigned int iCut = 0; iCut<vCut.size(); iCut++){
  
   double error = 0;
 
   numEvents_Normalized[iName][iCut] = histo_Normalized[iName][iCut]->IntegralAndError(0,histo_Normalized[iName][iCut]->GetNbinsX()+1,error);
   hTrend_Normalized[iName]->SetBinContent(iCut+1,numEvents_Normalized[iName][iCut]);
   hTrend_Normalized[iName]->SetBinError(iCut+1,error);

   TString nameBin = Form("%d",iCut);
   hTrend_Normalized[iName]->GetXaxis()->SetBinLabel(iCut+1,nameBin);
   if(treeEffVect[iName]!=0)
   { treeEffVect[iName]->SetBranchAddress("numEntriesBefore",&numEntriesBefore);
     treeEffVect[iName]->GetEntry(0);
   }
    std::cout << ">>>  numEvents[" << iName << "," << reduced_name_samples.at(iName) << "][" << iCut << "] = " << numEvents_Normalized[iName][iCut] <<
   " , " << histo_Normalized[iName][iCut]->GetEntries() << " , " << histo_Normalized[iName][iCut]->GetEffectiveEntries()<<" , "
   <<" error: "<<sqrt(numEvents_Normalized[iName][iCut])<<std::endl;
     
    if(iCut!=0) std::cout<<  ">>> CutEfficiency[" << iName << "," <<  reduced_name_samples.at(iName) 
     << "][" << iCut << "] = " << numEvents_Normalized[iName][iCut]/(numEvents_Normalized[iName][iCut-1]) << " , " 
     << histo_Normalized[iName][iCut]->GetEntries()/(histo_Normalized[iName][iCut-1]->GetEntries()) << " , "
     << histo_Normalized[iName][iCut]->GetEffectiveEntries()/(histo_Normalized[iName][iCut-1]->GetEffectiveEntries()) 
     << " error: " << sqrt(numEvents_Normalized[iName][iCut])/(numEvents_Normalized[iName][iCut-1]) << std::endl; 
   
   
  }
  
 }
 }
 
 std::cout << " ****************************************************************************** " << std::endl;
 std::cout << " *************************** Table Like Output ******************************** " << std::endl;
 std::cout << " ****************************************************************************** " << std::endl;
 
 std::string nameOutDataCard = "Table_Result.txt";
 
 std::string mass = "160";
 try {
  mass = gConfigParser -> readStringOption("Input::mass");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
 }
 std::cout << ">>>>> input::mass  " << mass << std::endl;  
 
 
 ///==== output - txt file name ====
 try {
  nameOutDataCard = gConfigParser -> readStringOption("Output::TableResult");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
 }
 
 std::ofstream myfile;
 myfile.open (nameOutDataCard.c_str(),std::ios::app);
 std::cout << "Writing to: " << nameOutDataCard << std::endl;
 std::cout << std::endl;

 
 for (unsigned int iName=0; iName<reduced_name_samples.size(); iName++){
    
 for (unsigned int iCut = 0; iCut<vCut.size(); iCut++){
   
  myfile<<"  |  "<< mass <<
  "  |  " << 
   numEvents[iName][iCut]/(numEntriesBefore_reduced[iName][iCut]) <<"  |  " <<std::setprecision(8) << 
   (1./(numEntriesBefore_reduced[iName][iCut]*numEntriesBefore_reduced[iName][iCut]))*
   sqrt(numEvents[iName][iCut]*numEntriesBefore_reduced[iName][iCut]*numEntriesBefore_reduced[iName][iCut]+numEvents[iName][iCut]*numEvents[iName][iCut]
   *numEntriesBefore_reduced[iName][iCut]) << "  |  "<<
   XSection_reduced[iName][iCut]<<"  |  "<<
   XSectionErrorUp_reduced[iName][iCut] <<"  |  "<<
   XSectionErrorDown_reduced[iName][iCut] <<"  |  "<<
   XSection_reduced[iName][iCut] * numEvents[iName][iCut]/(numEntriesBefore_reduced[iName][iCut])
   << "  |  "<<std::endl;
   }
 
   
}
  

    
 std::cerr << " ******************************************* end *******************************************" << std::endl;
 end = clock();
 std::cout <<"Time = " <<  ((double) (end - start)) << " (a.u.)" << std::endl;  
 
 
 
 ///==== save output ====
 outFile.cd();
 
 
}






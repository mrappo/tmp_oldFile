// global TMVA style settings                                                                                                                                                             
#ifndef TMVA_TMVAGLOB
#define TMVA_TMVAGLOB

#include <iostream>
#include <vector>
#include <istream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSystem.h"
#include "TImage.h"
#include "TKey.h"
#include "TH1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TClass.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"

#include "RVersion.h"

// default definition for ROC curve color, style and width                                                                                                                                

static double color []      = {1,2,210,6,7,4,12,95};
static double linestyle []  = {1,1,9,7,1,9,1,9};
static double linewidth []  = {2.6,2.6,2.5,2,2,2,2,2};

static std::vector<double> vec_color(color, color + sizeof(color)/sizeof(double));
static std::vector<double> vec_linewidth(linewidth, linewidth + sizeof(linewidth)/sizeof(double));
static std::vector<double> vec_linestyle(linestyle, linestyle + sizeof(linestyle)/sizeof(double));

class TMVAGlob {

  public:  

  TMVAGlob();
  ~TMVAGlob();
 
  enum TypeOfPlot { kId = 0,
		    kNorm,
		    kDecorrelated,
		    kPCA,
		    kGaussDecorr,
		    kNumOfMethods };


  void openFileInput( const TString& fin );
  void openFileInput( const std::vector<std::string> & fin );

  TFile* GetInputFile (const TString& fin );
  std::vector<TFile*> GetInputFile();

  void CreateCanvasandFrame(TFile* inputFile, const double & ptMin, const double & ptMax, const std::string & outputPlotDirectory);

  TKey *NextKey( TIter & keyIter, TString className);
  int GetListOfKeys( TList& keys, TString inherits, TDirectory *dir = 0 );

  void DestroyCanvases(); 
  void Initialize( Bool_t useTMVAStyle = kTRUE );
  void PrintImage( TDirectory* dir = 0 , const std::string & outputPlotDirectory = "");    
  void NormalizeHist( TH1* h );
  void NormalizeHists( TH1* sig, TH1* bkg = 0 );
  void GetMethodName( TString & name, TKey * mkey );
  void GetMethodTitle( TString & name, TKey * ikey );
  void GetMethodName( TString & name, TDirectory * mdir );
  void GetMethodTitle( TString & name, TDirectory * idir );

  int GetNumberOfTargets( TDirectory *dir );
  int GetNumberOfInputVariables( TDirectory *dir );  
  int GetNumberOfInputVariablesMultiClass( TDirectory *dir );

  std::vector<TString> GetInputVariableNames(TDirectory *dir );
  std::vector<TString> GetClassNames(TDirectory *dir );

  TKey* FindMethod( TString name, TDirectory *dir = 0 );
  bool ExistMethodName( TString name, TDirectory *dir = 0 );

  int GetListOfMethods( TList & methods, TDirectory *dir = 0 );
  int GetListOfJobs( TFile* file, TList& jobdirs);
  int GetListOfTitles( TDirectory *rfdir, TList & titles );
  int GetListOfTitles( TString & methodName, TList & titles, TDirectory *dir = 0 );

  TDirectory *GetInputVariablesDir( TMVAGlob::TypeOfPlot type, TDirectory *dir = 0 );
  TDirectory *GetCorrelationPlotsDir( TMVAGlob::TypeOfPlot type, TDirectory *dir = 0 );

  void plotEfficiency (TFile* inputFile, TDirectory* dir, const double & minPTbin = 200, const double & maxPTbin = 1000, const std::string & outputPlotDirectory = "");

  void banner4Plot (const bool & isLabel, const float & ptMin, const float & ptMax);

  void SetMethodName (const std::vector<std::string> & MethodName);

  std::vector<TFile*> inputFiles_ ;
  TCanvas* cROC_ ;
  TH2F*    frameROC_ ;
  TLegend* leg_ ;

  int color_index ;
  int method_index;
  std::vector<std::string> inputMethodName_ ;

};



#endif
 



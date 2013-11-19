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
#include "TMath.h"
#include "TNamed.h"
#include "TFormula.h"

#include "RVersion.h"
#include "TPaletteAxis.h"

// default definition for ROC curve color, style and width                                                                                                                                

static double color []      = {1,2,210,6,7,4,12,95};
static double linestyle []  = {1,1,9,7,1,9,1,9};
static double linewidth []  = {2.6,2.6,2.5,2,2,2,2,2};

static std::vector<double> vec_color(color, color + sizeof(color)/sizeof(double));
static std::vector<double> vec_linewidth(linewidth, linewidth + sizeof(linewidth)/sizeof(double));
static std::vector<double> vec_linestyle(linestyle, linestyle + sizeof(linestyle)/sizeof(double));

class significanceBox {
 
 public:

  significanceBox(){
    Signal_ = NULL ;
    Background_ = NULL ;
    efficiencySignal_ = NULL ;
    efficiencyBackground_ = NULL ;  
    significance_ = NULL;
    maxSig_    = 0; 
    maxSigErr_ = 0; 
    maxbin_    = 0;
  }
  ~significanceBox(){
    if(Signal_!=NULL) delete Signal_;
    if(efficiencySignal_!=NULL) delete efficiencySignal_;
    if(Background_!=NULL) delete Background_;
    if(efficiencyBackground_!=NULL) delete efficiencyBackground_;
    if(significance_!=NULL) delete significance_ ;
  } 

  TH1F* Signal_;
  TH1F* Background_;
  TH1F* efficiencySignal_;
  TH1F* efficiencyBackground_;
  TH1F* significance_;  
  TString methodName_;
  TString methodTitle_;

  double maxSig_ ;
  double maxSigErr_ ;
  int    maxbin_ ;

};

// class for global manipulation of TMVA output + ROC, correlation, output and significance plots

class TMVAGlob {

  public:  
  
  // default constructor
  TMVAGlob(); 
  // default deconstructor
  ~TMVAGlob();
 
  enum TypeOfPlot { kId = 0, kNorm, kDecorrelated, kPCA, kGaussDecorr, kNumOfMethods };

  enum HistType { MVAType = 0, ProbaType = 1, RarityType = 2, CompareType = 3 };

  enum SignificanceType { SoverB = 0, SoverSqrtB = 1, SoverSqrtSB = 2, Pvalue = 3};

  // open one input file given the name
  void openFileInput( const TString& fin );
  // open a set of input file and store them in the class
  void openFileInput( const std::vector<std::string> & fin );

  // get back a single file from the name
  TFile* GetInputFile (const TString& fin );
  // get back the whole set of file considered
  std::vector<TFile*> GetInputFile();


  // create canvas and frame for the ROC plot
  void Initialize( Bool_t useTMVAStyle = kTRUE );
  void CreateCanvasandFrameROC(TFile* inputFile, const double & ptMin, const double & ptMax, const std::string & outputPlotDirectory);
  void DestroyCanvases(); 
  void PrintImageROC( TDirectory* dir = 0 , const std::string & outputPlotDirectory = "");    
  void PrintImage   ( TCanvas* c = 0, const std::string & outputPlotDirectory = "");    


  void NormalizeHist( TH1* h );
  void NormalizeHists( TH1* sig, TH1* bkg = 0 );

  void GetMethodName( TString & name, TKey * mkey );
  void GetMethodTitle( TString & name, TKey * ikey );
  void GetMethodName( TString & name, TDirectory * mdir );
  void GetMethodTitle( TString & name, TDirectory * idir );

  int GetNumberOfTargets( TDirectory *dir );
  int GetNumberOfInputVariables( TDirectory *dir );  
  int GetNumberOfInputVariablesMultiClass( TDirectory *dir );

  int GetListOfMethods( TList & methods, TDirectory *dir = 0 );
  int GetListOfJobs( TFile* file, TList& jobdirs);
  int GetListOfTitles( TDirectory *rfdir, TList & titles );
  int GetListOfTitles( TString & methodName, TList & titles, TDirectory *dir = 0 );

  std::vector<TString> GetInputVariableNames(TDirectory *dir );
  std::vector<TString> GetClassNames(TDirectory *dir );

  // given the iterator on the file and a class name, get the next key
  TKey *NextKey( TIter & keyIter, TString className);
  TKey* FindMethod( TString name, TDirectory *dir = 0 );

  int  GetListOfKeys( TList& keys, TString inherits, TDirectory *dir = 0 );
  bool ExistMethodName( TString name, TDirectory *dir = 0 );

  TDirectory *GetInputVariablesDir( TMVAGlob::TypeOfPlot type, TDirectory *dir = 0 );
  TDirectory *GetCorrelationPlotsDir( TMVAGlob::TypeOfPlot type, TDirectory *dir = 0 );

  void banner4Plot (const bool & isLabel, const float & ptMin, const float & ptMax);
  void SetMethodName (const std::vector<std::string> & MethodName);

  // Methods for significance study

  void ReadHistograms(TFile* file);
  void SetFormula(const TString & f) { fFormula_ = f; }

  TString GetFormula();
  TString GetFormulaString() { return fFormula_; }
  TString GetLatexFormula();

  void SetSignalType    (const bool & type = false) { signalType_ = type; }
  void SetBackgroundType(const bool & type = false) { backgroundType_ = type; }


  // plot ROC curve
  void plotEfficiency (TFile* inputFile, TDirectory* dir, const double & minPTbin = 200, const double & maxPTbin = 1000, const std::string & outputPlotDirectory = "");
  // plot correlation Matrix
  void plotCorrelationMatrix (TFile* inputFile = 0, const int & iFile = 0, const std::string & outputPlotDirectory = "");
  // plot output distribution
  void plotMVAs( TFile*inputFile = 0, HistType htype = MVAType, const std::string & outputPlotDirectory = "");
  // plot signficance for different formula using or not the expected number of signal and background events
  void plotSignificance (TFile* inputFile = 0, const int & iFile = 0, SignificanceType stype = SoverB, const double & numberSignalEvents = 0., const double & numberBackgroundEvents = 0.,
                         const bool & UseSignalEfficiency = false, const bool & UseBakgroundEfficiency = false, const std::string & outputPlotDirectory = "");

 private:

  int color_index ;
  int method_index;
  int mvas_index;
  int thisMethod_;

  std::vector<TFile*> inputFiles_ ;
  std::vector<std::string> inputMethodName_ ;
  std::vector<std::string> originalMethodName_ ;

  TCanvas* cROC_ ;
  TCanvas* cCorrelationSignal_;
  TCanvas* cCorrelationBackground_;
  TCanvas* cMVAs_ ;
  TCanvas* cSignificance_ ;
  TH2F*    frameROC_ ;
  TLegend* legROC_ ;

  TString fFormula_ ;
  bool    signalType_ ;
  bool    backgroundType_ ;

  std::vector<significanceBox*>*  fInfoList_;

  TH1F* histoSignal_;
  TH1F* histoBackground_;
  TH1F* effBackground_;
  TH1F* effSignal_;

};

#endif
 



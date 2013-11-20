#include <iomanip>
#include "TMVAGlob.h"

TMVAGlob::TMVAGlob(){

  cROC_        = NULL;
  frameROC_    = NULL;
  legROC_      = NULL;

  color_index  = 0;
  method_index = 0;
  mvas_index   = 0;
  thisMethod_   = 0 ;
  cMVAs_       = NULL ;

  fInfoList_ = NULL ;

  histoSignal_     = NULL ;
  histoBackground_ = NULL ;
  effSignal_       = NULL ;
  effBackground_   = NULL ;

  signalType_     = false ;
  backgroundType_ = false ;

}

TMVAGlob::~TMVAGlob(){

  if(cROC_!=NULL)       delete cROC_;
  if(frameROC_!=NULL)   delete frameROC_;
  if(legROC_!=NULL)     delete legROC_;
  if(fInfoList_!=NULL){ fInfoList_->clear();  delete fInfoList_; }
  if(histoSignal_!=NULL) delete histoSignal_ ;
  if(histoBackground_!=NULL) delete histoBackground_;
  if(effSignal_!=NULL) delete effSignal_ ;
  if(effBackground_!=NULL) delete effBackground_;
  inputFiles_.clear() ;
  inputMethodName_.clear() ;
  
}

void TMVAGlob::DestroyCanvases(){
  TList* loc = (TList*)gROOT->GetListOfCanvases();
  TListIter itc(loc);
  TObject *o(0);
  while ((o = itc())) delete o;
}

// set style and remove existing canvas'
void TMVAGlob::Initialize( Bool_t useTMVAStyle){
 (*this).DestroyCanvases();
 // set style
 if (!useTMVAStyle) {
     gROOT->SetStyle("Plain");
     gStyle->SetOptStat(0);
     return;
 }
}

// checks if file with name "fin" is already open, and if not opens one
void TMVAGlob::openFileInput( const TString& fin ){
 TFile* file = gDirectory->GetFile();
 if (file==0 || fin != file->GetName()) {
   if (file != 0) {gROOT->cd();
                   file->Close();
   }
  std::cout << "--- Opening root file " << fin << " in read mode" << std::endl;
  file = TFile::Open( fin, "READ" );
 }
 else 
  file = gDirectory->GetFile();
      
 inputFiles_.push_back(file);
 return ;

}

void TMVAGlob::openFileInput( const std::vector<std::string> & fin ){
 for(size_t iFile = 0 ; iFile < fin.size(); iFile++){
  std::cout << "--- Opening root file " << fin.at(iFile) << " in read mode" << std::endl;
  inputFiles_.push_back(new TFile( fin.at(iFile).c_str(), "READ"));
 }

}


TFile* TMVAGlob::GetInputFile(const TString& fin){
  for(size_t iFile = 0 ; iFile < inputFiles_.size(); iFile++){
    if(inputFiles_.at(iFile)->GetName() == fin) 
      return inputFiles_.at(iFile) ;
  }

  return 0 ;
}

std::vector<TFile*> TMVAGlob::GetInputFile(){
  return inputFiles_ ;
}



TKey *TMVAGlob::NextKey( TIter & keyIter, TString className) {

 TKey *key=(TKey *)keyIter.Next();
 TKey *rkey=0;
 Bool_t loop=(key!=0);

 while (loop) {

   TClass *cl = gROOT->GetClass(key->GetClassName());
   if (cl->InheritsFrom(className.Data())) { loop = kFALSE;
                                             rkey = key;
   } 
   else {
          key = (TKey *)keyIter.Next();
          if (key==0) loop = kFALSE;
         }
   }

return rkey;

}


// used to create output file for canvas
void TMVAGlob::PrintImageROC(TDirectory* dir, const std::string & outputPlotDirectory){

  TString fname = outputPlotDirectory + "ROCcurve";

  if (TString(dir->GetName()).Contains("multicut")){
    TString fprepend(dir->GetName());
    fprepend.ReplaceAll("multicutMVA_","");
    fname = outputPlotDirectory + fprepend + "_" + "ROCcurve";
  }

  if (NULL == (*this).cROC_) std::cout << "*** Error in TMVAGlob::imgconv: canvas is NULL" << std::endl;
  else {
       // create directory if not existing
       TString f = fname;
       TString dir = f.Remove( f.Last( '/' ), f.Length() - f.Last( '/' ) );
       gSystem->mkdir( dir );

       TString pngName = fname + ".png";
       TString pdfName = fname + ".pdf";
       (*this).cROC_->cd();

       (*this).cROC_->Print(pdfName);
       (*this).cROC_->Print(pngName);
 }
}

// function to normalize a histo also if it has weights
void TMVAGlob::NormalizeHist( TH1F* h ) { 
 if (h==0) return;
 if (h->GetSumw2N() == 0) h->Sumw2();
 if(h->GetSumOfWeights()!=0) {
      Float_t dx = (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin())/h->GetNbinsX();
      h->Scale( 1.0/h->GetSumOfWeights()/dx );
  }
}

void TMVAGlob::NormalizeHists( TH1F* sig, TH1F* bkg) {
 if (sig->GetSumw2N() == 0) sig->Sumw2();
 if (bkg && bkg->GetSumw2N() == 0) bkg->Sumw2();
      
 if(sig->GetSumOfWeights()!=0) {
     Float_t dx = (sig->GetXaxis()->GetXmax() - sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
     sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
 }
 if (bkg != 0 && bkg->GetSumOfWeights()!=0) {
     Float_t dx = (bkg->GetXaxis()->GetXmax() - bkg->GetXaxis()->GetXmin())/bkg->GetNbinsX();
     bkg->Scale( 1.0/bkg->GetSumOfWeights()/dx );
 }
}

// the following are tools to help handling different methods and titles
void TMVAGlob::GetMethodName( TString & name, TKey * mkey ) {
   if (mkey==0) return;
   name = mkey->GetName();
   name.ReplaceAll("Method_","");
}

void TMVAGlob::GetMethodTitle( TString & name, TKey * ikey ) {
   if (ikey==0) return;
   name = ikey->GetName();
}

void TMVAGlob::GetMethodName( TString & name, TDirectory * mdir ) {
   if (mdir==0) return;
   name = mdir->GetName();
   name.ReplaceAll("Method_","");
}

void TMVAGlob::GetMethodTitle( TString & name, TDirectory * idir ) {
   if (idir==0) return;
   name = idir->GetName();
}


// get a list of keys with a given inheritance
// the list contains TKey objects
int TMVAGlob::GetListOfKeys( TList& keys, TString inherits, TDirectory *dir){

 if (dir==0) dir = gDirectory;
 TIter mnext(dir->GetListOfKeys());
 TKey *mkey;
 keys.Clear();
 keys.SetOwner(kFALSE);
 UInt_t ni=0;
 
 while ((mkey = (TKey*)mnext())) {
         
   // make sure, that we only look at TDirectory with name Method_<xxx>
   TClass *cl = gROOT->GetClass(mkey->GetClassName());
   if (cl->InheritsFrom(inherits)) {
      keys.Add(mkey);
      ni++;
   }
 }
 
 return ni;
}

int TMVAGlob::GetNumberOfTargets( TDirectory *dir ){
  if (!dir) {
             std::cout << "tmvaglob::GetNumberOfTargets is called with *dir==NULL :( " << std::endl;
             return 0;
  }
  
  TIter next(dir->GetListOfKeys());
  TKey* key    = 0;
 Int_t noTrgts = 0;
      
   while ((key = (TKey*)next())) {
     if (key->GetCycle() != 1) continue;        
     if (TString(key->GetName()).Contains("__Regression_target")) noTrgts++;
   }
   return noTrgts;
}


int TMVAGlob::GetNumberOfInputVariables( TDirectory *dir ){

 TIter next(dir->GetListOfKeys());
 TKey* key    = 0;
 Int_t noVars = 0;
   
 while ((key = (TKey*)next())) {
   if (key->GetCycle() != 1) continue;
   // count number of variables (signal is sufficient), exclude target(s)
   if (TString(key->GetName()).Contains("__Signal") || (TString(key->GetName()).Contains("__Regression") && !(TString(key->GetName()).Contains("__Regression_target")))) noVars++;
 }
      
 return noVars;
}


std::vector<TString> TMVAGlob::GetInputVariableNames(TDirectory *dir ){

 TIter next(dir->GetListOfKeys());
 TKey* key = 0;
 std::vector<TString> names;
      
 while ((key = (TKey*)next())) {
   if (key->GetCycle() != 1) continue;
   TClass *cl = gROOT->GetClass(key->GetClassName());
   if (!cl->InheritsFrom("TH1")) continue;
   TString name(key->GetName());
   Int_t pos = name.First("__");
   name.Remove(pos);
   Bool_t hasname = false;
   std::vector<TString>::const_iterator iter = names.begin();
   while(iter != names.end()){
     if(name.CompareTo(*iter)==0)
        hasname=true;
     iter++;
   }
   if(!hasname) names.push_back(name);  
 }    
 return names;
}


int TMVAGlob::GetNumberOfInputVariablesMultiClass( TDirectory *dir ){
  std::vector<TString> names(GetInputVariableNames(dir));
  return names.end() - names.begin();
}
   
std::vector<TString> TMVAGlob::GetClassNames(TDirectory *dir ){      
      
 TIter next(dir->GetListOfKeys());
 TKey* key = 0;
 std::vector<TString> names;
      
 while ((key = (TKey*)next())) {
   if (key->GetCycle() != 1) continue;
   TClass *cl = gROOT->GetClass(key->GetClassName());
   if (!cl->InheritsFrom("TH1")) continue;
   TString name(key->GetName());
   name.ReplaceAll("_Deco","");
   name.ReplaceAll("_Gauss","");
   name.ReplaceAll("_PCA","");
   name.ReplaceAll("_Id","");
   name.ReplaceAll("_vs_","");
   char c = '_';
   Int_t pos = name.Last(c);
   name.Remove(0,pos+1);
         
   Bool_t hasname = false;
   std::vector<TString>::const_iterator iter = names.begin();
   while(iter != names.end()){
       if(name.CompareTo(*iter)==0)
          hasname=true;
       iter++;
   }
   if(!hasname) names.push_back(name);
 }
 return names;
}


// find the key for a method
TKey* TMVAGlob::FindMethod( TString name, TDirectory *dir){
  
 if (dir==0) dir = gDirectory;
 TIter mnext(dir->GetListOfKeys());
 TKey *mkey;
 TKey *retkey=0;
 Bool_t loop=kTRUE;
 while (loop) {
   mkey = (TKey*)mnext();
   if (mkey==0) loop = kFALSE;
   else {
         TString clname = mkey->GetClassName();
         TClass *cl = gROOT->GetClass(clname);
         if (cl->InheritsFrom("TDirectory")) {
             TString mname = mkey->GetName(); // method name
             TString tname = "Method_"+name;  // target name
             if (mname==tname) { // target found!
                  loop = kFALSE;
                  retkey = mkey;
             }  
         }
   }
 }
 
 return retkey;
}


// find the key for a method
bool TMVAGlob::ExistMethodName( TString name, TDirectory *dir){
  
 if (dir==0) dir = gDirectory;
 TIter mnext(dir->GetListOfKeys());
 TKey *mkey;
 Bool_t loop=kTRUE;
 while (loop) {
    mkey = (TKey*)mnext();
    if (mkey==0) loop = kFALSE;
    else {
          TString clname  = mkey->GetClassName();
          TString keyname = mkey->GetName();
          TClass *cl = gROOT->GetClass(clname);
          if (keyname.Contains("Method") && cl->InheritsFrom("TDirectory")) {
               TDirectory* d_ = (TDirectory*)dir->Get( keyname );
               if  (!d_) {
	 	  std::cout << "HUUUGE TROUBLES IN TMVAGlob::ExistMethodName() --> contact authors" << std::endl;
                  return kFALSE;
               }
               TIter mnext_(d_->GetListOfKeys());
               TKey *mkey_;
               while ((mkey_ = (TKey*)mnext_())) {
                  TString clname_ = mkey_->GetClassName();
                  TClass *cl_ = gROOT->GetClass(clname_);
                  if (cl_->InheritsFrom("TDirectory")) {
                     TString mname = mkey_->GetName(); // method name
                     if (mname==name) return kTRUE;
		  }
	       }
	  }
    }
 }
   
 return kFALSE;
   
}


int TMVAGlob::GetListOfMethods( TList & methods, TDirectory *dir){

 // get a list of methods
 // the list contains TKey objects
 if (dir==0) dir = gDirectory;
 TIter mnext(dir->GetListOfKeys());
 TKey *mkey;
 methods.Clear();
 methods.SetOwner(kFALSE);
 UInt_t ni=0;
 while ((mkey = (TKey*)mnext())) {
     // make sure, that we only look at TDirectory with name Method_<xxx>
     TString name = mkey->GetClassName();
     TClass *cl = gROOT->GetClass(name);
     if (cl->InheritsFrom("TDirectory")) {
      if (TString(mkey->GetName()).BeginsWith("Method_")) {
         methods.Add(mkey);
         ni++;
      }
     }
 }
      
 std::cout << "--- Found " << ni << " classifier types" << std::endl;
 return ni;
}


// get a list of all jobs in all method directories
// based on ideas by Peter and Joerg found in macro deviations.C
int TMVAGlob::GetListOfJobs( TFile* file, TList& jobdirs){

  TIter next(file->GetListOfKeys());
  TKey *key(0);   
  while ((key = (TKey*)next())) {         
   if (TString(key->GetName()).BeginsWith("Method_")) {
    if (gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) {
        TDirectory* mDir = (TDirectory*)key->ReadObj();
        TIter keyIt(mDir->GetListOfKeys());
        TKey *jobkey;
        while ((jobkey = (TKey*)keyIt())) {
          if (!gROOT->GetClass(jobkey->GetClassName())->InheritsFrom("TDirectory")) continue;
          TDirectory *jobDir = (TDirectory *)jobkey->ReadObj();
	  std::cout << "jobdir name  " << jobDir->GetName() << std::endl;
          jobdirs.Add(jobDir);
	}
    }
   }
  }
  
  return jobdirs.GetSize();
}

 // get a list of titles (i.e TDirectory) given a method dir
int TMVAGlob::GetListOfTitles( TDirectory *rfdir, TList & titles ){
 UInt_t ni=0;
 if (rfdir==0) return 0;
 TList *keys = rfdir->GetListOfKeys();
 if (keys==0) {
   std::cout << "+++ Directory '" << rfdir->GetName() << "' contains no keys" << std::endl;
   return 0;
 }
 
 TIter rfnext(rfdir->GetListOfKeys());
 TKey *rfkey;
 titles.Clear();
 titles.SetOwner(kFALSE);
 while ((rfkey = (TKey*)rfnext())) {
    // make sure, that we only look at histograms
    TClass *cl = gROOT->GetClass(rfkey->GetClassName());
    if (cl->InheritsFrom("TDirectory")) {
        titles.Add(rfkey);
         ni++;
    }
 }
 std::cout << "--- Found " << ni << " instance(s) of the method " << rfdir->GetName() << std::endl;
 return ni;
}

// get the list of all titles for a given method
// if the input dir is 0, gDirectory is used
// returns a list of keys
int TMVAGlob::GetListOfTitles( TString & methodName, TList & titles, TDirectory *dir){

  UInt_t ni=0;
  if (dir==0) dir = gDirectory;
  TDirectory* rfdir = (TDirectory*)dir->Get( methodName );
  if (rfdir==0) {
    std::cout << "+++ Could not locate directory '" << methodName << std::endl;
    return 0;
  }

  return (*this).GetListOfTitles( rfdir, titles );

  TList *keys = rfdir->GetListOfKeys();
  if (keys==0) {
    std::cout << "+++ Directory '" << methodName << "' contains no keys" << std::endl;
    return 0;
  }
      
  TIter rfnext(rfdir->GetListOfKeys());
  TKey *rfkey;
  titles.Clear();
  titles.SetOwner(kFALSE);
   while ((rfkey = (TKey*)rfnext())) {
    // make sure, that we only look at histograms
    TClass *cl = gROOT->GetClass(rfkey->GetClassName());
    if (cl->InheritsFrom("TDirectory")) {
       titles.Add(rfkey);
        ni++;
    }
  }
  std::cout << "--- Found " << ni << " instance(s) of the method " << methodName << std::endl;
  return ni;
}


TDirectory *TMVAGlob::GetInputVariablesDir( TMVAGlob::TypeOfPlot type, TDirectory *dir){

 // get the InputVariables directory
 const TString directories[TMVAGlob::kNumOfMethods] = { "InputVariables_Id",
                                                        "InputVariables_Deco",
                                                        "InputVariables_PCA",
                                                        "InputVariables_Gauss_Deco" };
 if (dir==0) dir = gDirectory;

 // get top dir containing all hists of the variables
 dir = (TDirectory*)gDirectory->Get( directories[type] );
 if (dir==0) {
   std::cout << "+++ Could not locate input variable directory '" << directories[type] << std::endl;
   return 0;
 }
 
 return dir;

}

TDirectory *TMVAGlob::GetCorrelationPlotsDir( TMVAGlob::TypeOfPlot type, TDirectory *dir){
  // get the CorrelationPlots directory
  if (dir==0) dir = (*this).GetInputVariablesDir(type);
  if (dir==0) return 0;
  TDirectory* corrdir = (TDirectory*)dir->Get( "CorrelationPlots" );
  if (corrdir==0) {
    std::cout << "+++ Could not find CorrelationPlots directory 'CorrelationPlots'" << std::endl;
    return 0;
  }
  return corrdir;
}


void TMVAGlob::banner4Plot (const bool & isLabel, const float & ptMin, const float & ptMax){

  TPaveText* pt = new TPaveText(.76,0.71,.83,.88,"NDC");
  pt->AddText("CA R = 0.8");
  TString BoostLegend ; BoostLegend.Form("%d < p_{T} < %d GeV",int(ptMin),int(ptMax));
  pt->AddText(BoostLegend.Data());
  pt->AddText("p_{T} > 200 GeV");
  pt->AddText("|#eta|<2.4");
  pt->AddText("60 < m_{j} < 100 GeV");

  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetTextSize(0.027);
  pt->SetFillStyle(0);
  pt->SetLineColor(1);
  pt->SetLineStyle(1);
  pt->SetLineWidth(0);
  pt->SetMargin(0);
  pt->SetTextFont(42);
  pt->SetShadowColor(0);
  pt->Draw();

}


void TMVAGlob::CreateCanvasandFrameROC(TFile *inputFile, const double & minPTbin, const double & maxPTbin, const std::string & outputPlotDirectory){
                               

  cROC_ = new TCanvas( (std::string("cROC")+std::string(inputFile->GetName())).c_str(), "cROC",180,52,550,550);
  cROC_->cd();
  cROC_->SetTicks();
  cROC_->SetFillColor(0);
  cROC_->SetBorderMode(0);
  cROC_->SetBorderSize(2);
  cROC_->Range(-0.128266,-0.1538462,1.059382,1.128205);
  
  cROC_->SetTickx(1);
  cROC_->SetTicky(1);
  cROC_->SetRightMargin(0.05);
  cROC_->SetBottomMargin(0.12);
  cROC_->SetFrameBorderMode(0);


  frameROC_ = new TH2F((std::string("frameROC")+std::string(inputFile->GetName())).c_str(),"",500,0,1,500,0,1);
  frameROC_->SetLineWidth(2);
  frameROC_->SetMarkerStyle(21);
  frameROC_->SetMarkerSize(0.3);
  frameROC_->GetXaxis()->SetNdivisions(405);
  frameROC_->GetYaxis()->SetNdivisions(405);
  frameROC_->GetXaxis()->SetTitle("#epsilon_{sig}");
  frameROC_->GetXaxis()->SetLabelOffset(0.012);
  frameROC_->GetXaxis()->SetLabelSize(0.042);
  frameROC_->GetXaxis()->SetTitleSize(0.05);
  frameROC_->GetXaxis()->SetTitleOffset(1.05);
  frameROC_->GetYaxis()->SetTitle("1 - #epsilon_{bkg}");
  frameROC_->GetYaxis()->SetLabelOffset(0.012);
  frameROC_->GetYaxis()->SetLabelSize(0.042);
  frameROC_->GetYaxis()->SetTitleSize(0.05);
  frameROC_->GetYaxis()->SetTitleOffset(1.25);
  frameROC_->Draw("");

  banner4Plot(false,minPTbin,maxPTbin);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(21); // align right                                                                                                                                                  
  latex.SetTextSize(0.033);
  latex.DrawLatex(0.6,0.92,Form("CMS Preliminary Simulation,#sqrt{s} = 8 TeV, W+jets"));


  legROC_ = new TLegend(0.20,0.171,0.6,0.511,NULL,"brNDC");
  legROC_->SetBorderSize(0);
  legROC_->SetTextSize(0.033);
  legROC_->SetTextFont(42);
  legROC_->SetLineColor(1);
  legROC_->SetLineStyle(1);
  legROC_->SetLineWidth(2);
  legROC_->SetFillColor(19);
  legROC_->SetFillStyle(0);
}

void TMVAGlob::plotEfficiency (TFile* inputFile, TDirectory* dir, const double & minPTbin, const double & maxPTbin, const std::string & outputPlotDirectory){

  // Plot the ROC curve with a proper style from root file originated by TMVA                                                                                                         

  if((*this).cROC_==NULL) (*this).CreateCanvasandFrameROC(inputFile,minPTbin,maxPTbin,outputPlotDirectory); 
      
  TList TrainingMethods;
  TList hists;
  int res = (*this).GetListOfMethods(TrainingMethods);
  TIter next(&TrainingMethods);
  TKey *key = 0, *hkey = 0;
  inputFile->cd();
  (*this).cROC_->cd();

  // loop over all methods                                                                                                                                                                
  while ((key = (TKey*)next())) {
    TDirectory * myDir = (TDirectory*)key->ReadObj();
    TList Titles;
    int nTitles = (*this).GetListOfTitles(myDir,Titles);
    TIter nextTitle(&Titles);
    TKey *titkey = 0;
    TDirectory *titDir = 0;
    while ((titkey = (*this).NextKey(nextTitle,"TDirectory"))) {
      titDir = (TDirectory*)titkey->ReadObj();
      TString methodTitle;
      (*this).GetMethodTitle(methodTitle,titDir);
      TIter nextKey( titDir->GetListOfKeys() );
      while ((hkey = (*this).NextKey(nextKey,"TH1"))) {
        TH1F *h = (TH1F*)hkey->ReadObj();
        TString hname = h->GetName();
        if (hname.Contains("rejBvsS") && hname.BeginsWith("MVA_")) {
          if(size_t((*this).color_index) <= vec_color.size()){
	    h->SetLineWidth(vec_linewidth[(*this).color_index]);
	    h->SetLineColor(vec_color[(*this).color_index]);
	    h->SetLineStyle(vec_linestyle[(*this).color_index]);
	    h->Draw("csame");
	    hists.Add(h);
            (*this).color_index = (*this).color_index+1;
          }
          else{
	    h->SetLineWidth(vec_linewidth[(*this).color_index-vec_color.size()]);
	    h->SetLineColor(vec_color[(*this).color_index-vec_color.size()]);
	    h->SetLineStyle(vec_linestyle[(*this).color_index-vec_color.size()]);
	    h->Draw("csame");
	    hists.Add(h);
            (*this).color_index = (*this).color_index+1;
          }
        }
      }
    }
  }
  
  /// Loop on the different histos                                                                                                                                                        
  while (hists.GetSize()) {
    TListIter hIt(&hists); // define an iterator                                                                                                                                          
    TH1F* hist(0);
    Double_t largestInt=-1;
    TH1F* histWithLargestInt(0);
    while ((hist = (TH1F*)hIt())!=0) {
      Double_t integral = hist->Integral(1,hist->GetNbinsX());
      if (integral>largestInt) {
        largestInt = integral;
        histWithLargestInt = hist;
      }
    }
    if (histWithLargestInt == 0) {
      std::cout << "ERROR - unknown hist \"histWithLargestInt\" --> serious problem in ROOT file" << std::endl;
      break;
    }    
 
   if(TString(histWithLargestInt->GetTitle()).Contains("Cuts")){
     (*this).legROC_->AddEntry(histWithLargestInt,inputMethodName_.at((*this).method_index).c_str(),"l");
     (*this).method_index ++ ;
   }
   else if (TString(histWithLargestInt->GetTitle()).Contains("Likelihood"))
    (*this).legROC_->AddEntry(histWithLargestInt,"Likelihood","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("LD"))
    (*this).legROC_->AddEntry(histWithLargestInt,"Linear Discriminant","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("BDT") && ! TString(histWithLargestInt->GetTitle()).Contains("BDTG") && !TString(histWithLargestInt->GetTitle()).Contains("BDTF"))
    (*this).legROC_->AddEntry(histWithLargestInt,"Boosted Decision Tree (BDT)","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("BDTG"))
    (*this).legROC_->AddEntry(histWithLargestInt,"Gradient BDT (BDTG)","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("MLP"))
    (*this).legROC_->AddEntry(histWithLargestInt,"Multi-Layer Perceptron (MLP)","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("PDEFoam"))
    (*this).legROC_->AddEntry(histWithLargestInt,"PDEFoam","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("SVN"))
    (*this).legROC_->AddEntry(histWithLargestInt,"Supported Vector Machine (SVM)","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("Fisher"))
    (*this).legROC_->AddEntry(histWithLargestInt,"Fisher Discriminant","l");
   else if (TString(histWithLargestInt->GetTitle()).Contains("Fisher"))
     (*this).legROC_->AddEntry(histWithLargestInt,TString(histWithLargestInt->GetTitle()).ReplaceAll("MVA_",""),"l");   
   hists.Remove(histWithLargestInt);
  }
  
  (*this).legROC_->Draw("same");
  (*this).cROC_->Update();
  
  return;

}


void TMVAGlob::SetMethodName(const std::vector<std::string> & SetMethodName){

  (*this).inputMethodName_ = SetMethodName;
  (*this).originalMethodName_ = SetMethodName;

  for(size_t iMethod = 0 ; iMethod < (*this).inputMethodName_.size() ; iMethod++){
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll(":_:","#"));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll(":__:","_{"));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll(":___:","}"));
    inputMethodName_.at(iMethod) = std::string(TString(inputMethodName_.at(iMethod)).ReplaceAll("//"," "));
  }       

  for(size_t iMethod = 0 ; iMethod < (*this).originalMethodName_.size() ; iMethod++){
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll(":_:","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll(":__:","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll(":___:","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll("//","_"));
    originalMethodName_.at(iMethod) = std::string(TString(originalMethodName_.at(iMethod)).ReplaceAll("/","_"));
  }       
}


void TMVAGlob::plotCorrelationMatrix(TFile* inputFile, const int & iFile, const std::string & outputPlotDirectory){
  
  std::string nameCorrelationSignal     = "CorrelationMatrixS";
  std::string nameCorrelationBackground = "CorrelationMatrixB";

  TH2F* hSignal     = (TH2F*) inputFile->Get(nameCorrelationSignal.c_str());
  TH2F* hBackground = (TH2F*) inputFile->Get(nameCorrelationBackground.c_str());
  hSignal->SetName(std::string(Form("%s_%d",nameCorrelationSignal.c_str(),iFile)).c_str());
  hBackground->SetName(std::string(Form("%s_%d",nameCorrelationBackground.c_str(),iFile)).c_str());
  
  if(hSignal == 0 || hBackground == 0){ std::cerr<<" Null Pointer for correlation Matrix --> exit without plot "<<std::endl;  return ; }

  (*this).cCorrelationSignal_ = new TCanvas(std::string(Form("c%s_%d",nameCorrelationSignal.c_str(),iFile)).c_str(),Form("Correlation Matrix Signal"),180,52,550,550);
  float newMargin1 = 0.13;
  float newMargin2 = 0.15;
  float newMargin3 = 0.20;

  (*this).cCorrelationSignal_->SetGrid();
  (*this).cCorrelationSignal_->SetTicks();
  (*this).cCorrelationSignal_->SetLeftMargin(newMargin3);
  (*this).cCorrelationSignal_->SetBottomMargin(newMargin2);
  (*this).cCorrelationSignal_->SetRightMargin(newMargin1);
  (*this).cCorrelationSignal_->SetTopMargin(newMargin1);
  
  gStyle->SetPaintTextFormat("3g");

  hSignal->SetMarkerSize(1.5);
  hSignal->SetMarkerColor(0);
  hSignal->GetXaxis()->SetLabelSize(0.035);
  hSignal->GetYaxis()->SetLabelSize(0.035);
  hSignal->LabelsOption("d");
  hSignal->SetLabelOffset(0.011);// label offset on x axis                                                                                                                                
  
  hBackground->SetMarkerSize(1.5);
  hBackground->SetMarkerColor(0);
  hBackground->GetXaxis()->SetLabelSize(0.035);
  hBackground->GetYaxis()->SetLabelSize(0.035);
  hBackground->LabelsOption("d");
  hBackground->SetLabelOffset( 0.011 );// label offset on x axis                                                                                                                          
 
  // Plot correlation between signal      
  (*this).cCorrelationSignal_->cd();
  hSignal->Draw("colz"); // color pads                                                                                                                                                    
  hSignal->Draw("textsame");  // add text                                                                                                                                  
  // add comment                                                                                                                                                                       
  TText* text = new TText( 0.53, 0.88, "Linear correlation coefficients in %" );
  text->SetNDC();
  text->SetTextSize( 0.026 );
  text->AppendPad();
  (*this).cCorrelationSignal_->Update();

  nameCorrelationSignal = std::string(Form("%s_%d",nameCorrelationSignal.c_str(),iFile));
  nameCorrelationBackground = std::string(Form("%s_%d",nameCorrelationBackground.c_str(),iFile));

  (*this).PrintImage((*this).cCorrelationSignal_,outputPlotDirectory+"/"+nameCorrelationSignal);
    
  // Background correlation
  (*this).cCorrelationBackground_ = new TCanvas(std::string(Form("c%s_%d",nameCorrelationBackground.c_str(),iFile)).c_str(),Form("Correlation Matrix Signal"),180,52,550,550);
  (*this).cCorrelationBackground_->SetGrid();
  (*this).cCorrelationBackground_->SetTicks();
  (*this).cCorrelationBackground_->SetLeftMargin(newMargin3);
  (*this).cCorrelationBackground_->SetBottomMargin(newMargin2);
  (*this).cCorrelationBackground_->SetRightMargin(newMargin1);
  (*this).cCorrelationBackground_->SetTopMargin(newMargin1);

  (*this).cCorrelationBackground_->cd();
  hBackground->Draw("colz"); // color pads                                                                                                                             
  hBackground->Draw("textsame");  // add text                                                                                                                             
  // add comment                                                                                                                                                                       
  text->AppendPad();
  (*this).cCorrelationBackground_->Update();
  (*this).PrintImage((*this).cCorrelationBackground_,outputPlotDirectory+"/"+nameCorrelationBackground);

  if(hSignal!=0) delete hSignal;
  if(hBackground!=0) delete hBackground;
  if(text!=0) delete text;

}

void TMVAGlob::PrintImage(TCanvas* c, const std::string & fname){

  TString pngName = fname+".png";
  TString pdfName = fname+".pdf";

  c->Print(pdfName);
  c->Print(pngName);

}

void TMVAGlob::plotMVAs(TFile* inputFile, HistType htype, const std::string & outputPlotDirectory){


  // search for the right histograms in full list of keys                                                                                                                                 
  TIter next(inputFile->GetListOfKeys());
  TKey *key(0);

  while ((key = (TKey*)next())) {

    if (!TString(key->GetName()).BeginsWith("Method_")) continue;
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

    TString methodName;
    (*this).GetMethodName(methodName,key);

    TDirectory* mDir = (TDirectory*)key->ReadObj();

    TIter keyIt(mDir->GetListOfKeys());
    TKey *titkey;

    while ((titkey = (TKey*)keyIt())) {

      if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;

      TDirectory *titDir = (TDirectory *)titkey->ReadObj();
      TString methodTitle;
      (*this).GetMethodTitle(methodTitle,titDir);
      std::cout << "--- Found directory for method: " << methodName << "::" << methodTitle << std::flush;

      TString hname = "MVA_" + methodTitle;
      if      (htype == ProbaType  ) hname += "_Proba";
      else if (htype == RarityType ) hname += "_Rarity";

      (*this).histoSignal_ = dynamic_cast<TH1*>(titDir->Get(hname+"_S" ));
      (*this).histoBackground_ = dynamic_cast<TH1*>(titDir->Get(hname+"_B"));

      if ((*this).histoSignal_==0 || (*this).histoBackground_==0) {
	if (htype == MVAType)
	  std::cout << ":\t mva distribution not available (this is normal for Cut classifier)" << std::endl;
	else if(htype == ProbaType)
	  std::cout << ":\t probability distribution not available" << std::endl;
	else if(htype == RarityType)
	  std::cout << ":\t rarity distribution not available" << std::endl;
	else if(htype == CompareType)
	  std::cout << ":\t overtraining check not available" << std::endl;
	else std::cout << std::endl;
	continue;
      }

      std::cout << " containing " << hname << "_S/_B" << std::endl;

      if((methodTitle).Contains("Likelihood"))
         methodTitle = Form("Likelhood");
      else if ((methodTitle).Contains("LD"))
         methodTitle = Form("Linear Discriminant");
      else if ((methodTitle).Contains("BDT") && !(methodTitle).Contains("BDTG") && !(methodTitle).Contains("BDTF"))
         methodTitle = Form("Boosted Decision Tree (BDT)");
      else if ((methodTitle).Contains("BDTG")) 
         methodTitle = Form("Gradient BDT (BDTG)");
      else if ((methodTitle).Contains("MLP"))
         methodTitle = Form("Multi-Layer Perceptron (MLP)");
      else if ((methodTitle).Contains("PDEFoam")) 
         methodTitle = Form("PDEFoam");
      else if ((methodTitle).Contains("SVN")) 
         methodTitle = Form("Supported Vector Machine (SVM)");
      else if ((methodTitle).Contains("Fisher")) 
         methodTitle = Form("Fisher Discriminant");
        
      (*this).histoSignal_->SetTitle( Form("Response for classifier: %s", methodTitle.Data()) );
      if(htype == ProbaType)
	(*this).histoSignal_->SetTitle( Form("Probability for classifier: %s", methodTitle.Data()) );
      else if (htype == RarityType)
	(*this).histoSignal_->SetTitle( Form("Rarity for classifier: %s", methodTitle.Data()) );
      else if (htype == CompareType)
	(*this).histoSignal_->SetTitle( Form("Overtraining check for classifier: %s", methodTitle.Data()) );                                                                                


      TString cCanvasTitle = ((htype == MVAType)     ? Form("Response %s",methodTitle.Data()) :
			      (htype == ProbaType)   ? Form("Probability %s",methodTitle.Data()) :
			      (htype == CompareType) ? Form("Comparison %s",methodTitle.Data()) :
			      Form("Rarity %s",methodTitle.Data()));


      (*this).cMVAs_ = new TCanvas( Form("cMVAs_%d",(*this).mvas_index), cCanvasTitle,180,100,550,550);
      (*this).cMVAs_->cd();
      (*this).cMVAs_->SetTicks();
      (*this).cMVAs_->SetFillColor(0);
      (*this).cMVAs_->SetBorderMode(0);
      (*this).cMVAs_->SetBorderSize(2);
  
      (*this).cMVAs_->SetTickx(1);
      (*this).cMVAs_->SetTicky(1);
      (*this).cMVAs_->SetRightMargin(0.05);
      (*this).cMVAs_->SetBottomMargin(0.12);
      (*this).cMVAs_->SetFrameBorderMode(0);

      // set the histogram style                                                                                                                                                        
      (*this).histoSignal_->SetFillColor(4);
      (*this).histoSignal_->SetLineColor(4);
      (*this).histoSignal_->SetLineWidth(2);
      (*this).histoSignal_->SetFillStyle(1001);

      if(htype == CompareType) (*this).histoSignal_->SetFillStyle(3001);
         
      (*this).histoBackground_->SetFillColor(2);
      (*this).histoBackground_->SetLineColor(2);
      (*this).histoBackground_->SetLineWidth(2);
      (*this).histoBackground_->SetFillStyle(3005);
    
      // normalise both signal and background                                                                                                                                           
      (*this).NormalizeHists( (TH1F*)(*this).histoSignal_, (TH1F*)(*this).histoBackground_);

      // frame limits (choose judicuous x range)                                                                                                                                        
      Float_t nrms = 10;
      std::cout << "--- Mean and RMS (S): " << (*this).histoSignal_->GetMean() << ", " << (*this).histoSignal_->GetRMS() << std::endl;
      std::cout << "--- Mean and RMS (B): " << (*this).histoBackground_->GetMean() << ", " << (*this).histoBackground_->GetRMS() << std::endl;

      Float_t xmin = TMath::Max( TMath::Min((*this).histoSignal_->GetMean() - nrms*(*this).histoSignal_->GetRMS(),(*this).histoBackground_->GetMean() - nrms*(*this).histoBackground_->GetRMS() ),(*this).histoSignal_->GetXaxis()->GetXmin() );
      Float_t xmax = TMath::Min( TMath::Max((*this).histoSignal_->GetMean() + nrms*(*this).histoSignal_->GetRMS(),(*this).histoBackground_->GetMean() + nrms*(*this).histoBackground_->GetRMS() ),(*this).histoSignal_->GetXaxis()->GetXmax() );
      Float_t ymin = 0;
      Float_t maxMult = (htype == CompareType) ? 1.3 : 1.2;
      Float_t ymax = TMath::Max( (*this).histoSignal_->GetMaximum(), (*this).histoBackground_->GetMaximum() )*maxMult;

      // build a frame                                                                                                                                                                  
      Int_t nb = 500;
      TString hFrameName(TString("frame") + methodTitle);
      TObject *o = gROOT->FindObject(hFrameName);
      if(o) delete o;
      TH2F* frame = new TH2F( hFrameName, (*this).histoSignal_->GetTitle(),nb, xmin, xmax, nb, ymin, ymax );
      frame->GetXaxis()->SetTitle( methodTitle + ((htype == MVAType || htype == CompareType) ? " response" : "") );
      if      (htype == ProbaType  ) frame->GetXaxis()->SetTitle( "Signal probability" );
      else if (htype == RarityType ) frame->GetXaxis()->SetTitle( "Signal rarity" );
      frame->GetYaxis()->SetTitle("(1/N) dN/dx");


      frame->SetLineWidth(2);
      frame->SetMarkerStyle(21);
      frame->SetMarkerSize(0.3);
      frame->GetXaxis()->SetNdivisions(405);
      frame->GetYaxis()->SetNdivisions(405);
      frame->GetXaxis()->SetLabelOffset(0.012);
      frame->GetXaxis()->SetLabelSize(0.03);
      frame->GetXaxis()->SetTitleSize(0.035);
      frame->GetXaxis()->SetTitleOffset(1.05);
      frame->GetYaxis()->SetLabelOffset(0.012);
      frame->GetYaxis()->SetLabelSize(0.03);
      frame->GetYaxis()->SetTitleSize(0.035);
      frame->GetYaxis()->SetTitleOffset(1.10);
      frame->Draw("");
   
      // Draw legend                                                                                                                                                                    
      TLegend *legend= new TLegend( (*this).cMVAs_->GetLeftMargin()+0.03, 1-(*this).cMVAs_->GetTopMargin()-0.15,
				    (*this).cMVAs_->GetLeftMargin()+(htype == CompareType ? 0.40 : 0.3)+0.15, 1-(*this).cMVAs_->GetTopMargin()-0.03);
      legend->SetFillStyle(0);
      legend->SetFillColor(0);
      legend->SetTextSize (0.042);
      if(htype == CompareType) legend->SetTextSize (0.032); 
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetLineColor(1);
      legend->SetLineStyle(1);
      legend->SetLineWidth(2);

      legend->AddEntry((*this).histoSignal_,TString("Signal")     + ((htype == CompareType) ? " (testing)" : ""), "F");
      legend->AddEntry((*this).histoBackground_,TString("Background") + ((htype == CompareType) ? " (testing)" : ""), "F");

      legend->SetMargin( (htype == CompareType ? 0.2 : 0.3) );
      legend->Draw("same");

      (*this).histoSignal_->Draw("samehist");
      (*this).histoBackground_->Draw("samehist");

      TLegend *legend2= new TLegend( 1-(*this).cMVAs_->GetRightMargin()-0.4, 1-(*this).cMVAs_->GetTopMargin()-0.15,
	        		     1-(*this).cMVAs_->GetRightMargin()-0.02,1-(*this).cMVAs_->GetTopMargin()-0.03);

      if (htype == CompareType) {
	// if overtraining check, load additional histograms                                                                                                                           
	TH1* sigOv = 0;
	TH1* bgdOv = 0;

	TString ovname = hname += "_Train";
	sigOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_S" ));
	bgdOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_B" ));

	// normalise both signal and background                                                                                                                                        

	if (sigOv == 0 || bgdOv == 0) {
	  std::cout << "+++ Problem in \"mvas.C\": overtraining check histograms do not exist" << std::endl; return;
	}
	
	(*this).NormalizeHists((TH1F*)sigOv,(TH1F*)bgdOv );

	Int_t col = (*this).histoSignal_->GetLineColor();
	sigOv->SetMarkerColor(col);
	sigOv->SetMarkerSize(0.7);
	sigOv->SetMarkerStyle(20);
	sigOv->SetLineWidth(1);
	sigOv->SetLineColor(col);
	sigOv->Draw("e1same");

	col = (*this).histoBackground_->GetLineColor();
	bgdOv->SetMarkerColor(col);
	bgdOv->SetMarkerSize(0.7);
	bgdOv->SetMarkerStyle(20);
	bgdOv->SetLineWidth(1);
	bgdOv->SetLineColor(col);
	bgdOv->Draw("e1same");

	ymax = TMath::Max( ymax, float(TMath::Max(sigOv->GetMaximum(),bgdOv->GetMaximum())*maxMult));
	frame->GetYaxis()->SetLimits(0,ymax);

        std::cout << "--- Found comparison histograms for overtraining check" << std::endl;  

        legend2->SetFillStyle(0);
        legend2->SetFillColor(0);
        legend2->SetTextSize (0.03);
        legend2->SetBorderSize(0);
        legend2->SetTextFont(42);
        legend2->SetLineColor(1);
        legend2->SetLineStyle(1);
        legend2->SetLineWidth(2);
 	legend2->AddEntry(sigOv,"Signal (training)","P");
 	legend2->AddEntry(bgdOv,"Background (training)","P");
 	legend2->SetMargin( 0.1 );
 	legend2->Draw("same");

	// perform K-S test                                                                                                                                                            
	std::cout << "--- Perform Kolmogorov-Smirnov tests" << std::endl;
	Double_t kolS = (*this).histoSignal_->KolmogorovTest(sigOv);
	Double_t kolB = (*this).histoBackground_->KolmogorovTest(bgdOv);
	std::cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << std::endl;

	TString probatext = Form("Kolmogorov-Smirnov test:");
        TString probatext2 = Form("signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
	TLatex* tt  = new TLatex(0.25,0.72,probatext);
	TLatex* tt2 = new TLatex(0.25,0.69,probatext2);
	tt->SetNDC(); 
        tt->SetTextSize(0.025); 
        tt->AppendPad();
	tt2->SetNDC(); 
        tt2->SetTextSize(0.025); 
        tt2->AppendPad();

	// chi2 test                                                                                                                                            
	std::cout << "--- Perform Chi2 tests" << std::endl;
	Double_t chi2S = (*this).histoSignal_->Chi2Test(sigOv,"WW CHI2/NDF" );
	Double_t chi2B = (*this).histoBackground_->Chi2Test(bgdOv,"WW CHI2/NDF" );
	std::cout << "--- Goodness of signal (background) consistency: " << chi2S << " (" << chi2B << ")" << std::endl;

	TString probatext3 = Form("#Chi^{2} test:");
        TString probatext4 = Form("signal (background) #Chi^{2}/ndf = %5.3g (%5.3g)", chi2S, chi2B );
	TLatex* tt3  = new TLatex(0.25,0.65,probatext3);
	TLatex* tt4  = new TLatex(0.25,0.62,probatext4);
	tt3->SetNDC(); 
        tt3->SetTextSize(0.025); 
        tt3->AppendPad();
	tt4->SetNDC(); 
        tt4->SetTextSize(0.025); 
        tt4->AppendPad();

      }

      // redraw axes                                                                                                                                                                    
      frame->Draw("sameaxis");

      // text for overflows                                                                                                                                                             
      Int_t    nbin = (*this).histoSignal_->GetNbinsX();
      Double_t dxu  = (*this).histoSignal_->GetBinWidth(0);
      Double_t dxo  = (*this).histoSignal_->GetBinWidth(nbin+1);
      TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%",
			     (*this).histoSignal_->GetBinContent(0)*dxu*100, (*this).histoBackground_->GetBinContent(0)*dxu*100,
			     (*this).histoSignal_->GetBinContent(nbin+1)*dxo*100, (*this).histoBackground_->GetBinContent(nbin+1)*dxo*100 );

      TText* t = new TText( 0.975, 0.115, uoflow );
      t->SetNDC();
      t->SetTextSize( 0.030 );
      t->SetTextAngle( 90 );
      t->AppendPad();

      // update canvas                                                                                                                                                                  
      TLatex latex;
      latex.SetNDC();
      latex.SetTextAlign(21); // align right                                                                                                                                 
      latex.SetTextSize(0.033);
      latex.DrawLatex(0.6,0.92,Form("CMS Preliminary Simulation,#sqrt{s} = 8 TeV, W+jets"));
      (*this).cMVAs_->Update();

      methodTitle.ReplaceAll(" ","_");
      methodTitle.ReplaceAll("(","_");
      methodTitle.ReplaceAll(")","_");
      methodTitle.ReplaceAll("[","_");
      methodTitle.ReplaceAll("]","_");
      methodTitle.ReplaceAll("{","_");
      methodTitle.ReplaceAll("}","_");

      if      (htype == MVAType)     (*this).PrintImage((*this).cMVAs_, std::string(Form("%s/mva_output_%s",outputPlotDirectory.c_str(),methodTitle.Data())));
      else if (htype == ProbaType)   (*this).PrintImage((*this).cMVAs_, std::string(Form("%s/probability_%s",outputPlotDirectory.c_str(),methodTitle.Data())));
      else if (htype == CompareType) (*this).PrintImage((*this).cMVAs_, std::string(Form("%s/overtraining_%s",outputPlotDirectory.c_str(),methodTitle.Data())));
      else                           (*this).PrintImage((*this).cMVAs_, std::string(Form("%s/rarity_%s",outputPlotDirectory.c_str(),methodTitle.Data())));
    
      (*this).mvas_index = (*this).mvas_index  +1 ;

      if(frame!=0)   delete frame;
      if(legend!=0)  delete legend;
      if(t!=0)       delete t;
      if(legend2!=0) delete legend2 ;
	       
    }

  }	       

}

void TMVAGlob::ReadHistograms (TFile* inputFile){
  
  if ((*this).fInfoList_!=NULL) {
    delete (*this).fInfoList_;
    (*this).fInfoList_ = NULL ;
  }
  
  // search for the right histograms in full list of keys                                                                                                                                 
  (*this).fInfoList_ = new std::vector<significanceBox*> ;
  TIter next(inputFile->GetListOfKeys());
  TKey *key(0);
  while( (key = (TKey*)next()) ) {

    if(!TString(key->GetName()).BeginsWith("Method_")) continue;
    if(!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory") ) continue;

    std::cout << "--- Found directory: " << ((TDirectory*)key->ReadObj())->GetName() << std::endl;
    TDirectory* mDir = (TDirectory*)key->ReadObj();
    
    TIter keyIt(mDir->GetListOfKeys());
    TKey *titkey;
    int maxLenTitle = 0 ;
    while((titkey = (TKey*)keyIt())) {
      significanceBox* significance_ = new significanceBox;
      if(!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory") ) continue;
      TDirectory* titDir = (TDirectory *)titkey->ReadObj();      
      (*this).GetMethodName(significance_->methodName_,key);
      (*this).GetMethodTitle(significance_->methodTitle_,titDir);
      if (significance_->methodTitle_.Length() > maxLenTitle) maxLenTitle = significance_->methodTitle_.Length();
      TString hname = "MVA_" + significance_->methodTitle_;

      std::cout << "--- Classifier: " << significance_->methodTitle_ << std::endl;
      
      if((*this).histoBackground_ == NULL ) significance_->Background_ = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));
      else significance_->Background_ = (*this).histoBackground_ ;

      if((*this).histoSignal_ ==NULL)  significance_->Signal_ = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
      else significance_->Signal_ = (*this).histoSignal_ ;

      if((*this).effBackground_ == NULL ) significance_->efficiencyBackground_ = dynamic_cast<TH1*>(titDir->Get( hname + "_effB" ));
      else significance_->efficiencyBackground_ = (*this).effBackground_ ;

      if((*this).effSignal_ ==NULL)  significance_->efficiencySignal_      = dynamic_cast<TH1*>(titDir->Get( hname + "_effS" ));
      else significance_->efficiencySignal_ =(*this).effSignal_ ;

     
      if (significance_->efficiencyBackground_ == NULL || significance_->efficiencySignal_ == NULL){ delete significance_; continue; }
      (*this).fInfoList_->push_back(significance_);
    }
  }
  
 return;
}

TString TMVAGlob::GetFormula(){

  TString f = fFormula_;

  f.ReplaceAll("epsilonS","x");
  f.ReplaceAll("S","x");
  f.ReplaceAll("epsilonB","y");
  f.ReplaceAll("B","y"); 

  return f;
}

TString TMVAGlob::GetLatexFormula(){
  TString f = fFormula_;

  if (signalType_ == true) 
   f.ReplaceAll("S","#epsilon_{S}");
  if (backgroundType_ == true) 
   f.ReplaceAll("B","#epsilon_{B}");   

  f.ReplaceAll("(","{");
  f.ReplaceAll(")","}");
  f.ReplaceAll("sqrt","#sqrt");

  return f;
}


void TMVAGlob::plotSignificance (TFile* inputFile, const int & iFile, SignificanceType stype, const double & numberSignalEvents, const double & numberBackgroundEvents,
  		                 const bool & UseSignalEfficiency, const bool & UseBackgroundEfficiency, const std::string & outputPlotDirectory){


  if(inputFile == 0) {std::cerr<<" empty file --> exit "<<std::endl; return ; }

  (*this).ReadHistograms(inputFile); 
  
  if(stype == 0)       (*this).SetFormula("S/B"); 
  else if(stype == 1)  (*this).SetFormula("S/sqrt(B)");  
  else if(stype == 2)  (*this).SetFormula("S/sqrt(S+B)");  
  else if(stype == 3)  (*this).SetFormula("2*(sqrt(S+B)-sqrt(B))");  
  else { std::cerr<<" Not known formula --> exit "<<std::endl; return; }

  (*this).SetSignalType(UseSignalEfficiency);
  (*this).SetBackgroundType(UseBackgroundEfficiency);

  TFormula significanceFormula("significanceFormula",(*this).GetFormula());
  
  std::vector<significanceBox*>::iterator itList =  (*this).fInfoList_->begin();

  TString cname = "Classifier";
  int maxLenTitle = 0 ;
  if (cname.Length() >  maxLenTitle)  maxLenTitle = cname.Length();
  TString str = Form( "%*s   (  signal, backgr.)  Optimal-cut  %s      NSig      NBkg   EffSig   EffBkg",
		      maxLenTitle, cname.Data(), GetFormulaString().Data() );
  std::cout << "--- " << std::setfill('=') << std::setw(str.Length()) << "" << std::setfill(' ') << std::endl;
  std::cout << "--- " << str << std::endl;
  std::cout << "--- " << std::setfill('-') << std::setw(str.Length()) << "" << std::setfill(' ') << std::endl;
  Double_t maxSig    = -1;
  Double_t maxSigErr = -1;
  
  for ( ;itList!=(*this).fInfoList_->end(); ++itList) {
    // Cycle on the signal efficiency bin    
    if(signalType_)
      (*itList)->significance_ = new TH1F(Form("significance_%s_file%d_stype%d",(*itList)->methodTitle_.Data(),iFile,stype),"",(*itList)->efficiencySignal_->GetNbinsX(),(*itList)->efficiencySignal_->GetBinLowEdge(1),(*itList)->efficiencySignal_->GetBinLowEdge((*itList)->efficiencySignal_->GetNbinsX()+1));
    else
      (*itList)->significance_ = new TH1F(Form("significance_eff_%s_%d_stype%d",(*itList)->methodTitle_.Data(),iFile,stype),"",(*itList)->efficiencySignal_->GetNbinsX(),(*itList)->efficiencySignal_->GetBinLowEdge(1),(*itList)->efficiencySignal_->GetBinLowEdge((*itList)->efficiencySignal_->GetNbinsX()+1));

    for (Int_t iBin=1; iBin<=(*itList)->efficiencySignal_->GetNbinsX(); iBin++) {

      // as a function of the type choose to 
      Float_t S = 0;
      if(signalType_) S = (*itList)->efficiencySignal_->GetBinContent(iBin) * numberSignalEvents;
      else S = (*itList)->efficiencySignal_->GetBinContent(iBin) ;      

      Float_t B = 0;
      if(backgroundType_) B = (*itList)->efficiencyBackground_->GetBinContent(iBin) * numberBackgroundEvents;
      else B = (*itList)->efficiencyBackground_->GetBinContent(iBin) ;

      Double_t significance = 0. ;
      Double_t threshold = 0. ;

      if(signalType_) threshold = 0.005;
      else threshold = 0.01 ;

      if(S > threshold && B > threshold){
 
       if(stype == 0 && B != 0 )        significance = significanceFormula.Eval(S,B);
       else if(stype == 1 && B != 0 )   significance = significanceFormula.Eval(S,B);
       else if(stype == 2 && S+B != 0 ) significance = significanceFormula.Eval(S,B);
       else significance = significanceFormula.Eval(S,B);

                  
       if (significance > maxSig) {
	maxSig    = significance;
	if ((*this).GetFormulaString() == "S/B") maxSigErr = sqrt(S/(B*B)+S*S/(B*B*B));	
 	else if ((*this).GetFormulaString() == "S/sqrt(B)") maxSigErr = significance * sqrt( 1./S + 1./(2.*B));	
 	else if ((*this).GetFormulaString() == "S/sqrt(S+B)") maxSigErr = sqrt(S*(TMath::Power(1-0.5/sqrt(S+B),2))*1/(S+B)+B*S*S/(4*(S+B)));	
 	else if ((*this).GetFormulaString() == "2*(sqrt(S+B)-sqrt(B))") maxSigErr = sqrt(S*TMath::Power(1/sqrt(S+B),2)+B*TMath::Power(1/sqrt(S+B)-1/sqrt(B),2));	
       }
      }
      (*itList)->significance_->SetBinContent(iBin,significance);
      (*itList)->significance_->SetBinError(iBin,maxSigErr);
    }

    Int_t maxbin = (*itList)->significance_->GetMaximumBin();
    
    (*itList)->significance_->Scale(1/(*itList)->significance_->GetMaximum());
    (*itList)->maxSig_    = maxSig ;
    (*itList)->maxSigErr_ = maxSigErr ;      
    (*itList)->maxbin_    = maxbin ;
    
    TString opt = Form( "%%%is:  (%%8.4g,%%8.4g)    %%9.4g   %%10.6g  %%8.7g  %%8.7g %%8.4g %%8.4g",maxLenTitle);     
    std::cout << "--- "<< Form( opt.Data(), (*itList)->methodTitle_.Data(), numberSignalEvents, numberBackgroundEvents, (*itList)->significance_->GetXaxis()->GetBinCenter(maxbin),maxSig,
    			       (*itList)->efficiencySignal_->GetBinContent(maxbin)*numberSignalEvents, (*itList)->efficiencyBackground_->GetBinContent( maxbin )*numberBackgroundEvents,
                               (*itList)->efficiencySignal_->GetBinContent(maxbin), (*itList)->efficiencyBackground_->GetBinContent(maxbin) ) <<std::endl;

    
  }

  
  /// Plot the results  
  std::vector<significanceBox*>::iterator itInfoList = (*this).fInfoList_->begin();
  for ( ; itInfoList!=(*this).fInfoList_->end(); ++itInfoList ){
 
  // create new canvas                                                                                                                                                                 
   if(signalType_)
     (*this).cSignificance_ = new TCanvas( Form("cSignificance_%s_file_%d_type_%d",(*itInfoList)->methodTitle_.Data(),iFile,stype),Form("Efficiencies Classifier : %s",(*itInfoList)->methodTitle_.Data()),180,52,550,550);
   else
     (*this).cSignificance_ = new TCanvas( Form("cSignificance_eff_%s_%d_type_%d",(*itInfoList)->methodTitle_.Data(),iFile,stype),Form("Efficiencies Classifier : %s",(*itInfoList)->methodTitle_.Data()),180,52,550,550);

   (*this).cSignificance_->cd();
   (*this).cSignificance_->SetTicks();
   (*this).cSignificance_->SetFillColor(0);
   (*this).cSignificance_->SetBorderMode(0);
   (*this).cSignificance_->SetBorderSize(2);
  
   (*this).cSignificance_->SetTickx(1);
   (*this).cSignificance_->SetTicky(1);
   (*this).cSignificance_->SetLeftMargin(0.1);
   (*this).cSignificance_->SetRightMargin(0.12);
   (*this).cSignificance_->SetBottomMargin(0.12);
   (*this).cSignificance_->SetFrameBorderMode(0);

   // and the signal purity and quality                                                                                                                                                 
   (*itInfoList)->efficiencySignal_->SetTitle("Efficiencies and Optimal Cut");
   if ((*itInfoList)->methodTitle_.Contains("Cuts")) 
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "Signal Efficiency (#epsilon_{sig})" );    
   else if ((*itInfoList)->methodTitle_.Contains("Likelihood"))
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "Likelihood output" );    
   else if ((*itInfoList)->methodTitle_.Contains("LD"))
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "Linear Discriminant output" );    
   else if ((*itInfoList)->methodTitle_.Contains("BDT") && !(*itInfoList)->methodTitle_.Contains("BDTG"))
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "BDT output" );    
   else if ((*itInfoList)->methodTitle_.Contains("MLP"))
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "MLP output" );    
   else if ((*itInfoList)->methodTitle_.Contains("PDEFoam"))
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "PDEFoam output" );    
   else if ((*itInfoList)->methodTitle_.Contains("SVM"))
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "SVM output" );    
   else if ((*itInfoList)->methodTitle_.Contains("Fisher"))
      (*itInfoList)->significance_->GetXaxis()->SetTitle( "Fisher output" );    
   else
       (*itInfoList)->significance_->GetXaxis()->SetTitle( (*itInfoList)->methodTitle_ + " output" );

   (*itInfoList)->significance_->GetYaxis()->SetTitle("Efficiency");

   (*itInfoList)->significance_->GetXaxis()->SetTitleSize(0.035);
   (*itInfoList)->significance_->GetXaxis()->SetLabelSize(0.035);
   (*itInfoList)->significance_->GetXaxis()->SetTitleOffset(1.05);
   (*itInfoList)->significance_->GetYaxis()->SetTitleSize(0.035);
   (*itInfoList)->significance_->GetYaxis()->SetLabelSize(0.035);
   (*itInfoList)->significance_->GetYaxis()->SetTitleOffset(1.05);

   (*itInfoList)->significance_->GetYaxis()->SetRangeUser(0.,(*itInfoList)->significance_->GetMaximum()*1.3);

   (*itInfoList)->efficiencySignal_->SetLineColor(kBlue);
   (*itInfoList)->efficiencySignal_->SetLineWidth(2);

   (*itInfoList)->efficiencyBackground_->SetLineColor(kRed);
   (*itInfoList)->efficiencyBackground_->SetLineWidth(2);

   (*itInfoList)->significance_->SetLineColor(kBlack);
   (*itInfoList)->significance_->SetLineWidth(2);

   (*itInfoList)->significance_->Draw("histl");
   (*itInfoList)->efficiencySignal_->Draw("samehistl");   
   (*itInfoList)->efficiencyBackground_->Draw("samehistl");
   (*itInfoList)->significance_->Draw("sameaxis");

   
   // Draw legend                                                                                                                                                                       
   TLegend *legend1= new TLegend( cSignificance_->GetLeftMargin()+0.05, 1-cSignificance_->GetTopMargin()-0.17, cSignificance_->GetLeftMargin()+0.25, 1-cSignificance_->GetTopMargin()-0.02);
   legend1->SetFillStyle(0);
   legend1->SetFillColor(0);
   legend1->SetTextFont(42);
   legend1->SetTextSize(0.032);
   legend1->SetBorderSize(0.);
   legend1->SetMargin(0.3);

   legend1->AddEntry((*itInfoList)->efficiencySignal_,"Signal efficiency","L");
   legend1->AddEntry((*itInfoList)->efficiencyBackground_,"Background efficiency","L");
   legend1->Draw("same");
   
   TLegend *legend2= new TLegend(cSignificance_->GetLeftMargin()+0.4,1-cSignificance_->GetTopMargin()-0.09,1-cSignificance_->GetRightMargin()-0.2,1-cSignificance_->GetTopMargin()-0.02);
   legend2->SetFillStyle(0);
   legend2->SetFillColor(0);
   legend2->SetBorderSize(0.);
   legend2->SetTextFont(42);
   legend2->SetTextSize(0.032);
   legend2->SetMargin(0.3);

   legend2->AddEntry((*itInfoList)->significance_,Form("Significance %s",(*this).GetLatexFormula().Data()),"L");
   legend2->Draw("same");

   // line to indicate maximum efficiency                                                                                                                                               
   TLine* effline = new TLine( (*itInfoList)->significance_->GetXaxis()->GetXmin(), 1, (*itInfoList)->significance_->GetXaxis()->GetXmax(), 1 );
   effline->SetLineWidth(3);
   effline->SetLineStyle(7);
   effline->SetLineColor(210);
   effline->Draw("same");

   (*itInfoList)->significance_->Draw("samehistl");

   (*this).cSignificance_->Update();

   TGaxis* rightAxis = new TGaxis(cSignificance_->GetUxmax(),cSignificance_->GetUymin(),cSignificance_->GetUxmax(),cSignificance_->GetUymax(),0,1.1*(*itInfoList)->significance_->GetMaximum(),510,"+L");

   rightAxis->SetLineColor (kBlack);
   rightAxis->SetLabelColor(kBlack);
   rightAxis->SetTitleColor(kBlack);

   rightAxis->SetLabelSize(0.035);
   rightAxis->SetTitleSize(0.035);
   rightAxis->SetTitleOffset(1.22);
   rightAxis->SetTextFont(42);
   rightAxis->SetLabelFont(42);
   rightAxis->SetTitle("Significance");
   rightAxis->Draw();
   
   (*this).cSignificance_->Update();

   TLatex latex;
   latex.SetNDC();
   latex.SetTextAlign(21); // align right                                                                                                                                                 
   latex.SetTextSize(0.033);
   latex.DrawLatex(0.56,0.92,Form("CMS Preliminary Simulation,#sqrt{s} = 8 TeV, W+jets"));
   latex.Delete();

   // print comments                                                                                                                                                                    
   TString name = Form("For %.4g signal and %.4g background",numberSignalEvents,numberBackgroundEvents);
   TPaveText* line1 = new TPaveText(0.22,0.21,0.5,0.26,"NDC");
   line1->SetBorderSize(0);
   line1->SetFillColor(0);
   line1->SetFillStyle(0);
   line1->AddText(name.Data());
   line1->SetTextSize(0.027);
   line1->SetTextFont(42);
   line1->Draw("same");
  
   TPaveText* line2 = new TPaveText(0.22,0.15,0.5,0.2,"NDC");
   name = Form("max significance at %0.4g, cut at %0.4g",(*itInfoList)->maxSig_,(*itInfoList)->significance_->GetXaxis()->GetBinCenter((*itInfoList)->maxbin_));
   line2->AddText(name);
   line2->SetBorderSize(0);
   line2->SetFillColor(0);
   line2->SetFillStyle(0);
   line2->SetTextSize(0.027);
   line2->SetTextFont(42);
   line2->Draw("same");

   TString baseName ;
   if(UseSignalEfficiency && UseBackgroundEfficiency){
     if((*this).thisMethod_ >= (*this).method_index) continue ;
     if ((*itInfoList)->methodTitle_.Contains("Cuts"))      
       baseName  = Form("mva_significance_eff_%s_file%d",(*itInfoList)->methodTitle_.Data(),iFile);
     else
       baseName  = Form("mva_significance_eff_%s",(*itInfoList)->methodTitle_.Data());

     if(stype == 0)      (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_S_over_B_file",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 1) (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_S_over_sqrtB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 2) (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_S_over_sqrtSB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 3) (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_pval",outputPlotDirectory.c_str(),baseName.Data())));

     if ((*itInfoList)->methodTitle_.Contains("Cuts")) (*this).thisMethod_++;

   }
   else if( !UseSignalEfficiency && !UseBackgroundEfficiency){
     if((*this).thisMethod_ > (*this).method_index) continue ;
     if ((*itInfoList)->methodTitle_.Contains("Cuts"))      
       baseName  = Form("mva_significance_%s_file%d",(*itInfoList)->methodTitle_.Data(),iFile);
     else
       baseName  = Form("mva_significance_%s",(*itInfoList)->methodTitle_.Data());

     if(stype == 0) (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_S_over_B",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 1) (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_S_over_sqrtB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 2) (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_S_over_sqrtSB",outputPlotDirectory.c_str(),baseName.Data())));
     else if(stype == 3) (*this).PrintImage((*this).cSignificance_, std::string(Form("%s/%s_pval",outputPlotDirectory.c_str(),baseName.Data())));

     if ((*itInfoList)->methodTitle_.Contains("Cuts")) (*this).thisMethod_++;

   }

   if((*this).thisMethod_ == (*this).method_index) (*this).thisMethod_ = 0;
   
   //   delete rightAxis ;
   delete cSignificance_;
   legend1->Delete() ;   
   legend2->Delete() ;
   line1->Delete() ;
   line2->Delete() ;
   effline->Delete() ;

  }
}

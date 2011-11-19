//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jul  2 02:30:39 2011 by ROOT version 5.27/06b
// from TTree ntu/ntu
// found on file: /data1/dimatteo/Calibration/Ntuples/Run2011A/WZAnalysisSingleXtal/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1.root
//////////////////////////////////////////////////////////

#ifndef FastCalibratorWeight_h
#define FastCalibratorWeight_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <vector>
#include "hChain.h"
#include "h2Chain.h"
#include <TGraphErrors.h>

#include <TLorentzVector.h>
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"


class FastCalibratorWeight {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
    Int_t           runId;
    Int_t           lumiId;
    Int_t           isW;
    Int_t           isZ;
    std::vector<float>   *ele1_recHit_E;
    std::vector<int>     *ele1_recHit_hashedIndex;
    std::vector<int>     *ele1_recHit_ietaORix;
    std::vector<int>     *ele1_recHit_iphiORiy;
  
    Float_t         ele1_scERaw;
    Float_t         ele1_scE;
    Float_t         ele1_es;
    Float_t         ele1_e3x3;
    Float_t         ele1_tkP;
    Float_t         ele1_fbrem;
    Float_t         ele1_EOverP;
    Int_t           ele1_isEB;
    Int_t           ele1_isEBEEGap;
    Int_t           ele1_isEBEtaGap;
    Int_t           ele1_isEBPhiGap;
    Int_t           ele1_isEEDeeGap;
    Int_t           ele1_isEERingGap;
    std::vector<float>   *ele2_recHit_E;
    std::vector<int>     *ele2_recHit_hashedIndex;
    std::vector<int>     *ele2_recHit_iphiORiy;
    std::vector<int>     *ele2_recHit_ietaORix;
    
    Float_t         ele2_scERaw;
    Float_t         ele2_scE;
    Float_t         ele2_es;
    Float_t         ele2_e3x3;
    Float_t         ele2_tkP;
    Float_t         ele2_fbrem;
    Float_t         ele2_EOverP;
    Int_t           ele2_isEB;
    Int_t           ele2_isEBEEGap;
    Int_t           ele2_isEBEtaGap;
    Int_t           ele2_isEBPhiGap;
    Int_t           ele2_isEEDeeGap;
    Int_t           ele2_isEERingGap;

   // List of branches
    TBranch        *b_runId;   //!
    TBranch        *b_lumiId;   //!
    TBranch        *b_isW;   //!
    TBranch        *b_isZ;   //!
    TBranch        *b_ele1_recHit_E;   //!
    TBranch        *b_ele1_recHit_hashedIndex;
    TBranch        *b_ele1_recHit_iphiORiy;
    TBranch        *b_ele1_recHit_ietaORix;
       //!
    TBranch        *b_ele1_scERaw;   //!
    TBranch        *b_ele1_scE;   //!
    TBranch        *b_ele1_es;   //!
    TBranch        *b_ele1_e3x3;   //!
    TBranch        *b_ele1_tkP;   //!
    TBranch        *b_ele1_fbrem;   //!
    TBranch        *b_ele1_EOverP;   //!
    TBranch        *b_ele1_isEB;   //!
    TBranch        *b_ele1_isEBEEGap;   //!
    TBranch        *b_ele1_isEBEtaGap;   //!
    TBranch        *b_ele1_isEBPhiGap;   //!
    TBranch        *b_ele1_isEEDeeGap;   //!
    TBranch        *b_ele1_isEERingGap;   //!
    TBranch        *b_ele2_recHit_E;   //!
    TBranch        *b_ele2_recHit_hashedIndex;
    TBranch        *b_ele2_recHit_iphiORiy;
    TBranch        *b_ele2_recHit_ietaORix;   //!
    TBranch        *b_ele2_scERaw;   //!
    TBranch        *b_ele2_scE;   //!
    TBranch        *b_ele2_es;   //!
    TBranch        *b_ele2_e3x3;   //!
    TBranch        *b_ele2_tkP;   //!
    TBranch        *b_ele2_fbrem;   //!
    TBranch        *b_ele2_EOverP;   //!
    TBranch        *b_ele2_isEB;   //!
    TBranch        *b_ele2_isEBEEGap;   //!
    TBranch        *b_ele2_isEBEtaGap;   //!
    TBranch        *b_ele2_isEBPhiGap;   //!
    TBranch        *b_ele2_isEEDeeGap;   //!
    TBranch        *b_ele2_isEERingGap;   //!

    FastCalibratorWeight(TTree *tree=0);
    virtual ~FastCalibratorWeight();
    virtual void     bookHistos(int);
    virtual void     saveHistos(TFile *f1);
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop(int, int, int, int, int);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual void     printOnTxt(std::string outputTxtFile);
    virtual void     BuildEoPeta_ele1(int,int,int,int,std::vector<float>);
    virtual void     BuildEoPeta_ele2(int,int,int,int,std::vector<float>); 
    virtual void     saveEoPeta(TFile *f1);

    hChain     *hC_EoP_eta_ele1;
    hChain     *hC_EoP_eta_ele2;

    hChain     *hC_IntercalibValues;
    hChain     *hC_EoP;
    hChain     *hC_PullFromScalib;
    h2Chain    *hC_scale_EB;
    TH1F       *h_Occupancy_hashedIndex;
    TH2F       *h_occupancy;
    TProfile   *p_IntercalibValues_iEta;
    TH2F       *h_scalib_EB;
    TH2F       *h_scale_EB;
    TH2F       *h_scale_EB_meanOnPhi;
    TH1F       *h_scale_EB_hashedIndex;
    TH1F       *h_IntercalibSpread_iEta;
    TH1F       *h_IntercalibValues_test;
    TH1F       *h_Init_IntercalibValues;

    std::vector<int>   IetaValues;
    std::vector<int>   IphiValues;
    std::vector<float> ICValues;
    std::vector<float> meanICforPhiRingValues;
    
    TGraphErrors *g_ICmeanVsLoop;
    TGraphErrors *g_ICrmsVsLoop;

};

#endif

#include "FastCalibratorEE.h"
#include "GetHashedIndexEE.h"
#include "EERings.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <TRandom.h>



FastCalibratorEE::FastCalibratorEE(TTree *tree,TString outEPDistribution):
outEPDistribution_p(outEPDistribution)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data1/dimatteo/Calibration/Ntuples/Run2011A/WZAnalysisSingleXtal/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1.root");
    if (!f) {
      f = new TFile("/data1/dimatteo/Calibration/Ntuples/Run2011A/WZAnalysisSingleXtal/WZAnalysis_SingleElectron_Run2011A-WElectron-May10ReReco-v1.root");
    }
    tree = (TTree*)gDirectory->Get("ntu");

  }
   

  SumIC_Ring_EEP.assign(40,0);
  SumIC_Ring_EEM.assign(40,0);
  Sumxtal_Ring_EEP.assign(40,0);
  Sumxtal_Ring_EEM.assign(40,0);
  
  Init(tree);
}

FastCalibratorEE::~FastCalibratorEE()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t FastCalibratorEE::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t FastCalibratorEE::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void FastCalibratorEE::Init(TTree *tree)
{
   // Set object pointer
  ele1_recHit_E = 0;
  ele1_recHit_hashedIndex = 0;
  ele1_recHit_iphiORiy = 0;
  ele1_recHit_ietaORix =0 ;
  ele1_recHit_flag =0 ;

  ele2_recHit_E = 0;
  ele2_recHit_hashedIndex = 0;
  ele2_recHit_iphiORiy = 0;
  ele2_recHit_ietaORix = 0;
  ele2_recHit_flag =0 ;
   // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("runId", &runId, &b_runId);
  fChain->SetBranchAddress("lumiId", &lumiId, &b_lumiId);
  fChain->SetBranchAddress("isW", &isW, &b_isW);
  fChain->SetBranchAddress("isZ", &isZ, &b_isZ);
  fChain->SetBranchAddress("ele1_recHit_E", &ele1_recHit_E, &b_ele1_recHit_E);
  fChain->SetBranchAddress("ele1_recHit_hashedIndex", &ele1_recHit_hashedIndex, &b_ele1_recHit_hashedIndex);
  fChain->SetBranchAddress("ele1_recHit_iphiORiy", &ele1_recHit_iphiORiy, &b_ele1_recHit_iphiORiy);
  fChain->SetBranchAddress("ele1_recHit_ietaORix", &ele1_recHit_ietaORix, &b_ele1_recHit_ietaORix);
  fChain->SetBranchAddress("ele1_recHit_flag", &ele1_recHit_flag, &b_ele1_recHit_flag);
 
  fChain->SetBranchAddress("ele1_scERaw", &ele1_scERaw, &b_ele1_scERaw);
  fChain->SetBranchAddress("ele1_scE", &ele1_scE, &b_ele1_scE);
  fChain->SetBranchAddress("ele1_es", &ele1_es, &b_ele1_es);
  fChain->SetBranchAddress("ele1_e3x3", &ele1_e3x3, &b_ele1_e3x3);
  fChain->SetBranchAddress("ele1_tkP", &ele1_tkP, &b_ele1_tkP);
  fChain->SetBranchAddress("ele1_fbrem", &ele1_fbrem, &b_ele1_fbrem);
  fChain->SetBranchAddress("ele1_EOverP", &ele1_EOverP, &b_ele1_EOverP);
  fChain->SetBranchAddress("ele1_isEB", &ele1_isEB, &b_ele1_isEB);
  fChain->SetBranchAddress("ele1_isEBEEGap", &ele1_isEBEEGap, &b_ele1_isEBEEGap);
  fChain->SetBranchAddress("ele1_isEBEtaGap", &ele1_isEBEtaGap, &b_ele1_isEBEtaGap);
  fChain->SetBranchAddress("ele1_isEBPhiGap", &ele1_isEBPhiGap, &b_ele1_isEBPhiGap);
  fChain->SetBranchAddress("ele1_isEEDeeGap", &ele1_isEEDeeGap, &b_ele1_isEEDeeGap);
  fChain->SetBranchAddress("ele1_isEERingGap", &ele1_isEERingGap, &b_ele1_isEERingGap);
  fChain->SetBranchAddress("ele2_recHit_E", &ele2_recHit_E, &b_ele2_recHit_E);
  fChain->SetBranchAddress("ele2_recHit_hashedIndex", &ele2_recHit_hashedIndex, &b_ele2_recHit_hashedIndex);
  fChain->SetBranchAddress("ele2_recHit_iphiORiy", &ele2_recHit_iphiORiy, &b_ele2_recHit_iphiORiy);
  fChain->SetBranchAddress("ele2_recHit_ietaORix", &ele2_recHit_ietaORix, &b_ele2_recHit_ietaORix);
  fChain->SetBranchAddress("ele2_recHit_flag", &ele2_recHit_flag, &b_ele2_recHit_flag);
 
  fChain->SetBranchAddress("ele2_scERaw", &ele2_scERaw, &b_ele2_scERaw);
  fChain->SetBranchAddress("ele2_scE", &ele2_scE, &b_ele2_scE);
  fChain->SetBranchAddress("ele2_es", &ele2_es, &b_ele2_es);
  fChain->SetBranchAddress("ele2_e3x3", &ele2_e3x3, &b_ele2_e3x3);
  fChain->SetBranchAddress("ele2_tkP", &ele2_tkP, &b_ele2_tkP);
  fChain->SetBranchAddress("ele2_fbrem", &ele2_fbrem, &b_ele2_fbrem);
  fChain->SetBranchAddress("ele2_EOverP", &ele2_EOverP, &b_ele2_EOverP);
  fChain->SetBranchAddress("ele2_isEB", &ele2_isEB, &b_ele2_isEB);
  fChain->SetBranchAddress("ele2_isEBEEGap", &ele2_isEBEEGap, &b_ele2_isEBEEGap);
  fChain->SetBranchAddress("ele2_isEBEtaGap", &ele2_isEBEtaGap, &b_ele2_isEBEtaGap);
  fChain->SetBranchAddress("ele2_isEBPhiGap", &ele2_isEBPhiGap, &b_ele2_isEBPhiGap);
  fChain->SetBranchAddress("ele2_isEEDeeGap", &ele2_isEEDeeGap, &b_ele2_isEEDeeGap);
  fChain->SetBranchAddress("ele2_isEERingGap", &ele2_isEERingGap, &b_ele2_isEERingGap);
  Notify();
}

Bool_t FastCalibratorEE::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  return kTRUE;
}

void FastCalibratorEE::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t FastCalibratorEE::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  return 1;
}

void FastCalibratorEE::bookHistos(int nLoops)
{

  //service histos
  h_scale_hashedIndex_EE = new TH1F ("h_scale_hashedIndex_EE","h_scale_hashedIndex_EE", kEEhalf*2, 0, kEEhalf*2 - 1 );
  h_occupancy_hashedIndex_EE = new TH1F ("h_occupancy_hashedIndex_EE","h_occupancy_hashedIndex_EE", kEEhalf*2, 0, kEEhalf*2 - 1 );
  hC_EoP = new hChain ("EoP", "EoP", 500,0.2,3.0, nLoops);
 
  //EE+

  hC_IntercalibValues_EEP = new hChain ("IntercalibValues_EEP", "IntercalibValues_EEP", 400,0.2,1.9, nLoops);
  hC_PullFromScalib_EEP = new hChain ("hC_PullFromScalib_EEP", "hC_PullFromScalib_EEP", 2000,-0.5,0.5, nLoops);
  hC_scale_EEP = new h2Chain("hC_scale_EEP", "hC_scale_EEP", 100,1, 101, 100, 1, 101, nLoops );
  
  h_scale_EEP = new TH2F("h_scale_EEP", "h_scale_EEP", 100,1, 101, 100, 1, 101 );
  h_occupancy_EEP = new TH2F("h_occupancy_EEP", "h_occupancy_EEP", 100,1, 101, 100, 1, 101 );
  h_scalib_EEP = new TH2F("h_scalib_EEP", "h_scalib_EEP", 100,1, 101, 100, 1, 101);
  h_map_Dead_Channels_EEP = new TH2F("h_map_Dead_Channels_EEP","h_map_Dead_Channels_EEP",100,1,101,100,1,101);
  h_scale_meanOnring_EEP = new TH2F ("h_scale_meanOnring_EEP", "h_scale_meanOnring_EEP",  100,1, 101, 100, 1, 101);
  
  g_ICmeanVsLoop_EEP = new TGraphErrors();
  g_ICmeanVsLoop_EEP -> SetName("g_ICmeanVsLoop_EEP");
  g_ICmeanVsLoop_EEP -> SetTitle("g_ICmeanVsLoop_EEP");
  
  g_ICrmsVsLoop_EEP = new TGraphErrors();
  g_ICrmsVsLoop_EEP -> SetName("g_ICrmsVsLoop_EEP");
  g_ICrmsVsLoop_EEP -> SetTitle("g_ICrmsVsLoop_EEP");

 
  
  //EE-
  hC_IntercalibValues_EEM = new hChain ("IntercalibValues_EEM", "IntercalibValues_EEM", 400,0.2,1.9, nLoops);
  hC_PullFromScalib_EEM = new hChain ("hC_PullFromScalib_EEM", "hC_PullFromScalib_EEM", 2000,-0.5,0.5, nLoops);
  hC_scale_EEM = new h2Chain("hC_scale_EEM", "hC_scale_EEM", 100,1, 101, 100, 1, 101, nLoops );
  
  h_scale_EEM = new TH2F("h_scale_EEM", "h_scale_EEM", 100,1, 101, 100, 1, 101 );
  h_occupancy_EEM = new TH2F("h_occupancy_EEM", "h_occupancy_EEM", 100,1, 101, 100, 1, 101 );
  h_scalib_EEM = new TH2F("h_scalib_EEM", "h_scalib_EEM", 100,1, 101, 100, 1, 101);
  h_map_Dead_Channels_EEM = new TH2F("h_map_Dead_Channels_EEM","h_map_Dead_Channels_EEM",100,1,101,100,1,101);
  h_scale_meanOnring_EEM = new TH2F ("h_scale_meanOnring_EEM", "h_scale_meanOnring_EEM",  100,1, 101, 100, 1, 101);

  g_ICmeanVsLoop_EEM = new TGraphErrors();
  g_ICmeanVsLoop_EEM -> SetName("g_ICmeanVsLoop_EEM");
  g_ICmeanVsLoop_EEM -> SetTitle("g_ICmeanVsLoop_EEM");
  
  g_ICrmsVsLoop_EEM = new TGraphErrors();
  g_ICrmsVsLoop_EEM -> SetName("g_ICrmsVsLoop_EEM");
  g_ICrmsVsLoop_EEM -> SetTitle("g_ICrmsVsLoop_EEM");


  return;
}




///===== Build E/p for electron 1 and 2

void FastCalibratorEE::BuildEoPeta_ele(int iLoop, int nentries , int useW, int useZ, std::vector<float> theScalibration,bool       isSaveEPDistribution, bool isR9selection)
{
  if(iLoop ==0)  
  {
   TString name = Form ("hC_EoP_eta_%d",iLoop);
   hC_EoP_ir_ele = new hChain (name,name, 250,0.1,3.0,41);
  }
  else{
          hC_EoP_ir_ele -> Reset();
          TString name = Form ("hC_EoP_eta_%d",iLoop);
          hC_EoP_ir_ele = new hChain (name,name, 250,0.1,3.0,41);
      }

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   
   nbytes += nb;
   if (!(jentry%1000000))std::cerr<<"building E/p distribution ----> "<<jentry<<" vs "<<nentries<<std::endl;

   float pIn, pSub, FdiEta;

   ///=== electron tight W or Z only Endcaps
   if ( ele1_isEB == 0 && (( useW == 1 && isW == 1 ) ||  ( useZ== 1 && isZ == 1 ))) {

    FdiEta = ele1_scE/(ele1_scERaw+ele1_es);
   
    float thisE = 0;
    int   iseed = 0 ;
    int seed_hashedIndex = 0;
    float E_seed = 0;
    float thisE3x3 = 0;
    // Cycle on the all the recHits of the Event: to get the old IC and the corrected SC energy
    for (unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
          

            if(ele1_recHit_E -> at(iRecHit) > E_seed && ele1_recHit_flag -> at(iRecHit) < 4 )
            {
              seed_hashedIndex=ele1_recHit_hashedIndex -> at(iRecHit);
              iseed=iRecHit;
              E_seed=ele1_recHit_E -> at(iRecHit);

            }
    
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
            
            if(ele1_recHit_flag -> at(iRecHit) < 4)
            thisE += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
          
     }

    for (unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
         
            if(fabs(ele1_recHit_ietaORix->at(iRecHit)-ele1_recHit_ietaORix->at(iseed))<=1 && 
               fabs(ele1_recHit_iphiORiy->at(iRecHit)-ele1_recHit_iphiORiy->at(iseed))<=1 &&
               ele1_recHit_flag -> at(iRecHit) < 4)
              thisE3x3+=theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
           }

          
     int ix_seed = GetIxFromHashedIndex(seed_hashedIndex);
     int iy_seed = GetIyFromHashedIndex(seed_hashedIndex);
     int iz_seed = GetZsideFromHashedIndex(seed_hashedIndex);
     int ir_seed = EERings(ix_seed,iy_seed,iz_seed);
 
     pSub = 0.; //NOTALEO : test dummy
     pIn = ele1_tkP;
     bool skipElectron = false;
     if ( fabs(thisE3x3/thisE) < 0.9 && isR9selection == true) skipElectron = true;
     if(!skipElectron)    hC_EoP_ir_ele -> Fill(ir_seed,thisE/(ele1_tkP-ele1_es));
     
  
  }
  ///=== Second medium electron from Z only Endcaps
  if ( ele2_isEB == 0 && (( useW == 1 && isW == 1 ) || ( useZ == 1 && isZ == 1 )) ){

    FdiEta = ele2_scE/(ele2_scERaw+ele2_es);
      // Electron energy
    float thisE = 0;
    int   iseed = 0 ;
    int seed_hashedIndex = 0;
    float E_seed = 0;
    float thisE3x3 = 0;
  
    // Cycle on the all the recHits of the Event: to get the old IC and the corrected SC energy
    for (unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);

            if(ele2_recHit_E -> at(iRecHit) > E_seed && ele2_recHit_flag -> at(iRecHit) < 4 )
            {
              seed_hashedIndex=ele2_recHit_hashedIndex -> at(iRecHit);
              iseed=iRecHit;
              E_seed=ele2_recHit_E -> at(iRecHit);

            }
    
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
            
            if(ele2_recHit_flag -> at(iRecHit) < 4)
            thisE += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;
             
     }

     for (unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
         
            if(fabs(ele2_recHit_ietaORix->at(iRecHit)-ele2_recHit_ietaORix->at(iseed))<=1 && 
               fabs(ele2_recHit_iphiORiy->at(iRecHit)-ele2_recHit_iphiORiy->at(iseed))<=1 &&
               ele2_recHit_flag -> at(iRecHit) < 4)
              thisE3x3+=theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;
           }

  
     int ix_seed = GetIxFromHashedIndex(seed_hashedIndex);
     int iy_seed = GetIyFromHashedIndex(seed_hashedIndex);
     int iz_seed = GetZsideFromHashedIndex(seed_hashedIndex);
     int ir_seed = EERings(ix_seed,iy_seed,iz_seed);
 
     
     pSub = 0.; //NOTALEO : test dummy
     pIn = ele2_tkP;
 
     bool skipElectron = false;
     if ( fabs(thisE3x3/thisE) < 0.9 && isR9selection==true) skipElectron = true;
     if(!skipElectron) hC_EoP_ir_ele -> Fill(ir_seed,thisE/(ele2_tkP-ele2_es));
  
  }
  
  
 }
 
 for(unsigned int ir=0 ; ir < hC_EoP_ir_ele->Size() ; ir++)
 {
     hC_EoP_ir_ele->Normalize(ir);
 }
 
 if(isSaveEPDistribution == true) 
 {
   TFile *f2 = new TFile(outEPDistribution_p.Data(),"UPDATE");
   saveEoPeta(f2);
 }

}



void FastCalibratorEE::Loop(int nentries, int useZ, int useW, int splitStat, int nLoops, bool isMiscalib,bool isSaveEPDistribution,
                                bool isEPselection,bool isR9selection)
{
   if (fChain == 0) return;
   
   // Define the number of crystal you want to calibrate
   int m_regions = kEEhalf;
   
   std::cout << "m_regions " << m_regions << std::endl;
  
     // build up scalibration map
   std::vector<float> theScalibration(m_regions*2, 0.);
   TRandom genRand;
   for ( int iIndex = 0; iIndex < m_regions*2; iIndex++ ){
     bool isDeadXtal = false ;
     if(DeadXtal_HashedIndex.at(0)!=-9999) isDeadXtal = CheckDeadXtal(GetIxFromHashedIndex(iIndex), GetIyFromHashedIndex(iIndex),GetZsideFromHashedIndex(iIndex));
     if(isDeadXtal == true ) {
     theScalibration[iIndex]=0;
    
     if(GetZsideFromHashedIndex(iIndex)>0)
     h_map_Dead_Channels_EEP->Fill(GetIxFromHashedIndex(iIndex),GetIyFromHashedIndex(iIndex));
     else h_map_Dead_Channels_EEM->Fill(GetIxFromHashedIndex(iIndex),GetIyFromHashedIndex(iIndex));
     }
     else{
     
         if(isMiscalib==true) theScalibration[iIndex] = genRand.Gaus(1.,0.05);
         if(isMiscalib == false) theScalibration[iIndex] = 1.;
         if(GetZsideFromHashedIndex(iIndex)>0)
         h_scalib_EEP -> Fill ( GetIxFromHashedIndex(iIndex), GetIyFromHashedIndex(iIndex), theScalibration[iIndex] );
         else  h_scalib_EEM-> Fill ( GetIxFromHashedIndex(iIndex), GetIyFromHashedIndex(iIndex), theScalibration[iIndex] );

     }
   }
  
  
   /// ----------------- Calibration Loops -----------------------------//
   for ( int iLoop = 0; iLoop < nLoops; iLoop++ ) {
    // loop over events
    std::cout << "Starting iteration " << iLoop + 1 << std::endl;
  
    std::vector<float> theNumerator_EEP(m_regions*2+1, 0.);
    std::vector<float> theDenominator_EEP(m_regions*2+1, 0.);
    std::vector<float> theNumerator_EEM(m_regions+1, 0.);
    std::vector<float> theDenominator_EEM(m_regions+1, 0.);

    BuildEoPeta_ele(iLoop,nentries,useW,useZ,theScalibration,isSaveEPDistribution,isR9selection); ///==== build E/p distribution ele 1 and 2
    
    // loop over events
    std::cout << "Number of analyzed events = " << nentries << std::endl;
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
        if (!(jentry%10000))std::cerr<<jentry;
        if (!(jentry%1000)) std::cerr<<".";
      
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   
        nbytes += nb;
              
        // Event cut:
        //if ( runId != ilrunchevuoi ) continue;
  
        // Electron quantities used for the calibration
        float pIn, pSub, FdiEta;
      
        std::map<int,double> map;
        bool skipElectron=false;
        // Only W only Endcap
        if ( ele1_isEB == 0 && (( useW == 1 && isW == 1 ) || ( useZ == 1 && isZ == 1 )) ) {
                  
          // SCL energy containment correction
          FdiEta = ele1_scE/(ele1_scERaw+ele1_es);
          // Electron energy
          float thisE = 0;
          float thisE3x3 =0 ;
          int iseed = 0 ;
          int seed_hashedIndex = 0 ;
          float E_seed = 0;

         
          // Cycle on the all the recHits of the Event: to get the old IC and the corrected SC energy
          for (unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
  
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
            
            if(ele1_recHit_flag -> at(iRecHit) < 4)
            thisE += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
     
            if(ele1_recHit_E -> at(iRecHit) > E_seed && ele1_recHit_flag -> at(iRecHit) < 4 )
             {
              E_seed=ele1_recHit_E -> at(iRecHit);
              iseed=iRecHit;
              seed_hashedIndex=ele1_recHit_hashedIndex -> at(iRecHit);
 
             }
          
          }
          
          // Cycle on the all the recHits of the Event
          for (unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
    
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
            
            if(fabs(ele1_recHit_ietaORix->at(iRecHit)-ele1_recHit_ietaORix->at(iseed))<=1 && 
               fabs(ele1_recHit_iphiORiy->at(iRecHit)-ele1_recHit_iphiORiy->at(iseed))<=1 &&
                ele1_recHit_flag -> at(iRecHit) < 4)
              thisE3x3+=theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
                 
              }
            
          pSub = 0.; //NOTALEO : test dummy
          pIn = ele1_tkP;
                
          // find the zside
         int thisCaliBlock = -1;
         if (GetZsideFromHashedIndex(ele1_recHit_hashedIndex -> at(iseed)) < 0) thisCaliBlock = 0;
         else thisCaliBlock = 1;
 
         int ix_seed = GetIxFromHashedIndex(seed_hashedIndex);
         int iy_seed = GetIyFromHashedIndex(seed_hashedIndex);
         int iz_seed = GetZsideFromHashedIndex(seed_hashedIndex);
         int ir_seed = EERings(ix_seed,iy_seed,iz_seed);
      
         TH1F* EoPHisto = hC_EoP_ir_ele->GetHisto(ir_seed);
       
         if ( fabs(thisE/(ele1_tkP-ele1_es) - 1) > 0.7 && isEPselection==true) skipElectron = true;
         if ( fabs(thisE3x3/thisE) < 0.9 && isR9selection==true) skipElectron = true;
         if ( thisE/(ele1_tkP-ele1_es) < EoPHisto->GetXaxis()->GetXmin() || thisE/(ele1_tkP-ele1_es) > EoPHisto->GetXaxis()->GetXmax()) skipElectron=true;
 
         if ( !skipElectron ) {
                  
          for ( unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
           
           if(ele1_recHit_flag -> at(iRecHit) >= 4) continue;
         
           int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
           float thisIC = 1.;
         
           if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
  
               // Fill the occupancy map JUST for the first Loop
           if ( iLoop == 0 ) {
                     h_occupancy_hashedIndex_EE -> Fill(thisIndex);
                     if ( GetZsideFromHashedIndex(thisIndex) < 0 )
                       h_occupancy_EEM -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
                     else h_occupancy_EEP -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
                          
             }

            
            if ( splitStat == 0) {
             
                if(thisCaliBlock == 0) {
                
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele1_es));
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele1_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
 
                
                }
                
                if(thisCaliBlock == 1) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele1_es));
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele1_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
              }
             
              }
             
             
          // use evens    
          if ( splitStat == 1 && jentry%2 == 0 ) {
                  
                if(thisCaliBlock == 0) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele1_es));
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele1_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
               }
                
                if(thisCaliBlock == 1) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele1_es));
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele1_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
              }
            }
             
           // use odd    
           if ( splitStat == 1 && jentry%2 != 0 ) {
                  
                if(thisCaliBlock == 0) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele1_es));
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele1_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);}
                
                if(thisCaliBlock == 1) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele1_es));
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele1_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
              }
            }
            }
          }
          //Fill EoP
           if (thisCaliBlock != -1) hC_EoP -> Fill(iLoop, thisE/(pIn-ele1_es));
        
        }  
         skipElectron = false;     
        /// Fill the map with the ele (if any) from the second Z leg
        // Only Z only Barrel
        if ( ele2_isEB == 0 && ( useZ == 1 && isZ == 1 ) ) {
          
          // SCL energy containment correction
          FdiEta = ele2_scE/(ele2_scERaw+ele2_es);
          // Electron energy
          float thisE = 0;
          float thisE3x3 =0 ;
          int iseed = 0 ;
          int seed_hashedIndex = 0;
          float E_seed = 0;

         
          // Cycle on the all the recHits of the Event: to get the old IC and the corrected SC energy
          for (unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
            
            if( ele2_recHit_flag -> at(iRecHit) < 4 )
            thisE += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;

              
            if(ele2_recHit_E -> at(iRecHit) > E_seed && ele2_recHit_flag -> at(iRecHit) < 4)
             {
              E_seed=ele2_recHit_E -> at(iRecHit);
              iseed=iRecHit;
              seed_hashedIndex=ele2_recHit_hashedIndex -> at(iRecHit);
 
             }
          
          }
          
          // Cycle on the all the recHits of the Event
          for (unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);
            
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
            
            if(fabs(ele2_recHit_ietaORix->at(iRecHit)-ele2_recHit_ietaORix->at(iseed))<=1 && 
               fabs(ele2_recHit_iphiORiy->at(iRecHit)-ele2_recHit_iphiORiy->at(iseed))<=1 &&
               ele2_recHit_flag -> at(iRecHit) < 4)
              thisE3x3+=theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;
                 
              }
            
          pSub = 0.; //NOTALEO : test dummy
          pIn = ele2_tkP;
                
          // find the zside
          int thisCaliBlock = -1;
          if (GetZsideFromHashedIndex(ele2_recHit_hashedIndex -> at(iseed)) < 0) thisCaliBlock = 0;
          else thisCaliBlock = 1;
 
          int ix_seed = GetIxFromHashedIndex(seed_hashedIndex);
          int iy_seed = GetIyFromHashedIndex(seed_hashedIndex);
          int iz_seed = GetZsideFromHashedIndex(seed_hashedIndex);
          int ir_seed = EERings(ix_seed,iy_seed,iz_seed);
 
          TH1F* EoPHisto = hC_EoP_ir_ele->GetHisto(ir_seed);
          
          if ( fabs(thisE/(ele2_tkP-ele2_es) - 1) > 0.7 && isEPselection==true) skipElectron = true;
          if ( fabs(thisE3x3/thisE) < 0.9 && isR9selection==true) skipElectron = true;
          if ( thisE/(ele2_tkP-ele2_es) < EoPHisto->GetXaxis()->GetXmin() || thisE/(ele2_tkP-ele2_es) > EoPHisto->GetXaxis()->GetXmax()) skipElectron=true;
 
         if ( !skipElectron ) {
                  
          for ( unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
           
           if(ele2_recHit_flag -> at(iRecHit) >= 4) continue;
         
           int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);
           float thisIC = 1.;
              
           if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
  
               // Fill the occupancy map JUST for the first Loop
           if ( iLoop == 0 ) {
                     h_occupancy_hashedIndex_EE -> Fill(thisIndex);
                     if ( GetZsideFromHashedIndex(thisIndex) < 0 ) 
                       h_occupancy_EEM -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
                     else h_occupancy_EEP -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
             }

            
            if ( splitStat == 0) {
             
                if(thisCaliBlock == 0) {
                
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele2_es));
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele2_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
 
                
                }
                
                if(thisCaliBlock == 1) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele2_es));
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele2_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
              }
             
              }
             
             
          // use evens    
          if ( splitStat == 1 && jentry%2 == 0 ) {
                  
                if(thisCaliBlock == 0) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele2_es));
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele2_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
               }
                
                if(thisCaliBlock == 1) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele2_es));
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele2_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
              }
            }
             
           // use odd    
           if ( splitStat == 1 && jentry%2 != 0 ) {
                  
                if(thisCaliBlock == 0) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele2_es));
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele2_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);}
                
                if(thisCaliBlock == 1) {
                int EoPbin = EoPHisto->FindBin(thisE/(pIn-pSub-ele2_es));
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub-ele1_es)/thisE*EoPHisto->GetBinContent(EoPbin);
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*EoPHisto->GetBinContent(EoPbin);
              }
            }
            }
          }
          //Fill EoP
           if (thisCaliBlock != -1) hC_EoP -> Fill(iLoop, thisE/(pIn-ele2_es));
          
        }
   }  
    // Cycle on the events
       
    //Now solve everything
    std::cout << ">>>>> [L3][endOfLoop] entering..." << std::endl;
  
    TH1F auxiliary_IC_EEM("auxiliary_IC_EEM","auxiliary_IC_EEM",50,0.2,1.9);
    TH1F auxiliary_IC_EEP("auxiliary_IC_EEP","auxiliary_IC_EEP",50,0.2,1.9);

    //Fill the histo of IntercalibValues before the solve
    for ( int iIndex = 0; iIndex < kEEhalf*2; iIndex++ ){

      if ( h_occupancy_hashedIndex_EE -> GetBinContent(iIndex+1) > 0 ){
   
        int thisCaliBlock = -1;
        if (GetZsideFromHashedIndex(iIndex) < 0) thisCaliBlock = 0;
        else thisCaliBlock = 1;
        
        float thisIntercalibConstant =1. ;

        if(thisCaliBlock == 0 && theDenominator_EEM[iIndex] != 0.){ 
            thisIntercalibConstant = theNumerator_EEM[iIndex]/theDenominator_EEM[iIndex];}
        else{
             if(thisCaliBlock == 1 && theDenominator_EEP[iIndex] != 0.)
             thisIntercalibConstant = theNumerator_EEP[iIndex]/theDenominator_EEP[iIndex];}

        float oldIntercalibConstant = 1.;
        if ( iLoop > 0 ) oldIntercalibConstant = h_scale_hashedIndex_EE -> GetBinContent (iIndex+1);
        h_scale_hashedIndex_EE -> SetBinContent(iIndex+1, thisIntercalibConstant*oldIntercalibConstant);
              
        if ( thisCaliBlock == 0 ) {
          hC_IntercalibValues_EEM -> Fill (iLoop, thisIntercalibConstant);
          hC_PullFromScalib_EEM -> Fill(iLoop,(thisIntercalibConstant*oldIntercalibConstant-1./theScalibration[iIndex]));
          hC_scale_EEM -> Fill(iLoop, GetIxFromHashedIndex(iIndex), GetIyFromHashedIndex(iIndex),thisIntercalibConstant*oldIntercalibConstant);
     
          auxiliary_IC_EEM.Fill(thisIntercalibConstant);
          
        }
        else {
                 if( thisCaliBlock == 1)
                 {
                  hC_IntercalibValues_EEP -> Fill (iLoop, thisIntercalibConstant);
                  hC_PullFromScalib_EEP -> Fill(iLoop,(thisIntercalibConstant*oldIntercalibConstant-1./theScalibration[iIndex]));
                  hC_scale_EEP -> Fill(iLoop, GetIxFromHashedIndex(iIndex), GetIyFromHashedIndex(iIndex),thisIntercalibConstant*oldIntercalibConstant);
  
                  auxiliary_IC_EEP.Fill(thisIntercalibConstant);}
         }


        }
      
      }
   
    g_ICmeanVsLoop_EEM -> SetPoint(iLoop, iLoop, auxiliary_IC_EEM.GetMean());
    g_ICmeanVsLoop_EEM -> SetPointError(iLoop, 0.,auxiliary_IC_EEM.GetMeanError());
    
    g_ICrmsVsLoop_EEM -> SetPoint(iLoop, iLoop, auxiliary_IC_EEM . GetRMS());
    g_ICrmsVsLoop_EEM -> SetPointError(iLoop, 0., auxiliary_IC_EEM . GetRMSError());

    g_ICmeanVsLoop_EEP -> SetPoint(iLoop, iLoop, auxiliary_IC_EEP . GetMean());
    g_ICmeanVsLoop_EEP -> SetPointError(iLoop, 0., auxiliary_IC_EEP . GetMeanError());
    
    g_ICrmsVsLoop_EEP -> SetPoint(iLoop, iLoop, auxiliary_IC_EEP . GetRMS());
    g_ICrmsVsLoop_EEP -> SetPointError(iLoop, 0., auxiliary_IC_EEP . GetRMSError());
    
   }// Calibration Loops
      
   //Fill the histo of IntercalibValues after the loops at last step
   for ( int iIndex = 0; iIndex < kEEhalf*2; iIndex++ ){
           
     if ( h_occupancy_hashedIndex_EE -> GetBinContent(iIndex+1) > 0 ){
        
       int thisCaliBlock = -1;
       if (GetZsideFromHashedIndex(iIndex) < 0) thisCaliBlock = 0;
       else thisCaliBlock = 1;
       
       int thisIx = GetIxFromHashedIndex(iIndex);
       int thisIy = GetIyFromHashedIndex(iIndex);
       int thisIz = GetZsideFromHashedIndex(iIndex);

       float thisIntercalibConstant = h_scale_hashedIndex_EE -> GetBinContent (iIndex+1);
       if ( thisCaliBlock == 0 ) 
         h_scale_EEM -> Fill (thisIx, thisIy, thisIntercalibConstant);
       else
         h_scale_EEP -> Fill (thisIx, thisIy, thisIntercalibConstant);

       if ( thisCaliBlock == 0 )
       {
            //Vectors
         IxValues_EEM.push_back(thisIx);
         IyValues_EEM.push_back(thisIy);
         ICValues_EEM.push_back(thisIntercalibConstant);
       }
       else{
             IxValues_EEP.push_back(thisIx);
             IyValues_EEP.push_back(thisIy);
             ICValues_EEP.push_back(thisIntercalibConstant);
           }

       int thisIr = EERings(thisIx,thisIy,thisIz);
       if(thisIz >0)
       {
        SumIC_Ring_EEP.at(thisIr) = SumIC_Ring_EEP.at(thisIr) + thisIntercalibConstant;
        Sumxtal_Ring_EEP.at(thisIr) = Sumxtal_Ring_EEP.at(thisIr) + 1;
       }
       else{
              SumIC_Ring_EEM.at(thisIr) = SumIC_Ring_EEM.at(thisIr) + thisIntercalibConstant;
              Sumxtal_Ring_EEM.at(thisIr) = Sumxtal_Ring_EEM.at(thisIr) + 1;
           }
               
       }

      
   }
   for ( int iIndex = 0; iIndex < kEEhalf*2; iIndex++ ){
   
    if ( h_occupancy_hashedIndex_EE -> GetBinContent(iIndex+1) > 0 ){
        
       int thisCaliBlock = -1;
       if (GetZsideFromHashedIndex(iIndex) < 0) thisCaliBlock = 0;
       else thisCaliBlock = 1;
       
       int thisIx = GetIxFromHashedIndex(iIndex);
       int thisIy = GetIyFromHashedIndex(iIndex);
       int thisIz = GetZsideFromHashedIndex(iIndex);

       int thisIr = EERings(thisIx,thisIy,thisIz);

       float thisIntercalibConstant = h_scale_hashedIndex_EE -> GetBinContent (iIndex+1);
     
       
       if(thisIz > 0)
       {
          if(Sumxtal_Ring_EEP.at(thisIr) != 0 && SumIC_Ring_EEP.at(thisIr)!= 0)
          h_scale_meanOnring_EEP->Fill(thisIx,thisIy,thisIntercalibConstant/(SumIC_Ring_EEP.at(thisIr)/Sumxtal_Ring_EEP.at(thisIr)));
       }
       else{
            if(Sumxtal_Ring_EEM.at(thisIr) != 0 && SumIC_Ring_EEM.at(thisIr) != 0)
            h_scale_meanOnring_EEM->Fill(thisIx,thisIy,thisIntercalibConstant/(SumIC_Ring_EEM.at(thisIr)/Sumxtal_Ring_EEM.at(thisIr)));
           }
       }
   }
    
}

  
void FastCalibratorEE::saveHistos(TFile * f1)
{

  f1->cd();
  
  // EE+
  hC_IntercalibValues_EEP-> Write(*f1);
  hC_PullFromScalib_EEP->Write(*f1);
  hC_EoP->Write(*f1);
  hC_scale_EEP->Write("",*f1);
  
  h_occupancy_EEP->Write();
  h_scale_EEP->Write();
  
  h_scale_hashedIndex_EE->Write();
  h_occupancy_hashedIndex_EE->Write();
  h_map_Dead_Channels_EEP->Write();
    
  g_ICmeanVsLoop_EEP->Write();
  g_ICrmsVsLoop_EEP->Write();
  h_scale_meanOnring_EEP->Write("h_scale_map_EEP");
  
  // EE-
  hC_IntercalibValues_EEM-> Write(*f1);
  hC_scale_EEM->Write("",*f1);
  hC_PullFromScalib_EEM->Write(*f1);
  h_occupancy_EEM->Write();
  h_scale_EEM->Write();
  h_map_Dead_Channels_EEP->Write();
    
  g_ICmeanVsLoop_EEM->Write();
  g_ICrmsVsLoop_EEM->Write();
  h_scale_meanOnring_EEM->Write("h_scale_map_EEM");
  

  f1->Close();

  return;
}


void FastCalibratorEE::saveEoPeta(TFile * f2)
{
 f2->cd();
 hC_EoP_ir_ele ->Write(*f2);
 f2->Close(); 
}

void FastCalibratorEE::printOnTxt(TString outputTxtFile)
{
  std::ofstream outTxt (outputTxtFile.Data(),std::ios::out);

  outTxt << "---------------------------------------------------------------" << std::endl;
  outTxt << "--- ix ---- iy ------ IC value  --- EE+ ------" << std::endl;
  outTxt << "---------------------------------------------------------------" << std::endl;

  for (unsigned int iIndex = 0; iIndex < ICValues_EEP.size(); iIndex++)
  {
    int thisIr = EERings(IxValues_EEP.at(iIndex),IyValues_EEP.at(iIndex),1);
    if(Sumxtal_Ring_EEP.at(thisIr) == 0 || SumIC_Ring_EEP.at(thisIr) == 0) continue;
   
    outTxt << "  " << std::fixed << std::setw(4)  << std::right << IxValues_EEP.at(iIndex) 
           << std::fixed << std::setw(10)  << std::right << IyValues_EEP.at(iIndex) 
           << "          " << (float) ICValues_EEP.at(iIndex)/(SumIC_Ring_EEP.at(thisIr)/Sumxtal_Ring_EEP.at(thisIr))<< std::endl;
  }
  
  outTxt << "---------------------------------------------------------------" << std::endl;
  outTxt << "--- ix ---- iy ------ IC value  --- EE- -----" << std::endl;
  outTxt << "---------------------------------------------------------------" << std::endl;

  for (unsigned int iIndex = 0; iIndex < ICValues_EEM.size(); iIndex++)
  {
   int thisIr = EERings(IxValues_EEM.at(iIndex),IyValues_EEM.at(iIndex),-1);
   if(Sumxtal_Ring_EEM.at(thisIr) == 0 || SumIC_Ring_EEM.at(thisIr) == 0) continue;
       
    outTxt << "  " << std::fixed << std::setw(4)  << std::right << IxValues_EEM.at(iIndex) 
           << std::fixed << std::setw(10)  << std::right << IyValues_EEM.at(iIndex) 
           << "          " << (float) ICValues_EEM.at(iIndex)/(SumIC_Ring_EEM.at(thisIr)/Sumxtal_Ring_EEM.at(thisIr)) << std::endl;
  }
  
}

void FastCalibratorEE::AcquireDeadXtal(TString inputDeadXtal)
{
  if(inputDeadXtal!="NULL")
  {
   std::ifstream DeadXtal (inputDeadXtal.Data(),std::ios::binary);
   
   std::string buffer;
   int iX, iY ,iZ;
  

   while(!DeadXtal.eof())
   {
    getline(DeadXtal,buffer);
    std::stringstream line( buffer );
    line >> iX >> iY >>iZ ;
    DeadXtal_HashedIndex.push_back(GetHashedIndexEE(iX,iY,iZ)) ;
   
   }

  sort(DeadXtal_HashedIndex.begin(), DeadXtal_HashedIndex.end());
  }
  else{
       DeadXtal_HashedIndex.push_back(-9999);
      }

}

bool FastCalibratorEE::CheckDeadXtal(const int & iX, const int & iY, const int & iZ)
{
  int hashed_Index;
  hashed_Index = GetHashedIndexEE(iX,iY,iZ);
  
  std::vector<int>::iterator iter = find(DeadXtal_HashedIndex.begin(),DeadXtal_HashedIndex.end(),hashed_Index);

  if(iter !=DeadXtal_HashedIndex.end())
     return true;
  else return false;
}
   

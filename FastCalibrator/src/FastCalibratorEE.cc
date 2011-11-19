#include "FastCalibratorEE.h"
#include "VEcalCalibBlock.h"
#include "L3CalibBlock.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <TRandom.h>

const int FastCalibratorEE::kxf[] =  {
    41,  51,  41,  51,  41,  51,  36,  51,  36,  51,
    26,  51,  26,  51,  26,  51,  21,  51,  21,  51,
    21,  51,  21,  51,  21,  51,  16,  51,  16,  51,
    14,  51,  14,  51,  14,  51,  14,  51,  14,  51,
    9,  51,   9,  51,   9,  51,   9,  51,   9,  51,
    6,  51,   6,  51,   6,  51,   6,  51,   6,  51,
    6,  51,   6,  51,   6,  51,   6,  51,   6,  51,
    4,  51,   4,  51,   4,  51,   4,  51,   4,  56,
    1,  58,   1,  59,   1,  60,   1,  61,   1,  61,
    1,  62,   1,  62,   1,  62,   1,  62,   1,  62,
    1,  62,   1,  62,   1,  62,   1,  62,   1,  62,
    1,  61,   1,  61,   1,  60,   1,  59,   1,  58,
    4,  56,   4,  51,   4,  51,   4,  51,   4,  51,
    6,  51,   6,  51,   6,  51,   6,  51,   6,  51,
    6,  51,   6,  51,   6,  51,   6,  51,   6,  51,
    9,  51,   9,  51,   9,  51,   9,  51,   9,  51,
    14,  51,  14,  51,  14,  51,  14,  51,  14,  51,
    16,  51,  16,  51,  21,  51,  21,  51,  21,  51,
    21,  51,  21,  51,  26,  51,  26,  51,  26,  51,
    36,  51,  36,  51,  41,  51,  41,  51,  41,  51
} ;

const int FastCalibratorEE::kdi[] =  {
  0,   10,   20,   30,   40,   50,   60,   75,   90,  105,
  120,  145,  170,  195,  220,  245,  270,  300,  330,  360,
  390,  420,  450,  480,  510,  540,  570,  605,  640,  675,
  710,  747,  784,  821,  858,  895,  932,  969, 1006, 1043,
  1080, 1122, 1164, 1206, 1248, 1290, 1332, 1374, 1416, 1458,
  1500, 1545, 1590, 1635, 1680, 1725, 1770, 1815, 1860, 1905,
  1950, 1995, 2040, 2085, 2130, 2175, 2220, 2265, 2310, 2355,
  2400, 2447, 2494, 2541, 2588, 2635, 2682, 2729, 2776, 2818,
  2860, 2903, 2946, 2988, 3030, 3071, 3112, 3152, 3192, 3232,
  3272, 3311, 3350, 3389, 3428, 3467, 3506, 3545, 3584, 3623,
  3662, 3701, 3740, 3779, 3818, 3857, 3896, 3935, 3974, 4013,
  4052, 4092, 4132, 4172, 4212, 4253, 4294, 4336, 4378, 4421,
  4464, 4506, 4548, 4595, 4642, 4689, 4736, 4783, 4830, 4877,
  4924, 4969, 5014, 5059, 5104, 5149, 5194, 5239, 5284, 5329,
  5374, 5419, 5464, 5509, 5554, 5599, 5644, 5689, 5734, 5779,
  5824, 5866, 5908, 5950, 5992, 6034, 6076, 6118, 6160, 6202,
  6244, 6281, 6318, 6355, 6392, 6429, 6466, 6503, 6540, 6577,
  6614, 6649, 6684, 6719, 6754, 6784, 6814, 6844, 6874, 6904,
  6934, 6964, 6994, 7024, 7054, 7079, 7104, 7129, 7154, 7179,
  7204, 7219, 7234, 7249, 7264, 7274, 7284, 7294, 7304, 7314
} ;


FastCalibratorEE::FastCalibratorEE(TTree *tree)
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
  ele2_recHit_E = 0;
  ele2_recHit_hashedIndex = 0;
  ele2_recHit_iphiORiy = 0;
  ele2_recHit_ietaORix = 0;
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

  //EE+

  hC_IntercalibValues_EEP = new hChain ("IntercalibValues_EEP", "IntercalibValues_EEP", 400,0.2,1.9, nLoops);
  hC_EoP_EEP = new hChain ("EoP_EEP", "EoP_EEP", 400,0.2,1.9, nLoops);
  hC_PullFromScalib_EEP = new hChain ("hC_PullFromScalib_EEP", "hC_PullFromScalib_EEP", 2000,-0.5,0.5, nLoops);
  hC_scale_EEP = new h2Chain("hC_scale_EEP", "hC_scale_EEP", 100,1, 101, 100, 1, 101, nLoops );
  
  h_scale_EEP = new TH2F("h_scale_EEP", "h_scale_EEP", 100,1, 101, 100, 1, 101 );
  h_occupancy_EEP = new TH2F("h_occupancy_EEP", "h_occupancy_EEP", 100,1, 101, 100, 1, 101 );
  
  g_ICmeanVsLoop_EEP = new TGraphErrors();
  g_ICrmsVsLoop_EEP = new TGraphErrors();
  
  //EE-
  hC_IntercalibValues_EEM = new hChain ("IntercalibValues_EEM", "IntercalibValues_EEM", 400,0.2,1.9, nLoops);
  hC_EoP_EEM = new hChain ("EoP_EEM", "EoP_EEM", 100,0.2,1.9, nLoops);
  hC_PullFromScalib_EEM = new hChain ("hC_PullFromScalib_EEM", "hC_PullFromScalib_EEM", 2000,-0.5,0.5, nLoops);
  hC_scale_EEM = new h2Chain("hC_scale_EEM", "hC_scale_EEM", 100,1, 101, 100, 1, 101, nLoops );
  
  h_scale_EEM = new TH2F("h_scale_EEM", "h_scale_EEM", 100,1, 101, 100, 1, 101 );
  h_occupancy_EEM = new TH2F("h_occupancy_EEM", "h_occupancy_EEM", 100,1, 101, 100, 1, 101 );
  
  g_ICmeanVsLoop_EEM = new TGraphErrors();
  g_ICrmsVsLoop_EEM = new TGraphErrors();

  return;
}

int FastCalibratorEE::GetHashedIndexEE(int ix, int iy, int zside)
{
  int jx = ix ;
  int jd =  2*( iy - 1 ) + ( jx - 1 )/50  ;
  return (  ( zside < 0 ? 0 : kEEhalf ) + kdi[jd] + jx - kxf[jd] ) ;
}

int FastCalibratorEE::GetIxFromHashedIndex(int Index)
{
  int di = Index%kEEhalf ;
  int ii = ( std::upper_bound( kdi, kdi+(2*IY_MAX), di ) - kdi ) - 1 ;
  return ( kxf[ii] + di - kdi[ii] ) ;
}

int FastCalibratorEE::GetIyFromHashedIndex(int Index)
{
  int di = Index%kEEhalf ;
  int ii = ( std::upper_bound( kdi, kdi+(2*IY_MAX), di ) - kdi ) - 1 ;
  return ( 1 + ii/2 ) ;
}

int FastCalibratorEE::GetZsideFromHashedIndex(int Index)
{
  return ( Index < kEEhalf ? -1 : 1) ;
}

void FastCalibratorEE::Loop(int nentries, int useZ, int useW, int splitStat, int nLoops)
{
   if (fChain == 0) return;
   
   // Define the number of crystal you want to calibrate
   int m_regions = kEEhalf;
   
   std::cout << "m_regions " << m_regions << std::endl;
  
     // build up scalibration map
   std::vector<float> theScalibration(m_regions, 0.);
   TRandom genRand;
   for ( int iIndex = 0; iIndex < m_regions; iIndex++ )  
//      theScalibration[iIndex] = genRand.Gaus(1.,0.05);
     theScalibration[iIndex] = 1.;

  
  
   /// ----------------- Calibration Loops -----------------------------//
   for ( int iLoop = 0; iLoop < nLoops; iLoop++ ) {
    // loop over events
    std::cout << "Starting iteration " << iLoop + 1 << std::endl;
  
    std::vector<float> theNumerator_EEP(m_regions*2, 0.);
    std::vector<float> theDenominator_EEP(m_regions*2, 0.);
    std::vector<float> theNumerator_EEM(m_regions*2, 0.);
    std::vector<float> theDenominator_EEM(m_regions*2, 0.);


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
        
        // Prepare the map : 
        // index = hashIndex = (MAX_IETA + (positiveZ() ? ietaAbs()-1 : -ietaAbs()) )*MAX_IPHI+ iphi()-1 
        // index = hashIndex = (MAX_IETA + ieta() [ieta-1 if ieta >= 0] )*MAX_IPHI+ iphi()-1 
        // MIN_HASH =  0; // always 0 ...
        // MAX_HASH =  2*MAX_IPHI*MAX_IETA-1; = 2*360*85-1 = 61199
        std::map<int,double> map;
        bool skipElectron=false;
        // Only W only Endcap
        if ( ele1_isEB == 0 && (( useW == 1 && isW == 1 ) || ( useZ == 1 && isZ == 1 )) ) {
                  
          // SCL energy containment correction
          FdiEta = ele1_scE/ele1_scERaw;
          // Electron energy
          float thisE = 0;
          float thisE3x3 =0 ;
          int iseed = 0 ;
          float E_seed = 0;


          for (unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
            
            if(ele1_recHit_E -> at(iRecHit) > E_seed)
             {
              E_seed=ele1_recHit_E -> at(iRecHit);
              iseed=iRecHit;
             }
          }
        
          
          // Cycle on the all the recHits of the Event
          for (unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
            
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);
            thisE += ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
            
            if(fabs(ele1_recHit_ietaORix->at(iRecHit)-ele1_recHit_ietaORix->at(iseed))<=1 && 
               fabs(ele1_recHit_iphiORiy->at(iRecHit)-ele1_recHit_iphiORiy->at(iseed))<=1)
              thisE3x3+=theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
                 
              }
            
          pSub = 0.; //NOTALEO : test dummy
          pIn = ele1_tkP;
                
          // find the zside
          int thisCaliBlock = -1;
          if (GetZsideFromHashedIndex(ele1_recHit_hashedIndex -> at(iseed)) < 0) thisCaliBlock = 0;
          else thisCaliBlock = 1;
 
          // discard electrons with bad E/P
            if ( fabs(thisE/ele1_tkP - 1) > 0.3 ) skipElectron = true;
//            if ( fabs(thisE3x3/thisE) < 0.9 ) skipElectron = true;
          
          if ( !skipElectron ) {
                  
              for ( unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
  
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
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E->at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
                }
                
                if(thisCaliBlock == 1) {
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E->at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }
             
              }
             
             
          // use evens    
          if ( splitStat == 1 && jentry%2 == 0 ) {
                  
                if(thisCaliBlock == 0) {
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
                }
                
                if(thisCaliBlock == 1) {
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
         
              }
            }
             
           // use odd    
           if ( splitStat == 1 && jentry%2 != 0 ) {
                  
                if(thisCaliBlock == 0) {
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
                }
                
                if(thisCaliBlock == 1) {
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
         
              }
            }
            }
          }
          //Fill EoP
           if (thisCaliBlock == 0) hC_EoP_EEM -> Fill(iLoop, thisE/pIn);
           else hC_EoP_EEP -> Fill(iLoop, thisE/pIn);
        
        }  
         skipElectron = false;     
        /// Fill the map with the ele (if any) from the second Z leg
        // Only Z only Barrel
        if ( ele2_isEB == 1 && ( useZ == 1 && isZ == 1 ) ) {
          
           // SCL energy containment correction
          FdiEta = ele2_scE/ele2_scERaw;
          // Electron energy
          float thisE = 0;
          float thisE3x3 =0 ;
          float E_seed = 0;
          float iseed = 0 ;
          float thisIC = 0 ;

          for (unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
             
            if(ele2_recHit_E -> at(iRecHit) > E_seed)
             {
              E_seed = ele2_recHit_E -> at(iRecHit) ;
              iseed=iRecHit;
             }
          }
        

          // Cycle on the all the recHits of the Event: to get the old IC and the corrected SC energy
          for (unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);

            if ( iLoop == 0 ) {
              h_occupancy_hashedIndex_EE -> Fill(thisIndex);
              if ( GetZsideFromHashedIndex(thisIndex) < 0 ) 
                h_occupancy_EEM -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
              else h_occupancy_EEP -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
            
            }

            
            thisE += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;

           if(fabs(ele2_recHit_ietaORix->at(iRecHit)-ele2_recHit_ietaORix->at(iseed))<=1 && 
              fabs(ele2_recHit_iphiORiy->at(iRecHit)-ele2_recHit_iphiORiy->at(iseed))<=1)
             thisE3x3+=theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;
             
 
          }

            
          pSub = 0.; //NOTALEO : test dummy
          pIn = ele2_tkP;
          
          // find the zside
          int thisCaliBlock = -1;
          if (GetZsideFromHashedIndex(iseed) < 0) thisCaliBlock = 0;
          else thisCaliBlock = 1;
          
          // discard electrons with bad E/P
           if ( fabs(thisE/ele1_tkP - 1) > 0.3 ) skipElectron = true;
//            if ( fabs(thisE3x3/thisE) < 0.9 ) skipElectron = true;
        
          if ( !skipElectron ) {
                  
              for ( unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
  
              int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);
              float thisIC = 1.;

              if (iLoop > 0) thisIC = h_scale_hashedIndex_EE -> GetBinContent(thisIndex+1);

             if ( iLoop == 0 ) {
              h_occupancy_hashedIndex_EE -> Fill(thisIndex);
              if ( GetZsideFromHashedIndex(thisIndex) < 0 ) 
                h_occupancy_EEM -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
              else h_occupancy_EEP -> Fill(GetIxFromHashedIndex(thisIndex), GetIyFromHashedIndex(thisIndex) );
            
            }
   
  
              if ( splitStat == 0) {
             
                if(thisCaliBlock == 0) {
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E->at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
                }
                
                if(thisCaliBlock == 1) {
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E->at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }
             
              }
             
             
          // use evens    
           if ( splitStat == 1 && jentry%2 == 0 ) {
                  
                if(thisCaliBlock == 0) {
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
                }
                
                if(thisCaliBlock == 1) {
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
         
              }
            }
             
           // use odd    
           if ( splitStat == 1 && jentry%2 != 0 ) {
                  
                if(thisCaliBlock == 0) {
                theNumerator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEM[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
                }
                
                if(thisCaliBlock == 1) {
                theNumerator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator_EEP[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
         
              }
            }
           }
          
          
          //Fill EoP
           if (thisCaliBlock == 0) hC_EoP_EEM -> Fill(iLoop, thisE/pIn);
           else hC_EoP_EEP -> Fill(iLoop, thisE/pIn);
         
        }
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
        
        float thisIntercalibConstant = 1.;
        
        if(thisCaliBlock == 0 && theDenominator_EEM[iIndex] != 0.) 
        thisIntercalibConstant = theNumerator_EEM[iIndex]/theDenominator_EEM[iIndex];

        if(thisCaliBlock == 1 && theDenominator_EEP[iIndex] != 0.)
        thisIntercalibConstant = theNumerator_EEP[iIndex]/theDenominator_EEP[iIndex];

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
          hC_IntercalibValues_EEP -> Fill (iLoop, thisIntercalibConstant);
          hC_PullFromScalib_EEP -> Fill(iLoop,(thisIntercalibConstant*oldIntercalibConstant-1./theScalibration[iIndex]));
          hC_scale_EEP -> Fill(iLoop, GetIxFromHashedIndex(iIndex), GetIyFromHashedIndex(iIndex),thisIntercalibConstant*oldIntercalibConstant);
  
          auxiliary_IC_EEP.Fill(thisIntercalibConstant);
         }


        }
      
      }
     
    g_ICmeanVsLoop_EEM -> SetPoint(iLoop, iLoop, auxiliary_IC_EEM . GetMean());
    g_ICmeanVsLoop_EEM -> SetPointError(iLoop, 0., auxiliary_IC_EEM . GetMeanError());
    
    g_ICrmsVsLoop_EEM -> SetPoint(iLoop, iLoop, auxiliary_IC_EEM . GetRMS());
    g_ICrmsVsLoop_EEM -> SetPointError(iLoop, 0., auxiliary_IC_EEM . GetRMSError());

    g_ICmeanVsLoop_EEP -> SetPoint(iLoop, iLoop, auxiliary_IC_EEP . GetMean());
    g_ICmeanVsLoop_EEP -> SetPointError(iLoop, 0., auxiliary_IC_EEP . GetMeanError());
    
    g_ICrmsVsLoop_EEP -> SetPoint(iLoop, iLoop, auxiliary_IC_EEP . GetRMS());
    g_ICrmsVsLoop_EEP -> SetPointError(iLoop, 0., auxiliary_IC_EEP . GetRMSError());
    
   }// Calibration Loops
      
   //Fill the histo of IntercalibValues after the loops
   for ( int iIndex = 0; iIndex < kEEhalf*2; iIndex++ ){
           
     if ( h_occupancy_hashedIndex_EE -> GetBinContent(iIndex+1) > 0 ){
        
       int thisCaliBlock = -1;
       if (GetZsideFromHashedIndex(iIndex) < 0) thisCaliBlock = 0;
       else thisCaliBlock = 1;
       
       int thisIx = GetIxFromHashedIndex(iIndex);
       int thisIy = GetIyFromHashedIndex(iIndex);

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
         
       }

      
     }

}

  
void FastCalibratorEE::saveHistos(TFile * f1)
{

  f1->cd();
  
  h_scale_hashedIndex_EE->Write();
  h_occupancy_hashedIndex_EE->Write();

  // EE+
  hC_IntercalibValues_EEP-> Write(*f1);
  hC_EoP_EEP->Write(*f1);
  hC_scale_EEP->Write("",*f1);
  hC_PullFromScalib_EEP->Write(*f1);
  h_occupancy_EEP->Write();
  h_scale_EEP->Write();
    
  g_ICmeanVsLoop_EEP->Write();
  g_ICrmsVsLoop_EEP->Write();

  // EE-
  hC_IntercalibValues_EEM-> Write(*f1);
  hC_EoP_EEM->Write(*f1);
  hC_scale_EEP->Write("",*f1);
  hC_PullFromScalib_EEM->Write(*f1);

  h_occupancy_EEM->Write();
  h_scale_EEM->Write();
    
  g_ICmeanVsLoop_EEM->Write();
  g_ICrmsVsLoop_EEM->Write();

  f1->Close();

  return;
}

void FastCalibratorEE::printOnTxt(std::string outputTxtFile)
{
  std::ofstream outTxt (outputTxtFile.c_str(),std::ios::out);

  outTxt << "---------------------------------------------------------------" << std::endl;
  outTxt << "--- ix ---- iy ------ IC value  --- EE+ ------" << std::endl;
  outTxt << "---------------------------------------------------------------" << std::endl;

  for (unsigned int iIndex = 0; iIndex < ICValues_EEP.size(); iIndex++)
    outTxt << "  " << std::fixed << std::setw(4)  << std::right << IxValues_EEP.at(iIndex) 
           << std::fixed << std::setw(10)  << std::right << IyValues_EEP.at(iIndex) 
           << "          " << (float) ICValues_EEP.at(iIndex) << std::endl;
  
  outTxt << "---------------------------------------------------------------" << std::endl;
  outTxt << "--- ix ---- iy ------ IC value  --- EE- -----" << std::endl;
  outTxt << "---------------------------------------------------------------" << std::endl;

  for (unsigned int iIndex = 0; iIndex < ICValues_EEM.size(); iIndex++)
    outTxt << "  " << std::fixed << std::setw(4)  << std::right << IxValues_EEM.at(iIndex) 
           << std::fixed << std::setw(10)  << std::right << IyValues_EEM.at(iIndex) 
           << "          " << (float) ICValues_EEM.at(iIndex) << std::endl;
  
}

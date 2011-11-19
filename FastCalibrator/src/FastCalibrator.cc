#include "FastCalibrator.h"
#include "GetHashedIndexEB.h"
#include "VEcalCalibBlock.h"
#include "L3CalibBlock.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <iostream>
#include <fstream>
#include <vector>


FastCalibrator::FastCalibrator(TTree *tree)
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

FastCalibrator::~FastCalibrator()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t FastCalibrator::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t FastCalibrator::LoadTree(Long64_t entry)
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

void FastCalibrator::Init(TTree *tree)
{
  // Set object pointer
  ele1_recHit_E = 0;
  ele1_recHit_hashedIndex = 0;
  ele1_recHit_iphiORiy = 0;
  ele1_recHit_ietaORix = 0;
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
  fChain->SetBranchAddress("ele1_recHit_ietaORix", &ele1_recHit_ietaORix, &b_ele1_recHit_ietaORix);
  fChain->SetBranchAddress("ele1_recHit_iphiORiy", &ele1_recHit_iphiORiy, &b_ele1_recHit_iphiORiy);

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

Bool_t FastCalibrator::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  return kTRUE;
}

void FastCalibrator::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t FastCalibrator::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  return 1;
}

void FastCalibrator::bookHistos(int nLoops)
{
  hC_IntercalibValues = new hChain ("IntercalibValues", "IntercalibValues", 2000,0.5,1.5, nLoops);
  hC_PullFromScalib = new hChain ("hC_PullFromScalib", "hC_PullFromScalib", 2000,-0.5,0.5, nLoops);
  hC_EoP = new hChain ("EoP", "EoP", 100,0.2,1.9, nLoops);
  hC_scale_EB = new h2Chain("hC_scale_EB", "hC_scale_EB", 360,1, 361, 171, -85, 86, nLoops );
  h_scalib_EB = new TH2F("h_scalib_EB", "h_scalib_EB", 360,1, 361, 171, -85, 86 );
  h_Occupancy_hashedIndex = new TH1F ("h_Occupancy_hashedIndex", "h_Occupancy_hashedIndex", 61201,-0.5,61199.5);
  p_IntercalibValues_iEta = new TProfile ("p_IntercalibValues_iEta","p_IntercalibValues_iEta", 171, -85, 86, -0.1, 2.1);
  h_IntercalibSpread_iEta = new TH1F ("h_IntercalibSpread_iEta", "h_IntercalibSpread_iEta", 171, -85, 86);
  h_IntercalibValues_test = new TH1F ("h_IntercalibValues_test", "h_IntercalibValues_test", 400, -1, 1);
  h_scale_EB_hashedIndex = new TH1F("h_scale_EB_hashedIndex", "h_scale_EB_hashedIndex", 61201,-0.5,61999.5 );
  h_Init_IntercalibValues = new TH1F("h_Init_IntercalibValues","h_Init_IntercalibValues",2000,0.5,1.5);

  g_ICmeanVsLoop = new TGraphErrors();
  g_ICmeanVsLoop -> SetName("g_ICmeanVsLoop");
  g_ICmeanVsLoop -> SetTitle("g_ICmeanVsLoop");

  g_ICrmsVsLoop  = new TGraphErrors();
  g_ICrmsVsLoop -> SetName("g_ICrmsVsLoop");
  g_ICrmsVsLoop ->SetTitle("g_ICrmsVsLoop");
  // essential plot for validation
  h_scale_EB = new TH2F("h_scale_EB", "h_scale_EB", 360,1, 361, 171, -85, 86 );
  h_scale_EB_meanOnPhi = new TH2F("h_scale_EB_meanOnPhi", "h_scale_EB_meanOnPhi", 360,1, 361, 171, -85, 86 );
  h_occupancy = new TH2F("h_occupancy", "h_occupancy", 360,1, 361, 171, -85, 86 );
  
  return;
}


void FastCalibrator::Loop(int nentries, int useZ, int useW, int splitStat, int nLoops)
{
   if (fChain == 0) return;
   
   // Define the number of crystal you want to calibrate
   int m_regions = 0;
 

   // Define useful numbers
   static const int MIN_IETA = 1;
   static const int MIN_IPHI = 1;
   static const int MAX_IETA = 85;
   static const int MAX_IPHI = 360;
   //In this example set everything at 1 for EB
   for ( int iabseta = MIN_IETA; iabseta <= MAX_IETA; iabseta++ ){
     for ( int iphi = MIN_IPHI; iphi <= MAX_IPHI; iphi++ ){
       for ( int theZside = -1; theZside < 2; theZside = theZside+2 ){
         
         m_regions++;
         
       }
     }
   }
   
   std::cout << "m_regions " << m_regions << std::endl;
   
   //Prepare the calibration blocks
   int eventWeight = 2; //Pres8 prenscription

   // build up scalibration map
   std::vector<float> theScalibration(m_regions, 0.);
   TRandom genRand;
   for ( int iIndex = 0; iIndex < m_regions; iIndex++ )  
    { theScalibration[iIndex] = genRand.Gaus(1.,0.05);
//      theScalibration[iIndex] = 1.;
       h_Init_IntercalibValues->Fill(theScalibration[iIndex]);
       h_scalib_EB -> Fill ( GetIphiFromHashedIndex(iIndex), GetIetaFromHashedIndex(iIndex), theScalibration[iIndex] );
   }

   /// ----------------- Calibration Loops -----------------------------//
   for ( int iLoop = 0; iLoop < nLoops; iLoop++ ) {

    std::cout << "Starting iteration " << iLoop + 1 << std::endl;
    
    // prepare the numerator and denominator for each Xtal
    std::vector<float> theNumerator(m_regions, 0.);
    std::vector<float> theDenominator(m_regions, 0.);
    
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
       
       
        // Only W only Barrel
        if ( ele1_isEB == 1 && (( useW == 1 && isW == 1 ) || ( useZ == 1 && isZ == 1 )) ) {
                  
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
          /// Fill the map with the ele from W or the first Z leg
            
          // Cycle on the all the recHits of the Event: to get the old IC and the corrected SC energy
          for (unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
            
            float thisIC = 1.;
            int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
            // IC obtained from previous Loops
            if (iLoop > 0) thisIC = h_scale_EB_hashedIndex -> GetBinContent(thisIndex+1);
            
            thisE += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
            if(fabs(ele1_recHit_ietaORix->at(iRecHit)-ele1_recHit_ietaORix->at(iseed))<=1 && 
               fabs(ele1_recHit_iphiORiy->at(iRecHit)-ele1_recHit_iphiORiy->at(iseed))<=1)
              thisE3x3+=theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC;
              
          }
          
         
 
          pSub = 0.; //NOTALEO : test dummy
          pIn = ele1_tkP;
//           pIn = ele1_scE;
          
          // discard electrons with bad E/P
            if ( fabs(thisE/ele1_tkP - 1) > 0.1 ) skipElectron = true;
//             if ( fabs(thisE3x3/thisE) < 0.8 ) skipElectron = true;
        
          if ( !skipElectron ) {
          
            // Now cycle on the all the recHits and update the numerator and denominator
            for ( unsigned int iRecHit = 0; iRecHit < ele1_recHit_E->size(); iRecHit++ ) {
  
              int thisIndex = ele1_recHit_hashedIndex -> at(iRecHit);
              float thisIC = 1.;
              if (iLoop > 0) thisIC = h_scale_EB_hashedIndex -> GetBinContent(thisIndex+1);
  
              // Fill the occupancy map JUST for the first Loop
              if ( iLoop == 0 ) {
                h_Occupancy_hashedIndex -> Fill(thisIndex);  
                h_occupancy -> Fill(GetIphiFromHashedIndex(thisIndex), GetIetaFromHashedIndex(thisIndex));
              }
  
              // use full statistics
              if ( splitStat == 0 ) {
                theNumerator[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }
              // use evens    
              else if ( splitStat == 1 && jentry%2 == 0 ) {
                theNumerator[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }  
              // use odds
              else if ( splitStat == -1 && jentry%2 != 0 ) {
                theNumerator[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator[thisIndex] += theScalibration[thisIndex]*ele1_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }
            
            }
          
          }

          
          //Fill EoP
          hC_EoP -> Fill(iLoop, thisE/pIn);
        
        }  
              
        /// Fill the map with the ele (if any) from the second Z leg
        skipElectron = false;
      
        // Only Z only Barrel
        if ( ele2_isEB == 1 && (( useW == 1 && isW == 1 ) || ( useZ == 1 && isZ == 1 )) ) {
  
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
            if (iLoop > 0) thisIC = h_scale_EB_hashedIndex -> GetBinContent(thisIndex+1);

            if ( iLoop == 0 ) {
              h_Occupancy_hashedIndex -> Fill(thisIndex);
              h_occupancy -> Fill(GetIphiFromHashedIndex(thisIndex), GetIetaFromHashedIndex(thisIndex));
            }

            
            thisE += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;

           if(fabs(ele2_recHit_ietaORix->at(iRecHit)-ele2_recHit_ietaORix->at(iseed))<=1 && 
              fabs(ele2_recHit_iphiORiy->at(iRecHit)-ele2_recHit_iphiORiy->at(iseed))<=1)
             thisE3x3+=theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC;
             
 
          }

          pSub = 0.; //NOTALEO : test dummy
          pIn = ele2_tkP;
//           pIn = ele2_scE;
          
          // discard electrons with bad E/P
            if ( fabs(thisE/ele2_tkP - 1) > 0.1 ) skipElectron = true;
//            if ( fabs(thisE3x3/thisE) < 0.8 ) skipElectron = true;
    
          if ( !skipElectron ) {
          
            // Now cycle on the all the recHits and update the numerator and denominator
            for ( unsigned int iRecHit = 0; iRecHit < ele2_recHit_E->size(); iRecHit++ ) {
  
              int thisIndex = ele2_recHit_hashedIndex -> at(iRecHit);
              float thisIC = 1.;
              if (iLoop > 0) thisIC = h_scale_EB_hashedIndex -> GetBinContent(thisIndex+1);
  
              // Fill the occupancy map JUST for the first Loop
              if ( iLoop == 0 ) {
                h_Occupancy_hashedIndex -> Fill(thisIndex);  
                h_occupancy -> Fill(GetIphiFromHashedIndex(thisIndex), GetIetaFromHashedIndex(thisIndex));
              }
  
              // use full statistics
              if ( splitStat == 0 ) {
                theNumerator[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }
              // use evens    
              else if ( splitStat == 1 && jentry%2 == 0 ) {
                theNumerator[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }  
              // use odds
              else if ( splitStat == -1 && jentry%2 != 0 ) {
                theNumerator[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE*(pIn-pSub)/thisE;
                theDenominator[thisIndex] += theScalibration[thisIndex]*ele2_recHit_E -> at(iRecHit)*FdiEta*thisIC/thisE;
              }
            
            }
          
          }

          
          //Fill EoP
          hC_EoP -> Fill(iLoop, thisE/pIn);

        }
          
    }// Cycle on the events
       
    //Now solve everything
    std::cout << ">>>>> [L3][endOfLoop] entering..." << std::endl;

    TH1F auxiliary_IC("auxiliary_IC","auxiliary_IC",50,0.2,1.9);

    //Fill the histo of IntercalibValues before the solve
    for ( int iIndex = 0; iIndex < 61200; iIndex++ ){
      
     
      if ( h_Occupancy_hashedIndex -> GetBinContent(iIndex+1) > 0 ){
	
        float thisIntercalibConstant = 1.;
        // Solve the cases where the recHit energy is always 0 (dead Xtal?)
        if ( theDenominator[iIndex] != 0. ) thisIntercalibConstant = theNumerator[iIndex]/theDenominator[iIndex];
        float oldIntercalibConstant = 1.;
        if ( iLoop > 0 ) oldIntercalibConstant = h_scale_EB_hashedIndex -> GetBinContent (iIndex+1);
        
	h_scale_EB_hashedIndex -> SetBinContent(iIndex+1, thisIntercalibConstant*oldIntercalibConstant);
	hC_IntercalibValues -> Fill(iLoop, thisIntercalibConstant);
        hC_PullFromScalib -> Fill(iLoop,(thisIntercalibConstant*oldIntercalibConstant-1./theScalibration[iIndex]));
        hC_scale_EB -> Fill(iLoop, GetIphiFromHashedIndex(iIndex), GetIetaFromHashedIndex(iIndex),thisIntercalibConstant*oldIntercalibConstant);
        
        //Save the new IC coefficient
        auxiliary_IC.Fill(thisIntercalibConstant);

      }

    }

    g_ICmeanVsLoop -> SetPoint(iLoop, iLoop, auxiliary_IC . GetMean());
    g_ICmeanVsLoop -> SetPointError(iLoop, 0., auxiliary_IC . GetMeanError());
    
    g_ICrmsVsLoop -> SetPoint(iLoop, iLoop, auxiliary_IC . GetRMS());
    g_ICrmsVsLoop -> SetPointError(iLoop, 0., auxiliary_IC . GetRMSError());
   }// end calibration loop
      
   
   int myPhiIndex = 0;

   //Fill the histo of IntercalibValues after the solve: Cycle on iEtaiPhi
   for ( int iabseta = MIN_IETA; iabseta < MAX_IETA + 1; iabseta++ ){
     for ( int theZside = -1; theZside < 2; theZside = theZside+2 ){

     //Setup the histo for fit
     TH1F histoAuxiliary("histoAuxiliary","histoAuxiliary",400, 0.2, 1.9);
     TF1 f1("f1","gaus",0.2,1.9);
     
     int totIphi = 0;
     float meanICforPhiRing = 0.;

       for ( int iphi = MIN_IPHI; iphi <= MAX_IPHI; iphi++ ){
         

         int thisHashedIndex = GetHashedIndexEB(iabseta*theZside, iphi, theZside);
	 if ( h_Occupancy_hashedIndex -> GetBinContent(thisHashedIndex+1) == 0 ) continue;
	 float thisIntercalibConstant = h_scale_EB_hashedIndex -> GetBinContent (thisHashedIndex+1);
         h_scale_EB -> Fill(iphi, iabseta*theZside, thisIntercalibConstant);
         
         if (GetIetaFromHashedIndex(thisHashedIndex) == 85) h_IntercalibValues_test -> Fill (thisIntercalibConstant);
         p_IntercalibValues_iEta -> Fill(GetIetaFromHashedIndex(thisHashedIndex), thisIntercalibConstant);
         
         histoAuxiliary . Fill (thisIntercalibConstant);
         
         //Vectors
         IetaValues.push_back(iabseta*theZside);
         IphiValues.push_back(iphi);
         ICValues.push_back(thisIntercalibConstant);
         
         meanICforPhiRing += thisIntercalibConstant;
         totIphi++;

       }
       
       for ( int myPhiIndex = 0; myPhiIndex < totIphi; myPhiIndex++ )
         meanICforPhiRingValues.push_back(meanICforPhiRing/totIphi);


       for ( int iphi = MIN_IPHI; iphi <= MAX_IPHI; iphi++ ){

         int thisHashedIndex = GetHashedIndexEB(iabseta*theZside, iphi, theZside);
	 if ( h_Occupancy_hashedIndex -> GetBinContent(thisHashedIndex+1) == 0 ) continue;

         h_scale_EB_meanOnPhi -> Fill(iphi, iabseta*theZside, ICValues.at(myPhiIndex)/meanICforPhiRingValues.at(myPhiIndex));
         myPhiIndex++;

       }
       f1.SetParameters(histoAuxiliary.GetEntries(),histoAuxiliary.GetMean(),histoAuxiliary.GetRMS());
       f1.SetRange(1-5*histoAuxiliary.GetRMS(), 1+5*histoAuxiliary.GetRMS());
       histoAuxiliary . Fit("f1","QR");
       
       if ( f1.GetParError(2) > 0.5 ) continue;
       h_IntercalibSpread_iEta -> SetBinContent( iabseta*theZside + 85 + 1, f1.GetParameter(2) );
       h_IntercalibSpread_iEta -> SetBinError( iabseta*theZside + 85 + 1, f1.GetParError(2) );

     }
 
   }

}

void FastCalibrator::saveHistos(TFile * f1)
{

  f1->cd();
  hC_IntercalibValues -> Write(*f1);
  hC_PullFromScalib -> Write(*f1);
  hC_EoP -> Write(*f1);
  hC_scale_EB -> Write("",*f1);
  h_scalib_EB -> Write();

  h_IntercalibValues_test -> Write();
  h_Occupancy_hashedIndex -> Write();
  p_IntercalibValues_iEta -> Write();
  h_Init_IntercalibValues -> Write(); 

  h_IntercalibSpread_iEta -> Write();
  h_scale_EB -> Write();
  h_scale_EB_meanOnPhi -> Write("h_scale_map");
  h_scale_EB_hashedIndex -> Write(); 
  
  h_occupancy -> Write();
  
  g_ICmeanVsLoop -> Write();
  g_ICrmsVsLoop -> Write();


  f1->Close();

  return;
}

void FastCalibrator::printOnTxt(std::string outputTxtFile)
{
  std::ofstream outTxt (outputTxtFile.c_str(),std::ios::out);

  outTxt << "---------------------------------------------------------------" << std::endl;
  outTxt << "--- iEta ---- iPhi ------ IC value (mean on phi ring) ---------" << std::endl;
  outTxt << "---------------------------------------------------------------" << std::endl;

  for (unsigned int iIndex = 0; iIndex < ICValues.size(); iIndex++)
    outTxt << "  " << std::fixed << std::setw(4)  << std::right << IetaValues.at(iIndex) 
           << std::fixed << std::setw(10)  << std::right << IphiValues.at(iIndex) 
           << "          " << (float) ICValues.at(iIndex)/meanICforPhiRingValues.at(iIndex) << std::endl;

}

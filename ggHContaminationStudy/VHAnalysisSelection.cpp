#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TMath.h"

/// Selection on Leptons
#define PtLeadingCut  20.
#define PtTrailingCut 10.
#define EtaLeptonCut  2.5

/// Selection for cleaning
#define DeltaRJetLep  0.3
#define PtThresholdForCleaning  5.
#define JetNumber 5
#define JetTreshold 15.

/// Selection for jets
#define TagJetPtCut     30.
#define JetAcceptance   4.7
#define TagJetEtaCut    2.5
#define DEtaTagJetCut   2.1 
#define MjjTagJetCutMin 65.
#define MjjTagJetCutMax 115.
  
 
// Compile with g++ -Wall -o VHAnalysisSelection `root-config --glibs --cflags` VHAnalysisSelection.cpp

// Global Variables for Branch Settings
float f1V1pt, f1V1eta, f1V1phi, f1V1pdgId;
float f2V1pt, f2V1eta, f2V1phi, f2V1pdgId;
float f1V2pt, f1V2eta, f1V2phi, f1V2pdgId;
float f2V2pt, f2V2eta, f2V2phi, f2V2pdgId;

float jet1pt, jet1eta, jet1phi;
float jet2pt, jet2eta, jet2phi;
float jet3pt, jet3eta, jet3phi;
float jet4pt, jet4eta, jet4phi;
float jet5pt, jet5eta, jet5phi;


// External Functions
void InitializeTree (TTree * tree);
void CleaningJetLep (std::map<int,bool> *JetMap, float LeptonPt, float LeptonEta, float LeptonPhi);

//**************** Main Programme ********************//

int main (int argc, char *argv){

  // Input information 

  std::string TreeName         = "latino/latino";
  std::string BaseDirPath      = "/data2/amassiro/VBF/Summer12/23Nov2012/CMSSW_5_3_5/src/WWAnalysis/AnalysisStep/test/crab_mcDumper_53X/";
  std::string ListOfFile       = "ListOfFile.txt";
  std::string CrossSectionFile = "ggHCrossSection_Powheg8TeV_Summer12.txt"; 
  std::string TableOutputName  = "PowhegResult_VHSelection.txt";

  /// Cross Section Definition

  std::ifstream InputCrossSectionFile (CrossSectionFile.c_str());
  std::map <int,std::pair<float,std::pair<float,float> > > ggHCrossSectionMap ;
  std::string   buffer; 

  if(InputCrossSectionFile.fail()){ std::cerr<<" Fail to open Input XS File --> exit "<<std::endl; return 1; }

  std::cout<<" ########## Open Cross Section File ...... ######### "<<std::endl;
  std::cout<<std::endl;

  while(!InputCrossSectionFile.eof()){
  
    getline(InputCrossSectionFile,buffer);
    if(buffer=="" || !buffer.find("#") || buffer.size()==0) continue;
    std::stringstream line( buffer );      
    int mass ;
    line >> mass;                                     // Mass value
    line >> ggHCrossSectionMap[mass].first ;          // XS Value 
    line >> ggHCrossSectionMap[mass].second.first  ;  // XS Upper Error 
    line >> ggHCrossSectionMap[mass].second.second ;  // XS Lower Error
    std::cout<<" Higgs Mass "<<mass<<" XS = "<<ggHCrossSectionMap[mass].first<<" ErrorUp "<<ggHCrossSectionMap[mass].second.first
             <<" ErrorDown "<<ggHCrossSectionMap[mass].second.second << std::endl;
  }
  std::cout<<std::endl;
  if(ggHCrossSectionMap.size()==0){ std::cerr<<" No input xs, check the input XS file --> exit "<<std::endl; return 1; }
  InputCrossSectionFile.close();

  std::cout<<" ########## Open Cross Section File --> Ok ######### "<<std::endl;
  std::cout<<std::endl;
  // Make List of file to be executed for each WP
 
  std::string command     = "ls -l "+BaseDirPath+" | grep .root | grep ggH | awk '{print \""+BaseDirPath+"\"$9}'  > "+ListOfFile ;

  int status = system (command.c_str());
  if(status == -1) { std::cerr<< " Failing in build LisOfFile --> exit "<<std::endl; return 1;}
  
  // Open each File and Run the Analysis

  std::ifstream InputFileList (ListOfFile.c_str());
  std::vector <std::string> NameOfInputFile ;

  if(InputFileList.fail()){ std::cerr<<" Fail to open Input File List --> exit "<<std::endl; return 1; }

  std::cout<<" ########## Open Input File List ...... ######### "<<std::endl;  
  std::cout<<std::endl;

  while(!InputFileList.eof()){
  
    getline(InputFileList,buffer);
    if(buffer!="" && buffer.find("#") && buffer.size()!=0 ){ NameOfInputFile.push_back(buffer); std::cout<<buffer<<std::endl; }
  }

  if(NameOfInputFile.size()==0) { std::cerr<<" Empty input File List, size == 0; Check the Path --> exit"<<std::endl; return 1;}
  InputFileList.close();

  std::cout<<std::endl;
  std::cout<<" ########## Open Input File List --> Ok ######### "<<std::endl;  

  // OutputFile with format compatible with the mcfm code for the analysis
  std::ofstream OutputTable;
  command = " if [ -f "+TableOutputName+" ] ; then rm "+TableOutputName+" ; fi " ;
  system(command.c_str());

  
  OutputTable.open (TableOutputName.c_str(),std::ios::app);
  std::cout<<std::endl;
  std::cout << "Writing Output File to: " << TableOutputName << std::endl;
  std::cout << std::endl;

  std::cout << " ******************************************* start *******************************************" << std::endl;
  Long64_t start = clock();

  // Cycle on each input File --> One root file for each WP for xs evaluation
  
  std::map <int,std::pair<float,std::pair<float,float> > >::const_iterator  itMapggH = ggHCrossSectionMap.begin() ; 
 
  for( unsigned int iFile = 0; iFile <NameOfInputFile.size(); iFile ++){
 
    TFile *f = TFile::Open (NameOfInputFile.at(iFile).c_str()) ;
    if(f==0) { std::cerr<< " File not allocated ; skip "<<NameOfInputFile.at(iFile)<<std::endl; continue;}
    
    TTree * inputTree = (TTree*) f->Get(TreeName.c_str());
    if(inputTree==0) { std::cerr<< " Tree not allocated, wrong input name ; skip "<<std::endl; continue; }

    int NumberEntriesBefore = inputTree -> GetEntries();

    if(ggHCrossSectionMap.size()<(iFile+1)) { std::cerr<<" Lack of info in the Cross Section Map, Add info for: "<<NameOfInputFile.at(iFile) <<std::endl;
                                              continue; 
                                            }
                                   
    InitializeTree (inputTree);

    std::cout<<" Starting Analyze : "<<NameOfInputFile.at(iFile)<<std::endl;

    int NumberEntriesAfterSelection = 0 ; float efficiency =0. ; float efficiencyError = 0.;

    for( int iEntry = 0 ; iEntry < NumberEntriesBefore ; iEntry ++){

      if(iEntry%10000 == 0 ) std::cout << ">>> reading entry " << iEntry << " / " << inputTree->GetEntries() << "\r" << std::flush;

      inputTree -> GetEntry(iEntry);
   
      // Lepton Selection of the Analysis 

      float LeadingPt=-99., LeadingEta=-99., LeadingPhi=-99.;
      float TrailingPt=-99., TrailingEta=-99., TrailingPhi=-99.;
      
      // No ambiguity in the id, only two leptons in the genration at GenLevel and no difference between ele and mu

      // f1V1 = ele o mu , f1V2 = ele o mu
      if( (( fabs(f1V1pdgId) == 11 || fabs(f1V1pdgId) == 13 || fabs(f1V1pdgId) == 15)  && 
           ( fabs(f2V1pdgId) != 11 && fabs(f2V1pdgId) != 13 && fabs(f2V1pdgId) != 15)) &&
      	  (( fabs(f1V2pdgId) == 11 || fabs(f1V2pdgId) == 13 || fabs(f1V2pdgId) == 15)  &&
           ( fabs(f2V2pdgId) != 11 && fabs(f2V2pdgId) != 13 && fabs(f2V2pdgId) != 15)) ){

       if(f1V1pt >= f1V2pt){ LeadingPt=f1V1pt;   TrailingPt=f1V2pt; LeadingEta=f1V1eta; TrailingEta=f1V2eta; 
	                       LeadingPhi=f1V1phi; TrailingPhi=f1V2phi;
                             }

          else               { LeadingPt=f1V2pt;   TrailingPt=f1V1pt; LeadingEta=f1V2eta; TrailingEta=f1V1eta; 
	                       LeadingPhi=f1V2phi; TrailingPhi=f1V1phi;
                             }
      }
  
      // f2V1 = ele o mu , f1V2 = ele o mu
      if( (( fabs(f2V1pdgId) == 11 || fabs(f2V1pdgId) == 13 || fabs(f2V1pdgId) == 15)  &&
           ( fabs(f1V1pdgId) != 11 && fabs(f1V1pdgId) != 13 && fabs(f1V1pdgId) != 15)) &&
          (( fabs(f1V2pdgId) == 11 || fabs(f1V2pdgId) == 13 || fabs(f1V2pdgId) == 15)  &&
           ( fabs(f2V2pdgId) != 11 && fabs(f2V2pdgId) != 13 && fabs(f2V2pdgId) != 13)) ){
       
	if(f2V1pt >= f1V2pt){ LeadingPt=f2V1pt;   TrailingPt=f1V2pt; LeadingEta=f2V1eta; TrailingEta=f1V2eta; 
	                      LeadingPhi=f2V1phi; TrailingPhi=f1V2phi;
                            }
        
        else                { LeadingPt=f1V2pt;   TrailingPt=f2V1pt; LeadingEta=f1V2eta; TrailingEta=f2V1eta; 
	                      LeadingPhi=f1V2phi; TrailingPhi=f2V1phi;
                            }
      }

      // f1V1 = ele o mu , f2V2 = ele o mu
      if( (( fabs(f1V1pdgId) == 11 || fabs(f1V1pdgId) == 13 || fabs(f1V1pdgId) == 15)  &&
           ( fabs(f2V1pdgId) != 11 && fabs(f2V1pdgId) != 13 && fabs(f2V1pdgId) != 15)) &&
          (( fabs(f2V2pdgId) == 11 || fabs(f2V2pdgId) == 13 || fabs(f2V2pdgId) == 15 ) &&
           ( fabs(f1V2pdgId) != 11 && fabs(f1V2pdgId) != 13 && fabs(f1V2pdgId) != 15)) ){
       
	if(f1V1pt >= f2V2pt){ LeadingPt=f1V1pt;   TrailingPt=f2V2pt; LeadingEta=f1V1eta; TrailingEta=f2V2eta; 
	                      LeadingPhi=f1V1phi; TrailingPhi=f2V2phi;
                            }
        
        else                { LeadingPt=f2V2pt;   TrailingPt=f1V1pt; LeadingEta=f2V2eta; TrailingEta=f1V1eta; 
	                      LeadingPhi=f2V2phi; TrailingPhi=f1V1phi;
                            }
      }

      // f2V1 = ele o mu , f2V2 = ele o mu
      if( (( fabs(f2V1pdgId) == 11 || fabs(f2V1pdgId) == 13 || fabs(f2V1pdgId) == 15 )  &&
           ( fabs(f1V1pdgId) != 11 && fabs(f1V1pdgId) != 13 && fabs(f1V1pdgId) != 15 )) &&
          (( fabs(f2V2pdgId) == 11 || fabs(f2V2pdgId) == 13 || fabs(f2V2pdgId) == 15 )  && 
           ( fabs(f1V2pdgId) != 11 && fabs(f1V2pdgId) != 13 && fabs(f1V2pdgId) != 15 )) ){
       
	if(f2V1pt >= f2V2pt){ LeadingPt=f2V1pt;   TrailingPt=f2V2pt; LeadingEta=f2V1eta; TrailingEta=f2V2eta; 
	                      LeadingPhi=f2V1phi; TrailingPhi=f2V2phi;
                            }
        
        else                { LeadingPt=f2V2pt;   TrailingPt=f2V1pt; LeadingEta=f2V2eta; TrailingEta=f2V1eta; 
	                      LeadingPhi=f2V2phi; TrailingPhi=f2V1phi;
                            }
    


     }
   
     // map for jet lep cleaning
     std::map<int,bool> *JetLepCleaningMap = new std::map<int,bool> ();
     for(int iJet = 1; iJet<=JetNumber ; iJet++){ (*JetLepCleaningMap)[iJet]=true;}

     if(LeadingPt  > PtThresholdForCleaning  && fabs(LeadingEta)  < EtaLeptonCut ) CleaningJetLep(JetLepCleaningMap,LeadingPt,LeadingEta,LeadingPhi);
     if(TrailingPt > PtThresholdForCleaning  && fabs(TrailingEta) < EtaLeptonCut ) CleaningJetLep(JetLepCleaningMap,TrailingPt,TrailingEta,TrailingPhi);
     
     // At least two cleaned jet

     std::map <float,std::pair<float,float> > *SortedJetCleaned = new std::map <float,std::pair<float,float> > ();   

     for( unsigned int iJet = 0; iJet<JetLepCleaningMap->size() ; iJet++){

       if((*JetLepCleaningMap)[iJet]!=false){
	   
	 std::pair<float,float> *pairTemp                   = new std::pair<float,float> ();

	 if(iJet==1 && fabs(jet1eta) <= JetAcceptance ) { (*pairTemp).first = jet1eta ; (*pairTemp).second = jet1phi;
	                                                  (*SortedJetCleaned)[jet1pt] = *pairTemp ;
	                                                }
	 if(iJet==2 && fabs(jet2eta) <= JetAcceptance ) { (*pairTemp).first = jet2eta ; (*pairTemp).second = jet2phi;
	                                                  (*SortedJetCleaned)[jet2pt] = *pairTemp ;
                                                        }
         if(iJet==3 && fabs(jet3eta) <= JetAcceptance ) { (*pairTemp).first = jet3eta ; (*pairTemp).second = jet3phi;
	                                                  (*SortedJetCleaned)[jet3pt] = *pairTemp ;
                                                        }
         if(iJet==4 && fabs(jet4eta) <= JetAcceptance ) { (*pairTemp).first = jet4eta ; (*pairTemp).second = jet4phi;
	                                                  (*SortedJetCleaned)[jet4pt] = *pairTemp ;
                                                        }
         if(iJet==5 && fabs(jet5eta) <= JetAcceptance ) { (*pairTemp).first = jet5eta ; (*pairTemp).second = jet5phi;
	                                                  (*SortedJetCleaned)[jet5pt] = *pairTemp ;
                                                        }
	 delete pairTemp; 
	      
	 }
     }

     if(SortedJetCleaned->size()<2) continue; // skip event if less than two cleaned jet

     // Start with the analysis kinematical selections
    
     if( LeadingPt   < PtLeadingCut ) continue;
     if( TrailingPt  < PtTrailingCut) continue;
     if( fabs(LeadingEta)  > EtaLeptonCut ) continue;
     if( fabs(TrailingEta) > EtaLeptonCut ) continue;

     std::map <float,std::pair<float,float> >::const_reverse_iterator itMapLeadingJet = SortedJetCleaned->rbegin();
     std::map <float,std::pair<float,float> >::const_reverse_iterator itMapTrailingJet = SortedJetCleaned->rbegin(); itMapTrailingJet++;
         
     if(itMapLeadingJet->first  < TagJetPtCut || fabs(itMapLeadingJet->second.first)  > TagJetEtaCut ) continue;
     if(itMapTrailingJet->first < TagJetPtCut || fabs(itMapTrailingJet->second.first) > TagJetEtaCut ) continue;

     // Delta eta jj selection
     float DEtajj = 0.; float Mjj = 0.;
   
     DEtajj = fabs(itMapLeadingJet->second.first-itMapTrailingJet->second.first);
    
     if( DEtajj > DEtaTagJetCut ) continue ;
    
     // Invariant mass of dijet system in the hypothesis of particle without rest mass --> relativistic limit
     Mjj = sqrt(2.* itMapLeadingJet->first * itMapTrailingJet->first * (TMath::CosH(itMapLeadingJet->second.first-
	   itMapTrailingJet->second.first) - TMath::Cos(itMapLeadingJet->second.second - itMapTrailingJet->second.second)));

     if( Mjj > MjjTagJetCutMax || Mjj < MjjTagJetCutMin ) continue ;

     // Jet Veto Selection 
     bool doJetVeto = false ;
     itMapTrailingJet++;

     for( ;itMapTrailingJet != SortedJetCleaned->rend(); itMapTrailingJet++){
       if(itMapTrailingJet->first > TagJetPtCut ) doJetVeto = true;
     }
     if(doJetVeto==true) continue ;
     NumberEntriesAfterSelection = NumberEntriesAfterSelection +1 ;

     delete JetLepCleaningMap;
     delete SortedJetCleaned;
   }
    
   efficiency = float(NumberEntriesAfterSelection)/float(NumberEntriesBefore) ; 
   efficiencyError = 1./(float(NumberEntriesBefore))* sqrt(float(NumberEntriesAfterSelection));

   OutputTable <<"    "<< itMapggH->first                <<std::setprecision(3)
               <<"    "<< efficiency                     <<std::setprecision(7)       
               <<"    "<< efficiencyError                <<std::setprecision(7)
               <<"    "<< itMapggH->second.first         <<std::setprecision(3)
               <<"    "<< itMapggH->second.second.first  <<std::setprecision(4)
	       <<"    "<< itMapggH->second.second.second <<std::setprecision(4)
   << std::endl;
  
   itMapggH ++ ;   
  
   delete inputTree;
   delete f;
  
  }
  OutputTable.close() ;

  std::cout << " ******************************************* end *******************************************" << std::endl;
  Long64_t end = clock();
  std::cout <<"Time = " <<  ((double) (end - start)) << " (a.u.)" << std::endl;
    
  return 0;
  
}
    

void InitializeTree (TTree * tree){

  tree -> SetBranchStatus("*",0);

  tree -> SetBranchStatus("f1V1pt",1);      tree -> SetBranchAddress("f1V1pt",&f1V1pt);
  tree -> SetBranchStatus("f2V1pt",1);      tree -> SetBranchAddress("f2V1pt",&f2V1pt);
  tree -> SetBranchStatus("f1V2pt",1);      tree -> SetBranchAddress("f1V2pt",&f1V2pt);
  tree -> SetBranchStatus("f2V2pt",1);      tree -> SetBranchAddress("f2V2pt",&f2V2pt);

  tree -> SetBranchStatus("f1V1eta",1);     tree -> SetBranchAddress("f1V1eta",&f1V1eta);
  tree -> SetBranchStatus("f2V1eta",1);     tree -> SetBranchAddress("f2V1eta",&f2V1eta);
  tree -> SetBranchStatus("f1V2eta",1);     tree -> SetBranchAddress("f1V2eta",&f1V2eta);
  tree -> SetBranchStatus("f2V2eta",1);     tree -> SetBranchAddress("f2V2eta",&f2V2eta);

  tree -> SetBranchStatus("f1V1phi",1);     tree -> SetBranchAddress("f1V1phi",&f1V1phi);
  tree -> SetBranchStatus("f2V1phi",1);     tree -> SetBranchAddress("f2V1phi",&f2V1phi);
  tree -> SetBranchStatus("f1V2phi",1);     tree -> SetBranchAddress("f1V2phi",&f1V2phi);
  tree -> SetBranchStatus("f2V2phi",1);     tree -> SetBranchAddress("f2V2phi",&f2V2phi);

  tree -> SetBranchStatus("f1V1pdgId",1);   tree -> SetBranchAddress("f1V1pdgId",&f1V1pdgId);
  tree -> SetBranchStatus("f2V1pdgId",1);   tree -> SetBranchAddress("f2V1pdgId",&f2V1pdgId);
  tree -> SetBranchStatus("f1V2pdgId",1);   tree -> SetBranchAddress("f1V2pdgId",&f1V2pdgId);
  tree -> SetBranchStatus("f2V2pdgId",1);   tree -> SetBranchAddress("f2V2pdgId",&f2V2pdgId);

  tree -> SetBranchStatus("jet1pt",1);      tree -> SetBranchAddress("jet1pt",&jet1pt);
  tree -> SetBranchStatus("jet2pt",1);      tree -> SetBranchAddress("jet2pt",&jet2pt);
  tree -> SetBranchStatus("jet3pt",1);      tree -> SetBranchAddress("jet3pt",&jet3pt);
  tree -> SetBranchStatus("jet4pt",1);      tree -> SetBranchAddress("jet4pt",&jet4pt);
  tree -> SetBranchStatus("jet5pt",1);      tree -> SetBranchAddress("jet5pt",&jet5pt);

  tree -> SetBranchStatus("jet1eta",1);     tree -> SetBranchAddress("jet1eta",&jet1eta);
  tree -> SetBranchStatus("jet2eta",1);     tree -> SetBranchAddress("jet2eta",&jet2eta);
  tree -> SetBranchStatus("jet3eta",1);     tree -> SetBranchAddress("jet3eta",&jet3eta);
  tree -> SetBranchStatus("jet4eta",1);     tree -> SetBranchAddress("jet4eta",&jet4eta);
  tree -> SetBranchStatus("jet5eta",1);     tree -> SetBranchAddress("jet5eta",&jet5eta);

  tree -> SetBranchStatus("jet1phi",1);     tree -> SetBranchAddress("jet1phi",&jet1phi);
  tree -> SetBranchStatus("jet2phi",1);     tree -> SetBranchAddress("jet2phi",&jet2phi);
  tree -> SetBranchStatus("jet3phi",1);     tree -> SetBranchAddress("jet3phi",&jet3phi);
  tree -> SetBranchStatus("jet4phi",1);     tree -> SetBranchAddress("jet4phi",&jet4phi);
  tree -> SetBranchStatus("jet5phi",1);     tree -> SetBranchAddress("jet5phi",&jet5phi);

  return;

}


void CleaningJetLep (std::map<int,bool> *JetMap, float LeptonPt, float LeptonEta, float LeptonPhi){
  
  for(unsigned int iJet = 1 ;  iJet<=JetMap->size() ; iJet++){

    if((iJet ==1 && jet1pt > JetTreshold && sqrt(TMath::Power((jet1eta-LeptonEta),2)+TMath::Power((jet1phi-LeptonPhi),2))< DeltaRJetLep) || jet1pt<0) 
      (*JetMap)[iJet]=false;

    if((iJet ==2 && jet2pt > JetTreshold && sqrt(TMath::Power((jet2eta-LeptonEta),2)+TMath::Power((jet2phi-LeptonPhi),2))< DeltaRJetLep) || jet2pt<0) 
      (*JetMap)[iJet]=false;

    if((iJet ==3 && jet3pt > JetTreshold && sqrt(TMath::Power((jet3eta-LeptonEta),2)+TMath::Power((jet3phi-LeptonPhi),2))< DeltaRJetLep) || jet3pt<0)
      (*JetMap)[iJet]=false;

    if((iJet ==4 && jet4pt > JetTreshold && sqrt(TMath::Power((jet4eta-LeptonEta),2)+TMath::Power((jet4phi-LeptonPhi),2))< DeltaRJetLep) || jet4pt<0) 
      (*JetMap)[iJet]=false;

    if((iJet ==5 && jet5pt > JetTreshold && sqrt(TMath::Power((jet5eta-LeptonEta),2)+TMath::Power((jet5phi-LeptonPhi),2))< DeltaRJetLep) || jet5pt<0) 
      (*JetMap)[iJet]=false;

  }
  
  return;
}

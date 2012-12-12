// Code for mcfm powheg result comparison
// Complie with g++ -Wall `root-config --libs --glibs --cflags` -o ComparisonAnalysis_VH ComparisonAnalysis_VH.cpp

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>

#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"

#define NumberOfBranchingRatio 9
#define HistoScaleMax    30.
#define HistoScaleMin    0.
#define HistoNbin        1000
#define FractHiggsMass4CentralValue 1.
#define HiggsMassFractionMin 0.25
#define ScaleRatioMin 0.5
#define ScaleRatioMax 2

// PDF Error : set values for each desidered Higgs mass from yellow report
  /*
  std::vector< std::pair< float, float> > PDF_Error;
  PDF_Error.push_back(std::pair<int,int> (7.6,7.0));
  PDF_Error.push_back(std::pair<int,int> (7.6,7.2));
  PDF_Error.push_back(std::pair<int,int> (7.6,7.3));
  PDF_Error.push_back(std::pair<int,int> (7.6,7.5));
  PDF_Error.push_back(std::pair<int,int> (7.5,7.6));
  PDF_Error.push_back(std::pair<int,int> (7.5,7.8));
  PDF_Error.push_back(std::pair<int,int> (7.5,7.8));
  PDF_Error.push_back(std::pair<int,int> (7.5,7.8));
  PDF_Error.push_back(std::pair<int,int> (7.6,7.8));
  PDF_Error.push_back(std::pair<int,int> (7.5,7.9));
  PDF_Error.push_back(std::pair<int,int> (7.6,7.9));
  PDF_Error.push_back(std::pair<int,int> (7.7,8.0));
  PDF_Error.push_back(std::pair<int,int> (7.8,8.1));
  PDF_Error.push_back(std::pair<int,int> (8.0,8.3));
  */

// Basic struct to identify a MCFM job output
struct Sample { float Higgs_Mass;
                float Fact_Scale;
                float Rin_Scale;
                int   Run;
              };

// Information of cross section for each group of mcfm job 
struct Sample_xs { float Higgs_Mass;
                   float Fact_Scale;
                   float Rin_Scale;
                   TH1F* XSec_Virt;
                   TH1F* XSec_Real;
                   float err_prop_Virt;
                   float err_prop_Real;
                   float xs_value_Run;
                   float xs_error_Run;
                   float xs_error_RMS_Run;
           
                   bool operator==(const Sample_xs & a){
                    if(Higgs_Mass == a.Higgs_Mass && Fact_Scale==a.Fact_Scale && Rin_Scale==a.Rin_Scale)
                    return true;
                    else return false;
                   }

                   bool operator!=(const Sample_xs & a){
                   if(Higgs_Mass != a.Higgs_Mass || Fact_Scale!=a.Fact_Scale || Rin_Scale!=a.Rin_Scale)
                   return true;
                   else return false;
                  }
               
                  void operator=(const Sample_xs & a) {
                   Higgs_Mass       = a.Higgs_Mass;
                   Fact_Scale       = a.Fact_Scale;
                   Rin_Scale        = a.Rin_Scale;
                   err_prop_Virt    = a.err_prop_Virt;
                   err_prop_Real    = a.err_prop_Real;
                   XSec_Virt        = (TH1F*) a.XSec_Virt->Clone(XSec_Virt->GetName());
                   XSec_Real        = (TH1F*) a.XSec_Real->Clone(XSec_Real->GetName());
                   xs_value_Run     = a.xs_value_Run;
                   xs_error_Run     = a.xs_error_Run;
                   xs_error_RMS_Run = a.xs_error_RMS_Run;
               
                 }
             };

// Struct useful to calculate mcfm Cross section scale dependence
struct Sample_Scale { float Higgs_Mass;
                      float Fact_Scale_Max;
                      float Fact_Scale_Min;
                      float Rin_Scale_Max;
                      float Rin_Scale_Min;

                      bool operator==(const Sample_Scale & a){
                        if(Higgs_Mass == a.Higgs_Mass)
                        return true;
                        else return false;
                      }
                      
                      bool operator!=(const Sample_Scale & a){
                          if(Higgs_Mass != a.Higgs_Mass)
                          return true;
                          else return false;
                      }
               
                      void operator=(const Sample_Scale & a){
                        Higgs_Mass = a.Higgs_Mass;
                        Fact_Scale_Max = a.Fact_Scale_Max;
                        Fact_Scale_Min  = a.Fact_Scale_Min;
                        Rin_Scale_Max = a.Rin_Scale_Max;
                        Rin_Scale_Min = a.Rin_Scale_Min;
                      }

                   };

// Struct useful to store powheg information
struct Sample_powheg { float HiggsMass;
                       float efficiency;
                       float err_efficiency;
                       float xs;
                       float err_xs_up;
                       float err_xs_down;
                       float xs_eff;
                       float err_xs_eff_up;
                       float err_xs_eff_down;
                     };

// Set of Functions used in the code
float GetMaximum_RS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass);

float GetMinimum_RS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass);

float GetMaximum_FS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass);

float GetMinimum_FS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass);

// Main Function             

int main (int argc, char** argv){

  //  if(argc!=3){ std::cerr<<" Wrong number of parameter passed on the command line --> Exit "<<std::endl; return 1;}

  std::string VirtualFileDirectory = "VHAnalysis_JOB/VHWW_Virtual_NEW/" ;
  std::string RealFileDirectory    = "VHAnalysis_JOB/VHWW_Real_NEW/" ;
  std::string ListOfVirtualFile    = "VHWW_VirtualList.txt" ;
  std::string ListOfRealFile       = "VHWW_RealList.txt" ;
  std::string VirtualFileNameBase  = "HWW2jt_virt_mstw8nl";
  std::string RealFileNameBase     = "HWW2jt_real_mstw8nl";

  std::string OutputROOTFileName   = argv[1] ;
  std::string PowhegInputFile      = argv[2] ;
  std::string OutputDumpFile       = argv[3] ;

  // Create the output File

  TFile *outputFile = new TFile(OutputROOTFileName.c_str(),"RECREATE");
  outputFile->cd();

  // Creaion of virtual and  real mcfm job file list
 
  std::cout<<" ############## Generating Virtual and Real file List ############### "<<std::endl;
  std::cout<<"  "<<std::endl;

  std::string command = "if [ -f "+ListOfVirtualFile+" ] ; then rm "+ListOfVirtualFile+" ; fi ";

  system(command.c_str());
 
  command = "ls "+VirtualFileDirectory+" | grep virt | grep .dat  | sed -e s%vhHww%%g  | sed -e s%dat%%g | tr \"_\" \" \" "+
            " | awk '{print \"echo @vhHww_\"$6\"_\"$7\"_\"$8\"_\"$9\"_\"$10\" \"$7\" \"$8\" \"$9\" \"$10\" @ >> "+ListOfVirtualFile+" \"}'"+
            " | tr \"@\" \"\\\"\" | /bin/sh" ;

  int staus = system(command.c_str());
  if(staus==-1){ std::cerr<< " Error in creation of List of Virtual jobs --> Exit "<<std::endl; return 1; }


  command = "if [ -f "+ListOfRealFile+" ] ; then rm "+ListOfRealFile+" ; fi ";

  system(command.c_str());

  command = "ls "+RealFileDirectory+" | grep real | grep .dat  | sed -e s%vhHww%%g  | sed -e s%dat%%g | tr \"_\" \" \" "+
            " | awk '{print \"echo @vhHww_\"$6\"_\"$7\"_\"$8\"_\"$9\"_\"$10\" \"$7\" \"$8\" \"$9\" \"$10\" @ >> "+ListOfRealFile+" \"}'"+
            " | tr \"@\" \"\\\"\" | /bin/sh" ;
 
  staus = system(command.c_str());
  if(staus==-1){ std::cerr<< " Error in creation of List of Virtual jobs --> Exit "<<std::endl; return 1; }
  std::cout<<" ############## Virtual and Real file List ---> Ok ############### "<<std::endl;
  std::cout<<"  "<<std::endl;

  // Acquisition from Virtual file list  --> build virtual map

  std::cout<<" ############## Fill Virtual File Map  ############## "<<std::endl;
  std::cout<<"  "<<std::endl;

  std::ifstream FileVirtual(ListOfVirtualFile.c_str());
  if(FileVirtual.fail()){ std::cerr<<" Fail to Open Virtual File List --> Check it "<<std::endl; return 1;}

  std::string buffer;
  std::map <std::string, Sample > map_Virt;
  
  while (!FileVirtual.eof()){  
  
     getline(FileVirtual,buffer);
     Sample* VirtualSampleInfo = new Sample ();
     std::string name_temp;

     if (buffer == "" || !buffer.find('#')) continue;
     std::stringstream line( buffer );
     line >> name_temp;
     line >> VirtualSampleInfo->Higgs_Mass;
     line >> VirtualSampleInfo->Fact_Scale;
     line >> VirtualSampleInfo->Rin_Scale;
     line >> VirtualSampleInfo->Run;
     map_Virt[name_temp]= *(VirtualSampleInfo);
  
     delete VirtualSampleInfo ;
  }
  
  FileVirtual.close();
  if(map_Virt.size()==0){ std::cerr<<" No virtual sample information --> Exit from the code "<<std::endl; return 1;}

  std::cout<<" ############## Fill Virtual Map --> Ok  ############## "<<std::endl;
  std::cout<<"  "<<std::endl;

  // Acquisition from Real file list  --> build Real map

  std::cout<<" ############## Fill Real File Map  ############## "<<std::endl;
  std::cout<<"  "<<std::endl;

 
  std::ifstream FileReal(ListOfRealFile.c_str());
  if(FileReal.fail()){ std::cerr<<" Fail to Open Real File List --> Check it "<<std::endl; return 1;}

  std::map <std::string, Sample > map_Real;
  
  while (!FileReal.eof()){

     getline(FileReal,buffer);
     Sample* RealSampleInfo = new Sample ();
     std::string name_temp;
  
     if (buffer == "" || !buffer.find('#')) continue;

     std::stringstream line( buffer );
     line >> name_temp;
     line >> RealSampleInfo->Higgs_Mass;
     line >> RealSampleInfo->Fact_Scale;
     line >> RealSampleInfo->Rin_Scale;
     line >> RealSampleInfo->Run;
     map_Real[name_temp]= *(RealSampleInfo);

     delete RealSampleInfo ;     
  }

  FileReal.close();
  if(map_Real.size()==0){ std::cerr<<" No real sample information --> Exit from the code "<<std::endl; return 1;}

  std::cout<<" ############## Fill Real Map --> Ok  ############## "<<std::endl;
  std::cout<<"  "<<std::endl;

  // Fill virtual cross section value for each run
  std::cout<<" ############## Build Information and table for Virtual Result  ############## "<<std::endl;
  std::cout<<"  "<<std::endl;

  std::vector<Sample_xs> Vector_xs_sample ;

  for(std::map <std::string, Sample >::const_iterator itMap = map_Virt.begin(); itMap != map_Virt.end() ; itMap++){
    
    TString nameFileIn ;
    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale<100. )
      nameFileIn = Form("%s/%s_%d__%d__%sdat",VirtualFileDirectory.c_str(),VirtualFileNameBase.c_str(),int(itMap->second.Rin_Scale),
			int(itMap->second.Fact_Scale),itMap->first.c_str());
    
    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale>=100. )
      nameFileIn = Form("%s/%s_%d_%d_%sdat",VirtualFileDirectory.c_str(),VirtualFileNameBase.c_str(),int(itMap->second.Rin_Scale),
			int(itMap->second.Fact_Scale),itMap->first.c_str());

    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale>=100. )
      nameFileIn = Form("%s/%s_%d_%d__%sdat",VirtualFileDirectory.c_str(),VirtualFileNameBase.c_str(),int(itMap->second.Rin_Scale),
			int(itMap->second.Fact_Scale),itMap->first.c_str());
    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale<100. )
      nameFileIn = Form("%s/%s_%d__%d_%sdat",VirtualFileDirectory.c_str(),VirtualFileNameBase.c_str(),int(itMap->second.Rin_Scale),
			int(itMap->second.Fact_Scale),itMap->first.c_str());

    std::ifstream FileIn (nameFileIn.Data());
    if(FileIn.fail()){ std::cerr<<" Wrong file virt: "<<nameFileIn<<" check the name ---> Skip "<<std::endl; continue;}

    // Read xs value from mcfm .dat output file
    std::string tempString;
    float xsec, errxsec;
 
    getline(FileIn,buffer);
    std::stringstream line( buffer );
    line>>tempString>>tempString;
      
    if(tempString == "Intermediate"){
    getline(FileIn,buffer);
    getline(FileIn,buffer);
    getline(FileIn,buffer);
    }
    else{
       getline(FileIn,buffer);
       getline(FileIn,buffer);
     }
    
    Sample_xs *var = new Sample_xs();
    var->Higgs_Mass = itMap->second.Higgs_Mass;
    var->Fact_Scale = itMap->second.Fact_Scale;
    var->Rin_Scale  = itMap->second.Rin_Scale;
 
    line << buffer ;
    line >> tempString >> tempString;
    line >> xsec;
    line >> tempString;
    line >> errxsec;

    std::vector<Sample_xs>::iterator pos = find(Vector_xs_sample.begin(),Vector_xs_sample.end(),*(var));
    
    if(pos==Vector_xs_sample.end() || (*pos).XSec_Virt==0) {
                 
        var->XSec_Virt = new TH1F (itMap->first.c_str(),itMap->first.c_str(),HistoNbin,HistoScaleMin,HistoScaleMax);
        std::string name_tmp = itMap->first+"_Real";
        var->XSec_Real = new TH1F (name_tmp.c_str(),name_tmp.c_str(),HistoNbin,-fabs(HistoScaleMax),-fabs(HistoScaleMin));
        var->XSec_Virt->Fill(xsec);
        var->err_prop_Virt = errxsec*errxsec;
        Vector_xs_sample.push_back(*(var));
    }
    else{
          (*pos).XSec_Virt->Fill(xsec);
          (*pos).err_prop_Virt=(*pos).err_prop_Virt+errxsec*errxsec;
	}

    delete var;
    FileIn.close();
  }
 
  // Fill real cross section value for each run
  std::cout<<" ############## Build Information and table for Real Result  ############## "<<std::endl;
  std::cout<<"  "<<std::endl;
  

  for(std::map <std::string, Sample >::const_iterator itMap = map_Real.begin(); itMap != map_Real.end() ; itMap++){
    
    TString nameFileIn;
    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale<100. )
      nameFileIn = Form("%s/%s_%d__%d__%sdat",RealFileDirectory.c_str(),RealFileNameBase.c_str(),
                         int(itMap->second.Rin_Scale),int(itMap->second.Fact_Scale),itMap->first.c_str());

    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale>=100. )
      nameFileIn = Form("%s/%s_%d_%d_%sdat",RealFileDirectory.c_str(),RealFileNameBase.c_str(),
                         int(itMap->second.Rin_Scale),int(itMap->second.Fact_Scale),itMap->first.c_str());
    
    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale>=100. )
      nameFileIn = Form("%s/%s_%d_%d__%sdat",RealFileDirectory.c_str(),RealFileNameBase.c_str(),
                         int(itMap->second.Rin_Scale),int(itMap->second.Fact_Scale),itMap->first.c_str());

    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale<100. )
      nameFileIn = Form("%s/%s_%d__%d_%sdat",RealFileDirectory.c_str(),RealFileNameBase.c_str(),
                         int(itMap->second.Rin_Scale),int(itMap->second.Fact_Scale),itMap->first.c_str());
 
     std::ifstream FileIn (nameFileIn.Data());
     if(FileIn.fail()){ std::cerr<<" Wrong file real: "<<nameFileIn<<" check the name ---> Skip "<<std::endl; continue;}
  
     std::string tempString;
     float xsec, errxsec;
 
     getline(FileIn,buffer);
     std::stringstream line( buffer );
     line>>tempString>>tempString;

     if(tempString == "Intermediate"){
       getline(FileIn,buffer);
       getline(FileIn,buffer);
       getline(FileIn,buffer);
     }
     else{
          getline(FileIn,buffer);
          getline(FileIn,buffer);
     }

    Sample_xs *var = new Sample_xs();
    var->Higgs_Mass = itMap->second.Higgs_Mass;
    var->Fact_Scale = itMap->second.Fact_Scale;
    var->Rin_Scale  = itMap->second.Rin_Scale;
 
    line << buffer ;
    line >> tempString >> tempString;
    line >> xsec;
    line >> tempString;
    line >> errxsec;

    //    std::cout<<" xsec "<<xsec<<" errxsec "<<errxsec<<std::endl; 
    std::vector<Sample_xs>::iterator pos = find(Vector_xs_sample.begin(),Vector_xs_sample.end(),*(var));
    (*pos).XSec_Real->Fill(xsec);
    (*pos).err_prop_Real=(*pos).err_prop_Real+errxsec*errxsec;

    delete var;
    FileIn.close();
      
  }

  std::ofstream OutputDump;
  command = " if [ -f "+OutputDumpFile+" ] ; then rm "+OutputDumpFile+" ; fi " ;
  system(command.c_str());

  OutputDump.open (OutputDumpFile.c_str(),std::ios::app);

  
  OutputDump<<"--------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  OutputDump<<"------------------------------------  Value obtained for each scale and Higgs Mass  --------------------------------------------"<<std::endl;
  OutputDump<<"--------------------------------------------------------------------------------------------------------------------------------"<<std::endl;

  OutputDump<<" Mass -- Fact Scale -- Rin Scale -- XSec Virt -- Err XSec -- Err XSec Prop -- XSec Real -- Err XSec -- Err XSec Propagation ----"<<std::endl;
  OutputDump<<"-------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  

 
  for( std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){
   
    OutputDump<<(*itVec).Higgs_Mass                                        <<"        "<<std::setprecision(4)<<
                (*itVec).Fact_Scale                                        <<"        "<<std::setprecision(4)<<
                (*itVec).Rin_Scale                                         <<"        "<<std::setprecision(4)<<
                (*itVec).XSec_Virt->GetMean()*NumberOfBranchingRatio       <<"        "<<std::setprecision(5)<<
                (*itVec).XSec_Virt->GetMeanError()*NumberOfBranchingRatio  <<"        "<<std::setprecision(5)<<
                sqrt((*itVec).err_prop_Virt/(*itVec).XSec_Virt->GetEntries())*NumberOfBranchingRatio <<"       "<<std::setprecision(5)<<
                (*itVec).XSec_Real->GetMean()*NumberOfBranchingRatio       <<"        "<<std::setprecision(5)<<
                (*itVec).XSec_Real->GetMeanError()*NumberOfBranchingRatio  << "       "<<std::setprecision(5)<<
                sqrt((*itVec).err_prop_Real/(*itVec).XSec_Real->GetEntries())*NumberOfBranchingRatio <<std::setprecision(5)<<std::endl;
  }

  OutputDump<<std::endl;
  OutputDump<<"--------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  OutputDump<<"--------------------------------  Value obtained from mean for each scale and Higgs Mass  --------------------------------------"<<std::endl;
  OutputDump<<"--------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  OutputDump<<"  Mass -- Fact Scale -- Rin Scale -- XSec -- Err XSec---Err XSec Propagation                                                    "<<std::endl;
  OutputDump<<"--------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  

  std::vector<float> Higgs_Mass_Vect;

  for( std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){

   if(Higgs_Mass_Vect.size()==0) Higgs_Mass_Vect.push_back((*itVec).Higgs_Mass);
   else{

         std::vector<float>::iterator pos = find (Higgs_Mass_Vect.begin(),Higgs_Mass_Vect.end(),(*itVec).Higgs_Mass);
         if(pos == Higgs_Mass_Vect.end()) Higgs_Mass_Vect.push_back((*itVec).Higgs_Mass);                            
   }
  
  
   (*itVec).xs_value_Run       = (*itVec).XSec_Virt->GetMean()+(*itVec).XSec_Real->GetMean();

   (*itVec).xs_error_Run       = sqrt(TMath::Power((*itVec).err_prop_Virt/(*itVec).XSec_Virt->GetEntries(),2)+
                                      TMath::Power((*itVec).err_prop_Real/(*itVec).XSec_Real->GetEntries(),2));

   (*itVec).xs_error_RMS_Run   =  sqrt((*itVec).XSec_Virt->GetMeanError()*(*itVec).XSec_Virt->GetMeanError()+
                                       (*itVec).XSec_Real->GetMeanError()*(*itVec).XSec_Real->GetMeanError());

   OutputDump<< (*itVec).Higgs_Mass                               <<"        "<<std::setprecision(5)<< 
               (*itVec).Fact_Scale                               <<"        "<<std::setprecision(5)<<
               (*itVec).Rin_Scale                                <<"        "<<std::setprecision(5)<<
               (*itVec).xs_value_Run*NumberOfBranchingRatio      <<"        "<<std::setprecision(5)<<
               (*itVec).xs_error_RMS_Run*NumberOfBranchingRatio  <<"        "<<std::setprecision(5)<<
               (*itVec).xs_error_Run*NumberOfBranchingRatio      <<std::setprecision(5)<<std::endl;
  }

    
  // Save Real and virtual xs distribution + Error in the output file
  for(unsigned int i=0; i<Vector_xs_sample.size(); i++){  

       Vector_xs_sample.at(i).XSec_Virt->Write();
       Vector_xs_sample.at(i).XSec_Real->Write();
  }
   

  /// Scale Uncertainty Study

  std::vector<Sample_Scale> map_Scale;

  for( std::vector<float>::iterator itVec =Higgs_Mass_Vect.begin();  itVec !=Higgs_Mass_Vect.end(); ++itVec){

      float max_RS = GetMaximum_RS(Vector_xs_sample,*itVec);
      float min_RS = GetMinimum_RS(Vector_xs_sample,*itVec);
      float max_FS = GetMaximum_FS(Vector_xs_sample,*itVec);
      float min_FS = GetMinimum_FS(Vector_xs_sample,*itVec);

      Sample_Scale *var = new Sample_Scale();
      var->Higgs_Mass     = (*itVec);
      var->Fact_Scale_Max = max_FS; 
      var->Fact_Scale_Min = min_FS;
      var->Rin_Scale_Max  = max_RS; 
      var->Rin_Scale_Min  = min_RS;
      map_Scale.push_back(*(var));
      delete var;
  }
 
  std::map<float,std::vector<float> > XSec_Scale, XSec_Scale_Virt, XSec_Scale_Real;

  for(std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){
   
      Sample_Scale* var = new Sample_Scale();
      var->Higgs_Mass = (*itVec).Higgs_Mass;

      std::vector<Sample_Scale>::iterator pos = find (map_Scale.begin(),map_Scale.end(),*(var));
      delete var ;

      // Skip extreme scale values --> recipe from yellow report
      //  if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
      //  if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      //  if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      //  if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
      
      std::map<float,std::vector<float> >::iterator itMap = XSec_Scale.find((*itVec).Higgs_Mass);
      if(itMap == XSec_Scale.end()) XSec_Scale[(*itVec).Higgs_Mass].assign(3,0);
     
      if((*itVec).Fact_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue && (*itVec).Rin_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue)
       XSec_Scale[(*itVec).Higgs_Mass].at(0)=(*itVec).xs_value_Run;
    
      std::map<float,std::vector<float> >::iterator itMap_Virt = XSec_Scale_Virt.find((*itVec).Higgs_Mass);
      if(itMap_Virt == XSec_Scale_Virt.end()) XSec_Scale_Virt[(*itVec).Higgs_Mass].assign(3,0);
     
      if((*itVec).Fact_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue && (*itVec).Rin_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue)
       XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0)=(*itVec).XSec_Virt->GetMean();

      std::map<float,std::vector<float> >::iterator itMap_Real = XSec_Scale_Real.find((*itVec).Higgs_Mass);
      if(itMap_Real == XSec_Scale_Real.end()) XSec_Scale_Real[(*itVec).Higgs_Mass].assign(3,0);
     
      if((*itVec).Fact_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue && (*itVec).Rin_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue)
       XSec_Scale_Real[(*itVec).Higgs_Mass].at(0)=(*itVec).XSec_Real->GetMean();
       
  }

  for(std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){
   
      Sample_Scale *var = new Sample_Scale () ;
      var->Higgs_Mass = (*itVec).Higgs_Mass;
      std::vector<Sample_Scale>::iterator pos = find (map_Scale.begin(),map_Scale.end(),*(var));
      delete var;
      
      // Skip extreme scale values --> recipe from yellow report
      // if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
      // if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      // if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      // if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
     
     
      if((*itVec).Fact_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue && (*itVec).Rin_Scale ==(*itVec).Higgs_Mass*FractHiggsMass4CentralValue) continue;
      if(((*itVec).Fact_Scale/(*itVec).Rin_Scale)<ScaleRatioMin || ((*itVec).Fact_Scale/(*itVec).Rin_Scale)> ScaleRatioMax) continue ;
      if((*itVec).Fact_Scale ==(*itVec).Higgs_Mass*HiggsMassFractionMin || (*itVec).Rin_Scale==(*itVec).Higgs_Mass*HiggsMassFractionMin) continue;

      if(XSec_Scale[(*itVec).Higgs_Mass].at(1) < ((*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0)))
         XSec_Scale[(*itVec).Higgs_Mass].at(1) =  (*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0);
       
      if(XSec_Scale[(*itVec).Higgs_Mass].at(2) > ((*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0))) 
         XSec_Scale[(*itVec).Higgs_Mass].at(2) =  (*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0);
            
      if(XSec_Scale_Virt[(*itVec).Higgs_Mass].at(1) < ((*itVec).XSec_Virt->GetMean()-XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0))) 
         XSec_Scale_Virt[(*itVec).Higgs_Mass].at(1) =  (*itVec).XSec_Virt->GetMean()-XSec_Scale[(*itVec).Higgs_Mass].at(0);
      
      if(XSec_Scale_Virt[(*itVec).Higgs_Mass].at(2) > ((*itVec).XSec_Virt->GetMean()-XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0))) 
         XSec_Scale_Virt[(*itVec).Higgs_Mass].at(2) =  (*itVec).XSec_Virt->GetMean()-XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0);
      
      if(XSec_Scale_Real[(*itVec).Higgs_Mass].at(1) < ((*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0))) 
         XSec_Scale_Real[(*itVec).Higgs_Mass].at(1) =  (*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0);
         
      if(XSec_Scale_Real[(*itVec).Higgs_Mass].at(2) > ((*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0))) 
         XSec_Scale_Real[(*itVec).Higgs_Mass].at(2) =  (*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0);
      
  }
         
  // ****************************************************************************
  // ********************* PLOT MCFM RESULT *************************************
  // ****************************************************************************

  TGraphAsymmErrors* mcfm_xs = new TGraphAsymmErrors();
  mcfm_xs->SetMarkerStyle(20);
  mcfm_xs->SetMarkerSize(1);
  mcfm_xs->SetMarkerColor(kBlue+2);
  mcfm_xs->SetLineColor(kBlue+2);
  mcfm_xs->SetLineWidth(2);

  TGraphAsymmErrors* mcfm_xs_stat = new TGraphAsymmErrors();
  mcfm_xs_stat->SetMarkerStyle(20);
  mcfm_xs_stat->SetMarkerSize(1);
  mcfm_xs_stat->SetMarkerColor(kGreen);
  mcfm_xs_stat->SetLineColor(kGreen);
  mcfm_xs_stat->SetLineWidth(2);

  /*  TGraphAsymmErrors* mcfm_xs_pdf = new TGraphAsymmErrors();
  mcfm_xs_pdf->SetMarkerStyle(20);
  mcfm_xs_pdf->SetMarkerSize(1);
  mcfm_xs_pdf->SetMarkerColor(kYellow);
  mcfm_xs_pdf->SetLineColor(kYellow);
  mcfm_xs_pdf->SetLineWidth(2);
  */

  OutputDump<<"----------------------------------------------------------------------------"<<std::endl;
  OutputDump<<"------------  Scale dependence for each Higgs Mass evaluation  -------------"<<std::endl;
  OutputDump<<"----------------------------------------------------------------------------"<<std::endl;
  OutputDump<<"Mass ---------- XSec --------- Err XSec Up -------- Err XSec Down"<<std::endl;
  OutputDump<<"----------------------------------------------------------------------------"<<std::endl;
  
  int iPoint =0;

  for(std::map<float,std::vector<float> >::iterator itVec =XSec_Scale.begin();  itVec !=XSec_Scale.end(); ++itVec){
  
    if(itVec->second.at(1)>0 && itVec->second.at(2)<0 ){
      OutputDump<<itVec->first                                     <<"         "<<std::setprecision(4)
               <<itVec->second.at(0)*NumberOfBranchingRatio       <<"         "<<std::setprecision(5)
               <<fabs(itVec->second.at(1))*NumberOfBranchingRatio <<"         "<<std::setprecision(5)
               <<fabs(itVec->second.at(2))*NumberOfBranchingRatio<< std::setprecision<<std::endl;

      mcfm_xs->SetPoint(iPoint,itVec->first,itVec->second.at(0)*NumberOfBranchingRatio);
      mcfm_xs->SetPointError(iPoint,0.,0.,fabs(itVec->second.at(2))*NumberOfBranchingRatio,fabs(itVec->second.at(1))*NumberOfBranchingRatio);
      iPoint++; continue;
    }
  
  }

  // Stat and PDF Error
  iPoint=0;
  for( std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec) {
 
    if((*itVec).Higgs_Mass!=(*itVec).Fact_Scale*FractHiggsMass4CentralValue || (*itVec).Higgs_Mass!=(*itVec).Rin_Scale*FractHiggsMass4CentralValue) continue;

    mcfm_xs_stat->SetPoint(iPoint+1,(*itVec).Higgs_Mass,((*itVec).XSec_Virt->GetMean()+(*itVec).XSec_Real->GetMean())*NumberOfBranchingRatio);

    mcfm_xs_stat->SetPointError(iPoint+1,0.,0.,sqrt(TMath::Power((*itVec).err_prop_Virt/(*itVec).XSec_Virt->GetEntries(),2)+
                                                    TMath::Power((*itVec).err_prop_Real/(*itVec).XSec_Real->GetEntries(),2))*NumberOfBranchingRatio,
				               sqrt(TMath::Power((*itVec).err_prop_Virt/(*itVec).XSec_Virt->GetEntries(),2)+
                                                     TMath::Power((*itVec).err_prop_Real/(*itVec).XSec_Real->GetEntries(),2))*NumberOfBranchingRatio);
    iPoint++;
  }

  /*
  iPoint=0;
  for(std::map<float,std::vector<float> >::iterator itVec =XSec_Scale.begin();  itVec !=XSec_Scale.end(); ++itVec){
      
	mcfm_xs_pdf->SetPoint(iPoint+1,itVec->first,itVec->second.at(0)*NumberOfBranchingRatio);
	mcfm_xs_pdf->SetPointError(iPoint+1,0.,0.,sqrt(fabs(itVec->second.at(2))*fabs(itVec->second.at(2))+fabs(itVec->second.at(0)*
                                                  PDF_Error.at(iPoint).second/100.)*fabs(itVec->second.at(0)*PDF_Error.at(iPoint).second/100.))*
                                                  NumberOfBranchingRatio,NumberOfBranchingRatio*sqrt(fabs(itVec->second.at(1))*fabs(itVec->second.at(1))+
                                                  fabs(itVec->second.at(0)*PDF_Error.at(iPoint).first/100.)*fabs(itVec->second.at(0)*
                                                  PDF_Error.at(iPoint).first/100.)));
	iPoint++;
    }
  */

  // Powheg Analysis

  std::ifstream FilePowheg(PowhegInputFile.c_str());
  if(FilePowheg.fail()){ std::cerr<<" Fail to open powheg file --> Exit "<<std::endl; return 1;}

  std::map <float, Sample_powheg > map_powheg;

  while (!FilePowheg.eof()){
  
    Sample_powheg *PowhegInfo = new Sample_powheg ();

     getline(FilePowheg,buffer);
     if(buffer=="" || !buffer.find('#')) continue;
     std::stringstream line( buffer );
   
     line >> PowhegInfo->HiggsMass;      
     line >> PowhegInfo->efficiency; 
     line >> PowhegInfo->err_efficiency;  
     float temp;
     line >> temp ; PowhegInfo->xs          = temp *1000; 
     line >> temp ; PowhegInfo->err_xs_up   = temp *1000;  
     line >> temp ; PowhegInfo->err_xs_down = temp *1000;   

     PowhegInfo->xs_eff          = PowhegInfo->efficiency * PowhegInfo->xs;
     PowhegInfo->err_xs_eff_up   = sqrt(PowhegInfo->xs*PowhegInfo->xs*PowhegInfo->err_efficiency*PowhegInfo->err_efficiency+
                                        PowhegInfo->efficiency*PowhegInfo->efficiency*PowhegInfo->err_xs_up*PowhegInfo->err_xs_up );
     PowhegInfo->err_xs_eff_down = sqrt(PowhegInfo->xs*PowhegInfo->xs*PowhegInfo->err_efficiency*PowhegInfo->err_efficiency+
                                        PowhegInfo->efficiency*PowhegInfo->efficiency*PowhegInfo->err_xs_down*PowhegInfo->err_xs_down );

    
     map_powheg[PowhegInfo->HiggsMass]= *(PowhegInfo);
     delete PowhegInfo ;
     
  }

 
  TGraphAsymmErrors* mcfm_powheg_xs_ratio = new TGraphAsymmErrors();
  mcfm_powheg_xs_ratio->SetMarkerStyle(20);
  mcfm_powheg_xs_ratio->SetMarkerSize(2);
  mcfm_powheg_xs_ratio->SetMarkerColor(kBlue);
  mcfm_powheg_xs_ratio->SetLineColor(kBlue);
  mcfm_powheg_xs_ratio->SetLineWidth(2);

  TGraphAsymmErrors* powheg_mcfm_xs_ratio = new TGraphAsymmErrors();
  powheg_mcfm_xs_ratio->SetMarkerStyle(20);
  powheg_mcfm_xs_ratio->SetMarkerSize(2);
  powheg_mcfm_xs_ratio->SetMarkerColor(kRed);
  powheg_mcfm_xs_ratio->SetLineColor(kRed);
  powheg_mcfm_xs_ratio->SetLineWidth(2);

  iPoint   =0;

  OutputDump<<" ---------------------------------------------------------------------------------------------------------------------- "<<std::endl;
  OutputDump<<" ---------------------------------- Scale Factor Value and Error ------------------------------------------------------ "<<std::endl;
  OutputDump<<" ---------------------------------------------------------------------------------------------------------------------- "<<std::endl;
  OutputDump<<" -- Mass -- XS_mcfm -- Error up -- Error dw -- Xs_powheg -- Error up -- Error dw -- Scale Factor -- err_up -- err-dw -- "<<std::endl;
  OutputDump<<" ---------------------------------------------------------------------------------------------------------------------- "<<std::endl;
  OutputDump<<std::endl;

  for(std::map<float,Sample_powheg>::iterator itVec =map_powheg.begin();  itVec !=map_powheg.end(); ++itVec){

   double x_val      = (*itVec).first;
   double powh_val   = (*itVec).second.xs_eff;
   double powh_err_d = (*itVec).second.err_xs_eff_down;
   double powh_err_u = (*itVec).second.err_xs_eff_up;
   
   std::map<float,std::vector<float> >::iterator itMap = XSec_Scale.find(x_val);
   double mcfm_val=0.; double mcfm_err_d=0.; double mcfm_err_u=0.;
   double ratio=0.; double err_u=0.; double err_d=0.;

   if (itMap!=XSec_Scale.end()) { 
                                  mcfm_val   = itMap->second.at(0)*NumberOfBranchingRatio;
                                  mcfm_err_d = fabs(itMap->second.at(2))*NumberOfBranchingRatio;
                                  mcfm_err_u = fabs(itMap->second.at(1))*NumberOfBranchingRatio;

   }

   if (mcfm_val!=0 && powh_val!=0){
                     ratio = powh_val / mcfm_val;
                     err_u = sqrt((powh_err_u*powh_err_u)/(mcfm_val*mcfm_val)+(mcfm_err_d*mcfm_err_d)*(powh_val/mcfm_val/mcfm_val)*(powh_val/mcfm_val/mcfm_val));
                     err_d = sqrt((powh_err_d*powh_err_d)/(mcfm_val*mcfm_val)+(mcfm_err_u*mcfm_err_u)*(powh_val/mcfm_val/mcfm_val)*(powh_val/mcfm_val/mcfm_val));
                     powheg_mcfm_xs_ratio->SetPoint(iPoint,x_val,ratio);
                     powheg_mcfm_xs_ratio->SetPointError(iPoint,0.,0.,err_d,err_u);
    
                     ratio = mcfm_val / powh_val;
                     err_u = sqrt(mcfm_err_u * mcfm_err_u /(powh_val*powh_val));
                     err_d = sqrt(mcfm_err_d * mcfm_err_d /(powh_val*powh_val));

                     mcfm_powheg_xs_ratio->SetPoint(iPoint,x_val,ratio);
                     mcfm_powheg_xs_ratio->SetPointError(iPoint,0.,0.,err_d,err_u);
                     iPoint++;
       
		     OutputDump<<"   "<<x_val<<"     "<<mcfm_val<<"    "<<mcfm_err_u<<"    "<<mcfm_err_d<<"    "<<powh_val<<"    "<<powh_err_u<<"    "
                               <<powh_err_d  <<"     "<<ratio   <<"    "<<err_u     <<"    "<<err_d<<std::endl;
                     } 
				  
    
  }

  
  TGraphAsymmErrors* powheg_xs = new TGraphAsymmErrors();
  powheg_xs->SetMarkerStyle(20);
  powheg_xs->SetMarkerSize(1);
  powheg_xs->SetMarkerColor(kRed);
  powheg_xs->SetLineColor(kRed);
  powheg_xs->SetLineWidth(2);

  
  iPoint =0;
  for(std::map<float,Sample_powheg>::iterator itVec =map_powheg.begin();  itVec !=map_powheg.end(); ++itVec){

   powheg_xs->SetPoint(iPoint,(*itVec).first,(*itVec).second.xs_eff);
   powheg_xs->SetPointError(iPoint,0.,0.,(*itVec).second.err_xs_eff_down,(*itVec).second.err_xs_eff_up);
   iPoint++;
  }

  TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);

  leg->AddEntry(mcfm_xs,  "MCFM   @ VHWW 8TeV ", "LP");
  leg->AddEntry(powheg_xs,"Powheg @ VHWW 8TeV ", "LP");
  
  TCanvas *c = new TCanvas("powheg_vs_mcfm","powheg_vs_mcfm");
  c->cd();
  c->SetGridx();
  c->SetGridy();
  c->SetFillColor(0);

  powheg_xs->GetHistogram()->GetYaxis()-> SetRangeUser(HistoScaleMin,HistoScaleMax);
  powheg_xs->GetHistogram()->GetXaxis()-> SetRangeUser(110.,200.);
  powheg_xs->GetHistogram()->GetYaxis()-> SetTitle("#sigma (fb)");
  powheg_xs->GetHistogram()->GetXaxis()-> SetTitle("m_{H}");
  powheg_xs->SetFillColor(0);
  powheg_xs->SetFillStyle(3001);
  powheg_xs->Draw("ap30");

  mcfm_xs  ->SetFillColor(0);
  mcfm_xs  ->SetFillStyle(3002);
  mcfm_xs  ->Draw("p30same");

  leg->Draw("same");

  c->Write();

  powheg_xs ->Write("powheg_xs");
  mcfm_xs   ->Write("mcfm_xs");

  
  
  TLegend * leg2 = new TLegend(0.6,0.7,0.89, 0.89);
  leg2->SetFillColor(kWhite);
  leg2->SetFillStyle(0);

  leg2->AddEntry(powheg_mcfm_xs_ratio,"powheg / mcfm", "LP");
  
  TCanvas *cp = new TCanvas("powheg_vs_mcfm_ratio","powheg_vs_mcfm_ratio");
  cp->cd();
  cp->SetGridx();
  cp->SetGridy();
  cp->SetFillColor(0);
  powheg_mcfm_xs_ratio->GetHistogram()->GetXaxis()-> SetRangeUser(110.,200.);
  powheg_mcfm_xs_ratio->GetHistogram()->GetYaxis()-> SetTitle("#frac{#sigma_{POWHEG}}{#sigma_{mcfm}} @ 8TeV");
  powheg_mcfm_xs_ratio->GetHistogram()->GetXaxis()-> SetTitle("m_{H}");
  powheg_mcfm_xs_ratio->SetFillColor(0);

  powheg_mcfm_xs_ratio->Draw("ap30");
  leg2->Draw("same");

  cp->Write();

   
  
  TLegend * leg3 = new TLegend(0.6,0.7,0.89, 0.89);
  leg3->SetFillColor(kWhite);
  leg3->SetFillStyle(0);
 
  leg3->AddEntry(mcfm_powheg_xs_ratio,"mcfm / powheg", "LP");
  
  TCanvas *cp3 = new TCanvas("mcfm_powheg_xs_ratio","mcfm_powheg_xs_ratio");
  cp3->cd();
  cp3->SetGridx();
  cp3->SetGridy();
  cp3->SetFillColor(0);
  mcfm_powheg_xs_ratio->GetHistogram()->GetXaxis()-> SetRangeUser(110.,200.);
  mcfm_powheg_xs_ratio->GetHistogram()->GetYaxis()-> SetTitle("#frac{#sigma_{mcfm}}{#sigma_{POWHEG}} @ 8TeV");
  mcfm_powheg_xs_ratio->GetHistogram()->GetXaxis()-> SetTitle("m_{H}");
  mcfm_powheg_xs_ratio->SetFillColor(0);

  mcfm_powheg_xs_ratio->Draw("ap30");
  leg3->Draw("same");

  cp3->Write();

  
    
  TCanvas *c4 = new TCanvas("MCFM Error","MCFM Error");
  c4->cd();
  c4->SetGridx();
  c4->SetGridy();
  c4->SetFillColor(0);
  TLegend * leg4 = new TLegend(0.6,0.7,0.89, 0.89);
  leg4->SetFillColor(0);
  leg4->AddEntry(mcfm_xs_stat,"MCFM +/- Statistical Error", "LP");
  leg4->AddEntry(mcfm_xs,"MCFM +/- Scale Uncertainty", "LP");
  // leg4->AddEntry(mcfm_xs_pdf,"MCFM +/- Scale + PDF Error", "LP");

  mcfm_xs_stat->GetHistogram()->GetYaxis()-> SetRangeUser(0.,10.);
  mcfm_xs_stat->GetHistogram()->GetXaxis()-> SetRangeUser(110.,200.);
  mcfm_xs_stat->GetHistogram()->GetYaxis()-> SetTitle("#sigma (fb)");
  mcfm_xs_stat->GetHistogram()->GetXaxis()-> SetTitle("m_{H}");
  mcfm_xs_stat->SetFillColor(0);
  mcfm_xs_stat->SetFillColor(0);
  mcfm_xs_stat->Draw("ap30");
  mcfm_xs->Draw("p30same");
  //  mcfm_xs_pdf->Draw("p3same");

  leg4->Draw("same");
  c4->Write();
  outputFile->Close();

  return 0;
 
}  

//**************************************************************************************

float GetMaximum_RS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass){

   float max=-1.;
   for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){
    
      if((*itVec).Higgs_Mass != Higgs_Mass) continue;
      if((*itVec).Rin_Scale >= max) max=(*itVec).Rin_Scale;
   }
 
   if(max ==-1.) return -1.;
   else return max;
}   


float GetMinimum_RS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass){

  float min=100000000.;
  for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){    
     if((*itVec).Higgs_Mass != Higgs_Mass) continue;
     if((*itVec).Rin_Scale <= min) min=(*itVec).Rin_Scale;
  }
 
  if(min ==100000000.) return -1.;
  else return min;
}           

float GetMaximum_FS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass){

  float max=-1.;
  for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){
    
    if((*itVec).Higgs_Mass != Higgs_Mass) continue;
    if((*itVec).Fact_Scale >= max) max=(*itVec).Fact_Scale;
  }
 
  if(max ==-100.) return -1.;
  else return max;
}  

float GetMinimum_FS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass){

  float min=100000000.;
  for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec){
    
    if((*itVec).Higgs_Mass != Higgs_Mass) continue;
    if((*itVec).Fact_Scale <= min) min=(*itVec).Fact_Scale;
  }
 
  if(min ==100000000.) return -1.;
  else return min;
}           

/// Generate File with job information for the analysis
//ls Virtual/ | grep HWW2jt_virt_mstw8nl_ | grep .dat | tr "w" " " | tr "d" " " | tr "_" " " | awk '{print "echo @w"$7"_"$8"_"$9"_"$10"_"$11"_"$12" "$9" "$10" "$11" "$12"@ >> Virtual.txt"}' | tr "@" "\""  | /bin/sh

//ls Real/ | grep HWW2jt_real_mstw8nl_ | grep .dat | tr "w" " " | tr "d" " " | tr "_" " " | awk '{print "echo @w"$7"_"$8"_"$9"_"$10"_"$11"_"$12" "$9" "$10" "$11" "$12"@ >> Real.txt"}' | tr "@" "\""  | /bin/sh

// Complie with g++ -Wall `root-config --libs --glibs --cflags` -o Analysis Analysis.cxx
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include "TString.h"
#include "TH1F.h"
#include <algorithm>
#include "TFile.h"
#include <cmath>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TCanvas.h"


using namespace std;

struct Sample {
               float Higgs_Mass;
               float Fact_Scale;
               float Rin_Scale;
               int   Run;

              };

struct Sample_xs {
               
                  float Higgs_Mass;
                  float Fact_Scale;
                  float Rin_Scale;
                  TH1F* XSec_Virt;
                  TH1F* XSec_Real;
                  float err_prop_Virt;
                  float err_prop_Real;
                  float xs_value_Run;
                  float xs_error_Run;
                  float xs_error_RMS_Run;
           
                  bool operator==(const Sample_xs & a)
                  {
                    if(Higgs_Mass == a.Higgs_Mass && Fact_Scale==a.Fact_Scale && Rin_Scale==a.Rin_Scale)
                    return true;
                    else return false;
                  }

                  bool operator!=(const Sample_xs & a)
                  {
                   if(Higgs_Mass != a.Higgs_Mass || Fact_Scale!=a.Fact_Scale || Rin_Scale!=a.Rin_Scale)
                   return true;
                   else return false;
                  }
               
                  void operator=(const Sample_xs & a)
                  {
                   Higgs_Mass = a.Higgs_Mass;
                   Fact_Scale = a.Fact_Scale;
                   Rin_Scale  = a.Rin_Scale;
                   err_prop_Virt = a.err_prop_Virt;
                   err_prop_Real = a.err_prop_Real;
                   XSec_Virt = (TH1F*) a.XSec_Virt->Clone(XSec_Virt->GetName());
                   XSec_Real = (TH1F*) a.XSec_Real->Clone(XSec_Real->GetName());
                   xs_value_Run = a.xs_value_Run;
                   xs_error_Run = a.xs_error_Run;
                   xs_error_RMS_Run = a.xs_error_RMS_Run;
               
                 }
               
             
             };


struct Sample_Scale {
                      float Higgs_Mass;
                      float Fact_Scale_Max;
                      float Fact_Scale_Min;
                      float Rin_Scale_Max;
                      float Rin_Scale_Min;

                      bool operator==(const Sample_Scale & a)
                      {
                        if(Higgs_Mass == a.Higgs_Mass)
                        return true;
                        else return false;
                      }
                      
                      bool operator!=(const Sample_Scale & a)
                      {
                          if(Higgs_Mass != a.Higgs_Mass)
                          return true;
                          else return false;
                      }
               
                      void operator=(const Sample_Scale & a)
                      {
                        Higgs_Mass = a.Higgs_Mass;
                        Fact_Scale_Max = a.Fact_Scale_Max;
                        Fact_Scale_Min  = a.Fact_Scale_Min;
                        Rin_Scale_Max = a.Rin_Scale_Max;
                        Rin_Scale_Min = a.Rin_Scale_Min;
                      }

                   };

struct Sample_powheg {
                       
                       float efficiency;
                       float err_efficiency;
                       float xs;
                       float err_xs_up;
                       float err_xs_down;
                       float xs_eff;
                       float err_xs_eff_up;
                       float err_xs_eff_down;
                     };

float GetMaximum_RS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass)
{
   float max=-1.;
   for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
   {    
      if((*itVec).Higgs_Mass != Higgs_Mass) continue;
      if((*itVec).Rin_Scale >= max) max=(*itVec).Rin_Scale;
   }
 
   if(max ==-1.) return -1.;
   else return max;
}   

float GetMinimum_RS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass)
{
  float min=100000000.;
  for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
  {    
     if((*itVec).Higgs_Mass != Higgs_Mass) continue;
     if((*itVec).Rin_Scale <= min) min=(*itVec).Rin_Scale;
  }
 
  if(min ==100000000.) return -1.;
  else return min;
}           
               

float GetMaximum_FS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass)
{
  float max=-100.;
  for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
  {    
    if((*itVec).Higgs_Mass != Higgs_Mass) continue;
    if((*itVec).Fact_Scale >= max) max=(*itVec).Fact_Scale;
  }
 
  if(max ==-100.) return -1.;
  else return max;
}  

float GetMinimum_FS(const std::vector<Sample_xs> & Vector_xs_sample, float Higgs_Mass)
{
  float min=100000000.;
  for( std::vector<Sample_xs>::const_iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
  {    
    if((*itVec).Higgs_Mass != Higgs_Mass) continue;
    if((*itVec).Fact_Scale <= min) min=(*itVec).Fact_Scale;
  }
 
  if(min ==100000000.) return -1.;
  else return min;
}           
            

int main (int argc, char** argv)

{
  std::ifstream FileVirtual(argv[1]);
  std::string buffer;
  std::map <std::string, Sample > map_Virt;
  
  /// Bulding map for virtual jobs

  while (!FileVirtual.eof())
  {  
     std::string name_temp;
     float temp;
     Sample infos;

     getline(FileVirtual,buffer);
     if (buffer != ""){ ///---> save from empty line at the end!
     if (buffer.at(0) != '#'){
     std::stringstream line( buffer );
     line >> name_temp;
     line >> temp;
     infos.Higgs_Mass = temp;
     line >> temp;
     infos.Fact_Scale = temp;
     line >> temp;
     infos.Rin_Scale = temp;
     line >> temp;
     infos.Run = temp;
     map_Virt[name_temp]= infos;
     }
    }
  }

  std::ifstream FileReal(argv[2]);
  std::map <std::string, Sample > map_Real;

  /// Buolding map for real job
  
  while (!FileReal.eof())
  {  
     std::string name_temp;
     float temp;
     Sample infos;

     getline(FileReal,buffer);
     if (buffer != ""){ ///---> save from empty line at the end!
     if (buffer.at(0) != '#'){
     std::stringstream line( buffer );
     line >> name_temp;
     line >> temp;
     infos.Higgs_Mass = temp;
     line >> temp;
     infos.Fact_Scale = temp;
     line >> temp;
     infos.Rin_Scale = temp;
     line >> temp;
     infos.Run = temp;
     map_Real[name_temp]= infos;
     }
    }
  }

  std::vector<Sample_xs> Vector_xs_sample;


 /// Virtual cross section distribution for different run

  for(std::map <std::string, Sample >::const_iterator itMap = map_Virt.begin();
      itMap != map_Virt.end() ; itMap++)
  {
    
    TString nameFileIn;
    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale<100. )
     nameFileIn = Form("Virtual/HWW2jt_virt_mstw8nl_%d__%d__%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale>=100. )
     nameFileIn = Form("Virtual/HWW2jt_virt_mstw8nl_%d_%d_%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
    
    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale>=100. )
     nameFileIn = Form("Virtual/HWW2jt_virt_mstw8nl_%d_%d__%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale<100. )
     nameFileIn = Form("Virtual/HWW2jt_virt_mstw8nl_%d__%d_%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
 
     std::ifstream FileIn (nameFileIn.Data());
  
     getline(FileIn,buffer);
     std::stringstream line_head( buffer );
     float xsec, errxsec;
     std::string temp;

     line_head>>temp>>temp;
     if(temp == "Intermediate") continue;
 
     getline(FileIn,buffer);
     getline(FileIn,buffer);
 
     std::stringstream line (buffer);
     line>>temp>>temp>>temp>>temp;
     line>>xsec;
     line>>temp;
     line>>errxsec;
     
     Sample_xs var;
     var.Higgs_Mass = itMap->second.Higgs_Mass;
     var.Fact_Scale = itMap->second.Fact_Scale;
     var.Rin_Scale = itMap->second.Rin_Scale;
     
     std::vector<Sample_xs>::iterator pos = find(Vector_xs_sample.begin(),Vector_xs_sample.end(),var);

     if(pos==Vector_xs_sample.end() || (*pos).XSec_Virt==0) {
                 var.XSec_Virt = new TH1F (itMap->first.c_str(),itMap->first.c_str(),1000.,0.,3.);
                 std::string name_tmp = itMap->first+"_Real";
                 var.XSec_Real = new TH1F (name_tmp.c_str(),name_tmp.c_str(),1000.,-3.,0.);
           
                 var.XSec_Virt->Fill(xsec);
                 var.err_prop_Virt = errxsec*errxsec;
                 Vector_xs_sample.push_back(var);
                 }
      else{
            (*pos).XSec_Virt->Fill(xsec);
            (*pos).err_prop_Virt=(*pos).err_prop_Virt+errxsec*errxsec;}
  }


/// Real cross section distribution for different run

  for(std::map <std::string, Sample >::const_iterator itMap = map_Real.begin();
      itMap != map_Real.end() ; itMap++)
  {
    
    TString nameFileIn;
    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale<100. )
     nameFileIn = Form("Real/HWW2jt_real_mstw8nl_%d__%d__%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale>=100. )
     nameFileIn = Form("Real/HWW2jt_real_mstw8nl_%d_%d_%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
    
    if(itMap->second.Fact_Scale<100. && itMap->second.Rin_Scale>=100. )
     nameFileIn = Form("Real/HWW2jt_real_mstw8nl_%d_%d__%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
    if(itMap->second.Fact_Scale>=100. && itMap->second.Rin_Scale<100. )
     nameFileIn = Form("Real/HWW2jt_real_mstw8nl_%d__%d_%sdat",int(itMap->second.Rin_Scale),
                               int(itMap->second.Fact_Scale),itMap->first.c_str());
 
     std::ifstream FileIn (nameFileIn.Data());
  
     getline(FileIn,buffer);
     std::stringstream line_head( buffer );
     float xsec, errxsec;
     std::string temp;

     line_head>>temp>>temp;
     if(temp == "Intermediate") continue;
 
     getline(FileIn,buffer);
     getline(FileIn,buffer);
 
     std::stringstream line (buffer);
     line>>temp>>temp>>temp>>temp;
     line>>xsec;
     line>>temp;
     line>>errxsec;
     Sample_xs var;
     var.Higgs_Mass = itMap->second.Higgs_Mass;
     var.Fact_Scale = itMap->second.Fact_Scale;
     var.Rin_Scale = itMap->second.Rin_Scale;
     
     std::vector<Sample_xs>::iterator pos = find(Vector_xs_sample.begin(),Vector_xs_sample.end(),var);
     (*pos).XSec_Real->Fill(xsec);
     (*pos).err_prop_Real=(*pos).err_prop_Real+errxsec*errxsec;
      
  }

  cout<<"------------------------------------------------------------------------------------------"<<endl;
  cout<<"-------------------Value obtained from mean for each scale and Higgs Mass-----------------"<<endl;
  cout<<"------------------------------------------------------------------------------------------"<<endl;
  cout<<"Mass -- Fact Scale -- Rin Scale -- XSec Virt -- Err XSec Virt -- Err XSec Propagation Virt --"
  " XSec Real -- Err XSec Real ---Err XSec Propagation Real"<<endl;
  cout<<"------------------------------------------------------------------------------------------"<<endl;
  

 
  for( std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
  {
   
   cout<<(*itVec).Higgs_Mass<<"        "<<(*itVec).Fact_Scale<<"         "<<(*itVec).Rin_Scale<<"          "<<(*itVec).XSec_Virt->GetMean()<<"       "<<(*itVec).XSec_Virt->GetMeanError()<<"           "<<sqrt((*itVec).err_prop_Virt/(*itVec).XSec_Virt->GetEntries())<<"       "<<
    (*itVec).XSec_Real->GetMean()<<"        "<<(*itVec).XSec_Real->GetMeanError()<< "       "<<sqrt((*itVec).err_prop_Real/(*itVec).XSec_Real->GetEntries())<<endl;
  }


  cout<<"------------------------------------------------------------"<<endl;
  cout<<"---Value obtained from mean for each scale and Higgs Mass---"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Mass -- Fact Scale -- Rin Scale -- XSec -- Err XSec---Err XSec Propagation---"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  

  std::vector<float> Higgs_Mass_Vect;

  for( std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
  {
   if(Higgs_Mass_Vect.size()==0) Higgs_Mass_Vect.push_back((*itVec).Higgs_Mass);
   else
   {
    std::vector<float>::iterator pos = find (Higgs_Mass_Vect.begin(),Higgs_Mass_Vect.end(),(*itVec).Higgs_Mass);
    if(pos == Higgs_Mass_Vect.end()) Higgs_Mass_Vect.push_back((*itVec).Higgs_Mass); 
                                   
   }
  
  
   (*itVec).xs_value_Run = (*itVec).XSec_Virt->GetMean()+(*itVec).XSec_Real->GetMean();
   (*itVec).xs_error_Run  = sqrt((*itVec).err_prop_Virt/(*itVec).XSec_Virt->GetEntries()+(*itVec).err_prop_Real/(*itVec).XSec_Real->GetEntries());
   (*itVec).xs_error_RMS_Run   =  sqrt((*itVec).XSec_Virt->GetMeanError()*(*itVec).XSec_Virt->GetMeanError()+(*itVec).XSec_Real->GetMeanError()*(*itVec).XSec_Real->GetMeanError());
   cout<<(*itVec).Higgs_Mass<<"        "<<(*itVec).Fact_Scale<<"      "<<(*itVec).Rin_Scale<<"     "<<(*itVec).xs_value_Run<<"    "<<(*itVec).xs_error_RMS_Run<<"      "<<(*itVec).xs_error_Run<<endl;
  }

  TFile *output = new TFile(argv[3],"RECREATE");
  output->cd();
  
  
  for(unsigned int i=0; i<Vector_xs_sample.size(); i++)
  {    Vector_xs_sample.at(i).XSec_Virt->Write();
       Vector_xs_sample.at(i).XSec_Real->Write();
  }
   
 

  /// Scale study

 std::vector<Sample_Scale> map_Scale;

 for( std::vector<float>::iterator itVec =Higgs_Mass_Vect.begin();  itVec !=Higgs_Mass_Vect.end(); ++itVec)
  {
      float max_RS = GetMaximum_RS(Vector_xs_sample,*itVec);
      float min_RS = GetMinimum_RS(Vector_xs_sample,*itVec);
      float max_FS = GetMaximum_FS(Vector_xs_sample,*itVec);
      float min_FS = GetMinimum_FS(Vector_xs_sample,*itVec);
      Sample_Scale var;
      var.Higgs_Mass = (*itVec);
      var.Fact_Scale_Max = max_FS; 
      var.Fact_Scale_Min = min_FS;
      var.Rin_Scale_Max = max_RS; 
      var.Rin_Scale_Min = min_RS;
      map_Scale.push_back(var);
  }
 
  
  std::map<float,std::vector<float> > XSec_Scale,XSec_Scale_Virt,XSec_Scale_Real;

  for(std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
  {   
      Sample_Scale var ;
      var.Higgs_Mass = (*itVec).Higgs_Mass;
      std::vector<Sample_Scale>::iterator pos = find (map_Scale.begin(),map_Scale.end(),var);
      if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
      if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
      
      std::map<float,std::vector<float> >::iterator itMap = XSec_Scale.find((*itVec).Higgs_Mass);
      if(itMap == XSec_Scale.end()) XSec_Scale[(*itVec).Higgs_Mass].assign(3,0);
     
      if((*itVec).Fact_Scale ==0.5*(*itVec).Higgs_Mass && (*itVec).Rin_Scale ==0.5*(*itVec).Higgs_Mass)
      XSec_Scale[(*itVec).Higgs_Mass].at(0)=(*itVec).xs_value_Run;
    
      std::map<float,std::vector<float> >::iterator itMap_Virt = XSec_Scale_Virt.find((*itVec).Higgs_Mass);
      if(itMap_Virt == XSec_Scale_Virt.end()) XSec_Scale_Virt[(*itVec).Higgs_Mass].assign(3,0);
     
      if((*itVec).Fact_Scale ==0.5*(*itVec).Higgs_Mass && (*itVec).Rin_Scale ==0.5*(*itVec).Higgs_Mass)
      XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0)=(*itVec).XSec_Virt->GetMean();

      std::map<float,std::vector<float> >::iterator itMap_Real = XSec_Scale_Real.find((*itVec).Higgs_Mass);
      if(itMap_Real == XSec_Scale_Real.end()) XSec_Scale_Real[(*itVec).Higgs_Mass].assign(3,0);
     
      if((*itVec).Fact_Scale ==0.5*(*itVec).Higgs_Mass && (*itVec).Rin_Scale ==0.5*(*itVec).Higgs_Mass)
      XSec_Scale_Real[(*itVec).Higgs_Mass].at(0)=(*itVec).XSec_Real->GetMean();
  }

  for(std::vector<Sample_xs>::iterator itVec =Vector_xs_sample.begin();  itVec !=Vector_xs_sample.end(); ++itVec)
  {   
      Sample_Scale var ;
      var.Higgs_Mass = (*itVec).Higgs_Mass;
      std::vector<Sample_Scale>::iterator pos = find (map_Scale.begin(),map_Scale.end(),var);
      if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
      if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      if((*pos).Fact_Scale_Max == (*itVec).Fact_Scale && (*pos).Rin_Scale_Min == (*itVec).Rin_Scale) continue;
      if((*pos).Fact_Scale_Min == (*itVec).Fact_Scale && (*pos).Rin_Scale_Max == (*itVec).Rin_Scale) continue;
     
     
      if((*itVec).Fact_Scale ==0.5*(*itVec).Higgs_Mass && (*itVec).Rin_Scale ==0.5*(*itVec).Higgs_Mass) continue;
      
      if(XSec_Scale[(*itVec).Higgs_Mass].at(1)<((*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0))) 
      XSec_Scale[(*itVec).Higgs_Mass].at(1)=(*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0);
       
      if(XSec_Scale[(*itVec).Higgs_Mass].at(2)>((*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0))) XSec_Scale[(*itVec).Higgs_Mass].at(2)=(*itVec).xs_value_Run-XSec_Scale[(*itVec).Higgs_Mass].at(0);
            
      if(XSec_Scale_Virt[(*itVec).Higgs_Mass].at(1)<(*itVec).XSec_Virt->GetMean()-XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0)) XSec_Scale_Virt[(*itVec).Higgs_Mass].at(1)=(*itVec).XSec_Virt->GetMean()-XSec_Scale[(*itVec).Higgs_Mass].at(0);
      
      if(XSec_Scale_Virt[(*itVec).Higgs_Mass].at(2)>(*itVec).XSec_Virt->GetMean()-XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0)) XSec_Scale_Virt[(*itVec).Higgs_Mass].at(2)=(*itVec).XSec_Virt->GetMean()-XSec_Scale_Virt[(*itVec).Higgs_Mass].at(0);
      
      if(XSec_Scale_Real[(*itVec).Higgs_Mass].at(1)<(*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0)) XSec_Scale_Real[(*itVec).Higgs_Mass].at(1)=(*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0);
         
      if(XSec_Scale_Real[(*itVec).Higgs_Mass].at(2)>(*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0)) XSec_Scale_Real[(*itVec).Higgs_Mass].at(2)=(*itVec).XSec_Real->GetMean()-XSec_Scale_Real[(*itVec).Higgs_Mass].at(0);
      
  }
         
 


  TGraphAsymmErrors* mcfm_xs = new TGraphAsymmErrors();
  mcfm_xs->SetMarkerStyle(20);
  mcfm_xs->SetMarkerSize(1);
  mcfm_xs->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors* mcfm_virtual = new TGraphAsymmErrors();
  mcfm_virtual->SetMarkerStyle(20);
  mcfm_virtual->SetMarkerSize(1);
  mcfm_virtual->SetMarkerColor(kBlue+2);

  TGraphAsymmErrors* mcfm_real = new TGraphAsymmErrors();
  mcfm_real->SetMarkerStyle(20);
  mcfm_real->SetMarkerSize(1);
  mcfm_real->SetMarkerColor(kBlue+2);


  int iPoint=0;

  cout<<"------------------------------------------------------------"<<endl;
  cout<<"-------Scale dependence for each Higgs Mass evaluation Virtual ------"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Mass ------- XSec ------- Err XSec Up ----- Err XSec Down"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  for(std::map<float,std::vector<float> >::iterator itVec =XSec_Scale_Virt.begin();  itVec !=XSec_Scale_Virt.end(); ++itVec)
  {  cout<<itVec->first<<"         "<<itVec->second.at(0)<<"       "<<itVec->second.at(1)<<"    "<<itVec->second.at(2)<<endl;
     mcfm_virtual->SetPoint(iPoint+1,itVec->first,itVec->second.at(0));
     mcfm_virtual->SetPointError(iPoint+1,0.,0.,fabs(itVec->second.at(2)),itVec->second.at(1));
     iPoint++;
  }

  cout<<"------------------------------------------------------------"<<endl;
  cout<<"-------Scale dependence for each Higgs Mass evaluation Real------"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Mass ------- XSec ------- Err XSec Up ----- Err XSec Down"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  iPoint =0;

  for(std::map<float,std::vector<float> >::iterator itVec =XSec_Scale_Real.begin();  itVec !=XSec_Scale_Real.end(); ++itVec)
  {  cout<<itVec->first<<"         "<<itVec->second.at(0)<<"       "<<itVec->second.at(1)<<"    "<<itVec->second.at(2)<<endl;
     mcfm_real->SetPoint(iPoint+1,itVec->first,itVec->second.at(0));
     mcfm_real->SetPointError(iPoint+1,0.,0.,fabs(itVec->second.at(2)),itVec->second.at(1));
     iPoint++;
  }


  cout<<"------------------------------------------------------------"<<endl;
  cout<<"-------Scale dependence for each Higgs Mass evaluation------"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Mass ------- XSec ------- Err XSec Up ----- Err XSec Down"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  iPoint =0;

  for(std::map<float,std::vector<float> >::iterator itVec =XSec_Scale.begin();  itVec !=XSec_Scale.end(); ++itVec)
  {  cout<<itVec->first<<"         "<<itVec->second.at(0)<<"       "<<itVec->second.at(1)<<"    "<<itVec->second.at(2)<<endl;
     mcfm_xs->SetPoint(iPoint+1,itVec->first,itVec->second.at(0));
     mcfm_xs->SetPointError(iPoint+1,0.,0.,fabs(itVec->second.at(2)),itVec->second.at(1));
     iPoint++;
  }
     
 
  
  /// Comparison with powheg + yellow report result

  std::ifstream FilePowheg(argv[4]);
  std::map <float, Sample_powheg > map_powheg;
  
  /// Bulding map for virtual jobs

  while (!FilePowheg.eof())
  {  
     float higgs;
     float temp;
     Sample_powheg infos;

     getline(FilePowheg,buffer);
     if (buffer != ""){ ///---> save from empty line at the end!
     if (buffer.at(0) != '#'){
     std::stringstream line( buffer );
     line >> higgs;
     line >> temp;
     infos.efficiency = temp;
     line >> temp;
     infos.err_efficiency = temp;
     line >> temp;
     infos.xs = temp*1000.;
     line >> temp;
     infos.err_xs_up = temp*1000.;
     line >> temp;
     infos.err_xs_down = temp*1000.;
     infos.xs_eff = infos.efficiency * infos.xs;
     infos.err_xs_eff_up = sqrt(infos.xs*infos.xs*infos.err_efficiency*infos.err_efficiency+infos.efficiency*infos.efficiency*
                                infos.err_xs_up*infos.err_xs_up );
     infos.err_xs_eff_down = sqrt(infos.xs*infos.xs*infos.err_efficiency*infos.err_efficiency+infos.efficiency*infos.efficiency*
                                infos.err_xs_down*infos.err_xs_down );


     map_powheg[higgs]= infos;
     }
    }
  }

  TGraphAsymmErrors* powheg_xs = new TGraphAsymmErrors();
  powheg_xs->SetMarkerStyle(20);
  powheg_xs->SetMarkerSize(1);
  powheg_xs->SetMarkerColor(kRed);

  iPoint =0;
  for(std::map<float,Sample_powheg>::iterator itVec =map_powheg.begin();  itVec !=map_powheg.end(); ++itVec)
  {
   powheg_xs->SetPoint(iPoint+1,(*itVec).first,(*itVec).second.xs_eff);
   powheg_xs->SetPointError(iPoint+1,0.,0.,(*itVec).second.err_xs_eff_down,(*itVec).second.err_xs_eff_up);
   iPoint++;
  }
  

  TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
  leg->SetFillColor(0);
  leg->AddEntry(mcfm_xs,"mcfm", "LP");
  leg->AddEntry(powheg_xs,"powheg_xs", "LP");
  
  output->cd();
  TCanvas *c = new TCanvas("powheg_vs_mcfm","powheg_vs_mcfm");
  c->cd();
  c->SetGridx();
  c->SetGridy();
  c->SetFillColor(0);
  powheg_xs->GetHistogram()->SetTitle(" Powheg vs mcfm ");
  powheg_xs->GetHistogram()->GetYaxis()-> SetRangeUser(0.,3.);
  powheg_xs->GetHistogram()->GetXaxis()-> SetRangeUser(100.,300.);
  powheg_xs->GetHistogram()->GetYaxis()-> SetTitle("#sigma");
  powheg_xs->GetHistogram()->GetXaxis()-> SetTitle("m_{H}");
  powheg_xs->SetFillColor(0);
  mcfm_xs->SetFillColor(0);
  powheg_xs->Draw("ap");
  mcfm_xs->Draw("psame");
  leg->Draw("same");
  c->Write();


  TCanvas *c2 = new TCanvas("mcfm virtual","mcfm virtual");
  c2->cd();
  c2->SetGridx();
  c2->SetGridy();
  c2->SetFillColor(0);
  mcfm_virtual->GetHistogram()->SetTitle(" mcfm  virtual ");
  mcfm_virtual->GetHistogram()->GetYaxis()-> SetRangeUser(0.,3.);
  mcfm_virtual->GetHistogram()->GetXaxis()-> SetRangeUser(100.,300.);
  mcfm_virtual->GetHistogram()->GetYaxis()-> SetTitle("#sigma");
  mcfm_virtual->GetHistogram()->GetXaxis()-> SetTitle("m_{H}");
  mcfm_virtual->SetFillColor(0);
  mcfm_virtual->SetFillColor(0);
  mcfm_virtual->Draw("ap");
  c2->Write();
  
  TCanvas *c3 = new TCanvas("mcfm real","mcfm real");
  c3->cd();
  c3->SetGridx();
  c3->SetGridy();
  c3->SetFillColor(0);
  mcfm_real->GetHistogram()->SetTitle(" mcfm  real ");
  mcfm_real->GetHistogram()->GetYaxis()-> SetRangeUser(-3.,0.);
  mcfm_real->GetHistogram()->GetXaxis()-> SetRangeUser(100.,300.);
  mcfm_real->GetHistogram()->GetYaxis()-> SetTitle("#sigma");
  mcfm_real->GetHistogram()->GetXaxis()-> SetTitle("m_{H}");
  mcfm_real->SetFillColor(0);
  mcfm_real->SetFillColor(0);
  mcfm_real->Draw("ap");
  c3->Write();


 output->Close();

return 0;
}  
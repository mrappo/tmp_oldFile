#ifndef stdHisto_h
#define stdHisto_h

#include "treeReader.h"
#include "hFactory.h"
#include "ntpleUtils.h"

#include <iostream>
#include <vector>

#include "Math/GenVector/VectorUtil.h"



class stdHisto
{
 public:
  
  //!ctor
  stdHisto(const int& nStep,
           const std::string& outFileName,
           treeReader* reader = NULL);
  
  //!dtor
  ~stdHisto();
  
  
  //! methods
  
  // add histograms
  void Add1(const std::string& histoName,
            const int& nStep,
            const bool& doZeppPlots = false);
  
  void Add2(const std::string& histoName,
            const int& nStep,
            const bool& doZeppPlots = false);
  
  // fill histograms
  void Fill1(const std::string& histoName,
             const std::string& branchName,
             const int& nStep,
             std::vector<int>* selectionIt = NULL,
             const float* eta1_tag = NULL,
             const float* eta2_tag = NULL);
  
  void Fill1(const std::vector<ROOT::Math::XYZTVector>& vet,
             const std::string& histoName,
             const int& step,
             const float* eta1_tag = NULL,
             const float* eta2_tag = NULL);
  
  void Fill1(const ROOT::Math::XYZTVector& p,
             const std::string& histoName,
             const int& step,
             const float* eta1_tag = NULL,
             const float* eta2_tag = NULL);
  
  void Fill2(const std::string& histoName,
             const std::string& branchName,
             const int& nStep,
             const int& it1, const int& it2,
             const float* eta1_tag = NULL,
             const float* eta2_tag = NULL);
  
  void Fill2(const ROOT::Math::XYZTVector& v1,
             const ROOT::Math::XYZTVector& v2, 
             const std::string& histoName,
             const int& step,
             const float* eta1_tag = NULL,
             const float* eta2_tag = NULL);
  
  
  
 private:
  
  int m_nStep;
  treeReader* m_reader; 
  std::string m_outFileName;
  
  hFactory* m_hFactory;
  
};

#endif

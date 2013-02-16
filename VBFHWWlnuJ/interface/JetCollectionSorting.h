#ifndef JetCollectionSorting_h
#define JetCollectionSorting_h

#include <iostream>
#include <vector>
#include <functional>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lorentzVector ;


// two jets vbf with highest Deta first , then W jets with higher pT
std::vector<const lorentzVector *> Jet_Deta_pT       (const lorentzVector * input, int njets);

std::vector<const lorentzVector *> Jet_FNAL_Criteria (const lorentzVector * input, int njets);

// two jets vbf with highest Mjj first , then W jets with higher pT
std::vector<const lorentzVector *> Jet_Mjj_pT (const lorentzVector * input, int njets);

// lorentzVector pointers sorting algos
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#ifndef TLVP_EtaSort_h
#define TLVP_EtaSort_h

class TLVP_EtaSort : public std::binary_function<int,int,bool>{

 public: 
  
  // default constructor
  TLVP_EtaSort(){};

  // default de-constructor
  ~TLVP_EtaSort(){};

  // Operator ()
  bool operator() (const lorentzVector * x, const lorentzVector * y){
    return x->Eta () < y->Eta () ;
  };

} ;

#endif

#ifndef TLVP_AbsEtaSort_h
#define TLVP_AbsEtaSort_h

class TLVP_AbsEtaSort : public std::binary_function<int,int,bool> {

 public: 
  
  // default constructor
  TLVP_AbsEtaSort(){};

  // default de-constructor
  ~TLVP_AbsEtaSort(){};

  // Operator ()
  bool operator() (const lorentzVector * x, const lorentzVector * y){
    return fabs (x->Eta ()) < fabs (y->Eta ()) ;
  };

} ;

#endif

#ifndef TLVP_PtSort_h
#define TLVP_PtSort_h

class TLVP_PtSort : public std::binary_function<int,int,bool>{

 public :
  // default constructor
  TLVP_PtSort(){};

  // default de-constructor
  ~TLVP_PtSort(){};

  // Operator ()
  bool operator() (const lorentzVector * x, const lorentzVector * y){
      return x->Pt () < y->Pt () ;
  };
  
} ;

#endif

#endif

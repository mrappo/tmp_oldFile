#ifndef JetCollectionSorting_h
#define JetCollectionSorting_h

#include <iostream>
#include <vector>
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

struct TLVP_EtaSort{
  bool operator() (const lorentzVector * x, const lorentzVector * y){
    return x->Eta () < y->Eta () ;
  }
} ;


struct TLVP_AbsEtaSort{
  bool operator() (const lorentzVector * x, const lorentzVector * y){
    return fabs (x->Eta ()) < fabs (y->Eta ()) ;
  }
} ;

struct TLVP_PtSort{
  bool operator() (const lorentzVector * x, const lorentzVector * y){
    return x->Pt () < y->Pt () ;
  }
} ;



#endif

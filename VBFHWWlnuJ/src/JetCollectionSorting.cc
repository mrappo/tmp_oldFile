#include "JetCollectionSorting.h"


// default constructor
JetAK5::JetAK5(){}


// constructor
JetAK5::JetAK5(const int & position, const std::string & NameCollection, TLorentzVector & Vect){

  position_ = position ;
  NameCollection_ = NameCollection ;
  Momentum_ = Vect ;

}

// default de-constructor                                                                                                                                                                       
JetAK5::~JetAK5(){}


// Operator ()                                                                                                                                                                                  
bool TLVP_EtaSort::operator() (const JetAK5 & x, const JetAK5 & y){
  return x.Momentum_.Eta () < y.Momentum_.Eta () ;
};

// Operator ()                                                                                                                                                                                  
bool TLVP_AbsEtaSort::operator() (const JetAK5 & x, const JetAK5 & y){
  return fabs (x.Momentum_.Eta ()) < fabs (y.Momentum_.Eta ()) ;
};

// Operator ()                                                                                                                                                                                  
bool TLVP_PtSort::operator() (const JetAK5 & x, const JetAK5 & y){
  return x.Momentum_.Pt () < y.Momentum_.Pt () ;
};

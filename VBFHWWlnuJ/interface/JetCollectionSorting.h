#ifndef JetCollectionSorting_h
#define JetCollectionSorting_h

#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>


#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#ifndef JetAK5_h
#define JetAK5_h

class JetAK5 {

 public :

  // default constructor
  JetAK5();

  JetAK5(const int & position, const std::string & NameCollection, TLorentzVector & Vect) ;

  // default de-constructor
  ~JetAK5();

  int position_ ;
  std::string NameCollection_ ;
  TLorentzVector Momentum_ ;

  
} ;

#endif


#ifndef TLVP_EtaSort_h
#define TLVP_EtaSort_h

class TLVP_EtaSort : public std::binary_function<int,int,bool>{

 public: 
  
  // default constructor
  TLVP_EtaSort(){};

  // default de-constructor
  ~TLVP_EtaSort(){};

  // Operator ()
  bool operator() (const JetAK5 & x, const JetAK5 & y) ;

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
  bool operator() (const JetAK5 & x, const JetAK5 & y);

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
  bool operator() (const JetAK5 & x, const JetAK5 & y);
  
} ;

#endif


#endif

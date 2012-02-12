///---- class for PU MC scaling ----
#ifndef __PUclass__
#define __PUclass__

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include <vector>
#include <iostream>
#include <fstream>


class PUclass : public TObject {
 
 public:
  PUclass();
  ~PUclass();
  
  std::vector<double> PUWeight;
   
  double getPUWeight(int);
  void setPUWeight(const std::vector<double>&);
  void push_back(double);
  void Write(std::string nameoutfile);
  
  ClassDef(PUclass, 1);
  
};


#endif


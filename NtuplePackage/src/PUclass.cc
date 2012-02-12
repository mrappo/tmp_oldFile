#include "PUclass.h"

ClassImp(PUclass)

PUclass::PUclass(){};

PUclass::~PUclass(){};

double PUclass::getPUWeight(int it){
 if (it > (PUWeight.size()-1)) return PUWeight.at(PUWeight.size()-1);
 else return PUWeight.at(it);
};


void PUclass::setPUWeight(const std::vector<double>& in){
 for (int i=0; i< in.size(); i++){
  push_back(in.at(i));
 }
};


void PUclass::push_back(double value){
 PUWeight.push_back(value);
};


void PUclass::Write(std::string nameoutfile){
 ofstream myfile;
 size_t found;
 found = nameoutfile.find(".");
  
 std::string nameoutfileAll = nameoutfile;
 if (found != std::string::npos){
  nameoutfile.erase (found,std::string::npos);
 }
 myfile.open (nameoutfileAll.c_str());
 myfile << "double " << nameoutfile << "(int PU){" << std::endl;
 
 for (int i=0; i<PUWeight.size(); i++){
  myfile << " if (PU == " << i << ") return (" << getPUWeight(i) << ");" << std::endl;  
 }
 myfile << " return 1; " << std::endl;
 myfile << "}" << std::endl;
 
 myfile.close(); 
};





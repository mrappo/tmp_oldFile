#ifndef TKinFitterUtils_h
#define TKinFitterUtils_h

#include <iostream>
#include <cmath>
#include <cstdlib>


double ErrEt (double Et, double Eta); 
double ErrEta(double Et, double Eta); 
double ErrPhi(double Et, double Eta);

bool electronResolution(const double et, const double eta, double& etRes, double& etaRes, double& phiRes);
bool muonResolution    (const double et, const double eta, double& etRes, double& etaRes, double& phiRes);
bool PFMETResolution   (const double et, double& etRes, double& phiRes);

inline double square(const double x) {return x*x;}

#endif

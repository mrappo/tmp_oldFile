#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>



//
// find kinematic quantities for W+W -> u+v + j + j
//

// find the cross product of two vectors and sin(theta) between them
void dg_cross2(TVector3 &cross, double &sintheta, 
               TLorentzVector p1, TLorentzVector p2);

// inverse euler angle rotations on daughter b from direction of parent a
TLorentzVector dgieuler( TLorentzVector parent, TLorentzVector daughter);



////////////////////////////////

// does lorentz trans by beta, gamma (sense ikey)
// on p vector (px,py,pz,e)
TLorentzVector  dgloren( TLorentzVector p, double b, 
                         double g, double ikey);


/// routine to find cm decay angle of p1 in CM system p1+p2
// returns the cosine of the Jackson angle
double JacksonAngle( TLorentzVector p1, TLorentzVector p2);


void dg_kin_Wuv_Wjj( TLorentzVector pu, TLorentzVector pv, 
                     TLorentzVector pj1, TLorentzVector pj2, 
		     float &cosphipl, float &ctuv, float &ctjj);

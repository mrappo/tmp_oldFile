#ifndef METzCalculator_h
#define METzCalculator_h


#include<iostream>

#include"TMath.h"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

template<typename T> 
 class METzCalculator {
 
    public:

      /// constructor
     
      METzCalculator();
 
      virtual ~METzCalculator();

      void SetMET(T MET) {
	  MET_ = MET;
	}

      void SetLepton(T lepton) {
	  lepton_ = lepton;
	}

      void SetLeptonType(std::string leptonName) {
	  if(leptonName == "muon")      leptonMass_ = 0.105658367;
          if(leptonName == "electron")  leptonMass_ = 0.00051099891;
          if(leptonName == "tau")       leptonMass_ = 1.77682;
	}

        /// Calculate METz
	/// options to choose roots from quadratic equation:
	/// type = 0 (defalut): if real roots, pick the one nearest to
	///                     the lepton Pz except when the Pz so chosen
	///                     is greater than 300 GeV in which case pick
	///                     the most central root.
	/// type = 1: if real roots, choose the one closest to the lepton Pz
	///           if complex roots, use only the real part.
	/// type = 2: if real roots, choose the most central solution.
	///           if complex roots, use only the real part.
	/// type = 3: if real roots, pick the largest value of the cosine*

      double Calculate(int type=0) ;
 
      /// check for complex root

      bool IsComplex() const { return isComplex_; };
      double getOther() const { return otherSol_; };
      double getPtneutrino(int option=1) const { 
	  if ( option == 1 ) return newPtneutrino1_;
	  else return newPtneutrino2_;
	}
      void Print() {
		std::cout << " METzCalculator: pxmu = " << lepton_.Px() << " pzmu= " << lepton_.Pz() << std::endl;
		std::cout << " METzCalculator: pxnu = " << MET_.Px() << " pynu= " << MET_.Py() << std::endl;
	}
	
  private:
	
	bool isComplex_;
        T lepton_;
	T MET_;
	double otherSol_;
	double leptonMass_;
	double newPtneutrino1_;
	double newPtneutrino2_;
};

#endif

template <typename T>
METzCalculator<T>::METzCalculator(){
  
  isComplex_ = false;
  otherSol_ = 0.;
  leptonMass_ = 0.105658367;
  newPtneutrino1_ = -1;
  newPtneutrino2_ = -1;
}

template <typename T>
METzCalculator<T>::~METzCalculator(){};

template <typename T>
double METzCalculator<T>::Calculate(int type){


  double M_W  = 80.4;
  double M_mu =  leptonMass_;
  double emu = lepton_.E();
  double pxmu = lepton_.Px(); 
  double pymu = lepton_.Py();
  double pzmu = lepton_.Pz();
  double pxnu = MET_.Px();
  double pynu = MET_.Py();
  double pznu = 0.;
  otherSol_ = 0.;

  double a = M_W*M_W - M_mu*M_mu + 2.0*pxmu*pxnu + 2.0*pymu*pynu;
  double A = 4.0*(emu*emu - pzmu*pzmu);
  double B = -4.0*a*pzmu;
  double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;

  double tmproot = B*B - 4.0*A*C;

  if (tmproot<0) {
        isComplex_= true;
        pznu = - B/(2*A); // take real part of complex roots                                                                                                                             
        otherSol_ = pznu;
        double pnu = MET_.E();
        double Delta = (M_W*M_W - M_mu*M_mu);
        double alpha = (pxmu*pxnu/pnu + pymu*pynu/pnu);
        double ptnu = TMath::Sqrt( pxnu*pxnu + pynu*pynu); // old                                                                                                                        
        double AA = 4.*pzmu*pzmu - 4*emu*emu + 4*alpha*alpha;
        double BB = 4.*alpha*Delta;
        double CC = Delta*Delta;

        double tmpdisc = BB*BB - 4.0*AA*CC;
        double tmpsolpt1 = (-BB + TMath::Sqrt(tmpdisc))/(2.0*AA);
        double tmpsolpt2 = (-BB - TMath::Sqrt(tmpdisc))/(2.0*AA);

        if ( fabs( tmpsolpt1 - ptnu ) < fabs( tmpsolpt2 - ptnu) ) { newPtneutrino1_ = tmpsolpt1; newPtneutrino2_ = tmpsolpt2;}
         else { newPtneutrino1_ = tmpsolpt2; newPtneutrino2_ = tmpsolpt1; }

  }
  else {
         isComplex_ = false;
         double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
         double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);

         if (type == 0 ) {
                             // two real roots, pick the one closest to pz of muon                                                                                                            
                             if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2; otherSol_ = tmpsol1;}
                             else { pznu = tmpsol1; otherSol_ = tmpsol2; }
                             // if pznu is > 300 pick the most central root                                                                                                                   
                             if ( pznu > 300. ) {
                               if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pznu = tmpsol1; otherSol_ = tmpsol2; }
                               else { pznu = tmpsol2; otherSol_ = tmpsol1; }
                             }
         }
         if (type == 1 ) {
                           // two real roots, pick the one closest to pz of muon                                                                                                            
                           if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2; otherSol_ = tmpsol1; }
                           else {pznu = tmpsol1; otherSol_ = tmpsol2; }
         }
         if (type == 2 ) {
                            // pick the most central root.                                                                                                                                   
                            if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pznu = tmpsol1; otherSol_ = tmpsol2; }
                            else { pznu = tmpsol2; otherSol_ = tmpsol1; }
         }
         if (type == 3 ) {
                             // pick the largest value of the cosine                                                                                                                          
                             TVector3 p3w, p3mu;
                             p3w.SetXYZ(pxmu+pxnu, pymu+pynu, pzmu+ tmpsol1);
                             p3mu.SetXYZ(pxmu, pymu, pzmu );

                             double sinthcm1 = 2.*(p3mu.Perp(p3w))/M_W;
                             p3w.SetXYZ(pxmu+pxnu, pymu+pynu, pzmu+ tmpsol2);
                             double sinthcm2 = 2.*(p3mu.Perp(p3w))/M_W;

                             double costhcm1 = TMath::Sqrt(1. - sinthcm1*sinthcm1);
                             double costhcm2 = TMath::Sqrt(1. - sinthcm2*sinthcm2);

                             if ( costhcm1 > costhcm2 ) { pznu = tmpsol1; otherSol_ = tmpsol2; }
                             else { pznu = tmpsol2;otherSol_ = tmpsol1; }
        }

      }
        
  return pznu;
}

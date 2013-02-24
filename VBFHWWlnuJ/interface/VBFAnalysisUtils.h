#ifndef VBFAnalysisUtils_h
#define VBFAnalysisUtils_h

#include "ntpleUtils.h"

#include "TFitConstraintMGaus.h"
#include "TFitConstraintEp.h"
#include "TFitParticleEtEtaPhi.h"
#include "TFitParticleCart.h"
#include "TKinFitter.h"

#include "TLorentzVector.h"
#include "Resolution.h"

void  doKinematicFit(int fflage, const TLorentzVector & mup, const TLorentzVector & nvp, const TLorentzVector & ajp, const TLorentzVector & bjp,
				 TLorentzVector & fit_mup, TLorentzVector & fit_nvp, TLorentzVector & fit_ajp, TLorentzVector & fit_bjp, float & fit_chi2,
                		  int  & fit_NDF, int & fit_status, const std::string & LeptonName);


void calculateAngles  (TLorentzVector& thep4M11, TLorentzVector& thep4M12, TLorentzVector& thep4M21, TLorentzVector& thep4M22, double& costheta1, double& costheta2,
                       double& phi, double& costhetastar, double& phistar1, double& phistar2);

#endif

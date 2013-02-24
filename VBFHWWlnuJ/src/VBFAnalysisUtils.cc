#include "VBFAnalysisUtils.h"

#include "TVector3.h"
#include "TMatrixD.h"

void doKinematicFit(int fflage, const TLorentzVector  & mup, const TLorentzVector & nvp, const TLorentzVector & ajp, const TLorentzVector & bjp,
			        TLorentzVector  & fit_mup, TLorentzVector & fit_nvp, TLorentzVector & fit_ajp, TLorentzVector & fit_bjp,
		    float  & fit_chi2, int & fit_NDF, int & fit_status, const std::string & LeptonName ){

  fit_status = false;
  Resolution* resolution      = new Resolution();

  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  TMatrixD m3(3,3);
  TMatrixD m4(3,3);
  m1.Zero();
  m2.Zero();
  m3.Zero();
  m4.Zero();

  double etRes, etaRes, phiRes;

  // lepton resolution                                                                                                                                                                       
  const TLorentzVector lepton   = mup;
  if(LeptonName == "electron") {
    fit_status = resolution->electronResolution(lepton.Et(), lepton.Eta(), etRes, etaRes, phiRes);
    if(!fit_status) return ;
  } else {
    fit_status = resolution->muonResolution(    lepton.Et(), lepton.Eta(), etRes, etaRes, phiRes);
    if(!fit_status) return ;
  }

  m1(0,0) = resolution->square(etRes);
  m1(1,1) = resolution->square(etaRes);
  m1(2,2) = resolution->square(phiRes);

  // MET resolution                                                                                                                                                                          
  fit_status = resolution->PFMETResolution(     nvp.Et(),            etRes, etaRes, phiRes);
  if(!fit_status) return ;
  m2(0,0) = resolution->square(etRes);
  m2(1,1) = 0.01; // resolution->square(etaRes)                                                                                                                                              
  m2(2,2) = resolution->square(phiRes);

  // Leading Jet resolution                                                                                                                                                                  
  fit_status = resolution->udscPFJetResolution( ajp.Et(), ajp.Eta(), etRes, etaRes, phiRes);
  if(!fit_status) return ;
  m3(0,0) = resolution->square(etRes);
  m3(1,1) = resolution->square(etaRes);
  m3(2,2) = resolution->square(phiRes);

  // Trailing Jet resolution                                                                                                                                                                  
  fit_status = resolution->udscPFJetResolution( bjp.Et(), bjp.Eta(), etRes, etaRes, phiRes);
  if(!fit_status) return ;
  m4(0,0) = resolution->square(etRes);
  m4(1,1) = resolution->square(etaRes);
  m4(2,2) = resolution->square(phiRes);

  TLorentzVector tmp_mup = mup;
  TLorentzVector tmp_nvp = nvp;
  TLorentzVector tmp_ajp = ajp;
  TLorentzVector tmp_bjp = bjp;

  // Fit Particle                                                                                                                                                                            
  TFitParticleEtEtaPhi* particle1 = new TFitParticleEtEtaPhi( "Lepton",   "Lepton",   &tmp_mup,    &m1 );
  TFitParticleEtEtaPhi* particle2 = new TFitParticleEtEtaPhi( "Neutrino", "Neutrino", &tmp_nvp,    &m2 );
  TFitParticleEtEtaPhi* particle3 = new TFitParticleEtEtaPhi( "Jeta",     "Jeta",     &tmp_ajp,    &m3 );
  TFitParticleEtEtaPhi* particle4 = new TFitParticleEtEtaPhi( "Jetb",     "Jetb",     &tmp_bjp,    &m4 );

  // Leptonic WMass Constraint 
  TFitConstraintMGaus* mCons1 = new TFitConstraintMGaus( "W1MassConstraint", "W1Mass-Constraint", 0, 0 , 80.399, 2.085);
  mCons1->addParticles1( particle1, particle2 );

  // Hadronic WMass Constraint 
  TFitConstraintMGaus* mCons2 = new TFitConstraintMGaus( "W2MassConstraint", "W2Mass-Constraint", 0, 0 , 80.399, 2.085);
  mCons2->addParticles1( particle3, particle4 );

  // Conservation of transverse momentum constraint
  TFitConstraintEp *pxCons = new TFitConstraintEp( "PxConstraint", "Px-Constraint", 0, TFitConstraintEp::pX , (mup+nvp+ajp+bjp).Px() );
  pxCons->addParticles( particle1, particle2, particle3, particle4 );

  TFitConstraintEp *pyCons = new TFitConstraintEp( "PyConstraint", "Py-Constraint", 0, TFitConstraintEp::pY , (mup+nvp+ajp+bjp).Py() );
  pyCons->addParticles( particle1, particle2, particle3, particle4 );

  TKinFitter* fitter = new TKinFitter("fitter", "fitter");

  // Only W Lep + W Had  Mass contraint
  if(fflage == 1 ){

    fitter->addMeasParticle( particle1 );
    fitter->addMeasParticle( particle2 );
    fitter->addMeasParticle( particle3 );
    fitter->addMeasParticle( particle4 );
    fitter->addConstraint( mCons1 );
    fitter->addConstraint( mCons2 );

  }
  else if(fflage == 2 ){ // W Lep + W Had + Pt conservation Constraint

    fitter->addMeasParticle( particle1 );
    fitter->addMeasParticle( particle2 );
    fitter->addMeasParticle( particle3 );
    fitter->addMeasParticle( particle4 );
    fitter->addConstraint( pxCons );
    fitter->addConstraint( pyCons );
    fitter->addConstraint( mCons1 );
    fitter->addConstraint( mCons2 );
 }
 else if(fflage == 3 ){ // Only W Had constraint

      fitter->addMeasParticle( particle3 );
      fitter->addMeasParticle( particle4 );
      fitter->addConstraint( mCons2 );

 }else return ;

  //Set convergence criteria                                                                                                                                                                 

  fitter->setMaxNbIter( 50 );
  fitter->setMaxDeltaS( 1e-2 );
  fitter->setMaxF( 1e-1 );
  fitter->setVerbosity(1);
  fitter->fit();

  //Return the kinematic fit results                                                                                                                                                         
  fit_status   = fitter->getStatus();
  fit_chi2     = fitter->getS();
  fit_NDF      = fitter->getNDF();
  fit_mup      = *(particle1->getCurr4Vec());
  fit_nvp      = *(particle2->getCurr4Vec());
  fit_ajp      = *(particle3->getCurr4Vec());
  fit_bjp      = *(particle4->getCurr4Vec());

  if(fitter->getStatus() == 0) { fit_status = true; } else { fit_status = false; }
  delete resolution;
  delete particle1;
  delete particle2;
  delete particle3;
  delete particle4;
  delete mCons1;
  delete mCons2;
  delete pxCons;
  delete pyCons;
  delete fitter;

}


void calculateAngles( TLorentzVector& thep4M11, TLorentzVector& thep4M12, TLorentzVector& thep4M21, TLorentzVector& thep4M22, double& costheta1, double& costheta2,
                      double& phi, double& costhetastar, double& phistar1, double& phistar2){


  TLorentzVector thep4H = thep4M11 + thep4M12 + thep4M21 + thep4M22;
  TLorentzVector thep4Z1 = thep4M11 + thep4M12;
  TLorentzVector thep4Z2 = thep4M21 + thep4M22;

  double norm;

  TVector3 boostX = -(thep4H.BoostVector());

  TLorentzVector thep4Z1inXFrame( thep4Z1 );
  TLorentzVector thep4Z2inXFrame( thep4Z2 );

  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );

  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );

  // Calculate cos theta star
  costhetastar = theZ1X_p3.CosTheta();

  ///////////////////////////////////////////////                                                                                                                                            
  // check for z1/z2 convention, redefine all 4 vectors with convention                                                                                                                      
  ///////////////////////////////////////////////                                                                                                                                            

  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;

  // now helicity angles................................                                                                                                                                     
  // ...................................................                                                                                                                                     

  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);

  //find the decay axis                                                                                                                                                                      
  /////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);                                                                                                                                               
  TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;

  //boost daughters of z2                                                                                                                                                                    
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);

  //create z and y axes                                                                                                                                                                      
  /////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));                                                                                                                    
  TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
  TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
  TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);

  //caculate theta1                                                                                                                                                                          
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
  TVector3 unitM11 = p3M11.Unit();
  double x_m11 = unitM11.Dot(unitx_1); double y_m11 = unitM11.Dot(unity_1); double z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);

  costheta1 = M11_Z1frame.CosTheta();

  phi = M11_Z1frame.Phi();

  //set axes for other system                                                                                                                                                                
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
  norm = 1/(unitx_2.Mag());
  unitx_2*=norm;

  //boost daughters of z2                                                                                                                                                                    
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
  TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
  TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);

  //calcuate theta2                                                                                                                                                                          
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
  TVector3 unitM21 = p3M21.Unit();
  double x_m21 = unitM21.Dot(unitx_2); double y_m21 = unitM21.Dot(unity_2); double z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);

  costheta2 = M21_Z2frame.CosTheta();

  // calculate phi                                                                                                                                                                           
  //calculating phi_n                                                                                                                                                                        

  TLorentzVector n_p4Z1inXFrame( p4Z1 );
  TLorentzVector n_p4M11inXFrame( p4M11 );
  n_p4Z1inXFrame.Boost( boostX );
  n_p4M11inXFrame.Boost( boostX );
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();
  TVector3 n_unitz_1( n_p4Z1inXFrame_unit );


  TVector3 n_unity_1 = n_unitz_1.Cross( n_p4M11inXFrame_unit );
  TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );

  TLorentzVector n_p4M21inXFrame( p4M21 );
  n_p4M21inXFrame.Boost( boostX );
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
  //rotate into other plane                                                                                                                                                                  
  TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );

  TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
  TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );

  // negative sign is for arrow convention in paper                                                                                                                                          
  phistar1 = (n_p4PartoninXFrame_unitprime.Phi());

  // and the calculate phistar2                                                                                                                                                              
  TLorentzVector n_p4Z2inXFrame( p4Z2 );
  n_p4Z2inXFrame.Boost( boostX );
  TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
  TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
  TVector3 n_unity_2 = n_unitz_2.Cross( n_p4M21inXFrame_unit );
  TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );

  TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
  phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());

}


#include "CalibrationUtils.h"
/////////////////////////////////////////////////////////// 

bool CheckxtalIC_EB (TH2F* h_scale_EB,int iPhi, int iEta ){
  if(h_scale_EB->GetBinContent(iPhi,iEta) ==0) return false;
  
  int bx= h_scale_EB->GetNbinsX();
  int by= h_scale_EB->GetNbinsY();

  if((iPhi<bx && h_scale_EB->GetBinContent(iPhi+1,iEta) ==0) || (h_scale_EB->GetBinContent(iPhi-1,iEta)==0 && iPhi>1)) return false;

  if((iEta<by && h_scale_EB->GetBinContent(iPhi,iEta+1) ==0 && iEta!=85 ) || (h_scale_EB->GetBinContent(iPhi,iEta-1)==0 && iEta>1 && iEta!=87)) return false;

  if((iPhi<bx && h_scale_EB->GetBinContent(iPhi+1,iEta+1) ==0 && iEta!=85 && iEta<by) || ( h_scale_EB->GetBinContent(iPhi-1,iEta-1)==0 && iEta>1 && iEta!=87 && iPhi>1)) return false;
 
  if((h_scale_EB->GetBinContent(iPhi+1,iEta-1) ==0 && iEta>1 && iEta!=87 && iPhi<bx) || ( h_scale_EB->GetBinContent(iPhi-1,iEta+1)==0 && iPhi>1 && iEta!=85 && iEta<by )) return false;

  return true;
}

/////////////////////////////////////////////////////////

bool CheckxtalTT_EB (int iPhi, int iEta, const std::vector<std::pair<int,int> >& TT_centre ){
 for(unsigned int k =0; k<TT_centre.size(); k++){
   if(fabs(iPhi-TT_centre.at(k).second)<5 && fabs(iEta-86-TT_centre.at(k).first)<5) return false;
 }
 return true;
}

//////////////////////////////////////////////////////////////////
void InitializeDeadTT_EB(std::vector<std::pair<int,int> >& TT_centre){
 TT_centre.push_back(std::pair<int,int> (58,49));
 TT_centre.push_back(std::pair<int,int> (53,109));
 TT_centre.push_back(std::pair<int,int> (8,114));
 TT_centre.push_back(std::pair<int,int> (83,169));
 TT_centre.push_back(std::pair<int,int> (53,174));
 TT_centre.push_back(std::pair<int,int> (63,194));
 TT_centre.push_back(std::pair<int,int> (83,224));
 TT_centre.push_back(std::pair<int,int> (73,344));
 TT_centre.push_back(std::pair<int,int> (83,358));
 TT_centre.push_back(std::pair<int,int> (-13,18));
 TT_centre.push_back(std::pair<int,int> (-18,23));
 TT_centre.push_back(std::pair<int,int> (-8,53));
 TT_centre.push_back(std::pair<int,int> (-3,63));
 TT_centre.push_back(std::pair<int,int> (-53,128));
 TT_centre.push_back(std::pair<int,int> (-53,183));
 TT_centre.push_back(std::pair<int,int> (-83,193));
 TT_centre.push_back(std::pair<int,int> (-74,218));
 TT_centre.push_back(std::pair<int,int> (-8,223));
 TT_centre.push_back(std::pair<int,int> (-68,303));
 TT_centre.push_back(std::pair<int,int> (-43,328));
}

///////////////////////////////////////////////////////////////////
void NormalizeIC_EB(TH2F* h_scale_EB, TH2F* hcmap, const std::vector< std::pair<int,int> > & TT_centre, bool skip){

 /// Mean over phi corrected skipping dead channel 
  for (int iEta = 1 ; iEta < h_scale_EB->GetNbinsY()+1; iEta ++){
   float SumIC = 0;
   int numIC = 0;
   
   for(int iPhi = 1 ; iPhi < h_scale_EB->GetNbinsX()+1 ; iPhi++){
    bool isGood = CheckxtalIC_EB(h_scale_EB,iPhi,iEta);
    bool isGoodTT = CheckxtalTT_EB(iPhi,iEta,TT_centre);
 
     if(isGood && isGoodTT){
      SumIC = SumIC + h_scale_EB->GetBinContent(iPhi,iEta);
      numIC ++ ;
     }
   }

   ///fede: skip bad channels and bad TTs
   for (int iPhi = 1; iPhi< h_scale_EB->GetNbinsX()+1  ; iPhi++){ 
     if(numIC==0 || SumIC==0) continue;
     if(!skip){ hcmap->SetBinContent(iPhi,iEta,h_scale_EB->GetBinContent(iPhi,iEta)/(SumIC/numIC)); continue;}
     bool isGood = CheckxtalIC_EB(h_scale_EB,iPhi,iEta);
     bool isGoodTT = CheckxtalTT_EB(iPhi,iEta,TT_centre);
     if (!isGood || !isGoodTT) continue;
     hcmap->SetBinContent(iPhi,iEta,h_scale_EB->GetBinContent(iPhi,iEta)/(SumIC/numIC));
   }
  }
}

///////////////////////////////////////////////////////////////////////
void BookSpreadHistos_EB(TH2F* hcmap, TH1F **hspreadEtaFold, const int & ringGroupSize,const int & nEtaRing){

 char hname[100];
 int nStep = 0;
 int nbins = 500;

 /// Spread histos folding EB+ and EB-
 for (int jbin = 1; jbin < hcmap-> GetNbinsY()+1; jbin++){
  if (jbin < nEtaRing+1 && (jbin-1)%ringGroupSize == 0 ) {
      nStep++;
      sprintf(hname,"hspread_ringGroup_ietaFolded%02d",nStep);
      hspreadEtaFold[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
   }
   if (jbin > nEtaRing+1 && (jbin-2)%ringGroupSize == 0 ) {
      nStep++;
   }

   for (int ibin = 1; ibin < hcmap-> GetNbinsX()+1; ibin++){
      float ic = hcmap->GetBinContent(ibin,jbin);
   if (ic>0 && ic<2 )    {
        if (nStep <= nEtaRing) hspreadEtaFold[nStep-1]->Fill(ic);
        else                   hspreadEtaFold[nEtaRing*2-nStep]->Fill(ic);
      }
    }
  }
}  

/////////////////////////////////////////////////////////////////////////
void BookSpreadStatHistos_EB(TH2F* hcmap2,TH2F* hcmap3,TH1F **hstatprecisionEtaFold,const int & ringGroupSize,const int & nEtaRing){

 char hname[100];
 int nStep = 0;
 int nbins = 500;

  for (int jbin = 1; jbin < hcmap2-> GetNbinsY()+1; jbin++){
    if (jbin < nEtaRing+1 && (jbin-1)%ringGroupSize == 0 ) {
       nStep++;
       sprintf(hname,"hstatprecision_ringGroup_ietaFolded%02d",nStep);
       hstatprecisionEtaFold[nStep-1]= new TH1F(hname, hname, nbins,-0.5,0.5);
    }
    if (jbin > nEtaRing+1 && (jbin-2)%ringGroupSize == 0 ) {
       nStep++;
    }

    for (int ibin = 1; ibin < hcmap2-> GetNbinsX()+1; ibin++){
       float ic1 = hcmap2->GetBinContent(ibin,jbin);
       float ic2 = hcmap3->GetBinContent(ibin,jbin);
    if (ic1>0 && ic1<2 && ic2>0 && ic2 <2 )    {
        if (nStep <= nEtaRing) hstatprecisionEtaFold[nStep-1]->Fill((ic1-ic2)/(ic1+ic2));
        else             hstatprecisionEtaFold[nEtaRing*2-nStep]->Fill((ic1-ic2)/(ic1+ic2));
       }
     }
   }

}
///////////////////////////////////////////////////////////////////////////////
void PhiProfile(TGraphErrors *phiProjection, TGraphErrors **MomentumScale, TH2F* hcmap){

  std::vector<double> vectSum;
  std::vector<double> vectCounter;
 
  vectCounter.assign(MomentumScale[0]->GetN(),0.);
  vectSum.assign(MomentumScale[0]->GetN(),0.);

  for(int iPhi =1; iPhi<hcmap->GetNbinsX()+1; iPhi++){
   for(int iEta =1; iEta<hcmap->GetNbinsY()+1; iEta++){

    if(hcmap->GetBinContent(iPhi,iEta)==0) continue;
    unsigned int Phi = int((iPhi-1)/(360./MomentumScale[0]->GetN()));
    if(Phi != vectCounter.size()){
      vectSum.at(Phi)=vectSum.at(Phi)+hcmap->GetBinContent(iPhi,iEta);
      vectCounter.at(Phi)=vectCounter.at(Phi)+1;}
    else{
           vectSum.at(0)=vectSum.at(0)+hcmap->GetBinContent(iPhi,iEta);
           vectCounter.at(0)=vectCounter.at(0)+1;}
           
   }
  }

 for(unsigned int i=0; i<vectCounter.size();i++){
  phiProjection->SetPoint(i,i,vectSum.at(i)/vectCounter.at(i));
  phiProjection->SetPointError(i,0.,0.002);
 }

}

////////////////////////////////////////////////////////////////////////////
void ResidualSpread (TGraphErrors *statprecision, TGraphErrors *Spread, TGraphErrors *Residual){

 for (int i= 0; i < statprecision-> GetN(); i++){
      double spread, espread;
      double stat, estat;
      double residual, eresidual;
      double xdummy,ex;
      Spread->GetPoint(i, xdummy, spread);
      espread = Spread-> GetErrorY(i);
      statprecision->GetPoint(i, xdummy, stat);
      estat = statprecision-> GetErrorY(i);
      ex = statprecision-> GetErrorX(i);
      if (spread > stat ){
	residual  = sqrt( spread*spread - stat*stat );
	eresidual = sqrt( pow(spread*espread,2) + pow(stat*estat,2))/residual;
      }
      else {
	residual = 0;
	eresidual = 0;
      }
      Residual->SetPoint(i,xdummy, residual);
      Residual->SetPointError(i,ex,eresidual);
    }
 }

////////////////////////////////////////////////////////////// 

bool CheckxtalIC_EE(TH2F* h_scale_EE,int ix, int iy, int ir){
  if(h_scale_EE->GetBinContent(ix,iy) ==0) return false;
  
  int bx= h_scale_EE->GetNbinsX();
  int by= h_scale_EE->GetNbinsY();

  if((ix<bx && h_scale_EE->GetBinContent(ix+1,iy) ==0 && (ir!=0 || ir<33)) || (h_scale_EE->GetBinContent(ix-1,iy)==0 && ix>1 && (ir!=0 || ir<33))) return false;

  if((iy<by && h_scale_EE->GetBinContent(ix,iy+1) ==0 && (ir!=0 || ir<33)) || (h_scale_EE->GetBinContent(ix,iy-1)==0 && iy>1 && (ir!=0 || ir<33))) return false;

  if((ix<bx && h_scale_EE->GetBinContent(ix+1,iy+1) ==0 && iy<by && (ir!=0 || ir<33)) || ( h_scale_EE->GetBinContent(ix-1,iy-1)==0 && iy>1 && ix>1 && (ir!=0 || ir<33))) return false;
 
  if((h_scale_EE->GetBinContent(ix+1,iy-1) ==0 && iy>1 && ix<bx && (ir!=0 || ir<33)) || ( h_scale_EE->GetBinContent(ix-1,iy+1)==0 && ix>1 && iy<by && (ir!=0 || ir<33))) return false;

  return true;

}

////////////////////////////////////////////////////////////////
bool CheckxtalTT_EE(int ix, int iy, int ir, const std::vector<std::pair<int,int> >& TT_centre ){
 for( unsigned int k =0; k<TT_centre.size(); k++){
   if(fabs(ix-TT_centre.at(k).first)<5 && fabs(iy-TT_centre.at(k).second)<5) return false;

 }
 return true;
}

//////////////////////////////////////////////////////////////////
void InitializeDeadTTEEP(std::vector<std::pair<int,int> >& TT_centre){
 TT_centre.push_back(std::pair<int,int> (78,78));
 TT_centre.push_back(std::pair<int,int> (83,28));
 TT_centre.push_back(std::pair<int,int> (83,23));
}

/////////////////////////////////////////////////////////////////////
void InitializeDeadTTEEM(std::vector<std::pair<int,int> >& TT_centre){

TT_centre.push_back(std::pair<int,int> (53,28)); 

}


/////////////////////////////////////////////////////////////////////////////
void BookSpreadHistos_EE(TH2F** hcmap, TH1F ***hspread, TH1F **hspreadAll,  TEndcapRings *eRings){

 char hname[100];
 int nbins = 200;

 for (int k = 0; k < 2; k++){         
  for (int iring = 0; iring < 40 ; iring++){
      if (k==0){
       sprintf(hname,"hspreadAll_ring%02d",iring);
       hspreadAll[iring] = new TH1F(hname, hname, nbins,0.,2.);
       sprintf(hname,"hspreadEEM_ring%02d",iring);
       hspread[k][iring] = new TH1F(hname, hname, nbins,0.,2.);
      }
      else{ sprintf(hname,"hspreadEEP_ring%02d",iring);
            hspread[k][iring] = new TH1F(hname, hname, nbins,0.,2.); }
  }
 }
  
 /// spread all distribution, spread for EE+ and EE- and comparison with the MC truth
 for (int k = 0; k < 2 ; k++){
   for (int ix = 1; ix < hcmap[k]->GetNbinsX()+1; ix++){
      for (int iy = 1; iy < hcmap[k]->GetNbinsY()+1; iy++){
	int mybin = hcmap[k] -> FindBin(ix,iy);
        int ring;
	if(k==0)  ring = eRings->GetEndcapRing(ix,iy,-1);
	else      ring = eRings->GetEndcapRing(ix,iy,1);
	float ic = hcmap[k]->GetBinContent(mybin);

 	if (ic>0){
	  hspread[k][ring]->Fill(ic);
          hspreadAll[ring]->Fill(ic);}
      }
    }
  }
}

/////////////////////////////////////////////////////////////////
void BookSpreadStatHistos_EE(TH2F** hcmap2,TH2F** hcmap3, TH1F ***hstatprecision, TH1F **hstatprecisionAll,  TEndcapRings *eRings){

 char hname[100];
 int nbins = 200;

 /// stat precision histos for each EE ring
 for (int k = 0; k < 2; k++){
  for (int iring = 0; iring < 40 ; iring ++){
     if (k==0){ sprintf(hname,"hstatprecisionAll_ring%02d",iring);
               hstatprecisionAll[iring] = new TH1F(hname, hname, nbins,-1.3,1.3);
               sprintf(hname,"hstatprecisionEEM_ring%02d",iring);
               hstatprecision[k][iring] = new TH1F(hname, hname, nbins,-1.3,1.3);}
	else {sprintf(hname,"hstatprecisionEEP_ring%02d",iring);
              hstatprecision[k][iring] = new TH1F(hname, hname, nbins,-1.3,1.3);}
       }
     }
//     
 for (int k = 0; k < 2 ; k++){
   for (int ix = 1; ix < hcmap2[k]->GetNbinsX()+1; ix++){
	for (int iy = 1; iy < hcmap2[k]->GetNbinsY()+1; iy++){
 	  int mybin = hcmap2[k] -> FindBin(ix,iy);
          int ring;
	  if(k==0)  ring = eRings->GetEndcapRing(ix,iy,-1);
	  else      ring = eRings->GetEndcapRing(ix,iy,1);
	  float ic1 = hcmap2[k]->GetBinContent(mybin);
	  float ic2 = hcmap3[k]->GetBinContent(mybin);
          if (ic1>0 && ic2 >0){
	    hstatprecision[k][ring]->Fill((ic1-ic2)/(ic1+ic2)); /// sigma (diff/sum) gives the stat. precision on teh entire sample
            hstatprecisionAll[ring]->Fill((ic1-ic2)/(ic1+ic2));
	  }
	}
      }
    }
 
}

///////////////////////////////////////////////////////////////////////
void PhiProfileEE(TGraphErrors *phiProjection, TGraphErrors **MomentumScale, TH2F* hcmap,TEndcapRings *eRings, const int & iz){

std::vector<double> vectSum;
std::vector<double> vectCounter;
 
vectCounter.assign(MomentumScale[0]->GetN(),0.);
vectSum.assign(MomentumScale[0]->GetN(),0.);

 for(int ix=1; ix<hcmap->GetNbinsX()+1;ix++){
   for(int iy=1; iy<hcmap->GetNbinsY()+1;iy++){
    if(hcmap->GetBinContent(ix,iy)==0) continue;
      int iPhi = int(eRings->GetEndcapIphi(ix,iy,iz)/(360./MomentumScale[0]->GetN()));
      vectSum.at(iPhi)=vectSum.at(iPhi)+hcmap->GetBinContent(ix,iy);
      vectCounter.at(iPhi)=vectCounter.at(iPhi)+1;
  }
 }


 for(unsigned int i=0; i<vectCounter.size();i++)
  phiProjection->SetPoint(i,int(i*(360./MomentumScale[0]->GetN())),vectSum.at(i)/vectCounter.at(i));
}

/////////////////////////////////////////////////////////////////////////
void NormalizeIC_EE(TH2F** hcmap, TH2F** hcmap2, const std::vector< std::pair<int,int> > & TT_centre_EEP,const  std::vector< std::pair<int,int> > & TT_centre_EEM, TEndcapRings *eRings){

 std::vector<float> SumIC_Ring_EEP,SumIC_Ring_EEM,Sumxtal_Ring_EEP,Sumxtal_Ring_EEM;
  
 SumIC_Ring_EEP.assign(40,0);
 SumIC_Ring_EEM.assign(40,0);
 Sumxtal_Ring_EEP.assign(40,0);
 Sumxtal_Ring_EEM.assign(40,0);
 
 /// Mean over phi corrected skipping dead channel 
  
 for(int k=0; k<2 ; k++){
    for( int ix = 0; ix < hcmap[k]->GetNbinsX()+1 ; ix++ ){
     for(int iy = 0; iy < hcmap[k]->GetNbinsY()+1 ; iy++ ){
           
       int ring;
       if(k==0)  ring = eRings->GetEndcapRing(ix,iy,-1);
       else      ring = eRings->GetEndcapRing(ix,iy,1);
       
       bool isGood = CheckxtalIC_EE(hcmap[k],ix,iy,ring);
       bool isGoodTT;
       if(k==0) isGoodTT = CheckxtalTT_EE(ix,iy,ring,TT_centre_EEM);
       else isGoodTT = CheckxtalTT_EE(ix,iy,ring,TT_centre_EEP);
 
       if(k!=0 && isGoodTT && isGood ){
        SumIC_Ring_EEP.at(ring) = SumIC_Ring_EEP.at(ring) + hcmap[k]->GetBinContent(ix,iy);
        Sumxtal_Ring_EEP.at(ring) = Sumxtal_Ring_EEP.at(ring) + 1.;
       }
       if(k==0 && isGoodTT && isGood){
              SumIC_Ring_EEM.at(ring) = SumIC_Ring_EEM.at(ring) + hcmap[k]->GetBinContent(ix,iy);
              Sumxtal_Ring_EEM.at(ring) = Sumxtal_Ring_EEM.at(ring) + 1.;
        }
      }
    }
  }

 for(int k=0; k<2 ; k++){
   for( int ix = 0; ix < hcmap[k]->GetNbinsX()+1 ; ix++ ){
     for(int iy = 0; iy < hcmap[k]->GetNbinsY()+1 ; iy++ ){

       int ring;
       if(k==0)  ring = eRings->GetEndcapRing(ix,iy,-1);
       else      ring = eRings->GetEndcapRing(ix,iy,1);
      
       if(k!=0){
          if(ring>33){ hcmap2[k]->Fill(ix,iy,0.);
                     continue;}
          if(Sumxtal_Ring_EEP.at(ring) != 0 && SumIC_Ring_EEP.at(ring)!= 0)
          hcmap2[k]->Fill(ix,iy,hcmap[k]->GetBinContent(ix,iy)/(SumIC_Ring_EEP.at(ring)/Sumxtal_Ring_EEP.at(ring)));
       }
       else{
            if(ring>33){hcmap2[k]->Fill(ix,iy,0.);
                      continue;
                      }
            if(Sumxtal_Ring_EEM.at(ring) != 0 && SumIC_Ring_EEM.at(ring) != 0)
            hcmap2[k]->Fill(ix,iy,hcmap[k]->GetBinContent(ix,iy)/(SumIC_Ring_EEM.at(ring)/Sumxtal_Ring_EEM.at(ring)));
           }
       }
    }
  }
}

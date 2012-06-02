// default constructor, reading the map from file
#include"TEndcapRings.h"

TEndcapRings::TEndcapRings() {
  FILE *fRing;
  fRing = fopen("./CommonTools/eerings.dat","r");
  std::cout << "Inizializing endcap geometry from: eerings.dat" << std::endl;
  int ix,iy,iz,ir;
  while(fscanf(fRing,"(%d,%d,%d) %d \n",&ix,&iy,&iz,&ir) !=EOF ) {
    if (iz<0) iz=0; 
    iEndcapRing[ix][iy][iz] = ir;
  }
  return; 
}

TEndcapRings::~TEndcapRings() { return;}

Int_t TEndcapRings::GetEndcapRing(Int_t ix, Int_t iy, Int_t iz){
  return iEndcapRing[ix][iy][iz];
}

Int_t TEndcapRings::GetEndcapIeta(Int_t ix, Int_t iy, Int_t iz){
  Int_t iSide = iz; 
  if (iSide<0) iSide=0; 
  Int_t iEtaOffset = 86*iz; 
  Int_t iEta = iEtaOffset + iz*iEndcapRing[ix][iy][iSide];
  return iEta;
}

Int_t TEndcapRings::GetEndcapIphi(Int_t ix, Int_t iy, Int_t iz){

  Int_t iSide = iz; 
  Double_t iX = ix-50.;
  Double_t iY = iy-50.;
  Int_t iPhi;
  if(iY>=0) iPhi = 90.+int (TMath::ATan(iY/iX)*360./(2.*TMath::Pi()));
  else iPhi = 270.+int (TMath::ATan(iY/iX)*360./(2.*TMath::Pi()));
  return iPhi;

}

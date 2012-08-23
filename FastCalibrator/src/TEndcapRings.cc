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

int TEndcapRings::GetEndcapRing(int ix, int iy, int iz){
  return iEndcapRing[ix][iy][iz];
}

int TEndcapRings::GetEndcapIeta(int ix, int iy, int iz){
  int iSide = iz; 
  if (iSide<0) iSide=0; 
  int iEtaOffset = 86*iz; 
  int iEta = iEtaOffset + iz*iEndcapRing[ix][iy][iSide];
  return iEta;
}

int TEndcapRings::GetEndcapIphi(int ix, int iy, int iz){

  int iSide = iz; 
  double iX = ix-50.;
  double iY = iy-50.;
  int iPhi;
  if(iY>=0) iPhi = 90 + int(TMath::ATan(iY/iX)*360./(2.*TMath::Pi()));
  else iPhi = 270 + int(TMath::ATan(iY/iX)*360./(2.*TMath::Pi()));
  return iPhi;
}

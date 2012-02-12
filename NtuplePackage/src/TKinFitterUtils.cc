#include "TKinFitterUtils.h"



// pfjet resolutions. taken from AN-2010-371
double ErrEt( double Et, double Eta) {
  
  double InvPerr2;

  double N, S, C, m;
  if(fabs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( fabs(Eta) < 1. ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( fabs(Eta) < 1.5 ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( fabs(Eta) < 2. ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( fabs(Eta) < 2.5 ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
  } else if( fabs(Eta) < 3. ) {
    N = -3.33814;
    S = 0.73360;
    C = 0.;
    m = 0.08264;
  } else if( fabs(Eta) < 5. ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }
 else {
    std::cout << "Not implemented Et Err for eta > 5. Exiting." << std::endl;
    exit(1123);
  }

  // this is the absolute resolution (squared), not sigma(pt)/pt
  // so have to multiply by pt^2, thats why m+1 instead of m-1
  InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;


  return InvPerr2;

}



//pfjet position resolutions taken from AN 2010-104
double ErrEta( double Et, double Eta) {

  double InvPerr2;

  double a,b,c,d,e;
  if( fabs(Eta) < 0.5 ) {
    a = 487.702;
    b = 1.72868;
    c = 0.0297405;
    d = -487.197;
    e = -1.16389;
  } else if( fabs(Eta) < 1.0 ) {
    a = 277.114;
    b = 1.31746;
    c = 0.0232343;
    d = -276.588;
    e = -1.07289;
  } else if( fabs(Eta) < 1.5 ) {
    a = 19.7603;
    b = 0.406775;
    c = 0.0056006;
    d = -19.1144;
    e = -1.24005;
  } else if( fabs(Eta) < 2.0 ) {
    a = 41.55;
    b = 0.556349;
    c = 0.0094941;
    d = -40.8018;
    e = -1.39179;
  } else if( fabs(Eta) < 2.5 ) {
    a = 0.833363;
    b = 0.0786743;
    c = 0.0036199;
    d = 0.0507317;
    e = -1.5492;
  } else if( fabs(Eta) < 3. ) {
    a = 3.4712;
    b = 0.383594;
    c = 2.6831e-7;
    d = -2.93791;
    e = -0.259687;
  } else {
    std::cout << "Not implemented Eta Err for eta > 3. Exiting." << std::endl;
    exit(1123);
  }
  
  InvPerr2 = sqrt(a*a/(Et * Et) + b*b/Et + c*c) + d/Et + e/pow(Et,1.5);
  InvPerr2 *= InvPerr2;

  return InvPerr2;

}



double ErrPhi( double Et, double Eta) {

  double InvPerr2;

  double a,b,c,d,e;
  if( fabs(Eta) < 0.5 ) {
    a = 926.978;
    b = 2.52747;
    c = 0.0304001;
    d = -926.224;
    e = -1.94117;
  } else if( fabs(Eta) < 1.0 ) {
    a = 3.3251e-6;
    b = 0.063941;
    c = 0.0038759;
    d = 0.301932;
    e = -0.825352;
  } else if( fabs(Eta) < 1.5 ) {
    a = 0.38469;
    b = 0.0755727;
    c = 0.0044353;
    d = 0.453887;
    e = -1.8947;
  } else if( fabs(Eta) < 2.0 ) {
    a = 2.9200e-7;
    b = 0.0718389;
    c = 0.0038557;
    d = 0.403668;
    e = -0.62698;
  } else if( fabs(Eta) < 2.5 ) {
    a = 0.0033663;
    b = 0.0880209;
    c = 0.0023084;
    d = 0.214304;
    e = -0.416353;
  } else if( fabs(Eta) < 3. ) {
    a = 11.1957;
    b = 0.643236;
    c = 0.0071142;
    d = -10.7613;
    e = 0.280927;
  } else {
    std::cout << "Not implemented Phi Err for eta > 3. Exiting." << std::endl;
    exit(1123);
  }
  
  InvPerr2 = sqrt(a*a/(Et * Et) + b*b/Et + c*c) + d/Et + e/pow(Et,1.5);
  InvPerr2 *= InvPerr2;

  return InvPerr2;

}



bool electronResolution(const double et, const double eta, double& etRes, double& etaRes, double& phiRes)
{
  // Check that eta is in range

  if(fabs(eta)>2.5)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0;
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0;

  // Set the coefficients according to the eta interval
  // If no eta interval qualifies, return false to signal failure.
  if(0.000<=fabs(eta) && fabs(eta)<0.174) {
    aEt  = 0.01188;   bEt  = 0.045;     cEt  = 0.29;
    aEta = 0.0004763; bEta = 0.00059;   cEta = 0.0;
    aPhi = 0.0;       bPhi = 0.0014437; cPhi = 0.0;
  } else if(0.174<=fabs(eta) && fabs(eta)<0.261) {
    aEt  = 0.01256;   bEt  = 0.0564;   cEt  = 0.0;
    aEta = 0.0003963; bEta = 0.000848; cEta = 0.0;
    aPhi = 8.8e-05;   bPhi = 0.001193; cPhi = 0.0041;
  } else if(0.261<=fabs(eta) && fabs(eta)<0.348) {
    aEt  = 0.01129;  bEt  = 0.0703;   cEt  = 0.0;
    aEta = 0.000348; bEta = 0.00091;  cEta = 0.0;
    aPhi = 9.5e-05;  bPhi = 0.001192; cPhi = 0.00437;
  } else if(0.348<=fabs(eta) && fabs(eta)<0.435) {
    aEt  = 0.01275;   bEt  = 0.0621;  cEt  = 0.0;
    aEta = 0.0003152; bEta = 0.00096; cEta = 0.0;
    aPhi = 5.5e-05;   bPhi = 0.00143; cPhi = 0.00293;
  } else if(0.435<=fabs(eta) && fabs(eta)<0.522) {
    aEt  = 0.01256;   bEt  = 0.0678;   cEt  = 0.0;
    aEta = 0.0003111; bEta = 0.00093;  cEta = 0.0;
    aPhi = 7.4e-05;   bPhi = 0.001391; cPhi = 0.00326;
  } else if(0.522<=fabs(eta) && fabs(eta)<0.609) {
    aEt  = 0.01139;   bEt  = 0.0729;   cEt  = 0.0;
    aEta = 0.0003167; bEta = 0.00088;  cEta = 0.0;
    aPhi = 0.000114;  bPhi = 0.001294; cPhi = 0.00392;
  } else if(0.609<=fabs(eta) && fabs(eta)<0.696) {
    aEt  = 0.01285;   bEt  = 0.0599;   cEt  = 0.0;
    aEta = 0.0003251; bEta = 0.00102;  cEta = 0.0;
    aPhi = 7.8e-05;   bPhi = 0.001452; cPhi = 0.00304;
  } else if(0.696<=fabs(eta) && fabs(eta)<0.783) {
    aEt  = 0.01147;   bEt  = 0.0784;   cEt  = 0.0;
    aEta = 0.0003363; bEta = 0.001;    cEta = 0.0;
    aPhi = 0.000108;  bPhi = 0.001513; cPhi = 0.00293;
  } else if(0.783<=fabs(eta) && fabs(eta)<0.870) {
    aEt  = 0.01374;  bEt  = 0.0761;   cEt  = 0.0;
    aEta = 0.000324; bEta = 0.00106;  cEta = 0.0;
    aPhi = 0.000127; bPhi = 0.001556; cPhi = 0.00294;
  } else if(0.870<=fabs(eta) && fabs(eta)<0.957) {
    aEt  = 0.01431;   bEt  = 0.0754;  cEt  = 0.0;
    aEta = 0.0003081; bEta = 0.001;   cEta = 0.0;
    aPhi = 0.000164;  bPhi = 0.00149; cPhi = 0.00411;
  } else if(0.957<=fabs(eta) && fabs(eta)<1.044) {
    aEt  = 0.01196;   bEt  = 0.1066;   cEt  = 0.0;
    aEta = 0.0003212; bEta = 0.001;    cEta = 0.0;
    aPhi = 0.0001111; bPhi = 0.001933; cPhi = 0.0;
  } else if(1.044<=fabs(eta) && fabs(eta)<1.131) {
    aEt  = 0.01613;   bEt  = 0.1164; cEt   = 0.0;
    aEta = 0.0003348; bEta = 0.0011; cEta  = 0.0;
    aPhi = 0.000164;  bPhi = 0.00195; cPhi = 0.0022;
  } else if(1.131<=fabs(eta) && fabs(eta)<1.218) {
    aEt  = 0.0227;    bEt  = 0.1091;  cEt  = 0.0;
    aEta = 0.0003474; bEta = 0.00109; cEta = 0.0;
    aPhi = 0.000191;  bPhi = 0.00216; cPhi = 0.0026;
  } else if(1.218<=fabs(eta) && fabs(eta)<1.305) {
    aEt  = 0.0158;    bEt  = 0.1718;  cEt  = 0.0;
    aEta = 0.0003354; bEta = 0.00102; cEta = 0.0;
    aPhi = 0.000274;  bPhi = 0.00208; cPhi = 0.0028;
  } else if(1.305<=fabs(eta) && fabs(eta)<1.392) {
    aEt  = 0.0176;   bEt  = 0.1718;   cEt  = 0.0;
    aEta = 0.000332; bEta = 0.00109;  cEta = 0.0;
    aPhi = 0.000253; bPhi = 0.002472; cPhi = 0.0;
  } else if(1.392<=fabs(eta) && fabs(eta)<1.479) {
    aEt  = 0.0077;   bEt  = 0.2288;   cEt  = 0.0;
    aEta = 0.000317; bEta = 0.001049; cEta = 0.0;
    aPhi = 0.000285; bPhi = 0.00255;  cPhi = 0.003;
  } else if(1.479<=fabs(eta) && fabs(eta)<1.653) {
    aEt  = 0.047;     bEt  = 0.158;   cEt  = 0.0;
    aEta = 0.0003479; bEta = 0.0;     cEta = 0.0036;
    aPhi = 0.000333;  bPhi = 0.00277; cPhi = 0.0;
  } else if(1.653<=fabs(eta) && fabs(eta)<1.740) {
    aEt  = 0.0;       bEt  = 0.2;     cEt  = 0.0;
    aEta = 0.0003390; bEta = 0.0004;  cEta = 0.0027;// Values interpolated from neighbors.
    aPhi = 0.00038;   bPhi = 0.00282; cPhi = 0.0;
  } else if(1.740<=fabs(eta) && fabs(eta)<1.830) {
    aEt  = 0.04019;  bEt  = 0.0;     cEt  = 0.0;
    aEta = 0.00033;  bEta = 0.0009;  cEta = 0.0019;
    aPhi = 0.000269; bPhi = 0.00324; cPhi = 0.0;
  } else if(1.830<=fabs(eta) && fabs(eta)<1.930) {
    aEt  = 0.039;    bEt  = 0.048;   cEt  = 0.0;// Values interpolated from neighbors.
    aEta = 0.000348; bEta = 0.00096; cEta = 0.0016;
    aPhi = 0.000271; bPhi = 0.00369; cPhi = 0.0;
  } else if(1.930<=fabs(eta) && fabs(eta)<2.043) {
    aEt  = 0.038;     bEt  = 0.096;  cEt  = 0.0;
    aEta = 0.0003786; bEta = 0.0;    cEta = 0.00424;
    aPhi = 0.00028;   bPhi = 0.0031; cPhi = 0.0;
  } else if(2.043<=fabs(eta) && fabs(eta)<2.172) {
    aEt  = 0.0382;   bEt  = 0.076;   cEt  = 0.28;
    aEta = 0.000389; bEta = 0.00106; cEta = 0.0;
    aPhi = 0.000401; bPhi = 0.0025;  cPhi = 0.0114;
  } else if(2.172<=fabs(eta) && fabs(eta)<2.322) {
    aEt  = 0.035;    bEt  = 0.11;    cEt  = 0.0;
    aEta = 0.000486; bEta = 0.0002;  cEta = 0.0052;
    aPhi = 0.0;      bPhi = 0.00432; cPhi = 0.0088;
  } else if(2.322<=fabs(eta) && fabs(eta)<2.500) {
    aEt  = 0.0354;   bEt  = 0.123; cEt  = 0.1;
    aEta = 0.000568; bEta = 0.0;   cEta = 0.00734;
    aPhi = 0.000671; bPhi = 0.0;   cPhi = 0.0158;
  } else {
    return false;
  }
  
  etRes  = et * (sqrt(square(aEt)  + square(bEt/sqrt(et))  + square(cEt/et)));
  etaRes =       sqrt(square(aEta) + square(bEta/sqrt(et)) + square(cEta/et));
  phiRes =       sqrt(square(aPhi) + square(bPhi/sqrt(et)) + square(cPhi/et));
  return true;
}

bool muonResolution(const double et, const double eta, double& etRes, double& etaRes, double& phiRes)
{ 
  // Check that eta is in range

  if(fabs(eta)>2.4)  return false;

  double aEt  = 0.0, bEt  = 0.0, cEt  = 0.0;
  double aEta = 0.0, bEta = 0.0, cEta = 0.0; 
  double aPhi = 0.0, bPhi = 0.0, cPhi = 0.0; 

  // Set the coefficients according to the eta interval
  // If no eta interval qualifies, return false to signal failure.
  if(0.000<=fabs(eta) && fabs(eta)<0.100) {
    aEt  = 0.00475;   bEt  = 0.0002365; cEt  = 0.0;
    aEta = 0.0004348; bEta = 0.001063;  cEta = 0.0;
    aPhi = 6.28e-05;  bPhi = 0.0;       cPhi = 0.004545;
  } else if(0.100<=fabs(eta) && fabs(eta)<0.200) {
    aEt  = 0.00509;   bEt  = 0.0002298; cEt  = 0.0;
    aEta = 0.0004348; bEta = 0.001063;  cEta = 0.0;
    aPhi = 5.53e-05;  bPhi = 0.0;       cPhi = 0.004763;
  } else if(0.200<=fabs(eta) && fabs(eta)<0.300) {
    aEt  = 0.005942;  bEt  = 0.0002138; cEt  = 0.0;
    aEta = 0.0003412; bEta = 0.000857;  cEta = 0.00147;
    aPhi = 5.39e-5;   bPhi = 0.0;       cPhi = 0.004842;
  } else if(0.300<=fabs(eta) && fabs(eta)<0.400) {
    aEt  = 0.006989;  bEt  = 0.0002003; cEt  = 0.0;
    aEta = 0.0003208; bEta = 0.000604;  cEta = 0.00187;
    aPhi = 5.63e-5;   bPhi = 0.0;       cPhi = 0.00494;
  } else if(0.400<=fabs(eta) && fabs(eta)<0.500) {
    aEt  = 0.007227;  bEt  = 0.0001996; cEt  = 0.0;
    aEta = 0.0002908; bEta = 0.000733;  cEta = 0.00151;
    aPhi = 5.58e-5;   bPhi = 0.0;       cPhi = 0.00501;
  } else if(0.500<=fabs(eta) && fabs(eta)<0.600) {
    aEt  = 0.007528; bEt  = 0.0001935; cEt  = 0.0;
    aEta = 0.000289; bEta = 0.00076;   cEta = 0.00154;
    aPhi = 5.65e-5;  bPhi = 0.0;       cPhi = 0.005082;
  } else if(0.600<=fabs(eta) && fabs(eta)<0.700) {
    aEt  = 0.007909; bEt  = 0.0001863; cEt  = 0.0;
    aEta = 0.000309; bEta = 0.000667;  cEta = 0.00194;
    aPhi = 5.58e-5;  bPhi = 0.0;       cPhi = 0.005241;
  } else if(0.700<=fabs(eta) && fabs(eta)<0.800) {
    aEt  = 0.008298;  bEt  = 0.000185; cEt  = 0.0;
    aEta = 0.0002887; bEta = 0.000876; cEta = 0.00179;
    aPhi = 5.97e-5;   bPhi = 0.0;      cPhi = 0.005085;
  } else if(0.800<=fabs(eta) && fabs(eta)<0.900) {
    aEt  = 0.00918;   bEt  = 0.0001911; cEt  = 0.0;
    aEta = 0.0002956; bEta = 0.000752;  cEta = 0.00208;
    aPhi = 5.9e-5;    bPhi = 0.0;       cPhi = 0.005506;
  } else if(0.900<=fabs(eta) && fabs(eta)<1.000) {
    aEt  = 0.01096;   bEt  = 0.0001899; cEt  = 0.0;
    aEta = 0.0002734; bEta = 0.000967;  cEta = 0.00134;
    aPhi = 7.48e-5;   bPhi = 0.0;       cPhi = 0.005443;
  } else if(1.000<=fabs(eta) && fabs(eta)<1.100) {
    aEt  = 0.01262;   bEt  = 0.0001614; cEt  = 0.0;
    aEta = 0.0002831; bEta = 0.000968;  cEta = 0.00166;
    aPhi = 7.81e-5;   bPhi = 0.0;       cPhi = 0.005585;
  } else if(1.100<=fabs(eta) && fabs(eta)<1.200) {
    aEt  = 0.01379;  bEt  = 0.0001618; cEt  = 0.0;
    aEta = 0.000293; bEta = 0.000942;  cEta = 0.002;
    aPhi = 8.19e-5;  bPhi = 0.0;       cPhi = 0.005921;
  } else if(1.200<=fabs(eta) && fabs(eta)<1.300) {
    aEt  = 0.01485;   bEt  = 0.0001574; cEt  = 0.0;
    aEta = 0.0002907; bEta = 0.000832;  cEta = 0.002;
    aPhi = 7.89e-5;   bPhi = 0.00039;   cPhi = 0.00593;
  } else if(1.300<=fabs(eta) && fabs(eta)<1.400) {
    aEt  = 0.0152;    bEt  = 0.0001719; cEt  = 0.0;
    aEta = 0.0002937; bEta = 0.000839;  cEta = 0.00232;
    aPhi = 5.9e-5;    bPhi = 0.000724;  cPhi = 0.005664;
  } else if(1.400<=fabs(eta) && fabs(eta)<1.500) {
    aEt  = 0.01471;   bEt  = 0.0001828; cEt  = 0.0;
    aEta = 0.0002999; bEta = 0.000864;  cEta = 0.00229;
    aPhi = 4.7e-5;    bPhi = 0.000834;  cPhi = 0.00527;
  } else if(1.500<=fabs(eta) && fabs(eta)<1.600) {
    aEt  = 0.01337;   bEt  = 0.0002375; cEt  = 0.0;
    aEta = 0.0003035; bEta = 0.000746;  cEta = 0.00258;
    aPhi = 8.16e-5;   bPhi = 0.000757;  cPhi = 0.005558;
  } else if(1.600<=fabs(eta) && fabs(eta)<1.700) {
    aEt  = 0.01308;   bEt  = 0.000285; cEt  = 0.0;
    aEta = 0.0002967; bEta = 0.000798; cEta = 0.00263;
    aPhi = 6.2e-5;    bPhi = 0.001025; cPhi = 0.00523;
  } else if(1.700<=fabs(eta) && fabs(eta)<1.800) {
    aEt  = 0.01302;   bEt  = 0.0003797; cEt  = 0.0;
    aEta = 0.0003063; bEta = 0.000776;  cEta = 0.00278;
    aPhi = 0.000107;  bPhi = 0.001011;  cPhi = 0.00554;
  } else if(1.800<=fabs(eta) && fabs(eta)<1.900) {
    aEt  = 0.0139;    bEt  = 0.000492; cEt  = 0.0;
    aEta = 0.0003285; bEta = 0.00077;  cEta = 0.00292;
    aPhi = 0.000119 ; bPhi = 0.001163; cPhi = 0.00519;
  } else if(1.900<=fabs(eta) && fabs(eta)<2.000) {
    aEt  = 0.01507;   bEt  = 0.000581; cEt  = 0.0;
    aEta = 0.0003365; bEta = 0.00084;  cEta = 0.00323;
    aPhi = 0.000193;  bPhi = 0.00067;  cPhi = 0.00613;
  } else if(2.000<=fabs(eta) && fabs(eta)<2.100) {
    aEt  = 0.01711;   bEt  = 0.000731; cEt  = 0.0;
    aEta = 0.0003504; bEta = 0.00078;  cEta = 0.00365;
    aPhi = 0.000217;  bPhi = 0.00121;  cPhi = 0.00558;
  } else if(2.100<=fabs(eta) && fabs(eta)<2.200) {
    aEt  = 0.01973;  bEt  = 0.000823; cEt  = 0.0;
    aEta = 0.000381; bEta = 0.00088; cEta = 0.00369;
    aPhi = 0.000283; bPhi = 0.00082;  cPhi = 0.00608;
  } else if(2.200<=fabs(eta) && fabs(eta)<2.300) {
    aEt  = 0.02159;  bEt  = 0.0001052; cEt  = 0.0;
    aEta = 0.00042;  bEta = 0.00097;   cEta = 0.00393;
    aPhi = 0.000304; bPhi = 0.00149;   cPhi = 0.00549;
  } else if(2.300<=fabs(eta) && fabs(eta)<2.400) {
    aEt  = 0.02155;  bEt  = 0.001346; cEt  = 0.0;
    aEta = 0.000403; bEta = 0.00153;  cEta = 0.00403;
    aPhi = 0.000331; bPhi = 0.00183;  cPhi = 0.00585;
  } else {
    return false;
  }
  etRes  = et * (aEt + bEt * et);
  etaRes = sqrt(square(aEta) + square(bEta/sqrt(et)) + square(cEta/et));
  phiRes = sqrt(square(aPhi) + square(bPhi/sqrt(et)) + square(cPhi/et));
  return true;
}



bool PFMETResolution(const double et, double& etRes,double& phiRes)
{
  etRes  = et * (sqrt(square(0.05469) + square(10.549/et)));
  phiRes = sqrt(square(0.164/sqrt(et)) + square(11.068/et));
  return true;
}

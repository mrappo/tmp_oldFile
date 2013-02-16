#include "JetCollectionSorting.h"

std::vector<const lorentzVector *> Jet_Deta_pT (const lorentzVector input [], const int & njets){

  std::vector<const lorentzVector *> input_p ;
  std::vector<const lorentzVector *> output_p ;

  for (int iJet = 0 ; iJet < njets ; ++iJet){

    if(input[iJet].Pt()<=0) continue ;
    input_p.push_back (&(input[iJet])) ;
  }

  if(input_p.size()<2) return output_p;

  std::sort (input_p.begin(), input_p.end(), TLVP_EtaSort() ) ;

  output_p.push_back (input_p.front() ) ;
  output_p.push_back (input_p.back() ) ;

  input_p.erase (input_p.begin(), input_p.begin() + 1) ; // erase first element
  input_p.erase (input_p.end()-1, input_p.end() ) ;

  if(input_p.size()<2) return output_p;

  std::sort (input_p.begin(), input_p.end(), TLVP_PtSort()) ;

  output_p.push_back (input_p.at (0)) ;
  output_p.push_back (input_p.at (1)) ;

  return output_p ;
}


// - FNAL criterion, with embedded cuts:
//     |etaj|<4.5
//     etaj1 x etaj2 < 0
//     |Deta|>3.5
//     mj1j2 > 300

std::vector<const lorentzVector *> JS_FNAL (const lorentzVector input [], int njets){

  std::vector<const lorentzVector *> input_p ;
  std::vector<const lorentzVector *> output_p ;

  for (int iJet = 0 ; iJet < njets ; ++iJet){

    if(input[iJet].Pt()<=0) continue ;
    input_p.push_back (&(input[iJet])) ;
  }

  if(input_p.size()<2) return output_p;

  float maxDeta = 0. ;
  int iJ1 = 0 ;
  int iJ2 = 0 ;

  for (int iJet = 0 ; iJet < njets - 1 ; ++iJet){

    if (fabs (input_p.at(iJet)->Eta()) > 4.7 ) continue ;

      for (int jJet = iJet + 1 ; jJet < njets ; ++jJet) {

          if (fabs (input_p.at(jJet)->Eta()) > 4.7) continue ;

          if (input_p.at(iJet)->Eta() * input_p.at(jJet)->Eta() > 0) continue ;

          float deta = fabs (input_p.at(iJet)->Eta() - input_p.at(jJet)->Eta()) ;

          if (deta < 3.5) continue ;

          lorentzVector sum = *(input_p.at(iJet)) + *(input_p.at(jJet)) ;
          float mjj = sum.M2() ;

          if (mjj < 300) continue ;
          if (deta > maxDeta){
              maxDeta = deta ;
              iJ1 = iJet ;
              iJ2 = jJet ;
            }
       }
  }

  if (maxDeta < 0.1 || iJ1 == iJ2) return output_p ;

  output_p.push_back (input_p.at(iJ1)) ; 
  output_p.push_back (input_p.at(iJ2)) ; 

  if (iJ1 < iJ2) std::swap (iJ1, iJ2) ;

  input_p.erase (input_p.begin() + iJ1) ;
  input_p.erase (input_p.begin() + iJ2) ;

  if(input_p.size()<2) return output_p;

  std::sort (input_p.begin(), input_p.end(), TLVP_PtSort()) ;

  output_p.push_back (input_p.at(0)) ;
  output_p.push_back (input_p.at(1)) ;

  return output_p ;
}

// - two jets with highest Mjj are the VBF ones, then two jets with highest pT are the W ones
std::vector<const lorentzVector *> Jet_Mjj_pT (const lorentzVector input [], int njets){

  std::vector<const lorentzVector *> input_p ;
  std::vector<const lorentzVector *> output_p ;

  for (int iJet = 0 ; iJet < njets ; ++iJet){

   if(input[iJet].Pt()<=0) continue ;
   input_p.push_back (&(input[iJet])) ;

  }

  if(input_p.size()<2) return output_p;

  float maxMjj = 0. ;
  int iJ1 = 0 ;
  int iJ2 = 0 ;

  for (int iJet = 0 ; iJet < njets - 1 ; ++iJet){
      for (int jJet = iJet + 1 ; jJet < njets ; ++jJet) {

	lorentzVector sum = *(input_p.at(iJet)) + *(input_p.at(jJet)) ;
          float mjj = sum.M2() ;

          if (mjj > maxMjj){
              maxMjj = mjj ;
              iJ1 = iJet ;
              iJ2 = jJet ;
          }
      }
  }

  output_p.push_back (input_p.at(iJ1)) ; 
  output_p.push_back (input_p.at(iJ2)) ; 

  if (iJ1 < iJ2) std::swap (iJ1, iJ2) ;

  input_p.erase (input_p.begin() + iJ1) ;
  input_p.erase (input_p.begin() + iJ2) ;

  if(input_p.size()<2) return output_p;

  std::sort (input_p.begin(), input_p.end(), TLVP_PtSort()) ;
 
  output_p.push_back (input_p.at(0)) ;
  output_p.push_back (input_p.at(1)) ;

  return output_p ;
}


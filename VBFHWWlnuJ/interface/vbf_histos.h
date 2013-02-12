/*
collection of histograms for the vbf part of the analysis for the 4lepton and 3lepton cases
*/

#ifndef vbf_histos_h
#define vbf_histos_h

//PG VBF histograms set for 2jW 2jVBF
//PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

//PG assumes the four jets to be passed as vbf, vbf, W, W, not sorted in pt
struct H_VBF_4j
{
  H_VBF_4j (TString baseN) :
    baseName (baseN),
    ntuple (new TNtuple (baseName + "_4j_nt", "", "pt1:pt2:mjj:deta")),
    vbf_deta    (baseName + "_4j_vbf_deta"    , "", 60, 0, 10),
    vbf_dphi    (baseName + "_4j_vbf_dphi"    , "", 20, -3.15, 3.15),
    vbf_mjj     (baseName + "_4j_vbf_mjj"     , "", 100, 0, 1000),
    vbf_detaMjj (baseName + "_4j_vbf_detaMjj" , "", 100, 0, 1000, 60, 0, 6),
    W_deta      (baseName + "_4j_W_deta"      , "", 60, 0, 10),
    W_dphi      (baseName + "_4j_W_dphi"      , "", 20, -3.15, 3.15),
    W_mjj       (baseName + "_4j_W_mjj"       , "", 100, 0, 300),
    ptj1        (baseName + "_4j_ptj1"        , "", 100, 0, 500),
    ptj2        (baseName + "_4j_ptj2"        , "", 100, 0, 500),
    ptj3        (baseName + "_4j_ptj3"        , "", 100, 0, 500),
    ptj4        (baseName + "_4j_ptj4"        , "", 100, 0, 500),
    etaj1       (baseName + "_4j_etaj1"       , "", 60, -5, 5),
    etaj2       (baseName + "_4j_etaj2"       , "", 60, -5, 5),
    etaj3       (baseName + "_4j_etaj3"       , "", 60, -5, 5),
    etaj4       (baseName + "_4j_etaj4"       , "", 60, -5, 5),
    etavbfjout  (baseName + "_4j_etavbfjout"  , "", 60, -5, 5),
    etavbfjin   (baseName + "_4j_etavbfjin"   , "", 60, -5, 5),
    etatj       (baseName + "_4j_etatj"       , "", 60, -5, 5),
    etaWj       (baseName + "_4j_etaWj"       , "", 60, -5, 5)
    {
      cout << "created histos set with name " << baseName << endl ;
    }

  void fill (vector<const lorentzVector *> input)
    {
      int seq[4] = {0, 1, 2, 3} ;
      if (input.at (0)->Pt () < input.at (1)->Pt ()) 
        { seq[0] = 1 ; seq[1] = 0 ;}
      if (input.at (2)->Pt () < input.at (3)->Pt ()) 
        { seq[2] = 3 ; seq[3] = 2 ;}
      ptj1.Fill (input.at (seq[0])->Pt ()) ;
      ptj2.Fill (input.at (seq[1])->Pt ()) ;
      ptj3.Fill (input.at (seq[2])->Pt ()) ;
      ptj4.Fill (input.at (seq[3])->Pt ()) ;

      etaj1.Fill (input.at (seq[0])->Eta ()) ;
      etaj2.Fill (input.at (seq[1])->Eta ()) ;
      etaj3.Fill (input.at (seq[2])->Eta ()) ;
      etaj4.Fill (input.at (seq[3])->Eta ()) ;

      etatj.Fill (input.at (seq[0])->Eta ()) ;
      etatj.Fill (input.at (seq[1])->Eta ()) ;
      etaWj.Fill (input.at (seq[2])->Eta ()) ;
      etaWj.Fill (input.at (seq[3])->Eta ()) ;

      if (fabs (input.at (0)->Eta ()) > fabs (input.at (1)->Eta ()))
        {      
          etavbfjout.Fill (input.at (0)->Eta ()) ;
          etavbfjin.Fill  (input.at (1)->Eta ()) ;
        }
      else
        {      
          etavbfjout.Fill (input.at (1)->Eta ()) ;
          etavbfjin.Fill  (input.at (0)->Eta ()) ;
        }

      vbf_deta.Fill (fabs (input.at (0)->Eta () - input.at (1)->Eta ())) ;
      vbf_dphi.Fill (deltaPhi (input.at (0)->Phi (), input.at (1)->Phi ())) ;
      lorentzVector vbf_sum = *input.at (0) + *input.at (1) ;
      vbf_mjj.Fill (sqrt (vbf_sum.M2 ())) ;
      vbf_detaMjj.Fill (sqrt (vbf_sum.M2 ()), fabs (input.at (0)->Eta () - input.at (1)->Eta ())) ;
      ntuple->Fill (input.at (seq[0])->Pt (), input.at (seq[1])->Pt (), 
                    sqrt (vbf_sum.M2 ()), fabs (input.at (0)->Eta () - input.at (1)->Eta ())) ;
      W_deta.Fill (fabs (input.at (2)->Eta () - input.at (3)->Eta ())) ;
      W_dphi.Fill (deltaPhi (input.at (2)->Phi (), input.at (3)->Phi ())) ;
      lorentzVector W_sum = *input.at (2) + *input.at (3) ;
      W_mjj.Fill (sqrt (W_sum.M2 ())) ;

    }

  void write (TFile & output)
    {
      output.cd () ;
      ntuple->Write () ;
      vbf_deta.Write () ;
      vbf_dphi.Write () ;
      vbf_mjj.Write () ;
      vbf_detaMjj.Write () ;
      W_deta.Write () ;
      W_dphi.Write () ;
      W_mjj.Write () ;
      ptj1.Write () ;
      ptj2.Write () ;
      ptj3.Write () ;
      ptj4.Write () ;
      etaj1.Write () ;
      etaj2.Write () ;
      etaj3.Write () ;
      etaj4.Write () ;
      etavbfjout.Write () ;
      etavbfjin.Write () ;
      etatj.Write () ;
      etaWj.Write () ;
    }

  TString baseName ;
  TNtuple * ntuple ;
  TH1F vbf_deta ;
  TH1F vbf_dphi ;
  TH1F vbf_mjj ;
  TH2F vbf_detaMjj ;
  TH1F W_deta ;
  TH1F W_dphi ;
  TH1F W_mjj ;
  TH1F ptj1 ;
  TH1F ptj2 ;
  TH1F ptj3 ;
  TH1F ptj4 ;
  TH1F etaj1 ;
  TH1F etaj2 ;
  TH1F etaj3 ;
  TH1F etaj4 ;
  TH1F etavbfjout ;
  TH1F etavbfjin  ;
  TH1F etatj ;
  TH1F etaWj ;
} ;


//PG VBF histograms set for 1jW 2jVBF
//PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


struct H_VBF_3j
{
  H_VBF_3j (TString baseN) :
    baseName (baseN),
    ntuple (new TNtuple (baseName + "_3j_nt", "", "pt1:pt2:mjj:deta")),
    vbf_deta    (baseName + "_3j_vbf_deta"    , "", 60, 0, 10),
    vbf_dphi    (baseName + "_3j_vbf_dphi"    , "", 20, -3.15, 3.15),
    vbf_mjj     (baseName + "_3j_vbf_mjj"     , "", 100, 0, 1000),
    vbf_detaMjj (baseName + "_3j_vbf_detaMjj" , "", 100, 0, 1000, 60, 0, 6),
    W_mj        (baseName + "_3j_W_mj"        , "", 100, 0, 300),
    ptj1        (baseName + "_3j_ptj1"        , "", 100, 0, 500),
    ptj2        (baseName + "_3j_ptj2"        , "", 100, 0, 500),
    ptj3        (baseName + "_3j_ptj3"        , "", 100, 0, 500),
    etaj1       (baseName + "_3j_etaj1"       , "", 60, -5, 5),
    etaj2       (baseName + "_3j_etaj2"       , "", 60, -5, 5),
    etaj3       (baseName + "_3j_etaj3"       , "", 60, -5, 5),
    etavbfjout  (baseName + "_3j_etavbfjout"  , "", 60, -5, 5),
    etavbfjin   (baseName + "_3j_etavbfjin"   , "", 60, -5, 5),
    etatj       (baseName + "_3j_etatj"       , "", 60, -5, 5),
    etaWj       (baseName + "_3j_etaWj"       , "", 60, -5, 5)
    {
      cout << "created histos set with name " << baseName << endl ;
    }

  void fill (vector<const lorentzVector *> input)
    {
      int seq[4] = {0, 1, 2} ;
      if (input.at (0)->Pt () < input.at (1)->Pt ()) 
        { seq[0] = 1 ; seq[1] = 0 ;}
      ptj1.Fill (input.at (seq[0])->Pt ()) ;
      ptj2.Fill (input.at (seq[1])->Pt ()) ;
      ptj3.Fill (input.at (seq[2])->Pt ()) ;

      etaj1.Fill (input.at (seq[0])->Eta ()) ;
      etaj2.Fill (input.at (seq[1])->Eta ()) ;
      etaj3.Fill (input.at (seq[2])->Eta ()) ;

      etatj.Fill (input.at (seq[0])->Eta ()) ;
      etatj.Fill (input.at (seq[1])->Eta ()) ;
      etaWj.Fill (input.at (seq[2])->Eta ()) ;

      if (fabs (input.at (0)->Eta ()) > fabs (input.at (1)->Eta ()))
        {      
          etavbfjout.Fill (input.at (0)->Eta ()) ;
          etavbfjin.Fill  (input.at (1)->Eta ()) ;
        }
      else
        {      
          etavbfjout.Fill (input.at (1)->Eta ()) ;
          etavbfjin.Fill  (input.at (0)->Eta ()) ;
        }

      vbf_deta.Fill (fabs (input.at (0)->Eta () - input.at (1)->Eta ())) ;
      vbf_dphi.Fill (deltaPhi (input.at (0)->Phi (), input.at (1)->Phi ())) ;
      lorentzVector vbf_sum = *input.at (0) + *input.at (1) ;
      vbf_mjj.Fill (sqrt (vbf_sum.M2 ())) ;
      vbf_detaMjj.Fill (sqrt (vbf_sum.M2 ()), fabs (input.at (0)->Eta () - input.at (1)->Eta ())) ;
      ntuple->Fill (input.at (seq[0])->Pt (), input.at (seq[1])->Pt (), 
                    sqrt (vbf_sum.M2 ()), fabs (input.at (0)->Eta () - input.at (1)->Eta ())) ;
      W_mj.Fill (sqrt (input.at (2)->M2 ())) ;

    }

  void write (TFile & output)
    {
      output.cd () ;
      ntuple->Write () ;
      vbf_deta.Write () ;
      vbf_dphi.Write () ;
      vbf_mjj.Write () ;
      vbf_detaMjj.Write () ;
      W_mj.Write () ;
      ptj1.Write () ;
      ptj2.Write () ;
      ptj3.Write () ;
      etaj1.Write () ;
      etaj2.Write () ;
      etaj3.Write () ;
      etavbfjout.Write () ;
      etavbfjin.Write () ;
      etatj.Write () ;
      etaWj.Write () ;
    }

  TString baseName ;
  TNtuple * ntuple ;
  TH1F vbf_deta ;
  TH1F vbf_dphi ;
  TH1F vbf_mjj ;
  TH2F vbf_detaMjj ;
  TH1F W_mj ;
  TH1F ptj1 ;
  TH1F ptj2 ;
  TH1F ptj3 ;
  TH1F etaj1 ;
  TH1F etaj2 ;
  TH1F etaj3 ;
  TH1F etavbfjout ;
  TH1F etavbfjin  ;
  TH1F etatj ;
  TH1F etaWj ;
} ;


#endif
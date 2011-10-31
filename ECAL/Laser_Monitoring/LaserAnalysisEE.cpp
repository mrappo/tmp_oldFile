#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TDirectory.h"
#include "TFile.h"



#include "./LaserDataAnalysis.h"
#include "./StepRecognition.h"

using namespace std;

int main(int argc, char ** argv)
{
  int selected_fed = atoi(argv[1]);
  
  if ( selected_fed >= 610 &&  selected_fed <= 645){
    cout << ">>>>>  Fed  " << selected_fed << " is not in EE !!!"  << endl;
    cout << "       Exiting ...  "  << endl;
    return 0;
  }
    
  bool saveAllChannels = false ;

  int * wl = NULL, nwl = 0;
  TChain * tx = new TChain("x");
  tx->Add(argv[2]);
  
  init_ttree(tx, &x);
  tx->SetBranchStatus("*",0); //disable all branches
  // only read needed branches
  tx->SetBranchStatus("run",1);
  tx->SetBranchStatus("fed",1);
  tx->SetBranchStatus("seq",1);
  tx->SetBranchStatus("detId",1);
  tx->SetBranchStatus("apdpnAB",1);
  tx->SetBranchStatus("apdpnA",1);
  tx->SetBranchStatus("apdpnB",1);
  tx->SetBranchStatus("time",1);
  tx->SetBranchStatus("eta",1);
  tx->SetBranchStatus("lumi",1);
  tx->SetBranchStatus("field", 1);
  tx->SetBranchStatus("wl",1);
  tx->SetBranchStatus("harness",1);
  tx->SetBranchStatus("elecId",1);
  tx->SetBranchStatus("tmax",1);
  tx->SetBranchStatus("l_ampli",1);
  tx->SetBranchStatus("l_rise_time",1);
  tx->SetBranchStatus("l_fwhm",1);
  tx->SetBranchStatus("l_prepulse",1);

  cout << " Number of entries " << tx->GetEntries() << endl;
  
  // 18 --> number of feds in EE
  // 850 : max number of channels in the fed
  TGraph* gvptpnlas[850];
  TGraph* gvptpnled[850];
  TGraph* gratio[850];
  char gname[100];
  int fed = 0;

    for (int ixtal = 0; ixtal < 850; ixtal++){
      sprintf(gname,"g_vptpnlas_fed%d_xtal%d", selected_fed, ixtal);
      gvptpnlas[ixtal] = new TGraph();
      gvptpnlas[ixtal] -> SetTitle(gname);
      gvptpnlas[ixtal] -> SetName(gname);
      gvptpnlas[ixtal] -> SetMarkerColor(kBlue);
      gvptpnlas[ixtal] -> SetMarkerStyle(20);
      gvptpnlas[ixtal] -> SetMarkerSize(0.5);
      
      sprintf(gname,"g_vptpnled_fed%d_xtal%d", selected_fed, ixtal);
      gvptpnled[ixtal] = new TGraph();
      gvptpnled[ixtal] -> SetTitle(gname);
      gvptpnled[ixtal] -> SetName(gname);
      gvptpnled[ixtal] -> SetMarkerColor(kCyan+2);
      gvptpnled[ixtal] -> SetMarkerStyle(20);
      gvptpnled[ixtal] -> SetMarkerSize(0.5);
      
      sprintf(gname,"g_ratio_fed%d_xtal%d", selected_fed, ixtal);
      gratio[ixtal] = new TGraph();
      gratio[ixtal] -> SetTitle(gname);
      gratio[ixtal] -> SetName(gname);
      gratio[ixtal] -> SetMarkerColor(kBlack);
      gratio[ixtal] -> SetMarkerStyle(20);
      gratio[ixtal] -> SetMarkerSize(0.5);
    }
  
  
  // Plot results:
  TTimeStamp dateMin(2011, 1,  01, 0, kTRUE, 0); 
  TTimeStamp dateMax(2011, 10, 31, 0, kTRUE, 0); 
  
  TTimeStamp tfake(2011, 9, 28, 0, kTRUE, 0);

  // profile averaging on one harness
  TProfile *p_vptpn_las[20];
  TProfile *p_vptpn_led[20];
  TProfile *p_ratio[20];
  TProfile *p_ratio_corr[20];
  TProfile *p_dT[20];
  
  // other quantities per harness 
  TProfile *p_ratiopn_las[20];
  TProfile *p_ratiopn_led[20];
  TProfile *p_tmax_las[20];
  TProfile *p_tmax_led[20];
  
  TProfile *p_matacqampli_las[20];
  TProfile *p_matacqrisetime_las[20];
  TProfile *p_matacqfwhm_las[20];
  TProfile *p_matacqprepulse_las[20];
  TProfile *p_matacqtmax_las[20];

  // TH2
  TH2F *h_ratio_vs_matacqfwhm[20];
  TH2F *h_ratio_vs_matacqampli[20];
  //TH3F *g2[20];
  //TGraph2D
  TGraph2D *g2[20];
   	  	   
    for (int ih = 0; ih < 20; ih++){
      
      sprintf(gname,"p_vptpn_las_fed%d_harness%d", selected_fed, ih+1);
      p_vptpn_las[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_vptpn_las[ih] ->SetLineColor(kBlue);
      p_vptpn_las[ih] ->SetMarkerColor(kBlue);
      p_vptpn_las[ih] ->SetMarkerStyle(20);
      p_vptpn_las[ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_vptpn_led_fed%d_harness%d", selected_fed, ih+1);
      p_vptpn_led[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_vptpn_led[ih] ->SetLineColor(kCyan+2);
      p_vptpn_led[ih] ->SetMarkerColor(kCyan+2);
      p_vptpn_led[ih] ->SetMarkerStyle(20);
      p_vptpn_led[ih] ->SetMarkerSize(0.5);

      sprintf(gname,"p_ratio_fed%d_harness%d", selected_fed, ih+1);
      p_ratio[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_ratio[ih] ->SetLineColor(kBlack);
      p_ratio[ih] ->SetMarkerColor(kBlack);
      p_ratio[ih] ->SetMarkerStyle(20);
      p_ratio[ih] ->SetMarkerSize(0.5);   
      
      sprintf(gname,"p_ratio_corr_fed%d_harness%d", selected_fed, ih+1);
      p_ratio_corr[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_ratio_corr[ih] ->SetLineColor(kBlack);
      p_ratio_corr[ih] ->SetMarkerColor(kBlack);
      p_ratio_corr[ih] ->SetMarkerStyle(20);
      p_ratio_corr[ih] ->SetMarkerSize(0.5);   
      
      sprintf(gname,"p_dT_fed%d_harness%d", fed, ih+1);
      p_dT[ih] = new TProfile(gname,gname,10000,dateMin, dateMax, -5000,5000);
      p_dT[ih] ->SetLineColor(kBlack);
      p_dT[ih] ->SetMarkerColor(kBlack);
      p_dT[ih] ->SetMarkerStyle(20);
      p_dT[ih] ->SetMarkerSize(0.5);


      sprintf(gname,"p_ratiopn_las_fed%d_harness%d", selected_fed, ih+1);
      p_ratiopn_las[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_ratiopn_las[ih] ->SetLineColor(kBlue);
      p_ratiopn_las[ih] ->SetMarkerColor(kBlue);
      p_ratiopn_las[ih] ->SetMarkerStyle(20);
      p_ratiopn_las[ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_ratiopn_led_fed%d_harness%d", selected_fed, ih+1);
      p_ratiopn_led[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_ratiopn_led[ih] ->SetLineColor(kCyan+2);
      p_ratiopn_led[ih] ->SetMarkerColor(kCyan+2);
      p_ratiopn_led[ih] ->SetMarkerStyle(20);
      p_ratiopn_led[ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_tmax_las_fed%d_harness%d", selected_fed, ih+1);
      p_tmax_las[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_tmax_las[ih] ->SetLineColor(kBlue);
      p_tmax_las[ih] ->SetMarkerColor(kBlue);
      p_tmax_las[ih] ->SetMarkerStyle(20);
      p_tmax_las[ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_tmax_led_fed%d_harness%d", selected_fed, ih+1);
      p_tmax_led[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_tmax_led[ih] ->SetLineColor(kCyan+2);
      p_tmax_led[ih] ->SetMarkerColor(kCyan+2);
      p_tmax_led[ih] ->SetMarkerStyle(20);
      p_tmax_led[ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqampli_las_fed%d_harness%d", selected_fed, ih+1);
      p_matacqampli_las[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqampli_las[ih] ->SetLineColor(kBlue);
      p_matacqampli_las[ih] ->SetMarkerColor(kBlue);
      p_matacqampli_las[ih] ->SetMarkerStyle(20);
      p_matacqampli_las[ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqrisetime_las_fed%d_harness%d", selected_fed, ih+1);
      p_matacqrisetime_las[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqrisetime_las[ih] ->SetLineColor(kBlue);
      p_matacqrisetime_las[ih] ->SetMarkerColor(kBlue);
      p_matacqrisetime_las[ih] ->SetMarkerStyle(20);
      p_matacqrisetime_las[ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqfwhm_las_fed%d_harness%d", selected_fed, ih+1);
      p_matacqfwhm_las[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqfwhm_las[ih] ->SetLineColor(kBlue);
      p_matacqfwhm_las[ih] ->SetMarkerColor(kBlue);
      p_matacqfwhm_las[ih] ->SetMarkerStyle(20);
      p_matacqfwhm_las[ih] ->SetMarkerSize(0.5);   

      sprintf(gname,"p_matacqprepulse_las_fed%d_harness%d", selected_fed, ih+1);
      p_matacqprepulse_las[ih] = new TProfile(gname,gname,8000,dateMin, dateMax, 0,100);
      p_matacqprepulse_las[ih] ->SetLineColor(kBlue);
      p_matacqprepulse_las[ih] ->SetMarkerColor(kBlue);
      p_matacqprepulse_las[ih] ->SetMarkerStyle(20);
      p_matacqprepulse_las[ih] ->SetMarkerSize(0.5);  
      
      sprintf(gname,"h_ratio_vs_matacqfwhm_fed%d_harness%d", fed, ih+1);
      h_ratio_vs_matacqfwhm[ih] = new TH2F(gname,gname, 100,0,50, 200, 0.9,1.1);
      h_ratio_vs_matacqfwhm[ih] ->GetXaxis()->SetTitle("fwhm (ns)");
      h_ratio_vs_matacqfwhm[ih] ->GetYaxis()->SetTitle("laser/LED");
 	  	 
      sprintf(gname,"h_ratio_vs_matacqampli_fed%d_harness%d", fed, ih+1); 
      h_ratio_vs_matacqampli[ih] = new TH2F(gname,gname,500,0,1000, 200, 0.9,1.1);
      h_ratio_vs_matacqampli[ih] ->GetXaxis()->SetTitle("amplitude (ADC)");
      h_ratio_vs_matacqampli[ih] ->GetYaxis()->SetTitle("laser/LED");
 	  	 
      sprintf(gname,"g2_fed%d_harness%d", fed, ih+1);
      g2[ih]= new TGraph2D();
      g2[ih]->SetTitle(gname);
      g2[ih]->SetName(gname);
      //g2[ih]= new TH3F(gname,gname,500,0.,50., 5000,0.,5000.,200,0.,2.);
 	  	 

    }
  
   TTimeStamp tmin1(2011, 1, 1, 0, kTRUE, 0);
   TTimeStamp tmax1(2011, 3, 20, 0, kTRUE, 0);
 	  	 
 	  	 
   TTimeStamp tmin2(2011, 6, 1, 0, kTRUE, 0);
   TTimeStamp tmax2(2011, 6, 6, 0, kTRUE, 0);
    
  int evtlas[850] = {0};
  int evtled[850] = {0};
  int evt[850] = {0};
  int tAve;

  float apdpnlas0[850], apdpnled0[850], tmaxlas0[850], tmaxled0[850] ;
  float las_ampli0[850], las_risetime0[850], las_fwhm0[850] , las_prepulse0[850];

  float apdpntemp[850] = {-999.};
  int tledtemp[850] = {0};
  float apdpnlednew;

  for (int ientry = 0; ientry < tx->GetEntries(); ientry++){
    tx->GetEntry(ientry);
    
    if (ientry%1000000 == 0 ) cout << "Analyzing entry " << ientry << endl;
    
    // select only EE
    if ( isEB(x.fed)) continue;
    
    int harn    = x.harness;
    int fed     = x.fed;
    int ism     = iSM(fed);
    float apdpnlas = x.apdpnAB[0];
    float apdpnled = x.apdpnAB[1];
    int tlas       = x.time[0];
    int tled       = x.time[1];
    int xtal       = x.elecId;      

       
    // check only one fed
    if ( fed != selected_fed ) continue;

    
    if (apdpnlas<=0) continue;
    if (apdpnled<=0) continue;
   
    // normalize to first point 
    if (evtlas[xtal] == 0) {
      apdpnlas0[xtal]     = apdpnlas;
      tmaxlas0[xtal]      = x.tmax[0];
      las_ampli0[xtal]    = x.l_ampli[0];
      las_risetime0[xtal] = x.l_rise_time[0];
      las_fwhm0[xtal]     = x.l_fwhm[0];
      las_prepulse0[xtal] = x.l_prepulse[0];
    }

    if (evtled[xtal] == 0) {
      apdpnled0[xtal]     = apdpnled;
      tmaxled0[xtal]      = x.tmax[1];
    }

    // average time between two measurements
    tAve = tlas + (tled - tlas)/ 2;
    if ( tlas > tled ) tAve = tled + (tlas - tled)/ 2;
    
    // extrapolate led measurement to laser time
         if ( apdpntemp[xtal] > 0 && tledtemp[xtal] > 0) 
           apdpnlednew = (apdpnled - apdpntemp[xtal] )/(tled - tledtemp[xtal]) * (tlas-tled) + apdpnled;
         else apdpnlednew = apdpnled;
    
    
    //apply a fake step
    //if ( tlas > int(tfake.GetSec())) apdpnlas+=0.005;

    float ratio = (apdpnlas/apdpnlas0[xtal])/(apdpnled/apdpnled0[xtal]);

    // --- profile plots  (average on one harness)
    p_vptpn_las[harn-1]->Fill(tlas, apdpnlas/apdpnlas0[xtal]);
    p_vptpn_led[harn-1]->Fill(tled, apdpnled/apdpnled0[xtal]);
    p_ratio[harn-1]   ->Fill(tAve, (apdpnlas/apdpnlas0[xtal])/(apdpnled/apdpnled0[xtal]));
    p_ratio_corr[harn-1] ->Fill(tlas, apdpnlas/apdpnlednew);
    p_dT[harn-1]       ->Fill(tlas, tled-tlas);
    

    if (x.apdpnB[0]>0 ) p_ratiopn_las[harn-1]->Fill(tlas, x.apdpnA[0]/x.apdpnB[0]);
    if (x.apdpnB[1]>0 ) p_ratiopn_led[harn-1]->Fill(tled, x.apdpnA[1]/x.apdpnB[1]);  
    
     //if ( (tlas > tmin1.GetSec() && tlas < tmax1.GetSec()) || (tlas > tmin2.GetSec() && tlas < tmax2.GetSec()) ){
     if ( tlas > tmin1.GetSec() && tlas < tmax1.GetSec() ){
     if ( fabs(ratio-1)<0.1 ){
 	  	               h_ratio_vs_matacqfwhm[harn-1] ->Fill( x.l_fwhm[0],ratio );
 	  	               h_ratio_vs_matacqampli[harn-1]->Fill( x.l_ampli[0],ratio );
 	  	               g2[harn-1]->SetPoint(evt[harn],x.l_fwhm[0],x.l_ampli[0], ratio);
	  	               //g2[harn-1]->Fill(x.l_fwhm[0],x.l_ampli[0],ratio);
 	  	               evt[harn]++;
 	  	 
 	  	             }
     }

    p_tmax_las[harn-1]->Fill(tlas, x.tmax[0]/tmaxlas0[xtal]);
    p_tmax_led[harn-1]->Fill(tled, x.tmax[1]/tmaxled0[xtal]);

    p_matacqampli_las[harn-1]->Fill(tlas, x.l_ampli[0]/las_ampli0[xtal]);

    p_matacqrisetime_las[harn-1]->Fill(tlas, x.l_rise_time[0]/las_risetime0[xtal]);

    p_matacqfwhm_las[harn-1]->Fill(tlas, x.l_fwhm[0]/las_fwhm0[xtal]);

    p_matacqprepulse_las[harn-1]->Fill(tlas, x.l_prepulse[0]/las_prepulse0[xtal]);

    // plot by channel
    if (saveAllChannels){
      gvptpnlas[xtal]->SetPoint(evtlas[xtal], tlas , apdpnlas);
      gvptpnled[xtal]->SetPoint(evtled[xtal], tled , apdpnled);
      gratio[xtal]->SetPoint(evtled[xtal], tAve , apdpnlas/apdpnled);
      //gratio[xtal]->SetPoint(evtled[xtal], tlas , apdpnlas/apdpnlednew);
    }
    
    evtlas[xtal]++;
    evtled[xtal]++;
    
    apdpntemp[xtal]= apdpnled;
    tledtemp[xtal] = tled;

    
  }// end loop over entries
   

     
  char fname[100];
  sprintf(fname,argv[3],selected_fed);

  TFile *fout = new TFile(fname,"recreate");
   
  if (saveAllChannels){
    TDirectory *mydir      =   fout -> mkdir("by_channel","by_channel");
    mydir->cd();
      for (int ixtal = 0; ixtal < 850; ixtal++){
	if ( gvptpnlas[ixtal]-> GetN() > 0) gvptpnlas[ixtal]-> Write();
	if ( gvptpnled[ixtal]-> GetN() > 0) gvptpnled[ixtal]-> Write();
	if ( gratio[ixtal]-> GetN()    > 0) gratio[ixtal]-> Write();
      }
  }
 
  fout->cd();

    for (int ih = 0; ih < 20; ih++){
      if ( p_vptpn_las[ih]-> GetEntries() > 0)  	p_vptpn_las[ih]->Write();
      if ( p_vptpn_led[ih]-> GetEntries() > 0)  	p_vptpn_led[ih]->Write();
      if ( p_ratio[ih]-> GetEntries()    > 0)  	p_ratio[ih]->Write();
      if ( p_ratio_corr[ih]-> GetEntries()    > 0)  	p_ratio_corr[ih]->Write();
      if ( p_dT[ih]       -> GetEntries() > 0)     p_dT[ih]->Write();
      if ( p_ratiopn_las[ih]-> GetEntries() > 0) p_ratiopn_las[ih]->Write();
      if ( p_ratiopn_led[ih]-> GetEntries() > 0) p_ratiopn_led[ih]->Write();

      if ( p_tmax_las[ih]-> GetEntries() > 0) p_tmax_las[ih]->Write();
      if ( p_tmax_led[ih]-> GetEntries() > 0) p_tmax_led[ih]->Write();
      
      if ( p_matacqampli_las[ih]-> GetEntries() > 0) p_matacqampli_las[ih]->Write();
      if ( p_matacqrisetime_las[ih]-> GetEntries() > 0) p_matacqrisetime_las[ih]->Write();
      if ( p_matacqfwhm_las[ih]-> GetEntries() > 0) p_matacqfwhm_las[ih]->Write();
      if ( p_matacqprepulse_las[ih]-> GetEntries() > 0) p_matacqprepulse_las[ih]->Write();
      if ( h_ratio_vs_matacqfwhm[ih] -> GetEntries() > 0 ) h_ratio_vs_matacqfwhm[ih]->Write();
      if ( h_ratio_vs_matacqampli[ih]-> GetEntries() > 0 ) h_ratio_vs_matacqampli[ih]->Write();
     if ( g2[ih]                    -> GetN() > 0       ) g2[ih]->Write();
      
    }  
  fout->Close();

   for(int i=0; i<20; i++)
   { 
    if(p_matacqampli_las[i]->GetEntries()==0) continue;

    TString output_dumper_steps=Form(argv[4],i+1); 

    Step_ID(p_matacqampli_las[i],output_dumper_steps);
  }
 
 
return 0;
  
}

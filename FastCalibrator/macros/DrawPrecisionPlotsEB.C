//
// Macro to produce ECAL single electron calibration plots
//


void DrawPrecisionPlotsEB(
			    Char_t* infile2 = "/data1/rgerosa/L3_Weight/MC_WJets/noEP_Z/Odd_WZAnalysis_SingleEle_WJetsToLNu_Z_noEP.root ",
			    Char_t* infile3 = "/data1/rgerosa/L3_Weight/MC_WJets/noEP_Z/Even_WZAnalysis_SingleEle_WJetsToLNu_Z_noEP.root",
			    int evalStat = 1,
                            int inputLoops = 25,
			    Char_t* fileType = "png", 
			    Char_t* dirName = ".")
{

  /// Draw plots in different modality : xtal, TT, SM
  TString modality = "xtal";
  bool printPlots = false;
  const int nLoops = inputLoops;

  // Set style options
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(1); 
  gStyle->SetOptStat(1110); 
  gStyle->SetOptFit(0); 
  gStyle->SetFitFormat("6.3g"); 
  gStyle->SetPalette(1); 
 
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.05);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.05);
  gStyle->SetTitleXOffset(0.8);
  gStyle->SetTitleYOffset(1.1);
  gROOT->ForceStyle();


  if ( evalStat && (!infile2 || !infile3 )){
    cout << " No input files to evaluate statistical precision specified !" << endl;
    return;
}

  cout << "Making calibration plots for: " << infile2 << endl;
  
  char hname[100];
  char htitle[100];
  TF1 *fgaus = new TF1("fgaus","gaus",-10,10);
  TF1 *pol0_0 = new TF1("pol0_0","pol1",0,20);
  TF1 *pol0_1 = new TF1("pol0_1","pol1",20,40);
  TF1 *pol0_2 = new TF1("pol0_2","pol1",40,60);
  TF1 *pol0_3 = new TF1("pol0_3","pol1",60,85);

  if (evalStat){
    
    TFile *f2 = new TFile(infile2);
    TFile *f3 = new TFile(infile3);
    TH2F *hcmap2[nLoops];
    TH2F *hcmap3[nLoops];
    TGraphErrors *statprecision_vs_ieta[nLoops];
    TGraphErrors *statprecision_vs_loop[4];

    int ipoint = 0;

    for ( int ietaregion=0; ietaregion<4; ietaregion++){
      statprecision_vs_loop[ietaregion] = new TGraphErrors();
      statprecision_vs_loop[ietaregion]->SetMarkerStyle(20);
      statprecision_vs_loop[ietaregion]->SetMarkerSize(1);
      statprecision_vs_loop[ietaregion]->SetMarkerColor(kBlue+2);
      if (ietaregion == 0) statprecision_vs_loop[ietaregion]->SetTitle("i#eta < 20");
      if (ietaregion == 1) statprecision_vs_loop[ietaregion]->SetTitle("20 < i#eta < 40");
      if (ietaregion == 2) statprecision_vs_loop[ietaregion]->SetTitle("40 < i#eta < 60");
      if (ietaregion == 3) statprecision_vs_loop[ietaregion]->SetTitle("60 < i#eta < 85");
    }

        
    TH2F *hcmap2_TT[nLoops];
    TH2F *hcmap3_TT[nLoops];

    for ( int iLoop = 0; iLoop < nLoops; iLoop++ ) {
    
      sprintf(hname,"h2_%d_hC_scale_EB",iLoop);
      hcmap2[iLoop] = (TH2F*)f2->Get(hname); 
      hcmap3[iLoop] = (TH2F*)f3->Get(hname); 
  
      if ( modality == "TT" ) {
        sprintf(hname,"h2_%d_hC_scale_EB_TT",iLoop);
        hcmap2_TT[iLoop] = new TH2F(hname, hname, 72,1, 72, 35, -17, 18);
        sprintf(hname,"h3_%d_hC_scale_EB_TT",iLoop);
        hcmap3_TT[iLoop] = new TH2F(hname, hname, 72,1, 72, 35, -17, 18);
        
        for ( int iTTphi = 1; iTTphi < 72; iTTphi++ ) {
          
          for ( int iTTeta = -17; iTTeta < 18; iTTeta++ ) {
            if ( iTTeta == 0 ) continue;
            float theICsum = hcmap2[iLoop] -> Integral(iTTphi*5, iTTphi*5+4, (iTTeta+18)*5, (iTTeta+18)*5+4)/25.;
            hcmap2_TT[iLoop] -> SetBinContent (iTTphi, iTTeta+18, theICsum);
            theICsum = hcmap3[iLoop] -> Integral(iTTphi*5, iTTphi*5+4, (iTTeta+18)*5, (iTTeta+18)*5+4)/25.;
            hcmap3_TT[iLoop] -> SetBinContent (iTTphi, iTTeta+18, theICsum);

          }
          
        }
        
      }
            
      TH1F *hstatprecision[171];
      
      if ( modality == "xtal" ) {
      
        for (int jbin = 1; jbin < hcmap2[iLoop]-> GetNbinsY()+1; jbin++){
          //int etaring = -85+(jbin-1);
          float etaring = hcmap2[iLoop]-> GetYaxis()->GetBinCenter(jbin);
          sprintf(hname,"hstatprecision_ring_ieta%02d",etaring);
          //hstatprecision[jbin-1] = new TH1F(hname, hname, 250,-0.1,0.1);
          hstatprecision[jbin-1] = new TH1F(hname, hname, 500,-0.5,0.5);
          for (int ibin = 1; ibin < hcmap2[iLoop]-> GetNbinsX()+1; ibin++){

            float ic1 = hcmap2[iLoop]->GetBinContent(ibin,jbin);
            float ic2 = hcmap3[iLoop]->GetBinContent(ibin,jbin);
            if (ic1>0 && ic1<2 && ic1!=1 && ic2>0 && ic2 <2 && ic2!=1){
              hstatprecision[jbin-1]->Fill((ic1-ic2)/(ic1+ic2)); // sigma (diff/sum) gives the stat. precision on teh entire sample

            }
          }
        }
      
      }
      
      else if ( modality == "TT" ) {
        
        for (int jbin = 1; jbin < hcmap2_TT[iLoop]-> GetNbinsY()+1; jbin++){
          //int etaring = -85+(jbin-1);
          float etaring = hcmap2_TT[iLoop]-> GetYaxis()->GetBinCenter(jbin);
          sprintf(hname,"hstatprecision_ring_ieta%02d",etaring);
          hstatprecision[jbin-1] = new TH1F(hname, hname, 250,-0.1,0.1);
          for (int ibin = 1; ibin < hcmap2_TT[iLoop]-> GetNbinsX()+1; ibin++){
            float ic1 = hcmap2_TT[iLoop]->GetBinContent(ibin,jbin);
            float ic2 = hcmap3_TT[iLoop]->GetBinContent(ibin,jbin);
            if (ic1>0 && ic1<2 && ic1!=1 && ic2>0 && ic2 <2 && ic2!=1){
              hstatprecision[jbin-1]->Fill((ic1-ic2)/(ic1+ic2)); // sigma (diff/sum) gives the stat. precision on teh entire sample
            }
          }
        }
        
      }
      else if ( modality == "SM" ) {
        
        for (int ibin = 1; ibin < hcmap2[iLoop]-> GetNbinsX()+1; ibin++){
          // Get the SM number
          int iSM = (ibin-1)/20 + 1;
          if ( (ibin-1)%20 == 0 ) {
            std::cout << iSM << std::endl;
            sprintf(hname,"hstatprecision_ring_ieta%02d",iSM);
            hstatprecision[iSM-1] = new TH1F(hname, hname, 250,-0.1,0.1);
            sprintf(hname,"hstatprecision_ring_ieta%02d",iSM+18);
            hstatprecision[iSM-1+18] = new TH1F(hname, hname, 250,-0.1,0.1);
          }
          for (int jbin = 1; jbin < hcmap2[iLoop]-> GetNbinsY()+1; jbin++){
            float ic1 = hcmap2[iLoop]->GetBinContent(ibin,jbin);
            float ic2 = hcmap3[iLoop]->GetBinContent(ibin,jbin);
            if (ic1>0 && ic1<2 && ic1!=1 && ic2>0 && ic2 <2 && ic2!=1){
              if ( jbin < 86 ) //EE-
               hstatprecision[iSM-1+18]->Fill((ic1-ic2)/(ic1+ic2)); // sigma (diff/sum) gives the stat. precision on teh entire sample
              else //EE+
               hstatprecision[iSM-1]->Fill((ic1-ic2)/(ic1+ic2)); 
            }
          }
        }
        
      }
  
      
      statprecision_vs_ieta[iLoop] = new TGraphErrors();
      statprecision_vs_ieta[iLoop]->SetMarkerStyle(20);
      statprecision_vs_ieta[iLoop]->SetMarkerSize(1);
      statprecision_vs_ieta[iLoop]->SetMarkerColor(kRed+2);
    
      int n = 0;
      
      if ( modality == "xtal" ) {

        for (int i = 1; i < hcmap2[iLoop]-> GetNbinsY()+1; i++){
          etaring = hcmap2[iLoop]-> GetYaxis()->GetBinCenter(i);
          if (etaring==0) continue;
          if ( hstatprecision[i-1]->GetEntries() == 0) continue;
          float e     = 0.5*hcmap2[iLoop]-> GetYaxis()->GetBinWidth(i);
          fgaus->SetParameter(1,1);
          fgaus->SetParameter(2,hstatprecision[i-1]->GetRMS());
          fgaus->SetRange(-5*hstatprecision[i-1]->GetRMS(),5*hstatprecision[i-1]->GetRMS());
          hstatprecision[i-1]->Fit("fgaus","QR");
          statprecision_vs_ieta[iLoop]-> SetPoint(n,etaring,fgaus->GetParameter(2));
          statprecision_vs_ieta[iLoop]-> SetPointError(n,e,fgaus->GetParError(2));
          n++;
        }
       
       statprecision_vs_ieta[iLoop]->Fit("pol0_0","QR");
       statprecision_vs_ieta[iLoop]->Fit("pol0_1","QR");
       statprecision_vs_ieta[iLoop]->Fit("pol0_2","QR");
       statprecision_vs_ieta[iLoop]->Fit("pol0_3","QR");

        statprecision_vs_loop[0]->SetPoint(ipoint,iLoop+1,pol0_0->GetParameter(0)+pol0_0->GetParameter(1)*10);
       statprecision_vs_loop[0]->SetPointError(ipoint,0.5,sqrt(pol0_0->GetParError(0)*pol0_0->GetParError(0)+pol0_0->GetParError(1)*pol0_0->GetParError(1)));
       statprecision_vs_loop[1]->SetPoint(ipoint,iLoop+1,pol0_1->GetParameter(0)+pol0_1->GetParameter(1)*30);
       statprecision_vs_loop[1]->SetPointError(ipoint,0.5,sqrt(pol0_1->GetParError(0)*pol0_1->GetParError(0)+pol0_1->GetParError(1)*pol0_1->GetParError(1)));
       statprecision_vs_loop[2]->SetPoint(ipoint,iLoop+1,pol0_2->GetParameter(0)+pol0_2->GetParameter(1)*50);
       statprecision_vs_loop[2]->SetPointError(ipoint,0.5,sqrt(pol0_0->GetParError(0)*pol0_2->GetParError(0)+pol0_2->GetParError(1)*pol0_2->GetParError(1)));
       statprecision_vs_loop[3]->SetPoint(ipoint,iLoop+1,pol0_3->GetParameter(0)+pol0_3->GetParameter(1)*72.5);
       statprecision_vs_loop[3]->SetPointError(ipoint,0.5,sqrt(pol0_3->GetParError(0)*pol0_3->GetParError(0)+pol0_3->GetParError(1)*pol0_3->GetParError(1)));
      ipoint++;

      }
      
      else if ( modality == "TT" ) {
        
        for (int i = 1; i < hcmap2_TT[iLoop]-> GetNbinsY()+1; i++){
          etaring = hcmap2_TT[iLoop]-> GetYaxis()->GetBinCenter(i);
          if (etaring==0) continue;
          if ( hstatprecision[i-1]->GetEntries() == 0) continue;
          float e     = 0.5*hcmap2_TT[iLoop]-> GetYaxis()->GetBinWidth(i);
          fgaus->SetParameter(1,1);
          fgaus->SetParameter(2,hstatprecision[i-1]->GetRMS());
          fgaus->SetRange(-5*hstatprecision[i-1]->GetRMS(),5*hstatprecision[i-1]->GetRMS());
          hstatprecision[i-1]->Fit("fgaus","QR");
          statprecision_vs_ieta[iLoop]-> SetPoint(n,etaring,fgaus->GetParameter(2));
          statprecision_vs_ieta[iLoop]-> SetPointError(n,e,fgaus->GetParError(2));
          n++;
        }
        
      }
      
      else if ( modality == "SM" ) {
        
        for (int i = 1; i < 36+1; i++){

          if ( hstatprecision[i-1]->GetEntries() == 0) continue;
          float e     = 0.5;
          fgaus->SetParameter(1,1);
          fgaus->SetParameter(2,hstatprecision[i-1]->GetRMS());
          fgaus->SetRange(-5*hstatprecision[i-1]->GetRMS(),5*hstatprecision[i-1]->GetRMS());
          hstatprecision[i-1]->Fit("fgaus","QR");
          statprecision_vs_ieta[iLoop]-> SetPoint(n,i,fgaus->GetParameter(2));
          statprecision_vs_ieta[iLoop]-> SetPointError(n,e,fgaus->GetParError(2));
          n++;
  
          
        }
        
        
      }
      
//       delete hstatprecision;
      
  }
}



  //------------------------------------------------------------------------


    
  
  //-----------------------------------------------------------------
  //--- Draw plots
  //-----------------------------------------------------------------
  TCanvas *c[nLoops];
  TCanvas *c2[4];
  //TFile * out = new TFile ("StatPrec_MC_noEP.root","RECREATE");
  TFile * out = new TFile ("StatPrec.root","RECREATE");
  // --- plot 5 : statistical precision vs ieta
  if (evalStat){
    
    for ( int iLoop = 0; iLoop < nLoops; iLoop++ ) {
      sprintf(hname,"cstat_%d",iLoop);
      c[iLoop] = new TCanvas(hname,hname);
      c[iLoop]->SetGridx();
      c[iLoop]->SetGridy();
      statprecision_vs_ieta[iLoop]->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.07);
      statprecision_vs_ieta[iLoop]->GetHistogram()->GetXaxis()-> SetRangeUser(-85,85);
      statprecision_vs_ieta[iLoop]->GetHistogram()->GetYaxis()-> SetTitle("#sigma((c_{P}-c_{D})/(c_{P}+c_{D}))");
      statprecision_vs_ieta[iLoop]->GetHistogram()->GetXaxis()-> SetTitle("i#eta");
      if ( modality == "SM" ) statprecision_vs_ieta[iLoop]->GetHistogram()->GetXaxis()-> SetTitle("iSM");
      statprecision_vs_ieta[iLoop]->Draw("ap");
      if(iLoop == nLoops -1)
      {
       out->cd();
       statprecision_vs_ieta[iLoop]->SetName("gr_stat_prec");
       statprecision_vs_ieta[iLoop]->Write();
      }
       
    }


  // --- plot 6 : statistical precision vs loop
  if (evalStat){
    
    for ( int ietaregion = 0; ietaregion < 4; ietaregion++ ) {
      sprintf(hname,"ietaregion_%d",ietaregion);
      c2[ietaregion] = new TCanvas(hname,hname);
      c2[ietaregion]->SetGridx();
      c2[ietaregion]->SetGridy();
      statprecision_vs_loop[ietaregion]->GetHistogram()->GetYaxis()-> SetRangeUser(0.,0.04);
      statprecision_vs_loop[ietaregion]->GetHistogram()->GetXaxis()-> SetRangeUser(0,nLoops+1);
      statprecision_vs_loop[ietaregion]->GetHistogram()->GetYaxis()-> SetTitle("Statistical precision");
      statprecision_vs_loop[ietaregion]->GetHistogram()->GetXaxis()-> SetTitle("n#circ iteration");
      statprecision_vs_loop[ietaregion]->Draw("ap");
    }
  }
}

  //-----------------------------------------------------------------
  //--- Print plots
  //-----------------------------------------------------------------
  
  if (printPlots){

    //gStyle->SetOptStat(1110);
    c[0]->Print("IC_map.png",fileType);
    c[1]->Print("IC_precision_vs_ieta.png",fileType);
    c[2]->Print("IC_scale_vs_ieta.png",fileType);
    c[3]->Print("IC_spread.png",fileType);
    c[4]->Print("occupancy_map.png",fileType);
}

}

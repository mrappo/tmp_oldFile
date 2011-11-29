//
// Macro to produce ECAL single electron calibration plots
//

void DrawFoldedPlotsEB(     Char_t* infile1 = "/data1/rgerosa/L3_Weight/PromptSkim/WZAnalysis_SingleEle_Run2011AB-PromptSkim_Z_noEP.root",
			    Char_t* infile2 = "/data1/rgerosa/L3_Weight/PromptSkim/Even_WZAnalysis_SingleEle_Run2011AB-PromptSkim_Z_noEP.root",
			    Char_t* infile3 = "/data1/rgerosa/L3_Weight/PromptSkim/Odd_WZAnalysis_SingleEle_Run2011AB-PromptSkim_Z_noEP.root",
			    int evalStat = 1,
			    Char_t* fileType = "png", 
			    Char_t* dirName = ".")
{

  bool  printPlots = false;

  // by xtal
  int nbins = 500;

  // Set style options
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
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

  if ( !infile1 ) {
    cout << " No input file specified !" << endl;
    return;
  }


  if ( evalStat && (!infile2 || !infile3 )){
    cout << " No input files to evaluate statistical precision specified !" << endl;
    return;
  }

  cout << "Making calibration plots for: " << infile1 << endl;
  
  TFile *f = new TFile(infile1);
  TH2F *hcmap = (TH2F*)f->Get("h_scale_map");

  //-----------------------------------------------------------------
  //--- Build the precision vs ieta plot starting from the TH2F of IC
  //-----------------------------------------------------------------
  
  TH1F *hspreadEtaFold[85];

  char hname[100];
  char htitle[100];
  
  int ringGroupSize = 1;
  int nStep = 0;

  for (int jbin = 1; jbin < hcmap-> GetNbinsY()+1; jbin++){
   if (jbin < 86 && (jbin-1)%ringGroupSize == 0 ) {
      nStep++;
      sprintf(hname,"hspread_ringGroup_ietaFolded%02d",nStep);
      hspreadEtaFold[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
   }
   if (jbin > 86 && (jbin-2)%ringGroupSize == 0 ) {
      nStep++;
   }

   for (int ibin = 1; ibin < hcmap-> GetNbinsX()+1; ibin++){
      float ic = hcmap->GetBinContent(ibin,jbin);
   if (ic>0 && ic<2 && ic!=1)    {
        if (nStep <= 85) hspreadEtaFold[nStep-1]->Fill(ic);
        else             hspreadEtaFold[170-nStep]->Fill(ic);
      }
    }
  }

  TGraphErrors *sigma_vs_EtaFold = new TGraphErrors();
  sigma_vs_EtaFold->SetMarkerStyle(20);
  sigma_vs_EtaFold->SetMarkerSize(1);
  sigma_vs_EtaFold->SetMarkerColor(kBlue+2);

  TF1 *fgaus = new TF1("fgaus","gaus",-10,10);
  TF1 *fgaus2 = new TF1("fgaus2","gaus",-100,100);

  int np = 0;

  for (int i = 1; i < 86; i++){
    float etaring = hcmap->GetYaxis()->GetBinCenter((ringGroupSize*i + ringGroupSize*(i-1))/2 + 1);
    float e     = 0.5*ringGroupSize;
    fgaus->SetParameter(1,1);
    fgaus->SetParameter(2,hspreadEtaFold[i-1]->GetRMS());
    fgaus->SetRange(1-5*hspreadEtaFold[i-1]->GetRMS(),1+5*hspreadEtaFold[i-1]->GetRMS());
    hspreadEtaFold[i-1]->Fit("fgaus","QR");

    sigma_vs_EtaFold-> SetPoint(np,fabs(etaring),fgaus->GetParameter(2));
    sigma_vs_EtaFold-> SetPointError(np,e,fgaus->GetParError(2));
    np++;
  }


  if (evalStat){
    TFile *f2 = new TFile(infile2);
    TFile *f3 = new TFile(infile3);

    TH2F *hcmap2 = (TH2F*)f2->Get("h_scale_map");
    TH2F *hcmap3 = (TH2F*)f3->Get("h_scale_map");

    TH1F *hstatprecisionEtaFold[85];

    nStep = 0;
    for (int jbin = 1; jbin < hcmap-> GetNbinsY()+1; jbin++){
    if (jbin < 86 && (jbin-1)%ringGroupSize == 0 ) {
       nStep++;
       sprintf(hname,"hstatprecision_ringGroup_ietaFolded%02d",nStep);
       hstatprecisionEtaFold[nStep-1]= new TH1F(hname, hname, nbins,-0.5,0.5);
    }
    if (jbin > 86 && (jbin-2)%ringGroupSize == 0 ) {
       nStep++;
    }

    for (int ibin = 1; ibin < hcmap2-> GetNbinsX()+1; ibin++){
       float ic1 = hcmap2->GetBinContent(ibin,jbin);
       float ic2 = hcmap3->GetBinContent(ibin,jbin);
    if (ic1>0 && ic1<2 && ic1!=1 && ic2>0 && ic2 <2 && ic2!=1)    {
        if (nStep <= 85) hstatprecisionEtaFold[nStep-1]->Fill((ic1-ic2)/(ic1+ic2));
        else             hstatprecisionEtaFold[170-nStep]->Fill((ic1-ic2)/(ic1+ic2));
       }
     }
   }

    
    TGraphErrors *statprecision_vs_EtaFold = new TGraphErrors();
    statprecision_vs_EtaFold->SetMarkerStyle(20);
    statprecision_vs_EtaFold->SetMarkerSize(1);
    statprecision_vs_EtaFold->SetMarkerColor(kRed+2);


    np = 0;

    for (int i = 1; i < 86; i++){
      float etaring = hcmap2->GetYaxis()->GetBinCenter((ringGroupSize*i + ringGroupSize*(i-1))/2 + 1);
      float e     = 0.5*ringGroupSize;
      fgaus->SetParameter(1,1);
      fgaus->SetParameter(2,hstatprecisionEtaFold[i-1]->GetRMS());
      fgaus->SetRange(-5*hstatprecisionEtaFold[i-1]->GetRMS(),5*hstatprecisionEtaFold[i-1]->GetRMS());
      hstatprecisionEtaFold[i-1]->Fit("fgaus","QR");

      statprecision_vs_EtaFold-> SetPoint(np,fabs(etaring),fgaus->GetParameter(2));
      statprecision_vs_EtaFold-> SetPointError(np,e,fgaus->GetParError(2));
      np++;
    }


    TGraphErrors *residual_vs_EtaFold = new TGraphErrors();
    residual_vs_EtaFold->SetMarkerStyle(20);
    residual_vs_EtaFold->SetMarkerSize(1);
    residual_vs_EtaFold->SetMarkerColor(kGreen+2);


    for (int i= 0; i < statprecision_vs_EtaFold-> GetN(); i++){
      double spread, espread;
      double stat, estat;
      double residual, eresidual;
      double xdummy,ex;
      sigma_vs_EtaFold-> GetPoint(i, xdummy, spread);
      espread = sigma_vs_EtaFold-> GetErrorY(i);
      statprecision_vs_EtaFold-> GetPoint(i, xdummy, stat);
      estat = statprecision_vs_EtaFold-> GetErrorY(i);
      ex = statprecision_vs_EtaFold-> GetErrorX(i);
      if (spread > stat ){
	residual  = sqrt( spread*spread - stat*stat );
	eresidual = sqrt( pow(spread*espread,2) + pow(stat*estat,2))/residual;
      }
      else {
	residual = 0;
	eresidual = 0;
      }
      residual_vs_EtaFold->SetPoint(i,xdummy, residual);
      residual_vs_EtaFold->SetPointError(i,ex,eresidual);
    }

  ///////////////////////////////////////////////////////////

    TGraphErrors *ic_vs_PhiFold = new TGraphErrors();
    ic_vs_PhiFold->SetMarkerStyle(20);
    ic_vs_PhiFold->SetMarkerSize(1);
    ic_vs_PhiFold->SetMarkerColor(kOrange+2);

    TH1F* hspreadPhiFold[18];
    int nStep =0;
    
    for(int jbin = 1; jbin < hcmap-> GetNbinsX()+1; jbin++){
    if ((jbin-1)%20 == 0 ) {
      nStep++;
      sprintf(hname,"hspread_iphiFolded%02d",nStep);
      hspreadPhiFold[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
    }
   
    for(int ibin = 1; ibin < hcmap-> GetNbinsY()+1; ibin++){
      float ic = hcmap->GetBinContent(jbin,ibin);
    if (ic>0 && ic<2 && ic!=1) {
      hspreadPhiFold[nStep-1]->Fill(ic);
      }
    }
  }
  
   np=0;
   for(int i=1; i<=18; i++)
  {
      float phibin = i;
      fgaus2->SetParameter(1,hspreadPhiFold[i-1]->GetMean());
      fgaus2->SetParameter(2,hspreadPhiFold[i-1]->GetRMS());
      fgaus2->SetRange(hspreadPhiFold[i-1]->GetMean()-5*hspreadPhiFold[i-1]->GetRMS(),hspreadPhiFold[i-1]->GetMean()+5*hspreadPhiFold[i-1]->GetRMS());
      hspreadPhiFold[i-1]->Fit("fgaus2","QR");
      ic_vs_PhiFold-> SetPoint(np,phibin,fgaus2->GetParameter(1));
      ic_vs_PhiFold-> SetPointError(np,0,fgaus2->GetParError(1));
      np++;
  

}

  ////////////////////////////////////////////////////////////////////

    TGraphErrors *ic_vs_PhiFold_crack_EBp = new TGraphErrors();
    ic_vs_PhiFold_crack_EBp->SetMarkerStyle(20);
    ic_vs_PhiFold_crack_EBp->SetMarkerSize(1);
    ic_vs_PhiFold_crack_EBp->SetMarkerColor(kRed);
    
    TGraphErrors *ic_vs_PhiFold_crack_EBm = new TGraphErrors();
    ic_vs_PhiFold_crack_EBm->SetMarkerStyle(20);
    ic_vs_PhiFold_crack_EBm->SetMarkerSize(1);
    ic_vs_PhiFold_crack_EBm->SetMarkerColor(kBlue);


    TH1F* hspreadPhiFold_crack_EBp[20];
    TH1F* hspreadPhiFold_crack_EBm[20];
    int nStep =0;
    
    for(int jbin = 1; jbin < hcmap-> GetNbinsX()+1; jbin++){
      if ((jbin-1)<20 ) {
        nStep++;
        sprintf(hname,"hspread_iphiFolded_crack_EBp%02d",nStep);
        hspreadPhiFold_crack_EBp[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);
        sprintf(hname,"hspread_iphiFolded_crack_EBm%02d",nStep);
        hspreadPhiFold_crack_EBm[nStep-1]= new TH1F(hname, hname, nbins/2,0.5,1.5);

      }
   
    for(int ibin = 1; ibin < hcmap-> GetNbinsY()+1; ibin++){
     float ic = hcmap->GetBinContent(jbin,ibin);
     if (ic>0 && ic<2 && ic!=1) {
       if((jbin-1)<20)
       {
          if ((ibin-85)<0 ) hspreadPhiFold_crack_EBm[jbin-1]->Fill(ic);
          if ((ibin-85)>=0 ) hspreadPhiFold_crack_EBp[jbin-1]->Fill(ic);
       }
       else{
           int kbin = jbin%20 ;
           if ((ibin-85)<0 ) hspreadPhiFold_crack_EBm[kbin]->Fill(ic);
           if ((ibin-85)>=0 ) hspreadPhiFold_crack_EBp[kbin]->Fill(ic);
   
         }
      }
    }
  }

   np=0;
   for(int i=1; i<=20; i++)
  {
      float phibin = i;
      fgaus2->SetParameter(1,hspreadPhiFold_crack_EBp[i-1]->GetMean());
      fgaus2->SetParameter(2,hspreadPhiFold_crack_EBp[i-1]->GetRMS());
      fgaus2->SetRange(hspreadPhiFold_crack_EBp[i-1]->GetMean()-5*hspreadPhiFold_crack_EBp[i-1]->GetRMS(),hspreadPhiFold_crack_EBp[i-1]->GetMean()+5*hspreadPhiFold_crack_EBp[i-1]->GetRMS());
      hspreadPhiFold_crack_EBp[i-1]->Fit("fgaus2","QR");
      ic_vs_PhiFold_crack_EBp-> SetPoint(np,phibin-1,fgaus2->GetParameter(1));
//       ic_vs_PhiFold_crack-> SetPoint(np,phibin,hspreadPhiFold_crack[i-1]->GetMean()); 
      ic_vs_PhiFold_crack_EBp-> SetPointError(np,0,fgaus2->GetParError(1));
      
      fgaus2->SetParameter(1,hspreadPhiFold_crack_EBm[i-1]->GetMean());
      fgaus2->SetParameter(2,hspreadPhiFold_crack_EBm[i-1]->GetRMS());
      fgaus2->SetRange(hspreadPhiFold_crack_EBm[i-1]->GetMean()-5*hspreadPhiFold_crack_EBm[i-1]->GetRMS(),hspreadPhiFold_crack_EBm[i-1]->GetMean()+5*hspreadPhiFold_crack_EBm[i-1]->GetRMS());
      hspreadPhiFold_crack_EBm[i-1]->Fit("fgaus2","QR");
      ic_vs_PhiFold_crack_EBm-> SetPoint(np,phibin-1,fgaus2->GetParameter(1));
//       ic_vs_PhiFold_crack-> SetPoint(np,phibin,hspreadPhiFold_crack[i-1]->GetMean()); 
      ic_vs_PhiFold_crack_EBm-> SetPointError(np,0,fgaus2->GetParError(1));
      
      np++;
  

}



}
  //------------------------------------------------------------------------

  //-----------------------------------------------------------------
  //--- Draw plots
  //-----------------------------------------------------------------
  TCanvas *c[10];

  c[0] = new TCanvas("csigmaFold","csigmaFold");
  c[0]->SetGridx();
  c[0]->SetGridy();
  sigma_vs_EtaFold->GetHistogram()->GetYaxis()-> SetRangeUser(0.00,0.10);
  sigma_vs_EtaFold->GetHistogram()->GetXaxis()-> SetRangeUser(0,85);
  sigma_vs_EtaFold->GetHistogram()->GetYaxis()-> SetTitle("#sigma_{c}");
  sigma_vs_EtaFold->GetHistogram()->GetXaxis()-> SetTitle("|i#eta|");
  sigma_vs_EtaFold->Draw("ap");
  if (evalStat){
    statprecision_vs_EtaFold->Draw("psame");
    sigma_vs_EtaFold->Draw("psame");
    TLegend * leg = new TLegend(0.6,0.7,0.89, 0.89);
    leg->SetFillColor(0);
    leg->AddEntry(statprecision_vs_EtaFold,"statistical precision", "LP");
    leg->AddEntry(sigma_vs_EtaFold,"spread", "LP");
    leg->Draw("same");
  }


  c[1] = new TCanvas("cresidualFold","cresidualFold");
  c[1]->SetGridx();
  c[1]->SetGridy();
  residual_vs_EtaFold->GetHistogram()->GetYaxis()-> SetRangeUser(0.0001,0.05);
  residual_vs_EtaFold->GetHistogram()->GetXaxis()-> SetRangeUser(0,85);
  residual_vs_EtaFold->GetHistogram()->GetYaxis()-> SetTitle("residual spread");
  residual_vs_EtaFold->GetHistogram()->GetXaxis()-> SetTitle("|i#eta|");
  residual_vs_EtaFold->Draw("ap");



   c[2] = new TCanvas("cphimeanfold","cphimeanfold");
   c[2]->SetGridx();
   c[2]->SetGridy();
   ic_vs_PhiFold->GetHistogram()->GetYaxis()-> SetRangeUser(0.99.,1.01);
   ic_vs_PhiFold->GetHistogram()->GetXaxis()-> SetRangeUser(0,18);
   ic_vs_PhiFold->GetHistogram()->GetYaxis()-> SetTitle("mean IC");
   ic_vs_PhiFold->GetHistogram()->GetXaxis()-> SetTitle("Phi%20");
   ic_vs_PhiFold->Draw("ap");

   c[3] = new TCanvas("cphimeanfold_crack_EB+","cphimeanfold_crack_EB+");
   c[3]->SetGridx();
   c[3]->SetGridy();
   ic_vs_PhiFold_crack_EBp->SetTitle(" Mean IC EB+");
   ic_vs_PhiFold_crack_EBp->GetHistogram()->GetYaxis()-> SetRangeUser(0.99.,1.01);
   ic_vs_PhiFold_crack_EBp->GetHistogram()->GetXaxis()-> SetRangeUser(0,20);
   ic_vs_PhiFold_crack_EBp->GetHistogram()->GetYaxis()-> SetTitle("mean IC");
   ic_vs_PhiFold_crack_EBp->GetHistogram()->GetXaxis()-> SetTitle("Phi%20");
   ic_vs_PhiFold_crack_EBp->Draw("ap");

   c[4] = new TCanvas("cphimeanfold_crack_EB-","cphimeanfold_crackEB-");
   c[4]->SetGridx();
   c[4]->SetGridy();
   ic_vs_PhiFold_crack_EBm->SetTitle(" Mean IC EB-");
   ic_vs_PhiFold_crack_EBm->GetHistogram()->GetYaxis()-> SetRangeUser(0.99.,1.01);
   ic_vs_PhiFold_crack_EBm->GetHistogram()->GetXaxis()-> SetRangeUser(0,20);
   ic_vs_PhiFold_crack_EBm->GetHistogram()->GetYaxis()-> SetTitle("mean IC");
   ic_vs_PhiFold_crack_EBm->GetHistogram()->GetXaxis()-> SetTitle("Phi%20");
   ic_vs_PhiFold_crack_EBm->Draw("ap");

}

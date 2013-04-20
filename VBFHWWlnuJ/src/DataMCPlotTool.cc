#include "DataMCPlotTool.h"

void DrawStackError(THStack* hs, const std::string & Labels,  const std::map<int,double> & SystematicErrorMap,  const bool & isLog, const bool & isLabel, const double & syst){

  TObjArray* histos = hs->GetStack () ;
  if (histos) {

    Int_t number = histos->GetEntries();
    TH1F* last = (TH1F*) histos->At (number-1) ;

    last->GetYaxis()->SetTitle("entries");
    last->GetYaxis()->SetTitleSize(0.045);
    last->GetYaxis()->SetLabelSize(0.035);
    last->GetYaxis()->SetTitleOffset(0.95);

    if(isLog){ last->SetMaximum(last->GetBinContent(last->GetMaximumBin())*10000);
               last->SetMinimum(0.001) ;
    }
    else { last->SetMaximum(last->GetBinContent(last->GetMaximumBin())*1.6);
           last->SetMinimum(0) ;
    }

    if(isLabel){ last->GetXaxis()->SetTitle(Labels.c_str());
                 last->GetXaxis()->SetTitleSize(0.04);
                 last->GetXaxis()->SetTitleOffset(0.95);
                 last->GetXaxis()->SetLabelSize(0.03);

                 last->GetYaxis()->SetTitleSize(0.035);
                 last->GetYaxis()->SetTitleOffset(1.25);
                 last->GetYaxis()->SetLabelSize(0.03);
    }
    else {
                 last->GetXaxis()->SetTitleSize(0.);
                 last->GetXaxis()->SetLabelSize(0.);
    }

    last->SetFillStyle(3001);
    last->SetLineColor(kBlack);
    last->DrawClone ("hist") ;

    std::vector <double> vErrSys (last->GetNbinsX(),0.) ;
    std::vector <double> vErrStat (last->GetNbinsX(),0.) ;

    for (int i = number-2 ; i >=0 ; --i) {
      TH1F * histo = (TH1F*) histos->At (i) ;
      histo->GetXaxis()->SetTitle(Labels.c_str());
      histo->SetFillStyle(3001);
      histo->SetLineColor(kBlack);
      histo->Draw ("same hist") ;
      for(int iBin = 0 ; iBin < histo->GetNbinsX(); iBin++){
	vErrSys.at(iBin) = vErrSys.at(iBin) + histo->GetBinContent(iBin+1)*histo->GetBinContent(iBin+1)*SystematicErrorMap.at(i)*SystematicErrorMap.at(i);
        vErrStat.at(iBin) = vErrStat.at(iBin) + histo->GetBinError(iBin+1)*histo->GetBinError(iBin+1);  
      }
    }
    
    //last->SetFillStyle(3005);
    //  last->SetFillColor(kBlack);
    // last->SetMarkerSize(0);
    
    for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
      double additionalError = sqrt(vErrSys.at(iBin) + last->GetBinContent(iBin+1) * last->GetBinContent(iBin+1) * syst * syst );
      last->SetBinError(iBin+1,sqrt(additionalError*additionalError + vErrStat.at(iBin)));
    }
    
    //    last->DrawClone ("sameE2") ;

  }
}


void DrawStackError(THStack* hs, const std::string & Labels, const TH1F* dataHist, const std::map<int,double> & SystematicErrorMap, const bool & isLog, const bool & isLabel, const double & syst){

  TObjArray* histos = hs->GetStack () ;
  if (histos) {

    Int_t number = histos->GetEntries();
    TH1F* last = (TH1F*) histos->At (number-1) ;

    last->GetYaxis()->SetTitle("entries");
    last->GetYaxis()->SetTitleSize(0.045);
    last->GetYaxis()->SetLabelSize(0.035);
    last->GetYaxis()->SetTitleOffset(0.95);

    if(isLog){ last->SetMaximum(std::max(last->GetBinContent(last->GetMaximumBin())*10000,dataHist->GetBinContent(dataHist->GetMaximumBin())*10000));
               last->SetMinimum(0.001) ;
    }
    else { last->SetMaximum(std::max(last->GetBinContent(last->GetMaximumBin())*1.6,dataHist->GetBinContent(dataHist->GetMaximumBin())*1.6));
           last->SetMinimum(0.) ;
    }

    if(isLabel){ last->GetXaxis()->SetTitle(Labels.c_str());
                 last->GetXaxis()->SetTitleSize(0.04);
                 last->GetXaxis()->SetTitleOffset(0.95);
                 last->GetXaxis()->SetLabelSize(0.03);

                 last->GetYaxis()->SetTitleSize(0.035);
                 last->GetYaxis()->SetTitleOffset(1.25);
                 last->GetYaxis()->SetLabelSize(0.03);
    }
    else {
                 last->GetXaxis()->SetTitleSize(0.);
                 last->GetXaxis()->SetLabelSize(0.);

    }
   
    last->SetFillStyle(3001);
    last->SetLineColor(kBlack);
    last->DrawClone ("hist") ;

    std::vector <double> vErrSys (last->GetNbinsX(),0.) ;
    std::vector <double> vErrStat (last->GetNbinsX(),0.) ;

    for (int i = number-2 ; i >=0 ; --i) {
      TH1F * histo = (TH1F*) histos->At (i) ;
      histo->GetXaxis()->SetTitle(Labels.c_str());
      histo->SetLineColor(kBlack);
      histo->SetFillStyle(3001);
      histo->Draw ("same hist") ;
      for(int iBin = 0 ; iBin < histo->GetNbinsX(); iBin++){
	vErrSys.at(iBin) = vErrSys.at(iBin) + histo->GetBinContent(iBin+1)*histo->GetBinContent(iBin+1)*SystematicErrorMap.at(i)*SystematicErrorMap.at(i);
        vErrStat.at(iBin) = vErrStat.at(iBin) + histo->GetBinError(iBin+1)*histo->GetBinError(iBin+1);  
      }
    }
    

    //    last->SetFillStyle(3005);
    //    last->SetFillColor(kBlack);
    //    last->SetMarkerSize(0);
    
    for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
      double additionalError = sqrt(vErrSys.at(iBin) + last->GetBinContent(iBin+1) * last->GetBinContent(iBin+1) * syst * syst );
      last->SetBinError(iBin+1,sqrt(additionalError*additionalError + vErrStat.at(iBin)));
    }
    
    //    last->DrawClone ("sameE2") ;
  }
}

void LatexCMS (const double & lumi, const std::string & LeptonType, const bool & isLabel){

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(21); // align right                                                                                                                        

  if(isLabel)   latex.SetTextSize(0.032);
  else   latex.SetTextSize(0.04);

  if(LeptonType == "muon" || LeptonType == "mu") latex.DrawLatex(0.55,0.962,Form("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, mu+jets",(float)lumi/1000));
  else if(LeptonType == "electron" || LeptonType == "el") latex.DrawLatex(0.55,0.962,Form("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, el+jets",(float)lumi/1000));

}


void SetTotalSystematicVector( std::vector<double> & SysError, THStack* hs, const std::map<int,double> & SystematicErrorMap, const double & syst){

 TObjArray* histos = hs->GetStack () ;
  if (histos) {

    Int_t number = histos->GetEntries();
    TH1F* last = (TH1F*) histos->At (number-1) ;

    SysError.assign(last->GetNbinsX(),0.);
    for (int i = number-2 ; i >=0 ; --i) {
      TH1F * histo = (TH1F*) histos->At (i) ;
      for(int iBin = 0 ; iBin < histo->GetNbinsX(); iBin++)
	SysError.at(iBin) = SysError.at(iBin) + histo->GetBinContent(iBin+1)*histo->GetBinContent(iBin+1)*SystematicErrorMap.at(i)*SystematicErrorMap.at(i) + 
                 	    last->GetBinContent(iBin+1) * last->GetBinContent(iBin+1) * syst * syst;
    }
  }

}

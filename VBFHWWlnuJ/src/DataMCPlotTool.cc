#include "DataMCPlotTool.h"

void DrawStackError(THStack* hs,const std::string & Labels, const std::map<int,double> & SystematicErrorMap, const double & syst){ 
  TObjArray* histos = hs->GetStack () ;
  if (histos) {

    Int_t number = histos->GetEntries();
    TH1F* last = (TH1F*) histos->At (number-1) ;

    last->GetXaxis()->SetTitle(Labels.c_str());
    last->GetXaxis()->SetTitleSize(0.04);
    last->GetYaxis()->SetTitle("Entries");
    last->GetYaxis()->SetTitleSize(0.04);
    last->SetMaximum(last->GetBinContent(last->GetMaximumBin())*1.8) ;
    last->SetMinimum(1) ;


    last->DrawClone ("hist") ;

    std::vector <double> vErrSys (last->GetNbinsX(),0.) ;
    std::vector <double> vErrStat (last->GetNbinsX(),0.) ;

    for (int i = number-2 ; i >=0 ; --i) {
      TH1F * histo = (TH1F*) histos->At (i) ;
      histo->GetXaxis()->SetTitle(Labels.c_str());
      histo->Draw ("same hist") ;
      for(int iBin = 0 ; iBin < histo->GetNbinsX(); iBin++){
	vErrSys.at(iBin) = vErrSys.at(iBin) + histo->GetBinContent(iBin+1)*histo->GetBinContent(iBin+1)*SystematicErrorMap.at(i)*SystematicErrorMap.at(i);
        vErrStat.at(iBin) = vErrStat.at(iBin) + histo->GetBinError(iBin+1)*histo->GetBinError(iBin+1);  
      }
    }
    
    Style_t origStyleLast = last->GetFillStyle ();
    Color_t origColorLast = last->GetFillColor ();
    last->SetFillStyle(3005);
    last->SetFillColor(kBlack);
    last->SetMarkerSize(0);
    
    for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
      double additionalError = sqrt(vErrSys.at(iBin) + last->GetBinContent(iBin+1) * last->GetBinContent(iBin+1) * syst * syst );
      last->SetBinError(iBin+1,sqrt(additionalError*additionalError + vErrStat.at(iBin)));
    }
    
    last->DrawClone ("sameE2") ;

  }
}


void DrawStackError(THStack* hs, const std::string & Labels, const TH1F*  dataHist, const std::map<int,double> & SystematicErrorMap, const double & syst){ 
  TObjArray* histos = hs->GetStack () ;
  if (histos) {

    Int_t number = histos->GetEntries();
    TH1F* last = (TH1F*) histos->At (number-1) ;

    last->GetXaxis()->SetTitle(Labels.c_str());
    last->GetXaxis()->SetTitleSize(0.04);
    last->GetYaxis()->SetTitle("Entries");
    last->GetYaxis()->SetTitleSize(0.04);
    last->SetMaximum(std::max(last->GetBinContent(last->GetMaximumBin())*1.3,dataHist->GetBinContent(dataHist->GetMaximumBin())*1.3)) ;
    last->SetMinimum(1) ;
 
    last->DrawClone ("hist") ;

    std::vector <double> vErrSys (last->GetNbinsX(),0.) ;
    std::vector <double> vErrStat (last->GetNbinsX(),0.) ;

    for (int i = number-2 ; i >=0 ; --i) {
      TH1F * histo = (TH1F*) histos->At (i) ;
      histo->GetXaxis()->SetTitle(Labels.c_str());
      histo->Draw ("same hist") ;
      for(int iBin = 0 ; iBin < histo->GetNbinsX(); iBin++){
	vErrSys.at(iBin) = vErrSys.at(iBin) + histo->GetBinContent(iBin+1)*histo->GetBinContent(iBin+1)*SystematicErrorMap.at(i)*SystematicErrorMap.at(i);
        vErrStat.at(iBin) = vErrStat.at(iBin) + histo->GetBinError(iBin+1)*histo->GetBinError(iBin+1);  
      }
    }
    
    Style_t origStyleLast = last->GetFillStyle ();
    Color_t origColorLast = last->GetFillColor ();
    last->SetFillStyle(3005);
    last->SetFillColor(kBlack);
    last->SetMarkerSize(0);
    
    for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
      double additionalError = sqrt(vErrSys.at(iBin) + last->GetBinContent(iBin+1) * last->GetBinContent(iBin+1) * syst * syst );
      last->SetBinError(iBin+1,sqrt(additionalError*additionalError + vErrStat.at(iBin)));
    }
    
    last->DrawClone ("sameE2") ;

  }
}


void LatexCMS (const double & lumi, const std::string & LeptonType){

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);

  latex.SetTextAlign(31); // align right                                                                                                                        
  latex.DrawLatex(0.68,0.962,"#sqrt{s} = 8 TeV");
  latex.SetTextAlign(31); // align right                                                                                                                        
  latex.DrawLatex(0.5,0.962,Form("#int #font[12]{L} dt = %.1f fb^{-1}", (float)lumi/1000));
  latex.SetTextAlign(31); // align right                                                                                                                        
  if(LeptonType=="Muon" || LeptonType=="muon" || LeptonType == "mu" ) latex.DrawLatex(0.85,0.962,Form("mu+jets"));
  else  latex.DrawLatex(0.85,0.962,Form("el+jets"));

  latex.SetTextAlign(11); // align left                                                                                                                         
  //  latex.DrawLatex(0.15,0.93,"CMS,  #sqrt{s} = 7 TeV");//preliminary 2011");                                                                                 
  latex.DrawLatex(0.15,0.962,"CMS preliminary");

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


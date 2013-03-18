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
      std::cout<<" additional "<<additionalError<<" stat "<<sqrt(vErrStat.at(iBin))<<std::endl;
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
    last->SetMaximum(std::max(last->GetBinContent(last->GetMaximumBin())*1.1,dataHist->GetBinContent(dataHist->GetMaximumBin())*1.1)) ;
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
      std::cout<<" additional "<<additionalError<<" stat "<<sqrt(vErrStat.at(iBin))<<std::endl;
      last->SetBinError(iBin+1,sqrt(additionalError*additionalError + vErrStat.at(iBin)));
    }
    
    last->DrawClone ("sameE2") ;

  }
}


void LatexCMS (double lumi){

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);

  latex.SetTextAlign(31); // align right                                                                                                                        
  latex.DrawLatex(0.80,0.962,"#sqrt{s} = 8 TeV");
  latex.SetTextAlign(31); // align right                                                                                                                        
  latex.DrawLatex(0.66,0.962,Form("#int #font[12]{L} dt = %.1f pb^{-1}", (float)lumi));

  latex.SetTextAlign(11); // align left                                                                                                                         
  //  latex.DrawLatex(0.15,0.93,"CMS,  #sqrt{s} = 7 TeV");//preliminary 2011");                                                                                 
  latex.DrawLatex(0.15,0.962,"CMS preliminary");

}

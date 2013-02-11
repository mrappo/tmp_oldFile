#include "DataMCPlotTool.h"

void DrawStackError(THStack* hs, double syst){ 
  TObjArray* histos = hs->GetStack () ;
  if (histos) {
    Int_t number = histos->GetEntries();
    TH1F* last = (TH1F*) histos->At (number-1) ;
    last->DrawClone ("hist") ;
    for (int i = number-2 ; i >=0 ; --i) {
      TH1F * histo = (TH1F*) histos->At (i) ;
      //    Style_t origStyle = histo->GetFillStyle ();
      //    Color_t origColor = histo->GetFillColor ();
      //    TH1F* dummy = (TH1F*) histo->Clone () ;
      //    dummy->SetFillStyle (1001) ; 
      //    dummy->SetFillColor (10) ;        
      //    dummy->Draw ("same hist") ;
      histo->Draw ("same hist") ;
    }
    Style_t origStyleLast = last->GetFillStyle ();
    Color_t origColorLast = last->GetFillColor ();
    last->SetFillStyle(3005);
    last->SetFillColor(kBlack);
    last->SetMarkerSize(0);
  
    std::vector <double> vErr ;
    for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
      double additionalError = last->GetBinContent(iBin+1) * syst;
      vErr.push_back(last->GetBinError(iBin+1));
      last->SetBinError(iBin+1,sqrt(additionalError*additionalError + last->GetBinError(iBin+1) * last->GetBinError(iBin+1)) );
    }
    last->DrawClone ("sameE2") ;
    //---- restore hist ----
    last->SetFillStyle(origStyleLast);
    last->SetFillColor(origColorLast);
    for (int iBin = 0 ; iBin < last->GetNbinsX(); iBin++) {
      last->SetBinError(iBin+1, vErr.at(iBin));
    }
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

#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TAxis.h"

void plotVsNvertex(TString vbfFileName="", TString ggHFileName=""){
  
  gStyle->SetOptStat(00000);

  std::string treeName = "otree";
  
  std::string variable_reco  = "abs(vbf_maxpt_j1_eta-vbf_maxpt_j2_eta)";
  std::string variable_gen   = "abs(vbf_maxpt_j1_eta_gen-vbf_maxpt_j2_eta_gen)";
  std::string variable_quark = "abs(genTagQuark1eta-genTagQuark2eta)";

  std::string cut_reco = "issignal && v_pt > 200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBin >= 2 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679";

  std::string cut_reco_nvtx10 = "issignal && v_pt > 200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBin >= 2 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && nPV <=10";

  std::string cut_reco_nvtx20 = "issignal && v_pt > 200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBin >= 2 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && nPV >10 && nPV<=20";

  std::string cut_reco_nvtx30 = "issignal && v_pt > 200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBin >= 2 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && nPV >20";

  std::string cut_gen = "issignal && v_pt > 200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBinGen >= 2 && vbf_maxpt_j1_bDiscriminatorCSV_gen <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV_gen <=0.679";

  std::string cut_quark = "issignal && v_pt > 200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBinGen >= 2 && vbf_maxpt_j1_bDiscriminatorCSV_gen <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV_gen <=0.679";
  
  //Import File
  TFile *fvbf = new TFile(vbfFileName.Data(),"READ");
  TFile *fggH = new TFile(ggHFileName.Data(),"READ");
  
  //Import Tree
  TTree *tvbf = (TTree*) fvbf->Get(treeName.c_str());
  TTree *tggH = (TTree*) fggH->Get(treeName.c_str());
  
  TH1F* histo_reco_vbf  = new TH1F("histo_reco_vbf","",100,0,9);
  TH1F* histo_reco_vbf_nvtx10  = new TH1F("histo_reco_vbf_nvtx","",100,0,9);
  TH1F* histo_reco_vbf_nvtx20  = new TH1F("histo_reco_vbf_nvtx20","",100,0,9);
  TH1F* histo_reco_vbf_nvtx30  = new TH1F("histo_reco_vbf_nvtx30","",100,0,9);
  TH1F* histo_gen_vbf   = new TH1F("histo_gen_vbf","",100,0,9);
  TH1F* histo_quark_vbf = new TH1F("histo_quark_vbf","",100,0,9);
  TH1F* histo_reco_ggH  = new TH1F("histo_reco_ggH","",100,0,9);
  TH1F* histo_gen_ggH   = new TH1F("histo_gen_ggH","",100,0,9); 
  

  tvbf->Draw((variable_reco+">> histo_reco_vbf").c_str(),cut_reco.c_str(),"goff");
  tvbf->Draw((variable_reco+">> histo_reco_vbf_nvtx").c_str(),cut_reco_nvtx10.c_str(),"goff");
  tvbf->Draw((variable_reco+">> histo_reco_vbf_nvtx20").c_str(),cut_reco_nvtx20.c_str(),"goff");
  tvbf->Draw((variable_reco+">> histo_reco_vbf_nvtx30").c_str(),cut_reco_nvtx30.c_str(),"goff");
  tggH->Draw((variable_reco+">> histo_reco_ggH").c_str(),cut_reco.c_str(),"goff");

  tvbf->Draw((variable_gen+">> histo_gen_vbf").c_str(),cut_gen.c_str(),"goff");
  tggH->Draw((variable_gen+">> histo_gen_ggH").c_str(),cut_gen.c_str(),"goff");

  tvbf->Draw((variable_quark+">> histo_quark_vbf").c_str(),cut_quark.c_str(),"goff");
  
  histo_reco_vbf->Scale(1/histo_reco_vbf->Integral());
  histo_reco_vbf_nvtx10->Scale(1/histo_reco_vbf_nvtx10->Integral());
  histo_reco_vbf_nvtx20->Scale(1/histo_reco_vbf_nvtx20->Integral());
  histo_reco_vbf_nvtx30->Scale(1/histo_reco_vbf_nvtx30->Integral());
  histo_gen_vbf->Scale(1/histo_gen_vbf->Integral());
  histo_quark_vbf->Scale(1/histo_quark_vbf->Integral());

  histo_reco_ggH->Scale(1/histo_reco_ggH->Integral());
  histo_gen_ggH->Scale(1/histo_gen_ggH->Integral());
		     
  // Plots ggH and vbf for each
  TCanvas* c1 = new TCanvas("c1","",500,515);
  c1->cd();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  histo_reco_vbf->SetMaximum(histo_reco_vbf->GetMaximum()*100);
  histo_reco_vbf->SetLineColor(kBlue);
  histo_reco_ggH->SetLineColor(kRed);
  histo_reco_vbf->SetLineWidth(2);
  histo_reco_ggH->SetLineWidth(2);
  histo_reco_vbf->GetXaxis()->SetTitle("tag jet #Delta#eta_{jj}"); 
  histo_reco_vbf->GetYaxis()->SetTitle("Normalized Unit"); 
  histo_reco_vbf->Draw(""); 
  histo_reco_ggH->Draw("same"); 

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(histo_reco_vbf,"qqH mH=600 GeV","l");
  leg->AddEntry(histo_reco_ggH,"ggH mH=600 GeV","l");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.032);
  leg->Draw("same"); 

  // Plots ggH and vbf for each
  TCanvas* c2 = new TCanvas("c2","",500,515);
  c2->cd();
  c2->SetLogy();
  c2->SetGridx();
  c2->SetGridy();
  histo_gen_vbf->SetMaximum(histo_gen_vbf->GetMaximum()*100);
  histo_gen_vbf->SetLineColor(kBlue);
  histo_gen_ggH->SetLineColor(kRed);
  histo_gen_vbf->SetLineWidth(2);
  histo_gen_ggH->SetLineWidth(2);
  histo_gen_vbf->GetXaxis()->SetTitle("tag jet #Delta#eta_{jj}"); 
  histo_gen_vbf->GetYaxis()->SetTitle("Normalized Unit"); 
  histo_gen_vbf->Draw(""); 
  histo_gen_ggH->Draw("same"); 

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(histo_reco_vbf,"qqH mH=600 GeV","l");
  leg->AddEntry(histo_reco_ggH,"ggH mH=600 GeV","l");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.032);
  leg->Draw("same"); 
  
  // Plots ggH and vbf for each
  TCanvas* c3 = new TCanvas("c3","",500,515);
  c3->cd();
  c3->SetLogy();
  c3->SetGridx();
  c3->SetGridy();
  histo_quark_vbf->SetMaximum(histo_quark_vbf->GetMaximum()*100);
  histo_quark_vbf->SetLineColor(kBlue);
  histo_quark_vbf->SetLineWidth(2);
  histo_quark_vbf->GetXaxis()->SetTitle("tag jet #Delta#eta_{jj}"); 
  histo_quark_vbf->GetYaxis()->SetTitle("Normalized Unit"); 
  histo_quark_vbf->Draw(""); 

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(histo_quark_vbf,"qqH mH=600 GeV","l");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.032);
  leg->Draw("same"); 

  // All together vbf 
  TCanvas* c4 = new TCanvas("c4","",500,515);
  c4->cd();
  c4->SetLogy();
  c4->SetGridx();
  c4->SetGridy();
  histo_reco_vbf->SetMaximum(histo_reco_vbf->GetMaximum());
  histo_reco_vbf->SetLineColor(kBlue);
  histo_gen_vbf_2 = (TH1F*) histo_gen_vbf->Clone("histo_gen_vbf_2");
  histo_gen_vbf_2->SetLineColor(kRed);
  histo_gen_vbf_2->SetLineWidth(2);
  histo_quark_vbf->SetLineColor(kBlack);
  histo_quark_vbf->SetLineWidth(2);
  histo_reco_vbf->GetXaxis()->SetTitle("tag jet #Delta#eta_{jj}"); 
  histo_reco_vbf->GetYaxis()->SetTitle("Normalized Unit"); 
  histo_reco_vbf->Draw(""); 
  histo_gen_vbf_2->Draw("same"); 
  histo_quark_vbf->Draw("same"); 

  TLegend *leg2 = new TLegend(0.5,0.6,0.9,0.9);
  leg2->AddEntry(histo_reco_vbf,"reco qqH mH=600 GeV","l");
  leg2->AddEntry(histo_gen_vbf_2,"gen qqH mH=600 GeV","l");
  leg2->AddEntry(histo_quark_vbf,"quark qqH mH=600 GeV","l");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.032);
  leg2->Draw("same");

  // All together vbf 
  TCanvas* c5 = new TCanvas("c5","",500,515);
  c5->cd();
  c5->SetLogy();
  c5->SetGridx();
  c5->SetGridy();
  histo_reco_vbf_nvtx20->SetMaximum(histo_reco_vbf_nvtx20->GetMaximum()*100);
  histo_reco_vbf_nvtx20->SetMinimum(0.0001);
  histo_reco_vbf_nvtx10->SetLineColor(kBlue);
  histo_reco_vbf_nvtx10->SetLineWidth(2);
  histo_reco_vbf_nvtx20->SetLineColor(kRed);
  histo_reco_vbf_nvtx20->SetLineWidth(2);
  histo_reco_vbf_nvtx30->SetLineColor(kBlack);
  histo_reco_vbf_nvtx30->SetLineWidth(2);
  histo_reco_vbf_nvtx20->GetXaxis()->SetTitle("tag jet #Delta#eta_{jj}"); 
  histo_reco_vbf_nvtx20->GetYaxis()->SetTitle("Normalized Unit"); 
  histo_reco_vbf_nvtx20->Draw(""); 
  histo_reco_vbf_nvtx10->Draw("same"); 
  histo_reco_vbf_nvtx30->Draw("same"); 

  TLegend *leg2 = new TLegend(0.5,0.6,0.9,0.9);
  leg2->AddEntry(histo_reco_vbf_nvtx10,"qqH mH=600 GeV N_{vxt}#in[0,10]","l");
  leg2->AddEntry(histo_reco_vbf_nvtx20,"qqH mH=600 GeV N_{vxt}#in[10,20]","l");
  leg2->AddEntry(histo_reco_vbf_nvtx30,"qqH mH=600 GeV N_{vxt}#in[20,50]","l");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.032);
  leg2->Draw("same");
  

}

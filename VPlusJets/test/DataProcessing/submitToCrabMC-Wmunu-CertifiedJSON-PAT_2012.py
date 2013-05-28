import os,sys
import string, re
from time import gmtime, localtime, strftime

isQCD    = "isQCD = False"
isHEEPID = "isHEEPID = True"
isTransverseMassCut = "isTransverseMassCut = False"

##------ Please set ONLY one of the four flags to True -------
physMode   = "WmunuJets_"
ConfigFile = "../WmunuJetsAnalysisPAT_cfg.py"
dataset    = [#"/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/jdamgov-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/ilyao-SQWaT_PAT_53X_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_V2-402ec72c9517e76360adc5ca179e6efb/USER",
              #"/W1JetsToLNu_TuneZ2Star_8TeV-madgraph/ilyao-SQWaT_PAT_53X_W1JetsToLNu_TuneZ2Star_8TeV-madgraph-402ec72c9517e76360adc5ca179e6efb/USER",
              #"/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/ilyao-SQWaT_PAT_53X_W2JetsToLNu_TuneZ2Star_8TeV-madgraph-402ec72c9517e76360adc5ca179e6efb/USER",
              #"/W3JetsToLNu_TuneZ2Star_8TeV-madgraph/ilyao-SQWaT_PAT_53X_W3JetsToLNu_TuneZ2Star_8TeV-madgraph-402ec72c9517e76360adc5ca179e6efb/USER",
              #"/W4JetsToLNu_TuneZ2Star_8TeV-madgraph/ilyao-SQWaT_PAT_53X_W4JetsToLNu_TuneZ2Star_8TeV-madgraph-402ec72c9517e76360adc5ca179e6efb/USER",
              "/WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph/ntran-SQWaT_PAT_WJetsPt100MG_53x-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/dimatteo-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/WZ_TuneZ2star_8TeV_pythia6_tauola/goodell-SQWaT_PAT_43X_Fall12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/ZZ_TuneZ2star_8TeV_pythia6_tauola/ajkumar-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/T_t-channel_TuneZ2star_8TeV-powheg-tauola/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/goodell-SQWaT_PAT_43X_Fall12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/T_s-channel_TuneZ2star_8TeV-powheg-tauola/goodell-SQWaT_PAT_43X_Fall12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/goodell-SQWaT_PAT_53X_Fall12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/goodell-SQWaT_PAT_43X_Fall12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/goodell-SQWaT_PAT_53X_Fall12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/WW_TuneZ2star_8TeV_pythia6_tauola/ntran-SQWaT_PAT_WW_53x-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/LQ-ggh600_GEN_53X/ntran-SQWaT_PAT_ggH600_53x-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/LQ-ggh700_GEN_53X/shuai-SQWaT_PAT_53X_ggH700-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/GluGluToHToWWToLAndTauNuQQ_M-800_8TeV-powheg-pythia6/zixu-SQWaT_PAT_53X_ggH800-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/GluGluToHToWWToLAndTauNuQQ_M-900_8TeV-powheg-pythia6/shuai-SQWaT_PAT_53X_ggH900-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/GluGluToHToWWToLAndTauNuQQ_M-1000_8TeV-powheg-pythia6/zixu-SQWaT_PAT_53X_ggH1000-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/VBF_HToWWToLAndTauNuQQ_M-600_8TeV-powheg-pythia6/weizou-SQWaT_PAT_53X_Summer12_v2-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/VBF_HToWWToLAndTauNuQQ_M-700_8TeV-powheg-pythia6/weizou-SQWaT_PAT_53X_Summer12_v2-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/VBF_HToWWToLAndTauNuQQ_M-800_8TeV-powheg-pythia6/weizou-SQWaT_PAT_53X_Summer12_v2-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/VBF_HToWWToLAndTauNuQQ_M-900_8TeV-powheg-pythia6/weizou-SQWaT_PAT_53X_Summer12_v2-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/VBF_HToWWToLAndTauNuQQ_M-1000_8TeV-powheg-pythia6/weizou-SQWaT_PAT_53X_Summer12_v2-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/WJetsToLNu_PtW-100_8TeV-herwigpp/ajkumar-SQWaT_PAT_53X_Summer12_v2-829f288d768dd564418efaaf3a8ab9aa/USER",  
              #"/WJetsToLNu_matchingdown_8TeV-madgraph-tauola/aperloff-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/WJetsToLNu_matchingup_8TeV-madgraph-tauola/aperloff-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/WJetsToLNu_scaledown_8TeV-madgraph-tauola/aperloff-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              #"/WJetsToLNu_scaleup_8TeV-madgraph-tauola/aperloff-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/TTJets_matchingup_TuneZ2star_8TeV-madgraph-tauola/ajkumar-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/TTJets_scaledown_TuneZ2star_8TeV-madgraph-tauola/ajkumar-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/TTJets_matchingdown_TuneZ2star_8TeV-madgraph-tauola/jdamgov-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/TTJets_scaleup_TuneZ2star_8TeV-madgraph-tauola/ajkumar-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/TT_CT10_TuneZ2star_8TeV-powheg-tauola/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/WJetsToLNu_PtW-180_TuneZ2star_8TeV-madgraph-tarball/ajkumar-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl02_M-1000_TuneZ2star_8TeV-pythia6/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl02_M-1500_TuneZ2star_8TeV-pythia6/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl02_M-2000_TuneZ2star_8TeV-pythia6/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl01_M-1000_Tune23_8TeV-herwigpp/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl01_M-1500_Tune23_8TeV-herwigpp/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl01_M-1000_TuneZ2star_8TeV-pythia6/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl01_M-1500_TuneZ2star_8TeV-pythia6/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/RSGravitonToWW_kMpl01_M-2000_TuneZ2star_8TeV-pythia6/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M600-JHU-v2/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M700-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M800-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M900-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1000-JHU-v3/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1100-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1200-JHU-v2/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1300-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1400-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1500-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1600-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1700-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1800-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M1900-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M2000-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M2100-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M2200-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M2300-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M2400-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",
              "/BulkG_WW_lvjj_c0p2_M2500-JHU-v1/custodio-SQWaT_PAT_53X_Summer12_v1-829f288d768dd564418efaaf3a8ab9aa/USER",              

    ]
channels   = [
    #"WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_part1",
    #"WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_part2",
    #"W1JetsToLNu_TuneZ2Star_8TeV-madgraph",    
    #"W2JetsToLNu_TuneZ2Star_8TeV-madgraph",
    #"W3JetsToLNu_TuneZ2Star_8TeV-madgraph",            
    #"W4JetsToLNu_TuneZ2Star_8TeV-madgraph",
    "WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph",
    "DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball", 
    "TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola",
    "WZ_TuneZ2star_8TeV_pythia6_tauola",
    "ZZ_TuneZ2star_8TeV_pythia6_tauola",
    "T_t-channel_TuneZ2star_8TeV-powheg-tauola",
    "Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola",
    "T_s-channel_TuneZ2star_8TeV-powheg-tauola",
    "Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola",    
    "T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola",
    "Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola",
    "WW_TuneZ2star_8TeV_pythia6_tauola",
    #"GluGluToHToWWToLAndTauNuQQ_M-600_8TeV-powheg-pythia6",
    #"GluGluToHToWWToLAndTauNuQQ_M-700_8TeV-powheg-pythia6",
    #"GluGluToHToWWToLAndTauNuQQ_M-800_8TeV-powheg-pythia6",
    #"GluGluToHToWWToLAndTauNuQQ_M-900_8TeV-powheg-pythia6",
    #"GluGluToHToWWToLAndTauNuQQ_M-1000_8TeV-powheg-pythia6",
    #"qqToHToWWToLAndTauNuQQ_M-600_8TeV-powheg-pythia6",
    #"qqToHToWWToLAndTauNuQQ_M-700_8TeV-powheg-pythia6",
    #"qqToHToWWToLAndTauNuQQ_M-800_8TeV-powheg-pythia6",
    #"qqToHToWWToLAndTauNuQQ_M-900_8TeV-powheg-pythia6",
    #"qqToHToWWToLAndTauNuQQ_M-1000_8TeV-powheg-pythia6"
    "WJetsToLNu_PtW-100_8TeV-herwigpp",
    #"WJetsToLNu_matchingdown_8TeV-madgraph-tauola",
    #"WJetsToLNu_matchingup_8TeV-madgraph-tauola",
    #"WJetsToLNu_scaledown_8TeV-madgraph-tauola",
    #"WJetsToLNu_scaleup_8TeV-madgraph-tauola",
    "TTJets_matchingup_TuneZ2star_8TeV-madgraph-tauola",
    "TTJets_scaledown_TuneZ2star_8TeV-madgraph-tauola",
    "TTJets_matchingdown_TuneZ2star_8TeV-madgraph-tauola",
    "TTJets_scaleup_TuneZ2star_8TeV-madgraph-tauola",
    "TT_CT10_TuneZ2star_8TeV-powheg-tauola",
    "WJetsToLNu_PtW-180_TuneZ2star_8TeV-madgraph-tarball",
    "RSGravitonToWW_kMpl02_M-1000_TuneZ2star_8TeV-pythia6",
    "RSGravitonToWW_kMpl02_M-1500_TuneZ2star_8TeV-pythia6",
    "RSGravitonToWW_kMpl02_M-2000_TuneZ2star_8TeV-pythia6",
    "RSGravitonToWW_kMpl01_M-1000_Tune23_8TeV-herwigpp",
    "RSGravitonToWW_kMpl01_M-1500_Tune23_8TeV-herwigpp",
    "RSGravitonToWW_kMpl01_M-1000_TuneZ2star_8TeV-pythia6",
    "RSGravitonToWW_kMpl01_M-1500_TuneZ2star_8TeV-pythia6",
    "RSGravitonToWW_kMpl01_M-2000_TuneZ2star_8TeV-pythia6",    
    "BulkG_WW_lvjj_c0p2_M600-JHU",
    "BulkG_WW_lvjj_c0p2_M700-JHU",
    "BulkG_WW_lvjj_c0p2_M800-JHU",
    "BulkG_WW_lvjj_c0p2_M900-JHU",
    "BulkG_WW_lvjj_c0p2_M1000-JHU",
    "BulkG_WW_lvjj_c0p2_M1100-JHU",
    "BulkG_WW_lvjj_c0p2_M1200-JHU",
    "BulkG_WW_lvjj_c0p2_M1300-JHU",
    "BulkG_WW_lvjj_c0p2_M1400-JHU",
    "BulkG_WW_lvjj_c0p2_M1500-JHU",
    "BulkG_WW_lvjj_c0p2_M1600-JHU",
    "BulkG_WW_lvjj_c0p2_M1700-JHU",
    "BulkG_WW_lvjj_c0p2_M1800-JHU",
    "BulkG_WW_lvjj_c0p2_M1900-JHU",
    "BulkG_WW_lvjj_c0p2_M2000-JHU",
    "BulkG_WW_lvjj_c0p2_M2100-JHU",
    "BulkG_WW_lvjj_c0p2_M2200-JHU",
    "BulkG_WW_lvjj_c0p2_M2300-JHU",
    "BulkG_WW_lvjj_c0p2_M2400-JHU",
    "BulkG_WW_lvjj_c0p2_M2500-JHU",
    
]
condor   = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
MyResilientArea = "/uscms_data/d3/rgerosa1/CMSSW_5_3_5/src/ElectroWeakAnalysis/VPlusJets/test/DataProcessing/"+physMode


def changeMainConfigFile(trigpath,isQCD,isHEEPID,isTransverseMassCut):
    fin  = open(ConfigFile)
    pset_cfg      = physMode+"py_" + trigpath + ".py"
    outfile_root  = physMode + trigpath + ".root"
    fout = open(pset_cfg,"w")
    for line in fin.readlines():
        if  line.find("isMC = False")!=-1:
            line=line.replace("isMC = False", "isMC = True")

        if  line.find("isQCD = True")!=-1:
            line=line.replace("isQCD = True",isQCD)
        if  line.find("isQCD = False")!=-1:
            line=line.replace("isQCD = False",isQCD)

        if  line.find("isHEEPID = True")!=-1:
            line=line.replace("isHEEPID = True",isHEEPID)
        if  line.find("isHEEPID = False")!=-1:
            line=line.replace("isHEEPID = False",isHEEPID)

        if  line.find("isTransverseMassCut = True")!=-1:
            line=line.replace("isTransverseMassCut = True",isTransverseMassCut)
        if  line.find("isTransverseMassCut = False")!=-1:
            line=line.replace("isTransverseMassCut = False",isTransverseMassCut)

        if  line.find("WmunuJetAnalysisntuple.root")!=-1:
            line=line.replace("WmunuJetAnalysisntuple.root",outfile_root)
        fout.write(line)
    print pset_cfg + " has been written.\n"

                                                                                    

def changeCrabTemplateFile(outfile, index):
    fin  = open("crabTemplateMC.cfg")
    pset_cfg      = physMode+"py_" + outfile + ".py"
    pset_crab     = physMode+"cb_" + outfile + ".cfg"
    outfile_root  = physMode + outfile + ".root"
    fout = open(pset_crab,"w")
    for line in fin.readlines():
        if  line.find("mydataset")!=-1:
            line=line.replace("mydataset",dataset[index])
            fout.write("\n")
        if  line.find("total_number_of_lumis")!=-1:
            line=line.replace("total_number_of_lumis=-1","total_number_of_events=-1")
        if  line.find("lumis_per_job")!=-1:
            line=line.replace("lumis_per_job = 100","events_per_job = 5000")     
        if line.find("myanalysis")!=-1:
            line=line.replace("myanalysis",pset_cfg)    
        if  line.find("myrootfile")!=-1:
            line=line.replace("myrootfile",outfile_root)
            fout.write("dbs_url = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_01/servlet/DBSServlet"+"\n")
        if  line.find("myresilient")!=-1:
            line=line.replace("myresilient",MyResilientArea+outfile)    
        if line.find("glite")!=-1 and condor[index]!=0:
            line=line.replace("glite", "condor")
        if line.find("outputDir")!=-1 :
            line=line.replace("outputDir",outfile)
        fout.write(line)        
    if condor[index]!=0:
        fout.write("ce_white_list = cmssrm.fnal.gov")
      
    print pset_crab + " has been written.\n"


    
###################
for i in range(len(channels)):
    changeMainConfigFile(channels[i],isQCD,isHEEPID,isTransverseMassCut)
    changeCrabTemplateFile(channels[i],i)

for i in range(len(channels)):
    submitcommand = "crab -create -cfg " + physMode + "cb_" + channels[i] + ".cfg"    
    child   = os.system(submitcommand)

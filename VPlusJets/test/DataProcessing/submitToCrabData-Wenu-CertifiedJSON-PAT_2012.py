import os,sys
import string, re
from time import gmtime, localtime, strftime

physMode   = "WenuJets_"
ConfigFile = "../WenuJetsAnalysisPAT_cfg.py"
DefTrig    = "'HLT_Ele80_CaloIdVT_*'"
isMC       = "isMC = False"
isQCD      = "isQCD = False"
isHEEPID   = "isHEEPID = True"
isTransverseMassCut = "isTransverseMassCut = False"

dataset    = ["/SingleElectron/jdamgov-SQWaT_PAT_53X_Run2012A-recover-06Aug2012-v1-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/ajkumar-SQWaT_PAT_53X_2012A-13Jul2012-v1-dee4c99a1b5d294b4043f483391f854a/USER",
              "/SingleElectron/dimatteo-SQWaT_PAT_53X_2012B-13Jul2012-v1_part1-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/custodio-SQWaT_PAT_53X_2012B-13Jul2012-v1_part2-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/dimatteo-SQWaT_PAT_53X_2012C-24Aug2012-v1-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/ajkumar-SQWaT_PAT_53X_2012C-PromptReco-v2-dee4c99a1b5d294b4043f483391f854a/USER",
              "/SingleElectron/ajkumar-SQWaT_PAT_53X_2012C-PromptReco-v2_pt2-dee4c99a1b5d294b4043f483391f854a/USER",
              "/SingleElectron/ilyao-SQWaT_PAT_53X_SingleElectron_Run2012C-PromptReco-v2-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/jdamgov-SQWaT_PAT_53X_2012D_pt1-PromptReco-v1-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/jdamgov-SQWaT_PAT_53X_2012D_pt2-PromptReco-v1-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/jdamgov-SQWaT_PAT_53X_2012D_pt3-PromptReco-v1-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/custodio-SQWaT_PAT_53X_2012D-PromptReco-v1_p4-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/custodio-SQWaT_PAT_53X_2012D-PromptReco-v1_p5-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/ajkumar-SQWaT_PAT_53X_2012D-PromptReco-v1_p6-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/custodio-SQWaT_PAT_53X_2012D-PromptReco-v1_p7-3e4086321697e2c39c90dad08848274b/USER",
              "/SingleElectron/custodio-SQWaT_PAT_53X_2012D-PromptReco-v1_p8-3e4086321697e2c39c90dad08848274b/USER"
              ]
channels   = ["SingleElectron_Run2012A-recover-06Aug2012-v1",
              "SingleElectron_Run2012A-13Jul2012-v1",
              "SingleElectron_Run2012B-13Jul2012-v1_part1",
              "SingleElectron_Run2012B-13Jul2012-v1_part2",
              "SingleElectron_Run2012C-24Aug2012-v1",
              "SingleElectron_Run2012C-PromptReco-v2_part1",
              "SingleElectron_Run2012C-PromptReco-v2_part2",
              "SingleElectron_Run2012C-PromptReco-v2_part3",
              "SingleElectron_Run2012D-PromptReco-v1_part1",
              "SingleElectron_Run2012D-PromptReco-v1_part2",
              "SingleElectron_Run2012D-PromptReco-v1_part3",
              "SingleElectron_Run2012D-PromptReco-v1_part4",
              "SingleElectron_Run2012D-PromptReco-v1_part5",
              "SingleElectron_Run2012D-PromptReco-v1_part6",
              "SingleElectron_Run2012D-PromptReco-v1_part7",
              "SingleElectron_Run2012D-PromptReco-v1_part8"
              ]
trigname   = [
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'",
               "'HLT_Ele80_CaloIdVT_*'"]

RunRange   = [
        "190782-190949",
        "190456-196531",
        "193834-195182",
        "195183-196531",
        "190645-208686",
        "190645-208686",
        "200519-202016",
        "194631-203002",
        "203894-205618",
        "205620-206539",
        "205339-206088",
        "200961-206940",
        "206744-207469",
        "207477-207898",
        "207889-208357",
        "190645-208686",
]

JSON       = [
        "json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt",
        "json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt",
        "json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt",
        "json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt",
        "json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_200519-202016_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_194631-203002_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_203894-205618_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_205620-206539_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_205339-206088_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_200961-206940_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_206744-207469_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_207477-207898_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_207889-208357_8TeV_PromptReco_Collisions12_JSON.txt",
        "json/Cert_194480-208686_8TeV_PromptReco_Collisions12_JSON.txt",
        ]


condor     = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] # Total jobs 13 now
MyResilientArea = "/uscms_data/d3/rgerosa1/CMSSW_5_3_3_patch3/src/ElectroWeakAnalysis/VPlusJets/test/DataProcessing/"+physMode

## ------------------------------------------


def changeMainConfigFile(trigpath,nowtrigname,MCFlag, isQCD, isHEEPID,isTransverseMassCut):
    fin  = open(ConfigFile)
    pset_cfg      = physMode + "py_" + trigpath + ".py"
    outfile_root  = physMode + trigpath + ".root"
    fout = open(pset_cfg,"w")
    for line in fin.readlines():
        if  line.find("isMC = True")!=-1:
            line=line.replace("isMC = True",MCFlag)
            
        if  line.find("isQCD = True")!=-1: 
            line=line.replace("isQCD = True",isQCD)
        if  line.find("isQCD = False")!=-1:
            line=line.replace("isQCD = False",isQCD)

        if  line.find("isHEEPID = True")!=-1: 
            line=line.replace("isHEEPID = True",isHEEPID)
        if  line.find("isHEEPID = False")!=-1:
            line=line.replace("isHEEPID = False",isHEEPID)
            
        if  line.find("isTransverseMassCut = True")!=-1: 
            line=line.replace("isTransverseMassCut = True", isTransverseMassCut)
        if  line.find("isTransverseMassCut = False")!=-1:
            line=line.replace("isTransverseMassCut = False", isTransverseMassCut)

        if line.find("WenuJetAnalysisntuple.root")!=-1:
            line=line.replace("WenuJetAnalysisntuple.root",outfile_root)
        if  line.find(DefTrig)!=-1:
            line=line.replace(DefTrig,nowtrigname)
        fout.write(line)
    print pset_cfg + " has been written.\n"


def changeCrabTemplateFile(outfile, index, nowJSON, nowRang):
    fin  = open("crabTemplateData.cfg")
    pset_cfg      = physMode + "py_" + outfile + ".py"
    pset_crab     = physMode + "cb_" + outfile + ".cfg"
    outfile_root  = physMode + outfile + ".root"
    fout = open(pset_crab,"w")
    for line in fin.readlines():
        if  line.find("mydataset")!=-1:
            line=line.replace("mydataset",dataset[index])
            fout.write("\n")
            fout.write("runselection="+RunRange[index]+"\n")
            fout.write("lumi_mask="+nowJSON+"\n")
        if line.find("myanalysis")!=-1:
            line=line.replace("myanalysis",pset_cfg)    
            fout.write("dbs_url = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet"+"\n")
        if  line.find("myrootfile")!=-1:
            line=line.replace("myrootfile",outfile_root)
        if  line.find("myresilient")!=-1:
            line=line.replace("myresilient",MyResilientArea+nowRang+"_"+outfile)    
        if line.find("glite")!=-1 and condor[index]!=0:
            line=line.replace("glite", "condor")
        if line.find("outputDir")!=-1 :
            line=line.replace("outputDir",outfile)            
        fout.write(line)        
    if condor[index]!=0:
       fout.write("ce_white_list = cmssrm.fnal.gov")
        
    print pset_crab + " has been written.\n"

    
####################
## move along channels vector

for i in range(len(channels)):
    changeMainConfigFile(channels[i],trigname[i],isMC,isQCD, isHEEPID,isTransverseMassCut)
    changeCrabTemplateFile(channels[i],i,JSON[i],RunRange[i])

for i in range(len(channels)):
    submitcommand = "crab -create -cfg " + physMode + "cb_" + channels[i] + ".cfg"
    child   = os.system(submitcommand)

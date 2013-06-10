import os
import glob
import math
import time

import ROOT
from ROOT import *
from ROOT import gROOT, gStyle, gSystem, TLatex
import subprocess

from subprocess import Popen
from optparse import OptionParser

#python herculesUtils.py -b --batchMode --createTree --sampleToProcess all  --channel mu --queque longcms
############################################################
############################################
#            Job steering                  #
############################################
parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')

parser.add_option('--createTree', action='store_true', dest='createTree', default=False,
                  help='no X11 windows')

# submit jobs to condor
parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False,
                  help='no X11 windows')
parser.add_option('--sampleToProcess',action="store",type="string",dest="sampleToProcess",default=None)

parser.add_option('--channel',action="store",type="string",dest="channel",default="mu")

parser.add_option('--queque',action="store",type="string",dest="queque",default="longcms")


(options, args) = parser.parse_args()
############################################################
############################################################


def submitBatchJob( command, fn ):

    currentDir = os.getcwd();


    # create a dummy bash/csh
    outScript=open(fn+".sh","w");

    outScript.write('#!/bin/bash');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+command);
    outScript.close();

    # submit the condor job 

#    print "qsub -V -d "+currentDir+" -q shortcms "+currentDir+"/"+fn+".sh " ;
    Lunch = "lancia";
    outLunch=open(Lunch+".sh","a");
    outLunch.write("qsub -V -d "+currentDir+" -q "+options.queque+" "+currentDir+"/"+fn+".sh \n") ;
    
#    time.sleep(0.5);
#    os.system("mv "+"condor_"+fn+" "+fn+".sh"+" condorTmp/.");


if __name__ == '__main__':

    all = [
#           "data",
#            "data_194712_559_405404310", 
#            "data_204544_102_152417542",
#            "data_195930_83_62880579",
#            "data_208487_118_209515269",
            "data_195013_114_117238404",
            "data_195948_317_517532532",
            "data_202016_935_952022882",
            "data_195397_899_1123713321",
            "data_201278_532_72515017",
            "data_202973_241_260228320",              
            "data_195655_60_73510965",
            "data_201602_497_68268345",
#            "data_190703_11_8519680",  
#            "data_196438_586_479961054",
#            "data_191834_99_120189151",
#            "data_196452_193_232654780",
#            "data_gr2_194199_190_149891910",
#            "data_194199_258_229967752",
#            "data_198522_104_75971527",
#            "data_gr2_194897_50_91029885",
#            "data_194699_80_97577634", 
#            "data_199428_457_557329806", 
#            "data_gr2_195397_558_772053348", 
#            "data_199428_86_8009087", 
#            "data_gr2_199608_1123_1247460423", 
#            "data_194778_190_248720484", 
#            "data_199961_177_187360166", 
#            "data_gr2_200091_1678_1762413892", 
#            "data_194912_390_657997761",
#            "data_200188_104_142020231", 
#            "data_gr2_202299_249_336447785", 
#            "data_195013_303_44790264", 
#            "data_203002_1517_1719001770",
#            "data_gr2_202314_128_132710041",
#            "data_195099_171_245228003",
#            "data_gr2_206476_128_148892985",
#            "data_195398_899_752941172",
#            "data_204564_359_400306553",
#            "data_gr2_206859_764_1094226833",
#            "data_195399_159_130828536",
#            "data_207372_380_547252863",
#            "data_gr2_207231_920_1260547040",
#            "data_207490_80_87321589",
#             "data_201191_317_488053419",
#           "ggH600",
#           "ggH700",
#           "ggH800",
#           "ggH900",           
#           "ggH1000",
#           "vbfH600",
#           "vbfH700",
#           "vbfH800",
#           "vbfH900",           
#           "vbfH1000",
#           "rsg1000_kMpl01_py",
#           "rsg1000_kMpl01_hw",
#           "rsg1500_kMpl01_py",
#           "rsg1500_kMpl01_hw",
#           "rsg2000_kmpl01_py",
#           "BulkG_c0p2_M600",
#           "BulkG_c0p2_M700",
#           "BulkG_c0p2_M800",
#           "BulkG_c0p2_M900",
#           "BulkG_c0p2_M1000",
#           "BulkG_c0p2_M1100",
#           "BulkG_c0p2_M1200",
#           "BulkG_c0p2_M1300",
#           "BulkG_c0p2_M1400",
#           "BulkG_c0p2_M1500",
#           "BulkG_c0p2_M1600",
#           "BulkG_c0p2_M1700",
#           "BulkG_c0p2_M1800",
#           "BulkG_c0p2_M1900",
#           "BulkG_c0p2_M2000",
#           "BulkG_c0p2_M2100",
#           "BulkG_c0p2_M2200",
#           "BulkG_c0p2_M2400",
#           "BulkG_c0p2_M2500",
#           "WJets_Pythia",
#           "WJets_Herwig",
#           "WJets_Pythia180",
#           "ZJets",
#           "TTbar",
#           "TTbar_Powheg",
#           "TTbar_matchDn",           
#           "TTbar_matchUp",                      
#           "TTbar_scaleDn",           
#           "TTbar_scaleUp",                      
#           "WW",           
#           "WZ",                      
#           "ZZ",
#           "tch",
#           "tWch",
#           "sch",
#           "tch_bar",
#           "tWch_bar",
#           "sch_bar"           
           ]


    os.system("rm lancia.sh");
    os.system("touch lancia.sh");
        
    if options.sampleToProcess == "all" and options.batchMode and options.createTree:
    
        for i in range(len(all)):
            
            cmmd = "python runAnalysis.py -b --createTrees --channel "+options.channel+" --sampleToProcess "+ all[i]
            print cmmd
            fn = "createTreeScript_%s_%s"%(options.channel,all[i]);
            submitBatchJob( cmmd, fn );

    
    elif options.createTree and options.batchMode and not options.sampleToProcess == None:

        cmmd = "python runAnalysis.py -b --createTrees --channel "+options.channel+" --sampleToProcess "+ options.sampleToProcess
        print cmmd
        fn = "createTreeScript_%s_%s"%(options.channel,options.sampleToProcess);
        submitBatchJob( cmmd, fn );

    else:
        print "do nothing"







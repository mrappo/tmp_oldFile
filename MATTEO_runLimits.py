#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT

from optparse import OptionParser
from subprocess import Popen
from ROOT import gROOT, gStyle, gSystem, TLatex, TGaxis, TPaveText, TH2D, TColor, gPad, TGraph2D, TLine,TGraph,TList,TPad
import ROOT as rt

ROOT.gStyle.SetPadRightMargin(0.16);

from collections import defaultdict

############################################
#            Job steering                  #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-f', '--treeFolder',              action="store", type="string", dest="treeFolder",              default="./")

### basic allowed methods
parser.add_option('--makeCards',     action='store_true', dest='makeCards',         default=False, help='options to produce datacards and workspaces via fitting analysis')
parser.add_option('--computeLimits', action='store_true', dest='computeLimits',     default=False, help='basic option in order to compute asymptotic limits')
parser.add_option('--plotLimits',    action='store_true', dest='plotLimits',        default=False, help='basic option to plot asymptotic and pvalue')
parser.add_option('--computePvalue', action='store_true', dest='computePvalue',     default=False, help='basic option to plot asymptotic and pvalue')
parser.add_option('--computeSignif', action='store_true', dest='computeSignif',     default=False, help='basic option to plot asymptotic and pvalue')
parser.add_option('--biasStudy',     action='store_true', dest='biasStudy',         default=False, help='basic option to perform bias study with our own tool')
parser.add_option('--maximumLikelihoodFit', action='store_true', dest='maximumLikelihoodFit', default=False, help='basic option to run max Likelihood fit inside combination tool')
parser.add_option('--asymptotic', action='store_true', dest='asymptotic', default=False, help='basic option to run asymptotic CLS inside combination tool')
parser.add_option('--fullCLs', action='store_true', dest='fullCLs', default=False, help='basic option to run full CLS inside combination tool')
parser.add_option('--generateOnly',  action='store_true', dest='generateOnly', default=False, help='basic option to run the only generation with combiner')
parser.add_option('--makeLikelihoodScan', action='store_true', dest='makeLikelihoodScan', default=False, help='basic option to run the likelihood scan')

##### submit jobs to condor, lxbatch and hercules 
parser.add_option('--batchMode',      action='store_true', dest='batchMode',      default=False, help='to run jobs on condor fnal')
parser.add_option('--lxbatchCern',    action='store_true', dest='lxbatchCern',    default=False, help='run jobs on lxbatch system at cern, default is condor fnal')
parser.add_option('--herculesMilano', action='store_true', dest='herculesMilano', default=False, help='run jobs on hercules system at Milano, default is condor fnal')

##### other basci options for all the methods 
parser.add_option('--datacardDIR', action="store", type="string", dest="datacardDIR", default="")
parser.add_option('--category',    action="store", type="string", dest="category",    default="HP")
parser.add_option('--channel',     action="store", type="string", dest="channel",     default="em")
parser.add_option('--queque',      action="store", type="string", dest="queque",      default="")
parser.add_option('--pseudodata',  action="store", type="int",    dest="pseudodata",  default=0)
parser.add_option('--systematics', action="store", type="int",    dest="systematics", default=1)
parser.add_option('--massPoint',   action="store", type="int",    dest="massPoint",   default=-1)
parser.add_option('--closuretest', action="store", type="int",    dest="closuretest", default=0)
parser.add_option('--cPrime',      action="store", type="int",    dest="cPrime",      default=-1)
parser.add_option('--brNew',       action="store", type="int",    dest="brNew",       default=-1)
parser.add_option('--odir',        action="store", type="string", dest="odir",        default=".")
parser.add_option('--sigChannel',  action="store", type="string", dest="sigChannel",  default="")
parser.add_option('--jetBin',      action="store", type="string", dest="jetBin",      default="")
parser.add_option('--skipJetSystematics',   action="store", type="int", dest="skipJetSystematics",    default=0)
parser.add_option('--turnOnAnalysis',       action="store", type="int",  dest="turnOnAnalysis",       default=0)
parser.add_option('--injectSingalStrenght', action="store", type=float,  dest="injectSingalStrenght", default=0., help='inject a singal in the toy generation')

parser.add_option('--higgsCombination', action="store", type="int", dest="higgsCombination", default=0)
parser.add_option('--interferenceModel', action="store", type="string", dest="interferenceModel", default="3")


###### options for Bias test in the combination tool
parser.add_option('--nToys',        action="store", type="int",    dest="nToys",       default=0)
parser.add_option('--crossedToys',  action="store", type="int",    dest="crossedToys", default=0)
parser.add_option('--inputGeneratedDataset', action="store", type="string",    dest="inputGeneratedDataset", default="")
parser.add_option('--outputTree',   action="store", type="int",    dest="outputTree",  default=0)


##### options specific for the bias tool 
parser.add_option('--shapetest',           action="store", type="int",    dest="shapetest",           default=0)
parser.add_option('--ttbarcontrolregion',  action="store", type="int",    dest="ttbarcontrolregion",  default=0)
parser.add_option('--mlvjregion',          action="store", type="string", dest="mlvjregion",          default="_sb_lo")
parser.add_option('--fitjetmass',          action="store", type="int",    dest="fitjetmass",          default=0)
parser.add_option('--onlybackgroundfit',   action="store", type="int",    dest="onlybackgroundfit",   default=0, help='run only background fit')
parser.add_option('--inflatejobstatistic', action="store", type="int",    dest="inflatejobstatistic", default=1, help='enlarge the generated statistics in the fit')
parser.add_option('--scalesignalwidth',    action="store", type="int",    dest="scalesignalwidth",    default=1, help='reduce the signal width by a factor x')

##### final plot options
parser.add_option('--makeSMLimitPlot',       action="store", type="int",    dest="makeSMLimitPlot",        default=0)
parser.add_option('--makeBSMLimitPlotMass',  action="store", type="int",    dest="makeBSMLimitPlotMass",   default=0)
parser.add_option('--makeBSMLimitPlotBRnew', action="store", type="int",    dest="makeBSMLimitPlotBRnew",  default=0)
parser.add_option('--makeBSMLimitPlot2D',    action="store", type="int",    dest="makeBSMLimitPlot2D",     default=0)
parser.add_option('--blindObservedLine',     action="store", type="int",    dest="blindObservedLine",      default=0)
parser.add_option('--plotPValue',            action="store", type="int",    dest="plotPValue",             default=0)
parser.add_option('--plotxsec',              action="store", type="int",    dest="plotxsec",             default=0)
parser.add_option('--plotSignalStrenght',    action="store", type="int",    dest="plotSignalStrenght",     default=0)
parser.add_option('--plotLikelihoodScan',    action="store", type="int",    dest="plotLikelihoodScan",     default=0)

###############################################
####### MATTEO ADDED
###############################################

#parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--ntuple', action="store",type="string",dest="ntuple",default="WWTree_22sep_jecV7_lowmass")
#parser.add_option('--category', action="store",type="string",dest="category",default="HP")
parser.add_option('--sample', action="store",type="string",dest="sample",default="BulkGraviton")
#parser.add_option('--jetalgo', action="store",type="string",dest="jetalgo",default="jet_mass_pr")
#parser.add_option('--interpolate', action="store_true",dest="interpolate",default=False)
#parser.add_option('--batchMode', action="store_true",dest="batchMode",default=False)
parser.add_option('--vbf', action="store_true",dest="VBF_process",default=False)




(options, args) = parser.parse_args()

#############################################
######### Get Some Global Variables #########
#############################################
'''
##########################
### ORIGINAL CODE
##########################
#mass = [1000,2000,3000,4000]
#mass = [4000]
#mass = [4000]
mass = [800]
#mass = [750]
ccmlo = [600]
#ccmhi = [1500,2500,3500,4500]
ccmhi = [1500]
mjlo = [40]
mjhi = [150]
#mjlo = [35,35,35,35]
#mjhi = [110,110,110,110]
mlo = [600]
#mlo = [1000,1000,1000,1000]
#mhi = [1500,2500,3500,4500]
mhi = [1500]
#shape = ["ExpN","ExpN","ExpN","ExpN","ExpN","ExpN"]
shape = ["Exp"]
#shape = ["ExpN","ExpN","ExpN","ExpN"]
shapeAlt = ["ExpTail"]
#shape = ["Exp","Exp","Exp","Exp"]
#shapeAlt = ["Pow","Pow","Pow","Pow"]
'''

#####################################
##### MATTEO CHANGED
#####################################
sample_name=options.sample


if options.VBF_process:
   if sample_name.find('BulkGraviton') !=-1:
      mass=[600,800,1000]
      ccmlo = [400,600,800]
      ccmhi = [1500,2000,2500]
      mjlo = [40,40,40]
      mjhi = [150,150,150]
      mlo = [400,600,800]
      mhi = [1500,2000,2500]
      shape = ["Exp","Exp","Exp"]
      shapeAlt = ["ExpTail","ExpTail","ExpTail"]
      xsec_value = [0.010939,0.0021845,0.0006584]
      xsec_corr_value = [1.,1.,1.] 

     
       
   if sample_name.find('Higgs') !=-1:
      mass=[650,1000]
      ccmlo = [400,600]
      ccmhi = [1500,2000]
      mjlo = [40,40]
      mjhi = [150,150]
      mlo = [400,600]
      mhi = [1500,2000]
      shape = ["Exp","Exp"]
      shapeAlt = ["ExpTail","ExpTail"]
      xsec_value = [0.067953,0.023859]
      xsec_corr_value = [1.,1.,1.]
   
   
   
else:   
   if sample_name.find('BulkGraviton') !=-1:
      mass=[600,800,1000]
      ccmlo = [400,600,800]
      ccmhi = [1500,2000,2500]
      mjlo = [40,40,40]
      mjhi = [150,150,150]
      mlo = [200,600,800]
      mhi = [1500,2000,2500]
      shape = ["Exp","Exp","Exp"]
      shapeAlt = ["Exp","ExpTail","ExpTail"]
      xsec_value = [0.40683,0.07605,0.020499]
      xsec_corr_value = [1.,1.,1.] 

     
       
   if sample_name.find('Higgs') !=-1:
      mass=[650,1000]
      ccmlo = [200,200]
      ccmhi = [1500,3000]
      mjlo = [40,40]
      mjhi = [150,150]
      mlo = [200,200]
      mhi = [1500,3000]
      shape = ["Exp","Exp"]
      shapeAlt = ["ExpTail","ExpTail"]
      xsec_value = [0.19579,0.05041]
      xsec_corr_value = [1.,1.,1.]

points = []  
for p in range(1,2):
   points+=[float(p/10.)]
   points+=[float(p/10.+0.05)]
   points+=[float(p/1.)]
   points+=[float(p/1.+0.5)]
   points+=[float(p*10.)]
   points+=[float(p*10.+5.)]

'''
mass = [2000]#,2000,3000,4000]
ccmlo = [700]#,800,800,800]
ccmhi = [2500]#,4800,4800,4800]
mjlo = [35]#,40,40,40]
mjhi = [140]#,130,130,130]
mlo = [700]#,700,700,700]
mhi = [2500]#,5000,5000,5000]
shape = ["ExpN"]#,"Exp","Exp","Exp"]
shapeAlt = ["ExpTail"]#,"Pow","Pow","Pow"]
'''
################## options turnOn Analysis
'''
if options.turnOnAnalysis :

 mlo      =  [ 400, 400, 400, 400, 400] ## min mlvj cut 
 mhi      =  [1500,1500,1500,1500,1500] ## max mlvj cut

else: 

 mlo      =  [ 550, 550, 550, 550, 550] ## min mlvj cut
 mhi      =  [1500,1500,1500,1500,1500] ## max mlvj cut
'''

################## options for makeCards
''' 
if options.makeCards and options.turnOnAnalysis:
 
 shapeAlt =  ["ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1"] ## basic shape
 shape    =  ["ErfPow2_v1", "ErfPow2_v1","ErfPow2_v1","ErfPow2_v1","ErfPow2_v1"] ## alternate one

elif not options.turnOnAnalysis and options.makeCards:

 shape      =  ["Exp","Exp","Exp","Exp","Exp"] ## basic shape
 shapeAlt   =  ["Pow","Pow","Pow","Pow","Pow"] ## alternate one
'''
################## options for bias Study

if options.biasStudy:

 if not options.turnOnAnalysis and not options.fitjetmass:

#  shape_gen = ["Exp","Exp","Exp","Exp","Exp"]    
  shape_fit = ["Exp","Exp","Exp","Exp","Exp"]
  shape_gen = ["Pow2","Pow2","Pow2","Pow2","Pow2"]    
#  shape_fit = ["Pow2","Pow2","Pow2","Pow2","Pow2"]
#  shape_gen = ["Pow","Pow","Pow","Pow","Pow"]    
#  shape_fit = ["Pow","Pow","Pow","Pow","Pow"]

 elif options.turnOnAnalysis and not options.fitjetmass:

  shape_gen = ["ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1"];    
  shape_fit = ["ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1"];
#  shape_gen = ["ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1"];    
#  shape_fit = ["ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1"];
#  shape_gen = ["ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1"];    
#  shape_fit = ["ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1"];

 elif options.fitjetmass:

#  shape_gen = ["ErfExp","ErfExp","ErfExp","ErfExp","ErfExp"];    
#  shape_fit = ["ErfExp","ErfExp","ErfExp","ErfExp","ErfExp"];
#  shape_gen = ["User1","User1","User1","User1","User1"];    
  shape_fit = ["User1","User1","User1","User1","User1"];
  shape_gen = ["ErfPow","ErfPow","ErfPow","ErfPow","ErfPow"];    
#  shape_fit = ["ErfPow","ErfPow","ErfPow","ErfPow","ErfPow"];

 isMC      = [0,0,0,0,0];

#cprime = [01,02,03,05,07,10]
#BRnew  = [00,01,02,03,04,05]
cprime = [10]
BRnew = [00]

####CMS lumi

cmsText     = "CMS";
cmsTextFont   = 61  

writeExtraText = True
extraText   = "Preliminary"
extraTextFont = 52 

lumiTextSize     = 0.6
lumiTextOffset   = 0.2

cmsTextSize      = 0.75
cmsTextOffset    = 0.1

relPosX    = 0.045
relPosY    = 0.035
relExtraDY = 1.2

extraOverCmsTextSize  = 0.76

lumi_13TeV = "2.3 fb^{-1}"
lumi_8TeV  = "19.7 fb^{-1}" 
lumi_7TeV  = "5.1 fb^{-1}"
lumi_sqrtS = ""

drawLogo      = False

def CMS_lumi(pad,  iPeriod,  iPosX ):
    outOfFrame    = False
    if(iPosX/10==0 ): outOfFrame = True

    alignY_=3
    alignX_=2
    if( iPosX/10==0 ): alignX_=1
    if( iPosX==0    ): alignY_=1
    if( iPosX/10==1 ): alignX_=1
    if( iPosX/10==2 ): alignX_=2
    if( iPosX/10==3 ): alignX_=3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025

    pad.cd()

    lumiText = ""
    if( iPeriod==1 ):
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif ( iPeriod==2 ):
        lumiText += lumi_8TeV
        lumiText += " (8 TeV)"

    elif( iPeriod==3 ):      
        lumiText = lumi_8TeV 
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif ( iPeriod==4 ):
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
    elif ( iPeriod==7 ):
        if( outOfFrame ):lumiText += "#scale[0.85]{"
        lumiText += lumi_13TeV 
        lumiText += " (13 TeV)"
        lumiText += " + "
        lumiText += lumi_8TeV 
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
        if( outOfFrame): lumiText += "}"
    elif ( iPeriod==12 ):
        lumiText += "8 TeV"
    elif ( iPeriod==0 ):
        lumiText += lumi_sqrtS
            
    print lumiText

    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)    
    
    extraTextSize = extraOverCmsTextSize*cmsTextSize
    
    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)    

    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)

    if( outOfFrame ):
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11) 
        latex.SetTextSize(cmsTextSize*t)    
        latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText)
  
    pad.cd()

    posX_ = 0
    if( iPosX%10<=1 ):
        posX_ =   l + relPosX*(1-l-r)
    elif( iPosX%10==2 ):
        posX_ =  l + 0.5*(1-l-r)
    elif( iPosX%10==3 ):
        posX_ =  1-r - relPosX*(1-l-r)

    posY_ = 1-t - relPosY*(1-t-b)

    if( not outOfFrame ):
        if( drawLogo ):
            posX_ =   l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            CMS_logo = rt.TASImage("CMS-BW-label.png")
            pad_logo =  rt.TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 )
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()          
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if( writeExtraText ) :
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize*t)
                latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText)
    elif( writeExtraText ):
        if( iPosX==0):
            posX_ =   l +  relPosX*(1-l-r)
            posY_ =   1-t+lumiTextOffset*t

        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)      

    pad.Update()

########################################
###### Make Asymptotic Limit Plot ######
########################################

def getAsymLimits(file):
        
 f = ROOT.TFile(file);
 t = f.Get("limit"); 
 entries = t.GetEntries();
 lims = [0,0,0,0,0,0];

 limit_entries = 0;  
 print "\n\n------------------ MATTEO CHECK ----------------------\n\n"   
 for i in range(entries):
        
  t.GetEntry(i);
  t_quantileExpected = t.quantileExpected;
  t_limit = t.limit;
        
#  if t_quantileExpected == 0.5:  print "entry: ",i," limit: ", t_limit, ", quantileExpected: ",t_quantileExpected;        

  if t_quantileExpected == -1.: lims[0] += t_limit; 
  elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] += t_limit;
  elif t_quantileExpected >= 0.15 and t_quantileExpected <= 0.17: lims[2] += t_limit;            
  elif t_quantileExpected == 0.5: lims[3] += t_limit; limit_entries += 1 ; print "limit_entries modificato:\t contatore: %f \t %f"%(i,limit_entries)           
  elif t_quantileExpected >= 0.83 and t_quantileExpected <= 0.85: lims[4] += t_limit;
  elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] += t_limit;
  else: print "Unknown quantile!"
  print "Counter: %f"%i

# print lims[0]," / ",limit_entries, " = ",lims[0]/limit_entries;
 ## MATTEO CHANGED
 #if limit_entries==0:
    #limit_entries=1.;
 
  #### MATTEO CHANGED
 
 print "Entries:\t%f\n"%entries
 print "Range entries: \t"
 print range(entries) 
 print "\nlimit Entries: \t"
 print limit_entries
 print "\n\n------------------------------------------------------\n\n\n"
 lims[0] = lims[0]/limit_entries ;
 lims[1] = lims[1]/limit_entries ;
 lims[2] = lims[2]/limit_entries ;
 lims[3] = lims[3]/limit_entries ;
 lims[4] = lims[4]/limit_entries ;
 lims[5] = lims[5]/limit_entries ;
    
 return lims;

########################################
###### Submit batch job for cards ######
########################################

def submitBatchJob( command, fn ):
    
 currentDir = os.getcwd();
 CMSSWDir = currentDir+"/../";
    
 # create a dummy bash/csh
 outScript = open(fn+".sh","w");
 
 if not options.lxbatchCern and not options.herculesMilano :
  outScript.write('#!/bin/bash');
  outScript.write("\n"+'date');
  outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
  outScript.write("\n"+'echo "condor dir: " ${_CONDOR_SCRATCH_DIR}');    
  outScript.write("\n"+'cd '+currentDir);
  outScript.write("\n"+'eval `scram runtime -sh`');
  outScript.write("\n"+'cd -');    
  outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
  outScript.write("\n"+'echo ${PATH}');
  outScript.write("\n"+'ls'); 
  outScript.write("\n"+command);  
  outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *');    
  outScript.close();

  condorScript = open("condor_"+fn,"w");
  condorScript.write('universe = vanilla')
  condorScript.write("\n"+"Executable = "+fn+".sh")
  condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
  condorScript.write("\n"+'Should_Transfer_Files = YES')
  condorScript.write("\n"+'Transfer_Input_Files = doFit_class_higgs.py, BiasStudy/do_fitBias_higgs.py, AutoDict_std__map_std__string_std__string__cxx.so')    
  condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
  condorScript.write("\n"+'Output = out_$(Cluster).stdout')
  condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
  condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
  condorScript.write("\n"+'Log    = out_$(Cluster).log')
  condorScript.write("\n"+'Notification    = Error')
  condorScript.write("\n"+'Queue 1')
  condorScript.close();

  os.system("condor_submit "+"condor_"+fn);

 elif options.lxbatchCern and not options.herculesMilano:
  outScript.write('#!/bin/bash');
  outScript.write("\n"+'cd '+CMSSWDir);
  outScript.write("\n"+'eval `scram runtime -sh`');
  outScript.write("\n"+'cd -');
  outScript.write("\necho $PWD");
  outScript.write("\nls");
  outScript.write("\ncp "+currentDir+"/"+options.datacardDIR+"/ww* ./")
#  outScript.write("\n"+'cd '+currentDir);
#  outScript.write("\n"+'eval `scram runtime -sh`');
#  outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
#  outScript.write("\n"+'echo ${PATH}');
#  outScript.write("\n"+'ls');  
  outScript.write("\n"+command);
  outScript.write("\n"+'rm *.out');  
  outScript.close();
         
  os.system("chmod 777 "+currentDir+"/"+fn+".sh");
  if options.queque != "":
   os.system("bsub -q "+options.queque+" -cwd "+currentDir+" "+fn+".sh");
  else:
   os.system("bsub -q cmscaf1nd -cwd "+currentDir+" "+fn+".sh");
      
 elif not options.lxbatchCern and options.herculesMilano:

  outScript.write('#!/bin/bash');
  outScript.write("\n"+'cd '+currentDir);
  outScript.write("\n"+'eval `scram runtime -sh`');
  outScript.write("\n"+'cd -');
  outScript.write("\n"+'cp '+currentDir+'/BiasStudy/do_fitBias_higgs.py ./');
  outScript.write("\n"+'cp '+currentDir+'/doFit_class_higgs.py ./');
  outScript.write("\n"+'ls');  
  outScript.write("\n"+"unbuffer "+command+" > /gwteray/users/brianza/output"+fn+".txt");
  outScript.close();
         
  os.system("chmod 777 "+currentDir+"/"+fn+".sh");
  
  if options.queque != "":
   os.system("qsub -V -d "+currentDir+" -q "+options.queque+" "+currentDir+"/"+fn+".sh");
  else:
   os.system("qsub -V -d "+currentDir+" -q longcms "+currentDir+"/"+fn+".sh");
      

##########################################
###### Submit batch job for combine ######
##########################################

def submitBatchJobCombine( command, fn, mass, cprime, BRnew ):
    
    currentDir = os.getcwd();
    CMSSWDir = currentDir+"/../";
    
    if options.sigChannel !="": 
     SIGCH = options.jetBin+"_"+options.sigChannel;
    else:
     SIGCH = options.jetBin;
#    currentDir = os.getcwd();
    
    # create a dummy bash/csh
    outScript = open(fn+".sh","w");

    file1 = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP_unbin.txt"%(mass,SIGCH);
    file2 = "wwlvj_BulkGraviton_newxsec%03d_mu%s_HP_workspace.root"%(mass,SIGCH);
    file3 = "wwlvj_BulkGraviton_newxsec%03d_el%s_HP_workspace.root"%(mass,SIGCH);

    if not options.lxbatchCern and not options.herculesMilano :   
     outScript.write('#!/bin/bash');
     outScript.write("\n"+'date');
     outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
     outScript.write("\n"+'cd '+currentDir);
     outScript.write("\n"+'eval `scram runtime -sh`');
     outScript.write("\n"+'cd -');
     outScript.write("\n"+'ls');    
     outScript.write("\n"+command);
     if options.inputGeneratedDataset != "" and options.generateOnly:
      outScript.write("\n "+"mv higgsCombine* "+currentDir+"/"+options.inputGeneratedDataset);
      
     outScript.write("\n "+"rm rootstats* ");
     outScript.close();
    
     # link a condor script to your shell script
     condorScript = open("subCondor_"+fn,"w");
     condorScript.write('universe = vanilla')
     condorScript.write("\n"+"Executable = "+fn+".sh")
     condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
     condorScript.write("\n"+'Should_Transfer_Files = YES')
     condorScript.write("\n"+'Transfer_Input_Files = '+file1+', '+file2+', '+file3)    
     condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
     condorScript.write("\n"+'Output = out_$(Cluster).stdout')
     condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
     condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
     condorScript.write("\n"+'Log    = out_$(Cluster).log')
     condorScript.write("\n"+'Notification    = Error')
     condorScript.write("\n"+'Queue 1')
     condorScript.close();

     # submit the condor job     
     os.system("condor_submit "+"subCondor_"+fn)
 
    elif options.lxbatchCern and not options.herculesMilano: 

     outScript.write('#!/bin/bash');
     outScript.write("\n"+'cd '+CMSSWDir);
     outScript.write("\n"+'eval `scram runtime -sh`');
     outScript.write("\n"+'cd -');
     outScript.write("\necho $PWD");
     outScript.write("\nls");
#     outScript.write('#!/bin/bash');
#     outScript.write("\n"+'date');
#     outScript.write("\n"+'cd '+currentDir);
#     outScript.write("\n"+'eval `scram runtime -sh`');
#     outScript.write("\n"+'cd -');
#     outScript.write("\n"+'ls');    
    
     file1 = "wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt"%(mass,options.channel,SIGCH);
     file2 = "wwlvj_BulkGraviton_newxsec%03d_mu%s_HP_workspace.root"%(mass,SIGCH);
     file3 = "wwlvj_BulkGraviton_newxsec%03d_el%s_HP_workspace.root"%(mass,SIGCH);
     file4 = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP_workspace.root"%(mass,SIGCH);

     outScript.write("\ncp "+currentDir+"/"+file1+" ./")
     outScript.write("\ncp "+currentDir+"/"+file2+" ./")
     outScript.write("\ncp "+currentDir+"/"+file3+" ./")
     outScript.write("\ncp "+currentDir+"/"+file4+" ./")
     outScript.write("\ncp -r "+currentDir+"/"+options.inputGeneratedDataset+" ./")
     outScript.write("\ncp "+currentDir+"/"+"list* ./")

#     outScript.write("\n"+'cp '+currentDir+'/'+file1+' ./');
#     outScript.write("\n"+'cp '+currentDir+'/'+file2+' ./');
#     outScript.write("\n"+'cp '+currentDir+'/'+file3+' ./');
#     outScript.write("\n"+'cp '+currentDir+'/'+file4+' ./');
     outScript.write("\n"+command);
     if options.inputGeneratedDataset != "" and options.generateOnly:
#      outScript.write("\n "+"cp higgsCombine* "+currentDir+"/"+options.inputGeneratedDataset);
      outScript.write("\n "+"cp higgsCombine* "+currentDir+"/"+options.datacardDIR);
     outScript.write("\n "+"cp higgsCombine* "+currentDir+"/"+options.inputGeneratedDataset);
     outScript.write("\n "+"cp mlfit* "+currentDir+"/"+options.inputGeneratedDataset);
     outScript.write("\n "+"rm rootstats* ");
     outScript.close();

     os.system("chmod 777 "+currentDir+"/"+fn+".sh");
     if options.queque!="" :
      os.system("bsub -q "+options.queque+" -cwd "+currentDir+" "+fn+".sh");
     else: 
      os.system("bsub -q cmscaf1nd -cwd "+currentDir+" "+fn+".sh");

    elif not options.lxbatchCern and options.herculesMilano: 

     outScript.write('#!/bin/bash');
     outScript.write("\n"+'date');
     outScript.write("\n"+'cd '+currentDir);
     outScript.write("\n"+'eval `scram runtime -sh`');
     outScript.write("\n"+'cd -');
     outScript.write("\n"+'ls');    
    
     file1 = "wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt"%(mass,options.channel,SIGCH);
     file2 = "wwlvj_BulkGraviton_newxsec%03d_mu%s_HP_workspace.root"%(mass,SIGCH);
     file3 = "wwlvj_BulkGraviton_newxsec%03d_el%s_HP_workspace.root"%(mass,SIGCH);
     file4 = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP_workspace.root"%(mass,SIGCH);

     outScript.write("\n"+'cp '+currentDir+'/'+file1+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file2+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file3+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file4+' ./');
     outScript.write("\n"+command);
     if options.inputGeneratedDataset != "" and options.generateOnly:
      outScript.write("\n "+"mv higgsCombine* "+currentDir+"/"+options.inputGeneratedDataset);
     outScript.write("\n "+"rm rootstats* ");
     outScript.close();

     os.system("chmod 777 "+currentDir+"/"+fn+".sh");

     if options.queque != "":
      os.system("qsub -V -d "+currentDir+" -q "+options.queque+" "+currentDir+"/"+fn+".sh");
     else:
      os.system("qsub -V -d "+currentDir+" -q longcms "+currentDir+"/"+fn+".sh");
         

#####################################
###### definition of the style ######
#####################################

def setStyle():

  gStyle.SetPadBorderMode(0);
  gStyle.SetFrameBorderMode(0);
  gStyle.SetPadBottomMargin(0.12);
  gStyle.SetPadLeftMargin(0.12);
  gStyle.SetCanvasColor(ROOT.kWhite);
  gStyle.SetCanvasDefH(600); #Height of canvas                                                                                                                                            
  gStyle.SetCanvasDefW(600); #Width of canvas                                                                                                                                             
  gStyle.SetCanvasDefX(0);   #POsition on screen
  gStyle.SetCanvasDefY(0);
  gStyle.SetPadTopMargin(0.05);
  gStyle.SetPadBottomMargin(0.15);#0.13);
  gStyle.SetPadLeftMargin(0.15);#0.16);
  gStyle.SetPadRightMargin(0.05);#0.02);                                                                                                                                                  
  # For the Pad:                                                                                                                                                                        
  gStyle.SetPadBorderMode(0);
  gStyle.SetPadColor(ROOT.kWhite);
  gStyle.SetPadGridX(ROOT.kFALSE);
  gStyle.SetPadGridY(ROOT.kFALSE);
  gStyle.SetGridColor(0);
  gStyle.SetGridStyle(3);
  gStyle.SetGridWidth(1);

  # For the frame:                                                                                                                                                                       
  gStyle.SetFrameBorderMode(0);
  gStyle.SetFrameBorderSize(1);
  gStyle.SetFrameFillColor(0);
  gStyle.SetFrameFillStyle(0);
  gStyle.SetFrameLineColor(1);
  gStyle.SetFrameLineStyle(1);
  gStyle.SetFrameLineWidth(1);

  gStyle.SetAxisColor(1, "XYZ");
  gStyle.SetStripDecimals(ROOT.kTRUE);
  gStyle.SetTickLength(0.03, "XYZ");
  gStyle.SetNdivisions(505, "XYZ");
  gStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame                                                                                                        
  gStyle.SetPadTickY(1);
  gStyle.SetGridColor(0);
  gStyle.SetGridStyle(3);
  gStyle.SetGridWidth(1);

  gStyle.SetTitleColor(1, "XYZ");
  gStyle.SetTitleFont(42, "XYZ");
  gStyle.SetTitleSize(0.05, "XYZ");
  gStyle.SetTitleXOffset(1.15);#0.9);                                                                                                                                                  
  gStyle.SetTitleYOffset(1.3); # => 1.15 if exponents                                                                                                                                   
  gStyle.SetLabelColor(1, "XYZ");
  gStyle.SetLabelFont(42, "XYZ");
  gStyle.SetLabelOffset(0.007, "XYZ");
  gStyle.SetLabelSize(0.045, "XYZ");

  gStyle.SetPadBorderMode(0);
  gStyle.SetFrameBorderMode(0);
  gStyle.SetTitleTextColor(1);
  gStyle.SetTitleFillColor(10);
  gStyle.SetTitleFontSize(0.05);

  gStyle.SetOptStat(0);
  gStyle.SetOptTitle(0)
  gStyle.SetOptFit(1)

  NRGBs = 5
  NCont = 255
  stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
  red   = [ 0.00, 0.00, 0.87, 1.00, 0.51 ]
  green = [ 0.00, 0.81, 1.00, 0.20, 0.00 ]
  blue  = [ 0.51, 1.00, 0.12, 0.00, 0.00 ]
  stopsArray = array('d', stops)
  redArray   = array('d', red)
  greenArray = array('d', green)
  blueArray  = array('d', blue)
  TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont)
  gStyle.SetNumberContours(NCont)

####################################################
### Get PValue from combine -M ProfileLikelihood ###
####################################################

def getPValueFromCard(file,observed):

 f = ROOT.TFile(file);
 t = f.Get("limit");
 entries = t.GetEntries();

 lims = 0;

 for i in range(entries):
  t.GetEntry(i);
  lims += t.limit ;

 return lims/entries;

##############################
#### Make SM Limits Plots ####  
##############################  

def makeSMLimitPlot(SIGCH,cprime = 10, brnew = 00):

    nPoints = len(mass);
    xsec = xsec_value; ### MATTEO CHANGED #,0.165094354,0.110500282,0.076058293,0.0382161,0.020499693];
    xsec_corr = xsec_corr_value; ### MATTEO CHANGED #,1.,1.,1.,1.,1.];

    xbins     = array('f', []); xbins_env = array('f', []);
    ybins_exp = array('f', []); ybins_obs = array('f', []);
    ybins_1s  = array('f', []); ybins_2s  = array('f', []);

    setStyle();
     
    for i in range(len(mass)):
	curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
        print curFile
	curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
	xbins_env.append( mass[i] );
        ybins_exp.append( curAsymLimits[3] );
        ybins_obs.append( curAsymLimits[0] );
        ybins_2s.append( curAsymLimits[1] );
        ybins_1s.append( curAsymLimits[2] );
        
    for i in range( len(mass)-1, -1, -1 ):
	curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
	curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5] );
        ybins_1s.append( curAsymLimits[4] );
                                                                    
    
    curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
    curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
    curGraph_1s  = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_1s);
    curGraph_2s  = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_2s);


    curGraph_obs.SetMarkerStyle(20);
    curGraph_obs.SetLineWidth(3);
    curGraph_obs.SetLineStyle(1);
    curGraph_obs.SetMarkerSize(1.6);
    curGraph_exp.SetMarkerSize(1.3);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_exp.SetLineStyle(2);
    curGraph_exp.SetLineWidth(3);
    curGraph_exp.SetMarkerSize(2);
    curGraph_exp.SetMarkerStyle(24);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_1s.SetFillColor(ROOT.kGreen);
    curGraph_1s.SetFillStyle(1001);
    curGraph_1s.SetLineStyle(ROOT.kDashed);
    curGraph_1s.SetLineWidth(3);

    curGraph_2s.SetFillColor(ROOT.kYellow);
    curGraph_2s.SetFillStyle(1001);
    curGraph_2s.SetLineStyle(ROOT.kDashed);
    curGraph_2s.SetLineWidth(3);
                               
    oneLine = ROOT.TF1("oneLine","1",599,1001);
    oneLine.SetLineColor(ROOT.kRed);
    oneLine.SetLineWidth(3);

    setStyle();
    
    can_SM = ROOT.TCanvas("can_SM","can_SM",600,650); 

    hrl_SM = can_SM.DrawFrame(599,0.01,1001,1000);#ROOT.TMath.MaxElement(curGraph_2s.GetN(),curGraph_2s.GetY())*1.2);
    hrl_SM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{theory}");
    hrl_SM.GetYaxis().SetTitleOffset(1.35);
    hrl_SM.GetYaxis().SetTitleSize(0.045);
    hrl_SM.GetYaxis().SetTitleFont(42);

    hrl_SM.GetXaxis().SetTitle("M_{G} (GeV)");
    hrl_SM.GetXaxis().SetTitleSize(0.045);
    hrl_SM.GetXaxis().SetTitleFont(42);

#    hrl_SM.GetYaxis().SetNdivisions(505);
    can_SM.SetGridx(1);
    can_SM.SetGridy(1);
    ROOT.gPad.SetLogy();
                   
    curGraph_2s.Draw("F");
    curGraph_1s.Draw("Fsame");
#    if not options.blindObservedLine : curGraph_obs.Draw("PCsame");
    curGraph_obs.Draw("PCsame");
    curGraph_exp.Draw("Csame");
    oneLine.Draw("same");

    leg2 = ROOT.TLegend(0.3,0.72,0.75,0.9);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);

    leg2.AddEntry(curGraph_exp,"Asympt. CL_{S} Expected","L")
    leg2.AddEntry(curGraph_1s, "Asympt. CL_{S} Expected #pm 1#sigma","LF")
    leg2.AddEntry(curGraph_2s, "Asympt. CL_{S} Expected #pm 2#sigma","LF")
                                       
#    if not options.blindObservedLine:     leg2.AddEntry(curGraph_obs,"Asympt. CL_{S} Observed","LP")
    leg2.AddEntry(curGraph_obs,"Asympt. CL_{S} Observed","LP")

    can_SM.Update();
    can_SM.RedrawAxis();
    can_SM.RedrawAxis("g");
    can_SM.Update();

    leg2.Draw();

    banner = TPaveText( 0.145, 0.953, 0.56, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.03);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.5,0.953,0.96,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.03);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("W #rightarrow l#nu, L = 2.3 fb^{-1} at #sqrt{s} = 13 TeV");
    label_sqrt.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    
    can_SM.SaveAs("limitFigs/SMLim_%s_HP.png"%(options.channel));
    can_SM.SaveAs("limitFigs/SMLim_%s_HP.pdf"%(options.channel));

#############################
#### Make limit in xsec #####
#############################
def makeSMXsecPlot(SIGCH,cprime = 10, brnew = 00):

    nPoints = len(mass);
    xsec = xsec_value;### MATTEO CHANGED #,0.165094354,0.110500282,0.076058293,0.0382161,0.020499693];
    xsec_corr = xsec_corr_value; ### MATTEO CHANGED#,1,1,1,1,1 ];
    
    xbins     = array('f', []); xbins_env = array('f', []);
    ybins_exp = array('f', []); ybins_obs = array('f', []);
    ybins_1s  = array('f', []); ybins_2s  = array('f', []);
    ybins_xsec= array('f', []);

    setStyle();
     
    for i in range(len(mass)):
	curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
        print curFile
	curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
	xbins_env.append( mass[i] );

        ybins_exp.append( curAsymLimits[3]*xsec[i]*xsec_corr[i] );
        print "expected: ",curAsymLimits[3]*xsec[i]*xsec_corr[i];
        ybins_obs.append( curAsymLimits[0]*xsec[i]*xsec_corr[i] );
        print "observed: ",curAsymLimits[0]*xsec[i]*xsec_corr[i];
        ybins_2s.append( curAsymLimits[1]*xsec[i]*xsec_corr[i] );
        ybins_1s.append( curAsymLimits[2]*xsec[i]*xsec_corr[i] );
#        ybins_xsec.append(25.*xsec[i]/0.43924356);
        ybins_xsec.append(xsec[i]);

    for i in range( len(mass)-1, -1, -1 ):
	curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
	curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5]*xsec[i]*xsec_corr[i]  );
        ybins_1s.append( curAsymLimits[4]*xsec[i]*xsec_corr[i]  );
                                                                    
    
    curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
    curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
    curGraph_1s  = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_1s);
    curGraph_2s  = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_2s);
    curGraph_xsec= ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_xsec);


    curGraph_obs.SetMarkerStyle(20);
    curGraph_obs.SetLineWidth(3);
    curGraph_obs.SetLineStyle(1);
    curGraph_obs.SetMarkerSize(1.6);
    curGraph_exp.SetMarkerSize(1.3);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_exp.SetLineStyle(2);
    curGraph_exp.SetLineWidth(3);
    curGraph_exp.SetMarkerSize(2);
    curGraph_exp.SetMarkerStyle(24);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_1s.SetFillColor(ROOT.kGreen);
    curGraph_1s.SetFillStyle(1001);
    curGraph_1s.SetLineStyle(ROOT.kDashed);
    curGraph_1s.SetLineWidth(3);

    curGraph_2s.SetFillColor(ROOT.kYellow);
    curGraph_2s.SetFillStyle(1001);
    curGraph_2s.SetLineStyle(ROOT.kDashed);
    curGraph_2s.SetLineWidth(3);

#    curGraph_obs.SetMarkerStyle(20);
    curGraph_xsec.SetLineWidth(3);
    curGraph_xsec.SetLineStyle(1);
#    curGraph_xsec.SetMarkerSize(1.6);
    curGraph_xsec.SetLineColor(ROOT.kRed);
                               
    oneG = ROOT.TGraph(4);
    print "length: ", len(xsec)
    for i in range(len(xsec)):
        oneG.SetPoint(i, mass[i], xsec[i]);
    oneLine = ROOT.TSpline3("oneLine",oneG);    
    oneLine.SetLineColor(ROOT.kRed);
    oneLine.SetLineWidth(3);

    setStyle();
    
    can_SM = ROOT.TCanvas("can_SM","can_SM",600,650); 
#    can_SM.cd();

    hrl_SM = can_SM.DrawFrame(599,0.001,1001,100);#ROOT.TMath.MaxElement(curGraph_2s.GetN(),curGraph_2s.GetY())*1.2);
#    upperPad    = ROOT.TPad("upperPad", "upperPad", .005, .005, .995, .950);
#    upperPad.SetLeftMargin(0.18);
#    upperPad.Draw();
#    upperPad.cd();

    hrl_SM.GetYaxis().SetTitle("#sigma_{95%} x BR(G_{Bulk} #rightarrow WW)(pb)");
    hrl_SM.GetYaxis().SetTitleOffset(1.35);
    hrl_SM.GetYaxis().SetTitleSize(0.045);
    hrl_SM.GetYaxis().SetTitleFont(42);

    hrl_SM.GetXaxis().SetTitle("M_{G} (GeV)");
    hrl_SM.GetXaxis().SetTitleSize(0.045);
    hrl_SM.GetXaxis().SetTitleFont(42);

#    hrl_SM.GetYaxis().SetNdivisions(505);
    can_SM.SetGridx(1);
    can_SM.SetGridy(1);
    ROOT.gPad.SetLogy();
                   
    curGraph_2s.Draw("F");
    curGraph_1s.Draw("Fsame");
    if not options.blindObservedLine : curGraph_obs.Draw("PCsame");
    curGraph_exp.Draw("Csame");
    curGraph_xsec.Draw("Csame");
#    oneLine.Draw("same");

#    leg2 = ROOT.TLegend(0.3,0.72,0.75,0.9);
    leg2 = ROOT.TLegend(0.4,0.72,0.9,0.9);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);

    leg2.AddEntry(curGraph_exp,"Asympt. CL_{S} Expected","L")
#    leg2.AddEntry(curGraph_1s, "Asympt. CL_{S} Expected #pm 1#sigma","LF")
#    leg2.AddEntry(curGraph_2s, "Asympt. CL_{S} Expected #pm 2#sigma","LF")
    leg2.AddEntry(curGraph_1s, "Asympt. CL_{S} Expected #pm 1 s.d.","LF")
    leg2.AddEntry(curGraph_2s, "Asympt. CL_{S} Expected #pm 2 s.d.","LF")
    leg2.AddEntry(curGraph_xsec, "#sigma_{TH} #times BR_{G_{Bulk}#rightarrow WW}, k = 0.5","L")
                                       
    if not options.blindObservedLine:     leg2.AddEntry(curGraph_obs,"Asympt. CL_{S} Observed","LP")

    CMS_lumi(can_SM,4,11);

    can_SM.Update();
    can_SM.RedrawAxis();
    can_SM.RedrawAxis("g");
    can_SM.Update();

    leg2.Draw();

    banner = TPaveText( 0.145, 0.953, 0.56, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.03);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
#    banner.AddText(leftText);
#    banner.Draw();

    label_sqrt = TPaveText(0.23,0.253,0.33,0.275, "brNDC");
#    label_sqrt = TPaveText(0.2,0.803,0.3,0.825, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.03);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                   

    if (options.channel=="el"):
        label_sqrt.AddText("W #rightarrow e#nu");
    elif (options.channel=="mu"):
        label_sqrt.AddText("W #rightarrow #mu#nu");
    else:
        label_sqrt.AddText("W #rightarrow l#nu");

    label_sqrt.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    
    can_SM.SaveAs("limitFigs/SMXsec_%s_HP.png"%(options.channel));
    can_SM.SaveAs("limitFigs/SMXsec_%s_HP.pdf"%(options.channel));

##############################
#### Make SM PValue Plots ####  
##############################  

def makeSMPValuePlot(SIGCH,cprime = 10,brnew = 00):

    nPoints = len(mass);
    
    xbins_obs     = array('f', []);
    xbins_exp     = array('f', []);

    ybins_obs     = array('f', []);
    ybins_exp     = array('f', []);

    setStyle();

    for i in range(len(mass)):
	curFile_obs = "higgsCombinewwlvj_pval_obs_BulkGraviton_newxsec%03d_%s%s_HP_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
	curFile_exp = "higgsCombinewwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
        xbins_obs.append(mass[i]); 
        xbins_exp.append(mass[i]); 
        print "higgsCombinewwlvj_pval_obs_BulkGraviton_newxsec%03d_%s%s_HP_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
        ybins_obs.append(getPValueFromCard(curFile_obs,1));
        print getPValueFromCard(curFile_obs,1)
        print "higgsCombinewwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]); 
        ybins_exp.append(getPValueFromCard(curFile_exp,0));
        print getPValueFromCard(curFile_exp,0)

    gr_obs = ROOT.TGraphAsymmErrors(nPoints,xbins_obs,ybins_obs);
    gr_exp = ROOT.TGraphAsymmErrors(nPoints,xbins_exp,ybins_exp);
                    
    gr_obs.SetLineColor(1); gr_obs.SetMarkerColor(1); gr_obs.SetMarkerStyle(20); gr_obs.SetLineWidth(3);gr_obs.SetMarkerSize(1.6);
    gr_exp.SetLineColor(2); gr_exp.SetMarkerColor(2); gr_exp.SetMarkerStyle(20); gr_exp.SetLineWidth(3);gr_exp.SetMarkerSize(1.6);
    gr_exp.SetLineStyle(8);

    oneSLine = ROOT.TF1("oneSLine","1.58655253931457074e-01",mass[0],mass[len(mass)-1]);
    oneSLine.SetLineColor(ROOT.kRed); oneSLine.SetLineWidth(2); oneSLine.SetLineStyle(2);
    twoSLine = ROOT.TF1("twoSLine","2.27501319481792155e-02",mass[0],mass[len(mass)-1]);
    twoSLine.SetLineColor(ROOT.kRed); twoSLine.SetLineWidth(2); twoSLine.SetLineStyle(2);
    threeSLine = ROOT.TF1("threeSLine","1.34989803163009588e-03",mass[0],mass[len(mass)-1]);
    threeSLine.SetLineColor(ROOT.kRed); threeSLine.SetLineWidth(2); threeSLine.SetLineStyle(2);
    fourSLine = ROOT.TF1("fourSLine","3.16712418331199785e-05",mass[0],mass[len(mass)-1]);
    fourSLine.SetLineColor(ROOT.kRed); fourSLine.SetLineWidth(2); fourSLine.SetLineStyle(2);
    
    banner = TLatex(0.32,0.955,("CMS Preliminary, 1 fb^{-1} at #sqrt{s}=13 TeV"));
    banner.SetNDC(); banner.SetTextSize(0.035);

    ban1s = TLatex(950,1.58655253931457074e-01,("1 #sigma"));
    ban1s.SetTextSize(0.028); ban1s.SetTextColor(1)
    ban2s = TLatex(950,2.27501319481792155e-02,("2 #sigma"));
    ban2s.SetTextSize(0.028); ban2s.SetTextColor(1)
    ban3s = TLatex(950,1.34989803163009588e-03,("3 #sigma"));
    ban3s.SetTextSize(0.028); ban3s.SetTextColor(1);
    ban4s = TLatex(950,3.16712418331199785e-05,("4 #sigma"));
    ban4s.SetTextSize(0.028); ban4s.SetTextColor(1)

    leg2 = ROOT.TLegend(0.5,0.2,0.9,0.35);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(1);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    
    if options.jetBin == "":
     leg2.AddEntry( gr_obs, "obs signif, %s "%options.channel, "pl" );
     leg2.AddEntry( gr_exp, "exp signif, Bulk Graviton", "pl" );
    else:
     leg2.AddEntry( gr_obs, "obs signif, VBF e+mu", "pl" );
     leg2.AddEntry( gr_exp, "exp signif, SM Higgs", "pl" );

    zeroLine = ROOT.TF1("zeroLine","0.5",mass[0],mass[len(mass)-1]);
    zeroLine.SetLineColor(ROOT.kBlue); oneSLine.SetLineWidth(2);

    can = ROOT.TCanvas("can","can",600,650);
    hrl = can.DrawFrame(mass[0],1e-6,mass[len(mass)-1],1.);
    hrl.GetYaxis().SetTitle("p-value");
    hrl.GetXaxis().SetTitle("m_{X} (GeV)");
    can.SetGrid();
    ROOT.gPad.SetLogy();
    if options.blindObservedLine==1:
        gr_exp.Draw("PL");
    else:
        gr_obs.Draw("PL");
        gr_exp.Draw("PLsame");
    zeroLine.Draw("same");
    oneSLine.Draw("same");
    twoSLine.Draw("same");
    threeSLine.Draw("same");
    fourSLine.Draw("same");
    banner.Draw();
    leg2.Draw();
    ban1s.Draw();
    ban2s.Draw();
    ban3s.Draw();
    ban4s.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can.SaveAs("limitFigs/SMPvals_%s_HP.pdf"%(options.channel),"pdf");
    can.SaveAs("limitFigs/SMPvals_%s_HP.png"%(options.channel),"png");

####################################
#### Make Signal Strenght Plots ####  
####################################  

def makeSignalStrenghtPlot(SIGCH,cprime = 10,brnew = 00):

    nPoints = len(mass);

    xbins_mu         = array('f', []);
    xbins_mu_err_up  = array('f', []);
    xbins_mu_err_dn  = array('f', []);
    ybins_mu         = array('f', []);
    ybins_mu_err_up  = array('f', []);
    ybins_mu_err_dn  = array('f', []);
    ybins_mu_err_up_2s  = array('f', []);
    ybins_mu_err_dn_2s  = array('f', []);
    ybins_mu_ratio   = array('f', []);
    xbins_obs     = array('f', [0.]);
    ybins_obs     = array('f', [0.]);

    setStyle();
        
    for i in range(len(mass)):
        
     curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.MaxLikelihoodFit.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
 
     f = ROOT.TFile(curFile);
     t = f.Get("limit");
     entries = t.GetEntries();

     mu = 0;
     effect_entry = 0 ;
     mu_err_up = 0 ; mu_err_dn = 0 ;
     mu_err_up_2s = 0 ; mu_err_dn_2s = 0 ;
     
     mu_temp = 0 ;  
     for ientry in range(t.GetEntries()):
       t.GetEntry(ientry);
       if t.quantileExpected == 0.5 : mu += t.limit ;  effect_entry += 1; mu_temp = t.limit ; continue ; 
       if t.quantileExpected > 0.15 and t.quantileExpected < 0.17 : mu_err_dn += abs(t.limit-mu_temp) ; continue ;
       if t.quantileExpected > 0.83 and t.quantileExpected < 0.85 : mu_err_up += abs(t.limit-mu_temp) ; continue ;
       if t.quantileExpected >= 0.024 and t.quantileExpected <= 0.026: mu_err_up_2s += abs(t.limit-mu_temp) ; continue ;
       if t.quantileExpected >= 0.974 and t.quantileExpected <= 0.976: mu_err_dn_2s += abs(t.limit-mu_temp) ; continue ;
                         
     mu = mu/effect_entry;
     mu_err_dn = mu_err_dn/effect_entry;
     mu_err_up = mu_err_up/effect_entry;
     mu_err_dn_2s = mu_err_dn_2s/effect_entry;
     mu_err_up_2s = mu_err_up_2s/effect_entry;
     
     xbins_mu.append(mass[i]); 
     xbins_mu_err_up.append(0.); 
     xbins_mu_err_dn.append(0.); 

     ybins_mu.append(mu);
     ybins_mu_err_up.append(mu_err_up);
     ybins_mu_err_dn.append(mu_err_dn);
     ybins_mu_err_up_2s.append(mu_err_up_2s);
     ybins_mu_err_dn_2s.append(mu_err_dn_2s);

     if mu >= 0 : ybins_mu_ratio.append(mu/mu_err_dn);
     else:        ybins_mu_ratio.append(mu/mu_err_up);

     if options.plotPValue == 1:
	curFile_obs = "higgsCombinewwlvj_pval_obs_BulkGraviton_newxsec%03d_%s%s_HP_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);
        xbins_obs.append(mass[i]); 
        ybins_obs.append(getPValueFromCard(curFile_obs,1));

                    
    gr_mu_1 = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu);
    gr_mu_1.SetLineColor(1); gr_mu_1.SetMarkerColor(1); gr_mu_1.SetMarkerStyle(20); gr_mu_1.SetLineWidth(3);gr_mu_1.SetMarkerSize(1.6);

    gr_mu_ratio = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu_ratio);
    gr_mu_ratio.SetLineColor(1); gr_mu_ratio.SetMarkerColor(1); gr_mu_ratio.SetMarkerStyle(20); gr_mu_ratio.SetLineWidth(3);gr_mu_ratio.SetMarkerSize(1.6);

    gr_mu_1_1s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn,ybins_mu_err_up);
    gr_mu_1_1s.SetLineColor(1); gr_mu_1_1s.SetMarkerColor(1); gr_mu_1_1s.SetMarkerStyle(20); gr_mu_1_1s.SetLineWidth(5); gr_mu_1_1s.SetMarkerSize(1.6);

    gr_mu_1_2s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn_2s,ybins_mu_err_up_2s);
    gr_mu_1_2s.SetLineColor(ROOT.kBlue); gr_mu_1_2s.SetMarkerColor(ROOT.kBlue); gr_mu_1_2s.SetMarkerStyle(20); gr_mu_1_2s.SetLineWidth(3); gr_mu_1_2s.SetMarkerSize(1.6);

    gr_mu_zero_1s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,xbins_mu_err_up,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn,ybins_mu_err_up);
    gr_mu_zero_2s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,xbins_mu_err_up,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn_2s,ybins_mu_err_up_2s);
    
    gr_mu_zero_1s.SetFillStyle(3001);
    gr_mu_zero_1s.SetFillColor(ROOT.kRed);
    gr_mu_zero_2s.SetFillStyle(3002);
    gr_mu_zero_2s.SetFillColor(210);

    gr_obs = ROOT.TGraphAsymmErrors(nPoints+1,xbins_obs,ybins_obs);                    
    gr_obs.SetLineColor(ROOT.kBlue); gr_obs.SetMarkerColor(ROOT.kBlue); gr_obs.SetMarkerStyle(28); gr_obs.SetLineWidth(3);gr_obs.SetMarkerSize(1.6);

    oneSLine0 = ROOT.TF1("oneSLine","0",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine0.SetLineColor(ROOT.kRed); oneSLine0.SetLineWidth(2);
    oneSLine1 = ROOT.TF1("oneSLine","1",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine1.SetLineColor(ROOT.kRed); oneSLine1.SetLineWidth(2); oneSLine1.SetLineStyle(2);
    oneSLine2 = ROOT.TF1("oneSLine","2",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine2.SetLineColor(ROOT.kRed); oneSLine2.SetLineWidth(2); oneSLine2.SetLineStyle(2);
    oneSLine3 = ROOT.TF1("oneSLine","3",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine3.SetLineColor(ROOT.kRed); oneSLine3.SetLineWidth(2); oneSLine3.SetLineStyle(2);
    
    banner = TLatex(0.32,0.955,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV"));
    banner.SetNDC(); banner.SetTextSize(0.035);

    ban1s = TLatex(950,1.,("#mu SM"));
    ban1s.SetTextSize(0.028); ban1s.SetTextColor(1)

    leg2 = ROOT.TLegend(0.25,0.2,0.6,0.35);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(1);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    
    if options.jetBin == "":
     leg2.AddEntry( gr_mu_1, "signal strenght %s "%options.channel, "pl" );
    else:
     leg2.AddEntry( gr_mu_1, "signal strenght, VBF e+#mu", "pl" );

    can = ROOT.TCanvas("can","can",600,650);
    hrl = can.DrawFrame(mass[0],-5,mass[len(mass)-1],5);
    hrl.GetYaxis().SetTitle("signal strenght");
    hrl.GetXaxis().SetTitle("m_{H} (GeV)");
    can.SetGrid();
    gr_mu_1.Draw("P");
    gr_mu_zero_2s.Draw("e3same");
    gr_mu_zero_1s.Draw("e3same"); 
    gr_mu_1.Draw("Psame");
    oneSLine0.Draw("same");
   
    banner.Draw();
    leg2.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can.SaveAs("limitFigs/signal_strenght_band_%s%s_HP.pdf"%(options.channel,SIGCH),"pdf");
    can.SaveAs("limitFigs/signal_strenght_band_%s%s_HP.png"%(options.channel,SIGCH),"png");

    can_2 = ROOT.TCanvas("can_2","can_2",600,650);
    hrl_2 = can_2.DrawFrame(mass[0]-50,-5,mass[len(mass)-1]+50,5);
    hrl_2.GetYaxis().SetTitle("signal strenght");
    hrl_2.GetXaxis().SetTitle("m_{H} (GeV)");
    can_2.SetGrid();
    gr_mu_1_2s.Draw("P");
    gr_mu_1_1s.Draw("P"); 
    oneSLine0.Draw("same");
    oneSLine1.Draw("same");
    oneSLine2.Draw("same");
    oneSLine3.Draw("same");
   
    banner.Draw();
    leg2.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can_2.SaveAs("limitFigs/signal_strenght_%s%s_HP.pdf"%(options.channel,SIGCH),"pdf");
    can_2.SaveAs("limitFigs/signal_strenght_%s%s_HP.png"%(options.channel,SIGCH),"png");


    can_ratio = ROOT.TCanvas("can_ratio","can_ratio",600,650);
    hrl_ratio = can_ratio.DrawFrame(mass[0]-50,-5,mass[len(mass)-1]+50,5);
    hrl_ratio.GetYaxis().SetTitle("#mu/#sigma_{#mu}");
    hrl_ratio.GetXaxis().SetTitle("m_{H} (GeV)");
    can_ratio.SetGrid();
    leg2.Clear()
    if options.jetBin == "":
     leg2.AddEntry( gr_mu_ratio, "#sigma_{#mu}/#mu %s "%options.channel, "pl" );
    else:
     leg2.AddEntry( gr_mu_ratio, "#sigma_{#mu}/#mu, VBF e+#mu", "pl" );
    gr_mu_ratio.Draw("P");
    if options.plotPValue == 1 : gr_obs.Draw("P"); leg2.AddEntry( gr_obs, "local significance", "pl" );
    
    oneSLine0.Draw("same");
    oneSLine1.Draw("same");
    oneSLine2.Draw("same");
    oneSLine3.Draw("same");
    banner.Draw();
    leg2.Draw();
    
    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can_ratio.SaveAs("limitFigs/signal_strenght_ratio_%s%s_HP.pdf"%(options.channel,SIGCH),"pdf");
    can_ratio.SaveAs("limitFigs/signal_strenght_ratio_%s%s_HP.png"%(options.channel,SIGCH),"png");

#################################
### make likelihood scan plot ###
#################################    

def makeLikelihoodScanPlot(SIGCH):

    nPoints = len(mass);
    
    setStyle();

    banner = TLatex(0.32,0.955,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV"));
    banner.SetNDC(); banner.SetTextSize(0.035);

    leg2 = ROOT.TLegend(0.55,0.65,0.8,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(1);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    
    for i in range(len(mass)):

	curFile = "higgsCombinewwlvj_LikelihoodScan_BulkGraviton_newxsec%03d_%s%s_HP_unbin.MultiDimFit.mH%02d.root"%(mass[i],options.channel,SIGCH,mass[i]);
        f = ROOT.TFile(curFile);
        t = f.Get("limit");
        
        xbins_mu = array('f',[]); 
        ybins_mu = array('f',[]); 

        for ientry in range(t.GetEntries()):
         if ientry == 0 : continue ;   
         t.GetEntry(ientry);   
         xbins_mu.append(t.r); 
         ybins_mu.append(2*t.deltaNLL); 

        gr_mu = ROOT.TGraph(t.GetEntries()-1,xbins_mu,ybins_mu);
        gr_mu.Sort();
                                    
        gr_mu.SetLineWidth(2);
        gr_mu.SetLineColor(2);
        gr_mu.SetMarkerStyle(20);

        if options.jetBin == "":
         leg2.AddEntry(gr_mu, "#mu Likelihood scan %s "%options.channel, "pl" );
        else:
         leg2.AddEntry(gr_mu, "#mu Likelihood scan, VBF e+#mu", "pl" );

        can = ROOT.TCanvas("can_%d"%(mass[i]),"can_%d"%(mass[i]),600,650);
        gr_mu.GetYaxis().SetTitle("-2#Delta ln(L)");
        gr_mu.GetXaxis().SetTitle("signal strenght");
        can.SetGrid();
        gr_mu.Draw("AC");
        banner.Draw();
        leg2.Draw();

        os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
        can.SaveAs("limitFigs/LikelihoodScan_%03d_%s%s_HP.pdf"%(mass[i],options.channel,SIGCH),"pdf");
        can.SaveAs("limitFigs/LikelihoodScan_%03d_%s%s_HP.png"%(mass[i],options.channel,SIGCH),"png");
        leg2.Clear();
    
                                           
##############################################
### Make the BSM Limit plot vs mass and c' ###
##############################################    
def makeBSMLimitPlotMass(SIGCH,brnew):

    print "module ===> makeBSMLimits_vsMass";
    curcolors = [1,2,210,4,6,12,7];
    nPoints = len(mass);

    massCS  = [];
    if SIGCH == "" or SIGCH == "_2jet":
        massCS.append((0.5230+0.09688));
        massCS.append((0.2290+0.06330));
        massCS.append((0.1097+0.04365));
        massCS.append((0.0571+0.03164));
        massCS.append((0.0320+0.02399));
    elif SIGCH == "_BulkGraviton_newxsec":
        massCS.append((0.5230));
        massCS.append((0.2290));
        massCS.append((0.1097));
        massCS.append((0.0571));
        massCS.append((0.0320));
    elif SIGCH == "_vbfH":
        massCS.append((0.09688));
        massCS.append((0.06330));
        massCS.append((0.04365));
        massCS.append((0.03164));
        massCS.append((0.02399));
    else:
        print "problem!"
        
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];

    gridMaxExp    = -999;
    gridMaxSigExp = -999;
    gridMaxObs    = -999;
    gridMaxSigObs = -999;

    tGraphs_exp = [];
    tGraphs_obs = [];
    tGraphs_csXbr_exp = [];
    tGraphs_csXbr_obs = [];
    tGraphs_csXbr_th = [];
    
    for j in range(len(cprime)):

        if (cprime[j]*0.1) > (1-brnew*0.1): continue ;
            
        xbins           = array('f', []); ybins_exp       = array('f', []);
        ybins_obs       = array('f', []); ybins_csXbr_exp = array('f', []);
        ybins_csXbr_obs = array('f', []); ybins_csXbr_th  = array('f', []);

        for i in range(len(mass)):

            curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,mass[i]);

            curAsymLimits = getAsymLimits(curFile);
            xbins.append( mass[i] );
            ybins_exp.append( curAsymLimits[3] );
            ybins_obs.append( curAsymLimits[0] );
            ybins_csXbr_exp.append(curAsymLimits[3]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );
            ybins_csXbr_obs.append(curAsymLimits[0]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );
            ybins_csXbr_th.append(1.*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );

            if gridMaxExp < curAsymLimits[3]:
                gridMaxExp = curAsymLimits[3];

            if gridMaxObs < curAsymLimits[0]:
                gridMaxObs = curAsymLimits[0];

            cscur = ( curAsymLimits[3]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );

            if gridMaxSigExp  < ( curAsymLimits[3]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] ):
                gridMaxSigExp = ( curAsymLimits[3]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );

            if gridMaxSigObs  < ( curAsymLimits[0]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] ):
                gridMaxSigObs = ( curAsymLimits[0]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );

        curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
        curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
        curGraph_obs.SetMarkerStyle(20);
        curGraph_obs.SetLineWidth(3);
        curGraph_obs.SetLineStyle(1);
        curGraph_obs.SetMarkerSize(1);
        curGraph_exp.SetMarkerSize(1.3);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(3);
        curGraph_exp.SetMarkerSize(1);
        curGraph_exp.SetMarkerStyle(24);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);

        curGraph_csXbr_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_exp);
        curGraph_csXbr_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_obs);
        curGraph_csXbr_th  = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_th);
        curGraph_csXbr_exp.SetLineStyle(10);
        curGraph_csXbr_exp.SetLineWidth(2);
        curGraph_csXbr_obs.SetLineWidth(2);
        curGraph_csXbr_exp.SetMarkerSize(2);
        curGraph_csXbr_obs.SetMarkerSize(2);
        curGraph_csXbr_th.SetLineWidth(2);

        tGraphs_exp.append( curGraph_exp );
        tGraphs_obs.append( curGraph_obs );
        tGraphs_csXbr_exp.append( curGraph_csXbr_exp );
        tGraphs_csXbr_obs.append( curGraph_csXbr_obs );
        tGraphs_csXbr_th.append( curGraph_csXbr_th );


    setStyle();

    can_BSM = ROOT.TCanvas("can_BSM","can_BSM",630,600);
    hrl_BSM = can_BSM.DrawFrame(599,0.0,1001,max(gridMaxExp,gridMaxObs)*1.3);

    hrl_BSM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_BSM.GetYaxis().SetTitleOffset(1.35);
    hrl_BSM.GetYaxis().SetTitleSize(0.045);
    hrl_BSM.GetYaxis().SetTitleFont(42);

    hrl_BSM.GetXaxis().SetTitle("M_{H} (GeV)");
    hrl_BSM.GetXaxis().SetTitleSize(0.045);
    hrl_BSM.GetXaxis().SetTitleFont(42);

    hrl_BSM.GetYaxis().SetNdivisions(505);
    can_BSM.SetGridx(1);
    can_BSM.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    for k in range(len(tGraphs_exp)):
        tGraphs_exp[k].SetLineStyle(1);
        tGraphs_exp[k].SetLineColor(curcolors[k]);
        tGraphs_exp[k].SetMarkerColor(curcolors[k]);
        tGraphs_obs[k].SetLineColor(curcolors[k]);
        tGraphs_obs[k].SetMarkerColor(curcolors[k]);
        tGraphs_exp[k].SetLineWidth(2);
        tGraphs_obs[k].SetLineWidth(2);
        if not options.blindObservedLine:
            tGraphs_exp[k].SetLineStyle(7);
            tGraphs_exp[k].Draw("PL");            
            tGraphs_obs[k].Draw("PL");
        else:
            tGraphs_exp[k].Draw("PL");            
            
        tmplabel = "exp., C'^{ 2} = %1.1f"%(float((cprime[k])/10.));
        leg2.AddEntry(tGraphs_exp[k],tmplabel,"L");
        tmplabel = "obs., C'^{ 2} = %1.1f"%(float((cprime[k])/10.));
        if not options.blindObservedLine: leg2.AddEntry(tGraphs_obs[k],tmplabel,"L");


    banner2 = TLatex(0.17,0.91,("BR_{new} = 0"));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    can_BSM.Update();
    can_BSM.RedrawAxis();
    can_BSM.RedrawAxis("g");
    can_BSM.Update();
    
    banner = TPaveText( 0.145, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.036);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.5,0.953,0.96,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.036);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    label_sqrt.Draw();

    leg2.Draw();
    banner2.Draw();
    can_BSM.SaveAs("limitFigs/BSMLim%s_%s_vsMass_brNew_%s.png"%(SIGCH,options.channel,brnew));
    can_BSM.SaveAs("limitFigs/BSMLim%s_%s_vsMass_brNew_%s.pdf"%(SIGCH,options.channel,brnew));

    can_BSMsig = ROOT.TCanvas("can_BSMsig","can_BSMsig",630,600);
    hrl_BSMsig = can_BSMsig.DrawFrame(599,0.001,1001,max(gridMaxSigExp,gridMaxSigObs)*4);

    hrl_BSMsig.GetYaxis().SetTitle("#sigma #times BR_{WW}");
    hrl_BSMsig.GetYaxis().SetTitleOffset(1.35);
    hrl_BSMsig.GetYaxis().SetTitleSize(0.045);
    hrl_BSMsig.GetYaxis().SetTitleFont(42);

    hrl_BSMsig.GetXaxis().SetTitle("M_{H} (GeV)");
    hrl_BSMsig.GetXaxis().SetTitleSize(0.045);
    hrl_BSMsig.GetXaxis().SetTitleFont(42);

    hrl_BSMsig.GetYaxis().SetNdivisions(505);
    can_BSMsig.SetGridx(1);
    can_BSMsig.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.70,0.75,0.88);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    can_BSMsig.Update();
    can_BSMsig.RedrawAxis();
    can_BSMsig.RedrawAxis("g");
    can_BSMsig.Update();
    ROOT.gPad.SetLogy();
    
    for k in range(len(tGraphs_csXbr_exp)):
        tGraphs_csXbr_exp[k].SetLineStyle(1);
        tGraphs_csXbr_exp[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_obs[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_th[k].SetLineStyle(2);
        tGraphs_csXbr_th[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_exp[k].SetLineWidth(2);
        tGraphs_csXbr_obs[k].SetLineWidth(2);
        tGraphs_csXbr_exp[k].Draw("PL");
        if not options.blindObservedLine :        
            tGraphs_csXbr_exp[k].SetLineStyle(7);
            tGraphs_csXbr_exp[k].Draw("PL");
            tGraphs_csXbr_obs[k].Draw("PL");
        else:
            tGraphs_csXbr_exp[k].Draw("PL");
            
        tGraphs_csXbr_th[k].Draw("PL");

        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprime[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_exp[k],tmplabel,"L")
        if not options.blindObservedLine:
         tmplabel = "th., C'^{ 2} = %1.1f"%( float((cprime[k])/10.) )
         leg2.AddEntry(tGraphs_csXbr_th[k],tmplabel,"L");

    leg2.Draw();
    banner.Draw();
    label_sqrt.Draw();
    banner2.Draw();
    
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_%s_vsMass_brNew_%s_Sigma.png"%(SIGCH,options.channel,brnew));
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_%s_vsMass_brNew_%s_Sigma.pdf"%(SIGCH,options.channel,brnew));


###############################################
### Make the BSM Limit plot vs c' and brNew ###
###############################################
    
def makeBSMLimitPlotBRnew(SIGCH,mass):

    print "module ===> makeBSMLimits_vsBRnew";

    curcolors = [1,2,210,4,6,12,7];
    BRnew_x   = [0,0.1,0.2,0.3,0.4,0.5];
    massindex = {600:0,700:1,800:2,900:3,1000:4}
    massXS    = [];

    if SIGCH == "" or SIGCH == "_2jet":
	massXS.append((0.5230+0.09688));
        massXS.append((0.2288+0.06330));
        massXS.append((0.1095+0.04365));
        massXS.append((0.05684+0.03164));
        massXS.append((0.03163+0.02399));
    elif SIGCH == "_BulkGraviton_newxsec":
        massXS.append((0.5230));
        massXS.append((0.2288));
        massXS.append((0.1095));
        massXS.append((0.05684));
        massXS.append((0.03163));
    elif SIGCH == "_vbfH":
        massXS.append((0.09688));
        massXS.append((0.06330));
        massXS.append((0.04365));
        massXS.append((0.03164));
        massXS.append((0.02399));
    else:
        print "problem!"
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];

    gridMaxExp = -999;
    gridMaxSigExp = -999;
    gridMaxObs = -999;
    gridMaxSigObs = -999;

    tGraphs_exp = [];
    tGraphs_obs = [];
    tGraphs_csXbr_exp = [];
    tGraphs_csXbr_obs = [];
    tGraphs_csXbr_th = [];

    for j in range(len(cprime)):

        if cprime[j]*0.1 == 1.0 : continue ; 

        xbins           = array('f', []);
        ybins_exp       = array('f', []);
        ybins_obs       = array('f', []);
        ybins_csXbr_exp = array('f', []);
        ybins_csXbr_obs = array('f', []);
        ybins_csXbr_th  = array('f', []);
        nPoints = 0 ;
        for i in range(len(BRnew)):
            if (cprime[j]*0.1) > (1-BRnew[i]*0.1): continue ;
            nPoints = nPoints+1 ;
            curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass,options.channel,SIGCH,mass);
            curAsymLimits = getAsymLimits(curFile);
            xbins.append(BRnew_x[i]);
            ybins_exp.append( curAsymLimits[3] );
            ybins_obs.append( curAsymLimits[0] );
            ybins_csXbr_exp.append( curAsymLimits[3]*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );
            ybins_csXbr_obs.append( curAsymLimits[0]*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );
            ybins_csXbr_th.append( 1.*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );

            if gridMaxExp < curAsymLimits[3]: gridMaxExp = curAsymLimits[3];
            if gridMaxObs < curAsymLimits[0]: gridMaxObs = curAsymLimits[0];
            cscur = ( curAsymLimits[3]*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );
            if gridMaxSigExp < cscur: gridMaxSigExp = cscur;
            cscurObs = ( curAsymLimits[0]*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );
            if gridMaxSigObs < cscurObs: gridMaxSigObs = cscurObs;

        if not xbins or not ybins_exp or not ybins_obs : continue ;

        curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
        curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
        curGraph_obs.SetMarkerStyle(20);
        curGraph_obs.SetLineWidth(3);
        curGraph_obs.SetLineStyle(1);
        curGraph_obs.SetMarkerSize(1);
        curGraph_exp.SetMarkerSize(1.3);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(3);
        curGraph_exp.SetMarkerSize(1);
        curGraph_exp.SetMarkerStyle(24);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);

        curGraph_csXbr_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_exp);
        curGraph_csXbr_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_obs);
        curGraph_csXbr_th  = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_th);
        curGraph_obs.SetMarkerStyle(20);
        curGraph_obs.SetLineWidth(3);
        curGraph_obs.SetLineStyle(1);
        curGraph_obs.SetMarkerSize(1);
        curGraph_exp.SetMarkerSize(1.3);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(3);
        curGraph_exp.SetMarkerSize(1);
        curGraph_exp.SetMarkerStyle(24);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);

        curGraph_csXbr_exp.SetMarkerSize(1);
        curGraph_csXbr_exp.SetMarkerStyle(24);
        curGraph_csXbr_exp.SetMarkerColor(ROOT.kBlack);

        tGraphs_exp.append(curGraph_exp);
        tGraphs_obs.append(curGraph_obs);
        tGraphs_csXbr_exp.append(curGraph_csXbr_exp);
        tGraphs_csXbr_obs.append(curGraph_csXbr_obs);
        tGraphs_csXbr_th.append(curGraph_csXbr_th);

    setStyle();

    can_BSM = ROOT.TCanvas("can_BSM_BR","can_BSM",630,600);
    hrl_BSM = can_BSM.DrawFrame(0.0,0.0,0.5,max(gridMaxExp,gridMaxObs)*1.5);

    hrl_BSM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_BSM.GetYaxis().SetTitleOffset(1.35);
    hrl_BSM.GetYaxis().SetTitleSize(0.045);
    hrl_BSM.GetYaxis().SetTitleFont(42);

    hrl_BSM.GetXaxis().SetTitle("BR_{new}");
    hrl_BSM.GetXaxis().SetTitleSize(0.045);
    hrl_BSM.GetXaxis().SetTitleFont(42);

    hrl_BSM.GetYaxis().SetNdivisions(505);
    can_BSM.SetGridx(1);
    can_BSM.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    for k in range(len(tGraphs_exp)): 
        tGraphs_exp[k].SetLineStyle(1);
        tGraphs_exp[k].SetLineColor(curcolors[k]);
        tGraphs_exp[k].SetMarkerColor(curcolors[k]);
        tGraphs_obs[k].SetLineColor(curcolors[k]);
        tGraphs_obs[k].SetMarkerColor(curcolors[k]);
        tGraphs_exp[k].SetLineWidth(2);
        tGraphs_obs[k].SetLineWidth(2);
        if not options.blindObservedLine:
            tGraphs_exp[k].SetLineStyle(7);
            tGraphs_exp[k].Draw("PL");
            tGraphs_obs[k].Draw("PL");
        else:
            tGraphs_exp[k].Draw("PL");
            
        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprime[k])/10.) )
        leg2.AddEntry(tGraphs_exp[k],tmplabel,"L")
        if not options.blindObservedLine:
         tmplabel = "obs., C'^{ 2} = %1.1f    "%( float((cprime[k])/10.) )
         leg2.AddEntry(tGraphs_obs[k],tmplabel,"L")


    can_BSM.Update();
    can_BSM.RedrawAxis();
    can_BSM.RedrawAxis("g");
    can_BSM.Update();

    banner = TPaveText( 0.145, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.038);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.5,0.953,0.96,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.038);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");

    banner2 = TLatex(0.17,0.91,("Higgs Mass, %i GeV"%(mass)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    leg2.Draw();
    banner.Draw();
    banner2.Draw();
    label_sqrt.Draw();

    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu_vsBRnew_%i.png"%(SIGCH,mass));
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu_vsBRnew_%i.pdf"%(SIGCH,mass));

    can_BSMsig = ROOT.TCanvas("can_BSM_Sig","can_BSM",630,600);
    hrl_BSMSig = can_BSMsig.DrawFrame(0.,0.0,0.5,max(gridMaxSigExp,gridMaxSigObs)*1.5);

    hrl_BSMSig.GetYaxis().SetTitle("#sigma #times BR_{WW}");
    hrl_BSMSig.GetYaxis().SetTitleOffset(1.35);
    hrl_BSMSig.GetYaxis().SetTitleSize(0.045);
    hrl_BSMSig.GetYaxis().SetTitleFont(42);

    hrl_BSMSig.GetXaxis().SetTitle("BR_{new}");
    hrl_BSMSig.GetXaxis().SetTitleSize(0.045);
    hrl_BSMSig.GetXaxis().SetTitleFont(42);

    hrl_BSMSig.GetYaxis().SetNdivisions(505);
    can_BSMsig.SetGridx(1);
    can_BSMsig.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    for k in range(len(tGraphs_csXbr_exp)):
        tGraphs_csXbr_exp[k].SetLineStyle(1);
        tGraphs_csXbr_exp[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_obs[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_th[k].SetLineStyle(2);
        tGraphs_csXbr_th[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_exp[k].SetLineWidth(2);
        tGraphs_csXbr_obs[k].SetLineWidth(2);
        tGraphs_csXbr_exp[k].SetMarkerColor(curcolors[k]);
        tGraphs_csXbr_exp[k].Draw("PL");
        if not options.blindObservedLine:
            tGraphs_csXbr_exp[k].SetLineStyle(7);
            tGraphs_csXbr_exp[k].Draw("PL");
            tGraphs_csXbr_obs[k].Draw("PL");

        tmplabel = "exp., C'^{ 2} = %1.1f"%( float((cprime[k])/10.) );
        leg2.AddEntry(tGraphs_csXbr_exp[k],tmplabel,"L")

        if not options.blindObservedLine:
         tmplabel = "obs., C'^{ 2} = %1.1f"%( float((cprime[k])/10.) );
         leg2.AddEntry(tGraphs_csXbr_obs[k],tmplabel,"L")

    leg2.Draw();
    banner.Draw();
    banner2.Draw();
    label_sqrt.Draw();
     
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma_vsBRnew_%i.png"%(SIGCH,mass));
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma_vsBRnew_%i.pdf"%(SIGCH,mass));

#############################################
### Make the BSM 2D Scane vs c' and brNew ###
###################B##########################

def makeBSMLimitPlot2D( SIGCH, mass, contourListMassExp=0, contourListMassObs=0):


    stylePath = os.getenv("ROOTStyle");
    #gROOT.ProcessLine(".x "+stylePath+"/rootPalette.C");
    
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];            
    print "module ===> makeBSMLimits_2D";

    massindex = {600:0,700:1,800:2,900:3,1000:4}
    mass_XS  = [];
    BRnew_y  = [-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6];
    cprime_x = [-0.1,0.0,0.1,0.2,0.3,0.5,0.7,0.8];
    setStyle();
    
    if SIGCH == "" or SIGCH == "_2jet":
        mass_XS.append((0.5230+0.09688));
        mass_XS.append((0.2288+0.06330));
        mass_XS.append((0.1095+0.04365));
        mass_XS.append((0.05684+0.03164));
        mass_XS.append((0.03163+0.02399));
    elif SIGCH == "_BulkGraviton_newxsec":
        mass_XS.append((0.5230));
        mass_XS.append((0.2288));
        mass_XS.append((0.1095));
        mass_XS.append((0.05684));
        mass_XS.append((0.03163));
    elif SIGCH == "_vbfH":
        mass_XS.append((0.09688));
        mass_XS.append((0.06330));
        mass_XS.append((0.04365));
        mass_XS.append((0.03164));
        mass_XS.append((0.02399));
    else:
        print "problem!"

    setStyle();

    h2d_exp = TH2D("h2d_exp_%d"%(mass),"",200,cprime_x[0],cprime_x[len(cprime_x)-1],300,BRnew_y[0],BRnew_y[len(BRnew_y)-1]);
    h2d_obs = TH2D("h2d_obs_%d"%(mass),"",200,cprime_x[0],cprime_x[len(cprime_x)-1],300,BRnew_y[0],BRnew_y[len(BRnew_y)-1]);
    h2d_csXbr_exp = TH2D("h2d_csXbr_exp_%d"%(mass),"",200,cprime_x[0],cprime_x[len(cprime_x)-1],300,BRnew_y[0],BRnew_y[len(BRnew_y)-1]);
    h2d_csXbr_obs = TH2D("h2d_csXbr_obs_%d"%(mass),"",200,cprime_x[0],cprime_x[len(cprime_x)-1],300,BRnew_y[0],BRnew_y[len(BRnew_y)-1]);

    for j in range(len(cprime)):
        for i in range(len(BRnew)):
         gamFactor = (float(cprime[j])/10.)/(1-(float(BRnew[i])/10.));
         if gamFactor > 1 : continue ;
         binX = h2d_exp.GetXaxis().FindBin(cprime[j]*0.1);
         binY = h2d_exp.GetYaxis().FindBin(BRnew[i]*0.1);
                
         curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass,options.channel,SIGCH,mass); 
         curAsymLimits = getAsymLimits(curFile);
         h2d_exp.SetBinContent(binX+1,binY+1,curAsymLimits[3]);
         h2d_csXbr_exp.SetBinContent(binX+1,binY+1,curAsymLimits[3]*mass_XS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]]);
         h2d_obs.SetBinContent(binX+1,binY+1,curAsymLimits[0]);
         h2d_csXbr_obs.SetBinContent(binX+1,binY+1,curAsymLimits[0]*mass_XS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]]);

    g2d_obs = TGraph2D(h2d_obs);
    g2d_exp = TGraph2D(h2d_exp);
    
    for binX in range(h2d_exp.GetNbinsX()):
     for binY in range(h2d_exp.GetNbinsY()):
      if h2d_exp.GetBinContent(binX+1,binY+1) != 0 : continue ;   
      g2d_exp.Interpolate(h2d_exp.GetXaxis().GetBinCenter(binX+1),h2d_exp.GetYaxis().GetBinCenter(binY+1));

    for binX in range(h2d_obs.GetNbinsX()):
     for binY in range(h2d_obs.GetNbinsY()):
      if h2d_obs.GetBinContent(binX+1,binY+1) != 0 : continue ;   
      g2d_obs.Interpolate(h2d_obs.GetXaxis().GetBinCenter(binX+1),h2d_obs.GetYaxis().GetBinCenter(binY+1));

    h2d_exp = g2d_exp.GetHistogram();
    h2d_obs = g2d_obs.GetHistogram();
    
    for binX in range(h2d_exp.GetNbinsX()):
        for binY in range(h2d_exp.GetNbinsY()):
            if float(h2d_exp.GetXaxis().GetBinCenter(binX+1))/(1-(float(h2d_exp.GetYaxis().GetBinCenter(binY+1)))) > 1 :
               h2d_exp.SetBinContent(binX+1,binY+1,0);                                                       

    for binX in range(h2d_obs.GetNbinsX()):
        for binY in range(h2d_obs.GetNbinsY()):
            if float(h2d_obs.GetXaxis().GetBinCenter(binX+1))/(1-(float(h2d_obs.GetYaxis().GetBinCenter(binY+1)))) > 1 :
               h2d_obs.SetBinContent(binX+1,binY+1,0);                                                      

        
    h2d_obs.GetXaxis().SetTitle("C^{'2}");
    h2d_csXbr_exp.GetXaxis().SetTitle("C^{'2}");
    h2d_csXbr_obs.GetXaxis().SetTitle("C^{'2}");

    h2d_exp.GetYaxis().SetTitle("BR_{new}");
    h2d_obs.GetYaxis().SetTitle("BR_{new}");
    h2d_csXbr_exp.GetYaxis().SetTitle("BR_{new}");
    h2d_csXbr_obs.GetYaxis().SetTitle("BR_{new}");

    h2d_csXbr_obs.GetZaxis().SetLimits(0.,h2d_csXbr_obs.GetMaximum());
    h2d_csXbr_exp.GetZaxis().SetLimits(0.,h2d_csXbr_exp.GetMaximum());
    h2d_exp.GetZaxis().SetLimits(0,h2d_exp.GetMaximum());
    h2d_obs.GetZaxis().SetLimits(0,h2d_obs.GetMaximum());

    h2d_exp.SetLineStyle(0);
    h2d_exp.SetMarkerStyle(20);
    h2d_exp.GetXaxis().SetTitle("C^{'2}");
    h2d_exp.GetXaxis().SetNdivisions(504);
    h2d_exp.GetXaxis().SetLabelFont(42);
    h2d_exp.GetXaxis().SetLabelOffset(0.007);
    h2d_exp.GetXaxis().SetLabelSize(0.036);
    h2d_exp.GetXaxis().SetTitleSize(0.045);
    h2d_exp.GetXaxis().SetTitleOffset(1.02);
    h2d_exp.GetXaxis().SetTitleFont(42);
    h2d_exp.GetYaxis().SetTitle("BR_{new}");
    h2d_exp.GetYaxis().SetNdivisions(9);
    h2d_exp.GetYaxis().SetLabelFont(42);
    h2d_exp.GetYaxis().SetLabelOffset(0.007);
    h2d_exp.GetYaxis().SetLabelSize(0.036);
    h2d_exp.GetYaxis().SetTitleSize(0.045);
    h2d_exp.GetYaxis().SetTitleOffset(1.35);
    h2d_exp.GetYaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitle("signal strenght excluded at 95% C.L.");
    h2d_exp.GetZaxis().SetLabelFont(42);
    h2d_exp.GetZaxis().CenterTitle();
    h2d_exp.GetZaxis().SetLabelSize(0.025);
    h2d_exp.GetZaxis().SetTitleOffset(0.85);
    h2d_exp.GetZaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitleSize(0.05);

    h2d_obs.SetLineStyle(0);
    h2d_obs.SetMarkerStyle(20);
    h2d_obs.GetXaxis().SetTitle("C^{'2}");
    h2d_obs.GetXaxis().SetNdivisions(504);
    h2d_obs.GetXaxis().SetLabelFont(42);
    h2d_obs.GetXaxis().SetLabelOffset(0.007);
    h2d_obs.GetXaxis().SetLabelSize(0.036);
    h2d_obs.GetXaxis().SetTitleSize(0.045);
    h2d_obs.GetXaxis().SetTitleOffset(1.02);
    h2d_obs.GetXaxis().SetTitleFont(42);
    h2d_obs.GetYaxis().SetTitle("BR_{new}");
    h2d_obs.GetYaxis().SetNdivisions(9);
    h2d_obs.GetYaxis().SetLabelFont(42);
    h2d_obs.GetYaxis().SetLabelOffset(0.007);
    h2d_obs.GetYaxis().SetLabelSize(0.036);
    h2d_obs.GetYaxis().SetTitleSize(0.045);
    h2d_obs.GetYaxis().SetTitleOffset(1.25);
    h2d_obs.GetYaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitle("signal strenght excluded at 95% C.L.");
    h2d_obs.GetZaxis().SetLabelFont(42);
    h2d_obs.GetZaxis().SetLabelSize(0.026);
    h2d_obs.GetZaxis().SetTitleOffset(0.95);
    h2d_obs.GetZaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitleSize(0.05);

    banner = TPaveText( 0.155, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.033);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.45,0.953,0.89,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.033);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    
    banner2 = TLatex(0.2,0.9,("Higgs mass, %i GeV"%(mass)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);


    can1_BSM2D = ROOT.TCanvas("can1_BSM2D","can1_BSM2D",1,1,600,600);
    can1_BSM2D.SetHighLightColor(2);
    can1_BSM2D.SetFillColor(0);
    can1_BSM2D.SetBorderMode(0);
    can1_BSM2D.SetBorderSize(2);
    can1_BSM2D.SetTickx(1);
    can1_BSM2D.SetTicky(1);
    can1_BSM2D.SetLeftMargin(0.15);
    can1_BSM2D.SetRightMargin(0.17);
    can1_BSM2D.SetTopMargin(0.05);
    can1_BSM2D.SetBottomMargin(0.1);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);
    can1_BSM2D.SetGrid();    
    
    h2d_exp.GetXaxis().SetLimits(0.,0.7);
    h2d_exp.GetYaxis().SetLimits(0.,BRnew[len(BRnew)-1]*0.1);
    h2d_exp.GetXaxis().SetNdivisions(510);
    h2d_exp.GetYaxis().SetNdivisions(510); 
    h2d_exp.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());    
    counturLevel = [3.];
#    if mass == 1000 :
#     counturLevel[0] = 6.;
#    elif mass == 900 :
#     counturLevel[0] = 6.;
#    elif mass == 800 :
#     counturLevel[0] = 4.;
#    elif mass == 700 :
#     counturLevel[0] = 3.;
        
    h2d_exp.SetContour(1,array('d',counturLevel));    
    h2d_exp.Draw("cont,list");
    can1_BSM2D.Update();
    arrayList = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
    counturList = ROOT.TList();
    for i in range(arrayList.GetSize()):
        List_tmp = arrayList.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 if j == 0: counturList.Add(gr1.Clone());
                     
    h2d_exp.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());      

    banner3 = ROOT.TLegend(0.2,0.7,0.4,0.8);
    banner3.AddEntry(counturList.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
    banner3.SetTextSize(0.032);
    banner3.SetFillStyle(0);
    banner3.SetFillColor(0);
    banner3.SetShadowColor(0);
    banner3.SetTextFont(42);
    banner3.SetBorderSize(0);

    counturLevel = [4.];
#    if mass == 1000 :
#     counturLevel[0] = 7.;
#    elif mass == 900 :
#     counturLevel[0] = 7.;
#    elif mass == 800 :
#     counturLevel[0] = 5.;
#    elif mass == 700 :
#     counturLevel[0] = 4.;

    h2d_exp.SetContour(1,array('d',counturLevel));    
    h2d_exp.Draw("cont,list");
    can1_BSM2D.Update();
    arrayList_2 = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
    counturList_2 = ROOT.TList();
    for i in range(arrayList_2.GetSize()):
        List_tmp = arrayList_2.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 gr1.SetLineStyle(4);                 
                 if j == 0:
                     counturList_2.Add(gr1.Clone());
                     
    h2d_exp.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());      

    h2d_exp.Draw("colz");
    counturList.Draw("lsame");
    banner.Draw();
    label_sqrt.Draw();
    banner2.Draw();
    counturList_2.Draw("lsame");
    banner3.AddEntry(counturList_2.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
    banner3.Draw("same");

    gPad.Update();
    
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_%i.png"%(SIGCH,mass));
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_%i.pdf"%(SIGCH,mass));


    contourListMassExp[mass].append(counturList);
    

    if not options.blindObservedLine:

     can2_BSM2D = ROOT.TCanvas("can2_BSM2D","can2_BSM2D",1,1,600,600);
     can2_BSM2D.SetHighLightColor(2);
     can2_BSM2D.SetFillColor(0);
     can2_BSM2D.SetBorderMode(0);
     can2_BSM2D.SetBorderSize(2);
     can2_BSM2D.SetTickx(1);
     can2_BSM2D.SetTicky(1);
     can2_BSM2D.SetLeftMargin(0.15);
     can2_BSM2D.SetRightMargin(0.17);
     can2_BSM2D.SetTopMargin(0.05);
     can2_BSM2D.SetBottomMargin(0.1);
     can2_BSM2D.SetFrameFillStyle(0);
     can2_BSM2D.SetFrameBorderMode(0);
     can2_BSM2D.SetFrameFillStyle(0);
     can2_BSM2D.SetFrameBorderMode(0);
     can2_BSM2D.SetGrid();    
    
     h2d_obs.GetXaxis().SetLimits(0.,0.7);
     h2d_obs.GetYaxis().SetLimits(0.,BRnew[len(BRnew)-1]*0.1);
     h2d_obs.GetXaxis().SetNdivisions(510);
     h2d_obs.GetYaxis().SetNdivisions(510); 
     h2d_obs.SetContour(h2d_obs.GetNbinsX()*h2d_obs.GetNbinsY());
     counturLevel = [3.];
#     if mass == 1000 :
#      counturLevel[0] = 6.;
#     elif mass == 900 :
#      counturLevel[0] = 6.;
#     elif mass == 800 :
#      counturLevel[0] = 5.;
#     elif mass == 700 :
#      counturLevel[0] = 4.;

     h2d_obs.SetContour(1,array('d',counturLevel));    
     h2d_obs.Draw("cont,list");
     can2_BSM2D.Update();
     arrayList = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
     counturList = ROOT.TList();
     for i in range(arrayList.GetSize()):
        List_tmp = arrayList.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 if j == 0: counturList.Add(gr1.Clone());
                     
     h2d_obs.SetContour(h2d_obs.GetNbinsX()*h2d_obs.GetNbinsY());      
 
     banner3 = ROOT.TLegend(0.2,0.7,0.4,0.8);
     banner3.AddEntry(counturList.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
     banner3.SetTextSize(0.032);
     banner3.SetFillStyle(0);
     banner3.SetFillColor(0);
     banner3.SetShadowColor(0);
     banner3.SetTextFont(42);
     banner3.SetBorderSize(0);

     counturLevel = [4.];
#     if mass == 1000 :
#      counturLevel[0] = 7.;
#     elif mass == 900 :
#      counturLevel[0] = 7.;
#     elif mass == 800 :
#      counturLevel[0] = 6.;
#     elif mass == 700 :
#      counturLevel[0] = 4.;
     h2d_obs.SetContour(1,array('d',counturLevel));    
     h2d_obs.Draw("cont,list");
     can2_BSM2D.Update();
     arrayList_2 = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
     counturList_2 = ROOT.TList();
     for i in range(arrayList_2.GetSize()):
        List_tmp = arrayList_2.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 gr1.SetLineStyle(4);                 
                 if j == 0:
                     counturList_2.Add(gr1.Clone());
                     
     h2d_obs.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());      
 
     h2d_obs.Draw("colz");
     counturList.Draw("lsame");
     banner.Draw();
     label_sqrt.Draw();
     banner2.Draw();
     counturList_2.Draw("lsame");
     banner3.AddEntry(counturList_2.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
     banner3.Draw("same");
     gPad.Update();
    
     can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_%i.png"%(SIGCH,mass));
     can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_%i.pdf"%(SIGCH,mass));

     contourListMassObs[mass].append(counturList);

#############################################
### Make the BSM 2D Scane vs c' and brNew ###
#############################################

def makeBSMLimitPlot2DBRnew( SIGCH, brNew, contourListBrNewExp =0, counturListBrNewObs =0 ):


    stylePath = os.getenv("ROOTStyle");
#    gROOT.ProcessLine(".x "+stylePath+"/rootPalette.C");

    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];            

    print "module ===> makeBSMLimits_2DMass";

    massindex = {600:0,700:1,800:2,900:3,1000:4}
    mass_XS  = [];
    cprime_x = array('d',[-0.1,0.0,0.1,0.2,0.3,0.5,0.7,0.8]);
    mass_x   = array('d',[580,600,700,800,900,1000,1020]);

    setStyle();

    if SIGCH == "" or SIGCH == "_2jet":
        mass_XS.append((0.5230+0.09688));
        mass_XS.append((0.2288+0.06330));
        mass_XS.append((0.1095+0.04365));
        mass_XS.append((0.05684+0.03164));
        mass_XS.append((0.03163+0.02399));
    elif SIGCH == "_BulkGraviton_newxsec":
        mass_XS.append((0.5230));
        mass_XS.append((0.2288));
        mass_XS.append((0.1095));
        mass_XS.append((0.05684));
        mass_XS.append((0.03163));
    elif SIGCH == "_vbfH":
        mass_XS.append((0.09688));
        mass_XS.append((0.06330));
        mass_XS.append((0.04365));
        mass_XS.append((0.03164));
        mass_XS.append((0.02399));
    else:
        print "problem!"


    h2d_exp = TH2D("h2d_exp_brNew_%02d"%(brNew),"",200,mass_x[0],mass_x[len(mass_x)-1],300,cprime_x[0],cprime_x[len(cprime_x)-1]);
    h2d_obs = TH2D("h2d_obs_brNew_%02d"%(brNew),"",200,mass_x[0],mass_x[len(mass_x)-1],300,cprime_x[0],cprime_x[len(cprime_x)-1]);
    h2d_csXbr_exp = TH2D("h2d_csXbr_exp_brNew_%02d"%(brNew),"",200,mass_x[0],mass_x[len(mass_x)-1],300,cprime_x[0],cprime_x[len(cprime_x)-1]);
    h2d_csXbr_obs = TH2D("h2d_csXbr_obs_brNew_%02d"%(brNew),"",200,mass_x[0],mass_x[len(mass_x)-1],300,cprime_x[0],cprime_x[len(cprime_x)-1]);

    for j in range(len(mass)):
        for i in range(len(cprime)):
            gamFactor = (float(cprime[i])/10.)/(1-(float(brNew*0.1)));
            if gamFactor > 1 : continue ;
            binX = h2d_exp.GetXaxis().FindBin(mass[j]);
            binY = h2d_exp.GetYaxis().FindBin(cprime[i]*0.1);
                           
            curFile = "higgsCombinewwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.Asymptotic.mH%03d.root"%(mass[j],options.channel,SIGCH,cprime[i],brNew,mass[j]);
            curAsymLimits = getAsymLimits(curFile);
            h2d_exp.SetBinContent(binX+1,binY+1,curAsymLimits[3]);
            h2d_obs.SetBinContent(binX+1,binY+1,curAsymLimits[0]);
            h2d_csXbr_exp.SetBinContent(binX+1,binY+1,curAsymLimits[3]*mass_XS[massindex[mass[j]]]*cprime[i]*0.1*(1-brNew*0.1)*massBRWW[massindex[mass[j]]]);
            h2d_csXbr_obs.SetBinContent(binX+1,binY+1,curAsymLimits[0]*mass_XS[massindex[mass[j]]]*cprime[i]*0.1*(1-brNew*0.1)*massBRWW[massindex[mass[j]]]);

    
    g2d_obs = TGraph2D(h2d_obs);
    g2d_exp = TGraph2D(h2d_exp);
    
    for binX in range(h2d_exp.GetNbinsX()):
     for binY in range(h2d_exp.GetNbinsY()):
      if h2d_exp.GetBinContent(binX+1,binY+1) != 0 : continue ;   
      g2d_exp.Interpolate(h2d_exp.GetXaxis().GetBinCenter(binX+1),h2d_exp.GetYaxis().GetBinCenter(binY+1));

    for binX in range(h2d_obs.GetNbinsX()):
     for binY in range(h2d_obs.GetNbinsY()):
      if h2d_obs.GetBinContent(binX+1,binY+1) != 0 : continue ;   
      g2d_obs.Interpolate(h2d_obs.GetXaxis().GetBinCenter(binX+1),h2d_obs.GetYaxis().GetBinCenter(binY+1));

    h2d_exp = g2d_exp.GetHistogram();
    h2d_obs = g2d_obs.GetHistogram();
    
    for binX in range(h2d_exp.GetNbinsX()):
        for binY in range(h2d_exp.GetNbinsY()):
            if float(h2d_exp.GetYaxis().GetBinCenter(binY+1))/(1-(brNew*0.1)) > 1 :
               h2d_exp.SetBinContent(binX+1,binY+1,0);                                                       

    for binX in range(h2d_obs.GetNbinsX()):
        for binY in range(h2d_obs.GetNbinsY()):
            if float(h2d_obs.GetYaxis().GetBinCenter(binY+1))/(1-(brNew*0.1)) > 1 :
               h2d_obs.SetBinContent(binX+1,binY+1,0);                                                      
    
    h2d_obs.GetYaxis().SetTitle("C^{'2}");
    h2d_csXbr_exp.GetYaxis().SetTitle("C^{'2}");
    h2d_csXbr_obs.GetYaxis().SetTitle("C^{'2}");

    h2d_exp.GetXaxis().SetTitle("m_{H} (GeV)");
    h2d_obs.GetXaxis().SetTitle("m_{H} (GeV)");
    h2d_csXbr_exp.GetXaxis().SetTitle("m_{H} (GeV)");
    h2d_csXbr_obs.GetXaxis().SetTitle("m_{H} (GeV)");

    h2d_csXbr_obs.GetZaxis().SetLimits(0.,h2d_csXbr_obs.GetMaximum());
    h2d_csXbr_exp.GetZaxis().SetLimits(0.,h2d_csXbr_exp.GetMaximum());
    h2d_exp.GetZaxis().SetLimits(0,h2d_exp.GetMaximum());
    h2d_obs.GetZaxis().SetLimits(0,h2d_obs.GetMaximum());

    h2d_exp.SetLineStyle(0);
    h2d_exp.SetMarkerStyle(20);
    h2d_exp.GetYaxis().SetTitle("C^{'2}");
    h2d_exp.GetXaxis().SetNdivisions(504);
    h2d_exp.GetXaxis().SetLabelFont(42);
    h2d_exp.GetXaxis().SetLabelOffset(0.007);
    h2d_exp.GetXaxis().SetLabelSize(0.036);
    h2d_exp.GetXaxis().SetTitleSize(0.045);
    h2d_exp.GetXaxis().SetTitleOffset(1.02);
    h2d_exp.GetXaxis().SetTitleFont(42);
    h2d_exp.GetXaxis().SetTitle("m_{H} (GeV)");
    h2d_exp.GetYaxis().SetNdivisions(9);
    h2d_exp.GetYaxis().SetLabelFont(42);
    h2d_exp.GetYaxis().SetLabelOffset(0.007);
    h2d_exp.GetYaxis().SetLabelSize(0.036);
    h2d_exp.GetYaxis().SetTitleSize(0.045);
    h2d_exp.GetYaxis().SetTitleOffset(1.35);
    h2d_exp.GetYaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitle("signal strenght excluded at 95% C.L.");
    h2d_exp.GetZaxis().SetLabelFont(42);
    h2d_exp.GetZaxis().CenterTitle();
    h2d_exp.GetZaxis().SetLabelSize(0.025);
    h2d_exp.GetZaxis().SetTitleOffset(0.85);
    h2d_exp.GetZaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitleSize(0.05);

    h2d_obs.SetLineStyle(0);
    h2d_obs.SetMarkerStyle(20);
    h2d_obs.GetYaxis().SetTitle("C^{'2}");
    h2d_obs.GetXaxis().SetNdivisions(504);
    h2d_obs.GetXaxis().SetLabelFont(42);
    h2d_obs.GetXaxis().SetLabelOffset(0.007);
    h2d_obs.GetXaxis().SetLabelSize(0.036);
    h2d_obs.GetXaxis().SetTitleSize(0.045);
    h2d_obs.GetXaxis().SetTitleOffset(1.02);
    h2d_obs.GetXaxis().SetTitleFont(42);
    h2d_obs.GetXaxis().SetTitle("m_{H} (GeV)");
    h2d_obs.GetYaxis().SetNdivisions(9);
    h2d_obs.GetYaxis().SetLabelFont(42);
    h2d_obs.GetYaxis().SetLabelOffset(0.007);
    h2d_obs.GetYaxis().SetLabelSize(0.036);
    h2d_obs.GetYaxis().SetTitleSize(0.045);
    h2d_obs.GetYaxis().SetTitleOffset(1.25);
    h2d_obs.GetYaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitle("signal strenght excluded at 95% C.L.");
    h2d_obs.GetZaxis().SetLabelFont(42);
    h2d_obs.GetZaxis().SetLabelSize(0.026);
    h2d_obs.GetZaxis().SetTitleOffset(0.95);
    h2d_obs.GetZaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitleSize(0.05);

    banner = TPaveText( 0.155, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.033);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.45,0.953,0.89,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.033);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    
    banner2 = TLatex(0.2,0.2,("BR_{new} = %.1f "%(brNew*0.1)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    can1_BSM2D = ROOT.TCanvas("can1_BSM2D","can1_BSM2D",1,1,600,600);
    can1_BSM2D.SetHighLightColor(2);
    can1_BSM2D.SetFillColor(0);
    can1_BSM2D.SetBorderMode(0);
    can1_BSM2D.SetBorderSize(2);
    can1_BSM2D.SetTickx(1);
    can1_BSM2D.SetTicky(1);
    can1_BSM2D.SetLeftMargin(0.15);
    can1_BSM2D.SetRightMargin(0.17);
    can1_BSM2D.SetTopMargin(0.05);
    can1_BSM2D.SetBottomMargin(0.1);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);
    can1_BSM2D.SetGrid();    

    h2d_exp.GetXaxis().SetLimits(mass[0],mass[len(mass)-1]);
    h2d_exp.GetYaxis().SetLimits(0.,cprime[len(cprime)-2]*0.1);
    h2d_exp.GetXaxis().SetNdivisions(510);
    h2d_exp.GetYaxis().SetNdivisions(510); 
    h2d_exp.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());
    counturLevel = [3.];
#    if brNew*0.1 == 0.2:
#     counturLevel[0] = 3.;
#    elif brNew*0.1 == 0.3:
#     counturLevel[0] = 4.;
#    elif brNew*0.1 == 0.4:
#     counturLevel[0] = 5.;
#    elif brNew*0.1 == 0.5:
#     counturLevel[0] = 6.;
        
    h2d_exp.SetContour(1,array('d',counturLevel));    
    h2d_exp.Draw("cont,list");
    can1_BSM2D.Update();
    arrayList = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
    counturList = ROOT.TList();
    for i in range(arrayList.GetSize()):
        List_tmp = arrayList.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 if j == 0: counturList.Add(gr1.Clone());
                     
    h2d_exp.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());      
    banner3 = ROOT.TLegend(0.45,0.15,0.6,0.3);
    banner3.AddEntry(counturList.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
    banner3.SetTextSize(0.032);
    banner3.SetFillStyle(0);
    banner3.SetFillColor(0);
    banner3.SetShadowColor(0);
    banner3.SetTextFont(42);
    banner3.SetBorderSize(0);

    counturLevel = [4.];
#    if brNew*0.1 == 0.2:
#     counturLevel[0] = 4.;
#    elif brNew*0.1 == 0.3:
#     counturLevel[0] = 6.;
#    elif brNew*0.1 == 0.4:
#     counturLevel[0] = 7.;
#    elif brNew*0.1 == 0.5:
#     counturLevel[0] = 8.;

    h2d_exp.SetContour(1,array('d',counturLevel));    
    h2d_exp.Draw("cont,list");
    can1_BSM2D.Update();
    arrayList_2 = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
    counturList_2 = ROOT.TList();
    for i in range(arrayList_2.GetSize()):
        List_tmp = arrayList_2.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 gr1.SetLineStyle(4);                 
                 if j == 0:
                     counturList_2.Add(gr1.Clone());
                     
    h2d_exp.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());      
    h2d_exp.Draw("colz");
    counturList.Draw("lsame");
    banner.Draw();
    label_sqrt.Draw();
    banner2.Draw();
    counturList_2.Draw("lsame");
    banner3.AddEntry(counturList_2.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
    banner3.Draw("same");

    gPad.Update();

    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_brNew_%0.1f.png"%(SIGCH,brNew*0.1));
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_brNew_%0.1f.pdf"%(SIGCH,brNew*0.1));

    contourListBrNewExp[brNew*0.1].append(counturList) ;

    if not options.blindObservedLine:

     can2_BSM2D = ROOT.TCanvas("can2_BSM2D","can2_BSM2D",1,1,600,600);
     can2_BSM2D.SetHighLightColor(2);
     can2_BSM2D.SetFillColor(0);
     can2_BSM2D.SetBorderMode(0);
     can2_BSM2D.SetBorderSize(2);
     can2_BSM2D.SetTickx(1);
     can2_BSM2D.SetTicky(1);
     can2_BSM2D.SetLeftMargin(0.15);
     can2_BSM2D.SetRightMargin(0.17);
     can2_BSM2D.SetTopMargin(0.05);
     can2_BSM2D.SetBottomMargin(0.1);
     can2_BSM2D.SetFrameFillStyle(0);
     can2_BSM2D.SetFrameBorderMode(0);
     can2_BSM2D.SetFrameFillStyle(0);
     can2_BSM2D.SetFrameBorderMode(0);
     can2_BSM2D.SetGrid();    
    
     h2d_obs.GetXaxis().SetLimits(mass[0],mass[len(mass)-1]);
     h2d_obs.GetYaxis().SetLimits(cprime[0]*0.1,cprime[len(cprime)-2]*0.1);
     h2d_obs.GetXaxis().SetNdivisions(510);
     h2d_obs.GetYaxis().SetNdivisions(510); 
     h2d_obs.SetContour(h2d_obs.GetNbinsX()*h2d_obs.GetNbinsY());
     counturLevel = [3.];
#    if brNew*0.1 == 0.2:
#     counturLevel[0] = 3.;
#    elif brNew*0.1 == 0.3:
#     counturLevel[0] = 4.;
#    elif brNew*0.1 == 0.4:
#     counturLevel[0] = 5.;
#    elif brNew*0.1 == 0.5:
#     counturLevel[0] = 6.;
     h2d_obs.SetContour(1,array('d',counturLevel));    
     h2d_obs.Draw("cont,list");
     can2_BSM2D.Update();
     arrayList = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
     counturList = ROOT.TList();
     for i in range(arrayList.GetSize()):
        List_tmp = arrayList.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 if j == 0: counturList.Add(gr1.Clone());
                     
     h2d_obs.SetContour(h2d_obs.GetNbinsX()*h2d_obs.GetNbinsY());      
 
     banner3 = ROOT.TLegend(0.45,0.15,0.6,0.3);
     banner3.AddEntry(counturList.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
     banner3.SetTextSize(0.032);
     banner3.SetFillStyle(0);
     banner3.SetFillColor(0);
     banner3.SetShadowColor(0);
     banner3.SetTextFont(42);
     banner3.SetBorderSize(0);

     counturLevel = [4.];
#    if brNew*0.1 == 0.2:
#     counturLevel[0] = 3.;
#    elif brNew*0.1 == 0.3:
#     counturLevel[0] = 4.;
#    elif brNew*0.1 == 0.4:
#     counturLevel[0] = 5.;
#    elif brNew*0.1 == 0.5:
#     counturLevel[0] = 6.;

     h2d_obs.SetContour(1,array('d',counturLevel));    
     h2d_obs.Draw("cont,list");
     can2_BSM2D.Update();
     arrayList_2 = ROOT.gROOT.GetListOfSpecials().FindObject("contours");
     counturList_2 = ROOT.TList();
     for i in range(arrayList_2.GetSize()):
        List_tmp = arrayList_2.At(i);
        for j in range(List_tmp.GetSize()):
                 gr1 = TGraph (List_tmp.At(j));
                 gr1.SetLineColor(ROOT.kBlack);
                 gr1.SetLineWidth(2);
                 gr1.SetLineStyle(4);                 
                 if j == 0:
                     counturList_2.Add(gr1.Clone());
                     
     h2d_obs.SetContour(h2d_exp.GetNbinsX()*h2d_exp.GetNbinsY());      
 
     h2d_obs.Draw("colz");
     counturList.Draw("lsame");
     banner.Draw();
     label_sqrt.Draw();
     banner2.Draw();
     counturList_2.Draw("lsame");
     banner3.AddEntry(counturList_2.At(0),"%d #times #sigma_{Th} 95%s C.L. Limit"%(counturLevel[0],"%"),"l");
     banner3.Draw("same");
     gPad.Update();
    
     can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_brNew_%0.1f.png"%(SIGCH,brNew*0.1));
     can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_brNew_%0.1f.pdf"%(SIGCH,brNew*0.1));

     contourListBrNewObs[brNew*0.1].append(counturList);

def makeContourPlotMass(contourListMassExp,contourListMassObs):

  setStyle();

  can1_BSM = ROOT.TCanvas("can1_BSM","can1_BSM",1,1,600,600);
  can1_BSM.SetHighLightColor(2);
  can1_BSM.SetFillColor(0);
  can1_BSM.SetBorderMode(0);
  can1_BSM.SetBorderSize(2);
  can1_BSM.SetTickx(1);
  can1_BSM.SetTicky(1);
  can1_BSM.SetLeftMargin(0.15);
  can1_BSM.SetRightMargin(0.095);
  can1_BSM.SetTopMargin(0.05);
  can1_BSM.SetBottomMargin(0.10);
  can1_BSM.SetFrameFillStyle(0);
  can1_BSM.SetFrameBorderMode(0);
  can1_BSM.SetFrameFillStyle(0);
  can1_BSM.SetFrameBorderMode(0); 

  hrl = can1_BSM.DrawFrame(0.,BRnew[0]*0.1,cprime[len(cprime)-2]*0.1,BRnew[len(BRnew)-1]*0.1+0.3);
  hrl.GetYaxis().SetTitle("BR_{new}");
  hrl.GetXaxis().SetTitle("C^{'2}");
  hrl.GetXaxis().SetNdivisions(504);
  hrl.GetXaxis().SetLabelFont(42);
  hrl.GetXaxis().SetLabelOffset(0.007);
  hrl.GetXaxis().SetLabelSize(0.036);
  hrl.GetXaxis().SetTitleSize(0.045);
  hrl.GetXaxis().SetTitleOffset(1.02);
  hrl.GetXaxis().SetNdivisions(510);
  hrl.GetYaxis().SetNdivisions(510);
  hrl.GetXaxis().SetTitleFont(42);
  hrl.GetYaxis().SetLabelFont(42);
  hrl.GetYaxis().SetLabelOffset(0.007);
  hrl.GetYaxis().SetLabelSize(0.036);
  hrl.GetYaxis().SetTitleSize(0.045);
  hrl.GetYaxis().SetTitleOffset(1.35);
  hrl.GetYaxis().SetTitleFont(42);
  hrl.GetZaxis().SetTitle("signal strenght excluded at 95% C.L.");
  hrl.GetZaxis().SetLabelFont(42);
  hrl.GetZaxis().CenterTitle();
  hrl.GetZaxis().SetLabelSize(0.025);
  hrl.GetZaxis().SetTitleOffset(0.85);
  hrl.GetZaxis().SetTitleFont(42);
  hrl.GetZaxis().SetTitleSize(0.05);
  hrl.Draw("z");

  banner3 = ROOT.TLegend(0.243,0.63,0.69,0.92);

  list_temp = TList(); graph_temp = TGraph();
  list_temp = contourListMassExp[600];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineStyle(7);
  graph_temp.SetLineColor(ROOT.kBlack); 
  graph_temp.Draw("l") 
  banner3.AddEntry(graph_temp,"Expected 3 #times #sigma_{Th} m_{H}=600 GeV","l");

  list_temp = contourListMassObs[600];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineColor(ROOT.kBlack); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Observed 3 #times #sigma_{Th} m_{H}=600 GeV","l");

  list_temp = contourListMassExp[700];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineStyle(7);
  graph_temp.SetLineColor(ROOT.kRed); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Expected 3 #times #sigma_{Th} m_{H}=700 GeV","l");

  list_temp = contourListMassObs[700];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineColor(ROOT.kRed); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Observed 3 #times #sigma_{Th} m_{H}=700 GeV","l");

  list_temp = contourListMassExp[800];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineStyle(7);
  graph_temp.SetLineColor(ROOT.kBlue); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Expected 3 #times #sigma_{Th} m_{H}=800 GeV","l");

  list_temp = contourListMassObs[900];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineColor(ROOT.kBlue); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Observed 3 #times #sigma_{Th} m_{H}=800 GeV","l");

  banner = TPaveText(0.185, 0.953, 0.66, 0.975, "brNDC");
  banner.SetFillColor(ROOT.kWhite);
  banner.SetTextSize(0.033);
  banner.SetTextAlign(11);
  banner.SetTextFont(62);
  banner.SetBorderSize(0);
  leftText = "CMS Preliminary";
  banner.AddText(leftText);
  banner.Draw();

  label_sqrt = TPaveText(0.45,0.953,0.84,0.975, "brNDC");
  label_sqrt.SetFillColor(ROOT.kWhite);
  label_sqrt.SetBorderSize(0);
  label_sqrt.SetTextSize(0.033);
  label_sqrt.SetTextFont(62);
  label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
  label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    
  banner.Draw();
  label_sqrt.Draw();

  line = ROOT.TF1("line","0.5",0.,cprime[len(cprime)-2]*0.1);
  line.SetLineWidth(2);
  line.SetLineColor(ROOT.kMagenta);
  
  banner3.SetTextSize(0.032);
  banner3.SetFillStyle(0);
  banner3.SetFillColor(0);
  banner3.SetShadowColor(0);
  banner3.SetTextFont(42);
  banner3.SetBorderSize(0);
  banner3.Draw("same");
  line.Draw("same")

  can1_BSM.SaveAs("limitFigs/ContourMass.png");
  can1_BSM.SaveAs("limitFigs/ContourMass.pdf");
  can1_BSM.SaveAs("limitFigs/ContourMass.root");

def makeContourPlotBrNew(contourListMassExp,contourListMassObs):

  setStyle();

  can1_BSM = ROOT.TCanvas("can1_BSM","can1_BSM",1,1,600,600);
  can1_BSM.SetHighLightColor(2);
  can1_BSM.SetFillColor(0);
  can1_BSM.SetBorderMode(0);
  can1_BSM.SetBorderSize(2);
  can1_BSM.SetTickx(1);
  can1_BSM.SetTicky(1);
  can1_BSM.SetLeftMargin(0.15);
  can1_BSM.SetRightMargin(0.095);
  can1_BSM.SetTopMargin(0.05);
  can1_BSM.SetBottomMargin(0.10);
  can1_BSM.SetFrameFillStyle(0);
  can1_BSM.SetFrameBorderMode(0);
  can1_BSM.SetFrameFillStyle(0);
  can1_BSM.SetFrameBorderMode(0); 

  hrl = can1_BSM.DrawFrame(mass[0],0.,mass[len(mass)-1],cprime[len(cprime)-2]*0.1+0.4);
  hrl.GetYaxis().SetTitle("C^{'2}");
  hrl.GetXaxis().SetTitle("m_{H} (GeV)");
  hrl.GetXaxis().SetNdivisions(504);
  hrl.GetXaxis().SetLabelFont(42);
  hrl.GetXaxis().SetLabelOffset(0.007);
  hrl.GetXaxis().SetLabelSize(0.036);
  hrl.GetXaxis().SetTitleSize(0.045);
  hrl.GetXaxis().SetTitleOffset(1.02);
  hrl.GetXaxis().SetNdivisions(510);
  hrl.GetYaxis().SetNdivisions(510);
  hrl.GetXaxis().SetTitleFont(42);
  hrl.GetYaxis().SetLabelFont(42);
  hrl.GetYaxis().SetLabelOffset(0.007);
  hrl.GetYaxis().SetLabelSize(0.036);
  hrl.GetYaxis().SetTitleSize(0.045);
  hrl.GetYaxis().SetTitleOffset(1.35);
  hrl.GetYaxis().SetTitleFont(42);
  hrl.GetZaxis().SetTitle("signal strenght excluded at 95% C.L.");
  hrl.GetZaxis().SetLabelFont(42);
  hrl.GetZaxis().CenterTitle();
  hrl.GetZaxis().SetLabelSize(0.025);
  hrl.GetZaxis().SetTitleOffset(0.85);
  hrl.GetZaxis().SetTitleFont(42);
  hrl.GetZaxis().SetTitleSize(0.05);
  hrl.Draw("z");

  banner3 = ROOT.TLegend(0.22,0.66,0.80,0.91);

  list_temp = TList(); graph_temp = TGraph();
  list_temp = contourListMassExp[00];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineStyle(7);
  graph_temp.SetLineColor(ROOT.kBlack); 
  graph_temp.Draw("l") 
  banner3.AddEntry(graph_temp,"Expected 3 #times #sigma_{Th} BR_{new}=0","l");

  list_temp = contourListMassObs[0.0];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineColor(ROOT.kBlack); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Observed 3 #times #sigma_{Th} BR_{new}=0","l");

  list_temp = contourListMassExp[0.1];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineStyle(7);
  graph_temp.SetLineColor(ROOT.kRed); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Expected 3 #times #sigma_{Th} BR_{new}=0.1","l");

  list_temp = contourListMassObs[0.1];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineColor(ROOT.kRed); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Observed 3 #times #sigma_{Th} BR_{new}=0.1","l");

  list_temp = contourListMassExp[0.2];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineStyle(7);
  graph_temp.SetLineColor(ROOT.kBlue); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Expected 3 #times #sigma_{Th} BR_{new}=0.2","l");

  list_temp = contourListMassObs[0.2];
  list_temp = list_temp[0];
  graph_temp = list_temp.At(0);
  graph_temp.SetLineWidth(2);
  graph_temp.SetLineColor(ROOT.kBlue); 
  graph_temp.Draw("lsame") 
  banner3.AddEntry(graph_temp,"Observed 3 #times #sigma_{Th} BR_{new}=0.2","l");

  banner = TPaveText(0.185, 0.953, 0.66, 0.975, "brNDC");
  banner.SetFillColor(ROOT.kWhite);
  banner.SetTextSize(0.033);
  banner.SetTextAlign(11);
  banner.SetTextFont(62);
  banner.SetBorderSize(0);
  leftText = "CMS Preliminary";
  banner.AddText(leftText);
  banner.Draw();

  label_sqrt = TPaveText(0.45,0.953,0.84,0.975, "brNDC");
  label_sqrt.SetFillColor(ROOT.kWhite);
  label_sqrt.SetBorderSize(0);
  label_sqrt.SetTextSize(0.033);
  label_sqrt.SetTextFont(62);
  label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
  label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    
  banner.Draw();
  label_sqrt.Draw();

  banner3.SetTextSize(0.032);
  banner3.SetFillStyle(0);
  banner3.SetFillColor(0);
  banner3.SetShadowColor(0);
  banner3.SetTextFont(42);
  banner3.SetBorderSize(0);
  banner3.Draw("same");

  line = ROOT.TF1("line","0.7",mass[0],mass[len(mass)-1]);
  line.SetLineWidth(2);
  line.SetLineColor(ROOT.kMagenta);
  line.Draw("same");

  can1_BSM.SaveAs("limitFigs/ContourBrNew.png");
  can1_BSM.SaveAs("limitFigs/ContourBrNew.pdf");
  can1_BSM.SaveAs("limitFigs/ContourBrNew.root");

##################################
########### Main Code ############
##################################    

if __name__ == '__main__':

        
    CHAN = options.channel;
    DIR  = options.datacardDIR;
    if options.sigChannel !="": 
     SIGCH = options.jetBin+"_"+options.sigChannel;
    else:
     SIGCH = options.jetBin;
        
    moreCombineOpts = "";
    print "channel ",CHAN," directiory ",DIR," signal channel ",SIGCH," more options ",moreCombineOpts ;
    
    if options.makeCards and DIR != "":
        if not os.path.isdir(DIR):
            os.system("mkdir "+DIR);
        else: 
            print "Directory "+DIR+" already exists...";

    if options.computeLimits or options.plotLimits: os.chdir(DIR);

    nMasses = len(mass);
    mLo = 0;
    mHi = nMasses;
    nCprimes = len(cprime);
    cpLo = 0;
    cpHi = nCprimes;
    nBRnews = len(BRnew);
    brLo = 0;
    brHi = nBRnews;
    
    if options.massPoint > 0:   
        curIndex = mass.index(options.massPoint);
        mLo = curIndex;
        mHi = curIndex+1;
    if options.cPrime > 0:   
        curIndex = cprime.index(options.cPrime);
        cpLo = curIndex;
        cpHi = curIndex+1;
        nCprimes = 1;
    if options.brNew >= 0:   
        curIndex = BRnew.index(options.brNew);
        brLo = curIndex;
        brHi = curIndex+1;
        nBRnews = 1;    


    if options.injectSingalStrenght == 0 : rMin = -10 ; rMax = 10 ;
    else: rMin = - 50 ; rMax = 50

    #######################################################
    # ======= make the datacards running the fit =======  #
    #######################################################
    
    if options.makeCards:
        if not os.path.isdir("log"): os.system("mkdir log" );
        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                for k in range(brLo,brHi):

                    print "--------------------------------------------------";                
                    print "--------------------------------------------------";                
                    print "R U N N I N G   F I T S" 
                    print "mass = ",mass[i],", cprime = ",cprime[j],", brnew = ",BRnew[k],", channel: ",options.channel," pseudodata ",options.pseudodata
                    print "--------------------------------------------------";                
                    print "--------------------------------------------------";  
                    
                    time.sleep(0.3);
                    
                    command = "python doFit_class_run2exo.py %s BulkGraviton_newxsec%03d %02d %02d %02d %02d %02d %02d %s %s -b -s --cprime %01d --BRnew %01d --inPath %s/ --jetBin %s --channel %s --pseudodata %d --closuretest %d --skipJetSystematics %d --interferenceModel %s -f %s --category %s"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k], os.getcwd(), options.jetBin, options.channel,options.pseudodata,options.closuretest,options.skipJetSystematics,options.interferenceModel,options.treeFolder,options.category);
                    print command ;

                    if options.batchMode :
                        fn = "fitScript_%s_%03d_HP_%s"%(options.channel,mass[i],shape[i]);
                        submitBatchJob( command, fn );
                    if not options.batchMode: 
                        os.system(command);

    #########################################
    # ===================================== #
    #########################################

    if options.computeLimits:

        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                for k in range(brLo,brHi):

                    print "--------------------------------------------------";
                    print "analyzing card: wwlvj_BulkGraviton_newxsec%03d_em%s_HP_unbin.txt"%(mass[i],SIGCH);
                    print "--------------------------------------------------";                

                    ###############################################
                    ####### Combine single channel cards ##########
                    ###############################################


                    if options.channel =="em" and options.jetBin != "_2jet" :
                     combineCmmd = "combineCards.py wwlvj_BulkGraviton_newxsec%03d_el%s_HP_unbin.txt wwlvj_BulkGraviton_newxsec%03d_mu%s_HP_unbin.txt > wwlvj_BulkGraviton_newxsec%03d_em%s_HP_unbin.txt"%(mass[i],SIGCH,mass[i],SIGCH,mass[i],SIGCH);
                     print "combineCmmd ",combineCmmd; 
                     os.system(combineCmmd);

                    if options.higgsCombination == 1:
                     options.channel = "combo" ; 
                     combineCmmd = "combineCards.py wwlvj_BulkGraviton_newxsec%03d_mu%s_HP_unbin.txt wwlvj_BulkGraviton_newxsec%03d_el%s_HP_unbin.txt wwlvj_BulkGraviton_newxsec%03d_em_2jet%s_HP_unbin.txt > wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt"%(mass[i],SIGCH,mass[i],SIGCH,mass[i],SIGCH,mass[i],options.channel,SIGCH);
                     print "combineCmmd ",combineCmmd; 
                     os.system(combineCmmd);
                        

                    ########################################
                    ####### Generate only options ##########
                    ########################################
                    ### ORIGINAL
                    ## os.system("mkdir "+options.inputGeneratedDataset);
                    #### MATTEO CHANGED
                    if options.inputGeneratedDataset:

                       os.system("mkdir "+options.inputGeneratedDataset);
            
                       print "\n\n----------------- MATTEO CHECK------------\n\n"
                       print "OS COMMAND"
                       print options.inputGeneratedDataset
                       print "\n--------------------------------------------\n\n"
                     
                    if options.generateOnly == 1 :
                      if options.outputTree == 0 :
                          for iToy in range(options.nToys):                              
                           if options.systematics == 1:
                               runCmmd =  "combine -M GenerateOnly --saveToys -s -1 -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -t 1 --expectSignal=%f "%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.injectSingalStrenght);
                           else:
                               runCmmd =  "combine -M GenerateOnly --saveToys -s -1 -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -t 1 -S 0 --expectSignal=%f "%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.injectSingalStrenght);
                           print "runCmmd ",runCmmd;
                           if options.batchMode:
                              fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,iToy);
                              cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                           else: 
                              os.system(runCmmd);
                              os.system("mv higgsCombine* "+options.inputGeneratedDataset);   
                           continue ;

                      else:
                           runCmmd =  "combine -M GenerateOnly --saveToys -s -1 -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -t %d --expectSignal=%f  --saveNormalizations --saveWithUncertainties"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.nToys,options.injectSingalStrenght);
                           print "runCmmd ",runCmmd;
                           if options.batchMode:
                              fn = "combineScript_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                              cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                           else: 
                              os.system(runCmmd);
                              os.system("mv higgsCombine* "+options.inputGeneratedDataset);   
                           continue ;

                    ######################################
                    ####### Maximum Likelihood Fits ######
                    ######################################   
    
                    elif options.maximumLikelihoodFit == 1:
                        
                       ################################################# 
                       #### just one fit with the defined datacards ####
                       #################################################
                        
                       if options.nToys == 0 and options.crossedToys == 0 : 
                        runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties  -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2  --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);                     
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                           fn = "combineScript_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                           submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                        else:   
                         os.system(runCmmd);
                         
                       ######################################################## 
                       #### run many toys and fits with the same datacards  ###
                       ########################################################

                       elif options.nToys != 0 and options.crossedToys == 0 :
                          if options.outputTree == 0:  
                           for iToy in range(options.nToys):
                             if options.systematics == 1:  
                                 runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveToys --saveWithUncertainties --toysNoSystematics -s -1 -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin_%d -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -t 1 --expectSignal=%f --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,iToy,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.injectSingalStrenght);                     
                             else:  
                                 runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveToys --saveWithUncertainties --toysNoSystematics -s -1 -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin_%d -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -t 1 --expectSignal=%f --robustFit=1 -S 0 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,iToy,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.injectSingalStrenght);                     
                             print "runCmmd ",runCmmd;
                             if options.batchMode:
                              fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,iToy);
                              cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                             else: 
                              os.system(runCmmd);
                           continue ;
                          else:
                             runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties  --toysNoSystematics --saveToys -s -1 -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -t %d --expectSignal=%f --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.nToys,options.injectSingalStrenght);                     
                             print "runCmmd ",runCmmd;
                             if options.batchMode:
                              fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,options.nToys);
                              cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                             else: 
                              os.system(runCmmd);
                             continue ;

                       ##################################
                       ###### Make the crossed toys ##### 
                       ##################################  

                       elif options.nToys != 0 and options.crossedToys == 1  and not options.asymptotic==1:

                          command = "ls "+options.inputGeneratedDataset+" | grep root | grep higgsCombine | grep GenerateOnly | grep BulkGraviton_newxsec"+str(mass[i])+" > list_temp.txt";
#                          command = "ls "+os.getcwd()+" | grep root | grep higgsCombine | grep BulkGraviton_newxsec"+str(mass[i])+" > list_temp.txt";
                          os.system(command);
                          print command;
#                          os.system("ls "+options.inputGeneratedDataset+" | grep root | grep higgsCombine | grep BulkGraviton_newxsec"+str(mass[i])+" > list_temp.txt"); 
                          iToy = 0 ;
                          if options.outputTree == 0:  
                           with open("list_temp.txt") as input_list:
                            print "qui";   
                            for line in input_list:
                             print line;   
                             for name in line.split():
                                if iToy >= options.nToys: continue ; 
                                runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin_%d -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -s -1 -t 1  --toysFile %s/%s --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,iToy,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.inputGeneratedDataset,name);
                                iToy = iToy + 1 ;
                                print "runCmmd ",runCmmd;                                
                                if options.batchMode:
                                  fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,iToy);
                                  cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                                  submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                                else: 
                                  os.system(runCmmd);

                          else:
                           with open("list_temp.txt") as input_list:
                              for line in input_list:
                               for name in line.split():                                                                       
                                runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -s -1 -t %d --toysFile %s/%s --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.nToys,options.inputGeneratedDataset,name);                     
                                print "runCmmd ",runCmmd;                                
                                if options.batchMode:
                                  fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,iToy);
                                  cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                                  submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                                else: 
                                  os.system(runCmmd);

#                          os.system("rm list_temp.txt")
                          continue ;



                       elif options.nToys != 0 and options.crossedToys == 1 and options.asymptotic==1:

                          command = "ls "+options.inputGeneratedDataset+" | grep root | grep higgsCombine | grep GenerateOnly | grep BulkGraviton_newxsec"+str(mass[i])+" > list_temp.txt";
#                          command = "ls "+os.getcwd()+" | grep root | grep higgsCombine | grep BulkGraviton_newxsec"+str(mass[i])+" > list_temp.txt";
                          os.system(command);
                          print command;
#                          os.system("ls "+options.inputGeneratedDataset+" | grep root | grep higgsCombine | grep BulkGraviton_newxsec"+str(mass[i])+" > list_temp.txt"); 
                          iToy = 0 ;
                          if options.outputTree == 0:  
                           with open("list_temp.txt") as input_list:
                            print "qui";   
                            for line in input_list:
                             print line;   
                             for name in line.split():
                                if iToy >= options.nToys: continue ; 
                                runCmmd =  "combine -M Asymptotic --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin_%d -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -s -1 -t 1  --toysFile %s/%s"%(rMin,rMax,mass[i],options.channel,SIGCH,iToy,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.inputGeneratedDataset,name);
                                iToy = iToy + 1 ;
                                print "runCmmd ",runCmmd;                                
                                if options.batchMode:
                                  fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,iToy);
                                  cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                                  submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                                else: 
                                  os.system(runCmmd);

                          else:
                           with open("list_temp.txt") as input_list:
                              for line in input_list:
                               for name in line.split():                                                                       
                                runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -s -1 -t %d --toysFile %s/%s --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts,options.nToys,options.inputGeneratedDataset,name);                     
                                print "runCmmd ",runCmmd;                                
                                if options.batchMode:
                                  fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,iToy);
                                  cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                                  submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                                else: 
                                  os.system(runCmmd);

#                          os.system("rm list_temp.txt")
                          continue ;

                    #full CLs part    

                    elif options.systematics == 1 and not options.computePvalue == 1 and not options.computeSignif == 1 and not options.makeLikelihoodScan == 1 and options.fullCLs == 1:

                       for p in range(len(points)):
                          point = points[p];
                          runCmmd = "combine -M HybridNew --frequentist --clsAcc 0 -T 50 -i 30 -s -1 --saveHybridResult --saveToys --singlePoint %0.2f -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2"%(point,mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                          print "runCmmd ",runCmmd ;

                          if options.batchMode:
                             fn = "combineScript_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                             cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                             submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                          else: 
                             os.system(runCmmd);

                    ###############################
                    #### Asymptotic Limit part  ###
                    ###############################

                    elif options.systematics == 1 and not options.computePvalue == 1 and not options.computeSignif == 1 and not options.makeLikelihoodScan == 1:
                       runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                       print "runCmmd ",runCmmd ;
                       print "\n\n---------- MATTEO CHECK -------------\n\n"

                       if options.batchMode:
                        fn = "combineScript_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                        cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                        submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                       else: 
                        os.system(runCmmd);

                    elif options.systematics == 0 and not options.computePvalue == 1 and not options.computeSignif == 1 and not options.makeLikelihoodScan == 1:
                       runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -S 0"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                       print "runCmmd ",runCmmd ;

                       if options.batchMode:
                        fn = "combineScript_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                        cardStem = "wwlvj_BulkGraviton_newxsec%03d_em%s_HP"%(mass[i],SIGCH);
                        submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                       else: 
                        os.system(runCmmd);
                      

                    elif not options.computePvalue == 1 and not options.computeSignif == 1 and not options.makeLikelihoodScan == 1:

                       #############################################
                       ###### run Asymptotic on the final card ##### 
                       #############################################  

                       if options.nToys == 0: 
                        runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);                                        
                        print "runCmmd ",runCmmd;

                        if options.batchMode:
                         fn = "combineScript_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                         submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                        else: 
                         os.system(runCmmd);

                       else:

                       #############################################
                       ###### run Asymptotic making many toys  ##### 
                       #############################################  

                        for iToy in range(options.nToys):   
                         runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin_%d -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 -t 1 -s -1"%(mass[i],options.channel,SIGCH,iToy,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);                                        
                         print "runCmmd ",runCmmd;

                         if options.batchMode:
                          fn = "combineScript_%s_%03d%s_HP_iToy%d"%(options.channel,mass[i],SIGCH,iToy);
                          submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                         else: 
                          os.system(runCmmd);
                           
                      
                    elif options.computePvalue == 1 and not options.makeLikelihoodScan == 1: 

                       ##################################################
                       ###### run the observed and expected pvalue  ##### 
                       ##################################################  

                        runCmmd = "combine -M ProfileLikelihood --signif --pvalue -n wwlvj_pval_obs_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                        print "runCmmd ",runCmmd;

                        if options.batchMode:
                         fn = "combineScript_ProfileLikelihood_obs_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                         submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                        else:
                         os.system(runCmmd);

                        if options.nToys==-1: 
                            if options.systematics==1:
                                runCmmd = "combine -M ProfileLikelihood --signif --pvalue -n wwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 --toysFreq -t -1 --expectSignal=1 "%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                            else:
                                runCmmd = "combine -M ProfileLikelihood --signif --pvalue -n wwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 --toysFreq -t -1 -S 0 --expectSignal=1 "%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                            os.system(runCmmd);
                        else:        
                         for iToy in range(options.nToys):
                          runCmmd = "combine -M ProfileLikelihood --signif --pvalue -n wwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin_%d -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 --toysFreq -t 1 --expectSignal=1 -s -1 "%(mass[i],options.channel,SIGCH,iToy,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                          print "runCmmd ",runCmmd;

                          if options.batchMode:
                           fn = "combineScript_ProfileLikelihood_exp_%s_%03d%s_HP_%d"%(options.channel,mass[i],SIGCH,iToy);
                           submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                          else:
                           os.system(runCmmd);
                        print "runCmmd ",runCmmd;


                    elif options.computeSignif == 1 and not options.makeLikelihoodScan == 1: 
                        runCmmd = "combine -M ProfileLikelihood --signif -n wwlvj_pval_obs_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                        print "runCmmd ",runCmmd;

                        if options.batchMode:
                         fn = "combineScript_ProfileLikelihood_obs_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                         submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                        else:
                         os.system(runCmmd);

                        if options.nToys==-1: 
                            if options.systematics==1:
                                runCmmd = "combine -M ProfileLikelihood --signif -n wwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 --toysFreq -t -1 --expectSignal=1 "%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                            else:
                                runCmmd = "combine -M ProfileLikelihood --signif -n wwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 --toysFreq -t -1 -S 0 --expectSignal=1"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                            os.system(runCmmd);
                        else:        
                         for iToy in range(options.nToys):
                          runCmmd = "combine -M ProfileLikelihood --signif -n wwlvj_pval_exp_BulkGraviton_newxsec%03d_%s%s_HP_unbin_%d -m %03d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt %s -v 2 --toysFreq -t 1 --expectSignal=1 -s -1 "%(mass[i],options.channel,SIGCH,iToy,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                          print "runCmmd ",runCmmd;

                          if options.batchMode:
                           fn = "combineScript_ProfileLikelihood_exp_%s_%03d%s_HP_%d"%(options.channel,mass[i],SIGCH,iToy);
                           submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                          else:
                           os.system(runCmmd);
                        print "runCmmd ",runCmmd;


                    elif options.makeLikelihoodScan == 1:
                        
                         runCmmd = "combine -M MultiDimFit -n wwlvj_LikelihoodScan_BulkGraviton_newxsec%03d_%s%s_HP_unbin -m %03d -d wwlvj_BulkGraviton_newxsec%03d_%s%s_HP_unbin.txt  --algo=grid --points=150 --setPhysicsModelParameterRanges r=-1,5 %s"%(mass[i],options.channel,SIGCH,mass[i],mass[i],options.channel,SIGCH,moreCombineOpts);
                         print "runCmmd ",runCmmd;
 
                         if options.batchMode:
                          fn = "combineScript_LikelihoodScan_%s_%03d%s_HP"%(options.channel,mass[i],SIGCH);
                          submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                         else:
                          os.system(runCmmd);
                                                                                                                               
    ####################################################
    # =================== Bias Analysis ============== #
    ####################################################
    
    if options.biasStudy:

        for i in range(mLo,mHi): 
            print "--------------------------------------------------";                
            print "--------------------------------------------------";                
            print "B I A S  S T U D Y   F I T S" 
            print "mass = ",mass[i]," channel: ",options.channel," pseudodata ",options.pseudodata
            print "--------------------------------------------------";                
            print "--------------------------------------------------";  

            if options.jetBin == "_2jet" :
             command = "python BiasStudy/do_fitBias_higgs.py BulkGraviton_newxsec%d %d %d %d %d -b --pseudodata %d --fgen %s --fres %s --nexp %d --isMC %d --storeplot %d --channel %s --inPath %s --ttbarcontrolregion %d --fitjetmass %d --mlvjregion %s --onlybackgroundfit %d --inflatejobstatistic %d --scalesignalwidth %0.1f --injectSingalStrenght %0.1f --jetBin %s"%(mass[i],mlo[i],mhi[i],mjlo[i],mjhi[i],options.pseudodata,shape_gen[i],shape_fit[i],options.nToys,isMC[i],1,options.channel,os.getcwd(),options.ttbarcontrolregion,options.fitjetmass,options.mlvjregion,options.onlybackgroundfit,options.inflatejobstatistic,options.scalesignalwidth,options.injectSingalStrenght,options.jetBin);
            else:
             command = "python BiasStudy/do_fitBias_higgs.py BulkGraviton_newxsec%d %d %d %d %d -b --pseudodata %d --fgen %s --fres %s --nexp %d --isMC %d --storeplot %d --channel %s --inPath %s --ttbarcontrolregion %d --fitjetmass %d --mlvjregion %s --onlybackgroundfit %d --inflatejobstatistic %d --scalesignalwidth %0.1f --injectSingalStrenght %0.1f "%(mass[i],mlo[i],mhi[i],mjlo[i],mjhi[i],options.pseudodata,shape_gen[i],shape_fit[i],options.nToys,isMC[i],1,options.channel,os.getcwd(),options.ttbarcontrolregion,options.fitjetmass,options.mlvjregion,options.onlybackgroundfit,options.inflatejobstatistic,options.scalesignalwidth,options.injectSingalStrenght);
                
            print command ;
            if options.batchMode:
             suffix = options.channel;
             if options.jetBin == "_2jet" : suffix = suffix+"_2jet";
             if options.ttbarcontrolregion : suffix = suffix+"_ttbar";
             if options.fitjetmass         : suffix = suffix+"_jetmass";
             if options.turnOnAnalysis     : suffix = suffix+"_turnOn";
             if options.scalesignalwidth !=1 : suffix = suffix+ ("_width_%0.1f")%(options.scalesignalwidth);
             if options.injectSingalStrenght !=0 : suffix = suffix+"_SB";
             else :suffix = suffix+"_B";
             if options.onlybackgroundfit  : suffix = suffix+"_B";
             else: suffix = suffix+"_SB";
             fn = "biasScript_BulkGraviton_newxsec%03d_%s_%s%s"%(mass[i],shape_gen[i],shape_fit[i],suffix);
             submitBatchJob( command, fn );
            else: 
             os.system(command);

    #############################################################
    # =================== Plot of the Limit  ================== #
    #############################################################
    
    if options.plotLimits:

      if options.makeSMLimitPlot == 1:
          makeSMLimitPlot(SIGCH,10,00);
          makeSMXsecPlot(SIGCH,10,00);
          '''
          makeSMLimitPlot(SIGCH,01,00);
          makeSMLimitPlot(SIGCH,01,01);
          makeSMLimitPlot(SIGCH,01,02);
          makeSMLimitPlot(SIGCH,01,03);
          makeSMLimitPlot(SIGCH,01,04);
          makeSMLimitPlot(SIGCH,01,05);
          makeSMLimitPlot(SIGCH,02,00);
          makeSMLimitPlot(SIGCH,02,01);
          makeSMLimitPlot(SIGCH,02,02);
          makeSMLimitPlot(SIGCH,02,03);
          makeSMLimitPlot(SIGCH,02,04);
          makeSMLimitPlot(SIGCH,02,05);
          makeSMLimitPlot(SIGCH,03,00);
          makeSMLimitPlot(SIGCH,03,01);
          makeSMLimitPlot(SIGCH,03,02);
          makeSMLimitPlot(SIGCH,03,03);
          makeSMLimitPlot(SIGCH,03,04);
          makeSMLimitPlot(SIGCH,03,05);
          makeSMLimitPlot(SIGCH,05,00);
          makeSMLimitPlot(SIGCH,05,01);
          makeSMLimitPlot(SIGCH,05,02);
          makeSMLimitPlot(SIGCH,05,03);
          makeSMLimitPlot(SIGCH,05,04);
          makeSMLimitPlot(SIGCH,05,05);
          makeSMLimitPlot(SIGCH,07,00);
          makeSMLimitPlot(SIGCH,07,01);
          makeSMLimitPlot(SIGCH,07,02);
          makeSMLimitPlot(SIGCH,07,03);
          '''

      if options.plotPValue == 1:
              makeSMPValuePlot(SIGCH,10,00);
              '''
              makeSMPValuePlot(SIGCH,10,01);
              makeSMPValuePlot(SIGCH,10,02);
              makeSMPValuePlot(SIGCH,10,03);
              makeSMPValuePlot(SIGCH,10,04);
              makeSMPValuePlot(SIGCH,10,05);

              makeSMPValuePlot(SIGCH,01,00);
              makeSMPValuePlot(SIGCH,01,01);
              makeSMPValuePlot(SIGCH,01,02);
              makeSMPValuePlot(SIGCH,01,03);
              makeSMPValuePlot(SIGCH,01,04);
              makeSMPValuePlot(SIGCH,01,05);
              makeSMPValuePlot(SIGCH,02,00);
              makeSMPValuePlot(SIGCH,02,01);
              makeSMPValuePlot(SIGCH,02,02);
              makeSMPValuePlot(SIGCH,02,03);
              makeSMPValuePlot(SIGCH,02,04);
              makeSMPValuePlot(SIGCH,02,05);
              makeSMPValuePlot(SIGCH,03,00);
              makeSMPValuePlot(SIGCH,03,01);
              makeSMPValuePlot(SIGCH,03,02);
              makeSMPValuePlot(SIGCH,03,03);
              makeSMPValuePlot(SIGCH,03,04);
              makeSMPValuePlot(SIGCH,03,05);
              makeSMPValuePlot(SIGCH,05,00);
              makeSMPValuePlot(SIGCH,05,01);
              makeSMPValuePlot(SIGCH,05,02);
              makeSMPValuePlot(SIGCH,05,03);
              makeSMPValuePlot(SIGCH,05,04);
              makeSMPValuePlot(SIGCH,05,05);
              makeSMPValuePlot(SIGCH,07,00);
              makeSMPValuePlot(SIGCH,07,01);
              makeSMPValuePlot(SIGCH,07,02);
              makeSMPValuePlot(SIGCH,07,03);
              '''

      if options.makeBSMLimitPlotMass == 1:
          makeBSMLimitPlotMass(SIGCH,00);
          makeBSMLimitPlotMass(SIGCH,01);
          makeBSMLimitPlotMass(SIGCH,02);
          makeBSMLimitPlotMass(SIGCH,03);
          makeBSMLimitPlotMass(SIGCH,04);
          makeBSMLimitPlotMass(SIGCH,05);
                 
      if options.makeBSMLimitPlotBRnew == 1:
          
          makeBSMLimitPlotBRnew(SIGCH,600);
          makeBSMLimitPlotBRnew(SIGCH,700);
          makeBSMLimitPlotBRnew(SIGCH,800);
          makeBSMLimitPlotBRnew(SIGCH,900);
          makeBSMLimitPlotBRnew(SIGCH,1000);

      if options.makeBSMLimitPlot2D == 1:

          contourListMassExp = defaultdict(list);          
          contourListMassObs = defaultdict(list);          

          makeBSMLimitPlot2D(SIGCH,600,contourListMassExp,contourListMassObs);
          makeBSMLimitPlot2D(SIGCH,700,contourListMassExp,contourListMassObs);
          makeBSMLimitPlot2D(SIGCH,800,contourListMassExp,contourListMassObs);
          makeBSMLimitPlot2D(SIGCH,900,contourListMassExp,contourListMassObs);
          makeBSMLimitPlot2D(SIGCH,1000,contourListMassExp,contourListMassObs);

          contourListBrNewExp = defaultdict(list);          
          contourListBrNewObs = defaultdict(list);          

          makeBSMLimitPlot2DBRnew(SIGCH,00,contourListBrNewExp,contourListBrNewObs);
          makeBSMLimitPlot2DBRnew(SIGCH,01,contourListBrNewExp,contourListBrNewObs);
          makeBSMLimitPlot2DBRnew(SIGCH,02,contourListBrNewExp,contourListBrNewObs);
          makeBSMLimitPlot2DBRnew(SIGCH,03,contourListBrNewExp,contourListBrNewObs);
          makeBSMLimitPlot2DBRnew(SIGCH,04,contourListBrNewExp,contourListBrNewObs);
          makeBSMLimitPlot2DBRnew(SIGCH,05,contourListBrNewExp,contourListBrNewObs);

          makeContourPlotMass(contourListMassExp,contourListMassObs);
          makeContourPlotBrNew(contourListBrNewExp,contourListBrNewObs);

      if options.plotSignalStrenght == 1:
          makeSignalStrenghtPlot(SIGCH,10,00);

      if options.plotLikelihoodScan == 1:
          makeLikelihoodScanPlot(SIGCH,10,00);
          

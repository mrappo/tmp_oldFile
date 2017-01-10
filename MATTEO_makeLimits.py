import os,commands
import sys
from optparse import OptionParser
import subprocess

parser = OptionParser()

parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--ntuple', action="store",type="string",dest="ntuple",default="WWTree_22sep_jecV7_lowmass")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")
#parser.add_option('--type', action="store",type="string",dest="type",default="")
#parser.add_option('--jetalgo', action="store",type="string",dest="jetalgo",default="jet_mass_pr")
#parser.add_option('--interpolate', action="store_true",dest="interpolate",default=False)
#parser.add_option('--batchMode', action="store_true",dest="batchMode",default=False)
parser.add_option('--vbf', action="store_true",dest="VBF_process",default=False)
parser.add_option('--pseudodata', action="store_true",dest="pseudodata",default=False)
parser.add_option('--copyDC', action="store_true",dest="copyDC",default=False)
(options, args) = parser.parse_args()

currentDir = os.getcwd();


Ntuple_dir_name="Ntuple_%s"%(options.ntuple)

if not os.path.isdir(Ntuple_dir_name):
       os.system("mkdir "+Ntuple_dir_name);

#Ntuple_Path_lxplus="/afs/cern.ch/user/l/lbrianza/work/public/%s/"%options.ntuple
samples=["BulkGraviton","Higgs"]
lumi_float_true=2197.96;
luminosities=[lumi_float_true]


if options.VBF_process:
   tmp_vbf_name="_VBF_";
   
else:
   tmp_vbf_name="_";


channel_in=options.channel
   
      
      
for lumi_float_value in luminosities:
    
    
    
    lumi_str=str("%.0f"%lumi_float_value);
    
    
       
    if options.pseudodata:
          
          pseudodata_dir=Ntuple_dir_name+"/pseudoData"
          if not os.path.isdir(pseudodata_dir):
                 os.system("mkdir "+pseudodata_dir);
                     
          if options.VBF_process: 
      
             datacards_dir_in="../../../CMSSW_5_3_13/src/EXOVVFitter_mr/pseudoData/Ntuple_%s/Lumi_%s_VBF/cards_%s_%s"%(options.ntuple,lumi_str,options.channel,options.category,sample);
             lumi_dir=Ntuple_dir_name+"/pseudoData/Lumi_%s_VBF"%lumi_str;
             
   
          else:
             datacards_dir_in="../../../CMSSW_5_3_13/src/EXOVVFitter_mr/Ntuple_%s/pseudoData/Lumi_%s/cards_%s_%s"%(options.ntuple,lumi_str,options.channel,options.category);
             lumi_dir=Ntuple_dir_name+"/pseudoData/Lumi_%s"%lumi_str;
             





    else:

          truedata_dir=Ntuple_dir_name+"/trueData"
          if not os.path.isdir(truedata_dir):
                 os.system("mkdir "+truedata_dir);
          
          if options.VBF_process:
             datacards_dir_in="../../../CMSSW_5_3_13/src/EXOVVFitter_mr/Ntuple_%s/trueData/Lumi_%s_VBF/cards_%s_%s"%(options.ntuple,lumi_str,options.channel,options.category);
             lumi_dir=Ntuple_dir_name+"/trueData/Lumi_%s_VBF"%lumi_str;
             
          else:
             datacards_dir_in="../../../CMSSW_5_3_13/src/EXOVVFitter_mr/Ntuple_%s/trueData/Lumi_%s/cards_%s_%s"%(options.ntuple,lumi_str,options.channel,options.category);
             lumi_dir=Ntuple_dir_name+"/trueData/Lumi_%s"%lumi_str;
             
    
    
    datacards_dir_out=lumi_dir
    if not options.copyDC:
           if not os.path.isdir(datacards_dir_out):
                  print datacards_dir_in  
                  print "\n\nDatacard Mancanti. Uscita. \n\n"
                  exit();
    
    if options.copyDC:
       if not os.path.isdir(datacards_dir_out):
              os.system("mkdir "+datacards_dir_out);
   
       p1 = subprocess.Popen(['cp','-r',datacards_dir_in,datacards_dir_out])
       p1.wait()
    
    
    for sample in samples:
        
        if sample.find('BulkGraviton') !=-1:
           masses=[600,800,1000]
       
       
        if sample.find('Higgs') !=-1:
           masses=[650,1000]
        
        for m in masses:
            mass=str(m);
            '''
            if (options.ntuple=="WWTree_18feb_jecV7_lowmass" and sample=="BulkGraviton"):
                   datacard_file_in=datacards_dir_out+"/cards_%s_%s/%s/wwlvj%s%s_newxsec%s_%s_%s_lumi_%s_unbin.txt"%(options.channel,options.category,sample,tmp_vbf_name,sample,mass,options.channel,options.category,lumi_str)
                   datacard_file_out=datacards_dir_out+"/cards_%s_%s/%s/wwlvj_BulkGraviton_newxsec%s_%s_HP_unbin.txt"%(options.channel,options.category,sample,mass,channel_in)
                   p2 = subprocess.Popen(['cp',datacard_file_in,datacard_file_out])
                   p2.wait()
            
            else:
            '''
            datacard_file_in=datacards_dir_out+"/cards_%s_%s/%s/wwlvj%s%s%s_%s_%s_lumi_%s_unbin.txt"%(options.channel,options.category,sample,tmp_vbf_name,sample,mass,options.channel,options.category,lumi_str)
            datacard_file_out=datacards_dir_out+"/cards_%s_%s/%s/wwlvj_BulkGraviton_newxsec%s_%s_HP_unbin.txt"%(options.channel,options.category,sample,mass,channel_in)
            p2 = subprocess.Popen(['cp',datacard_file_in,datacard_file_out])
            p2.wait()
            
            
        
        if options.VBF_process:
           
           cards_dir=datacards_dir_out+"/cards_%s_%s/%s/"%(options.channel,options.category,sample)
           cmd="python MATTEO_runLimits.py -b --computeLimits --channel %s --datacardDIR %s --makeSMLimitPlot 1 --plotLimits 1 --systematics 1 --sample %s --vbf TRUE "%(options.channel,cards_dir,sample)
        
           #p2 = subprocess.Popen([cmd])
           #p2.wait()
           print cmd
           os.system(cmd)
        
        else:
           cards_dir=datacards_dir_out+"/cards_%s_%s/%s/"%(options.channel,options.category,sample)
           cmd="python MATTEO_runLimits.py -b --computeLimits --channel %s --datacardDIR %s --makeSMLimitPlot 1 --plotLimits 1 --systematics 1 --sample %s"%(options.channel,cards_dir,sample)
        
           #p2 = subprocess.Popen([cmd])
           #p2.wait()
           print cmd
           os.system(cmd)
        
    
    
    
    
    
  
    
  
    
    
    
    
    
       

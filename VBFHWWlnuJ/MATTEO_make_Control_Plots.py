import os,commands
import sys
from optparse import OptionParser
import subprocess

parser = OptionParser()

parser.add_option('--channel', action="store", type="string", dest="channel", default="mu")
#parser.add_option('--sample', action="store", type="string", dest="sample", default="BulkGraviton")
parser.add_option('--ntuple', action="store", type="string", dest="ntuple", default="22sep")
parser.add_option('--lumi', action="store", type="string", dest="lumi", default="2300")
parser.add_option('--mass', action="store", type="string", dest="mass", default="2300")
parser.add_option('--sample', action="store", type="string", dest="sample", default="BulkGraviton")
(options, args) = parser.parse_args()
currentDir = os.getcwd();
        

process_name=options.sample+options.mass


    
cfg_file_BG="cfg/DataMCComparison_InputCfgFile/Run2_DataMCComparisonPlot_76x_mu_22sep_%s.cfg"%process_name
p1 = subprocess.Popen(['./bin/DataMCComparisonPlot.exe',cfg_file_BG]) 
p1.wait()
        
 
        
path_dir_in="output/run2/MCDATAComparisonPlot_mu_22sep_%s/Run2_MCDataComparisonRSGraviton2000_mu_plot"%process_name
path_dir_out="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s"%(options.ntuple,options.lumi,options.channel,process_name)
p2 = subprocess.Popen(['cp','-r',path_dir_in,path_dir_out])
p2.wait()



data_in="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Run2_MCDataComparisonRSGraviton2000_mu_plot"%(options.ntuple,options.lumi,options.channel,process_name)
data_out="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Data"%(options.ntuple,options.lumi,options.channel,process_name)
p3 = subprocess.Popen(['mv',data_in,data_out])
p3.wait()

root_in="output/run2/MCDATAComparisonPlot_mu_22sep_%s/Run2_MCDataComparisonRSGraviton2000_mu.root"%process_name
root_out="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Data/Root_ControlPlots_out_%s.root"%(options.ntuple,options.lumi,options.channel,process_name,process_name)

p4 = subprocess.Popen(['cp',root_in,root_out])
p4.wait()

if os.path.isdir(data_in):
   p5=subprocess.Popen(['rm','-r',data_in])
   p5.wait()


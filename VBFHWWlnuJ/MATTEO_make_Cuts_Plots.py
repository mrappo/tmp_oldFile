import os,commands
import sys
from optparse import OptionParser
import subprocess

parser = OptionParser()

parser.add_option('--mass', action="store", type="string", dest="mass", default="600")
parser.add_option('--channel', action="store", type="string", dest="channel", default="mu")
parser.add_option('--sample', action="store", type="string", dest="sample", default="BulkGraviton")
parser.add_option('--ntuple', action="store", type="string", dest="ntuple", default="22sep")
parser.add_option('--lumi', action="store", type="string", dest="lumi", default="600")
(options, args) = parser.parse_args()

currentDir = os.getcwd();

process_name=options.sample+options.mass
cfg_file_input="cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_22sep_%s.cfg"%process_name
p1 = subprocess.Popen(['./bin/VBFOptimizeSelections.exe',cfg_file_input]) 
p1.wait()
#####os.system("./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGraviton1000_data2015.cfg")

p9 = subprocess.Popen(['cp','-r','output/outputTMVATraining_BulkG1000/TMVATrainingResult_BulkG1000_mu_PTBin_0_5000.root','TMVA.root'])
p9.wait()
#os.system("cp output/outputTMVATraining_BulkG1000/TMVATrainingResult_BulkG1000_mu_PTBin_0_5000.root TMVA.root")
#####os.system("root -l macros/TMVAMacro/mvaeffs.C")
#####p2 = subprocess.Popen(['root','-l','macros/TMVAMacro/mvaeffs.C'])
p2 = subprocess.Popen(['root','-l','macros/TMVAMacro/MATTEO_mvaeff.C'])
p2.wait()

#command1 = "--sample %s --channel %s --ntuple %s --mass %s"%(options.sample,options.channel,options.ntuple,options.mass)
p3 = subprocess.Popen(['python','MATTEO_read_file.py','--sample',options.sample,'--channel',options.channel,'--ntuple',options.ntuple,'--mass',options.mass])
p3.wait()
#os.system(command1)

plots_out_dir="output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/plots"%(options.ntuple,options.lumi,options.channel,process_name)
p4 = subprocess.Popen(['cp','-r','plots/',plots_out_dir])
p4.wait()


#os.system("cp -r plots/ output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/plots"%(options.ntuple,options.lumi,options.channel,process_name))

data_out_dir="output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Data"%(options.ntuple,options.lumi,options.channel,process_name)
p5 = subprocess.Popen(['cp','-r','output/outputTMVATraining_BulkG1000/TMVAWeight_CutsMC_mu_PTBin_0_5000/',data_out_dir])
p5.wait()
#os.system("cp -r output/outputTMVATraining_BulkG1000/TMVAWeight_CutsMC_mu_PTBin_0_5000/ output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Data"%(options.ntuple,options.lumi,options.channel,process_name))
tmva_out="output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Data/TMVA.root"%(options.ntuple,options.lumi,options.channel,process_name)
p6 = subprocess.Popen(['cp','-r','TMVA.root',tmva_out])
p6.wait()

#os.system("cp -r TMVA.root output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Data/TMVA.root"%(options.ntuple,options.lumi,options.channel,process_name))
p7 = subprocess.Popen(['rm','-r','plots/'])
p7.wait()
p8 = subprocess.Popen(['rm','-r','TMVA.root'])
p8.wait()
#os.system("rm -r plots/")
#os.system("rm -r TMVAWeight_CutsMC_mu_PTBin_0_5000/")
#os.system("rm -r TMVA.root")





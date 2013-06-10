#!/usr/bin/env python

########################################
##     
##       Author: Wei Zou
##       
##       Email: weizou.pku@gmail.com
#######################################

from ROOT import *

class Samples:
      
      def __init__(self, CHANNEL):
          
          self.filenames = {}
          self.filepath = ""
          self.luminosity = 0.
          self.treename = ""
          self.channel = CHANNEL

      def SetFilePath(self,path):

          self.filepath = path

      def SetLumi(self,lumi):

          self.luminosity = lumi
      
      def SetTreeName(self,tree):

          self.treename = tree

      def SetFileNames(self):
          if self.channel == "mu":
              self.filenames["data"] = self.filepath + "RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root"
              self.filenames["data_194712_559_405404310"]      = self.filepath + "RD_WmunuJetAnalysisntuple_194712_559_405404310.root"
              self.filenames["data_204544_102_152417542"]      = self.filepath + "RD_WmunuJetAnalysisntuple_204544_102_152417542.root"
              self.filenames["data_195930_83_62880579"]        = self.filepath + "RD_WmunuJetAnalysisntuple_195930_83_62880579.root"
              self.filenames["data_208487_118_209515269"]      = self.filepath + "RD_WmunuJetAnalysisntuple_208487_118_209515269.root"
              self.filenames["data_190703_11_8519680"]         = self.filepath + "RD_WmunuJetAnalysisntuple_190703_11_8519680.root"
              self.filenames["data_196438_586_479961054"]      = self.filepath + "RD_WmunuJetAnalysisntuple_196438_586_479961054.root"
              self.filenames["data_191834_99_120189151"]       = self.filepath + "RD_WmunuJetAnalysisntuple_191834_99_120189151.root"
              self.filenames["data_196452_193_232654780"]      = self.filepath + "RD_WmunuJetAnalysisntuple_196452_193_232654780.root"
              self.filenames["data_gr2_194199_190_149891910"]  = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_194199_190_149891910.root"
              self.filenames["data_194199_258_229967752"]      = self.filepath + "RD_WmunuJetAnalysisntuple_194199_258_229967752.root"
              self.filenames["data_198522_104_75971527"]       = self.filepath + "RD_WmunuJetAnalysisntuple_198522_104_75971527.root"
              self.filenames["data_gr2_194897_50_91029885"]    = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_194897_50_91029885.root"
              self.filenames["data_194699_80_97577634"]        = self.filepath + "RD_WmunuJetAnalysisntuple_194699_80_97577634.root"
              self.filenames["data_199428_457_557329806"]      = self.filepath + "RD_WmunuJetAnalysisntuple_199428_457_557329806.root"
              self.filenames["data_gr2_195397_558_772053348"]  = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_195397_558_772053348.root"
              self.filenames["data_199428_86_8009087"]         = self.filepath + "RD_WmunuJetAnalysisntuple_199428_86_8009087.root"
              self.filenames["data_gr2_199608_1123_1247460423"] = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_199608_1123_1247460423.root"
              self.filenames["data_194778_190_248720484"]       = self.filepath + "RD_WmunuJetAnalysisntuple_194778_190_248720484.root"
              self.filenames["data_199961_177_187360166"]       = self.filepath + "RD_WmunuJetAnalysisntuple_199961_177_187360166.root"
              self.filenames["data_gr2_200091_1678_1762413892"] = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_200091_1678_1762413892.root"
              self.filenames["data_194912_390_657997761"]       = self.filepath + "RD_WmunuJetAnalysisntuple_194912_390_657997761.root"
              self.filenames["data_200188_104_142020231"]       = self.filepath + "RD_WmunuJetAnalysisntuple_200188_104_142020231.root"
              self.filenames["data_gr2_202299_249_336447785"]   = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_202299_249_336447785.root"
              self.filenames["data_195013_303_44790264"]        = self.filepath + "RD_WmunuJetAnalysisntuple_195013_303_44790264.root"
              self.filenames["data_203002_1517_1719001770"]     = self.filepath + "RD_WmunuJetAnalysisntuple_203002_1517_1719001770.root"
              self.filenames["data_gr2_202314_128_132710041"]   = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_202314_128_132710041.root"
              self.filenames["data_195099_171_245228003"]       = self.filepath + "RD_WmunuJetAnalysisntuple_195099_171_245228003.root"
              self.filenames["data_gr2_206476_128_148892985"]   = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_206476_128_148892985.root"
              self.filenames["data_195398_899_752941172"]       = self.filepath + "RD_WmunuJetAnalysisntuple_195398_899_752941172.root"
              self.filenames["data_204564_359_400306553"]       = self.filepath + "RD_WmunuJetAnalysisntuple_204564_359_400306553.root"
              self.filenames["data_gr2_206859_764_1094226833"]  = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_206859_764_1094226833.root"
              self.filenames["data_195399_159_130828536"]       = self.filepath + "RD_WmunuJetAnalysisntuple_195399_159_130828536.root"
              self.filenames["data_207372_380_547252863"]       = self.filepath + "RD_WmunuJetAnalysisntuple_207372_380_547252863.root"
              self.filenames["data_gr2_207231_920_1260547040"]  = self.filepath + "RD_WmunuJetAnalysisntuple_gr2_207231_920_1260547040.root"
              self.filenames["data_207490_80_87321589"]         = self.filepath + "RD_WmunuJetAnalysisntuple_2DD07490_80_87321589.root"

#              self.filenames["data"] = "extraScripts/RDclone.root" 
#              self.filenames["data"] = "/eos/uscms/store/user/lnujj/Moriond2013/ReducedTrees/RD_WmunuJets_DataAll_GoldenJSON_5p3invfb.root" 
              self.filenames["TTbar"] = self.filepath + "RD_mu_TTbar_CMSSW532.root"
              self.filenames["TTbar_matchDn"] = self.filepath + "RD_mu_TTbar_matchingdown_CMSSW532.root"
              self.filenames["TTbar_matchUp"] = self.filepath + "RD_mu_TTbar_matchingup_CMSSW532.root"
              self.filenames["TTbar_Powheg"] = self.filepath + "RD_mu_TTbar_powheg_CMSSW532.root"
              self.filenames["TTbar_scaleDn"] = self.filepath + "RD_mu_TTbar_scaleup_CMSSW532.root"
              self.filenames["TTbar_scaleUp"] = self.filepath + "RD_mu_TTbar_scaledown_CMSSW532.root"
              self.filenames["WJets_Pythia"] = self.filepath + "RD_mu_WpJPt100_CMSSW532.root"
              self.filenames["WJets_Pythia180"] = self.filepath + "RD_mu_WpJPt180_CMSSW532.root"              
              self.filenames["WJets_Pythia180_higgs"] = self.filepath + "RD_mu_WpJ_PT180_CMSSW532_higgs.root"
              self.filenames["WJets_Pythia180_newid"] = self.filepath + "RD_mu_WpJ_PT180_CMSSW532_newid.root"              
#              self.filenames["WJets_Pythia"] = self.filepath + "RD_el_WpJPt100_CMSSW532.root"
              self.filenames["WJets_Herwig"] = self.filepath + "RD_mu_WpJPt100_herwig_CMSSW532.root";
              self.filenames["ZJets"] = self.filepath + "RD_mu_ZpJ_CMSSW532.root"
              self.filenames["tch"] = self.filepath + "RD_mu_STopT_T_CMSSW532.root"
              self.filenames["tWch"] = self.filepath + "RD_mu_STopTW_T_CMSSW532.root"
              self.filenames["sch"] = self.filepath + "RD_mu_STopS_T_CMSSW532.root"
              self.filenames["tch_bar"] = self.filepath + "RD_mu_STopT_Tbar_CMSSW532.root"
              self.filenames["tWch_bar"] = self.filepath + "RD_mu_STopTW_Tbar_CMSSW532.root"
              self.filenames["sch_bar"] = self.filepath + "RD_mu_STopS_Tbar_CMSSW532.root"
              self.filenames["WW"] = self.filepath + "RD_mu_WW_CMSSW532.root"
              self.filenames["WZ"] = self.filepath + "RD_mu_WZ_CMSSW532.root"
              self.filenames["ZZ"] = self.filepath + "RD_mu_ZZ_CMSSW532.root"
              self.filenames["ggH600"] = self.filepath + "RD_mu_HWWMH600_CMSSW532_private.root"
              self.filenames["ggH700"] = self.filepath + "RD_mu_HWWMH700_CMSSW532_private.root"
              self.filenames["ggH800"] = self.filepath + "RD_mu_HWWMH800_CMSSW532_private.root"
              self.filenames["ggH900"] = self.filepath + "RD_mu_HWWMH900_CMSSW532_private.root"
              self.filenames["ggH1000"] = self.filepath + "RD_mu_HWWMH1000_CMSSW532_private.root"
              self.filenames["vbfH600"] = self.filepath + "RD_mu_VBFHWWMH600_CMSSW532_private.root"
              self.filenames["vbfH700"] = self.filepath + "RD_mu_VBFHWWMH700_CMSSW532_private.root"
              self.filenames["vbfH800"] = self.filepath + "RD_mu_VBFHWWMH800_CMSSW532_private.root"
              self.filenames["vbfH900"] = self.filepath + "RD_mu_VBFHWWMH900_CMSSW532_private.root"
              self.filenames["vbfH1000"] = self.filepath + "RD_mu_VBFHWWMH1000_CMSSW532_private.root"
              self.filenames["rsg1000_kMpl01_py"] = self.filepath + "RD_mu_RSGravitonToWW_kMpl01_M-1000_pythia_CMSSW532_private.root"          
              self.filenames["rsg1000_kMpl01_hw"] = self.filepath + "RD_mu_RSGravitonToWW_kMpl01_M-1000_herwig_CMSSW532_private.root"          
              self.filenames["rsg1500_kMpl01_py"] = self.filepath + "RD_mu_RSGravitonToWW_kMpl01_M-1500_pythia_CMSSW532_private.root"          
              self.filenames["rsg1500_kMpl01_hw"] = self.filepath + "RD_mu_RSGravitonToWW_kMpl01_M-1500_herwig_CMSSW532_private.root"          
              self.filenames["rsg2000_kMpl01_py"] = self.filepath + "RD_mu_RSGravitonToWW_kMpl01_M-2000_pythia_CMSSW532_private.root"                    
              self.filenames["BulkG_c0p2_M600"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M600.root"                    
              self.filenames["BulkG_c0p2_M700"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M700.root"                    
              self.filenames["BulkG_c0p2_M800"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M800.root"                    
              self.filenames["BulkG_c0p2_M900"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M900.root"                    
              self.filenames["BulkG_c0p2_M1000"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1000.root"                    
              self.filenames["BulkG_c0p2_M1100"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1100.root"                    
              self.filenames["BulkG_c0p2_M1200"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1200.root"                    
              self.filenames["BulkG_c0p2_M1300"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1300.root"                    
              self.filenames["BulkG_c0p2_M1400"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1400.root"                    
              self.filenames["BulkG_c0p2_M1500"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1500.root"                    
              self.filenames["BulkG_c0p2_M1600"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1600.root"                    
              self.filenames["BulkG_c0p2_M1700"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1700.root"                    
              self.filenames["BulkG_c0p2_M1800"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1800.root"                    
              self.filenames["BulkG_c0p2_M1900"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M1900.root"                    
              self.filenames["BulkG_c0p2_M2000"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M2000.root"                    
              self.filenames["BulkG_c0p2_M2100"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M2100.root"                    
              self.filenames["BulkG_c0p2_M2200"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M2200.root"                    
              self.filenames["BulkG_c0p2_M2300"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M2300.root"                    
              self.filenames["BulkG_c0p2_M2400"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M2400.root"                    
              self.filenames["BulkG_c0p2_M2500"] = self.filepath + "RD_mu_BulkG_WW_lvjj_c0p2_M2500.root"                    

          elif self.channel == "el":
              self.filenames["data"] = self.filepath + "RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root"
              self.filenames["data_195013_114_117238404"]  = self.filepath + "RD_WenuJetAnalysisntuple_195013_114_117238404.root"
              self.filenames["data_195948_317_517532532"]   = self.filepath + "RD_WenuJetAnalysisntuple_195948_317_517532532.root"
              self.filenames["data_201602_497_68268345"]   = self.filepath + "RD_WenuJetAnalysisntuple_201602_497_68268345.root"
              self.filenames["data_195397_899_1123713321"] = self.filepath + "RD_WenuJetAnalysisntuple_195397_899_1123713321.root"
              self.filenames["data_201191_317_488053419"]  = self.filepath + "RD_WenuJetAnalysisntuple_201191_317_488053419.root"
              self.filenames["data_202016_935_952022882"]  = self.filepath + "RD_WenuJetAnalysisntuple_202016_935_952022882.root"
              self.filenames["data_195655_60_73510965"]    = self.filepath + "RD_WenuJetAnalysisntuple_195655_60_73510965.root"
              self.filenames["data_201278_532_72515017"]   = self.filepath + "RD_WenuJetAnalysisntuple_201278_532_72515017.root"
              self.filenames["data_202973_241_260228320"]  = self.filepath + "RD_WenuJetAnalysisntuple_202973_241_260228320.root"
              self.filenames["TTbar"] = self.filepath + "RD_el_TTbar_CMSSW532.root";
              self.filenames["TTbar_matchDn"] = self.filepath + "RD_el_TTbar_matchingdown_CMSSW532.root"
              self.filenames["TTbar_matchUp"] = self.filepath + "RD_el_TTbar_matchingup_CMSSW532.root"
              self.filenames["TTbar_Powheg"] = self.filepath + "RD_el_TTbar_powheg_CMSSW532.root"
              self.filenames["TTbar_scaleDn"] = self.filepath + "RD_el_TTbar_scaleup_CMSSW532.root"
              self.filenames["TTbar_scaleUp"] = self.filepath + "RD_el_TTbar_scaledown_CMSSW532.root"              
              self.filenames["WJets_Pythia"] = self.filepath + "RD_el_WpJPt100_CMSSW532.root";
              self.filenames["WJets_Pythia180"] = self.filepath + "RD_el_WpJ_PT180_CMSSW532.root"                            
              self.filenames["WJets_Pythia180_higgs"] = self.filepath + "RD_el_WpJ_PT180_CMSSW532_higgs.root"              
              self.filenames["WJets_Pythia180_newid"] = self.filepath + "RD_el_WpJ_PT180_CMSSW532_newid.root"              
#              self.filenames["WJets_Pythia"] = self.filepath + "RD_el_WpJPt100_CMSSW532.root";              
              self.filenames["WJets_Herwig"] = self.filepath + "RD_el_WpJPt100_herwig_CMSSW532.root";
              self.filenames["ZJets"] = self.filepath + "RD_el_ZpJ_CMSSW532.root"
              self.filenames["tch"] = self.filepath + "RD_el_STopT_T_CMSSW532.root"
              self.filenames["tWch"] = self.filepath + "RD_el_STopTW_T_CMSSW532.root"
              self.filenames["sch"] = self.filepath + "RD_el_STopS_T_CMSSW532.root"
              self.filenames["tch_bar"] = self.filepath + "RD_el_STopT_Tbar_CMSSW532.root"
              self.filenames["tWch_bar"] = self.filepath + "RD_el_STopTW_Tbar_CMSSW532.root"
              self.filenames["sch_bar"] = self.filepath + "RD_el_STopS_Tbar_CMSSW532.root"
              self.filenames["WW"] = self.filepath + "RD_el_WW_CMSSW532.root"
              self.filenames["WZ"] = self.filepath + "RD_el_WZ_CMSSW532.root"
              self.filenames["ZZ"] = self.filepath + "RD_el_ZZ_CMSSW532.root"
              self.filenames["ggH600"] = self.filepath + "RD_el_HWWMH600_CMSSW532_private.root"
              self.filenames["ggH700"] = self.filepath + "RD_el_HWWMH700_CMSSW532_private.root"
              self.filenames["ggH800"] = self.filepath + "RD_el_HWWMH800_CMSSW532_private.root"
              self.filenames["ggH900"] = self.filepath + "RD_el_HWWMH900_CMSSW532_private.root"
              self.filenames["ggH1000"] = self.filepath + "RD_el_HWWMH1000_CMSSW532_private.root"
              self.filenames["vbfH600"] = self.filepath + "RD_el_VBFHWWMH600_CMSSW532_private.root"
              self.filenames["vbfH700"] = self.filepath + "RD_el_VBFHWWMH700_CMSSW532_private.root"
              self.filenames["vbfH800"] = self.filepath + "RD_el_VBFHWWMH800_CMSSW532_private.root"
              self.filenames["vbfH900"] = self.filepath + "RD_el_VBFHWWMH900_CMSSW532_private.root"
              self.filenames["vbfH1000"] = self.filepath + "RD_el_VBFHWWMH1000_CMSSW532_private.root"
              self.filenames["rsg1000_kMpl01_py"] = self.filepath + "RD_el_RSGravitonToWW_kMpl01_M-1000_pythia_CMSSW532_private.root"          
              self.filenames["rsg1000_kMpl01_hw"] = self.filepath + "RD_el_RSGravitonToWW_kMpl01_M-1000_herwig_CMSSW532_private.root"          
              self.filenames["rsg1500_kMpl01_py"] = self.filepath + "RD_el_RSGravitonToWW_kMpl01_M-1500_pythia_CMSSW532_private.root"          
              self.filenames["rsg1500_kMpl01_hw"] = self.filepath + "RD_el_RSGravitonToWW_kMpl01_M-1500_herwig_CMSSW532_private.root"          
              self.filenames["rsg2000_kMpl01_py"] = self.filepath + "RD_el_RSGravitonToWW_kMpl01_M-2000_pythia_CMSSW532_private.root"                    

              self.filenames["BulkG_c0p2_M600"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M600.root"                    
              self.filenames["BulkG_c0p2_M700"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M700.root"                    
              self.filenames["BulkG_c0p2_M800"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M800.root"                    
              self.filenames["BulkG_c0p2_M900"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M900.root"                    
              self.filenames["BulkG_c0p2_M1000"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1000.root"                    
              self.filenames["BulkG_c0p2_M1100"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1100.root"                    
              self.filenames["BulkG_c0p2_M1200"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1200.root"                    
              self.filenames["BulkG_c0p2_M1300"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1300.root"                    
              self.filenames["BulkG_c0p2_M1400"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1400.root"                    
              self.filenames["BulkG_c0p2_M1500"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1500.root"                    
              self.filenames["BulkG_c0p2_M1600"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1600.root"                    
              self.filenames["BulkG_c0p2_M1700"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1700.root"                    
              self.filenames["BulkG_c0p2_M1800"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1800.root"                    
              self.filenames["BulkG_c0p2_M1900"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M1900.root"                    
              self.filenames["BulkG_c0p2_M2000"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M2000.root"                    
              self.filenames["BulkG_c0p2_M2100"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M2100.root"                    
              self.filenames["BulkG_c0p2_M2200"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M2200.root"                    
              self.filenames["BulkG_c0p2_M2300"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M2300.root"                    
              self.filenames["BulkG_c0p2_M2400"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M2400.root"                    
              self.filenames["BulkG_c0p2_M2500"] = self.filepath + "RD_el_BulkG_WW_lvjj_c0p2_M2500.root"                    


      def GetLumiScaleFactor(self,txtfile,keyname):
          
          multiplicitylabel = 1.0
          scalefactor = 1.0
          SFfile = open(txtfile)
          for sfline in SFfile:
              if sfline.find("#")!=-1: continue
              if(sfline.find(keyname) != -1):
                 scalefactor = float(sfline.split()[1])
                 if len(sfline.split()) > 2:
                    multiplicitylabel = float(sfline.split()[2])
                 scalefactor = scalefactor * multiplicitylabel
                 break 
          SFfile.close()
          return scalefactor

      def GetFileNames(self):
         
          return self.filenames

      def GetLumi(self):
          
          return self.luminosity

      def GetTreeName(self):
          
          return self.treename
      

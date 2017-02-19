import string
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

f1 = open('MATTEO_bin_file.txt', 'r')
bin = f1.readline()
print bin
f1.close()

f2 = open('output/outputTMVATraining_BulkG1000/TMVAWeight_CutsMC_mu_PTBin_0_5000/otree_mu_PTBin_0_5000_CutsMC.weights.xml','r')
counter=0
for line in f2:
    counter = counter +1
    if line.find(bin) != -1:
       #print line
       #print counter
       nline=int(counter)
f2.close()
intestation="%s\t\t%s%s\t%s\n"%(options.ntuple,options.sample,options.mass,options.channel)
print"\n\n Ricerca Limiti: \n\n"
print intestation
f3= open('output/outputTMVATraining_BulkG1000/TMVAWeight_CutsMC_mu_PTBin_0_5000/otree_mu_PTBin_0_5000_CutsMC.weights.xml','r')
filename= "MATTEO_Cut_result_%s.txt"%options.ntuple
f4= open(filename,'a')
#print
lines=[nline-1,nline]
i=0

f4.write(intestation)
for lin in f3:
    
    if i in lines:
       print lin
       f4.write(lin)
    i+=1
f4.write("\n\n")
f3.close()
f4.close()


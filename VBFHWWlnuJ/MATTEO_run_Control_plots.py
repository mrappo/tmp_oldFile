import os,commands
import sys
from optparse import OptionParser
import subprocess

parser = OptionParser()

parser.add_option('--channel', action="store", type="string", dest="channel", default="mu")
#parser.add_option('--sample', action="store", type="string", dest="sample", default="BulkGraviton")
parser.add_option('--ntuple', action="store", type="string", dest="ntuple", default="22sep")
#parser.add_option('--lumi', action="store", type="float", dest="lumi", default="2300")
(options, args) = parser.parse_args()
currentDir = os.getcwd();

sample=["BulkGraviton","Higgs"]


lumi_true_value=2248.51

luminosity=[lumi_true_value]
for lumi in luminosity:


    for s in sample:

        if s.find('BulkGraviton') !=-1:
           masses=[600,800,1000]
       
       
        if s.find('Higgs') !=-1:
           masses=[650,1000]
       
    
        for m in masses:
           
            mass_str=str(m)
            process_name=s+mass_str
            lumi_str=str("%.0f"%lumi)
            
            sys.stdout.write("\n\n--------------------------------\n\n")
            sys.stdout.write("STARTING\t\t")
            sys.stdout.write(process_name)
            sys.stdout.write("\n\n--------------------------------\n\n")

            if not os.path.isdir("output/Ntuple_%s"%(options.ntuple)):
                   os.system("mkdir output/Ntuple_%s"%(options.ntuple));

            if not os.path.isdir("output/Ntuple_%s/Lumi_%s"%(options.ntuple,lumi_str)):
                   os.system("mkdir output/Ntuple_%s/Lumi_%s"%(options.ntuple,lumi_str));


            if not os.path.isdir("output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s"%(options.ntuple,lumi_str,options.channel,process_name)):
                   os.system("mkdir output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s"%(options.ntuple,lumi_str,options.channel,process_name));
        
            if not os.path.isdir("output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log"%(options.ntuple,lumi_str,options.channel,process_name)):
                   os.system("mkdir output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log"%(options.ntuple,lumi_str,options.channel,process_name));
        
            #if not os.path.isdir("output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Data"%(options.ntuple,lumi_str,options.channel,process_name)):
                   #os.system("mkdir output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Data"%(options.ntuple,lumi_str,options.channel,process_name));
            
            LogFile1="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log/Log_Control_Plots_%s_%s.log"%(options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name)
            LogFile2="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log/Error_Log_Control_Plots_%s_%s.log"%(options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name)         
            output_Log=open(LogFile1,'w+')
            output_Log2=open(LogFile2,'w+')
            #output_Log="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log/Log_Control_Plots_%s_%s.log"%(options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name)
            p3 = subprocess.Popen(['python','MATTEO_make_Control_Plots.py','--lumi',lumi_str,'--mass',mass_str,'--channel',options.channel,'--sample',s,'--ntuple',options.ntuple],stdout=output_Log,stderr=subprocess.PIPE)
            
            for line in p3.stderr:
                sys.stdout.write(line)
                output_Log2.write(line)
            p3.wait()
            output_Log.close()
            output_Log2.close()
            
            sys.stdout.write("\n\n--------------------------------\n\n")
            sys.stdout.write("ENDED\t\t")
            sys.stdout.write(process_name)
            sys.stdout.write("\n\n--------------------------------\n\n")
            

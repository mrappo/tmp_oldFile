import os,commands
import sys
from optparse import OptionParser
import subprocess

parser = OptionParser()

#parser.add_option('--mass', action="store", type="string", dest="mass", default="600")
parser.add_option('--channel', action="store", type="string", dest="channel", default="mu")
#parser.add_option('--sample', action="store", type="string", dest="sample", default="BulkGraviton")
parser.add_option('--ntuple', action="store", type="string", dest="ntuple", default="WWTree_22sep_jecV7_lowmass")
(options, args) = parser.parse_args()

currentDir = os.getcwd();

lumi_true_value=2197.96

luminosity=[lumi_true_value]



sample=["BulkGraviton","Higgs"]

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
           
           if not os.path.isdir("output/Ntuple_%s"%(options.ntuple)):
                  os.system("mkdir output/Ntuple_%s"%(options.ntuple));

           if not os.path.isdir("output/Ntuple_%s/Lumi_%s"%(options.ntuple,lumi_str)):
                  os.system("mkdir output/Ntuple_%s/Lumi_%s"%(options.ntuple,lumi_str));


           if not os.path.isdir("output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s"%(options.ntuple,lumi_str,options.channel,process_name)):
                  os.system("mkdir output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s"%(options.ntuple,lumi_str,options.channel,process_name));
       
           if not os.path.isdir("output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/plots"%(options.ntuple,lumi_str,options.channel,process_name)):
                  os.system("mkdir output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/plots"%(options.ntuple,lumi_str,options.channel,process_name));

           if not os.path.isdir("output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Data"%(options.ntuple,lumi_str,options.channel,process_name)):
                  os.system("mkdir output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Data"%(options.ntuple,lumi_str,options.channel,process_name));

           if not os.path.isdir("output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Log"%(options.ntuple,lumi_str,options.channel,process_name)):
                  os.system("mkdir output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Log"%(options.ntuple,lumi_str,options.channel,process_name));
           
           
           os.system("python MATTEO_make_Cuts_Plots.py --mass %s --channel %s --sample %s --ntuple %s --lumi %s > output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Log/logfile_%s_%s.log"%(mass_str,options.channel,s,options.ntuple,lumi_str,options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name))
           
           
           
           log_output_dir="output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Log/"%(options.ntuple,lumi_str,options.channel,process_name)
           
           log_file=log_output_dir+"logfile_%s_%s.log"%(options.channel,process_name)
           output_Log=open(log_file,'a')
           
           #p3 = subprocess.Popen(['python','MATTEO_make_Cuts_Plots.py','--mass',mass_str,'--channel',options.channel,'--sample',s,'--ntuple',options.ntuple,'--lumi',lumi_str],stdout=output_Log,stderr=subprocess.PIPE)
           #p3.wait()
           
           #for line in p3.stderr:
           #     sys.stdout.write(line)
                #output_Log2.write(line)
           #p3.wait()
           
            #output_Log2.close()
            
           sys.stdout.write("\n\n--------------------------------\n\n")
           sys.stdout.write("ENDED\t\t")
           sys.stdout.write(process_name)
           sys.stdout.write("\n\n--------------------------------\n\n")
           
           output_Log.write("\n\n--------------------------------\n\n")
           output_Log.write("ENDED\t\t")
           output_Log.write(process_name)
           output_Log.write("\n\n--------------------------------\n\n")
           
           output_Log.close()
           #LogFile1="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log/Log_Control_Plots_%s_%s.log"%(options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name)
           # LogFile2="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log/Error_Log_Control_Plots_%s_%s.log"%(options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name)         
           # output_Log=open(LogFile1,'w+')
           # output_Log2=open(LogFile2,'w+')
            #output_Log="output/Ntuple_%s/Lumi_%s/Control_Plots_%s_%s/Log/Log_Control_Plots_%s_%s.log"%(options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name)
           # p3 = subprocess.Popen(['python','MATTEO_make_Control_Plots.py','--lumi',lumi_str,'--mass',mass_str,'--channel',options.channel,'--sample',s,'--ntuple',options.ntuple],stdout=output_Log,stderr=subprocess.PIPE)
            
           #for line in p3.stderr:
           #     sys.stdout.write(line)
           #     output_Log2.write(line)
           #p3.wait()
           # output_Log.close()
           # output_Log2.close()
            
           # sys.stdout.write("\n\n--------------------------------\n\n")
           # sys.stdout.write("ENDED\t\t")
           # sys.stdout.write(process_name)
           # sys.stdout.write("\n\n--------------------------------\n\n")
           
           
           
           #os.system("python MATTEO_make_Cuts_Plots.py --mass %s --channel %s --sample %s --ntuple %s --lumi %s > output/Ntuple_%s/Lumi_%s/BKG_Cuts_Evaluation_%s_%s/Log/logfile_%s_%s.log"%(mass_str,options.channel,s,options.ntuple,lumi_str,options.ntuple,lumi_str,options.channel,process_name,options.channel,process_name))
           
           

filename="MATTEO_Cut_result_%s.txt"%options.ntuple
output_dir="output/Ntuple_%s/Lumi_%s/"%(options.ntuple,lumi_str)
output_file_name1 =output_dir+"Cut_result.txt"
output_file_name2 =output_dir+"max_significance_file.txt"
p1 = subprocess.Popen(['cp','-r',filename,output_file_name1]) 
p1.wait()

p6 = subprocess.Popen(['cp','-r','MATTEO_max_significance_file.txt',output_file_name2]) 
p6.wait()

p2 = subprocess.Popen(['rm','-r',filename]) 
p2.wait()

p4 = subprocess.Popen(['rm','-r','MATTEO_bin_file.txt']) 
p4.wait()

p5 = subprocess.Popen(['rm','-r','MATTEO_max_significance_file.txt']) 
p5.wait()



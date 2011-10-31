#!/usr/bin/perl

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------


#PG lettura dei parametri da cfg file
#PG --------------------------------
print "reading ".$ARGV[0]."\n" ;

open (USERCONFIG,$ARGV[0]) ;

while (<USERCONFIG>)
  {
    chomp; 
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
#    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences{$var} = $value;
  }

$BASEDir                       = $User_Preferences{"BASEDir"} ; 
$INPUTFile = $User_Preferences{"INPUTFile"} ;
$OUTPUTFile = $User_Preferences{"OUTPUTFile"} ;

$EXEName =          $User_Preferences{"EXEName"};
$MONTHList= $User_Preferences{"MONTHList"};
@MONTHS = split(/,/, $MONTHList);
$EXECUTEName = $EXEName.".cpp " ;


print "BASEDir = "          .$BASEDir."\n" ;
print "INPUTFile = "    .$INPUTFile."\n" ;
print "OUTPUTFile = "   .$OUTPUTFile."\n" ;
print "EXEName= ".$EXEName."\n" ;
print "MONTHList= ".$MONTHList."\n" ;
print "EXECUTEName= ".$EXECUTEName."\n" ;

system("rm ".$EXEName) ;

$Compile="g++ -Wall -o ".$EXEName." `root-config --glibs --libs --cflags` ".$EXECUTEName ;

print "Complile= ".$Compile."\n" ;

system ($Compile) ; 

$jobNumber=0;

for($imonth=0; $imonth<@MONTHS; $imonth++)
{
  $jobNumber = $jobNumber+1;
}
  
print "NumberOfJobs = ".$jobNumber."\n";

#####################################################
# make jobs
#####################################################

 for($imonth=1; $imonth<=@MONTHS ; $imonth++)
    {
     $currDir = `pwd` ;
     
     $jobDir = $BASEDir."Jobs_Skim/".$EXEName.@MONTHS[$imonth-1] ;  
     system ("rm -r ".$jobDir." \n") ;
     system ("mkdir ".$jobDir." \n") ;
     

     ######################
     # make job files
     ######################    
    
     $tempBjob = $jobDir."/bjob_".@MONTHS[$imonth-1].".sh" ;
    
     $command = "touch ".$tempBjob ;             system ($command) ;
     $command = "chmod 777 ".$tempBjob ;         system ($command) ;
     $command = "cd ".$BASEDir ;                 system ("echo ".$command." >> ".$tempBjob) ;
     $command = "df -h";                         system ("echo ".$command." >> ".$tempBjob) ;

     $command = "./".$EXEName." ".@MONTHS[$imonth-1]." "."\\\"".$INPUTFile."\\\""." "."\\\"".$OUTPUTFile."\\\"";  system ("echo ".$command." >> ".$tempBjob); 
;
   ############
   # submit job."
   ############
   $command = "qsub -V -d ".$jobDir." -q longcms ".$tempBjob."\n" ;      
   print ($command."\n");
   system ($command);
    
   print "\n" ;
}  



   


    

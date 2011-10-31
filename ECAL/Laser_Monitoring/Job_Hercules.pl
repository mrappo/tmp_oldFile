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
$INPUTSAVEPath        = $User_Preferences{"INPUTSAVEPath"} ;
$OUTPUTSAVEPath    = $User_Preferences{"OUTPUTSAVEPath"} ;
$INPUTFile = $User_Preferences{"INPUTFile"} ;
$OUTPUTFile = $User_Preferences{"OUTPUTFile"} ;

$EXEName =          $User_Preferences{"EXEName"};
$FEDList = $User_Preferences{"FEDList"};
$MONTHList= $User_Preferences{"MONTHList"};
@FEDS = split(/,/, $FEDList);
@MONTHS = split(/,/, $MONTHList);
$EXECUTEName = $EXEName.".cpp " ;


print "BASEDir = "          .$BASEDir."\n" ;
print "INPUTSAVEPath = "    .$INPUTSAVEPath."\n" ;
print "OUTPUTSAVEPath = "   .$OUTPUTSAVEPath."\n" ;
print "INPUTFile = "    .$INPUTFile."\n" ;
print "OUTPUTFile = "   .$OUTPUTFile."\n" ;
print "EXEName= ".$EXEName."\n" ;
print "FEDList= ".$FEDList."\n" ;
print "MONTHList= ".$MONTHList."\n" ;
print "EXECUTEName= ".$EXECUTEName."\n" ;

system("rm ".$EXEName) ;

$Compile="g++ -Wall -o ".$EXEName." `root-config --glibs --libs --cflags` ".$EXECUTEName ;

print "Complile= ".$Compile."\n" ;

system ($Compile) ; 

$jobNumber=0;

for($imonth=0; $imonth<@MONTHS; $imonth++)
{
  
  for($ifed=0; $ifed<@FEDS ;  $ifed++ )
   {$jobNumber = $jobNumber+1;}
}
  
print "NumberOfJobs = ".$jobNumber."\n";

#####################################################
# make jobs
#####################################################

for($ifed = 1; $ifed <=  @FEDS; ++$ifed)
  { 
    for($imonth=0; $imonth<@MONTHS ; $imonth++)
    {
     $currDir = `pwd` ;
     
     $jobDir = $BASEDir."Jobs/".$EXEName.@FEDS[$ifed-1].@MONTHS[$imonth] ;  
     system ("rm -r ".$jobDir." \n") ;
     system ("mkdir ".$jobDir." \n") ;
     system ("cd ".$OUTPUTSAVEPath."\n");
    
     system ("rm -r ".$OUTPUTSAVEPath.$EXEName.@FEDS[$ifed-1]."\n");
    
     system ("mkdir ".$OUTPUTSAVEPath.$EXEName.@FEDS[$ifed-1]."\n");
     

    ######################
    # make job files
    ######################    
    
    $tempBjob = $jobDir."/bjob_".@FEDS[$ifed-1].@MONTHS[$imonth].".sh" ;
    
    $command = "touch ".$tempBjob ;             system ($command) ;
    $command = "chmod 777 ".$tempBjob ;         system ($command) ;
    $command = "cd ".$BASEDir ;                 system ("echo ".$command." >> ".$tempBjob) ;
    $command = "df -h";                         system ("echo ".$command." >> ".$tempBjob) ;

    $INPUT = sprintf($INPUTFile,@MONTHS[$imonth]) ;
    $OUTPUT = $OUTPUTSAVEPath.$EXEName.@FEDS[$ifed-1]."/".sprintf($OUTPUTFile,@FEDS[$ifed-1],@MONTHS[$imonth]) ;
    $DUMPEROutput = $OUTPUTSAVEPath.$EXEName.@FEDS[$ifed-1]."/Dumper_output_".@FEDS[$ifed-1].@MONTHS[$imonth]."_%d" ;
    $command = "./".$EXEName." ".@FEDS[$ifed-1]." "."\\\"".$INPUT."\\\""." "."\\\"".$OUTPUT."\\\""." "."\\\"".$DUMPEROutput.".txt"."\\\"";  system ("echo ".$command." >> ".$tempBjob); 

   ############
   # submit job."
   ############
   $command = "qsub -V -d ".$jobDir." -q longcms ".$tempBjob."\n" ;      
   print ($command."\n");
   system ($command);
    
   print "\n" ;
}  

}

   


    

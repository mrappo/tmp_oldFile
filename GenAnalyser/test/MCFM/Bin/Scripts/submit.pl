#!/usr/bin/perl
use POSIX ;

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------

#RG acquisition infos from CFG file
#RG --------------------------------
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

$HiggsList = $User_Preferences{"Higgs_Mass"};
@Higgs_mass = split(/,/, $HiggsList);
$run = $User_Preferences{"NRun"};
$FacScale = $User_Preferences{"Factorization_scale"};
@Factorization_scale = split(/,/, $FacScale);
$RinScale = $User_Preferences{"Rinormalization_scale"};
@Rinormalization_scale = split(/,/, $RinScale);
$CFG_Template = $User_Preferences{"CFG_Template"};
$baseName = $User_Preferences{"BaseName"};
$Mode = $User_Preferences{"DiagramPart"};
@DiagramMode = split(/,/, $Mode);
$BaseDir = $User_Preferences{"BaseDir"};
$icallReal = $User_Preferences{"Real_iteration"};
$icallVirt = $User_Preferences{"Virtual_iteration"};

# creation of .DAT and submission on LSF
$sampleJobListFile = "./lancia.sh";
open(SAMPLEJOBLISTFILE, ">", $sampleJobListFile);

for ($irun = 1; $irun <=$run ; $irun++) 
 {
  for($imass = 0; $imass <@Higgs_mass ; $imass++)
  {
    
    for($iFS = 0 ; $iFS< @Factorization_scale ; $iFS++)
    {
      
      for($iRS = 0 ; $iRS < @Rinormalization_scale ; $iRS++)
        {
       
         my $random_number = floor (rand () * 9999) ;
         $FSValue = @Factorization_scale[$iFS]*@Higgs_mass[$imass];
         $RSValue = @Rinormalization_scale[$iRS]*@Higgs_mass[$imass];

        
         if( @DiagramMode[0] == "real" && @DiagramMode[1] == "virt"  )
         {
          $fileName=$baseName."_".@DiagramMode[0]."_".@Higgs_mass[$imass]."_".$FSValue."_".$RSValue."_".$irun.".DAT" ;
           $command = "cat ".$BaseDir.$CFG_Template."   | sed -e s%SEEDFIXME%".$random_number.                                                                                    "%g | sed -e s%DiagramType%".@DiagramMode[0].
                                                   "%g | sed -e s%HIMASS%".@Higgs_mass[$imass].
                                                   "%g | sed -e s%FACTS%".$FSValue.
                                                   "%g | sed -e s%RINS%".$RSValue.
                                                   "%g | sed -e s%RUNNS%".$irun.
                                                   "%g | sed -e s%NrunCall%".$icallReal.
                                                   "%g | sed -e s%baseName%".$baseName.
                                                   "%g > ".$BaseDir.$fileName ;
         print $command."\n" ;
         system ($command) ;
        
          $fileName2=$baseName."_".@DiagramMode[1]."_".@Higgs_mass[$imass]."_".$FSValue."_".$RSValue."_".$irun.".DAT" ;
         
         $command2 =  "cat ".$BaseDir.$CFG_Template."   | sed -e s%SEEDFIXME%".$random_number.                                                                                      "%g | sed -e s%DiagramType%".@DiagramMode[1].
                                                   "%g | sed -e s%HIMASS%".@Higgs_mass[$imass].
                                                   "%g | sed -e s%FACTS%".$FSValue.
                                                   "%g | sed -e s%RINS%".$RSValue.
                                                   "%g | sed -e s%RUNNS%".$irun.
                                                   "%g | sed -e s%NrunCall%".$icallVirt.
                                                   "%g | sed -e s%baseName%".$baseName.
                                                   "%g > ".$BaseDir.$fileName2 ;
          print $command2."\n" ;
          system ($command2) ;

          $command = "bsub -cwd ./LSF -q 2nd \"cd \$PWD ; ./mcfm ".$BaseDir."/".$fileName."\"\n"; 
          print SAMPLEJOBLISTFILE  $command;
    
          print $command ;
#          system ($command) ;

          $command_2 = "bsub -cwd ./LSF -q 2nd \"cd \$PWD ; ./mcfm ".$BaseDir."/".$fileName2."\"\n" ; 
         
          print $command_2 ;
 #         system ($command_2) ;
	  print SAMPLEJOBLISTFILE  $command_2;

  
         }
       
      
         else{
         
         if(@DiagramMode[0] == "virt" && @DiagramMode[1] == "real" )
         {
          $fileName=$baseName."_".@DiagramMode[0]."_".@Higgs_mass[$imass].."_".$FSValue."_".$RSValue."_".$irun.".DAT" ;
         
          $command =  "cat ".$BaseDir.$CFG_Template."   | sed -e s%SEEDFIXME%".$random_number.                                                                                      "%g | sed -e s%DiagramType%".@DiagramMode[0].
                                                   "%g | sed -e s%HIMASS%".@Higgs_mass[$imass].
                                                   "%g | sed -e s%FACTS%".$FSValue.
                                                   "%g | sed -e s%RINS%".$RSValue.
                                                   "%g | sed -e s%RUNNS%".$irun.
                                                   "%g | sed -e s%NrunCall%".$icallVirt.
                                                   "%g | sed -e s%baseName%".$baseName.
                                                   "%g > ".$BaseDir.$fileName ;
          print $command."\n" ;
          system ($command) ;

          $fileName2=$baseName."_".@DiagramMode[1]."_".@Higgs_mass[$imass]."_".$FSValue."_".$RSValue."_".$irun.".DAT" ;

          $command2=  "cat ".$BaseDir.$CFG_Template."   | sed -e s%SEEDFIXME%".$random_number.                                                                                    "%g | sed -e s%DiagramType%".@DiagramMode[1].
                                                   "%g | sed -e s%HIMASS%".@Higgs_mass[$imass].
                                                   "%g | sed -e s%FACTS%".$FSValue.
                                                   "%g | sed -e s%RINS%".$RSValue.
                                                   "%g | sed -e s%RUNNS%".$irun.
                                                   "%g | sed -e s%NrunCall%".$icallReal.
                                                   "%g | sed -e s%baseName%".$baseName.
                                                   "%g > ".$BaseDir.$fileName2 ;

         print $command2."\n" ;
         system ($command2) ;
         
 
         system ($command);
  
	  $command = "bsub -cwd ./LSF  -q 2nd \"cd \$PWD ; ./mcfm ".$BaseDir."/".$fileName."\"\n";
          
         print $command ;
      #   system ($command) ;
         print SAMPLEJOBLISTFILE  $command;
  
	  $command_2 = "bsub -cwd ./LSF  -q 2nd \"cd \$PWD ; ./mcfm ".$BaseDir."/".$fileName_2."\"\n";
         
         print $command_2 ;
       #  system ($command_2) ;
	  print SAMPLEJOBLISTFILE  $command_2;


          }

          }
         
         
       } 
    }  
  }
 }

###### Job without full iteration set complete

####### Virtual
#grep -H "Intermediate"  output_ggH/HWW2jt_virt_mstw8nl_* | grep virt | grep .dat | tr "da" " " | tr "_" " " | awk '{print "bsub -cwd ./LSF -q 1nd \" cd $PWD ; ./mcfm REAL_TEST/"$7"_"$8"_"$9"_"$10"_"$11"_"$12"DAT\""}'
#in order to cancel output that are spourius derived from exited jobs


##### Real 
#grep -H "Intermediate"  output_ggH/HWW2jt_virt_mstw8nl_* | grep real | grep .dat | tr "da" " " | tr "_" " "  | awk '{if($20<9) print "bsub -cwd ./LSF -q 2nd \" cd $PWD ; ./mcfm REAL_TEST/"$7"_"$8"_"$9"_"$10"_"$11"_"$12".DAT\""}'

#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir+"/../";
ReducedTreeDir = "";

command = [
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_RSGraviton1000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_RSGraviton2000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_RSGraviton3000.cfg",
"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_RSGraviton4000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGraviton1000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGraviton2000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGraviton3000.cfg",
"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGraviton4000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_Wprime1000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_Wprime2000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_Wprime3000.cfg",
"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_Wprime4000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGWprime1000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGWprime2000.cfg",
#"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGWprime3000.cfg",
"./bin/VBFOptimizeSelections.exe cfg/VBFOptimizeSelections_InputCfgFile/Run2OptimizeSelection_BulkGWprime4000.cfg"
]

name = [
#"RSG1000",
#"RSG2000",
#"RSG3000",
"RSG4000",
#"BulkG1000",
#"BulkG2000",
#"BulkG3000",
"BulkG4000",
#"Wprime1000",
#"Wprime2000",
#"Wprime3000",
"Wprime4000",
#"BulkGWprime1000",
#"BulkGWprime2000",
#"BulkGWprime3000",
"BulkGWprime4000"
]

#MC
for a in range(len(command)):
        fn = "Job/Job_"+name[a];
        outScript = open(fn+".sh","w");
        print command[a];
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+command[a]);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
        command2 = "qsub -V -d "+currentDir+" -q shortcms "+currentDir+"/"+fn+".sh";
        os.system(command2);
        print command2

{
 ///==== cfg creator ====
 ///==== lunch with: root -l cfgCreator_MC_Local_Gen.cxx
 
 #include <iomanip>
 #include <string>
 #include "test/Read.cc"
 
 char outputDirectory[1000] = {"/home/raffaele/Programmi/RGerosa/AnalysisPackage_qqHWWlnulnu/output/out_NtupleProducer"};
 char    inputDirectory[1000] = {"/home/raffaele/Programmi/RGerosa/AnalysisPackage_qqHWWlnulnu/data/data_1011"};  
  
 char *nameSample[1000];
 char *nameHumanReadable[1000];
 char * xsection[1000];
 char * xsectionErrorUp[1000];
 char * xsectionErrorDown[1000];
 
 char nameFileIn[1000] = {"test/Latinos/sample_GenLevel.txt"};
 
 int numberOfSamples = ReadFileXSection(nameFileIn,nameSample,nameHumanReadable,xsection, 
		                         xsectionErrorUp,xsectionErrorDown);
 
 TString toDoShell;
 
 toDoShell= Form ("mkdir test/Latinos/dir_cfg_skimmed_MC_GenLevel");
 system (toDoShell.Data());
  
 toDoShell = Form ("mkdir %s", outputDirectory);
 system (toDoShell.Data());
 
 for (int iSample = 0; iSample < numberOfSamples; iSample++){
   std::string name_samples = nameHumanReadable[iSample];
  if (name_samples  == "DATA") continue;
  
  std::ofstream myfile;
  char nameFile[1000];
  sprintf(nameFile,"test/Latinos/dir_cfg_skimmed_MC_GenLevel/NtupleProducerNT_Gen_%s.cfg",nameSample[iSample]);
  myfile.open (nameFile);
  
  myfile << std::fixed;
  myfile << std::setprecision(20);
  myfile << "[Input]" << std::endl;
  myfile << "treeName = SimpleNtupleGen/SimpleNtupleGen " << std::endl;

  myfile << "inputFile = " << inputDirectory << "/" << nameSample[iSample] << ".root" << std::endl;
  
  char toDo[1000];
  sprintf(toDo,"double xsec = %s",xsection[iSample]);
  gROOT->ProcessLine(toDo);
  myfile << "histoNameEvents = events" << std::endl;
  myfile << std::endl;
  myfile << "entryMIN = 0 " << std::endl;
  myfile << "entryMAX = -1 " << std::endl;
  myfile << "entryMOD = 1000 " << std::endl;
  myfile << "nStepToDo = 6 " << std::endl;
  myfile << std::endl;
  
  myfile << "[Options]" << std::endl;
  myfile << "XSection = "<<xsection[iSample]<< std::endl;
  myfile << "XSectionErrorUp = "<<xsectionErrorUp[iSample]<< std::endl;
  myfile << "XSectionErrorDown = "<<xsectionErrorDown[iSample]<< std::endl;
  myfile << std::endl;
  
  myfile <<  "[Selection]" << std::endl;
  myfile << std::endl;
  myfile << "[Output]" << std::endl;
    
  myfile << "outFileName = " << outputDirectory << "/out_NtupleProducer_" << nameSample[iSample] << ".root" << std::endl;
  myfile << std::endl;
  myfile.close();
 
 }
 
}



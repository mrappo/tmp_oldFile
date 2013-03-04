#include "ReadInputFile.h"

int ReadInputSampleFile(const std::string & InputSampleList, std::vector<std::string> & NameSample, std::vector<std::string> & NameReducedSample, std::vector<int> & ColorSample, 
                        std::vector<double>& SampleCrossSection, std::vector<int> & NumEntriesBefore){


  std::ifstream inputFile (InputSampleList.c_str());
  std::string buffer;

  if(inputFile.fail()) return -1; 

  while(!inputFile.eof()){
  
    getline(inputFile,buffer);

    if(buffer.empty() || !buffer.find("#") ) continue ;
    std::stringstream line(buffer);

    std::string  NameSampleTemp;
    std::string  NameReducedSampleTemp;
    std::string  ColorSampleTemp;
    std::string  SampleCrossSectionTemp;
    std::string  NumEntresBeforeTemp;

    line >> NameSampleTemp >> NameReducedSampleTemp >> ColorSampleTemp >> SampleCrossSectionTemp >> NumEntresBeforeTemp ;

    NameSample.push_back(NameSampleTemp);   NameReducedSample.push_back(NameReducedSampleTemp);
    ColorSample.push_back(atoi(ColorSampleTemp.c_str())); SampleCrossSection.push_back(atof(SampleCrossSectionTemp.c_str())) ; NumEntriesBefore.push_back(atoi(NumEntresBeforeTemp.c_str()));

  }

  return NameSample.size() ;

}


int ReadInputVariableFile( const std::string & InputVariableList , std::vector<std::string> & Variables, std::vector<int> & VariablesNbin,
                           std::vector<double> & VariablesMinValue, std::vector<double> & VariablesMaxValue, std::vector<std::string> & VariablesTitle ){


  std::ifstream inputFile (InputVariableList.c_str());
  std::string buffer;

  if(inputFile.fail()) return -1; 

  while(!inputFile.eof()){
  
    getline(inputFile,buffer);

    if(buffer.empty() || !buffer.find("#") ) continue ;
    std::stringstream line(buffer);

    std::string  VariablesTemp;
    std::string  VariablesNbinTemp;
    std::string  VariablesMinValueTemp;
    std::string  VariablesMaxValueTemp;
    std::string  VariablesTitleTemp;

    line >> VariablesTemp >> VariablesNbinTemp >> VariablesMinValueTemp >> VariablesMaxValueTemp >> VariablesTitleTemp ;

    Variables.push_back(VariablesTemp); VariablesNbin.push_back(atoi(VariablesNbinTemp.c_str())) ; VariablesMinValue.push_back(atof(VariablesMinValueTemp.c_str())); 
    VariablesMaxValue.push_back(atof(VariablesMaxValueTemp.c_str())); VariablesTitle.push_back(VariablesTitleTemp);

  }

  return Variables.size() ;

}


int ReadInputCutFile( const std::string & InputCutList , std::vector<std::string> & CutList){


  std::ifstream inputFile (InputCutList.c_str());
  std::string buffer;

  if(inputFile.fail()) return -1; 

  while(!inputFile.eof()){
  
    getline(inputFile,buffer);

    if(buffer.empty() || !buffer.find("#") ) continue ;

    CutList.push_back(buffer);
  }

  return CutList.size() ;

}

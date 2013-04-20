#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>
#include <map>
#include <algorithm>


#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TLatex.h"



void DrawStackError(THStack* hs, const std::string & Labels, const TH1F* dataHist, const std::map<int,double> & SystematicErrorMap, const bool & isLog = false, const bool & isLabel = false, const double & syst = 0.044);

void DrawStackError(THStack* hs, const std::string & Labels,  const std::map<int,double> & SystematicErrorMap,  const bool & isLog = false, const bool & isLabel = false, const double & syst = 0.044);

void SetTotalSystematicVector( std::vector<double> & SysError, THStack* hs, const std::map<int,double> & SystematicErrorMap, const double & syst = 0.044);

void LatexCMS (const double & lumi, const std::string & LeptonType, const bool & isLabel = false);

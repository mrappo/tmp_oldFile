#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>

#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TLatex.h"



void DrawStackError(THStack* hs, const double & syst, const std::string & Labels, const TH1F* dataHist);

void DrawStackError(THStack* hs, const double & syst, const std::string & Labels);

void LatexCMS (double lumi);

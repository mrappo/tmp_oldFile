// Select the 2011 Ntuples as a function of the month
// compile with:
// g++ `root-config --libs --cflags` SkimNtuple_month.cpp -o SkimNtuple_month
// usage:
// ./SkimNtuple_month number        (e.g: ./SkimNtuple_month 3   for March 2011)        

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TDirectory.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TTimeStamp.h"


#include "./LaserDataAnalysis.h"

using namespace std;

int main(int argc, char ** argv)
{
  // select fed
  int month = atoi(argv[1]);
  string month_name [12] = {"Gen","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

  if (month < 0 || month > 12) {
   cout << "Give as argument a number between 1 and 12 corresponding to the chosen month" << endl;
   return 0;
  }

  if (month > 10) {
   cout << "No data corresponding to the chosen month..yet!" << endl;
   return 0;
  }


  //Get old tree 
  TChain *tx = new TChain("x");
  tx->Add(argv[2]);
 
  init_ttree(tx, &x);

  TTree *oldtree = (TTree*)tx;

  Long64_t nentries = oldtree->GetEntries();
  cout<< "Number of entries in the tree : " << nentries << endl;


  //Create a new file + a clone of old tree in new file
  char fname[1000];
  sprintf(fname,argv[3],(month_name[month - 1]).c_str());
  TFile *newfile = new TFile(fname,"recreate");
  TTree *newtree = oldtree->CloneTree(0);

  int lowerRun, upperRun = 0;

  if (month == 1) {
    lowerRun = 153941;
    upperRun = 156054;
  }

  if (month == 2) {
    lowerRun = 156225;
    upperRun = 159130;
  }

  if (month == 3) {
    lowerRun = 159248;
    upperRun = 161732;
  }

  if (month == 4) {
    lowerRun = 161846;
    upperRun = 163726;
  }

  if (month == 5) {
    lowerRun = 163762;
    upperRun = 166139;
  }

  if (month == 6) {
    lowerRun = 166243;
    upperRun = 168148;
  }

  if (month == 7) {
    lowerRun = 168266;
    upperRun = 172266;
  }

  if (month == 8) {
    lowerRun = 172319;
    upperRun = 174912;
  }

  if (month == 9) {
    lowerRun = 175118;
    upperRun = 177519;
  }

  if (month == 10) {
    lowerRun = 177624;
    upperRun = 179998;
  }


  for (int i=0; i<nentries; i++) {

    if (i%10000000==0) cout << "Analyzing entry " << i << endl;
    oldtree->GetEntry(i);

    if ( x.run > lowerRun && x.run < upperRun ) newtree->Fill();
  }
  
  newtree->Print();
  newtree->AutoSave();
  
  delete newfile;


}

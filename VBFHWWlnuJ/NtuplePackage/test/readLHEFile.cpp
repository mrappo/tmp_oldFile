#include "LHEReader.h"

#include <cmath>
#include <cstdlib>






int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there                                                                                                                               
  if(argc != 2)
  {
    std::cerr << ">>>>> readLHEFile.cpp::usage: " << argv[0] << " LHEFileName" << std::endl ;
    return 1;
  }
  
  int entryMAX = -1;
  int entryMODULO = 10000;
  LHEReader myReader(argv[1]);
  
  
  
  
  
  
  //*********************
  // LOOP OVER THE EVENTS
  std::vector<int> n_plus(5);
  std::vector<int> n_minus(5);
  
  
  
  int entry = 0;
  while( myReader.GetNextEvent() == true )
  {
    if( entry == entryMAX ) break;
    if( (entry%entryMODULO) == 0 ) std::cout << ">>>>> readLHEFile::GetEntry " << entry << std::endl;
    //myReader.PrintEvent();
    
    
    int nup = myReader.GetNUP();
    std::vector<int> idup = myReader.GetIDUP();
    std::vector<int> istup = myReader.GetISTUP();
    
    
    bool isPlus = false;
    bool isMinus = false;
    int nJets = 0;
    for(int i = 0; i < nup; ++i)
    {
      if( istup.at(i) != 1 ) continue;
      
      if( (idup.at(i) == 11) || (idup.at(i) == 13) || (idup.at(i) == 15) )
        isMinus = true;
      if( (idup.at(i) == -11) || (idup.at(i) == -13) || (idup.at(i) == -15) )
        isPlus = true;
      if( (abs(idup.at(i)) == 1) ||
          (abs(idup.at(i)) == 2) ||
          (abs(idup.at(i)) == 3) ||
          (abs(idup.at(i)) == 4) ||
          (abs(idup.at(i)) == 5) ||
          (abs(idup.at(i)) == 6) ||
          (abs(idup.at(i)) == 21) )
        ++nJets;
    }
    
    if( isPlus == true )  n_plus.at(nJets) += 1;
    if( isMinus == true ) n_minus.at(nJets) += 1;
    
    
    ++entry;
  }
  
  
  
  for(int i = 0; i < 5; ++i)
  {
    std::cout << "nJets = " << i << "\tn+ + n- = " << std::fixed << std::setw(6) << n_plus.at(i) + n_minus.at(i)
                                 << "\tn+ = "    << std::fixed << std::setw(5) << n_plus.at(i)
                                 << "\tn- = "    << std::fixed << std::setw(5) << n_minus.at(i)
                                 << "\tn+/n- = " << std::fixed << std::setprecision(4) << std::setw(6) << 1.*n_plus.at(i)/n_minus.at(i)
                                 << " +/- " << std::fixed << std::setprecision(3) << std::setw(6) << 1.*n_plus.at(i)/n_minus.at(i) * sqrt( 1./n_plus.at(i) + 1./n_minus.at(i) )
                                 << std::endl;
  }
  
  return 0;
}

// Versione 0.0.1

#ifndef StepRecognition_h
#define StepRecognition_h

#include "TProfile.h"
#include "TTimeStamp.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


// ----- Define Thresholds -----

#define THR_AMP       0.11
#define THR_dR        0.10
#define THR_TIMESTEP  30
#define THR_TIMEFLUC  10

// -----------------------------


// -----------------------------------------------------------------------------------

//Return a true value if the two numbers on input have the same sign
bool isSameSign (double a, double b)
{
  if ( (a * b) > 0 )
    return true;
  
  else
    return false;
 
}


// -----------------------------------------------------------------------------------

//Evaluate the separation between two points as their (signed) difference
double EvalDelta (double y2, double y1)
{
  return y2 - y1;
}


// -----------------------------------------------------------------------------------

//Evaluate the (signed) increment between a point (x1,y1) and a point (x2,y2) as dR := (y2 - y1)/(x2 -x1)
double EvalIncrement (double y2, double y1, double x2, double x1)
{
  return ( (y2 - y1)/(x2 -x1) );
}


// -----------------------------------------------------------------------------------

//Return a true value if there's a possible candidate step between two points of the TProfile passed as input
bool isStep (TProfile *profile, int i, int j,std::ofstream & outFile)
{
  
  //No step by default
  bool step = false;
  
  //The step condition could be given both by the fact that the relative difference is above threshold (THR_AMP) 
  //or by the fact that the increment is above THR_dR.
  //The second condition is especially useful when the two points are close to each other (< 5 bins)
  
  if ( (fabs(EvalDelta(profile -> GetBinContent(i),profile -> GetBinContent(j))) > THR_AMP) ||
       (fabs(j - i) < 5 && fabs(EvalIncrement(profile -> GetBinContent(j),profile -> GetBinContent(i),j,i)) > THR_dR))
   {
      step = true;
      
      //outFile << "Applying a Delta/Amp test to points " << i << " and " << j << std::endl;
   }
  
  return step;
}


// -----------------------------------------------------------------------------------

//Takes as input an IndexVector which, at first, contains only two points of the TProfile. Then cycles on the following ones 
//and checks if the conditions for a step candidate are preserved. When the points get stable or the end of the TProfile is
//reached, return a true value
bool isStepLong (std::vector<int> &IndexVector, TProfile *profile, std::ofstream & outFile )
{
    
  //Increment between the first two points of the vector
  float firstdR = EvalIncrement (profile->GetBinContent(IndexVector.at(IndexVector.size()-1)),
			         profile->GetBinContent(IndexVector.at(IndexVector.size()-2)),
			         IndexVector.at(IndexVector.size()-1), IndexVector.at(IndexVector.size()-2)); 
  
  //Vector containing the following increments
  std::vector<float> dR;
  
  
  //The analyzed points "i" go from the last one of the vector given as input to (ideally) the end of the TProfile 
  for (int i = IndexVector.at(IndexVector.size()-1); i < profile -> GetNbinsX(); i++)
  {
    //Skip the empty bins
    if (profile->GetBinContent(i+1) == 0)   continue;
    
    //Add to the dR vector the increment between a point "i+1" and the previous one (i.e. the last of IndexVector)
    dR.push_back( EvalIncrement (profile->GetBinContent(i+1),
				 profile->GetBinContent(IndexVector.at(IndexVector.size()-1)), 
				 i+1, IndexVector.at(IndexVector.size()-1)) );
   
    
    
    //If there is an increment above threshold and in the same direction, add the point to IndexVector
    if ( fabs(dR.at(dR.size()-1)) > THR_dR && isSameSign(firstdR, dR.at(dR.size()-1)) )
    {
       IndexVector.push_back(i+1);
       //outFile << "Cand. Step rising at point " << i+1 << std::endl;
    }  
    
    
    //If not, analyze the two possible options:
    else
    {
      
      // 1) Big excursion but with a sign opposite to the previous increment: check if it's a fluctuation
      if ( fabs(dR.at(dR.size()-1)) > THR_dR && !isSameSign(firstdR, dR.at(dR.size()-1)) )
      {
        
        IndexVector.push_back(i+1);
        //outFile << "Cand. Step decreasing at point " << i+1 << std::endl;

        //Cycle on the points following the decrease 
        int k = i+1;

        while (k < profile -> GetNbinsX())
        {

          //If empty: increment "k" and continue with the for cycle
          if (profile->GetBinContent(k+1) == 0) {
             k++;
             continue;
          }

          //Evaluate the increment between "k" and the last point added to IndexVector
          float tmpdR = EvalIncrement (profile->GetBinContent(k+1), 
                                       profile->GetBinContent(IndexVector.at(IndexVector.size()-1)),
				       k+1, IndexVector.at(IndexVector.size()-1));
          
          
          //If the points "k" keep on decreasing, add them to the IndexVector
          if ( !isSameSign(tmpdR,firstdR) )
          {
            IndexVector.push_back(k);
            k++;

            //outFile<< "The point " << k << " has the opposite sign" << std::endl;
          } 


          //As soon as there is a new inversion, break the while on "k" and check the
          //options A or B below
          else break; //exit from the while cycle  
  
        }


        IndexVector.push_back(k);


        // A) At least three points added to the decrease, but they reach an amplitude compatible
        // with the starting one: is a fake step 
        if (fabs(k - i) > 3 && fabs(EvalDelta(profile->GetBinContent(k),profile->GetBinContent(IndexVector.at(0)))) < THR_AMP)  
          {
            //outFile << "Fake step between points " << i << " and " << k << std::endl;
            return false;
          }

        // B) Even if some points are decreasing, the overall structure is compatible with a step
        else 
          { 
            //outFile << "Cand.Step getting stable between " << i+1 << " and " << k << " exiting.." << std::endl;
            return true;          
          }
          
      }
      
      
      // 2) The points are under threshold (getting stable)
      else 
      {
        IndexVector.push_back(i+1);
        //outFile << "Cand.Step getting stable at point " << i+1 << " exiting.." << std::endl;

        return true;
      }
      
    }
    
  } //end of for cycle  
  
 return true;
 
 //N.B: The IndexVector now contains all the points from the first one to the last one of the step candidate

}


// -----------------------------------------------------------------------------------

//Count the flat plateau of the step which, at first, contain only one (starting) point. The points are added to the StepSize vector.
//Return false if the step is finished before the end ot the TProfile, true if not.
bool CountStepExtension (std::vector<int> & StepSize, TProfile *profile, std::ofstream & outFile)
{
  
  bool isFinish = true;
  
  for (int i = StepSize.at(0); i < profile -> GetNbinsX(); i++)
  {
    
    //Skip the empty bins
    if (profile->GetBinContent(i+1) == 0)
    {
       //outFile << "Zero Bin found" << std::endl;
       continue;
    }
     
    
    //If the points remain stable on the plateau, add them to the StepSize vector and continue with the for cycle
    if ( (fabs(EvalDelta(profile -> GetBinContent(i+1),profile -> GetBinContent(StepSize.at(StepSize.size()-1)))) < THR_AMP) ||
         (fabs(i+1 - StepSize.at(StepSize.size()-1)) < 5 && 
          fabs(EvalIncrement(profile -> GetBinContent(i+1),profile->GetBinContent(StepSize.at(StepSize.size()-1)),
                             i+1,StepSize.at(StepSize.size()-1))) < THR_dR) )   
      
    {
       //outFile << "Add point " << i+1 << " to the Cand. Step" << std::endl;
       StepSize.push_back(i+1);      
       continue; 
    } 
    
    
    //The plateau condition is not satisfied anymore, and the TProfile is not finished yet: return false
    if ((i+1) <= profile -> GetNbinsX())
    {
       //outFile << "End of the Cand.Step" << std::endl;
       isFinish = false;

       return isFinish;
    }
      
  } //end of for cycle  

  //The plateau condition is not satisfied anymore, and the TProfile is finished: return true
  isFinish = true;
  
  return isFinish;
  
}


// -----------------------------------------------------------------------------------

//Basic tool for the step identification in the MATACQ amplitude TProfiles
void Step_ID(TProfile *amplitude, TString & output_dumper_steps)
{
   
   std::vector<int> IndexVector;
   std::vector<int> StepSize;
   
   TTimeStamp date_step_rise(2011, 2, 22, 0, kFALSE, 0);
   TTimeStamp date_step_inf (2011, 2, 22, 0, kFALSE, 0);
   TTimeStamp date_step_sup (2011, 2, 22, 0, kFALSE, 0);
   
   std::ofstream outFile(output_dumper_steps.Data(),std::ios::out);   //Save the output on a text file
   
   for(int it = 0; it < amplitude -> GetNbinsX(); it++)
   {  
     
     //Skip the empty bins
     if(amplitude -> GetBinContent(it) == 0)
     { 
       //outFile << "Empty Bin" << std::endl;
       continue;
     }
     
     
     //Look for the first non-empty bin following "it"
     int j = it + 1;
     
     while (amplitude->GetBinContent(j) == 0 && j <= amplitude -> GetNbinsX())
     { 
       j++;
     }

     //You have reached the end of the TProfile: exit
     if(j >= amplitude -> GetNbinsX())
     {
       outFile << "End Sample!" << std::endl;
       break;
     }

     //No candidate step between "it" and "j": move forward!
     if (isStep(amplitude,j,it,outFile) == false) 
     {
       //outFile << "Non Step Candidate between " << it << " and " << j << std::endl;
       
       //Move "it" to the point BEFORE "j" and let it be incremented by the for cycle
       it = j-1;
       continue;
     } 

     //There is a candidate step between "it" and "j": cycle on the bins after "j" and look at their behaviour
     else 
     {
       
       //outFile << "Step Candidate at point " << it << std::endl;
       IndexVector.push_back(it);
       IndexVector.push_back(j);
      
       
       //If there is no step or it finds something not compatible with a step (like a fluctuation)
       if ( !isStepLong(IndexVector,amplitude,outFile) ) {
          
           //outFile << "Step Size Long method has failed to recognize a step structure" << std::endl;
           
           //Move the iterator to the last point of the IndexVector
           it = IndexVector.at(IndexVector.size()-1);

           IndexVector.clear();

           continue;
         }
      
       //If it finds a step-like structure: counts the extension of that step
       else {
	  
           StepSize.push_back(IndexVector.at(IndexVector.size()-1));

           //Call the CountStepExtension methof
           bool resultCS = CountStepExtension(StepSize,amplitude,outFile);
           
           //outFile << "Good Step Size Long at point " << IndexVector.at(IndexVector.size()-1) << std::endl;
          
          
           //End of the step found AND it is extended for a reasonable amount of time
           if ( resultCS == false && EvalDelta(StepSize.at(StepSize.size()-1),StepSize.at(0)) > THR_TIMESTEP) 
           {
              float step_rise = amplitude -> GetBinCenter(IndexVector.at(0));
              float step_inf =  amplitude -> GetBinCenter(StepSize.at(0));
              float step_sup =  amplitude -> GetBinCenter(StepSize.at(StepSize.size()-1));
             
              date_step_rise.SetSec(int (step_rise));
              date_step_inf.SetSec(int (step_inf));
              date_step_sup.SetSec(int (step_sup));
             
              outFile << " ---------------------------------------" << std::endl;
              outFile << "           STEP IDENTIFICATION          " << std::endl;
              outFile << " ---------------------------------------" << std::endl;
              outFile << "riseBegin: " << date_step_rise.AsString( ) << std::endl;
	      outFile << "stepBegin: " << date_step_inf.AsString( )  << std::endl;
	      outFile << "stepEnd:   " << date_step_sup.AsString( )  << std::endl << std::endl;


              //Move the iterator to the penultimate point of the StepSize vector
              it = StepSize.at(StepSize.size()-1) - 1;
              
              //Clear the vectors
              IndexVector.clear();
              StepSize.clear();
              
              continue;

           }

          //End of the step found AND it is extended for a time compatible with a small fluctuation
          if ( resultCS == false && EvalDelta(StepSize.at(StepSize.size()-1),StepSize.at(0)) > THR_TIMEFLUC &&
	                            EvalDelta(StepSize.at(StepSize.size()-1),StepSize.at(0)) < THR_TIMESTEP ) 
          {
              float step_rise = amplitude -> GetBinCenter(IndexVector.at(0));
              float step_inf =  amplitude -> GetBinCenter(StepSize.at(0));
              float step_sup =  amplitude -> GetBinCenter(StepSize.at(StepSize.size()-1));
             
              date_step_rise.SetSec(int (step_rise));
              date_step_inf.SetSec(int (step_inf));
              date_step_sup.SetSec(int (step_sup));

              outFile << " ---------------------------------------" << std::endl;
              outFile << "       FLUCTUATION IDENTIFICATION       " << std::endl;
              outFile << " ---------------------------------------" << std::endl;
              outFile << "riseBegin: " << date_step_rise.AsString( ) << std::endl;
	      outFile << "flucBegin: " << date_step_inf.AsString( )  << std::endl;
	      outFile << "flucEnd:   " << date_step_sup.AsString( )  << std::endl << std::endl;


             //Move the iterator to the penultimate point of the StepSize vector
             it = StepSize.at(StepSize.size()-1) - 1;
             
             //Clear the vectors
             IndexVector.clear();
             StepSize.clear();
             
             continue;

          }
          
          
          //End of the step found BUT is NOT extended for a reasonable amount of time
          if ( resultCS == false && EvalDelta(StepSize.at(StepSize.size()-1),StepSize.at(0)) < THR_TIMEFLUC) 
          {  

              float step_rise = amplitude -> GetBinCenter(IndexVector.at(0));
              float step_inf =  amplitude -> GetBinCenter(StepSize.at(0));
              float step_sup =  amplitude -> GetBinCenter(StepSize.at(StepSize.size()-1));
             
              date_step_rise.SetSec(int (step_rise));
              date_step_inf.SetSec(int (step_inf));
              date_step_sup.SetSec(int (step_sup));

              outFile << " ---------------------------------------" << std::endl;
              outFile << "            FAKE FLUCTUATION            " << std::endl;
              outFile << "----------------------------------------" << std::endl;
              
              outFile << "riseBegin: " << date_step_rise.AsString( ) << std::endl;
	      outFile << "flucBegin: " << date_step_inf.AsString( )  << std::endl;
	      outFile << "flucEnd:   " << date_step_sup.AsString( )  << std::endl << std::endl;
 
 
              //Move the iterator to the penultimate point of the StepSize vector
              it = StepSize.at(StepSize.size()-1) - 1;
              
              //Clear the vectors
              IndexVector.clear();
              StepSize.clear();
              
              continue;
          }


          //You have arrived at the end of the TProfile
          if (resultCS == true)
          {  

              float step_rise = amplitude -> GetBinCenter(IndexVector.at(0));
              float step_inf =  amplitude -> GetBinCenter(StepSize.at(0));
              float step_sup =  amplitude -> GetBinCenter(StepSize.at(StepSize.size()-1));
             
              date_step_rise.SetSec(int (step_rise));
              date_step_inf.SetSec(int (step_inf));
              date_step_sup.SetSec(int (step_sup));
            
              outFile << " ---------------------------------------" << std::endl;
              outFile << "          END OF MONTH SAMPLE           " << std::endl;
              outFile << "----------------------------------------" << std::endl;
             
	      outFile << "riseBegin: " << date_step_rise.AsString( ) << std::endl;
	      outFile << "stepBegin: " << date_step_inf.AsString( )  << std::endl;
	      outFile << "stepEnd:   " << "End  Current Month"      << std::endl << std::endl;
            
            
              //Clear the vectors
              IndexVector.clear();
              StepSize.clear();
              
              outFile << "End Sample!" << std::endl;
              
              break;
          }
       }
      
     } 
      
   } //end of for cycle on "it"
    
  return;  
      
}      
      
#endif

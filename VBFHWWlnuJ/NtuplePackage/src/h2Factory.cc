#include "h2Factory.h"
#include <iostream>
#include <string>

h2Factory::h2Factory (std::string fileName, 
                    bool print) :
 m_fileName (fileName) ,
 m_print (print)                     
{}


//PG --------------------------------------------------------   


h2Factory::~h2Factory () 
{
//  if (m_H2content.size () && m_print)
//    for (std::map <std::string, h2Chain *>::iterator mapIt = m_H2content.begin () ;
//         mapIt != m_H2content.end () ;
//         ++mapIt)
//      {
//        mapIt->second->Print () ; 
//        //PG FIXME generalizzare le opzioni di stampa    
//      }
  if (m_H2content.size () && m_fileName != "no")
    {
      TFile histosFile (m_fileName.c_str (),"update") ;
      histosFile.cd () ;
      for (std::map <std::string, h2Chain *>::iterator mapIt = m_H2content.begin () ;
           mapIt != m_H2content.end () ;
           ++mapIt)
        {
          mapIt->second->Write (mapIt->first, histosFile) ; 
        }
      histosFile.Close () ;
    }  
}


//PG --------------------------------------------------------   


void
h2Factory::add_h2 (std::string baseName, std::string baseTitle, 
                   int nbinsx, double minx, double maxx,
                   int nbinsy, double miny, double maxy, int NUM)
{
  if (m_H2content.find (baseName) != m_H2content.end ())
    {
      std::cerr << "ERROR : histograms series " << baseName
                << " already existing, NOT replaced" << std::endl ;
      return ;                
    }
  h2Chain * dummy = new h2Chain (baseName, baseTitle, 
                                 nbinsx, minx, maxx,
                                 nbinsy, miny, maxy, NUM) ;
  m_H2content[baseName] = dummy ;     
  return ;                          
}
                 

//PG --------------------------------------------------------   


void 
h2Factory::Fill (const std::string & name, int i, double valx, double valy) 
  {
    if (m_H2content.find (name) != m_H2content.end ())
      m_H2content[name]->Fill (i,valx, valy) ;
    return ;
  }


//PG --------------------------------------------------------   


h2Chain * 
h2Factory::operator[] (const std::string& name)
  {
    if (m_H2content.find (name) == m_H2content.end ())
      {
        std::cerr << "ERROR : histograms series " << name
                  << " not existing, returning the first one" << std::endl ;
        return m_H2content.begin ()->second ;                
      }
    return m_H2content[name] ;
  
  }


//PG --------------------------------------------------------   


void 
h2Factory::Print (int isLog, int rebin)
{
  if (!m_H2content.size ()) return ;
  for (std::map <std::string, h2Chain *>::iterator mapIt = m_H2content.begin () ;
       mapIt != m_H2content.end () ;
       ++mapIt)
    {
      mapIt->second->Print (isLog, rebin) ; 
    }
}


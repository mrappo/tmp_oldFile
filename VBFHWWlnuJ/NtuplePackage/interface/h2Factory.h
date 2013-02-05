#ifndef h2Factory_h
#define h2Factory_h

#include <map>
#include <string>
#include <functional>
#include "hChain.h"
#include "h2Chain.h"

class h2Factory
{
  public :
    h2Factory (std::string fileName = "no", bool print = true) ;
    ~h2Factory () ;

    void add_h2 (std::string baseName, std::string baseTitle, 
                 int nbinsx, double minx, double maxx,
                 int nbinsy, double miny, double maxy, int NUM) ;

    template <class UnaryFunction>
    void applyToAll (UnaryFunction function) 
      {
        for (std::map <std::string, h2Chain *>::iterator mapIt = m_H2content.begin () ;
             mapIt != m_H2content.end () ;
             ++mapIt)
          {
            function (mapIt->second) ;
          }
      }

    void Fill (const std::string & name, int i, double valx, double valy) ;
    h2Chain * operator[] (const std::string& name) ;
    void Print (int isLog = 0, int rebin = 1) ;

  private :
  
    std::string m_fileName ;
    bool m_print ;
    std::map <std::string, h2Chain *> m_H2content ;

} ;

#endif

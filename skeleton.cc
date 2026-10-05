//////////////////////////////////////////////////////////////////////////////////////////////
// skeleton.cc
//
// This is the main skeleton program to use StocasticSolver class.
//
// There are two ways to run this code in Linux command line:
// 1. using root environment
//    $ root skeleton.cc
//
// 2. using standalone c++ (compile and run respectively)
//    $ g++ skeleton.cc -o skeleton -O3 `root-config --cflags --glibs`
//    $ ./skeleton
//
// Outputs: Energy eigenvalues {Etar} and {Ecal} and discrete radial potential values.
//
// Author: Ahmet.Bingul(at)cern.ch
// Date  ; Oct 2026
//////////////////////////////////////////////////////////////////////////////////////////////

#include "TApplication.h"
#include "StocasticSolver.h"

/* skeleton function */
void skeleton()
{
  StocasticSolver sc;   // generate an object

  sc.SetSeed(1981);     // determine seed
  sc.ShowPlot(true);    // draw true potential profile and reco potential points
  sc.LoopMC(5e6);       // start 5M Monte Carlo iterations
}

/* Actual main program for standalone c++ application */
int main(int argc, char** argv){
  // start root application
  TApplication *myApp = new TApplication("theApp",&argc, argv);
  // call the function above
  skeleton();
  // run root application
  myApp->Run();
  return 0.0;
}

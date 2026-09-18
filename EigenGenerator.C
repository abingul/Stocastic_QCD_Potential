en//////////////////////////////////////////////////////////////////////////////////////////////
// EigenGenerator.C
//
// This ROOT macro solves Schrodinger equation to obtain energy eigenvalues and wavefunctions
// for a simple non-relativistic q-qbar bound system. We define an effective potential as
//
//     Veff(r) = -alpha_s/r +  (sigma/mu)*(1-exp(-mu*r))
//
// here first term is called the coulomb part and second term is the confinement part
// which is in this case a screeing function. We assume that quantum numbers are L = S = 0.
//
// How to run in Linux command line:
//
// 1. using root environment
//    $ root EigenGenerator.C
//
// 2. using standalone c++ (compile and run respectively)
//    $ g++ EigenGenerator.C -o EigenGenerator -O3 `root-config --cflags --glibs`
//    $ ./EigenGenerator
//
// For the given configuration output are the graph of Vconf(r) and values
//
//    Eigenvalues of 12 states:
//    0.624983, 1.39276, 1.92537, ... ,4.37172,
//
// Author: Ahmet.Bingul@cern.ch
// Date  ; Sep 2026
//////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include "TCanvas.h"
#include "TGraph.h"
using namespace std;

#define eps        1.0e-12
#define tinyValue  0.01
#define K          1.0 // m/hbar^2 value
#define toler      1.0e-6

// Screening potential as a confinement part
double Vconf(double r){
  double sigma = 0.4, mu = 0.08;
  return sigma / mu * (1-exp(-mu*r));
}

// Effective potential.
// This part can be extended by adding angular momentum and spin dependent part.
double Veff(double r){
   double alpha_s = 0.4;
   return -alpha_s /(r + eps) + Vconf(r); // + other terms
}

// to be used in numerical solution
double F(double r, double E){
   return 2 * K * (E-Veff(r));
}

// the genrator function
void EigenGenerator(){
  const int  N = 2000;
  double rmin  = -0;
  double rmax  = +20;
  double h     = (rmax-rmin)/N;
  double E     = 0.0;
  double c     = h*h/12.0;
  int    n     = 12;

  double r[N+1] = {0.0}, psi[N+1] = {0.0}, con[N+1], pot[N+11];
  double psiold, dE;

  for(int i=0;i<N+1;i++){
      r[i] = rmin + i*h;
      pot[i] = Veff(r[i]);
      con[i] = Vconf(r[i]);
  }

  cout << "Eigenvalues of " << n << " states:" << endl;
  for (int j = 0; j<n; j++)
  {
    dE = 0.1;

    for(int i=0; i<N+1; i++) psi[i] = 0.0;
    psi[1] = tinyValue;

    // initial numerov steps
    for(int i=1; i<N; i++){
         psi[i+1] = (  2*(1-5*c*F(r[i],E))*psi[i] - (1+c*F(r[i-1],E))*psi[i-1] )
                  / ( 1 + c*F(r[i+1],E) );
    }
    psiold = psi[N];

    while(1){
      for(int i=0; i<N+1; i++) psi[i] = 0.0;
      psi[1] = tinyValue;
      E = E + dE;

      for(int i=1; i<N; i++){
         psi[i+1] = (  2*(1-5*c*F(r[i],E))*psi[i] - (1+c*F(r[i-1],E))*psi[i-1] )
                  / ( 1 + c*F(r[i+1],E) );
      }

      if (psi[N] * psiold < 0){
         E = E - dE;
         dE = dE / 2.0;

         if (dE < toler){
            cout << E <<  ", ";
            E = E + 0.1;
            break;
         }
      }
    } // end of while
  } // endof for j

  // draw confinement function
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TGraph *gr1 = new TGraph(N+1, r, con);
  gr1->Draw("APL");

  cout << "\nALL OK." << endl;
} // end of numeroc

// main function
int main(){
  EigenGenerator();
}

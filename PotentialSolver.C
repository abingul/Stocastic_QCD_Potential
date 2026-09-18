//////////////////////////////////////////////////////////////////////////////////////////////
// PotentialSolver.C
//
// For the given eigenvalues, this ROOT macro solves Schrodinger equation to obtain
// approximate form of the confinement part of a simple non-relativistic q-qbar bound system.
// We define an effective potential as
//
//     Veff(r) = -alpha_s/r +  Vconf(r)
//
// here first term is called the coulomb part and second term is the confinement part
// We assume that quantum numbers are L = S = 0.
//
// Generalized description can be found at:
// arXiv:2605.23896
//
// It is suggested that first run EigenGenerator.C before running this code
//
// How to run in Linux command line:
//
// 1. using root environment
//    $ root PotentialSolver.C
//
// 2. using standalone c++ (compile and run respectively)
//    $ g++ PotentialSolver.C -o PotentialSolver -O3 `root-config --cflags --glibs`
//    $ ./PotentialSolver
//
// Output utput are eigenvalues and potential values
//
// Author: Ahmet.Bingul@cern.ch
// Date  ; Sep 2026
//////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
using namespace std;


class StocasticSolver{
   public:

     TRandom3 *rnd;
     TCanvas  *c1;
     double   *Vr, *Vr_best, *rho, *c, *rval, *psi;
     //double   *Etar, *Ecal;
     double   eps, tinyValue, K, toler, rmin, rmax;
     double   dr, drho, cons, chi2max;
     bool     showPlot;
     int      N, NN, n, seed;
     vector<double> Etar, Ecal;

     //---------------------------------------------------------------------------------------------
     StocasticSolver(){
        /* *** options ******************* */
        n         = 12;                     // number of states
        eps       = 1.0e-10;                // epsilon
        tinyValue = 0.01;                   // to be used in numerical solution
        toler     = 1.0e-3;                 // toleance of the numerical solution
        N         = 2000;                   // number of sample points taken from wave function
        NN        = 9;                      // number of stocastic slopes
        rmin      = 0.0;                    // fm
        rmax      = 20.0;                   // fm
        showPlot  = false;                  // to show progress of the stocastic function
        chi2max   = 1e-3,                   // terminational condition of the program

        /* *** scalar or array objects *****/
        dr        = (rmax-rmin)/double(N);  // radial step size
        cons      = dr*dr/12.0;             // constant,, to be used in numerov method
        K         = 1.0;                    // m/hbar^2, value, to be used in numberov method
        c         = new double[NN];         // stocastic slopes
        Vr        = new double[NN+1];       // random potential values
        Vr_best   = new double[NN+1];       // best random potential values
        rho       = new double[NN+1];       // radial stocastic step points, fm
        drho      = (rmax-rmin)/double(NN); // radial stocastic step size
        rval      = new double[N+1];        // radial values
        psi       = new double[N+1];        // wavefunction

        if(showPlot) c1 = new TCanvas("c1","c1",800,600);

        rho[0] = 0.0;
        for(int i=0; i<NN; i++)
           rho[i+1] = rho[i] + drho;

        rval[0] = 0.0;
        for(int i=0; i<N; i++)
           rval[i+1] = rval[i] + dr;

        Etar.resize(n);
        Ecal.resize(n);
        // You should provide eigenvalues.
        // e.g. the following values are taken from the output of generator.C
        Etar = {0.624983, 1.39276, 1.92537, 2.34877, 2.70275, 3.00647, 3.27107, 3.50385, 3.71141, 3.90917, 4.12433, 4.37172};
     }
     //---------------------------------------------------------------------------------------------
     ~StocasticSolver(){
        if(showPlot){
          c1->Print(Form("output_%d.pdf",seed));
          c1->Print(Form("output_%d.png",seed));
        }
     }
     //---------------------------------------------------------------------------------------------
     // Cornell (linear) potental
     double Vcornell(double r){
        double alpha_s = 0.4, b = 0.8;
        return -alpha_s /(r + eps) + b*r;
     }
     //---------------------------------------------------------------------------------------------
     // Screeening potential (to test)
     double Vscreening(double r){
        double sigma = 0.4, mu = 0.08;
        return sigma / mu * (1-exp(-mu*r));
     }
     //---------------------------------------------------------------------------------------------
     // Random or stocastic potential
     double Vrand(double r){
        // c[i] values are set in LoopMC function
        Vr[0] = 0;
        for(int i=0;i<NN;i++)
           Vr[i+1] = Vr[i] + c[i] * drho;

        // caompute potential value at distance r using linear interpolation
        for(int i=0; i<NN; i++){
           if(r > rho[i] && r <= rho[i+1])
              return  Vr[i] + c[i] * (r - rho[i]);
        }
        return 0.0;
     }
     //---------------------------------------------------------------------------------------------
     // Effective potential = coulomb + confinement
     double Veff(double r){
        double alpha_s = 0.4;
        return -alpha_s /(r + eps) + Vrand(r); // + other terms if needed
     }
     //---------------------------------------------------------------------------------------------
     // This function is used in numerical solution
     double F(double r, double E){
        return 2 * K * (E - Veff(r));
     }
     //---------------------------------------------------------------------------------------------
     double Solve(){
       double psiold, dE, chi2 = 0.0;

       Ecal[0] = 0;
       for (int j = 0; j<n; j++)
       {
          dE = 0.1;
          for(int i=0; i<N+1; i++) psi[i] = 0.0;
          psi[1] = tinyValue;

         // initial numerov steps to get final value of wfn
         for(int i=1; i<N; i++){
             psi[i+1] = ( 2*(1-5*cons*F(rval[i],Ecal[j]))*psi[i] - (1+cons*F(rval[i-1],Ecal[j]))*psi[i-1] )
                      / ( 1 + cons*F(rval[i+1],Ecal[j]) );
         }
         psiold = psi[N];

         // recalculate wfn until it converges
         while(1){
             for(int i=0; i<N+1; i++) psi[i] = 0.0;
             psi[1] = tinyValue;
             Ecal[j] = Ecal[j] + dE;

             for(int i=1; i<N; i++){
                psi[i+1] = ( 2*(1-5*cons*F(rval[i],Ecal[j]))*psi[i] - (1+cons*F(rval[i-1],Ecal[j]))*psi[i-1] )
                         / ( 1 + cons*F(rval[i+1],Ecal[j]) );
             }

            if (psi[N] * psiold < 0){
                Ecal[j] = Ecal[j] - dE;
                dE = dE / 2.0;
                if (dE < toler) break;
             }
         } // end of while

         chi2 += (Ecal[j] - Etar[j]) * (Ecal[j] - Etar[j]);
         Ecal[j+1] = Ecal[j] + 0.1;

       } // end of for j

       return chi2;
     }
     //---------------------------------------------------------------------------------------------
     void ShowPlot(bool sp = false){
       showPlot = sp;
       if(showPlot) c1 = new TCanvas("c1","c1",800,600);
     }
     //---------------------------------------------------------------------------------------------
     // you should call tgis method at the beginning
     void SetSeed(int s = 42){
       rnd = new TRandom3(seed);
       seed = s;
     }
     //---------------------------------------------------------------------------------------------
     int GetSeed(){
       return seed;
     }
     //---------------------------------------------------------------------------------------------
     // Monte Carlo trials to form sample points of the appoximate continement poterntial function
     void LoopMC(int nTrials = 10000){
        double chi2best = 1e9, chi2;

        cout << "Starting Monte Carlo trials for nTrials =  " << nTrials << " and seed = " << GetSeed() << endl;

        for(int trial=0; trial < nTrials; trial++)
        {
          // set random slopes here
          for(int i = 0;i<NN; i++)
             c[i] = rnd->Uniform(0.0, 0.5);

          // compute chi2
          chi2 = Solve();

          // reasonable check
          if(chi2 < chi2best){
             chi2best = chi2;
             cout.precision(4);
             cout << fixed;

             cout <<"Trail# = " << trial << "   chi2 = " << chi2 << " seed = " << GetSeed() << endl;

             cout << "True  Energies  : ";
             for(int j=0; j<n; j++) cout << Etar[j] << " ";
             cout << endl;

             cout << "Calc. Energies  : ";
             for(int j=0; j<n; j++) cout << Ecal[j] << " ";
             cout << endl << endl;

             cout << "True  Potential : ";
             for(int j=0; j<NN+1; j++) {
               double r = rho[j];
               cout << Vscreening(r) << " ";
             }
             cout << endl;

             cout << "Calc. Potential : ";
             for(int j=0; j<NN+1; j++) {
               Vr_best[j] = Vr[j];
               cout << Vr[j] << " ";
             }
             cout << endl << endl << endl;

             if(showPlot) DrawPotential();
          }

          // terminational condition
          if(chi2 < chi2max) {
              cout << "Stop loop since chi2 < chi2max.\n";
              break;
          }

        } // trial
     }
     //---------------------------------------------------------------------------------------------
     void DrawPotential(){
        double *Vscr = new double[N+1];
        double *Vrnd = new double[NN+1];
        for(int i=0;i<N+1; i++) Vscr[i] = Vscreening(rval[i]);
        for(int i=0;i<NN+1;i++) Vrnd[i] = Vr_best[i];

        TGraph *gr1 = new TGraph(N+1,  rval, Vscr);
        TGraph *gr2 = new TGraph(NN+1, rho,  Vrnd);

        gr1->SetTitle(Form("Random potential for seed %d",seed));
        gr1->GetXaxis()->SetTitle("Radial distance, r [fm]");
        gr1->GetYaxis()->SetTitle("Confinement Potential, V_{conf}(r) [GeV]");

        gr1->SetLineWidth(2);
        gr1->SetMaximum(7);
        gr1->Draw("APL");

        gr2->SetMarkerStyle(20);
        gr2->SetMarkerSize(2);
        gr2->Draw("P");
        c1->Update();

        delete [] Vscr;
        delete [] Vrnd;
     }

}; // end of class

/*
    main functions
*/

void PotentialSolver(){

  StocasticSolver sc;
  sc.SetSeed(time(NULL)%1093);
  sc.ShowPlot(true);
  sc.LoopMC(5e6);

  cout << "ALL OK." << endl;
}

int main(){
   PotentialSolver();
}

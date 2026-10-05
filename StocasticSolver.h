//////////////////////////////////////////////////////////////////////////////////////////////
// StochasticSolver.h
//
// Supporting C++/ROOT class for the benchmark validation presented in Appendix A of:
// A. Bingül and A. Özpineci, "A Stochastic Approach for Determining the Quark Confinement
// Potential of Charmonia," arXiv:2605.23896 [hep-ph].
//
// Core Functionality:
// 1. Forward Solver: Computes target (exact) reference eigenvalues {Etar} for a
//    spin-independent screened interaction potential (L = S = 0) by numerically solving
//    the radial Schrödinger equation using the Numerov integration scheme.
//
// 2. Inverse Reconstruction: Reconstructs the spatial form of the confinement/screening
//    potential and reconstructed energies (Ecal) via stochastic Monte Carlo optimization,
//    fitting the candidate potential to reproduce the mock eigenvalues generated in Step 1.
//
// Interaction Potential Model:
//    V_eff(r) = - (4 * alpha_s) / (3 * r) + (sigma / mu) * (1 - exp(-mu * r))
//    where the first term represents one-gluon exchange (Coulomb) and the second term
//    models color screening / string breaking.
//
//
// Author : Ahmet.Bingul(at)cern.ch
// Date   : September 2026
//////////////////////////////////////////////////////////////////////////////////////////////

#ifndef STOCASTIC_SOLVER_H
#define STOCASTIC_SOLVER_H

// C++ headers
#include <iostream>
#include <cmath>

// ROOT headers
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
     double   eps, tinyValue, K, toler, rmin, rmax;
     double   dr, drho, cons, chi2max;
     double   alpha_s, sigma, mu;
     bool     showPlot, sortSlopes;
     int      N, NN, n, seed;
     vector<double> Etar, Ecal;

     //---------------------------------------------------------------------------------------------
     // constuctor method
     StocasticSolver(){
        /* *** options ******************* */
        n         = 12;                     // number of states to evaluate
        eps       = 1.0e-10;                // epsilon, a small number to eliminate division error, 1/(r+eps)
        tinyValue = 1.0e-6;                 // to be used in numerical solution
        toler     = 1.0e-4;                 // toleance of the numerical solution
        N         = 500;                    // number of sample points to be used in wavefunction integration
        NN        = 10;                     // number of stocastic slopes
        rmin      = 0.0;                    // fm, minumum radial distance
        rmax      = 20.0;                   // fm, maximum radial distance
        showPlot  = false;                  // to show progress of the stocastic function
        chi2max   = 5.0e-3;                 // to terminate program when chi2 < chi2max
        alpha_s   = 0.4;                    // QCD running coupling constant
        sigma     = 0.6;                    // string tension in screening potential
        mu        = 0.1;                    // screening scale
        sortSlopes= false;                  // true if you want to sort the slopes to get plateau at large distances

        /* *** scalar or array objects *****/
        dr        = (rmax-rmin)/double(N);  // radial step size for wfn calculation
        cons      = dr*dr/12.0;             // constant, to be used in Numerov method
        K         = 1.0;                    // m/hbar^2, value, to be used in Numerov method
        c         = new double[NN];         // stocastic slopes
        Vr        = new double[NN+1];       // random potential values
        Vr_best   = new double[NN+1];       // best random potential values
        rho       = new double[NN+1];       // radial stocastic step points, fm
        drho      = (rmax-rmin)/double(NN); // radial stocastic step size
        rval      = new double[N+1];        // radial values
        psi       = new double[N+1];        // wavefunction

        /* *** initilization ***************/
        if(showPlot) c1 = new TCanvas("c1","c1",800,600);

        rho[0] = 0.0;
        for(int i=0; i<NN; i++)
           rho[i+1] = rho[i] + drho;

        rval[0] = 0.0;
        for(int i=0; i<N; i++)
           rval[i+1] = rval[i] + dr;

        // size of both vectors is n
        Etar.resize(n+1); // target     eigenvalues
        Ecal.resize(n+1); // calculated eigenvalues

        cout << "StocasticSolver has been initilized." << endl;
        cout.precision(4);
        cout << fixed;

        // first compute target eigenvalues
        // and expectation value of <r>
        CalculateTargetEigenValues();
     }
     //---------------------------------------------------------------------------------------------
     // destructer method
     ~StocasticSolver(){
        cout << "StocasticSolver finalized." << endl;
        if(showPlot){
          c1->Print(Form("output_%d.pdf",seed));
          c1->Print(Form("output_%d.png",seed));
        }
     }
     //---------------------------------------------------------------------------------------------
     // The linear potential if required
     double Vlinear(double r){
        double b = 0.3;
        return b*r;
     }
     //---------------------------------------------------------------------------------------------
     // The Screeening potential
     double Vscreening(double r){
        return sigma / mu * ( 1-exp(-mu*r) );
     }
     //---------------------------------------------------------------------------------------------
     // The confinement part of the potential
     double Vconf(double r){
        //return Vlinear(r);
        return Vscreening(r);
     }
     //---------------------------------------------------------------------------------------------
     // Random or stocastic potential
     double Vrand(double r){
        // c[i] values are set in LoopMC function
        Vr[0] = 0;
        for(int i=0; i<NN; i++)
           Vr[i+1] = Vr[i] + c[i] * drho;

        // compute potential value at distance r using linear interpolation
        for(int i=0; i<NN; i++){
           if(r > rho[i] && r <= rho[i+1])
              return  Vr[i] + c[i] * (r - rho[i]);
        }
        return 0.0;
     }
     //---------------------------------------------------------------------------------------------
     // Reconstruced (effective) potential = coulomb + stocastic confinement
     double Veff(double r){
        return -4.0*alpha_s /(3.0*(r + eps)) + Vrand(r); // + other terms if needed
     }
     //---------------------------------------------------------------------------------------------
     // true (target) potential = coulomb + screening
     double Vtar(double r){
        return -4.0*alpha_s /(3.0*(r + eps)) + Vconf(r); // + other terms if needed
     }
     //---------------------------------------------------------------------------------------------
     // This function is used in numerical solution using stocastic potential
     double F(double r, double E){
        return 2 * K * (E - Veff(r));
     }
     //---------------------------------------------------------------------------------------------
     // This function is used in numerical solution using true potential
     double Ksqr(double r, double E){
        return 2 * K * (E - Vtar(r));
     }
     //---------------------------------------------------------------------------------------------
     // Calculates the global target (true) states, Etar[j] where j = 0,1,...,n using Numerov Method
     void CalculateTargetEigenValues(){
        double dE;
        int diverge = 0, last_diverge = 0, max_iter = 1000;

        Etar[0] = 0.0;
        for(int j=0; j<n; j++) // loop over n states
        {
           dE = 0.15;
           for(int it = 0;it < max_iter; it++){
             // numerov integration in backward direction
             psi[N] = 0;
             psi[N-1] = tinyValue;
             for (int i = N-1; i > 0; --i) {
                 double f_prev = 2.0 * (1.0 - 5.0*cons *  Ksqr(rval[i],Etar[j])) * psi[i]
                               - (1.0 + cons * Ksqr(rval[i+1],Etar[j])) * psi[i+1];
                 f_prev /= (1.0 + cons * Ksqr(rval[i-1],Etar[j]));
                 psi[i-1] = f_prev;
             }
             // check divergence
             if(psi[0] > 0.0)  diverge = +1;
             else              diverge = -1;
             if( diverge * last_diverge < 0 ) dE = -dE / 2.0;
             last_diverge = diverge;
             // update energy
             Etar[j] += dE;
             if (fabs(dE) < toler) break;
           }

           // go to next state
           Etar[j+1] = Etar[j] + 0.05;
        }

        cout << "Target eigenvalues are computed up to " << n  << " states "
             << "for sigma = " << sigma << " and mu = "  << mu << " as follows:"
             << endl;

        cout << "Energies   : ";
        for(int j=0; j<n; j++)
          cout << Etar[j] <<" ";
        cout << endl;
     }
     //---------------------------------------------------------------------------------------------
     // Normalize the current wavefunction via trapezoidal integration method
     void NormalizeWavefunction(){
        double sum = 0.0;
        for(int i=1; i<N-1; i++)
           sum += psi[i]*psi[i];

        double norm = sqrt(sum * dr);
        for(int i=0; i<N+1; i++)
           psi[i] /= norm;
     }
     //---------------------------------------------------------------------------------------------
     // Compute <r> using current wfn via trapezoidal integration method
     double GetMean_r(){
        NormalizeWavefunction();
        double sum = 0.0;
        for(int i=1; i<N-1; i++)
           sum += rval[i]*psi[i]*psi[i];

        return sum * dr;
     }
     //---------------------------------------------------------------------------------------------
     // Returns the current radial wavefunction if required
     vector<double> GetRadialWavefunction(){
        NormalizeWavefunction();
        vector<double> Rwfn(N+1);
        for(int i=0; i<N+1; i++)
           Rwfn[i] = psi[i]/(rval[i]+eps);
        return Rwfn;
     }
     //---------------------------------------------------------------------------------------------
     // Evaluates the global predicted (calculated) states, Ecal[j] where j = 0,1,...,n
     // using Numerov Method and returns chi2 value based on energy differences
     double GetChi2(){
       double psiold, dE, chi2 = 0.0;
       int diverge = 0, last_diverge = 0, max_iter = 1000;

       Ecal[0] = 0;
       for (int j = 0; j<n; j++)
       {
          dE = 0.2;
          for(int it = 0;it < max_iter; it++){
             // numerov integration in backward direction
             psi[N] = 0;
             psi[N-1] = tinyValue;
             for (int i = N-1; i > 0; --i) {
                 double f_prev = 2.0 * (1.0 - 5.0*cons *  F(rval[i],Ecal[j])) * psi[i]
                               - (1.0 + cons * F(rval[i+1],Ecal[j])) * psi[i+1];
                 f_prev /= (1.0 + cons * F(rval[i-1],Ecal[j]));
                 psi[i-1] = f_prev;
             }
             // check divergence
             if(psi[0] > 0.0)  diverge = +1;
             else              diverge = -1;
             if( diverge * last_diverge < 0 ) dE = -dE / 2.0;
             last_diverge = diverge;
             // update energy
             Ecal[j] += dE;
             if (fabs(dE) < toler) break;
          }
          // compute chi-square using calculated and target energies
          chi2 += (Ecal[j] - Etar[j]) * (Ecal[j] - Etar[j]);

          // initlize next energy state
          Ecal[j+1] = Ecal[j] + 0.1;

       } // end of for j

       return chi2;
     }
     //---------------------------------------------------------------------------------------------
     // to show uptated plots and save true potential profile and reconstructed data points
     void ShowPlot(bool sp = false){
         showPlot = sp;
         if(showPlot) c1 = new TCanvas("c1","c1",800,600);
     }
     //---------------------------------------------------------------------------------------------
     // Intilize the seed of the random number generator.
     // You should call this method at the beginning of the program
     void SetSeed(int s = 42){
         rnd = new TRandom3(seed);
         seed = s;
     }
     //---------------------------------------------------------------------------------------------
     // Get seed value (default is 42)
     int GetSeed(){
         return seed;
     }
     //---------------------------------------------------------------------------------------------
     // This function constructs a loop of Monte Carlo trials to form random (stochastic)
     // sample points of the appoximate confinement potential function
     void LoopMC(int nTrials = 10000){
        double chi2best = 1e9, chi2;

        // check if random number generator is initilized
        if(rnd == NULL) rnd = new TRandom3(42);

        cout << "Starting Monte Carlo for =  " << nTrials << " trials and seed = " << GetSeed() << endl;

        for(int trial=0; trial < nTrials; trial++)
        {
          // First set random slopes here
          for(int i = 0; i<NN; i++)
             c[i] = rnd->Uniform(0.0, 0.5);

          // Then sort slopes in decreasing order if necessary
          // This approach is useful for zero-slope plateau at large distances
          if(sortSlopes)
             sort(c, c + NN, std::greater<double>());

          // Then compute chi2
          chi2 = GetChi2();

          // Finally update the values if it is better and print out performance values
          if(chi2 < chi2best){
             chi2best = chi2;
             cout << "Trail# = " << trial << " chi2 = " << chi2 << " seed = " << GetSeed() << endl;
             cout << "Slopes          : ";
             for(int j=0; j<NN; j++) cout << c[j] << " ";
             cout << endl;

             cout << "True  Energies  : ";
             for(int j=0; j<n; j++) cout << Etar[j] << " ";
             cout << endl;

             cout << "Calc. Energies  : ";
             for(int j=0; j<n; j++) cout << Ecal[j] << " ";
             cout << endl;

             cout << "True  Potential : ";
             for(int j=0; j<NN+1; j++) cout << Vconf(rho[j]) << " ";
             cout << endl;

             cout << "Calc. Potential : ";
             for(int j=0; j<NN+1; j++) {
               Vr_best[j] = Vr[j];
               cout << Vr[j] << " ";
             }
             cout << endl << endl << endl;

             // you can watch the updated potential
             if(showPlot) DrawPotential();
          }

          // terminational condition
          if(chi2 < chi2max) {
              cout << "End of Monte Carlo trials since chi2 < chi2max.\n";
              return;
          }

        } // MC trial

        cout << "End of Monte Carlo trials. Note that chi2 > chi2max." << endl;
     }
     //---------------------------------------------------------------------------------------------
     // Plot the target confinement function and predicted sample points
     void DrawPotential(){
        double *Vcon = new double[N+1];
        double *Vrnd = new double[NN+1];
        for(int i=0;i<N+1; i++) Vcon[i] = Vconf(rval[i]);
        for(int i=0;i<NN+1;i++) Vrnd[i] = Vr_best[i];

        TGraph *gr1 = new TGraph(N+1,  rval, Vcon);
        TGraph *gr2 = new TGraph(NN+1, rho,  Vrnd);

        gr1->SetTitle(Form("Random potential for seed %d",GetSeed()));
        gr2->SetTitle(Form("Random potential for seed %d",GetSeed()));
        gr1->GetXaxis()->SetTitle("Radial distance, r [fm]");
        gr1->GetYaxis()->SetTitle("Confinement Potential [GeV]");

        gr1->SetLineWidth(2);
        gr2->SetMarkerStyle(20);
        gr2->SetMarkerSize(2);
        gr2->Draw("AP");
        gr1->Draw("PL");
        c1->Update();

        delete [] Vcon;
        delete [] Vrnd;
     }
     //---------------------------------------------------------------------------------------------

}; // end of class

#endif

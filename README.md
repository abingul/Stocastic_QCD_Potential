Here We provide two independent ROOT macro programs.

* EigenGenerator.C

  solves non relaativistic Schrodinger equation to obtain energy eigenvalues and wavefunctions
  for a simple q-qbar bound system where the effective potential is defined as
  
     Veff(r) = -alpha_s/r +  (sigma/mu)*(1-exp(-mu*r))

 Here first term is called the Coulomb part and second term is the confinement part
 which is in this case a screeing function. We assume that quantum numbers are L = S = 0.
 For the given configuration outputa this program are the graph of confinement function and 
 correspondig eigenvalues like:

    Eigenvalues of 12 states:
    0.624983, 1.39276, 1.92537, ... ,4.37172,

* PotentialSolver.C

 For the given eigenvalues, this ROOT macro solves Schrodinger equation to obtain
 approximate form of the confinement part of a simple non-relativistic q-qbar bound system.
 We define an effective potential as
 
     Veff(r) = -alpha_s/r +  Vconf(r)
     
 The aim is to extract Vconf(r) using a stocastic (Monte Carlo) approach if eigenvalues are given.

The details of the study can be found at: arXiv:2605.23896

Here We provide two independent ROOT macro programs.

* EigenGenerator.C

solves Schrodinger equation to obtain energy eigenvalues and wavefunctions
for a simple q-qbar bound system where the effective potential is defined as
```math
      V_{eff}(r) = -\frac{\alpha_s} {r} +  \frac{\sigma} {\mu} \left( 1-\exp(-\mu r) \right)
```
Here first term is called the Coulomb part and second term is the confinement part
which is in this case a screeing function. We assume that quantum numbers are $L = S = 0$.
For the given configuration outputs of this program are the graph of the confinement function and 
correspondig eigenvalues like:

    Eigenvalues of 12 states:
    0.624983, 1.39276, 1.92537, ... ,4.37172

* PotentialSolver.C

For the given eigenvalues, this ROOT macro solves Schrodinger equation to obtain
approximate form of the confinement part of a simple q-qbar bound system.
We define an effective potential as
```math
     V_{eff}(r) = -\frac{\alpha_s} {r} +  V_{conf}(r)
```
The aim is to extract $V_{conf}(r)$ using a stocastic (Monte Carlo) approach if eigenvalues are given.
Output utput are eigenvalues and potential values like

   . . .
   Trail# = 27798   chi2 = 0.0144 seed = 669
   True  Energies  : 0.6250 1.3928 1.9254 2.3488 2.7027 3.0065 3.2711 3.5038 3.7114 3.9092 4.1243 4.3717
   Calc. Energies  : 0.6219 1.3484 1.9156 2.3656 2.6688 3.0422 3.3313 3.4922 3.6938 3.9297 4.1688 4.4281

   True  Potential : 0.0000 0.8144 1.4961 2.0668 2.5445 2.9444 3.2792 3.5595 3.7941 3.9905
   Calc. Potential : 0.0000 0.8387 1.3887 2.2306 2.3856 3.1690 3.3479 3.3742 3.8925 4.6044
   . . .

The details of the study can be found at: arXiv:2605.23896

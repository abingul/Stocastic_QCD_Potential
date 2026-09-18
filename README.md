**The details of the study can be found at: arXiv:2605.23896**

Here We provide two independent ROOT macro programs.

* **EigenGenerator.C**

solves Schrodinger equation to obtain energy eigenvalues and wavefunctions
for a simple q-qbar bound system where the effective potential is defined as
```math
      V(r) = -\frac{\alpha_s} {r} +  \frac{\sigma} {\mu} \left( 1-\exp(-\mu r) \right)
```
Here first term is called the Coulomb part and second term is the confinement part
which is in this case a screeing function. We assume that quantum numbers are $L = S = 0$.
For the given configuration outputs of this program are the graph of the confinement function and 
correspondig eigenvalues like:

    Eigenvalues of 12 states:
    0.624983, 1.39276, 1.92537, 2.34877, 2.70275, 3.00647, 3.27107, 3.50385, 3.71141, 3.90917, 4.12433, 4.37172

These values can be input for PotentialSolver() function given below.

* **PotentialSolver.C**

For the given eigenvalues, this ROOT macro solves Schrodinger equation to obtain
approximate form of the confinement part of a simple q-qbar bound system.
We define an effective potential as
```math
     V(r) = -\frac{\alpha_s} {r} +  V_{conf}(r)
```
The aim is to extract $V_{conf}(r)$ using a stocastic (Monte Carlo) approach if eigenvalues are given.
Outputs of this program are true and predicted eigenvalues, true and predicted potential values and the graph of the sample points like

    . . .
    Trail# = 2882999   chi2 = 0.0042 seed = 669
    True  Energies  : 0.6250 1.3928 1.9254 2.3488 2.7027 3.0065 3.2711 3.5038 3.7114 3.9092 4.1243 4.3717
    Calc. Energies  : 0.6469 1.3969 1.9547 2.3891 2.6859 3.0031 3.2859 3.5031 3.6875 3.9156 4.1328 4.3750

    True  Potential : 0.0000 0.8144 1.4961 2.0668 2.5445 2.9444 3.2792 3.5595 3.7941 3.9905
    Calc. Potential : 0.0000 0.8754 1.4646 2.2177 2.4531 2.9272 3.3763 3.4309 3.8719 3.9300
    . . .

An example graph is shown below for seed = 669. The target confinement function (curve) and computed sampled 
stochastic potential (points) fit fairly well after about 2.9 million iterations.


<img width="796" height="572" alt="output_669" src="https://github.com/user-attachments/assets/67fb9d06-5928-4211-ae7d-3426a8268c32" />

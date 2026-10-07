# Stochastic QCD Potential Solver

[![arXiv](https://img.shields.io/badge/arXiv-2605.23896-b31b1b.svg)](https://arxiv.org/abs/2605.23896)
[![C++](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![ROOT](https://img.shields.io/badge/ROOT-6.xx-red.svg)](https://root.cern/)

This repository contains the C++/ROOT source code supporting the benchmark validation presented in **Appendix A** of the manuscript:

> **A Stochastic Approach for Determining the Quark Confinement Potential of Charmonia**  
> Ahmet Bingül and Altuğ Özpineci  
> *arXiv:2605.23896 [hep-ph]*

---

## Overview

The code provides a numerical suite to solve the radial Schrödinger equation for heavy quarkonium ($q\bar{q}$) bound states using the **Numerov integration method** and extracts the spatial form of the confinement potential via a **Monte Carlo stochastic reconstruction algorithm**.

For simplicity and to isolate radial dynamics, spin-dependent potential terms are turned off ($L = S = 0$).

---

## Repository Structure

### 1. `StochasticSolver.h`
A supporting C++/ROOT class providing forward and inverse solving tools:
* **Forward Calculation:** Solves the Schrödinger equation for a target screened baseline interaction:
  ```math
  V_{\text{eff}}(r) = -\frac{4}{3}\frac{\alpha_s}{r} + \frac{\sigma}{\mu} \left(1 - e^{-\mu r}\right)

where the first term is the Coulomb interaction (one-gluon exchange) and the second term models color screening / string breaking.

### 2. `skeleton.cc`
The main ROOT macro program that implements the stochastic reconstruction algorithm using StochasticSolver.h.
Outputs: Energy eigenvalues True and calculated energy eigenvalues and discrete radial potential values. Such as
Starting Monte Carlo for =  5000000 trials and seed = 4507

    Trail# = 0   chi2 = 0.3392 seed = 4507
    True  Energies  : 0.7340 1.7617 2.4564 3.0006 3.4500 3.8313 4.1594 4.4443 4.6932 4.9131 5.1223 5.3514
    Calc. Energies  : 0.7584 1.8641 2.4693 2.9297 3.4061 3.7945 4.1541 4.5025 4.8510 5.1793 5.4518 5.6859
    True  Potential : 0.0000 1.0876 1.9781 2.7071 3.3040 3.7927 4.1928 4.5204 4.7886 5.0082 5.1880
    Calc. Potential : 0.0000 1.1558 2.2218 2.4261 3.2600 3.7521 4.2756 5.0150 5.4116 5.6103 5.9779
    . . .
    Trail# = 93   chi2 = 0.0556 seed = 4507
    True  Energies  : 0.7340 1.7617 2.4564 3.0006 3.4500 3.8313 4.1594 4.4443 4.6932 4.9131 5.1223 5.3514
    Calc. Energies  : 0.7029 1.7871 2.4893 2.9799 3.4051 3.7316 4.0189 4.3199 4.6258 4.8939 5.1330 5.3693
    True  Potential : 0.0000 1.0876 1.9781 2.7071 3.3040 3.7927 4.1928 4.5204 4.7886 5.0082 5.1880
    Calc. Potential : 0.0000 1.0630 2.1238 2.6526 3.2496 3.6215 3.8805 4.5005 4.9084 5.1332 5.9414


## Prerequisites & Usage
### Dependencies

   * C++ Compiler: C++17 or higher
   * ROOT Framework: CERN ROOT (v6.xx recommended)

### Running the Code

    Execute the main script within ROOT interactively or in standalone c++:
    Run interactively in ROOT
      $ root skeleton.cc

    Compile and run in standalone c++
      $ g++ skeleton.cc -o skeleton -O3 `root-config --cflags --glibs`
      $ ./skeleton

## Citation
If you use this code or method in your research, please cite

    @article{Bingul:2026stochastic,
     author        = "Bing{\"u}l, Ahmet and {\"O}zpineci, Altu{\u{g}}",
     title         = "{A Stochastic Approach for Determining the Quark Confinement Potential of Charmonia}",
     journal       = "arXiv preprint arXiv:2605.23896",
     archivePrefix = "arXiv",
     eprint        = "2605.23896",
     primaryClass  = "hep-ph",
     year          = "2026"
    }

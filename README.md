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

Output: Generates reference target energy eigenvalues $\{E_n^{\text{tar}}\}$, and exact and reconstructed potential profile values.

### 2. `skeleton.cc`
he main ROOT macro program that implements the stochastic reconstruction algorithm using StochasticSolver.h.

## Prerequisites & Usage
### Dependencies

   * C++ Compiler: C++17 or higher
   * ROOT Framework: CERN ROOT (v6.xx recommended)

### Running the Code

    Execute the main script within ROOT interactively or in standalone c++:
    Run interactively in ROOT
      $ root -l skeleton.cc

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

# GME_Cpp
Welcome to the c++ implementation of the GME method. 

*Please be aware that this is a project still under construction. Practical tutorials will be added soon! ;-)*

# Introduction

The GME (guided mode expansion) method solves the photonic band structure and diffraction (radiation or out-of-plane) losses of photonic crystal slabs. The underlying method is based on the expansion of the photonic crystal Bloch modes in the guided modes of the homogeneous slab (simple sine and cosine functions). The resulting matrix elements of the Maxwell operator (in the guided mode basis) are analytical and can be easily computed. This makes GME one of the fastest techniques to solve band structures of dispersionless and transparent photonic crystal slabs. The radiation (out-of-plane) losses is computed via the photonic equivalent of Fermi's golden rule (the so-called photonic golden rule), where the imaginary part of the frequency is estimated through time dependent perturbation theory. Please refer to the original work by L.C. Andreani and D. Gerace for full details on the method https://doi.org/10.1103/PhysRevB.73.235114

The c++ GME code takes advantage of the very efficient LAPACK subroutines to invert the dielectric matrix of the system and diagonalize the Maxwell operator. Therefore, *column-major* order is adopted for arranging the matrices throughout the code. Please note that, this is different from the c++ standard (*row-major* order). The GME method relies on the guided mode solutions of the effective homogeneous slab, which are computed numerically by solving the implicit equations appearing after imposing the electromagnetic boundary conditions of the slab. These implicit equations are solved using the GSL library via bracketing algorithms. For disordered dielectric structures, GSL is also invoked to generate random numbers with gaussian probability distribution.

# Compilation

Since the c++ GME code uses the GSL and LAPACK libraries to compute the dispersion and losses, they must be installed before compiling the GME. Assuming the default location of gsl (`/usr/local/include/gsl` for the included .h files and `/usr/lib` for the installed libraries) and lapack (`/usr/lib` for the installed libraries) the code can be compiled as follows:
```
g++ main.cpp -o main.out -lgsl -llapack -lblas -lm
```
where the `-l` options link the GME program with the corresponding libraries. If the location of the libraries is not the default location the `-L $library-location$` option must be used. The compilation creates the executable file `main.out` which will perform the calculation.

If you have access to the c++ Intel compiler I strongly suggest you to use it to compile the code. This will speedup the computation by a factor of 3 or 4 (in Intel architectures) with respect to the standard compilation with g++.

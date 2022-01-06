# GME_Cpp
This is the C++ implementation of the GME method

# 1   Introduction

GME (guided mode expansion) method computes the photonic dispersion and the diffraction losses for photonic crystal slabs. The underlying method is based on the expansion of the photonic crystal modes in the basis of guided modes associated to the homogeneous slab. To compute the diffraction losses, the time dependent perturbation theory is used to obtain an equivalent expression to the Fermi's golden rule (the photonic golden rule) and estimate the imaginary part of the mode frequency. The details of the method are presented in the work https://doi.org/10.1103/PhysRevB.73.235114

The c++ GME code takes advantage of the very efficient LAPACK subroutines to execute the principal matrix procedures, therefore, it is adopted the *column-major* order for arranging the matrices throughout the code. This is different from the c++ standard (row-major order). The heart of the GME method is composed of the guided mode solutions of the effective homogeneous slab, which are computed through the numerical solution of the implicit equations coming from the electromagnetic boundary conditions. These implicit equations are solved using the GSL library via bracketing algorithms. For some disordered dielectric structures, GSL is also used to generate random numbers with gaussian probability distribution.

# 2   Compilation

Since the c++ GME code uses the GSL and LAPACK libraries to compute the dispersion and losses, they must be installed before compiling the GME. Assuming the default location of gsl (`/usr/local/include/gsl` for the included .h files and `/usr/lib` for the installed libraries) and lapack (`/usr/lib` for the installed libraries) the code can be compiled as follows:
```
g++ main.cpp -o main.out -lgsl -llapack -lblas -lm
```
where the `-l` options link the GME program with the corresponding libraries. If the location of the libraries is not the default location the `-L $library-location$` option must be used. The compilation creates the executable file `main.out` which will perform the calculation.

# 3   Main file and libraries

The following files are included in this package:

| File  | Function  |
|-------|-----------|
|`main.cpp`|Input parameters and main control sequencies| 
|`eps.h`|Dielectric structures|
|`gme.h`|Main functions and routines of the GME method|
|`matrix.h`|Matrix computations|
|`compile.h`|Compiles the code|

**3.1   `main.cpp`**

The file `main.ccp` contains the input parameters and the main control sequences of the GME method. Generally speaking, the file is composed of *Initial instructions*, *Photonic dispersion*, *Diffraction losses*, *Structure*, *Basis and symmetry*, *Writing epsilon*, *K list*, *Initial checks and rountines of the code*, *Main loop for each k*

**_Initial instructions_**





<!-- <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1"> -->

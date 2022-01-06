# Documentation

The following files are included in the GME_Cpp package:

| File  | Function  |
|-------|-----------|
|`main.cpp`|Input parameters and main control sequencies| 
|`eps.h`|Dielectric structures|
|`gme.h`|Main functions and routines of the GME method|
|`matrix.h`|Matrix computations|
|`compile.h`|Compiles the code|

# `main.cpp`

The file `main.ccp` contains the input parameters and the main control sequences of the GME method. Generally speaking, the file is composed of *Initial instructions*, *Photonic dispersion*, *Diffraction losses*, *Structure*, *Basis and symmetry*, *Writing epsilon*, *K list*, *Initial checks and rountines of the code*, *Main loop for each k*

**1. Initial instructions**

`set_pre`<br/>
Numerical precision of the output.

`gme.sym`<br/>
If 1 assumes the hermiticity of the H matrix (Maxwell operator), this means, only the upper triangular part of H is built. If different from 1 the complete H matrix is built and the its hermiticity is checked.

`gme.psym`<br/>
If 1 sets inversion symmetry of the lattice, which makes the dielectric and H matrices real and symmetric. This option saves memory and makes the computation faster. If different from 1, the dielectric and H matrices are complex and hermitian. Note that, the code does not verify if the selected structure has indeed inversion symmetry.

**2. Photonic dispersion**

`gme.stat_i` and `gme.stat_f`<br/>
Initial and final state in the eigensolver. Values start from zero.

`gme.nkp`<br/>
Number of k points in the reciprocal space.

**3. Diffraction losses**

`gme.cimw`<br/>
If 1 computes the imaginary part of the frequency, and the eigenvectors are also computed by the eigensolver.

`gme.lstat_i` and `gme.lstat_f`<br/>
Initial and final state to compute the imaginary part of the frequency. Values start from zero.

`gme.klstat_i` and `gme.klstat_f`<br/>
Initial and final k points to compute the imaginary part of the frequency. Values start from zero.

**4. Fields**

`flim.cfds`<br/>
If 1 computes the fields, and the eigenvectors are also computed by the eigensolver.

`flim.stat_i` and `flim.stat_f`<br/>
Initial and final state to compute the fields. Values start from zero.

`flim.kstat_i` and `flim.kstat_f`<br/>
Initial and final k points to compute the fields. Values start from zero.

`flim.Dx`, `flim.Dy`, `flim.Dz`, `flim.Hx`, `flim.Hy` and `flim.Hz`<br/>
If 1 computes the corresponding field component.


<!-- <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1"> -->

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

`flim.xi`, `flim.xf`, `flim.dx`, `flim.yi`, `flim.yf`, `flim.dy`, `flim.zi`, `flim.zf` and `flim.dz`<br/>  
Initial value, final value and step size in the corresponding coordinate axes to evaluate the fields. These input variables define a 3D rectangular grid. 

**5. Structure**

`gme.e1`<br/>
Dielectric constant of the lower cladding.

`gme.e3`<br/>
Dielectric constant of the upper cladding.

`eps.`<br/>
Struct with the dielectric data structures. See below in the *eps.h* file

**6. Basis and symmetry**

`Gmaxvsnpw(eps,Gmax_initial,Gmax_final,Gmax_step)`<br/>
Writes the file *Gmax_vs_npw.dat* with the number of plane waves (first column) as a function of the cutoff values *Gmax* (second column) from `Gmax_initial` to `Gmax_final` with step size `Gmax_step`.

`gme.Gmax`<br/>
Cutoff value which truncates the plane wave basis. The truncated set of plane waves is defined by a circle in the reciprocal space, with radius `gme.Gmax` and centered at the origin.

`gme.Ng`<br/>
Number of guided modes.

`gme.sigm`<br/>
Symmetry of the eigenmodes to be computed with respect to the reflection operator <img src="https://render.githubusercontent.com/render/math?math=\hat{\sigma}_{xy}">. If 1 the even modes (TE-like) are computed, and if -1 the odd modes (TM-like) are computed. If `gme.e1` is different from `gme.e3` the value of `gme.sigm` is ignored.

**7. Write epsilon**

`epslim.weps`<br/>
If 1 the 2D Fourier expansion of the in-plane dielectric function is written.

`epslim.Gmax`<br/>
Cutoff used to write the Fourier expansion of epsilon. This is independent of `gme.Gmax`.

`epslim.xi`, `epslim.xf`, `epslim.dx`, `epslim.yi`, `epslim.yf` and `epslim.dy`<br/>
Initial value, final value and step size in the corresponding coordinate axes to evaluate the dielectric function. These input variables define a 2D rectangular grid.

**8. k list**

`gme.K1[0]`, `gme.K1[1]`, `gme.K2[0]`, `gme.K2[1]`, `gme.K3[0]` and `gme.K3[1]`<br/>
Defines a point (only `gme.K1`), a line (`gme.K1` and `gme.K2`), or a triangle (`gme.K1`, `gme.K2` and `gme.K3`) in the reciprocal space to compute the photonic eigenmodes.

`klist_point(gme,flim)`<br/>
The system is solved at the point `gme.K1` only. Here, `gme.nkp`, `gme.K2` and `gme.K3` are ignored.

`klist_line(gme,flim)`<br/>
The photonic dispersion is computed along the straight line between `gme.K1` and `gme.K2`. There are a total of `gme.nkp` points in this line (`gme.K1` and `gme.K2` included). Here `gme.K3` is ignored.

`klist_tri(gme,flim)`<br/>
The photonic dispersion is computed along the perimiter of the triangle with vertices `gme.K1`, `gme.K2` and `gme.K3`. There are `gme.nkp` points along each perimeter line (total of 3*`gme.nkp`). Points `gme.K1`, `gme.K2` and `gme.K3` are not included in the simulation.

**9. Initial routines of the program**

`initcheck(gme,flim)`<br/>
Checks the compatibility of the *stat_* variables.

`write_eps(eps,epslim)`<br/>
Writes the dielectric function.

`set_npw(gme,eps)`<br/>
Sets the number of the plane waves.

`alloc_Ggs(gme)`<br/>
Allocates que G and g arrays.

`set_G(gme)`<br/>
Sets the reciprocal lattice vectors <img src="https://render.githubusercontent.com/render/math?math=\vec{G}">

`minfo(gme)`<br/>
Writes the size of basis and states to be computed.

`set_epsm(gme,eps)`<br/>
Sets the average dielectric constants.

`set_Metaeps1(gme,eps)`, `set_Metaeps2(gme,eps)` and `set_Metaeps3(gme,eps)`<br/>
Allocates and computes the dielectric matrices in the corresponding three regions.

`alloc_TETM(gme)`<br/>
Allocates the arrays where the guided mode solutions will be stored.

**10. Main loop over the k list**
This is the main loop of the code, which iterate over the k list. The function `dnsgetvi(gme.lk,gme.k,kbz,gme.dimlk,5)` extract the components of the kbz-th k-vector from the array `gme.lk` (where all k points are stored) and copy them in the `gme.k` array. The latter is one used by the functions within the loop. The function `kmess(gme,kbz)` writes the current value of k in the format (kx,ky).

**_Computing the guided mode basis_**

`set_g(gme)`<br/>
Sets the g vetors.

`TETMstates(gme)`<br/>
Computes the guided mode states.

`dimcheck(gme,flim)`<br/>
Checks the compatibility of the *stat_* variables with the actual dimension of the H matrix.

**_Computing the photonic dispersion_**

`MHf(gme)`<br/>
Allocates and builts the H matrix.

`eigensolver(gme,flim)`<br/>
Computes the eigenmodes of the H matrix.

`writedis(gme)`<br/>
Writes the photonic dispersion in the file *dispersion.dat* with the format (kx,ky,k_brillouin,k_mag,w.a/2.pi.c), where k_brillouin is the corresponding value in the horizontal axis of the band diagram, and k_mag is the magnitude of the k vector over 2.pi (useful to plot the light line of the slab). Note that, the frequency output w.a/2.pi.c = a/lambda is dimesionless with 'a' and 'c' representing the lattice parameter and speed of light, respectively.

`writedis_1p(gme,p1)`<br/>
Writes the photonic dispersion in the file *dispersion_1p.dat* with one parameter using the format (p1,kx,ky,k_brillouin,k_mag,w.a/2.pi.c), where p1 is one free parameter to be defined by the user. All other quantities are defined as in `writedis(gme)`.

`writedis_2p(gme,p1,p2)`<br/>
Writes the photonic dispersion in the file *dispersion_1p.dat* with two parameters using the format (p1,p2,kx,ky,k_brillouin,k_mag,w.a/2.pi.c), where p1 and p2 are two free parameters to be defined by the user. All other quantities are defined as in `writedis(gme)`.

**_Computing the diffraction losses_**

`alloc_imw(gme)`<br/>
Allocates the array where the `gme.lstat_f-gme.lstat_i+1` values of Im{w} will be stored.

`computeimw(gme,eps)`<br/>
Computes the imaginary part of the frequency Im{w}.

`writeloss(gme)`<br/>
Writes the diffraction losses in the file *losses.dat* with the format (kx,ky,k_brillouin,lstate,w.a/2.pi.c,Im{w.a/2.pi.c}), where k_brillouin is the corresponding value in the horizontal axis of the band diagram, and lstate is the state for which the imaginary part of the frequency is calculated. Note that, as in the case of the photonic dispersion, the frequency output w.a/2.pi.c = a/lambda is dimesionless with 'a' and 'c' representing the lattice parameter and speed of light, respectively.

`writeloss_1p(gme,p1)`<br/>
Writes the diffraction losses in the file *losses_1p.dat* with one parameter using the format (p1,kx,ky,k_brillouin,lstate,w.a/2.pi.c,Im{w.a/2.pi.c}), where p1 is one free parameter to be defined by the user. All other quantities are defined as in `writeloss(gme)`.

`writeloss_2p(gme,p1,p2)`<br/>
Writes the diffraction losses in the file *losses_2p.dat* with two parameters using the format (p1,p2,kx,ky,k_brillouin,lstate,w.a/2.pi.c,Im{w.a/2.pi.c}), where p1 and p2 are two free parameters to be defined by the user. All other quantities are defined as in `writeloss(gme)`.

`free_imw(gme)`<br/>
Deallocates the array where the values of Im{w} are stored.

**_Computing fields_**











<!-- <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1"> -->

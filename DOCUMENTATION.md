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

The file `main.ccp` contains the input parameters and the main control sequences of the GME method. This file is broken down as follows: *Initial instructions*, *Photonic dispersion*, *Diffraction losses*, *Structure*, *Basis and symmetry*, *Writing epsilon*, *K list*, *Initial checks and rountines of the code*, *Main loop for each k*

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

`writefields(gme,flim,eps,kbz)`<br/>
Writes the field in the file *Fc-k_kbz-band_state.dat*, where *F* refers to the electromagnetic field (*D* or *H*), *c* refers to the field component (x, y or z), kbz refers to the kstate and state refers to the photonic mode state (the band number starting from zero). The format in the output file is (x,y,z,Re{Fc},Im{Fc}).

**_Memory deallocation_**

`free_MHeigen(gme)`<br/>
Deallocates all arrays associated to the eigenvalue problem.

`free_TETM(gme)`<br/>
Deallocates the arrays associated to the guided mode solutions.

`free_Metaeps1(gme)`, `free_Metaeps2(gme)` and `free_Metaeps3(gme)`<br/>
Deallocates the dielectric matrices in the corresponding three regions.

`free_Ggs(gme)`<br/>
Deallocates que G and g arrays.


# The `eps.h` library

The library *eps.h* contains the dielectric data structures and the main routines associated to the computation of the dielectric matrices. When using a given dielectric structure, it must be initialized through the constructor `init()`. Some of the structures, also must be deallocated before introducing any change in the geometry or material by means of the function `free()`, and they have to be initialized again. The following structures are currently included in the GME_Cpp package:


**<br/>`eps.recreg`**

Rectangular lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`lx`: Cell length along 'x'.<br/>
`ly`: Cell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexreg`**

Regular hexagonal lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.Thexreg`**

Regular hexagonal lattice of triangles.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`et`: dielectric constant of the triangles.<br/>
`L`: triangle side length.<br/>
`theta`: triangle rotation.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.rech1`**

H1 defect in a rectangular lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`lxp`: Basic cell length along 'x'.<br/>
`lyp`: Basic cell length along 'y'.<br/>
`lx`: Supercell length along 'x' (integer multiple of `lxp`).<br/>
`ly`: Supercell length along 'y' (integer multiple of `lyp`).<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexh1`**

H1 defect in a hexagonal lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`l`: Supercell parameter (integer number).<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexh1r`**

H1 defect in a hexagonal lattice of cylinders with rectangular cell.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`lx`: Supercell length along 'x'.<br/>
`ly`: Supercell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.recl2`**

L2 defect in a rectangular lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`lxp`: Basic cell length along 'x'.<br/>
`lyp`: Basic cell length along 'y'.<br/>
`lx`: Supercell length along 'x' (integer multiple of `lxp`).<br/>
`ly`: Supercell length along 'y' (integer multiple of `lyp`).<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexl2`**

L2 defect in a hexagonal lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`l`: Supercell parameter (integer number).<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexl2r`**

L2 defect in a hexagonal lattice of cylinders with rectangular cell.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`lx`: Supercell length along 'x'.<br/>
`ly`: Supercell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.recl3`**

L3 defect in a rectangular lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`s`: outward displacement of the closest lateral holes.<br/>
`rs`: radii of the closest lateral holes.<br/>
`lxp`: Basic cell length along 'x'.<br/>
`lyp`: Basic cell length along 'y'.<br/>
`lx`: Supercell length along 'x' (integer multiple of `lxp`).<br/>
`ly`: Supercell length along 'y' (integer multiple of `lyp`).<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexl3`**

L3 defect in a hexagonal lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`s`: outward displacement of the closest lateral holes.<br/>
`rs`: radii of the closest lateral holes.<br/>
`l`: Supercell parameter (integer number).<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexl3r`**

L3 defect in a hexagonal lattice of cylinders with rectangular cell.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`s`: outward displacement of the closest lateral holes.<br/>
`rs`: radii of the closest lateral holes.<br/>
`lx`: Supercell length along 'x'.<br/>
`ly`: Supercell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.hexl3sr`**

L3 defect with outward displacement of 5 closest lateral holes in a hexagonal lattice of cylinders with rectangular cell.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`s1`: outward displacement of the closest lateral holes.<br/>
`rs1`: radii of the closest lateral holes.<br/>
`s2`: outward displacement of the second closest lateral holes.<br/>
`rs2`: radii of the second closest lateral holes.<br/>
`s3`: outward displacement of the third closest lateral holes.<br/>
`rs3`: radii of the third closest lateral holes.<br/>
`s4`: outward displacement of the fourth closest lateral holes.<br/>
`rs4`: radii of the fourth closest lateral holes.<br/>
`s5`: outward displacement of the fifth closest lateral holes.<br/>
`rs5`: radii of the fifth closest lateral holes.<br/>
`lx`: Supercell length along 'x'.<br/>
`ly`: Supercell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.recwg`**

Waveguide defect in a rectangular lattice of cylinders.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`ax`: Basic cell length along 'x'.<br/>
`lyp`: Basic cell length along 'y'.<br/>
`ly`: Supercell length along 'y' (integer multiple of `lyp`).<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.Thexwgr`**

Waveguide defect in a hexagonal lattice of triangles with rectangular cell.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`et`: dielectric constant of the triangles.<br/>
`L`: triangle side length.<br/>
`theta`: triangle rotation.<br/>
`ly`: Supercell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.
`free()`: deallocates the structure


**<br/>`eps.Phexwgr`**

Waveguide defect in a hexagonal lattice of polygons (N sides) with rectangular cell

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`ep`: dielectric constant of the polygons.<br/>
`N`: number of sides.<br/>
`l`: polygon side length.<br/>
`th`: polygon rotation.<br/>
`ly`: Supercell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.
`free()`: deallocates the structure


**<br/>`eps.hexwgr`**

Waveguide defect in a hexagonal lattice of cylinders with rectangular cell.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`ly`: Supercell length along 'y'.<br/>

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.


**<br/>`eps.dishexl3r`**

L3 defect in a hexagonal lattice of cylinders with disorder (hole radii only) and rectangular cell.

**_Parameters_**<br/>
`es`: dielectric constant of the slab.<br/>
`er`: dielectric constant of the cylinders.<br/>
`r`: radii of the cylinders.<br/>
`s`: outward displacement of the closest lateral holes.<br/>
`rs`: radii of the closest lateral holes.<br/>
`lx`: Supercell length along 'x'.<br/>
`ly`: Supercell length along 'y'.<br/>
`sigma`: Variance of the gaussian distribution.

**_Constructor and functions_**<br/>
`init()`: initializes the dielectric structure.
`free()`: deallocates the structure


# The `gme.h` library

The library ``gme.h'' contains the main functions and routines of the GME method, corresponding to the photonic dispersion and diffraction losses. The functions associated to the electromagnetic fields and the root solver are also in this library.


# The `matrix.h` library

The library ``matrix.h'' contains the LAPACK matrix routines used by the GME method.


# The `compile.sh` library

The *compile.sh* is an executable file with the compilation command.


# Comments on units

Every length in the c++ GME code is scaled with respect to the lattice parameter 'a'. For example, when setting $d=0.5$, it must be interpreted as 'd/a=0.5'. The wavenumber 'k' is also scaled with respect to the lattice parameter as 'ak, therefore, when setting 'k=2.pi', it must be interpreted as ka=2.pi. The frequency magnitude used in the code is 'w.a/2.pi.c' or 'a/lambda', which is a dimensionless quantity. Therefore, when the code shows a frequency (or its imaginary part) value of '0.32', this means 'w.a/2.pi.c=0.32' so that the actual angular frequency is 'w=0.32.2.pi/a'. 

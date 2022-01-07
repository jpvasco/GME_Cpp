static int set_pre=14; // sets the precision of the output 
static double t,ttot; // time variables 

#include "matrix.h"
#include "eps.h"
#include "gme.h"



int main (){

ttot=clock();
cout.precision(set_pre);


gme.sym=1; // if 1 assumes the hermiticity of the H matrix	
gme.psym=1; // if 1 sets point symmetry 

// PHOTONIC DISPERSION
gme.stat_i=100;	gme.stat_f=100;	gme.nkp=10;	// initial state; final state; number of kpoints  

// DIFFRACTION LOSSES 
gme.cimw=0;	gme.lstat_i=100;	gme.lstat_f=100;	// if 1 computes losses; initial state; final state 
gme.klstat_i=0;	gme.klstat_f=11;			// initial kstate; final kstate

// FIELDS
flim.cfds=1;	flim.stat_i=100;	flim.stat_f=100;	// if 1 computes fields; initial state; final state
flim.kstat_i=11;	flim.kstat_f=11;			// initial kstate; final kstate
flim.Dx=1;	flim.Dy=1;	flim.Dz=0;	// if 1 write Ex; if 1 write Ey; if 1 write Ez 	
flim.Hx=0;	flim.Hy=0;	flim.Hz=0;	// if 1 write Hx; if 1 write Hy; if 1 write Hz

flim.xi=0.0;	flim.xf=0.0;	flim.dx=0.1;	// initial x; final x; step in x
flim.yi=0.0;	flim.yf=0.0;	flim.dy=0.1*sqrt(3);	// initial y; final y; step in y
flim.zi=0.0;	flim.zf=0.0;	flim.dz=0.1;	// initial z; final z; step in z

// STRUCTURE
gme.e1=1;	gme.e3=1;	gme.d=0.4615384;	// eps lower cladding; eps upper cladding; core thickness

eps.hexl3r.es=11.6281;	eps.hexl3r.er=1.0;	eps.hexl3r.r=0.25;	eps.hexl3r.s=0.15;
eps.hexl3r.rs=0.25*0.8;	eps.hexl3r.lx=10;	eps.hexl3r.ly=5*sqrt(3);
eps.hexl3r.init();

// BASIS AND SYMMETRY
Gmaxvsnpw(eps,0,100.0,0.1);			// query Gmax vs number of planewaves (eps,G_initial,G_final,G_step)
gme.Gmax=25.0;	gme.Ng=1;	gme.sigm=1;	// cutoff; guided modes; sigmaxy

// WRITE EPS
epslim.weps=1;	epslim.Gmax=20.0;		// if 1 writes eps pattern; cutoff in the eps expansion

epslim.xi=-5.0;	epslim.xf=5.0;	epslim.dx=0.1;	// initial x; final x; step in x
epslim.yi=-2.5*sqrt(3);	epslim.yf=2.5*sqrt(3);	epslim.dy=0.1;	// initial y; final y; step in y

// K LIST
gme.K1[0]=0.001;		gme.K1[1]=0.0;			
gme.K2[0]=M_PI/eps.hexl3r.lx;		gme.K2[1]=0.0;	
gme.K3[0]=2*M_PI/3.0;	gme.K3[1]=2*M_PI/sqrt(3);	

klist_line(gme,flim);



// INITIAL ROUTINES OF THE PROGRAM
initcheck(gme,flim);
write_eps(eps,epslim);
set_npw(gme,eps);
alloc_Ggs(gme);
set_G(gme);
minfo(gme);
set_epsm(gme,eps);
set_Metaeps1(gme,eps);
set_Metaeps2(gme,eps);
set_Metaeps3(gme,eps);
alloc_TETM(gme);
breakline(2);


// MAIN LOOP, FOR EACH K VECTOR
for (int kbz=0;kbz<gme.dimlk;kbz++){
 dnsgetvi(gme.lk,gme.k,kbz,gme.dimlk,5);
 kmess(gme,kbz);
 
 // guided mode states
 set_g(gme);
 TETMstates(gme);
 dimcheck(gme,flim);
 
 // photonic dispersion
 MHf(gme);
 eigensolver(gme,flim);
 writedis(gme);
 
 // diffraction losses
 alloc_imw(gme);
 computeimw(gme,eps);
 writeloss(gme);
 free_imw(gme);

 // fields
 writefields(gme,flim,eps,kbz);

 // fields with data
 //writefields_data(gme,flim,eps,kbz);

 // deallocation principal loop
 free_MHeigen(gme);

 breakline(2);
}

// deallocation initial routines
free_TETM(gme);

free_Meta1(gme);
free_Meta2(gme);
free_Meta3(gme);

free_Ggs(gme);

ttot=(clock()-ttot)/CLOCKS_PER_SEC;
cout<<"Total calculation done in "<<ttot<<" seconds\n";
cout.flush();
}

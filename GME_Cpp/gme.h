#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace std;

/*******************************************************************************/
/*****************************    GME c++ code     *****************************/
/*******************************************************************************/

struct st_gme { // variables and arrays used throughout the program
 int npw,Ng,lte,ltm,sym,psym,cimw,stat_i,stat_f,lstat_i,lstat_f,nkp,dimlk,klstat_i,klstat_f,sigm,Nm[2];
 double Gmax,e1,e3,d,epsm1,epsm2,epsm3,wk,imw,K1[2],K2[2],K3[2],k[5],enorm,hnorm;
 double *gsm,*eta1,*eta2,*eta3,*Meps1,*Meps2,*Meps3,*G,*lk,*gs,*solte,*soltm,*imwv,*MH,*eigenw,*eigenvec,*vk;
 complex<double> *ceta2,*cMeps2,*cMH,*ceigenvec,*cvk;
}gme;
 
/*######################################################################################################*/

/*******************************************************************************/
/****************************    Basic functions    ****************************/
/*******************************************************************************/

// Wave numbers 

double qq (double w, double g, double epsm2){
  return sqrt(epsm2*w*w-g*g);
}

double Xx (double w, double g, double epsm13){
  return sqrt(g*g-epsm13*w*w);
}

complex<double> qx (double w, double g, double epsm123){
  return sqrt(complex<double>(epsm123*w*w-g*g,0));
}

// Unitary vectors 


void eg (double *gv, double *egc){
 double norm;
 egc[0]=-gv[1];
 egc[1]=gv[0];
 norm=sqrt(egc[0]*egc[0]+egc[1]*egc[1]);
 egc[0]=egc[0]/norm;
 egc[1]=egc[1]/norm;
}

void gu (double *gv){
 double norm=sqrt(gv[0]*gv[0]+gv[1]*gv[1]);
 gv[0]=gv[0]/norm;
 gv[1]=gv[1]/norm;
}

/* Number of modes 
 * Writes in the 2-element array Nm the number of te and tm modes (te,tm) for a given symmetry sim
 *
 * sig+: sim=1	sig-: sim=-1	non-symmetric: sim=0 
 */

void Nmod (int sim, int Ng, int *Nm){
 if (Ng%2==0){
   Nm[0]=Ng/2;
   Nm[1]=Ng/2;
 } 
 else{
  if (sim==-1){
   Nm[0]=(Ng-1)/2;
   Nm[1]=(Ng+1)/2;
  }  
  if (sim==0 || sim==1){
   Nm[0]=(Ng+1)/2;
   Nm[1]=(Ng-1)/2;
  }
 }
}

// Photonic density of states 

double roj (double w, double g, double epsm1, double epsm3, int j){
 if (j==1){
  double arg1=w*w-g*g/epsm1;
  if (arg1<=0)
   return 0;
  else 
   return sqrt(epsm1)/(4*M_PI*sqrt(arg1));
 }
 if (j==3){
  double arg3=w*w-g*g/epsm3;
  if (arg3<=0)
   return 0;
  else 
   return sqrt(epsm3)/(4*M_PI*sqrt(arg3));
 }
}

/*******************************************************************************/
/***************************    eps gme functions    ***************************/
/*******************************************************************************/

void set_epsm (st_gme &gme, st_eps &eps){
 gme.epsm1=gme.e1;
 gme.epsm2=real(Ceps(0,0,eps));
 gme.epsm3=gme.e3;
}

void set_Metaeps1 (st_gme &gme, st_eps &eps){
 cout<<"Computing dielectric matrix: lower cladding, ";
 cout.flush();
 gme.Meps1=new double[gme.npw*(gme.npw+1)/2]();
 gme.eta1=new double[gme.npw*(gme.npw+1)/2]();
 dDcp(gme.epsm1,gme.Meps1,gme.npw);
 dDcp(1/gme.epsm1,gme.eta1,gme.npw);
 cout<<"done\n"; 
 cout.flush();
 delete[] gme.Meps1;
}

void set_Metaeps3 (st_gme &gme, st_eps &eps){
 cout<<"Computing dielectric matrix: upper cladding, ";
 cout.flush();
 gme.Meps3=new double[gme.npw*(gme.npw+1)/2]();
 gme.eta3=new double[gme.npw*(gme.npw+1)/2]();
 dDcp(gme.epsm3,gme.Meps3,gme.npw);
 dDcp(1/gme.epsm3,gme.eta3,gme.npw);
 cout<<"done\n"; 
 cout.flush();
 delete[] gme.Meps3;
}

void set_Metaeps2 (st_gme &gme, st_eps &eps){
 cout<<"Computing dielectric matrix: core, ";
 cout.flush();
 if (gme.psym==1){ // real case
  gme.Meps2=new double[gme.npw*(gme.npw+1)/2]();
  gme.eta2=new double[gme.npw*(gme.npw+1)/2]();
  dseps(eps,gme.G,gme.Meps2,gme.npw);
  dacp(gme.Meps2,gme.eta2,gme.npw*(gme.npw+1)/2);
  dspinv(gme.eta2,gme.npw);
  delete[] gme.Meps2;
 }
 else{ // complex case
  gme.cMeps2=new complex<double>[gme.npw*(gme.npw+1)/2]();
  gme.ceta2=new complex<double>[gme.npw*(gme.npw+1)/2]();
  zheps(eps,gme.G,gme.cMeps2,gme.npw);
  zacp(gme.cMeps2,gme.ceta2,gme.npw*(gme.npw+1)/2);
  zhpinv(gme.ceta2,gme.npw);
  delete[] gme.cMeps2;
 }
 cout<<"done\n"; 
 cout.flush();
}

// dealloc dielectric matrices

void free_Meta1 (st_gme &gme){
 delete[] gme.eta1;
}

void free_Meta3 (st_gme &gme){
 delete[] gme.eta3;
} 

void free_Meta2 (st_gme &gme){
 if (gme.psym==1) // real case
  delete[] gme.eta2;
 else //complex case
  delete[] gme.ceta2;
}

/*######################################################################################################*/

/*******************************************************************************/
/***************************    Implicit equations   ***************************/
/*******************************************************************************/

struct rootpar {double g,epsm1,epsm2,epsm3,d;}; // Parameters

/**************************************************************/
/**************    Symmetric case epsm1=epsm3    **************/
/**************************************************************/

/* Asymptotes of the function tan(qd/2). 
 *
 * We use them to define the bounded region in the bracketing algorithm.
 * first asymptote n=0, second n=1....
 */

double wasymsym (int n, double g, double epsm2, double d){
  return sqrt((2*n+1)*(2*n+1)*M_PI*M_PI/(d*d*epsm2)+g*g/epsm2);
}

// TE implicit equation with symmetry sig+

double FTEp (double w, void *params){
  
  struct rootpar *p=(struct rootpar *)params;
  
  double g=p->g,epsm1=p->epsm1,epsm2=p->epsm2,d=p->d,q;
  q=qq(w,g,epsm2);
  
  return tan(q*d*0.5)-Xx(w,g,epsm1)/q;
}

// TE implicit equation with symmetry sig-

double FTEm (double w, void *params){
  
  struct rootpar *p=(struct rootpar *)params;
  
  double g=p->g,epsm1=p->epsm1,epsm2=p->epsm2,d=p->d,q;
  q=qq(w,g,epsm2);
  
  return tan(q*d*0.5)+q/Xx(w,g,epsm1);
}

// TM implicit equation with symmetry sig+

double FTMp (double w, void *params){
  
  struct rootpar *p=(struct rootpar *)params;
  
  double g=p->g,epsm1=p->epsm1,epsm2=p->epsm2,d=p->d,q;
  q=qq(w,g,epsm2);
  
  return tan(q*d*0.5)+(q*epsm1)/(Xx(w,g,epsm1)*epsm2);
}

// TM implicit equation with symmetry sig-

double FTMm (double w, void *params){
  
  struct rootpar *p=(struct rootpar *)params;
  
  double g=p->g,epsm1=p->epsm1,epsm2=p->epsm2,d=p->d,q;
  q=qq(w,g,epsm2);
  
  return tan(q*d*0.5)-(Xx(w,g,epsm1)*epsm2)/(q*epsm1);
}

/**************************************************************/
/***********    Non symmetric case epsm1!=epsm3    ************/
/**************************************************************/

/* Asymptotes of the function tan(qd)+q(X1+X3)/(X1X3-q^2). 
 *
 * We use them to define the bounded region in the bracketing algorithm.
 * first asymptote n=0, second n=1....
 */

double wasymnonsymTE (int n, double g, double epsm1, double epsm2, double epsm3, double d){
  double sqrtn,sqrtnm1,wasysTE;
  if(n==0)
   sqrtnm1=0;
  else 
   sqrtnm1=sqrt((2*n-1)*(2*n-1)*M_PI*M_PI/(4*d*d*epsm2)+g*g/epsm2);
  
  sqrtn=sqrt((2*n+1)*(2*n+1)*M_PI*M_PI/(4*d*d*epsm2)+g*g/epsm2);  
  wasysTE=sqrt(g*g*(2*epsm2-epsm1-epsm3)/(epsm2*epsm2-epsm1*epsm3));
  
  if (sqrtn<wasysTE)
    return sqrtn;
  if (sqrtnm1<wasysTE && wasysTE<sqrtn)
    return wasysTE;
  if (wasysTE<sqrtnm1)
    return sqrtnm1;
  else{
    cout<<"an error has occurred in function wasymnonsymTE when the root-solver was trying to find the roots \n";
    abort();
  }
}

/* Asymptotes of the function tan(qd)+q*epsm2(X1*epsm3+X3epsm1)/(X1X3*epsm2^2-epsm1*epsm3*q^2). 
 *
 * We use them to define the bounded region in the bracketing algorithm.
 * first asymptote n=0, second n=1....
 */

double wasymnonsymTM (int n, double g, double epsm1, double epsm2, double epsm3, double d){
  double sqrtn,sqrtnm1,wasysTM,a,b,c;
  if(n==0)
   sqrtnm1=0;
  else 
   sqrtnm1=sqrt((2*n-1)*(2*n-1)*M_PI*M_PI/(4*d*d*epsm2)+g*g/epsm2);

  sqrtn=sqrt((2*n+1)*(2*n+1)*M_PI*M_PI/(4*d*d*epsm2)+g*g/epsm2);
  a=epsm2*epsm2*epsm2*epsm2*epsm1*epsm3-epsm1*epsm1*epsm2*epsm2*epsm3*epsm3;
  b=g*g*(2*epsm2*epsm1*epsm3*epsm1*epsm3-epsm2*epsm2*epsm2*epsm2*(epsm1+epsm3));
  c=g*g*g*g*(epsm2*epsm2*epsm2*epsm2-epsm1*epsm3*epsm1*epsm3);
  wasysTM=sqrt((-b-sqrt(b*b-4*a*c))/(2*a));
  
  if (sqrtn<wasysTM)
    return sqrtn;
  if (sqrtnm1<wasysTM && wasysTM<sqrtn)
    return wasysTM;
  if (wasysTM<sqrtnm1)
    return sqrtnm1;
  else{
    cout<<"an error has occurred in function wasymnonsymTM when the root-solver was trying to find the roots \n";
    abort();
  }

}

/* TE implicit equation */

double FTE (double w, void *params){
  
  struct rootpar *p=(struct rootpar *)params;
  
  double g=p->g,epsm1=p->epsm1,epsm2=p->epsm2,epsm3=p->epsm3,d=p->d,q,Xx1,Xx3;
  Xx1=Xx(w,g,epsm1);
  q=qq(w,g,epsm2);
  Xx3=Xx(w,g,epsm3);
  
  return tan(q*d)+q*(Xx1+Xx3)/(Xx1*Xx3-q*q);
}

/* TM implicit equation */

double FTM (double w, void *params){
  
  struct rootpar *p=(struct rootpar *)params;
  
  double g=p->g,epsm1=p->epsm1,epsm2=p->epsm2,epsm3=p->epsm3,d=p->d,q,Xx1,Xx3;
  Xx1=Xx(w,g,epsm1);
  q=qq(w,g,epsm2);
  Xx3=Xx(w,g,epsm3);
  
  return tan(q*d)+q*epsm2*(Xx1*epsm3+Xx3*epsm1)/(Xx1*Xx3*epsm2*epsm2-epsm1*epsm3*q*q);
}

/*######################################################################################################*/


/*******************************************************************************/
/******************************    Root solver    ******************************/
/*******************************************************************************/
/*
 * TE case: tem=1		TM case: tem=-1
 * TE-like: sigm=1		TM-like: sigm=-1
 * ********************************************
 * If epsm1!=epsm3 the variable sigm is ignored
 */


void RootSolver (double epsm1, double epsm2, double epsm3, double d, int tem, int sigm, double g, int nroot, double *roots){
  
  int status,n=0,rn=0,iter,maxiter=100;
  double wmax,wlo,whi,gg=g,r,dw=pow(10,-10),pre=0.000001,cond;
      
  if(gg<0)
    g=-gg;
  
//  ofstream Log ("Log.txt",ios_base::app); 
    
   //Initialization of the solver and the function.
   const gsl_root_fsolver_type *T;
   gsl_root_fsolver *s;
   T=gsl_root_fsolver_brent;
   s=gsl_root_fsolver_alloc(T);
   gsl_function F;
   struct rootpar params={g,epsm1,epsm2,epsm3,d};
   F.params=&params;
   
   //Choosing the gme implicit function.
   if (epsm1==epsm3){
     if (sigm==1){
       if (tem==1){
	   F.function=&FTEp;
//           Log<<"Solving for TE modes sigma=1\n";
       } 
       else if (tem==-1){
	   F.function=&FTMp;
//           Log<<"Solving for TM modes sigma=1\n";
       }
       else{
//	   Log<<"invalid tem param\n";
	   cout<<"invalid tem param\n";
	   abort();
       }
     }
     else if (sigm==-1){
           if (tem==1){
	       F.function=&FTEm;
//               Log<<"Solving for TE modes sigma=-1\n";
           }
           else if (tem==-1){
	       F.function=&FTMm;
//               Log<<"Solving for TM modes sigma=-1\n";
           }
	   else{
//	       Log<<"invalid tem param\n";
	       cout<<"invalid tem param\n";
	       abort();
	     }
     }
     else{
//       Log<<"invalid tem param\n";
       cout<<"invalid sigm param\n";
       abort();
     }
   }
   else{
     if (tem==1){
        F.function=&FTE;
//        Log<<"Solving for TE modes\n";
     }
     else if (tem==-1){
	F.function=&FTM;
//        Log<<"Solving for TM modes\n";
     }
      else{
//	Log<<"invalid tem param\n";
	cout<<"invalid tem param\n";
	abort();
      }
   }   
  
 wmax=g/sqrt(GSL_MAX(epsm1,epsm3))-dw;
 wlo=g/sqrt(epsm2)+dw;

// Log<<"Using the "<<gsl_root_fsolver_name(s)<<" method for g="<<gg<<"\nConvergence condition: |diff|<"<<pre<<"min(w_lower,w_upper)\n";  

 do {
  if (n>0){
    if (epsm1==epsm3)
      wlo=wasymsym(n-1,g,epsm2,d)+dw;
    else{
      if (tem==1)
        wlo=wasymnonsymTE(n-1,g,epsm1,epsm2,epsm3,d)+dw;
      else
        wlo=wasymnonsymTM(n-1,g,epsm1,epsm2,epsm3,d)+dw;
    }
  }
  
  if (epsm1==epsm3)
    whi=wasymsym(n,g,epsm2,d)-dw;
  else{
    if (tem==1)
      whi=wasymnonsymTE(n,g,epsm1,epsm2,epsm3,d)-dw;
    else
      whi=wasymnonsymTM(n,g,epsm1,epsm2,epsm3,d)-dw;
  }
  
  if (wmax<whi)
    whi=wmax;
  
  if (GSL_SIGN(GSL_FN_EVAL(&F,wlo))==GSL_SIGN(GSL_FN_EVAL(&F,whi))){
    n++;
    if (epsm1==epsm3)
    cond=wasymsym(n-1,g,epsm2,d);
    else{
      if (tem==1)
        cond=wasymnonsymTE(n-1,g,epsm1,epsm2,epsm3,d);
      else
        cond=wasymnonsymTM(n-1,g,epsm1,epsm2,epsm3,d);
    }
    continue;
  }
  
  gsl_root_fsolver_set(s,&F,wlo,whi);
 
//  Log<<"iteration  [ w_lower , w_upper ]  root  diff\n";

  //Iteration
  iter=0;
  do {
    iter++;
    status=gsl_root_fsolver_iterate(s);
    r=gsl_root_fsolver_root(s);
    wlo=gsl_root_fsolver_x_lower(s);
    whi=gsl_root_fsolver_x_upper(s);
    status=gsl_root_test_interval(wlo,whi,0,pre);
   
    if (status==GSL_SUCCESS){
//      Log<<"Converged:\n";
      if(gg<0)
      	roots[rn]=-r;
      else
     	roots[rn]=r;
     	rn++;
    }
   
//    if (gg<0)
//     Log<<iter<<"  "<<-whi<<"  "<<-wlo<<"  "<<-r<<"  "<<whi-wlo<<"\n";
//    else 
//     Log<<iter<<"  "<<wlo<<"  "<<whi<<"  "<<r<<"  "<<whi-wlo<<"\n"; 
    
//    if (status==GSL_SUCCESS)
//     Log<<"\n";
  }
  while (status==GSL_CONTINUE && iter<maxiter);
  if (iter==maxiter){
//    Log<<"convergence not achieved in the root-solver after "<<maxiter<<" iterations\n";
    cout<<"convergence not achieved in the root-solver after "<<maxiter<<" iterations\n";
    abort();
  }
 
  n++;
  if (epsm1==epsm3)
    cond=wasymsym(n-1,g,epsm2,d);
  else{
    if (tem==1)
      cond=wasymnonsymTE(n-1,g,epsm1,epsm2,epsm3,d);
    else
      cond=wasymnonsymTM(n-1,g,epsm1,epsm2,epsm3,d);
  }
 }
 while (cond<wmax && rn<nroot);

 gsl_root_fsolver_free(s);
// Log.close();
}

/*######################################################################################################*/

/*******************************************************************************/
/*********************   Gg's and basis of guided modes    *********************/
/*******************************************************************************/

// number of plane waves

void set_npw (st_gme &gme, st_eps &eps){
 Gi(eps,gme.G,1,gme.npw,gme.Gmax,0);
}

// alloc Gg's

void alloc_Ggs(st_gme &gme){
 gme.gs=new double[gme.npw*2];
 gme.gsm=new double[gme.npw];
 gme.G=new double[gme.npw*4];
}

// dealloc Gg's

void free_Ggs(st_gme &gme){
 delete[] gme.gs;
 delete[] gme.gsm;
 delete[] gme.G;
}

// reciprocal lattice vectors

void set_G (st_gme &gme){
 int n=gme.npw; // this n, which is the total number of rows of the G array must be different
                // from gme.npw since the Gi routine recalculates the value of gme.npw
 Gi(eps,gme.G,0,gme.npw,gme.Gmax,n);
}

// g vectors, g=k+G

void set_g (st_gme &gme){
 for (int i=0;i<gme.npw;i++){
  gme.gs[i]=gme.k[0]+gme.G[i];
  gme.gs[i+gme.npw]=gme.k[1]+gme.G[i+gme.npw];
  gme.gsm[i]=sqrt(gme.gs[i]*gme.gs[i]+gme.gs[i+gme.npw]*gme.gs[i+gme.npw]);
 }
}

// TE and TM states

/* The maximum number of TE and TM states are npw*Nm[0] and npw*Nm[1], respectively.
 * But for some values of k the number of states is smaller. Nevertheless,
 * solte and soltm are allocated and initialized to zero in the main program with 
 * these maximum dimensions, as wrte and wrtm in the function below; in this way, for 
 * such values of k, where the dimension is smaller, the "latest entries" of these 
 * arrays are zero, and the real dimension of the system is lte+ltm.
 * We ignore such zeros in the conditions wrte[p]!=0 and wrtm[p]!=0 below!
 */

// alloc TE and TM

void alloc_TETM(st_gme &gme){
 if (gme.epsm1==gme.epsm3) // number of TE and TM modes  
  Nmod(gme.sigm,gme.Ng,gme.Nm);
 else               
  Nmod(0,gme.Ng,gme.Nm);
 gme.solte=new double[gme.npw*gme.Nm[0]*2]();
 gme.soltm=new double[gme.npw*gme.Nm[1]*2]();
}

// dealloc TE and TM

void free_TETM(st_gme &gme){
 delete[] gme.solte;
 delete[] gme.soltm;
}

// states

void TETMstates (st_gme &gme){
 int p;
 double g,*wrte,*wrtm;
 gme.lte=0;
 gme.ltm=0;
 
 for (int i=0;i<gme.npw;i++){ // base of the states TE and TM 
  g=gme.gsm[i];
  
  if (gme.Nm[0]==0){ // eigen-frequencies of the homogeneous waveguide 
   wrte=new double[gme.Nm[0]]();
   wrtm=new double[gme.Nm[1]]();
   RootSolver(gme.epsm1,gme.epsm2,gme.epsm3,gme.d,-1,gme.sigm,g,gme.Nm[1],wrtm);
  } 
  else if (gme.Nm[1]==0){
   wrte=new double[gme.Nm[0]]();
   wrtm=new double[gme.Nm[1]]();
   RootSolver(gme.epsm1,gme.epsm2,gme.epsm3,gme.d,1,gme.sigm,g,gme.Nm[0],wrte);
  }
  else{
   wrte=new double[gme.Nm[0]]();
   wrtm=new double[gme.Nm[1]]();
   RootSolver(gme.epsm1,gme.epsm2,gme.epsm3,gme.d,1,gme.sigm,g,gme.Nm[0],wrte);  
   RootSolver(gme.epsm1,gme.epsm2,gme.epsm3,gme.d,-1,gme.sigm,g,gme.Nm[1],wrtm);
  }
  p=0;  
  while (p<gme.Nm[0]){  // making the set of TE states 
   if (wrte[p]!=0){
    gme.solte[gme.lte]=i;
    gme.solte[gme.lte+gme.npw*gme.Nm[0]]=wrte[p];
    gme.lte++;
   }
   else
    break;
   p++;
  }
  p=0; 
  while (p<gme.Nm[1]){  // making the set of TM states 
   if (wrtm[p]!=0){
    gme.soltm[gme.ltm]=i;
    gme.soltm[gme.ltm+gme.npw*gme.Nm[1]]=wrtm[p];
    gme.ltm++;
   }
   else
    break;
   p++;
  }
 }

 delete[] wrte;
 delete[] wrtm; 
}

/*######################################################################################################*/

/*******************************************************************************/
/***************************    Field coefficients    **************************/
/***************************   Photonic dispersion    **************************/
/*******************************************************************************/

// TE coefficients, tec={B1,B2,A2,A3}

void TEc (double w, double g, double epsm1, double epsm2, double epsm3, double d, complex<double> &B1, complex<double> &B2, complex<double> &A2, complex<double> &A3){
  double A2p2,B2p2,A3p,B1inv2,A2B2,q,q2,X1,X12,X3,X32,vcos,vsin,g2;
  q=qq(w,g,epsm2);
  q2=q*q;
  X1=Xx(w,g,epsm1);
  X12=X1*X1;
  X3=Xx(w,g,epsm3);
  X32=X3*X3;
  vcos=cos(q*d);
  vsin=sin(q*d);
  g2=g*g;
  A2p2=0.25*(q2+X12)/q2;
  B2p2=A2p2;
  A3p=0.5*(q*(X3-X1)*vcos+(q2+X1*X3)*vsin)/(q*X3);
  A2B2=0.25*(2*(q2-X12)*vcos+4*q*X1*vsin)/q2;
  B1inv2=0.5*((X12+g2)/X1+(X32+g2)*A3p*A3p/X3)+d*((g2+q2)*(A2p2+B2p2)+(g2-q2)*A2B2*vsin/(q*d));
  B1=1.0/sqrt(B1inv2);
  B2=(0.5*B1/q)*complex<double>(q,X1)*exp(complex<double>(0,-0.5*q*d)); 
  A2=(0.5*B1/q)*complex<double>(q,-X1)*exp(complex<double>(0,0.5*q*d)); 
  A3=B1*A3p;                                                            
}

// TM coefficients, tmc={D1,D2,C2,C3}

void TMc (double w, double g, double epsm1, double epsm2, double epsm3, double d, complex<double> &D1, complex<double> &D2, complex<double> &C2, complex<double> &C3){
  double C2p2,D2p2,C3p,D1inv2,C2D2,q,q2,X1,X12,X3,vcos,vsin;
  q=qq(w,g,epsm2)/epsm2;
  q2=q*q;
  X1=Xx(w,g,epsm1)/epsm1;
  X12=X1*X1;
  X3=Xx(w,g,epsm3)/epsm3;
  vcos=cos(q*d*epsm2);
  vsin=sin(q*d*epsm2);
  C2p2=0.25*(q2+X12)/q2;
  D2p2=C2p2;
  C3p=0.5*(q*(X3-X1)*vcos+(q2+X1*X3)*vsin)/(q*X3);
  C2D2=0.25*(2*(q2-X12)*vcos+4*q*X1*vsin)/q2;
  D1inv2=0.5*(1/(epsm1*X1)+C3p*C3p/(X3*epsm3))+d*(C2p2+D2p2+C2D2*vsin/(q*d*epsm2));
  D1=1.0/sqrt(D1inv2);
  D2=(0.5*D1/q)*complex<double>(q,X1)*exp(complex<double>(0,-0.5*q*d*epsm2)); 
  C2=(0.5*D1/q)*complex<double>(q,-X1)*exp(complex<double>(0,0.5*q*d*epsm2)); 
  C3=D1*C3p;                                                                 
}

/*******************************************************************************/
/***************************    Field coefficients    **************************/
/***************************    Diffraction losses    **************************/
/*******************************************************************************/

/* TEr coefficients, tecr={X1,X2,X3,W1,W2,W3}
 *
 * outgoing in the lower cladding: j=1
 * outgoing in the upper cladding: j=3
 */

void TEcr (double w, double g, double epsm1, double epsm2, double epsm3, double d, complex<double> &X1, complex<double> &X2, complex<double> &X3, complex<double> &W1, complex<double> &W2, complex<double> &W3, int j){
 complex<double> q1,q2,q3,expm,expp,Mte21_11,Mte21_12,Mte21_21,Mte21_22,Mte32_11,Mte32_12,Mte32_21,Mte32_22,Mte31_11,Mte31_12,Mte31_21;
 q1=qx(w,g,epsm1);
 q2=qx(w,g,epsm2);
 q3=qx(w,g,epsm3);
 expm=exp(complex<double>(0,-0.5*d)*q2);
 expp=exp(complex<double>(0,0.5*d)*q2);
 Mte21_11=0.5*(q2+q1)*expp/q2;
 Mte21_12=0.5*(q2-q1)*expp/q2;
 Mte21_21=0.5*(q2-q1)*expm/q2;
 Mte21_22=0.5*(q2+q1)*expm/q2;
 Mte32_11=0.5*(q3+q2)*expp/q3;
 Mte32_12=0.5*(q3-q2)*expm/q3;
 Mte32_21=0.5*(q3-q2)*expp/q3;
 Mte32_22=0.5*(q3+q2)*expm/q3;
 Mte31_11=Mte32_11*Mte21_11+Mte32_12*Mte21_21;
 Mte31_12=Mte32_11*Mte21_12+Mte32_12*Mte21_22;
 Mte31_21=Mte32_21*Mte21_11+Mte32_22*Mte21_21;
 if (j==1){
 X1=1/sqrt(epsm1);
 W1=-X1*Mte31_12/Mte31_11;
 W2=Mte21_11*W1+Mte21_12*X1;
 X2=Mte21_21*W1+Mte21_22*X1;
 W3=0.0;
 X3=Mte32_21*W2+Mte32_22*X2;
 }
 else if (j==3){
 W3=1/sqrt(epsm3);
 X3=W3*Mte31_21/Mte31_11;
 W1=W3/Mte31_11;
 X1=0.0;
 W2=Mte21_11*W1;
 X2=Mte21_21*W1;
 }
}

/* TMr coefficients, tecr={Z1,Z2,Z3,Y1,Y2,Y3}
 *
 * outgoing in the lower cladding: j=1
 * outgoing in the upper cladding: j=3
 */

void TMcr (double w, double g, double epsm1, double epsm2, double epsm3, double d, complex<double> &Z1, complex<double> &Z2, complex<double> &Z3, complex<double> &Y1, complex<double> &Y2, complex<double> &Y3, int j){
 complex<double> q1,q2,q3,expm,expp,Mtm21_11,Mtm21_12,Mtm21_21,Mtm21_22,Mtm32_11,Mtm32_12,Mtm32_21,Mtm32_22,Mtm31_11,Mtm31_12,Mtm31_21;
 q1=qx(w,g,epsm1)/epsm1;
 q2=qx(w,g,epsm2)/epsm2;
 q3=qx(w,g,epsm3)/epsm3;
 expm=exp(complex<double>(0,-0.5*d)*q2*epsm2);
 expp=exp(complex<double>(0,0.5*d)*q2*epsm2);
 Mtm21_11=0.5*(q2+q1)*expp/q2;
 Mtm21_12=0.5*(q2-q1)*expp/q2;
 Mtm21_21=0.5*(q2-q1)*expm/q2;
 Mtm21_22=0.5*(q2+q1)*expm/q2;
 Mtm32_11=0.5*(q3+q2)*expp/q3;
 Mtm32_12=0.5*(q3-q2)*expm/q3;
 Mtm32_21=0.5*(q3-q2)*expp/q3;
 Mtm32_22=0.5*(q3+q2)*expm/q3;
 Mtm31_11=Mtm32_11*Mtm21_11+Mtm32_12*Mtm21_21;
 Mtm31_12=Mtm32_11*Mtm21_12+Mtm32_12*Mtm21_22;
 Mtm31_21=Mtm32_21*Mtm21_11+Mtm32_22*Mtm21_21;
 if (j==1){
 Z1=1.0;                               
 Y1=-Z1*Mtm31_12/Mtm31_11;        
 Y2=Mtm21_11*Y1+Mtm21_12*Z1; 
 Z2=Mtm21_21*Y1+Mtm21_22*Z1; 
 Y3=0.0;                              
 Z3=Mtm32_21*Y2+Mtm32_22*Z2; 
 }
 else if (j==3){
 Y3=1.0;                           
 Z3=Y3*Mtm31_21/Mtm31_11;       
 Y1=Y3/Mtm31_11;                
 Z1=0.0;                              
 Y2=Mtm21_11*Y1;                
 Z2=Mtm21_21*Y1;
 }                  
}

/*######################################################################################################*/

/*******************************************************************************/
/******************************     z integrals    *****************************/
/*******************************************************************************/

// Photonic dispersion, iz={I1,I2-,I2+,I3} 

void Iz (double w, double g, double wp, double gp, double epsm1, double epsm2, double epsm3, double d, double &I1, double &I2m, double &I2p, double &I3){
 double q,qp;
 q=qq(w,g,epsm2);
 qp=qq(wp,gp,epsm2);
 I1=1.0/(Xx(w,g,epsm1)+Xx(wp,gp,epsm1));
 if (q==qp)
  I2m=d;
 else
  I2m=2.0*sin(0.5*(q-qp)*d)/(q-qp);
 I2p=2.0*sin(0.5*(q+qp)*d)/(q+qp);
 I3=1.0/(Xx(w,g,epsm3)+Xx(wp,gp,epsm3));
}
 
// Diffraction losses, izr={I1-,I1+,I2-,I2+,I3-,I3+} 

void Izr (double w, double g, double wp, double gp, double epsm1, double epsm2, double epsm3, double d, complex<double> &I1m, complex<double> &I1p, complex<double> &I2m, complex<double> &I2p, complex<double> &I3m, complex<double> &I3p){
 double q;
 complex<double> q2r;
 q=qq(w,g,epsm2);
 q2r=qx(wp,gp,epsm2);
 I1m=1.0/(Xx(w,g,epsm1)-complex<double>(0,1)*qx(wp,gp,epsm1));
 I1p=1.0/(Xx(w,g,epsm1)+complex<double>(0,1)*qx(wp,gp,epsm1));
 if (q==q2r)
  I2m=d;
 else
  I2m=2.0*sin(0.5*(q-q2r)*d)/(q-q2r);
 I2p=2.0*sin(0.5*(q+q2r)*d)/(q+q2r);
 I3m=1.0/(Xx(w,g,epsm3)-complex<double>(0,1)*qx(wp,gp,epsm3));
 I3p=1.0/(Xx(w,g,epsm3)+complex<double>(0,1)*qx(wp,gp,epsm3));
}

/*######################################################################################################*/

/*******************************************************************************/
/****************************   Hamiltonian Matrix  ****************************/
/**************************    Photonic dispersion    **************************/
/*******************************************************************************/  

void MHf (st_gme &gme){ 
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *eta1d=(double *)gme.eta1,*eta3d=(double *)gme.eta3,*solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d;
 int dim=ltm+lte;
 double w,wp,g,gp,q,qp,degegp,I1,I2m,I2p,I3,eta1v,eta2vd,eta3v,factor,*MHd,*eta2d;
 complex<double> A2,A3,B1,B2,C2,C3,D1,D2,A2p,A3p,B1p,B2p,C2p,C3p,D1p,D2p,eta2vc,*MHc,*eta2c;

 if (gme.sym==1)
  cout<<"Assuming the hermiticity of the hamiltonian matrix\n";
 else
  cout<<"Building the complete hamiltonian matrix\n";
 cout<<"dimension of the hamiltonian matrix: "<<ltm+lte<<"\n";
 cout.flush();
 t=clock();

 if (gme.psym==1){ // real case
  gme.MH=new double[(lte+ltm)*(lte+ltm)]();
  MHd=(double *)gme.MH;
  eta2d=(double *)gme.eta2;

  //Assuming the symmetry of the hamiltonian matrix

  for (int i=0;i<ltm;i++){  // TM-TM block matrix  
   for (int j=i;j<ltm;j++){   
    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    q=qq(w,g,epsm2);
    wp=soltm[j+gme.npw*gme.Nm[1]];
    gp=gsm[int(soltm[j])];
    qp=qq(wp,gp,epsm2);  
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
    degegp=(gs[int(soltm[i])]*gs[int(soltm[j])]+gs[int(soltm[i])+npw]*gs[int(soltm[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>soltm[j])
     eta2vd=eta2d[iju(soltm[j],soltm[i])];
    else
     eta2vd=eta2d[iju(soltm[i],soltm[j])];
    MHd[i+j*dim]=real(eta2vd*((conj(C2)*C2p+conj(D2)*D2p)*(q*qp*degegp+g*gp)*I2m+(conj(C2)*D2p+conj(D2)*C2p)*(-q*qp*degegp+g*gp)*I2p));

    if (soltm[i]==soltm[j]){
     eta1v=eta1d[iju(soltm[i],soltm[j])];
     eta3v=eta3d[iju(soltm[i],soltm[j])];
     MHd[i+j*dim]=real(MHd[i+j*dim]+eta1v*conj(D1)*D1p*(Xx(w,g,epsm1)*Xx(wp,gp,epsm1)*degegp+g*gp)*I1+eta3v*conj(C3)*C3p*(Xx(w,g,epsm3)*Xx(wp,gp,epsm3)*degegp+g*gp)*I3);
    }
   }
  }

  for (int i=0;i<lte;i++){  // TE-TE block matrix 
   for (int j=i;j<lte;j++){ 
    w=solte[i+gme.npw*gme.Nm[0]];
    g=gsm[int(solte[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(gs[int(solte[i])]*gs[int(solte[j])]+gs[int(solte[i])+npw]*gs[int(solte[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (solte[i]>solte[j])
     eta2vd=eta2d[iju(solte[j],solte[i])];
    else
     eta2vd=eta2d[iju(solte[i],solte[j])];
    factor=w*w*wp*wp*degegp;
    MHd[i+ltm+(j+ltm)*dim]=real(factor*epsm2*epsm2*eta2vd*((conj(A2)*A2p+conj(B2)*B2p)*I2m+(conj(A2)*B2p+conj(B2)*A2p)*I2p)); 

    if (solte[i]==solte[j]){
     eta1v=eta1d[iju(solte[i],solte[j])];
     eta3v=eta3d[iju(solte[i],solte[j])];
     MHd[i+ltm+(j+ltm)*dim]=real(MHd[i+ltm+(j+ltm)*dim]+factor*(epsm1*epsm1*eta1v*conj(B1)*B1p*I1+epsm3*epsm3*eta3v*conj(A3)*A3p*I3));
    }
   }  
  }

  for (int i=0;i<ltm;i++){  // TM-TE block matrix  
   for (int j=0;j<lte;j++){   
    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(-gs[int(soltm[i])]*gs[int(solte[j])+npw]+gs[int(solte[j])]*gs[int(soltm[i])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>solte[j])
     eta2vd=eta2d[iju(solte[j],soltm[i])];
    else
     eta2vd=eta2d[iju(soltm[i],solte[j])];
    factor=wp*wp*degegp;
    MHd[i+(j+ltm)*dim]=real(factor*complex<double>(0,-epsm2)*eta2vd*qq(w,g,epsm2)*((-conj(C2)*A2p+conj(D2)*B2p)*I2m+(conj(D2)*A2p-conj(C2)*B2p)*I2p));

    if (soltm[i]==solte[j]){
     eta1v=eta1d[iju(soltm[i],solte[j])];
     eta3v=eta3d[iju(soltm[i],solte[j])];
     MHd[i+(j+ltm)*dim]=real(MHd[i+(j+ltm)*dim]+factor*(-epsm1*eta1v*conj(D1)*B1p*Xx(w,g,epsm1)*I1+epsm3*eta3v*conj(C3)*A3p*Xx(w,g,epsm3)*I3));
    }
   }
  }

  // The complete hamiltonian matrix is made 

  if (gme.sym!=1){

   for (int i=0;i<ltm;i++){  // TM-TM block matrix 
    for (int j=0;j<i;j++){   
     w=soltm[i+gme.npw*gme.Nm[1]];
     g=gsm[int(soltm[i])];
     q=qq(w,g,epsm2);
     wp=soltm[j+gme.npw*gme.Nm[1]];
     gp=gsm[int(soltm[j])];
     qp=qq(wp,gp,epsm2);  
     TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
     TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
     degegp=(gs[int(soltm[i])]*gs[int(soltm[j])]+gs[int(soltm[i])+npw]*gs[int(soltm[j])+npw])/(g*gp);
     Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
     if (soltm[i]>soltm[j])
      eta2vd=eta2d[iju(soltm[j],soltm[i])];
     else
      eta2vd=eta2d[iju(soltm[i],soltm[j])];
     MHd[i+j*dim]=real(eta2vd*((conj(C2)*C2p+conj(D2)*D2p)*(q*qp*degegp+g*gp)*I2m+(conj(C2)*D2p+conj(D2)*C2p)*(-q*qp*degegp+g*gp)*I2p));

     if (soltm[i]==soltm[j]){
      eta1v=eta1d[iju(soltm[i],soltm[j])];
      eta3v=eta3d[iju(soltm[i],soltm[j])];
      MHd[i+j*dim]=real(MHd[i+j*dim]+eta1v*conj(D1)*D1p*(Xx(w,g,epsm1)*Xx(wp,gp,epsm1)*degegp+g*gp)*I1+eta3v*conj(C3)*C3p*(Xx(w,g,epsm3)*Xx(wp,gp,epsm3)*degegp+g*gp)*I3);
     }
    }
   }

   for (int i=0;i<lte;i++){  // TE-TE block matrix  
    for (int j=0;j<i;j++){ 
     w=solte[i+gme.npw*gme.Nm[0]];
     g=gsm[int(solte[i])];
     wp=solte[j+gme.npw*gme.Nm[0]];
     gp=gsm[int(solte[j])];
     TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
     TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
     degegp=(gs[int(solte[i])]*gs[int(solte[j])]+gs[int(solte[i])+npw]*gs[int(solte[j])+npw])/(g*gp);
     Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
     if (solte[i]>solte[j])
      eta2vd=eta2d[iju(solte[j],solte[i])];
     else
      eta2vd=eta2d[iju(solte[i],solte[j])];
     factor=w*w*wp*wp*degegp;
     MHd[i+ltm+(j+ltm)*dim]=real(factor*epsm2*epsm2*eta2vd*((conj(A2)*A2p+conj(B2)*B2p)*I2m+(conj(A2)*B2p+conj(B2)*A2p)*I2p));

     if (solte[i]==solte[j]){
      eta1v=eta1d[iju(solte[i],solte[j])];
      eta3v=eta3d[iju(solte[i],solte[j])];
      MHd[i+ltm+(j+ltm)*dim]=real(MHd[i+ltm+(j+ltm)*dim]+factor*(epsm1*epsm1*eta1v*conj(B1)*B1p*I1+epsm3*epsm3*eta3v*conj(A3)*A3p*I3));
     }
    }
   }

   for (int i=0;i<lte;i++){  // TE-TM block matrix  
    for (int j=0;j<ltm;j++){   
     w=solte[i+gme.npw*gme.Nm[0]];
     g=gsm[int(solte[i])];
     wp=soltm[j+gme.npw*gme.Nm[1]];
     gp=gsm[int(soltm[j])];
     TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
     TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
     degegp=(-gs[int(soltm[j])]*gs[int(solte[i])+npw]+gs[int(solte[i])]*gs[int(soltm[j])+npw])/(g*gp);
     Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
     if (solte[i]>soltm[j])
      eta2vd=eta2d[iju(soltm[j],solte[i])];
     else
      eta2vd=eta2d[iju(solte[i],soltm[j])];
     factor=w*w*degegp;
     MHd[i+ltm+j*dim]=real(factor*complex<double>(0,epsm2)*eta2vd*qq(wp,gp,epsm2)*((-conj(A2)*C2p+conj(B2)*D2p)*I2m+(conj(A2)*D2p-conj(B2)*C2p)*I2p));

     if (solte[i]==soltm[j]){
      eta1v=eta1d[iju(solte[i],soltm[j])];
      eta3v=eta3d[iju(solte[i],soltm[j])];
      MHd[i+ltm+j*dim]=real(MHd[i+ltm+j*dim]+factor*(-epsm1*eta1v*conj(B1)*D1p*Xx(wp,gp,epsm1)*I1+epsm3*eta3v*conj(A3)*C3p*Xx(wp,gp,epsm3)*I3));
     }
    }
   }

   if (checksM(MHd,lte+ltm)==0){
    cout<<"the hamiltonian matrix is not hermitian\n";
    cout.flush();
    abort();
   }
   cout<<"hermiticity of the hamiltonian matrix: ok\n";
   cout.flush();
  }
 }
 else{ // complex case
  gme.cMH=new complex<double>[(lte+ltm)*(lte+ltm)]();
  MHc=(complex<double> *)gme.cMH;
  eta2c=(complex<double> *)gme.ceta2;

  //Assuming the hermiticity of the hamiltonian matrix

  for (int i=0;i<ltm;i++){  // TM-TM block matrix 
   for (int j=i;j<ltm;j++){   
    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    q=qq(w,g,epsm2);
    wp=soltm[j+gme.npw*gme.Nm[1]];
    gp=gsm[int(soltm[j])];
    qp=qq(wp,gp,epsm2);  
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
    degegp=(gs[int(soltm[i])]*gs[int(soltm[j])]+gs[int(soltm[i])+npw]*gs[int(soltm[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>soltm[j])
     eta2vc=conj(eta2c[iju(soltm[j],soltm[i])]);
    else
     eta2vc=eta2c[iju(soltm[i],soltm[j])];
    MHc[i+j*dim]=eta2vc*((conj(C2)*C2p+conj(D2)*D2p)*(q*qp*degegp+g*gp)*I2m+(conj(C2)*D2p+conj(D2)*C2p)*(-q*qp*degegp+g*gp)*I2p);

    if (soltm[i]==soltm[j]){
     eta1v=eta1d[iju(soltm[i],soltm[j])];
     eta3v=eta3d[iju(soltm[i],soltm[j])];
     MHc[i+j*dim]=MHc[i+j*dim]+eta1v*conj(D1)*D1p*(Xx(w,g,epsm1)*Xx(wp,gp,epsm1)*degegp+g*gp)*I1+eta3v*conj(C3)*C3p*(Xx(w,g,epsm3)*Xx(wp,gp,epsm3)*degegp+g*gp)*I3;
    }
   }
  }

  for (int i=0;i<lte;i++){  // TE-TE block matrix  
   for (int j=i;j<lte;j++){ 
    w=solte[i+gme.npw*gme.Nm[0]];
    g=gsm[int(solte[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(gs[int(solte[i])]*gs[int(solte[j])]+gs[int(solte[i])+npw]*gs[int(solte[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (solte[i]>solte[j])
     eta2vc=conj(eta2c[iju(solte[j],solte[i])]);
    else
     eta2vc=eta2c[iju(solte[i],solte[j])];
    factor=w*w*wp*wp*degegp;
    MHc[i+ltm+(j+ltm)*dim]=factor*epsm2*epsm2*eta2vc*((conj(A2)*A2p+conj(B2)*B2p)*I2m+(conj(A2)*B2p+conj(B2)*A2p)*I2p);

    if (solte[i]==solte[j]){
     eta1v=eta1d[iju(solte[i],solte[j])];
     eta3v=eta3d[iju(solte[i],solte[j])];
     MHc[i+ltm+(j+ltm)*dim]=MHc[i+ltm+(j+ltm)*dim]+factor*(epsm1*epsm1*eta1v*conj(B1)*B1p*I1+epsm3*epsm3*eta3v*conj(A3)*A3p*I3);
    }
   }
  }

  for (int i=0;i<ltm;i++){  // TM-TE block matrix  
   for (int j=0;j<lte;j++){   
    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(-gs[int(soltm[i])]*gs[int(solte[j])+npw]+gs[int(solte[j])]*gs[int(soltm[i])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>solte[j])
     eta2vc=conj(eta2c[iju(solte[j],soltm[i])]);
    else
     eta2vc=eta2c[iju(soltm[i],solte[j])];
    factor=wp*wp*degegp;
    MHc[i+(j+ltm)*dim]=factor*complex<double>(0,-epsm2)*eta2vc*qq(w,g,epsm2)*((-conj(C2)*A2p+conj(D2)*B2p)*I2m+(conj(D2)*A2p-conj(C2)*B2p)*I2p);

    if (soltm[i]==solte[j]){
     eta1v=eta1d[iju(soltm[i],solte[j])];
     eta3v=eta3d[iju(soltm[i],solte[j])];
     MHc[i+(j+ltm)*dim]=MHc[i+(j+ltm)*dim]+factor*(-epsm1*eta1v*conj(D1)*B1p*Xx(w,g,epsm1)*I1+epsm3*eta3v*conj(C3)*A3p*Xx(w,g,epsm3)*I3);
    }
   }
  }

  // The complete hamiltonian matrix is made 

  if (gme.sym!=1){

   for (int i=0;i<ltm;i++){  // TM-TM block matrix 
    for (int j=0;j<i;j++){   
     w=soltm[i+gme.npw*gme.Nm[1]];
     g=gsm[int(soltm[i])];
     q=qq(w,g,epsm2);
     wp=soltm[j+gme.npw*gme.Nm[1]];
     gp=gsm[int(soltm[j])];
     qp=qq(wp,gp,epsm2);  
     TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
     TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
     degegp=(gs[int(soltm[i])]*gs[int(soltm[j])]+gs[int(soltm[i])+npw]*gs[int(soltm[j])+npw])/(g*gp);
     Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
     if (soltm[i]>soltm[j])
      eta2vc=conj(eta2c[iju(soltm[j],soltm[i])]);
     else
      eta2vc=eta2c[iju(soltm[i],soltm[j])];
     MHc[i+j*dim]=eta2vc*((conj(C2)*C2p+conj(D2)*D2p)*(q*qp*degegp+g*gp)*I2m+(conj(C2)*D2p+conj(D2)*C2p)*(-q*qp*degegp+g*gp)*I2p);

     if (soltm[i]==soltm[j]){
      eta1v=eta1d[iju(soltm[i],soltm[j])];
      eta3v=eta3d[iju(soltm[i],soltm[j])];
      MHc[i+j*dim]=MHc[i+j*dim]+eta1v*conj(D1)*D1p*(Xx(w,g,epsm1)*Xx(wp,gp,epsm1)*degegp+g*gp)*I1+eta3v*conj(C3)*C3p*(Xx(w,g,epsm3)*Xx(wp,gp,epsm3)*degegp+g*gp)*I3;
     }
    }
   }

   for (int i=0;i<lte;i++){  // TE-TE block matrix  
    for (int j=0;j<i;j++){ 
     w=solte[i+gme.npw*gme.Nm[0]];
     g=gsm[int(solte[i])];
     wp=solte[j+gme.npw*gme.Nm[0]];
     gp=gsm[int(solte[j])];
     TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
     TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
     degegp=(gs[int(solte[i])]*gs[int(solte[j])]+gs[int(solte[i])+npw]*gs[int(solte[j])+npw])/(g*gp);
     Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
     if (solte[i]>solte[j])
      eta2vc=conj(eta2c[iju(solte[j],solte[i])]);
     else 
      eta2vc=eta2c[iju(solte[i],solte[j])];
     factor=w*w*wp*wp*degegp;
     MHc[i+ltm+(j+ltm)*dim]=factor*epsm2*epsm2*eta2vc*((conj(A2)*A2p+conj(B2)*B2p)*I2m+(conj(A2)*B2p+conj(B2)*A2p)*I2p);

     if (solte[i]==solte[j]){
      eta1v=eta1d[iju(solte[i],solte[j])];
      eta3v=eta3d[iju(solte[i],solte[j])];
      MHc[i+ltm+(j+ltm)*dim]=MHc[i+ltm+(j+ltm)*dim]+factor*(epsm1*epsm1*eta1v*conj(B1)*B1p*I1+epsm3*epsm3*eta3v*conj(A3)*A3p*I3);
     }
    }
   }

   for (int i=0;i<lte;i++){  // TE-TM block matrix 
    for (int j=0;j<ltm;j++){   
     w=solte[i+gme.npw*gme.Nm[0]];
     g=gsm[int(solte[i])];
     wp=soltm[j+gme.npw*gme.Nm[1]];
     gp=gsm[int(soltm[j])];
     TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
     TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
     degegp=(-gs[int(soltm[j])]*gs[int(solte[i])+npw]+gs[int(solte[i])]*gs[int(soltm[j])+npw])/(g*gp);
     Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
     if (solte[i]>soltm[j])
      eta2vc=conj(eta2c[iju(soltm[j],solte[i])]);
     else
      eta2vc=eta2c[iju(solte[i],soltm[j])];
     factor=w*w*degegp;
     MHc[i+ltm+j*dim]=factor*complex<double>(0,epsm2)*eta2vc*qq(wp,gp,epsm2)*((-conj(A2)*C2p+conj(B2)*D2p)*I2m+(conj(A2)*D2p-conj(B2)*C2p)*I2p);

     if (solte[i]==soltm[j]){
      eta1v=eta1d[iju(solte[i],soltm[j])];
      eta3v=eta3d[iju(solte[i],soltm[j])];
      MHc[i+ltm+j*dim]=MHc[i+ltm+j*dim]+factor*(-epsm1*eta1v*conj(B1)*D1p*Xx(wp,gp,epsm1)*I1+epsm3*eta3v*conj(A3)*C3p*Xx(wp,gp,epsm3)*I3);
     }
    }
   }

   if (checkhM(MHc,lte+ltm)==0){
    cout<<"the hamiltonian matrix is not hermitian\n";
    cout.flush();
    abort();
   }
   cout<<"hermiticity of the hamiltonian matrix: ok\n";
   cout.flush();
  }
 }

 t=(clock()-t)/CLOCKS_PER_SEC;
 cout<<"hamiltonian matrix made in "<<t<<" seconds\n";
 cout.flush();
}

/*######################################################################################################*/

/*******************************************************************************/
/***************************   Imaginary part of w   ***************************/
/**************************    Diffraction losses    ***************************/
/*******************************************************************************/

void Imw (st_gme &gme){
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *eta1d=(double *)gme.eta1,*eta3d=(double *)gme.eta3,wk=gme.wk,*solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d;
 int dim=ltm+lte;
 double w,wp,g,gp,Xx1,Xx3,q,degegp,simw2,eta1v,eta2vd,eta3v,*eta2d,*vkd;
 complex<double> sradte,sradtm,tem,q1r,q2r,q3r,Z1,Z2,Z3,Y1,Y2,Y3,X1,X2,X3,W1,W2,W3,A2,A3,B1,B2,C2,C3,D1,D2,I1m,I1p,I2m,I2p,I3m,I3p,eta2vc,*eta2c,*vkc;

 simw2=0;

 if (gme.psym==1){ // real case
  eta2d=(double *)gme.eta2;
  vkd=(double *)gme.vk;
 
  for (int j=1;j<4;j=j+2){ // computing Im(w^2) 
   for (int m=0;m<npw;m++){
    sradte=complex<double>(0,0);
    sradtm=complex<double>(0,0);

    if (roj(wk,gsm[m],epsm1,epsm3,j)!=0){
   //wp=wk
     gp=gsm[m];
     q1r=qx(wk,gp,epsm1);
     q2r=qx(wk,gp,epsm2);
     q3r=qx(wk,gp,epsm3);
     TMcr(wk,gp,epsm1,epsm2,epsm3,d,Z1,Z2,Z3,Y1,Y2,Y3,j);
     TEcr(wk,gp,epsm1,epsm2,epsm3,d,X1,X2,X3,W1,W2,W3,j);
    
     for (int i=0;i<ltm;i++) { // TM-TM
      w=soltm[i+gme.npw*gme.Nm[1]];
      g=gsm[int(soltm[i])];
      Xx1=Xx(w,g,epsm1);
      Xx3=Xx(w,g,epsm3);
      q=qq(w,g,epsm2);
      TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(gs[int(soltm[i])]*gs[m]+gs[int(soltm[i])+npw]*gs[m+npw])/(g*gp);
      if (soltm[i]>m)
       eta2vd=eta2d[iju(m,soltm[i])];
      else
       eta2vd=eta2d[iju(soltm[i],m)];
      tem=eta2vd*((conj(C2)*Y2+conj(D2)*Z2)*(g*gp+q*q2r*degegp)*I2m+(conj(C2)*Z2+conj(D2)*Y2)*(g*gp-q*q2r*degegp)*I2p);

      if (soltm[i]==m){
       eta1v=eta1d[iju(soltm[i],m)];
       eta3v=eta3d[iju(soltm[i],m)];
       tem=tem+eta1v*conj(D1)*((g*gp+complex<double>(0,1)*Xx1*q1r*degegp)*Y1*I1p+(g*gp-complex<double>(0,1)*Xx1*q1r*degegp)*Z1*I1m)+eta3v*conj(C3)*((g*gp-complex<double>(0,1)*Xx3*q3r*degegp)*Y3*I3m+(g*gp+complex<double>(0,1)*Xx3*q3r*degegp)*Z3*I3p);
      }

      sradtm=sradtm+vkd[i]*tem;
     }

     for (int i=0;i<lte;i++) { // TE-TM
      w=solte[i+gme.npw*gme.Nm[0]];
      g=gsm[int(solte[i])];
      TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(-gs[m]*gs[int(solte[i])+npw]+gs[int(solte[i])]*gs[m+npw])/(g*gp);
      if (solte[i]>m)
       eta2vd=eta2d[iju(m,solte[i])];
      else
       eta2vd=eta2d[iju(solte[i],m)];
      tem=epsm2*eta2vd*q2r*((-conj(A2)*Y2+conj(B2)*Z2)*I2m+(conj(A2)*Z2-conj(B2)*Y2)*I2p);

      if (solte[i]==m){
       eta1v=eta1d[iju(solte[i],m)];
       eta3v=eta3d[iju(solte[i],m)];
       tem=tem+epsm1*eta1v*q1r*conj(B1)*(-Y1*I1p+Z1*I1m)+epsm3*eta3v*q3r*conj(A3)*(-Y3*I3m+Z3*I3p);
      }
      tem=complex<double>(0,1)*w*w*degegp*tem;

      sradtm=sradtm+vkd[i+ltm]*tem;
     }

     for (int i=0;i<lte;i++) { // TE-TE
      w=solte[i+gme.npw*gme.Nm[0]];
      g=gsm[int(solte[i])];
      TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(gs[int(solte[i])]*gs[m]+gs[int(solte[i])+npw]*gs[m+npw])/(g*gp);
      if (solte[i]>m)
       eta2vd=eta2d[iju(m,solte[i])];
      else
       eta2vd=eta2d[iju(solte[i],m)];
      tem=epsm2*epsm2*eta2vd*((conj(A2)*W2+conj(B2)*X2)*I2m+(conj(A2)*X2+conj(B2)*W2)*I2p);

      if (solte[i]==m){
       eta1v=eta1d[iju(solte[i],m)];
       eta3v=eta3d[iju(solte[i],m)];
       tem=tem+epsm1*epsm1*eta1v*conj(B1)*(W1*I1p+X1*I1m)+epsm3*epsm3*eta3v*conj(A3)*(W3*I3m+X3*I3p);
      }
      tem=w*w*wk*degegp*tem;

      sradte=sradte+vkd[i+ltm]*tem;
     }

     for (int i=0;i<ltm;i++) { // TM-TE
      w=soltm[i+gme.npw*gme.Nm[1]];
      g=gsm[int(soltm[i])];
      Xx1=Xx(w,g,epsm1);
      Xx3=Xx(w,g,epsm3);
      q=qq(w,g,epsm2);
      TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(-gs[int(soltm[i])]*gs[m+npw]+gs[m]*gs[int(soltm[i])+npw])/(g*gp);
      if (soltm[i]>m)
       eta2vd=eta2d[iju(m,soltm[i])];
      else
       eta2vd=eta2d[iju(soltm[i],m)];
      tem=-complex<double>(0,1)*epsm2*eta2vd*q*((conj(D2)*X2-conj(C2)*W2)*I2m+(conj(D2)*W2-conj(C2)*X2)*I2p);

      if (soltm[i]==m){
       eta1v=eta1d[iju(soltm[i],m)];
       eta3v=eta3d[iju(soltm[i],m)];
       tem=tem-epsm1*eta1v*Xx1*conj(D1)*(W1*I1p+X1*I1m)+epsm3*eta3v*Xx3*conj(C3)*(W3*I3m+X3*I3p);
      }
      tem=wk*degegp*tem;

      sradte=sradte+vkd[i]*tem;
     }

     simw2=simw2+M_PI*(abs(sradte)*abs(sradte)+abs(sradtm)*abs(sradtm))*roj(wk,gp,epsm1,epsm3,j);
    }
   }
  }

  gme.imw=0.5*simw2/(2*M_PI*wk);
 }
 else{ // complex case
  eta2c=(complex<double> *)gme.ceta2;
  vkc=(complex<double> *)gme.cvk;
 
  for (int j=1;j<4;j=j+2){ // computing Im(w^2) 
   for (int m=0;m<npw;m++){
    sradte=complex<double>(0,0);
    sradtm=complex<double>(0,0);

    if (roj(wk,gsm[m],epsm1,epsm3,j)!=0){
     //wp=wk
     gp=gsm[m];
     q1r=qx(wk,gp,epsm1);
     q2r=qx(wk,gp,epsm2);
     q3r=qx(wk,gp,epsm3);
     TMcr(wk,gp,epsm1,epsm2,epsm3,d,Z1,Z2,Z3,Y1,Y2,Y3,j);
     TEcr(wk,gp,epsm1,epsm2,epsm3,d,X1,X2,X3,W1,W2,W3,j);
    
     for (int i=0;i<ltm;i++) { // TM-TM
      w=soltm[i+gme.npw*gme.Nm[1]];
      g=gsm[int(soltm[i])];
      Xx1=Xx(w,g,epsm1);
      Xx3=Xx(w,g,epsm3);
      q=qq(w,g,epsm2);
      TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(gs[int(soltm[i])]*gs[m]+gs[int(soltm[i])+npw]*gs[m+npw])/(g*gp);
      if (soltm[i]>m)
       eta2vc=conj(eta2c[iju(m,soltm[i])]);
      else
       eta2vc=eta2c[iju(soltm[i],m)];
      tem=eta2vc*((conj(C2)*Y2+conj(D2)*Z2)*(g*gp+q*q2r*degegp)*I2m+(conj(C2)*Z2+conj(D2)*Y2)*(g*gp-q*q2r*degegp)*I2p);

      if (soltm[i]==m){
       eta1v=eta1d[iju(soltm[i],m)];
       eta3v=eta3d[iju(soltm[i],m)];
       tem=tem+eta1v*conj(D1)*((g*gp+complex<double>(0,1)*Xx1*q1r*degegp)*Y1*I1p+(g*gp-complex<double>(0,1)*Xx1*q1r*degegp)*Z1*I1m)+eta3v*conj(C3)*((g*gp-complex<double>(0,1)*Xx3*q3r*degegp)*Y3*I3m+(g*gp+complex<double>(0,1)*Xx3*q3r*degegp)*Z3*I3p);
      }

      sradtm=sradtm+conj(vkc[i])*tem;
     }

     for (int i=0;i<lte;i++) { // TE-TM
      w=solte[i+gme.npw*gme.Nm[0]];
      g=gsm[int(solte[i])];
      TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(-gs[m]*gs[int(solte[i])+npw]+gs[int(solte[i])]*gs[m+npw])/(g*gp);
      if (solte[i]>m)
       eta2vc=conj(eta2c[iju(m,solte[i])]);
      else
       eta2vc=eta2c[iju(solte[i],m)];
      tem=epsm2*eta2vc*q2r*((-conj(A2)*Y2+conj(B2)*Z2)*I2m+(conj(A2)*Z2-conj(B2)*Y2)*I2p);

      if (solte[i]==m){
       eta1v=eta1d[iju(solte[i],m)];
       eta3v=eta3d[iju(solte[i],m)];
       tem=tem+epsm1*eta1v*q1r*conj(B1)*(-Y1*I1p+Z1*I1m)+epsm3*eta3v*q3r*conj(A3)*(-Y3*I3m+Z3*I3p);
      }
      tem=complex<double>(0,1)*w*w*degegp*tem;

      sradtm=sradtm+conj(vkc[i+ltm])*tem;
     }

     for (int i=0;i<lte;i++) { // TE-TE
      w=solte[i+gme.npw*gme.Nm[0]];
      g=gsm[int(solte[i])];
      TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(gs[int(solte[i])]*gs[m]+gs[int(solte[i])+npw]*gs[m+npw])/(g*gp);
      if (solte[i]>m)
       eta2vc=conj(eta2c[iju(m,solte[i])]);
      else
       eta2vc=eta2c[iju(solte[i],m)];
      tem=epsm2*epsm2*eta2vc*((conj(A2)*W2+conj(B2)*X2)*I2m+(conj(A2)*X2+conj(B2)*W2)*I2p);

      if (solte[i]==m){
       eta1v=eta1d[iju(solte[i],m)];
       eta3v=eta3d[iju(solte[i],m)];
       tem=tem+epsm1*epsm1*eta1v*conj(B1)*(W1*I1p+X1*I1m)+epsm3*epsm3*eta3v*conj(A3)*(W3*I3m+X3*I3p);
      }
      tem=w*w*wk*degegp*tem;

      sradte=sradte+conj(vkc[i+ltm])*tem;
     }

     for (int i=0;i<ltm;i++) { // TM-TE
      w=soltm[i+gme.npw*gme.Nm[1]];
      g=gsm[int(soltm[i])];
      Xx1=Xx(w,g,epsm1);
      Xx3=Xx(w,g,epsm3);
      q=qq(w,g,epsm2);
      TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
      Izr(w,g,wk,gp,epsm1,epsm2,epsm3,d,I1m,I1p,I2m,I2p,I3m,I3p);
      degegp=(-gs[int(soltm[i])]*gs[m+npw]+gs[m]*gs[int(soltm[i])+npw])/(g*gp);
      if (soltm[i]>m)
       eta2vc=eta2c[iju(m,soltm[i])];
      else
       eta2vc=eta2c[iju(soltm[i],m)];
      tem=-complex<double>(0,1)*epsm2*eta2vc*q*((conj(D2)*X2-conj(C2)*W2)*I2m+(conj(D2)*W2-conj(C2)*X2)*I2p);

      if (soltm[i]==m){
       eta1v=eta1d[iju(soltm[i],m)];
       eta3v=eta3d[iju(soltm[i],m)];
       tem=tem-epsm1*eta1v*Xx1*conj(D1)*(W1*I1p+X1*I1m)+epsm3*eta3v*Xx3*conj(C3)*(W3*I3m+X3*I3p);
      }
      tem=wk*degegp*tem;

      sradte=sradte+conj(vkc[i])*tem;
     }

     simw2=simw2+M_PI*(abs(sradte)*abs(sradte)+abs(sradtm)*abs(sradtm))*roj(wk,gp,epsm1,epsm3,j);
    }
   }
  }

  gme.imw=0.5*simw2/(2*M_PI*wk);
 }
}

/*######################################################################################################*/
  
/*******************************************************************************/
/*************************    Electromagnetic fields   *************************/
/*************************   From photonic dispersion   ************************/
/*******************************************************************************/  

/* Coefficients and wave numbers 
 *
 * MTEc dimensions: (lte,4); ith_row: (B1_i,B2_i,A2_i,A3_i)
 * MTMc dimensions: (ltm,4); ith_row: (D1_i,D2_i,C2_i,C3_i)
 * MXqn dimensions: (lte+ltm,3); ith_row: (Xx1_i,q_i,Xx3_i)
 */

void MTE (double *solte, int lte, int dimsolte, double *gsm, double epsm1, double epsm2, double epsm3, double d, complex<double> *MTEc){ 
 complex<double> B1,B2,A2,A3;
 for (int i=0;i<lte;i++){
  TEc(solte[i+dimsolte],gsm[int(solte[i])],epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
  MTEc[i]=B1;
  MTEc[i+lte]=B2;
  MTEc[i+2*lte]=A2;
  MTEc[i+3*lte]=A3;
 }
}

void MTM (double *soltm, int ltm, int dimsoltm, double *gsm, double epsm1, double epsm2, double epsm3, double d, complex<double> *MTMc){
 complex<double> D1,D2,C2,C3;
 for (int i=0;i<ltm;i++){
  TMc(soltm[i+dimsoltm],gsm[int(soltm[i])],epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
  MTMc[i]=D1;
  MTMc[i+ltm]=D2;
  MTMc[i+2*ltm]=C2;
  MTMc[i+3*ltm]=C3;
 }
}

void MXq (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gsm, double epsm1, double epsm2, double epsm3, double *MXqn){
 int dim=lte+ltm;
 double w,g;
 
 for (int i=0;i<ltm;i++){
  w=soltm[i+dimsoltm];
  g=gsm[int(soltm[i])];
  MXqn[i]=Xx(w,g,epsm1);
  MXqn[i+dim]=qq(w,g,epsm2);
  MXqn[i+2*dim]=Xx(w,g,epsm3);
 }
 for (int i=0;i<lte;i++){
  w=solte[i+dimsolte];
  g=gsm[int(solte[i])];
  MXqn[i+ltm]=Xx(w,g,epsm1);
  MXqn[i+ltm+dim]=qq(w,g,epsm2);
  MXqn[i+ltm+2*dim]=Xx(w,g,epsm3);
 }
}

/*******************************************************************************/
/**********************  Field components in each region   *********************/
/*******************************************************************************/  

// From real diagonalization

complex<double> Hxu_r (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hx upper cladding
 int dim=lte+ltm;
 complex<double> Hxs(0,0);

 for (int i=0;i<ltm;i++){
  Hxs=Hxs-vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*gs[int(soltm[i])+npw]*MTMc[i+3*ltm]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hxs=Hxs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*MTEc[i+3*lte]*MXqn[i+ltm+2*dim]*gs[int(solte[i])]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hxs;
}

complex<double> Hxc_r (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hx core
 int dim=lte+ltm;
 complex<double> Hxs(0,0);

 for (int i=0;i<ltm;i++){
  Hxs=Hxs-vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*gs[int(soltm[i])+npw]*(MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hxs=Hxs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*MXqn[i+ltm+dim]*gs[int(solte[i])]*(-MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*hfnorm);
 }
 return Hxs;
}

complex<double> Hxl_r (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hx lower cladding
 int dim=lte+ltm;
 complex<double> Hxs(0,0); 

 for (int i=0;i<ltm;i++){
  Hxs=Hxs-vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*gs[int(soltm[i])+npw]*MTMc[i]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hxs=Hxs-vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*MTEc[i]*MXqn[i+ltm]*gs[int(solte[i])]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hxs;
}

complex<double> Hyu_r (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hy upper cladding
 int dim=lte+ltm;
 complex<double> Hys(0,0);

 for (int i=0;i<ltm;i++){
  Hys=Hys+vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*gs[int(soltm[i])]*MTMc[i+3*ltm]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hys=Hys+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*MTEc[i+3*lte]*MXqn[i+ltm+2*dim]*gs[int(solte[i])+npw]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hys;
}

complex<double> Hyc_r (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hy core
 int dim=lte+ltm;
 complex<double> Hys(0,0);

 for (int i=0;i<ltm;i++){
  Hys=Hys+vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*gs[int(soltm[i])]*(MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hys=Hys+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*MXqn[i+ltm+dim]*gs[int(solte[i])+npw]*(-MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*hfnorm);
 }
 return Hys;
}

complex<double> Hyl_r (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hy lower cladding
 int dim=lte+ltm;
 complex<double> Hys(0,0); 

 for (int i=0;i<ltm;i++){
  Hys=Hys+vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*gs[int(soltm[i])]*MTMc[i]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hys=Hys-vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*MTEc[i]*MXqn[i+ltm]*gs[int(solte[i])+npw]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hys;
}

complex<double> Hzu_r (double *solte, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, double x, double y, double z){ // Hz upper cladding
 int dim=lte+ltm;
 complex<double> Hzs(0,0);

 for (int i=0;i<lte;i++){
  Hzs=Hzs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*MTEc[i+3*lte]*complex<double>(0,1)*gsm[int(solte[i])]/hfnorm;
 }
 return Hzs;
}

complex<double> Hzc_r (double *solte, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, double x, double y, double z){ // Hz core
 int dim=lte+ltm;
 complex<double> Hzs(0,0);

 for (int i=0;i<lte;i++){
  Hzs=Hzs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*gsm[int(solte[i])]*(MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/hfnorm;
 }
 return Hzs;
}

complex<double> Hzl_r (double *solte, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double hfnorm, double *MXqn, complex<double> *MTEc, double x, double y, double z){ // Hz lower cladding
 int dim=lte+ltm;
 complex<double> Hzs(0,0);

 for (int i=0;i<lte;i++){
  Hzs=Hzs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*MTEc[i]*complex<double>(0,1)*gsm[int(solte[i])]/hfnorm;
 }
 return Hzs;
}

complex<double> Exu_r (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ex upper cladding
 int dim=lte+ltm;
 complex<double> Exs(0,0);

 for (int i=0;i<ltm;i++){
  Exs=Exs+vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*complex<double>(0,1)*MTMc[i+3*ltm]*MXqn[i+2*dim]*gs[int(soltm[i])]/(gsm[int(soltm[i])]*epsm3*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Exs=Exs-vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])+npw]*MTEc[i+3*lte]/(gsm[int(solte[i])]*efnorm);
 }
 return Exs;
}

complex<double> Exc_r (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ex core
 int dim=lte+ltm;
 complex<double> Exs(0,0);

 for (int i=0;i<ltm;i++){
  Exs=Exs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*MXqn[i+dim]*gs[int(soltm[i])]*(-MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*epsm2*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Exs=Exs-vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])+npw]*(MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*efnorm);
 }
 return Exs;
}

complex<double> Exl_r (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ex lower cladding
 int dim=lte+ltm;
 complex<double> Exs(0,0); 

 for (int i=0;i<ltm;i++){
  Exs=Exs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*complex<double>(0,1)*MTMc[i]*MXqn[i]*gs[int(soltm[i])]/(gsm[int(soltm[i])]*epsm1*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Exs=Exs-vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])+npw]*MTEc[i]/(gsm[int(solte[i])]*efnorm);
 }
 return Exs;
}

complex<double> Eyu_r (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ey upper cladding
 int dim=lte+ltm;
 complex<double> Eys(0,0);

 for (int i=0;i<ltm;i++){
  Eys=Eys+vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*complex<double>(0,1)*MTMc[i+3*ltm]*MXqn[i+2*dim]*gs[int(soltm[i])+npw]/(gsm[int(soltm[i])]*epsm3*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Eys=Eys+vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])]*MTEc[i+3*lte]/(gsm[int(solte[i])]*efnorm);
 }
 return Eys;
}

complex<double> Eyc_r (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ey core
 int dim=lte+ltm;
 complex<double> Eys(0,0);

 for (int i=0;i<ltm;i++){
  Eys=Eys-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*MXqn[i+dim]*gs[int(soltm[i])+npw]*(-MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*epsm2*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Eys=Eys+vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])]*(MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*efnorm);
 }
 return Eys;
}

complex<double> Eyl_r (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ey lower cladding
 int dim=lte+ltm;
 complex<double> Eys(0,0); 

 for (int i=0;i<ltm;i++){
  Eys=Eys-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*complex<double>(0,1)*MTMc[i]*MXqn[i]*gs[int(soltm[i])+npw]/(gsm[int(soltm[i])]*epsm1*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Eys=Eys+vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])]*MTEc[i]/(gsm[int(solte[i])]*efnorm);
 }
 return Eys;
}

complex<double> Ezu_r (double *soltm, int lte, int ltm, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTMc, double x, double y, double z){ // Ez upper cladding
 int dim=lte+ltm;
 complex<double> Ezs(0,0);

 for (int i=0;i<ltm;i++){
  Ezs=Ezs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*MTMc[i+3*ltm]*gsm[int(soltm[i])]/(epsm3*soltm[i+dimsoltm]*efnorm);
 }
 return Ezs;
}

complex<double> Ezc_r (double *soltm, int lte, int ltm, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTMc, double x, double y, double z){ // Ez core
 int dim=lte+ltm;
 complex<double> Ezs(0,0);

 for (int i=0;i<ltm;i++){
  Ezs=Ezs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*gsm[int(soltm[i])]*(MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(epsm2*soltm[i+dimsoltm]*efnorm);
 }
 return Ezs;
}

complex<double> Ezl_r (double *soltm, int lte, int ltm, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, double *vk, double efnorm, double *MXqn, complex<double> *MTMc, double x, double y, double z){ // Ez lower cladding
 int dim=lte+ltm;
 complex<double> Ezs(0,0);

 for (int i=0;i<ltm;i++){
  Ezs=Ezs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*MTMc[i]*gsm[int(soltm[i])]/(epsm1*soltm[i+dimsoltm]*efnorm);
 }
 return Ezs;
}

// From complex diagonalization

complex<double> Hxu_c (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hx upper cladding
 int dim=lte+ltm;
 complex<double> Hxs(0,0);

 for (int i=0;i<ltm;i++){
  Hxs=Hxs-vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*gs[int(soltm[i])+npw]*MTMc[i+3*ltm]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hxs=Hxs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*MTEc[i+3*lte]*MXqn[i+ltm+2*dim]*gs[int(solte[i])]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hxs;
}

complex<double> Hxc_c (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hx core
 int dim=lte+ltm;
 complex<double> Hxs(0,0);

 for (int i=0;i<ltm;i++){
  Hxs=Hxs-vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*gs[int(soltm[i])+npw]*(MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hxs=Hxs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*MXqn[i+ltm+dim]*gs[int(solte[i])]*(-MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*hfnorm);
 }
 return Hxs;
}

complex<double> Hxl_c (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hx lower cladding
 int dim=lte+ltm;
 complex<double> Hxs(0,0); 

 for (int i=0;i<ltm;i++){
  Hxs=Hxs-vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*gs[int(soltm[i])+npw]*MTMc[i]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hxs=Hxs-vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*MTEc[i]*MXqn[i+ltm]*gs[int(solte[i])]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hxs;
}

complex<double> Hyu_c (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hy upper cladding
 int dim=lte+ltm;
 complex<double> Hys(0,0);

 for (int i=0;i<ltm;i++){
  Hys=Hys+vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*gs[int(soltm[i])]*MTMc[i+3*ltm]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hys=Hys+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*MTEc[i+3*lte]*MXqn[i+ltm+2*dim]*gs[int(solte[i])+npw]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hys;
}

complex<double> Hyc_c (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hy core
 int dim=lte+ltm;
 complex<double> Hys(0,0);

 for (int i=0;i<ltm;i++){
  Hys=Hys+vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*gs[int(soltm[i])]*(MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hys=Hys+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*MXqn[i+ltm+dim]*gs[int(solte[i])+npw]*(-MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*hfnorm);
 }
 return Hys;
}

complex<double> Hyl_c (double *solte, double *soltm, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Hy lower cladding
 int dim=lte+ltm;
 complex<double> Hys(0,0); 

 for (int i=0;i<ltm;i++){
  Hys=Hys+vk[i]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*gs[int(soltm[i])]*MTMc[i]/(gsm[int(soltm[i])]*hfnorm);
 }
 for (int i=0;i<lte;i++){
  Hys=Hys-vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*MTEc[i]*MXqn[i+ltm]*gs[int(solte[i])+npw]/(gsm[int(solte[i])]*hfnorm);
 }
 return Hys;
}

complex<double> Hzu_c (double *solte, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, double x, double y, double z){ // Hz upper cladding
 int dim=lte+ltm;
 complex<double> Hzs(0,0);

 for (int i=0;i<lte;i++){
  Hzs=Hzs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*MTEc[i+3*lte]*complex<double>(0,1)*gsm[int(solte[i])]/hfnorm;
 }
 return Hzs;
}

complex<double> Hzc_c (double *solte, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, double x, double y, double z){ // Hz core
 int dim=lte+ltm;
 complex<double> Hzs(0,0);

 for (int i=0;i<lte;i++){
  Hzs=Hzs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*gsm[int(solte[i])]*(MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/hfnorm;
 }
 return Hzs;
}

complex<double> Hzl_c (double *solte, int lte, int ltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double hfnorm, double *MXqn, complex<double> *MTEc, double x, double y, double z){ // Hz lower cladding
 int dim=lte+ltm;
 complex<double> Hzs(0,0);

 for (int i=0;i<lte;i++){
  Hzs=Hzs+vk[i+ltm]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*MTEc[i]*complex<double>(0,1)*gsm[int(solte[i])]/hfnorm;
 }
 return Hzs;
}

complex<double> Exu_c (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ex upper cladding
 int dim=lte+ltm;
 complex<double> Exs(0,0);

 for (int i=0;i<ltm;i++){
  Exs=Exs+vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*complex<double>(0,1)*MTMc[i+3*ltm]*MXqn[i+2*dim]*gs[int(soltm[i])]/(gsm[int(soltm[i])]*epsm3*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Exs=Exs-vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])+npw]*MTEc[i+3*lte]/(gsm[int(solte[i])]*efnorm);
 }
 return Exs;
}

complex<double> Exc_c (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ex core
 int dim=lte+ltm;
 complex<double> Exs(0,0);

 for (int i=0;i<ltm;i++){
  Exs=Exs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*MXqn[i+dim]*gs[int(soltm[i])]*(-MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*epsm2*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Exs=Exs-vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])+npw]*(MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*efnorm);
 }
 return Exs;
}

complex<double> Exl_c (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ex lower cladding
 int dim=lte+ltm;
 complex<double> Exs(0,0); 

 for (int i=0;i<ltm;i++){
  Exs=Exs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*complex<double>(0,1)*MTMc[i]*MXqn[i]*gs[int(soltm[i])]/(gsm[int(soltm[i])]*epsm1*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Exs=Exs-vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])+npw]*MTEc[i]/(gsm[int(solte[i])]*efnorm);
 }
 return Exs;
}

complex<double> Eyu_c (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ey upper cladding
 int dim=lte+ltm;
 complex<double> Eys(0,0);

 for (int i=0;i<ltm;i++){
  Eys=Eys+vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*complex<double>(0,1)*MTMc[i+3*ltm]*MXqn[i+2*dim]*gs[int(soltm[i])+npw]/(gsm[int(soltm[i])]*epsm3*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Eys=Eys+vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)-MXqn[i+ltm+2*dim]*(z-0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])]*MTEc[i+3*lte]/(gsm[int(solte[i])]*efnorm);
 }
 return Eys;
}

complex<double> Eyc_c (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ey core
 int dim=lte+ltm;
 complex<double> Eys(0,0);

 for (int i=0;i<ltm;i++){
  Eys=Eys-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*MXqn[i+dim]*gs[int(soltm[i])+npw]*(-MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(gsm[int(soltm[i])]*epsm2*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Eys=Eys+vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])]*(MTEc[i+2*lte]*exp(complex<double>(0,1)*MXqn[i+ltm+dim]*z)+MTEc[i+lte]*exp(complex<double>(0,-1)*MXqn[i+ltm+dim]*z))/(gsm[int(solte[i])]*efnorm);
 }
 return Eys;
}

complex<double> Eyl_c (double *solte, double *soltm, int lte, int ltm, int dimsolte, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTEc, complex<double> *MTMc, double x, double y, double z){ // Ey lower cladding
 int dim=lte+ltm;
 complex<double> Eys(0,0); 

 for (int i=0;i<ltm;i++){
  Eys=Eys-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*complex<double>(0,1)*MTMc[i]*MXqn[i]*gs[int(soltm[i])+npw]/(gsm[int(soltm[i])]*epsm1*soltm[i+dimsoltm]*efnorm);
 }
 for (int i=0;i<lte;i++){
  Eys=Eys+vk[i+ltm]*solte[i+dimsolte]*exp(complex<double>(0,1)*(gs[int(solte[i])]*x+gs[int(solte[i])+npw]*y)+MXqn[i+ltm]*(z+0.5*d))*complex<double>(0,1)*solte[i+dimsolte]*gs[int(solte[i])]*MTEc[i]/(gsm[int(solte[i])]*efnorm);
 }
 return Eys;
}

complex<double> Ezu_c (double *soltm, int lte, int ltm, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTMc, double x, double y, double z){ // Ez upper cladding
 int dim=lte+ltm;
 complex<double> Ezs(0,0);

 for (int i=0;i<ltm;i++){
  Ezs=Ezs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)-MXqn[i+2*dim]*(z-0.5*d))*MTMc[i+3*ltm]*gsm[int(soltm[i])]/(epsm3*soltm[i+dimsoltm]*efnorm);
 }
 return Ezs;
}

complex<double> Ezc_c (double *soltm, int lte, int ltm, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTMc, double x, double y, double z){ // Ez core
 int dim=lte+ltm;
 complex<double> Ezs(0,0);

 for (int i=0;i<ltm;i++){
  Ezs=Ezs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y))*gsm[int(soltm[i])]*(MTMc[i+2*ltm]*exp(complex<double>(0,1)*MXqn[i+dim]*z)+MTMc[i+ltm]*exp(complex<double>(0,-1)*MXqn[i+dim]*z))/(epsm2*soltm[i+dimsoltm]*efnorm);
 }
 return Ezs;
}

complex<double> Ezl_c (double *soltm, int lte, int ltm, int dimsoltm, double *gs, double *gsm, int npw, double epsm1, double epsm2, double epsm3, double d, complex<double> *vk, double efnorm, double *MXqn, complex<double> *MTMc, double x, double y, double z){ // Ez lower cladding
 int dim=lte+ltm;
 complex<double> Ezs(0,0);

 for (int i=0;i<ltm;i++){
  Ezs=Ezs-vk[i]*soltm[i+dimsoltm]*exp(complex<double>(0,1)*(gs[int(soltm[i])]*x+gs[int(soltm[i])+npw]*y)+MXqn[i]*(z+0.5*d))*MTMc[i]*gsm[int(soltm[i])]/(epsm1*soltm[i+dimsoltm]*efnorm);
 }
 return Ezs;
}

/*######################################################################################################*/

/*******************************************************************************/
/*********************    Electric field normalization   ***********************/
/*******************************************************************************/ 

void enorm (st_gme &gme, st_eps &eps){
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *eta1=(double *)gme.eta1,*eta3=(double *)gme.eta3,*solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d;
 int dim=ltm+lte;
 double w,wp,g,gp,q,qp,degegp,I1,I2m,I2p,I3,eta1v,eta2vd,eta3v,factor,fd,*eta2d,*vkd;
 complex<double> A2,A3,B1,B2,C2,C3,D1,D2,A2p,A3p,B1p,B2p,C2p,C3p,D1p,D2p,eta2vc,*eta2c,*vkc,efnorm(0,0);

 if (gme.psym==1){ // real case
  eta2d=(double *)gme.eta2;
  vkd=(double *)gme.vk;

  for (int i=0;i<ltm;i++){  // TM-TM integral  
   for (int j=i;j<ltm;j++){ 
    if (i==j)
     fd=0.5;
    else
     fd=1.0;
  
    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    q=qq(w,g,epsm2);
    wp=soltm[j+gme.npw*gme.Nm[1]];
    gp=gsm[int(soltm[j])];
    qp=qq(wp,gp,epsm2);  
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
    degegp=(gs[int(soltm[i])]*gs[int(soltm[j])]+gs[int(soltm[i])+npw]*gs[int(soltm[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>soltm[j])
     eta2vd=eta2d[iju(soltm[j],soltm[i])];
    else
     eta2vd=eta2d[iju(soltm[i],soltm[j])];
    efnorm=efnorm+vkd[i]*vkd[j]*(eta2vd*((conj(C2)*C2p+conj(D2)*D2p)*(q*qp*degegp+g*gp)*I2m+(conj(C2)*D2p+conj(D2)*C2p)*(-q*qp*degegp+g*gp)*I2p))*fd;

    if (soltm[i]==soltm[j]){
     eta1v=eta1[iju(soltm[i],soltm[j])];
     eta3v=eta3[iju(soltm[i],soltm[j])];
     efnorm=efnorm+vkd[i]*vkd[j]*(eta1v*conj(D1)*D1p*(Xx(w,g,epsm1)*Xx(wp,gp,epsm1)*degegp+g*gp)*I1+eta3v*conj(C3)*C3p*(Xx(w,g,epsm3)*Xx(wp,gp,epsm3)*degegp+g*gp)*I3)*fd;
    }
   }
  }

  for (int i=0;i<lte;i++){  // TE-TE integral 
   for (int j=i;j<lte;j++){ 
    if (i==j)
     fd=0.5;
    else
     fd=1.0;

    w=solte[i+gme.npw*gme.Nm[0]];
    g=gsm[int(solte[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(gs[int(solte[i])]*gs[int(solte[j])]+gs[int(solte[i])+npw]*gs[int(solte[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (solte[i]>solte[j])
     eta2vd=eta2d[iju(solte[j],solte[i])];
    else
     eta2vd=eta2d[iju(solte[i],solte[j])];
    factor=w*w*wp*wp*degegp;
    efnorm=efnorm+vkd[i+ltm]*vkd[j+ltm]*(factor*epsm2*epsm2*eta2vd*((conj(A2)*A2p+conj(B2)*B2p)*I2m+(conj(A2)*B2p+conj(B2)*A2p)*I2p))*fd; 

    if (solte[i]==solte[j]){
     eta1v=eta1[iju(solte[i],solte[j])];
     eta3v=eta3[iju(solte[i],solte[j])];
     efnorm=efnorm+vkd[i+ltm]*vkd[j+ltm]*(factor*(epsm1*epsm1*eta1v*conj(B1)*B1p*I1+epsm3*epsm3*eta3v*conj(A3)*A3p*I3))*fd;
    }
   }  
  }

  for (int i=0;i<ltm;i++){  // TM-TE integral  
   for (int j=0;j<lte;j++){
     fd=1.0;
   
    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(-gs[int(soltm[i])]*gs[int(solte[j])+npw]+gs[int(solte[j])]*gs[int(soltm[i])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>solte[j])
     eta2vd=eta2d[iju(solte[j],soltm[i])];
    else
     eta2vd=eta2d[iju(soltm[i],solte[j])];
    factor=wp*wp*degegp;
    efnorm=efnorm+vkd[i]*vkd[j+ltm]*factor*complex<double>(0,-epsm2)*eta2vd*qq(w,g,epsm2)*((-conj(C2)*A2p+conj(D2)*B2p)*I2m+(conj(D2)*A2p-conj(C2)*B2p)*I2p)*fd;

    if (soltm[i]==solte[j]){
     eta1v=eta1[iju(soltm[i],solte[j])];
     eta3v=eta3[iju(soltm[i],solte[j])];
     efnorm=efnorm+vkd[i]*vkd[j+ltm]*factor*(-epsm1*eta1v*conj(D1)*B1p*Xx(w,g,epsm1)*I1+epsm3*eta3v*conj(C3)*A3p*Xx(w,g,epsm3)*I3)*fd;
    }
   }
  }

 gme.enorm=sqrt(2*Acell(eps)*real(efnorm));
 }
 else{ // complex case
  eta2c=(complex<double> *)gme.ceta2;
  vkc=(complex<double> *)gme.cvk;

  for (int i=0;i<ltm;i++){  // TM-TM integral  
   for (int j=i;j<ltm;j++){ 
    if (i==j)
     fd=0.5;
    else
     fd=1.0;
  
    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    q=qq(w,g,epsm2);
    wp=soltm[j+gme.npw*gme.Nm[1]];
    gp=gsm[int(soltm[j])];
    qp=qq(wp,gp,epsm2);  
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TMc(wp,gp,epsm1,epsm2,epsm3,d,D1p,D2p,C2p,C3p);
    degegp=(gs[int(soltm[i])]*gs[int(soltm[j])]+gs[int(soltm[i])+npw]*gs[int(soltm[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>soltm[j])
     eta2vc=conj(eta2c[iju(soltm[j],soltm[i])]);
    else
     eta2vc=eta2c[iju(soltm[i],soltm[j])];
    efnorm=efnorm+conj(vkc[i])*vkc[j]*(eta2vc*((conj(C2)*C2p+conj(D2)*D2p)*(q*qp*degegp+g*gp)*I2m+(conj(C2)*D2p+conj(D2)*C2p)*(-q*qp*degegp+g*gp)*I2p))*fd;

    if (soltm[i]==soltm[j]){
     eta1v=eta1[iju(soltm[i],soltm[j])];
     eta3v=eta3[iju(soltm[i],soltm[j])];
     efnorm=efnorm+conj(vkc[i])*vkc[j]*(eta1v*conj(D1)*D1p*(Xx(w,g,epsm1)*Xx(wp,gp,epsm1)*degegp+g*gp)*I1+eta3v*conj(C3)*C3p*(Xx(w,g,epsm3)*Xx(wp,gp,epsm3)*degegp+g*gp)*I3)*fd;
    }
   }
  }

  for (int i=0;i<lte;i++){  // TE-TE integral 
   for (int j=i;j<lte;j++){ 
    if (i==j)
     fd=0.5;
    else
     fd=1.0;

    w=solte[i+gme.npw*gme.Nm[0]];
    g=gsm[int(solte[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TEc(w,g,epsm1,epsm2,epsm3,d,B1,B2,A2,A3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(gs[int(solte[i])]*gs[int(solte[j])]+gs[int(solte[i])+npw]*gs[int(solte[j])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (solte[i]>solte[j])
     eta2vc=conj(eta2c[iju(solte[j],solte[i])]);
    else
     eta2vc=eta2c[iju(solte[i],solte[j])];
    factor=w*w*wp*wp*degegp;
    efnorm=efnorm+conj(vkc[i+ltm])*vkc[j+ltm]*(factor*epsm2*epsm2*eta2vc*((conj(A2)*A2p+conj(B2)*B2p)*I2m+(conj(A2)*B2p+conj(B2)*A2p)*I2p))*fd; 

    if (solte[i]==solte[j]){
     eta1v=eta1[iju(solte[i],solte[j])];
     eta3v=eta3[iju(solte[i],solte[j])];
     efnorm=efnorm+conj(vkc[i+ltm])*vkc[j+ltm]*(factor*(epsm1*epsm1*eta1v*conj(B1)*B1p*I1+epsm3*epsm3*eta3v*conj(A3)*A3p*I3))*fd;
    }
   }  
  }

  fd=1.0;
  for (int i=0;i<ltm;i++){  // TM-TE integral  
   for (int j=0;j<lte;j++){

    w=soltm[i+gme.npw*gme.Nm[1]];
    g=gsm[int(soltm[i])];
    wp=solte[j+gme.npw*gme.Nm[0]];
    gp=gsm[int(solte[j])];
    TMc(w,g,epsm1,epsm2,epsm3,d,D1,D2,C2,C3);
    TEc(wp,gp,epsm1,epsm2,epsm3,d,B1p,B2p,A2p,A3p);
    degegp=(-gs[int(soltm[i])]*gs[int(solte[j])+npw]+gs[int(solte[j])]*gs[int(soltm[i])+npw])/(g*gp);
    Iz(w,g,wp,gp,epsm1,epsm2,epsm3,d,I1,I2m,I2p,I3);
    if (soltm[i]>solte[j])
     eta2vc=conj(eta2c[iju(solte[j],soltm[i])]);
    else
     eta2vc=eta2c[iju(soltm[i],solte[j])];
    factor=wp*wp*degegp;
    efnorm=efnorm+conj(vkc[i])*vkc[j+ltm]*factor*complex<double>(0,-epsm2)*eta2vc*qq(w,g,epsm2)*((-conj(C2)*A2p+conj(D2)*B2p)*I2m+(conj(D2)*A2p-conj(C2)*B2p)*I2p)*fd;

    if (soltm[i]==solte[j]){
     eta1v=eta1[iju(soltm[i],solte[j])];
     eta3v=eta3[iju(soltm[i],solte[j])];
     efnorm=efnorm+conj(vkc[i])*vkc[j+ltm]*factor*(-epsm1*eta1v*conj(D1)*B1p*Xx(w,g,epsm1)*I1+epsm3*eta3v*conj(C3)*A3p*Xx(w,g,epsm3)*I3)*fd;
    }
   }
  }

 gme.enorm=sqrt(2*Acell(eps)*real(efnorm));
 }
}

/*######################################################################################################*/

/*******************************************************************************/
/**********************  Field components in all regions   *********************/
/*******************************************************************************/ 

complex<double> Hx (st_gme &gme, double x, double y, double z){ // Hx
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,hfnorm=gme.hnorm;
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  if (z<-0.5*d)
   return Hxl_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return Hxc_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return Hxu_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,x,y,z);
}
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  if (z<-0.5*d)
   return Hxl_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return Hxc_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return Hxu_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,x,y,z);
}

 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

complex<double> Hy (st_gme &gme, double x, double y, double z){ // Hy
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,hfnorm=gme.hnorm;
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  if (z<-0.5*d)
   return Hyl_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return Hyc_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return Hyu_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,x,y,z);
}
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  if (z<-0.5*d)
   return Hyl_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return Hyc_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return Hyu_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,x,y,z);
}

 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

complex<double> Hz (st_gme &gme, double x, double y, double z){ // Hz
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,hfnorm=gme.hnorm;
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  if (z<-0.5*d)
   return Hzl_r(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return Hzc_r(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,x,y,z);
  else
   return Hzu_r(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,x,y,z);
}
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  if (z<-0.5*d)
   return Hzl_c(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return Hzc_c(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,x,y,z);
  else
   return Hzu_c(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,x,y,z);
}

 delete[] MTEc;
 delete[] MXqn;
}

complex<double> Dx (st_gme &gme, double x, double y, double z){ // Dx
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,efnorm=gme.enorm;
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  if (z<-0.5*d)
   return epsm1*Exl_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return epsm2*Exc_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return epsm3*Exu_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,x,y,z);
}
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  if (z<-0.5*d)
   return epsm1*Exl_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return epsm2*Exc_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return epsm3*Exu_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,x,y,z);
}

 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

complex<double> Dy (st_gme &gme, double x, double y, double z){ // Dy
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,efnorm=gme.enorm;
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  if (z<-0.5*d)
   return epsm1*Eyl_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return epsm2*Eyc_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return epsm3*Eyu_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,x,y,z);
}
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  if (z<-0.5*d)
   return epsm1*Eyl_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return epsm2*Eyc_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,x,y,z);
  else
   return epsm3*Eyu_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,x,y,z);
}

 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

complex<double> Dz (st_gme &gme, double x, double y, double z){ // Dz
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,efnorm=gme.enorm;
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,*vkc;

 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  if (z<-0.5*d)
   return epsm1*Ezl_r(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return epsm2*Ezc_r(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTMc,x,y,z);
  else
   return epsm3*Ezu_r(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTMc,x,y,z);
}
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  if (z<-0.5*d)
   return epsm1*Ezl_c(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTMc,x,y,z);
  else if (fabs(z)<=0.5*d)
   return epsm2*Ezc_c(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTMc,x,y,z);
  else
   return epsm3*Ezu_c(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTMc,x,y,z);
}

 delete[] MTMc;
 delete[] MXqn;
}

/*######################################################################################################*/

/*******************************************************************************/
/******************  Writing Electromagnetic field components  *****************/
/*******************************************************************************/ 

struct field_box_lim {  
 double xi; // initial x value
 double xf; // final x value
 double dx; // step in x
 double yi; // initial y value
 double yf; // final y value
 double dy; // step in y
 double zi; // initial z value
 double zf; // final z value
 double dz; // step in z
 double stat_i; // initial state to compute (begins from zero)
 double stat_f; // final state to compute
 double kstat_i; // initial k state to compute
 double kstat_f; // final k state to compute
 int cfds; // if equal to 1 computes the fields
 int Dx; // if equal to 1 writes the Dx component
 int Dy; // if equal to 1 writes the Dy component
 int Dz; // if equal to 1 writes the Dz component
 int Hx; // if equal to 1 writes the Hx component
 int Hy; // if equal to 1 writes the Hy component
 int Hz; // if equal to 1 writes the Hz component
}flim;

void writeHx (st_gme &gme, field_box_lim &flim, int kbz, int stat){ // writes the Hx component
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,hfnorm=gme.hnorm;

 std::stringstream s_s;
 s_s<<"Hx-k_"<<kbz<<"-band_"<<stat<<".dat";
 std::string s_e = s_s.str();
 const char *c_e=s_e.c_str();

 ofstream Hxout (c_e);
 Hxout.precision(set_pre);
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,Hxc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hxc=Hxl_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hxc)<<"	"<<imag(Hxc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hxc=Hxc_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hxc)<<"	"<<imag(Hxc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hxc=Hxu_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hxc)<<"	"<<imag(Hxc)<<"\n";
     }
    }
   }
  }
 }
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hxc=Hxl_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hxc)<<"	"<<imag(Hxc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hxc=Hxc_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hxc)<<"	"<<imag(Hxc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hxc=Hxu_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hxc)<<"	"<<imag(Hxc)<<"\n";
     }
    }
   }
  }
 }
 Hxout.close();
 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

void writeHy (st_gme &gme, field_box_lim &flim, int kbz, int stat){ // writes the Hy component
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,hfnorm=gme.hnorm;

 std::stringstream s_s;
 s_s<<"Hy-k_"<<kbz<<"-band_"<<stat<<".dat";
 std::string s_e = s_s.str();
 const char *c_e=s_e.c_str();

 ofstream Hyout (c_e);
 Hyout.precision(set_pre);
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,Hyc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hyc=Hyl_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hyc)<<"	"<<imag(Hyc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hyc=Hyc_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hyc)<<"	"<<imag(Hyc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hyc=Hyu_r(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hyc)<<"	"<<imag(Hyc)<<"\n";
     }
    }
   }
  }
 }
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hyc=Hyl_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hyc)<<"	"<<imag(Hyc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hyc=Hyc_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hyc)<<"	"<<imag(Hyc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hyc=Hyu_c(solte,soltm,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Hyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hyc)<<"	"<<imag(Hyc)<<"\n";
     }
    }
   }
  }
 }
 Hyout.close();
 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

void writeHz (st_gme &gme, field_box_lim &flim, int kbz, int stat){ // writes the Hz component

 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,hfnorm=gme.hnorm;

 std::stringstream s_s;
 s_s<<"Hz-k_"<<kbz<<"-band_"<<stat<<".dat";
 std::string s_e = s_s.str();
 const char *c_e=s_e.c_str();

 ofstream Hzout (c_e);
 Hzout.precision(set_pre);
 double *MXqn,*vkd;
 complex<double> *MTEc,Hzc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hzc=Hzl_r(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,ix,iy,iz);
      Hzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hzc)<<"	"<<imag(Hzc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hzc=Hzc_r(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,ix,iy,iz);
      Hzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hzc)<<"	"<<imag(Hzc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hzc=Hzu_r(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,hfnorm,MXqn,MTEc,ix,iy,iz);
      Hzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hzc)<<"	"<<imag(Hzc)<<"\n";
     }
    }
   }
  }
 }
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hzc=Hzl_c(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,ix,iy,iz);
      Hzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hzc)<<"	"<<imag(Hzc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hzc=Hzc_c(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,ix,iy,iz);
      Hzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hzc)<<"	"<<imag(Hzc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Hzc=Hzu_c(solte,lte,ltm,gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,hfnorm,MXqn,MTEc,ix,iy,iz);
      Hzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Hzc)<<"	"<<imag(Hzc)<<"\n";
     }
    }
   }
  }
 }
 Hzout.close();
 delete[] MTEc;
 delete[] MXqn;
}

void writeDx (st_gme &gme, field_box_lim &flim, int kbz, int stat){ // writes the Dx component
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,efnorm=gme.enorm;

 std::stringstream s_s;
 s_s<<"Dx-k_"<<kbz<<"-band_"<<stat<<".dat";
 std::string s_e = s_s.str();
 const char *c_e=s_e.c_str();

 ofstream Dxout (c_e);
 Dxout.precision(set_pre);
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,Dxc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dxc=epsm1*Exl_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dxc)<<"	"<<imag(Dxc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dxc=epsm2*Exc_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dxc)<<"	"<<imag(Dxc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dxc=epsm3*Exu_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dxc)<<"	"<<imag(Dxc)<<"\n";
     }
    }
   }
  }
 }
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dxc=epsm1*Exl_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dxc)<<"	"<<imag(Dxc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dxc=epsm2*Exc_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dxc)<<"	"<<imag(Dxc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dxc=epsm3*Exu_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dxout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dxc)<<"	"<<imag(Dxc)<<"\n";
     }
    }
   }
  }
 }
 Dxout.close();
 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

void writeDy (st_gme &gme, field_box_lim &flim, int kbz, int stat){ // writes the Dy component

 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,efnorm=gme.enorm;

 std::stringstream s_s;
 s_s<<"Dy-k_"<<kbz<<"-band_"<<stat<<".dat";
 std::string s_e = s_s.str();
 const char *c_e=s_e.c_str();

 ofstream Dyout (c_e);
 Dyout.precision(set_pre);
 double *MXqn,*vkd;
 complex<double> *MTEc,*MTMc,Dyc,*vkc;

 MTEc=new complex<double>[lte*4]();
 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTE(solte,lte,gme.npw*gme.Nm[0],gsm,epsm1,epsm2,epsm3,d,MTEc);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dyc=epsm1*Eyl_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dyc)<<"	"<<imag(Dyc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dyc=epsm2*Eyc_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dyc)<<"	"<<imag(Dyc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dyc=epsm3*Eyu_r(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dyc)<<"	"<<imag(Dyc)<<"\n";
     }
    }
   }
  }
 }
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dyc=epsm1*Eyl_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dyc)<<"	"<<imag(Dyc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dyc=epsm2*Eyc_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dyc)<<"	"<<imag(Dyc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dyc=epsm3*Eyu_c(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTEc,MTMc,ix,iy,iz);
      Dyout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dyc)<<"	"<<imag(Dyc)<<"\n";
     }
    }
   }
  }
 }
 Dyout.close();
 delete[] MTEc;
 delete[] MTMc;
 delete[] MXqn;
}

void writeDz (st_gme &gme, field_box_lim &flim, int kbz, int stat){ // writes the Dz component
 int lte=gme.lte,ltm=gme.ltm,npw=gme.npw;
 double *solte=(double *)gme.solte,*soltm=(double *)gme.soltm,*gs=(double *)gme.gs,*gsm=(double *)gme.gsm,epsm1=gme.epsm1,epsm2=gme.epsm2,epsm3=gme.epsm3,d=gme.d,efnorm=gme.enorm;

 std::stringstream s_s;
 s_s<<"Dz-k_"<<kbz<<"-band_"<<stat<<".dat";
 std::string s_e = s_s.str();
 const char *c_e=s_e.c_str();

 ofstream Dzout (c_e);
 Dzout.precision(set_pre);
 double *MXqn,*vkd;
 complex<double> *MTMc,Dzc,*vkc;

 MTMc=new complex<double>[ltm*4]();
 MXqn=new double[(lte+ltm)*3]();

 MXq(solte,soltm,lte,ltm,gme.npw*gme.Nm[0],gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,MXqn);
 MTM(soltm,ltm,gme.npw*gme.Nm[1],gsm,epsm1,epsm2,epsm3,d,MTMc);

 if (gme.psym==1){ // real case
  vkd=(double *)gme.vk; 

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dzc=epsm1*Ezl_r(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTMc,ix,iy,iz);
      Dzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dzc)<<"	"<<imag(Dzc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dzc=epsm2*Ezc_r(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTMc,ix,iy,iz);
      Dzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dzc)<<"	"<<imag(Dzc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dzc=epsm3*Ezu_r(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkd,efnorm,MXqn,MTMc,ix,iy,iz);
      Dzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dzc)<<"	"<<imag(Dzc)<<"\n";
     }
    }
   }
  }
 }
 else{ // complex case
  vkc=(complex<double> *)gme.cvk;

  for (double iz=flim.zi;iz<=flim.zf;iz=iz+flim.dz){
   if (iz<-0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dzc=epsm1*Ezl_c(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTMc,ix,iy,iz);
      Dzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dzc)<<"	"<<imag(Dzc)<<"\n";
     }
    }
   }
   else if (fabs(iz)<=0.5*d){
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dzc=epsm2*Ezc_c(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTMc,ix,iy,iz);
      Dzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dzc)<<"	"<<imag(Dzc)<<"\n";
     }
    }
   } 
   else{
    for (double iy=flim.yi;iy<=flim.yf;iy=iy+flim.dy){
     for (double ix=flim.xi;ix<=flim.xf;ix=ix+flim.dx){
      Dzc=epsm3*Ezu_c(soltm,lte,ltm,gme.npw*gme.Nm[1],gs,gsm,npw,epsm1,epsm2,epsm3,d,vkc,efnorm,MXqn,MTMc,ix,iy,iz);
      Dzout<<ix<<"	"<<iy<<"	"<<iz<<"	"<<real(Dzc)<<"	"<<imag(Dzc)<<"\n";
     }
    }
   }
  }
 }
 Dzout.close();
 delete[] MTMc;
 delete[] MXqn;
}

// Calling the writing field functions

void writefields (st_gme &gme, field_box_lim &flim, st_eps &eps, int kbz){
 if (gme.k[4]!=0 && flim.cfds!=0){
  cout<<"Writing fields:\n";
  cout.flush();
 
  for (int stat=flim.stat_i;stat<(flim.stat_f+1);stat++){
   cout<<"state "<<stat;
   cout.flush();
   t=clock();

   if (gme.psym==1){ // real case
    dgetvj(gme.eigenvec,gme.vk,stat-gme.stat_i,gme.lte+gme.ltm);
   }
   else{ // complex case
    zgetvj(gme.ceigenvec,gme.cvk,stat-gme.stat_i,gme.lte+gme.ltm);
   } 

   gme.hnorm=sqrt(Acell(eps)); 

   if (flim.Dx==1||flim.Dy==1||flim.Dz==1){
    t=clock();
    cout<<" [normalizing electric field... ";
    cout.flush();
    enorm(gme,eps);
    t=(clock()-t)/CLOCKS_PER_SEC;
    cout<<t<<" seconds] ";
   }

   t=clock();
   if (flim.Dx==1)
    writeDx(gme,flim,kbz,stat);
   if (flim.Dy==1)
    writeDy(gme,flim,kbz,stat);
   if (flim.Dz==1)
    writeDz(gme,flim,kbz,stat);
   if (flim.Hx==1)
    writeHx(gme,flim,kbz,stat);
   if (flim.Hy==1)
    writeHy(gme,flim,kbz,stat);
   if (flim.Hz==1)
    writeHz(gme,flim,kbz,stat);
       
   t=(clock()-t)/CLOCKS_PER_SEC;
   cout<<"  written in "<<t<<" seconds\n";
   cout.flush();
  }
 }
}


/*######################################################################################################*/

/*******************************************************************************/
/**************    Diagonalization of the hamiltonian matrix    ****************/
/*******************************************************************************/

void eigensolver (st_gme &gme, field_box_lim &flim){
 gme.eigenw=new double[gme.lte+gme.ltm]();
 if (gme.psym==1){ // real case
  if (gme.cimw!=1 && flim.cfds!=1){
   cout<<"Computing eigenvalues\n";
   cout.flush();
   dseigen(gme.MH,gme.eigenw,gme.lte+gme.ltm,gme.stat_i,gme.stat_f);
  }
  else{
   gme.eigenvec=new double[(gme.lte+gme.ltm)*(gme.lte+gme.ltm)]();
   cout<<"Computing eigenvalues and eigenvectors\n";
   cout.flush();
   dseigenv(gme.MH,gme.eigenw,gme.eigenvec,gme.lte+gme.ltm,gme.stat_i,gme.stat_f);
   gme.vk=new double[gme.lte+gme.ltm];
  }
 }
 else{ // complex case
  if (gme.cimw!=1 && flim.cfds!=1){
   cout<<"Computing eigenvalues\n";
   cout.flush();
   zheigen(gme.cMH,gme.eigenw,gme.lte+gme.ltm,gme.stat_i,gme.stat_f);
  }
  else{
   gme.ceigenvec=new complex<double>[(gme.lte+gme.ltm)*(gme.lte+gme.ltm)]();
   cout<<"Computing eigenvalues and eigenvectors\n";
   cout.flush();
   zheigenv(gme.cMH,gme.eigenw,gme.ceigenvec,gme.lte+gme.ltm,gme.stat_i,gme.stat_f);
   gme.cvk=new complex<double>[gme.lte+gme.ltm];
  }
 }
}

// dealloc MH and the matrices and vectors from the eigensystem solution

void free_MHeigen (st_gme &gme){
 delete[] gme.eigenw;
 if (gme.psym==1){ // real case
  delete[] gme.MH;
  delete[] gme.eigenvec;
  delete[] gme.vk;
 }
 else{ // complex case
  delete[] gme.cMH;
  delete[] gme.ceigenvec;
  delete[] gme.cvk;
 }
}

/*######################################################################################################*/

/*******************************************************************************/
/**************************    losses computation    ***************************/
/*******************************************************************************/

void alloc_imw (st_gme &gme){
 if (gme.k[3]!=0 && gme.cimw!=0)
  gme.imwv=new double[2*(gme.lstat_f-gme.lstat_i+1)]();
}
 
void computeimw (st_gme &gme, st_eps &eps){
 if (gme.k[3]!=0 && gme.cimw!=0){
  cout<<"Computing losses:\n";
  cout.flush();

  for (int stat=gme.lstat_i;stat<(gme.lstat_f+1);stat++){
   cout<<"state "<<stat;
   cout.flush();
    
   t=clock();
   gme.wk=sqrt(gme.eigenw[stat-gme.stat_i]);
   gme.imwv[stat-gme.lstat_i]=gme.wk;

   if (gme.psym==1){ // real case
    dgetvj(gme.eigenvec,gme.vk,stat-gme.stat_i,gme.lte+gme.ltm);
    Imw(gme);
   }
   else{ // complex case
    zgetvj(gme.ceigenvec,gme.cvk,stat-gme.stat_i,gme.lte+gme.ltm);
    Imw(gme);
   }
   gme.imwv[stat-gme.lstat_i+(gme.lstat_f-gme.lstat_i+1)]=gme.imw;
   t=(clock()-t)/CLOCKS_PER_SEC;
   cout<<"  computed in "<<t<<" seconds\n";
   cout.flush();
  }
 }
} 

void free_imw (st_gme &gme){
 if (gme.k[3]!=0 && gme.cimw!=0)
  delete[] gme.imwv;
}

/*######################################################################################################*/

/********************************************************************************/
/****************************    output messages    *****************************/
/********************************************************************************/ 

void initcheck (st_gme &gme, field_box_lim &flim){
 if (gme.cimw==0) 
  cout<<"\nCOMPUTING PHOTONIC DISPERSION"<<"\n\n";
 else
  cout<<"\nCOMPUTING PHOTONIC DISPERSION AND LOSSES"<<"\n\n";
  cout.flush();

 if (gme.cimw!=0 && (gme.lstat_i<gme.stat_i || gme.stat_f<gme.lstat_i || gme.lstat_f<gme.stat_i || gme.stat_f<gme.lstat_f)){
  cout<<"The sates chosen to calculate the losses are out of the range ["<<gme.stat_i<<","<<gme.stat_f<<"]\n";
  cout.flush();
  abort();
 }
 if (flim.cfds!=0 && (flim.stat_i<gme.stat_i || gme.stat_f<flim.stat_i || flim.stat_f<gme.stat_i || gme.stat_f<flim.stat_f)){
  cout<<"The sates chosen to compute the fields are out of the range ["<<gme.stat_i<<","<<gme.stat_f<<"]\n";
  cout.flush();
  abort();
 }
}

void minfo (st_gme &gme){
 cout<<"Planewaves: "<<gme.npw<<"\n";
 cout<<"Guided modes: "<<gme.Ng<<"\n\n";
 cout<<"Initial state in the eigensolver: "<<gme.stat_i<<"\n";
 cout<<"Final state in the eigensolver: "<<gme.stat_f<<"\n";
 cout<<"The number of bands to obtain is: "<<(gme.stat_f-gme.stat_i+1)<<"\n\n";
 cout.flush();
}

void kmess (st_gme, int kbz){
 cout<<"Computing for k_"<<kbz<<"=("<<gme.k[0]<<","<<gme.k[1]<<")\n";
}

void dimcheck (st_gme &gme, field_box_lim &flim){
 if ((gme.stat_f+1)>(gme.lte+gme.ltm)){
  cout<<"The final state gme.stat_f: "<<gme.stat_f<<"th to compute the photonic dispersion is higher than the dimension of the eigenvalue problem: "<<(gme.lte+gme.ltm)<<"\n";
  cout.flush();
  abort();
 }
 if (gme.cimw!=0){ 
  if (gme.lstat_f>=(gme.lte+gme.ltm)){ /* remember that gme.lstat must be in the interval from 0 to (gme.lte+gme.ltm)-1*/
   cout<<"The final state gme.lstat_f: "<<gme.lstat_f+1<<"th to calculate the losses is higher than the dimension of the eigenvalue problem: "<<(gme.lte+gme.ltm)<<"\n";
   cout.flush();
   abort();
  }
 }
 if (flim.cfds!=0){ 
  if (flim.stat_f>=(gme.lte+gme.ltm)){ /* remember that flim.stat must be in the interval from 0 to (gme.lte+gme.ltm)-1*/
   cout<<"The final state flim.stat_f: "<<flim.stat_f+1<<"th to write the fields is higher than the dimension of the eigenvalue problem :"<<(gme.lte+gme.ltm)<<"\n";
   cout.flush();
   abort();
  }
 }
}

void breakline (int nl){
 for (int i=0;i<nl;i++){
  cout<<"\n";
 }
}

/*######################################################################################################*/

/********************************************************************************/
/*************************    cutoff and plane waves   **************************/
/********************************************************************************/ 

/* This function calculates and write the number of plane waves as a 
 * function of the cutoff values Gmax
 */

void Gmaxvsnpw (st_eps &eps, double Ginit, double Gfinal, double dG){
 int npw;
 double Gval,*G;
 ofstream Gvsnpwout ("Gmax_vs_npw.dat");
 Gvsnpwout.precision(set_pre);

 for (Gval=Ginit;Gval<=Gfinal;Gval=Gval+dG){
  Gi(eps,G,1,npw,Gval,0);
  Gvsnpwout<<Gval<<"	"<<npw<<"\n";
 }

 Gvsnpwout.close();
}

/********************************************************************************/
/*********************************    k lists   **********************************/
/********************************************************************************/ 

/* This function sets the points k=(kx,ky) in the perimeter of the triangle 
 * determined by the three vertices K1, K2 and K3. The list is created in the direction
 * K1 to K2 to K3. There are np points between two consecutive vertices.
 * By default the points K1, K2 and K3 are not included, but it is possible to 
 * include also these points commenting the if conditions (p!=0)
 */

void klist_tri (st_gme &gme, field_box_lim &flim){
 double tl,Dl,m,Dx;
 int lp,f;

 gme.dimlk=3*gme.nkp; // If you consider also de vertices it is necessary to add 3 in the allocation

 if ((gme.klstat_f+1)>gme.dimlk){
  cout<<"The last k state to compute the losses: "<<gme.klstat_f+1<<"th is higher than the total number of k points to compute: "<<gme.dimlk<<"\n";
  cout.flush();
  abort();
 }

 if ((flim.kstat_f+1)>gme.dimlk){
  cout<<"The last k state to compute the fields: "<<gme.klstat_f+1<<"th is higher than the total number of k points to compute: "<<gme.dimlk<<"\n";
  cout.flush();
  abort();
 }

 gme.lk=new double[gme.dimlk*5](); 
 lp=0;
 tl=0;
  
 // K1 to K2
 Dl=sqrt((gme.K2[0]-gme.K1[0])*(gme.K2[0]-gme.K1[0])+(gme.K2[1]-gme.K1[1])*(gme.K2[1]-gme.K1[1]))/(gme.nkp+1);
 for (int p=0;p<(gme.nkp+1);p++){
  if (p!=0){ // if you want to set also the point K1 comment this if (if this you need to allocate gme.lk with the proper dimension)
   if (gme.K2[0]!=gme.K1[0]){
    m=(gme.K2[1]-gme.K1[1])/(gme.K2[0]-gme.K1[0]);
    Dx=Dl/sqrt(m*m+1);
    if (gme.K2[0]<gme.K1[0])
     f=-1;
    else
     f=1;
    gme.lk[lp]=gme.K1[0]+f*Dx*p;
    gme.lk[lp+gme.dimlk]=gme.K1[1]+m*((gme.K1[0]+f*Dx*p)-gme.K1[0]);
    gme.lk[lp+2*gme.dimlk]=tl;
   }
   else{
    if (gme.K2[1]<gme.K1[1])
     f=-1;
    else
     f=1;
    gme.lk[lp]=gme.K1[0];
    gme.lk[lp+gme.dimlk]=gme.K1[1]+f*Dl*p;
    gme.lk[lp+2*gme.dimlk]=tl;
   }
   if (gme.klstat_i<=lp && lp<=gme.klstat_f)
    gme.lk[lp+3*gme.dimlk]=1;
   if (flim.kstat_i<=lp && lp<=flim.kstat_f)
    gme.lk[lp+4*gme.dimlk]=1;
   lp++;
  } // final of the if associated to K1
  tl=tl+Dl;
 }

 // K2 to K3
 Dl=sqrt((gme.K3[0]-gme.K2[0])*(gme.K3[0]-gme.K2[0])+(gme.K3[1]-gme.K2[1])*(gme.K3[1]-gme.K2[1]))/(gme.nkp+1);
 for (int p=0;p<(gme.nkp+1);p++){
  if (p!=0){ // if you want to set also the point K2 comment this if (if this you need to allocate gme.lk with the proper dimension)
   if (gme.K3[0]!=gme.K2[0]){
    m=(gme.K3[1]-gme.K2[1])/(gme.K3[0]-gme.K2[0]);
    Dx=Dl/sqrt(m*m+1);
    if (gme.K3[0]<gme.K2[0])
     f=-1;
    else
     f=1;
    gme.lk[lp]=gme.K2[0]+f*Dx*p;
    gme.lk[lp+gme.dimlk]=gme.K2[1]+m*((gme.K2[0]+f*Dx*p)-gme.K2[0]);
    gme.lk[lp+2*gme.dimlk]=tl;
   }
   else{
    if (gme.K3[1]<gme.K2[1])
     f=-1;
    else
     f=1;
    gme.lk[lp]=gme.K2[0];
    gme.lk[lp+gme.dimlk]=gme.K2[1]+f*Dl*p;
    gme.lk[lp+2*gme.dimlk]=tl;
   }
   if (gme.klstat_i<=lp && lp<=gme.klstat_f)
    gme.lk[lp+3*gme.dimlk]=1;
   if (flim.kstat_i<=lp && lp<=flim.kstat_f)
    gme.lk[lp+4*gme.dimlk]=1;
   lp++;
  } // final of the if associated to K2
  tl=tl+Dl;
 }

 // K3 to K1
 Dl=sqrt((gme.K1[0]-gme.K3[0])*(gme.K1[0]-gme.K3[0])+(gme.K1[1]-gme.K3[1])*(gme.K1[1]-gme.K3[1]))/(gme.nkp+1);
 for (int p=0;p<(gme.nkp+1);p++){
  if (p!=0){ // if you want to set also the point K3 comment this if (if this you need to allocate gme.lk with the proper dimension)
   if (gme.K1[0]!=gme.K3[0]){
    m=(gme.K1[1]-gme.K3[1])/(gme.K1[0]-gme.K3[0]);
    Dx=Dl/sqrt(m*m+1);
    if (gme.K1[0]<gme.K3[0])
     f=-1;
    else
     f=1;
    gme.lk[lp]=gme.K3[0]+f*Dx*p;
    gme.lk[lp+gme.dimlk]=gme.K3[1]+m*((gme.K3[0]+f*Dx*p)-gme.K3[0]);
    gme.lk[lp+2*gme.dimlk]=tl;
   }
   else{
    if (gme.K1[1]<gme.K3[1])
     f=-1;
    else
     f=1;
    gme.lk[lp]=gme.K3[0];
    gme.lk[lp+gme.dimlk]=gme.K3[1]+f*Dl*p;
    gme.lk[lp+2*gme.dimlk]=tl;
   }
   if (gme.klstat_i<=lp && lp<=gme.klstat_f)
    gme.lk[lp+3*gme.dimlk]=1;
   if (flim.kstat_i<=lp && lp<=flim.kstat_f)
    gme.lk[lp+4*gme.dimlk]=1;
   lp++;
  } // final of the if associated to K3
  tl=tl+Dl;
 }
}

/* This function sets the points k=(kx,ky) in the line determined by 
 * the two points K1 and K2. Only K1 and K2 are considered. The list 
 * is created in the direction K1 to K2. There are np points between 
 * K1 and K2. The points K1 are K2 are included!!
 */

void klist_line (st_gme &gme, field_box_lim &flim){
 double tl,Dl,m,Dx;
 int lp,f;

 gme.dimlk=gme.nkp+2; 

 if ((gme.klstat_f+1)>gme.dimlk){
  cout<<"The last k state to compute the losses: "<<gme.klstat_f+1<<"th is higher than the total number of k points to compute: "<<gme.dimlk<<"\n";
  cout.flush();
  abort();
 }

 if ((flim.kstat_f+1)>gme.dimlk){
  cout<<"The last k state to compute the fields: "<<gme.klstat_f+1<<"th is higher than the total number of k points to compute: "<<gme.dimlk<<"\n";
  cout.flush();
  abort();
 }

 gme.lk=new double[gme.dimlk*5](); 
 lp=0;
 tl=0;

 // K1 to K2
 Dl=sqrt((gme.K2[0]-gme.K1[0])*(gme.K2[0]-gme.K1[0])+(gme.K2[1]-gme.K1[1])*(gme.K2[1]-gme.K1[1]))/(gme.nkp+1);
 for (int p=0;p<(gme.nkp+2);p++){
  if (gme.K2[0]!=gme.K1[0]){
   m=(gme.K2[1]-gme.K1[1])/(gme.K2[0]-gme.K1[0]);
   Dx=Dl/sqrt(m*m+1);
   if (gme.K2[0]<gme.K1[0])
    f=-1;
   else
    f=1;
   gme.lk[lp]=gme.K1[0]+f*Dx*p;
   gme.lk[lp+gme.dimlk]=gme.K1[1]+m*((gme.K1[0]+f*Dx*p)-gme.K1[0]);
   gme.lk[lp+2*gme.dimlk]=tl;
  }
  else{
   if (gme.K2[1]<gme.K1[1])
    f=-1;
   else
    f=1;
   gme.lk[lp]=gme.K1[0];
   gme.lk[lp+gme.dimlk]=gme.K1[1]+f*Dl*p;
   gme.lk[lp+2*gme.dimlk]=tl;
  }
  if (gme.klstat_i<=lp && lp<=gme.klstat_f)
   gme.lk[lp+3*gme.dimlk]=1;
  if (flim.kstat_i<=lp && lp<=flim.kstat_f)
   gme.lk[lp+4*gme.dimlk]=1;
  lp++;
 tl=tl+Dl;
 }
}

/* This function set the point k=(kx,ky). Only K1 is considered!! 
 * Here gme.lkstat_i, gme.lkstat_f, flim.kstat_i, flim.kstat_f and gme.nkp
 * are ignored, and the calculation of the losses and the fields
 * is determined by the gme.cimw and the flim.cfds varables, respectively!
 */

void klist_point (st_gme &gme, field_box_lim &flim){

 gme.lk=new double[5]();
 gme.dimlk=1;
 
 gme.lk[0]=gme.K1[0];
 gme.lk[1]=gme.K1[1];
 gme.lk[2]=0;
 gme.lk[3]=1;
 gme.lk[4]=1;
}

// dealloc klist (if you need)

void klist_free (st_gme &gme){
 delete[] gme.lk;
}

/*######################################################################################################*/

/********************************************************************************/
/******************************    output files   *******************************/
/********************************************************************************/

/* These functions have a default output, if you want to configure your
 * own output function it is necessary to create it here or in the main 
 * program
 */ 

// Photonic dispersion (kx,ky,k_brillouin,kmag,w)

void writedis (st_gme &gme){
 
 ofstream Disp ("dispersion.dat",ios_base::app);
 Disp.precision(set_pre);

 for (int i=0;i<(gme.stat_f-gme.stat_i+1);i++){ 
  Disp<<gme.k[0]<<" "<<gme.k[1]<<" "<<gme.k[2]<<" "<<sqrt(gme.k[0]*gme.k[0]+gme.k[1]*gme.k[1])/(2*M_PI)<<" "<<sqrt(gme.eigenw[i])/(2*M_PI)<<"\n";
 }
 Disp.flush();

 Disp.close();
}

// Photonic dispersion with one parameter (parameter,kx,ky,k_brillouin,kmag,w)

void writedis_1p (st_gme &gme, double param){
 
 ofstream Disp ("dispersion_1p.dat",ios_base::app);
 Disp.precision(set_pre);

 for (int i=0;i<(gme.stat_f-gme.stat_i+1);i++){ 
  Disp<<param<<" "<<gme.k[0]<<" "<<gme.k[1]<<" "<<gme.k[2]<<" "<<sqrt(gme.k[0]*gme.k[0]+gme.k[1]*gme.k[1])/(2*M_PI)<<" "<<sqrt(gme.eigenw[i])/(2*M_PI)<<"\n";
 }
 Disp.flush();

 Disp.close();
}

// Photonic dispersion with two parameters (parameter1,parameter2,kx,ky,k_brillouin,kmag,w)

void writedis_2p (st_gme &gme, double param1, double param2){
 
 ofstream Disp ("dispersion_2p.dat",ios_base::app);
 Disp.precision(set_pre);

 for (int i=0;i<(gme.stat_f-gme.stat_i+1);i++){ 
  Disp<<param1<<" "<<param2<<" "<<gme.k[0]<<" "<<gme.k[1]<<" "<<gme.k[2]<<" "<<sqrt(gme.k[0]*gme.k[0]+gme.k[1]*gme.k[1])/(2*M_PI)<<" "<<sqrt(gme.eigenw[i])/(2*M_PI)<<"\n";
 }
 Disp.flush();

 Disp.close();
}

// Diffraction losses (kx,ky,k_brillouin,w,imw)

void writeloss (st_gme &gme){
 if (gme.k[3]!=0 && gme.cimw!=0){
  ofstream Loss ("losses.dat",ios_base::app);
  Loss.precision(set_pre);

  for (int i=0;i<(gme.lstat_f-gme.lstat_i+1);i++){
   Loss<<gme.k[0]<<" "<<gme.k[1]<<" "<<gme.k[2]<<" "<<gme.lstat_i+i<<" "<<gme.imwv[i]/(2*M_PI)<<" "<<gme.imwv[i+(gme.lstat_f-gme.lstat_i+1)]<<"\n";
  }
  Loss.flush();
  Loss.close();
 }
}

// Diffraction losses with one parameter (parameter,kx,ky,k_brillouin,w,imw)

void writeloss_1p (st_gme &gme, double param){
 if (gme.k[3]!=0 && gme.cimw!=0){
  ofstream Loss ("losses_1p.dat",ios_base::app);
  Loss.precision(set_pre);

  for (int i=0;i<(gme.lstat_f-gme.lstat_i+1);i++){
   Loss<<param<<" "<<gme.k[0]<<" "<<gme.k[1]<<" "<<gme.k[2]<<" "<<gme.lstat_i+i<<" "<<gme.imwv[i]/(2*M_PI)<<" "<<gme.imwv[i+(gme.lstat_f-gme.lstat_i+1)]<<"\n";
  }
  Loss.flush();
  Loss.close();
 }
}

// Diffraction losses with two parameters (parameter1,parameter2,kx,ky,k_brillouin,w,imw)

void writeloss_2p (st_gme &gme, double param1, double param2){
 if (gme.k[3]!=0 && gme.cimw!=0){
  ofstream Loss ("losses_2p.dat",ios_base::app);
  Loss.precision(set_pre);

  for (int i=0;i<(gme.lstat_f-gme.lstat_i+1);i++){
   Loss<<param1<<" "<<param2<<" "<<gme.k[0]<<" "<<gme.k[1]<<" "<<gme.k[2]<<" "<<gme.lstat_i+i<<" "<<gme.imwv[i]/(2*M_PI)<<" "<<gme.imwv[i+(gme.lstat_f-gme.lstat_i+1)]<<"\n";
  }
  Loss.flush();
  Loss.close();
 }
}

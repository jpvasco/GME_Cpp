/*******************************************************************************/
/*************************    Dielectric structures   **************************/
/*******************************************************************************/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;

/*######################################################################################################*/

/*******************************************************************/
/******************    Mathematical functions    *******************/
/*******************************************************************/


// Complex exponential 

complex<double> cexp (double arge){
 return complex<double>(cos(arge),sin(arge));
}

// sinc function

double sinc (double x){
 if(x==0)
  return 1;
 else
  return sin(x)/x;
} 

/*######################################################################################################*/

/***************************************************************/
/*************    Reciprocal primitive vectors    **************/
/***************************************************************/

void bvec (double *a1, double *a2, double *b1, double *b2, double &A){
 A=a1[0]*a2[1]-a1[1]*a2[0]; 
 b1[0]=2*M_PI*a2[1]/A;
 b1[1]=-2*M_PI*a2[0]/A;
 b2[0]=-2*M_PI*a1[1]/A;
 b2[1]=2*M_PI*a1[0]/A;
 A=fabs(A); /* area of the unit cell */
}

/*######################################################################################################*/

/***************************************************************/
/**************    General Fourier coefficients   **************/
/***************************************************************/

/* Hexagonal lattice of cylinders with hexagonal cell*/

double Chex (int i, int j, double e_s, double e_r, double r, double *b1, double *b2, double A){
  if (i==0 && j==0)
   return e_s+2*M_PI*(e_r-e_s)*r*r/A;
  else{
    double G=sqrt((i*b1[0]+j*b2[0])*(i*b1[0]+j*b2[0])+(i*b1[1]+j*b2[1])*(i*b1[1]+j*b2[1]));
    return 4*M_PI*(e_r-e_s)*r*j1(r*G)/(A*G);
  }
}

/* Hexagonal lattice of cylinders with rectangular cell*/

double Chexr (int i, int j, double e_s, double e_r, double r, double *b1, double *b2, double A){
  if (i==0 && j==0)
    return e_s+2*M_PI*(e_r-e_s)*r*r/A;
  else{
    double G=sqrt(i*i*b1[0]*b1[0]+j*j*b2[1]*b2[1]);
    return 2*M_PI*(e_r-e_s)*r*(1+cos(i*0.5*b1[0])*cos(j*0.5*sqrt(3)*b2[1]))*j1(r*G)/(A*G);
  }
}

/* Hexagonal lattice of triangular rods with hexagonal cell*/

complex<double> Thex (int i, int j, double e_s, double e_t, double L, double theta, double *b1, double *b2, double A){
  double f,cost,sint,gxp,gyp,j01m,j01p,j02;
  complex<double> Ip,Im;

  f=sqrt(3)*L*L/(4*A);
  cost=cos(theta);
  sint=sin(theta);
  gxp=(i*(b1[0]*cost+b1[1]*sint)+j*(b2[0]*cost+b2[1]*sint))*L/2.0;
  gyp=(i*(b1[1]*cost-b1[0]*sint)+j*(-b2[0]*sint+b2[1]*cost))*sqrt(3)*L/2.0;

  if (gxp==0.0 && gyp==0.0){
    return f*e_t+(1-f)*e_s;
  }
  else{
    if (gyp==0.0){
      Ip=(complex<double>(1.0,0)-exp(complex<double>(0,-gxp)))/(gxp*gxp)-complex<double>(0,1.0/gxp);
      Im=(complex<double>(1.0,0)-exp(complex<double>(0,gxp)))/(gxp*gxp)-complex<double>(0,-1.0/gxp);
    }
    else{
      if ((gxp-gyp)==0.0)
	j01p=1;
      else
	j01p=2*sin(0.5*(gxp-gyp))/(gxp-gyp);
      if ((gxp+gyp)==0.0)
	j01m=1;
      else
	j01m=2*sin(0.5*(gxp+gyp))/(gxp+gyp);
      if (gxp==0.0)
	j02=1;
      else
	j02=2*sin(0.5*gxp)/gxp;
      
      Ip=complex<double>(0,1.0/gyp)*exp(complex<double>(0,gyp/3.0-0.5*gxp))*(exp(complex<double>(0,-0.5*gyp))*j01p-j02);
      Im=complex<double>(0,1.0/gyp)*exp(complex<double>(0,gyp/3.0+0.5*gxp))*(exp(complex<double>(0,-0.5*gyp))*j01m-j02);
    }
  return f*(e_t-e_s)*(Ip+Im);
  }
}

/* Fourier transform of triangular rods */

complex<double> Tri (int i, int j, double e_s, double e_t, double L, double theta, double *b1, double *b2, double A){
  double f,cost,sint,gxp,gyp,gxp2,j01m,j01p,j02,g;
  complex<double> Ip,Im,exp05gy;

  f=sqrt(3)*L*L/(4*A);
  cost=cos(theta);
  sint=sin(theta);
  gxp=(i*(b1[0]*cost+b1[1]*sint)+j*(b2[0]*cost+b2[1]*sint))*L/2.0;
  gyp=(i*(b1[1]*cost-b1[0]*sint)+j*(-b2[0]*sint+b2[1]*cost))*sqrt(3)*L/2.0;
  gxp2=gxp*gxp;
  g=sqrt(gxp2+gyp*gyp);

  if (g<0.000000000001){
    return f*e_t+(1-f)*e_s;
  }
  else{
    if (fabs(gyp)<0.000000000001){
      Ip=(complex<double>(1.0,0)-cexp(-gxp))/(gxp2)-complex<double>(0,1.0/gxp);
      Im=(complex<double>(1.0,0)-cexp(gxp))/(gxp2)+complex<double>(0,1.0/gxp);
    }
    else{
     j01p=sinc(0.5*(gxp-gyp));
     j01m=sinc(0.5*(gxp+gyp));
     j02=sinc(0.5*gxp);
     exp05gy=cexp(-0.5*gyp);    
     Ip=complex<double>(0,1.0/gyp)*cexp(gyp/3.0-0.5*gxp)*(exp05gy*j01p-j02);
     Im=complex<double>(0,1.0/gyp)*cexp(gyp/3.0+0.5*gxp)*(exp05gy*j01m-j02);
    }
  return f*(e_t-e_s)*(Ip+Im);
  }
}

/* Fourier transform of polygons rods of N sides*/

complex<double> Polgme (int i, int j, double e_s, double e_p, double l, int N, double th, double *b1, double *b2, double A){
 double ang,L,x1,y1,gxp,gyp;
 complex<double> coef,coeft(0,0);

 L=2.0*l*sin(M_PI/N);
 x1=sqrt(l*l-0.25*L*L);
 y1=0.5*L;

 for (double thi=0;thi<1.999999*M_PI;thi=thi+2.0*M_PI/N){
  ang=thi+0.5*M_PI+th;

  gxp=i*(b1[0]*cos(ang)+b1[1]*sin(ang))+j*(b2[0]*cos(ang)+b2[1]*sin(ang));
  gyp=i*(b1[1]*cos(ang)-b1[0]*sin(ang))+j*(-b2[0]*sin(ang)+b2[1]*cos(ang));

  coef=complex<double>(0,0);

 if (fabs(gyp)<0.000000000001 && fabs(gxp)<0.000000000001)
  coef=e_s/N+(e_p-e_s)*x1*y1/A;
 else if (fabs(gyp)<0.000000000001)
  coef=2.0*(e_p-e_s)*y1*(cexp(gxp*x1)*complex<double>(1,-gxp*x1)-complex<double>(1,0))/(x1*gxp*gxp*A);
 else if (fabs(gxp*gxp*x1*x1-gyp*gyp*y1*y1)<0.000000000001)
  coef=(e_p-e_s)*y1*(0.5*(complex<double>(1,0)-cexp(2.0*gxp*x1))+complex<double>(0,gxp*x1))/(gxp*gxp*x1*A);
 else
  coef=2*(e_p-e_s)*x1*(cexp(gxp*x1)*complex<double>(gyp*y1*cos(gyp*y1),-gxp*x1*sin(gyp*y1))-gyp*y1)/((gxp*x1-gyp*y1)*(gxp*x1+gyp*y1)*gyp*A);

 coeft=coeft+coef; 
 }
 
 return coeft;
}

/* Rectangular lattice of cylinders with rectangular cell*/

double Crec (int i, int j, double e_s, double e_r, double r, double *b1, double *b2, double A){
  if (i==0 && j==0)
    return e_s+M_PI*(e_r-e_s)*r*r/A;
  else{
    double G=sqrt(i*i*b1[0]*b1[0]+j*j*b2[1]*b2[1]);
    return 2*M_PI*(e_r-e_s)*r*j1(r*G)/(A*G);
  }
}

/*######################################################################################################*/

/***************************************************************/
/***************   Dielectric data structures   ****************/
/***************************************************************/

struct st_eps {

 struct st_recreg { /* regular rectangular lattice of cylinders */
  double es,er,r,lx,ly,a1[2],a2[2],b1[2],b2[2],A;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
  }
  st_recreg():es(0),er(0),r(0),lx(0),ly(0){}
 }recreg;

 struct st_hexreg { /* regular hexagonal lattice of cylinders */
  double es,er,r,a1[2],a2[2],b1[2],b2[2],A;
  void init(){ 
   a1[0]=1;a1[1]=0;a2[0]=0.5;a2[1]=0.5*sqrt(3);
   bvec(a1,a2,b1,b2,A);
   A=2*A;
  }
  st_hexreg():es(0),er(0),r(0){}
 }hexreg;
 
 struct st_Thexreg { /* regular hexagonal lattice of triangles */
  double es,et,L,theta,a1[2],a2[2],b1[2],b2[2],A;
  void init(){ 
   a1[0]=1;a1[1]=0;a2[0]=0.5;a2[1]=0.5*sqrt(3);
   bvec(a1,a2,b1,b2,A);
  }
 st_Thexreg():es(0),et(0),L(0),theta(0){}
 }Thexreg;
 
 struct st_rech1 { /* H1 defect in a rectangular lattice of cylinders */
  double es,er,r,lxp,lyp,lx,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=lxp;ap1[1]=0;ap2[0]=0;ap2[1]=lyp;
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_rech1():es(0),er(0),r(0),lxp(0),lyp(0),lx(0),ly(0){}
 }rech1;

 struct st_hexh1 { /* H1 defect in a hexagonal lattice of cylinders */
  double es,er,r,l,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=l;a1[1]=0;a2[0]=0.5*l;a2[1]=0.5*l*sqrt(3);
   bvec(a1,a2,b1,b2,A);
   A=2*A;
   ap1[0]=1;ap1[1]=0;ap2[0]=0.5;ap2[1]=0.5*sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
   Ap=2*Ap;
  }
  st_hexh1():es(0),er(0),r(0),l(0){}
 }hexh1;

 struct st_hexh1r { /* H1 defect in a hexagonal lattice of cylinders with rectangular cell*/
  double es,er,r,lx,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=1;ap1[1]=0;ap2[0]=0;ap2[1]=sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_hexh1r():es(0),er(0),r(0),lx(0),ly(0){}
 }hexh1r;

 struct st_recl2 { /* L2 defect in a rectangular lattice of cylinders */
  double es,er,r,lxp,lyp,lx,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=lxp;ap1[1]=0;ap2[0]=0;ap2[1]=lyp;
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_recl2():es(0),er(0),r(0),lxp(0),lyp(0),lx(0),ly(0){}
 }recl2;

 struct st_hexl2 { /* L2 defect in a hexagonal lattice of cylinders */
  double es,er,r,l,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=l;a1[1]=0;a2[0]=0.5*l;a2[1]=0.5*l*sqrt(3);
   bvec(a1,a2,b1,b2,A);
   A=2*A;
   ap1[0]=1;ap1[1]=0;ap2[0]=0.5;ap2[1]=0.5*sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
   Ap=2*Ap;
  }
  st_hexl2():es(0),er(0),r(0),l(0){}
 }hexl2;

 struct st_hexl2r { /* L2 defect in a hexagonal lattice of cylinders with rectangular cell*/
  double es,er,r,lx,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=1;ap1[1]=0;ap2[0]=0;ap2[1]=sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_hexl2r():es(0),er(0),r(0),lx(0),ly(0){}
 }hexl2r;

 struct st_recl3 { /* L3 defect in a rectangular lattice of cylinders */
  double es,er,r,s,rs,lxp,lyp,lx,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=lxp;ap1[1]=0;ap2[0]=0;ap2[1]=lyp;
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_recl3():es(0),er(0),r(0),s(0),rs(0),lxp(0),lyp(0),lx(0),ly(0){}
 }recl3;

 struct st_hexl3 { /* L3 defect in a hexagonal lattice of cylinders */
  double es,er,r,s,rs,l,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=l;a1[1]=0;a2[0]=0.5*l;a2[1]=0.5*l*sqrt(3);
   bvec(a1,a2,b1,b2,A);
   A=2*A;
   ap1[0]=1;ap1[1]=0;ap2[0]=0.5;ap2[1]=0.5*sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
   Ap=2*Ap;
  }
  st_hexl3():es(0),er(0),r(0),s(0),rs(0),l(0){}
 }hexl3;

 struct st_hexl3r { /* L3 defect in a hexagonal lattice of cylinders with rectangular cell*/
  double es,er,r,s,rs,lx,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=1;ap1[1]=0;ap2[0]=0;ap2[1]=sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_hexl3r():es(0),er(0),r(0),s(0),rs(0),lx(0),ly(0){}
 }hexl3r;

 struct st_hexl3sr { /* L3 defect with outward displacement of the lateral holes in a hexagonal lattice of cylinders with rectangular cell*/
  double es,er,r,s1,rs1,s2,rs2,s3,rs3,s4,rs4,s5,rs5,lx,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=1;ap1[1]=0;ap2[0]=0;ap2[1]=sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_hexl3sr():es(0),er(0),r(0),s1(0),rs1(0),s2(0),rs2(0),s3(0),rs3(0),s4(0),rs4(0),s5(0),rs5(0),lx(0),ly(0){}
 }hexl3sr;

 struct st_recwg { /* Waveguide defect in a rectangular lattice of cylinders */
  double es,er,r,lyp,ax,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=ax;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=ax;ap1[1]=0;ap2[0]=0;ap2[1]=lyp;
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_recwg():es(0),er(0),r(0),lyp(0),ax(0),ly(0){}
 }recwg;

struct st_Thexwgr { /* Waveguide defect in a hexagonal lattice of triangles with rectangular cell */
  int Nh;
  double es,et,L,theta,ly,a1[2],a2[2],b1[2],b2[2],A,*Mch;
  void init(){ 
   int py,p=0;
   a1[0]=1;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   py=(int)round(ly/sqrt(3));
   Nh=3*py-1;
   Mch=new double[3*Nh];

   for (double i=(-0.5*py+1);i<=(0.5*py-1);i++){
    Mch[p]=0.0;
    Mch[p+Nh]=i*sqrt(3);
    Mch[p+2*Nh]=et-es;
    p++;
   }
   Mch[p]=0.0;
   Mch[p+Nh]=-0.5*py*sqrt(3);
   Mch[p+2*Nh]=0.5*(et-es);
   p++;
   Mch[p]=0.0;
   Mch[p+Nh]=0.5*py*sqrt(3);
   Mch[p+2*Nh]=0.5*(et-es);
   p++;
   for (double i=1;i<=0.5*(py-1);i++){
    Mch[p]=0.5;
    Mch[p+Nh]=i*sqrt(3);
    Mch[p+2*Nh]=0.5*(et-es);
    p++;
    Mch[p]=-0.5;
    Mch[p+Nh]=i*sqrt(3);
    Mch[p+2*Nh]=0.5*(et-es);
    p++;
    Mch[p]=0.5;
    Mch[p+Nh]=-i*sqrt(3);
    Mch[p+2*Nh]=0.5*(et-es);
    p++;
    Mch[p]=-0.5;
    Mch[p+Nh]=-i*sqrt(3);
    Mch[p+2*Nh]=0.5*(et-es);
    p++;
   }
  }
  void free(){
   delete[] Mch;
  }
 st_Thexwgr():es(0),et(0),L(0),theta(0),ly(0){}
 }Thexwgr;

struct st_Phexwgr { /* Waveguide defect in a hexagonal lattice of polygons with rectangular cell */
  int Nh,N;
  double es,ep,l,th,ly,a1[2],a2[2],b1[2],b2[2],A,*Mch;
  void init(){ 
   int py,p=0;
   a1[0]=1;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   py=(int)round(ly/sqrt(3));
   Nh=3*py-1;
   Mch=new double[3*Nh];

   for (double i=(-0.5*py+1);i<=(0.5*py-1);i++){
    Mch[p]=0.0;
    Mch[p+Nh]=i*sqrt(3);
    Mch[p+2*Nh]=ep-es;
    p++;
   }
   Mch[p]=0.0;
   Mch[p+Nh]=-0.5*py*sqrt(3);
   Mch[p+2*Nh]=0.5*(ep-es);
   p++;
   Mch[p]=0.0;
   Mch[p+Nh]=0.5*py*sqrt(3);
   Mch[p+2*Nh]=0.5*(ep-es);
   p++;
   for (double i=1;i<=0.5*(py-1);i++){
    Mch[p]=0.5;
    Mch[p+Nh]=i*sqrt(3);
    Mch[p+2*Nh]=0.5*(ep-es);
    p++;
    Mch[p]=-0.5;
    Mch[p+Nh]=i*sqrt(3);
    Mch[p+2*Nh]=0.5*(ep-es);
    p++;
    Mch[p]=0.5;
    Mch[p+Nh]=-i*sqrt(3);
    Mch[p+2*Nh]=0.5*(ep-es);
    p++;
    Mch[p]=-0.5;
    Mch[p+Nh]=-i*sqrt(3);
    Mch[p+2*Nh]=0.5*(ep-es);
    p++;
   }
  }
  void free(){
   delete[] Mch;
  }
 st_Phexwgr():es(0),ep(0),l(0),N(0),th(0),ly(0){}
 }Phexwgr;

 struct st_hexwgr { /* Waveguide defect in a hexagonal lattice of cylinders with rectangular cell */
  double es,er,r,ly,a1[2],a2[2],b1[2],b2[2],A,ap1[2],ap2[2],bp1[2],bp2[2],Ap;
  void init(){ 
   a1[0]=1;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   ap1[0]=1;ap1[1]=0;ap2[0]=0;ap2[1]=sqrt(3);
   bvec(ap1,ap2,bp1,bp2,Ap);
  }
  st_hexwgr():es(0),er(0),r(0),ly(0){}
 }hexwgr;
 
 struct st_dishexl3r { /* L3 defect in a hexagonal lattice of cylinders with disorder and rectangular cell*/
  int Nh;
  double es,er,r,s,rs,lx,ly,sigma,a1[2],a2[2],b1[2],b2[2],A,*Mch;
  void init(){
   int px,py,p=0;
   double rr;
   a1[0]=lx;a1[1]=0;a2[0]=0;a2[1]=ly;
   bvec(a1,a2,b1,b2,A);
   px=(int)lx;
   py=(int)round(ly/sqrt(3));
   gsl_rng *rn;
   rn=gsl_rng_alloc(gsl_rng_mt19937);
   gsl_rng_set(rn,time(NULL));
   Nh=(px+1)*(py+1)+px*py-1-3;
   Mch=new double[4*Nh];

   if(px%2==py%2){
    cout<<"The values of lx and ly/sqrt(3) cannot be both odd or both even\n";
    abort();
   }
   
   // if you want to see the gaussian distribution uncomment the 4 lines below
   //ofstream gauss ("gaussian-dishexl3r.dat");
   //for (int i=0;i<100000;i++)
   // gauss<<gsl_ran_gaussian(rn,sigma)+r<<"\n";
   //gauss.close();
   
   ofstream gaussr ("gaussianr-dishexl3r.dat"); // random values of r  
   for (double i=-px/2.0;i<=px/2.0;i++){
    for (double j=(-py/2.0+0.5);j<=(py/2.0-0.5);j++){

     if ((fabs(j)<0.001 && fabs(i)<0.001) || (fabs(j)<0.001 && fabs(fabs(i)-1)<0.001) || (fabs(j)<0.001 && fabs(fabs(i)-2)<0.001))
      continue;
     else{
      if (fabs(fabs(i)-lx/2.0)<0.00000001 && fabs(fabs(j*sqrt(3))-ly/2.0)<0.00000001){
       Mch[p]=i;
       Mch[p+Nh]=j*sqrt(3);
       Mch[p+2*Nh]=r;
       Mch[p+3*Nh]=(er+3*es)/4.0;     
       p++;
      }
      else if (fabs(fabs(i)-lx/2.0)<0.00000001 || fabs(fabs(j*sqrt(3))-ly/2.0)<0.00000001){
       Mch[p]=i;
       Mch[p+Nh]=j*sqrt(3);
       Mch[p+2*Nh]=r;
       Mch[p+3*Nh]=(er+es)/2.0;     
       p++;
      }
      else{
       rr=r+gsl_ran_gaussian(rn,sigma);
       Mch[p]=i;
       Mch[p+Nh]=j*sqrt(3);
       Mch[p+2*Nh]=rr;
       Mch[p+3*Nh]=er;
       gaussr<<rr<<"\n";     
       p++;
      }
     }
    }
   }
   for (double i=(-px/2.0+0.5);i<=(px/2.0-0.5);i++){
    for (double j=-py/2.0;j<=py/2.0;j++){

     if ((fabs(j)<0.001 && fabs(i)<0.001) || (fabs(j)<0.001 && fabs(fabs(i)-1)<0.001) || (fabs(j)<0.001 && fabs(fabs(i)-2)<0.001))
      continue;
     else{
      if (fabs(fabs(i)-lx/2.0)<0.00000001 && fabs(fabs(j*sqrt(3))-ly/2.0)<0.00000001){
       Mch[p]=i;
       Mch[p+Nh]=j*sqrt(3);
       Mch[p+2*Nh]=r;
       Mch[p+3*Nh]=(er+3*es)/4.0;     
       p++;
      }
      else if (fabs(fabs(i)-lx/2.0)<0.00000001 || fabs(fabs(j*sqrt(3))-ly/2.0)<0.00000001){
       Mch[p]=i;
       Mch[p+Nh]=j*sqrt(3);
       Mch[p+2*Nh]=r;
       Mch[p+3*Nh]=(er+es)/2.0;     
       p++;
      }
      else{
       rr=r+gsl_ran_gaussian(rn,sigma);
       Mch[p]=i;
       Mch[p+Nh]=j*sqrt(3);
       Mch[p+2*Nh]=rr;
       Mch[p+3*Nh]=er;    
       gaussr<<rr<<"\n"; 
       p++;
      }
     }
    }
   }
   rr=rs+gsl_ran_gaussian(rn,sigma);
   Mch[p]=2+s;
   Mch[p+Nh]=0;
   Mch[p+2*Nh]=rr;
   Mch[p+3*Nh]=er;
   gaussr<<rr<<"\n";
   p++;
   rr=rs+gsl_ran_gaussian(rn,sigma);
   Mch[p]=-2-s;
   Mch[p+Nh]=0;
   Mch[p+2*Nh]=rr;
   Mch[p+3*Nh]=er; 
   gaussr<<rr<<"\n";
   p++;

   gaussr.close();
   gsl_rng_free(rn);
  }
  void free(){
   delete[] Mch;
  }
  st_dishexl3r():es(0),er(0),r(0),s(0),rs(0),lx(0),ly(0),sigma(0){}
 }dishexl3r;  

}eps;

/*######################################################################################################*/

/***************************************************************/
/********************   Area of the cell   *********************/
/***************************************************************/

double Acell (st_eps &eps){
 if (!(eps.recreg.es==0&&eps.recreg.er==0&&eps.recreg.r==0&&eps.recreg.lx==0&&eps.recreg.ly==0))
  return eps.recreg.A;

 else if (!(eps.hexreg.es==0&&eps.hexreg.er==0&&eps.hexreg.r==0))
  return eps.hexreg.A;
 
 else if (!(eps.Thexreg.es==0&&eps.Thexreg.et==0&&eps.Thexreg.L==0&&eps.Thexreg.theta==0))
  return eps.Thexreg.A;
 
 else if (!(eps.rech1.es==0&&eps.rech1.er==0&&eps.rech1.r==0&&eps.rech1.lxp==0&&eps.rech1.lyp==0&&eps.rech1.lx==0&&eps.rech1.ly==0))
  return eps.rech1.A;

 else if (!(eps.hexh1.es==0&&eps.hexh1.er==0&&eps.hexh1.r==0&&eps.hexh1.l==0))
  return eps.hexh1.A;

 else if (!(eps.hexh1r.es==0&&eps.hexh1r.er==0&&eps.hexh1r.r==0&&eps.hexh1r.lx==0&&eps.hexh1r.ly==0))
  return eps.hexh1r.A;

 else if (!(eps.recl2.es==0&&eps.recl2.er==0&&eps.recl2.r==0&&eps.recl2.lxp==0&&eps.recl2.lyp==0&&eps.recl2.lx==0&&eps.recl2.ly==0))
  return eps.recl2.A;

 else if (!(eps.hexl2.es==0&&eps.hexl2.er==0&&eps.hexl2.r==0&&eps.hexl2.l==0))
  return eps.hexl2.A;

 else if (!(eps.hexl2r.es==0&&eps.hexl2r.er==0&&eps.hexl2r.r==0&&eps.hexl2r.lx==0&&eps.hexl2r.ly==0))
  return eps.hexl2r.A;

 else if (!(eps.recl3.es==0&&eps.recl3.er==0&&eps.recl3.r==0&&eps.recl3.s==0&&eps.recl3.rs==0&&eps.recl3.lxp==0&&eps.recl3.lyp==0&&eps.recl3.lx==0&&eps.recl3.ly==0))
  return eps.recl3.A;

 else if (!(eps.hexl3.es==0&&eps.hexl3.er==0&&eps.hexl3.r==0&&eps.hexl3.s==0&&eps.hexl3.rs==0&&eps.hexl3.l==0))
  return eps.hexl3.A;

 else if (!(eps.hexl3r.es==0&&eps.hexl3r.er==0&&eps.hexl3r.r==0&&eps.hexl3r.s==0&&eps.hexl3r.rs==0&&eps.hexl3r.lx==0&&eps.hexl3r.ly==0))
  return eps.hexl3r.A;

 else if (!(eps.hexl3sr.es==0&&eps.hexl3sr.er==0&&eps.hexl3sr.r==0&&eps.hexl3sr.s1==0&&eps.hexl3sr.rs1==0&&eps.hexl3sr.s2==0&&eps.hexl3sr.rs2==0&&eps.hexl3sr.s3==0&&eps.hexl3sr.rs3==0&&eps.hexl3sr.s4==0&&eps.hexl3sr.rs4==0&&eps.hexl3sr.s5==0&&eps.hexl3sr.rs5==0&&eps.hexl3sr.lx==0&&eps.hexl3sr.ly==0))
  return eps.hexl3sr.A;

 else if (!(eps.recwg.es==0&&eps.recwg.er==0&&eps.recwg.r==0&&eps.recwg.ax==0&&eps.recwg.ly==0))
  return eps.recwg.A;

 else if (!(eps.Thexwgr.es==0&&eps.Thexwgr.et==0&&eps.Thexwgr.L==0&&eps.Thexwgr.theta==0&&eps.Thexwgr.ly==0))
  return eps.Thexwgr.A;

 else if (!(eps.Phexwgr.es==0&&eps.Phexwgr.ep==0&&eps.Phexwgr.l==0&&eps.Phexwgr.N==0&&eps.Phexwgr.th==0&&eps.Phexwgr.ly==0))
  return eps.Phexwgr.A;

 else if (!(eps.hexwgr.es==0&&eps.hexwgr.er==0&&eps.hexwgr.r==0&&eps.hexwgr.ly==0))
  return eps.hexwgr.A;
 
 else if (!(eps.dishexl3r.es==0&&eps.dishexl3r.er==0&&eps.dishexl3r.r==0&&eps.dishexl3r.s==0&&eps.dishexl3r.rs==0&&eps.dishexl3r.lx==0&&eps.dishexl3r.ly==0&&eps.dishexl3r.sigma==0))
  return eps.dishexl3r.A;
}

/*######################################################################################################*/

/***************************************************************/
/***************   Reciprocal lattice vector    ****************/
/***************************************************************/

void Gi (st_eps &eps, double *Gif, int qd, int &npw, double Gmax, int n){
 int imax,jmax;
 double nG;
 npw=0;

 if (!(eps.recreg.es==0&&eps.recreg.er==0&&eps.recreg.r==0&&eps.recreg.lx==0&&eps.recreg.ly==0)){
  if (eps.recreg.b1[0]>eps.recreg.b2[1]){
   imax=(int)(Gmax/eps.recreg.b1[0]+1);
   jmax=(int)(imax*eps.recreg.b1[0]/eps.recreg.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.recreg.b2[1]+1);
   imax=(int)(jmax*eps.recreg.b2[1]/eps.recreg.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.recreg.b1[0]*eps.recreg.b1[0]+j*j*eps.recreg.b2[1]*eps.recreg.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.recreg.b1[0];
      Gif[npw+n]=j*eps.recreg.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexreg.es==0&&eps.hexreg.er==0&&eps.hexreg.r==0)){
  imax=2*((int)(Gmax/sqrt(eps.hexreg.b1[0]*eps.hexreg.b1[0]+eps.hexreg.b1[1]*eps.hexreg.b1[1])+1));
  jmax=imax;
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt((i*eps.hexreg.b1[0]+j*eps.hexreg.b2[0])*(i*eps.hexreg.b1[0]+j*eps.hexreg.b2[0])+(i*eps.hexreg.b1[1]+j*eps.hexreg.b2[1])*(i*eps.hexreg.b1[1]+j*eps.hexreg.b2[1]));   
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=(i*eps.hexreg.b1[0]+j*eps.hexreg.b2[0]);
      Gif[npw+n]=(i*eps.hexreg.b1[1]+j*eps.hexreg.b2[1]);
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }
 
 else if (!(eps.Thexreg.es==0&&eps.Thexreg.et==0&&eps.Thexreg.L==0&&eps.Thexreg.theta==0)){
  imax=2*((int)(Gmax/sqrt(eps.Thexreg.b1[0]*eps.Thexreg.b1[0]+eps.Thexreg.b1[1]*eps.Thexreg.b1[1])+1));
  jmax=imax;
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt((i*eps.Thexreg.b1[0]+j*eps.Thexreg.b2[0])*(i*eps.Thexreg.b1[0]+j*eps.Thexreg.b2[0])+(i*eps.Thexreg.b1[1]+j*eps.Thexreg.b2[1])*(i*eps.Thexreg.b1[1]+j*eps.Thexreg.b2[1]));   
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=(i*eps.Thexreg.b1[0]+j*eps.Thexreg.b2[0]);
      Gif[npw+n]=(i*eps.Thexreg.b1[1]+j*eps.Thexreg.b2[1]);
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.rech1.es==0&&eps.rech1.er==0&&eps.rech1.r==0&&eps.rech1.lxp==0&&eps.rech1.lyp==0&&eps.rech1.lx==0&&eps.rech1.ly==0)){
  if (eps.rech1.b1[0]>eps.rech1.b2[1]){
   imax=(int)(Gmax/eps.rech1.b1[0]+1);
   jmax=(int)(imax*eps.rech1.b1[0]/eps.rech1.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.rech1.b2[1]+1);
   imax=(int)(jmax*eps.rech1.b2[1]/eps.rech1.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.rech1.b1[0]*eps.rech1.b1[0]+j*j*eps.rech1.b2[1]*eps.rech1.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.rech1.b1[0];
      Gif[npw+n]=j*eps.rech1.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexh1.es==0&&eps.hexh1.er==0&&eps.hexh1.r==0&&eps.hexh1.l==0)){
  imax=2*((int)(Gmax/sqrt(eps.hexh1.b1[0]*eps.hexh1.b1[0]+eps.hexh1.b1[1]*eps.hexh1.b1[1])+1));
  jmax=imax;
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt((i*eps.hexh1.b1[0]+j*eps.hexh1.b2[0])*(i*eps.hexh1.b1[0]+j*eps.hexh1.b2[0])+(i*eps.hexh1.b1[1]+j*eps.hexh1.b2[1])*(i*eps.hexh1.b1[1]+j*eps.hexh1.b2[1]));           
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=(i*eps.hexh1.b1[0]+j*eps.hexh1.b2[0]);
      Gif[npw+n]=(i*eps.hexh1.b1[1]+j*eps.hexh1.b2[1]);
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexh1r.es==0&&eps.hexh1r.er==0&&eps.hexh1r.r==0&&eps.hexh1r.lx==0&&eps.hexh1r.ly==0)){
  if (eps.hexh1r.b1[0]>eps.hexh1r.b2[1]){
   imax=(int)(Gmax/eps.hexh1r.b1[0]+1);
   jmax=(int)(imax*eps.hexh1r.b1[0]/eps.hexh1r.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.hexh1r.b2[1]+1);
   imax=(int)(jmax*eps.hexh1r.b2[1]/eps.hexh1r.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.hexh1r.b1[0]*eps.hexh1r.b1[0]+j*j*eps.hexh1r.b2[1]*eps.hexh1r.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.hexh1r.b1[0];
      Gif[npw+n]=j*eps.hexh1r.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.recl2.es==0&&eps.recl2.er==0&&eps.recl2.r==0&&eps.recl2.lxp==0&&eps.recl2.lyp==0&&eps.recl2.lx==0&&eps.recl2.ly==0)){
  if (eps.recl2.b1[0]>eps.recl2.b2[1]){
   imax=(int)(Gmax/eps.recl2.b1[0]+1);
   jmax=(int)(imax*eps.recl2.b1[0]/eps.recl2.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.recl2.b2[1]+1);
   imax=(int)(jmax*eps.recl2.b2[1]/eps.recl2.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.recl2.b1[0]*eps.recl2.b1[0]+j*j*eps.recl2.b2[1]*eps.recl2.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.recl2.b1[0];
      Gif[npw+n]=j*eps.recl2.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexl2.es==0&&eps.hexl2.er==0&&eps.hexl2.r==0&&eps.hexl2.l==0)){
  imax=2*((int)(Gmax/sqrt(eps.hexl2.b1[0]*eps.hexl2.b1[0]+eps.hexl2.b1[1]*eps.hexl2.b1[1])+1));
  jmax=imax;
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt((i*eps.hexl2.b1[0]+j*eps.hexl2.b2[0])*(i*eps.hexl2.b1[0]+j*eps.hexl2.b2[0])+(i*eps.hexl2.b1[1]+j*eps.hexl2.b2[1])*(i*eps.hexl2.b1[1]+j*eps.hexl2.b2[1]));           
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=(i*eps.hexl2.b1[0]+j*eps.hexl2.b2[0]);
      Gif[npw+n]=(i*eps.hexl2.b1[1]+j*eps.hexl2.b2[1]);
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexl2r.es==0&&eps.hexl2r.er==0&&eps.hexl2r.r==0&&eps.hexl2r.lx==0&&eps.hexl2r.ly==0)){
  if (eps.hexl2r.b1[0]>eps.hexl2r.b2[1]){
   imax=(int)(Gmax/eps.hexl2r.b1[0]+1);
   jmax=(int)(imax*eps.hexl2r.b1[0]/eps.hexl2r.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.hexl2r.b2[1]+1);
   imax=(int)(jmax*eps.hexl2r.b2[1]/eps.hexl2r.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.hexl2r.b1[0]*eps.hexl2r.b1[0]+j*j*eps.hexl2r.b2[1]*eps.hexl2r.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.hexl2r.b1[0];
      Gif[npw+n]=j*eps.hexl2r.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.recl3.es==0&&eps.recl3.er==0&&eps.recl3.r==0&&eps.recl3.s==0&&eps.recl3.rs==0&&eps.recl3.lxp==0&&eps.recl3.lyp==0&&eps.recl3.lx==0&&eps.recl3.ly==0)){
  if (eps.recl3.b1[0]>eps.recl3.b2[1]){
   imax=(int)(Gmax/eps.recl3.b1[0]+1);
   jmax=(int)(imax*eps.recl3.b1[0]/eps.recl3.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.recl3.b2[1]+1);
   imax=(int)(jmax*eps.recl3.b2[1]/eps.recl3.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.recl3.b1[0]*eps.recl3.b1[0]+j*j*eps.recl3.b2[1]*eps.recl3.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.recl3.b1[0];
      Gif[npw+n]=j*eps.recl3.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexl3.es==0&&eps.hexl3.er==0&&eps.hexl3.r==0&&eps.hexl3.s==0&&eps.hexl3.rs==0&&eps.hexl3.l==0)){
  imax=2*((int)(Gmax/sqrt(eps.hexl3.b1[0]*eps.hexl3.b1[0]+eps.hexl3.b1[1]*eps.hexl3.b1[1])+1));
  jmax=imax;
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt((i*eps.hexl3.b1[0]+j*eps.hexl3.b2[0])*(i*eps.hexl3.b1[0]+j*eps.hexl3.b2[0])+(i*eps.hexl3.b1[1]+j*eps.hexl3.b2[1])*(i*eps.hexl3.b1[1]+j*eps.hexl3.b2[1]));           
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=(i*eps.hexl3.b1[0]+j*eps.hexl3.b2[0]);
      Gif[npw+n]=(i*eps.hexl3.b1[1]+j*eps.hexl3.b2[1]);
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexl3r.es==0&&eps.hexl3r.er==0&&eps.hexl3r.r==0&&eps.hexl3r.s==0&&eps.hexl3r.rs==0&&eps.hexl3r.lx==0&&eps.hexl3r.ly==0)){
  if (eps.hexl3r.b1[0]>eps.hexl3r.b2[1]){
   imax=(int)(Gmax/eps.hexl3r.b1[0]+1);
   jmax=(int)(imax*eps.hexl3r.b1[0]/eps.hexl3r.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.hexl3r.b2[1]+1);
   imax=(int)(jmax*eps.hexl3r.b2[1]/eps.hexl3r.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.hexl3r.b1[0]*eps.hexl3r.b1[0]+j*j*eps.hexl3r.b2[1]*eps.hexl3r.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.hexl3r.b1[0];
      Gif[npw+n]=j*eps.hexl3r.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexl3sr.es==0&&eps.hexl3sr.er==0&&eps.hexl3sr.r==0&&eps.hexl3sr.s1==0&&eps.hexl3sr.rs1==0&&eps.hexl3sr.s2==0&&eps.hexl3sr.rs2==0&&eps.hexl3sr.s3==0&&eps.hexl3sr.rs3==0&&eps.hexl3sr.s4==0&&eps.hexl3sr.rs4==0&&eps.hexl3sr.s5==0&&eps.hexl3sr.rs5==0&&eps.hexl3sr.lx==0&&eps.hexl3sr.ly==0)){
  if (eps.hexl3sr.b1[0]>eps.hexl3sr.b2[1]){
   imax=(int)(Gmax/eps.hexl3sr.b1[0]+1);
   jmax=(int)(imax*eps.hexl3sr.b1[0]/eps.hexl3sr.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.hexl3sr.b2[1]+1);
   imax=(int)(jmax*eps.hexl3sr.b2[1]/eps.hexl3sr.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.hexl3sr.b1[0]*eps.hexl3sr.b1[0]+j*j*eps.hexl3sr.b2[1]*eps.hexl3sr.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.hexl3sr.b1[0];
      Gif[npw+n]=j*eps.hexl3sr.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.recwg.es==0&&eps.recwg.er==0&&eps.recwg.r==0&&eps.recwg.ax==0&&eps.recwg.ly==0)){
  imax=(int)(Gmax/eps.recwg.b1[0]+1);
  jmax=(int)(imax*eps.recwg.b1[0]/eps.recwg.b2[1]+1);
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.recwg.b1[0]*eps.recwg.b1[0]+j*j*eps.recwg.b2[1]*eps.recwg.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.recwg.b1[0];
      Gif[npw+n]=j*eps.recwg.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.Thexwgr.es==0&&eps.Thexwgr.et==0&&eps.Thexwgr.L==0&&eps.Thexwgr.theta==0&&eps.Thexwgr.ly==0)){
  imax=(int)(Gmax/eps.Thexwgr.b1[0]+1);
  jmax=(int)(imax*eps.Thexwgr.b1[0]/eps.Thexwgr.b2[1]+1);
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.Thexwgr.b1[0]*eps.Thexwgr.b1[0]+j*j*eps.Thexwgr.b2[1]*eps.Thexwgr.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.Thexwgr.b1[0];
      Gif[npw+n]=j*eps.Thexwgr.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.Phexwgr.es==0&&eps.Phexwgr.ep==0&&eps.Phexwgr.l==0&&eps.Phexwgr.N==0&&eps.Phexwgr.th==0&&eps.Phexwgr.ly==0)){
  imax=(int)(Gmax/eps.Phexwgr.b1[0]+1);
  jmax=(int)(imax*eps.Phexwgr.b1[0]/eps.Phexwgr.b2[1]+1);
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.Phexwgr.b1[0]*eps.Phexwgr.b1[0]+j*j*eps.Phexwgr.b2[1]*eps.Phexwgr.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.Phexwgr.b1[0];
      Gif[npw+n]=j*eps.Phexwgr.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

 else if (!(eps.hexwgr.es==0&&eps.hexwgr.er==0&&eps.hexwgr.r==0&&eps.hexwgr.ly==0)){
  imax=(int)(Gmax/eps.hexwgr.b1[0]+1);
  jmax=(int)(imax*eps.hexwgr.b1[0]/eps.hexwgr.b2[1]+1);
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.hexwgr.b1[0]*eps.hexwgr.b1[0]+j*j*eps.hexwgr.b2[1]*eps.hexwgr.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.hexwgr.b1[0];
      Gif[npw+n]=j*eps.hexwgr.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 } 

 else if (!(eps.dishexl3r.es==0&&eps.dishexl3r.er==0&&eps.dishexl3r.r==0&&eps.dishexl3r.s==0&&eps.dishexl3r.rs==0&&eps.dishexl3r.lx==0&&eps.dishexl3r.ly==0&&eps.dishexl3r.sigma==0)){
  if (eps.dishexl3r.b1[0]>eps.dishexl3r.b2[1]){
   imax=(int)(Gmax/eps.dishexl3r.b1[0]+1);
   jmax=(int)(imax*eps.dishexl3r.b1[0]/eps.dishexl3r.b2[1]+1);
  }
  else{
   jmax=(int)(Gmax/eps.dishexl3r.b2[1]+1);
   imax=(int)(jmax*eps.dishexl3r.b2[1]/eps.dishexl3r.b1[0]+1);
  }
  for (int i=-imax;i<=imax;i++){
   for (int j=-jmax;j<=jmax;j++){
    nG=sqrt(i*i*eps.dishexl3r.b1[0]*eps.dishexl3r.b1[0]+j*j*eps.dishexl3r.b2[1]*eps.dishexl3r.b2[1]);
    if (nG<=Gmax){
     if(qd!=1){
      Gif[npw]=i*eps.dishexl3r.b1[0];
      Gif[npw+n]=j*eps.dishexl3r.b2[1];
      Gif[npw+2*n]=i;
      Gif[npw+3*n]=j;
      npw++;
     }
     else
      npw++;
    }
   }
  }
 }

}

/*######################################################################################################*/

/***************************************************************/
/*****************    Fourier coefficients    ******************/
/***************************************************************/

complex<double> Ceps (int i, int j, st_eps &eps){

 if (!(eps.recreg.es==0&&eps.recreg.er==0&&eps.recreg.r==0&&eps.recreg.lx==0&&eps.recreg.ly==0))
  return complex<double>(Crec(i,j,eps.recreg.es,eps.recreg.er,eps.recreg.r,eps.recreg.b1,eps.recreg.b2,eps.recreg.A),0);
 
 else if (!(eps.hexreg.es==0&&eps.hexreg.er==0&&eps.hexreg.r==0))
  return complex<double>(Chex(i,j,eps.hexreg.es,eps.hexreg.er,eps.hexreg.r,eps.hexreg.b1,eps.hexreg.b2,eps.hexreg.A),0);
 
 else if (!(eps.Thexreg.es==0&&eps.Thexreg.et==0&&eps.Thexreg.L==0&&eps.Thexreg.theta==0))
   return Thex(i,j,eps.Thexreg.es,eps.Thexreg.et,eps.Thexreg.L,eps.Thexreg.theta,eps.Thexreg.b1,eps.Thexreg.b2,eps.Thexreg.A); 

 else if (!(eps.rech1.es==0&&eps.rech1.er==0&&eps.rech1.r==0&&eps.rech1.lxp==0&&eps.rech1.lyp==0&&eps.rech1.lx==0&&eps.rech1.ly==0)){
  double cbl=0;
  if (fmod(i,eps.rech1.lx/eps.rech1.lxp)==0&&fmod(j,eps.rech1.ly/eps.rech1.lyp)==0)
   cbl=Crec((int)round(i*eps.rech1.lxp/eps.rech1.lx),(int)round(j*eps.rech1.lyp/eps.rech1.ly),eps.rech1.es,eps.rech1.er,eps.rech1.r,eps.rech1.bp1,eps.rech1.bp2,eps.rech1.Ap);
  return complex<double>(cbl+Crec(i,j,0,eps.rech1.es-eps.rech1.er,eps.rech1.r,eps.rech1.b1,eps.rech1.b2,eps.rech1.A),0);
 }

 else if (!(eps.hexh1.es==0&&eps.hexh1.er==0&&eps.hexh1.r==0&&eps.hexh1.l==0)){
  double cbl=0;
  if (fmod(i,eps.hexh1.l)==0&&fmod(j,eps.hexh1.l)==0)
   cbl=Chex(i/eps.hexh1.l,j/eps.hexh1.l,eps.hexh1.es,eps.hexh1.er,eps.hexh1.r,eps.hexh1.bp1,eps.hexh1.bp2,eps.hexh1.Ap);
  return complex<double>(cbl+Chex(i,j,0,eps.hexh1.es-eps.hexh1.er,eps.hexh1.r,eps.hexh1.b1,eps.hexh1.b2,eps.hexh1.A),0);
 }

 else if (!(eps.hexh1r.es==0&&eps.hexh1r.er==0&&eps.hexh1r.r==0&&eps.hexh1r.lx==0&&eps.hexh1r.ly==0)){
  double cbl=0; 
  if (fmod(i,eps.hexh1r.lx)==0&&fmod(j,round(eps.hexh1r.ly/sqrt(3)))==0)
   cbl=Chexr(i/eps.hexh1r.lx,(int)round(j*sqrt(3)/eps.hexh1r.ly),eps.hexh1r.es,eps.hexh1r.er,eps.hexh1r.r,eps.hexh1r.bp1,eps.hexh1r.bp2,eps.hexh1r.Ap);
  return complex<double>(cbl+Crec(i,j,0,eps.hexh1r.es-eps.hexh1r.er,eps.hexh1r.r,eps.hexh1r.b1,eps.hexh1r.b2,eps.hexh1r.A),0);
 }

 else if (!(eps.recl2.es==0&&eps.recl2.er==0&&eps.recl2.r==0&&eps.recl2.lxp==0&&eps.recl2.lyp==0&&eps.recl2.lx==0&&eps.recl2.ly==0)){
  double cbl=0;
  if (fmod(i,eps.recl2.lx/eps.recl2.lxp)==0&&fmod(j,eps.recl2.ly/eps.recl2.lyp)==0)
   cbl=Crec((int)round(i*eps.recl2.lxp/eps.recl2.lx),(int)round(j*eps.recl2.lyp/eps.recl2.ly),eps.recl2.es,eps.recl2.er,eps.recl2.r,eps.recl2.bp1,eps.recl2.bp2,eps.recl2.Ap);
  return (cbl+(exp(complex<double>(0,2*M_PI*i*eps.recl2.lxp/eps.recl2.lx))+complex<double>(1,0))*Crec(i,j,0,eps.recl2.es-eps.recl2.er,eps.recl2.r,eps.recl2.b1,eps.recl2.b2,eps.recl2.A))*exp(complex<double>(0,-M_PI*i*eps.recl2.lxp/eps.recl2.lx));
 }

 else if (!(eps.hexl2.es==0&&eps.hexl2.er==0&&eps.hexl2.r==0&&eps.hexl2.l==0)){
  double cbl=0;
  if (fmod(i,eps.hexl2.l)==0&&fmod(j,eps.hexl2.l)==0)
   cbl=Chex(i/eps.hexl2.l,j/eps.hexl2.l,eps.hexl2.es,eps.hexl2.er,eps.hexl2.r,eps.hexl2.bp1,eps.hexl2.bp2,eps.hexl2.Ap);
  return (cbl+(exp(complex<double>(0,2*M_PI*i/eps.hexl2.l))+complex<double>(1,0))*Chex(i,j,0,eps.hexl2.es-eps.hexl2.er,eps.hexl2.r,eps.hexl2.b1,eps.hexl2.b2,eps.hexl2.A))*exp(complex<double>(0,-M_PI*i/eps.hexl2.l));
 }

 else if (!(eps.hexl2r.es==0&&eps.hexl2r.er==0&&eps.hexl2r.r==0&&eps.hexl2r.lx==0&&eps.hexl2r.ly==0)){
  double cbl=0;
  if (fmod(i,eps.hexl2r.lx)==0&&fmod(j,round(eps.hexl2r.ly/sqrt(3)))==0)
   cbl=Chexr(i/eps.hexl2r.lx,(int)round(j*sqrt(3)/eps.hexl2r.ly),eps.hexl2r.es,eps.hexl2r.er,eps.hexl2r.r,eps.hexl2r.bp1,eps.hexl2r.bp2,eps.hexl2r.Ap);
  return (cbl+(exp(complex<double>(0,2*M_PI*i/eps.hexl2r.lx))+complex<double>(1,0))*Crec(i,j,0,eps.hexl2r.es-eps.hexl2r.er,eps.hexl2r.r,eps.hexl2r.b1,eps.hexl2r.b2,eps.hexl2r.A))*exp(complex<double>(0,-M_PI*i/eps.hexl2r.lx));  
 }

 else if (!(eps.recl3.es==0&&eps.recl3.er==0&&eps.recl3.r==0&&eps.recl3.s==0&&eps.recl3.rs==0&&eps.recl3.lxp==0&&eps.recl3.lyp==0&&eps.recl3.lx==0&&eps.recl3.ly==0)){
  double cbl=0;
  if (fmod(i,eps.recl3.lx/eps.recl3.lxp)==0&&fmod(j,eps.recl3.ly/eps.recl3.lyp)==0)
   cbl=Crec((int)round(i*eps.recl3.lxp/eps.recl3.lx),(int)round(j*eps.recl3.lyp/eps.recl3.ly),eps.recl3.es,eps.recl3.er,eps.recl3.r,eps.recl3.bp1,eps.recl3.bp2,eps.recl3.Ap);
  return complex<double>(cbl+(1+2*(cos(2*M_PI*i*eps.recl3.lxp/eps.recl3.lx)+cos(4*M_PI*i*eps.recl3.lxp/eps.recl3.lx)))*Crec(i,j,0,eps.recl3.es-eps.recl3.er,eps.recl3.r,eps.recl3.b1,eps.recl3.b2,eps.recl3.A)+2*cos(2*M_PI*i*(2+eps.recl3.s)*eps.recl3.lxp/eps.recl3.lx)*Crec(i,j,0,eps.recl3.er-eps.recl3.es,eps.recl3.rs,eps.recl3.b1,eps.recl3.b2,eps.recl3.A),0);
 }

 else if (!(eps.hexl3.es==0&&eps.hexl3.er==0&&eps.hexl3.r==0&&eps.hexl3.s==0&&eps.hexl3.rs==0&&eps.hexl3.l==0)){
  double cbl=0;
  if (fmod(i,eps.hexl3.l)==0&&fmod(j,eps.hexl3.l)==0)
   cbl=Chex(i/eps.hexl3.l,j/eps.hexl3.l,eps.hexl3.es,eps.hexl3.er,eps.hexl3.r,eps.hexl3.bp1,eps.hexl3.bp2,eps.hexl3.Ap);
  return complex<double>(cbl+(1+2*(cos(2*M_PI*i/eps.hexl3.l)+cos(4*M_PI*i/eps.hexl3.l)))*Chex(i,j,0,eps.hexl3.es-eps.hexl3.er,eps.hexl3.r,eps.hexl3.b1,eps.hexl3.b2,eps.hexl3.A)+2*cos(2*M_PI*i*(2+eps.hexl3.s)/eps.hexl3.l)*Chex(i,j,0,eps.hexl3.er-eps.hexl3.es,eps.hexl3.rs,eps.hexl3.b1,eps.hexl3.b2,eps.hexl3.A),0);
 }

 else if (!(eps.hexl3r.es==0&&eps.hexl3r.er==0&&eps.hexl3r.r==0&&eps.hexl3r.s==0&&eps.hexl3r.rs==0&&eps.hexl3r.lx==0&&eps.hexl3r.ly==0)){
  double cbl=0;
  if (fmod(i,eps.hexl3r.lx)==0&&fmod(j,round(eps.hexl3r.ly/sqrt(3)))==0)
   cbl=Chexr(i/eps.hexl3r.lx,(int)round(j*sqrt(3)/eps.hexl3r.ly),eps.hexl3r.es,eps.hexl3r.er,eps.hexl3r.r,eps.hexl3r.bp1,eps.hexl3r.bp2,eps.hexl3r.Ap);
  return complex<double>(cbl+(1+2*(cos(2*M_PI*i/eps.hexl3r.lx)+cos(4*M_PI*i/eps.hexl3r.lx)))*Crec(i,j,0,eps.hexl3r.es-eps.hexl3r.er,eps.hexl3r.r,eps.hexl3r.b1,eps.hexl3r.b2,eps.hexl3r.A)+2*cos(2*M_PI*i*(2+eps.hexl3r.s)/eps.hexl3r.lx)*Crec(i,j,0,eps.hexl3r.er-eps.hexl3r.es,eps.hexl3r.rs,eps.hexl3r.b1,eps.hexl3r.b2,eps.hexl3r.A),0);
 }

 else if (!(eps.hexl3sr.es==0&&eps.hexl3sr.er==0&&eps.hexl3sr.r==0&&eps.hexl3sr.s1==0&&eps.hexl3sr.rs1==0&&eps.hexl3sr.s2==0&&eps.hexl3sr.rs2==0&&eps.hexl3sr.s3==0&&eps.hexl3sr.rs3==0&&eps.hexl3sr.s4==0&&eps.hexl3sr.rs4==0&&eps.hexl3sr.s5==0&&eps.hexl3sr.rs5==0&&eps.hexl3sr.lx==0&&eps.hexl3sr.ly==0)){
  double cbl=0;
  if (fmod(i,eps.hexl3sr.lx)==0&&fmod(j,round(eps.hexl3sr.ly/sqrt(3)))==0)
   cbl=Chexr(i/eps.hexl3sr.lx,(int)round(j*sqrt(3)/eps.hexl3sr.ly),eps.hexl3sr.es,eps.hexl3sr.er,eps.hexl3sr.r,eps.hexl3sr.bp1,eps.hexl3sr.bp2,eps.hexl3sr.Ap);
  return complex<double>(cbl+(1+2*(cos(2*M_PI*i/eps.hexl3sr.lx)+cos(4*M_PI*i/eps.hexl3sr.lx)+cos(6*M_PI*i/eps.hexl3sr.lx)+cos(8*M_PI*i/eps.hexl3sr.lx)+cos(10*M_PI*i/eps.hexl3sr.lx)+cos(12*M_PI*i/eps.hexl3sr.lx)))*Crec(i,j,0,eps.hexl3sr.es-eps.hexl3sr.er,eps.hexl3sr.r,eps.hexl3sr.b1,eps.hexl3sr.b2,eps.hexl3sr.A)+2*cos(2*M_PI*i*(2+eps.hexl3sr.s1)/eps.hexl3sr.lx)*Crec(i,j,0,eps.hexl3sr.er-eps.hexl3sr.es,eps.hexl3sr.rs1,eps.hexl3sr.b1,eps.hexl3sr.b2,eps.hexl3sr.A)+2*cos(2*M_PI*i*(3+eps.hexl3sr.s2)/eps.hexl3sr.lx)*Crec(i,j,0,eps.hexl3sr.er-eps.hexl3sr.es,eps.hexl3sr.rs2,eps.hexl3sr.b1,eps.hexl3sr.b2,eps.hexl3sr.A)+2*cos(2*M_PI*i*(4+eps.hexl3sr.s3)/eps.hexl3sr.lx)*Crec(i,j,0,eps.hexl3sr.er-eps.hexl3sr.es,eps.hexl3sr.rs3,eps.hexl3sr.b1,eps.hexl3sr.b2,eps.hexl3sr.A)+2*cos(2*M_PI*i*(5+eps.hexl3sr.s4)/eps.hexl3sr.lx)*Crec(i,j,0,eps.hexl3sr.er-eps.hexl3sr.es,eps.hexl3sr.rs4,eps.hexl3sr.b1,eps.hexl3sr.b2,eps.hexl3sr.A)+2*cos(2*M_PI*i*(6+eps.hexl3sr.s5)/eps.hexl3sr.lx)*Crec(i,j,0,eps.hexl3sr.er-eps.hexl3sr.es,eps.hexl3sr.rs5,eps.hexl3sr.b1,eps.hexl3sr.b2,eps.hexl3sr.A),0);
 }

 else if (!(eps.recwg.es==0&&eps.recwg.er==0&&eps.recwg.r==0&&eps.recwg.lyp==0&&eps.recwg.ax==0&&eps.recwg.ly==0)){
  double cbl=0;
  if (fmod(j,eps.recwg.ly/eps.recwg.lyp)==0)
   cbl=Crec(i,(int)round(j*eps.recwg.lyp/eps.recwg.ly),eps.recwg.es,eps.recwg.er,eps.recwg.r,eps.recwg.bp1,eps.recwg.bp2,eps.recwg.Ap);
  return complex<double>(cbl+Crec(i,j,0,eps.recwg.es-eps.recwg.er,eps.recwg.r,eps.recwg.b1,eps.recwg.b2,eps.recwg.A),0);
 }

 else if (!(eps.Thexwgr.es==0&&eps.Thexwgr.et==0&&eps.Thexwgr.L==0&&eps.Thexwgr.theta==0&&eps.Thexwgr.ly==0)){
  complex<double> cl(0,0);
  
  if(i==0&&j==0)
   cl=eps.Thexwgr.es;
  for (int p=0;p<eps.Thexwgr.Nh;p++){
   cl=cl+Tri(i,j,0.0,eps.Thexwgr.Mch[p+2*eps.Thexwgr.Nh],eps.Thexwgr.L,eps.Thexwgr.theta,eps.Thexwgr.b1,eps.Thexwgr.b2,eps.Thexwgr.A)*cexp(-2*M_PI*(eps.Thexwgr.Mch[p]*i+eps.Thexwgr.Mch[p+eps.Thexwgr.Nh]*j/eps.Thexwgr.ly));
  }
  return cl;
 }

else if (!(eps.Phexwgr.es==0&&eps.Phexwgr.ep==0&&eps.Phexwgr.l==0&&eps.Phexwgr.N==0&&eps.Phexwgr.th==0&&eps.Phexwgr.ly==0)){
  complex<double> cl(0,0);
  
  if(i==0&&j==0)
   cl=eps.Phexwgr.es;
  for (int p=0;p<eps.Phexwgr.Nh;p++){
   cl=cl+Polgme(i,j,0.0,eps.Phexwgr.Mch[p+2*eps.Phexwgr.Nh],eps.Phexwgr.l,eps.Phexwgr.N,eps.Phexwgr.th,eps.Phexwgr.b1,eps.Phexwgr.b2,eps.Phexwgr.A)*cexp(-2*M_PI*(eps.Phexwgr.Mch[p]*i+eps.Phexwgr.Mch[p+eps.Phexwgr.Nh]*j/eps.Phexwgr.ly));
  }
  return cl;
 }

 else if (!(eps.hexwgr.es==0&&eps.hexwgr.er==0&&eps.hexwgr.r==0&&eps.hexwgr.ly==0)){
  double cbl=0;
  if (fmod(j,round(eps.hexwgr.ly/sqrt(3)))==0)
   cbl=Chexr(i,(int)round(j*sqrt(3)/eps.hexwgr.ly),eps.hexwgr.es,eps.hexwgr.er,eps.hexwgr.r,eps.hexwgr.bp1,eps.hexwgr.bp2,eps.hexwgr.Ap);
  return complex<double>(cbl+Crec(i,j,0,eps.hexwgr.es-eps.hexwgr.er,eps.hexwgr.r,eps.hexwgr.b1,eps.hexwgr.b2,eps.hexwgr.A),0);
 } 

 else if (!(eps.dishexl3r.es==0&&eps.dishexl3r.er==0&&eps.dishexl3r.r==0&&eps.dishexl3r.s==0&&eps.dishexl3r.rs==0&&eps.dishexl3r.lx==0&&eps.dishexl3r.ly==0&&eps.dishexl3r.sigma==0)){
  complex<double> cl(0,0);

  if (i==0 && j==0)
   cl=eps.dishexl3r.es*(1-eps.dishexl3r.Nh);
  for (int p=0;p<eps.dishexl3r.Nh;p++){
   cl=cl+Crec(i,j,eps.dishexl3r.es,eps.dishexl3r.Mch[p+3*eps.dishexl3r.Nh],eps.dishexl3r.Mch[p+2*eps.dishexl3r.Nh],eps.dishexl3r.b1,eps.dishexl3r.b2,eps.dishexl3r.A)*exp(complex<double>(0,-2*M_PI*(eps.dishexl3r.Mch[p]*i/eps.dishexl3r.lx+eps.dishexl3r.Mch[p+eps.dishexl3r.Nh]*j/eps.dishexl3r.ly)));
  }
  return cl;
 }
 
}

/*######################################################################################################*/

/***************************************************************/
/*******************    Dielectric matrix    *******************/
/***************************************************************/

void dseps (st_eps &eps, double *G, double *meps, int n){ // real case
 t=clock(); 
 for (int j=0;j<n;j++){
   for (int i=0;i<=j;i++){
     meps[iju(i,j)]=real(Ceps(G[i+2*n]-G[j+2*n],G[i+3*n]-G[j+3*n],eps));
   }
 }
 t=(clock()-t)/CLOCKS_PER_SEC;
 cout<<"made in "<<t<<" seconds, ";
}

void zheps (st_eps &eps, double *G, complex<double> *meps, int n){ // complex case
 t=clock(); 
 for (int j=0;j<n;j++){
   for (int i=0;i<=j;i++){
     meps[iju(i,j)]=Ceps(G[i+2*n]-G[j+2*n],G[i+3*n]-G[j+3*n],eps);
   }
 }
 t=(clock()-t)/CLOCKS_PER_SEC;
 cout<<"made in "<<t<<" seconds, ";
}

/*######################################################################################################*/

/***************************************************************/
/***************    Writing fourier expansion    ***************/
/***************************************************************/

struct eps_box_lim {  
 double xi; // initial x value
 double xf; // final x value
 double dx; // step in x
 double yi; // initial y value
 double yf; // final y value
 double dy; // step in y
 double Gmax; // cutoff in the Fourier expansion
 double weps; //if equal to 1 writes the Fourier expansion of epsilon
}epslim;

void write_eps (st_eps &eps, eps_box_lim &epslim){
 if (epslim.weps==1){
  int npw;
  double *G;
  complex<double> *VCeps,feps;

  cout<<"Writing epsilon: ";
  cout.flush();
  t=clock(); 

  ofstream epsout ("epsilon.dat");
  epsout.precision(set_pre);

  Gi(eps,G,1,npw,epslim.Gmax,npw);
  G=new double[4*npw];
  Gi(eps,G,0,npw,epslim.Gmax,npw);
  VCeps=new complex<double>[npw];

  for (int i=0;i<npw;i++)
   VCeps[i]=Ceps(G[i+2*npw],G[i+3*npw],eps);

   for (double x=epslim.xi;x<=epslim.xf;x=x+epslim.dx){
    for (double y=epslim.yi;y<=epslim.yf;y=y+epslim.dy){
     feps=complex<double>(0,0); 
      for (int i=0;i<npw;i++){
       feps=feps+VCeps[i]*exp(complex<double>(0,G[i]*x+G[i+npw]*y));
      }
      epsout<<x<<"	"<<y<<"	"<<real(feps)<<"\n";
    }
   }
  epsout.close();

  delete[] G;
  delete[] VCeps;

  t=(clock()-t)/CLOCKS_PER_SEC;
  cout<<"written in "<<t<<" seconds\n\n\n"; 
  cout.flush();
 }
}


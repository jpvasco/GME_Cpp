/*******************************************************************************/
/**************    Matrix and vector operations calling lapack    **************/
/*******************************************************************************/

/* 
 * THIS LIBRARY CONSIDERS THE MATRICES IN A COLUMN-MAJOR ORDERED CONVENTION
 */


/* It is important to take into account that the lapack subroutines can be
 * used only if the lapack library is correctly linked in the compilation
 * process. To use all the lapack features the libraries are linked in the 
 * following form:
   
 *  g++ $$other-options$$ -L $library-location$ -llapack -lblas -lgfortran -lm
 
 * when $library-location$ is the default location it is not necessary to
 * put the -L option.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex>

//#include "matrix-cula.h" /* CUDA compatibility */

using namespace std;

/* fortran packed storage for the
 * upper triangle of the matrix
 * 0<=i<=j
 */

int iju (int i, int j){
  return ((i+1)+j*(j+1)/2)-1;
}


/******************************************************************/
/*********************    Matrix operations   *********************/
/******************************************************************/

/*######################################################################################################*/

/**********************************************************/
/** Addition of two matrixes in vector form y=y+alpha*x ***/
/**********************************************************/

/* DAXPY subroutine
 * DAXPY constant times a vector plus a vector.
 * uses unrolled loops for increments equal to one.
 */

extern "C" {
     void daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
}

/* ZAXPY subroutine
 * ZAXPY constant times a vector plus a vector.
 * uses unrolled loops for increments equal to one.
 */

extern "C" {
     void zaxpy_(int *n, complex<double> *da, complex<double> *dx, int *incx, complex<double> *dy, int *incy);
}


/* Matrix sum a=a+b;
 * real case
 */

void dApB (double *a, double *b, int n){
 int incx=1,incy=1,nn=n*n;
 double alpha=1.0;
 daxpy_(&nn,&alpha,b,&incx,a,&incy);
}

/* Matrix subtraction a=a-b;
 * real case
 */

void dAmB (double *a, double *b, int n){
 int incx=1,incy=1,nn=n*n;
 double alpha=-1.0;
 daxpy_(&nn,&alpha,b,&incx,a,&incy);
}

/* Matrix addition a=a+alpha*b;
 * real case
 */

void dApaB (double *a, double alpha, double *b, int n){
 int incx=1,incy=1,nn=n*n;
 daxpy_(&nn,&alpha,b,&incx,a,&incy);
}

/* Matrix sum a=a+b;
 * complex case
 */

void zApB (complex<double> *a, complex<double> *b, int n){
 int incx=1,incy=1,nn=n*n;
 complex<double> alpha(1.0,0);
 zaxpy_(&nn,&alpha,b,&incx,a,&incy);
}

/* Matrix subtraction a=a-b;
 * complex case
 */

void zAmB (complex<double> *a, complex<double> *b, int n){
 int incx=1,incy=1,nn=n*n;
 complex<double> alpha(-1.0,0);
 zaxpy_(&nn,&alpha,b,&incx,a,&incy);
}

/* Matrix addition a=a+alpha*b;
 * complex case
 */

void zApaB (complex<double> *a, complex<double> alpha, complex<double> *b, int n){
 int incx=1,incy=1,nn=n*n;
 zaxpy_(&nn,&alpha,b,&incx,a,&incy);
}

/*######################################################################################################*/

/**********************************************************/
/***************** Matrix times a constant ****************/
/**********************************************************/

/* DSCAL subroutine
 * DSCAL scales a vector by a constant.
 * uses unrolled loops for increment equal to one.
 */

extern "C" {
     void dscal_(int *n, double *da, double *dx, int *incx);
}

/* ZSCAL subroutine
 * ZSCAL scales a vector by a constant.
 */

extern "C" {
     void zscal_(int *n, complex<double> *da, complex<double> *dx, int *incx);
}

/* Matrix times a constant c*A   
 * real case
 */

void dcA (double c, double *a, int n){
 int incx=1,nn=n*n;
 dscal_(&nn,&c,a,&incx);
}

/* Matrix times a constant c*A   
 * complex case
 */

void zcA (complex<double> c, complex<double> *a, int n){
 int incx=1,nn=n*n;
 zscal_(&nn,&c,a,&incx);
}


/* Diagonal-Matrix times a constant c*A   
 * real case
 */

void dcDA (double c, double *a, int n){
 for (int i=0;i<n;i++){
  a[i+i*n]=c*a[i+i*n];
 }
}

void dcDAp (double c, double *a, int n){ /* packed storage */
 for (int i=0;i<n;i++){
  a[iju(i,i)]=c*a[iju(i,i)];
 }
}

/* Diagonal-Matrix times a constant c*A   
 * complex case
 */

void zcDA (complex<double> c, complex<double> *a, int n){
 for (int i=0;i<n;i++){
  a[i+i*n]=c*a[i+i*n];
 }
}

void zcDAp (complex<double> c, complex<double> *a, int n){ /* packed storage */
 for (int i=0;i<n;i++){
  a[iju(i,i)]=c*a[iju(i,i)];
 }
}

/*######################################################################################################*/

/**********************************************************/
/****************** Matrix-Matrix product *****************/
/**********************************************************/

/* DGEMM subroutine
 * DGEMM  performs one of the matrix-matrix operations
 *
 *    C := alpha*op( A )*op( B ) + beta*C,
 *
 * where  op( X ) is one of
 *
 *    op( X ) = X   or   op( X ) = X**T,
 *
 * alpha and beta are scalars, and A, B and C are matrices, with op( A )
 * an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 */

extern "C" {
     void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);      
}

/* ZGEMM subroutine
 * ZGEMM  performs one of the matrix-matrix operations
 *
 *    C := alpha*op( A )*op( B ) + beta*C,
 *
 * where  op( X ) is one of
 *
 *    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
 *
 * alpha and beta are scalars, and A, B and C are matrices, with op( A )
 * an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 */

extern "C" {
     void zgemm_(char *transa, char *transb, int *m, int *n, int *k, complex<double> *alpha, complex<double> *a, int *lda, complex<double> *b, int *ldb, complex<double> *beta, complex<double> *c, int *ldc);      
}


/* Square Matrix-Matrix product   
 * real case
 */

void dAxB (double *a, double *b, double *c, int n){
 char transa='N',transb='N';
 double alpha=1.0,beta=0.0;

 dgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}

/* Square Matrix-Matrix product   
 * complex case
 */

void zAxB (complex<double> *a, complex<double> *b, complex<double> *c, int n){
 char transa='N',transb='N';
 complex<double> alpha(1.0,0),beta(0.0,0.0);

 zgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}

/* Square alpha-Matrix-Matrix product   
 * real case
 */

void daAxB (double alpha, double *a, double *b, double *c, int n){
 char transa='N',transb='N';
 double beta=0.0;

 dgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}

/* Square alpha-Matrix-Matrix product   
 * complex case
 */

void zaAxB (complex<double> alpha, complex<double> *a, complex<double> *b, complex<double> *c, int n){
 char transa='N',transb='N';
 complex<double> beta(0.0,0.0);

 zgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}

/* Square alpha-Matrix-Matrix product plus beta-Matrix  
 * real case
 */

void daAxBpbC (double alpha, double *a, double *b, double beta, double *c, int n){
 char transa='N',transb='N';

 dgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}

/* Square alpha-Matrix-Matrix product plus beta-Matrix  
 * complex case
 */

void zaAxBpbC (complex<double> alpha, complex<double> *a, complex<double> *b, complex<double> beta, complex<double> *c, int n){
 char transa='N',transb='N';

 zgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}

/*######################################################################################################*/

/**********************************************************/
/*********************** Utilities ************************/
/**********************************************************/

/* check if the square matrix M is symmetric */

int checksM (double *M, int n){
 int res=1;
 double tol=0.000000000001,elemij,elemji;
 
 for(int i=0;i<n;i++){
  for(int j=0;j<n;j++){
   elemij=M[i+j*n];        
   elemji=M[j+i*n];
   if (fabs(elemij-elemji)>tol){
    res=0;
    break;
   }
  }
 }
 return res;
}

/* check if the square matrix M is hermitian */

int checkhM (complex<double> *M, int n){
 int res=1;
 double tol=0.000000000001;
 complex<double> elemij,elemji;
 
 for(int i=0;i<n;i++){
  for(int j=0;j<n;j++){
   elemij=M[i+j*n];        
   elemji=M[j+i*n];
   if (fabs(real(elemij)-real(elemji))>tol || fabs(imag(elemij)-imag(conj(elemji)))>tol){
    res=0;
    break;
   }
  }
 }
 return res;
}

/* makes a square matrix M with the square block matrices a, b, c and d.
   The M matrix will be {{a,b},{c,d}}
*/

/* real case */

void dblockM (double *a, double *b, double *c, double *d, int n, double *M){
 
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[i+j*2*n]=a[i+j*n];
  }
 }
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[i+(j+n)*2*n]=b[i+j*n];
  }
 }
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[(i+n)+j*2*n]=c[i+j*n];
  }
 }
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[(i+n)+(j+n)*2*n]=d[i+j*n];
  }
 } 
}

/* complex case */

void zblockM (complex<double> *a, complex<double> *b, complex<double> *c, complex<double> *d, int n, complex<double> *M){
 
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[i+j*2*n]=a[i+j*n];
  }
 }
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[i+(j+n)*2*n]=b[i+j*n];
  }
 }
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[(i+n)+j*2*n]=c[i+j*n];
  }
 }
 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   M[(i+n)+(j+n)*2*n]=d[i+j*n];
  }
 } 
}

/* build the square matrix a from the upper packed storage matrix ap
 * ap must be symmetric or hermitian
 */

/* real case */

void dsap2a (double *ap, double *a, int n){

 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   if (j>=i)
    a[i+j*n]=ap[iju(i,j)];
   else
    a[i+j*n]=ap[iju(j,i)];
  }
 }
}

/* complex hermitian case */

void zhap2a (complex<double> *ap, complex<double> *a, int n){

 for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
   if (j>=i)
    a[i+j*n]=ap[iju(i,j)];
   else
    a[i+j*n]=conj(ap[iju(j,i)]);
  }
 }
}

/* copy the array a in ac
 */
 
/* real case */

void dacp (double *a, double *ac, int la){

 for (int i=0;i<la;i++)
  ac[i]=a[i];
}

/* complex case */

void zacp (complex<double> *a, complex<double> *ac, int la){

 for (int i=0;i<la;i++)
  ac[i]=a[i];
}

/* Diagonal-Matrix with diagonal elements c  
 * real case
 */

void dDc (double c, double *a, int n){
 for (int i=0;i<n;i++){
  a[i+i*n]=c;
 }
}

void dDcp (double c, double *a, int n){ /* packed storage */
 for (int i=0;i<n;i++){
  a[iju(i,i)]=c;
 }
}

/* Diagonal-Matrix with diagonal elements c  
 * complex case
 */

void zDc (complex<double> c, complex<double> *a, int n){
 for (int i=0;i<n;i++){
  a[i+i*n]=c;
 }
}

void zDcp (complex<double> c, complex<double> *a, int n){ /* packed storage */
 for (int i=0;i<n;i++){
  a[iju(i,i)]=c;
 }
}

/* Copy column jth from the square matrix a in vj 
 * jth run from 0 to (n-1)
 */

/* real case */

void dgetvj (double *a, double *vj, int jth, int n){
 for (int i=0;i<n;i++){
  vj[i]=a[i+jth*n];
 }
}

/* complex case */

void zgetvj (complex<double> *a, complex<double> *vj, int jth, int n){
 for (int i=0;i<n;i++){
  vj[i]=a[i+jth*n];
 }
}

/* Copy row ith from the square matrix a in vi 
 * ith run from 0 to (n-1)
 */

/* real case */

void dgetvi (double *a, double *vi, int ith, int n){
 for (int j=0;j<n;j++){
  vi[j]=a[ith+j*n];
 }
}

/* complex case */

void zgetvi (complex<double> *a, complex<double> *vi, int ith, int n){
 for (int j=0;j<n;j++){
  vi[j]=a[ith+j*n];
 }
}

/* Copy column jth from the non-square matrix a in vj 
 * jth run from 0 to (nj-1). 
 * ni: number of rows, nj: number of columns
 */

/* real case */

void dnsgetvj (double *a, double *vj, int jth, int ni){
 for (int i=0;i<ni;i++){
  vj[i]=a[i+jth*ni];
 }
}

/* complex case */

void znsgetvj (complex<double> *a, complex<double> *vj, int jth, int ni){
 for (int i=0;i<ni;i++){
  vj[i]=a[i+jth*ni];
 }
}

/* Copy row ith from the non-square matrix a in vi 
 * ith run from 0 to (ni-1)
 * ni: number of rows, nj: number of columns
 */

/* real case */

void dnsgetvi (double *a, double *vi, int ith, int ni, int nj){
 for (int j=0;j<nj;j++){
  vi[j]=a[ith+j*ni];
 }
}

/* complex case */

void znsgetvi (complex<double> *a, complex<double> *vi, int ith, int ni, int nj){
 for (int j=0;j<nj;j++){
  vi[j]=a[ith+j*ni];
 }
}

/*######################################################################################################*/

/**********************************************************/
/********** Inverse of a real symmetric matrix ************/
/**********************************************************/
/* DSPTRF and DSPTRI subroutines 
 *
 * DSPTRF computes the factorization of a real symmetric matrix A stored
 * in packed format using the Bunch-Kaufman diagonal pivoting method:
 *
 *    A = U*D*U**T  or  A = L*D*L**T
 *
 * where U (or L) is a product of permutation and unit upper (lower)
 * triangular matrices, and D is symmetric and block diagonal with
 * 1-by-1 and 2-by-2 diagonal blocks.
 *
 * DSPTRI computes the inverse of a real symmetric indefinite matrix
 * A in packed storage using the factorization A = U*D*U**T or
 * A = L*D*L**T computed by DSPTRF.
 */

extern "C" {
	void dsptrf_(char *uplo, int *n, double *ap, int *ipiv, int *info);
}

extern "C" {
	void dsptri_(char *uplo, int *n, double *ap, int *ipiv, double *work, int *info);
}


void dspinv (double *ap, int n){ /* packed storage*/

 //lapack subroutine parameters
 char uplo='U';
 int *ipiv,info;
 double *work;

 ipiv=new int[n];
 work=new double[n];

 t=clock();
 dsptrf_(&uplo,&n,ap,ipiv,&info);
 dsptri_(&uplo,&n,ap,ipiv,work,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"inverted in "<<t<<" seconds, "; 
 else{
  cout<<"An error has occurred inverting the matrix\n";
  abort();
 }

 delete[] ipiv;
 delete[] work;
}

/**********************************************************/
/************* Inverse of a hermitian matrix **************/
/**********************************************************/
/* ZHPTRF and ZHPTRI subroutines 
 * 
 * ZHPTRF computes the factorization of a complex Hermitian packed
 * matrix A using the Bunch-Kaufman diagonal pivoting method:
 * 
 *    A = U*D*U**H  or  A = L*D*L**H
 * 
 * where U (or L) is a product of permutation and unit upper (lower)
 * triangular matrices, and D is Hermitian and block diagonal with
 * 1-by-1 and 2-by-2 diagonal blocks.
 *
 * ZHPTRI computes the inverse of a complex Hermitian indefinite matrix
 * A in packed storage using the factorization A = U*D*U**H or
 * A = L*D*L**H computed by ZHPTRF.
 */

extern "C" {
	void zhptrf_(char *uplo, int *n, complex<double> *ap, int *ipiv, int *info);
}

extern "C" {
	void zhptri_(char *uplo, int *n, complex<double> *ap, int *ipiv, complex<double> *work, int *info);
}


void zhpinv (complex<double> *ap, int n){ /* packed storage*/

 //lapack subroutine parameters
 char uplo='U';
 int *ipiv,info;
 complex<double> *work;

 ipiv=new int[n];
 work=new complex<double>[n];
 
 t=clock();
 zhptrf_(&uplo,&n,ap,ipiv,&info);
 zhptri_(&uplo,&n,ap,ipiv,work,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"inverted in "<<t<<" seconds, "; 
 else{
  cout<<"An error has occurred inverting the matrix\n";
  abort();
 }

 delete[] ipiv;
 delete[] work;
}

/**********************************************************/
/*********** Inverse of a real general matrix *************/
/**********************************************************/
/* DGETRF and DGETRI subroutines
 *
 * DGETRF computes an LU factorization of a general M-by-N matrix A
 * using partial pivoting with row interchanges.
 *
 * The factorization has the form
 *    A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapezoidal if m > n), and U is upper
 * triangular (upper trapezoidal if m < n).
 *
 * This is the right-looking Level 3 BLAS version of the algorithm.
 *
 * DGETRI computes the inverse of a matrix using the LU factorization
 * computed by DGETRF.
 *
 * This method inverts U and then computes inv(A) by solving the system
 * inv(A)*L = inv(U) for inv(A).
 */

extern "C" {
	void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
}

extern "C" {
	void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
}

void dgeinv (double *a, int n){

 int *ipiv,info;
 double *work;

 ipiv=new int[n];
 work=new double[n];

 t=clock();
 dgetrf_(&n,&n,a,&n,ipiv,&info);
 dgetri_(&n,a,&n,ipiv,work,&n,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"inverted in "<<t<<" seconds, "; 
 else{
  cout<<"An error has occurred inverting the matrix\n";
  abort();
 }

 delete[] ipiv;
 delete[] work;
}

/**********************************************************/
/********* Inverse of a complex general matrix ************/
/**********************************************************/
/* ZGETRF and ZGETRI subroutines
 *
 * ZGETRF computes an LU factorization of a general M-by-N matrix A
 * using partial pivoting with row interchanges.
 *
 * The factorization has the form
 *    A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapezoidal if m > n), and U is upper
 * triangular (upper trapezoidal if m < n).
 *
 * This is the Crout Level 3 BLAS version of the algorithm.
 *
 * ZGETRI computes the inverse of a matrix using the LU factorization
 * computed by ZGETRF.
 *
 * This method inverts U and then computes inv(A) by solving the system
 * inv(A)*L = inv(U) for inv(A).
 */

extern "C" {
	void zgetrf_(int *m, int *n, complex<double> *a, int *lda, int *ipiv, int *info);
}

extern "C" {
	void zgetri_(int *n, complex<double> *a, int *lda, int *ipiv, complex<double> *work, int *lwork, int *info);
}

void zgeinv (complex<double> *a, int n){

 int *ipiv,info;
 complex<double> *work;

 ipiv=new int[n];
 work=new complex<double>[n];

 t=clock();
 zgetrf_(&n,&n,a,&n,ipiv,&info);
 zgetri_(&n,a,&n,ipiv,work,&n,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"inverted in "<<t<<" seconds, "; 
 else{
  cout<<"An error has occurred inverting the matrix\n";
  abort();
 }

 delete[] ipiv;
 delete[] work;
}

/*######################################################################################################*/

/**********************************************************/
/*** Eigenvalues and eigenvectors of a symmetric matrix ***/
/**********************************************************/
/* DSYEVR subroutine
 *
 * DSYEVR computes selected eigenvalues and, optionally, eigenvectors
 * of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
 * selected by specifying either a range of values or a range of
 * indices for the desired eigenvalues.
 *
 * DSYEVR first reduces the matrix A to tridiagonal form T with a call
 * to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute
 * the eigenspectrum using Relatively Robust Representations.  DSTEMR
 * computes eigenvalues by the dqds algorithm, while orthogonal
 * eigenvectors are computed from various "good" L D L^T representations
 * (also known as Relatively Robust Representations). Gram-Schmidt
 * orthogonalization is avoided as far as possible. More specifically,
 * the various steps of the algorithm are as follows.
 *
 * For each unreduced block (submatrix) of T,
 *    (a) Compute T - sigma I  = L D L^T, so that L and D
 *        define all the wanted eigenvalues to high relative accuracy.
 *        This means that small relative changes in the entries of D and L
 *        cause only small relative changes in the eigenvalues and
 *        eigenvectors. The standard (unfactored) representation of the
 *        tridiagonal matrix T does not have this property in general.
 *    (b) Compute the eigenvalues to suitable accuracy.
 *        If the eigenvectors are desired, the algorithm attains full
 *        accuracy of the computed eigenvalues only right before
 *        the corresponding vectors have to be computed, see steps c) and d).
 *    (c) For each cluster of close eigenvalues, select a new
 *        shift close to the cluster, find a new factorization, and refine
 *        the shifted eigenvalues to suitable accuracy.
 *    (d) For each eigenvalue with a large enough relative separation compute
 *        the corresponding eigenvector by forming a rank revealing twisted
 *        factorization. Go back to (c) for any clusters that remain.
 *
 * The desired accuracy of the output can be specified by the input
 * parameter ABSTOL.
 *
 * For more details, see DSTEMR's documentation and:
 * - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
 *   to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
 *   Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
 * - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
 *   Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
 *   2004.  Also LAPACK Working Note 154.
 * - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
 *   tridiagonal eigenvalue/eigenvector problem",
 *   Computer Science Division Technical Report No. UCB/CSD-97-971,
 *   UC Berkeley, May 1997.
 *
 *
 * Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested
 * on machines which conform to the ieee-754 floating point standard.
 * DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and
 * when partial spectrum requests are made.
 *
 * Normal execution of DSTEMR may create NaNs and infinities and
 * hence may abort due to a floating point exception in environments
 * which do not handle NaNs and infinities in the ieee standard default
 * manner.
 */

extern "C" {
     void dsyevr_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);
}


/* eigenvalues only */

void dseigen (double *a, double *w, int n, int stat_i, int stat_f){

 //lapack subroutine parameters
 char jobz='N',range='I',uplo='U';
 int lda=n,il=stat_i+1,iu=stat_f+1,m=iu-il+1,ldz=1,*isuppz,lwork=26*n,liwork=10*n,*iwork,info;
 double vl=0,vu=1,abstol=0;
 double *z,*work;

 iwork=new int[liwork];
 work=new double[lwork];
 
 t=clock();
 dsyevr_(&jobz,&range,&uplo,&n,a,&lda,&vl,&vu,&il,&iu,&abstol,&m,w,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }

 delete[] iwork;
 delete[] work;
}

/* eigenvalues and eigenvectors */

void dseigenv (double *a, double *w, double *vw, int n, int stat_i, int stat_f){

 //lapack subroutine parameters
 char jobz='V',range='I',uplo='U';
 int lda=n,il=stat_i+1,iu=stat_f+1,m=iu-il+1,ldz=n,*isuppz,lwork=26*n,liwork=10*n,*iwork,info;
 double vl=0,vu=1,abstol=0;
 double *work;

 isuppz=new int[2*m];
 iwork=new int[liwork];
 work=new double[lwork];
 
 t=clock();
 dsyevr_(&jobz,&range,&uplo,&n,a,&lda,&vl,&vu,&il,&iu,&abstol,&m,w,vw,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }

 delete[] isuppz;  
 delete[] iwork;
 delete[] work;
} 

/**********************************************************/
/*** Eigenvalues and eigenvectors of a hermitian matrix ***/
/**********************************************************/
/* ZHEEVR subroutine
 *
 * ZHEEVR computes selected eigenvalues and, optionally, eigenvectors
 * of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can
 * be selected by specifying either a range of values or a range of
 * indices for the desired eigenvalues.
 * 
 * ZHEEVR first reduces the matrix A to tridiagonal form T with a call
 * to ZHETRD.  Then, whenever possible, ZHEEVR calls ZSTEMR to compute
 * eigenspectrum using Relatively Robust Representations.  ZSTEMR
 * computes eigenvalues by the dqds algorithm, while orthogonal
 * eigenvectors are computed from various "good" L D L^T representations
 * (also known as Relatively Robust Representations). Gram-Schmidt
 * orthogonalization is avoided as far as possible. More specifically,
 * the various steps of the algorithm are as follows.
 * 
 * For each unreduced block (submatrix) of T,
 *    (a) Compute T - sigma I  = L D L^T, so that L and D
 *        define all the wanted eigenvalues to high relative accuracy.
 *        This means that small relative changes in the entries of D and L
 *       cause only small relative changes in the eigenvalues and
 *        eigenvectors. The standard (unfactored) representation of the
 *        tridiagonal matrix T does not have this property in general.
 *    (b) Compute the eigenvalues to suitable accuracy.
 *        If the eigenvectors are desired, the algorithm attains full
 *        accuracy of the computed eigenvalues only right before
 *        the corresponding vectors have to be computed, see steps c) and d).
 *    (c) For each cluster of close eigenvalues, select a new
 *        shift close to the cluster, find a new factorization, and refine
 *        the shifted eigenvalues to suitable accuracy.
 *    (d) For each eigenvalue with a large enough relative separation compute
 *        the corresponding eigenvector by forming a rank revealing twisted
 *        factorization. Go back to (c) for any clusters that remain.
 * 
 * The desired accuracy of the output can be specified by the input
 * parameter ABSTOL.
 * 
 * For more details, see DSTEMR's documentation and:
 * - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
 *   to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
 *   Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
 * - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
 *   Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
 *   2004.  Also LAPACK Working Note 154.
 * - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
 *   tridiagonal eigenvalue/eigenvector problem",
 *   Computer Science Division Technical Report No. UCB/CSD-97-971,
 *   UC Berkeley, May 1997.
 * 
 * 
 * Note 1 : ZHEEVR calls ZSTEMR when the full spectrum is requested
 * on machines which conform to the ieee-754 floating point standard.
 * ZHEEVR calls DSTEBZ and ZSTEIN on non-ieee machines and
 * when partial spectrum requests are made.
 * 
 * Normal execution of ZSTEMR may create NaNs and infinities and
 * hence may abort due to a floating point exception in environments
 * which do not handle NaNs and infinities in the ieee standard default
 * manner. 
 */

extern "C" {
     void zheevr_(char *jobz, char *range, char *uplo, int *n, complex<double> *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, complex<double> *z, int *ldz, int *isuppz, complex<double> *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
}


/* eigenvalues only */

void zheigen (complex<double> *a, double *w, int n, int stat_i, int stat_f){

 //lapack subroutine parameters
 char jobz='N',range='I',uplo='U';
 int lda=n,il=stat_i+1,iu=stat_f+1,m=iu-il+1,ldz=1,*isuppz,lwork=2*n,lrwork=24*n,liwork=10*n,*iwork,info;
 double vl=0,vu=1,abstol=0,*rwork;
 complex<double> *z,*work;
 
 isuppz=new int[2*m];
 iwork=new int[liwork];
 rwork=new double[lrwork];
 work=new complex<double>[lwork];
 
 t=clock();
 zheevr_(&jobz,&range,&uplo,&n,a,&lda,&vl,&vu,&il,&iu,&abstol,&m,w,z,&ldz,isuppz,work,&lwork,rwork,&lrwork,iwork,&liwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }

 delete[] isuppz; 
 delete[] iwork;
 delete[] rwork;
 delete[] work;
}

/* eigenvalues and eigenvectors */

void zheigenv (complex<double> *a, double *w, complex<double> *vw, int n, int stat_i, int stat_f){

 //lapack subroutine parameters
 char jobz='V',range='I',uplo='U';
 int lda=n,il=stat_i+1,iu=stat_f+1,m=iu-il+1,ldz=n,*isuppz,lwork=2*n,lrwork=24*n,liwork=10*n,*iwork,info;
 double vl=0,vu=1,abstol=0,*rwork;
 complex<double> *work;

 isuppz=new int[2*m];
 iwork=new int[liwork];
 rwork=new double[lrwork];
 work=new complex<double>[lwork];
 
 t=clock();
 zheevr_(&jobz,&range,&uplo,&n,a,&lda,&vl,&vu,&il,&iu,&abstol,&m,w,vw,&ldz,isuppz,work,&lwork,rwork,&lrwork,iwork,&liwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }

 delete[] isuppz;  
 delete[] iwork;
 delete[] rwork;
 delete[] work;
} 

/***********************************************************/
/***** Eigenvalue problem of a real assymmetric matrix *****/
/***********************************************************/
/* DGEEV subroutine
 *
 * DGEEV computes for an N-by-N real nonsymmetric matrix A, the
 * eigenvalues and, optionally, the left and/or right eigenvectors.
 *
 * The right eigenvector v(j) of A satisfies
 *                  A * v(j) = lambda(j) * v(j)
 * where lambda(j) is its eigenvalue.
 * The left eigenvector u(j) of A satisfies
 *               u(j)**H * A = lambda(j) * u(j)**H
 * where u(j)**H denotes the conjugate-transpose of u(j).
 *
 * The computed eigenvectors are normalized to have Euclidean norm
 * equal to 1 and largest component real.
 */

extern "C" {
     void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
}	

/* eigenvalues only */

void dgeeigen (double *a, complex<double> *w, int n){

 //lapack subroutine parameters
 char jobvl='N',jobvr='N';
 int lda=n,ldvl=1,ldvr=1,lwork=3*n,info;
 double *wr,*wi,*vl,*vr,*work;

 wr=new double[n];
 wi=new double[n];
 work=new double[lwork];
 
 t=clock();
 dgeev_(&jobvl,&jobvr,&n,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }
 
 for (int i=0;i<n;i++)
  w[i]=complex<double>(wr[i],wi[i]);

 delete[] wr;
 delete[] wi;
 delete[] work;
}

/* eigenvalues and eigenvectors */

void dgeeigenv (double *a, complex<double> *w, complex<double> *vw, int n){

 //lapack subroutine parameters
 char jobvl='N',jobvr='V';
 int lda=n,ldvl=1,ldvr=n,lwork=4*n,info,*vj,pj;
 double *wr,*wi,*vl,*vr,*work;

 wr=new double[n];
 wi=new double[n];
 vr=new double[ldvr*n];
 work=new double[lwork];
 vj=new int[n]();
 
 t=clock();
 dgeev_(&jobvl,&jobvr,&n,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }
 
 for (int i=0;i<n;i++){
  w[i]=complex<double>(wr[i],wi[i]);
  if (wi[i]>0)
   vj[i]=1;
 }

 for (int i=0;i<n;i++){
  pj=0;
  do{
   if (vj[pj]==1){
    vw[i+pj*n]=complex<double>(vr[i+pj*n],vr[i+(pj+1)*n]);
    vw[i+(pj+1)*n]=complex<double>(vr[i+pj*n],-vr[i+(pj+1)*n]);
    pj=pj+2;
   }
   else{
    vw[i+pj*n]=complex<double>(vr[i+pj*n],0);
    pj++;
   }
  }while(pj<n);
 }

 delete[] wr;
 delete[] wi;
 delete[] vr;
 delete[] work;
 delete[] vj;
}

/***********************************************************/
/*** Eigenvalue problem of a complex assymmetric matrix ****/
/***********************************************************/
/* ZGEEV subroutine
 *
 * ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
 * eigenvalues and, optionally, the left and/or right eigenvectors.
 *
 * The right eigenvector v(j) of A satisfies
 *                  A * v(j) = lambda(j) * v(j)
 * where lambda(j) is its eigenvalue.
 * The left eigenvector u(j) of A satisfies
 *               u(j)**H * A = lambda(j) * u(j)**H
 * where u(j)**H denotes the conjugate transpose of u(j).
 *
 * The computed eigenvectors are normalized to have Euclidean norm
 * equal to 1 and largest component real.
 */

extern "C" {
     void zgeev_(char *jobvl, char *jobvr, int *n, complex<double> *a, int *lda, complex<double> *w, complex<double> *vl, int *ldvl, complex<double> *vr, int *ldvr, complex<double> *work, int *lwork, double *rwork, int *info);
}

/* eigenvalues only */

void zgeeigen (complex<double> *a, complex<double> *w, int n){

 //lapack subroutine parameters
 char jobvl='N',jobvr='N';
 int lda=n,ldvl=1,ldvr=1,lwork=2*n,info;
 double *rwork;
 complex<double> *vl,*vr,*work;

 work=new complex<double>[lwork];
 rwork=new double[lwork];
 
 t=clock();
 zgeev_(&jobvl,&jobvr,&n,a,&lda,w,vl,&ldvl,vr,&ldvr,work,&lwork,rwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }

 delete[] work;
 delete[] rwork;
}

/* eigenvalues and eigenvectors */

void zgeeigenv (complex<double> *a, complex<double> *w, complex<double> *vw, int n){

 //lapack subroutine parameters
 char jobvl='N',jobvr='V';
 int lda=n,ldvl=1,ldvr=n,lwork=2*n,info;
 double *rwork;
 complex<double> *vl,*work;

 work=new complex<double>[lwork];
 rwork=new double[lwork];
 
 t=clock();
 zgeev_(&jobvl,&jobvr,&n,a,&lda,w,vl,&ldvl,vw,&ldvr,work,&lwork,rwork,&info);
 t=(clock()-t)/CLOCKS_PER_SEC;
 
 if (info==0)
  cout<<"eigensystem solved in "<<t<<" seconds\n"; 
 else{
  cout<<"An error has occurred solving the eigensystem\n";
  abort();
 }

 delete[] work;
 delete[] rwork;
}

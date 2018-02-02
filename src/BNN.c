

#include <stdio.h>
#include <stdlib.h>
//#include <cstdlib>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rconfig.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

//#define CSTACK_DEFNS 7

#define pi 3.14159265

//#include "nrutil.h"

//
/*
double unif_rand()
{
  double un;
  GetRNGstate();
  un=unif_rand();
  PutRNGstate();
  return un;
}
*/

//begin: "nrutil.c"
int nrerror(error_text)
  char error_text[];
{
    //void exit();
    
    Rprintf("Numerical Recipes run-time error...\n");
    Rprintf("%s\n",error_text);
    Rprintf("...now exiting to system...\n");
    return 0;
}



float *vector(nl,nh)
  int nl,nh;
{
    float *v;
    
    v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl;
}

int *ivector(nl,nh)
  int nl,nh;
{
    int *v;
    
    v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl;
}

double *dvector(nl,nh)
  int nl,nh;
{
    double *v;
    
    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl;
}



float **matrix(nrl,nrh,ncl,nch)
  int nrl,nrh,ncl,nch;
{
    int i;
    float **m;
    
    m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m -= nrl;
    
    for(i=nrl;i<=nrh;i++) {
      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
      if (!m[i]) nrerror("allocation failure 2 in matrix()");
      m[i] -= ncl;
    }
    return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
  int nrl,nrh,ncl,nch;
{
    int i;
    double **m;
    
    m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m) nrerror("allocation failure 1 in dmatrix()");
    m -= nrl;
    
    for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
      if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
      m[i] -= ncl;
    }
    return m;
}

int **imatrix(nrl,nrh,ncl,nch)
  int nrl,nrh,ncl,nch;
{
    int i,**m;
    
    m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
    if (!m) nrerror("allocation failure 1 in imatrix()");
    m -= nrl;
    
    for(i=nrl;i<=nrh;i++) {
      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
      if (!m[i]) nrerror("allocation failure 2 in imatrix()");
      m[i] -= ncl;
    }
    return m;
}



float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
  float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
  int i,j;
  float **m;
  
  m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
  if (!m) nrerror("allocation failure in submatrix()");
  m -= newrl;
  
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;
  
  return m;
}



void free_vector(v,nl,nh)
  float *v;
int nl,nh;
{
  free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
  int *v,nl,nh;
{
    free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
  double *v;
int nl,nh;
{
  free((char*) (v+nl));
}



void free_matrix(m,nrl,nrh,ncl,nch)
  float **m;
int nrl,nrh,ncl,nch;
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
  double **m;
int nrl,nrh,ncl,nch;
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
  int **m;
int nrl,nrh,ncl,nch;
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}



void free_submatrix(b,nrl,nrh,ncl,nch)
  float **b;
int nrl,nrh,ncl,nch;
{
  free((char*) (b+nrl));
}



float **convert_matrix(a,nrl,nrh,ncl,nch)
  float *a;
int nrl,nrh,ncl,nch;
{
  int i,j,nrow,ncol;
  float **m;
  
  nrow=nrh-nrl+1;
  ncol=nch-ncl+1;
  m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
  if (!m) nrerror("allocation failure in convert_matrix()");
  m -= nrl;
  for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
  return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
  float **b;
int nrl,nrh,ncl,nch;
{
  free((char*) (b+nrl));
}

// end: nrutil.c

// begin: "lib.c"
#define TINY 1.0e-20;

void ludcmp(a,n,indx,d)
  int n,*indx;
double **a,*d;
{
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  double *vv;
  
  vv=dvector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big <= 0.0) nrerror("Singular matrix in routine LUDCMP");
      vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free_dvector(vv,1,n);
}

#undef TINY


void lubksb(a,n,indx,b)
  double **a,*b;
int n,*indx;
{
  int i,ii=0,ip,j;
  double sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}



double matrix_logdet(X, n)
  double **X;
int n;
{
  int j, *indx;
 
  double d, logdet;
  
  indx=ivector(1,n);
  ludcmp(X, n, indx, &d);
  for(logdet=0.0,j=1; j<=n; j++) logdet+=log(fabs(X[j][j]));
  
  free_ivector(indx,1,n);
  return logdet;
}


/* Y=inv(X), return d=log(det(X)) */ 
double matrix_inverse(X, Y, n)
  double **X, **Y;
int n;
{
  double d, *col;
  int i, j, *indx;
  
  double logdet;
  
  col=dvector(1,n);
  indx=ivector(1,n);
  
  ludcmp(X, n, indx, &d);
  
  for(logdet=0.0,j=1; j<=n; j++) logdet+=log(fabs(X[j][j]));
  /*
  if(d==0.0){  Rprintf("Singular matrix\n");  return; }
  */
  
  for(j=1; j<=n; j++){
    for(i=1; i<=n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(X, n, indx, col);
    for(i=1; i<=n; i++) Y[i][j]=col[i];
  }
  
  for(i=1; i<=n; i++)
    for(j=1; j<=n; j++) { Y[i][j]=(Y[i][j]+Y[j][i])*0.5; Y[j][i]=Y[i][j]; }
    
    free_dvector(col,1,n); free_ivector(indx,1,n);
    
    return logdet;
}


/* Y=inv(X), return d=log(det(X)) */
int matrix_inverse_diag(X, Y, diag, n)
  double **X, **Y, *diag;
int n;
{
  double d, *col;
  int i, j, *indx;
  
  
  
  col=dvector(1,n);
  indx=ivector(1,n);
  
  ludcmp(X, n, indx, &d);
  for(j=1; j<=n; j++) diag[j]=X[j][j];
  
  for(j=1; j<=n; j++){
    for(i=1; i<=n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(X, n, indx, col);
    for(i=1; i<=n; i++) Y[i][j]=col[i];
  }
  
  for(i=1; i<=n; i++)
    for(j=1; j<=n; j++) { Y[i][j]=(Y[i][j]+Y[j][i])*0.5; Y[j][i]=Y[i][j]; }
    
    free_dvector(col,1,n); free_ivector(indx,1,n);
    
    return 0;
}



double matrix_trace(A,p)
  double **A;
int p;
{
  int i;
  double sum;
  
  for(sum=0.0, i=1; i<=p; i++) sum+=A[i][i];
  
  return sum;
}


int matrix_sum(A,B,C,n,p)
  double **A, **B, **C;
int n, p;
{
  int i, j;
  for(i=1; i<=n; i++)
    for(j=1; j<=p; j++) C[i][j]=A[i][j]+B[i][j];
  return 0;
}


/* Matrix: A: n by p; B: p by m;  C: n by m */
int matrix_multiply(A,B,C,n,p,m)
  double **A, **B, **C;
int n, p, m;
{
  int i, j, k;
  for(i=1; i<=n; i++)
    for(j=1; j<=m; j++){
      C[i][j]=.0;
      for(k=1; k<=p; k++) C[i][j]+=A[i][k]*B[k][j];
    }
    return 0;
}


int matrix_vector_prod(A,b,d,n,p)
  double **A, *b, *d;
int n,p;
{
  int i,j;
  for(i=1; i<=n; i++){
    d[i]=0.0;
    for(j=1; j<=p; j++) d[i]+=A[i][j]*b[j];
  }
  return 0;
}

double vector_matrix_vector(a,X,b,m,n)
  double *a,**X,*b;
int m, n;
{
  double sum;
  int i, j;
  
  for(sum=0.0, i=1; i<=m; i++)
    for(j=1; j<=n; j++) sum+=a[i]*X[i][j]*b[j];
  
  return sum;
}


void copy_vector(a,b,p)
  double *a, *b;
int p;
{
  int i;
  for(i=1; i<=p; i++) b[i]=a[i];
}

void copy_matrix(a,b,n,p)
  double **a,**b;
int n, p;
{
  int i, j;
  for(i=1; i<=n; i++)
    for(j=1; j<=p; j++) b[i][j]=a[i][j];
}


int choldc(a, n, D)
  double **a;
int n;
double **D;
{
  int i, j, k;
  double sum, *p;
  
  p=dvector(1,n);
  for (i=1; i<=n; i++) {
    for (j=i; j<=n; j++) {
      for (sum=a[i][j], k=i-1; k>=1; k--) sum -= a[i][k]*a[j][k];
      if (i==j) {
        if (sum <=0.0){  Rprintf("choldc failed");   }
        p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }
  /* Transfer the lower triangular part of A to the lower triangular matrix D */
  for(i=1; i<=n; i++){
    D[i][i]=p[i];
    for(j=1; j<i; j++){ D[i][j]=a[i][j]; D[j][i]=0.0; }
  }
  free_dvector(p,1,n);
  return 0;
}


/* calculate log(Gamma(s))  */
double loggamma(xx)
  double xx;
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
                          -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
      x += 1.0;
      ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}


/* calculate log(k!) */
double logpum(k)
  int k;
{
    double value;
    int i;
    
    for(value=0.0, i=1; i<=k; i++) value+=log(1.0*i);
    
    return value;
}


/* generate the random variable form Gamma(a,b) */
double Rgamma(a,b)
  double a, b;
{
    int ok;
    double d, q, un, u1=0.0, y, z;
    
    if(a<=0.0 || b<=0.0) {  Rprintf("Gamma parameter error (<0.0)\n"); return 0.0; }
    
    if(a<1.0){  /* Ahrens, (*P).213 */
ok=0;
      while(ok==0){
        un=0.0;
        while(un<=0.0 || un>=1.0){
          
          un=unif_rand();
          
        }
        d=(2.718282+a)/2.718282;
        q=d*un;
        
        if(q<=1.0){
          z=exp(1.0/a*log(q));
          
          u1=unif_rand();
          if(u1<exp(-z)) ok=1;
        }
        else{
          z=-1.0*log((d-q)/a);
          
          u1=unif_rand();
          if(u1<exp((a-1)*log(z))) ok=1;
        }
      } /* end ok */
    }
    else {  /* a>=1.0 Fishman, (*P).214 */
ok=0;
      while(ok==0){
        un=0.0;
        while(un<=0.0 || un>=1.0){
          
          un=unif_rand();
          
        }
        y=-1.0*log(un);
        
        un=unif_rand();
        
        if(u1<exp((a-1)*(log(y)-(y-1)))) { z=a*y; ok=1; }
      }
    }
    
    return z/b;
}


/* Generate a random variable from Beta(1,k), where
the first parameter is 1, the second parameter is b */
double Rbeta(b)
  double b;
{
    double un;
    un=0.0;
    while(un<=0.0 || un>=1.0){
      
      un=unif_rand();
      
    }
    return 1.0-exp(1.0/b*log(un));
}


/* Generate deviates from Dirichlet(a1,a2,\ldots,a_k) */
int RDirichlet(w,a,k)
  double *w,*a;
int k;
{
  double sum;
  int i;
  for(sum=0.0,i=1; i<=k; i++){
    w[i]=Rgamma(a[i],1.0);
    sum+=w[i];
  }
  for(i=1; i<=k; i++) w[i]/=sum;
  return 0;
}


double gasdev()
{
  static int iset=0;
  static double gset;
  double fac,r,v1,v2;
  
  if  (iset == 0) {
    do {
      
      v1=unif_rand()*2.0-1.0;
      
      
      v2=unif_rand()*2.0-1.0;
      
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}


double Rgasdev(mean,variance)
  double mean,variance;
{
    static int iset=0;
    static double gset;
    double fac,r,v1,v2;
    
    if  (iset == 0) {
      do {
        
        v1=unif_rand()*2.0-1.0;
        
        
        v2=unif_rand()*2.0-1.0;
        
        r=v1*v1+v2*v2;
      } while (r >= 1.0);
      fac=sqrt(-2.0*log(r)/r);
      gset=v1*fac;
      iset=1;
      return v2*fac*sqrt(variance)+mean;
    } else {
      iset=0;
      return gset*sqrt(variance)+mean;
    }
}


int RNORM(x,mu,Sigma,p)
  double *x,*mu,**Sigma;
int p;
{
  int i, j;
  double **D, *z;
  
  D=dmatrix(1,p,1,p);
  z=dvector(1,p);
  
  choldc(Sigma,p,D);
  
  for(i=1; i<=p; i++) z[i]=gasdev();
  for(i=1; i<=p; i++){
    x[i]=mu[i];
    for(j=1; j<=i; j++) x[i]+=D[i][j]*z[j];
  }
  
  free_dmatrix(D,1,p,1,p);
  free_dvector(z,1,p);
  
  return 0;
}


int Rwishart(B,df,Sigma,p)
  double **B,df,**Sigma;
int p;
{
  double **Z, **A, *Y;
  int i, j, k;
  
  Z=dmatrix(1,p,1,p);
  A=dmatrix(1,p,1,p);
  Y=dvector(1,p);
  
  for(i=1; i<=p; i++)
    for(j=1; j<=p; j++) Z[i][j]=A[i][j]=B[i][j]=0.0;
  
  choldc(Sigma, p, A);
  
  for(j=1; j<=p; j++)
    for(i=1; i<=p; i++) Z[i][j]=gasdev();
  for(i=1; i<=p; i++) Y[i]=Rgamma(0.5*df,0.5);
  
  B[1][1]=Y[1];
  for(j=2; j<=p; j++){
    B[j][j]=Y[j];
    for(i=1; i<j; i++) B[j][j]+=Z[i][j]*Z[i][j];
    B[1][j]=Z[1][j]*sqrt(Y[1]);
  }
  for(j=2; j<=p; j++)
    for(i=2; i<j; i++){
      B[i][j]=Z[i][j]*sqrt(Y[i]);
      for(k=1; k<=i-1; k++) B[i][j]+=Z[k][i]*Z[k][j];
    }
    for(i=1; i<=p; i++)
      for(j=1; j<i; j++) B[i][j]=B[j][i];
  
  matrix_multiply(A,B,Z,p,p,p);
  
  for(i=1; i<=p; i++)
    for(j=1; j<i; j++){ A[j][i]=A[i][j]; A[i][j]=0.0; }
    
    matrix_multiply(Z,A,B,p,p,p);
  
  free_dmatrix(Z,1,p,1,p);
  free_dmatrix(A,1,p,1,p);
  free_dvector(Y,1,p);
  
  return 0;
}


/* calculated the log-density of  z~gamma(a,b) */
double dloggamma(x,a,b)
  double x,a,b;
{
    double logcon, den;
    logcon=loggamma(a);
    den=log(b)-b*x+(a-1)*log(b*x)-logcon;
    return den;
}


double dloggauss(z,mean,variance)
  double z,mean,variance;
{
    double sum;
    sum=-0.5*log(2.0*pi*variance);
    sum+=-0.5*(z-mean)*(z-mean)/variance;
    return sum;
}

double dlogstudent(z,k)
  double z;
int k;  /* the degree of freedom */
{
  double logprob;
  logprob=-0.5*(k+1)*log(1.0+z*z/k);
  return logprob;
}



double DLOGGAUSS(z,mean,variance,p)
  double *z, *mean, **variance;
int p;
{
  int i;
  double logdet, sum, **mat,*vect, *mu;
  
  mat=dmatrix(1,p,1,p);
  vect=dvector(1,p);
  mu=dvector(1,p);
  
  for(i=1; i<=p; i++) mu[i]=z[i]-mean[i];
  logdet=matrix_inverse(variance,mat,p);
  matrix_vector_prod(mat,mu,vect,p,p);
  
  for(sum=0.0,i=1; i<=p; i++) sum+=mu[i]*vect[i];
  sum*=-0.5;
  sum+=-0.5*logdet-0.5*p*log(2.0*pi);
  
  free_dmatrix(mat,1,p,1,p);
  free_dvector(vect,1,p);
  free_dvector(mu,1,p);
  
  return sum;
}


double Dlogwishart(D,df,Sigma,p)
  double **D,df,**Sigma;
int p;
{
  int i, j;
  double a, sum, logdet1, logdet2, **mt1, **mt2;
  
  mt1=dmatrix(1,p,1,p);
  mt2=dmatrix(1,p,1,p);
  
  for(i=1; i<=p; i++)
    for(j=1; j<=p; j++) mt1[i][j]=D[i][j];
  logdet1=matrix_inverse(mt1,mt2,p);
  
  for(i=1; i<=p; i++)
    for(j=1; j<=p; j++) mt1[i][j]=Sigma[i][j];
  logdet2=matrix_inverse(mt1,mt2,p);
  
  matrix_multiply(mt2,D,mt1,p,p,p);
  
  for(sum=0.0,i=1; i<=p; i++) sum+=mt1[i][i];
  sum*=-0.5;
  sum+=0.5*(df-p-1)*logdet1;
  sum+=-0.5*df*logdet2;
  sum+=-0.5*df*p*log(2.0);
  sum+=-0.25*p*(p-1)*log(pi);
  for(i=1; i<=p; i++){
    a=df-i+1;
    if(a<0.0001) a=0.0001;
    sum+=-loggamma(0.5*a);
  }
  
  free_dmatrix(mt1,1,p,1,p);
  free_dmatrix(mt2,1,p,1,p);
  
  return sum;
}


int uniform_direction(d, n)
  double *d;
int n;
{
  double sum;
  int k;
  
  for(sum=0, k=1; k<=n; k++){
    d[k]=gasdev();
    sum+=d[k]*d[k];
  }
  for(k=1; k<=n; k++)
    d[k]=d[k]/sqrt(sum);
  
  return 0;
}


int dmaxclass(z,n)
  double *z;
int n;
{
  int i, maxi;
  double maxd;
  
  maxd=z[1]; maxi=1;
  for(i=2; i<=n; i++)
    if(z[i]>maxd){ maxd=z[i]; maxi=i; }
    
    return maxi;
} 

double dmaximum(z,n)
  double *z;
int n;
{
  int i;
  double maxd;
  
  maxd=z[1];
  for(i=2; i<=n; i++)
    if(z[i]>maxd){ maxd=z[i]; }
    
    return maxd;
}


int imaxclass(z,n)
  int *z, n;
{
    int i, maxi;
    int maxd;
    
    maxd=z[1]; maxi=1;
    for(i=2; i<=n; i++)
      if(z[i]>maxd){ maxd=z[i]; maxi=i; }
      
      return maxi;
}


int binary_trans(k,l,d)
  int k, l,*d;
{
    int i, j;
    for(i=1; i<=l; i++) d[i]=0;
    j=l;
    while(k>0){
      d[j]=k%2;
      k=k/2;
      j--;
    }
    return 0;
}


double logsum(a,b)
  double a, b;
{
    double sum;
    if(a>b) sum=a+log(1.0+exp(b-a));
    else sum=b+log(1.0+exp(a-b));
    return sum;
}


double maxvector(x,n)
  double *x;
int n;
{
  double max;
  int i;
  max=x[1];
  for(i=2; i<=n; i++)
    if(x[i]>max) max=x[i];
    return max;
}


double minvector(x,n)
  double *x;
int n;
{
  double min;
  int i;
  min=x[1];
  for(i=2; i<=n; i++)
    if(x[i]<min) min=x[i];
    return min;
}



double sample_variance(x,n)
  double *x;
int n;
{
  int i;
  double sum1, sum2, mean, var;
  
  sum1=sum2=0.0;
  for(i=1; i<=n; i++){
    sum1+=x[i];
    sum2+=x[i]*x[i];
  }
  mean=sum1/n;
  var=(sum2-n*mean*mean)/(n-1);
  
  return var;
}



/* Return the value ln[Gamma(xx)] for xx>0 */
double gammln(xx)
  double xx;
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
                          -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
      x += 1.0;
      ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}


#define ITMAX 100
#define EPS 3.0e-7
void gser(gamser,a,x,gln)
  double a,x,*gamser,*gln;
{
    int n;
    double sum,del,ap;
    
    *gln=gammln(a);
    if (x <= 0.0) {
      if (x < 0.0) nrerror("x less than 0 in routine GSER");
      *gamser=0.0;
      return;
    } else {
      ap=a;
      del=sum=1.0/a;
      for (n=1;n<=ITMAX;n++) {
        ap += 1.0;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum)*EPS) {
          *gamser=sum*exp(-x+a*log(x)-(*gln));
          return;
        }
      }
      nrerror("a too large, ITMAX too small in routine GSER");
      return;
    }
}

#undef ITMAX
#undef EPS

#define ITMAX 500
#define EPS 1.0e-9
void gserln(gamser,a,x,gln)
  double a,x,*gamser,*gln;
{
    int n;
    double sum,del,ap;
    
    *gln=gammln(a);
    if (x <= 0.0) {
      if (x < 0.0) nrerror("x less than 0 in routine GSERLN");
      *gamser=0.0;
      return;
    } else {
      ap=a;
      del=sum=1.0/a;
      for (n=1;n<=ITMAX;n++) {
        ap += 1.0;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum)*EPS) {
          *gamser=log(sum)-x+a*log(x)-(*gln);
          return;
        }
      }
      nrerror("a too large, ITMAX too small in routine GSERLN");
      return;
    }
}
#undef ITMAX
#undef EPS



#define ITMAX 100
#define EPS 3.0e-7
void gcf(gammcf,a,x,gln)
  double a,x,*gammcf,*gln;
{
    int n;
    double gold=0.0,g,fac=1.0,b1=1.0;
    double b0=0.0,anf,ana,an,a1,a0=1.0;
    /* float gammln();
    void nrerror(); */
    
    *gln=gammln(a);
    a1=x;
    for (n=1;n<=ITMAX;n++) {
      an=(double) n;
      ana=an-a;
      a0=(a1+a0*ana)*fac;
      b0=(b1+b0*ana)*fac;
      anf=an*fac;
      a1=x*a0+anf*a1;
      b1=x*b0+anf*b1;
      if (a1) {
        fac=1.0/a1;
        g=b1*fac;
        if (fabs((g-gold)/g) < EPS) {
          *gammcf=exp(-x+a*log(x)-(*gln))*g;
          return;
        }
        gold=g;
      }
    }
    nrerror("a too large, ITMAX too small in routine GCF");
}
#undef ITMAX
#undef EPS

#define ITMAX 500
#define EPS 1.0e-9
void gcfln(gammcf,a,x,gln)
  double a,x,*gammcf,*gln;
{
    int n;
    double gold=0.0,g,fac=1.0,b1=1.0;
    double b0=0.0,anf,ana,an,a1,a0=1.0;
    /* float gammln();
    void nrerror(); */
    
    *gln=gammln(a);
    a1=x;
    for (n=1;n<=ITMAX;n++) {
      an=(double) n;
      ana=an-a;
      a0=(a1+a0*ana)*fac;
      b0=(b1+b0*ana)*fac;
      anf=an*fac;
      a1=x*a0+anf*a1;
      b1=x*b0+anf*b1;
      if (a1) {
        fac=1.0/a1;
        g=b1*fac;
        if (fabs((g-gold)/g) < EPS) {
          *gammcf=-x+a*log(x)-(*gln)+log(g);
          return;
        }
        gold=g;
      }
    }
    nrerror("a too large, ITMAX too small in routine GCFLN");
}
#undef ITMAX
#undef EPS



double gammp(a,x)
  double a,x;
{
    double gamser=1.0,gammcf,gln;
    /* void gser(),gcf(),nrerror(); */
    
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMP");
    if (x < (a+1.0)) {
      gser(&gamser,a,x,&gln);
      return gamser;
    } else {
      gcf(&gammcf,a,x,&gln);
      return 1.0-gammcf;
    }
}


double gammq(a,x)
  double a,x;
{
    double gamser=0.0,gammcf,gln;
    // void gcf(),gser(),nrerror();
    
    if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMQ");
    if (x < (a+1.0)) {
      gser(&gamser,a,x,&gln);
      return 1.0-gamser;
    } else {
      gcf(&gammcf,a,x,&gln);
      return gammcf;
    }
}



double gammpln(a,x,tail)
  double a,x;
int *tail;
{
  double gamserln=1.0,gammcfln,gln;
  /* void gser(),gcf(),nrerror(); */
  
  
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMP");
  if (x < (a+1.0)) {
    *tail=1;
    gserln(&gamserln,a,x,&gln);
    return gamserln;
  } else {
    gcfln(&gammcfln,a,x,&gln);
    *tail=-1;
    return gammcfln;
  }
}


/* Return the CDF of the standard normal distribution */
double gauss_cdf(x)
  double x;
{
    double s, prob;
    
    s=0.5*x*x;
    if(x>0) prob=0.5+0.5*gammp(0.5,s); 
    else prob=0.5-0.5*gammp(0.5,s);
    
    return prob;
}


/* Return the log(CDF) of the standard normal distribution */
double gauss_cdf_ln(x)
  double x;
{
    double s, a, b, logprob;
    int tail;
    
    s=0.5*x*x;
    a=log(0.5);
    b=gammpln(0.5,s,&tail)+a;
    
    if(tail==1){
      if(x>0.0) logprob=a+log(1+exp(b-a)); 
      else logprob=a+log(1-exp(b-a));
    }
    else{
      if(x<0.0) logprob=b;
      else logprob=b+log(exp(-b)-1.0);
    }
    
    return logprob;
}



/* return Gamma'(z)/Gamma(z)   */
/* Refer to "mathematics handbook pp.287" */
double diGamma(z)
  double z;
{
    int i;
    double sum, delta, epsilon=3.0e-7;
    
    sum=-1.0/z-0.5772156649;
    
    delta=1.0-1.0/(1+z);
    sum+=delta;
    i=1;
    while(delta>epsilon){
      i++;
      delta=1.0/i-1.0/(i+z);
      sum+=delta;
    }
    
    return sum;
}


/* return Gamma'(z)/Gamma(z)   */
/* Refer to "mathematics handbook pp.287" */
double derivative_gamma(z)
  double z;
{
    int i;
    double sum, delta, epsilon=3.0e-7;
    
    sum=-1.0/z-0.5772156649;
    
    delta=1.0-1.0/(1+z);
    sum+=delta;
    i=1;
    while(delta>epsilon){
      i++;
      delta=1.0/i-1.0/(i+z);
      sum+=delta;
    }
    
    return sum;
}

double dlogGnormal(z,mu,alpha,beta)
  double z, mu, alpha,beta;
{
    double sum;
    
    if(fabs(z-mu)<1.0e-10) sum=log(beta)-log(2.0*alpha)-gammln(1.0/beta);
    else sum=log(beta)-log(2.0*alpha)-gammln(1.0/beta)-exp(beta*log(fabs(z-mu)/alpha));
    
    return sum;
}

double correlation(z1,z2,p)
  double *z1, *z2;
int p;
{
  double ave1, ave2, sq1, sq2, sum;
  int i;
  
  ave1=ave2=0.0; sq1=sq2=0.0;
  for(i=1; i<=p; i++){
    ave1+=z1[i]; ave2+=z2[i];
    sq1+=z1[i]*z1[i]; sq2+=z2[i]*z2[i];
  }
  ave1/=p; ave2/=p;
  sq1=(sq1-p*ave1*ave1)/(p-1);
  sq2=(sq2-p*ave2*ave2)/(p-1);
  
  if(sq1<=0.0 || sq2<=0.0) return 0.0;
  else{
    for(sum=0.0,i=1; i<=p; i++) sum+=(z1[i]-ave1)*(z2[i]-ave2);
    sum/=p;
    sum/=sqrt(sq1*sq2);
    return sum;
  }
}


int permut_sample(sam,n)
  int *sam, n;
{
    int j,k,u,v,*b;
    
    b=ivector(1,n);
    
    for(j=1; j<=n; j++) b[j]=j;
    k=0;
    while(k<n){
      u=0;
      while(u<=0 || u>n-k){
        
        u=floor(unif_rand()*(n-k))+1;
        
      } 
      sam[k+1]=b[u];
      for(v=u; v<n-k; v++) b[v]=b[v+1];
      k++;
    }
    
    return 0;
}



int random_order(x,n)
  int *x, n;
{
    int i, j, k, m, *y;
    
    y=ivector(1,n);
    
    m=n;
    for(i=1; i<=m; i++) y[i]=i;
    for(i=1; i<=n; i++){
      j=0;
      while(j<1 || j>m){
        
        j=(int)(unif_rand()*m)+1;
        
      } 
      x[i]=y[j];
      for(k=j+1; k<=m; k++) y[k-1]=y[k];
      m--;
    }
    
    free_ivector(y,1,n);
    
    return 0;
}


/* Generate a subset sample of size M from the set 1:N */
/*
int subset_sample(x,z,M,N)
int *x, *z,M, N;
{
int i, j, k, m, *y;

y=ivector(1,N);

m=N;
for(i=1; i<=N; i++) y[i]=i;
for(i=1; i<=M; i++){
j=0;
while(j<1 || j>m) j=floor(unif_rand()*m)+1;
x[i]=y[j];
for(k=j+1; k<=m; k++) y[k-1]=y[k];
m--;
}

for(i=1; i<=N-M; i++) z[i]=y[i]; 

free_ivector(y,1,N);

return 0;
}
*/


void indexx(n,arrin,indx)
  int n,*indx;
double *arrin;
{
  int l,j,ir,indxt,i;
  double q;
  
  for (j=1;j<=n;j++) indx[j]=j;
  l=(n >> 1) + 1;
  ir=n;
  for (;;) {
    if (l > 1)
      q=arrin[(indxt=indx[--l])];
    else {
      q=arrin[(indxt=indx[ir])];
      indx[ir]=indx[1];
      if (--ir == 1) {
        indx[1]=indxt;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir) {
      if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
      if (q < arrin[indx[j]]) {
        indx[i]=indx[j];
        j += (i=j);
      }
      else j=ir+1;
    }
    indx[i]=indxt;
  }
}


void indexx_integer(n,arrin,indx)
  int n,*indx;
int *arrin;
{
  int l,j,ir,indxt,i;
  double q;
  
  for (j=1;j<=n;j++) indx[j]=j;
  l=(n >> 1) + 1;
  ir=n;
  for (;;) {
    if (l > 1)
      q=arrin[(indxt=indx[--l])];
    else {
      q=arrin[(indxt=indx[ir])];
      indx[ir]=indx[1];
      if (--ir == 1) {
        indx[1]=indxt;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir) {
      if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
      if (q < arrin[indx[j]]) {
        indx[i]=indx[j];
        j += (i=j);
      }
      else j=ir+1;
    }
    indx[i]=indxt;
  }
}


void indexx_convert_double(n,indx,x,y)
  int n, *indx;
double *x, *y;
{
  int i;
  for(i=1; i<=n; i++) y[indx[i]]=x[i];
}


void indexx_convert_integer(n,indx,x,y)
  int n, *indx;
int *x, *y;
{
  int i;
  for(i=1; i<=n; i++) y[indx[i]]=x[i];
}

/* coefficients for the rational approximants for the normal probit: */
#define a1	(-3.969683028665376e+01)
#define a2	( 2.209460984245205e+02)
#define a3	(-2.759285104469687e+02)
#define a4	( 1.383577518672690e+02)
#define a5	(-3.066479806614716e+01)
#define a6	( 2.506628277459239e+00)
#define b1	(-5.447609879822406e+01)
#define b2	( 1.615858368580409e+02)
#define b3	(-1.556989798598866e+02)
#define b4	( 6.680131188771972e+01)
#define b5	(-1.328068155288572e+01)
#define c1	(-7.784894002430293e-03)
#define c2	(-3.223964580411365e-01)
#define c3	(-2.400758277161838e+00)
#define c4	(-2.549732539343734e+00)
#define c5	( 4.374664141464968e+00)
#define c6	( 2.938163982698783e+00)
#define d1	( 7.784695709041462e-03)
#define d2	( 3.224671290700398e-01)
#define d3	( 2.445134137142996e+00)
#define d4	( 3.754408661907416e+00)
#define p_low	0.02425
#define logp_low -3.719339
#define p_high	(1.0 - p_low)
#define logp_high -0.02454887

/**
* Returns the probit value of the normal distribution CDF.  This is
* an implementation of the algorithm published at
* http://home.online.no/~pjacklam/notes/invnorm/
*/
double inverse_normal_cdf(double p) {
  double q, x;
  
  if(0.0 < p && p < p_low) {
    /* rational approximation for the lower region */
    q = sqrt(-2.0*log(p));
    x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  } else if(p_low <= p && p <= p_high) {
    double r;
    /* rational approximation for the central region */
    q = p - 0.5;
    r = q*q;
    x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
  } else /* if(p_high < p && p < 1.0) */ {
    /* rational approximation for the upper region */
    q = sqrt(-2.0*log(1.0-p));
    x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
  }
  
  if(0.0 < p && p < 1.0) {
    double u, e;
    e = 0.5 * erfc(-x/sqrt(2.0)) - p;
    u = e * sqrt(2.0*M_PI) * exp(x*x/2.0);
    x = x - u/(1.0 + x*u/2.0);
  }
  
  return x;
}

/* the input p is given in logairithm */
double inverse_normal_cdf_log(double logp) {
  double q, x;
  
  if(logp < logp_low) {
    /* rational approximation for the lower region */
    q = sqrt(-2.0*logp);
    x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  } else if(logp_low <= logp && logp <= logp_high) {
    double r;
    /* rational approximation for the central region */
    q = exp(logp) - 0.5;
    r = q*q;
    x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
  } else /* if(p_high < p && p < 1.0) */ {
    /* rational approximation for the upper region */
    q = sqrt(-2.0*(logp+log(exp(-logp)-1)));
    x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
  }
  
  if(logp < 0.0) {
    double u, e;
    e = 0.5 * erfc(-x/sqrt(2.0)) - exp(logp);
    u = e * sqrt(2.0*M_PI) * exp(x*x/2.0);
    x = x - u/(1.0 + x*u/2.0);
  }
  
  return x;
}

#undef a1 
#undef a2  
#undef a3     
#undef a4     
#undef a5    
#undef a6    
#undef b1     
#undef b2      
#undef b3   
#undef b4   
#undef b5   
#undef c1  
#undef c2    
#undef c3   
#undef c4    
#undef c5    
#undef c6     
#undef d1    
#undef d2    
#undef d3    
#undef d4   
#undef p_low 
#undef p_high
#undef logp_low
#undef logp_high

/**
* Returns the quantile for the two-tailed student's t distribution.
* This is an implementation of the algorithm in
* G. W. Hill. "Algorithm 396: Student's t-Quantiles." Communications
* of the ACM 13(10):619--620.  ACM Press, October, 1970.
*/
double inverse_t_cdf(double p, long n) {
  double a, b, c, d, x, y;
  
  if(n < 1) {
    /* you can put your own error handling here */
    Rprintf("tquantile(%f, %d): error: second argument must be >= 1 !", p, n);
    return 0.0;
  } else if(p > 1.0 || p <= 0.0) {
    /* you can put your own error handling here */
    Rprintf("tquantile(%f, %d): error: first argument must be in (0.0, 1.0] !", p, n);
    return 0.0;
  }
  
  if(n == 1) {
    /* special case */
    p *= M_PI_2;
    return cos(p) / sin(p);
  }
  
  a = 1.0 / (n-0.5);
  b = pow(48.0 / a, 2.0);
  c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
  d = ((94.5 / (b + c) - 3.0) / b + 1.0) * sqrt(a * M_PI_2) * (double)n;
  x = d * p;
  y = pow(x, 2.0/(double)n);
  if(y > 0.05 + a) {
    /* asymptotic inverse expansion about the normal */
    x = inverse_normal_cdf(p * 0.5);
    y = x * x;
    if(n < 5) {
      c += 0.3 * ((double)n - 4.5) * (x + 0.6);
      c = (((0.5 * d * x - 0.5) * x - 7.0) * x - 2.0) * x + b + c;
      y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
      y *= a * y;
      if(y > 0.002)
        y = exp(y) - 1.0;
      else
        y += 0.5 * y * y;
    }
  } else
    y = ((1.0/((((double)n + 6.0)/((double)n * y) - 0.089 * d - 0.822) * ((double)n+2.0) * 3.0) + 0.5 / ((double)n+4.0))*y - 1.0) * ((double)n + 1.0) / ((double)n + 2.0) + 1.0 / y;
  
  return sqrt((double)n * y);
}

int iminimum(d,n)
  int *d, n;
{
    int min, i;
    
    min=d[1];
    for(i=2; i<=n; i++){
      if(d[i]<min) min=d[i];
    }
    
    return min;
}

int imaximum(d,n)
  int *d, n;
{
    int max, i;
    
    max=d[1];
    for(i=2; i<=n; i++){
      if(d[i]>max) max=d[i];
    }
    
    return max;
}

int iminmax(d,n,min,max)
  int *d, n,*min,*max;
{
    int i;
    
    *min=d[1]; *max=d[1];
    for(i=2; i<=n; i++){
      if(d[i]<*min) *min=d[i];
      else if(d[i]>*max) *max=d[i];
    }
    
    return 0;
}

/* Generate n MCMC samples from generalized normal(mu,alpha,beta) and 
return the samples z  */
double MCMCGnormal(mu,alpha,beta,n,z)
  double mu, alpha, beta, *z;
int n;
{
  double x,y, fx, fy, un, r;
  int i, accept, warm=50;
  
  x=mu+0.5*alpha;
  fx=dlogGnormal(x,mu,alpha,beta);
  
  for(i=1; i<=n+warm; i++){
    
    y=x+gasdev()*alpha;
    fy=dlogGnormal(y,mu,alpha,beta);
    
    r=fy-fx;
    if(r>0.0) accept=1;
    else{
      un=0.0;
      while(un<=0.0){
        
        un=unif_rand();
        
      }
        
      if(un<exp(r)) accept=1;
      else accept=0;
    }
    
    if(accept==1){ x=y; fx=fy; }
    if(i>warm) z[i-warm]=x;
  }
  
  return 0;
}

// x is a vector of z-values, gamma is a parameter with a default value of 0.1, 
// the output is the estimated mean mu and standard deviation sigma 
//
int EstNull_fdr(x,n,gamma,mu,sigma)
  double *x,gamma,*mu,*sigma;
int n;
{
  double *t, *phiplus, *phiminus, *dphiplus,*dphiminus, *phi, *dphi;
  double gan, shat, uhat;
  double s, tt, a, b, c, da, db;
  int i, j;
  
  t=dvector(1,1000);
  phiplus=dvector(1,1000);
  phiminus=dvector(1,1000);
  dphiplus=dvector(1,1000);
  dphiminus=dvector(1,1000);
  phi=dvector(1,1000);
  dphi=dvector(1,1000);
  
  for(i=1; i<=1000; i++){
    t[i]=1.0*i/200;
    phiplus[i]=1;
    phiminus[i]=1;
    dphiplus[i]=1;
    phi[i]=1;
    dphi[i]=1;
  }
  gan=exp(-gamma*log(1.0*n));
  
  for(i=1; i<=1000; i++){
    s=t[i];
    for(phiplus[i]=0, j=1; j<=n; j++) phiplus[i]+=cos(s*x[j]);
    phiplus[i]/=n;
    
    for(phiminus[i]=0, j=1; j<=n; j++) phiminus[i]+=sin(s*x[j]);
    phiminus[i]/=n;
    
    for(dphiplus[i]=0, j=1; j<=n; j++) dphiplus[i]+=-x[j]*sin(s*x[j]);
    dphiplus[i]/=n; 
    
    for(dphiminus[i]=0, j=1; j<=n; j++) dphiminus[i]+=x[j]*cos(s*x[j]);
    dphiminus[i]/=n; 
    
    phi[i]=sqrt(phiplus[i]*phiplus[i]+phiminus[i]*phiminus[i]);
  }
  
  
  i=1;
  while(phi[i]>gan && i<1000) i++;   
  
  //  Rprintf("i=%d gan=%g phi=%g\n",i, gan, phi[i]);
  
  tt=t[i];
  a=phiplus[i];
  b=phiminus[i];
  da=dphiplus[i];
  db=dphiminus[i];
  c=phi[i];
  
  shat=sqrt(-(a*da+b*db)/(tt*c*c));
  uhat=-(da*b-db*a)/(c*c);
  
  
  *mu=uhat; 
  *sigma=shat; 
  
  
  free_dvector(t,1,1000);
  free_dvector(phiplus,1,1000);
  free_dvector(phiminus,1,1000);
  free_dvector(dphiplus,1,1000);
  free_dvector(dphiminus,1,1000);
  free_dvector(phi,1,1000);
  free_dvector(dphi,1,1000);
  
  return 0;
}


double EpsEst_fdr(x,n,mu,sigma)
  double *x, mu, sigma;
int n;
{
  double *z, tmax, epshat, epsest, *xi,*f,*w,*co,*tt;
  double t, sum1, sum2;
  int i, j, KK, l; 
  
  z=dvector(1,n);
  xi=dvector(0,100);
  f=dvector(0,100);
  w=dvector(0,100);
  co=dvector(0,100);
  
  
  for(i=1; i<=n; i++) z[i]=(x[i]-mu)/sigma;
  for(i=0; i<=100; i++) xi[i]=1.0*i/100;
  
  tmax=sqrt(log(1.0*n));
  
  KK=floor(tmax/0.1); 
  tt=dvector(0,KK);
  for(i=0; i<=KK; i++) tt[i]=0.1*i;
  
  epsest=0.0;
  for(j=0; j<=KK; j++){
    
    t=tt[j];
    for(i=0; i<=100; i++){
      f[i]=exp((t*xi[i])*(t*xi[i])/2);
      w[i]=1-fabs(xi[i]);
    }
    
    for(i=0; i<=100; i++){
      co[i]=0.0;
      for(l=1; l<=n; l++) co[i]+=cos(t*xi[i]*z[l]);
      co[i]/=n;
    }
    
    for(sum1=sum2=0.0, i=0; i<=100; i++){ 
      sum1+=w[i]*f[i]*co[i];
      sum2+=w[i];
    }
    
    epshat=1.0-sum1/sum2;
    if(epshat>epsest) epsest=epshat;
  }
  
  free_dvector(z,1,n);
  free_dvector(xi,0,100);
  free_dvector(f,0,100);
  free_dvector(w,0,100);
  free_dvector(co,0,100);
  free_dvector(tt,0,KK);
  
  return epsest;
}

double standard_deviation(x,n)
  double *x;
int n;
{
  int i;
  double ave, sumq, sd;
  
  for(ave=0.0, sumq=0.0, i=1; i<=n; i++){
    ave+=x[i]; 
    sumq+=x[i]*x[i];
  }
  
  sd=sqrt((sumq-ave*ave/n)/(n-1));
  
  return sd;
}

// end:"lib.c"

// begin: global variables for pratio function

int shortcut=1;
int WARM=1;  
double *refden, tau=0.6, rho=1.0, delta;
double lowE, maxE, bestE, range=50.0, scale=5.0, exL=50.0, exR=0.0;
int grid, extgrid=0, BLOCKSIZE=50;
double hightem=1.0, lowtem=1.0, seltem=1.0, EBICgamma=1.0; 
double mut_step=0.25;   /* 1.0 */
double alpha1=0.2, alpha2=0.2, alpha3=0.2;   /* =1/variance */
double prior_var_alpha=0.01, prior_var_beta=0.01;
double **data_mat, **test_mat,*wei, **data_y, **test_y;
int dim, connection_threshold=25;
double accept_mut,accept_change;
double total_mut, total_change;

// begin: "priorcost.c"
double priorcost(int *W, int hd1, int OUT_UNIT, int P, double lambda)
{
  int i, m1, m2, Okay=1;
  double sum; 
  
  // Okay=check(W);
  
  if(Okay==1){
    for(m1=0, i=1; i<=((P)+1)*(hd1); i++) m1+=W[i]; 
    
    if(shortcut==0){
      m2=0;
      for(i=((P)+1)*(hd1)+1; i<=((P)+1)*(hd1)+((hd1)+1)*(OUT_UNIT); i++) m2+=W[i]; 
    }
    else{
      m2=0;
      for(i=((P)+1)*(hd1)+1; i<=((P)+1)*(hd1)+((P)+(hd1)+1)*(OUT_UNIT); i++) m2+=W[i]; 
    }
    
    if(m1+m2>connection_threshold || m1+m2<1) return 1.0e+100;
    else sum=(m1+m2)*log((lambda)/(1.0-(lambda)));
  }
  else { return 1.0e+100; }
  
  return  -1.0*sum; 
}
// end priorcost.c


// begin "mutationPriorg.c"
int add_connection_prior(int *W,double *fz,double tem,int *region, double **hist, int hd1, int OUT_UNIT, int P, double lambda)
{
  int i,j,m,s,H,k,K,L,d,accept,k1,k2,**node,AA,B;
  double fnew, un, r;
  double p1, p2, pxy, pyx;
  int newregion,cold, cnew;
  
  
  
  node=imatrix(1,(hd1)+(OUT_UNIT),1,2);
  for(cold=0, i=1; i<=dim; i++) cold+=W[i];
  
  // sampling H, if(H=0) update shortcut connections  
  if(shortcut==1){
    H=(hd1)+(OUT_UNIT)+1;
    while(H>(hd1)+(OUT_UNIT)){
      
      H=floor(unif_rand()*((hd1)+(OUT_UNIT)))+1; 
      
    } 
  }  
  else{ 
    H=(hd1)+1;
    while(H>(hd1)){
      
      H=floor(unif_rand()*(hd1))+1;
      
    } 
  } 
  
  // Rprintf("H=%d\n", H);
  
  if(H>(hd1)){
    m=((P)+1)*(hd1)+(1+(hd1)+(P))*(H-(hd1)-1);
    L=W[m+1];
    for(s=1; s<=(P); s++) L+=W[m+1+(hd1)+s]; 
    node[H][1]=L;  node[H][2]=0; 
    
    if(node[H][1]==(P)+1){ total_mut++; free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 0; } 
    else{
      
      B=1;  // type 1 addiation 
      K=((P)+1)-node[H][1];
      k=K+1;
      while(k>K){
        
        k=floor(unif_rand()*K)+1;
        
      } 
      
      if(W[m+1]==0) j=1;
      else j=0; 
      
      // determine the connection to add 
      if(j==k) k1=m+1;
      else{ 
        s=m+1+(hd1); 
        while(j<k && s<m+1+(hd1)+(P)){
          s++;
          if(W[s]==0) j++;
        } 
        k1=s; 
      } 
      W[k1]=1;  
      
      // calculation of proposal probability 
      pxy=1.0/K;  pyx=1.0/(node[H][1]+1); 
    }
  } 
  else{
    
    for(L=0, s=((P)+1)*(H-1)+1; s<=((P)+1)*H; s++) L+=W[s];
    node[H][1]=L;
    
    if(shortcut==1){
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1)+(P))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)+(P)) L+=W[s]; }
    else{
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)) L+=W[s]; 
    }
    node[H][2]=L; 
    
    if((node[H][1]>0 && node[H][2]==0) || (node[H][1]==0 && node[H][2]>0))
    { Rprintf("network structure error\n");  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 1; }
    else{ //no structure error 
      
      if(node[H][1]==(P)+1 && node[H][2]==(OUT_UNIT)){ total_mut++; free_imatrix(node,1, (hd1)+(OUT_UNIT),1,2); return 0; } 
      else if(node[H][2]==(OUT_UNIT) && node[H][1]<(P)+1){
        B=1; // type 1 addition
        
        K=((P)+1)-node[H][1];
        k=K+1;
        while(k>K){
          
          k=floor(unif_rand()*K)+1;
          
        } 
        
        s=((P)+1)*(H-1); j=0;
        while(j<k && s<((P)+1)*H){
          s++;
          if(W[s]==0) j++;
        }
        k1=s;
        W[k1]=1;
        
        // calculation of proposal probability 
        pxy=1.0/K;  
        if((OUT_UNIT)==1) pyx=1.0/(node[H][1]+1); 
        else pyx=0.5/(node[H][1]+1); 
      }
      else if(node[H][1]==(P)+1 && node[H][2]<(OUT_UNIT)){ 
        
        B=1; // type 1 addition 
        
        K=(OUT_UNIT)-node[H][2];
        k=K+1;
        while(k>K){
          
          k=floor(unif_rand()*K)+1;
          
        } 
        
        j=0; d=0;
        while(j<k && d<(OUT_UNIT)){
          
          d++;
          if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H; 
          else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H; 
          
          if(W[s]==0) j++; 
        } 
        k1=s;
        W[k1]=1;
        
        // calculation of proposal probability 
        pxy=1.0/K;  pyx=0.5/(node[H][2]+1); 
      }
      else if(node[H][1]==0 && node[H][2]==0){
        B=2; // type 2 addition 
        
        k=(P)+2;
        while(k>(P)+1){
          
          k=floor(unif_rand()*((P)+1))+1;
          
        } 
        k1=((P)+1)*(H-1)+k;
        
        if((OUT_UNIT)==1) k2=((P)+1)*(hd1)+1+H;
        else{
          k=(OUT_UNIT)+1;
          while(k>(OUT_UNIT)){
            
            k=floor(unif_rand()*(OUT_UNIT))+1;
            
          } 
          
          if(shortcut==1) k2=((P)+1)*(hd1)+(1+(hd1)+(P))*(k-1)+1+H;
          else k2=((P)+1)*(hd1)+(1+(hd1))*(k-1)+1+H;
        }
        W[k1]=1; W[k2]=1;
        
        // calculation of proposal probability 
        pxy=1.0/((P)+1)/(OUT_UNIT);  pyx=1.0;
      }
      else{ // case node[H][1]<(P)+1 && node[H][2]<(OUT_UNIT)  
        
        B=1; // type 1 addition 
        
        
        un=unif_rand();
        
        if(un<0.5) AA=1;
        else AA=2; 
        
        if(AA==1){ 
          K=((P)+1)-node[H][1];
          k=K+1;
          while(k>K){
            
            k=floor(unif_rand()*K)+1;
            
          } 
          
          s=((P)+1)*(H-1); j=0;
          while(j<k && s<((P)+1)*H){
            s++;
            if(W[s]==0) j++;
          }
          k1=s; W[k1]=1;
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/(node[H][1]+1); 
        } 
        else{
          
          K=(OUT_UNIT)-node[H][2];
          k=K+1;
          while(k>K) k=floor(unif_rand()*K)+1;
          
          j=0; d=0;
          while(j<k && d<(OUT_UNIT)){
            d++;
            if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
            else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
            if(W[s]==0) j++;
          }
          k1=s; W[k1]=1; 
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/(node[H][2]+1);
        } 
      }
    } //end no structure error 
  } // end of H is not equal to 0 
  
  
  if(B==1) cnew=cold+1;
  else if(B==2) cnew=cold+2;
  
  if(cold==1) p1=1.0;
  else if(cold==connection_threshold) p1=0.0;
  else p1=1.0/2;
  if(cnew==connection_threshold) p2=1.0;
  else if(cnew==1) p2=0.0;
  else p2=1.0/2;
  
  
  if(B==1){
    fnew=priorcost(W, hd1, OUT_UNIT, P, lambda);
    newregion=cnew;
    
    // Rprintf("B=%d fnew=%g newregion=%d old region=%d, grid=%d\n", B, fnew, newregion,*region, grid);
    
    if(newregion>grid) accept=0;
    else{
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem;
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++; }
    else W[k1]=0;
    total_mut++;
  }
  else if(B==2){
    
    fnew=priorcost(W,hd1, OUT_UNIT, P, lambda);
    newregion=cnew;
    
    //    Rprintf("fnew=%g newregion=%d\n", fnew, newregion);
    
    if(newregion>grid) accept=0;
    else{
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem; 
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++; }
    else { W[k1]=0; W[k2]=0; }
    total_mut++;
  }
  
  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2);
  return 0;
  
}



int del_connection_prior(int *W,double *fz,double tem,int *region, double **hist, int hd1, int OUT_UNIT, int P, double lambda)
{
  int i,j,k,K,m,s,d,H,L,accept,k1=1,k2,AA,B=1;
  double fnew,un,r;
  double p1, p2, pxy=1.0, pyx=1.0;
  int newregion,cold, cnew, **node;
  
  
  node=imatrix(1,(hd1)+(OUT_UNIT),1,2);
  for(cold=0, i=1; i<=dim; i++) cold+=W[i];
  
  
  // sampling H, if(H=0) update shortcut connections  
  if(shortcut==1){
    H=(hd1)+(OUT_UNIT)+1;
    while(H>(hd1)+(OUT_UNIT)) H=floor(unif_rand()*((hd1)+(OUT_UNIT)))+1;
  }
  else{
    H=(hd1)+1;
    while(H>(hd1)) H=floor(unif_rand()*(hd1))+1;
  }
  
  // Rprintf("del H=%d\n", H);
  
  if(H>(hd1)){
    m=((P)+1)*(hd1)+(1+(hd1)+(P))*(H-(hd1)-1);
    L=W[m+1];
    for(s=1; s<=(P); s++) L+=W[m+1+(hd1)+s];
    node[H][1]=L;  node[H][2]=0;
    
    if(node[H][1]==0){ total_mut++; free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 0; }
    else{
      
      B=1;  // type 1 addiation 
      K=node[H][1];
      k=K+1;
      while(k>K) k=floor(unif_rand()*K)+1;
      
      if(W[m+1]==1) j=1;
      else j=0;
      
      // determine the connection to add 
      if(j==k) k1=m+1;
      else{
        s=m+1+(hd1);
        while(j<k && s<m+1+(hd1)+(P)){
          s++;
          if(W[s]==1) j++;
        }
        k1=s;
      }
      W[k1]=0;
      
      // calculation of proposal probability 
      pxy=1.0/K;  pyx=1.0/(((P)+1)-(node[H][1]-1));
    }
  } 
  else{ // H is not equal to 0 
    
    for(L=0, s=((P)+1)*(H-1)+1; s<=((P)+1)*H; s++) L+=W[s];
    node[H][1]=L;
    
    if(shortcut==1){
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1)+(P))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)+(P)) L+=W[s]; }
    else{
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)) L+=W[s];
    }
    node[H][2]=L; 
    
    if((node[H][1]>0 && node[H][2]==0) || (node[H][1]==0 && node[H][2]>0))
    { Rprintf("network structure error\n");  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 1; }
    else{ //no structure error 
      
      if(node[H][1]==0 && node[H][2]==0){ total_mut++; free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 0; }
      else if(node[H][1]==1 && node[H][2]==1){ 
        
        B=2; // type 2 deletion 
        
        s=((P)+1)*(H-1)+1; j=1;
        while(W[s]==0 && j<=(P)+1){ s++; j++; }
        if(j==(P)+2) Rprintf("network error\n"); 
        k1=s; 
        
        
        if((OUT_UNIT)==1) k2=((P)+1)*(hd1)+1+H;
        else{
          
          d=1; s=((P)+1)*(hd1)+1+H;
          while(W[s]==0 && d<(OUT_UNIT)){
            d++;
            if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
            else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
          }
          k2=s; 
        }
        
        W[k1]=0; W[k2]=0; 
        
        
        // calculation of proposal probability 
        pxy=1.0; pyx=1.0/((P)+1)/(OUT_UNIT);   
      } 
      else if(node[H][1]==1 && node[H][2]>1){ 
        
        B=1; //type 1 deletion 
        
        K=node[H][2];
        k=K+1;
        while(k>K) k=floor(unif_rand()*K)+1;
        
        j=0; d=0;
        while(j<k && d<(OUT_UNIT)){
          d++;
          if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
          else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
          
          if(W[s]==1) j++;
        }
        k1=s; W[k1]=0;
        
        // calculation of proposal probability 
        pxy=1.0/K;  pyx=0.5/((OUT_UNIT)-(node[H][2]-1)); 
      }
      else if(node[H][1]>1 && node[H][2]==1){
        
        B=1; 
        K=node[H][1];
        k=K+1;
        while(k>K) k=floor(unif_rand()*K)+1;
        
        s=((P)+1)*(H-1); j=0;
        while(j<k && s<((P)+1)*H){
          s++;
          if(W[s]==1) j++;
        }
        k1=s;  W[k1]=0;
        
        // calculation of proposal probability 
        pxy=1.0/K;  
        if((OUT_UNIT)==1) pyx=1.0/(((P)+1)-(node[H][1]-1)); 
        else pyx=0.5/(((P)+1)-(node[H][1]-1)); 
      }
      else if(node[H][1]>1 && node[H][2]>1){ 
        
        B=1; 
        
        un=unif_rand();
        if(un<0.5) AA=1;
        else AA=2;
        
        if(AA==1){
          K=node[H][1];
          k=K+1;
          while(k>K) k=floor(unif_rand()*K)+1;
          
          s=((P)+1)*(H-1); j=0;
          while(j<k && s<((P)+1)*H){
            s++;
            if(W[s]==1) j++;
          }
          k1=s; W[k1]=0;
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/(((P)+1)-(node[H][1]-1));
        }
        else{ 
          K=node[H][2];
          k=K+1;
          while(k>K) k=floor(unif_rand()*K)+1;
          
          j=0; d=0;
          while(j<k && d<(OUT_UNIT)){
            d++;
            if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
            else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
            if(W[s]==1) j++;
          }
          k1=s; W[k1]=0;
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/((OUT_UNIT)-(node[H][2]-1));
        }
      }
    } // end of no structure error 
  } // end H is not equal to 0 
  
  
  
  if(B==1) cnew=cold-1;
  else if(B==2) cnew=cold-2;
  if(cold==connection_threshold) p1=1.0;
  else if(cold==1) p1=0.0;
  else p1=1.0/2;
  if(cnew==1) p2=1.0;
  else if(cnew==connection_threshold) p2=0.0;
  else p2=1.0/2;
  
  
  if(B==1){
    fnew=priorcost(W,hd1, OUT_UNIT, P, lambda);
    newregion=cnew; 
    
    if(newregion>grid || cnew==0) accept=0; 
    else{
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem;
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++; }
    else W[k1]=1;
    total_mut++;
  }
  else{
    fnew=priorcost(W,hd1, OUT_UNIT, P, lambda);
    newregion=cnew;
    
    if(newregion>grid || cnew==0) accept=0;
    else{
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem;
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++;  }
    else { W[k1]=1; W[k2]=1; }
    total_mut++;
  }
  
  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2);
  return 0;
}

int mutation_SAMC_prior(int *WW,double *fvalue,double t,int *pos,double **hist, int hd1, int OUT_UNIT, int P, double lambda)
{
  int i, m, region;
  double un, fx;
  
  
  fx=*fvalue; 
  region=*pos;
  
  for(m=0, i=1; i<=dim; i++) m+=WW[i];
  
  un=unif_rand();
  
  // Rprintf("m=%d un=%g\n", m, un);
  
  if(m<=1) add_connection_prior(WW,&fx,t,&region,hist, hd1, OUT_UNIT, P, lambda);
  else if(m==connection_threshold) del_connection_prior(WW,&fx,t,&region,hist, hd1, OUT_UNIT, P, lambda);
  else{
    if(un<1.0/2) add_connection_prior(WW,&fx,t,&region,hist, hd1, OUT_UNIT, P, lambda);
    else  del_connection_prior(WW,&fx,t,&region,hist, hd1, OUT_UNIT, P, lambda);
  }
  
  *fvalue=fx;  *pos=region; 
  
  
  return 0;
}
//PutRNGstate();
// end: mutationPriorg.c


void pratio(int *OUT_UNIT, int *P, int *hd1, double *lambda,
              int *total_iteration, int *popN)
{

  
  GetRNGstate();
  
  
  double *fvalue,**hist, *t, total_weight, *locfreq, ave, max, sum;
  int i, j, k0, iter, m,s1,s2,dim0,**WW,*pos;
  double w, linearity, nonlinearity;
  FILE *ins;
  
  int warm2=ceil((double)*total_iteration/5.0),step=ceil((double)*total_iteration/100000.0),stepscale=ceil((double)*total_iteration/2000.0);
  if(shortcut==0){ dim=((*P)+1)*(*hd1)+(1+(*hd1))*(*OUT_UNIT); dim0=((*P)+1)*(*hd1); }
  else{
    dim0=((*P)+1)*(*hd1);
    dim=((*P)+1)*(*hd1)+(1+(*hd1)+(*P))*(*OUT_UNIT);
  }
  if(connection_threshold>dim) connection_threshold=dim;
  
  grid=connection_threshold;
  hist=dmatrix(1,grid,1,3);
  refden=dvector(1,grid);
  pos=ivector(1,(*popN));
  locfreq=dvector(1,grid);
  WW=imatrix(1,(*popN),1,dim);
  fvalue=dvector(1,(*popN));
  t=dvector(1,(*popN));
  
  /* 
  ins=fopen("bb.log", "a");
  fprintf(ins, "\nseed=%d HD=%d (*P)=%d\n", stime,(*hd1),(*P));
  fclose(ins);
  */
  /*
  t[1]=hightem;  t[(*popN)]=lowtem;
  q=(1.0/lowtem-1.0/hightem)/((*popN)-1);
  for(i=2; i<=(*popN)-1; i++) t[i]=1.0/(q+1.0/t[i-1]);
  */
  for(i=1; i<=(*popN); i++) t[i]=1.0;
  
  /*
  ins=fopen("bb.log", "a");
  fprintf(ins, "lambda=%g\n", (*lambda));
  fclose(ins);
  */

  // #pragma omp parallel for private(i,j,iter,k0,m) default(shared)
  for(i=1; i<=(*popN); i++){
    
    for(j=1; j<=dim; j++) WW[i][j]=0;
    k0=6;
    while (k0>5) k0=floor(unif_rand()*5)+1;
    for(j=1; j<=k0; j++){
      m=(*P)+2;
      while(m>(*P)+1) m=floor(unif_rand()*((*P)+1))+1;
      WW[i][m]=1;
    }
    WW[i][(*hd1)*((*P)+1)+2]=1;
    fvalue[i]=priorcost(WW[i],*hd1, *OUT_UNIT, *P, *lambda);
    for(pos[i]=0, j=1; j<=dim; j++) pos[i]+=WW[i][j];
    //Rprintf(" pos[%d]=%d fvalue=%g\n",i,pos[i],fvalue[i]);
  }
  
  
  /* Initialize the configurations */
  for(sum=0.0, i=1; i<=grid; i++){ refden[i]=1.0; sum+=refden[i]; }
  // for(i=0; i<=grid; i++) refden[i]/=sum;
  for(i=1; i<=grid; i++){
    hist[i][1]=i;
    hist[i][2]=0.0;
    hist[i][3]=0.0;
  }
  
  
  total_weight=0.0;  linearity=nonlinearity=0.0; 
  accept_mut=0; total_mut=0;
  /* 
  instr=fopen("bb", "a");
  if(instr==NULL){ Rprintf("can't write to file\n");  }
  */
  
  for(iter=1; iter<=(*total_iteration)+warm2; iter++){
    
    if(iter<=WARM*stepscale) delta=rho;
    else delta=rho*exp(-tau*log(1.0*(iter-(WARM-1)*stepscale)/stepscale));
    
    // #pragma omp parallel for private(i) default(shared)
    for(i=1; i<=(*popN); i++){
      
      mutation_SAMC_prior(WW[i],&(fvalue[i]),t[i],&(pos[i]),hist, *hd1,*OUT_UNIT,*P, *lambda);
      
       //if(iter%10==0) Rprintf("iter=%d pos[%d]=%d fvalue=%g\n",iter,i,pos[i],fvalue[i]);
    }
    
    
    for(j=1; j<=grid; j++) locfreq[j]=0;
    for(i=1; i<=(*popN); i++) locfreq[pos[i]]+=1.0; 
    
    for(j=1; j<=grid; j++){
      hist[j][2]+=delta*(1.0*locfreq[j]/(*popN)-1.0*refden[j]/grid);
      hist[j][3]+=locfreq[j];
    } 
    
    /* weight normalization */
    if(iter==warm2){
      for(sum=0.0,k0=0,i=1; i<=grid; i++)
        if(hist[i][3]<=0.0){ sum+=1.0*refden[i]/grid; k0++; }
        if(k0>0) ave=sum/k0;
        else ave=0.0;
        for(i=1; i<=grid; i++) hist[i][2]=hist[i][2]+log(refden[i]+ave);
        max=hist[1][2];
        for(i=2; i<=grid; i++)
          if(hist[i][2]>max) max=hist[i][2];
          for(sum=0.0, i=1; i<=grid; i++){ hist[i][2]-=max; sum+=exp(hist[i][2]); }
          for(i=1; i<=grid; i++) hist[i][2]=hist[i][2]-log(sum)+log(100.0);
    }
    
    
    if(iter>warm2 && iter%step==0){
      for(i=1; i<=(*popN); i++){
        
        // output the network
        /*
        if(iter%1000==0){ 
        fprintf(instr, "%g %g ", fvalue[i], hist[pos[i]][2]);
        for(j=1; j<=dim; j++){ if(WW[i][j]!=0) fprintf(instr, " %d", j); }
        fprintf(instr, "\n");
        }
        */
        
        w=exp(hist[pos[i]][2]); total_weight+=w;
        
        for(s1=s2=0, j=1; j<=dim; j++){
          if(WW[i][j]==1){
            if(j<=dim0) s1++;
            else if((j-dim0)%(1+(*hd1)+(*P))==1) s2++;
            else if((j-dim0)%(1+(*hd1)+(*P))>=2 && (j-dim0)%(1+(*hd1)+(*P))<=1+(*hd1)) s1++;
            else s2++;
          }
        }
        if(s1==0 && s2>0) linearity+=w; 
        else nonlinearity+=w; 
        
        // if(iter%10000==0) Rprintf("w=%g\n", w);
      }
  } /* end for iter>warm */
} /* end for iter */
        //   fclose(instr);
        
        /* 
  ins=fopen("bb.log", "a");
  fprintf(ins, "linearity prob=%g\n", linearity/total_weight);
  fprintf(ins, "nonlinearity prob=%g\n", nonlinearity/total_weight);
  fclose(ins);
  
  ins=fopen("bb.log", "a");
  fprintf(ins, "linearity prob=%g\n", linearity/total_weight);
  fprintf(ins, "nonlinearity prob=%g\n", nonlinearity/total_weight);
  fclose(ins);
  
  */
  //Rprintf(" %g %g %g\n", linearity/total_weight, nonlinearity/total_weight, nonlinearity/linearity);
  ins=fopen("bb.BF", "w");
  fprintf(ins, " %g %g %g\n", linearity/total_weight, nonlinearity/total_weight, nonlinearity/linearity);
  fclose(ins); 
  
  
  for(sum=0.0,k0=0,i=1; i<=grid; i++)
    if(hist[i][3]<=0.0){ sum+=1.0*refden[i]/grid; k0++; }
    if(k0>0) ave=sum/k0;
    else ave=0.0;
    for(i=1; i<=grid; i++) hist[i][2]=hist[i][2]+log(refden[i]+ave);
    max=hist[1][2];
    for(i=2; i<=grid; i++)
      if(hist[i][2]>max) max=hist[i][2];
      for(sum=0.0, i=1; i<=grid; i++){ hist[i][2]-=max; sum+=exp(hist[i][2]); }
      for(i=1; i<=grid; i++) hist[i][2]=hist[i][2]-log(sum)+log(100.0);
      /*
      ins=fopen("bb.est", "a");
      fprintf(ins, "delta=%g \n", delta);
      if(ins==NULL){ Rprintf("Can't write to file\n"); }
      for(i=1; i<=grid; i++){
      fprintf(ins, "%5d  %10.6f  %10.6f  %10.6f  %g\n",i,hist[i][1],exp(hist[i][2]),
              hist[i][3],hist[i][2]);
      }
      fclose(ins);
      */
      PutRNGstate();
      free_imatrix(WW,1,(*popN),1,dim);
      free_dvector(fvalue,1,(*popN));
      free_dvector(t,1,(*popN));
      free_dmatrix(hist,1,grid,1,3);
      free_dvector(refden,1,grid);
      free_ivector(pos,1,(*popN));
      free_dvector(locfreq,1,grid);
      
      // return 0;
      }





// begin: "cost.c"
#define LOGISTIC(z) 1.0/(1.0+exp(-z))
int fden(double *ox,double *z,int *W,double *ps2, int hd1, int OUT_UNIT, int P)
{
  int i, j, k;
  double ps1[(hd1)+1],ave; 
  
  /* calculate the output of the first hidden layer */
  for(i=1; i<=(hd1); i++){
    k=((P)+1)*(i-1);
    ps1[i]=z[k+1]*W[k+1];
    for(j=k+2; j<=k+(P)+1; j++) ps1[i]+=ox[j-k-1]*z[j]*W[j];
    ps1[i]=tanh(ps1[i]);
  }
  
  
  /* calculate the predicted mean from the mean unit */
  if(shortcut==0){
    for(i=1; i<=(OUT_UNIT); i++){
      k=((P)+1)*(hd1)+(1+(hd1))*(i-1);
      ave=z[k+1]*W[k+1];
      for(j=k+2; j<=k+(hd1)+1; j++) ave+=ps1[j-k-1]*z[j]*W[j];
      ps2[i]=ave;
    }
  }
  else{
    for(i=1; i<=(OUT_UNIT); i++){
      k=((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1);
      ave=z[k+1]*W[k+1];
      for(j=k+2; j<=k+(hd1)+1; j++) ave+=ps1[j-k-1]*z[j]*W[j];
      for(j=k+(hd1)+2; j<=k+(P)+(hd1)+1; j++) ave+=ox[j-k-(hd1)-1]*z[j]*W[j];
      ps2[i]=ave;
    }
  }
  
  return 0;
}

double cost(double *z,int *W,double *aic0,double *bic0,double *ebic0, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  double *ps2,sum,sum0,sum1,sum2,sum3,sum4; 
  int i, j,k, k0, m1, m2, OK;
  
  
  // the connection weights from hidden to output: k0+2, ..., k0+(hd1)+1  
  k0=((P)+1)*(hd1); 
  OK=1;  j=1; 
  if(z[k0+1+j]*W[k0+1+j]<0.0) OK=0;
  while(OK==1 && j<(hd1)){ 
    j++;
    if(z[k0+1+j]*W[k0+1+j]<0.0) OK=0; 
    else if(z[k0+j]*W[k0+j] <z[k0+1+j]*W[k0+1+j]) OK=0;
  }
  
  if(OK==1){
    for(k=0, j=1; j<=dim; j++) k+=W[j];
    if(k>connection_threshold || k<1) OK=0;
  }
  
  
  if(OK==0) return 1.0e+100;
  else{
    
    ps2=dvector(1,(OUT_UNIT));
    
    for(sum0=0.0, i=1; i<=(data_num); i++){
      fden(data_mat[i],z,W,ps2,hd1,OUT_UNIT,P);
      for(j=1; j<=(OUT_UNIT); j++) 
        sum0+=(data_y[i][j]-ps2[j])*(data_y[i][j]-ps2[j]); 
    }
    sum0=-(0.5*(data_num)+prior_var_alpha)*log(0.5*sum0+prior_var_beta); 
    
    
    /*Normal(0,1/a1) prior for the weights from input to the first hidden layer */
    sum1=0.0; m1=0;
    for(i=1; i<=((P)+1)*(hd1); i++){  sum1+=z[i]*z[i]*W[i];  m1+=W[i]; }
    sum1=-0.5*alpha1*sum1+0.5*m1*log(alpha1);
    
    
    /*Normal(0,1/a2) prior for the weights from the first layer to the second and
    the weights from the input to the second: Mean Node */
    if(shortcut==0){
      sum2=0.0;  m2=0; 
      for(i=((P)+1)*(hd1)+1; i<=((P)+1)*(hd1)+((hd1)+1)*(OUT_UNIT); i++){ 
        sum2+=z[i]*z[i]*W[i]; m2+=W[i]; }
    }
    else{
      sum2=0.0;  m2=0;
      for(i=((P)+1)*(hd1)+1; i<=((P)+1)*(hd1)+((P)+(hd1)+1)*(OUT_UNIT); i++){ 
        sum2+=z[i]*z[i]*W[i];  m2+=W[i]; }
    }
    sum2=-0.5*alpha2*sum2+0.5*m2*log(alpha2);
    sum3=sum1+sum2-0.5*(m1+m2)*log(2.0*pi);
    
    sum4=(m1+m2)*log((lambda)/(1.0-(lambda)));
    
    sum=sum0+sum3+sum4;
    
    *aic0=-sum0+1.0*(m1+m2); /* 1/2 AIC*/ 
    *bic0=-sum0+0.5*log(1.0*(data_num))*(m1+m2);  /* 1/2 BIC */
    *ebic0=-sum0+0.5*log(1.0*(data_num))*(m1+m2)+EBICgamma*(m1+m2)*log(1.0*dim);
    
    sum*=-1.0;
    if(sum>1.0e+100) sum=1.0e+100;
    
    free_dvector(ps2,1,(OUT_UNIT));
    return sum;
  }
}

// end: cost.c


// begin: "mutationMH.c"
int Metropolis_mut(double *z,double *fz,double tem,int *W,double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  double *y, *d, fnew, r,un,aicnew,bicnew, ebicnew;
  int i, k,accept, blocksize, tail, m1, m2;
  
  /*
  if(check(W)==1){Rprintf("network simution error 1\n"); exit; }
  */
  
  y=dvector(1,dim+1);  d=dvector(1,dim+1);
  
  for(i=1; i<=dim+1; i++) y[i]=z[i];
  for(tail=0, i=1; i<=dim; i++) tail+=W[i];
  
  k=0; 
  while(tail>0){
    
    blocksize=BLOCKSIZE;
    if(tail<blocksize) blocksize=tail;
    
    uniform_direction(d,blocksize);
    un=gasdev()*sqrt(tem)*mut_step; /* 2 6.0 9.0 */
  
  m1=k+1; 
  for(i=1; i<=blocksize; i++){
    k++;
    while(k<=dim && W[k]==0) k++;
    if(k<=dim && W[k]==1) y[k]=z[k]+un*d[i];
  }
  m2=k;
  
  fnew=cost(y, W,&aicnew,&bicnew,&ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
  r=-1.0*(fnew-(*fz))/tem;
  
  if(r>0.0) accept=1;
  else{
    un=0.0;
    while(un<=0.0) un=unif_rand();
    if(un<exp(r)) accept=1;
    else accept=0;
  }
  
  if(accept==1){
    for(i=m1; i<=m2; i++) z[i]=y[i];
    *fz=fnew; *aic=aicnew; *bic=bicnew; *ebic=ebicnew;
    accept_mut++;
  }
  else{
    for(i=m1; i<=m2; i++) y[i]=z[i];
  }
  total_mut++;
  
  tail-=blocksize;
  }
  
  free_dvector(y,1,dim+1);  free_dvector(d,1,dim+1);
  /*
  if(check(W)==1){ Rprintf("network simulation error 2\n"); exit; }
  */
  
  return 0;
}



int add_connection(double *z,double *fz,double tem,int *W,double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  int i,j,k,m,L,s,accept, k2=1,*D,**node,A,B;
  double fnew, un, mean, var, r,aicnew,bicnew, ebicnew;
  double p1, p2;
  int cold, cnew;
  
  
  /*
  if(check(W)==1){ Rprintf("network starting error\n"); exit; }
  */
  
  
  mean=var=0.0;
  for(m=0, i=1; i<=dim; i++){ m+=W[i]; mean+=z[i]*W[i]; var+=z[i]*z[i]*W[i]; }
  mean/=m;  var=(var/m-mean*mean)*m/(m-1); 
  cold=m;
  m=dim-cold;  /* the number of 0's in W, the possible connections to add */
  
  
  /* determine the connection to add: k */
  s=m+1;
  while(s>m) s=(int)(unif_rand()*m)+1;
  j=0; k=0;
  while(j<s && k<=dim){
    k++;
    if(W[k]==0) j++;
  }
  if(W[k]==1) { Rprintf("Mutation adding error 1\n"); return 1; }
  
  node=imatrix(1,(hd1),1,2);
  for(i=1; i<=(hd1); i++){
    for(L=0, s=(i-1)*((P)+1)+1; s<=i*((P)+1); s++) L+=W[s];
    node[i][1]=L;
    for(L=0,s=((P)+1)*(hd1)+1+i;s<=((P)+1)*(hd1)+(1+(hd1)+(P))*((OUT_UNIT)-1)+1+i;s+=1+(hd1)+(P))
      L+=W[s];
    node[i][2]=L;
    if((node[i][1]==0 && node[i][2]>0)||(node[i][1]>0 && node[i][2]==0)){
      Rprintf("node error mutation adding\n"); return 1; }
  }
  
  if(k>((P)+1)*(hd1)){
    B=0; i=(OUT_UNIT);
    while(i>=1 && B==0){
      if(k>((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1)+1+(hd1)) B=1;
      else if(k>((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1)+1){
        A=k-((P)+1)*(hd1)-(1+(hd1)+(P))*(i-1)-1;
        if(node[A][2]>=1) B=1;
        else B=2;
      }
      else if(k>((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1)) B=1;
      i--;
    }
  }
  else{
    if(k%((P)+1)==0) A=k/((P)+1);
    else A=(int)(1.0*k/((P)+1))+1;
    if(A>(hd1)) A=(hd1);
    
    if(node[A][2]==0 && node[A][1]==0) B=3;
    else B=1;
  }
  
  if(B==1){ k2=0; W[k]=1; }
  else if(B==2){
    k2=(P)+2;
    while(k2>(P)+1) k2=(int)(unif_rand()*((P)+1))+1;
    
    A=(k-((P)+1)*(hd1))%(1+(hd1)+(P))-1; 
    
    k2+=((P)+1)*(A-1);
    if(W[k2]==1) { Rprintf("Mutation adding error 2\n"); return 1; }
    else W[k2]=1;
    W[k]=1;
  }
  else if(B==3){
    if(k%((P)+1)==0) A=k/((P)+1);
    else A=(int)(1.0*k/((P)+1))+1;
    if(A>(hd1)) A=(hd1);
    
    k2=(OUT_UNIT)+1;
    while(k2>(OUT_UNIT)) k2=(int)(unif_rand()*(OUT_UNIT))+1;
    
    k2=((P)+1)*(hd1)+(k2-1)*(1+(hd1)+(P))+1+A;
    if(W[k2]==1) { Rprintf("Mutation adding error 3\n"); return 1; }
    else W[k2]=1;
    W[k]=1;
  }
  
  for(cnew=0, i=1; i<=dim; i++) cnew+=W[i];
  if(cnew>dim){ accept=0;  total_mut++;
  if(B==1) W[k]=0;
  else{ W[k]=0; W[k2]=0; }
  free_imatrix(node,1,(hd1),1,2);
  return 0;
  }
  
  
  for(i=1; i<=(hd1); i++){
    for(L=0, s=(i-1)*((P)+1)+1; s<=i*((P)+1); s++) L+=W[s];
    node[i][1]=L;
    for(L=0,s=((P)+1)*(hd1)+1+i;s<=((P)+1)*(hd1)+(1+(hd1)+(P))*((OUT_UNIT)-1)+1+i;s+=1+(hd1)+(P))
      L+=W[s];
    node[i][2]=L;
    if((node[i][1]==0 && node[i][2]>0)||(node[i][1]>0 && node[i][2]==0)){
      Rprintf("node error mutation adding\n"); return 1; }
  }
  
  D=ivector(1,dim);
  
  for(j=0, i=1; i<=((P)+1)*(hd1); i++)
    if(W[i]==1){
      if(i%((P)+1)==0) A=i/((P)+1);
      else A=(int)(1.0*i/((P)+1))+1;
      if(node[A][1]==1 && node[A][2]>1) A=0;
      else{ j++; D[j]=i; }
    }
    
    for(L=1; L<=(OUT_UNIT); L++){
      
      i=((P)+1)*(hd1)+(L-1)*(1+(hd1)+(P))+1;
      if(W[i]==1){ j++; D[j]=i; }
      
      for(i=((P)+1)*(hd1)+(L-1)*(1+(hd1)+(P))+2; i<=((P)+1)*(hd1)+(L-1)*(1+(hd1)+(P))+1+(hd1); i++)
        if(W[i]==1){
          A=i-((P)+1)*(hd1)-(L-1)*(1+(hd1)+(P))-1;
          if(node[A][2]>1){ j++; D[j]=i; }
          else if(node[A][2]==1 && node[A][1]==1){ j++; D[j]=i; }
        }
        
        for(i=((P)+1)*(hd1)+(L-1)*(1+(hd1)+(P))+1+(hd1)+1;i<=((P)+1)*(hd1)+(L-1)*(1+(hd1)+(P))+1+(hd1)+(P);i++)
          if(W[i]==1){ j++; D[j]=i; }
    }
    L=j;  /* the number of connections deletable */
  
  
  if(cold==3) p1=2.0/3;
  else if(cold==dim) p1=0.0;
  else p1=1.0/3;
  if(cnew==dim) p2=2.0/3;
  else if(cnew==3) p2=0.0;
  else p2=1.0/3;
  
  
  if(B==1){
    z[k]=gasdev()*sqrt(var)+mean;
    fnew=cost(z,W,&aicnew,&bicnew,&ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    
    r=-1.0*(fnew-(*fz))/tem-log(1.0*L)+log(1.0*m)-dloggauss(z[k],mean,var);
    r+=log(p2/p1);
    
    if(r>0.0) accept=1;
    else{
      un=0.0;
      while(un<=0.0) un=unif_rand();
      if(un<exp(r)) accept=1;
      else accept=0;
    }
    
    if(accept==1){ *fz=fnew; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew; }
    else W[k]=0;
    total_mut++;
  }
  else if(B>=2){
    z[k]=gasdev()*sqrt(var)+mean;
    z[k2]=gasdev()*sqrt(var)+mean;
    fnew=cost(z,W,&aicnew,&bicnew,&ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    
    r=-1.0*(fnew-(*fz))/tem+log(2.0/L)-log(1.0/m*(1.0/((P)+1)+1.0/(OUT_UNIT)))
      -dloggauss(z[k],mean,var)-dloggauss(z[k2],mean,var);
    r+=log(p2/p1);
    
    if(r>0.0) accept=1;
    else{
      un=0.0;
      while(un<=0.0) un=unif_rand();
      if(un<exp(r)) accept=1;
      else accept=0;
    }
    
    if(accept==1){ *fz=fnew; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew; }
    else { W[k]=0; W[k2]=0; }
    total_mut++;
  }
  
  
  /*
  if(check(W)==1){ Rprintf("network adding error happened\n"); exit; }
  */
  
  free_ivector(D,1,dim); 
  free_imatrix(node,1,(hd1),1,2);
  return 0;
}



int del_connection(double *z,double *fz,double tem,int *W,double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  int i,j,k,m,s,L,accept,k2=1,A,B,*D;
  double fnew,un,r,mean,var,aicnew,bicnew, ebicnew;
  double p1, p2;
  int cold, cnew, **node;
  
  
  
  /*
  if(check(W)==1){ Rprintf("network starting error del\n"); exit; }
  */
  
  for(cold=0, i=1; i<=dim; i++) cold+=W[i];
  
  D=ivector(1,dim);
  node=imatrix(1,(hd1),1,2);
  
  for(i=1; i<=(hd1); i++){
    for(L=0, s=(i-1)*((P)+1)+1; s<=i*((P)+1); s++) L+=W[s];
    node[i][1]=L;
    for(L=0,s=((P)+1)*(hd1)+1+i;s<=((P)+1)*(hd1)+(1+(hd1)+(P))*((OUT_UNIT)-1)+1+i;s+=1+(hd1)+(P))
      L+=W[s];
    node[i][2]=L;
    if((node[i][1]==0 && node[i][2]>0)||(node[i][1]>0 && node[i][2]==0)){
      Rprintf("node error mutation adding in del\n"); return 1; }
  }
  
  
  for(j=0, i=1; i<=((P)+1)*(hd1); i++)
    if(W[i]==1){ 
      if(i%((P)+1)==0) A=i/((P)+1);
      else A=(int)(1.0*i/((P)+1))+1;
      if(node[A][1]==1 && node[A][2]>1) A=0;
      else{ j++; D[j]=i; }
    }
    
    for(B=1; B<=(OUT_UNIT); B++){
      i=((P)+1)*(hd1)+(B-1)*(1+(hd1)+(P))+1;
      if(W[i]==1){ j++; D[j]=i; }
      
      for(i=((P)+1)*(hd1)+(B-1)*(1+(hd1)+(P))+2; i<=((P)+1)*(hd1)+(B-1)*(1+(hd1)+(P))+1+(hd1); i++)
        if(W[i]==1){
          A=i-((P)+1)*(hd1)-(B-1)*(1+(hd1)+(P))-1;
          if(node[A][2]>1){ j++; D[j]=i; }
          else if(node[A][2]==1 && node[A][1]==1){ j++; D[j]=i; }
        }
        
        for(i=((P)+1)*(hd1)+(B-1)*(1+(hd1)+(P))+1+(hd1)+1;i<=((P)+1)*(hd1)+(B-1)*(1+(hd1)+(P))+1+(hd1)+(P);i++)
          if(W[i]==1){ j++; D[j]=i; }
    }
    m=j;  /* the number of connections deletable */
  
  /* determine the connection to delete: k */
  L=m+1;
  while(L>m) L=(int)(unif_rand()*m)+1;
  k=D[L];
  if(W[k]==0) Rprintf("Mutation deleting error 1\n");
  
  
  if(k>((P)+1)*(hd1)){
    B=0; i=(OUT_UNIT);
    while(i>=1 && B==0){
      if(k>((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1)+1+(hd1)) B=1;
      else if(k>((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1)+1){
        A=k-((P)+1)*(hd1)-(1+(hd1)+(P))*(i-1)-1;
        if(node[A][2]>1) B=1;
        else if(node[A][2]==1 && node[A][1]==1) B=2;
      }
      else if(k>((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1)) B=1;
      i--;
    }
  }
  else{
    if(k%((P)+1)==0) A=k/((P)+1);
    else A=(int)(1.0*k/((P)+1))+1;
    if(A>(hd1)) A=(hd1);
    if(node[A][2]==1 && node[A][1]==1) B=3;
    else B=1;
  }
  
  if(B==1) W[k]=0;
  else if(B==2){
    A=(k-((P)+1)*(hd1))%(1+(hd1)+(P))-1;
    i=(A-1)*((P)+1)+1;
    while(W[i]==0 && i<=A*((P)+1)) i++;
    k2=i;
    if(W[k2]==0){ Rprintf("Mutation delete error 3\n"); return 1; }
    W[k]=0; W[k2]=0;
  }
  else if(B==3){
    if(k%((P)+1)==0) A=k/((P)+1);
    else A=(int)(1.0*k/((P)+1))+1;
    
    i=((P)+1)*(hd1)+1+A;
    while(W[i]==0 && i<=((P)+1)*(hd1)+((OUT_UNIT)-1)*(1+(hd1)+(P))+1+A) i+=1+(hd1)+(P);
    k2=i;
    if(W[k2]==0){ Rprintf("Mutation delete error 4\n"); return 1; }
    W[k]=0; W[k2]=0;
  }
  
  
  for(cnew=0, i=1; i<=dim; i++) cnew+=W[i];
  if(cnew<3){ accept=0;  total_mut++; 
  if(B==1) W[k]=1;
  else{ W[k]=1; W[k2]=1; }
  free_ivector(D,1,dim);
  free_imatrix(node,1,(hd1),1,2);
  return 0;
  }
  
  if(cold==dim) p1=2.0/3;
  else if(cold==3) p1=0.0;
  else p1=1.0/3;
  if(cnew==3) p2=2.0/3;
  else if(cnew==dim) p2=0.0;
  else p2=1.0/3;
  
  mean=0.0; var=0.0; L=0;
  for(i=1; i<=dim; i++){
    mean+=z[i]*W[i];
    var+=z[i]*z[i]*W[i];
    L+=W[i];
  }
  mean/=L;
  var=(var/L-mean*mean)*L/(L-1);
  
  if(B==1){
    fnew=cost(z,W,&aicnew,&bicnew,&ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    r=-1.0*(fnew-(*fz))/tem-log(1.0*(dim-L))+dloggauss(z[k],mean,var)+log(1.0*m);
    r+=log(p2/p1);
    
    if(r>0.0) accept=1;
    else{
      un=0.0;
      while(un<=0.0) un=unif_rand();
      if(un<exp(r)) accept=1;
      else accept=0;
    }
    
    if(accept==1){ *fz=fnew; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew; }
    else W[k]=1;
    total_mut++;
  }
  else{
    fnew=cost(z,W,&aicnew,&bicnew,&ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    r=-1.0*(fnew-(*fz))/tem+log(1.0/(dim-L)*(1.0/(OUT_UNIT)+1.0/((P)+1)))
      +dloggauss(z[k],mean,var)+dloggauss(z[k2],mean,var)-log(2.0/m);
    r+=log(p2/p1);
    
    if(r>0.0) accept=1;
    else{
      un=0.0;
      while(un<=0.0) un=unif_rand();
      if(un<exp(r)) accept=1;
      else accept=0;
    }
    
    if(accept==1){ *fz=fnew; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew; }
    else { W[k]=1; W[k2]=1; }
    total_mut++;
  }
  
  
  /*
  if(check(W)==1){ Rprintf("network end error del\n"); exit; }
  */
  
  
  free_ivector(D,1,dim);
  free_imatrix(node,1,(hd1),1,2);
  return 0;
}


int mutationMH(double *x,double *fvalue,double t,int *WW, int state, double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)

{
  int i, m;
  double un, fx, aicm, bicm, ebicm;
  
  
  fx=*fvalue; aicm=*aic; bicm=*bic; ebicm=*ebic;
  if(state==0) Metropolis_mut(x,&fx,t,WW,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
  else{
    for(m=0, i=1; i<=dim; i++) m+=WW[i];
    
    un=unif_rand();
    if(m<=3){
      if(un<2.0/3) add_connection(x,&fx,t,WW,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else Metropolis_mut(x,&fx,t,WW,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
    }
    else if(m==dim){
      if(un<2.0/3) del_connection(x,&fx,t,WW,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else Metropolis_mut(x,&fx,t,WW,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
    }
    else{
      if(un<1.0/3) 
        add_connection(x,&fx,t,WW,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else if(un<2.0/3) 
        del_connection(x,&fx,t,WW,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else Metropolis_mut(x,&fx,t,WW,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
    }
  } 
  *fvalue=fx; *aic=aicm; *bic=bicm; *ebic=ebicm;
  
  return 0;
}

// end:mutationMH.c

//begin:  "mutationSAMCg.c"
int Metropolis_mut_SAMC(double *z,double *fz,double tem,int *W,int *region,double **hist,double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  double *y, *d, fnew, r,un,aicnew,bicnew, ebicnew;
  int newregion,i, k,accept, blocksize, tail, m1, m2;
  
  /*
  if(check(W)==1){Rprintf("network simution error 1\n"); exit; }
  */
  
  // Rprintf("mut\n");
  
  grid=ceil((maxE-lowE)*scale);
  y=dvector(1,dim);  d=dvector(1,dim);
  
  for(i=1; i<=dim; i++) y[i]=z[i];
  for(tail=0, i=1; i<=dim; i++) tail+=W[i];
  
  k=0; 
  while(tail>0){
    
    blocksize=BLOCKSIZE;
    if(tail<blocksize) blocksize=tail;
    
    uniform_direction(d,blocksize);
    un=gasdev()*sqrt(tem)*mut_step; /* 2 6.0 9.0 */
  
  m1=k+1; 
  for(i=1; i<=blocksize; i++){
    k++;
    while(k<=dim && W[k]==0) k++;
    if(k<=dim && W[k]==1) y[k]=z[k]+un*d[i];
  }
  m2=k;
  
  fnew=cost(y, W,&aicnew,&bicnew, &ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
  if(fnew>maxE) newregion=grid;
  else if(fnew<lowE) newregion=0;
  else newregion=floor((fnew-lowE)*scale);
  
  r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem;
  
  if(newregion>=grid) accept=0;
  else{
    if(r>0.0) accept=1;
    else{
      un=0.0;
      while(un<=0.0) un=unif_rand();
      if(un<exp(r)) accept=1;
      else accept=0;
    }
  }
  
  if(accept==1){
    for(i=m1; i<=m2; i++) z[i]=y[i];
    *fz=fnew; *region=newregion; *aic=aicnew; *bic=bicnew; *ebic=ebicnew;
    accept_mut++;
  }
  else{
    for(i=m1; i<=m2; i++) y[i]=z[i];
  }
  total_mut++;
  
  tail-=blocksize;
  }
  
  free_dvector(y,1,dim);  free_dvector(d,1,dim);
  /*
  if(check(W)==1){ Rprintf("network simulation error 2\n"); exit; }
  */
  
  // Rprintf("endmut\n");
  
  return 0;
}



int add_connection_SAMC(double *z,double *fz,double tem,int *W,int *region,double **hist,double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  int i,j,m,s,H,k,K,L,d,accept,k1,k2,**node,AA,B;
  double fnew, un, mean=0.0, var=1.0, r,aicnew,bicnew, ebicnew;
  double p1, p2, pxy, pyx;
  int newregion,cold, cnew;
  
  
  grid=ceil((maxE-lowE)*scale);
  node=imatrix(1,(hd1)+(OUT_UNIT),1,2);
  for(cold=0, i=1; i<=dim; i++) cold+=W[i];
  
  // sampling H, if(H=0) update shortcut connections  
  if(shortcut==1){
    H=(hd1)+(OUT_UNIT)+1;
    while(H>(hd1)+(OUT_UNIT)) H=floor(unif_rand()*((hd1)+(OUT_UNIT)))+1; 
  }  
  else{ 
    H=(hd1)+1;
    while(H>(hd1)) H=floor(unif_rand()*(hd1))+1;
  } 
  
  // Rprintf("H=%d\n", H);
  
  if(H>(hd1)){
    m=((P)+1)*(hd1)+(1+(hd1)+(P))*(H-(hd1)-1);
    L=W[m+1];
    for(s=1; s<=(P); s++) L+=W[m+1+(hd1)+s]; 
    node[H][1]=L;  node[H][2]=0; 
    
    if(node[H][1]==(P)+1){ total_mut++; free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 0; } 
    else{
      
      B=1;  // type 1 addiation 
      K=((P)+1)-node[H][1];
      k=K+1;
      while(k>K) k=floor(unif_rand()*K)+1;
      
      if(W[m+1]==0) j=1;
      else j=0; 
      
      // determine the connection to add 
      if(j==k) k1=m+1;
      else{ 
        s=m+1+(hd1); 
        while(j<k && s<m+1+(hd1)+(P)){
          s++;
          if(W[s]==0) j++;
        } 
        k1=s; 
      } 
      W[k1]=1;  
      
      // calculation of proposal probability 
      pxy=1.0/K;  pyx=1.0/(node[H][1]+1); 
    }
  } 
  else{
    
    for(L=0, s=((P)+1)*(H-1)+1; s<=((P)+1)*H; s++) L+=W[s];
    node[H][1]=L;
    
    if(shortcut==1){
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1)+(P))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)+(P)) L+=W[s]; }
    else{
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)) L+=W[s]; 
    }
    node[H][2]=L; 
    
    if((node[H][1]>0 && node[H][2]==0) || (node[H][1]==0 && node[H][2]>0))
    { Rprintf("network structure error\n");  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 1; }
    else{ //no structure error 
      
      if(node[H][1]==(P)+1 && node[H][2]==(OUT_UNIT)){ total_mut++; free_imatrix(node,1, (hd1)+(OUT_UNIT),1,2); return 0; } 
      else if(node[H][2]==(OUT_UNIT) && node[H][1]<(P)+1){
        B=1; // type 1 addition
        
        K=((P)+1)-node[H][1];
        k=K+1;
        while(k>K) k=floor(unif_rand()*K)+1;
        
        s=((P)+1)*(H-1); j=0;
        while(j<k && s<((P)+1)*H){
          s++;
          if(W[s]==0) j++;
        }
        k1=s;
        W[k1]=1;
        
        // calculation of proposal probability 
        pxy=1.0/K;  
        if((OUT_UNIT)==1) pyx=1.0/(node[H][1]+1); 
        else pyx=0.5/(node[H][1]+1); 
      }
      else if(node[H][1]==(P)+1 && node[H][2]<(OUT_UNIT)){ 
        
        B=1; // type 1 addition 
        
        K=(OUT_UNIT)-node[H][2];
        k=K+1;
        while(k>K) k=floor(unif_rand()*K)+1;
        
        j=0; d=0;
        while(j<k && d<(OUT_UNIT)){
          
          d++;
          if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H; 
          else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H; 
          
          if(W[s]==0) j++; 
        } 
        k1=s;
        W[k1]=1;
        
        // calculation of proposal probability 
        pxy=1.0/K;  pyx=0.5/(node[H][2]+1); 
      }
      else if(node[H][1]==0 && node[H][2]==0){
        B=2; // type 2 addition 
        
        k=(P)+2;
        while(k>(P)+1) k=floor(unif_rand()*((P)+1))+1;
        k1=((P)+1)*(H-1)+k;
        
        if((OUT_UNIT)==1) k2=((P)+1)*(hd1)+1+H;
        else{
          k=(OUT_UNIT)+1;
          while(k>(OUT_UNIT)) k=floor(unif_rand()*(OUT_UNIT))+1;
          
          if(shortcut==1) k2=((P)+1)*(hd1)+(1+(hd1)+(P))*(k-1)+1+H;
          else k2=((P)+1)*(hd1)+(1+(hd1))*(k-1)+1+H;
        }
        W[k1]=1; W[k2]=1;
        
        // calculation of proposal probability 
        pxy=1.0/((P)+1)/(OUT_UNIT);  pyx=1.0;
      }
      else{ // case node[H][1]<(P)+1 && node[H][2]<(OUT_UNIT)  
        
        B=1; // type 1 addition 
        
        un=unif_rand();
        if(un<0.5) AA=1;
        else AA=2; 
        
        if(AA==1){ 
          K=((P)+1)-node[H][1];
          k=K+1;
          while(k>K) k=floor(unif_rand()*K)+1;
          
          s=((P)+1)*(H-1); j=0;
          while(j<k && s<((P)+1)*H){
            s++;
            if(W[s]==0) j++;
          }
          k1=s; W[k1]=1;
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/(node[H][1]+1); 
        } 
        else{
          
          K=(OUT_UNIT)-node[H][2];
          k=K+1;
          while(k>K) k=floor(unif_rand()*K)+1;
          
          j=0; d=0;
          while(j<k && d<(OUT_UNIT)){
            d++;
            if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
            else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
            if(W[s]==0) j++;
          }
          k1=s; W[k1]=1; 
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/(node[H][2]+1);
        } 
      }
    } //end no structure error 
  } // end of H is not equal to 0 
  
  
  if(B==1) cnew=cold+1;
  else if(B==2) cnew=cold+2;
  if(cold==1) p1=2.0/3;
  else if(cold==connection_threshold) p1=0.0;
  else p1=1.0/3;
  if(cnew==connection_threshold) p2=2.0/3;
  else if(cnew==1) p2=0.0;
  else p2=1.0/3;
  
  
  if(B==1){
    z[k1]=gasdev()*sqrt(var)+mean;
    
    fnew=cost(z,W,&aicnew,&bicnew,&ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    if(fnew>maxE) newregion=grid;
    else if(fnew<lowE) newregion=0;
    else newregion=floor((fnew-lowE)*scale);
    
    if(newregion>=grid) accept=0;
    else{
      
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem-dloggauss(z[k],mean,var);
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew; }
    else W[k1]=0;
    total_mut++;
  }
  else if(B==2){
    z[k1]=gasdev()*sqrt(var)+mean;
    z[k2]=gasdev()*sqrt(var)+mean;
    
    fnew=cost(z,W,&aicnew,&bicnew, &ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    
    if(fnew>maxE) newregion=grid;
    else if(fnew<lowE) newregion=0;
    else newregion=floor((fnew-lowE)*scale);
    
    if(newregion>=grid) accept=0;
    else{
      
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem-dloggauss(z[k1],mean,var)-dloggauss(z[k2],mean,var);
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew; }
    else { W[k1]=0; W[k2]=0; }
    total_mut++;
  }
  
  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2);
  
  return 0;
}



int del_connection_SAMC(double *z,double *fz,double tem,int *W,int *region,double **hist,double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  int i,j,k,K,m,s,d,H,L,accept,k1=1,k2,AA,B=1;
  double fnew,un,r,mean=0.0,var=1.0,aicnew,bicnew, ebicnew;
  double p1, p2, pxy=1.0, pyx=1.0;
  int newregion,cold, cnew, **node;
  
  
  
  grid=ceil((maxE-lowE)*scale);
  node=imatrix(1,(hd1)+(OUT_UNIT),1,2);
  for(cold=0, i=1; i<=dim; i++) cold+=W[i];
  
  // sampling H, if(H=0) update shortcut connections  
  if(shortcut==1){
    H=(hd1)+(OUT_UNIT)+1;
    while(H>(hd1)+(OUT_UNIT)) H=floor(unif_rand()*((hd1)+(OUT_UNIT)))+1;
  }
  else{
    H=(hd1)+1;
    while(H>(hd1)) H=floor(unif_rand()*(hd1))+1;
  }
  
  // Rprintf("del H=%d\n", H);
  
  if(H>(hd1)){
    m=((P)+1)*(hd1)+(1+(hd1)+(P))*(H-(hd1)-1);
    L=W[m+1];
    for(s=1; s<=(P); s++) L+=W[m+1+(hd1)+s];
    node[H][1]=L;  node[H][2]=0;
    
    if(node[H][1]==0){ total_mut++; free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 0; }
    else{
      
      B=1;  // type 1 addiation 
      K=node[H][1];
      k=K+1;
      while(k>K) k=floor(unif_rand()*K)+1;
      
      if(W[m+1]==1) j=1;
      else j=0;
      
      // determine the connection to add 
      if(j==k) k1=m+1;
      else{
        s=m+1+(hd1);
        while(j<k && s<m+1+(hd1)+(P)){
          s++;
          if(W[s]==1) j++;
        }
        k1=s;
      }
      W[k1]=0;
      
      // calculation of proposal probability 
      pxy=1.0/K;  pyx=1.0/(((P)+1)-(node[H][1]-1));
    }
  } 
  else{ // H is not equal to 0 
    
    for(L=0, s=((P)+1)*(H-1)+1; s<=((P)+1)*H; s++) L+=W[s];
    node[H][1]=L;
    
    if(shortcut==1){
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1)+(P))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)+(P)) L+=W[s]; }
    else{
      for(L=0,s=((P)+1)*(hd1)+1+H; s<=((P)+1)*(hd1)+(1+(hd1))*((OUT_UNIT)-1)+1+H; s+=1+(hd1)) L+=W[s];
    }
    node[H][2]=L; 
    
    if((node[H][1]>0 && node[H][2]==0) || (node[H][1]==0 && node[H][2]>0))
    { Rprintf("network structure error\n");  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 1; }
    else{ //no structure error 
      
      if(node[H][1]==0 && node[H][2]==0){ total_mut++; free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2); return 0; }
      else if(node[H][1]==1 && node[H][2]==1){ 
        
        B=2; // type 2 deletion 
        
        s=((P)+1)*(H-1)+1; j=1;
        while(W[s]==0 && j<=(P)+1){ s++; j++; }
        if(j==(P)+2) Rprintf("network error\n"); 
        k1=s; 
        
        
        if((OUT_UNIT)==1) k2=((P)+1)*(hd1)+1+H;
        else{
          
          d=1; s=((P)+1)*(hd1)+1+H;
          while(W[s]==0 && d<(OUT_UNIT)){
            d++;
            if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
            else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
          }
          k2=s; 
        }
        
        W[k1]=0; W[k2]=0; 
        
        
        // calculation of proposal probability 
        pxy=1.0; pyx=1.0/((P)+1)/(OUT_UNIT);   
      } 
      else if(node[H][1]==1 && node[H][2]>1){ 
        
        B=1; //type 1 deletion 
        
        K=node[H][2];
        k=K+1;
        while(k>K) k=floor(unif_rand()*K)+1;
        
        j=0; d=0;
        while(j<k && d<(OUT_UNIT)){
          d++;
          if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
          else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
          
          if(W[s]==1) j++;
        }
        k1=s; W[k1]=0;
        
        // calculation of proposal probability 
        pxy=1.0/K;  pyx=0.5/((OUT_UNIT)-(node[H][2]-1)); 
      }
      else if(node[H][1]>1 && node[H][2]==1){
        
        B=1; 
        K=node[H][1];
        k=K+1;
        while(k>K) k=floor(unif_rand()*K)+1;
        
        s=((P)+1)*(H-1); j=0;
        while(j<k && s<((P)+1)*H){
          s++;
          if(W[s]==1) j++;
        }
        k1=s;  W[k1]=0;
        
        // calculation of proposal probability 
        pxy=1.0/K;  
        if((OUT_UNIT)==1) pyx=1.0/(((P)+1)-(node[H][1]-1)); 
        else pyx=0.5/(((P)+1)-(node[H][1]-1)); 
      }
      else if(node[H][1]>1 && node[H][2]>1){ 
        
        B=1; 
        
        un=unif_rand();
        if(un<0.5) AA=1;
        else AA=2;
        
        if(AA==1){
          K=node[H][1];
          k=K+1;
          while(k>K) k=floor(unif_rand()*K)+1;
          
          s=((P)+1)*(H-1); j=0;
          while(j<k && s<((P)+1)*H){
            s++;
            if(W[s]==1) j++;
          }
          k1=s; W[k1]=0;
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/(((P)+1)-(node[H][1]-1));
        }
        else{ 
          K=node[H][2];
          k=K+1;
          while(k>K) k=floor(unif_rand()*K)+1;
          
          j=0; d=0;
          while(j<k && d<(OUT_UNIT)){
            d++;
            if(shortcut==1) s=(1+(P))*(hd1)+(1+(hd1)+(P))*(d-1)+1+H;
            else s=(1+(P))*(hd1)+(1+(hd1))*(d-1)+1+H;
            if(W[s]==1) j++;
          }
          k1=s; W[k1]=0;
          
          // calculation of proposal probability 
          pxy=0.5/K;  pyx=0.5/((OUT_UNIT)-(node[H][2]-1));
        }
      }
    } // end of no structure error 
  } // end H is not equal to 0 
  
  
  if(B==1) cnew=cold-1;
  else if(B==2) cnew=cold-2;
  
  if(cold==connection_threshold) p1=2.0/3;
  else if(cold==1) p1=0.0;
  else p1=1.0/3;
  if(cnew==1) p2=2.0/3;
  else if(cnew==connection_threshold) p2=0.0;
  else p2=1.0/3;
  
  
  if(B==1){
    
    fnew=cost(z,W,&aicnew, &bicnew, &ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    if(fnew>maxE) newregion=grid;
    else if(fnew<lowE) newregion=0;
    else{
      newregion=floor((fnew-lowE)*scale);
    }
    
    if(newregion>=grid || cnew==0) accept=0;
    else{
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem+dloggauss(z[k1],mean,var);
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew; }
    else W[k1]=1;
    total_mut++;
  }
  else{
    
    fnew=cost(z,W,&aicnew,&bicnew, &ebicnew, hd1, OUT_UNIT, P, lambda, data_num);
    if(fnew>maxE) newregion=grid;
    else if(fnew<lowE) newregion=0;
    else newregion=floor((fnew-lowE)*scale);
    
    if(newregion>=grid || cnew==0) accept=0;
    else{
      r=hist[*region][2]-hist[newregion][2]-1.0*(fnew-(*fz))/tem+dloggauss(z[k1],mean,var)+dloggauss(z[k2],mean,var);
      r+=log(p2/p1)+log(pyx/pxy);
      
      if(r>0.0) accept=1;
      else{
        un=0.0;
        while(un<=0.0) un=unif_rand();
        if(un<exp(r)) accept=1;
        else accept=0;
      }
    }
    
    if(accept==1){ *fz=fnew; *region=newregion; accept_mut++; *aic=aicnew; *bic=bicnew; *ebic=ebicnew;  }
    else { W[k1]=1; W[k2]=1; }
    total_mut++;
  }
  
  
  free_imatrix(node,1,(hd1)+(OUT_UNIT),1,2);
  
  return 0;
}

int mutation_SAMC(double *x,double *fvalue,double t,int *WW,int *pos,double **hist,int state,double *aic,double *bic,double *ebic, int hd1, int OUT_UNIT, int P, double lambda, int data_num)
{
  int i, m, region;
  double un, fx, aicm, bicm, ebicm;
  
  
  fx=*fvalue; aicm=*aic; bicm=*bic; ebicm=*ebic;
  region=*pos;
  // Rprintf("fx=%g region=%d\n",fx,region);
  
  if(state==0) Metropolis_mut_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
  else{
    for(m=0, i=1; i<=dim; i++) m+=WW[i];
    
    un=unif_rand();
    if(m<=1){
      if(un<2.0/3) add_connection_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else Metropolis_mut_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
    }
    else if(m==connection_threshold){
      if(un<2.0/3) del_connection_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else Metropolis_mut_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
    }
    else{
      if(un<1.0/3) 
        add_connection_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else if(un<2.0/3) 
        del_connection_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm,&ebicm, hd1, OUT_UNIT, P, lambda, data_num);
      else Metropolis_mut_SAMC(x,&fx,t,WW,&region,hist,&aicm,&bicm, &ebicm, hd1, OUT_UNIT, P, lambda, data_num);
    }
  }
  
  *fvalue=fx; *aic=aicm; *bic=bicm; *ebic=ebicm;
  *pos=region; 
  
  // Rprintf("****\n");
  
  return 0;
}

// end: mutationSAMCg.c


//begin: "fitting.c"
int fitting(double *ox,double *z,int *W,double *ps2, int hd1, int OUT_UNIT, int P)
{
  int i, j, k;
  double ps1[(hd1)+1], ave; 
  
  /* calculate the output of the first hidden layer */
  for(i=1; i<=(hd1); i++){
    k=((P)+1)*(i-1);
    ps1[i]=z[k+1]*W[k+1];
    for(j=k+2; j<=k+(P)+1; j++) ps1[i]+=ox[j-k-1]*z[j]*W[j];
    ps1[i]=tanh(ps1[i]);
  }
  
  /* calculate the predicted mean from the mean unit */
  if(shortcut==0){
    for(i=1; i<=(OUT_UNIT); i++){
      k=((P)+1)*(hd1)+(1+(hd1))*(i-1);
      ave=z[k+1]*W[k+1];
      for(j=k+2; j<=k+(hd1)+1; j++) ave+=ps1[j-k-1]*z[j]*W[j];
      ps2[i]=ave;
    }
  }
  else{
    for(i=1; i<=(OUT_UNIT); i++){
      k=((P)+1)*(hd1)+(1+(hd1)+(P))*(i-1);
      ave=z[k+1]*W[k+1];
      for(j=k+2; j<=k+(hd1)+1; j++) ave+=ps1[j-k-1]*z[j]*W[j];
      for(j=k+(hd1)+2; j<=k+(P)+(hd1)+1; j++) ave+=ox[j-k-(hd1)-1]*z[j]*W[j];
      ps2[i]=ave;
    }
  }
  
  return 0;
}
//end: fitting.c


void posratio(double DataX[], double DataY[], int *data_num, int *test_num, int *OUT_UNIT, int *P, int *hd1, double *lambda, int *total_iteration, int *popN, int *nCPUs)
{
  //R_CStackLimit=(uintptr_t)-1;
  int inipopN = (*popN);
  
  
#ifdef _OPENMP
  omp_set_num_threads(*nCPUs);
#endif
  
  GetRNGstate();
  
  
  
  double **x,**sx,*y,*fvalue,*sfvalue,*fx,*fy,*t,**total_mat,*ps2,**FIT,**PRED,**total_y,**hist;
  int count,grid0,maxpopN,i, j, k0, k, iter, m, s1,dim0,bestk=1,newgrid,**WW,**sWW,*W,*rem,*sel,*pos,*indx;
  double q,w,bias,predbias,mis1,mis2, linearity;
  double *mar, *prob, total_weight, *net, *locfreq, energy, ave, max, sum, maxhist;
  double  *AIC, *BIC, *EBIC, *sAIC, *sBIC, *aAIC, *aBIC, logIAIC, logIBIC, KLAIC, KLBIC, KLEBIC;
  FILE *ins;
  
  int warm0=ceil((double)*total_iteration/40.0),warm1=ceil((double)*total_iteration/40.0),warm2=ceil((double)*total_iteration/5.0), div=ceil((double)*total_iteration/10000.0),step=ceil((double)*total_iteration/10000.0),stepscale=ceil((double)*total_iteration/50.0);
  /*
  ins=fopen("aa.log", "a");
  fprintf(ins, "\nseed=%d HD=%d (*P)=%d\n", stime,(*hd1),(*P));
  fclose(ins);
  */
  if(shortcut==0){ dim=((*P)+1)*(*hd1)+(1+(*hd1))*(*OUT_UNIT); dim0=((*P)+1)*(*hd1); }
  else{  
    dim0=((*P)+1)*(*hd1); 
    dim=((*P)+1)*(*hd1)+(1+(*hd1)+(*P))*(*OUT_UNIT);
  }  
  if(connection_threshold>dim) connection_threshold=dim;
  
  if(inipopN>(*popN)) maxpopN=inipopN; 
  else maxpopN=(*popN);
  
  x=dmatrix(1,maxpopN,1,dim);
  sx=dmatrix(1,maxpopN,1,dim);
  y=dvector(1,dim);
  net=dvector(1,dim);
  WW=imatrix(1,maxpopN,1,dim);
  sWW=imatrix(1,maxpopN,1,dim);
  W=ivector(1,dim);
  fvalue=dvector(1,maxpopN);
  sfvalue=dvector(1,maxpopN);
  t=dvector(1,maxpopN);
  total_mat=dmatrix(1,((*data_num)+(*test_num)),1,(*P));
  total_y=dmatrix(1,((*data_num)+(*test_num)),1,(*OUT_UNIT));
  data_mat=dmatrix(1,(*data_num),1,(*P));
  test_mat=dmatrix(1,(*test_num),1,(*P));
  data_y=dmatrix(1,(*data_num),1,(*OUT_UNIT));
  test_y=dmatrix(1,(*test_num),1,(*OUT_UNIT));
  ps2=dvector(1,(*OUT_UNIT));
  sel=ivector(1,(*data_num));
  rem=ivector(1,((*data_num)+(*test_num)));
  FIT=dmatrix(1,(*data_num),1,(*OUT_UNIT));
  PRED=dmatrix(1,(*test_num),1,(*OUT_UNIT));
  wei=dvector(1,(*data_num));
  indx=ivector(1,1000);
  fx=dvector(1,1000);
  fy=dvector(1,1000);
  mar=dvector(1,(*P)+1);
  prob=dvector(1,(*P)+1);
  AIC=dvector(1,maxpopN);
  BIC=dvector(1,maxpopN);
  EBIC=dvector(1,maxpopN);
  sAIC=dvector(1,(*popN)*(*total_iteration));
  sBIC=dvector(1,(*popN)*(*total_iteration));
  aAIC=dvector(1,(*popN)*(*total_iteration));
  aBIC=dvector(1,(*popN)*(*total_iteration));
  
  
  t[1]=hightem;  t[(*popN)]=lowtem;
  q=(1.0/lowtem-1.0/hightem)/((*popN)-1);
  for(i=2; i<=(*popN)-1; i++) t[i]=1.0/(q+1.0/t[i-1]);
  
  
  /*
  if(v==0) lambda=0.01;
  else if(v==1) lambda=0.05;
  else if(v==2) lambda=0.1;
  else lambda=0.2;
  */
  /*
  ins=fopen("aa.log", "a");
  fprintf(ins, "lambda=%g\n", lambda);
  fclose(ins);
  */
  /* read the data from DataX, DataY */
  for(i=1; i<=(*data_num)+(*test_num); i++){
    for(j=1; j<=(*OUT_UNIT); j++){
      total_y[i][j] = DataY[(i-1)*(*OUT_UNIT)+(j-1)];
      //total_y[i][j] = sqrt(total_y[i][j]);
    }
    for(j=1; j<=(*P); j++){
      total_mat[i][j] = DataX[(i-1)*(*P)+(j-1)];
    } 
  }
  /*
  if(v==0){ 
  for(i=1; i<=(*test_num); i++) rem[i]=i;
  for(i=1; i<=(*data_num); i++) sel[i]=(*test_num)+i;
  }
  else if(v==1){ 
  for(i=1; i<=(*test_num); i++) sel[i]=i;
  for(i=(*test_num)+1; i<=2*(*test_num); i++) rem[i-(*test_num)]=i;
  for(i=2*(*test_num)+1; i<=(*data_num)+(*test_num); i++) sel[i-(*test_num)]=i; 
  }
  else if(v==2){
  for(i=1; i<=2*(*test_num); i++) sel[i]=i;
  for(i=2*(*test_num)+1; i<=3*(*test_num); i++) rem[i-2*(*test_num)]=i;
  for(i=3*(*test_num)+1; i<=(*data_num)+(*test_num); i++) sel[i-(*test_num)]=i; 
  }
  else if(v==3){
  for(i=1; i<=3*(*test_num); i++) sel[i]=i;
  for(i=3*(*test_num)+1; i<=4*(*test_num); i++) rem[i-3*(*test_num)]=i;
  for(i=4*(*test_num)+1; i<=(*data_num)+(*test_num); i++) sel[i-(*test_num)]=i;
  }
  else{
  for(i=1; i<=(*data_num); i++) sel[i]=i;
  for(i=1; i<=(*test_num); i++) rem[i]=(*data_num)+i;
  }
  */
  
  for(i=1; i<=(*data_num); i++){
    for(j=1; j<=(*OUT_UNIT); j++) data_y[i][j]=total_y[i][j];
    for(j=1; j<=(*P); j++) data_mat[i][j]=total_mat[i][j];
  }
  
  for(i=1; i<=(*test_num); i++){
    for(j=1; j<=(*OUT_UNIT); j++) test_y[i][j]=total_y[i+(*data_num)][j];
    for(j=1; j<=(*P); j++) test_mat[i][j]=total_mat[i+(*data_num)][j];
  }
  for(i=1; i<=(*data_num); i++) wei[i]=1.0;
  
  
#pragma omp parallel for private(i,j,iter,k0,m) default(shared)
  for(i=1; i<=inipopN; i++){
    
    for(j=1; j<=dim; j++) WW[i][j]=0;
    k0=6;
    while (k0>5) k0=floor(unif_rand()*5)+1;
    for(j=1; j<=k0; j++){
      m=(*P)+2;
      while(m>(*P)+1) m=floor(unif_rand()*((*P)+1))+1;
      WW[i][m]=1;
    }
    WW[i][(*hd1)*((*P)+1)+2]=1;
    
    fvalue[i]=1.0e+100;
    while(fvalue[i]>1.0e+99){
      for(j=1; j<=dim; j++) x[i][j]=gasdev()*0.1;
      fvalue[i]=cost(x[i],WW[i],&(AIC[i]), &(BIC[i]), &(EBIC[i]), *hd1, *OUT_UNIT, *P, *lambda, *data_num);
    }
    
    for(iter=1; iter<=warm0; iter++){
      if(iter<warm0/div) mutationMH(x[i],&(fvalue[i]),1.0,WW[i],0, &(AIC[i]), &(BIC[i]), &(EBIC[i]), *hd1, *OUT_UNIT, *P, *lambda, *data_num);
      else mutationMH(x[i],&(fvalue[i]),1.0,WW[i],1, &(AIC[i]), &(BIC[i]), &(EBIC[i]), *hd1, *OUT_UNIT, *P, *lambda, *data_num); 
      // if(iter%10000==0) Rprintf("inipop=%d iter=%d fvalue=%g\n", i, iter, fvalue[i]);
    }
  }
#pragma omp barrier    
  k=inipopN;
  for(i=1; i<=k; i++) fx[i]=fvalue[i];
  indexx(k,fx,indx); 
  k0=floor(k*0.9); 
  if(k0<1) k0=1;
  for(i=1; i<=k0; i++) fy[i]=fx[indx[i]];
  q=standard_deviation(fy,k0);
  lowE=fx[indx[1]]-exL*q;
  bestE=fx[indx[1]];
  maxE=bestE+range; 
  //  Rprintf("lowE=%g  maxE=%g\n", lowE, maxE);
  //  lowE=200.0; maxE=300.0;
  
  /*
  ins=fopen("aa.log","a");
  fprintf(ins, "lowE=%g  maxE=%g\n", lowE, maxE);
  fclose(ins);
  */
  
  grid=ceil((maxE-lowE)*scale);
  hist=dmatrix(0,grid,1,3);
  refden=dvector(0,grid);
  pos=ivector(0,grid);
  locfreq=dvector(0,grid);
  grid0=floor((bestE-lowE)*scale)-extgrid;
  if(grid0<0) grid0=0;
  
  /* Initialize the configurations */
  for(sum=0.0, i=0; i<grid; i++){ refden[i]=1.0; sum+=refden[i]; }
  // for(i=0; i<=grid; i++) refden[i]/=sum;
  for(i=0; i<=grid; i++){
    hist[i][1]=lowE+i*1.0/scale;
    hist[i][2]=0.0;
    hist[i][3]=0.0;
  }
  
  
  // run MH algorithm to get good starting points
  k=0;
  while(k<(*popN)){
    
#pragma omp parallel for private(i,j,k0,m,iter) default(shared)
    for(i=1; i<=(*popN); i++){ 
      
      for(j=1; j<=dim; j++) sWW[i][j]=0;
      k0=6;
      while (k0>5) k0=floor(unif_rand()*5)+1;
      for(j=1; j<=k0; j++){
        m=(*P)+2;
        while(m>(*P)+1) m=floor(unif_rand()*((*P)+1))+1;
        sWW[i][m]=1;
      }
      sWW[i][(*hd1)*((*P)+1)+2]=1;
      
      sfvalue[i]=1.0e+100;
      while(sfvalue[i]>1.0e+99){
        for(j=1; j<=dim; j++) sx[i][j]=gasdev()*0.1;
        sfvalue[i]=cost(sx[i],sWW[i],&(AIC[i]), &(BIC[i]), &(EBIC[i]), *hd1, *OUT_UNIT, *P, *lambda, *data_num);
      }
      
      for(iter=1; iter<=warm1; iter++){
        if(iter<warm1/div) mutationMH(sx[i],&(sfvalue[i]),t[i],sWW[i],0,&(AIC[i]),&(BIC[i]),&(EBIC[i]), *hd1, *OUT_UNIT, *P, *lambda, *data_num);
        else mutationMH(sx[i],&(sfvalue[i]),t[i],sWW[i],1,&(AIC[i]),&(BIC[i]),&(EBIC[i]), *hd1, *OUT_UNIT, *P, *lambda, *data_num); 
        // if(iter%10000==0) Rprintf("pop=%d iter=%d fvalue=%g\n", i, iter, sfvalue[i]);
      }
    } 
#pragma omp barrier      
    for(i=1; i<=(*popN); i++){
      if(k<(*popN) && sfvalue[i]<maxE){
        k++;
        fvalue[k]=sfvalue[i];
        for(j=1; j<=dim; j++){ x[k][j]=sx[i][j]; WW[k][j]=sWW[i][j]; }
        
        if(fvalue[k]>maxE) pos[k]=grid;
        else if(fvalue[k]<lowE) pos[k]=0;
        else pos[k]=floor((fvalue[k]-lowE)*scale);
      }
      if(k>0 && pos[k]>=grid) k--; 
    }
    // Rprintf("Initial points k=%d\n", k);
  }
  
  
  /*
  for(i=1; i<=(*popN); i++){
  Rprintf("pos[%d]=%d, %g\n",i,pos[i],fvalue[i]);
  un=cost(x[i],WW[i],&(AIC[i]),&(BIC[i]));
  Rprintf("un=%g\n", un);
  }
  */
  
  for(j=1; j<=dim; j++) net[j]=0.0;
  for(j=1; j<=(*P)+1; j++) prob[j]=0.0;
  total_weight=0.0; energy=0.0; bias=predbias=0.0;  KLAIC=KLBIC=KLEBIC=0.0;
  accept_mut=0; total_mut=1; linearity=0.0;  count=0;
  for(i=1; i<=(*data_num); i++)
    for(j=1; j<=(*OUT_UNIT); j++) FIT[i][j]=0.0;
  for(i=1; i<=(*test_num); i++)
    for(j=1; j<=(*OUT_UNIT); j++) PRED[i][j]=0.0;
  
  /* 
  instr=fopen("aa", "a");
  if(instr==NULL){ Rprintf("can't write to file\n"); }
  */  
  
  for(iter=1; iter<=warm2+(*total_iteration); iter++){
    
    if(iter<=WARM*stepscale) delta=rho;
    else delta=rho*exp(-tau*log(1.0*(iter-(WARM-1)*stepscale)/stepscale));
    
#pragma omp parallel for private(i) default(shared)
    for(i=1; i<=(*popN); i++){
      mutation_SAMC(x[i],&(fvalue[i]),t[i],WW[i],&(pos[i]),hist,1,&(AIC[i]),&(BIC[i]),&(EBIC[i]), *hd1, *OUT_UNIT, *P, *lambda, *data_num);
      
      // if(iter%1000==0) Rprintf("iter=%d pos[%d]=%d fvalue=%g AIC=%g BIC=%g EBIC=%g\n",iter,i,pos[i],fvalue[i],AIC[i],BIC[i],EBIC[i]);
    }
#pragma omp barrier      
    
    for(j=0; j<grid; j++) locfreq[j]=0;
    for(i=1; i<=(*popN); i++) locfreq[pos[i]]+=1.0; 
    
    for(i=1; i<=(*popN); i++){ if(fvalue[i]<bestE){ bestE=fvalue[i]; bestk=i; }}
    grid0=floor((bestE-lowE)*scale)-extgrid;
    if(grid0<0) grid0=0;
    
    for(j=grid0; j<grid; j++){
      hist[j][2]+=delta*(1.0*locfreq[j]/(*popN)-1.0*refden[j]/(grid-grid0));
      hist[j][3]+=locfreq[j];
    } 
    
    maxE=bestE+range;
    newgrid=ceil((maxE-lowE)*scale);
    if(newgrid<grid){
      for(maxhist=hist[newgrid][2], i=newgrid+1; i<grid; i++){
        if(hist[i][2]>maxhist) maxhist=hist[i][2];
        hist[newgrid][3]+=hist[i][3];
      }
      for(sum=0.0, i=newgrid; i<grid; i++) sum+=exp(hist[i][2]-maxhist);
      hist[newgrid][2]=log(sum)+maxhist;
      
      for(i=1; i<=(*popN); i++){
        if(fvalue[i]>maxE){
          fvalue[i]=fvalue[bestk];
          pos[i]=pos[bestk];
          for(j=1; j<=dim; j++){ WW[i][j]=WW[bestk][j]; x[i][j]=x[bestk][j]; }
        }
      }
      grid=newgrid;
    }
    
    /* weight normalization */
    if(iter==warm2){
      for(sum=0.0,k0=0,i=grid0; i<grid; i++)
        if(hist[i][3]<=0.0){ sum+=1.0*refden[i]/(grid-grid0); k0++; }
        if(k0>0) ave=sum/k0;
        else ave=0.0;
        for(i=grid0; i<grid; i++) hist[i][2]=hist[i][2]+log(refden[i]+ave);
        max=hist[grid0][2];
        for(i=grid0+1; i<grid; i++)
          if(hist[i][2]>max) max=hist[i][2];
          for(sum=0.0, i=grid0; i<grid; i++){ hist[i][2]-=max; sum+=exp(hist[i][2]); }
          for(i=grid0; i<grid; i++) hist[i][2]=hist[i][2]-log(sum)+log(100.0);
          /*
          max=hist[grid0][2];
          for(i=grid0+1; i<grid; i++)
          if(hist[i][2]>max) max=hist[i][2];
          for(i=grid0; i<=grid; i++) hist[i][2]-=max;
          */
    }
    
    
    if(iter>warm2 && iter%step==0){
      for(i=1; i<=(*popN); i++){
        
        if(pos[i]<grid){
          // output the network
          /*
          if(iter%1000==0){ 
          fprintf(instr, "%g %g ", fvalue[i], hist[pos[i]][2]);
          for(j=1; j<=dim; j++){ if(WW[i][j]!=0) fprintf(instr, " %d", j); }
          fprintf(instr, "\n");
          for(j=1; j<=dim; j++){ if(WW[i][j]!=0) fprintf(instr, " %g", x[i][j]); }
          fprintf(instr, "\n");
          }
          */ 
          
          w=exp(hist[pos[i]][2]); total_weight+=w;
          energy+=w*fvalue[i]; 
          
          count++;
          sAIC[count]=AIC[i]+hist[pos[i]][2]; 
          sBIC[count]=BIC[i]+hist[pos[i]][2]; 
          aAIC[count]=-AIC[i]+hist[pos[i]][2];
          aBIC[count]=-BIC[i]+hist[pos[i]][2];
          KLAIC+=-w*AIC[i];
          KLBIC+=-w*BIC[i];
          KLEBIC+=-w*EBIC[i];
          
          
          for(s1=0, j=1; j<=dim0; j++){
            net[j]+=w*WW[i][j];
            if(j%(1+(*P))!=1 && WW[i][j]==1) s1++;
          }
          if(s1==0) linearity+=w;
          for(j=dim0+1; j<=dim; j++) net[j]+=w*WW[i][j]; 
          
          // marginal inclusion prob.
          for(j=1; j<=(*P)+1; j++) mar[j]=0.0;
          for(k=1; k<=(*hd1); k++){
            for(j=(k-1)*((*P)+1)+1; j<=k*((*P)+1); j++) mar[j-(k-1)*((*P)+1)]+=WW[i][j]; }
          if(shortcut==1){
            mar[1]+=WW[i][(*hd1)*((*P)+1)+1];
            for(j=(*hd1)*((*P)+1)+(*hd1)+2; j<=dim; j++) mar[j-(*hd1)*((*P)+1)-(*hd1)]+=WW[i][j];
          }
          for(j=1; j<=(*P)+1; j++){ if(mar[j]>0.5) prob[j]+=w; }
          
          for(mis1=0, k=1; k<=(*data_num); k++){
            fitting(data_mat[k],x[i],WW[i],ps2, *hd1, *OUT_UNIT, *P);
            for(j=1; j<=(*OUT_UNIT); j++){ 
              mis1+=w*fabs(0.5*(data_y[k][j]-ps2[j]));
              FIT[k][j]+=w*ps2[j];
              // Rprintf(" %g %g\n", data_y[k][j], ps2[j]);
            }
          }
          bias+=mis1;
          
          for(mis2=0, k=1; k<=(*test_num); k++){
            fitting(test_mat[k],x[i],WW[i],ps2, *hd1, *OUT_UNIT, *P);
            for(j=1; j<=(*OUT_UNIT); j++){
              mis2+=w*fabs(0.5*(test_y[k][j]-ps2[j]));
              PRED[k][j]+=w*ps2[j];
            }
          }
          predbias+=mis2;
        }
    }
  } /* end for iter>warm */
} /* end for iter */
          //    fclose(instr);
          for(j=1; j<=dim; j++) net[j]/=total_weight;  
  for(j=1; j<=(*P)+1; j++) prob[j]/=total_weight;
  for(k=1; k<=(*data_num); k++) 
    for(j=1; j<=(*OUT_UNIT); j++) FIT[k][j]/=total_weight;
  for(k=1; k<=(*test_num); k++) 
    for(j=1; j<=(*OUT_UNIT); j++) PRED[k][j]/=total_weight;  
  bias/=total_weight; 
  predbias/=total_weight;
  KLAIC/=total_weight; KLBIC/=total_weight; KLEBIC/=total_weight;
  
  max=sAIC[1];
  for(j=2; j<=count; j++){
    if(sAIC[j]>max) max=sAIC[j]; }
  for(sum=0.0, j=1; j<=count; j++)  sum+=exp(sAIC[j]-max);  
  logIAIC=-1.0*(max+log(sum)-log(total_weight)); 
  
  max=sBIC[1];
  for(j=2; j<=count; j++){
    if(sBIC[j]>max) max=sBIC[j]; }
  for(sum=0.0, j=1; j<=count; j++)  sum+=exp(sBIC[j]-max);
  logIBIC=-1.0*(max+log(sum)-log(total_weight));
  
  max=aAIC[1];
  for(j=2; j<=count; j++){
    if(aAIC[j]>max) max=aAIC[j]; }
  for(sum=0.0, j=1; j<=count; j++)  sum+=exp(aAIC[j]-max);
  //aveAIC=max+log(sum)-log(total_weight);
  
  max=aBIC[1];
  for(j=2; j<=count; j++){
    if(aBIC[j]>max) max=aBIC[j]; }
  for(sum=0.0, j=1; j<=count; j++)  sum+=exp(aBIC[j]-max);
  //aveBIC=max+log(sum)-log(total_weight);
  
  KLAIC-=logIAIC; KLBIC-=logIBIC;
  /*
  ins=fopen("aa.s","a");
  fprintf(ins, " %g %g %g %g %g %g %g %g %g\n",logIAIC,logIBIC,KLAIC,KLBIC,KLEBIC,aveAIC,aveBIC,logIAIC+KLAIC,logIBIC+KLBIC);
  fclose(ins);
  
  ins=fopen("aa.log", "a");
  fprintf(ins, "linearity probability=%g\n", linearity/total_weight);
  fprintf(ins, "mutation: %g  ", 1.0*accept_mut/total_mut); 
  fprintf(ins, "bias=%g predbias=%g \n", bias, predbias);
  fprintf(ins, "logIAIC=%g logIBIC=%g KLAIC=%g KLBIC=%g KLEBIC=%g aveAIC=%g aveBIC=%g mAIC=%g mBIC=%g\n",logIAIC,logIBIC,KLAIC,KLBIC,KLEBIC, aveAIC,aveBIC,logIAIC+KLAIC,logIBIC+KLBIC);
  for(i=1; i<=dim; i++) fprintf(ins, " %g", net[i]);
  fprintf(ins, "\n");
  for(i=1; i<=(*P)+1; i++) fprintf(ins, " %g", prob[i]);
  fprintf(ins, "\n");
  fclose(ins); 
  */
  ins=fopen("aa.net","w");
  for(i=1; i<=dim; i++) fprintf(ins, " %g", net[i]);
  fprintf(ins, "\n");
  fclose(ins);
  
  ins=fopen("aa.BF", "w");
  fprintf(ins, " %g %g %g\n", linearity/total_weight, 1-linearity/total_weight, linearity/(total_weight-linearity));
  fclose(ins);
  
  ins=fopen("aa.mar", "w");
  for(i=1; i<=(*P)+1; i++) fprintf(ins, " %g", prob[i]);
  fprintf(ins, "\n\n");
  fclose(ins);
  
  for(sum=0.0,k0=0,i=grid0; i<grid; i++)
    if(hist[i][3]<=0.0){ sum+=1.0*refden[i]/(grid-grid0); k0++; }
    if(k0>0) ave=sum/k0;
    else ave=0.0;
    for(i=grid0; i<grid; i++) hist[i][2]=hist[i][2]+log(refden[i]+ave);
    // ignore the last subregion in normalization
    max=hist[grid0][2];
    for(i=grid0+1; i<grid; i++)
      if(hist[i][2]>max) max=hist[i][2];
      for(sum=0.0, i=grid0; i<grid; i++){ hist[i][2]-=max; sum+=exp(hist[i][2]); }
      for(i=grid0; i<grid; i++) hist[i][2]=hist[i][2]-log(sum)+log(100.0);
      /*
      ins=fopen("aa.est", "a");
      fprintf(ins, "delta=%g \n", delta);
      if(ins==NULL){ Rprintf("Can't write to file\n"); }
      for(i=grid0; i<grid; i++){
      fprintf(ins, "%5d  %10.6f  %10.6f  %10.6f  %g\n",i,hist[i][1],exp(hist[i][2]),
              hist[i][3],hist[i][2]);
      }
      fclose(ins);
      */ 
      
      ins=fopen("aa.fit", "w");
      if(ins==NULL){ Rprintf("Can't write to file\n");  }
      for(i=1; i<=(*data_num); i++){ 
        // fprintf(ins, " %d", i);
        for(j=1; j<=(*OUT_UNIT); j++) fprintf(ins, "%g",FIT[i][j]);
        fprintf(ins, "\n");
      }
      fclose(ins);
      
      ins=fopen("aa.pred","w");
      for(i=1; i<=(*test_num); i++){
        //  fprintf(ins, " %d", i);
        for(j=1; j<=(*OUT_UNIT); j++) fprintf(ins, "  %g",PRED[i][j]);
        fprintf(ins, "\n");
      }
      fclose(ins);
      
      
      PutRNGstate();
      free_dmatrix(x,1,maxpopN,1,dim);
      free_dmatrix(sx,1,maxpopN,1,dim);
      free_dvector(y,1,dim);
      free_imatrix(WW,1,maxpopN,1,dim);
      free_imatrix(sWW,1,maxpopN,1,dim);
      free_ivector(W,1,dim);
      free_dvector(fvalue,1,maxpopN);
      free_dvector(sfvalue,1,maxpopN);
      free_dvector(t,1,maxpopN);
      free_dmatrix(total_mat,1,((*data_num)+(*test_num)),1,(*P));
      free_dmatrix(total_y,1,((*data_num)+(*test_num)),1,(*OUT_UNIT));
      free_dmatrix(data_mat,1,(*data_num),1,(*P));
      free_dmatrix(test_mat,1,(*test_num),1,(*P));
      free_dmatrix(data_y,1,(*data_num),1,(*OUT_UNIT));
      free_dmatrix(test_y,1,(*test_num),1,(*OUT_UNIT));
      free_ivector(sel,1,(*data_num));
      free_ivector(rem,1,((*data_num)+(*test_num)));
      free_dmatrix(FIT,1,(*data_num),1,(*OUT_UNIT));
      free_dmatrix(PRED,1,(*test_num),1,(*OUT_UNIT));
      free_dvector(ps2,1,(*OUT_UNIT));
      free_dvector(wei,1,(*data_num)); 
      free_dmatrix(hist,0,grid,1,3);
      free_dvector(refden,0,grid);
      free_ivector(pos,0,grid);
      free_dvector(locfreq,0,grid);
      free_dvector(net,1,dim);
      free_ivector(indx,1,1000);
      free_dvector(fx,1,1000);
      free_dvector(fy,1,1000);
      free_dvector(prob,1,(*P)+1);
      free_dvector(mar,1,(*P)+1);
      free_dvector(AIC,1,maxpopN);
      free_dvector(BIC,1,maxpopN);
      free_dvector(EBIC,1,maxpopN);
      free_dvector(sAIC,1,(*popN)*(*total_iteration));
      free_dvector(sBIC,1,(*popN)*(*total_iteration));
      free_dvector(aAIC,1,(*popN)*(*total_iteration));
      free_dvector(aBIC,1,(*popN)*(*total_iteration));
      
      //return 0;
      }            
              
      
void
R_init_OpenMPController(DllInfo *info)
{
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}



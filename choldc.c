#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "allocate.h"
#define EPS 0.000000001

int choldc(double **L, double **a, int N);

/****************************************************************************
    Cholevsky decomposition of a symmetric and positive definite matrix
    A = LL^t, with L lower triangular
    Uses: L_ii = sqrt(A_ii- sum_{k<i} L_ik ^2)
    L_ji = (A_ij - sum_{k<i} L_ik*L_jk)/L_ii
    Input: A[N][N] (only upper diagonal needed)
    The matrix L is returned in the lower triangle of A

*****************************************************************************/
int choldc(double **L, double **a, int N)
{
  int i, j, k, ik, jk;
  double sum, diag[N];
  for(i=0; i<N; i++){
    for(j=i; j<N; j++){
      sum=a[i][j];
      for(k=i-1; k>=0; k--)sum-=L[i][k]*L[j][k];
      if(j==i){
	if(sum<EPS){
	  printf("WARNING, Choldc failed L_ii^2= %.2g i=%d\n", sum, i);
	  if(L==a){
	    // Ugo: Same matrix, restore lower diagonal
	    for(ik=0; ik<i; ik++){
	      a[ik][ik]=diag[ik];
	      for(jk=ik+1; jk<N; jk++)a[jk][ik]=a[ik][jk];
	    }
	  }
	  //exit(8);
	  return(-1);
	}
	diag[i]=a[i][i]; L[i][i]=sqrt(sum);
      }else{
	L[j][i]=sum/L[i][i]; // j<=i 
      }
    }
  }
  return(0);
}

void Forward_substitution(double *X, double **L, double *Y, int N){
  // Solve the matrix equation LX=Y, with L=lower diagonal (only i>=j)
  // X_i = (Y_i - sum_j<i L_ij X_j)/L_ii
  int i, j; double *XX=X;
  for(i=0; i<N; i++){
    double *XY=X, *LL=L[i];
    *XX=Y[i]; for(j=0; j<i; j++){(*XX)-=(*LL)*(*XY); XY++; LL++;}
    (*XX)/=(*LL);
    XX++;
  }
}

void Backward_substitution(double *X, double **L, double *Y, int N)
{
  // Solve the matrix equation L^tX=Y, with L=lower diagonal (only i>=j)
  // X_i =(Y_i - sum_j>i L_ji X_j)/L_ii
  int i, j;
  for(i=N-1; i>=0; i--){
    X[i]=Y[i];
    for(j=i+1; j<N; j++)X[i]-= X[j]*L[j][i];
    X[i]/=L[i][i];
  }
}

float **Cholevsky_inversion_d2f(double **L, int N){
  float **L_inv=Allocate_mat2_f(N, N);
  double X[N];
  int i, j, k;
  for(k=0; k<N; k++){
  // Solve the matrix equation L(i,j)X(j,k)=delta(i,k),
  // with L=lower diagonal (only i>=j)
  // X_jk = (delta_ik - sum_j<i L_ij X_jk)/L_ii
    for(i=0; i<N; i++){
      double *LL=L[i], *XY=X, sum=0; if(i==k)sum=1;
      for(j=0; j<i; j++){sum-=(*LL)*(*XY); XY++; LL++;}
      X[i]=sum/(*LL);
    }
    for(i=0; i<N; i++)L_inv[i][k]=X[i];
  }
  return(L_inv);
}

double **Cholevsky_inversion(double **L, int N){
  double **L_inv=Allocate_mat2_d(N, N), X[N];
  int i, j, k;
  for(k=0; k<N; k++){
  // Solve the matrix equation L(i,j)X(j,k)=delta(i,k),
  // with L=lower diagonal (only i>=j)
  // X_jk = (delta_ik - sum_j<i L_ij X_jk)/L_ii
    for(i=0; i<N; i++){
      double *LL=L[i], *XY=X, sum=0; if(i==k)sum=1;
      for(j=0; j<i; j++){sum-=(*LL)*(*XY); XY++; LL++;}
      X[i]=sum/(*LL); L_inv[i][k]=X[i];
    }
  }
  return(L_inv);
}

int choldc_f(float **L, float **a, int N)
{
  int i, j, k, ik, jk;
  double sum, diag[N];
  for(i=0; i<N; i++){
    for(j=i; j<N; j++){
      sum=a[i][j];
      for(k=i-1; k>=0; k--)sum-=L[i][k]*L[j][k];
      if(j==i){
	if(sum<EPS){
	  printf("WARNING, Choldc failed L_ii^2= %.2g i=%d\n", sum, i);
	  if(a==L){
	    // Ugo: Restore lower diagonal 
	    for(ik=0; ik<i; ik++){
	      a[ik][ik]=diag[ik];
	      for(jk=ik+1; jk<N; jk++)a[jk][ik]=a[ik][jk];
	    }
	  }
	  return(-1);
	}
	diag[i]=a[i][i]; L[i][i]=sqrt(sum);
      }else{
	L[j][i]=sum/L[i][i]; // j<=i 
      }
    }
  }
  return(0);
}

void Forward_substitution_f(double *X, float **L, double *Y, int N){
  // Solve the matrix equation LX=Y, with L=lower diagonal (only i>=j)
  // X_i = (Y_i - sum_j<i L_ij X_j)/L_ii
  int i, j; double *XX=X;
  for(i=0; i<N; i++){
    double *XY=X; float *LL=L[i];
    *XX=Y[i]; for(j=0; j<i; j++){(*XX)-=(*LL)*(*XY); XY++; LL++;}
    (*XX)/=(*LL);
    XX++;
  }
}


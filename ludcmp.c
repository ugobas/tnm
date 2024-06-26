#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define TINY 1.0e-20;

/*
Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above; indx[1..n] is an output vector that records the row permutation effected by the partial pivoting; d is output as ±1 depending on whether the number of row interchanges was even or odd, respectively. This routine is used in combination with lubksb to solve linear equations or invert a matrix.
 */

float ludcmp(double **a, int n, int *indx)
{
  float d=1.0;
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  double vv[n];
 
  for(i=0;i<n;i++) {
    big=0.0;
    for(j=0;j<n;j++){
      if ((temp=fabs(a[i][j])) > big) big=temp;
    }
    if(big == 0.0){
      printf("ERROR, singular matrix in routine LUDCMP\n"); exit(8);
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if(j != imax) {
      for (k=0; k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i< n;i++) a[i][j] *= dum;
    }
  }
  return(d);
}

#undef TINY

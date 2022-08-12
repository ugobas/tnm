#include "vector.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/************************ Codes  ****************************/

void Vector_product(float *w, float *v1, float *v2)
{
  // w = v1 X v2
  w[0]=v1[1]*v2[2]-v1[2]*v2[1];
  w[1]=v1[2]*v2[0]-v1[0]*v2[2];
  w[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

void Vector_product_d(double *w, double *v1, double *v2)
{
  // w = v1 X v2
  w[0]=v1[1]*v2[2]-v1[2]*v2[1];
  w[1]=v1[2]*v2[0]-v1[0]*v2[2];
  w[2]=v1[0]*v2[1]-v1[1]*v2[0];
}


float Scalar_product_3(float *x, float *y){
  return(x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
}

double Scalar_product_3_d(double *x, double *y){
  return(x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
}

void Rotate(float *r, float **rot){
  float tmp[3];
  tmp[0]=rot[0][0]*r[0]+rot[0][1]*r[1]+rot[0][2]*r[2];
  tmp[1]=rot[1][0]*r[0]+rot[1][1]*r[1]+rot[1][2]*r[2];
  tmp[2]=rot[2][0]*r[0]+rot[2][1]*r[1]+rot[2][2]*r[2];
  r[0]=tmp[0]; r[1]=tmp[1]; r[2]=tmp[2]; 
}

void Rotate_d(double *r, double **rot){
  double tmp[3];
  tmp[0]=rot[0][0]*r[0]+rot[0][1]*r[1]+rot[0][2]*r[2];
  tmp[1]=rot[1][0]*r[0]+rot[1][1]*r[1]+rot[1][2]*r[2];
  tmp[2]=rot[2][0]*r[0]+rot[2][1]*r[1]+rot[2][2]*r[2];
  r[0]=tmp[0]; r[1]=tmp[1]; r[2]=tmp[2]; 
}

void Subtract_vector_3(float *a, float *b, float *c){
  // a = b-c 
  a[0]=b[0]-c[0]; a[1]=b[1]-c[1]; a[2]=b[2]-c[2]; 
}

void Subtract_vector_3_d(double *a, double *b, double *c){
  // a = b-c 
  a[0]=b[0]-c[0]; a[1]=b[1]-c[1]; a[2]=b[2]-c[2]; 
}

void Sum_vector_3(float *a, float *b, float *c){
  // a + b-c 
  a[0]=b[0]+c[0]; a[1]=b[1]+c[1]; a[2]=b[2]+c[2]; 
}

void Sum_vector_3_d(double *a, double *b, double *c){
  // a + b-c 
  a[0]=b[0]+c[0]; a[1]=b[1]+c[1]; a[2]=b[2]+c[2]; 
}

void Matrix_multiplication(double *w, double **M, double *v, int n)
{
  int i, j;
  for(i=0; i<n; i++){
    double sum=0, *mi=M[i], *v0=v;
    for(j=0; j<n; j++){sum+=(*mi)*(*v0); mi++; v0++;}
    w[i]=sum;
  }
}

void Transform_matrix(float A[3][3], float **R){
  // Similarity transformation A'=RAR^T
  float tmp[3][3]; int i, j, k;
  for(i=0; i<3; i++){
    float *Ai=A[i];
    for(j=0; j<3; j++){
      float t=0, *Rj=R[j];
      for(k=0; k<3; k++)t+=Ai[k]*Rj[k];
      tmp[i][j]=t; 
    }
  }
  for(i=0; i<3; i++){
    float *Ri=R[i];
    for(j=0; j<3; j++){
      float t=0;
      for(k=0; k<3; k++)t+=Ri[k]*tmp[k][j];
      A[i][j]=t; 
    }
  }
}

void Transform_matrix_d(double **A, double **R){
  // Similarity transformation A'=RAR^T
  double tmp[3][3]; int i, j, k;
  for(i=0; i<3; i++){
    double *Ai=A[i];
    for(j=0; j<3; j++){
      double t=0, *Rj=R[j];
      for(k=0; k<3; k++)t+=Ai[k]*Rj[k];
      tmp[i][j]=t; 
    }
  }
  for(i=0; i<3; i++){
    double *Ri=R[i];
    for(j=0; j<3; j++){
      double t=0;
      for(k=0; k<3; k++)t+=Ri[k]*tmp[k][j];
      A[i][j]=t; 
    }
  }
}


float Scalar_product(float *xx, float *yy, int n){
  double q=0; float *x=xx, *y=yy; int i;
  for(i=0; i<n; i++){q+=(*x)*(*y); x++; y++;}
  return((float)q);
}

float Scalar_product_weighted(float *xx, float *yy, float *ww, int n){
  double q=0; float *x=xx, *y=yy, *w=ww; int i;
  for(i=0; i<n; i++){q+=(*x)*(*y)*(*w); x++; y++; w++;}
  return((float)q);
}

void Normalize_vector(float *vv, int n){
  int i; float *v;
  float q=Scalar_product(vv, vv, n); q=sqrt(q);
  v=vv; for(i=0; i<n; i++){*v/=q; v++;}
}

void Normalize_vector_weighted(float *vv, float *ww, int n){
  int i; float *v=vv, *w=ww; double q=0;
  for(i=0; i<n; i++){q+=(*v)*(*v)*(*w); v++; w++;} q=1./sqrt(q);
  v=vv; for(i=0; i<n; i++){(*v)*=q; v++;}
}

float Corr_coeff(float *xx, float *yy, int N,
		 float *slope, float *offset)
{
  double xy=0, x1=0, y1=0, x2=0, y2=0;
  float *x=xx, *y=yy; int i;

  for(i=0; i<N; i++){
    xy += (*x)*(*y);
    x1+=(*x); x2+= (*x)*(*x);
    y1+=(*y); y2+= (*y)*(*y);
    x++; y++;
  }
  x2=N*x2-x1*x1; y2=N*y2-y1*y1; xy=N*xy-x1*y1;
  if((x2==0.0)||(y2==0.0))return(0.0);
  if(slope)*slope=xy/x2;
  if(offset)*offset=(y1-(*slope)*x1)/N;
  //printf("N= %d xy=%.2g x2=%.2g y2=%.2g\n", N, xy, x2, y2);
  return(xy/sqrt(x2*y2));
}



float Corr_coeff_w(float *xx, float *yy, int *w, int N,
		   float *slope, float *offset, int *norm)
{
  double xy=0, x1=0, y1=0, x2=0, y2=0;
  float *x=xx, *y=yy; int i;

  (*norm)=0;
  for(i=0; i<N; i++){
    if(w[i]){
      xy += (*x)*(*y); (*norm)++;
      x1+=(*x); x2+= (*x)*(*x);
      y1+=(*y); y2+= (*y)*(*y);
    }
    x++; y++;
  }
  x2=(*norm)*x2-x1*x1; y2=(*norm)*y2-y1*y1; xy=(*norm)*xy-x1*y1;
  *slope=xy/x2; *offset=(y1-(*slope)*x1)/(*norm);
  //printf("N= %d xy=%.2g x2=%.2g y2=%.2g\n", N, xy, x2, y2);
  return(xy/sqrt(x2*y2));
}

float Collectivity_norm1(float *v, int N){
  double sum=0, Entropy=0;
  for(int i=0; i<N; i++){
    float w=v[i]; if(w<=0)continue; sum+=w; Entropy-=w*log(w);
  }
  if(sum==0){
    printf("ERROR in collectivity_norm1, sum= %.2g\n", sum);
    return(0);
  }
  Entropy=Entropy/sum+log(sum);
  float kappa=exp(Entropy);
  return(kappa);
}

float Collectivity_norm2(float *v, int N){
  float kappa, w; double sum=0, Entropy=0; int i;
  for(i=0; i<N; i++){
    w=v[i]*v[i]; if(w>0){sum+=w; Entropy-=w*log(w);}
  }
  if(sum>0){
    Entropy=Entropy/sum;
    Entropy+=log(sum);
  }
  kappa=exp(Entropy);
  if(isnan(kappa)){
    printf("Warning in Collectivity_norm2 ");
    printf("Entropy: %.2g sum: %.2g kappa: %.2g\n", Entropy, sum, kappa);
    printf("v: ");
    for(i=0; i<30; i++)printf(" %.2g", v[i]);
    printf(" ...\n");
    //exit(8);
  }
  return(kappa);
}

float Collectivity_norm2_outlier(float *v, int N, int *outlier){
  float kappa, w; double sum=0, Entropy=0; int i;
  for(i=0; i<N; i++){
    if(outlier && (outlier[i]))continue;
    w=v[i]*v[i]; if(w>0){sum+=w; Entropy-=w*log(w);}
  }
  Entropy=Entropy/sum+log(sum);
  kappa=exp(Entropy);
  return(kappa);
}


float Collectivity_Renyi_norm1(float *v, int N){
  double sum=0, sum2=0, SR; int i; float w;
  for(i=0; i<N; i++){
    w=v[i]; if(w<=0)continue; sum+=w; sum2+=w*w;
  }
  SR=(sum*sum)/sum2;
  return(SR);
}

float Collectivity_Renyi_norm2(float *v, int N){
  double sum=0, sum2=0, SR; int i; float w;
  for(i=0; i<N; i++){
    w=v[i]*v[i]; if(w<=0)continue; sum+=w; sum2+=w*w;
  }
  SR=(sum*sum)/sum2;
  return(SR);
}

float Area_norm1(float *v, int N){
  double sum=0, sum2=0, imean; int i; float w;
  for(i=0; i<N; i++){
    w=v[i]; if(w<=0)continue; sum+=w; sum2+=i*w;
  }
  imean=sum2/sum;
  return((N-imean-1)/(N-1));
}

float Area_norm2(float *v, int N){
  double sum=0, sum2=0, imean; int i; float w;
  for(i=0; i<N; i++){
    w=v[i]*v[i]; if(w<=0)continue; sum+=w; sum2+=i*w;
  }
  imean=sum2/sum;
  return((N-imean-1)/(N-1));
}


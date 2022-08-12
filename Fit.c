#include <math.h>
#include <stdio.h>

double Corr[2][2], corr_y[2], X_ave[2], X_dev[2];
double Y1=0, Y2=0; 

float Fit_2(float *Y, float *X0, float *X1, int N,
	    float *A0, float *A1, float *bias, int first)
{
  float mu1=0;
  // 0: do not allow offset 1: do not penalize offset 0.5: penalize

  // Fit Y=XA with X=(2*n) setting bias b=0
  // and evaluating it as b=<Y>-<X>A
  // A=<X^T X>^{-1}<X^T Y>
  float *x[2]; x[0]=X0; x[1]=X1;
  int amax=0, a, b, i;
  if(1){ // Compute also for variable X1, which is the same at all calls
    amax=1;
    Y1=0; Y2=0; for(i=0; i<N; i++){Y1+=Y[i]; Y2+=Y[i]*Y[i];}
    Y1/=N; Y2=(Y2/N-Y1*Y1); 
  }


  for(a=amax; a>=0; a--){
    float *xa=x[a];
    double s=0; for(i=0; i<N; i++){s+=xa[i];} X_ave[a]=s/N;
    //s=0; for(i=0; i<N; i++){s+=xa[i]*xa[i];}
    //X_dev[a]=sqrt(s/N-X_ave[a]*X_ave[a]);
    s=0; for(i=0; i<N; i++){s+=xa[i]*Y[i];}
    corr_y[a]=s/N-X_ave[a]*Y1;
    for(b=1; b>=a; b--){
      float *xb=x[b];
      s=0; for(i=0; i<N; i++){s+=xa[i]*xb[i];}
      Corr[a][b]=s/N-mu1*X_ave[a]*X_ave[b];
      if(b>a)Corr[b][a]=Corr[a][b];
    }
  }

  double Corr_inv[2][2];
  double det=Corr[0][0]*Corr[1][1]-Corr[0][1]*Corr[1][0];
  Corr_inv[0][0]=Corr[1][1]/det;
  Corr_inv[1][1]=Corr[0][0]/det;
  Corr_inv[0][1]=-Corr[1][0]/det;
  Corr_inv[1][0]=-Corr[0][1]/det;

  double A[2], r2=0; *bias=Y1;
  for(a=0; a<2; a++){
    A[a]=0; for(b=0; b<2; b++)A[a]+=Corr_inv[a][b]*corr_y[b];
    r2+=A[a]*corr_y[a];
    (*bias)-=X_ave[a]*A[a];
  }
  *A0=A[0];
  *A1=A[1];

  return(sqrt(r2/Y2));
}

float Fit_1(float *Y, float *X0, int N, float *A0, float *bias)
{
  float mu1=0;
  // 0: do not allow offset 1: do not penalize offset 0.5: penalize

  // Fit Y=XA with X=(1*n) setting bias b=0
  // and evaluating it as b=<Y>-<X>A
  // A=<X^T X>^{-1}<X^T Y>
  float *x[1]; x[0]=X0; int i;
  if(1){ // Compute also for variable X1, which is the same at all calls
    Y1=0; Y2=0; for(i=0; i<N; i++){Y1+=Y[i]; Y2+=Y[i]*Y[i];}
    Y1/=N; Y2=sqrt(Y2/N-Y1*Y1); 
  }

  int a=0;
  float *xa=x[a];
  double s=0; for(i=0; i<N; i++){s+=xa[i];} X_ave[a]=s/N;
  //s=0; for(i=0; i<N; i++){s+=xa[i]*xa[i];}
  //X_dev[a]=sqrt(s/N-X_ave[a]*X_ave[a]);
  s=0; for(i=0; i<N; i++){s+=xa[i]*Y[i];}
  corr_y[a]=s/N-X_ave[a]*Y1;
  s=0; for(i=0; i<N; i++){s+=xa[i]*xa[i];}
  Corr[a][a]=s/N-mu1*X_ave[a]*X_ave[a];

  double A=corr_y[a]/Corr[0][0];
  double r2=A*corr_y[a];
  (*bias)=Y1-X_ave[a]*A;
  *A0=A;

  return(sqrt(r2/Y2));
}

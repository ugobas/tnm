#include "allocate.h"
#include "diagonalize.h"
#include "EC.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define VERBOSE 1
#define EPS 0.00001
#define EPS2 0.0000000001
#define IT_MAX 100
#define W_MIN 0.000001 // Weight below which eigenvectors are discarded
#define X_MIN -100
#define X_MAX 10000


// Global variables
float *EC, *eval, **evec, *w_ev;
float Lambda, c1c2;
int N;


float C1sqr_over_C2(double **Corr, int N, int l);
static int Find_zero_Newton(float *x0, float X_min, float X_max);
static int Find_zero_quad(float *x0);
static float Compute_score(float Lambda);
static float Compute_deriv(float Lambda, float *f);
static void Compute_EC(float *EC, int N, float *eval, float **evec,
		       float *w_ev, float Lambda, float AA);
static float Score_c1c2(float Lambda, float c1c2, int N,
			float *eval, float **evec, float *w_ev);

static float Derivative_c1c2(float *f0, float Lambda, float c1c2, int N,
			     float *eval, float **evec, float *w_ev);

float *EC_profile(float *PE, double **Mat, int N_in, char *name)
{
  N=N_in;
  float *EC_out=malloc(N*sizeof(float));
  EC=EC_out;
  int i,j;

  // Eigensystem
  eval=malloc(N*sizeof(float));
  evec=malloc(N*sizeof(float *));
  for(i=0; i<N; i++)evec[i]=malloc(N*sizeof(float));
  d_Diagonalize(N, Mat, eval, evec, 1);

  w_ev=malloc(N*sizeof(float));
  for(i=0; i<N; i++){
    float *v=evec[i];
    double c1=0; for(j=0; j<N; j++)c1+=v[j];
    c1/=N; for(j=0; j<N; j++)v[j]/=c1;
    w_ev[i]=N*c1*c1;
  }

  // c1c2 set equal to value for matrix
  float *ci=CV_profile(Mat, N);
  double cv1=0, cv2=0; float *c=ci;
  for(i=0; i<N; i++){cv1+=*c; cv2+=(*c)*(*c); c++;}
  free(ci);
  if(1){
    c1c2=(cv1*cv1)/(N*cv2);
  }else{
    // Choose variance such that the expected minimum component is zero
    double eps=1./(2*log((float)N)-log(6.283));
    c1c2=1./(1+eps);
  }
  printf("<EC_i>^2/<EC_i^2>= %.3f\n", c1c2);

  // Compute lambda;
  Lambda=eval[0]*1.05;
  if(Find_zero_Newton(&Lambda, X_MIN, X_MAX))Find_zero_quad(&Lambda);

  //float E=Compute_score(Lambda);
  Compute_EC(EC, N, eval, evec, w_ev, Lambda, 1);
  if(name){printf("EC of %s", name);}else{printf("EC");}
  printf(" computed, lambda1= %.3g Lambda= %.3g\n",eval[0], Lambda);

  // Normalize so that mean is mean of matrix
  double c1=0;
  for(i=0; i<N; i++)c1+=EC[i];
  for(i=0; i<N; i++)EC[i]*=(cv1/c1);

  if(PE!=NULL){
    float *v=evec[0]; for(i=0; i<N; i++)PE[i]=v[i];
  }
  Empty_matrix_f(evec, N); free(eval); free(w_ev);
  return(EC_out);
}

float Compute_score(float Lambda){
  return(Score_c1c2(Lambda, c1c2, N, eval, evec, w_ev));
}

float Compute_deriv(float Lambda, float *f){
  return(Derivative_c1c2(f, Lambda, c1c2, N, eval, evec, w_ev));
}

float C1sqr_over_C2(double **Corr, int N, int l)
{
  int i, j;
  double *ci=malloc(N*sizeof(double));
  for(i=0; i<N; i++){
    ci[i]=0; double *cc=Corr[i];
    for(j=0; j<=(i-l); j++){
      ci[i]+=cc[j];
      ci[j]+=cc[j];
    }
  }
  double c1=0, c2=0;
  for(i=0; i<N; i++){c1+=ci[i]; c2+=ci[i]*ci[i];}
  free(ci);
  return((c1*c1)/(N*c2));
}

int Find_zero_Newton(float *x0, float X_min, float X_max){
  float f1, f, d; int iter;
  float f_opt=1, x_opt=*x0;
  float x_min=*x0, x_max=*x0, f_min=0, f_max=1;
  int i_min=0, i_max=0;
  float x_ini=*x0;

  for(iter=0; iter<IT_MAX; iter++){
    f1=Compute_deriv((*x0), &f);
    if((i_min==0)&&(f<0)){i_min=1; f_min=f; x_min=*x0;}
    if((i_max==0)&&(f>0)){i_max=1; f_max=f; x_max=*x0;}
    if(VERBOSE)printf("%3d %.3g %.2g %.2g\n", iter, *x0, f, f1);
    if(fabs(f)<EPS)return(0);
    if(fabs(f)<fabs(f_opt)){f_opt=f; x_opt=*x0;}
    d = f/f1; //if(d>*x0)d=(*x0)*0.5; // Watch out!
    (*x0) -= d;
    if((*x0<X_min)||(*x0 > X_max))break;
  }
  *x0=x_opt;
  if((i_min)&&(i_max)){
    float x = x_min-f_min*(x_max-x_min)/(f_max-f_min);
    float f=Compute_score(x);
    if(fabs(f) < fabs(f_opt)){*x0=x; f_opt=f;}
  }
  if(fabs(f_opt)<EPS)return(0);
  if(VERBOSE)printf("%3d %.3g %.2g Not converged\n", iter, *x0, f_opt);
  *x0=x_ini;
  return(1);
}


int Find_zero_quad(float *x0)
{
  int iter=0;
  float x1=(*x0)*(0.98), x2=*x0*(1.02), x1_new, x2_new;
  float f_opt=1, x_opt=x1;
  float x_min=x1, x_max=x1, f_min=0, f_max=1;
  int i_min=0, i_max=0;
  double f1, f2, b;

  // Initialization

  f1=Compute_score(x1); f1*=f1;
  f2=Compute_score(x2); f2*=f2;
  for(iter=0; iter<IT_MAX; iter++){
    if(VERBOSE)printf("%3d %.3g %.3g %.2g %.2g\n", iter,x1,x2,f1,f2);
    if(f1<f_opt){f_opt=f1; x_opt=x1;}
    if(f2<f_opt){f_opt=f2; x_opt=x2;}
    if((i_min==0)&&(f1<0)){i_min=1; f_min=f1; x_min=x1;}
    if((i_min==0)&&(f2<0)){i_min=1; f_min=f2; x_min=x2;}
    if((i_max==0)&&(f1>0)){i_max=1; f_max=f1; x_max=x1;}
    if((i_max==0)&&(f2>0)){i_max=1; f_max=f2; x_max=x2;}
    b=f1/f2;
    if(b==1)break;
    x1_new=(x1-b*x2)/(1-b);
    x2_new=(x1+b*x2)/(1+b);
    f1=Compute_score(x1_new); f1*=f1;
    f2=Compute_score(x2_new); f2*=f2;
    if(f1 <EPS2){*x0=x1_new; return(0);}
    if((x1_new<X_MIN)||(x2_new<X_MIN))break;
    if((x1_new>X_MAX)||(x2_new>X_MAX))break;
    x1=x1_new; x2=x2_new;
  }
  // Not converged
  *x0=x_opt;
  if((i_min)&&(i_max)){
    float x = x_min-f_min*(x_max-x_min)/(f_max-f_min);
    float f=Compute_score(x), f2=f*f;
    if(f_opt > f2){*x0=x; f_opt=f2;}
  }
  if(VERBOSE)printf("%3d %.3f %.8f\n", iter, *x0, f_opt);
  if(f_opt <EPS2)return(0);
  return(1);
}

void Compute_EC(float *EC, int N, float *eval, float **evec, float *w_ev,
		float Lambda, float AA)
{
  int i, j; float eps=0.0000001;
  for(i=0; i<N; i++)EC[i]=0;
  for(j=0; j< N; j++){
    if(w_ev[j] < W_MIN)continue;
    float diff=Lambda-eval[j];
    if(fabs(diff)<eps){if(diff>0){diff=eps;}else{diff=-eps;}}
    float w=w_ev[j]/diff, *v=evec[j], *ec=EC;
    for(i=0; i<N; i++){*ec+=w*(*v); ec++; v++;}
  }
}

float Score_c1c2(float Lambda, float c1c2, int N,
		 float *eval, float **evec, float *w_ev)
{
  double S1=0, S2=0, eps=0.000001;
  int j;
  for(j=0; j< N; j++){
    if(w_ev[j] < W_MIN)continue;
    float diff=Lambda-eval[j];
    if(fabs(diff)<eps){if(diff>0){diff=eps;}else{diff=-eps;}}
    double w=w_ev[j]/diff; S1+=w;
    w/=diff; S2+=w;
  }
  return(c1c2 - (S1*S1)/S2);
}

float Derivative_c1c2(float *f0, float Lambda, float c1c2, int N,
		      float *eval, float **evec, float *w_ev)
{
  // Derivative of -N<c>^2/<c^2> with respect to Lambda 
  double S1=0, S2=0, S3=0, eps=0.000001;
  int j;
  for(j=0; j< N; j++){
    if(w_ev[j] < W_MIN)continue;
    float diff=Lambda-eval[j];
    if(fabs(diff)<eps){if(diff>0){diff=eps;}else{diff=-eps;}}
    double w=w_ev[j]/diff; S1+=w;
    w/=diff; S2+=w;
    w/=diff; S3+=w;
  }
  (*f0)=c1c2 - (S1*S1)/S2;
  double f=2*S1*(1-(S1*S3)/(S2*S2));
  return(f);
}

float *CV_profile(double **Corr, int N)
{
  float *CV=malloc(N*sizeof(float));
  for(int i=0; i<N; i++){
    double sum=0, *c=Corr[i];
    for(int j=0; j<N; j++)sum+=c[j];
    CV[i]=sum;
  }
  return(CV);
}

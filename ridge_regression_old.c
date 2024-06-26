#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "diagonalize.h"
#include "ridge_regression.h"
#define VBS 1

static void Score_computation(float *E, float *Cv, float *RR, float *GCV,
			      float *A2, float *Cv2, float *Sv1,
			      float Lambda, float *lambda, float *ya2,
			      float E0, float Y2, int P, int N);
static float Get_A(float Lambda, float *yalpha, float **evec, float *eval,
		   float *D, int Npar, int a0);
static float Compute_error(float Lambda, float *lambda, float *ya2,
			   float Y2, int P);
static float Make_fit(float *A, float Lambda, float *lambda, float **evec,
		      float *yalpha, float *D, int Npar, float **X, float *Y,
		      int N);
static float Print_make_fit(struct ridge_fit *F, float Lambda, char *what,
			    FILE *file_out, float *eval, float **evec,
			    float *yalpha, float *D, float **X, float *Y,
			    float *ya2, float E0, float Y2, int Npar, int N);
extern float Corr_coeff(float *xx, float *yy, int N,
			float *slope, float *offset);
void Print_fit(struct ridge_fit *F, double *sum_x, double sum_y,
	       int Npar, int N, char *nameout, char *what,
	       float **X, int VERBOSE);

float Find_Lambda(float *Ec, float *Cvc, float *RRc, float *GCVc,
		  float *A2c, float *Cv2c, float *Sv1c,
		  float L1, float L2, float C1, float C2,
		  float *eval, float *ya2, float E0, float Y2,
		  int Npar, int N);

float Ridge_regression(struct ridge_fit *Sv1_fit,
		      char *nameout, float **X, float *Y, int N, int Npar)
{
  int i, a, b;

  // Compute the covar
  printf("Computing and diagonalizing covariance matrix\n");
  float **C=malloc(Npar*sizeof(float *));
  for(a=0; a<Npar; a++){
    C[a]=malloc(Npar*sizeof(float));
    for(b=0; b<Npar; b++)C[a][b]=0;
  }
  for(i=0; i<N; i++){
    float *x=X[i];
    for(a=0; a<Npar; a++){
      for(b=0; b<=a; b++)C[a][b]+=x[a]*x[b];
    }
  }
  float *D=malloc(Npar*sizeof(float));
  for(a=0; a<Npar; a++){
    D[a]=sqrt(C[a][a]);
    for(b=0; b<=a; b++){
      C[a][b]/=(D[a]*D[b]);
      C[b][a]=C[a][b];
    }
  }
  
  float *eval=malloc(Npar*sizeof(float));
  float **evec=malloc(Npar*sizeof(float *));
  for(a=0; a<Npar; a++)evec[a]=malloc(Npar*sizeof(float));
  f_Diagonalize(Npar, C, eval, evec, 1, 0.000001, NULL, 0);
  printf("lambda= ");
  for(a=0; a<Npar; a++)printf(" %.4g", eval[a]);
  printf("\n");

  // Sequence of values for lambda
  int K=7, Kp=50, Nl=K*(Npar+Kp+1)+1, j=0, k;
  float *lambda=malloc(Nl*sizeof(float));
  float lambda_min=-2*eval[0], delta, la;
  for(a=0; a<Npar; a++){
    delta=(-eval[a]-lambda_min)/(K+1);
    la=lambda_min; lambda_min=-eval[a];
    for(k=0; k<K; k++){
      la+=delta; lambda[j]=la; j++;
    }
  }
  delta=-lambda_min/(K+1);
  la=lambda_min; lambda_min=0;
  for(k=0; k<K; k++){
    la+=delta; lambda[j]=la; j++;
  }
  lambda[j]=0; j++;
  delta=1.2*eval[0]/(Kp*K+1);
  la=lambda_min;
  for(k=0; k<Kp*K; k++){
    la+=delta; lambda[j]=la; j++;
  }


  // Projection of Y over principal components
  double sum;
  float *ya=malloc(Npar*sizeof(float));
  for(a=0; a<Npar; a++){
    sum=0; for(i=0; i<N; i++)sum+=X[i][a]*Y[i];
    ya[a]=sum;
  }
  float *yalpha=malloc(Npar*sizeof(float));
  float *ya2=malloc(Npar*sizeof(float));
  for(a=0; a<Npar; a++){
    sum=0; for(b=0; b<Npar; b++)sum+=evec[a][b]*ya[b]/D[b];
    yalpha[a]=sum;
    ya2[a]=sum*sum;
  }


  // Compute
  int ini_Cv2=0, ini_Sv=0;
  float Y2=0; for(i=0; i<N; i++)Y2+=Y[i]*Y[i];
  float E0=Compute_error(0, eval, ya2, Y2, Npar);
  float E, Cv, RR, GCV, A2, A0, Cv2, Sv1, LS;
  float Ec, Cvc, RRc, GCVc, A2c, Cv2c, Sv1c, Lc;
  float Sv1_old=1000, C2_old=1000; 
  float min_RR=1000, min_GCV=0, max_Cv=0, max_Cv2=0, max_Sv=0;
  float L_RR=-1000, L_GCV=-1000, L_Cv=-1000, L_Cv2=-1000,L_Sv=-1000; 

  char namefile[200];
  sprintf(namefile, "%s.dat", nameout);
  FILE *file_out=fopen(namefile, "w");
  printf("Writing %s\n", namefile); 
  if(VBS){
    fprintf(file_out,
	    "# Lambda RangeRisk GenCrossVal Error dE/dL d^3E/dL^3 Sv1 A0\n");
    fprintf(file_out,
	    "# N= %d Npar= %d\n", N, Npar);
  }
  // Loop of Lambda
  for(j=0; j<Nl; j++){
    float Lambda=lambda[j];
    Score_computation(&E, &Cv, &RR, &GCV, &A2, &Cv2, &Sv1,
		      Lambda, eval, ya2, E0, Y2, Npar, N);

    if((j==0)||(RR < min_RR)){min_RR=RR; L_RR=Lambda;}
    if((j==0)||(GCV< min_GCV)){min_GCV=GCV; L_GCV=Lambda;}
    if((Lambda>=0)&&(Cv>max_Cv)){max_Cv=Cv; L_Cv=Lambda;}
    // Cv2
    if((C2_old <0)&&(Cv2 >0)){ // Zero of Cv2=d^2C/DL^2

      Lc=Find_Lambda(&Ec, &Cvc, &RRc, &GCVc, &A2c, &Cv2c, &Sv1c,
		     lambda[j-1], Lambda, C2_old, Cv2,
		     eval, ya2, E0, Y2, Npar, N);
      LS=Lc*A2c;
      if((ini_Cv2==0)||(LS > max_Cv2)){max_Cv2=LS; L_Cv2=Lc;}
      if(RRc < min_RR){min_RR=RRc; L_RR=Lc;}
      if(GCVc < min_GCV){min_GCV=GCVc; L_GCV=Lc;}
      if((Lc>=0)&&(Cvc > max_Cv)){max_Cv=Cvc; L_Cv=Lc;}
      if(VBS){
	A0=Get_A(Lc, yalpha, evec, eval, D, Npar, 0);
	fprintf(file_out,
		"%5.3f\t%7.4g\t%7.4g\t%7.4g\t%9.3g\t%9.3g\t%9.3g\t%7.4g\n",
		Lc, RRc, GCVc, Ec/N, Cvc/N, Cv2c/N, Sv1c/N, A0); //A2/N, 
      }
      ini_Cv2=1;
    }
    C2_old=Cv2;

    // Sv1  // Zero of Sv1=d(Lambda S)/dLambda
    if((Sv1_old >0)&&(Sv1 <0)){
      Lc=(lambda[j-1]*Sv1-Lambda*Sv1_old)/(Sv1-Sv1_old);
      Score_computation(&Ec, &Cvc, &RRc, &GCVc, &A2c, &Cv2c, &Sv1c,
			Lc, eval, ya2, E0, Y2, Npar, N);
      LS=Lc*A2c;
      if((ini_Sv==0)||(LS > max_Sv)){max_Sv=LS; L_Sv=Lc;}
      if(RRc < min_RR){min_RR=RRc; L_RR=Lc;}
      if(GCVc < min_GCV){min_GCV=GCVc; L_GCV=Lc;}
      if((Lc>=0)&&(Cvc > max_Cv)){max_Cv=Cvc; L_Cv=Lc;}
      if(VBS){
	A0=Get_A(Lc, yalpha, evec, eval, D, Npar, 0);
	fprintf(file_out,
		"%5.3f\t%7.4g\t%7.4g\t%7.4g\t%9.3g\t%9.3g\t%9.3g\t%7.4g\n",
		Lc, RRc, GCVc, Ec/N, Cvc/N, Cv2c/N, Sv1c/N, A0); //A2/N, 
      }
      ini_Sv=1;
    }
    Sv1_old=Sv1;

    if(VBS){
      A0=Get_A(Lambda, yalpha, evec, eval, D, Npar, 0);
      fprintf(file_out,
	      "%5.3f\t%7.4g\t%7.4g\t%7.4g\t%9.3g\t%9.3g\t%9.3g\t%7.4g\n",
	      Lambda, RR, GCV, E/N, Cv/N, Cv2/N, Sv1/N, A0); //A2/N, 
    }
  }
  LS=lambda[Nl-1]*A2;
  if((Cv2<0)&&(LS > max_Cv2)){max_Cv2=LS; L_Cv2=lambda[Nl-1];} 
  if((Sv1>0)&&(LS > max_Sv )){max_Sv=LS; L_Sv=lambda[Nl-1];} 


  // Print results
  float r;
  struct ridge_fit OLS_fit; OLS_fit.A=malloc(Npar*sizeof(float));
  r=Print_make_fit(&OLS_fit, 0.00, "Ordinary Least Squares (OLS)",
		   file_out,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);

  struct ridge_fit RR_fit;  RR_fit.A=malloc(Npar*sizeof(float));
  r=Print_make_fit(&RR_fit, L_RR, "Range Risk (RR)",
		   file_out,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);

  struct ridge_fit GCV_fit; GCV_fit.A=malloc(Npar*sizeof(float));
  r=Print_make_fit(&GCV_fit, L_GCV, "Generalized Cross Validation (GCV)",
		   file_out,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);

  struct ridge_fit Cv_fit;  Cv_fit.A=malloc(Npar*sizeof(float));
  r=Print_make_fit(&Cv_fit, L_Cv, "dE/dLambda (Cv)",
		   file_out,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);

  struct ridge_fit Cv2_fit; Cv2_fit.A=malloc(Npar*sizeof(float));
  r=Print_make_fit(&Cv2_fit, L_Cv2, "d^2 Cv/dLambda^2 (Cv2)",
		   file_out,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);

  r=Print_make_fit(Sv1_fit, L_Sv, "d(Lambda S)/dLambda (entropy)",
		   file_out,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);

  // Compute fractions
  double *sum_x=malloc(Npar*sizeof(double)), sum_y=0, sum_y2=0;
  for(a=0; a<Npar; a++)sum_x[a]=0;
  for(i=0; i<N; i++){
    sum_y+=Y[i]; sum_y2+=Y[i]*Y[i];
    for(a=0; a<Npar; a++)sum_x[a]+=X[i][a];
  }
  sum_y2=(sum_y2-sum_y*sum_y/N);

  // Print results of the fits
  fprintf(file_out, "# Normalization of the error: %.2f\n", sum_y);
  fprintf(file_out, "# Force_const Lambda Rel.error ");
  fprintf(file_out, "Internal_frac Trans_frac Rot_frac fit\n");

  // Only B_ENM and translation:
  float *B_ENM=malloc(N*sizeof(float));
  for(i=0; i<N; i++)B_ENM[i]=X[i][0];
  float a0, a1, r_B=Corr_coeff(B_ENM, Y, N, &a0, &a1);
  fprintf(file_out, "#\t%.1f\t%.2f\t%.4f\t%.3f\t%.3f\t%.3f\t%s\n",
	  1./a0, 0.0, 1-r_B*r_B, a0*sum_x[0]/sum_y, a1*N/sum_y, 0.0,
	  "Norot_fit");
  fclose(file_out);

  //  Ordinary least squares
  Print_fit(&OLS_fit, sum_x, sum_y2, Npar, N, nameout, "OLS_fit", X, VBS);

  // Range risk
  Print_fit(&RR_fit, sum_x, sum_y2, Npar, N, nameout, "RR_fit", X, VBS);

  // Generalized Cross Validation
  Print_fit(&GCV_fit, sum_x, sum_y2, Npar, N, nameout, "GCV_fit", X, VBS);

  // dE/dLambda
  Print_fit(&Cv_fit, sum_x, sum_y2, Npar, N, nameout, "Cv_fit", X, VBS);

  // d^2Cv/dLambda^2
  Print_fit(&Cv2_fit, sum_x, sum_y2, Npar, N, nameout, "Cv2_fit", X, VBS);

  // d(Lambda S)/dLambda=0
  Print_fit(Sv1_fit, sum_x, sum_y2, Npar, N, nameout, "entropy_fit", X, VBS);

  float fact, Lambda; char name[80];
  fact=1.5; Lambda=fact*L_Sv; sprintf(name, "%.1fentro_fit", fact);

  struct ridge_fit Sv2; Sv2.A=malloc(Npar*sizeof(float));
  Print_make_fit(&Sv2, Lambda, "fact*(Lambda S)/dLambda (entropy)",
		 NULL,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);
  Print_fit(&Sv2, sum_x, sum_y2, Npar, N, nameout, name, X, VBS);
  
  fact=2.0; Lambda=fact*L_Sv; sprintf(name, "%.1fentro_fit", fact);
  Print_make_fit(&Sv2, Lambda, "fact*(Lambda S)/dLambda (entropy)",
		 NULL,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);
  Print_fit(&Sv2, sum_x, sum_y2, Npar, N, nameout, name, X, VBS);

  fact=10; Lambda=fact*L_Sv; sprintf(name, "%.1fentro_fit", fact);
  Print_make_fit(&Sv2, Lambda, "fact*(Lambda S)/dLambda (entropy)",
		 NULL,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);
  Print_fit(&Sv2, sum_x, sum_y2, Npar, N, nameout, name, X, VBS);

  fact=100; Lambda=fact*L_Sv; sprintf(name, "%.1fentro_fit", fact);
  Print_make_fit(&Sv2, Lambda, "fact*(Lambda S)/dLambda (entropy)",
		 NULL,eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N);
  Print_fit(&Sv2, sum_x, sum_y2, Npar, N, nameout, name, X, VBS);


  free(OLS_fit.A); free(RR_fit.A); free(GCV_fit.A); free(Cv_fit.A);
  free(Cv2_fit.A); free(Sv2.A);
  free(sum_x);


  free(lambda); free(eval);
  for(a=0; a<Npar; a++)free(evec[a]); free(evec);
  for(a=0; a<Npar; a++)free(C[a]); free(C); free(D);
  free(ya); free(yalpha); free(ya2);
  return(r);
}

float Compute_error(float Lambda, float *lambda, float *ya2, float Y2, int P)
{
  double E=0; int a;
  for(a=0; a<P; a++){
    float laLa=lambda[a]+Lambda;
    E+=ya2[a]*(laLa+Lambda)/(laLa*laLa);
  }
  return(Y2-E);
}

void Score_computation(float *E, float *Cv, float *RR, float *GCV,
		       float *A2, float *Cv2, float *Sv1,
		       float Lambda, float *lambda, float *ya2,
		       float E0, float Y2, int P, int N)
{
  int a;
  *E=Compute_error(Lambda, lambda, ya2, Y2, P);
  double yL2=0, yL3=0, layL3=0, yL5=0, L1=0; //LN=0, 
  for(a=0; a<P; a++){
    float laLa=lambda[a]+Lambda;
    float ylaLak=ya2[a]/(laLa*laLa);
    L1+=1./laLa;
    //LN+=1./(lambda[a]+N*Lambda);
    yL2+=ylaLak;
    ylaLak/=laLa;
    yL3+=ylaLak;
    layL3+=ylaLak*lambda[a];
    yL5+=ylaLak*(Lambda-lambda[a])/(laLa*laLa);
  }
  *Cv=Lambda*yL3;
  *Sv1=layL3-*Cv;
  *A2=yL2;
  *Cv2=yL5;
  float Trace=(N-P+Lambda*L1);
  //float Trace=(N-P+N*Lambda*LN);
  *GCV=N*(*E)/(Trace*Trace);
  *RR=(N-P)*(*E)/E0-2*Trace; *RR/=N;
}

float Make_fit(float *A, float Lambda, float *lambda, float **evec,
	       float *yalpha, float *D, int Npar, float **X, float *Y, int N)
{
  // Compute A
  int a, b, i; double sum;
  float *coeff=malloc(Npar*sizeof(float));
  for(a=0; a<Npar; a++)coeff[a]=yalpha[a]/(lambda[a]+Lambda);
  for(a=0; a<Npar; a++){
    sum=0; for(b=0; b<Npar; b++)sum+=coeff[b]*evec[b][a];
    A[a]=sum/D[a];
  }
  free(coeff);

  double xy=0, x2=0, y2=0, x1=0, y1=0;
  for(i=0; i<N; i++){
    float *Xi=X[i]; double Y_pred=0;
    for(a=0; a<Npar; a++)Y_pred+=Xi[a]*A[a];
    x1+=Y[i]; x2+=Y[i]*Y[i]; xy+=Y[i]*Y_pred;
    y1+=Y_pred; y2+=Y_pred*Y_pred;
  }
  x2=x2-x1*x1/N;
  y2=y2-y1*y1/N;
  xy=(xy-x1*y1/N)/sqrt(x2*y2);
  return(xy);
}

float Get_A(float Lambda, float *yalpha, float **evec, float *eval, float *D,
	    int Npar, int a0)
{
  double A0=0; int a;
  for(a=0; a<Npar; a++)A0+=yalpha[a]*evec[a][a0]/(eval[a]+Lambda);
  return(A0/D[a0]);
}

float Find_Lambda(float *Ec, float *Cvc, float *RRc, float *GCVc,
		  float *A2c, float *Cv2c, float *Sv1c,
		  float L1, float L2, float C1, float C2,
		  float *eval, float *ya2, float E0, float Y2,
		  int Npar, int N)
{
  float Lc, La1=L1, La2=L2, Cv1=C1, Cv2=C2; int i;
  for(i=0; i<30; i++){
    Lc=(La1*Cv2-La2*Cv1)/(Cv2-Cv1);
    Score_computation(Ec, Cvc, RRc, GCVc, A2c, Cv2c, Sv1c, Lc,
		      eval, ya2, E0, Y2, Npar, N);
    if(fabs(*Cv2c)<0.001)break;
    if(*Cv2c < 0){Cv1=*Cv2c; La1=Lc;}else{Cv2=*Cv2c; La2=Lc;}
  }
  return(Lc);
}

float Print_make_fit(struct ridge_fit *F, float Lambda, char *what,
		    FILE *file_out, float *eval, float **evec, float *yalpha,
		    float *D, float **X, float *Y, float *ya2, float E0,
		    float Y2, int Npar, int N)
{
  float Cv, RR, GCV, A2, Cv2, Sv1;
  F->Lambda=Lambda;
  Score_computation(&(F->E), &Cv, &RR, &GCV, &A2, &Cv2, &Sv1,
		    F->Lambda, eval, ya2, E0, Y2, Npar, N);
  float r=Make_fit(F->A, F->Lambda, eval, evec, yalpha, D, Npar, X, Y, N);
  if(file_out){
    fprintf(file_out, "# Optimal fit %s:\n# ", what);
    fprintf(file_out,
	    "%5.3f\t%7.4g\t%7.4g\t%7.4g\t%9.3g\t%9.3g\t%9.3g\t%7.4g\n",
	    F->Lambda, RR, GCV, F->E/N, Cv/N, Cv2/N, Sv1/N, F->A[0]);
  }
  return(r);
}

void Print_fit(struct ridge_fit *F, double *sum_x, double var_y,
	       int Npar, int N, char *nameout, char *what,
	       float **X, int VERBOSE)
{
  char namefile[200]; sprintf(namefile, "%s.dat", nameout);
  FILE *file_out=fopen(namefile, "a");
  int a; double sum=0; for(a=0; a<Npar; a++)sum+=F->A[a]*sum_x[a];
  fprintf(file_out, "#\t%.1f\t%.2f\t%.4f\t%.3f\t%.3f\t%.3f\t%s\n",
	  1./F->A[0], F->Lambda, F->E/var_y, F->A[0]*sum_x[0]/sum,
	  F->A[1]*N/sum, 1-(F->A[0]*sum_x[0]+F->A[1]*N)/sum, what);
  fclose(file_out);
  if(VERBOSE){
    sprintf(namefile, "%s_%s.dat", nameout, what);
    file_out=fopen(namefile, "w"); int i, a;
    printf("Writing %s fit in %s\n", what, namefile);
    for(i=0; i<N; i++){
      double Y=0;
      for(a=0; a<Npar; a++)Y+=F->A[a]*X[i][a];
      fprintf(file_out, "%.3f\n", Y);
    }
    fclose(file_out);
  }
}

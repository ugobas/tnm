#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "diagonalize.h"
#include "ridge_regression.h"
#include "vector.h"
#include "allocate.h"
#define VBS 0
#define VB2 0

// char TYPE='C'; // C: return(Cv_fit) M: return(cis_mu) O: OLS L=Use Lambda

int P_FIT=0; // Print fit parameters?
int Nprint=10;
float csi_inf;

static void Score_computation(struct ridge_fit *fit,
			      float *lambda, float **evec,
			      float *yalpha, float *ya2, float *Da,
			      float E0, float Y2, float *A_ref, float *A_inf,
			      int P, int N);
/*static float Error(float *A, float **X, float *Y, int Npar, int N);
static float Get_A(float Lambda, float *yalpha, float **evec, float *eval,
		   float *D, int Npar, int a0);
static float Get_zero_quad(float *x1, float *y1, float *x2, float *y2,
			   float *eval, float *ya2, int Npar);
*/
static float Compute_error(float Lambda, float *lambda, float *ya2,
			   float Y2, int P);
static void Copy_fit(struct ridge_fit *return_fit, struct ridge_fit *fit,
		     int Npar);
float Make_fit(float *A, float Lambda, float *lambda, float **evec,
	       float *yalpha, float *D, int Npar, float **X, float *Y,
	       int N);
static void Print_make_fit(struct ridge_fit *F, float Lambda, char *what,
		    FILE *file_out, char *nameout, float *eval, float **evec,
		    float *yalpha, float *D, float **X, float *Y, float *ya2,
		    float E0, float Y2, float *A_ref, float *A_inf,
		    int Npar, int N, double *sum_x, double var_y);
extern float Corr_coeff(float *xx, float *yy, int N,
			float *slope, float *offset);

static void Print_fit_2(float a0, float a1, int N, char *nameout, char *what,
			float *X);
static void Get_max(float *Max, float *L, float x, float Lambda, int j);
static float New_lambda(float Lambda_old, float *eval, float *ya2,
			int Npar, float e);
static float Lambda_e(float e, float *eval, float *ya2, int Npar);
static float Critical_Lambda(float *eval, float *ya2, int Npar);
static float Old_Cv(float *eval, float *ya2, int Npar);
static float df_dL(float Lambda, float e, float *eval, float *ya2, int Npar);
static float Get_zero_lin(float *x1, float *y1, float *x2, float *y2,
			  float *eval, float *ya2, int Npar);

float Ridge_regression(struct ridge_fit *return_fit, float *Y_pred,
		       float *D_out, char *nameout, float **X, float *Y,
		       int N, int Npar, char TYPE)
{
  int i, a, b;
  if((TYPE!='C')&&(TYPE!='M')&&(TYPE!='O')&&(TYPE!='L')
     &&(TYPE!='0')&&(TYPE!='1')){
    printf("Error, ridge criterion %c is not defined, ", TYPE);
    printf("valid options are C (specific heat), M (max penalty mu) ");
    printf(" O (ordinary) L (give Lambda) and 0 (NoRigid)\n");
    exit(8);
  }
  if(Nprint>Npar)Nprint=Npar; // Number of eigenvalues and D to be printed
  // Compute the covar
  if(VBS)printf("Computing and diagonalizing covariance matrix\n");
  double **C=malloc(Npar*sizeof(double *));
  for(a=0; a<Npar; a++){
    C[a]=malloc(Npar*sizeof(double));
    for(b=0; b<Npar; b++)C[a][b]=0;
  }
  for(i=0; i<N; i++){
    float *x=X[i];
    for(a=0; a<Npar; a++){
      for(b=0; b<=a; b++)C[a][b]+=x[a]*x[b];
    }
  }
  float D[Npar];
  for(a=0; a<Npar; a++){
    D[a]=sqrt(C[a][a]);
    for(b=0; b<=a; b++){
      C[a][b]/=(D[a]*D[b]);
      C[b][a]=C[a][b];
    }
  }

  int nan=0;
  for(a=0; a<Npar; a++){
    if(VBS)printf("Var%d sqrt(C_aa)=%.3g\n", a, D[a]);
    if(D[a]<=0){
      printf("ERROR, variable %d has zero or negative variance\n", a+1);
      nan=1;
    }
  }
  if(nan){
    Empty_matrix_d(C, Npar);
    return(-100);
  }

  float eval[Npar], **evec=malloc(Npar*sizeof(float *));
  for(a=0; a<Npar; a++)evec[a]=malloc(Npar*sizeof(float));
  //f_Diagonalize(Npar, C, eval, evec);
  d_Diagonalize(Npar, C, eval, evec, 1);
  if(VBS){
    printf("lambda= ");
    for(a=0; a<Nprint; a++)printf(" %.4g", eval[a]);
    printf("\n");
  }

  for(a=0; a<Npar; a++)if(isnan(eval[a]))nan++;
  if(nan){
    printf("ERROR, %d eigenvalues are nan\n", nan);
    Empty_matrix_d(C, Npar);
    return(-100);
  }

  // Projection of Y over principal components
  float ya[Npar]; // y_a=sum_i X_ia Y_i/D_a = (x^T Y) where x_ia=X_ia/D_a
  float Y2=0; for(i=0; i<N; i++)Y2+=Y[i]*Y[i];
  double sum;
  for(a=0; a<Npar; a++){
    sum=0; for(i=0; i<N; i++)sum+=X[i][a]*Y[i];
    ya[a]=sum/D[a];
  }
  float yalpha[Npar]; // y_alpha= (x^T Y, u^alpha)
  float ya2[Npar];   // y_alpha2 = (x^T Y, u^alpha)^2
  double yy=0;
  for(a=0; a<Npar; a++){
    sum=0; for(b=0; b<Npar; b++)sum+=evec[a][b]*ya[b];
    yalpha[a]=sum;
    ya2[a]=sum*sum;
    yy+=ya2[a];
  }
  if(VBS)printf("Y2/N= %.3g yy/N= %.3g\n", Y2/N, yy/N);

  // Reference fit parameters
  float A_ref[Npar]; for(a=0; a<Npar; a++)A_ref[a]=ya[a];
  double H2=0, L2=0;
  for(a=0; a<Npar; a++){H2+=ya2[a]; L2+=ya2[a]*eval[a];}
  csi_inf=H2/L2;
  float A_inf[Npar]; for(a=0; a<Npar; a++)A_inf[a]=csi_inf*A_ref[a];

  // Compute fractions, normalize error
  double *sum_x=malloc(Npar*sizeof(double)), sum_y=0, sum_y2=0;
  for(a=0; a<Npar; a++)sum_x[a]=0;
  for(i=0; i<N; i++){
    sum_y+=Y[i]; sum_y2+=Y[i]*Y[i];
    for(a=0; a<Npar; a++)sum_x[a]+=X[i][a];
  }
  sum_y2=(sum_y2-sum_y*sum_y/N);

  // Prepare output
  char namefile[200];
  sprintf(namefile, "%s.dat", nameout);
  FILE *file_out=fopen(namefile, "w");
  if(VBS)printf("Writing %s\n", namefile);
  fprintf(file_out,
	  "# Lambda GenCrossVal Error dE/dL d(LS)/dL Sh LF(csi) LF(inf)\n");
  fprintf(file_out, "# N= %d Npar= %d\n", N, Npar); //RangeRisk
 
  struct ridge_fit fit; fit.A=malloc(Npar*sizeof(float));
  float E0=Compute_error(0, eval, ya2, Y2, Npar);

  // Sequence of values for lambda
  int K=7, Kp=120, Nl=K*(Npar+Kp+1)+1, j=0, k;
  float lambda[Nl], delta, la;
  int NEGLAM=0;
  if(NEGLAM){
    float lambda_min=-2*eval[0];
    for(a=0; a<Npar; a++){
      delta=(-eval[a]-lambda_min)/(K+1);
      la=lambda_min; lambda_min=-eval[a];
      for(k=0; k<K; k++){
	la+=delta; lambda[j]=la; j++;
      }
    }
    delta=-lambda_min/(K+1);
    la=lambda_min;
    for(k=0; k<K; k++){
      la+=delta; lambda[j]=la; j++;
    }
  }
  la=0; lambda[j]=0; j++;
  delta=4*eval[0]/(Kp*K+1);
  for(k=0; k<Kp*K; k++){
    la+=delta; lambda[j]=la; j++;
  }
  Nl=j;

  // Compute
  // Fixed Lambda
  if(TYPE=='L'){
    Print_make_fit(&fit, return_fit->Lambda, "Fix Lambda",
		   file_out, nameout, eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf, Npar,N,sum_x,sum_y2);
    Copy_fit(return_fit, &fit, Npar); goto end;
  }

  // Ordinary Least Squares
  Print_make_fit(&fit, 0.00, "OLS_fit", file_out, nameout,
		 eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		 A_ref,A_inf, Npar,N,sum_x,sum_y2);
  if(TYPE=='O'){Copy_fit(return_fit, &fit, Npar); goto end;}

  //int ini_Sv=0; //ini_Cv2=0, 
  float Sv1_old=-1000; // C2_old=1000; 
  float DUMM=-100;
  float L_GCV=DUMM, L_Cv=DUMM, L_Sv=DUMM, L_R2=DUMM, L_Sh=DUMM; //L_RR=DUMM, 
  float min_GCV=0, max_Cv=0, min_R2=1000, max_Sh=0; ; //max_Cv2=0,min_RR=1000
  float LF_csi_old=0, L_LF_csi=0;
  float max_LF_inf=0, L_LF_inf=0;
  float max_LF_mu=0,  L_LF_mu=0;

  // Loop of Lambda
  for(j=0; j<Nl; j++){
    fit.Lambda=lambda[j];
    //if(fit.Lambda<ev){R_ptr=NULL; S_ptr=NULL;}else{R_ptr=&R2; S_ptr=&Sh;}
    Score_computation(&fit,eval,evec,yalpha,ya2,D,E0,Y2,A_ref,A_inf,Npar,N);

    //if((j==0)||(RR < min_RR)){min_RR=RR; L_RR=Lambda;}
    Get_max(&min_GCV,&L_GCV, -fit.GCV, fit.Lambda, j);
    Get_max(&max_Cv, &L_Cv, fit.Cv, fit.Lambda, j);
    Get_max(&max_Sh, &L_Sh, fit.Sh, fit.Lambda, j);
    Get_max(&min_R2, &L_R2, -fit.R2, fit.Lambda, j);
    Get_max(&max_LF_inf, &L_LF_inf, fit.LF_inf, fit.Lambda, j);
    Get_max(&max_LF_mu, &L_LF_mu, fit.LF_mu, fit.Lambda, j);
    // Select last maximum of LF_csi;
    if(fit.LF_csi > LF_csi_old)L_LF_csi=fit.Lambda; LF_csi_old=fit.LF_csi;

    if(VBS){
      fprintf(file_out,
	      "%5.3f\t%7.4g\t%7.4g\t%9.3g\t%9.3g\t%.4f\t%.4f\t%7.4g\n",
	      fit.Lambda, fit.GCV, fit.E/N, fit.Cv/N,
	      fit.Sv1/N, fit.Sh, fit.LF_csi, fit.LF_inf);
    }

    // Sv1  // Zero of Sv1=d(Lambda S)/dLambda
    if((fit.Lambda>0)&&(Sv1_old <0)&&(fit.Sv1 >0)){
      fit.Lambda=(lambda[j-1]*fit.Sv1-fit.Lambda*Sv1_old)/(fit.Sv1-Sv1_old);
      Score_computation(&fit,eval,evec,yalpha,ya2,D,E0,Y2,A_ref,A_inf,Npar,N);
      L_Sv=fit.Lambda; //LS=Lc*A2c; 
      if(VBS){
	fprintf(file_out,
		"%5.3f\t%7.4g\t%7.4g\t%9.3g\t%9.3g\t%.4f\t%.4f\t%7.4g\n",
		fit.Lambda, fit.GCV, fit.E/N, fit.Cv/N,
		fit.Sv1/N, fit.Sh, fit.LF_csi, fit.LF_inf);
      }
      //ini_Sv=1;
    }
    Sv1_old=fit.Sv1;
  }

  //LS=lambda[Nl-1]*A2;
  //if((Cv2<0)&&(LS > max_Cv2)){max_Cv2=LS; L_Cv2=lambda[Nl-1];} 
  //if((Sv1>0)&&(LS > max_Sv )){max_Sv=LS; L_Sv=lambda[Nl-1];} 

  // Print results of the fits
  fprintf(file_out, "# Normalization of the error: %.2f\n", sum_y2/N);
  fprintf(file_out, "# force=1/a(0) Lambda Rel.error ");
  fprintf(file_out, "Renyi2 Shannon Lambda*F(csi) Lambda*F(inf) fit");
  if(P_FIT)for(a=0; a<Npar; a++)fprintf(file_out, "\tvar%d", a+1);
  fprintf(file_out, "\n");

  if(0){
    // No rigid body
    float X0[N]; for(i=0; i<N; i++)X0[i]=X[i][0];
    float nu=ya[0]/D[0], err=Y2-nu*ya[0]*D[0];
    fprintf(file_out,
	    "#\t%.3g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%s",
	    1./nu, 0.0, err/sum_y2, 1.00, 0.00, 0.00, 0.00, "X1_fit");
    if(P_FIT){
      fprintf(file_out,"\t%.3f", nu);
      for(a=1; a<Npar; a++)fprintf(file_out, "\t0.0");
    }
    fprintf(file_out, "\n");
    if(VB2)Print_fit_2(nu, 0, N, nameout, "X1_fit", X0);
    if(TYPE=='0'){
      fit.Lambda=0; fit.nu=nu; fit.E=err;
      fit.A[0]=nu; for(int a=1; a<Npar; a++)fit.A[a]=0;
      Copy_fit(return_fit, &fit, Npar);
    }
    
    // No rotations
    float a0, a1, r_B=Corr_coeff(X0, Y, N, &a0, &a1);
    double norm=fabs(a0)+fabs(a1);
    double S=fabs(a0)*log(fabs(a0))+fabs(a1)*log(fabs(a1));
    fprintf(file_out,
	    "#\t%.3g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%s",
	    1./a0, 0.0, 1-r_B*r_B,(a0*a0+a1*a1)/(norm*norm), -S/norm+log(norm),
	    0.,0., "X0_X1_fit");
    if(P_FIT){
      fprintf(file_out,"\t%.3f\t%.3f", a0, a1);
      for(a=2; a<Npar; a++)fprintf(file_out, "\t0.0");
    }
    fprintf(file_out, "\n");
    if(VB2)Print_fit_2(a0, a1, N, nameout, "X0_X1_fit", X0);
    if(TYPE=='1'){
      fit.Lambda=0; fit.nu=nu; fit.E=err;
      fit.A[0]=a0; fit.A[1]=a1; for(int a=2; a<Npar; a++)fit.A[a]=0;
      Copy_fit(return_fit, &fit, Npar);
    }
  }

  // Ridge regression fits

  /*// Range risk
  Print_make_fit(&fit, L_RR, "RR_fit", file_out, nameout,
		 eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		 A_ref,A_inf,Npar,N,sum_x,sum_y2);*/

  // Generalized Cross Validation
  Print_make_fit(&fit, L_GCV, "GCV_fit", file_out, nameout,
		 eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		 A_ref,A_inf, Npar,N,sum_x,sum_y2);

 // dE/dLambda
  Print_make_fit(&fit, L_Cv, "Cv_fit", file_out, nameout,
		 eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		 A_ref,A_inf, Npar,N,sum_x,sum_y2);
  if(TYPE=='C'){Copy_fit(return_fit, &fit, Npar); goto end;}

  /*// d^2Cv/dLambda^2
    Print_make_fit(&fit, L_Cv2, "Cv2_fit", file_out, nameout, //d^2 Cv/dLambda^2
    eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N,sum_x,sum_y2);
  */

  // Old Cv
  float Lambda;
  if(0){
    Lambda=Old_Cv(eval, ya2, Npar);
    Print_make_fit(&fit, Lambda, "Oldcv_fit", file_out, nameout,
		   eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf,Npar,N, sum_x, sum_y2);
    Get_max(&max_LF_inf, &L_LF_inf, fit.LF_inf, fit.Lambda, 1);
  }

  // Critical Lambda
  if(0){
    Lambda=Critical_Lambda(eval, ya2, Npar);
    Print_make_fit(&fit, Lambda, "critical_fit", file_out, nameout,
		   eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf,Npar,N, sum_x, sum_y2);
    Get_max(&max_LF_inf, &L_LF_inf, fit.LF_inf, fit.Lambda, 1);
  }

  // Old entropy fit //d(Lambda F)/dLambda
  if(0){
    Print_make_fit(&fit, L_Sv, "entropy_fit",file_out, nameout, 
		   eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf,Npar,N,sum_x,sum_y2);
  }

  if(0){
    int Nfact=10;
    for(i=0; i<Nfact; i++){
      float Lambda=pow(2,i)*0.01*eval[0];
      char name[80]; sprintf(name, "2E%dlmax_fit", i);
      Print_make_fit(&fit, Lambda, name, file_out, nameout,
		     eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		     A_ref,A_inf,Npar,N, sum_x, sum_y2);
      if(fit.Sh > max_Sh){L_Sh=Lambda; max_Sh=fit.Sh;}
      if(fit.R2 < min_R2){L_R2=Lambda; min_R2=fit.R2;}
      //Get_max(&max_LF_csi, &L_LF_csi, fit.LF_csi, fit.Lambda, 1);
      Get_max(&max_LF_inf, &L_LF_inf, fit.LF_inf, fit.Lambda, 1);
      Get_max(&max_LF_mu, &L_LF_mu, fit.LF_mu, fit.Lambda, 1);
    }
  }

  // Entropy csi
  if(0){
    Print_make_fit(&fit, L_LF_csi, "LF_csi_fit", file_out, nameout,
		   eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf,Npar,N,sum_x,sum_y2);
  }

  // Entropy mu
  Print_make_fit(&fit, L_LF_mu, "LF_mu_fit", file_out, nameout,
		 eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		 A_ref,A_inf,Npar,N,sum_x,sum_y2);
  if(TYPE=='M'){Copy_fit(return_fit, &fit, Npar); goto end;}

  // Entropy inf
  if(0){
    Print_make_fit(&fit, L_LF_inf, "LF_inf_fit", file_out, nameout,
		   eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf,Npar,N,sum_x,sum_y2);
  }



  /*/ <l>L,e=4 = L
  float e=4.0, Lambda_old=eval[0]; int it;
  sprintf(name, "%.1fexp_fit", e);
  for(it=0; it<50; it++){
    Lambda=New_lambda(Lambda_old, eval, ya2, Npar, e);
    if(fabs(Lambda-Lambda_old)<0.0001)break;
    Lambda_old=Lambda;
  }
  Print_make_fit(&fit, Lambda, name, file_out, nameout,
		 eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		 A_ref,A_inf,Npar,N, sum_x, sum_y2);
  Get_max(&max_LF_inf, &L_LF_inf, fit.LF_inf, fit.Lambda, 1); */

  if(0){
    // max of the Reniy entropy
    Print_make_fit(&fit, L_R2, "R2_fit", file_out, nameout,
		   eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf,Npar,N, sum_x, sum_y2);
    
    // max of the Shannon entropy
    Print_make_fit(&fit, L_Sh, "Sh_fit", file_out, nameout,
		   eval,evec,yalpha,D,X,Y,ya2,E0,Y2,
		   A_ref,A_inf,Npar,N,sum_x,sum_y2);    
  }

  /*
  int Nfact=7, eprint=0; float l_1=lambda[Nl-1];
  float fact[7]={0.3, 0.5, 0.8, 2.0, 4.0, 10.0, 1000};
  for(i=0; i<Nfact; i++){
    if((eprint==0)&&(fact[i]>1)){
      // entropy fit
      Print_make_fit(Sv1_fit, L_Sv, "entropy_fit",file_out, nameout,
      eval,evec,yalpha,D,X,Y,ya2,E0,Y2,Npar,N,sum_x,sum_y2);
      eprint=1;
    }
    int prn=0; if(VBS &&((i==0)||(i==(Nfact-1))))prn=1;
    float f=fact[i], Lambda=f*L_Sv;
    char name[80]; sprintf(name, "%.2gentr_fit", f);
    Print_make_fit(&fit, Lambda, name, file_out, nameout,
    eval, evec, yalpha, D, X, Y, ya2, E0, Y2,Npar,N,sum_x,sum_y2);
    if(fit.Sh > max_Sh){L_Sh=Lambda; max_Sh=fit.Sh;}
    if(fit.R2 < min_R2){L_R2=Lambda; min_R2=fit.R2;}
  }
  */

  /*// Lambda=sum_a w_a(e) lambda_a, we_a=(y_a)^2/(Lambda+lambda_a)^e/sum
  float e, Lambda_old=0, Lambda; int it;
  for(e=9; e>=-1; e-=2){
    char name[80]; sprintf(name, "%.1fexp_fit", e);
    for(it=0; it<50; it++){
      Lambda=New_lambda(Lambda_old, eval, ya2, Npar, e);
      if(fabs(Lambda-Lambda_old)<0.001)break;
      Lambda_old=Lambda;
    }
    Print_make_fit(&fit, Lambda, name, file_out, nameout,
    eval, evec, yalpha, D, X, Y, ya2, E0, Y2,Npar,N, sum_x, sum_y2);
    if(fit.Sh > max_Sh){L_Sh=Lambda; max_Sh=fit.Sh;}
    if(fit.R2 < min_R2){L_R2=Lambda; min_R2=fit.R2;}
  }
  */

  // Make prediction
 end:
  for(i=0; i<N; i++){
    float *Xi=X[i]; double Yp=0;
    for(a=0; a<Npar; a++)Yp+=Xi[a]*return_fit->A[a];
    Y_pred[i]=Yp;
  }
  float slope, offset, r=Corr_coeff(Y, Y_pred, N, &slope, &offset);
  printf("Ridge regression fit of type %c, Lambda= %.2g\n",
	 TYPE, return_fit->Lambda);
  printf("r(Y_pred, Y_obs)= %.3f slope=%.4g off= %.3g", r, slope, offset);
  printf(" (%d var ridge regr %d samples)\n", Npar, N);
  //printf("Fit error= %.3f\n", return_fit->E/sum_y2);
  for(a=0; a<Npar; a++)D_out[a]=D[a];
  if(VBS){
    printf("D= ");
    for(a=0; a<Nprint; a++)printf(" %.3g", D_out[a]);
    printf("\n");
  }
  return_fit->E/=sum_y2; return_fit->Cv/=sum_y2;

  free(fit.A);
  free(sum_x);
  for(a=0; a<Npar; a++)free(evec[a]); free(evec);
  for(a=0; a<Npar; a++)free(C[a]); free(C);
  fclose(file_out);

  return(r);
}

float Compute_error(float Lambda, float *lambda, float *ya2, float Y2, int P)
{
  double H1=0, H2l=0; int a;
  for(a=0; a<P; a++){
    float laLa=lambda[a]+Lambda;
    float yl=ya2[a]/laLa; H1+=yl;
    H2l+=yl*lambda[a]/laLa;
  }
  return(Y2-(H1*H1/H2l));
}

void Score_computation(struct ridge_fit *fit,
		       float *lambda, float **evec,
		       float *yalpha, float *ya2, float *D,
		       float E0, float Y2, float *A_ref, float *A_inf,
		       int P, int N)
{
  int a;
  double H1=0, H2=0, H3=0, H2l=0, H3l=0, L1=0; 
  for(a=0; a<P; a++){
    float laLa=lambda[a]+fit->Lambda;   
    L1+=1./laLa;
    float yl=ya2[a]/laLa;  H1+=yl;
    yl/=laLa; H2+=yl;
    H2l+=yl*lambda[a];
    yl/=laLa; H3+=yl;
    H3l+=yl*lambda[a];
  }
  float csi=H2/H2l;
  //fit->nu=1.+fit->Lambda*csi; 
  fit->nu=H1/H2l; 
  fit->E=Y2-H1*fit->nu;
  fit->Cv=fit->Lambda*H1*(H1*H3-H2*H2)/(H2l*H2l);
  fit->A2=H2*fit->nu*fit->nu;
  fit->Sv1=fit->Lambda*H3-H3l;

  //*Cv2=yL5;
  //double E_nu=Y2-(H1+Lambda*H2);
  float Trace=(N-P+fit->Lambda*L1);
  float E_nu=Y2-2*H1+H2l;
  fit->GCV=N*E_nu/(Trace*Trace);
  fit->RR=((N-P)*E_nu/E0-2*Trace)/N;

  float A_ns[P]; 
  for(a=0; a<P; a++)A_ns[a]=yalpha[a]/(lambda[a]+fit->Lambda);
  for(a=0; a<P; a++){
    double sum=0; for(int b=0; b<P; b++)sum+=A_ns[b]*evec[b][a];
    fit->A[a]=fit->nu*sum;
  }

  double sum_csi=0, sum_inf=0;
  for(a=0; a<P; a++){
    float d=fit->A[a]-csi*A_ref[a]; sum_csi+=d*d;
    d=fit->A[a]-A_inf[a]; sum_inf+=d*d;
  }
  fit->LF_csi=fit->Lambda*sum_csi;
  fit->LF_inf=fit->Lambda*sum_inf;
  float mu_inv=1./(1+2*fit->Lambda*(csi-csi_inf));
  fit->LF_mu=fit->LF_inf*mu_inv;

  double norm=0, R2=0, Sh=0;
  for(a=0; a<P; a++){
    float x=fabs(fit->A[a]); R2+=x*x; if(x)Sh+=x*log(x); norm+=x;
  }
  fit->R2=R2/(norm*norm);
  fit->Sh=-Sh/norm+log(norm);

}

float Make_fit(float *A, float Lambda, float *lambda, float **evec,
	       float *yalpha, float *D, int Npar, float **X, float *Y, int N)
{
  // Compute A
  int a, b, i; double sum;
  float coeff[Npar];
  for(a=0; a<Npar; a++)coeff[a]=yalpha[a]/(lambda[a]+Lambda);
  double H1=0, H2l=0;
  for(a=0; a<Npar; a++){
    float yl=yalpha[a]*coeff[a];
    H1+=yl; H2l+=yl*lambda[a]/(lambda[a]+Lambda);
  }
  float nu=H1/H2l;
  for(a=0; a<Npar; a++){
    sum=0; for(b=0; b<Npar; b++)sum+=coeff[b]*evec[b][a];
    A[a]=nu*sum/D[a];
  }

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

/*float Get_A(float Lambda, float *yalpha, float **evec, float *lambda,
              float *D, int Npar, int a0)
{
  double H1=0, H2l=0;  int a;
  for(a=0; a<Npar; a++){
    float yl=yalpha[a]*yalpha[a]/(lambda[a]+Lambda);
    H1+=yl; H2l+=yl*lambda[a]/(lambda[a]+Lambda);
  }
  float nu=H1/H2l; double A0=0;
  for(a=0; a<Npar; a++)A0+=yalpha[a]*evec[a][a0]/(lambda[a]+Lambda);
  return(nu*A0);
  }
*/


void Print_make_fit(struct ridge_fit *F, float Lambda, char *what,
		    FILE *file_out, char *nameout, float *eval, float **evec,
		    float *yalpha, float *D, float **X, float *Y, float *ya2,
		    float E0, float Y2, float *A_ref, float *A_inf,
		    int Npar, int N, double *sum_x, double var_y)
{
  F->Lambda=Lambda;
  Score_computation(F, eval,evec,yalpha,ya2,D,E0,Y2,A_ref,A_inf,Npar,N);
  int a; for(a=0; a<Npar; a++)F->A[a]/=D[a];
  /*double sum=0; for(a=0; a<Npar; a++)sum+=F->A[a]*sum_x[a];
  F->dof[0]=F->A[0]*sum_x[0]/sum; // Internal fraction
  F->dof[1]=F->A[1]*N/sum;  // Translations
  F->dof[2]=1-F->dof[0]-F->dof[1]; // Rotations*/
  // Print
  fprintf(file_out,
	  "#\t%.3g\t%.3g\t%.4g\t%.3g\t%.4g\t%.4g\t%.4g\t%12s",
	  1./F->A[0], F->Lambda, F->E/var_y, F->R2, F->Sh,
	  F->LF_csi, F->LF_inf, what);
  if(P_FIT)for(a=0; a<Npar; a++)fprintf(file_out, "\t%.3g", F->A[a]);
  fprintf(file_out, "\n");

  if(VB2){
    char namefile[200];
    sprintf(namefile, "%s_%s.dat", nameout, what);
    FILE *file_o=fopen(namefile, "w"); int i, a;
    printf("Writing %s fit in %s\n", what, namefile);
    fprintf(file_o, "#A= ");
    for(i=0; i<Npar; i++)fprintf(file_o, " %.3g",F->A[i]);
    fprintf(file_o, "\n# Y ");
    for(i=0; i<Npar; i++)fprintf(file_o, " X%d", i+1);
    fprintf(file_o, "\n");
    fprintf(file_o, "#nu= %.3f\n", F->nu);
    fprintf(file_o, "#eval= ");
    for(i=0; i<Npar; i++)fprintf(file_o, " %.3g",eval[i]);
    fprintf(file_o, "\n");
    fprintf(file_o, "#Complete only_internal\n");
    for(i=0; i<N; i++){
      double Y=0;
      for(a=0; a<Npar; a++)Y+=F->A[a]*X[i][a];
      fprintf(file_o, "%.3f %.3f\n", Y, F->A[0]*X[i][0]);
    }
    fclose(file_o);
  }
}

void Print_fit_2(float a0, float a1, int N, char *nameout, char *what,
		 float *X)
{
  char namefile[200]; sprintf(namefile, "%s_%s.dat", nameout, what);
  FILE *file_out=fopen(namefile, "w"); int i;
  printf("Writing %s fit in %s\n", what, namefile);
  fprintf(file_out, "# slope= %.3g offset= %.3g\n", a0, a1);
  fprintf(file_out, "#Complete only_internal\n");
  for(i=0; i<N; i++){
    float Y=a0*X[i];
    fprintf(file_out, "%.3f %.3f\n", Y+a1, Y);
  }
  fclose(file_out);
 }

/*float Error(float *A, float **X, float *Y, int Npar, int N)
{
  double E=0; int i, a;
  for(i=0; i<N; i++){
    float *Xi=X[i]; double Yp=0;
    for(a=0; a<Npar; a++)Yp+=Xi[a]*A[a];
    Yp-=Y[i]; E+=Yp*Yp;
  }
  return(E);
  }*/

float New_lambda(float Lambda_old, float *eval, float *ya2, int Npar, float e)
{
  double Lambda=0, norm=0; int a; 
  for(a=0; a<Npar; a++){
    float w=ya2[a]/pow(eval[a]+Lambda_old, e);
    Lambda+=eval[a]*w; norm+=w;
  }
  return(Lambda/norm);
}

float Old_Cv(float *eval, float *ya2, int Npar){
  float Lambda_old=0, Lambda; int it;
  for(it=0; it<50; it++){
    Lambda=0.5*New_lambda(Lambda_old, eval, ya2, Npar, 4.0);
    if(fabs(Lambda-Lambda_old)<0.001)return(Lambda);
    Lambda_old=Lambda;
  }
  printf("WARNING, Lambda not found in 50 steps, e= 4.0\n");
  return(Lambda);
}

float Critical_Lambda(float *eval, float *ya2, int Npar){
  // Lambda=sum_a w_a(e) lambda_a, we_a=(y_a)^2/(Lambda+lambda_a)^e/sum
  // Find the e where dLambda/de diverges => df/dLambda=0

  double FACT=1.1; float e=3; int it;
  float Lambda=Lambda_e(e, eval, ya2, Npar);
  float dfdL=df_dL(Lambda, e, eval, ya2, Npar);
  if(VBS)printf("df_dL= %.3f  e= %.4f L= %.3f\n", dfdL, e, Lambda);
  float e1=e, Lambda_ret=Lambda, dfdL_1=dfdL;

  for(it=0; it<50; it++){
    e=e1*FACT; Lambda=Lambda_e(e, eval, ya2, Npar);
    dfdL=df_dL(Lambda, e, eval, ya2, Npar);
    if(VBS)
      printf("df_dL= %.3f  e= %.4f e1= %.3f L= %.3f\n", dfdL,e,e1,Lambda);
    if(dfdL <= 0){
      Lambda_ret=Get_zero_lin(&e1, &dfdL_1, &e, &dfdL, eval, ya2, Npar);
      FACT=sqrt(FACT);
    }else if(dfdL < dfdL_1){
      e1=e; dfdL_1=dfdL; Lambda_ret=Lambda;
    }else{
      FACT=sqrt(FACT);
      if(FACT < 1.00000001)return(Lambda_ret);
    }
    if(fabs(dfdL_1) < 0.01)return(Lambda_ret);
  }
  printf("WARNING, critical e not found in 50 steps\n");
  return(Lambda_ret);
}

float Lambda_e(float e, float *eval, float *ya2, int Npar){
  float Lambda_old=0, Lambda; int it, Nt=100;
  for(it=0; it<Nt; it++){
    Lambda=New_lambda(Lambda_old, eval, ya2, Npar, e);
    if(fabs(Lambda-Lambda_old)<0.001)return(Lambda);
    Lambda_old=Lambda;
  }
  if(VBS)printf("WARNING, Lambda not found in %d steps, e= %.3f\n", Nt, e);
  return(Lambda);
}

float df_dL(float Lambda, float e, float *eval, float *ya2, int Npar)
{
  float e1=e+1; double f1=0, f0=0; int a; 
  for(a=0; a<Npar; a++){
    float w=ya2[a]/pow(eval[a]+Lambda, e1);
    f1+=eval[a]*w; f0+=w;
  }
  float df_dL=f1*e1-Lambda*f0*(e-1);
  return(df_dL);
}

float Get_zero_lin(float *x1, float *y1, float *x2, float *y2,
		   float *eval, float *ya2, int Npar)
{
  float e=*x1-(*y1)*((*x1)-(*x2))/((*y1)-(*y2));
  float Lambda=Lambda_e(e, eval, ya2, Npar);
  float y=df_dL(Lambda, e, eval, ya2, Npar);
  if((y < *y1)&&(y > 0)){*x1=e; *y1=y;}
  if(VBS)printf("df_dL= %.2f  e= %.3f e1= %.3f L= %.3f\n", y, e, *x1, Lambda);
  return(Lambda);
}

/*float Get_zero_quad(float *x1, float *y1, float *x3, float *y3,
		    float *eval, float *ya2, int Npar)
{
  float x2=((*x1)+(*x3))/2;
  float Lambda=Lambda_e(x2, eval, ya2, Npar);
  float y2=df_dL(Lambda, x2, eval, ya2, Npar);
  if(VBS)printf("df_dL= %.2f  e= %.3f e1= %.3f L= %.3f\n", y2,x2,*x1,Lambda);
  if(y2 <=0)return(Get_zero_lin(x1, y1, &x2, &y2, eval, ya2, Npar));
  float dy2=(*y1-y2)/(*x1-x2);
  float dy3=(*y1-*y3)/(*x1-*x3);
  float a=(dy2-dy3)/(x2-*x3);
  float b=dy2-a*(*x1+x2);
  float c=*y1-a*(*x1)*(*x1)-b*(*x1);
  float x0=(-b+sqrt(b*b-4*a*c))/(2*a);
  Lambda=Lambda_e(x0, eval, ya2, Npar);
  float y0=df_dL(Lambda, x0, eval, ya2, Npar);
  if(VBS)printf("df_dL= %.2f  e0= %.3f e1= %.3f L= %.3f\n", y0,x0,*x1,Lambda);
  if(y0 <=0)return(Get_zero_lin(x1, y1, &x0, &y0, eval, ya2, Npar));
  if(y0 < *y1){*y1=y0; *x1=x0;}
  return(Lambda);
}
*/

void Get_max(float *Max, float *L, float x, float Lambda, int j){
    if((j==0)||(x > *Max)){*L=Lambda; *Max=x;}
}

void Copy_fit(struct ridge_fit *return_fit, struct ridge_fit *fit, int N){
  float *A_tmp=return_fit->A; int i;
  *return_fit=*fit; return_fit->A=A_tmp;
  for(i=0; i<N; i++)return_fit->A[i]=fit->A[i];
}

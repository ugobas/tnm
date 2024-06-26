#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ridge_regression.h"
#include "vector.h"
#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "nma_para.h"
#include "Fit_B.h"


int VBS=0;
float Compute_correlation(float **Cart_mode, float *sigma2,
			  int N_modes, int i1, int i2);
void Print_data(float **X, float *Y, int N, int Npar, char *nameout);
int Find_outliers_Z(int *outlier, float *y, int N, float sig);
int Find_outliers_mean(int *outlier, float *y, int N, float sig);


float Fit_fluctuations(int FIT_B, float *B_TNM,
		       float *RMSD_NM, char *out,
		       float *Temp_comp, float Temp,
		       int anharmonic, struct Normal_Mode NM,
		       struct Reference Ref_kin, char *REF,
		       struct residue *seq1, int nres1, int nmr1,
		       atom *atoms1, char *nameout1,
		       float *factor_B, int *outlier_B)
{

  float *sigma2; char nameout[100];
  if(anharmonic==0){
    sigma2=NM.sigma2;
    sprintf(nameout, "%s", nameout1);
  }else{
    sigma2=NM.sigma2_anhar;
    sprintf(nameout, "%s.anharmonic", nameout1);
  }

  int i, k, N_ref=Ref_kin.N_ref;


  double bsum=0, bsum2=0; // mass_CA=0; 
  double pi=3.1415, norm_B=8.0*pi*pi/3;
  char BSEL[4]="CA";
  if(strcmp(REF, "CB")==0)strcpy(BSEL, "CB");
  int Na=0, j=0, Bref[nres1];
  float Xref[(3*nres1)];
  float *B_exp=malloc(N_ref*sizeof(float)), *B_pred_all=NULL;
  for(i=0; i<N_ref; i++){
    atom *atm=atoms1+Ref_kin.atom_num[i];
    if((strncmp(atm->name, BSEL, 2)!=0)&&
       ((seq1[atm->res].amm!='G')||(strncmp(atm->name, "CA", 2)!=0)))continue;
    Bref[Na]=Ref_kin.atom_num[i];
    B_TNM[Na]=Compute_correlation(NM.Cart,sigma2,NM.N,i,i);
    //if(ANM)B_TNM[Na]/=Ref_kin.mass_atom[i];
    //mass_CA+=Ref_kin.mass_atom[i];
    B_exp[Na]=atm->B_factor/norm_B; //
    bsum+=B_exp[Na]; bsum2+=B_exp[Na]*B_exp[Na];
    for(k=0; k<3; k++){Xref[j]=atm->r[k]; j++;}
    Na++;
  }
  bsum/=Na; bsum2=(bsum2/Na-bsum*bsum);
  float RMSD_crys=sqrt(bsum), RMSD_crys_int=-1;
  printf("%d %s atoms selected in %d residues\n", Na, BSEL, nres1);
  printf("RMSD of fluctuations (B fact tot)= %.3f\n", RMSD_crys);

  char ridge='C'; //'M'
  float bsmall=0.01, r_B=0, slope_B=0, slope_NoRot=0, f_dof[3]; 
  f_dof[0]=0; f_dof[1]=0; f_dof[2]=0;
  for(i=0; i<nres1; i++)outlier_B[i]=0; //modyves
  if((bsum < bsmall)||(bsum2<bsmall)||(nmr1)){
    //YYY this is a weird contstraint... ditches many cases..
    printf("WARNING, no Bfact fit bsum= %.2g var= %.2g nmr=%d\n",
	   bsum,bsum2,nmr1);
    free(B_exp);
    B_exp=NULL;
    RMSD_crys=0;
    // FIT_B=-1; // YYY this is probably necessary ?? Ugo: This is wrong
  }else{
    // Fit B factors
    // Fit without rotations
    float offset=0, r=Corr_coeff(B_TNM, B_exp, Na, &slope_NoRot, &offset);
    printf("r(B_pred,B_exp)= %.3f offset= %.3f force const factor=%.4g",
	   r, offset, 1./slope_NoRot);
    printf(" (2 var fit)\n");
    // Fit with rotations
    B_pred_all=malloc(nres1*sizeof(float));
    slope_B=Bfactors_fit(&r_B, outlier_B, f_dof, B_pred_all, B_exp, B_TNM,
			 Xref, Na, nameout1, ridge);
    printf("B-factor fit with rotations, force const factor= %.3g\n", slope_B);
    RMSD_crys_int=sqrt(RMSD_crys*RMSD_crys*f_dof[0]);
    if((r_B<R_MIN)||(slope_B<SLOPE_MIN)){
      printf("WARNING, r_B= %.2f or Bfact slope= %.2g too small\n",
	     r_B,slope_B); 
      if(FIT_B>0)FIT_B=0;
    }
  }

  /*if(FIT_B && ANISOU){
    float **aniso_pred=malloc(N_ref*sizeof(float **));
    Compute_anisou(aniso_pred, N_modes, N_ref, eigen_B, Cart_mode);
    int weight[N_ref];
    float dot_aniso=0, ov_aniso=0, delta_aniso=2.;
    float **aniso_exp=malloc(N_ref*sizeof(float **));
    int n_aniso=Compare_anisou(aniso_pred, aniso_exp, nca, model, inter,
    file_aniso, prot_name, N, weight,
    &dot_aniso, &ov_aniso, &delta_aniso);
    fprintf(file_sum, " %3d %.3f %6.3f %.4f",
    n_aniso, ov_aniso, dot_aniso, delta_aniso);
    }*/

  // factor for simulating ensemble
  if(slope_B){*factor_B=sqrt(slope_NoRot/slope_B);}
  else{*factor_B=-1;}

  double sum_sigma2=0;
  for(i=0; i<NM.N; i++)if(sigma2[i]>0)sum_sigma2+=sigma2[i];
  float kappa=1;
  if((FIT_B>0)&&(slope_B>0)){
    kappa=1./slope_B;
    printf("Force constant: %.4g x %.4g = %.4g    -> K0 = %.4g"
	   "  (Bfactor fit)\n", kappa,KAPPA,kappa*KAPPA,kappa*K0);
  }else if(FIT_B<0){
    kappa= sum_sigma2/(Ref_kin.mass_sum*RMSD_EXP*RMSD_EXP);
    printf("Force constant: %.4g x %.4g = %.4g     -> K0 = %.4g"
	   "   (<RMSD>= %.2f)\n", kappa,KAPPA,kappa*KAPPA,kappa*K0,RMSD_EXP);
  }else { //(FIT_B==0)||(slope_B<=0)
    kappa=1;
    printf("Force constant: %.4g x %.4g = %.4g    -> K0 = %.4g"
	   "    (Default value)\n", kappa,KAPPA,kappa*KAPPA,K0);
  }
  *RMSD_NM=sqrt(sum_sigma2/(Ref_kin.mass_sum*kappa));
  printf("RMSD of normal modes= %.3f (pred)\n", *RMSD_NM);
  if(FIT_B>0 && Temp>0){
    *Temp_comp=Temp;
  }else{
    //*Temp_comp=TEMP_DEF*kappa*KAPPA/KAPPA_DEF;
    *Temp_comp=TEMP_DEF;
  }

  if(anharmonic){sprintf(out,   "#B-factors, anharmonic analysis\n");}
  else{sprintf(out,   "#B-factors, harmonic analysis\n");}

  char tmp[60];
  sprintf(tmp, "r(B_pred,B_exp)       %.3f\n", r_B); strcat(out, tmp);
  sprintf(tmp, "Force constant factor %.3f\n", kappa); strcat(out, tmp);
  sprintf(tmp, "RMSD_crystal          %.3f (all)\n", RMSD_crys);
  strcat(out, tmp);
  sprintf(tmp, "RMSD_crystal          %.3f (internal)\n", RMSD_crys_int);
  strcat(out, tmp);
  sprintf(tmp, "RMSD_Norm.Modes       %.3f\n", *RMSD_NM); strcat(out, tmp);
  /*sprintf(tmp, "Estimated temperature %.1f (kappa= %.3g)\n",
   *Temp_comp, kappa*KAPPA); strcat(out, tmp);*/
  sprintf(tmp, "Internal_fluct        %.3f\n", f_dof[0]); strcat(out, tmp);
  sprintf(tmp, "Translation_fluct     %.3f\n", f_dof[1]); strcat(out, tmp);
  sprintf(tmp, "Rotation_fluct        %.3f\n", f_dof[2]); strcat(out, tmp);

  //  Print B factors
  // Rescale predicted B factors
  for(i=0; i<nres1; i++)B_TNM[i]/=kappa;
  Print_B_fact(B_TNM, B_pred_all, B_exp, Na, atoms1, Bref,
	      seq1, nameout, "Bfact", r_B, slope_B, f_dof, ridge);
  if(B_exp)free(B_exp);
  if(B_pred_all)free(B_pred_all);
  printf("FIT_B= %d kappa= %.3g\n", FIT_B, kappa);

  return(kappa);
}

float Bfactors_fit(float *r, int *outlier, float *dof,
		   float *B_pred_all, float *B_PDB, float *B_ENM,
		   float *R_eq, int N, char *name, char type)
{
  /* It fits B factors as the combination of the internal motion
     predicted by the ENM and the rigid body motion dr_i = t+I X r_i,
     which leads to
     B^exp_i ~ A0 B^ENM_i + A1(translations) +
               A2*r_ix+A3*r_iy+A+A4*r_iz (rototranslation) +
	       A5*r_ix*r_ix+...A10*r_iz*r_iz (rotations)
     The 11 free parameters A0... A10 are fitted through ridge regression.

     OUTPUT file <name>_Bfact_fit.dat
     For each examined Lambda we report:
     Lambda, RangeRisk, GenCrossVal(Golub_etal_1979), Error, dE/dLambda, A0
     The program outputs results of three kinds of fit:
     (1) B^exp_i ~ A0 B^ENM_i + A1 (only translations, fit for Lambda=0)
     (2) 11 parameter fit for Lambda that minimizes RangeRisk (Mallows, 1973,
     as reported in Golub, Heath & Wahba 1979);
     (3) 11 parameters, for Lambda that maximizes the "specific heat"
     dE/dLambda (Bastolla, unpublished)
     For each fit (1), (2), (3) we report:
     Force constant (A0), Lambda, Relative error of the fit,
     Fraction of B^exp atributed to internal motion, fraction attributed to
     rotation and rototranslation.
   */


  int Npar=11, i, j=0;
  float **Xall=malloc(N*sizeof(float *));
  for(i=0; i<N; i++){
    Xall[i]=malloc(Npar*sizeof(float));
    float *x=Xall[i];
    x[0]=B_ENM[i]; x[1]=1;
    float rx=R_eq[j],ry=R_eq[j+1],rz=R_eq[j+2]; j+=3;
    x[2]=rx; x[3]=ry; x[4]=rz;
    x[5]=rx*rx; x[6]=rx*ry; x[7]=rx*rz;
    x[8]=ry*ry; x[9]=ry*rz; x[10]=rz*rz;
  }
  char nameout[100];
  sprintf(nameout, "%s.Bfact_fit", name);
  if(VBS)Print_data(Xall, B_PDB, N, Npar, name);

  // Eliminate outliers
  for(i=0; i<N; i++)outlier[i]=0;
  int num=N; float sigma=2.5;
  // num-=
  // Removed UGO 03/02/2021
  // It is dangerous to remove outliers in the normal modes because this may
  // prevent from fitting some modes that only move outliers
  num-=Find_outliers_mean(outlier, B_PDB, N, sigma);
  /*float sigma=3;
  num-=Find_outliers_Z(outlier, B_ENM, N, sigma); 
  num-=Find_outliers_Z(outlier, B_PDB, N, sigma);*/

  printf("%d outliers eliminated for B factors fit\n", N-num);
  float Y[num]; int k=0;
  float **X=malloc(num*sizeof(float *));
  for(i=0; i<N; i++){
    if(outlier[i])continue;
    X[k]=malloc(Npar*sizeof(float));
    float *x=X[k], *z=Xall[i];
    for(int a=0; a<Npar; a++)x[a]=z[a];
    Y[k]=B_PDB[i]; k++;
  }

  // Complete fit:
  // char type='M'; C= specific heat M= Max penalty
  printf("Fitting B factors with internal and rigid body motions\n");
  float D_out[Npar];
  struct ridge_fit fit; fit.A=malloc(Npar*sizeof(float));
  *r=Ridge_regression(&fit, B_pred_all, D_out, nameout,
		      X, Y, num, Npar, type);

  // Returned results
  int a; double sum_x[Npar];
  for(a=0; a<Npar; a++)sum_x[a]=0;
  for(i=0; i<N; i++){
    float *Xi=Xall[i]; double Yp=0; 
    for(a=0; a<Npar; a++){
      Yp+=Xi[a]*fit.A[a]; sum_x[a]+=Xi[a];
    }
    B_pred_all[i]=Yp;
  }
  float slope=fit.A[0];
  // Degrees of freedom
  double sum=0; for(a=0; a<Npar; a++)sum+=fit.A[a]*sum_x[a];
  dof[0]=fit.A[0]*sum_x[0]/sum; // Internal fraction
  dof[1]=fit.A[1]*N/sum;  // Translations
  dof[2]=1-dof[0]-dof[1]; // Rotations*/
  free(fit.A);
  for(i=0; i<N; i++){free(Xall[i]);} free(Xall);
  for(i=0; i<num; i++){free(X[i]);} free(X);
  return(slope);
}

void Print_data(float **X, float *Y, int N, int Npar, char *name){
  char nameout[100];
  sprintf(nameout, "%s_Bfact_data_fits.dat", name);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  int a,i;
  /*float *norm=malloc(Npar*sizeof(float)); 
  for(a=0; a<Npar; a++){
    double sum=0;
    for(i=0; i<N; i++)sum+=X[i][a]*X[i][a];
    norm[a]=sqrt(sum/N);
    }*/
  fprintf(file_out, "#Y ENM const x y z xx xy xz yy yz zz\n");
  for(i=0; i<N; i++){
    fprintf(file_out, "%.3f\t",Y[i]);
    for(a=0; a<Npar; a++)fprintf(file_out, "%.3f\t",X[i][a]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
  //free(norm);
}

int Find_outliers_Z(int *outlier, float *y, int N, float sig){
  double y1=0, y2=0; float *yy=y; int i;
  for(i=0; i<N; i++){
    y1+=*yy; y2+=(*yy)*(*yy); yy++;
  }
  y1/=N; y2=y2-N*y1*y1; if(N>1)y2=sqrt(y2/(N-1));
  float max=y1+sig*y2, min=y1-sig*y2; int n=0;
  for(i=0; i<N; i++){
    if(outlier[i])continue;
    if((y[i]>max)||(y[i]<min)){outlier[i]=1; n++;}
  }
  return(n);
}

int Find_outliers_mean(int *outlier, float *y, int N, float sig){
  double y1=0; int i;
  for(i=0; i<N; i++)y1+=y[i];
  float max=sig*y1/N; int n=0;
  for(i=0; i<N; i++){
    if(outlier[i])continue;
    if(y[i]>max){outlier[i]=1; n++;}
  }
  return(n);
}

float Fit_no_outliers(int *outlier, float *XX, float *YY, int N,
		      float *slope, float *offset)
{
  int k=0; float X[N], Y[N];
  for(int i=0; i<N; i++){
    if(outlier[i])continue;
    X[k]=XX[i]; Y[k]=YY[i]; k++;
  }
  float r=Corr_coeff(X, Y, k, slope, offset);
  return(r);
}

int Filter_modes_outliers(struct Normal_Mode NM, struct Reference Ref_kin,
			  float t_out, float t_mode)
{
  int discard_tot=0;
  int i, N_ref=Ref_kin.N_ref;
  float B_TNM[N_ref];
  int outlier[N_ref]; 

  for(int iter=0; iter<10; iter++){
    for(i=0; i<N_ref; i++){
      B_TNM[i]=Compute_correlation(NM.Cart,NM.sigma2,NM.N,i,i);
    }
    for(i=0; i<N_ref; i++)outlier[i]=0;
    int nout=Find_outliers_mean(outlier, B_TNM, N_ref, t_out);
    printf("%d iter %d outliers with fluctuations > %.2f of the mean\n",
	   iter, nout, t_out);
    int discard=0;
    for(int k=0; k<NM.N; k++){
      if((NM.select[k]==0)||(NM.sigma2[k]==0))continue;
      double B_out=0, B_no_out=0;
      for(int j=0; j<3*N_ref; j++){
	float v2=NM.Cart[k][j]; v2*=v2; j++;
	i=j/3; if(outlier[i]){B_out+=v2;}else{B_no_out+=v2;}
      }
      if(B_out>t_mode*(B_out+B_no_out)){
	printf("Discard: mode %d, sigma2= %.3g Collectivity= %.3g ",
	       k, NM.sigma2[k], NM.Cart_coll[k]);
	printf("Outliers= %.3g No outliers= %.3g\n", B_out, B_no_out);
	NM.select[k]=0; NM.sigma2[k]=0; discard++;
      }
    }
    discard_tot+=discard;
    printf("%d iter %d modes discarded for outliers represent > %.2f motion\n",
	   iter, discard_tot, t_mode);
    if(discard==0)break;
  }
  return(discard_tot);
}

float Compute_correlation(float **Cart_mode, float *sigma2,
			  int N_modes, int i1, int i2)
{
  int ia, k1=3*i1, k2=3*i2;
  double corr=0;
  for(ia=0; ia< N_modes; ia++){
    if(sigma2[ia]==0)continue;
    float *x1=Cart_mode[ia]+k1;
    float *x2=Cart_mode[ia]+k2;
    double xx=x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
    corr += xx*sigma2[ia];
  }
  return(corr);
}

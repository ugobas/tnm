#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "nma_para.h"
#include "allostery.h"
#include "mutation.h"
#include "contacts.h"
#include "vector.h"
#include "atom_numb.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "allocate.h"
#include "contacts.h"
#include "optimization.h"
#include "externals.h"
#include "ridge_regression.h"


float RC=3.25; // length scale for rescaling mutation parameters
float REXP=1;  // REXP=RC^EXP_F C_size=C_SIZE*KAPPA*REXP
static int OPT_COEFF=0;
static char OPT_SCORE='C';
// C= Cartesian correlation; A=absolute value correlation;
//static float Lambda=0.10; // Penalization in optimization
static float EXP_FORCE=0, Exp;
float COLL_THR_MUT=0.01; // Minimal collectivity of normal modes to predict mut
//float C_SIZE=0.015, C_STAB=1.4, C_DIST=1.0;
float C_SIZE=0.36, C_STAB=0.22, C_DIST=4.81;
int PRINT_FLANK=1;
float scale_C;

int ncmax=100; 
int mut_cont[100], wt_cont[100];
int mut_res[100], wt_res[100];

float Econt_norm[21][21], D_Res_norm[21][21], natm_norm[21];
int ini_fit=1;
float Lambda_RRR=0;
int amax=3;

// Econt D_Residues atom_numb

void Read_pred_para(char *name_in);
void Normalize_parameters();
float Pred_mut_str(float C_size, float C_stab, float C_dist,
		   float *DE, float *MSD_mut, float *MSD_cross,
		   float *A_mut, float *A_nomut, float *A_cross, float *offset,
		   float *pred_mut_3, float *pred_mut_2,
		   float *pred_mut_tot_2,
		   int Na, int *n_cont_mut, float *Force, float *Force_coeff,
		   float *pred_mut_Cart, int N_Cart, 
		   float *pred_mut_Tors, int N_axes,
		   char *AAwt, int *Posmut, char *AAmut, int Nmut,
		   int **clist, int **cnum, struct interaction *Int_list,
		   struct residue *seq, int Nres, atom *atoms, int *resatom,
		   struct Reference Ref_kin,
		   int *kin_atom, struct Normal_Mode NM,
		   float *sigma2, float **B3, float *B_pred,
		   float *str_diff_3, float *str_diff_2);
float Optimize_coeff(float *C_size, float *C_stab, float *C_dist,
		     float *MSD_mut, float *MSD_cross,
		     float *A_mut,float *A_nomut,float *A_cross,float *offset,
		     float *pred_mut_3, float *pred_mut_2,
		     float *pred_mut_tot_2,
		     int Na, int *n_cont_mut, float *Force, float *Force_coeff,
		     float *pred_mut_Cart, int N_Cart, 
		     float *pred_mut_Tors, int N_axes,
		     char *AAwt, int *Posmut, char *AAmut, int Nmut,
		     int **clist, int **cnum, struct interaction *Int_list,
		     struct residue *seq, int Nres, atom *atoms, int *resatom,
		     struct Reference Ref_kin,
		     int *kin_atom, struct Normal_Mode NM, float *sigma2,
		     float **B3, float *B_pred,
		     float *str_diff_3, float *str_diff_2); 

int Mutation_force(float *Force, int N_Cart,
		   float C_size, float C_stab, float C_dist, 
		   char *AAwt, int *Posmut, char *AAmut, int Nmut,
		   int **clist, int **cnum, struct interaction *Int_list,
		   struct residue *seq, int nres, atom *atoms,
		   struct Reference Ref_kin, int *kin_atom);
int Mut_neighbor(int *mut_neigh, int Na,
		 char *AAwt, int *Posmut, char *AAmut, int Nmut,
		 int **clist, int **cnum, struct interaction *Int_list,
		 atom *atoms, int nres);
int Extract_reference_atoms(int *resatom, atom *atoms, int natoms,
			    struct residue *seq, int Nres);

void Compute_mut_def(float *pred_mut_Cart, int N_Cart,
		     float *pred_mut_Tors, int N_axes,
		     float *DE, float *Force_coeff, float *Force,
		     struct Normal_Mode NM, float *sigma2, float mass);

float *Convert_to_tors(float *Cart, int N3, float **J_ar, int N_axes);
void  Match_atoms(int *ref_atom, atom *atoms, int natoms,
		  int *atom_ref, int N_ref);
float Cosine(float *scale, float *xx, float *yy, int n);
extern float Fit_2(float *Y, float *X0, float *X1, int N,
		   float *A0, float *A1, float *b, int first);
extern float Fit_1(float *Y, float *X, int N, float *A, float *b);
extern void f_sort(float d[], int n, int *i_rank);
extern void f_Diagonalize(int N, float **MATRIX, float *eigen_values,
			  float **eigen_vector, int SIGN);

int Mutation(float *mut_Tors_out, int N_axes, float *mut_CC_out,
	     float KAPPA, struct Reference Ref_kin, 
	     int Nmut, char *AAwt, int *Posmut, char *AAmut,
	     //int Force_type, 0=size
	     struct Normal_Mode NM, atom *atoms, int natoms,
	     struct interaction *Int_list, int N_int,
	     struct residue *seq,
	     float *B_pred,
	     float *Confchange_Cart, struct Reference Ref_CC,
	     float *Confchange_Tors, float *Tors_fluct,
	     char *nameout1, int anhar, char *Mut_para)
{
  // Determine mutations
  if(Nmut==0){
    printf("Sorry, no mutation is present to predict its effect\n");
    return(0);
  }

  // Parameters of the prediction
  /**************************************************************************
           Normalize size, distance and energy param. for calculations
           Set parameters for gaps
  ****************************************************************************/
  Read_pred_para(Mut_para); //"Mutation_para.in"
  Normalize_parameters();

  /****************************************************
           Find mutated residues
 ******************************************************/
  printf("Examining the effect of mutations: ");
  // Interactions directed along representative atoms
  int Nres=atoms[Ref_kin.atom_num[Ref_kin.N_ref-1]].res+1;
  int i, imut=1;
  char mut_string[1000]="", tmp[200];
  for(i=0; i<Nmut; i++){
    int pos=Posmut[i];
    if((pos<0)||(pos>=Nres)||(Code_AA(AAmut[i])<0)||(Code_AA(AAmut[i])>20)){
      printf("\nWARNIG mutation does not exist (only %d residues)\n", Nres);
      imut=0;
    }
    sprintf(tmp, "%c%d%c", AAwt[i], pos, AAmut[i]); strcat(mut_string, tmp);
    if(i < (Nmut-1))strcat(mut_string, "_");
  }
  printf(" %s\n", mut_string);
  if(imut==0)return(0);

  // Anharmonicity
  float *sigma2; char name_har[80];
  if(anhar==0){sigma2=NM.sigma2; strcpy(name_har, ""); }
  else{sigma2=NM.sigma2_anhar; strcpy(name_har, "_anharmonic");}
  
  /* Three nested sets of atoms:
     - All (compute interactions)
     - kin (Cartesian normal modes) -> kin_atom
     - CC  (Conformation change)    -> cc_kin
     - res (one per residue) -> B_fact
     The Cartesian force is computed for kin atoms and
     projected onto normal modes to obtain the torsional force
  */

  /**************************************************************************
           Prepare reference atoms and determine their displacement str_diff_3
 ****************************************************************************/

  // Match any atom to closest kinetic atom
  int N_kin=Ref_kin.N_ref, N_Cart=3*N_kin, kin_atom[natoms]; 
  Match_atoms(kin_atom, atoms, natoms, Ref_kin.atom_num, N_kin);

  // Structural deformation profile for reference atoms
  // Extract one reference atom per residue
  char SEL[4]="CA"; if(strcmp(REF_CC, "CB")==0)strcpy(SEL, "CB");

  // resatom[i]= Reference atom for residue i
  // Center of mass and inertia tensor
  //double **corr_sum=Allocate_mat2_d(3, 3);
  //double **inertia=Allocate_mat2_d(3, 3);
  int Na=0, j, resatom[Nres]; double m_tot=0;
  float str_diff_3[3*Nres], *str=str_diff_3;
  float str_diff_ave[3]; for(j=0; j<3; j++)str_diff_ave[j]=0;
  //float *xr= x_ref;
  for(i=0; i<Ref_CC.N_ref; i++){
    atom *atm=atoms+Ref_CC.atom_num[i];
    if((strncmp(atm->name, SEL, 2)!=0)&&
       ((seq[atm->res].amm!='G')||(strncmp(atm->name, "CA", 2)!=0)))continue;
    resatom[Na]=Ref_CC.atom_num[i];
    float m=Ref_CC.mass_atom[i];
    float *x=Confchange_Cart+3*i;
    //float *r=atoms[resatom[Na]].r;
    for(j=0; j<3; j++){
      *str=*x; str_diff_ave[j]+=m*(*x); str++; x++; // r++;
    }
    m_tot+=m; Na++;
  }
  if(Na!=Nres)printf("WARNING, %d residues expected, %d found\n", Nres, Na);
  if(Na > Nres){printf("Leaving\n"); return(0);}

  // Subtract center of mass
  for(j=0; j<3; j++)str_diff_ave[j]/=m_tot;
  float str_diff_2[Na]; //, str_diff_abs[Na];
  double RMSD=0;
  str=str_diff_3;
  for(i=0; i<Na; i++){
    double d2=0;
    for(j=0; j<3; j++){*str-=str_diff_ave[j]; d2+=(*str)*(*str); str++;}
    str_diff_2[i]=d2;
    //str_diff_abs[i]=sqrt(d2);
    RMSD+=d2;
  }
  // Predicted fluctuations B1
  double B_sum=0; for(i=0; i<Na; i++)B_sum+=B_pred[i];
  float MSD_nomut=(B_sum/Na);
  //double B_norm=sqrt(B_sum/RMSD);
  RMSD=sqrt(RMSD/Na);
  printf("RMSD(%s atoms)= %.2f\n", SEL, RMSD);

  // Contact list
  int *nc, **clist, **cnum;
  Get_contact_list(&nc, &clist, &cnum, Nres, atoms, Int_list, N_int);

  int mut_neigh[Na],
    n_neigh=Mut_neighbor(mut_neigh, Na, AAwt, Posmut, AAmut, Nmut,
			 clist, cnum, Int_list, atoms, Nres);


  /******************** Compute excess confchange ********************/

  float ave_diff=0;
  for(i=0; i<Na; i++)ave_diff+=str_diff_2[i];
  float ratio=B_sum/ave_diff; //j=0;
  float excess[Na], max_exc_mut=0, max_exc_nomut=0;
  for(i=0; i<Na; i++){
    excess[i]=ratio*str_diff_2[i]/B_pred[i];
    if(mut_neigh[i]){
      if(excess[i]>max_exc_mut)max_exc_mut=excess[i];
    }else if(excess[i]>max_exc_nomut){
      max_exc_nomut=excess[i];
    }
    //for(int k=0; k<3; k++){exc_3[j]=str_diff_3[j]/B1[i]; j++;}
  }

  int N3=3*Na;
  float *B3[3];
  for(j=0; j<3; j++){
    B3[j]=malloc(N3*sizeof(float));
    for(i=0; i<N3; i++)B3[j][i]=0;
  }
  int B_EVEC=0;
  /* B_EVEC=1: amax=3, B3[3][n] B3_ji is the j Cartesian component of the
  // eigenvector of sum_a s^2_a va_ik va_il for the i atom
  // B_EVEC=0: amax=1, B3[1][n] B3_0i = sum_a s_a va_i
  */
  if(B_EVEC){
    amax=3;
    float *B_mat[3]; for(i=0; i<3; i++)B_mat[i]=malloc(3*sizeof(float));
    float **e_vec=Allocate_mat2_f(3, 3), e_val[3];
    for(i=0; i<Na; i++){
      int jj=3*kin_atom[resatom[i]], k;
      for(j=0; j<3; j++)for(k=0; k<=j; k++)B_mat[j][k]=0;
      for(int a=0; a<NM.N; a++){
	if(sigma2[a]<=0)continue;
	float *va=NM.Cart[a]+jj;
	for(j=0; j<3; j++){
	  float sv=sigma2[a]*va[j];
	  for(k=0; k<=j; k++)B_mat[j][k]+=sv*va[k];
	}
      } // end normal modes, symmetrize:
      for(j=0; j<3; j++){
	for(k=j+1; k<3; k++)B_mat[j][k]=B_mat[k][j];
      }
      f_Diagonalize(3, B_mat, e_val, e_vec, 1);
      jj=3*i;
      for(int a=0; a<3; a++){
	float lam=sqrt(e_val[a]), *ev=e_vec[a], *B=B3[a];
	for(j=0; j<3; j++){B[jj]=lam*ev[j]; jj++;}
      }
    } // end i
  }else{
    amax=1;
    for(int a=0; a<NM.N; a++){
      if(sigma2[a]<=0)continue;
      float om_minus_1=sqrt(sigma2[a]), *va=NM.Cart[a]; int k=0;
      for(i=0; i<Na; i++){
	int jj=3*kin_atom[resatom[i]];
	for(j=0; j<3; j++){B3[0][k]+=om_minus_1*va[jj]; jj++; k++;}
      }
    }
  }

  /************************************************************************
           Prepare output
  *************************************************************************/
  // Header
  char header[2000];
  sprintf(header,
	  "# RMSD(%s atoms)= %.2f N= %d residues n= %d kinetic atoms\n",
	  SEL, RMSD, Nres, Ref_kin.N_ref);
  if(Nres!=Na){
    sprintf(tmp, "# gap= %d residues\n", Nres-Na); strcat(header, tmp);
  }
  sprintf(tmp,
	  "# RMSD predicted (no mutation): %.3g "
	  "with force constant kappa= %.3g\n", sqrt(MSD_nomut), KAPPA);
  strcat(header, tmp);

  //Correlation predicted_fluctuations confchange (Cartesian):
  float slope, offset, r_CC_B;
  /*r_CC_B=Corr_coeff(B1, str_diff_abs, Na, &slope, &offset);
  sprintf(tmp, "# Correlation predicted_fluct confchange (abs): %.3f"
  " slope= %.3g offset= %.3g\n", r_CC_B, slope, offset); 
  strcat(header, tmp); */
  r_CC_B=Corr_coeff(B_pred, str_diff_2, Na, &slope, &offset);
  sprintf(tmp, "# Correlation predicted_fluct confchange (^2) : %.3f"
	  " slope= %.3g offset= %.3g\n", r_CC_B, slope, offset);
  strcat(header, tmp);

  //Correlation predicted_fluctuations confchange (torsional):
  float Tors_fluct_obs[N_axes];
  for(i=0; i<N_axes; i++)
    Tors_fluct_obs[i]=Confchange_Tors[i]*Confchange_Tors[i];
  float r=Corr_coeff(Tors_fluct, Tors_fluct_obs, N_axes, &slope, &offset);
  strcat(header, "# Correlation predicted_Tors_fluct Tors_confchange:");
  sprintf(tmp, " %.3f slope=%.2g offset=%.2g\n", r,slope,offset);
  strcat(header, tmp);
  sprintf(tmp, "# Mut: Max.conf_change^2(obs/pred)= %.3g\n"
	  "# No-mut: Max.conf_change^2(obs/pred)= %.3g\n", 
	  max_exc_mut, max_exc_nomut);
  strcat(header, tmp);

  char nameout[200];
  sprintf(nameout, "%s%s_EXPF%.1f_pred_mut.dat",
	  nameout1, name_har, EXP_FORCE); //, OPT_SCORE
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "%s", header);


  fprintf(file_out, "# Input value of force constant: %.4g\n", KAPPA);
  fprintf(file_out, "# Prediction based on %d mutations: %s\n",
	  Nmut, mut_string);
  fprintf(file_out, "# Predicted Mut force based on combination of changes "
	  " of size, stability and optimal distance\n");
  fprintf(file_out, "# Mut.force directed along axes between contacts,"
	  "each weighted with r^-%.1f\n", EXP_FORCE);

  /**************************************************************************
           Compute predicted mutations
 ****************************************************************************/

  // Exponent of the force for computation
  Exp=(1.+EXP_FORCE)/2;
  REXP=pow(RC, EXP_FORCE);
  scale_C=KAPPA*REXP;

  // prediction based on combination of size, stability, distance 
  // Coefficients for the computation of the force
  /*if(C_SIZE==0){C_SIZE=1;}
  if(C_DIST==0){C_DIST=1;}
  if(C_STAB==0){C_STAB=1;}*/
  float C_size=C_SIZE*scale_C, C_stab=C_STAB*scale_C, C_dist=C_DIST*scale_C;

  float Force[N_Cart], Force_coeff[NM.N]; 
  float pred_mut_Cart[N_Cart], pred_mut_Tors[N_axes],
    pred_mut_3[Na], pred_mut_tot_2[Na], pred_mut_2[N3]; //, pred_mut_abs[Na] 
  float sigma2_sc[NM.N], B_pred_sc[Na], *B3_sc[3]; // B1_sc[Na];
  for(i=0; i<NM.N; i++)sigma2_sc[i]=sigma2[i];
  for(i=0; i<Na; i++){
    B_pred_sc[i]=B_pred[i];
    //B1_sc[i]=B1[i];
  }
  for(j=0; j<3; j++){
    B3_sc[j]=malloc(N3*sizeof(float));
    for(i=0; i<N3; i++)B3_sc[j][i]=B3[j][i];
  }
  float Kappa_sc=KAPPA, MSD_nomut_sc=MSD_nomut;

  // Loop on predictions
  int N_force=1; if(OPT_COEFF)N_force++;
  float cc_opt=0; float DE=0;
  for(int iforce=0; iforce<N_force; iforce++){ 
    ini_fit=1; 
    // Compute the force due to the mutation
    if(iforce)Kappa_sc=KAPPA;
    scale_C=Kappa_sc*REXP;

    float cc, MSD_mut, MSD_cross, A_mut=0,A_nomut=0,A_cross=0,offset,scale=1;
    int n_cont_mut;
    if(OPT_COEFF && iforce==1){
      // Optimize coefficients
      cc=Optimize_coeff(&C_size, &C_stab, &C_dist, &MSD_mut, &MSD_cross,
			&A_mut, &A_nomut, &A_cross, &offset,
			pred_mut_3, pred_mut_2, pred_mut_tot_2, Na,
			&n_cont_mut, Force, Force_coeff,
			pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
			AAwt, Posmut, AAmut, Nmut,
			clist, cnum, Int_list, seq, Nres,
			atoms, resatom, Ref_kin, kin_atom,
			NM, sigma2_sc, B3_sc, B_pred_sc, //, B1_sc
			str_diff_3, str_diff_2); //, str_diff_abs
    }else{
      cc=Pred_mut_str(C_size, C_stab, C_dist, &DE, &MSD_mut, &MSD_cross,
		      &A_mut, &A_nomut, &A_cross, &offset,
		      pred_mut_3, pred_mut_2, pred_mut_tot_2, Na,
		      &n_cont_mut, Force, Force_coeff,
		      pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
		      AAwt, Posmut, AAmut, Nmut,
		      clist, cnum, Int_list, seq, Nres,
		      atoms, resatom, Ref_kin, kin_atom,
		      NM, sigma2_sc, B3_sc, B_pred_sc,
		      str_diff_3, str_diff_2); //, pred_mut_abs[Na]
    }
    printf("Mut_pred_score: %.4g  Weights: %.3g %.3g %.3g ",
	   cc, C_size/scale_C, C_stab/scale_C, C_dist/scale_C);
    printf("Predicted RMSD due to mutation= %.3g\n", sqrt(MSD_mut)); 

    if(cc>cc_opt){cc_opt=cc;}

    /***********************************************************
                Print and store the prediction
    *******************************************************/
 
    // Print
    if(iforce==0){
      // Columns printed in the output file for each prediction type
      fprintf(file_out, "# %d contacts of mutated residues:\n#", n_neigh);
      //fprintf(file_out, "# Prediction based on %d contacts: ", n_cont_mut);
      for(i=0; i<n_cont_mut; i++){
	fprintf(file_out, " %c%d-%c%d",
		AA_code[wt_cont[i]], wt_res[i],
		AA_code[mut_cont[i]], mut_res[i]);
	if(i==ncmax)break;
      }
      fprintf(file_out, "\n");
      fprintf(file_out,
	      "# MSD_tot=A_mut*MSD_mut+A_nomut*MSD_nomut+A_cross*MSD_cross\n");
      fprintf(file_out, "#1=MSD 2=MSD_pred(tot) "
	      "3=Correlation(confchange^2,pred_tot^2) "
	      "4=MSD_pred(mut) 5=MSD_pred(no mut) 6=MSD_cross 7=offset "
	      "8=A_mut 9=A_nomut 10=A_cross ");
	      //"11=Cosine(confchange,pred_mut) 12=scale_cosine");
      //" 10=Corr(pred_mut,confchange)_tors 11=slope 12=offset\n"); 
      fprintf(file_out, "\n#\n");
    }

    for(int iter=0; iter<5; iter++){
      if(iter){ // rescale mutation parameters and force constant
	if(iforce==0)continue; 
	printf("Rescaling weights and Kappa, iter= %d\n", iter);
	ini_fit=1;
	cc=Pred_mut_str(C_size, C_stab, C_dist, &DE, &MSD_mut, &MSD_cross,
			&A_mut, &A_nomut, &A_cross, &offset,
			pred_mut_3, pred_mut_2, pred_mut_tot_2, Na,
			&n_cont_mut, Force, Force_coeff,
			pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
			AAwt, Posmut, AAmut, Nmut,
			clist, cnum, Int_list, seq, Nres,
			atoms, resatom, Ref_kin, kin_atom,
			NM, sigma2_sc, B3_sc, B_pred_sc,
			str_diff_3, str_diff_2); //, str_diff_abs
	printf("Mut_pred_score: %.4g  Weights: %.3g %.3g %.3g ",
	       cc, C_size/scale_C, C_stab/scale_C, C_dist/scale_C);
	printf("Predicted RMSD due to mutation= %.3g\n", sqrt(MSD_mut)); 
      }

      scale_C=REXP*Kappa_sc;
      fprintf(file_out,"# Force constant: %.4g ", Kappa_sc);
      fprintf(file_out," C_SIZE= %.3g C_STAB= %.3g C_DIST= %.3g ",
	      C_size/scale_C, C_stab/scale_C, C_dist/scale_C);
      if(iforce==1 && OPT_COEFF){fprintf(file_out," optimized");}
      else{fprintf(file_out, " input");}
      fprintf(file_out, " weights");
      if(iter)fprintf(file_out, ", rescaled by %.3g", scale);
      fprintf(file_out, "\n");

      fprintf(file_out, "%.3f", RMSD*RMSD);
      float MSD_all=A_mut*MSD_mut+A_nomut*MSD_nomut_sc+A_cross*MSD_cross;
      fprintf(file_out, "\t%.3f", MSD_all);
      fprintf(file_out, "\t%.4g", cc);
      fprintf(file_out, "\t%.3f", MSD_mut);
      fprintf(file_out, "\t%.3f", MSD_nomut_sc);
      fprintf(file_out, "\t%.3f", MSD_cross);
      fprintf(file_out, "\t%.3g\t%.3g\t%.3g\t%.3g",
	      offset, A_mut, A_nomut, A_cross);
      //float c_scale, cos=Cosine(&c_scale, pred_mut_3, str_diff_3, N3);
      //fprintf(file_out, "\t%.3g\t%.3g", cos, c_scale);
      //Correlation predicted_mutation confchange (torsional):
      //for(i=0; i<N_axes; i++)pred_mut_Tors[i]*=pred_mut_Tors[i];
      //r=Corr_coeff(pred_mut_Tors, Confchange_Tors, N_axes, &slope, &offset);
      //fprintf(file_out, "\t%.3f\t%.2g\t%.2g", r, slope, offset);
      fprintf(file_out, "\n");

      float scale_K2=1;
      if(iforce){
	// Rescale parameters for next iteration
	if(A_nomut>0)scale_K2=A_nomut;
	scale=1./scale_K2;
	if(A_mut>0)scale*=sqrt(A_mut);
	C_size*=scale; C_stab*=scale; C_dist*=scale;
	for(i=0; i<NM.N; i++)sigma2_sc[i]*=scale_K2;
	float scale_K1=sqrt(scale_K2);
	for(i=0; i<Na; i++){
	  B_pred_sc[i]*=scale_K2;
	  for(j=0; j<3; j++)B3_sc[j][i]*=scale_K1;
	}
	Kappa_sc/=scale_K2;
	MSD_nomut_sc*=scale_K2;
	fprintf(file_out,"# Next force constant: %.4g ", Kappa_sc);
	if(A_nomut<=0){ // Do not rescale
	  fprintf(file_out,"(Not rescaled, A_nomut=%.2g<=0)\n", A_nomut);
	}else{
	  fprintf(file_out,"(Kappa/A_nomut=%.3g)\n", A_nomut);
	}
	fprintf(file_out,"# Next parameters: "
		"C_SIZE= %.3g C_STAB= %.3g C_DIST= %.3g ",
		C_size/scale_C, C_stab/scale_C, C_dist/scale_C);
	if(scale==1){ // Do not rescale
	  fprintf(file_out,"(Not rescaled, A=%.2g)\n#\n", A_nomut);
	}else{
	  fprintf(file_out,"(C*scale=%.3g)\n#\n", scale);
	}
      }
      if(fabs(scale-1)<0.05 && fabs(scale_K2-1)<0.05)break;

    } // end iter

    if(0){
      // Print contribution of normal modes to the prediction
      int N_print=30, a; char nameout[200];
      sprintf(nameout, "%s%s_mutation_%d.dat", nameout1, name_har, iforce);
      FILE *file_out=fopen(nameout, "w");
      printf("Writing %s\n", nameout);
      float kappa=Collectivity_norm2(Force, N_Cart);
      fprintf(file_out, "# Collectivity of force: %.2f\n", kappa/3);
      fprintf(file_out,"#mod sigma2 contr2cc contr2force collectivity\n");
      double sum=0;
      for(a=0; a<NM.N_relevant; a++){
	Force_coeff[a]*=Force_coeff[a]; sum+=Force_coeff[a];
      }
      for(a=0; a<N_print; a++){
	fprintf(file_out, "%d\t%.4f\t%.4f\t%.4f\t%.3f\n", //
		a, sigma2[a], NM.confchange2[a],
		Force_coeff[a]/sum, NM.Cart_coll[a]);
      }
      fclose(file_out);
    }
    
  } // end iforce
  printf("Writing %s\n", nameout);  
  fclose(file_out);
  
  // Output predicted deformation of reference atoms mut_CC
  float *m_CC=mut_CC_out;
  for(i=0; i<Ref_CC.N_ref; i++){
    float *x=pred_mut_Cart+3*kin_atom[Ref_CC.atom_num[i]];
    for(int j=0; j<3; j++){*m_CC=*x; x++; m_CC++;}
  }
  for(i=0; i<N_axes; i++)mut_Tors_out[i]=pred_mut_Tors[i];
  
  
  if(0){
    // Print on screen (disabled)
    printf("Force: ");
    for(i=0; i<N_Cart; i++)printf(" %.2f", Force[i]);
    printf("\n");
    printf("mut_Cart: ");
    for(i=0; i<N_Cart; i++)printf(" %.2f",pred_mut_Cart[i]);
    printf("\n");
    printf("mut_CC: ");
    for(i=0; i<Ref_CC.N_Cart; i++)printf(" %.2f", mut_CC_out[i]);
    printf("\n");
    printf("mut_Tors: ");
    for(i=0; i<N_axes; i++)printf(" %.2f", mut_Tors_out[i]);
    printf("\n");
  }

  /******************************************************************/
  if(1){
    // Write excess_confchange and predicted mutations for every residue
    sprintf(nameout, "%s%s_EXPF%.1f_mutation.dat", nameout1,name_har,EXP_FORCE);
    file_out=fopen(nameout, "w");
    fprintf(file_out, "%s", header);
    fprintf(file_out,"#1=conf_ch^2 2=pred_tot^2 3=pred_mut^2 4=mut_neighbor");
    //int k; for(k=1; k<=4; k++)fprintf(file_out, " %d=pred_mut_%d", k+4, k); 
    fprintf(file_out, " pos aa\n");
    for(i=0; i<Na; i++){
      fprintf(file_out, "%.3g\t%.3g\t%.4g\t%d",
	      str_diff_2[i], pred_mut_tot_2[i], pred_mut_2[i], mut_neigh[i]);
      //for(k=0; k<N_TYPES; k++)
      //fprintf(file_out, "\t%.3g", pred_mut_abs_f[k][i]);
      int a=atoms[resatom[i]].res;
      fprintf(file_out, "\t%s\t%c\n", seq[a].pdbres, seq[a].amm);
    }
    fclose(file_out);
    printf("Writing %s\n", nameout);
  }

  /******************************************************************/
  // Clean memory
  Empty_matrix_i(clist, Nres);
  Empty_matrix_i(cnum, Nres); free(nc);
  for(j=0; j<3; j++){
    free(B3[j]); free(B3_sc[j]);
  }
  return(1);
}

float Pred_mut_str(float C_size, float C_stab, float C_dist,
		   float *DE, float *MSD_mut, float *MSD_cross,
		   float *A_mut, float *A_nomut, float *A_cross, float *offset,
		   float *pred_mut_3, float *pred_mut_2, float *pred_mut_tot_2,
		   int Na, int *n_cont_mut, float *Force, float *Force_coeff,
		   float *pred_mut_Cart, int N_Cart, 
		   float *pred_mut_Tors, int N_axes,
		   char *AAwt, int *Posmut, char *AAmut, int Nmut,
		   int **clist, int **cnum, struct interaction *Int_list,
		   struct residue *seq, int Nres, atom *atoms, int *resatom,
		   struct Reference Ref_kin,
		   int *kin_atom, struct Normal_Mode NM,
		   float *sigma2, float **B3, float *B_pred, 
		   float *str_diff_3, float *str_diff_2)
{
  // Compute force
  *n_cont_mut=
    Mutation_force(Force, N_Cart, C_size, C_stab, C_dist, 
		   AAwt, Posmut, AAmut, Nmut, clist, cnum, Int_list,
		   seq, Nres, atoms, Ref_kin, kin_atom);

  // Compute deformation mut_Cart and mut_Tors
  Compute_mut_def(pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes,
		  DE, Force_coeff, Force, NM, sigma2, Ref_kin.mass_sum);
      
  // Deformation of resatom mut_3, mut_abs2 and MSD_mut
  float pred_mut_ave[3], *str=pred_mut_3; int i, j;
  for(j=0; j<3; j++)pred_mut_ave[j]=0;
  for(i=0; i<Na; i++){
    float *x=pred_mut_Cart+3*kin_atom[resatom[i]];
    for(j=0; j<3; j++){*str=*x; pred_mut_ave[j]+=*x; x++; str++;}
  }
  // Subtract center of mass
  for(j=0; j<3; j++)pred_mut_ave[j]/=Na;
  *MSD_mut=0; str=pred_mut_3;
  for(i=0; i<Na; i++){
    double d2=0;
    for(j=0; j<3; j++){*str-=pred_mut_ave[j]; d2+=(*str)*(*str); str++;}
    *MSD_mut+=d2;
    pred_mut_2[i]=d2;
  }
  *MSD_mut/=Na;

  //for(i=0; i<20; i++){printf("%.3g ",pred_mut_2[i]);} printf("\n");
  //printf("Predicted MSD: %.4g\n", *MSD_mut);

  *MSD_cross=0;
  // No mutant structure is present for computing observed differences
  if((str_diff_3==NULL)||(str_diff_2==NULL))return(0);


  // pred_mut_cross[i]= 2 pred_mut[i] X B3[i]
  // B3_0i = sum_a s_a va_i (B_EVEC==0)
  // Does it make sense? va_i and -va_i are equivalent,
  // the dot product should be zero... 
  float pred_mut_cross[Na];
  for(i=0; i<Na; i++){pred_mut_cross[i]=0;}
  for(int a=0; a<amax; a++){ // amax=1 (B_EVEC==0) or 3 (B_EVEC==1)
    float *bb=B3[a]; str=pred_mut_3;
    // Scalar product pred_mut_3*B3
    for(i=0; i<Na; i++){
      double d_cross=0;
      for(j=0; j<3; j++){d_cross+=(*str)*(*bb); str++; bb++;}
      pred_mut_cross[i]+=d_cross*d_cross;
    }
  }
  (*MSD_cross)=0;
  for(i=0; i<Na; i++){
    pred_mut_cross[i]=4*sqrt(pred_mut_cross[i]); // Note the factor 4!
    (*MSD_cross)+=pred_mut_cross[i];
  }
  *MSD_cross/=Na;

  // Fit observed deformation profile squared str_diff_2[i] 
  float score;
  int Npar=3;
  float D_out[Npar], Y_pred[Na], *X[Na];
  for(i=0; i<Na; i++){
    X[i]=malloc(Npar*sizeof(float));
    X[i][0]=pred_mut_2[i];
    X[i][1]=B_pred[i];
    X[i][2]=pred_mut_cross[i];
  }

  struct ridge_fit fit; fit.A=malloc(Npar*sizeof(float)); char type;
  if(ini_fit){type='C';} // Specific-heat fit
  else{type='L'; fit.Lambda=Lambda_RRR;}
  score=Ridge_regression(&fit, Y_pred, D_out, "ridge_regression",
			 X, str_diff_2, Na, Npar, type);
  printf("RRR3: A_mut= %.3g A_nomut= %.3g A_cross= %.3g\n",
	 fit.A[0], fit.A[1], fit.A[2]);

  if(fit.A[0]<=0 || fit.A[1]<=0){ //*A_mut<=0 *A_nomut<0
    // Fit without A_cross
    fit.A[2]=0; Npar=2;
    score=Ridge_regression(&fit, Y_pred, D_out, "ridge_regression",
			   X, str_diff_2, Na, Npar, type);
    printf("RRR2: A_mut= %.3g A_nomut= %.3g A_cross= %.3g\n",
	   fit.A[0], fit.A[1], fit.A[2]);
  }

  Lambda_RRR=fit.Lambda;
  *A_mut=fit.A[0];
  *A_nomut=fit.A[1];
  *A_cross=fit.A[2];
  *offset=0;
  for(i=0; i<Na; i++){
    (*offset)+=str_diff_2[i]-Y_pred[i];
    pred_mut_tot_2[i]=Y_pred[i];
    free(X[i]);
  }
  (*offset)/=Na;
  free(fit.A);

  ini_fit=0;

  if(*A_nomut<=0){
    // Fit only with the predicted effect of mutation
    score=Fit_1(str_diff_2, pred_mut_2, Na, A_mut,offset);
    *A_nomut=0; *A_cross=0;
    printf("Fit1: A_mut= %.3g A_nomut= %.3g A_cross= %.3g\n",
	   *A_mut, *A_nomut, *A_cross);
  }else if(*A_mut<=0){
    // Fit only with thermal fluctuations
    score=Fit_1(str_diff_2, B_pred, Na, A_nomut,offset);
    *A_mut=0; *A_cross=0;
    printf("Fit1: A_mut= %.3g A_nomut= %.3g A_cross= %.3g\n",
	   *A_mut, *A_nomut, *A_cross);
  }

  return(score);
}


int Mutation_force(float *Force, int N_Cart,
		   float C_size, float C_stab, float C_dist, 
		   char *AAwt, int *Posmut, char *AAmut, int Nmut,
		   int **clist, int **cnum, struct interaction *Int_list,
		   struct residue *seq, int nres, atom *atoms,
		   struct Reference Ref_kin, int *kin_atom)
{

  int am=-1, ai, j, nc=0, iter;
  for(j=0; j<N_Cart; j++)Force[j]=0;
  for(int kmut=0; kmut<Nmut; kmut++){
    int pmut=Posmut[kmut];
    int aaw=Code_AA(AAwt[kmut]);   // wild-type aa
    int aam=Code_AA(AAmut[kmut]);  // mutated aa
    double dsize=(natm_norm[aam]-natm_norm[aaw]);
    if((aaw==20)||(aam==20)){iter=2;} // gap!
    else{iter=1;}
    /* Since atoms in the gap are removed, we have to analyze
       the flanking residues */
    int p1=-1;
    for(int it=0; it<iter; it++){
      if(iter==2){
	if(it==0){
	  while((pmut>=0)&&(clist[pmut][0]<=0))pmut--;
	  if(pmut<0)continue;
	  p1=pmut;
	}else if(it==1){
	  pmut=Posmut[kmut]+1; while((pmut<nres)&&(clist[pmut][0]<=0))pmut++;
	  if(pmut>=nres)continue;
	}
      }
      int *ic=clist[pmut]; // residue in contact with pmut
      int *ii=cnum[pmut];  // interaction label of the contact
      while(*ic >= 0){
	//if(abs(pmut-(*ic))<3)goto next; // short range contact discarded
	int i1=Int_list[*ii].i1, i2=Int_list[*ii].i2;
	if(atoms[i1].res==pmut){
	  am=i1; ai=i2;
	}else if(atoms[i2].res==pmut){
	  am=i2; ai=i1;
	}else{
	  printf("ERROR in mutation, interaction %d res %d - %d pmut= %d\n",
		 *ii, atoms[i1].res, atoms[i2].res, pmut); exit(8);
	}
	int aai=Code_AA(seq[*ic].amm);
	ic++; ii++;
	
	int res_i=atoms[ai].res;
	//if((abs(res_i-pmut)<3))goto next; //&&(FORCE_TYPE==1)

	// contacts of main chain atoms of mutant residue are discarded (?)
	if(PRINT_FLANK &&(res_i==p1)){
	  printf("Gap flanking residues %d and %d\n", p1, pmut);
	  PRINT_FLANK=0;
	  //}else if(Mainchain(atoms[am].name)){
	  //continue; //&& FORCE_TYPE!=3
	}

	float f =C_size*dsize;
	f+=C_dist*(D_Res[aam][aai]-D_Res[aaw][aai]);
	float dE=(Econt[aam][aai]-Econt[aaw][aai]);
	if(((aaw==20)&&(dE<0))|| // Insertion, attracting
	   ((aam==20)&&(dE>0)))continue; // Deletion, repelling
	f+=C_stab*dE;

	// Record contacts considered for force computation
	if(nc<ncmax){
	  mut_cont[nc]=aai; mut_res[nc]=res_i;
	  wt_cont[nc]=aaw; wt_res[nc]=pmut;
	}
	nc++;


	double r[3], r2=0, *rm=atoms[am].r, *ri=atoms[ai].r;
	for(j=0; j<3; j++){
	  r[j]=rm[j]-ri[j]; r2+=r[j]*r[j];
	}
	f/=pow(r2, Exp);

	int km=3*kin_atom[am], ki=3*kin_atom[ai]; 
	for(j=0; j<3; j++){
	  float fj=f*r[j];
	  Force[km+j]+=fj; //?
	  Force[ki+j]-=fj;
	}
      } // end contact
    } // end iter (in case of gap)
  }
  //printf("Force computation based on %d contacts\n", nc);
  return(nc);
}


void Compute_mut_def(float *pred_mut_Cart, int N_Cart,
		     float *pred_mut_Tors, int N_axes,
		     float *DE, float *Force_coeff, float *Force,
		     struct Normal_Mode NM, float *sigma2, float mass)
{
  int a, i;
  for(i=0; i<N_Cart; i++)pred_mut_Cart[i]=0;
  for(i=0; i<N_axes; i++)pred_mut_Tors[i]=0;

  /*double F=0;
  for(i=0; i<N_Cart; i++){F+=Force[i]*Force[i];}
  printf("RMS of mutational force: %.3f\n", sqrt(3*F/N_Cart)); */
  double DE_tmp=0;

  // Compute predicted deformation of kinetic atoms (Cart)
  for(a=0; a<NM.N; a++){
    if((NM.Cart_coll[a]<COLL_THR_MUT)||(sigma2[a]<=0)){
      Force_coeff[a]=0; continue;
    }
    float *f=Force, *x=NM.Cart[a];
    double xf=0; // Projection of force on normal mode
    for(i=0; i<N_Cart; i++){xf+=(*x)*(*f); x++; f++;}
    Force_coeff[a]=xf;
    //(*DE)+=xf*xf/sigma2[a];
    DE_tmp+=xf*xf*sigma2[a]; // OK
    xf*=sigma2[a];
    float *m=pred_mut_Cart; x=NM.Cart[a];
    for(i=0; i<N_Cart; i++){
      (*m)+=(*x)*xf; x++; m++;
    }
    m=pred_mut_Tors; x=NM.Tors[a];
    for(i=0; i<N_axes; i++){
      (*m)+=(*x)*xf; x++; m++;
    }
  }
  *DE=DE_tmp;
}

int Mainchain(char *atom_name){
  if(strncmp(atom_name,"N ",2)==0)return(1);
  if(strncmp(atom_name,"CA",2)==0)return(1);
  if(strncmp(atom_name,"C ",2)==0)return(1);
  if(strncmp(atom_name,"O ",2)==0)return(1);
  return(0);
}

float *Convert_to_tors(float *Cart, int N3, float **J_ar, int N_axes)
{
  // Torsional force = J^t f = sum_i dr_i/dphi_a f_i
  float *Tors=malloc(N_axes*sizeof(float)); int a, i;
  for(a=0; a<N_axes; a++){
    double sum=0; float *J=J_ar[a], *f=Cart;
    for(i=0; i<N3; i++){sum+=(*J)*(*f); J++; f++;}
    Tors[a]=sum;
  }
  return(Tors);
}

void  Match_atoms(int *ref_atom, atom *atoms, int natoms,
		  int *atom_ref, int N_ref)
{
  /* For every atom i=1,natoms ref_atom[i] is the index of the reference
     atom that is closest to i */
  int i; for(i=0; i<natoms; i++)ref_atom[i]=-1;
  for(i=0; i<N_ref; i++)ref_atom[atom_ref[i]]=i;
  int Nr1=N_ref-1, ini_r=0, n=0;
  for(i=0; i<natoms; i++){
    if(ref_atom[i]>=0)continue;
    n++;
    atom *atom=atoms+i, *atom2;
    while((ini_r<Nr1)&&(atoms[atom_ref[ini_r]].res<atom->res))ini_r++;
    if(atoms[atom_ref[ini_r]].res != atom->res){
      printf("ERROR unmatch_atom in residue %c%d\n", atom->aa, atom->res);
      exit(8);
    }
    float d2_min=1000; int ir=ini_r, ir_min=-1;
    while(ir<N_ref){
      atom2=atoms+atom_ref[ir];
      if(atom2->res==atom->res){
	float d2=Distance_square(atom->r, atom2->r);
	if((d2<d2_min)||(ir_min<0)){d2_min=d2; ir_min=ir;}
      }
      ir++;
    }
    if(ir_min<0){
      printf("ERROR, reference atom %d in atom %d residue %c%d not valid\n",
	     ir_min, i, atom->aa, atom->res);
    }
    ref_atom[i]=ir_min;
  }
  if((N_ref+n)!=natoms){
    printf("ERROR, N_ref+%d=%d atoms matched over %d\n",
	   n, N_ref+n, natoms); exit(8);
  }
  if(0){
    for(i=0; i<natoms; i++){
      atom *atom=atoms+i;
      printf("%d %s(%c%d) ", i, atom->name, atom->aa, atom->res);
      int j=atom_ref[ref_atom[i]]; atom=atoms+j;
      printf(" ref: %d  %s(%c%d)\n", j, atom->name, atom->aa, atom->res);
    }
    exit(8);
  }
}

float Cosine(float *scale, float *xx, float *yy, int n)
{
  double sum_xy=0, sum_x2=0, sum_y2=0;
  float *x=xx, *y=yy;
  for(int i=0; i<n; i++){
    sum_xy+=(*x)*(*y);
    sum_x2+=(*x)*(*x);
    sum_y2+=(*y)*(*y);
    x++; y++;
  }
  if(sum_x2){
    *scale=sqrt(sum_y2/sum_x2);
    return(sum_xy/sqrt(sum_x2*sum_y2));
  }else{
    *scale=0;
    return(0);
  }
}

int Mut_neighbor(int *mut_neigh, int Na,
		 char *AAwt, int *Posmut, char *AAmut, int Nmut,
		 int **clist, int **cnum, struct interaction *Int_list,
		 atom *atoms, int nres)
{
  int i, iter;
  for(i=0; i<Na; i++)mut_neigh[i]=0;
  for(int kmut=0; kmut<Nmut; kmut++){
    int pmut=Posmut[kmut];
    int aaw=Code_AA(AAwt[kmut]);   // wild-type aa
    int aam=Code_AA(AAmut[kmut]);  // mutated aa
    if((aaw==20)||(aam==20)){pmut--; iter=2;} //
    else{iter=1;}
    for(int it=0; it<iter; it++){
      if(pmut<0)continue;
      if(it==1){pmut=Posmut[kmut]+1; if(pmut>=nres)continue;}
      int *ii=cnum[pmut];  // interaction label of the contact
      int *ic=clist[pmut]; // residue in contact with pmut
      while(*ic >= 0){
	//if(abs(pmut-(*ic))<3)goto next; // short range contact discarded
	int res1=atoms[Int_list[*ii].i1].res, res2=atoms[Int_list[*ii].i2].res;
	if((res1 !=pmut)&&(res2!=pmut)){
	  printf("ERROR in mutation, interaction %d res %d - %d pmut= %d\n",
		 *ii, res1, res2, pmut); exit(8);
	}
	mut_neigh[res1]=1;
	mut_neigh[res2]=1;
	ic++; ii++;
      }
    }
  }
  int n=0; for(i=0; i<Na; i++)if(mut_neigh[i])n++;
  return(n-1);
}

void Read_pred_para(char *name_in){
  FILE *file_in=fopen(name_in, "r");
  if(file_in==NULL){
    printf("WARNING, parameter file %s not found\n", name_in);
    return;
  }
  char string[100];
  printf("Reading weights of mutation model in %s: ", name_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "EXP_FORCE", 9)==0){
      sscanf(string+10, "%f", &EXP_FORCE);
    }else if(strncmp(string, "OPT_COEFF", 9)==0){
      sscanf(string+10, "%d", &OPT_COEFF);
    }else if(strncmp(string, "OPT_SCORE", 9)==0){
      char c;
      sscanf(string+10, "%c", &c);
      if((c!='C')&&(c!='A')){
	printf("WARNING, %c is not an allowed score type\n", c);
	printf("Allowed types are C (Cartesian) and A (Absolute). ");
	printf("Using default %c\n", OPT_SCORE);
      }else{
	OPT_SCORE=c;
      }
    }else if(strncmp(string, "C_SIZE", 6)==0){
      sscanf(string+7, "%f", &C_SIZE);
    }else if(strncmp(string, "C_DIST", 6)==0){
      sscanf(string+7, "%f", &C_DIST);
    }else if(strncmp(string, "C_STAB", 6)==0){
      sscanf(string+7, "%f", &C_STAB);
    }else{
      printf("WARNING, unknown record %s\n", string);
    }
  }
  printf("size,dist,stab= %.3g %.3g %.3g\n", C_SIZE, C_DIST, C_STAB);
  fclose(file_in);
}

float Optimize_coeff(float *C_size_opt, float *C_stab_opt, float *C_dist_opt,
		     float *MSD_mut, float *MSD_cross,
		     float *A_mut, float *A_nomut,float *A_cross,float *offset,
		     float *pred_mut_3, float *pred_mut_2,
		     float *pred_mut_tot_2, 
		     int Na, int *n_cont_mut, float *Force, float *Force_coeff,
		     float *pred_mut_Cart, int N_Cart, 
		     float *pred_mut_Tors, int N_axes,
		     char *AAwt, int *Posmut, char *AAmut, int Nmut,
		     int **clist, int **cnum, struct interaction *Int_list,
		     struct residue *seq, int Nres, atom *atoms, int *resatom,
		     struct Reference Ref_kin,
		     int *kin_atom, struct Normal_Mode NM, float *sigma2,
		     float **B3, float *B_pred,
		     float *str_diff_3, float *str_diff_2)
{
  float y_max=-1000; float DE=0;
  ini_fit=1;

  // Initialize parameters
  float C_size=*C_size_opt, C_stab=*C_stab_opt, C_dist=*C_dist_opt;

  printf("Optimizing mutation coefficients\n");
  int ifit=0, j, n_fail=0;
  for(int it=0; it<200; it++){
    float x[3], y[3], yy, *C_move, *C_move_opt;
    if(ifit==0){C_move=&C_size; C_move_opt=C_size_opt; ifit=1;}
    else{C_move=&C_dist; C_move_opt=C_dist_opt; ifit=0;}
    if(it==0){
      y_max=Pred_mut_str(C_size, C_stab, C_dist, &DE, MSD_mut, MSD_cross,
			 A_mut, A_nomut, A_cross, offset,
			 pred_mut_3, pred_mut_2, pred_mut_tot_2, Na,
			 n_cont_mut, Force, Force_coeff,
			 pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
			 AAwt,Posmut,AAmut,Nmut,clist,cnum,Int_list,seq,Nres,
			 atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
			 B3, B_pred, str_diff_2, str_diff_3);
      printf("%d cc: %.4g Weights: %.3g %.3g %.3g\n",
	     -1,y_max,C_size/scale_C,C_stab/scale_C,C_dist/scale_C);
    }
    x[1]=*C_move; y[1]=y_max; int j_max=1;
    float x_min=0, x_max=x[1]*100;
    for(j=0; j<=2; j+=2){
      if(j==0){x[0]=x[1]*0.85;}
      else if(j==2){x[2]=x[1]*1.15;}
      *C_move=x[j];
      y[j]=Pred_mut_str(C_size, C_stab, C_dist, &DE, MSD_mut, MSD_cross,
			A_mut, A_nomut, A_cross, offset,
			pred_mut_3, pred_mut_2, pred_mut_tot_2, Na,
			n_cont_mut, Force, Force_coeff,
			pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
			AAwt,Posmut,AAmut,Nmut,clist,cnum,Int_list,seq,Nres,
			atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
			B3, B_pred, str_diff_3, str_diff_2);
      if(y[j] > y_max){y_max=y[j]; j_max=j;}
    }
    // Determine new value of x and test its score
    *C_move=Find_max_quad(x[0], x[1], x[2], y[0], y[1], y[2], x_min, x_max);
    if(*C_move==x[j_max]){yy=y_max;}
    else{
      yy=Pred_mut_str(C_size, C_stab, C_dist, &DE, MSD_mut, MSD_cross, 
		      A_mut, A_nomut, A_cross, offset,
		      pred_mut_3, pred_mut_2, pred_mut_tot_2, Na,
		      n_cont_mut, Force, Force_coeff,
		      pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
		      AAwt,Posmut,AAmut,Nmut,clist,cnum,Int_list,seq,Nres,
		      atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
		      B3, B_pred, str_diff_3, str_diff_2);
    }
    printf("%d cc: %.4g Weights: %.3g %.3g %.3g\n",
	   it,yy,C_size/scale_C,C_stab/scale_C,C_dist/scale_C);
    if(yy > y_max){
      n_fail=0; y_max=yy;
    }else if(j_max!=1){
      n_fail=0; *C_move=x[j_max];
    }else{
      // Score does not improve. If two consecutive failures, exit.
      n_fail++; *C_move=x[j_max];
    }
    *C_move_opt=*C_move;
    if(n_fail>=2){
      printf("%d consecutive decreases, exiting\n", n_fail); break;
    }
  } // end iterations
  // Recompute
  y_max=Pred_mut_str(*C_size_opt, *C_stab_opt, *C_dist_opt, &DE,
		     MSD_mut, MSD_cross,
		     A_mut, A_nomut, A_cross, offset,
		     pred_mut_3, pred_mut_2, pred_mut_tot_2, Na,
		     n_cont_mut, Force, Force_coeff,
		     pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
		     AAwt,Posmut,AAmut,Nmut,clist,cnum,Int_list,seq,Nres,
		     atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
		     B3, B_pred, str_diff_3, str_diff_2);
  printf("mut_pred_score: %.4g  Weights: %.3g %.3g %.3g\n",
  	 y_max, *C_size_opt/scale_C, *C_stab_opt/scale_C, *C_dist_opt/scale_C);

  return(y_max);
}


void Normalize_parameters()
{
  // a.a. 20 = gap
  double E1=0, E2=0, D1=0, D2=0, na1=0, na2=0;
  int a, b, n=21;
  for(a=0; a<=20; a++){
    na1+=atom_numb[a];
    na2+=atom_numb[a]*atom_numb[a];
  }
  na1/=n; na2=sqrt(na2/n-na1*na1);
  for(a=0; a<=20; a++){
    natm_norm[a]=(atom_numb[a]-na1)/na2;
  }
  // E(gap)=0; D(gap)=min D_Res
  n=230;
  double D_min=1000;
  for(a=0; a<20; a++){
    for(b=0; b<=a; b++){
      E1+=Econt[a][b];
      E2+=Econt[a][b]*Econt[a][b];
      D1+=D_Res[a][b];
      D2+=D_Res[a][b]*D_Res[a][b];
      if(D_Res[a][b]<D_min)D_min=D_Res[a][b];
    }
  }
  E1/=n; E2=sqrt(E2/n-E1*E1);
  D1+=D_min*20; D2+=D_min*D_min*20;
  D1/=n; D2=sqrt(D2/n-D1*D1);
  float E_gap=-E1/E2;
  float D_gap=(D_min-D1)/D2;
  for(a=0; a<20; a++){
    for(b=0; b<20; b++){
      Econt_norm[a][b]=(Econt[a][b]-E1)/E2;
      D_Res_norm[a][b]=(D_Res[a][b]-D1)/D2;
    }
    Econt_norm[a][20]=E_gap;
    Econt_norm[20][a]=E_gap;
    D_Res_norm[a][20]=D_gap;
    D_Res_norm[20][a]=D_gap;
  }
  printf("Parameters of mutation model type, ave, s.d.:\n");
  printf("SIZE natom: %.3g %.3g\n\n", na1, na2);
  printf("DIST D_Res: %.3g %.3g\n", D1, D2);
  printf("STAB Econt: %.3g %.3g\n", E1, E2);
}

void Predict_mutations(struct Normal_Mode NM, float KAPPA,
		       atom *atoms, int natoms,
		       int N_axes, struct Reference Ref_kin, 
		       struct interaction *Int_list, int N_int,
		       struct residue *seq, int Nres,
		       char *nameout1, char *Mut_para,
		       int PRED_MUT, int IWT)
{
  int i_rank[20]; // For RMSD[a_wt] //IWT=3, 
  int ALLPAIR=0; if(PRED_MUT==2)ALLPAIR=1;

  // Exponent of the force for computation
  Read_pred_para(Mut_para);
  Exp=(1.+EXP_FORCE)/2;
  REXP=pow(RC, EXP_FORCE); scale_C=KAPPA*REXP;
  float C_size=C_SIZE*scale_C, C_stab=C_STAB*scale_C, C_dist=C_DIST*scale_C;
  Normalize_parameters();

  /* For all positions in the protein and all possible amino acid changes,
     predict the structural effect of the mutation */

  float *sigma2=NM.sigma2;
  /**************************************************************************
           Prepare reference atoms and determine their displacement str_diff_3
  ****************************************************************************/

  // Match any atom to closest kinetic atom
  int N_kin=Ref_kin.N_ref, N_Cart=3*N_kin, kin_atom[natoms]; 
  Match_atoms(kin_atom, atoms, natoms, Ref_kin.atom_num, N_kin);

  // Extract one reference atom per residue
  int resatom[Nres];
  int Na=Extract_reference_atoms(resatom, atoms, natoms, seq, Nres);

  // Contact list
  int *nc, **clist, **cnum;
  Get_contact_list(&nc, &clist, &cnum, Nres, atoms, Int_list, N_int);

  // Output files
  char para[200], head[800],tmp[50];
  sprintf(para, "# Mutation par.: C_SIZE= %.3g C_STAB= %.3g C_DIST=%.3g IWT=%d",
	  C_SIZE, C_STAB, C_DIST, IWT);
  strcpy(head, "#Mut\tncont");
  for(int a=0; a<20; a++){
    sprintf(tmp, "\t%c", AA_code[a]); strcat(head, tmp);
  }
  strcat(head,"\n");

  char name[200]; sprintf(name,"%s.mut_RMSD.dat", nameout1);
  char what[200]="Predicted RMSD of all mutations from WT";
  FILE *file_rmsd=fopen(name, "w");
  fprintf(file_rmsd, "# %s\n%s\n%s", what, para, head);
  printf("Writing %s in %s\n", what,name);

  sprintf(name,"%s.mut_DE.dat", nameout1);
  strcpy(what, "Predicted DE of all mutations from WT");
  FILE *file_de=fopen(name, "w");
  fprintf(file_de, "# %s\n%s\n%s", what,para,head);
  printf("Writing %s in %s\n",what, name);

  FILE *file_de_all=NULL, *file_rmsd_all=NULL;
  if(ALLPAIR){
    strcpy(head, "#pos");
    for(int a=0; a<20; a++){
      for(int b=a+1; b<20; b++){
	sprintf(tmp, "\t%c%c", AA_code[a], AA_code[b]);
	strcat(head, tmp);
      }
    }
    strcat(head,"\n");

    strcpy(what, "Predicted RMSD of all pairs of mutations");
    sprintf(name,"%s.mut_RMSD_all.dat", nameout1);
    file_rmsd_all=fopen(name, "w");
    fprintf(file_rmsd_all, "# %s\n%s\n%s", what,para,head);
    printf("Writing %s in %s\n",what, name);

    strcpy(what, "Predicted DE of all pairs of mutations");
    sprintf(name,"%s.mut_DE_all.dat", nameout1);
    file_de_all=fopen(name, "w");
    fprintf(file_de_all, "# %s\n%s\n%s", what,para,head);
    printf("Writing %s in %s\n",what, name);
  }

  FILE *file3=NULL;
  if(0){
    strcpy(what, "Predicted RMSD and DE of all mutations from WT");
    sprintf(name,"%s.mut_all.dat", nameout1);
    fopen(name, "w");
    printf("Writing %s in %s\n",what, name);
    fprintf(file3, "# %s\nRMSD DE\n", what);
  }

  FILE *file4=NULL;
  if(0){
    sprintf(what, "Predicted force of all mutations from WT");
    sprintf(name,"%s.mut_force.dat", nameout1);
    file4=fopen(name, "w");
    printf("Writing %s in %s\n",what, name);
    fprintf(file4, "# %s\n%s\n%s", what, para, head);
  }

  // Define memory
  float Force[N_Cart], Force_coeff[NM.N];
  int N3=3*Na;
  float pred_mut_Cart[N_Cart], pred_mut_Tors[N_axes],
    pred_mut_3[N3], pred_mut_2[Na]; // pred_mut_tot_2[Na];

  int Nmut=1, Posmut[1]; char AAwt[1], AAmut[1];
  float MSD_mut[20][20], MSD_cross=0, MSD[20];
  float DE_pred[20][20], DE[20];
  double RMSD_prof[Na]; for(int i=0; i<Na; i++)RMSD_prof[i]=0;
  for(int ipos=0; ipos<Nres; ipos++){
    printf("%d ",ipos);
    Posmut[0]=ipos; AAwt[0]=seq[ipos].amm; AAmut[0]='A';
    int wt=seq[ipos].i_aa;
    int mut_neigh[Na], n_cont_mut;
    Mut_neighbor(mut_neigh, Na, AAwt, Posmut, AAmut, Nmut,
		 clist, cnum, Int_list, atoms, Nres);
    char out[80]; int a, b;
    sprintf(out,"%c%s\t%d",AAwt[0],seq[ipos].pdbres, nc[ipos]); //ipos+1
    fprintf(file_rmsd,"%s",out);
    fprintf(file_de,"%s",out);
    if(file4)fprintf(file4,"%s",out);
    for(a=0; a<20; a++){
      if(a==wt){DE[a]=1000; MSD[a]=1000; continue;}
      AAmut[0]=AA_code[a];
      Pred_mut_str(C_size, C_stab, C_dist,
		   DE+a, MSD+a, &MSD_cross,
		   NULL, NULL, NULL, NULL,
		   pred_mut_3, pred_mut_2, NULL, Na,
		   &n_cont_mut, Force, Force_coeff,
		   pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
		   AAwt, Posmut, AAmut, Nmut,
		   clist, cnum, Int_list, seq, Na,
		   atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
		   NULL, NULL, NULL, NULL);
      for(int i=0; i<Na; i++)RMSD_prof[i]+=pred_mut_2[i];
      if(file4){
	double F=0;
	if(AAmut[0]!=AAwt[0]){
	  for(int i=0; i<Na; i++)F+=Force[i]*Force[i];
	  F=sqrt(F/Na);
	}
	fprintf(file4,"\t%.2g",F);
      }
    }
    // Compute mutation to same a.a. as mean of IWT lowest ones
    if(IWT>0){
      DE[wt]=10000; f_sort(DE, 20, i_rank);
      double Y=0; for(a=0; a<IWT; a++)Y+=DE[i_rank[19-a]]; DE[wt]=Y/IWT;
      MSD[wt]=10000; f_sort(MSD, 20, i_rank);
      Y=0; for(a=0; a<IWT; a++)Y+=MSD[i_rank[19-a]]; MSD[wt]=Y/IWT;
    }else{
      double sum=0, norm=0; float *x=DE, p;
      for(a=0; a<20; a++){
	if(a!=wt){p=exp(-x[a]); sum+=x[a]*p; norm+=p;}
      }
      x[wt]=sum/norm;
      sum=0; norm=0; x=MSD;
      for(a=0; a<20; a++){
	if(a!=wt){p=exp(-x[a]); sum+=x[a]*p; norm+=p;}
      }
      x[wt]=sum/norm;
    }
    for(a=0; a<20; a++){
      float RMSD_mut=sqrt(MSD[a]);
      fprintf(file_rmsd,"\t%.2g",RMSD_mut);
      fprintf(file_de,"\t%.3g",DE[a]);
      if(file3)fprintf(file3,"%.2g %.3g\n", RMSD_mut, DE[a]);
    }
    fprintf(file_rmsd,"\n");
    fprintf(file_de,"\n");
    if(file3)fprintf(file3,"\n");
    if(file4)fprintf(file4,"\n");

    if(ALLPAIR){
      // Compute all pairs
      for(a=0; a<20; a++){
	if(a==wt)continue;
	DE_pred[a][wt]=DE[a];
	DE_pred[wt][a]=DE[a];
	MSD_mut[a][wt]=MSD[a];
	MSD_mut[wt][a]=MSD[a];
	AAwt[0]=AA_code[a];
	for(b=a+1; b<20; b++){
	  if(b==wt){continue;}
	  AAmut[0]=AA_code[b];
	  Pred_mut_str(C_size, C_stab, C_dist,
		       &(DE_pred[a][b]), &(MSD_mut[a][b]), &MSD_cross,
		       NULL, NULL, NULL, NULL,
		       pred_mut_3, pred_mut_2, NULL, Na,
		       &n_cont_mut, Force, Force_coeff,
		       pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes, 
		       AAwt, Posmut, AAmut, Nmut,
		       clist, cnum, Int_list, seq, Na,
		       atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
		       NULL, NULL, NULL, NULL);
	}
      } // end pairs
      // print all pairs
      char out[20];
      sprintf(out,"%c%s",AA_code[wt],seq[ipos].pdbres); //ipos+1
      fprintf(file_rmsd_all,"%s", out);
      fprintf(file_de_all,"%s", out);
      for(a=0; a<20; a++){
	for(b=a+1; b<20; b++){
	  fprintf(file_rmsd_all,"\t%.2g",sqrt(MSD_mut[a][b]));
	  fprintf(file_de_all,"\t%.3g",DE_pred[a][b]);
	}
      }
      fprintf(file_rmsd_all,"\n");
      fprintf(file_de_all,"\n");
    } // end ALLPAIR
  }
  fclose(file_rmsd);
  fclose(file_de);
  if(file3)fclose(file3);
  if(file4)fclose(file4);
  if(file_de_all)fclose(file_de_all);
  if(file_rmsd_all)fclose(file_rmsd_all);

  // Print profile RMSD
  int NMUT=19*Nres;
  sprintf(name,"%s.mut_prof_RMSD.dat",nameout1);
  FILE *file1=fopen(name,"w");
  printf("Printing profile of structural RMSD in %s\n",name);
  fprintf(file1,"# Profile of structural deviations");
  fprintf(file1," averaged over %d simulated point mutations\n",NMUT);
  fprintf(file1,"# 1=res 2=n_contact 3=predicted RMSD from w.t.\n");
  for(int i=0; i<Na; i++){
    fprintf(file1,"%c%d\t%d\t%.3f\n",
	    seq[i].amm, i+1, nc[i], sqrt(RMSD_prof[i]/NMUT));
  }
  fclose(file1);
  printf("\n");
}

int Extract_reference_atoms(int *resatom, atom *atoms, int natoms,
			    struct residue *seq, int Nres)
{
  char SEL[4]="CA"; if(strcmp(REF_CC, "CB")==0)strcpy(SEL, "CB");
  int Na=0; 
  // resatom[i]= Reference atom for residue i
  // Center of mass and inertia tensor

  for(int i=0; i<natoms; i++){
    atom *atm=atoms+i;
    if((strncmp(atm->name, SEL, 2)!=0)&&
       ((seq[atm->res].amm!='G')||(strncmp(atm->name, "CA", 2)!=0)))continue;
    resatom[Na]=i; Na++;
  }
  if(Na!=Nres)printf("WARNING, %d residues expected, %d found\n", Nres, Na);
  if(Na > Nres){printf("Leaving\n"); return(0);}
  return(Na);
}

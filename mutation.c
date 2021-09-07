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

int OPT_COEFF=0;
char OPT_SCORE='C'; // C= Cartesian correlation; A=absolute value correlation;
float Lambda=0.10; // Penalization in optimization
float EXP_FORCE=0, Exp;
float COLL_THR_MUT=0.01; // Minimal collectivity of normal modes to predict mut
//float C_SIZE=0.015, C_STAB=1.4, C_DIST=1.0;
float C_SIZE=40, C_STAB=300, C_DIST=52;
int PRINT_FLANK=1;

int ncmax=100; 
int mut_cont[100], wt_cont[100];
int mut_res[100], wt_res[100];

float Econt_norm[21][21], D_Res_norm[21][21], natm_norm[21];

// Econt D_Residues atom_numb

void Read_pred_para(char *name_in);
void Normalize_parameters();
float Pred_mut_str(float C_size, float C_stab, float C_dist,
		   double *DE, float *RMSD_pred, int *n_cont_mut,
		   float *Force, float *Force_coeff,
		   float *pred_mut_Cart, int N_Cart, 
		   float *pred_mut_Tors, int N_axes,
		   float *pred_mut_3, float *pred_mut_abs, int Na,
		   char *AAwt, int *Posmut, char *AAmut, int Nmut,
		   int **clist, int **cnum, struct interaction *Int_list,
		   struct residue *seq, int Nres, atom *atoms, int *resatom,
		   struct Reference Ref_kin,
		   int *kin_atom, struct Normal_Mode NM,
		   float *sigma2, float *str_diff_abs, float *str_diff_3);
int Extract_reference_atoms(int *resatom, atom *atoms, int natoms,
			    struct residue *seq, int Nres);
//float *exc_3
float Optimize_coeff(float *C_size, float *C_stab, float *C_dist,
		     int N_Cart,  int N_axes, int Na,
		     char *AAwt, int *Posmut, char *AAmut, int Nmut,
		     int **clist, int **cnum, struct interaction *Int_list,
		     struct residue *seq, int Nres, atom *atoms, int *resatom,
		     struct Reference Ref_kin,
		     int *kin_atom, struct Normal_Mode NM, float *sigma2,
		     float *str_diff_abs, float *str_diff_3); //, float *exc_3

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


void Compute_mut_def(float *pred_mut_Cart, int N_Cart,
		     float *pred_mut_Tors, int N_axes,
		     double *DE, float *Force_coeff, float *Force,
		     struct Normal_Mode NM, float *sigma2, float mass);

float *Convert_to_tors(float *Cart, int N3, float **J_ar, int N_axes);
void  Match_atoms(int *ref_atom, atom *atoms, int natoms,
		  int *atom_ref, int N_ref);
float Cosine(float *scale, float *xx, float *yy, int n);

int Mutation(float *mut_Tors_out, int N_axes, float *mut_CC_out,
	     struct Reference Ref_kin, 
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

  /****************************************************
           Find mutated residues
 ******************************************************/
  printf("Examining the effect of mutations: ");
  // Interactions directed along representative atoms
  int Nres=atoms[Ref_kin.atom_num[Ref_kin.N_ref-1]].res+1;
  int i, imut=1;
  char mut_string[200]="";
  for(i=0; i<Nmut; i++){
    int pos=Posmut[i];
    if((pos<0)||(pos>=Nres)||(Code_AA(AAmut[i])<0)||(Code_AA(AAmut[i])>20)){
      printf("\nWARNIG mutation does not exist (only %d residues)\n", Nres);
      imut=0;
    }
    sprintf(mut_string, "%s%c%d%c", mut_string, AAwt[i], pos, AAmut[i]);
    if(i < (Nmut-1))sprintf(mut_string, "%s_", mut_string);
  }
  printf(" %s\n", mut_string);
  if(imut==0)return(0);

  // Anharmonicity
  float *sigma2; char name_har[80];
  if(anhar==0){sigma2=NM.sigma2; strcpy(name_har, ""); }
  else{sigma2=NM.sigma2_anhar; strcpy(name_har, "_anharmonic");}
  
  /**************************************************************************
           Normalize size, distance and energy param. for calculations
           Set parameters for gaps
 ****************************************************************************/
  Normalize_parameters();

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
  float str_diff_3[3*Nres], str_diff_abs[Nres];

  // Match any atom to closest kinetic atom
  int N_kin=Ref_kin.N_ref, N_Cart=3*N_kin, kin_atom[natoms]; 
  Match_atoms(kin_atom, atoms, natoms, Ref_kin.atom_num, N_kin);

  // Structural deformation profile for reference atoms
  // Extract one reference atom per residue
  char SEL[4]="CA"; if(strcmp(REF_CC, "CB")==0)strcpy(SEL, "CB");
  int Na=0, j, resatom[Nres]; double m_tot=0;
  // resatom[i]= Reference atom for residue i
  // Center of mass and inertia tensor
  float str_diff_ave[3]; for(j=0; j<3; j++)str_diff_ave[j]=0;
  //double **corr_sum=Allocate_mat2_d(3, 3);
  //double **inertia=Allocate_mat2_d(3, 3);

  float *str=str_diff_3; 
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
  double RMSD=0;
  str=str_diff_3;
  for(i=0; i<Na; i++){
    double d2=0;
    for(j=0; j<3; j++){
      *str-=str_diff_ave[j]; d2+=(*str)*(*str); str++;
    }
    str_diff_abs[i]=sqrt(d2); RMSD+=d2;
  }
  // Predicted fluctuations B1
  double B_norm=0; for(i=0; i<Na; i++)B_norm+=B_pred[i];
  B_norm=sqrt(B_norm/RMSD);
  float B1[Na];
  for(i=0; i<Na; i++)B1[i]=sqrt(B_pred[i])/B_norm;
  RMSD=sqrt(RMSD/Na);
  printf("RMSD(%s atoms)= %.2f\n", SEL, RMSD);

  /******************** Compute excess confchange ********************/
  float excess[Na], exc_3[3*Na], max_exc_mut=0, max_exc_nomut=0;

  // Contact list
  int *nc, **clist, **cnum;
  Get_contact_list(&nc, &clist, &cnum, Nres, atoms, Int_list, N_int);

  int mut_neigh[Na],
    n_neigh=Mut_neighbor(mut_neigh, Na, AAwt, Posmut, AAmut, Nmut,
			 clist, cnum, Int_list, atoms, Nres);
  j=0;
  for(i=0; i<Na; i++){
    excess[i]=str_diff_abs[i]/B1[i];
    if(mut_neigh[i]){
      if(excess[i]>max_exc_mut)max_exc_mut=excess[i];
    }else if(excess[i]>max_exc_nomut){
      max_exc_nomut=excess[i];
    }
    for(int k=0; k<3; k++){exc_3[j]=str_diff_3[j]/B1[i]; j++;}
  }

  /************************************************************************
           Prepare output
  *************************************************************************/
  // Header
  char header[1000];
  sprintf(header, "# RMSD(%s atoms)= %.2f N= %d residues n= %d kinetic atoms\n",
	  SEL, RMSD, Nres, Ref_kin.N_ref);
  sprintf(header, "%s# Prediction based on %d mutations: %s\n",
	  header, Nmut, mut_string);
  if(Nres!=Na)sprintf(header, "%s# gap= %d residues\n", header, Nres-Na);

  //Correlation predicted_fluctuations confchange (Cartesian):
  float slope, offset;
  float r_CC_B=Corr_coeff(B1, str_diff_abs, Na, &slope, &offset);
  sprintf(header, "%s# Correlation predicted_fluct confchange: %.3f",
	  header, r_CC_B);
  sprintf(header, "%s slope= %.3g offset= %.3g\n", header, slope, offset);

  //Correlation predicted_fluctuations confchange (torsional):
  float r=Corr_coeff(Tors_fluct, Confchange_Tors, N_axes, &slope, &offset);
  sprintf(header, "%s# Correlation predicted_Tors_fluct Tors_confchange:",
	  header);
  sprintf(header, "%s %.3f slope=%.2g offset=%.2g\n", header, r,slope,offset);

  sprintf(header, "%s# %d neighbors of mutated residues\n",
	  header, n_neigh);
  sprintf(header, "%s# Mut: Max.conf_change(obs/pred)= %.3f\n",
	  header, max_exc_mut);
  sprintf(header, "%s# No-mut: Max.conf_change(obs/pred)= %.3f\n",
	  header, max_exc_nomut);

  // Parameters of the prediction
  Read_pred_para(Mut_para); //"Mutation_para.in"
  char nameout[200];
  sprintf(nameout, "%s%s_EXPF%.2f_SC_%c_pred_mut.dat",
	  nameout1, name_har, EXP_FORCE, OPT_SCORE);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "%s", header);
  fprintf(file_out, "# Power-law of the force: r^-%.1f\n", EXP_FORCE);

  // Columns printed in the output file for each prediction type
  char head[800];
  sprintf(head, "# Pred.force directed along axes between kinetic atoms\n");
  sprintf(head, "%s#1=RMSD 2=RMSD_pred ", head);
  sprintf(head, "%s3=Cosine(pred_mut,confchange) 4=scale ", head);
  sprintf(head, "%s5=Corr(|pred_mut|,|confchange|) 6=slope 7=offset ", head);
  sprintf(head, "%s8=Cosine(pred_mut,excess) 9=scale ", head);
  sprintf(head, "%s10=Corr(|pred_mut|,excess) 11=slope 12=offset ", head);
  sprintf(head, "%s13=Corr(pred_mut,confchange)_tors\n", head); 
  //14=slope 15=offset

  /**************************************************************************
           Compute predicted mutations
 ****************************************************************************/

  // Exponent of the force for computation
  Exp=(1.+EXP_FORCE)/2;

  // 1 type of prediction based on combination of size, stability, distance 
  int N_TYPES=0; float pred_mut_abs_f[N_TYPES+1][Na];
  // Coefficients for the computation of the force
  float C_stab=0, C_size=0, C_dist=0;
  if(C_SIZE==0){C_SIZE=1;}
  if(C_DIST==0){C_DIST=1;}
  if(C_STAB==0){C_STAB=1;}

  int N3=3*Na;
  float Force[N_Cart], Force_coeff[NM.N_relevant]; 
  float pred_mut_Cart[N_Cart], pred_mut_Tors[N_axes],
    pred_mut_abs[Na], pred_mut_3[N3];

  // Loop on predictions
  char pred_type[400];
  float cc_opt=0; double DE=0;
  int N_force=N_TYPES+1; //if(OPT_COEFF==0)N_force++;
  for(int iforce=0; iforce<N_force; iforce++){
    sprintf(pred_type, "Predicted mut based on changes of");
    if(iforce>=N_TYPES){
      sprintf(pred_type,
	      "%s combination of size, stability and optimal distance",
	      pred_type);
      C_size=C_SIZE;  C_stab=C_STAB; C_dist=C_DIST;
    }else if(iforce==0){
      sprintf(pred_type, "%s size, only s.c. contacts", pred_type);
      C_size=C_SIZE;  C_stab=0; C_dist=0;
    }else if(iforce==1){
      sprintf(pred_type, "%s stability", pred_type);
      C_size=0;  C_stab=C_STAB; C_dist=0;
    }else if(iforce==2){
      sprintf(pred_type, "%s opt_distance", pred_type);
      C_size=0;  C_stab=0; C_dist=C_DIST;
    }
    printf("%s\n", pred_type);
    
    // Compute the force due to the mutation
    float RMSD_pred, scale, cc, cc_all; int n_cont_mut;
    if(OPT_COEFF && (iforce==N_TYPES)){
      // Optimize coefficients
      fprintf(file_out, "# Starting optimization, ");
      fprintf(file_out, "C_SIZE= %.3g C_STAB= %.3g C_DIST= %.3g\n",
	      C_size, C_stab, C_dist);
      cc_opt=Optimize_coeff(&C_size, &C_stab, &C_dist,
			    N_Cart, N_axes, Na,
			    AAwt, Posmut, AAmut, Nmut,
			    clist, cnum, Int_list, seq, Nres,
			    atoms, resatom, Ref_kin, kin_atom, NM,
			    sigma2, str_diff_abs, str_diff_3);
    }
    cc=Pred_mut_str(C_size, C_stab, C_dist, &DE, &RMSD_pred, &n_cont_mut,
		    Force, Force_coeff, pred_mut_Cart, N_Cart,
		    pred_mut_Tors, N_axes, pred_mut_3, pred_mut_abs, Na,
		    AAwt, Posmut, AAmut, Nmut,
		    clist, cnum, Int_list, seq, Nres,
		    atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
		    str_diff_abs, str_diff_3);  //, exc_3
    if(cc>cc_opt){cc_opt=cc;}
    cc_all=Cosine(&scale, pred_mut_3, str_diff_3, N3);
    if(OPT_COEFF){
      // Normalize predictions by fitted scale
      if(C_size){C_size*=scale; C_SIZE=C_size;}
      if(C_stab){C_stab*=scale; C_STAB=C_stab;}
      if(C_dist){C_dist*=scale; C_DIST=C_dist;}
      for(i=0; i<N3; i++)pred_mut_3[i]*=scale;
      for(i=0; i<Na; i++)pred_mut_abs[i]*=scale;
      for(i=0; i<N_axes; i++)pred_mut_Tors[i]*=scale;
      RMSD_pred*=scale;
    }
    printf("Mut_pred_score: %.4g  C_size=%.3g C_stab=%.3g C_dist=%.3g\n",
	   cc, C_size, C_stab, C_dist);

    /***********************************************************
                Print and store the prediction
    *******************************************************/
    // Store
    int kforce=iforce; //if(iforce==4)kforce=3;
    for(i=0; i<Na; i++)pred_mut_abs_f[kforce][i]=pred_mut_abs[i];

    // Print
    printf("Predicted RMSD= %.2f\n", RMSD_pred); 
    // Print contacts
    if(iforce==0){
      fprintf(file_out, "# Prediction based on %d contacts: ", n_cont_mut);
      for(i=0; i<n_cont_mut; i++){
	fprintf(file_out, " %c%d-%c%d",
		AA_code[wt_cont[i]], wt_res[i],
		AA_code[mut_cont[i]], mut_res[i]);
	if(i==ncmax)break;
      }
      fprintf(file_out, "\n");
      fprintf(file_out, "%s", head);
    }

    // Print type of prediction and parameters
    fprintf(file_out, "# %s\n", pred_type);
    if(OPT_COEFF)
      fprintf(file_out, "# Parameters were rescaled to fit the RMSD\n");
    fprintf(file_out, "# C_SIZE= %.3g C_STAB= %.3g C_DIST= %.3g\n",
	    C_size, C_stab, C_dist);
    fprintf(file_out, "%.2f\t%.2f", RMSD, RMSD_pred);
    fprintf(file_out, "\t%.3g\t%.2g", cc_all, scale);

    float slope, offset, r_CC_A=
      Corr_coeff(pred_mut_abs, str_diff_abs, Na, &slope, &offset);
    fprintf(file_out, "\t%.3f\t%.2g\t%.2g", r_CC_A, slope, offset);

    /*
      float r_A_B=Corr_coeff(B_pred, pred_mut_abs2, Na, &slope, &offset);
      fprintf(file_out,"\t%.3f", r_A_B);
      //Partial correlation predicted_mut confchange| B_pred:
      float rp=(r_CC_A-r_CC_B*r_A_B)/sqrt((1-r_CC_B*r_CC_B)*(1-r_A_B*r_A_B));
      fprintf(file_out, "\t%.3f", rp);
    */
  
    // Cosine predicted-excess confchange
    float scale_exc, cc_exc=Cosine(&scale_exc, pred_mut_3, exc_3, N3); 
    fprintf(file_out, "\t%.3g\t%.2g", cc_exc, scale_exc);
  
    float r=Corr_coeff(pred_mut_abs, excess, Na, &slope, &offset);
    fprintf(file_out, "\t%.3f\t%.2g\t%.2g", r, slope, offset);
    printf("Correlation predicted_mutation excess (abs): %.3f\n",r);
  
    //Correlation predicted_mutation confchange (torsional):
    for(i=0; i<N_axes; i++)pred_mut_Tors[i]*=pred_mut_Tors[i];
    r=Corr_coeff(pred_mut_Tors, Confchange_Tors, N_axes, &slope, &offset);
    fprintf(file_out, "\t%.3f", r); //\t%.2g\t%.2g slope, offset
    fprintf(file_out, "\n");


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
    for(i=0; i<N_Cart; i++)printf(" %.2f", Force[i]); printf("\n");
    printf("mut_Cart: ");
    for(i=0; i<N_Cart; i++)printf(" %.2f",pred_mut_Cart[i]);printf("\n");
    printf("mut_CC: ");
    for(i=0; i<Ref_CC.N_cart; i++)printf(" %.2f", mut_CC_out[i]);
    printf("\n");
    printf("mut_Tors: ");
    for(i=0; i<N_axes; i++)printf(" %.2f", mut_Tors_out[i]); printf("\n");
  }


  fclose(file_out);

  /******************************************************************/
  // Write excess_confchange and predicted mutations for every residue
  sprintf(nameout, "%s%s_EXPF%.2f_SC_%c_mutation.dat",
	  nameout1, name_har, EXP_FORCE, OPT_SCORE);
  file_out=fopen(nameout, "w");
  fprintf(file_out, "%s", header);
  fprintf(file_out,"#1=conf_ch 2=pred_fluct 3=conf_ch/pred 4=mut_neighbor");
  int k; for(k=1; k<=4; k++)fprintf(file_out, " %d=pred_mut_%d", k+4, k); 
  fprintf(file_out, " pos aa\n");
  for(i=0; i<Na; i++){
    fprintf(file_out, "%.3g\t%.3g\t%.4g\t%d",
	    str_diff_abs[i], B1[i], excess[i], mut_neigh[i]);
    for(k=0; k<N_TYPES; k++)fprintf(file_out, "\t%.3g", pred_mut_abs_f[k][i]);
    int a=atoms[resatom[i]].res;
    fprintf(file_out, "\t%s\t%c\n", seq[a].pdbres, seq[a].amm);
  }
  fclose(file_out);
  printf("Writing %s\n", nameout);

  /******************************************************************/
  // Clean memory
  Empty_matrix_i(clist, Nres);
  Empty_matrix_i(cnum, Nres); free(nc);
  return(1);
}

float Pred_mut_str(float C_size, float C_stab, float C_dist,
		   double *DE, float *RMSD_pred, int *n_cont_mut,
		   float *Force, float *Force_coeff,
		   float *pred_mut_Cart, int N_Cart, 
		   float *pred_mut_Tors, int N_axes,
		   float *pred_mut_3, float *pred_mut_abs, int Na,
		   char *AAwt, int *Posmut, char *AAmut, int Nmut,
		   int **clist, int **cnum, struct interaction *Int_list,
		   struct residue *seq, int Nres, atom *atoms, int *resatom,
		   struct Reference Ref_kin,
		   int *kin_atom, struct Normal_Mode NM,
		   float *sigma2,
		   float *str_diff_abs, float *str_diff_3) //, float *exc_3
{
  // Compute force
  *n_cont_mut=
    Mutation_force(Force, N_Cart, C_size, C_stab, C_dist, 
		   AAwt, Posmut, AAmut, Nmut, clist, cnum, Int_list,
		   seq, Nres, atoms, Ref_kin, kin_atom);

  // Compute deformation mut_Cart and mut_Tors
  Compute_mut_def(pred_mut_Cart, N_Cart, pred_mut_Tors, N_axes,
		  DE, Force_coeff, Force, NM, sigma2, Ref_kin.mass_tot);
      
  // Deformation of resatom mut_3, mut_abs2 and RMSD_pred
  float pred_mut_ave[3], *str=pred_mut_3; int i, j;
  for(j=0; j<3; j++)pred_mut_ave[j]=0;
  for(i=0; i<Na; i++){
    float *x=pred_mut_Cart+3*kin_atom[resatom[i]];
    for(j=0; j<3; j++){*str=*x; pred_mut_ave[j]+=*x; x++; str++;}
  }
  // Subtract center of mass
  for(j=0; j<3; j++)pred_mut_ave[j]/=Na;
  *RMSD_pred=0; str=pred_mut_3;
  for(i=0; i<Na; i++){
    double d2=0;
    for(j=0; j<3; j++){*str-=pred_mut_ave[j]; d2+=(*str)*(*str); str++;}
    *RMSD_pred+=d2; pred_mut_abs[i]=sqrt(d2);
  }
  *RMSD_pred=sqrt(*RMSD_pred/Na);

  if((str_diff_3==NULL)||(str_diff_abs==NULL))return(0);

  float score=0; //-Lambda*(C_size*C_size+C_stab*C_stab+C_dist*C_dist);
  float scale, slope, offset;
  if(OPT_SCORE=='C'){
    score+=Cosine(&scale, pred_mut_3, str_diff_3, 3*Na);
    //score+=Cosine(&scale, pred_mut_3, exc_3, N3);
  }else{
    score+=Corr_coeff(pred_mut_abs, str_diff_abs, Na, &slope, &offset);
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
	  if(pmut<0)continue; p1=pmut;
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
		     double *DE, float *Force_coeff, float *Force,
		     struct Normal_Mode NM, float *sigma2, float mass)
{
  int a, i;
  for(i=0; i<N_Cart; i++)pred_mut_Cart[i]=0;
  for(i=0; i<N_axes; i++)pred_mut_Tors[i]=0;

  /*double F=0;
  for(i=0; i<N_Cart; i++){F+=Force[i]*Force[i];}
  printf("RMS of mutational force: %.3f\n", sqrt(3*F/N_Cart)); */
  *DE=0;

  // Compute predicted deformation of kinetic atoms (Cart)
  for(a=0; a<NM.N_relevant; a++){
    if((NM.Cart_coll[a]<COLL_THR_MUT)||(sigma2[a]<=0)){
      Force_coeff[a]=0; continue;
    }
    float *f=Force, *x=NM.Cart[a];
    double xf=0; // Projection of force on normal mode
    for(i=0; i<N_Cart; i++){xf+=(*x)*(*f); x++; f++;}
    //(*DE)+=xf*xf/sigma2[a];
    (*DE)+=xf*xf*sigma2[a];
    xf*=sigma2[a];
    Force_coeff[a]=xf;
    float *m=pred_mut_Cart; x=NM.Cart[a];
    for(i=0; i<N_Cart; i++){
      (*m)+=(*x)*xf; x++; m++;
    }
    m=pred_mut_Tors; x=NM.Tors[a];
    for(i=0; i<N_axes; i++){
      (*m)+=(*x)*xf; x++; m++;
    }
  }
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
    if(ref_atom[i]>=0)continue; n++;
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
  printf("Reading parameters in %s\n", name_in);
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
  fclose(file_in);
}

float Optimize_coeff(float *C_size_opt, float *C_stab_opt, float *C_dist_opt,
		     int N_Cart,  int N_axes, int Na,
		     char *AAwt, int *Posmut, char *AAmut, int Nmut,
		     int **clist, int **cnum, struct interaction *Int_list,
		     struct residue *seq, int Nres, atom *atoms, int *resatom,
		     struct Reference Ref_kin,
		     int *kin_atom, struct Normal_Mode NM, float *sigma2,
		     float *str_diff_abs, float *str_diff_3) //, float *exc_3
{
  float cc_opt=-1000; double DE=0;
  int N3=3*Na;
  float Force[N_Cart], pred_mut_Cart[N_Cart], pred_mut_Tors[N_axes],
    pred_mut_abs[Na], pred_mut_3[N3];
  float Force_coeff[NM.N_relevant];
  float RMSD_pred; int n_cont_mut;

  // Initialize parameters
  float C_size=*C_size_opt, C_stab=*C_stab_opt, C_dist=*C_dist_opt;

  float x[3], y[3], yy, x_min, x_max;
  float *C_move, y_max=cc_opt;
  int ifit=0, j, ini=0, n_fail=0;
  printf("Optimizing mutation coefficients\n");
  for(int it=0; it<200; it++){
    if(ifit==0){C_move=&C_size; x[1]=*C_size_opt; ifit=1;}
    else{C_move=&C_dist; x[1]=*C_dist_opt; ifit=0;}
    x_min=x[1]*0.1; x_max=x[1]*10;
    for(j=0; j<3; j++){
      if(j==0){x[0]=x[1]*0.85;}
      else if((j==1)&& ini){y[1]=y_max; continue;}
      else{x[2]=x[1]*1.15;}
      *C_move=x[j];
      y[j]=Pred_mut_str(C_size, C_stab, C_dist, &DE, &RMSD_pred, &n_cont_mut,
			Force, Force_coeff, pred_mut_Cart, N_Cart,
			pred_mut_Tors, N_axes, pred_mut_3, pred_mut_abs, Na,
			AAwt,Posmut,AAmut,Nmut, clist,cnum,Int_list,seq,Nres,
			atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
			str_diff_abs, str_diff_3); // exc_3
      if((ini==0)&&(j==1)){
	printf("mut_pred_score: %.4g    ifit=%d it=%d %.3g %.3g %.3g\n",
	       y[j], ifit, it, C_size, C_stab, C_dist);
	ini=1;
      }
    }
    // Determine new value of x and test its score
    *C_move=Find_max_quad(x[0], x[1], x[2], y[0], y[1], y[2], x_min, x_max);
    yy=Pred_mut_str(C_size, C_stab, C_dist, &DE, &RMSD_pred, &n_cont_mut,
		    Force, Force_coeff, pred_mut_Cart, N_Cart,
		    pred_mut_Tors, N_axes, pred_mut_3, pred_mut_abs, Na,
		    AAwt, Posmut, AAmut, Nmut, clist, cnum, Int_list,seq,Nres,
		    atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
		    str_diff_abs, str_diff_3); // exc_3
    if(yy>y_max){
      y_max=yy; n_fail=0;
      if(yy>cc_opt){
	cc_opt=yy;
	*C_size_opt=C_size; *C_stab_opt=C_stab; *C_dist_opt=C_dist;
      }
    }else{
      // Score does not improve. If two consecutive failures, exit.
      n_fail++;
      if(n_fail>=2){
	printf("%d consecutive decreases, exiting\n", n_fail); break;
      }
    }
    printf("mut_pred_score: %.4g    ifit=%d it=%d %.3g %.3g %.3g\n",
	   yy, ifit, it, C_size, C_stab, C_dist);
  } // end iterations
  printf("mut_pred_score: %.4g  C_size=%.3g C_stab=%.3g C_dist=%.3g\n",
	 cc_opt, *C_size_opt, *C_stab_opt, *C_dist_opt);
  return(cc_opt);
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
  printf("Econt: %.3g %.3g\n", E1, E2);
  printf("D_Res: %.3g %.3g\n", D1, D2);
  printf("natom: %.3g %.3g\n", na1, na2);
}

void Predict_mutations(struct Normal_Mode NM, atom *atoms, int natoms,
		       int N_axes, struct Reference Ref_kin, 
		       struct interaction *Int_list, int N_int,
		       struct residue *seq, int Nres,
		       char *nameout1, char *Mut_para)
{
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

  // Exponent of the force for computation
  Exp=(1.+EXP_FORCE)/2;
  Read_pred_para(Mut_para);
  float C_size=C_SIZE, C_stab=C_STAB, C_dist=C_DIST;
  Normalize_parameters();

  // Output files
  char head[400];
  sprintf(head, "# Mutation par.: C_SIZE= %.3g C_STAB= %.3g C_DIST=%.3g\n",
	  C_size, C_stab, C_dist);
  sprintf(head, "%s#Mut\tncont", head);
  for(int a=0; a<20; a++)sprintf(head, "%s\t%c", head, AA_code[a]);
  sprintf(head,"%s\n", head);
  char name1[200]; sprintf(name1,"%s.mut_RMSD.dat", nameout1);
  char what1[200]="Predicted RMSD of all possible mutations";
  FILE *file1=fopen(name1, "w");
  printf("Writing %s in %s\n",what1,name1);
  fprintf(file1, "# %s\n%s", what1,head);
  char name2[200]; sprintf(name2,"%s.mut_DE.dat", nameout1);
  char what2[200]="Predicted DE of all possible mutations";
  FILE *file2=fopen(name2, "w");
  printf("Writing %s in %s\n",what2, name2);
  fprintf(file2, "# %s\n%s", what2,head);

  FILE *file3=NULL;
  if(0){
    char what3[200]="Predicted RMSD and DE of all possible mutations";
    char name3[200]; sprintf(name3,"%s.mut_all.dat", nameout1);
    fopen(name3, "w");
    printf("Writing %s in %s\n",what3, name3);
    fprintf(file3, "# %s\nRMSD DE\n", what3);
  }

  // Define memory
  float Force[N_Cart], Force_coeff[NM.N_relevant];
  float pred_mut_Cart[N_Cart], pred_mut_Tors[N_axes],
    pred_mut_3[3*Na], pred_mut_abs[Na];

  int Nmut=1, Posmut[1]; char AAwt[1], AAmut[1];
  float RMSD_pred=0; double DE_pred=0;
  double RMSD_prof[Na]; for(int i=0; i<Na; i++)RMSD_prof[i]=0;
  for(int ipos=0; ipos<Nres; ipos++){
    printf("%d ",ipos);
    Posmut[0]=ipos; AAwt[0]=seq[ipos].amm; AAmut[0]='A';
    int mut_neigh[Na], n_cont_mut;
    Mut_neighbor(mut_neigh, Na, AAwt, Posmut, AAmut, Nmut,
		 clist, cnum, Int_list, atoms, Nres);
    char out[80]; sprintf(out,"%c%d\t%d",AAwt[0],ipos+1, nc[ipos]);
    fprintf(file1,"%s",out);
    fprintf(file2,"%s",out);
    for(int a=0; a<20; a++){
      AAmut[0]=AA_code[a];
      if(AAmut[0]==AAwt[0]){RMSD_pred=0; DE_pred=0;}
      else{
	Pred_mut_str(C_size, C_stab, C_dist, &DE_pred, &RMSD_pred, &n_cont_mut,
		     Force, Force_coeff, pred_mut_Cart, N_Cart,
		     pred_mut_Tors, N_axes, pred_mut_3, pred_mut_abs, Na,
		     AAwt, Posmut, AAmut, Nmut,
		     clist, cnum, Int_list, seq, Na,
		     atoms, resatom, Ref_kin, kin_atom, NM, sigma2,
		     NULL, NULL);
	for(int i=0; i<Na; i++)RMSD_prof[i]+=pred_mut_abs[i];
      }
      fprintf(file1,"\t%.2f",RMSD_pred);
      fprintf(file2,"\t%.3g",DE_pred);
      if(file3)fprintf(file3,"%.2f %.3g\n", RMSD_pred, DE_pred);
    }
    fprintf(file1,"\n");
    fprintf(file2,"\n");
  }
  fclose(file1);
  fclose(file2);
  if(file3)fclose(file3);
  // Print profile RMSD
  int NMUT=19*Nres;
  sprintf(name1,"%s_mut_prof_RMSD.dat",nameout1);
  file1=fopen(name1,"w");
  printf("Printing profile of structural RMSD in %s\n",name1);
  fprintf(file1,"# Profile of structural deviations");
  fprintf(file1," averaged over %d simulated point mutations\n",NMUT);
  fprintf(file1,"# 1=res 2=n_contact 3=predicted RMSD from w.t.\n");
  for(int i=0; i<Na; i++){
    fprintf(file1,"%c%d\t%d\t%.3f\n",
    seq[i].amm, i+1, nc[i], RMSD_prof[i]/NMUT);
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

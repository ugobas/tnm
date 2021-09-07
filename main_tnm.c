#define PRG "tnm"
int LABEL=0; // Print model parameters in file name?

/********** Some of these parameters are not set from input file **********/
// Fit of B factors
float KAPPA_DEF=218.4; // Force const. for parameters C=4.5 E=6 if B not fitted

// Selection of modes
#define COLL_THR_DEF 30   
float COLL_THR= COLL_THR_DEF;  // Select mode if exp(-S_cart) > COLL_THR
// E_MIN

// Analysis of conformation change
int SWITCH_BONDS=1; // Select only bonds with aligned atoms (1)
int MIN_ALI_LEN_DEF=20; // Chains are aligned if L>MIN_ALI_LEN, if not peptides
float RMSD_MAX=100;   // max. RMSD for analyzing conformation change
float RMSD_MIN=0.5;   // min. RMSD for analyzing conformation change
#define RMSD_THR_DEF 0.1  // If confchange projections smaller, discard
#define SEQID_THR_DEF 50  // If Seq.Id<THR, align structures with Mammoth
int NMUT=-2; /* Mutations between the two proteins required for analysis.
		NMUT=-2: No analysis. NMUT=-1: Always analysis.
		NMUT=0: Sequences must be identical.
		NMUT=1: Exactly one mutation */
int PRED_MUT=0; // Predict mutations?
//int OPT_MUT=0; // Optimize coefficients to fit RMSD of mutation?
char Mut_para[100]="Mutation_para.in";

// Simulation of conformation changes
#define AMPLITUDE_DEF    1.0  // NM amplitude for simulated structures
#define E_THR_DEF        20.0 // Threshold for stopping motion
#define E_THR1_DEF       1.0  // Threshold for deciding privileged direction
#define D_REP_DEF        2.5  // Threshold for repulsions
#define MAX_ANGLE_DEF    0.4  // Max. angle for torsional deformations, radiants
float s0=0.07;       // Entropy per residue, used for simulations

// Analysis performed
int ANHARMONIC=0;
int SITE_DYNAMICS=0; // Examine dynamics of binding sites


// Output
#define PRINT_PDB_DEF 1
#define PRINT_MODE_SUMM_DEF 1
#define PRINT_FORCE_DEF 0
#define PDB_STEP_DEF  0.4  // Step for printing PDB


// Dynamical couplings
int PRINT_DEF_COUPLING=0;  // Print deformation coupling?
int PRINT_DIR_COUPLING=0;  // Print directionality coupling?
int PRINT_COV_COUPLING=0;  // Print covariance coupling?
int PRINT_COORD_COUPLING=0; // Print coordination coupling?
int PRINT_SIGMA_DIJ=0;  // Print variance of DIJ distanc
int PRINT_CMAT=0;       // Print contact matrix
int STRAIN=0; // Compute strain profile, similar to PNAS 106:14253 2009 ?
int ALL_PAIRS=0; // Print couplings for all pairs?
float SIGMA=1.0; // Print only couplings > SIGMA*std.dev.
char PROF_TYPE='A';
char SITES[100]="";

// Default parameters
// Model parameters
#define REF_DEF "ALL"      // Ref. atoms Allowed: CA BB(backbone) EB ALL
#define INT_TYPE_DEF "MIN" // Interaction type. Allowed: CA CB ALL MIN
#define THR_CA 9.0         // For CA contacts
#define THR_CB 9.0         // For CB contacts
#define THR_ALL 4.5        // For ALL and MIN atoms contacts
#define EXP_HESSIAN_DEF 6  // force constant ~ r^(-EXP_HESSIAN)
#define POW_DEF 1          // Power-law or exponential decay of force constant?
#define K_OMEGA_DEF 0.5      // Force constant for omega torsions
#define K_PSI_DEF 0.2      // Force constant for psi torsions
#define K_PHI_DEF 0.2      // Force constant for phi torsions
#define K_CHI_DEF 0.4      // Force constant for chi torsions
#define K_BA_DEF 0.5       // Force constant for bond angles
#define K_BL_DEF 0.5       // Force constant for bond length
#define MIN_INT_DEF 1      // Freeze unconstrained degrees of freedom


// Elastic constants defined in nma_para.h

#ifndef __USE_GNU   // this is to use the get_current_dir_name
#define __USE_GNU
#endif

#include "coord.h"
#include "tnm.h"
#include "nma_para.h"
#include "nma.h"
#include "read.h"
#include "dof_tnm.h"
#include "kinetic_tnm.h"
#include "interactions_tnm.h"
#include "vector.h"
#include "buildup.h"
#include "allocate.h"
#include "McLachlan.h"
#include "align_tnm.h"
#include "output_tnm.h"
#include "contacts.h"
#include "energy_BKV.h"
#include "force2confchange.h"
#include "simulation.h"
#include "allostery.h"
#include "mutation.h"
#include "diagonalize.h"
#include "Fit_B.h"
#include "Residues_distances.h"
#include "Residues_propensity.h"

// #include "rmsd.h"
// #include "rotation.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


// Global variables
double mass_sum=0;


extern float Fit_fluctuations(double *mass_sum, float *B_TNM,
			      float *RMSD_NM, char *out, int anharmonic,
			      struct Normal_Mode NM,
			      struct Reference Ref_kin, char *REF,
			      struct residue *seq1, int nres1, int nmr1,
			      atom *atoms1, char *nameout1,
			      float *factor_B, int *outlier_B);
extern float RMSD_CA(float *coord_1, float *coord_2, struct Reference Ref,
	      atom *atoms, int nres, struct ali_atoms ali_atoms);


/*********************************************************************/
/**********************  Memory handling  ****************************/
void Allocate_memory(struct Normal_Mode *NM,
		     struct Jacobian *J,
		     double ***Hessian,
		     int N_axes, int N_ref, int N_cart,
		     int N_modes, int N_res);
void Clean_memory(struct Normal_Mode NM,
		  struct Jacobian J,
		  double **Hessian,
		  int N_axes, int N_ref,
		  int N_cart, int N_modes,
		  int ANM);

/***************************  Printing  *********************************/

void Print_summary(char *name_out, char *name1, char *parameters, int Nchain,
		   struct Normal_Mode NM, char *out_B1, char *out_B2,
		   float Tors_fluct,struct axe *axes, int N_axes, int N_main,
		   int Nskip, int N_atoms, int N_res, int N_diso, int N_int,
		   int N_inter, float Cont_overlap,
		   int N_disc_coll, int N_disc_freq, int N_disc_out,
		   float mass_sum);
void Print_thermal(int anharmonic, FILE *file_out,
		   struct Normal_Mode NM, char *out_B,
		   struct axe *axes, int N_main, int N_res, int Nchain);


/*****************************  Input  *********************************/
int getArgs(int argc, char **argv,
	    char *file_pdb1, char *chain1,
            char *file_pdb2, char *chain2,
	    int *ANM, char *REF, int *LABEL,
	    int *SIDECHAINS, int *OMEGA, int *PSI,
	    double *K_OMEGA, double *K_PSI,
	    double *K_PHI, double *K_CHI,
	    double *K_BA, double *K_BL,
	    float *E_MIN, float *COLL_THR,
	    int *MIN_INT, char *INT_TYPE, float *thr,
	    int *S_TYPE, float *S_THR,
	    int *ONEINT, int *N_RESRES,
	    int *N_MODE, char *outdir,
	    int *PRINT_CONFCHANGE, int *PRINT_FORCE,
	    int *PRINT_PDB, int *PRINT_MODE_SUMM,
	    char **FILE_FORCE, int *ALLOSTERY,
	    float *KAPPA, int *FIT_B, float *RMSD_EXP, 
	    float *RMSD_MIN, int *NMUT, int *PRED_MUT,
	    char *Mut_para,
	    struct Para_simul *Para_simul,
	    struct Para_confchange *Para_confchange,
	    float *SEQID_THR,
	    char *parameters,
	    int *PRINT_COV_COUPLING,
	    int *PRINT_DEF_COUPLING,
	    int *PRINT_DIR_COUPLING,
	    int *PRINT_COORD_COUPLING,
	    int *PRINT_SIGMA_DIJ,
	    char *PROF_TYPE,
	    int *ALL_PAIRS, float *SIGMA,
	    int *STRAIN, char *SITES,
	    int *PRINT_CMAT, int *ANHARMONIC,
	    int *PRINT_PDB_ANHARM, int *N_PDB_print);

int Read_para(char *filename,
	      char *file_pdb1, char *chain1,
	      char *file_pdb2, char *chain2,
	      int *ANM, char *REF, int *LABEL,
	      int *SIDECHAINS, int *OMEGA, int *PSI,
	      double *K_OMEGA, double *K_PSI,
	      double *K_PHI, double *K_CHI,
	      double *K_BA, double *K_BL,
	      float *E_MIN, float *COLL_THR, int *MIN_INT,
	      char *INT_TYPE, float *thr,
	      int *S_TYPE, float *S_THR,
	      int *ONEINT, int *N_RESRES,
	      int *N_MODE, char *outdir,
	      int *PRINT_CONFCHANGE, int *PRINT_FORCE,
	      int *PRINT_PDB, int *PRINT_MODE_SUMM,
	      char **FILE_FORCE, int *ALLOSTERY,
	      float *KAPPA, int *FIT_B, float *RMSD_EXP,
	      float *RMSD_MIN, int *NMUT, int *PRED_MUT,
	      char *Mut_para,
	      struct Para_simul *Para_simul,
	      struct Para_confchange *Para_confchange,
	      float *SEQID_THR,
	      int *PRINT_COV_COUPLING,
	      int *PRINT_DEF_COUPLING,
	      int *PRINT_DIR_COUPLING,
	      int *PRINT_COORD_COUPLING,
	      int *PRINT_SIGMA_DIJ,
	      char *PROF_TYPE,
	      int *ALL_PAIRS, float *SIGMA,
	      int *STRAIN, char *SITES,
	      int *PRINT_CMAT, int *ANHARMONIC,
	      int *PRINT_PDB_ANHARM, int *N_PDB_print);
void help (void);

/******************************* Computations ****************************/
int Torsional_Hessian(double **Hessian, int N_modes,
		      struct axe *axe,
		      float K_OMEGA, float K_PSI, float K_PHI, float K_CHI,
		      float K_BA, float K_BL);
int Relevant_modes(float *sigma2, int *select, int N, float N_atoms);
double Rescale_force(struct Normal_Mode NM,
		     struct interaction *Int_list, int N_int,
		     float kappa, int anharmonic);

extern void Binding_site_dynamics(struct Normal_Mode NM, struct Reference Ref,
				  struct residue *seq, atom *atoms,
				  char *nameout1,
				  char *pdb, char *chain, int Nchain,
				  int Nres, char *SITES);
extern void Center_atoms(atom *atoms, int N_atoms, struct Reference Ref);
extern void Reorient_axes(atom *atoms, int N_atoms, struct Reference Ref);
extern int Periodic_angles(float *dphi, struct axe *axe, int n);

/****************************** Main code  ********************************/

int main(int argc , char *argv[]){

  
  printf("\n************************************************\n");
  printf("*********** Get arguments/parameters ***********\n");
  printf("************************************************\n\n");
  
  /****************************  Variables  *******************************/
  // Degrees of freedom
  int PSI=1;        // Use psi angle or freeze it?
  int SIDECHAINS=0; // Use side chains as degrees of freedom?
  int OMEGA=0;      // Use also omega as degree of freedom?
                    // 0=No 1=Yes -1=Yes if Domega>D_omega_thr
  D_omega_thr=90;
  int MIN_INT_MAIN=1; // Minimum number of interactions per degree of freedom
  int MIN_INT_SIDE=1; // Minimum number of interactions per degree of freedom

  // Control
  int ANM=0;
  char INT_TYPE[40];
  //char REF[20]; // In interactions_tnm.h
  int N_MODE_PRINT=0;
  int ALLOSTERY=0;

  
  //float S_THR=0.20;
  //int S_TYPE=0; // type of shadow interaction

  // Struct 1
  int ANISOU1=0;
  int nres1=0, nmr1=0, natoms1=0;
  atom *atoms1=NULL;
  char  file_pdb1[150], chain1[CHMAX], pdbid1[100], nameprot1[100];
  struct residue *seq1;
  float *coord_1=NULL;
  struct axe *axe1=NULL;
  int naxe1=0;  // total number of degrees of freedom
  int nmain=0;  // main axes degrees of freedom
  int N_diso1=0;// Number of disordered axis for protein 1
  int nskip=0;  // skipped main chains

  int n_lig1=0, na_lig1=0;
  struct residue *ligres1=NULL;
  atom *ligatom1=NULL;

  // Chains
  int Nchain1=0, Nchain;
  struct chain *chains1=NULL;

  // Struct 2
  int nres2=0, natoms2=0, Nchain2=0;
  int ANISOU2=0, nmr2=0;
  int N_diso2=0; // Number of disordered axis for protein 2
  atom *atoms2=NULL;
  char  file_pdb2[150], chain2[CHMAX], pdbid2[100], nameprot2[100];
  struct residue *seq2;
  float *coord_2=NULL;
  struct chain *chains2=NULL;
  int n_lig2=0, na_lig2=0;
  struct residue *ligres2=NULL;
  atom *ligatom2=NULL;

  // Comparison
  int CONF_CHANGE=0;
  int last_ali_res;
  struct Ali_score ali;
  double rmsd=0.0;
  struct Para_confchange Para_confchange;
  float SEQID_THR=50;

  // Input force
  char *FILE_FORCE=NULL;

  // Simulation
  struct Para_simul Para_simul;

  /* Global (nma_para.h):
     int ALL_AXES;
     float EXP_HESSIAN;
     int PRINT_LAMBDA;
     int PRINT_AXES;
   */

  // GO Model
  double **Hessian=NULL;
  int N_int=0;
  float Cont_overlap=-1;

  // Kinematics
  struct Jacobian J;

  // Normal modes
  struct Normal_Mode NM;

  // Dummies
  // int i,j,k,l,jj, ia, ib, na, nc, a;
  int i,ia,k;

  // time
  double nbtops = CLOCKS_PER_SEC, t0=clock();

  // File names
  #define NCH 150
  char outdir[NCH]="", nameout1[NCH]="", summary1[NCH]=""; //directory[150],
  char fullmodel[NCH]="", MODEL[80]="";
  char name2[NCH]="", nameout2[NCH]="", summary2[NCH]="";
  char parameters[1000]="";

  /*****************************  INPUT  **********************************/
  int PRINT_CONFCHANGE=0;
  int PRINT_FORCE=PRINT_FORCE_DEF;
  int PRINT_PDB= PRINT_PDB_DEF;
  int PRINT_MODE_SUMM=PRINT_MODE_SUMM_DEF;

  outdir[0]='\0';
  ALL_AXES=1;
  KINETIC=1;
  E_MIN=0.00000001; //E_MIN_DEF;
  ONEINT=1;
  N_RESRES=1;
  HYD=0; HYD_REF=0; // Read hydrogen atoms or not?
  KAPPA=KAPPA_DEF;
  FIT_B=1; R_MIN=0; SLOPE_MIN=0.00; RMSD_EXP=0.5;
  K_OMEGA=K_OMEGA_DEF; K_PSI=K_PSI_DEF; K_PHI=K_PHI_DEF;
  K_CHI=K_CHI_DEF; K_BA=K_BA_DEF; K_BL=K_BL_DEF;
	 
  getArgs(argc,argv,file_pdb1,chain1,file_pdb2, chain2,
	  &ANM, REF, &LABEL, &SIDECHAINS, &OMEGA, &PSI,
	  &K_OMEGA, &K_PSI, &K_PHI, &K_CHI, &K_BA, &K_BL,
	  &E_MIN, &COLL_THR, &MIN_INT_SIDE, INT_TYPE, &C_THR, &S_TYPE,
	  &S_THR, &ONEINT, &N_RESRES, &N_MODE_PRINT, outdir,
	  &PRINT_CONFCHANGE, &PRINT_FORCE, &PRINT_PDB, &PRINT_MODE_SUMM,
	  &FILE_FORCE, &ALLOSTERY, &KAPPA, &FIT_B, &RMSD_EXP, &RMSD_MIN,
	  &NMUT, &PRED_MUT, Mut_para,
	  &Para_simul, &Para_confchange, &SEQID_THR, parameters,
	  &PRINT_COV_COUPLING, &PRINT_DEF_COUPLING,
	  &PRINT_DIR_COUPLING, &PRINT_COORD_COUPLING,
	  &PRINT_SIGMA_DIJ, &PROF_TYPE, &ALL_PAIRS, &SIGMA, &STRAIN, SITES,
	  &PRINT_CMAT, &ANHARMONIC, &PRINT_PDB_ANHARM, &N_PDB_print);
  MIN_INT_MAIN=MIN_INT_SIDE;
  // Rescale torsion springs
  K_PSI*=KAPPA; K_PHI*=KAPPA; K_OMEGA*=KAPPA; K_CHI*=KAPPA;
  K_BA*=KAPPA; K_BL*=KAPPA;

  // Control parameters
  if(HNM){
    strcpy(MODEL, "HNM");
    strcpy(REF, "CA");
    strcpy(INT_TYPE, "CA");
    printf("WARNING, Hinsen ANM (HNM) chosen, ");
    printf("setting interaction types and ref. atoms to CA\n");
    ANM=1;
  }
  if(ANM){
    char REFT[4]="CA";
    strcpy(MODEL, "ANM");
    if((strncmp(INT_TYPE, "MIN", 3)==0)||
       (strncmp(INT_TYPE, "ALL", 3)==0)){
      strcpy(INT_TYPE, "SCA");
      printf("WARNING, all atoms contacts but interacting atoms set to CA.");
    }else if((strncmp(INT_TYPE, "CB", 2)==0)||
	     (strncmp(INT_TYPE, "SCB", 3)==0)){
      strcpy(REFT, "CB");
    }else if(strncmp(INT_TYPE, "CA", 2)!=0){
      printf("WARNING, interaction type %s not allowed with ANM\n", INT_TYPE);
      strcpy(INT_TYPE, "CA");
    }
    if(strcmp(REF, REFT)!=0){
      printf("WARNING, %s not allowed reference atoms with ANM.", REF);
      printf("Reference atoms set to %s\n", REFT);
      strcpy(REF, REFT);
    }
    printf("ANM (Cartesian d.o.f.) ");
    //E_MIN*=10;
  }else{
    // TNM
    strcpy(MODEL, "TNM");
    if(REF[0]=='\0')strcpy(REF, REF_DEF);
    char REFT[4]; strcpy(REFT, REF);
    if((strncmp(INT_TYPE, "CA", 2)==0)||
       (strncmp(INT_TYPE, "SCA", 3)==0)){
	strcpy(REFT, "CA");
    }else if((strncmp(INT_TYPE, "CB", 2)==0)||
	     (strncmp(INT_TYPE, "SCB", 3)==0)){
	strcpy(REFT, "CB");
    }else if(strncmp(INT_TYPE, "HYD", 3)==0){
      HYD=1; HYD_REF=0;
    }
    if((strcmp(REF,REFT)!=0)&&(strcmp(REF,"ALL")!=0)&&(strcmp(REF,"EB")!=0)){
      printf("WARNING, %s not allowed reference atoms with %s contacts",
	     REF, INT_TYPE);
      printf(" Setting reference atoms to %s\n", REFT);
      strcpy(REF, REFT);
    }
  }
  printf("Reference atoms: %s Contact type: %s\n", REF, INT_TYPE);
  NM.ANM=ANM;

  AA_code=AA_BKV;
  Econt=Allocate_mat2_f(20, 20);
  for(i=0; i<20; i++)for(k=0; k<20; k++)Econt[i][k]=E_BKV[i][k];
  //for(i=0; i<20; i++)for(k=0; k<20; k++)Econt[i][k]=U_Residues[i][k];
  for(i=0; i<20; i++)for(k=0; k<20; k++)U_Res[i][k]=U_Residues[i][k];
  for(i=0; i<20; i++)for(k=0; k<20; k++)D_Res[i][k]=D_Residues[i][k];


  printf("\n************************************************\n");
  printf("*********** Read protein structure 1 ***********\n");
  printf("************************************************\n\n");


  /********************  Protein structure 1  ***********************/
  /* Read protein structure 1*/
  Read_PDB(&nres1, &seq1, &n_lig1, &ligres1, chain1, &ANISOU1, &nmr1,
	   &natoms1, &atoms1, &na_lig1, &ligatom1, &chains1, &Nchain1,
	   pdbid1, file_pdb1);
  Nchain=Nchain1;
  // Collectivity
  //if(COLL_THR > natoms1*0.15)COLL_THR = natoms1*0.15;
  COLL_THR *=3;
  //N_TNM=naxe1;
  

  /*******************  Protein structure 2, if any  ********************/
  /* Read protein structure 2 */
  int Nmut=0, *Posmut=NULL; char *AAmut=NULL, *AAwt=NULL;
  if(file_pdb2[0]!='\0'){
  
    printf("\n************************************************\n");
    printf("*********** Read/align structure 2 *************\n");
    printf("************************************************\n\n");
    
    Read_PDB(&nres2, &seq2, &n_lig2, &ligres2, chain2, &ANISOU2,
	     &nmr2, &natoms2, &atoms2, &na_lig2, &ligatom2,
	     &chains2, &Nchain2, pdbid2, file_pdb2);
    if(nres2<=0){
      printf("WARNING, no residues found in file %s\n", file_pdb2);
      printf("Performing single structure computation\n");
      for(i=0;i<natoms1;i++)atoms1[i].ali=i;
      goto end_ali;
    }
    CONF_CHANGE=1;
    MIN_ALI_LEN=MIN_ALI_LEN_DEF;
    Posmut=malloc(nres2*sizeof(int));
    AAmut=malloc(nres2*sizeof(char));
    AAwt=malloc(nres2*sizeof(char));
    Align_chains(&Nmut, AAwt, Posmut, AAmut, &ali, &last_ali_res,
		 &Nchain, SEQID_THR,
		 chains1, Nchain1, seq1, pdbid1,
		 chains2, Nchain2, seq2, pdbid2);
    printf("%d mutations found: ", Nmut);
    for(i=0; i<Nmut; i++)printf(" %c%d%c", AAwt[i], Posmut[i], AAmut[i]);
    printf("\n");
    if((NMUT>=0)&&(Nmut!=NMUT)){
      printf("WARNING, %d mutations found, %d expected. ", Nmut, NMUT);
      printf("If you want that the analysis is run, please ");
      printf("set the parameter NMUT=-1\nExiting the program\n");
      exit(8);
    }
    int N_ali_atom=0;
    for(i=0; i<Nchain1; i++){
      struct chain *ch1= chains1+i, *ch2= chains2+ch1->match;
      if(i<Nchain){
	N_ali_atom+=Align_atoms(atoms1, atoms2, ch1, ch2);
      }else{
	for(int j=ch1->ini_atom; j<(ch1->ini_atom+ch1->natoms); j++){
	  atoms1[j].ali=-1;
	}
      }
    }
    printf("%d aligned atoms out of %d and %d in %d chains\n",
	   N_ali_atom, natoms1, natoms2, Nchain);

    // Eliminate not aligned atoms before setting bonds and reference atoms
    // Also maintain not aligned atoms that belong to aligned residues
    //if(SWITCH_BONDS){
    Purge_not_aligned(chains1, atoms1, &natoms1,
		      chains2, atoms2, &natoms2);
    printf("Prot 1: %d Prot 2: %d\n", natoms1, natoms2);
    // From this point on, atoms1[i].ali=i

  }else{
    for(i=0;i<natoms1;i++)atoms1[i].ali=i;
  }
 
 end_ali:
  if(CONF_CHANGE)printf("%d atoms found in structure 2\n", natoms2);


  printf("\n************************************************\n");
  printf("*********** Ref atoms Kinetic Energy ***********\n");
  printf("************************************************\n\n");

  /*************  Reference atoms for kinetic energy  ***********************/
  struct Reference Ref_kin;
  int N_ref=Set_reference(&Ref_kin, 0, REF, atoms1, 0, natoms1);
  int N_cart=3*N_ref;
  // Check that reference atoms exist
  if((natoms1==0)||(N_ref==0)){
    printf("ERROR, no atoms found in file %s", file_pdb1);
    printf(" natoms= %d N_ref= %d\n", natoms1, N_ref);
    exit(8);
  }
  printf("%d reference atoms\n", N_ref);

  printf("\n************************************************\n");
  printf("*********** Center coordinates *****************\n");
  printf("************************************************\n\n");
  /*******************************************************************/
  // Change the coordinate system
  // Put center of mass at the origin
  Center_atoms(atoms1, natoms1, Ref_kin);
  // Reorient the Cartesian axes along the principal axes
  Reorient_axes(atoms1, natoms1, Ref_kin);

  printf("\n************************************************\n");
  printf("*********** Prepare output files ***************\n");
  printf("************************************************\n\n");

  /************************  Output files ***************************/
  sprintf(nameprot1, "%s", pdbid1);
  for(i=0; i<Nchain; i++){
    if((chains1[i].label!=' ')&&(chains1[i].label!='\0'))
      sprintf(nameprot1, "%s%c", nameprot1, chains1[i].label);
  }
  strcpy(nameout1, nameprot1);

  if(strncmp(INT_TYPE, "SCR", 3)==0){
    if(ONEINT==0){
      sprintf(fullmodel, "%s_d%.2f_t%.1f", INT_TYPE, S_THR, C_THR);
    }else{
      sprintf(fullmodel, "%s_MIN_d%.2f_t%.1f", INT_TYPE, S_THR, C_THR);
    }
  }else if(strncmp(INT_TYPE, "SHA", 3)==0){
    if(ONEINT==0){
      sprintf(fullmodel, "%s_%d_d%.2f_t%.1f", INT_TYPE, S_TYPE, S_THR, C_THR);
    }else{
      sprintf(fullmodel, "%s_MIN_%d_d%.2f_t%.1f",
	      INT_TYPE, S_TYPE, S_THR, C_THR);
    }
  }else if((strncmp(INT_TYPE, "MIN", 3)==0)&&(N_RESRES>1)){
    sprintf(fullmodel, "%s%.1f_M%d", INT_TYPE, C_THR, N_RESRES);
  }else{
    sprintf(fullmodel, "%s%.1f", INT_TYPE, C_THR);
  }
  sprintf(fullmodel, "%s_%s", fullmodel, REF);
  sprintf(fullmodel, "%s_PHI", fullmodel);
  if(PSI)sprintf(fullmodel, "%sPSI", fullmodel);
  if(OMEGA)sprintf(fullmodel, "%sOME", fullmodel);
  if(OMEGA<0)sprintf(fullmodel, "%sCISTR", fullmodel);
  if(SIDECHAINS)sprintf(fullmodel, "%sSCH", fullmodel);
  if(LABEL)sprintf(nameout1,  "%s_%s", nameout1, fullmodel);
  // UUU: Restore the option of the full name

  sprintf(summary1, "%s.summary.dat", nameout1);   
  if(Check_make_dir(outdir)){
    sprintf(summary1, "%s/%s", outdir, summary1);
    sprintf(nameout1,  "%s/%s", outdir, nameout1);
  }

  printf("output file (nameout1): %s\n",nameout1);
  printf("output file (summary1): %s\n",summary1);

  if(CONF_CHANGE){
    sprintf(nameprot2, "%s", pdbid2);
    for(i=0; i<Nchain; i++){
      if((chains2[i].label!=' ')&&(chains2[i].label!='\0'))
	sprintf(nameprot2, "%s%c", nameprot2, chains2[i].label);
    }
    sprintf(name2, "%s_%s", nameprot1, nameprot2);
    //sprintf(name2, "%s-%s", nameprot1, nameprot2);
    //sprintf(name2, "%s_%s", pdbid1, pdbid2); // YYY don't write chain names 


    sprintf(nameout2, "%s_%s", nameprot1, nameprot2); 
    if(LABEL)sprintf(nameout2, "%s_%s", nameout2, fullmodel);

    sprintf(summary2, "%s.summary.dat", nameout2);
    if(Check_make_dir(outdir)){
      sprintf(summary2, "%s/%s", outdir, summary2);
      sprintf(nameout2, "%s/%s", outdir, nameout2);
    }
    printf("output file (nameout2): %s\n",nameout2);
    printf("output file (summary2): %s\n",summary2);
  }
  
  printf("\n************************************************\n");
  printf("*********** Topology? ***************\n");
  printf("************************************************\n\n");

  // bonds are made only for atoms with ali=1!
  struct bond *bonds=Set_bonds_topology(natoms1, atoms1, seq1);
  // Test bonds
  int naxe=3*nres1;
  Test_buildup(bonds, atoms1, natoms1, naxe);
  // New structure modified by build-up
  for(i=0; i<natoms1; i++)for(int j=0; j<3; j++)atoms1[i].r[j]=bonds[i].r[j];

  struct Reference Ref1;
  struct ali_atoms ali_atoms;
  if(CONF_CHANGE){
    // Check RMSD of conformational change
    // Set reference atoms only if aligned
    Set_reference(&Ref1, 0, REF, atoms1, 0, natoms1); //"EB"
    int N_ali=Align_references(&ali_atoms, Ref1, atoms1);
    printf("%d ref. atoms for conf. change, aligned: %d\n",Ref1.N_ref, N_ali);
    // Compute RMSD
    int N_cart=ali_atoms.N_cart;
    coord_1=malloc(N_cart*sizeof(float));
    Write_ref_coord_atom(coord_1, N_ali, atoms1, ali_atoms.ali1);
    coord_2=malloc(N_cart*sizeof(float));
    Write_ref_coord_atom(coord_2, N_ali, atoms2, ali_atoms.ali2);
    float rmsd=RMSD_CA(coord_1, coord_2, Ref1, atoms1, nres1, ali_atoms);
    printf("RMSD_CA(%s,%s)=%.2f\n", pdbid1, pdbid2, rmsd);
    if((rmsd < RMSD_MIN)||(rmsd > RMSD_MAX)){
      printf("Too small RMSD, exiting the program ");
      printf("allowed: %.2f-%.2f\n", RMSD_MIN, RMSD_MAX);
      exit(8);
    }
    //rmsd=rmsd_mclachlan_f(coord_1, coord_2, ali_atoms.mass, N_ali);
    if((OMEGA<0)){ //SWITCH_BONDS || 
      // bond->d_phi is computed here
      Switch_bonds(bonds, naxe, name2, ali_atoms,
		   atoms1, natoms1, seq1,
		   atoms2, natoms2, seq2);
    }
    
  }else{ // No conformation change
    last_ali_res=nres1;
    for(i=0; i<natoms1; i++)atoms1[i].ali=i;
  }
 
  printf("\n************************************************\n");
  printf("*********** ENM contacts/interactions **********\n");
  printf("************************************************\n\n");

  /*******************************************************************
                          ENM Computations
  ********************************************************************/
  // Consider covalent contacts in ANM but not in TNM
  NOCOV=1; if(ANM)NOCOV=0;

  /*************************  Interactions  *************************/
  struct interaction *Int_list; 
  Compute_interactions(&N_int, &Int_list, INT_TYPE,
		       atoms1, natoms1, nres1, nameprot1);
  if(PRINT_CMAT){
    printf("Printing contact matrix\n");
    Print_contact_matrix(Int_list,N_int,atoms1,seq1,fullmodel,nameprot1);
  }
  printf("%d interactions\n", N_int);
  int N_inter=0;  // Count interchain interactions
  if(Nchain > 1){
    for(k=0; k<N_int; k++){
      if(atoms1[Int_list[k].i1].chain!=atoms1[Int_list[k].i2].chain)
	N_inter++;
    }
    printf("%d Interchain contacts\n", N_inter);
  }

  if((strncmp(INT_TYPE, "SCR", 3)==0)||(strncmp(INT_TYPE, "SHA", 3)==0)){ 
      // Screened interactions
    int alignres[nres1], N_int3;
    for(i=0; i<nres1; i++)alignres[i]=i;
    struct interaction *Int_list3;
    Compute_interactions(&N_int3, &Int_list3, "MIN",
			 atoms1, natoms1, nres1, nameprot1);
    Cont_overlap=Contact_overlap(Int_list3,N_int3,Int_list,N_int,alignres);
    printf("MIN: %d interactions q=%.3f (%s-MIN)\n",
	   N_int3, Cont_overlap, INT_TYPE);
    free(Int_list3);
  }
  printf("Contact matrix computation.Time= %.2lf sec.\n",(clock()-t0)/nbtops);
  t0=clock();


  printf("\n************************************************\n");
  printf("*********** ENM degrees of freedom *************\n");
  printf("************************************************\n\n");
  
  /***********************  Degrees of freedom ************************/
  axe1=Set_DegofFreed(&naxe1, &nmain, &nskip, &N_diso1, bonds,
		      atoms1, natoms1, Ref_kin.atom_num, N_ref,
		      seq1, last_ali_res, chains1, Nchain,
		      Int_list, N_int, MIN_INT_MAIN, MIN_INT_SIDE,
		      OMEGA, SIDECHAINS, PSI);
  printf("Number of degrees of freedom: %d\n", naxe1);
  if(naxe1!=nmain)printf("Side chains dof: %d\n", naxe1-nmain);

  // Standardize bond lengths and bond angles
  int STANDARDIZE=0;
  if(STANDARDIZE){
    Standardize_bonds(bonds, atoms1, natoms1, nameprot1,
		      axe1, naxe1, seq1, nres1, chains1,
		      Nchain, 0.1, Para_simul);
  }

   printf("\n************************************************\n");
  printf("*********** ENM memory allocation **************\n");
  printf("************************************************\n\n");

  /************************ Prepare memory ***************************/
  int N_modes;
  if(ANM){N_modes=N_cart;}else{N_modes=naxe1;}
  Allocate_memory(&NM, &J, &Hessian, naxe1, N_ref, N_cart, N_modes, nres1);
  printf("Memory allocated\n");


   printf("\n************************************************\n");
  printf("*********** ENM kinetic en, Eckart cdt *********\n");
  printf("************************************************\n\n");


  /********* Eckart conditions, Jacobian, Kinetic energy  ************/
  J.N_kin=
    Compute_kinetic(&J, axe1, naxe1, atoms1, natoms1, Ref_kin, 1);
  if((ANM==0)&&(KINETIC)){N_modes=J.N_kin; NM.N=N_modes;}
  printf("Kinetic energy computation. Time= %.2lf sec.\n",(clock()-t0)/nbtops);
  t0=clock();

  /************************** Debugging *****************************/
  if(DEBUG){
    for(i=0; i<Nchain; i++){
      struct chain *chp=chains1+i;
      printf("### Chain %c\n", chp->label);
      printf("%3d axes: %d - %d\n", chp->mainaxes,
	     chp->ini_main, chp->ini_main+chp->mainaxes-1);
      printf("%4d atoms: %d - %d\n", chp->natoms,chp->ini_atom,
	     chp->natoms+chp->ini_atom-1);
      printf("%3d residues: %d - %d\n", chp->nres,chp->ini_res,
	     chp->nres+chp->ini_res-1);
    }
  }


   printf("\n************************************************\n");
  printf("*********** ENM Hessian / normal modes  ********\n");
  printf("************************************************\n\n");

  /*************************  Normal modes TNM *************************/
  if(ANM==0){ // TNM
    printf("Torsional network model\n");
    if(KINETIC==0)printf("WARNING, kinetic energy not considered\n\n");
    Compute_Hessian_TNM(Hessian, J.T_sqrt, J.T_sqrt_inv, Int_list, N_int,
			atoms1, natoms1, axe1, nmain, naxe1, N_modes,
			chains1, Nchain, KINETIC,
			K_OMEGA, K_PSI, K_PHI, K_CHI, K_BA, K_BL);

    printf("Hessian computation. Time= %.2lf sec.\n",(clock()-t0)/nbtops);
    t0=clock();

    d_Diagonalize(N_modes, Hessian, NM.omega2, NM.MW_Tors, -1);
    printf("Hessian diagonalization. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
    t0=clock();

    // Allocate only relevant modes
    //NM.N_relevant=Relevant_modes(NM.sigma2, NM.select, NM.N, NM.N_cart/40);
    NM.N_relevant=NM.N;
    for(ia=0; ia<NM.N; ia++){
      Transform_tors_modes(NM.Tors[ia], NM.MW_Tors[ia],
			   J.T_sqrt_tr, J.T_sqrt_inv_tr,
			   N_modes, naxe1);
      Convert_torsion2cart(NM.Cart[ia], atoms1, NM.Tors[ia],
			   axe1, naxe1, Ref_kin, ia);
      /* Outdated
	 Convert_torsion2cart_old(NM.Cart[ia], atoms1, NM.Tors[ia], axe1,
			       naxe1, Ref1.atom_num, N_ref, &J);
	 Normalize_vector_weighted(NM.Cart[ia], Ref1.mass_coord, N_cart);
      */
    }

  /*************************  Normal modes ANM *************************/
  }else{ // ANM
    if(HNM){printf("Hinsen network model");}
    else{printf("Anisotropic network model\n");}

    Compute_Hessian_ANM(Hessian, N_ref, Ref_kin.atom_num,
			Int_list, N_int, atoms1, natoms1);
    printf("Hessian computation. Time= %.2lf sec.\n",(clock()-t0)/nbtops);

    d_Diagonalize(N_modes, Hessian, NM.omega2, NM.Cart, -1);
    printf("Hessian diagonalization. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
    t0=clock();
    
    //  Convert normal modes
    struct Tors *D_NM=malloc(N_modes*sizeof(struct Tors));
    for(ia=0; ia<NM.N; ia++){
      //if(ia<10)printf("om2= %.2g\n", NM.omega2[ia]);
      Normalize_vector_weighted(NM.Cart[ia], Ref_kin.mass_coord, N_cart);
      D_NM[ia].Cart=NM.Cart[ia];
      D_NM[ia].Tors=NM.Tors[ia];
      D_NM[ia].MW_Tors=NM.MW_Tors[ia];
      D_NM[ia].coeff=malloc(1*sizeof(float));
      D_NM[ia].N_axes=naxe1;
      D_NM[ia].N_cart=N_cart;
      Convert_cart2torsion(&(D_NM[ia]), Ref_kin, &J);
      NM.Tors_frac[ia]=Tors_fraction(&(D_NM[ia]), Ref_kin.mass_coord);
    }
  }
  printf("Normal mode computation. Time= %.2lf sec.\n",
	 (clock()-t0)/nbtops);
  t0=clock();



  printf("\n************************************************\n");
  printf("*********** ENM normal mode coll. **************\n");
  printf("************************************************\n\n");



  /********************  Collectivity of normal modes *******************/
  // Select modes based on eigenvalues and collectivity
  int N_disc_freq=0, N_disc_coll=0;

  // Select modes based on collectivity
  float inv_sq_mass[N_cart];
  for(i=0; i<N_cart; i++)inv_sq_mass[i]=1./sqrt(Ref_kin.mass_coord[i]);

  double E_ave=0; for(ia=0; ia<NM.N; ia++)E_ave+=NM.omega2[ia]; E_ave/=NM.N;
  printf("Average eigenvalue omega^2 prior to norm: %.2g\n", E_ave);
  double E_THR=E_MIN*E_ave;
  
  if(1) // YYY additional output in log file (related to discarded modes)
  	{
  	 printf("(discard) AV_OME2     %10.4g   | average EV (or omega^2) before kappa fit\n", E_ave);
 	 double E_logave=0.0;
  	 double E_xmin=pow(10,10);
  	 double E_xmax=(-1)*pow(10,10);
  	 double E_normsig2=0.0;
  	 for(ia=0; ia<NM.N; ia++)
  		{E_logave+=log(NM.omega2[ia]);
  	 	 E_normsig2+=1.0/(NM.omega2[ia]);
  	 	 if(NM.omega2[ia]>E_xmax)E_xmax=NM.omega2[ia];
  	 	 if(NM.omega2[ia]<E_xmin)E_xmin=NM.omega2[ia];
  		}
  	 E_logave=exp(E_logave/NM.N);
  	 printf("(discard) LOGAV_OME2  %10.4g   | exp(average log(omega^2))\n",
		E_logave);
 	 printf("(discard) MIN_OME2    %10.4g   | min omega^2\n", E_xmin);
  	 printf("(discard) MAX_OME2    %10.4g   | max omega^2\n", E_xmax);
  	 printf("(discard) NORM_SIG2   %10.4g   | ", E_normsig2);
	 printf("sum(1/omega^2) --> CONTR_TH = (1/omega^2)/NORM_SIG2\n");
  	 printf("(discard) THR_OME2    %10.4g   | = EMIN * AV_OME2\n", E_THR);
  	 printf("(discard) THR_COLL    %10.4g   |", COLL_THR);
	 printf("= COLL * 3 (cart coord per atom), max= 0.15*3*natoms\n");
  	 printf("(discard) N_cart      %10d   | ", N_cart);
	 printf("= natoms * 3 --> collect = xxx/N_cart\n");
  	}
  	
  for(ia=0; ia<NM.N; ia++){

    NM.Max_dev[ia]=
      Compute_Max_dev(NM.Cart[ia], NM.omega2[ia], N_cart, inv_sq_mass);
    NM.Tors_coll[ia]=Collectivity_norm2(NM.Tors[ia], naxe1)/naxe1;
    NM.MW_Tors_coll[ia]=Collectivity_norm2(NM.MW_Tors[ia],naxe1)/naxe1;
    NM.Cart_coll[ia]=Collectivity_norm2(NM.Cart[ia], N_cart); //_Renyi
    if((NM.omega2[ia]<=E_THR)||(NM.Cart_coll[ia] < COLL_THR)){
      printf("WARNING, discarding mode %d ", ia);
      printf(" omega^2= %.3g Coll= %.0f RMSD= %.2g\n",
	     NM.omega2[ia], NM.Cart_coll[ia], NM.Max_dev[ia]);
      if(NM.Cart_coll[ia] < COLL_THR){N_disc_coll++;}
      else{N_disc_freq++;}
      NM.select[ia]=0; NM.sigma2[ia]=0;
    }else{
      NM.select[ia]=1;
      NM.sigma2[ia]=1./NM.omega2[ia];
    }
    NM.Cart_coll[ia]/=N_cart;
  }
  printf("%d modes selected\n", NM.N-(N_disc_coll+N_disc_freq));
  printf("%d modes discarded for low collectivity ", N_disc_coll);
  printf("and %d for low frequency\n", N_disc_freq);
  int amin=0; while(NM.omega2[amin]<0)amin++;
  printf("Smallest positive eigenvalue: %.2g a=%d ",NM.omega2[amin],amin);
  printf("Threshold for selection: %.2g*<omega^2>\n", E_MIN);

  // Ugo_mod 05/21
  int N_disc_out=0;
  /* float thr_outlier=2.5, thr_mode=0.5;
     printf("Remove modes that move > %.2f outlier dof (rmsf > %.1f)\n",
     thr_mode, thr_outlier);
     N_disc_out=Filter_modes_outliers(NM, Ref_kin, thr_outlier, thr_mode);
     printf("%d modes removed because they mostly move outliers\n",
     N_disc_out);*/
  printf("Normal mode selection. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
  t0=clock();

  /************************ Debugging *****************************/
  if(DEBUG &&(ANM==0)){
    int ib;
    float **ev=NM.MW_Tors;
    printf("Scalar product between mass-weighted TNM modes\n");
    for(ia=0; ia<20; ia++){
      for(ib=0; ib<=ia; ib++)
	printf(" %.3f", Scalar_product(ev[ia],ev[ib],naxe1));
      printf("\n");
    }
    printf("Scalar product between mass weighted Cartesian TNM modes\n");
    for(ia=0; ia<20; ia++){
      if(NM.select[ia]==0)continue;
      if(ia>=NM.N_relevant)continue;
      for(ib=0; ib<=ia; ib++){
	float *Ca=NM.Cart[ia], *Cb=NM.Cart[ib]; double q=0;
	for(i=0; i<N_cart; i++)q+=Ca[i]*Cb[i]*Ref_kin.mass_coord[i];
	printf(" %.3f", q);
      }
      printf("\n");
    }
  }



  printf("\n**************************************************\n");
  printf("**** Rescale force constant with B-factors ********\n");
  printf("***************************************************\n\n");
  
  char out_B1[1000], out_B2[1000]; int outlier_B[nres1];
  int anharmonic=0;
  float *B_TNM[2]; B_TNM[0]=malloc(N_ref*sizeof(float)); B_TNM[1]=NULL;
  float RMSD_NM, factor_B, kappa=
    Fit_fluctuations(&mass_sum, B_TNM[0], &RMSD_NM, out_B1, anharmonic,
		     NM, Ref_kin, REF, seq1, nres1, nmr1,
		     atoms1, nameout1, &factor_B, outlier_B);

  // Rescale force constant
  double sum_sigma2=Rescale_force(NM, Int_list, N_int, kappa, anharmonic);

  printf("Fitting B factors. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
  printf("(discard) KAPPA   %10.4g    | ratio between OME2 above and in _Modes file\n",kappa); // YYY added output in log file
  t0=clock();


  printf("************************************************\n");
  printf("*************** Anharmonicity ******************\n");
  printf("************************************************\n\n");

  int NA=20; // amino acids or atom types
  int aa_seq[nres1]; for(i=0; i<nres1; i++)aa_seq[i]=seq1[i].i_aa;
  struct interaction *Int_KB[NA];
  for(i=0; i<NA; i++)Int_KB[i]=malloc(NA*sizeof(struct interaction));
  Assign_interactions_KB(Int_KB, NA, POW, EXP_HESSIAN, KAPPA);
  Assign_interactions_SB(Int_list, N_int, Int_KB, aa_seq);

  int n_har=1; 
  if(ANHARMONIC){
    // New
    n_har=2;
    // If atom types, you have to input the atom type instead of aa_seq!
    Print_interactions(Int_list, N_int, Int_KB, aa_seq, AA_code,NA,nameout1);
    double sum=0.0, sum_thr=0.90*sum_sigma2; //0.95
    printf("Anharmonicity corrections\n");
    int MIN_ANH=10, k=0;
    for(i=0; i<NM.N; i++){
      NM.sigma2_anhar[i]=NM.sigma2[i];
      if((sum<sum_thr)||(k<MIN_ANH)){
	if(NM.sigma2[i]>0){
	  Anharmonicity(NM, i, atoms1, natoms1, axe1, naxe1, bonds,
			seq1, nres1, Int_list, N_int, Int_KB, 20, nameout1);
	  k++;
	}
	sum+=NM.sigma2[i];
	if((i<=10)||(i/10)*10==i)
	  printf("%d modes sum= %.3g\n", i, sum/sum_sigma2);
      }else{
	printf("%d modes over %d sum= %.3g exiting\n",
	       i, NM.N, sum/sum_sigma2); break;
      }
    }
    //printf("----> %4d modes, %4d anha | sum_sigma2 = %10.5f  sum_ok = %10.5f  sum_anh = %10.5f\n",NM.N,ix,sum_sigma2,sum,sum_anh);
    // YYY correct the last modes
    int ix; double sum_anh=0; sum=0; 
    for(ix=i-1; ix>=i-10; ix--){
      if(ix<0)break; sum+=NM.sigma2[ix]; sum_anh+=NM.sigma2_anhar[ix];
    }
    if(sum>0){
      sum_anh/=sum;
      for(ix=i; ix<NM.N; ix++)NM.sigma2_anhar[ix]=NM.sigma2[ix]*sum_anh;
    }
    anharmonic=1;
    B_TNM[1]=malloc(N_ref*sizeof(float));
    kappa=Fit_fluctuations(&mass_sum, B_TNM[1], &RMSD_NM, out_B2, anharmonic,
			   NM, Ref_kin, REF, seq1, nres1, nmr1,
			   atoms1, nameout1, &factor_B, outlier_B);

    // Rescale force constant
    //double sum_sigma2_anhar=
    Rescale_force(NM, Int_list, N_int, kappa, anharmonic);
    printf("Anharmonicity analysis. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
    t0=clock();
  }

  printf("************************************************\n");
  printf("*********** Torsional fluctuations *************\n");
  printf("************************************************\n\n");

 /*********************  RMSF of torsion angles ********************/
  Predict_fluctuations(&NM, NM.sigma2);
  double Tors_fluct=0;
  for(i=0; i<NM.N_axes; i++){float x=NM.tors_fluct[i]; Tors_fluct+=x*x;}
  Tors_fluct=sqrt(Tors_fluct/NM.N_axes);
  int outlier_tors[naxe1];
  Tors_outliers(outlier_tors, naxe1, outlier_B, nres1, axe1, chains1, Nchain);

  /************************ Print_summary *****************************/

  float Mass_sqrt=sqrt(mass_sum);
  char *summary, *nameprot;
  if(CONF_CHANGE){
    summary=summary2; nameprot=name2;
  }else{
    summary=summary1; nameprot=nameprot1;
  }
  // Single protein properties
  // Mean_Cart_coll and Coll_thr_cc are computed here
  Print_summary(summary, nameprot, parameters, Nchain, NM,
		out_B1, out_B2, Tors_fluct, axe1, naxe1, nmain, nskip, natoms1,
		nres1, N_diso1, N_int, N_inter, Cont_overlap,
		N_disc_coll, N_disc_freq, N_disc_out, mass_sum);
  //if(CONF_CHANGE==0){
    Print_mode_summary(nameout1,"Modes",NM,Mass_sqrt,anharmonic,kappa);
  //}
  
  printf("************************************************\n");
  printf("************ Conformational Change *************\n");
  printf("************************************************\n\n");

  /*********************  Conformation change  **********************/
  char name_cc[400];
  float *Confchange=NULL;

  if(CONF_CHANGE){

    strcpy(name_cc,name2); // YYY change output file name
    //sprintf(name_cc, "%s_%s", name2, REF);
    //if(OMEGA>0){sprintf(name_cc, "%s_OMEGA", name_cc);}
    //else if(OMEGA<0){sprintf(name_cc, "%s_CISTRANS", name_cc);}
    //if(PSI==0)sprintf(name_cc, "%s_NOPSI", name_cc);

    // Store conformation change over kinetic energy atoms for allostery
    Confchange=malloc(Ref1.N_cart*sizeof(float));
    rmsd=rmsd_mclachlan_f(coord_1, coord_2, Ref1.mass_atom, Ref1.N_ref);
    for(i=0; i<N_ref; i++)Confchange[i]=coord_2[i]-coord_1[i];
    
    struct Tors Diff;
    Allocate_tors(&Diff, naxe1, Ref1.N_cart, N_modes);

    // int n_har=1, anharm; if(ANHARMONIC)n_har=2;
    for(int anharm=0; anharm<n_har; anharm++){
      rmsd=
	Examine_confchange(&Diff, bonds, axe1, chains1, Nchain, summary2,
			   nameout2, name_cc,
			   // YYY change output file name: nameprot2 -> nameout2
			   atoms1, coord_1, seq1, nres1, N_diso1, natoms1, 
			   atoms2, coord_2, seq2, nres2, N_diso2, natoms2, 
			   Ref1, ali_atoms, Int_list, N_int, INT_TYPE, s0,
			   N_modes, NM, 1, ali, // J_CC, 
			   Para_simul, //diff_phi, greedy_phi,
			   0, B_TNM[anharm], outlier_tors, nameout2,
			   PRINT_CONFCHANGE, anharm);
      if(anharm==(n_har-1)){
	Print_modes_confchange(name2,nameout2,NM,Diff,ANM,
			       Mass_sqrt,rmsd, anharm, kappa);
      }
    }

    //Print_diff_fluct(&Diff, &NM, outlier_tors, nameout2);
    
    printf("Writing %s\n", summary);
    printf("Output conf.change. Time= %.2lf sec.\n\n",
	   (clock()-t0)/nbtops); t0=clock();
    
    if((NMUT!=-2)&&(NMUT!=0)){
      // Analysis of mutation
      float mut_phi[naxe1], mut_CC[N_cart];
      for(int anharm=0; anharm<n_har; anharm++){
	int mut_ok=Mutation(mut_phi, naxe1, mut_CC, // output
			    Ref_kin, Nmut, AAwt, Posmut, AAmut,
			    NM, atoms1, natoms1, Int_list, N_int, seq1,
			    B_TNM[anharm], Diff.Cart, Ref1, Diff.Tors,
			    NM.tors_fluct, nameout2, anharm, Mut_para);
	if(0 && mut_ok){
	  Simulate_confchange(30, mut_phi, mut_CC, 1, RMSD_NM, coord_1, coord_2,
			      rmsd, Para_simul, atoms1, natoms1, axe1, naxe1,
			      bonds, seq1, nres1, N_diso1, NM, J, //*J_CC
			      Ref1, Int_list, N_int, INT_TYPE, s0, nameout2,
			      Tors_fluct, anharm);
	}
      }
      printf("Analysis of mutation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
    }

    /* Computing and printing force producing confchange
       The force is computed for all atoms of kynetic energy (Ref_kin),
       not only those representative of conformation change (Ref1) */
    struct Tors Force;
    Compute_force(&Force, Diff.coeff, NM);
    float k_Force_cart=Collectivity_norm2(Force.Cart, 3*N_ref);
    FILE *file_out=fopen(summary2, "a");
    fprintf(file_out, "Coll(force_Lambda)     %.0f\n", k_Force_cart);
    fclose(file_out);
    if(PRINT_FORCE){
      float F_module[N_ref];
      int j=0;
      for(i=0; i<N_ref; i++){
	F_module[i]=
	  Force.Cart[j]*Force.Cart[j]+
	  Force.Cart[j+1]*Force.Cart[j+1]+
	  Force.Cart[j+2]*Force.Cart[j+2];
	j+=3;
      }
      char nameforce[200]; sprintf(nameforce, "%s_force.pdb", nameout2);
      Print_cart_fluct(Ref_kin.atom_num, N_ref, F_module, atoms1, seq1,
		       nameforce,"FORCE");
      Print_force_confchange(NM, Diff, atoms1, axe1, naxe1, seq1, nres1,
			     nameout2);
    }
    printf("Computing and printing force. Time= %.2lf sec.\n",
	   (clock()-t0)/nbtops); t0=clock();
    Empty_tors(Force);
    Empty_tors(Diff); // double free, why?
    //if(Change_J)Empty_Jacobian(*J_CC);
    //if(diff_phi)free(diff_phi);
    //if(greedy_phi)free(greedy_phi);
  }
  /************************** End of conformation change **************/

  Print_tors_fluct2(NM, axe1, naxe1, seq1, nameout1);

  /********************** Print all modes *****************************/

  int N_print=N_MODE_PRINT;
  if(N_print>NM.N_relevant)N_print=NM.N_relevant;
  if(N_print){
    Print_modes(N_print, nameout1, "Cart", NM.select, N_cart,
		NM.Cart, NM.Cart_coll, NM.sigma2, sum_sigma2, axe1, naxe1,
		atoms1, natoms1, seq1, nres1, Ref_kin.atom_num, N_ref);
    Print_modes(N_print, nameout1, "Tors", NM.select, naxe1,
		NM.Tors, NM.Tors_coll, NM.sigma2, sum_sigma2, axe1, naxe1,
		atoms1, natoms1, seq1, nres1, Ref_kin.atom_num, N_ref);
    Print_modes(N_print, nameout1, "MWTors", NM.select, naxe1,
		NM.MW_Tors, NM.MW_Tors_coll, NM.sigma2, sum_sigma2, axe1,
		naxe1, atoms1, natoms1, seq1, nres1, Ref_kin.atom_num, N_ref);
    printf("Printing normal modes. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
    t0=clock();
  }

  // Print modes in PDB format
  if((PRINT_PDB)&&(N_MODE_PRINT)){
    //double RMSD=1.00, SDEV_SIM=RMSD*RMSD*Ref_kin.mass_tot;
    //double Coll_ave=0; for(ia=0; ia<NM.N; ia++)Coll_ave+=NM.Cart_coll[ia];
    //Coll_ave/=NM.N;
    int N_STEP=40, ip=0;
    for(ia=0; ia<NM.N; ia++){
      if((NM.sigma2[ia]==0)||(NM.select[ia]==0))continue;
      Print_mode_PDB(atoms1, natoms1, axe1, naxe1, bonds, seq1,
		     NM.Tors[ia], NM.omega[ia], N_STEP,
		     Para_simul, nameout1, ia, ip);
      ip++; if(ip==N_print)break;
    }
    Set_bonds_measure(bonds, natoms1, atoms1);
    printf("Printing modes as PDB. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
    t0=clock();
  }

  /**********************  Predict RMSD of mutations ***********************/
  if(PRED_MUT){
    Predict_mutations(NM, atoms1, natoms1, naxe1, Ref_kin,
		      Int_list, N_int, seq1, nres1, nameout1, Mut_para);
  }

  /**********************  Binding sites dynamics *************************/
  if(SITE_DYNAMICS){
    Binding_site_dynamics(NM, Ref_kin, seq1, atoms1, nameout1,
			  file_pdb1, chain1, Nchain, nres1, SITES);
  }

  /**********************  Dynamical couplings ****************************/
  if(ALLOSTERY){
    for(int anharm=0; anharm<n_har; anharm++){
      printf("Predicting dynamical couplings, anharmonicity=%d\n", anharm);
      Predict_allostery(NM, atoms1, Ref_kin, Int_list, N_int, seq1, nameout1,
			Confchange, file_pdb1, chain1, nres1, SITES, anharm);
    }
    printf("Dynamical coupling. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
    t0=clock();
  }

 /******************** Simulations  ************************/
  if(Para_simul.N_SIMUL){
    // To compute the ensemble NoRot, set strcpy(name_ens,"NoRot");
    // factor=factor_B*Para_simul.AMPLITUDE
    float fact; int N_ens=0;
    for(fact=Para_simul.AMPL_MIN; fact<=Para_simul.AMPL_MAX;
	fact*=Para_simul.AMPL_FACT){
      for(int anhar=0; anhar<=1; anhar++){
	if((ANHARMONIC==0)&&anhar)continue;
	char name_ens[40]; sprintf(name_ens,"Ampl_%.2g", fact);
	if(anhar==0){sprintf(name_ens,"%s_harmonic",name_ens);}
	else{sprintf(name_ens,"%s_anharmonic",name_ens);}
	if(Simulate_ensemble(Para_simul.N_SIMUL, fact, name_ens,
			     Para_simul, Mass_sqrt, atoms1, natoms1,
			     axe1, naxe1, bonds, seq1, nres1, N_diso1, NM, J,
			     Ref_kin, Int_list, N_int, Int_KB, INT_TYPE, s0,
			     nameout1, anhar)<0)break;
	N_ens++;
      }
    }
    printf("Drawing %d ensembles of conformations. Time= %.2lf sec.\n",
	   N_ens, (clock()-t0)/nbtops); t0=clock();
  }

  /******************** Example of Build up  ************************/
  /*if(BUILDUP){ // Example of application of build-up
    float fact=1.; int imode=0, nmove=10;
    Trajectory(atoms1,natoms1,axe1,naxe1,NM.Tors,imode,fact,nmove,bonds);
    }*/

  /**************  Computing conformation change given a force ********/
  if((FILE_FORCE!=NULL)&& (ANM==0)){
    float coord_ini[N_cart], coord_end[N_cart];
    Write_ref_coord_atom(coord_ini, N_ref, atoms1, Ref_kin.atom_num);
    char chain3=chain1[0];
    Force2Confchange(coord_end, FILE_FORCE, naxe1, N_ref, atoms1,
		     natoms1, Ref_kin.atom_num, J.Jacobian_ar, NM.omega2,
		     NM.Tors, NM.select, axe1, &chain3, nres1, seq1, bonds);
    rmsd=rmsd_mclachlan_f(coord_ini, coord_end, Ref_kin.mass_atom, N_ref);
    printf("rms force perturbed structure = %8.3f\n",rmsd);
  }

  /* Anharmonicity analysis
     int STEP_MAX=100;   // Number of steps for Move_struct or Move_all_modes
  if(0){ // It has to be corrected! 4/2/2013
    int step, moves;
    for(step=2; step<STEP_MAX; step++){
      moves=Move_all_modes(nameout1,axe1,naxe1,&NM,atoms1,natoms1,seq1,nres1,
			   atoms2, N_ref, Ref1.atom_num, Ref2.atom_num, step,
			   Para_simul, bonds);
      printf("%3d moves of all modes\n", moves);
      if(moves==0)break;
    }
    printf("Moving structure, all modes. Time= %.5lf sec.\n",
	   (clock()-t0)/nbtops); t0=clock();
	   }*/

  /*if(0 && CONF_CHANGE && PRINT_CONFCHANGE){
    printf("Simulating conformation change with normal modes\n");
    int MODE_MAX=50;    // Number of modes used for Move_struct
    float RMSD_STOP=0.1; // When Move_struct is stopped

    // Anharmonicity is computed here
    float rmsd_fin=rmsd;
    atom *atoms_sim=malloc(natoms1*sizeof(atom));
    for(i=0; i<natoms1; i++)atoms_sim[i]=atoms1[i];
    for(i=0; i<NM.N_relevant; i++)NM.Max_factor[i]=0;

    for(k=0; k<STEP_MAX; k++){
      int move=0; float rmsd_k=rmsd_fin;
      if(k){
	// WARNING: The matrix J is updated here
	J.N_kin=Compute_kinetic(&J,axe1,naxe1,atoms_sim,natoms1,Ref1,1);
      }
      for(i=0; i<MODE_MAX; i++){
	if(i>=NM.N_relevant)break;
	if(CONF_CHANGE){ia=NM.sort[i];}else{ia=i;}
	if((NM.sigma2[ia]==0)||(ia>=NM.N_relevant))continue;
	if(k){
	  Transform_tors_modes(NM.Tors[ia], NM.MW_Tors[ia],
			       J.T_sqrt_tr, J.T_sqrt_inv_tr,
			       N_modes, naxe1);
	  NM.Max_factor[ia]=0;
	}
	move+=Move_struct_confchange(atoms_sim, natoms1, &rmsd_fin,
				     nameout1, ia, 1, axe1, naxe1,
				     &NM, seq1, nres1, atoms2, N_ref,
				     Ref1.atom_num, Ref2.atom_num,
				     Para_simul, bonds);
	if(rmsd_fin < Para_simul.PDB_STEP)break;
      }
      printf("%2d round, %3d moves rmsd= %.2f\n", k+1, move, rmsd_fin);
      if((move==0)||(rmsd_fin < Para_simul.PDB_STEP))break;
      if(rmsd_k-rmsd_fin < RMSD_STOP)break;
    }
    free(atoms_sim);
    printf("Making conformation change. Time= %.5lf sec.\n",
	   (clock()-t0)/nbtops); t0=clock();
	   }*/

  /******************** Clean memory *************************/
   //modyves: several things below were never freed
  free(Int_list);						  
  free(bonds);						  
  free(atoms1);
  free(seq1);
  free(axe1);
  if(file_pdb2[0]!='\0'){
    free(atoms2);
    free(seq2);
    for(i=0;i<Nchain2;i++){
      free(chains2[i].seq);
      free(chains2[i].seqres);
      free(chains2[i].ali_seqres);
    }
    free(chains2);
    for(i=0;i<Nchain1;i++)free(chains1[i].alignres);
    free(ali.alignres);
    free(Posmut);
    free(AAmut);
  }
  for(i=0;i<Nchain1;i++){
    free(chains1[i].seq);
    free(chains1[i].seqres);
    free(chains1[i].ali_seqres);
  }
  free(chains1);
  if(CONF_CHANGE){
    free(coord_1);
    free(coord_2);
    free(Ref1.mass_atom);
    free(Ref1.mass_coord);
    free(Ref1.atom_num);
    free(ali_atoms.ali1);
    free(ali_atoms.ali2);
    free(ali_atoms.mass);
  }
  Clean_memory(NM, J, Hessian, naxe1, N_ref, N_cart, N_modes, ANM);
  Empty_Ref(&Ref_kin);
  return(0);
  
}

/******************  Other routines  ***************************/

double Rescale_force(struct Normal_Mode NM,
		     struct interaction *Int_list, int N_int,
		     float kappa, int anharmonic)
{
  double sum_sigma2=0; int i;
  printf("Rescaling force constant by factor %.3g\n", kappa);
  KAPPA*=kappa;
  K_PSI*=kappa; K_PHI*=kappa; K_OMEGA*=kappa; K_CHI*=kappa;
  K_BA*=kappa; K_BL*=kappa;
  // rescale interactions
  for(i=0; i<N_int; i++)Int_list[i].sec_der*=kappa;

  // rescale modes
  for(i=0; i<NM.N; i++){
    if(anharmonic==0){
      NM.omega2[i]*=kappa;
      NM.sigma2[i]/=kappa;
      if(NM.select[i]==0){NM.omega[i]=0;}
      else{NM.omega[i]=sqrt(NM.omega2[i]);}
      sum_sigma2+=NM.sigma2[i];
    }else{
      NM.sigma2_anhar[i]/=kappa;
      sum_sigma2+=NM.sigma2_anhar[i];
    }
  }
  return(sum_sigma2);
}

/*********************** Handling memory  **************************/

void Allocate_memory(struct Normal_Mode *NM,
		     struct Jacobian *J,
		     double ***Hessian,
		     int N_axes, int N_ref, int N_cart,
		     int N_modes, int N_res)
{

  // Hessian
  *Hessian=Allocate_mat2_d(N_modes, N_modes);
  printf("Hessian allocated\n");

  // Kinematics:
  Allocate_Jacobian(J, N_axes, N_cart);
  printf("Jacobian allocated\n");

  // Normal modes
  Allocate_Normal_modes(NM, N_modes, N_axes, N_cart);
  printf("Normal modes allocated\n");

}

void Allocate_Normal_modes(struct Normal_Mode *NM,
			   int N_modes, int N_axes, int N_cart)
{
  NM->N=N_modes;
  NM->N_kin=N_modes;
  NM->N_axes=N_axes;
  NM->N_cart=N_cart;
  int i;

  NM->select=malloc(N_modes*sizeof(int));
  NM->omega=malloc(N_modes*sizeof(float));
  NM->omega2=malloc(N_modes*sizeof(float));
  NM->sigma2=malloc(N_modes*sizeof(float));
  NM->Cart=Allocate_mat2_f(N_modes, N_cart);
  NM->Tors=Allocate_mat2_f(N_modes, N_axes);
  NM->MW_Tors=Allocate_mat2_f(N_modes, N_axes);
  NM->Cart_coll=malloc(N_modes*sizeof(float));
  NM->Tors_coll=malloc(N_modes*sizeof(float));
  NM->MW_Tors_coll=malloc(N_modes*sizeof(float));
  NM->Max_dev=malloc(N_modes*sizeof(float));
  NM->sort=malloc(N_modes*sizeof(int));
  NM->Tors_frac=malloc(N_modes*sizeof(float));
  for(i=0; i<N_modes; i++)NM->Tors_frac[i]=1.0;
  NM->tors_fluct=malloc(NM->N_axes*sizeof(float));
  NM->cart_fluct=malloc(NM->N_cart*sizeof(float));
  NM->confchange2=malloc(N_modes*sizeof(float));
  NM->sigma2_anhar=malloc(N_modes*sizeof(float));
  NM->d_KL=malloc(N_modes*sizeof(float));
  NM->Anharmonicity=malloc(N_modes*sizeof(float));
  NM->Anharm_struct=malloc(N_modes*sizeof(float));
  for(i=0; i<N_modes; i++){
    NM->Anharmonicity[i]=0; NM->Anharm_struct[i]=0;
  }
  NM->Max_factor=malloc(N_modes*sizeof(float));
  NM->Max_RMSD=malloc(N_modes*sizeof(float));
  for(i=0; i<N_modes; i++)NM->Max_RMSD[i]=-1;
}


void Empty_Normal_modes(struct Normal_Mode NM)
{

  // Cleaning normal modes
  Empty_matrix_f(NM.Cart, NM.N);
  Empty_matrix_f(NM.Tors, NM.N);
  Empty_matrix_f(NM.MW_Tors, NM.N);
  free(NM.omega2);
  free(NM.omega);
  free(NM.select); //modyves: was never freed
  free(NM.sigma2);
  free(NM.Cart_coll);
  free(NM.Tors_coll);
  free(NM.MW_Tors_coll);
  free(NM.Tors_frac);
  free(NM.sort);
  free(NM.Anharmonicity);
  free(NM.Anharm_struct);
  free(NM.Max_factor);
  free(NM.Max_RMSD);
  free(NM.Max_dev); 
  free(NM.cart_fluct);
  free(NM.tors_fluct);
  free(NM.confchange2);

}

void Clean_memory(struct Normal_Mode NM,
		  struct Jacobian J,
		  double **Hessian,
		  int N_axes, int N_ref, int N_cart,
		  int N_modes, int ANM)
{

  // Free structure 1
  // free(axe1);
  // free(atoms1);

  // Free Hessian
  Empty_matrix_d(Hessian, N_modes);

  // Normal modes
  Empty_Normal_modes(NM);

  // Cleaning reference set for kinetic energy
  Empty_Jacobian(J);

}


/*************************   Printing  *************************/

void Print_summary(char *name_out, char *name1, char *parameters, int Nchain,
		   struct Normal_Mode NM,
		   char *out_B1, char *out_B2, float Tors_fluct,
		   struct axe *axes, int N_axes, int N_main,
		   int Nskip, int N_atoms, int N_res, int N_diso,
		   int N_int, int N_inter, float Cont_overlap,
		   int N_disc_coll, int N_disc_freq, int N_disc_out,
		   float mass_sum)
{
  FILE *file_out;

  printf("Writing %s\n", name_out);
  file_out=fopen(name_out, "w");

  fprintf(file_out, "%s\n", parameters);
  fprintf(file_out, "protein               %s\n", name1);
  fprintf(file_out, "Nchain                %d\n", Nchain);
  fprintf(file_out, "Nres                  %d\n", N_res);
  fprintf(file_out, "Mass                  %.0f\n", mass_sum);
  fprintf(file_out, "Degrees of freedom    %d\n", NM.N);
  fprintf(file_out, "Side chain            %d\n", N_axes-N_main);
  fprintf(file_out, "Skipped main chain    %d\n", Nskip);
  fprintf(file_out, "Discarded modes coll. %d\n", N_disc_coll);
  fprintf(file_out, "Discarded modes freq. %d\n", N_disc_freq);
  fprintf(file_out, "Discarded modes outl. %d\n", N_disc_out);
  fprintf(file_out, "Disordered gaps       %d\n", N_diso);
  fprintf(file_out, "Native interactions   %d\n", N_int);
  if(Nchain>1)fprintf(file_out, "Interchain contacts    %d\n", N_inter);
  if(Cont_overlap >= 0)
    fprintf(file_out, "Overlap with MIN int. %.3f\n", Cont_overlap);

  int anharmonic=0;
  Print_thermal(anharmonic, file_out, NM, out_B1, axes, N_main, N_res, Nchain);
  fprintf(file_out, "Torsional_RMSD_Harm   %.3f\n", Tors_fluct);
  fprintf(file_out, "#\n");

  if(ANHARMONIC){
    anharmonic=1;
    Print_thermal(anharmonic, file_out, NM, out_B2, axes,N_main,N_res,Nchain);
    fprintf(file_out, "#\n");
  }
  fclose(file_out);
}

void Print_thermal(int anharmonic, FILE *file_out,
		   struct Normal_Mode NM, char *out_B,
		   struct axe *axes, int N_main,
		   int N_res, int Nchain)
{
  float *sigma2; char name[80]="";
  if(anharmonic){
    sigma2=NM.sigma2_anhar; strcpy(name, "anhar");
    fprintf(file_out, "# Anharmonic correction:\n");
  }else{
    sigma2=NM.sigma2;
    fprintf(file_out, "#\n# Harmonic approximation:\n");
  }

  double k_Therm=Collectivity_norm1(sigma2, NM.N); //Renyi
  fprintf(file_out, "Recp.Coll(therm)      %.0f\n", k_Therm);
  /*fprintf(file_out, "Area(therm)           %.3f\n",
    Area_norm1(NM.sigma2, NM.N));*/

  double Mean_Cart_coll=0, sum_w=0; int ia, i;
  for(ia=0; ia<NM.N; ia++){
    float w=sigma2[ia];
    Mean_Cart_coll+=NM.Cart_coll[ia]*w; sum_w+=w;
  }
  if(sum_w)Mean_Cart_coll/=sum_w;
  fprintf(file_out, "Ave_collectivity_NM   %.3f\n", Mean_Cart_coll);
  /*
  Coll_thr_cc=Mean_Cart_coll*0.6666;
  int N_modes_coll=0;
  for(ia=0; ia<NM.N; ia++){
    if((NM.select[ia])&&(NM.Cart_coll[ia] >= Coll_thr_cc))N_modes_coll++;
  }
  fprintf(file_out, "Num.collective_modes  %d\n", N_modes_coll);
  fprintf(file_out, "Threshold_collectivity %.3f\n", Coll_thr_cc);
  */

  {  /* Harmonic entropy quantum
	/h= 1.055 e-34 J.sec
	k_B T at 293K = 1.381 e-23 *293= 4.05e-21 J
	/h/k_BT = 0.260e-13 s^(-1)
	omega in units of s^(-1): sqrt(k_B T/M )(1/Dr)
	where 1/Dr = (omega) 10^10 omega
	sqrt(k_B T/M)= sqrt(1/M) in units of atomic masses *
	sqrt(4.05e-21 J/1.66e-27 Kg)=1562
	Time unit in s^-1: 1.562*10^13/sqrt(M) in atomic masses
     */
    double beta_h=0.26e-13; 
    double omega_scale=1.562e+13/sqrt(mass_sum);
    double beta_omega=beta_h*omega_scale;
    double entropy=0;
    for(i=0; i<NM.N; i++){
      if((NM.select[i]==0)||(sigma2[i]<=0))continue;
      float omega=1./sqrt(sigma2[i]);
      double x=beta_omega*omega, f=exp(-x);
      if(f<1)entropy+=(-log(1.-f)+x*f/(1.-f));
    }
    if(anharmonic==0)
      fprintf(file_out,"kT//h=                %.3g\n", 1./beta_omega);
    fprintf(file_out, "Entropy per res.      %.3g\n",  entropy/N_res);
  }

  fprintf(file_out, "%s", out_B);

  // Rigid body degrees of freedom
  if(Nchain>1){
    double tra=0, rot=0, all=0;
    for(ia=0; ia<NM.N; ia++){
      float *v=NM.MW_Tors[ia];
      double tra_ia=0, rot_ia=0;
      for(i=0; i<N_main; i++){
	if(axes[i].rigid){
	  if(axes[i].type=='l'){tra_ia+=(*v)*(*v);}
	  else{rot_ia += (*v)*(*v);}
	}
	v++;
      }
      tra+=tra_ia*sigma2[ia];
      rot+=rot_ia*sigma2[ia];
      all+=sigma2[ia];
    }
    fprintf(file_out, "Fraction_internal_tra %.3f\n", tra/all);
    fprintf(file_out, "Fraction_internal_rot %.3f\n", rot/all);
  }

}

/***************************  Input **********************************/

int getArgs(int argc, char **argv,
	    char *file_pdb1, char *chain1,
            char *file_pdb2, char *chain2,
	    int *ANM, char *REF, int *LABEL,
	    int *SIDECHAINS, int *OMEGA, int *PSI,
	    double *K_OMEGA, double *K_PSI,
	    double *K_PHI, double *K_CHI,
	    double *K_BA, double *K_BL,
	    float *E_MIN, float *COLL_THR,
	    int *MIN_INT, char *INT_TYPE,
	    float *C_THR, int *S_TYPE, float *S_THR,
	    int *ONEINT, int *N_RESRES,
	    int *N_MODE, char *outdir,
	    int *PRINT_CONFCHANGE, int *PRINT_FORCE,
	    int *PRINT_PDB, int *PRINT_MODE_SUMM,
	    char **FILE_FORCE, int *ALLOSTERY,
	    float *KAPPA, int *FIT_B, float *RMSD_EXP,
	    float *RMSD_MIN, int *NMUT, int *PRED_MUT,
	    char *Mut_para,
	    struct Para_simul *Para_simul,
	    struct Para_confchange *Para_confchange,
	    float *SEQID_THR,
	    char *parameters,
	    int *PRINT_COV_COUPLING,
	    int *PRINT_DEF_COUPLING,
	    int *PRINT_DIR_COUPLING,
	    int *PRINT_COORD_COUPLING,
	    int *PRINT_SIGMA_DIJ,
	    char *PROF_TYPE,
	    int *ALL_PAIRS, float *SIGMA,
	    int *STRAIN, char *SITES,
	    int *PRINT_CMAT, int *ANHARMONIC,
	    int *PRINT_PDB_ANHARM, int *N_PDB_print)
{
  int i; //p1=0, p2=0, out=0


  /************** Initialize ****************/
  // Input
  file_pdb1[0]='\0'; strcpy(chain1, "");
  file_pdb2[0]='\0'; strcpy(chain2, "");
  // Model variables
  *ANM=0;  *C_THR=THR_ALL;
  INT_TYPE[0]='\0'; REF[0]='\0';
  EXP_HESSIAN=EXP_HESSIAN_DEF;
  POW=POW_DEF;
  *MIN_INT=MIN_INT_DEF;
  *K_PSI=K_PSI_DEF;
  *K_PHI=K_PHI_DEF;
  *K_CHI=K_CHI_DEF;
  *K_OMEGA=K_OMEGA_DEF;
  DA_THR=DA_THR_DEF;
  COS_DAA=COS_DAA_DEF;
  ENE_HB=ENE_HB_DEF;
  KINETIC=1;
  // Output
  *N_MODE=0;
  *PRINT_PDB_ANHARM=0; *N_PDB_print=5;

  // Analysis of conformation change
  Para_confchange->RMSD_THR=RMSD_THR_DEF; // Minimum confchange per mode
  Para_confchange->COLL_THR=COLL_THR_DEF; // Threshold on collectivity
  *SEQID_THR=SEQID_THR_DEF;

  // Simulations
  Para_simul->N_SIMUL=0;
  Para_simul->AMPLITUDE=AMPLITUDE_DEF;
  Para_simul->AMPL_MIN=AMPLITUDE_DEF;
  Para_simul->AMPL_MAX=AMPLITUDE_DEF;
  Para_simul->AMPL_FACT=1.0;
  Para_simul->D_REP=D_REP_DEF;
  Para_simul->E_THR=E_THR_DEF;
  Para_simul->PDB_STEP=PDB_STEP_DEF;
  Para_simul->MAX_ANGLE=MAX_ANGLE_DEF;
  Para_simul->RESET=0;
  Para_simul->SELECT_ENE=0;
  // Conformation changes
  Para_simul->NSTEPS=100;
  Para_simul->STEP_MAX=0.5;
  Para_simul->STEP_MIN=0.001;
  Para_simul->ANGLE=0.05;

  if(argc<2)help();
  int infile=Read_para(argv[1], file_pdb1, chain1,file_pdb2, chain2,
		       ANM, REF, LABEL, SIDECHAINS, OMEGA, PSI,
		       K_OMEGA, K_PSI, K_PHI, K_CHI, K_BA, K_BL,
		       E_MIN, COLL_THR, MIN_INT, INT_TYPE, C_THR,
		       S_TYPE, S_THR, ONEINT, N_RESRES,
		       N_MODE, outdir, PRINT_CONFCHANGE, PRINT_FORCE,
		       PRINT_PDB, PRINT_MODE_SUMM, FILE_FORCE, ALLOSTERY,
		       KAPPA, FIT_B, RMSD_EXP, RMSD_MIN,
		       NMUT, PRED_MUT, Mut_para,
		       Para_simul, Para_confchange, SEQID_THR,
		       PRINT_COV_COUPLING, PRINT_DEF_COUPLING,
		       PRINT_DIR_COUPLING, PRINT_COORD_COUPLING,
		       PRINT_SIGMA_DIJ, PROF_TYPE, ALL_PAIRS, SIGMA,
		       STRAIN, SITES, PRINT_CMAT, ANHARMONIC,
		       PRINT_PDB_ANHARM, N_PDB_print);
  if(infile)goto inform;
  printf("WARNING, input file not specified or absent (%s)\n", argv[1]);
  printf("WARNING, reading parameters from command line\n");

  for(i=1; i<argc; i++){
      // Proteins:
    if (strncmp(argv[i],"-p2",3)==0){
      i++; if(i>=argc)continue;
      strcpy(file_pdb2,argv[i]); //p2=1;
      printf("file_pdb2=%s\n",file_pdb2);
    }else if (strncmp(argv[i],"-c2",3)==0){
      i++; if(i>=argc)continue;
      if((strncmp(argv[i], "all", 3)==0)||(strncmp(argv[i], "ALL",3)==0))
	{*chain2='*';}
      else{strcpy(chain2,argv[i]);}
      printf("chain2=%s\n",chain2);
    }else if (strncmp(argv[i],"-p1",3)==0){
      i++; if(i>=argc)continue;
      strcpy(file_pdb1,argv[i]); //p1=1;
      printf("file_pdb1=%s\n",file_pdb1);
    }else if (strncmp(argv[i],"-c1",3)==0){
      i++; if(i>=argc)continue;
      if((strncmp(argv[i], "all",3)==0)||(strncmp(argv[i], "ALL",3)==0))
	{*chain1='*';}
      else{strcpy(chain1,argv[i]);}
      printf("chain(s) to read: %s\n",chain1);

      // Model:
    }else if (strncmp(argv[i],"-anm",4)==0){
      *ANM=1;
    }else if (strncmp(argv[i],"-hnm",4)==0){
      *ANM=1; HNM=1;
    }else if (strncmp(argv[i],"-sc",3)==0){
      *SIDECHAINS=1;
    }else if (strncmp(argv[i],"-omega",6)==0){
      *OMEGA=1;
    }else if (strncmp(argv[i],"-ref",4)==0){
      i++; if(i>=argc)continue;
      strcpy(REF, argv[i]);
    }else if (strncmp(argv[i],"-cont_type",10)==0){
      i++; if(i>=argc)continue;
      strcpy(INT_TYPE, argv[i]);
    }else if (strncmp(argv[i],"-expo", 5)==0){
      i++; if(i>=argc)continue;
      sscanf(argv[i], "%f", &EXP_HESSIAN);
    }else if (strncmp(argv[i],"-cont_thr",9)==0){
      i++; if(i>=argc)continue;
      sscanf(argv[i], "%f", C_THR);
      printf("Threshold distance= %.1f\n", *C_THR);
    }else if (strncmp(argv[i],"-nokin",6)==0){
      KINETIC=0;

      // Force
    }else if (strncmp(argv[i],"-force",6)==0){
      i++; if(i>=argc)continue;
      *FILE_FORCE=malloc(400*sizeof(char));
      strcpy(*FILE_FORCE, argv[i]);
      // Simulation
    }else if (strncmp(argv[i],"-simul",6)==0){
      i++; if(i>=argc)continue;
      sscanf(argv[i], "%d", &(Para_simul->N_SIMUL));
      printf("Simulations to be done= %d\n", Para_simul->N_SIMUL);
      // Output:
    }else if (strncmp(argv[i],"-outdir",7)==0){
      i++; if(i>=argc)continue;
      strcpy(outdir, argv[i]); //out=1;
      printf("output directory= %s\n",outdir);
    }else if (strncmp(argv[i],"-modes",6)==0){
      i++; if(i>=argc)continue;
      *N_MODE=atoi(argv[i]);  //sscanf(argv[i], "%d", N_MODE);
      printf("%d normal modes to be printed\n",*N_MODE);
    }else if (strncmp(argv[i],"-print_pdb",10)==0){
      *PRINT_PDB=1;
    }else if (strncmp(argv[i],"-allostery",8)==0){
      *ALLOSTERY=1;
    }else if (strncmp(argv[i],"-print_confchange",13)==0){
      *PRINT_CONFCHANGE=1;
    }else if (strncmp(argv[i],"-print_force",12)==0){
      *PRINT_FORCE=1;
    }else if (strncmp(argv[i],"-print_summ",11)==0){
      *PRINT_MODE_SUMM=1;
    }else if (strncmp(argv[i],"-debug",6)==0){
      DEBUG=1;
    }else if (strncmp(argv[i],"-h",2)==0){
      help();
    }else if (argv[i][0]=='-'){
      printf("WARNING, argument %s does not exist\n", argv[i]);
    }
  }

 inform:
  if(*SEQID_THR > 1.0)(*SEQID_THR)/=100;
  /*if(*SEQID_THR > 1.0)
    printf("WARNING, SEQID_THR=%.2f > 1, performing Mammoth alignment\n",
    *SEQID_THR);*/

  if(file_pdb1[0]=='\0'){
    printf("ERROR, at least one PDB file is mandatory\n");
    help();
  }
  if(*OMEGA>0){
    printf("Omega angle taken as degree of freedom.\n");
  }else if(*OMEGA<0){
    printf("Omega angle taken as dof if Domega>%d degrees.\n", D_omega_thr);
  }else{
    *K_OMEGA=0;
  }
  if(*PSI==0){
    printf("Psi angle frozen\n"); *K_PSI=0;
  }
  if(*SIDECHAINS){
    printf("Side chains degrees of freedom used\n");
    if((strncmp(INT_TYPE, "MIN", 3)!=0)||(strncmp(INT_TYPE, "SCR", 3)!=0)){
      printf("WARNING, interaction type %s not allowed with side chains\n",
	     INT_TYPE);
      strcpy(INT_TYPE, "MIN");
    }
    if((strncmp(REF, "ALL", 3)!=0)){
      printf("WARNING, reference atoms %s not allowed with side chains\n",
	     REF);
      strcpy(REF, "ALL");
    }
    printf("Reference atoms set to %s, interaction type set to %s\n",
	   REF, INT_TYPE);
    printf("Minimum number of interactions for accepting a degree of freedom");
    printf(": %d\n", *MIN_INT);
  }else{
    *K_CHI=0;
  }
  if((strncmp(INT_TYPE, "CA", 2)!=0)&&
     (strncmp(INT_TYPE, "CB", 2)!=0)&&
     (strncmp(INT_TYPE, "SCA", 3)!=0)&&
     (strncmp(INT_TYPE, "SCB", 3)!=0)&&
     (strncmp(INT_TYPE, "MIN", 3)!=0)&&
     (strncmp(INT_TYPE, "HYD", 3)!=0)&&
     (strncmp(INT_TYPE, "ALL", 3)!=0)&&
     (strncmp(INT_TYPE, "SCR", 3)!=0)&&
     (strncmp(INT_TYPE, "SHA", 3)!=0)&&
     (strncmp(INT_TYPE, "HB", 2)!=0)){
    if(INT_TYPE[0]!='\0')
      printf("WARNING, interaction type %s does not exist\n", INT_TYPE);
    printf("Interaction type set to default %s\n", INT_TYPE_DEF);
    strcpy(INT_TYPE, INT_TYPE_DEF);
  }else{
    printf("Contact type set to %s\n", INT_TYPE);
    if(strncmp(INT_TYPE, "HB", 2)==0){
      printf("Estimating hydrogen bonds, ");
      printf("Min.dist= %.3f max.cos(DAA')= %.3f Force factor= %.3f\n",
	     DA_THR, COS_DAA, ENE_HB);
    }
  }

  if(REF[0]!='\0'){
    if((strncmp(REF, "CA", 2)!=0)&&
       (strncmp(REF, "CB", 2)==0)&&
       (strncmp(REF, "BB", 2)==0)&&
       (strncmp(REF, "EB", 2)==0)&&
       (strncmp(REF, "ALL", 3)==0)){
      printf("WARNING, %s are not allowed reference atoms\n", REF);
      REF[0]='\0';
    }
  }
  if(EXP_HESSIAN<0){
    printf("ERROR, exponent must be positive (%.3f), setting it to zero\n",
	   EXP_HESSIAN); EXP_HESSIAN=0;
  }else if (EXP_HESSIAN){
    printf("Exponent of force constant = %.3f\n", EXP_HESSIAN);
  }


  if(*C_THR==0){
    if(strncmp(INT_TYPE, "CA", 2)==0){
      *C_THR=THR_CA;
    }else if(strncmp(INT_TYPE, "CB", 2)==0){
      *C_THR=THR_CB;
    }else if(strncmp(INT_TYPE, "SCB", 3)==0){
      *C_THR=THR_ALL;
    }else if(strncmp(INT_TYPE, "ALL", 3)==0){
      *C_THR=THR_ALL;
    }else if(strncmp(INT_TYPE, "MIN", 3)==0){
      *C_THR=THR_ALL;
    }else if(strncmp(INT_TYPE, "HYD", 3)==0){
      *C_THR=THR_ALL;
    }else if(strncmp(INT_TYPE, "SCR", 3)==0){
      *C_THR=THR_CA;
    }else if(strncmp(INT_TYPE, "HB", 2)==0){
      *C_THR=THR_ALL;
    }else{
      printf("ERROR, wrong interaction type %s\n", INT_TYPE);
      exit(8);
    }
    printf("Contact threshold set to default %.1f\n", *C_THR);
  }

  if(*ANM){
    printf("Degrees of freedom: Cartesian (ANM)\n");
    if(REF[0]!='\0')
      printf("WARNING, reference atoms can not be assigned in the ANM\n");
  }else{
    printf("Degrees of freedom: Torsional (TNM, default)\n");
  }

  char tmpstr[500];
  //modyves: some changes to not sprintf parameters into parameters, which gave me some warnings.
  sprintf(parameters, "Para: Psi=%d  Omega=%d Sidechain=%d Min_int=%d  REF=%s",
	  *PSI, *OMEGA, *SIDECHAINS, *MIN_INT, REF);
  sprintf(tmpstr," E_MIN=%.2g COLL_THR=%.0f", *E_MIN, *COLL_THR);
  strcat(parameters,tmpstr);
  sprintf(tmpstr," K_O=%.2g K_F=%.2g K_P=%.2g K_C=%.2g K_BA=%.2g K_BL=%.2g",
	  *K_OMEGA, *K_PHI, *K_PSI, *K_CHI, *K_BA, *K_BL);
  strcat(parameters,tmpstr);
  sprintf(tmpstr," CONT=%s ", INT_TYPE); strcat(parameters,tmpstr);
  if(strncmp(INT_TYPE, "SCR", 3)==0){
    sprintf(tmpstr, " Tolerance=%.2f THR=%.2f", *S_THR, *C_THR);
  }else if(strncmp(INT_TYPE, "SHA", 3)==0){
    sprintf(tmpstr, " type=%d d=%.2f THR=%.1f", *S_TYPE, *S_THR, *C_THR);
  }else{
    sprintf(tmpstr, " thr=%.2f", *C_THR);
  }
  strcat(parameters,tmpstr);
  sprintf(tmpstr, " EXPO=%.1f ", EXP_HESSIAN); strcat(parameters,tmpstr);
  if(*FIT_B==0){sprintf(tmpstr, " KAPPA=%.1f", *KAPPA);}
  else{sprintf(tmpstr, " KAPPA obtained from fit");}
  strcat(parameters,tmpstr);
  fflush(stdout);
  return (0);
}

void help(void)
{
  fprintf(stderr, "Program %s\n", PRG);
  fprintf(stderr, "author Ugo Bastolla <ubastolla@cbm.csic.es> ");
  fprintf(stderr, "Centro de Biologia Molecular Severo Ochoa (CSIC-UAM)\n");
  fprintf(stderr, "Computes normal modes of proteins with the TNM or ANM.\n");
  fprintf(stderr, "If two structures are provided, it projects the ");
  fprintf(stderr, "experimental conformation change over normal modes\n");
  fprintf(stderr, "Inputs can be given either by file or by command line\n");
  fprintf(stderr, "   USAGE:");
  fprintf(stderr, "   %s <inputile>  or\n", PRG);
  fprintf(stderr, "   %s -p1 <pdbfile1>", PRG);
  fprintf(stderr, "   (reference structure for computing normal modes)\n");
  fprintf(stderr, "   OPTIONS:\n");
  fprintf(stderr, "       -h prints this help\n");
  fprintf(stderr, "       -c1 <chain_id1> <all>:reading all chains, A, AB..\n");
  fprintf(stderr, "       -p2 <pdbfile2> for conformation change\n");
  fprintf(stderr, "       -c2 <chain_id2> <all>:reading all chains, A, AB..\n");
  fprintf(stderr, "       -ref Reference atoms. Allowed: CA CB BB EB ALL. ");
  fprintf(stderr, "       -omega  Use also omega angle as degree of freedom\n");
  fprintf(stderr, "Default: %s\n", REF_DEF);
  fprintf(stderr, "       -anm use ANM d.o.freedom (by default TNM is used)\n");
  fprintf(stderr, "Interaction model:\n");
  fprintf(stderr, "       -cont_type Contact type CA CB, ALL, MIN HYD or HB. ");
  fprintf(stderr, "Default: %s\n", INT_TYPE_DEF);
  fprintf(stderr, "        if -cont_type=HB hydrogen bonds are estimated\n");
  fprintf(stderr, "       -cont_thr Distance threshold default %.1f\n",THR_CB);
  fprintf(stderr, "       -expo <e> k~r^(-e) Default 0\n");
  fprintf(stderr, "       -hnm HNM (ref=CA, cont_type=CA, e=6, covalent)\n");
  fprintf(stderr, "       -debug  print debugging information\n");
  fprintf(stderr, "       -outdir <output_directory>\n");
  fprintf(stderr, "       -modes <number of modes printed>\n");
  fprintf(stderr, "       -print_pdb    Print modes in PDB format\n");
  fprintf(stderr, "       -print_summ   Print modes summary\n");
  fprintf(stderr, "       -print_confchange Print conformation change\n");
  fprintf(stderr, "       -print_force  Print force induced by c.change\n");
  fprintf(stderr, "       -force        PDB file with force as coordinates\n");
  fprintf(stderr, "       -simul        <Num. simulated structures>\n");
  fprintf(stderr, "       -allostery Predict allostery\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "FORMAT of input file:\n");
  fprintf(stderr, "####   INPUT:\n");
  fprintf(stderr, "PDB1= /data/ortizg/databases/pdb/1usg.pdb ");
  fprintf(stderr, "! Reference structure\n");
  fprintf(stderr, "CH1=  A                                   ");
  fprintf(stderr, "! Chain\n");
  fprintf(stderr, " Conformation change (optional):\n");
  fprintf(stderr, "PDB2= /data/ortizg/database/pdb/1usk.pdb");
  fprintf(stderr, "! Conformation change\n");
  fprintf(stderr, "CH2=  A                                  ");
  fprintf(stderr, "! Chain\n");
  fprintf(stderr, "#\n");
  fprintf(stderr, "####  Model parameters (reccomended: do not change)\n");
  fprintf(stderr,
	  "#===============================================================\n");
  fprintf(stderr, "DOF=  TORS  (CART)		      ");
  fprintf(stderr, "! Internal or Cartesian deg. of freed.\n");
  fprintf(stderr, "REF=  ALL (allowed: ALL EB CB  CA BB)     ");
  fprintf(stderr, "! Reference atoms\n");
  fprintf(stderr, "PSI=1                              ");
  fprintf(stderr, "! Allow PSI angles to rotate?\n");
  fprintf(stderr, "OMEGA=0                              ");
  fprintf(stderr, "! Allow OMEGA angles to rotate?\n");
  fprintf(stderr, "SIDECHAIN=0                              ");
  fprintf(stderr, "! Allow sidechain angles to rotate?\n");
  fprintf(stderr,
	  "# The model with fewer degrees of freedom is faster and more robust w.r.t. parameter values\n");
  fprintf(stderr,
	  "# More degrees of freedom may introduce non collective motions that worsen performances.\n");
  fprintf(stderr, "MIN_INT=1                              ");
  fprintf(stderr, "! Min. number of interactions to accept a dof\n");
  fprintf(stderr, "E_MIN=0.00001                              ");
  fprintf(stderr, "! Min. eigenvalue/<evalue> of a normal mode\n");
  fprintf(stderr, "COLL_THR=30                                ");
  fprintf(stderr, "! Discard modes that move < COLL_THR atoms\n");
  fprintf(stderr, "K_BL=%.2f                                ", K_BA);
  fprintf(stderr, "! Elastic constant for bond lengths\n");
  fprintf(stderr, "K_BA=%.2f                                ", K_BA);
  fprintf(stderr, "! Elastic constant for bond angles\n");
  fprintf(stderr, "K_OMEGA=%.2f                             ", K_OMEGA);
  fprintf(stderr, "! Elastic constant for angle omega\n");
  fprintf(stderr, "K_PHI=%.2f                                ", K_PHI);
  fprintf(stderr, "! Elastic constant for angle phi\n");
  fprintf(stderr, "K_PSI=%.2f                                ", K_PSI);
  fprintf(stderr, "! Elastic constant for angle psi\n");
  fprintf(stderr, "K_CHI=%.1f                                ", K_CHI);
  fprintf(stderr, "! Elastic constant for angle Chi\n");
  fprintf(stderr, "CONT_TYPE= MIN ");
  fprintf(stderr, "! Interaction model (MIN SCR SHA SCB CB CA HB HNM)\n");
  fprintf(stderr, "CONT_THR=  4.5            ");
  fprintf(stderr, "! Threshold for contacts\n");
  /*fprintf(stderr, "S_THR= 0.2                          ");
  fprintf(stderr, "! Parameter for screened or shadow interactions\n");
  fprintf(stderr, "S_TYPE= 0.2                             ");
  fprintf(stderr, "! Type of shadow interactions\n");
  fprintf(stderr, "# 0=angle 1=axis distance 2=S-ball-axis-distance\n");
  fprintf(stderr, "ONEINT= 1                             ");
  fprintf(stderr, "! Only one interaction per residue pair?\n");
  fprintf(stderr, "N_RESRES= 1                           ");
  fprintf(stderr, "! Max num interactions per residue pair\n");*/
  fprintf(stderr, "POW= 1                                   ");
  fprintf(stderr, "! Force constant power-law (1) or exponential (0)\n");
  fprintf(stderr, "EXP_HESSIAN= 0                           ");
  fprintf(stderr, "! Force constant k~r^(-e)\n");
  fprintf(stderr, "FIT_B= 1                            ");
  fprintf(stderr, "! Determine force constant by fitting B factors if present? (default YES)\n");
  fprintf(stderr, "KAPPA= 218.4                        ");
  fprintf(stderr, "! force constant when not fitting B factors\n");
  fprintf(stderr,
	  "#==============================================================\n");
  fprintf(stderr, "############### OUTPUT ################\n");
  fprintf(stderr, "#\n### General output:\n");
  fprintf(stderr, "ANHARMONIC= 0  ");
  fprintf(stderr, "! Examine anharmonicity of each mode?\n");
  fprintf(stderr, "PRINT_PDB_ANHARM=0   ");
  fprintf(stderr, "! Print PDB while computing the anharmonic energy?\n");
  fprintf(stderr, "N_PDB_print= 4        ! Number of PDB printed ");
  fprintf(stderr, "for each norm.mode (if PRINT_PDB_ANHARM)\n");
  fprintf(stderr, "LABEL= 0 or 1			    ");
  fprintf(stderr, "! Print Mod. par. in file names\n");
  fprintf(stderr, "PRINT_CONT_MAT=0                         ");
  fprintf(stderr, "! Print Contact matrix?");
  fprintf(stderr, "PRINT_SUMM= 1                            ");
  fprintf(stderr, "! Print summary of results\n");
  fprintf(stderr, "DEBUG= 0                                 ");
  fprintf(stderr, "! Print debugging information\n");

  fprintf(stderr, "#\n### Normal modes:\n");
  fprintf(stderr, "NMODES= 5                                ");
  fprintf(stderr, "! Number of modes to print\n");
  fprintf(stderr, "PRINT_PDB= 0                             ");
  fprintf(stderr, "! Print modes as PDB files?\n");
 fprintf(stderr, "AMPLITUDE=1.0             ! Amplitude of printed modes w.r.t. thermal fluctuations\n");
  fprintf(stderr, "PDB_STEP=0.5              ! Minimum Rmsd between structures printed as PDB\n");
  fprintf(stderr, "E_THR=5.0                 ! Maximum allowed energy of generated structures\n");
  fprintf(stderr, "D_REP=2.5                 ! Maximum distance for computing repulsion energy\n");

 fprintf(stderr, "#\n#### Ensembles of alternative structures\n");
 
  fprintf(stderr, "SIMUL=20		     ! Number of simulated structures\n");
  fprintf(stderr, "AMPL_MIN=1.0              ! Minimum amplitude of simulated structures\n");
  fprintf(stderr, "AMPL_MAX=8.0              ! Maximum amplitude of simulated structures\n");
  fprintf(stderr, "AMPL_FACT=2.0             ! Factor multiplying conecutive amplitudes\n");

  fprintf(stderr, "#\n###  Analysis of conformation change\n");
  fprintf(stderr, "RMSD_MIN= 0.5                            ");
  fprintf(stderr, "! Min. RMSD for analyzing conformation change\n");
  fprintf(stderr, "RMSD_THR= 0.2                            ");
  fprintf(stderr, "! Do not analyze modes that move < RMSD_THR\n");
  fprintf(stderr, "PRINT_CHANGE= 0                          ");
  fprintf(stderr, "! Print conformation change?\n");
  fprintf(stderr, "STEP_MAX= 0.5                            ");
  fprintf(stderr, "! Max. step in confchange\n");
  fprintf(stderr, "STEP_MIN= 0.001                          ");
  fprintf(stderr, "! Smallest step in confchange\n");
  fprintf(stderr, "NSTEPS= 100                              ");
  fprintf(stderr, "! Number of steps in confchange\n");
  fprintf(stderr, "ANGLE= 0.1                               ");
  fprintf(stderr, "! Angular step in confchange\n");
  fprintf(stderr, "PRINT_FORCE= 0                           ");
  fprintf(stderr, "! Print force from linear response\n");
  fprintf(stderr, "#FILE_FORCE= \n");
  fprintf(stderr, "! File with force, compute deformation\n");

  fprintf(stderr, "#\n#### Prediction and analysis of mutational effects\n");
  fprintf(stderr, "PRED_MUT=0                           ");
  fprintf(stderr, "! Predict RMSD of all possible mutations?\n");
  fprintf(stderr, "MUT_PARA=Mutation_para.in           ");
  fprintf(stderr, "! File with mutation parameters\n");
  fprintf(stderr, "NMUT= -2                            ");
  fprintf(stderr, "! Analyze conf.change as mutation if number of mut==NMUT\n");
  fprintf(stderr, "                                    ");
  fprintf(stderr, "! NMUT=-1: Analyze always NMUT=-2: Mutation not analyzed\n");

  fprintf(stderr, "#\n### Dynamical couplings\n");
  fprintf(stderr, "ALLOSTERY=0                              ");
  fprintf(stderr, "! Compute dynamical couplings\n");
  fprintf(stderr, "ALL_PAIRS=0         ! Print couplings for all pairs? (1=YES)\n");
  fprintf(stderr, "SIGMA=1.0           ! If ALL_PAIRS=0, couplings are printed\n");
  fprintf(stderr, "# only for pairs with Coupling_ij > SIGMA*Std.dev\n");
  fprintf(stderr, "PROF_TYPE= C        ! Type of profile of the couplings p_i({C_kl})\n");
  fprintf(stderr, "# A=average E=effective connectivity P=principal eigenvector\n");
  fprintf(stderr, "PRINT_COV_COUPLING=0                    ");
  fprintf(stderr, "! Print covariance couplings\n");
  fprintf(stderr, "PRINT_DEF_COUPLING=0                    ");
  fprintf(stderr, "! Print deformation couplings\n");
  fprintf(stderr, "# Allosteric coupling A_ij is the deformations produced on site j by an\n");
  fprintf(stderr, "# unitary force applied in j in the direction that maximizes the deformation.\n");
  fprintf(stderr, "# The coupling is output only when it is > SIGMA standard deviations above the mean value.\n");
  fprintf(stderr, "# Output file: <>_deformation_coupling.dat\n");
  fprintf(stderr, "PRINT_DIR_COUPLING= 0	! Print directionality couplings\n");
  fprintf(stderr, "# Directionality coupling D_ij is the Boltzmann average of the scalar\n");
  fprintf(stderr, "# product of the direction of motion of the residues i and j. If it is\n");
  fprintf(stderr, "# positive the two residues tend to move in the same direction.\n");
  fprintf(stderr, "# The coupling is output only when it is > SIGMA standard deviations above\n");
  fprintf(stderr, "# Output file: <>_directionality_coupling.dat and\n");
  fprintf(stderr, "# <>_directionality_coupling_neg.dat (< -SIGMA standard deviations)\n");
  fprintf(stderr, "# It is output in the format of a N*N matrix for all pairs of sites\n");
  fprintf(stderr, "# Output file: <>_directionality_matrix.dat\n");
  fprintf(stderr, "PRINT_COORD_COUPLING=0			  ! Print coordination couplings\n");
  fprintf(stderr, "# Coordination coupling C_ij is the Boltzmann average of the squared\n");
  fprintf(stderr, "# fluctuations of the distance d_ij with respect to the equilibrium value.\n");
  fprintf(stderr, "# If it is small the two residues maintain an almost fixed distance during\n");
  fprintf(stderr, "# their dynamics.\n");
  fprintf(stderr, "# The coupling is output only when it is > SIGMA standard deviations above\n");
  fprintf(stderr, "# the mean value.\n");
  fprintf(stderr, "# Output file: <>_coordination_coupling.dat\n");
  fprintf(stderr, "PRINT_SIGMA_DIJ=0   ! Print variance of interatomic distance (only CA atoms)\n");
  fprintf(stderr, "# Note that it is the square of the coordination coupling\n");
  fprintf(stderr, "# Output file: <>_interatomic_distance_variance.dat");
  fprintf(stderr, "SITES=leut_sites.in	! File where an active site is read.\n");
  fprintf(stderr, "# If not specified, read from the SITE record in the PDB file.\n");
  fprintf(stderr, "# format: SITE AA (3 letter) CHAIN RESNUM (PDB)\n");
  fprintf(stderr, "# The mean coupling of pairs of residues in the active site is computed.\n");
  fprintf(stderr, "#\n");
  exit(1);
}

void GetPdbId(char *pdb_file_in, char *pdbid){
     /* This subroutine pretends to get the
        PDB id from a pdb file name, and ressembles
	quite a lot my "old-and-dirty-Perl" days */

  int start=0, end=0, i,j; //end2=0

  for(i=strlen(pdb_file_in)-1;i>=0;i--){
    if (pdb_file_in[i]=='.'){
      end=i-1;
    }else if (pdb_file_in[i]=='/'){
      start=i+1; //end2=i-1;
      break;
    }
  }
  j=0;
  for (i=start;i<=end;i++){
    pdbid[j]=pdb_file_in[i];
    j++;
  }
  pdbid[j]='\0';
}

void Copy_vector(float *xx, float *yy, int n){
  int i; float *x=xx, *y=yy;
  for(i=0; i<n; i++){*x=*y; x++; y++;}
}

void Write_ref_coord_d(double *coords, int N_ref, atom *atoms, int *atom_num)
{
  double *coord=coords;
  for(int i=0; i<N_ref; i++){
    double *r=atoms[atom_num[i]].r;
    for(int j=0; j<3; j++){*coord=*r; coord++; r++;}
  }
}

int Read_para(char *filename,
	      char *file_pdb1, char *chain1,
	      char *file_pdb2, char *chain2,
	      int *ANM, char *REF, int *LABEL,
	      int *SIDECHAINS, int *OMEGA, int *PSI,
	      double *K_OMEGA, double *K_PSI,
	      double *K_PHI, double *K_CHI,
	      double *K_BA, double *K_BL,
	      float *E_MIN, float *COLL_THR, int *MIN_INT,
	      char *INT_TYPE, float *C_THR,
	      int *S_TYPE, float *S_THR,
	      int *ONEINT, int *N_RESRES,
	      int *N_MODE, char *outdir,
	      int *PRINT_CONFCHANGE, int *PRINT_FORCE,
	      int *PRINT_PDB, int *PRINT_MODE_SUMM,
	      char **FILE_FORCE, int *ALLOSTERY,
	      float *KAPPA, int *FIT_B, float *RMSD_EXP,
	      float *RMSD_MIN, int *NMUT, int *PRED_MUT, char *Mut_para,
	      struct Para_simul *Para_simul,
	      struct Para_confchange *Para_confchange,
	      float *SEQID_THR,
	      int *PRINT_COV_COUPLING,
	      int *PRINT_DEF_COUPLING,
	      int *PRINT_DIR_COUPLING,
	      int *PRINT_COORD_COUPLING,
	      int *PRINT_SIGMA_DIJ,
	      char *PROF_TYPE,
	      int *ALL_PAIRS, float *SIGMA,
	      int *STRAIN, char *SITES, int *PRINT_CMAT,
	      int *ANHARMONIC, int *PRINT_PDB_ANHARM, int *N_PDB_print)
{
  FILE *file_in=fopen(filename, "r");
  if(file_in==NULL){
    printf("ERROR, TNM input file %s not found\n", filename); 
    return(0);
  }
  char string[1000], dumm[80];
  printf("Reading parameters in %s\n", filename);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "PDB1", 4)==0){
      sscanf(string+5, "%s", file_pdb1);
    }else if(strncmp(string, "CH1",3 )==0){
      sscanf(string+4, "%s", dumm);
      if((strncmp(dumm,"ALL", 3)==0)||(strncmp(dumm,"all", 3)==0)){
	*chain1='*';
      }else{
	strcpy(chain1,dumm);
      }
    }else if(strncmp(string, "PDB2", 4)==0){
      sscanf(string+5, "%s", file_pdb2);
    }else if(strncmp(string, "CH2",3 )==0){
      sscanf(string+4, "%s", dumm);
      if((strncmp(dumm,"ALL", 3)==0)||(strncmp(dumm,"all", 3)==0)){
	*chain2='*';
      }else{
	strcpy(chain2,dumm);
      }
    }else if(strncmp(string, "DOF", 3)==0){
      sscanf(string+4, "%s", dumm);
      if(strncmp(dumm, "CART", 4)==0){*ANM=1;}
      else if(strncmp(dumm, "TORS", 4)!=0){
	printf("WARNING, not allowed degrees of freedom: %s\n", dumm);
      }
    }else if(strncmp(string, "REF", 3)==0){
      sscanf(string+4, "%s", REF);
    }else if(strncmp(string, "LABEL", 5)==0){
      sscanf(string+6, "%d", LABEL);
    }else if(strncmp(string, "SIDECHAIN", 9)==0){
      sscanf(string+10, "%d", SIDECHAINS);
    }else if(strncmp(string, "MIN_INT", 7)==0){
      sscanf(string+8, "%d", MIN_INT);
    }else if(strncmp(string, "E_MIN", 5)==0){
      sscanf(string+6, "%f", E_MIN);
    }else if(strncmp(string, "COLL_THR", 8)==0){
      sscanf(string+9, "%f", COLL_THR);
      printf("COLL_THR= %.1f\n", *COLL_THR);
    }else if(strncmp(string, "PSI", 3)==0){
      sscanf(string+4, "%d", PSI);
      if(*PSI==0)printf("Psi angle frozen\n");
    }else if(strncmp(string, "OMEGA", 5)==0){
      sscanf(string+6, "%d", OMEGA);
      if(*OMEGA)printf("Omega angle used as degree of freedom\n");
      if(*OMEGA<0)
	printf("Omega angle is dof if Domega>%d degrees.\n", D_omega_thr);
    }else if(strncmp(string, "K_OMEGA", 7)==0){
      sscanf(string+8, "%lf", K_OMEGA);
      if(*OMEGA){
	printf("Elastic constant for omega angle: %.3f\n", *K_OMEGA);
      }else if(*K_OMEGA){
	printf("WARNING, elastic constant for omega set but omega is frozen\n");
      }
    }else if(strncmp(string, "K_BA", 4)==0){
      sscanf(string+5, "%lf", K_BA);
      printf("Elastic constant for bond angles: %.2f\n", *K_BA);
    }else if(strncmp(string, "K_BL", 4)==0){
      sscanf(string+5, "%lf", K_BL);
      printf("Elastic constant for bond lengths: %.2f\n", *K_BL);
    }else if(strncmp(string, "K_PSI", 5)==0){
      sscanf(string+6, "%lf", K_PSI);
      if(*K_PSI)printf("Elastic constant for psi angle:    %.3f\n", *K_PSI);
    }else if(strncmp(string, "K_PHI", 5)==0){
      sscanf(string+6, "%lf", K_PHI);
      if(*K_PHI)printf("Elastic constant for phi angle:    %.3f\n", *K_PHI);
    }else if(strncmp(string, "K_CHI", 5)==0){
      sscanf(string+6, "%lf", K_CHI);
      if(*K_PHI)printf("Elastic constant for chi angle:    %.3f\n", *K_CHI);
    }else if(strncmp(string, "CONT_TYPE", 9)==0){
      sscanf(string+10, "%s", INT_TYPE);
      if(strncmp(INT_TYPE, "HNM", 3)==0)*ANM=1;
    }else if(strncmp(string, "CONT_THR", 8)==0){
      sscanf(string+ 9, "%f", C_THR);
    }else if(strncmp(string, "S_THR", 5)==0){
      sscanf(string+6, "%f", S_THR);
    }else if(strncmp(string, "S_TYPE", 6)==0){
      sscanf(string+7, "%d", S_TYPE);
    }else if(strncmp(string, "ONEINT", 6)==0){
      sscanf(string+7, "%d", ONEINT);
    }else if(strncmp(string, "N_RESRES", 8)==0){
      sscanf(string+9, "%d", N_RESRES);
      if(*N_RESRES <1){
	printf("WARNING, not allowed value of N_RESRES: %d\n", *N_RESRES);
	printf("Setting to default value 1\n"); *N_RESRES=1;
      }
    }else if(strncmp(string, "POW", 3)==0){
      sscanf(string+4,  "%d", &POW);
    }else if(strncmp(string, "EXP_HESSIAN", 11)==0){
      sscanf(string+12, "%f", &EXP_HESSIAN);
    }else if(strncmp(string, "ANHARMONIC", 10)==0){
      sscanf(string+11, "%d", ANHARMONIC);
    }else if(strncmp(string, "KAPPA", 5)==0){
      sscanf(string+6, "%f",  KAPPA);
    }else if(strncmp(string, "FIT_B", 5)==0){
      sscanf(string+6, "%d",  FIT_B);
    }else if(strncmp(string, "RMSD_EXP", 8)==0){
      sscanf(string+9, "%f",  RMSD_EXP);
    }else if(strncmp(string, "RMSD_MIN", 8)==0){
      sscanf(string+9, "%f",  RMSD_MIN);
    }else if(strncmp(string, "PRINT_PDB_ANHARM", 16)==0){
      sscanf(string+17, "%d", PRINT_PDB_ANHARM);
    }else if(strncmp(string, "N_PDB_print", 11)==0){
      sscanf(string+12, "%d", N_PDB_print);
    }else if(strncmp(string, "NMUT", 4)==0){
      sscanf(string+5, "%d",  NMUT);
    }else if(strncmp(string, "PRED_MUT", 8)==0){
      sscanf(string+9, "%d",  PRED_MUT);
    }else if(strncmp(string, "MUT_PARA", 8)==0){
      sscanf(string+9, "%s",  Mut_para);
    }else if(strncmp(string, "RMSD_THR", 8)==0){
      sscanf(string+9, "%f",  &(Para_confchange->RMSD_THR));
    }else if(strncmp(string, "SEQID_THR", 9)==0){
      sscanf(string+10, "%f", SEQID_THR);
    }else if(strncmp(string, "DEBUG", 5)==0){
      sscanf(string+6, "%d", &DEBUG);
    }else if(strncmp(string, "FILE_FORCE", 10)==0){
      *FILE_FORCE=malloc(400*sizeof(char));
      sscanf(string+11, "%s", *FILE_FORCE);
    }else if(strncmp(string, "SIMUL", 5)==0){
      sscanf(string+6, "%d", &(Para_simul->N_SIMUL));
    }else if(strncmp(string, "STEP_MAX", 8)==0){
      sscanf(string+9, "%f", &(Para_simul->STEP_MAX));
    }else if(strncmp(string, "STEP_MIN", 8)==0){
      sscanf(string+9, "%f", &(Para_simul->STEP_MIN));
    }else if(strncmp(string, "NSTEPS", 6)==0){
      sscanf(string+7, "%d", &(Para_simul->NSTEPS));
    }else if(strncmp(string, "ANGLE", 5)==0){
      sscanf(string+6, "%f", &(Para_simul->ANGLE));
    }else if(strncmp(string, "OUTDIR", 6)==0){
      sscanf(string+7, "%s", outdir);
    }else if(strncmp(string, "NMODES", 6)==0){
      sscanf(string+7, "%d", N_MODE);
    }else if(strncmp(string, "PRINT_PDB", 9)==0){
      sscanf(string+10, "%d", PRINT_PDB);
    }else if(strncmp(string, "PRINT_CHANGE", 12)==0){
      sscanf(string+13, "%d", PRINT_CONFCHANGE);
    }else if(strncmp(string, "PRINT_FORCE", 11)==0){
      sscanf(string+12, "%d", PRINT_FORCE);
    }else if(strncmp(string, "PRINT_SUMM", 10)==0){
      sscanf(string+11, "%d", PRINT_MODE_SUMM);
    }else if(strncmp(string, "ALLOSTERY", 9)==0){
      sscanf(string+10, "%d", ALLOSTERY);
    }else if(strncmp(string, "SITES", 5)==0){
      sscanf(string+6,  "%s", SITES);
    }else if(strncmp(string, "PRINT_COV_COUPLING", 18)==0){
      sscanf(string+19, "%d", PRINT_COV_COUPLING);
    }else if(strncmp(string, "PRINT_DEF_COUPLING", 18)==0){
      sscanf(string+19, "%d", PRINT_DEF_COUPLING);
    }else if(strncmp(string, "PRINT_DIR_COUPLING", 18)==0){
      sscanf(string+19, "%d", PRINT_DIR_COUPLING);
    }else if(strncmp(string, "PRINT_COORD_COUPLING", 20)==0){
      sscanf(string+21, "%d", PRINT_COORD_COUPLING);
    }else if(strncmp(string, "PRINT_SIGMA_DIJ", 15)==0){
      sscanf(string+16, "%d", PRINT_SIGMA_DIJ);
    }else if(strncmp(string, "PRINT_CONT_MAT", 14)==0){
      sscanf(string+15, "%d", PRINT_CMAT);
    }else if(strncmp(string, "PROF_TYPE", 9)==0){
      sscanf(string+10, "%s", dumm);
      if((dumm[0]=='A')||(dumm[0]=='C')||(dumm[0]=='P')){
	*PROF_TYPE=dumm[0];
      }else{
	printf("WARNING, profile type %s does not exist\n", dumm);
      }
    }else if(strncmp(string, "ALL_PAIRS", 9)==0){
      sscanf(string+10, "%d", ALL_PAIRS);
    }else if(strncmp(string, "SIGMA", 5)==0){
      sscanf(string+6, "%f", SIGMA);
    }else if(strncmp(string, "STRAIN", 6)==0){
      sscanf(string+7, "%d", STRAIN);
    }else if(strncmp(string, "AMPLITUDE", 9)==0){
      sscanf(string+10, "%f", &(Para_simul->AMPLITUDE));
    }else if(strncmp(string, "AMPL_MIN", 8)==0){
      sscanf(string+9, "%f", &(Para_simul->AMPL_MIN));
    }else if(strncmp(string, "AMPL_MAX", 8)==0){
      sscanf(string+9, "%f", &(Para_simul->AMPL_MAX));
    }else if(strncmp(string, "AMPL_FACT", 9)==0){
      sscanf(string+10, "%f", &(Para_simul->AMPL_FACT));
    }else if(strncmp(string, "E_THR", 5)==0){
      sscanf(string+6, "%f", &(Para_simul->E_THR));
    }else if(strncmp(string, "D_REP", 5)==0){
      sscanf(string+6, "%f", &(Para_simul->D_REP));
    }else if(strncmp(string, "PDB_STEP", 8)==0){
      sscanf(string+9, "%f", &(Para_simul->PDB_STEP));

      //}else if(strncmp(string, "", )==0){
    }else if(string[0]!='\n'){
      printf("WARNING, unrecognized line:\n%s", string);
    }
  }
  fclose(file_in);
  if((*K_PHI==0)&&(*K_PSI))*K_PHI=*K_PSI;
  if((*K_PSI==0)&&(*K_PHI))*K_PSI=*K_PHI;
  if((*K_CHI==0)&&(*K_PHI))*K_CHI=*K_PHI;
  if((*S_TYPE < 0)||(*S_TYPE >2)){
    printf("WARNING, not allowed Shadow interaction type %d\n", *S_TYPE);
    *S_TYPE=0;
    printf("Setting to default %d (angle kij, d_ij<d_ik<d_jk)\n", *S_TYPE);
  }
  if(strncmp(INT_TYPE, "SHA", 3)==0){
    if(*S_TYPE==0){ // Thereshold on cosine(kij)
      if((*S_THR<0)||(*S_THR>1)){
	printf("ERROR Shadow type %d, shadow parameter %f out of range. ",
	       *S_TYPE, *S_THR); exit(8);
      }
    }else if (*S_TYPE==1){ // Threshold on distance from axis
      if(*S_THR<0){
	printf("ERROR Shadow type %d, shadow parameter %f out of range. ",
	       *S_TYPE, *S_THR); exit(8);
      }
    }else if (*S_TYPE==2){ // Atom radius
      if(*S_THR<0){
	printf("ERROR Shadow type %d, shadow parameter %f out of range. ",
	       *S_TYPE, *S_THR); exit(8);
      }
    }
  }
  return(1);
}

int Relevant_modes(float *sigma2, int *select, int N, float N_atoms)
{
  double sum=0; int i;
  for(i=0; i<N; i++)if(select[i])sum+=sigma2[i];
  float eps; if(N_atoms>100){eps=1./N_atoms;}else{eps=1./100;}
  float thr=sum*(1.-eps); sum=0;
  for(i=0; i<N; i++){
    if(select[i])sum+=sigma2[i];
    if(sum>thr)break;
  }
  return(i);
}


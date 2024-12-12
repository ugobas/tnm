
/*
  Units used in the TNM program:
  energy kBT= 1 -> kBT at 293K = 1.381 e-23 *293= 4.05e-21 J
  length = 1A   -> e-10 m
  mass = atomic mass -> 1.66e-27 kg
  time unit =  0.64 e-13 = 6.4 e-12 = 6.4 ps
  frequency = Freq_unit = 1/time_unit = 0.156 ps^(-1)
  frequency omega_q = kBT/h/ = (0.1636 ps)^(-1)= 6.1 ps^(-1) = 39 Freq_unit
  1/(omega_q)^2= 0.00065 in internal units

  Normal modes:
  omega^2 represents both frequency^2 in Freq_unit^2 and
  omega^2/RT in units of mass^-1 length^-2
  (RT/M*omega^2) has units of length
  kinetic energy matrix T has dimension m*l^2
  Mass_weighted torsional eigenvector u^alpha has dimension 0
  Torsional eigenvector v^alpha = T^(-1/2)u^alpha has dimension m^(-1/2)l^(-1)
  Cartesian eigenvector x^alpha = Jv^alpha has dimension m^(-1/2)
  Cartesian eigenvector r^alpha = m^(1/2)x^alpha has dimension 0 
  displacement d_i=c_alpha x^\alpha has dimension of length, with
  dimension(c_alpha)=m^(1/2)*l
  <c_alpha>^2 = (RT/omega_alpha^2)
  c_alpha*x^alpha has dimension l
  c_alpha*r^alpha has dimension m^(1/2)*l
  c_alpha*v^alpha has dimension 0 (torsion angle)
  c_alpha*u^alpha has dimension m^(1/2)*l
  MSD = sum_i m_i*|d_i|^2/M = sum_alpha(c_alpha)^2/M
  = sum_alpha(RT/M*omega_alpha^2)
 */

#define PRG "tnm"
int LABEL=0; // Print model parameters in file name?

char DOF_LABEL[30];
// Free energy variables global
double G_NM_holo=0;
double DG_NM_internal=0;
double DG_NM_rigid=0;
double E_cont_bind=0;
int N_rigid_inter;

/********** Some of these parameters are not set from input file **********/
// Force constant
//float KAPPA_DEF=218.4; // Force const. if B not fitted
float KAPPA_DEF_PHI=5; //15; // Force const. if B not fitted
float KAPPA_DEF_PHIPSI=8.2;  //24.5; // Force const. if B not fit
float KAPPA_DEF_OMEGA=15;   // 69.2; // Force const. if B not fit
float KAPPA_DEF_SCHAIN=14;  //42.1; // Force const. if B not fit
// In all above cases, C=4.5 E=6 K_PHI=K_PSI=K_OM=K_SCH=0.2 T=293K r0=3.5A

float TEMPERATURE=-1;
float Temp_comp=-1;
float KAPPA_DEF;
float TEMP_DEF=293;
float Freq_unit=0.156; // internal units in units of ps^(-1)

// Selection of modes
#define COLL_THR_DEF 30   
float COLL_THR= COLL_THR_DEF;  // Select mode if exp(-S_cart) > COLL_THR
double E_THR;
// E_MIN

// Analysis of conformation change
int SWITCH_BONDS=1; // Select only bonds with aligned atoms (1)
int MIN_ALI_LEN_DEF=20; // Chains are aligned if L>MIN_ALI_LEN, if not peptides
float RMSD_MAX=100;   // max. RMSD for analyzing conformation change
float RMSD_MIN=0.5;   // min. RMSD for analyzing conformation change
#define RMSD_THR_DEF 0.1  // If confchange projections smaller, discard
#define SEQID_THR_DEF 0.3  // If Seq.Id<THR, align structures with Mammoth
int NMUT=-2; /* Mutations between the two proteins required for analysis.
		NMUT=-2: No analysis. NMUT=-1: Always analysis.
		NMUT=0: Sequences must be identical.
		NMUT=1: Exactly one mutation */
int PRED_MUT=0; // Predict mutations?
int IWT=0; // Average rmsd of WT mut over the lowest IWT ones
//int OPT_MUT=0; // Optimize coefficients to fit RMSD of mutation?
char Mut_para[100]="Mutation_para.in";

// Simulation of conformation changes
#define AMPLITUDE_DEF    1.0  // NM amplitude for simulated structures
#define E_THR_DEF        20.0 // Threshold for stopping motion
#define E_THR1_DEF       1.0  // Threshold for deciding privileged direction
#define D_REP_DEF        2.5  // Threshold for repulsions
#define MAX_ANGLE_DEF    0.4  // Max. angle for torsional deformations, radiants
float s0=0.07;       // Entropy per residue, used for simulations
float cont_inf=3.9; // Asymptotic value of number of contacts
int N_int_long;

// Unfolding
float K_tors_unfold=0.2; // Use torsional springs in unfolding calculations

// Analysis performed
int ANHARMONIC=0;
int SITE_DYNAMICS=0; // Examine dynamics of binding sites

// Output
#define PRINT_PDB_DEF 1
#define PRINT_MODE_SUMM_DEF 1
#define PRINT_FORCE_DEF 0
#define PDB_STEP_DEF  0.4  // Step for printing PDB


// Dynamical couplings
int PRINT_DEF_COUPLING=1;  // Print deformation coupling?
int PRINT_DIR_COUPLING=1;  // Print directionality coupling?
int PRINT_COV_COUPLING=1;  // Print covariance coupling?
int PRINT_COORD_COUPLING=1; // Print coordination coupling?
int PRINT_SIGMA_DIJ=0;  // Print variance of DIJ distanc
int PRINT_CMAT=0;       // Print contact matrix
int STRAIN=0; // Compute strain profile, similar to PNAS 106:14253 2009 ?
int ALL_PAIRS=0; // Print couplings for all pairs?
float SIGMA=1.5; // Print only couplings > SIGMA*std.dev.
char PROF_TYPE='C';
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
#define K_OMEGA_DEF 0.5    // Force constant for omega torsions
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
//#include "unfolding.h"
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
double nbtops = CLOCKS_PER_SEC;

// External variables
char *AA_code;
int HYD=0;
int HYD_REF=0;

// atoms
int num_BB=5; //=5
char BB_name[8][3]; // N CA C CB O
//                       C-N    N-CA   CA-C   CA-CB  C-O  
float BB_length[8]; //={1.3301,1.4572,1.5236,1.5319,1.2316};
//                    CA-C-N C-N-CA N-CA-C N-CA-CB CA-C-O  
float BB_angle[8]; //={1.1091,1.0196,1.2018,1.2162,1.0359};

// Computation of native interactions and force constant in
// interactions_tnm.c and anharmonic_tnm.c
// k_ij = KAPPA*(r_ij/C0)^-EXP_HESSIAN if r<C1 0 if r>=C_THR (if POW)
int POW=1;      // Power law (1) or exponential (0)?
float C_THR=0;  // Threshold for contacts
float EXP_HESSIAN=0; 
float C1_THR=0; // Start of interpolation region
float C0=0;     // Reference distance for force constant
float KAPPA=0;  // Force constant at r=C0
double K0=0; // K0=KAPPA*C0^EXP;
double K1=0; // K1=K0*C1^-EXP/(C_THR-C1);
double K_OMEGA=0, K_PSI=0, K_PHI=0, K_CHI=0, K_BA=0, K_BL=0;
// Force constants for internal variables
int D_omega_thr=0; // Omega unfrozen if D_omega>D_omega_thr

// Fit of B factors
int FIT_B=0;
float R_MIN=0;     // Minimum correlation above which B-factors are fitted
float SLOPE_MIN=0; // Minimum slope above which B-factors are fitted
float RMSD_EXP=0;

// Energy function
float **Econt; 
float lambda=0;
float E_repulsion=0;
float D_Res[20][20], U_Res[20][20];
float ENE_HB=2.0; // Ratio between HB interactions and other types
float DA_THR=3.5; // Max donor-acceptor dist in h bond
float COS_DAA=-0.008;   // Max cos(A'AD) (min. 105 degrees)
float THR_NO=4;   // Threshold for NO interactions
int S_TYPE=0;        // Type of shadow interaction
int ONEINT=0;        // Only one interaction per atom vs. residue
float S_THR=0;       // Screening parameter
int HNM=0;           // Calpha interactions with Hissen parameters
int NOCOV=1;         // Covalent neighbors do not interact
int N_RESRES=0;      // Number of interactions retained for each residue pair

// General
char REF[20];      // Reference atoms
int N_MODES=0;
float E_MIN=0;  // Minimum eigenvalue allowed for normal modes
int KINETIC=1;      // Considering kinetic energy in tnm model
int MIN_INT_MAIN=MIN_INT_DEF; // Min. numb. of interactions per degree of freedom

// Input
int L_MAX;
int MIN_ALI_LEN;  // Chains are aligned if L>MIN_ALI_LEN, if not peptides


// Output
int Verbose=0;
int PRINT_LAMBDA=0; // Printing list of eigenvalues
int ANISOU=0;     // Examine anisotropic temperature factors?
int PRINT_AXES=0;
int N_PDB_print=0;
int PRINT_PDB_ANHARM=0;
char REF_CC[5];

// Other
float sqrarg=0;

extern float Fit_fluctuations(int FIT_B, float *B_TNM,
			      float *RMSD_NM, char *out,
			      float *Temp_comp, float TEMPERATURE,
			      int anharmonic, struct Normal_Mode NM,
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

double Print_summary(char *name_out, char *name1, char *parameters, int Nchain,
		     struct Normal_Mode NM, char *out_B1, char *out_B2,
		     float Tors_fluct,struct axe *axes, int N_axes, int N_main,
		     int Nskip, int N_atoms, int N_res, int N_diso, int N_int,
		     int N_inter, float Cont_overlap,
		     int N_disc_coll, int N_disc_freq, int N_disc_out,
		     float mass_sum,
		     struct Normal_Mode NM_unfold,
		     int N_int_unfold, double E_cont_unfold, int FOLDING,
		     struct Normal_Mode NM_apo, int N_int_inter,
		     double E_cont_holo, double E_cont_bind,
		     int BINDING,
		     float *mass_coord, int *Cart_interface, int N_Cart_inter,
		     int *rigid_inter, int N_rigid);
double Print_thermal(int anharmonic, FILE *file_out,
		     struct Normal_Mode NM, char *out_B,
		     struct axe *axes, int N_main, int N_res, int Nchain,
		     struct Normal_Mode NM_unfold,
		     int N_int_unfold, int N_int,
		     double E_cont_unfold, int FOLDING,
		     struct Normal_Mode NM_apo, int N_int_inter,
		     double E_cont_bind, int BINDING,
		     float *mass_coord, int *Cart_interface, int N_Cart_inter,
		     int *rigid_inter, int N_rigid);
double Compute_entropy(double *Free_Ene, double *Free_Ene_rigid,
		       struct Normal_Mode NM, float *sigma2,
		       double beta_omega, float *mass_coord,
		       int *Cart_interface, int N_Cart_inter,
		       int *rigid_inter, int N_rigid);

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
	    float *RMSD_MIN, int *NMUT,
	    int *PRED_MUT, int *IWT, char *Mut_para,
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
	    int *PRINT_PDB_ANHARM, int *N_PDB_print,
	    int *FOLDING, int *BINDING, int *UNFOLDING,
	    int *n_apo, char *apo_chains);

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
	      float *RMSD_MIN, int *NMUT,
	      int *PRED_MUT, int *IWT, char *Mut_para,
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
	      int *PRINT_PDB_ANHARM, int *N_PDB_print,
	      int *FOLDING, int *BINDING, int *UNFOLDING,
	      int *n_apo, char *apo_chains);
void help (void);

/******************************* Computations ****************************/
int Torsional_Hessian(double **Hessian, int N_modes,
		      struct axe *axe,
		      float K_OMEGA, float K_PSI, float K_PHI, float K_CHI,
		      float K_BA, float K_BL);
double Rescale_Hessian(struct Normal_Mode NM,
		       struct interaction *Int_list, int N_int,
		       float kappa, int anharmonic);
void Rescale_para(float kappa);

void Get_modes_tors(struct Normal_Mode NM, int N_modes,
		    double **Hessian, struct Jacobian J,
		    struct axe *axe1, int naxe1,
		    struct Reference Ref_kin, atom *atoms1);
void Get_modes_Cart(struct Normal_Mode NM, int N_modes,
		    double **Hessian, struct Jacobian J,
		    struct axe *axe1, int naxe1, struct Reference Ref_kin);
int Compute_sigma2(int *N_disc_freq, int *N_disc_out, struct Normal_Mode NM,
		   float *mass_sqrt, int naxe1, int N_cart,
		   float COLL_THR);
struct interaction *Apo_interactions(int *N_int_apo,
				     struct interaction **Int_list_inter,
				     int *N_int_inter,
				     struct interaction *Int_list, int N_int,
				     atom *atoms,int *apo_chain_num,
				     int n_apo, int Nchain);
struct interaction *Unfold_interactions(int *N_int_unfold,
					struct interaction *Int_list,
					int N_int, atom *atoms);
int N_res_1=0, N_res_2=0, N_cont_long_1=0, N_cont_long_2=0;
int Get_interface(int *Cart_interface,
		  struct interaction *Int_list, int N_int,
		  atom *atoms, int *apo_chain_num,
		  int n_apo, struct chain *chains, int Nchain,
		  int *atom_num, int N_ref, int natoms, int nres);
int Get_chain_num(int *apo_chain_num, char *apo_chains, int n_apo,
		  struct chain *chains, int Nchain);
int Belong_to_apo(int chain, int *apo_chain_num, int n_apo);

extern void Binding_site_dynamics(struct Normal_Mode NM, struct Reference Ref,
				  struct residue *seq, atom *atoms,
				  char *nameout1,
				  char *pdb, char *chain, int Nchain,
				  int Nres, char *SITES);

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
  char  file_pdb1[150], chain1[CHMAX], pdbid1[100];
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
  char  file_pdb2[150], chain2[CHMAX], pdbid2[100];
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
  float SEQID_THR=SEQID_THR_DEF;

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
  double t0=clock(), t_all=t0;

  // File names
  #define NCH 150
  char outdir[50]="", nameout1[300]="", summary1[400]=""; //directory[150],
  char fullmodel[100]="", MODEL[80]="";
  char name2[201]="", nameout2[300]="", summary2[400]="";
  char nameprot1[100], nameprot2[100];
  char parameters[1000]="";

  /*****************************  INPUT  **********************************/
  int PRINT_CONFCHANGE=0;
  int PRINT_FORCE=PRINT_FORCE_DEF;
  int PRINT_PDB= PRINT_PDB_DEF;
  int PRINT_MODE_SUMM=PRINT_MODE_SUMM_DEF;

  /***************** Binding free energy ********************************/
  int BINDING=0, n_apo=1; char apo_chains[100]; 
  /**************** Unfolding entropy ***********************************/
  int FOLDING=0;
  int UNFOLDING=0;

  outdir[0]='\0';
  ALL_AXES=1;
  KINETIC=1;
  E_MIN=0.00000001; //E_MIN_DEF;
  ONEINT=1;
  N_RESRES=1;
  HYD=0; HYD_REF=0; // Read hydrogen atoms or not?
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
	  &NMUT, &PRED_MUT, &IWT, Mut_para,
	  &Para_simul, &Para_confchange, &SEQID_THR, parameters,
	  &PRINT_COV_COUPLING, &PRINT_DEF_COUPLING,
	  &PRINT_DIR_COUPLING, &PRINT_COORD_COUPLING,
	  &PRINT_SIGMA_DIJ, &PROF_TYPE, &ALL_PAIRS, &SIGMA, &STRAIN, SITES,
	  &PRINT_CMAT, &ANHARMONIC, &PRINT_PDB_ANHARM, &N_PDB_print,
	  &FOLDING, &BINDING, &UNFOLDING, &n_apo, apo_chains);
  MIN_INT_MAIN=MIN_INT_SIDE;
  if(BINDING){
    if(ALLOSTERY==0){
      printf("WARNING, dynamical couplings are needed for "
	     "estimating binding affinity, setting it.\n");
    }
    ALLOSTERY=1; //PRINT_DEF_COUPLING=1;
    PRINT_DIR_COUPLING=1; PRINT_COORD_COUPLING=1;
    if(FIT_B==0 && KAPPA){
      printf("WARNING, Binding computation must be run with default KAPPA\n"); 
      KAPPA=0;
    }
  }
  if(FOLDING){
    if(FIT_B){
      printf("WARNING, Unfolding computation must be run with FITB=0\n");
      FIT_B=0;
    }
    if(KAPPA){
      printf("WARNING, Unfolding computation must be run with default KAPPA\n");
	     KAPPA=0;
    }
  }

  sprintf(DOF_LABEL,"PHI");
  if(PSI)strcat(DOF_LABEL, "PSI");
  if(OMEGA>0)strcat(DOF_LABEL, "OME");
  if(OMEGA<0)strcat(DOF_LABEL, "CISTR");
  if(SIDECHAINS)strcat(DOF_LABEL, "SCH");

  // Rescale torsion springs
  if(PSI==0 && OMEGA==0 && SIDECHAINS==0){
    KAPPA_DEF=KAPPA_DEF_PHI; 
  }else if(PSI && OMEGA==0 && SIDECHAINS==0){
    KAPPA_DEF=KAPPA_DEF_PHIPSI;
  }else if(PSI && OMEGA && SIDECHAINS==0){
    KAPPA_DEF=KAPPA_DEF_OMEGA;
  }else if(PSI && OMEGA && SIDECHAINS){
    KAPPA_DEF=KAPPA_DEF_SCHAIN;
  }else{
    printf("ERROR, not recommended set of degrees of freedom\n");
    printf("PSI= %d OMEGA= %d SIDECHAINS= %d\n",
	   PSI,OMEGA,SIDECHAINS); exit(8);
  }
  if(KAPPA==0){
    KAPPA=KAPPA_DEF;
    printf("Setting KAPPA to default %.4g\n", KAPPA);
  }

  /*char tmp[20]; sprintf(tmp, " KAPPA=%.1f", KAPPA);
    strcat(parameters,tmp);
  K_PSI*=KAPPA; K_PHI*=KAPPA; K_OMEGA*=KAPPA; K_CHI*=KAPPA;
  K_BA*=KAPPA; K_BL*=KAPPA;*/

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
	   &TEMPERATURE, pdbid1, file_pdb1);
  Nchain=Nchain1;

  //N_TNM=naxe1;
  

  /*******************  Protein structure 2, if any  ********************/
  /* Read protein structure 2 */
  int Nmut=0, *Posmut=NULL; char *AAmut=NULL, *AAwt=NULL;
  char *str_mut=NULL;
  if(file_pdb2[0]!='\0'){
  
    printf("\n************************************************\n");
    printf("*********** Read/align structure 2 *************\n");
    printf("************************************************\n\n");
    float TEMP2;
    Read_PDB(&nres2, &seq2, &n_lig2, &ligres2, chain2, &ANISOU2,
	     &nmr2, &natoms2, &atoms2, &na_lig2, &ligatom2,
	     &chains2, &Nchain2, &TEMP2, pdbid2, file_pdb2);
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
    str_mut=malloc(8*Nmut*sizeof(char));
    strcpy(str_mut, "");
    for(i=0; i<Nmut; i++){
      char mmut[8];
      sprintf(mmut,"%c%d%c ", AAwt[i], Posmut[i], AAmut[i]);
      strcat(str_mut, mmut);
    }
    printf("%d mutations found: %s\n", Nmut, str_mut);
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
  int N_Cart=3*N_ref;
  // Check that reference atoms exist
  if((natoms1==0)||(N_ref==0)){
    printf("ERROR, no atoms found in file %s", file_pdb1);
    printf(" natoms= %d N_ref= %d\n", natoms1, N_ref);
    exit(8);
  }
  printf("%d reference atoms\n", N_ref);
  // Mass of the protein
  mass_sum=Ref_kin.mass_sum;

  {
    char tmp[20]; sprintf(tmp, " KAPPA=%.1f", KAPPA);
    strcat(parameters,tmp);
    K_PSI*=KAPPA; K_PHI*=KAPPA; K_OMEGA*=KAPPA; K_CHI*=KAPPA;
    K_BA*=KAPPA; K_BL*=KAPPA;
  }

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
  char tmp[120];
  sprintf(nameprot1, "%s", pdbid1);
  for(i=0; i<Nchain; i++){
    if((chains1[i].label!=' ')&&(chains1[i].label!='\0')){
      sprintf(tmp, "%c", chains1[i].label);
      strcat(nameprot1, tmp);
    }
  }
  strcpy(nameout1, nameprot1);

  if(strncmp(INT_TYPE, "SCR", 3)==0){
    if(ONEINT==0){
      sprintf(tmp, "%s_d%.2f_t%.1f", INT_TYPE, S_THR, C_THR);
    }else{
      sprintf(tmp, "%s_MIN_d%.2f_t%.1f", INT_TYPE, S_THR, C_THR);
    }
  }else if(strncmp(INT_TYPE, "SHA", 3)==0){
    if(ONEINT==0){
      sprintf(tmp, "%s_%d_d%.2f_t%.1f", INT_TYPE, S_TYPE, S_THR, C_THR);
    }else{
      sprintf(tmp, "%s_MIN_%d_d%.2f_t%.1f",INT_TYPE, S_TYPE, S_THR, C_THR);
    }
  }else if((strncmp(INT_TYPE, "MIN", 3)==0)&&(N_RESRES>1)){
    sprintf(tmp, "%s%.1f_M%d", INT_TYPE, C_THR, N_RESRES);	      
  }else{
    sprintf(tmp, "%s%.1f", INT_TYPE, C_THR);
  }
  strcat(fullmodel, tmp);

  strcat(fullmodel,"_"); strcat(fullmodel,REF);
  strcat(fullmodel,"_"); strcat(fullmodel,DOF_LABEL);
  if(PSI)strcat(fullmodel, "PSI");
  if(OMEGA)strcat(fullmodel, "OME");
  if(OMEGA<0)strcat(fullmodel, "CISTR");
  if(SIDECHAINS)strcat(fullmodel, "SCH");
  if(LABEL){strcat(nameout1,"_"); strcat(nameout1,fullmodel);}
  // UUU: Restore the option of the full name

  sprintf(summary1, "%s.summary.dat", nameout1);   
  if(Check_make_dir(outdir)){
    char tmp2[400]; sprintf(tmp2, "%s", summary1);
    sprintf(summary1, "%s/", outdir); strcat(summary1, tmp2);
    sprintf(tmp2, "%s", nameout1);
    sprintf(nameout1, "%s/", outdir); strcat(nameout1, tmp2);
  }
  printf("output file (nameout1): %s\n",nameout1);
  printf("output file (summary1): %s\n",summary1);

  if(CONF_CHANGE){
    sprintf(nameprot2, "%s", pdbid2);
    for(i=0; i<Nchain; i++){
      if((chains2[i].label!=' ')&&(chains2[i].label!='\0')){
	sprintf(tmp, "%c", chains2[i].label);
	strcat(nameprot2, tmp);
      }
    }
    //sprintf(name2, "%s_%s", nameprot1, nameprot2);
    //sprintf(name2, "%s-%s", nameprot1, nameprot2);
    //sprintf(name2, "%s_%s", pdbid1, pdbid2); // YYY don't write chain names 

    sprintf(nameout2, "%s_%s", nameprot1, nameprot2); 
    if(LABEL){strcat(nameout2,"_"); strcat(nameout2, fullmodel);}

    sprintf(summary2, "%s.summary.dat", nameout2);
    if(Check_make_dir(outdir)){
      char tmp2[400]; sprintf(tmp2, "%s", summary2);
      sprintf(summary2, "%s/", outdir); strcat(summary2, tmp2);
      sprintf(tmp2, "%s", nameout2);
      sprintf(nameout2, "%s/", outdir); strcat(nameout2, tmp2);
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
    int N_Cart=ali_atoms.N_Cart;
    coord_1=malloc(N_Cart*sizeof(float));
    Write_ref_coord_atom(coord_1, N_ali, atoms1, ali_atoms.ali1);
    coord_2=malloc(N_Cart*sizeof(float));
    Write_ref_coord_atom(coord_2, N_ali, atoms2, ali_atoms.ali2);
    float rmsd=RMSD_CA(coord_1, coord_2, Ref1, atoms1, nres1, ali_atoms);
    printf("RMSD_CA(%s,%s)= %.2f\n", pdbid1, pdbid2, rmsd);
    if((rmsd < RMSD_MIN)||(rmsd > RMSD_MAX)){
      char namefile[100]; sprintf(namefile, "%s.RMSD", nameout2);
      FILE *file_out=fopen(namefile, "w"); char out[400], tmp[80];
      sprintf(out, "%d mutations", Nmut);
      if(str_mut){sprintf(tmp, ": %s", str_mut); strcat(out, tmp);}
      sprintf(out, "\nRMSD= %.3f\n"
	      "Too small RMSD< %.2f\nExiting the program\n",rmsd, RMSD_MIN);
      strcat(out, tmp);
      fprintf(file_out, "%s", out);
      fclose(file_out);
      printf("%sRMSD printed in %s\n", out, namefile);
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
  if(ANM)NOCOV=0;

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
  struct rigid *Rigid_dof=NULL; int N_rigid=0;
  axe1=Set_DegofFreed(&naxe1, &nmain, &nskip, &N_diso1, &Rigid_dof, &N_rigid,
		      bonds, atoms1, natoms1, Ref_kin.atom_num, N_ref,
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
  if(ANM){N_modes=N_Cart;}else{N_modes=naxe1;}
  Allocate_memory(&NM, &J, &Hessian, naxe1, N_ref, N_Cart, N_modes, nres1);


  // Unfolding entropy
  struct Normal_Mode NM_unfold; double **Hessian_unfold=NULL;
  struct interaction *Int_list_unfold=NULL; int N_int_unfold=0;
  double E_cont_holo=Contact_energy(Int_list, N_int, seq1), E_cont_unfold=0;
  if(FOLDING){
    Int_list_unfold=Unfold_interactions(&N_int_unfold, Int_list, N_int,atoms1);
    E_cont_unfold= Contact_energy(Int_list_unfold, N_int_unfold, seq1);
    printf("Found %d contacts in folded and %d in unfolded str.\n",
	   N_int,  N_int_unfold);
    printf("Contact energy: %.3g (folded) %.3g (unfolded)\n",
	   E_cont_holo, E_cont_unfold);
    Allocate_Normal_modes(&NM_unfold, N_modes, naxe1, N_Cart);
    NM_unfold.ANM=ANM; NM_unfold.N=N_modes;
    Hessian_unfold=Allocate_mat2_d(N_modes, N_modes);
  }

  // Interface
  double E_cont_apo;
  struct Normal_Mode NM_apo; double **Hessian_apo=NULL;
  struct interaction *Int_list_apo=NULL; int N_int_apo=0;
  struct interaction *Int_list_inter=NULL; int N_int_inter=0;
  int Cart_interface[N_Cart], N_Cart_inter=0;
  int rigid_inter[N_rigid]; for(i=0; i<N_rigid; i++)rigid_inter[i]=-1;
  if(BINDING){
    if(Nchain1<=n_apo){
      printf("WARNING, you requested to compute the binding free energy "
	     "of %d chains with the others but we found only %d chains\n",
	     n_apo, Nchain1); BINDING=0; goto No_bind;
    }
    int apo_chain_num[n_apo];
    if(Get_chain_num(apo_chain_num, apo_chains, n_apo, chains1, Nchain1)){
      Int_list_apo=Apo_interactions(&N_int_apo, &Int_list_inter, &N_int_inter,
				    Int_list,N_int,
				    atoms1,apo_chain_num,n_apo,Nchain1);
      N_Cart_inter=Get_interface(Cart_interface, Int_list,N_int,
				 atoms1,apo_chain_num,n_apo,chains1,Nchain1,
				 Ref_kin.atom_num, N_ref, natoms1, nres1);
      N_rigid_inter=0;
      for(i=0; i<N_rigid; i++){
	int apo1=Belong_to_apo(Rigid_dof[i].chain1, apo_chain_num, n_apo);
	int apo2=Belong_to_apo(Rigid_dof[i].chain2, apo_chain_num, n_apo);
	if(apo1 != apo2){rigid_inter[i]=Rigid_dof[i].iaxe; N_rigid_inter++;}
      }

      if(N_int_apo==N_int){
	printf("WARNING, no contacts found between %s and the other chains\n"
	       "I cannot compute binding free energy\n",
	       apo_chains);
	BINDING=0; goto No_bind;
      }
      Allocate_Normal_modes(&NM_apo, N_modes, naxe1, N_Cart);
      NM_apo.ANM=ANM; NM_apo.N=N_modes;
      Hessian_apo=Allocate_mat2_d(N_modes, N_modes);
      E_cont_apo= Contact_energy(Int_list_apo, N_int_apo, seq1);
      E_cont_bind=E_cont_holo-E_cont_apo; //N_int_inter=N_int-N_int_apo;
      printf("Found %d contacts in holo str. and %d in apo str.\n",
	      N_int,  N_int_apo);
      printf("Contact energy: %.3g (holo) %.3g (apo) %.3g (bind)\n",
	     E_cont_holo, E_cont_apo, E_cont_bind);
      printf("Computing binding free energy with normal modes\n");
    }else{
      printf("WARNING, not all chains in %s were found in PDB file. "
	     "I cannot perform binding free energy computation\n",
	     apo_chains); BINDING=0;
    }
  }
 No_bind:
  printf("Memory allocated\n");

  printf("\n************************************************\n");
  printf("*********** ENM kinetic en, Eckart cdt *********\n");
  printf("************************************************\n\n");


  /********* Eckart conditions, Jacobian, Kinetic energy  ************/
  J.N_kin=
    Compute_kinetic(&J, axe1, naxe1, atoms1, natoms1, Ref_kin, 1);
  if((ANM==0)&&(KINETIC)){
    N_modes=J.N_kin; NM.N=N_modes; NM_apo.N=N_modes; NM_unfold.N=N_modes;
  }
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
    Get_modes_tors(NM, N_modes, Hessian, J, axe1, naxe1, Ref_kin, atoms1);

    printf("Normal modes computation. Time= %.2lf sec.\n",(clock()-t0)/nbtops);
    t0=clock();
    if(FOLDING){
      Compute_Hessian_TNM(Hessian_unfold, J.T_sqrt, J.T_sqrt_inv,
			  Int_list_unfold, N_int_unfold,
			  atoms1, natoms1, axe1, nmain, naxe1, N_modes,
			  chains1, Nchain, KINETIC,
			  K_OMEGA, K_PSI, K_PHI, K_CHI, K_BA, K_BL);
      printf("Hessian unfold computation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
      Get_modes_tors(NM_unfold, N_modes, Hessian_unfold, J, axe1, naxe1,
		     Ref_kin, atoms1);
      printf("Normal modes unfold computation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
    }
    if(BINDING){
      Compute_Hessian_TNM(Hessian_apo, J.T_sqrt, J.T_sqrt_inv,
			  Int_list_apo, N_int_apo,
			  atoms1, natoms1, axe1, nmain, naxe1, N_modes,
			  chains1, Nchain, KINETIC,
			  K_OMEGA, K_PSI, K_PHI, K_CHI, K_BA, K_BL);
      printf("Hessian apo computation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
      Get_modes_tors(NM_apo, N_modes, Hessian_apo,
		     J, axe1, naxe1, Ref_kin, atoms1);
      printf("Normal modes apo computation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
    }

  /*************************  Normal modes ANM *************************/
  }else{ // ANM
    if(HNM){printf("Hinsen network model");}
    else{printf("Anisotropic network model\n");}

    t0=clock();
    Compute_Hessian_ANM(Hessian, N_ref, Ref_kin.atom_num,
			Int_list, N_int, atoms1, natoms1);
    printf("Hessian computation. Time= %.2lf sec.\n",(clock()-t0)/nbtops);
    t0=clock();
    Get_modes_Cart(NM, N_modes, Hessian, J, axe1, naxe1, Ref_kin);
    printf("ANM modes computation. Time= %.2lf sec.\n",(clock()-t0)/nbtops);
    t0=clock();
    if(FOLDING){
      Compute_Hessian_ANM(Hessian_unfold, N_ref, Ref_kin.atom_num,
			  Int_list_unfold, N_int_unfold, atoms1, natoms1);
      printf("Hessian unfold computation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
      Get_modes_Cart(NM_unfold,N_modes,Hessian_unfold,J,axe1,naxe1,Ref_kin);
      printf("Normal modes unfold computation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
    }
    if(BINDING){
      Compute_Hessian_ANM(Hessian_apo, N_ref, Ref_kin.atom_num,
			  Int_list_apo, N_int_apo, atoms1, natoms1);
      printf("Hessian apo computation. Time= %.2lf sec.\n",(clock()-t0)/nbtops);
      t0=clock();
      Get_modes_Cart(NM_apo, N_modes, Hessian_apo, J, axe1, naxe1, Ref_kin);
      printf("Normal modes apo computation. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
    }
  }
  printf("Normal mode computation finished\n");


  printf("\n************************************************\n");
  printf("*********** ENM normal mode coll. **************\n");
  printf("************************************************\n\n");



  /********************  Collectivity of normal modes *******************/
  // Select modes based on eigenvalues and collectivity

  // Collectivity
  if(COLL_THR > natoms1*0.05)COLL_THR = natoms1*0.05;
  COLL_THR *=3;

  double E_ave=0; for(ia=0; ia<NM.N; ia++)E_ave+=NM.omega2[ia]; E_ave/=NM.N;
  printf("Average omega^2 prior to norm: %.2g\n", E_ave);
  E_THR=E_MIN*E_ave;

  if(1){ // YYY additional output in log file (related to discarded modes)
    printf("(discard) AV_OME2     %10.4g   | average EV (or omega^2) "
	   "before kappa fit\n", E_ave);
    double E_logave=0.0;
    double E_xmin=pow(10,10);
    double E_xmax=(-1)*pow(10,10);
    double E_normsig2=0.0;
    for(ia=0; ia<NM.N; ia++){
      E_logave+=log(NM.omega2[ia]);
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
    printf("(discard) N_Cart      %10d   | ", N_Cart);
    printf("= natoms * 3 --> collect = xxx/N_Cart\n");
  }


  // Select modes based on frequency collectivity and outliers
  printf("Selecting modes of full interaction matrix\n");
  int N_disc_freq=0, N_disc_out=0, N_disc_coll=
    Compute_sigma2(&N_disc_freq, &N_disc_out, NM,
		   Ref_kin.mass_sqrt, naxe1, N_Cart, COLL_THR);

  if(FOLDING){
    int N_disc_unfold_f=0, N_disc_unfold_o=0;
    printf("Selecting modes of unfolded interaction matrix\n");
    Compute_sigma2(&N_disc_unfold_f, &N_disc_unfold_o, NM_unfold,
		   Ref_kin.mass_sqrt, naxe1, N_Cart, COLL_THR);
  }
  if(BINDING){
    int N_disc_apo_f=0, N_disc_apo_o=0;
    printf("Selecting modes of unbound interaction matrix\n");
    Compute_sigma2(&N_disc_apo_f, &N_disc_apo_o, NM_apo,
		   Ref_kin.mass_sqrt, naxe1, N_Cart, COLL_THR);
  }

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
	for(i=0; i<N_Cart; i++)q+=Ca[i]*Cb[i]*Ref_kin.mass_coord[i];
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
    Fit_fluctuations(FIT_B, B_TNM[0], &RMSD_NM, out_B1, &Temp_comp,
		     TEMPERATURE, anharmonic, NM, Ref_kin, REF, seq1, nres1,
		     nmr1, atoms1, nameout1, &factor_B, outlier_B);

  // Rescale force constant
  Rescale_para(kappa);
  double sum_sigma2= Rescale_Hessian(NM, Int_list, N_int, kappa, anharmonic);
  if(FOLDING){
    Rescale_Hessian(NM_unfold, Int_list_unfold, N_int_unfold, kappa,anharmonic);
  }
  if(BINDING){
    Rescale_Hessian(NM_apo, Int_list_apo, N_int_apo, kappa, anharmonic);
  }


  printf("Fitting B factors. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
  printf("(discard) KAPPA   %10.4g    | "
	 "ratio between OME2 above and in .Modes file\n\n",kappa);
  // YYY added output in log file
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
      if(ix<0)break;
      sum+=NM.sigma2[ix]; sum_anh+=NM.sigma2_anhar[ix];
    }
    if(sum>0){
      sum_anh/=sum;
      for(ix=i; ix<NM.N; ix++)NM.sigma2_anhar[ix]=NM.sigma2[ix]*sum_anh;
    }
    anharmonic=1;
    B_TNM[1]=malloc(N_ref*sizeof(float));
    kappa=Fit_fluctuations(FIT_B, B_TNM[1], &RMSD_NM, out_B2,
			   &Temp_comp, TEMPERATURE, anharmonic, NM, Ref_kin,
			   REF, seq1, nres1, nmr1,
			   atoms1, nameout1, &factor_B, outlier_B);

    // Rescale force constant
    Rescale_para(kappa);
    //double sum_sigma2_anhar=
    Rescale_Hessian(NM, Int_list, N_int, kappa, anharmonic);
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

  float mass_sum_sqrt=sqrt(mass_sum);
  char *summary, *nameprot;
  if(CONF_CHANGE){
    summary=summary2; nameprot=name2;
  }else{
    summary=summary1; nameprot=nameprot1;
  }
  // Single protein properties
  // Mean_Cart_coll and Coll_thr_cc are computed here
  N_int_long=0;
  for(int i=0; i<N_int; i++){
    if(atoms1[Int_list[i].i2].res-atoms1[Int_list[i].i1].res>2)N_int_long++;
  }

  G_NM_holo=
    Print_summary(summary, nameprot, parameters, Nchain, NM,
		  out_B1, out_B2, Tors_fluct, axe1,naxe1,nmain,nskip, natoms1,
		  nres1, N_diso1, N_int, N_inter, Cont_overlap,
		  N_disc_coll, N_disc_freq, N_disc_out, mass_sum,
		  NM_unfold, N_int_unfold, E_cont_unfold, FOLDING,
		  NM_apo, N_int_inter, E_cont_holo, E_cont_bind, BINDING,
		  Ref_kin.mass_coord, Cart_interface, N_Cart_inter,
		  rigid_inter, N_rigid);
  //if(CONF_CHANGE==0){
  Print_mode_summary(nameout1,"Modes",NM,mass_sum_sqrt,anharmonic,kappa);
    //}
    
  printf("************************************************\n");
  printf("************ Unfolding simulation  *************\n");
  printf("************************************************\n\n");
  
  /*if(UNFOLDING){
    Unfolding(bonds, atoms1, natoms1, axe1, naxe1, nmain,
	      chains1, Nchain, seq1, nres1,
	      Int_list, N_int, Int_KB, NA, K_tors_unfold,
	      file_pdb1, chain1);
	      }*/


  /*********************  Conformation change  **********************/

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
    Confchange=malloc(Ref1.N_Cart*sizeof(float));
    rmsd=rmsd_mclachlan_f(coord_1, coord_2, Ref1.mass_atom, Ref1.N_ref);
    for(i=0; i<N_ref; i++)Confchange[i]=coord_2[i]-coord_1[i];
    
    struct Tors Diff;
    Allocate_tors(&Diff, naxe1, Ref1.N_Cart, N_modes);

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
			       mass_sum_sqrt,rmsd, anharm, kappa);
      }
    }

    //Print_diff_fluct(&Diff, &NM, outlier_tors, nameout2);
    
    printf("Writing %s\n", summary);
    printf("Output conf.change. Time= %.2lf sec.\n\n",
	   (clock()-t0)/nbtops); t0=clock();
    
    if((NMUT!=-2)&&(NMUT!=0)){
      // Analysis of mutation
      float mut_phi[naxe1], mut_CC[N_Cart];
      for(int anharm=0; anharm<n_har; anharm++){
	int mut_ok=Mutation(mut_phi, naxe1, mut_CC, // output
			    KAPPA, Ref_kin, Nmut, AAwt, Posmut, AAmut,
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
      char nameforce[400]; sprintf(nameforce, "%s_force.pdb", nameout2);
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
    Print_modes(N_print, nameout1, "Cart", NM.select, N_Cart,
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
    //double Coll_ave=0; for(ia=0; ia<NM.N; ia++)Coll_ave+=NM.Cart_coll[ia];
    //Coll_ave/=NM.N;
    int N_STEP=40, ip=0;
    for(ia=0; ia<NM.N; ia++){
      if((NM.sigma2[ia]==0)||(NM.select[ia]==0))continue;
      Print_mode_PDB(atoms1, natoms1, axe1, naxe1, bonds, seq1,
		     NM.Tors[ia], NM.omega[ia], N_STEP,
		     Para_simul, nameout1, ia, ip);
      ip++; if(ip==N_MODE_PRINT)break;
    }
    Set_bonds_measure(bonds, natoms1, atoms1);
    printf("Printing modes as PDB. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
    t0=clock();
  }

  /**********************  Predict RMSD of mutations ***********************/
  if(PRED_MUT){
    Predict_mutations(NM, KAPPA, atoms1, natoms1, naxe1, Ref_kin,
		      Int_list, N_int, seq1, nres1, nameout1,
		      Mut_para, PRED_MUT, IWT);
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
			Confchange, file_pdb1, chain1, nres1, SITES, anharm,
			Int_list_inter, N_int_inter);
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
	if(anhar==0){strcat(name_ens,"_harmonic");}
	else{strcat(name_ens,"_anharmonic");}
	if(Simulate_ensemble(Para_simul.N_SIMUL, fact, name_ens,
			     Para_simul, mass_sum_sqrt, atoms1, natoms1,
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
    float coord_ini[N_Cart], coord_end[N_Cart];
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
  Clean_memory(NM, J, Hessian, naxe1, N_ref, N_Cart, N_modes, ANM);
  Empty_Ref(&Ref_kin);

  printf("Total Time= %.2lf sec.\n",(clock()-t_all)/nbtops);
  return(0);
  
}

/******************  Other routines  ***************************/

void Rescale_para(float kappa)
{
  printf("Rescaling force constant by factor %.4g\n", kappa);
  KAPPA*=kappa;
  K_PSI*=kappa; K_PHI*=kappa; K_OMEGA*=kappa; K_CHI*=kappa;
  K_BA*=kappa; K_BL*=kappa;
  K0=KAPPA*pow(C0,EXP_HESSIAN);
}

double Rescale_Hessian(struct Normal_Mode NM,
		       struct interaction *Int_list, int N_int, 
		       float kappa, int anharmonic)
{

  // rescale interactions
  for(int i=0; i<N_int; i++)Int_list[i].sec_der*=kappa;

  double sum_sigma2=0; int i;
  printf("Rescaling omega2 by factor %.3g\n", kappa);

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
		     int N_axes, int N_ref, int N_Cart,
		     int N_modes, int N_res)
{

  // Hessian
  *Hessian=Allocate_mat2_d(N_modes, N_modes);
  printf("Hessian allocated\n");

  // Kinematics:
  Allocate_Jacobian(J, N_axes, N_Cart);
  printf("Jacobian allocated\n");

  // Normal modes
  Allocate_Normal_modes(NM, N_modes, N_axes, N_Cart);
  printf("Normal modes allocated\n");

}


void Clean_memory(struct Normal_Mode NM,
		  struct Jacobian J,
		  double **Hessian,
		  int N_axes, int N_ref, int N_Cart,
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

double Print_summary(char *name_out, char *name1, char *parameters, int Nchain,
		     struct Normal_Mode NM,
		     char *out_B1, char *out_B2, float Tors_fluct,
		     struct axe *axes, int N_axes, int N_main,
		     int Nskip, int N_atoms, int N_res, int N_diso,
		     int N_int, int N_inter, float Cont_overlap,
		     int N_disc_coll, int N_disc_freq, int N_disc_out,
		     float mass_sum,
		     struct Normal_Mode NM_unfold,
		     int N_int_unfold, double E_cont_unfold, int FOLDING,
		     struct Normal_Mode NM_apo, int N_int_inter,
		     double Econt_holo, double Econt_bind, int BINDING,
		     float *mass_coord, int *Cart_interface, int N_Cart_inter,
		     int *rigid_inter, int N_rigid)
{
  FILE *file_out;

  printf("Writing %s\n", name_out);
  file_out=fopen(name_out, "w");

  //float r02=pow(N_res, 2/3);
  //float c=(float)N_int/(float)N_res;
  //float Surf=(cont_inf -c)*r02;
  float Surf=N_res-N_int_long/cont_inf;

  fprintf(file_out, "%s\n", parameters);
  fprintf(file_out, "protein               %s\n", name1);
  fprintf(file_out, "Nchain                %d\n", Nchain);
  fprintf(file_out, "Nres                  %d\n", N_res);
  fprintf(file_out, "Mass                  %.0f\n", mass_sum);
  fprintf(file_out, "Degrees of freedom    %d\n", NM.N);
  fprintf(file_out, "Rigid Deg of freed    %d\n", N_rigid);
  fprintf(file_out, "Side chain            %d\n", N_axes-N_main);
  fprintf(file_out, "Skipped main chain    %d\n", Nskip);
  fprintf(file_out, "Discarded modes coll. %d\n", N_disc_coll);
  fprintf(file_out, "Discarded modes freq. %d\n", N_disc_freq);
  fprintf(file_out, "Discarded modes outl. %d\n", N_disc_out);
  fprintf(file_out, "Disordered gaps       %d\n", N_diso);
  fprintf(file_out, "Native interactions   %d\n", N_int);
  if(Cont_overlap >= 0)
    fprintf(file_out, "Overlap with MIN int. %.3f\n", Cont_overlap);
  if(Nchain>1)fprintf(file_out, "Interchain contacts    %d\n", N_inter);
  fprintf(file_out, "Force constant at %.1fA %.4g\n", C0, KAPPA);
  if(TEMPERATURE>0)
    fprintf(file_out, "Temperature:          %.3g K\n", TEMPERATURE);
  fprintf(file_out, "Estimated surface     %.4g\n", Surf);
  fprintf(file_out, "Contact energy        %.4g\n", Econt_holo);

  int anharmonic=0;
  G_NM_holo=
    Print_thermal(anharmonic, file_out, NM, out_B1, axes, N_main, N_res,Nchain,
		  NM_unfold, N_int_unfold, N_int, E_cont_unfold, FOLDING,
		  NM_apo, N_int_inter, Econt_bind, BINDING, 
		  mass_coord,Cart_interface,N_Cart_inter,rigid_inter,N_rigid);
  fprintf(file_out, "Torsional_RMSD_Harm   %.3f\n", Tors_fluct);
  fprintf(file_out, "#\n");

  if(ANHARMONIC){
    anharmonic=1;
    Print_thermal(anharmonic, file_out, NM, out_B2, axes,N_main,N_res,Nchain,
		  NM_unfold, N_int_unfold, N_int, E_cont_unfold, FOLDING,
		  NM_apo, N_int_inter, Econt_bind, BINDING, 
		  mass_coord,Cart_interface,N_Cart_inter,rigid_inter,N_rigid);
    fprintf(file_out, "#\n");
  }
  fclose(file_out);
  return(G_NM_holo);
}

double Print_thermal(int anharmonic, FILE *file_out,
		     struct Normal_Mode NM, char *out_B,
		     struct axe *axes, int N_main,
		     int N_res, int Nchain,
		     struct Normal_Mode NM_unfold,
		     int N_int_unfold, int N_int, 
		     double E_cont_unfold, int FOLDING,
		     struct Normal_Mode NM_apo, int N_int_inter,
		     double Econt_bind, int BINDING,
		     float *mass_coord, int *Cart_interface, int N_Cart_inter,
		     int *rigid_inter, int N_rigid)
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
  double Free_Ene_holo=0, Free_Ene_holo_rigid=0;

  {  /* Harmonic entropy quantum
	h/= 1.055 e-34 J.sec
	h = 6.626 e-34 J.sec
	k_B T at 293K = 1.381 e-23 *293= 4.05e-21 J
	/h/k_BT = 0.260e-13 s^(-1) = 2.6 ps
	h/k_BT = 1.637e-13 s^(-1) = 0.1636 ps
	omega_q = kBT/h/ = 6.1 ps^(-1)
	(agrees with K. Hinsen)

	Our internal units of time are:
	[t]=[m l^2/E]^(1/2)
	Energy unit: kT= 4.05e-21 J at T=293
	Mass unit: Atomic mass 1.66e-27 kg
	mass/Energy: 0.4099e-6 sqrt: 0.64 e-3
	//sqrt(k_B T/M)= sqrt(4.05e-21 J/1.66e-27 Kg)=1562
	length unit = e-10 m
	time unit: 0.64 e-13 = 6.4 e-12 = 6.4 ps
	Frequency unit: 0.156 ps^(-1)
	omega_q = 39 internal units
	1/(omega_q)^2= 0.00065 in internal units 

	Other constants:
	Avogadro number NA= 6.02e+23
	Gas constant R= 8.314 J/(K*mol) 

     */

    double omega_q=6.1*(Temp_comp/293)/Freq_unit;
    double beta_omega=1./omega_q;
    double entropy=
      Compute_entropy(&Free_Ene_holo, &Free_Ene_holo_rigid,
		      NM, sigma2, beta_omega, NULL, NULL, 0, NULL, 0);
    if(anharmonic==0){
      fprintf(file_out,"kT/h=                %.4g\n", omega_q);
    }
    int m=0, k; for(k=0; k<NM.N; k++){if(NM.omega[k]>omega_q)break; m++;}
    fprintf(file_out, "Quantum modes:        %d / %d\n", NM.N-m, NM.N);
    fprintf(file_out, "Entropy per res       %.4g", entropy/N_res);
    if(BINDING){fprintf(file_out, " (holo)");}
    fprintf(file_out, "\n");
    fprintf(file_out, "Free energy per res   %.4g", Free_Ene_holo/N_res);
    if(BINDING){fprintf(file_out, " (holo)");}
    fprintf(file_out, "\n");

    if(FOLDING){
      float *s2_unfold;
      if(anharmonic){s2_unfold=NM_unfold.sigma2_anhar;}
      else{s2_unfold=NM_unfold.sigma2;}
      double Free_Ene_unfold=0, Free_Ene_unfold_rigid=0;
      double entropy_unfold= 
	Compute_entropy(&Free_Ene_unfold, &Free_Ene_unfold_rigid, NM,
			s2_unfold, beta_omega, NULL, NULL, 0, NULL, 0);
      fprintf(file_out, "Unfolded contacts:   %d\n", N_int_unfold);
      fprintf(file_out, "Contact energy unfold %.4g\n", E_cont_unfold);
      fprintf(file_out, "Entropy per res       %.4g (unfolded)\n",
	      (entropy_unfold)/N_res);
      fprintf(file_out, "Free energy per res   %.4g (unfolded)\n",
	      Free_Ene_unfold/N_res);
      double DG_NM_unfold=(Free_Ene_holo-Free_Ene_unfold)/N_res;
      fprintf(file_out, "Unfold free energy NM  %.4g\n", DG_NM_unfold);
      fprintf(file_out, "Unfold entropy NM      %.4g\n",
	      (entropy_unfold-entropy)/N_res);


      //float r0=pow(N_res, 1/3);
      //float c=(float)N_int/(float)N_res;
      //float Surf=(cont_inf -c)*r0*r0;
      float Surf=N_res-N_int_long/cont_inf;


      //float c_1=(float)N_int_unfold/(float)N_res;
      //float Surf_1=(3.94 -c_1)*r0*r0;

      fprintf(file_out, "SASA diff unfold       %.4g\n",
	      N_res-Surf);
    }
  
    if(BINDING){
      float *s2_apo;
      if(anharmonic){s2_apo=NM_apo.sigma2_anhar;}
      else{s2_apo=NM_apo.sigma2;}
      entropy= //double entropy_holo=
	Compute_entropy(&Free_Ene_holo, &Free_Ene_holo_rigid, NM, sigma2,
			beta_omega, mass_coord, Cart_interface, N_Cart_inter,
			rigid_inter, N_rigid);
      double Free_Ene_apo=0, Free_Ene_apo_rigid=0;
      entropy=
	Compute_entropy(&Free_Ene_apo, &Free_Ene_apo_rigid, NM_apo, s2_apo,
			beta_omega, mass_coord, Cart_interface, N_Cart_inter,
			rigid_inter, N_rigid);
      fprintf(file_out, "Interface contacts:   %d\n", N_int_inter);
      fprintf(file_out, "Entropy per res       %.4g (apo, interface)\n",
	      (entropy)/N_res);
      fprintf(file_out, "Free energy per res   %.4g (apo, intern)\n",
	      Free_Ene_apo/N_res);
      fprintf(file_out, "Free energy rigid     %.4g (apo, /nrigid)\n",
	      Free_Ene_apo_rigid/N_rigid);
      DG_NM_internal=(Free_Ene_holo-Free_Ene_apo)/N_res;
      DG_NM_rigid=(Free_Ene_holo_rigid-Free_Ene_apo_rigid)/N_res;
      //N_rigid_inter;
      Free_Ene_holo/=N_res;

      fprintf(file_out, "Binding free energy NM  %.4g (internal, /Nres)\n",
	      DG_NM_internal);
      fprintf(file_out, "Binding free energy NM  %.4g (rigid body, /nrigid)\n",
	      DG_NM_rigid);

      //fprintf(file_out, "Binding entropy NM      %.4g\n",
      //	      entropy_holo-entropy);
      fprintf(file_out, "Binding free energy cont %.4g\n", Econt_bind);

      int N_int=N_cont_long_1+N_cont_long_2+N_int_inter;
      float Surf=N_res_1+N_res_2-N_int/cont_inf;
      float Surf_1=N_res_1-N_cont_long_1/cont_inf;
      float Surf_2=N_res_2-N_cont_long_2/cont_inf;

      fprintf(file_out, "SASA difference          %.4g\n",
	      Surf_1+Surf_2-Surf);
    }
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
    fprintf(file_out, "Fraction_RB_tra       %.3g\n", tra/all);
    fprintf(file_out, "Fraction_RB_rot       %.3g\n", rot/all);
  }
  return(Free_Ene_holo);
}

double Compute_entropy(double *Free_Ene, double *Free_Ene_rigid,
		       struct Normal_Mode NM, float *sigma2,
		       double beta_omega, float *mass_coord,
		       int *Cart_interface, int N_Cart_inter,
		       int *rigid_inter, int N_rigid)
{
  double entropy=0; *Free_Ene=0;
  for(int i=0; i<NM.N; i++){
    if((NM.select[i]==0)||(sigma2[i]<=0))continue;
    double w_interface=1;
    if(Cart_interface){
      w_interface=0; float *xi=NM.Cart[i];
      for(int k=0; k<N_Cart_inter; k++){
	int j=Cart_interface[k]; w_interface+=mass_coord[j]*xi[j]*xi[j];
      }
    }
    double w_internal=1;
    if(rigid_inter){
      float *xi=NM.MW_Tors[i];
      for(int k=0; k<N_rigid; k++){
	int a=rigid_inter[k]; if(a>=0)w_internal-=xi[a]*xi[a];
      }
    }

    float omega=1./sqrt(sigma2[i]);
    double x=beta_omega*omega, f=exp(-x);
    if(f<1){
      double lf= log(1.-f), G=w_interface*(x/2+lf);
      entropy+=w_internal*w_interface*(-lf+x*f/(1.-f));
      (*Free_Ene)+=w_internal*G;
      (*Free_Ene_rigid)+=(1-w_internal)*G;
    }else{
      printf("WARNING, exp(-h/omega/kT)= %.2g\n", f);
    }
  }
  return(entropy);
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
	    float *RMSD_MIN, int *NMUT,
	    int *PRED_MUT, int *IWT, char *Mut_para,
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
	    int *PRINT_PDB_ANHARM, int *N_PDB_print,
	    int *FOLDING, int *BINDING, int *UNFOLDING,
	    int *n_apo, char *apo_chains)
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
  printf("Reading parameters from command line\n");

  char input[80]="";
  for(i=0; i<argc; i++){
    if(strncmp(argv[i],"-file",3)==0){
      i++; if(i<argc){strcpy(input,argv[i]);}
      break;
    }
  }

  if(input[0]!='\0' && input[0]!='-'){
    printf("Reading parameters from file %s\n", input);
    Read_para(input, file_pdb1, chain1,file_pdb2, chain2,
		       ANM, REF, LABEL, SIDECHAINS, OMEGA, PSI,
		       K_OMEGA, K_PSI, K_PHI, K_CHI, K_BA, K_BL,
		       E_MIN, COLL_THR, MIN_INT, INT_TYPE, C_THR,
		       S_TYPE, S_THR, ONEINT, N_RESRES,
		       N_MODE, outdir, PRINT_CONFCHANGE, PRINT_FORCE,
		       PRINT_PDB, PRINT_MODE_SUMM, FILE_FORCE, ALLOSTERY,
		       KAPPA, FIT_B, RMSD_EXP, RMSD_MIN,
		       NMUT, PRED_MUT, IWT, Mut_para,
		       Para_simul, Para_confchange, SEQID_THR,
		       PRINT_COV_COUPLING, PRINT_DEF_COUPLING,
		       PRINT_DIR_COUPLING, PRINT_COORD_COUPLING,
		       PRINT_SIGMA_DIJ, PROF_TYPE, ALL_PAIRS, SIGMA,
		       STRAIN, SITES, PRINT_CMAT, ANHARMONIC,
		       PRINT_PDB_ANHARM, N_PDB_print,
		       FOLDING, BINDING, UNFOLDING, n_apo, apo_chains);
  }else{
    printf("WARNING, input file not specified or absent\n");
  }

  printf("Reading parameters from command line\n");
  for(i=0; i<argc; i++){
    if(strncmp(argv[i],"-file",3)==0){
      i++;
      // Proteins:
    }else if (strcmp(argv[i],"-p2")==0 || strcmp(argv[i],"-pdb2")==0){
      i++; if(i>=argc)continue;
      strcpy(file_pdb2,argv[i]); //p2=1;
      printf("file_pdb2=%s\n",file_pdb2);
    }else if (strncmp(argv[i],"-c2",3)==0){
      i++; if(i>=argc)continue;
      if((strncmp(argv[i], "all", 3)==0)||(strncmp(argv[i], "ALL",3)==0))
	{*chain2='*';}
      else{strcpy(chain2,argv[i]);}
      printf("chain2=%s\n",chain2);
    }else if (strcmp(argv[i],"-p1")==0 ||
	      strcmp(argv[i],"-pdb1")==0 || strcmp(argv[i],"-pdb")==0){
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
    }else if (strncmp(argv[i],"-pred_mut",9)==0){
      *PRED_MUT=1;
    }else if (strncmp(argv[i],"-mut_para",9)==0){
      i++; strcpy(Mut_para, argv[i]);
    }else if (strncmp(argv[i],"-couplings",8)==0){
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


  //modyves: some changes to not sprintf parameters into parameters, which gave me some warnings.
  char tmp[100];
  sprintf(parameters, "Para: Psi=%d  Omega=%d Sidechain=%d Min_int=%d  REF=%s",
	  *PSI, *OMEGA, *SIDECHAINS, *MIN_INT, REF);
  sprintf(tmp," E_MIN=%.2g COLL_THR=%.0f", *E_MIN, *COLL_THR);
  strcat(parameters,tmp);
  sprintf(tmp," K_O=%.2g K_F=%.2g K_P=%.2g K_C=%.2g K_BA=%.2g K_BL=%.2g",
	  *K_OMEGA, *K_PHI, *K_PSI, *K_CHI, *K_BA, *K_BL);
  strcat(parameters,tmp);
  sprintf(tmp," CONT=%s ", INT_TYPE); strcat(parameters,tmp);
  if(strncmp(INT_TYPE, "SCR", 3)==0){
    sprintf(tmp, " Tolerance=%.2f THR=%.2f", *S_THR, *C_THR);
  }else if(strncmp(INT_TYPE, "SHA", 3)==0){
    sprintf(tmp, " type=%d d=%.2f THR=%.1f", *S_TYPE, *S_THR, *C_THR);
  }else{
    sprintf(tmp, " thr=%.2f", *C_THR);
  }
  strcat(parameters,tmp);
  sprintf(tmp, " EXPO=%.1f ", EXP_HESSIAN); strcat(parameters,tmp);
  if(*FIT_B==0){sprintf(tmp, " KAPPA=%.1f", *KAPPA);}
  else{sprintf(tmp, " KAPPA obtained from fit");}
  strcat(parameters,tmp);
  fflush(stdout);
  return (0);
}

void help(void)
{
  fprintf(stderr,
	  "Program %s\n"
	  "author Ugo Bastolla <ubastolla@cbm.csic.es> "
	  "Centro de Biologia Molecular Severo Ochoa (CSIC-UAM)\n"
	  "Computes normal modes of proteins with the TNM or ANM.\n"
	  "If two structures are provided, it projects the "
	  "experimental conformation change over normal modes\n", PRG);

  fprintf(stderr,
	  "FORMAT of input file (lines with # are comments):\n"
	  "####   INPUT:\n"
	  "PDB1= /data/ortizg/databases/pdb/1usg.pdb ! Reference structure\n"
	  "CH1=  A                                   ! Chain\n"
	  "# Conformation change (optional):\n"
	  "PDB2= /data/ortizg/database/pdb/1usk.pdb  ! Conformation change\n"
	  "CH2=  A                                   ! Chain\n"
	  "#\n"
	  "BINDING=  XY     ! Compute binding free energy between\n"
	  "#                ! chains XY and the other chains (optional)\n"
	  "FOLDING=   0,1   ! Compute unfolding entropy\n"
	  "UNFOLDING= 0,1   ! Perform unfolding simulation\n"
	  "####  Model parameters (reccomended: do not change)\n"
	  "#===============================================================\n"
	  "DOF= TORS     ! Internal/Cartesian coordinates Allowed: TORS CART\n"
	  "REF= ALL      ! Reference atoms allowed: ALL EB CB  CA BB  \n"
	  "PSI=1         ! Allow PSI angles to rotate?\n"
	  "OMEGA=0       ! Allow OMEGA angles to rotate?\n"
	  "SIDECHAIN=0   ! Allow sidechain angles to rotate?\n"
	  "# The model with fewer degrees of freedom is faster and"
	  " more robust w.r.t. parameter values\n"
	  "# More degrees of freedom introduce non collective motions"
	  " that may worsen performances.\n"
	  "MIN_INT=%d    ! Min. number of interactions to accept a dof\n"
	  "E_MIN=%.1g  ! Min. eigenvalue/<evalue> of a normal mode\n"
	  "COLL_THR=%.0f   ! Discard modes that move < COLL_THR atoms\n"
	  "K_BL=%.2f     ! Elastic constant for bond lengths\n"
	  "K_BA=%.2f     ! Elastic constant for bond angles\n"
	  "K_OMEGA=%.2f  ! Elastic constant for angle omega\n"
	  "K_PHI=%.2f    ! Elastic constant for angle phi\n"
	  "K_PSI=%.2f    ! Elastic constant for angle psi\n"
	  "K_CHI=%.1f    ! Elastic constant for angle Chi\n"
	  "CONT_TYPE= MIN ! Interaction model (MIN SCR SHA SCB CB CA HB HNM)\n"
	  "CONT_THR= %.1f  ! Threshold for contacts\n"
	  /*"S_THR= 0.2  ! Parameter for screened or shadow interactions\n"
	    "S_TYPE= 2   ! Type of shadow interactions\n"
	    "# 0=angle 1=axis distance 2=S-ball-axis-distance\n"
	    "ONEINT= 1  ! Only one interaction per residue pair?\n"
	    "N_RESRES= 1  ! Max num interactions per residue pair\n");*/
	  "POW= 1         ! Force constant power-law (1) or exponential (0)\n"
	  "EXP_HESSIAN= 0 ! Force constant k~r^(-e)\n"
	  "FIT_B= 1  ! Set force constant by fitting B factors if present?\n"
	  "KAPPA= %.2f ! force constant when not fitting B factors\n",
	  MIN_INT_MAIN, E_MIN, COLL_THR,
	  K_BL, K_BA, K_OMEGA, K_PHI, K_PSI, K_CHI, THR_ALL, KAPPA);


  fprintf(stderr,
	  "#==============================================================\n"
	  "############### OUTPUT ################\n"
	  "#\n### General output:\n"
	  "LABEL= 0    ! Print Model parameters in file names\n"
	  "PRINT_SUMM= 1      ! Print summary of results\n"
	  "PRINT_CONT_MAT=0   ! Print Contact matrix?"
	  "DEBUG= 0           ! Print debugging information\n"
	  "ANHARMONIC= 0      ! Examine anharmonicity of each mode? (slow)\n"
	  "PRINT_PDB_ANHARM=0 ! Print PDB computing the anharmonic energy?\n"
	  "N_PDB_print= 4     ! Number of PDB printed "
	  "for each norm.mode (if PRINT_PDB_ANHARM)\n");

  fprintf(stderr,
	  "#\n### Normal modes:\n"
	  "NMODES= 5          ! Number of printed normal modes\n"
	  "PRINT_PDB= 0       ! Print modes as PDB files?\n"
	  "AMPLITUDE=1.0      ! Max. amplitude of printed modes"
	  " w.r.t. thermal fluctuations\n"
	  "PDB_STEP=0.2       ! Minimum RMSD between printed PDBs\n"
	  "E_THR=5.0          ! Maximum energy of printed structures\n"
	  "D_REP=2.5          ! Maximum distance for computing repulsion ene\n"
	  "#\n#### Ensembles of alternative structures\n"
	  "SIMUL=20	      ! Numb. simulated structures for each amplitude\n"
	  "AMPL_MIN=1.0       ! Minimum amplitude of simulated structures\n"
	  "AMPL_MAX=8.0       ! Maximum amplitude of simulated structures\n"
	  "AMPL_FACT=2.0      ! Factor multiplying conecutive amplitudes\n"
	  "#\n###  Analysis of conformation change\n"
	  "RMSD_MIN= 0.5      ! Min. RMSD for analyzing conformation change\n"
	  "PRINT_CHANGE= 0    ! Print conformation change as PDB?\n"
	  "RMSD_THR= 0.2      ! Do not analyze modes that move < RMSD_THR\n"
	  "STEP_MAX= 0.5      ! Max. step in simulated confchange\n"
	  "STEP_MIN= 0.001    ! Smallest step in simulated confchange\n"
	  "NSTEPS= 100        ! Number of steps in simulated confchange\n"
	  "ANGLE= 0.1         ! Angular step in confchange (radiants)\n"
	  "PRINT_FORCE= 0     ! Print linear resp. force that produces"
	  " the confchange\n"
	  "# FILE_FORCE=     ! Input file with force, compute deformation\n"
	  "#\n#### Prediction and analysis of mutational effects\n"
	  "PRED_MUT=0         ! Predict RMSD of all possible mutations?\n"
	  "MUT_PARA=Mutation_para.in  ! File with mutation parameters\n"
	  "NMUT= -2           ! Analyze mutation as conf.change\n"
	  "# NMUT>0: analyze if numb.mutations=NMUT NMUT=-1: Analyze always"
	  " NMUT=-2: Do not analyze\n");
  fprintf(stderr,
	  "#\n### Dynamical couplings\n"
	  "ALLOSTERY=0    ! Compute dynamical couplings\n"
	  "ALL_PAIRS=0    ! Print couplings for all pairs? (1=YES)\n"
	  "SIGMA=1.0      ! If ALL_PAIRS=0, print if |coupling|>SIGMA*Std.dev\n"
	  "PRINT_COV_COUPLING=0    ! Print covariance couplings\n"
	  "# Output file: <>_covariance_coupling.dat\n"
	  "PRINT_DIR_COUPLING= 0   ! Print directionality couplings\n"
	  "# The coupling Dir_ij is the Boltzmann average of the scalar\n"
	  "# product of the direction of motion of residues i and j. If it is\n"
	  "# positive the two residues tend to move in similar directions.\n"
	  "# Output files: <>_directionality_coupling.dat and\n"
	  "# <>_directionality_coupling_neg.dat (<0)\n"
	  "PRINT_COORD_COUPLING=0  ! Print coordination couplings\n"
	  "# Output file: <>_coordination_coupling.dat\n"
	  "# The coupling Coord_ij is the Boltzmann average of the squared\n"
	  "# fluctuations of the distance d_ij with respect to the equilibrium value.\n"
	  "# If it is small the two residues maintain an almost fixed distance\n"
	  "# during their dynamics.\n"
	  "PRINT_DEF_COUPLING=0    ! Print deformation couplings\n"
	  "# Output file: <>_deformation_coupling.dat\n"
	  "# The coupling Def_ij is the deformation produced on site j by an\n"
	  "# unitary force applied at j in the direction that maximizes the"
	  " deformation.\n"
	  "PROF_TYPE= C   ! Type of profile of the couplings p_i({C_kl})\n"
	  "# A=average E=effective connectivity P=principal eigenvector\n"
	  "PRINT_SIGMA_DIJ=0   ! Print variance of the distance between CA atoms\n"
	  "# (it is the square of the coordination coupling)\n"
	  "# Output file: <>_interatomic_distance_variance.dat"
	  "SITES=leut_sites.in	! Input file with a binding site\n"
	  "# format: SITE AA (3 letter) CHAIN RESNUM (PDB)\n"
	  "# If not specified, read from the SITE record in the PDB file.\n"
	  "# The program computes the mean coupling of pairs of residues"
	  " in each binding site\n");

  fprintf(stderr,
	  "#\n#\n"
	  "Inputs can be given both by file and by command line\n"
	  "   USAGE:"
	  "   %s -file <input_file>  or\n"
	  "   %s -pdb1 <pdbfile1>   (reference pdb structure)\n"
	  "   OPTIONS:\n"
	  "       -h prints this help\n"
	  "       -c1 <chain_id1>  all:read all chains, A, AB..\n"
	  "       -pdb2 <pdbfile2> for conformation change\n"
	  "       -c2 <chain_id2>  all:read all chains, A, AB..\n"
	  "       -ref Reference atoms Allowed: ALL (default), CA CB BB EB\n"
	  "       -omega  Use also omega angle as degree of freedom\n"
	  "       -anm use ANM d.o.freedom (by default TNM is used)\n"
	  "       -cont_type Interaction model: MIN (default) CA CB ALL HYD HB\n"
	  "       -cont_thr Distance threshold default: %.1f\n"
	  "       -expo <e> k~r^(-e) Default: %.0f\n"
	  "       -print_pdb  Print modes in PDB format\n"
	  "       -modes <number of printed modes\n"
	  "       -pred_mut  Print RMSD and DE produced by all possible mutations\n"
	  "       -mut_para <file with parameters of mutation model>\n"
	  "       -print_confchange Print conformation change as PDB\n"
	  "       -print_force  Print force producing the c.change\n"
	  "       -force     Input PDB file with force as coordinates\n"
	  "       -couplings Compute dynamical couplings\n"
	  "       -simul    <Num. simulated structures>\n"
	  "       -hnm HNM (ref=CA, cont_type=CA, e=6, covalent)\n"
	  "       -debug  print debugging information\n"
	  "\n", PRG, PRG, THR_ALL, EXP_HESSIAN);

  exit(1);
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
	      float *RMSD_MIN, int *NMUT,
	      int *PRED_MUT, int *IWT, char *Mut_para,
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
	      int *ANHARMONIC, int *PRINT_PDB_ANHARM, int *N_PDB_print,
	      int *FOLDING, int *BINDING, int *UNFOLDING,
	      int *n_apo, char *apo_chains)
{
  FILE *file_in=fopen(filename, "r");
  if(file_in==NULL){
    printf("ERROR, TNM input file %s not found\n", filename); 
    return(0);
  }
  char string[1000], dumm[80]; int itmp;
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
    }else if(strncmp(string, "IWT", 3)==0){
      sscanf(string+4, "%d",  &itmp);
      if(itmp<0 || itmp >=20) {
	printf("WARNING, IWT=%d not allowed, using default %d\n",itmp,*IWT);
      }else{
	*IWT=itmp;
      }
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

    }else if(strncmp(string, "FOLDING", 7)==0){
      sscanf(string+8, "%d", FOLDING);
    }else if(strncmp(string, "UNFOLDING",9)==0){
      sscanf(string+10, "%d", UNFOLDING);
    }else if(strncmp(string, "BINDING", 7)==0){
      *n_apo=0; char *s=string+8;
      while(*s!='\n' && *s!='!'){
	if(*s!=' '){apo_chains[*n_apo]=*s; (*n_apo)++;}
	s++;
      }
      if(*n_apo>0 && *n_apo<100){
	*BINDING=1; apo_chains[*n_apo]='\0';
	printf("Computing binding free energy between chains ");
	for(int j=0; j<*n_apo; j++)printf("%c", apo_chains[j]);
	printf(" and the other chains\n"); 
      }else{
	printf("WARNING, %d is not an allowed value of n_apo\n", *n_apo);
	printf("Read string: %s", string);
      }

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

void Get_modes_tors(struct Normal_Mode NM, int N_modes,
		    double **Hessian, struct Jacobian J,
		    struct axe *axe1, int naxe1, struct Reference Ref_kin,
		    atom *atoms1)
{
  float t0=clock();

  d_Diagonalize(N_modes, Hessian, NM.omega2, NM.MW_Tors, -1);
  printf("Hessian diagonalization. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
  t0=clock();

  // UGO: Rescale omega^2 with mass 09/08/2022 It may induce mistakes
  //float omega_norm=Ref_kin.mass_sum/Ref_kin.N_ref;
  //for(int i=0; i<N_modes; i++)NM.omega2[i]*=omega_norm;
  
  // Allocate only relevant modes
  NM.N_relevant=NM.N;
  for(int ia=0; ia<NM.N; ia++){
    Transform_tors_modes(NM.Tors[ia], NM.MW_Tors[ia],
			 J.T_sqrt_tr, J.T_sqrt_inv_tr,
			 N_modes, naxe1);
    Convert_torsion2cart(NM.Cart[ia], atoms1, NM.Tors[ia],
			 axe1, naxe1, Ref_kin, ia);
  }
  printf("Torsional and Cartesian modes. Time= %.2lf sec.\n",
	 (clock()-t0)/nbtops);
}


void Get_modes_Cart(struct Normal_Mode NM, int N_modes,
		    double **Hessian, struct Jacobian J,
		    struct axe *axe1, int naxe1, struct Reference Ref_kin)
{
  float t0=clock();
  d_Diagonalize(N_modes, Hessian, NM.omega2, NM.Cart, -1);
  printf("Hessian diagonalization. Time= %.2lf sec.\n", (clock()-t0)/nbtops);
  t0=clock();
    
  //  Convert normal modes
  struct Tors D_NM;
  D_NM.N_axes=naxe1;
  D_NM.N_Cart=NM.N_Cart;
  D_NM.coeff=malloc(1*sizeof(float));
  for(int ia=0; ia<NM.N; ia++){
    // Normalize_vector_weighted(NM.Cart[ia], Ref_kin.mass_coord, N_modes);
    D_NM.Cart=NM.Cart[ia];
    D_NM.Tors=NM.Tors[ia];
    D_NM.MW_Tors=NM.MW_Tors[ia];
    Convert_cart2torsion(&D_NM, Ref_kin, &J);
    NM.Tors_frac[ia]=Tors_fraction(&D_NM, Ref_kin.mass_coord);
  }
  printf("Torsional and Cartesian modes. Time= %.2lf sec.\n",
	 (clock()-t0)/nbtops);
}

int Get_chain_num(int *apo_chain_num, char *apo_chains, int n_apo,
		   struct chain *chains, int Nchain){
  int success=1;
  for(int i=0; i<n_apo; i++){
    int j;
    for(j=0; j<Nchain; j++)if(chains[j].label==apo_chains[i])break;
    if(j==Nchain){
      printf("ERROR, chain %c not found in PDB\n",apo_chains[i]);
      success=0;
    }
    apo_chain_num[i]=j;
  }
  return(success);
}

int Belong_to_apo(int chain, int *apo_chain_num, int n_apo){
  for(int j=0; j<n_apo; j++)if(apo_chain_num[j]==chain)return(1);
  return(0);
}

struct interaction *Apo_interactions(int *N_int_apo,
				     struct interaction **Int_list_inter,
				     int *N_int_inter,
				     struct interaction *Int_list, int N_int,
				     atom *atoms, int *apo_chain_num,
				     int n_apo, int Nchain)
{
  struct interaction *Int_list_apo=
    malloc(N_int*sizeof(struct interaction)), *tmp=Int_list;
  *N_int_apo=0; *N_int_inter=0;
  for(int k=0; k<N_int; k++){
    if(atoms[tmp->i1].chain!=atoms[tmp->i2].chain){
      int apo1=Belong_to_apo(atoms[tmp->i1].chain, apo_chain_num, n_apo);
      int apo2=Belong_to_apo(atoms[tmp->i2].chain, apo_chain_num, n_apo);
      if((apo1 && apo2==0) || (apo1==0 && apo2)){
	(*N_int_inter)++; goto next_int;
      }
    }
    Int_list_apo[*N_int_apo]=*tmp; (*N_int_apo)++; 
  next_int:
    tmp++;
  }
  printf("%d interactions belong to the interface "
	 "between %d chains and %d chains\n",
	 N_int-*N_int_apo, n_apo, Nchain-n_apo);

  *Int_list_inter=malloc(*N_int_inter*sizeof(struct interaction));
  int n=0; tmp=Int_list; struct interaction *store=*Int_list_inter;
  for(int k=0; k<N_int; k++){
    if(atoms[tmp->i1].chain!=atoms[tmp->i2].chain){
      int apo1=Belong_to_apo(atoms[tmp->i1].chain, apo_chain_num, n_apo);
      int apo2=Belong_to_apo(atoms[tmp->i2].chain, apo_chain_num, n_apo);
      if((apo1 && apo2==0) || (apo1==0 && apo2)){
	*store=*tmp; store++; n++; 
      }
    }
    tmp++;
  }
  if(n!=*N_int_inter){printf("ERROR of Interface interactions\n"); exit(8);}

  return(Int_list_apo);
}

int Get_interface(int *Cart_interface,
		  struct interaction *Int_list, int N_int,
		  atom *atoms, int *apo_chain_num,
		  int n_apo, struct chain *chains, int Nchain,
		  int *atom_num, int N_ref, int natoms, int nres)
{
  N_res_1=0; N_res_2=0; N_cont_long_1=0; N_cont_long_2=0;

  int atom_ref[natoms], i;
  for(i=0; i<natoms; i++)atom_ref[i]=-1;
  for(i=0; i<N_ref; i++)atom_ref[atom_num[i]]=i;
  int res_inter[nres]; for(i=0; i<nres; i++)res_inter[i]=0;
  for(i=0; i<Nchain; i++){
    if(Belong_to_apo(i, apo_chain_num, n_apo)){N_res_1+=chains[i].nres;}
    else{N_res_2+=chains[i].nres;}
  }

  struct interaction *tmp=Int_list;
  for(int k=0; k<N_int; k++){
    if(atoms[tmp->i1].chain!=atoms[tmp->i2].chain){
      int apo1=Belong_to_apo(atoms[tmp->i1].chain, apo_chain_num, n_apo);
      int apo2=Belong_to_apo(atoms[tmp->i2].chain, apo_chain_num, n_apo);
      if(apo1 && apo2){
	if(atoms[tmp->i2].res-atoms[tmp->i1].res>2)N_cont_long_1++;
      }else if(apo1==0 && apo2==0){
	if(atoms[tmp->i2].res-atoms[tmp->i1].res>2)N_cont_long_2++;
      }else if((apo1 && apo2==0)||(apo1==0 && apo2)){
	res_inter[atoms[tmp->i1].res]=1;
	res_inter[atoms[tmp->i2].res]=1;
      }
    }
    tmp++;
  }
  int nres_inter=0, natom_inter=0, N_Cart_inter=0;
  for(i=0; i<nres; i++)if(res_inter[i])nres_inter++;
  for(i=0; i<natoms; i++){
    if(res_inter[atoms[i].res] && atom_ref[i]>=0){
      natom_inter++; int k=3*atom_ref[i];
      for(int j=0; j<3; j++){
	Cart_interface[N_Cart_inter]=k; N_Cart_inter++; k++;
      }
    }
  }

  printf("Interface between %d chains and %d chains:\n", n_apo, Nchain-n_apo);
  printf("%d residues out of %d "
	 "%d atoms out of %d "
	 "%d Cartesian coordinates out of %d\n",
	 nres_inter, nres, natom_inter, natoms, N_Cart_inter, 3*N_ref);
  return(N_Cart_inter);
}


int Compute_sigma2(int *N_disc_freq, int *N_disc_out,
		   struct Normal_Mode NM,
		   float *mass_sum_sqrt, int naxe1, int N_Cart,
		   float COLL_THR)
{
  int N_disc_coll=0; *N_disc_freq=0;

  for(int ia=0; ia<NM.N; ia++){

    NM.Max_dev[ia]=
      Compute_Max_dev(NM.Cart[ia], NM.omega2[ia], N_Cart, mass_sum_sqrt);
    NM.Tors_coll[ia]=Collectivity_norm2(NM.Tors[ia], naxe1)/naxe1;
    NM.MW_Tors_coll[ia]=Collectivity_norm2(NM.MW_Tors[ia],naxe1)/naxe1;
    NM.Cart_coll[ia]=Collectivity_norm2(NM.Cart[ia], N_Cart); //_Renyi
    if((NM.omega2[ia]<=E_THR)||(NM.Cart_coll[ia] < COLL_THR)){
      printf("WARNING, discarding mode %d ", ia);
      printf(" omega^2= %.3g Coll= %.0f RMSD= %.2g\n",
	     NM.omega2[ia], NM.Cart_coll[ia], NM.Max_dev[ia]);
      if(NM.Cart_coll[ia] < COLL_THR){N_disc_coll++;}
      else{(*N_disc_freq)++;}
      NM.select[ia]=0; NM.sigma2[ia]=0;
    }else{
      NM.select[ia]=1;
      NM.sigma2[ia]=1./NM.omega2[ia];
    }
    NM.Cart_coll[ia]/=N_Cart;
  }
  printf("%d modes selected\n", NM.N-(N_disc_coll+*N_disc_freq));
  printf("%d modes discarded for low collectivity ", N_disc_coll);
  printf("and %d for low frequency\n", *N_disc_freq);
  int amin=0; while(NM.omega2[amin]<0)amin++;
  printf("Smallest positive eigenvalue: %.2g a=%d ",NM.omega2[amin],amin);
  printf("Threshold for selection: %.2g*<omega^2>\n", E_MIN);

  // Ugo_mod 05/21
  /* float thr_outlier=2.5, thr_mode=0.5;
     printf("Remove modes that move > %.2f outlier dof (rmsf > %.1f)\n",
     thr_mode, thr_outlier);
     *N_disc_out=Filter_modes_outliers(NM, Ref_kin, thr_outlier, thr_mode);
     printf("%d modes removed because they mostly move outliers\n",
     *N_disc_out);*/

  return(N_disc_coll);
}

struct interaction *Unfold_interactions(int *N_int_unfold,
					struct interaction *Int_list,
					int N_int, atom *atoms)
{
  struct interaction *Int_list_unfold=
    malloc(N_int*sizeof(struct interaction)), *tmp=Int_list;
  int DIJ_MIN=2, n=0;
  for(int k=0; k<N_int; k++){
    if(atoms[tmp->i1].chain==atoms[tmp->i2].chain &&
       abs(atoms[tmp->i1].res-atoms[tmp->i2].res)<=DIJ_MIN){
      Int_list_unfold[n]=*tmp; n++; 
    }
    tmp++;
  }
  *N_int_unfold=n;
  return(Int_list_unfold);
}

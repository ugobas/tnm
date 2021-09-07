#include <stdio.h>
extern int DEBUG,ALL_AXES;
extern char SEL[10];
extern int READ_RESTRAINT;
extern int ini_print;
extern char chain;
// float dist_CA[L_MAX]; // Squared distance from the center of mass
extern char file_aniso[200];
extern int INT_MAX;

/******************** Shared variables ****************************/
#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWHX"
#define AA3 "ALA GLU GLN ASP ASN LEU GLY LYS SER VAL ARG THR PRO ILE MET PHE TYR CYS TRP HIS "

/****************************************************************
                      SHARED ROUTINES
*****************************************************************/

/**************************  Input  ******************************/
// files: read.c, arguments_nma.c

int Read_command_line_nma(int argc, char **argv, char *file_pdb, char *output,
			  char *chain, char *file_para, char *prot_name,
			  int *N_MODES, char *FILE_RESTR);
void GetPdbId(char *pdb_file_in, char *pdbid);

/*********************************  INM ************************************/
void Normal_modes_INM(int N,           //degrees of freedom
		      struct interaction *Int_list, int N_int,
		      atom *atoms, int N_atoms,
		      double **Hessian, float **eigen_vector,
		      float *eigen_value, float *eigen_B,
		      float *B_CA, float *Cart_collectivity);

/*********************************  ANM ************************************/
void Compute_Hessian_ANM(double **Hessian, int N_ref, int *atom_ref,
			 struct interaction *Int_list, int N_int,
			 atom *atoms, int natoms);
void Normal_modes_ANM(// OUTPUT:
		      float *eigen_value, float **eigen_vector,
		      double **Hessian,
		      // float *eigen_B, float *B_CA, float *Cart_collectivity
		      // INPUT:
		      int N,           //degrees of freedom
		      struct interaction *Int_list, int N_int,
		      atom *atoms, int N_atoms);
int Compute_Bfact_ANM(float *B_CA, float *eigen_B, float *eigen_value,
		      float **eigen_vector, int N, int N_atoms);

/**************************** B factors ************************************/

void Predict_fluctuations(struct Normal_Mode *NM, float *sigma2);
void Print_tors_fluct2(struct Normal_Mode NM, struct axe *axe, int naxe,
		       struct residue *res, char *nameout1);
float Compute_correlation(float **Cart_mode, float *sigma2,
			  int N_modes, int i1, int i2);
void Compute_anisou(float ***aniso_pred, int N_modes, int N_ref,
		    float *sigma2, float **Cart_mode);


/*********************************  TNM ************************************/
void Normal_modes_TNM(// OUTPUT:
		      float *eigen_value, float **eigen_vector,
		      float **Masswtd_evector, float **Cart_mode,
		      double **Hessian,
		      // INPUT:
		      int N,           //degrees of freedom
		      struct interaction *Int_list, int N_int,
		      atom *atoms, int N_atoms,
		      struct axe *axes, int N_axes,
		      struct chain *chains, int Nchain,
		      float *mass_coord, int N_ref);
void Compute_Hessian_TNM(double **Hessian, double **T_sqrt, float **T_sqrt_inv,
			 struct interaction *Int_list, int N_int,
			 atom *atoms, int N_atoms, struct axe *axes,
			 int nmain, int N_axes, int N_modes,
			 struct chain *chains, int Nchain, int kinetic,
			 float K_OMEGA, float K_PSI, float K_PHI,
			 float K_CHI, float K_BA, float K_BL);
float **Dr_daxe_int(struct axe *axes, int N, atom *atoms, int N_atoms,
		    int *atom_kyn, float *mass_kyn, int N_atom_kyn);
void Transform_tors_modes(float *Tors_mode, float *Masswtd_Tors_mode,
			  float **T_sqrt_tr, float **T_sqrt_inv_tr,
			  int N_modes, int N_axes);
void Transform_tors_modes_old(float *Tors_mode, float *Masswtd_Tors_mode,
			      double **T_sqrt, double **T_sqrt_inv,
			      int N_modes, int N_axes);
void Convert_torsion2cart(float *Cart_mode, atom *atoms, float *u_a,
			  struct axe *axes, int N_axes,
			  struct Reference Ref, int ia);
void Convert_torsion2cart_old(float *Cart_mode, atom *atoms, float *u_a,
			  struct axe *axes, int N_axes, int *atom_ref,
			      int N_ref, struct Jacobian *J);
int Convert_cart2torsion(struct Tors *Diff, struct Reference Ref,
			 struct Jacobian *J);
int Convert_cart2torsion_fit(struct Tors *D, struct Reference Ref,
			     struct Jacobian *J, char *nameout,
			     char type, float *Lambda);
int Convert_cart2torsion1(float *Tors_dev, float *Masswtd_Tors_dev,
			  float *Cart_dev, struct Reference Ref,
			  struct Jacobian *J);
void Compute_MW(float *MW_Tors, float *Tors, struct Jacobian *J);
float Tors_fraction(struct Tors *Diff, float *mass);
float Convert_cart2torsion_old(float *Tors_dev, float *Masswtd_Tors_dev,
			       float *Cart_dev, int N_axes, int N_cart,
			       double **T_sqrt, double **T_sqrt_inv,
			       double **Jacobian_ar, float *mass);
void Rescale_eigenvector(float *eigen_vector, float *Masswtd_evector,
			 int N, double **T_sqrt);
float Compute_Max_dev(float *Cart_mode, float omega2, int N_cart,
		      float *inv_sq_mass);
int Compute_Bfact_TNM(float *B_CA, float *eigen_B, float *eigen_value,
		      float **Cart_mode, int N, int N_CA, float mass);

/*// Setting axes
int Set_axes_atoms(struct residue *seq, atom *atoms, int Nres,
		   struct axe *axes, int *N_atoms, int *Nchain);
int Set_chains(struct chain *chains,
	       struct axe *axes, int N_axes,
	       atom *atoms, int N_atoms,
	       struct residue *res, int Nchain);
void Set_atom_axe(atom *atoms, struct chain *chains, int Nchain,
		  struct axe *axes);
*/
// Transforming axes
void Principal_axis_frame(atom *atoms, int N_atoms, int *atom_ref,
			  int N_atom_ref, struct axe *axes, int N_axes,
			  int ANISOU);
// Inertial set
void Eckart_main(struct axe *axes, int N_axes, atom *atoms, int N_atoms,
		 struct Reference Ref);
		 //int *atom_ref, double *mass_ref, int N_atom_ref);
/*int Set_reference_ali(// Output:
		      int *atom_ref1, float *mass_atom, int *atom_ref2,
		      // Input:
		      int ini_ref, char *SEL, int *alignres,
		      atom *atoms1, int ini1, int natoms1, int ini_res1,
		      atom *atoms2, int ini2, int natoms2, int ini_res2);
*/
int Set_reference(// Output:
		  struct Reference *Ref,
		  // Input:
		  int ini_ref, char *SEL, atom *atoms, int ini1, int N_atoms);
void Empty_Ref(struct Reference *Ref);
int Align_references(struct ali_atoms *ali_atoms,
		     struct Reference Ref1, atom *atoms1);

int Internal_Jacobian(// Output:
		      struct Jacobian J,
		      // Input:
		      struct Reference Ref,
		      atom *atoms, int N_atoms,
		      struct axe *axe, int N_axes);
//int Internal_Jacobian(// Output:
//		      double **Jacobian_ar, double **Rot_ra, double **Shift_ra,
		      // Input:
//		      int *atom_ref, double *mass_atom, int N_atom_ref,
//		      atom *atoms, int N_atoms, struct axe *axes, int N_axe);
int Kinetic_energy_all(struct Jacobian J, struct Reference Ref, int N_axes);
//int Kinetic_energy_all(double **T_mat, double **Jacobian_ar,
//		       float *mass_coord, int N_axes, int N_ref);
int Kinetic_energy(double ***T_mat, float *mass_coord,
		   struct axe *axe, int ini_axe, int N_axes,
		   atom *atoms, int *atom_num, int N_ref);
int Kinetic_sqrt(struct Jacobian *J, int *N_kin,
		 int N, int CONTROL, float E_MIN);
//int Kinetic_sqrt(double **T_sqrt, double ***T_sqrt_inv, int *N_kin,
//		 int N, int CONTROL, float E_MIN);

void Allocate_Jacobian(struct Jacobian *J, int N_axes, int N_cart);
void Empty_Jacobian(struct Jacobian J);


float Distance_square(double *r1, double *r2);


/************************ Auxiliary computation ****************************/
// atoms
int Find_atom(atom *atoms, int *i1, int res, int N_atoms, char *atom_type);
int Select_atoms(int *atom_num, atom *atoms, int N_atoms, char *SEL);
void Match_references(int *kin_cc, struct Reference Ref1, struct Reference Ref);
float Mass(atom *atom);
float Mass_atom(char *atom_name);

void gaussj0(float **a, int n,float **b, int m);
// Tests
void Compare_modes(float **Cart_mode_ANM, float **Cart_mode_TNM, int N_cart,
		   float **Tors_mode_ANM, float **Tors_mode_TNM, int N_tors,
		   float *mass_coord, char *prot_name, char *inter);
int Test_Eckart(float *Cart_dev, atom *atoms, int N_atoms,
		struct Reference Ref, int mode, struct axe *axe);
//void Test_Eckart(float *Cart_dev, atom *atoms, int N_atoms,
//		 int *atom_ref, double *mass_ref, int N_ref, int mode);

/**************************** B factors *****************************/
int Set_B_exp(atom *atom, int N_atoms, int N_res, float *B_exp, char *REF);
void Set_anisou(atom *atom, int N_atoms, int nca,
		float ***anisou, float ***aniso_pred, char *REF);
float Rigid_fraction(float *f_tras, int n_RB, struct Normal_Mode NM);

/********************************  Output **********************************/
float Output_modes(int N,           //degrees of freedom
		   char *model,     // INM, ANM or TNM
		   char *REF_ATM,   // Reference atoms
		   atom *atoms, int N_atoms,
		   struct residue *seq, int N_res,
		   struct axe *axes, int N_axes,
		   float *B_CA, float *B_CA_exp,
		   int *weight, float *dist_CA,
		   float ***aniso_pred, float ***aniso_exp,
		   int N_MODES, float cutoff,
		   char *inter, char *prot_name,
		   FILE *file_B, FILE *file_sum,
		   float **eigen_vector, float *eigen_value, float *eigen_B,
		   float *Cart_collectivity, float *Tors_collectivity,
		   float **Cart_mode, int *atom_ref, int N_ref);
void Print_PDB_3(float *Cart_mode, float eigen_value, float eigen_B,
		 atom *atoms, int N_atom_ref, int *atom_ref,
		 struct residue *seq, char *file_name, int ia);

void Print_B(float *B, int nca, FILE *file_out, char *model, float cc);
void Name3(char *aaname3, int i_aa);
void Print_B_fact(float *B_TNM, float *B_pred_all, float *B_exp,
		  int N, atom *atoms,  int *atom_ref, struct residue *seq,
		  char *name, char *what,float cc,float slope,
		  float *dof, char ridge);
/* Other
extern void f_Diagonalize(int N, float **MATRIX, float *eigen_values,
			  float **eigen_vector, int SIGN, float E_MIN);
extern void d_Diagonalize(int N, double **MATRIX, float *eigen_values,
			  float **eigen_vector, int SIGN, float E_MIN);
extern void dd_Diagonalize(int N, double **MATRIX, double *eigen_values,
			   double **eigen_vector, int SIGN, float E_MIN);
*/

float Force_constant(float r);
void Compute_force(struct Tors *Force, float *cc, struct Normal_Mode NM);

void Torsional_force(struct Tors Force, struct Tors *Diff,
		     double **Hessian, int N_axes,
		     int N_ref, int *atom_num,
		     float **Jacobian_ar, int N_int,
		     struct interaction *Int_list,
		     atom *atoms, int natoms);

void Cartesian_force(float *Cart_force, float *atom_diff,
		     double **Hessian,int N_cart);

//double Mean_Cart_coll; // Mean collectivity of normal modes
//double Coll_thr_cc;    // Minimal collectivity for conformation change
char REF_CC[5]; // Reference atoms for computing conformation change
void Tors_outliers(int *outlier_tors, int naxe, int *outlier, int Na,
		   struct axe *axe, struct chain *chains, int Nchain);

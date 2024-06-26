#ifndef ___TNM_H___12345
#define ___TNM_H___12345

extern int N_PDB_print;
extern int PRINT_PDB_ANHARM;
extern float KAPPA_DEF;
extern float TEMP_DEF;

/*struct Molecule{
  int N_atoms, N_axes, N_res;
  struct residue *seq;
  struct axe *axe;
  char chain;
  char pdbid[10];
  char nameprot[100];
};

struct Dof{
  int N_DOF;
  int N_main;
  struct axe *main;
  int N_rigid;
  struct axe *rigid;
  int N_side;
  struct axe *side;
  };*/

struct Reference{
  int N_ref, N_Cart;
  float mass_sum;
  float *mass_atom;
  float *mass_coord;
  float *mass_sqrt;
  //float *inv_sq_mass;
  int *atom_num;
};

struct Normal_Mode{
  int N, N_kin;
  int N_axes, N_Cart;
  int N_relevant;
  float *omega;
  float *omega2;
  float *sigma2;
  float *confchange2;
  float **Cart;
  float **Tors;
  float **MW_Tors;
  float *tors_fluct;  // Predicted fluctuations of torsional dofs
  float *cart_fluct;  // Predicted fluctuations of cartesian dofs
  float *Tors_frac;
  float *Cart_coll;
  float *Tors_coll;
  float *MW_Tors_coll;
  float *Max_dev;
  // Anharmonicity
  float *sigma2_anhar;
  float *d_KL;
  float *Anharmonicity;
  float *Anharm_struct;
  float *Max_factor;
  float *Max_RMSD;
  int *sort;
  int *select;
  int ANM;
};

struct Jacobian{
  int N_axes, N_Cart, N_kin;
  float **Jacobian_ar;
  float **Jtilde_ar;
  double **T;
  double **T_sqrt;
  float **T_sqrt_inv;
  float **T_sqrt_inv_tr;
  float **T_sqrt_tr;
  //float **Rot_ra;
  //float **Shift_ra;
};

struct Tors{
  int N_axes, N_Cart;
  float *Cart;
  float *MW_Tors;
  float *Tors;
  float *coeff;
  float Tors_frac;
  float RMSD_NoTors;
  float RMSD;
  float RMSD_W;
  float RMSD_Tors;
  float M;
};

struct Para_simul{
  int N_SIMUL;
  float STEP_MAX;
  float STEP_MIN;
  float ANGLE;
  int NSTEPS;
  int RESET;
  int SELECT_ENE;
  float PDB_STEP;
  float AMPLITUDE;
  float AMPL_MIN;
  float AMPL_MAX;
  float AMPL_FACT;
  float D_REP;
  float MAX_ANGLE;
  float E_THR;
};

struct Para_confchange{
  float RMSD_THR;
  float COLL_THR;
};


void Allocate_tors(struct Tors *X, int N_axes, int N_Cart, int N_modes);
void Empty_tors(struct Tors X);


void Print_force_confchange(struct Normal_Mode NM, struct Tors Diff,
			    atom *atoms, struct axe *axes, int N_axes,
			    struct residue *seq, int N_res,
			    char *nameout2);

int Find_atom(atom *atoms, int *i1, int res, int N_atoms, char *atom_type);
double Energy_anharmonic(float *r, atom *atoms, int natoms,
			 struct residue *seq, int nres,
			 struct interaction *Int_list, int N_int,
			 struct interaction **Int_KB, int NA,
			 struct axe *axe, int naxes, double *delta_phi);
float Energy_clashes(float *atom_coord, int N);
float Anharmonicity(struct Normal_Mode NM, int ia,
		    atom *atoms, int natoms,
		    struct axe *axe, int naxes,
		    struct bond *bonds,
		    struct residue *seq, int nres,
		    struct interaction *Int_list, int N_int,
		    struct interaction **Int_KB, int NA,
		    char *nameout);

void Anharmonicity_analysis(float *Anharmonicity, int *direction,
			    float *Cart, float omega,
			    atom *atoms, int natoms, struct residue *seq,
			    int nres, int *atom_num, int N_ref);

void Empty_Normal_modes(struct Normal_Mode NM);

#endif

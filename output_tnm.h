
int Check_make_dir(char *outdir);
void Print_modes(int N_print, char *nameout,
		 char *label, int *select,
		 int N, float **Mode, float *collectivity,
		 float *sigma2, double sum_sigma2,
		 struct axe *axe, int naxe,
		 atom *atoms1, int natoms,
		 struct residue *seq, int nres,
		 int *atom_num, int N_ref);
int Print_PDB_mode(char *nameout, int ia, float *Cart_mode,
		   float Amplitude, float Cart_collectivity,
		   float Tors_collectivity, float MW_Tors_collectivity,
		   float eigen_value, float eigen_B,
		   atom *atoms, int natoms, struct residue *seq, int nres,
		   int *atom_num, int N_ref);
void Print_mode_summary(char *nameout, char *label, struct Normal_Mode NM,
			float M_sqrt, int anharmonic, float xkappa);
int Print_Bfactors(int N,           //degrees of freedom
		   char *model,     // INM, ANM or TNM
		   char *REF_ATM,   // Reference atoms
		   atom *atoms, int N_atoms,
		   struct residue *seq, int N_res,
		   struct axe *axes, int N_axes,
		   float *B_CA, float *B_CA_exp,
		   int N_MODES, float cutoff,
		   char *inter, char *prot_name,
		   FILE *file_B, FILE *file_sum,
		   float **eigen_vector, float *eigen_value, float *eigen_B,
		   float *Cart_collectivity, float *Tors_collectivity,
		   float **Cart_mode, int *atom_ref, int N_ref);
void Print_tors_fluct(struct axe *axe, int N_axes, float *fluct,
		      atom *atoms, struct residue *seq,
		      char *name, char *type);
void Print_cart_fluct(int *atom_num, int N_atom, float *fluct,
		      atom *atoms, struct residue *seq,
		      char *name, char *type);
void Print_structures(char *pdbout,
		      double *atom_str1, atom *atoms1,
		      struct residue *seq1, char chain1,
		      double *atom_str2, atom *atoms2,
		      struct residue *seq2, char chain2,
		      int *atom_num, int N_ref, int N_cart,
		      float **Cart_mode, float *coeff,
		      int *sort, int N_MODE_PRINT, int N_modes);
int Write_coord(char *name_out, double *atoms_ref,
		atom *atoms, int *atom_num, int N_atom_ref,
		struct residue *seq, char chain, int i_model);
int Print_change(float *Fluct_pred, float *Fluct_obs, int N,
		char *nameout, char *what);
void Print_atom(atom *atom, int *num, float B, struct residue *seq, 
		FILE *file_out);
void Print_diff_fluct(struct Tors *Diff, struct Normal_Mode *NM,
		      int *outlier_tors, char *nameout2);

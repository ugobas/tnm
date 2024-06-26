void Write_ref_coord(float *coord_ref, int N_ref, float *coord_all,
		     int *atom_num);
void Write_ref_coord_atom(float *coord_ref, int N_ref, atom *atoms,
                          int *atom_num);
float Optimize_amplitude(struct bond *bonds_min, float *f_min,
			 float *Diff_Tors, int naxes,
			 struct bond *bonds, int natoms1,
			 struct Reference Ref1, float *coord_ref_2);
void Set_masses(float *mass, atom *atoms, int natoms);
float Generate_coords(float *coord_new, float *RMSD_target,
		      struct Tors *Diff, int naxes, struct bond *bonds,
		      atom *atoms1, int natoms1, struct Reference Ref,
		      float *coord_old, float *coord_target);
float RMSD_LIN(float *RMSD_target, float *RMSD_full,
	       struct Jacobian *J, float *Delta_phi,
	       struct Reference Ref,
	       float *coord_old, float *coord_target,
	       float *coord_full);
int Simulate_ensemble(int N_struct, float factor, char *name,
		      struct Para_simul Para_simul, float Mass_sqrt,
		      atom *atoms, int natoms,
		      struct axe *axe, int naxes,
		      struct bond *bonds,
		      struct residue *seq, int nres,
		      int N_diso, struct Normal_Mode NM,
		      struct Jacobian J, struct Reference Ref,
		      struct interaction *Int_list, int N_int,
		      struct interaction **Int_KB,
		      char *INT_TYPE, float s0, char *nameout,
		      int anhar);
void Simulate_ensemble_Cart(int N_struct, float factor, char *name,
			    struct Para_simul Para_simul,
			    atom *atoms, int natoms,
			    struct axe *axe, int naxes,
			    struct bond *bonds,
			    struct residue *seq, int nres,
			    int N_diso, struct Normal_Mode NM,
			    struct Jacobian J, struct Reference Ref,
			    struct interaction *Int_list, int N_int,
			    char *INT_TYPE, float s0, char *nameout,
			    int anhar);
void Simulate_confchange(int N_frames, float *tors_dir, float *d_Cart,
			 int Print_simul_PDB, float RMSD_NM,
			 float *coord_1, float *coord_2, float rmsd,
			 struct Para_simul Para_simul,
			 atom *atoms, int natoms,
			 struct axe *axe, int naxes,
			 struct bond *bonds,
			 struct residue *seq, int nres,
			 int N_diso, struct Normal_Mode NM,
			 struct Jacobian J, struct Reference Ref,
			 struct interaction *Int_list, int N_int,
			 char *INT_TYPE, float s0, char *nameout,
			 double tors2, int anhar);
float Make_step(atom *atoms_sim, float *d_phi, struct bond *bonds,
		atom *atoms, int natoms, int naxes,
		struct Reference Ref, double *coord_ref);
void Torsional_confchange_RRR(double *diff_phi, char *diff_type,
			      struct bond *BONDS,
			      struct Jacobian *J, struct Reference Ref1,
			      struct ali_atoms ali_a, double t2,
			      atom *atoms1, int natoms, char *nameout,
			      struct axe *axes, int naxes,
			      struct residue *seq, int nres,
			      struct chain *chains, int Nchain, 
			      atom *atoms2, int natoms2,
			      struct Para_simul Para_simul,
			      struct Normal_Mode *NM, int NMODES);
void Torsional_confchange(double *diff_phi, char *diff_type,
			  struct bond *bonds_ini,
			  struct Jacobian *J, struct Reference Ref1,
			  struct ali_atoms ali_a,
			  double Tors_fluct,
			  atom *atoms1, int natoms1, char *nameout,
			  struct axe *axes, int naxes,
			  struct residue *seq, int nres,
			  struct chain *chains, int Nchain, 
			  atom *atoms2, int natoms2,
			  struct Para_simul Para_simul,
			  struct Normal_Mode *NM,
			  char *mode_name, int NMODES, float rmsd_thr);
int  Change_internal(double *diff_phi, int naxe,
		     struct bond *bonds, struct bond *bonds2,
		     int natoms, int type, int same);
int Put_coord(float *coord, struct bond *bonds, int natoms);
int Periodic_angles(double *dphi, struct axe *axe, int n);

int Print_mode_PDB(atom *atoms, int natoms,
		   struct axe *axe, int naxes,
		   struct bond *bonds, struct residue *seq,
		   float *Tors, float omega,
		   int N_STEP, struct Para_simul Para_simul,
		   char *nameout, int ia, int ip);
void Step(atom *atoms_sim, struct bond *bonds,
	  float *d_phi, atom *atoms, int natoms, int naxes);
void Write_all_coord(double *coord, atom *atoms, int natoms);
void Copy_bonds(struct bond *bonds1, struct bond *bonds, int natoms);
int Move_struct_print(char *nameout, int ia, float Amplitude, int Nstep,
		      int PRINT, int DIR, struct axe *axes, int naxes,
		      struct Normal_Mode NM,
		      atom *atoms, int natoms,
		      struct residue *seq, int nres,
		      struct bond *bonds);
int Move_struct_confchange(atom *atoms, int natoms,
			   float *rmsd_fin, char *nameout, int ia, 
		            int PRINT, struct axe *axes, int naxes,
		            struct Normal_Mode *NM,
		            struct residue *seq, int nres,
                            atom *atoms2, int N_ref,
	                    int *atom_num1, int *atom_num2,
			    struct Para_simul Para_simul,
			    struct bond *bonds);
int Move_all_modes(char *nameout,
                   struct axe *axes, int naxes, struct Normal_Mode *NM,
                   atom *atoms, int natoms, struct residue *seq, int nres,
		   atom *atoms2, int N_ref,
		   int *atom_num1, int *atom_num2,
		   int nstep, struct Para_simul Para_simul,
		   struct bond *bonds);
void Print_PDB(FILE *file_out, atom *atoms, int natoms, float *atom_str3,
	       struct residue *seq, int k, float rmsd); //
void Write_ref_coord_d(double *coord, int N_ref, atom *atoms, int *atom_num);
float Test_buildup(struct bond *bonds, atom *atoms, int natoms, int naxes);
float Switch_bonds(struct bond *bonds, int naxes, char *prots,
		   struct ali_atoms ali_a,
		   atom *atoms1, int natoms1, struct residue *seq1,
		   atom *atoms2, int natoms2, struct residue *seq2);
void Change_kin(int *first_kin, int *last_kin, struct axe *axes,
		int naxes, int natoms1, struct Reference Ref1);
void Change_kin_back(struct axe *axes, int naxes,
		     int *first_kin, int *last_kin);
void Standardize_bonds(struct bond *bonds, atom *atoms,
		       int natoms, char *nameout1,
			struct axe *axes, int naxes,
		       struct residue *seq, int nres, 
		       struct chain *chains, int Nchain, 
		       double Tors_fluct,
		       struct Para_simul Para_simul);

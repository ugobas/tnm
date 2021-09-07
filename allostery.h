void Predict_allostery(struct Normal_Mode NM, atom *atoms,
		       struct Reference Ref,
		       struct interaction *Int_list, int N_int,
		       struct residue *seq, char *nameout,
		       float *Confchange,
		       //int mode,
		       char *pdb, char *chain,
		       int Nres, char *SITES, int anharmonic);
void Get_interaction_list(int **clist, int **cnum, int *nc,
			  atom *atoms, int *atomres, int Na,
			  struct interaction *Int_list, int N_int, int NCMAX);
int Extract_atoms(int *iref, int *iatom, int *atomres,
		  struct residue *seq, float *mass,
		  atom *atoms, struct Reference Ref, char *ANAME);

extern int ALL_PAIRS;  // Print all pairs?
extern float SIGMA; // if(ALL_PAIRS==0) couplings are printed only for
// Coupling_ij > SIGMA*std_dev; default SIGMA=1.0
extern int STRAIN; // Compute strain profile, inspired by
  // Balabin, Yang, Beratan PNAS 106:14253 2009 ?
extern int PRINT_COV_COUPLING; // Print covariance coupling?
extern int PRINT_DEF_COUPLING; // Print deformation coupling?
extern int PRINT_DIR_COUPLING; // Print directionality coupling?
extern int PRINT_COORD_COUPLING; // Print coordination coupling?
extern int PRINT_SIGMA_DIJ; // Print variance of DIJ distance
extern char PROF_TYPE;     // Type of profile A, E or P

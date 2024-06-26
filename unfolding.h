struct inter_new{
  int i1, i2;
  float f0; // v(r)
  float f1; // v'(r)/r
  float f2; // v''(r)
  float f2c; //(v''-v'/r)/r^2
  float r; // |rj-ri|
  float rji[3]; // rj-ri
  float crji[3]; // rj X ri
  int nat;
};
extern struct interaction *Int_list_all;
extern int N_int_all;
extern float d_CA_high;


int Unfolding(struct bond *bonds_in, atom *atoms_in, int natoms_in,
	      struct axe *axes_in, int naxes_in, int nmain_in,
	      struct chain *chains_in, int Nchain_in,
	      struct residue *seq_in, int nres_in,
	      struct interaction *Int_list_in, int N_int_in,
	      struct interaction **Int_KB_in, int NA_in,
	      float K_TORS, char *file_pdb_in);
float Energy_anharmonic_new(struct inter_new *New_inter,
			    int *N_int_new, int N_int_max,
			    float *coord, atom *atoms, int natoms,
			    struct residue *seq, int nres,
			    struct interaction *Int_list, int N_int,
			    struct interaction **Int_KB, int NA,
			    struct axe *axe, int naxes,
			    double *delta_phi, int K_tors, int sec_der);
void Move_atoms(atom *atoms, int natoms, float *coord);

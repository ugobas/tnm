// Parameters for interactions

int S_TYPE;        // Type of shadow interaction
int ONEINT;        // Only one interaction per atom vs. residue
float S_THR;       // Screening parameter
int HNM;           // Calpha interactions with Hissen parameters
int NOCOV;         // Covalent neighbors do not interact
char REF[20];      // Reference atoms
int N_RESRES;      // Number of interactions retained for each residue pair

void Assign_interactions_SB(struct interaction *Int_list, int N_int,
			    struct interaction **Int_KB, int *aa_seq);
void Assign_interactions_KB(struct interaction **Int_KB, int Na,
			    int POW, float EXPO, float k0);
void Print_interactions(struct interaction *Int_list, int N_int,
			struct interaction **Int_KB,
			int *aa_seq, char *AA_code, int NA, char *name);
void Compute_sec_der(struct interaction *Int_list, int N_int, char *nameout);

/*******************************  Go model  *********************************/
int Interactions(struct interaction *Int_list, char *atom_type, float thr,
		 atom *atoms, int N_atoms, int N_res);
int Interactions_hnm(struct interaction *Int_list, char *atom_type, float thr,
		     atom *atoms, int N_atoms, int N_res);
int Interactions_all(struct interaction *Int_list, float thr,
		     atom *atoms, int N_atoms, int N_res);
int Interactions_all_CB(struct interaction *Int_list, float thr,
			atom *atoms, int N_atoms, int N_res);
int Interactions_HB(struct interaction *Int_list, float thr,
		    atom *atoms, int N_atoms, int N_res,
		    float thr_HB, float cos_HB, float ene_HB);

int Backbone(char *name);
void Compute_interactions(int *N_int, struct interaction **Int_list,
			  char *INT_TYPE, atom *atoms1, int natoms1, int nres1, char *nameout);
/*int Interactions_shadow(struct interaction *Int_list, float THR,
			atom *atoms, int N_atoms,
			int TYPE, float S_THR, int NOCOV);*/
int Screened_Interactions(struct interaction *Int_list,
			  float THR, atom *atoms, int N_atoms, int nres,
			  int TYPE, float DELTA_SCR);

void Print_contact_matrix(struct interaction *Int_list, int N_int,
			  atom *atoms, struct residue *res,
			  char *name, char *pdb);
int Covalent(int res1, int res2, atom *atom1, atom *atom2);
int Bond(int res3, int res1, float d2_13, int res2, float d2_23);
void Scale_interactions(char *INT_TYPE, atom *atoms, int natoms, int nres);
struct interaction *Check_pair(int *new1, int res1, int *res2, int nres,
			       struct interaction **ini,
			       struct interaction *end);

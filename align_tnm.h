#ifndef ALIGN_TNM
#define ALIGN_TNM

extern int MIN_ALI_LEN; // Chains are aligned if L>MIN_ALI_LEN, if not peptides

struct Ali_score{
  float seq_id, aligned, psi, mamm_score;
  int ngaps;
  int *alignres;
  int mammoth;
};


int Align_chains(int *Nmut, char *AAwt, int *Posmut, char *AAmut,
	         struct Ali_score *ali, int *last_ali_res,
		 int *Nchain, float SEQID_THR,
		 struct chain *chains1, int Nchain1,
		 struct residue *seq1, char *pdbid1,
		 struct chain *chains2, int Nchain2,
		 struct residue *seq2, char *pdbid2);
int Align_atoms(atom *atoms1, atom *atoms2, // Output:
		struct chain *ch1, struct chain *ch2); 	// Input:
void Purge_not_aligned(struct chain *chains1, atom *atoms1, int *natoms1,
		       struct chain *chains2, atom *atoms2, int *natoms2);
float Examine_confchange(struct Tors *Diff,
			 struct bond *bonds, struct axe *axe,
			 struct chain *chains, int Nchain,
			 char *nameout, char *name1, char *namepdb,
			 atom *atoms1, float *coord_ref1,
			 struct residue *seq1, int nres1,
			 int N_diso1, int natoms1,
			 atom *atoms2, float *coord_ref2,
			 struct residue *seq2, int nres2,
			 int N_diso2, int natoms2,
			 struct Reference Ref, //struct Jacobian *J,
			 struct ali_atoms ali_atoms,
			 struct interaction *Int_list1, int N_int1,
			 char *INT_TYPE, float S0, int N_modes,
			 struct Normal_Mode NM, int nprint,
			 struct Ali_score ali,
			 struct Para_simul Para_simul,
			 //float *diff_phi,  float *greedy_phi,
			 int nstruct, float *B_TNM, int *outlier_tors,
			 char *name2, int PRINT_CHANGE,
			 int anharmonic);
void Print_modes_confchange(char *name2, char *confchange,
			    struct Normal_Mode NM, struct Tors Diff,
			    int ANM, float M_sqrt, float rmsd, int anharmonic,
			    float xkappa);
void Allocate_Normal_modes(struct Normal_Mode *NM,
			   int N_modes, int N_axes, int N_Cart);

#endif

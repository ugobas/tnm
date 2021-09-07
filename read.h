int Read_coord(char *pdb_name, int *nmr, struct residue *seq, atom *atoms,
	       char *chain_to_read, int *ANISOU);
int Count_residues(char *pdb_name, char *chain_to_read);
int Get_amm(char *amm, int *exo, char *res_type_old, int hetatm, int n_exo);
void Delete_tmp_file();
FILE *Open_compressed_file(char *name, int *Compression);
int Get_compression(char *name);
void Code_3_2(char *res2, char *res3);
char Code_3_1(char *res);
void Align_seqres(int *ali_seqres, char *seqres, int N_seqres,
		  char *seq, int nres);
int Read_seqres(char **seqres, int *N_seqres, int NCHAIN,
		char *file_pdb, char *chain_to_read);
void Find_PDB_order(atom *atoms, int N_atoms);

int L_MAX;


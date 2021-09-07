
int Read_sites(int **site_res, int *nres_site, char **chain_res,
	       char **description, char *file_site,
	       int NMAX, char *SITES, char *pdb, char *chain, int Nchain,
	       struct residue *seq, int Nres, int PRINT_SITE);
int Read_sites_PDB(int **site_res, char **chain_res,
		   int *nres_site, char **description,
		   char *pdb, char *chain, struct residue *seq,
		  int Nres, int NMAX, int SMAX);
int Read_sites_file(int **site_res, char **chain_res,
		    int *nres_site, char **description,
		    char *file, char *chain, struct residue *seq,
		    int Nres, int NMAX, int SMAX);
int Read_site(char *resnum, char *resnam, char *reschain,
	      char *string, char *chain);
int Site_overlap(int *isite, int n1, int *jsite, int n2);

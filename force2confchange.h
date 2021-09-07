int Force2Confchange(float *coord_force, char *FILE_FORCE,
		     int naxe1, int N_ref, atom *atoms1, int natoms1,
		     int *atom_num1, float **Jacobian_ar, float *eigen_value,
		     float **Tors_mode, int *select, struct axe *axe1,
		     char *chain1, int nres1, struct residue *seq1,
		     struct bond *bonds);

int Interactions_CA(struct interaction *Int_list, float thr,
		    atom *atoms, int natoms, int N_res, float *coord);
void Initialize_energy(struct interaction *Int_list, int N_int,
		       float *r, atom *atoms, int natoms, int nres);

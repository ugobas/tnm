int Compute_kinetic(struct Jacobian *J,
		    struct axe *axe1, int naxe1,
		    atom *atoms1, int natoms1,
		    struct Reference Ref, int SQR);
void  Set_rot_shift(struct axe *axes, int N_axes, atom *atoms, int N_atoms);
void Inertia_tensor(double **inertia, double **corr_sum, double *r_sum,
		    double M_sum);
float Sum_inertia(double *mr_tot, double **corr_sum,
		  atom *atoms, struct Reference Ref,
		  int ini_ref, int end_ref);
float Empty_inertia(double *mr_tot, double **corr_sum);
void Center_atoms(atom *atoms, int N_atoms, struct Reference Ref);
void Reorient_axes(atom *atoms, int N_atoms, struct Reference Ref);
void Eckart_main(struct axe *axes, int N_axes, atom *atoms, int N_atoms,
		 struct Reference Ref);

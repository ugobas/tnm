int Compute_kinetic(struct Jacobian *J,
		    struct axe *axe1, int naxe1,
		    atom *atoms1, int natoms1,
		    struct Reference Ref, int SQR);
void  Set_rot_shift(struct axe *axes, int N_axes, atom *atoms, int N_atoms);

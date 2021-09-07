
struct bond *Set_bonds_topology(int N_atoms, atom *atoms, struct residue *seq);
int Set_bonds_measure(struct bond *bonds, int N_atoms, atom *atoms);
int Set_bonds_prot(struct bond *bonds, int N_atoms,
		   atom *atoms, struct axe *axes, int naxes);
void Build_up(struct bond *bonds, int N_atoms,
	      float *d_phi, int N_axes);
void Internal_coordinates(double *int_coord, int N_axes,
			  struct bond *bonds, int N_atoms);
void Trajectory(atom *atoms, int natoms, struct axe *axes, int naxe,
		float **Tors_mode, int imode, float fact, int nmove,
		struct bond *bonds);
struct axe *Set_DegofFreed(int *naxe,     // Number of degrees of freedom
			   int *nmain,    // Number of mainaxes
			   int *nskip,    // Skipped mainaxes
			   int *N_diso,   // Number of disorder gaps
			   struct bond *bonds,  // Covalent topology
			   atom *atoms,   // Pointer to protein atoms
			   int natoms,    // Number of protein atoms
			   int *atom_ref, // Reference atoms
			   int N_ref,     // Pointer to ligand atoms
			   struct residue *res,  // residues
			   int nres,      // Number of residues
			   struct chain *chains, // chains
			   int Nchain,    // Number of chains
			   struct interaction *Int_list, // interactions
			   int N_int,     // Number of interactions
			   int MIN_INT_MAIN, // Min. interactions per dof
			   int MIN_INT_SIDE, // Min. interactions per dof
			   int OMEGA,     // Are Omega angles stored?
			   int SCHAIN,   // Are side chains used?
			   int PSI);     // Are psi angles used?
int D_omega_thr; // Omega unfrozen if D_omega>D_omega_thr


/*void Set_atom_axe(atom *atoms, struct chain *chains, int Nchain,
  struct axe *axes);*/

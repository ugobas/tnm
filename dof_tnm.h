struct axe *Set_DegofFreed(int *naxe,     // Number of degrees of freedom
			   int *nmain,    // Number of mainaxes
			   int *nskip,    // Skipped mainaxes
			   int *N_diso,   // Number of disorder gaps
			   struct rigid **Rigid_dof, // Chains, axe index
			   int *N_rigid,  // N. rigid body dofs
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

void Read_PDB(int *nres,   // Number of amino acid residues
	      struct residue **seq, // Where ligand and aa res are stored
	      int *n_lig,  // Number of ligand residues
	      struct residue **ligres, // Pointer to first ligand residue
	      char *chain, // Chain identifiers
	      int *ANISOU, // Anisotropic B factors?
	      int *nmr,    // NMR structure?
	      int *natoms, // Number of protein atoms
	      atom **atoms,// Pointer to protein atoms
	      int *na_lig, // Number of ligand atoms
	      atom **ligatoms, // Pointer to ligand atoms
	      struct chain **chains, // Pointer to protein chains
	      int *Nchain, // Number of protein chains (1-Nchain)
	      float *TEMP, // Temperature of the experiment, if given
	      char *pdbid, // PDB identifier
	      char *file_pdb); // INPUT: Path to PDB file


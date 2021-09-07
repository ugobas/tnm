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
	      char *pdbid, // PDB identifier
	      char *file_pdb); // INPUT: Path to PDB file


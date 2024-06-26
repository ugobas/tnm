#ifndef __INI_BUILDUP
#define __INI_BUILDUP

struct bond *Set_bonds_topology(int N_atoms, atom *atoms, struct residue *seq);
int Set_bonds_measure(struct bond *bonds, int N_atoms, atom *atoms);
int Set_bonds_prot(struct bond *bonds, int N_atoms,
		   atom *atoms, struct axe *axes, int naxes);
void Build_up(struct bond *bonds, int N_atoms,
	      double *d_phi, int N_axes);
void Internal_coordinates(double *int_coord, int N_axes,
			  struct bond *bonds, int N_atoms);
void Trajectory(atom *atoms, int natoms, struct axe *axes, int naxe,
		float **Tors_mode, int imode, float fact, int nmove,
		struct bond *bonds);

extern int D_omega_thr; // Omega unfrozen if D_omega>D_omega_thr


/*void Set_atom_axe(atom *atoms, struct chain *chains, int Nchain,
  struct axe *axes);*/

#endif

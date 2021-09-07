// #define L_MAX 40000       /* Maximal length of a chain */

#define CHMAX 800          // Maximum number of chains

typedef struct {
  // Read from PDB
  double r[3];
  char name[4];
  float B_factor;
  float anisou[3][3];
  float occupancy;
  // main chain or side chain?
  int main;
  // Connections
  int i_next, i_num;
  int ali;
  // Degrees of freedom
  //int first_dof_main; // first  main-chain DoF
  int last_dof_main;  // last  main-chain DoF
  int last_dof_side;  // last  side-chain DoF
  //int first_dof_side; // first side-chain DoF
  // Interactions
  int N_int;
  // Residue
  int res;
  char aa;
  int chain;
} atom;

struct bond{
  atom *atom;
  double r[3];           // Coordinates of representative atom
  int i_atom;
  //char *name;
  //int chain;  
  //int res;
  //char type; // m= main chain, s= side chain, r= rigid body
  // Degrees of freedom
  struct bond *previous;
  int terminal;       // Is the bond at terminus of the molecule?
  int len;           // DoF associated with bond length
  int angle;         // DoF associated with bond angle
  int torsion;       // DoF associated with torsion, previous rot_axe
  //int stored;        // Rotation around the axis already considered? 
  // Internal coordinates
  double dr[3];          // dr - previous->dr
  double l;             // bond length |dr|;
  double r_cos_theta;   // l*cos(theta)
  double phi;           // arccos(cos_phi)
  short d_phi;          // Difference of phi in degrees
  double dy[3];          // cross product dr X previous->dr
  double dx[3];          // norm to dr in plane dr dr1
  double dz[3];          // Direction of dr
  //float ax, ay, az;     // Coordinates
};

struct axe{
  char type;                  // l=length a=angle other=torsion
  int rigid;
  struct bond *bond;          // Bond modified by the dof
  //struct bond *axe;           // axe corresponding to the dof
  //atom *atom;               // Last atom
  //int atom1, atom2;         // atoms that define the axe
  //int res;                    // Residue of last atom
  //int chain;                  // chain of last atom
  int first_atom;             // first and last atoms moved by the axe
  int last_atom;
  int first_kin;              // first reference atom moved by the axe
  int last_kin;
  int first_main;             // First main chain degree of freedom
  int last_main;              // Last main chain degree of freedom
  int first_side;             // First side chain degree of freedom
  int previous;
  double offset[3];           // origin of the axis
  double rot[3];              // rotation
  double shift[3];            // translation (remember, minus sign)
  double vers[3];             // axe versor
  double global_shift[3];     // Rigid body translation
  double global_rot[3];       // Rigid body rotation
  double local_shift[3];      // Rb translation - shift
  double local_rot[3];        // Rb rotation + rot
  double mass;
  double MR[3];
  double I[3][3];
};

struct interaction{
  int i1, i2; // atom2 > atom1
  int res1, res2;
  float r1, r2, r3; // thresholds
  double r0, sec_der;
  double A0, A1, A2, A4; // for r<r1; A4=repulsion A2=attraction
  double B0, B1, B2, B4; // for r2<r<r3;
  double a0, a1, a2, a3, a4, a5; // r1<r<r2
  double b3, b4, b5; // r2<r<r3;
  double E_min;
  int type;
};

struct residue{
  atom *atom_ptr;
  int n_atom;
  char chain;
  char pdbres[6];
  char amm;
  short i_aa;
  short exo;
  short type;     // Poly: AA=1, RNA=2, DNA=3, other=0, Single: *-1
  short n_cont;
  float asa;
  //float PE;
  //float B_factor;
  //float B_NM;
};

struct chain{
  int ini_atom;
  int natoms;
  int ini_main;
  int mainaxes;
  int ini_side;
  int nsides;
  int ini_res;
  int nres;
  //int nref;
  int ini_cart;
  int ncart;
  int type;   // Polymer: AA=1, RNA=2, DNA=3, other=0; Monomer: -1
  int last_atom;
  int last_ref;
  char label;
  int *alignres;
  int match;
  struct residue *res;
  char *seq;     // Only structured residues
  char *seqres;  // SEQRES record
  int N_seqres;
  int *ali_seqres;
};

struct molecule{
  atom *atoms;
  int natoms;
  double *atom_str;
  int N_ref, N_cart;
  struct residue *seq;
  int nres;
  int N_diso;
  struct chain *chains;
  char chain[CHMAX];
  int Nchain;
  struct interaction *Int_list;
  int N_int;
  struct Reference *Ref;
  struct axe *axe;
  int naxe;
  char name[100];
  char pdbid[100];
  int ANISOU, nmr;
};

struct ali_atoms{
  int *ali1;
  int *ali2;
  float *mass;
  int N_ref;
  int N_cart;
};

#define FILE_MSG   "message.out"


int Verbose;
char cont_type;
float cont_thr, cont_thr2;


short **Compute_map(int N_res, struct residue *seq, int *N_cont);
void Average_B(struct residue *seq, int N);
int Is_atom_main(char *name, int seq_type);
int Code_AA(char res);
char *AA_code;

#include "coord.h"
#include "buildup.h"
#include "nma_para.h"
#include "tnm.h"
#include "nma.h"
#include "dof_tnm.h"
#include "read.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

// Arrays limiters
// #define ATOM_MAX 50000
#define ATOM_MAX 50 // Number of atoms per residue
#define DEBUG 0

/*********************************************************************/
// Setting degrees of freedom
int Set_dof(char type, struct axe **axe, int *num_dof,
	    struct bond *bond, int *ndof, int *oldref, int iref);
int Set_side_chains(struct axe *axes, struct bond *bonds, atom *atoms,
		    int natoms, int *refatom, int *atom_ref,
		    int nmain, int MIN_INT);
int Order_atoms(atom *atoms, struct residue *seq, int Nres,
		int *Nchain, int natoms);
int Count_interactions(atom *atoms, int natoms,
		       struct interaction *Int_list, int N_int);
int Set_chains(struct chain *chains,
	       atom *atoms, int N_atoms,
	       struct residue *res, int Nchain);
int Copy_atom(atom **atom2, atom *atom1, int ires, int i_pdb, int Nc, int i);
int Store_atoms(atom **atom2, char BRANCH,
		atom *res_atom, int n_atom,
		int ires, int Nchain, int N_pdb);
struct bond *Find_bond_l(struct axe *axe, struct bond *bonds);
struct bond *Find_bond_a(struct axe *axe, struct bond *bonds);
struct bond *Find_bond_t(struct axe *axe, struct bond *bonds, int *i_bond);

/*********************************************************************/
// Auxiliary
void GetPdbId(char *pdb_file_in, char *pdbid);
float Distance_square(float *r1, float *r2);
extern int Find_atom(atom *atoms, int *i1, int res, int Natoms, char *type);

/*********************************************************************/
/*********************************************************************/
// Start:

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
	      char *file_pdb) // INPUT: Path to PDB file
{
  /****************** Reading atoms: ************************/
  int numres=Count_residues(file_pdb, chain); L_MAX=numres;
  int numatoms=ATOM_MAX*numres;
  atom *atom_read=malloc(numatoms*sizeof(atom));
  *seq=malloc(numres*sizeof(struct residue));
  *nres=Read_coord(file_pdb, nmr, *seq, atom_read, chain, ANISOU);
  if(*nres<=0){
    printf("WARNING, file %s, no residue found\n", file_pdb); return;
  }
  GetPdbId(file_pdb,pdbid); if(*chain==' ')*chain='_';

  // *n_lig=Read_ligand(file_pdb, *seq+*nres, atom_read+*natoms, na_lig);
  // Ligand atoms are stored at the end

  /****************** Storing atoms: ************************/
  int i, j;
  *natoms=0;
  for(i=0; i<*nres; i++)*natoms+=(*seq)[i].n_atom;
  *atoms=malloc((*natoms+*na_lig)*sizeof(atom));
  printf("PDB %s chain %s nres=%d natoms=%d\n",
	 pdbid,chain,*nres,*natoms);

  // Store protein atoms
  *natoms=Order_atoms(*atoms, *seq, *nres, Nchain, *natoms);
  // Store ligand atoms
  atom *atom1=atom_read+*natoms, *atom2=*atoms+*natoms;
  for(j=*natoms; j<(*natoms+*na_lig); j++){
    *atom2=*atom1; atom2->chain=-1;
    atom2->i_num=j; atom2->i_next=j+1;
    atom1++; atom2++;
  }

  /**************** Set chains *********************/
  *chains=malloc((*Nchain)*sizeof(struct chain));
  (*chains)->ini_atom=0;
  Set_chains(*chains, *atoms, *natoms, *seq, *Nchain);
  for(i=0; i<*natoms; i++)(*atoms)[i].aa=(*seq)[(*atoms)[i].res].amm;
  printf("%d chains found\n", *Nchain);

  return;
}

struct axe *Set_DegofFreed(int *naxe,     // Number of degrees of freedom
			   int *nmain,    // Number of mainaxes
			   int *nskip,    // Skipped mainaxes
			   int *N_diso,   // Number of disorder gaps
			   struct bond *bonds,  // Covalent topology
			   int natoms,    // Number of protein atoms
			   atom *atoms,   // Pointer to protein atoms
			   int *atom_ref, // Reference atoms
			   int N_ref,     // Pointer to ligand atoms
			   int nres,      // Number of residues
			   struct residue *res,  // residues
			   struct chain *chains, // chains
			   struct interaction *Int_list, // interactions
			   int N_int,     // Number of interactions
			   int MIN_INT_MAIN, // Min. interactions per dof
			   int MIN_INT_SIDE, // Min. interactions per dof
			   int OMEGA,     // Are Omega angles stored?
			   int SCHAIN)    // Are side chains used?
{
  // Allocate degrees of freedom
  int N_axe_max=3*nres;
  if(OMEGA)N_axe_max+=nres;
  if(SCHAIN)N_axe_max+=3*nres;
  struct axe *axes=malloc(N_axe_max*sizeof(struct axe));

  // By default, degrees of freedom are frozen
  int ibond; struct bond *bond=bonds;
  for(ibond=0; ibond<natoms; ibond++){
    bond->len=-1; bond->angle=-1; bond->torsion=-1; bond++;
  }

  // Set types of bonds (aa, dna, rna, ligand)
  for(ibond=0; ibond<natoms; ibond++){
    bond->type=res[atoms[bond->i_atom].res].type;
  }

  // Prepare reference atoms
  int *refatom=malloc(natoms*sizeof(int)), i, iref=0;
  int refmax=N_ref-1;
  for(i=0; i<natoms; i++){
    if(atom_ref[iref]<i)iref++;
    if(iref<=refmax){refatom[i]=iref;}
    else{refatom[i]=-1;}
  }

  // Store axes
  /* A new degree of freedom is created only if there is at least one
     reference atom for 3 new degrees of freedom, and if the number
     of interactions is at least
  */
  int N_axes=0,r_ini=-1, r_last=nres-1;
  int oldref=1, ndof=0;
  struct axe *axe=axes;
  int move_atom; //last_move_atom=-1;

  for(ibond=1; ibond<(natoms-1); ibond++){
    bond=bonds+ibond;
    if(bond->i_atom<0)break;
    if((bond->previous==NULL)||(bond->terminal))continue;
    if((bond->type!=1)&&(bond->chain==bond->previous->chain))continue;
    if(r_ini<0)r_ini=bond->res;
    char type='q';
    if(strncmp(bond->name,"N ", 2)==0){
      if((strncmp(bond->previous->name,"C ", 2)==0)&&
	 (bond->res==(bond->previous->res+1))){
	if(OMEGA){
	  type='o';   // omega angle
	}else if(bond->chain==bond->previous->chain){continue;}
      }else{type='q';}         // Pseudo-axe after disordered loop
    }else if(strncmp(bond->name,"CA", 2)==0){
      if((strncmp(bond->previous->name,"N ", 2)==0)&&
	 (bond->res==bond->previous->res)){
	type='f'; // phi
	if(bond->res==r_ini)continue;
	if((bond->res==r_last)&&(ALL_AXES==0))continue;
      }else{type='q';}
    }else if(strncmp(bond->name,"C ", 2)==0){
      if((strncmp(bond->previous->name,"CA", 2)==0)&&
	 (bond->res==bond->previous->res)){
	type='p'; // psi
	if(bond->res==r_last)continue;
	if((bond->res==r_ini)&&(ALL_AXES==0))continue;
      }else{type='q';}
    }else{
      continue;
    }

    if((bond->chain!=bond->previous->chain)||
       (bond->type!=bond->previous->type)){
      type='c';
    }else{
      if(Distance_square(bond->r, bond->previous->r) > 5)type='q';
      if(type=='q'){
	(*N_diso)++;  // Disordered residue
      }
    }

    if((type=='q')||(type=='c')){
      iref=refatom[bond->i_atom]; if(iref<0)continue;
      Set_dof('l', &axe, &N_axes, bond, &ndof, &oldref, iref);
      Set_dof('a', &axe, &N_axes, bond, &ndof, &oldref, iref);
      //last_move_atom=bond->i_atom;
      //if(type=='c')continue;
      type='t';
    }
    move_atom=bond->i_atom+1;
    iref=refatom[move_atom];
    if(iref<0)continue;

    Set_dof(type, &axe, &N_axes, bond, &ndof, &oldref, iref);

  }


  // Exclude dof with few interactions
  int n_skipped=0;
  if(MIN_INT_MAIN>0){
    Count_interactions(atoms, natoms, Int_list, N_int);
    int n=0, last=N_axes-1, last_atom, a, num_int=0;
    for(a=0; a<N_axes; a++){
      axe=axes+a; num_int=0;
      if(a<last){
	struct axe *next=axe+1;
	if(next->type=='l'){if(a+1<last)next++;}
	else if(axe->type=='a'){if(a+1<last)next++;}
	last_atom=next->first_atom;
      }else{
	last_atom=natoms;
      }
      for(i=axe->first_atom; i< last_atom; i++)num_int+=atoms[i].N_int;
      if(num_int >= MIN_INT_MAIN){ //||(axe->type=='l')||(axe->type=='a')){
	if(a!=n)axes[n]=*axe;
	n++;
      }else{
	n_skipped++;
      }
    }
    N_axes=n;
  }
  (*nmain)=N_axes;
  (*naxe)=N_axes;
  (*nskip)=n_skipped;
  printf("%d axes over %d skipped since <= %d interactions\n",
	 n_skipped, N_axes+n_skipped, MIN_INT_MAIN);

  // Print degrees of freedom
  int ig=1, chain=atoms[0].chain;
  for(i=0; i<N_axes; i++){
    axe=axes+i;
    if(atoms[axe->atom2].chain!=chain){
      printf("-ChainBreak-"); chain=atoms[axe->atom2].chain; ig=0;
    }else if((axe->type=='l')||((axe->type=='a')&&ig)){
      printf("-Gap %2s%d-%2s%d-",
	     atoms[axe->atom1].name, atoms[axe->atom1].res,
	     atoms[axe->atom2].name, atoms[axe->atom2].res);
      ig=0;
    }else{
      ig=1;
    }
    printf("%c", axe->type);
  }
  printf("\n");


  // Set last atom for each main axe
  axe=axes;
  int last_atom=natoms-1, last_ref=N_ref-1;
  for(i=0; i<N_axes; i++){
    axe->last_atom=last_atom;
    axe->last_kin=last_ref;
    axe->first_side=-1;
    axe++;
  }

  // First and last degree of freedom for each atom
  int iatom=axes[0].first_atom, nmax=N_axes-1; ndof=-1;
  for(i=0; i<natoms; i++){
    if((i >= iatom)&&(ndof < nmax)){
      ndof++; if(ndof<nmax)iatom=axes[ndof+1].first_atom;
    }
    atoms[i].last_dof_main=ndof;
    atoms[i].last_dof_side=-1;
  }

  // Association between main chain degrees of freedom and bonds
  for(i=0; i<N_axes; i++){
    axe=axes+i; 
    if(axe->type=='l'){
      bond=Find_bond_l(axe, bonds);
      bond->len=i;
    }else if(axe->type=='a'){
      bond=Find_bond_a(axe, bonds);
      bond->angle=i;
    }else{
      int i_bond=-1;
      while(1){
	bond=Find_bond_t(axe, bonds, &i_bond);
	if(bond==NULL)break;
	bond->torsion=i;
      }
    }
  }

  if(DEBUG){
    axe=axes;
    for(i=0; i<N_axes; i++){
      printf("%c  %3d %s %3d %s  %3d %3d\n", axe->type,  axe->atom1,
	     atoms[axe->atom1].name, axe->atom2, atoms[axe->atom2].name,
	     axe->first_atom, axe->last_atom);  axe++;
    }
    //exit(8);
  }

  // axes, n_axe for each chain
  //~ struct chain *chp=chains;
  chains[0].ini_axe=0;  chains[0].mainaxes=N_axes;
  /*int ichain=0, n=0, chain=res[0].chain;
    for(i=0; i<N_axes; i++)axes[i].chain=ichain;*/


  if(SCHAIN){
    (*naxe)+=Set_side_chains(axes+N_axes, bonds, atoms, natoms,
			     refatom, atom_ref, *naxe, MIN_INT_SIDE);
  }
  free(refatom);

  return(axes);
}

int Set_dof(char type, struct axe **axe, int *num_dof,
	    struct bond *bond, int *ndof, int *oldref, int iref)
{
  /* Before updating the reference atom, maximum 3 degrees of freedom
     can be accepted   */
  if((*ndof >= 2)&&(iref <= *oldref)){
    printf("="); return(0);
  }
  if(iref > *oldref){
    *oldref=iref; *ndof=0;
  }else{
    (*ndof)++;
  }
  //if(type=='l')printf("*");
  (*axe)->type=type;
  (*axe)->atom1=bond->previous->i_atom;
  (*axe)->atom2=bond->i_atom;
  (*axe)->chain=bond->chain;
  (*axe)->res=bond->res;
  if(type=='l'){
    (*axe)->first_atom=bond->i_atom;
  }else{
    (*axe)->first_atom=bond->i_atom+1;
  }
  (*axe)->first_kin=iref;
  (*axe)++;
  (*num_dof)++;
  return(1);
}

int Order_atoms(atom *atoms, struct residue *seq, int Nres,
		int *Nchain, int natoms)
{
  int i, ires, N_pdb=0;
  char chain=seq[0].chain;
  int N_atoms=0;

  *Nchain=0;
  for(ires=0; ires<Nres; ires++){
    atom *res_atom=seq[ires].atom_ptr;
    int n_atom=seq[ires].n_atom;

    if(chain!=seq[ires].chain){
      chain=seq[ires].chain; (*Nchain)++;
    }

    // Store main chain before side chain
    int i_N=-1, i_CA=-1, i_C=-1, i_O=-1, i_OT=-1;
    int i_HN[3], i_HA[3], n_HN=0, n_HA=0, j;
    if(HYD)for(j=0; j<3; j++){i_HN[j]=-1; i_HA[j]=-1;}
    atom *atom1=res_atom;
    for(i=0; i<n_atom; i++){
      if(strncmp(atom1->name, "N ", 2)==0){i_N=i;}
      else if(strncmp(atom1->name, "CA", 2)==0){i_CA=i;}
      else if(strncmp(atom1->name, "C ", 2)==0){i_C=i;}
      else if(strncmp(atom1->name, "O ", 2)==0){i_O=i;}
      else if((strncmp(atom1->name, "OT", 2)==0)||
	      (strncmp(atom1->name, "OX", 2)==0)){i_OT=i;}
      else if(atom1->name[0]=='H'){
	   if((atom1->name[1]==' ')||(atom1->name[1]=='1')||
	      (atom1->name[1]=='2')||(atom1->name[1]=='3')){
	     atom1->name[2]=atom1->name[1]; atom1->name[1]='N';
	     i_HN[n_HN]=i; n_HN++;
	   }else if(atom1->name[1]=='A'){
	     i_HA[n_HA]=i; n_HA++;
	   }
      }
      atom1++;
    }

    // Store main chain in the order N (HN) CA (HA)
    atom *atom2=atoms+N_atoms;
    int N_atom_res=0;
    N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i_N);
    if(HYD){
      for(j=0; j<3; j++)
	N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i_HN[j]);
    }
    N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i_CA);
    if(HYD){
      for(j=0; j<3; j++)
	N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i_HA[j]);
    }

    // Store side chain before C-O
    int N_side=n_atom;
    if(i_C>=0)N_side--; if(i_O>=0)N_side--; 
    if(N_atom_res==N_side)goto CO;
    N_atom_res+=Store_atoms(&atom2, 'B', res_atom,n_atom,ires,N_pdb,*Nchain);
    if(N_atom_res==N_side)goto CO;
    N_atom_res+=Store_atoms(&atom2, 'G', res_atom,n_atom,ires,N_pdb,*Nchain);
    if(N_atom_res==N_side)goto CO;
    N_atom_res+=Store_atoms(&atom2, 'D', res_atom,n_atom,ires,N_pdb,*Nchain);
    if(N_atom_res==N_side)goto CO;
    N_atom_res+=Store_atoms(&atom2, 'E', res_atom,n_atom,ires,N_pdb,*Nchain);
    if(N_atom_res==N_side)goto CO;
    N_atom_res+=Store_atoms(&atom2, 'Z', res_atom,n_atom,ires,N_pdb,*Nchain);
    if(N_atom_res==N_side)goto CO;
    N_atom_res+=Store_atoms(&atom2, 'H', res_atom,n_atom,ires,N_pdb,*Nchain);
    if(N_atom_res==N_side)goto CO;

    // Store other side chain atoms before C-O
    atom1=res_atom;
    for(i=0; i<n_atom; i++){
      if((strncmp(atom1->name, "N ", 2)!=0)&&
	 (strncmp(atom1->name, "CA", 2)!=0)&&
	 (strncmp(atom1->name, "C ", 2)!=0)&&
	 (strncmp(atom1->name, "O ", 2)!=0)&&
	 (strncmp(atom1->name, "OT", 2)!=0)&&
	 (strncmp(atom1->name, "OX", 2)!=0)&&
	 (atom1->name[1]!='N')&&
	 (atom1->name[1]!='A')&&
	 (atom1->name[1]!='B')&&
	 (atom1->name[1]!='G')&&
	 (atom1->name[1]!='D')&&
	 (atom1->name[1]!='E')&&
	 (atom1->name[1]!='Z')&&
	 (atom1->name[1]!='H')){
	N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i);
      }
      atom1++;
    }

    // Store CO
  CO:
    N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i_C);
    N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i_O);
    N_atom_res+=Copy_atom(&atom2, res_atom, ires, N_pdb, *Nchain, i_OT);

    N_atoms+=N_atom_res;
    N_pdb+=n_atom;
  }
  (*Nchain)++;

  if(N_atoms!=natoms){
    printf("ERROR, only %d atoms over %d have been stored\n",
	   N_atoms, natoms); exit(8);
  }

  if(DEBUG){
    for(i=0; i<N_atoms; i++)printf("%s ",atoms[i].name); printf("\n");
    //exit(8);
  }

  // Find next atom in original PDB order, i_next
  Find_PDB_order(atoms, N_atoms);
  return(N_atoms);
}

int Store_atoms(atom **atom2, char BRANCH,
		atom *res_atom, int n_atom,
		int ires, int Nchain, int N_pdb)
{
  atom *atom1=res_atom; int i, n=0;
  for(i=0; i<n_atom; i++){
    if(atom1->name[1]==BRANCH)
      n+=Copy_atom(atom2, res_atom, ires, N_pdb, Nchain, i);
    atom1++;
  }
  return(n);
}

int Copy_atom(atom **atom2, atom *atom1, int ires, int N_pdb, int Nc, int i)
{
  if(i<0)return(0);
  *(*atom2)=*(atom1+i);
  (*atom2)->res=ires;
  (*atom2)->i_num=N_pdb+i;
  (*atom2)->chain=Nc;
  (*atom2)++;
  return(1);
}

void Find_PDB_order(atom *atoms, int N_atoms)
{
  int i, i_next, res=atoms[0].res, ares=0, i2=0, i_min;
  int i_max=N_atoms+1000;
  atom *atom1=atoms, *atom2;
  for(i=0; i<N_atoms; i++){
    atom1=atoms+i;
    if(atom1->res > res){ares=i; res=atom1->res;}
    i_min=i_max; i_next=atom1->i_num+1;
    atom2=atoms+ares; i2=ares;
    while((i2<N_atoms)&&(atom2->res==res)){
      if((atom2->i_num>=i_next)&&(atom2->i_num<i_min)){
	i_min=i2; if(i_min==i_next)break;
      }
      atom2++; i2++;
    }
    if(i_min < i_max){atom1->i_next=i_min;}
    else if(i2<N_atoms){atom1->i_next=i2;}
    else{atom1->i_next=-1;}
  }
}

int Set_chains(struct chain *chains,
	       atom *atoms, int N_atoms,
	       struct residue *res, int Nchain)
{
  int i, ichain=0, n=0, ires=0, resold=-1, nres=0;
  char chain=res[0].chain;
  struct chain *chp=chains;

  for(i=0; i<Nchain; i++)chains[i].alignres=NULL;

  // ini_atom, ini_res etc.
  chp=chains;
  chp->ini_atom=0;
  chp->ini_res=0;
  chp->label=chain;
  for(i=0; i<N_atoms; i++){
    ires=atoms[i].res;
    if(ires!=resold){
      if(res[ires].chain!=chain){
	chain=res[ires].chain;
	chp->natoms=n; n=0;
	chp->nres=nres; nres=0;
	chp++; ichain++;
	chp->ini_atom=i;
	chp->ini_res=ires;
	chp->label=chain;
      }
      resold=ires; nres++;
    }
    atoms[i].chain=ichain; n++;
  }
  chp->natoms=n;
  chp->nres=nres;

  return(0);
}

float Distance_square(float *r1, float *r2){
  float d0=r1[0]-r2[0], d1=r1[1]-r2[1], d2=r1[2]-r2[2];
  return(d0*d0+d1*d1+d2*d2);
}

int Set_side_chains(struct axe *axes,
		    struct bond *bonds,
		    atom *atoms,
		    int natoms,
		    int *refatom,
		    int *atom_ref,
		    int Nmain,
		    int MIN_INT)
{
  int N_axes=0, ndof=0, numdof=0, sidedof=Nmain;
  int first_atom=-1, last_atom=0, last_ref=0; //first_ref=0
  int ibond, ires=-1, i;
  char type[1], aa=' ';
  struct axe *axe=axes;
  struct bond *bond;
  int iref, oldref=1;
  int last_main=0, first_side=0;
  int move_atom; // last_move_atom=0;

  for(ibond=1; ibond<(natoms-1); ibond++){
    bond=bonds+ibond;
    if(bond->i_atom<0)break;
    if(bond->previous==NULL)continue;

    if(bond->res!=ires){
      // Residue start
      ires=bond->res;
      aa=atoms[bond->i_atom].aa;
      if((aa=='G')||(aa=='A')||(aa=='P'))continue;
      // degrees of freedom
      first_side=sidedof;
      ndof=0; numdof=0;
      // Find first and last atom in side chain
      first_atom=-1;
      for(i=bond->i_atom; i<natoms; i++){
	atom *atom_i=atoms+i;
	if((atom_i->res>ires)||
	   //(strncmp(atom_i->name,"C ", 2)==0)||
	   //(strncmp(atom_i->name,"O ", 2)==0)||
	   (strncmp(atom_i->name,"OT", 2)==0)||
	   (strncmp(atom_i->name,"OX", 2)==0))
	  break;
	if((strncmp(atom_i->name,"N ", 2)==0)||
	   (strncmp(atom_i->name,"CA", 2)==0)||
	   (strncmp(atom_i->name,"CB", 2)==0))continue;
	if(first_atom<0)first_atom=i;
      }

      if(first_atom<0)continue;
      //last_move_atom=first_atom;

      last_atom=i-1;
      /*printf(" Res %c first_atom= %s %d last_atom= %s %d\n",
	     aa, atoms[first_atom].name, first_atom,
	     atoms[last_atom].name, last_atom);*/
      last_main=atoms[first_atom].last_dof_main;
      if(last_main<0)continue;
      //first_ref=refatom[first_atom];
      last_ref=refatom[last_atom];
      while(atoms[atom_ref[last_ref]].res>ires)last_ref--;
    }

    if((first_atom<0)||
       (last_main<0)||
       (bond->i_atom < first_atom)||
       (bond->i_atom > last_atom)){
      continue;
    }else if((aa=='G')||(aa=='A')||(aa=='P')){
      continue;     // exit if gly, ala or pro
    }else if((strncmp(bond->name,"CB", 2)==0)||
	     (strncmp(bond->name,"CG", 2)==0)){
      goto test_ref;
    }else if((aa=='F')||(aa=='W')||
	     (aa=='Y')||(aa=='H')){
      continue;  // if aromatic ring, exit
    }else if((bond->terminal)||
	     (strncmp(bond->name,"N ", 2)==0)||
	     (strncmp(bond->name,"CA", 2)==0)||
	     (strncmp(bond->name,"C ", 2)==0)||
	     (strncmp(bond->name,"O ", 2)==0)){
      continue;
    }

  test_ref:
    move_atom=bond->i_atom+1;
    iref=refatom[move_atom];
    if((iref<0)||(iref>last_ref)||
       (atoms[atom_ref[iref]].res!=ires))
      continue;

    sprintf(type, "%1d", numdof);
    if(Set_dof(type[0], &axe, &N_axes, bond, &ndof, &oldref, iref)){
      struct axe *axe1=axe-1;
      (axe1)->type='t';
      (axe1)->last_atom=last_atom;
      (axe1)->last_kin=last_ref;
      (axe1)->last_main=last_main;
      (axe1)->first_side=first_side;
      sidedof++;
      numdof++;
    }
  }

  if(MIN_INT>0){
    // Test number of interactions
    int n=0, a;
    for(a=0; a<N_axes; a++){
      axe=axes+a;
      int num_int=0;
      for(i=axe->first_atom; i< natoms; i++){
	if(atoms[i].res!=axe->res)break;
	num_int+=atoms[i].N_int;
      }
      if(num_int >= MIN_INT){
	if(a!=n)axes[n]=*axe; n++;
      }
    }
    N_axes=n;
  }

  // First and last degree of freedom for each atom
  int iatom=axes[0].first_atom, nmax=N_axes-1; ndof=0;
  for(i=0; i<natoms; i++){
    if((i >= iatom)&&(ndof < nmax)){
      ndof++; if(ndof<nmax)iatom=axes[ndof+1].first_atom;
    }
    if(axes[ndof].res==atoms[i].res)atoms[i].last_dof_side=ndof;
  }

  // Association between main chain degrees of freedom and bonds


  // Print degrees of freedom of side chains
  ires=-1;
  for(i=0; i<N_axes; i++){
    axe=axes+i;
    if(axe->res!=ires){
      ires=atoms[axe->atom1].res; numdof=0;
      printf("-%c%d:", atoms[axe->atom1].aa, ires+1);
    }
    numdof++;
    printf("%1d", numdof);
  }
  printf("\n");


  return(N_axes);
}

int Count_interactions(atom *atoms, int natoms,
		       struct interaction *Int_list, int N_int)
{
  int i;
  for(i=0; i<natoms; i++)atoms[i].N_int=0;
  for(i=0; i<N_int; i++){
    atoms[Int_list[i].i1].N_int++;
    atoms[Int_list[i].i2].N_int++;
  }
  return i;
}

struct bond *Find_bond_l(struct axe *axe, struct bond *bonds)
{
  // Bond length of bond associated with dof axe
  // Find the bond whose defining atom coincides with  axe->atom2
  int i_atom=axe->atom2;
  if(bonds[i_atom].i_atom==i_atom)return(bonds+i_atom);
  struct bond *bond=bonds;
  while(bond->i_atom>=0){
    if(bond->i_atom==i_atom)return(bond); bond++;
  }
  printf("ERROR, axe %d %d %c not found\n",
	 axe->atom1, axe->atom2, axe->type);
  exit(8);
  return(NULL);
}

struct bond *Find_bond_a(struct axe *axe, struct bond *bonds)
{
  // Bond angle of bond associated with dof axe
  // Find the min. bond whose previous defining atom coincides with axe->atom2
  int i_atom=axe->atom2;

  int i=i_atom+1; struct bond *bond=bonds+i;
  while(bond->i_atom>=0){
    if((bond->previous)&&((bond->previous)->i_atom==i_atom))
      return(bond);
    bond++;
  }
  printf("ERROR, axe %d %d %c not found\n",
	 axe->atom1, axe->atom2, axe->type);
  exit(8);
  return(NULL);
}

struct bond *Find_bond_t(struct axe *axe, struct bond *bonds, int *i_bond)
{
  // Torsion angle of bond associated with dof axe
  // Find all bonds whose previous defining atom coincides with axe->atom2
  int i_atom=axe->atom2;

  int i=*i_bond+1, res=axe->res+1; struct bond *bond=bonds+i;
  while(bond->i_atom>=0){
    if(bond->res > res)break;
    if((bond->previous)&&((bond->previous)->i_atom==i_atom)){
      *i_bond=i; return(bond);
    }
    bond++; i++;
  }
  if(*i_bond<0){
    printf("ERROR, axe %d %d %c not found\n",
	   axe->atom1, axe->atom2, axe->type);
    exit(8);
  }
  return(NULL);
}

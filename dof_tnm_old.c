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

char NONE='-';

struct side_chain{
  char aa;
  int res;
  int first_dof;
  int first_atom, first_ref;
  int last_atom, last_ref;
  int ndof, oldref;
  int numdof;
};



#define DEBUG 0

/*********************************************************************/
// Setting degrees of freedom
int Set_dof(char type, struct axe **axe, int *num_dof, struct bond *bond,
	    int *ndof, int *oldref, int *refatom, int natoms,
	    int last_atom, int last_ref);
char Axe_prot_main(struct bond *bond, int r_last);
char Axe_nuc_main(struct bond *bond, int r_last);
int Find_ref(int i_atom, int natoms, int *refatom);
void Axe_not_found(struct axe *axe);
void Axe_prot_side(struct axe *axe_side, int *N_side, struct side_chain *sc,
		   struct bond *bond, atom *atoms, int natoms, int *refatom);
int Set_side_chains(struct axe *axes, struct bond *bonds, atom *atoms,
		    int natoms, int *refatom, int *atom_ref,
		    int nmain, int MIN_INT);

int Count_interactions(atom *atoms, int natoms,
		       struct interaction *Int_list, int N_int);


float D2_GAP=6.5;
int Count_gaps(struct bond *bonds, int natoms);
int Select_axes_int(struct axe *axes, int *N_axes, int start,
		    atom *atoms, int natoms, int MIN_INT);
int Last_axe(int *aa, struct axe *axes, int ini, int n, int i);

/*********************************************************************/
// Auxiliary

float Distance_square(double *r1, double *r2);

/*********************************************************************/
/*********************************************************************/
// Start:

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
			   int MIN_INT, // Min. interactions per dof
			   int MIN_INT_SIDE, // Min. interactions per dof
			   int OMEGA,     // Are Omega angles stored?
			   int SCHAIN,    // Are side chains used?
			   int PSI)       // Are psi angles used?
{
  printf("\nSetting degrees of freedom for %d chains\n", Nchain);
  for(int i=0; i<Nchain; i++){
    printf("Chain %d type %d\n", i+1, chains[i].type);
  }

  // Allocate degrees of freedom
  int ndof_prot=3, ndof_dna=7, ndof_rb=6*(Nchain-1), i;
  int N_axe_max=ndof_rb; 
  if(OMEGA)ndof_prot=3;
  if(SCHAIN){ndof_prot+=4; ndof_dna+=1;}
  for(i=0; i<Nchain; i++){
    if((chains[i].type==1)||(chains[i].type==-1)){
      N_axe_max+=ndof_prot*chains[i].nres;
    }else if((chains[i].type==2)||(chains[i].type==-2)||
	     (chains[i].type==3)||(chains[i].type==-3)){
      N_axe_max+=ndof_dna*chains[i].nres;
    }
  }
  //int Ngap=Count_gaps(bonds, natoms);
  //N_axe_max+=3*Ngap;

  // Prepare reference atoms
  int refatom[natoms];
  for(i=0; i<natoms; i++)refatom[i]=-1;
  for(i=0; i<N_ref; i++)refatom[atom_ref[i]]=i;
  printf("%d reference atoms, last: %d\n", N_ref, atom_ref[N_ref-1]);
  /*printf("Reference atoms: ");
  for(i=0; i<natoms; i++)
    if(refatom[i]<0){printf("0");}else{printf("1");}
    printf("\n");*/

  // Prepare chains
  int ichain; struct chain *chain=chains;
  for(ichain=0; ichain<Nchain; ichain++){
    chain->ini_main=-1;
    chain->last_atom=chain->ini_atom+chain->natoms-1;
    printf("ichain= %d atoms=%d-%d\n",
	   ichain, chain->ini_atom, chain->last_atom);
    for(i=chain->last_atom; i>=chain->ini_atom; i--){
      if(refatom[i]>=0){chain->last_ref=refatom[i]; break;}
    }
    if(chain->last_ref==0){
      printf("ERROR, chain %d last ref %d %d last_atom %d\n",
	     ichain, i, refatom[i], chain->last_atom);
    }
    chain++;
  }

  // By default, degrees of freedom are frozen and not stored
  int ibond; struct bond *bond=bonds;
  for(ibond=0; ibond<natoms; ibond++){
    //bond->ibond=ibond;
    bond->len=-1; bond->angle=-1; bond->torsion=-1;
    bond->stored=0; bond++;
  }

  // Store axes
  /* A new degree of freedom is created only if there is at least one
     reference atom for 3 new degrees of freedom, and if the number
     of interactions is at least MIN_INT
  */
  struct axe *axes=malloc(N_axe_max*sizeof(struct axe)), *axe=axes;
  int N_axes=0, nrigid=0;
  struct axe axe_side[3*nres]; int N_side=0;
  struct side_chain sc; sc.res=-1;
  // Store main chain and side chain dofs for polypeptide and nuc. acids
  char type;  // Bond type: lat fpo 123456 s

  // Set first chain
  int achain=2; chain=chains; ichain=0;
  int last_r=chain->ini_res+chain->nres-1;
  int last_a=chain->last_atom;
  int last_ref=chain->last_ref, oldref=0;
  int ctype=chain->type, ndof=0;
  for(ibond=3; ibond<natoms; ibond++){
    bond=bonds+ibond;
    if(bond->atom==NULL)break;
    if(bond->atom->ali<0)continue;
    if(bond->previous==NULL)continue;
    if((bond->atom->chain!=bond->previous->atom->chain)||(ichain<0)){
      ichain++; chain=chains+ichain;
      last_r=chain->ini_res+chain->nres-1;
      last_a=chain->last_atom;
      last_ref=chain->last_ref;
      ctype=chain->type;
      achain=0;
    }
    if((achain==0)&&(ichain)){
      Set_dof('l', &axe, &N_axes, bond, &ndof,
	      &oldref, refatom, natoms, last_a, last_ref);
      Set_dof('a', &axe, &N_axes, bond, &ndof,
	      &oldref, refatom, natoms, last_a, last_ref);
      Set_dof('t', &axe, &N_axes, bond, &ndof,
	      &oldref, refatom, natoms, last_a, last_ref);
      (axe-3)->rigid=1; (axe-2)->rigid=1; (axe-1)->rigid=1; 
      nrigid+=3;
    }else if((achain==1)&&(ichain)){
      Set_dof('a', &axe, &N_axes, bond, &ndof,
	      &oldref, refatom, natoms, last_a, last_ref);
      Set_dof('t', &axe, &N_axes, bond, &ndof,
	      &oldref, refatom, natoms, last_a, last_ref);
      (axe-2)->rigid=1; (axe-1)->rigid=1; nrigid+=2;
    }else if((achain==2)&&(ichain)){
      Set_dof('t', &axe, &N_axes, bond, &ndof,
	      &oldref, refatom, natoms, last_a, last_ref);
      (axe-1)->rigid=1; nrigid+=1;
    }else{
      //if(bond->terminal)continue;
      if(strncmp(bond->atom->name,"CB", 2)==0)continue;
      if(bond->previous->stored)continue;
      if(ichain>=Nchain)break;
      if(ctype==0){
	continue;
      }else if(ctype==1){
	type=Axe_prot_main(bond->previous, last_r);
      }else{
	type=Axe_nuc_main(bond->previous, last_r);
      }
      if(type=='s'){
	if((SCHAIN)&&(bond->atom->res>0))
	  Axe_prot_side(axe_side, &N_side, &sc, bond->previous,
			atoms, natoms, refatom);
	continue;
      }
      if((Distance_square(bond->r, bond->previous->r) > D2_GAP)||
	 (atoms[bond->previous->i_atom].ali<0)){
	(*N_diso)++;
	Set_dof('l', &axe, &N_axes, bond, &ndof,
		&oldref, refatom, natoms, last_a, last_ref);
	Set_dof('a', &axe, &N_axes, bond, &ndof,
		&oldref, refatom, natoms, last_a, last_ref);
	Set_dof('t', &axe, &N_axes, bond, &ndof,
		&oldref, refatom, natoms, last_a, last_ref);
      }else if(type==NONE){
	continue;
      }else if(type=='o'){ // omega angle: set for transitions with D_omega>thr
	if((OMEGA>0)||((OMEGA<0)&&((bond+1)->d_phi>D_omega_thr))){
	  Set_dof(type, &axe, &N_axes, bond, &ndof,
		  &oldref, refatom, natoms, last_a, last_ref);
	  if(OMEGA<0){
	    atom *atom=(bond+1)->atom;
	    printf("Omega %s %c%d d=%d\n", atom->name,
		   atom->aa, atom->res, (bond+1)->d_phi);
	  }
	}
      }else{
	if((PSI==0)&&(type=='p'))continue;
	Set_dof(type, &axe, &N_axes, bond, &ndof,
		&oldref, refatom, natoms, last_a, last_ref);
      }
    }
    achain++; // Number of amino acid in the chain
  }
  printf("\n");
  ichain++;
  if(ichain != Nchain){
    printf("WARNING, expected %d aligned chains, found %d\n", Nchain, ichain);
  }
  printf(" %d rigid degrees of freedom found (expected %d for %d chains)\n",
	 nrigid, ndof_rb, Nchain);

  // Exclude dof with few interactions (except rigid body)
  *nskip=0; int nskip_side=0;
  if(MIN_INT>0){
    Count_interactions(atoms, natoms, Int_list, N_int);
    printf("Main axes skipped since < %d interactions:\n", MIN_INT);
    (*nskip)=
      Select_axes_int(axes, &N_axes, 0, atoms, natoms, MIN_INT);
    printf("\nSkipped %d axes over %d (main)\n", *nskip, N_axes+*nskip);
    if(N_side){
      printf("Side axes skipped since < %d interactions:\n", MIN_INT_SIDE);
      nskip_side=
	Select_axes_int(axe_side, &N_side, 0, atoms, natoms, MIN_INT_SIDE);
      printf("\nSkipped %d axes over %d (side)\n",
	     nskip_side, N_side+nskip_side);
    }
  }
  (*nmain)=N_axes;
  (*naxe)=N_axes;
  if(N_side){
    struct axe *axe1=axes+N_axes, *axe2=axe_side; int a;
    for(a=0; a<N_side; a++){*axe1=*axe2; axe1++; axe2++;}
    N_axes+=N_side; *naxe=N_axes; (*nskip)+=nskip_side;
  }

  // Print degrees of freedom
  int ig=1; ichain=atoms[0].chain;
  for(i=0; i<N_axes; i++){
    axe=axes+i;
    if(axe->bond->atom->chain!=ichain){
      if(axe->rigid)printf("-NewChain-");
      ichain=axe->bond->atom->chain; ig=0;
    }else if((axe->type=='l') && ig){
      printf("-Gap %2s%d-%2s%d-",
	     axe->bond->previous->atom->name,
	     axe->bond->previous->atom->res,
	     axe->bond->atom->name, axe->bond->atom->res);
      ig=0;
    }else{
      ig=1;
    }
    printf("%c%c", axe->type, chains[ichain].label);
    if(axe->type=='s')printf("%d",axe->bond->atom->res);
  }
  printf("\n");
  exit(8);

  // First and last degree of freedom for each atom
  atom *atom1=atoms;
  int am=0, as=*nmain; // main and side axe
  //printf("Last axe of atoms\n");
  for(i=0; i<natoms; i++){
    atom1->last_dof_main=Last_axe(&am, axes, 0, *nmain, i);
    if(N_side){
      atom1->last_dof_side=Last_axe(&as, axes, *nmain, N_axes, i);
    }else{
      atom1->last_dof_side=-1;
    }
    if(DEBUG){
      if(atom1->chain != (atom1-1)->chain)printf("\nChain %d ", atom1->chain);
      printf("%d ", atom1->last_dof_main);
      if(0){
	printf("%s %d axes: %d-%d %d  am: %d %d  as: %d %d %d\n",
	       atom1->name, atom1->res,
	       axes[atom1->last_dof_main].first_main, atom1->last_dof_main,
	       atom1->last_dof_side, am, axes[am].first_atom,
	       as, axes[as].first_atom, axes[as].bond->atom->res);
      }
    }
    atom1++;
  }
  if(DEBUG)printf("\n");

  // First degree of freedom of each chain and links between dofs
  for(ichain=0; ichain<Nchain; ichain++){
    chains[ichain].ini_main=-1; chains[ichain].mainaxes=-1;
  }
  printf("%d Axes: %d main, %d side\n", N_axes, *nmain, N_axes-*nmain);
  int a, first_main=-1; ichain=-1;
  printf("Sequence of dofs: ");
  for(a=0; a<*nmain; a++){
    axe=axes+a;
    axe->first_side=-1;
    if(axe->bond->atom->chain != ichain){
      ichain=axe->bond->atom->chain;
      chains[ichain].ini_main=a;
      axe->previous=-1;
      first_main=a;
      if(a)printf("%d", a-1);
      printf("\nChain %d: ", axe->bond->atom->chain+1);
      if(ichain){
	chains[ichain-1].mainaxes=a-chains[ichain-1].ini_main;
      }
    }else{
      axe->previous=a-1;
    }
    axe->first_main=first_main;
    axe->last_main=a;
    printf("%d ", axe->previous);
  }
  printf("%d\n", a-1);
  for(a=*nmain; a<N_axes; a++){
    axe=axes+a;
    if(axe->bond->atom->res == (axe-1)->bond->atom->res)axe->previous=a-1;
    int last=atoms[axe->first_atom].last_dof_main;
    if(last>=0){
      axe->last_main=last;
      axe->first_main=axes[last].first_main;
      if(axes[last].first_side<0){
	axes[last].first_side=a;
	axe->previous=last;
      }
    }
    printf("%d ", axe->previous);
  }
  if(N_axes>*nmain)printf("%d\n", a-1);

  ichain=Nchain-1;
  if(chains[ichain].ini_main<0)chains[ichain].ini_main=*nmain;
  chains[ichain].mainaxes=*nmain-chains[ichain].ini_main;
  for(i=0; i<Nchain; i++){
    struct chain *ch=chains+i;
    printf("chain %d (%c): %3d res %3d axes %d-%d\n",
	   i, ch->label, ch->nres, ch->mainaxes,
	   ch->ini_main, ch->ini_main+ch->mainaxes-1);
  }

  // Association between main chain degrees of freedom and bonds
  printf("Degrees of freedom associated to bonds:\n");
  for(i=0; i<N_axes; i++){
    axe=axes+i;
    printf("%c %s%d ", axe->type, axe->bond->atom->name, axe->bond->atom->res);
    if(axe->type=='l'){
       axe->bond->len=i;
    }else if(axe->type=='a'){
      axe->bond->angle=i;
    }else{
      axe->bond->torsion=i;
      // Find all bonds whose previous bond coincides with axe->axe
      struct bond *bond=axe->axe+1;
      int res_max=bond->atom->res+1;
      while((bond->atom)&&(bond->atom->res <= res_max)){
	if((bond!=axe->bond)&&(bond->previous==axe->axe)){
	  bond->torsion=i;
	  printf("%s%d ", bond->atom->name, bond->atom->res);
	}
	bond++;
      }
    }
  }
  printf("\n");

  if(0){
    axe=axes;
    printf("#Axes: type a1 a2 first_atom last_atom first_kin last_kin\n");
    for(i=0; i<N_axes; i++){
      printf("%c  ", axe->type);
      atom *a=axe->axe->previous->atom;
      printf(" %s%d", a->name, a->res);
      a=axe->axe->atom; printf(" %s%d", a->name, a->res);
      a=atoms+axe->first_atom; printf(" %s%d", a->name, a->res);
      a=atoms+axe->last_atom;  printf(" %s%d", a->name, a->res);
      a=atoms+axe->first_kin;  printf(" %s%d", a->name, a->res);
      a=atoms+axe->last_kin;   printf(" %s%d\n", a->name, a->res);
      axe++;
    }
  }

  // axes, n_axe for each chain
  //~ struct chain *chp=chains;
  //chains[0].ini_main=0;  chains[0].mainaxes=N_axes;
  // Change: already assigned

  /*int ichain=0, n=0, chain=res[0].chain;
    for(i=0; i<N_axes; i++)axes[i].chain=ichain;*/

  /*if(SCHAIN){
    (*naxe)+=Set_side_chains(axes+N_axes, bonds, atoms, natoms,
			     refatom, atom_ref, *naxe, MIN_INT_SIDE);
			     }*/
  return(axes);
}

int Set_dof(char type, struct axe **axe, int *num_dof, struct bond *bond,
	    int *ndof, int *oldref, int *refatom, int natoms,
	    int last_atom, int last_ref)
{
  /* Before updating the reference atom, maximum 3 degrees of freedom
     can be accepted   */
  int first_atom=bond->i_atom;
  int iref=Find_ref(first_atom, natoms, refatom);
  if((iref<0)||(last_ref<0)||(last_ref<iref)){
    atom *atom1=bond->previous->atom, *atom2=bond->atom;
    printf("\nWARNING, axe %s/%c%d/%d %s/%c%d/%d type %c ",
	   atom1->name, atom1->aa, atom1->res, atom1->chain,
	   atom2->name, atom2->aa, atom2->res, atom2->chain, type);
    printf("doesn't move any atom\n");
    printf("First: %d Last: %d first_ref: %d last_ref: %d\n",
	   first_atom, last_atom, iref, last_ref);
    printf("Refatom[first]= %d Refatom[last]= %d\n",
	   refatom[first_atom], refatom[last_atom]); 
    atom2=(bond+1)->atom;
    printf("Next bond: %s%d%c\n", atom2->name, atom2->res, atom2->aa);
    return(0);
  }

  if((*ndof >= 2)&&(iref <= *oldref)){
    printf("=(%d-%d)",iref, *oldref); return(0);
  }
  if(iref > *oldref){
    *oldref=iref; *ndof=0;
  }else{
    (*ndof)++;
  }

  //printf("%c%d", type, bond->atom->chain);
  printf("%c", type);
  (*axe)->type=type; //l,a,t=(p,f,o)(1,2,3,4,5,6)(s)
  (*axe)->bond=bond;
  // rotation axe associated to dof
  if((type=='l')||(type=='a')){(*axe)->axe=bond;}
  else{(*axe)->axe=bond->previous; (*axe)->axe->stored=1;} // Torsion
  (*axe)->first_atom=first_atom;
  (*axe)->last_atom=last_atom;
  (*axe)->first_kin=iref;
  (*axe)->last_kin=last_ref;
  (*axe)->first_side=-1;
  (*axe)->rigid=0;
  (*axe)++;
  (*num_dof)++;
  return(1);
}

int Count_gaps(struct bond *bonds, int natoms){
  int Ngap=0, ichain=-1, ires=-1;
  for(int i=0; i<natoms; i++){
    struct bond *bond=bonds+i;
    if(bond->atom->chain!=ichain){
      ichain=bond->atom->chain; ires=bond->atom->res;
      continue;
    }
    if(bond->atom->res!=ires){
      if((ires>=0)&&(bond->previous!=NULL)){
	if(Distance_square(bond->r, bond->previous->r) > D2_GAP)Ngap++;
      }
      ires=bond->atom->res;
    }
  }
  printf("%d gaps found in protein\n", Ngap);
  return(Ngap);
}


float Distance_square(double *r1, double *r2){
  double d0=r1[0]-r2[0], d1=r1[1]-r2[1], d2=r1[2]-r2[2];
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
    if(bond->atom==NULL)break;
    if(bond->previous==NULL)continue;
    if(bond->atom->res!=ires){
      // Residue start
      ires=bond->atom->res;
      aa=bond->atom->aa;
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
    }else if((strncmp(bond->atom->name,"CB", 2)==0)||
	     (strncmp(bond->atom->name,"CG", 2)==0)){
      goto test_ref;
    }else if((aa=='F')||(aa=='W')||
	     (aa=='Y')||(aa=='H')){
      continue;  // if aromatic ring, exit
    }else if((bond->terminal)||
	     (strncmp(bond->atom->name,"N ", 2)==0)||
	     (strncmp(bond->atom->name,"CA", 2)==0)||
	     (strncmp(bond->atom->name,"C ", 2)==0)||
	     (strncmp(bond->atom->name,"O ", 2)==0)){
      continue;
    }

  test_ref:
    move_atom=bond->i_atom+1;
    iref=refatom[move_atom];
    if((iref<0)||(iref>last_ref)||
       (atoms[atom_ref[iref]].res!=ires))
      continue;

    sprintf(type, "%1d", numdof);
    if(Set_dof(type[0], &axe, &N_axes, bond, &ndof,
	       &oldref, refatom, natoms, last_atom, last_ref)){
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
	if(atoms[i].res!=axe->bond->atom->res)break;
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
    if(axes[ndof].bond->atom->res==atoms[i].res)atoms[i].last_dof_side=ndof;
  }

  // Association between main chain degrees of freedom and bonds
  for(i=0; i<natoms; i++){
    bond=bonds+i;
    if(bond->atom==NULL)break;
    atom *atom=bond->atom;
    int n=atom->last_dof_side; if(n<0)continue;
    if((bond->previous!=NULL)&&(axes[n].axe==bond->previous)){
      if(axes[n].type=='t')bond->torsion=n;
    }
  }

  // Print degrees of freedom of side chains
  ires=-1;
  for(i=0; i<N_axes; i++){
    axe=axes+i;
    if(axe->bond->atom->res!=ires){
      ires=axe->bond->atom->res; numdof=0;
      printf("-%c%d:", axe->bond->atom->aa, ires+1);
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

void Axe_not_found(struct axe *axe){
  printf("ERROR, %c axe %s%d-%s%d not found\n", axe->type,
	 axe->bond->atom->name, axe->bond->atom->res,
	 axe->bond->previous->atom->name, axe->bond->previous->atom->res);
  exit(8);
}

int Find_ref(int i_atom, int natoms, int *refatom){
   int i=i_atom, iref=refatom[i]; 
   while(iref<0){
      i++; if(i>=natoms)break; iref=refatom[i]; 
    }
   return(iref);
 }

char Axe_prot_main(struct bond *bond, int r_last)
{
  char *name=bond->atom->name, *previous_name=bond->previous->atom->name;
  if(strncmp(name,"N ", 2)==0){
    if((strncmp(previous_name,"C ", 2)==0)&&
       (bond->atom->res==(bond->previous->atom->res+1))){
      return('o'); // omega 
    }else{return('q');}         // Pseudo-axe after disordered loop
  }else if(strncmp(name,"CA", 2)==0){
    if((strncmp(previous_name,"N ", 2)==0)&&
       (bond->atom->res==bond->previous->atom->res)){
      return('f'); // phi
    }else{return('q');}
  }else if(strncmp(name,"C ", 2)==0){
    if((strncmp(previous_name,"CA", 2)==0)&&
       (bond->atom->res==bond->previous->atom->res)){
      if(bond->atom->res==r_last)return(NONE);
      return('p'); // psi
    }else{return('q');}
  }else{
    return('s');
  }
  return(NONE);
}

char Axe_nuc_main(struct bond *bond, int r_last)
{
  char *name=bond->atom->name, *previous_name=bond->previous->atom->name;
  if(strncmp(name,"O5'", 3)==0){
    if(strncmp(previous_name,"P ", 3)==0)return('1');
  }else if(strncmp(name,"C5'", 3)==0){
    if(strncmp(previous_name,"O5'", 3)==0)return('2');
  }else if(strncmp(name,"C4'", 3)==0){
    if(strncmp(previous_name,"C5'", 3)==0)return('3');
  }else if(strncmp(name,"C3'", 3)==0){
    if(strncmp(previous_name,"C4'", 3)==0)return('4');
  }else if(strncmp(name,"O3'", 3)==0){
    if(strncmp(previous_name,"C3'", 3)==0)return('5');
  }else if(strncmp(name,"P  ", 3)==0){
    if(strncmp(previous_name,"O3'", 3)==0)return('6');
    /*}else if(strncmp(name,"C1'", 3)==0){
      if(strncmp(previous_name,"N", 1)==0)return('s');*/ // side chain!
  }else{
    return(NONE);
  }
  return(NONE);
  //return('q'); // Disordered main axis
}

void Axe_prot_side(struct axe *axe_side, int *N_side, struct side_chain *sc,
		   struct bond *bond, atom *atoms, int natoms, int *refatom)
{
  if(bond->atom->res!=sc->res){
    // Residue start
    sc->res=bond->atom->res;
    sc->aa= bond->atom->aa;
    if((sc->aa=='G')||(sc->aa=='A')||(sc->aa=='P'))return;
    
    // Find first and last atom in side chain
    sc->first_atom=-1;
    sc->last_atom=-1;
    sc->first_ref=-1;
    sc->last_ref=-1;
    int i;
    for(i=bond->i_atom; i<natoms; i++){
      atom *atom_i=atoms+i;
      if((atom_i->res>sc->res)||
	 (strncmp(atom_i->name,"OT", 2)==0)||
	 (strncmp(atom_i->name,"OX", 2)==0))
	break;
      if((strncmp(atom_i->name,"N ", 2)==0)||
	 (strncmp(atom_i->name,"CA", 2)==0)||
	 (strncmp(atom_i->name,"CB", 2)==0))continue;
      if((strncmp(atom_i->name,"C ", 2)!=0)&&
	 (strncmp(atom_i->name,"O ", 2)!=0)){
	sc->last_atom=i;
	if(sc->first_atom<0)sc->first_atom=i;
	int iref=refatom[i];
	if(iref>=0){
	  sc->last_ref=iref; if(sc->first_ref<0)sc->first_ref=iref;
	}
      }
    }
    sc->ndof=0; sc->numdof=0; sc->oldref=0;
  }
  
  if((bond->terminal)||
     (strncmp(bond->atom->name,"N ", 2)==0)||
     (strncmp(bond->atom->name,"CA", 2)==0)||
     (strncmp(bond->atom->name,"C ", 2)==0)||
     (strncmp(bond->atom->name,"O ", 2)==0))
    return;
    
  if((sc->first_atom<0)||(sc->first_ref<0)||
     (sc->aa=='G')||(sc->aa=='A')||(sc->aa=='P'))return;
  
  if(((strncmp(bond->atom->name,"CB", 2)!=0)&&
      (strncmp(bond->atom->name,"CG", 2)==0))&&
     ((sc->aa=='F')||(sc->aa=='W')||(sc->aa=='Y')||(sc->aa=='H'))){
    return;  // if aromatic ring, exit
  }

  //sprintf(type, "%1d", sc->numdof);
  struct axe *axe=axe_side+*N_side;
  if(Set_dof('s', &axe, N_side, bond, &sc->ndof,
	     &sc->oldref, refatom, natoms, sc->last_atom, sc->last_ref)){
    sc->numdof++;
  }
  return;
}


int Select_axes_int(struct axe *axes, int *N_axes, int start,
		    atom *atoms, int natoms, int MIN_INT)
{
  int n_skipped=0;
  int n=start, last=*N_axes-1, last_a, a, i;
  for(a=start; a<*N_axes; a++){
    struct axe *axe=axes+a;
    if(a<last){
      struct axe *next=axe+1;
      if(((next->type=='l')||(axe->type=='a'))&&((a+1)<last)){
	next++;
      }
      last_a=next->first_atom;
    }else{
      last_a=natoms;
    }
    int num_int=0;
    for(i=axe->first_atom; i< last_a; i++)num_int+=atoms[i].N_int;
    if((num_int >= MIN_INT)||(axe->rigid)){
      if(a!=n)axes[n]=*axe;
      n++;
    }else{
      printf("%d ", a);
      n_skipped++;
    }
  }
  *N_axes=n;
  return(n_skipped);
}

int Last_axe(int *aa, struct axe *axes, int ini, int n, int i)
{
  int a=*aa; struct axe *axe=axes+(*aa);
  while(a<n){
    if(axe->first_atom>i)break; a++; axe++;
  }
  if(a==n){a--; axe--;}
  if(a && axe->first_atom>i && (axe-1)->first_atom<=i){a--; axe--;}
  if(a<ini){*aa=ini;}else{*aa=a;}
  if((axe->first_atom<=i))return(a); //&&(axe->last_atom>=i)
  if(i>10){
    printf("WARNING, last_axe<0: atom %d axe %d n= %d first %d last %d\n",
	   i, a, n, axe->first_atom, axe->last_atom); //exit(8);
  }
  return(-1);
}

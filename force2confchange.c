#include "coord.h"
#include "buildup.h"
#include "force2confchange.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "read.h"

static atom *Find_atom_seq(struct residue *seq, int nres, atom atom1);
int Read_force(float *Force_input_cart, atom *atoms1,
		      int *atom_num1, int N_ref, char *FILE_FORCE);
static int Find_atom_seq1(int *iref, char *name, char *aa1, int res,
			  atom *atoms1, int *atom_num1, int N_ref);

int Force2Confchange(float *atom_move, char *FILE_FORCE,
		     int naxe1, int N_ref, atom *atoms1, int natoms1,
		     int *atom_num1, float **Jacobian_ar, float *eigen_value,
		     float **Tors_mode, int *select, struct axe *axe1,
		     char *chain1, int nres1, struct residue *seq1,
		     struct bond *bonds)
{
  int N_cart=3*N_ref, i, ia, j, k;
  double Delta_theta[naxe1];
  float Force_input_cart[N_cart];
  float Force_input_tors[naxe1];
  float Force_input_coef[naxe1];
  double sum;

  //i=Read_force(Force_input_cart, atoms1, atom_num1, N_ref, FILE_FORCE);

  // Read force
  int nres, nmr, ANISOU;
  //~ int natoms;
  atom atom_read[natoms1], *atom;
  struct residue seq[nres1];
  printf("Reading force in file %s\n", FILE_FORCE);
  nres=Read_coord(FILE_FORCE, &nmr, seq, atom_read, chain1, &ANISOU, NULL);
  if(nres<0)return(-1); // file not found
  if(nres!=nres1){
    printf("ERROR, wrong number of residues in force file:");
    printf(" %d instead of %d\n", nres, nres1); return(-1);
  }
  for(i=0; i<nres; i++){
    if(seq[i].amm!=seq1[i].amm){
      printf("ERROR, wrong amino acid in force file\n");
      printf("Position %d, %c instead of %c\n", i, seq[i].amm, seq1[i].amm);
      return(-1);
    }
  }
  for(i=0; i<3*N_ref; i++)Force_input_cart[i]=0;
  for(i=0; i<N_ref; i++){
    atom=Find_atom_seq(seq, nres, atoms1[atom_num1[i]]);
    if(atom==NULL)continue;
    k=3*i; for(j=0; j<3; j++){Force_input_cart[k]=atom->r[j]; k++;}
  }
  // End read force

  for(i=0; i<naxe1; i++){  // Torsional force F^tors_a = sum_i F_i*J_ia
    sum=0; for(j=0; j<N_cart; j++)sum+=Jacobian_ar[i][j]*Force_input_cart[j];
    Force_input_tors[i]=sum;
  }
  for(ia=0; ia<naxe1; ia++){  // c^alpha =(F^tors,u^alpha)
    if(select[ia]==0)continue;
    float *ua=Tors_mode[ia];
    sum=0; for(i=0; i<naxe1; i++)sum+=ua[i]*Force_input_tors[i];
    Force_input_coef[ia]=sum;
  }
  for(i=0; i<naxe1; i++)Delta_theta[i]=0;
  for(ia=0; ia<naxe1; ia++){  // Delta_theta= H^{-1} F^tors
    if(select[ia]==0)continue;
    float *ua=Tors_mode[ia];
    float w=Force_input_coef[ia]/eigen_value[ia];
    for(i=0; i<naxe1; i++)Delta_theta[i]+=w*ua[i];
  }

  // Initialize build up
  Set_bonds_measure(bonds, natoms1, atoms1);
  // Build up Cartesian coordinates from Delta_theta
  Build_up(bonds, natoms1, Delta_theta, naxe1);
  k=3;
  for(i=1;i<N_ref;i++){
    struct bond *bond=bonds+atom_num1[i];
    for(j=0; j<3; j++){
      atom_move[k]=bond->r[j]; k++;
    }
  }
  return(0);
}

int Read_force(float *Force_input_cart, atom *atoms1,
	       int *atom_num1, int N_ref, char *FILE_FORCE)
{
  FILE *file_in=fopen(FILE_FORCE, "r");
  int iref=0, res; float f1, f2, f3;
  char name[20], aa1[5], string[1000];
  if(file_in==NULL){
    printf("ERROR, I could not find file with input force %s\n",
	   FILE_FORCE); return(-1);
  }
  while(fgets(string, sizeof(string), file_in)!=NULL){
    sscanf(string, "%s%s%d%f%f%f", name, aa1, &res, &f1, &f2, &f3);
    if(Find_atom_seq1(&iref, name, aa1, res, atoms1,atom_num1, N_ref)<0){
      printf("ERROR, reference number not found for %s", string);
      return(-1);
    }
    Force_input_cart[3*iref]=f1;
    Force_input_cart[3*iref+1]=f2;
    Force_input_cart[3*iref+2]=f3;
  }
  return(0);
}

atom *
Find_atom_seq(struct residue *seq, int nres, atom atom1)
{
  struct residue *seq_i=seq+atom1.res;
  if(seq_i->amm != atom1.aa){
    printf("ERROR, wrong residue type\n"); return(NULL);
  }
  atom *atom=seq_i->atom_ptr; int i;
  for(i=0; i<seq_i->n_atom; i++){
    if(strncmp(atom->name, atom1.name, 3)==0)return(atom);
    atom++;
  }
  return(NULL);
}


int Find_atom_seq1(int *iref, char *name, char *aa1, int res, atom *atoms1,
		   int *atom_num1, int N_ref)
{
  if(1)return(0);
  return(-1);
}

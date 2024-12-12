#include "nma_para.h"
#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "interactions_tnm.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MINT 1000 // Maximum number of screened interactions
                  // for any residue pair

/*******************************  Go model  *********************************/
int Interactions(struct interaction *Int_list, char *atom_type, float thr,
		 atom *atoms, int N_atoms, int N_res);
int Interactions_min(struct interaction *Int_list, float thr,
		     atom *atoms, int N_atoms, int N_res);
int Interactions_min_M(struct interaction *Int_list, float thr,
		       atom *atoms, int N_atoms, int N_res);
int Interactions_all(struct interaction *Int_list, float thr,
		     atom *atoms, int N_atoms, int N_res);
int Interactions_hnm(struct interaction *Int_list, char *atom_type, float thr,
		     atom *atoms, int N_atoms, int N_res);
int Interactions_HB(struct interaction *Int_list, float thr,
		    atom *atoms, int N_atoms, int N_res,
		    float thr_HB, float cos_HB, float ene_HB);
int Interactions_all_ref(struct interaction *Int_list, float thr,
			 atom *atoms, int N_atoms, int N_res, char *REF);

void Compute_sec_der(struct interaction *Int_list, int N_int, float expo);
void Purge_interactions_atom_res(int *N_int, struct interaction *Int_list,
				 int nres);
void Purge_interactions_res_pair(int *N_int, struct interaction *Int_list,
				 int nres);
void Select_pairs(int *selected, struct interaction **res12, int n12);

int Find_atom(atom *atoms, int *i1, int res, int N_atoms, char *atom_type);
float Dist2(float *r1, float *r2);
static float Bond_angle(float *r1, float *r2, float *r3);
int Check_store_HB(struct interaction *Int_list, int *N_int,
		   atom *atoms, int i_AA, int i_A, int i_D, int rA,
		   float d2);

/*****************************************************************************/
void Compute_interactions(int *N_int, struct interaction **Int_list,
			  char *INT_TYPE, atom *atoms, int natoms, int nres)
{
  INT_MAX=50*natoms;
  struct interaction *Int_tmp=malloc(INT_MAX*sizeof(struct interaction));

  if(HNM){
    *N_int=Interactions_hnm(Int_tmp, "CA", 20, atoms, natoms, nres);
  }else if((strncmp(INT_TYPE, "CA", 2)==0)||(strncmp(INT_TYPE, "CB", 2)==0)){
    *N_int=Interactions(Int_tmp, INT_TYPE, C_THR, atoms, natoms, nres);
  }else if((strncmp(INT_TYPE, "SCA", 3)==0)||
	   (strncmp(INT_TYPE, "SCB", 3)==0)){
    *N_int=Interactions_all_ref(Int_tmp, C_THR, atoms, natoms, nres, REF);
  }else if(strncmp(INT_TYPE, "MIN", 3)==0){
    if(N_RESRES==1){
      *N_int=Interactions_min(Int_tmp, C_THR, atoms, natoms, nres);
    }else{
      *N_int=Interactions_min_M(Int_tmp, C_THR, atoms, natoms, nres);
    }
  }else if(strncmp(INT_TYPE, "MIN", 3)==0){
    *N_int=Interactions_all(Int_tmp, C_THR, atoms, natoms, nres);
  }else if(strncmp(INT_TYPE, "SCR", 3)==0){
    *N_int=Screened_Interactions(Int_tmp,C_THR,atoms,natoms,3,S_THR,NOCOV);
  }else if(strncmp(INT_TYPE, "SHA", 3)==0){
    *N_int=Screened_Interactions(Int_tmp,C_THR,atoms,natoms,
				 S_TYPE,S_THR,NOCOV);
  }else if(strncmp(INT_TYPE, "HB", 2)==0){
    *N_int=Interactions_HB(Int_tmp, C_THR, atoms, natoms, nres,
			   DA_THR, COS_DAA, ENE_HB);
  }else{
    printf("WARNING, incorrect interaction type %s\n", INT_TYPE);
    printf("Setting to default MIN\n"); strcpy(INT_TYPE, "MIN");
    *N_int=Interactions_min(Int_tmp, C_THR, atoms, natoms, nres);
  }

  int k;
  for(k=0; k<*N_int; k++){
    Int_tmp[k].res1=atoms[Int_tmp[k].i1].res;
    Int_tmp[k].res2=atoms[Int_tmp[k].i2].res;
  }


  // Each atoms only interacts once per residue
  if((ONEINT)&&
     ((strncmp(INT_TYPE, "SCR", 3)==0)||(strncmp(INT_TYPE, "SHA", 3)==0))){
    //Purge_interactions_atom_res(N_int, Int_list, nres);
    Purge_interactions_res_pair(N_int, Int_tmp, nres);
  }

  if(EXP_HESSIAN>0){
    Compute_sec_der(Int_tmp, *N_int, EXP_HESSIAN);
    // Normalize
    double t=0;
    for(k=0; k<*N_int; k++)t+=Int_tmp[k].sec_der; t/=(*N_int);
    for(k=0; k<*N_int; k++)Int_tmp[k].sec_der/=t;
  }

  *Int_list=malloc((*N_int)*sizeof(struct interaction));
  for(k=0; k<*N_int; k++)(*Int_list)[k]=Int_tmp[k];
  free(Int_tmp);

}

/*******************************  Go model  *********************************/
int Interactions(struct interaction *Int_list, char *atom_type, float thr,
		 atom *atoms, int N_atoms, int N_res)
{
  int N_int=0, i1=0, i2;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr;
  //~ float t_max2=t_max*t_max;
  float d, d2;
  int res1, res2;

  for(res1=0; res1<N_res; res1++){

    if(Find_atom(atoms, &i1, res1, N_atoms, atom_type)< 0)continue;
    atom1=atoms+i1; i2=i1+1;

    for(res2=res1+1; res2<N_res; res2++){
      if(Find_atom(atoms, &i2, res2, N_atoms, atom_type)< 0)continue;
      atom2=atoms+i2;

      // Interatomic distance
      d=atom2->r[0]-atom1->r[0]; d2 =d*d; if(d2>thr2)continue;
      d=atom2->r[1]-atom1->r[1]; d2+=d*d; if(d2>thr2)continue;
      d=atom2->r[2]-atom1->r[2]; d2+=d*d; if(d2>thr2)continue;

      // store interaction
      Int_list[N_int].i1=i1; Int_list[N_int].i2=i2;
      Int_list[N_int].r2=d2; Int_list[N_int].sec_der=1;
      Int_list[N_int].rmin=d2;
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    }
  }
  return(N_int);
}


int Interactions_hnm(struct interaction *Int_list, char *atom_type, float thr,
		     atom *atoms, int N_atoms, int N_res)
{
  int N_int=0, i1=0, i2;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr;
  float d, d2;
  int res1, res2;

  for(res1=0; res1<N_res; res1++){

    if(Find_atom(atoms, &i1, res1, N_atoms, atom_type)< 0)continue;
    atom1=atoms+i1; i2=i1+1;

    for(res2=res1+1; res2<N_res; res2++){
      if(Find_atom(atoms, &i2, res2, N_atoms, atom_type)< 0)continue;
      atom2=atoms+i2;

      // Interatomic distance
      d=atom2->r[0]-atom1->r[0]; d2 =d*d; if(d2>thr2)continue;
      d=atom2->r[1]-atom1->r[1]; d2+=d*d; if(d2>thr2)continue;
      d=atom2->r[2]-atom1->r[2]; d2+=d*d; if(d2>thr2)continue;

      // store interaction
      Int_list[N_int].i1=i1; Int_list[N_int].i2=i2;
      Int_list[N_int].r2=d2; Int_list[N_int].rmin=d2;
      if(d2<=4){
	Int_list[N_int].sec_der=205*sqrt(d2)-571.2;
      }else{
	Int_list[N_int].sec_der=305.9*pow(d2, -3);
      }
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    }
  }
  return(N_int);
}

int Find_atom(atom *atoms, int *i1, int res, int N_atoms, char *atom_type){
  int i; atom *atom;

  atom=atoms+*i1;
  for(i=*i1; i<N_atoms; i++){
    if(atom->res > res)break;
    if(atom->res < res)goto next1;
    if(strncmp(atom->name, atom_type, 2)==0){*i1=i; return(0);}
  next1: atom++;
  }
  // Not found? Look for CA atom
  atom=atoms+*i1;
  for(i=*i1; i<N_atoms; i++){
    if(atom->res > res)break;
    if(atom->res < res)goto next2;
    if(strncmp(atom->name, "CA", 2)==0){*i1=i; return(0);}
  next2: atom++;
  }
  return(-1);
}

int Interactions_min(struct interaction *Int_list, float thr,
		     atom *atoms, int N_atoms, int N_res)
{
  // All pairs of atoms tested, only the pair at minimal distance interact
  int N_int=0, i1, i2, i1_start=0, i2_start=0;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr, t_max=4*thr;
  float t_max2=t_max*t_max;
  float d, d2, d2_min;
  int res1, res2;
  int i1_min=-1, i2_min=-1;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++; i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      while(atoms[i2_start].res< res2)i2_start++;
      d2_min=10000;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;

	  // Check if covalent
	  if((NOCOV) && Covalent(res1, res2, atom1, atom2))continue;

	  // Interatomic distance
	  d=atom2->r[0]-atom1->r[0];
	  d2 =d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[1]-atom1->r[1];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[2]-atom1->r[2];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  if(d2 < d2_min){d2_min=d2; i1_min=i1; i2_min=i2;}
	}
      }

      // store interaction
      if(d2_min < thr2){
	Int_list[N_int].i1=i1_min; Int_list[N_int].i2=i2_min;
	Int_list[N_int].r2=d2_min; Int_list[N_int].sec_der=1;
	Int_list[N_int].rmin=d2_min;
	N_int++;
	if(N_int >= INT_MAX){
	  printf("ERROR, too many interactions %d\n", N_int); exit(8);
	}
      }
    next2: continue;
    }
  }
  return(N_int);
}

int Interactions_min_M(struct interaction *Int_list, float thr,
		     atom *atoms, int N_atoms, int N_res)
{
  //  For each pair of residues only N_RESRES atoms
  //  with smallest distance interact
  int N_int=0, i1, i2, i1_start=0, i2_start=0;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr, t_max=4*thr;
  float t_max2=t_max*t_max;
  float d, d2;
  int res1, res2, k;

  float *d2_min=malloc(N_RESRES*sizeof(float));
  int *i1_min=malloc(N_RESRES*sizeof(int));
  int *i2_min=malloc(N_RESRES*sizeof(int));

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++; i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      while(atoms[i2_start].res< res2)i2_start++;
      int n_int=0;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;

	  // Check if covalent
	  if((NOCOV) && Covalent(res1, res2, atom1, atom2))continue;

	  // Interatomic distance
	  d=atom2->r[0]-atom1->r[0];
	  d2 =d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[1]-atom1->r[1];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[2]-atom1->r[2];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  if((n_int<N_RESRES)||(d2 < d2_min[n_int-1])){
	    int j, n;
	    for(k=0; k<n_int; k++)if(d2<d2_min[k])break;
	    if(n_int<N_RESRES){n=n_int-1; n_int++;}else{n=N_RESRES-1;}
	    for(j=n; j>k; j--){
	      d2_min[j]=d2_min[j-1];
	      i1_min[j]=i1_min[j-1];
	      i2_min[j]=i2_min[j-1];
	    }
	    d2_min[k]=d2; i1_min[k]=i1; i2_min[k]=i2;
	  }
	}
      }

      // store interaction
      for(k=0; k<n_int; k++){
	Int_list[N_int].i1=i1_min[k]; Int_list[N_int].i2=i2_min[k];
	Int_list[N_int].r2=d2_min[k]; Int_list[N_int].sec_der=1;
	Int_list[N_int].rmin=d2_min[k];
	N_int++;
	if(N_int >= INT_MAX){
	  printf("ERROR, too many interactions %d\n", N_int); exit(8);
	}
      }
    next2: continue;
    }
  }
  free(d2_min); free(i1_min); free(i2_min);
  return(N_int);
}


int Covalent(int res1, int res2, atom *atom1, atom *atom2){
  if(res2==(res1+1)){
    if((strncmp(atom2->name, "N ", 2)==0)&&
       ((strncmp(atom1->name, "CA", 2)==0)||
	(strncmp(atom1->name, "C ", 2)==0)||
	(strncmp(atom1->name, "O ", 2)==0)))return(1);
    if((strncmp(atom2->name, "CA", 2)==0)&&
       (strncmp(atom1->name, "C ", 2)==0)
       //||(strncmp(atom1->name, "O ", 2)==0)
       )return(1);
  }
  return(0);
}

int Bond(int res3, int res1, float d2_13, int res2, float d2_23)
{
  if((abs(res3-res2)<=1)&&(d2_23<5))return(1);
  if((abs(res3-res1)<=1)&&(d2_13<5))return(1);
  return(0);
}

int Interactions_all(struct interaction *Int_list, float thr,
		     atom *atoms, int N_atoms, int N_res)
{
  // All pairs of atoms below the thresholds interact
  int N_int=0, i1, i2, i1_start=0, i2_start=0;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr, t_max=4*thr;
  float t_max2=t_max*t_max;
  float d, d2;
  int res1, res2;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++; i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      while(atoms[i2_start].res< res2)i2_start++;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;

	  // Check if covalent
	  if((NOCOV) && Covalent(res1, res2, atom1, atom2))continue;

	  // Interatomic distance
	  d=atom2->r[0]-atom1->r[0];
	  d2 =d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[1]-atom1->r[1];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[2]-atom1->r[2];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  if(d2 < thr2){
	    Int_list[N_int].i1=i1; Int_list[N_int].i2=i2;
	    Int_list[N_int].r2=d2; Int_list[N_int].sec_der=1;
	    Int_list[N_int].rmin=d2;
	    N_int++;
	    if(N_int >= INT_MAX){
	      printf("ERROR, too many interactions %d\n", N_int);
	      exit(8);
	    }
	  } // End store interaction
	}
      }
    next2: continue;
    }
  }
  return(N_int);
}


int Interactions_all_ref(struct interaction *Int_list, float thr,
			atom *atoms, int N_atoms, int N_res, char *REF)
{
  int N_int=0, i1, i2, i1_start=0, i2_start=0;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr, t_max=4*thr;
  float t_max2=t_max*t_max;
  float d, d2, d2_min;
  int res1, res2;
  int CB1=0, CB2=0; //i1_min, i2_min, 
  int *notfound=malloc(N_res*sizeof(int));

  for(res1=0; res1<N_res; res1++)notfound[res1]=0;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++; CB1=i1_start;
    if(Find_atom(atoms, &CB1, res1, N_atoms, REF)<0){
      if(notfound[res1]==0)
	printf("WARNING, %s atom not found, res= %d\n", REF, res1);
      continue;
    }

    i2_start=i1_start+1;
    for(res2=res1+1; res2<N_res; res2++){
      while(atoms[i2_start].res< res2)i2_start++;  CB2=i2_start;
      if(Find_atom(atoms, &CB2, res2, N_atoms, REF)<0){
	if(notfound[res2]==0)
	  printf("WARNING, %s atom not found, res= %d\n", REF, res2);
	notfound[res2]=1;
	continue;
      }
      d2_min=10000;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;

	  // Check if covalent
	  if((NOCOV) && Covalent(res1, res2, atom1, atom2))continue;

	  // Interatomic distance
	  d=atom2->r[0]-atom1->r[0];
	  d2 =d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[1]-atom1->r[1];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[2]-atom1->r[2];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  if(d2 < d2_min)d2_min=d2; //i1_min=i1; i2_min=i2;}
	}
      }

      // store interaction
      if(d2_min < thr2){
	// Interacting atoms set as CB (CA if Gly)
	Int_list[N_int].i1=CB1; Int_list[N_int].i2=CB2;
	Int_list[N_int].r2=Dist2(atoms[CB1].r, atoms[CB2].r);
	Int_list[N_int].sec_der=1; Int_list[N_int].rmin=d2_min;
	N_int++;
	if(N_int >= INT_MAX){
	  printf("ERROR, too many interactions %d\n", N_int); exit(8);
	}
      }
    next2: continue;
    }
  }
  free(notfound);
  return(N_int);
}

int Interactions_HB(struct interaction *Int_list, float thr,
		    atom *atoms, int N_atoms, int N_res,
		    float thr_NO, float cos_HB, float ene_HB)
{
  // Hydrogen bonds, in the form of N-O pairs, are separated
  int N_int=0, i1, i2, i1_start=0, i2_start=0;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr, t_max=4*thr;
  float thr2_NO=thr_NO*thr_NO;
  float t_max2=t_max*t_max;
  float d, d2, d2_min;
  int res1, res2;
  int i1_min=-1, i2_min=-1, i_hb;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++; i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      while(atoms[i2_start].res< res2)i2_start++;
      d2_min=10000; i_hb=0;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	if((i_hb)&&(atom1->name[0]!='N')&&(atom1->name[0]!='O'))continue;
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;
	  if((i_hb)&&(atom2->name[0]!='N')&&(atom2->name[0]!='O'))continue;

	  // Check if covalent
	  if((NOCOV) && Covalent(res1, res2, atom1, atom2))continue;

	  // Interatomic distance
	  d=atom2->r[0]-atom1->r[0];
	  d2 =d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[1]-atom1->r[1];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[2]-atom1->r[2];
	  d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	  if(d2 < d2_min){d2_min=d2; i1_min=i1; i2_min=i2;}
	  if(d2 < thr2_NO){
	    if((atom1->name[0]=='N')&&(atom2->name[0]=='O')){
	      int i_AA=i2-1;
	      i_hb+=Check_store_HB(Int_list,&N_int,atoms,i_AA,i2,i1,1,d2);
	      if(i_hb==2)goto next2;
	    }else if((atom1->name[0]=='O')&&(atom2->name[0]=='N')){
	      int i_AA=i1-1;
	      i_hb+=Check_store_HB(Int_list,&N_int,atoms,i_AA,i1,i2,2,d2);
	      if(i_hb==2)goto next2;
	    }
	  }
	}
      }

      // store interaction
      if((i_hb)||(d2_min >= thr2))goto next2;
      Int_list[N_int].type=1;
      Int_list[N_int].i1=i1_min; Int_list[N_int].i2=i2_min;
      Int_list[N_int].r2=d2_min; Int_list[N_int].sec_der=1;
      Int_list[N_int].rmin=d2_min;
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    next2: continue;
    }
  }
  return(N_int);
}


void Compute_sec_der(struct interaction *Int_list, int N_int, float expo)
{
  int n;
  for(n=0; n<N_int; n++){
    if(Int_list[n].rmin<=0){
      printf("ERROR, d2= %.2f\n", Int_list[n].rmin); exit(8);
    }
    Int_list[n].sec_der*=pow(Int_list[n].rmin, -expo);
  }
}

float Dist2(float *r1, float *r2){
  return(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]);
}

int Check_store_HB(struct interaction *Int_list, int *N_int,
		   atom *atoms, int i_AA, int i_A, int i_D, int rA,
		   float d2)
{
  float cos_tau=Bond_angle(atoms[i_AA].r, atoms[i_A].r, atoms[i_D].r);
  if(cos_tau > COS_DAA)return(0);
  float d_DAA=Dist2(atoms[i_AA].r, atoms[i_D].r);
  Int_list[*N_int].type=1;   Int_list[*N_int].sec_der=1;
  Int_list[*N_int].r2=d_DAA; Int_list[*N_int].rmin=d_DAA;
  if(rA==1){
    Int_list[*N_int].i1=i_AA; Int_list[*N_int].i2=i_D;
  }else{
    Int_list[*N_int].i1=i_D;  Int_list[*N_int].i2=i_AA;
  }
  (*N_int)++;
  Int_list[*N_int].type=2; Int_list[*N_int].sec_der=ENE_HB;
  if(rA==1){
    Int_list[*N_int].i1=i_A; Int_list[*N_int].i2=i_D;
  }else{
    Int_list[*N_int].i1=i_D; Int_list[*N_int].i2=i_A;
  }
  Int_list[*N_int].r2=d2; Int_list[*N_int].rmin=d2;
  (*N_int)++;
  return(1);
}

float Bond_angle(float *r1, float *r2, float *r3){
  float ux, uy, uz, vx, vy, vz;
  float scal_prod;

  if((r1==NULL)||(r2==NULL)||(r3==NULL))return(0);
  ux=(r1[0]-r2[0]); vx=(r3[0]-r2[0]);
  uy=(r1[1]-r2[1]); vy=(r3[1]-r2[1]);
  uz=(r1[2]-r2[2]); vz=(r3[2]-r2[2]);
  scal_prod=ux*vx+uy*vy+uz*vz;
  scal_prod/=sqrt((ux*ux+uy*uy+uz*uz)*(vx*vx+vy*vy+vz*vz));
  return(scal_prod);
}

void Print_contact_matrix(struct interaction *Int_list, int N_int,
			  atom *atoms, struct residue *res,
			  char *name, char *pdb)
{
  char nameout[100]; sprintf(nameout, "%s_%s_Cont_Mat.txt", pdb,name);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing contact matrix in %s\n", nameout);
  int i, SCR=0, ALL=1;
  if((strncmp(name,"SCR",3)==0))SCR=1;
  if((strncmp(name,"CA",2)==0)||(strncmp(name,"CB",2)==0)||
     (strncmp(name,"SCA",3)==0)||(strncmp(name,"SCB",3)==0))ALL=0;
  fprintf(file_out, "#res1\tres2\td");
  if(SCR)fprintf(file_out, "\tw");
  if(ALL)fprintf(file_out, "\tatom1\tatom2");
  fprintf(file_out, "\tcontacts=%d\n",N_int);
  for(i=0; i<N_int; i++){
    fprintf(file_out, "%s\t%s\t%.2f",
	    res[atoms[Int_list[i].i1].res].pdbres,
	    res[atoms[Int_list[i].i2].res].pdbres,
	    sqrt(Int_list[i].r2));
    if(SCR)fprintf(file_out, "\t%.4f", Int_list[i].sec_der);
    if(ALL)fprintf(file_out, "\t%s\t%s",
		   atoms[Int_list[i].i1].name,
		   atoms[Int_list[i].i2].name);
    fprintf(file_out, "\n");
  }
  fclose(file_out);

}

void Purge_interactions_atom_res(int *N_int, struct interaction *Int_list,
				 int nres)
{
  int Nnew=0, res1, res2, k; 
  int selected[MINT];
  struct interaction *res12[MINT];
  struct interaction *ptr, *ini=Int_list, *end=Int_list+*N_int;
  struct interaction *tmp=malloc(*N_int*sizeof(struct interaction));

  for(res1=0; res1<nres; res1++){
    for(res2=res1+1; res2<nres; res2++){
      // Make a list of res1-res2 contacts.
      // Find the next value of res2
      int n12=0, res2min=nres;
      for(ptr=ini; ptr<end; ptr++){
	if(ptr->res1 > res1)break;
	if(ptr->res1 == res1){
	  if(ptr->res2 == res2){
	    res12[n12]=ptr; n12++;
	    if(n12 >= MINT){
	      printf("ERROR, more than %d interactions ", MINT);
	      printf("for residues %d and %d\n", res1, res2); exit(8);
	    }
	  }else if((ptr->res2 < res2min)&&(ptr->res2 > res2)){
	    res2min=ptr->res2;
	  }
	}
      } // End list

      if(n12){
	Select_pairs(selected, res12, n12);
	for(k=0; k<n12; k++){
	  if(selected[k]==1){tmp[Nnew]=*res12[k]; Nnew++;}
	}
      }

      // Prepare next value
      if(res2min==nres){ini=ptr; break;} // end res1
      else{res2=res2min-1;}

    }
  }

  for(k=0; k<Nnew; k++){
    Int_list[k]=tmp[k];
  }
  *N_int=Nnew;
  free(tmp);

}

void Purge_interactions_res_pair(int *N_int, struct interaction *Int_list,
				 int nres)
{
  // Maximum one interaction per residue pair
  int Nnew=0, res1, res2, k; 
  struct interaction *res12;
  struct interaction *ptr, *ini=Int_list, *end=Int_list+*N_int;
  struct interaction *tmp=malloc(*N_int*sizeof(struct interaction));

  for(res1=0; res1<nres; res1++){
    for(res2=res1+1; res2<nres; res2++){
      // Find strongest res1-res2 contact
      // Find the next value of res2
      int res2min=nres; res12=NULL;
      for(ptr=ini; ptr<end; ptr++){
	if(ptr->res1 > res1)break;
	if(ptr->res1 == res1){
	  if(ptr->res2 == res2){
	    if((res12==NULL)||(ptr->sec_der>res12->sec_der)||
	       ((ptr->sec_der==res12->sec_der)&&(ptr->rmin < res12->rmin))){
	      res12=ptr;
	    }
	  }else if((ptr->res2 < res2min)&&(ptr->res2 > res2)){
	    res2min=ptr->res2;
	  }
	}
      } // End pair

      if(res12){tmp[Nnew]=*res12; Nnew++;} // store

      // Prepare next value
      if(res2min==nres){ini=ptr; break;} // end res2
      else{res2=res2min-1;}

    } // end res2
  } // end res1

  for(k=0; k<Nnew; k++)Int_list[k]=tmp[k];
  *N_int=Nnew;
  free(tmp);

}


void Select_pairs(int *selected, struct interaction **res12, int n12)
{
  int k, kopt=0, nsel=0; float w=0;
  //selected: 0=not yet -1=discarded 1=selected
  for(k=0; k<n12; k++)selected[k]=0;
  while(kopt>=0){
    kopt=-1;
    for(k=0; k<n12; k++){
      if(selected[k])continue;
      if(res12[k]->sec_der==1){
	if((kopt<0)||(res12[k]->rmin < w)){
	  kopt=k; w=res12[k]->rmin; 
	}
      }else{
	if((kopt<0)||(res12[k]->sec_der > w)){
	  kopt=k; w=res12[k]->sec_der;
	}
      } 
    }
    if(kopt<0)return;
    selected[kopt]=1; nsel++;
    int i1=res12[kopt]->i1, i2=res12[kopt]->i2;
    for(k=0; k<n12; k++){
      if(selected[k])continue;
      if((res12[k]->i1==i1)||(res12[k]->i2==i2))selected[k]=-1;
    }
  }
  return;
}

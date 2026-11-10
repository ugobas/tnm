#include "nma_para.h"
#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "interactions_tnm.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Var_dist.h"

extern int INT_ATOM;
extern float K_THR; // Where the force constant vanishes
extern float K1_THR; // < K_THR, start of interpolatin region
float C0_CA=6.0; // distance for normalizing contact interactions (CA)
float C0_MIN=3.5; //3.5; // distance for normalizing contact interactions (MIN)
//int IJ_MIN=2; // Minimum distance for computing backbone atom contacts
#define MINT 1000 // Maximum number of screened interactions
                  // for any residue pair

/*******************************  Go model  *********************************/
int Cont_comparison(struct interaction c1, struct interaction c2);
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
void Heavy_atoms_interactors(struct interaction *Int_list, int N_int,
			     atom *atoms);
int Heavy_atom(int *iatom, atom *atoms);

void Compute_for_cons_old(struct interaction *Int_list, int N_int, float expo);
void Purge_interactions_atom_res(int *N_int, struct interaction *Int_list,
				 int nres);
void Purge_interactions_res_pair(int *N_int, struct interaction *Int_list,
				 int nres);
void Purge_interactions_res_pair_new(int *N_int, struct interaction *Int_list,
				     int nres);
void Select_pairs(int *selected, struct interaction **res12, int n12);

float Dist2(double *r1, double *r2);
static float Bond_angle(double *r1, double *r2, double *r3);
int Check_store_HB(struct interaction *Int_list, int *N_int,
		   atom *atoms, int i_AA, int i_A, int i_D, int rA,
		   float d2);

/*****************************************************************************/
void Compute_interactions(int *N_int, struct interaction **Int_list,
			  char *INT_TYPE, atom *atoms, int natoms,
			  int nres, char *nameout)
{
  INT_MAX=100*nres;
  struct interaction *Int_tmp=malloc(INT_MAX*sizeof(struct interaction));
  if(Int_tmp==NULL){
    printf("ERROR, no memory available for storing %d interactions\n",
	   INT_MAX); exit(8);
  }
  if(C_THR<1){
    printf("ERROR, Rc= %.2g\n", C_THR); exit(8);
  }

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
      //*N_int=Interactions_all(Int_tmp, C_THR, atoms, natoms, nres);
      *N_int=Interactions_min_M(Int_tmp, C_THR, atoms, natoms, nres);
    }
  }else if(strncmp(INT_TYPE, "HYD", 3)==0){
    *N_int=Interactions_min(Int_tmp, C_THR, atoms, natoms, nres);
    //Heavy_atoms_interactors(Int_tmp, *N_int,atoms);
  }else if(strncmp(INT_TYPE, "SCR", 3)==0){ // S_TYPE=3
    *N_int=Screened_Interactions(Int_tmp,atoms,natoms,nres,3,S_TOL,SCR_THR);
    K_THR=SCR_THR;
  }else if(strncmp(INT_TYPE, "SHA", 3)==0){ // S_TYPE=0 (
    *N_int=Screened_Interactions(Int_tmp,atoms,natoms,nres,
				 S_TYPE,S_TOL,SCR_THR);
    K_THR=SCR_THR;
  }else if(strncmp(INT_TYPE, "HB", 2)==0){
    *N_int=Interactions_HB(Int_tmp, C_THR, atoms, natoms, nres,
			   DA_THR, COS_DAA, ENE_HB);
  }else{
    printf("WARNING, incorrect interaction type %s\n", INT_TYPE);
    printf("Setting to default MIN\n"); strcpy(INT_TYPE, "MIN");
    *N_int=Interactions_min(Int_tmp, C_THR, atoms, natoms, nres);
  }
  printf("%d interactions of type %s found\n", *N_int, INT_TYPE);

  // Each atoms only interacts once per residue
  if((ONEINT)&&
     ((strncmp(INT_TYPE, "SCR", 3)==0)||(strncmp(INT_TYPE, "SHA", 3)==0))){
    //Purge_interactions_atom_res(N_int, Int_list, nres);
    Purge_interactions_res_pair(N_int, Int_tmp, nres);
    printf("%d interactions after keeping 1 per residue pair\n", *N_int);
  }

  int k;
  for(k=0; k<*N_int; k++){
    Int_tmp[k].res1=atoms[Int_tmp[k].i1].res;
    Int_tmp[k].res2=atoms[Int_tmp[k].i2].res;
  }

  if(UNION || INTERSECT){
    struct interaction *Int_tmp1=malloc(*N_int*sizeof(struct interaction));
    for(k=0; k<*N_int; k++){Int_tmp1[k]=Int_tmp[k];}
    struct interaction *Int_tmp2=malloc(INT_MAX*sizeof(struct interaction));
    int N_int2=Interactions_min(Int_tmp2, C_THR, atoms, natoms, nres);
    for(k=0; k<N_int2; k++){
      Int_tmp2[k].res1=atoms[Int_tmp2[k].i1].res;
      Int_tmp2[k].res2=atoms[Int_tmp2[k].i2].res;
    }

    if(UNION){printf("Union ");}
    if(INTERSECT){printf("Intersection ");}
    printf("of %d screened (%.1f) and %d MIN (%.1f) interactions:\n",
	   *N_int, SCR_THR, N_int2, C_THR);
    int m=0, j, k2=0, N_intersect=0, N_end=0;
    for(k=0; k<*N_int; k++){
      int intersect=0;
      for(j=k2; j<N_int2; j++){ 
	int cc=Cont_comparison(Int_tmp2[j], Int_tmp1[k]);
	if(cc<0){
	  if(UNION){Int_tmp[m]=Int_tmp2[j]; m++;}
	}else if(cc==0){
	  Int_tmp[m]=Int_tmp1[k]; m++; intersect=1; N_intersect++;
	  k2=j+1; break;
	}else{
	  k2=j; break;
	}
      } // end j
      if(UNION && intersect==0){Int_tmp[m]=Int_tmp1[k]; m++;}
    } // end k
    if(UNION){// Last contacts
      for(j=k2; j<N_int2; j++){Int_tmp[m]=Int_tmp2[j]; m++; N_end++;}
    }
    free(Int_tmp1);
    free(Int_tmp2);
    if(m>INT_MAX){
      printf("Error, too many contacts %d UNION= %d INTERSECT= %d\n",
	     m, UNION, INTERSECT); exit(8);
    }
    *N_int=m; // Corrected, 5/10/2026
    printf("%d interactions after UNION or INTERSECTION\n", *N_int);
    printf("%d intersections, %d terminal pairs\n", N_intersect, N_end);
  } // end UNION || INTERSECT

  C0=C0_MIN;
  if((strncmp(INT_TYPE,"CA",2)==0)||(strncmp(INT_TYPE,"CB",2)==0)||HNM)
    C0=C0_CA;

  Compute_for_cons(Int_tmp, *N_int, nameout);

  //Compute_for_cons(Int_tmp, *N_int, POW, EXP_HESSIAN, C0, C1_THR, C_THR);
  //Compute_for_cons_old(Int_tmp, *N_int, EXP_HESSIAN);

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
  //float t_max2=t_max+22; t_max2*=t_max2;
  float d, d2;
  int res1, res2;

  for(res1=0; res1<N_res; res1++){

    if(Find_atom(atoms, &i1, res1, N_atoms, atom_type)< 0)continue;
    atom1=atoms+i1; i2=i1+1;

    for(res2=res1+1; res2<N_res; res2++){
      if(Find_atom(atoms, &i2, res2, N_atoms, atom_type)< 0)continue;
      int d12=res2-res1;
      // Interatomic distance
      atom2=atoms+i2;
      d=atom2->r[0]-atom1->r[0]; d2 =d*d; if(d12>INT_STIFF && d2>thr2)continue;
      d=atom2->r[1]-atom1->r[1]; d2+=d*d; if(d12>INT_STIFF && d2>thr2)continue;
      d=atom2->r[2]-atom1->r[2]; d2+=d*d; if(d12>INT_STIFF && d2>thr2)continue;
      
      // store interaction
      struct interaction *Int=Int_list+N_int;
      Int->i1=i1;
      Int->i2=i2;
      Int->r0=sqrt(d2);
      Int->force_aa=1;
      Int->for_cons=1;
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    }
  }
  return(N_int);
}

int Cont_comparison(struct interaction c1, struct interaction c2){
  if(c1.res1 < c2.res1){return(-1);}
  else if(c1.res1 == c2.res1){
    if(c1.res2 < c2.res2){return(-1);}
    else if(c1.res2 == c2.res2){return(0);}
    else{return(1);}
  }else{
    return(1);
  }
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
      int d12=res2-res1;
      atom2=atoms+i2;
      // Interatomic distance
      d=atom2->r[0]-atom1->r[0]; d2 =d*d; if(d12>INT_STIFF && d2>thr2)continue;
      d=atom2->r[1]-atom1->r[1]; d2+=d*d; if(d12>INT_STIFF && d2>thr2)continue;
      d=atom2->r[2]-atom1->r[2]; d2+=d*d; if(d12>INT_STIFF && d2>thr2)continue;

      // store interaction
      struct interaction *Int=Int_list+N_int;
      Int->i1=i1;
      Int->i2=i2;
      Int->r0=sqrt(d2);
      if(d2<=4){
	Int->for_cons=205*sqrt(d2)-571.2;
      }else{
	Int->for_cons=305.9*pow(d2, -3);
      }
      Int->force_aa=Int->for_cons;
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
  float thr2=thr*thr;
  //float t_max2=thr+22; t_max2*=t_max2;
  float d, d2, d2_min;
  int res1, res2, ij;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++;
    i2_start=i1_start+1; ij=0;

    for(res2=res1+1; res2<N_res; res2++){ // Omit i+1 because covalent
      int d12=res2-res1;
      while(atoms[i2_start].res< res2)i2_start++;

      int i1_min=i1_start, i2_min=i2_start, ncheck=0;
      d2_min=10000; ij++;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	//if(atom1->name[0]=='H')continue;
	//int ib; if(ij<IJ_MIN)ib=Backbone(atom1->name);
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;
	  //if(atom2->name[0]=='H')continue;

	  // Check if covalent
	  if(NOCOV && Covalent(res1, res2, atom1, atom2)){continue;}
	  //if((ij<IJ_MIN)&&(ib)&&(Backbone(atom2->name)))continue;
	  ncheck++;

	  // Interatomic distance
	  d=atom2->r[0]-atom1->r[0];
	  d2 =d*d; // if(d2>thr2){continue;}//if(d2 > t_max2)goto next2;  
	  d=atom2->r[1]-atom1->r[1];
	  d2+=d*d; // if(d2>thr2){continue;} //if(d2 > t_max2)goto next2; 
	  d=atom2->r[2]-atom1->r[2];
	  d2+=d*d; // if(d2>thr2){continue;} //if(d2 > t_max2)goto next2; 
	  if(d2 < d2_min){d2_min=d2; i1_min=i1; i2_min=i2;}
	}
      }

      // store interaction
      if(d2_min < thr2 || d12<= INT_STIFF){
	struct interaction *Int=Int_list+N_int;
	Int->i1=i1_min;
	Int->i2=i2_min;
	Int->r0=sqrt(d2_min);
	Int->for_cons=1;
	Int->force_aa=1;
	N_int++;
	if(N_int >= INT_MAX){
	  printf("ERROR, too many interactions %d\n", N_int); exit(8);
	}
      }
      //next2: continue;
    }
  }
  return(N_int);
}

void Heavy_atoms_interactors(struct interaction *Int_list, int N_int,
			     atom *atoms)
{
  int k, h=0;
  for(k=0; k<N_int; k++){
    h+=Heavy_atom(&(Int_list[k].i1), atoms);
    h+=Heavy_atom(&(Int_list[k].i2), atoms);
    if(1){
      Int_list[k].r0=
	sqrt(Distance_square(atoms[Int_list[k].i1].r,atoms[Int_list[k].i2].r));
    }

  }
  printf("%d hydrogen atoms changed to heavy atoms\n", h);
}

int Heavy_atom(int *iatom, atom *atoms)
{
  if(atoms[*iatom].name[0]!='H')return(0); // Not hydrogen atoms
  atom *atom1=atoms+*iatom;
  int i, imin=*iatom, res=atom1->res;
  double dmin=100, *r=atom1->r;
  char branch=atom1->name[1];
  if((branch=='N')||(branch=='1')||(branch=='2')||(branch=='3'))branch=' ';
  for(i=*iatom-1; i>=0; i--){
    atom1--; if(atom1->res!=res)break;
    if((atom1->name[0]!='H')&&(atom1->name[1]==branch)){
      double x=r[0]-atom1->r[0], d=x*x;
      x=r[1]-atom1->r[1]; d+=x*x;
      x=r[2]-atom1->r[2]; d+=x*x;
      if(d<dmin){dmin=d; imin=i;}
    }
  }

  if(*iatom==imin){
    printf("WARNING, atom %s %d, no heavy atom found\n",
	   atoms[*iatom].name, *iatom); return(0);
  }else if(dmin > 4){
    printf("WARNING, atom %s %d, heavy atom %s %d has d=%.3f\n",
	   atoms[*iatom].name, *iatom, atoms[imin].name, imin, dmin);
  }
  *iatom=imin; return(1);
}


int Interactions_min_M(struct interaction *Int_list, float thr,
		     atom *atoms, int N_atoms, int N_res)
{
  //  For each pair of residues only N_RESRES atoms
  //  with smallest distance interact
  int N_int=0, i1, i2, i1_start=0, i2_start=0;
  atom *atom1=atoms, *atom2;
  float thr2=thr*thr;
  float t_max2=thr+22; t_max2*=t_max2;
  float d, d2;
  int res1, res2, k;

  float *d2_min=malloc(N_RESRES*sizeof(float));
  int *i1_min=malloc(N_RESRES*sizeof(int));
  int *i2_min=malloc(N_RESRES*sizeof(int));

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++;
    i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      while(atoms[i2_start].res< res2)i2_start++;
      int n_int=0;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;

	  // Check if covalent
	  if(NOCOV && Covalent(res1, res2, atom1, atom2)){continue;}

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
	struct interaction *Int=Int_list+N_int;
	Int->i1=i1_min[k];
	Int->i2=i2_min[k];
	Int->r0=sqrt(d2_min[k]);
	Int->for_cons=1;
	Int->force_aa=1;
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
  if(res2!=(res1+1)){return(0);}
   if(strncmp(atom2->name, "N ", 2)==0 &&
      strncmp(atom1->name, "C ", 2)==0){return(1);}
  return(0);
}

/*
int Covalent_old(int res1, int res2, atom *atom1, atom *atom2){
      // Old version: 
  if(res2==(res1+1))return(0);
  if((strncmp(atom2->name, "N ", 2)==0)&&
     ((strncmp(atom1->name, "O ", 2)==0)
      ||(strncmp(atom1->name, "CA", 2)==0)
      ||(strncmp(atom1->name, "C ", 2)==0)
      ))return(1);
  if((strncmp(atom2->name, "CA", 2)==0)&&
     (strncmp(atom1->name, "C ", 2)==0)
     //||(strncmp(atom1->name, "O ", 2)==0)
     )return(1);
  return(0);
}
*/

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
  float thr2=thr*thr;
  float t_max2=thr+22; t_max2*=t_max2;
  float d, d2;
  int res1, res2;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++;
    i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      int d12=res2-res1;
      while(atoms[i2_start].res< res2)i2_start++;

      // Check if contact
      for(i1=i1_start; i1<N_atoms; i1++){
	atom1=atoms+i1; if(atom1->res != res1)break;
	for(i2=i2_start; i2<N_atoms; i2++){
	  atom2=atoms+i2; if(atom2->res != res2)break;

	  // Check if covalent
	  if(d12==1 && NOCOV &&  Covalent(res1, res2, atom1, atom2)){continue;}

	  // Interatomic distance
	  d=atom2->r[0]-atom1->r[0]; d2 =d*d;
	  if(d2>thr2 && d12>INT_STIFF){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[1]-atom1->r[1]; d2+=d*d;
	  if(d2>thr2 && d12>INT_STIFF){if(d2 > t_max2)goto next2; continue;}
	  d=atom2->r[2]-atom1->r[2]; d2+=d*d;
	  if(d2>thr2 && d12>INT_STIFF){if(d2 > t_max2)goto next2; continue;}
	  if(d2 < thr2 || d12<=INT_STIFF){
	    struct interaction *Int=Int_list+N_int;
	    Int->i1=i1;
	    Int->i2=i2;
	    Int->r0=sqrt(d2);
	    Int->for_cons=1;
	    Int->force_aa=1;
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
  float thr2=thr*thr;
  float t_max2=thr+22; t_max2*=t_max2;
  float d, d2, d2_min;
  int res1, res2;
  int CB1=0, CB2=0; //i1_min, i2_min, 
  int *notfound=malloc(N_res*sizeof(int));

  for(res1=0; res1<N_res; res1++)notfound[res1]=0;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++;
    CB1=i1_start;
    if(Find_atom(atoms, &CB1, res1, N_atoms, REF)<0){
      if(notfound[res1]==0)
	printf("WARNING, %s atom not found, res= %d\n", REF, res1);
      continue;
    }
    i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      int d12=res2-res1;
      while(atoms[i2_start].res< res2)i2_start++;
      CB2=i2_start;
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
	  if(NOCOV && Covalent(res1, res2, atom1, atom2)){continue;}
	  
	  // Interatomic distance
	  if(d12>INT_STIFF){
	    d=atom2->r[0]-atom1->r[0];
	    d2 =d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	    d=atom2->r[1]-atom1->r[1];
	    d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	    d=atom2->r[2]-atom1->r[2];
	    d2+=d*d; if(d2>thr2){if(d2 > t_max2)goto next2; continue;}
	    if(d2 < d2_min)d2_min=d2; //i1_min=i1; i2_min=i2;}
	  }
	}
      }

      // store interaction
      if(d2_min < thr2 || d12<=INT_STIFF){
	// Interacting atoms set as CB (CA if Gly)
	struct interaction *Int=Int_list+N_int;
	Int->i1=CB1;
	Int->i2=CB2;
	//Int->r0=sqrt(Dist2(atoms[CB1].r, atoms[CB2].r));
	Int->r0=sqrt(d2_min);
	Int->for_cons=1;
	Int->force_aa=1;
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
  float thr2=thr*thr;
  float thr2_NO=thr_NO*thr_NO;
  float t_max2=thr+22; t_max2*=t_max2;
  float d, d2, d2_min;
  int res1, res2;
  int i1_min=-1, i2_min=-1, i_hb;

  for(res1=0; res1<N_res; res1++){

    while(atoms[i1_start].res< res1)i1_start++;
    i2_start=i1_start+1;

    for(res2=res1+1; res2<N_res; res2++){
      int d12=res2-res1;
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
	  if(NOCOV && Covalent(res1, res2, atom1, atom2)){continue;}

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
      if((i_hb)||((d2_min >= thr2)&& d12>INT_STIFF))goto next2;
      struct interaction *Int=Int_list+N_int;
      Int->type=1;
      Int->i1=i1_min;
      Int->i2=i2_min;
      Int->r0=sqrt(d2_min);
      Int->for_cons=1;
      Int->force_aa=1;
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    next2: continue;
    }
  }
  return(N_int);
}


void Compute_for_cons(struct interaction *Int_list, int N_int, char *nameout)
{
  // Normalize such that the average contact interaction does
  // not depend on C and E

  // Set parameters for force constant computation
  printf("Computing %d force constants\n", N_int);

  K1_THR=K_THR; //-0.025; // Lower threshold for interpolation
  if(K_THR>K1_THR){K1=KAPPA/(K_THR-K1_THR);}
  else{K1=0;}

  if(POW){ // power law: C0^-e
    K0=KAPPA*pow(C0,EXP_HESSIAN); 
    K1*=pow((K1_THR/C0),-EXP_HESSIAN);
    printf("Power law k= %.3g (r/%.1f)^-%.1f\n", KAPPA, C0, EXP_HESSIAN);
  }else{
    K0=KAPPA*exp(C0*EXP_HESSIAN);
    K1*=exp((C0-K1_THR)*EXP_HESSIAN);
    printf("Exponential law k= %.3g exp(-(r/%.1f)*%.1f\n",
	   KAPPA, C0, EXP_HESSIAN);
  }
  printf("Interpolating to zero between %.2f and %.2f\n", K1_THR, K_THR);

  char name_out[200];
  if(nameout){
    sprintf(name_out,"%s.Inter_Force_constant.dat",nameout);
    // YYY change output file name
    //if(POW){sprintf(name_out, "%s_POW_%.1f.dat", name_out, EXP_HESSIAN);}
    //else{sprintf(name_out, "%s_EXP_%.1f.dat", name_out, EXP_HESSIAN);}
    FILE *file_out=fopen(name_out, "w");

    if(INT_ATOM){
      fprintf(file_out, "# Force constant inversely proportional to variance "
	      "of atom types distance\n");
      printf("# Force constant inversely proportional to variance "
	     "of atom types distance\n");
    }    
    fprintf(file_out, "# POW= %d EXPO= %.2f C0= %.2f KAPPA= %.2f "
	    "K_THR= %.2f K1_THR= %.2f %d interactions\n",
	    POW, EXP_HESSIAN, C0, KAPPA, K_THR, K1_THR, N_int);

    printf("# Force constant: POW= %d EXPO= %.2f C0= %.2f KAPPA= %.2f "
	   "K_THR= %.2f K1_THR= %.2f\n",
	   POW, EXP_HESSIAN, C0, KAPPA, K_THR, K1_THR);


    if(POW){ // power law: C0^-e
      fprintf(file_out, "# f(r)= %.2f (%.2f/r)^%.2f if r< %.2f\n",
	      KAPPA, C0, EXP_HESSIAN, K1_THR);
    }else{ // Exponential
      fprintf(file_out, "# f(r)= %.2f exp((%.2f-r)%.2f) if r< %.2f\n",
	      KAPPA, C0, EXP_HESSIAN, K1_THR);
    }
    fprintf(file_out, "# f(r)= 0 r>= %.2f\n", K_THR);
    fprintf(file_out, "# f(r)= linear interpolation %.2f < r < %.2f\n",
	    K1_THR, K_THR);
    
    for(float r=2.0; r<K_THR; r+=0.2){
      float kappa=Force_constant(r);
      fprintf(file_out, "%.2f %.3g\n", r, kappa);
    }
    fclose(file_out);
  }

  int n_err=0;
  struct interaction *Int=Int_list;
  for(int n=0; n<N_int; n++){
    if(Int->r0<=0){
      printf("ERROR at pair %d, d2= %.2f\n", n, Int->r0);
      n_err++; Int++; continue;
    }
    Int->for_cons*=Force_constant(Int->r0);
    Int++;
  }
  if(n_err){exit(8);}
  if(0){
    printf("\nExpo= %.2f\nV''(r^PDB):\n", EXP_HESSIAN);
    for(int n=0; n<N_int; n++){
      printf("%d-%d %.2g %.2g\n", Int_list[n].res1, Int_list[n].res2,
	     Int_list[n].r0, Int_list[n].for_cons);
    }
    exit(8);
  }
}


float Force_constant(float r)
{
  if(r <= K1_THR){
    if(POW){return(K0*pow(r, -EXP_HESSIAN));}
    else{return(K0*exp(-r*EXP_HESSIAN));}
  }else if(r < K_THR){
    return(K1*(K_THR-r));
  }else{
    return(0);
  }
}

void Type_dependent_force(struct interaction *Int_list, int N_int,
			  atom *atoms, int natoms,
			  struct residue *res, int nres)
{
  if(INT_ATOM==0)return;
  double sum=0;
  struct interaction *Int=Int_list;
  for(int n=0; n<N_int; n++){
    float var=0; int t1=-1, t2=-1;
    if(INT_ATOM>0){
      int i1=Int->i1, i2=Int->i2;
      if(i1<0 || i2<0 || i1>=natoms || i2>=natoms){
	printf("ERROR, wrong interacting atoms %d %d, n=%d\n", i1,i2,natoms);
	exit(8);
      }
      t1=atoms[i1].type; t2=atoms[i2].type; 
      var=Var_dist_atom[t1][t2];
    }else if(INT_ATOM<0){
      int r1=Int->r1, r2=Int->r2;
      if(r1<0 || r2<0 || r1>=nres || r2>=nres){
	printf("ERROR, wrong interacting res %d %d, n=%d\n", r1,r2,nres);
	exit(8);
      }
      t1=res[r1].i_aa; t2=res[r2].type;
      var=Var_dist_res[t1][t2];
    }
    if(var >0){
      //printf("%.3f %.3f\n", Int->for_cons, Int->for_cons/var);
      Int->force_aa=1/var;
      sum+=Int->force_aa;
    }else{
      printf("ERROR, wrong force constant "
	     "INT_ATOM= %d t1=%d t2=%d var= %.2g\n",
	     INT_ATOM, t1, t2, var);
      exit(8);
    }
    Int++;
  }
  double factor=N_int/sum;
  Int=Int_list;
  for(int n=0; n<N_int; n++){
    Int->force_aa*=factor;
    Int->for_cons*=Int->force_aa;
    Int++;
  }
}

float Dist2(double *r1, double *r2){
  return(r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]);
}

int Check_store_HB(struct interaction *Int_list, int *N_int,
		   atom *atoms, int i_AA, int i_A, int i_D, int rA,
		   float d2)
{
  float cos_tau=Bond_angle(atoms[i_AA].r, atoms[i_A].r, atoms[i_D].r);
  if(cos_tau > COS_DAA)return(0);
  float d_DAA=Dist2(atoms[i_AA].r, atoms[i_D].r);
  struct interaction *Int=Int_list+*N_int;
  Int->type=1;
  Int->force_aa=1;
  Int->for_cons=1;
  Int->r0=sqrt(d_DAA); 
  if(rA==1){
    Int->i1=i_AA;
    Int->i2=i_D;
  }else{
    Int->i1=i_D;
    Int->i2=i_AA;
  }
  (*N_int)++;
  Int++;
  Int->type=2;
  Int->force_aa=ENE_HB;
  Int->for_cons=ENE_HB;
  if(rA==1){
    Int->i1=i_A; Int->i2=i_D;
  }else{
    Int->i1=i_D; Int->i2=i_A;
  }
  Int->r0=sqrt(d2);
  (*N_int)++;
  return(1);
}

float Bond_angle(double *r1, double *r2, double *r3)
{
  if((r1==NULL)||(r2==NULL)||(r3==NULL))return(0);
  double ux=(r1[0]-r2[0]), vx=(r3[0]-r2[0]);
  double uy=(r1[1]-r2[1]), vy=(r3[1]-r2[1]);
  double uz=(r1[2]-r2[2]), vz=(r3[2]-r2[2]);
  double scal_prod=ux*vx+uy*vy+uz*vz;
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
    int ij=abs(atoms[Int_list[i].i1].res-atoms[Int_list[i].i2].res);
    if(ij<=2)continue;
    fprintf(file_out, "%s\t%s\t%.2f",
	    res[atoms[Int_list[i].i1].res].pdbres,
	    res[atoms[Int_list[i].i2].res].pdbres,
	    Int_list[i].r0);
    if(SCR)fprintf(file_out, "\t%.4f", Int_list[i].for_cons);
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
	    if((res12==NULL)||(ptr->for_cons>res12->for_cons)||
	       ((ptr->for_cons==res12->for_cons)&&(ptr->r0 < res12->r0))){
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

void Purge_interactions_res_pair_new(int *N_int, struct interaction *Int_list,
				     int nres)
{
  // Maximum one interaction per residue pair
  int Nnew=0, res1, res2, k; 
  struct interaction *ini=Int_list, *end=Int_list+*N_int;
  struct interaction *tmp=malloc(*N_int*sizeof(struct interaction));

  for(res1=0; res1<nres; res1++){
    for(res2=res1+1; res2<nres; res2++){
      // Find strongest res1-res2 contact
      // Find the next value of res2
      int new1;
      struct interaction *res12=
	Check_pair(&new1, res1, &res2, nres, &ini, end);

      if(res12){tmp[Nnew]=*res12; Nnew++;} // store
      if(new1)break;

    } // end res2
  } // end res1

  for(k=0; k<Nnew; k++)Int_list[k]=tmp[k];
  *N_int=Nnew;
  free(tmp);
}

struct interaction *Check_pair(int *new1, int res1, int *res2, int nres,
			       struct interaction **ini,
			       struct interaction *end)
{
  struct interaction *res12=NULL, *ptr;
  int res2min=nres; *new1=0;
  for(ptr=*ini; ptr<end; ptr++){
    if(ptr->res1 > res1)break;
    if(ptr->res1 == res1){
      if(ptr->res2 == *res2){
	if((res12==NULL)||(ptr->for_cons>res12->for_cons)||
	   ((ptr->for_cons==res12->for_cons)&&(ptr->r0 < res12->r0))){
	  res12=ptr;
	}
      }else if((ptr->res2 < res2min)&&(ptr->res2 > *res2)){
	res2min=ptr->res2;
      }
    }
  } // End pair
  // Prepare next value
  if(res2min==nres){*ini=ptr; *new1=1;} // end res2
  else{*res2=res2min-1;}
  return(res12);
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
      if(res12[k]->for_cons==1){
	if((kopt<0)||(res12[k]->r0 < w)){
	  kopt=k; w=res12[k]->r0; 
	}
      }else{
	if((kopt<0)||(res12[k]->for_cons > w)){
	  kopt=k; w=res12[k]->for_cons;
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

int Backbone(char *name){
  if((strncmp(name, "N ", 2)==0)||
     (strncmp(name, "CA", 2)==0)||
     (strncmp(name, "C ", 2)==0))//||
     //(strncmp(name, "O ", 2)==0)||
     //(strncmp(name, "H ", 2)==0)||
     //(strncmp(name, "HN", 2)==0))
    return(1);
  return(0);
}

void Compute_for_cons_old(struct interaction *Int_list, int N_int, float expo)
{
  double t=0; int n;
    for(n=0; n<N_int; n++){
      if(Int_list[n].r0<=0){
	printf("ERROR, d2= %.2f\n", Int_list[n].r0); exit(8);
      }
      if(expo)Int_list[n].for_cons*=pow(Int_list[n].r0, -expo);
      t+=Int_list[n].for_cons;
    }
    t/=N_int;
    for(n=0; n<N_int; n++)Int_list[n].for_cons/=t;
}

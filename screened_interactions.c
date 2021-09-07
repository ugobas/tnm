#include "coord.h"
#include "tnm.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interactions_tnm.h"
#include "nma.h"

#define RELAX_COV 0  // Less screening from covalently bound atoms

static int RMAX=2000; // Maximum number of stored neighbors 
//static float DTHR=4.6; // Minimum weight: 2*exp(-DTHR)=0.02
static int PROTECT=0; // Protect from screening MIN interactions?


struct store{
  int i;
  float x;
  struct store *next;
};
struct ranked{
  int num;
  struct store *store;
  struct store *first;
  struct store *last;
};

int **Protect_interactions(struct ranked *rank, atom *atoms,
			   int N_atoms, int nres);
int Screening_distance_old(float *e2_nb, float *e2_bo,
			   int res1,struct store *store1,
			   atom *atom2,atom *atom);
int Screening_distance(float *d2_nb, float *d2_bo, float D2_NB, float D2_BO,
		       struct store *store1, float d2_ij,
		       atom *atom_i, atom *atom_j, atom *atoms);
int Shadow_cosine(float *d_nb, float *d_bo, float S_NB, float S_BO,
		   struct store *store1, float d2_ij,
		   atom *atom_i, atom *atom_j, atom *atoms);
int Axis_distance(float *d_nb, float *d_bo, float S2_NB, float S2_BO,
		   struct store *store1, float d2_ij,
		   atom *atom_i, atom *atom_j, atom *atoms);
int Axis_distance_S(float *d_nb, float *d_bo, float S_THR, float S_BOND,
		     struct store *store1, float d2_ij,
		     atom *atom_i, atom *atom_j, atom *atoms);

static float Dist2(double *ri, double *rj);
static float Scalar(double *ri, double *rj, double *rk);
static void  Rank(struct ranked *rank, int i, float x, int Rmax);

int Screened_Interactions(struct interaction *Int_list,
			  float THR, atom *atoms, int N_atoms, int nres,
			  int TYPE, float S_THR)
{
  // For each pair of atoms ij, determine e_ij such that
  // e_ij=min_k(max(d_ik,d_jk) (minimum of the maximum distance).
  // C_ij=exp(e_ij-d_ij)/(1+exp(e_ij-d_ij))-> 0 if e_ij << d_ij
  // Fast algorithm: for each i store the first RMAX neighbors d<thr
  // For each i and j in list of neighbors, find e_ij

  int i1, i2, i;
  atom *atom1=atoms, *atom2;
  float Thr=THR, thr2=Thr*Thr, d, d2;

  // Allocate and rank RMAX neighbors, including bonded
  struct ranked *rank=malloc(N_atoms*sizeof(struct ranked));
  for(i=0; i<N_atoms; i++){
    rank[i].num=0;
    rank[i].store=malloc(RMAX*sizeof(struct store));
  }
  for(i1=0; i1<N_atoms; i1++){
    atom2=atom1+1;
    struct ranked *rank1=rank+i1;
    for(i2=i1+1; i2<N_atoms; i2++){
      // Interatomic distance
      d=atom2->r[0]-atom1->r[0]; d2 =d*d; if(d2>thr2)goto next;
      d=atom2->r[1]-atom1->r[1]; d2+=d*d; if(d2>thr2)goto next;
      d=atom2->r[2]-atom1->r[2]; d2+=d*d; if(d2>thr2)goto next;
      Rank(rank1,  i2, d2, RMAX);
      Rank(rank+i2,i1, d2, RMAX);
    next:
      atom2++;
    }
    atom1++;
  }

  // Thresholds
  float S_BOND=S_THR;
  if(RELAX_COV){
    if(TYPE==0){S_BOND=S_THR+(1.-S_THR)/2;} // cosine
    else if(TYPE==1){S_BOND=0.5*S_THR;} // distance from axis
    else if(TYPE==2){S_BOND=0.5*S_THR;} // Shadow
    else if(TYPE==3){S_BOND=2*S_THR;}   // SCR
  }
  float S2_NB=S_THR*S_THR, S2_BO=S_BOND*S_BOND; 
  //float delta_nb=DTHR*S_THR;
  //float delta_bo=DTHR*S_BOND;
  printf("Screened interactions type %d D=%.3f (nobond) %.3f (bond)\n",
	 TYPE, S_THR, S_BOND);
  if((TYPE<0)||(TYPE>3)){
    printf("ERROR, undefined type of screened interaction %d\n", TYPE);
    exit(8);
  }

  // Protect from screening interactions that satisfy the MIN criterion
  int **protect=NULL;
  if(PROTECT)protect=Protect_interactions(rank, atoms, N_atoms, nres);
  // WARNING! store is a linked list !!!!!

  int N_int=0, tested=0, Nprotect=0;
  float e2=1;
  for(i1=0; i1<N_atoms; i1++){
    atom1=atoms+i1;
    int res1=atom1->res, n=-1;
    struct store *store=rank[i1].first;
    while(store != NULL){
      if(protect)n++;
      i2=store->i;
      if(i2<=i1)goto No_int;
      atom2=atoms+i2;
      int res2=atom2->res;
      if(res2==res1)goto No_int;
      if(NOCOV && Covalent(res1, res2, atom1, atom2))goto No_int;
      if((protect)&&(protect[i1][n])){Nprotect++; goto store;}
      tested++;

      float d2_ij= store->x; 
      float d2_bo=100, d2_nb=d2_bo;  // Bonded and non-bonded, squared

      if(TYPE==3){
	// SCR interactions
	float d_ij=sqrt(d2_ij);
	float D2_NB=d_ij-S_THR;  D2_NB*=D2_NB;
	float D2_BO=d_ij-S_BOND; D2_BO*=D2_BO;
	//Screening_distance(&e2_nb,&e2_bo,res2,rank[i2].first, atom1, atoms);
	// Change 27/3/2015: weight SCR=0,1
	if(Screening_distance(&d2_nb, &d2_bo, D2_NB, D2_BO,
			      rank[i1].first, d2_ij, atom1, atom2, atoms))
	  goto No_int;

	/* 
	e2_nb=sqrt(d2_ij)-sqrt(e2_nb);
	if(e2_nb > delta_nb)goto No_int;
	e2_bo=sqrt(d2_ij)-sqrt(e2_bo);
	if(e2_bo > delta_bo)goto No_int;
	if(S_THR > 0){
	  e2_nb=exp(-e2_nb/S_THR);
	  e2_bo=exp(-e2_bo/S_BOND);
	  if(e2_nb<e2_bo){e2=e2_nb/(1.+e2_nb);}
	  else{e2=e2_bo/(1.+e2_bo);}
	} */

      }else if(TYPE==0){
	// Screening based on cosine
	float cos2_bo=0, cos2_nb=0;
	if(Shadow_cosine(&cos2_nb, &cos2_bo, S2_NB, S2_BO,
			 rank[i1].first, d2_ij, atom1, atom2, atoms))
	  goto No_int;

      }else if(TYPE==1){
	// Screening based on distance from axis
	if(Axis_distance(&d2_nb, &d2_bo, S2_NB, S2_BO,
			 rank[i1].first, d2_ij, atom1, atom2, atoms))
	  goto No_int;
	//w=exp(-d_min/S_THR); w=(1-w)/(1+w); if(w<0.005)goto No_int;

      }else if(TYPE==2){
	// Screening based on Onuchic radius
	if(Axis_distance_S(&d2_nb, &d2_bo, S_THR, S_BOND,
			   rank[i1].first, d2_ij, atom1, atom2, atoms))
	  goto No_int;
      }

    store:
      Int_list[N_int].i1=i1; Int_list[N_int].i2=store->i;
      Int_list[N_int].r0=sqrt(store->x); Int_list[N_int].sec_der=e2;
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    No_int:
      store=store->next;
    }
  }

  // INform
  printf("%d pairs tested, %d screened interactions type %d found,",
	 tested, N_int, TYPE);
  printf(" screening par=%.4f cut-off distance=%.1f", S_THR, THR);
  //if(TYPE==3)printf("cut-off for screening=%.2g", 2*exp(-DTHR));
  printf("\n");
  printf("%d pairs protected from screening because of MIN\n", Nprotect);

  for(i=0; i<N_atoms; i++)free(rank[i].store);
  free(rank);

  return(N_int);
}

int Screening_distance_old(float *d2_nb, float *d2_bo,
			   int res1, struct store *store1,
			   atom *atom2, atom *atoms)
{
  int res2=atom2->res;
  while(store1 != NULL){
    atom *atom3=atoms+store1->i;
    if(atom3!=atom2){
      float d_13=store1->x, d_23=Dist2(atom2->r, atom3->r);
      int bond=Bond(atom3->res, res1, d_13, res2, d_23);
      if(d_23 > d_13){
	if(bond){if(d_23 < *d2_bo)*d2_bo=d_23;}
	else{if(d_23 < *d2_nb)*d2_nb=d_23;}
      }else{
	if(bond){if(d_13 < *d2_bo)*d2_bo=d_13;}
	else{if(d_13 < *d2_nb)*d2_nb=d_13;}
      }
    }
    store1=store1->next;
  }
  return(0);
}

int Screening_distance(float *d2_nb, float *d2_bo, float D2_NB, float D2_BO,
		       struct store *store1, float d2_ij,
		       atom *atom_i, atom *atom_j, atom *atoms)
{
  int resi=atom_i->res, resj=atom_j->res;
  while(store1 != NULL){
    if(store1->x > d2_ij)break; // store1->x= d2_ik
    atom *atom_k=atoms+store1->i;
    if(atom_k!=atom_j){
      float d2, d2_jk=Dist2(atom_j->r, atom_k->r);
      if(d2_jk >= d2_ij)goto next3;
      if(d2_jk > store1->x){
	d2=d2_jk;
      }else{
	d2=store1->x;
      }
      if(d2 < D2_NB){
	int bond=Bond(atom_k->res, resi, store1->x, resj, d2_jk);
	if(bond){
	  if(d2 < D2_BO){*d2_bo=d2; return(1);}
	}else{
	  if(d2 < D2_NB){*d2_nb=d2; return(1);}
	}
      }
    }
  next3:
    store1=store1->next;
  }
  return(0);
}

int Shadow_cosine(float *c_nb, float *c_bo, float S2_NB, float S2_BO,
		   struct store *store1, float d2_ij,
		   atom *atom_i, atom *atom_j, atom *atoms)
{
  int resi=atom_i->res, resj=atom_j->res;
  float cos2_a, d2_jk;
  while(store1 != NULL){
    if(store1->x > d2_ij)break; // store1->x= d2_ik
    atom *atom_k=atoms+store1->i;
    if(atom_k!=atom_j){
      d2_jk=Dist2(atom_j->r, atom_k->r);
      if(d2_jk >= d2_ij)goto next;
      if(d2_jk > store1->x){
	cos2_a=Scalar(atom_j->r, atom_i->r, atom_k->r);
	cos2_a=(cos2_a*cos2_a)/(d2_ij*d2_jk); 
      }else{
	cos2_a=Scalar(atom_i->r, atom_j->r, atom_k->r);
	cos2_a=(cos2_a*cos2_a)/(d2_ij*store1->x);
      }
      if(cos2_a > S2_NB){
	int bond=Bond(atom_k->res, resi, store1->x, resj, d2_jk);
	if(bond){
	  if(cos2_a > S2_BO){*c_bo=cos2_a; return(1);}
	}else{
	  if(cos2_a > S2_NB){*c_nb=cos2_a; return(1);}
	}
      }
    }
  next:
    store1=store1->next;
  }
  return(0);
}

int Axis_distance(float *d2_nb, float *d2_bo, float S2_NB, float S2_BO,
		  struct store *store1, float d2_ij,
		  atom *atom_i, atom *atom_j, atom *atoms)
{
  int resi=atom_i->res, resj=atom_j->res;
  float d2, d2_jk, d2_c;
  while(store1 != NULL){
    if(store1->x > d2_ij)break;
    atom *atom_k=atoms+store1->i;
    if(atom_k!=atom_j){
      d2_jk=Dist2(atom_j->r, atom_k->r);
      if(d2_jk >= d2_ij)goto next1;
      if(d2_jk > store1->x){
	d2_c=Scalar(atom_j->r, atom_i->r, atom_k->r);
	d2=d2_jk-d2_c*d2_c/d2_ij; 
      }else{
	d2_c=Scalar(atom_i->r, atom_j->r, atom_k->r);
	d2=store1->x-d2_c*d2_c/d2_ij;
      }
      if(d2 < S2_NB){
	int bond=Bond(atom_k->res, resi, store1->x, resj, d2_jk);
	if(bond){
	  if(d2 < S2_BO){*d2_bo=d2; return(1);}
	}else{
	  if(d2 < S2_NB){*d2_nb=d2; return(1);}
	}
      }
    }
  next1:
    store1=store1->next;
  }
  return(0);
}

int Axis_distance_S(float *d2_nb, float *d2_bo, float S_NB, float S_BO,
		    struct store *store1, float d2_ij,
		    atom *atom_i, atom *atom_j, atom *atoms)
{
  int resi=atom_i->res, resj=atom_j->res;
  float d2, d2_jk, d2_k, d2_c, S, S2;
  float S2_NB=S_NB*S_NB;
  float S2_BO=S_BO*S_BO;
  float S4_NB=4*S2_NB;
  float S4_BO=4*S2_BO;
  float d_ij=sqrt(d2_ij);
  while(store1 != NULL){
    if(store1->x > d2_ij)break;
    atom *atom_k=atoms+store1->i;
    if(atom_k!=atom_j){
      d2_jk=Dist2(atom_j->r, atom_k->r);
      if(d2_jk >= d2_ij)goto next2;
      if(d2_jk > store1->x){
	d2_k=d2_jk;
	d2_c=Scalar(atom_j->r, atom_i->r, atom_k->r);
      }else{
	d2_k=store1->x;
	d2_c=Scalar(atom_i->r, atom_j->r, atom_k->r);
      }
      d2=d2_k-d2_c*d2_c/d2_ij;
      // d2 is the squared distance between atom k and the axis ij
      if(d2 > S4_NB)goto next2; // No screening if d > 2*S
      int bond=Bond(atom_k->res, resi, store1->x, resj, d2_jk);
      if(bond){
	if(d2 > S4_BO)goto next2;
	if(d2 < S2_BO){*d2_bo=0; return(1);} // There is shadow on the axis
	S=S_BO; S2=S2_BO;
      }else{
	if(d2 < S2_NB){*d2_nb=0; return(1);}
	S=S_NB; S2=S2_NB;
      }
      float ang=asin(sqrt(d2/d2_k)); if(ang<0)ang=-ang;
      float ang1=asin(S/d_ij); if(ang1<0)ang1=-ang1;
      float ang2=ang-ang1, s2;
      if(ang2<=0){*d2_nb=0; return(1);}
      else{s2=sin(ang2); s2*=s2;} // s2=sin(ang2)^2
      if(d2_jk > store1->x){d2=s2*d2_jk;}
      else{d2=s2*store1->x;}
      if(d2<S2){*d2_nb=0; return(1);}

      /*float d2_ic=d2_k-d2;       // Sq. dist. from atom i to projection of k
      // Sq. distance from bottom of atom k to axis
      if(bond){d2=(d2+S2_BO-2*S_BO*sqrt(d2));}
      else{d2=(d2+S2_NB-2*S_NB*sqrt(d2));}
      float d2a_ij=d2*d2_ij/d2_ic; // Sq. dist. from shadow to axis at j
      if(bond){if(d2a_ij < *d2_bo)*d2_bo=d2a_ij;}
      else{if(d2a_ij < *d2_nb)*d2_nb=d2a_ij;}*/
    }
  next2:
    store1=store1->next;
  }
  return(0);
}

float Dist2(double *ri, double *rj){
  double d=ri[0]-rj[0], d2=d*d;
  d=ri[1]-rj[1]; d2+=d*d;
  d=ri[2]-rj[2]; d2+=d*d;
  return(d2);
}

float Scalar(double *ri, double *rj, double *rk){
  double dijk=0; int a;
  for(a=0; a<3; a++){
    double dij=ri[a]-rj[a], dik=ri[a]-rk[a];
    dijk+=dij*dik;
  } 
  return(dijk);
}

void Rank(struct ranked *rank,  int i, float x, int Rmax)
{
  // Rank in decreasing order (lowest = 0)
  if((rank->num >= Rmax)&&(x >= rank->last->x))return;

  struct store *store;
  if(rank->num==0){
    store=rank->store;
    store->i=i; store->x=x;
    store->next=NULL;
    rank->num=1;
    rank->first=store;
    rank->last= store;
    return;
  }

  // Find position in the list
  struct store *k=rank->first, *k1=k->next;
  while(k1!=NULL){
    if(k1->x >= x)break;
    k=k1; k1=k1->next;
  }

  if(rank->num < Rmax){
    store=rank->store+rank->num;
    (rank->num)++;
    if(x > rank->last->x)rank->last=store;
  }else if(k1==NULL){
    store=k; rank->last=k;
  }else{
    store=rank->last;
    if(k1==rank->last){
      k1=NULL;
    }else{
      struct store *k2=k1, *k3=k1->next->next;
      while(k3!=NULL){k2=k2->next; k3=k3->next;}
      rank->last=k2; k2->next=NULL;
    }
  }
  if(x < rank->first->x){
    rank->first=store; k1=k;
  }else{
    k->next=store;
  }
  store->i=i; store->x=x; store->next=k1;

  /*printf("%d: ", rank->num);
  struct store *k2=rank->first; float y=0;
  while(k2!=NULL){
    printf("%.2f ", k2->x); 
    if(k2->x < y){
      printf("\nERROR bad ordered list\n");
      exit(8);
    }
    y=k2->x; k2=k2->next;
  }
  printf("\n");*/

  return;
}

int **Protect_interactions(struct ranked *rank, atom *atoms,
			   int N_atoms, int nres)
{
  int **protect=malloc(N_atoms*sizeof(int *));
  int N_protect=0, i1=0, i;
  int *ii2=malloc(N_atoms*sizeof(int));
  int *i1start=malloc(nres*sizeof(int));
  for(i=0; i<nres; i++)i1start[i]=N_atoms;
  for(i1=0; i1<N_atoms; i1++){
    int n=rank[i1].num;
    protect[i1]=malloc(n*sizeof(int));
    for(i=0; i<n; i++)protect[i1][i]=0;
    if(i1<i1start[atoms[i1].res])i1start[atoms[i1].res]=i1;
    ii2[i1]=0;
  }

  float dmin=100;
  int i1min=-1, imin=-1, res1, res2;
  for(res1=0; res1<nres; res1++){
    for(res2=res1+1; res2<nres; res2++){
      // Find strongest res1-res2 contact
      // Find the next value of res2
      int res2next=nres; i1min=-1;
      for(i1=i1start[res1]; i1<N_atoms; i1++){
	if(atoms[i1].res > res1)break;
	if(atoms[i1].res == res1){
	  struct store *store=rank[i1].first; i=0;
	  while(store!=NULL){
	    int res2atom=atoms[store->i].res;
	    if(res2atom == res2){
	      if(NOCOV && Covalent(res1, res2, atoms+i1, atoms+store->i))
		goto next;
	      if((i1min<0)|| (store->x< dmin)){
		i1min=i1; dmin=store->x; imin=i;
	      }
	    }else if(res2atom > res2){
	      if(res2atom  < res2next)res2next=res2atom;
	      //ii2[i1]=i; break;
	    }
	  next:
	    store=store->next; i++;
	  }
	}
      } // end pair
      if(i1min >=0){protect[i1min][imin]=1; N_protect++;}

      // Prepare next value
      if(res2next<nres)res2=res2next-1;
    }
  }
  printf("%d contacts protected from screening (MIN)\n", N_protect);
  free(ii2); free(i1start);
  return(protect);
}

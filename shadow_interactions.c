#include "coord.h"
#include "tnm.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nma.h"
#include "interactions_tnm.h"

#define RELAX_COV 1  // Less screening from covalently bound atoms

static int RMAX=300; // Maximum number of stored neighbors 
static float DTHR=5.3; // Minimum weight: 2*exp(-DTHR)=0.005

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

void Shadow_cosine(float *d_nb, float *d_bo,
		   struct store *store1, float d2_ij,
		   atom *atom_i, atom *atom_j, atom *atoms);
void Axis_distance(float *d_nb, float *d_bo,
		   struct store *store1, float d2_ij,
		   atom *atom_i, atom *atom_j, atom *atoms);
void Axis_distance_S(float *d_nb, float *d_bo, float S_THR,
		     struct store *store1, float d2_ij,
		     atom *atom_i, atom *atom_j, atom *atoms);
static float Scalar(float *ri, float *rj, float *rk);
static float Dist2(float *ri, float *rj);
static void Rank(struct ranked *rank, int i, float x, int Rmax);

int Interactions_shadow(struct interaction *Int_list, float THR, 
			atom *atoms, int N_atoms,
			int TYPE, float S_THR, int NOCOV)
{
  // For each pair of atoms ij, determine e_ij such that
  // e_ij=min_k(max(d_ik,d_jk) (minimum of the maximum distance).
  // C_ij=exp(e_ij-d_ij)/(1+exp(e_ij-d_ij))-> 0 if e_ij << d_ij
  // Fast algorithm: for each i store the first RMAX neighbors d<thr
  // For each i and j in list of neighbors, find e_ij

  int i1, i2, i;
  atom *atom1=atoms, *atom2;
  float Thr=THR, thr2=Thr*Thr, d=0, d2;
  //float EPS=0.005, A_THR=-log((1-EPS)/(1+EPS))*S_THR;

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

  // Threshold
  float S_BOND=S_THR; // More tolerant threshold for covalent bonds
  if(RELAX_COV){
    if(TYPE==0){S_BOND=2*S_THR; if(S_BOND > 0.99)S_BOND=0.99;}
    else if(TYPE==1){S_BOND=0.5*S_THR;}
    else if(TYPE==2){S_BOND=0.5*S_THR;}
    else{printf("ERROR, undefined type of shadow interaction %d\n", TYPE);
      exit(8);
    }
  }
  printf("Testing screened interactions, type %d\n", TYPE);

  int N_int=0, tested=0;
  for(i1=0; i1<N_atoms; i1++){
    atom1=atoms+i1;
    int res1=atom1->res;
    struct store *store=rank[i1].first;
    while(store != NULL){
      i2=store->i;
      if(i2<=i1)goto No_int;
      atom2=atoms+i2;
      if(atom2->res==res1)goto No_int;
      if(NOCOV && Covalent(res1, atom2->res, atom1, atom2))goto No_int;
      tested++;
      float d2_ij= store->x, w=1;
      if(TYPE==0){
	// Screening based on cosine
	float cos_bo=0, cos_nb=0;
	Shadow_cosine(&cos_nb, &cos_bo,
		      rank[i1].first, d2_ij, atom1, atom2, atoms);
	//printf("Shadow cosine= %.3f (nobond) %.3f (bond)\n",cos_nb,cos_bo);
	if(cos_nb > S_THR)goto No_int;
	if(cos_bo > S_BOND)goto No_int;
      }else if(TYPE==1){
	// Screening based on distance from axis
	float d_bo=10, d_nb=d_bo;
	Axis_distance(&d_nb, &d_bo,
		      rank[i1].first, d2_ij, atom1, atom2, atoms);
	//printf("d=  %.3f (nobond) %.3f (bond)\n", d_nb, d_bo);
	if(d_nb < S_THR)goto No_int;
	if(d_bo < S_BOND)goto No_int;
	//w=exp(-d_min/S_THR); w=(1-w)/(1+w); if(w<0.005)goto No_int;
      }else if(TYPE==2){
	// Screening based on Onuchic radius
	float d_bo=2*S_THR, d_nb=d_bo;
	Axis_distance_S(&d_nb, &d_bo,
			S_THR, rank[i1].first, d2_ij, atom1, atom2, atoms);
	//printf("d=  %.3f (nobond) %.3f (bond)\n", d_nb, d_bo);
	if(d_nb < S_THR)goto No_int;
	if(d_bo < S_BOND)goto No_int;
      }
      Int_list[N_int].i1=i1; Int_list[N_int].i2=i2;
      Int_list[N_int].r2=d2_ij; Int_list[N_int].sec_der=w;
      Int_list[N_int].rmin=d2_ij;
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    No_int:
      store=store->next;
    }
  }


  // INform
  printf("%d pairs tested, %d shadow interactions type %d found,",
	 tested, N_int, TYPE);
  printf(" parameter=%.3g cut-off distance=%.1f cut-off for screening=%.2g\n",
	 S_THR, THR, 2*exp(-DTHR));

  for(i=0; i<N_atoms; i++)free(rank[i].store);
  free(rank);

  return(N_int);
}

void Shadow_cosine(float *c_nb, float *c_bo,
		   struct store *store1, float d2_ij,
		   atom *atom_i, atom *atom_j, atom *atoms)
{
  int resi=atom_i->res, resj=atom_j->res;
  float cos_a, d2_jk;
  while(store1 != NULL){
    if(store1->x > d2_ij)break;
    atom *atom_k=atoms+store1->i;
    if(atom_k!=atom_j){
      d2_jk=Dist2(atom_j->r, atom_k->r);
      if(d2_jk >= d2_ij)goto next;
      if(d2_jk > store1->x){
	cos_a=Scalar(atom_j->r, atom_i->r, atom_k->r)/sqrt(d2_ij*d2_jk); 
      }else{
	cos_a=Scalar(atom_i->r, atom_j->r, atom_k->r)/sqrt(d2_ij*store1->x);
      }
      int bond=Bond(atom_k->res, resi, store1->x, resj, d2_jk);
      if(bond){if(cos_a > *c_nb)*c_bo=cos_a;}
      else{if(cos_a > *c_nb)*c_nb=cos_a;}
    }
  next:
    store1=store1->next;
  }
  *c_nb=sqrt(*c_nb);
  *c_bo=sqrt(*c_bo);
}

float Dist2(float *ri, float *rj){
  float d=ri[0]-rj[0];
  float d2=d*d;
  d=ri[1]-rj[1]; d2+=d*d;
  d=ri[2]-rj[2]; d2+=d*d;
  return(d2);
}

float Scalar(float *ri, float *rj, float *rk){
  float dijk=0; int a;
  for(a=0; a<3; a++){
    float dij=ri[a]-rj[a], dik=ri[a]-rk[a];
    dijk+=dij*dik;
  } 
  return(dijk);
}

void Axis_distance(float *d_nb, float *d_bo,
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
      int bond=Bond(atom_k->res, resi, store1->x, resj, d2_jk);
      if(bond){if(d2 < *d_bo)*d_bo=d2;}
      else{if(d2 < *d_nb)*d_nb=d2;}
    }
  next1:
    store1=store1->next;
  }
  *d_bo=sqrt(*d_bo);
  *d_nb=sqrt(*d_nb);
}

void Axis_distance_S(float *d_nb, float *d_bo, float S_THR,
		     struct store *store1, float d2_ij,
		     atom *atom_i, atom *atom_j, atom *atoms)
{
  int resi=atom_i->res, resj=atom_j->res;
  float d2, d2_jk, d2_k, d2_c;
  float S2=S_THR*S_THR;
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
      if(d2 < S2)return;           // There is shadow on the axis 
      d2=sqrt(d2)-S_THR; d2*=d2;   // Sq. distance bottom of k - axis 
      float d2_ic=d2_k-d2;         // Sq. dist. from atom i to projection of k
      float d2a_ij=d2*d2_ij/d2_ic; // Sq. dist. from shadow to axis at j
      int bond=Bond(atom_k->res, resi, store1->x, resj, d2_jk);
      if(bond){if(d2a_ij < *d_bo)*d_bo=d2a_ij;}
      else{if(d2a_ij < *d_nb)*d_nb=d2a_ij;}
    }
  next2:
    store1=store1->next;
  }
  *d_bo=sqrt(*d_bo);
  *d_nb=sqrt(*d_nb);
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
    struct store *k2=k1, *k3=k1->next;
    while(k3!=NULL){k2=k3; k3=k2->next;}
    rank->last=k2;
  }
  if(x < rank->first->x){
    rank->first=store; k1=k;
  }else{
    k->next=store;
  }
  store->i=i; store->x=x; store->next=k1;

  /*printf("%d: ", rank->num);
  struct store *k2=rank->first;
  while(k2!=NULL){printf("%.2f ", k2->x); k2=k2->next;}
  printf("\n");
  if(rank->num == 10){
    k=rank->first; float y=0;
    printf("Ordered list: %d\n", rank->num);
    while(k!=NULL){
      printf("%.2f ", k->x);
      if(k->x <= y){
	printf("\nERROR bad ordered list\n");
	exit(8);
      }
      y=k->x; k=k->next;
    }
    printf("\n");
    exit(8);
    }*/

  return;
}

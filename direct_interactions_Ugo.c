#include "coord.h"
#include "tnm.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "direct_interactions.h"
#include "nma.h"

int RMAX=100;
float Thr=15;

struct store{
  int i;
  float x;
  struct store *next;
};
struct ranked{
  int num;
  float max, min;
  struct store *store;
  struct store *first;
};

static float Dist2(atom *atom1, atom *atom2);
static void Rank(struct ranked *rank, int i, float x, int Rmax);

int Interactions_direct(struct interaction *Int_list,
			atom *atoms, int N_atoms)
{
  // For each pair of atoms ij, determine e_ij such that
  // e_ij=min_k(max(d_ik,d_jk) (minimum of the maximum distance).
  // C_ij=exp(e_ij-d_ij)/(1+exp(e_ij-d_ij))-> 0 if e_ij << d_ij
  // Fast algorithm: for each i store the first RMAX neighbors d<thr
  // For each i and j in list of neighbors, find e_ij

  int i1, i2, i;
  atom *atom1=atoms, *atom2;
  float thr2=Thr*Thr, d, d2;

  // Allocate and rank RMAX neighbors, including bonded
  struct ranked *rank=malloc(N_atoms*sizeof(struct ranked));
  for(i=0; i<N_atoms; i++){
    rank[i].num=0; rank[i].max=0; rank[i].min=10000;
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

  float delta=5*DELTA_DIR;
  int N_int=0, tested=0;
  for(i1=0; i1<N_atoms; i1++){
    atom1=atoms+i1;
    struct store *store=rank[i1].first;
    while(store != NULL){
      i2=store->i;
      if(i2<=i1)goto No_direct;
      atom2=atoms+i2;
      if((atom2->res==atom1->res)&&(store->x < 5))goto No_direct;
      if(atom2->res==(atom1->res+1)){
	if((strncmp(atom2->name, "N ", 2)==0)&&
	   ((strncmp(atom1->name, "C ", 2)==0)||
	    (strncmp(atom1->name, "O ", 2)==0)))goto No_direct;
	if((strncmp(atom2->name, "CA", 2)==0)&&
	   (strncmp(atom1->name, "C ", 2)==0))goto No_direct;
      }
      tested++;
      float d2_ij= store->x;
      float e2_ij=thr2, e2;
      struct store *store1=rank[i1].first;
      while((store1 != NULL)&&(store1->x < thr2)){
	if(store1!=store){
	  e2=Dist2(atom2, atoms+store1->i);
	  if(e2 > store1->x){
	    if(e2 < e2_ij)e2_ij=e2;
	  }else{
	    if(store1->x < e2_ij)e2_ij=store1->x;
	  }
	}
	store1=store1->next;
      }
      e2=sqrt(d2_ij)-sqrt(e2_ij);
      if(e2 > delta)goto No_direct;
      e2=exp(-e2/DELTA_DIR);
      e2=e2/(1.+e2);
      // Direct interaction
      Int_list[N_int].i1=i1; Int_list[N_int].i2=i2;
      Int_list[N_int].r2=d2; Int_list[N_int].sec_der=e2;
      Int_list[N_int].rmin=d2;
      N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int); exit(8);
      }

    No_direct:
      store=store->next;
    }
  }


  // INform
  printf("%d pairs tested, %d direct interactions found,",tested, N_int);
  printf(" tolerance=%.4f cut-off=%.4f\n", DELTA_DIR,  delta);

  for(i=0; i<N_atoms; i++)free(rank[i].store);
  free(rank);

  return(N_int);
}

float Dist2(atom *atom1, atom *atom2){
  float d2, d;
  d=atom2->r[0]-atom1->r[0]; d2 =d*d;
  d=atom2->r[1]-atom1->r[1]; d2+=d*d;
  d=atom2->r[2]-atom1->r[2]; d2+=d*d;
  return(d2);
}

void Rank(struct ranked *rank,  int i, float x, int Rmax)
{
  // Rank in decreasing order (lowest = 0)
  if((rank->num >= Rmax)&&(x >= rank->max))return;

  struct store *store;
  if(rank->num==0){
    store=rank->store;
    store->i=i; store->x=x;
    store->next=NULL;
    rank->num=1;
    rank->first=store;
    rank->max=x;
    rank->min=x;
    return;
  }

  // Find position in the list
  struct store *k=rank->first, *k1=k->next;
  while(k1!=NULL){
    if(k1->x >= x)break;
    k=k1; k1=k1->next;
  }

  if(rank->num < Rmax){
    if(x > rank->max)rank->max=x;
    store=rank->store+rank->num;
    (rank->num)++;
  }else if(k1==NULL){
    store=k; rank->max=x;
  }else{
    store=k1; // store1 store store2
    struct store *store1=k, *store2=store->next;
    while(store2!=NULL){
      store1=store; store=store2; store2=store2->next;
    }
    if(x > store1->x){rank->max=x; k1=NULL;}
    else{rank->max=store1->x; store1->next=NULL;}
  }
  if(x < rank->min){
    rank->min=x;
    k1=(rank->first);
    rank->first=store;
  }else{
    k->next=store;
  }
  store->i=i; store->x=x; store->next=k1;

  /*k=rank->first; float y=0;
  while(k!=NULL){
    printf("%f ", k->x);
    if(k->x <= y){
      printf("\nERROR %f %f %f\n", k->x, store->x, k1->x);
      exit(8);
    }
    y=k->x; k=k->next;
  }
  printf("\n");*/

  return;
}

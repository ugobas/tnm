#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "coord.h"
#include "mammoth_interface.h"
#include "mammoth_ali.h"
#include "Allocate_fortran.h"
#define VERBOSE 1

int Sequence_identity(int *nali, char *aseq1, char *aseq2,
		      int *alignres, int n1, int n2);
void  Set_coord(int *nca, float ***xca, float ***vec, char **aseq,
		struct residue *seq, struct chain *chp, int n);
float **Set_CA_vectors(float **XCA, int n);
void Empty_coord(float **xca, float **vec, char *aseq, int n);


float Mammoth_ali(int *aligned, int *seq_id, float *psi, int *alignres, 
		  struct chain *chp1, struct residue *seq1, 
		  struct chain *chp2, struct residue *seq2)
{
  float lne;
  float **xca1, **vec1; char *aseq1; int n1=chp1->nres, nca1;
  float **xca2, **vec2; char *aseq2; int n2=chp2->nres, nca2;

  Set_coord(&nca1, &xca1, &vec1, &aseq1, seq1, chp1, n1);
  Set_coord(&nca2, &xca2, &vec2, &aseq2, seq2, chp2, n2);

  lne=Call_mammoth(alignres, psi, VERBOSE,
		   n1, nca1, xca1, vec1, aseq1,
		   n2, nca2, xca2, vec2, aseq2);
  *seq_id=Sequence_identity(aligned, aseq1, aseq2, alignres, n1, n2);
  Empty_coord(xca1, vec1, aseq1, n1);
  Empty_coord(xca2, vec2, aseq2, n2);
  return(lne);
}

void  Set_coord(int *nca, float ***xca, float ***vec, char **aseq,
		struct residue *seq, struct chain *chp, int n)
{
  int i, n1=n-1, nn=0;
  struct residue *res=seq+chp->ini_res;
  *xca=Allocate_mat2_f_fortran(n,3);
  *aseq=malloc(n*sizeof(char));
  for(i=0; i<n; i++){
    (*aseq)[i]=res->amm;
    atom *atom1=res->atom_ptr;
    int n_atom=res->n_atom, j, k;
    float *ca=(*xca)[nn];
    for(j=0; j<n_atom; j++){
      if(strncmp(atom1->name, "CA", 2)==0){
	for(k=0; k<3; k++)ca[k]=atom1->r[k];
	nn++; break;
      }
      atom1++;
    }
    res++;
  }
  *vec=Set_CA_vectors(*xca, n);
  *nca=nn;
}

float **Set_CA_vectors(float **XCA, int n){
  int nv=n-1, i, k;
  float **vec, *xca=XCA[0], *xca1=XCA[1], *v;

  vec=Allocate_mat2_f_fortran(nv,3);
  for(i=0; i<nv; i++){
    double vmod=0;
    v=vec[i];
    for(k=0; k<3; k++){v[k]=xca1[k]-xca[k]; vmod+=v[k]*v[k];}
    vmod=sqrt(vmod);
    for(k=0; k<3; k++)v[k]/=vmod;
    xca=xca1; xca1=XCA[i+2];
  }
  return(vec);
}

void Empty_coord(float **xca, float **vec, char *aseq, int n){
  Empty_matrix_f_fortran(n-1, vec);
  Empty_matrix_f_fortran(n-1, xca);
  free(aseq); 
}

int Sequence_identity(int *Nali, char *aseq1, char *aseq2,
		      int *alignres, int n1, int n2)
{
  int i, Nid=0; *Nali=0;
  for(i=0; i<n1; i++){
    if(alignres[i]>=0){
      if(aseq1[i]==aseq2[alignres[i]])Nid++;
      (*Nali)++;
    }
  }
  return(Nid);
  //return((float)Nid/(float)Nali);
}

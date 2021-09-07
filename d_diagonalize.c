#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "diagonalize.h"

/* nrutil.h/nrutil.c */
extern double *dvector( long, long);
extern double **dmatrix( long, long, long, long);
extern void free_dmatrix( double **, long, long, long, long);
extern void free_dvector( double *, long, long);
/* tred2.c (float -> double) */
extern void tred2( double **, int, double *, double *);
/* tqli.c (float -> double, bolean return value)  */
extern int tqli( double *, double *, int, double **);
/* jacobi.c (float -> double)  */
extern void jacobi( double **, int, double *, double **, int *);
/* eigsrt.c (float -> double)  */
void eigsrt( double *, double **, int);
void d_sort(double d[], int n, int *i_rank);


void f_Diagonalize(int N, float **MATRIX, float *eigen_values,
		   float **eigen_vector, int SIGN)
{
  double **eigen_vt, *a; float *m;
  int jacobi_rotations=0, i, j, k, ik, ir;
  int *i_rank=malloc(N*sizeof(int));

  double *e_vector= dvector( (long)1, (long)N);
  double *eigen_va= dvector( (long)1, (long)N);
  double **a_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);
  double **v_matrix= NULL;

  /* fill matrix a_matrix. Indeces from 1 to N! */
  for(i=0; i<N; i++){
    a=a_matrix[i+1]; m=MATRIX[i];
    //for(j=0; j<N; j++)a[j+1]= m[j];
    if(SIGN<0){for(j=0; j<N; j++)a[j+1]=-m[j];}
    else{for(j=0; j<N; j++)a[j+1]= m[j];}
  }

  // obtain eigensystem with tridimensional reduction
  tred2( a_matrix, (int)N, eigen_va, e_vector);
  if( tqli(eigen_va, e_vector, (int)N, a_matrix)){
    eigen_vt=a_matrix;
  }else{
    // tqli failed 
    // reconstruct matrix a_matrix (destroyed by tred2/tqli!)
    printf("tqli routine has failed, trying jacobi\n");
    for(i=0; i<N; i++){
      a=a_matrix[i+1]; m=MATRIX[i];
      //for(j=0; j<N; j++)a[j+1]= m[j];
      if(SIGN<0){for(j=0; j<N; j++)a[j+1]=-m[j];}
      else{for(j=0; j<N; j++)a[j+1]= m[j];}
    }
    // obtain eigensystem with jacobi routine
    v_matrix=dmatrix( (long)1, (long)N, (long)1, (long)N);
    jacobi( a_matrix, N, eigen_va, v_matrix, &jacobi_rotations);
    eigen_vt=v_matrix;
  }

  // Sort:
  d_sort(eigen_va+1, N, i_rank);
  if(SIGN<0) for(j=1; j<=N; j++)eigen_va[j]=-eigen_va[j];
  /* d_sort(eigen_va+1, N, i_rank);
     if(SIGN<0) for(j=1; j<=N; j++)eigen_va[j]=-eigen_va[j];*/
  // eigsrt(eigen_va, eigen_vt, N);

  // Transpose and shift indices. Omit lambda < E_THR
  ik=0; // int ilow=N-1;
  for(i=0; i<N; i++){
    //ir=i+1; k=ik; ik++;
    ir=i_rank[i]+1;
    /*if(select && (select[ir]==0)){
      k=ilow; ilow--; if(selected)selected[k]=0;
      }else{*/
    k=ik; ik++;
    //}
    eigen_values[k]=eigen_va[ir];
    for(j=1; j<=N; j++){
      eigen_vector[k][j-1]=eigen_vt[j][ir];
    }
  }

  free_dvector(e_vector, 1, N);
  free_dvector(eigen_va, 1, N);
  free_dmatrix(a_matrix, 1, N, 1, N);
  if(v_matrix)free_dmatrix(v_matrix, 1, N, 1, N);
  free(i_rank);
  return;
}

void d_Diagonalize(int N, double **MATRIX, float *eigen_values,
		   float **eigen_vector, int SIGN)
{
  double **eigen_vt, *a; double *m;
  int jacobi_rotations=0, i, j, k, ik, ilow, ir;
  int *i_rank=malloc(N*sizeof(int));

  double *e_vector= dvector( (long)1, (long)N);
  double *eigen_va= dvector( (long)1, (long)N);
  double **a_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);
  double **v_matrix= NULL;

  /* fill matrix a_matrix. Indeces from 1 to N! */
  for(i=0; i<N; i++){
    a=a_matrix[i+1]; m=MATRIX[i];
    for(j=0; j<N; j++)a[j+1]= m[j];
    /* if(SIGN<0){for(j=0; j<N; j++)a[j+1]=-m[j];}
       else{for(j=0; j<N; j++)a[j+1]= m[j];}*/
  }

  /* obtain eigensystem */
  tred2( a_matrix, (int)N, eigen_va, e_vector);
  if( tqli(eigen_va, e_vector, (int)N, a_matrix)){
    eigen_vt=a_matrix;
  }else{
    /* tqli failed */
    /* reconstruct matrix a_matrix (destroyed by tred2/tqli!) */
    for(i=0; i<N; i++){
      a=a_matrix[i+1]; m=MATRIX[i];
      for(j=0; j<N; j++)a[j+1]= m[j];
      /*if(SIGN<0){for(j=0; j<N; j++)a[j+1]=-m[j];}
	else{for(j=0; j<N; j++)a[j+1]= m[j];}*/
    }

    /* obtain eigensystem */
    v_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);
    jacobi( a_matrix, N, eigen_va, v_matrix, &jacobi_rotations);
    eigen_vt=v_matrix;
  }

  // Sort:
  if(SIGN<0) for(j=1; j<=N; j++)eigen_va[j]=-eigen_va[j];
  d_sort(eigen_va+1, N, i_rank);
  if(SIGN<0) for(j=1; j<=N; j++)eigen_va[j]=-eigen_va[j];
  // eigsrt(eigen_va, eigen_vt, N);

  // Change sign of eigenvectors
  for(j=1; j<=N; j++){
    double sum=0;
    for(i=1; i<=N; i++)sum+=eigen_vt[i][j];
    if(sum<0)for(i=1; i<=N; i++)eigen_vt[i][j]=-eigen_vt[i][j];
  }

  // Transpose and shift indices. Omit lambda < E_THR
  float E_THR=0;
  ik=0; ilow=N-1;
  for(i=0; i<N; i++){
    //ir=i+1; k=ik; ik++;
    ir=i_rank[i]+1;
    if(fabs(eigen_va[ir])<=E_THR){k=ilow; ilow--;}
    else{k=ik; ik++;}
    eigen_values[k]=eigen_va[ir];
    for(j=1; j<=N; j++){
      eigen_vector[k][j-1]=eigen_vt[j][ir];
    }
  }

  free_dvector(e_vector, 1, N);
  free_dvector(eigen_va, 1, N);
  free_dmatrix(a_matrix, 1, N, 1, N);
  if(v_matrix)free_dmatrix(v_matrix, 1, N, 1, N);
  free(i_rank);

  return;
}

void dd_Diagonalize(int N, double **MATRIX, double *eigen_values,
		    double **eigen_vector, int SIGN)
{
  double **eigen_vt, *a; double *m;
  int jacobi_rotations=0, i, j, k, ik, ilow, ir;
  int *i_rank=malloc(N*sizeof(int));

  double *e_vector= dvector( (long)1, (long)N);
  double *eigen_va= dvector( (long)1, (long)N);
  double **a_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);
  double **v_matrix= NULL;

  /* fill matrix a_matrix. Indeces from 1 to N! */
  for(i=0; i<N; i++){
    a=a_matrix[i+1]; m=MATRIX[i];
    for(j=0; j<N; j++)a[j+1]= m[j];
    /* if(SIGN<0){for(j=0; j<N; j++)a[j+1]=-m[j];}
       else{for(j=0; j<N; j++)a[j+1]= m[j];}*/
  }

  /* obtain eigensystem */
  tred2( a_matrix, (int)N, eigen_va, e_vector);
  if( tqli(eigen_va, e_vector, (int)N, a_matrix)){
    eigen_vt=a_matrix;
  }else{
    /* tqli failed */
    /* reconstruct matrix a_matrix (destroyed by tred2/tqli!) */
    for(i=0; i<N; i++){
      a=a_matrix[i+1]; m=MATRIX[i];
      for(j=0; j<N; j++)a[j+1]= m[j];
      /* if(SIGN<0){for(j=0; j<N; j++)a[j+1]=-m[j];}
	 else{for(j=0; j<N; j++)a[j+1]= m[j];}*/
    }

    /* obtain eigensystem */
    v_matrix=dmatrix( (long)1, (long)N, (long)1, (long)N);
    jacobi( a_matrix, N, eigen_va, v_matrix, &jacobi_rotations);
    eigen_vt=v_matrix;
  }

  // Sort:
  if(SIGN<0) for(j=1; j<=N; j++)eigen_va[j]=-eigen_va[j];
  d_sort(eigen_va+1, N, i_rank);
  if(SIGN<0) for(j=1; j<=N; j++)eigen_va[j]=-eigen_va[j];
  // eigsrt(eigen_va, eigen_vt, N);

  // Transpose and shift indices. Omit lambda < E_THR
  float E_THR=0;
  ik=0; ilow=N-1;
  for(i=0; i<N; i++){
    //ir=i+1;  k=ik; ik++;
    ir=i_rank[i]+1;
    if(fabs(eigen_va[ir])<=E_THR){k=ilow; ilow--;}
    else{k=ik; ik++;}
    eigen_values[k]=eigen_va[ir];
    for(j=1; j<=N; j++){
      eigen_vector[k][j-1]=eigen_vt[j][ir];
    }
  }

  free_dvector(e_vector, 1, N);
  free_dvector(eigen_va, 1, N);
  free_dmatrix(a_matrix, 1, N, 1, N);
  if(v_matrix)free_dmatrix(v_matrix, 1, N, 1, N);
  free(i_rank);

  return;
}

void d_sort(double d[], int n, int *i_rank)
{
  // Sort the vector d from large to small.
  // Returns i_rank[i]= index of object at rank i
  int k,j,i, jmax, *not_ranked=malloc(n*sizeof(int));
  double d_max;

  for (i=0; i<n; i++){not_ranked[i]=1; i_rank[i]=i;}
  for (k=0; k<n; k++) {
    for(i=0; i<n; i++)if(not_ranked[i])break;
    d_max=d[jmax=i];
    for (j=i+1;j<n; j++){
      if(not_ranked[j]&&(d[j] > d_max))d_max=d[jmax=j];
    }
    i_rank[k]=jmax; not_ranked[jmax]=0;
  }
  free(not_ranked);
}


void eigsrt(double d[], double **v, int n)
{
  int k,j,i;
  double p;
  
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

void Select_modes(int *selected, double *eigen_va, float E_MIN, int N)
{
  /* printf("Discarding modes with e_val/<e_val> < %.2g ", E_MIN*0.001);
     printf(" or e_val/<e_val> < %.2g, coll < %.f\n",E_MIN, COLL_THR);
     int n=0; */
  float E_THR=E_MIN; int i;
  double e=0; for(i=1; i<=N; i++)e+=eigen_va[i];
  if(e<0)e=-e; E_THR*=(e/N);
  for(i=1; i<=N; i++){
    if(eigen_va[i]<E_THR){
      selected[i]=0;
    }else{
      selected[i]=1;
    }
  }
}

void Select_modes_coll(int *selected, double *eigen_va, double **eigen_vt,
		       float E_MIN, float COLL_THR, int N)
{
  /* printf("Discarding modes with e_val/<e_val> < %.2g ", E_MIN*0.001);
     printf(" or e_val/<e_val> < %.2g, coll < %.f\n",E_MIN, COLL_THR);
     int n=0; */
  float E_THR=E_MIN; int i, j;
  double e=0; for(j=1; j<=N; j++)e+=eigen_va[j];
  if(e<0)e=-e; E_THR*=(e/N);
  float E_THR2=E_THR*0.001;
  for(i=1; i<=N; i++){
    selected[i]=1;
    if(eigen_va[i]<E_THR){
      if(eigen_va[i]<E_THR2){
	selected[i]=0;
      }else{
	double sum=0, Entropy=0;
	for(j=1; j<=N; j++){
	  float v=eigen_vt[j][i];
	  if(v){v*=v; sum+=v; Entropy-=v*log(v);}
	}
	float c=sum*exp(Entropy/sum);
	if(c<COLL_THR){
	  selected[i]=0;
	  //printf("%.1f %f\n", c, eigen_va[i]); n++;
	}
      }
    }
  }
}

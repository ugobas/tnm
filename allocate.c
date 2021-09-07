#include <stdlib.h>
float  **Allocate_mat2_f(int n1, int n2);
double **Allocate_mat2_d(int n1, int n2);

void Empty_matrix_f(float **matrix, int N){
  int i;
  if(matrix==NULL)return;
  for(i=0; i<N; i++)free(matrix[i]);
  free(matrix);
}

void Empty_matrix_d(double **matrix, int N){
  int i;
  if(matrix==NULL)return;
  for(i=0; i<N; i++)free(matrix[i]);
  free(matrix);
}


void Empty_matrix_i(int **matrix, int N){
  int i;
  if(matrix==NULL)return;
  for(i=0; i<N; i++)free(matrix[i]);
  free(matrix);
}

float **Allocate_mat2_f(int n1, int n2){
  int i, j;
  float **matrix=malloc(n1*sizeof(float *));
  for(i=0; i<n1; i++){
    matrix[i]=malloc(n2*sizeof(float));
    for(j=0; j<n2; j++)matrix[i][j]=0;
  }
  return(matrix);
}


int ** Allocate_mat2_i(int n1, int n2){
  int i, j;
  int **matrix=malloc(n1*sizeof(int *));
  for(i=0; i<n1; i++){
    matrix[i]=malloc(n2*sizeof(int));
    for(j=0; j<n2; j++)matrix[i][j]=0;
  }
  return(matrix);
}

double **Allocate_mat2_d(int n1, int n2){
  int i, j;
  double **matrix=malloc(n1*sizeof(double *));
  for(i=0; i<n1; i++){
    matrix[i]=malloc(n2*sizeof(double));
    for(j=0; j<n2; j++)matrix[i][j]=(double)0;
  }
  return(matrix);
}

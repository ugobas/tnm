/************************ List of routines ****************************/
// 3D Vector operations with float
float Scalar_product_3(float *x, float *y);
void Vector_product(float *w, float *v1, float *v2);
void Rotate(float *r, float **rot);
void Subtract_vector_3(float *a, float *b, float *c);
void Sum_vector_3(float *a, float *b, float *c);
void Transform_matrix(float A[3][3], float **R);
//  3D Vector operations with double
double Scalar_product_3_d(double *x, double *y);
void Vector_product_d(double *w, double *v1, double *v2);
void Rotate_d(double *r, double **rot);
void Subtract_vector_3_d(double *a, double *b, double *c);
void Sum_vector_3_d(double *a, double *b, double *c);
void Transform_matrix_d(double **A, double **R);
// generic vector operations
float Scalar_product(float *x, float *y, int n);
float Scalar_product_weighted(float *x, float *y, float *w, int n);
void Normalize_vector(float *vv, int n);
void Normalize_vector_weighted(float *vv, float *w, int n);
void Matrix_multiplication(double *w, double **M, double *v, int n);
// Correlations
float Corr_coeff_w(float *x, float *y, int *w, int N,
		   float *slope, float *offset, int *norm);
float Corr_coeff(float *x, float *y, int N, float *slope, float *offset);
// Collectivity
float Collectivity_norm1(float *v, int N);
float Collectivity_norm2(float *v, int N);
float Collectivity_norm2_outlier(float *v, int N, int *outlier);
float Collectivity_Renyi_norm1(float *v, int N);
float Collectivity_Renyi_norm2(float *v, int N);
float Area_norm1(float *v, int N);
float Area_norm2(float *v, int N);

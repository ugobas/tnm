void f_Diagonalize(int N, float **MATRIX, float *eigen_values,
		   float **eigen_vector, int SIGN);
void d_Diagonalize(int N, double **MATRIX, float *eigen_values,
		   float **eigen_vector, int SIGN);
void dd_Diagonalize(int N, double **MATRIX, double *eigen_values,
		    double **eigen_vector, int SIGN);
void Select_modes(int *selected, double *eigen_va, float E_MIN, int N);
void Select_modes_coll(int *selected, double *eigen_va, double **eigen_vt,
		       float E_MIN, float COLL_THR, int N);

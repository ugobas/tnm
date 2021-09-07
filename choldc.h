int choldc(double **L, double **a, int N);
void Forward_substitution(double *X, double **L, double *Y, int N);
void Backward_substitution(double *X, double **L, double *Y, int N);
int choldc_f(float **L, float **a, int N);
float **Cholevsky_inversion_d2f(double **L, int N);
double **Cholevsky_inversion(double **L, int N);

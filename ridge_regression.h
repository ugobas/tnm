#ifndef RIDGE_REGRESSION
#define RIDGE_REGRESSION

struct ridge_fit{
  float *A;      // Fitted parameters
  float nu;      // Scale of the fit
  float Lambda;  // Ridge parameter
  float E;       // Error of the fit
  float Cv;      // dE/dLambda
  float LF_csi;  // Lambda*sum (A-A_ref(csi))^2
  float LF_inf;  // Lambda*sum (A-A_ref_inf)^2
  float LF_mu;   // Lambda*(1-mu1)*sum (A-A_ref_inf)^2
  float A2;      // sum A^2 ("entropy")
  float Sv1;     // d(Lambda*A2)/dLambda;
  float R2;      // Renyi entropy
  float Sh;      // Shannon entropy
  float GCV;     // Generalized cross validation
  float RR;      // Range risk, outdated
  float dof[3];
};

float Ridge_regression(struct ridge_fit *return_fit, float *Y_pred,
		       float *D_out, char *nameout, float **X, float *Y,
		       int N, int Npar, char TYPE);
#endif

float Fit_no_outliers(int *outlier, float *XX, float *YY, int N,
		      float *slope, float *offset);
float Bfactors_fit(float *r, int *outlier, float *f_dof,
		   float *B_pred_all, float *B_exp, float *B_ENM,
		   float *R_eq, int N, char *name, char type);
int Filter_modes_outliers(struct Normal_Mode NM, struct Reference Ref_kin,
			  float t_out, float t_mode);

#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "allostery.h"
#include "read.h"
#include "vector.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "allocate.h"
#include "diagonalize.h"
#include "interactions_tnm.h"
#include "random3.h"
#include "EC.h"
#include <time.h>

int PRINT_DIFF=1;  // Print difference in Conformation per residue?
int PRINT_SITE=1;  // Print active site read in the PDB?
#define VERBOSE 0
#define DIR_AVE 0 // Compute directionality profile as average besides as PE
#define SIGMA 2.0 // Number of standard deviations to print coupling
int ini_ran=0;

float *PE_profile(double **Corr, int N, int IJMIN, int invert);
float Eigenvalue_main(double **Corr, int Na);
float Singvalue_main(double **A, int Na, int Nb);
void Normalize_profile(double *v, int N);
int Read_sites_PDB(int **site_res, char **chain_res, int *nres_site,
		   char *pdb, char *chain, struct residue *seq,
		   int Nres, int NMAX, int SMAX);
int Read_sites_file(int **site_res, char **chain_res, int *nres_site,
		    char *pdb, char *chain, struct residue *seq,
		    int Nres, int NMAX, int SMAX);
int *Random_sites(int nres, int *site_res, int N);
unsigned long randomgenerator(void);
int Read_site(char *resnum, char *resnam, char *reschain,
	      char *string, char *chain);
int Find_residue(char *resnum, char chain,
		 struct residue *seq, int Nres, int i_ini);
int Print_links(double **Coupling, int N,
		char **pdbres, char *amm, char *chain,
		char *nameout1, char *ext, int sign);
float *Matrix_average(double **Mat, int N, int IJMIN, char *name);
int Test_profile(float *prof_1, float *prof_2, float *pair_1, float *pair_2,
		 float *profile, double **coupling, int Na,
		 int ***site_res, int *nres_site, int N_sites,
		 struct residue *seq, int *atomres,
		 char **pdbres, char *amm, char *chain,
		 char *nameout1, char *filename,
		 char *sitename, char *profname, int num);
double *Cont_ave(float *profile, int N, int **clist, int *nc);
int Extract_atoms(int *iref, int *atomres, struct residue *seq, float *mass,
		  atom *atoms, struct Reference Ref, char *ANAME);
void Get_interaction_list(int **clist, int *nc,
			  atom *atoms, int *atomres, int Na,
			  struct interaction *Int_list, int N_int);
void Pairwise_distances(float ***r0ijk, atom *atoms, int *iref, int Na);
void Fill_low_diag(double **Corr, int Na);
int Print_matrix(double **Coupling, int N, char **pdbres, char *amm,
		 char *nameout1, char *ext);
void Allosteric_propagation(double **Allosteric_coupling, int Na,
			    char **pdbres, char *amm, char *chain,
			    char *nameout);
void Print_sites(int **site_res, int *nres_site, int nsites,
		 struct residue *seq, char *pdb, char *chain);
float Normalized_allostery(float *strdiff,float *Allostery,float *Bpa,int Na);
int Get_zeta(double **Coupling, int N);
double **Negexp_distance(atom *atoms, int *iref, int Na, float D);

/****************************************************************

      Main routine

*****************************************************************/

void Predict_allostery(struct Normal_Mode NM, atom *atoms,
		       struct Reference Ref,
		       struct interaction *Int_list, int N_int,
		       struct residue *seq, char *nameout1, //int mode,
		       float *Confchange,
		       char *pdb, char *chain,
		       int Nres, char *SITES,
		       float *B_pred)
{

  // Extract c_alpha
  float mass=0;
  int *iref=malloc(Ref.N_ref*sizeof(int));
  int *atomres=malloc(Ref.N_ref*sizeof(int));
  char SEL[4]="CA"; if(strcmp(REF, "CB")==0)strcpy(SEL, "CB");
  int Na=Extract_atoms(iref, atomres, seq, &mass, atoms, Ref, SEL);

  // Interaction list
  int *nc=malloc(Na*sizeof(int));
  int **clist=malloc(Na*sizeof(int *));
  Get_interaction_list(clist, nc, atoms, atomres, Na, Int_list, N_int);

  // Pairwise distances
  float ***r0ijk=malloc(Na*sizeof(float **));
  Pairwise_distances(r0ijk, atoms, iref, Na);

  // Allocate
  double **Corr_dr=NULL, **Corr_dir=NULL, **Corr_dr_dir=NULL;
  double **Corr_Bahar=NULL, **Corr_str=NULL, **Corr_str_dir=NULL;
  double **Corr_coord=NULL;
  Corr_dir=Allocate_mat2_d(Na, Na);
  Corr_coord=Allocate_mat2_d(Na, Na);
  if(0)Corr_dr=Allocate_mat2_d(Na, Na);
  if(0)Corr_dr_dir=Allocate_mat2_d(Na, Na);
  if(0)Corr_Bahar=Allocate_mat2_d(Na, Na);
  if(STRAIN){
    Corr_str=Allocate_mat2_d(Na, Na);
    Corr_str_dir=Allocate_mat2_d(Na, Na);
  }

  float *abs_dr=malloc(Na*sizeof(float));
  float *abs_dr2=malloc(Na*sizeof(float));
  float *strain=malloc(Na*sizeof(float));
  float **dr=Allocate_mat2_f(Na, 3);
  float **dir=Allocate_mat2_f(Na, 3);

  // For each normal mode
  int i, j, ia, ja, a, n;
  double norm_w=0;
  for(a=0; a<NM.N_relevant; a++){
    if(NM.contr2fluct[a]==0)continue;
    float *Mode=NM.Cart[a];
    //float w=NM.contr2fluct[a];
    float w=1./(NM.omega2[a]);
    norm_w+=w;

    // Individual atoms
    for(ia=0; ia<Na; ia++){
      int i3=3*iref[ia], j3;
      double r2=0, x;
      for(j=0; j<3; j++){x=Mode[i3+j]; dr[ia][j]=x; r2+=x*x;}
      abs_dr2[ia]=r2;
      r2=sqrt(r2); abs_dr[ia]=r2;
      for(j=0; j<3; j++)dir[ia][j]=dr[ia][j]/r2;
      if(STRAIN){  // Compute strain
	double s=0;
	for(n=0; n<nc[ia]; n++){
	  j3=3*clist[ia][n];
	  for(j=0; j<3; j++){x=Mode[i3+j]-Mode[j3+j]; s+=x*x;}
	}
	strain[ia]=sqrt(s/nc[ia]);
      }
    }

    // Pairwise computations
    for(ia=0; ia<Na; ia++){
      float **r0i=r0ijk[ia];
      for(ja=0; ja<=ia; ja++){
	float dri_drj=0;
	for(j=0; j<3; j++)dri_drj+=dr[ia][j]*dr[ja][j];
	if(Corr_dr_dir)Corr_dr_dir[ia][ja]+=w*dri_drj;       // <ri*rj>
	if(Corr_dr)Corr_dr[ia][ja]+=w*abs_dr[ia]*abs_dr[ja]; // <|ri|*|rj|>
	if(Corr_Bahar)
	  Corr_Bahar[ia][ja]+=w*(abs_dr2[ia]+abs_dr2[ja]-2*dri_drj);
	double sum=0;
	for(j=0; j<3; j++)sum+=r0i[ja][j]*(dr[ia][j]-dr[ja][j]);
	Corr_coord[ia][ja]+=w*(sum*sum);
	float dd=0; for(j=0; j<3; j++)dd+=dir[ia][j]*dir[ja][j];
	Corr_dir[ia][ja]+=w*dd;                 // d_i*d_j   d=r/|r|
	if(STRAIN){
	  float ss=w*strain[ia]*strain[ja];
	  Corr_str[ia][ja]+=ss;                   // strain_i*strain_j
	  Corr_str_dir[ia][ja]+=dd*ss;            // (str_i)d_i*(str_j)d_j
	}
      }
    }
  }

  // Normalize
  for(ia=0; ia<Na; ia++){
    for(ja=0; ja<=ia; ja++){
      Corr_dir[ia][ja]/=norm_w;
      //Corr_coord[ia][ja]=sqrt(Corr_coord[ia][ja]);
    }
    if(Corr_Bahar){
      for(ja=0; ja<ia; ja++)Corr_Bahar[ia][ja]=sqrt(Corr_Bahar[ia][ja]);
    }
  }

  // Fill lower diagonal
  Fill_low_diag(Corr_dir, Na);
  Fill_low_diag(Corr_coord, Na);
  if(Corr_dir)Fill_low_diag(Corr_dir, Na);
  if(Corr_dr_dir)Fill_low_diag(Corr_dr_dir, Na);
  if(Corr_str)Fill_low_diag(Corr_str, Na);
  if(Corr_str_dir)Fill_low_diag(Corr_str_dir, Na);
  if(Corr_Bahar)Fill_low_diag(Corr_Bahar, Na);

  /********************* Print couplings ******************/
  char **pdbres=malloc(Na*sizeof(char *));
  char *amm=malloc(Na*sizeof(char));
  char *ch=malloc(Na*sizeof(char));
  for(ia=0; ia<Na; ia++){
    struct residue *res=seq+(atoms+Ref.atom_num[iref[ia]])->res;
    pdbres[ia]=malloc(6*sizeof(char));
    strcpy(pdbres[ia], res->pdbres);
    amm[ia]=res->amm;
    ch[ia]=res->chain;
  }
  if(PRINT_SIGMA_DIJ){
    Print_links(Corr_coord, Na, pdbres, amm, ch,
		nameout1, "_interatomic_distance_variance.dat", 0);
  }
  
  // Get_zeta(Corr_coord, Na); // Transform into Z score!
  if(PRINT_COORD_COUPLING){
    Print_links(Corr_coord, Na, pdbres, amm, ch,  
		nameout1, "_coordination_coupling.dat", -1); //_zeta
  }
  //Get_zeta(Corr_dir, Na); // Transform into Z score!
  if(PRINT_DIR_COUPLING){
    Print_links(Corr_dir, Na, pdbres, amm, ch,    
		nameout1, "_directionality_coupling.dat", 1); //_zeta
    Print_links(Corr_dir, Na, pdbres, amm, ch, 
		nameout1, "_directionality_coupling_all.dat", 0); //_zeta
  }


  if(0){
    // Broadcaster profile with force that maximizes deformation
    // produced by residue i
    double *Broadcast_profile=malloc(Na*sizeof(double));
    printf("Computing Broadcast profile\n");
    for(ia=0; ia<Na; ia++){
      double **Def=Allocate_mat2_d(3,3);
      int i3=3*iref[ia], i, j;
      for(a=0; a<NM.N_relevant; a++){
	float w=NM.contr2fluct[a]; w*=w;
	float *Mode=NM.Cart[a]+i3;
	for(i=0; i<3; i++)for(j=0; j<=i; j++)Def[i][j]+=w*Mode[i]*Mode[j];
      }
      for(i=0; i<3; i++)for(j=i+1; j<3; j++)Def[i][j]=Def[j][i];
      Broadcast_profile[ia]=Eigenvalue_main(Def, 3);
      Empty_matrix_d(Def, 3);
    }
    Normalize_profile(Broadcast_profile, Na);
  }

  // Receiver profile with force that maximizes deformation
  // suffered by residue i
  float norm_F=mass;
  printf("Computing Allosteric coupling\n");
  double **Allosteric_coupling=Allocate_mat2_d(Na,Na);
  for(ia=0; ia<Na; ia++){
    int i3=3*iref[ia], i, j, k;
    for(ja=0; ja<=ia; ja++){
      int j3=3*iref[ja];
      double **F_matrix=Allocate_mat2_d(3, 3);
      double **Receiver_matrix=Allocate_mat2_d(3, 3);
      for(a=0; a<NM.N_relevant; a++){
	if(NM.omega2[a]<=0)continue;
	float w=1./NM.omega2[a]; //w*=w;
	float *Mode_i=NM.Cart[a]+i3;
	float *Mode_j=NM.Cart[a]+j3;
	for(i=0; i<3; i++)
	  for(j=0; j<3; j++)
	    F_matrix[i][j]+=w*Mode_i[i]*Mode_j[j];
      }
      for(i=0; i<3; i++){
	for(j=0; j<=i; j++){
	  for(k=0; k<3; k++)
	    Receiver_matrix[i][j]+=F_matrix[k][i]*F_matrix[k][j];
	  if(j!=i)Receiver_matrix[j][i]=Receiver_matrix[i][j];
	}
      }
      //float Lambda=Singvalue_main(Receiver_matrix, 3, 3)*norm_F;
      float Lambda=Eigenvalue_main(Receiver_matrix, 3);
      // Rij and Rji have the same eigenvalues, since:
      // Fji=(Fij)t Rij=(Fij)t(Fij) Rji=(Fij)(Fij)t
      // Lambda/=abs_dr2[ia];
      // Optimal forced deformation at i divided by thermal deformation
      if(Lambda<0){
	printf("WARNING, i=%d j=%d Lambda=%.3f\n", ia, ja, Lambda);
	Lambda=0;
      }
      Lambda=sqrt(Lambda)*norm_F; // Module, and not square displacement
      Allosteric_coupling[ia][ja]=Lambda;
      Allosteric_coupling[ja][ia]=Lambda;
      Empty_matrix_d(F_matrix, 3);
      Empty_matrix_d(Receiver_matrix, 3);
    }
  }
  /********* Print allosteric coupling ************/
  Get_zeta(Allosteric_coupling, Na);
  if(PRINT_ALLO_COUPLING){
    Print_links(Allosteric_coupling, Na, pdbres, amm, ch,
	       nameout1, "_allo_coupling_zeta.dat", 1);
  }
  if(PRINT_ALLO_MATRIX){
    Print_matrix(Allosteric_coupling, Na, pdbres, amm,
		 nameout1, "_allosteric_matrix_zeta.dat");
  }


  // Profiles
  int IJMIN=2;
  float *prof_coord=NULL, *prof_dir=NULL, *prof_dr=NULL, *prof_dr_dir=NULL;
  float *prof_Bahar=NULL, *prof_str=NULL, *prof_str_dir=NULL;
  float *prof_dir_ave=NULL;

  // Weighted_flux(Corr_dr, Na, IJMIN, "Corr_dr")  // CV_profile
  // Matrix_average(Corr_dr, Na, IJMIN, "Corr_dr");
  /* prof_dr=EC_profile(NULL, Corr_dr, Na, IJMIN, "Corr_dr",0);
    char profname[5]="EC";*/

  printf("Diagonalizing Corr_coord\n");
  prof_coord=PE_profile(Corr_coord, Na, IJMIN,1);
  printf("Diagonalizing Corr_dir\n");
  prof_dir=PE_profile(Corr_dir, Na, IJMIN, 0);
  if(Corr_dr){
    printf("Diagonalizing Corr_dr\n");
    prof_dr=PE_profile(Corr_dr, Na, IJMIN, 0);
  }
  if(Corr_dr_dir){
    printf("Diagonalizing Corr_dr_dir\n");
    prof_dr_dir=PE_profile(Corr_dr_dir, Na, IJMIN, 0);
  }
  if(Corr_Bahar){
    printf("Diagonalizing Corr_Bahar\n");
    prof_Bahar=PE_profile(Corr_Bahar, Na, IJMIN, 1);
  }
  if(Corr_str){
    printf("Diagonalizing Corr_str\n");
    prof_str=PE_profile(Corr_str, Na, IJMIN, 0);
  }
  if(Corr_str_dir){
    printf("Diagonalizing Corr_str_dir\n");
    prof_str_dir=PE_profile(Corr_str_dir, Na, IJMIN, 0);
  }
  char profname[5]="PE";

  // Average profile
  if(DIR_AVE){
    prof_dir_ave=Matrix_average(Corr_dir, Na, IJMIN, "Corr_dir");
  }

  //float *Receiver_profile=Get_profile(Allosteric_coupling, Na);

    /*
  // Compute the allosteric matrix
  // M_ij = |sum_a w_a ri^a rj^a*fj|^2
  // where fj=sum_a w_a^2 rj/|fj|
  double **Allostery=Allocate_mat2_d(Na, Na);
  float **force_mode=Allocate_mat2_f(Na, NM.N);
  //  Compute force that maximizes deformation
  for(ia=0; ia<Na; ia++){
    double f[3]; for(j=0; j<3; j++)f[j]=0;
    int i3=3*iref[ia];
    for(a=0; a<NM.N; a++){
      float w=NM.contr2fluct[a]; //w*=w;
      float *Mode=NM.Cart[a]+i3;
      for(j=0; j<3; j++)f[j]+=w*Mode[j];
    }
    double f2=f[0]*f[0]+f[1]*f[1]+f[2]*f[2]; f2=sqrt(f2);
    for(j=0; j<3; j++)f[j]/=f2;
    for(a=0; a<NM.N; a++){
      float *Mode=NM.Cart[a]+i3; double ff=0;
      for(j=0; j<3; j++)ff+=f[j]*Mode[j];
      force_mode[ia][a]=ff;
    }
  }
  for(ia=0; ia<Na; ia++){
    int i3=3*iref[ia];
    for(ja=0; ja<Na; ja++){
      double Dr[3]; for(j=0; j<3; j++)Dr[j]=0;
      for(a=0; a<NM.N; a++){
	float *Mode=NM.Cart[a]+i3;
	float ww=NM.contr2fluct[a]*force_mode[ja][a];
	for(j=0; j<3; j++)Dr[j]+=ww*Mode[j];
      }
      Allostery[ia][ja]=Dr[0]*Dr[0]+Dr[1]*Dr[1]+Dr[2]*Dr[2];
    }
  }
  Empty_matrix_f(force_mode, Na);
  float *prof_Allostery=Get_profile(Allostery, Na);
    */


  // Use directionality profile as a weight
  float min_w=10;   double norm_weight=0;
  float *weight=malloc(Na*sizeof(float));
  for(ia=0; ia<Na; ia++){
    if(prof_dir[ia]<min_w)min_w=prof_dir[ia];
    norm_weight+=prof_dir[ia];
  }
  norm_weight-=Na*min_w;
  for(ia=0; ia<Na; ia++){
    weight[ia]=(prof_dir[ia]-min_w)/norm_weight;
  }

  // Allosteric profiles
  float *Allosteric_profile=malloc(Na*sizeof(float));
  for(ia=0; ia<Na; ia++){
    double sum=0; double *Ca=Allosteric_coupling[ia];
    for(ja=0; ja<Na; ja++)sum+=Ca[ja];
    Allosteric_profile[ia]=sum/(Na-1);
  }

  float *Allosteric_profile_weight=malloc(Na*sizeof(float));
  for(ia=0; ia<Na; ia++){
    double sum=0; double *Ca=Allosteric_coupling[ia];
    for(ja=0; ja<Na; ja++)sum+=Ca[ja]*weight[ja];
    Allosteric_profile_weight[ia]=sum/(1.-weight[ia]);
  }

  float *Allosteric_distance_weight=malloc(Na*sizeof(float));
  for(ia=0; ia<Na; ia++){
    double sum=0, norm=0; double *Ca=Allosteric_coupling[ia];
    for(ja=0; ja<Na; ja++){
      int d=ia-ja; if(d<0)d=-d;
      sum+=Ca[ja]*d; norm+=d;
    }
    Allosteric_distance_weight[ia]=sum/norm;
  }

  //double *Allosteric_ave=Cont_ave(Allosteric_profile, Na, clist, nc);
  printf("Diagonalizing Allosteric coupling\n");
  float *Allosteric_PE=PE_profile(Allosteric_coupling, Na, IJMIN, 0);
  printf("Computing EC of allosteric coupling\n");
  float *Allosteric_EC=
    EC_profile(NULL, Allosteric_coupling, Na, IJMIN, "Allostery",0);


  // Conformation change
  float *strdiff=NULL;
  if(Confchange){
    strdiff=malloc(Na*sizeof(float));
    for(ia=0; ia<Na; ia++){
      double d=0; float *x=Confchange+3*iref[ia];
      for(j=0; j<3; j++){d+=(*x)*(*x); x++;}
      strdiff[ia]=d;
    }
  }

  // Print
  char nameout[200];
  sprintf(nameout, "%s_allostery.dat", nameout1);
  FILE *file_out=fopen(nameout, "w");

  printf("Writing %s\n", nameout);
  fprintf(file_out, "#res\taa");
  fprintf(file_out, "\tAllosteric_profile");
  fprintf(file_out, "\tAllosteric_profile_weighted");
  //fprintf(file_out, " Allosteric_profile_contact_ave");
  fprintf(file_out, "\tAllosteric_profile_distance_weight");
  fprintf(file_out, "\tAllosteric_PE");
  if(prof_dir)fprintf(file_out, "\t%s(<di*dj>)", profname);
  if(prof_coord)fprintf(file_out, "\t%s(<rij-r0ij>)", profname);
  fprintf(file_out, "\tAllosteric_EC");
  if(DIR_AVE)fprintf(file_out, "\tAve(<di*dj>)");
  if(prof_dr_dir)fprintf(file_out, "\t%s(<ri*rj>)", profname);
  if(prof_dr)fprintf(file_out, "\t%s(<|ri|*|rj|>)", profname);
  if(prof_Bahar)fprintf(file_out, "\t%s(<|ri-rj|^2>,Bahar)", profname);
  if(prof_str)fprintf(file_out, "\t%s(<s_i*s_j>)", profname);
  if(prof_str_dir)fprintf(file_out, "\t%s(<s_i d_i*s_j d_j>)", profname);
  if(strdiff && PRINT_DIFF)fprintf(file_out, "\tConfchange");
  fprintf(file_out, "\n");

  if(STRAIN)fprintf(file_out, "# with di=ri/|ri|, s=strain");
  if(strncmp(profname, "PE", 2)==0){
    fprintf(file_out, "# PE=principal eigenvector\n");
  }else if(strncmp(profname, "EC", 2)==0){
    fprintf(file_out, "# EC=effective connectivity\n");
  }else  if(strncmp(profname, "WF", 2)==0){
    fprintf(file_out, "# WF=weighted flux\n");
  }
  for(ia=0; ia<Na; ia++){
    fprintf(file_out, "%3s\t%c", pdbres[ia], amm[ia]);
    fprintf(file_out, "\t%8.3g", Allosteric_profile[ia]);
    fprintf(file_out, "\t%8.3g", Allosteric_profile_weight[ia]);
    //fprintf(file_out, "  %8.3g", Allosteric_ave[ia]);
    fprintf(file_out, "\t%8.3g", Allosteric_distance_weight[ia]);
    fprintf(file_out, "\t%8.3g", Allosteric_PE[ia]);
    if(prof_dir)fprintf(file_out, "\t%6.3f", prof_dir[ia]);
    if(prof_coord)fprintf(file_out, "\t%6.3f", prof_coord[ia]);
    fprintf(file_out, "  %8.3g", Allosteric_EC[ia]);
    if(DIR_AVE)fprintf(file_out, "\t%6.3f", prof_dir_ave[ia]);
    if(prof_dr_dir)fprintf(file_out, "\t%6.3f", prof_dr_dir[ia]);
    if(prof_dr)fprintf(file_out, "\t%6.3f", prof_dr[ia]);
    if(prof_Bahar)fprintf(file_out, "\t%6.3f", prof_Bahar[ia]);
    if(prof_str)fprintf(file_out, "\t%6.3f", prof_str[ia]);
    if(prof_str_dir)fprintf(file_out, "\t%6.3f", prof_str_dir[ia]);
    if(strdiff && PRINT_DIFF)fprintf(file_out, "\t%6.3f", strdiff[ia]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);

  //
  Allosteric_propagation(Allosteric_coupling, Na, pdbres, amm, ch, nameout1);

  /********************************************************/
  // Correlation between allostery, conformation change, normal modes
  float *Allostery=malloc(Na*sizeof(float));
  for(ia=0; ia<Na; ia++)Allostery[ia]=Allosteric_profile[ia];
  float *fluctuation=malloc(NM.N_relevant*sizeof(float));
  int imax=0, NMODE=15; if(NMODE>NM.N_relevant)NMODE=NM.N_relevant;
  float *cc=malloc(NMODE*sizeof(float)), cmax=-2;
  float *dev=malloc(Na*sizeof(float));
  for(i=0; i<NMODE; i++){
    float *Cart=NM.Cart[i]; double sum=0;
    for(ia=0; ia<Na; ia++){
      double d=0; float *x=Cart+3*iref[ia];
      for(j=0; j<3; j++){d+=(*x)*(*x); x++;}
      dev[ia]=d; sum+=d;
    }
    cc[i]=Corr_coeff(dev, Allostery, Na, NULL, NULL);
    if(cc[i]>cmax){cmax=cc[i]; imax=i;}
    fluctuation[i]=sqrt(sum/Na)/NM.omega[i];
  }

  // Writing allostery_summary
  sprintf(nameout, "%s_allostery_summary.dat", nameout1);
  file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  fprintf(file_out,
	  "# Maximum correlation allosteric prof. normal mode: %.3f mode=%d\n",
	  cmax, imax);

  if(Confchange){
    float r=Corr_coeff(strdiff, Allostery, Na, NULL, NULL);
    fprintf(file_out,
	    "# Correlation allosteric prof. confchange: %.3f\n", r);
    float *Bpa=malloc(Na*sizeof(float));
    for(ia=0; ia<Na; ia++)Bpa[ia]=B_pred[iref[ia]];
    r=Corr_coeff(strdiff, Bpa, Na, NULL, NULL);
    fprintf(file_out,
	    "# Correlation predicted B_factors confchange: %.3f\n", r);
    r=Normalized_allostery(strdiff, Allostery, Bpa, Na);
    fprintf(file_out,
	    "# Correlation Allostery/B_pred confchange/B_pred: %.3f\n", r);
    free(Bpa);

  }
  fprintf(file_out, "# mode\tmean_fluct\tcorr(allostery)\n");
  for(i=0; i<NMODE; i++){
    fprintf(file_out, "%d\t%.4f\t%.3f\n", i, fluctuation[i], cc[i]);
  }
  fclose(file_out);

  // Read active sites
  char file_site[200];
  int NMAX=80, SMAX=50, N_sites=0;
  int *nres_site=malloc(NMAX*sizeof(int));
  char **chain_res=malloc(NMAX*sizeof(char *));
  for(i=0; i<NMAX; i++)chain_res[i]=malloc(SMAX*sizeof(char));
  int ***site_res=malloc(1*sizeof(int **));
  site_res[0]=Allocate_mat2_i(NMAX, SMAX);
  if(SITES[0]!='\0'){
    N_sites=Read_sites_file(site_res[0], chain_res, nres_site,
			    SITES, chain, seq, Nres, NMAX, SMAX);
    if(N_sites)strcpy(file_site, SITES);
  }
  if(N_sites==0){
    N_sites=Read_sites_PDB(site_res[0], chain_res, nres_site,
			   pdb, chain, seq, Nres, NMAX, SMAX);
    if(N_sites)strcpy(file_site, pdb);
    // Print site
    if(PRINT_SITE && N_sites)
      Print_sites(site_res[0], nres_site, N_sites, seq, pdb, chain);
  }


  // Testing active sites
  if(N_sites){
    int numran=100, k;
    int ***site_res_ran=malloc(numran*sizeof(int **));
    for(k=0; k<numran; k++){
      site_res_ran[k]=malloc(N_sites*sizeof(int *));
      for(i=0; i<N_sites; i++){
	site_res_ran[k][i]=Random_sites(nres_site[i], site_res[0][i], Nres);
      }
    }

    float prof_1, prof_2, coup_1, coup_2;
    Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		 Allosteric_profile, Allosteric_coupling,
		 Na, site_res_ran, nres_site, N_sites,
		 seq, atomres, pdbres, amm, ch, nameout1, file_site,
		 "random_site", "allosteric_profile", numran);
    Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		 Allosteric_profile, Allosteric_coupling,
		 Na, site_res, nres_site, N_sites,
		 seq, atomres, pdbres, amm, ch, nameout1, file_site,
		 "active_site", "allosteric_profile", 1);

    Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		 prof_dir, Corr_dir, Na, site_res_ran, nres_site, N_sites,
		 seq, atomres, pdbres, amm, ch, nameout1, file_site,
		 "random_site", "directionality_EC", numran);
    Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		 prof_dir, Corr_dir, Na, site_res, nres_site, N_sites,
		 seq, atomres, pdbres, amm, ch, nameout1, file_site,
		 "active_site", "directionality_EC", 1);

   Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		prof_coord, Corr_coord, Na, site_res_ran, nres_site, N_sites,
		seq, atomres, pdbres, amm, ch, nameout1, file_site,
		"random_site", "coordination_EC", numran);
    Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		 prof_coord, Corr_coord, Na, site_res, nres_site, N_sites,
		 seq, atomres, pdbres, amm, ch, nameout1, file_site,
		 "active_site", "coordination_EC", 1);

    double **d_inv_ij=Negexp_distance(atoms, iref, Na, 4);
    float *prof_d=PE_profile(d_inv_ij, Na, 0, 0);
    Test_profile(&prof_1, &prof_2, &coup_1, &coup_2, prof_d, d_inv_ij,
		 Na, site_res_ran, nres_site, N_sites,
		 seq, atomres, pdbres, amm, ch, nameout1, file_site,
		 "random_site", "negexp_dist_profile", numran);
    Test_profile(&prof_1, &prof_2, &coup_1, &coup_2, prof_d, d_inv_ij,
		 Na, site_res, nres_site, N_sites,
		 seq, atomres, pdbres, amm, ch, nameout1, file_site,
		 "active_site", "negexp_dist_profile", 1);
    Empty_matrix_d(d_inv_ij, Na); free(prof_d);
    
    if(DIR_AVE){
      Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		   prof_dir_ave, Corr_dir, Na, site_res, nres_site, N_sites,
		   seq, atomres, pdbres, amm, ch, nameout1, file_site,
		   "active_site", "directionality_ave", 1);
      Test_profile(&prof_1, &prof_2, &coup_1, &coup_2,
		   prof_dir_ave, Corr_dir, Na, site_res_ran, nres_site, N_sites,
		   seq, atomres, pdbres, amm, ch, nameout1, file_site,
		   "random_site", "directionality_ave", numran);
    }
    for(k=0; k<numran; k++)Empty_matrix_i(site_res_ran[k], N_sites);
    free(site_res_ran);
  }



  // Empty
  free(cc); free(dev); if(strdiff)free(strdiff);
  free(fluctuation); free(weight);
  free(Allostery);
  free(nres_site);
  Empty_matrix_i(site_res[0], NMAX);
  free(site_res);
  free(amm); free(ch);
  for(i=0; i<Na; i++)free(pdbres[i]); free(pdbres);

  for(i=0; i<Na; i++)Empty_matrix_f(r0ijk[i], Na); free(r0ijk);
  free(iref); free(atomres); free(nc);
  free(abs_dr); free(abs_dr2); free(strain);
  Empty_matrix_d(Allosteric_coupling, Na);
  Empty_matrix_f(dr, Na);
  Empty_matrix_i(clist, Na);
  if(Corr_dr) {Empty_matrix_d(Corr_dr, Na); free(prof_dr);}
  if(Corr_dir){Empty_matrix_d(Corr_dir, Na); free(prof_dir);}
  if(Corr_dr_dir){Empty_matrix_d(Corr_dr_dir, Na); free(prof_dr_dir);}
  if(Corr_Bahar) {Empty_matrix_d(Corr_Bahar, Na);  free(prof_Bahar);}
  if(Corr_coord) {Empty_matrix_d(Corr_coord, Na); free(prof_coord);}
  if(Corr_str)   {Empty_matrix_d(Corr_str, Na); free(prof_str);}
  if(Corr_str_dir){Empty_matrix_d(Corr_str_dir, Na); free(prof_str_dir);}
  //Empty_matrix_d(Allostery, Na); free(prof_Allostery);
  free(Allosteric_profile);
  //free(Allosteric_ave);
  free(Allosteric_profile_weight);
  free(Allosteric_distance_weight);
}

float *PE_profile(double **Corr_in, int N, int IJMIN, int invert)
{
  float E_min=0.0000; int i, j;
  double **Corr=malloc(N*sizeof(double *)), Cmin=0.00001;
  for(i=0; i<N; i++){
    Corr[i]=malloc(N*sizeof(double));
    if(invert){
      for(j=0; j<N; j++){
	if(fabs(i-j)>=IJMIN){
	  if(Corr_in[i][j]>Cmin){
	    Corr[i][j]=1./Corr_in[i][j];
	  }else{
	    Corr[i][j]=1./Cmin;
	  }
	}else{
	  Corr[i][j]=0;
	}
      }
    }else{
      for(j=0; j<N; j++){
	if(fabs(i-j)>=IJMIN){Corr[i][j]=Corr_in[i][j];}
	else{Corr[i][j]=0;}
      }
    }
  }

  float *eigenval=malloc(N*sizeof(float));
  float **eigenvec=malloc(N*sizeof(float *));
  for(i=0; i<N; i++)eigenvec[i]=malloc(N*sizeof(float));
  d_Diagonalize(N, Corr, eigenval, eigenvec, 1, E_min);
  float *PE=eigenvec[0];
  double sum=0, norm=0; int num=0;
  for(i=0; i<N; i++){
    sum+=PE[i];
    for(j=0; j<=i; j++){
      if((i-j)>=IJMIN){norm+=Corr_in[i][j]; num++;}
    }
  }
  norm=(N*norm)/(sum*num);
  float *prof=malloc(N*sizeof(float));
  for(i=0; i<N; i++)prof[i]=norm*PE[i];
  //for(i=0; i<N; i++)prof[i]=eigenval[0]*PE[i]*PE[i];
  Empty_matrix_f(eigenvec, N); free(eigenval);
  Empty_matrix_d(Corr, N);
  return(prof);
}

float Eigenvalue_main(double **Corr, int Na)
{
  float E_min=0.00000000001; int i;
  float *eigenval=malloc(Na*sizeof(double));
  float **eigenvec=malloc(Na*sizeof(float *));
  for(i=0; i<Na; i++)eigenvec[i]=malloc(Na*sizeof(float));
  d_Diagonalize(Na, Corr, eigenval, eigenvec, 1, E_min);
  float Ev=eigenval[0];
  Empty_matrix_f(eigenvec, Na); free(eigenval);
  return(Ev);
}

float Singvalue_main(double **A, int Na, int Nb)
{
  float E_min=0.00; int i, j, k, N1, N2;
  if(Na<Nb){N1=Na; N2=Nb;}else{N1=Nb; N2=Na;}
  float *eigenval=malloc(N1*sizeof(float));
  float **eigenvec=Allocate_mat2_f(N1, N1);
  double **Mat=Allocate_mat2_d(N1, N1);
  for(i=0; i<N1; i++){
    for(j=0; j<N1; j++){
      for(k=0; k<N2; k++)Mat[i][j]+=A[i][k]*A[j][k];
    }
  }

  d_Diagonalize(N1, Mat, eigenval, eigenvec, 1, E_min);
  float Ev=sqrt(eigenval[0]);
  Empty_matrix_d(Mat, Na);
  Empty_matrix_f(eigenvec, Na);
  //for(i=0; i<N1; i++)printf("%.3g ", eigenval[i]); printf("\n");
  free(eigenval);
  return(Ev);
}


void Normalize_profile(double *v, int N){
  double sum=0; int i;
  for(i=0; i<N; i++)sum+=v[i]*v[i]; sum=sqrt(sum/N);
  for(i=0; i<N; i++)v[i]/=sum;
}


int Read_sites_PDB(int **site_res, char **chain_res, int *nres_site,
		   char *pdb, char *chain, struct residue *seq,
		   int Nres, int NMAX, int SMAX)
{
  int Compression=0;
  printf("Reading sites in file %s, chain %s %d residues\n",pdb, chain, Nres);
  FILE *file_in=Open_compressed_file(pdb, &Compression);
  if(file_in==NULL){
    printf("WARNING, file %s not found\n", pdb); return(0);
  }
  int nsites=-1, sres=0, l, i;
  for(i=0; i<NMAX; i++)nres_site[i]=0;
  char string[1000], resnum[6], resnam[5];
  char snam[5], snam_old[5]="   ", ch; i=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "SITE", 4)==0){
      sscanf(string+11, "%s", snam);
      if(strncmp(snam, snam_old, 3)!=0){
	if(nsites>=0)nres_site[nsites]=sres;
	strcpy(snam_old, snam); sres=0; nsites++;
	if(nsites >= NMAX){
	  printf("ERROR, too many sites > %d\n", NMAX); nsites--; break;
	}
      }
      for(l=18; l<60; l+=11){
	i=Read_site(resnum, resnam, &ch, string+l, chain);
	if(i==-2){break;}else if(i==-1){continue;}
	i=Find_residue(resnum, ch, seq, Nres, i+1);
	if(i<0){
	  printf("WARNING, residue %s %s not found\n", resnum, resnam);
	  continue;
	}
	if(seq[i].amm!=Code_3_1(resnam)){
	  printf("WARNING, residue %d is %c while it is %s in site\n",
		 i, seq[i].amm, resnam);
	}
	site_res[nsites][sres]=i;
	chain_res[nsites][sres]=ch;
	sres++;
	if(sres > SMAX){
	  printf("ERROR, too many residues > %d\n", SMAX); goto end;
	}
      }
    }else{
      if((nsites>=0)||(strncmp(string, "ATOM", 4)==0)){break;}
      else{continue;}
    }
  }
 end:
  fclose(file_in);
  if(Compression)Delete_tmp_file();
  if(nsites>=0)nres_site[nsites]=sres;
  nsites++;
  int num=0, msites=0;
  for(i=0; i<nsites; i++){
    num+=nres_site[i];
    if(nres_site[i]>0)msites++;
  }
  printf("%d active residues read in %d active sites out of %d\n",
	 num, msites, nsites);
  return(nsites);
}

int Read_sites_file(int **site_res, char **chain_res, int *nres_site,
		    char *file, char *chain, struct residue *seq,
		    int Nres, int NMAX, int SMAX)
{
  printf("Reading sites in file %s, chain %s %d residues\n",file, chain, Nres);
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("WARNING, active site file %s not found\n", file); return(0);
  }
  int nsites=-1, sres=0, i, res;
  for(i=0; i<NMAX; i++)nres_site[i]=0;
  char string[1000], resnum[6]="", resnam[5], icode=' ';
  char snam[5], snam_old[5]="   ", ch[6]; i=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      continue;
    }else{
      sscanf(string, "%s%s%s%d%c", snam, resnam, ch, &res, &icode);
      if(strcmp(snam, snam_old)!=0){
	if(nsites>=0)nres_site[nsites]=sres;
	sres=0; nsites++;
	if(nsites >= NMAX){
	  printf("ERROR, too many sites > %d\n", NMAX); break;
	}
	strcpy(snam_old, snam);
      }
      if(icode=='\n')icode=' ';
      sprintf(resnum, "%4d%c", res, icode);
      i=Find_residue(resnum, ch[0], seq, Nres, i+1);
      if(i<0){
	printf("Residue %s %s %c discarded, chain: %s\n",
	       resnam, resnum,ch[0],chain); continue;
      }
      if(seq[i].amm!=Code_3_1(resnam)){
	printf("WARNING, residue %d is %c while it is %s in site %s\n",
	       i, seq[i].amm, resnam, snam);
      }
      site_res[nsites][sres]=i;
      chain_res[nsites][sres]=ch[0];
      sres++;
      if(sres > SMAX){
	printf("ERROR, too many residues > %d\n", SMAX); goto end;
      }
    }
  }
 end:
  fclose(file_in);
  if(nsites>=0)nres_site[nsites]=sres;
  nsites++;
  int num=0, msites=0;
  for(i=0; i<nsites; i++){
    num+=nres_site[i];
    if(nres_site[i]>0)msites++;
  }
  printf("file %s, reading %d active residues in %d sites out of %d\n",
	 file, num, msites, nsites);
  return(nsites);
}

int Read_site(char *resnum, char *resnam, char *reschain,
	      char *string, char *chain){
  if(string[0]==' ')return(-2);
  char *c=chain; int inum;
  sscanf(string, "%s", resnam);
  *reschain=string[4];
  sscanf(string+5, "%d", &inum);
  sprintf(resnum, "%4d%c", inum, string[9]);
  if(*c=='*')return(0);
  while(*c!='\0'){if(*c==*reschain)return(0); c++;}
  if(VERBOSE)
    printf("Residue %s %s %c discarded, chain: %s\n",
	   resnam,resnum,*reschain,chain);
  return(-1);
}

int Find_residue(char *resnum, char chain,
		 struct residue *seq, int Nres, int i_ini)
{
  int i=0; if(i_ini>0){i=i_ini;}else{i=0;}
  while(i<Nres){
    if((strcmp(seq[i].pdbres, resnum)==0)&&(seq[i].chain==chain))return(i);
    i++;
  }
  i=0;
  while(i<i_ini){
    if((strcmp(seq[i].pdbres, resnum)==0)&&(seq[i].chain==chain))return(i);
    i++;
  }
  return(-1);
}

int Print_links(double **Coupling, int N,
		char **pdbres, char *amm, char *chain,
		char *nameout1, char *ext, int sign)
{
  char nameout[200]; int i, j;
  sprintf(nameout, "%s%s", nameout1, ext);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);

  // Mean and standard deviation
  double sum1=0, sum2=0; int norm=0;
  for(i=0; i<N; i++){
    double *c=Coupling[i];
    for(j=i+2; j<N; j++){
      sum1+=c[j]; sum2+=c[j]*c[j]; norm++;
    }
  }
  sum1/=norm; sum2=sqrt((sum2-norm*sum1*sum1)/(norm-1));
  double thr=sum1;
  if(sign>0){
    thr+=SIGMA*sum2;
    fprintf(file_out, "# mean= %.3f s.d.= %.3f thr=mean+%.2f*s.d.= %.3f\n",
	    sum1, sum2, SIGMA, thr);
  }else if(sign<0){
    thr-=SIGMA*sum2;
    fprintf(file_out, "# mean= %.3f s.d.= %.3f thr=mean-%.2f*s.d.= %.3f\n",
	    sum1, sum2, SIGMA, thr);
  }
  for(i=0; i<N; i++){
    double *c=Coupling[i];
    for(j=i+1; j<N; j++){
      if((sign==0)||
	 ((sign>0)&&(c[j]>=thr)&&(c[j]>=0))||
	 ((sign<0)&&(c[j]<=thr))){
	fprintf(file_out, "%c\t%s\t%c\t%c\t%s\t%c\t%.3f\n",
		amm[i], pdbres[i], chain[i],
		amm[j], pdbres[j], chain[j], c[j]);
      }
    }
  }
  fclose(file_out);
  return(0);
}

int Get_zeta(double **Coupling, int N)
{
  // Mean and standard deviation depending on l=|i-j|
  // Compute z=(c_{i,i+l}-<c>_l)/sigma_l
  // If |z| > 1, amplify c -> c*z
  int i, j, l, ll=-1, num_l=0;;
  float *c1=malloc(N*sizeof(double));
  float *c2=malloc(N*sizeof(double));
  double sum1_l=0, sum2_l=0;
  for(l=0; l<N; l++){
    double sum1=0, sum2=0; int num=N-l;
    for(i=0; i<(N-l); i++){
      float c=Coupling[i][i+l];
      sum1+=c; sum2+=c*c;
    }
    if(l==0){
      c1[l]=sum1/num;
      c2[l]=sqrt((sum2-num*c1[l]*c1[l])/(num-1));
    }else if((ll<0) && (l<N-9) && (fabs(c1[l]-c1[l-1]) > 3*c2[l-1]/sqrt(num))){
      c1[l]=sum1/num;
      c2[l]=sqrt((sum2-sum1*sum1/num)/(num-1));
    }else{
      if(ll<0)ll=l; sum1_l+=sum1; sum2_l+=sum2; num_l+=num;
    }
  }
  sum1_l/=num_l;
  sum2_l=sqrt((sum2_l-num_l*sum1_l*sum1_l)/(num_l-1));
  for(l=ll; l<N; l++){
    c1[l]=sum1_l;
    c2[l]=sum2_l;
  }

  for(i=0; i<N; i++){
    double *c=Coupling[i];
    for(j=i+1; j<N; j++){
      float z=(c[j]-c1[j-i])/c2[j-i];
      if(z<0)z=-z;
      if(z>1){c[j]*=z; Coupling[j][i]=c[j];}
    }
  }
  free(c1); free(c2);
  return(0);
}


int Test_profile(float *prof_1, float *prof_2, float *pair_1, float *pair_2,
		 float *profile, double **coupling, int Na,
		 int ***site_res, int *nres_site, int N_sites,
		 struct residue *seq, int *atomres,
		 char **pdbres, char *amm, char *chain,
		 char *nameout1, char *filename,
		 char *sitename, char *profname,
		 int num)
{
  char nameout[200];
  sprintf(nameout, "%s_%s_%s_summary.dat", nameout1, sitename, profname);
  FILE *file_out=NULL;
  if(num==1){
    file_out=fopen(nameout, "w");
    printf("Writing %s\n", nameout);
    fprintf(file_out, "# read active site in file %s\n", filename);
    fprintf(file_out, "# profile %s\n", profname);
    if(num >1){
      fprintf(file_out, "# site\tnres\t<prof>\ts.dev.\t<coup>\ts.dev.\n");
    }else{
      fprintf(file_out, "# site\tnres\t<prof>\tZ_score\t<coup>\tZ_score\n");
    }
  }
  int *atom_site=malloc(Na*sizeof(int)), i, i1, i2, j1, j2;
  for(i=0; i<Na; i++)atom_site[i]=0;

  int isite;
  for(isite=0; isite<N_sites; isite++){
    if(nres_site[isite]<=0)continue;
    //fprintf(file_out, "Site %d\n", isite+1);
    int nr=nres_site[isite], kran;
    double tot_pair1=0, tot_prof1=0;
    double tot_pair2=0, tot_prof2=0;
    // Sum over num random sites 
    for(kran=0; kran<num; kran++){
      double sum_prof=0, sum_pair=0; int npair=0;
      int *site=site_res[kran][isite];
      for(i1=0; i1<nr; i1++){
	j1=atomres[site[i1]];
	if((j1<0)||(j1>=Na)){
	  int r=site[i1];
	  printf("WARNING, reference atom of res. %d  %c%s%c does not exist",
		 r, seq[r].amm, seq[r].pdbres, seq[r].chain);
	  printf(" (%d dis in str.2?)\n", j1);
	  nr--; continue;
	}
	//struct residue *r=seq+site_res[k][i1];
	//fprintf(file_out, "%c %s %c\n", r->amm, r->pdbres, r->chain);
	atom_site[j1]=1;
	sum_prof+=profile[j1];
	for(i2=0; i2<i1; i2++){
	  j2=atomres[site[i2]];
	  if((j2<0)||(j2>=Na))continue;
	  sum_pair+=coupling[j1][j2];
	  npair++;
	}
      } // end of the random site kran
      sum_prof/=nr;
      tot_prof1+=sum_prof;
      if(num > 1)tot_prof2+=sum_prof*sum_prof;
      if(npair){
	sum_pair/=npair;
	tot_pair1+=sum_pair;
	if(num > 1)tot_pair2+=sum_pair*sum_pair;
      }
    } // end of the active site or all the random sites
    if(num > 1){
      tot_prof1/=num;
      tot_pair1/=num;
      tot_prof2=(tot_prof2-num*tot_prof1*tot_prof1);
      tot_prof2=sqrt(tot_prof2/(num-1));
      tot_pair2=(tot_pair2-num*tot_pair1*tot_pair1);
      tot_pair2=sqrt(tot_pair2/(num-1));
      *prof_1=tot_prof1; *prof_2=tot_prof2;
      *pair_1=tot_pair1; *pair_2=tot_pair2;
      //fprintf(file_out, "%d\t%2d\t%.3f\t%.2g\t%.3f\t%.2g\n",
      //      k+1,nr,tot_prof1,tot_prof2,tot_pair1,tot_pair2);
    }else{
      float Z_prof=(tot_prof1-*prof_1)/(*prof_2);
      float Z_pair=(tot_pair1-*pair_1)/(*pair_2);
      fprintf(file_out, "%d\t%2d\t%.3f\t%.2g\t%.3f\t%.2g\n",
	      isite+1,nr,tot_prof1,Z_prof,tot_pair1,Z_pair);
    }
  }

  // Profile for active residues
  if(num==1){
    fclose(file_out);
    sprintf(nameout, "%s_%s_%s.dat", nameout1, sitename, profname);
    file_out=fopen(nameout, "w");
    printf("Writing %s\n", nameout);
    fprintf(file_out, "#res\taa\tchain\tprofile\n");
    for(i=0; i<Na; i++){
      if(atom_site[i]){
	fprintf(file_out, "%3s\t%c\t%c\t%.3f\n",
		pdbres[i], amm[i], chain[i], profile[i]);
      }
    }

    fclose(file_out);
  }

  free(atom_site);
  return(0);
}


int *Random_sites(int nres, int *site_res, int N)
{
  int *site_ran=malloc(nres*sizeof(int)), i;
  //float RANFACTOR=pow(12.0,1/3.0);
  unsigned long iran=randomgenerator();
  if(ini_ran==0){
    InitRandom((RANDOMTYPE)iran); ini_ran=1;
  }
  int imin=site_res[0], imax=imin;
  for(i=1; i<nres; i++){
    if(site_res[i]<imin){imin=site_res[i];}
    else if(site_res[i]>imax){imax=site_res[i];}
  }
  int shiftmin=-imin, shiftmax=N-1-imax, shift, thr=8; i=0;
  while(i<30){
    shift=RandomFloating()*(shiftmax-shiftmin)+shiftmin;
    if((shift>thr)||(shift<-thr))break; i++;
  }
  if(i<30){
    for(i=0; i<nres; i++)site_ran[i]=site_res[i]+shift;
  }else{
    for(i=0; i<nres; i++){
      site_ran[i]=RandomFloating()*N;
    }
  }
  return(site_ran);
}


double *Cont_ave(float *v, int N, int **clist, int *nc)
{
  double *Ave=malloc(N*sizeof(double)); int i, j;
  for(i=0; i<N; i++){
    Ave[i]=v[i];
    for(j=0; j<nc[i]; j++)Ave[i]+=v[clist[i][j]];
    Ave[i]/=(nc[i]+1);
  }
  return(Ave);
}

int Print_matrix(double **Coupling, int N, char **pdbres, char *amm,
		 char *nameout1, char *ext)
{
  int IJ_MIN=0;
  char nameout[200]; int i, j;
  sprintf(nameout, "%s%s", nameout1, ext);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);

  // Mean and standard deviation
  double sum1=0, sum2=0; int norm=0;
  for(i=0; i<N; i++){
    double *c=Coupling[i];
    for(j=i+2; j<N; j++){
      sum1+=c[j]; sum2+=c[j]*c[j]; norm++;
    }
  }
  sum1/=norm; sum2=sqrt((sum2-norm*sum1*sum1)/(norm-1));
  double thrmax=sum1+SIGMA*sum2, thrmin=sum1-SIGMA*sum2;
  //fprintf(file_out, "# mean= %.3f s.d.= %.3f thr=mean+%.2f*s.d.= %.3f\n",
  //	  sum1, sum2, SIGMA, thrmax);
  int *signi=malloc(N*sizeof(int));
  for(i=0; i<N; i++)signi[i]=0;
  for(i=0; i<N; i++){
    if(signi[i])continue;
    double *c=Coupling[i]; int s=0;
    for(j=0; j<N; j++){
      if(((c[j]>=thrmax)||(c[j]<=thrmin))&&(abs(i-j)>IJ_MIN)){s=1; break;}
    }
    if(s){signi[i]=1; if(j<N)signi[j]=1;}
  }
  // Write matrix
  char **nameres=malloc(N*sizeof(char *));
  fprintf(file_out, "Residues");
  for(i=0; i<N; i++){
    char *name=malloc(7*sizeof(char)); nameres[i]=name;
    name[0]=amm[i];
    int k1=0, k2=1;
    for(k1=0; k1<6; k1++){
      if(pdbres[i][k1]!=' '){name[k2]=pdbres[i][k1]; k2++;}
    }
    while(k2<7){name[k2]='\0'; k2++;}
    if(signi[i])fprintf(file_out, "\t%s", nameres[i]);
  }
  fprintf(file_out, "\n");
  for(i=0; i<N; i++){
    if(signi[i]==0)continue;
    double *c=Coupling[i];
    fprintf(file_out, "%s", nameres[i]);
    for(j=0; j<N; j++)if(signi[j])fprintf(file_out, "\t%.3f",c[j]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
  for(i=0; i<N; i++)free(nameres[i]); free(nameres);
  free(signi);
  return(0);
}


float *Matrix_average(double **Mat, int N, int IJMIN, char *name)
{
  float *mean=malloc(N*sizeof(float));
  int i, j, num;
  if(IJMIN){num=N-2*IJMIN+1;}
  else{num=N;}
  for(i=0; i<N; i++){
    mean[i]=0; double *m=Mat[i];
    if(IJMIN){
      for(j=0; j<=(i-IJMIN); j++)mean[i]+=m[j];
      for(j=i+IJMIN; j<N; j++)mean[i]+=m[j];
    }else{
      for(j=0; j<N; j++)mean[i]+=m[j];
    }
  }
  for(i=0; i<N; i++)mean[i]/=num;
  return(mean);
}

int Extract_atoms(int *iref, int *atomres, struct residue *seq, float *mass,
		  atom *atoms, struct Reference Ref, char *ANAME)
{
  int ia, Na=0;
  for(ia=0; ia<Ref.N_ref; ia++)atomres[ia]=-1;
  for(ia=0; ia<Ref.N_ref; ia++){
    atom *atom_i=atoms+Ref.atom_num[ia];
    if((strncmp(atom_i->name, ANAME, 2)==0)||
       ((seq[atom_i->res].amm=='G')&&(strncmp(atom_i->name, "CA", 2)==0))){
      if(*mass==0)*mass=Ref.mass_atom[ia];
      if(atom_i->res >= Ref.N_ref){
	printf("ERROR, too many residues in extract atoms\n"); exit(8);
      }
      iref[Na]=ia;
      atomres[atom_i->res]=Na; Na++;
    }
  }
  /*for(i=0; i<Na; i++){
    int r=(atoms+Ref.atom_num[iref[i]])->res;
    printf("%d %c %s\n", i, seq[r].amm, seq[r].pdbres);
    }*/
  return(Na);
}

void Get_interaction_list(int **clist, int *nc,
			  atom *atoms, int *atomres, int Na,
			  struct interaction *Int_list, int N_int)
{
  int ia, ja, n;
  for(ia=0; ia<Na; ia++)nc[ia]=0;
  for(n=0; n<N_int; n++){
    ia=atomres[(atoms+Int_list[n].i1)->res];
    ja=atomres[(atoms+Int_list[n].i2)->res];
    if((ia<0)||(ja<0))continue;
    nc[ia]++; nc[ja]++;
  }
  for(ia=0; ia<Na; ia++){
    clist[ia]=malloc((nc[ia])*sizeof(int));
    nc[ia]=0;
  }
  for(n=0; n<N_int; n++){
    ia=atomres[(atoms+Int_list[n].i1)->res];
    ja=atomres[(atoms+Int_list[n].i2)->res];
    if((ia<0)||(ja<0))continue;
    clist[ia][nc[ia]]=ja;
    clist[ja][nc[ja]]=ia; nc[ia]++;
    clist[ja][nc[ja]]=ia; nc[ja]++;
  }
}

void Pairwise_distances(float ***r0ijk, atom *atoms, int *iref, int Na)
{
  int ia, ja, k;
  for(ia=0; ia<Na; ia++){
    float *ri=atoms[iref[ia]].r;
    float **rij=Allocate_mat2_f(Na,3);
    r0ijk[ia]=rij;
    for(ja=0; ja<ia; ja++){
      float *rj=atoms[iref[ja]].r, *d=rij[ja];
      double sum=0;
      for(k=0; k<3; k++){
	*d=ri[k]-(*rj); sum+=(*d)*(*d); d++; rj++;
      }
      sum=sqrt(sum);
      d=rij[ja]; for(k=0; k<3; k++){(*d)/=sum; d++;}
    }
  }
}


double **Negexp_distance(atom *atoms, int *iref, int Na, float D)
{
  double **eij=Allocate_mat2_d(Na, Na);
  int i, j, k;
  for(i=0; i<Na; i++){
    float *ri=atoms[iref[i]].r;
    eij[i][i]=1;
    for(j=0; j<i; j++){
      float *rj=atoms[iref[j]].r;
      double sum=0;
      for(k=0; k<3; k++){
	float d=ri[k]-(*rj); sum+=d*d; rj++;
      }
      sum=exp(-sqrt(sum)/D);
      eij[i][j]=sum; eij[j][i]=sum;
    }
  }
  return(eij);
}

void Fill_low_diag(double **Corr, int Na){
  int ia, ja;
  for(ia=0; ia<Na; ia++){
    for(ja=0; ja<ia; ja++)Corr[ja][ia]=Corr[ia][ja];
  }
}

void Allosteric_propagation(double **Allosteric_coupling, int Na,
			    char **pdbres, char *amm, char *chain,
			    char *nameout)
{
  char namefile[200];
  sprintf(namefile, "%s_perturbation_propagation.dat", nameout);
  FILE *file_out=fopen(namefile, "w");
  printf("Writing %s\n", namefile);
  fprintf(file_out, "#AA res chain Allo prop_length C.c.\n");
  float *d=malloc(Na*sizeof(float));
  float *y=malloc(Na*sizeof(float));
  float slope, offset, r; int i, j;
  for(i=0; i<Na; i++){
    double *Ci=Allosteric_coupling[i];
    for(j=0; j<Na; j++){
      if(i>j){d[j]=i-j;}else{d[j]=j-i;}
      y[j]=log(Ci[j]);
    }
    r=Corr_coeff(d, y, Na, &slope, &offset);
    fprintf(file_out, "%c\t%s\t%c\t%.3f\t%.3f\t%.3f\n",
	    amm[i], pdbres[i], chain[i], exp(offset), -1/slope, r);
  }
  fclose(file_out);
  free(d); free(y);
}

void Print_sites(int **site_res, int *nres_site, int nsites,
		 struct residue *seq, char *pdb, char *chain)
{
  // Printing
  char pdbid[80], nameout[200], amm[5]; int k, j;
  GetPdbId(pdb,pdbid);
  sprintf(nameout, "%s%s_sites.in", pdbid, chain);
  printf("Writing active sites in %s\n", nameout);
  FILE *file_out=fopen(nameout, "w");
  for(k=0; k<nsites; k++){
    for(j=0; j<nres_site[k]; j++){
      struct residue *res=seq+site_res[k][j];
      Name3(amm, res->i_aa);
      fprintf(file_out, "%d\t%s\t%c\t%s\n",
	      k+1, amm, res->chain, res->pdbres);
    }
  }	  
  fclose(file_out);
}

float Normalized_allostery(float *strdiff, float *Allostery, float *Bpa, int Na)
{
  float *x=malloc(Na*sizeof(float));
  float *y=malloc(Na*sizeof(float));
  int i;
  for(i=0; i<Na; i++){
    x[i]=Allostery[i]/Bpa[i];
    y[i]=strdiff[i]/Bpa[i];
  }
  float r=Corr_coeff(x, y, Na, NULL, NULL);
  free(x); free(y);
  return(r);
}

#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "allostery.h"
#include "read.h"
#include "vector.h"
#include "contacts.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "allocate.h"
#include "diagonalize.h"
#include "interactions_tnm.h"
#include "random3.h"
#include "EC.h"
#include "binding_sites.h"

#include <time.h>

double SCALE_COORD=2; // Coordination_ij=exp(-SCALE_COORD*<(dij-<dij>)^2)
double SCALE_D=2;      // Proximity_ij=exp(-d_ij/SCALE_D);

int PRINT_RAN=0; // Print results for random sites?
int PRINT_DIFF=1;  // Print difference in Conformation per residue?
int PRINT_SITE=1;  // Print active site read in the PDB?
#define VERBOSE 0
int ini_ran=0;

struct site_results{
  float prof, P_p, Z_p; // profile
  float coup, P_c, Z_c;  // coupling
  float p_ran1, p_ran2;  // profile_random
  float c_ran1, c_ran2;  // coupling random
};

float *Make_profile(float **Mat, int N, char PROF_TYPE, char *what);
double **Get_matrix_positive(float **Mat, int N);
float PE_profile(float *PE, double **Mat, int N);
float Eigenvalue_main(double **Corr, int Na);
float Singvalue_main(double **A, int Na, int Nb);
void Normalize_profile(double *v, int N);
int Print_links(float **Coupling, int N, char *nameout, char *ext, int sign);
void Print_names(int N, char **pdbres, char *amm, char *chain, char *nameout);

extern char DOF_LABEL[];
extern double G_NM_holo;
extern double DG_NM;
extern double DG_NM_internal;
extern double DG_NM_rigid;
extern double E_cont_bind;

// Binding sites
void Binding_sites(int *site,
		   float *Prof_dir, float **Coup_dir,
		   float *Prof_coord, float **Coup_coord,
		   float *Prof_proximity, float **Coup_proximity, 
		   float *Prof_deformation, float **Coup_deformation,
		   //float *Prof_deformation_dist,float **Coup_deformation_dist,
		   float *Prof_dr_dir, float **Coup_dr_dir,
		   int *atomres, char **pdbres, char *amm, 
		   int Na, char *chain, int Nchain, struct residue *seq,
		   int Nres, char *pdb, char *nameout, char *SITES);
/*
int Read_sites_PDB(int **site_res, char **chain_res,
		   int *nres_site, char **description,
		   char *pdb, char *chain, struct residue *seq,
		   int Nres, int NMAX, int SMAX);
int Read_sites_file(int **site_res, char **chain_res,
		    int *nres_site, char **description,
		    char *pdb, char *chain, struct residue *seq,
		    int Nres, int NMAX, int SMAX); 
int Site_overlap(int *isite, int n1, int *jsite, int n2);
*/
void Random_sites(int *site_ran, int nres, int *site_res,
		  int *ichain, int *min_chain, int *max_chain,
		  int *all_chain, int nchain,
		  int *ini_chain, int *end_chain);
int Count_chains(int *ichain, int *min_chain, int *max_chain, char *chain_lab,
		 char *chain_res, int *site_res, int nr);
void Identify_chains(int *all_chain, char *chain_lab, int nchain,
		     char *lab_chain, int Nchain);

unsigned long randomgenerator(void);
int Read_site(char *resnum, char *resnam, char *reschain,
	      char *string, char *chain);
int Find_residue(char *resnum, char chain,
		 struct residue *seq, int Nres, int i_ini);
float Ave_site(double *sigprof, double *ave_pair,
	       double *sig_pair, int *n_notfound,
	       int *site, int *atomres, int nr, int *atom_site,
	       float *profile, float **coupling, int Na);
int Test_site(struct site_results *site_results, int *atom_site,
	      float *profile, float **coupling, int Na,
	      int *site_res, int **site_res_ran, int numran,
	      int nres, int *atomres);
int Test_pairs(struct site_results *site_results, int *atom_site,
	       float *profile, float **coupling, int Na,
	       int *coupled_site, int **coupled_site_ran, int numran,
	       int ni, int nij, int *atomres);
float Pair_site(double *sigpair, int *n_notfound,
		int *site, int n1, int n2,
		int *atomres, int *atom_site,
		float **coupling, int Na);

int Extract_atoms(int *iref, int *iatom, int *atomres,
		  struct residue *seq, float *mass,
		  atom *atoms, struct Reference Ref, char *ANAME);
void Pairwise_distances(float **d02ij, float ***r0ijk,
			atom *atoms, int *iatom, int Na);
void Fill_low_diag(float **Corr, int Na);
void Deformation_propagation(float **Deformation_coupling, int Na,
			    char **pdbres, char *amm, char *chain,
			    char *nameout);
void Print_sites(int **site_res, int *nres_site,
		 char **description, int nsites, 
		 struct residue *seq, char *pdb, char *chain);
int Get_zeta(float **Coupling, int N);
float **Negexp_distance(float **d02ij, int Na, float D);
void Ave_St_dev(double *sum1, double *sum2, int n);
void Coupling_chains(double **Chain_dir, double **Chain_coord,
		     double **Chain_def, double **Chain_cov,
		     double **norm_chain, char *ch, int Na,
		     char *chain_label, int Nchain, char *nameout1);
double **Chain_coup(float **Coupling, double **norm_chain, float **Proximity,
		    int Nchain, char *ch, int Na);
void Print_interface_coupling(struct interaction *Int_list_inter,
			      int N_int_inter,
			      struct interaction *Int_list, int N_int,
			      atom *atoms, int Nchain, char *nameout,
			      float **Coup_dir, float**Coup_coord,
			      float **Coup_deformation,
			      double **Chain_dir, double **Chain_coord,
			      double **Chain_def, double **norm_chain, int Na);

/****************************************************************

      Main routine

*****************************************************************/

void Predict_allostery(struct Normal_Mode NM, atom *atoms,
		       struct Reference Ref,
		       struct interaction *Int_list, int N_int,
		       struct residue *seq, char *nameout1, //int mode,
		       float *Confchange,
		       char *pdb, char *chain,
		       int Nres, char *SITES, int anharmonic,
		       struct interaction *Int_list_inter, int N_int_inter)
{
  // Predicted fluctuations
  float *sigma2; char nameout[200];
  if(anharmonic==0){
    sigma2=NM.sigma2; sprintf(nameout, "%s", nameout1);
  }else{
    sigma2=NM.sigma2_anhar; sprintf(nameout, "%s.anharmonic", nameout1);
  }

  // Extract c_alpha
  char SEL[4]="CA"; if(strcmp(REF, "CB")==0)strcpy(SEL, "CB");
  int Nref=Ref.N_ref, iref[Nref], iatom[Nref], atomres[Nref];
  float mass=0;
  int Na=Extract_atoms(iref, iatom, atomres, seq, &mass, atoms, Ref, SEL);

  // Interaction list
  int *nc, **clist, **cnum;
  Get_contact_list(&nc, &clist, &cnum, Na, atoms, Int_list, N_int);

  // Pairwise differences
  float **d02ij=malloc(Na*sizeof(float *));
  float ***r0ijk=malloc(Na*sizeof(float **));
  Pairwise_distances(d02ij, r0ijk, atoms, iatom, Na);

  // Allocate
  float abs_dr[Na], abs_dr2[Na], strain[Na];
  float **dr=Allocate_mat2_f(Na, 3);
  float **dir=Allocate_mat2_f(Na, 3);

  ////////////////////////////////////////////////////////////
  // DIRECTIONALITY COUPLING and COORDINATION COUPLING

  float **Coup_coord=NULL;
  float **Coup_dr_dir=NULL;
  if(PRINT_COORD_COUPLING ||PRINT_SIGMA_DIJ || PRINT_COV_COUPLING){
    Coup_coord=Allocate_mat2_f(Na, Na);
    Coup_dr_dir=Allocate_mat2_f(Na, Na);
  }

  float **Coup_dir=NULL;
  if(PRINT_DIR_COUPLING)
    Coup_dir=Allocate_mat2_f(Na, Na);
  float **Coup_dr=NULL;
  if(0)Coup_dr=Allocate_mat2_f(Na, Na);
  float **Coup_Bahar=NULL;
  if(0)Coup_Bahar=Allocate_mat2_f(Na, Na);
  float **Coup_str=NULL;
  float **Coup_str_dir=NULL;
  if(STRAIN){
    Coup_str=Allocate_mat2_f(Na, Na);
    Coup_str_dir=Allocate_mat2_f(Na, Na);
  }

  int i, j, ia, ja, a, n;
  double dr2[Na], norm_w=0;
  for(ia=0; ia<Na; ia++)dr2[ia]=0;
  for(a=0; a<NM.N; a++){
    // For each normal mode
    if((NM.select[a]==0)||(sigma2[a]==0))
      continue;
    float *Mode=NM.Cart[a];
    float w=sigma2[a];
    norm_w+=w;

    // Individual atoms
    for(ia=0; ia<Na; ia++){
      int i3=3*iref[ia], j3;
      double r2=0, x;
      for(j=0; j<3; j++){x=Mode[i3+j]; dr[ia][j]=x; r2+=x*x;}
      dr2[ia]+=w*r2;
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
	double dri_drj=0;
	for(j=0; j<3; j++)dri_drj+=dr[ia][j]*dr[ja][j];
	if(Coup_dr_dir)Coup_dr_dir[ia][ja]+=w*dri_drj;       // <ri*rj>
	if(Coup_dr)Coup_dr[ia][ja]+=w*abs_dr[ia]*abs_dr[ja]; // <|ri|*|rj|>
	if(Coup_coord){
	  double sum=0;
	  for(j=0; j<3; j++)sum+=r0i[ja][j]*(dr[ia][j]-dr[ja][j]);
	  Coup_coord[ia][ja]+=w*(sum*sum);
	}
	if(Coup_Bahar)
	  Coup_Bahar[ia][ja]+=w*(abs_dr2[ia]+abs_dr2[ja]-2*dri_drj);
	if(Coup_dir|| STRAIN){
	  float dd=0; for(j=0; j<3; j++)dd+=dir[ia][j]*dir[ja][j];
	  if(Coup_dir)Coup_dir[ia][ja]+=w*dd;   // d_i*d_j   d=r/|r|
	  if(STRAIN){
	    float ss=w*strain[ia]*strain[ja];
	    Coup_str[ia][ja]+=ss;                   // strain_i*strain_j
	    Coup_str_dir[ia][ja]+=dd*ss;            // (str_i)d_i*(str_j)d_j
	  }
	}
      }
    }
  }


  // Empty 1
  for(i=0; i<Na; i++)Empty_matrix_f(r0ijk[i], Na);
  free(r0ijk);
  Empty_matrix_i(clist, Na);
  Empty_matrix_i(cnum, Na);
  Empty_matrix_f(dr, Na);
  Empty_matrix_f(dir, Na);

  // Normalize
  if(Coup_dir){
    for(ia=0; ia<Na; ia++){
      for(ja=0; ja<=ia; ja++)Coup_dir[ia][ja]/=norm_w;
    }
  }
  if(Coup_Bahar){
    for(ia=0; ia<Na; ia++){
      for(ja=0; ja<=ia; ja++)Coup_Bahar[ia][ja]=sqrt(Coup_Bahar[ia][ja]);
    }
  }

  // Fill lower diagonal
  if(Coup_dir)Fill_low_diag(Coup_dir, Na);
  if(Coup_coord)Fill_low_diag(Coup_coord, Na);
  if(Coup_dr_dir)Fill_low_diag(Coup_dr_dir, Na);
  if(Coup_str)Fill_low_diag(Coup_str, Na);
  if(Coup_str_dir)Fill_low_diag(Coup_str_dir, Na);
  if(Coup_Bahar)Fill_low_diag(Coup_Bahar, Na);

  /*
  if(Coup_dir){
    for(ia=0; ia<Na; ia++){
      Coup_dir[ia][ia]=1;
      for(ja=0; ja<ia; ja++){
	Coup_dir[ia][ja]=0.5*(1+Coup_dir[ia][ja]);
	// Make directionality positive: negative values -> 0
	if(Coup_dir[ia][ja] < 0)Coup_dir[ia][ja]=0;
	Coup_dir[ja][ia]=Coup_dir[ia][ja];
      }
    }
    } */
  //Get_zeta(Coup_dir, Na); // Transform into Z score!


  ////////////////////////////////////////////////////////////
  // DEFORMATION COUPLING
  // Receiver profile with force that maximizes deformation
  // suffered by residue i
  float **Coup_deformation=NULL;
  if(PRINT_DEF_COUPLING){
    printf("Computing Deformation coupling\n");
    Coup_deformation=Allocate_mat2_f(Na,Na);
    for(ia=0; ia<Na; ia++){
      int i3=3*iref[ia], i, j, k;
      for(ja=0; ja<=ia; ja++){
	int j3=3*iref[ja];
	double **F_matrix=Allocate_mat2_d(3, 3);
	double **Receiver_matrix=Allocate_mat2_d(3, 3);
	for(a=0; a<NM.N; a++){
	  if((NM.select[a]==0)||(sigma2[a]<=0))continue;
	  float w=sigma2[a]; //w*=w;
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
	//float Lambda=Singvalue_main(Receiver_matrix, 3, 3);
	float Lambda=Eigenvalue_main(Receiver_matrix, 3);
	// Rij and Rji have the same eigenvalues, since:
	// Fji=(Fij)t Rij=(Fij)t(Fij) Rji=(Fij)(Fij)t
	// Lambda/=abs_dr2[ia];
	// Optimal forced deformation at i divided by thermal deformation
	if(Lambda<0){
	  printf("WARNING, i=%d j=%d Lambda=%.3f\n", ia, ja, Lambda);
	  Lambda=0;
	}
	// Lambda=sqrt(Lambda)*norm_F; // Module, and not square displacement
	Lambda=sqrt(Lambda/mass); // Module, and not square displacement
	Coup_deformation[ia][ja]=Lambda;
	Coup_deformation[ja][ia]=Lambda;
	Empty_matrix_d(F_matrix, 3);
	Empty_matrix_d(Receiver_matrix, 3);
      }
    }
  }

  /************************************************************
                   Print couplings and profiles
  *************************************************************/

  // Print names
  char *pdbres[Na], amm[Na], ch[Na], ch_old='A';
  char *nameres[Na], chain_label[Na]; int Nchain=1;
  for(ia=0; ia<Na; ia++){
    struct residue *res=seq+(atoms+iatom[ia])->res;
    nameres[ia]=malloc(12*sizeof(char));
    int nres; sscanf(res->pdbres, "%d", &nres);
    sprintf(nameres[ia], "%c_%d_%c", res->amm, nres, res->chain);
    pdbres[ia]=malloc(6*sizeof(char));
    strcpy(pdbres[ia], res->pdbres);
    amm[ia]=res->amm;
    ch[ia]=res->chain;
    if(ia==0){
      ch_old=ch[0]; chain_label[0]=ch_old;
    }else if(ch[ia]!=ch_old){
      ch_old=ch[ia]; chain_label[Nchain]=ch[ia]; Nchain++; 
    }
  }
  Print_names(Na, pdbres, amm, ch, nameout1);

  ////////////////////////////////////////////////////////////
  // COORDINATION 
  if(PRINT_SIGMA_DIJ){
    Print_links(Coup_coord, Na, nameout,".DVAR.dat", 0);
    // YYY change output file name
  }
		
  // Transform coordination into similarity
  /*if(Coup_coord){
    for(ia=0; ia<Na; ia++){
      Coup_coord[ia][ia]=1;
      for(ja=0; ja<ia; ja++){
	Coup_coord[ia][ja]=1-sqrt(Coup_coord[ia][ja]); //d02ij[ia][ja]
	Coup_coord[ja][ia]=Coup_coord[ia][ja];
      }
    }
    }*/
  if(Coup_coord){ // similarity
    double rii=1, scale=SCALE_COORD; //1./0.1;
    for(ia=0; ia<Na; ia++){
      Coup_coord[ia][ia]=rii;
      for(ja=0; ja<ia; ja++){
	double rij2=dr2[ia]+dr2[ja]-2*Coup_dr_dir[ia][ja];
	//Coup_coord[ia][ja]=rii-scale*sqrt(rij2);
	Coup_coord[ia][ja]=exp(-scale*sqrt(rij2));
	// (1-sqrt(<|ri-rj|^2>))*scale
	//d02ij[ia][ja];
	Coup_coord[ja][ia]=Coup_coord[ia][ja];
      }
    }
  }
  // Get_zeta(Coup_coord, Na); // Transform into Z score!

  //////////////////////////////////////////////////////////////
  //  Coordination
  float *Prof_coord=NULL;
  if(PRINT_COORD_COUPLING){
    Print_links(Coup_coord, Na, nameout,
		".coordination_coupling.dat", 1); //_zeta
    Prof_coord=Make_profile(Coup_coord, Na, PROF_TYPE, "Coord");
  }

  ////////////////////////////////////////////////////////////
  //   Deformation
  float *Prof_deformation=NULL;
  if(PRINT_DEF_COUPLING){
    Print_links(Coup_deformation,Na,nameout,".deformation_coupling.dat",1);
    Prof_deformation=
      Make_profile(Coup_deformation, Na, PROF_TYPE, "Deformation");
  }
  ////////////////////////////////////////////////////////////
  //   Covariance
  //float *Prof_cov=NULL;
  if(PRINT_COV_COUPLING){
    Print_links(Coup_dr_dir, Na, nameout,
		".covariance_coupling.dat", 1); //_zeta
    //Prof_dir=Make_profile(Coup_dr_dir, Na, PROF_TYPE, "Covariance");
  }else if(Coup_dr_dir){
    Empty_matrix_f(Coup_dr_dir, Na);
    Coup_dr_dir=NULL;
  }


  ////////////////////////////////////////////////////////////
  //   Directionality
  float *Prof_dir=NULL;
  if(PRINT_DIR_COUPLING){
    Print_links(Coup_dir, Na, nameout,
		".directionality_coupling.dat", 1); //_zeta
    Prof_dir=Make_profile(Coup_dir, Na, PROF_TYPE, "Directionality");
  }

  ////////////////////////////////////////////////////////////
  //    Proximity
  float **Coup_proximity=NULL, *Prof_proximity=NULL;
  if(PRINT_COORD_COUPLING){
    float D=SCALE_D; // Scale for distances
    Coup_proximity=Negexp_distance(d02ij, Na, D);
    Print_links(Coup_proximity, Na, nameout,
		".proximity_coupling.dat", 1); //_zeta
    Prof_proximity=Make_profile(Coup_proximity, Na, PROF_TYPE, "proximity");
  }
  for(i=0; i<Na; i++)free(d02ij[i]);
  free(d02ij);

  // Other profiles
  float *Prof_dr=NULL;
  if(Coup_dr){
    Prof_dr=Make_profile(Coup_dr, Na, PROF_TYPE, "Dr");
  }

  float *Prof_dr_dir=NULL;
  if(PRINT_COV_COUPLING){
    Prof_dr_dir=Make_profile(Coup_dr_dir, Na, PROF_TYPE, "Dr_dir");
  }

  float *Prof_Bahar=NULL;
  if(Coup_Bahar){
    Prof_Bahar=Make_profile(Coup_Bahar, Na, PROF_TYPE, "Bahar");
  }

  float *Prof_str=NULL;
  if(Coup_str){
    Prof_str=Make_profile(Coup_str, Na, PROF_TYPE, "Coup_str");
  }
  float *Prof_str_dir=NULL;
  if(Coup_str_dir){
    Prof_str_dir=Make_profile(Coup_str_dir, Na, PROF_TYPE, "Coup_str_dir");
  }

  // Weighted_flux(Coup_dr, Na, "Coup_dr")  // CV_profile

  ////////////////////////////////////////////////////////////
  // BROADCASTER PROFILE
  if(0){
    // Broadcaster profile with force that maximizes deformation
    // produced by residue i
    double *Broadcast_profile=malloc(Na*sizeof(double));
    printf("Computing Broadcast profile\n");
    for(ia=0; ia<Na; ia++){
      double **Def=Allocate_mat2_d(3,3);
      int i3=3*iref[ia], i, j;
      for(a=0; a<NM.N; a++){
	float w=sigma2[a]; w*=w;
	float *Mode=NM.Cart[a]+i3;
	for(i=0; i<3; i++)for(j=0; j<=i; j++)Def[i][j]+=w*Mode[i]*Mode[j];
      }
      for(i=0; i<3; i++)for(j=i+1; j<3; j++)Def[i][j]=Def[j][i];
      Broadcast_profile[ia]=Eigenvalue_main(Def, 3);
      Empty_matrix_d(Def, 3);
    }
    Normalize_profile(Broadcast_profile, Na);
  }


  // Use directionality profile as a weight
  float *Prof_deformation_dir=NULL;
  if(0){
    Prof_deformation_dir=malloc(Na*sizeof(float));
    float weight_dir[Na], min_w=10;
    double norm_weight=0;
    for(ia=0; ia<Na; ia++){
      if(Prof_dir[ia]<min_w)min_w=Prof_dir[ia];
      norm_weight+=Prof_dir[ia];
    }
    norm_weight-=Na*min_w;
    for(ia=0; ia<Na; ia++){
      weight_dir[ia]=(Prof_dir[ia]-min_w)/norm_weight;
    }
    for(ia=0; ia<Na; ia++){
      double sum=0; float *Ca=Coup_deformation[ia];
      for(ja=0; ja<Na; ja++)sum+=Ca[ja]*weight_dir[ja];
      Prof_deformation_dir[ia]=sum/(1.-weight_dir[ia]);
    }
  }

  float *Prof_deformation_dist=NULL;
  float **Coup_deformation_dist=NULL;
  if(0){
    Prof_deformation_dist=malloc(Na*sizeof(float));
    Coup_deformation_dist=Allocate_mat2_f(Na, Na);
    for(ia=0; ia<Na; ia++){
      for(ja=ia+1; ja<Na; ja++){
	Coup_deformation_dist[ia][ja]=Coup_deformation[ia][ja]*(ja-ia);
	Coup_deformation_dist[ja][ia]=Coup_deformation_dist[ia][ja];
      }
    }
    Prof_deformation_dist=
      Make_profile(Coup_deformation_dist, Na, PROF_TYPE, "Deformation_dist");
  }
  
  //Deformation_propagation(Coup_deformation,Na,pdbres,amm,ch,nameout);
  
  /********************************************************/
  char nameout2[200]; FILE *file_out;
  // Correlation between allostery, conformation change, normal modes
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

  if(Prof_deformation && 0){
    // Correlation between normal mode and allostery
    float fluctuation[NM.N], dev[Na];
    int imax=0, NMODE=15; if(NMODE>NM.N)NMODE=NM.N;
    float cc[NMODE], cmax=-2;
    for(i=0; i<NMODE; i++){
      if(sigma2[i]<=0)continue;
      float *Cart=NM.Cart[i]; double sum=0;
      for(ia=0; ia<Na; ia++){
	double d=0; float *x=Cart+3*iref[ia];
	for(j=0; j<3; j++){d+=(*x)*(*x); x++;}
	dev[ia]=d; sum+=d;
      }
      cc[i]=Corr_coeff(dev, Prof_deformation, Na, NULL, NULL);
      if(cc[i]>cmax){cmax=cc[i]; imax=i;}
      fluctuation[i]=sqrt(sum/Na)*sqrt(sigma2[i]);
    }
    sprintf(nameout2, "%s.allostery_summary.dat", nameout);
    file_out=fopen(nameout2, "w");
    printf("Writing %s\n", nameout2);
    fprintf(file_out,
	    "# Maximum correlation def. prof. normal mode: %.3f mode=%d\n",
	    cmax, imax);
    fprintf(file_out, "# mode\tmean_fluct\tcorr(deformation)\n");
    for(i=0; i<NMODE; i++){
      fprintf(file_out, "%d\t%.4f\t%.3f\n", i, fluctuation[i], cc[i]);
    }
    fclose(file_out);
  }


  /*********************************************************
                  Binding sites
  **********************************************************/
  if(Prof_dir || Prof_coord || Prof_proximity || 
     Prof_deformation || Prof_dr_dir){

    int site[Na]; for(i=0; i<Na; i++)site[i]=0;
    Binding_sites(site,
		  Prof_dir, Coup_dir,
		  Prof_coord, Coup_coord,
		  Prof_proximity, Coup_proximity, 
		  Prof_deformation, Coup_deformation,
		  //Prof_deformation_dist, Coup_deformation_dist,
		  Prof_dr_dir, Coup_dr_dir,
		  atomres, pdbres, amm,
		  Na, chain, Nchain, seq, Nres, pdb, nameout, SITES);

    /***************************************************************
                               Interface
    ***************************************************************/
    
    double **norm_ch=NULL, **Chain_cov=NULL, **Chain_dir=NULL,
      **Chain_def=NULL, **Chain_coord=NULL;
    if(Nchain > 1){
      norm_ch=Allocate_mat2_d(Nchain, Nchain);
      //Chain_cov=Chain_coup(Coup_dr_dir, norm, Coup_proximity, Nchain, ch, Na);
      Chain_dir=Chain_coup(Coup_dir, norm_ch, Coup_proximity,Nchain,ch,Na);
      Chain_def=Chain_coup(Coup_deformation,norm_ch,Coup_proximity,Nchain,ch,Na);
      Chain_coord=Chain_coup(Coup_coord, norm_ch,Coup_proximity,Nchain,ch,Na);
      Coupling_chains(Chain_dir, Chain_coord, Chain_def, Chain_cov,
		      norm_ch, ch, Na, chain_label, Nchain, nameout);
    }

    if(Int_list_inter){
      int *nc_inter, **clist_inter, **cnum_inter;
      Get_contact_list(&nc_inter, &clist_inter, &cnum_inter,
		       Na, atoms, Int_list_inter, N_int_inter);
      Print_interface_coupling(Int_list_inter, N_int_inter,
			       Int_list, N_int, atoms, Nchain, nameout,
			       Coup_dir, Coup_coord, Coup_deformation,
			       Chain_dir, Chain_coord, Chain_def,
			       norm_ch, Na);
    }

    if(Chain_cov)Empty_matrix_d(Chain_cov, Nchain);
    if(Chain_dir)Empty_matrix_d(Chain_dir, Nchain);
    if(Chain_def)Empty_matrix_d(Chain_def, Nchain);
    if(Chain_coord)Empty_matrix_d(Chain_coord, Nchain);
    if(norm_ch)Empty_matrix_d(norm_ch, Nchain);

    // Print matrix of profiles
    strcat(nameout, ".profiles.dat");
    file_out=fopen(nameout, "w");
    printf("Writing %s\n", nameout);
    
    char profname[15];
    if(PROF_TYPE=='A'){sprintf(profname, "Ave");}
    else if(PROF_TYPE=='P'){sprintf(profname, "PE");}
    else{sprintf(profname, "EC");}
    fprintf(file_out, "# Profiles of type %s N=%d\n", profname, Na);
    fprintf(file_out, "### res AA binding");
    char hd[80]; sprintf(hd, "### 0 0 0"); int np=4;
    if(Prof_proximity){
      strcat(hd, " 1"); fprintf(file_out, "\t%d=proximity", np); np++;
    }
    if(Prof_coord){
      strcat(hd, " 1");
      fprintf(file_out, "\t%d=coord=exp(-sigma^2(d_ij)/D)", np); np++;
    }
    if(Prof_dir){
      strcat(hd, " 1");
      fprintf(file_out, "\t%d=dir=<dri*drj/|dri||drj|>", np); np++;
    }
    if(Prof_deformation){
      strcat(hd, " 1");
      fprintf(file_out, "\t%d=deformation=dri/dfj", np); np++;
    }
    if(Prof_deformation_dist){
      strcat(hd, " 1");
      fprintf(file_out, "\t%d=deformation_dist=dri/dfj(|i-j|)", np); np++;
    }
    if(Prof_dr_dir){
      strcat(hd, " 1"); fprintf(file_out, "\t%d=cov=<dri*drj>", np); np++;
    }
    if(Prof_dr){
      strcat(hd, " 1"); fprintf(file_out, "\t%d=<|dri||drj|>", np); np++;
    }
    if(Prof_Bahar){
      strcat(hd, " 1"); fprintf(file_out, "\t%d=<|ri-rj|^2>Bahar", np); np++;
    }
    if(Prof_str){
      strcat(hd, " 1"); fprintf(file_out, "\t%d=<s_i*s_j>", np); np++;
    }
    if(Prof_str_dir){
      strcat(hd, " 1"); fprintf(file_out, "\t%d=<s_i d_i*s_j d_j>", np); np++;
    }
    if(strdiff){
      strcat(hd, " 1"); fprintf(file_out, "\t%d=Confchange", np); np++;
    }
    fprintf(file_out, "\n%s\n", hd);
    
    for(i=0; i<Na; i++){
      fprintf(file_out, "%s\t%c\t%d",pdbres[i], amm[i], site[i]);
      if(Prof_proximity)fprintf(file_out, "\t%.3g", Prof_proximity[i]);
      if(Prof_coord)fprintf(file_out, "\t%.3g", Prof_coord[i]);
      if(Prof_dir)fprintf(file_out, "\t%.3g", Prof_dir[i]);
      if(Prof_deformation)fprintf(file_out, "\t%.3g", Prof_deformation[i]);
      if(Prof_deformation_dist)
	fprintf(file_out, "\t%.3g", Prof_deformation_dist[i]);
      
      if(Prof_dr_dir)fprintf(file_out, "\t%.3g", Prof_dr_dir[i]);
      if(Prof_dr)fprintf(file_out, "\t%.3g", Prof_dr[i]);
      if(Prof_Bahar)fprintf(file_out, "\t%.3g", Prof_Bahar[i]);
      if(Prof_str)fprintf(file_out, "\t%.3g", Prof_str[i]);
      if(Prof_str_dir)fprintf(file_out, "\t%.3g", Prof_str_dir[i]);
      if(strdiff)fprintf(file_out, "\t%.3g", strdiff[i]);
      fprintf(file_out, "\n");

    }
    fclose(file_out);
  }

  // Empty
  if(strdiff)free(strdiff);
 
  if(Coup_proximity){
    Empty_matrix_f(Coup_proximity, Na); free(Prof_proximity);
  }
  if(Coup_dr) {Empty_matrix_f(Coup_dr, Na); free(Prof_dr);}
  if(Coup_dir){Empty_matrix_f(Coup_dir, Na); free(Prof_dir);}
  if(Coup_dr_dir){Empty_matrix_f(Coup_dr_dir, Na); free(Prof_dr_dir);}
  if(Coup_Bahar) {Empty_matrix_f(Coup_Bahar, Na);  free(Prof_Bahar);}
  if(Coup_coord) {Empty_matrix_f(Coup_coord, Na); free(Prof_coord);}
  if(Coup_str)   {Empty_matrix_f(Coup_str, Na); free(Prof_str);}
  if(Coup_str_dir){Empty_matrix_f(Coup_str_dir, Na); free(Prof_str_dir);}
  if(Coup_deformation){
    Empty_matrix_f(Coup_deformation, Na); free(Prof_deformation);
  }
  if(Coup_deformation_dist){
    Empty_matrix_f(Coup_deformation_dist, Na); free(Prof_deformation_dist);
  }
  for(ia=0; ia<Na; ia++){
    free(pdbres[ia]); free(nameres[ia]);
  }
}

float PE_profile(float *PE, double **Mat, int N)
{
  float eigenval[N];
  float **eigenvec=malloc(N*sizeof(float *));
  for(int i=0; i<N; i++)eigenvec[i]=malloc(N*sizeof(float));
  d_Diagonalize(N, Mat, eigenval, eigenvec, 1);
  float *v=eigenvec[0];
  for(int i=0; i<N; i++)PE[i]=v[i];
  Empty_matrix_f(eigenvec, N);
  return(eigenval[0]);
}

float Eigenvalue_main(double **Corr, int Na)
{
  int i;
  float *eigenval=malloc(Na*sizeof(double));
  float **eigenvec=malloc(Na*sizeof(float *));
  for(i=0; i<Na; i++)eigenvec[i]=malloc(Na*sizeof(float));
  d_Diagonalize(Na, Corr, eigenval, eigenvec, 1);
  float Ev=eigenval[0];
  Empty_matrix_f(eigenvec, Na); free(eigenval);
  return(Ev);
}

float Singvalue_main(double **A, int Na, int Nb)
{
  int i, j, k, N1, N2;
  if(Na<Nb){N1=Na; N2=Nb;}else{N1=Nb; N2=Na;}
  float *eigenval=malloc(N1*sizeof(float));
  float **eigenvec=Allocate_mat2_f(N1, N1);
  double **Mat=Allocate_mat2_d(N1, N1);
  for(i=0; i<N1; i++){
    for(j=0; j<N1; j++){
      for(k=0; k<N2; k++)Mat[i][j]+=A[i][k]*A[j][k];
    }
  }

  d_Diagonalize(N1, Mat, eigenval, eigenvec, 1);
  float Ev=sqrt(eigenval[0]);
  Empty_matrix_d(Mat, Na);
  Empty_matrix_f(eigenvec, Na);
  //for(i=0; i<N1; i++)printf("%.3g ", eigenval[i]); printf("\n");
  free(eigenval);
  return(Ev);
}


void Normalize_profile(double *v, int N){
  double sum=0; int i;
  for(i=0; i<N; i++)sum+=v[i]*v[i];
  sum=sqrt(sum/N);
  for(i=0; i<N; i++)v[i]/=sum;
}



int Count_chains(int *ichain, int *min_chain, int *max_chain, char *chain_lab,
		 char *chain_res, int *site_res, int nr)
{
  int nc=0, j, k;
  for(int i=0; i<nr; i++){
    for(j=0; j<nc; j++){
      if(chain_lab[j]==chain_res[i])break;
    }
    ichain[i]=j;
    //num_chain[j]++;
    if(site_res){k=site_res[i];}
    else{k=i;}
    if(j==nc){
      chain_lab[nc]=chain_res[i];
      min_chain[j]=k;
      max_chain[j]=k;
      nc++;
    }else{
      if(k<min_chain[j])min_chain[j]=k;
      if(k>max_chain[j])max_chain[j]=k;
    }
  }
  return(nc);
}

void Identify_chains(int *all_chain, char *chain_lab, int nchain,
		     char *lab_chain, int Nchain)
{
  for(int nc=0; nc<nchain; nc++){
    all_chain[nc]=-1;
    for(int k=0; k<Nchain; k++){
      if(chain_lab[nc]==lab_chain[k]){all_chain[nc]=k; break;}
    }
    if(all_chain[nc]<0){
      printf("ERROR, chain %c not found among the %d chains\n",
	     chain_lab[nc], Nchain); exit(8);
    }
  }
}

int Print_links(float **Coupling, int N, char *nameout1, char *ext, int sign)
{
  char nameout[200]; int i, j;
  sprintf(nameout, "%s%s", nameout1, ext);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);

  // Mean and standard deviation
  double sum1=0, sum2=0; int norm=0;
  float max=0, min=1000;
  for(i=0; i<N; i++){
    float *c=Coupling[i]+i+2;
    for(j=i+2; j<N; j++){
      sum1+=*c; sum2+=(*c)*(*c); norm++;
      if(*c>max){max=*c;}
      if(*c<min){min=*c;} c++;
    }
  }
  sum1/=norm; sum2=sqrt((sum2-norm*sum1*sum1)/(norm-1));
  if(sum1==0 || sum2==0 || isnan(sum1)|| isnan(sum2)){
    printf("ERROR in Print_links, coupling= %s\n", ext);
    for(i=0; i<N; i++){
      float *c=Coupling[i];
      for(j=i+1; j<N; j++){printf("%.3g ",c[j]);} printf("\n");
    }
    printf("ERROR in Print_links, coupling= %s\n", ext);
    printf("Mean= %.3g s.d.= %.3g\n", sum1, sum2); exit(8);
  }

  double thrmax=sum1+SIGMA*sum2, thrmin=sum1-SIGMA*sum2;
  if(thrmax>=max){thrmax=(max+sum1)/2;}
  if(thrmin<=min){thrmin=(min+sum1)/2;}
  fprintf(file_out, "# mean= %.3g s.d.= %.3g thr=%.2g %.2g SIGMA= %.3g\n",
	  sum1, sum2, thrmin, thrmax, SIGMA);
  fprintf(file_out, "#Residue names written in %s.names.dat\n", nameout1);
  fprintf(file_out, "#resi resj coupling\n");
  for(i=0; i<N; i++){
    float *c=Coupling[i];
    for(j=i+1; j<N; j++){
      if((sign==0)||
	 ((sign>0)&&(c[j]>=thrmax))||
	 ((sign<0)&&(c[j]<=thrmin))){
	fprintf(file_out, "%d %d %.3g\n", i, j, c[j]);
      }
    }
  }
  fclose(file_out);
  return(0);
}

void Print_names(int N, char **pdbres, char *amm, char *chain,
		 char *nameout1)
{
  int i;  char nameout[200];
  sprintf(nameout, "%s.names.dat", nameout1);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  for(i=0; i<N; i++){
    fprintf(file_out, "%d %c_%s_%c\n", i, amm[i], pdbres[i], chain[i]);
  }
  fclose(file_out);
}

int Get_zeta(float **Coupling, int N)
{
  // Mean and standard deviation depending on l=|i-j|
  // Compute z=(c_{i,i+l}-<c>_l)/sigma_l
  // If |z| > 1, amplify c -> c*z
  int i, j, l, ll=-1, num_l=0;;
  float *c1=malloc(N*sizeof(float));
  float *c2=malloc(N*sizeof(float));
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
      if(ll<0){ll=l;} sum1_l+=sum1; sum2_l+=sum2; num_l+=num;
    }
  }
  sum1_l/=num_l;
  sum2_l=sqrt((sum2_l-num_l*sum1_l*sum1_l)/(num_l-1));
  for(l=ll; l<N; l++){
    c1[l]=sum1_l;
    c2[l]=sum2_l;
  }

  for(i=0; i<N; i++){
    float *c=Coupling[i];
    for(j=i+1; j<N; j++){
      float z=(c[j]-c1[j-i])/c2[j-i];
      if(z<0)z=-z;
      if(z>1){c[j]*=z; Coupling[j][i]=c[j];}
    }
  }
  free(c1); free(c2);
  return(0);
}

int Test_site(struct site_results *res, int *atom_site,
	      float *profile, float **coupling, int Na,
	      int *site_res, int **site_res_ran, int numran,
	      int nr, int *atomres)
{

  int n_not=0;

  // Observed site
  double p_ave, p_dev, c_ave, c_dev;
  int n_notfound=0;
  p_ave=Ave_site(&p_dev, &c_ave, &c_dev, &n_notfound,
		 site_res, atomres, nr, atom_site,
		 profile, coupling, Na);
  if(n_notfound){
    printf("WARNING, %d reference atoms of binding site not found\n",
	   n_notfound);
  } 

  // Sum over num random sites 
  double ran_prof1=0, ran_prof2=0;
  double ran_pair1=0, ran_pair2=0;
  double ran_sig1=0,  ran_sig2=0;
  int num_p_gt=0, num_c_gt=0;
  for(int kran=0; kran<numran; kran++){
    double ave_prof, sigprof, ave_pair, sigpair;
    ave_prof=Ave_site(&sigprof, &ave_pair, &sigpair, &n_notfound,
		      site_res_ran[kran], atomres, nr,
		      NULL, profile, coupling, Na);
    if(ave_prof > p_ave)num_p_gt++;
    if(ave_pair > c_ave)num_c_gt++;
    ran_prof1+=ave_prof;
    ran_prof2+=ave_prof*ave_prof;
    ran_sig1+=sigprof;
    ran_sig2+=sigprof*sigprof;
    ran_pair1+=ave_pair;
    ran_pair2+=ave_pair*ave_pair;
    n_not+=n_notfound;
  } 
  Ave_St_dev(&ran_prof1, &ran_prof2, numran);
  Ave_St_dev(&ran_pair1, &ran_pair2, numran);
  Ave_St_dev(&ran_sig1, &ran_sig2, numran);
  if(n_not)
    printf("WARNING, %d atoms not found in %d random sites\n", n_not, numran);
  
  float Z_p=(p_ave-ran_prof1)/(ran_prof2);
  float Z_c=(c_ave-ran_pair1);
  if(ran_pair2)Z_c/=ran_pair2;
  //float Z_sig= (p_dev - ran_sig1)/(ran_sig2);
  res->prof=p_ave; res->Z_p=Z_p; res->P_p= (float)num_p_gt/numran;
  res->coup=c_ave; res->Z_c=Z_c; res->P_c= (float)num_c_gt/numran;
  res->p_ran1=ran_prof1; res->p_ran2=ran_prof2;
  res->c_ran1=ran_pair1; res->c_ran2=ran_pair2;

  return(0);
}
 
float Ave_site(double *sigprof, double *ave_pair,
	       double *sigpair, int *n_notfound,
	       int *site, int *atomres, int nr, int *atom_site,
	       float *profile, float **coupling, int Na)
{
  double sum_prof=0; 
  *sigprof=0; *ave_pair=0; *sigpair=0; *n_notfound=0;
  int npair=0;
  for(int i1=0; i1<nr; i1++){
    int j1=atomres[site[i1]];
    if((j1<0)||(j1>=Na)){
      if(atom_site)atom_site[j1]=-site[i1];
      (*n_notfound)++;
      continue;
    }
    if(atom_site)atom_site[j1]=1;
    sum_prof+=profile[j1];
    *sigprof+=profile[j1]*profile[j1];
    for(int i2=0; i2<i1; i2++){
      int j2=atomres[site[i2]];
      if((j2<0)||(j2>=Na)||(j2==j1))continue;
      float x=coupling[j1][j2];
      *ave_pair+=x; *sigpair+=x*x;
      npair++;
    }
  } // end of the site
  Ave_St_dev(&sum_prof, sigprof, nr-*n_notfound);
  Ave_St_dev(ave_pair, sigpair, npair);
  return(sum_prof);
}

int Test_pairs(struct site_results *res, int *atom_site,
	       float *profile, float **coupling, int Na,
	       int *coupled_site, int **coupled_site_ran, int numran,
	       int ni, int nij, int *atomres)
{

  int n_not=0;

  // Observed site
  double c_ave, c_dev;
  int n_notfound=0;
  c_ave=Pair_site(&c_dev, &n_notfound,
		  coupled_site, ni, nij, atomres, atom_site, coupling, Na);
  if(n_notfound){
    printf("WARNING in Test_pairs, %d reference atoms not found\n",
	   n_notfound);
  }
      
  // Sum over num random sites 
  double ran_pair1=0, ran_pair2=0;
  int num_c_gt=0;
  for(int kran=0; kran<numran; kran++){
    double ave_pair, sigpair;
    ave_pair=Pair_site(&sigpair, &n_notfound,
		       coupled_site_ran[kran], ni, nij,
		       atomres, NULL, coupling, Na);
    if(ave_pair > c_ave)num_c_gt++;
    ran_pair1+=ave_pair;
    ran_pair2+=ave_pair*ave_pair;
    n_not+=n_notfound;
  }
  if(n_not)
    printf("WARNING, %d atoms not found in %d random sites\n",
	   n_not, numran);
  
  Ave_St_dev(&ran_pair1, &ran_pair2, numran);
  float Z_c=(c_ave-ran_pair1);
  if(ran_pair2)Z_c/=ran_pair2;
  res->coup=c_ave; res->Z_c=Z_c; res->P_c= (float)num_c_gt/numran;
  return(0);
}

float Pair_site(double *sigpair, int *n_notfound,
		int *coupled_site, int n1, int n12,
		int *atomres, int *atom_site,
		float **coupling, int Na)
{
  double ave_pair=0; *sigpair=0; *n_notfound=0;
  int npair=0;
  for(int i1=0; i1<n1; i1++){
    int j1=atomres[coupled_site[i1]];
    if((j1<0)||(j1>=Na)){
      if(atom_site)atom_site[j1]=-coupled_site[i1];
      (*n_notfound)++;
      continue;
    }
    if(atom_site)atom_site[j1]=1;
    for(int i2=n1; i2<n12; i2++){
      int j2=atomres[coupled_site[i2]];
      if((j2<0)||(j2>=Na)){
	if(atom_site)atom_site[j1]=-coupled_site[i2];
	(*n_notfound)++;
	continue;
      }
      if(j2==j1)continue;
      float x=coupling[j1][j2];
      ave_pair+=x; *sigpair+=x*x;
      npair++;
    }
  } // end of the site
  Ave_St_dev(&ave_pair, sigpair, npair);
  return(ave_pair);
}

 
void Random_sites(int *site_ran, int nres, int *site_res,
		  int *ichain, int *min_chain, int *max_chain,
		  int *all_chain, int nchain,
		  int *ini_chain, int *end_chain)
{
  int THR=6, RMAX=40, i, ic;
  //float RANFACTOR=pow(12.0,1/3.0);
  unsigned long iran=randomgenerator();
  if(ini_ran==0){
    InitRandom((RANDOMTYPE)iran); ini_ran=1;
  }
  for(i=0; i<nres; i++)site_ran[i]=-1;
  for(ic=0; ic<nchain; ic++){
    // Shift all chains independently
    int c=all_chain[ic];
    int shiftmin=ini_chain[c]-min_chain[ic];
    int shiftmax=end_chain[c]-max_chain[ic];
    int shift; i=0;
    while(i<RMAX){
      shift=RandomFloating()*(shiftmax-shiftmin)+shiftmin;
      if((shift>THR)||(shift<-THR))break;
      i++;
    }
    if(i==RMAX){
      printf("WARNING, no allowed random site found in %d attempts\n",
	     RMAX); shift=0;
    }
    for(i=0; i<nres; i++){
      if(ichain[i]==ic){
	site_ran[i]=site_res[i]+shift;
	// Check
	if((site_ran[i]<ini_chain[c])||(site_ran[i]>end_chain[c])){
	    printf("ERROR, incorrect random site %d = %d\n", i,site_ran[i]);
	    printf("Chain %d in site, %d in prot\n", ic, c);
	    printf("Site: ini %d end %d Protein: ini %d end %d\n",
		   min_chain[ic], max_chain[ic], ini_chain[c], end_chain[c]);
	    exit(8);
	}
      }
    }
  }
  for(i=0; i<nres; i++){
    if(site_ran[i]<0){
      printf("ERROR, residue %d not set\n", i); exit(8);
    }
  }
  return;
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

float *Make_profile(float **Mat, int N, char PROF_TYPE, char *what)
{
  printf("Computing profiles of %s type %c\n", what, PROF_TYPE);
  float *prof=malloc(N*sizeof(float)); int i;
  float weight[N], sum_weight=0;
  if(PROF_TYPE=='A'){
    sum_weight=N; for(i=0; i<N; i++)weight[i]=1;
  }else{
    double **Mat2=Get_matrix_positive(Mat, N);
    if(PROF_TYPE=='P'){
      PE_profile(weight, Mat2, N);
    }else if(PROF_TYPE=='C'){
      float *c=EC_profile(NULL, Mat2, N, what);
      for(i=0; i<N; i++)weight[i]=c[i];
      free(c);
    }
    for(i=0; i<N; i++)sum_weight+=weight[i];
    Empty_matrix_d(Mat2, N);
  }
  for(i=0; i<N; i++){
    double sum=0; float *Ca=Mat[i];
    for(int j=0; j<N; j++)if(j!=i)sum+=Ca[j]*weight[j];
    prof[i]=sum/(sum_weight-weight[i]);
  }
  return(prof);
}

double **Get_matrix_positive(float **Mat, int N)
{
  double **Mat2=malloc(N*sizeof(double *)); int i, j;
  double min=0;
  for(i=0; i<N; i++){
    Mat2[i]=malloc(N*sizeof(double));
    for(j=0; j<N; j++){
      Mat2[i][j]=Mat[i][j];
      if(Mat2[i][j]<min)min=Mat2[i][j];
    }
  }
  // Transform into positive
  if(min<0){
    for(i=0; i<N; i++){
      for(j=0; j<N; j++)Mat2[i][j]-=min;
    }
  }
  return(Mat2);
}

int Extract_atoms(int *iref, int *iatom, int *atomres, // Output
		  struct residue *seq, float *mass,
		  atom *atoms, struct Reference Ref, char *ANAME)
{
  /* 
     iref[i]= representative atom (kin) of representative atom i (allostery)
     iatom=   atom number of representative atom i (allostery)
     atomres= representative atom (allostery) of residue
   */
  int ia, Na=0, Nres=atoms[Ref.atom_num[Ref.N_ref-1]].res+1;
  for(ia=0; ia<Nres; ia++)atomres[ia]=-1;
  for(ia=0; ia<Ref.N_ref; ia++){
    atom *atom_i=atoms+Ref.atom_num[ia];
    if((strncmp(atom_i->name, ANAME, 2)==0)||
       ((seq[atom_i->res].amm=='G')&&(strncmp(atom_i->name, "CA", 2)==0))){
      if(*mass==0)*mass=Ref.mass_atom[ia];
      iref[Na]=ia;
      iatom[Na]=Ref.atom_num[ia];
      //ires[Na]=atom_i->res;
      atomres[atom_i->res]=Na; Na++;
    }
  }
  /*for(ia=0; ia<Na; ia++){
    int r=(atoms+Ref.atom_num[iref[ia]])->res;
    printf("%d %c %s %s\n", ia, seq[r].amm, seq[r].pdbres,
    (atoms+Ref.atom_num[iref[ia]])->name);
   }
   exit(8);*/
  return(Na);
}


void Pairwise_distances(float **d02ij, float ***r0ijk, 
			atom *atoms, int *iatom, int Na)
{
  int ia, ja, k;
  for(ia=0; ia<Na; ia++){
    double *ri=atoms[iatom[ia]].r;
    float *dij=malloc(Na*sizeof(float));
    d02ij[ia]=dij;
    float **rij=Allocate_mat2_f(Na,3);
    r0ijk[ia]=rij;
    for(ja=0; ja<ia; ja++){
      double *rj=atoms[iatom[ja]].r, sum=0;
      float *d=rij[ja];
      for(k=0; k<3; k++){
	*d=ri[k]-(*rj); sum+=(*d)*(*d); d++; rj++;
      }
      dij[ja]=sum;
      sum=sqrt(sum);
      d=rij[ja]; for(k=0; k<3; k++){(*d)/=sum; d++;}
    }
  }
}


float **Negexp_distance(float **d02ij, int Na, float D)
{
  float **eij=Allocate_mat2_f(Na, Na);
  for(int i=0; i<Na; i++){
    eij[i][i]=1;
    for(int j=0; j<i; j++){
      eij[i][j]=exp(-sqrt(d02ij[i][j])/D);
      eij[j][i]=eij[i][j];
    }
  }
  return(eij);
}

void Fill_low_diag(float **Corr, int Na){
  int ia, ja;
  for(ia=0; ia<Na; ia++){
    for(ja=0; ja<ia; ja++)Corr[ja][ia]=Corr[ia][ja];
  }
}

void Deformation_propagation(float **Deformation_coupling, int Na,
			    char **pdbres, char *amm, char *chain,
			    char *nameout)
{
  char namefile[200];
  sprintf(namefile, "%s.perturbation_propagation.dat", nameout);
  FILE *file_out=fopen(namefile, "w");
  printf("Writing %s\n", namefile);
  fprintf(file_out, "# Predicted deformation produced on sites j by a ");
  fprintf(file_out, "perturbation applied at site i:\n");
  fprintf(file_out, "# Def(i->j)~Allo_i exp(-|i-j|/prop_length_i)\n");
  fprintf(file_out, "#AA res chain Allo prop_length C.c.\n");
  float *d=malloc(Na*sizeof(float));
  float *y=malloc(Na*sizeof(float));
  float slope, offset, r; int i, j;
  for(i=0; i<Na; i++){
    float *Ci=Deformation_coupling[i];
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

void Print_sites(int **site_res, int *nres_site,
		 char **description, int nsites, 
		 struct residue *seq, char *pdb, char *chain)
{
  // Printing
  char pdbid[80], nameout[200], amm[5]; int k, j;
  GetPdbId(pdb,pdbid);
  sprintf(nameout, "%s%s_sites.in", pdbid, chain);
  printf("Writing binding sites in %s\n", nameout);
  FILE *file_out=fopen(nameout, "w");
  for(k=0; k<nsites; k++){
    if(description[k])
      fprintf(file_out, "SITE_DESCRIPTION: %s", description[k]);
    for(j=0; j<nres_site[k]; j++){
      struct residue *res=seq+site_res[k][j];
      Name3(amm, res->i_aa);
      fprintf(file_out, "%d\t%s\t%c\t%s\n",
	      k+1, amm, res->chain, res->pdbres);
    }
  }
  fclose(file_out);
}

void Ave_St_dev(double *sum1, double *sum2, int n){
  if(n<=0)return;
  if(n==1){*sum2=0; return;}
  (*sum1)/=n;
  (*sum2)=*sum2-n*(*sum1)*(*sum1);
  (*sum2)=sqrt((*sum2)/(n-1));
}

void Binding_sites(int *site,
		   float *Prof_dir, float **Coup_dir,
		   float *Prof_coord, float **Coup_coord,
		   float *Prof_proximity, float **Coup_proximity, 
		   float *Prof_deformation, float **Coup_deformation,
		   //float *Prof_deformation_dist,float **Coup_deformation_dist,
		   float *Prof_dr_dir, float **Coup_dr_dir,
		   int *atomres, char **pdbres, char *amm, 
		   int Na, char *chain, int Nchain,
		   struct residue *seq, int Nres,
		   char *pdb, char *nameout1, char *SITES)
{
  int numran=1000; int i, j, k;
  
  int NMAX=1000;
  int *site_res[NMAX], nres_site[NMAX];
  char *chain_res[NMAX], *description[NMAX], file_site[200];
  int PRINT_SITE=1;  // Print active site read in the PDB?
  printf("Reading binding sites\n");
  if(SITES[0]!='\0'){printf("Binding site file: %s\n", SITES);}
  int N_sites=
    Read_sites(site_res, nres_site, chain_res, description, file_site,
	       NMAX, SITES, pdb, chain, Nchain, seq, Nres,PRINT_SITE);
  printf("%d binding sites read\n", N_sites);
  if(N_sites==0)return;

  /*
  int NMAX=1000, SMAX=80, N_sites=0;
  int nres_site[NMAX]; char *chain_res[NMAX];
  for(i=0; i<NMAX; i++)chain_res[i]=malloc(SMAX*sizeof(char));
  char *description[NMAX]; for(i=0; i<NMAX; i++)description[i]=NULL;
  int **site_res=Allocate_mat2_i(NMAX,SMAX);
  char file_site[200];
  if(SITES[0]!='\0'){
    N_sites=Read_sites_file(site_res, chain_res, nres_site, description,
			    SITES, chain, seq, Nres, NMAX, SMAX);
    if(N_sites)strcpy(file_site, SITES);
  }
  if(N_sites==0){
    N_sites=Read_sites_PDB(site_res, chain_res, nres_site, description,
			   pdb, chain, seq, Nres, NMAX, SMAX);
    if(N_sites)strcpy(file_site, pdb);
    // Print site
    if(PRINT_SITE && N_sites)
      Print_sites(site_res, nres_site, description, N_sites,
		  seq, pdb, chain);
  }
  if(N_sites==0)return;
  for(i=0; i<N_sites; i++){
    if(description[i]==NULL){
      description[i]=malloc(10*sizeof(char));
      sprintf(description[i], "site%d", i+1);
    }else{
      char *ptr=description[i];
      while(*ptr!='\0'){
	if(*ptr=='\n'){*ptr='\0'; break;} ptr++;
      } 
    }
  }
  int S_sites=0;
  for(i=0; i<N_sites; i++){
    if(nres_site[i]>0)S_sites++;
  }
  if(S_sites==0)return;
  */

  // Testing active sites
  int numcoup=5, icoup=0;
  float *prof[numcoup], **coup[numcoup];
  char name_p[numcoup][40];
  if(Coup_dir){
    prof[icoup]=Prof_dir;
    coup[icoup]=Coup_dir;
    strcpy(name_p[icoup], "directionality");
    icoup++;
  }
  if(Coup_coord){
    prof[icoup]=Prof_coord;
    coup[icoup]=Coup_coord;
    strcpy(name_p[icoup], "coordination");
    icoup++;
  }
  if(Coup_proximity){
    prof[icoup]=Prof_proximity;
    coup[icoup]=Coup_proximity;
    strcpy(name_p[icoup], "proximity");
    icoup++;
  }
  if(Coup_deformation){
    prof[icoup]=Prof_deformation;
    coup[icoup]=Coup_deformation;
    strcpy(name_p[icoup], "deformation");
    icoup++;
  }
  if(Coup_dr_dir){
    prof[icoup]=Prof_dr_dir;
    coup[icoup]=Coup_dr_dir;
    strcpy(name_p[icoup], "covariance");
    icoup++;
  }

  /*if(Coup_deformation_dist){
    prof[icoup]=Prof_deformation_dist;
    coup[icoup]=Coup_deformation_dist;
    strcpy(name_p[icoup], "deformation_dist");
    icoup++;
    }*/

  int *site_res_ran[numran];
  for(j=0; j<numran; j++)
    site_res_ran[j]=malloc(Nres*sizeof(int *));

  char nameout[200];
  sprintf(nameout, "%s.binding_sites.dat", nameout1);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  fprintf(file_out, "# read binding sites in file %s\n", file_site);
  fprintf(file_out, "# profile type %c\n", PROF_TYPE);
  fprintf(file_out, "# N= %d, %d random sites\n", Na, numran);
  fprintf(file_out, "### 1=nres\t2=nchain"); j=3;
  for(k=0; k<icoup; k++){
    fprintf(file_out, "\t%d=<%s_prof>\t%d=Z(%s_prof)\t%d=P>(%s_prof)",
	    j, name_p[k], j+1, name_p[k], j+2, name_p[k]); j+=3;
    if(PRINT_RAN){
      fprintf(file_out, "\t%d=E(%s_p_ran>)\t%d=sigma(%s_p_ran>)",
	      j, name_p[k], j+1, name_p[k]); j+=2;
    }
    fprintf(file_out, "\t%d=<%s_coup>\t%d=Z(%s_coup)\t%d=P>(%s_coup)",
	    j, name_p[k], j+1, name_p[k], j+2, name_p[k]); j+=3;
    if(PRINT_RAN){
      fprintf(file_out, "\t%d=E(%s_p_ran>)\t%d=sigma(%s_p_ran>)",
	      j, name_p[k], j+1, name_p[k]); j+=2;
    }
  }
  fprintf(file_out, "\n");

  // Set chains
  char res_chain[Nres]; for(i=0; i<Nres; i++)res_chain[i]=seq[i].chain;
  char lab_chain[Nres]; int ini_chain[Nres], end_chain[Nres], ichain[Nres];
  //int Nchain=
  Count_chains(ichain, ini_chain, end_chain, lab_chain,res_chain, NULL, Nres);

  struct site_results site_results[numcoup], *res;
  int atom_site[Na]; for(i=0; i<Na; i++)atom_site[i]=0;
  int min_chain[Na], max_chain[Na], all_chain[Na];
  char chain_lab[Na];
  for(i=0; i<N_sites; i++){
    int nr=nres_site[i]; //if(nr<=1)continue;
    int *isite=site_res[i];
    int nchain=Count_chains(ichain, min_chain, max_chain, chain_lab,
			    chain_res[i], isite, nr);
    Identify_chains(all_chain, chain_lab, nchain, lab_chain, Nchain);
    fprintf(file_out, "## %s\n", description[i]);
    fprintf(file_out, "%d\t%d", nr, nchain);
    for(j=0; j<numran; j++){
      Random_sites(site_res_ran[j], nr, isite,
		   ichain, min_chain, max_chain, all_chain, nchain,
		   ini_chain, end_chain);
    }
    res=site_results;
    for(k=0; k<icoup; k++){
      Test_site(res, atom_site, prof[k], coup[k], Na,
		isite, site_res_ran, numran, nr, atomres);
      fprintf(file_out, "\t%.3g\t%.1f\t%.3f",
	      res->prof, res->Z_p, res->P_p);
      if(PRINT_RAN)
	fprintf(file_out, "\t%.3g\t%.1f", res->p_ran1, res->p_ran2);
      fprintf(file_out, "\t%.3g\t%.1f\t%.3f",
	      res->coup, res->Z_c, res->P_c); //->c_ran1
      if(PRINT_RAN)
	fprintf(file_out, "\t%.3g\t%.1f", res->c_ran1, res->c_ran2);
      res++;
    }
    fprintf(file_out, "\n");
  }
  fclose(file_out);

  // Profile for binding residues
  /*float ave_p[icoup];
  for(k=0; k<icoup; k++){
    double sum=0; float *p=prof[k];
    for(i=0; i<Na; i++)sum+=p[i];
    ave_p[k]=sum/Na;
    }*/

  sprintf(nameout, "%s.binding_residues.dat", nameout1);
  file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  fprintf(file_out, "#res\taa\tchain");
  for(k=0; k<icoup; k++){
    fprintf(file_out, "\t%s_p", name_p[k]);
  }
  fprintf(file_out, "\n");
  for(i=0; i<Na; i++){
    site[i]=0;
    if(atom_site[i]==0)continue;
    if(atom_site[i]<0){
      int r=-atom_site[i];
      printf("WARNING, reference atom of res. %d  %c%s%c does not exist",
	     r, seq[r].amm, seq[r].pdbres, seq[r].chain);
      continue;
    }
    site[i]=1;
    //fprintf(file_out, "%s\t%c\t%c", pdbres[i], amm[i], chain[i]);
    fprintf(file_out, "%s\t%c\t%c", seq[i].pdbres, seq[i].amm, seq[i].chain);
    for(k=0; k<icoup; k++){
      fprintf(file_out, "\t%.3g", prof[k][i]); ///ave_p[k]
    }
    fprintf(file_out, "\n");
  }
  fclose(file_out);

  // Test pairs of sites
  if(N_sites<2){
    for(j=0; j<numran; j++)free(site_res_ran[j]);
    return;
  }
  sprintf(nameout, "%s.pairs_sites.dat", nameout1);
  printf("Writing %s\n", nameout);
  file_out=fopen(nameout, "w");
  fprintf(file_out, "# read binding sites in file %s\n", file_site);
  fprintf(file_out, "# profile type %c\n", PROF_TYPE);
  fprintf(file_out, "# N= %d, %d random sites\n", Na, numran);
  fprintf(file_out, "###1=nres1\t2=nres2\t3=overlap"); j=4;
  for(k=0; k<icoup; k++){
    fprintf(file_out, "\t%d=<%s_coup>\t%d=Z(%s_coup)\t%d=P>(%s_coup)",
	    j, name_p[k], j+1, name_p[k], j+2, name_p[k]); j+=3;
  }
  fprintf(file_out, "\n");

  int coupled_site[Nres];
  char chain_site[Nres];
  for(i=0; i<N_sites; i++){
    int n1=nres_site[i]; if(n1<=0)continue;
    int *isite=site_res[i];
    for(k=0; k<n1; k++){
      coupled_site[k]=isite[k]; chain_site[k]=chain_res[i][k];
    }
    for(j=i+1; j<N_sites; j++){
      int n2=nres_site[j]; if(n2<=0)continue;
      int *jsite=site_res[j], l=n1;
      for(k=0; k<n2; k++){
	coupled_site[l]=jsite[k]; chain_site[l]=chain_res[j][k]; l++;
      }
      int nchain=Count_chains(ichain, min_chain, max_chain, chain_lab,
			      chain_site, coupled_site, l);
      Identify_chains(all_chain, chain_lab, nchain, lab_chain, Nchain);
      for(k=0; k<numran; k++){
	Random_sites(site_res_ran[k], l, coupled_site,
		     ichain, min_chain, max_chain, all_chain, nchain,
		     ini_chain, end_chain);
      }
      fprintf(file_out, "## %s\t%s", description[i], description[j]);
      if(Nchain>1)fprintf(file_out, "\t%d chains", nchain);
      fprintf(file_out, "\n");
      int q=Site_overlap(isite, n1, jsite, n2);
      fprintf(file_out, "%d\t%d\t%d", n1, n2, q);
      res=site_results;
      for(k=0; k<icoup; k++){
	Test_pairs(res, atom_site, prof[k], coup[k], Na,
		   coupled_site, site_res_ran, numran, n1, l, atomres);
	fprintf(file_out, "\t%.3g\t%.1f\t%.3f",
		res->coup, res->Z_c, res->P_c);
	res++;
      }
      fprintf(file_out, "\n");
    }
  }
  fclose(file_out);
  for(j=0; j<numran; j++)free(site_res_ran[j]);

}

void Coupling_chains(double **Chain_dir, double **Chain_coord,
		     double **Chain_def, double **Chain_cov,
		     double **norm_chain, char *ch, int Na,
		     char *chain_label, int Nchain, char *nameout1)
{
  printf("Printing average couplings of %d chains\n", Nchain);


  char nameout[200];
  //sprintf(nameout, "%s.chain%c_couplings.dat", nameout1,chain_label[i]);
  sprintf(nameout, "%s.chain_couplings.dat", nameout1);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);

  char header[200]; sprintf(header, "#chain1 chain2 ");
  if(Chain_dir)strcat(header,  " Directionality");
  if(Chain_coord)strcat(header," Coordination");
  if(Chain_cov)strcat(header,  " Covariance");
  if(Chain_def)strcat(header,  " Deformation");
  strcat(header,  " Proximity_weight");
  fprintf(file_out, "%s\n", header);

  for(int i=0; i<Nchain; i++){
    for(int j=i; j<Nchain; j++){
      fprintf(file_out, "%c %c ", chain_label[i], chain_label[j]);
      if(Chain_dir)fprintf(file_out, " %.4g", Chain_dir[i][j]);
      if(Chain_coord)fprintf(file_out, "\t%.4g", Chain_coord[i][j]);
      if(Chain_cov)fprintf(file_out, "\t%.4g", Chain_cov[i][j]);
      if(Chain_def)fprintf(file_out, "\t%.4g", Chain_def[i][j]);
      fprintf(file_out, "\t%.4g", norm_chain[i][j]);
      fprintf(file_out, "\n");
    }
  }
  fclose(file_out);
}

double **Chain_coup(float **Coupling, double **norm, float **Proximity,
		    int Nchain, char *ch, int Na)
{
  if(Coupling==NULL)return(NULL);
  double **Chain_coup=Allocate_mat2_d(Nchain, Nchain);
  int ichain, jchain;
  for(ichain=0; ichain<Nchain; ichain++){
    for(jchain=0; jchain<=ichain; jchain++)norm[ichain][jchain]=0;
  }

  char ch_old_i=ch[0]; ichain=0;
  for(int ia=0; ia<Na; ia++){
    if(ch[ia]!=ch_old_i){ch_old_i=ch[ia]; ichain++;}
    char ch_old_j=ch[0]; jchain=0;
    for(int ja=0; ja<ia; ja++){
      if(ch[ja]!=ch_old_j){ch_old_j=ch[ja]; jchain++;}
      float w=Proximity[ia][ja];
      Chain_coup[ichain][jchain]+=w*Coupling[ia][ja];
      norm[ichain][jchain]+=w;
    }
  }
  for(ichain=0; ichain<Nchain; ichain++){
    for(jchain=0; jchain<=ichain; jchain++){
      Chain_coup[ichain][jchain]/=norm[ichain][jchain];
      if(jchain!=ichain){
	Chain_coup[jchain][ichain]=Chain_coup[ichain][jchain];
	norm[jchain][ichain]=norm[ichain][jchain];
      }
    }
  }
  return(Chain_coup);
}

void Print_interface_coupling(struct interaction *Int_list_inter,
			      int N_int_inter,
			      struct interaction *Int_list, int N_int,
			      atom *atoms, int Nchain, char *nameout,
			      float **Coup_dir, float**Coup_coord,
			      float **Coup_deformation, 
			      double **Chain_dir, double **Chain_coord,
			      double **Chain_def, double **norm_chain, int Na)
{
  char name2[200]; sprintf(name2, "%s.interface.dat", nameout);
  FILE *file_out=fopen(name2, "w");
  printf("Writing %s\n", name2);

  int ia, ja, n, npair=N_int_inter; float y;
  double Ave_dir=0, Ave_coord=0, Ave_def=0;
  double Min_dir=100, Min_coord=100, Max_def=0;
  int mpair=0, mpair_max=Nchain*(Nchain-1)/2, m;
  int c1[mpair_max], c2[mpair_max];
  for(n=0; n<N_int_inter; n++){
    ia=(atoms+Int_list_inter[n].i1)->res;
    ja=(atoms+Int_list_inter[n].i2)->res;
    if((ia>=Na)||(ja>=Na)){npair--; continue;}
    {
      int ic1=(atoms+Int_list_inter[n].i1)->chain;
      int ic2=(atoms+Int_list_inter[n].i2)->chain;
      if(ic1!=ic2){
	if(ic1>ic2){int tmp=ic2; ic2=ic1; ic1=tmp;}
	for(m=0; m<mpair; m++)if(ic1==c1[m] && ic2==c2[m])break;
	if(m==mpair){c1[m]=ic1; c2[m]=ic2; mpair++;}
      }
    }
    if(Coup_dir){
      y=Coup_dir[ia][ja]; Ave_dir+=y; if(y<Min_dir)Min_dir=y;
    }
    if(Coup_coord){
      y=Coup_coord[ia][ja]; Ave_coord+=y; if(y<Min_coord)Min_coord=y;
    }
    if(Coup_deformation){
      y=Coup_deformation[ia][ja]; Ave_def+=y; if(y>Max_def)Max_def=y;
    }
  }
  if(npair){Ave_dir/=npair; Ave_coord/=npair; Ave_def/=npair;}

  int numran=200;
  unsigned long iran=randomgenerator();
  if(ini_ran==0){
    InitRandom((RANDOMTYPE)iran); ini_ran=1;
  }
  double Ave_dir_ran=0, Ave_coord_ran=0, Ave_def_ran=0;
  double p_dir_ran=0, p_coord_ran=0, p_def_ran=0;
  for(int kran=0; kran<numran; kran++){
    npair=N_int_inter;
    double Ave_dir_tmp=0, Ave_coord_tmp=0, Ave_def_tmp=0;
    for(int k=0; k<N_int_inter; k++){
      n=RandomFloating()*N_int;
      if(n>=N_int){
	printf("WARNING in RandomFloating\n"); n=N_int-1;
      }
      ia=(atoms+Int_list[n].i1)->res;
      ja=(atoms+Int_list[n].i2)->res;
      if((ia>=Na)||(ja>=Na)){npair--; continue;}
      if(Coup_dir)Ave_dir_tmp+=Coup_dir[ia][ja];
      if(Coup_coord)Ave_coord_tmp+=Coup_coord[ia][ja];
      if(Coup_deformation)Ave_def_tmp+=Coup_deformation[ia][ja];
    }
    Ave_dir_tmp/=npair; Ave_dir_ran+=Ave_dir_tmp;
    if(Ave_dir_tmp>Ave_dir)p_dir_ran++;
    Ave_coord_tmp/=npair; Ave_coord_ran+=Ave_coord_tmp;
    if(Ave_coord_tmp>Ave_coord)p_coord_ran++;
    Ave_def_tmp/=npair; Ave_def_ran+=Ave_def_tmp;
    if(Ave_def_tmp<Ave_def)p_def_ran++;
  }
  
  fprintf(file_out,"# %d interface contacts\n", N_int_inter);
  fprintf(file_out, "Directionality_coupling: "
	  "Ave= %.4g Min= %.4g Random= %.2g p> = %.3g\n",
	  Ave_dir, Min_dir, Ave_dir_ran/numran, p_dir_ran/numran);
  fprintf(file_out, "Coordination_coupling:   "
	  "Ave= %.4g Min= %.4g Random= %.2g p> = %.3g\n",
	  Ave_coord, Min_coord, Ave_coord_ran/numran, p_coord_ran/numran);
  fprintf(file_out, "Deformation_coupling:    "
	  "Ave= %.4g Max= %.4g Random= %.2g p< = %.3g\n",
	  Ave_def, Max_def, Ave_def_ran/numran, p_def_ran/numran);

  double Ave_chain_coord=0, Ave_chain_dir=0, norm=0;
  fprintf(file_out, "# %d pairs of chains: ", mpair);
  for(m=0; m<mpair; m++){
    int ic1= c1[m], ic2= c2[m];
    fprintf(file_out, " %d-%d", ic1, ic2);
    Ave_chain_coord+=Chain_coord[ic1][ic2]*norm_chain[ic1][ic2];
    Ave_chain_dir+=Chain_dir[ic1][ic2]*norm_chain[ic1][ic2];
    norm+=norm_chain[ic1][ic2];
  }
  fprintf(file_out, "\n");
  Ave_chain_coord/=norm; Ave_chain_dir/=norm; 
  
  float Chain_dir_coeff=0, Chain_coord_coeff=0, Ave_coord_coeff=0,
    Ave_dir_coeff=0, bias=0, G_NM_coeff=0, DG_NM_coeff=0, DG_cont_coeff=0;
  if(strcmp(DOF_LABEL,"PHI")==0){ 
    // Param: KAPPA=5.0, SCALE_COORD=2, SCALE_D=2 Err=0.57 bias=0.37
    G_NM_coeff=0.0; DG_NM_coeff=31.3; DG_cont_coeff=0.27;
    Ave_coord_coeff=-14.2; Ave_dir_coeff=-3.61;
    Chain_coord_coeff=-13.3; Chain_dir_coeff=-4.16; bias=0.37;
  }else if(strcmp(DOF_LABEL,"PHIPSI")==0){
    // Param: KAPPA=8.2, SCALE_COORD=2, SCALE_D=2 Err=0.59 bias=0.44
    G_NM_coeff=0; DG_NM_coeff=26.2; DG_cont_coeff=0.325;
    Ave_coord_coeff=-13.8; Ave_dir_coeff=-3.7;
    Chain_coord_coeff=-13.8; Chain_dir_coeff=-4.58; bias=0.44;
  }else if(strcmp(DOF_LABEL,"PHIPSIOME")==0){
    // Param: KAPPA=15, SCALE_COORD=2, SCALE_D=2 Err=0.63 bias=0.50
    G_NM_coeff=0; DG_NM_coeff=21.2; DG_cont_coeff=0.335;
    Ave_coord_coeff=-12.5; Ave_dir_coeff=-3.64;
    Chain_coord_coeff=-13.0; Chain_dir_coeff=-4.66; bias=0.50;
  }else if(strcmp(DOF_LABEL,"PHIPSIOMESCH")==0){
    // Param: KAPPA=14, SCALE_COORD=2, SCALE_D=2 Err=0.60 bias=24.7
    G_NM_coeff=1.08; DG_NM_coeff=0; DG_cont_coeff=0.255;
    Ave_coord_coeff=-10.8; Ave_dir_coeff=-5.22;
    Chain_coord_coeff=-10.5; Chain_dir_coeff=-5.16; bias=24.7;
  }
  double DG_bind=
    G_NM_coeff*G_NM_holo+
    DG_NM_coeff*DG_NM_internal+
    DG_cont_coeff*E_cont_bind+
    Ave_coord_coeff*Ave_coord+
    Ave_dir_coeff*Ave_dir+
    Chain_coord_coeff*Ave_chain_coord+
    Chain_dir_coeff*Ave_chain_dir+
    bias;
  fprintf(file_out,"# Predicted binding free energy: %.4g\n", DG_bind);
  fprintf(file_out,"# Calculated as: %.4g*G_NM + %.4g*DG_NM + %.4g*DE_cont"
	  " + %.4g*Ave_coord_interface + %.4g*Ave_dir_interface"
	  " + %.4g*Ave_coord_CHAIN + %.4g*Ave_dir_chain + %.4g\n",
	  G_NM_coeff, DG_NM_coeff, DG_cont_coeff,
	  Ave_coord_coeff, Ave_dir_coeff, Chain_coord_coeff, Chain_dir_coeff,
	  bias);
  fprintf(file_out, "# Coordination computed as Coord_ij="
	  "exp(-%.3g*<(dij-<dij>)^2)\n", SCALE_COORD);
  fprintf(file_out, "# Interchain contacts weighted with "
	  "Proximity_ij=exp(-d_ij/%.3g)\n", SCALE_D);
  fprintf(file_out, "# G_NM/nres  = %.4g\n", G_NM_holo);
  fprintf(file_out, "# DG_NM_internal/nres = %.4g\n", DG_NM_internal);
  fprintf(file_out, "# DG_NM_rigid/nres = %.4g\n", DG_NM_rigid);
  //fprintf(file_out, "# DG_NM/nres = %.4g\n", DG_NM_rigid+DG_NM_internal);
  fprintf(file_out, "# DE_cont = %.4g\n", E_cont_bind);
  fprintf(file_out, "# Ave_inter_coord = %.4g\n", Ave_coord);
  fprintf(file_out, "# Ave_inter_dir   = %.4g\n", Ave_dir);
  fprintf(file_out, "# Ave_chain_coord = %.4g\n", Ave_chain_coord);
  fprintf(file_out, "# Ave_chain_dir   = %.4g\n", Ave_chain_dir);

  fclose(file_out);
  }

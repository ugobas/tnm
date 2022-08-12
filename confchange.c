#include "coord.h"
//#include "Parameters.h"
#include "tnm.h"
#include "interactions_tnm.h"
#include "contacts.h"
#include "nma.h"
#include "vector.h"
#include "McLachlan.h"
#include "align_tnm.h"
#include "output_tnm.h"
#include "buildup.h"
#include "simulation.h"
#include "kinetic_tnm.h"
#include "random3.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int ALL_TYPES=0;
/* If ALL_TYPES, examine all ways to determine the torsion angles, namely: Diff_Tors, OLS fit, RRR CV, RRR MP, Greedy. If ALL_TYPES=0, adopt only RRR Cv.
 */

#define TEST_NULL 0  // Test best parameter of null model

//#define DIFF_THR 0.25
#define COLL_THR 2   // Select if exp(-S_cart) > COLL_THR

/*extern float Torsional_superimposition(double *x1, double *x2,
				       int N_axes, struct Reference Ref,
				       struct Jacobian J);
extern Test_Diff(struct Tors *Diff, struct Reference Ref, struct Jacobian J);
*/
struct Jacobian *J_CC=NULL;
extern float Freq_unit;

void Print_diff(int *num_dphi, char **name_dphi,
		float **diff_angles, float **MW_diff_angles,
		float **c2_alpha, int *outlier_tors,
		struct Tors *Diff,
		struct Jacobian *J,
		struct Normal_Mode *NM, float M_sqrt, 
		struct bond *bonds, atom *atoms1, int natoms1,
		struct Reference Ref, float *coord_old, float *coord_target,
		struct residue *seq,
		struct axe *axe,
		int nstruct, char *name,
		FILE *file_out, char *name_pdb,
		int anharmonic, char *out, int whichmode);
void Print_angular_differences(float **diff_angles, char *type,
			       char **name_dphi, int num_dphi,
			       char *name_para, struct axe *axe,
			       int N_axes, struct residue *seq1);

int Mode_confchange_fit(struct Tors *Diff, struct Normal_Mode *NM,
			struct bond *bonds, atom *atoms, int natoms,
			struct Reference Ref, struct residue *seq,
			float *coord_1, float *coord_2, float M_sqrt,
			char *nameout, char *name_pdb, int type, int NMODES);
int Mode_confchange_old(struct Tors *Diff, struct Normal_Mode *NM,
			struct bond *bonds, atom *atoms, int natoms,
			struct Reference Ref, struct residue *seq,
			float *coord_1, float *coord_2, float M_sqrt,
			char *nameout, char *name_pdb);
float RMSD_confchange(float cc, float *mode, float *d_phi,
		      float *coord_all, float *coord_new,
		      float *d_phi_opt, int naxes, struct bond *bonds,
		      atom *atoms, int natoms, struct Reference Ref,
		      float *coord_2);
int Lasso_confchange_fit(struct Tors *Diff, struct Normal_Mode *NM,
			 struct bond *bonds, atom *atoms1, int natoms1,
			 struct Reference Ref, float *coord_1, float *coord_2,
			 char *nameout);
float RMSD_Lasso_confchange(float lambda, float *cnew, float *d_phi,
			    float *coord_all, float *coord_new,
			    float *c, double c2_sum,
			    struct Normal_Mode *NM, struct bond *bonds,
			    atom *atoms, int natoms,
			    struct Reference Ref, float *coord_2);
extern float Find_max_quad(float x1, float x2, float x3,
			   float y1, float y2, float y3,
			   float MIN, float MAX);
float Fermi_function(float c2, float c2_thr, float S);
float Test_Renyi(float *ene, float *P, int N, char *name, FILE *file);
float Test_Null(float *omega_minus2, float *P, int N, char *name, FILE *file);
int Select_cc2(float *dx2, int nres1, float *Dcart, atom *atoms1,
	       struct Reference Ref, char *name);
int Periodic_angles(float *dphi, struct axe *axe, int n);
static float RMS(float *v, int n, int *outlier);
unsigned long randomgenerator(void);
static int ini_ran=0;
static long idum;

struct Normal_Mode NM_kin;
extern void d_Diagonalize(int N, double **MATRIX, float *eigen_values,
			  float **eigen_vector, int SIGN);


/************************************************************************/

float Examine_confchange(struct Tors *Diff,
			 struct bond *bonds, struct axe *axe,
			 struct chain *chains, int Nchain,
			 char *namesumm, char *name1, char *namepdb,
			 atom *atoms1, float *coord_1,
			 struct residue *seq1, int nres1,
			 int N_diso1, int natoms1,
			 atom *atoms2, float *coord_2,
			 struct residue *seq2, int nres2,
			 int N_diso2, int natoms2,
			 struct Reference Ref1, //struct Jacobian *J,
			 struct ali_atoms ali_a,
			 struct interaction *Int_list1, int N_int1,
			 char *INT_TYPE, float s0,
			 int N_modes,  struct Normal_Mode NM,
			 int nprint, struct Ali_score ali, 
			 struct Para_simul Para_simul,
			 //float *diff_phi, float *greedy_phi,
			 int nstruct, float *B_TNM,
			 int *outlier_tors,
			 char *nameout, int PRINT_CHANGE,
			 int anharmonic)
// nstruct = 0 if experimental confchange, = k if simulated structures
{
  char namenew[100]; float *sigma2;
  if(anharmonic==0){
    sprintf(namenew, "%s", nameout);
    sigma2=NM.sigma2;
  }else{
    sprintf(namenew, "%s.anharmonic", nameout);
    sigma2=NM.sigma2_anhar;
  }

  /*
  // Prepare reference atoms
  printf("Preparing to examine conformation change\n");
  if(nstruct==0){
    // Set reference atoms
    Set_reference(&Ref1, 0, "EB", atoms1, 0, natoms1);
    Align_references(Ref1, atoms1);
  }else if(Ref1.atom_num==NULL){
    Set_reference(&Ref1, 0, "EB", atoms1, 0, natoms1);
  }
  // Compute RMSD
  int N_Cart=Ref1.N_Cart; float coord_1[N_Cart], coord_2[N_Cart];
  Write_ref_coord_atom(coord_1, Ref1.N_ref, atoms1, Ref1.atom_num);
  Write_ref_coord_atom(coord_2, Ref1.N_ref, atoms2, Ref2.atom_num);
  */

  int N_Cart=ali_a.N_Cart, i;
  float rmsd=rmsd_mclachlan_f(coord_1, coord_2, ali_a.mass, ali_a.N_ref);
  for (i=0;i<N_Cart;i++)Diff->Cart[i]=coord_2[i]-coord_1[i];

  printf("Examining confchange, RMSD= %.2f\n", rmsd);
  printf("Changing reference atoms from %d (kin) to %d (CC)\n",
	 Ref1.N_Cart/3, ali_a.N_ref);

  // Jacobian matrix
  int naxe=NM.N_axes;
  if(J_CC==NULL){
    J_CC=malloc(1*sizeof(struct Jacobian));
    Allocate_Jacobian(J_CC, naxe, Ref1.N_Cart);
    int first_kin[naxe], last_kin[naxe];
    Change_kin(first_kin, last_kin, axe, naxe, natoms1, Ref1);
    J_CC->N_kin=Compute_kinetic(J_CC, axe, naxe, atoms1, natoms1,Ref1,1);
    Change_kin_back(axe, naxe, first_kin, last_kin);
  }
  // Mass modes
  int MMASS=0;
  if(anharmonic==0){
    Allocate_Normal_modes(&NM_kin, NM.N, NM.N_axes, NM.N_Cart);
    d_Diagonalize(NM.N_axes, J_CC->T, NM_kin.omega2, NM_kin.Tors, 1);
    FILE *file_mode=NULL;
    if(MMASS){
      char out[100]; sprintf(out, "%s.mass_modes.dat", nameout);
      file_mode=fopen(out, "w"); printf("Writing %s\n", out);
      fprintf(file_mode,"# Eigenvectors of the mass matrix\n");
      fprintf(file_mode,"# n 1/<K> Coll_cart Coll_tors most_sim_nor_mod c2\n");
    }
    for(i=0; i<NM.N; i++){
      NM_kin.select[i]=1;
      NM_kin.sigma2[i]=NM_kin.omega2[i];
      NM_kin.omega[i]=1./sqrt(NM_kin.omega2[i]);
      // Normalize so that Tors is orthonormal wrt kinetic energy
      float norm=NM_kin.omega[i], *mode=NM_kin.Tors[i];
      for(int a=0; a<naxe; a++)mode[a]*=norm;
      Compute_MW(NM_kin.MW_Tors[i], NM_kin.Tors[i], J_CC);
      Convert_torsion2cart(NM_kin.Cart[i], atoms1, NM_kin.Tors[i],
			   axe, naxe, Ref1, i);
      NM_kin.Tors_coll[i]=Collectivity_norm2(NM_kin.Tors[i], naxe)/naxe ;
      NM_kin.MW_Tors_coll[i]=Collectivity_norm2(NM_kin.MW_Tors[i],naxe)/naxe;
      NM_kin.Cart_coll[i]=
	Collectivity_norm2(NM_kin.Cart[i], NM_kin.N_Cart)/N_Cart;
      //Predict_fluctuations(&NM_kin, NM_kin.sigma2);
      if(MMASS){
	double cmax=0, c; int imax=0, j; float *mode=NM_kin.Tors[i];
	double sigma2=0;
	for(j=0; j<NM.N; j++){
	  c=Scalar_product(mode, NM.Tors[j], NM.N_axes); c*=c;
	  sigma2+=c*NM.omega2[j];
	  if(c>cmax){cmax=c; imax=j;}
	} 
	NM_kin.sigma2[i]=1./(sigma2*NM_kin.omega2[i]*NM_kin.omega2[i]);
	if(NM_kin.Cart_coll[i]<0.01)NM_kin.sigma2[i]=0.1;
	fprintf(file_mode, "%d\t%.3g\t%.3f\t%.3f\t%d\t%.3g\n",
		i,NM_kin.sigma2[i],NM_kin.Cart_coll[i],NM_kin.Tors_coll[i],
		imax,cmax);
      }
    }
    if(MMASS)fclose(file_mode);
  }


  // Open file for printing, one chain
  FILE *file_out;
  if(nprint){
    file_out=fopen(namesumm, "a");
  }else{
    file_out=fopen(namesumm, "w");
    printf("Writing %s\n", namesumm);
    if(anharmonic==0){
      fprintf(file_out, "protein               %s\n", name1);
      fprintf(file_out, "Nres                  %d\n", nres1);
      fprintf(file_out, "Degrees of freedom    %d\n", N_modes);
      fprintf(file_out, "Disordered gaps       %d\n", N_diso1);
      fprintf(file_out, "#\n");
    }
  }

  // Printing sequence comparison
  if((nstruct==0)&&(anharmonic==0)){
    fprintf(file_out, "# Conformation change\n");
    fprintf(file_out, "seq_id(percent)       %.1f\n", ali.seq_id);
    fprintf(file_out, "align                 %.3f\n", ali.aligned);
    fprintf(file_out, "ngaps                 %d\n", ali.ngaps);
    if(ali.mammoth){
      fprintf(file_out,"psi                   %.1f\n", ali.psi);
      fprintf(file_out,"Mammoth score         %.1f\n", ali.mamm_score);
    }
  }

  // Contact matrix and contact energy
  int N_int2=0; struct interaction *Int_list2;
  Compute_interactions(&N_int2, &Int_list2, INT_TYPE,
		       atoms2, natoms2, nres2, namenew);
  float Cont_overlap=
    Contact_overlap(Int_list1, N_int1, Int_list2, N_int2, ali.alignres);
  float E1=Contact_energy(Int_list1, N_int1, seq1);
  float E2=Contact_energy(Int_list2, N_int2, seq2);
  float DE=E2-E1-(nres1-nres2)*s0;
  free(Int_list2);

  // Printing properties of conformation change
  if(nstruct){
    fprintf(file_out, "Structure             %d\n", nstruct);
  }else if(anharmonic==0){
    fprintf(file_out, "Disordered axis2      %d\n", N_diso2);
  }
  if(anharmonic==0){
    fprintf(file_out, "Contact_overlap       %.3f\n", Cont_overlap);
    fprintf(file_out, "(ECont2-ECont1)       %.3f\n", DE);
    // float rmsd=
    //Torsional_superimposition(atom_str1,atom_str2,NM.N_axes,Ref,J);
    fprintf(file_out, "RMSD %s atoms:      %.2f\n", REF, rmsd);
  }

  if(B_TNM){
    // Compute correlation between Cartesian displacements and conf. change
    float dx2[nres1];
    int n=Select_cc2(dx2, nres1, Diff->Cart, atoms1, Ref1, "CA");
    if(n){
      float slope, off, r=Corr_coeff(B_TNM, dx2, n, &slope, &off);
      fprintf(file_out, "r(therm,change)_Cart_");
      if(anharmonic){fprintf(file_out, "anhar");}
      else{fprintf(file_out, "har  ");}
      fprintf(file_out, " %.3f\n",r);
    }
  }

  /************************* Torsional fits ****************************/
  char name_pdb[100];
  if(nstruct==0){
    sprintf(name_pdb, "%s.Confchange.pdb", nameout);
    FILE *file_pdb=fopen(name_pdb, "w"); fclose(file_pdb); 
  }

  float Lambda=0; //thr=Para_confchange.RMSD_THR, 
  double M_tot=0, M_sqrt;
  for(i=0; i<Ref1.N_Cart; i+=3)M_tot+=Ref1.mass_coord[i];
  M_sqrt=sqrt(M_tot);

  // Different methods to obtain diff_phi
  int max_dphi=5; int istore=2, jstore=-1; // select RRR Cv
  if((nstruct)||(PRINT_CHANGE==0))max_dphi=4; // Do not generate confchange
  int num_dphi=0;
  float *diff_angles[max_dphi], *MW_diff_angles[max_dphi], *c2_alpha[max_dphi];
  char *name_dphi[max_dphi], name[80];
  float coeff[NM.N]; 

  for(int i_dphi=0; i_dphi<max_dphi; i_dphi++){ //

    if(nstruct)continue; 
    if((ALL_TYPES==0)&&((i_dphi==0)||(i_dphi==3)))continue;

    if(i_dphi==1){
      continue;
      /*sprintf(name, "Lasso regression");
      Lasso_confchange_fit(Diff, &NM, bonds, atoms1, natoms1, Ref1,
      coord_1, coord_2, name1);*/
      sprintf(name, "Modes fit"); 
      int type=0; // rank by frequency (0) or by overlap (1)?
      int NMODES=10; //10;
      printf("Lasso regression of torsion angle changes\n");
      Mode_confchange_fit(Diff, &NM, bonds, atoms1, natoms1, Ref1, seq1,
			  coord_1, coord_2, M_sqrt, name1, name_pdb,
			  type, NMODES);
      continue;
    }else if(i_dphi<3){
      char type='X'; // Rescaled ridge regression fits
      if(i_dphi==0){type='O';} // Ordinary least squares
      else if(i_dphi==1){type='M'; continue;} // Substituted with Lasso
      else if(i_dphi==2){type='C';}
      sprintf(name, "Fit");
      if(i_dphi)strcat(name, " RRR");
      char tmp[40];
      sprintf(tmp, " type %c", type); strcat(name, tmp);
      printf("Ridge regression type %c of torsion angle changes\n",type);
      if(Convert_cart2torsion_fit(Diff, Ref1, J_CC, name1, type, &Lambda)<0){
	fprintf(file_out, "WARNING, ridge regression has failed\n");
	printf("WARNING, ridge regression has failed\n"); continue;
      }

    }else if(i_dphi==3){
      // Difference of internal coordinates
      sprintf(name, "Differences of torsion angles different bonds");
      struct bond *bonds2=Set_bonds_topology(natoms2, atoms2, seq2);
      Set_bonds_measure(bonds, natoms1, atoms1);
      Set_bonds_measure(bonds2,natoms2, atoms2);
      Change_internal(Diff->Tors, naxe, bonds, bonds2, natoms1, 0, 0);
      free(bonds2);
    }else  if(i_dphi==4){
      fflush(file_out);
      printf("Interpolating torsional conformation change\n");
      sprintf(name, "Differences of torsion angles, iterative");
      double nbtops = CLOCKS_PER_SEC, t0=clock();       // time
      // Torsional_confchange_RRR
      int NMODES=NM.N; if(NMODES>NM.N)NMODES=NM.N;
      struct Normal_Mode *Modes; char mode_name[20];
      float cc_thr=0.0; // Filter modes contributing < cc_thr A
      if(1||MMASS==0){Modes=&NM; strcpy(mode_name, "normal");}
      else{Modes=&NM_kin; strcpy(mode_name, "mass");}
      Torsional_confchange(Diff->Tors, "greedy", bonds, J_CC, Ref1, ali_a,
			   0.1, atoms1, natoms1, namenew, axe, naxe,seq1,
			   nres1, chains, Nchain, atoms2, natoms2, 
			   Para_simul, Modes, mode_name, NMODES, cc_thr);
      printf("\nMaking conformation change. Time= %.2lf sec.\n",
	     (clock()-t0)/nbtops); t0=clock();
    }
    fprintf(file_out, "### Projection on normal modes");
    Print_diff(&num_dphi, name_dphi, diff_angles, MW_diff_angles, c2_alpha,
	       outlier_tors, Diff, J_CC, &NM, M_sqrt,
	       bonds, atoms1, natoms1, Ref1, coord_1, coord_2,
	       seq1, axe, nstruct, name, file_out, name_pdb,
	       anharmonic, nameout, 0);
    if(MMASS && anharmonic==0){
      fprintf(file_out, "### Projection on mass modes");
      Print_diff(&num_dphi, name_dphi, diff_angles, MW_diff_angles, c2_alpha,
		 outlier_tors, Diff, J_CC, &NM_kin, M_sqrt,
		 bonds, atoms1, natoms1, Ref1, coord_1, coord_2,
		 seq1, axe, nstruct, name, file_out, name_pdb,
		 anharmonic, nameout, 1);
    }
    if(i_dphi==istore){
      jstore=num_dphi-1;
      for(i=0; i<NM.N; i++)coeff[i]=Diff->coeff[i];
    }

  }
  fclose(file_out);

  if(jstore>=0){
    printf("Storing confchange method %s\n", name_dphi[jstore]);
    for(i=0; i<NM.N; i++){
      Diff->coeff[i]=coeff[i];
      NM.confchange2[i]=c2_alpha[jstore][i];
    }
    for(i=0; i<naxe; i++)Diff->Tors[i]=diff_angles[jstore][i];

  }

  if(nstruct==0){

    // Print list of angles and correlations with normal modes
    if(nameout && num_dphi>1){
      /*Print_angular_differences(diff_angles, "", name_dphi, num_dphi,
				namenew, axe, NM.N_axes, seq1);
      Print_angular_differences(MW_diff_angles, "MW", name_dphi, num_dphi,
      namenew, axe, NM.N_axes, seq1);*/
      
      char nameout2[200]; int a;
      sprintf(nameout2, "%s.Mode_projections.dat", namenew);
      printf("Writing mode projections in %s\n", nameout2);
      file_out=fopen(nameout2, "w");
      fprintf(file_out, "# %d modes %d methods for deriving angles\n",
	      NM.N, num_dphi);
      for(i=0; i<num_dphi; i++)
	fprintf(file_out, "# Method %d: %s\n", i+1, name_dphi[i]);
      fprintf(file_out, "#mode ");
      for(i=0; i<num_dphi; i++)fprintf(file_out, " %d=dphi%d", i+2, i+1);
      fprintf(file_out, " %d=1/omega^2", i+2);
      if(anharmonic==0)fprintf(file_out, " %d=mass^2", i+3);
      fprintf(file_out, "\n");
      for(a=0; a<NM.N; a++){
	fprintf(file_out, "%d", a);
	for(i=0; i<num_dphi; i++)fprintf(file_out, "\t%.5f", c2_alpha[i][a]);
	fprintf(file_out, "\t%.3g", sigma2[a]);
	if(anharmonic==0)fprintf(file_out, "\t%.3g", NM_kin.sigma2[a]);
	fprintf(file_out, "\n");
      }
      fclose(file_out);
    }

    for(i=0; i<num_dphi; i++){
      free(c2_alpha[i]);
      free(diff_angles[i]);
      free(MW_diff_angles[i]); //modyves: was never freed
      free(name_dphi[i]);
    }
  }

  Empty_Jacobian(*J_CC);
  // should probably not depend on nstruct==0, since alloc doesn't
  free(J_CC); J_CC=NULL;
  if(anharmonic)Empty_Normal_modes(NM_kin);

  return(rmsd);
}

float RMSD_CA(float *coord_1, float *coord_2, struct Reference Ref,
	      atom *atoms, int nres, struct ali_atoms ali_atoms)
{
  int n3=3*nres; float coord_CA1[n3], coord_CA2[n3], mass[nres];
  int k=0, n=0;
  for(int i=0; i<ali_atoms.N_ref; i++){
    if(strncmp(atoms[Ref.atom_num[i]].name, "CA", 2))continue;
    int ii=3*i; mass[n]=ali_atoms.mass[i]; n++;
    if(n>nres){printf("ERROR, too many CA atoms\n"); exit(8);}
    for(int j=0; j<3; j++){
      coord_CA1[k]=coord_1[ii]; 
      coord_CA2[k]=coord_2[ii];
      ii++; k++;
    }
  }
  float rmsd=rmsd_mclachlan_f(coord_CA1, coord_CA2, mass, n);
  return(rmsd);
}

void Print_diff_fluct(struct Tors *Diff, struct Normal_Mode *NM,
		      int *outlier_tors, char *nameout2)
{
  FILE *file_out; char nameout[200];

  sprintf(nameout, "%s_Tors_diff.dat", nameout2);
  file_out=fopen(nameout, "w"); printf("Writing %s\n", nameout);
  fprintf(file_out, "# Squared fluctuations of torsion angles\n");
  if(outlier_tors)fprintf(file_out, "# Excluding outliers\n");
  fprintf(file_out, "# Thermal Conf.change\n");
  int a, ma=0;
  for(a=0; a<NM->N_axes; a++){
    if(outlier_tors && (outlier_tors[a])){ma++; continue;}
    fprintf(file_out,"%.4g\t%.4g\n", NM->tors_fluct[a],
	    Diff->Tors[a]*Diff->Tors[a]);
  }
  if(outlier_tors)fprintf(file_out, "# %d angles excluded\n", ma);
  fclose(file_out);
  
 sprintf(nameout, "%s_Cart_diff.dat", nameout2);
  file_out=fopen(nameout, "w"); printf("Writing %s\n", nameout);
  fprintf(file_out, "# Squared fluctuations of Cartesian coordinates\n");
  fprintf(file_out, "# Thermal Conf.change\n");
  for(a=0; a<NM->N_Cart; a++){
    fprintf(file_out,"%.4g\t%.4g\n", NM->cart_fluct[a],
	    Diff->Cart[a]*Diff->Cart[a]);
  }
  fclose(file_out);

}

void Print_diff(int *num_dphi, char **name_dphi,
		float **diff_angles, float **MW_diff_angles,
		float **c2_alpha, int *outlier_tors,
		struct Tors *Diff,
		struct Jacobian *J,
		struct Normal_Mode *NM, float M_sqrt, //float thr, 
		struct bond *bonds, atom *atoms1, int natoms1,
		struct Reference Ref, float *coord_1, float *coord_2,
		struct residue *seq1, struct axe *axe,
		int nstruct, char *name,
		FILE *file_out, char *name_pdb,
		int anharmonic, char *out, int whichmode)
{
  // RMSD
  double norm_c2=
    Scalar_product_weighted(Diff->Cart,Diff->Cart,Ref.mass_coord,Ref.N_Cart);
  Diff->RMSD=sqrt(norm_c2)/M_sqrt; Diff->M=M_sqrt*M_sqrt;

  // Normal modes
  int naxes=NM->N_axes;
  int N_modes=NM->N, i, k, a;
  float *sigma2;
  if(anharmonic==0){
    sigma2=NM->sigma2;
    fprintf(file_out, "\n");
  }else{
    sigma2=NM->sigma2_anhar;
    fprintf(file_out, " anharmonic\n");
    Predict_fluctuations(NM, sigma2);
  }

  int REG2=1; // Regularization: 0=L1 (thredhold on rmsd), 1=L2 (limit DE)
  float c_thr=0, Lambda=0, rmsd_thr=0;
  if(REG2==0){
    rmsd_thr=0.2;  // Threshold on mode RMSD
    float fmin=0.05; if(rmsd_thr<Diff->RMSD*fmin)rmsd_thr=Diff->RMSD*fmin;
    c_thr=rmsd_thr*M_sqrt;
  }else{
    float Coeff=0.1;
    sprintf(name, "RR DE C=%.2g", Coeff);
    double sum_sigma2=0; // Regularization on DE=sum_a (c_a/w_a)^2
    for (i=0; i<N_modes; i++)sum_sigma2+=sigma2[i];
    Lambda=Coeff*sum_sigma2/N_modes;
  }

   Periodic_angles(Diff->Tors, axe, naxes); // Remove 2pi
   Compute_MW(Diff->MW_Tors, Diff->Tors, J);
   // Initialize angles, they will be recomputed
   for(a=0; a<naxes; a++)Diff->Tors[a]=0;

  // Print name
  name_dphi[*num_dphi]=malloc(80*sizeof(char));
  strcpy(name_dphi[*num_dphi], name);
  fprintf(file_out, "##### %s\n", name);

  double sum_c2=0, sum_sigma2=0, Delta_E=0;
  float *c2=NM->confchange2, c; int nmodes=0;
  for (i=0; i<N_modes; i++){
    if(sigma2[i]==0){c2[i]=0; continue;}
    c=Scalar_product(Diff->MW_Tors, NM->MW_Tors[i], NM->N_axes);
    if(REG2==0){
      if(fabs(c)<c_thr){c=0;}else{nmodes++;}
      /*}else{
      c=Scalar_product_weighted(Diff->Cart,NM->Cart[i],
				Ref.mass_coord,Ref.N_Cart);
				c/=(1+Lambda/sigma2[i]);*/
    }
    float *tt=NM->Tors[i];
    for(a=0; a<naxes; a++)Diff->Tors[a]+=c*tt[a];
    Diff->coeff[i]=c;
    c2[i]=c*c;
    float de=c2[i]/sigma2[i];
    if(REG2){float Lom2=Lambda/sigma2[i]; de/=((1+Lom2)*(1+Lom2));}
    Delta_E+=de; 
    if(fabs(c)>0)sum_sigma2+=sigma2[i]; // sum only selected modes
    sum_c2+=c2[i];
  }
  for(i=0; i<N_modes;i++)c2[i]/=sum_c2; //norm_c2; // Normalize
  Periodic_angles(Diff->Tors, axe, naxes); // Remove 2pi
  Compute_MW(Diff->MW_Tors, Diff->Tors, J);

  // Scale for equating random and real confchange
  double scale2=sum_c2/sum_sigma2;

  //fprintf(file_out, "# Total RMSD: %.3g %d modes\n", Diff->RMSD, N_modes);
  if(REG2){
    nmodes=N_modes; // L2 regularization, all modes are selected
    fprintf(file_out, "# L2 regularization:\n"
	    "# Mode coordinate c'_a = c_a/(1+Lambda*w_a^2)\n"
	    "# Lambda=M*RMSD^2/(n*A^2)=sum_a 1/w_a^2/n = %.3g\n", Lambda);
  }else{
    fprintf(file_out, "# L1 regularization:\n"
	    "# Eliminate modes that contribute < %.2g A: RMSD %.3g %d modes\n"
	    "# Reduction factor: RMSD^2 %.3f Modes %.3f\n",
	    rmsd_thr, sqrt(sum_c2)/M_sqrt, nmodes,
	    sum_c2/(norm_c2), nmodes/(float)N_modes);
    fprintf(file_out, "# Scale factor: %.3g (ini) ", sqrt(scale2));
    // Select modes with individual contribution > threshold
    int IT_MAX=20, iter;
    for(iter=0; iter<IT_MAX; iter++){
      float om_thr=sqrt(scale2)/c_thr;
      double sum_sigma2_new=0;
      for(i=0; i<N_modes; i++){
	sum_sigma2_new+=sigma2[i];
	if(sigma2[i] && NM->omega[i]>om_thr)break;
      }
      float scale2_new=sum_c2/sum_sigma2_new;
      if(fabs(scale2_new-scale2)<0.00001*scale2){
	scale2=scale2_new; break;
      }
      scale2=scale2_new;
    }
    fprintf(file_out, "%.3g (end) iter=%d\n", sqrt(scale2), iter);
    if(iter==IT_MAX){
      fprintf(file_out,
	      "# WARNING, could not set the scale self-consistently\n");
    }
    nmodes=0; float om_thr=3*sqrt(scale2)/c_thr;
    for(i=0; i<N_modes; i++){
      if(sigma2[i] && NM->omega[i]>om_thr){nmodes=i; break;}
    }
    //int nm_thr=10; if(nmodes<nm_thr)nmodes=nm_thr;
    char tmp[100];
    sprintf(tmp, "# Number of modes with max. thermal rmsd > %.2g: %d\n",
	    rmsd_thr, nmodes); 
    printf("%s", tmp); fprintf(file_out, "%s", tmp);
  }
 
  // Generate new conformation from all torsion angles
  float coord_all[3*natoms1], coord_new[3*Ref.N_ref];
  float RMSD=0, RMSD_target=0;
  Set_bonds_measure(bonds, natoms1, atoms1);
  int OPTIMIZE_AMPL=1;
  if(OPTIMIZE_AMPL){
    struct bond bonds_min[natoms1]; float f_min=1;
    Copy_bonds(bonds_min, bonds, natoms1);
    RMSD_target=Optimize_amplitude(bonds_min, &f_min, Diff->Tors, NM->N_axes,
				   bonds, natoms1, Ref, coord_2);
    printf("RMSD to target= %.2f\n", RMSD_target);
    Copy_bonds(bonds, bonds_min, natoms1);
  }else{
    Build_up(bonds, natoms1, Diff->Tors, NM->N_axes);
  }
  Put_coord(coord_all, bonds, natoms1);
  Write_ref_coord(coord_new, Ref.N_ref, coord_all, Ref.atom_num);
  RMSD_target=rmsd_mclachlan_f(coord_2, coord_new, Ref.mass_atom, Ref.N_ref);
  RMSD=rmsd_mclachlan_f(coord_1, coord_new, Ref.mass_atom, Ref.N_ref);
  printf("RMSD to target= %.2f\n", RMSD_target);
  float rmsd1_pdb=RMSD, rmsd2_pdb=RMSD_target;

  // Generate new conformation from linear approximation
  float RMSD_lin=0, RMSD_lin_target=0, RMSD_lin_full=0;
  if(1){
    RMSD_lin=RMSD_LIN(&RMSD_lin_target, &RMSD_lin_full, J, Diff->Tors,
		      Ref, coord_1, coord_2, coord_all);
  }


  // Conformation change
  int nn=0, n_outlier=0; float THR=0.5;
  for(a=0; a<NM->N_axes; a++){
    if(outlier_tors && (outlier_tors[a])){n_outlier++; continue;}
    if(fabs(Diff->Tors[a])>THR)nn++;
  }
  Diff->Tors_frac=Tors_fraction(Diff, Ref.mass_coord);
  Diff->RMSD_Tors=RMS(Diff->Tors, Diff->N_axes, outlier_tors);

  if((nstruct==0)&&(anharmonic==0)&&(whichmode==0)){
    fprintf(file_out, "RMSD_Cart             %.3f\n", Diff->RMSD);
    fprintf(file_out, "RMSD_WTors            %.3f\n", Diff->RMSD_W);
    fprintf(file_out, "Non-tors_RMSD         %.3f\n", Diff->RMSD_NoTors);
    fprintf(file_out, "RMSD_reconstruct_1    %.3f\n", RMSD);
    fprintf(file_out, "RMSD_reconstruct_2    %.3f\n", RMSD_target);
    fprintf(file_out, "RMSD_lin_1            %.3f\n", RMSD_lin);
    fprintf(file_out, "RMSD_lin_2            %.3f\n", RMSD_lin_target);
    fprintf(file_out, "RMSD_lin_full         %.3f\n", RMSD_lin_full);
    fprintf(file_out, "Tors_frac             %.3f\n", Diff->Tors_frac);
    if(n_outlier)
      fprintf(file_out, "Eliminating %d outlier angles\n", n_outlier);
    fprintf(file_out, "Angular RMSD          %.3f\n", Diff->RMSD_Tors);
    fprintf(file_out, "Num(dtheta>%.1f)       %d\n", THR, nn);

    float Coll=Collectivity_norm2(Diff->Cart, Ref.N_Cart)/Ref.N_Cart;
    fprintf(file_out, "Coll.(cart)           %.3f\n", Coll);
    Coll=Collectivity_norm2(Diff->MW_Tors, NM->N_axes)/NM->N_axes;
    fprintf(file_out, "Coll.(MWtors)         %.3f\n", Coll);
    //Coll=Collectivity_norm2(Diff->Tors, NM->N_axes)/NM->N_axes;
    Coll=Collectivity_norm2_outlier(Diff->Tors, NM->N_axes, outlier_tors);
    Coll/=(NM->N_axes-n_outlier);
    fprintf(file_out, "Coll.(tors)           %.3f\n", Coll);
    fprintf(file_out, "#\n");
  }

  //////////////////////////////////////////////////////
  // Print PDB
  if(name_pdb && anharmonic==0 && whichmode==0){
    FILE *file_pdb=fopen(name_pdb, "a");
    fprintf(file_pdb, "MODEL %d  %s\n", *num_dphi+1, name);
    fprintf(file_pdb, "REMARK  RMSD from target: %.2f ",rmsd2_pdb);
    fprintf(file_pdb, "RMSD from starting structure: %.2f\n",rmsd1_pdb);
    float coord_all_1[3*natoms1], *c=coord_all_1, mass[natoms1]; int i, j;
    for(i=0; i<natoms1; i++){
      mass[i]=1;
      double *r=atoms1[i].r; for(j=0; j<3; j++){*c=*r; c++; r++;}
    }
    float rmsd1=rmsd_mclachlan_f(coord_all_1, coord_all, mass, natoms1);
    Print_PDB(file_pdb, atoms1, natoms1, coord_all, seq1, *num_dphi, rmsd1);
    fclose(file_pdb);
  }

  // Correlation angular fluctuations - conformation change
  Diff->RMSD_Tors=RMS(Diff->Tors, naxes, outlier_tors);
  printf("RMSD(Tors)= %.3f (without outliers)\n", Diff->RMSD_Tors);
  float r, slope, offset;
  float dt[NM->N_axes], tors_fluct[NM->N_axes]; int ma=0;
  for(a=0; a<NM->N_axes; a++){
    if(outlier_tors && (outlier_tors[a]))continue;
    dt[ma]=Diff->Tors[a]*Diff->Tors[a];
    tors_fluct[ma]=NM->tors_fluct[a];
    ma++;
  }
  r=Corr_coeff(tors_fluct, dt, ma, &slope, &offset);
  fprintf(file_out, "r(therm,change)_Tors  %.3f\n", r);

  // Prepare projections and compute RMSD


  /* Random conf changes with regularization */
  unsigned long iran=randomgenerator();
  if(ini_ran==0){
    InitRandom((RANDOMTYPE)iran); ini_ran=1;
  }
  double Delta_E0_ave=0, m_ave=0;
  float pval=0; int nran=200;
  float c2_thr=c_thr*c_thr/scale2;
  for(int j=0; j<nran; j++){
    double DE=0, sum_c2_ran=0;
    for (i=0; i<nmodes; i++){
      if(sigma2[i]==0)continue;
      float z2=Norm_var_2();
      float c2i=z2*sigma2[i];
      if(REG2){
	float reg=1+Lambda/sigma2[i]; z2/=(reg*reg);
      }else{
	if(c2i<c2_thr)continue;
      }
      m_ave++;
      DE+=z2;
      sum_c2_ran+=c2i;
    } // end modes
    if(sum_c2_ran)DE*=(sum_c2/sum_c2_ran);
    Delta_E0_ave+=DE;
    if(DE<Delta_E)pval++;
  }
  pval/=nran; Delta_E0_ave/=nran; m_ave/=nran;

 // Projection on normal modes and frequency
  /*fprintf(file_out, "Ratio_rel_modes_%.2f  %.3g  %d %d\n",
    sum_thr, nmodes_w2/(float)nmodes_c2, nmodes_w2, nmodes_c2);*/
  //(sum_c2*Delta_E0)/(sum_sigma2*Delta_E));
  fprintf(file_out, "Null_model_rmsd       %.3g\n", sqrt(sum_c2)/M_sqrt);
  fprintf(file_out, "Delta_E_ran/Delta_E   %.3g\n",
	  Delta_E0_ave/Delta_E); 
  fprintf(file_out, "P(Delta_E_ran<Delta_E) %.3g\n",pval);
  //fprintf(file_out, "Selected_modes        %d\n", nmodes); //Delta_E
  fprintf(file_out, "Mean_num. modes       %.1f\n", m_ave);
  Delta_E/=(2.0*M_sqrt*M_sqrt);
  fprintf(file_out, "Delta_E(kBT)          %.3g\n", Delta_E);
  fprintf(file_out, "Delta_E(kBT)/Dr^2     %.3g\n",
	  Delta_E/(Diff->RMSD*Diff->RMSD));

  // Equivalent threshold for normal modes
  //double s2_thr=sum_sigma2*sum_c2/sum_c2;

  // Rank normal modes on frequency
  /*int sort_w[N_modes], sorted[N_modes];
  for(i=0;i<N_modes; i++){sorted[i]=0; sort_w[i]=-1;}
  for(k=0; k<N_modes; k++){
    int ik=-1; float smax=-1; 
    for (i=0;i<N_modes; i++){
      if((sorted[i]==0)&&(sigma2[i]>=smax)&&(NM->select[i])){
	smax=sigma2[i]; ik=i;
      }
    }
    if(ik<0)break;
    sort_w[k]=ik; sorted[ik]=1;
    } */
 
  // Number of modes contributing more than sum_thr to sum_c2 and sum_sigma2
  /* float sum_thr=0.85, sum_thr_w2=sum_thr*sum_sigma2;
  int nmodes_c2=0, nmodes_w2=0;
  double tmp_s2=0, tmp_c2=0;
  for (k=0; k<N_modes; k++){
    i=sort_w[k]; if(i<0)break;
    if(tmp_c2<=sum_thr){tmp_c2+=c2[i]; nmodes_c2++;}
    if(tmp_s2<=sum_thr_w2){tmp_s2+=sigma2[i]; nmodes_w2++;}
  }
  printf("number of modes such that sum< %.2f total: %d (c2), %d (w2) %.3f\n",
	 sum_thr, nmodes_c2, nmodes_w2, nmodes_w2/(float)nmodes_c2);

  int selmodes=0; sum_sigma2=0;
  for (k=0; k<N_modes; k++){
    if(sum_sigma2>s2_thr)break;
    i=sort_w[k]; if(i<0)break;
    if(sigma2[i]>0){
      sum_sigma2+=sigma2[i]; selmodes++;
    }
  }
  printf("selmode= %d sum_sigma2= %.3g sum_c2= %.2g s2_thr= %.3g\n",
	 selmodes, sum_sigma2, sum_c2/sum_d2, s2_thr);
  */


  // Store angles and projections
  float *Dphi=malloc(NM->N_axes*sizeof(float));
  for(int a=0; a<NM->N_axes; a++)Dphi[a]=Diff->Tors[a];
  diff_angles[*num_dphi]=Dphi;
  
  Dphi=malloc(NM->N_axes*sizeof(float));
  for(int a=0; a<NM->N_axes; a++)Dphi[a]=Diff->MW_Tors[a];
  MW_diff_angles[*num_dphi]=Dphi;
  
  float *cc=malloc(NM->N*sizeof(float));
  for (i=0;i<N_modes;i++)cc[i]=c2[i];
  c2_alpha[*num_dphi]=cc;

  (*num_dphi)++;

  //rho1
  float x[N_modes], y[N_modes];  int imax=-1;
  float cmax=-1; k=0; 
  for(i=0; i<nmodes; i++){
    if(NM->select[i]==0)continue;
    //if(NM->Cart_coll[i] < Coll_thr_cc)continue;
    x[k]=sigma2[i]; y[k]=c2[i]; k++;
    if(c2[i]>cmax){cmax=c2[i]; imax=i;}
  }
  r= Corr_coeff(x, y, k, &slope, &offset);
  //*Confchange_mode=imax;

  //if(REG2)
  //fprintf(file_out, "# Warning, regularization violates the null model\n");
  fprintf(file_out, "r[c^2,1/w^2]          %.3f\n", r);
  fprintf(file_out, "slope                 %.2g\n", slope);
  fprintf(file_out, "offset                %.2g\n", offset);
  float Coll=Collectivity_norm1(y, k);
  //Coll=Collectivity_norm1(c2, NM->N_relevant);
  fprintf(file_out, "Recp.Coll(cc)         %.1f\n", Coll);

  /*if(nstruct==0){
    int na=Ref.N_Cart/3;
    fprintf(file_out, "Recp.Coll.(cc)/Natm   %.3f\n", Coll/na);
    fprintf(file_out, "R.Coll.(cc)Renyi/Natm %.3f\n",
      Collectivity_Renyi_norm1(c2,N_modes)/na);
      }*/
  /*fprintf(file_out, "Area(cc)              %.3f\n",
    Area_norm1(c2,N_modes));*/

  if(nstruct==0){
    fprintf(file_out, "Most_contr.mode(cc):  %d %.3f %.3g  %.2f %.2f %.2f\n",
	    imax, cmax, sigma2[imax], NM->Cart_coll[imax],
	    NM->MW_Tors_coll[imax], NM->Tors_coll[imax]);
    fprintf(file_out, "(mode cc 1/w^2 coll_cart coll_Wtors coll_tors)\n");
  }

  // Test of null model
  if(TEST_NULL){
    //float q=Test_Renyi(x, y, k, name1, file_out);
    //~ float q=Test_Null(x, y, k, name1, file_out);
  }

  // rho1_loglog
  /* k=0; cmax=-1;
  for(i=0; i<N_modes; i++){
   if(NM->select[i]){
     x[k]=log(sigma2[i]); y[k]=log(c2[i]); k++;
     if(c2[i]>cmax){cmax=c2[i]; imax=i;}
   }
   }*/

  if(nstruct==0)fprintf(file_out, "#\n");

  // rho
  k=0; cmax=-1;
  //~ float min_pred=-2*mu/slope;
  for(i=0; i<nmodes; i++){
    if(NM->select[i]==0)continue;
    //if(NM->Cart_coll[i] < Coll_thr_cc)continue;
    x[k]=sigma2[i];
    y[k]=(c2[i])/sigma2[i];
    if(y[k]>cmax){cmax=y[k]; imax=i;}  k++;
  }
  r= Corr_coeff(x, y, k, &slope, &offset);
  //~ float k_Therm=Collectivity_norm1(sigma2, N_modes);
  if(nstruct==0)
    fprintf(file_out, "Most_contr.mode(c*w): %d %.3f %.3g  %.2f %.2f %.2f\n",
	    imax, sqrt(cmax), sigma2[imax], NM->Cart_coll[imax],
	    NM->MW_Tors_coll[imax], NM->Tors_coll[imax]);
  fprintf(file_out, "r[(c*w)^2,1/w^2]      %.3f\n", r);
  fprintf(file_out, "Significance(Z)       %.3f nsel=%d\n", r*sqrt(k), k);
  //fprintf(file_out, "Significance(Z)       %.3f\n",r*sqrt(k_Therm));
  Coll=Collectivity_norm1(y, k);
  fprintf(file_out, "Recp.Coll((cc*w)^2)   %.1f\n", Coll);

  // Correlation conformation change and anharmonicity
  if(NM->Anharmonicity[0]){
    float P2[N_modes], norm=norm_c2/sum_c2;
    for(i=0; i<N_modes; i++){
      if(NM->select[i]){P2[i]=c2[i]*norm;}
      else{P2[i]=0;}
    }
    double Anharm_ene=0, Anharm_struct=0;
    for (i=0; i<N_modes; i++){
      if(((NM->Anharmonicity[i]>0)&&(Diff->coeff[i]>0))||
	 ((NM->Anharmonicity[i]<0)&&(Diff->coeff[i]<0)))
	Anharm_ene+=P2[i];
      if(((NM->Anharm_struct[i]>0)&&(Diff->coeff[i]>0))||
	 ((NM->Anharm_struct[i]<0)&&(Diff->coeff[i]<0)))
	Anharm_struct+=P2[i];
    }
    fprintf(file_out,"Frac.conf_change directed as anharmonicity (ene)=%.4f\n",
	    Anharm_ene);
    fprintf(file_out,"Frac.conf_change directed as anharmonicity (str)=%.4f\n",
	    Anharm_struct);
  }

  // Projecting force on normal modes
  double sum=0; k=0; float F_coeff[N_modes];
  for(i=0; i<N_modes; i++){
    if((NM->select[i])&&(sigma2[i])){
      float c=Diff->coeff[i]/sigma2[i], cc=c*c;
      F_coeff[i]=cc;
      sum+=cc;
      //x[k]=log(sigma2[i]); y[k]=log(cc);
      if(cc>cmax){cmax=cc; imax=i;} k++;
    }else{
      F_coeff[i]=0;
    }
  }
  //Corr_coeff(x, y, k, &slope, &offset);
  //fprintf(file_out, "r_ll(therm,force)     %.3f\n", r);
  //fprintf(file_out, "expo(therm,force)     %.2g\n", slope);

  for(i=0; i<NM->N; i++)F_coeff[i]/=sum;
  float k_Force_mode=Collectivity_norm1(F_coeff, N_modes);
  fprintf(file_out, "Recp.Coll.(Force)     %.1f\n", k_Force_mode);
  fprintf(file_out, "Most_contr.mode(force) %d %.3g %.3g  %.2f %.2f %.2f\n",
	  imax, cmax/sum, sigma2[imax], NM->Cart_coll[imax],
	  NM->MW_Tors_coll[imax], NM->Tors_coll[imax]);

  fprintf(file_out, "#\n");

}

void Print_force_confchange(struct Normal_Mode NM, struct Tors Diff,
			    atom *atoms, struct axe *axes, int N_axes,
			    struct residue *seq, int N_res,
			    char *nameout2)
{
  // Squared conformation changes
  struct Tors Diff2;   int i;
  Allocate_tors(&Diff2, NM.N_axes, NM.N_Cart, NM.N);
  for(i=0; i<NM.N_Cart;  i++){
    Diff2.Cart[i]=Diff.Cart[i]*Diff.Cart[i];
  }
  for(i=0; i<NM.N_axes;  i++){
    Diff2.Tors[i]=Diff.Tors[i]*Diff.Tors[i];
    Diff2.MW_Tors[i]=Diff.MW_Tors[i]*Diff.MW_Tors[i];
  }

  Print_change(NM.cart_fluct, Diff2.Cart, N_res,  nameout2, "Cart");
  Print_change(NM.tors_fluct, Diff2.Tors, N_axes, nameout2, "Tors");
  char namepdb[200];
  sprintf(namepdb, "%s1.pdb", nameout2);
  Print_tors_fluct(axes,N_axes,Diff2.Tors,atoms,seq,namepdb,"CONFCHANGE");
  sprintf(namepdb, "%s2.pdb", nameout2);
  Print_tors_fluct(axes,N_axes, NM.tors_fluct,atoms,seq,namepdb,"FLUCT");

  Empty_tors(Diff2);
}

void Print_modes_confchange(char *name2, char *nameout2,
			    struct Normal_Mode NM, struct Tors Diff,
			    int ANM, float M_sqrt, float rmsd,
			    int anharmonic,float xkappa)
{
  /***************** Normal modes confchange ********************/

  char nameout[400];
  sprintf(nameout, "%s.Modes.dat", nameout2);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  fprintf(file_out, "# Conf.change %s\n", name2);
  fprintf(file_out, "# kappa=  %.3f\n",xkappa);

  /*float *sigma2;
  if(anharmonic){
    fprintf(file_out, "# Omega2 corrected for anharmonicity\n");
    sigma2=NM.sigma2_anhar;
  }else{
    sigma2=NM.sigma2;
    }*/


  // Compute confchange and copy on normal modes

  // Sort on c2
  int sorted [NM.N], i, k;
  for(i=0;i<NM.N; i++)sorted[i]=0;
  float *c2=NM.confchange2;
  for(k=0; k<NM.N; k++){
    int ik=-1; float cmax=-1; 
    for (i=0;i<NM.N; i++){
      if((sorted[i]==0)&&(c2[i]>=cmax)&&(NM.select[i])){
	cmax=c2[i]; ik=i;
      }
    }
    if(ik>=0){
      NM.sort[k]=ik; sorted[ik]=1;
    }else{ // nothing found
      ik=k;
      for (i=0;i<NM.N; i++){
	if(NM.select[i]==0){NM.sort[ik]=i; ik++; if(ik==NM.N)break;}
      }
      break;
    }
  }

  // Print lines
  fprintf(file_out, "# Frequency reported in internal units: %.3f ps^(-1),"
	  " kT/h/ is %.3g internal units (om-2)= %.3g\n",
	  Freq_unit, 0.385/Freq_unit, Freq_unit*Freq_unit/0.1482);
  fprintf(file_out,"# mode pc_confch cumul  ");
  fprintf(file_out,"om-2(harm) ");
  if(anharmonic)fprintf(file_out, "om-2(corr) om-2(anha) ");
  fprintf(file_out,"   RMSD  Coll_cart");
  if(NM.MW_Tors_coll)fprintf(file_out, "  Coll_MW");
  if(NM.Tors_coll)fprintf(file_out, "  Coll_tors");

  //fprintf(file_out, "# mode\tContr_to_confchange\tCumulative");
  //if(anharmonic)fprintf(file_out, "\tOmega^(-2)_anharm");
  //fprintf(file_out, "\tOmega^(-2) RMSD ");
  //fprintf(file_out, "\tColl.(cart)\tColl.(MWtors)\tColl.(tors)");
  if(ANM)fprintf(file_out, "\tTors_frac"); //force_coeff
  if(NM.Anharmonicity[0]){
    fprintf(file_out, "\tAnharm_Str Str_confchange");
    fprintf(file_out, "\tAnharm_Ene\tEne_confchange\tMax_RMSD");
  }
  fprintf(file_out, "\n");

  float N_sqrt=sqrt(NM.N_Cart/3);
  double sum_coeff_C=0; //sum_B=0;
  for(k=0; k<NM.N; k++){
    int ik=NM.sort[k]; if(ik<0)continue;
    sum_coeff_C+=c2[ik];
    //sum_B+=NM.sigma2[ik];
    if(NM.select[ik]==0)fprintf(file_out,"#");

    fprintf(file_out, "%6d  %8.3g %5.3f     ",ik, c2[ik], sum_coeff_C);
    fprintf(file_out, "%7.4g  ",NM.sigma2[ik]);
    if(anharmonic)fprintf(file_out,"%7.4g  %7.4g  ",
			  xkappa*NM.sigma2_anhar[ik],NM.sigma2_anhar[ik]);
    fprintf(file_out, "  %8.3g", 1./(N_sqrt*NM.omega[ik]));
    fprintf(file_out, "  %5.3f", NM.Cart_coll[ik]);
    fprintf(file_out, "  %5.3f", NM.MW_Tors_coll[ik]);
    fprintf(file_out, "  %5.3f", NM.Tors_coll[ik]);
	    
    if(ANM)fprintf(file_out,"\t%.3f",NM.Tors_frac[ik]);
    if(NM.Anharmonicity[0]){
      int dir;
      if(((NM.Anharm_struct[ik]>0)&&(Diff.coeff[ik]>0))||
	 ((NM.Anharm_struct[ik]<0)&&(Diff.coeff[ik]<0))){dir=1;}
      else{dir=0;}
      fprintf(file_out,"\t%.3f\t%d",
	      fabs(NM.Anharm_struct[ik]),dir);
      if(((NM.Anharmonicity[ik]>0)&&(Diff.coeff[ik]>0))||
	 ((NM.Anharmonicity[ik]<0)&&(Diff.coeff[ik]<0))){dir=1;}
      else{dir=0;}
      fprintf(file_out,"\t%.3f\t%d\t%.2f",
	      fabs(NM.Anharmonicity[ik]),dir, NM.Max_RMSD[ik]);
    }
    fprintf(file_out, "\n");
  }
  //free(c2);
  fclose(file_out);
}

float Fermi_function(float c2, float c2_thr, float S){
  float x=exp(-S*(fabs(c2)-c2_thr));
  return(1./(1.+x));
}

float Test_Renyi(float *wminus2, float *P, int N, char *name, FILE *file1)
{
  /* Tests the Renyi probabilities P^(q-1)\propto -ene -mu
     for various values of the Reniy parameter q
     Returns the value of q that yields best correlation */
  int Newfile=0, i;
  float *y=malloc(N*sizeof(float));
  float q, q1=0, r, slope, offset, qmax=0, rmax=0;
  FILE *file2; char filename[80];
  if(Newfile){
    sprintf(filename, "Renyi_%s.dat", name);
    file2=fopen(filename, "w");
    fprintf(file2, "# Correlation between w_a^-2 and P_q(a)=c_a^2(q-1)\n");
    fprintf(file2, "# q Corr\n");
  }else{
    fprintf(file1, "# P_q(a)=c_a^2(q-1), d_KL r[P_q,1/w^2]\n");
  }

  float wmin=1000, mu;
  for(i=0; i<N; i++)if(wminus2[i]<wmin)wmin=wminus2[i];
  wmin-=0.0001;
  for(q=1; q<=3.5; q+=0.25){
    double d_KL=0, q2=0, Z=0, lnp;
    if(q==1){for(i=0; i<N; i++)y[i]=log(P[i]);}  // q=1: Shannon
    else{for(i=0; i<N; i++)y[i]=pow(P[i], q1);}
    r= Corr_coeff(wminus2, y, N, &slope, &offset);
    if(q==1){
      for(i=0; i<N; i++){
	d_KL+=P[i]*(y[i]-slope*wminus2[i]);
	Z+=exp(slope*wminus2[i]);
      }
    }else{
      q1=q-1; q2=1/q1; 
      if(offset<0){mu=-wmin;}else{mu=offset/slope;}
      for(i=0; i<N; i++){
	lnp=q2*log(wminus2[i]+mu); Z+=exp(lnp);
	d_KL+=P[i]*(log(P[i])-lnp);
      }
    }
    d_KL+=log(Z);
    if(r > rmax){rmax=r; qmax=q;}
    if(Newfile){fprintf(file2, "%.3f  %.3f  %.3f\n", q, d_KL, r);}
    else{fprintf(file1, "q=%.2f                %.3f  %.3f\n",q,d_KL,r);}
  }
  fprintf(file1, "Renyi parameter qopt  %.3f\n", qmax);
  fprintf(file1, "r[P_qopt,1/w^2]       %.3f\n#\n", rmax);
  if(Newfile){
    fclose(file2);
    printf("Writing %s\n", filename);
  }
  free(y);
  return(qmax);
}

float Test_Null(float *wminus2, float *P, int N, char *name, FILE *file1)
{
  /* Tests the null model P0(q) \propto w^-q varying q
     Returns the value of q that yields best correlation */

  int Newfile=0, i;
  float *x=malloc(N*sizeof(float));
  FILE *file2; char filename[80];
  if(Newfile){
    sprintf(filename, "Null_model_%s.dat", name);
    file2=fopen(filename, "w");
    fprintf(file2, "# Correlation between c_a^2 and P0_q(a)=w_a^(-q)\n");
    fprintf(file2, "# q Corr\n");
  }else{
    fprintf(file1, "# P_q(a)=w^(-q) d_KL(w^-q,c_a^2) d2 r(w^q*c_a^2,w^-2)\n");
  }

  float q, rho, slope, offset;
  float dmin=1000, q_dmin=0, rho_opt=0;
  //float q_rho=0;
  double P1=0, P2=0, dP=0;
  for(i=0; i<N; i++){P1+=P[i]; P2+=P[i]*P[i]; dP+=P[i]*log(P[i]);}
  for(i=0; i<N; i++)P[i]/=P1;
  P2=P2/(P1*P1)-1./N; dP=dP/P1-log(P1);
  //printf("P2= %.3f -S[P]=%.3f N=%d\n", P2, dP, N);
  for(q=0.5; q<=8; q+=0.25){
    double  q1=q/2., x1=0, d_KL=0, d2=0, d;
    for(i=0; i<N; i++){x[i]=pow(wminus2[i], q1); x1+=x[i];}
    for(i=0; i<N; i++)d_KL+=P[i]*log(x[i]);
    d_KL=dP-(d_KL-log(x1));
    if(d_KL < dmin){dmin=d_KL; q_dmin=q;}
    for(i=0; i<N; i++){d=P[i]-x[i]/x1; d2+=d*d;} d2/=P2;
    for(i=0; i<N; i++)x[i]=P[i]/x[i];
    rho= Corr_coeff(x, wminus2, N, &slope, &offset);
    if(fabs(rho)<rho_opt)rho_opt=fabs(rho); //q_rho=q;
    if(Newfile){
      fprintf(file2, "%.3f\t%.3f\t%.3f\t%.4f\n", q, d_KL, d2, rho);
    }else{
      fprintf(file1, "q=%.2f\t%.3f\t%.3f\t%.4f\n",q,d_KL,d2,rho);
    }
  }
  fprintf(file1, "Null model qopt       %.3f\n", q_dmin);
  fprintf(file1, "r^2[P_qopt,P_obs]     %.3f\n#\n", 1-dmin);
  if(Newfile){
    fclose(file2);
    printf("Writing %s\n", filename);
  }
  free(x);
  return(q_dmin);
}

int Select_cc2(float *dx2, int nres1, float *Dcart, atom *atoms1,
	       struct Reference Ref, char *name)
{
  int n=0, i, j;
  for(i=0; i<Ref.N_ref; i++){
    if(strncmp(atoms1[Ref.atom_num[i]].name, name, 2)!=0)continue;
    float *Dc=Dcart+3*i; double d2=0;
    for(j=0; j<3; j++){d2+=(*Dc)*(*Dc); Dc++;}
    dx2[n]=d2; n++;
    if(n > nres1){
      printf("ERROR in Select_cc2, too many %s atoms > %d residues\n",
	     name, nres1); return(0);
    }
  }
  if(n!=nres1){
    printf("WARNING in Select_cc2, %d %s atoms but %d residues\n",
	   n, name, nres1);
  }
  return(n);
}

int Periodic_angles(float *dphi, struct axe *axe, int n)
{
  float pi=acos(-1), twopi=2*pi; int m=0;
  for(int i=0; i<n; i++){
    if(axe[i].type=='l')continue;
    if((dphi[i]<=pi)&&(dphi[i]>=-pi))continue;
    while(dphi[i]>pi)dphi[i]-=twopi;
    while(dphi[i]<-pi)dphi[i]+=twopi;
    m++;
  }
  printf("%d angles > pi=%.5f changed\n", m, pi);
  return(m);
}

float RMS(float *v, int n, int *outlier){
  double RMS=0; int m=n;
  for(int i=0; i<n; i++){
    if(outlier && (outlier[i])){m--; continue;}
    RMS+=v[i]*v[i];
  }
  RMS=sqrt(RMS/m);
  return(RMS);
}

/////////////////////////////////////////////////////////////////////
void Print_angular_differences(float **diff_angles, char *type,
			       char **name_dphi, int num_dphi,
			       char *name_para, struct axe *axe,
			       int N_axes, struct residue *seq1)
{
  char nameout[100]; int a, i;
  sprintf(nameout, "%s.Angular_diff%s.dat", name_para, type);
  printf("Writing angular differences in %s\n", nameout);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "# %d angles %d methods for deriving angles\n",
	  N_axes, num_dphi);
  for(i=0; i<num_dphi; i++)
    fprintf(file_out, "# Method %d: %s\n", i, name_dphi[i]);	     
  fprintf(file_out, "#angle ");
  for(i=0; i<num_dphi; i++)fprintf(file_out, " %d=dphi%d", i+2, i+1);
  fprintf(file_out, "\n");
  for(a=0; a<N_axes; a++){
    struct residue *r=seq1+(axe[a].bond)->atom->res;
    fprintf(file_out, "%c%s%c_%c", axe[a].type,r->pdbres,r->chain,r->amm);
    for(i=0; i<num_dphi; i++){
      fprintf(file_out, " %.4f", diff_angles[i][a]);
    }
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}

int Mode_confchange_fit(struct Tors *Diff, struct Normal_Mode *NM,
			struct bond *bonds, atom *atoms, int natoms,
			struct Reference Ref, struct residue *seq,
			float *coord_1, float *coord_2, float M_sqrt,
			char *nameout, char *name_pdb, int type, int NMODES)
{
  int nround=1, nreduce=0, itmax, i, j, k, num_mod=0;
  float STEP_MAX=0.15, STEP_MIN=0.001, TOL=0.1, TOL2=0.000;
  float c_max=STEP_MAX*M_sqrt, c_thr=STEP_MIN*M_sqrt;

  idum=randomgenerator();
  InitRandom((RANDOMTYPE)idum);
  float gasdev(long *idum);
  float rmsd_ini=rmsd_mclachlan_f(coord_1,coord_2,Ref.mass_atom,Ref.N_ref);

  char out[100]; sprintf(out, "%s.mode_confchange.dat", nameout);
  printf("Writing %s\n", out);
  FILE *file_out=fopen(out, "w");
  sprintf(out, "%s.mode_confchange.pdb", nameout);
  printf("Writing %s\n", out);
  FILE *file_pdb=fopen(out, "w");
  char what[100];
  sprintf(what, "Torsional conformations, norm. modes ranked by"); 
  if(type){strcat(what, " overlap");}
  else{strcat(what, " frequency");}
  fprintf(file_out, "# %s\n", what);
  fprintf(file_out,"#1=rmsd_min 2=rmsd_tar 3=mode 4=c_alpha\n");
  fprintf(file_out, "%.3g\t%.3g\t0\t0\n",rmsd_ini,rmsd_ini);

  int naxes=NM->N_axes, ncoord=NM->N_Cart;
  float d_phi_opt[naxes], d_phi[naxes], d_phi_best[naxes];
  float coord_new[ncoord], coord_step[ncoord], coord_all[3*natoms];
  float D_coord[ncoord];   
  float sigma[NM->N]; for(int i=0; i<NM->N; i++)sigma[i]=sqrt(NM->sigma2[i]);
  for(i=0; i<naxes; i++)d_phi_opt[i]=0;
  for(i=0; i<ncoord; i++){
    coord_step[i]=coord_1[i];
    D_coord[i]=coord_2[i]-coord_step[i];
  }
  
  int jmode[NM->N]; float c_lin[NM->N];
  if(type==0){
    for(k=0; k<NMODES; k++){
      jmode[k]=k;
      c_lin[k]=Scalar_product_weighted(D_coord, NM->Cart[k],
				       Ref.mass_coord, ncoord);
    }
  }else{
    for (i=0; i<NM->N; i++){
      c_lin[i]=Scalar_product_weighted(D_coord, NM->Cart[i],
				       Ref.mass_coord, ncoord);
    }
    // Rank the modes by the correlation with conf.change
    int done[NM->N]; for(k=0; k<NM->N; k++)done[k]=0;
    for(i=0; i<NMODES; i++)jmode[i]=-1;
    for(i=0; i<NMODES; i++){
      int kj=-1; float cmax=0;
      for(k=0; k<NM->N; k++){
	if(done[k])continue;
	if(kj<0 || fabs(c_lin[k])>fabs(cmax)){kj=k; cmax=c_lin[k];}
      }
      if(cmax==0 || kj<0)break;
      done[kj]=1;
      jmode[i]=kj;
    }
  }

  int NMOD1=NMODES-1, accept=1, random=0;
  float score_opt=-1000, score_old=-1000, score_best=-1000,
    tol, rfact=0.5, yy, cc=0;

  for(int round=0; round<nround; round++){

    int done[NMODES], kj;
    if(random){tol=TOL; itmax=NMODES*0.5;}
    else{tol=0; for(i=0; i<NMODES; i++)done[i]=0; itmax=10*NMODES;}

    for(int iter=0; iter<=itmax; iter++){

      // Full reconstruction, optimize the coefficient
      // Determine mode (kj) and coefficient (cc)
      if(random){
	int j=NMODES*RandomFloating(); if(j>NMOD1)j=NMOD1;
	kj=jmode[j]; if(kj<0)continue;
	cc=gasdev(&idum)*sigma[kj]*rfact;
      }else{ // Minimization
	kj=-1; cc=0;
	for(j=0; j<NMODES; j++){
	  k=jmode[j]; if(k<0 || done[k])continue;
	  if(kj<0 || fabs(c_lin[k])>fabs(cc)){kj=k; cc=c_lin[k];}
	}
	if(kj<0 || fabs(cc)<c_thr){
	  fprintf(file_out,"# Too small max. rmsd: %.2g mode=%d\n",
		  STEP_MIN, kj);
	  break;
	}
      }

      // Compute RMSD
      if(fabs(cc)>c_max){if(cc>0){cc=c_max;}else{cc=-c_max;}}
      for(int j=0; j<=nreduce; j++){
	yy=-RMSD_confchange(cc, NM->Tors[kj], d_phi, coord_all,
			    coord_new, d_phi_opt, naxes, bonds,atoms,
			    natoms,Ref,coord_2);
	fprintf(file_out, "%.3g\t%.3g\t%d\t%.3g\n",-score_opt,-yy,kj,cc);

	// Test if improvement
	if(yy>score_opt-tol){
	  accept=1; score_old=score_opt; score_opt=yy; 
	  for(i=0; i<naxes; i++)d_phi_opt[i]=d_phi[i];
	  break;
	}else{accept=0; cc/=2; done[kj]=1;}

      }

      if(accept==0)continue;
      
      // Update structure
      for(i=0; i<ncoord; i++){
	coord_step[i]=coord_new[i];
	D_coord[i]=coord_2[i]-coord_new[i];
      }
      for(i=0; i<NMODES; i++){
	k=jmode[i];
	c_lin[k]=Scalar_product_weighted(D_coord, NM->Cart[k],
					 Ref.mass_coord, ncoord);
      } 
      if(yy>score_best){
	score_best=yy; 
	for(i=0; i<naxes; i++)d_phi_best[i]=d_phi[i];
      }

      if(fabs(score_old-score_opt)<TOL2 && random==0){
	fprintf(file_out, "# RMSD does not improve, %.4g -> %.4g\n",
		score_old, score_opt);
	done[kj]=1;
      }
    } // end iter
    if(random){random=0; rfact/=2;}else{random=1;}
  } // end round


  yy=-RMSD_confchange(0, NM->Tors[0], d_phi, coord_all, coord_new,
		      d_phi_best, naxes, bonds,atoms, natoms,Ref,coord_2);
  fprintf(file_out, "&\n");
  fprintf(file_out, "# %d modes min. RMSD=%.3f\n", NMODES, -yy);
    
  fclose(file_out);

  num_mod++;
  fprintf(file_pdb, "MODEL %d\n", num_mod);
  fprintf(file_pdb, "REMARK %s\n", what);
  fprintf(file_pdb, "REMARK RMSD from target: %.2f\n", -score_opt);
  fprintf(file_pdb, "REMARK  The fitted conformation change is ");
  fprintf(file_pdb, "obtained with %d normal modes\n", NMODES);
  float coord_all_1[3*natoms], *c=coord_all_1, mass[natoms];
  for(i=0; i<natoms; i++){
    mass[i]=1;
    double *r=atoms[i].r; for(k=0; k<3; k++){*c=*r; c++; r++;}
  }
  float rmsd=rmsd_mclachlan_f(coord_all_1, coord_all, mass, natoms);
  rmsd=rmsd_mclachlan_f(coord_1,coord_new,Ref.mass_atom,Ref.N_ref);
  Print_PDB(file_pdb, atoms, natoms, coord_all, seq, num_mod, rmsd);
  fclose(file_pdb);

  return(0);
}

int Mode_confchange_old(struct Tors *Diff, struct Normal_Mode *NM,
			struct bond *bonds, atom *atoms, int natoms,
			struct Reference Ref, struct residue *seq,
			float *coord_1, float *coord_2, float M_sqrt,
			char *nameout, char *name_pdb)
{
  int NMODES_CC=10, nround=1, itmax=10, nreduce=4, i, j, k, num_mod=0;
  float STEP_MAX=0.4, c_thr=0.02*M_sqrt, TOL=0.04, TOL2=0.01;
  float c_max=STEP_MAX*M_sqrt;
  char out[100]; sprintf(out, "%s.mode_confchange.dat", nameout);
  printf("Writing %s\n", out);
  FILE *file_out=fopen(out, "w");
  sprintf(out, "%s.mode_confchange.pdb", nameout);
  printf("Writing %s\n", out);
  FILE *file_pdb=fopen(out, "w");

  float rmsd_ini=rmsd_mclachlan_f(coord_1,coord_2,Ref.mass_atom,Ref.N_ref);
  int naxes=NM->N_axes, ncoord=NM->N_Cart;
  float d_phi_opt[naxes], d_phi[naxes];
  float coord_new[ncoord], coord_step[ncoord], coord_all[3*natoms];
  float D_coord[ncoord];   int done[NM->N];
  float gasdev(long *idum);
  float sigma[NM->N]; for(int i=0; i<NM->N; i++)sigma[i]=sqrt(NM->sigma2[i]);
  idum=randomgenerator();
  InitRandom((RANDOMTYPE)idum);

  for(int type=0; type<1; type++){
    for(int lin=0; lin<2; lin++){
      if(lin==1 && type==1)continue;

      char what[100];
      if(lin){sprintf(what, "Linear ");}
      else{sprintf(what, "Torsional ");}
      strcat(what, " conformations, norm. modes ranked by"); 
      if(type){strcat(what, " overlap");}
      else{strcat(what, " frequency");}
      fprintf(file_out, "# %s\n", what);
      fprintf(file_out,"#1=n.modes 2=rmsd_min 3=rmsd 4=mode 5=c_alpha\n");
      fprintf(file_out, "0\t%.3g\t%.3g\t0\t0\n",rmsd_ini,rmsd_ini);

      if(type)for(k=0; k<NM->N; k++)done[k]=0;
      for(i=0; i<naxes; i++)d_phi_opt[i]=0;
      for(i=0; i<ncoord; i++){
	coord_step[i]=coord_1[i];
	D_coord[i]=coord_2[i]-coord_step[i];
      }
      float c_lin[NM->N];
      for (i=0; i<NM->N; i++){
	c_lin[i]=Scalar_product_weighted(D_coord, NM->Cart[i],
					 Ref.mass_coord, ncoord);
      }

      int jmode[NM->N], imode, accept=1, kj=0;
      float score_opt=-1000, score_old=-1000, yy, cc=0;

      for(imode=0; imode<NM->N; imode++){

	if(type==0){
	  kj=imode; // Rank the modes by the frequency
	}else{ 
	  // Rank the modes by the correlation with conf.change
	  kj=-1; float cmax=0;
	  for(k=0; k<NM->N; k++){
	    if(done[k])continue;
	    if(kj<0 || fabs(c_lin[k])>fabs(cmax)){kj=k; cmax=c_lin[k];}
	  }
	  if(cmax==0 || kj<0)break;
	  done[kj]=1;
	}
	jmode[imode]=kj;

	int random=0; float tol, rfact=0.5;
	for(int round=0; round<nround; round++){

	  int done2[NM->N];
	  if(random){tol=TOL;}
	  else{tol=0; for(i=0; i<NM->N; i++)done2[i]=0;}

	  for(int iter=0; iter<=itmax; iter++){

	    if(lin){
	      cc=c_lin[kj]; float *cart=NM->Cart[kj];
	      for(i=0; i<ncoord; i++)coord_new[i]=coord_step[i]+cc*cart[i];
	      yy=-rmsd_mclachlan_f(coord_2,coord_new,Ref.mass_atom,Ref.N_ref);
	      if(yy>=score_opt){
		accept=1; score_opt=yy;
		float *mode=NM->Tors[kj];
		for(i=0; i<naxes; i++)d_phi_opt[i]+=cc*mode[i];
	      }else{accept=0;}
	      goto print_rmsd;
	    }

	    // Full reconstruction, optimize the coefficient
	    // Determine mode (kj) and coefficient (cc)
	    if(random){
	      if(iter>imode)break;
	      int iran=imode*RandomFloating();
	      if(iran>imode)iran=imode;
	      kj=jmode[iran];
	      cc=gasdev(&idum)*sigma[kj]*rfact;
	    }else{ // Minimization
	      kj=-1; cc=0;
	      for(j=imode; j>=0; j--){
		k=jmode[j]; if(done2[k])continue;
		if(kj<0 || fabs(c_lin[k])>fabs(cc)){kj=k; cc=c_lin[k];}
	      }
	      if(kj<0 || fabs(cc)<c_thr)break;
	    }

	    // Compute RMSD
	    if(fabs(cc)>c_max){if(cc>0){cc=c_max;}else{cc=-c_max;}}
	    for(int j=0; j<nreduce; j++){
	      yy=-RMSD_confchange(cc, NM->Tors[kj], d_phi, coord_all,
				  coord_new, d_phi_opt, naxes, bonds,atoms,
				  natoms,Ref,coord_2);
	      // Test if improvement
	      if(yy>score_opt-tol){
		accept=1; score_old=score_opt; score_opt=yy; 
		for(i=0; i<naxes; i++)d_phi_opt[i]=d_phi[i];
		break;
	      }else{accept=0; cc/=2; done2[kj]=1;}
	    }
	  print_rmsd:
	    fprintf(file_out, "%d\t%.3g\t%.3g\t%d\t%.3g\n",
		    imode,-score_opt,-yy,kj,cc);
	    if(accept==0)continue;
	    
	    for(i=0; i<ncoord; i++){
	      coord_step[i]=coord_new[i];
	      D_coord[i]=coord_2[i]-coord_new[i];
	    }
	    for (i=0; i<NM->N; i++){
	      c_lin[i]=Scalar_product_weighted(D_coord, NM->Cart[i],
					       Ref.mass_coord, ncoord);
	    } 

	    if(lin)goto print_pdb;
	    if(fabs(score_old-score_opt)<TOL2)break;
	  } // end iter
	  if(random){random=0; rfact/=2;}else{random=1;}
	} // end round
	
      print_pdb:
	if(imode==NMODES_CC){ //&& sel
	  num_mod++;
	  fprintf(file_pdb, "MODEL %d\n", num_mod);
	  fprintf(file_pdb, "REMARK %s\n", what);
	  fprintf(file_pdb, "REMARK RMSD from target: %.2f\n", -score_opt);
	  fprintf(file_pdb, "REMARK  The fitted conformation change is ");
	  fprintf(file_pdb, "obtained with %d normal modes\n", NMODES_CC);
	  float coord_all_1[3*natoms], *c=coord_all_1, mass[natoms];
	  for(i=0; i<natoms; i++){
	    mass[i]=1;
	    double *r=atoms[i].r; for(k=0; k<3; k++){*c=*r; c++; r++;}
	  }
	  float rmsd=rmsd_mclachlan_f(coord_all_1, coord_all, mass, natoms);
	  rmsd=rmsd_mclachlan_f(coord_1,coord_new,Ref.mass_atom,Ref.N_ref);
	  Print_PDB(file_pdb, atoms, natoms, coord_all, seq, num_mod, rmsd);
	}
    
      } // end imode
  
      if(type==1 &&  lin==0) // Store angles
	for(i=0; i<naxes; i++)Diff->Tors[i]=d_phi_opt[i];
      
      yy=-RMSD_confchange(0, NM->Tors[0], d_phi, coord_all, coord_new,
			  d_phi_opt,naxes,bonds,atoms,natoms,Ref,coord_2);
      fprintf(file_out, "#%d modes RMSD= %.3f RMSD tors. reconstr.= %.3f\n",
	      imode, -score_opt, -yy);
      fprintf(file_out, "&\n");
      printf("%d best modes RMSD=%.3f\n", imode, -score_opt);
      
      /*fprintf(file_out,"# %d non-zero coefficients of normal modes:\n",imode);
	for(i=0; i<NM->N; i++)if(cnew[i])
	fprintf(file_out,"%d %.3g\n",i,cnew[i]);*/
    } // end lin
  }// end type
  fclose(file_out);
  fclose(file_pdb);
  return(0);
}

int Lasso_confchange_fit(struct Tors *Diff, struct Normal_Mode *NM,
			 struct bond *bonds, atom *atoms, int natoms,
			 struct Reference Ref, float *coord_1, float *coord_2,
			 char *nameout)
{
  char out[100]; sprintf(out, "%s_Lasso.dat", nameout);
  printf("Writing %s\n", out);
  FILE *file_out=fopen(out, "w");
  
  float c[NM->N]; double c2_sum=0, sum_c1=0; int i;
  for (i=0; i<NM->N; i++){
    c[i]=Scalar_product_weighted(Diff->Cart, NM->Cart[i],
				 Ref.mass_coord, Ref.N_Cart);
    c2_sum+=c[i]*c[i]; sum_c1+=fabs(c[i]);
  }
  sum_c1/=NM->N; // average projection
  fprintf(file_out, "# %d modes <|c|>=%.4g\n#lambda rmsd(conf,str2)\n",
	  NM->N, sum_c1);

  float cnew[NM->N], coord_new[3*Ref.N_ref], coord_all[3*natoms];
  float d_phi[NM->N_axes];

  float x[3], y[3], lambda=sum_c1, lambda_step=lambda*0.2;
  y[1]=-RMSD_Lasso_confchange(0, cnew, d_phi, coord_all, coord_new,
			      c,c2_sum,NM,bonds,atoms,natoms,Ref,coord_2);
  fprintf(file_out, "%.3f %.3f\n", 0.0, -y[1]);
  x[1]=lambda;
  y[1]=-RMSD_Lasso_confchange(lambda, cnew, d_phi, coord_all, coord_new,
			      c,c2_sum,NM,bonds,atoms,natoms,Ref,coord_2);
  fprintf(file_out, "%.3f %.3f\n", x[1], -y[1]);
  int n_down=0; float score_opt=y[1], lambda_opt=lambda;
  
  for(int iter=0; iter<100; iter++){
    if(n_down==1)break;
    x[0]=x[1]-lambda_step; if(x[0]<0)x[0]=0;
    x[2]=x[1]+lambda_step;
    for(i=0; i<=2; i+=2){
      y[i]=-RMSD_Lasso_confchange(x[i], cnew, d_phi, coord_all, coord_new,
				  c,c2_sum,NM,bonds,atoms,natoms,Ref,coord_2);
      fprintf(file_out, "%.3f %.3f\n", x[i], -y[i]);
    }
    lambda=Find_max_quad(x[0],x[1],x[2],y[0],y[1],y[2],0,100);
    float yy=-RMSD_Lasso_confchange(lambda, cnew, d_phi, coord_all, coord_new,
				    c,c2_sum,NM,bonds,atoms,natoms,Ref,coord_2);
    x[1]=lambda; y[1]=yy;
    fprintf(file_out, "%.3f %.3f\n", x[1], -y[1]);
    if(yy>score_opt){
      lambda_opt=lambda;
      if(fabs(score_opt-yy)<0.01)break;
      score_opt=yy; n_down=0;
    }else{
      n_down++;
    }
  }

  RMSD_Lasso_confchange(lambda_opt, cnew, d_phi, coord_all, coord_new,
			c,c2_sum,NM,bonds,atoms,natoms,Ref,coord_2);
  for(i=0; i<NM->N_axes; i++)Diff->Tors[i]=d_phi[i];

  fclose(file_out);
  return(0);
}

float RMSD_Lasso_confchange(float lambda, float *cnew, float *d_phi,
			    float *coord_all, float *coord_new,
			    float *c, double c2_sum,
			    struct Normal_Mode *NM, struct bond *bonds,
			    atom *atoms, int natoms,
			    struct Reference Ref, float *coord_2)
{
  int naxes=NM->N_axes, k, i; double norm=0;
  for(k=0; k<NM->N; k++){
    if(fabs(c[k])<lambda){cnew[k]=0; continue;}
    else{cnew[k]=c[k]*(1-lambda/fabs(c[k]));}
    norm+=cnew[k]*cnew[k];
  }
  norm=c2_sum/norm; for(k=0; k<NM->N; k++)if(cnew[k])cnew[k]/=norm;
  // angular distance
  for(i=0; i<naxes; i++)d_phi[i]=0;
  for(k=0; k<NM->N; k++){
    if(cnew[k]==0)continue;
    float *mode=NM->Tors[k];
    for(i=0; i<naxes; i++)d_phi[i]+=cnew[k]*mode[i];
  }
  Set_bonds_measure(bonds, natoms, atoms);
  Build_up(bonds, natoms, d_phi, naxes);
  Put_coord(coord_all, bonds, natoms);
  Write_ref_coord(coord_new, Ref.N_ref, coord_all, Ref.atom_num);
  float rmsd2=rmsd_mclachlan_f(coord_2, coord_new, Ref.mass_atom, Ref.N_ref);
  return(rmsd2);
}

float RMSD_confchange(float cc, float *mode, float *d_phi,
		      float *coord_all, float *coord_new,
		      float *d_phi_opt, int naxes, struct bond *bonds,
		      atom *atoms, int natoms, struct Reference Ref,
		      float *coord_2)
{
  // angular distance
  for(int i=0; i<naxes; i++)d_phi[i]=d_phi_opt[i]+cc*mode[i];
  Set_bonds_measure(bonds, natoms, atoms);
  Build_up(bonds, natoms, d_phi, naxes);
  Put_coord(coord_all, bonds, natoms);
  Write_ref_coord(coord_new, Ref.N_ref, coord_all, Ref.atom_num);
  float rmsd2=rmsd_mclachlan_f(coord_2, coord_new, Ref.mass_atom, Ref.N_ref);
  return(rmsd2);
}

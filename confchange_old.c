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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TEST_NULL 0  // Test best parameter of null model

//#define DIFF_THR 0.25
#define COLL_THR 2   // Select if exp(-S_cart) > COLL_THR

/*extern float Torsional_superimposition(double *x1, double *x2,
				       int N_axes, struct Reference Ref,
				       struct Jacobian J);
extern Test_Diff(struct Tors *Diff, struct Reference Ref, struct Jacobian J);
*/

void Print_diff(int *num_dphi, char **name_dphi,
		float **diff_angles, float **MW_diff_angles,
		float **c2_alpha, int *outlier_tors,
		struct Tors *Diff,
		struct Jacobian *J,
		struct Normal_Mode *NM, float M_sqrt, float thr, 
		struct bond *bonds, atom *atoms1, int natoms1,
		struct Reference Ref, float *coord_old, float *coord_target,
		struct residue *seq,
		struct axe *axe,
		int nstruct, char *name,
		FILE *file_out, FILE *file_pdb);
void Print_angular_differences(float **diff_angles, char *type,
			       char **name_dphi, int num_dphi,
			       char *name_para, struct axe *axe,
			       int N_axes, struct residue *seq1);

float Fermi_function(float c2, float c2_thr, float S);
float Test_Renyi(float *ene, float *P, int N, char *name, FILE *file);
float Test_Null(float *omega_minus2, float *P, int N, char *name, FILE *file);
int Select_cc2(float *dx2, int nres1, float *Dcart, atom *atoms1,
	       struct Reference Ref, char *name);
int Periodic_angles(float *dphi, struct axe *axe, int n);
static float RMS(float *v, int n, int *outlier);

/************************************************************************/

float Examine_confchange(struct Tors *Diff,
			 struct bond *bonds, struct axe *axe,
			 char *nameout, char *name1, char *namepdb,
			 atom *atoms1, float *coord_ref1,
			 struct residue *seq1, int nres1,
			 int N_diso1, int natoms1,
			 atom *atoms2, float *coord_ref2,
			 struct residue *seq2, int nres2,
			 int N_diso2, int natoms2,
			 struct interaction *Int_list1, int N_int1,
			 char *INT_TYPE, float s0, int N_modes,
			 struct Reference Ref, struct Jacobian *J,
			 struct Normal_Mode NM, int nprint,
			 struct Ali_score ali, 
			 struct Para_confchange Para_confchange,
			 float *diff_phi, float *greedy_phi,
			 int nstruct, float *B_TNM,
			 int *outlier_tors,
			 char *name_para)
// nstruct = 0 if experimental confchange, = k if simulated structures
{
  int i;
  // Only modes that produce a conformation change > DIFF_THR are used
  printf("Computing confchange, %d modes %d relevant\n",
	 N_modes, NM.N_relevant);

  // Open file for printing, one chain
  FILE *file_out;
  if(nprint){
    file_out=fopen(nameout, "a");
  }else{
    file_out=fopen(nameout, "w");
    printf("Writing %s\n", nameout);
    fprintf(file_out, "protein               %s\n", name1);
    fprintf(file_out, "Nres                  %d\n", nres1);
    fprintf(file_out, "Degrees of freedom    %d\n", N_modes);
    fprintf(file_out, "Disordered gaps       %d\n", N_diso1);
    fprintf(file_out, "#\n");
  }

  // Printing sequence comparison
  if(nstruct==0){
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
  Compute_interactions(&N_int2, &Int_list2, INT_TYPE, atoms2, natoms2, nres2);
  float Cont_overlap=
    Contact_overlap(Int_list1, N_int1, Int_list2, N_int2, ali.alignres);
  float E1=Contact_energy(Int_list1, N_int1, seq1);
  float E2=Contact_energy(Int_list2, N_int2, seq2);
  float DE=E2-E1-(nres1-nres2)*s0;
  free(Int_list2);

  // Printing properties of conformation change
  if(nstruct==0){
    fprintf(file_out, "Disordered axis2      %d\n", N_diso2);
  }else{
    fprintf(file_out, "Structure             %d\n", nstruct);
  }
  fprintf(file_out, "Contact_overlap       %.3f\n", Cont_overlap);
  fprintf(file_out, "(ECont2-ECont1)       %.3f\n", DE);

  // Superimpose str2 on str1 based on rmsd
  float rmsd=rmsd_mclachlan_f(coord_ref1, coord_ref2, Ref.mass_atom, Ref.N_ref);
  for (i=0;i<Ref.N_cart;i++)Diff->Cart[i]=coord_ref2[i]-coord_ref1[i];
  // float rmsd=Torsional_superimposition(atom_str1,atom_str2,NM.N_axes,Ref,J);

  fprintf(file_out, "RMSD based on atoms   %s\n", REF_CC);
  fprintf(file_out, "RMSD                  %.2f\n", rmsd);

  if(B_TNM){
    // Compute correlation between Cartesian displacements and conf. change
    float dx2[nres1];
    int n=Select_cc2(dx2, nres1, Diff->Cart, atoms1, Ref, "CA");
    if(n){
      float slope, off, r=Corr_coeff(B_TNM, dx2, n, &slope, &off);
      fprintf(file_out, "r(therm,change)_Cart  %.3f\n", r);
    }
  }
  fprintf(file_out, "#\n");

  /************************* Torsional fits ****************************/
  // Projection on normal modes
  float thr=Para_confchange.RMSD_THR, Lambda=0;
  int num_dphi=0, max_dphi=5;
  float *diff_angles[max_dphi], *MW_diff_angles[max_dphi], *c2_alpha[max_dphi];
  char *name_dphi[max_dphi];

  double M_tot=0, M_sqrt;
  for(i=0; i<Ref.N_cart; i+=3)M_tot+=Ref.mass_coord[i];
  M_sqrt=sqrt(M_tot);

  /* Exclude modes contributing little to thermal motion
    Not a good idea!
  float mean=1./(float)N_modes;
  for(i=0; i<N_modes; i++){
    if(NM.contr2fluct[i]<mean)selected[i]=0;
    }*/

  FILE *file_pdb=NULL;
  if(nstruct==0){
    char name[100]; sprintf(name, "%s_confchange.pdb", namepdb);
    file_pdb=fopen(name, "w");
  }

  // Compute the differences between the two superimposed structures
  // and fit them with torsion angles, Lambda=0
  //Convert_cart2torsion(Diff, Ref, J); // There is some mistake in this...
  Convert_cart2torsion_fit(Diff, Ref, J, name1, 'O', &Lambda);
  Print_diff(&num_dphi, name_dphi, diff_angles, MW_diff_angles, c2_alpha,
	     outlier_tors, Diff, J, &NM, M_sqrt, thr,
	     bonds, atoms1, natoms1, Ref, coord_ref1, coord_ref2, seq1, axe,
	     nstruct, "Fit Lambda=0", file_out, file_pdb);

  // Fits with ridge regression M and C
  if(nstruct==0){ 
    for(i=0; i<2; i++){
      char type='M'; if(i)type='C';
      char name[40]; sprintf(name, "Fit RRR type %c", type);
      if(Convert_cart2torsion_fit(Diff, Ref, J, name1, type, &Lambda)<0){
	fprintf(file_out, "WARNING, ridge regression has failed\n");
	printf("WARNING, ridge regression has failed\n");
	continue;
      }
      Print_diff(&num_dphi, name_dphi, diff_angles, MW_diff_angles, c2_alpha,
		 outlier_tors, Diff, J, &NM, M_sqrt, thr,
		 bonds, atoms1, natoms1, Ref, coord_ref1, coord_ref2,
		 seq1, axe, nstruct, name, file_out, file_pdb);
    }
  }

  // Compute the differences of torsion angles
  if((diff_phi)||(greedy_phi)){
    struct Tors Diff_phi; float *diff;
    Allocate_tors(&Diff_phi, NM.N_axes, Ref.N_cart, NM.N);
    for(i=0; i<Ref.N_cart; i++)Diff_phi.Cart[i]=Diff->Cart[i];

    for(int k=0; k<=1; k++){
      char name[80]; sprintf(name, "Differences of torsion angles");
      if(k==0){
	diff=diff_phi; sprintf(name, "%s, different bonds", name);
      }else{
	diff=greedy_phi; sprintf(name, "%s, iterative", name);
      }
      if(diff==NULL)continue;
      for(i=0; i<NM.N_axes; i++)Diff_phi.Tors[i]=diff[i];
      //Compute_MW(&Diff_phi, J);
      Print_diff(&num_dphi, name_dphi, diff_angles, MW_diff_angles, c2_alpha,
		 outlier_tors, &Diff_phi, J, &NM, M_sqrt, thr,
		 bonds, atoms1, natoms1, Ref, coord_ref1, coord_ref2,
		 seq1, axe, nstruct, name, file_out, file_pdb);
      float a0, a1,
	r=Corr_coeff(Diff->Tors, diff, NM.N_axes, &a0, &a1);
      fprintf(file_out, "Corr.coeff(dphi_RRR, dphi)= %.3f\n", r);
    }
    Empty_tors(Diff_phi);
  }
  fclose(file_out);


  // Print list of angles and correlations with normal modes
  if(name_para){
    Print_angular_differences(diff_angles, "", name_dphi, num_dphi,
			      name_para, axe, NM.N_axes, seq1);
    Print_angular_differences(MW_diff_angles, "MW", name_dphi, num_dphi,
			      name_para, axe, NM.N_axes, seq1);

    char nameout[100]; int a;
    sprintf(nameout, "%s_mode_projections.dat", name_para);
    printf("Writing mode projections in %s\n", nameout);
    file_out=fopen(nameout, "w");
    fprintf(file_out, "# %d modes %d methods for deriving angles\n",
	    NM.N, num_dphi);
    for(i=0; i<num_dphi; i++)
      fprintf(file_out, "# Method %d: %s\n", i, name_dphi[i]);
    fprintf(file_out, "#mode ");
    for(i=0; i<num_dphi; i++)fprintf(file_out, " %d=dphi%d", i+2, i+1);
    fprintf(file_out, " 1/omega^2\n");
    for(a=0; a<NM.N; a++){
      fprintf(file_out, "%d ", a);
      for(i=0; i<num_dphi; i++){
	fprintf(file_out, " %.4f", c2_alpha[i][a]);
      }
      fprintf(file_out, " %.4f", 1./NM.omega2[a]);
      fprintf(file_out, "\n");
    }
    fclose(file_out);

  }
  for(i=0; i<num_dphi; i++){
    free(c2_alpha[i]);
    free(diff_angles[i]);
    free(name_dphi[i]);
  }

  return(rmsd);
}

void Print_diff(int *num_dphi, char **name_dphi,
		float **diff_angles, float **MW_diff_angles,
		float **c2_alpha, int *outlier_tors,
		struct Tors *Diff,
		struct Jacobian *J,
		struct Normal_Mode *NM, float M_sqrt, float thr, 
		struct bond *bonds, atom *atoms1, int natoms1,
		struct Reference Ref, float *coord_ref_1, float *coord_ref_2,
		struct residue *seq1, struct axe *axe,
		int nstruct, char *name,
		FILE *file_out, FILE *file_pdb)
{
  name_dphi[*num_dphi]=malloc(80*sizeof(char));
  strcpy(name_dphi[*num_dphi], name);
  fprintf(file_out, "##### %s\n", name);

  // Generate new conformation from torsion angles
  Set_bonds_measure(bonds, natoms1, atoms1);
  Build_up(bonds, natoms1, Diff->Tors, NM->N_axes);
  float coord_all[3*natoms1];
  Put_coord(coord_all, bonds, natoms1);
  float coord_new[3*Ref.N_ref], RMSD=0, RMSD_target=0;
  Write_ref_coord(coord_new, Ref.N_ref, coord_all, Ref.atom_num);
  RMSD_target=
    rmsd_mclachlan_f(coord_ref_2, coord_new, Ref.mass_atom, Ref.N_ref);
  RMSD=rmsd_mclachlan_f(coord_ref_1, coord_new, Ref.mass_atom, Ref.N_ref);
  printf("RMSD= %.2f\n", RMSD_target);

  // Generate new conformation from linear approximation
  float RMSD_lin=0, RMSD_lin_target=0, RMSD_lin_full=0;
  RMSD_lin=RMSD_LIN(&RMSD_lin_target, &RMSD_lin_full,
		    J, Diff->Tors, Ref,
		    coord_ref_1, coord_ref_2, coord_all);

  // Remove 2pi
  Periodic_angles(Diff->Tors, axe, NM->N_axes);
  Compute_MW(Diff, J);


  // Print PDB
  if(file_pdb){
    fprintf(file_pdb, "MODEL %d  %s\n", *num_dphi+1, name);
    float coord_1[3*natoms1], mass[natoms1]; int i, j; float *c=coord_1;
    for(i=0; i<natoms1; i++){
      mass[i]=1; double *r=atoms1[i].r;
      for(j=0; j<3; j++){*c=*r; c++; r++;}
    }
    rmsd_mclachlan_f(coord_1, coord_all, mass, natoms1);
    Print_PDB(file_pdb, atoms1, natoms1, coord_all, seq1, *num_dphi, RMSD);
  }

  int N_modes=NM->N, i, k, a;

  // Project on normal modes
  float *c2=NM->confchange2;
  double sum=0, Delta_E=0;
  for (i=0;i<N_modes;i++){
    float d=Scalar_product(Diff->MW_Tors, NM->MW_Tors[i], NM->N_axes);
    //Scalar_product_weighted(Diff->Cart,NM->Cart[i],Ref.mass_coord,Ref.N_cart);
    // Penalize modes that contribute to conf.change as little as noise
    //d*=Fermi_function(d/M_sqrt, thr, 50.0);
    float d2=d*d;
    c2[i]=d2; sum+=d2;
    Delta_E+=NM->omega2[i]*d2;   
    if(NM->select[i]==0){Diff->coeff[i]=0;}
    else{Diff->coeff[i]=d;}
  }
  for(i=0;i<N_modes;i++)c2[i]/=sum; //norm_c2; // Normalize 


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

  double norm_c2=
    Scalar_product_weighted(Diff->Cart,Diff->Cart,Ref.mass_coord,Ref.N_cart);
  float M_tot=M_sqrt*M_sqrt, rmsd=sqrt(norm_c2/M_tot);

  // Conformation change
  int nn=0, n_outlier=0; float THR=0.5;
  for(a=0; a<NM->N_axes; a++){
    if(outlier_tors && (outlier_tors[a])){n_outlier++; continue;}
    if(fabs(Diff->Tors[a])>THR)nn++;
  }
  Diff->Tors_frac=Tors_fraction(Diff, Ref.mass_coord);
  Diff->RMSD_Tors=RMS(Diff->Tors, Diff->N_axes, outlier_tors);

  if(nstruct==0){
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
  }
  float Coll=Collectivity_norm2(Diff->Cart, Ref.N_cart)/Ref.N_cart;
  fprintf(file_out, "Coll.(cart)           %.3f\n", Coll);
  Coll=Collectivity_norm2(Diff->MW_Tors, NM->N_axes)/NM->N_axes;
  fprintf(file_out, "Coll.(MWtors)         %.3f\n", Coll);
  //Coll=Collectivity_norm2(Diff->Tors, NM->N_axes)/NM->N_axes;
  Coll=Collectivity_norm2_outlier(Diff->Tors, NM->N_axes, outlier_tors);
  Coll/=(NM->N_axes-n_outlier);
  fprintf(file_out, "Coll.(tors)           %.3f\n", Coll);
  fprintf(file_out, "#\n");

  // Correlation angular fluctuations - conformation change
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

  // Projection on normal modes and frequency
  Delta_E/=(2.0*M_sqrt*M_sqrt);
  fprintf(file_out, "Delta_E(kBT)          %.3g\n", Delta_E);
  fprintf(file_out, "Delta_E(kBT)/Dr^2     %.3g\n", Delta_E/(rmsd*rmsd));

  //rho1
  float x[N_modes], y[N_modes];  int imax=-1;
  float cmax=-1; k=0; 
  for(i=0; i<N_modes; i++){
    if(NM->select[i]==0)continue;
    if(NM->Cart_coll[i] < Coll_thr_cc)continue;
    x[k]=NM->contr2fluct[i]; y[k]=c2[i]; k++;
    if(c2[i]>cmax){cmax=c2[i]; imax=i;}
  }
  r= Corr_coeff(x, y, k, &slope, &offset);
  //*Confchange_mode=imax;

  fprintf(file_out, "r[c^2,1/w^2]          %.3f\n", r);
  fprintf(file_out, "slope                 %.2f\n", slope);
  fprintf(file_out, "offset                %.4f\n", offset);
  Coll=Collectivity_norm1(y, k);
  //Coll=Collectivity_norm1(c2, NM->N_relevant);
  fprintf(file_out, "Recp.Coll(cc)         %.1f\n", Coll);

  /*if(nstruct==0){
    int na=Ref.N_cart/3;
    fprintf(file_out, "Recp.Coll.(cc)/Natm   %.3f\n", Coll/na);
    fprintf(file_out, "R.Coll.(cc)Renyi/Natm %.3f\n",
      Collectivity_Renyi_norm1(c2,N_modes)/na);
      }*/
  /*fprintf(file_out, "Area(cc)              %.3f\n",
    Area_norm1(c2,N_modes));*/

  if(nstruct==0){
    fprintf(file_out, "Most_contr.mode(cc):  %d %.3f %.3f  %.2f %.2f %.2f\n",
	    imax, cmax, NM->contr2fluct[imax], NM->Cart_coll[imax],
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
     x[k]=log(NM.contr2fluct[i]); y[k]=log(c2[i]); k++;
     if(c2[i]>cmax){cmax=c2[i]; imax=i;}
   }
   }*/

  if(nstruct==0)fprintf(file_out, "#\n");

  // rho
  k=0; cmax=-1;
  //~ float min_pred=-2*mu/slope;
  for(i=0; i<N_modes; i++){
    if((NM->select[i]==0)||(NM->contr2fluct[i]==0))continue;
    if(NM->Cart_coll[i] < Coll_thr_cc)continue;
    x[k]=NM->contr2fluct[i];
    y[k]=(c2[i])/NM->contr2fluct[i];
    if(y[k]>cmax){cmax=y[k]; imax=i;}  k++;
  }
  r= Corr_coeff(x, y, k, &slope, &offset);
  //~ float k_Therm=Collectivity_norm1(NM->contr2fluct, N_modes);
  if(nstruct==0)
    fprintf(file_out, "Most_contr.mode(c*w): %d %.3f %.3f  %.2f %.2f %.2f\n",
	    imax, sqrt(cmax), NM->contr2fluct[imax], NM->Cart_coll[imax],
	    NM->MW_Tors_coll[imax], NM->Tors_coll[imax]);
  fprintf(file_out, "r[(c*w)^2,1/w^2]      %.3f\n", r);
  fprintf(file_out, "Significance(Z)       %.3f nsel=%d\n", r*sqrt(k), k);
  //fprintf(file_out, "Significance(Z)       %.3f\n",r*sqrt(k_Therm));
  Coll=Collectivity_norm1(y, k);
  fprintf(file_out, "Recp.Coll((cc*w)^2)   %.1f\n", Coll);

  // Probability of native fluctuation
  /*double logP=0;
  for(i=0; i<N_modes; i++){
    if((NM->select[i])&&(c2[i])){
      logP+=NM->omega2[i]*Diff->coeff[i]*Diff->coeff[i];
    }
  }
  logP/=(2*M_tot);
  fprintf(file_out, "log(native probab.)   %.4f\n", logP);*/

  // Correlation conformation change and anharmonicity
  if(NM->Anharmonicity[0]){
    float P2[N_modes], norm=norm_c2/sum;
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
  sum=0; k=0; float F_coeff[N_modes];
  for(i=0; i<N_modes; i++){
    if((NM->select[i])&&(NM->contr2fluct[i])){
      float c=Diff->coeff[i]*NM->omega2[i], cc=c*c;
      F_coeff[i]=cc;
      sum+=cc;
      //x[k]=log(NM->contr2fluct[i]); y[k]=log(cc);
      if(cc>cmax){cmax=cc; imax=i;} k++;
    }else{
      F_coeff[i]=0;
    }
  }
  //Corr_coeff(x, y, k, &slope, &offset);
  //fprintf(file_out, "r_ll(therm,force)     %.3f\n", r);
  //fprintf(file_out, "expo(therm,force)     %.2f\n", slope);

  for(i=0; i<NM->N; i++)F_coeff[i]/=sum;
  float k_Force_mode=Collectivity_norm1(F_coeff, N_modes);
  fprintf(file_out, "Recp.Coll.(Force)     %.1f\n", k_Force_mode);
  fprintf(file_out, "Most_contr.mode(force) %d %.3g %.2g  %.2f %.2f %.2f\n",
	  imax, cmax/sum, NM->contr2fluct[imax], NM->Cart_coll[imax],
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
  Allocate_tors(&Diff2, NM.N_axes, NM.N_cart, NM.N);
  for(i=0; i<NM.N_cart;  i++){
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
			    int ANM, float M_sqrt, float rmsd)
{
  /***************** Normal modes confchange ********************/

  char nameout[400];
  sprintf(nameout, "%s_Modes.dat", nameout2);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  fprintf(file_out, "# Conf.change %s\n", name2);

  int rank, ik,  k; 

  // Compute confchange and copy on normal modes

  // Sort
  float *c2=NM.confchange2;
  /*float c2[NM.N_relevant]; double norm=0;
  for (k=0;k<NM.N_relevant; k++){
    float c=Diff.coeff[k]; c2[k]=c*c; norm+=c2[k];
  }
  for (k=0;k<NM.N_relevant; k++)c2[k]/=norm;*/

  int *sorted=malloc(NM.N*sizeof(int));
  for (k=0;k<NM.N_relevant; k++)sorted[k]=0;
  for(rank=0; rank<NM.N_relevant; rank++){
    ik=-1; float cmax=-1;
    for (k=0;k<NM.N_relevant; k++){
      if((NM.select[k])&&(sorted[k]==0)&&(c2[k]>=cmax)){
	cmax=c2[k]; ik=k;
      }
    }
    NM.sort[rank]=ik;
    if(ik>=0)sorted[ik]=1;
  }
  free(sorted);

  // Print lines
  fprintf(file_out, "# mode\tContr_to_confchange\tCumulative");
  fprintf(file_out, "\tOmega^(-2) RMSD ");
  fprintf(file_out, "\tColl.(cart)\tColl.(MWtors)\tColl.(tors)");
  if(ANM)fprintf(file_out, "\tTors_frac"); //force_coeff
  if(NM.Anharmonicity[0]){
    fprintf(file_out, "\tAnharm_Str Str_confchange");
    fprintf(file_out, "\tAnharm_Ene\tEne_confchange\tMax_RMSD");
  }
  fprintf(file_out, "\n");

  double sum_coeff_C=0; //sum_B=0;
  for(rank=0; rank<NM.N_relevant; rank++){
    ik=NM.sort[rank]; if(ik<0)continue;
    sum_coeff_C+=c2[ik];
    //sum_B+=NM.contr2fluct[ik];
    fprintf(file_out,
	    "%d\t%8.3g\t%.4f\t%7.4g\t%8.3g\t%.3f\t%.3f\t%.3f",
	    ik, c2[ik], sum_coeff_C, //NM.contr2fluct[ik],
	    1./NM.omega2[ik],
	    1./(M_sqrt*NM.omega[ik]),
	    NM.Cart_coll[ik], NM.MW_Tors_coll[ik],
	    NM.Tors_coll[ik]); //Force.coeff[ik], sum_B,
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
  float q, q1, r, slope, offset, qmax=0, rmax=0;
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
  for(i=0; i<N; i++)if(wminus2[i]<wmin)wmin=wminus2[i]; wmin-=0.0001;
  for(q=1; q<=3.5; q+=0.25){
    if(q==1){for(i=0; i<N; i++)y[i]=log(P[i]);}  // q=1: Shannon
    else{q1=q-1; for(i=0; i<N; i++)y[i]=pow(P[i], q1);}
    r= Corr_coeff(wminus2, y, N, &slope, &offset);
    double d_KL=0, q2=1/q1, Z=0, lnp;
    if(q==1){
      for(i=0; i<N; i++){
	d_KL+=P[i]*(y[i]-slope*wminus2[i]);
	Z+=exp(slope*wminus2[i]);
      }
    }else{
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
  for(i=0; i<N; i++)P[i]/=P1; P2=P2/(P1*P1)-1./N; dP=dP/P1-log(P1);
  //printf("P2= %.3f -S[P]=%.3f N=%d\n", P2, dP, N);
  for(q=0.5; q<=8; q+=0.25){
    double  q1=q/2., x1=0, d_KL=0, d2=0, d;
    for(i=0; i<N; i++){x[i]=pow(wminus2[i], q1); x1+=x[i];}
    for(i=0; i<N; i++)d_KL+=P[i]*log(x[i]); d_KL=dP-(d_KL-log(x1));
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
void Print_angular_differences(float **diff_angles, char *type,
			       char **name_dphi, int num_dphi,
			       char *name_para, struct axe *axe,
			       int N_axes, struct residue *seq1)
{
  char nameout[100]; int a, i;
  sprintf(nameout, "%s_%sangular_differences.dat", name_para, type);
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

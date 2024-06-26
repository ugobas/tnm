#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "nma_para.h"
#include "buildup.h"
#include "simulation.h"
#include "atoms.h"
#include "vector.h"
#include "random3.h"
#include "McLachlan.h"
#include "align_tnm.h"
#include "kinetic_tnm.h"
#include "dof_tnm.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "choldc.h"
#include "allocate.h"
#include "optimization.h"
#include "externals.h"
#include "interactions_tnm.h"
#include "anharmonic_tnm.h"
#include "unfolding.h"

int N_int_max;

float Unfold(// Modified in the routine:
	     float *Ene_T, float *Rg2_T,double **L_Hess_u,
	     int sec_der, // L_Hess is modified onlyn if sec_der=1
	     struct bond *bonds, atom *atoms, int natoms,
	     float *coord, int N_Cart,
	     double *dphi_all_opt, struct axe *axe,
	     // Given from input and not changed:
	     int naxes, int nmain,
	     struct chain *chains, int Nchain,
	     struct residue *seq, int nres,
	     struct interaction *Int_list, int N_int,
	     struct interaction **Int_KB, int NA, int K_tors,
	     struct Reference Ref,
	     float *coord_ini, float *mass, float Mass_tot,
	     FILE *file_fold, float Temp);

void Torsional_derivatives(double **Hessian, double *Force,
			   struct inter_new *New_inter, int N_int_new,
			   atom *atoms, int natoms, struct axe *axe,
			   int nmain, int naxes, double *dphi,
			   struct chain *chains, int Nchain, int sec_der);
void Unfolding_derivatives(double **Hessian_u, double *Force_u,
			   float *mass, atom *atoms, int natoms,
			   struct axe *axe, int naxes,
			   float Temp, float Mass_tot, int sec_der);
void Sum_matrix_d(double **M_sum, double **M, int n1, int n2);
void Sum_vector_d(double *v_sum, double *v, int n);

float Compute_Rg2(float *coord, float *mass, int natoms, float Mass_tot);
void Mass_axes(struct axe *axes, int naxes, atom *atoms, int natoms,
	       struct Reference Ref);
void Move_atoms(atom *atoms, int natoms, float *coord);
float RMS(double *ff, int n);
extern float Empty_inertia(double *mr_tot, double **corr_sum);
extern void Print_seqres(FILE *file_out, struct chain *chains, int Nchain);

int Unfolding(struct bond *bonds_in, atom *atoms_in, int natoms_in,
	      struct axe *axes_in, int naxes_in, int nmain_in,
	      struct chain *chains_in, int Nchain_in,
	      struct residue *seq_in, int nres_in,
	      struct interaction *Int_list_in, int N_int_in,
	      struct interaction **Int_KB_in, int NA_in,
	      float K_TORS, char *file_pdb_in)  
{
  printf("Simulating torsional unfolding\n");
  int K_tors=K_TORS;

  // Coordinates
  atom *atoms; int natoms; int i;
  struct chain *chains; int Nchain;
  struct residue *seq; int nres;
  char pdbid[20];
  if(atoms_in){
    atoms=atoms_in; natoms=natoms_in;
    seq=seq_in; nres=nres_in; chains=chains_in; nres=nres_in;
    GetPdbId(file_pdb_in, pdbid);
  }else{
    int ANISOU=0, nmr=0, n_lig=0, na_lig=0;
    atom *ligatom=NULL;  struct residue *ligres=NULL;
    float TEMPERATURE=-1; char chain[40];
    Read_PDB(&nres, &seq, &n_lig, &ligres, chain, &ANISOU, &nmr,
	     &natoms, &atoms, &na_lig, &ligatom, &chains, &Nchain,
	     &TEMPERATURE, pdbid, file_pdb_in);
  }

  // File names
  char nameout[80], nameprot[80];
  sprintf(nameprot, "%s", pdbid);
  {char tmp[10];
    for(i=0; i<Nchain; i++){
      if((chains[i].label!=' ')&&(chains[i].label!='\0')){
	sprintf(tmp, "%c", chains[i].label);
	strcat(nameprot, tmp);
      }
    }
  }
  strcpy(nameout, nameprot);

  // Bonds
  struct bond *bonds_ini;
  if(bonds_in==NULL){
    bonds_ini=Set_bonds_topology(natoms, atoms, seq);
  }else{
    bonds_ini=bonds_in; Set_bonds_measure(bonds_ini, natoms, atoms);
  }
  struct bond bonds[natoms];
  Copy_bonds(bonds, bonds_ini, natoms);

  /*************************  Interactions  *************************/
  struct interaction *Int_list; int N_int;
  if(Int_list_in==NULL){
    KAPPA=8.2; C_THR=4.5; POW=1; EXP_HESSIAN=6; C0=3.5; 
    Compute_interactions(&N_int, &Int_list, "MIN",
			 atoms, natoms, nres, nameprot);
    K_PHI=K_TORS; K_PSI=K_TORS; K_OMEGA=0; K_CHI=0; K_BA=0; K_BL=0;
  }else{
    Int_list=Int_list_in; N_int=N_int_in;
  }
  // Knowledge-based interactions
  struct interaction **Int_KB; int NA; // amino acids or atom types
  if(Int_KB_in){
    NA=NA_in; Int_KB=Int_KB_in;
  }else{
    NA=20;
    Int_KB=malloc(NA*sizeof(struct interaction *));
    for(i=0; i<NA; i++)Int_KB[i]=malloc(NA*sizeof(struct interaction));
    Assign_interactions_KB(Int_KB, NA, POW, EXP_HESSIAN, KAPPA);
    int aa_seq[nres]; for(i=0; i<nres; i++)aa_seq[i]=seq[i].i_aa;
    Assign_interactions_SB(Int_list, N_int, Int_KB, aa_seq);
  }

  // Reference atoms: all heavy atoms
  struct Reference Ref;
  int N_ref=Set_reference(&Ref, 0, "ALL", atoms, 0, natoms);
  if(N_ref!=natoms){
    printf("ERROR, discordant number of atoms N_ref=%d natoms=%d\n",
	   N_ref, natoms); exit(8);
  }


  // Internal degrees of freedom
  struct axe *axes_ini; int naxes, nmain;
  if(axes_in){
    axes_ini=axes_in; naxes=naxes_in; nmain=nmain_in;
  }else{
    struct rigid *Rigid_dof=NULL; int N_rigid=0, nskip=0, N_diso=0;
    int OMEGA=0, SIDECHAINS=0, PSI=1, MIN_INT_MAIN=1, MIN_INT_SIDE=1;
    int last_ali_res=nres;
    axes_ini=Set_DegofFreed(&naxes, &nmain,&nskip,&N_diso,&Rigid_dof,&N_rigid,
		      bonds, atoms, natoms, Ref.atom_num, N_ref,
		      seq, last_ali_res, chains, Nchain,
		      Int_list, N_int, MIN_INT_MAIN, MIN_INT_SIDE,
		      OMEGA, SIDECHAINS, PSI);
    printf("Number of degrees of freedom: %d\n", naxes);
    if(naxes!=nmain)printf("Side chains dof: %d\n", naxes-nmain);
  }
  struct axe axe[naxes]; for(i=0; i<naxes; i++){axe[i]=axes_ini[i];}
  // Torsion angles
  double phi_ini[naxes]; //, phi_barrier[naxes];
  printf("Internal coordinates of initial structure\n");
  Internal_coordinates(phi_ini, naxes, bonds_ini, natoms);
  double dphi_sum[naxes], dphi[naxes];
  for(i=0; i<naxes; i++){dphi[i]=0; dphi_sum[i]=0;}

  // Cartesian degrees of freedom
  int N_Cart=3*natoms;
  float coord_ini[N_Cart], coord[N_Cart], coord_pdb[N_Cart];
  Put_coord(coord_ini, bonds_ini, natoms);
  for(i=0; i<N_Cart; i++){
    coord[i]=coord_ini[i]; coord_pdb[i]=coord[i];
  }
  double Mass_tot=0; float mass[natoms];
  for(i=0; i<natoms; i++){mass[i]=Mass(atoms+i); Mass_tot+=mass[i];}

  // Allocate
  double **Hessian=Allocate_mat2_d(naxes, naxes);
  double **Hessian_u=Allocate_mat2_d(naxes, naxes);
  double **L_Hess_u=Allocate_mat2_d(naxes, naxes);
  double Force[naxes], Force_u[naxes];

  // Equilibrium values
  int Ntemp=100; // Ntemp temperature values with step dTemp
  float dTemp=0.02, Temp=0;
  float Rg2_T[Ntemp]; Rg2_T[0]=Compute_Rg2(coord_ini, mass, natoms, Mass_tot);
  float Ene_T[Ntemp]; Ene_T[0]=0;
  //float G_T[Ntemp]; G_T[0]=Ene_T[0]-Temp*Rg2_T[0];
  float conv_T[Ntemp]; conv_T[0]=1.0;
  float RMSD_T[Ntemp]; RMSD_T[0]=0;
  float Temp_T[Ntemp]; Temp_T[0]=0;
 
  // Print equilibrium
  char namefold[200]; sprintf(namefold, "%s.unfolding.dat", nameout);
  FILE *file_fold=fopen(namefold, "w");
  printf("Opening %s\n", namefold);
  fprintf(file_fold,"# Unfolding the protein with linear torsion response\n");
  fprintf(file_fold,"# 1=T 2=RMSD 3=E 4=Rg^2 5=G=E-TRg^2 "
	  "6=convergence 7=dE/dT 8=dRg^2/dT\n");

  // Print angles
  char namephi[200]; sprintf(namephi, "%s.torsion.dat", nameout);
  FILE *file_phi=fopen(namephi, "w");
  printf("Opening %s\n", namephi);
  fprintf(file_phi, "#1=index 2=delta_phi 3=type\n");

  // Opening file for printing PDB
  char pdbout[200];
  sprintf(pdbout, "%s_unfolding.pdb", nameout);
  FILE *file_pdb=fopen(pdbout, "w");
  printf("Writing trajectory in PDB format in %s\n", pdbout);
  Print_seqres(file_pdb, chains, Nchain);

  // Initialization of energy parameters
  N_int_max=30*nres; int N_int_new=0;
  struct inter_new New_inter[N_int_max];

  int sec_der=0; // Recompute second derivative at each step?
  float RMSD_THR=0.35;
  // Step after which PDB files are printed and Hessian is recomputed
  int it; int k_pdb=0; float rmsd_pdb=0;
  for(it=0; it<(Ntemp-1); it++){

    if(it==0 || rmsd_pdb>RMSD_THR){
      k_pdb++;
      for(i=0; i<N_Cart; i++)coord_pdb[i]=coord[i];
      Print_PDB(file_pdb,atoms,natoms,coord_pdb,seq,k_pdb,RMSD_T[it]);
      N_int_all=
	Interactions_CA(Int_list_all,d_CA_high,atoms,natoms,nres,coord);

      if(sec_der==0){
	printf("Computing axes, interactions and Hessian\n");
	Ene_T[it]=Energy_anharmonic_new(New_inter, &N_int_new, N_int_max,
					coord, atoms, natoms, seq, nres,
					Int_list, N_int, Int_KB, NA,
					axe, naxes, dphi, K_tors, 1);
	Torsional_derivatives(Hessian, Force, New_inter, N_int_new,
			      atoms, natoms, axe, nmain, naxes, dphi,
			      chains, Nchain, 1);
	Unfolding_derivatives(Hessian_u, Force_u, mass,
			      atoms, natoms, axe, naxes, Temp, Mass_tot, 1);
	Sum_matrix_d(Hessian_u, Hessian, naxes, naxes);
	Sum_vector_d(Force_u, Force, naxes);
	choldc(L_Hess_u, Hessian_u, naxes);
      }
    }

    Temp+=dTemp; printf("Unfolding at T=%.3f\n", Temp);

    int tt=it+1; // Compute next T
    for(i=0; i<naxes; i++)dphi[i]=0;
    conv_T[tt]=Unfold(Ene_T+tt, Rg2_T+tt, L_Hess_u, sec_der,
		      bonds, atoms, natoms, coord, N_Cart, dphi, axe,
		      naxes, nmain, chains, Nchain, seq, nres,
		      Int_list, N_int, Int_KB, NA, K_tors,
		      Ref, coord_ini, mass, Mass_tot, file_fold, Temp);
    RMSD_T[tt]=rmsd_mclachlan_f(coord_ini, coord, mass, natoms);
    rmsd_pdb  =rmsd_mclachlan_f(coord_pdb, coord, mass, natoms);

    float dE, dR2;
    if(it){
      dE=(Ene_T[tt]-Ene_T[it-1])/2; dR2=(Rg2_T[tt]-Rg2_T[it-1])/2;
    }else{
      dE=(Ene_T[tt]-Ene_T[it]); dR2=(Rg2_T[tt]-Rg2_T[it]);
    }
    fprintf(file_fold,"%.3f\t%.3f\t%.3g\t%.3g\t%.2g\t%.3g\t%.3g\t%.3g\n",
	    Temp_T[it], RMSD_T[it], Ene_T[it], Rg2_T[it],
	    Ene_T[it]-Temp_T[it]*Rg2_T[it], conv_T[it], dE/dTemp, dR2/dTemp);
    fprintf(file_phi, "# Temp= %.3f\n", Temp_T[it]);
    for(i=0; i<naxes; i++){
      fprintf(file_phi, "%d %.4f %c\n", i, dphi_sum[i], axe[i].type);
      dphi_sum[i]+=dphi[i];
    }
    fprintf(file_phi, "&\n");
  }

  float dE=(Ene_T[it]-Ene_T[it-2])/2, dR2=(Rg2_T[it]-Rg2_T[it-2])/2;
  fprintf(file_fold,"%.3f\t%.3f\t%.3g\t%.3g\t%.2g\t%.3g\t%.3g\t%.3g\n",
	  Temp_T[it], RMSD_T[it], Ene_T[it], Rg2_T[it],
	  Ene_T[it]-Temp_T[it]*Rg2_T[it], conv_T[it], dE/dTemp, dR2/dTemp);
  fprintf(file_phi, "# Temp= %.3f\n", Temp_T[it]);
  for(i=0; i<naxes; i++){
    fprintf(file_phi, "%d %.4f %c\n", i, dphi_sum[i], axe[i].type);
  }
  fprintf(file_phi, "&\n");

  fclose(file_phi);
  fclose(file_pdb);
  fclose(file_fold);
  printf("Writing PDB files in %s\n", pdbout);
  printf("Writing equilibrium values in %s\n", namefold);
  printf("Writing torsion angles in %s\n", namephi);

  // restore previous
  Move_atoms(atoms, natoms, coord_ini);
  Set_rot_shift(axe, naxes, atoms, natoms);
  Mass_axes(axe, naxes, atoms, natoms, Ref);

  Empty_matrix_d(Hessian, naxes);
  Empty_matrix_d(Hessian_u, naxes);
  Empty_matrix_d(L_Hess_u, naxes);
  return(0);
}

float Compute_Rg2(float *coord, float *mass, int natoms, float Mass_tot){
  double Rg2=0; float *r=coord;
  for(int i=0; i<natoms; i++){
    float r2=0; for(int j=0; j<3; j++){r2+=(*r)*(*r); r++;}
    Rg2+=mass[i]*r2;
  }
  return(Rg2/Mass_tot);
}

float Unfold(// Modified in the routine:
	     float *Ene_T, float *Rg2_T,double **L_Hess_u,
	     int sec_der, // L_Hess is modified onlyn if sec_der=1
	     struct bond *bonds, atom *atoms, int natoms,
	     float *coord, int N_Cart,
	     double *dphi_all_opt, struct axe *axe,
	     // Given from input and not changed:
	     int naxes, int nmain,
	     struct chain *chains, int Nchain,
	     struct residue *seq, int nres,
	     struct interaction *Int_list, int N_int,
	     struct interaction **Int_KB, int NA, int K_tors,
	     struct Reference Ref,
	     float *coord_ini, float *mass, float Mass_tot,
	     FILE *file_fold, float Temp)
{
  int IT_MAX=100;
  float A_MAX=0.1;
  int nfail_max=3;

  float Ene, Rg2, FF, F_opt, conv=100, conv_opt;
  int it_opt=-1, nfail=0, i;
  double dphi[naxes], dphi_all[naxes];
  for(i=0; i<naxes; i++){dphi[i]=0; dphi_all[i]=0;}
  struct bond bonds_tmp[natoms];
  Copy_bonds(bonds_tmp, bonds, natoms);

  double Force[naxes], Force_u[naxes];
  double **Hessian=NULL, **Hessian_u=NULL;
  if(sec_der){
    Hessian=Allocate_mat2_d(naxes, naxes);
    Hessian_u=Allocate_mat2_d(naxes, naxes);
  }

  int N_int_new=0;
  struct inter_new New_inter[N_int_max];

  int IT1=IT_MAX-1;
  for(int iter=0; iter<IT_MAX; iter++){

    Ene=Energy_anharmonic_new(New_inter, &N_int_new, N_int_max,
			      coord, atoms, natoms, seq, nres,
			      Int_list, N_int, Int_KB, NA,
			      axe, naxes, dphi_all, K_tors, sec_der);
    Rg2=Compute_Rg2(coord, mass, natoms, Mass_tot);
    FF=Ene-Temp*Rg2;
    float RMSD=rmsd_mclachlan_f(coord_ini, coord, mass, natoms);

    fprintf(file_fold,"%.3f\t%.3f\t%.3g\t%.3g\t%.2g\t%.3g\t%.3g\t%.3g\n",
	    Temp, RMSD, Ene, Rg2, FF, conv, 0.0, 0.0);
    fflush(file_fold);
    
    if(it_opt<0 || FF<F_opt){
      it_opt=iter; nfail=0; F_opt=FF;
      *Ene_T=Ene; *Rg2_T=Rg2; conv_opt=conv;
      for(i=0; i<naxes; i++){dphi_all_opt[i]=dphi_all[i];}
    }else{
      nfail++;
      if(nfail>=nfail_max)break;
    }
    if(iter==IT1)break;

    Torsional_derivatives(Hessian, Force, New_inter, N_int_new,
			  atoms, natoms, axe, nmain, naxes, dphi,
			  chains, Nchain, sec_der);
    Unfolding_derivatives(Hessian_u, Force_u, mass,
			  atoms, natoms, axe, naxes, Temp, Mass_tot, sec_der);

    conv=RMS(Force_u, naxes);

    // Compute dphi=(H)^(-1)Force
    // Symmetric matrix inversion: Solve (L L^t)dphi = Force
    if(sec_der)choldc(L_Hess_u, Hessian_u, naxes);
    double Z[naxes];
    Forward_substitution(Z, L_Hess_u, Force_u, naxes);
    Backward_substitution(dphi, L_Hess_u, Z, naxes);

    // Set scale
    float dphi_max=0;
    for(i=0; i<naxes; i++)
      if(fabs(dphi[i])>dphi_max)dphi_max=fabs(dphi[i]);
    if(dphi_max>A_MAX){
      float f=A_MAX/dphi_max; for(i=0; i<naxes; i++)dphi[i]*=f;
    }

    // Move coordinates
    for(i=0; i<naxes; i++)dphi_all[i]+=dphi[i];
    Build_up(bonds_tmp, natoms, dphi, naxes);
    Put_coord(coord, bonds_tmp, natoms);

    Move_atoms(atoms, natoms, coord);
    Set_rot_shift(axe, naxes, atoms, natoms);
    Mass_axes(axe, naxes, atoms, natoms, Ref);
  }

  if(it_opt!=IT1){
    for(i=0; i<naxes; i++)dphi[i]=dphi_all_opt[i];
    Build_up(bonds, natoms, dphi, naxes);
    Put_coord(coord, bonds, natoms);
    Move_atoms(atoms, natoms, coord);
    Set_rot_shift(axe, naxes, atoms, natoms);
    Mass_axes(axe, naxes, atoms, natoms, Ref);
  }else{
    Copy_bonds(bonds, bonds_tmp, natoms);
  }

  if(Hessian){
    Empty_matrix_d(Hessian, naxes);
    Empty_matrix_d(Hessian_u, naxes);
  }

  return(conv_opt);
}

void Torsional_derivatives(double **Hessian, double *Force,
			   struct inter_new *New_inter, int N_int_new,
			   atom *atoms, int natoms, struct axe *axe,
			   int nmain, int naxes, double *dphi,
			   struct chain *chains, int Nchain, int sec_der){

}

void Unfolding_derivatives(double **Hessian_u, double *Force_u,
			   float *mass, atom *atoms, int natoms,
			   struct axe *axe, int naxes,
			   float Temp, float Mass_tot, int sec_der){

}

void Sum_matrix_d(double **M_sum, double **M, int n1, int n2){
  for(int i=0; i<n1; i++){
    double *M1=M_sum[i], *M2=M[i];
    for(int j=0; j<n2; j++){*M1+=(*M2); M1++; M2++;}
  }
}

void Sum_vector_d(double *v_sum, double *v, int n){
  double *v1=v_sum, *v2=v;
  for(int i=0; i<n; i++){*v1+=(*v2); v1++; v2++;}
}

void Mass_axes(struct axe *axes, int naxes, atom *atoms, int natoms,
	       struct Reference Ref)
{
  // Compute inertia tensor and center of mass for all atoms
  int i, j;
  double **corr_sum=Allocate_mat2_d(3,3);
  double mr_tot[3]; for(i=0; i<3; i++)mr_tot[i]=0;
  double M_tot=Sum_inertia(mr_tot, corr_sum, atoms, Ref, 0, natoms-1);
  double **inertia_tot=Allocate_mat2_d(3,3);
  Inertia_tensor(inertia_tot, corr_sum, mr_tot, M_tot);

  // Partial inertia tensor for all axes
  // Exploits the fact that axes are nested
  double M_part=0, mr_part[3], **inertia_part=Allocate_mat2_d(3,3);
  Empty_inertia(mr_part, corr_sum);

  struct axe *axe=axes+naxes-1;
  int last_kin=-1, i_last=0;
  for(int k=naxes-1; k>=0; k--){
    // Compute partial center of mass and inertia tensor
    if(axe->last_kin!=last_kin){
      // New chain, initialize counters
      M_part=Empty_inertia(mr_part, corr_sum);
      last_kin=axe->last_kin;
      i_last=last_kin;
    }
    M_part+=Sum_inertia(mr_part, corr_sum, atoms, Ref, axe->first_kin, i_last);
    Inertia_tensor(inertia_part, corr_sum, mr_part, M_part);
    i_last=axe->first_kin-1;

    // Store kinetic variables
    axe->mass=M_part;
    for(i=0; i<3; i++){
      axe->MR[i]=mr_part[i];
      for(j=0; j<3; j++)axe->I[i][j]=inertia_part[i][j];
    }

    // Test that at least some atom is moved
    if(M_part==0){
      printf("WARNING, axis %d does not move any reference atom!\n", k);
      printf("Reference atoms: %d - %d\n", axe->first_kin, axe->last_kin);
      printf("All atoms: %d - %d\n", axe->first_atom, axe->last_atom);
      printf("dof: %d-%d (main) ", axe->first_main, axe->last_main);
      if(axe->first_side>=0)printf(" %d-%d (side)", axe->first_side, k);
      printf("\n");
    }

    // New axe
    axe--;
  }
  // Cleaning
  Empty_matrix_d(corr_sum, 3);
  Empty_matrix_d(inertia_part, 3);
  Empty_matrix_d(inertia_tot, 3);
}


void Move_atoms(atom *atoms, int natoms, float *coord){
  float *r=coord;
  for(int i=0; i<natoms; i++){
    double *ri=atoms[i].r;
    for(int j=0; j<3; j++){*ri=*r; r++; ri++;}
  }
}

float RMS(double *ff, int n){
  double sum=0, *f=ff;
  for(int i=0; i<n; i++){sum+=(*f)*(*f); f++;}
  return(sqrt(sum/n));
}

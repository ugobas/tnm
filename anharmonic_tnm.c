#include "coord.h"
#include "nma_para.h"
#include "tnm.h"
#include "buildup.h"
#include "McLachlan.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allocate.h"
#include "diagonalize.h"
#include "choldc.h"
#include "simulation.h"
#include "interactions_tnm.h"
#include "Residues_force_constant.h"
//#include "Residues_propensity.h"

int printed_modes=0;
int ENEG=0;

// Global parameters
int INI_ANHAR=0;
int INI_ENERGY=0;
int **nat_res1_res2;
double lambda_4, lambda_2;
int NPAR=2; // Number of parameters: 2 (only repulsion), 3 (barrier)
float C22;
double Ene_anhar_ini;
struct interaction *Int_list_ini; int N_int_ini;
struct interaction *Int_list_tmp; int N_int_tmp;
float *coord_ini;
// Covalent energy
int INI_COV=0;
double E_COV_FACT=2; // Factor between bonded and not bonded energy
double Ene_cov_ini;
struct interaction *Int_list_cov; int N_cov;

static int INT_MAX;        // Maximum number of interactions tested
float d_CA_high=22; // Maximum distance in stored interactions
float MAX_ANGLE=0.05; // Maximum allowed angular deviation
float RMSD_STEP=1.0;  // Interaction list computed when RMSD > RMSD_STEP

// Scale of energies in units of V''
float Scale_Cont_POW=282.3663;
float Scale_Prop_POW=120.8132;
float Scale_Cont_EXP=32.3058;
float Scale_Prop_EXP=13.8224;

// Fit
//float dE_max=15;    // Maximum value of (omega*c)^2
//float dE_step=0.01; // Step value of (omega*c)^2
double **Corr_inv;
double **pow_ia;

// Routines:

double Energy_anharmonic(float *r, atom *atoms, int natoms, int nres,
			 struct interaction *Int_list, int N_int,
			 struct interaction **Int_KB, int NA,
			 struct residue *seq,
			 struct axe *axe, int naxes,
			 double *delta_phi);
double Energy_debug(float *r, atom *atoms, int natoms, int nres,
		    struct interaction *Int_list, int N_int,
		    struct interaction **Int_KB, int NA,
		    struct residue *seq);
double Energy_covalent(float *r, atom *atoms, int natoms,
		       struct bond *bonds);

float Compute_V(struct interaction *Int, float r);
float Compute_V1(struct interaction *Int, float r);
float Compute_V2(struct interaction *Int, float r, int nomin);

int Determine_Ak(struct interaction *Int_ptr, double V0, double K, int n);
void Interpolation(struct interaction *Int, struct interaction *Int2,
		   float r1, float r2);
void Print_potential(FILE *file_out, struct interaction *Int,
		     int a, int b, char *what, char *AA_code);


int **Set_interactions(struct interaction *Int_list, int N_int,
		       atom *atoms, int Nres);
int Interactions_CA(struct interaction *Int_list, float thr,
		    atom *atoms, int natoms, int N_res, float *coord);
int Find_CA(atom *atoms, int *i1, int res, int natoms);
int Start_res(atom *atoms, int i, int res);
float Min_res_dist(float *r, atom *atoms, int *i1, int *i2, int n, int natm);
int Dist_thr(struct interaction *Int_cov, atom *atoms);

// Fit
double **Compute_corr_inv(int npara, int nsam);
double **Compute_pow_ia(int npara, int nsam);
void Get_fit_para(double *para, double *y,
		  double **Corr_inv, int npara, int nsam, float scale);
float Polynomial(float x, double *para, int n);

//////////////////////////////////////////////////////////////////////////
float Anharmonicity(struct Normal_Mode NM, int ia,
		    atom *atoms, int natoms,
		    struct axe *axe, int naxes,
		    struct bond *bonds,
		    struct residue *seq, int nres,
		    struct interaction *Int_list, int N_int,
		    struct interaction **Int_KB, int NA,
		    char *nameout)
{
  // Fit of energy
  /* The normal mode coordinate c is moved at steps of c=x/omega
     equal to x_step/omega, until x<=x_max, i.e. E_har<=0.5*x_max^2,
     //or E_anhar <= E_thr
     NSAM-1 energy differences are fitted with NPARA-1 parameters
     and each interval is interpolated with it_max internal points
     If the interpolated energy is negative, the energy is computed exactly
   */
  if((NM.select[ia]==0)||(NM.omega[ia]==0))return(-1);

  int C_LIN=1; // Compute conformations with linear approximation?

  float Pmin=0.00015; // Stop when the Boltzmann factor reaches Pmin
  int step_min=100;   // Minimum number of steps in one direction
  // Determine maximum possible factor for mode ia
  float MAX_T=0.1; // Maximum torsion angle in radiants
  float tmax=0, *Tors=NM.Tors[ia]; int i;
  for(i=0; i<naxes; i++)if(fabs(Tors[i])>tmax)tmax=fabs(Tors[i]);
  double omega=NM.omega[ia], omega2=omega*omega, omega22=omega2/2;
  double c2_max=-2.1*log(Pmin)/omega2, c_max=sqrt(c2_max), c_step=MAX_T/tmax; 
  int nstep=c_max/c_step;
  if(nstep<step_min){c_step=c_max/step_min; nstep=step_min;}
  int PDB_step=1000; if(N_PDB_print)PDB_step=nstep/N_PDB_print;
  double c2_all=0, c_all=0, dc=fabs(c_step);
  nstep=0;

  int MAKE_FIT=0;  // Interpolate points through fit?
  int NPARA=4, NSAM=6, it_max=10;
  double c0=0, c_step_fit=c_step/it_max, dcf=fabs(c_step_fit);
  // Prepare fits
  double para_Ene[NPARA], para_RMSD[NPARA];
  double E_fit[NSAM], RMSD_fit[NSAM];
  para_Ene[0]=0; para_RMSD[0]=0; E_fit[0]=0; RMSD_fit[0]=0;
  int i_sam=1, n_max=(NSAM-1)*it_max;

  int N_MODE_PRINT=10;
  int PRINT_MODE=(printed_modes < N_MODE_PRINT);
  FILE *file_out=NULL; char out_all[200];
  sprintf(out_all, "%s.anharmonic.dat", nameout);

  // Angular displacements
  float d_phi[naxes]; double delta_phi[naxes];
  for(i=0; i<naxes; i++){delta_phi[i]=0; d_phi[i]=c_step*Tors[i];}

  // Initialization of coordinates
  Set_bonds_measure(bonds, natoms, atoms);
  int n3=3*natoms;

  // Initialization of global parameters
  if(INI_ANHAR==0){
    INI_ANHAR=1;
    printed_modes=0;

    // Initial coordinates
    coord_ini=malloc(n3*sizeof(float));
    int num_atom=Put_coord(coord_ini, bonds, natoms);
    if(num_atom!=natoms){
      printf("ERROR, different number of atoms in bonds (%d) and atoms (%d)\n",
	     num_atom, natoms); exit(8);
    }
    // Initial energy
    Ene_anhar_ini=Energy_anharmonic(coord_ini,atoms,natoms,
				    nres,Int_list,N_int,Int_KB,NA,seq,
				    axe, naxes, delta_phi);
    printf("Writing %s\n", out_all);
    file_out=fopen(out_all, "w");
    fprintf(file_out, "#1=Mode 2=Coll 3=1/omega 4=E_harmonic ");
    fprintf(file_out,
	    "5=1/omega_anh 6=E_anh 7=<c_anh>/sqrt(<c_anh^2>) 8=dKL 9=rmsd_anh");
    fprintf(file_out, " 10=t_max 11=c_step 12=nstep 13=RMSD_skew");
    if(C_LIN){
      fprintf(file_out," 14=1/omega_anhar_lin");
      //fprintf(file_out," 15=E_anharm_lin 16=dKL_lin 17=rmsd_lin");
    }
    fprintf(file_out, "\n");
    fclose(file_out);

    // Initialize_fit
    if(MAKE_FIT){
      Corr_inv=Compute_corr_inv(NPARA, NSAM);
      pow_ia=Compute_pow_ia(NPARA, NSAM);
    }
  } 

  // RMSD
  float mass[natoms], coord_new[n3], coord_ref[n3];
  float rmsd_ref=0, rmsd_max=0, rmsd_ini=0;;
  for(i=0; i<natoms; i++)mass[i]=0;
  Set_masses(mass,atoms,natoms);
  for(i=0; i<n3; i++){
    coord_ref[i]=coord_ini[i];
    coord_new[i]=coord_ini[i];
  }

  // Initialize energy
  ENEG=0;
  N_int_tmp=N_int_ini;
  for(i=0; i<N_int_tmp; i++)Int_list_tmp[i]=Int_list_ini[i];
  double Ene_anhar=0, Ene_har=0;

  // Output file
  if(PRINT_MODE){
    printed_modes++;
    char out[200];
    sprintf(out, "%s.Mode%d_energy.dat", nameout, ia);
    file_out=fopen(out, "w");
    printf("Writing harmonic and anharmonic energy in %s\n", out);
    fprintf(file_out, "#Mode %d, 1/omega= %.3g ", ia, 1./omega);
    fprintf(file_out, "Max.ampl.= %.3g step=%.3g max.angle= %.2g\n",
	    c_max, c_step, c_step*tmax);
    fprintf(file_out, "# %d native contacts E_ini= %.3g\n",
	    N_int, Ene_anhar_ini);
    fprintf(file_out, "#c_alpha Ene_harmonic rmsd(ini,tor) Ene_noharmonic");
    if(C_LIN){
      fprintf(file_out," rmsd(ini,lin) Ene_lin rmsd(tor,lin)");
      fprintf(file_out," Ene_lin+Ene_cov_lin"); //Ene_cov_tors
    }
    fprintf(file_out,"\n");
    //Print_PDB(file_out, atoms, natoms, coord1, seq, 0, 0.0);
    fprintf(file_out, "%.4g\t%.2g\t%.4g\t%.4g",
	    c_all, rmsd_ini, Ene_har, Ene_anhar_ini-Ene_anhar_ini);
    if(C_LIN)fprintf(file_out, "\t%.3g\t%.4g\t%.3g\t%.3g",
		     0.0, 0.0, 0.0, 0.0); //\t%.3g, 0.0
    fprintf(file_out,"\n");
  }else{
    file_out=NULL;
  }

  FILE *file_pdb=NULL; int i_step=0, pdb=0;
  if(PRINT_PDB_ANHARM && PRINT_MODE){
    char name_pdb[100]; sprintf(name_pdb, "%s.mode%d.pdb", nameout, ia);
    file_pdb=fopen(name_pdb, "w");
    printf("Printing %d PDB files in %s\n", 2*N_PDB_print, name_pdb);
    fprintf(file_pdb, "REMARK RMSD=0 C_MODE=0\n");
    Print_PDB(file_pdb, atoms, natoms, coord_ini, seq, 0, 0.0);
  }

  // Counters
  int turn=0;
  double Z, Z0, E_ave_har=0, E_ave_anhar=0, Ene_max=0;
  double RMSD_skew=0; 
  double c1_ave_anhar=0, c2_ave_anhar=0, msd_ave_tors=0;
  if(MAKE_FIT){Z0=Z=c_step_fit;}
  else{Z0=Z=c_step;}

  // Linear approximation (if C_LIN)
  double Z_lin=Z, E_ave_lin=0, c1_ave_lin=0, c2_ave_lin=0, msd_ave_lin=0;
  double rmsd_lin, Ene_lin, rmsd2, Ene_cov_lin; //, Ene_cov_tors
  float coord_lin[n3], *Cart=NM.Cart[ia];

  while(1){

    // Compute coordinates, RMSD and energy and print
    nstep++; c_all+=c_step; c2_all=c_all*c_all;
    Ene_har=c2_all*omega22;

    for(i=0; i<naxes; i++)delta_phi[i]=c_all*Tors[i];

    Build_up(bonds, natoms, d_phi, naxes);
    Put_coord(coord_new, bonds, natoms);
    rmsd_ini= rmsd_mclachlan_f(coord_ini, coord_new, mass, natoms);
    rmsd_ref= rmsd_mclachlan_f(coord_ref, coord_new, mass, natoms);
    if(rmsd_ref > RMSD_STEP){
      for(i=0; i<n3; i++)coord_ref[i]=coord_new[i];
      N_int_tmp= 
	Interactions_CA(Int_list_tmp,d_CA_high,atoms,natoms,nres,coord_ref);
    }

    if(file_pdb){
      i_step++; if(i_step<PDB_step)continue;
      fprintf(file_pdb, "REMARK RMSD= %.2f C_MODE= %.4g\n",rmsd_ini, c_all);
      pdb++; Print_PDB(file_pdb, atoms, natoms, coord_new, seq, pdb, rmsd_ini);
      i_step=0;
    }

    Ene_anhar=Energy_anharmonic(coord_new,atoms,natoms,
				nres,Int_list,N_int,Int_KB,NA,seq,
				axe, naxes, delta_phi)
      - Ene_anhar_ini;
    //if(C_LIN)Ene_cov_tors=Energy_covalent(coord_new,atoms,natoms,bonds);
    if((Ene_anhar < -0.1)&&(ENEG==0)){
      printf("ERROR, mode %d has neg. energy",ia);
      printf(" E= %.4g E_har= %.2g c= %.2f rmsd= %.2g\n",
	     Ene_anhar, Ene_har, c_all, rmsd_ini);
      printf("Native energy: ");
      Energy_debug(coord_ini,atoms,natoms,nres,Int_list,N_int,Int_KB,NA,seq);
      printf("Actual energy: ");
      Energy_debug(coord_new,atoms,natoms,nres,Int_list,N_int,Int_KB,NA,seq);
      ENEG=1;
      //exit(8);
    }

    if(C_LIN){
      for(i=0; i<n3; i++)coord_lin[i]=coord_ini[i]+c_all*Cart[i];
      rmsd2   = rmsd_mclachlan_f(coord_new, coord_lin, mass, natoms);
      rmsd_lin= rmsd_mclachlan_f(coord_ini, coord_lin, mass, natoms);
      Ene_lin=
	Energy_anharmonic(coord_lin,atoms,natoms,nres,
			  Int_list,N_int,Int_KB,NA,seq,axe,naxes,delta_phi)
	- Ene_anhar_ini;
      Ene_cov_lin=Energy_covalent(coord_lin,atoms,natoms,bonds);
      if((Ene_lin < -0.02)&&(ENEG==0)){
	printf("ERROR, mode %d has neg. energy",ia);
	printf(" E= %.4g E_har= %.2g c= %.2f rmsd= %.2g\n",
	       Ene_lin, Ene_har, c_all, rmsd_lin);
	printf("Native energy: ");
	Energy_debug(coord_ini,atoms,natoms,nres,Int_list,N_int,Int_KB,NA,seq);
	printf("Actual energy: ");
	Energy_debug(coord_lin,atoms,natoms,nres,Int_list,N_int,Int_KB,NA,seq);
	ENEG=1;
      }
    }

    if(PRINT_MODE){
      fprintf(file_out, "%.4g\t%.4g\t%.3g\t%.4g",
	      c_all, Ene_har, rmsd_ini, Ene_anhar);
      if(C_LIN){
	fprintf(file_out, "\t%.3g\t%.4g\t%.3g",rmsd_lin,Ene_lin,rmsd2);
	fprintf(file_out, "\t%.4g",Ene_lin+Ene_cov_lin); //, Ene_cov_tors
      }
      fprintf(file_out, "\n");
    }

    // Store for fit
    if(MAKE_FIT==0){
      // Compute integral
      double p0=exp(-Ene_har)*dc, p=exp(-Ene_anhar)*dc;
      if(p0<Pmin)turn=1;
      E_ave_har+=p0*Ene_har; E_ave_anhar+=p0*Ene_anhar;
      Z+=p; Z0+=p0; c2_ave_anhar+=p*c2_all; c1_ave_anhar+=p*c_all;
      msd_ave_tors+=p*rmsd_ini*rmsd_ini;
      RMSD_skew+=rmsd_ini/c_all;
      if(Ene_anhar > Ene_max)Ene_max=Ene_anhar;
      if(C_LIN){
	Ene_lin+=Ene_cov_lin;
	double pl=exp(-Ene_lin)*dc;
	Z_lin+=pl; E_ave_lin+=p0*Ene_lin;
	c1_ave_lin+=pl*c_all; c2_ave_lin+=pl*c2_all;
	msd_ave_lin+=pl*rmsd_lin*rmsd_lin;
      }

    }else{
      E_fit[i_sam]=Ene_anhar-E_fit[0];
      RMSD_fit[i_sam]=rmsd_ini-RMSD_fit[0];
      i_sam++;

      if(i_sam==NSAM){
	// Make fits
	float scale=dc;
	Get_fit_para(para_Ene, E_fit, Corr_inv, NPARA, NSAM, scale);
	Get_fit_para(para_RMSD,RMSD_fit, Corr_inv, NPARA, NSAM, scale);
	if(PRINT_MODE){
	  fprintf(file_out, "# para_Ene:  %.2g %.2g %.2g %.2g scale= %.2g\n",
		  para_Ene[1], para_Ene[2], para_Ene[3], para_Ene[4], scale);
	  fprintf(file_out, "# para_RMSD: %.2g %.2g %.2g %.2g scale= %.2g\n",
		  para_RMSD[1],para_RMSD[2],para_RMSD[3],para_RMSD[4],scale);
	}

	double x=0, c, E, d, Eh;
	for(int it=1; it<=n_max; it++){
	  x+=c_step_fit; c=c0+x;
	  Eh=omega22*c*c;
	  d=RMSD_fit[0]+Polynomial(x, para_RMSD, NPARA);
	  E=E_fit[0]+Polynomial(x, para_Ene, NPARA);
	  if(E<0){
	    // Compute the energy without interpolation
	    printf("Computing E without interpol. mode %d it=%d\n",ia,it);
	    struct bond bonds_tmp[natoms];
	    Copy_bonds(bonds_tmp, bonds, natoms);
	    Set_bonds_measure(bonds_tmp, natoms, atoms);
	    float d_phi_tmp[naxes];
	    double delta_phi_tmp[naxes];
	    for(i=0; i<naxes; i++){
	      delta_phi_tmp[i]=Tors[i]*c; d_phi_tmp[i]=delta_phi_tmp[i];
	    }
	    Build_up(bonds_tmp, natoms, d_phi_tmp, naxes);
	    Put_coord(coord_new, bonds_tmp, natoms);
	    float d_old=d, E_old= E;
	    d=rmsd_mclachlan_f(coord_ini, coord_new, mass, natoms);
	    E=Energy_anharmonic(coord_new,atoms,natoms,nres,Int_list,N_int,
			      Int_KB,NA,seq,axe,naxes,delta_phi_tmp)
	      - Ene_anhar_ini;
	    if(fabs(d_old-d)>0.15){
	      printf("ERROR anhar mode %d c= %.2f interpol. d= %.2f d=%.2f\n",
		     ia, c, d_old, d);
	      FILE *file_tmp=fopen(out_all, "a");
	      fprintf(file_tmp,
		      "#WRONG interpolation, mode %d 1/omega= %.2g c= %.2f",
		      ia, 1./omega, c);
	      fprintf(file_tmp," interpolated d= %.2f d=%.2f\n", d_old, d);
	      fprintf(file_tmp," E_int= %.2g E= %.2g", E_old, E);
	      fclose(file_tmp);
	    }
	  } // end of the energy computation

	  if(PRINT_MODE){
	    fprintf(file_out, "%.4g\t%.4g\t%.3g\t%.4g\n", c, Eh, d, E);
	  }
	  // Compute integral
	  double p0=exp(-Eh)*dcf, p=exp(-E)*dcf;
	  if(p0<Pmin)turn=1;
	  Z+=p; Z0+=p0; E_ave_har+=p0*Eh; E_ave_anhar+=p0*E;
	  c2_ave_anhar+=p*c*c; c1_ave_anhar+=p*c;
	  msd_ave_tors+=p*rmsd_ini*rmsd_ini;
	  if(E > Ene_max)Ene_max=E;
	} // End interpolation loop
	if(PRINT_MODE && MAKE_FIT){
	  fprintf(file_out,"# Fits printed c_step= %.2g c0= %.2g\n",dc,c0);
	}
	i_sam=1;
	E_fit[0]=Ene_anhar;
	RMSD_fit[0]=rmsd_ini;
	c0=c_all;
      }// end interpolation check
    } // end if(MAKE_FIT)

    // Next step
    if(turn ||(c2_all > c2_max)){ //||(E>E_thr)
      if(c_all<0)break;
      turn=0;
      // c>0, change direction
      c_all=0; c0=0; for(i=0; i<n3; i++)coord_ref[i]=coord_ini[i];
      c_step=-c_step; c_step_fit=-c_step_fit; i_step=0;
      for(i=0; i<naxes; i++){delta_phi[i]=0; d_phi[i]=c_step*Tors[i];}
      Set_bonds_measure(bonds, natoms, atoms);
      N_int_tmp=N_int_ini;
      for(i=0; i<N_int_tmp; i++)Int_list_tmp[i]=Int_list_ini[i];
      rmsd_max=rmsd_ini;
      E_fit[0]=0; RMSD_fit[0]=0;
    }

  } // End computations
  if(rmsd_ini>rmsd_max)rmsd_max=rmsd_ini;
  E_ave_har/=Z0; E_ave_anhar/=Z0;
  c2_ave_anhar/=Z; c1_ave_anhar/=Z;
  msd_ave_tors/=Z;
  double c2_ave_har=E_ave_har/omega22;
  float d_KL=E_ave_anhar-E_ave_har+log(Z)-log(Z0);
  float d_KL_lin=0;
  if(C_LIN){
    E_ave_lin/=Z0; 
    d_KL_lin=E_ave_lin-E_ave_har+log(Z_lin)-log(Z0);
  }

  NM.sigma2_anhar[ia]=NM.sigma2[ia]*(c2_ave_anhar/c2_ave_har);
  NM.d_KL[ia]=d_KL;
  if(PRINT_MODE){
    fprintf(file_out, "# <E_har>= %.4g analytical= 0.500 ", E_ave_har);
    fprintf(file_out, "<E_anhar>= %.4g dKL= %.4g", E_ave_anhar, d_KL);
    if(C_LIN)
      fprintf(file_out,
	      " <E_lin_appr>= %.4g dKL_lin= %.4g", E_ave_lin, d_KL_lin);
    fprintf(file_out, "\n");
    fclose(file_out);
  }
  // Print mode results
  file_out=fopen(out_all, "a");
  fprintf(file_out, "%d\t%.3g\t%.4g\t%.3g",
	  ia, NM.Cart_coll[ia], 1./omega, E_ave_har);
  float ca=sqrt(c2_ave_anhar);
  fprintf(file_out, "\t%.4g\t%.3g\t%.2g\t%.2g\t%.3g",
	  sqrt(NM.sigma2_anhar[ia]), E_ave_anhar, 
	  c1_ave_anhar/ca, d_KL, sqrt(msd_ave_tors));
  fprintf(file_out, "\t%.2g\t%.2g\t%d", dc*tmax, dc, nstep);
  fprintf(file_out, "\t%.4f", (RMSD_skew*2/nstep)*ca);
  if(C_LIN){
    c2_ave_lin/=Z_lin;
    float f2=NM.sigma2[ia]*(c2_ave_lin/c2_ave_har);
    fprintf(file_out, "\t%.4g", sqrt(f2));
    /*c1_ave_lin/=Z_lin; msd_ave_lin/=Z_lin;
    fprintf(file_out, "\t%.3g\t%.3g\t%.2g\t%.3g",
	    E_ave_lin, c1_ave_lin/sqrt(c2_ave_lin),
	    d_KL_lin, sqrt(msd_ave_lin));*/
  }
  fprintf(file_out, "\n");
  fclose(file_out);
  if(file_pdb)fclose(file_pdb);

  return(d_KL);
}

//////////////////////////////////////////////////////////////////////////
// returns a matrix with the index of the interaction data (in Int_list) from the two residues ID, or -1 if no interaction
int **Set_interactions(struct interaction *Int_list, int N_int,
		       atom *atoms, int nres)
{
  int **matrix=malloc(nres*sizeof (int *)); int i1, i2;
  for(i1=0; i1<nres; i1++){
    matrix[i1]=malloc(nres*sizeof(int));
    for(i2=0; i2<nres; i2++)matrix[i1][i2]=-1;
  }
  struct interaction *Int=Int_list;
  for(int n=0; n<N_int; n++){
    i1=atoms[Int->i1].res;
    i2=atoms[Int->i2].res;
    if((i1<0)||(i1>=nres)||(i2<0)||(i2>=nres)){
      printf("ERROR in Set_interactions i1=%d i2= %d nres=%d n= %d\n",
	     i1, i2, nres, n); exit(8);
    }
    matrix[i1][i2]=n;
    matrix[i2][i1]=n;
    Int++;
  }
  printf("\n%d native interactions stored\n", N_int);
  return(matrix);
}

//////////////////////////////////////////////////////////////////////////
int Interactions_CA(struct interaction *Int_list, float thr,
		    atom *atoms, int natoms, int N_res, float *coord)
{
  int N_int=0, i1=0;
  float thr2=thr*thr, d, d2;

  for(int res1=0; res1<N_res; res1++){

    if(Find_CA(atoms, &i1, res1, natoms)< 0)continue;
    float *r1=coord+3*i1; atom *atom1=atoms+i1;
    int i2=i1+1;

    for(int res2=res1+1; res2<N_res; res2++){
      if(Find_CA(atoms, &i2, res2, natoms)< 0)continue;
      float *r2=coord+3*i2; atom *atom2=atoms+i2;

      int n=nat_res1_res2[atom1->res][atom2->res];

      // Interatomic distance
      d=r2[0]-r1[0]; d2 =d*d; if((d2>thr2)&&(n<0))continue;
      d=r2[1]-r1[1]; d2+=d*d; if((d2>thr2)&&(n<0))continue;
      d=r2[2]-r1[2]; d2+=d*d; if((d2>thr2)&&(n<0))continue;

      // store interaction
      Int_list[N_int].i1=i1; Int_list[N_int].i2=i2; N_int++;
      if(N_int >= INT_MAX){
	printf("ERROR, too many interactions %d\n", N_int);
	printf("%d residues, Threshold= %.2f\n", N_res, thr);
	exit(8);
      }
    }
  }
  return(N_int);
}

//////////////////////////////////////////////////////////////////////////
int Find_CA(atom *atoms, int *i1, int res, int natoms){
  atom *atom=atoms+*i1; int i=*i1;
  while(atom->res <= res){
    if((atom->res == res)&&(strncmp(atom->name, "CA", 2)==0)){
      *i1=i; return(0);
    }
    atom++; i++; if(i>=natoms)break;
  }
  if(*i1)(*i1)--;
  if(atoms[*i1].res==res)return(0);
  return(-1);
}

//////////////////////////////////////////////////////////////////////////
float Min_res_dist(float *r, atom *atoms, int *i1, int *i2, int n, int natm)
{
  int res1=atoms[*i1].res, j1=Start_res(atoms, *i1, res1);
  // YY ??? think because it's the ID of CA, so goes -1 to start with N.
  int res2=atoms[*i2].res, j2=Start_res(atoms, *i2, res2);

  float d, d2, d2_min=1000;
  atom *atom1=atoms+j1; float *r1=r+3*j1;
  // YY coord of atom j1 (r beginning of array with all atom coord)
  while(atom1->res==res1){
    atom *atom2=atoms+j2; float *r2=r+3*j2; int jj2=j2;
    while(atom2->res==res2){
      if((NOCOV) && Covalent(res1, res2, atom1, atom2))goto next;
      d=r1[0]-r2[0]; d2 =d*d; if(d2>d2_min)goto next;
      d=r1[1]-r2[1]; d2+=d*d; if(d2>d2_min)goto next;
      d=r1[2]-r2[2]; d2+=d*d; if(d2>d2_min)goto next;
      if(d2<d2_min){d2_min=d2; *i1=j1; *i2=jj2;}
    next:
      r2+=3; atom2++; jj2++; if(jj2>=natm)break;
    }
    r1+=3; atom1++; j1++; if(j1>=natm)break;
  }
  if((d2_min<=C22)||(n>=0))return(sqrt(d2_min));
  return(100);
}

//////////////////////////////////////////////////////////////////////////
int Start_res(atom *atoms, int i, int res){
  int k=i, j=i-1;
  while(j>=0 && atoms[j].res==res){k=j; j--;}
  return(k);
}

//////////////////////////////////////////////////////////////////////////
double Energy_anharmonic(float *r, atom *atoms, int natoms, int nres,
			 struct interaction *Int_list, int N_int,
			 struct interaction **Int_KB, int NA,
			 struct residue *seq,
			 struct axe *axe, int naxes,
			 double *delta_phi)
{
  if(INI_ENERGY==0){
    INI_ENERGY=1;
    E_repulsion=Int_list[0].A4;
    C22=C_THR*C_THR;
    lambda_2=2*lambda;
    lambda_4=4*lambda;

    // Auxiliary matrices for energy computation: interaction list
    INT_MAX=150*nres;
    nat_res1_res2=Set_interactions(Int_list, N_int, atoms, nres);
    // YY matrix with the index of the interaction data (in Int_list)
    Int_list_ini=malloc(INT_MAX*sizeof(struct interaction));
    Int_list_tmp=malloc(INT_MAX*sizeof(struct interaction));
    N_int_ini=
      Interactions_CA(Int_list_ini, d_CA_high, atoms, natoms, nres, r);
    N_int_tmp=N_int_ini;
    for(int i=0; i<N_int_tmp; i++)Int_list_tmp[i]=Int_list_ini[i];
    // YY List of interactions with d(CA-CA)<22
  }

  double E_nat=0, E_nonat=0;
  for(int i=0; i<N_int_tmp; i++){
    int i1=Int_list_tmp[i].i1, i2=Int_list_tmp[i].i2;
    int n=nat_res1_res2[atoms[i1].res][atoms[i2].res];
    float d=Min_res_dist(r, atoms, &i1, &i2, n, natoms);
    if(d>=C_THR)continue;
    if(d<=0){
      printf("ERROR, d(%d,%d)= %.2g n=%d\n", i1, i2, d, n);
      exit(8);
    }
    if(n>=0){
      E_nat+=Compute_V(Int_list+n, d); //-Int_list[n].E_min;
    }else{
      int a=seq[atoms[i1].res].i_aa, b=seq[atoms[i2].res].i_aa;
      E_nonat+=Compute_V(Int_KB[a]+b, d); //-Int_KB[a][b].E_min;
    }
  }

  // Compute_torsional energy
  double E_phi=0, E_psi=0, E_omg=0, E_chi=0, E_BA=0, E_BL=0, E_BT=0;
  struct axe *ax=axe;
  for(int i=0; i<naxes; i++){
    float d2=delta_phi[i]*delta_phi[i];
    if(ax->type=='f'){E_phi+=d2;}
    else if(ax->type=='p'){E_psi+=d2;}
    else if(ax->type=='l'){E_BL+=d2;}
    else if(ax->type=='a'){E_BA+=d2;}
    else if(ax->type=='t'){E_BT+=d2;}
    else if(ax->type=='o'){E_omg+=d2;}
    else if(ax->type=='s'){E_chi+=d2;}
    ax++;
  }
  double E_tors=K_PHI*E_phi;
  if(E_psi)E_tors+=K_PSI*E_psi;
  if(E_omg)E_tors+=K_OMEGA*E_omg;
  if(E_chi)E_tors+=K_CHI*E_chi;
  if(E_BL)E_tors+=K_BL*E_BL;
  if(E_BA)E_tors+=K_BA*E_BA;
  if(E_BT)E_tors+=K_PHI*E_BT;
  return(E_nat+E_nonat+0.5*E_tors);
}

double Energy_debug(float *r, atom *atoms, int natoms, int nres,
		    struct interaction *Int_list, int N_int,
		    struct interaction **Int_KB, int NA,
		    struct residue *seq)
{
  double E_nat=0, E_nonat=0; int n_nat=0, n_not=0;
  double d_min=100, d_nat_ave=0; int i1_min=-1, i2_min=-1;
  for(int i=0; i<N_int_tmp; i++){
    int i1=Int_list_tmp[i].i1, i2=Int_list_tmp[i].i2;
    int n=nat_res1_res2[atoms[i1].res][atoms[i2].res];
    float d=Min_res_dist(r, atoms, &i1, &i2, n, natoms);
    if(d>=C_THR)continue;
    if(d<=0){
      printf("ERROR, d(%d,%d)= %.2g n=%d\n", i1, i2, d, n);
      exit(8);
    }
    if(n>=0){
      float e=Compute_V(Int_list+n, d)-Int_list[n].E_min;
      if(e<0)printf("Neg SB de= %.2g n=%d d= %.3f d0=%.3f\n",
		     e, n, d, Int_list[n].r0);
      E_nat+=e;
      float x=fabs(d-Int_list[n].r0);
      if(x>0.01){
	atom *a1=atoms+i1, *a2=atoms+i2;
	printf("ERROR in interaction %d, %s r%d - %s r%d, d= %.2f\n",
	       n, a1->name, a1->res, a2->name, a2->res, d);
	a1=atoms+Int_list[n].i1; a2=atoms+Int_list[n].i2; 
	printf("Contact matrix: %s r%d - %s r%d, d= %.2f\n",
	       a1->name, a1->res, a2->name, a2->res, d); exit(8);
      }
      d_nat_ave+=x;
      n_nat++;
    }else{
      int a=seq[atoms[i1].res].i_aa, b=seq[atoms[i2].res].i_aa;
      float e=Compute_V(Int_KB[a]+b, d)-Int_KB[a][b].E_min;
      if(e<0){
	printf("Non-native int. d=%.5g\n", d);
	printf("Neg KB energy= %.3g ab=%d-%d d= %.3f\n", e, a, b, d);
	printf("A0= %.2g A1= %.2g A2= %.2g A4= %.2g\n",
	       Int_KB[a][b].A0, Int_KB[a][b].A1,
	       Int_KB[a][b].A2, Int_KB[a][b].A4);
	int j=3*i1;
	printf("a1: %s %d r%d\tr= %.3g %.3g %.3g\t",
	       atoms[i1].name, i1, atoms[i1].res,
	       r[j], r[j+1], r[j+2]);
	//atoms[i1].r[0], atoms[i1].r[1], atoms[i1].r[2]);
	j=3*i2;
	printf("a2: %s %d r%d\tr= %.3g %.3g %.3g\n",
	       atoms[i2].name, i2, atoms[i2].res,
	       r[j], r[j+1], r[j+2]); 
	//atoms[i2].r[0], atoms[i2].r[1], atoms[i2].r[2]);
      }
      E_nonat+=e; n_not++;
      if(d<d_min){d_min=d; i1_min=i1; i2_min=i2;}
    }
  }
  printf("%d tested interactions over %d pairs, ", N_int_tmp, nres*(nres-1)/2);
  printf("%d native over %d, %d non-native d_thr= %.2f\n",
	 n_nat, N_int, n_not, C_THR);
  printf("Energy: native %.4g |r0-d|=%.3g non-native %.4g\n",
	 E_nat, d_nat_ave/n_nat, E_nonat);
  if(i1_min>=0){
    printf("d_min_nonat=%.5g ", d_min);
    printf("a1: %s %d r %d  ",atoms[i1_min].name, i1_min, atoms[i1_min].res);
    printf("a2: %s %d r %d\n",atoms[i2_min].name, i2_min, atoms[i2_min].res);
  }
  return(E_nat+E_nonat);
}

double Energy_covalent(float *r, atom *atoms, int natoms,
		       struct bond *bonds)
{
  int i, j;
  if(INI_COV==0){
    INI_COV=1;
    int nmax=2*natoms;
    Int_list_cov=malloc(nmax*sizeof(struct interaction));
    struct interaction *Int_cov=Int_list_cov; int n=0;
    struct bond *bond=bonds;
    for(i=0; i<natoms; i++){
      if(bond->previous && bond->len<0){
	Int_cov->i1=bond->i_atom; Int_cov->i2=bond->previous->i_atom;
	if(Dist_thr(Int_cov, atoms)){
	  Int_cov->sec_der=E_COV_FACT; Int_cov++; n++;
	}
	if(bond->previous->previous && bond->angle<0){
	  Int_cov->i1=bond->i_atom;
	  Int_cov->i2=bond->previous->previous->i_atom;
	  if(Dist_thr(Int_cov, atoms)){
	    Int_cov->sec_der=E_COV_FACT; Int_cov++; n++;
	  }
	}
	if(n>nmax){printf("ERROR too many covalent bonds\n"); exit(8);}
      }
      bond++;
    }
    printf("%d covalent bonds stored for %d atoms\n", n, natoms);
    N_cov=n;
    Compute_sec_der(Int_list_cov, N_cov, NULL);
    double E_rep_tmp=E_repulsion; E_repulsion*=E_COV_FACT;
    Assign_interactions_SB(Int_list_cov, N_cov, NULL, NULL);
    E_repulsion=E_rep_tmp;

    Ene_cov_ini=0; Int_cov=Int_list_cov;
    for(i=0; i<N_cov; i++){
      Ene_cov_ini+=Compute_V(Int_cov, Int_cov->r0);
      Int_cov++;
    }
  } // End initialization of the energy

  double Ene_cov=0; struct interaction *Int_cov=Int_list_cov;
  for(i=0; i<N_cov; i++){
    double d2=0, d;
    float *r1=r+3*Int_cov->i1, *r2=r+3*Int_cov->i2;
    for(j=0; j<3; j++){d=r1[j]-r2[j]; d*=d; d2+=d;}
    Ene_cov+=Compute_V(Int_cov, sqrt(d2));
    Int_cov++;
  }
  return(Ene_cov-Ene_cov_ini);
}

int Dist_thr(struct interaction *Int_cov, atom *atoms){
  if(Int_cov->i1>=0 && Int_cov->i2>=0){   
    double *r1=atoms[Int_cov->i1].r, *r2=atoms[Int_cov->i2].r;
    double d2=0, d; for(int j=0; j<3; j++){d=r1[j]-r2[j]; d*=d; d2+=d;}
    Int_cov->r0=sqrt(d2);
    if(Int_cov->r0 < C_THR){return(1);}
    else{return(0);}
  }
  return(0);
}

//////////////////////////////////////////////////////////////////////////
void Assign_interactions_KB(struct interaction **Int_KB, int Na,
			    int POW, float EXPO, float k0)
{
  // Only repulsion
  printf("Assigning knowledge-based interactions\n");
  printf("--> NPAR = %8d        number of parameters: ", NPAR);
  printf("2 (only repulsion), 3 (barrier)\n");
  printf("--> Na   = %8d        amino acids or atom types\n",Na);
  printf("--> POW  = %8d        1 power law, 0 exponential (kappa)\n", POW);
  printf("--> EXPO = %8.2f        kappa exponent\n",EXPO);
  printf("--> k0   = %8.2f        kappa0 (default value or Bfactor fit)\n",k0);
  printf("--> OUTPUT: Int_KB[%d][%d]\n",Na,Na);
  printf("\n");

  double A4, Scale;
  if(POW){
    lambda=EXPO/4-0.5;
    A4=k0*pow(C0,EXPO); ///(8*lambda*lambda);
    Scale=Scale_Prop_POW/(lambda*lambda);
  }else{
    lambda=EXPO/4;
    A4=k0*exp(C0*EXPO); ///(4*lambda*lambda);
    Scale=Scale_Prop_EXP/(lambda*lambda);
  }
  A4/=(8*lambda*lambda);
  lambda_2=2*lambda;
  lambda_4=4*lambda;
  E_repulsion=A4;
  double A2=2*A4;

  int a, b, n=0, n_err=0; 
  for(a=0; a<Na; a++){
    for(b=a; b<Na; b++){
      double fc2, r3;
      struct interaction *Int_ptr=Int_KB[a]+b;
      Int_ptr->r1=100;
      if(NPAR==2){
	Int_ptr->r0=C_THR;
	Int_ptr->E_min=0;
	Int_ptr->A4=A4;
	if(POW){
	  fc2=pow(C_THR, -lambda_2);
	}else{
	  fc2=exp(-lambda_2*C_THR);
	}
	Int_ptr->A2=-A2*fc2;
	Int_ptr->A0=-Int_ptr->A4*fc2*fc2-Int_ptr->A2*fc2;
	Int_ptr->A1=0;
	r3=C1_THR;
      }else{ // Three parameters fit
	if(D_Res[a][b]>=C_THR){
	  printf("ERROR, too large KB distance %.2f (> %.2f)",
		 D_Res[a][b], C_THR);
	  printf("for residues %d and %d\n", a, b); exit(8);
	}
	Int_ptr->r0=D_Res[a][b];
	double V0=U_Res[a][b]*Scale;
	if(Determine_Ak(Int_ptr, V0, K_Residues[a][b], n))n_err++;
	r3=(Int_ptr->r0+C_THR)*0.5;
	if(r3>C1_THR)r3=C1_THR;
	Int_ptr->E_min=Compute_V(Int_ptr, Int_ptr->r0);
      }
      Int_ptr->r1=r3;
      Int_ptr->r2=r3;
      Int_ptr->sec_der=Compute_V2(Int_ptr, Int_ptr->r0, 0);
      Int_ptr->a0=0; Int_ptr->a1=0; Int_ptr->a2=0;
      Int_ptr->a3=0; Int_ptr->a4=0; Int_ptr->a5=0;
      //Interpolation(Int_ptr, NULL, r3, C_THR);
      if(b!=a)Int_KB[b][a]=*Int_ptr;
      n++;
    }
  }
  if(n_err){
    printf("ERRORS in KB in %d pairs of residues/atoms out of %d\n",
	   n_err, n); exit(8);
  }
}

//////////////////////////////////////////////////////////////////////////
int Determine_Ak(struct interaction *Int, double V0, double KK, int n)
{
  int err=0;
  double f, fc, V2=KK/(lambda*lambda);
  Int->r1=100;
  if(POW){
    f=pow(Int->r0, -lambda);
    fc=pow(C_THR, -lambda);
    V2*=(Int->r0*Int->r0);
  }else{
    f=exp(-lambda*Int->r0);
    fc=exp(-lambda*C_THR);
  }
  double alpha=fc/f, a2=alpha*alpha, aa=1.-alpha;
  double f2=f*f, fc2=fc*fc;
  Int->A4=(V0+0.5*V2*aa*aa)/(f2*f2*(3.0-a2*a2+6*a2-8.0*alpha)); // Condition on V
  if(Int->A4<0){
    printf("ERROR A4= %.3g <0 n= %d", Int->A4, n);
    printf(" alpha=%.2g V0= %.2g V2= %.2g\n", alpha, V0, V2); 
    err=-1;
  }
  Int->A2=-6.0*Int->A4*f2+0.5*V2/f2; // Condition on V''
  Int->A1=8.0*Int->A4*f2*f-V2/f;     // Condition on V'
  Int->A0=-Int->A4*fc2*fc2-Int->A2*fc2-Int->A1*fc; // Condition at C_THR
  { // Test
    float EPS=0.01;
    float V0_1=Compute_V(Int, Int->r0);
    float V1_1=Compute_V1(Int, Int->r0);
    float V2_1=Compute_V2(Int, Int->r0, 0);
    if(err ||(fabs(V0_1-V0)>EPS)||(fabs(V1_1)>EPS)||(fabs(V2_1-KK)>EPS)){
      printf("ERROR in Determine_Ak n= %d r0= %.2f r1= %.2f A4= %.2g\n",
	     n, Int->r0, Int->r1, Int->A4);
      printf("V(Rc)= %.4g\n",  Compute_V(Int, C_THR));
      printf("V(r0)  = %.4g (old) %.4g (new)\n", V0, V0_1);
      printf("V'(r0) = 0 (old) %.4g (new)\n", V1_1);
      printf("V''(r0)= %.4g (old) %.4g (new)\n", KK, V2_1);
      return(-1);
    }
  }
  return(err);

}

//////////////////////////////////////////////////////////////////////////
void Assign_interactions_SB(struct interaction *Int_list, int N_int,
			    struct interaction **Int_KB, int *aa_seq)
{
  printf("Assigning structure-based interactions\n");
  printf("--> NPAR  = %8d        number of parameters: ",NPAR);
  printf("2 (Morse) 3 (barrier)\n");
  printf("--> N_int = %8d        number of native interactions\n",N_int);
  printf("--> Int_list[]         list of all native interactions\n");
  printf("--> Int_KB[][]         knowledge-based interactions\n");
  printf("--> aa_seq[]           a.a. type ID of each residue in the seq\n");
  printf("--> OUTPUT: Int_list gets filled up with energy parameters\n");

  int n, n_err=0; struct interaction *Int=Int_list, *I_KB=NULL;
  for(n=0; n<N_int; n++){
    double V0, r3;
    Int->r1=100;
    if(NPAR==2 || aa_seq==NULL || Int_KB==NULL){
      Int->A1=0;
      Int->A4=E_repulsion;
      double fc;
      if(POW){
	Int->A2=-2*Int->A4*pow(Int->r0, -lambda_2);
	fc=pow(C_THR, -lambda);
      }else{
	Int->A2=-2*Int->A4*exp(-lambda_2*Int->r0);
	fc=exp(-lambda*C_THR);
      }
      double fc2=fc*fc;
      Int->A0=-Int->A4*fc2*fc2-Int->A2*fc2;
      V0=Compute_V(Int, Int->r0);
    }else{
      int a=aa_seq[Int->res1], b=aa_seq[Int->res1];
      I_KB=Int_KB[a]+b;
      V0=Compute_V(I_KB, Int->r0);
      if(Determine_Ak(Int, V0, Int->sec_der, n))n_err++;
      if(Int->A4<0){
	printf("POWER= %d lambda= %.2f (rSB/rKB)= %.3g\n",
	       POW, lambda, Int->r0/I_KB->r0);
      }
    }
    Int->E_min=V0;
    double kk=Compute_V2(Int, Int->r0, 0);
    if(fabs(1-Int->sec_der/kk)>0.001){
      printf("ERROR in Assign_int_SB n=%d r= %.2f ", n, Int->r0);
      printf("k= %.4g (old) %.4g (new) /= %.3f  ", 
	     Int->sec_der, kk, Int->sec_der/kk);
      printf("V0= %.3g\n", Compute_V(Int, Int->r0));
      if(I_KB){
	printf("Knowledge-based: NPAR=%d r0= %.3g r1=%.2g r2=%.2g r3=%.2g\n",
	       NPAR, I_KB->r0, I_KB->r1, I_KB->r2, I_KB->r3);
	printf("b= %.3g %.3g %.3g\n", I_KB->b3, I_KB->b4, I_KB->b5);
	printf("V_KB=  %.3g\n",Compute_V(I_KB, Int->r0));
      }
      exit(8);
    }
    if(NPAR==2){
      // Interpolation is not a good idea
      //r3=(Int->r0+C_THR)*0.5; if(r3>C1_THR)r3=C1_THR;
      r3=C_THR;
      Int->r1=r3; Int->r2=r3;
    }else{
      // Interpolate structure based and knowledge based
      float r1, r2;
      if(Int->r0<I_KB->r0){r1=Int->r0; r2=I_KB->r0;}
      else{r1=I_KB->r0; r2=Int->r0;}
      //r3=(r2+C_THR)*0.5; if(r3>C1_THR)r3=C1_THR;
      r3=C_THR;
      float d=(r2-r1)/3; r1+=d; r2-=d;
      Interpolation(Int, I_KB, r1, r2);
      if(isnan(Int->a0) || isnan(Int->a1) || isnan(Int->a2) ||
	 isnan(Int->a3) || isnan(Int->a4) || isnan(Int->a5)){
	printf("ERROR in Assign_SB, n=%d KB= %d\n", n, NPAR);
	printf("Interpolation r1= %.2f r2= %.2f\n", r1, r2);
	printf("a= %.3g %.3g %.3g %.3g %.3g %.3g\n",
	       Int->a0, Int->a1, Int->a2, Int->a3, Int->a4, Int->a5);
	exit(8);
      }
    }
    Int->r3=r3;
    /*Interpolation(Int, NULL, r3, C_THR);
    if(isnan(Int->b3) || isnan(Int->b4) || isnan(Int->b5)){
      printf("ERROR in Assign_SB, n=%d KB= %d\n", n, NPAR);
      printf("Interpolation r3= %.2f r4= %.2f\n", r3, C_THR);
      printf("b= %.3g %.3g %.3g\n", Int->b3, Int->b4, Int->b5);
      exit(8);
      }*/
    Int++;
  }
  printf("%d SB interactions assigned\n", N_int);
  if(n_err){
    printf("ERRORS in SB in %d pairs of residues/atoms\n", n_err);
    //exit(8);
  }
}

//////////////////////////////////////////////////////////////////////////
void Print_interactions(struct interaction *Int_list, int N_int,
			struct interaction **Int_KB,
			int *aa_seq, char *AA_code, int NA, char *name)
{
  // Prepare for printing
  char nameout[100];
  sprintf(nameout,"%s.Inter_SB_KB.dat",name);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  int printed[4]; for(int i=0; i<4; i++)printed[i]=0;
  fprintf(file_out, "# NPAR= %d\n#r V(r)\n", NPAR);

  int n, nprint=0; struct interaction *Int=Int_list;
  for(n=0; n<N_int; n++){
    int a=aa_seq[Int->res1], b=aa_seq[Int->res1];
    struct interaction *I_KB=Int_KB[a]+b;
    int type;
    if(Int->r0>I_KB->r0){type=0;}else{type=2;}
    if(I_KB->E_min>0)type++;
    if(printed[type]==0){
      Print_potential(file_out, Int, a, b, "Structure-based", AA_code);
      Print_potential(file_out, I_KB, a,b, "Knowledge-based", AA_code);
      printed[type]=1; nprint++;
      if(nprint==4)break;
    }
    Int++;
  }
  fclose(file_out);
}

void Print_interactions_KB(struct interaction **Int_KB,
			   char *AA_code, int NA, char *name)
{
  // Prepare for printing
  char nameout[100];
  sprintf(nameout,"%s.Inter_KB.dat",name);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  fprintf(file_out, "# NPAR= %d\n#r V(r)\n", NPAR);
  for(int a=0; a<NA; a++){
    for(int b=0; a<NA; a++){
      struct interaction *I_KB=Int_KB[a]+b;
      Print_potential(file_out, I_KB, a,b, "Knowledge-based", AA_code);
    }
  }
  fclose(file_out);
}


//////////////////////////////////////////////////////////////////////////
void Print_potential(FILE *file_out, struct interaction *Int,
		     int a, int b, char *what, char *AA_code)
{
  fprintf(file_out, "# %s potential %c-%c, r0= %.2f Vmin= %.3f V''= %.3f\n",
	  what, AA_code[a], AA_code[b], Int->r0, Int->E_min, Int->sec_der);
  fprintf(file_out, "# A4= %.3g A2= %.3g A1= %.3g A0= %.3g\n",
	  Int->A4, Int->A2, Int->A1, Int->A0);
  fprintf(file_out, "# r1= %.2f r2= %.2f r3= %.2f\n", Int->r1, Int->r2, Int->r3);
  int n=100; float r=1.5, dr=(5.0-r)/n;
  for(int k=0; k<n; k++){
    fprintf(file_out, "%.2f %.4g\n", r, Compute_V(Int, r)); r+=dr;
  }
  fprintf(file_out, "&\n");
} 


//////////////////////////////////////////////////////////////////////////
void Assign_interactions_SB_old(struct interaction *Int_list, int N_int,
				int POW, float EXPO, float k0)
{
  float A2; int n;
  if(POW){
    lambda=EXPO/4-0.5; lambda_2=2*lambda;
    A2=k0*pow(C0,EXPO)/(lambda_2*lambda_2);
    for(n=0; n<N_int; n++){
      Int_list[n].A2=A2*pow(Int_list[n].r0,-lambda_2);
    }
  }else{
    lambda=EXPO/4; lambda_2=2*lambda;
    A2=k0*exp(C0*EXPO)/(lambda_2*lambda_2);
    for(int n=0; n<N_int; n++){
      Int_list[n].A2=A2*exp(-Int_list[n].r0*lambda_2);
    }
  }
  float A4=A2/2;
  E_repulsion=A4;
  float kk0=A2*lambda_2*lambda_2, kk, EPS=0.005; int err=0;
  struct interaction *Int=Int_list;
  for(n=0; n<N_int; n++){
    Int->A4=A4;
    Int->A1=0;
    Int->A0=0;
    Int->r1=(Int->r0+C_THR)/2;
    Int->A0=-Compute_V(Int, C_THR);
    Int->E_min=Compute_V(Int, Int->r0);
    // Interpolation(Int_list+n);
    // Check force constant
    if(POW){kk=kk0*pow(Int->r0,-EXPO);}
    else{kk=kk0*exp(-Int->r0*EXPO);}
    if(fabs(kk-Int->sec_der)>EPS){
      err++;
      if(err <= 10){
	printf("ERROR in Assign_int n=%d r= %.2f ", n, Int->r0);
	printf("k= %.4g (old) %.4g (new) /= %.3f\n", 
	       Int->sec_der, kk, Int->sec_der/kk);
      }else{
	printf(".");
      }
    }
    Int++;
  }
  if(err){
    printf("\n%d errors in Assign_int, POW= %d EXPO= %.1f K0= %.4g r0= %.1f\n",
	   err, POW, EXPO, k0, C0); exit(8);
  }

  // Print energy versus r0
  char name_out[200]="Anharmonic_energy_TNM.dat";
  FILE *file_out=fopen(name_out, "w");
  printf("Writing %s\n", name_out);
  if(POW){ // power law: C0^-e
    fprintf(file_out, "# f(r)= %.2f r^(-%.2f)\n", -A4, EXPO);
  }else{ // Exponential
    fprintf(file_out, "# f(r)= %.2f exp(-%.2f*r)\n", -A4, EXPO);
  }
  for(float r=2.0; r<4.6; r+=0.005){
    float f;
    if(POW){f=-A4*pow(r,-EXPO);}
    else{f=-A4*exp(-EXPO*r);}
    fprintf(file_out, "%.2f %.4g\n", r, f);
  }
  fclose(file_out);
}

//////////////////////////////////////////////////////////////////////////
void Interpolation(struct interaction *Int,
		   struct interaction *Int2,
		   float r1, float r2)
{
  // Interpolate between r1 and r2 so that V=0 at r2,
  // V+=V-, V1+=V1-, V2+=V2- at r1
  float Dr=r1-r2, Dr2=Dr*Dr, Dr3=Dr2*Dr;
  float V0, V1, V2;  // Constraints at r1
  struct interaction *I1, *I2;
  if(Int2){
    if(Int2->r0<Int->r0){I1=Int2; I2=Int;}
    else{I1=Int; I2=Int2;}
  }else{
    I1=Int; I2=NULL;
  }
  V0=Compute_V(I1, r1);
  V1=Compute_V1(I1, r1);
  V2=Compute_V2(I1, r1, 1);
  V1*=Dr; V2*=Dr2;

  if(Int2==NULL){
    // Interpolate the final part
    Int->r3=r1;
    Int->b3=(10*V0-4*V1+0.5*V2)/Dr3;
    Int->b4=(-15*V0+7*V1-V2)/(Dr3*Dr);
    Int->b5=(6*V0-3*V1+0.5*V2)/(Dr3*Dr2);
  }else{ // Interpolate the median part
    Int->r1=r1; Int->r2=r2;
    if(Int2->E_min<Int->E_min)
      Int->E_min=Int2->E_min;
    // Constraints at r2
    float V0_2=Compute_V(I2, r2);
    float V1_2=Compute_V1(I2, r2);
    float V2_2=Compute_V2(I2, r2, 1);
    Int->a0=V0_2;
    Int->a1=V1_2;
    Int->a2=V2_2*0.5;
    V1_2*=Dr; V2_2*=Dr2; 
    V0-=(V0_2+V1_2+0.5*V2_2);
    V1-=(V1_2+V2_2);
    V2-=V2_2;
    Int->a3=(10*V0-4*V1+0.5*V2)/Dr3;
    Int->a4=(-15*V0+7*V1-V2)/(Dr3*Dr);
    Int->a5=(6*V0-3*V1+0.5*V2)/(Dr3*Dr2);

    if(Int2->r0<Int->r0){
      // The KB potential has minimum before the SB potential
      double A0=Int->A0, A1=Int->A1, A2=Int->A2, A4=Int->A4;
      Int->A0=Int2->A0; Int->A1=Int2->A1;
      Int->A2=Int2->A2; Int->A4=Int2->A4;
      Int->B0=A0; Int->B1=A1;
      Int->B2=A2; Int->B4=A4;
    }else{
      Int->B0=Int2->A0; Int->B1=Int2->A1;
      Int->B2=Int2->A2; Int->B4=Int2->A4;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
float Compute_V(struct interaction *Int, float r)
{
  if(Int==NULL)return(0);
  double V=0, f;
  if(r<=Int->r1){
    // Knowledge based zone for SB or structure based if r0^SB < r0^KB
    if(POW){f=pow(r,-lambda);}else{f=exp(-lambda*r);}
    if(Int->A0)V=Int->A0;
    if(Int->A1)V+=Int->A1*f; 
    f*=f; if(Int->A2)V+=Int->A2*f;
    f*=f; if(Int->A4)V+=Int->A4*f;
  }else if(r<=Int->r2){
    double d=(r-Int->r2), d2=d*d, d3=d2*d;
    V=Int->a0+Int->a1*d+Int->a2*d2+Int->a3*d3+Int->a4*d3*d+Int->a5*d3*d2;
  }else if(r<=Int->r3){
    if(POW){f=pow(r,-lambda);}else{f=exp(-lambda*r);}
    if(Int->B0)V=Int->B0;
    if(Int->B1)V+=Int->B1*f; 
    f*=f; if(Int->B2)V+=Int->B2*f;
    f*=f; if(Int->B4)V+=Int->B4*f;
  }else if(r<C_THR){
    // Final interpolation region
    double d=(r-C_THR), d2=d*d, d3=d2*d;
    V=Int->b3*d3+Int->b4*d3*d+Int->b5*d3*d2;
  }
  return(V);
}

//////////////////////////////////////////////////////////////////////////
float Compute_V1(struct interaction *Int, float r)
{
  if(Int==NULL)return(0);
  double V=0, f, f1;
  if(r<=Int->r1){
    if(POW){f=pow(r,-lambda); f1=lambda*f/r;}
    else{f=exp(-lambda*r); f1=lambda*f;}
    if(Int->A1)V+=Int->A1; 
    if(Int->A2)V+=2.0*Int->A2*f;
    if(Int->A4)V+=4.0*Int->A4*f*f*f;
    V*=-f1;
  }else if(r<=Int->r2){
    double d=(r-Int->r2), d2=d*d;
    V=Int->a1+2*Int->a2*d+3*Int->a3*d2+4*Int->a4*d2*d+5*Int->a5*d2*d2;
  }else if(r<=Int->r3){
    if(POW){f=pow(r,-lambda); f1=lambda*f/r;}
    else{f=exp(-lambda*r); f1=lambda*f;}
    if(Int->B1)V+=Int->B1; 
    if(Int->B2)V+=2*Int->B2*f;
    if(Int->B4)V+=4*Int->B4*f*f*f;
    V*=-f1;
  }else if(r<C_THR){
    // Final interpolation region
    double d=(r-C_THR), d2=d*d;
    V=3*Int->b3*d2+4*Int->b4*d2*d+5*Int->b5*d2*d2;
  }
  return(V);
}

//////////////////////////////////////////////////////////////////////////
float Compute_V2(struct interaction *Int, float r, int no_min)
{
  if(Int==NULL)return(0);
  double V=0, f, f1;
  if(r<=Int->r1){
    if(POW){f=pow(r,-lambda); f1=lambda*f/r;}
    else{f=exp(-lambda*r); f1=lambda*f;}
    V=2.0*Int->A2+12.0*Int->A4*f*f;
    V*=(f1*f1);
    if(no_min){
      double f2=f;
      if(POW){f2*=lambda*lambda;}else{f2*=lambda*(lambda+1)/(r*r);}
      V+=(Int->A1+2*Int->A2*f+4*Int->A4*f*f*f)*f2;
    }
  }else if(r<=Int->r2){
    double d=(r-Int->r2), d2=d*d;
    V=2*Int->a2+6*Int->a3*d+12*Int->a4*d2+20*Int->a5*d2*d;
  }else if(r<=Int->r3){
    if(POW){f=pow(r,-lambda); f1=lambda*f/r;}
    else{f=exp(-lambda*r); f1=lambda*f;}
    if(Int->B2)V+=2.0*Int->B2;
    if(Int->B4)V+=12.0*Int->B4*f*f;
    V*=(f1*f1);
    if(no_min){
      double f2=f;
      if(POW){f2*=lambda*lambda;}else{f2*=lambda*(lambda+1)/(r*r);}
      V+=(Int->B1+2*Int->B2*f+4*Int->B4*f*f*f)*f2;
    }
  }else if(r< C_THR){
    // Final interpolation region
    double d=(r-C_THR), d2=d*d;
    V=6*Int->b3*d+12*Int->b4*d2+20*Int->b5*d2*d;
  }
  return(V);
}

/****************************************************
                       Fits
*****************************************************/
/* Linear fit: x_i=x0+d*i E_i~E(x0)+sum_a=1^4 (p_a*d^a)*i^a i=1...4
   d can also be < 0
   Logarithmic fit: log(|xi|)=log(|x0|)+d*i
   log(E_i)~log(E(x0))+sum_a=1^4 (p_a*d^a)*i^a
   Correlation matrix:
   C_ab= sum_i=1^4 i^a*i^b = sum_i=1^4 i^(a+b)
   p_a*d^a = sum_b C^(-1)_ab sum_i i^b*y_i = C^(-1) Y
 */

double **Compute_corr_inv(int npara, int nsam)
{
  //C_ab= sum_i=1^4 i^a*i^b = sum_i=1^4 i^(a+b)
  double **Corr=Allocate_mat2_d(npara, npara); int a,b,k;
  for(a=0; a<npara; a++){
    for(b=0; b<=a; b++){
      double C=0;
      for(k=1; k<nsam; k++)C+=pow(k,a+b);
      Corr[a][b]=C; Corr[b][a]=C;
    }
  }
  Corr[0][0]=nsam+1;
  double **Corr_inv=Allocate_mat2_d(npara, npara);
  if(choldc(Corr, Corr, npara)==0){ //
    printf("Inverting correlation through Cholevsky decomposition\n");
    double **L_inv=Cholevsky_inversion(Corr, npara);
    for(k=0; k<npara; k++){
      double *v=L_inv[k];
      for(a=0; a<npara; a++)for(b=0; b<=a; b++)Corr_inv[a][b]+=v[a]*v[b];
    }
    Empty_matrix_d(L_inv, npara);
  }else{
    printf("Inverting correlation through diagonalization\n");
    float **eigen_vector=Allocate_mat2_f(npara, npara), eigen_value[npara];
    d_Diagonalize(npara, Corr, eigen_value, eigen_vector, 1);
    int num=0;
    for(k=0; k<npara; k++){
      if(eigen_value[k]<=0){num++; continue;}
      double t=1./eigen_value[k];
      float *v=eigen_vector[k];
      for(a=0; a<npara; a++){
	for(b=0; b<=a; b++)Corr_inv[a][b]+=t*v[a]*v[b];
      }
    }
    if(num)printf("WARNING, correlation has %d null modes\n", num);
    Empty_matrix_f(eigen_vector, npara);
  }
  for(a=0; a<npara; a++)for(b=0; b<a; b++)Corr_inv[b][a]=Corr_inv[a][b];
  Empty_matrix_d(Corr, npara);
  return(Corr_inv);
}

double **Compute_pow_ia(int npara, int nsam){
  double **pow_ia=Allocate_mat2_d(nsam, npara);
  for(int a=0; a<npara; a++){
    for(int i=1; i<nsam; i++)pow_ia[i][a]=pow(i,a);
  }
  return(pow_ia);
}

void Get_fit_para(double *para, double *y, double **Corr_inv,
		  int npara, int nsam, float scale)
{
  // p_a*d^a = sum_b C^(-1)_ab sum_i i^b*y_i = C^(-1) Y
  double YX[npara]; int a, b, i;
  for(a=0; a<npara; a++){
    YX[a]=0; for(i=1; i<nsam; i++)YX[a]+=pow_ia[i][a]*y[i]; //pow(i,a)
  }
  for(a=0; a<npara; a++){
    para[a]=0;
    for(b=0; b<npara; b++)para[a]+=Corr_inv[a][b]*YX[b];
    para[a]*=pow(scale, -a);
  }
  if(0){
    printf("Y for fit: "); for(a=1; a<npara; a++)printf(" %.2g",y[a]);
    printf("\n");
    printf("YX for fit: "); for(a=0; a<npara; a++)printf(" %.2g",YX[a]);
    printf("\n");
    for(a=0; a<npara; a++){if(isnan(YX[a]))exit(8);}
  }
}

float Polynomial(float x, double *para, int n){
  // E_i~E(x0)+sum_a=1^4 (p_a*d^a)*i^a i=1...4
  double y=0;
  for(int a=1; a<n; a++)y+=para[a]*pow(x,a);
  return(y);
}

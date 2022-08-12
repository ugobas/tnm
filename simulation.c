//float STEP_MAX=0.6;
//float STEP_MIN=0.01;
int ITER_MAX=200;
float RMSD_TOL=0.20; // Tolerance RMSD for confchange

//~ #define DEBUG 1
#define TOL 0.15  // Tollerence of RMSD in simulations
#define EXP_FORCE 1.5 // Maximum force constant = E0*E_THR/RMSD(0,1)^EXP_FORCE
#define IDUM -45731    // Seed for simulations

//#include "Parameters.h"
#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "buildup.h"
#include "simulation.h"
#include "atoms.h"
#include "vector.h"
#include "random3.h"
#include "McLachlan.h"
#include "align_tnm.h"
#include "kinetic_tnm.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "choldc.h"
#include "allocate.h"
#include "optimization.h"
#include "externals.h"

float Energy_clashes(float *r, int N);
int Metropolis(float Ene_tmp, float Ene, float T);
void Extract_dphi(float *d_phi, float **Tors_mode, float *sigma,
		  int N_modes, int N_axes, float factor, long *idum);
void Extract_dx(float *dx, float **Cart_mode, float *sigma,
		  int N_modes, int N_Cart, float factor, long *idum);
unsigned long randomgenerator(void);
void Copy_atoms_bonds(atom *atoms, struct bond *bonds, int natoms);
int Decide_direction(struct bond *bonds, float *d_phi, float *atom_str,
		     atom *atoms, int naxes, int natoms);
void Print_PDB(FILE *file_out, atom *atoms, int natoms, float *atom_str3,
	       struct residue *seq, int k, float rmsd);
void Print_seqres(FILE *file_out, struct chain *chains, int Nchain);
int Put_atoms(atom *atoms, float *coord, int natoms);
void Write_all_coord_atom(float *coord, atom *atoms, int natoms);
void Write_ref_coord_atom(float *coord, int N_ref, atom *atoms,
			  int *atom_num);
void Write_ref_coord(float *coord_ref, int N_ref, float *coord_all,
		     int *atom_num);
int Make_conformations(char type, int N_frames, float *d_phi, float *d_Cart,
		       float factor, float coeff, float rmsd0,
		       float *coord_all_f, float *coord_opt_f,
		       float *coord_ref, float *coord_ref2,
		       int PRINT_SIMUL_PDB, FILE *file_out, FILE *file_pdb,
		       float *coord_str1, float *coord_str2,
		       struct Para_simul Para_simul,
		       atom *atoms, int natoms,
		       struct axe *axe, int naxes,
		       struct bond *bonds,
		       struct residue *seq, int nres,
		       int N_diso, struct Normal_Mode NM,
		       struct Jacobian J, struct Reference Ref,
		       struct interaction *Int_list, int N_int,
		       char *INT_TYPE, float s0, char *nameout, int anhar);
static int Check_nan_f(float *x, int N, char *name);
void  Internal_differences(struct bond *bonds, struct bond *bonds2,
			   int natoms, char *prots);
int Test_standard(struct bond *bonds, int natoms);
void Torsional_confchange_1angle_old(struct bond *bonds, int natoms1,
				     int naxes, struct Reference Ref1,
				     float *coord_ref_1, float *coord_ref_2,
				     int N_coord, atom *atoms1,FILE *file_rmsd,
				     struct Para_simul Para_simul);
int Torsional_confchange_1angle(float *diff_tors, int naxes, float angle,
				struct bond *bonds, int natoms1,
				struct ali_atoms ali_a,
				float *coord_ref_1, float *coord_ref_2,
				float STEP_MIN, float RMSD_min, atom *atoms1,
				float *coord_all);
static void Ave_se(double *sum1, double *sum2, int n);


float Cos_sin(float *s, float c);
extern int Periodic_angles(float *dphi, struct axe *axe, int n);
extern void Compare_build_up(struct bond *bonds, int N_atoms,
			     float *d_phi, int N_axes, struct bond *bonds2);

float Compute_RMSD(float f, float *Diff_Tors, int naxes,
		   struct bond *bonds_ini, struct bond *bonds, int natoms1,
		   struct Reference Ref1, float *coord_ref_2);
int Write_RMSD(float *reduce_step, float *RMSD_opt, float *RMSD_kf,
	       float STEP_MIN, float STEP_MAX, float factor,
	       int iter, FILE *file_rmsd,
	       struct bond *bonds_min, float *coord_all, int natoms1,
	       struct ali_atoms ali_a, float *coord_new, float *coord_old,
	       float *coord_ref_1, float *coord_ref_2);
void Copy_ali(struct ali_atoms *ali_a, struct Reference Ref);

int Tors_step(struct Tors *Diff, double **T_Lambda, float **Jacobian_ar,
	      float *mass_coord, int naxes, int N_Cart);

/*static void Initialize(int Ndir, struct axe *axes, int naxes, int N_modes,
		       atom *atoms, int natoms,
		       atom *atoms2, int N_ref,
		       int *atom_num1, int *atom_num2,
		       struct Para_simul Para_simul,
		       struct bond *BONDS);
void Extract_dphi_old(float *d_phi, float **Tors_mode, float *sigma,
		      int N_modes, int N_axes, float SDEV_SIM);
*/

// Required for examine_confchange
struct Tors Diff_s;
struct Ali_score ali_sim;
int Ini_ch=1; 

int INI_MOVE=0;
int Ini_dthr=0;

static long idum; 
int PRINT_ENE=1;
int INIRAN;
//float RANFACTOR;
int INI_ALL;
int INI_ONE_MODE;
/*struct bond *bonds[2];
struct bond *bonds_act=NULL;
float *MASS_ALL=NULL;
float *MASS_REF=NULL;
float *ATOM_ALL1=NULL;
float *ATOM_REF1=NULL;
float *ATOM_REF2=NULL;
float *ATOM_DIFF=NULL;
float *ATOM_PDB=NULL;
float *ATOM_ALL[2];
float *ATOM_REF[2];
float *d_phi[2];
float RMSDif, RMS_NORM;
float Ene_ini, Ene_act, Ene_thr_step, Ene_thr_all, E_thr1, T_all, T_step;
float k_thr;
FILE *file_rmsd=NULL;
int N_modes_active;
int k_accept;
float rmsd_ini_prev=0;
float rmsd_print=0;
char pdbout[200];
char namermsd[200];*/

int Print_mode_PDB(atom *atoms, int natoms,
		   struct axe *axe, int naxes,
		   struct bond *bonds, struct residue *seq,
		   float *Tors, float omega,
		   int N_STEP, struct Para_simul Para_simul,
		   char *nameout, int ia, int is)
{
  // Output file
  char pdbout[200];
  sprintf(pdbout, "%s_mode%d.pdb", nameout, is);
  FILE *file_out=fopen(pdbout, "w");
  printf("Writing normal modes in PDB format in %s\n", pdbout);

  // Coordinates
  Set_bonds_measure(bonds, natoms, atoms);
  int n3=3*natoms, i;
  float coord1[n3], coord2[n3], coord_ini[n3];
  float *coord_old=coord1, *coord_new=coord2;
  float mass[natoms]; Set_masses(mass,atoms,natoms);
  int num_atom=Put_coord(coord_ini, bonds, natoms);
  for(i=0; i<3*num_atom; i++)coord_old[i]=coord_ini[i];
  float Ene=Energy_clashes(coord_ini, num_atom);
  float Ene_thr=100*(Ene+0.1);

  // Determine maximum possible factor for mode ia
  //float factor=Para_simul.AMPLITUDE/(omega*N_STEP);
  float AMAX=16, factor=AMAX/(omega*N_STEP); // Max amplitude = AMAX

  int nstep=N_STEP;
  float MAX_ANGLE=0.15, max_v=0;
  for(i=0; i<naxes; i++)if(fabs(Tors[i])>max_v)max_v=fabs(Tors[i]);
  if(factor*max_v > MAX_ANGLE){
    printf("WARNING, max. d_theta allowed for mode %d = %.3f > %.3f\n",
	   ia, factor*max_v, MAX_ANGLE);
    factor=(MAX_ANGLE/max_v);
    nstep=AMAX/(omega*factor);
  }
  float d_phi[naxes];
  for(i=0; i<naxes; i++)d_phi[i]=Tors[i]*factor;

  Print_PDB(file_out, atoms, natoms, coord1, seq, 0, 0.00);
  Set_bonds_measure(bonds, natoms, atoms);

  int step=0, forward=1, pdb=0, maxstep=0;
  float Ene_max=0, rmsd_max=0;
  while(1){
    if(forward){step++;}else{step--; if(step==0)break;}
    Build_up(bonds, natoms, d_phi, naxes);
    Put_coord(coord_new, bonds, num_atom);
    Ene=Energy_clashes(coord_new, num_atom);
    if(Ene>Ene_max)Ene_max=Ene;
    if((step==nstep)||(forward && (Ene>Ene_thr))){
      // change direction
      forward=0; for(i=0; i<naxes; i++)d_phi[i]=-d_phi[i];
      maxstep=step;
    }

    float rmsd=rmsd_mclachlan_f(coord_old, coord_new, mass, num_atom);
    if((rmsd > Para_simul.PDB_STEP)||((forward==0)&&(pdb==0))){
      pdb++;
      rmsd=rmsd_mclachlan_f(coord_ini, coord_new, mass, num_atom);
      Print_PDB(file_out, atoms, natoms, coord_new, seq, pdb, rmsd);
      float *tmp=coord_old; coord_old=coord_new; coord_new=tmp;
      if(rmsd>rmsd_max)rmsd_max=rmsd;
    }
  }
  fclose(file_out);
  printf("Mode %d, 1/omega= %.3g max.rmsd=%.2f Clash=%.3g snapshots=%d\n",
	 ia, 1./omega, rmsd_max, Ene_max, pdb);
  printf("Max.amplitude= %.3g step=%.3g max.angle= %.3f nstep=%d\n",
	 factor*maxstep, factor, factor*max_v, maxstep);
  return 0;
}


int Metropolis(float Ene_tmp, float Ene, float T)
{
  // Initialize random numbers
  if(INIRAN==0){
     INIRAN=1;
     unsigned long iran=randomgenerator();
     InitRandom( (RANDOMTYPE)iran);
     //RANFACTOR=pow(12.0,1/3.0);
  }
  if(Ene_tmp < Ene)return(1);
  float fact=exp((Ene-Ene_tmp)/T);
  //float ran=RandomFloating()*RANFACTOR;
  //if(ran < fact)return(1);
  if(RandomFloating() < fact)return(1);
  return(0);
}

void Set_masses(float *mass, atom *atoms, int natoms){
  int n=0, i;
  for(i=0; i<natoms; i++){
    //ALIGN
    if(atoms[i].ali==0)continue;
    mass[n] = Mass(atoms+i); n++;
  }
}

void Extract_dphi(float *d_phi, float **Tors_mode, float *sigma,
		  int N_modes, int N_axes, float factor, long *idum)
{
  float gasdev(long *idum);
  int i, k; float z, coeff;

  for(i=0; i<N_axes; i++)d_phi[i]=0;
  for(k=0; k<N_modes; k++){
    z=gasdev(idum);
    coeff = z*sigma[k];
    float *phi=d_phi, *mode=Tors_mode[k];
    for(i=0; i<N_axes; i++){*phi += coeff*(*mode); phi++; mode++;}
  }
  for(i=0; i<N_axes; i++)d_phi[i]*=factor;
}

void Extract_dx(float *dx, float **Cart_mode, float *sigma,
		int N_modes, int N_Cart, float factor, long *idum)
{
  float gasdev(long *idum);
  int i, k; float z, coeff;

  for(i=0; i<N_Cart; i++)dx[i]=0;
  for(k=0; k<N_modes; k++){
    z=gasdev(idum);
    coeff = z*sigma[k];
    float *x=dx, *mode=Cart_mode[k];
    for(i=0; i<N_Cart; i++){*x += coeff*(*mode); x++; mode++;}
  }
  for(i=0; i<N_Cart; i++)dx[i]*=factor;
}

unsigned long randomgenerator(void){

     unsigned long tm;
     time_t seconds;

     time(&seconds);
     srand((unsigned)(seconds % 65536));
     do   /* waiting time equal 1 second */
       tm= clock();
     while (tm/CLOCKS_PER_SEC < (unsigned)(1));
     return((unsigned long) (rand()));
}

void Print_PDB(FILE *file_out, atom *atoms, int natoms, float *atom_str3,
	       struct residue *seq, int k, float rmsd)
{
  int i, ij=0, j=0; atom *atom1=atoms;
  char aaname3[10]; for(i=0; i<4; i++)aaname3[i]='\0';
  if(k>=0)fprintf(file_out, "MODEL %d RMSD from PDB= %.3f\n", k, rmsd);
  for(i=0; i<natoms; i++){
    //ALIGN
    if(atom1->ali<0){atom1++; continue;}
    int res=atom1->res; Name3(aaname3,seq[res].i_aa);
    if(j < 99999)j++;
    fprintf(file_out,
	    "ATOM  %5d  %3s %3s %c%4s    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	    j,atom1->name,aaname3,seq[res].chain, seq[res].pdbres,
	    atom_str3[ij],atom_str3[ij+1],atom_str3[ij+2],1.0,
	    atom1->B_factor);
    ij+=3;
    atom1++;
  }
  fprintf(file_out, "ENDMDL\n");
}

void Print_seqres(FILE *file_out, struct chain *chains, int Nchain)
{
  char aaname3[8]; for(int i=0; i<8; i++)aaname3[i]='\0';
  int line=1;
  for(int i=0; i<Nchain; i++){
    struct chain *ch=chains+i;
    int lines=ch->N_seqres/13+1, k=0;
    for(int l=0; l<lines; l++){
      fprintf(file_out, "SEQRES %3d %c %4d ",
	      line, ch->label, ch->N_seqres);
      for(int j=0; j<13; j++){
	if(k<ch->N_seqres){
	  Name3(aaname3, Code_AA(ch->seqres[k]));
	  fprintf(file_out, " %3s", aaname3);
	}else{
	  fprintf(file_out, "    ");
	}
	k++;
      }
      fprintf(file_out, "\n"); line++;
    }
  }
}

void Write_ref_coord_atom(float *coord_ref, int N_ref, atom *atoms,
                          int *atom_num)
{
  float *coord=coord_ref;
  for(int i=0;i<N_ref;i++){
    double *r=atoms[atom_num[i]].r;
    for(int j=0; j<3; j++){*coord=*r; coord++; r++;}
  }
}

void Write_ref_coord(float *coord_ref, int N_ref, float *coord_all,
		     int *atom_num)
{
  float *coord=coord_ref;
  for(int i=0;i<N_ref;i++){
    float *r=coord_all+3*atom_num[i];
    for(int j=0; j<3; j++){*coord=*r; coord++; r++;}
  }
}

int Put_coord(float *coord_all, struct bond *bonds, int natoms)
{
  float *coord=coord_all; int i;
  for(i=0; i<natoms; i++){
    if(bonds[i].i_atom<0)break;
    for(int j=0; j<3; j++){*coord=bonds[i].r[j]; coord++;}
  }
  return(i);
}

int Put_atoms(atom *atoms, float *coords, int natoms){
  atom *atom=atoms; float *coord=coords;
  for(int i=0; i<natoms; i++){
    for(int j=0; j<3; j++){atom->r[j]=*coord; coord++;}
    atom++;
  }
  return(0);
}


float Energy_clashes(float *r, int N){
  float dthr=1.2, dthr2=dthr*dthr;
  float dthr_large=7*dthr, dthr2_large=dthr_large*dthr_large;
  int large_skip=10, r_skip=3*large_skip;
  double E_sum=0;
  int i, j, k;
  float *ri=r;
  for(i=0; i<N; i++){
    float *rj=r;
    for(j=0; j<i-3; j++){
      float dd=0, d;
      for(k=0; k<3; k++){
	d=ri[k]-rj[k];
	if(fabs(d)> dthr_large){
	  j+=large_skip; rj+=r_skip; goto next_j;
	}
	if(fabs(d)> dthr)goto next_j;
	dd += d*d;
	if(dd > dthr2_large){
	  j+=large_skip; rj+=r_skip; goto next_j;
	}
	if(dd > dthr2)goto next_j;
      }
      if(dd==0){
	printf("ERROR, zero distance for atoms %d and %d\n", i, j);
	exit(8);
      }else{
	E_sum+=1./(dd*dd*dd);
      }
    next_j: rj+=3;
    }
    ri+=3;
  }
  return(E_sum/N);
}


void Write_all_coord_atom(float *coords, atom *atoms, int natoms)
{
  int i, j; float *coord=coords;
  for(i=0; i<natoms; i++){
    double *r=atoms[i].r;
    for(j=0; j<3; j++){*coord=*r; coord++; r++;}
  }
}

int Simulate_ensemble(int N_struct, float factor, char *name,
		      struct Para_simul Para_simul, float Mass_sqrt,
		      atom *atoms, int natoms,
		      struct axe *axe, int naxes,
		      struct bond *bonds,
		      struct residue *seq, int nres,
		      int N_diso, struct Normal_Mode NM,
		      struct Jacobian J, struct Reference Ref,
		      struct interaction *Int_list, int N_int,
		      struct interaction **Int_KB,
		      char *INT_TYPE, float s0, char *nameout,
		      int anhar)
{
  float E_thr=factor*factor*log(4*N_struct);
  float RMSD_mod_thr=0.02;
  float RMSD_min_fact=0.1;
  float sigma_thr=Mass_sqrt*RMSD_mod_thr/factor;
  float c_thr=Mass_sqrt*RMSD_mod_thr;
  float gasdev(long *idum);

  float *sigma2; int i, k;
  if(anhar){sigma2=NM.sigma2_anhar;}
  else{sigma2=NM.sigma2;}
  float sigma[NM.N];
  for(k=0; k<NM.N; k++){
    if(sigma2[k]>0){sigma[k]=sqrt(sigma2[k]);}
    else{sigma[k]=0;}
  }
  double RMSD_all=0;
  for(k=0; k<NM.N; k++){
    if(sigma[k]<sigma_thr)break;
    RMSD_all+=sigma2[k];
  }
  RMSD_all=sqrt(RMSD_all)/Mass_sqrt;
  float RMSD_thr=factor*RMSD_all*RMSD_min_fact;

  idum=randomgenerator();
  int n3=3*natoms;
  float mass[natoms]; for(i=0; i<natoms; i++)mass[i]=Mass(atoms+i);
  float coord_all_1[n3], coord_all[n3];
  Write_all_coord_atom(coord_all_1, atoms, natoms);

  // Open files
  char name_pdb[100], name_out[200];
  sprintf(name_pdb, "%s_%s_ensemble.pdb", nameout, name);
  FILE *file_pdb=fopen(name_pdb, "w");
  sprintf(name_out, "%s_%s_ensemble.dat", nameout, name);
  FILE *file_out=fopen(name_out, "w");
  fprintf(file_out, "#struct rmsd energy modes accept\n");
  if(anhar){
    fprintf(file_out, "# Omega2 corrected for anharmonicity\n");
    fprintf(file_pdb, "REMARK Omega2 corrected for anharmonicity\n");
  }

  float d_phi[naxes];
  double delta_phi[naxes]; for(i=0; i<naxes; i++)delta_phi[i]=0;
  float Ene_anhar,
    Ene_anhar_ini=Energy_anharmonic(coord_all_1,atoms,natoms,
				    nres,Int_list,N_int,Int_KB,20,seq,
				    axe, naxes, delta_phi);
  // Write initial PDB
  //fprintf(file_pdb, "MODEL 0  RMSD from PDB= 0.0\n");
  //Print_PDB(file_pdb, atoms, natoms, coord_all_1, seq, 0, 0.00);
  fprintf(file_out, "%d %.3f %.3g 0 1\n", 0, 0.0, 0.0);

  int nsel=0, ndisc=0, ndisc_max=50, nconf=1;
  struct bond *bonds2=Set_bonds_topology(natoms, atoms, seq);
  //struct ali_atoms ali_a; Copy_ali(&ali_a, Ref);
  while(nsel<N_struct){
    if(ndisc>=ndisc_max)break;
    Set_bonds_measure(bonds, natoms, atoms);
    int nmode=0;
    for(k=0; k<NM.N; k++){
      // Change the structure with each relevant normal mode
      if(sigma[k]<sigma_thr)break;
      float c = gasdev(&idum)*factor*sigma[k], *mode=NM.Tors[k];
      if(c<c_thr)continue;
      nmode++;
      for(i=0; i<naxes; i++){d_phi[i]=c*mode[i]; delta_phi[i]=d_phi[i];}
      Copy_bonds(bonds2, bonds, natoms);
      Build_up(bonds2, natoms, d_phi, naxes);
      Put_coord(coord_all, bonds2, natoms);
      Ene_anhar=Energy_anharmonic(coord_all,atoms,natoms,
				  nres,Int_list,N_int,Int_KB,20,seq,
				  axe, naxes, delta_phi);
      Ene_anhar-= Ene_anhar_ini;
      if(Ene_anhar<=E_thr)Copy_bonds(bonds, bonds2, natoms); // accept
    }
    Put_coord(coord_all, bonds, natoms);
    float rmsd=rmsd_mclachlan_f(coord_all_1, coord_all, mass, natoms);
    Ene_anhar=Energy_anharmonic(coord_all,atoms,natoms,
				nres,Int_list,N_int,Int_KB,20,seq,
				axe, naxes, delta_phi);
    Ene_anhar-= Ene_anhar_ini;
    fprintf(file_out, "%d %.3f %.3g %d", nconf, rmsd, Ene_anhar, nmode);
    nconf++;
    if(Ene_anhar>E_thr || rmsd<RMSD_thr){
      ndisc++; fprintf(file_out, " 0\n"); continue;
    }
    nsel++; fprintf(file_out, " 1\n");
    Print_PDB(file_pdb, atoms, natoms, coord_all, seq, nsel, rmsd);
  }
  printf("Distances printed in %s\n", name_out); fclose(file_out);
  if(nsel)printf("%d simulated structures printed in %s\n", nsel, name_pdb);
  fclose(file_pdb);
  Set_bonds_measure(bonds, natoms, atoms);
  //if(Ini_ch==0){Empty_tors(Diff_s); free(ali_sim.alignres); Ini_ch=1;}
  if(ndisc>=ndisc_max){
    printf("WARNING, amplitude factor= %.1f", factor);
    printf(" only %d over %d structure printed in ensemble\n",nsel,N_struct);
    printf("Rejection ratio= %.1f\n", (float)nconf/nsel-1);
    return(-1);
  }
  return(nsel);
}

void Simulate_ensemble_Cart(int N_struct, float factor, char *name,
			    struct Para_simul Para_simul,
			    atom *atoms, int natoms,
			    struct axe *axe, int naxes,
			    struct bond *bonds,
			    struct residue *seq, int nres,
			    int N_diso, struct Normal_Mode NM,
			    struct Jacobian J, struct Reference Ref,
			    struct interaction *Int_list, int N_int,
			    char *INT_TYPE, float s0, char *nameout,
			    int anhar)
{
  float *sigma2;
  if(anhar){sigma2=NM.sigma2_anhar;}
  else{sigma2=NM.sigma2;}

  if(natoms!=Ref.N_ref){
    printf("WARNING, cannot simulate Cartesian ensemble\n"); return;
  }
  idum=IDUM;
  int i, k, n3=3*natoms, N_Cart=Ref.N_Cart, N_ref=Ref.N_ref;

  atom atoms_sim[natoms]; float mass[natoms];
  for(i=0; i<natoms; i++){
    atoms_sim[i]=atoms[i];
    mass[i]=Mass(atoms+i);
  }
  float coord_all_1[n3], coord_all[n3];
  Write_all_coord_atom(coord_all_1, atoms, natoms);
  float coord_ref_1[N_Cart]; // coord_ref[N_Cart];
  Write_ref_coord_atom(coord_ref_1, N_ref, atoms, Ref.atom_num);

  // Open files
  char name_pdb[100], name_out[200];
  sprintf(name_pdb, "%s_%s_ensemble_Cart.pdb", nameout, name);
  FILE *file_pdb=fopen(name_pdb, "w");
  sprintf(name_out, "%s_%s_ensemble_Cart.dat", nameout, name);
  FILE *file_out=fopen(name_out, "w");
  fprintf(file_out, "#struct rmsd\n");
  if(anhar){
    fprintf(file_out, "# Omega2 corrected for anharmonicity\n");
    fprintf(file_pdb, "REMARK Omega2 corrected for anharmonicity\n");
  }

  // Initialization for Examine_confchange
  if(Ini_ch){
    Allocate_tors(&Diff_s, naxes, N_Cart, NM.N);
    ali_sim.alignres=malloc(nres*sizeof(int));
    Ini_ch=0;
  }
  ali_sim.seq_id=100; ali_sim.mammoth=0;
  for(i=0; i<nres; i++)ali_sim.alignres[i]=i;
  char summary3[400];
  sprintf(summary3, "Summary_%s_ensemble.dat", nameout);

  struct ali_atoms ali_a; Copy_ali(&ali_a, Ref);
  float dx[N_Cart], sigma[NM.N];
  for(k=0; k<NM.N; k++){
    if(sigma2[k]>0){sigma[k]=sqrt(sigma2[k]);}
    else{sigma[k]=0;}
  }
  // Write initial PDB
  Print_PDB(file_pdb, atoms, natoms, coord_all_1, seq, 0, 0.00);
  fprintf(file_out, "%d %.3f\n", 0, 0.0);
  for(k=1; k<=N_struct; k++){
    Extract_dx(dx, NM.Cart, sigma, NM.N, N_Cart, factor, &idum);
    for(i=0; i<N_Cart; i++){
      coord_all[i]=coord_ref_1[i]+dx[i];
    }
    float rmsd=rmsd_mclachlan_f(coord_all_1, coord_all, mass, natoms);
    fprintf(file_out, "%d %.3f\n", k, rmsd);
    Print_PDB(file_pdb, atoms, natoms, coord_all, seq, k, rmsd);
    Put_atoms(atoms_sim, coord_all, natoms);
    /*Examine_confchange(&Diff_s, bonds, axe, NULL, 0, summary3, nameout, "",
		       atoms, coord_ref_1, seq, nres, N_diso, natoms,
		       atoms_sim, coord_ref, seq, nres, N_diso, natoms,  //&J, 
		       Ref, ali_a, Int_list, N_int, INT_TYPE, s0, NM.N, NM,
		       k, ali_sim, Para_simul, k, NULL, NULL, NULL, 0, anhar);
    */
  }
  printf("Distances printed in %s\n", name_out); fclose(file_out);
  printf("Simulated structures printed in %s\n", name_pdb);
  fclose(file_pdb);
  if(Ini_ch==0){Empty_tors(Diff_s); free(ali_sim.alignres); Ini_ch=1;}
}

void Torsional_confchange(float *diff_phi, char *diff_type,
			  struct bond *bonds_ini,
			  struct Jacobian *J, struct Reference Ref1,
			  struct ali_atoms ali_a,
			  double Tors_fluct,
			  atom *atoms1, int natoms1, char *nameout,
			  struct axe *axes, int naxes,
			  struct residue *seq, int nres,
			  struct chain *chains, int Nchain, 
			  atom *atoms2, int natoms2,
			  struct Para_simul Para_simul,
			  struct Normal_Mode *NM,
			  char *mode_name, int NMODES, float rmsd_thr)
{
  int i, a;
  int IT_MAX=Para_simul.NSTEPS; // 100 
  float STEP_MAX=Para_simul.STEP_MAX, STEP_MIN=Para_simul.STEP_MIN; //0.3 0.001
  float A_MAX=Para_simul.ANGLE, PDB_STEP=0.33;
  printf("Interpolating conformation change with %d %s modes\n",
	 NMODES, mode_name);
  
  // Initial bonds
  Set_bonds_measure(bonds_ini, natoms1, atoms1);
  struct bond bonds[natoms1], bonds_tmp[natoms1];
  Copy_bonds(bonds, bonds_ini, natoms1);
  printf("Internal coordinates of initial structure\n");
  double phi_ini[naxes], phi_opt[naxes]; 
  Internal_coordinates(phi_ini, naxes, bonds_ini, natoms1);
  float diff_phi0[naxes]; for(i=0; i<naxes; i++)diff_phi0[i]=0;

  // Set references 
  /*char SEL[4]="EB"; // "EB" "ALL"
  struct Reference Ref1, Ref2;
  int N_ref=Set_reference(&Ref1, 0, SEL, atoms1, 0, natoms1);
  printf("Change reference atoms to %s n=%d\n",SEL, Ref1.N_ref);
  int N_Cart=3*N_ref; Ref1.N_Cart=N_Cart;
  int first_kin[naxes], last_kin[naxes];
  Change_kin(first_kin, last_kin, axes, naxes, natoms1, Ref1);
  struct Jacobian J; Allocate_Jacobian(&J, naxes, Ref1.N_Cart);
  */
  //int N_ref=Ref1.N_ref;
  int N_ali=ali_a.N_ref, N_Cart=3*N_ali;
  double Mass_tot=0; for(i=0; i<Ref1.N_ref; i++)Mass_tot+=Ref1.mass_atom[i];

  float coord_ref_1[N_Cart], coord_ref_2[N_Cart];
  Write_ref_coord_atom(coord_ref_1, N_ali, atoms1, ali_a.ali1);
  Write_ref_coord_atom(coord_ref_2, N_ali, atoms2, ali_a.ali2);
  float coord1[N_Cart], coord2[N_Cart], coord_pdb[N_Cart]; //coord_kin[N_Cart]; 
  for(i=0; i<N_Cart; i++){
    coord1[i]=coord_ref_1[i]; coord2[i]=coord1[i];
    coord_pdb[i]=coord1[i]; //coord_kin[i]=coord1[i];
  }


  float *coord_new=coord1, *coord_old=coord2, *mass_ref=Ref1.mass_atom;
  float coord_all[3*natoms1]; Put_coord(coord_all, bonds_ini, natoms1);
  float Ene=Energy_clashes(coord_all, natoms1);

  // Print rmsd
  char namermsd[200];
  sprintf(namermsd, "%s_%s_rmsd.dat", nameout, diff_type);
  FILE *file_rmsd=fopen(namermsd, "w");
  printf("Opening %s\n", namermsd);
  fprintf(file_rmsd,
	  "# Interpolating conformation change with %d %s modes\n",
	  NMODES, mode_name);
  fprintf(file_rmsd,"# RMSD_step: max= %.2g min= %.2g max. steps %d\n",
	  STEP_MAX, STEP_MIN, IT_MAX);
  fprintf(file_rmsd,"# E_clash(PDB)= %.2g\n", Ene);
  // Superimpose second structure to first one
  float RMSD_kf=rmsd_mclachlan_f(coord_ref_1, coord_ref_2, ali_a.mass, N_ali);
  float RMSD_ki=0, RMSD_step=0, RMSD_opt=RMSD_kf;
  printf("Initial RMSD: %.2f\n\n", RMSD_kf);
  fprintf(file_rmsd,"# RMSD= %.3f\n", RMSD_kf);
  fprintf(file_rmsd,"# RMSD_ki RMSD_kf RMSD_step E_clash factor accept\n");
  fprintf(file_rmsd, "%.3f\t%.3f\t%.3f\t%.2g\t%.2g\t%d\n",
	  RMSD_ki, RMSD_kf, RMSD_step, Ene, 1.0, 1);

  // Opening file for printing PDB
  FILE *file_pdb=NULL;
  float coord_all_1[3*natoms1], mass[natoms1], rmsd;
  for(i=0; i<natoms1; i++)mass[i]=Mass(atoms1+i);
  if(strcmp(diff_type, "standard")!=0){
    char pdbout[200];
    sprintf(pdbout, "%s_confchange_%s.pdb", nameout, diff_type);
    file_pdb=fopen(pdbout, "w");
    printf("Writing trajectory in PDB format in %s\n", pdbout);
    Print_seqres(file_pdb, chains, Nchain);
    Put_coord(coord_all_1, bonds_ini, natoms1);
    Print_PDB(file_pdb, atoms1, natoms1, coord_all_1, seq, 0, 0.0);
  }

  float c_amax[NMODES];
  for (i=0; i<NMODES; i++){
    float *mode=NM->Tors[i], max=fabs(mode[0]); 
    for(a=1; a<naxes; a++)if(fabs(mode[a])>max)max=fabs(mode[a]);
    c_amax[i]=A_MAX/max; // Maximum coefficient of mode i for phi<A_MAX
  }

  int accept=1, iter, exit=0, k_pdb=0;
  //float FRAC=0.0;
  float Mass_sqrt=sqrt(Mass_tot), c_thr=rmsd_thr*Mass_sqrt;
  double skip_tot=0; int nskip=0;
  for(iter=0; iter<IT_MAX; iter++){

    /*// Update Jacobian matrix WRONG!!!!
      rmsd=rmsd_mclachlan_f(coord_kin, coord_old, mass_ref, N_ali);
    if((rmsd>RMSD_TOL)){
      printf("Computing kinetic energy\n");
      // Kinetic energy must be computed with new atoms!!!!
      for(i=0;i<N_Cart;i++)coord_kin[i]=coord_old[i];
      atom atoms_t[natoms1]; Copy_atoms_bonds(atoms_t, bonds, natoms1);
      Compute_kinetic(J, axes, naxes, atoms_t, natoms1, Ref1, 0);
      }*/

    // Difference of coordinates
    float Diff_Cart[N_Cart];
    for(i=0;i<N_Cart;i++)Diff_Cart[i]=coord_ref_2[i]-coord_old[i];
    float tors_coeff[naxes];
    for(i=0; i<naxes; i++){
	tors_coeff[i]=
	  Scalar_product_weighted(Diff_Cart, J->Jacobian_ar[i],
				  Ref1.mass_coord, N_Cart);
      }

    // Compute angles
    float diff_step[naxes], diff_all[naxes];
    for(a=0; a<naxes; a++)diff_step[a]=0;
    double rstep=0, factor=1, excess=1;
    float c_lin[NMODES];
    for (i=0; i<NMODES; i++){
      float *mode=NM->Tors[i]; 
      double c=0; for(a=0; a<naxes; a++)c+=tors_coeff[a]*mode[a];
      if(fabs(c)>excess*c_amax[i])excess=fabs(c)/c_amax[i];
      c_lin[i]=c; rstep+=c*c;
    }
    rstep=1.1*sqrt(rstep/Mass_tot);
    if(rstep<STEP_MIN){exit=1; goto print;}
    factor=rstep/STEP_MAX;
    if(excess>factor)factor=excess;

    //float thr=FRAC*RMSD_kf*Mass_sqrt; if(thr<c_thr)thr=c_thr;
    int skip=0;
    for (i=0; i<NMODES; i++){
      float c=c_lin[i]; if(factor>1)c/=factor;
      if(fabs(c)<c_thr){skip++; continue;}
      float *mode=NM->Tors[i];
      for(a=0; a<naxes; a++)diff_step[a]+=c*mode[a];
    }
    skip_tot+=skip; nskip++;
    float max=0;
    for(i=0; i<naxes; i++)if(fabs(diff_step[i])>max)max=fabs(diff_step[i]);
    if(max>A_MAX){
      float f=A_MAX/max;
      printf("WARNING, too large angle= %.2f ", max);
      printf("Reducing d_phi by a factor %.2g\n", f);
      rstep/=f; for(i=0; i<naxes; i++)diff_step[i]*=f;
    }
    if(rstep<STEP_MIN){exit=1; goto print;}

    // Move atoms
    float reduce_step=1; factor=1;
    while(1){ // loop on amplitudes
      for(i=0; i<naxes; i++)diff_all[i]=diff_phi0[i]+diff_step[i];
      Copy_bonds(bonds_tmp, bonds_ini, natoms1);
      Build_up(bonds_tmp, natoms1, diff_all, naxes);
      accept=Write_RMSD(&reduce_step, &RMSD_opt, &RMSD_kf,
			STEP_MIN, STEP_MAX, factor, iter, file_rmsd,
			bonds_tmp, coord_all, natoms1, ali_a,
			coord_new,coord_old,coord_ref_1,coord_ref_2);
      if(accept)break;
      factor*=reduce_step; rstep*=reduce_step;
      if(rstep<STEP_MIN){
	printf("Too small amplitude (%.2g), exiting\n", factor);
	accept=0; break;
      }
      printf("Reducing amplitude by factor %.2g\n", reduce_step);
      for(i=0; i<naxes; i++)diff_step[i]*=reduce_step;
    } // End loop amplitude 


    if(accept==0){
      exit=1;
    }else{
      exit=(RMSD_opt < RMSD_TOL);
      // Store conformation
      Copy_bonds(bonds, bonds_tmp, natoms1);
      Internal_coordinates(phi_opt, naxes, bonds, natoms1);
      for(i=0; i<naxes; i++)diff_phi0[i]=phi_opt[i]-phi_ini[i];
    }

    // Print PDB
  print:
    if(file_pdb){
      rmsd=rmsd_mclachlan_f(coord_pdb, coord_new, mass_ref, N_ali);
      if((rmsd > PDB_STEP)|| exit){
	k_pdb++;
	for(i=0; i<N_Cart; i++)coord_pdb[i]=coord_new[i];
	RMSD_ki=rmsd_mclachlan_f(coord_all_1, coord_all, mass, natoms1);
	Print_PDB(file_pdb,atoms1,natoms1,coord_all,seq,k_pdb,RMSD_ki);
      }
    }
    printf("\n");
    if(exit){
      printf("Exiting iterations, rmsd= %.2f\n",RMSD_opt);
      break;
    }

    float *tmp=coord_old; coord_old=coord_new; coord_new=tmp;

  }

  // Take final coordinates from bonds
  fprintf(file_rmsd, "# Final RMSD:\n");
  float reduce_step;
  Write_RMSD(&reduce_step, &RMSD_opt, &RMSD_kf, STEP_MIN, 100, // = STEP_MAX
	     1.00, iter, file_rmsd, bonds, coord_all, natoms1, ali_a,
	     coord_new,coord_ref_1,coord_ref_1,coord_ref_2);
  printf("Final RMSD= %.2f\n", RMSD_kf);

  printf("Storing internal coordinates\n");
  Internal_coordinates(phi_opt, naxes, bonds, natoms1);
  for(i=0; i<naxes; i++)diff_phi[i]=phi_opt[i]-phi_ini[i];
  Periodic_angles(diff_phi, axes, naxes);

  Copy_bonds(bonds, bonds_ini, natoms1);
  Build_up(bonds, natoms1, diff_phi, naxes);
  fprintf(file_rmsd, "# Reconstructed RMSD:\n");
  Write_RMSD(&reduce_step, &RMSD_opt, &RMSD_kf, STEP_MIN, 100,
	     1.00, iter, file_rmsd, bonds, coord_all, natoms1,
	     ali_a, coord_new, coord_ref_1,coord_ref_1,coord_ref_2);
  printf("%d %s modes, reconstructed RMSD= %.2f\n", NMODES,mode_name,RMSD_kf);
  printf("Mean number of skipped modes out of %d: %.1f\n",
	 NMODES, skip_tot/nskip);
  fprintf(file_rmsd, "# Mean number of skipped modes out of %d: %.1f",
	  NMODES, skip_tot/nskip);
  fprintf(file_rmsd, " threshold= %.3f A\n", rmsd_thr);


  fclose(file_rmsd);
  if(file_pdb)fclose(file_pdb);
  //Empty_Ref(&Ref1);
  //Change_kin_back(axes, naxes, first_kin, last_kin);
  //Empty_Jacobian(J);
  if(0)Compute_kinetic(J, axes, naxes, atoms1, natoms1, Ref1,0);
  return;
}

void Torsional_confchange_RRR(float *diff_phi, char *diff_type,
			  struct bond *bonds_ini,
			  struct Jacobian *J, struct Reference Ref1,
			  struct ali_atoms ali_a,
			  double Tors_fluct,
			  atom *atoms1, int natoms1, char *nameout,
			  struct axe *axes, int naxes,
			  struct residue *seq, int nres,
			  struct chain *chains, int Nchain, 
			      atom *atoms2, int natoms2,
			      struct Para_simul Para_simul,
			      struct Normal_Mode *NM, int NMODES)
{
  int IT_MAX=Para_simul.NSTEPS; int i;
  float STEP_MAX=Para_simul.STEP_MAX, STEP_MIN=Para_simul.STEP_MIN;
  float A_MAX=Para_simul.ANGLE;
  int FAIL_MAX=4;
  printf("Interpolating conformation change with steps\n");
  
  // Initial bonds
  Set_bonds_measure(bonds_ini, natoms1, atoms1);
  struct bond bonds[natoms1], bonds_tmp[natoms1];
  Copy_bonds(bonds, bonds_ini, natoms1);
  printf("Internal coordinates of initial structure\n");
  double phi_ini[naxes], phi_opt[naxes]; 
  Internal_coordinates(phi_ini, naxes, bonds_ini, natoms1);
  float diff_phi0[naxes]; for(i=0; i<naxes; i++)diff_phi0[i]=0;

  // Set references 
  /*char SEL[4]="EB"; // "EB" "ALL"
  struct Reference Ref1, Ref2;
  int N_ref=Set_reference(&Ref1, 0, SEL, atoms1, 0, natoms1);
  printf("Change reference atoms to %s n=%d\n",SEL, Ref1.N_ref);
  int N_Cart=3*N_ref; Ref1.N_Cart=N_Cart;
  int first_kin[naxes], last_kin[naxes];
  Change_kin(first_kin, last_kin, axes, naxes, natoms1, Ref1);
  struct Jacobian J; Allocate_Jacobian(&J, naxes, Ref1.N_Cart);
  */
  //int N_ref=Ref1.N_ref;
  int N_ali=ali_a.N_ref, N_Cart=3*N_ali;

  float coord_ref_1[N_Cart], coord_ref_2[N_Cart];
  Write_ref_coord_atom(coord_ref_1, N_ali, atoms1, ali_a.ali1);
  Write_ref_coord_atom(coord_ref_2, N_ali, atoms2, ali_a.ali2);
  float coord1[N_Cart], coord2[N_Cart], coord_pdb[N_Cart];//coord_kin[N_Cart];
  for(i=0; i<N_Cart; i++){
    coord1[i]=coord_ref_1[i]; coord2[i]=coord1[i];
    coord_pdb[i]=coord1[i]; //coord_kin[i]=coord1[i];
  }

  struct Tors Diff; 
  Allocate_tors(&Diff, naxes, N_Cart, naxes); //N_modes

  float *coord_new=coord1, *coord_old=coord2, *mass_ref=Ref1.mass_atom;
  float coord_all[3*natoms1]; Put_coord(coord_all, bonds_ini, natoms1);
  float Ene=Energy_clashes(coord_all, natoms1);

  // Print rmsd
  char namermsd[200];
  sprintf(namermsd, "%s_%s_rmsd.dat", nameout, diff_type);
  FILE *file_rmsd=fopen(namermsd, "w");
  printf("Opening %s\n", namermsd);
  fprintf(file_rmsd,
	  "# Interpolating conformation change with torsional steps\n");
  fprintf(file_rmsd,"# RMSD_step: max= %.2g min= %.2g N.steps %d\n",
	  STEP_MAX, STEP_MIN, IT_MAX);
  //fprintf(file_rmsd,"# max. angular change (end) %.2g\n", D_MAX_MIN);
  fprintf(file_rmsd,"# E_clash(PDB)= %.2g\n", Ene);
  // Superimpose second structure to first one
  float RMSD_kf=rmsd_mclachlan_f(coord_ref_1, coord_ref_2, ali_a.mass, N_ali);
  float RMSD_ki=0, RMSD_step=0, RMSD_opt=RMSD_kf;
  printf("Initial RMSD: %.2f\n\n", RMSD_kf);
  fprintf(file_rmsd,"# RMSD= %.3f Tors_RMSD= %.3f\n", RMSD_kf, Tors_fluct);
  fprintf(file_rmsd,"# RMSD_ki RMSD_kf RMSD_step E_clash factor accept\n");
  fprintf(file_rmsd, "%.3f\t%.3f\t%.3f\t%.2g\t%.2g\t%d\n",
	  RMSD_ki, RMSD_kf, RMSD_step, Ene, 0.0, 1);

  // Opening file for printing PDB
  FILE *file_pdb=NULL;
  float coord_all_1[3*natoms1], mass[natoms1];
  for(i=0; i<natoms1; i++)mass[i]=Mass(atoms1+i);
  if(strcmp(diff_type, "standard")!=0){
    char pdbout[200];
    sprintf(pdbout, "%s_confchange_%s.pdb", nameout, diff_type);
    file_pdb=fopen(pdbout, "w");
    printf("Writing trajectory in PDB format in %s\n", pdbout);
    Print_seqres(file_pdb, chains, Nchain);
    Put_coord(coord_all_1, bonds_ini, natoms1);
    Print_PDB(file_pdb, atoms1, natoms1, coord_all_1, seq, 0, 0.0);
  }
 
  double Mass_tot=0; for(i=0; i<Ref1.N_ref; i++)Mass_tot+=Ref1.mass_atom[i];
  int accept=1, imax=-1, iter, exit, k_fail=0, k_pdb=0,
    fit_fail=0, regression=1; 
  float Lambda=0; int l_Lambda=1; //float rmsd=100
  for(iter=0; iter<IT_MAX; iter++){

    float diff_step[naxes], diff_all[naxes];
    if(regression==1){
      //if(iter)rmsd=rmsd_mclachlan_f(coord_kin, coord_old, mass_ref, N_ali);
      if((iter==0)){//||(rmsd>RMSD_TOL)){
	printf("Computing kinetic energy\n");
	// Kinetic energy must be computed with new atoms!!!!
	//for(i=0;i<N_Cart;i++)coord_kin[i]=coord_old[i];
	atom atoms_t[natoms1]; Copy_atoms_bonds(atoms_t, bonds, natoms1);
	Compute_kinetic(J, axes, naxes, atoms_t, natoms1, Ref1, 0);
	l_Lambda=1;
      }
      for(i=0;i<N_Cart;i++)Diff.Cart[i]=coord_ref_2[i]-coord_old[i];
      if(iter==0){
	if(Convert_cart2torsion_fit(&Diff, Ref1, J, nameout, 'C', &Lambda)<0){
	  printf("WARNING, ridge regression failed\n");
	  accept=0; k_fail++; fit_fail++; goto test;
	}
      }else{
	if(l_Lambda){
	  printf("Computing T+Lambda, Lambda= %.3f\n", Lambda);
	  l_Lambda=0; for(i=0;i<naxes; i++)J->T[i][i]*=(1+Lambda);
	}
	if(Tors_step(&Diff,J->T,J->Jacobian_ar,Ref1.mass_coord,naxes,N_Cart)){
	  printf("WARNING, fast regression failed\n");
	  accept=0; k_fail++; fit_fail++; regression=0; continue;
	}	
      }
      for(i=0; i<naxes; i++)diff_step[i]=Diff.Tors[i];
    }else{
      goto test;
      // Move angles one by one
      imax=Torsional_confchange_1angle(diff_step, naxes, Para_simul.ANGLE,
				       bonds, natoms1, ali_a,
				       coord_ref_1, coord_ref_2,
				       STEP_MIN, RMSD_opt,
				       atoms1, coord_all);
      if(imax<0){
	accept=0; k_fail++; goto test;
      }
      fprintf(file_rmsd, "# 1angle, angle= %.2g i=%d %c:\n",
	      Para_simul.ANGLE, imax, axes[imax].type);
      printf("# angle= %.2g i=%d %c:\n",
	     Para_simul.ANGLE, imax, axes[imax].type);
    }

    float max=0;
    for(i=0; i<naxes; i++)if(fabs(diff_step[i])>max)max=fabs(diff_step[i]);
    if(max>A_MAX){
      float f=A_MAX/max;
      printf("WARNING, too large angle= %.2f ", max);
      printf("Reducing d_phi by a factor %.2g\n", f);
      for(i=0; i<naxes; i++)diff_step[i]*=f;
    }

    float factor=1., reduce_step=1;
    /**/
    while(1){ // loop on amplitudes
      
      // Move atoms
      for(i=0; i<naxes; i++)diff_all[i]=diff_phi0[i]+diff_step[i];
      Copy_bonds(bonds_tmp, bonds_ini, natoms1);
      Build_up(bonds_tmp, natoms1, diff_all, naxes);
      accept=Write_RMSD(&reduce_step, &RMSD_opt, &RMSD_kf,
			STEP_MIN, STEP_MAX, factor, iter, file_rmsd,
			bonds_tmp, coord_all, natoms1, ali_a,
			coord_new,coord_old,coord_ref_1,coord_ref_2);
      if(accept)break;
      factor*=reduce_step;
      if((factor<0.005)||(reduce_step<0.0001)){
	printf("Too small amplitude (%.2g), exiting\n", factor);
	accept=0; break;
      }
      printf("Reducing amplitude by factor %.2g\n", reduce_step);
      if(regression){
	for(i=0; i<naxes; i++)diff_step[i]*=reduce_step;
      }else{
	diff_step[imax]*=reduce_step;
      }

    } // End loop amplitude

      /**/
    /*
      Copy_bonds(bonds_tmp, bonds_ini, natoms1);
      RMSD_kf=Optimize_amplitude(bonds_tmp,&factor,diff_tors,naxes,
				 bonds,natoms1,Ref1, coord_ref_2);
      accept=Write_RMSD(&reduce_step, &RMSD_opt, &RMSD_kf,
			STEP_MIN, STEP_MAX, factor, iter, file_rmsd,
			bonds, coord_all, natoms1, ali_a,
			coord_new,coord_old,coord_ref_1,coord_ref_2);
    
      if(accept && (RMSD_step>STEP_MAX)){ // Reduce factor
	factor*=(STEP_MAX/RMSD_step);
	for(i=0; i<naxes; i++)diff_step[i]*=factor;
	Build_up(bonds, natoms1, diff_tors, naxes);
	accept=Write_RMSD(&reduce_step,&RMSD_opt, &RMSD_kf,
			  STEP_MIN, STEP_MAX, factor, iter, file_rmsd,
			  bonds, coord_all, natoms1, ali_a,
			  coord_new,coord_old,coord_ref_1,coord_ref_2);
      }
      // Switch coordinates
      if(accept)Copy_bonds(bonds_old, bonds, natoms1);
    }
    */
      /**/
  test:
    exit=(RMSD_opt < RMSD_TOL);
    if(accept==0){
      k_fail++;
      if(k_fail>=2)exit=1;
      if(regression){
	fit_fail++; regression=0;
      }else if(fit_fail<FAIL_MAX){
	regression=1;
      }
      if(exit==0)continue;
    }else{
      // Store conformation
      Copy_bonds(bonds, bonds_tmp, natoms1);
      Internal_coordinates(phi_opt, naxes, bonds, natoms1);
      for(i=0; i<naxes; i++)diff_phi0[i]=phi_opt[i]-phi_ini[i];
      k_fail=0;
      if(regression){fit_fail=0;}
      else if(fit_fail<FAIL_MAX){regression=1;}
    }

    // Print PDB
    if(file_pdb){
      float rmsd=rmsd_mclachlan_f(coord_pdb, coord_new, mass_ref, N_ali);
      if((rmsd > Para_simul.PDB_STEP)|| exit){
	k_pdb++;
	for(i=0; i<N_Cart; i++)coord_pdb[i]=coord_new[i];
	RMSD_ki=rmsd_mclachlan_f(coord_all_1, coord_all, mass, natoms1);
	Print_PDB(file_pdb,atoms1,natoms1,coord_all,seq,k_pdb,RMSD_ki);
      }
    }
    printf("\n");
    if(exit){
      printf("Exiting iterations, k_fail= %d fit_fail= %d\n",k_fail,fit_fail);
      break;
    }

    float *tmp=coord_old; coord_old=coord_new; coord_new=tmp;

  }

  // Take final coordinates from bonds
  fprintf(file_rmsd, "# Final RMSD:\n");
  float reduce_step;
  Write_RMSD(&reduce_step, &RMSD_opt, &RMSD_kf, STEP_MIN, 100, // = STEP_MAX
	     1.00, iter, file_rmsd,
	     bonds, coord_all, natoms1, ali_a,
	     coord_new,coord_ref_1,coord_ref_1,coord_ref_2);
  printf("Final RMSD= %.2f\n", RMSD_kf);

  printf("Storing internal coordinates\n");
  Internal_coordinates(phi_opt, naxes, bonds, natoms1);
  //float *diff_phi=malloc(naxes*sizeof(float));
  for(i=0; i<naxes; i++)diff_phi[i]=phi_opt[i]-phi_ini[i];
  Periodic_angles(diff_phi, axes, naxes);

  Copy_bonds(bonds, bonds_ini, natoms1);
  Build_up(bonds, natoms1, diff_phi, naxes);
  fprintf(file_rmsd, "# Reconstructed RMSD:\n");
  Write_RMSD(&reduce_step, &RMSD_opt, &RMSD_kf, STEP_MIN, 100,
	     1.00, iter, file_rmsd, bonds, coord_all, natoms1,
	     ali_a, coord_new, coord_ref_1,coord_ref_1,coord_ref_2);
  printf("Reconstructed RMSD= %.2f\n", RMSD_kf);

  fclose(file_rmsd);
  if(file_pdb)fclose(file_pdb);
  Empty_tors(Diff);
  //Empty_Ref(&Ref1);
  //Change_kin_back(axes, naxes, first_kin, last_kin);
  //Empty_Jacobian(J);
  Compute_kinetic(J, axes, naxes, atoms1, natoms1, Ref1,0);
  return;
}

int Torsional_confchange_1angle(float *diff_tors, int naxes, float angle,
				struct bond *bonds, int natoms1,
				struct ali_atoms ali_a,
				float *coord_ref_1, float *coord_ref_2,
				float STEP_MIN, float RMSD_min, atom *atoms1,
				float *coord_all)
{
  printf("Changing %d internal coordinates one by one\n", naxes);

  // Initialize random numbers
  if(INIRAN==0){
     INIRAN=1;
     unsigned long iran=randomgenerator();
     InitRandom( (RANDOMTYPE)iran);
  }

  int i, i_min=-1, dir_min=0;
  for(i=0; i<naxes; i++)diff_tors[i]=0;
  int tested[naxes]; for(i=0; i<naxes; i++)tested[i]=0;
  int N_to_test=naxes;

  struct bond bonds2[natoms1];
  int N_ref=ali_a.N_ref;
  float coord_new[ali_a.N_Cart];  
  float *mass_ref=ali_a.mass;
  float RMSD_kf, RMSD_step;
  for(int k=0; k<naxes; k++){
    // Change internal coordinate i in positive and negative direction
    int ran=RandomFloating()*N_to_test, iran=0;
    if(ran>=N_to_test)ran=N_to_test-1;
    for(i=0; i<naxes; i++){
      if(tested[i])continue;
      if(iran==ran)break;
      iran++;
    }
    if(i==naxes)break;

    // Test i
    tested[i]=1; N_to_test--;
    for(int dir=-1; dir<=1; dir+=2){
      diff_tors[i]=dir*angle;
      Copy_bonds(bonds2, bonds, natoms1);
      Build_up(bonds2, natoms1, diff_tors, naxes);
      Put_coord(coord_all, bonds2, natoms1);
      Write_ref_coord(coord_new, N_ref, coord_all, ali_a.ali1);
      RMSD_kf=rmsd_mclachlan_f(coord_ref_2, coord_new, mass_ref, N_ref);
      if(RMSD_kf<RMSD_min){
	RMSD_step=rmsd_mclachlan_f(coord_ref_1, coord_new, mass_ref, N_ref);
	if(RMSD_step > STEP_MIN){
	  i_min=i; dir_min=dir; RMSD_min=RMSD_kf; break;
	}
      }
    }
    diff_tors[i]=0.0;
  }

  if(i_min<0){
    printf("\nERROR, internal coordinate %d does not exist\n", i_min);
    return(-1);
  }
  diff_tors[i_min]=dir_min*angle;
  return(i_min);
}



void Torsional_confchange_1angle_old(struct bond *bonds, int natoms1,
				     int naxes, struct Reference Ref1,
				     float *coord_ref_1, float *coord_ref_2,
				     int N_Cart, atom *atoms1,FILE *file_rmsd,
				     struct Para_simul Para_simul)
{
  printf("Changing %d internal coordinates one by one\n", naxes);
  fprintf(file_rmsd, "# Changing %d internal coordinates one by one\n",naxes);

  float coord_all[3*natoms1];
  Put_coord(coord_all, bonds, natoms1);
  struct bond bonds2[natoms1];
  int N_ref=Ref1.N_ref;
  float coord1[N_Cart], coord2[N_Cart];
  Write_ref_coord(coord1, N_ref, coord_all, Ref1.atom_num);
  Write_ref_coord(coord2, N_ref, coord_all, Ref1.atom_num);
  float *coord_new=coord1, *coord_old=coord2;
  
  float *mass_ref=Ref1.mass_atom;
  float RMSD_opt=rmsd_mclachlan_f(coord_ref_2, coord_new, mass_ref, N_ref);
  float RMSD_kf, RMSD_ki, RMSD_step, Ene;

  for(int iter=0; iter<Para_simul.NSTEPS; iter++){ //ITER_MAX_1
    
    int i, j;
    float diff_tors[naxes]; for(j=0; j<naxes; j++)diff_tors[j]=0;
    float RMSD_min=10000000;
    int i_min=-1, dir_min=0, accept=1;

    for(i=0; i<naxes; i++){

      // Change internal coordinate i in positive direction
      diff_tors[i]=Para_simul.ANGLE;
      Copy_bonds(bonds2, bonds, natoms1);
      Build_up(bonds2, natoms1, diff_tors, naxes);
      Put_coord(coord_all, bonds2, natoms1);
      Write_ref_coord(coord_new, N_ref, coord_all, Ref1.atom_num);
      RMSD_kf=rmsd_mclachlan_f(coord_ref_2, coord_new, mass_ref, N_ref);
      if(RMSD_kf<RMSD_min){i_min=i; dir_min=1; RMSD_min=RMSD_kf;}

      // Change internal coordinate i in negative direction
      diff_tors[i]=-Para_simul.ANGLE;
      Copy_bonds(bonds2, bonds, natoms1);
      Build_up(bonds2, natoms1, diff_tors, naxes);
      Put_coord(coord_all, bonds2, natoms1);
      Write_ref_coord(coord_new, N_ref, coord_all, Ref1.atom_num);
      RMSD_kf=rmsd_mclachlan_f(coord_ref_2, coord_new, mass_ref, N_ref);
      if(RMSD_kf<RMSD_min){i_min=i; dir_min=-1; RMSD_min=RMSD_kf;}

      diff_tors[i]=0;
    }

    if(i_min<0){
      printf("\nERROR, internal coordinate %d does not exist\n", i_min);
      return;
    }

    float d_max=Para_simul.ANGLE;
    RMSD_kf=1000;
    while(1){
      diff_tors[i_min]=dir_min*d_max;
      Copy_bonds(bonds2, bonds, natoms1);
      Build_up(bonds2, natoms1, diff_tors, naxes);
      Put_coord(coord_all, bonds2, natoms1);
      Write_ref_coord(coord_new, N_ref, coord_all, Ref1.atom_num);
      RMSD_kf=rmsd_mclachlan_f(coord_ref_2, coord_new, mass_ref, N_ref);
      RMSD_step=rmsd_mclachlan_f(coord_old, coord_new, mass_ref, N_ref);
      RMSD_ki=rmsd_mclachlan_f(coord_ref_1, coord_new, mass_ref, N_ref);
      printf("iter= %d RMSD_kf= %.3f step= %.3f angle= %.2g\n",
	     iter+1, RMSD_kf, RMSD_step, d_max);
      Ene=Energy_clashes(coord_all, natoms1);
      fprintf(file_rmsd, "%.3f\t%.3f\t%.3f\t%.2g\t%.2g\t%d\n",
	      RMSD_ki, RMSD_kf, RMSD_step, d_max, Ene, accept);
      fflush(file_rmsd);
      if(isnan(RMSD_step)||isnan(RMSD_kf)||isnan(RMSD_ki))return;

      if(RMSD_kf < RMSD_opt){
	if(RMSD_step<Para_simul.STEP_MAX){
	  RMSD_opt=RMSD_kf; accept=1;
	  Copy_bonds(bonds, bonds2, natoms1);
	  break;
	}else{
	  accept=0; d_max*=0.98*(Para_simul.STEP_MAX/RMSD_step);
	  if(RMSD_step<Para_simul.STEP_MIN)break;
	}
      }else{
	if(RMSD_step<Para_simul.STEP_MIN)break;
	accept=0; d_max/=2;  // Reduce step
	printf("Reducing max. step to %.2g \n", d_max);
      }
      if(d_max<0.0001)break;
    }

    if((accept==0)||(RMSD_kf < RMSD_TOL))break;


  }
}

void Copy_atoms_bonds(atom *atoms, struct bond *bonds, int natoms){
  atom *atom=atoms; struct bond *bond=bonds;
  for(int i=0; i<natoms; i++){
    double *r1=atom->r, *r2=bond->r;
    for(int j=0; j<3; j++){*r1=*r2; r1++; r2++;}
    atom++; bond++;
  }
}

void Copy_bonds(struct bond *bonds1, struct bond *bonds, int natoms)
{
  struct bond *bond1=bonds1, *bond=bonds;
  bond1->previous=NULL;
  for(int i=0; i<natoms; i++){
    *bond1=*bond;
    //if(bond->i_atom<0)break;
    if(i){
      bond1->previous=bonds1+(int)(bond->previous-bonds);
    }
    bond++; bond1++;
  }
}

void Simulate_confchange(int N_frames, float *tors_dir, float *d_Cart,
			 int PRINT_SIMUL_PDB, float RMSD_NM,
			 float *coord_1, float *coord_2, float rmsd,
			 struct Para_simul Para_simul,
			 atom *atoms, int natoms,
			 struct axe *axe, int naxes,
			 struct bond *bonds,
			 struct residue *seq, int nres,
			 int N_diso, struct Normal_Mode NM,
			 struct Jacobian J, struct Reference Ref,
			 struct interaction *Int_list, int N_int,
			 char *INT_TYPE, float s0, char *nameout,
			 double tors2, int anhar)
{
  int i, Ini_ch=1, n3=3*natoms;
  float coord_all[n3], coord_opt[n3];
  float coord_ref[Ref.N_Cart], coord_ref2[Ref.N_Cart];
  for(i=0; i<Ref.N_Cart; i++)coord_ref2[i]=coord_1[i];

  // Compute factor for simulations
  double factor=1.00;
  if(rmsd >0 ){
    factor=sqrt(rmsd/RMSD_NM);
    printf("RMSD= %.2f (observed) factor=%.2f\n", rmsd, factor);
  }


  // Open files
  FILE *file_pdb=NULL, *file_out; char name_pdb[100], name_out[200];
  if(tors_dir==NULL){
      if(PRINT_SIMUL_PDB)sprintf(name_pdb, "%s_simul.pdb", nameout);
      sprintf(name_out, "%s_simul.dat", nameout);
  }else{
    if(PRINT_SIMUL_PDB)sprintf(name_pdb, "%s_pred_mut.pdb", nameout);
    sprintf(name_out, "%s_pred_mut.dat", nameout);
  }
  if(PRINT_SIMUL_PDB)file_pdb=fopen(name_pdb, "w");
  file_out=fopen(name_out, "w");
  fprintf(file_out, "#step rmsd_previous rmsd_ini rmsd_end E_clash\n");

  // Simulations
  float rmsd0=-1;
  if(coord_2){
    rmsd0=rmsd_mclachlan_f(coord_1, coord_2, Ref.mass_atom, Ref.N_ref);
  }

  char type='R'; // Random ensemble
  float coeff=3*factor/N_frames;
  float d_phi[naxes];
  if(tors_dir){   // Move along a predefined torsional direction   
    for(i=0; i<naxes; i++)d_phi[i]=tors_dir[i]*coeff;
    type='T'; // torsional displacement
  }
  Make_conformations(type, N_frames, d_phi, NULL, factor, coeff, rmsd0,
		     coord_all, coord_opt, coord_ref, coord_ref2, 
		     PRINT_SIMUL_PDB, file_out, file_pdb,
		     coord_1, coord_2, Para_simul,
		     atoms, natoms, axe, naxes, bonds, seq, nres, N_diso,
		     NM, J, Ref, Int_list, N_int, INT_TYPE, s0, nameout,
		     anhar);
  if(d_Cart){
    Make_conformations('C', N_frames, NULL, d_Cart, factor, coeff, rmsd0,
		       coord_all, coord_opt, coord_ref, coord_ref2, 
		       PRINT_SIMUL_PDB, file_out, file_pdb,
		       coord_1, coord_2, Para_simul,
		       atoms, natoms, axe, naxes, bonds, seq, nres, N_diso,
		       NM, J, Ref, Int_list, N_int, INT_TYPE, s0, nameout,
		       anhar);
  }
  printf("Distances printed in %s\n", name_out); fclose(file_out);
  if(PRINT_SIMUL_PDB){
    printf("Simulated structures printed in %s\n", name_pdb);
    fclose(file_pdb);
  }
  if(Ini_ch==0){Empty_tors(Diff_s); free(ali_sim.alignres); Ini_ch=1;}
}

int Make_conformations(char type, int N_frames, float *d_phi, float *d_Cart,
		       float factor, float coeff, float rmsd0,
		       float *coord_all_f, float *coord_opt_f,
		       float *coord_ref, float *coord_ref2,
		       int PRINT_SIMUL_PDB, FILE *file_out, FILE *file_pdb,
		       float *coord_str1, float *coord_str2,
		       struct Para_simul Para_simul,
		       atom *atoms, int natoms,
		       struct axe *axe, int naxes,
		       struct bond *bonds,
		       struct residue *seq, int nres,
		       int N_diso, struct Normal_Mode NM,
		       struct Jacobian J, struct Reference Ref,
		       struct interaction *Int_list, int N_int,
		       char *INT_TYPE, float s0, char *nameout,
		       int anhar)
{
  /* type: C= Cartesian interpolation T= torsional interpolation R= random */
  float *sigma2;
  if(anhar){sigma2=NM.sigma2_anhar;}
  else{sigma2=NM.sigma2;}

  float Ene=Energy_clashes(coord_str1, Ref.N_ref), Ene_thr=40*Ene;

  // Initialization

  char summary3[400]; sprintf(summary3, "Summary_%s_simul.dat", nameout);
  int k, i; 
  Set_bonds_measure(bonds, natoms, atoms);
  Put_coord(coord_all_f, bonds, natoms);
  Write_ref_coord(coord_ref2, Ref.N_ref, coord_all_f, Ref.atom_num);

  float *sigma=NULL;
  if(type=='R'){
    if(Ini_ch){
      Allocate_tors(&Diff_s, naxes, NM.N_Cart, NM.N);
      ali_sim.alignres=malloc(nres*sizeof(int));
      Ini_ch=0;
      ali_sim.seq_id=100; ali_sim.mammoth=0;
      for(i=0; i<nres; i++)ali_sim.alignres[i]=i;
      idum=IDUM;
    }
    sigma=malloc(NM.N*sizeof(float));
    for(k=0; k<NM.N; k++){
      if(sigma2[k]>0){sigma[k]=sqrt(sigma2[k]);}
      else{sigma[k]=0;}
    }
  }else if(type=='T'){
    float d_thr=0.01, d_max=0;
    for(i=0; i<naxes; i++)if(fabs(d_phi[i])>d_max)d_max=fabs(d_phi[i]);
    if(d_max>d_thr){
      float f=d_thr/d_max; for(i=0; i<naxes; i++)d_phi[i]*=f;
      printf("Reducing angular changes by a factor %.3f\n", f);
    }
  }


  int k_opt=0, nn=N_frames; if(type!='R')nn=1000;
  float rmsd_sim=0, rmsd1=0, rmsd2=rmsd0, rmsd_min=1000, rmsd1_min=0;

  for(k=0; k<=nn; k++){

    if(k==0)goto Print;
    if(type=='C'){
      float ck=k*coeff;
      for(i=0; i<Ref.N_Cart; i++)coord_ref[i]=coord_str1[i]+ck*d_Cart[i];
    }else{
      if(type=='R'){ // Random simulation starting from PDB
	Extract_dphi(d_phi, NM.Tors, sigma, NM.N, naxes, factor, &idum);
      }
      Build_up(bonds, natoms, d_phi, naxes);
      Put_coord(coord_all_f, bonds, natoms);
      Write_ref_coord(coord_ref, Ref.N_ref, coord_all_f, Ref.atom_num);
    }
    Ene=Energy_clashes(coord_ref, Ref.N_ref);

    // RMSD computation
    rmsd_sim=rmsd_mclachlan_f(coord_ref2, coord_ref, Ref.mass_atom, Ref.N_ref);
    for(i=0; i<Ref.N_Cart; i++)coord_ref2[i]=coord_ref[i];
    rmsd1=rmsd_mclachlan_f(coord_str1, coord_ref, Ref.mass_atom, Ref.N_ref);
    if(coord_str2){
      rmsd2=rmsd_mclachlan_f(coord_str2, coord_ref, Ref.mass_atom, Ref.N_ref);
      if(rmsd2<rmsd_min){
	Put_coord(coord_opt_f, bonds, natoms); 
	rmsd_min=rmsd2; rmsd1_min=rmsd1; k_opt=k;
      }
    }

    /*if(type=='R'){     // Projection of conformation change on normal modes
      atom atoms_sim[natoms]; int cc;
      for(i=0; i<natoms; i++)atoms_sim[i]=atoms[i];
      Examine_confchange(&Diff_s, bonds, axe, NULL, 0, summary3, nameout, "", 
			 atoms, coord_str1, seq, nres, N_diso, natoms,
			 atoms_sim, coord_ref, seq, nres, N_diso, natoms,
			 Int_list, N_int, INT_TYPE, s0, NM.N, Ref, NM, //&J, 
			 k, ali_sim, Para_simul, k, 
			 NULL, NULL, NULL, 0, anhar);
			 }*/

  Print:
    fprintf(file_out, "%d %.3f %.3f %.4f %.3g\n",
	    k, rmsd_sim, rmsd1, rmsd2, Ene);
    if((PRINT_SIMUL_PDB)&&((type=='R')||((k==0)&&(type=='T')))){
      if(Ene < Ene_thr){
	Print_PDB(file_pdb, atoms, natoms, coord_all_f, seq, k, rmsd_sim);
      }
    }
    if((type!='R')&&(rmsd2>rmsd_min))break;
  }
  char str[1000];
  sprintf(str, "# Min.RMSD= %.2f Starting_RMSD= %.2f"
	  " Opt_coeff=%.3f factor= %.3f type=%c\n",
	  rmsd_min, rmsd0, coeff*k_opt, factor, type);
  fprintf(file_out, "%s", str); printf("%s", str);
  fprintf(file_out, "&\n");
  if((PRINT_SIMUL_PDB)&&(type!='R')&&(k_opt>=0)){
    if(type=='C'){
      fprintf(file_pdb, "REMARK: predicted Cartesian mutation\n");
    }else{
      fprintf(file_pdb, "REMARK: predicted Torsional mutation\n");
    }
    Print_PDB(file_pdb, atoms, natoms, coord_opt_f, seq, k, rmsd1_min);
  }
  if(sigma)free(sigma);
  return(0);
}

void Change_kin_back(struct axe *axes, int naxes,
		     int *first_kin, int *last_kin){
  int *first=first_kin, *last=last_kin; struct axe *axe=axes;
  for(int i=0; i<naxes; i++){
    axe->first_kin=*first; axe->last_kin=*last; axe++; first++; last++;
  }
}

void Change_kin(int *first_kin, int *last_kin, struct axe *axes,
		int naxes, int natoms1, struct Reference Ref1)
{
  int refatom[natoms1], i; for(i=0; i<natoms1; i++)refatom[i]=-1;
  for(i=0; i<Ref1.N_ref; i++)refatom[Ref1.atom_num[i]]=i;
  int *first=first_kin, *last=last_kin;
  struct axe *axe=axes;
  //printf("#first_atom first_kin last_atom last_kin\n");
  printf("Changing reference atoms of axes\n");
  for(int a=0; a<naxes; a++){
    *first=axe->first_kin; axe->first_kin=-1;
    for(i=axe->first_atom; i<=axe->last_atom; i++){
      if((axe->first_kin<0)&&(refatom[i]>=0)){
	axe->first_kin=refatom[i]; break;
      }
    }
    *last=axe->last_kin; axe->last_kin=-1;
    for(i=axe->last_atom; i>=axe->first_atom; i--){
      if((axe->last_kin<0)&&(refatom[i]>=0)){
	axe->last_kin=refatom[i]; break;
      }
    }
    if((axe->first_kin<0)||(axe->last_kin<0)||(axe->first_kin>axe->last_kin)){
      printf("ERROR, axe %d out of %d\n", a+1, naxes);
      printf("first_kin= %d last_kin= %d first_atom= %d last_atom= %d\n",
	     axe->first_kin, axe->last_kin, axe->first_atom, axe->last_atom);
      printf("Last atom of axe: %s %d %c\n",
	     axe->bond->atom->name, axe->bond->atom->res, axe->bond->atom->aa);
      exit(8);
    }
    axe++; first++; last++;
  }
}


int Check_nan_f(float *x, int N, char *name){
  float *y=x;
  for(int i=0; i<N; i++){
    if(isnan(*y)){
      printf("ERROR, variable %s i= %d over %d is nan\n", name, i, N);
      printf("Reference atom: %d over %d\n", i/3, N/3);
      return(1);
    }
    y++;
  }
  return(0);
}

float RMSD_LIN(float *RMSD_target, float *RMSD_full,
	       struct Jacobian *J, float *Delta_phi,
	       struct Reference Ref,
	       float *coord_ref_1, float *coord_ref_2,
	       float *coord_full)
{
  // Compute new conformation
  int N_Cart=3*Ref.N_ref;
  float RMSD=0, coord_lin[N_Cart];
  for(int i=0; i<N_Cart; i++){
    double X=coord_ref_1[i];
    for(int a=0; a<J->N_axes; a++)
      X+=J->Jacobian_ar[a][i]*Delta_phi[a];
    coord_lin[i]=X;
  }
  if(coord_ref_2){
    *RMSD_target=
      rmsd_mclachlan_f(coord_ref_2, coord_lin, Ref.mass_atom, Ref.N_ref);
  }
  if(coord_ref_1){
    RMSD=rmsd_mclachlan_f(coord_ref_1, coord_lin, Ref.mass_atom, Ref.N_ref);
  }
  if(coord_full){
    float coord_ref[N_Cart];
    Write_ref_coord(coord_ref, Ref.N_ref, coord_full, Ref.atom_num);
    *RMSD_full=
      rmsd_mclachlan_f(coord_ref, coord_lin, Ref.mass_atom, Ref.N_ref);
  }
  //printf("RMSD_lin= %.2f (ini) %.2f (end) %.2f (torsional)\n",
  //	 RMSD, *RMSD_target, *RMSD_full);
  return(RMSD);
}


float Test_buildup(struct bond *bonds, atom *atoms, int natoms, int naxes)
{
  printf("Testing build up\n");
  // Compute new conformation
  float coord_old[3*natoms], coord_new[3*natoms];
  Set_bonds_measure(bonds, natoms, atoms);
  Put_coord(coord_old, bonds, natoms);
  Build_up(bonds, natoms, NULL, naxes);
  Put_coord(coord_new, bonds, natoms);
  float mass[natoms]; for(int i=0; i<natoms; i++)mass[i]=1;
  float RMSD=rmsd_mclachlan_f(coord_old, coord_new, mass, natoms);
  printf("RMSD with no torsion change: %.2g\n", RMSD);
  if(RMSD > 0.5){
    printf("ERROR in Test_buildup, RMSD= %.2f\n", RMSD);
    exit(8);
  }
  return(RMSD);
}

float Switch_bonds(struct bond *bonds, int naxes, char *prots,
		   struct ali_atoms ali_a,
		   atom *atoms1, int natoms1, struct residue *seq1,
		   atom *atoms2, int natoms2, struct residue *seq2)
{
  // Compute chimera conformations switching internal degrees of freedom
  Set_bonds_measure(bonds, natoms1, atoms1);

  struct bond *bonds2=Set_bonds_topology(natoms2, atoms2, seq2);
  Set_bonds_measure(bonds2, natoms2, atoms2);

  int i, k, N_ref=ali_a.N_ref;
  printf("Switching internal coordinates for %d (all) and %d (ref) atoms ",
	 natoms1, N_ref);
  printf("to match %d (all) atoms of str.2\n", natoms2);

  // Compute differences of internal coordinates
  Internal_differences(bonds, bonds2, natoms1, prots);

  // RMSD all atoms
  int ncart1=3*natoms1, ncart2=3*natoms2;
  float coord1[ncart1]; Put_coord(coord1, bonds, natoms1);
  float coord2[ncart2]; Put_coord(coord2, bonds2, natoms2);
  float coord_new[ncart1];
  int natom=natoms1; if(natoms2<natom)natom=natoms2;
  float mass[natom]; for(i=0; i<natom; i++)mass[i]=1;
  float RMSD=rmsd_mclachlan_f(coord1, coord2, mass, natom);

  // Reference coordinates
  int Ncart=3*N_ref;
  float coord_ref1[Ncart], coord_ref2[Ncart], coord_ref[Ncart];
  Write_ref_coord(coord_ref2, N_ref, coord2, ali_a.ali2);
  Write_ref_coord(coord_ref1, N_ref, coord1, ali_a.ali1);
  RMSD=rmsd_mclachlan_f(coord_ref1, coord_ref2, ali_a.mass, N_ref);
  printf("RMSD ini end: %.3f\n", RMSD);

  char name_pdb[100]; sprintf(name_pdb, "Chimera_%s.pdb", prots);
  FILE *file_pdb=fopen(name_pdb, "w");
  Print_PDB(file_pdb,atoms1,natoms1,coord1,seq1,0,0.00);

  int same=1;
  for(k=1; k<=3; k++){
    printf("Switching internal coordinates of type <= %d\n", k);
    Set_bonds_measure(bonds, natoms1, atoms1);
    Change_internal(NULL, 0, bonds, bonds2, natoms1, k, same);
    Build_up(bonds, natoms1, NULL, naxes);
    Put_coord(coord_new, bonds, natoms1);
    Write_ref_coord(coord_ref, N_ref, coord_new, ali_a.ali1);
    RMSD=rmsd_mclachlan_f(coord_ref1, coord_ref, mass, N_ref);
    printf("RMSD with struct 1: %.3f\n", RMSD);
    RMSD=rmsd_mclachlan_f(coord_ref2, coord_ref, mass, N_ref);
    printf("RMSD with struct 2: %.3f\n", RMSD);
    rmsd_mclachlan_f(coord1, coord_new, mass, natom);
    Print_PDB(file_pdb, atoms1, natoms1, coord_new, seq1, k, RMSD);
  }
  //Print_PDB(file_pdb, atoms2, natoms2, coord2, seq2, k, 0.00);
  printf("Writing %s\n", name_pdb);
  fclose(file_pdb);

  //Compare_build_up(bonds, natoms1, NULL, naxes, bonds2);
  free(bonds2);
  return(RMSD);
}

int Change_internal(float *diff_phi, int naxe,
		    struct bond *bonds, struct bond *bonds2,
		    int natoms, int type, int same)
{
  // type=0: only dofs active in the TNM. Write differences without changing.
  // type=1: change torsion 2: also change angles 3: also change length
  // Pseudo-bonds coordinates (disordered loops) are always changed
  // Same: bonds2 corresponds to same structure, without gaps

  if(diff_phi){
    for(int i=0; i<naxe; i++)diff_phi[i]=0;
  }else if(type==0){
    printf("WARNING in Change_internal, vector for writing differences ");
    printf("has not been provided\n"); return(-1);
  }

  struct bond *bond1=bonds, *b1, *bond2=bonds2; int n=0, m=0;
  for(int i=0; i<natoms; i++){
    bond1=bonds+i;
    if((same==0) && bond1->atom->ali<0)continue;
    //if(bond1->atom->main==0)continue;
    b1=bond1->previous;
    if(same){bond2=bonds2+i;}
    else{bond2=bonds2+bond1->atom->ali;}
    if(b1 && (same || b1->atom->ali>=0)){
      m++;
      if((type >= 3)||(bond1->len >=0)){ // bond length
	if(type){
	  bond1->l=bond2->l; n++;
	}else if(bond1->len >=0){
	  diff_phi[bond1->len]=bond2->l-bond1->l; n++;
	}
      }
      b1=b1->previous;
      if(b1 && (same || b1->atom->ali>=0)){
	m++;
	if((type >= 2)||(bond1->angle >=0)){   // bond angle
	  if(type){
	    bond1->r_cos_theta=bond2->r_cos_theta*(bond1->l/bond2->l);
	    n++;
	  }else if(bond1->angle >=0){
	    diff_phi[bond1->angle]=
	      acos(bond2->r_cos_theta/bond2->l)-
	      acos(bond1->r_cos_theta/bond1->l);
	    n++;
	  }
	}
	b1=b1->previous;
	if(b1 && (same || b1->atom->ali>=0)){
	  m++;
	  if((type)||(bond1->torsion >=0)){  // torsion angle
	    if(type){
	      bond1->phi=bond2->phi; n++;
	    }else if(bond1->torsion >=0){
	      diff_phi[bond1->torsion]=bond2->phi-bond1->phi; n++;
	    }
	  }
	}
      }
    }
  }
  printf("%d out of %d internal coordinates of type<=%d and pseudobonds changed\n",
	 n, m, type);
  return(0);
}


void  Internal_differences(struct bond *bonds, struct bond *bonds2,
			   int natoms, char *prots)
{
  int PRINT_ALL=0;
  if(PRINT_ALL){
    char name[80]; sprintf(name, "Internal_coordinates_%s.dat", prots);
    FILE *file=fopen(name, "w"); // Write internal coordinates
    fprintf(file,
	    "#1=l1 2=l2 3=cos(th1) 4=cos(th2) 5=phi1 6=phi2 atom res\n");
    for(int i=0; i<natoms; i++){
      struct bond *bond1=bonds+i; if(bond1->atom==NULL)break;
      if(bond1->atom->main==0)continue;
      if(bond1->atom->ali<0)continue;
      struct bond *bond2=bonds2+bond1->atom->ali;
      struct bond *b1=bond1->previous;
      if(b1 &&(b1->atom->ali>=0)){
	b1=b1->previous;
	if(b1 &&(b1->atom->ali>=0)){
	  fprintf(file, "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%d\n", 
		  bond1->l, bond2->l,
		  bond1->r_cos_theta/bond1->l, bond2->r_cos_theta/bond2->l,
		  bond1->phi, bond2->phi, bond1->atom->name, bond1->atom->res);
	}
      }
    }
    fclose(file); printf("Writing %s\n", name);
  }

  int n[3], k, a;
  double diff1[3], diff2[3], diff1_aa[3][20], diff2_aa[3][20], n_aa[3][20];
  //float diff[3][natoms];  //diff[k][n[k]]=d;
  for(k=0; k<3; k++){  
    n[k]=0; diff1[k]=0; diff2[k]=0;
    for(a=0; a<20; a++){
      n_aa[k][a]=0; diff1_aa[k][a]=0; diff2_aa[k][a]=0;
    }
  }
  Atom_type("CA"); // Needed to initialize num_BB
  double internal[3][num_BB], internal_2[3][num_BB];
  int n_int[3][num_BB];
  for(k=0; k<3; k++){
    for(a=0; a<num_BB; a++){
      internal[k][a]=0; internal_2[k][a]=0; n_int[k][a]=0;
    }
  }

  float pi=3.1415925, deg=180.0/pi, twopi=2*pi;
  float d, m;
  struct bond *bond1=bonds, *bond2=bonds2, *b1;
  float omega2[natoms];
  for(int i=0; i<natoms; i++){
    bond1=bonds+i;
    int na=Atom_type(bond1->atom->name);
    if(na<0)continue;
    if((bond1->len>=0)||(bond1->angle>=0))continue;  // exclude pseudo-bonds
    if((bond1->atom->ali<0))continue;
    a=Code_AA(bond1->atom->aa);
    bond2=bonds2+bond1->atom->ali; b1=bond1->previous;
    if(b1 &&(b1->atom->ali>=0)){
      k=0; // bond length
      d=bond1->l; // l1
      internal[k][na]+=d; internal_2[k][na]+=d*d;
      d=bond2->l; // l2
      internal[k][na]+=d; internal_2[k][na]+=d*d; n_int[k][na]+=2;
      d= bond1->l-bond2->l; // difference
      diff1[k]+=d; diff2[k]+=d*d; n[k]++;
      if((a>=0)&&(a<20)){
	diff1_aa[k][a]+=d; diff2_aa[k][a]+=d*d; n_aa[k][a]++;
      }
      b1=b1->previous;
      if(b1 && (b1->atom->ali>=0)){
	k=1;  // bond angle
	float s1, c1=Cos_sin(&s1, bond1->r_cos_theta/bond1->l);
	float s2, c2=Cos_sin(&s2, bond2->r_cos_theta/bond2->l);

	d=acos(c1); // theta1
	internal[k][na]+=d; internal_2[k][na]+=d*d;
	d=acos(c2); // theta2
	internal[k][na]+=d; internal_2[k][na]+=d*d; n_int[k][na]+=2;
	d=asin(s2*c1-s1*c2); // difference
	diff1[k]+=d; diff2[k]+=d*d; n[k]++;
	if((a>=0)&&(a<20)){
	  diff1_aa[k][a]+=d; diff2_aa[k][a]+=d*d; n_aa[k][a]++;
	}
	b1=b1->previous;
	if(b1 && (b1->atom->ali>=0)){
	  k=2;  // torsion angle
	  d=bond1->phi; // phi1
	  internal[k][na]+=d; internal_2[k][na]+=d*d;
	  d=bond2->phi; // phi2
	  internal[k][na]+=d; internal_2[k][na]+=d*d; n_int[k][na]+=2;
	  omega2[i]=bond2->phi;
	  d=bond1->phi-bond2->phi; // difference
	  if(d>pi){d-=twopi;}else if(d<-pi){d+=twopi;}
	  bond1->d_phi=fabs(d)*deg; // difference in degrees
	  diff1[k]+=d; diff2[k]+=d*d; n[k]++;
	  if((a>=0)&&(a<20)){
	    diff1_aa[k][a]+=d; diff2_aa[k][a]+=d*d; n_aa[k][a]++;
	  }
	}
      }
    }
  }

  char name[100]; sprintf(name, "Internal_differences_%s.dat", prots);
  FILE *file_out=fopen(name, "w"); printf("Writing %s\n", name);
  fprintf(file_out, "# atom\tn\t3=len s.e.\t5=angle s.e.\t7=torsion\n"); 
  for(a=0; a<num_BB; a++){
    fprintf(file_out, "# %s\t%d ", BB_name[a], n_int[1][a]);
    for(k=0; k<3; k++){
      Ave_se(&(internal[k][a]), &(internal_2[k][a]), n_int[k][a]);
      fprintf(file_out, "\t%.3f %.3f", internal[k][a], internal_2[k][a]);
    }
    fprintf(file_out, "\n");
  }
  fprintf(file_out, "#\n");
  for(k=0; k<3; k++){
    if(k==0){     fprintf(file_out, "# Bond lengths:   ");}
    else if(k==1){fprintf(file_out, "# Bond angles:    ");}
    else if(k==2){fprintf(file_out, "# Torsion angles: ");}
    diff1[k]/=n[k]; diff2[k]/=n[k];
    d=sqrt((diff2[k]-diff1[k]*diff1[k])/(n[k]-1));
    fprintf(file_out, "<d>= %7.4f s.e.= %.4f t= %5.1f RMSD= %.3f n= %.0f\n",
	    diff1[k], d, diff1[k]/d, sqrt(diff2[k]), (float)n[k]);
  }
  fprintf(file_out, "#aa Bond lengths: 2=mean 3=RMSD ");
  fprintf(file_out, "Bond angles: 4=mean 5=RMSD ");
  fprintf(file_out, "Torsion angles: 6=mean 7=RMSD 8=n\n");
  for(a=0; a<20; a++){
    fprintf(file_out, "%c", AA_code[a]);
    for(k=0; k<3; k++){
      if(n_aa[k][a]){
	m=diff1_aa[k][a]/n_aa[k][a];
	d=sqrt(diff2_aa[k][a]/n_aa[k][a]);
      }else{
	m=0; d=0;
      }
      fprintf(file_out, "\t%.3f\t%.3f", m, d);
    }
    fprintf(file_out, "\t%.0f\n", n_aa[0][a]);
  }

  printf("# List of residues with D_omega > %d\n", D_omega_thr);
  printf("#Atom Res nres d_phi\n");
  fprintf(file_out, "# List of residues with D_omega > %d\n", D_omega_thr);
  fprintf(file_out, "#Atom Res nres d_omega omega1 omega2\n");
  bond1=bonds;
  for(int i=0; i<natoms; i++){
    atom *atom=bond1->atom;
    if((strncmp(atom->name, "CA", 2)==0)&&(bond1->d_phi > D_omega_thr)){
      fprintf(file_out, "%s %c%dch%d  %d %.0f %.0f\n",
	      atom->name, atom->aa, atom->res, atom->chain+1,
	      bond1->d_phi, bond1->phi*deg, omega2[i]*deg);
      printf("%s %c%d %d\n", atom->name, atom->aa, atom->res, bond1->d_phi);
    }
    bond1++;
  }

  fclose(file_out);
}

void Ave_se(double *sum1, double *sum2, int n){
  *sum1/=n;
  *sum2=sqrt((*sum2/n-(*sum1)*(*sum1))/(n-1));
}

float Cos_sin(float *s, float c){
  if(c>1){*s=0; return(1);}
  else if(c<-1){*s=0; return(-1);}
  *s=sqrt(1.-c*c); return(c);
}


void Standardize_bonds(struct bond *bonds, atom *atoms,
		       int natoms, char *nameout1,
		       struct axe *axes, int naxes,
		       struct residue *seq, int nres,
		       struct chain *chains, int Nchain, 
		       double Tors_fluct,
		       struct Para_simul Para_simul)
{
  printf("Standardizing bond lengths and bond angles\n");

  int PR_STD=1; char name[80]="tmp.dat"; FILE *file_out;
  if(PR_STD){
    file_out=fopen(name, "w");
    fprintf(file_out, "#atm delta_r delta_cos_theta\n");
  }

  // Change bond angles and bond lengths in bonds
  Set_bonds_measure(bonds, natoms, atoms);
  int nc=0, i;
  for(i=0; i<natoms; i++){
    struct bond *bond=bonds+i;
    if(bond->previous==NULL)continue;
    if((bond->len>=0)||(bond->angle>=0))continue;  // exclude pseudo-bonds
    int na=Atom_type(bond->atom->name);
    if((na<0)||(na>4))continue;
    if(PR_STD){
      fprintf(file_out, "%s %.3f %.3f\n",
	      bond->atom->name, bond->l-BB_length[na],
	      bond->r_cos_theta/bond->l-cos(BB_angle[na]));
    }
    if(bond->previous->previous){
      bond->r_cos_theta =  BB_length[na]*cos(BB_angle[na]);
    }else{
      bond->r_cos_theta =  BB_length[na];
    }
    bond->l=BB_length[na];
    nc+=2;
  }
  if(PR_STD){fclose(file_out);
  printf("%d internal coordinates changed\n", nc);
  printf("Coordinates written in %s\n", name);}

  // Compute coordinates and RMSD
  float tors[naxes]; for(i=0; i<naxes; i++)tors[i]=0;
  Build_up(bonds, natoms, tors, naxes);

  atom atoms2[natoms], *atom2=atoms2, *atom1=atoms;
  for(i=0; i<natoms; i++){
    *atom2=*atom1; atom2->ali=i; atom1++; atom2++;
  }
  Copy_atoms_bonds(atoms2, bonds, natoms);

  // Adjust torsion angles
  float diff_phi[naxes];
  struct Reference Ref1;
  Set_reference(&Ref1, 0, "EB", atoms, 0, natoms);
  struct ali_atoms ali_a; Copy_ali(&ali_a, Ref1);
  struct Jacobian J;
  Allocate_Jacobian(&J, naxes, Ref1.N_Cart);
  int first_kin[naxes], last_kin[naxes];
  Change_kin(first_kin, last_kin, axes, naxes, natoms, Ref1);
  Compute_kinetic(&J, axes, naxes, atoms, natoms, Ref1, 0);
  Change_kin_back(axes, naxes, first_kin, last_kin);
  Torsional_confchange_RRR(diff_phi, "standard", bonds, &J, Ref1, ali_a,
			   Tors_fluct, atoms2, natoms, nameout1, axes, naxes,
			   seq, nres, chains, Nchain,atoms,natoms,Para_simul,
			   NULL,0);
  Empty_Ref(&Ref1);
  Empty_Jacobian(J);

  // Check
  Test_standard(bonds, natoms);
}

int Test_standard(struct bond *bonds, int natoms){
  int n=0;
  for(int i=0; i<natoms; i++){
    struct bond *bond=bonds+i;
    if((bond->len>=0)||(bond->angle>=0))continue;  // exclude pseudo-bonds
    if(bond->previous==NULL)continue;
    int na=Atom_type(bond->atom->name);
    if((na<0)||(na>4))continue;
    float d1=bond->l-BB_length[na], d2;
    if(bond->previous->previous){
      d2=bond->r_cos_theta/bond->l-cos(BB_angle[na]);
    }else{
      d2=0;
    }
    if((fabs(d1)>0.01)||(fabs(d2)>0.01)){
      printf("WARNING %c%d %s %.3f %.3f\n",
	     bond->atom->aa, bond->atom->res+1, bond->atom->name, d1, d2);
      n++;
    }
  }
  if(n){
    printf("WARNING, %d errors in standardize bonds\n", n); return(n);
  }
  printf("Standard bonds successfully tested\n");
  return(0);
}

int Write_RMSD(float *reduce_step, float *RMSD_opt, float *RMSD_kf,
	       float STEP_MIN, float STEP_MAX, float factor, int iter,
	       FILE *file_rmsd,
	       struct bond *bonds, float *coord_all, int natoms1,
	       struct ali_atoms ali_a, float *coord_new, float *coord_old,
	       float *coord_ref_1, float *coord_ref_2)
{
  Put_coord(coord_all, bonds, natoms1);
  Write_ref_coord(coord_new, ali_a.N_ref, coord_all, ali_a.ali1);
  if(Check_nan_f(coord_new, ali_a.N_Cart, "ref_coord")){
    printf("ERROR nan!!!!\n");
    fprintf(file_rmsd, "# nan\n"); return(0);   
  }
  
  //  Compute rmsd and energy
  float RMSD_ki, RMSD_step;
  RMSD_ki=rmsd_mclachlan_f(coord_ref_1, coord_new, ali_a.mass, ali_a.N_ref);
  RMSD_step=rmsd_mclachlan_f(coord_old, coord_new, ali_a.mass, ali_a.N_ref);
  *RMSD_kf=rmsd_mclachlan_f(coord_ref_2, coord_new, ali_a.mass, ali_a.N_ref);
  
  // Check if RMSD improves
  int accept; //float RMSD_TAR=STEP_MAX*0.3333;
  if(RMSD_step<STEP_MIN){
    if((*RMSD_kf<*RMSD_opt)&&(factor==1.0)){
      accept=0; *reduce_step=2*(STEP_MIN/RMSD_step);
    }else{
      printf("Too small step\n");
      accept=0; *reduce_step=0.0;
    }
  }else if(*RMSD_kf<*RMSD_opt){
    /*if((RMSD_step<RMSD_TAR)&&(factor==1.0)){
      accept=0; *reduce_step=1.1*(RMSD_TAR/RMSD_step);
      }*/
    if((RMSD_step>STEP_MAX)&&(factor==1.0)){
      accept=0; *reduce_step=0.95*(STEP_MAX/RMSD_step);
    }else{
      accept=1; *reduce_step=1.;
    }
  }else{ // RMSD increases
    accept=0; *reduce_step=0.5;
  }
  if(accept)*RMSD_opt=*RMSD_kf;   

  printf("%2d round, rmsd= %.2f %.2f step= %.3f f= %.5f accept=%d\n",
	 iter+1, RMSD_ki, *RMSD_kf, RMSD_step, factor, accept);
  float Ene_new=Energy_clashes(coord_all, natoms1);
  if(accept==0)fprintf(file_rmsd, "#");
  fprintf(file_rmsd, "%.3f\t%.3f\t%.3f\t%.2g\t%.2g\t%d\n",
	  RMSD_ki, *RMSD_kf, RMSD_step, Ene_new, factor, accept);
  fflush(file_rmsd);
  return(accept);
}

float Optimize_amplitude(struct bond *bonds_min, float *f_min,
			 float *Diff_Tors, int naxes,
			 struct bond *bonds,
			 int natoms1, struct Reference Ref1,
			 float *coord_ref_2)
{
  float EPS=0.005, RMSD_min=10000, RMSD_RRR=0; *f_min=-1;
  struct bond bonds_tmp[natoms1];
  float f[3], y[3], fmin=0.01, fmax=20; int i;
  f[0]=1.3; f[1]=1; f[2]=0.7;
  printf("Optimizing amplitude f of torsion angle change:\n");
  for(i=0; i<3; i++){
    y[i]=Compute_RMSD(f[i],Diff_Tors,naxes,
		      bonds,bonds_tmp,natoms1,Ref1,coord_ref_2);
    printf("RMSD=%.5f f=%.5f\n", y[i], f[i]);
    if(f[i]==1)RMSD_RRR=y[i];
    if(y[i]<RMSD_min){
      RMSD_min=y[i]; *f_min=f[i]; Copy_bonds(bonds_min, bonds_tmp, natoms1);
    }
    y[i]=exp(-y[i]);
  }

  for(int iter=0; iter<20; iter++){
    float ff=Find_max_quad(f[0],f[1],f[2],y[0],y[1],y[2],fmin,fmax);
    float yy=Compute_RMSD(ff, Diff_Tors, naxes, bonds, bonds_tmp, natoms1,
			  Ref1, coord_ref_2);
    printf("RMSD=%.5f f=%.5f\n", yy, ff);
    if(yy<RMSD_min){
      RMSD_min=yy; *f_min=ff; Copy_bonds(bonds_min, bonds_tmp, natoms1);
    }else{
      if(iter>10)break;
    }
    yy=exp(-yy);
    if(ff<f[0]){
      if((f[0]-ff)<EPS)break;
      f[2]=f[1]; y[2]=y[1]; f[1]=f[0]; y[1]=y[0]; f[0]=ff; y[0]=yy;
    }else if(ff < f[1]){
      if(((f[1]-ff)<EPS)||((ff-f[0])<EPS))break;
      f[2]=f[1]; y[2]=y[1]; f[1]=ff; y[1]=yy;
    }else if(ff < f[2]){
      if(((f[2]-ff)<EPS)||((ff-f[1])<EPS))break;
      f[0]=f[1]; y[0]=y[1]; f[1]=ff; y[1]=yy;
    }else{
      if((ff-f[2])<EPS)break;
      f[0]=f[1]; y[0]=y[1]; f[1]=f[2]; y[1]=y[2]; f[2]=ff; y[2]=yy;
    }
  }
  printf("RMSD RRR= %.4f factor= %.4g\n", RMSD_RRR, 1.0);
  printf("Minimal RMSD= %.4f factor= %.4g\n", RMSD_min, *f_min);
  return(RMSD_min);
}


float Compute_RMSD(float f, float *Diff_Tors, int naxes,
		   struct bond *bonds_ini, struct bond *bonds, int natoms1,
		   struct Reference Ref1, float *coord_ref_2)
{
  float diff_tors[naxes];
  for(int i=0; i<naxes; i++)diff_tors[i]=f*Diff_Tors[i];
  Copy_bonds(bonds, bonds_ini, natoms1);
  Build_up(bonds, natoms1, diff_tors, naxes);
  float coord_all[3*natoms1];
  Put_coord(coord_all, bonds, natoms1);
  float coord_new[Ref1.N_Cart];
  Write_ref_coord(coord_new, Ref1.N_ref, coord_all, Ref1.atom_num);
  float RMSD=
    rmsd_mclachlan_f(coord_ref_2, coord_new, Ref1.mass_atom, Ref1.N_ref);
  return(RMSD);
}

void Copy_ali(struct ali_atoms *ali_a, struct Reference Ref)
{    
  ali_a->N_ref=Ref.N_ref; ali_a->N_Cart=Ref.N_Cart;
  ali_a->ali1=Ref.atom_num;
  ali_a->ali2=Ref.atom_num;
  ali_a->mass=Ref.mass_atom;
}

int Tors_step(struct Tors *Diff, double **T_Lambda, float **Jacobian_ar,
	      float *mass_coord, int naxes, int N_Cart)
{
  double d_phi[naxes]; int fail=0;
  // Delta_phi=(T+Lambda I)^(-1)(M J^t Delta R)
  double Y[naxes]; int i, a;
  for(a=0; a<naxes; a++){
    Y[a]=0;
    for(i=0; i<N_Cart; i++)
      Y[a]+=mass_coord[i]*Jacobian_ar[a][i]*Diff->Cart[i];
    //Y[a]/=sqrt(T_Lambda[a][a]);
  }
  //(T+Lambda I)= LL^t Lower tridiagonal matrix L[i][j] with j<=i
  double **L=Allocate_mat2_d(naxes, naxes), X[naxes];
  if(choldc(L, T_Lambda, naxes)==0){
    //Solve L L^t d_phi = Y as: LX=Y (forward), L^t d_phi=X (backward)
    Forward_substitution(X, L, Y, naxes);
    Backward_substitution(d_phi, L, X, naxes);
    // The mean error is of order 1e-11
  }else{
    printf("Cholc decomposition failed\n");
    fail=1; for(a=0; a<naxes; a++)Diff->Tors[a]=0;
  }
  Empty_matrix_d(L, naxes);
  if(fail)return(1);
  /*double Yphi=0, phi2=0;
  for(a=0; a<naxes; a++){Yphi+=Y[a]*d_phi[a]; phi2+=d_phi[a]*d_phi[a];}
  float scale=Yphi/phi2;
  printf("Fast torsional step, scale= %.2g Yphi=%.2g phi2=%.2g\n",
	 scale, Yphi, phi2);
	 for(a=0; a<naxes; a++)Diff->Tors[a]=scale*d_phi[a];*/
  printf("Fast torsional step\n");
  for(a=0; a<naxes; a++)Diff->Tors[a]=d_phi[a];

  return(0);
}

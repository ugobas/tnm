#include "nma_para.h"
#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "vector.h"
#include "allocate.h"
#include "nrutil.h"
#include "interactions_tnm.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "diagonalize.h"
#include "ridge_regression.h"

int INT_MAX;
int DBG_NMA=0;
int DEBUG;
int ALL_AXES;
char SEL[10];
int READ_RESTRAINT;
int ini_print;
char chain;
// float dist_CA[L_MAX]; // Squared distance from the center of mass
char file_aniso[200];

/*********************************  INM  ************************************/
void Normal_modes_INM(int N,           //degrees of freedom
		      struct interaction *Int_list, int N_int,
		      atom *atoms, int N_atoms,
		      double **Hessian, float **eigen_vector,
		      float *eigen_value, float *eigen_B,
		      float *B_CA, float *Cart_collectivity);
void Compute_Hessian_INM (double **Hessian, int N,
			  struct interaction *Int_list, int N_int,
			  atom *atoms);
int Compute_Bfact_INM(float *B_CA, float *eigen_B, int N,
		      float *eigen_value, float **eigen_vector);

/*********************************  ANM  ************************************/
void Normal_modes_ANM(// OUTPUT:
		      float *eigen_value, float **eigen_vector,
		      double **Hessian,
		      // INPUT:
		      int N,           //degrees of freedom
		      struct interaction *Int_list, int N_int,
		      atom *atoms, int N_atoms);
void Compute_Hessian_ANM(double **Hessian, int N_ref, int *atom_ref,
			 struct interaction *Int_list, int N_int,
			 atom *atoms, int natoms);
int Compute_Bfact_ANM(float *B_CA, float *eigen_B, float *eigen_value,
		      float **eigen_vector, int N, int N_CA);

/*********************************  TNM  ************************************/
void Normal_modes_TNM(// OUTPUT:
		      float *eigen_value, float **eigen_vector,
		      float **Masswtd_evector, float **Cart_mode,
		      double **Hessian,
		      // INPUT:
		      int N,           //degrees of freedom
		      struct interaction *Int_list, int N_int,
		      atom *atoms, int N_atoms,
		      struct axe *axes, int N_axes,
		      struct chain *chains, int Nchain,
		      float *mass_coord, int N_ref);
extern void Forward_substitution(double *X, double **L, double *Y, int N);
extern void Forward_substitution_f(double *X, float **L, double *Y, int N);
extern void Backward_substitution(double *X, double **L, double *Y, int N);
int Compute_Bfact_TNM(float *B_CA, float *eigen_B, float *eigen_value,
		      float **Cart_mode, int N, int N_CA, float mass);
void Convert_torsion2cart(float *Cart_mode, atom *atoms, float *u_a,
			  struct axe *axes, int N_axes,
			  struct Reference Ref, int ia);
void Convert_torsion2cart_old(float *Cart_mode, atom *atoms, float *u_a,
			      struct axe *axes, int N_axes, int *atom_ref,
			      int N_ref, struct Jacobian *J);
extern int Compute_Jtilde(struct Jacobian *J);
void Axes_path(int *main1, int *main2, int *side1, int *side2,
	       struct axe *axes, int i, int nmain, int naxes);

// Transforming axes
void Principal_axis_frame(atom *atoms, int N_atoms, int *atom_ref,
			  int N_atom_ref, struct axe *axes, int N_axes,
			  int ANISOU);
void Get_rotation(double **inertia, double **rotation);
double Center_of_mass(double *r_sum, atom *atoms, int *atom_ref, int N_ref);
int Inertia_tensor_old(atom *atoms, int *atom_ref, int N, double *r_ave,
		       double **corr, double **inertia);

/************************ Auxiliary computation ****************************/
// Select atoms
int Select_atoms(int *atom_num, atom *atoms, int N_atoms, char *SEL);
  // if ali=1 select only aligned atoms, otherwise select all
int Selection (char *atom_name, atom *atom, char *SEL);
float Mass(atom *atom);
float Mass_residue(atom *atom);
// Other
void Name3(char *aaname3, int i_aa);
int Bivariate_fit(int n, float *z, float *x, float *y, int *w,
		  float *a, float *b, float *c, float *r, float *chi);
void f_sort(float d[], int n, int *i_rank);

/**************************** B factors *****************************/
int Set_B_exp(atom *atom, int N_atoms, int N_res, float *B_expp, char *REF);
void Rescale(float *x, int N, float slope, float offset);
void Set_anisou(atom *atom, int N_atoms, int nca,
		float ***anisou, float ***aniso_predp, char *REF);
int Compare_anisou(float ***aniso_pred, float ***aniso_exp, int nca,
		   char *model, char *inter, char *file_name, char *prot_name,
		   int N, int *w,
		   float *dot_aniso, float *ov_aniso, float *delta_aniso);

/********************************  Output **********************************/
float Output_modes(int N,           //degrees of freedom
		   char *model,     // INM, ANM or TNM
		   char *REF_ATM,   // Reference atoms
		   atom *atoms, int N_atoms,
		   struct residue *seq, int N_res,
		   struct axe *axes, int N_axes,
		   float *B_CA, float *B_CA_exp,
		   int *weight, float *dist_CA,
		   float ***aniso_pred, float ***aniso_exp,
		   int N_MODES, float cutoff,
		   char *inter, char *prot_name,
		   FILE *file_B, FILE *file_sum,
		   float **eigen_vector, float *eigen_value, float *eigen_B,
		   float *Cart_collectivity, float *Tors_collectivity,
		   float **Cart_mode, int *atom_ref, int N_ref);
void Print_B(float *B, int nca, FILE *file_out, char *model, float cc);
void Print_PDB_3(float *Cart_mode, float eigen_value, float eigen_B,
		 atom *atoms, int N_atom_ref, int *atom_ref,
		 struct residue *seq, char *file_name, int ia);
void Compare_modes(float **Cart_mode_ANM, float **Cart_mode_TNM, int N_Cart,
		   float **Tors_mode_ANM, float **Tors_mode_TNM, int N_tors,
		   float *mass_coord, char *prot_name, char *inter);
void Print_modes_old(float *eigen_value, float **eigen_vector, float *eigen_B,
		     float *Cart_collectivity, float *Tors_collectivity,
		     int N, int N_MODES, atom *atoms, int N_atoms, int nca,
		     struct residue *seq, struct axe *axes, int N_axes,
		     char *file_name, char *deg_of_freedom,
		     float **Cart_mode, int *atom_ref, int N_ref);
int Get_tors(struct axe *axe, int naxe, int res, int a_ini);

/***************************************************************************
                               CODES
****************************************************************************/


/*********************************  INM ************************************/
void Normal_modes_INM(int N,           //degrees of freedom
		      struct interaction *Int_list, int N_int,
		      atom *atoms, int N_atoms,
		      double **Hessian, float **eigen_vector,
		      float *eigen_value, float *eigen_B,
		      float *B_CA, float *Cart_collectivity)
{
  int i;

  Compute_Hessian_INM(Hessian, N, Int_list, N_int, atoms);
  d_Diagonalize(N, Hessian, eigen_value, eigen_vector, -1);
  Compute_Bfact_INM(B_CA, eigen_B, N, eigen_value, eigen_vector);
  for(i=0; i<N; i++)
    Cart_collectivity[i]=Collectivity_norm2(eigen_vector[i], N)/N;

}

void Compute_Hessian_INM (double **Hessian, int N,
			  struct interaction *Int_list, int N_int,
			  atom *atoms)
{
  int i, j, n;

  for(i=0; i<N; i++)for(j=0; j<N; j++)Hessian[i][j]=0;

  for(n=0; n<N_int; n++){
    float h=Int_list[n].sec_der;
    i=(atoms+Int_list[n].i1)->res;
    j=(atoms+Int_list[n].i2)->res;
    Hessian[i][j]=-h; Hessian[j][i]=-h;
    Hessian[i][i]+=h; Hessian[j][j]+=h;
  }
  return;
}

int Compute_Bfact_INM(float *B_CA, float *eigen_B, int N,
		      float *eigen_value, float **eigen_vector)
{
  int ia, i; double B, B_ia;
  for(i=0; i<N; i++)eigen_B[i]=0;

  for(i=0; i<N; i++){
    B=0;
    for(ia=0; ia<N; ia++){
      if(eigen_value[ia] < E_MIN)continue;
      B_ia= eigen_vector[ia][i]*eigen_vector[ia][i]/eigen_value[ia];
      B += B_ia; eigen_B[ia]+=B_ia;
    }
    B_CA[i]=B;
  }
  return(0);
}


/*********************************  ANM ************************************/

void Compute_Hessian_ANM  (double **Hessian, int N_ref, int *atom_ref,
			  struct interaction *Int_list, int N_int,
			   atom *atoms, int natoms)
{
  int N=3*N_ref;
  int n, i, j;
  int *i_ref=malloc(natoms*sizeof(int));
  for(i=0; i<natoms; i++)i_ref[i]=-1;
  for(i=0; i<N_ref; i++)i_ref[atom_ref[i]]=i;

  // Computing  elements in double precision
  for(i=0; i<N; i++)for(j=0; j<N; j++)Hessian[i][j]=0;
  for(n=0; n<N_int; n++){
    int i1_pdb=Int_list[n].i1, i2_pdb=Int_list[n].i2;
    int i1=3*i_ref[i1_pdb], i2=3*i_ref[i2_pdb], j1, j2;
    if((i1<0)||(i2<0))continue;
    atom *atom1= atoms+i1_pdb, *atom2=atoms+i2_pdb;
    int i1j1, i1j2, i2j1, i2j2;
    double r21[3], h, norm;
    r21[0]=atom2->r[0]-atom1->r[0];
    r21[1]=atom2->r[1]-atom1->r[1];
    r21[2]=atom2->r[2]-atom1->r[2];
    norm=Int_list[n].sec_der/Scalar_product_3_d(r21, r21);

    for(j1=0; j1<3; j1++){
      i1j1 = i1+j1; i2j1=i2+j1;
      for(j2=j1; j2<3; j2++){
	h = r21[j1]*r21[j2]*norm;
	i2j2 =i2+j2; i1j2=i1+j2;
        Hessian[i1j1][i2j2] =-h;
        Hessian[i2j2][i1j1] =-h;
        Hessian[i1j1][i1j2] +=h;
        Hessian[i2j2][i2j1] +=h;
        if(j1!=j2){
          Hessian[i1j2][i2j1] =-h;
          Hessian[i2j1][i1j2] =-h;
          Hessian[i1j2][i1j1]+= h;
          Hessian[i2j1][i2j2]+= h;
        }
      }
    }
  }
  free(i_ref);
  return;
}


int Compute_Bfact_ANM(float *B_CA, float *eigen_B, float *eigen_value,
		      float **eigen_vector, int N, int N_atoms)
{
  int ia, i_atom, j, k;
  double B, B_ia; float *ev;

  for(ia=0; ia< N; ia++)eigen_B[ia]=0;
  for(i_atom=0; i_atom<N_atoms; i_atom++){
    k = 3*i_atom; B=0;
    for(ia=0; ia< N; ia++){
      if(eigen_value[ia] < E_MIN)continue;
      ev=eigen_vector[ia]+k; B_ia=0;
      for(j=0; j<3; j++){B_ia += (*ev)*(*ev); ev++;}
      B_ia /= eigen_value[ia];
      B += B_ia; eigen_B[ia]+=B_ia;
    }
    B_CA[i_atom]=B;
  }
  return(0);
}


/*********************************  TNM ************************************/
void Transform_tors_modes(float *Tors_mode, float *Masswtd_Tors_mode,
			  float **T_sqrt_tr, float **T_sqrt_inv_tr,
			  int N_modes, int N_axes)
{
  int ib, j;
  float *Mw_ev=Masswtd_Tors_mode;
  float *ev=Tors_mode;

  if(KINETIC==0){
    for(ib=0; ib<N_modes; ib++){
      *ev=*Mw_ev; ev++; Mw_ev++;
    }
  }else if(T_sqrt_inv_tr){
    // ev= T_sqrt_inv Masswtd_ev
    for(ib=0; ib<N_axes; ib++){
      double sum=0; float *t=T_sqrt_inv_tr[ib];
      Mw_ev=Masswtd_Tors_mode;
      for(j=0; j<N_modes; j++){
	sum+=(*Mw_ev)*(*t); Mw_ev++; t++;
      }
      *ev=sum; ev++;
    }
  }else{
    // Solve (L^t)EV = Mw_EV by backward_substitution
    // x_i =(v_i - sum_j>i L_ji x_j)/L_ii
    double x[N_modes];
    int nn=N_modes-1;
    double *x1=x+nn; Mw_ev=Masswtd_Tors_mode+nn; ev=Tors_mode+nn;
    for(ib=nn; ib>=0; ib--){
      *x1=*Mw_ev; double *x2=x1+1; float *t=T_sqrt_tr[ib]+ib+1;
      for(j=ib+1; j<N_modes; j++){*x1-= (*x2)*(*t); t++; x2++;}
      (*x1)/=T_sqrt_tr[ib][ib];
      *ev=*x1; ev--; Mw_ev--; x1--;
    }
  }
}

void Transform_tors_modes_old(float *Tors_mode, float *Masswtd_Tors_mode,
			      double **T_sqrt, double **T_sqrt_inv,
			      int N_modes, int N_axes)
{
  int ib, j;
  float *Mw_ev=Masswtd_Tors_mode;
  float *ev=Tors_mode;

  if(KINETIC==0){
    for(ib=0; ib<N_modes; ib++)ev[ib]=Mw_ev[ib];
  }else if(T_sqrt_inv){
    // ev= T_sqrt_inv Masswtd_ev
    for(ib=0; ib<N_axes; ib++){
      double sum=0;
      for(j=0; j<N_modes; j++)sum+=Mw_ev[j]*T_sqrt_inv[j][ib];
      ev[ib]=sum;
    }
  }else{
    // Solve (L^t)EV = M_EV by backward_substitution
    // x_i =(v_i - sum_j>i L_ji x_j)/L_ii
    double *x=malloc(N_modes*sizeof(double));
    for(ib=N_modes-1; ib>=0; ib--){
      x[ib]=Mw_ev[ib];
      for(j=ib+1; j<N_modes; j++)x[ib]-= x[j]*T_sqrt[j][ib];
      x[ib]/=T_sqrt[ib][ib];
      ev[ib]=x[ib];
    }
    free(x);
  }
}


void Compute_Hessian_TNM(double **Hessian,  // Output
			 double **T_sqrt, float **T_sqrt_inv,
			 struct interaction *Int_list, int N_int,
			 atom *atoms, int N_atoms, struct axe *axes,
			 int nmain, int N_axes, int N_modes,
			 struct chain *chains, int Nchain, int kinetic,
			 float K_OMEGA, float K_PSI,
			 float K_PHI, float K_CHI,
			 float K_BA, float K_BL)
{


  // Compute transformed jacobian
  // Jtilde_ia = sum_b J_ib T_sqrt_inv_ab = eta_ia X r_i - tau_ia
  // eta_ak = sum_b chi_kb v_b T_sqrt_inv_ab
  // tau_ak = sum_b chi_kb shift_b T_sqrt_inv_ab
  //
  int i, j, k, l, n, a;
  for(i=0; i<N_axes; i++)for(j=0; j<N_axes; j++)Hessian[i][j]=0;

  double **eta_iaj[N_axes], **tau_iaj[N_axes];
  for(i=0; i<N_axes; i++){
    eta_iaj[i]=Allocate_mat2_d(N_modes, 3);
    tau_iaj[i]=Allocate_mat2_d(N_modes, 3);
  }

  if(kinetic==0){
    for(a=0; a<N_axes; a++){
      struct axe *a_axe=axes+a;
      for(i=a; i<N_axes; i++){
	double *eta=eta_iaj[i][a], *tau=tau_iaj[i][a];
	for(j=0; j<3; j++){
	  eta[j]=a_axe->rot[j];
	  tau[j]=a_axe->shift[j];
	}
      }
    }

  }else if(T_sqrt_inv){

    for(i=0; i<N_axes; i++){
      struct axe *i_axe=axes+i;
      // eta_ia = sum_{ini_chain<=b<=i} v_b T_sqrt_inv_ab
      for(a=0; a<N_modes; a++){
	double *eta=eta_iaj[i][a], *tau=tau_iaj[i][a];
	if(i_axe->previous >=0){
	  double *eta_prev=eta_iaj[i_axe->previous][a];
	  double *tau_prev=tau_iaj[i_axe->previous][a];
	  for(j=0; j<3; j++){
	    eta[j]=eta_prev[j]; tau[j]=tau_prev[j];
	  }
	}else{
	  for(j=0; j<3; j++){eta[j]=0; tau[j]=0;}
	}
	for(j=0; j<3; j++){
	  eta[j]+=(i_axe->rot[j])*T_sqrt_inv[a][i];
	  tau[j]+=(i_axe->shift[j])*T_sqrt_inv[a][i];
	}
      }
    }
  }else{   // Cholevsky decomposition has been successfully performed
    // Solve L eta_ia = v(i) by forward substitution
    // Note that v_b(i) =0 for b>rotable_i
    printf("Computing Hessian with Cholevsky decomposition\n");
    double X[N_axes], Y[N_axes];
    int main1=-1, main2=-1, side1=-1, side2=-1;
    for(i=0; i<N_axes; i++){

      for(l=0; l<N_axes; l++)Y[l]=0;
      Axes_path(&main1, &main2, &side1, &side2, axes, i, nmain, N_axes);
      double **eta=eta_iaj[i], **tau=tau_iaj[i];

      for(j=0; j<3; j++){
	for(l=main1; l<=main2; l++)Y[l]=axes[l].rot[j];
	if(side1>=0)for(l=side1; l<=side2; l++)Y[l]=axes[l].rot[j];
	Forward_substitution(X, T_sqrt, Y, N_axes);
	for(a=0; a<N_axes; a++)eta[a][j]=X[a];
	//
	for(l=main1; l<=main2; l++)Y[l]=axes[l].shift[j];
	if(side1>=0)for(l=side1; l<=side2; l++)Y[l]=axes[l].shift[j];
	Forward_substitution(X, T_sqrt, Y, N_axes);
	for(a=0; a<N_axes; a++)tau[a][j]=X[a];
      }
    }
  }

  // Compute Hessian in double precision (only upper diagonal)

  double **zero=Allocate_mat2_d(N_modes, 3);
  int N_cont[Nchain], n_inter=0;
  for(i=0; i<Nchain; i++)N_cont[i]=0;

  for(n=0; n<N_int; n++){
    atom *atom1=atoms+Int_list[n].i1;
    atom *atom2=atoms+Int_list[n].i2;

    if(atom1->chain==atom2->chain){
      N_cont[atom1->chain]++; // Intrachain contact
    }else{
      n_inter++;  // Interchain contact
    }

    // Determine set of rotable axes
    int ilast;
    double **eta_1, **tau_1, **eta_2, **tau_2;
    double eta[3], tau[3];
    if(atom1->last_dof_side >= 0){
      ilast=atom1->last_dof_side;
    }else{
      ilast=atom1->last_dof_main;
    }
    if(ilast>=0){
      eta_1=eta_iaj[ilast];
      tau_1=tau_iaj[ilast];
    }else{
      eta_1=zero; tau_1=zero;
    }

    if(atom2->last_dof_side >= 0){
      ilast=atom2->last_dof_side;
    }else{
      ilast=atom2->last_dof_main;
    }
    if(ilast>=0){
      eta_2=eta_iaj[ilast];
      tau_2=tau_iaj[ilast];
    }else{
      eta_2=zero; tau_2=zero;
    }

    // Compute r_ij/|r_ij|
    double r21[3], cr21[3], r1[3], r2[3];
    for(j=0; j<3; j++){
      r1[j]=atom1->r[j];
      r2[j]=atom2->r[j];
      r21[j]=r2[j]-r1[j];
    }
    Vector_product_d(cr21, r2, r1); //cr21= r2 X r1

    double h=Int_list[n].sec_der;
    h/=(r21[0]*r21[0]+r21[1]*r21[1]+r21[2]*r21[2]);

    double dk[N_modes];
    for(k=0; k< N_modes; k++){
      for(j=0; j<3; j++){
	eta[j]=eta_2[k][j]-eta_1[k][j];
	tau[j]=tau_2[k][j]-tau_1[k][j];
      }
      dk[k]=Scalar_product_3_d(cr21, eta)-
	Scalar_product_3_d(r21, tau);
      double Hk=h*dk[k];
      for(l=0; l<=k; l++){
	Hessian[l][k] += Hk*dk[l];
      }
    }
  }

  for(i=0; i<Nchain; i++)
    printf("chain %d (%c): %3d contacts\n", i, chains[i].label, N_cont[i]);
  if(Nchain>1)printf("%d interchain contacts\n", n_inter);

  // Hessian of the torsional potential K(a)*(theta_a-theta0_a)^2
  printf("Computing torsional potential\n");
  if((K_OMEGA==0)&&(K_PHI==0)&&(K_PSI==0)&&(K_CHI==0)&&(K_BA==0)&&(K_BL==0))
    goto upper_diagonal;
  if(T_sqrt_inv){
    for(k=0; k<N_modes; k++){
      float *tk=T_sqrt_inv[k];
      for(l=0; l<=k; l++){
	float *tl=T_sqrt_inv[l];
	double H_Phi=0, H_Psi=0, H_Ome=0, H_Chi=0, H_BA=0, H_BL=0, H_BT=0;
	struct axe *ax=axes;
	for(i=0; i<N_axes; i++){
	  if(ax->type=='f'){H_Phi+=tk[i]*tl[i];}
	  else if(ax->type=='p'){H_Psi+=tk[i]*tl[i];}
	  else if(ax->type=='l'){H_BL+=tk[i]*tl[i];}
	  else if(ax->type=='a'){H_BA+=tk[i]*tl[i];}
	  else if(ax->type=='t'){H_BT+=tk[i]*tl[i];}
	  else if(ax->type=='o'){H_Ome+=tk[i]*tl[i];}
	  else if(ax->type=='s'){H_Chi+=tk[i]*tl[i];}
	  ax++;
	}
	double HH=K_PHI*H_Phi;
	if(H_Psi)HH+=K_PSI*H_Psi;
	if(H_Ome)HH+=K_OMEGA*H_Ome;
	if(H_Chi)HH+=K_CHI*H_Chi;
	if(H_BL)HH+=K_BL*H_BL;
	if(H_BA)HH+=K_BA*H_BA;
	if(H_BT)HH+=K_PHI*H_BT;
	Hessian[l][k]+=HH;
      }
    }
  }else{ // Cholevsky decomposition
    double X[N_axes], Y[N_axes], *Y2[N_axes];
    for(i=0; i<N_axes; i++){
      Y[i]=0; Y2[i]=malloc(N_axes*sizeof(double));
      for(k=0; k<N_axes; k++)Y2[i][k]=0;
    }
    for(i=0; i<N_axes; i++){
      char type=axes[i].type;
      if(type=='f'){Y[i]=K_PHI;}
      else if(type=='p'){Y[i]=K_PSI;}
      else if(type=='l'){Y[i]=K_BL;}
      else if(type=='a'){Y[i]=K_BA;}
      else if(type=='t'){Y[i]=K_PHI;}
      else if(type=='o'){Y[i]=K_OMEGA;}
      else if(type=='s'){Y[i]=K_CHI;}
      if(Y[i]){
	Forward_substitution(X, T_sqrt, Y, N_axes);
	for(k=0; k<N_axes; k++)Y2[k][i]=X[k];
	Y[i]=0;
      }
    }
    for(k=0; k<N_axes; k++){
      Forward_substitution(X, T_sqrt, Y2[k], N_axes);
      for(l=0; l<=k; l++)Hessian[l][k]+=X[l];
    }
    for(i=0; i<N_axes; i++)free(Y2[i]);
  }

  // Upper diagonal; store to float
 upper_diagonal:
  for(k=0; k<N_modes; k++){
    for(l=0; l<k; l++){
      Hessian[k][l]=Hessian[l][k];
    }
  }

  if(DBG_NMA){
    // Print
    if(T_sqrt_inv==NULL){printf("Cholevsky decomposition\n");}
    else{printf("Diagonalization of kinetic energy\n");}
    printf("V''(r^PDB): ");
    for(n=0; n<N_int; n++)
      printf(" %.2g:%.2g", Int_list[n].r0,Int_list[n].sec_der);
    //printf(" %d-%d",Int_list[n].i1, Int_list[n].i2);
    printf("\n");
    int ini=N_modes-20;
    printf("###  Kinetic energy matrix:\n");
    for(k=ini; k<N_modes; k++){
      if(T_sqrt_inv==0)
	for(j=k+1; j<N_axes; j++)T_sqrt[k][j]=0;
      for(l=ini; l<=k; l++){
	double sum=0;
	for(j=0; j<N_axes; j++)
	  sum+=T_sqrt[k][j]*T_sqrt[l][j];
	printf(" %.4g", sum);
      }
      printf("\n");
    }
    if(T_sqrt_inv){
      printf("###  T_sqrt T_sqrt_inv:\n");
      for(k=ini; k<N_modes; k++){
	for(l=ini; l<=k; l++){
	  double sum=0;
	  for(j=0; j<N_axes; j++)
	    sum+=T_sqrt_inv[k][j]*T_sqrt[j][l];
	  printf(" %.4f", sum);
	}
	printf("\n");
      }
    }
    printf("###  |eta_ai|\n");
    for(k=ini; k<N_modes; k++){
      for(i=ini; i<N_modes; i++){
	double sum=0, *eta=eta_iaj[i][k];
	for(j=0; j<3; j++)sum+=eta[j]*eta[j];
	printf("%.2g ",sqrt(sum));
      }
      printf("\n");
    }
    printf("### |tau_ai|   (shift)\n");
    for(k=ini; k<N_modes; k++){
      for(i=ini; i<N_modes; i++){
	double sum=0, *tau=tau_iaj[i][k];
	for(j=0; j<3; j++)sum+=tau[j]*tau[j];
	printf("%.2g ",sqrt(sum));
      }
      printf("  %.3g", sqrt(Scalar_product_3_d(axes[k].shift, axes[k].shift)));
      printf("\n");
    }
    printf("###  Potential energy matrix H~:\n");
    for(k=ini; k<N_modes; k++){
      for(l=ini; l<=k; l++)printf(" %.4g", Hessian[k][l]);
      printf("\n");
    }
    printf("Atoms last dof (mainchain): ");
    for(i=0; i<100; i++)printf("%d ", atoms[i].last_dof_main);
    printf("\n");
    exit(8);
  }

  // Cleaning
  Empty_matrix_d(zero, N_modes);
  for(i=0; i<N_axes; i++){
    Empty_matrix_d(eta_iaj[i], N_modes);
    Empty_matrix_d(tau_iaj[i], N_modes);
  }

}

/*int Set_reference_ali(// Output:
		      int *atom_ref1, float *mass_atom, int *atom_ref2,
		      // Input:
		      int ini_ref, char *SEL, int *alignres,
		      atom *atoms1, int ini_atom1,int natoms1,int ini_res1,
		      atom *atoms2, int ini_atom2,int natoms2,int ini_res2)
{
  int N_ref=0, i, i2=0, iali, n;

  // Set reference atoms for structure 1
  int *aref1=malloc(natoms1*sizeof(int));
  int Nref1=Select_atoms(aref1, atoms1+ini_atom1, natoms1, SEL);

  // Look if aligned
  atom *atom2=atoms2+ini_atom2;
  int *sel1=malloc(natoms1*sizeof(int));
  int *sel2=malloc(natoms2*sizeof(int));
  for(i=0; i<natoms1; i++)sel1[i]=0;
  for(i=0; i<natoms2; i++)sel2[i]=0;
  for(i=0; i<Nref1; i++){
    int i1=aref1[i];
    atom *atom1=atoms1+ini_atom1+i1;
    int res2=alignres[atom1->res-ini_res1];
    if(res2<0)continue; iali=-1; res2+=ini_res2;
    while(i2<natoms2){
      if(atom2->res > res2){  // If CB does not exist, use CA!
	if((iali>=0)&&(sel2[iali]==0)){sel1[i1]=1; sel2[iali]=1;}
	break;
      }else if(atom2->res==res2){
	if(strcmp(atom2->name,atom1->name)==0){
	  sel1[i1]=1; sel2[i2]=1; break;
	}else if((iali<0)&&(strncmp(atom2->name,atom1->name,1)==0)){
	  iali=i2;
	}
      }
      i2++; atom2++;
    }
  }

  n=ini_ref;
  for(i=0; i<natoms1; i++){
    if(sel1[i]==0)continue;
    atom_ref1[n]=ini_atom1+i;
    mass_atom[n]=Mass(atoms1+ini_atom1+i);
    n++;
  }
  n=ini_ref;
  for(i=0; i<natoms2; i++){
    if(sel2[i]==0)continue;
    atom_ref2[n]=ini_atom2+i; n++;
  }
  if(DEBUG){
    printf("Reference atoms:\n");
    for(i=ini_ref; i<n; i++){
      atom *atom1=atoms1+atom_ref1[i];
      atom *atom2=atoms2+atom_ref2[i];
      printf("%s %3d %s %3d\n",atom1->name,atom1->res,
	     atom2->name,atom2->res);
    }
  }


  N_ref=n-ini_ref;
  printf("%d atoms out of %d selected for kinetic energy (%s)\n",
	 N_ref, natoms1, SEL);
  free(aref1); free(sel1); free(sel2);
  return(N_ref);
}
*/

int Set_reference(// Output:
		  struct Reference *Ref,
		  // Input:
		  int ini_ref, char *SEL, atom *atoms,
		  int ini_atom, int natoms)
{
  // Set reference atoms
  Ref->atom_num=malloc(natoms*sizeof(int));
  int N_ref=
    Select_atoms(Ref->atom_num+ini_ref, atoms+ini_atom, natoms, SEL);
  Ref->N_ref=N_ref;
  Ref->N_Cart=3*N_ref;
  printf("%d atoms out of %d selected for kinetic energy (%s)\n",
	 N_ref, natoms, SEL);
  Ref->mass_atom=malloc(N_ref*sizeof(float));
  Ref->mass_coord=malloc(Ref->N_Cart*sizeof(float));
  Ref->mass_sqrt=malloc(Ref->N_Cart*sizeof(float));

  float *ma=Ref->mass_atom, *mc=Ref->mass_coord, *ms=Ref->mass_sqrt;
  int *a_num=Ref->atom_num;
  double mass_sum=0;
  for(int i=ini_ref; i<ini_ref+N_ref; i++){
    (*a_num)+=ini_atom;
    (*ma) = Mass(atoms+(*a_num));
    mass_sum+=(*ma);
    for(int j=0; j<3; j++){*mc = *ma; *ms=sqrt(*ma); mc++; ms++;}
    ma++; a_num++;
  }
  Ref->mass_sum=mass_sum;
  return(N_ref);
}

int Align_references(struct ali_atoms *ali_atoms,
		     struct Reference Ref1, atom *atoms1)
{
  int N_ref=Ref1.N_ref;
  ali_atoms->ali1=malloc(N_ref*sizeof(int));
  ali_atoms->ali2=malloc(N_ref*sizeof(int));
  ali_atoms->mass=malloc(N_ref*sizeof(float));
  int iali=0;
  for(int i=0; i<N_ref; i++){
    int j=atoms1[Ref1.atom_num[i]].ali;
    if(j>=0){
      ali_atoms->ali1[iali]=Ref1.atom_num[i];
      ali_atoms->ali2[iali]=j;
      ali_atoms->mass[iali]=Ref1.mass_atom[i];
      iali++;
    }
  }
  printf("%d aligned reference atoms\n", iali);
  ali_atoms->N_ref=iali; ali_atoms->N_Cart=3*iali;
  return(iali);
}

void Match_references(int *kin_cc, struct Reference Ref1, struct Reference Ref)
{
  // kin_cc[i] = Kinetic atom corresponding to ref atom i in Ref1
  int j=0;
  for(int i=0; i<Ref1.N_ref; i++){
    while(Ref.atom_num[j]!=Ref1.atom_num[i])j++;
    if(j>=Ref.N_ref){
      printf("ERROR, unmatched reference %d atom %d out of %d %d\n",
	     i, Ref1.atom_num[i], Ref1.N_ref, Ref.N_ref); exit(8);
    }
    kin_cc[i]=j;
  }
}

void Empty_Ref(struct Reference *Ref){
  free(Ref->atom_num);
  free(Ref->mass_atom);
  free(Ref->mass_coord);
  free(Ref->mass_sqrt);
}

void Convert_torsion2cart(float *Cart_mode, atom *atoms, float *u_a,
			  struct axe *axes, int N_axes,
			  struct Reference Ref, int ia)
{
  // Computing derivative:
  // dr_ia ~ (sum_k u_a,k*rot_ik) X r_i + (sum_k u_a,k*shift_ik)
  int  i, j, k, jj;

  double dr_daxe_rot_global[3], dr_daxe_shift_global[3];
  for(j=0; j<3; j++){
    dr_daxe_rot_global[j]=0;
    dr_daxe_shift_global[j]=0;
  }
  float *u_ak=u_a;
  struct axe *axe=axes;
  for(k=0; k<N_axes; k++){
    for(j=0; j<3; j++){
      dr_daxe_shift_global[j]+=(*u_ak)*axe->global_shift[j];
      dr_daxe_rot_global[j]  +=(*u_ak)*axe->global_rot[j];
    }
    u_ak++; axe++;
  }

  double dr_daxe_rot[Ref.N_Cart], dr_daxe_shift[Ref.N_Cart];
  for(j=0; j<Ref.N_Cart; j++){
    dr_daxe_rot[j]=0;
    dr_daxe_shift[j]=0;
  }
  u_ak=u_a; axe=axes;
  for(k=0; k<N_axes; k++){
    jj=3*axe->first_kin;
    for(i=axe->first_kin; i<=axe->last_kin; i++){
      for(j=0; j<3; j++){
	dr_daxe_rot[jj]  +=(*u_ak)*axe->rot[j];
	dr_daxe_shift[jj]+=(*u_ak)*axe->shift[j];
	jj++;
      }
    }
    u_ak++; axe++;
  }

  double r[3], dd[3], rot[3]; jj=0;
  for(i=0; i<Ref.N_ref; i++){
    atom *atom=atoms+Ref.atom_num[i];
    for(j=0; j<3; j++){
      r[j]= atom->r[j];
      rot[j]=dr_daxe_rot[jj+j]+dr_daxe_rot_global[j];
    }
    Vector_product_d(dd, rot, r);
    for(j=0; j<3; j++){
      Cart_mode[jj]=dd[j]+dr_daxe_shift[jj]+dr_daxe_shift_global[j];
      jj++;
    }
  }
  Normalize_vector_weighted(Cart_mode, Ref.mass_coord, Ref.N_Cart);

  // UGO 10/08/22 r_i=sqrt(m_i)*x_i;
  float norm=sqrt(Ref.mass_sum/Ref.N_ref);
  for(j=0; j<Ref.N_Cart; j++)Cart_mode[j]*=norm; //Ref.mass_sqrt[j];

  if(0){
    // Test Eckart conditions. 15/03/19 I verified that they are fulfilled
    jj=0; int violation=0; float eps=0.001;
    double d[3], M_tot=0, Md_sum[3]; for(j=0; j<3; j++)Md_sum[j]=0;
    double MRd[3], MRd_sum[3]; for(j=0; j<3; j++)MRd_sum[j]=0;
    for(i=0; i<Ref.N_ref; i++){
      double *ri=atoms[Ref.atom_num[i]].r;
      float m= Ref.mass_atom[i];
      M_tot+=m;
      for(j=0; j<3; j++){
	r[j]=ri[j]; d[j]=Cart_mode[jj];
	Md_sum[j]+=m*Cart_mode[jj]; jj++;	
      }
      Vector_product_d(MRd, r, d);
      for(j=0; j<3; j++){MRd_sum[j]+=m*MRd[j];}
    }
    float eps_m=eps*M_tot;
    for(j=0; j<3; j++){
      if((fabs(Md_sum[j])>eps_m)||(fabs(MRd_sum[j])>eps_m)){
	violation++;
	printf("Eckart violation for mode %d! Md_sum= %.4f MRd_sum= %.4f\n",
	       ia, Md_sum[j]/M_tot, MRd_sum[j]/M_tot);
      }
    }
    if(violation){
      printf("Exiting the program\n"); exit(8);
    }else{
      printf("Mode %d ok\n", ia);
    }
  }
}

void Convert_torsion2cart_old(float *Cart_mode, atom *atoms, float *u_a,
			  struct axe *axes, int N_axes, int *atom_ref,
			  int N_ref, struct Jacobian *J)
{
  // Computing derivative:
  // dr_ia ~ (sum_k u_a,k*rot_ik) X r_i + (sum_k u_a,k*shift_ik)

  int  i, j, k, jj;

  double dr_daxe_rot[J->N_Cart], dr_daxe_shift[J->N_Cart];
  for(j=0; j<J->N_Cart; j++){
    dr_daxe_rot[j]=0; dr_daxe_shift[j]=0;
  }
  float *u_ak=u_a;
  struct axe *axe=axes;
  for(k=0; k<N_axes; k++){
    jj=0;
    for(i=0; i<N_ref; i++){
      if((i>=axe->first_kin)&&(i<=axe->last_kin)){
	for(j=0; j<3; j++){
	  dr_daxe_rot[jj]  +=(*u_ak)*axe->local_rot[j];
	  dr_daxe_shift[jj]+=(*u_ak)*axe->local_shift[j];
	  jj++;
	}
      }else{
	for(j=0; j<3; j++){
	  dr_daxe_rot[jj]  +=(*u_ak)*axe->global_rot[j];
	  dr_daxe_shift[jj]+=(*u_ak)*axe->global_shift[j];
	  jj++;
	}
      }
    }
    u_ak++; axe++;
  }
  
  double r[3], dd[3]; jj=0;
  for(i=0; i<N_ref; i++){
    atom *atom=atoms+atom_ref[i];
    for(j=0; j<3; j++)r[j]= atom->r[j];
    Vector_product_d(dd, dr_daxe_rot+jj, r);
    for(j=0; j<3; j++){
      Cart_mode[jj]=dr_daxe_shift[jj]+dd[j]; jj++;
    }
  }
}

int Convert_cart2torsion1(float *Tors_dev,
			  float *Masswtd_Tors_dev,
			  float *Cart_dev,
			  struct Reference Ref,
			  struct Jacobian *J)
{
  // Formula: (LL^t)Delta_theta = JM Delta_r  with:
  // L=T_sqrt (lower diagonal matrix if obtained from Cholevsky decomposition)
  // T=JMJ^t = LL^t
  // J is the jacobian matrix dr_int_i/d_theta_a
  // M_i is the diagonal mass matrix.
  // ==> L^t Delta_theta = L^{-1} J M Delta_r = Jtilde M Delta_r
  // To solve the equation:
  // Compute x=Jtilde MDelta_r =Masswtd_Tors
  // solve x=L^t Delta_theta (by backward substitution if T_sqrt lower diagonal)
  // Return the torsional fraction,
  // sum_i m_i (J d_theta)^2/sum_i m_i Dr_i^2

  if(J->Jtilde_ar==NULL)Compute_Jtilde(J);

  double sum=0; int i, j, a;

  // Compute Masswtd_Tors_dev = Jtilde M Delta_r
  for(a=0; a<J->N_axes; a++){
    float *Jt=J->Jtilde_ar[a]; sum=0;
    float *r=Cart_dev, *m=Ref.mass_coord;
    for(j=0; j<J->N_Cart; j++)sum+=m[j]*r[j]*Jt[j];
    Masswtd_Tors_dev[a]=sum;
  }
  if(J->T_sqrt_inv){
    // T_sqrt is not lower diagonal!
    for(a=0; a<J->N_axes; a++){
      sum=0;
      for(i=0; i<J->N_axes; i++)sum+=J->T_sqrt_inv[i][a]*Masswtd_Tors_dev[i];
      Tors_dev[a]=sum;
    }
  }else{
    // Exploit that T_sqrt is lower diagonal
    // Backward_substitution
    // theta_i =(x_i - sum_j>i L_ji x_j)/L_ii
    for(i=J->N_axes-1; i>=0; i--){
      Tors_dev[i]=Masswtd_Tors_dev[i];
      for(j=i+1; j<J->N_axes; j++)Tors_dev[i]-=Tors_dev[j]*J->T_sqrt[j][i];
      Tors_dev[i]/=J->T_sqrt[i][i];
    }
  }
  return(0);
}

int Convert_cart2torsion(struct Tors *D, struct Reference Ref,
			 struct Jacobian *J)
{
  // Formula: (LL^t)Delta_theta = JM Delta_r  with:
  // L=T_sqrt (lower diagonal matrix if obtained from Cholevsky decomposition)
  // T=JMJ^t = LL^t
  // J is the jacobian matrix dr_int_i/d_theta_a
  // M_i is the diagonal mass matrix.
  // ==> L^t Delta_theta = L^{-1} J M Delta_r = Jtilde M Delta_r
  // To solve the equation:
  // Compute x=Jtilde MDelta_r =Masswtd_Tors
  // solve x=L^t Delta_theta (by backward substitution if T_sqrt lower diagonal)
  // Return the torsional fraction,
  // sum_i m_i (J d_theta)^2/sum_i m_i Dr_i^2

  double sum=0; int i, j, a;
  if(J->Jtilde_ar==NULL)Compute_Jtilde(J);

  // Compute Masswtd_Tors_dev = Jtilde M Delta_r
  for(a=0; a<J->N_axes; a++){
    float *Jt=J->Jtilde_ar[a], *r=D->Cart, *m=Ref.mass_coord;
    sum=0; for(j=0; j<J->N_Cart; j++)sum+=m[j]*r[j]*Jt[j];
    D->MW_Tors[a]=sum;
  }
  if(J->T_sqrt_inv){
    // T_sqrt is not lower diagonal!
    for(a=0; a<J->N_axes; a++){
      sum=0;
      for(i=0; i<J->N_axes; i++)sum+=J->T_sqrt_inv[i][a]*D->MW_Tors[i];
      D->Tors[a]=sum;
    }
  }else{
    // Exploit that T_sqrt is lower diagonal
    // Backward_substitution
    // theta_i =(x_i - sum_j>i L_ji x_j)/L_ii
    for(i=J->N_axes-1; i>=0; i--){
      D->Tors[i]=D->MW_Tors[i];
      for(j=i+1; j<J->N_axes; j++)D->Tors[i]-=D->Tors[j]*J->T_sqrt[j][i];
      D->Tors[i]/=J->T_sqrt[i][i];
    }
  }
  return(0);
}


int Convert_cart2torsion_fit(struct Tors *D, struct Reference Ref,
			     struct Jacobian *J, char *nameout,
			     char type, float *Lambda)
{
  /* type of fit: C= specific heat M= Max penalty O= Ordinary
     L: known value of Lambda */

  /* Minimize |sum_a J_ia Delta_phi_a - sqrt(m_i) Delta_r_i|^2 where:
     J is the jacobian matrix dr_int_i/d_phi_a
     M_i is the diagonal mass matrix 
     The fit is performed through the ridge regression approach
     Lambda=0: ==> L^t Delta_phi = L^{-1} M Delta_r = Jtilde M Delta_r
     L=T_sqrt (lower diagonal matrix if obtained from Cholevsky decomposition)
     T=JMJ^t = LL^t
  */
  
  int Npar=J->N_axes, Nsam=J->N_Cart, a;
  float **X=malloc(Nsam*sizeof(float *));
  float D_out[Npar], Y[Nsam];
  for(int i=0; i<Nsam; i++){
    float sm=sqrt(Ref.mass_coord[i]);
    Y[i]=sm*D->Cart[i];
    X[i]=malloc(Npar*sizeof(float));
    float *Xi=X[i];
    for(a=0; a<Npar; a++){*Xi=sm*J->Jacobian_ar[a][i]; Xi++;}
  }

  //printf("Fitting %d Cartesian dofs with %d torsion angles", Nsam, Npar);
  //printf(" using ridge regression of type %c\n", type);
  char out[200]; sprintf(out, "%s.Ridge_regression", nameout);
  struct ridge_fit fit; fit.A=malloc(Npar*sizeof(float));
  fit.Lambda=*Lambda;
  float y_pred[Nsam];
  float r=Ridge_regression(&fit, y_pred, D_out, out, X, Y, Nsam, Npar, type);
  Empty_matrix_f(X, Nsam);
  *Lambda=fit.Lambda;
  for(a=0; a<Npar; a++)D->Tors[a]=fit.A[a];
  free(fit.A);

  if(r<-1)return(-1);
  // Compute Masswtd_Tors_dev = T_sqrt Delta_theta
  //Compute_MW(D, J); // WARNING: The J matrix has not the correct T,T_sqrt
  return(0);
}

void Compute_MW(float *MW_Tors, float *Tors, struct Jacobian *J){
  // Compute Masswtd_Tors_dev = T_sqrt Delta_theta
  int N=J->N_axes, a, b;
  for(a=0; a<N; a++){
    double W=0;
    if(J->T_sqrt_inv){ // T_sqrt is not lower diagonal! DMW =L^t Dtheta
      for(b=0; b<N; b++)W+=J->T_sqrt[b][a]*Tors[b];
    }else{  // b<a
      for(b=a; b<N; b++)W+=J->T_sqrt[b][a]*Tors[b];
    }
    MW_Tors[a]=W;
  }
}

float Tors_fraction(struct Tors *Diff, float *mass)
{
  double MRR=0, MRT=0, MRW=0, mtot=0; int i;
  // Torsional part / Total
  if(Diff->RMSD==0){
    for(i=0; i<Diff->N_Cart; i++){
      MRR+=mass[i]*Diff->Cart[i]*Diff->Cart[i]; mtot+=mass[i];
    }
    Diff->RMSD=sqrt(MRR/mtot);
    Diff->M=mtot;
  }else{
    mtot=Diff->M; MRR=Diff->M*Diff->RMSD*Diff->RMSD;
  }
  for(i=0; i<Diff->N_axes; i++){
    MRW += Diff->MW_Tors[i]*Diff->MW_Tors[i];
    MRT += Diff->Tors[i]*Diff->Tors[i];
  }
  Diff->RMSD_W=sqrt(MRW/mtot);
  Diff->RMSD_Tors=sqrt(MRT/Diff->N_axes);
  Diff->RMSD_NoTors=sqrt((MRR-MRW)/mtot);
  return(MRW/MRR);
}

/************************ Auxiliary computation 1 ****************************/

void Principal_axis_frame(atom *atoms, int N_atoms, int *atom_ref,
			  int N_atom_ref,struct axe *axes, int N_axes,
			  int ANISOU)
{
  double **inertia=Allocate_mat2_d(3,3);
  double **rotation=Allocate_mat2_d(3,3);
  double **corr=Allocate_mat2_d(3,3);
  double r_ave[3], ri[3];
  int i, j;

  /******************************************************************/
  /*              Transforming to principal axis frame              */
  /******************************************************************/
  // Compute center of mass and inertia tensor for reference atoms
  Inertia_tensor_old(atoms, atom_ref, N_atom_ref, r_ave, corr, inertia);

  // Transform all atoms and axes to principal axis frame of references
  Get_rotation(inertia, rotation);
  for(i=0; i<N_atoms; i++){
    double *r=atoms[i].r;
    for(j=0; j<3; j++)ri[j]=r[j];
    Subtract_vector_3_d(ri, ri, r_ave);
    Rotate_d(ri, rotation);
    for(j=0; j<3; j++)r[j]=ri[j];
  }
  for(i=0; i<N_axes; i++){
    Rotate_d(axes[i].rot, rotation);
    Subtract_vector_3_d(axes[i].offset, axes[i].offset, r_ave);
    Rotate_d(axes[i].offset, rotation);
  }

  // Transfrom anisotropic structure factors
  if(ANISOU){
    float rot[3][3];
    for(i=0; i<3; i++)for(j=0; j<3; j++)rot[i][j]=rotation[i][j];
    for(i=0; i<N_atoms; i++){
      Transform_matrix(atoms[i].anisou, (float **)rot);
    }
  }

  // Cleaning
  Empty_matrix_d(inertia, 3);
  Empty_matrix_d(corr, 3);
  Empty_matrix_d(rotation, 3);
}

int Inertia_tensor_old(atom *atoms, int *atom_ref, int N, double *r_ave,
		       double **corr, double **inertia)
{
  int i, j, k, i_atom;
  double norm_m=0, diag=0;

  for(i=0; i<3; i++){
    r_ave[i]=0;
    for(j=0; j<3; j++)corr[i][j]=0;
  }

  for(i=0; i<N; i++){
    i_atom=atom_ref[i];
    double m=Mass(atoms+i_atom); norm_m+=m;
    double *r=atoms[i_atom].r;
    for(j=0; j<3; j++){
      double mr = m*r[j];   r_ave[j]+= mr;    // center of mass
      for(k=j; k<3; k++)corr[j][k]+=mr*r[k];   // correlation matrix
    }
  }
  for(j=0; j<3; j++)r_ave[j]/=norm_m;

  for(j=0; j<3; j++){
    for(k=j; k<3; k++){
      corr[j][k] -= norm_m*r_ave[j]*r_ave[k];
    }
    diag += corr[j][j];
  }
  for(j=0; j<3; j++){
    inertia[j][j]=diag-corr[j][j];
    for(k=j+1; k<3; k++){
      inertia[j][k]=-corr[j][k];
      inertia[k][j]=inertia[j][k];
    }
  }

  // Printing
  printf("# Inertia tensor:\n");
  for(j=0; j<3; j++){
    for(k=0; k<3; k++)printf("%6.0f ", inertia[j][k]);
    printf("\n");
  }
  return(0);
}



void Get_rotation(double **inertia, double **rotation)
{
  int i, j;
  double **e_vec=Allocate_mat2_d(3,3), e_val[3], det;

  //svdcmp(inertia, 3, 3, e_val, e_vec);
  dd_Diagonalize(3, inertia, e_val, e_vec, 1);

  // Test determinant
  det=
    e_vec[0][0]*(e_vec[1][1]*e_vec[2][2]-e_vec[1][2]*e_vec[2][1])-
    e_vec[0][1]*(e_vec[1][0]*e_vec[2][2]-e_vec[1][2]*e_vec[2][0])+
    e_vec[0][2]*(e_vec[1][0]*e_vec[2][1]-e_vec[1][1]*e_vec[2][0]);
  if(det<0){
    printf("Determinant = %.3f\n", det);
    for(j=0; j<3; j++)e_vec[0][j]=-e_vec[0][j];
  }

  for(i=0; i<3; i++){
    for(j=0; j<3; j++)rotation[i][j]=e_vec[i][j];
  }

  // Clean
  Empty_matrix_d(e_vec, 3);
}

double Center_of_mass(double *r_sum, atom *atoms, int *atom_ref, int N_ref)
{
  int i, j, i_atom; double m_tot=0;
  for(j=0; j<3; j++)r_sum[j]=0;
  for(i_atom=0; i_atom<N_ref; i_atom++){
    i=atom_ref[i_atom];
    double m=Mass(atoms+i), *r=atoms[i].r;
    m_tot+=m; for(j=0; j<3; j++)r_sum[j]+=m*r[j];
  }
  for(j=0; j<3; j++)r_sum[j]/=m_tot;
  return(m_tot);
}



/************************ Auxiliary computation 2 ****************************/
void Name3(char *aaname3, int i_aa){
  int j, i=4*i_aa; for(j=0; j<3; j++){aaname3[j]=AA3[i]; i++;}
  aaname3[3]='\0';
}

int Select_atoms(int *atom_num,  atom *atoms, int N_atoms, char *SEL)
{
  // if ali=1 select only aligned atoms, otherwise select all
  int i_atom, n=0; atom *atom1=atoms;
  for(i_atom=0; i_atom<N_atoms; i_atom++){
    if(Selection(atom1->name, atom1, SEL)){
	atom_num[n]=i_atom; n++;
    }
    atom1++;
  }
  return(n);
}

int Selection (char *atom_name, atom *atom, char *SEL){
  if(strncmp(SEL, "ALL", 2)==0){               // All atoms
    if((atom_name[0]=='H')&&(HYD_REF==0)){return(0);}
    else{return(1);}
  }else if(strncmp(SEL, "CB", 2)==0){          // CB
    if(strncmp(atom_name,"CB",2)==0)return(1);
    if((atom->aa=='G')&&(strncmp(atom_name,"CA",2)==0))
      return(1);
  }else if(strncmp(atom_name,"CA",2)==0){      // CA or BB or EB
    return(1);
  }else if(strncmp(SEL, "BB", 2)==0){          // BB: N-CA-C
    if(strncmp(atom_name,"N ",2)==0)return(1);
    if(strncmp(atom_name,"C ",2)==0)return(1);
  }else if(strncmp(SEL, "EB", 2)==0){          // EB: N-CA-CB-C-O
    if(strncmp(atom_name,"N ",2)==0)return(1);
    if(strncmp(atom_name,"C ",2)==0)return(1);
    if(strncmp(atom_name,"O ",2)==0)return(1);
    if(strncmp(atom_name,"CB",2)==0)return(1);
  }
  return(0);
}

float Mass_residue(atom *atom){
  char *atom_name=atom->name, aa=atom->aa;

  //modyves: removed some double parentheses below, which gave me some warnings
  if(atom_name[0]=='C'){
    if(atom_name[1]==' '){
      return(12); // 0H
    }else if(atom_name[1]=='A'){
      if(aa!='G')return(13); // Adding 1H mass
    }else if(atom_name[1]=='B'){
      if(aa=='A')return(15);              // Adding 3H mass
      if((aa=='V')||(aa=='T'))return(13); // Adding 1H mass
    }else if((aa=='W')||(aa=='F')||(aa=='Y')||(aa=='H')){
      if(atom_name[1]=='G')return(12); // 0H
      if((atom_name[1]=='Z')&&(aa=='Y'))return(12); // 0H
      return(13); // Adding 1H mass
    }else if(atom_name[1]=='G'){
      if((aa=='I')&&(atom_name[2]=='2'))return(15);  // Adding 3H mass
      if(aa=='V')return(15);      // Adding 3H mass
      if(aa=='L')return(13);      // Adding 1H mass
      if((aa=='N')||(aa=='D'))return(12);  // 0H
    }else if(atom_name[1]=='D'){
      if((aa=='L')||(aa=='I'))return(15);    // Adding 3H mass
      if((aa=='Q')||(aa=='E'))return(12);    // 0H
    }else if(atom_name[1]=='E'){
      if(aa=='M')return(15);      // Adding 3H mass
    }else if(atom_name[1]=='Z'){
      if(aa=='R')return(12);      // 0H mass
    }
    return(14);                 // Adding 2H mass
  }

  if(atom_name[0]=='N'){
    if(atom_name[1]==' '){
      if(aa=='P')return(14);           // 0H
      return(15);                      // Adding 1H
    }else if(atom_name[1]=='D'){
      if(aa=='H')return(14);           // 0H
    }else if(atom_name[1]=='E'){
      if((aa=='H')||(aa=='R'))return(15); // Adding 1H
    }else if(atom_name[1]=='Z'){
      if(aa=='K')return(17); // Adding 3H
    }
    return(16);  // Adding 2H mass
  }

  if(atom_name[0]=='O'){
    if(atom_name[1]==' ')return(16); // 0H
    if((atom_name[1]=='G')||(atom_name[1]=='H'))return(17); // 1H
    return(16); // 0H
  }

  if(atom_name[0]=='S')return(32);  // Adding H mass
  if(atom_name[0]=='H')return(1);
  if(atom_name[0]=='P')return(31);
  if(atom_name[0]=='E')return(79); // Selenium
  printf("WARNING, unusual atom %s\n", atom_name);
  return(14);
}


float Mass(atom *atom){
  char *atom_name=atom->name;
  if(atom_name[0]=='C')return(12);
  if(atom_name[0]=='N')return(14); // Adding one hydrogen atom mass
  if(atom_name[0]=='O')return(16);
  if(atom_name[0]=='H')return(1);
  if(atom_name[0]=='P')return(31);
  if(strncmp(atom_name, "SE", 2)==0)return(79); // Selenium
  if(atom_name[0]=='S')return(32);
  printf("WARNING, unusual atom %s\n", atom_name);
  return(12);
}


int Bivariate_fit(int n, float *z, float *x, float *y, int *w,
		  float *a, float *b, float *c, float *r, float *chi)
{
  // fit: z_i = a*x_i + b*y_i + c
  // Returns a, b, c, correlation coefficient and normalized sum of residuals

  int i, norm=0;
  float xy_ave=0, xz_ave=0, yz_ave=0;
  float x2_ave=0, y2_ave=0, z2_ave=0;
  float x_ave=0, y_ave=0, z_ave=0, chi1;

  if(n<3){
    printf("Too few points (%d), the fit can not be performed!\n\n", n);
    return(0);
  }

  for(i=0; i<n; i++){
    if(w[i]){
      x_ave+=x[i]; y_ave+=y[i]; z_ave+=z[i]; norm++;
      x2_ave+=x[i]*x[i]; y2_ave+=y[i]*y[i]; z2_ave+=z[i]*z[i];
      xy_ave+=x[i]*y[i]; xz_ave+=x[i]*z[i]; yz_ave+=y[i]*z[i];
    }
  }

  xy_ave=(float)norm*xy_ave-(x_ave*y_ave);
  xz_ave=(float)norm*xz_ave-(x_ave*z_ave);
  yz_ave=(float)norm*yz_ave-(y_ave*z_ave);
  x2_ave=(float)norm*x2_ave-(x_ave*x_ave);
  y2_ave=(float)norm*y2_ave-(y_ave*y_ave);
  z2_ave=(float)norm*z2_ave-(z_ave*z_ave);

  *a=(y2_ave*xz_ave-xy_ave*yz_ave)/(x2_ave*y2_ave-xy_ave*xy_ave);
  *b=(x2_ave*yz_ave-xy_ave*xz_ave)/(x2_ave*y2_ave-xy_ave*xy_ave);
  *c=z_ave-(*a)*x_ave-(*b)*y_ave; *c/=(float)norm;
  *r=((*a)*xz_ave+(*b)*yz_ave)/
    sqrt(z2_ave*((*a)*(*a)*x2_ave+(*b)*(*b)*y2_ave+ 2*(*a)*(*b)*xy_ave));
  (*chi)=0;
  for(i=0; i<n; i++){
    chi1= z[i]-(*a)*x[i]-(*b)*y[i]-(*c);
    (*chi) += chi1*chi1;
  }
  (*chi)/=z2_ave;
  return(norm);
}


/**************************** B factors *****************************/
void Compute_anisou(float ***aniso_pred, int N_modes, int N_ref,
		    float *one_over_omega2, float **Cart_mode)
{
  // Predicted anisotropic temperature factors
  int ia, i_atom, i, j;
  float *Cmode, Cmode1[3], Cmode2[3];

  for(i_atom=0; i_atom<N_ref; i_atom++){
    int k = 3*i_atom;
    float **anisou=Allocate_mat2_f(3,3);
    aniso_pred[i_atom]=anisou;

    for(ia=0; ia< N_modes; ia++){
      if(one_over_omega2[ia] == 0)continue;
      Cmode=Cart_mode[ia];
      for(i=0; i<3; i++){
	Cmode1[i]=Cmode[k+i];
	Cmode2[i]=Cmode1[i]*one_over_omega2[ia];
	for(j=0; j<=i; j++)anisou[i][j]+=Cmode1[i]*Cmode2[j];
      }
    }
    for(i=0; i<3; i++)for(j=i+1; j<3; j++)anisou[i][j]=anisou[j][i];
  }
  return;
}

void Predict_fluctuations(struct Normal_Mode *NM, float *sigma2)
{
  int naxe=NM->N_axes, ncart=NM->N_Cart, i, a;
  double tors_fluct[naxe];
  for(a=0; a<naxe; a++)tors_fluct[a]=0;
  double cart_fluct[ncart];
  for(i=0; i<ncart; i++)cart_fluct[i]=0;
  for(int ia=0; ia< NM->N; ia++){
    if(NM->sigma2[ia]==0)continue;
    float w=sigma2[ia]; //NM->
    float *u=NM->Tors[ia];
    for(a=0; a<naxe; a++){
      tors_fluct[a]+=w*(*u)*(*u); u++;
    }
    float *x=NM->Cart[ia];
    for(i=0; i<ncart; i++){
      cart_fluct[i]+=w*(*x)*(*x); x++;
    }
  }
  // Copy in struct NM
  for(a=0; a<naxe; a++)NM->tors_fluct[a]=tors_fluct[a];
  for(a=0; a<ncart; a++)NM->cart_fluct[a]=cart_fluct[a];
}

void Print_tors_fluct2(struct Normal_Mode NM, struct axe *axe, int naxe,
		       struct residue *res, char *nameout1)
{
  // Print in file
  char namefile[200]="";
  sprintf(namefile, "%s.MSF_tors.dat", nameout1);
  printf("Writing %s\n", namefile);
  FILE *file_out=fopen(namefile, "w");
  fprintf(file_out, "# dof(type, res, chain, _a.a.)\tfluct\n");
  for(int a=0; a<naxe; a++){
    struct residue *r=res+(axe[a].bond)->atom->res;
    fprintf(file_out, "%c%s%c_%c\t%.3g\n",
	    axe[a].type, r->pdbres, r->chain, r->amm, NM.tors_fluct[a]);
  }
  fclose(file_out);
}


int Set_B_exp(atom *atom, int N_atoms, int N_res, float *B_exp, char *REF){
  int i, nca=0;
  float pi=3.1415, norm=8.0*pi*pi/3.0;
  for(i=0; i<N_atoms; i++){
    if(Selection(atom->name, atom, REF)){
      B_exp[nca]=atom->B_factor/norm; nca++;
      if(nca==N_res)break;
    }
    atom++;
  }
  printf("B factors: %d %s atoms selected\n", nca, REF);
  return(nca);
}

void Set_anisou(atom *atom, int N_atoms, int nca,
		float ***anisou, float ***aniso_pred, char *REF)
{
  int i_atom, i, j, ica=0;

  for(i_atom=0; i_atom<N_atoms; i_atom++){
    if(Selection(atom->name, atom, REF)){
      anisou[ica]=malloc(3*sizeof(float *));
      aniso_pred[ica]=malloc(3*sizeof(float *));
      for(i=0; i<3; i++){
	aniso_pred[ica][i]=malloc(3*sizeof(float));
	anisou[ica][i]=malloc(3*sizeof(float));
	for(j=0; j<3; j++)anisou[ica][i][j]=atom->anisou[i][j];
      }
      ica++; if(ica==nca)break;
    }
    atom++;
  }
  return;
}

int Compare_anisou(float ***aniso_pred, float ***aniso_exp, int nca,
		   char *model, char *inter, char *file_name, char *prot_name,
		   int N, int *weight, float *dot_aniso, float *ov_aniso,
		   float *delta_aniso)
{
  //FILE *file_out=fopen(file_name, "a");
  int ica, i, j, n_sum=0;
  float sum_overlap=0, sum_dot=0, sum_delta=0;
  float delta, overlap, dot, dd;
  float **e_vec1=Allocate_mat2_f(3,3), e_val1[3];
  float **e_vec2=Allocate_mat2_f(3,3), e_val2[3];
  float **e_vec12=Allocate_mat2_f(3,3), e_val12[3];
  float **u_12=Allocate_mat2_f(3,3);

  for(ica=0; ica<nca; ica++){
    if(weight[ica]){
      float **u1=aniso_exp[ica],  norm1=0;
      float **u2=aniso_pred[ica], norm2=0;
      float prod_lambda=1.0, prod_lambda_12=1.0;
      float isotropy;

      // Exclude quasi-spherical distributions, isotropy >= 0.5
      f_Diagonalize(3, u1, e_val1, e_vec1, 1);
      isotropy=e_val1[2]/e_val1[0]; //printf("%.2f\n", anisotropy);
      if(isotropy >= 0.5)continue;
      n_sum++;

      /****************  Overlap *********************/
      // Normalize by the trace
      for(i=0; i<3; i++)norm1+=u1[i][i];
      for(i=0; i<3; i++)norm2+=u2[i][i];
      for(i=0; i<3; i++)for(j=i; j<3; j++){
	  u1[i][j]/=norm1; u2[i][j]/=norm2;
	  u_12[i][j]=u1[i][j]+u2[i][j];
	}
      for(i=0; i<3; i++)for(j=0; j<i; j++){
	  u1[i][j]=u1[j][i]; u2[i][j]=u2[j][i]; u_12[i][j]=u_12[j][i];
	}

      for(i=0; i<3; i++)e_val1[i]/=norm1;
      f_Diagonalize(3, u2, e_val2, e_vec2, 1);
      f_Diagonalize(3, u_12, e_val12, e_vec12, 1);
      for(i=0; i<3; i++){
	prod_lambda*=(e_val1[i]*e_val2[i]);
	prod_lambda_12*=e_val12[i];
      }
      overlap=sqrt(8.*sqrt(prod_lambda)/prod_lambda_12);
      sum_overlap += overlap;
      // printf("overlap: %.3f\n", overlap);

      /*****   dot product   **************/
      /*
	dot=0; dot_norm=0;
	for(i=0; i<3; i++){
	float w=e_val1[i]*e_val2[i];
	dot+=w*Scalar_product_3(e_vec1[i], e_vec2[i]);
	dot_norm+=w;
	}
	dot /= dot_norm;
      */
      dot=Scalar_product_3(e_vec1[0], e_vec2[0]);
      if(dot<0){dot=-dot;} sum_dot += dot;

      /***********   delta   **************/
      delta=0;           // delta = sqrt(sum_ij(U1_ij-U2_ij)^2/9)
      for(i=0; i<3; i++){
	dd=u1[i][i]-u2[i][i]; delta+=dd*dd;
	for(j=0; j<i; j++){dd=u1[i][j]-u2[i][j]; delta+=2*dd*dd;}
      }
      sum_delta+=delta;

      /*
      //fprintf(file_out, "%6.3f %6.3f\n", overlap, dot);
      for(i=0; i<3; i++)for(j=0; j<=i; j++){
      fprintf(file_out, "%6.3f %6.3f ", u1[i][j], u2[i][j]);
      }
      fprintf(file_out, "\n");
      */
    }
  }
  *ov_aniso=sum_overlap/n_sum;
  *dot_aniso=sum_dot/n_sum;
  *delta_aniso=sqrt(sum_delta/(9*n_sum));

  /*
    fprintf(file_out, "%6.3f  %s %8s %6.3f %2.0f %8s %3d\n",
	  sum_overlap/n_sum, model, inter, sum_dot/n_sum, n_sum, prot_name, N);
  fclose(file_out);
  */

  // Clean
  Empty_matrix_f(e_vec1, 3);
  Empty_matrix_f(e_vec2, 3);
  Empty_matrix_f(e_vec12, 3);
  Empty_matrix_f(u_12, 3);

  return(n_sum);
}

void Rescale(float *x, int N, float slope, float offset)
{
  int i; for(i=0; i<N; i++)x[i]=x[i]*slope+offset;
}



/***************************   Eigensystem  *********************************/


void f_sort(float d[], int n, int *i_rank)
{
  // Sort the vector d from large to small.
  // Returns i_rank[i]= index of object at rank i
  int k,j,i, jmax, *not_ranked=malloc(n*sizeof(int));
  float d_max;

  for (i=0; i<n; i++)not_ranked[i]=1;
  for (k=0; k<n; k++) {
    for(i=0; i<n; i++)if(not_ranked[i])break;
    d_max=d[jmax=i];
    for (j=i+1;j<n; j++){
      if(not_ranked[j]&&(d[j] >= d_max))d_max=d[jmax=j];
    }
    i_rank[k]=jmax; not_ranked[jmax]=0;
  }
  free(not_ranked);
}


/********************************  Output **********************************/
float Output_modes(int N,           //degrees of freedom
		   char *model,     // INM, ANM or TNM
		   char *REF_ATM,   // Reference atoms
		   atom *atoms, int N_atoms,
		   struct residue *seq, int N_res,
		   struct axe *axes, int N_axes,
		   float *B_CA, float *B_CA_exp,
		   int *weight, float *dist_CA,
		   float ***aniso_pred, float ***aniso_exp,
		   int N_MODES, float cutoff,
		   char *inter, char *prot_name,
		   FILE *file_B, FILE *file_sum,
		   float **eigen_vector, float *eigen_value, float *eigen_B,
		   float *Cart_collectivity, float *Tors_collectivity,
		   float **Cart_mode, int *atom_ref, int N_ref)
{
  char outfile[200], string[100]; int n_aniso=0;
  float r, slope, offset, dot_aniso=0, ov_aniso=0, delta_aniso=2.;
  float slope_calc, slope_dist, offset2, r2, chi;
  int neg=0, i;
  int nca=N_res;

  // B factors comparison
  // Model 1: B_exp = slope*B_calc + offset
  // r=Corr_coeff(B_CA, B_CA_exp, weight, nca, &slope, &offset, &newnorm);
  r=Corr_coeff(B_CA, B_CA_exp, nca, &slope, &offset);
  printf("cc(B, %s %s)=  %6.3f   ref: %s\n", model, inter, r, REF_ATM);

  fprintf(file_sum, "%3s\t%3s\t%4.1f", model, inter, cutoff);
  fprintf(file_sum, "\t%6.3f\t%8.5f\t%5.3f", r, slope, offset);
  if(FIT_ROT){
    // Model 2: B_exp = slope*B_calc + offset
    Bivariate_fit(nca, B_CA_exp, B_CA, dist_CA, weight,
		  &slope_calc, &slope_dist, &offset2, &r2, &chi);
    fprintf(file_sum, "\t%6.3f\t%8.5f\t%8.5f\t%5.3f ",
	    r2, slope_calc, slope_dist, offset2);
  }
  //fprintf(file_sum, "%8s  %3d  %4.1f", prot_name, N, cutoff);
  // Anisotropic B factors
  if(ANISOU){
    if(strncmp(model, "INM", 3)!=0){
      n_aniso=Compare_anisou(aniso_pred, aniso_exp, nca, model, inter,
			     file_aniso, prot_name, N, weight,
			     &dot_aniso, &ov_aniso, &delta_aniso);
    }
    fprintf(file_sum, "\t%3d\t%.3f\t%6.3f\t%.4f",
	    n_aniso, ov_aniso, dot_aniso, delta_aniso);
  }
  fprintf(file_sum, "\t%s\n", REF_ATM);
  fflush(file_sum);

  // Print isotropic B factors
  sprintf(string, "%s_%s", model, inter);
  Rescale(B_CA, nca, slope, offset);
  Print_B(B_CA, nca, file_B, string, r);

  // Normalize eigenvalues
  if(slope>0)for(i=0; i<N; i++)eigen_value[i]/=slope;
  // for(i=0; i<N; i++)eigen_value[i]*=100;
  // Print negative eigenvalues, if any
  for(i=0; i<N; i++){
    if(eigen_value[i]<-0.001){
      if(neg==0){printf("Negative eigenvalues:\n"); neg=1;}
      printf("%.4f\n", eigen_value[i]);
    }
  }

  // Print normal modes
  sprintf(outfile, "%s_%s_%s", prot_name, inter, model);
  Print_modes_old(eigen_value, eigen_vector, eigen_B,
		  Cart_collectivity, Tors_collectivity, N, N_MODES,
		  atoms, N_atoms, nca, seq, axes, N_axes, outfile, model,
		  Cart_mode, atom_ref, N_ref);

  return(r);

}


void Compare_modes(float **Cart_mode_ANM, float **Cart_mode_TNM, int N_Cart,
		   float **Tors_mode_ANM, float **Tors_mode_TNM, int N_tors,
		   float *mass_coord, char *prot_name, char *inter)
{
  // Compare torsional and cartesian normal modes
  int ia, i;
  float q;
  float *Ovlp_cart_max_TNM=malloc(N_tors*sizeof(float));
  float *Ovlp_tors_max_TNM=malloc(N_tors*sizeof(float));
  float *Ovlp_cart_max_ANM=malloc(N_Cart*sizeof(float));
  float *Ovlp_tors_max_ANM=malloc(N_Cart*sizeof(float));
  int *mode_cart_max_TNM=malloc(N_tors*sizeof(int));
  int *mode_tors_max_TNM=malloc(N_tors*sizeof(int));
  int *mode_cart_max_ANM=malloc(N_Cart*sizeof(int));
  int *mode_tors_max_ANM=malloc(N_Cart*sizeof(int));

  FILE *file_out; char outfile[200];

  sprintf(outfile, "%s_%s_compare.dat", prot_name, inter);
  file_out=fopen(outfile, "w");
  printf("Writing %s\n", outfile);

  for(i=0; i<N_tors; i++){
    Ovlp_cart_max_TNM[i]=0;
    Ovlp_tors_max_TNM[i]=0;
  }
  for(i=0; i<N_Cart; i++){
    Ovlp_cart_max_ANM[i]=0;
    Ovlp_tors_max_ANM[i]=0;
  }
  // Normalize
  for(i=0; i<N_tors; i++){
    Normalize_vector_weighted(Cart_mode_TNM[i], mass_coord, N_Cart);
    Normalize_vector(Tors_mode_TNM[i], N_tors);
  }
  for(i=0; i<N_Cart; i++){
    Normalize_vector_weighted(Cart_mode_ANM[i], mass_coord, N_Cart);
    Normalize_vector(Tors_mode_ANM[i], N_tors);
  }

  for(ia=0; ia<N_tors; ia++){
    float *Ca=Cart_mode_TNM[ia];
    for(i=0; i<N_Cart; i++){
      q=Scalar_product_weighted(Ca, Cart_mode_ANM[i], mass_coord, N_Cart);
      q*=q;
      if(q > Ovlp_cart_max_TNM[ia]){
	Ovlp_cart_max_TNM[ia]=q; mode_cart_max_TNM[ia]=i;
      }
      if(q > Ovlp_cart_max_ANM[i] ){
	Ovlp_cart_max_ANM[i]=q; mode_cart_max_ANM[i]=ia;
      }
    }
  }
  for(ia=0; ia<N_tors; ia++){
    float *Ta=Tors_mode_TNM[ia];
    for(i=0; i<N_Cart; i++){
      q=Scalar_product(Ta, Tors_mode_ANM[i],N_tors); q*=q;
      if(q > Ovlp_tors_max_TNM[ia]){
	Ovlp_tors_max_TNM[ia]=q; mode_tors_max_TNM[ia]=i;
      }
      if(q > Ovlp_tors_max_ANM[i] ){
	Ovlp_tors_max_ANM[i]=q; mode_tors_max_ANM[i]=ia;
      }
    }
  }

  fprintf(file_out, "# Max. coeff Cartesian TNM - ANM modes\n");
  fprintf(file_out, "# TNM_mode  Max.coeff^2 ANM_mode\n");
  for(i=0; i<N_tors; i++)
    fprintf(file_out, "%3d  %.3f %3d\n",
	    i, Ovlp_cart_max_TNM[i], mode_cart_max_TNM[i]);
  fprintf(file_out, "&\n");

  fprintf(file_out, "# Max. coeff Torsional TNM - ANM modes\n");
  fprintf(file_out, "# TNM_mode  Max.coeff^2 ANM_mode\n");
  for(i=0; i<N_tors; i++)
    fprintf(file_out, "%3d\t%.3f\t%3d\n",
	    i, Ovlp_tors_max_TNM[i], mode_tors_max_TNM[i]);
  fprintf(file_out, "&\n");

  fprintf(file_out, "# Max. coeff Cartesian ANM - TNM modes\n");
  fprintf(file_out, "# ANM_mode  Max.coeff^2  TNM_mode\n");
  for(i=0; i<N_Cart; i++)
    fprintf(file_out, "%3d\t%.3f\t%3d\n",
	    i, Ovlp_cart_max_ANM[i], mode_cart_max_ANM[i]);
  fprintf(file_out, "&\n");

  fprintf(file_out, "# Max. coeff Torsional ANM - TNM modes\n");
  fprintf(file_out, "# ANM_mode  Max.coeff^2  TNM_mode\n");
  for(i=0; i<N_Cart; i++)
    fprintf(file_out, "%3d\t%.3f\t%3d\n",
	    i, Ovlp_tors_max_ANM[i], mode_tors_max_ANM[i]);
  fprintf(file_out, "&\n");

  fclose(file_out);
  free(Ovlp_cart_max_TNM); free(mode_cart_max_TNM);
  free(Ovlp_tors_max_TNM); free(mode_tors_max_TNM);
  free(Ovlp_cart_max_ANM); free(mode_cart_max_ANM);
  free(Ovlp_tors_max_ANM); free(mode_tors_max_ANM);

}

void Print_B_fact(float *B_TNM, float *B_pred_all, float *B_exp,
		  int N, atom *atoms, int *atom_ref, struct residue *seq,
		  char *name, char *what, float cc, float slope,
		  float *dof, char ridge)
{
  char nameout[200];  int i;
  atom *atom1=NULL;
  FILE *file_out;

  //sprintf(nameout, "%s_%s.dat", name, what);
  sprintf(nameout, "%s.MSF.dat", name);
  file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  if(B_pred_all && B_exp) {
    fprintf(file_out, "# Ridge regression of type %c\n", ridge);
    fprintf(file_out, "# slope=%.4g corr.=%.3f (%s)\n", slope, cc, what);
    fprintf(file_out, "# Internal_fluct= %.3f", dof[0]);
    fprintf(file_out, " Translation_fluct= %.3f", dof[1]);
    fprintf(file_out, " Rotation_fluct= %.3f\n", dof[2]);
    fprintf(file_out, "# aa_pos_chain pred_msd_int pred_msd_all B/(8pi^2/3)\n");
  }else{
    fprintf(file_out, "# aa_pos_chain pred_msd_internal\n");
  }
  
  for(i=0; i<N; i++){atom1=atoms+atom_ref[i];
    struct residue *r=seq+atom1->res;
    int j;
    for(j=(int)sizeof(r->pdbres)-1; j>=0; j--){
      if(r->pdbres[j]==' '){r->pdbres[j]='\0';}
      else{break;}
    }
    for(j=0; j<(int)sizeof(r->pdbres); j++)
    	if(r->pdbres[j]!=' ')break;
    char name[20];
    sprintf(name, "%c%s%c", r->amm, r->pdbres+j, r->chain);
    fprintf(file_out, "%s\t%.4g", name, B_TNM[i]);
    if(B_pred_all && B_exp)
      fprintf(file_out, "\t%7.3f\t%7.3f", B_pred_all[i], B_exp[i]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}

void Print_B(float *B, int nca, FILE *file_out, char *model, float cc)
{
  int i;
  fprintf(file_out, "# B factors: %s, corr.=%.3f\n", model, cc);
  for(i=0; i<nca; i++)fprintf(file_out, "%3d\t%7.3f\n", i, B[i]);
  fprintf(file_out, "&\n"); fflush(file_out);
}

void Print_PDB_3(float *Cart_mode, float eigen_value, float eigen_B,
		 atom *atoms, int N_atom_ref, int *atom_ref,
		 struct residue *seq, char *file_name, int ia)
{
  // Print PDB file
  FILE *file_out;  char outfile[200];
  int i, j, k; atom *atom1=atoms;
  char aaname3[10];
  // Amplitude factors
  float DTHETA=0;
  float AMPL_FACTOR=16; // Amplification with respect to thermal fluctuations
  //~ float AMAX=120;      // New amplitude factor
  float step, dtheta, r[3];
  int N_STEP=10; // Number of conformations per normal mode
  int i_step;

  for(i=0; i<10; i++)aaname3[i]='\0';

  // Opening file
  sprintf(outfile, "%s_%d.pdb", file_name, ia+1);
  file_out=fopen(outfile, "w");
  printf("Writing normal modes in PDB format in %s\n", outfile);
  fprintf(file_out, "REMARK  Normal mode %3d  Percent fluctuation= %.2f",
	  ia, eigen_B*100.0);

  // Compute amplitude
  DTHETA=1./sqrt(eigen_value);
  fprintf(file_out, " sqrt<theta^2>= %.4f\n", DTHETA);
  DTHETA*=AMPL_FACTOR;
  step=DTHETA/N_STEP; dtheta=0;

  for(i_step=0; i_step<N_STEP; i_step++){
    if(i_step==0){
      fprintf(file_out, "MODEL: ground state\n");
    }else{
      fprintf(file_out, "MODEL %d\n", i_step);
    }

    // Compute and print
    int jatom=0;
    for(i=0; i<N_atom_ref; i++){
      int m=3*i;
      atom1=atoms+atom_ref[i];
      for(j=0; j<3; j++)r[j]=atom1->r[j]+Cart_mode[m+j]*dtheta;
      k=atom1->res; Name3(aaname3, seq[k].i_aa);
      if(jatom <99999)jatom++;
      fprintf(file_out,
	      "ATOM  %5d%4s  %3s %c%4s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	      jatom, atom1->name, aaname3, chain,
	      seq[k].pdbres,  r[0], r[1], r[2], 1.0, atom1->B_factor);
      atom1=atoms+atom1->i_next;
    }
    fprintf(file_out, "ENDMDL\n");
    dtheta+=step;
  }
  fclose(file_out);
}



void Print_modes_old(float *eigen_value, float **eigen_vector, float *eigen_B,
		     float *Cart_collectivity, float *Tors_collectivity,
		     int N, int N_MODES, atom *atoms, int N_atoms, int nca,
		     struct residue *seq, struct axe *axes, int N_axes,
		     char *file_name, char *deg_of_freedom,
		     float **Cart_mode, int *atom_ref, int N_ref)
{
  // Notice that eigenvectors are defined between 0 and N-1 !!

  FILE *file_out; int ia, i, n_max=N_MODES, tnm=0;
  char outfile[200];

  if(PRINT_LAMBDA){
    double sum=0; float kappa;
    sprintf(outfile, "%s_lambda.dat", file_name);
    file_out=fopen(outfile, "w");
    printf("Writing %s\n", outfile);
    kappa=Collectivity_norm1(eigen_B, N);
    fprintf(file_out, "# Reciprocal collectivity of fluctuations=  %.1f\n",
	    kappa);
    fprintf(file_out, "# mode Contr_to_fluct Cumulative lambda Coll.(cart)");
    if(Tors_collectivity)fprintf(file_out, " Coll.(tors)");
    fprintf(file_out, "\n");

    for(i=0; i<N; i++){
      if(eigen_value[i] < E_MIN)fprintf(file_out, "#");
      sum+=eigen_B[i];
      fprintf(file_out, "%3d\t%.4f\t%.3f\t%7.4g\t%.3f",
	      i, eigen_B[i], sum, eigen_value[i], Cart_collectivity[i]);
      if(Tors_collectivity)fprintf(file_out, "\t%.3f", Tors_collectivity[i]);
      fprintf(file_out, "\n");
    }
    fclose(file_out);
  }


  if(n_max==0)return;
  if(n_max > N)n_max=N;

  sprintf(outfile, "%s.dat", file_name);
  file_out=fopen(outfile, "w");
  printf("Writing normal modes to file %s\n", outfile);
  fprintf(file_out, "# %d degrees of freedom: %s\n", N, deg_of_freedom);
  fprintf(file_out, "# Contribution to fluctuation:");
  for(ia=0; ia<n_max; ia++)fprintf(file_out, "\t%.3g", eigen_B[ia]);
  fprintf(file_out, "\n");
  fprintf(file_out, "# e_values:                    ");
  for(ia=0; ia<n_max; ia++)fprintf(file_out, "\t%.3g", eigen_value[ia]);
  fprintf(file_out, "\n");

  if(strncmp(deg_of_freedom, "TNM", 3)==0)tnm=1;
  for(i=0; i<N; i++){
    fprintf(file_out, "%4d", i);
    for(ia=0; ia<n_max; ia++)fprintf(file_out, "\t%.6g", eigen_vector[ia][i]);
    if(tnm){
      fprintf(file_out, "\t%c%s", axes[i].type,
	      seq[(axes[i].bond)->atom->res].pdbres);
    }
    fprintf(file_out, "\n");
  }

  /* // Check global rotation eigenvector (only ANM)
  if(N==nca*3){
    float *ev, *r, sum; int ica;
    for(n=1; n<=n_max; n++){
      // scalar product sum_i vi*ri
      sum=0; ica=0; ev=eigen_vector[N-n];
      for(i=0; i<N_atoms; i++){
	if(strncmp(atoms[i].name, "CA", 2)!=0)continue; r=atoms[i].r;
	sum += ev[ica]*r[0]+ev[ica+1]*r[1]+ev[ica+2]*r[2];
	ica+=3;
      }
      fprintf(file_out, "# ev %d, lambda= %.4f sum(ri*vi)= %6.3f\n",
	      n, eigen_value[N-n], sum);
    }
    }*/

  fclose(file_out);

  if(tnm){
    for(ia=0; ia<n_max; ia++){
      if(eigen_value[ia] < E_MIN)continue;
      Print_PDB_3(Cart_mode[ia], eigen_value[ia], eigen_B[ia], atoms,
		  N_ref, atom_ref, seq, file_name, ia);
    }
  }
}

void Compute_force(struct Tors *Force, float *cc, struct Normal_Mode NM)
{
  // Compute force from projection of conformation change on normal modes
  int i, a;
  Force->Cart=malloc(NM.N_Cart*sizeof(float));
  for(i=0; i<NM.N_Cart; i++)Force->Cart[i]=0;
  Force->Tors=malloc(NM.N_axes*sizeof(float));
  Force->MW_Tors=malloc(NM.N_axes*sizeof(float));
  for(i=0; i<NM.N_axes; i++){Force->Tors[i]=0; Force->MW_Tors[i]=0;}
  Force->coeff=malloc(NM.N*sizeof(float));
  for(a=0; a<NM.N_relevant; a++){
    if((NM.select[a]==0)||(NM.sigma2[a]==0))continue;
    float c=NM.omega2[a]*cc[a];
    float *x=NM.Cart[a], *u=NM.Tors[a], *w=NM.MW_Tors[a];
    for(i=0; i<NM.N_Cart; i++){Force->Cart[i]+=c*(*x); x++;}
    for(i=0; i<NM.N_axes; i++){
      Force->Tors[i]+=c*(*u); u++;
      Force->MW_Tors[i]+=c*(*w); w++;
    }
  }
}

void Cartesian_force(float *Cart_force, float *atom_diff,
		     double **Hessian,int N_Cart)
{
  int i, j;
  for(i=0; i<N_Cart; i++){
    double sum=0;
    for(j=0; j<N_Cart; j++)sum+=Hessian[i][j]*atom_diff[j];
    Cart_force[i]=sum;
  }
}

void Torsional_force(struct Tors Force, struct Tors *Diff, double **Hessian,
		     int N_axes, int N_ref, int *atom_num,
		     float **Jacobian_ar, int N_int,
		     struct interaction *Int_list, atom *atoms, int natoms)
{
  int N_Cart=3*N_ref, i, j, a; double sum;
  for(i=0; i<N_axes; i++){
    sum=0;
    for(j=0; j<N_axes; j++)sum += Hessian[i][j]*Diff->MW_Tors[j];
    Force.Tors[i]=sum;
  }
  // Cartesian force fr= HrJDelta_theta
  float dr[N_Cart];
  for(i=0; i<N_Cart; i++){
    sum=0; for(a=0; a<N_axes; a++)sum+=Jacobian_ar[a][i]*Diff->Tors[a];
    dr[i]=sum;
  }
  double **H=Allocate_mat2_d(N_Cart, N_Cart);
  Compute_Hessian_ANM(H,N_ref,atom_num,Int_list,N_int,atoms,natoms);
  for(i=0; i<N_Cart; i++){
    sum=0; for(j=0; j<N_Cart; j++)sum += H[i][j]*dr[j];
    Force.Cart[i]=sum;
  }
  Empty_matrix_d(H, N_Cart);
}

float Compute_Max_dev(float *Cart_mode, float omega2, int N_Cart,
		      float *mass_sqrt)
{
  float Max_dev=0; int i;
  for(i=0; i< N_Cart; i++){
    float x=fabs(Cart_mode[i])*mass_sqrt[i]; //;
    if(x>Max_dev)Max_dev=x;
  }
  //return(3*Max_dev/sqrt(omega2));
  return(Max_dev/sqrt(omega2));
}

void Axes_path(int *main1, int *main2, int *side1, int *side2,
	       struct axe *axes, int i, int nmain, int naxes)
{
  *main1=-1; *main2=-1; *side1=-1; *side2=-1;
  if(i < nmain){
    *main1=axes[i].first_main;
    *main2=i;
  }else{
    *main1=axes[i].first_main;
    *main2=axes[i].last_main;
    *side1=axes[i].first_side;
    *side2=i;
  }
  if((*main2 >= nmain) ||
     (*main2 < *main1)||
     (*side2 < *side1)||
     (*side2>=naxes)){
    printf("ERROR, wrong axes path\n");
    printf("i= %d res= %d main1= %d res=%d main2= %d res=%d",
	   i, (axes[i].bond)->atom->res,
	   *main1, (axes[*main1].bond)->atom->res,
	   *main2, (axes[*main2].bond)->atom->res);
    if(*side1>=0)
      printf(" side1= %d res=%d side2= %d res=%d\n",
	     *side1, (axes+*side1)->bond->atom->res,
	     *side2, (axes[*side2].bond)->atom->res);
    exit(8);
  }
}

void Tors_outliers(int *outlier_tors, int naxe, int *outlier, int Na,
		   struct axe *axe, struct chain *chains, int Nchain)
{
  int a=0, a_ini=0;
  for(a=0; a<naxe; a++)outlier_tors[a]=0;
  for(int c=0; c<Nchain; c++){
    struct chain *ch=chains+c;
    int i=ch->ini_res, end_res=ch->ini_res+ch->nres-1;
    while(1){
      if((outlier[i]==0)&&(outlier[i+1]==0))break;
      while(1){
	int a=Get_tors(axe, naxe, i, a_ini);
	if(a>=0){outlier_tors[a]=1; a_ini=a+1;}
	else{break;}
      }
      i++; if(i>=end_res)break;
    }
    i=end_res;
    while(1)
    	{if((outlier[i]==0)&&(outlier[i-1]==0))break;
      	 while(1){
			int a=Get_tors(axe, naxe, i, a_ini);
			if(a>=0){outlier_tors[a]=1; a_ini=a+1;}
			else{break;}
      		}
      	i--;
      	if(i<=ch->ini_res)break;
    	}
  }
}

int Get_tors(struct axe *axe, int naxe, int res, int a_ini){
  for(int a=a_ini; a<naxe; a++){
    if((axe[a].bond)->atom->res == res){
      return(a);
    }else if((axe[a].bond)->atom->res > res){
      return(-1);
    }
  }
  return(-1);
}

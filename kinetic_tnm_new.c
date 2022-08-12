int DBG=0;
#define T_COMPUTE 0 // 0=Cholevsky 1=Diagonalization
#include "nma_para.h"
#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "vector.h"
#include "allocate.h"
#include "kinetic_tnm.h"
#include "buildup.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// Rigid body degrees of freedom
void Center_atoms(atom *atoms, int N_atoms, struct Reference Ref);
void Reorient_axes(atom *atoms, int N_atoms, struct Reference Ref);
void Eckart_main(struct axe *axes, int N_axes, atom *atoms, int N_atoms,
		 struct Reference Ref);
void Eckart_test(float **Jacobian_ar, int N_axes,
		 atom *atoms, int N_atoms, struct Reference Ref);
int Test_Eckart(float *Cart_dev, atom *atoms, int N_atoms,
		 struct Reference Ref, int mode);

void Inertia_tensor(double **inertia, double **corr_sum, double *r_sum,
		    double M_sum);
float Sum_inertia(double *mr_tot, double **corr_sum,
		  atom *atoms, struct Reference Ref,
		  int ini_ref, int end_ref);
float Empty_inertia(double *mr_tot, double **corr_sum);

int Kinetic_energy(double ***T_mat, float *mass_coord,
		   struct axe *axe, int ini_axe, int N_axes,
		   atom *atoms, int *atom_num, int N_ref);
int Kinetic_energy_Jac(struct Jacobian J, struct Reference Ref, int N_axes);
int Kinetic_energy_fast(struct Jacobian *J, struct axe *axe, int N_axes,
			atom *atoms, int N_atoms, struct Reference Ref);
int Kinetic_sqrt(struct Jacobian *J, int *N_kin,
		 int N, int CONTROL, float E_MIN);
int Internal_Jacobian(// Output:
		      struct Jacobian J,
		      // Input:
		      struct Reference Ref,
		      atom *atoms, int N_atoms,
		      struct axe *axe, int N_axes);
int Correlate_Jacobian(struct Jacobian J,
		       struct Reference Ref,
		       struct axe *axe, int N_axes);
int Compute_Jtilde(struct Jacobian *J);
int Compute_Rot_Shift(struct Jacobian *J, struct axe *axe);
float **Transpose_matrix(float **M, int n1, int n2);
float **Transpose_matrix_d2f(double **M, int n1, int n2);

// Other
extern void f_Diagonalize(int N, float **MATRIX, float *eigen_values,
			  float **eigen_vector, int SIGN);
extern void d_Diagonalize(int N, double **MATRIX, float *eigen_values,
			  float **eigen_vector, int SIGN);
extern void dd_Diagonalize(int N, double **MATRIX, double *eigen_values,
			   double **eigen_vector, int SIGN);
extern int choldc(double **L, double **a, int N);
extern int choldc_f(float **L, float **a, int N);
extern float **Cholevsky_inversion_d2f(double **L, int N);
extern void Forward_substitution(double *X, double **L, double *Y, int N);
extern void Forward_substitution_f(double *X, float **L, double *Y, int N);
extern void Backward_substitution(double *X, double **L, double *Y, int N);
void Normal_vector(double *v, struct bond *bond, double *vers);
int Dof_overlap(struct axe *axe_k, struct axe *axe_l);
static float Module(double *v);

int Compute_kinetic(struct Jacobian *J,
		    struct axe *axe1, int naxe1,
		    atom *atoms1, int natoms1,
		    struct Reference Ref)
{
  int N_kin=0;
  int test=(J->Jacobian_ar[0][0]==0);
  // Put center of mass at the origin
  Center_atoms(atoms1, natoms1, Ref);
  // Reorient the Cartesian axes along the principal axes
  //Reorient_axes(atoms1, natoms1, Ref);

  Set_rot_shift(axe1, naxe1, atoms1, natoms1);
  Eckart_main(axe1, naxe1, atoms1, natoms1, Ref);
  Internal_Jacobian(*J, Ref, atoms1, natoms1, axe1, naxe1);
  if(test){
    Eckart_test(J->Jacobian_ar, naxe1, atoms1, natoms1, Ref);
    Correlate_Jacobian(*J, Ref, axe1, naxe1);
  }
  if(J->T_sqrt_inv){
    Empty_matrix_f(J->T_sqrt_inv, naxe1); J->T_sqrt_inv=NULL;
  }

  Kinetic_energy_fast(J, axe1, naxe1, atoms1, natoms1, Ref);
  //Kinetic_energy_Jac(*J, Ref, naxe1);

  Kinetic_sqrt(J, &(N_kin), naxe1, T_COMPUTE, E_MIN);

  // Transpose matrix to accelerate computation
  if(J->T_sqrt_inv_tr!=NULL){
    Empty_matrix_f(J->T_sqrt_inv_tr, naxe1); J->T_sqrt_inv_tr=NULL;
  }
  if(J->T_sqrt_tr!=NULL){
    Empty_matrix_f(J->T_sqrt_tr, N_kin); J->T_sqrt_tr=NULL;
  }
  if(J->T_sqrt_inv){
    J->T_sqrt_inv_tr=Transpose_matrix(J->T_sqrt_inv, N_kin, naxe1);
  }else{
    J->T_sqrt_tr=Transpose_matrix_d2f(J->T_sqrt, N_kin, N_kin);
  }
  J->N_kin=N_kin;

  // Empty matrices
  if(J->Jtilde_ar!=NULL){
    Empty_matrix_f(J->Jtilde_ar, naxe1); J->Jtilde_ar=NULL;
  }

  // Show results
  if(DBG){
    int i, j, m=10;
    printf("Masses: ");
    for(i=0; i<3*m; i+=3)printf("%.0f ", Ref.mass_coord[i]);
    printf("\nJacobian: ");
    for(i=0; i<m; i++)printf("%.1f ", J->Jacobian_ar[i][0]);
    printf("\nKinetic energy matrix. ");
    printf("%d dofs atoms= %d cart= %d Nkin= %d\n",
	   naxe1, natoms1, Ref.N_cart, N_kin);
    for(i=0; i<m; i++){
      for(j=0; j<m; j++)printf("%.2g\t", J->T_sqrt[i][j]);
      printf("\n");
    }
    printf("Axis properties\n");
    struct axe *axe=axe1;
    for(i=0; i<30; i++){
      printf("%c%d\t%s-%s\t%d-%d\t%d-%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
	     axe->type, i,
	     axe->bond->previous->atom->name, axe->bond->atom->name,
	     axe->first_atom, axe->last_atom,
	     axe->first_kin, axe->last_kin,
	     axe->rot[0],
	     axe->local_shift[0], axe->local_rot[0],
	     axe->global_shift[0], axe->global_rot[0]);
      axe++;
    }
  }
  //exit(8);

  return(N_kin);
}

int Kinetic_energy(double ***T_mat, float *mass_coord,
		   struct axe *axe, int ini_axe, int N_axes,
		   atom *atoms, int *atom_num, int N_ref)
{
  // Compute kinetic energy tensor T=T_sqrt T_sqrt^T (only internal dof)
  // Lower tridiagonal matrix T_sqrt[i][j] with j<=i
  int k, l, i, j, ini_ref, N_cart=3*N_ref;
  float **J_ar=Allocate_mat2_f(N_axes, N_cart), *J, dr[3];
  *T_mat=Allocate_mat2_d(N_axes, N_axes);

  // Local Jacobian (only atoms upstream of k)
  for(l=0; l<N_axes; l++){
    struct axe *l_axe=axe+ini_axe+l;
    float offset[3], rot[3];
    for(j=0; j<3; j++){
      offset[j]=l_axe->offset[j]; rot[j]=l_axe->rot[j];
    }
    ini_ref=l_axe->first_kin;
    J=J_ar[l]+3*ini_ref;
    for(i=ini_ref; i<N_ref; i++){
      for(j=0; j<3; j++)dr[j]=atoms[atom_num[i]].r[j]-offset[j];
      Vector_product(J, rot, dr); J+=3;
    }
  }

  for(l=0; l<N_axes; l++){
    float *Jl=J_ar[l];
    //~ struct axe *l_axe=axe+ini_axe+l;
    int end_cart=3*(axe[ini_axe+l].first_kin);

    for(k=0; k<=l; k++){
      double sum=0; float *Jk=J_ar[k];
      int ini_cart=3*(axe[ini_axe+k].first_kin);
      for(j=ini_cart; j<end_cart; j++)sum+=Jk[j]*Jl[j]*mass_coord[j];
      (*T_mat)[k][l]=sum;
    }
  }
  for(l=0; l<N_axes; l++){
    for(k=l+1; k<N_axes; k++)(*T_mat)[k][l]=(*T_mat)[l][k];
  }
  if(DEBUG){
    for(l=0; l<10; l++)
      for(k=0; k<=l; k++)
	printf("%d %d %d %d\n", l, k,
	       axe[ini_axe+l].first_kin, axe[ini_axe+k].first_kin);
      printf("\n");
    for(l=0; l<20; l++){
      for(k=0; k<=l; k++)printf("%.2g ", (*T_mat)[l][k]);
      printf("\n");
    }
  }

  Empty_matrix_f(J_ar, N_axes);
  return(0);
}

int Kinetic_energy_fast(struct Jacobian *J, struct axe *axe, int N_axes,
			atom *atoms, int N_atoms, struct Reference Ref)
{
  // The computation assumes that the center of mass is at zero and
  // the global rotations and shift have been computed for all axes
  // in the routine Eckart_main

  int k, l, i, j, iref;

  // Compute inertia tensor and center of mass for all atoms
  double **corr_sum=Allocate_mat2_d(3, 3);
  double **inertia=Allocate_mat2_d(3, 3);
  double M_tot=0, mr_tot[3];
  for(i=0; i<3; i++)mr_tot[i]=0;
  for(iref=0; iref< Ref.N_ref; iref++){
    double *r=atoms[Ref.atom_num[iref]].r;
    double  m= Ref.mass_atom[iref];
    M_tot+=m;
    for(i=0; i<3; i++){
      float mri=m*r[i]; mr_tot[i]+=mri;
      for(j=i; j<3; j++)corr_sum[i][j]+=mri*r[j];
    }
  }
  Inertia_tensor(inertia, corr_sum, mr_tot, M_tot);

  if(0){
    printf("Total: M= %.0f MR= %.4f %.4f %.4f\n",
	   M_tot, mr_tot[0], mr_tot[1], mr_tot[2]);
    printf("Inertia tensor:\n");
    printf("%8.0f %8.0f %8.0f\n%8.0f %8.0f %8.0f\n%8.0f %8.0f %8.0f\n",
	   inertia[0][0], inertia[0][1], inertia[0][2],
	   inertia[1][0], inertia[1][1], inertia[1][2],
	   inertia[2][0], inertia[2][1], inertia[2][2]);
  }

  // Compute kinetic energy tensor (only internal dof)
  // Lower tridiagonal matrix T_sqrt[l][k] with k<=l

  double I_tot_rot_l[3], I_part_rot_l[3];
  double **inertia_part=Allocate_mat2_d(3,3);
  for(l=0; l<N_axes; l++){
    struct axe *axe_l=axe+l;

    // Inertia tensor must NOT be in center of mass frame
    double *mr_l=axe_l->MR;
    double r2=Scalar_product_3_d(mr_l, mr_l)/axe_l->mass;
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
        inertia_part[i][j]=axe_l->I[i][j]-mr_l[i]*mr_l[j]/axe_l->mass;
      }
      inertia_part[i][i]+=r2;
    }
    /*for(i=0; i<3; i++){
      for(j=0; j<3; j++)inertia_part[i][j]=axe_l->I[i][j];
      }*/
    Matrix_multiplication(I_part_rot_l, inertia_part, axe_l->rot, 3);
    Matrix_multiplication(I_tot_rot_l, inertia, axe_l->global_rot, 3);

    for(k=0; k<=l; k++){
      struct axe *axe_k=axe+k;
      double T=0, tv_kl[3], tv_lk[3];
 
      // Intersection of k and l
      // We use the fact that axes are nested: axe l>=k moves a
      // subset of atoms moved by axe k, unless side chain
      if(Dof_overlap(axe_k, axe_l)){
	T+=Scalar_product_3_d(I_part_rot_l, axe_k->rot);
	T+=axe_l->mass*Scalar_product_3_d(axe_l->shift, axe_k->shift);
	Vector_product_d(tv_kl, axe_k->shift, axe_l->rot);
	Vector_product_d(tv_lk, axe_l->shift, axe_k->rot);
	double tv[3]; for(j=0; j<3; j++)tv[j]=tv_kl[j]+tv_lk[j];
	T+=Scalar_product_3_d(tv, axe_l->MR);
      }

      // If Eckart conditions have been imposed, it is sufficient
      // to add the following terms:
      T-=M_tot*Scalar_product_3_d(axe_l->global_shift,axe_k->global_shift);
      T-=Scalar_product_3_d(I_tot_rot_l, axe_k->global_rot);
      J->T_sqrt[l][k]=T;
      if(l!=k)J->T_sqrt[k][l]=T;
    }
  }

  if(DEBUG){
    int ini=N_axes-20;
    printf("###  Kinetic energy matrix (calc):\n");
    for(k=ini; k<ini+20; k++){
      for(l=ini; l<=k; l++)printf(" %.3g", J->T_sqrt[k][l]);
      printf("\n");
    }
  }

  Empty_matrix_d(inertia_part, 3);
  Empty_matrix_d(corr_sum, 3);
  Empty_matrix_d(inertia, 3);

  return(0);
}

int Kinetic_energy_Jac(struct Jacobian J, struct Reference Ref, int N_axes)
{
  // Compute kinetic energy tensor T=T_sqrt T_sqrt^T (only internal dof)
  // Lower tridiagonal matrix T_sqrt[i][j] with j<=i
  int k, l, j, diff=0;
  for(l=0; l<N_axes; l++){
    float *Jl=J.Jacobian_ar[l];
    for(k=0; k<=l; k++){
      float *Jk=J.Jacobian_ar[k]; double sum=0;
      for(j=0; j<Ref.N_cart; j++)sum+=Jk[j]*Jl[j]*Ref.mass_coord[j];
      if(fabs(sum-J.T_sqrt[l][k])>0.0001*fabs(J.T_sqrt[l][l])){
	printf("%d %d %.6f %.6f\n", l, k, sum, J.T_sqrt[l][k]);
	diff++;
      }
      J.T_sqrt[l][k]=sum;
      if(l!=k)J.T_sqrt[k][l]=J.T_sqrt[l][k];
    }
  }
  if(diff){
    printf("ERROR, %d components of the kinetic energy are different\n",
	   diff); exit(8);
  }else{
    printf("All components of the kinetic energy are equal\n"); exit(8);
  }
  if(DEBUG){
    int ini=N_axes-20;
    printf("###  Kinetic energy matrix (calc):\n");
    for(k=ini; k<ini+20; k++){
      for(l=ini; l<=k; l++)printf(" %.4g", J.T_sqrt[k][l]);
      printf("\n");
    }
  }
  return(0);
}


int Kinetic_sqrt(struct Jacobian *J, int *N_kin,
		 int N, int CONTROL, float E_MIN)
{
  // Compute kinetic energy tensor T=T_sqrt T_sqrt^T (only internal dof)
  // Lower tridiagonal matrix T_sqrt[i][j] with j<=i
  int k, l;
  (*N_kin)=N;

  if(CONTROL==0){
    printf("Kinetic_sqrt with Cholevsky decomposition\n");
    if(choldc(J->T_sqrt, J->T_sqrt, N)==0){
      //J->T_sqrt_inv=Cholevsky_inversion_d2f(J->T_sqrt, N);
      for(k=0; k<N; k++){
	double *Lk=J->T_sqrt[k]+k+1;
	for(l=k+1; l<N; l++){*Lk=0; Lk++;}
      }
      return(0);
    }
  }

  // New decomposition diagonalizing kinetic energy
  printf("Square root of kinetic energy through diagonalization\n");
  {
    float **eigen_vector=Allocate_mat2_f(N, N);
    float eigen_value[N], t, t_inv;
    int num=0, ik;
    if(J->T_sqrt_inv)Empty_matrix_f(J->T_sqrt_inv, N);
    J->T_sqrt_inv=Allocate_mat2_f(N, N);
    d_Diagonalize(N, J->T_sqrt, eigen_value, eigen_vector, 1);
    for(k=0; k<N; k++){
      if(eigen_value[k]<E_MIN){t=0; num++; continue;}
      else{t=sqrt(eigen_value[k]); t_inv=1./t; ik=k-num;}
      for(l=0; l<N; l++){
	J->T_sqrt[l][ik]=t*eigen_vector[k][l];
	if(t)J->T_sqrt_inv[ik][l]=eigen_vector[k][l]*t_inv;
      }
    }
    if(num)printf("WARNING, kinetic energy has %d null modes\n", num);
    (*N_kin)-=num;
    Empty_matrix_f(eigen_vector, N);
  }

  return(0);
}

void  Set_rot_shift(struct axe *axes, int N_axes, atom *atoms, int N_atoms)
{
  int k, j;
  double v[3], offset[3], shift[3];
  // Axes vector and shift
  for(k=0; k<N_axes; k++){
    struct axe *axe=axes+k;
    double *r1=axe->axe->previous->atom->r;
    double *r2=axe->axe->atom->r;
    double v2=0;
    for(j=0; j<3; j++){
      axe->offset[j]=r1[j];               // Get offset as r1
      offset[j]=axe->offset[j];
      v[j]=r2[j]-r1[j]; v2+=v[j]*v[j];    // Get axe as (r1-r2)/|r1-r2|
    }
    v2=1./sqrt(v2);
    for(j=0; j<3; j++)axe->vers[j]=v[j]*v2;

    if(axe->type=='l'){ // bond length
      for(j=0; j<3; j++){
	// Bonds used as dof have a standard length LBOND=1
	axe->shift[j]=axe->vers[j];
	axe->rot[j]=0;
      }
    }else if(axe->type=='a'){ // bond angle
      Normal_vector(v, axe->bond, axe->vers);
      Vector_product_d(shift, v, offset);
      for(j=0; j<3; j++){
	axe->shift[j]=-shift[j];
	axe->rot[j]=v[j];
      }
    }else{   // torsion
      Vector_product_d(shift, axe->vers, offset);
      for(j=0; j<3; j++){
	axe->shift[j]=-shift[j]; axe->rot[j]=axe->vers[j];
      }
    }
  }
}

void Center_atoms(atom *atoms, int N_atoms, struct Reference Ref)
{
  double r_tot[3], M_tot=0; int iref, i, j;
  for(j=0; j<3; j++)r_tot[j]=0;
  for(iref=0; iref< Ref.N_ref; iref++){
    double *r=atoms[Ref.atom_num[iref]].r;
    double m = Ref.mass_atom[iref];
    M_tot+=m; for(j=0; j<3; j++)r_tot[j]+=m*r[j];
  }
  for(j=0; j<3; j++)r_tot[j]/=M_tot;
  for(i=0; i<N_atoms; i++){
    for(j=0; j<3; j++)atoms[i].r[j]-=r_tot[j];
  }
}

void Eckart_main(struct axe *axes, int N_axes, atom *atoms, int N_atoms,
		 struct Reference Ref)
{
  /******************************************************************/
  /*   Computing internal rotations and translations (Eckart frame)  */
  /******************************************************************/
  // Initialize
  int i, j, k;
  double **corr_sum=Allocate_mat2_d(3,3), mr_tot[3];
  double M_tot=Sum_inertia(mr_tot, corr_sum, atoms, Ref, 0, Ref.N_ref-1);
  double r_ave[3]; for(i=0; i<3; i++)r_ave[i]=mr_tot[i]/M_tot;
  double **inertia_tot=Allocate_mat2_d(3,3);
  Inertia_tensor(inertia_tot, corr_sum, mr_tot, M_tot);

  // Cholevsky decomposition of global inertia tensor
  double **L_mat=Allocate_mat2_d(3,3);
  choldc(L_mat, inertia_tot, 3);

  // Partial inertia tensor for all axes
  // Exploits the fact that axes are nested
  double mr_part[3], M_part=0, **inertia_part=Allocate_mat2_d(3,3);
  struct axe *axe=axes+N_axes-1;
  int last_kin=-1, i_last=0;
  for(int k=N_axes-1; k>=0; k--){
    // Compute partial center of mass and inertia tensor
    if(axe->last_kin!=last_kin){
      // New chain, initialize counters
      M_part=Empty_inertia(mr_part, corr_sum);
      last_kin=axe->last_kin;
      i_last=last_kin;
    }
    M_part+=Sum_inertia(mr_part, corr_sum, atoms, Ref, axe->first_kin, i_last);
    i_last=axe->first_kin-1;
    Inertia_tensor(inertia_part, corr_sum, mr_part, M_part);

    // Test that at least some atom is moved
    if(M_part==0){
      printf("WARNING, axis %d does not move any reference atom!\n", k);
      printf("Reference atoms: %d - %d\n", axe->first_kin, axe->last_kin);
      printf("All atoms: %d - %d\n", axe->first_atom, axe->last_atom);
      printf("dof: %d-%d (main) ", axe->first_main, axe->last_main);
      if(axe->first_side>=0)printf(" %d-%d (side)", axe->first_side, k);
      printf("\n");
      //exit(8);
    }
 
    /*
      ECKART ROTATION AND SHIFT OPERATORS:
      Vector global rotation about axe k:
      Iner global_rot[k]=-iner_part[k]v[k]+M_k(R-Rk)X(axe[k]X Rk+tau[k])}
      Vector global shift:
      global_shift[k]= - (axe[k]X Rk+tau[k])Mk/M - global_rot[k] X R
     */

    // Center of mass displacement
    double delta_cm[3], r_shift[3], rot_tot[3], X[3], Y[3];
    Vector_product_d(delta_cm, axe->rot, mr_part);
    for(i=0; i<3; i++)delta_cm[i]=delta_cm[i]/M_part+axe->shift[i];
    for(i=0; i<3; i++)r_shift[i]=M_part*r_ave[i]-mr_part[i];
    Vector_product_d(Y, r_shift, delta_cm);
    for(i=0; i<3; i++){
      for(j=0; j<3; j++)Y[i]-=inertia_part[i][j]*axe->rot[j];
    }

    // Downstream rotataion: Solve (L L^t)A_a = Y_a
    Forward_substitution(X, L_mat, Y, 3);
    Backward_substitution(axe->global_rot, L_mat, X, 3);

    // Downstream shift
    Vector_product_d(rot_tot, axe->global_rot, r_ave);
    for(i=0; i<3; i++){
      axe->global_shift[i]=-delta_cm[i]*M_part/M_tot-rot_tot[i]; //
    }

    // Downstream rotation and shift
    for(i=0; i<3; i++){
      axe->local_rot[i]  =axe->global_rot[i] + axe->rot[i];
      axe->local_shift[i]=axe->global_shift[i]+axe->shift[i];
    }

    // Store kinetic variables
    axe->mass=M_part;
    for(i=0; i<3; i++){
      axe->MR[i]=mr_part[i];
      for(j=0; j<3; j++)axe->I[i][j]=inertia_part[i][j];
    }

    // New axe
    axe--;
  }

  if(DEBUG){
    double *v, q;
    printf("Axis res local_rot local_shift global_rot global_shift");
    printf(" first_ref first_rotable_atom type\n");
    for(k=0; k<N_axes; k++){
      axe=axes+k;
      printf("%3d %3d ", k, axe->axe->atom->res);
      v=axe->local_rot;
      q=Scalar_product_3_d(v,v); printf(" %.3f ", sqrt(q));
      v=axe->local_shift;
      q=Scalar_product_3_d(v,v); printf("%4.1f ", sqrt(q));
      v=axe->global_rot;
      q=Scalar_product_3_d(v,v); printf(" %.3f ", sqrt(q));
      v=axe->global_shift;
      q=Scalar_product_3_d(v,v); printf("%4.1f ", sqrt(q));
      printf("%3d %3d %c\n",axe->first_kin,axe->first_atom,axe->type);
    }
  }

  // Cleaning
  Empty_matrix_d(corr_sum, 3);
  Empty_matrix_d(inertia_part, 3);
  Empty_matrix_d(inertia_tot, 3);
  Empty_matrix_d(L_mat, 3);
}

int Compute_Jtilde(struct Jacobian *J)
{
  // Jtilde = J L^t(-1)
  // Jtilde_ia = sum_b J_ib L^t(-1)_ba =  sum_b J_ib L(-1)_ab
  // Jtilde_ai = L(-1)_ab sum_b J_bi => LJtilde = J
  int a, i, b;

  if(J->Jtilde_ar==NULL)
    J->Jtilde_ar=Allocate_mat2_f(J->N_axes, J->N_cart);
  if(J->T_sqrt_inv!=NULL){
    // T_sqrt is not lower diagonal!
    for(a=0; a<J->N_kin; a++){
      for(i=0; i<J->N_cart; i++){
	double sum=0;
	for(b=0; b<J->N_axes; b++)
	  sum+=J->T_sqrt_inv[a][b]*J->Jacobian_ar[b][i];
	J->Jtilde_ar[a][i]=sum;
      }
    }
  }else{
    // Exploit that T_sqrt is lower diagonal
    // T_sqrt Jtilde = J
    double *X=malloc(J->N_axes*sizeof(double));
    double *Y=malloc(J->N_axes*sizeof(double));
    for(i=0; i<J->N_cart; i++){
      for(a=0; a<J->N_axes; a++)Y[a]=J->Jacobian_ar[a][i];
      Forward_substitution(X, J->T_sqrt, Y, J->N_axes);
      for(a=0; a<J->N_axes; a++)J->Jtilde_ar[a][i]=X[a];
    }
    free(X); free(Y);
  }
  return(0);
}

int Correlate_Jacobian(struct Jacobian J,
		       struct Reference Ref,
		       struct axe *axe, int N_axes)
{
  // Jacobian: derivatives with respect to torsions for all atoms
  // (only internal degrees of freedom)
  char nameout[100]="Jacobian_correlations.dat";
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);
  fprintf(file_out, "#r(J[k], index) loc_shift glob_shift loc_rot glob_rot\n");
  int N=Ref.N_ref, i, k;
  float D2[N], index[N];
  for(i=0; i<N; i++)index[i]=i;
  for(k=0; k<N_axes; k++){
    float *J_k=J.Jacobian_ar[k];
    for(i=0; i<N; i++){
      double D=0;
      D+=(*J_k)*(*J_k); J_k++;
      D+=(*J_k)*(*J_k); J_k++; 
      D+=(*J_k)*(*J_k); J_k++;
      D2[i]=D;
    }
    float slope, offset, r=Corr_coeff(index, D2, N, &slope, &offset);
    struct axe *axe_k=axe+k;
    fprintf(file_out, "%.3f %c", r, axe_k->type);
    float v2=Module(axe_k->local_shift); fprintf(file_out, " %.3f", v2);
    v2=Module(axe_k->global_shift);fprintf(file_out, " %.3f", v2);
    v2=Module(axe_k->local_rot);   fprintf(file_out, " %.3f", v2);
    v2=Module(axe_k->global_rot);  fprintf(file_out, " %.3f\n", v2);
  }
  fclose(file_out); //exit(8);
  return(0);
}

float Module(double *v){
  return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

int Internal_Jacobian(// Output:
		      struct Jacobian J,
		      // Input:
		      struct Reference Ref,
		      atom *atoms, int N_atoms,
		      struct axe *axe, int N_axes)
{
  // Jacobian: derivatives with respect to torsions for all atoms
  // (only internal degrees of freedom)
  int i, j, k;
  for(k=0; k<N_axes; k++){
    float *Jacobian_a=J.Jacobian_ar[k]; int jj=0;
    struct axe *axe_k=axe+k;
    for(i=0; i<Ref.N_ref; i++){
      double *rot,  *shift, dk[3];
      if((i>=axe_k->first_kin)&&(i<=axe_k->last_kin)){
	rot=axe[k].local_rot;  shift=axe[k].local_shift;
      }else{
	rot=axe[k].global_rot; shift=axe[k].global_shift;
      }
      double *ri=atoms[Ref.atom_num[i]].r;
      double r[3]; for(j=0; j<3; j++)r[j]=ri[j];
      Vector_product_d(dk, rot, r);
      for(j=0; j<3; j++){
	Jacobian_a[jj]=(dk[j]+shift[j]); jj++;
      }
    }
  }

  return(0);
}

void Eckart_test(float **Jacobian_ar, int N_axes,
		 atom *atoms, int N_atoms, struct Reference Ref)
{
  printf("Testing Eckart conditions for each axis (%d)\n", N_axes);
  int EV=0;
  for(int k=0; k<N_axes; k++){
    EV+=Test_Eckart(Jacobian_ar[k], atoms, N_atoms, Ref, k);
  }
  if(EV){
    printf("ERROR, %d violations of the Eckart conditions\n", EV);
    exit(8);
  }
}

int Test_Eckart(float *Cart_dev, atom *atoms, int N_atoms,
		 struct Reference Ref, int mode)
{
  // r_sum = sum_i m_i Delta_r_i
  // r_rot = sum_i r_i X Delta_r_i
  double r_sum[3], r_rot_sum[3], norm=0;
  double r_rot[3], r1=0, r2=0, x[3];
  int iref, j, jj=0;
  for(j=0; j<3; j++){r_sum[j]=0; r_rot_sum[j]=0;}
  for(iref=0; iref< Ref.N_ref; iref++){
    atom *atom1=atoms+Ref.atom_num[iref];
    double m = Ref.mass_atom[iref];
    for(j=0; j<3; j++){x[j]=Cart_dev[jj]; r_sum[j]+=m*x[j]; jj++;}
    Vector_product_d(r_rot, atom1->r, x);
    for(j=0; j<3; j++)r_rot_sum[j]+=m*r_rot[j];
    norm+=m;
  }
  for(j=0; j<3; j++){
    r1+= r_sum[j]*r_sum[j];
    r2+= r_rot_sum[j]*r_rot_sum[j];
  }
  r1=sqrt(r1)/norm; r2=sqrt(r2)/norm;
  if((r1>0.001)||(r2>0.001)){
    printf("%3d <d>= %.8f <r_X_d>= %.8f\n", mode, r1, r2);
    return(1);
  }
  return(0);
}


void Inertia_tensor(double **inertia, double **corr_sum,
	            double *r_sum, double M_sum)
{
  // Only upper diagonal of corr_sum is needed
  int i, j; double r2;
  for(i=0; i<3; i++){
    for(j=i; j<3; j++)
      inertia[i][j]=corr_sum[i][j]-r_sum[i]*r_sum[j]/M_sum;
  }
  r2=inertia[0][0]+inertia[1][1]+inertia[2][2];
  for(i=0; i<3; i++){
    inertia[i][i]=r2-inertia[i][i];
    for(j=i+1; j<3; j++){
      inertia[i][j]=-inertia[i][j];
      inertia[j][i]= inertia[i][j];
    }
  }
}

float **Transpose_matrix(float **M, int n1, int n2){
  float **M_tr=Allocate_mat2_f(n2, n1); int i, j;
  for(i=0; i<n2; i++){
    for(j=0; j<n1; j++)M_tr[i][j]=M[j][i];
  }
  return(M_tr);
}

float **Transpose_matrix_d2f(double **M, int n1, int n2){
  float **M_tr=Allocate_mat2_f(n2, n1); int i, j;
  for(i=0; i<n2; i++){
    for(j=0; j<n1; j++)M_tr[i][j]=M[j][i];
  }
  return(M_tr);
}

void Normal_vector(double *v, struct bond *bond, double *vers)
{
  double *r0=NULL; int j;
  if((bond)&&(bond->previous)&&(bond->previous->previous)&&
     (bond->previous->previous->atom)){
    r0=bond->previous->previous->atom->r;
  }
  if(r0==NULL){
    printf("ERROR, no previous atom found for bond angle\n");
    goto error_normal;
  }
  double *r1=bond->previous->atom->r, dr[3];
  for(j=0; j<3; j++)dr[j]=r1[j]-r0[j];
  Vector_product_d(v, dr, vers);
  // Normalize:
  double v2=0; for(j=0; j<3; j++)v2+=v[j]*v[j];
  if(v2<=0){
    printf("ERROR for bond angle, zero norm (%.6f) of rotation axis\n",v2);
    goto error_normal;
  }
  v2=1./sqrt(v2); for(j=0; j<3; j++)v[j]*=v2;
  return;
 error_normal:
  printf("atoms: %s %d %s %d prev. bond: %d prev. prev.: %d\n",
	 bond->previous->atom->name, bond->previous->atom->res,
	 bond->atom->name, bond->atom->res,
	 bond->previous->i_atom, bond->previous->previous->i_atom);
  exit(8);
}

int Dof_overlap(struct axe *axe_k, struct axe *axe_l)
{
  // dof k<=l. Returns 1 if they move a common subset of atoms
  if((axe_l->first_kin>=axe_k->first_kin)&&
     (axe_l->last_kin <=axe_k->last_kin)){
    return(1);
  }else if((axe_k->last_kin < axe_l->first_kin)||
	   (axe_l->last_kin < axe_k->first_kin)){
    return(0);
  }
  printf("ERROR, not nested axes\n");
  printf("k: type %c first_kin: %d last_kin: %d\n",
	 axe_k->type, axe_k->first_kin, axe_k->last_kin);
  printf("l: type %c first_kin: %d last_kin: %d\n",
	 axe_l->type, axe_l->first_kin, axe_l->last_kin);
  exit(8);
}

void Allocate_Jacobian(struct Jacobian *J, int N_axes, int N_cart)
{
  J->N_axes=N_axes; J->N_cart=N_cart;
  J->T_sqrt=Allocate_mat2_d(N_axes,N_axes);
  J->Jacobian_ar=Allocate_mat2_f(N_axes,N_cart);
  J->Jtilde_ar=NULL;  // Allocate_mat2_f(N_axes,N_cart);
  J->T_sqrt_tr=NULL;
  J->T_sqrt_inv=NULL;
  J->T_sqrt_inv_tr=NULL;
}

void Empty_Jacobian(struct Jacobian J)
{
  if(J.Jacobian_ar!=NULL)Empty_matrix_f(J.Jacobian_ar, J.N_axes);
  if(J.Jtilde_ar!=NULL)Empty_matrix_f(J.Jtilde_ar, J.N_axes);
  if(J.T_sqrt!=NULL)Empty_matrix_d(J.T_sqrt, J.N_axes);
  if(J.T_sqrt_inv!=NULL)Empty_matrix_f(J.T_sqrt_inv, J.N_axes);
}

float Empty_inertia(double *mr_tot, double **corr_sum)
{
  for(int i=0; i<3; i++){
    mr_tot[i]=0; for(int j=0; j<3; j++)corr_sum[i][j]=0;
  }
  return(0.);
}

float Sum_inertia(double *mr_tot, double **corr_sum,
		  atom *atoms, struct Reference Ref,
		  int ini_ref, int end_ref)
{
  double M_tot=0; int i, j;
  for(i=0; i<3; i++)for(j=i; j<3; j++)corr_sum[i][j]=0;
  for(int iref=ini_ref; iref<=end_ref; iref++){
    double *r=atoms[Ref.atom_num[iref]].r;
    double m = Ref.mass_atom[iref];
    M_tot+=m;
    for(i=0; i<3; i++){
      float mri=m*r[i]; mr_tot[i]+=mri;
      for(j=i; j<3; j++)corr_sum[i][j]+=mri*r[j];
    }
  }
  return(M_tot);
}

void Reorient_axes(atom *atoms, int N_atoms, struct Reference Ref)
{
  double **inertia_tot=Allocate_mat2_d(3,3);
  double **corr_sum=Allocate_mat2_d(3,3);
  double r_ave[3], mr_tot[3]; int i, k, l;

  // Compute inertia tensor and center of mass for all atoms
  double M_tot=Sum_inertia(mr_tot, corr_sum, atoms, Ref, 0, Ref.N_ref-1);
  Inertia_tensor(inertia_tot, corr_sum, mr_tot, M_tot);

  for(i=0; i<3; i++)r_ave[i]=mr_tot[i]/M_tot;
  printf("Center of mass: ");
  for(i=0; i<3; i++)printf("  %.7f", r_ave[i]);
  printf("  M=%.1f\n", M_tot);

  float eigen_value[3], **rot=Allocate_mat2_f(3, 3);
  d_Diagonalize(3, inertia_tot, eigen_value, rot, 1);
  printf("Inertia_values: ");
  for(k=0; k<3; k++)printf(" %.0f", eigen_value[k]); printf("\n");
  printf("Principal axes:\n");
  for(k=0; k<3; k++){
    float *v=rot[k];
    for(l=0; l<3; l++)printf(" %.3f", v[l]); printf("\n");
  }

  // Rotate atoms coordinates
  for(i=0; i<N_atoms; i++){
    double rr[3], *rk=rr, *r=atoms[i].r;
    for(k=0; k<3; k++){
      *rk=0; for(l=0; l<3; l++)*rk+=rot[k][l]*r[l]; rk++;
    }
    for(k=0; k<3; k++)r[k]=rr[k];
  }

  if(0){
    M_tot=Sum_inertia(mr_tot, corr_sum, atoms, Ref, 0, Ref.N_ref-1);
    Inertia_tensor(inertia_tot, corr_sum, mr_tot, M_tot);
    d_Diagonalize(3, inertia_tot, eigen_value, rot, 1);
    printf("Principal axes:\n");
    for(k=0; k<3; k++){
      float *v=rot[k];
      for(l=0; l<3; l++)printf(" %.3f", v[l]); printf("\n");
    }
    exit(8);
  }

  Empty_matrix_d(corr_sum, 3);
  Empty_matrix_d(inertia_tot, 3);
  Empty_matrix_f(rot, 3);
}

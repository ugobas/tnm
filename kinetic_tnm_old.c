int DBG=0;
#define T_COMPUTE 1 // 0=Cholevsky 1=Diagonalization
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
void Center_atoms(atom *atoms, int N_atoms, struct axe *axes, int N_axes,
		  struct Reference Ref);
void Eckart_main(struct axe *axes, int N_axes, atom *atoms, int N_atoms,
		 struct Reference Ref);
void Inertia_tensor(double **inertia, double **corr_sum, double *r_sum,
		    double M_sum);
float Sum_inertia(double *mr_tot, double **corr_sum,
		  atom *atoms, struct Reference Ref,
		  int ini_ref, int end_ref);
float Empty_inertia(double *mr_tot, double **corr_sum);

int Kinetic_energy(double ***T_mat, float *mass_coord,
		   struct axe *axe, int ini_axe, int N_axes,
		   atom *atoms, int *atom_num, int N_ref);
int Kinetic_energy_all(struct Jacobian J, struct Reference Ref, int N_axes);
int Kinetic_energy_fast(struct Jacobian *J, struct axe *axe,
			int N_axes, int mainaxes, int rigidaxes, 
			atom *atoms, int N_atoms,
			struct Reference Ref);
int Kinetic_sqrt(struct Jacobian *J, int *N_kin,
		 int N, int CONTROL, float E_MIN);
int Internal_Jacobian(// Output:
		      struct Jacobian J,
		      // Input:
		      struct Reference Ref,
		      atom *atoms, int N_atoms,
		      struct axe *axe, int N_axes);
void Test_Eckart(float *Cart_dev, atom *atoms, int N_atoms,
		 struct Reference Ref, int mode);
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
void Normal_vector(double *v, struct axe *axe, atom *atoms, int k);
int Dof_overlap(int k, int l, struct axe *axe_k, struct axe *axe_l,
		int rigidaxes, int mainaxes, struct axe *axe);

int Compute_kinetic(struct Jacobian *J,
		    struct axe *axe1,
		    int naxe1, int mainaxes, int rigidaxes,
		    atom *atoms1, int natoms1,
		    struct Reference Ref)
{
  int N_kin=0;
  Center_atoms(atoms1, natoms1, axe1, naxe1, Ref);
  Set_rot_shift(axe1, naxe1, atoms1, natoms1);
  Eckart_main(axe1, naxe1, atoms1, natoms1, Ref);
  Internal_Jacobian(*J, Ref, atoms1, natoms1, axe1, naxe1);
  if(J->T_sqrt_inv){
    Empty_matrix_f(J->T_sqrt_inv, naxe1); J->T_sqrt_inv=NULL;
  }
  //if(rigidaxes==0){
    Kinetic_energy_fast(J,axe1,naxe1,mainaxes,rigidaxes,atoms1,natoms1,Ref);
  //}else{
  //  Kinetic_energy_all(*J, Ref, naxe1);
  //}
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
  if(J->Rot_ra!=NULL){
    Empty_matrix_f(J->Rot_ra, 3*Ref.N_ref); J->Rot_ra=NULL;
  }
  if(J->Shift_ra!=NULL){
    Empty_matrix_f(J->Shift_ra, 3*Ref.N_ref); J->Shift_ra=NULL;
  }

  // Show results
  if(DBG){
    int i, j, m=10;
    printf("Masses: ");
    for(i=0; i<3*m; i+=3)printf("%.0f ", Ref.mass_coord[i]);
    printf("\nJacobian: ");
    for(i=0; i<m; i++)printf("%.1f ", J->Jacobian_ar[i][0]);
    printf("\nKinetic energy matrix. ");
    printf("%d dofs rigid= %d main= %d atoms= %d cart= %d Nkin= %d\n",
	   naxe1, rigidaxes, mainaxes, natoms1, Ref.N_cart, N_kin);
    for(i=0; i<m; i++){
      for(j=0; j<m; j++)printf("%.2g\t", J->T_sqrt[i][j]);
      printf("\n");
    }
    printf("Axis properties\n");
    struct axe *axe=axe1;
    for(i=0; i<30; i++){
      printf("%c%d\t%d-%d\t%d-%d\t%d-%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
	     axe->type, i, axe->atom1, axe->atom2,
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
  float **J_ar=Allocate_mat2_f(N_axes, N_cart), *J, r[3];
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
      Subtract_vector_3(r, atoms[atom_num[i]].r, offset);
      Vector_product(J, rot, r); J+=3;
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
			int mainaxes, int rigidaxes, atom *atoms, int N_atoms,
			struct Reference Ref)
{
  // The computation assumes that the center of mass is at zero and
  // the global rotations and shift have been computed for all axes
  // in the routine Eckart_main

  int k, l, i, j, iref;

  // Compute inertia tensor and center of mass for all atoms
  double **corr_sum=Allocate_mat2_d(3, 3);
  double **inertia=Allocate_mat2_d(3, 3);
  double M_tot=0, r_tot[3], mri;
  for(i=0; i<3; i++)r_tot[i]=0;
  for(iref=0; iref< Ref.N_ref; iref++){
    float  *r=atoms[Ref.atom_num[iref]].r;
    double  m= Ref.mass_atom[iref];
    M_tot+=m;
    for(i=0; i<3; i++){
      mri=m*r[i]; r_tot[i]+=mri;
      for(j=i; j<3; j++)corr_sum[i][j]+=mri*r[j];
    }
  }
  Inertia_tensor(inertia, corr_sum, r_tot, M_tot);

  if(0){
    printf("Total: M= %.0f MR= %.4f %.4f %.4f\n",
	   M_tot, r_tot[0], r_tot[1], r_tot[2]);
    printf("Inertia tensor:\n");
    printf("%8.0f %8.0f %8.0f\n%8.0f %8.0f %8.0f\n%8.0f %8.0f %8.0f\n",
	   inertia[0][0], inertia[0][1], inertia[0][2],
	   inertia[1][0], inertia[1][1], inertia[1][2],
	   inertia[2][0], inertia[2][1], inertia[2][2]);
  }

  // Compute kinetic energy tensor (only internal dof)
  // Lower tridiagonal matrix T_sqrt[l][k] with k<=l

  double I_rot_l[3], **I_part_rot=Allocate_mat2_d(N_axes, 3);
  double **inertia_part=Allocate_mat2_d(3,3), r2, *mr_l;
  for(l=0; l<N_axes; l++){
    struct axe *axe_l=axe+l;

    // Inertia tensor must NOT be in center of mass frame
    mr_l=axe_l->MR;
    r2=Scalar_product_3_d(mr_l, mr_l)/axe_l->mass;
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
        inertia_part[i][j]=axe_l->I[i][j]-mr_l[i]*mr_l[j]/axe_l->mass;
      }
      inertia_part[i][i]+=r2;
    }
    Matrix_multiplication(I_rot_l, inertia, axe_l->global_rot, 3);
    Matrix_multiplication(I_part_rot[l], inertia_part, axe_l->rot, 3);

    if(0){
      printf("axe %d: M= %.0f MR= %.0f %.0f %.0f\n",
	     l, axe_l->mass, axe_l->MR[0], axe_l->MR[1], axe_l->MR[2]);
      printf("Inertia tensor:\n");
      printf("%8.0f %8.0f %8.0f\n%8.0f %8.0f %8.0f\n%8.0f %8.0f %8.0f\n",
	     inertia_part[0][0], inertia_part[0][1], inertia_part[0][2],
	     inertia_part[1][0], inertia_part[1][1], inertia_part[1][2],
	     inertia_part[2][0], inertia_part[2][1], inertia_part[2][2]);
    }


    for(k=0; k<=l; k++){
      struct axe *axe_k=axe+k;
      double tv1[3], tv2[3], T=0;

      // Intersection of k and l
      // We use the fact that axes are nested: axe l>=k moves a
      // subset of atoms moved by axe k, unless side chain

      if(Dof_overlap(k, l, axe_k, axe_l, rigidaxes, mainaxes, axe)){
	T+=Scalar_product_3_d(I_part_rot[l], axe_k->rot);
	T+=axe_l->mass*Scalar_product_3_d(axe_l->shift, axe_k->shift);
	Vector_product_d(tv1, axe_k->shift, axe_l->rot);
	Vector_product_d(tv2, axe_l->shift, axe_k->rot);
	for(j=0; j<3; j++)tv1[j]+=tv2[j];
	T+=Scalar_product_3_d(tv1, axe_l->MR);
      }

      // If Eckart conditions have been imposed, it is sufficient
      // to add the following two terms
      T-=M_tot*Scalar_product_3_d(axe_l->global_shift,axe_k->global_shift);
      T-=Scalar_product_3_d(I_rot_l, axe_k->global_rot);


      /*// All atoms
      T+=Scalar_product_3_d(I_rot_l, axe_k->global_rot);
      T+=M_tot*Scalar_product_3_d(axe_l->global_shift,axe_k->global_shift);
      Vector_product_d(tv1, axe_k->global_shift, axe_l->global_rot);
      Vector_product_d(tv2, axe_l->global_shift, axe_k->global_rot);
      for(j=0; j<3; j++)tv1[j]+=tv2[j];
      T+=Scalar_product_3_d(tv1, r_tot);

      // Only l
      T+=Scalar_product_3_d(I_part_rot[l], axe_k->global_rot);
      T+=axe_l->mass*Scalar_product_3_d(axe_l->shift, axe_k->global_shift);
      Vector_product_d(tv1, axe_k->global_shift, axe_l->rot);
      Vector_product_d(tv2, axe_l->shift, axe_k->global_rot);
      for(j=0; j<3; j++)tv1[j]+=tv2[j];
      T+=Scalar_product_3_d(tv1, axe_l->MR);

      // Only k
      T+=Scalar_product_3_d(I_part_rot[k], axe_l->global_rot);
      T+=axe_k->mass*Scalar_product_3_d(axe_l->global_shift, axe_k->shift);
      Vector_product_d(tv1, axe_l->global_shift, axe_k->rot);
      Vector_product_d(tv2, axe_k->shift, axe_l->global_rot);
      for(j=0; j<3; j++)tv1[j]+=tv2[j];
      T+=Scalar_product_3_d(tv1, axe_k->MR);*/

      /*int L=N_axes-6; T=T+T0;
      if(((l<5)&&(k<5))||((l>L)&&(k>L))){
	printf("%3d %3d %.0f %.0f %.0f\n",
	l, k, J->T_sqrt[l][k], T, T0-T1);
      }*/

      J->T_sqrt[l][k]=T;
    }
  }
  for(l=0; l<N_axes; l++){
    for(k=l+1; k<N_axes; k++)J->T_sqrt[l][k]=J->T_sqrt[k][l];
  }
  if(DEBUG){
    int ini=N_axes-20;
    printf("###  Kinetic energy matrix (calc):\n");
    for(k=ini; k<ini+20; k++){
      for(l=ini; l<=k; l++)printf(" %.3g", J->T_sqrt[k][l]);
      printf("\n");
    }
  }

  Empty_matrix_d(I_part_rot, N_axes);
  Empty_matrix_d(inertia_part, 3);
  Empty_matrix_d(corr_sum, 3);
  Empty_matrix_d(inertia, 3);

  return(0);
}


int Kinetic_energy_all(struct Jacobian J, struct Reference Ref, int N_axes)
{
  // Compute kinetic energy tensor T=T_sqrt T_sqrt^T (only internal dof)
  // Lower tridiagonal matrix T_sqrt[i][j] with j<=i
  int k, l, j;
  for(l=0; l<N_axes; l++){
    float *Jl=J.Jacobian_ar[l];
    for(k=0; k<=l; k++){
      float *Jk=J.Jacobian_ar[k]; double sum=0;
      for(j=0; j<Ref.N_cart; j++)sum+=Jk[j]*Jl[j]*Ref.mass_coord[j];
      J.T_sqrt[l][k]=sum;
    }
  }
  for(l=0; l<N_axes; l++){
    for(k=l+1; k<N_axes; k++)J.T_sqrt[l][k]=J.T_sqrt[k][l];
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
    float *r1=atoms[axe->atom1].r;
    float *r2=atoms[axe->atom2].r;
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
      Normal_vector(v, axe, atoms, k);
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

void Center_atoms(atom *atoms, int N_atoms, struct axe *axes, int N_axes,
		  struct Reference Ref)
{
  double r_tot[3], M_tot=0; int iref, i, j;
  for(j=0; j<3; j++)r_tot[j]=0;
  for(iref=0; iref< Ref.N_ref; iref++){
    float *r=atoms[Ref.atom_num[iref]].r;
    double m = Ref.mass_atom[iref];
    M_tot+=m; for(j=0; j<3; j++)r_tot[j]+=m*r[j];
  }
  for(j=0; j<3; j++)r_tot[j]/=M_tot;
  for(i=0; i<N_atoms; i++){
    for(j=0; j<3; j++)atoms[i].r[j]-=r_tot[j];
  }
  for(i=0; i<N_axes; i++){
    for(j=0; j<3; j++)axes[i].offset[j]-=r_tot[j];
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
  double r_tot[3], mr_tot[3], mr_part[3];
  double M_tot=0, M_part=0;
  double **corr_sum=Allocate_mat2_d(3,3);
  double **inertia_tot=Allocate_mat2_d(3,3);
  double **inertia_part=Allocate_mat2_d(3,3);
  double **L_mat=Allocate_mat2_d(3,3);

  // Compute inertia tensor and center of mass for all atoms
  M_tot+=Sum_inertia(mr_tot, corr_sum, atoms, Ref, 0, Ref.N_ref-1);
  for(i=0; i<3; i++)r_tot[i]=mr_tot[i]/M_tot;
  Inertia_tensor(inertia_tot, corr_sum, mr_tot, M_tot);
  // Cholevsky decomposition of global inertia tensor
  choldc(L_mat, inertia_tot, 3);

  // Partial inertia tensor for all axes
  // Exploits the fact that axes are nested
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
      printf("Atom last dof: %d (rigid) %d (main), %d (side)\n",
	     (atoms+Ref.atom_num[axe->first_kin])->last_dof_rigid,
	     (atoms+Ref.atom_num[axe->first_kin])->last_dof_main,
	     (atoms+Ref.atom_num[axe->first_kin])->last_dof_side);
      //exit(8);
    }
 
    /*
      ECKART ROTATION AND SHIFT OPERATORS:
      Vector global rotation about axe k:
      Iner global_rot[k]=-iner_part[k]axe[k]+M_k(R-Rk)X(axe[k]X Rk+tau[k])}
      Vector global shift:
      global_shift[k]= - (axe[k]X Rk+tau[k])Mk/M - global_rot[k] X R
     */

    // Center of mass displacement
    double delta_cm[3], rot_tot[3];
    double X[3], Y[3], r_shift[3];
    for(i=0; i<3; i++)r_shift[i]=M_part*r_tot[i]-mr_part[i];
    Vector_product_d(delta_cm, axe->rot, mr_part);
    for(i=0; i<3; i++)delta_cm[i]=delta_cm[i]/M_part+axe->shift[i];
    Vector_product_d(Y, r_shift, delta_cm);
    for(i=0; i<3; i++){
      for(j=0; j<3; j++)Y[i]-=inertia_part[i][j]*axe->rot[j];
    }
    // Downstream rotataion: Solve (L L^t)A_a = Y_a
    Forward_substitution(X, L_mat, Y, 3);
    Backward_substitution(axe->global_rot, L_mat, X, 3);

    // Downstream shift
    Vector_product_d(rot_tot, axe->global_rot, r_tot);
    for(i=0; i<3; i++){
      axe->global_shift[i]=-rot_tot[i]-delta_cm[i]*(M_part/M_tot);
    }

    // Upstream rotation and shift
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
      printf("%3d %3d ", k, atoms[axe->atom1].res);
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

int Internal_Jacobian(// Output:
		      struct Jacobian J,
		      // Input:
		      struct Reference Ref,
		      atom *atoms, int N_atoms,
		      struct axe *axe, int N_axes)
{
  int i, j, k;

  // Jacobian: derivatives with respect to torsions for all atoms
  // (only internal degrees of freedom)
  for(k=0; k<N_axes; k++){
    int jj=0;
    float *Jacobian_a=J.Jacobian_ar[k];
    double *rot,  *shift, r[3], dk[3];
    int ini_ref=axe[k].first_kin;
    int last_ref=axe[k].last_kin;
    for(i=0; i<Ref.N_ref; i++){
      atom *atom1=atoms+Ref.atom_num[i];
      for(j=0; j<3; j++)r[j]=atom1->r[j];

      if((i>=ini_ref)&&(i<=last_ref)){
	rot=axe[k].local_rot;  shift=axe[k].local_shift;
      }else{
	rot=axe[k].global_rot; shift=axe[k].global_shift;
      }
      Vector_product_d(dk, rot, r);
      for(j=0; j<3; j++){
	Jacobian_a[jj]=(dk[j]+shift[j]); jj++;
      }
    }
  }

  if(DEBUG){
    float *Cart_dev=malloc(Ref.N_cart*sizeof(float));
    printf("Testing Eckart condition for each axis (%d)\n", N_axes);
    for(k=0; k<N_axes; k++){
      float *Jacobian_a=J.Jacobian_ar[k];
      for(j=0; j<J.N_cart; j++)Cart_dev[j]=Jacobian_a[j];
      Test_Eckart(Cart_dev,atoms,N_atoms,Ref, k);
    }
    free(Cart_dev);
  }

  return(0);
}

int Compute_Rot_Shift(// Output:
		      struct Jacobian *J,
		      // Input:
		      struct axe *axe)
{
  if(J->Rot_ra==NULL){
    //printf("Allocating Rot Shift\n");
    J->Rot_ra=Allocate_mat2_f(J->N_cart, J->N_axes);
    J->Shift_ra=Allocate_mat2_f(J->N_cart, J->N_axes);
  }

  // Jacobian: derivatives with respect to torsions for all atoms
  // (only internal degrees of freedom)
  int i, j, k;
  for(k=0; k<J->N_axes; k++){
    int jj=0;
    double *rot,  *shift;
    int ini_ref=axe[k].first_kin;
    int last_ref=axe[k].last_kin;
    for(i=0; i<(J->N_cart/3); i++){
      if((i>=ini_ref)&&(i<=last_ref)){
	rot=axe[k].local_rot;   shift=axe[k].local_shift;
      }else{
	rot=axe[k].global_rot; shift=axe[k].global_shift;
      }
      for(j=0; j<3; j++){
	J->Rot_ra[jj][k]=rot[j];
	J->Shift_ra[jj][k]=shift[j];
	jj++;
      }
    }
  }
  return(0);
}

void Test_Eckart(float *Cart_dev, atom *atoms, int N_atoms,
		 struct Reference Ref, int mode)
{
  // r_sum = sum_i m_i Delta_r_i
  // r_rot = sum_i r_i X Delta_r_i
  double r_sum[3], r_rot_sum[3], norm=0;
  float r_rot[3], r1=0, r2=0;
  int iref, j, jj=0;
  for(j=0; j<3; j++){r_sum[j]=0; r_rot_sum[j]=0;}
  for(iref=0; iref< Ref.N_ref; iref++){
    atom *atom1=atoms+Ref.atom_num[iref];
    double m = Ref.mass_atom[iref];
    for(j=0; j<3; j++)r_sum[j]+=m*Cart_dev[jj+j];
    Vector_product(r_rot, atom1->r, Cart_dev+jj);
    for(j=0; j<3; j++)r_rot_sum[j]+=m*r_rot[j];
    norm+=m;
    jj+=3;
  }
  for(j=0; j<3; j++){
    r1+= r_sum[j]*r_sum[j];
    r2+= r_rot_sum[j]*r_rot_sum[j];
  }
  r1=sqrt(r1)/norm; r2=sqrt(r2)/norm;
  printf("%3d <d>= %.3f <r_X_d>= %.3f\n", mode, r1, r2);
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

void Normal_vector(double *v, struct axe *axe, atom *atoms, int k)
{
  float *r0; int i0=-1, j;
  if((axe->bond)&&(axe->bond->previous)&&
     (axe->bond->previous->previous)){
    i0=axe->bond->previous->previous->i_atom;
    if(i0>=0)r0=atoms[i0].r;
  }
  if(i0<0){
    printf("ERROR, no previous atom of axe %d (a) found\n", k);
    goto error_normal;
  }
  float *r1=atoms[axe->atom1].r;
  double dr[3]; for(j=0; j<3; j++)dr[j]=r1[j]-r0[j];
  Vector_product_d(v, dr, axe->vers);
  // Normalize:
  double v2=0; for(j=0; j<3; j++)v2+=v[j]*v[j];
  if(v2<=0){
    printf("ERROR for dof %d, zero norm (%.6f) of rotation axis\n",k,v2);
    goto error_normal;
  }
  v2=1./sqrt(v2); for(j=0; j<3; j++)v[j]*=v2;
  return;
 error_normal:
  printf("atoms: %d %s %d %d %s %d previous bond: %d previous previous: %d\n",
	 axe->atom1, atoms[axe->atom1].name, atoms[axe->atom1].res,
	 axe->atom2, atoms[axe->atom2].name, atoms[axe->atom2].res,
	 axe->bond->previous->i_atom,
	 axe->bond->previous->previous->i_atom);
  exit(8);
}

int Dof_overlap(int k, int l, struct axe *axe_k, struct axe *axe_l,
		int rigidaxes, int mainaxes, struct axe *axe)
{
  // dof k<=l. Returns 1 if they move a common subset of atoms
  if(l<rigidaxes){
    return(1);
  }else if(l<mainaxes){
    if((axe_l->atom->chain==axe_k->atom->chain)||
       (k<=axe_l->last_rigid))return(1);
  }else{  //l>mainaxes = side chain
    if(k<=axe[axe_l->last_main].last_rigid){
      return(1);
    }else if((axe_k->atom->chain==axe_l->atom->chain)&&(k<=axe_l->last_main)){
      return(1);
    }else if(axe_k->atom->res==axe_l->atom->res){
      return(1);
    }
  }
  return(0);
}

void Allocate_Jacobian(struct Jacobian *J, int N_axes, int N_cart)
{
  J->N_axes=N_axes; J->N_cart=N_cart;
  J->T_sqrt=Allocate_mat2_d(N_axes,N_axes);
  J->Jacobian_ar=Allocate_mat2_f(N_axes,N_cart);
  J->Jtilde_ar=NULL;  // Allocate_mat2_f(N_axes,N_cart);
  J->Rot_ra=NULL;     // Allocate_mat2_f(N_cart, N_axes);
  J->Shift_ra=NULL;   // Allocate_mat2_f(N_cart, N_axes);
  J->T_sqrt_tr=NULL;
  J->T_sqrt_inv=NULL;
  J->T_sqrt_inv_tr=NULL;
}

void Empty_Jacobian(struct Jacobian J)
{
  if(J.Jacobian_ar!=NULL)Empty_matrix_f(J.Jacobian_ar, J.N_axes);
  if(J.Jtilde_ar!=NULL)Empty_matrix_f(J.Jtilde_ar, J.N_axes);
  if(J.Shift_ra!=NULL)Empty_matrix_f(J.Shift_ra, J.N_cart);
  if(J.Rot_ra!=NULL)Empty_matrix_f(J.Rot_ra, J.N_cart);
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
  for(int iref=ini_ref; iref<=end_ref; iref++){
    float *r=atoms[Ref.atom_num[iref]].r;
    double m = Ref.mass_atom[iref];
    M_tot+=m;
    for(i=0; i<3; i++){
      float mri=m*r[i]; mr_tot[i]+=mri;
      for(j=i; j<3; j++)corr_sum[i][j]+=mri*r[j];
    }
  }
  return(M_tot);
}


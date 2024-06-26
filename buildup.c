#include "coord.h"
#include "buildup.h"
#include "vector.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define DEBUG 0
float BONDTHR1=2.2; //Threshold for bond
float BONDTHR1_2;
int Ini_buildup=0;
float TWOPI=6.2831853;
float  PI=3.1415926;

struct bond *Find_previous(struct bond *bond, int n, int N_atoms,
			   struct residue *seq);
int Make_frame(struct bond *bond, struct bond *prev, double *pz);
static void Get_previous_nuc(char *apre, struct bond *bond, char res);
static int Get_previous_AA(char *pre, struct bond *bond, int res);

void Trajectory(atom *atoms, int natoms, struct axe *axes, int naxe,
		float **Tors_mode, int imode, float fact, int nmove,
		struct bond *bonds)
{
  int iter, i, iatom, j;

  Set_bonds_measure(bonds, natoms, atoms);
  double d_phi[naxe];
  for(i=0; i<naxe; i++)d_phi[i]=Tors_mode[imode][i]*fact;

  atom *atommove=malloc(natoms*sizeof(atom));
  for(iter=0; iter<nmove; iter++){
    Build_up(bonds, natoms, d_phi, naxe);
    double r2=0; float x;
    for(i=0; i<natoms; i++){
      struct bond *bond=bonds+i;
      if(bond->i_atom<0)break;
      iatom=bond->i_atom;
      for(j=0; j<3; j++){
	atommove[i].r[j]=bonds[i].r[j];
	x=atommove[i].r[j]-atoms[iatom].r[j];
	r2+=x*x;
      }
      // printf("%3d %s %.4f\n", i, bond->name, x);
    }
    printf("mode %2d step %3d d2=%.4f\n", imode, iter, r2);
  }
  free(atommove);
}

struct bond *Set_bonds_topology(int N_atoms, atom *atoms, struct residue *seq)
{
  struct bond *bonds=malloc((N_atoms+1)*sizeof(struct bond));
  BONDTHR1_2=BONDTHR1*BONDTHR1;

  if(DEBUG)printf("Setting bond topology\n");
  int i;
  char chain=atoms[0].chain;
  atom *atom=atoms;
  struct bond *bond=bonds, *previous;
  for(i=0; i<N_atoms; i++){
    //ALIGN
    if(atom->ali<0)continue; // for 2 str, not aligned atoms are omitted
    atom->main=Is_atom_main(atom->name, seq->type);
    bond->i_atom=i;
    bond->atom=atom;
    for(int j=0; j<3; j++)bond->r[j]=atom->r[j];
    bond->terminal=1;
    bond->d_phi=0; bond->len=-1; bond->angle=-1; bond->torsion=-1;
    // The first atom of a chain is joined to the third one
    // of the first chain (n=3)
    if(atom->chain!=chain){
      chain=atom->chain;
      previous=bonds+2;
    }else if(i){
      previous=Find_previous(bond, i, N_atoms, seq+bond->atom->res);
      if(previous==NULL){
	printf("ERROR, atom %d %s res %d has no previous atom\n",
	       i, atom->name, atom->res); exit(8);
      }
    }else{previous=NULL;}
    bond->previous=previous;
    if(previous)previous->terminal=0;
    bond++; atom++;
  }
  bonds[i].atom=NULL; // Last bond
  printf("%d atoms, %d bonds for building the molecule\n", N_atoms, i);
  return(bonds);
}

int Set_bonds_measure(struct bond *bonds, int N_atoms, atom *atoms)
{
  int i, j; atom *atom=atoms;
  struct bond *bond=bonds, *previous;
  //~ struct bond *bond2;

  if(DEBUG)printf("Setting bond lengths and angles\n");

  bonds[0].previous=NULL;
  for(i=0; i<N_atoms; i++){
    bond=bonds+i;
    if(bond->i_atom<0)break;
    atom=atoms+bond->i_atom;
    for(j=0; j<3; j++)bond->r[j]=atom->r[j];
    previous=bond->previous;

    if(previous==NULL){ // First atom
      if(i){
	printf("ERROR, atom %d %s res %d has no previous atom\n",
	       i, atom->name, atom->res); exit(8);
      }
      double d=0;
      for(j=0; j<3; j++){
	bond->dx[j]=0; bond->dy[j]=0;
	float dz=(atom+1)->r[j]-atom->r[j];
	bond->dz[j]=dz; d+=dz*dz;
      }
      d=sqrt(d); for(j=0; j<3; j++)bond->dz[j]/=d;
      continue;
    }

    // Set r
    double dr_dr=0;
    for(j=0; j<3; j++){
      double dr=bond->r[j]-previous->r[j];
      bond->dr[j]=dr; dr_dr+=dr*dr;
    }
    bond->l=sqrt(dr_dr);
    for(j=0; j<3; j++){
      bond->dz[j]=bond->dr[j]/bond->l;
    }
    // Initialize dx
    for(j=0; j<3; j++)bond->dx[j]=0;

    if(previous->previous==NULL){
      bond->r_cos_theta=-bond->l; bond->phi=0;
      double dx[3];
      for(j=0; j<3; j++)dx[j]=-(atom+1)->r[j]+atom->r[j];
      if(Make_frame(bond, NULL, dx)<0){
	printf("ERROR in Make_frame\n"); exit(8);
      }
      continue;
    }

    // i>0,1 set angle and torsion

    // Make orthogonal frame for next atom
    if(bond->terminal==0){
      if(Make_frame(bond, previous, previous->dz)<0){
	printf("ERROR in Make_frame\n"); exit(8);
      }
    }

    // Set theta
    bond->r_cos_theta=-Scalar_product_3_d(bond->dr, previous->dz);
    if(bond->r_cos_theta > bond->l){bond->r_cos_theta = bond->l;}
    else if(bond->r_cos_theta < -bond->l){bond->r_cos_theta = -bond->l;}
    double r_sin_theta= sqrt(dr_dr-bond->r_cos_theta*bond->r_cos_theta);

    // Set phi
    double ax=Scalar_product_3_d(bond->dr, previous->dx);
    double cos_phi=ax/r_sin_theta;
    if(cos_phi>1){cos_phi=1;}
    else if(cos_phi<-1){cos_phi=-1;}
    double ay=Scalar_product_3_d(bond->dr, previous->dy);
    if(ay >= 0){
      bond->phi=acos(cos_phi);
    }else{
      bond->phi=-acos(cos_phi);
    }

    if(DEBUG){
      double dr[3], dx[3];
      double az=-bond->r_cos_theta;
      for(j=0; j<3; j++){
	dr[j]=
	  az*previous->dz[j]+
	  ax*previous->dx[j]+
	  ay*previous->dy[j];
      }
      dx[0]=Scalar_product_3_d(bond->dr, previous->dx)-
	Scalar_product_3_d(dr, previous->dx);
      dx[1]=Scalar_product_3_d(bond->dr, previous->dy)-
	Scalar_product_3_d(dr, previous->dy);
      dx[2]=Scalar_product_3_d(bond->dr, previous->dz)-
	Scalar_product_3_d(dr, previous->dz);
      printf("%3d %4s %4s %8.2g %8.2g %8.2g\n", bond->i_atom,
	     atom->name, atoms[bond->previous->i_atom].name,
	     dx[0], dx[1], dx[2]);
    }
  }
  return(0);
}

void Build_up(struct bond *bonds, int N_atoms, double *d_phi, int N_axes)
{
  float COSMAX=0.999; // Maximum value of cos(theta) to avoid collinear bonds
  struct bond *bond, *previous; int i, j;

  for(i=1; i<N_atoms; i++){
    bond=bonds+i;
    if(bond->i_atom<0)break;
    previous=bond->previous;
    if(previous==NULL){
      printf("ERROR, build up can not be applied at bond %d\n", i);
      return;
    }
    if((previous->previous==NULL)||
     (previous->previous->previous==NULL))continue;
    if((d_phi)&&(bond->len >= 0)){
      double factor=1.+d_phi[bond->len]/bond->l; //
      if(factor <= 0.1){
	printf("WARNING, bond length too short with buildup! ");
	printf("atoms: %s%d%c  %s%d%c f= %.2g\n",
	       bond->atom->name, bond->atom->res, bond->atom->aa,
	       bond->previous->atom->name, bond->previous->atom->res,
	       bond->previous->atom->aa, factor);
	factor=0.1;
	printf("Setting shrink factor to %.2g\n", factor);
      }
      bond->l*=factor;
      bond->r_cos_theta*=factor;
      double lcosmax=bond->l*COSMAX;
      if(bond->r_cos_theta>lcosmax){bond->r_cos_theta=lcosmax;}
      else if(bond->r_cos_theta<-lcosmax){bond->r_cos_theta=-lcosmax;}
    }

    double r_sin_theta=
      sqrt(bond->l*bond->l-bond->r_cos_theta*bond->r_cos_theta);
    if((d_phi)&&(bond->angle >= 0)){
      // Cos(theta)->Cos(theta+d_phi)
      double cos_a=cos(d_phi[bond->angle]);
      double sin_a=sin(d_phi[bond->angle]);
      bond->r_cos_theta=bond->r_cos_theta*cos_a-r_sin_theta*sin_a;
      double lcosmax=bond->l*COSMAX;
      if(bond->r_cos_theta>lcosmax){bond->r_cos_theta=lcosmax;}
      else if(bond->r_cos_theta<-lcosmax){bond->r_cos_theta=-lcosmax;}
      r_sin_theta=
	sqrt(bond->l*bond->l-bond->r_cos_theta*bond->r_cos_theta);
    }

    if((d_phi)&&(bond->torsion >= 0)){
      bond->phi += d_phi[bond->torsion];
      while(bond->phi > TWOPI){bond->phi-=TWOPI;}
      while(bond->phi <-TWOPI){bond->phi+=TWOPI;}
    }
    double ax= r_sin_theta*cos(bond->phi);
    double ay= r_sin_theta*sin(bond->phi);
    double az= -bond->r_cos_theta;

    for(j=0; j<3; j++){
      double dr=
	az*previous->dz[j]+
	ax*previous->dx[j]+
	ay*previous->dy[j];
      bond->dr[j]=dr;
      bond->r[j]=previous->r[j]+dr;
    }
    if(isnan(bond->r[0])||isnan(bond->r[1])||isnan(bond->r[2])){
      printf("ERROR nan, atom %s%d ", bond->atom->name, i);
      printf("len: %d ", bond->len);
      if(bond->len>=0)printf(" %.2g", d_phi[bond->len]);
      printf("  angle: %d ", bond->angle);
      if(bond->angle>=0){
	printf(" %.3f cos=%.3f",d_phi[bond->angle], cos(d_phi[bond->angle]));
      }
      printf("  torsion: %d ", bond->torsion);
      if(bond->torsion>=0)printf(" %.3f", d_phi[bond->torsion]);
      printf("\n");
      printf("r= %.3f", bond->l);
      printf(" r*cos(theta)= %.3f", bond->r_cos_theta);
      printf(" r*sin(theta)= %.3f", r_sin_theta);
      printf(" phi= %.3f", bond->phi);
      double *v=previous->dx; printf(" x: %.3f %.3f %.3f",v[0],v[1],v[2]);
      v=previous->dy; printf(" y: %.3f %.3f %.3f",v[0],v[1],v[2]);
      v=previous->dz; printf(" z: %.3f %.3f %.3f",v[0],v[1],v[2]);
      printf("\n"); 
      //for(i=0; i<N_axes; i++)printf("%.2f ",d_phi[i]); printf("\n"); 
      //exit(8);
    }
    double dr_dr=Scalar_product_3_d(bond->dr, bond->dr);
    bond->l=sqrt(dr_dr);
    for(j=0; j<3; j++){
      bond->dz[j]=bond->dr[j]/bond->l;
    }

    // Make new orthogonal frame
    if(bond->terminal==0 && Make_frame(bond, previous, previous->dz)<0){
      printf("Could not make last frame, exiting\n");
      return;
    }
  }
  return;
}

void Compare_build_up(struct bond *bonds, int N_atoms,
		      double *d_phi, int N_axes,
		      struct bond *bonds2)
{
  struct bond *bond, *previous; int i, j, error=0;
  double dr_dr;

  printf("Testing differences atom by atom\n");

  for(i=0; i<N_atoms; i++){
    bond=bonds+i;
    if(bond->i_atom<0)break;
    previous=bond->previous;
    if((previous==NULL)||(previous->previous==NULL)||
       (previous->previous->previous==NULL)){
      if(bonds2){
	for(j=0; j<3; j++){
	  bond->r[j]=bonds2[i].r[j];
	  bond->dr[j]=bonds2[i].dr[j];
	  bond->dz[j]=bonds2[i].dz[j];
	}
      }
      if(previous){
	goto test;
      }else{
	continue;
      }
    }
    if((d_phi)&&(bond->len >= 0)){
      double factor=1.+d_phi[bond->len]/bond->l; //
      if(factor <= 0.1){
	printf("WARNING, bond length too short with buildup! ");
	printf("atoms: %s%d%c  %s%d%c f= %.2g\n",
	       bond->atom->name, bond->atom->res, bond->atom->aa,
	       bond->previous->atom->name, bond->previous->atom->res,
	       bond->previous->atom->aa, factor);
	factor=0.1;
	printf("Setting shrink factor to %.2g\n", factor);
      }
      bond->l*=factor;
      bond->r_cos_theta*=factor;
    }
    double r_sin_theta=
      sqrt(bond->l*bond->l-bond->r_cos_theta*bond->r_cos_theta);
    if((d_phi)&&(bond->angle >= 0)){
      double cos_a=cos(d_phi[bond->angle]);
      double sin_a=-sin(d_phi[bond->angle]);
      bond->r_cos_theta=bond->r_cos_theta*cos_a-r_sin_theta*sin_a;
      r_sin_theta=
	sqrt(bond->l*bond->l-bond->r_cos_theta*bond->r_cos_theta);
    }
    if((d_phi)&&(bond->torsion >= 0)){
      bond->phi += d_phi[bond->torsion];
    }
    double ax= r_sin_theta*cos(bond->phi);
    double ay= r_sin_theta*sin(bond->phi);
    double az= -bond->r_cos_theta;

    for(j=0; j<3; j++){
      double dr=
	az*previous->dz[j]+
	ax*previous->dx[j]+
	ay*previous->dy[j];
      bond->dr[j]=dr;
      bond->r[j]=previous->r[j]+dr;
    }

  test:
    
    dr_dr=Scalar_product_3_d(bond->dr, bond->dr);
    bond->l=sqrt(dr_dr);
    for(j=0; j<3; j++){
      bond->dz[j]=bond->dr[j]/bond->l;
    }

    // Make new orthogonal frame
    if(bond->terminal==0)Make_frame(bond, previous, previous->dz);

    if(bond->atom->main==0)continue;
    if(bonds2){
      double d2=0;
      for(j=0; j<3; j++){
	double d=bond->r[j]-bonds2[i].r[j]; d2+=d*d;
      }
      if(d2>0.005){
	printf("ERROR in Compare_build_up, bond %d atom %s %c%d chain %d",
	       i, bond->atom->name, bond->atom->aa, bond->atom->res,
	       bond->atom->chain);
	printf("d=%.4g\n", sqrt(d2));
	if(bond->len>=0)
	  printf("Len: %.2f %d\n", d_phi[bond->len], bond->len);
	if(bond->angle>=0)
	  printf("Angle: %.2f %d\n", d_phi[bond->angle], bond->angle);
	if(bond->torsion>=0)
	  printf("Torsion: %.2f %d\n", d_phi[bond->torsion], bond->torsion);
	error++; exit(8);
      }
    }

  }
  if(error)exit(8);
  return;
}

void Internal_coordinates(double *int_coord, int N_axes,
			  struct bond *bonds, int N_atoms)
{
  struct bond *bond; int i, k=0;
  for(i=0; i<N_axes; i++)int_coord[i]=0;
  for(i=1; i<N_atoms; i++){
    bond=bonds+i;
    if(bond->i_atom<0)break;
    if(bond->len >= 0){
      if(bond->len>=N_axes)goto error;
      int_coord[bond->len]=bond->l; k++;
      if(isnan(int_coord[bond->len]))goto error;
    }
    if(bond->angle >= 0){
      if(bond->angle>=N_axes)goto error;
      double c=bond->r_cos_theta/bond->l;
      if(c>1){c=1;}else if(c<-1){c=-1;}
      int_coord[bond->angle]=acos(c); k++;
      if(isnan(int_coord[bond->angle]))goto error;
    }
    if((bond->torsion >= 0)&&(bond->atom->main)){
      if(bond->torsion>=N_axes)goto error;
      double c=bond->phi;
      while(c > TWOPI){c-=TWOPI;}
      while(c <-TWOPI){c+=TWOPI;}
      int_coord[bond->torsion]=c; k++;
      if(isnan(int_coord[bond->torsion]))goto error;
    }
  }
  if(0)printf("%d internal coordinates over %d stored\n", k, N_axes);
  return;
 error:
  printf("ERROR in internal_coordinates, bond %d atom %s %c%d chain %d\n",
	 i, bond->atom->name, bond->atom->aa,
	 bond->atom->res, bond->atom->chain);
  printf("len: %d angle: %d torsion: %d Total_dofs= %d\n",
	 bond->len, bond->angle, bond->torsion, N_axes);
  if(bond->len>=0)printf("Len: %.2f\n", int_coord[bond->len]);
  if(bond->angle>=0)printf("Angle: %.2f\n", int_coord[bond->angle]);
  if(bond->torsion>=0)printf("Torsion: %.2f\n", int_coord[bond->torsion]);
  printf("r= %.2f r_cos_theta %.4f phi= %.4f\n",
	  bond->l, bond->r_cos_theta, bond->phi);
  exit(8);
}


int Make_frame(struct bond *bond, struct bond *previous, double *pz)
{
  double dr_dr1=
    bond->dz[0]*pz[0]+bond->dz[1]*pz[1]+bond->dz[2]*pz[2];

  double norm=0; int j;
  for(j=0; j<3; j++){
    double dx=-pz[j]+dr_dr1*bond->dz[j]; // Standard definition
    //double dx=pz[j]-dr_dr1*bond->dz[j];
    bond->dx[j]=dx; norm+=dx*dx;
  }
  if(norm>0){ 
    norm=1./sqrt(norm);
    for(j=0; j<3; j++)bond->dx[j]*=norm;
  }else{  // Collinear bonds, theta=0 or pi. Use previous axis
    //double dx3[3]; for(j=0; j<3; j++)dx3[j]=bond->dx[j];
    if(previous==NULL){printf("ERROR no previous bond\n"); exit(8);}
    for(j=0; j<3; j++)bond->dx[j]=previous->dx[j];
  } 
  
  Vector_product_d(bond->dy, bond->dz, bond->dx);
  norm=0; for(j=0; j<3; j++)norm+=bond->dy[j]*bond->dy[j];
  if(norm>0){
    norm= 1./sqrt(norm);
    for(j=0; j<3; j++)bond->dy[j]*=norm;
  }else{
    printf("ERROR in Make_frame norm_y= %.2g atom= %s res=%d\n"
	   "dy=%.2g %.2g %.2g dx= %.2g %.2g %.2g dz= %.2g %.2g %.2g\n",
	   norm, bond->atom->name, bond->atom->res,
	   bond->dy[0], bond->dy[1], bond->dy[2],
	   bond->dx[0], bond->dx[1], bond->dx[2],
	   bond->dz[0], bond->dz[1], bond->dz[2]);
    return(-1);
  }
  return(0);
}

struct bond *Find_previous(struct bond *bond, int n, int N_atoms,
			   struct residue *seq)
{
  // Find previous atom in same chain
  struct bond *previous=bond-1;
  int res=previous->atom->res;
  if((res+1)<bond->atom->res){
    printf("Gap %d-%d ", res, bond->atom->res);
  }

  int ichar=-1; char name_p[4]="   "; // Ex. B in CB
  if((seq->type==1)||(seq->type==-1)){ // Amino-acid
    ichar=Get_previous_AA(name_p, bond, res);
  }else if((seq->type==2)||(seq->type==3)||
	   (seq->type==-2)||(seq->type==-3)){ // Nucleic acid
    Get_previous_nuc(name_p, bond, seq->amm);
  }

  struct bond *prev=NULL, *pname=NULL;
  double dmin=1000000, dname=100000; int i;
  for(i=n; i>0; i--){
    if(previous->atom->res<(res-1))break;
    if((bond->atom->main)&&(previous->atom->main==0))goto next;
    if((((ichar>=0)&&(previous->atom->name[ichar]==name_p[0]))||
	(strncmp(previous->atom->name, name_p, 3)==0)||
	(strncmp(previous->atom->name, name_p, 2)==0))&&
       (previous->atom->name[0]!='H')){
      double d2=0, x; int j;
      for(j=0; j<3; j++){
	x=previous->r[j]-bond->r[j]; d2+=x*x;
      }
      if(d2< BONDTHR1_2)return(previous);
      if(d2<dname){pname=previous; dname=d2;}
    }else{ // Not the name I am looking for!
      double d2=0, x; int j;
      for(j=0; j<3; j++){
	x=previous->r[j]-bond->r[j]; d2+=x*x;
      }
      if(d2<dmin){prev=previous; dmin=d2;}
    }
  next:
    previous--;
  }
  /*// Look ahead
  previous=bond+1;
  for(i=n+1; i<N_atoms; i++){
    if(previous->atom->res>res)break;
    if((bond->atom->main)&&(previous->atom->main==0))goto next2;
    if((((ichar>=0)&&(previous->atom->name[ichar]==name_p[0]))||
	(strncmp(previous->atom->name, name_p, 3)==0)||
	(strncmp(previous->atom->name, name_p, 2)==0))&&
       (previous->atom->name[0]!='H')){
      float d2=0, x; int j;
      for(j=0; j<3; j++){
	x=previous->r[j]-bond->r[j]; d2+=x*x;
      }
      if(d2< BONDTHR1_2)return(previous);
      if(d2<dname){pname=previous; dname=d2;}
    }else{ // Not the name I am looking for!
      float d2=0, x; int j;
      for(j=0; j<3; j++){
	x=previous->r[j]-bond->r[j]; d2+=x*x;
      }
      if(d2<dmin){prev=previous; dmin=d2;}
    }
  next2:
    previous++;
    }*/


  if(pname!=NULL){
    prev=pname;
  }else{
    printf("WARNING, no proper previous atom of %s ", bond->atom->name);
    printf("res %d chain %d type %d %c\n", bond->atom->res,
	   bond->atom->chain, seq->type, seq->amm);
    printf("Expected previous: %s %d %c\n", name_p, ichar, name_p[0]);
  }

  if(prev==NULL){ 
    printf("ERROR, no previous atom of %s found\n", bond->atom->name);
    printf("n= %d i= %d ", n, i);
    printf("d(atom: %s%d%c previous: %s %d%c)= %.1f\n",
	 bond->atom->name, bond->atom->res, bond->atom->aa,
	 previous->atom->name, previous->atom->res, previous->atom->aa,
	   sqrt(dmin));
    exit(8);
  }
  printf("d(atom: %s%d%c previous: %s%d%c)= %.1f\n",
	 bond->atom->name, bond->atom->res, bond->atom->aa,
	 prev->atom->name, prev->atom->res, prev->atom->aa, sqrt(dmin));
  return(prev);
}

int Get_previous_AA(char *pre, struct bond *bond, int res){
  // returns ichar position of character that should match *pre
  if(res<bond->atom->res){
    strcpy(pre, "C "); return(-1);
  }
  char *name=bond->atom->name, index=name[1];
  if(name[0]=='H'){
    if((index=='N')||(index=='1')||(index=='2')||(index=='3')){*pre=' ';}
    else{*pre=index;}
  }else if(index==' '){
    if(name[0]=='N'){
      strcpy(pre, "C "); return(-1);
    }else if(name[0]=='C'){
      strcpy(pre, "CA"); return(-1); // CA
    }else if(name[0]=='O'){
      strcpy(pre, "C "); return(-1); // C
    }
  }else if(index=='A'){
    strcpy(pre, "N "); return(-1);
  }else if(index=='B'){
    strcpy(pre, "CA"); return(-1);
  }else if(index=='G'){
    *pre='B';
  }else if(index=='D'){
    *pre='G';
  }else if(index=='E'){
    *pre='D';
  }else if(index=='Z'){
    *pre='E';
  }else if(index=='H'){
    *pre='Z';
  }else if((strncmp(name,"OXT",3)==0)||
	   (strncmp(name,"OT" ,2)==0)){
    strcpy(pre, "C "); return(-1);
  }else{
    printf("WARNING, unknown atom name: %s res %d chain %d\n",
	   name, bond->atom->res, bond->atom->chain);
    strcpy(pre, "  "); return(-1);
    //exit(8);
  }
  return(1);
}

void Get_previous_nuc(char *apre, struct bond *bond, char res){
  // Backbone
  char *name=bond->atom->name;
  if(strncmp(name,"P  ",3)==0){
    strcpy(apre, "O3'");
  }else if(strncmp(name,"OP1",3)==0){
    strcpy(apre, "P  ");
  }else if(strncmp(name,"OP2",3)==0){
    strcpy(apre, "P  ");
  }else if(strncmp(name,"O5'",3)==0){
    strcpy(apre, "P  ");
  }else if(strncmp(name,"C5'",3)==0){
    strcpy(apre, "O5'");
  }else if(strncmp(name,"C4'",3)==0){
    strcpy(apre, "C5'");
  }else if(strncmp(name,"O4'",3)==0){
    strcpy(apre, "C4'");
  }else if(strncmp(name,"C3'",3)==0){
    strcpy(apre, "C4'");
  }else if(strncmp(name,"O3'",3)==0){
    strcpy(apre, "C3'");
  }else if(strncmp(name,"C2'",3)==0){
    strcpy(apre, "C3'");
  }else if(strncmp(name,"O2'",3)==0){ // RNA!
    strcpy(apre, "C2'");
  }else if(strncmp(name,"C1'",3)==0){
    strcpy(apre, "C2'");
    // Bases: Purines T,C
  }else if((res=='T')||(res=='U')||(res=='C')){
    if(strncmp(name,"N1 ",3)==0){
      strcpy(apre, "C1'");
    }else if(strncmp(name,"C2 ",3)==0){
      strcpy(apre, "N1 ");
    }else if(strncmp(name,"O2 ",3)==0){
      strcpy(apre, "C2 ");
    }else if(strncmp(name,"N3 ",3)==0){
      strcpy(apre, "C2 ");
    }else if(strncmp(name,"C4 ",3)==0){
      strcpy(apre, "N3 ");
    }else if(strncmp(name,"O4 ",3)==0){
      strcpy(apre, "C4 ");
    }else if(strncmp(name,"N4 ",3)==0){
      strcpy(apre, "C4 ");
    }else if(strncmp(name,"C5 ",3)==0){
      strcpy(apre, "C4 ");
    }else if(strncmp(name,"C7 ",3)==0){
      strcpy(apre, "C5 ");
    }else if(strncmp(name,"C6 ",3)==0){
      strcpy(apre, "N1 ");
    }
    // Pyrimidine
  }else if((res=='G')||(res=='A')){
    if(strncmp(name,"N9",2)==0){
      strcpy(apre, "C1'");
    }else if(strncmp(name,"C8",2)==0){
      strcpy(apre, "N9");
    }else if(strncmp(name,"N7",2)==0){
      strcpy(apre, "C8");
    }else if(strncmp(name,"C5",2)==0){
      strcpy(apre, "N7");
    }else if(strncmp(name,"C6",2)==0){
      strcpy(apre, "C5");
    }else if(strncmp(name,"O6",2)==0){
      strcpy(apre, "C6");
    }else if(strncmp(name,"N6",2)==0){
      strcpy(apre, "C6");
    }else if(strncmp(name,"N1",2)==0){
      strcpy(apre, "C6");
    }else if(strncmp(name,"C2",2)==0){
      strcpy(apre, "N1");
    }else if(strncmp(name,"N2",2)==0){
      strcpy(apre, "C2");
    }else if(strncmp(name,"N3",2)==0){
      strcpy(apre, "C2");
    }else if(strncmp(name,"C4",2)==0){
      strcpy(apre, "N9");
    }
  }else{
    printf("WARNING, unknown atom name: %s res %d chain %d\n",
	   bond->atom->name, bond->atom->res, bond->atom->chain);
  }
}

int Is_atom_main(char *name, int seq_type){
  if((seq_type==1)||(seq_type==-1)){ // Amino-acid
    if((strncmp(name, "N ", 2)==0)||
       (strncmp(name, "CA", 2)==0)||
       (strncmp(name, "C ", 2)==0)){return(1);}
    else{return(0);}
  }else if((seq_type== 2)||(seq_type== 3)||
	   (seq_type==-2)||(seq_type==-3)){ // Nucleic acid
    if((strncmp(name,"P ", 3)==0)||
       (strncmp(name,"O5'", 3)==0)||
       (strncmp(name,"C5'", 3)==0)||
       (strncmp(name,"C4'", 3)==0)||
       (strncmp(name,"C3'", 3)==0)||
       (strncmp(name,"O3'", 3)==0)){return(1);}
    else{return(0);}
  }else{
    return(1);
  }
}
    


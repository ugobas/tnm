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

struct bond *Find_previous(struct bond *bond, int n, int N_atoms,
			   struct residue *seq);
void Make_frame(struct bond *bond, float *ddz);
static void Get_previous_nuc(char *apre, struct bond *bond, char res);
static int Get_previous_AA(char *pre, struct bond *bond, int res);

void Trajectory(atom *atoms, int natoms, struct axe *axes, int naxe,
		float **Tors_mode, int imode, float fact, int nmove,
		struct bond *bonds)
{
  int iter, i, iatom, j;

  Set_bonds_measure(bonds, natoms, atoms);
  float *d_phi=malloc(naxe*sizeof(float));
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
  free(d_phi);
  free(atommove);
}

struct bond *Set_bonds_topology(int N_atoms, atom *atoms, struct residue *seq)
{
  if(DEBUG)printf("Setting bond topology\n");
  struct bond *bonds=malloc((N_atoms+1)*sizeof(struct bond));

  int achain=1, n=0, lastchain=-1, i, j;
  char chain=' '; atom *atom=atoms;
  bonds[0].previous=NULL;
  BONDTHR1_2=BONDTHR1*BONDTHR1;
  for(i=0; i<N_atoms; i++){
    struct bond *bond=bonds+n, *bond1=NULL;
    bond->i_atom=i;
    atom=atoms+i;
    bond->atom=atom;
    atom->main=Is_atom_main(atom->name, seq->type);
    for(j=0; j<3; j++)bond->r[j]=atom->r[j];
    bond->terminal=1;

    // The first atom of a chain is joined to the third one
    // of the previous chain 
    if(achain<4)achain++;
    if(achain==3)lastchain=n; 
    if(atom->chain!=chain){
      chain=atom->chain; achain=1;
      bond1=bonds+lastchain;
    }else if(n){
      bond1=Find_previous(bond, n, N_atoms, seq+bond->atom->res);
      if(bond1==NULL){
	printf("ERROR, atom %d %s res %d has no previous atom\n",
	       i, atom->name, atom->res); exit(8);
      }
    }
    n++; if(n==1)continue;
    bond->previous=bond1;
    bond1->terminal=0;
  }
  bonds[n].atom=NULL; // Last bond
  printf("%d bonds for building the molecule\n", n);
  return(bonds);
}

int Set_bonds_measure(struct bond *bonds, int N_atoms, atom *atoms)
{
  int i, j; atom *atom=atoms;
  struct bond *bond=bonds, *previous;

  if(DEBUG)printf("Setting bond lengths and angles\n");

  bonds[0].previous=NULL;
  for(i=0; i<N_atoms; i++){
    bond=bonds+i;
    if(bond->i_atom<0)break;
    atom=atoms+bond->i_atom;
    for(j=0; j<3; j++)bond->r[j]=atom->r[j];
    previous=bond->previous;

    if(previous==NULL){
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
    double d2=0;
    for(j=0; j<3; j++){
      double dr=bond->r[j]-previous->r[j];
      bond->dr[j]=dr; d2+=dr*dr;
    }
    bond->l=sqrt(d2);
    for(j=0; j<3; j++){
      bond->dz[j]=bond->dr[j]/bond->l;
    }

    if(previous->previous==NULL){
      bond->r_cos_theta=bond->l; bond->phi=0;
      double d=0;
      for(j=0; j<3; j++){
	float dx=(atom+1)->r[j]-atom->r[j];
	bond->dx[j]=dx; d+=dx*dx;
      }
      d=sqrt(d); for(j=0; j<3; j++)bond->dx[j]/=d;
      Vector_product(bond->dy, bond->dz, bond->dx);
      d=Scalar_product_3(bond->dy, bond->dy);
      d= sqrt(d);
      for(j=0; j<3; j++)bond->dy[j]/=d;
      continue;
    }

    // i>0,1 set angle and torsion
    // Make orthogonal frame for next atom
    if(bond->terminal==0)Make_frame(bond, previous->dz);

    // Set theta
    bond->r_cos_theta=Scalar_product_3(bond->dr, previous->dz);
    if(bond->r_cos_theta > bond->l){bond->r_cos_theta = bond->l;}
    else if(bond->r_cos_theta < -bond->l){bond->r_cos_theta = -bond->l;}
    float r_sin_theta= sqrt(d2-bond->r_cos_theta*bond->r_cos_theta);

    // Set phi
    if(previous->previous->previous==NULL)continue;
    float ax=Scalar_product_3(bond->dr, previous->dx);
    float cos_phi=ax/r_sin_theta;
    if(fabs(cos_phi)>1){
      //printf("WARNING cos_phi= %.3f\n", cos_phi);
      if(cos_phi>1){cos_phi=1;}
      else{cos_phi=-1;}
    }
    float ay=Scalar_product_3(bond->dr, previous->dy);
    if(ay >= 0){
      bond->phi=acos(cos_phi);
    }else{
      bond->phi=-acos(cos_phi);
    }

    if(DEBUG){
      float dr[3], dx[3];
      float az=bond->r_cos_theta;
      for(j=0; j<3; j++){
	dr[j]=
	  az*previous->dz[j]+
	  ax*previous->dx[j]+
	  ay*previous->dy[j];
      }
      dx[0]=Scalar_product_3(bond->dr, previous->dx)-
	Scalar_product_3(dr, previous->dx);
      dx[1]=Scalar_product_3(bond->dr, previous->dy)-
	Scalar_product_3(dr, previous->dy);
      dx[2]=Scalar_product_3(bond->dr, previous->dz)-
	Scalar_product_3(dr, previous->dz);
      printf("%3d %4s %4s %8.2g %8.2g %8.2g\n", bond->i_atom,
	     atom->name, atoms[previous->i_atom].name,
	     dx[0], dx[1], dx[2]);
    }
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
  float dmin=1000000, dname=100000; int i;
  for(i=n; i>0; i--){
    if(previous->atom->res<res)break;
    if((bond->atom->main)&&(previous->atom->main==0))goto next;
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
  next:
    previous--;
  }

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
    exit(8);
  }
  printf("d(atom: %s%d%c previous: %s%d%c)= %.1f\n",
	 bond->atom->name, bond->atom->res, bond->atom->aa,
	 prev->atom->name, prev->atom->res, prev->atom->aa, sqrt(dmin));
  return(prev);
}


void Build_up(struct bond *bonds, int N_atoms, float *d_phi, int N_axes)
{
  struct bond *bond, *previous; int i, j;
  double dr;

  for(i=1; i<N_atoms; i++){
    bond=bonds+i;
    if(bond->i_atom<0)break;
    previous=bond->previous;
    if(previous==NULL){
      printf("ERROR, build up cannot be applied at bond %d\n", i);
      exit(8);
    }
    if((previous->previous==NULL)||
     (previous->previous->previous==NULL))continue;
    if((d_phi)&&(bond->len >= 0)){
      double factor=1.+d_phi[bond->len]/bond->l;
      if(factor <= 0.5){
	printf("WARNING, bond length too short with buildup! ");
	printf("atoms: %s%d%c  %s%d%c f= %.2g\n",
	       bond->atom->name, bond->atom->res, bond->atom->aa,
	       previous->atom->name, previous->atom->res,
	       previous->atom->aa, factor);
	factor=0.5;
	printf("Setting shrink factor to %.2g\n", factor);
      }
      bond->l*=factor;
      bond->r_cos_theta*=factor;
    }

    float r_sin_theta=
      sqrt(bond->l*bond->l-bond->r_cos_theta*bond->r_cos_theta);
    if((d_phi)&&(bond->angle >= 0)){
      double cos_a=cos(d_phi[bond->angle]);
      double sin_a=sin(d_phi[bond->angle]);
      double rcos_th=bond->r_cos_theta*cos_a-r_sin_theta*sin_a;
      bond->r_cos_theta=rcos_th;
      double rsin_th=bond->r_cos_theta*sin_a+r_sin_theta*cos_a;
      if(rsin_th < 0){r_sin_theta=-rsin_th;}
      else{r_sin_theta=rsin_th;}
    }
    if((d_phi)&&(bond->torsion >= 0)){
      //bond->phi += d_phi[bond->torsion]; change
      bond->phi -= d_phi[bond->torsion];
    }
    float ax= r_sin_theta*cos(bond->phi);
    float ay= r_sin_theta*sin(bond->phi);
    float az= bond->r_cos_theta;

    double dr_dr=0;
    for(j=0; j<3; j++){
      dr= az*previous->dz[j]+
	  ax*previous->dx[j]+
	  ay*previous->dy[j];
      bond->dr[j]=dr; dr_dr+=dr*dr;
      bond->r[j]=previous->r[j]+dr;
    }
    bond->l=sqrt(dr_dr);
    for(j=0; j<3; j++){
      bond->dz[j]=bond->dr[j]/bond->l;
    }

    // Make new orthogonal frame
    if(bond->terminal==0)Make_frame(bond, previous->dz);

  }

  return;
}

void Make_frame(struct bond *bond, float *pdz)
{
  // pdz= previous->dz;
  double dr_dr1=
    bond->dz[0]*pdz[0]+bond->dz[1]*pdz[1]+bond->dz[2]*pdz[2];

  double norm=0; int j;
  for(j=0; j<3; j++){
    double dx=-pdz[j]+dr_dr1*bond->dz[j]; // Standard definition
    //double dx=pdz[j]-dr_dr1*bond->dz[j];
    bond->dx[j]=dx; norm+=dx*dx;
  }
  norm=1./sqrt(norm);
  for(j=0; j<3; j++)bond->dx[j]*=norm;

  //Vector_product(bond->dy, bond->dz, pdz);
  Vector_product(bond->dy, bond->dz, bond->dx);
  norm=Scalar_product_3(bond->dy, bond->dy);
  norm= 1./sqrt(norm);
  for(j=0; j<3; j++)bond->dy[j]*=norm;
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
    printf("ERROR, unknown atom name: %s res %d chain %d\n",
	   name, bond->atom->res, bond->atom->chain);
    exit(8);
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
    

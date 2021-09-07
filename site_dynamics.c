#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "read.h"
#include "vector.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "allocate.h"
#include "diagonalize.h"
#include "binding_sites.h"

#define VERBOSE 0
extern void Print_sites(int **site_res, int *nres_site,
			char **description, int nsites, 
			struct residue *seq, char *pdb, char *chain);
extern int Extract_atoms(int *iref, int *iatom, int *atomres,
			 struct residue *seq, float *mass,
			 atom *atoms, struct Reference Ref, char *ANAME);
void Eigen_directions(float *d2, float *lambda, float **v,
		     int *site_res, int nres, char *chain_res,
		      struct Normal_Mode NM, int *atomres, int *iref);

int Find_residue(char *resnum, char chain,
		 struct residue *seq, int Nres, int i_ini);

int New_site(int *nsite,int *sres,int *nres_site,char **description,int NMAX);

//////////////////////////////////////////////////////////////////
void Binding_site_dynamics(struct Normal_Mode NM, struct Reference Ref,
			   struct residue *seq, atom *atoms,
			   char *nameout1, //int mode,
			   char *pdb, char *chain, int Nchain,
			   int Nres, char *SITES)
{
  int NMAX=1000; // Maximum number of sites
  int *site_res[NMAX], nres_site[NMAX];
  char *chain_res[NMAX], *description[NMAX], file_s[200];
  int PRINT_SITE=1;  // Print active site read in the PDB?
  printf("Reading binding sites\n");
  int N_sites=Read_sites(site_res, nres_site, chain_res, description, file_s,
			 NMAX, SITES, pdb, chain, Nchain, seq, Nres,PRINT_SITE);
  printf("%d binding sites read\n", N_sites);
  if(N_sites==0)return;

  // Extract c_alpha
  char SEL[4]="CA"; float mass=0;
  int Nref=Ref.N_ref, iref[Nref], iatom[Nref], atomres[Nref];
  int Na=Extract_atoms(iref, iatom, atomres, seq, &mass, atoms, Ref, SEL);
  if(Na!=Nres){
    printf("WARNING, %d residues but %d representative atoms", Nres, Na);
    printf(" found in Binding_site_dynamics\n");
  }

  char nameout[200]; int i, a, k;
  sprintf(nameout, "%s_sites_dynamics.dat", nameout1);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "# read binding sites in file %s\n", file_s);
  float lambda[3], d2[3];
  float **v=Allocate_mat2_f(3, 3);
  fprintf(file_out, "# 1=<d>^2 2=<d>^2_x 3=<d>^2_y 4=<d>^2_z "); k=5;
  for(a=0; a<3; a++){
    fprintf(file_out,
	    " %d=projection_%d %d=v%d_x %d=v%d_y %d=v%d_z",
	    k,a,k+1,a,k+2,a,k+3,a); k+=4;
  }
  fprintf(file_out, "\n");
  for(i=0; i<N_sites; i++){
    Eigen_directions(d2, lambda, v,
		     site_res[i], nres_site[i], chain_res[i], NM,
		     atomres, iref);
    fprintf(file_out, "# site %s\n", description[i]);
    float sum=d2[0]+d2[1]+d2[2];
    fprintf(file_out, "%.3g  %.3f %.3f %.3f  ",
	    sum, d2[0]/sum, d2[1]/sum, d2[2]/sum);
    for(a=0; a<3; a++){
      fprintf(file_out, "%.3g  %.3f %.3f %.3f    ",
	      lambda[a], v[a][0], v[a][1], v[a][2]);
    }
    fprintf(file_out, "\n");
  }
  printf("Writing %s\n", nameout);

  Empty_matrix_f(v,3); //modyves: v was never freed
  for(i=0; i<N_sites; i++){
    free(site_res[i]);
    free(chain_res[i]);
    if(description[i])free(description[i]);
  }
}

//////////////////////////////////////////////////////////////////
int Read_sites(int **site_res, int *nres_site, char **chain_res,
	       char **description, char *file_site,
	       int NMAX, char *SITES, char *pdb, char *chain, int Nchain,
	       struct residue *seq, int Nres, int PRINT_SITE)
{
  int SMAX=400, N_sites=0, i;
  for(i=0; i<NMAX; i++){
    site_res[i]=malloc(SMAX*sizeof(int));
    chain_res[i]=malloc(SMAX*sizeof(char));
    description[i]=NULL;
  }
  if(SITES[0]!='\0'){
    N_sites=Read_sites_file(site_res, chain_res, nres_site, description,
			    SITES, chain, seq, Nres, NMAX, SMAX);
    if(N_sites)strcpy(file_site, SITES);
  }
  if(N_sites==0){
    N_sites=Read_sites_PDB(site_res, chain_res, nres_site, description,
			   pdb, chain, seq, Nres, NMAX, SMAX);
    if(N_sites)strcpy(file_site, pdb);
    // Print site
    if(PRINT_SITE && N_sites)
      Print_sites(site_res, nres_site, description, N_sites,
		  seq, pdb, chain);
  }
  if(N_sites==0)return(0);
  for(i=0; i<N_sites; i++){
    if(description[i]==NULL){
      description[i]=malloc(10*sizeof(char));
      sprintf(description[i], "site%d", i+1);
    }else{
      char *ptr=description[i];
      while(*ptr!='\0'){
	if(*ptr=='\n'){*ptr='\0'; break;} ptr++;
      } 
    }
  }

  // Eliminate sites that are not in the selected chains 
  int S=0;
  for(i=0; i<N_sites; i++){
    if(nres_site[i]<=0)continue;
    int n=0;
    for(int j=0; j<nres_site[i]; j++){
      char chj=chain_res[i][j]; int nc;
      for(nc=0; nc<Nchain; nc++){
	if(chj==chain[nc]){n++; break;}
      }
    }
    if(n==nres_site[i]){
      if(S!=i){
	nres_site[S]=nres_site[i];
	strcpy(description[S], description[i]);
	for(int j=0; j<nres_site[i]; j++){
	  site_res[S][j]=site_res[i][j];
	  chain_res[S][j]=chain_res[i][j];
	}
      }
      S++;
    }
  }
  return(S);
}

//////////////////////////////////////////////////////////////////
int Read_sites_PDB(int **site_res, char **chain_res,
		   int *nres_site, char **description,
		   char *pdb, char *chain, struct residue *seq,
		   int Nres, int NMAX, int SMAX)
{
  int Compression=0;
  printf("Reading sites in file %s, chain %s %d residues\n",pdb, chain, Nres);
  FILE *file_in=Open_compressed_file(pdb, &Compression);
  if(file_in==NULL){
    printf("WARNING, file %s not found\n", pdb); return(0);
  }
  int nsites=-1, dsite=0, sres=0, l, i;
  for(i=0; i<NMAX; i++)nres_site[i]=0;
  char string[1000], resnum[6], resnam[5];
  char snam[5], snam_old[5]="   ", ch; i=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "SITE", 4)==0){
      sscanf(string+11, "%s", snam);
      if(strncmp(snam, snam_old, 3)!=0){
	if(nsites>=0)nres_site[nsites]=sres;
	strcpy(snam_old, snam); sres=0; nsites++;
	if(nsites >= NMAX){
	  printf("ERROR, too many sites > %d\n", NMAX); nsites--; break;
	}
      }
      for(l=18; l<60; l+=11){
	i=Read_site(resnum, resnam, &ch, string+l, chain);
	if(i==-2){break;}else if(i==-1){continue;}
	i=Find_residue(resnum, ch, seq, Nres, i);
	if(i<0){
	  printf("WARNING, residue %s_%c_%s not found\n",
		 resnum, ch, resnam); continue;
	}
	if(seq[i].amm!=Code_3_1(resnam)){
	  printf("WARNING, residue %d is %c while it is %s in site\n",
		 i, seq[i].amm, resnam);
	}
	site_res[nsites][sres]=i;
	chain_res[nsites][sres]=ch;
	sres++;
	if(sres > SMAX){
	  printf("ERROR, too many residues > %d\n", SMAX); goto end;
	}
      }
    }else if(strncmp(string+11, "SITE_DESCR", 10)==0){
      description[dsite]=malloc(100*sizeof(char));
      strcpy(description[dsite], string+28);
      dsite++;
    }else{
      if((nsites>=0)||(strncmp(string, "ATOM", 4)==0)){break;}
      else{continue;}
    }
  }
 end:
  fclose(file_in);
  if(Compression)Delete_tmp_file();
  if(nsites>=0)nres_site[nsites]=sres;
  nsites++;
  int num=0, msites=0;
  for(i=0; i<nsites; i++){
    num+=nres_site[i];
    if(nres_site[i]>0)msites++;
  }
  printf("%d residues read in %d binding sites out of %d\n",
	 num, msites, nsites);
  return(nsites);
}

//////////////////////////////////////////////////////////////////
int Read_sites_file(int **site_res, char **chain_res,
		    int *nres_site, char **description,
		    char *file, char *chain, struct residue *seq,
		    int Nres, int NMAX, int SMAX)
{
  printf("Reading sites in file %s, chain %s %d residues\n",file, chain, Nres);
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("WARNING, active site file %s not found\n", file); return(0);
  }
  int nsite=-1, sres=0, i, res; 
  char string[1000], resnum[6]="\0", resnam[5]="\0"; //, icode='\0';
  char snam[5], snam_old[5]="   ", ch[6]; i=-1;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      continue;
    }else if(strncmp(string, "SITE_DESCRIPTION", 10)==0){
      if(New_site(&nsite, &sres, nres_site, description, NMAX))break;
      strcpy(description[nsite], string+18);
      strcpy(snam_old, "   ");
    }else{
      sscanf(string, "%s%s%s%d", snam, resnam, ch, &res); //, &icode
      if(strcmp(snam, snam_old)!=0){
	if(strcmp(snam_old, "   ")!=0){
	  if(New_site(&nsite, &sres, nres_site, description, NMAX))break;
	  sprintf(description[nsite], "SITE_%s", snam);
	}
	strcpy(snam_old, snam);
      }
      //if(icode=='\n')icode='\0';
      //sprintf(resnum, "%d%c", res, icode);
      sprintf(resnum, "%d", res);
      i=Find_residue(resnum, ch[0], seq, Nres, i+1);
      if(i<0){
	printf("Residue %s_%s_%c discarded, chain: %s L=%d\n",
	       resnam, resnum,ch[0],chain,Nres); continue;
      }
      if(seq[i].amm!=Code_3_1(resnam)){
	printf("WARNING, residue %d is %c while it is %s in site %s\n",
	       i, seq[i].amm, resnam, snam);
      }
      site_res[nsite][sres]=i;
      chain_res[nsite][sres]=ch[0];
      sres++;
      if(sres > SMAX){
	printf("ERROR, too many residues > %d\n", SMAX); goto end;
      }
    }
  }
 end:
  fclose(file_in);
  if(nsite>=0)nres_site[nsite]=sres;
  nsite++;
  int num=0, msites=0;
  for(i=0; i<nsite; i++){
    num+=nres_site[i];
    if(nres_site[i]>0)msites++;
  }
  printf("file %s, reading %d residues in %d binding sites out of %d\n",
	 file, num, msites, nsite);
  return(nsite);
}

int Read_site(char *resnum, char *resnam, char *reschain,
	      char *string, char *chain){
  if(string[0]==' ')return(-2);
  char *c=chain; int inum;
  sscanf(string, "%s", resnam);
  *reschain=string[4];
  sscanf(string+5, "%d", &inum);
  //sprintf(resnum, "%4d%c", inum, string[9]); 
  sprintf(resnum, "%d", inum);
  if(string[9]!=' ')sprintf(resnum, "%s%c", resnum, string[9]);
  if(*c=='*')return(0);
  while(*c!='\0'){if(*c==*reschain)return(0); c++;}
  if(VERBOSE)
    printf("Residue %s %s_%c discarded, chain: %s\n",
	   resnam,resnum,*reschain,chain);
  return(-1);
}

int Find_residue(char *resnum, char chain,
		 struct residue *seq, int Nres, int i_ini)
{
  int i, ini=0; if(i_ini>0){ini=i_ini;}else{ini=0;}
  for(i=ini; i<Nres; i++){
    if((strcmp(seq[i].pdbres, resnum)==0)&&(seq[i].chain==chain))return(i);
  }
  for(i=0; i<ini; i++){
    if((strcmp(seq[i].pdbres, resnum)==0)&&(seq[i].chain==chain))return(i);
  }
  return(-1);
}

int Site_overlap(int *isite, int n1, int *jsite, int n2){
  int nc=0;
  for(int i=0; i<n1; i++){
    int i1=isite[i];
    for(int j=0; j<n2; j++){
      if(i1==jsite[j]){nc++; break;}
    }
  }
  return(nc);
}

void Eigen_directions(float *d2, float *lambda, float **v,
		      int *site_res, int nres, char *chain_res,
		      struct Normal_Mode NM, int *atomres, int *iref)
{
  /* Compute the matrix M_kl=sum_a 1/w^2_a x^a_ik x^a_jl
     Returns the eigenvalues lambda and eigenvectors v */

  double **matrix=Allocate_mat2_d(3,3);
  int atom_ref[nres], i, j, a, k, l;
  for(i=0; i<nres; i++){
    if(atomres[site_res[i]]>=0){
      atom_ref[i]=iref[atomres[site_res[i]]];
    }else{
      atom_ref[i]=-1;
    }
  }

  for(a=0; a<NM.N_relevant; a++){
    // For each normal mode
    if((NM.select[a]==0)||(NM.sigma2[a]==0))continue;
    float *Mode=NM.Cart[a];
    float w=1./(NM.omega2[a]);

    for(i=0; i<nres; i++){
      if(atom_ref[i]<0)continue;
      float *Mode_i=Mode+3*atom_ref[i];
      for(k=0; k<3; k++){
	float wx=w*Mode_i[k];
	for(j=0; j<nres; j++){
	  if(atom_ref[j]<0)continue;
	  float *Mode_j=Mode+3*atom_ref[j];
	  for(l=0; l<=k; l++){
	    matrix[k][l]+=wx*Mode_j[l];
	  }
	}
      }
    }
  }
  float norm=nres*nres;
  for(k=0; k<3; k++){
    d2[k]=matrix[k][k]/norm;
    for(l=k+1; l<3; l++)matrix[k][l]=matrix[l][k];
  }
  d_Diagonalize(3, matrix, lambda, v, 1);
  for(k=0; k<3; k++)lambda[k]/=norm;
  Empty_matrix_d(matrix, 3);
  return;
}

int New_site(int *nsite, int *sres, int *nres_site, char **description,
	      int NMAX)
{
  if(*nsite>=0)nres_site[*nsite]=*sres; *sres=0; (*nsite)++;
  if(*nsite > NMAX){
    printf("ERROR, too many sites > %d\n", NMAX); return(-1);
  }
  description[*nsite]=malloc(100*sizeof(char));
  return(0);
}


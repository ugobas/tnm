#include "coord.h"
#include "tnm.h"
#include "nma_para.h"
#include "nma.h"
#include "simulation.h"
#include "read.h"
#include "dof_tnm.h"
#include "kinetic_tnm.h"
#include "interactions_tnm.h"
#include "buildup.h"
#include "allocate.h"
#include "energy_BKV.h"

// #include "screened_interactions.h"
// #include "vector.h"
// #include "McLachlan.h"
// #include "align_tnm.h"
// #include "output_tnm.h"
// #include "contacts.h"
// #include "energy_BKV.h"
// #include "force2confchange.h"
// #include "allostery.h"
// #include "mutation.h"
// #include "diagonalize.h"
// #include "Fit_B.h"
// #include "rmsd.h"
// #include "rotation.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

float THR=4.5;
char INT_TYPE[40]="MIN";

/*****************************  Input  *********************************/
int getArgs(char *parameters, int argc, char **argv,
	    char *file_pdb, char *chain,
	    char *REF, int *SIDECHAINS, int *OMEGA, int *PSI,
	    int *MIN_INT, char *INT_TYPE, float *thr,
	    struct Para_simul *Para_simul);

int Read_para(char *filename,
	      char *file_pdb, char *chain,
	      char *REF, int *SIDECHAINS, int *OMEGA, int *PSI,
	      int *MIN_INT, char *INT_TYPE, float *thr,
	      struct Para_simul *Para_simul);
void help (char *name);

/***************************************************************************/
int main(int argc , char *argv[]){

  // Default parameters in .h file
  strcpy(REF, "EB"); C_THR=THR; N_RESRES=1; 

  // Parameters
  int SIDECHAINS=0; // Use side chains as degrees of freedom?
  int OMEGA=0;      // Use also omega as degree of freedom?
                    // 0=No 1=Yes -1=Yes, if D_Omega>90 degrees
  int PSI=1;        // Use psi angles?
  int MIN_INT=1;    // Minimum number of interactions per degree of freedom
  char  file_pdb[150], chain[200];
  struct Para_simul Para_simul;
  
  char parameters[1000];
  getArgs(parameters, argc, argv, file_pdb, chain, REF,
	  &SIDECHAINS, &OMEGA, &PSI, &MIN_INT, INT_TYPE, &C_THR, 
	  &Para_simul);

  /********************  Read protein   ***********************/
  int nres=0, nmr=0, natoms=0, ANISOU=0;
  atom *atoms=NULL;
  char  pdbid[100];
  struct residue *seq;

  // Ligand
  int n_lig=0, na_lig=0;
  struct residue *ligres=NULL;
  atom *ligatom=NULL;

  // Chains
  int Nchain=0;
  struct chain *chains=NULL;
  AA_code=AA_BKV;
  Read_PDB(&nres, &seq, &n_lig, &ligres, chain, &ANISOU, &nmr,
	   &natoms, &atoms, &na_lig, &ligatom, &chains, &Nchain,
	   pdbid, file_pdb);

  // Name
  char nameprot[100];
  sprintf(nameprot, "%s", pdbid);
  for(int i=0; i<Nchain; i++){
    if((chains[i].label!=' ')&&(chains[i].label!='\0'))
      sprintf(nameprot, "%s%c", nameprot, chains[i].label);
  }

  /*************************************/
  // Reference atoms
  struct Reference Ref_kin;
  int N_ref=Set_reference(&Ref_kin, 0, REF, atoms, 0, natoms);

  // Interactions
  struct interaction *Int_list; int N_int; 
  Compute_interactions(&N_int, &Int_list, INT_TYPE, atoms, natoms, nres, nameprot);
  printf("%d interactions of type %s, thr= %.2f\n", N_int, INT_TYPE, C_THR);

  // Set bonds
  struct bond *bonds=Set_bonds_topology(natoms, atoms, seq);

  // Degrees of freedom
  int last_ali_res=nres;
  int naxe=0, nmain=0, nrigid=0, N_diso=0, nskip=0;
  struct axe *axe=
    Set_DegofFreed(&naxe, &nmain, &nrigid, &nskip, &N_diso, bonds,
		   atoms, natoms, Ref_kin.atom_num, N_ref,
		   seq, last_ali_res, chains, Nchain,
		   Int_list, N_int, MIN_INT, MIN_INT,
		   OMEGA, SIDECHAINS, PSI);

  // Standardize bonds
  Standardize_bonds(bonds, atoms, natoms, nameprot,
		    axe, naxe, nmain, nrigid,
		    seq, nres, chains, Nchain, 0.1, Para_simul);
}

/***************************  Input **********************************/

int getArgs(char *parameters, int argc, char **argv,
	    char *file_pdb, char *chain,
	    char *REF, int *SIDECHAINS, int *OMEGA, int *PSI,
	    int *MIN_INT, char *INT_TYPE, float *C_THR,
	    struct Para_simul *Para_simul)
{
  int i; //p1=0, p2=0, out=0


  /************** Initialize ****************/
  // Input
  file_pdb[0]='\0'; strcpy(chain, "");

  // Simulations
  Para_simul->N_SIMUL=0;
  Para_simul->RESET=1;      // Reset native structure after build-up
  Para_simul->SELECT_ENE=0; // Build structures based on anharnomicity
  Para_simul->AMPLITUDE=10; // Max. amplitude for normal modes deformations
  Para_simul->D_REP=2.5;    // Threshold for repulsion
  Para_simul->RMSD_STEP=0.1;// If confchange projections smaller, discard
  Para_simul->MAX_ANGLE=0.7;// Max. angle for torsional deformations, radiants
  Para_simul->E_THR=5.0;    // Threshold for stopping motion

  if(argc<2)help(argv[0]);
  int infile=Read_para(argv[1], file_pdb, chain, REF, SIDECHAINS, OMEGA, PSI,
		       MIN_INT, INT_TYPE, C_THR, Para_simul);
  if(infile)goto inform;
  printf("WARNING, input file not specified or absent (%s)\n", argv[1]);
  printf("WARNING, reading parameters from command line\n");

  for(i=1; i<argc; i++){
    // Proteins:
    if (strncmp(argv[i],"-p",3)==0){
      i++; if(i>=argc)continue;
      strcpy(file_pdb ,argv[i]);
      printf("file_pdb=%s\n", file_pdb);
    }else if (strncmp(argv[i],"-ch",3)==0){
      i++; if(i>=argc)continue;
      if((strncmp(argv[i],"ALL", 3)==0)||(strncmp(argv[i],"all", 3)==0)){
	*chain='*';
      }else{
	strcpy(chain, argv[i]);
      }
      printf("chain=%s\n", chain);
    }else if (strncmp(argv[i],"-sc",3)==0){
      *SIDECHAINS=1;
    }else if (strncmp(argv[i],"-omega",6)==0){
      *OMEGA=1;
    }else if (strncmp(argv[i],"-ref",4)==0){
      i++; if(i>=argc)continue;
      strcpy(REF, argv[i]);
    }else if (strncmp(argv[i],"-cont_type",10)==0){
      i++; if(i>=argc)continue;
      strcpy(INT_TYPE, argv[i]);
    }else if (strncmp(argv[i],"-expo", 5)==0){
      i++; if(i>=argc)continue;
      sscanf(argv[i], "%f", &EXP_HESSIAN);
    }else if (strncmp(argv[i],"-cont_thr",9)==0){
      i++; if(i>=argc)continue;
      sscanf(argv[i], "%f", C_THR);
      printf("Threshold distance= %.1f\n", *C_THR);
    }else if (argv[i][0]=='-'){
      printf("WARNING, argument %s does not exist\n", argv[i]);
    }
  }

 inform:
  if(file_pdb[0]=='\0'){
    printf("ERROR, PDB file is mandatory\n"); help(argv[0]);
  }
  if(*OMEGA)printf("Omega angle taken as degree of freedom.\n");
  if(*SIDECHAINS){
    printf("Side chains degrees of freedom used\n");
    if(strncmp(INT_TYPE, "MIN", 3)!=0){
      printf("WARNING, interaction type %s not allowed with side chains\n",
	     INT_TYPE); strcpy(INT_TYPE, "MIN");
    }
    if((strncmp(REF, "ALL", 3)!=0)){
      printf("WARNING, reference atoms %s not allowed with side chains\n",
	     REF); strcpy(REF, "ALL");
    }
    printf("Reference atoms set to %s, interaction type set to %s\n",
	   REF, INT_TYPE);
    printf("Minimum number of interactions for accepting a degree of freedom");
    printf(": %d\n", *MIN_INT);
  }
  if((strncmp(INT_TYPE, "CA", 2)!=0)&&
     (strncmp(INT_TYPE, "CB", 2)!=0)&&
     (strncmp(INT_TYPE, "MIN", 3)!=0)){
    if(INT_TYPE[0]!='\0')
      printf("WARNING, interaction type %s does not exist\n", INT_TYPE);
    printf("Interaction type set to default %s\n", "MIN");
    strcpy(INT_TYPE, "MIN");
  }else{
    printf("Contact type set to %s\n", INT_TYPE);
  }
  
  if(REF[0]!='\0'){
    if((strncmp(REF, "CA", 2)!=0)&&
       (strncmp(REF, "CB", 2)==0)&&
       (strncmp(REF, "BB", 2)==0)&&
       (strncmp(REF, "EB", 2)==0)&&
       (strncmp(REF, "ALL", 3)==0)){
      printf("WARNING, %s are not allowed reference atoms\n", REF);
      REF[0]='\0';
    }
  }
  if(*C_THR==0){
    *C_THR=THR;
    printf("Contact threshold set to default %.1f\n", *C_THR);
  }
  sprintf(parameters, "Para: REF=%s Sidechain=%d Min_int=%d",
	  REF, *SIDECHAINS, *MIN_INT);
  sprintf(parameters, "%s Omega=%d ", parameters, *OMEGA);
  sprintf(parameters, "%s Psi=%d ", parameters, *PSI);
  sprintf(parameters,"%s CONT=%s ", parameters, INT_TYPE);
  sprintf(parameters, "%s thr=%.2f", parameters, *C_THR);
  fflush(stdout);
  return (0);
}

void help(char *PRG)
{
  fprintf(stderr, "Program %s\n", PRG);
  fprintf(stderr, "author Ugo Bastolla <ubastolla@cbm.csic.es> ");
  fprintf(stderr, "Centro de Biologia Molecular Severo Ochoa.\n");
  fprintf(stderr, "It standardizes bond lengths and bond angles of a protein structure and adjusts the torsion angles to minimize RMSD\n");
  fprintf(stderr, "   USAGE:");
  fprintf(stderr, "   %s <inputile>  or\n", PRG);
  fprintf(stderr, "   %s -p <pdbfile>", PRG);
  fprintf(stderr, "   OPTIONS:\n");
  fprintf(stderr, "       -h prints this help\n");
  fprintf(stderr, "       -ch <chain_id> <ALL>:reading all chains, A, AB..\n");
  fprintf(stderr, "       -ref Reference atoms. Allowed: CA CB BB EB ALL. ");
  fprintf(stderr, "       -omega  Use also omega angle as degree of freedom\n");
  fprintf(stderr, "Default: %s\n", REF);
  fprintf(stderr, "Interaction model:\n");
  fprintf(stderr, "       -cont_type Contact type CA CB, MIN ");
  fprintf(stderr, "Default: %s\n", INT_TYPE);
  fprintf(stderr, "       -cont_thr Distance threshold default %.1f\n",THR);
  fprintf(stderr, "\n");
  fprintf(stderr, "FORMAT of input file:\n");
  fprintf(stderr, "####   INPUT:\n");
  fprintf(stderr, "PDB= 1usg.pdb ");
  fprintf(stderr, "! Input structure\n");
  fprintf(stderr, "CH=  A                                   ");
  fprintf(stderr, "! Chain\n");
  fprintf(stderr, "#\n");
  fprintf(stderr, "####  Model parameters (reccomended: do not change)\n");
  fprintf(stderr,
	  "#===============================================================\n");
  fprintf(stderr, "REF=  EB CB  CA BB ALL                   ");
  fprintf(stderr, "! Reference atoms\n");
  fprintf(stderr, "SIDECHAIN=0                              ");
  fprintf(stderr, "! Use sidechains as degree of freedom?\n");
  fprintf(stderr,
	  "# This option allows using side chains chi angles as additional degree of\n");
  fprintf(stderr,
	  "# freedom, but introduces non collective motions that worsen performances.\n");
  fprintf(stderr, "MIN_INT=1                              ");
  fprintf(stderr, "! Min. number of interactions to accept a dof\n");
  fprintf(stderr, "OMEGA=0                                  ");
  fprintf(stderr, "! Use omega angle as degree of freedom?\n");
  fprintf(stderr, "PSI=1                                    ");
  fprintf(stderr, "! Use psi angle as degree of freedom?\n");

  exit(1);
}

void GetPdbId(char *pdb_file_in, char *pdbid){
     /* This subroutine pretends to get the
        PDB id from a pdb file name, and ressembles
	quite a lot my "old-and-dirty-Perl" days */

  int start=0, end=0, i,j; //end2=0

  for(i=strlen(pdb_file_in)-1;i>=0;i--){
    if (pdb_file_in[i]=='.'){
      end=i-1;
    }else if (pdb_file_in[i]=='/'){
      start=i+1; //end2=i-1;
      break;
    }
  }
  j=0;
  for (i=start;i<=end;i++){
    pdbid[j]=pdb_file_in[i];
    j++;
  }
  pdbid[j]='\0';
}

int Read_para(char *filename,
	      char *file_pdb, char *chain,
	      char *REF, int *SIDECHAINS, int *OMEGA, int *PSI,
	      int *MIN_INT, char *INT_TYPE, float *C_THR,
	      struct Para_simul *Para_simul)
{
  FILE *file_in=fopen(filename, "r");
  if(file_in==NULL){
    printf("ERROR, TNM input file %s not found\n", filename); 
    return(0);
  }
  char string[1000], dumm[80];
  printf("Reading parameters in %s\n", filename);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "PDB1 ", 4)==0){
      sscanf(string+5, "%s", file_pdb);
    }else if(strncmp(string, "CH",2 )==0){
      sscanf(string+4, "%s", dumm);
      if((strncmp(dumm,"ALL", 3)==0)||(strncmp(dumm,"all", 3)==0)){
	*chain='*';
      }else{
	strcpy(chain,dumm);
      }
    }else if(strncmp(string, "REF", 3)==0){
      sscanf(string+4, "%s", REF);
    }else if(strncmp(string, "SIDECHAIN", 9)==0){
      sscanf(string+10, "%d", SIDECHAINS);
    }else if(strncmp(string, "MIN_INT", 7)==0){
      sscanf(string+8, "%d", MIN_INT);
   }else if(strncmp(string, "PSI", 3)==0){
      sscanf(string+4, "%d", PSI);
    }else if(strncmp(string, "OMEGA", 5)==0){
      sscanf(string+6, "%d", OMEGA);
      if(*OMEGA>0){printf("Omega angle used as degree of freedom\n");}
      else if(*OMEGA<0){printf("Omega angle used as dof if D_omega>90\n");}
    }else if(strncmp(string, "CONT_TYPE", 9)==0){
      sscanf(string+10, "%s", INT_TYPE);
    }else if(strncmp(string, "CONT_THR", 8)==0){
      sscanf(string+ 9, "%f", C_THR);
    }else if(string[0]!='\n'){
      printf("WARNING, unrecognized line:\n%s", string);
    }
  }
  fclose(file_in);
  return(1);
}


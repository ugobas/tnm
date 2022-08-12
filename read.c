#include "coord.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "read.h"
#include "nma_para.h"
#include "externals.h"

int SC_BEFORE=1; // Put side chain before CO
char CHAIN_CODE[200];

#define AMIN_CODE "AEQDNLGKSVRTPIMFYCWHX"
#define ASTPATH    "/data/ortizg/databases/astral_40/"
#define PDBPATH    "/data/ortizg/databases/pdb/"
#define PDBEXT     ".pdb" /*.ent.Z*/
#define PDBCAT     "zcat"  /* zcat */
#define PDBTMP     "pdb.dat"


// Arrays limiters
// #define ATOM_MAX 50000
#define ATOM_MAX 50 // Number of atoms per residue
int MIN_ATOM=2; // Minimum number of atoms to set a chain
#define DEBUG 0

static int n_atom;
static char res_exo[400][5], res_std[400][5];
char modres[20][100], stdres[20][5]; //modyves: these were allocated inside a function, but never freed.
int n_modres[20], ini_modres=0;
char mod_0[]="CSD MDO ORN DBZ NAL LAL HAC AYA NCB "; // Ala
char mod_1[]="CYG CY3 CSU CSP CYM CAS CSB CSR CME CMH CSS CSX CSW CSO ALS SMC CEA OCS CCS SCY YCM"; // Cys 
char mod_2[]="BFD SNN ASA BHD SUI "; // Asp
char mod_3[]="PCA AR4 GLQ GMA ";    // Glu 
char mod_4[]="DPN"; // Phe
char mod_5[]="GL3 CHG GLZ ACY FGL "; //  Gly 
char mod_6[]="MHS NEP HIC "; //  His
char mod_7[]="IIL "; // Ile
char mod_8[]="LLP KCX LYZ ALY LCX MCL LYN "; // Lys 
char mod_9[]="LEF DLE "; //  Leu 
char mod_10[]="MSE FME MME MSO FOR NRQ CH6 "; // Met
char mod_11[]="IAS MEN "; // Asn
char mod_12[]="HYP DPR "; // Pro 
char mod_13[]="MGN 5HP "; // Gln 
char mod_14[]="AGM BOR ARM ACL DAR "; // Arg
char mod_15[]="OSE DSN HSL DHL SAC SEP "; // Ser
char mod_16[]="TPO "; //  Thr
char mod_17[]="DVA "; // Val
char mod_18[]="TRO TRQ TRW TRN "; //  Trp
char mod_19[]="TPQ PTR STY YOF TYI TYC TYS "; // Tyr

void GetPdbId(char *pdb_file_in, char *pdbid);

// Residues
void Ini_modres(); //modyves: parameters are global static variables
static int Next_residue(int *N_res, int *start, struct residue *seq,
			atom *first_atom, short *i_atom,
			char *res_type_old, int *res_num_old, char *icode_old,
			char *chain_old, int *hetatm_old, int n_exo,
			char *res_type, int res_num, char icode, char chain,
			int hetatm);
static short Write_residue(char *res_type_old, atom *first_atom, short i_atom,
			   struct residue *ptr_tmp, int n_exo, int res_num,
			   char icode, char chain, int hetatm);
void Code_3_2(char *res2, char *res3);
extern int Code_AA(char res);

// Chains
int Set_chains(struct chain **chains,
	       atom **atoms, int *N_atoms,
	       struct residue *res, int *nres);
void Store_chain(struct chain *chp, int *ichain, int n, int *nr);
int Check_chain(int *ichain, int *sel, int *na, int ini_r, int nres);

// Atoms
extern int Find_atom(atom *atoms, int *i1, int res, int Natoms, char *type);
int Order_atoms(atom *atoms, int natoms, struct residue *seq, int Nres);
int Order_atoms_aa(atom *res_atom, int n_atom);
int Store_branch(atom **atom2, char BRANCH,
		 atom *res_atom, int n_atom);
int Copy_atom(atom **atom2, atom *atom1, int i);
#define NA_RES_MAX 500
atom atom_tmp[NA_RES_MAX];

// Start codes
//////////////////////////////////////////////////////////////////////////////////////
void Read_PDB(int *nres,   // Number of amino acid residues
	      struct residue **seq, // Where ligand and aa res are stored
	      int *n_lig,  // Number of ligand residues
	      struct residue **ligres, // Pointer to first ligand residue
	      char *chain, // Chain identifiers
	      int *ANISOU, // Anisotropic B factors?
	      int *nmr,    // NMR structure?
	      int *natoms, // Number of protein atoms
	      atom **atoms,// Pointer to protein atoms
	      int *na_lig, // Number of ligand atoms
	      atom **ligatoms, // Pointer to ligand atoms
	      struct chain **chains, // Pointer to protein chains
	      int *Nchain, // Number of protein chains (1-Nchain)
	      float *TEMP, // Temperature of the experiment, if given
	      char *pdbid, // PDB identifier
	      char *file_pdb) // INPUT: Path to PDB file
{
  strcpy(CHAIN_CODE,"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz:=-_ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz:=-_");

  /****************** Reading atoms: ************************/
  int i; 
  int numres=Count_residues(file_pdb, chain); L_MAX=numres;
  int numatoms=ATOM_MAX*numres;
  atom *atom_read=malloc(numatoms*sizeof(atom));
  *seq=malloc(numres*sizeof(struct residue)); 
  //for(i=0;i<numres;i++){(*seq)[i].chain='-';}
  
  *nres=Read_coord(file_pdb, nmr, *seq, atom_read, chain, ANISOU, TEMP);
  if(*nres<=0){
    printf("WARNING, file %s, no residue found\n", file_pdb); return;
  }
  GetPdbId(file_pdb,pdbid); //if(*chain==' ')*chain='_';

  // *n_lig=Read_ligand(file_pdb, *seq+*nres, atom_read+*natoms, na_lig);
  // Ligand atoms are stored at the end

  
  /****************** Storing atoms: ************************/
  *natoms=0;  
  for(i=0; i<*nres; i++)*natoms+=(*seq)[i].n_atom;
  *atoms=malloc((*natoms+*na_lig)*sizeof(atom));
  printf("PDB %s chain %s nres=%d natoms=%d\n",pdbid, chain, *nres, *natoms);

  /**************** Set chains *********************/
  // Eliminate chains with fewer than MIN_ATOM atoms and their residues
  // Allocate and store chains
  // Allocate and store atoms
  *Nchain=Set_chains(chains, atoms, natoms, *seq, nres);
  printf("%d chains found: ", *Nchain);
  for(i=0; i<*Nchain; i++)printf("%c", (*chains)[i].label);
  printf(".\n");
  free(atom_read); //modyves: this was never freed

  int N_seqres[*Nchain];
  char *seqres[*Nchain], chainlabel[*Nchain+1];
  for(i=0; i<*Nchain; i++)chainlabel[i]=(*chains)[i].label;
  chainlabel[*Nchain]='\0';
  if(Read_seqres(seqres, N_seqres, *Nchain, file_pdb, chainlabel)<0){
    // Seqres not found, copy from atoms
    for(i=0; i<*Nchain; i++){
      struct chain *ch=*chains+i;
      int L=ch->nres;
      N_seqres[i]=L;
      seqres[i]=malloc(L*sizeof(char)); char *seq=seqres[i];
      for(int j=0; j<L; j++)seq[j]=ch->res[j].amm;
    }
  }
  for(i=0; i<*Nchain; i++){
    struct chain *ch=*chains+i;
    ch->seqres=seqres[i];
    ch->N_seqres=N_seqres[i];
    ch->seq=malloc(ch->nres*sizeof(char));
    for(int j=0; j<ch->nres; j++)ch->seq[j]=ch->res[j].amm;
    ch->ali_seqres=malloc(N_seqres[i]*sizeof(int));
    Align_seqres(ch->ali_seqres, ch->seqres, ch->N_seqres, ch->seq, ch->nres);
  }
  
  if(SC_BEFORE){
    // Change the order of protein atoms setting side chain before CO
    // and testing that CA is before C
    Order_atoms(*atoms, *natoms, *seq, *nres);
    // Find next atom in original PDB order, i_next
    Find_PDB_order(*atoms, *natoms);
  }

  return;
}

/////////////////////////////////////////////////////////////////////////
int Count_residues(char *pdb_name, char *chain_to_read){
  int N_res=0, Compression;
  FILE *file_in=Open_compressed_file(pdb_name, &Compression);
  if(file_in==NULL)return(0);
  if(ini_modres==0)Ini_modres(); //modyves: parameters were global static variables

  char string[200];
  int hetatm=0;
  int i, res_num, res_num_old=10000;
  char  chain=' ', chain_old='\0';
  char res_type[5], res_type_old[5], icode=' ', icode_old='.';
  float x, y, z;
  int nchains=0, ic, ini_chain=0, readall=0, *rchain=NULL;
  char *xchain=NULL;//modyves: there were two things called "chain"

  char *chain_pdb=CHAIN_CODE;

  // Count chains to read
  if((chain_to_read[0]=='*')||
     ((chain_to_read[0]=='0')&&(chain_to_read[1]=='1'))){
    readall=1;
    printf("Counting residues in ALL chains\n");
  }else{
    i=0; while(chain_to_read[i]!='\0')i++;
    nchains=i;  if(nchains==0)nchains=1;
    rchain=(int *)malloc(nchains*sizeof(int));
    for(i=0; i<nchains; i++)rchain[i]=0;
    printf("Counting residues in %d chains %s\n",
	   nchains, chain_to_read);
  }

  strcpy(res_type_old,"xxx");
  while(fgets(string, sizeof(string), file_in)!=NULL){

    if(strncmp(string,"ATOM", 4)==0){
      /* Standard residue or DNA basis */
      hetatm=0; if(ini_chain==0)ini_chain=1;
    }else if(strncmp(string,"HETATM", 6)==0){
      if(ini_chain==0)continue;
      hetatm=1;
      /* Cofactor or exotic residue */
    }else if((strncmp(string,"TER",3)==0)&&(N_res>0)){
      chain_pdb++; ini_chain=0; continue;
    }else if(strncmp(string,"ENDMDL", 6)==0){
      break;                                    /* end model */
    }else{
      continue;
    }

    chain=string[21];
    if(chain==' ')chain=*chain_pdb;
    if(readall){
      if(chain!=chain_old){
	//chain_to_read[nchains]=chain;
	nchains++; chain_old=chain;
      }
    }else{
      if((chain_to_read[0]==' ')||(chain_to_read[0]=='\0')){
	chain_to_read[0]=chain; chain_to_read[1]='\0';
      }
      for(ic=0; ic<nchains; ic++){
	if(chain==chain_to_read[ic]){rchain[ic]=1; goto read;}
      }
      continue;
    }

  read:    /* Omit hydrogen if HYD=0 */
    if(HYD==0){
      if((string[13]=='H')||(string[12]=='H') ||
	 (string[13]=='D') ||(string[12]=='D'))continue;
    }

    /* Read residue; check if water molecule */
    res_type[0]=string[17]; res_type[1]=string[18]; res_type[2]=string[19];
    res_type[3]='\0';
    if((hetatm==1)&&((strncmp(res_type,"HOH",3)==0)||
		     (strncmp(res_type,"DOD",3)==0)))continue;
    icode=string[26]; string[26]=' ';
    sscanf(string+22,"%d %f %f %f", &res_num, &x, &y, &z);

    /* New residue */
    if((res_num!=res_num_old)||(icode!=icode_old)||
       (strncmp(res_type, res_type_old, 3)!=0)){
      N_res++;
      strcpy(res_type_old, res_type); icode_old=icode; res_num_old=res_num;
    }
  }
  //N_res++;
  fclose(file_in);
  if(Verbose)printf("%3d residues\n", N_res);
  if(Compression)Delete_tmp_file();
  if(readall==0){
    // Count how many chains have been read
    int nc=0;
    xchain=(char *)malloc(nchains*sizeof(char)); //modyves
    for(i=0; i<nchains; i++){
      if(rchain[i]){xchain[nc]=chain_to_read[i]; nc++;}//modyves
      else{printf("WARNING, chain %c not found\n", chain_to_read[i]);}
    }
    for(i=0; i<nc; i++)chain_to_read[i]=xchain[i];
    nchains=nc;//modyves
  }
  chain_to_read[nchains]='\0';
  printf("Counting %3d residues in %d chains: %s\n",N_res, nchains, chain_to_read);
  if(readall==0)//modyves
  	{free(xchain);xchain=NULL; 
   	 free(rchain);rchain=NULL;
  	}
  return(N_res);
}

/////////////////////////////////////////////////////////////////////////
int Read_coord(char *pdb_name, int *nmr,
	       struct residue *seq, atom *atoms,
	       char *chain_to_read, int *ANISOU, float *TEMP)
{
  /* Open file */
  int Compression;
  FILE *file_in=Open_compressed_file(pdb_name, &Compression);
  if(file_in==NULL)return(0);

  int N_res=0, n_exo=0, start=0;
  char string[200];
  atom *atom_ptr=atoms, *first_atom=atoms;

  short i_atom=0, alternative=0;
  int hetatm=0, hetatm_old=0;
  int i, j, res_num, res_num_old=10000;
  char altloc, altloc_sel=' ';
  char chain_old='\0', chain=' ';
  char res_type[5], res_type_old[5], icode=' ', icode_old='U';
  float x, y, z;
  int nchains=0, ic, ichain=0, ini_chain=0, *rchain=NULL, readall=0;
  char *chain_pdb=CHAIN_CODE;

  n_atom=0;
  if(HYD)printf("Considering hydrogen atoms\n");

  // Count chains to read
  if(chain_to_read[0]=='*'){
    readall=1;
     printf("Reading atoms in all chains\n");
  }else{
    i=0;
    while(chain_to_read[i]!='\0')i++;
    nchains=i; 
    if(nchains==0)nchains=1;
    rchain=(int *)malloc(nchains*sizeof(int));
    for(i=0; i<nchains; i++)rchain[i]=0;
    printf("%d chains to read\n", nchains);
  }
  
  *nmr=0;
  strcpy(res_type_old,"xxx");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string,"ATOM", 4)==0){
      /* Standard residue or DNA basis */
      hetatm=0;
      if(ini_chain==0)ini_chain=1; 	 
    }else if(strncmp(string,"HETATM", 6)==0){
      if(ini_chain==0)continue;
      hetatm=1;
      /* Cofactor or exotic residue */
      }else if(TEMP && strncmp(string,"REMARK 200", 10)==0){
      char word1[40], word2[40], word3[40], word4[40];
      sscanf(string+12, "%s%s%s%s", word1, word2, word3, word4);
      if((strncmp(word1,"TEMPERATURE", 11)==0)&&
	 (strncmp(word4,"NULL", 4)!=0)){
	sscanf(word4, "%f", TEMP);
	if(strncmp(word2,"(KELVIN)", 8)!=0){*TEMP+=273;}
	}
      continue;
    }else if(strncmp(string,"EXPDTA", 6)==0){
      char word1[80], word2[80];
      sscanf(string+8, "%s%s", word1, word2);
      if(strncmp(word1, "NMR", 3)==0){*nmr=1;}
      else if(strncmp(word2, "NMR", 3)==0){*nmr=1;}
      continue;
      /* NMR structure */
    }else if((strncmp(string,"TER",3)==0)&&(N_res>0)){
      chain_pdb++;
      //printf("-NR- |%3d|%c|%c|\n",N_res,chain_old,chain);
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
      ini_chain=0;
      continue;
    }else if(strncmp(string,"MODRES", 6)==0){
      res_exo[n_exo][0]=string[12]; res_std[n_exo][0]=string[24];
      res_exo[n_exo][1]=string[13]; res_std[n_exo][1]=string[25];
      res_exo[n_exo][2]=string[14]; res_std[n_exo][2]=string[26];
      for(j=0; j<n_exo; j++){
	if(strncmp(res_exo[j],res_exo[n_exo],3)==0)break;
      }
      if(j==n_exo)n_exo++;
      continue;
    }else if(strncmp(string,"ENDMDL", 6)==0){
      break;                                    /* end model */
    }else if(strncmp(string,"ANISOU", 6)==0){
      // Anisotropic structure factor
      int i, j;
      if(*ANISOU==0)*ANISOU=1;
      //for(i=0; i<3; i++)atom_ptr->anisou[i]=malloc(3*sizeof(float));
      sscanf(string+28, "%f %f %f %f %f %f",
	     &(atom_ptr->anisou[0][0]), &(atom_ptr->anisou[1][1]),
	     &(atom_ptr->anisou[2][2]), &(atom_ptr->anisou[0][1]),
	     &(atom_ptr->anisou[0][2]), &(atom_ptr->anisou[1][2]));
      for(i=0; i<3; i++)
	{for(j=i; j<3; j++)
	    {atom_ptr->anisou[i][j]*=0.0001;
	      if(j!=i)atom_ptr->anisou[j][i]=atom_ptr->anisou[i][j];
	    }
	}
      // Check
      //for(i=0; i<3; i++)B+=aniso[i][i]; B*=26.319; // 8pi^2/2
      //printf("B= %.2f %.2f %d\n", atom_ptr->B_factor, B, n_atom);
      continue;
      /* }else if(strncmp(string,"HELIX ", 6)==0){
	 Read_sec_str(string, chain, 'H'); continue;
	 }else if(strncmp(string,"SHEET ", 6)==0){
	 Read_sec_str(string, chain, 'E'); continue;
	 }else if(strncmp(string,"TURN ", 5)==0){
	 Read_sec_str(string, chain, 'T'); continue;
	 }*/
    }else{
      continue;
    }

    chain=string[21];
    if(chain==' ')chain=*chain_pdb;
    if(readall)
    	{if(chain!=chain_old)
    		{printf("Reading chain %c\n", chain);
		 	 chain_to_read[nchains]=chain; nchains++; chain_old=chain;
      		}
    	}
    else
    	{if((chain_to_read[0]==' ')||(chain_to_read[0]=='\0'))
    		{chain_to_read[0]=chain; chain_to_read[1]='\0';
      		}
      	 for(ic=0; ic<nchains; ic++)
      	 	{if(chain==chain_to_read[ic])
      	 		{rchain[ic]=1; goto read;}
      		}
      	 continue;
    	}

    /*if(*chain_to_read!='*'){
      if((*chain_to_read==' ')||(*chain_to_read=='\0'))*chain_to_read=chain;
      for(ic=0; ic<nchains; ic++)if(chain==chain_to_read[ic])goto read;
      continue;
      }*/




  read:
    // Hydrogen atoms
    if((string[13]=='H')||(string[12]=='H') || (string[13]=='D') ||(string[12]=='D'))
    	{if(HYD==0){continue;}   // Omit hydrogen
      	 else if((string[13]=='H')&&(string[12]!=' '))
      	 	{// if name 1HB1 change to HB11
		  	 char c=string[12];
			 if((c=='1')||(c=='2')||(c=='3'))
			 	{string[12]='H'; string[13]=string[14];
	  			 string[14]=string[15]; string[15]=c; 
				}
      		}
    	}
  
    /* Read residue; check if water molecule */
    res_type[0]=string[17]; res_type[1]=string[18]; res_type[2]=string[19];
    res_type[3]='\0';
    if((hetatm==1)&&((strncmp(res_type,"HOH",3)==0)||(strncmp(res_type,"DOD",3)==0)))continue;

    icode=string[26]; string[26]=' ';

    /* Read coordinates */
    sscanf(&string[22],"%d %f %f %f", &res_num, &x, &y, &z); //modyves

    /* Check if alternative conformation */
    if((string[72]=='A')&&(string[73]=='L')&&(string[74]=='T')&&(string[75]!='1')&&(string[75]!=' '))continue;

    altloc=string[16];
    if(altloc!=' '){
      if(altloc_sel==' ')altloc_sel=altloc;
      if(altloc!=altloc_sel)continue;
    }

    if((icode!=icode_old)&&(res_num==res_num_old)){
      if(alternative==1){
	continue;
      }else{
	atom *atom_old=first_atom;
	float dx, dy, dz;
	dx=x-atom_old->r[0]; dy=y-atom_old->r[1]; dz=z-atom_old->r[2];
	if((dx*dx+dy*dy+dz*dz)<.5){alternative=1; continue;}
      }
    }

    /* New residue */
    if((res_num!=res_num_old)||(icode!=icode_old)||
       (strncmp(res_type, res_type_old, 3)!=0)){
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
    }
    atom_ptr=atoms+n_atom;

    if(i_atom==0)first_atom=atom_ptr;
    i_atom++;
    n_atom++;
    atom_ptr->r[0] = x;
    atom_ptr->r[1] = y;
    atom_ptr->r[2] = z;
    
    // Copy name
    if(string[12]!=' ')
    	{ // Hydrogen atoms
      	 if(strncmp(string+12, "SE", 2)==0)
      	 	{// Selenium in Seleniomethionine is at delta position
			 string[13]='D';
      		}
      	 for(ic=0; ic<3; ic++)atom_ptr->name[ic]=string[12+ic];
    	}
    else if(string[13]!=' ')
     	{for(ic=0; ic<3; ic++)atom_ptr->name[ic]=string[13+ic];
   		}
    else
    	{for(ic=0; ic<3; ic++)atom_ptr->name[ic]=string[14+ic];
    	}
    atom_ptr->name[3]='\0';
    
    if(string[62]!=' ')
    	{sscanf(string+60, "%f", &atom_ptr->B_factor);
     	}
     else
     	{atom_ptr->B_factor=0;
    	}
    if(string[60]!=' ')string[60]=' ';
    if(string[56]!=' ')
    	{sscanf(string+54, "%f", &atom_ptr->occupancy);
    	}
    else
    	{atom_ptr->occupancy=0;
    	}

    atom_ptr->chain=ichain;
  }
  Next_residue(&N_res, &start, seq, first_atom, &i_atom,
	       res_type_old, &res_num_old, &icode_old,
	       &chain_old, &hetatm_old, n_exo,
	       res_type, res_num, icode, chain, hetatm);
	       
  if(N_res > L_MAX)
  	{printf("\n ERROR, more than %d residues found\n", L_MAX); exit(8);
  	}
  fclose(file_in);
  printf("Reading %3d residues in %d chains: %s\n",N_res, nchains, chain_to_read);
  chain_to_read[nchains]='\0';
  
  
  if(readall==0){free(rchain);rchain=NULL;} //modyves: this was never freed
  if(Compression)Delete_tmp_file();
  return(N_res);
}

//////////////////////////////////////////////////////////////////////////////////
int Read_ligand(char *pdb_name, int *nmr, struct residue *seq, atom *atoms,
	       char *chain_to_read, int *ANISOU)
{
  int N_res=0, n_exo=0, Compression=0, start=0;
  FILE *file_in;
  char string[200], command[200];
  atom *atom_ptr=atoms, *first_atom=NULL;

  short i_atom=0, alternative=0;
  int hetatm=0, hetatm_old=0;
  int i, j, res_num, res_num_old=10000;
  char altloc, altloc_sel=' ';
  char chain_old='#', chain=' ';
  char res_type[5], res_type_old[5], icode=' ', icode_old='U';
  float x, y, z;
  char file_name[500];
  int nchains=1, ic, ichain=0, ini_chain=0;
  n_atom=0;

  /* Open file */
  Compression=Get_compression(pdb_name);
  if(Compression){
    sprintf(command, "%s %s > %s\n", PDBCAT, pdb_name, PDBTMP);
    system(command); strcpy(file_name, PDBTMP);
  }else{
    sprintf(file_name, "%s", pdb_name);
  }
  file_in=fopen(file_name, "r");
  if(file_in==NULL){
    printf("\nERROR, %s not found\n", file_name); return(0);
  }

  // Count chains to read
  for(i=0; i<sizeof(chain_to_read); i++){
    nchains=i; if(chain_to_read[i]=='\0')break;
  }
  nchains++; if(nchains==0)nchains=1;

  *nmr=0;
  strcpy(res_type_old,"xxx");
  while(fgets(string, sizeof(string), file_in)!=NULL){

    if(strncmp(string,"ATOM", 4)==0){
      /* Standard residue or DNA basis */
      hetatm=0; if(ini_chain==0)ini_chain=1;
    }else if(strncmp(string,"HETATM", 6)==0){
      if(ini_chain==0)continue;
      hetatm=1;
      /* Cofactor or exotic residue */

    }else if(strncmp(string,"EXPDTA", 6)==0){
      if(strncmp(string+10, "NMR", 3)==0)*nmr=1;
      continue;
      /* NMR structure */

    }else if((strncmp(string,"TER",3)==0)&&(N_res>0)){
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
      ini_chain=0;
      continue;

    }else if(strncmp(string,"MODRES", 6)==0){
      res_exo[n_exo][0]=string[12]; res_std[n_exo][0]=string[24];
      res_exo[n_exo][1]=string[13]; res_std[n_exo][1]=string[25];
      res_exo[n_exo][2]=string[14]; res_std[n_exo][2]=string[26];
      for(j=0; j<n_exo; j++){
	if(strncmp(res_exo[j],res_exo[n_exo],3)==0)break;
      }
      if(j==n_exo)n_exo++;
      continue;

    }else if(strncmp(string,"ENDMDL", 6)==0){
      break;                                    /* end model */
    }else if(strncmp(string,"ANISOU", 6)==0){
      // Anisotropic structure factor
      int i, j;
      if(*ANISOU==0)*ANISOU=1;
      sscanf(string+28, "%f %f %f %f %f %f",
	     &(atom_ptr->anisou[0][0]), &(atom_ptr->anisou[1][1]),
	     &(atom_ptr->anisou[2][2]), &(atom_ptr->anisou[0][1]),
	     &(atom_ptr->anisou[0][2]), &(atom_ptr->anisou[1][2]));
      for(i=0; i<3; i++){
	for(j=i; j<3; j++){
	  atom_ptr->anisou[i][j]*=0.0001;
	  if(j!=i)atom_ptr->anisou[j][i]=atom_ptr->anisou[i][j];
	}
      }
      // Check
      //for(i=0; i<3; i++)B+=aniso[i][i]; B*=26.319; // 8pi^2/2
      //printf("B= %.2f %.2f %d\n", atom_ptr->B_factor, B, n_atom);
      continue;

      /*}else if(strncmp(string,"HELIX ", 6)==0){
	Read_sec_str(string, chain, 'H'); continue;
	}else if(strncmp(string,"SHEET ", 6)==0){
	Read_sec_str(string, chain, 'E'); continue;
	}else if(strncmp(string,"TURN ", 5)==0){
	Read_sec_str(string, chain, 'T'); continue;*/
    }else{
      continue;
    }

    chain=string[21];
    if(*chain_to_read!='*'){
      if((*chain_to_read==' ')||(*chain_to_read=='\0'))*chain_to_read=chain;
      for(ic=0; ic<nchains; ic++)if(chain==chain_to_read[ic])goto read;
      continue;
    }

  read:
    /* Read atom name */
    if((string[13]=='H')||(string[12]=='H') ||
       (string[13]=='D') ||(string[12]=='D'))continue;

    /* Read residue; check if water molecule */
    res_type[0]=string[17]; res_type[1]=string[18]; res_type[2]=string[19];
    res_type[3]='\0';
    if((hetatm==1)&&((strncmp(res_type,"HOH",3)==0)||
		     (strncmp(res_type,"DOD",3)==0)))continue;

    icode=string[26]; string[26]=' ';


    /* Read coordinates */
    sscanf(string+22,"%d %f %f %f", &res_num, &x, &y, &z);

    /* Check if alternative conformation */
    if((string[72]=='A')&&(string[73]=='L')&&(string[74]=='T')&&
       (string[75]!='1')&&(string[75]!=' '))continue;

    altloc=string[16];
    if(altloc!=' '){
      if(altloc_sel==' ')altloc_sel=altloc;
      if(altloc!=altloc_sel)continue;
    }

    if((icode!=icode_old)&&(res_num==res_num_old)){
      if(alternative==1){
	continue;
      }else{
	atom *atom_old=first_atom;
	float dx, dy, dz;
	dx=x-atom_old->r[0]; dy=y-atom_old->r[1]; dz=z-atom_old->r[2];
	if((dx*dx+dy*dy+dz*dz)<.5){
	  alternative=1; continue;
	}
      }
    }

    /* New residue */
    if((res_num!=res_num_old)||(icode!=icode_old)||
       (strncmp(res_type, res_type_old, 3)!=0)){
      Next_residue(&N_res, &start, seq, first_atom, &i_atom,
		   res_type_old, &res_num_old, &icode_old,
		   &chain_old, &hetatm_old, n_exo,
		   res_type, res_num, icode, chain, hetatm);
    }
    atom_ptr=atoms+n_atom;

    if(i_atom==0)first_atom=atom_ptr;
    i_atom++; n_atom++;
    atom_ptr->r[0] = x;
    atom_ptr->r[1] = y;
    atom_ptr->r[2] = z;
    if(string[12]!=' '){ // Hydrogen atoms
      for(ic=0; ic<4; ic++)atom_ptr->name[ic]=string[12+ic];
    }else{
      for(ic=0; ic<3; ic++)atom_ptr->name[ic]=string[13+ic];
      atom_ptr->name[3]=' ';
    }
    sscanf(string+56, "%f", &atom_ptr->occupancy);
    sscanf(string+60, "%f", &atom_ptr->B_factor);
    atom_ptr->chain=ichain;
  }
  Next_residue(&N_res, &start, seq, first_atom, &i_atom,
	       res_type_old, &res_num_old, &icode_old,
	       &chain_old, &hetatm_old, n_exo,
	       res_type, res_num, icode, chain, hetatm);
  if(N_res > L_MAX){
    printf("\n ERROR, more than %d residues found\n", L_MAX); exit(8);
  }
  fclose(file_in);
  if(Verbose)printf("%3d residues\n", N_res);
  if(Compression){
    sprintf(command, "rm -f %s\n", PDBTMP); system(command);
  }
  return(N_res);
}

////////////////////////////////////////////////////////////////////////////////
static short Write_residue(char *res_type_old, atom *first_atom, short i_atom,
			   struct residue *ptr_tmp, int n_exo, int res_num,
			   char icode, char xchain, int hetatm)
{
  /* Returns type =1 (a.a.), -1 (a.a. as ligand), 2/-2: DNA 3/-3: RNA */
  char amm[5]; //modyves: this was the cause of the overwriting problems I had (that were not visible for you, though something else might have gotten overwritten)
  int exo=0, type;
  type = Get_amm(amm, &exo, res_type_old, hetatm, n_exo);

  	
  /* Check backbone */
  if((type==0) && (i_atom >=3))
  	{// If backbone atoms exist: Modified residue
     int i_N=0, i_CA=0, i_C=0;
     atom *atom_ptr=first_atom;
    
     for(int i=0; i<i_atom; i++)
     	{if(strncmp(atom_ptr->name, "N ", 2)==0){i_N=1;}
      	 else if(strncmp(atom_ptr->name, "CA", 2)==0){i_CA=1;}
      	 else if(strncmp(atom_ptr->name, "C ", 2)==0){i_C=1;}
      	 if(i_N && i_CA && i_C){type=1; exo=1; break;}
      	 atom_ptr++;
    	}
  	}

  ptr_tmp->type=type;
  ptr_tmp->chain=xchain;

  strcpy(ptr_tmp->pdbres,"     "); //modyves: this was sometimes uninitialized
  sprintf(ptr_tmp->pdbres, "%d", res_num);
  char tmp[40];
  if(icode!=' '){sprintf(tmp, "%c", icode); strcat(ptr_tmp->pdbres, tmp);}
  ptr_tmp->amm=amm[1]; ptr_tmp->exo=exo;
  if((type==1)||(type==-1)){
    ptr_tmp->i_aa=Code_AA(amm[1]);
    if(ptr_tmp->i_aa<0){
      printf("Unknown residue %s %d%c\n", res_type_old, res_num, icode);
      ptr_tmp->i_aa=0;
    }
  }
  ptr_tmp->atom_ptr=first_atom; 
  ptr_tmp->n_atom=i_atom;

  if(type==0)
    printf("Group %s %d%c  %s (%d atoms) not a residue\n",
	   res_type_old, res_num, icode, amm, i_atom);

  return(type);
}

////////////////////////////////////////////////////////////////////////////////
int Get_amm(char *amm, int *exo, char *res_type_old, int hetatm, int n_exo)
{
  /* Check amino acid type */
  Code_3_2(amm, res_type_old);
  int type=0, exo_check=0;

 find_type:
  if(amm[0] == ' '){
    // Standard amino acid residue
    if(hetatm==0 || exo_check){type=1;}
    else{type=-1;} // Standard res. and HETATM => ligand
  }else if(amm[0]=='D'){
    // deoxyribo-nucleotide
    if(hetatm==0 || exo_check){type=2;}
    else{type=-2;} // Standard nuc and HETATM => ligand
  }else if(amm[0]=='R'){
    // ribo-nucleotide
    if(hetatm==0 || exo_check){type=3;}
    else{type=-3;} // Standard nuc and HETATM => ligand
  }else if(exo_check==0){
    exo_check=1;
    // Non-standard residue
    for(int i=0; i<20; i++){
      for(int j=0; j<n_modres[i]; j++){
	if(strncmp(res_type_old, modres[i]+4*j, 3)==0){
	  Code_3_2(amm, stdres[i]); *exo=1; goto find_type;
	}
      }
    }
    if(n_exo){
      for(int i=n_exo-1; i>=0; i--){
	if(strncmp(res_type_old,res_exo[i],3)==0){
	  Code_3_2(amm, res_std[i]); *exo=1; goto find_type;
	}
      }
    }
  }
  //if(type<0)type=-type;
  return(type);
}

////////////////////////////////////////////////////////////////////////////////
int Next_residue(int *N_res, int *start, struct residue *seq,
		 atom *first_atom, short *i_atom,
		 char *res_type_old, int *res_num_old, char *icode_old,
		 char *chain_old, int *hetatm_old, int n_exo,
		 char *res_type, int res_num, char icode, char chain,
		 int hetatm)
{
  if((*start)==0){*start=1;}
  else if(*i_atom)
  	{Write_residue(res_type_old, first_atom, *i_atom, seq+*N_res,
		  n_exo, *res_num_old, *icode_old, *chain_old, *hetatm_old);
     //if(type==1){(*N_res)++;}else{n_atom-=(*i_atom);}
     (*N_res)++;
     (*i_atom)=0;
  	}
  	
  strcpy(res_type_old,res_type);
  *icode_old=icode;
  *res_num_old=res_num;
  *chain_old=chain;
  *hetatm_old=hetatm;

  return(0);
}



////////////////////////////////////////////////////////////////////////////////
void Code_3_2(char *res2, char *res){

  if(strncmp(res,"ALA",3)==0){ sprintf(res2," A");
  }else if(strncmp(res,"GLU",3)==0){ sprintf(res2," E");
  }else if(strncmp(res,"GLN",3)==0){ sprintf(res2," Q");
  }else if(strncmp(res,"ASP",3)==0){ sprintf(res2," D");
  }else if(strncmp(res,"ASN",3)==0){ sprintf(res2," N");
  }else if(strncmp(res,"LEU",3)==0){ sprintf(res2," L");
  }else if(strncmp(res,"GLY",3)==0){ sprintf(res2," G");
  }else if(strncmp(res,"LYS",3)==0){ sprintf(res2," K");
  }else if(strncmp(res,"SER",3)==0){ sprintf(res2," S");
  }else if(strncmp(res,"VAL",3)==0){ sprintf(res2," V");
  }else if(strncmp(res,"ARG",3)==0){ sprintf(res2," R");
  }else if(strncmp(res,"THR",3)==0){ sprintf(res2," T");
  }else if(strncmp(res,"PRO",3)==0){ sprintf(res2," P");
  }else if(strncmp(res,"ILE",3)==0){ sprintf(res2," I");
  }else if(strncmp(res,"MET",3)==0){ sprintf(res2," M");
  }else if(strncmp(res,"PHE",3)==0){ sprintf(res2," F");
  }else if(strncmp(res,"TYR",3)==0){ sprintf(res2," Y");
  }else if(strncmp(res,"CYS",3)==0){ sprintf(res2," C");
  }else if(strncmp(res,"TRP",3)==0){ sprintf(res2," W");
  }else if(strncmp(res,"HIS",3)==0){ sprintf(res2," H");
  }else if(strncmp(res,"ASX",3)==0){ sprintf(res2," N");
  }else if(strncmp(res,"GLX",3)==0){ sprintf(res2," Q");
  }else if(strncmp(res," DG",3)==0){ sprintf(res2,"DG");
  }else if(strncmp(res," DC",3)==0){ sprintf(res2,"DC");
  }else if(strncmp(res," DA",3)==0){ sprintf(res2,"DA");
  }else if(strncmp(res," DT",3)==0){ sprintf(res2,"DT");
  }else if(strncmp(res," DI",3)==0){ sprintf(res2,"DI");
  }else if(strncmp(res,"  G",3)==0){ sprintf(res2,"RG");
  }else if(strncmp(res,"  C",3)==0){ sprintf(res2,"RC");
  }else if(strncmp(res,"  A",3)==0){ sprintf(res2,"RA");
  }else if(strncmp(res,"  U",3)==0){ sprintf(res2,"RT");
  }else if(strncmp(res,"  I",3)==0){ sprintf(res2,"RI");
  }else{ sprintf(res2,"00");
  }
  return;
}

/*int Code_AA(char res){
  short i; char r;
  i=(int)res; if(i>96){r=(char)(i-32);}else{r=res;}
  for(i=0; i<20; i++)if(r==AMIN_CODE[i])return(i);
  if(res=='X')return(0);
  if((res!='-')&&(res!='.')&&(res!='*'))
    printf("Warning, wrong aa type %c\n", res);
  return(-1);
}
*/

FILE *Open_compressed_file(char *name, int *Compression)
{
  FILE *file_in;
  printf("Reading %s ", name);
  *Compression=Get_compression(name);
  char command[200];
  if(*Compression){
    sprintf(command, "%s %s > %s\n", PDBCAT, name, PDBTMP);
    system(command);
    file_in=fopen(PDBTMP, "r");
  }else{
    file_in=fopen(name, "r");
  }
  if(file_in==NULL){
    printf("\nERROR, %s not found\n", name); return(NULL);
  }
  return(file_in);
}

int Get_compression(char *name){
  char *tmp=name;
  while(*tmp!='\0'){
    if((*tmp=='.')&&(*(tmp+1)=='g')&&(*(tmp+2)=='z')){
      return(1);
    }
    tmp++;
  }
  return(0);
}

void Delete_tmp_file(){
  char command[200];
  sprintf(command, "rm -f %s\n", PDBTMP);
  system(command);
}

char Code_3_1(char *res){
  char code;

  if(strncmp(res,"ALA",3)==0){ code='A';
  }else if(strncmp(res,"GLU",3)==0){ code='E';
  }else if(strncmp(res,"GLN",3)==0){ code='Q';
  }else if(strncmp(res,"ASP",3)==0){ code='D';
  }else if(strncmp(res,"ASN",3)==0){ code='N';
  }else if(strncmp(res,"LEU",3)==0){ code='L';
  }else if(strncmp(res,"GLY",3)==0){ code='G';
  }else if(strncmp(res,"LYS",3)==0){ code='K';
  }else if(strncmp(res,"SER",3)==0){ code='S';
  }else if(strncmp(res,"VAL",3)==0){ code='V';
  }else if(strncmp(res,"ARG",3)==0){ code='R';
  }else if(strncmp(res,"THR",3)==0){ code='T';
  }else if(strncmp(res,"PRO",3)==0){ code='P';
  }else if(strncmp(res,"ILE",3)==0){ code='I';
  }else if(strncmp(res,"MET",3)==0){ code='M';
  }else if(strncmp(res,"PHE",3)==0){ code='F';
  }else if(strncmp(res,"TYR",3)==0){ code='Y';
  }else if(strncmp(res,"CYS",3)==0){ code='C';
  }else if(strncmp(res,"TRP",3)==0){ code='W';
  }else if(strncmp(res,"HIS",3)==0){ code='H';
  }else if(strncmp(res,"ASX",3)==0){ code='N';
  }else if(strncmp(res,"GLX",3)==0){ code='Q';
  }else{ return('X');
  }
  return(code);
}

//modyves: parameters were global static variables
void Ini_modres(){
  ini_modres=0;
  int i;
 // for(i=0; i<20; i++){ 									// modyves: no need to allocate given definition on top of file
 //   modres[i]=malloc(100*sizeof(char));
 //   stdres[i]=malloc(4*sizeof(char));
 // }
  strcpy(modres[0], mod_0); strcpy(stdres[0],"ALA");   // A  //modyves: use strcpy
  strcpy(modres[1], mod_1); strcpy(stdres[1],"CYS");   // C 
  strcpy(modres[2], mod_2); strcpy(stdres[2],"ASP");   // D
  strcpy(modres[3], mod_3); strcpy(stdres[3],"GLU");   // E 
  strcpy(modres[4], mod_4); strcpy(stdres[4],"PHE");   // F
  strcpy(modres[5], mod_5); strcpy(stdres[5],"GLY");   // G 
  strcpy(modres[6], mod_6); strcpy(stdres[6],"HIS");   // H
  strcpy(modres[7], mod_7); strcpy(stdres[7],"ILE");   // I
  strcpy(modres[8], mod_8); strcpy(stdres[8],"LYS");   // K 
  strcpy(modres[9], mod_9); strcpy(stdres[9],"LEU");   // L 
  strcpy(modres[10],mod_10); strcpy(stdres[10],"MET"); // M
  strcpy(modres[11],mod_11); strcpy(stdres[11],"ASN"); // N
  strcpy(modres[12],mod_12); strcpy(stdres[12],"PRO"); // P 
  strcpy(modres[13],mod_13); strcpy(stdres[13],"GLN"); // Q 
  strcpy(modres[14],mod_14); strcpy(stdres[14],"ARG"); // R
  strcpy(modres[15],mod_15); strcpy(stdres[15],"SER"); // S
  strcpy(modres[16],mod_16); strcpy(stdres[16],"THR"); // T
  strcpy(modres[17],mod_17); strcpy(stdres[17],"VAL"); // V
  strcpy(modres[18],mod_18); strcpy(stdres[18],"TRP"); // W
  strcpy(modres[19],mod_19); strcpy(stdres[19],"TYR"); // Y
  for(i=0; i<20; i++){
    int k=0; while(modres[i][k]!='\0')k++; n_modres[i]=k/4;
    //printf("%c %d\n", AANAME1[i], n_modres[i]);
  }
}  

////////////////////////////////////////////////////////////////////
int Set_chains(struct chain **chains, atom **atoms, int *N_atoms, struct residue *seq, int *nres)
{
  int Nchain=0, Nc_del=0, ires, ini_r=0, na=0;
  int *sel;
  char chain=seq[0].chain;

  sel=(int*)malloc((*nres)*sizeof(int));
  
  // Eliminate chains with < MIN_ATOM atoms
  for(ires=0; ires<(*nres); ires++)
  	{sel[ires]=1;
     if(chain!=seq[ires].chain)
     	{Nc_del+=Check_chain(&Nchain, sel, &na, ini_r, ires);
         chain=seq[ires].chain;
      	 ini_r=ires; 
    	}
     na+=seq[ires].n_atom;
  	}
  Nc_del+=Check_chain(&Nchain, sel, &na, ini_r, ires);

  // Eliminate residues in empty chains
  int natoms=0,ires2=0;
  struct residue *res1=seq, *res2=seq;
  for(ires=0; ires<(*nres); ires++)
  	{if(sel[ires])
  		{if(ires2!=ires)*res2=*res1;
      	 natoms+=res1->n_atom;
      	 res2++; ires2++;
    	}
     res1++;
  	}
  
  if((*nres!=ires2)||(Nc_del))
  	{printf("%d residues eliminated in %d chains with < %d atoms",*nres-ires2, Nc_del, MIN_ATOM);
  	}
  printf("Storing %d atoms in %d chains\n", natoms, Nchain);

  *nres=ires2;
 // *atoms=malloc(natoms*sizeof(atom)); //modyves: already allocated
  *N_atoms=natoms;
  *chains=malloc(Nchain*sizeof(struct chain)); int i;
  for(i=0; i<Nchain; i++)(*chains)[i].alignres=NULL;

  // ini_atom, ini_res etc.
  //int ichain=-1, resold=-1, nr=0;
  int ichain=0, nr=0;
  struct chain *chp=*chains;
  chain=seq[0].chain;
  natoms=0;
  
  for(ires=0; ires<*nres; ires++){
    struct residue *res=seq+ires;
    if((ires==0)||(res->chain!=chain)){
      chain=res->chain;
      if(nr){Store_chain(chp, &ichain, natoms, &nr); chp++;}
      chp->ini_atom=natoms;
      chp->ini_res=ires;
      chp->type=res->type; // Set type: aa=1, nuc=2,3, other=0
      chp->label=chain;
      chp->res=res;
    }
    nr++;
    atom *atom1=res->atom_ptr, *atom2=*atoms+natoms;
    for(i=0; i<res->n_atom; i++){
      *(atom2)=*(atom1);      
      atom2->res=ires;
      atom2->i_num=natoms+i;
      atom2->chain=ichain;
      atom2->aa=res->amm;
      atom2++; atom1++; 
    }
    res->atom_ptr=*atoms+natoms;
    natoms+=res->n_atom;
  }

  if(nr)Store_chain(chp, &ichain, natoms, &nr);
  printf("Storing %d chains: ", ichain);
  for(i=0; i<ichain; i++)printf("%c",(*chains)[i].label);
  printf("\n");

  *N_atoms=natoms;
  free(sel);sel=NULL; //modyves: this was never freed
  return(ichain);
}

////////////////////////////////////////////////////////////////////
int Check_chain(int *ichain, int *sel, int *na, int ini_r, int nres)
{
  int del;
  if(*na < MIN_ATOM){
    int j; for(j=ini_r; j< nres; j++)sel[j]=0;
    del=1;
  }else{
    (*ichain)++; del=0;
  }
  *na=0;
  return(del);
}

////////////////////////////////////////////////////////////////////
void Store_chain(struct chain *chp, int *ichain, int n, int *nr)
{
  (*ichain)++; 
  chp->natoms=n-chp->ini_atom;
  chp->nres=*nr; 
  *nr=0;
}

////////////////////////////////////////////////////////////////////
void Find_PDB_order(atom *atoms, int N_atoms)
{
  int i, i_next, res=atoms[0].res, ares=0, i2=0, i_min;
  int i_max=N_atoms+1000;
  atom *atom1=atoms, *atom2;
  for(i=0; i<N_atoms; i++){
    atom1=atoms+i;
    if(atom1->res > res){ares=i; res=atom1->res;}
    i_min=i_max; i_next=atom1->i_num+1;
    atom2=atoms+ares; i2=ares;
    while((i2<N_atoms)&&(atom2->res==res)){
      if((atom2->i_num>=i_next)&&(atom2->i_num<i_min)){
	i_min=i2; if(i_min==i_next)break;
      }
      atom2++; i2++;
    }
    if(i_min < i_max){atom1->i_next=i_min;}
    else if(i2<N_atoms){atom1->i_next=i2;}
    else{atom1->i_next=-1;}
  }
}

////////////////////////////////////////////////////////////////////
int Order_atoms(atom *atoms, int natoms, struct residue *seq, int Nres)
{
  int ires;
  for(ires=0; ires<Nres; ires++){
    struct residue *res=seq+ires;
    if((res->type==1)||(res->type==-1)){
      // Store side chain before CO
      Order_atoms_aa(res->atom_ptr, res->n_atom);
    }
  }
  return(1);
}

////////////////////////////////////////////////////////////////////
int Order_atoms_aa(atom *res_atom, int n_atom)
{
  int i_N=-1, i_CA=-1, i_C=-1, i_O=-1, i_OT=-1;
  int i_HN[3], i_HA[3], n_HN=0, n_HA=0, i, j;
  if(HYD)for(j=0; j<3; j++){i_HN[j]=-1; i_HA[j]=-1;}
  atom *atom1=res_atom;
  for(i=0; i<n_atom; i++){
    if(strncmp(atom1->name, "N ", 2)==0){i_N=i;}
    else if(strncmp(atom1->name, "CA", 2)==0){i_CA=i;}
    else if(strncmp(atom1->name, "C ", 2)==0){i_C=i;}
    else if(strncmp(atom1->name, "O ", 2)==0){i_O=i;}
    else if((strncmp(atom1->name, "OT", 2)==0)||
	    (strncmp(atom1->name, "OX", 2)==0)){i_OT=i;}
    else if(atom1->name[0]=='H'){
      if((atom1->name[1]==' ')||(atom1->name[1]=='1')||
	 (atom1->name[1]=='2')||(atom1->name[1]=='3')){
	atom1->name[2]=atom1->name[1]; atom1->name[1]='N';
	i_HN[n_HN]=i; n_HN++;
      }else if(atom1->name[1]=='A'){
	i_HA[n_HA]=i; n_HA++;
      }
    }
    atom1++;
  }
  
  // Store main chain in the order N (HN) CA (HA)
  atom *atom2=atom_tmp;
  int N_atom_res=0;
  N_atom_res+=Copy_atom(&atom2, res_atom, i_N);
  if(HYD){
    for(j=0; j<3; j++)
      N_atom_res+=Copy_atom(&atom2, res_atom, i_HN[j]);
  }
  N_atom_res+=Copy_atom(&atom2, res_atom, i_CA);
  if(HYD){
    for(j=0; j<3; j++)
      N_atom_res+=Copy_atom(&atom2, res_atom, i_HA[j]);
  }
  
  // Store side chain before C-O
  int N_side=n_atom;
  if(i_C>=0)N_side--;
  if(i_O>=0)N_side--; 
  if(N_atom_res==N_side)goto CO;
  N_atom_res+=Store_branch(&atom2, 'B', res_atom, n_atom);
  if(N_atom_res==N_side)goto CO;
  N_atom_res+=Store_branch(&atom2, 'G', res_atom, n_atom);
  if(N_atom_res==N_side)goto CO;
  N_atom_res+=Store_branch(&atom2, 'D', res_atom, n_atom);
  if(N_atom_res==N_side)goto CO;
  N_atom_res+=Store_branch(&atom2, 'E', res_atom, n_atom);
  if(N_atom_res==N_side)goto CO;
  N_atom_res+=Store_branch(&atom2, 'Z', res_atom, n_atom);
  if(N_atom_res==N_side)goto CO;
  N_atom_res+=Store_branch(&atom2, 'H', res_atom, n_atom);
  if(N_atom_res==N_side)goto CO;
  
  // Store other side chain atoms before C-O
  atom1=res_atom;
  for(i=0; i<n_atom; i++){
    if((strncmp(atom1->name, "N ", 2)!=0)&&
       (strncmp(atom1->name, "CA", 2)!=0)&&
       (strncmp(atom1->name, "C ", 2)!=0)&&
       (strncmp(atom1->name, "O ", 2)!=0)&&
       (strncmp(atom1->name, "OT", 2)!=0)&&
       (strncmp(atom1->name, "OX", 2)!=0)&&
       (atom1->name[1]!='N')&&
       (atom1->name[1]!='A')&&
       (atom1->name[1]!='B')&&
       (atom1->name[1]!='G')&&
       (atom1->name[1]!='D')&&
       (atom1->name[1]!='E')&&
       (atom1->name[1]!='Z')&&
       (atom1->name[1]!='H')){
      N_atom_res+=Copy_atom(&atom2, res_atom, i);
    }
    atom1++;
  }
  
  // Store CO
 CO:
  N_atom_res+=Copy_atom(&atom2, res_atom, i_C);
  N_atom_res+=Copy_atom(&atom2, res_atom, i_O);
  N_atom_res+=Copy_atom(&atom2, res_atom, i_OT);

  if(N_atom_res!=n_atom){
    printf("ERROR changing the order of protein atoms\n");
    printf("residue %d %c cain %d has %d residues but %d stored\n",
	   res_atom->res, res_atom->aa, res_atom->chain, n_atom, N_atom_res);
    exit(8);
  }

  for(i=0; i<n_atom; i++){
    res_atom[i]=atom_tmp[i];
  }
  return(N_atom_res);
}

int Store_branch(atom **atom2, char BRANCH,
		 atom *res_atom, int n_atom)
{
  atom *atom1=res_atom; int i, n=0;
  for(i=0; i<n_atom; i++){
    if(atom1->name[1]==BRANCH){
      n+=Copy_atom(atom2, res_atom, i);
    }
    atom1++;
  }
  return(n);
}

int Copy_atom(atom **atom2, atom *atom1, int i)
{
  if(i<0)return(0);
  *(*atom2)=*(atom1+i);
  (*atom2)++;
  return(1);
}

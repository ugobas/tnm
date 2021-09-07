#include "coord.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "read.h"

int Read_seqres(char **seqres, int *N_seqres, int NCHAIN,
		char *file_pdb, char *chain_to_read)
{
  /* Open file */
  int Compression;
  FILE *file_in=Open_compressed_file(file_pdb, &Compression);
  if(file_in==NULL)return(0);

  int nchain=0, read=0, L=0, n=0;
  char string[200], chain_old='-', *seq=NULL;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string,"ATOM", 4)==0)break;
    if(strncmp(string,"SEQRES", 6)==0){
      char chain=string[11];
      if(chain!=chain_old){
	if(n){
	  printf("%d residues read %d expected for chain %d\n",
		 n, L, nchain); n=0;
	}
	chain_old=chain; int ic=0;
	while((chain_to_read[ic]!='\0')&&(chain!=chain_to_read[ic]))ic++;
	if(chain_to_read[ic]!=chain){read=0; continue;}
	read=1; sscanf(string+12, "%d", &L);
	printf("Reading chain %d %c %d residues\n", nchain, chain, L);
	seqres[nchain]=malloc(L*sizeof(char)); seq=seqres[nchain];
	N_seqres[nchain]=L;
	nchain++;
      }
      if(read==0)continue;
      char amm[3], aa3[4], *str=string+19;
      for(int k=0; k<13; k++){
	sscanf(str, "%s", aa3); str+=4;
	int exo=0; Get_amm(amm, &exo, aa3, 0, 0);
	*seq=amm[1]; seq++; n++; if(n==L)break;
      }
    }
  }
  if(n){
    printf("%d residues read %d expected for chain %d\n",
	   n, L, nchain); n=0;
  }
  if(nchain!=NCHAIN){
    printf("WARNING in SEQRES, %d chains expected (%s), %d found\n",
	   NCHAIN, chain_to_read, nchain);
    if(nchain==0){printf("Getting residues from atoms\n"); return(-1);}
    else{printf("Exiting\n"); exit(8);}
  }
  return(0);
}

void Align_seqres(int *ali_seqres, char *seqres, int N_seqres,
		  char *seq, int nres)
{
  int IDMIN=7;
  int i=0, k=0; char *s1=seqres, *s2=seq;
  for(i=0; i<N_seqres; i++){
    if((*s1==*s2)&&(k<nres)){
      int ali=1;
      if((i==0)||(ali_seqres[i-1]<0)){
	// gap, look forward
	for(int i1=i; i1<N_seqres; i1++){
	  char *t1=seqres+i1, *t2=s2; int id=0;
	  while(*t2==*t1){
	    id++; if(id==IDMIN)break; t1++; t2++;
	  }
	  if(id==IDMIN){if(i1!=i)ali=0; break;}
	}
      }
      if(ali){ali_seqres[i]=k; k++; s2++;} // Align here
      else{ali_seqres[i]=-1;} // Align forward, not here
    }else{
      ali_seqres[i]=-1;
    }
    s1++;
  }
  if(0){
    printf("%d aligned residues over %d disordered: %d\n",
	   k, N_seqres, N_seqres-k);
    printf("SEQRES:  ");
    for(i=0; i<N_seqres; i++)printf("%c", seqres[i]); printf("\n");
    printf("ATOMRES: ");
    for(i=0; i<N_seqres; i++){
      if(ali_seqres[i]<0){printf("-");}
      else{printf("%c", seq[ali_seqres[i]]);}
    }
    printf("\n");
    exit(8);
  }
}

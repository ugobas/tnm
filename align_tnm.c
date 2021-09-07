#include "coord.h"
#include "tnm.h"
#include "align_tnm.h"
#include "NeedlemanWunsch.h"
#include "mammoth_interface.h"
//#include "needle.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define DEBUG 0

int PDB_ali(int *seq_id, int nali,
	    char *ali1, struct chain *chp1, int n1,
	    char *ali2, struct chain *chp2, int n2);
int Max_length(struct chain *chains, int Nchain);
void Match_chains_by_length(int *match1,
			    struct chain *chains1, int Nchain1,
			    struct chain *chains2, int Nchain2);
void Best_match(int *match1, int *free2,
		struct chain *chains1, int Nchain1,
		struct chain *chains2, int Nchain2);

int Align_chains(int *Nmut, char *AAwt, int *Posmut, char *AAmut,
		 struct Ali_score *ali, int *last_ali_res,
		 int *Nchain, float SEQID_THR, 
		 struct chain *chains1, int Nchain1,
		 struct residue *seq1, char *pdbid1,
		 struct chain *chains2, int Nchain2,
		 struct residue *seq2, char *pdbid2)
{
  // Sequence alignment parameters
  int VBS=0;  // Verbose:
  int IDE=1;  // Use identity to score alignment
  int GAP=7;  // Gap opening penalty

  *Nchain=Nchain1;
  if(Nchain2!=Nchain1){
    if(Nchain2<*Nchain)*Nchain=Nchain2;
    printf("WARNING, different number of chains: ");
    printf("%d (%s) and %d (%s)\n", Nchain1, pdbid1, Nchain2, pdbid2);
    printf("Considering only the %d best matching chains\n", *Nchain);
    //exit(8);
  }

  ali->mammoth=0;
  ali->aligned=0;
  ali->ngaps=0;
  ali->seq_id=0;
  *Nmut=0;

  int n2max=Max_length(chains2, Nchain2);
  int matched[Nchain2]; for(int i=0; i<Nchain2; i++)matched[i]=0;

  //int match_chain[Nchain1];
  //Match_chains_by_length(match_chain, chains1, Nchain1, chains2, Nchain2);

  int NMUT=1000;
  int nres1=0, nres2=0;
  int ichain, jchain, i;
  for(ichain=0; ichain< Nchain1; ichain++){
    struct chain *chp1=chains1+ichain;
    int n1=chp1->nres;
    chp1->alignres=malloc(n1*sizeof(int));
    for(i=0; i<n1; i++)chp1->alignres[i]=-1;
    int jopt=-1, dL_opt=1000, seq_id_opt=0, nali_opt=0;
    char ali1_opt[(n1+n2max)], ali2_opt[(n1+n2max)];

    for(jchain=0; jchain< Nchain2; jchain++){
      if(matched[jchain])continue;
      struct chain *chp2=chains2+jchain;
      int n2=chp2->nres;

      printf("Aligning chains %d %c prot1 L=%d and %d %c prot2 L=%d\n",
	     ichain, chp1->label, n1, jchain, chp2->label, n2);

      // Make alignment
      //for(i=0; i<n1; i++)aligntmp[i]=-1;
      //float psi;
      /*if(SEQID_THR > 1.0){
	ali->mamm_score=
	Mammoth_ali(&aligned, &seq_id, &psi, chp1->alignres,
	chp1, seq1, chp2, seq2);
	ali->mammoth=1;
	ali->psi+=sqrt(n1*n2)*psi;
	goto update_chain;
	}*/	
      // Alignment based on sequence
      //ali->seq_id=
      //NeedlemanWunsch(chp1->seq,chp2->seq,n1,n2,"BLOSUM62",chp1->alignres);

      char ali1[(n1+n2)], ali2[(n1+n2)];
      int seq_id=0, nali=0;
      int al=alignNW(chp1->seqres, chp1->N_seqres,
		     chp2->seqres, chp2->N_seqres,
		     VBS,IDE,GAP, ali1, ali2, &nali);
      if(al==0){printf("ERROR, alignment failed\n"); exit(8);}

      nali=PDB_ali(&seq_id, nali, ali1, chp1, n1, ali2, chp2, n2);
      printf("Seq.id= %d %d total residues\n", seq_id, nali);

      // Choose matching chain based on alignment or length
      if(seq_id > 0.3){
	if(seq_id > seq_id_opt){
	  seq_id_opt=seq_id; nali_opt=nali;
	  jopt=jchain;
	  for(i=0; i<nali; i++){
	    ali1_opt[i]=ali1[i];
	    ali2_opt[i]=ali2[i];
	  }
	  if(seq_id_opt > 0.9*n1)break;
	}
      }else if(jopt<0){
	int dL=abs(n1-n2);
	if(dL < dL_opt){jopt=jchain; dL_opt=dL;}
      }
    } // end jchain

    if(jopt<0){
      printf("ERROR, could not match chain %d\n", ichain);
      break;
    }

    matched[jopt]=1;
    chp1->match=jopt;
    struct chain *chp2=chains2+jopt;
    int n2=chp2->nres;
    chp2->match=ichain;
    nres1+=n1; nres2+=n2;

    if((n1 < MIN_ALI_LEN)||(n2 < MIN_ALI_LEN)){
      // Trivial alignment for peptides or nucleic acids
      nali_opt=n1; if(n2>nali_opt)nali_opt=n2;
      for(i=0; i<nali_opt; i++){
	if(i<n1){ali1_opt[i]=chp1->seq[i];}
	else{ali1_opt[i]='-';}
	if(i<n2){ali2_opt[i]=chp2->seq[i];}
	else{ali2_opt[i]='-';}
      }
    }

    // Print alignment
    printf("Protein_1: ");
    for(i=0; i<nali_opt; i++)printf("%c", ali1_opt[i]); printf("\n");
    printf("Protein_2: ");
    for(i=0; i<nali_opt; i++)printf("%c", ali2_opt[i]); printf("\n");

    // Store alignment
    int seq_id=0, aligned=0;
    int i1=0, i2=0, gap=1;
    for(i=0; i<nali_opt; i++){
      if((ali1_opt[i]!='-')||(ali2_opt[i]!='-')){
	if((ali1_opt[i]==ali2_opt[i])&&(ali1_opt[i]!='-')){
	  seq_id++;
	}else{
	  Posmut[*Nmut]=i1; AAmut[*Nmut]=ali2_opt[i];
	  if(*Nmut<NMUT)AAwt[*Nmut]=ali1_opt[i];
	  (*Nmut)++;
	} 
	aligned++;
      }
      if(ali1_opt[i]!='-'){
	if(gap)gap=0;
	if(ali2_opt[i]!='-'){chp1->alignres[i1]=i2;}
	else{chp1->alignres[i1]=-1;}
	i1++; 
      }else if(gap==0){
	ali->ngaps++; gap=1;
      }
      if(ali2_opt[i]!='-')i2++;
    }

      /* float ss=(float)seq_id/(aligned);
       if(ss < SEQID_THR){
       printf("WARNING! Low sequence identity %.1f\%\n", ss*100);
       printf("Performing structure alignment with Mammoth\n");
       ali->mamm_score=
       Mammoth_ali(&aligned, &seq_id, &psi, chp1->alignres,
       chp1, seq1, chp2, seq2);
       ali->mammoth=1;
       ali->psi+=sqrt(n1*n2)*psi;
       }*/

    if(gap)ali->ngaps--; 
    ali->seq_id+=seq_id;
    ali->aligned+=aligned;

  }
  printf("%d mutations: ", *Nmut);
  for(i=0; i<*Nmut; i++)
    if(*Nmut<NMUT)printf(" %c%s%c", AAwt[i],seq1[Posmut[i]].pdbres,AAmut[i]);
  printf("\n");
  
  ali->seq_id/=ali->aligned; ali->seq_id*=100;
  float norm=sqrt((float)nres1*nres2);
  ali->psi/=norm;
  //ali->aligned/=norm;
  ali->aligned/=nres1;

  ali->alignres=malloc(nres1*sizeof(int)); i=0;
  for(ichain=0; ichain< *Nchain; ichain++){
    struct chain *chp1=chains1+ichain;
    int n1=chp1->nres, i1, ini2=(chains2+chp1->match)->ini_res;
    for(i1=0; i1<n1; i1++){
      if(chp1->alignres[i1]>=0){
	ali->alignres[i]=chp1->alignres[i1]+ini2;
      }else{
	ali->alignres[i]=-1;
      }
      i++;
    }
  }

  // last aligned residue
  for(i=nres1-1; i>=0; i--){
    if(ali->alignres[i]>=0){
      *last_ali_res=i+1; break;
    }
  }
  return 0;
}

void Match_chains_by_length(int *match1,
			    struct chain *chains1, int Nchain1,
			    struct chain *chains2, int Nchain2)
{
  int match2[Nchain2], free1[Nchain1], free2[Nchain2], i, j;
  for(i=0; i<Nchain1; i++){match1[i]=-1; free1[i]=1;}
  for(i=0; i<Nchain2; i++){match2[i]=-1; free2[i]=1;}
  Best_match(match1, free1, chains1, Nchain1, chains2, Nchain2);
  Best_match(match2, free2, chains2, Nchain2, chains1, Nchain1);
  for(i=0; i<Nchain1; i++){
    j=match1[i];
    if(match2[j]==i){free2[j]=0;} // Freeze reciprocal best matches
    else{match1[i]=-1;}
  }
  Best_match(match1, free1, chains1, Nchain1, chains2, Nchain2);  

}

void Best_match(int *match1, int *free2,
		struct chain *chains1, int Nchain1,
		struct chain *chains2, int Nchain2)
{
  for(int i=0; i<Nchain1; i++){
    if(match1[i]>=0)continue;
    int imin=-1, d_min=100000, n1=chains1[i].nres;
    if((i<Nchain2)&&(free2[i])){
      d_min=abs(n1-chains2[i].nres); imin=i;
    }
    for(int j=0; j<Nchain2; j++){
      if((j==i)||(free2[j]==0))continue;
      int d=abs(n1-chains2[j].nres);
      if(d<d_min){d_min=d; imin=j;}
    }
    match1[i]=imin;
  }
}

void Purge_not_aligned(struct chain *chains1, atom *atoms1, int *natoms1,
		       struct chain *chains2, atom *atoms2, int *natoms2)
{
  // Not aligned atoms are purged, except side chain atoms in aligned residues
  int k1=0, k2=0, j_old=0, i, j, ichain=-1;
  int notali=0, notalimax=1000; atom notali_atom[notalimax];
  printf("Not aligned atoms: ");
  struct chain *ch1=NULL, *ch2=NULL;
  for(i=0; i<*natoms1; i++){
    atom *a1=atoms1+i;
    if(a1->chain!=ichain){
      ichain=a1->chain;
      ch1=chains1+ichain; ch1->ini_atom=-1;
      ch2=chains2+ch1->match; ch2->ini_atom=-1;
      printf("Chain %d ",ichain);
    }
    j=a1->ali;
    if(j>=0){
      // if((j>=0)||
      //  ((ch1->alignres[ch1->inires+a1->res]>=0)&&(a1->main==0))){ //
      // aligned atom or side chain atom in aligned residue
      if(ch1->ini_atom<0)ch1->ini_atom=k1;
      if(i!=k1)atoms1[k1]=*a1;
      // if(j>=0){
      if(ch2->ini_atom<0)ch2->ini_atom=k2;
      if(j!=k2){
	for(int k=j_old+1; k<j; k++){
	  if(notali>=notalimax)break;
	  notali_atom[notali]=atoms2[k]; notali++;
	}
	atoms2[k2]=atoms2[j];
      }
      atoms1[k1].ali=k2; j_old=j; k2++;
      // }else{printf("Saved: %s/%c%d", a1->name, a1->aa, a1->res);}
      k1++;
    }else{
      // not aligned
      printf(" %s/%c%d", a1->name, a1->aa, a1->res);
      ch1->natoms--;
    }
  }
  printf("\nn1= %d n2= %d\n", k1, k2);
  if(k1 != *natoms1){
    printf("WARNING, eliminating %d not aligned atoms of prot 1 (%d)\n",
	   *natoms1-k1, *natoms1); *natoms1=k1;
  }
  if(k2 != *natoms2){
    printf("WARNING, eliminating %d not aligned atoms of prot 2 (%d)\n",
	   *natoms2-k2, *natoms2); *natoms2=k2;
    atom *a1=notali_atom;
    for(int k=0; k<notali; k++){
      printf(" %s/%c%d/%d", a1->name, a1->aa, a1->res, a1->chain); a1++;
    }
    printf("\n");
  }
}

int Align_atoms(atom *atoms1, atom *atoms2,          // Output
		struct chain *ch1, struct chain *ch2) // Input
{
  int ini1=ch1->ini_atom, natoms1=ch1->natoms, ini_res1=ch1->ini_res;
  int ini2=ch2->ini_atom, natoms2=ch2->natoms, ini_res2=ch2->ini_res;
  int *alignres=ch1->alignres;
  int n=0, i, i1, i2max=ini2+natoms2;
  printf("Aligning atoms in chain %d %c %d residues and %d %c %d residues\n",
	 ch1->match,ch1->label,ch1->nres, ch2->match,ch2->label,ch2->nres);

  if(0){
    for(i1=0; i1<ch1->nres; i1++){
      printf("%c ", ch1->seq[i1]);
      int i2=ch1->alignres[i1];
      if(i2>=0){printf("%c\n", ch2->seq[i2]);}
      else{printf("-\n");}
    }
    printf("chains %c %c\n", ch1->label, ch2->label);
    exit(8);
  }

  // Initialize
  for(i=ini1; i<(ini1+natoms1); i++)atoms1[i].ali=-1;
  for(i=ini2; i<(ini2+natoms2); i++)atoms2[i].ali=-1;
  // Look if aligned
  int ares2=ini2;
  atom *atom1;
  for(i1=ini1; i1<(ini1+natoms1); i1++){
    atom1=atoms1+i1;
    int res2=alignres[atom1->res-ini_res1];
    if(res2<0)continue; res2+=ini_res2;
    int iali=-1, i2=ares2;
    atom *atom2=atoms2+ares2; ;
    while((atom2->res<res2)&&(i2<i2max)){atom2++; i2++;}
    ares2=i2;
    while((i2<i2max)&&(atom2->res==res2)){
      if((atom2->ali<0)&&(atom2->name[1]==atom1->name[1])){
	if((atom2->name[0]==atom1->name[0])&&
	   (atom2->name[2]==atom1->name[2])){
	  atom1->ali=i2; atom2->ali=i1; n++; break; // Same name
	}else if(iali<0){
	  iali=i2; // Same side-chain index
	}
      }
      i2++; atom2++;
    }
    if(atom1->ali>=0)continue;
    // For main chain atoms we require that the name coincides exactly
    // For side chain atoms, the same side-chain index is enough
    atom1->main=Is_atom_main(atom1->name, ch1->type);
    if((atom1->main==0)&&(iali>=0)){
      atom1->ali=iali; (atoms2+iali)->ali=i1; n++; continue;
    }
    // Not aligned atom
    if(iali>=0){atom2=atoms2+iali;}
    else{atom2=atoms2+ares2;}
    printf("Not aligned: 1: %s%d%c %c 2: %s%d%c %c",
	   atom1->name, atom1->res, atom1->aa, ch1->label,
	   atom2->name, atom2->res, atom2->aa, ch2->label);
    printf(" Aligned res: %d iali: %d\n", ares2, iali);
  }

  if(DEBUG){
    printf("%d Aligned atoms:\n", n);
    for(i=ini1; i<(ini1+30); i++){
      atom *atom1=atoms1+i;
      if(atom1->ali<0)continue;
      atom *atom2=atoms2+atom1->ali;
      printf("%s %3d %d - %s %3d %d\n",atom1->name,atom1->res,i,
	     atom2->name,atom2->res, atom1->ali);
    }
    //exit(8);
    }
  return(n);
}

int PDB_ali(int *seq_id, int nali,
	    char *ali1, struct chain *chp1, int n1,
	    char *ali2, struct chain *chp2, int n2)
{
  int i1=0, i2=0, j1, j2, k=0; *seq_id=0;
  for(int i=0; i<nali; i++){
    // seq.1
    if(ali1[i]!='-'){j1=chp1->ali_seqres[i1]; i1++;}
    else{j1=-1;}
    if(j1>=0){ali1[k]=chp1->seq[j1];}
    else{ali1[k]='-';}
    // seq.2
    if(ali2[i]!='-'){j2=chp2->ali_seqres[i2]; i2++;}
    else{j2=-1;}
    if(j2>=0){ali2[k]=chp2->seq[j2];}
    else{ali2[k]='-';}
    if((ali1[k]!='-')&&(ali1[k]==ali2[k]))(*seq_id)++;
    if((j1>=0)||(j2>=0)){k++;}
  }
  return(k);
}

int Max_length(struct chain *chains, int Nchain)
{
  int nmax=0;
  for(int i=0; i<Nchain; i++)if(chains[i].nres>nmax)nmax=chains[i].nres;
  return(nmax);
}

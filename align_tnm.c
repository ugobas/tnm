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

int Invert_dimer=1; // Align A with B and B with A in dimers


static void Superimpose_cluster(struct chain *chains1, int Nchain1,
				struct chain *chains2, int Nchain2,
				int *i1cl, int n1cl, int *i2cl, int n2cl,
				int ***ali_all, int **nali);
static void Superimpose_centers(int **ichain2, int ncompl,
				struct chain *chains1, int Nchain1,
				struct chain *chains2, int Nchain2,
				int *i1cl, int *i2cl, 
				int ***ali_all, int **nali);
static float Superimpose_complex(int **ichain2, int ncompl,
				 struct chain *chains1, int Nchain1,
				 struct chain *chains2, int Nchain2,
				 int *i1cl, int *i2cl, 
				 int ***ali_all, int **nali);

static void Test_chains(int *icl, int ncl,
			struct chain *chains, int Nchain, int n);

static int PDB_ali(int *ali, int *seq_id, int nali,
	    char *ali1, struct chain *chp1, int n1,
	    char *ali2, struct chain *chp2, int n2);
static int Max_length(struct chain *chains, int Nchain);
void Match_chains_by_length(int *match1,
			    struct chain *chains1, int Nchain1,
			    struct chain *chains2, int Nchain2);
static void Best_match(int *match1, int *free2,
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
  float seqid_min=0.4, si_tol=0.05; // seqid_max=0.9;

  // Sequence alignment parameters
  int VBS=0;  // Verbose:
  int IDE=1;  // Use identity to score alignment
  int GAP=7;  // Gap opening penalty

  *Nchain=Nchain1;
  if(Nchain2<Nchain1){
    printf("WARNING, different number of chains: ");
    printf("%d (%s) and %d (%s)\n", Nchain1, pdbid1, Nchain2, pdbid2);
    if(Nchain2<Nchain1){printf("Exiting the program\n"); exit(8);}
    printf("Considering only the %d best matching chains\n", *Nchain);
  }

  ali->mammoth=0;
  ali->aligned=0;
  ali->ngaps=0;
  ali->seq_id=0;

  //Cluster chains_by_sequence identity
  // Alignments
  int nmax=Max_length(chains1, Nchain1)+Max_length(chains2, Nchain2);
  int **ali_all[Nchain1], *nali[Nchain1];
  int seqid[Nchain1][Nchain2];

  // Clusters
  int ic1, ic2, i;
  int ncl=0, n1cl[Nchain1], n2cl[Nchain1],*i1cl[Nchain1], *i2cl[Nchain1];
  for(ic1=0; ic1<Nchain1; ic1++){
    ali_all[ic1]=malloc(Nchain2*sizeof(int *));
    nali[ic1]=malloc(Nchain2*sizeof(int));
    i1cl[ic1]=malloc(Nchain1*sizeof(int));
    i2cl[ic1]=malloc(Nchain2*sizeof(int));
  }
  int cl2[Nchain2]; //cl1[Nchain1], 
  for(i=0; i<Nchain1; i++){n1cl[i]=0; n2cl[i]=0;}  //cl1[i]=-1; 
  for(i=0; i<Nchain2; i++)cl2[i]=-1;

  for(ic1=0; ic1< Nchain1; ic1++){
    struct chain *chp1=chains1+ic1;
    int n1=chp1->nres;
    chp1->match=-1;

    int seqid_opt=0, ic2_opt=-1;      
    for(ic2=0; ic2< Nchain2; ic2++){
      struct chain *chp2=chains2+ic2;
      int n2=chp2->nres, na;

      printf("Aligning chains %d %c prot1 L=%d and %d %c prot2 L=%d\n",
	     ic1, chp1->label, n1, ic2, chp2->label, n2);

      int nmax_seqres=chp1->N_seqres+chp1->N_seqres;
      char ali1[nmax_seqres], ali2[nmax_seqres];

      if(n1 <= MIN_ALI_LEN){
	// Trivial alignment for peptides or nucleic acids
	printf("Trivial alignment for peptide or nucleic acid\n");
	na=n1;
	for(i=0; i<n1; i++){
	  ali1[i]=chp1->seq[i];
	  if(i<n2){ali2[i]=chp2->seq[i];}
	  else{ali2[i]='-';}
	}
      }else{
	int al=alignNW(chp1->seqres, chp1->N_seqres,
		       chp2->seqres, chp2->N_seqres,
		       VBS,IDE,GAP, ali1, ali2, &na);
	if(al==0){printf("ERROR, alignment failed\n"); exit(8);}
      }

      // Transform the alignment of seqres into an alignment of
      // structured residues (seq)
      ali_all[ic1][ic2]=malloc(nmax*sizeof(int));
      nali[ic1][ic2]=PDB_ali(ali_all[ic1][ic2], &(seqid[ic1][ic2]),
			     na, ali1, chp1, n1, ali2, chp2, n2);
 
      printf("Seq.id= %d alignment length %d\n",
	     seqid[ic1][ic2],nali[ic1][ic2]);

      // Best matching chain
      if(seqid[ic1][ic2] > seqid_opt){
	seqid_opt=seqid[ic1][ic2]; ic2_opt=ic2;
      }
    } // end ic2

    printf("Best match: chain %d seqid= %d n=%d chain %d",
	   ic1, chp1->nres, seqid_opt, ic2_opt);
    if(ic2_opt>=0)printf(" n=%d\n", chains2[ic2_opt].nres);
    if(ic2_opt<0 || (seqid_opt<seqid_min && n1>MIN_ALI_LEN)){
      // If short chain (peptide) accept the match even for low identity
      printf("\nERROR, could not reliably match chain %d\n", ic1);
      exit(8);
    }

    /* float ss=(float)seqid_opt/n1;
       if(ss < SEQID_THR){
       printf("WARNING! Low sequence identity %.1f\%\n", ss*100);
       printf("Performing structure alignment with Mammoth\n");
       ali->mamm_score=
       Mammoth_ali(&aligned, &seqid_opt, &psi, chp1->alignres,
       chp1, seq1, chp2, seq2);
       ali->mammoth=1;
       ali->psi+=sqrt(n1*n2)*psi;
       }*/
    
    int icl, seqid_thr;
    if(cl2[ic2_opt]>=0){icl=cl2[ic2_opt];} // already clustered
    else{icl=ncl; ncl++;}
    if(seqid_opt>=seqid_min){seqid_thr=seqid_opt-si_tol;}
    else{seqid_thr=seqid_min;}

    // Set clusters
    i1cl[icl][n1cl[icl]]=ic1; n1cl[icl]++; //cl1[ic1]=icl; 
    for(ic2=0; ic2< Nchain2; ic2++){
      if((seqid[ic1][ic2]>=seqid_thr)&&(cl2[ic2]!=icl)){
	cl2[ic2]=icl; i2cl[icl][n2cl[icl]]=ic2; n2cl[icl]++;
      }
    }
  } // end ic1

  if(ncl==0){
    printf("ERROR, no sequence clusters could be found\n"); exit(8);
  }

  for(int ic=0; ic<ncl; ic++){
    Superimpose_cluster(chains1, Nchain1, chains2, Nchain2,
			i1cl[ic], n1cl[ic], i2cl[ic], n2cl[ic],
			ali_all, nali);
  }


  *Nmut=0;
  int nres1=0, nres2=0;
  // Store alignments
  for(ic1=0; ic1<Nchain1; ic1++){
    struct chain *chp1=chains1+ic1;
    int n1=chp1->nres;
    ic2=chp1->match;
    if(ic2<0){
      printf("ERROR, chain %d %c n=%d could not be superimposed\n",
	     ic1, chp1->label, n1); exit(8);
    }
    
    chp1->alignres=malloc(n1*sizeof(int));
   
    struct chain *chp2=chains2+ic2;
    int n2=chp2->nres;
    chp2->match=ic1;
    nres1+=n1; nres2+=n2;
    int *al=ali_all[ic1][ic2];
   
    // Print alignment
    printf("Protein_1: ");
    for(i=0; i<n1; i++)printf("%c", chp1->seq[i]);
    printf("\nProtein_2: ");
    for(i=0; i<n1; i++){
      if(al[i]>=0){printf("%c", chp2->seq[al[i]]);}
      else{printf("-");}
    }
    printf("\n");
    
    // Store alignment, count mutations
    int seq_id=0, ngap=0, gap=1;
    printf("Alignment: ");
    for(i=0; i<n1; i++){
      if(al[i]>=0){
	if(gap)gap=0;
	if(chp1->seq[i]==chp2->seq[al[i]]){
	  seq_id++;
	}else{ // Mutation
	  Posmut[*Nmut]=i;
	  AAwt[*Nmut]=chp1->seq[i];
	  AAmut[*Nmut]=chp2->seq[al[i]];
	  (*Nmut)++;
	} 
	chp1->alignres[i]=al[i];
	printf("%c",chp2->seq[al[i]]);
      }else{
	chp1->alignres[i]=-1;
	printf("-");
	ngap++;
	if(gap==0){ali->ngaps++; gap=1;}
      }
    }
    printf("\n");

    if(gap)ali->ngaps--; 
    ali->seq_id+=seq_id;
    int aligned=n1-ngap;
    ali->aligned+=aligned;
  } // end ic1

  printf("%d mutations: ", *Nmut); int NMUT=1000;
  for(i=0; i<*Nmut; i++){
    if(*Nmut==NMUT)break;
    printf(" %c%s%c", AAwt[i],seq1[Posmut[i]].pdbres,AAmut[i]);
  }
  printf("\n");
  
  ali->seq_id/=ali->aligned; ali->seq_id*=100;
  float norm=sqrt((float)nres1*nres2);
  ali->psi/=norm;
  //ali->aligned/=norm;
  ali->aligned/=nres1;

  ali->alignres=malloc(nres1*sizeof(int)); i=0;
  for(ic1=0; ic1< *Nchain; ic1++){
    struct chain *chp1=chains1+ic1;
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

  for(ic1=0; ic1<Nchain1; ic1++){
    for(i=0; i<Nchain2; i++)free(ali_all[ic1][i]);
    free(ali_all[ic1]);
    free(nali[ic1]);
    free(i1cl[ic1]);
    free(i2cl[ic1]);
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
    if(res2<0){continue;} res2+=ini_res2;
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


int Max_length(struct chain *chains, int Nchain)
{
  int nmax=0;
  for(int i=0; i<Nchain; i++)if(chains[i].nres>nmax)nmax=chains[i].nres;
  return(nmax);
}

int PDB_ali(int *ali, int *seq_id, int nali,
	    char *ali1, struct chain *chp1, int n1,
	    char *ali2, struct chain *chp2, int n2)
{
  int i1=0, i2=0, j1, j2, k=0; *seq_id=0;
  for(int i=0; i<n1; i++)ali[i]=-1;
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
    if(j1>=0){
      ali[j1]=j2;
      if(ali1[k]==ali2[k])(*seq_id)++;
    }
    if((j1>=0)||(j2>=0)){k++;}
  }
  return(k);
}

void Superimpose_cluster(struct chain *chains1, int Nchain1,
			 struct chain *chains2, int Nchain2,
			 int *i1cl, int n1cl, int *i2cl, int n2cl,
			 int ***ali_all, int **nali)
{
  if(n1cl!=n2cl){
    printf("WARNING, different number of chains in clusters of PDB1 (%d) "
	   "and PDB2 (%d)\n", n1cl, n2cl);
    if(n1cl!=1){printf("Exiting\n"); exit(8);}
  }
  Test_chains(i1cl, n1cl, chains1, Nchain1, 1);
  Test_chains(i2cl, n2cl, chains2, Nchain2, 2);

  if(n1cl==1){chains1[i1cl[0]].match=i2cl[0]; return;}
  if(Invert_dimer && n1cl==2 && n2cl==2){
    chains1[i1cl[0]].match=i2cl[1];
    chains1[i1cl[1]].match=i2cl[0]; return;
  }
  for(int i=0; i<n1cl; i++)chains1[i1cl[i]].match=i2cl[i];
  return;

  int ncompl=n1cl, *ichain2[ncompl], i;
  for(i=0; i<ncompl; i++)ichain2[i]=malloc(ncompl*sizeof(int));

  if(ncompl==2 || ncompl>=7){
    // Circular permutations
    for(i=0; i<ncompl; i++){
      int ij=i, j;
      for(j=0; j<ncompl; j++){
	ichain2[i][j]=i2cl[ij]; ij++; if(ij==ncompl)ij=0;
	if(ichain2[i][j]<0 || ichain2[i][j]>=Nchain2){
	  printf("ERROR in circular perm., wrong chain %d\n",
		 ichain2[i][j]); exit(8);
	}
      }
    }
  }else{
    Superimpose_centers(ichain2, ncompl, chains1, Nchain1, chains2, Nchain2,
			i1cl, i2cl, ali_all, nali);
  }

  // Find largest RMSD
  float RMSD_max=-1; int imax=-1;
  for(i=0; i<ncompl; i++){
    float RMSD=Superimpose_complex(ichain2, ncompl, chains1, Nchain1,
				   chains2,Nchain2,i1cl,i2cl, ali_all, nali);
    if(RMSD>RMSD_max){RMSD_max=RMSD; imax=i;}
  }
  if(imax<0){
    printf("ERROR, optimal superimposition could not be found\n");
    exit(8);
  }

  for(i=0; i<n1cl; i++){
    chains1[i1cl[i]].match=ichain2[imax][i];
  }
  for(i=0; i<ncompl; i++)free(ichain2[i]);

  return;
}

void Test_chains(int *icl, int ncl,struct chain *chains, int Nchain, int n)
{
  int i, j;
  for(j=0; j<ncl; j++){
    i=icl[j];
    if(i<0 || i>=Nchain){
      printf("ERROR, wrong index of chain in PDB%d: %d not in [0,%d]\n",
	     n, i, Nchain); exit(8);
    }
  }
}

void Superimpose_centers(int **ichain2, int ncompl,
			 struct chain *chains1, int Nchain1,
			 struct chain *chains2, int Nchain2,
			 int *i1cl, int *i2cl, 
			 int ***ali_all, int **nali)
{

}

float Superimpose_complex(int **ichain2, int ncompl,
			  struct chain *chains1, int Nchain1,
			  struct chain *chains2, int Nchain2,
			  int *i1cl, int *i2cl, 
			  int ***ali_all, int **nali)
{
  float RMSD=0;
  return(RMSD);
}

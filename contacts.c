#include "nma_para.h"
#include "coord.h"
#include "contacts.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static void Write_contact_list(short **cont_i, short **cont_j,
			       struct interaction *Int_list, int *N_int);

float Weighted_overlap(struct interaction *Int_list1, int N_int1,
		       struct interaction *Int_list2, int N_int2)
{
  double overlap=0;

  return(overlap);
}

float Contact_overlap(struct interaction *Int_list1, int N_int1,
		      struct interaction *Int_list2, int N_int2, int *alires)
{
  // Make contact matrix 1
  short *cont1_i, *cont1_j; int Nc1=N_int1;
  Write_contact_list(&cont1_i, &cont1_j, Int_list1, &Nc1);
  short *cont2_i, *cont2_j; int Nc2=N_int2;
  Write_contact_list(&cont2_i, &cont2_j, Int_list2, &Nc2);
  printf("Initial number of contacts: 1 %d 2 %d\n", N_int1, N_int2);
  printf("After identifying residues: 1 %d 2 %d\n", Nc1, Nc2);
  int k;
  for(k=0; k<Nc1; k++){
    cont1_i[k]=alires[cont1_i[k]];
    cont1_j[k]=alires[cont1_j[k]];
  }
  int overlap=0, k2=0;
  for(k=0; k<Nc1; k++){
    if((cont1_i[k]<0)||(cont1_j[k]<0))continue;
    while((k2<Nc2)&&(cont2_i[k2] < cont1_i[k]))k2++;
    while((k2<Nc2)&&(cont2_i[k2]== cont1_i[k])&&
	  (cont2_j[k2] < cont1_j[k]))k2++;
    if((cont1_i[k]==cont2_i[k2])&&(cont1_j[k]==cont2_j[k2])){
      overlap++; k2++;
    }
  }
  free(cont1_i);free(cont1_j); //modyves: was never freed
  free(cont2_i);free(cont2_j); //modyves: was never freed
  return((float)overlap/sqrt((float)(Nc1*Nc2)));
}

float Contact_energy(struct interaction *Int_list, int N_int,
		     struct residue *seq){
  double E=0; int i, ai, aj, res=-1;
  for(i=0; i<N_int; i++){
    ai=Code_AA(seq[Int_list[i].res1].amm);
    if(ai<0){
      if(res!=Int_list[i].res1){
        res=Int_list[i].res1;
	printf("WARNING, unknown a.a. %c residue %d\n",seq[res].amm, res);
      }
      continue;
    }
    aj=Code_AA(seq[Int_list[i].res2].amm);
    if(aj<0){
      if(res!=Int_list[i].res2){
        res=Int_list[i].res2;
	printf("WARNING, unknown a.a. %c residue %d\n",seq[res].amm, res);
      }
      continue;
    }
    E+=Econt[ai][aj];
  }
  return(E);
}

void Get_contact_list(int **nc1, int ***clist1, int ***cnum1, int Nres,
		      atom *atoms, struct interaction *Int_list, int N_int)
{
  int **clist=(int **)malloc(Nres*sizeof(int *));
  int **cnum=(int **)malloc(Nres*sizeof(int *));
  int *nc=(int *)malloc(Nres*sizeof(int));
  *nc1=nc; *clist1=clist; *cnum1=cnum;
  int ia, ja, n;

  for(ia=0; ia<Nres; ia++)nc[ia]=0;
  for(n=0; n<N_int; n++){
    ia=(atoms+Int_list[n].i1)->res;
    ja=(atoms+Int_list[n].i2)->res;
    if((ia>=Nres)||(ja>=Nres))continue;
    nc[ia]++; nc[ja]++;
  }

  int Ncmax=0;
  for(ia=0; ia<Nres; ia++){
    if(nc[ia]>Ncmax)Ncmax=nc[ia];
    clist[ia]=malloc((nc[ia]+1)*sizeof(int));
    cnum[ia]=malloc((nc[ia]+1)*sizeof(int));
    nc[ia]=0;
  }
  printf("Max number contact per residue: %d\n", Ncmax);

  for(n=0; n<N_int; n++){
    ia=(atoms+Int_list[n].i1)->res;
    ja=(atoms+Int_list[n].i2)->res;
    if((ia>=Nres)||(ja>=Nres))continue;
    clist[ia][nc[ia]]=ja; cnum[ia][nc[ia]]=n; nc[ia]++;
    clist[ja][nc[ja]]=ia; cnum[ja][nc[ja]]=n; nc[ja]++;
  }    
  for(ia=0; ia<Nres; ia++){
    clist[ia][nc[ia]]=-1;
    cnum[ia][nc[ia]]=-1;
  }
}


void Write_contact_list(short **cont_i, short **cont_j,
			struct interaction *Int_list, int *N_int)
{
  int k, nc=0, res1=-1, res2=-1;
  for(k=0; k<*N_int; k++){
    if((Int_list[k].res1!=res1)||(Int_list[k].res2!=res2)){
      res1=Int_list[k].res1; res2=Int_list[k].res2; nc++;
    } 
  }

  *cont_i=malloc(nc*sizeof(int));
  *cont_j=malloc(nc*sizeof(int));
  nc=0; res1=-1; res2=-1;
  for(k=0; k<*N_int; k++){
    if((Int_list[k].res1!=res1)||(Int_list[k].res2!=res2)){
      res1=Int_list[k].res1; res2=Int_list[k].res2;
      (*cont_i)[nc]=res1;
      (*cont_j)[nc]=res2;
      nc++;
    } 
  }
  *N_int=nc;
}

int Code_AA(char res){
  int i; for(i=0; i<20; i++)if(res==AA_code[i])return(i);
  if(res=='-')return(20); // gap
  i=(int)res; char r; if(i>96){r=(char)(i-32);}else{r=res;}
  for(i=0; i<20; i++)if(r==AA_code[i])return(i);
  //printf("WARNING, wrong aa type %c\n", res);
  return(-1);
}

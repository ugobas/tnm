#include "atoms.h"
#include <string.h>
int init;


int Atom_type(char *name){
  if(init==0){
    num_BB=8;
    strcpy(BB_name[0],"N ");
    strcpy(BB_name[1],"CA");
    strcpy(BB_name[2],"C ");
    strcpy(BB_name[3],"CB");
    strcpy(BB_name[4],"O ");
    strcpy(BB_name[5],"CG");
    strcpy(BB_name[6],"OG");
    strcpy(BB_name[7],"SG");

    BB_length[0]=1.3301; BB_angle[0]=1.1091; // C-N,  CA-C-N
    BB_length[1]=1.4572; BB_angle[1]=1.0196; // N-CA, C-N-CA
    BB_length[2]=1.5236; BB_angle[2]=1.2018; // CA-C, N-CA-C
    BB_length[3]=1.5319; BB_angle[3]=1.2162; // CA-CB,N-CA-CB
    BB_length[4]=1.2316; BB_angle[4]=1.0359; // C-O,  CA-C-O


    init=1;
  }
  for(int i=0; i<num_BB; i++){
    if((name[0]==BB_name[i][0])&&(name[1]==BB_name[i][1]))
      return(i);
  }
  return(-1);
}


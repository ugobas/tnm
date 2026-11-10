#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom_types.h"

int get_atom_type_nuc(char *atom_name, char *res_name, FILE *file_out)	
{
  char at_name[7];
  if(strncmp(res_name, "  U", 3)==0 || strncmp(res_name, " DU", 3)==0){
    if(strncmp(atom_name, "N1", 2)==0){
      strcpy(at_name, "N3H0u");
    }else if(strncmp(atom_name, "C2", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "O2", 2)==0){
      strcpy(at_name, "O1H0u");
    }else if(strncmp(atom_name, "N3", 2)==0){
      strcpy(at_name, "N3H1s");
    }else if(strncmp(atom_name, "C4", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "O4", 2)==0){
      strcpy(at_name, "O1H0u");
    }else if(strncmp(atom_name, "C5", 2)==0){
      strcpy(at_name, "C3H1s");
    }else if(strncmp(atom_name, "C6", 2)==0){
      strcpy(at_name, "C3H1s");
    }else if(strncmp(atom_name, "C1'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C2'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C3'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C4'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "O2'", 3)==0){
      strcpy(at_name, "O2H1u");
    }else if(strncmp(atom_name, "O3'", 3)==0){
      strcpy(at_name, "O2H0b");
    }else if(strncmp(atom_name, "O4'", 3)==0){
      strcpy(at_name, "O2H0s");
    }else if(strncmp(atom_name, "C5'", 3)==0){
      strcpy(at_name, "C4H2s");
    }else if(strncmp(atom_name, "O5'", 3)==0){
      strcpy(at_name, "O2H0b");
    }else if(strncmp(atom_name, "P  ", 3)==0){
      strcpy(at_name, "P4H0u");
    }else if(strncmp(atom_name, "O1P", 3)==0){
      strcpy(at_name, "O1H0u");
    }else if(strncmp(atom_name, "O2P", 3)==0){
      strcpy(at_name, "O1H0u");
    }else{
      printf("ERROR, not standard atom %s in residue %s\n",atom_name,res_name);
      if(file_out){fprintf(file_out, "%s %s\n", res_name, atom_name);}
      return(ATYPES-1);
    }
  }else if(strncmp(res_name, "  A", 3)==0 || strncmp(res_name, " DA", 3)==0){

   if(strncmp(atom_name, "N1", 2)==0){
      strcpy(at_name, "N3H0u");
    }else if(strncmp(atom_name, "C2", 2)==0){
      strcpy(at_name, "C3H1s");
   }else if(strncmp(atom_name, "N3", 2)==0){
     strcpy(at_name, "N3H0u");
   }else if(strncmp(atom_name, "C4", 2)==0){
     strcpy(at_name, "C3H0s");
   }else if(strncmp(atom_name, "C5", 2)==0){
     strcpy(at_name, "C3H0s");
   }else if(strncmp(atom_name, "C6", 2)==0){
     strcpy(at_name, "C3H0s");
   }else if(strncmp(atom_name, "N6", 2)==0){
     strcpy(at_name, "N3H2u");
   }else if(strncmp(atom_name, "N7", 2)==0){
     strcpy(at_name, "N3H1b");
   }else if(strncmp(atom_name, "C8", 2)==0){
     strcpy(at_name, "C3H1s");
   }else if(strncmp(atom_name, "N9", 2)==0){
     strcpy(at_name, "N3H0u");
   }else if(strncmp(atom_name, "C1'", 3)==0){
     strcpy(at_name, "C4H1s");
   }else if(strncmp(atom_name, "C2'", 3)==0){
     strcpy(at_name, "C4H1ss");
   }else if(strncmp(atom_name, "C3'", 3)==0){
     strcpy(at_name, "C4H1s");
   }else if(strncmp(atom_name, "C4'", 3)==0){
     strcpy(at_name, "C4H1s");
   }else if(strncmp(atom_name, "O2'", 3)==0){
     strcpy(at_name, "O2H1u");
   }else if(strncmp(atom_name, "O3'", 3)==0){
     strcpy(at_name, "O2H0b");
   }else if(strncmp(atom_name, "O4'", 3)==0){
     strcpy(at_name, "O2H0s");
   }else if(strncmp(atom_name, "C5'", 3)==0){
     strcpy(at_name, "C4H2s");
   }else if(strncmp(atom_name, "O5'", 3)==0){
     strcpy(at_name, "O2H0b");
   }else if(strncmp(atom_name, "P  ", 3)==0){
     strcpy(at_name, "P4H0u");
   }else if(strncmp(atom_name, "O1P", 3)==0){
     strcpy(at_name, "O1H0u");
   }else if(strncmp(atom_name, "O2P", 3)==0){
     strcpy(at_name, "O1H0u");
   }else{
     printf("ERROR, not standard atom %s in residue %s\n",atom_name, res_name);
     if(file_out){fprintf(file_out, "%s %s\n", res_name, atom_name);}
     return(ATYPES-1);
   }
   
  }else if(strncmp(res_name, "  G", 3)==0 || strncmp(res_name, " DG", 3)==0){
    
    if(strncmp(atom_name, "N1", 2)==0){
      strcpy(at_name, "N3H1s");
    }else if(strncmp(atom_name, "C2", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "N2", 2)==0){
      strcpy(at_name, "N3H2u");
    }else if(strncmp(atom_name, "N3", 2)==0){
      strcpy(at_name, "N3H1b");
    }else if(strncmp(atom_name, "C4", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "C5", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "C6", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "O6", 2)==0){
      strcpy(at_name, "O1H0u");
    }else if(strncmp(atom_name, "N7", 2)==0){
      strcpy(at_name, "N3H1b");
    }else if(strncmp(atom_name, "C8", 2)==0){
      strcpy(at_name, "C3H1s");
    }else if(strncmp(atom_name, "N9", 2)==0){
      strcpy(at_name, "N3H0u");
    }else if(strncmp(atom_name, "C1'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C2'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C3'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C4'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "O2'", 3)==0){
      strcpy(at_name, "O2H1u");
    }else if(strncmp(atom_name, "O3'", 3)==0){
      strcpy(at_name, "O2H0b");
    }else if(strncmp(atom_name, "O4'", 3)==0){
      strcpy(at_name, "O2H0s");
    }else if(strncmp(atom_name, "C5'", 3)==0){
      strcpy(at_name, "C4H2s");
    }else if(strncmp(atom_name, "O5'", 3)==0){
      strcpy(at_name, "O2H0b");
    }else if(strncmp(atom_name, "P  ", 3)==0){
      strcpy(at_name, "P4H0u");
    }else if(strncmp(atom_name, "O1P", 3)==0){
      strcpy(at_name, "O1H0u");
    }else if(strncmp(atom_name, "O2P", 3)==0){
      strcpy(at_name, "O1H0u");
    }else{
     printf("ERROR, not standard atom %s in residue %s\n",atom_name, res_name);
     if(file_out){fprintf(file_out, "%s %s\n", res_name, atom_name);}
     return(ATYPES-1);
    }
    
  }else if(strncmp(res_name, "  C", 3)==0 || strncmp(res_name, " DC", 3)==0){
    
    if(strncmp(atom_name, "N1", 2)==0){
      strcpy(at_name, "N3H0u");
    }else if(strncmp(atom_name, "C2", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "O2", 2)==0){
      strcpy(at_name, "O1H0u");
    }else if(strncmp(atom_name, "N3", 2)==0){
      strcpy(at_name, "N3H1s");
    }else if(strncmp(atom_name, "C4", 2)==0){
      strcpy(at_name, "C3H0s");
    }else if(strncmp(atom_name, "N4", 2)==0){
      strcpy(at_name, "N3H2u");
    }else if(strncmp(atom_name, "C5", 2)==0){
      strcpy(at_name, "C3H1s");
    }else if(strncmp(atom_name, "C6", 2)==0){
      strcpy(at_name, "C3H1s");
    }else if(strncmp(atom_name, "C1'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C2'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C3'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "C4'", 3)==0){
      strcpy(at_name, "C4H1s");
    }else if(strncmp(atom_name, "O2'", 3)==0){
      strcpy(at_name, "O2H1u");
    }else if(strncmp(atom_name, "O3'", 3)==0){
      strcpy(at_name, "O2H0b");
    }else if(strncmp(atom_name, "O4'", 3)==0){
      strcpy(at_name, "O2H0s");
    }else if(strncmp(atom_name, "C5'", 3)==0){
      strcpy(at_name, "C4H2s");
    }else if(strncmp(atom_name, "O5'", 3)==0){
      strcpy(at_name, "O2H0b");
    }else if(strncmp(atom_name, "P  ", 3)==0){
      strcpy(at_name, "P4H0u");
    }else if(strncmp(atom_name, "O1P", 3)==0){
      strcpy(at_name, "O1H0u");
    }else if(strncmp(atom_name, "O2P", 3)==0){
      strcpy(at_name, "O1H0u");
    }else{
      printf("ERROR, not standard atom %s in residue %s\n",atom_name,res_name);
      if(file_out){fprintf(file_out, "%s %s\n", res_name, atom_name);}
      return(ATYPES-1);
    }
    
  }else{
    printf("ERROR, not standard atom %s in residue %s\n",atom_name,res_name);
    if(file_out){fprintf(file_out, "%s %s\n", res_name, atom_name);}
    return(ATYPES-1);
  }

  for(int i=0; i<ATYPES; i++){
    if(strncmp(atom_names[i], at_name, 3)==0){return(i);}
  }
  printf("ERROR for atom %s in residue %s\n",atom_name,res_name);
  if(file_out){
    fprintf(file_out, "%s %s %s\n", res_name, atom_name, at_name);
  }
  printf("ERROR, atom %s does not exist\n",at_name);
  exit(8);
}

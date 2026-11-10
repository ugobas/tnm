#include "atom_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char atom_names[ATYPES][6]={
  "C3H0s", // Carbonyl carbons with branching (C not gly)
  "C3H0b", // Carboxyl and carbonyl atom w/o branching (s.c. and gly) and aromatic w/o hydrogen
  "C4H1s", // Aliphatic carbons with 1 hydrogen branching from all three heavy atoms (CA)
  "C4H1b", // Aliphatic carbons with 1H and no branching through at least 1 heavy atom bond (Ala, Pro)
  "C3H1s", // Aromatic carbons with 1 hydrogen, small
  "C3H1b", // Aromatic carbons with 1 hydrogen, big
  "C4H2s", // Aliphatic carbons with 2 hydrogens, small (CA Gly)
  "C4H2b", // Aliphatic carbons with 2 hydrogens, big
  "C4H3u", // Aliphatic carbons with 3 hydrogens
  "N3H0u", // Imide hydrogens (PRO)
  "N3H1s", // Amide nitrogens with 1 hydrogen (N mainchain)
  "N3H1b", // Amide nitrogens with 1 hydrogen (side chains)
  "N3H2u", // Amide nitrogens with 2 hydrogens (side chains)
  "N4H3u", // Amide nitrogen charged, 3 hydrogens
  "O1H0u", // All oxygens in carbonyl or carboxyl groups (O, OXT)
  "O2H1u", // All hydroxyl atoms (SER, THR)
  "S2H0u", // Sulfurs with no hydrogens (MET)
  "S2H1u", // Sulfurs with 1 hydrogen (CYS)
  "S=O__", // sulfinyl, sulfonyl, .. any S=0 sulfur atom
  "SEL__", // Selenium
  "P____", // Phosphorous in Phosphate groups. Also Arsenium atoms
  "OP___", // Oxygen atoms from a phosphate group
  "OS___", // sufinyl, sulfonyll or Se=O oxigen atoms
  "X____", // any halogen atom (Cl, F, I)
  "AEth_", // Alkaline earth metals (BE)
  "XXX__"};


int init_pair=0;
int num_pair[ATYPES];

int Ini_atom_types=0;
int ATOM_PAIRS=0;			    
int Pair_index[ATYPES][ATYPES];

void Atom_name(char *name, int n){
  strcpy(name, atom_names[n]);
}

int Get_atom_pair(int a1, int a2)
{
  if(init_pair==0){
    init_pair=1;
    int n=ATYPES, sum=n;
    for(int i=0; i<ATYPES; i++){
      num_pair[i]=sum;
      n--; sum+=n;
    } 
  }
  int amin, amax, p;
  if(a1 < a2){amin=a1; amax=a2;}else{amin=a2; amax=a1;}
  if(amin==0){p=amax;}
  else{p=num_pair[amin-1]+(amax-amin);}
  return(p);
}
	
int find_pair_indx(int indx1, int indx2)
{
  /* This subroutine is analogous to the find_dicam_indx,
     except that, here we input the indices 1 and 2 for
     the atom types 1 and 2 respectively */

   if(Ini_atom_types==0){
     ATOM_PAIRS=ATYPES*ATYPES;
     int k=0;
     for(int i=0; i<ATYPES; i++){
       for(int j=i; j<ATYPES; j++){
	 Pair_index[i][j]=k; k++;
       }
     }
     Ini_atom_types=1;
   }
   if((indx1>=ATYPES)||(indx2>=ATYPES)){
     printf("ERROR, types %d %d not defined (max. is %d)\n",
	    indx1, indx2, ATYPES-1); exit(8);
   }
   if(indx1<=indx2)return(Pair_index[indx1][indx2]);
   return(Pair_index[indx2][indx1]);
}
		    		    
int get_atom_type(char *atom_name, char *res_name, FILE *file_out)	
{
  /* This subroutine assigns atom type according to the classification
     by Tsai et al.  Bioinformatics 2001;17:949-956 */
  
  /**************** Carboxyl carbon **********************/
  if(strncmp(atom_name,"C ",2)==0){
      // Carbonyl carbons with branching (except Gly)
    if(strncmp(res_name,"GLY",3)==0 ||
       strncmp(res_name,"IAS",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"SUI",3)==0 ||
       strncmp(res_name,"KOR",3)==0){
      return(1);
    }else{
      /* Carboxyl and carbonyl atom w/o branching (side chain 
	 and glycine's) and aromatic carbons w/o hydrogen */
      if(strncmp(res_name,"ALA",3)==0 ||
	 strncmp(res_name,"ARG",3)==0 ||
	 strncmp(res_name,"ASN",3)==0 ||
	 strncmp(res_name,"ASP",3)==0 ||
	 strncmp(res_name,"CYS",3)==0 ||
	 strncmp(res_name,"GLN",3)==0 ||
	 strncmp(res_name,"GLU",3)==0 ||
	 strncmp(res_name,"HIS",3)==0 ||
	 strncmp(res_name,"ILE",3)==0 ||
	 strncmp(res_name,"LEU",3)==0 ||
	 strncmp(res_name,"LYS",3)==0 ||
	 strncmp(res_name,"MET",3)==0 ||
	 strncmp(res_name,"PHE",3)==0 ||
	 strncmp(res_name,"PRO",3)==0 ||
	 strncmp(res_name,"SER",3)==0 ||
	 strncmp(res_name,"THR",3)==0 ||
	 strncmp(res_name,"TRP",3)==0 ||
	 strncmp(res_name,"TYR",3)==0 ||
	 strncmp(res_name,"VAL",3)==0 ||
	 // Modified residues:
	 // Ala
	 strncmp(res_name,"DDZ",3)==0 ||
	 // Asn
	 strncmp(res_name,"AHB",3)==0 ||
	 // Asp
	 strncmp(res_name,"BFD",3)==0 ||
	 strncmp(res_name,"CSS",3)==0 ||
	 strncmp(res_name,"SNC",3)==0 ||
	 strncmp(res_name,"KYN",3)==0 ||
	 strncmp(res_name,"PHD",3)==0 ||
	 // Cys
	 strncmp(res_name,"CSO",3)==0 ||
	 strncmp(res_name,"CME",3)==0 ||
	 strncmp(res_name,"OCS",3)==0 ||
	 strncmp(res_name,"CSX",3)==0 ||
	 strncmp(res_name,"CAS",3)==0 ||
	 strncmp(res_name,"CSW",3)==0 ||
	 strncmp(res_name,"YCM",3)==0 ||
	 strncmp(res_name,"CAF",3)==0 ||
	 strncmp(res_name,"CSD",3)==0 ||
	 strncmp(res_name,"CSE",3)==0 ||
	 strncmp(res_name,"BCS",3)==0 ||
	 strncmp(res_name,"CMT",3)==0 ||
	 strncmp(res_name,"QPA",3)==0 ||
	 strncmp(res_name,"CS4",3)==0 ||
	 strncmp(res_name,"CZZ",3)==0 ||
	 // Gln
	 strncmp(res_name,"DGN",3)==0 ||
	 // His
	 strncmp(res_name,"HIC",3)==0 ||
	 strncmp(res_name,"NEP",3)==0 ||
	 // Leu
	 strncmp(res_name,"NLE",3)==0 ||
	 // Lys
	 strncmp(res_name,"MLY",3)==0 ||
	 strncmp(res_name,"KCX",3)==0 ||
	 strncmp(res_name,"MLZ",3)==0 ||
	 strncmp(res_name,"ALY",3)==0 ||
	 strncmp(res_name,"M3L",3)==0 ||
	 strncmp(res_name,"LLP",3)==0 ||
	 // Met
	 strncmp(res_name,"SME",3)==0 ||
	 strncmp(res_name,"MSE",3)==0 ||
	 strncmp(res_name,"MHO",3)==0 ||
	 strncmp(res_name,"MSO",3)==0 ||
	 strncmp(res_name,"FME",3)==0 ||
	 // Phe
	 strncmp(res_name,"DAH",3)==0 ||
	 // Pro
	 strncmp(res_name,"PCA",3)==0 ||
	 strncmp(res_name,"HYP",3)==0 ||
	 // Ser
	 strncmp(res_name,"DSN",3)==0 ||
	 strncmp(res_name,"SAC",3)==0 ||
	 strncmp(res_name,"SEP",3)==0 ||
	 strncmp(res_name,"SEB",3)==0 ||
	 strncmp(res_name,"MIR",3)==0 ||
	 strncmp(res_name,"SGB",3)==0 ||
	 strncmp(res_name,"MIS",3)==0 ||
	 // Thr
	 strncmp(res_name,"TPO",3)==0 ||
	 // Trp
	 strncmp(res_name,"TRO",3)==0 ||
	 strncmp(res_name,"6CW",3)==0 ||
	 // Tyr
	 strncmp(res_name,"TYI",3)==0 ||
	 strncmp(res_name,"IYR",3)==0 ||
	 strncmp(res_name,"PTR",3)==0 ||
	 strncmp(res_name,"OTY",3)==0 ||
	 // Val
	 strncmp(res_name,"ABA",3)==0 ||
	 strncmp(res_name,"AIB",3)==0 ||
	 strncmp(res_name,"PM3",3)==0 ||
	 // ?
	 strncmp(res_name,"UF0",3)==0 ||
	 strncmp(res_name,"CGN",3)==0 ||
	 strncmp(res_name,"SEC",3)==0 ||
	 strncmp(res_name,"CXM",3)==0){
	return(0);
      }   
    }  
    /******************* Alpha carbon ****************************/
  }else if(strncmp(atom_name,"CA ",3)==0){
    if(strncmp(res_name,"GLY",3)==0 ||
       strncmp(res_name,"UF0",3)==0){
      return(6);
      /* Aliphatic carbons with 2 hydrogen, small */
      
    }else if(strncmp(res_name,"ARG",3)==0 ||
	     strncmp(res_name,"ASN",3)==0 ||
	     strncmp(res_name,"ASP",3)==0 ||
	     strncmp(res_name,"CYS",3)==0 ||
	     strncmp(res_name,"GLN",3)==0 ||
	     strncmp(res_name,"GLU",3)==0 ||
	     strncmp(res_name,"HIS",3)==0 ||
	     strncmp(res_name,"ILE",3)==0 ||
	     strncmp(res_name,"LEU",3)==0 ||
	     strncmp(res_name,"LYS",3)==0 ||
	     strncmp(res_name,"MET",3)==0 ||
	     strncmp(res_name,"PHE",3)==0 ||
	     strncmp(res_name,"SER",3)==0 ||
	     strncmp(res_name,"THR",3)==0 ||
	     strncmp(res_name,"TRP",3)==0 ||
	     strncmp(res_name,"TYR",3)==0 ||
	       //
	     strncmp(res_name,"AHB",3)==0 ||
	     strncmp(res_name,"BFD",3)==0 ||
	     strncmp(res_name,"PHD",3)==0 ||
	     strncmp(res_name,"CSS",3)==0 ||
	     strncmp(res_name,"SNC",3)==0 ||
	     strncmp(res_name,"KYN",3)==0 ||
	     strncmp(res_name,"CSO",3)==0 ||
	     strncmp(res_name,"CME",3)==0 ||
	     strncmp(res_name,"OCS",3)==0 ||
	     strncmp(res_name,"CSX",3)==0 ||
	     strncmp(res_name,"CAS",3)==0 ||
	     strncmp(res_name,"CSW",3)==0 ||
	     strncmp(res_name,"YCM",3)==0 ||
	     strncmp(res_name,"CAF",3)==0 ||
	     strncmp(res_name,"CSD",3)==0 ||
	     strncmp(res_name,"CSE",3)==0 ||
	     strncmp(res_name,"BCS",3)==0 ||
	     strncmp(res_name,"CMT",3)==0 ||
	     strncmp(res_name,"KOR",3)==0 ||
	     strncmp(res_name,"QPA",3)==0 ||
	     strncmp(res_name,"CS4",3)==0 ||
	     strncmp(res_name,"CZZ",3)==0 ||
	     strncmp(res_name,"HIC",3)==0 ||
	     strncmp(res_name,"NEP",3)==0 ||
	     strncmp(res_name,"NLE",3)==0 ||
	     strncmp(res_name,"MLY",3)==0 ||
	     strncmp(res_name,"KCX",3)==0 ||
	     strncmp(res_name,"MLZ",3)==0 ||
	     strncmp(res_name,"ALY",3)==0 ||
	     strncmp(res_name,"M3L",3)==0 ||
	     strncmp(res_name,"LLP",3)==0 ||
	     strncmp(res_name,"MSE",3)==0 ||
	     strncmp(res_name,"MHO",3)==0 ||
	     strncmp(res_name,"SME",3)==0 ||
	     strncmp(res_name,"CXM",3)==0 ||
	     strncmp(res_name,"MSO",3)==0 ||
	     strncmp(res_name,"FME",3)==0 ||
	     strncmp(res_name,"DAH",3)==0 ||
	     strncmp(res_name,"SAC",3)==0 ||
	     strncmp(res_name,"SEP",3)==0 ||
	     strncmp(res_name,"SEB",3)==0 ||
	     strncmp(res_name,"MIR",3)==0 ||
	     strncmp(res_name,"SGB",3)==0 ||
	     strncmp(res_name,"MIS",3)==0 ||
	     strncmp(res_name,"TPO",3)==0 ||
	     strncmp(res_name,"TRO",3)==0 ||
	     strncmp(res_name,"6CW",3)==0 ||
	     strncmp(res_name,"TYI",3)==0 ||
	     strncmp(res_name,"IYR",3)==0 ||
	     strncmp(res_name,"PTR",3)==0 ||
	     strncmp(res_name,"OTY",3)==0 ||
	     strncmp(res_name,"VAL",3)==0 ||
	     strncmp(res_name,"IAS",3)==0 ||
	     strncmp(res_name,"ABA",3)==0 ||
	     strncmp(res_name,"AIB",3)==0 ||
	     strncmp(res_name,"PM3",3)==0 ||
	     strncmp(res_name,"SEC",3)==0 ||
	     strncmp(res_name,"DDZ",3)==0){
      return(2);
      /* Aliphatic carbons with 1 hydrogen and branching from 
	 all three heavy atoms */
    }else if(strncmp(res_name,"ALA",3)==0 ||
	     strncmp(res_name,"PRO",3)==0 ||
	     strncmp(res_name,"PCA",3)==0 ||
	     strncmp(res_name,"HYP",3)==0 ||
	     strncmp(res_name,"CGN",3)==0 ){
      return(3);
      /* Aliphatic carbons with 1 hydrogen and no branching through
	 at least 1 heavy atom bond */
    }
    
    /******************* Nitrogen ****************************/
  }else if(strncmp(atom_name,"N ",2)==0){
    
    if(strncmp(res_name,"ALA",3)==0 ||
       strncmp(res_name,"ARG",3)==0 ||
       strncmp(res_name,"ASN",3)==0 ||
       strncmp(res_name,"ASP",3)==0 ||
       strncmp(res_name,"CYS",3)==0 ||
       strncmp(res_name,"GLN",3)==0 ||
       strncmp(res_name,"GLU",3)==0 ||
       strncmp(res_name,"GLY",3)==0 ||
       strncmp(res_name,"HIS",3)==0 ||
       strncmp(res_name,"ILE",3)==0 ||
       strncmp(res_name,"LEU",3)==0 ||
       strncmp(res_name,"LYS",3)==0 ||
       strncmp(res_name,"MET",3)==0 ||
       strncmp(res_name,"PHE",3)==0 ||
       strncmp(res_name,"SER",3)==0 ||
       strncmp(res_name,"THR",3)==0 ||
       strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TYR",3)==0 ||
       //
       strncmp(res_name,"DDZ",3)==0 ||
       strncmp(res_name,"AHB",3)==0 ||
       strncmp(res_name,"BFD",3)==0 ||
       strncmp(res_name,"PHD",3)==0 ||
       strncmp(res_name,"CSS",3)==0 ||
       strncmp(res_name,"SNC",3)==0 ||
       strncmp(res_name,"KYN",3)==0 ||
       strncmp(res_name,"CSO",3)==0 ||
       strncmp(res_name,"CME",3)==0 ||
       strncmp(res_name,"OCS",3)==0 ||
       strncmp(res_name,"CSX",3)==0 ||
       strncmp(res_name,"CAS",3)==0 ||
       strncmp(res_name,"CSW",3)==0 ||
       strncmp(res_name,"YCM",3)==0 ||
       strncmp(res_name,"CSD",3)==0 ||
       strncmp(res_name,"CSE",3)==0 ||
       strncmp(res_name,"BCS",3)==0 ||
       strncmp(res_name,"CMT",3)==0 ||
       strncmp(res_name,"QPA",3)==0 ||
       strncmp(res_name,"CS4",3)==0 ||
       strncmp(res_name,"CZZ",3)==0 ||
       strncmp(res_name,"KOR",3)==0 ||
       strncmp(res_name,"HIC",3)==0 ||
       strncmp(res_name,"CAF",3)==0 ||
       strncmp(res_name,"NEP",3)==0 ||
       strncmp(res_name,"NLE",3)==0 ||
       strncmp(res_name,"MLY",3)==0 ||
       strncmp(res_name,"KCX",3)==0 ||
       strncmp(res_name,"MLZ",3)==0 ||
       strncmp(res_name,"ALY",3)==0 ||
       strncmp(res_name,"M3L",3)==0 ||
       strncmp(res_name,"LLP",3)==0 ||
       strncmp(res_name,"MSE",3)==0 ||
       strncmp(res_name,"SME",3)==0 ||
       strncmp(res_name,"MHO",3)==0 ||
       strncmp(res_name,"CXM",3)==0 ||
       strncmp(res_name,"MSO",3)==0 ||
       strncmp(res_name,"FME",3)==0 ||
       strncmp(res_name,"DAH",3)==0 ||
       strncmp(res_name,"SAC",3)==0 ||
       strncmp(res_name,"SEP",3)==0 ||
       strncmp(res_name,"SEB",3)==0 ||
       strncmp(res_name,"MIR",3)==0 ||
       strncmp(res_name,"SGB",3)==0 ||
       strncmp(res_name,"MIS",3)==0 ||
       strncmp(res_name,"TPO",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"PTR",3)==0 ||
       strncmp(res_name,"OTY",3)==0 ||
       strncmp(res_name,"VAL",3)==0 ||
       strncmp(res_name,"IAS",3)==0 ||
       strncmp(res_name,"ABA",3)==0 ||
       strncmp(res_name,"SUI",3)==0 ||
       strncmp(res_name,"IYR",3)==0 ||
       strncmp(res_name,"AIB",3)==0 ||
       strncmp(res_name,"SEC",3)==0 ||
       strncmp(res_name,"PM3",3)==0){
      return(10);
      /* Amide nitrogens with 1 hydrogen (all other mainchain N's) */

    }else if(strncmp(res_name,"PRO",3)==0 ||
	     strncmp(res_name,"PCA",3)==0 ||
	     strncmp(res_name,"SNN",3)==0 ||
	     strncmp(res_name,"CSY",3)==0 ||
	     strncmp(res_name,"HYP",3)==0 ||
	     strncmp(res_name,"CGN",3)==0 ){
      /* Imide nitrogens */
      return(9);
    }else if(strncmp(res_name,"UF0",3)==0){
      /* Amide nitrogens with 1 hydrogen (on side chains) */
      return(11);
    }
    /******************** Carboxyl Oxygen **********************/
  }else if(strncmp(atom_name,"O ",2)==0 ||
	   strncmp(atom_name,"OXT",3)==0){
    if(strncmp(res_name,"ALA",3)==0 ||
       strncmp(res_name,"ARG",3)==0  ||
       strncmp(res_name,"ASN",3)==0  ||
       strncmp(res_name,"GLN",3)==0  ||
       strncmp(res_name,"GLU",3)==0  ||
       strncmp(res_name,"GLY",3)==0  ||
       strncmp(res_name,"HIS",3)==0  ||
       strncmp(res_name,"LYS",3)==0  ||
       strncmp(res_name,"ILE",3)==0  ||
       strncmp(res_name,"LEU",3)==0  ||
       strncmp(res_name,"MET",3)==0  ||
       strncmp(res_name,"PHE",3)==0  ||
       strncmp(res_name,"PRO",3)==0  ||
       strncmp(res_name,"SER",3)==0  ||
       strncmp(res_name,"THR",3)==0  ||
       strncmp(res_name,"TRP",3)==0  ||
       strncmp(res_name,"TYR",3)==0  ||
       //
       strncmp(res_name,"DDZ",3)==0  ||
       strncmp(res_name,"AHB",3)==0  ||
       strncmp(res_name,"ASP",3)==0  ||
       strncmp(res_name,"BFD",3)==0  ||
       strncmp(res_name,"CSS",3)==0  ||
       strncmp(res_name,"CYS",3)==0  ||
       strncmp(res_name,"CSO",3)==0  ||
       strncmp(res_name,"CME",3)==0  ||
       strncmp(res_name,"OCS",3)==0  ||
       strncmp(res_name,"CSX",3)==0  ||
       strncmp(res_name,"CAS",3)==0  ||
       strncmp(res_name,"CSW",3)==0  ||
       strncmp(res_name,"YCM",3)==0  ||
       strncmp(res_name,"CAF",3)==0  ||
       strncmp(res_name,"CSD",3)==0  ||
       strncmp(res_name,"CSE",3)==0  ||
       strncmp(res_name,"BCS",3)==0  ||
       strncmp(res_name,"CMT",3)==0  ||
       strncmp(res_name,"CZZ",3)==0  ||
       strncmp(res_name,"HIC",3)==0  ||
       strncmp(res_name,"NEP",3)==0  ||
       strncmp(res_name,"NLE",3)==0  ||
       strncmp(res_name,"MLY",3)==0  ||
       strncmp(res_name,"KCX",3)==0  ||
       strncmp(res_name,"MLZ",3)==0  ||
       strncmp(res_name,"ALY",3)==0  ||
       strncmp(res_name,"M3L",3)==0  ||
       strncmp(res_name,"LLP",3)==0  ||
       strncmp(res_name,"MSE",3)==0  ||
       strncmp(res_name,"SME",3)==0  ||
       strncmp(res_name,"MHO",3)==0 ||
       strncmp(res_name,"CXM",3)==0  ||
       strncmp(res_name,"MSO",3)==0  ||
       strncmp(res_name,"DAH",3)==0  ||
       strncmp(res_name,"PCA",3)==0  ||
       strncmp(res_name,"HYP",3)==0  ||
       strncmp(res_name,"SAC",3)==0  ||
       strncmp(res_name,"SEP",3)==0  ||
       strncmp(res_name,"SEB",3)==0  ||
       strncmp(res_name,"MIR",3)==0  ||
       strncmp(res_name,"SGB",3)==0  ||
       strncmp(res_name,"TPO",3)==0  ||
       strncmp(res_name,"MIS",3)==0  ||
       strncmp(res_name,"TRO",3)==0  ||
       strncmp(res_name,"6CW",3)==0  ||
       strncmp(res_name,"TYI",3)==0  ||
       strncmp(res_name,"PTR",3)==0  ||
       strncmp(res_name,"VAL",3)==0  ||
       strncmp(res_name,"IAS",3)==0  ||
       strncmp(res_name,"ABA",3)==0  ||
       strncmp(res_name,"IYR",3)==0  ||
       strncmp(res_name,"AIB",3)==0  ||
       strncmp(res_name,"SNC",3)==0  ||
       strncmp(res_name,"EYG",3)==0  ||
       strncmp(res_name,"KYN",3)==0  ||
       strncmp(res_name,"KOR",3)==0  ||
       strncmp(res_name,"AOP",3)==0  ||
       strncmp(res_name,"UF0",3)==0  ||
       strncmp(res_name,"QPA",3)==0  ||
       strncmp(res_name,"CS4",3)==0  ||
       strncmp(res_name,"PHD",3)==0  ||
       strncmp(res_name,"OTY",3)==0  ||
       strncmp(res_name,"C12",3)==0  ||
       strncmp(res_name,"SEC",3)==0 ||
       strncmp(res_name,"CGN",3)==0){
      return(14);	
      /* All oxygens in carbonyl or carboxyl groups (no distinction
	 made between oxygens in carbonyl in carboxyl group) */
    }
    
    
    /****************** Beta carbon ********************/
  }else if(strncmp(atom_name,"CB ",3)==0){
    if(strncmp(res_name,"ILE",3)==0 ||
       strncmp(res_name,"THR",3)==0 ||
       strncmp(res_name,"VAL",3)==0 ||
       strncmp(res_name,"TPO",3)==0 ||
       strncmp(res_name,"AHB",3)==0 ||
       strncmp(res_name,"DDZ",3)==0){
      
      /* Aliphatic carbons with 1 hydrogen and no branching through
	 at least 1 heavy atom bond */
      return(3);
      
    }else if(strncmp(res_name,"ARG",3)==0 ||
	     strncmp(res_name,"ASN",3)==0 ||
	     strncmp(res_name,"ASP",3)==0 ||
	     strncmp(res_name,"GLN",3)==0 ||
	     strncmp(res_name,"GLU",3)==0 ||
	     strncmp(res_name,"HIS",3)==0 ||
	     strncmp(res_name,"LEU",3)==0 ||
	     strncmp(res_name,"LYS",3)==0 ||
	     strncmp(res_name,"MET",3)==0 ||
	     strncmp(res_name,"MSE",3)==0 ||
	     strncmp(res_name,"PHE",3)==0 ||
	     strncmp(res_name,"SER",3)==0 ||
	     strncmp(res_name,"TRP",3)==0 ||
	     strncmp(res_name,"TYR",3)==0 ||
	     //
	     strncmp(res_name,"BFD",3)==0 ||
	     strncmp(res_name,"KYN",3)==0 ||
	     strncmp(res_name,"PHD",3)==0 ||
	     strncmp(res_name,"HIC",3)==0 ||
	     strncmp(res_name,"NEP",3)==0 ||
	     strncmp(res_name,"NLE",3)==0 ||
	     strncmp(res_name,"MLY",3)==0 ||
	     strncmp(res_name,"KCX",3)==0 ||
	     strncmp(res_name,"MLZ",3)==0 ||
	     strncmp(res_name,"ALY",3)==0 ||
	     strncmp(res_name,"M3L",3)==0 ||
	     strncmp(res_name,"LLP",3)==0 ||
	     strncmp(res_name,"SME",3)==0 ||
	     strncmp(res_name,"CXM",3)==0 ||
	     strncmp(res_name,"MHO",3)==0 ||
	     strncmp(res_name,"MSO",3)==0 ||
	     strncmp(res_name,"FME",3)==0 ||
	     strncmp(res_name,"DAH",3)==0 ||
	     strncmp(res_name,"CSY",3)==0 ||
	     strncmp(res_name,"SAC",3)==0 ||
	     strncmp(res_name,"SEP",3)==0 ||
	     strncmp(res_name,"SEB",3)==0 ||
	     strncmp(res_name,"MIR",3)==0 ||
	     strncmp(res_name,"SGB",3)==0 ||
	     strncmp(res_name,"MIS",3)==0 ||
	     strncmp(res_name,"TRO",3)==0 ||
	     strncmp(res_name,"6CW",3)==0 ||
	     strncmp(res_name,"TYI",3)==0 ||
	     strncmp(res_name,"IAS",3)==0 ||
	     strncmp(res_name,"ABA",3)==0 ||
	     strncmp(res_name,"IYR",3)==0 ||
	     strncmp(res_name,"UF0",3)==0 ||
	     strncmp(res_name,"PTR",3)==0 ||
	     strncmp(res_name,"OTY",3)==0 ||
	     strncmp(res_name,"SEC",3)==0 ||
	     strncmp(res_name,"PM3",3)==0){
      
      /* Aliphatic carbons with 2 hydrogen, small */
      return(6);
      
    }else if(strncmp(res_name,"CYS",3)==0 ||
	     strncmp(res_name,"PRO",3)==0 ||
	     strncmp(res_name,"CSS",3)==0 ||
	     strncmp(res_name,"SNC",3)==0 ||
	     strncmp(res_name,"CSO",3)==0 ||
	     strncmp(res_name,"CME",3)==0 ||
	     strncmp(res_name,"OCS",3)==0 ||
	     strncmp(res_name,"CSX",3)==0 ||
	     strncmp(res_name,"CAS",3)==0 ||
	     strncmp(res_name,"CSW",3)==0 ||
	     strncmp(res_name,"YCM",3)==0 ||
	     strncmp(res_name,"CAF",3)==0 ||
	     strncmp(res_name,"CSD",3)==0 ||
	     strncmp(res_name,"CSE",3)==0 ||
	     strncmp(res_name,"BCS",3)==0 ||
	     strncmp(res_name,"CMT",3)==0 ||
	     strncmp(res_name,"CZZ",3)==0 ||
	     strncmp(res_name,"PCA",3)==0 ||
	     strncmp(res_name,"KOR",3)==0 ||
	     strncmp(res_name,"QPA",3)==0 ||
	     strncmp(res_name,"CS4",3)==0 ||
	     strncmp(res_name,"HYP",3)==0 ||
	     strncmp(res_name,"CGN",3)==0 ){
      
      /* Aliphatic carbons with 2 hydrogens, big */ 	    
      return(7);
      
    }else if(strncmp(res_name,"ALA",3)==0){
      
      /* Aliphatic carbons with three hydrogens, i.e. methyl groups */
      return(8);
      
    }else if(strncmp(res_name,"ASN",3)==0 ||
	     strncmp(res_name,"ASP",3)==0 ||
	     strncmp(res_name,"HIS",3)==0 ||
	     strncmp(res_name,"PHE",3)==0 ||
	     strncmp(res_name,"TRP",3)==0 ||
	     strncmp(res_name,"TYR",3)==0 ||
	     //
	     strncmp(res_name,"AHB",3)==0 ||
	     strncmp(res_name,"KYN",3)==0 ||
	     strncmp(res_name,"BFD",3)==0 ||
	     strncmp(res_name,"PHD",3)==0 ||
	     strncmp(res_name,"CSY",3)==0 ||
	     strncmp(res_name,"HIC",3)==0 ||
	     strncmp(res_name,"NEP",3)==0 ||
	     strncmp(res_name,"6CW",3)==0 ||
	     strncmp(res_name,"TYI",3)==0 ||
	     strncmp(res_name,"SUI",3)==0 ||
	     strncmp(res_name,"PTR",3)==0 ||
	     strncmp(res_name,"OTY",3)==0 ||
	     strncmp(res_name,"DAH",3)==0 ||
	     strncmp(res_name,"PM3",3)==0){
      
      /* Carboxyl and carbonyl atom w/o branching (side chain 
	 and glycine's) and aromatic carbons w/o hydrogen */
      
      return(1);
    }
    
    
    /********************** Gamma carbon *************************/
  }else if(strncmp(atom_name,"CG ",3)==0){
    if(strncmp(res_name,"IAS",3)==0){
      /*Carbonyl carbons with branching (mainchain carbonyls 
	from residues with a Cp so no gly carbon) */       
      return(0);
      
    }else if(strncmp(res_name,"ASN",3)==0 ||
	     strncmp(res_name,"ASP",3)==0 ||
	     strncmp(res_name,"HIS",3)==0 ||
	     strncmp(res_name,"PHE",3)==0 ||
	     strncmp(res_name,"TYR",3)==0 ||
	     strncmp(res_name,"TRP",3)==0 ||
	     strncmp(res_name,"AHB",3)==0 ||
	     strncmp(res_name,"KYN",3)==0 ||
	     strncmp(res_name,"BFD",3)==0 ||
	     strncmp(res_name,"PHD",3)==0 ||
	     strncmp(res_name,"CSY",3)==0 ||
	     strncmp(res_name,"HIC",3)==0 ||
	     strncmp(res_name,"NEP",3)==0 ||
	     strncmp(res_name,"6CW",3)==0 ||
	     strncmp(res_name,"TYI",3)==0 ||
	     strncmp(res_name,"SUI",3)==0 ||
	     strncmp(res_name,"PTR",3)==0 ||
	     strncmp(res_name,"TRO",3)==0 ||
	     strncmp(res_name,"OTY",3)==0 ||
	     strncmp(res_name,"DAH",3)==0 ||
	     strncmp(res_name,"PM3",3)==0){
      
      /* Carboxyl and carbonyl atom w/o branching (side chain 
	 and glycine's) and aromatic carbons w/o hydrogen */	    
      return(1);
    }else if(strncmp(res_name,"LEU",3)==0 ||
	     strncmp(res_name,"HYP",3)==0){
      
      /* Aliphatic carbons with 1 hydrogen and no branching through
	 at least 1 heavy atom bond */
      return(3);
      
    }else if(strncmp(res_name,"IYR",3)==0){
      /* Small aromatic carbons with 1 hydrogen */	    
      return(4);
    }else if(strncmp(res_name,"ARG",3)==0 ||
	     strncmp(res_name,"GLN",3)==0 ||
	     strncmp(res_name,"GLU",3)==0 ||
	     strncmp(res_name,"LYS",3)==0 ||
	     strncmp(res_name,"NLE",3)==0 ||
	     strncmp(res_name,"MLY",3)==0 ||
	     strncmp(res_name,"KCX",3)==0 ||
	     strncmp(res_name,"MLZ",3)==0 ||
	     strncmp(res_name,"ALY",3)==0 ||
	     strncmp(res_name,"M3L",3)==0 ||
	     strncmp(res_name,"LLP",3)==0 ){
      /* Aliphatic carbons with 2 hydrogen, small */
      return(6);
    }else if(strncmp(res_name,"MET",3)==0 ||
	     strncmp(res_name,"PRO",3)==0 ||
	     strncmp(res_name,"MSE",3)==0 ||
	     strncmp(res_name,"SME",3)==0 ||
	     strncmp(res_name,"CXM",3)==0 ||
	     strncmp(res_name,"MHO",3)==0 ||
	     strncmp(res_name,"MSO",3)==0 ||
	     strncmp(res_name,"FME",3)==0 ||
	     strncmp(res_name,"PCA",3)==0 ||
	     strncmp(res_name,"EYG",3)==0 ||
	     strncmp(res_name,"CGN",3)==0){
      
      /* Aliphatic carbons with 2 hydrogens, big */ 
      return(7);
    }else if(strncmp(res_name,"ABA",3)==0){
      /* Aliphatic carbons with three hydrogens, i.e. methyl groups */
      return(8);
    }
    
  }else if(strncmp(atom_name,"CG1",3)==0){
    if(strncmp(res_name,"ILE",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0) return(7);
    if(strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"VAL",3)==0 ||
       strncmp(res_name,"CRO",3)==0) return(8);

  }else if(strncmp(atom_name,"CG2",3)==0){
    if(strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0) return(1);
    if(strncmp(res_name,"ILE",3)==0 ||
       strncmp(res_name,"THR",3)==0 ||
       strncmp(res_name,"TPO",3)==0 ||
       strncmp(res_name,"VAL",3)==0) return(8);

    /************************** Delta carbon **************************/
  }else if(strncmp(atom_name,"CD ",3)==0){
    if(strncmp(res_name,"GLN",3)==0 ||
       strncmp(res_name,"GLU",3)==0 ||
       strncmp(res_name,"PCA",3)==0 ) return(1);
    if(strncmp(res_name,"ARG",3)==0 ||
       strncmp(res_name,"LYS",3)==0 ||
       strncmp(res_name,"PRO",3)==0 ||
       strncmp(res_name,"NLE",3)==0 ||
       strncmp(res_name,"MLY",3)==0 ||
       strncmp(res_name,"KCX",3)==0 ||
       strncmp(res_name,"MLZ",3)==0 ||
       strncmp(res_name,"ALY",3)==0 ||
       strncmp(res_name,"M3L",3)==0 ||
       strncmp(res_name,"LLP",3)==0 ||
       strncmp(res_name,"HYP",3)==0 ||
       strncmp(res_name,"YCM",3)==0 ) return(6);
    if(strncmp(res_name,"CGN",3)==0 ) return(0);
    if(strncmp(res_name,"IYR",3)==0 ) return(4);
    if(strncmp(res_name,"BCS",3)==0 ) return(7);
    
  }else if(strncmp(atom_name,"CD1",3)==0){

    if(strncmp(res_name,"PHE",3)==0 ||
       strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TYR",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"PTR",3)==0 ||
       strncmp(res_name,"DAH",3)==0 ||
       strncmp(res_name,"PM3",3)==0 )return(4);
    if(strncmp(res_name,"ILE",3)==0 ||
       strncmp(res_name,"LEU",3)==0) return(8);
    if(strncmp(res_name,"KYN",3)==0)return(5);
    if(strncmp(res_name,"OTY",3)==0)return(1);
    
  }else if(strncmp(atom_name,"CD2",3)==0){
    if(strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0 ||
       strncmp(res_name,"KYN",3)==0 ) return(1);
    if(strncmp(res_name,"HIS",3)==0 ||
       strncmp(res_name,"TYR",3)==0 ||
       strncmp(res_name,"HIC",3)==0 ||
       strncmp(res_name,"NEP",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"PTR",3)==0 ||
       strncmp(res_name,"OTY",3)==0 ||
       strncmp(res_name,"DAH",3)==0 ||
       strncmp(res_name,"PM3",3)==0 )return(4);
    if(strncmp(res_name,"PHE",3)==0 )return(5);
    if(strncmp(res_name,"LEU",3)==0 )return(8);

    /************* epsilon carbon ************/
  }else if(strncmp(atom_name,"CE ",3)==0){
    if(strncmp(res_name,"LYS",3)==0 ||
       strncmp(res_name,"CME",3)==0 ||
       strncmp(res_name,"MLY",3)==0 ||
       strncmp(res_name,"KCX",3)==0 ||
       strncmp(res_name,"MLZ",3)==0 ||
       strncmp(res_name,"ALY",3)==0 ||
       strncmp(res_name,"M3L",3)==0 ||
       strncmp(res_name,"LLP",3)==0 ) return(7);
    if(strncmp(res_name,"MET",3)==0 ||
       strncmp(res_name,"NLE",3)==0 ||	     
       strncmp(res_name,"MSE",3)==0 ||
       strncmp(res_name,"SME",3)==0 ||
       strncmp(res_name,"MLY",3)==0 ||
       strncmp(res_name,"CXM",3)==0 ||
       strncmp(res_name,"MHO",3)==0 ||
       strncmp(res_name,"MSO",3)==0 ||
       strncmp(res_name,"FME",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ) return(8);
    if(strncmp(res_name,"YCM",3)==0 ||
       strncmp(res_name,"BCS",3)==0 ||
       strncmp(res_name,"IYR",3)==0 ) return(1);
    if(strncmp(res_name,"SEB",3)==0 ) return(6);
    
  }else if(strncmp(atom_name,"CE1",3)==0){
    if(strncmp(res_name,"HIS",3)==0 ||
       strncmp(res_name,"TYR",3)==0 ||
       strncmp(res_name,"HIC",3)==0 ||
       strncmp(res_name,"NEP",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"PTR",3)==0 ||
       strncmp(res_name,"OTY",3)==0 ||
       strncmp(res_name,"DAH",3)==0 ||
       strncmp(res_name,"PM3",3)==0 ) return(4);
    if(strncmp(res_name,"PHE",3)==0 ||
       strncmp(res_name,"KYN",3)==0 ) return(5);
    if(strncmp(res_name,"CAS",3)==0 ||
       strncmp(res_name,"CAF",3)==0 ) return(8);

  }else if(strncmp(atom_name,"CE2",3)==0){
    if(strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0 ||
       strncmp(res_name,"DAH",3)==0 )return(1);
    if(strncmp(res_name,"TYR",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"PTR",3)==0 ||
       strncmp(res_name,"OTY",3)==0 ||
       strncmp(res_name,"PM3",3)==0 ) return(4);
    if(strncmp(res_name,"PHE",3)==0 ||
       strncmp(res_name,"KYN",3)==0 ) return(5);
    if(strncmp(res_name,"CAS",3)==0 ||
       strncmp(res_name,"CAF",3)==0 ) return(8);
    
  }else if(strncmp(atom_name,"CE3",3)==0){
    if(strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0) return(5);
    
    /************* zeta carbon ************/
  }else if(strncmp(atom_name,"CZ ",3)==0){
    if(strncmp(res_name,"ARG",3)==0 ||
       strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"TYR",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"SEB",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"PTR",3)==0 ||
       strncmp(res_name,"OTY",3)==0 ||
       strncmp(res_name,"DAH",3)==0 ||
       strncmp(res_name,"PM3",3)==0 ) return(1);
    if(strncmp(res_name,"PHE",3)==0 ||
       strncmp(res_name,"KYN",3)==0 ) return(5);
    if(strncmp(res_name,"CME",3)==0 ) return(6);
    if(strncmp(res_name,"HIC",3)==0 ) return(8);

  }else if(strncmp(atom_name,"CZ2",3)==0){
    if(strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0 ||
       strncmp(res_name,"BCS",3)==0 ) return(5); 
    
  }else if(strncmp(atom_name,"CZ3",3)==0){
    if(strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0 ) return(5);
    
    /************* eta carbon ************/
  }else if(strncmp(atom_name,"CH ",3)==0){
    if(strncmp(res_name,"ALY",3)==0 )return(1);
    if(strncmp(res_name,"IYR",3)==0 )return(4);
    if(strncmp(res_name,"BCS",3)==0 )return(5);
    
  }else if(strncmp(atom_name,"CH1",3)==0){
    if(strncmp(res_name,"SEB",3)==0 ) return(4);
    if(strncmp(res_name,"MLY",3)==0 ) return(8);
    
  }else if(strncmp(atom_name,"CH2",3)==0){
    if(strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0 ||
       strncmp(res_name,"SEB",3)==0 ) return(5);
    if(strncmp(res_name,"MLY",3)==0 ) return(8);
    
    /*********** Side chain Oxygen ****************/
  }else if(strncmp(atom_name,"OG1",3)==0){
    if(strncmp(res_name,"THR",3)==0 ||
       strncmp(res_name,"TPO",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"DDZ",3)==0 ) return(15);
    
  }else if(strncmp(atom_name,"OG2",3)==0){
    if(strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"DDZ",3)==0) return(15);
    
  }else if(strncmp(atom_name,"OG",2)==0){
    if(strncmp(res_name,"SER",3)==0 ||
       strncmp(res_name,"SEB",3)==0 || 
       strncmp(res_name,"SAC",3)==0 ||
       strncmp(res_name,"SEP",3)==0 ||
       strncmp(res_name,"MIR",3)==0 ) return(15);
    if(strncmp(res_name,"SGB",3)==0 ||
       strncmp(res_name,"MIS",3)==0 ) return(21);
    
  }else if(strncmp(atom_name,"OD1",3)==0){
    if(strncmp(res_name,"ASN",3)==0 ||
       strncmp(res_name,"ASP",3)==0 ||
       strncmp(res_name,"AHB",3)==0 ||
       strncmp(res_name,"IAS",3)==0 ||
       strncmp(res_name,"AHB",3)==0 ||
       strncmp(res_name,"BFD",3)==0) return(14);
    if(strncmp(res_name,"CYS",3)==0 ||
       strncmp(res_name,"MET",3)==0 ||
       strncmp(res_name,"MHO",3)==0 ||
       strncmp(res_name,"OCS",3)==0 ||
       strncmp(res_name,"CSD",3)==0 )return(22);
    if(strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"HYP",3)==0)return(15);
    if(strncmp(res_name,"PHD",3)==0)return(21);
    
  }else if(strncmp(atom_name,"OD2",3)==0){
    if(strncmp(res_name,"ASP",3)==0 ||
       strncmp(res_name,"BFD",3)==0 ||
       strncmp(res_name,"OCS",3)== 0)return(14);
    if(strncmp(res_name,"CSD",3)==0)return(22);
    
  }else if(strncmp(atom_name,"OD3",3)==0){
    if(strncmp(res_name,"OCS",3)== 0)return(22);
    
  }else if(strncmp(atom_name,"OD",2)==0){
    if(strncmp(res_name,"CYS",3)==0 ||
       strncmp(res_name,"CSX",3)==0 )return(22);
    
  }else if(strncmp(atom_name,"OE1",3)==0){
    if(strncmp(res_name,"GLN",3)==0 ||
       strncmp(res_name,"GLU",3)==0 ||
       strncmp(res_name,"CGN",3)==0) return(14);
    if(strncmp(res_name,"K1R",3)==0) return(22);
    
  }else if(strncmp(atom_name,"OE2",3)==0){
    if(strncmp(res_name,"GLU",3)==0)return(14);
    if(strncmp(res_name,"DAH",3)==0)return(15);
    
  }else if(strncmp(atom_name,"OE ",3)==0){
    if(strncmp(res_name,"PCA",3)==0 ||
       strncmp(res_name,"SNC",3)==0)return(14);
    if(strncmp(res_name,"SME",3)==0 ||
       strncmp(res_name,"MSO",3)==0)return(22);
    
  }else if(strncmp(atom_name,"OF ",3)==0){
    if(strncmp(res_name,"IYR",3)==0)return(15);
    
  }else if(strncmp(atom_name,"OZ ",3)==0){
    if(strncmp(res_name,"DAH",3)==0)return(15);
    
  }else if(strncmp(atom_name,"OZ1",3)==0){
    if(strncmp(res_name,"YCM",3)==0)return(14);
    
  }else if(strncmp(atom_name,"OH ",3)==0){
    if(strncmp(res_name,"TYR",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"TYI",3)==0 ||
       strncmp(res_name,"CME",3)==0 ||
       strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"OTY",3)==0 )return(15);
    if(strncmp(res_name,"ALY",3)==0 )return(14);
    if(strncmp(res_name,"PTR",3)==0)return(21);
    
  }else if(strncmp(atom_name,"O1 ",3)==0){
    if(strncmp(res_name,"CYS",3)==0)return(22);
    
  }else if(strncmp(atom_name,"O2 ",3)==0){
    if(strncmp(res_name,"SNN",3)==0 ||
       strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"K1R",3)==0) return(14);
    if(strncmp(res_name,"C12",3)==0)return(15);
    if(strncmp(res_name,"SGB",3)==0)return(21);
    
    /************** Nitrogen ********************/
  }else if(strncmp(atom_name,"ND2",3)==0){
    if(strncmp(res_name,"ASN",3)==0 ||
       strncmp(res_name,"AHB",3)==0 )return(12);
  }else if(strncmp(atom_name,"ND1",3)==0){
    if(strncmp(res_name,"HIS",3)==0 ||
       strncmp(res_name,"HIC",3)==0 ||
       strncmp(res_name,"NEP",3)==0)return(11);
    
  }else if(strncmp(atom_name,"NE2",3)==0){
    if(strncmp(res_name,"HIS",3)==0 ||
       strncmp(res_name,"HIC",3)==0 ||
       strncmp(res_name,"NEP",3)==0)return(11);
    if(strncmp(res_name,"GLN",3)==0)return(12);
    
  }else if(strncmp(atom_name,"NE ",3)==0){
    if(strncmp(res_name,"ARG",3)==0)return(11);
    
  }else if(strncmp(atom_name,"NE1",3)==0){
    if(strncmp(res_name,"TRP",3)==0 ||
       strncmp(res_name,"TRO",3)==0 ||
       strncmp(res_name,"6CW",3)==0)return(11);
    
  }else if(strncmp(atom_name,"NZ ",3)==0){
    if(strncmp(res_name,"LYS",3)==0 ||
       strncmp(res_name,"MLY",3)==0 ||
       strncmp(res_name,"KCX",3)==0 ||
       strncmp(res_name,"MLZ",3)==0 ||
       strncmp(res_name,"M3L",3)==0 ||
       strncmp(res_name,"LLP",3)==0 )return(13);
    if(strncmp(res_name,"ALY",3)==0 )return(11);
    
  }else if(strncmp(atom_name,"NZ2",3)==0){
    if(strncmp(res_name,"YCM",3)==0)return(12);
    
  }else if(strncmp(atom_name,"NH",2)==0){
    if(strncmp(res_name,"ARG",3)==0 )return(12);
    
  }else if(strncmp(atom_name,"N1",2)==0){
    if(strncmp(res_name,"SNN",3)==0)return(9);
    if(strncmp(res_name,"CSY",3)==0 ||
       strncmp(res_name,"C12",3)==0 ||
       strncmp(res_name,"CR2",3)==0 ||
       strncmp(res_name,"CRO",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0)return(10);
    if(strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"LLP",3)==0)return(11);
    if(strncmp(res_name,"KYN",3)==0 ||
       strncmp(res_name,"K1R",3)==0 )return(13);
    
    /************** Sulfur ****************/
  }else if(strncmp(atom_name,"SD ",3)==0){
    if(strncmp(res_name,"MET",3)==0 ||
       strncmp(res_name,"MSE",3)==0 ||
       strncmp(res_name,"MHO",3)==0 ||
       strncmp(res_name,"CSS",3)==0 ||
       strncmp(res_name,"CME",3)==0 ||
       strncmp(res_name,"CXM",3)==0 ||
       strncmp(res_name,"NRQ",3)==0 ||
       strncmp(res_name,"CH6",3)==0 ||
       strncmp(res_name,"EYG",3)==0 ||
       strncmp(res_name,"FME",3)==0 )return(16);
    if(strncmp(res_name,"SEB",3)==0 )return(18);
    
  }else if(strncmp(atom_name,"SG ",3)==0){
    if(strncmp(res_name,"CSS",3)==0 ||
       strncmp(res_name,"SNC",3)==0 ||
       strncmp(res_name,"CME",3)==0 ||
       strncmp(res_name,"CSO",3)==0 ||
       strncmp(res_name,"CAS",3)==0 ||
       strncmp(res_name,"YCM",3)==0 ||
       strncmp(res_name,"CAF",3)==0 ||  
       strncmp(res_name,"K1R",3)==0 ||
       strncmp(res_name,"KOR",3)==0 ||
       strncmp(res_name,"QPA",3)==0 ||
       strncmp(res_name,"CS4",3)==0 ||
       strncmp(res_name,"CZZ",3)==0 )return(16);
    if(strncmp(res_name,"CYS",3)==0 ||
       strncmp(res_name,"BCS",3)==0 ||
       strncmp(res_name,"CMT",3)==0 )return(17);
    if(strncmp(res_name,"OCS",3)==0 ||
       strncmp(res_name,"CSX",3)==0 ||
       strncmp(res_name,"CSW",3)==0 ||
       strncmp(res_name,"CSD",3)==0 )return(18);
    
    /******** Other atoms ********************/
  }else if(strncmp(atom_name,"SE",2)==0){     
    /* Selenium atom */
    return(19);
    
  }else if(strncmp(atom_name,"P ",2)==0 || strncmp(atom_name,"AS",2)==0){
    /* Phosphorous in Phosphate groups. Also Arsenium atoms*/	 
    return(20);
    
  }else if(strncmp(atom_name,"I1",2)==0 || strncmp(atom_name,"I2",2)==0 ||
	   strncmp(atom_name,"F1",2)==0 || strncmp(atom_name,"F2",2)==0 ||
	   strncmp(atom_name,"F3",2)==0 ||
	   strncmp(atom_name,"IE",2)==0 ||
	   strncmp(atom_name,"CL1D",4)==0 || strncmp(atom_name,"CL1E",4)==0 ||
	   strncmp(atom_name,"CLL",3)==0) { 
    return(23); // halogen atom
    
  }else if(strncmp(atom_name,"BE",2)==0){
    /* Alkaline earth metals */
    return(24);
    
  }

  // Look by residue name for rare atom names
  if(strncmp(res_name,"BCS",3)==0){
    if(strncmp(atom_name,"CZ2",3)==0 ||
       strncmp(atom_name,"CT1",3)==0 ||
       strncmp(atom_name,"CT2",3)==0 ||
       strncmp(atom_name,"CH",2)==0 ) return(5);
    if(strncmp(atom_name,"CZ1",3)==0)return(4);
    
  }else if(strncmp(res_name,"CAF",3)==0){
    if(strncmp(atom_name,"O1",2)==0)return(14);
    if(strncmp(atom_name,"SG",2)==0)return(16);
    
  }else if(strncmp(res_name,"CSY",3)==0){
    if(strncmp(atom_name,"C2",2)==0 ||
       strncmp(atom_name,"C3",2)==0 ||
       strncmp(atom_name,"CG",2)==0 ||
       strncmp(atom_name,"CZ",2)==0)return(1);
    
  }else if(strncmp(res_name,"CS4",3)==0){
    if(strncmp(atom_name,"C1N",3)==0 ||
       strncmp(atom_name,"C1Q",3)==0 ||
       strncmp(atom_name,"C1P",3)==0 ||
       strncmp(atom_name,"C1O",2)==0)return(1);
    if(strncmp(atom_name,"C1G",3)==0 ||
       strncmp(atom_name,"C1H",3)==0 ||
       strncmp(atom_name,"C1I",3)==0)return(4);
    if(strncmp(atom_name,"C1K",3)==0)return(6);
    if(strncmp(atom_name,"C1J",3)==0)return(7);
    if(strncmp(atom_name,"O1C",3)==0)return(14);
    if(strncmp(atom_name,"CL1",3)==0)return(23);

  }else if(strncmp(res_name,"CXM",3)==0){
    if(strncmp(atom_name,"C ",2)==0 ||
       strncmp(atom_name,"CN",2)==0)return(0);
    if(strncmp(atom_name,"ON1",3)==0  ||
       strncmp(atom_name,"ON2",3)==0)return(14);
    if(strncmp(atom_name,"SD ",2)==0)return(16);
    
  }else if(strncmp(res_name,"CZZ",3)==0){
    if(strncmp(atom_name,"O1 ",3)==0)return(15);
    if(strncmp(atom_name,"SG",2)==0)return(16);
    
  }else if(strncmp(res_name,"C12",3)==0){
    if(strncmp(atom_name,"C1",2)==0 ||
       strncmp(atom_name,"CA2",3)==0||
       strncmp(atom_name,"C3",2)==0)return(0);
    if(strncmp(atom_name,"C2 ",2)==0 ||
       strncmp(atom_name,"CZ ",2)==0 ||
       strncmp(atom_name,"CG2",3)==0)return(1);
    if(strncmp(atom_name,"CA1",3)==0)return(2);
    if(strncmp(atom_name,"CB1",3)==0)return(3);
    if(strncmp(atom_name,"CD1",3)==0  ||
       strncmp(atom_name,"CD2",3)==0  ||
       strncmp(atom_name,"CE1",3)==0  ||
       strncmp(atom_name,"CE2",3)==0)return(4);
    if(strncmp(atom_name,"CB2",3)==0 ||
       strncmp(atom_name,"CA3",3)==0)return(6);
    if(strncmp(atom_name,"CG1",3)==0)return(8);
    if(strncmp(atom_name,"N1 ",2)==0)return(10);
    if(strncmp(atom_name,"N3",2)==0  ||
       strncmp(atom_name,"N2",2)==0)return(11);
    if(strncmp(atom_name,"O3",2)==0)return(14);
    if(strncmp(atom_name,"OG1",3)==0 ||
       strncmp(atom_name,"O2",2)==0  ||
       strncmp(atom_name,"OH",2)==0)return(15);

  }else if(strncmp(res_name,"FME",3)==0){
    if(strncmp(atom_name,"CN ",2)==0)return(1);
    if(strncmp(atom_name,"O1 ",3)==0)return(14);
    if(strncmp(atom_name,"SD ",2)==0)return(16);
    
  }else if(strncmp(res_name,"K1R",3)==0){
    if(strncmp(atom_name,"CA1",3)==0 ||
       strncmp(atom_name,"C2", 2)==0)return(1);
    if(strncmp(atom_name,"CA2",3)==0)return(2);
    if(strncmp(atom_name,"CB1",3)==0)return(6);
    if(strncmp(atom_name,"CG1",3)==0 ||
       strncmp(atom_name,"CB2",3)==0 )return(7);
    if(strncmp(atom_name,"N2", 2)==0)return(10);
    if(strncmp(atom_name,"NAF",3)==0)return(11);
    if(strncmp(atom_name,"N1", 2)==0)return(13);
    if(strncmp(atom_name,"O11",3)==0 ||
       strncmp(atom_name,"O12",3)==0 ||
       strncmp(atom_name,"O2", 2)==0)return(14);
    if(strncmp(atom_name,"SG2",3)==0)return(16);
    if(strncmp(atom_name,"SAG",3)==0)return(18);
    if(strncmp(atom_name,"OE1",3)==0)return(22);
    
  }else if(strncmp(res_name,"KYN",3)==0){
    if(strncmp(atom_name,"C1 ",2)==0 ||
       strncmp(atom_name,"CD2",3)==0 ||
       strncmp(atom_name,"CG ",2)==0)return(1);
    if(strncmp(atom_name,"O2 ",2)==0)return(14);
    
  }else if(strncmp(res_name,"LLP",3)==0){
    if(strncmp(atom_name,"OP1",3)==0 ||
       strncmp(atom_name,"OP2",3)==0 ||
       strncmp(atom_name,"OP3",2)==0 ||
       strncmp(atom_name,"OP4",2)==0 ) return(21);
    if(strncmp(atom_name,"C2",2)==0 ||
       strncmp(atom_name,"C3",2)==0 ) return(1);
    if(strncmp(atom_name,"C5",2)==0 ) return(6);
    if(strncmp(atom_name,"C6",2)==0 ) return(4);
    if(strncmp(atom_name,"C4",2)==0 ) return(7);
    if(strncmp(atom_name,"O3",2)==0 ) return(15);

  }else if(strncmp(res_name,"MIR",3)==0){
    if(strncmp(atom_name,"O1P",3)==0 || 
       strncmp(atom_name,"O2P",3)==0 ||
       strncmp(atom_name,"O3P",3)==0) return(21);
    if(strncmp(atom_name,"C1",2)==0) return(6);
    if(strncmp(atom_name,"C2",2)==0) return(8);

  }else if(strncmp(res_name,"MIS",3)==0){
    if(strncmp(atom_name,"OG",2)==0 ||
       strncmp(atom_name,"O1P",3)==0 ||
       strncmp(atom_name,"O2P",3)==0 ||
       strncmp(atom_name,"O3P",3)==0 ) return(21);
    if(strncmp(atom_name,"C1",2)==0 ) return(3);
    if(strncmp(atom_name,"C2",2)==0 ||
       strncmp(atom_name,"C3",2)==0 ) return(8);
    
  }else if(strncmp(res_name,"NEP",3)==0){ // Phospho?
    if(strncmp(atom_name,"O1P",3)==0 ||
       strncmp(atom_name,"O2P",3)==0 ||
       strncmp(atom_name,"O3P",3)==0 ) return(21);
    if(strncmp(atom_name,"P  ",3)==0 ) return(20);

  }else if(strncmp(res_name,"PHD",3)==0){
    if(strncmp(atom_name,"OD1",3)==0 ||
       strncmp(atom_name,"OP1",3)==0 ||
       strncmp(atom_name,"OP2",3)==0 ||
       strncmp(atom_name,"OP3",3)==0 ) return(21);
    
  }else if(strncmp(res_name,"PM3",3)==0){
    if(strncmp(atom_name,"O3",2)==0 || 
       strncmp(atom_name,"O1",3)==0 ||
       strncmp(atom_name,"O2",3)==0 ) return(21);
    if(strncmp(atom_name,"CH4",3)==0 )return(6);
    
  }else if(strncmp(res_name,"PTR",3)==0){
    if(strncmp(atom_name,"OH",2)==0 ||
       strncmp(atom_name,"O1P",3)==0 ||
       strncmp(atom_name,"O2P",3)==0 ||
       strncmp(atom_name,"O3P",3)==0 ) return(21);
    
  }else if(strncmp(res_name,"QPA",3)==0){
    if(strncmp(atom_name,"C10",3)==0)return(1);
    if(strncmp(atom_name,"C8 ",2)==0)return(3);
    if(strncmp(atom_name,"C9 ",2)==0)return(8);
    if(strncmp(atom_name,"O9",2)==0  ||
       strncmp(atom_name,"O10",3)==0)return(14);
    if(strncmp(atom_name,"P2",2)==0)return(20);
    if(strncmp(atom_name,"O11",3)==0 ||
       strncmp(atom_name,"O12",3)==0 ||
       strncmp(atom_name,"O13",3)==0 ||
       strncmp(atom_name,"O14",3)==0)return(21);
    
  }else if(strncmp(res_name,"SAC",3)==0){
    if(strncmp(atom_name,"C1A",3)==0)return(1);
    if(strncmp(atom_name,"C2A",3)==0)return(8);
    if(strncmp(atom_name,"OAC",3)==0)return(14);
    
  }else if(strncmp(res_name,"SEP",3)==0){ // Phosphoserine
    if(strncmp(atom_name,"O1P",3)==0 ||
       strncmp(atom_name,"O2P",3)==0 ||
       strncmp(atom_name,"O3P",3)==0 ) return(21);
    if(strncmp(atom_name,"P  ",3)==0 ) return(20);

  }else if(strncmp(res_name,"TPO",3)==0){ // Phospho?
    if(strncmp(atom_name,"O1P",3)==0 ||
       strncmp(atom_name,"O2P",3)==0 ||
       strncmp(atom_name,"O3P",3)==0 ) return(21);
    if(strncmp(atom_name,"P  ",3)==0 ) return(20);

  }else if(strncmp(res_name,"SGB",3)==0){
    if(strncmp(atom_name,"P1",2)==0) return(20);
    if(strncmp(atom_name,"OG",2)==0 ||
       strncmp(atom_name,"O1",2)==0 ||
       strncmp(atom_name,"O2",2)==0) return(21);
    if(strncmp(atom_name,"C2",2)==0) return(3);
    if(strncmp(atom_name,"C1",2)==0  ||
       strncmp(atom_name,"C4",2)==0  ||
       strncmp(atom_name,"C3",2)==0) return(8);
    
  }else if(strncmp(res_name,"SME",3)==0){
    if(strncmp(atom_name,"S ",2)==0)return(16);
    
  }else if(strncmp(res_name,"SUI",3)==0){
    if(strncmp(atom_name,"C1 ",3)==0)return(1);
    if(strncmp(atom_name,"C1A",3)==0)return(3);
    if(strncmp(atom_name,"C2A",3)==0)return(6);
    if(strncmp(atom_name,"N2",2)==0)return(9);
    if(strncmp(atom_name,"O1 ",3)==0)return(14);
    
  }else if(strncmp(res_name,"UF0",3)==0 ){
    if(strncmp(atom_name,"C1E",3)==0 ||
       strncmp(atom_name,"C1J",3)==0 ||
       strncmp(atom_name,"C1U",2)==0)return(0);
    if(strncmp(atom_name,"C1C",3)==0)return(2);
    if(strncmp(atom_name,"C1D",3)==0)return(3);
    if(strncmp(atom_name,"C1B",3)==0 ||
       strncmp(atom_name,"C1H",3)==0||
       strncmp(atom_name,"C1I",3)==0||
       strncmp(atom_name,"C1M",3)==0||
       strncmp(atom_name,"C1N",3)==0||
       strncmp(atom_name,"C1V",3)==0)return(6);
    if(strncmp(atom_name,"C1P",3)==0  ||
       strncmp(atom_name,"C1Q",2)==0)return(8);
    if(strncmp(atom_name,"N1F",3)==0 ||
       strncmp(atom_name,"N1K",3)==0 ||
       strncmp(atom_name,"N1O",3)==0)return(10);
    if(strncmp(atom_name,"N ",2)==0)return(11);
    if(strncmp(atom_name,"O1G",3)==0  ||
       strncmp(atom_name,"O1L",3)==0  ||
       strncmp(atom_name,"O1W",3)==0)return(14);
    if(strncmp(atom_name,"O1X",3)==0  ||
       strncmp(atom_name,"O1",2)==0)return(15);
    if(strncmp(atom_name,"P1R",3)==0)return(20);
    if(strncmp(atom_name,"O1S",3)==0 ||
       strncmp(atom_name,"O1T",3)==0 ||
       strncmp(atom_name,"O1A",2)==0 ||
       strncmp(atom_name,"O1Y",2)==0)return(21);

  }
  // Residue name not found but standard atom name
  if(strncmp(atom_name,"N ",2)==0)return(10);
  if(strncmp(atom_name,"CA ",3)==0)return(2);
  if(strncmp(atom_name,"C ",2)==0)return(0);
  if(strncmp(atom_name,"O ",2)==0)return(14);
  if(strncmp(atom_name,"CB ",3)==0)return(6);
  
  printf("ERROR, atom %s in residue %s does not exist\n", atom_name, res_name);
  if(file_out){fprintf(file_out, "%s %s\n", res_name, atom_name);}
  return(ATYPES-1);
}


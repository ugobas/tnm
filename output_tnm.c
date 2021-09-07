#include "nma_para.h"
//#include "Parameters.h"
#include "coord.h"
#include "tnm.h"
#include "nma.h"
#include "vector.h"
#include "buildup.h"
#include "output_tnm.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void Print_modes(int N_print, char *nameout,
		 char *label, int *select,
		 int N, float **Mode, float *collectivity,
		 float *sigma2, double sum_sigma2,
		 struct axe *axe, int naxe,
		 atom *atoms, int natoms,
		 struct residue *seq, int nres,
		 int *atom_num, int N_ref)
{
  int PRINT_ALL=0;
  FILE *file_out; int ia, i=0, cart=0, k=0, j=0, ik=0, ires;
  char outfile[200], xyz[4];

  strcpy(xyz, "xyz");
  sprintf(outfile, "%s.Modes_%s.dat", nameout, label);
  file_out=fopen(outfile, "w");
  printf("Writing normal modes to file %s\n", outfile);

  fprintf(file_out, "# %d degrees of freedom: %s\n", N, label);
  fprintf(file_out, "# Contribution to fluctuation: ");
  for(ia=0; ia<N_print; ia++)
    if((select[ia])||PRINT_ALL)
      fprintf(file_out, "\t%.3g", sigma2[ia]/sum_sigma2);
  fprintf(file_out, "\n");
  fprintf(file_out, "# e_values:                    ");
  for(ia=0; ia<N_print; ia++)
    if((select[ia])||PRINT_ALL)
      fprintf(file_out, "\t%.3g", 1./sigma2[ia]);
  fprintf(file_out, "\n");
  fprintf(file_out, "# collectivity:                ");
  for(ia=0; ia<N_print; ia++)
    if((select[ia])||PRINT_ALL)
      fprintf(file_out, "\t%.3g", collectivity[ia]);
  fprintf(file_out, "\n");

  if(strncmp(label, "Cart", 4)==0)cart=1;
  for(i=0; i<N; i++){
    fprintf(file_out, "%4d", i);
    for(ia=0; ia<N_print; ia++)
      if((select[ia])||PRINT_ALL)
	fprintf(file_out, "\t%11.5g", Mode[ia][i]);
    if(cart==0){
      if(axe[i].type=='f'){fprintf(file_out, "\tphi");}
      else if(axe[i].type=='p'){fprintf(file_out, "\tpsi");}
      else if(axe[i].type=='l'){fprintf(file_out, "\tlen");}
      else if(axe[i].type=='a'){fprintf(file_out, "\tban");}
      else if(axe[i].type=='t'){fprintf(file_out, "\ttan");}
      else{ fprintf(file_out, "\t%c", axe[i].type);}
      sscanf(seq[axe[i].bond->atom->res].pdbres, "%d", &ires);
      fprintf(file_out, "%d\n", ires);
    }else{
      k=i/3; j=i-3*k; ik=atom_num[k];
      sscanf(seq[atoms[ik].res].pdbres, "%d", &ires);
      fprintf(file_out, "\t%s%d%c\n", atoms[ik].name, ires, xyz[j]);
    }
  }

  fclose(file_out);

}


int Print_PDB_mode_old(char *nameout, int ia, float *Cart_mode,
		       float Amplitude, float Cart_collectivity,
		       float Tors_collectivity, float MW_Tors_collectivity,
		       float eigen_value, float eigen_B,
		       atom *atoms, int natoms, struct residue *seq, int nres,
		       int *atom_num, int N_ref)
{
//Print_PDB_3(Cart_mode[ia], eigen_value[ia], eigen_B[ia], atoms,
//		  N_ref, atom_ref, seq, file_name, ia);
  FILE *file_out;  char outfile[200];
  //~ FILE * CHAIN[2];
  int i, j, k; atom *atom1=atoms;
  char aaname3[10];
  // Amplitude factors
  //~ float DTHETA=0;
  //~ float AMPL_FACTOR=16; // Amplification with respect to thermal fluctuations
  //~ float AMAX=120;      // New amplitude factor
  //~ float step, dtheta;
  float r[3];
  //~ int N_STEP=10; // Number of conformations per normal mode
  //~ int i_step;
  //~ double norm;

  for(i=0; i<10; i++)aaname3[i]='\0';

  // Opening file
  if(ia<1){
    sprintf(outfile, "%s_modes.pdb", nameout);
    file_out=fopen(outfile, "w");
    printf("Writing normal modes in PDB format in %s\n", outfile);
    fprintf(file_out, "MODEL 0:  Native structure\n");
  }else{
    file_out=fopen(outfile, "a");
    fprintf(file_out, "MODEL %d\n", ia+1);
    fprintf(file_out, "REMARK  Normal mode %3d  Percent fluctuation= %.2f",
	    ia, eigen_B*100.0);
    fprintf(file_out,
	    "REMARK Collectivity: Cartesian %.3f Torsional %.3f MW_Torsional %.3f\n",
	    Cart_collectivity, Tors_collectivity, MW_Tors_collectivity);
  }

  /*
  // Compute amplitude
  DTHETA=1./sqrt(eigen_value);
  fprintf(file_out, " sqrt<theta^2>= %.4f\n", DTHETA);
  DTHETA*=AMPL_FACTOR;
  step=DTHETA/N_STEP; dtheta=0;
  */

  // Compute and print
  int jatom=0;
  for(i=0; i<N_ref; i++){
    int m=3*i;
    atom1=atoms+atom_num[i];
    for(j=0; j<3; j++){
      r[j]=atom1->r[j];
      if(ia>=0)r[j]+=Cart_mode[m+j]*Amplitude;
    }
    k=atom1->res; Name3(aaname3, seq[k].i_aa);
    if(jatom <99999)jatom++;
    fprintf(file_out,
	    "ATOM  %5d%4s  %3s %c%4s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	    jatom, atom1->name, aaname3, chain,
	    seq[k].pdbres,  r[0], r[1], r[2], 1.0, atom1->B_factor);
    atom1=atoms+atom1->i_next;
  }
  fprintf(file_out, "ENDMDL\n");
  fclose(file_out);
  return(0);
}


void Print_mode_summary(char *nameout, char *label, struct Normal_Mode NM,float M_sqrt,int anharmonic,float xkappa)
{ 
  int i;
  FILE *file_out; char outfile[200];
  sprintf(outfile, "%s.%s.dat", nameout, label);
  file_out=fopen(outfile, "w");
  printf("Writing %s\n", outfile);
  float kappa=Collectivity_norm1(NM.sigma2, NM.N);
  fprintf(file_out, "# Reciprocal collectivity of fluctuations=  %.1f\n",kappa);
  fprintf(file_out, "# kappa=  %.3f\n",xkappa);

  double norm=0; for(i=0; i<NM.N; i++)norm+=NM.sigma2[i];
  if(anharmonic){
    double Anhar_ene=0, Anhar_str=0;
    for(i=0; i<NM.N; i++){
      if(NM.sigma2[i]==0)continue;
      Anhar_ene+=NM.Anharmonicity[i]*NM.sigma2[i];
      Anhar_str+=NM.Anharm_struct[i]*NM.sigma2[i];
    }
    fprintf(file_out, "# Frequency weighted anharmonicity (ene)= %.3f\n",
	  Anhar_ene/norm);
    fprintf(file_out, "# Frequency weighted anharmonicity (str)= %.3f\n",
	    Anhar_str/norm);
  }

  fprintf(file_out,"#mode pc_therm cumul ");
  fprintf(file_out,"om-2(harm) ");
  if(anharmonic)fprintf(file_out, "om-2(corr) om-2(anha) ");
  fprintf(file_out," RMSD  C_cart");
  if(NM.MW_Tors_coll)fprintf(file_out, " C_MW");
  if(NM.Tors_coll)fprintf(file_out, " C_tors");
  fprintf(file_out, " Max_dev_atom");
  /*fprintf(file_out, " Anharmonicity_(ene) Anharmonicity_(str)");
    fprintf(file_out, " Max_RMSD(DE<E_THR*E_nat)"); //, E_THR*/
  fprintf(file_out, "\n");

  double sum=0;
  for(i=0; i<NM.N; i++){
    if(NM.select[i]==0)fprintf(file_out, "#");
    sum+=NM.sigma2[i];
    fprintf(file_out, "%5d  %6.4f %5.3f  ",i, NM.sigma2[i]/norm, sum/norm);
    fprintf(file_out, "%7.4g ",NM.sigma2[i]);
    if(anharmonic)
      fprintf(file_out,"%7.4g %7.4g  ",
	      xkappa*NM.sigma2_anhar[i],NM.sigma2_anhar[i]);
    fprintf(file_out, " %7.3g", 1./(M_sqrt*sqrt(NM.omega2[i])));
    fprintf(file_out, "  %5.3f", NM.Cart_coll[i]);
    if(NM.MW_Tors_coll)fprintf(file_out, " %5.3f", NM.MW_Tors_coll[i]);
    if(NM.Tors_coll)fprintf(file_out, " %5.3f", NM.Tors_coll[i]);
    fprintf(file_out, "  %5.3f", NM.Max_dev[i]);
    /*fprintf(file_out, " %.3f", NM.Anharmonicity[i]);
    fprintf(file_out, " %.3f", NM.Anharm_struct[i]);
    fprintf(file_out, " %.3f", NM.Max_RMSD[i]);*/
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}

int Check_make_dir(char *outdir){
  char name_out[1000]; FILE *file_out;
  if(outdir[0]=='\0')return(0);
  sprintf(name_out, "%s/%s", outdir, "tmp");
  file_out=fopen(name_out, "w");
  if(file_out!=NULL){fclose(file_out); return(1);}
  sprintf(name_out, "mkdir -p %s\n", outdir);
  if(system(name_out)==0)return(1);
  return(0);
}

void Print_cart_fluct(int *atom_num, int N_atom, float *fluct,
		      atom *atoms, struct residue *seq,
		      char *name, char *type)
{
  float MAX_B=50.0;
  int i, ia=0, n=1;
  FILE *file_out=fopen(name, "w");
  printf("Writing %s\n", name);
  fprintf(file_out, "REMARK %d B_factor = %s\n",n, type);

  { // Normalize such that max_fluctuation = MAX_B
    float max_f=0, norm;
    for(i=0; i<N_atom; i++)if(fluct[i]>max_f)max_f=fluct[i];
    norm=MAX_B/max_f; for(i=0; i<N_atom; i++)fluct[i]*=norm;
  }
  for(i=0; i<N_atom; i++){
    Print_atom(atoms+atom_num[i], &ia, fluct[i], seq, file_out);
  }
  fprintf(file_out, "TER\n");
  fclose(file_out);
}


void Print_tors_fluct(struct axe *axe, int N_axes, float *fluct,
		      atom *atoms, struct residue *seq,
		      char *name, char *type)
{
  float MAX_B=100.0, f; int i, ia=0, n=1;
  FILE *file_out=fopen(name, "w");
  printf("Writing %s\n", name);
  fprintf(file_out, "MODEL %d  B_factor = %s\n", n, type);

  { // Normalize such that max_fluctuation = MAX_B
    float max_f=0, norm;
    for(i=0; i<N_axes; i++)if(fluct[i]>max_f)max_f=fluct[i];
    norm=MAX_B/max_f;
    for(i=0; i<N_axes; i++)fluct[i]*=norm;
  }

  for(i=0; i<N_axes; i++){
    if(axe[i].type=='f'){
      Print_atom(axe[i].bond->previous->atom, &ia, fluct[i], seq, file_out);
      f=fluct[i]; if(i+1<N_axes)f=0.5*(f+fluct[i+1]);
      Print_atom(axe[i].bond->atom, &ia, f, seq, file_out);
    }else if(axe[i].type=='p'){
      if(ia==0)
	Print_atom(axe[i].bond->previous->atom, &ia, fluct[i], seq, file_out);
      Print_atom(axe[i].bond->atom, &ia, fluct[i], seq, file_out);
    }
  }
  fprintf(file_out, "ENDMDL\n");
  fclose(file_out);
}

void Print_structures(char *pdbout,
		      double *atom_str1, atom *atoms1,
		      struct residue *seq1, char chain1,
		      double *atom_str2, atom *atoms2,
		      struct residue *seq2, char chain2,
		      int *atom_num, int N_ref, int N_cart,
		      float **Cart_mode, float *coeff,
		      int *sort, int N_MODE_PRINT, int N_modes)
{
  double *atom_tmp;
  int imod=0, kmod=0, nmod, i;
  if(N_MODE_PRINT<=0)return;
  if(N_MODE_PRINT < N_modes){nmod = N_MODE_PRINT;}else{nmod= N_modes;}

  // Print reference structure
  //Write_coord(pdbout,atom_str1,atoms1,atom_num,N_ref,seq1,chain1,kmod);
  // kmod++;
  atom_tmp=malloc(N_cart*sizeof(double));
  for(i=0; i<N_cart; i++)atom_tmp[i]=atom_str1[i];

  for(imod=0; imod< nmod; imod++){
    int ik=sort[imod];
    float *mode=Cart_mode[ik];
    /* Update intermediate structure */
    for(i=0;i<N_cart;i++)atom_tmp[i]+=coeff[ik]*mode[i];
    Write_coord(pdbout,atom_tmp,atoms1,atom_num,N_ref,seq1,chain1,kmod);
    kmod++;
  }

  // Print experimental structure 2
  //Write_coord(pdbout,atom_str2,atoms2,atom_num,N_ref,seq2,chain2,-1);
  printf("Writing %s\n", pdbout);
  free(atom_tmp);
}

int Write_coord(char *name_out, double *atoms_ref,
		atom *atoms, int *atom_num, int N_atom_ref,
		struct residue *seq, char chain, int i_model)
{
  int i, ini=0, m=0, ires=0, iatom=0, print_bb=1, nprint_B, nprint_S;
  int *printed=malloc(N_atom_ref*sizeof(int));
  struct residue *res;
  char aaname3[4]; FILE *fh;
  atom *atom1=NULL;

  // printf("Printing PDB %d\n", i_model);
  if(i_model==0){
    fh=fopen(name_out, "w");
  }else{
    fh=fopen(name_out, "a");
  }
  fprintf(fh, "MODEL %d: ", i_model+1);
  if(i_model==0){
    fprintf(fh, "STR.1\n");
  }else if(i_model > 0){
    fprintf(fh, "%d normal modes\n", i_model);
  }else if(i_model<0){
    fprintf(fh, "STR.2\n");
  }
  for(i=0; i<4; i++)aaname3[i]='\0';


  int jatom=0;
  ires=0; ini=0; print_bb=1; iatom=0;
  for (i=ini;i<N_atom_ref;i++)printed[i]=0;
  while(iatom < N_atom_ref){

    nprint_B=0;  nprint_S=0;
    for (i=ini;i<N_atom_ref;i++){
      if(printed[i])continue;
      atom1=atoms+atom_num[i];
      if(atom1->res > ires)break;
      if((strncmp(atom1->name, "N ", 2)==0)||
	 (strncmp(atom1->name, "CA", 2)==0)||
	 (strncmp(atom1->name, "C ", 2)==0)||
	 (strncmp(atom1->name, "O ", 2)==0)){
	if(print_bb){nprint_B=1; break;}
      }else{
	if(print_bb==0){nprint_S=1; break;}
      }
    }
    if((print_bb) && (nprint_B==0)){print_bb=0; continue;}
    if((print_bb==0) && (nprint_S==0)){
      ires=ires+1; ini=i; print_bb=1; continue;
    }

    printed[i]=1; iatom++;
    if(jatom < 99999)jatom++;
    res=seq+atom1->res;
    m=3*i;

    Name3(aaname3, res->i_aa);
    if (strlen(atom1->name)==4){
      fprintf(fh,"%-6s%5d %-4s%1s%3s %c%4d%1s   %8.3f%8.3f%8.3f%6.2f\n",
	      "ATOM", jatom, atom1->name," ", aaname3,chain,
	      //res->pdbres," ",atoms_ref[m],atoms_ref[m+1],atoms_ref[m+2],
	      atom1->res+1," ",atoms_ref[m],atoms_ref[m+1],atoms_ref[m+2],
	      atom1->occupancy);
    }else{
      fprintf(fh,"%-6s%5d  %-3s%1s%3s %c%4d%s   %8.3f%8.3f%8.3f%6.2f\n",
	      "ATOM", jatom, atom1->name," ",aaname3,chain,
	      //res->pdbres," ",atoms_ref[m],atoms_ref[m+1],atoms_ref[m+2],
	      atom1->res+1," ",atoms_ref[m],atoms_ref[m+1],atoms_ref[m+2],
	      atom1->occupancy);
    }
  }
  fprintf(fh,"TER\n");
  fclose(fh); free(printed);
  return (0);
}
void Print_atom(atom *atom, int *num, float B, struct residue *seq,
		FILE *file_out)
{
  struct residue *s=seq+atom->res;
  char aaname3[4]; Name3(aaname3, s->i_aa);
  (*num)++; int jatom=*num; if(jatom>99999)jatom=99999;
  fprintf(file_out,
	  "ATOM  %5d  %4s%3s %c%4d   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	  jatom, atom->name, aaname3, 'A', atom->res+1,
	  atom->r[0], atom->r[1], atom->r[2], 1.0, B);
}

int Print_change(float *Fluct_pred, float *Fluct_obs, int N,
		char *nameout, char *what)
{
  char namefile[200]; FILE *file_out; int i;
  sprintf(namefile, "%s_%s_dev.dat", nameout, what);
  file_out=fopen(namefile, "w");
  printf("Writing %s\n", namefile);

  fprintf(file_out, "# Squared deviation\n# Pred Obs\n");
  for(i=0; i<N; i++)
    fprintf(file_out, "%.4f %.4f\n", Fluct_pred[i], Fluct_obs[i]);
  fclose(file_out);
  return(0);
}

void Initialize(int Ndir, struct axe *axes, int naxes, int N_modes,
                atom *atoms, int natoms,
                atom *atoms2, int N_ref,
                int *atom_num1, int *atom_num2,
		struct Para_simul Para_simul,
		struct bond *BONDS)
{
  int N_coord=3*N_ref, i, dir;
   printf("Initialization of bonds\n");
   Set_bonds_measure(BONDS, natoms, atoms);
   //for(i=0; i<natoms; i++)printf("%.0f ", BONDS[i].r[0]); printf("\n");
   ATOM_ALL1=malloc(3*natoms*sizeof(float));
   Put_coord(ATOM_ALL1, BONDS, natoms);
   ATOM_REF1=malloc(N_coord*sizeof(float));
   Write_ref_coord_float(ATOM_REF1, N_ref, ATOM_ALL1, atom_num1);
   ATOM_PDB=malloc(N_coord*sizeof(float));
   for(i=0; i<N_coord; i++)ATOM_PDB[i]=ATOM_REF1[i];
   MASS_ALL=malloc(natoms*sizeof(float));
   for(i=0; i<natoms; i++)MASS_ALL[i]=Mass(atoms+i);
   MASS_REF=malloc(N_ref*sizeof(float));
   for(i=0; i<N_ref; i++)MASS_REF[i] = MASS_ALL[atom_num1[i]];
   if(Ini_dthr==0){
     dthr=Para_simul.D_REP; dthr2=dthr*dthr; Ini_dthr=1;
   }
   Ene_ini=Energy_clashes(ATOM_ALL1, natoms);
   printf("Initial energy= %.3g\n", Ene_ini);
   Ene_act=Ene_ini;
   T_all=Para_simul.E_THR*Ene_ini;
   T_step=T_all/N_modes;
   Ene_thr_step =Ene_ini+T_step; Ene_thr_all=Ene_ini+T_all;
   E_thr1=Ene_thr_step;
   //E_thr1=E_THR1*Ene_ini/N_modes;
   N_modes_active=N_modes;
   for(dir=0; dir<Ndir; dir++){
     d_phi[dir]=malloc(naxes*sizeof(float));
     bonds[dir]=malloc((natoms+1)*sizeof(struct bond));
     ATOM_ALL[dir]=malloc(3*natoms*sizeof(float));
     ATOM_REF[dir]=malloc(N_coord*sizeof(float));
   }
   bonds_act=malloc((natoms+1)*sizeof(struct bond));
   Copy_bonds(bonds_act, BONDS, natoms);
   if(atoms2!=NULL){
     ATOM_REF2=malloc(N_coord*sizeof(float));
     Write_ref_coord_atom(ATOM_REF2, N_ref, atoms2, atom_num2);
     RMSDif=rmsd_mclachlan_f(ATOM_REF1, ATOM_REF2, MASS_REF, N_ref);
     k_thr=(Para_simul.E_THR*Ene_ini)/pow(RMSDif, EXP_FORCE);
     RMS_NORM=sqrt((float)N_ref)*RMSDif;
     ATOM_DIFF=malloc(N_coord*sizeof(float));
     for(i=0; i<N_coord; i++)ATOM_DIFF[i]=ATOM_REF2[i]-ATOM_REF1[i];
   }
}


int Move_struct_confchange(atom *atoms, int natoms,
			   float *rmsd_fin,
			   char *nameout, int ia,
			   int PRINT,
			   struct axe *axes, int naxes,
			   struct Normal_Mode *NM,
			   struct residue *seq, int nres,
			   atom *atoms2, int N_ref,
			   int *atom_num1, int *atom_num2,
			   struct Para_simul Para_simul,
			   struct bond *BONDS)
{
  // WARNING: atoms is updated here!
  if(NM->contr2fluct[ia]<=0)return(0);

  int i, j, N_coord=3*N_ref, dir;
  FILE *file_out=NULL;
  
  // Initialize native structure and output files
  if(INI_MOVE==0){
    Initialize(2, axes, naxes, NM->N, atoms, natoms, atoms2, N_ref,
	       atom_num1, atom_num2, Para_simul, BONDS);
    INI_MOVE=1;
    if(PRINT){        // Opening file for printing PDB
      sprintf(pdbout, "%s_confchange_all_modes.pdb", nameout);
      file_out=fopen(pdbout, "w");
      printf("Writing trajectory in PDB format in %s\n", pdbout);
      Print_PDB(file_out, atoms, natoms, ATOM_ALL1, seq, 0, 0.0);
      fclose(file_out); k_accept=1;
    }
    rmsd_ini_prev=0;
    if(atoms2!=NULL)*rmsd_fin=RMSDif;
    sprintf(namermsd, "%s_rmsd.dat", nameout);
    file_rmsd=fopen(namermsd, "w");
    printf("Opening %s\n", namermsd);
    fprintf(file_rmsd,"# E_THR= %.2f Eini= %.2g k_thr= %.3g",
	    Para_simul.E_THR*Ene_ini, Ene_ini, k_thr);
     if(atoms2!=NULL)fprintf(file_rmsd," RMSD= %.3f", RMSDif);
     fprintf(file_rmsd,"\n# Force_const (1,-1) RMSD_ini (1,-1)");
     if(atoms2!=NULL){
       fprintf(file_rmsd," RMSD_fin (1,-1)");
       fprintf(file_rmsd," RMSD_force\n");
     }
     fclose(file_rmsd);
  }
  if(Ini_dthr==0){
    dthr=Para_simul.D_REP; dthr2=dthr*dthr; Ini_dthr=1;
  }


  // Determine maximum possible factor for mode ia
  float *mode=NM->Tors[ia];
  float factor=NM->Max_factor[ia];
  if(factor==0){
    factor=Para_simul.AMPLITUDE/NM->omega[ia];
    float max_v=0;
    for(i=0; i<naxes; i++)if(fabs(mode[i])>max_v)max_v=fabs(mode[i]);
    float Max_factor=Para_simul.MAX_ANGLE/max_v;
    if(factor > Max_factor){
      printf("WARNING, max. amplitude allowed for mode %d = %.3f < %.3f\n",
	     ia, Max_factor, factor);
      factor=Max_factor;
    }
    NM->Max_factor[ia]=factor;
  }
  /*if(atoms2 !=NULL){
    float c=Scalar_product(ATOM_DIFF, NM->Cart[ia], N_coord)/RMS_NORM;
    if(factor > fabs(c))factor=fabs(c);
    }*/

  // Move in both directions
  float Ene[2], rmsd_i[2], rmsd_f[2], force[2], *dphi, *A_REF;
  for(dir=0; dir<2; dir++){
    dphi=d_phi[dir]; A_REF=ATOM_REF[dir];
    for(i=0; i<naxes; i++)dphi[i]=mode[i]*factor;
    Copy_bonds(bonds[dir], bonds_act, natoms);
    Build_up(bonds[dir], natoms, dphi, naxes);
    Put_coord(ATOM_ALL[dir], bonds[dir], natoms);
    Ene[dir]=Energy_clashes(ATOM_ALL[dir], natoms);
    Write_ref_coord_float(A_REF,N_ref, ATOM_ALL[dir], atom_num1);
    rmsd_i[dir]=rmsd_mclachlan_f(ATOM_REF1, A_REF, MASS_REF, N_ref);
    if(atoms2 !=NULL){
      rmsd_f[dir]=rmsd_mclachlan_f(ATOM_REF2, A_REF, MASS_REF, N_ref);
    }
    force[dir]=(Ene[dir]-Ene_ini)/(rmsd_i[dir]*rmsd_i[dir]);
    factor=-factor;
  }

  // Store anharmonicity
  if(NM->Anharmonicity[ia]==0){
    NM->Anharmonicity[ia]=(force[0]-force[1])/fabs(Ene_ini);
  }
  if(NM->Anharm_struct[ia]==0){
    NM->Anharm_struct[ia]=(rmsd_i[0]-rmsd_i[1])/(rmsd_i[0]+rmsd_i[1]);
  }

  // Write energy and rmsd
  if(PRINT_ENE){
    fprintf(file_rmsd,"# mode %d contr2fluc= %.3f factor=%.2f\n",
	    ia, NM->contr2fluct[ia], NM->Max_factor[ia]);
    file_rmsd=fopen(namermsd, "a");
    fprintf(file_rmsd,"%8.3g\t%8.3g\t%.2f\t%.2f",
	    force[0], force[1], rmsd_i[0], rmsd_i[1]);
    if(atoms2 != NULL){
      fprintf(file_rmsd,"\t%.2f\t%.2f", rmsd_f[0], rmsd_f[1]);
      if(((force[0] < force[1])&&(rmsd_f[0] < rmsd_f[1]))||
	 ((force[1] < force[0])&&(rmsd_f[1] < rmsd_f[0]))){
	fprintf(file_rmsd,"\t1");
      }else{
	fprintf(file_rmsd,"\t0");
      }
    }
    fprintf(file_rmsd,"\n");
    fclose(file_rmsd);
  }


  // Select direction
  if(atoms2 !=NULL){ // Based on distance
    if(rmsd_f[0] < rmsd_f[1]){
      dir=0;
    }else{
      dir=1;
    }
  }else{             // Based on energy
    if(force[1] > force[0]){dir=0;}
    else{dir=1;}
  }

  // Accept move based on energy ?
  if(force[dir] > k_thr)return(0);
  // Accept move based on rmsd ?
  if((rmsd_f[dir] > *rmsd_fin)&&(rmsd_i[dir] < rmsd_ini_prev))return(0);
  // Store
  Copy_bonds(bonds_act, bonds[dir], natoms); Ene_act=Ene[dir];
  *rmsd_fin=rmsd_f[dir]; rmsd_ini_prev=rmsd_i[dir];

  // Update coordinates
  float *atom_coord=ATOM_ALL[dir]; int k=0;
  for(i=0; i<natoms; i++){
    for(j=0; j<3; j++){atoms[i].r[j]=atom_coord[k]; k++;}
  }

  // Print PDB
  if(PRINT){
    float *A_REF=ATOM_REF[dir];
    float rmsd=rmsd_mclachlan_f(ATOM_PDB, A_REF, MASS_REF, N_ref);
    if(rmsd > Para_simul.RMSD_STEP){
      for(i=0; i<N_coord; i++)ATOM_PDB[i]=A_REF[i];
      file_out=fopen(pdbout, "a");
      rmsd_mclachlan_f(ATOM_ALL1, ATOM_ALL[dir], MASS_ALL, natoms);
      Print_PDB(file_out,atoms,natoms,ATOM_ALL[dir],seq,k_accept,rmsd_i[dir]);
      fclose(file_out);
      k_accept++;
    }
  }

  return(1);
}

int Move_struct_print(char *nameout, int ia, float Amplitude,
		      int Nstep, int PRINT, int DIR,
		      struct axe *axes, int naxes,
		      struct Normal_Mode NM,
		      atom *atoms, int natoms,
		      struct residue *seq, int nres,
		      struct bond *bonds)
{
//Print_PDB_3(Cart_mode[ia], eigen_value[ia], eigen_B[ia], atoms,
//		  N_ref, atom_ref, seq, file_name, ia);

  if(NM.contr2fluct[ia]<=0)return(-1);
  FILE *file_out=NULL; int i, k;

  // Opening file
  if(PRINT){
    char outfile[200];
    sprintf(outfile, "%s_mode%d.pdb", nameout, ia);
    file_out=fopen(outfile, "w");
    printf("Writing normal modes in PDB format in %s\n", outfile);
    fprintf(file_out, "REMARK  Normal mode %3d  Percent fluctuation= %.2f\n",
	    ia, NM.contr2fluct[ia]*100.0);
    fprintf(file_out,
	    "REMARK Collectivity: Cartes. %.3f Tors. %.3f MW_Tors. %.3f\n",
	    NM.Cart_coll[ia], NM.Tors_coll[ia], NM.MW_Tors_coll[ia]);
  }

  Set_bonds_measure(bonds, natoms, atoms);
  float factor=DIR*Amplitude/NM.omega[ia];
  float *d_phi=malloc(naxes*sizeof(float));
  for(i=0; i<naxes; i++)d_phi[i]=NM.Tors[ia][i]*factor;

  float *atom_str1=malloc(3*natoms*sizeof(float));
  Put_coord(atom_str1, bonds, natoms);

  float *atom_str3=malloc(3*natoms*sizeof(float));
  float *mass_atom=Set_masses(atoms,natoms);

  // Decide direction to go ???
  // int pos=Decide_direction(bonds, d_phi, atom_str3, atoms, naxes,  natoms);

  if(PRINT)Print_PDB(file_out, atoms, natoms, atom_str1, seq, 0, 0.0);
  for(k=1; k<=2*Nstep; k++){
    Build_up(bonds, natoms, d_phi, naxes); //if(k>1)
    Put_coord(atom_str3, bonds, natoms);
    float rmsd=rmsd_mclachlan_f(atom_str1, atom_str3, mass_atom, natoms);
    if(PRINT)Print_PDB(file_out, atoms, natoms, atom_str3, seq, k, rmsd);
    if(k==Nstep)for(i=0; i<naxes; i++)d_phi[i]=-d_phi[i];
  }
  free(atom_str1); free(atom_str3); free(mass_atom);
  free(d_phi); if(PRINT)fclose(file_out);
  return(0);
}

float Simulation_old(atom *atoms_sim, atom *atoms, int natoms,
		 struct axe *axes, int naxes,
		 float **Tors_mode, float *omega_inverse, int N_modes,
		 double SDEV_INI, double RMSD,
		 struct Para_simul Para_simul,
		 struct bond *bonds)
{
  float SDEV=SDEV_INI/13, rmsd;
  float *d_phi=malloc(naxes*sizeof(float));
  float RMSD_MIN=RMSD*(1.-TOL), RMSD_MAX=RMSD*(1.+TOL);
  int MAX_RMSD_INTER = 10, rmsd_iter = 0, i, j;

  float *atom1=malloc(3*natoms*sizeof(float));
  float *atom2=malloc(3*natoms*sizeof(float));
  int num_atom=Put_coord(atom1, bonds, natoms);
  float *mass=Set_masses(atoms, natoms);
  Set_bonds_measure(bonds, natoms, atoms);

  while(rmsd_iter < MAX_RMSD_INTER){
    Extract_dphi_old(d_phi, Tors_mode, omega_inverse, N_modes, naxes, SDEV);
    Build_up(bonds, num_atom, d_phi, naxes);
    // Evaluate rmsd
    Put_coord(atom2, bonds, natoms);
    rmsd=rmsd_mclachlan_f(atom1, atom2, mass, num_atom);
    printf("Simulation: rmsd= %5.2f A=%7.0f\n", rmsd, SDEV);
    if(rmsd<RMSD_MIN){
      SDEV*=1.4;
    }else if(rmsd > RMSD_MAX){
      SDEV/=1.25;
    }else{
      printf("-> confchange\n"); break;
    }
    rmsd_iter++;
  }

  int n=0;
  for(i=0; i<natoms; i++){
    if(atoms[i].ali==0)continue;
    for(j=0; j<3; j++)atoms_sim[n].r[j]=bonds[n].r[j];
    n++;
  }

  free(d_phi); free(mass);
  free(atom1); free(atom2);
  return(rmsd);
}

void Write_all_coord(double *coord, atom *atoms, int natoms)
{
  int i, j, k=0;
  for(i=0; i<natoms; i++){
    for(j=0; j<3; j++){coord[k]=atoms[i].r[j]; k++;}
  }
}

void Extract_dphi_old(float *d_phi, float **Tors_mode, float *omega_inverse,
		      int N_modes, int N_axes, float SDEV_SIM)
{
  int i, k;

  // Initialize random numbers
  if(INIRAN==0){
    INIRAN=1;
    unsigned long iran=randomgenerator();
    InitRandom( (RANDOMTYPE)iran);
    RANFACTOR=pow(12.0,1/3.0);
  }

  double sum=0;
  for(i=0; i<N_axes; i++)d_phi[i]=0;
  for(k=0; k<N_modes; k++){
    float ran=RandomFloating()-0.5;
    float coeff = ran*omega_inverse[k];
    //if(DEBUG)printf("%.3f\n", ran);
    for(i=0; i<N_axes; i++)d_phi[i]+=coeff*Tors_mode[k][i];
    sum+=coeff*coeff;
  }
  float norm=sqrt(SDEV_SIM/sum);
  for(i=0; i<N_axes; i++)d_phi[i]*=norm;

  //if(DEBUG)
  //printf("sum c**2= %.3f, N*rmsd**2=%.3f\n", sum*norm*norm, SDEV_SIM);

}


int Move_all_modes(char *nameout,
		   struct axe *axes, int naxes, struct Normal_Mode *NM,
		   atom *atoms, int natoms, struct residue *seq, int nres,
		   atom *atoms2, int N_ref,
		   int *atom_num1, int *atom_num2,
		   int nstep,
		   struct Para_simul Para_simul,
		   struct bond *bonds)
{
  FILE *file_out;
  int i, ia;
  //int N_coord=3*N_ref;
  char outfile_pdb[200], outfile_rmsd[200];

  /*if(INI_MOVE==0){
     printf("Initialization of Move_all_modes\n");
     Initialize(2, axes, naxes, NM->N, atoms, natoms, atoms2, N_ref,
     		atom_num1, atom_num2, Para_simul, BONDS);
     INI_MOVE=1;
     }*/

  // Opening file
  if(INI_ALL==0){
    Ene_act=Ene_ini;
    if(file_rmsd!=NULL)fclose(file_rmsd);
    sprintf(outfile_rmsd, "%s_Allmodes_rmsd.dat", nameout);
    file_rmsd=fopen(outfile_rmsd, "w");
    printf("Writing energy and RMSD in %s\n", outfile_rmsd);
    fprintf(file_rmsd,
            "# E_THR=%.3f E_thr1= %.3f T(E-E_ini)= %.3g  T(E-E')= %.3g\n",
            Para_simul.E_THR, E_thr1, T_all, T_step);
    fprintf(file_rmsd,"#RMSD(C,ini)\tEne");
    if(atoms2!=NULL)fprintf(file_rmsd,"\tRMSD(C,fin)");
    fprintf(file_rmsd,"\tmode\tdir\n");
    fprintf(file_rmsd,"%.3f\t%8.3g", 0.0,  0.0);
    if(atoms2!=NULL)fprintf(file_rmsd,"\t%.3f", RMSDif);
    fprintf(file_rmsd,"\t-1\t-1\n");
    //
    sprintf(outfile_pdb, "%s_Allmodes.pdb", nameout);
    file_out=fopen(outfile_pdb, "w");
    //printf("Writing displacements in PDB format in %s\n", outfile_pdb);
    //Print_PDB(file_out, atoms, natoms, ATOM_ALL1, seq, 0, 0.0);
    INI_ALL=1;
  }else{
    file_out=fopen(outfile_pdb, "a");
    file_rmsd=fopen(outfile_rmsd, "a");
  }

  // Simulation
  float factor=Para_simul.AMPLITUDE;
  if(Para_simul.RESET)factor*=nstep;
  if(Ini_dthr==0){
    dthr=Para_simul.D_REP; dthr2=dthr*dthr; Ini_dthr=1;
  }
  if(PRINT_ENE){
    fprintf(file_rmsd, "# Normal mode multiplier= %.2f/omega\n", factor);
  }
  float rmsd_i[2], Ene[2], Emax, *A_REF;
  int dir;
  k_accept=0;
  if(Para_simul.RESET){
    Ene_act=Ene_ini; Ene_thr_step=Ene_act+T_step;
    Copy_bonds(bonds_act, BONDS, natoms);
  }
  for(dir=0; dir<2; dir++)
    Copy_bonds(bonds[dir], bonds_act, natoms);

  for(ia=0; ia<NM->N; ia++){
    if(NM->contr2fluct[ia]==0)continue;

    // Compute factor
    //if(NM->Max_factor[ia] <AMPLITUDE)factor=NM->Max_factor[ia];
    float factor_ia=factor/NM->omega[ia];

    for(dir=0; dir<2; dir++){
      float *dphi=d_phi[dir], *mode=NM->Tors[ia];
      A_REF=ATOM_REF[dir];
      for(i=0; i<naxes; i++)dphi[i]=mode[i]*factor_ia;
      Build_up(bonds[dir], natoms, dphi, naxes);
      Put_coord(ATOM_ALL[dir], bonds[dir], natoms);
      Ene[dir]=Energy_clashes(ATOM_ALL[dir], natoms);
      Write_ref_coord_float(A_REF, N_ref, ATOM_ALL[dir], atom_num1);
      rmsd_i[dir]=rmsd_mclachlan_f(ATOM_REF1,A_REF,MASS_REF,N_ref);
      factor_ia=-factor_ia;
    }

    if(Ene[0]<Ene[1]){
      dir=0; Emax=Ene[1];
    }else{
      dir=1; Emax=Ene[0];
    }

    if((NM->Anharmonicity[ia]==0)&&(Emax > E_thr1)){
      NM->Anharmonicity[ia]=fabs(Ene[0]-Ene[1])/fabs(Ene_ini);
    }
    if(NM->Anharm_struct[ia]==0){
      NM->Anharm_struct[ia]=(rmsd_i[0]-rmsd_i[0])/(rmsd_i[0]+rmsd_i[0]);
    }


    // Evaluate energy
    //if(Metropolis(Ene[dir], Ene_act, T_step)&&
    //   Metropolis(Ene[dir], Ene_ini, T_all)){
    if((Ene[dir]< Ene_thr_step)&&(Ene[dir] < Ene_thr_all)){
      k_accept++;
      Ene_act=Ene[dir]; Ene_thr_step=Ene_act+T_step;
      Copy_bonds(bonds_act, bonds[dir], natoms);
      if(dir==0){Copy_bonds(bonds[1], bonds_act, natoms);}
      else{Copy_bonds(bonds[0], bonds_act, natoms);}
      if(PRINT_ENE)fprintf(file_rmsd, "%8.3g\t%.3f", Ene_act, rmsd_i[dir]);
      if(atoms2 !=NULL){
	float rmsd_f=rmsd_mclachlan_f(ATOM_REF2,ATOM_REF[dir],MASS_REF,N_ref);
	if(PRINT_ENE)fprintf(file_rmsd, "\t%.3f", rmsd_f);
      }
      if(PRINT_ENE)fprintf(file_rmsd, "\t%3d\n", ia);
      //Print_PDB(file_out,atoms,natoms,ATOM_ALL[dir],seq,k_accept,rmsd);
    }else{
      for(dir=0; dir<2; dir++)Copy_bonds(bonds[dir], bonds_act, natoms);
    }
  }
  fclose(file_out); fclose(file_rmsd);
  return(k_accept);
}


void Binding_sites(int *site,
		   float *Prof_dir, float **Coup_dir,
		   float *Prof_coord, float **Coup_coord,
		   float *Prof_proximity, float **Coup_proximity, 
		   float *Prof_deformation, float **Coup_deformation,
		   //float *Prof_deformation_dist,float **Coup_deformation_dist,
		   float *Prof_dr_dir, float **Coup_dr_dir,
		   int *atomres, char **pdbres, char *amm, 
		   int Na, char *chain, struct residue *seq,
		   int Nres, char *pdb, char *nameout1, char *SITES);

  /*** Check if binding site exists ***/
  {
    char SEL[4]="CA";
    int Nref=Ref_kin.N_ref, iref[Nref], iatom[Nref], atomres[Nref];
    float mass=0;
    int Na=Extract_atoms(iref,iatom,atomres,seq1,&mass,atoms1,Ref_kin,SEL);
    int site[Na];
    char *pdbres[Na], amm[Na];
    for(int ia=0; ia<Na; ia++){
      struct residue *res=seq1+(atoms1+iatom[ia])->res;
      pdbres[ia]=malloc(6*sizeof(char));
      strcpy(pdbres[ia], res->pdbres);
      amm[ia]=res->amm;
    }
    Binding_sites(site,
		  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
		  atomres, pdbres, amm, Na,
		  chain1, seq1, nres1, file_pdb1, nameout1, SITES);
  }

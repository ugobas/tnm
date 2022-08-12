#ifndef ATOMS_H
#define ATOMS_H

extern int num_BB; //=5
extern char BB_name[8][3]; // N CA C CB O
//                       C-N    N-CA   CA-C   CA-CB  C-O  
extern float BB_length[8]; //={1.3301,1.4572,1.5236,1.5319,1.2316};
//                    CA-C-N C-N-CA N-CA-C N-CA-CB CA-C-O  
extern float BB_angle[8]; //={1.1091,1.0196,1.2018,1.2162,1.0359};

int Atom_type(char *name);

/*
  CYS  SG
  ASP  CG  OD1 OD2
  GLU  CG  CD  OE1 OE2
  PHE  CG  CD1 CD2 CE1 CE2=CZ  (CZ connected to both CE1, CE2)
  HIS  CG  ND1 CD2 CE1=NE2 (CE1 connected to NE2)
  ILE  CG1 CG2 ND1
  LYS  CG  CD  CE  NZ
  LEU  CG  CD1 CD2
  MET  CG  SD  CE
  ASN  CG  OD1 ND2
  PRO  CG  CD  (CD is connected to N)
  GLN  CG  CD  OE1 NE2
  ARG  CG  CD  NE  CZ  NH1  NH2
  SER  OG
  THR  OG1 CG2
  VAL  CG1 CG2
  TRP  CG  CD1 CD2 NE1=CE2 CE3 CZ2 CZ3=CH2 (CH2 connected to both CZ2,CZ3)   
  TYR  CG  CD1 CD2 CE1 CE2=CZ  OH  (CZ connected to both CE1, CE2)
 */

#endif

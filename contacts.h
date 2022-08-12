#ifndef __INI_CONTACTS
#define __INI_CONTACTS


float Contact_overlap(struct interaction *Int_list1, int N_int1,
		      struct interaction *Int_list2, int N_int2, int *alires);
float Contact_energy(struct interaction *Int_list, int N_int,
		     struct residue *seq);
void Get_contact_list(int **nc1, int ***clist1, int ***cnum1, int Nres,
		      atom *atoms, struct interaction *Int_list, int N_int);
extern float **Econt; 

#endif

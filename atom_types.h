#define ATYPES 26
#define MAX_CLOSEST_NEIGHBORS 5000
#define DICAM_PARAM 4
#define DICAM_CLOSEST_NEIGHBORS 300
#include <stdio.h>

extern char atom_names[ATYPES][6];
extern int ATOM_PAIRS;	    

int get_atom_type(char *atom_name, char *res_name, FILE *file_out);
int find_pair_indx(int, int);
int Get_atom_pair(int a1, int a2);
int read_dicam_library(char *libfile, float **dicam, char **names);
void Atom_name(char *name, int n);

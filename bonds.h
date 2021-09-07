// bond lengths and elastic contacts for AMBER force field
// Based on ff03 CA=C  C=C2 N OH
// Second column based on http://www.cryst.bbk.ac.uk/PPS95/course/3_geometry/peptide2.html 1995
/*
#Bonds
"C N " 1.45 1.32
"N CA" 1.45 1.47
"CAC " 1.53 1.53
"CACB" 1.53 1.53
"C O " 1.41 1.24
# Angles
"CA C  N " 116.6 114
"C  N  CA" 121.9 123
"N  CA C " 116.6 123
"N  CA CB" 116.6 ?
"CA C  O " 121   121
"O  C  N " = 360 - (CA-C-O + CA-C-N) = 122
"CB CA C " ?
*/
/* Fixed degrees of freedom of backbone = 5atoms*3-3torsions = 12
   5 bond lengths, 2 torsions of CB and O, 5 bond angles
 */
int  num_BB=5;
char BB_name[]="N  CA C  CB O  "; //N=0 CA=1 C=2 CB=3 O=4
float BB_bond[num_BB]={1.45,1.45,1.53,1.53,1.41};
float BB_angle[num_BB]={116.6,121.9,116.6,116.6,121};




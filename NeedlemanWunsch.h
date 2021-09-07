int alignNW(char *seq1, 
	    int  length1, 
	    char *seq2, 
	    int  length2, 
	    short verbose, 
	    short identity, 
	    int  penalty, 
	    char *align1, 
	    char *align2,
	    int  *align_len);
int CalcMDMScore(char resa, char resb);
int CalcMDMScoreUC(char resa, char resb);

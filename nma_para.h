// Reference atoms
#define ALL_AXES_DEF 1 // If 0, first phi and last psi angles excluded
#define SEL_DEF "CA"
char REF[20];

int HYD; // Are hydrogen atoms considered?
int HYD_REF;


// Computation of native interactions and force constant in
// interactions_tnm.c and anharmonic_tnm.c
// k_ij = KAPPA*(r_ij/C0)^-EXP_HESSIAN if r<C1 0 if r>=C_THR (if POW)
int POW;      // Power law (1) or exponential (0)?
float C_THR;  // Threshold for contacts
float EXP_HESSIAN; 
float C1_THR; // Start of interpolation region
float C0;     // Reference distance for force constant
float KAPPA;  // Force constant at r=C0
double K0; // K0=KAPPA*C0^EXP;
double K1; // K1=K0*C1^-EXP/(C_THR-C1);
double K_OMEGA, K_PSI, K_PHI, K_CHI, K_BA, K_BL;
// Force constants for internal variables

// Fit of B factors
int FIT_B;
float R_MIN;     // Minimum correlation above which B-factors are fitted
float SLOPE_MIN; // Minimum slope above which B-factors are fitted
float RMSD_EXP;

// Energy function
float lambda;
float E_repulsion;
float D_Res[20][20], U_Res[20][20];


// defaulta values for interaction thresholds
#define THR_A_MIN   9.0   //7.0
#define THR_A_MAX   9.0    //14.0
#define THR_A_STEP  1.0
#define THR_B_MIN   9.0   // 7.0
#define THR_B_MAX   9.0   // 14.0
#define THR_B_STEP  1.0
#define THR_C_MIN   5.0   // 4.0
#define THR_C_MAX   5.0   // 7.0
#define THR_C_STEP  0.5

#define THR_A 9.0
#define THR_B 9.0
#define THR_C 4.5

#define N_MODES_DEF 0 // Lowest frequency modes printed
int N_MODES;

#define E_MIN_DEF   0.00000001  // Minimum relative eigenvalue 
float E_MIN;  // Minimum eigenvalue allowed for normal modes

#define VERBOSE 1      /* Set to one for more output */
#define KINETIC_DEF 1 // If(KINETIC), compute T_inv
int KINETIC;      // Considering kinetic energy in tnm model

#define PRINT_LAMBDA_DEF 0
int PRINT_LAMBDA; // Printing list of eigenvalues

#define ANISOU_DEF 0 // Examine anisotropic temperature factors?
int ANISOU;

// Hydrogen bonds
#define ENE_HB_DEF 2.0  // Ratio between HB interactions and other types
float ENE_HB;
#define DA_THR_DEF    3.5    // Max donor-acceptor dist in h bond
float DA_THR;
#define COS_DAA_DEF  -0.08   // Max cos(A'AD) (min. 105 degrees)
float COS_DAA;

#define THR_NO_DEF 4.0   // Threshold for NO interactions
float THR_NO;

#define N_INT_TYPE 4  // if 4, computes also hydrogen bonds
// Interactions types: 1=CA, 2=CB, 3=all, 4=HB

#define FIT_ROT 0 // Include global rotations in B factors fit? (3 param)

#define PRINT_AXES_DEF 1
int PRINT_AXES;

int COMPARE;

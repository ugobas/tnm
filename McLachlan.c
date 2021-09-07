/*************************************************************************

   File: McLachlan.c
         
	 This files contains a subroutine to perform McLachlan's least
	 square fit, as implemented in the program Profit by 
*************************************************************************/

/* Includes
*/
#include "McLachlan.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"
/************************************************************************/
/* Defines and macros
*/

#define SMALL  1.0e-20     /* Convergence cutoffs                       */
#define SMALSN 1.0e-10
#define DBG 0

#ifndef PI
#define PI (4.0 * atan(1.0))
#endif

#ifndef _MATHTYPE_H
#define _MATHTYPE_H
#endif

/* This is for compilers running on machines such as Amigas, Macs and
   older Sun workstations using 680X0 series processors with maths
   coprocessors. This assumes that the symbol _M68881 is defined when
   the compiler is run to use the maths coprocessor and that a file
   called m68881.h is to be included to make full use of the coprocessor
*/
#ifdef _M68881
#include <m68881.h>
#endif

/************************************************************************/
/* Prototypes
*/

static short matfit(double *x1, double *x2, double rm[3][3],
             int n, double *wt1, short column); 
static void qikfit(double umat[3][3], double rm[3][3], short column);
static void center_of_mass(double *ret, double *v, int N, double *mass);
static void TranslateVector(double *v, double transl[3], int N);
static void Apply3Matrix2dot(double *pt, double matrix[3][3]);
static double rmsd_w (double *alignv1, double *alignv2, int N, double *mass);

/************************************************************************/
/*>double rmsd_maclachlan(double *x1, double *x2, double *wt1, int n)
   ------------------------------------------------------------------
   Input:   double *x1         First (fixed) array of coordinates
            double *x2         Second (mobile) array of coordinates
	    double *wt1        Weight array or NULL
	    int    n           Number of atoms to be fitted
   
   Returns: double             rmsd
   
   
   This subroutine returns the rmsd after least squares superposition of
   x2 onto x1 following McLachlan's algorithm (1982) Acta Cryst. 
   A38 871-873; as implemented in the program Profit
			       
*/			       


double rmsd_mclachlan(double *x1, double *x2, double *wt1, int n)
{
    
  int i=0,j=0,m=0, ncoord=3*n;
  double cent1[3], cent2[3], rm[3][3];
  double tmpcent[3];
  double rmsd;

  double tmpatom1[ncoord], *x=x1, *t=tmpatom1;
  for(i=0; i<ncoord; i++){*t=*x; x++; t++;}

  if(DBG){
    double sum=0;
    for(i=0; i<n; i++)sum+=wt1[i];
    printf("Average mass: %.2f\n", sum/n);
  }

  // compute centre of masses for the reference
  // and the mobile proteins    
  center_of_mass(cent1,x1,n,wt1);
  center_of_mass(cent2,x2,n,wt1);
  
  if(DBG){
    printf("Centre of Mass for str1\n");
    printf("\t%.4f %.4f %.4f\n",cent1[0],cent1[1],cent1[2]);
    printf("Centre of Mass for str2\n");
    printf("\t%.4f %.4f %.4f\n",cent2[0],cent2[1],cent2[2]);
  }
  
  // centre both structures in their centres of masses
  tmpcent[0] = -cent1[0];
  tmpcent[1] = -cent1[1];
  tmpcent[2] = -cent1[2];
  
  TranslateVector(tmpatom1,tmpcent,n);
  
  tmpcent[0] = -cent2[0];
  tmpcent[1] = -cent2[1];
  tmpcent[2] = -cent2[2];
  
  TranslateVector(x2,tmpcent,n);
    
  // compute the superposition
  matfit(tmpatom1,x2,rm,n,wt1,1);
  
  // Print rotation matrix
  if(DBG){
    printf("Rotation Matrix\n");
    for (i=0;i<3;i++){
      printf("\t");
      for (j=0;j<3;j++)printf("%.4f ",rm[i][j]);
      printf("\n");
    }
  }
 
  for (i=0;i<n;i++){  
    Apply3Matrix2dot(x2+m,rm);
    m+=3;
  }
  TranslateVector(x2,cent1,n);
    
  if(DBG){
    printf("Translation vector (between Centres of masses)\n");
    printf("\t%.4f %.4f %.4f\n",cent1[0]-cent2[0],
	   cent1[1]-cent2[1],cent1[2]-cent2[2]);
  }
  
  // computes weighted rmsd
  rmsd=rmsd_w(x1,x2,n,wt1);
  
  if(DBG)printf("\t mass-weighted RMS: %.4f\n", rmsd);
  
  free(tmpatom1);
  
  return (rmsd);
  
}

/************************************************************************/
/*>short matfit(double *x1, double *x2, double rm[3][3], int n,
               double *wt1, short column)
   -----------------------------------------------------
   Input:   double  *x1         First (fixed) array of coordinates
            double  *x2         Second (mobile) array of coordinates
            int     n           Number of atoms
            double  *wt1        Weight array or NULL
            short   column      TRUE: Output a column-wise matrix (as used
                                 by FRODO)
                                FALSE: Output a standard row-wise matrix.
   Output:  double  rm[3][3]    Returned rotation matrix
   Returns: double              TRUE:  success
                                FALSE: error

   Fit coordinate array x2 to x1 both centred around the origin and of 
   length n. Optionally weighted with the wt1 array if wt1 is not NULL.
   If column is set the matrix will be returned column-wise rather 
   than row-wise.

*/
short matfit(double    *x1,        /* First coord array    */
             double    *x2,        /* Second coord array   */
             double    rm[3][3],   /* Rotation matrix      */
             int       n,          /* Number of points     */
             double    *wt1,       /* Weight array         */
             short     column)     /* Column-wise output   */
{
   int  i,j,m=0;
   double umat[3][3];

   
   if(n<2){return(0);}

   if(wt1){
      for(i=0;i<3;i++){
	for(j=0;j<3;j++) umat[i][j] = 0.0;
	m=0;
	for(j=0;j<n;j++){
	  switch(i){
	  case 0:
	    umat[i][0] += wt1[j] * x1[m] * x2[m];
	    umat[i][1] += wt1[j] * x1[m] * x2[m+1];
	    umat[i][2] += wt1[j] * x1[m] * x2[m+2];
	    break;
	  case 1:
	    umat[i][0] += wt1[j] * x1[m+1] * x2[m];
	    umat[i][1] += wt1[j] * x1[m+1] * x2[m+1];
	    umat[i][2] += wt1[j] * x1[m+1] * x2[m+2];
	    break;
	  case 2:
	    umat[i][0] += wt1[j] * x1[m+2] * x2[m];
	    umat[i][1] += wt1[j] * x1[m+2] * x2[m+1];
	    umat[i][2] += wt1[j] * x1[m+2] * x2[m+2];
	    break;
	  }
	  m+=3;
	}
      }
   }else{
     for(i=0;i<3;i++){
       for(j=0;j<3;j++) umat[i][j] = 0.0;
       for(j=0;j<n;j++){
	 switch(i)
	   {
	   case 0:
	     umat[i][0] += x1[m] * x2[m];
	     umat[i][1] += x1[m] * x2[m+1];
	     umat[i][2] += x1[m] * x2[m+2];
	     break;
	   case 1:
	     umat[i][0] += x1[m+1] * x2[m];
	     umat[i][1] += x1[m+1] * x2[m+1];
	     umat[i][2] += x1[m+1] * x2[m+2];
	     break;
	   case 2:
	     umat[i][0] += x1[m+2] * x2[m];
	     umat[i][1] += x1[m+2] * x2[m+1];
	     umat[i][2] += x1[m+2] * x2[m+2];
	     break;
	   }
	 m+=3;
       }
     }
   }
   qikfit(umat,rm,column);
   return(1);
}
/************************************************************************/
/*>static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column)
   ---------------------------------------------------------------
   Input:   double  umat[3][3]     The U matrix
            short   column         TRUE: Create a column-wise matrix
                                  (other way round from normal).
   Output:  double  rm[3][3]       The output rotation matrix
  
*/
static void qikfit(double  umat[3][3],
                   double  rm[3][3],
                   short   column)
{
   
   double  rot[3][3],
     turmat[3][3],
     c[3][3],
     coup[3],
     dir[3],
     step[3],
     v[3],
     rtsum,rtsump,
     //rsum,
     stp,stcoup,
     ud,tr,ta,cs,sn,ac,
     delta,//deltap,
     gfac,
     cle,clep=0.0;
   int i,j,k,l,m,
     jmax,
     ncyc,
     nsteep,
     nrem;

   /* Rotate repeatedly to reduce couple about initial direction to zero.
      Clear the rotation matrix
   */
   for(l=0;l<3;l++)
   {
      for(m=0;m<3;m++)
         rot[l][m] = 0.0;
      rot[l][l] = 1.0;
   }

   /* Copy vmat[][] (sp) into umat[][] (dp)                             */
   jmax = 30;
   rtsum = umat[0][0] + umat[1][1] + umat[2][2];
   delta = 0.0;
   for(i=0;i<3;i++) step[i]=0.0;

   for(ncyc=0;ncyc<jmax;ncyc++)
   {
      /* Modified CG. For first and every NSTEEP cycles, set previous
         step as zero and do an SD step
      */
      nsteep = 3;
      nrem = ncyc-nsteep*(int)(ncyc/nsteep);

      if(!nrem)
      {
         for(i=0;i<3;i++) step[i]=0.0;
         clep = 1.0;
      }
      
      /* Couple                                                         */
      coup[0] = umat[1][2]-umat[2][1];
      coup[1] = umat[2][0]-umat[0][2];
      coup[2] = umat[0][1]-umat[1][0];
      cle     = sqrt(coup[0]*coup[0] + coup[1]*coup[1] + coup[2]*coup[2]);

      /* Gradient vector is now -coup                                   */
      gfac = (cle/clep)*(cle/clep);

      /* Value of rtsum from previous step                              */
      rtsump = rtsum;
      //deltap = delta;
      clep   = cle;
      if(cle < SMALL) break;

      /* Step vector conjugate to  previous                             */
      stp = 0.0;
      for(i=0;i<3;i++)
      {
         step[i]=coup[i]+step[i]*gfac;
         stp   += (step[i] * step[i]);
      }
      stp = 1.0/sqrt(stp);
         
      /* Normalised step                                                */
      for(i=0;i<3;i++) dir[i] = stp*step[i];

      /* Couple resolved along step direction                           */
      stcoup = coup[0]*dir[0] + coup[1]*dir[1] + coup[2]*dir[2];

      /* Component of UMAT along direction                              */
      ud = 0.0;
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            ud += umat[l][m]*dir[l]*dir[m];


      tr = umat[0][0]+umat[1][1]+umat[2][2]-ud;
      ta = sqrt(tr*tr + stcoup*stcoup);
      cs=tr/ta;
      sn=stcoup/ta;
         
      /* If cs<0 then posiiton is unstable, so don't stop               */
      if((cs>0.0) && (ABS(sn)<SMALSN)) break;
            
      /* Turn matrix for correcting rotation:

         Symmetric part
      */
      ac = 1.0-cs;
      for(l=0;l<3;l++)
      {
         v[l] = ac*dir[l];
         for(m=0;m<3;m++)
            turmat[l][m] = v[l]*dir[m];
         turmat[l][l] += cs;
         v[l]=dir[l]*sn;
      }

      /* Asymmetric part                                                */
      turmat[0][1] -= v[2];
      turmat[1][2] -= v[0];
      turmat[2][0] -= v[1];
      turmat[1][0] += v[2];
      turmat[2][1] += v[0];
      turmat[0][2] += v[1];

      /* Update total rotation matrix                                   */
      for(l=0;l<3;l++)
      {
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*rot[k][m];
         }
      }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            rot[l][m] = c[l][m];

      /* Update umat tensor                                             */
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*umat[k][m];
         }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            umat[l][m] = c[l][m];

      rtsum = umat[0][0] + umat[1][1] + umat[2][2];
      delta = rtsum - rtsump;

      /* If no improvement in this cycle then stop                      */
      if(ABS(delta)<SMALL) break;

      /* Next cycle                                                     */
   }

   //rsum = rtsum;

   /* Copy rotation matrix for output                                   */
   if(column)
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[j][i] = rot[i][j];
   }
   else
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[i][j] = rot[i][j];
   }
}

static void center_of_mass(double *ret, double *v, int N, double *mass)
{
  int i, m=0; double norm=0;

  for(i=0;i<3;i++)ret[i] = 0;
  for(i=0;i<N;i++){
    double w=mass[i];
    ret[0] += w*v[m];
    ret[1] += w*v[m+1];
    ret[2] += w*v[m+2];
    norm+=w; m+=3;
  }
  for(i=0;i<3;i++)ret[i]/= norm;
}
static void TranslateVector(double *v, double *transl, int N)
{ 
    int i,m=0;
    for (i=0;i<N;i++){
       v[m]+=transl[0];
       v[m+1]+=transl[1];
       v[m+2]+=transl[2];
       m+=3;
    } 
}
static double rmsd_w (double *alignv1, double *alignv2, int N, double *mass){
  
  double answer =0.0, norm=0.0;
  double d, *x1=alignv1, *x2=alignv2;
  int i, j; 

  for(i=0;i<N;i++){
    double w=mass[i]; norm+=w;
    for(j=0; j<3; j++){
      d = *x1-*x2; x1++; x2++;
      answer += w*d*d;
    }
  }
  return sqrt(answer/norm);
}

static void Apply3Matrix2dot(double *pt, double matrix[3][3])
{
  double v[3]; //={0.0,0.0,0.0};   int i,j;
  for(int i=0;i<3;i++){
    //v[i]=0.0; for (j=0;j<3;j++){v[i]+=matrix[i][j] * pt[j];}
    v[i]=matrix[i][0]*pt[0]+matrix[i][1]*pt[1]+matrix[i][2]*pt[2];
  }
  pt[0]=v[0];
  pt[1]=v[1];
  pt[2]=v[2];
  
}

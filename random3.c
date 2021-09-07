/*
	File: random3.c
	
	Deals with rumber numbers (four-tap shift-register).
	
	Revision: Jan 17, 2007
*/

#include <stdio.h>

#include "random3.h"
#include "math.h"

/* limit */
#define RANDOMLIMIT				(RANDOMTYPE)2147483647
#define INVRANDOMLIMITFLOATING	(double)(1.0/2147483648.0)

/* multiplier */
#define MULTIPLIER				(RANDOMTYPE)16807

/* the array of random numbers */
RANDOMTYPE random_array[16384];

/* next random number */
RANDOMTYPE *random_which1, *random_which2, *random_which3, *random_which4, *random_which5;

/* end of the array */
RANDOMTYPE *random_end;

# define TWOPI 6.283185307

void InitRandom( RANDOMTYPE seed)
{
  int i, j;
  RANDOMTYPE random, composite;
  
  /* check size of RANDOMTYPE */
  if( sizeof( RANDOMTYPE) != 8)
    {
      fprintf( stderr, "Wrong size for RANDOMTYPE.\n");
    }
  
  /* create array of random numbers */
  random= seed;
  for( i= 0; i < 16384; i++)
    {
      composite= (RANDOMTYPE)0;
      for( j= 0; j < 31; j++)
	{
	  /* overflow ? */
	  if( (random *= MULTIPLIER) > RANDOMLIMIT){
	    random -= (RANDOMLIMIT+1);
	  }	  
	  /* use only most significant bit */
	  composite= (composite << 1) | ((random >> 30) & 1);
	}
      random_array[i]= composite;
    }
  
  /* initialize the pointer */
  random_which1= random_array+10000;
  random_which2= random_array+10000-471;
  random_which3= random_array+10000-1586;
  random_which4= random_array+10000-6988;
  random_which5= random_array+10000-9689;
  random_end= random_array+16384;
}

RANDOMTYPE Random()
{
  /* create random number */
  if( ++random_which1 == random_end)
    {
      random_which1= random_array;
    }
  if( ++random_which2 == random_end)
    {
      random_which2= random_array;
    }
  if( ++random_which3 == random_end)
    {
      random_which3= random_array;
    }
  if( ++random_which4 == random_end)
    {
      random_which4= random_array;
    }
  if( ++random_which5 == random_end)
    {
      random_which5= random_array;
    }
  
  return( *random_which1=
	  *random_which2 ^ *random_which3 ^ *random_which4 ^ *random_which5);
}

double RandomFloating()
{
  /* create random number */
  if( ++random_which1 == random_end)
    {
      random_which1= random_array;
    }
  if( ++random_which2 == random_end)
    {
      random_which2= random_array;
    }
  if( ++random_which3 == random_end)
    {
      random_which3= random_array;
    }
  if( ++random_which4 == random_end)
    {
      random_which4= random_array;
    }
  if( ++random_which5 == random_end)
    {
      random_which5= random_array;
    }
  
  return( INVRANDOMLIMITFLOATING*
	  (double)(*random_which1=
		   *random_which2 ^ *random_which3 ^ 
		   *random_which4 ^ *random_which5));
}

float Norm_var_2(){
  // Draws squared Gaussian variable z^2 with Box-Muller transform method
  float U1=RandomFloating();
  float U2=RandomFloating();
  float c=cos(TWOPI*U2);
  return(-2*log(U1)*c*c);
}

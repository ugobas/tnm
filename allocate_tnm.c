#include "coord.h"
#include "tnm.h"
#include <stdlib.h>

void Allocate_tors(struct Tors *X, int N_axes, int N_cart, int N_modes)
{
  X->N_axes=N_axes;
  X->N_cart=N_cart;
  X->Cart=malloc(N_cart*sizeof(float));
  X->Tors=malloc(N_axes*sizeof(float));
  X->MW_Tors=malloc(N_axes*sizeof(float));
  X->coeff=malloc(N_modes*sizeof(float));
}

void Empty_tors(struct Tors X){
  free(X.Tors);
  free(X.MW_Tors);
  free(X.Cart);
  free(X.coeff);
}

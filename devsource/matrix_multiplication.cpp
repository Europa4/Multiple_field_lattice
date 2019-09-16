#include "Prot.h"

void matrix_multiplication(dcomp delta[Ntot], dcomp Mat[Ntot][Ntot], dcomp eta[Ntot])
{
  for (int c = 0; c < Ntot; ++c)
  {
    delta[c] = 0;
    for (int r = 0; r < Ntot; ++r)
    {
      delta[c] += Mat[c][r]*eta[r];
    }
  }
}

void matrix_multiplication(dcomp delta[Ntot], dcomp eta[Ntot], dcomp Mat[Ntot][Ntot])
{
  for (int c = 0; c < Ntot; ++c)
  {
    delta[c] = 0;
    for (int r = 0; r < Ntot; ++r)
    {
      delta[c] += Mat[r][c]*eta[r];
    }
  }
}

dcomp dot_product(dcomp v1[Ntot], dcomp v2[Ntot])
{
  dcomp output = 0;
  for (int i = 0; i < Ntot; ++i)
  {
    output += v1[i]*v2[i];
  }
  return output;
}

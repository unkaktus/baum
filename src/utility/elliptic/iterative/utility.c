/* utility.c */
/* Bernd Bruegmann 01/03 */

#include "bam.h"
#include "iterative.h"




void minus(tVarList *vlw, tVarList *vlu, tVarList *vlv)
{
  int i, j;
  double *u, *v, *w;

  for (j = 0; j < vlu->n; j++) {
    u = VLPtr(vlu, j);
    v = VLPtr(vlv, j);
    w = VLPtr(vlw, j);

    forallinner(vlu->level, i) 
      w[i] = u[i] - v[i];
  }
}




double dot(tVarList *vlu, tVarList *vlv)
{
  int i, j;
  double *u, *v;
  double sum = 0;

  return bampi_allreduce_alldot(vlu, vlv);

  for (j = 0; j < vlu->n; j++) {
    u = vlu->level->v[vlu->index[j]];
    v = vlv->level->v[vlv->index[j]];

    forinner1(vlu->level, i) 
      sum += u[i] * v[i];
  }
  return sum;
}




double norm2(tVarList *u)
{
  return bampi_allreduce_allnorm2(u);
}



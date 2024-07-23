/* jacobi.c */
/* Bernd Bruegmann 01/00 */

#include "bam.h"
#include "iterative.h"


/* Jacobi iteration
   NOTE: this version only works for the flat Laplace operator
         we should pass in dL/du to generalize
*/
int jacobi(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *coeffs,
	   int itmax, double tol, double *normres,
	   void (*lop)(tL *, tVarList *, tVarList *), 
	   void (*precon)(tL *, tVarList *, tVarList *))
{
  int i, j, k;
  int pr = Getv("ell_verbose", "yes");
  double c = - 1 / (2/(l->dx*l->dx) + 2/(l->dy*l->dy) + 2/(l->dz*l->dz));
  
  if (pr) printf("jacobi:  itmax %d, tol %e\n", itmax, tol);

  for (k = 0; k < itmax; k++) {
    lop(l, r, x);
    minus(r, b, r);

    for (j = 0; j < x->n; j++) {
      double *xp = VLPtr(x, j);
      double *rp = VLPtr(r, j);

      forallinner(l, i)
	xp[i] += c * rp[i];
    }

    *normres = bampi_allreduce_allnorm2(r);
    if (pr) printf("%5d  %10.3e %10.3e  residuum\n",
		   k, log10(*normres), *normres);
    if (*normres <= tol) break;
  }

  return 0;
}




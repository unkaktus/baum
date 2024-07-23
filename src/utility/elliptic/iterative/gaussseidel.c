/* gaussseidel.c */
/* Bernd Bruegmann 4/03 */

#include "bam.h"
#include "iterative.h"




int gaussseidel(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *coeffs,
		int itmax, double tol, double *normres,
		void (*lop)(tL *, tVarList *, tVarList *), 
		void (*gs)(tL *, tVarList *, tVarList *))
{
  int i, j, k;
  int pr = Getv("ell_verbose", "yes");
  tVarList *y  = AddDuplicateEnable(x, "_iterative_y");
  
  if (pr) printf("gaussseidel:  itmax %d, tol %e\n", itmax, tol);

  for (k = 0; k < itmax; k++) {
    lop(l, r, x);
    minus(r, b, r);
    *normres = bampi_allreduce_allnorm2(r);
    if (pr) printf("%5d  residuum %10.3e\n", k, *normres);
    if (*normres <= tol) break;

    //prvare(l, VarName(x->index[0]));
    //prdivider(0);
    gs(l, y, x);
    vlcopy(x, y);
    //prvare(l, VarName(x->index[0]));
  }

  return 0;
}




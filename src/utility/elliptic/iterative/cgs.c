/* cgs.c */
/* Bernd Bruegmann 01/00 */

/* CGS: preconditioned conjugate gradient squared method 
   cmp. "Templates for the solution of linear systems", Barett et al
        http://www.netlib.org/templates
*/

#include "bam.h"
#include "iterative.h"




int cgs(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c,
	int itmax, double tol, double *res,
	void (*lop)(tL *, tVarList *, tVarList *), 
	void (*precon)(tL *, tVarList *, tVarList *))
{
  errorexit(
    "cgs has not been adapted to the new version of bam, use bicgstab.");

#if 0
  double *u, *p, *q, *uh, *ph, *rt;
  double alpha, beta, bnorm2;
  double rho = 0, rho1 = 1, rhotol = 1e-50;
  int i, ii;
  int pr = Getv("ell_verbose", "yes");

  /* check */
  if (pr) printf("cgs:  itmax %d, tol %e\n", itmax, tol);
  nvector = l->nnodes;

  /* temporary storage */
  u = dmalloc(6*nvector);
  p =  u+1*nvector;
  q =  u+2*nvector;
  uh = u+3*nvector;
  ph = u+4*nvector;
  rt = u+5*nvector;

  /* compute initial residual rt = r = b - A x */
  lop(l, r, x);
  forvector(i) rt[i] = r[i] = b[i] - r[i];
  *res = norm2(r);
  if (pr) printf("%5d  %10.3e\n", 0, *res);
  if (0) prvare(l, "r");
  if (*res <= tol) return;

  /* norm of rhs for weighting */
  bnorm2 = norm2(b);
  if (bnorm2 == 0) bnorm2 = 1;

  /* cgs iteration */
  for (ii = 0; ii < itmax; ii++) {

    rho = dot(rt, r);
    if (fabs(rho) < rhotol) break;

    /* compute direction vectors p and u */
    if (ii == 0) {
      forvector(i) p[i] = u[i] = r[i];
    } else {
      beta = rho/rho1;
      forvector(i) {
	u[i] = r[i] + beta * q[i];
	p[i] = u[i] + beta * (q[i] + beta * p[i]);
      }
    }

    /* compute direction adjusting scalar alpha */
    precon(l, ph, p);
    lop(l, uh, ph);
    alpha = rho/dot(rt, uh);
    forvector(i) {
      q[i] = u[i] - alpha * uh[i];
      ph[i] = u[i] + q[i];
    }

    /* compute direction adjusting vector */
    precon(l, uh, ph);
    lop(l, ph, uh);
    forvector(i) {
      x[i] += alpha * uh[i];
      r[i] -= alpha * ph[i];
    }

    /* are we done? */
    *res = norm2(r);
    if (0) printf("%5d  %10.3e\n", ii+1, *res);
    if (pr) printf("%5d  %10.3e  %10.3e  %10.3e  %10.3e\n", 
		  ii+1, *res, alpha, beta, rho);
    if (*res <= tol) break;
    rho1 = rho;
  }

  /* free temporary storage */
  free(u);

  /* iteration failed */
  if (ii > itmax) return 1;

  /* breakdown */
  if (fabs(rho) < rhotol) return -10;
#endif
  /* success! */
  return 0;
}



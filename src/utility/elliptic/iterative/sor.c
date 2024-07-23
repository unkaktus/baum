/* sor.c */
/* Bernd Bruegmann 01/00 */

/* SOR: successive overrelaxation method 
   cmp. "Templates for the solution of linear systems", Barett et al
        http://www.netlib.org/templates
*/

#include "bam.h"
#include "iterative.h"




int sor(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c,
	int itmax, double tol, double *res,
	void (*lop)(tL *, tVarList *, tVarList *), 
	void (*precon)(tL *, tVarList *, tVarList *))
{
#if 0
  double omega, rho2;
  int i, it, j;
  int pr = Getv("ell_verbose", "yes");
  double ot, normrt, *rt, *xt, min;

  /* check */
  if (pr) printf("sor:  itmax %d, tol %e\n", itmax, tol);
  nvector = l->nnodes;

  /* temporary storage */
  xt = dmalloc(2*nvector);
  rt = xt + nvector;

  /* compute initial residual r = b - A x */
  lop(l, r, x);
  forvector(i) r[i] = b[i] - r[i];
  *res = norm2(r);
  if (pr) printf("%5d  %10.3e\n", 0, *res);
  if (0) prvare(l, "r");
  if (*res <= tol) return;

  /* diagonal preconditioning is identical to Gauss-Seidel weighting */
  DPflag = 1;
  precon(l, r, r);

  /* inital overrelaxation parameters */
  rho2 = cos(PI*l->ds);
  rho2 *= rho2;
  omega = 2 - 2*PI*l->ds;
  omega = 2 / (1 + sqrt(1 - rho2));
  omega = 1.0;
  if (omega < 1) omega = 1;
  if (omega > 2) omega = 2;

  /* sor iteration */
  for (it = 0; it < itmax; it++) {

    /* update */
    forvector(i) x[i] += omega * r[i];
    if (0) prvare(l, "u"); 

    /* are we done? */
    lop(l, r, x);
    forvector(i) r[i] = b[i] - r[i];
    *res = norm2(r);
    precon(l, r, r);

    /* awfully crude hack to find next omega */
    if (0) {
      min = 1e50;
      ot = omega - l->ds/10;
      for (j = 0; j < 4; j++) {
	forvector(i) xt[i] = x[i] + ot * r[i];
	lop(l, rt, xt);
	forvector(i) rt[i] = b[i] - rt[i];
	normrt = norm2(rt);
	if (normrt < min) {min = normrt; omega = ot;}
	if (0) printf("ot %e  normrt %e\n", ot, normrt);
	ot += l->ds/10;
      }
    }

    if (pr) printf("%5d  %10.3e   %10.3e\n", it+1, *res, omega);
    if (*res <= tol) break;
  }

  /* free temporary storage */
  free(xt);

  /* iteration failed */
  if (it > itmax) return 1;
#endif
  /* success! */
  return 0;
}





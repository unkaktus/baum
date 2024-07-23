/* bicgstab.c */
/* Bernd Bruegmann 01/00 */

/* BiCGSTAB: preconditioned biconjugate gradient stabilized method 
   cmp. "Templates for the solution of linear systems", Barett et al
        http://www.netlib.org/templates
*/

#include "bam.h"
#include "iterative.h"




int bicgstab(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c,
	     int itmax, double tol, double *normres,
	     void (*lop)(tL *, tVarList *, tVarList *), 
	     void (*precon)(tL *, tVarList *, tVarList *))
{
  double alpha = 0, beta = 0;
  double rho = 0, rho1 = 1, rhotol = 1e-50;
  double omega = 0, omegatol = 1e-50;
  int i, ii, j;
  int pr = Getv("ell_verbose", "yes");

  /* temporary storage */
  tVarList *p  = AddDuplicateEnable(x, "_iterative_p");
  tVarList *ph = AddDuplicateEnable(x, "_iterative_ph");
  tVarList *rt = AddDuplicateEnable(x, "_iterative_rt");
  tVarList *s  = AddDuplicateEnable(x, "_iterative_s");
  tVarList *sh = AddDuplicateEnable(x, "_iterative_sh");
  tVarList *t  = AddDuplicateEnable(x, "_iterative_t");
  tVarList *v  = AddDuplicateEnable(x, "_iterative_v");
     
#if 0
  /* hack for scalar only: variables passed to lop need boundary info */
  if (VarFallOff(x->index[0])) {
    VarNameSetBoundaryInfo("iterative_ph", 0, 1, 0);
    VarNameSetBoundaryInfo("iterative_sh", 0, 1, 0);
  }
#endif

  /* check */
  if (pr) printf("bicgstab:  itmax %d, tol %e\n", itmax, tol);

  /* compute initial residual rt = r = b - A x */
  lop(l, r, x);
  for (j = 0; j < r->n; j++) {
    double *prt = VLPtr(rt, j);
    double *pr  = VLPtr(r,  j);
    double *pb  = VLPtr(b,  j);
    forallinner(l, i) prt[i] = pr[i] = pb[i] - pr[i];
  }
  *normres = norm2(r);
  if (pr) printf("bicgstab: %5d  %10.3e\n", 0, *normres);
  if (0) prvare(l, "r");
  if (*normres <= tol) return 0;

  /* cgs iteration */
  for (ii = 0; ii < itmax; ii++) {
    rho = dot(rt, r);
    if (fabs(rho) < rhotol) break;

    /* compute direction vector p */
    if (ii == 0)
      vlcopy(p, r);
    else {
      beta = (rho/rho1)*(alpha/omega);
      for (j = 0; j < r->n; j++) {
	double *pp = VLPtr(p, j);
	double *rp = VLPtr(r, j);
	double *vp = VLPtr(v, j);
	forallinner(l, i) 
	  pp[i] = rp[i] + beta * (pp[i] - omega * vp[i]);
      }
    }

    /* compute direction adjusting vector ph and scalar alpha */
    precon(l, ph, p);
    lop(l, v, ph);
    alpha = rho/dot(rt, v);
    for (j = 0; j < r->n; j++) {
      double *sp = VLPtr(s, j);
      double *rp = VLPtr(r, j);
      double *vp = VLPtr(v, j);
      forallinner(l, i) 
	sp[i] = rp[i] - alpha * vp[i];
    }

    /* early check of tolerance */
    *normres = norm2(s);
    if (*normres <= tol) {
      for (j = 0; j < r->n; j++) {
	double *xp  = VLPtr(x,  j);
	double *php = VLPtr(ph, j);
	forallinner(l, i) 
	  xp[i] += alpha * php[i];
      }
      if (pr) printf("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n", 
		     ii+1, *normres, alpha, beta, omega);
      break;
    }

    /* compute stabilizer vector sh and scalar omega */
    precon(l, sh, s);
    lop(l, t, sh);
    omega = dot(t, s) / dot (t, t);

    /* compute new solution approximation */
    for (j = 0; j < r->n; j++) {
      double *rp = VLPtr(r, j);
      double *sp = VLPtr(s, j);
      double *tp = VLPtr(t, j);
      double *xp = VLPtr(x, j);
      double *php = VLPtr(ph, j);
      double *shp = VLPtr(sh, j);
      forallinner(l, i) {
	xp[i] += alpha * php[i] + omega * shp[i];
	rp[i] = sp[i] - omega * tp[i];
      }
    }

    /* are we done? */
    *normres = norm2(r);
    if (0) printf("bicgstab: %5d  %10.3e\n", ii+1, *normres);
    if (pr) printf("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n", 
		  ii+1, *normres, alpha, beta, omega);
    if (*normres <= tol) break;
    rho1 = rho;
    if (fabs(omega) < omegatol) break;
  }

  /* free temporary storage */
  VLDisableFree(p);
  VLDisableFree(ph);
  VLDisableFree(rt);
  VLDisableFree(s);
  VLDisableFree(sh);
  VLDisableFree(t);
  VLDisableFree(v);

  /* iteration failed */
  if (ii > itmax) return -1;

  /* breakdown */
  if (fabs(rho) < rhotol) return -10;
  if (fabs(omega) < omegatol) return -11;

  /* success! */
  return ii+1;
}




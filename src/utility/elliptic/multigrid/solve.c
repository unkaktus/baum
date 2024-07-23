/* solve.c */
/* Bernd Bruegmann 4/97, 3/03 */

#include "bam.h"
#include "multigrid.h"


double tolres2, tolresi, factau2, factaui;




/* Solve by performing full V multigrid cycles. 
   I.e. call multigrid solve for coarsest level,
   then prolong to one level finer and call multigrid for that finer level,
   then prolong to one level finer and so forth until finest is reached.
   This is called full multigrid. 
   There also is an option to start with the finest level.
*/
void fullVdriver(tL *level, double tol)
{
  int vnmin = Geti("mg_vcycle_nmin");
  int vnmax = Geti("mg_vcycle_nmax");
  int wnmin = Geti("mg_wcycle_nmin");
  int wnmax = Geti("mg_wcycle_nmax");
  int npre  = Geti("mg_relax_npre");
  int npost = Geti("mg_relax_npost");
  int l, lmin, lmax, n;
  int lvcycle, nvcycle;
  double res2, resi, tau2, taui;

  /* set tolerances (not fully implemented yet in new bam) */
  tolres2 = tol;
  tolresi = 0;
  factau2 = 0;
  factaui = 0;

  /* cycle up? */
  lmax = level->grid->lmax;
  lmin = Getv("mg_typeofcycles", "fullV") ? 0 : lmax;

  /* for all levels in multigrid stack */
  for (l = lmin; l <= lmax; l++) {
    lvcycle = l;
    
    /* this is the place where in the old BAM new refinements may be created */
    /* let's not get into that now */

    /* perform V cycles */
    for (n = 0; n < vnmax; n++) {
      nvcycle = n;

      /* are we done? */
      if (n >= vnmin) {
	
	/* check errors */
	res2 = resi = tau2 = taui = 0;
	restau(l, (tolres2)?&res2:0, (tolresi)?&resi:0, 
	          (factau2)?&tau2:0, (factaui)?&taui:0, 0);
	if (res2 <= tolres2 && 
	    resi <= tolresi &&
	    (res2*factau2 <= tau2) &&
	    (resi*factaui <= taui))
	  break;
      }

      /* print */
      if (vb == 2) {
	prresiduals(l); /* if factau not 1, will be incorrect on uH */
	printf("----- level %d, V-cycle %d -------------------------\n\n", 
	       l, n+1);
      }

      /* solve */
      solve(l, lmax, npre, npost, wnmin, wnmax);

      /* for l == 0, solve is exact */
      if (l == 0) break;
    }

    /* print */
    if (vb == 1 && l == lmax) printf(
      "------------------------ after %d V-cycles ---------------------\n",n);
    if (vb == 2 || vb == 1 && l == lmax) prresiduals(l);

    /* for full multigrid, prolong u for refined initial guess,
       f should be zero from stack initialization
    */
    if (l < lmax)
      prolong_u(l, l+1);
  }
}





/* solve 
   recursive V or W cycles
*/
void solve(int l, int lmax, int npre, int npost, int wnmin, int wnmax) 
{
  int ressmallenough;
  double res2, resi;
  int nwcycle = 0;
  int pr = 0;

  if (pr) printf(".............. entering solve level %d ............\n", l);

  /* on coarsest level */
  if (l == 0) {

    /* solve "exactly" */
    exactsolve(l);
  }

  /* if not on the coarsest level */
  else {

    /* relax */
    relax(l, npre);

    /* transfer from fine to coarse */
    ressmallenough = finetocoarse(l, l-1);

    /* perform W cycles if requested */
    while (1) {
      
      /* solve on coarser level */
      solve(l-1, lmax, npre, npost, wnmin, wnmax);

      /* 0 W cycles if we don't want W cycles at all (the default) */
      if (nwcycle >= wnmax) break;

      /* perform minimal number of cycles without any check */
      if (nwcycle < wnmin) continue;

      /* check errors to find out whether we are done with W cycling */
      res2 = resi = 0;
      restau(l-1, (tolres2)?&res2:0, (tolresi)?&resi:0, 0, 0, 0);
      if (vb==2) {
	printf(". . . . . . .  W-cycle %d  . . . . . . .\n", nwcycle);
        printf("l %d   res2 %e   tolres2 %e\n\n", l-1, res2, tolres2);
      }
      if (0) printf("resi %e   tolresi %e\n", resi, tolresi);

      if (res2 <= tolres2 && resi <= tolresi)
	break;
      
      nwcycle++;
    }

    /* transfer from coarse to fine, perform coarse grid correction */
    coarsetofine(l-1, l);

    /* relax */
    relax(l, npost);
  }

  if (pr) printf("..............  leaving solve level %d ............\n\n",l);
}



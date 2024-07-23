/* exactsolve.c */
/* Bernd Bruegmann 4/97, 3/03 */

#include "bam.h"
#include "multigrid.h"

#define PR 0



/* exactly solve on this level by sufficiently many relaxation sweeps
   applicable only to very small grids
   if you really want a headache, use PETSc here
   -> but now we have the new bam and can use bicgstab etc

   needs clean up
*/
void exactsolve(int lh) 
{
  int i, u, v, w, f;
  tL *gh = setlevelh(lh);
  int nmin = Geti("mg_exact_nmin");
  int nmax = Geti("mg_exact_nmax");
  double normmax = tolres2/10;
  double normres, normres0;
  int pr = PR;
  int pr2 = 0;
  //int o = (problem == VEC4) ? 1 : 0;
  //extern int lvcycle, nvcycle;
  int lvcycle = 0, nvcycle = 0;
  
  //forlevel(l, g) {
  //shrinkrange(g);
  if (gh) {
    normres0 = findres(0, gh);
    if (pr) printf("exactsolve %4d   normres %.15e\n", 0, normres0);
                               if (pr2) prprim("exact bottom: f", fh);
    for (i = 1; i <= nmax; i++) {
      //setrangeall(g);
                               if (pr2&&i<5) prprim("exact bottom: u", uh);
      //shrinkrange(g);
      rb_gauss_seidel(uh, vh, fh);
      if (i >= nmin) {
	normres = findres(0, gh);
#ifdef HAVE_ISNAN
	if (isnan(normres)) 
	  errorexit("exactsolve(): nan");
#endif
	if (pr&&i%nmin==0) 
	  printf("exactsolve %4d   normres %.15e\n", i, normres);
	if (normres < normmax)       /* if error small enough */
	  break;
      }
      if (1 && i < nmin) {
	if (pr) printf("exactsolve %4d   normres %.15e\n", i, findres(0, gh));
      }
    }
    //setrangeall(g);
                               if (pr2) prprim("after exact : u", uh);
                               if (pr2) prprim("after exact : f", fh);
    if (1) {
      int *ib = gh->ibbox;
      //shrinkrange(g);
      normres = findres(0, gh);
      if (vb==2||pr) printf("exactsolve: swept %3d for res %10.3e, V %d, ", 
		     i-1, normres, (ib[1]-2)*(ib[3]-2)*(ib[5]-2));
      if (vb==2||pr) printf("vcycle l=%d, n=%d\n\n", lvcycle, nvcycle);
#ifdef HAVE_ISNAN
      if (isnan(normres)) errorexit("exactsolve(): nan");
#endif
      if (0 && normres >= normmax) {
	prprim("after inexact : u", uh);
	prprim("after inexact : f", fh);
      }
    }
  }
}





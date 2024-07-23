/* inject.c */
/* Bernd Bruegman 12/99, 3/03 */

#include "bam.h"
#include "amr.h"

tCAMR CAMR;

/* check whether parent level exists and is aligned in time */
int parentaligned(tL *level) 
{
  tG *g = level->grid;
  int l = level->l;

  if (l <= g->lmin)
    return 0;

  if (dequal(level->time, g->level[l-1]->time))
    return 1;

  return 0;
}




/* helper function for ghost parent initialization
   if ghost parent is used, make sure it is initialized 
   in particular, variables may be added and/or enabled at any time
*/
void enableparent(tL *lf, tL **plc)
{
  tL *lp = lf->prlocal;
  tL *lc = *plc;
  
  if (lp && lc != lp) {
    if (0) printf("enableparent: lp %p  lc %p  l%d  %d %d\n", 
		    lp, lc, lc->l, lc->nvariables, lp->nvariables);
    //lp->time = lc->time;
    //lp->iteration = lc->iteration;
    if (lc->nvariables != lp->nvariables)
      realloc_levelvariables(lp, lc->nvariables);
    enablesamevars(lc, lp);
    *plc = lp;
  }
}








/****************************************************************************/
/* restrict and prolong using flagrestrict and flagprolong
   respects buffer zone
   uses 4th order polynomial using 4^3 cube
   assumes all cubes exist
*/




/* wrapper for RP operation */
void restrict_prolong_varlist(tG *g, int lcoarse, int lfine, 
			      tVarList *uc, tVarList *uf, int nbuffer)
{
  
  if (0) printf(" do a restrict prolong:\n");
  
  timer_start(0, "restrict_prolong");

  /* restrict prolong between shells and box in a special way */
  if (g->level[lcoarse]->shells) {

    timer_start(0, "restrict_prolong_shells");

    restrict_prolong_shells(g, lfine, lcoarse, uf, uc);
    
    timer_stop(0, "restrict_prolong_shells");

  /* general box to box restrict prolong call */
  } else {

    timer_start(0, "restrict_varlist");

    restrict_varlist(g, lfine, lcoarse, uf, uc,  nbuffer);

    timer_stop(0, "restrict_varlist");

      if (GetvLax("conservative_amr","yes")) {

         timer_start(0, "correct_varlist_camr"); 
      
         CAMR.correct_varlist(g, lfine, lcoarse, uf, uc, nbuffer);     
         
         timer_stop(0, "correct_varlist_camr"); 
       }
    
    timer_start(0, "prolong_varlist");
    
    prolong_varlist(g, lcoarse, lfine, uc, uf,  nbuffer);

    timer_stop(0, "prolong_varlist");
    
  }
  
  timer_stop(0, "restrict_prolong");
}




/* wrapper for applying restriction and prolongation to evolution variables */
void restrict_prolong_evolve(tG *g, int lcoarse, int lfine) 
{ 
  tVarList *u, *up, *upp;

  if (0) printf("restrict_prolong_evolve  %d %d\n", lcoarse, lfine);
  
  /* get list of variables involved in evolution */
  evolve_vlretrieve(&u, &up, &upp);

  if (parentaligned(g->level[lfine])) upp = u;

  /* do it */

  restrict_prolong_varlist(g, lcoarse, lfine, upp, u, 3);

 // restrict_prolong_varlist(g, lcoarse, lfine, upp, u, Geti("amr_nbuffer"));
}




/* fill all levels by R/P as done during evolution */
void restrict_prolong_grid(tG *g, tVarList *u)
{
  int l;
  
  for (l = g->lmax; l > g->lmin; l--) {
     restrict_prolong_varlist(g, l-1, l, u, u, 3);
    //   restrict_prolong_varlist(g, l-1, l, u, u, Geti("amr_nbuffer"));
  }
}





/* apply restriction and prolongation for refinements 
   to be called from fine to coarse
*/
void set_boundary_refinement(tL *level, tVarList *u) 
{ 
  tG *g = level->grid;
  int lmin = level->l;
  int l;

  /* only apply if next coarser level exists and is time aligned */
  if (!parentaligned(level)) return;

  for (l = g->lmax; l > lmin; l--) {
    if (0) printf("set_boundary_refinement c%d f%d\n", l-1, l); 
      restrict_prolong_varlist(g, l-1, l, u, u, 3);
    // restrict_prolong_varlist(g, l-1, l, u, u, Geti("amr_nbuffer"));  
  } 
}




/* set boundary data for refinements so that best available data is used 
   - set boundary of fine by prolongation from coarse
   - set interior of coarse by restriction from fine
   call from finest to coarsest in a stack (not coarsest to finest)
*/
void set_boundary_refinement_cf(tG *g, int lcoarse, int lfine, tVarList *u)
{ 
  if (lfine != lcoarse + 1)
    errorexit("set_boundary_refinement_cf: watch your levels");
  if (lfine > g->lmax || lcoarse < g->lmin) return;
  if (!u || !u->n) return;
  if (0) printf("set_boundary_refinement_cf c%d f%d\n", lcoarse, lfine); 

  restrict_prolong_varlist(g, lcoarse, lfine, u, u, 3);
 //   restrict_prolong_varlist(g, lcoarse, lfine, u, u, Geti("amr_nbuffer"));
}









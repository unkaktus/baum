/* fasmg.c */
/* Bernd Bruegmann, 4/97, 3/03 */

/* Full Approximation Storage MultiGrid */
/* Here we perform the fine to coarse to fine operations. */

#include "bam.h"
#include "multigrid.h"

#define PR 0



/* memory usage is close to optimal in that only one auxiliary field is used
   a single scratch field should always be available anyway
   makes rbGS much simpler
   however, the algorithm would look simpler if auxiliary fields were
   available for the truncation error and the residuum
   ... already added one more field :)
*/

/* optimization:
   - make special versions of all restrictors and prolongators that
     add their result to a given field, subtract from and morf ...
*/




/* prepare coarse grid equation:  
     LH uH = fH
     uH = R uh
     fH = R fh + tauHh,  tauHh = L R uh - R L uh
   control: 
     no sense to make residuum smaller than estimate of
     true truncation error, which for refinement factor 2 is
     approximated by tauHh/3
   return 1 if |dh| < |tauHh|/3 
   return 0 else
*/
int finetocoarse(int lh, int lH) 
{
  double normres, normtau;
  double epsrestau = -1;
  int pr = PR;

  setlevelHh(lH, lh);

  /* operate on fine seperately because this may be expensive */ 
                                 if (pr) prprim("uh", uh);
  lop(uh, vh);
                                 if (pr) prprim("vh = L(uh)", vh);

  /* for all the overlap between fine and coarse */
  mg_restrict(vh, vH, 2);
                                 if (pr) prprim("vH = R(vh)", vH);
  subtract(vh, fh, vh);
                                 if (pr) prprim("dh = vh = vh - fh", vh);
                                 if (pr) normres = l2norm(vh);
                                
                                 if (pr) prprim("fh", fh);
                                 if (pr) prprim("uh", uh);

  // unclear: what is correct range and order? BB 29.5.04
  mg_restrict(uh, uH, 2); 
  setbound(uH);            // should not be needed after restrict 

                                 if (pr) prprim("uH = R(uh)", uH);
  lop(uH, fH);
                                 if (pr) prprim("fH = LR(uh)", fH);
  subtract(fH, vH, fH);
                                 if (pr) prprim("tau = LRu-RLu = fH-vH", fH);

  if (pr) normtau = l2norm(fH);
  if (0 && normres < normtau * epsrestau) 
    return 1;

  mg_restrict(fh, vH, 2);
                                 if (pr) prprim("fH", vH);
  add(vH, fH, fH);
  //setghosts(fH);
                                 if (pr) prprim("fH = R(fh) + tau", fH);
  copyall(uH, wH);
                                 if (pr) prprim("uHold = uH", wH);
  return 0;
}






/* perform coarse grid correction
     uh = uh + P (uH - R uh)
   one of several beautiful instances where doing something harmless
     like setting P R uh = uh spoils everything
*/
void coarsetofine(int lH, int lh)
{
  int pr = PR;

  /* for all the overlap between fine and coarse */
  setlevelHh(lH, lh);

  /* correct fine u */
  /* (!) don't set vH and vh to zero to allow update of ref bound */
  subtract(uH, wH, vH);
  //setbound(vH);
                                     if (pr) prprim("vH = uH - uoldH", vH);
  mg_prolong(vH, vh, 2);
                                     if (pr) prprim("vh = P(vH)", vh);
  add(uh, vh, uh);
  //setbound(uh);
                                     if (pr) prprim("uh = uh + vh", uh);

  if (lh > uh->level->grid->ltop) {
    if (pr) printf("prolongating into refinement boundary, l=%d, ltop=%d\n", 
		   lh, uh->level->grid->ltop);
    mg_prolong(uH, uh, 4);
  }
}





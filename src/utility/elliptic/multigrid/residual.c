/* residual.c */
/* Bernd Bruegmann 4/97, 3/03 */

#include "bam.h"
#include "multigrid.h"

#define PR 0




/* find residual for interior of a given box, uses v as temporary storage */
double findres(int p, tL *gh)
{
  double normres = 0;

  setlevelh(gh->l);
  lop(uh, vh);

  /*
  if (p) prprim("\nfindres uh", uh);
  if (p) prprim("findres fh", fh);
  if (p) prprim("findres luh", vh);
  */
  subtract(vh, fh, vh);
  normres = l2norm(vh);
  
  //if (p) prprim("findres res", vh);
  if (p) printf(" %2dx%2dx%2d  %6.3e\n", 
		gh->ibbox[1]+1, gh->ibbox[3]+1, gh->ibbox[5]+1, normres);
  
  return normres;
}  





/* for diagnostics, compute residuals and truncation errors for one level 
   does not leave fields unchanged:  vh = resh, uH = Ruh, wH = tauH,
*/
void restau(int lh, double *res2, double *resi, double *tau2, double *taui, 
	    int p) 
{
  int n;
  double x;
  int pr = PR;

  setlevelh(lh);
  if (pr) printf(".............. prnormres l%d ............\n", lh);

  /* residuum, only good for printing since not cumulative */
  if (resi || res2) {    
        
    //                         if (pr) prprime("uh", uh);    
    lop(uh, vh);
    //                           if (pr) prprime("vh = L(uh)", vh);    
    subtract(vh, fh, vh);
    //                           if (pr) prprime("dh = vh = vh - fh", vh);
    if (resi) *resi = linorm(vh);
    if (res2) *res2 = l2norm(vh);
  }

  /* truncation error
     to be independent of nesting property, average over boxes */
  if (taui || tau2) {
    if (taui) *taui = 0;
    if (tau2) *tau2 = 0;
    n = 0;
    if (lh > 0) {
      setlevelH(lh-1);

      mg_restrict(vh, vH, 2);
      //if (pr) prprim("vH = R(vh)", vH);    

      //                         if (pr) prprim("uh", uh);
      mg_restrict(uh, uH, 2);

      //                         if (pr) prprim("uH = R(uh)", uH);
      lop(uH, wH);
      //                         if (pr) prprim("wH = LR(uh)", wH);
      subtract(wH, vH, wH);
      //                         if (pr) prprim("tau = LRu-RLu = fH-vH", wH);
      if (taui) {			 
       x = linorm(wH);
      if (*taui < x) *taui = x;
      } 
      if (tau2) {
	*tau2 += l2norm(wH);
	n++; /* unused */
      }
    }

    /* turn tau^H_h into tau^H, the relative local truncation error */
    if (taui) *taui /= 3;
    if (tau2) *tau2 /= 3;

    /* print */
    if (p) {
      int *ibbox = levelof(lh)->ibbox;
      printf("  l%-2d  %3dx%3dx%3d   ",
	     lh, ibbox[1]+1, ibbox[3]+1, ibbox[5]+1);
      if (1) {
	printf(" %.9e %.9e  ", *res2, *resi);
	if (n) printf(" %.9e %.9e", *tau2, *taui);
      } else {
	printf(" %.3e %.3e  ", *res2, *resi);
	if (n) printf(" %.3e %.3e", *tau2, *taui);
      }
    }
  }
  if (p) printf("\n");
}




/* for diagnostics, compute and print residuals per grid for one level 
   does not leave fields unchanged:  vh = resh, uH = Ruh, wH = tauH,
*/
void prnormres(int lh) 
{
  double res2, resi, tau2, taui;

  restau(lh, &res2, &resi, &tau2, &taui, 1);
}




/* for diagnostics, print residuals for all grids on all levels 
   does not leave u and f on l=0..lmax-1 unchanged
*/
void prresiduals(int lmax)
{
  int l;  

  printf("level      size        l2(res)   li(res)     l2(tau)   li(tau)\n");
  for (l = 0; l <= lmax; l++)
    prnormres(l);
  printf("\n");
  fflush(stdout);
}

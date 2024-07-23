/* multigrid.c */
/* Bernd Bruegmann 11/02 */

#include "bam.h"
#include "multigrid.h"



int vb;
tG *stack;
tVarList *uh, *vh, *wh, *fh;
tVarList *uH, *vH, *wH, *fH;

void (*multigrid_lop)(tL *, tVarList *, tVarList *);
void (*multigrid_lgs)(tL *, tVarList *, tVarList *);

/* v = L u */
void lop(tVarList *u, tVarList *v)
{
  multigrid_lop(u->level, v, u);
}

/* v = u + (f-Lu)/Lii */
void lgs(tVarList *u, tVarList *v)
{
  multigrid_lgs(u->level, v, u);
}




/* solve elliptic equation with multigrid method
   this is the standard bam interface routine which is called from outside
*/
int multigrid(tL *level, 
	      tVarList *u, tVarList *f, tVarList *v, tVarList *c, 
	      int itmax, double tol, double *normres,
	      void (*lop)(tL *, tVarList *, tVarList *), 
	      void (*lgs)(tL *, tVarList *, tVarList *))
{
  int l;
  tG *grid = level->grid;
  tVarList *all, *w;

  if (Geti("bampi_nghosts") < 2) 
    errorexit("multigrid: need nghosts >= 2");

  vb = 2*Getv("mg_verbose", "yes") + Getv("mg_verbose", "some");
  if (vb) printf("Calling multigrid solver\n");

  /* register functions */
  multigrid_lop = lop;
  multigrid_lgs = lgs;
  if(Getv("mg_setbound", "set_boundary_elliptic"))
  {
    multigrid_setbound = set_boundary_elliptic;
  }
  else
  {
    printf("multigrid: using function pointer: "
           "void (*multigrid_setbound)(tL *, tVarList *)\n"
           "to set boundaries. "
           "Caller has to set multigrid_setbound to a "
           "user supplied boundary routine!\n");
    printf("multigrid_setbound = %p\n",multigrid_setbound);
    if(multigrid_setbound==NULL)
      errorexit("multigrid_setbound is not set correctly (it is NULL)!");
  }

  /* temporary storage */
  w  = AddDuplicateEnable(u, "_multigrid_w");

  /* create and set global variable lists */
  uH = vlduplicate(uh = u);
  vH = vlduplicate(vh = v);
  wH = vlduplicate(wh = w);
  fH = vlduplicate(fh = f);

  /* create variable list with all the required variables */
  all = vlalloc(level);
  if (c) vlpushvl(all, c);
  vlpushvl(all, u);
  vlpushvl(all, v);
  vlpushvl(all, w);
  vlpushvl(all, f);

  /* create multigrid stack and initialize it */
  makestack(level);
  initializestack(level, all);
  stack = level->grid;
  
  for (l = grid->lmin; l <= grid->lmax; l++) {
    vh->level = grid->level[l];
    vlsetconstant(vh,0);
    wh->level = grid->level[l];
    vlsetconstant(wh,0);
    fh->level = grid->level[l];
    vlsetconstant(fh,0);
   }

  /* solve */
  fullVdriver(level, tol);

  /* provide a one number check */
  if (vb > 1) {
    printf("After multigrid solve\n");
    printvarinfo(level, u->index[0]);
    bampi_vlsynchronize(v);
    set_boundary_elliptic(level, v);
    printvarinfo(level, v->index[0]);
  }

  /* clean up */
  vlfree(uH);
  vlfree(vH);
  vlfree(wH);
  vlfree(fH);
  vlfree(all);
  VLDisableFree(w);
  removestack(level);
  if (vb) printf("Leaving multigrid solver\n");
  return 0;
}






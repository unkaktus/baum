/* makestack.c */
/* Bernd Bruegmann 9/97, 3/03 */

/* Make multigrid stack of coarser grids. 

   The strategy is to call the standard top level grid routine because then
   we do not even have to think about parallelization. Afterwards we
   sort out the parent/child relations.

   The trick is to allow only certain number of grid points in each direction
   which are small multiples of powers of 2 such that all the interlevel
   operations are indeed processor local.

   MUCH cleaner than the old PUGH_makestack. You have to see that one to 
   believe it.
*/

#include "bam.h"
#include "multigrid.h"


void findparentchild(tL *c, tL *f);
void mg_setRPflags(tL *c, tL *f); 






/* check whether coarsening a level is allowed and return new size

   currently one has to add two points for two Robins by hand in parameter
   the new size is therefore new = (old - 2)/2 + 2 = old/2 + 1
   and note that old and old-2 have the same parity 
*/
int coarsegridsize(int nlevels, int *m, int *n, int *o)
{
  int nghosts = Geti("bampi_nghosts");
  int mnew, nnew, onew;
  int i;

  /* check parameter that limits the number of levels allowed */
  i = Geti("mg_levels_nmax");
  if (i != -1 && nlevels > i)
    return 0;

  /* fake it for Cartoon: 3 points in Cartoon direction, 4/2 + 1 = 3 */
  // if (Getv("grid", "1d"))
  //   *m = *n = 4;
  
  /* check whether number of points is suitably even */
  if ((*m-2)%4 || (*n-2)%4 || (*o-2)%4)
    return 0;

  /* compute new size */
  mnew = (*m)/2 + 1;
  nnew = (*n)/2 + 1;
  onew = (*o)/2 + 1;

  /* check parameter that limits the size of global grid */
  i = Geti("mg_size_nmin");
  if (mnew < i || nnew < i || onew < i)
    return 0;

  /* check whether there would be at least 2 points per processor */
  i = Geti("mg_size_nmin_local");
  if ((mnew/Geti("bampi_xsize")) < i || 
      (nnew/Geti("bampi_ysize")) < i ||
      (onew/Geti("bampi_zsize")) < i) {
    printf("Grid cannot be further coarsened: %d %d %d\n", *m, *n, *o);
    return 0;
  }

  /* ok, go for it */
  *m = mnew;
  *n = nnew;
  *o = onew;
  return 1;
}




/* make stack */
void makestack(tL *level)
{
  int l, m, n, o, flag;
  double dx, dy, dz;
  tG *grid = level->grid;
  tL *newlevel;

  /* parameters of the given level */
  l = grid->nlevels;
  m = Geti("nx");
  n = Geti("ny");
  o = Geti("nz");
  dx = level->dx;
  dy = level->dy;
  dz = level->dz;

  /* add coarsened levels as long as allowed */
  while (coarsegridsize(l, &m, &n, &o)) {
    if (vb > 1) prdivider(0);
    if (vb > 1) printf("Adding level %d to multigrid stack\n", -l);

    /* make new level 
       turn off automatic gridsize
    */
    dx *= 2;
    dy *= 2;
    dz *= 2;
    flag = Getv("grid_fournplustwo", "yes");
    if (flag) Sets("grid_fournplustwo", "no");
    newlevel = make_level(dx, dy, dz, m, n, o, vb>1);
    if (flag) Sets("grid_fournplustwo", "yes");

    /* add level to stack */
    push_level(grid, newlevel);
    l++;

    /* find parent/child relationships */
    findparentchild(grid->level[0], grid->level[1]);

    /* set flags */
    mg_setRPflags(grid->level[0], grid->level[1]);
  }

  /* print info */
  if (vb > 1) {
    prdivider(0);
    printgrid(grid);
  }
}




/* initialize stack */
void initializestack(tL *level, tVarList *fine)
{
  tG *grid = level->grid;
  tVarList *coarse = vlduplicate(fine);
  int i, l;

  if (0) printf("initializestack l %d, lmin %d, lmax %d\n",
		level->l, grid->lmin, grid->lmax);

  /* make sure we have storage on every level */
  for (l = grid->lmin; l <= grid->lmax; l++) {
    fine->level = grid->level[l];
    vlenable(fine);
  }

  /* initialize finest to coarse by restriction */
  for (l = grid->lmax; l > grid->lmin; l--) {
    fine->level   = grid->level[l];
    coarse->level = grid->level[l-1]; 

    mg_restrict(fine, coarse, 2);
    
    if (Getv("boundary", "robin")) find_robin_normal(coarse->level);
    /* boundary is untouched, but
       - for u et al we set boundary later
       - for coefficients we do not need boundary values */

    if (Getv("boundary", "excision"))
      set_boundary_flags_excision(coarse->level);
  }

  /* debug */
  if (0) {
    prvare(grid->level[0], "punctures_b");
    prvare(grid->level[1]->prlocal, "punctures_b");
    prvare(grid->level[1], "punctures_b");
    prvare(grid->level[0], "x");
    prvare(grid->level[1]->prlocal, "x");
    prvare(grid->level[1], "x");
  }
  if (0)
    for (l = grid->lmin; l <= grid->lmax; l++)
      prvare(grid->level[l], "punctures_b");
  /*
  if (0) {
    prolong(coarse, fine);
    prvare(fine->level, "punctures_b");
  } */

  /* initialize red-black flag on every level (misnomer) */
  initializeredblack(level);

  /* clean up */
  vlfree(coarse);
  if (vb > 1) prdivider(0);
}




/* remove stack */
void removestack(tL *level)
{
  int i;

  if (vb > 1) {
    prdivider(0);
    printf("Removing stack: level->l = %d\n", level->l);
  }

  /* remove and free levels that were added as multigrid stack */
  while (level->l > 0)
    remove_top_level(level->grid);

  /* the new top level no longer has any parents */
  level->box[0]->pr = 0;

  /* get your own if you want one elsewhere in bam ... (should change) */
  disablevar(level, Ind("mg_redblack"));
}




/* find parent/child relationships of two overlapping levels */
/* assumes Robin ghost at outer boundary
   note that processor boundaries are allowed to shift between levels
   (this makes some non-power-two configurations possible)
*/
void findparentchild(tL *coarse, tL *fine)
{

  /* preliminary: assume just one box */
  fine->box[0]->pr = coarse->box[0];
  
  // FIXME: is there something missing???
}





/* set R and P flags
   we cannot use the standard AMR function setRPflags because
   the multigrid stack implements the Robin ghosts with dangling nodes
   (parents may have incomplete set of children)

   default:
   course, Rflag:    3   3   3   3   2   1   0   0
   fine,   Pflag:   0 0 0 0 0 0 0 0 0 3 2 1

   mg stack:
   course, Rflag:    1   1   1   1   1   1   0
   fine,   Pflag:   1 1 1 1 1 1 1 1 1 1 1 1 0

   unfinished: standard R/P operators are not compatible yet
*/
void mg_setRPflags(tL *lc, tL *lf) 
{
  //double *p = PtrEnable(lf, "flagprolong");
  double *r = PtrEnable(lc, "flagrestrict");
  int i;

  /* preliminary */
  forallpoints(lc, i)
    r[i] =  boundaryflag(lc, i) ? 0 : 1;

}














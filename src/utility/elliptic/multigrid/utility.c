/* utility.c */
/* Bernd Bruegmann 4/97, 3/03 */

#include "bam.h"
#include "multigrid.h"

#define PR 0




/* get level pointer at a given level index */
tL *levelof(int l)
{
  if (l < stack->lmin || l > stack->lmax)
    errorexit("levelof: index out of range");
  
  return stack->level[l];
}

void debuglevelof(int l)
{
  printf("l=%d stack=%p ", l, stack);
  //setlevelH(l);
  printf("level=%p l=%d, uH->n = %d\n", levelof(l), levelof(l)->l, uH->n);
}




/* set levels in global variable lists */
tL *setlevelh(int lh)
{
  tL *level = levelof(lh);
  uh->level = level;
  vh->level = level;
  wh->level = level;
  fh->level = level;
  return level;
}

tL *setlevelH(int lH)
{
  tL *level = levelof(lH);
  uH->level = level;
  vH->level = level;
  wH->level = level;
  fH->level = level;
  return level;
}

void setlevelHh(int lH, int lh)
{
  setlevelH(lH);
  setlevelh(lh);
}



/* prolongate u */
void prolong_u(int lH, int lh)
{
  tL *gH = setlevelH(lH);
  tL *gh = setlevelh(lh);
  int pr = 0;

  if (vb==2 || pr) 
    printf(".............. prolongating u %d -> %d .............\n", lH, lh);

  /* prolongate everywhere */
  if (pr) prprime("prolong_uH", uH);
  if (pr) prprime("prolong_uh", uh);
  
  mg_prolong(uH, uh, -4);
  
  
  if (pr) prprime("prolong_u: uh after PHh", uh);
  setbound(uh);
  if (pr) prprime("prolong_u: uh after setbound", uh);
  //setboundmask(gh, u, 1);
  //                         if (pr) prprim("prolong_u: uh bound set", gh,u);
}









/* set boundary */
void setbound(tVarList *u)
{
  multigrid_setbound(u->level, u);
}




/* set all ghosts: set outer and symmetry boundary, synchronize ghosts */
void setghosts(tVarList *u)
{
  bampi_vlsynchronize(u);
  setbound(u);
}




/* copy: v = u */
void copyall(tVarList *u, tVarList *v)
{
  vlcopy(v, u);
}




/* setconstant: w = c */
void setconst(tVarList *vlu, double c)
{
  int i, j;
  double *u;

  for (j = 0; j < vlu->n; j++) {
    u = VLPtr(vlu, j);

    forallpoints(vlu->level, i) 
      u[i] = c;
  }
}




/* add: w = u + v 
   it does not work to loop over all points, so fix ghosts/boundaries later
*/
void add(tVarList *vlu, tVarList *vlv, tVarList *vlw)
{
  int i, j;
  double *u, *v, *w;

  for (j = 0; j < vlu->n; j++) {
    u = VLPtr(vlu, j);
    v = VLPtr(vlv, j);
    w = VLPtr(vlw, j);

    forallinner(vlu->level, i) 
      w[i] = u[i] + v[i];
  }

  if (1) setghosts(vlw);  // not always needed, could be optimized
}




/* subtract: w = u - v 
   it does not work to loop over all points, so fix ghosts/boundaries later
*/
void subtract(tVarList *vlu, tVarList *vlv, tVarList *vlw)
{
  int i, j;
  double *u, *v, *w;

  for (j = 0; j < vlu->n; j++) {
    u = VLPtr(vlu, j);
    v = VLPtr(vlv, j);
    w = VLPtr(vlw, j);

    forallinner(vlu->level, i) 
      w[i] = u[i] - v[i];
  }

  if (1) setghosts(vlw);  // not always needed, could be optimized
}




/* norm */
double l2norm(tVarList *u)
{
  return bampi_allreduce_allnorm2(u);
}

double linorm(tVarList *u)
{
  return bampi_allreduce_allnormInf(u);
}




/* precompute red-black flag */
int mg_redblack = 0;
void initializeredblack(tL *level)
{
  tG *grid = level->grid;
  double *fp, *xp, *yp, *zp;
  double xmin, ymin, zmin;
  int ccc, i, j, k, l;

  for (l = grid->lmax; l >= 0; l--) {
    level = grid->level[l];
    fp = PtrEnable(level, "mg_redblack");
    xp = Ptr(level, "x");
    yp = Ptr(level, "y");
    zp = Ptr(level, "z");

    /* this is important: we use the origin of the GLOBAL bounding box 
       to make the result independent of parallelization
    */
    xmin = level->bbox[0];
    ymin = level->bbox[2];
    zmin = level->bbox[4];

    forallpoints(level, ccc) {
      i = (xp[ccc] - xmin)/level->dx + 0.5;
      j = (yp[ccc] - ymin)/level->dy + 0.5;
      k = (zp[ccc] - zmin)/level->dz + 0.5;
      fp[ccc] = ((i + j + k) % 2) ? 0 : 1;
    }

    if (0) prvar01(level, "mg_redblack");
  }
}


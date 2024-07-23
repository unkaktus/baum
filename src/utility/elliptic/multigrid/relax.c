/* relax.c */
/* Bernd Bruegmann 4/97, 3/03 */

#include "bam.h"
#include "multigrid.h"

#define PR 0


void rb2_gauss_seidel(tVarList *u, tVarList *v, tVarList *f);
void uni_gauss_seidel(tVarList *u, tVarList *v, tVarList *f);
int rb_decouple, bam_relax;

int solve_nmin, solve_nmax;
double solve_normmax;


/* perform a number of relaxation sweeps */
void relax(int l, int nsweeps)
{
  int i;
  tL *g = setlevelh(l);
  int pr = 0;

  if (g) {
    //forlevel(l, g) {
    //shrinkrange(g);
    if (pr) printf("before %d sweeps\n", nsweeps);
    if (pr&&l>1) findres(1, g);

    for (i = 1; i <= nsweeps; i++) {
      rb_gauss_seidel(uh, vh, fh);
      if (0) printf("relax %4d   normres %10.3e\n", i, findres(0, g));
      if (0)       /* if error small enough, or convergence too slow, or ... */
	break;
    }

    if (pr) printf("after %d sweeps\n", nsweeps);
    if (pr&&l>1) findres(0, g);
    //setrangeall(g);
                                  if (pr) prprim("after relax", uh);
  }
}





/* red black non-linear Gauss Seidel relaxation
   note that stencils with mixed derivatives overlap, but 
     stencil separation is not strictly required if loss of symmetry allowed,
     tried red, black, green, yellow, ..., didn't really matter
   v is passed in case an intermediate layer is needed
*/
void rb_gauss_seidel(tVarList *u, tVarList *v, tVarList *f)
{
  int pr = 0;
  int i, j, k, ijk, nc;
  //double lu[MGCMAX], lii[MGCMAX];
  //int o = octoffset + bitoffset + ((cartoon)?oxdo+oydo:0);

  /* may be we don't want red-black at all? */
  if (1 /* bam_relax == 1 */) {
    uni_gauss_seidel(u, v, f);
    return;
  }

#if 0
  /* may be we want to decouple? will be slower */
  if (rb_decouple) {
    rb2_gauss_seidel(g, u, v, f);
    return;
  }

  /* boundary */
  //if (1||!cartoon) setboundrobin(g, u);

  /* for odd points */
  //setrangePUGH(g, 1);
  forrange(g) {
    if ((i+j+k+o)%2) {
      stencil(g, u, ijk, lu, lii);
      for (nc = 0; nc < g->ncomps; nc++)
	g->store[u+nc][ijk] += (g->store[f+nc][ijk] - lu[nc])/lii[nc];
    }
  }
  PUGHsyncvec(g, u);
  if (masks || (cartoon && !cartoon_stencil)) setboundrobin(g, u);
                               if (pr&&g->l==1) prprim("rbGS after odd", g, u);


  /* for even points */
  setrangePUGH(g, 1);
  forrange(g) {
    if ((i+j+k+o+1)%2) {
      stencil(g, u, ijk, lu, lii);
      for (nc = 0; nc < g->ncomps; nc++)
	g->store[u+nc][ijk] += (g->store[f+nc][ijk] - lu[nc])/lii[nc];
    }
  }
  PUGHsyncvec(g, u);
  setboundrobin(g, u);
                              if (pr&&g->l==1) prprim("rbGS after even", g, u);
#endif
}




/* red black non-linear Gauss Seidel relaxation
   decouple by using two levels
*/
void rb2_gauss_seidel(tVarList *u, tVarList *v, tVarList *f)
{
#if 0
  int pr = 0;
  int i, j, k, ijk, nc;
  double lu[MGCMAX], lii[MGCMAX];
  int o = octoffset + bitoffset + ((cartoon)?oxdo+oydo:0);

  /* boundary */
  setboundrobin(g, u);

  /* make copy */
  setrangeall(g);
  copy(g, u, v);
  shrinkrange(g);

  /* for odd points */
  setrangePUGH(g, 1);
  forrange(g) {
    if ((i+j+k+o)%2) {
      stencil(g, u, ijk, lu, lii);
      for (nc = 0; nc < g->ncomps; nc++)
	g->store[v+nc][ijk] += (g->store[f+nc][ijk] - lu[nc])/lii[nc];
    }
  }
  PUGHsyncvec(g, v);
  if (masks || (cartoon && !cartoon_stencil)) setboundrobin(g, v);
                               if (pr&&g->l==1) prprim("rbGS after odd", g, v);

  /* for even points */
  setrangePUGH(g, 1);
  forrange(g) {
    if ((i+j+k+o+1)%2) {
      stencil(g, v, ijk, lu, lii);
      for (nc = 0; nc < g->ncomps; nc++)
	g->store[u+nc][ijk] += (g->store[f+nc][ijk] - lu[nc])/lii[nc];
    } else {
      for (nc = 0; nc < g->ncomps; nc++)
	g->store[u+nc][ijk] = g->store[v+nc][ijk];
    }
  }
  PUGHsyncvec(g, u);

  /* boundary */
  setboundrobin(g, u);
                              if (pr&&g->l==1) prprim("rbGS after even", g, u);

  /* may have to worry about cleanup */
  setrangeall(g);
  subtract(g, v, v, v);
  shrinkrange(g);
#endif
}





/* uniform non-linear Gauss Seidel relaxation
   decoupled by using two levels
*/
void uni_gauss_seidel(tVarList *u, tVarList *v, tVarList *f)
{
  int pr = 0;
  int l = u->level->l;

  copyall(u, v);
  lgs(v, u);
  setconst(v, 0);

  /* boundary */
  //setboundrobin(g, u);

  /* make copy */
  //setrangeall(g);
  //copy(u, v);
  //shrinkrange(g);

  /* for all points */
  //setrangePUGH(g, 1);

  /* u = u + (f - Lu) / (dLu/du) */
  //forrange(g) {
  //  stencil(g, v, ijk, lu, lii);
  //  for (nc = 0; nc < g->ncomps; nc++)
  //    g->store[u+nc][ijk] += (g->store[f+nc][ijk] - lu[nc])/lii[nc];
  //}
  //PUGHsyncvec(g, u);
  //                          if (pr&&l==0) prprim("rbGS after uni", v);

  /* boundary */
  //setboundrobin(g, u);
  //                          if (pr&&g->l==1) prprim("rbGS after even", u);

  /* may have to worry about cleanup */
  //setrangeall(g);
  //subtract(v, v, v);
  //shrinkrange(g);

  //#endif
}






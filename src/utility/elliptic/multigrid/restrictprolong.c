/* makestack.c */
/* Bernd Bruegmann 
   mth 05/2012
*/


#include "bam.h"
#include "multigrid.h"

#define PR 0






/****************************************************************************/
/* restrict and prolong with simple trilinear operators */

/* restrict single variable */
void mg_restrict_var(tL *lf, int nvarf, tL *lc, int nvarc, int doboundary)
{
  
  
  // here is a quick hack, all this has to be cleaned up
  double *flagrestrict = Ptr(lc, "flagrestrict");
  double *f = lf->v[nvarf];
  double *c = lc->v[nvarc];
  tB *fbox;  
  int fijk, fdi, fdj, fdk;
  
  if (PR) printf("mg_restrict_var:  %s  %p %p %p\n",VarName(nvarf),flagrestrict,f,c );

  if (!doboundary) errorexit("BOX: implement doboundary = 0");

  /* for all nodes on coarse level that need restriction */
  forallpoints_ijk(lc) if (flagrestrict[ijk]) {
    fbox = box_containing_ijk(lf, box, i, j, k);

    fdi = fbox->di;
    fdj = fbox->dj;
    fdk = fbox->dk;
     
    fijk = box_child_ijk(fbox, box, i, j, k);
    
    c[ijk] = 0.125 * (f[fijk] + f[fijk+fdi+fdj+fdk] +
        f[fijk+fdi] + f[fijk+fdj] + f[fijk+fdk] +
        f[fijk+fdi+fdj] + f[fijk+fdi+fdk] +f[fijk+fdj+fdk]);
  } endfor_ijk;
}

/* restrict */
void mg_restrict(tVarList *f, tVarList *c, int order)
{
  if (PR) printf("mg_restrict:  %d -> %d\n",f->level->l, c->level->l);
  
  tG *g = f->level->grid;
  int lfine = f->level->l;
  int lcoarse = c->level->l;
  tVarList *t;
  int i;

  /* check */
  if (abs(lfine-lcoarse) != 1)
    errorexit("restrict: watch your levels");
  if (lfine < lcoarse) {t = f; f = c; c = t;}
  if (f->n != c->n)
    errorexit("restrict: variable lists of different length");

  /* call appropriate restriction operator */
  if (lfine <= g->ltop || order == 2) {
    enableparent(f->level, &c->level);
    for (i = 0; i < f->n; i++)
      mg_restrict_var(f->level, f->index[i], c->level, c->index[i], 1);
    bampi_syncparent_send123(f->level, c);
  } else {
    errorexit("not implemented");
  }
}




/***************************************************************************/

/* prolong single variable */
void mg_prolong_var(tL *lc, int nvarc, tL *lf, int nvarf)
{

  // here is a quick hack, all this has to be cleaned up
  double *f = lf->v[nvarf];
  double *c = lc->v[nvarc];
  tB *cbox;  
  int ci, cj, ck, cijk;
  int cdi, cdj, cdk;
  double kx0, ky0, kz0, kx1, ky1, kz1;
  const int mgstack = (lc->l < lc->grid->ltop);

  /* for all nodes on fine level that are not in the mgstack/robin boundary */
  forallpoints_ijk(lf) if (!mgstack || mgstack && !boundary[ijk]) {
    cbox = box->pr;
    cdi = cbox->di;
    cdj = cbox->dj;
    cdk = cbox->dk;
    box_point_indices(cbox, box, &ci, &cj, &ck, i, j, k);
    cijk = ijkofbox(cbox, ci, cj, ck);
    kx1 = (xofbox(box, i) - xofbox(cbox, ci))/lc->dx;
    ky1 = (yofbox(box, j) - yofbox(cbox, cj))/lc->dy;
    kz1 = (zofbox(box, k) - zofbox(cbox, ck))/lc->dz;
    kx0 = 1-kx1;
    ky0 = 1-ky1;
    kz0 = 1-kz1;
    f[ijk] = kx0*ky0*kz0 * c[cijk] +
        kx1*ky0*kz0 * c[cijk+cdi] +
        kx0*ky1*kz0 * c[cijk+cdj] +
        kx1*ky1*kz0 * c[cijk+cdi+cdj] +
        kx0*ky0*kz1 * c[cijk+cdk] +
        kx1*ky0*kz1 * c[cijk+cdi+cdk] +
        kx0*ky1*kz1 * c[cijk+cdj+cdk] +
        kx1*ky1*kz1 * c[cijk+cdi+cdj+cdk];
  } endfor_ijk;
}

/* higher order prolong, this is ONLY an adapted version form 9.06 
   normally this has to be done better, however we only need 4th order
*/
void mg_prolong_varlist4(tG *g, int lcoarse, int lfine, 
                         tVarList *uc, tVarList *uf, int nbuffer)
{
  tL *lf = g->level[lfine];
  tL *lc = g->level[lcoarse];
  double *flagprolong = Ptr(lf, "flagprolong");
  
  /* check */
  if (lfine != lcoarse + 1)
    errorexit("prolong: watch your levels");
  if (!uc || !uc->n || !uf || !uf->n) return;

  /* if ghost parent is used, make sure it is initialized */
  enableparent(lf, &lc);

  /* synchronize parent */
  bampi_syncparent_recv0123(lf, uc);
  
  /* for all nodes on fine level that need prolongation */
  forallpoints_ijk(lf) {
    if (boundaryflag(lf, ijk) == SYMBOUND) continue;
    if (boundaryflag(lf, ijk) == GHOBOUND && nbuffer < 0) continue;
    if (!(nbuffer >= 3 && flagprolong[ijk]  > 0 ||  // 3 or more buffer nodes
          nbuffer == 1 && flagprolong[ijk] == 1 ||  // 1 buffer node 
          nbuffer < 0))                             // everywhere
      continue;
    
    /* interpolate from coarse to fine level */
    interpolate_TriN_varlist(lc, lf, uc, uf, box->pr, box, i, j, k, 4);

  } endfor_ijk;
 
  /* synchronize fine level */
  uf->level = lf;
  bampi_vlsynchronize(uf);

  /* some boundary conditions need to be applied here */
  set_boundary_symmetry(lf, uf);
}

/* prolong */
void mg_prolong(tVarList *c, tVarList *f, int order)
{
  tG *g = f->level->grid;
  int lfine = f->level->l;
  int lcoarse = c->level->l;
  tVarList *t;
  int i;

  if (abs(lfine-lcoarse) != 1)
    errorexit("prolong: watch your levels");
  if (lfine < lcoarse) {t = f; f = c; c = t;}
  if (f->n != c->n)
    errorexit("restrict: variable lists of different length");

  if (lfine <= g->ltop || order == 2) {
    enableparent(f->level, &c->level);
    bampi_syncparent_recv0123(f->level, c);
    for (i = 0; i < f->n; i++) 
      mg_prolong_var(c->level, c->index[i], f->level, f->index[i]);
    bampi_vlsynchronize(f);
    set_boundary_symmetry(f->level, f);
  } else {
    if (order == 4) 
      mg_prolong_varlist4(c->level->grid, lcoarse, lfine, c, f, 1); 
    else if (order == -4) 
      mg_prolong_varlist4(c->level->grid, lcoarse, lfine, c, f, -1); 
    else
      errorexit("not implemented");
  }

}














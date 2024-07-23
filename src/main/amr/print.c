/* print.c */
/* Bernd Bruegmann, 12/99, 3/03 */

#include "bam.h"
#include "amr.h"


/* print information */

void printlocalgrid(tG *g)
{
  tL *level;
  int ii,jj;
  int m1,n1,o1,m2,n2,o2;

  printf("local boxes:\n");
  for (ii = g->lmin; ii <= g->lmax; ii++) {
    level = g->level[ii];  

    forallboxes(level) {
      printf("l=%d: box=%d   [%7.3f,%7.3f]x[%7.3f,%7.3f]x[%7.3f,%7.3f]  %dx%dx%d = %d\n", 
          level->l, nbox, 
          box->x0, box->x0+(box->m-1)*box->dx, 
          box->y0, box->y0+(box->n-1)*box->dy, 
          box->z0, box->z0+(box->o-1)*box->dz,
          box->m,box->n,box->o,
          box->m*box->n*box->o);
    } endforboxes;
  }
  
  printf("local inner boxes:\n");
  for (ii = g->lmin; ii <= g->lmax; ii++) {
      level = g->level[ii];
      forallboxes(level) {
        m1=-1;
        forinnerpoints_boxijk(level,box) {
          if (m1==-1) {
              m1 = i;
              n1 = j;
              o1 = k;
          }
          m2 = i;
          n2 = j;
          o2 = k;
        } endfor_ijk;
        printf("l=%d: box=%d   [%7.3f,%7.3f]x[%7.3f,%7.3f]x[%7.3f,%7.3f]  %dx%dx%d = %d\n", 
          level->l, nbox, 
          box->x0+m1*box->dx,box->x0+m2*box->dx,
          box->y0+n1*box->dy,box->y0+n2*box->dy,
          box->z0+o1*box->dz,box->z0+o2*box->dz,
          m2-m1+1, n2-n1+1, o2-o1+1,
          (m2-m1+1)*(n2-n1+1)*(o2-o1+1));
      } endforboxes;
  }
}



void printgrid(tG *g) 
{
  int i;
  tL *level;

  printf("grid %p: nlevels %d, lmin %d, lmax %d, nvariables %d\n",
	 g, g->nlevels, g->lmin, g->lmax, g->nvariables);
  for (i = g->lmin; i <= g->lmax; i++) {
    level = g->level[i];

    forallboxes(level) {
      printf("p=%dx%dx%d, ",
	     box->com->sizexyz[0], box->com->sizexyz[1], box->com->sizexyz[2]);
      printbbox(level, box->bbox, box->ibbox);
    } endforboxes;
  }
  
  printlocalgrid(g);
}




void printlevel(tL *l) 
{
  printf("level %p: grid %p, l %d, nboxes %d, box %p, v %p\n", 
	 l, l->grid, l->l, l->nboxes, l->box, l->v);
}







void printbbox(tL *level, double *bbox, int *ibbox)
{
  if (dless(bbox[1], bbox[0]) ||
      dless(bbox[3], bbox[2]) ||
      dless(bbox[5], bbox[4])) {
    printf("l=%d, h=%.3f,%.3f,%.3f: empty\n",
	   level->l, level->dx, level->dy, level->dz);
    return;
  }

  printf("l=%d, h=%.3f,%.3f,%.3f: [%7.3f,%7.3f]x[%7.3f,%7.3f]x[%7.3f,%7.3f] ", 
	 level->l, level->dx, level->dy, level->dz,
	 bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
  if (ibbox) {
    if (!(ibbox[0] || ibbox[2] || ibbox[4])) 
      printf(" %dx%dx%d = %d", 
	     ibbox[1]+1, ibbox[3]+1, ibbox[5]+1, 
	     (ibbox[1]+1)*(ibbox[3]+1)*(ibbox[5]+1));
    else
      printf(" [%d,%d]x[%d,%d]x[%d,%d] = %d", 
	     ibbox[0], ibbox[1], ibbox[2], ibbox[3], ibbox[4], ibbox[5], 
	     (ibbox[1]-ibbox[0]+1)*(ibbox[3]-ibbox[2]+1)
	     *(ibbox[5]-ibbox[4]+1));
  }
  printf("\n");
}







/* print list of all variables */
void printvariables(tL *level)
{
  int n = level->nvariables;
  int i, j;
  char *s;

  for (i = j = 0; i < n; i++) {
    if (level->v[i]) {
      s = "ON ";
      j++;
    } else
      s = "off";
    printf("level->v[%3d] = %s   %s\n", i, s, VarName(i));
  }
  printenabled(level);
}




/* print norm and value at origin for one variable */
void printvarinfo(tL *level, int i)
{
  double *norm, origin;
  tVarList *vl;

  if (!level) return;
  if (!level->v[i]) return;

  vl = vlalloc(level);
  vlpush(vl, i);
  bampi_allreduce_norm(vl, &norm);
  origin = interpolate_xyz_scalar(level, 0, 0, 0, i,  2,LAGRANGE);

  printf("%s:  norm %23.16e,  at origin %23.16e\n",
	 VarName(i), norm[0], origin);
    
  free(norm);
  vlfree(vl);
}













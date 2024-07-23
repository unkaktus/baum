/* interpolate_global.c */
/* Bernd Bruegmann, Jose Gonzalez 4/2006 */

/* interface to interpolation routines that use the best available
   information on an AMR grid
*/


#include "bam.h"
#include "interpolate.h"



/* find and fill interpolation cube in global grid 
   - criterion ...
*/

void interpolate_xyz_grid(tG *g, double x, double y, double z, 
                          int nv, int *iv, double *vinterp,
                          int order, int scheme)
{
  int i, j, k, i0, j0, k0, i1, j1, k1;
  int icorner, iowner;
  int a = order;
  int b = a/2;
  int flag, l;
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
    
  for (l = g->lmin; l <= g->lmax; l++) {
    flag = 0;
    forallboxes(g->level[l]) {
      if (!xyzinsidebbox(box->bbox, x, y, z)) continue;

      i = floor((x - box->com->bbox[0] + dequaleps)/box->dx);
      i0 = i - b + 1;
      i1 = i + b;
      if (i0 < box->com->ibbox[0] || i1 > box->com->ibbox[1]) continue;
      
      j = floor((y - box->com->bbox[2] + dequaleps)/box->dy);
      j0 = j - b + 1;
      j1 = j + b;
      if (j0 < box->com->ibbox[2] || j1 > box->com->ibbox[3]) continue;
      
      k = floor((z - box->com->bbox[4] + dequaleps)/box->dz);
      k0 = k - b + 1;
      k1 = k + b;
      if (k0 < box->com->ibbox[4] || k1 > box->com->ibbox[5]) continue;

      iowner = ijkofbox(box, i, j, k);
      if (boundaryflag(box->level, iowner)) continue;
      //if (flagprolong[iowner] && flagprolong[iowner] < nnn) continue;
  
      /*
      printf("found %7.3f %7.3f %7.3f in ", x, y, z);
      printbbox(box->level, box->bbox, 0);
      */

      /* fill index cube */
      icorner = ijkofbox(box, i0, j0, k0);
      for (k = 0; k < a; k++)
      for (j = 0; j < a; j++)
      for (i = 0; i < a; i++)
	uui[i][j][k] = icorner + i*box->di + j*box->dj + k*box->dk;

      /* for all variables, get data and interpolate */
      for (i = 0; i < nv; i++) {
	interpolate_filldatacube(box->level, iv[i], uui,uuu,order);
	vinterp[i] =
	  interpolate_tripolynomial(x, y, z,
	      xofbox(box,i0), yofbox(box,j0), zofbox(box,k0),
	      box->dx, box->dy, box->dz,
              uuu,order,scheme);
      }
      flag = 1;

    } endforboxes;
    if (flag) break;
  }
}








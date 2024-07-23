/* matter_fluxes_grmhd.c
   AN 2023 */

#include "bam.h"
#include "matter.h"

#include <stdarg.h>

#define PR 0

/* not finite fluxes are set to zero.
   do you want to stop if it happens ? */
#define StopFlxNOTFinite 0

/*stupid hack to make camr work with rk4*/
static int rhs_camr_counter = 1;

// main function which splits the 3d problem into several 1d problems.
// compute 1d fluxes and add them to rhs
void matter_fluxes_grmhd(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars)
{
  bampi_openmp_start

      tL *level = ucur->level;
  int i, j, k, v;
  int imax, jmax, kmax, ijkmax;

  double **fluxs1d = (double **)malloc(MATTER.NVf * sizeof(double *));

  // do for all boxes and all points
  forallboxes(level)
  {

    // save one point more in each direction in compare to normal-variables
    imax = box->m - 1;
    jmax = box->n - 1;
    kmax = box->o - 1;
    ijkmax = (box->m > box->n) ? box->m : box->n;
    ijkmax = (ijkmax > box->o) ? ijkmax : box->o;
    ijkmax += 1;

    // alloc memory
    for (v = 0; v < MATTER.NVf; v++)
      fluxs1d[v] = (double *)malloc((ijkmax + 1) * sizeof(double));

    // compute fluxes in each direction seperatly
    // and add to rhs

    // x-direction
    bampi_openmp_loop_collapse2 for (k = 0; k <= kmax; k++)
    {
      for (j = 0; j <= jmax; j++)
      {

        compute_flux1d_grmhd(1, box, j, k, unew, upre, c, ucur, primVars, otherVars, fluxs1d);
        add_flx1d_rhs(1, box, j, k, unew, upre, c, ucur, primVars, fluxs1d);
      }
    }

    // y-direction
    bampi_openmp_loop_collapse2 for (k = 0; k <= kmax; k++)
    {
      for (i = 0; i <= imax; i++)
      {

        compute_flux1d_grmhd(2, box, i, k, unew, upre, c, ucur, primVars, otherVars, fluxs1d);
        add_flx1d_rhs(2, box, i, k, unew, upre, c, ucur, primVars, fluxs1d);
      }
    }

    // z-direction
    bampi_openmp_loop_collapse2 for (j = 0; j <= jmax; j++)
    {
      for (i = 0; i <= imax; i++)
      {

        compute_flux1d_grmhd(3, box, i, j, unew, upre, c, ucur, primVars, otherVars, fluxs1d);
        add_flx1d_rhs(3, box, i, j, unew, upre, c, ucur, primVars, fluxs1d);
      }
    }

    for (v = 0; v < MATTER.NVf; v++)
      free(fluxs1d[v]);
  }
  endforboxes;

  free(fluxs1d);

  bampi_openmp_stop
}

int dir_perm_full(int direction, char d1, char d2)
{
  /* define: xx = 0, xy = 1, xz = 2, yx = 3, yy = 4, yz = 5, zx = 6, zy = 7, zz = 8 */
  if (direction == 1)
  {
    if (d1 == 'x')
    {
      if (d2 == 'x')
        return 0; // xx
      if (d2 == 'y')
        return 1; // xy
      if (d2 == 'z')
        return 2; // xz
    }
    else if (d1 == 'y')
    {
      if (d2 == 'x')
        return 3; // yx
      if (d2 == 'y')
        return 4; // yy
      if (d2 == 'z')
        return 5; // yz
    }
    else
    {
      if (d2 == 'x')
        return 6; // zx
      if (d2 == 'y')
        return 7; // zy
      if (d2 == 'z')
        return 8; // zz
    }
  }
  else if (direction == 2)
  {
    // x->y, y->z, z->x
    if (d1 == 'x')
    {
      if (d2 == 'x')
        return 4; // yy
      if (d2 == 'y')
        return 5; // yz
      if (d2 == 'z')
        return 3; // yx
    }
    else if (d1 == 'y')
    {
      if (d2 == 'x')
        return 7; // zy
      if (d2 == 'y')
        return 8; // zz
      if (d2 == 'z')
        return 6; // zx
    }
    else
    {
      if (d2 == 'x')
        return 1; // xy
      if (d2 == 'y')
        return 2; // xz
      if (d2 == 'z')
        return 0; // xx
    }
  }
  else
  {
    // x->z, y->x, z->y
    if (d1 == 'x')
    {
      if (d2 == 'x')
        return 8; // zz
      if (d2 == 'y')
        return 6; // zx
      if (d2 == 'z')
        return 7; // zy
    }
    else if (d1 == 'y')
    {
      if (d2 == 'x')
        return 2; // xz
      if (d2 == 'y')
        return 0; // xx
      if (d2 == 'z')
        return 1; // xy
    }
    else
    {
      if (d2 == 'x')
        return 5; // yz
      if (d2 == 'y')
        return 3; // yx
      if (d2 == 'z')
        return 4; // yy
    }
  }
}

/* build 1d arrays for metric, primitives, and conservatives
   extrapolate into ghosts
   call the 1d flux routine pointed by hrsc_flx1d */
void compute_flux1d_grmhd(int direction, tB *box, int dir1, int dir2,
                          tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars,
                          double **fluxs1d)
{
  tL *level = ucur->level;

  int l, m, n, p, ijk;
  double detginv;

  int Nghost = MATTER.HRSC_NGHOST;
  int Nmax = 0;
  int NVg = 10;
  int Nfields = NVg + MATTER.NVw + MATTER.NVq + MATTER.NVo + 1;
  int index_w = NVg;
  int index_q = NVg + MATTER.NVw;
  int index_o = NVg + MATTER.NVw + MATTER.NVq;

  // attention: direction = {1,2,3} NOT {0,1,2}
  int N = (direction == 1) * box->m + (direction == 2) * box->n + (direction == 3) * box->o;
  int Ntot = N + 2 * Nghost;
  int off = Nghost - 1;

  /* contains all informations */
  double **buffer = (double **)malloc(Nfields * sizeof(double *));
  for (l = 0; l < Nfields; l++)
    buffer[l] = (double *)malloc((Ntot) * sizeof(double));

  /* pointer to information ... for convenience */
  double **g1d = (double **)malloc(NVg * sizeof(double *));
  double **w1d = (double **)malloc(MATTER.NVw * sizeof(double *));
  double **q1d = (double **)malloc(MATTER.NVq * sizeof(double *));
  double **o1d = (double **)malloc(MATTER.NVo * sizeof(double *));

  // set the 1d storages to buffer and offset such that var[1] -> buffer[Nghost]
  for (l = 0; l < NVg; l++)
    g1d[l] = buffer[l] + off;
  for (m = 0; m < MATTER.NVw; m++)
    w1d[m] = buffer[index_w + m] + off;
  for (n = 0; n < MATTER.NVq; n++)
    q1d[n] = buffer[index_q + n] + off;
  for (p = 0; p < MATTER.NVo; p++)
    o1d[p] = buffer[index_o + p] + off;

  /* mask information */
  double *mask1d = buffer[Nfields - 1] + off;

  // copy physical data from 3d to 1d arrays ... [1..N]
  for (l = 0; l < N; l++)
  {
    if (direction == 1)
    {
      ijk = box->noffset + l + box->m * dir1 + box->m * box->n * dir2;
      // if (l == 0) printf(" %d \n", box->noffset);
    }
    else if (direction == 2)
      ijk = box->noffset + dir1 + box->m * l + box->m * box->n * dir2;
    else
      ijk = box->noffset + dir1 + box->m * dir2 + box->m * box->n * l;

    // mask
    if (MATTER.USEMASK)
      mask1d[l + 1] = level->v[MATTER.INDX_VAR_mask][ijk];

    // alpha + beta from evolved vars
    g1d[0][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_alp)[ijk];
    g1d[1][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_bet + dir_perm(direction, 'x', ' '))[ijk];
    g1d[2][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_bet + dir_perm(direction, 'y', ' '))[ijk];
    g1d[3][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_bet + dir_perm(direction, 'z', ' '))[ijk];

    // adm metric, which should be computed by the metric evolution project
    g1d[4][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'x', 'x')][ijk];
    g1d[5][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'x', 'y')][ijk];
    g1d[6][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'x', 'z')][ijk];
    g1d[7][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'y', 'y')][ijk];
    g1d[8][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'y', 'z')][ijk];
    g1d[9][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'z', 'z')][ijk];

    // primitives
    for (m = 0; m < MATTER.NVw; m++)
    {
      if (VarNComponents(primVars->index[m]) == 3)
      {
        w1d[m][l + 1] = vldataptr(primVars, m + dir_perm(direction, 'x', ' '))[ijk];
        w1d[m + 1][l + 1] = vldataptr(primVars, m + dir_perm(direction, 'y', ' '))[ijk];
        w1d[m + 2][l + 1] = vldataptr(primVars, m + dir_perm(direction, 'z', ' '))[ijk];
        m += 2;
      }
      else
      {
        w1d[m][l + 1] = vldataptr(primVars, m)[ijk];
      }
    }

    // conservatives
    for (m = 0; m < MATTER.NVq; m++)
    {
      if (VarNComponents(ucur->index[MATTER.INDX_VAR_q + m]) == 3)
      {
        q1d[m][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_q + m + dir_perm(direction, 'x', ' '))[ijk];
        q1d[m + 1][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_q + m + dir_perm(direction, 'y', ' '))[ijk];
        q1d[m + 2][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_q + m + dir_perm(direction, 'z', ' '))[ijk];
        m += 2;
      }
      else
      {
        q1d[m][l + 1] = vldataptr(ucur, MATTER.INDX_VAR_q + m)[ijk];
      }
    }

    // other
    for (m = 0; m < MATTER.NVo; m++)
    {
      if (VarNComponents(otherVars->index[m]) == 3)
      {
        o1d[m][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'x', ' '))[ijk];
        o1d[m + 1][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'y', ' '))[ijk];
        o1d[m + 2][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'z', ' '))[ijk];
        m += 2;
      }
      else if (VarNComponents(otherVars->index[m]) == 9)
      {
        o1d[m][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'x', 'x'))[ijk];
        o1d[m + 1][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'x', 'y'))[ijk];
        o1d[m + 2][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'x', 'z'))[ijk];
        o1d[m + 3][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'y', 'x'))[ijk];
        o1d[m + 3][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'y', 'y'))[ijk];
        o1d[m + 4][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'y', 'z'))[ijk];
        o1d[m + 5][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'z', 'x'))[ijk];
        o1d[m + 5][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'z', 'y'))[ijk];
        o1d[m + 5][l + 1] = vldataptr(otherVars, m + dir_perm_full(direction, 'z', 'z'))[ijk];
        m += 8;
      }
      else
      {
        o1d[m][l + 1] = vldataptr(otherVars, m)[ijk];
      }
    }
  }

  // extrapolate all the data to the ghost points
  for (l = -Nghost + 1; l <= N + Nghost; l++)
  {
    if (l == 1)
      l = N + 1;
    for (m = 0; m < Nfields; m++)
      extrapolate(buffer[m] + off, l, 1, N);
  }

  // finally call the flx 1d routine
  MATTER.hrsc_flx1d(g1d, w1d, q1d, o1d, mask1d, N, Nghost, direction, fluxs1d);

  /* free all fields */
  free(g1d);
  free(w1d);
  free(q1d);
  free(o1d);
  for (l = 0; l < Nfields; l++)
    free(buffer[l]);
  free(buffer);
}

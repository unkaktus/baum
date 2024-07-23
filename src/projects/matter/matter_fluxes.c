/* matter_fluxes.c
   sbernuz, mth 03/2012 */

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
void matter_fluxes(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars)
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

        compute_flux1d(1, box, j, k, unew, upre, c, ucur, primVars, otherVars, fluxs1d);
        add_flx1d_rhs(1, box, j, k, unew, upre, c, ucur, primVars, fluxs1d);
      }
    }

    // y-direction
    bampi_openmp_loop_collapse2 for (k = 0; k <= kmax; k++)
    {
      for (i = 0; i <= imax; i++)
      {

        compute_flux1d(2, box, i, k, unew, upre, c, ucur, primVars, otherVars, fluxs1d);
        add_flx1d_rhs(2, box, i, k, unew, upre, c, ucur, primVars, fluxs1d);
      }
    }

    // z-direction
    bampi_openmp_loop_collapse2 for (j = 0; j <= jmax; j++)
    {
      for (i = 0; i <= imax; i++)
      {

        compute_flux1d(3, box, i, j, unew, upre, c, ucur, primVars, otherVars, fluxs1d);
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

/* we want to compute 1D stuff -> tensors components have to be
   permuted, that we can assume any direction as the x-direction
   FIXME: this is slow and readable, can be accelerated by using
   fancy operations */
int dir_perm(int direction, char d1, char d2)
{
  // vector
  if (d2 == ' ')
  {

    if (direction == 1)
    {

      if (d1 == 'x')
        return 0;
      if (d1 == 'y')
        return 1;
      if (d1 == 'z')
        return 2;
    }
    else if (direction == 2)
    {

      if (d1 == 'x')
        return 1; // y
      if (d1 == 'y')
        return 2; // z
      if (d1 == 'z')
        return 0; // x
    }
    else
    {

      if (d1 == 'x')
        return 2; // z
      if (d1 == 'y')
        return 0; // x
      if (d1 == 'z')
        return 1; // y
    }

    // tensor
  }
  else
  {

    /* define: xx = 0, xy = 1, xz = 2, yy = 3, yz = 4, zz = 5 */
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
        if (d2 == 'y')
          return 3; // yy
        if (d2 == 'z')
          return 4; // yz
      }
      else
      {
        return 5; // zz
      }
    }
    else if (direction == 2)
    {

      // x->y, y->z, z->x
      if (d1 == 'x')
      {
        if (d2 == 'x')
          return 3; // yy
        if (d2 == 'y')
          return 4; // yz
        if (d2 == 'z')
          return 1; // yx
      }
      else if (d1 == 'y')
      {
        if (d2 == 'y')
          return 5; // zz
        if (d2 == 'z')
          return 2; // zx
      }
      else
      {
        return 0; // xx
      }
    }
    else
    {

      // x->z, y->x, z->y
      if (d1 == 'x')
      {
        if (d2 == 'x')
          return 5; // zz
        if (d2 == 'y')
          return 2; // zx
        if (d2 == 'z')
          return 4; // zy
      }
      else if (d1 == 'y')
      {
        if (d2 == 'y')
          return 0; // xx
        if (d2 == 'z')
          return 1; // xy
      }
      else
      {
        return 3; // yy
      }
    }
  }
  errorexit("something is wrong");
}

int dir_perm_b(int direction, char d1)
{

  if (direction == 1)
  {

    if (d1 == 'x')
      return 0;
    if (d1 == 'y')
      return 1;
    if (d1 == 'z')
      return 2;
  }
  else if (direction == 2)
  {

    if (d1 == 'x')
      return 1; // y
    if (d1 == 'y')
      return 2; // z
    if (d1 == 'z')
      return 0; // x
  }
  else
  {

    if (d1 == 'x')
      return 2; // z
    if (d1 == 'y')
      return 0; // x
    if (d1 == 'z')
      return 1; // y
  }
}

/* add fluxes to the RHS at the physical points */
void add_flx1d_rhs(int direction, tB *box, int dir1, int dir2,
                   tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars,
                   double **fluxs1d)
{
  tL *level = ucur->level;
  int addlinear = (c != 0.0l);
  double *mask = level->v[MATTER.INDX_VAR_mask];
  double *mask_A, *mask_B;
  double dtime = level->dt;
  int l, ijk, v;
  double *r[MATTER.NVf], deltaF[MATTER.NVf], *camr[MATTER.NVf];
  double camr_direction, camr_fact;
  int N = (direction == 1) * box->m + (direction == 2) * box->n + (direction == 3) * box->o;
  double h = (direction == 1) * level->dx + (direction == 2) * level->dy + (direction == 3) * level->dz;
  double ooh = 1. / (h);

  int f_ijk_notfinite, sdir;
  tVarList *evolve_vars = get_evolve_vlregister(level);

  if (MATTER.USECAMR)
  {
    mask_A = Ptr(level, "camr_mask_A");
    mask_B = Ptr(level, "camr_mask_B");
    if (addlinear)
      camr_fact = dtime * ooh; /* for euler integration*/
    else
      camr_fact = ooh; /* for all forms of RK, but not rk4g  */
    for (v = 0; v < MATTER.NVf; v++)
    {
      if (VarNComponents(ucur->index[MATTER.INDX_VAR_q + v]) == 3)
      {
        camr[v + 0] = vldataptr(unew, MATTER.INDX_VAR_camr + v + dir_perm_b(direction, 'x'));
        camr[v + 1] = vldataptr(unew, MATTER.INDX_VAR_camr + v + dir_perm_b(direction, 'y'));
        camr[v + 2] = vldataptr(unew, MATTER.INDX_VAR_camr + v + dir_perm_b(direction, 'z'));
        v += 2;
      }
      else
      {
        camr[v] = vldataptr(unew, MATTER.INDX_VAR_camr + v);
      }
    }
  }

  // pointers to rhs vars with permutation, because fluxs1d is computed
  // as it would be in x-dir
  for (v = 0; v < MATTER.NVq; v++)
  {
    if (VarNComponents(ucur->index[MATTER.INDX_VAR_q + v]) == 3)
    {
      r[v + 0] = vldataptr(unew, MATTER.INDX_VAR_q + v + dir_perm_b(direction, 'x'));
      r[v + 1] = vldataptr(unew, MATTER.INDX_VAR_q + v + dir_perm_b(direction, 'y'));
      r[v + 2] = vldataptr(unew, MATTER.INDX_VAR_q + v + dir_perm_b(direction, 'z'));
      v += 2;
    }
    else
    {
      r[v] = vldataptr(unew, MATTER.INDX_VAR_q + v);
    }
  }

  // compute RHS and add (in this case with a minus)
  // only for inner points like for finite differences
  //       dt q + dx f = s

  // Check if fluxs finite at pt l=0
  l = 0;
  f_ijk_notfinite = 0;
  for (v = 0; v < MATTER.NVf; v++)
  {
    if (!finite(fluxs1d[v][l]))
    {
      f_ijk_notfinite++;
#if (0)
      fluxs1d[v][l] = 0.; // set this to 0
#endif
    }
  }
#if (1)
  if (f_ijk_notfinite)
    for (v = 0; v < MATTER.NVf; v++)
      fluxs1d[v][l] = 0.; // set all to 0 (default)
#endif

  // Loop on 1d points
  for (l = 1; l <= N; l++)
  {

    // ijk index
    if (direction == 1)
      ijk = box->noffset + (l - 1) + box->m * dir1 + box->m * box->n * dir2;
    else if (direction == 2)
      ijk = box->noffset + dir1 + box->m * (l - 1) + box->m * box->n * dir2;
    else
      ijk = box->noffset + dir1 + box->m * dir2 + box->m * box->n * (l - 1);

    // Handle ATM pt
    if ((MATTER.USEMASK) && (mask[ijk] > 0.99))
    {
      // add zero to RHS
      for (v = 0; v < MATTER.NVf; v++)
        r[v][ijk] -= 0.;
      // go next point
      continue;
    }

    // Check if flx finite
    f_ijk_notfinite = 0;
    for (v = 0; v < MATTER.NVf; v++)
    {
      if (!finite(fluxs1d[v][l]))
      {
        f_ijk_notfinite++;
#if (0)
        fluxs1d[v][l] = 0.; // set this to 0
#endif
      }
    }
#if (1)
    if (f_ijk_notfinite)
      for (v = 0; v < MATTER.NVf; v++)
        fluxs1d[v][l] = 0.;
        // set all to 0 (default)
#endif

    // Compute div F
    for (v = 0; v < MATTER.NVf; v++)
      deltaF[v] = ooh * (fluxs1d[v][l] - fluxs1d[v][l - 1]);

    // Nans were found above and cured, what we want to do here?
    if (f_ijk_notfinite)
    {

      if (PR)
      {
        printf(" problem with nan's in add_flx1d_rhs\n");
        printf("  %e %e %e    l=%d\n", Ptr(level, "x")[ijk], Ptr(level, "y")[ijk], Ptr(level, "z")[ijk], level->l);
        printf("  var \t F_i F_i-1 \t deltaF_i \n");
        for (v = 0; v < MATTER.NVf; v++)
          printf(" %d \t %e %e %e\n", v, fluxs1d[v][l], fluxs1d[v][l - 1], deltaF[v]);
      }

      if (StopFlxNOTFinite || PR)
      {
        if (ijkinsidefinerlevel(box, ijk))
          printf(" hydro fluxes in rhs not finite: point is inside finer box or in symmetry area\n");
        else
        {
          printf(" hydro fluxes in rhs not finite: point is NOT inside finer box NOR in some symmetry area \n");
#if (StopFlxNOTFinite)
          errorexit(" I stop.");
#endif
        }
      }
    }

    /*              for (v=0; v<MATTER.NVf; v++)
                       camr[v][ijk] += deltaF[v];
  */
    // Add to RHS
    if (addlinear)
      for (v = 0; v < MATTER.NVf; v++)
        r[v][ijk] -= c * deltaF[v];
    else
      for (v = 0; v < MATTER.NVf; v++)
        r[v][ijk] -= deltaF[v];

    // save fluxes for the c_amr_XXX in type B-cells
    if ((MATTER.USECAMR) && (MATTER.CAMRACTIVE))
    {

      if (mask_B[ijk] > 0.75)
      {
        if ((direction == 1 && mask_B[ijk] == 2) ||
            (direction == 2 && mask_B[ijk] == 4) ||
            (direction == 3 && mask_B[ijk] == 6))
        {
          for (v = 0; v < MATTER.NVf; v++)
            camr[v][ijk] += camr_fact * fluxs1d[v][l - 1];
        }

        if ((direction == 1 && mask_B[ijk] == 1) ||
            (direction == 2 && mask_B[ijk] == 3) ||
            (direction == 3 && mask_B[ijk] == 5))
        {
          for (v = 0; v < MATTER.NVf; v++)
            camr[v][ijk] -= camr_fact * fluxs1d[v][l];
        }
      }

      // and in type A-cells
      if (mask_A[ijk] > 0.75)
      {
        if ((direction == 1 && mask_A[ijk] == 2) ||
            (direction == 2 && mask_A[ijk] == 4) ||
            (direction == 3 && mask_A[ijk] == 6))
        {
          for (v = 0; v < MATTER.NVf; v++)
            camr[v][ijk] -= camr_fact * fluxs1d[v][l - 1];
        }

        if ((direction == 1 && mask_A[ijk] == 1) ||
            (direction == 2 && mask_A[ijk] == 3) ||
            (direction == 3 && mask_A[ijk] == 5))
        {
          for (v = 0; v < MATTER.NVf; v++)
            camr[v][ijk] += camr_fact * fluxs1d[v][l];
        }
      }

    } // end loop of c-amr

  } // loop 1d pts
}

/* flat extrapolation for 1d buffer ghosts */
void extrapolate(double *f, int l, int lb, int rb)
{
  // ATTENTION: using linear extrapolation
  // leads to negative rho or a velocity larger than light ...
  // BUT ONLY at MPI or SYMMETRY boundaries
  // => just try 0. order extrapolation
  if (l < lb)
    f[l] = f[lb]; // + (double)(l-lb)*( f[lb+1]-f[lb] );
  if (l > rb)
    f[l] = f[rb]; // + (double)(l-rb)*( f[rb]-f[rb-1] );
}

/* linear/flat extrapolation for 1d buffer ghosts */
void extrapolate2(double *f, int l, int lb, int rb, const int ord)
{
  // ord == 0 : copy
  // ord == 1 : line extrapolation
  if (l < lb)
    f[l] = f[lb] + ord * ((double)(l - lb) * (f[lb + 1] - f[lb]));
  if (l > rb)
    f[l] = f[rb] + ord * ((double)(l - rb) * (f[rb] - f[rb - 1]));
}

/* build 1d arrays for metric, primitives, and conservatives
   extrapolate into ghosts
   call the 1d flux routine pointed by hrsc_flx1d */
void compute_flux1d(int direction, tB *box, int dir1, int dir2,
                    tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars,
                    double **fluxs1d)
{
  tL *level = ucur->level;

  int l, m, n, p, ijk;
  double detginv;

  int Nghost = MATTER.HRSC_NGHOST;
  int Nmax = 0;
  int NVg = 8;
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

    // adm metric, which should be computed by the metric evolution project
    g1d[2][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'x', 'x')][ijk];
    g1d[3][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'x', 'y')][ijk];
    g1d[4][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'x', 'z')][ijk];
    g1d[5][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'y', 'y')][ijk];
    g1d[6][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'y', 'z')][ijk];
    g1d[7][l + 1] = level->v[MATTER.INDX_VAR_gxx + dir_perm(direction, 'z', 'z')][ijk];

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
      else if (VarNComponents(otherVars->index[m]) == 6)
      {
        o1d[m][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'x', 'x'))[ijk];
        o1d[m + 1][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'x', 'y'))[ijk];
        o1d[m + 2][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'x', 'z'))[ijk];
        o1d[m + 3][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'y', 'y'))[ijk];
        o1d[m + 4][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'y', 'z'))[ijk];
        o1d[m + 5][l + 1] = vldataptr(otherVars, m + dir_perm(direction, 'z', 'z'))[ijk];
        m += 5;
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

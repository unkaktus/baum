/* hrsc_numflx1d.c
   generic routine for numerical fluxes 1d
   sbernuz 03/2012 */

#include "bam.h"
#include "matter.h"

#define TINY 1e-32

#define ZERO 0.0
#define HALF 0.5
#define ONE 1.0

// HLL flux
void hrsc_numflx1d_hll(double *uL, double *uR,
                       double *fL, double *fR,
                       double *lamL, double *lamR,
                       double *f, int nf)
{
  int v;

  double aplus = 0.;
  double amins = 0.;

  // max speeds
  for (v = 0; v < nf; v++)
  {
    aplus = DMAX(aplus, DMAX(lamR[v], lamL[v]));
    amins = DMIN(amins, DMIN(lamR[v], lamL[v]));
  }

  aplus = fabs(aplus);
  amins = fabs(amins);

  double oda = ONE / (aplus + amins + TINY);
  double apm = aplus * amins;

  // build fluxes
  for (v = 0; v < nf; v++)
    f[v] = oda * (aplus * fL[v] + amins * fR[v] - apm * (uR[v] - uL[v]));
}

// HLL flux grmhd
void hrsc_numflx1d_hll_grmhd(double *uL, double *uR,
                             double *fL, double *fR,
                             double *lamL, double *lamR,
                             double *f, int nf)
{
  int v;

  double aplus = 0.;
  double amins = 0.;

  // max speeds (excluding the eigenvalues for divergence cleaning)
  for (v = 0; v < nf - 2; v++)
  {
    aplus = DMAX(aplus, DMAX(lamR[v], lamL[v]));
    amins = DMIN(amins, DMIN(lamR[v], lamL[v]));
  }

  aplus = fabs(aplus);
  amins = fabs(amins);

  double oda = ONE / (aplus + amins + TINY);
  double apm = aplus * amins;

  // build fluxes (for divergence cleaning aplus = 1 and amins = -1)
  for (v = 0; v < nf; v++)
  {
    if ((MATTER.q_names_list[v] != "grmhd_Bx") && (MATTER.q_names_list[v] != "grmhd_By") &&
        (MATTER.q_names_list[v] != "grmhd_Bz") && (MATTER.q_names_list[v] != "grmhd_zeta"))
      f[v] = oda * (aplus * fL[v] + amins * fR[v] - apm * (uR[v] - uL[v]));
    else
      f[v] = ONE / (2 + TINY) * (fL[v] + fR[v] - (uR[v] - uL[v]));
  }
}

// LLF flux
void hrsc_numflx1d_llf(double *uL, double *uR,
                       double *fL, double *fR,
                       double *lamL, double *lamR,
                       double *f, int nf)
{
  int v;

  double amax;
  double aplus = 0.;
  double amins = 0.;

  // max speed
  /*
  for(v=0; v<nf; v++)
  amax = DMAX( amax, DMAX(fabs(lamR[v]), fabs(lamL[v])) );
  */

  // max speeds
  for (v = 0; v < nf; v++)
  {
    aplus = DMAX(aplus, DMAX(lamR[v], lamL[v]));
    amins = DMIN(amins, DMIN(lamR[v], lamL[v]));
  }

  aplus = fabs(aplus);
  amins = fabs(amins);
  amax = DMAX(amins, aplus);

  // build fluxes
  for (v = 0; v < nf; v++)
    f[v] = HALF * (fR[v] + fL[v] - amax * (uR[v] - uL[v]));
}

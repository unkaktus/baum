/* maximal.c */
/* Bernd Bruegmann, 7/03 */

#include "bam.h"
#include "Gauge.h"




/* wrapper for maximal_L */
void compute_maximal_L(tL *level, tVarList *vlLu, tVarList *vlu)
{
  tVarList *vlgi = VLPtrEnable1(level, "maximal_gixx");
  tVarList *vlG  = VLPtrEnable1(level, "maximal_Gx");
  tVarList *vlKK = VLPtrEnable1(level, "maximal_KK");

  /* apply boundary conditions */
  set_boundary_elliptic(level, vlu);

  /* interior */
  maximal_L(level, vlLu, vlu, vlgi, vlG, vlKK);

  /* synchronize */
  bampi_vlsynchronize(vlLu);

  /* cleanup */
  vlfree(vlgi);
  vlfree(vlG);
  vlfree(vlKK);
}




/* special operator for spherical symmetry 
   more efficient, in particular because it avoids application of 
   the expensive cartoon boundary condition
*/
void compute_maximal_L1d(tL *level, tVarList *vlLu, tVarList *vlu) 
{
  double *Lu = VLPtr(vlLu, 0);
  double *u  = VLPtr(vlu,  0);
  double *gixx = Ptr(level, "maximal_gixx");
  double *giyy = Ptr(level, "maximal_giyy");
  double *gizz = Ptr(level, "maximal_gizz");
  double *Gz   = Ptr(level, "maximal_Gz");
  double *KK   = Ptr(level, "maximal_KK");
  double *z    = Ptr(level, "z");
  double dz = level->dz;
  double oo2dz = 1/(2*dz);
  double oodz2 = 1/(dz*dz);
  double duz;
  double dduzz;

  set_boundary_elliptic_1d(level, vlu);

  forinner7(level) {
    duz = oo2dz*(-u[ccm] + u[ccp]);
    dduzz = oodz2*(-2.*u[ccc] + u[ccm] + u[ccp]);

    Lu[ccc] = ((gixx[ccc] + giyy[ccc])/z[ccc] - Gz[ccc]) * duz +
              gizz[ccc] * dduzz - KK[ccc]*u[ccc];
  } endfor;

  bampi_vlsynchronize(vlLu);
}




/* compute Gauss-Seidel sweep, calls maximal_L */
void compute_maximal_GS(tL *level, tVarList *vlv, tVarList *vlu)
{
  tVarList *vlgi = VLPtrEnable1(level, "maximal_gixx");
  tVarList *vlG  = VLPtrEnable1(level, "maximal_Gx");
  tVarList *vlKK = VLPtrEnable1(level, "maximal_KK");
  double *u = VLPtr(vlu, 0);
  double *v = VLPtr(vlv, 0);
  double *f    = Ptr(level, "maximal_f");
  double *gixx = Ptr(level, "maximal_gixx");
  double *giyy = Ptr(level, "maximal_giyy");
  double *gizz = Ptr(level, "maximal_gizz");
  double *KK   = Ptr(level, "maximal_KK");
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double lii;
  int i;

  /* interior */
  maximal_L(level, vlv, vlu, vlgi, vlG, vlKK);

  /* GS */
  forinner1(level, i) {
    lii = -KK[i] - 2*(gixx[i]/(dx*dx) + giyy[i]/(dy*dy) + gizz[i]/(dz*dz));
    v[i] = u[i] + (f[i] - v[i])/lii;
  }

  /* synchronize and set boundary */
  bampi_vlsynchronize(vlv);
  set_boundary_elliptic(level, vlv);

  /* cleanup */
  vlfree(vlgi);
  vlfree(vlG);
  vlfree(vlKK);
}




/***************************************************************************/
/* solve maximal slicing equation */
int maximal(tL *level)
{
  tG *g = level->grid;
  tL *lvl, *toplvl = g->level[g->lmin];
  tVarList *vlu = VLPtrEnable1(level, "alpha");
  tVarList *vlf = VLPtrEnable1(level, "maximal_f");
  tVarList *vlr = VLPtrEnable1(level, "maximal_r");
  tVarList *vlc = vlalloc(level);
  tVarList *vlgi = VLPtrEnable1(level, "maximal_gixx");
  tVarList *vlG  = VLPtrEnable1(level, "maximal_Gx");
  tVarList *vlKK = VLPtrEnable1(level, "maximal_KK");
  int pr = Getv("maximal_verbose", "yes");
  int i, l;
  double normres;
  double vgauge;
  double *alpha, *f, *KK;
  int i_alpha = Ind("alpha");
  int i_f     = Ind("maximal_f");
  int i_KK    = Ind("maximal_KK");
  static int firstcall = 1;

  /* make special robin boundary known */
  if (firstcall) {
    firstcall = 0;
    if (Getv("maximal_boundary", "robin"))
      Appends("boundary", "robin");
  }

  /* first experiment: wait for alignement (see main.c for real thing) */
  if (level->l < g->lmax) return 0;
  if (dless(level->time, toplvl->time - level->dt)) return 0;
  if (pr) printf("Computing lapse for maximal slicing\n");

  /* check how often we really want to solve maximal */
  /* fixme: this is not quite what one wants with Berger-Oliger AMR */
  if (toplvl->iteration % Geti("maximal_every")) return 0;

  /* prepare Robin boundary condition, set alpha at infinity to 0 */
  vgauge = VarPropSpeed(Ind("alpha"));
  if (Getv("maximal_boundary", "robin")) {
    if (level->iteration == 0) 
      for (l = g->lmin; l <= g->lmax; l++)
	find_robin_normal(g->level[l]);
    VarNameSetBoundaryInfo("alpha", 0, 1, 0);
  } else
    VarNameSetBoundaryInfo("alpha", 0, 0, 0);  

  /* make variable list for coefficients */
  vlpushvl(vlc, vlgi);
  vlpushvl(vlc, vlG);
  vlpushvl(vlc, vlKK);

  /* precompute coefficients */
  for (l = g->lmin; l <= g->lmax; l++) {
    lvl = g->level[l];
    vlenablelevel(lvl, vlc);
    vlenablelevel(lvl, vlf);
    vlenablelevel(lvl, vlr);
    maximal_init(lvl);

    /* temporary hack because the current Robin formula leads to an
       instability when u_infinity != 0 
    */
    alpha = lvl->v[i_alpha];
    f     = lvl->v[i_f];
    KK    = lvl->v[i_KK];
    forallpoints(lvl, i) {
      alpha[i] -= 1.0;
      f[i] = KK[i];
    }
  }


  /* solve */
  /* bicgstab */
  if (Getv("maximal_solver", "bicgstab")) {
    if (Getv("maximal_use1dstencil", "yes"))
      i = bicgstab(level, vlu, vlf, vlr, vlc, 
		   Geti("maximal_itmax"), Getd("maximal_tolerance"), &normres, 
		   compute_maximal_L1d, DPflatlinear);
    else
      i = bicgstab(level, vlu, vlf, vlr, vlc, 
		   Geti("maximal_itmax"), Getd("maximal_tolerance"), &normres, 
		   compute_maximal_L, DPflatlinear);
    if (i >= 0 && pr) 
      printf("maximal slicing, bicgstab: %5d  %10.3e\n", i, normres);
    if (i < 0)
      printf("maximal slicing, bicgstab failed: %5d  %10.3e\n", i, normres);
  }

  /* multigrid, call with coarsest level as argument */
  else if (Getv("maximal_solver", "multigrid")) {
    multigrid(toplvl, vlu, vlf, vlr, vlc, 
	      Geti("maximal_itmax"), Getd("maximal_tolerance"), &normres, 
	      compute_maximal_L, compute_maximal_GS);
  }
  
  /* default */
  else 
    errorexit("maximal: unknown solver");


  /* finish result */
  for (l = g->lmin; l <= g->lmax; l++) {
    lvl = g->level[l];

    /* set boundary for final result */
    set_boundary_elliptic(lvl, vlu);

    /* back to normal alpha */
    alpha = lvl->v[i_alpha];
    forallpoints(lvl, i) alpha[i] += 1.0;
  }

  /* restore on exit so that radiative works for intermediate 1+log, say */
  VarNameSetBoundaryInfo("alpha", 1, 1, vgauge);

  /* clean up */
  vlfree(vlc);
  vlfree(vlgi);
  vlfree(vlG);
  vlfree(vlKK);
  return 0;
}

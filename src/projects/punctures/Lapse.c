/* Lapse.c */
/* Bernd Bruegmann 2/03 */
/* set lapse appropriate for puncture data 
   in particular, compute maximal slicing
*/

#include "bam.h"
#include "punctures.h"



/* maximal slicing equation specialized to punctures initial data and 
   without singularities:
     Laplace v = (v+1) c,  v = 0 at boundary
   where
     v = psi alpha - 1    or    v = psi alpha - (1 - 1/a)     
     c = 7/8 psi^-8 A_ij A^ij
     cf = c * (1 - 1/a);
*/


/* set coefficient */
void SetPunctureMaxLapseCoefficient(tL *level)
{
  double *a   = Ptr(level, "punctures_a");
  double *b   = Ptr(level, "punctures_b");
  double *u   = Ptr(level, "punctures_u");
  double *psi = Ptr(level, "adm_psi");
  double *c   = PtrEnable(level, "punctures_c");
  double *cf  = PtrEnable(level, "punctures_cf");
  int i;
  double f = -Getd("punctures_lapse_at_puncture");

  /* compute coefficient c from a and b using psi = psi_0 + u */
  forallpoints(level, i) {
    c[i] = 7.0 * pow(psi[i]+u[i], -8.0) * pow(a[i], -7.0) * b[i];
    cf[i] = c[i] * (1 - f/a[i]);
  }
}




/* set lapse */
void SetPunctureMaxLapse(tL *level)
{
  double *alpha = Ptr(level, "alpha");
  double *a = Ptr(level, "punctures_a");
  double *u = Ptr(level, "punctures_u");
  double *v = Ptr(level, "punctures_v");
  double *psi = Ptr(level, "adm_psi");
  int i;
  double f = -Getd("punctures_lapse_at_puncture");

  forallpoints(level, i)
    alpha[i] = (v[i] + 1 - f/a[i]) / (psi[i] + u[i]);
}




/* apply linear elliptic operator:  lu = Laplace v - c v */
void LPunctureMaxLapse(tL *level, tVarList *vllu, tVarList *vlu)
{
  double *c  = Ptr(level, "punctures_c");
  double *u  = VLPtr(vlu, 0);
  double *lu = VLPtr(vllu, 0);
  double cx = 1/(level->dx*level->dx);
  double cy = 1/(level->dy*level->dy);
  double cz = 1/(level->dz*level->dz);
  double cc = -2*(cx + cy + cz);
  int i;
  
  /* apply boundary conditions */
  set_boundary_elliptic(level, vlu);

  /* interior */
  forinner7(level) {
    lu[ccc] = cc * u[ccc] + 
              cx * (u[mcc] + u[pcc]) +
              cy * (u[cmc] + u[cpc]) +
              cz * (u[ccm] + u[ccp]) - c[ccc] * u[ccc];
    /* produces significant round-off differences ?!
       lu[ccc] = (cc - c[ccc]) * u[ccc] +
                 cx * (u[mcc] + u[pcc]) +
                 cy * (u[cmc] + u[cpc]) +
                 cz * (u[ccm] + u[ccp]);
    */
  } endfor;

  /* synchronize */
  bampi_vlsynchronize(vllu);
}




/* linear Gauss-Seidel:  v = u + (f - Lu)/Lii */
/* elliptic operator:  Lu = Laplace v - c v */
void LPunctureMaxLapse_GS(tL *level, tVarList *vlv, tVarList *vlu)
{
  double *c = Ptr(level, "punctures_c");
  double *f = Ptr(level, "punctures_cf");
  double *u = VLPtr(vlu, 0);
  double *v = VLPtr(vlv, 0);
  double cx = 1/(level->dx*level->dx);
  double cy = 1/(level->dy*level->dy);
  double cz = 1/(level->dz*level->dz);
  double cc = -2*(cx + cy + cz);
  double lu, lii;
  
  /* interior */
  forinner7(level) {
    lu = cc * u[ccc] + 
         cx * (u[mcc] + u[pcc]) +
         cy * (u[cmc] + u[cpc]) +
         cz * (u[ccm] + u[ccp]) - c[ccc] * u[ccc];

    lii = cc - c[ccc];

    v[ccc] = u[ccc] + (f[ccc] - lu)/lii;
  } endfor;

  /* synchronize and set boundary */
  bampi_vlsynchronize(vlv);
  set_boundary_elliptic(level, vlv);
}




/* compute lapse for maximal slicing */
void PunctureMaximalSlicing(tL *level)
{
  tVarList *vlv  = VLPtrEnable1(level, "punctures_v");
  tVarList *vlc  = VLPtrEnable1(level, "punctures_c");
  tVarList *vlf  = VLPtrEnable1(level, "punctures_cf");
  tVarList *vlr  = VLPtrEnable1(level, "punctures_r");
  tVarList *vlcoeff = vlalloc(level);
  int itmax = Geti("punctures_itmax");
  double tol = Getd("punctures_tolerance");
  double normres;

  prdivider(0);
  printf("computing maximal slicing lapse for puncture data\n");

  /* set coefficient */
  SetPunctureMaxLapseCoefficient(level);

  /* make variable list for coefficients */
  vlpush(vlcoeff, Ind("punctures_c"));
  vlpush(vlcoeff, Ind("punctures_cf"));

  /* solve */
  /* multigrid */
  if (Getv("punctures_solver", "multigrid"))
    multigrid(level, vlv, vlf, vlr, vlc, itmax, tol, &normres, 
	      LPunctureMaxLapse, LPunctureMaxLapse_GS);

  /* bicgstab */
  else if (Getv("punctures_solver", "bicgstab"))
    bicgstab(level, vlv, vlf, vlr, 0, itmax, tol, &normres, 
	     LPunctureMaxLapse, DPflatlinear);

  /* unknown solver */
  else 
    errorexit("unknown elliptic solver in punctures/Lapse.c");

  /* set lapse */
  set_boundary_elliptic(level, vlv);
  SetPunctureMaxLapse(level);

  /* clean up */
  vlfree(vlcoeff);
  if (Getv("punctures_persist", "no")) {
    VLDisableFree(vlc);
    VLDisableFree(vlf);
    VLDisableFree(vlr);
  }
  printf("computed maximal slicing lapse for puncture data\n");
  prdivider(0);
}

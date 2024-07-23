/* adm_data.c */
/* Bernd Bruegmann 6/05 */
/* provide some simple types of analytic initial data */


#include "bam.h"
#include "adm.h"




/* set Minkowski data */
int adm_Minkowski(tL *level)
{
  tVarList *vl;
  
  /* make list of all variables */
  vl = vlalloc(level);
  vlpush(vl, Ind("adm_gxx"));
  vlpush(vl, Ind("adm_Kxx"));
  vlpush(vl, Ind("alpha"));
  vlpush(vl, Ind("betax"));
  vlpush(vl, Ind("adm_psi"));
  vlpush(vl, Ind("adm_dpsiopsix"));
  vlpush(vl, Ind("adm_ddpsiopsixx"));
  
  /* enable variable list, which initializes all variables to zero */
  vlenable(vl);
  vlfree(vl);

  /* make list of all variables that we will set to one */
  vl = vlalloc(level);
  vlpushone(vl, Ind("adm_gxx"));
  vlpushone(vl, Ind("adm_gyy"));
  vlpushone(vl, Ind("adm_gzz"));
  vlpushone(vl, Ind("alpha"));
  vlpushone(vl, Ind("adm_psi"));
  
  /* set these variables to one */
  vlsetconstant(vl, 1.0);
  vlfree(vl);

  /* set adm matter parts to zero */
  if (Getv("physics","matter")) {
    vl = vlalloc(level);
    vlpush(vl, Ind("adm_rho"));
    vlpush(vl, Ind("adm_Sx"));
    vlpush(vl, Ind("adm_SSxx"));
    vlpush(vl, Ind("adm_ST"));
    vlenable(vl);
    vlfree(vl);
  }
  
  return 0;
}




/* set single black hole puncture data */
int adm_puncture(tL *level)
{
  double *x = level->v[Ind("x")];
  double *y = level->v[Ind("y")];
  double *z = level->v[Ind("z")];
  double *gxx, *gyy, *gzz;
  double *psi, *dp1op, *dp2op, *dp3op;
  double *ddp11op, *ddp12op, *ddp13op;
  double *ddp22op, *ddp23op, *ddp33op;
  double *alpha;
  double m = Getd("mass1");
  int i;
  double p, dp1, dp2, dp3, ddp11, ddp12, ddp13, ddp22, ddp23, ddp33;
  double s1, s3, s5;
  double r, ri, xi, yi, zi;

  /* for simplicity, start by initializing with Minkowski */
  adm_Minkowski(level);

  /* set pointer if we want special initialization of the lapse */
  alpha = (Getv("adm_data", "psiBL^(-2)")) ? Ptr(level, "alpha") : 0;
  alpha = Ptr(level, "alpha");

  /* set non-conformal puncture data */
  if (Getv("adm_data", "nonconformal")) {
    gxx = Ptr(level, "adm_gxx");
    gyy = Ptr(level, "adm_gyy");
    gzz = Ptr(level, "adm_gzz");

    forallpoints(level, i) {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      r  = sqrt(xi*xi+yi*yi+zi*zi);
      if (dequal(r, 0.0)) r = dequaleps/1000.0;
      p = 1.0 + 0.5*m/r;
      gxx[i] = gyy[i] = gzz[i] = p*p*p*p;
      if (alpha) alpha[i] = 1/(p*p);
    }
    return 0;
  }

  /* set conformal factor and its derivatives 
     compare ConformalFactor.c in projects/punctures 
  */
  psi = Ptr(level, "adm_psi");
  dp1op = Ptr(level, "adm_dpsiopsix");
  dp2op = Ptr(level, "adm_dpsiopsiy");
  dp3op = Ptr(level, "adm_dpsiopsiz");
  ddp11op = Ptr(level, "adm_ddpsiopsixx");
  ddp12op = Ptr(level, "adm_ddpsiopsixy");
  ddp13op = Ptr(level, "adm_ddpsiopsixz");
  ddp22op = Ptr(level, "adm_ddpsiopsiyy");
  ddp23op = Ptr(level, "adm_ddpsiopsiyz");
  ddp33op = Ptr(level, "adm_ddpsiopsizz");

  forallpoints(level, i) {
    xi = x[i];
    yi = y[i];
    zi = z[i];

    r  = sqrt(xi*xi+yi*yi+zi*zi);
    if (dequal(r, 0.0)) r = dequaleps/1000.0;
    ri = 1/r;

    s1 = m*ri*(0.5);
    s3 = s1*ri*ri*(-1.0);
    s5 = s3*ri*ri*(-3.0);

    p = s1;

    dp1 = xi*s3; 
    dp2 = yi*s3; 
    dp3 = zi*s3; 

    ddp11 = xi*xi*s5 + s3;
    ddp12 = xi*yi*s5;
    ddp13 = xi*zi*s5;
    ddp22 = yi*yi*s5 + s3;
    ddp23 = yi*zi*s5;
    ddp33 = zi*zi*s5 + s3;
  
    p += 1.0;

    psi[i] = p;
    dp1op[i] = dp1/p;
    dp2op[i] = dp2/p;
    dp3op[i] = dp3/p;
    ddp11op[i] = ddp11/p;
    ddp12op[i] = ddp12/p;
    ddp13op[i] = ddp13/p;
    ddp22op[i] = ddp22/p;
    ddp23op[i] = ddp23/p;
    ddp33op[i] = ddp33/p;

    if (alpha) alpha[i] = 1/(p*p);
  }

  return 0;
}




/* gauge wave */
int adm_gauge_wave(tL *level)
{
    adm_Minkowski(level);
    
    double *alpha = level->v[Ind("alpha")];
    double *x = level->v[Ind("x")];
    double *y = level->v[Ind("y")];
    double *z = level->v[Ind("z")];

    double v1 = Getd("gauge_wave_amplitude");
    double v2 = Getd("gauge_wave_diameter");
    
    forallpoints_ijk(level) {
        alpha[ijk] = 1. + v1*exp(-v2*(x[ijk]*x[ijk]+y[ijk]*y[ijk]+z[ijk]*z[ijk]));
    } endfor_ijk;
    printf("  set lapse to one + gauge wave\n");
    
    return 0;
}














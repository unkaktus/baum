/* forMovingPunctures.c   */
/* Wolfgang Tichy 11/2005 */

#include "bam.h"
#include "punctures.h"


/* Absorb psi into 3-metric */
int Absorb_psi(tL *level)
{
  printf("Moving conformal factor psi into 3-metric...  ");
  Move_psi_intoMetric(level, Ind("adm_gxx"), 
                      Ind("adm_psi"), Ind("adm_dpsiopsix"), Ind("adm_ddpsiopsixx"));
  printf("Done.\n");

  return 0;
}


/* move psi into 3-metric */
void Move_psi_intoMetric(tL* level, int i_g, 
                         int i_psi, int i_dpsiopsi, int i_ddpsiopsi)
{
  double *g11 = level->v[i_g+0];
  double *g12 = level->v[i_g+1];
  double *g13 = level->v[i_g+2];
  double *g22 = level->v[i_g+3];
  double *g23 = level->v[i_g+4];
  double *g33 = level->v[i_g+5];
  double *psi = level->v[i_psi];
  double *dpsiopsi1 = level->v[i_dpsiopsi+0];
  double *dpsiopsi2 = level->v[i_dpsiopsi+1];
  double *dpsiopsi3 = level->v[i_dpsiopsi+2];
  double *ddpsiopsi11 = level->v[i_ddpsiopsi+0];
  double *ddpsiopsi12 = level->v[i_ddpsiopsi+1];
  double *ddpsiopsi13 = level->v[i_ddpsiopsi+2];
  double *ddpsiopsi22 = level->v[i_ddpsiopsi+3];
  double *ddpsiopsi23 = level->v[i_ddpsiopsi+4];
  double *ddpsiopsi33 = level->v[i_ddpsiopsi+5];
  double psi4;
  long ccc;

  forallpoints(level,ccc)
  {
    /* g = psi^4 gb */
    psi4 = pow(psi[ccc],4);
    g11[ccc] = psi4 * g11[ccc];
    g12[ccc] = psi4 * g12[ccc];
    g13[ccc] = psi4 * g13[ccc];
    g22[ccc] = psi4 * g22[ccc];
    g23[ccc] = psi4 * g23[ccc];
    g33[ccc] = psi4 * g33[ccc];

    psi[ccc] = 1.0;
    dpsiopsi1[ccc] = dpsiopsi2[ccc] = dpsiopsi3[ccc] = 0.0;
    ddpsiopsi11[ccc] = ddpsiopsi12[ccc] = ddpsiopsi13[ccc] = 0.0;
    ddpsiopsi22[ccc] = ddpsiopsi23[ccc] = ddpsiopsi33[ccc] = 0.0;
  }
}

/* set alpha=psi^(-2) */
void Set_alpha_psip(tL *level, int psipower)
{
  double *psi   = level->v[Ind("adm_psi")];
  double *alpha = level->v[Ind("alpha")];
  long ccc;

  printf("Setting alpha=psi^(psipower).\n");
  forallpoints(level,ccc)
    alpha[ccc] = pow(psi[ccc],psipower);
}

/* set alpha= 1/[ 1 + m1/(N r1) + m2/(N r2) ]^N */
void Set_alpha_rtoN_atPunc(tL *level)
{
  double *alpha = level->v[Ind("alpha")];
  double *x = level->v[Ind("x")];
  double *y = level->v[Ind("y")];
  double *z = level->v[Ind("z")];
  double N  = Getd("punctures_lapse_rPower_atPunc");
  double m1 = Getd("bhmass1");
  double m2 = Getd("bhmass2");
  double x1 = Getd("bhx1");
  double y1 = Getd("bhy1");
  double z1 = Getd("bhz1");
  double x2 = Getd("bhx2");
  double y2 = Getd("bhy2");
  double z2 = Getd("bhz2");
  double r1, r2;
  long ccc;

  printf("Setting alpha=1/[ 1 + m1/(%f r1) + m2/(%f r2) ]^%f.\n", N, N, N);
  forallpoints(level,ccc)
  {
    r1 = sqrt(  (x[ccc]-x1)*(x[ccc]-x1) + (y[ccc]-y1)*(y[ccc]-y1)
               +(z[ccc]-z1)*(z[ccc]-z1)  );
    r2 = sqrt(  (x[ccc]-x2)*(x[ccc]-x2) + (y[ccc]-y2)*(y[ccc]-y2)
               +(z[ccc]-z2)*(z[ccc]-z2)  );
    alpha[ccc] = pow(1.0 + m1/(N*r1) + m2/(N*r2), -N);
  }
}

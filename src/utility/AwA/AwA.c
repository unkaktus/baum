/* AwA.c */
/* mth 06/2011 */
/* dmh 06/2011 */
/* sb  05/2019 */
/* no  07/2019 */

#include "bam.h"
#include "AwA.h"

#include <time.h>

#define PR 0

/* Macros for wave profile */
#define  SINWAVE(A,d,x) ( (A)*      sin(2.*PI*(x)/(d) ) )
#define DSINWAVE(A,d,x) (-(A)*2.*PI*cos(2.*PI*(x)/(d) ) )

//=============== 3D Robust Stability Test ================
//=========================================================
int AwA_robust_stability(tL *level)
{
  /* Set ADM vars to flat spacetime */
  adm_Minkowski(level);

  /* Ptrs to ADM variables and gauge */
  double *alpha = Ptr(level, "alpha");
  double *betax = Ptr(level, "betax");
  double *betay = Ptr(level, "betay");
  double *betaz = Ptr(level, "betaz");
  double *gxx   = Ptr(level, "adm_gxx");
  double *gxy   = Ptr(level, "adm_gxy");
  double *gxz   = Ptr(level, "adm_gxz");
  double *gyy   = Ptr(level, "adm_gyy");
  double *gyz   = Ptr(level, "adm_gyz");
  double *gzz   = Ptr(level, "adm_gzz");
  double *Kxx   = Ptr(level, "adm_Kxx");
  double *Kxy   = Ptr(level, "adm_Kxy");
  double *Kxz   = Ptr(level, "adm_Kxz");
  double *Kyy   = Ptr(level, "adm_Kyy");
  double *Kyz   = Ptr(level, "adm_Kyz");
  double *Kzz   = Ptr(level, "adm_Kzz");

  /* Adding noise at inner points */
  const double eps = Getd("AwA_amplitude");
  double *var[16] = { alpha, betax, betay, betaz,
                      gxx, gxy, gxz, gyy, gyz, gzz,
                      Kxx, Kxy, Kxz, Kyy, Kyz, Kzz };
  int v;

  forinnerpoints_ijk(level) {
    for (v = 0; v < 16; v++) {
      double r = (double) ((rand() / (RAND_MAX + 1.0)));
      var[v][ijk] += eps * (2.*r - 1.);
    }
  } endfor_ijk;

  return 1;
}
//============ End of 3D Robust Stability Test ============

//================ Simple 3D Gauge Wave Test =================
//=========================================================
int AwA_simple_gauge_wave(tL *level)
{
    adm_Minkowski(level);

    double *alpha = level->v[Ind("alpha")];
    double *x = level->v[Ind("x")];
    double *y = level->v[Ind("y")];
    double *z = level->v[Ind("z")];

    double A = Getd("AwA_amplitude");
    double s = Getd("AwA_sigma");

    forallpoints_ijk(level) {
        //alpha[ijk] = 1. + A*exp(-(x[ijk]*x[ijk]+y[ijk]*y[ijk]+z[ijk]*z[ijk])/SQR(s));
        alpha[ijk] = 1. + A * pow(sin(2.*PI*x[ijk]),6) * pow(sin(2.*PI*y[ijk]),6) * pow(sin(2.*PI*z[ijk]),6);
    } endfor_ijk;
    printf("  set lapse to one + gauge wave\n");

    return 0;
}
//============= End of simple 3D Gauge Wave Test =============

//================= 1D Gauge Wave Test ====================
//=========================================================
int AwA_gauge_wave1(tL *level)
{
  /* Set ADM vars to flat spacetime */
  adm_Minkowski(level);
  
  /* Ptrs to AMD variables and gauge */
  double *alpha = Ptr(level, "alpha");
  double *gxx   = Ptr(level, "adm_gxx");
  double *Kxx   = Ptr(level, "adm_Kxx");

  const double *xp = Ptr(level, "x");
  double x, b, dot_b;

  const double A = Getd("AwA_amplitude");
  const double d = Getd("AwA_sigma");

  /* Set gauge wave at inner points along x*/
    
//    forinnerpoints_ijk(level) {
  forallpoints_ijk(level) {

    x     =  xp[i];
    b     =  SINWAVE(A,d,x);
    dot_b = DSINWAVE(A,d,x);

    gxx[ijk] -= b;
    Kxx[ijk] += 0.5*dot_b/sqrt(1. - b);
    alpha[ijk] = sqrt(1. - b);
      
  } endfor_ijk;
  
  return 1;
}
//=============== End of 1D Gauge Wave Test ===============

//================ 1D Linear Wave Test ====================
//=========================================================
int AwA_linear_wave(tL *level)
{
  /* Set ADM vars to flat spacetime */
  adm_Minkowski(level);

  /* Ptrs to AMD variables and gauge */
  double *alpha = Ptr(level, "alpha");
  double *gyy   = Ptr(level, "adm_gyy");
  double *gzz   = Ptr(level, "adm_gzz");
  double *Kyy   = Ptr(level, "adm_Kyy");
  double *Kzz   = Ptr(level, "adm_Kzz");

  const double *xp = Ptr(level, "x");
  double x, b, dot_b;

  const double A = Getd("AwA_amplitude");
  const double d = Getd("AwA_sigma");

  if (Getv("AwA", "linear_wave1")) {

  /* Set linear wave at inner points along x*/

//    forinnerpoints_ijk(level) {
    forallpoints_ijk(level) {

      x     =  xp[i];
      b     =  SINWAVE(A,d,x);
      dot_b = DSINWAVE(A,d,x);

      gyy[ijk] += b;
      gzz[ijk] -= b;
      Kyy[ijk] += 0.5*dot_b;
      Kzz[ijk] -= 0.5*dot_b;
      alpha[ijk] = 1.;

    } endfor_ijk;
  }

  return 1;
}
//============== End of 1D Linear Wave Test ===============

/* Old implementation of robust stability test, 
   Is it still working? Please test. */
int AwA_add_highfrequ(tL *level)
{
  printf("adding high frequency noise\n");
  if (IndLax("alpha")==-1)
    errorexit("no alpha defined, AwA does not work");
  
  srand(time(0));

  int nv;
  double *var,r;
  double eps;
  //double epsl_met  = Getd("AwA_amplitude_met");
  //double epsl_extr = Getd("AwA_amplitude_extr");
  double epsl_met  = Getd("AwA_amplitude");
  double epsl_extr = Getd("AwA_amplitude");

  tVarList *u = get_evolve_vlregister(level); 

  for (nv=0; nv<u->n; nv++) {
    var = level->v[ u->index[nv] ];

    eps = epsl_extr;

    /* determine metric components */
    if (strcmp( VarName(IndComponent0(u->index[nv])), "bssn_gxx" )==0) eps = epsl_met;
    if (strcmp( VarName(IndComponent0(u->index[nv])), "alpha" )==0) eps = epsl_met;
    if (strcmp( VarName(IndComponent0(u->index[nv])), "betax" )==0) eps = epsl_met;

    if (PR) printf("  %d %s\n", nv,VarName( u->index[nv] ));
    
    /* add high frequ */
    forallpoints_ijk(level) {
      r = (double) ((rand() / (RAND_MAX + 1.0)));
      var[ijk] += eps * (2.*r-1.);
    } endfor_ijk;
    
  }
  
  if (PR) {
    printf("  set metric   components to var = var + %e * [-1,+1]\n", epsl_met);
    printf("  set remaining variables to var = var + %e * [-1,+1]\n", epsl_extr);
  }
  
  /* set periodic boundary condition */
  tVarList *vl = vlalloc(level);
  vlpush(vl,Ind("bssn_gxx"));
  vlpush(vl,Ind("alpha"));
  vlpush(vl,Ind("betax"));
  set_boundary(vl,vl,1,vl);
  
  return 1;
}

/* Old implementation of Dplus norm
   Is it still working? Please test & generalize */
int AwA_compute_Dplus(tL *level)
{
  int l,m,n;
  double tmp;
  double delta[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  
  enablevar(level, Ind("Dplus"));
  enablevar(level, Ind("Dplusx"));
  double* Dplus  = level->v[Ind("Dplus")];
  double* Dplusx = level->v[Ind("Dplusx")];
  
  double* b[3];
  double* G[3];
  double* g[3][3];
  double* A[3][3];
  
  double* a = level->v[Ind("alpha")];
  double* c = level->v[Ind("bssn_chi")];
  double* K = level->v[Ind("bssn_K")];
  double* T = NULL;
  
  b[0]    = level->v[Ind("betax")];
  b[1]    = level->v[Ind("betay")];
  b[2]    = level->v[Ind("betaz")];
  
  G[0]    = level->v[Ind("bssn_Gx")];
  G[1]    = level->v[Ind("bssn_Gy")];
  G[2]    = level->v[Ind("bssn_Gz")];
  
  g[0][0] = level->v[Ind("bssn_gxx")];
  g[0][1] = g[1][0] = level->v[Ind("bssn_gxy")];
  g[0][2] = g[2][0] = level->v[Ind("bssn_gxz")];
  g[1][1] = level->v[Ind("bssn_gyy")];
  g[1][2] = g[2][1] = level->v[Ind("bssn_gyz")];
  g[2][2] = level->v[Ind("bssn_gzz")];
  
  A[0][0] = level->v[Ind("bssn_Axx")];
  A[0][1] = A[1][0] = level->v[Ind("bssn_Axy")];
  A[0][2] = A[2][0] = level->v[Ind("bssn_Axz")];
  A[1][1] = level->v[Ind("bssn_Ayy")];
  A[1][2] = A[2][1] = level->v[Ind("bssn_Ayz")];
  A[2][2] = level->v[Ind("bssn_Azz")];
  
  if (Getv("physics","Z4d") || Getv("physics","z4"))
    T = level->v[Ind("Z4d_Theta")];
  
  double ooh = 1./level->dx;
  int dijk[3];
  
  if (PR) printf("compute D+\n");    
  
  forinnerpoints_ijk(level) {
   
    dijk[0] = di;
    dijk[1] = dj;
    dijk[2] = dk;
    
    Dplus[ijk] = 0.;        
    
    // alpha
    tmp = a[ijk]-1.;
    Dplus[ijk] += tmp*tmp;
    
    // beta
    tmp = 0.;
    for (l=0; l<3; l++)
      tmp += b[l][ijk];
    Dplus[ijk] += tmp*tmp;
    
    // chi
    tmp = c[ijk]-1.;
    Dplus[ijk] += tmp*tmp;
    
    // metric
    tmp = 0.;
    for (l=0; l<3; l++)
      for (m=0; m<3; m++)
        tmp += (g[l][m][ijk]-delta[l][m]);
    Dplus[ijk] += tmp*tmp;
    
    // D+ alpha
    tmp = 0.;
    for (l=0; l<3; l++)
      tmp += (a[ijk+dijk[l]]-a[ijk])*ooh;
    Dplus[ijk] += tmp*tmp;
    
    // D+ beta
    tmp = 0.;
    for (l=0; l<3; l++)
      for (m=0; m<3; m++)
        tmp += (b[l][ijk+dijk[l]]-b[l][ijk])*ooh;
    Dplus[ijk] += tmp*tmp;
    
    // D+ chi
    tmp = 0.;
    for (l=0; l<3; l++)
      tmp += (c[ijk+dijk[l]]-c[ijk])*ooh;
    Dplus[ijk] += tmp*tmp;
    
    // D+ metric
    tmp = 0.;
    for (l=0; l<3; l++)
      for (m=0; m<3; m++)
        for (n=0; n<3; n++)
          tmp += (g[m][n][ijk+dijk[l]]-g[m][n][ijk]);
    Dplus[ijk] += tmp*tmp;
        
    // Khat
    tmp = K[ijk];
    Dplus[ijk] += tmp*tmp;
    
    // A
    tmp = 0.;
    for (l=0; l<3; l++)
      for (m=0; m<3; m++)
        tmp += A[l][m][ijk];
    Dplus[ijk] += tmp*tmp;
    
    // G
    tmp = 0.;
    for (l=0; l<3; l++)
      tmp += G[l][ijk];
    Dplus[ijk] += tmp*tmp;
    
    // Theta
    if (T!=NULL) {
      tmp = T[ijk];
      Dplus[ijk] += tmp*tmp;
    }    
    
  } endfor_ijk;
      
  // FIXME: this is not generic!!!!
  forallpoints_ijk(level) {
    
    // assume periodic has 2 ghosts
    if ((i<2) || (i>imax-2) || (j<2) || (j>jmax-2) || (k<2) || (k>imax-2))
      Dplus[ijk] = 0.;
    
    // now only for the y=z=0 line
    if ((j==(jmax+1)/2) && (k==(kmax+1)/2) )
      Dplusx[ijk] = Dplus[ijk];
    else 
      Dplusx[ijk] = 0.;
  
  } endfor_ijk;    
  
  if (PR) printf("finished compute D+\n");   
  return 1;
}








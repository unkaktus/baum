/* moving_punctures_spin.c */
/* Bernd Bruegmann, 2/2008 */

/* compute spin based on local puncture information

   UNFINISHED, work in progress, do not use
*/ 

#include "bam.h"
#include "Gauge.h"


/* macros for some common operations
   use with care, check with gcc -E 
*/
#define fori for (i = 0; i < 3; i++)
#define forj for (j = 0; j < 3; j++)
#define fork for (k = 0; k < 3; k++)
#define forl for (l = 0; l < 3; l++)
#define sumi(a,b)  {a = 0; fori a += b;}
#define sumj(a,b)  {a = 0; forj a += b;}




/* spin */
int moving_puncture_spin(tL *level)
{
  double *px = Ptr(level, "x");
  double *py = Ptr(level, "y");
  double *pz = Ptr(level, "z");
  double x[3];
  
  double *chi = Ptr(level, "bssn_chi");
  double psim2;
  double psipower = Getd("bssn_chi_psipower");

  double *gxx = Ptr(level, "bssn_gxx");
  double *gxy = Ptr(level, "bssn_gxy");
  double *gxz = Ptr(level, "bssn_gxz");
  double *gyy = Ptr(level, "bssn_gyy");
  double *gyz = Ptr(level, "bssn_gyz");
  double *gzz = Ptr(level, "bssn_gzz");
  double g[3][3];
  double gi[3][3];

  double *Axx = Ptr(level, "bssn_Axx");
  double *Axy = Ptr(level, "bssn_Axy");
  double *Axz = Ptr(level, "bssn_Axz");
  double *Ayy = Ptr(level, "bssn_Ayy");
  double *Ayz = Ptr(level, "bssn_Ayz");
  double *Azz = Ptr(level, "bssn_Azz");
  double A[3][3];

  double *rpsi2 = PtrEnable(level, "mprpsi2");
  double *sdh   = PtrEnable(level, "mpsdh");
  double *sx = PtrEnable(level, "mpsx");
  double *sy = PtrEnable(level, "mpsy");
  double *sz = PtrEnable(level, "mpsz");  
  double *tx = PtrEnable(level, "mptx");
  double *ty = PtrEnable(level, "mpty");
  double *tz = PtrEnable(level, "mptz");
  double *ux = PtrEnable(level, "mpux");
  double *uy = PtrEnable(level, "mpuy");
  double *uz = PtrEnable(level, "mpuz");
  double *vx = PtrEnable(level, "mpvx");
  double *vy = PtrEnable(level, "mpvy");
  double *vz = PtrEnable(level, "mpvz");

  double r, r2;
  double nx, ny, nz;
  double Ax, Ay, Az;

  double nn, rho, sqrtdeth;
  double As[3], f[3], n[3], s[3], u[3];
  double dxdX[3][3], h[3][3];

  int i, j, k, l, ijk;


  /* sanity check, could be relaxed */
  if (!Getv("bssn_moving_punctures", "yes")) return 0;
  if (!Getv("bssn_moving_punctures_type", "chi")) {
    printf("moving_puncture_spin: for now require punctures type chi\n");
    return 0;
  }

  /* for all points */
  forallpoints(level, ijk) {

    /* conformal factor:  chi = psi^N, psi^(-2) = chi^(-2/N) */
    /* consider psi^(-2) which is regular at r = 0 */
    psim2 = pow(chi[ijk], -2.0/psipower);

    /* flat approximation */
    if (1) {
      r = sqrt(px[ijk]*px[ijk] + py[ijk]*py[ijk] + pz[ijk]*pz[ijk]);
      if (dequal(r,0)) r = dequaleps;
      
      nx = px[ijk]/r;
      ny = py[ijk]/r;
      nz = pz[ijk]/r;
      
      Ax = Axx[ijk] * nx + Axy[ijk] * ny + Axz[ijk] * nz;
      Ay = Axy[ijk] * nx + Ayy[ijk] * ny + Ayz[ijk] * nz;
      Az = Axz[ijk] * nx + Ayz[ijk] * ny + Azz[ijk] * nz;
      
      /* ti = eijk nj Akl nl */
      tx[ijk] = ny * Az - nz * Ay;
      ty[ijk] = nz * Ax - nx * Az;
      tz[ijk] = nx * Ay - ny * Ax;

      /* ui = rescale ti by radius and conformal factor */
      nn = pow(r/psim2, 3.0) / 3.0;
      ux[ijk] = nn * tx[ijk];
      uy[ijk] = nn * ty[ijk];
      uz[ijk] = nn * tz[ijk];
    }
    
    /* read data at one grid point into local variables */
    x[0] = px[ijk];
    x[1] = py[ijk];
    x[2] = pz[ijk];
    g[0][0] = gxx[ijk];
    g[1][1] = gyy[ijk];
    g[2][2] = gzz[ijk];
    g[0][1] = g[1][0] = gxy[ijk];
    g[0][2] = g[2][0] = gxz[ijk];
    g[1][2] = g[2][1] = gyz[ijk];
    A[0][0] = Axx[ijk];
    A[1][1] = Ayy[ijk];
    A[2][2] = Azz[ijk];
    A[0][1] = A[1][0] = Axy[ijk];
    A[0][2] = A[2][0] = Axz[ijk];
    A[1][2] = A[2][1] = Ayz[ijk];

    /* gi = inverse metric */
    invg(g[0][0], g[0][1], g[0][2], g[1][1], g[1][2], g[2][2],
	 &gi[0][0], &gi[0][1], &gi[0][2], &gi[1][1], &gi[1][2], &gi[2][2]);
    gi[1][0] = gi[0][1];
    gi[2][0] = gi[0][2];
    gi[2][1] = gi[1][2];

    /* fixme: shift center depending on which puncture we are working on */

    /* r = radius, guard against r = 0 */
    sumi(r2, x[i]*x[i]);
    if (dequal(r2,0)) r = dequaleps;
    else              r = sqrt(r2);
   
    /* n[i] = cartesian normal covector */
    fori n[i] = x[i]/r;

    /* s[i] = normal vector to sphere normalized by conformal metric */
    fori sumj(u[i], gi[i][j] * n[j]);
    sumi(nn, n[i] * u[i]);
    nn = sqrt(fabs(nn));
    fori s[i] = u[i]/nn;

    /* As[i] = projection of Aij into one normal direction */
    fori sumj(As[i], A[i][j] * s[j]);

    /* f[i] = integrand without area factor = eijk nj Akl sl */
    f[0] = n[1] * As[2] - n[2] * As[1];
    f[1] = n[2] * As[0] - n[0] * As[2];
    f[2] = n[0] * As[1] - n[1] * As[0];

    /* conformal factor
       normalization of s[i] still needs r psi^2
       area element contains r^2 psi^4
       => overall (r psi^2)^3 = R0^3 assuming spherical moving puncture
    */
    fori f[i] *= pow(r/psim2, 3.0);

    /* fixme: add factor sqrt(det L L gtilde) */

    /* transformation to spherical coordinates
       see e.g. utility/coordtrans/CoordTrans.nb
    */
    rho = x[0]*x[0] + x[1]*x[1];
    if (dequal(rho,0)) rho = dequaleps; else rho = sqrt(rho);
    dxdX[0][0] = x[0]/r;          // dxdr
    dxdX[0][1] = x[1]/r;          // dydr
    dxdX[0][2] = x[2]/r;          // dzdr
    dxdX[1][0] = x[0]*x[2]/rho;   // dxdt
    dxdX[1][1] = x[1]*x[2]/rho;   // dydt
    dxdX[1][2] = -rho;            // dzdt
    dxdX[2][0] = -x[1];           // dxdp
    dxdX[2][1] = x[0];            // dydp
    dxdX[2][2] = 0;               // dzdp

    /* sqrtdeth = sqrt(det L L h) 
       3x3 for simplicity, we really just need the 2x2 metric
    */
    fori forj h[i][j] = 0;
    fori forj fork forl h[i][j] += dxdX[i][k] * dxdX[j][l] * g[k][l];
    sqrtdeth = sqrt(fabs(h[1][1]*h[2][2] - h[1][2]*h[2][1]));

    /* normalize for directly reading off spin on coordinate axis */
    fori f[i] /= 3;

    /* result with conformal factor and r, but without area element */
    sx[ijk] = f[0];
    sy[ijk] = f[1];
    sz[ijk] = f[2];
    rpsi2[ijk] = r/psim2;

    /* result with proper area element */
    /* remove sintheta = rho/r for reading off on coordinate axis */
    fori f[i] *= sqrtdeth/r2 * r/rho;
    vx[ijk] = f[0];
    vy[ijk] = f[1];
    vz[ijk] = f[2];
    sdh[ijk] = sqrtdeth/r2 * r/rho;
  }

  return 0;
}





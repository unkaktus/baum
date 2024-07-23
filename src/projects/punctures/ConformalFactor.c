/* ConformalFactor.c */
/* Bernd Bruegmann 6/02 */
/* based on brbr.c, Bernd Bruegmann, 11/97 */


#include "bam.h"
#include "punctures.h"



/* Misner-Wheeler-Brill-Lindquist conformal factor psi and derivatives over psi

   psi(x,y,z) = 1 + 1/2 Sum_{i} m_i / |r - c_i|

*/
double mwblpsi(double x, double y, double z)
{
  int i;
  double r, xi, yi, zi;
  double p = 1;

  if (fisheye) {coordtrans(&x, &y, &z);}
  
  for (i = 0; MBL[i] != 0.0 && i < N_MWBL; i++) {
    xi = x-CBL[i][0];
    yi = y-CBL[i][1];
    zi = z-CBL[i][2];

    r  = sqrt(xi*xi+yi*yi+zi*zi);
    if (dequal(r, 0.0)) r = EPSZERO;

    p += 0.5 * MBL[i] / r;
  }
  return p;
}




 /* psi and its derivatives, finite difference calculation */
void fdmwblpsiofxyz(double x, double y, double z, 
		    double *psi,
		    double *dp1op, double *dp2op, double *dp3op,
		    double *ddp11op, double *ddp12op, double *ddp13op,
		    double *ddp22op, double *ddp23op, double *ddp33op)
{
  double dx = fisheyedx, dy = fisheyedy, dz = fisheyedz;
  double ccc = mwblpsi(x   , y   , z   );
  double mcc = mwblpsi(x-dx, y   , z   );
  double pcc = mwblpsi(x+dx, y   , z   );
  double cmc = mwblpsi(x   , y-dy, z   );
  double cpc = mwblpsi(x   , y+dy, z   );
  double ccm = mwblpsi(x   , y   , z-dz);
  double ccp = mwblpsi(x   , y   , z+dz);
  double mmc = mwblpsi(x-dx, y-dy, z   );
  double mpc = mwblpsi(x-dx, y+dy, z   );
  double mcm = mwblpsi(x-dx, y   , z-dz);
  double mcp = mwblpsi(x-dx, y   , z+dz);
  double pmc = mwblpsi(x+dx, y-dy, z   );
  double ppc = mwblpsi(x+dx, y+dy, z   );
  double pcm = mwblpsi(x+dx, y   , z-dz);
  double pcp = mwblpsi(x+dx, y   , z+dz);
  double cmm = mwblpsi(x   , y-dy, z-dz);
  double cmp = mwblpsi(x   , y-dy, z+dz);
  double cpm = mwblpsi(x   , y+dy, z-dz);
  double cpp = mwblpsi(x   , y+dy, z+dz);
  
  *psi = ccc;
  *dp1op = (pcc - mcc)/(2*dx*ccc);
  *dp2op = (cpc - cmc)/(2*dy*ccc);
  *dp3op = (ccp - ccm)/(2*dz*ccc);
  *ddp11op = (pcc - 2*ccc + mcc)/(dx*dx*ccc);
  *ddp22op = (cpc - 2*ccc + cmc)/(dy*dy*ccc);
  *ddp33op = (ccp - 2*ccc + ccm)/(dz*dz*ccc);
  *ddp12op = (ppc - pmc - mpc + mmc)/(4*dx*dy*ccc);
  *ddp13op = (pcp - pcm - mcp + mcm)/(4*dx*dz*ccc);
  *ddp23op = (cpp - cpm - cmp + cmm)/(4*dy*dz*ccc);
}




/* psi and its derivatives, analytic calculation */
void mwblpsiofxyz(double x, double y, double z, double *psi,
		  double *dp1op, double *dp2op, double *dp3op,
		  double *ddp11op, double *ddp12op, double *ddp13op,
		  double *ddp22op, double *ddp23op, double *ddp33op)
{
  int i;
  double p, dp1, dp2, dp3, ddp11, ddp12, ddp13, ddp22, ddp23, ddp33;
  double s1, s3, s5;
  double r, ri, xi, yi, zi;

  if (fisheye) {
    fdmwblpsiofxyz(x, y, z, psi, dp1op, dp2op, dp3op,
		   ddp11op, ddp12op, ddp13op, ddp22op, ddp23op, ddp33op);
    return;
  }

  p = dp1 = dp2 = dp3 = ddp11 = ddp12 = ddp13 = ddp22 = ddp23 = ddp33 = 0;

  for (i = 0; MBL[i] != 0.0 && i < N_MWBL; i++) {
    xi = x-CBL[i][0];
    yi = y-CBL[i][1];
    zi = z-CBL[i][2];

    r  = sqrt(xi*xi+yi*yi+zi*zi);
    if (dequal(r, 0.0)) r = EPSZERO;
    ri = 1/r;

    s1 = MBL[i]*ri*(0.5);
    s3 = s1*ri*ri*(-1.0);
    s5 = s3*ri*ri*(-3.0);

    p += s1;

    dp1 += xi*s3; 
    dp2 += yi*s3; 
    dp3 += zi*s3; 

    ddp11 += xi*xi*s5 + s3;
    ddp12 += xi*yi*s5;
    ddp13 += xi*zi*s5;
    ddp22 += yi*yi*s5 + s3;
    ddp23 += yi*zi*s5;
    ddp33 += zi*zi*s5 + s3;
  }
  p += 1.0;

  *psi = p;
  *dp1op = dp1/p;
  *dp2op = dp2/p;
  *dp3op = dp3/p;
  *ddp11op = ddp11/p;
  *ddp12op = ddp12/p;
  *ddp13op = ddp13/p;
  *ddp22op = ddp22/p;
  *ddp23op = ddp23/p;
  *ddp33op = ddp33/p;
}





/* set MWBL, conformal metric: gab = dab, Kab = 0, psi and derivs */
void SetMWBL(tL *level)
{
  int i, j;
  double *x = level->v[Ind("x")];
  double *y = level->v[Ind("y")];
  double *z = level->v[Ind("z")];
  double *p[10];
  double *g11, *g12, *g13, *g22, *g23, *g33;
  double *K11, *K12, *K13, *K22, *K23, *K33;

  i = Ind("adm_gxx");
  enablevar(level, i);
  g11 = level->v[i++];
  g12 = level->v[i++];
  g13 = level->v[i++];
  g22 = level->v[i++];
  g23 = level->v[i++];
  g33 = level->v[i++];

  i = Ind("adm_Kxx");
  enablevar(level, i);
  K11 = level->v[i++];
  K12 = level->v[i++];
  K13 = level->v[i++];
  K22 = level->v[i++];
  K23 = level->v[i++];
  K33 = level->v[i++];

  
  i = Ind("adm_psi");
  for (j = 0; j < 10; j++) p[j] = level->v[i+j];

  enablevar(level, Ind("alpha"));
  enablevar(level, Ind("betax"));

  forallpoints(level, i) {

    mwblpsiofxyz(x[i], y[i], z[i], 
		 p[0]+i, p[1]+i, p[2]+i, p[3]+i, p[4]+i,
		 p[5]+i, p[6]+i, p[7]+i, p[8]+i, p[9]+i);

    g11[i] = g22[i] = g33[i] = 1.0;
    g12[i] = g13[i] = g23[i] = 0.0;

    K11[i] = K22[i] = K33[i] = 0.0;
    K12[i] = K13[i] = K23[i] = 0.0;

    
    if (fisheye) {
      indextrans_ddsym(x[i], y[i], z[i],
		       g11+i, g12+i, g13+i, g22+i, g23+i, g33+i);
    }
  }

}




/* transform with scalar function 
     u = v * psi^p
*/
void ConfTrans(tL *level, char *uname, char *vname, char *psiname, 
	       double power)
{
  int iu = Ind(uname);
  int iv = Ind(vname);
  int n = VarNComponents(iu);
  double *psi = Ptr(level, psiname);
  double **var = level->v;
  double c, *v;
  int i, j;

  if (n != VarNComponents(iv))
    errorexit("ConfTrans: the given tensors need same number of components");

  forallpoints(level, i) {
    c = pow(psi[i], power);
    for (j = 0; j < n; j++) 
      var[iu+j][i] = c * var[iv+j][i];
  }
}




/* given conformal data, set bam's adm data */
void ConfToPhys(tL *level, char *use_u)
{
  double *psi = Ptr(level, "adm_psi");
  int ig = Ind("adm_gxx");
  int iK = Ind("adm_Kxx");
  double *g, *K, *u, phi;
  int i, j;

  /* set metric and extrinsic curvature */
  if (*use_u) {
    u = Ptr(level, use_u);
    for (j = 0; j < 6; j++) {
      g = level->v[ig+j];
      K = level->v[iK+j];
      forallpoints(level, i) {
	phi = psi[i] + u[i];
	g[i] *= pow(phi/psi[i], 4.0);
	K[i] *= pow(phi, -2.0);
      }
    } 
  } else {
    for (j = 0; j < 6; j++) {
      K = level->v[iK+j];
      forallpoints(level, i)
	K[i] *= pow(psi[i], -2.0);
    }
  }

  /* hack to compare with BAM_Elliptic in Cactus */
  if (0 && *use_u) {
    forallpoints(level, i) u[i] += 1;
  }
}





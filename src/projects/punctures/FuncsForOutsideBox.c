/* FuncsForOutsideBox.c  */
/* Wolfgang Tichy 3/2003 */
/* Functions to integrate over with OutsideBox,
   in order to get masses outside the gird box */

#include "bam.h"
#include "punctures.h"




/* coeffs. of u and alpha and v */
static double C_u, C_alpha, C_v;
static double lapse_atpunc;

/* Set the coeffs. of u and alpha and v */
void SetFuncCoeffs(double c_u, double c_alpha, double c_v)
{
 C_u=c_u;
 C_alpha=c_alpha;
 C_v=c_v;
 lapse_atpunc=Getd("punctures_lapse_at_puncture");
}


/* RHS of Laplace u :  Laplace u = RHS_of_Lapl_punctures_u  */
double RHS_of_Lapl_punctures_u(double x, double y, double z)
{
 double a,ooa,b,u;
 double r = sqrt(x*x + y*y + z*z);   /* distance from (0,0,0)     */

 if (dequal(r, 0.0)) r = EPSZERO;
 
 ooa=OneOver_a(x, y, z); 
 a=1.0/ooa;                       
                       
 u=C_u/r;
 b=(1.0/8.0) * pow(a,7) * BY_KKofxyz(x, y, z); 
 return -b * pow( 1 + a*(1+u), -7); 
}


/* RHS of Laplace alphapsi :  Laplace v = RHS_of_Lapl_punctures_v  */
double RHS_of_Lapl_punctures_v(double x, double y, double z)
{
 double ooa, v, psi, u;
 double r = sqrt(x*x + y*y + z*z);   /* distance from (0,0,0)	  */

 if (dequal(r, 0.0)) r = EPSZERO;
 
 v=C_v/r;

 ooa = OneOver_a(x, y, z);
 u = C_u/r;
 psi = 1 + ooa + u;

 return (7.0/8.0) * (1.0+lapse_atpunc*ooa+v) * pow(psi,-8) * BY_KKofxyz(x, y, z);
}


/* compute the puncture 1/a from Bernd's puncture paper */ 
double OneOver_a(double x, double y, double z)
{
 double ooa;
 double xc,yc,zc, rc;                /* vector from BH to (x,y,z) */
 int n;
 
 ooa=0.0;
 for (n = 0; MBL[n] != 0.0 && n < N_MWBL; n++)
 {
   xc = x-CBL[n][0];
   yc = y-CBL[n][1];
   zc = z-CBL[n][2];
   rc  = sqrt( xc*xc + yc*yc + zc*zc);
   if (dequal(rc, 0.0)) rc = EPSZERO;
   ooa += MBL[n] / (2.0*rc);
 }
 return ooa;
}


/* RHS of Laplace alpha :  Laplace alpha = RHS_of_Lapl_alpha  */
double RHS_of_Lapl_alpha(double x, double y, double z)
{
 double alpha, ooa, psi;
 double r = sqrt(x*x + y*y + z*z);   /* distance from (0,0,0)	  */
 double xc,yc,zc, rc;                /* vector from BH to (x,y,z) */
 double nxc,nyc,nzc;
 double dooax,dooay,dooaz,  dpsix,dpsiy,dpsiz,  dalphax,dalphay,dalphaz;
 double RH,Derivs;
 int n;
 
 if (dequal(r, 0.0)) r = EPSZERO;
 
 ooa=0.0;
 dooax=0.0;
 dooay=0.0;
 dooaz=0.0;
 for (n = 0; MBL[n] != 0.0 && n < N_MWBL; n++)
 {
   xc = x-CBL[n][0];
   yc = y-CBL[n][1];
   zc = z-CBL[n][2];
   rc  = sqrt( xc*xc + yc*yc + zc*zc);
   if (dequal(rc, 0.0)) rc = EPSZERO;
   
   nxc = xc/rc;
   nyc = yc/rc;
   nzc = zc/rc;
   
   dooax += -( MBL[n] / (2.0*rc*rc) )*nxc;
   dooay += -( MBL[n] / (2.0*rc*rc) )*nyc;
   dooaz += -( MBL[n] / (2.0*rc*rc) )*nzc;   
   
   ooa += MBL[n] / (2.0*rc);
 }
 psi = 1 + ooa + C_u/r;
 dpsix = dooax - ( C_u/(r*r) )*(x/r);
 dpsiy = dooay - ( C_u/(r*r) )*(y/r);
 dpsiz = dooaz - ( C_u/(r*r) )*(z/r);
 
  
 alpha = C_alpha/r;
 dalphax = - ( C_alpha/(r*r) )*(x/r);
 dalphay = - ( C_alpha/(r*r) )*(y/r);
 dalphaz = - ( C_alpha/(r*r) )*(z/r);
   
 RH = pow(psi,-8) * alpha * BY_KKofxyz(x, y, z);
 Derivs = -(2.0/psi)*( dpsix*dalphax + dpsiy*dalphay + dpsiz*dalphaz );
 
 /* return RH; */
 return RH + Derivs;
}


/* generic_Gauge.c */
/* mth,  10/2010 */


#include "bam.h"
#include "Gauge.h"











void    generic_Gauge_boost(generic_Gauge_VARS)
{
    static double boost[4];
    
    
    if (*firstcall) {
        boost[1] = Getd("boost_x");
        boost[2] = Getd("boost_y");
        boost[3] = Getd("boost_z");
        *firstcall = 0;
    }
    
    *ralpha     = boost[1]*dadx + boost[2]*dady + boost[3]*dadz;
    
    *rbetax     = boost[1]*dbxdx + boost[2]*dbxdy + boost[3]*dbxdz;
    *rbetay     = boost[1]*dbydx + boost[2]*dbydy + boost[3]*dbydz;
    *rbetaz     = boost[1]*dbzdx + boost[2]*dbzdy + boost[3]*dbzdz;
    
    *rBx        = 0.;
    *rBy        = 0.;
    *rBz        = 0.;
    
}

void    generic_Gauge_onepluslog(generic_Gauge_VARS)
{
    static double eta,betaf,xi;
    
    if (*firstcall) {
        eta = Getd("Gauge_eta");
        xi    = Getd("Gauge_xi");
        betaf = Getd("Gauge_betaf");
        *firstcall = 0;
    }
    
    *ralpha     = advalpha - 2.*alpha*K;
    
    *rbetax     = advbetax + betaf*Bx;
    *rbetay     = advbetay + betaf*By;
    *rbetaz     = advbetaz + betaf*Bz;
    
    *rBx        = advBx    + xi*(rGx - advGx) - eta*Bx;
    *rBy        = advBy    + xi*(rGy - advGy) - eta*By;
    *rBz        = advBz    + xi*(rGz - advGz) - eta*Bz;
    
}

void    generic_Gauge_onepluslog_dookie(generic_Gauge_VARS)
{
    static double eta,xi;
    
    if (*firstcall) {
        eta = Getd("Gauge_eta");
        xi  = Getd("Gauge_xi");
        *firstcall = 0;
    }
    
    *ralpha     = advalpha - 2.*alpha*K;
    
    *rbetax     = advbetax + xi*Gx - eta*betax;
    *rbetay     = advbetay + xi*Gy - eta*betay;
    *rbetaz     = advbetaz + xi*Gz - eta*betaz;
    
    *rBx        = 0.;
    *rBy        = 0.;
    *rBz        = 0.;
    
}

void    generic_Gauge_harmonic(generic_Gauge_VARS)
{
    static double eta,xi,betaf;
    
    if (*firstcall) {
        eta   = Getd("Gauge_eta");
        xi    = Getd("Gauge_xi");
        betaf = Getd("Gauge_betaf");
        *firstcall = 0;
    }
    
    *ralpha     = advalpha - 2.*alpha*alpha*K;
    
    *rbetax     = advbetax + betaf*Bx;
    *rbetay     = advbetay + betaf*By;
    *rbetaz     = advbetaz + betaf*Bz;
    
    *rBx        = advBx    + alpha*alpha*xi*(rGx - advGx) - eta*Bx;
    *rBy        = advBy    + alpha*alpha*xi*(rGy - advGy) - eta*By;
    *rBz        = advBz    + alpha*alpha*xi*(rGz - advGz) - eta*Bz;
    
}






/* set the different gauge function which are directly 
   called inside bssn_rhs_noGauge.c */
void    set_generic_Gauge()
{
    if (Getv("Gauge", "boost"))
        generic_Gauge = generic_Gauge_boost;
    else if (Getv("Gauge", "boost_all"))
        generic_Gauge = generic_Gauge_boost;
    else if (Getv("Gauge", "1+log"))
        generic_Gauge = generic_Gauge_onepluslog;
    else if (Getv("Gauge", "1+log_dookie"))
        generic_Gauge = generic_Gauge_onepluslog_dookie;
    else if (Getv("Gauge", "harmonic"))
        generic_Gauge = generic_Gauge_harmonic;
    else
        errorexit("use a gauge!");
}




















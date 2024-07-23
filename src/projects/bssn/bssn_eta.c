/* bssn.c */
/* mth 03/12 */

#include "bam.h"
#include "bssn.h"



int etaswitch;
double m1,m2,M;
double xp1,yp1,zp1,xp2,yp2,zp2, Sep12;
double width,width1,width2 ,eta_power,eta_power0, et0,et1,epunc,einf;

void bssn_eta_init(tG* grid)
{
  if (!Getv("bssn_use_eta","yes")) return;
  
  etaswitch = Geti("bssn_eta_useswitch");
  m1 = Getd("mass1");
  m2 = Getd("mass2");
  M = m1 + m2;
  xp1 = grid->puncpos[0][0];
  yp1 = grid->puncpos[0][1];
  zp1 = grid->puncpos[0][2];
  xp2 = grid->puncpos[1][0];
  yp2 = grid->puncpos[1][1];
  zp2 = grid->puncpos[1][2];
  Sep12 = sqrt((xp1-xp2)*(xp1-xp2) +
      (yp1-yp2)*(yp1-yp2) + (zp1-zp2)*(zp1-zp2));
  width   = Getd("bssn_eta_width");
  width1  = Getd("bssn_eta_width1");
  width2  = Getd("bssn_eta_width2");
  eta_power   = Getd("bssn_eta_rpower");
  eta_power0  = Getd("bssn_eta_rpower0");
  et0     = Getd("bssn_eta_eta0");
  et1     = Getd("bssn_eta_eta1");
  epunc   = Getd("bssn_eta_epunc");
  einf    = Getd("bssn_eta_einf");
  
}

double bssn_eta_set(double xp,double yp,double zp, double absdpsim2, double psim2)
{
  double shiftDr;
  double R0 = 1.31;
  
  /* gauss */
  if (etaswitch == 5) { 

    if (Sep12 < .0001) {
      shiftDr = 2/M;
    } else {
      shiftDr = 2/M +
                (1/m1-2/M) * exp(- width *
                ((xp1-xp)*(xp1-xp)+(yp1-yp)*(yp1-yp)+(zp1-zp)*(zp1-zp))/(Sep12 * Sep12))
                +
                (1/m2-2/M)*exp(- width *
                ((xp2-xp)*(xp2-xp)+(yp2-yp)*(yp2-yp)+(zp2-zp)*(zp2-zp))/(Sep12 * Sep12));
    }

  /* rational fct*/
  } else if (etaswitch == 6) {

    if (Sep12 < .0001) {
      shiftDr = 2/M;
    } else {
      shiftDr = 2/M + (1/m1-2/M)/(1 + width1 *
                pow((xp1-xp)*(xp1-xp)+(yp1-yp)*(yp1-yp)+(zp1-zp)*(zp1-zp), eta_power)/(Sep12 * Sep12)) 
                + (1/m2-2/M)/(1 + width2 *
                pow((xp2-xp)*(xp2-xp)+(yp2-yp)*(yp2-yp)+(zp2-zp)*(zp2-zp), eta_power)/(Sep12 * Sep12));
    }

  /* rational fct, but proper scale in width; make this an extra switch for backwards compatibility */
  } else if (etaswitch == 8) {

    if (Sep12 < .0001) {
      shiftDr = 2/M;
    } else {
      double q = m2/m1;
      double Minv = 1/M;
      double Sep12invSq = 1/(Sep12 * Sep12);
      shiftDr = 2*Minv + (1/m1-2*Minv)/(1 + width1 * q *
                pow((xp1-xp)*(xp1-xp)+(yp1-yp)*(yp1-yp)+(zp1-zp)*(zp1-zp), eta_power) * Sep12invSq) 
                + (1/m2-2*Minv)/(1 + width2 * q *
                pow((xp2-xp)*(xp2-xp)+(yp2-yp)*(yp2-yp)+(zp2-zp)*(zp2-zp), eta_power) * Sep12invSq);
    }

  /* rational fct, but proper scale in width; make this an extra switch for backwards compatibility */
  } else if (etaswitch == 9) {
    
    if (Sep12 < .0001) {
      shiftDr = et0/M;
    } else {
      double Minv = 1/M;
      double Sep12invSq = 1/(Sep12 * Sep12);
      double etaiso = et1/et0 + (1-et1/et0)/(1+pow((xp*xp+yp*yp+zp*zp)/((Sep12+width)*(Sep12+width)),eta_power0));
      double etabase = et0*(Minv + (0.5/m1-Minv)/(1 + width1 *
                       pow((xp1-xp)*(xp1-xp)+(yp1-yp)*(yp1-yp)+(zp1-zp)*(zp1-zp), eta_power) * Sep12invSq)
                       + (0.5/m2-Minv)/(1 + width2 *
                       pow((xp2-xp)*(xp2-xp)+(yp2-yp)*(yp2-yp)+(zp2-zp)*(zp2-zp), eta_power) * Sep12invSq));
      shiftDr =  etabase * etaiso;
    }
  } else {
    
    shiftDr = R0*absdpsim2/pow(1 - psim2,2);
  }
  
  return shiftDr;
}









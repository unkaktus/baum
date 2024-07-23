/* scale_shift.c */
/* Wolfgang Tichy 12/2003 */

#include "bam.h"
#include "Gauge.h"

/* global vars containing lapse at bh on x-, y-, and z-axis */
extern double alpha_bh_x1;
extern double alpha_bh_y1;
extern double alpha_bh_z1;
extern double alpha_bh_x2;
extern double alpha_bh_y2;
extern double alpha_bh_z2;



/***************************************************************************/
/* try something simple: */
/* scale radial shift in order to keep alpha constant near bh */
int scale_radial_shift(tG *g)
{
  double shift_scale = Getd("corotation_shift_scale");
  double bh_lapse    = Getd("corotation_bh_lapse");

  double fac1x = bh_lapse / alpha_bh_x1;
  double fac1y = bh_lapse / alpha_bh_y1;
  double fac1z = bh_lapse / alpha_bh_z1;
  double fac1x2= fac1x*fac1x;
  double fac1y2= fac1y*fac1y;
  double fac1z2= fac1z*fac1z;

  double fac2x = bh_lapse / alpha_bh_x2;
  double fac2y = bh_lapse / alpha_bh_y2;
  double fac2z = bh_lapse / alpha_bh_z2;
  double fac2x2= fac2x*fac2x;
  double fac2y2= fac2y*fac2y;
  double fac2z2= fac2z*fac2z;

  double x1 = Getd("bhx1");
  double y1 = Getd("bhy1");
  double z1 = Getd("bhz1");
  double x2 = Getd("bhx2");
  double y2 = Getd("bhy2");
  double z2 = Getd("bhz2");
  

  int i, l;

  for (l = g->lmin; l <= g->lmax; l++)
  {
    tL *level = g->level[l];
    double *xp = Ptr(level, "x");
    double *yp = Ptr(level, "y");
    double *zp = Ptr(level, "z");
    double *betax = Ptr(level, "betax");
    double *betay = Ptr(level, "betay");
    double *betaz = Ptr(level, "betaz");
    double *betarotx = Ptr(level, "betarotx");
    double *betaroty = Ptr(level, "betaroty");
    double *betarotz = Ptr(level, "betarotz");
    double dbetx, dbety, dbetz;

    if(0) printf("scaling shift on level %d\n", level->l);

    forallpoints(level, i)
    {
      double r1x = xp[i]-x1;
      double r1y = yp[i]-y1;
      double r1z = zp[i]-z1;
      double r1x2= r1x*r1x;
      double r1y2= r1y*r1y;
      double r1z2= r1z*r1z;
      double r1  = sqrt(r1x2 + r1y2 + r1z2);

      double r2x = xp[i]-x2;
      double r2y = yp[i]-y2;
      double r2z = zp[i]-z2;
      double r2x2= r2x*r2x;
      double r2y2= r2y*r2y;
      double r2z2= r2z*r2z;
      double r2  = sqrt(r2x2 + r2y2 + r2z2);
      double fac1, fac2, fac;

      if(r1>0)
        fac1 = sqrt(fac1x2*r1x2 + fac1y2*r1y2 + fac1z2*r1z2)/r1;
      else
        fac1 = 1.0;

      if(r2>0)
        fac2 = sqrt(fac2x2*r2x2 + fac2y2*r2y2 + fac2z2*r2z2)/r2;
      else
   	fac2 = 1.0;
      
      fac1 = 1.0 + shift_scale*(fac1-1.0);
      fac2 = 1.0 + shift_scale*(fac2-1.0);            
      fac  = ( fac1*r2 + fac2*r1 )/(r1 + r2);

      dbetx = betax[i] - betarotx[i];
      dbety = betay[i] - betaroty[i];
      dbetz = betaz[i] - betarotz[i];
      
      betax[i] = fac*dbetx + betarotx[i];
      betay[i] = fac*dbety + betaroty[i];
      betaz[i] = fac*dbetz + betarotz[i];
    }
  } /* end for all levels */
  return 0;
}



/* shape of radial shift correction 
   used in addto_radial_shift and adjust_radial_shift */
double radial_shiftshape(double r, double m)
{
 return ( 1.0 - Attenuation01(r/(5.0*m), 2.0, 0.5) ) * 4.0 *
      	( 0.3 + r/m )/pow( r/m + 2, 3);
}


/* increase (decrease) the radial shift, if the lapse is too small (large). */
int addto_radial_shift(tG *g)
{
  double shift_factor = Getd("corotation_shift_additionFactor");
  double bh_lapse     = Getd("corotation_bh_lapse");

  double dif1x = bh_lapse - alpha_bh_x1;
  double dif1y = bh_lapse - alpha_bh_y1;
  double dif1z = bh_lapse - alpha_bh_z1;

  double dif2x = bh_lapse - alpha_bh_x2;
  double dif2y = bh_lapse - alpha_bh_y2;
  double dif2z = bh_lapse - alpha_bh_z2;

  double x1 = Getd("bhx1");
  double y1 = Getd("bhy1");
  double z1 = Getd("bhz1");
  double x2 = Getd("bhx2");
  double y2 = Getd("bhy2");
  double z2 = Getd("bhz2");
  double m1 = Getd("bhmass1");
  double m2 = Getd("bhmass2");
  int i, l;

  for (l = g->lmin; l <= g->lmax; l++)
  {
    tL *level = g->level[l];
    double *xp = Ptr(level, "x");
    double *yp = Ptr(level, "y");
    double *zp = Ptr(level, "z");
    double *betax = Ptr(level, "betax");
    double *betay = Ptr(level, "betay");
    double *betaz = Ptr(level, "betaz");
    double *betarotx = Ptr(level, "betarotx");
    double *betaroty = Ptr(level, "betaroty");
    double *betarotz = Ptr(level, "betarotz");
    double dbetx, dbety, dbetz;

    if(0) printf("scaling shift on level %d\n", level->l);

    forallpoints(level, i)
    {
      double r1x = xp[i]-x1;
      double r1y = yp[i]-y1;
      double r1z = zp[i]-z1;
      double r1x2= r1x*r1x;
      double r1y2= r1y*r1y;
      double r1z2= r1z*r1z;
      double r1r1= r1x2 + r1y2 + r1z2;
      double r1  = sqrt(r1r1);

      double r2x = xp[i]-x2;
      double r2y = yp[i]-y2;
      double r2z = zp[i]-z2;
      double r2x2= r2x*r2x;
      double r2y2= r2y*r2y;
      double r2z2= r2z*r2z;
      double r2r2= r2x2 + r2y2 + r2z2;
      double r2  = sqrt(r2r2);
      double dif1, dif2, shape1, shape2;

      if(r1>0)
        dif1 = (dif1x*r1x2 + dif1y*r1y2 + dif1z*r1z2)/(r1r1);
      else
        dif1 = 0.0;

      if(r2>0)
        dif2 = (dif2x*r2x2 + dif2y*r2y2 + dif2z*r2z2)/(r2r2);
      else
   	dif2 = 0.0;
      
      dif1 = shift_factor * dif1;
      dif2 = shift_factor * dif2;            

      shape1 = dif1 * radial_shiftshape(r1, m1);
      shape2 = dif2 * radial_shiftshape(r2, m2);

      dbetx = shape1 * r1x/r1 + shape2 * r2x/r2;
      dbety = shape1 * r1y/r1 + shape2 * r2y/r2;
      dbetz = shape1 * r1z/r1 + shape2 * r2z/r2;
      
      betarotx[i] += dbetx;
      betaroty[i] += dbety;
      betarotz[i] += dbetz;
       
      betax[i] += dbetx;
      betay[i] += dbety;
      betaz[i] += dbetz;
    }
  } /* end for all levels */
  return 0;
}


/* adjust the radial shift, i.e.                                            */
/* increase (decrease) the radial shift, if the lapse is too small (large). */
int adjust_radial_shift(tG *g)
{
  double shift_Adjust = Getd("corotation_shift_radialAdjustment");
  double bh_lapse     = Getd("corotation_bh_lapse");

  double x1 = Getd("bhx1");
  double y1 = Getd("bhy1");
  double z1 = Getd("bhz1");
  double x2 = Getd("bhx2");
  double y2 = Getd("bhy2");
  double z2 = Getd("bhz2");
  double m1 = Getd("bhmass1");
  double m2 = Getd("bhmass2");
  double lapse_radius = Getd("corotation_lapse_radius");
  double xoffset = Getd("commotion_offset_x");
  double yoffset = Getd("commotion_offset_y");
  tL *level_min = g->level[g->lmin];
  tL *level_max = g->level[g->lmax];
  tL *lapse_level1;
  tL *lapse_level2;
  double dx = level_max->dx;
  double dy = level_max->dy;
  double dz = level_max->dz;
  double h = pow(dx*dy*dz, 0.33333333333333333333333333333333);
  double r_ah1, r_ah2;
  int pr = 0;
  int i, l;

  /* center of spheres where we look at lapse */
  /* Note: this works only if we have equal mass puncs on y-axis */
  if(y1>0)
  {
    x1 += xoffset;    y1 += yoffset;
    x2 -= xoffset;    y2 -= yoffset;
  }
  else
  {
    x1 -= xoffset;    y1 -= yoffset;
    x2 += xoffset;    y2 += yoffset;
  }
  
  /* radius where we look at lapse */
  r_ah1 = lapse_radius*m1 + 2.0*h;
  r_ah2 = lapse_radius*m2 + 2.0*h;
  
  /* level on which this radius lies */
  lapse_level1 = NULL;
  lapse_level2 = NULL;
  for (l = g->lmin; l <= g->lmax; l++)
  {
    tL *level = g->level[l];
    double ymin = level->bbox[2] + 2.0*level->dy;
    double ymax = level->bbox[3] - 2.0*level->dy;
    if( (y1 + r_ah1 < ymax)&&(y1 - r_ah1 > ymin) )
      lapse_level1 = level;
    if( (y2 + r_ah2 < ymax)&&(y2 - r_ah2 > ymin) )
      lapse_level2 = level;
  }
  if( lapse_level1==NULL)  lapse_level1 = lapse_level2;
  if( lapse_level2==NULL)  lapse_level2 = lapse_level1;
  if( lapse_level1==NULL && lapse_level2==NULL )
  {
    lapse_level1 = level_min;
    lapse_level2 = level_min;
  }
  printf("adjust_radial_shift:\n");
  printf(" beta will be adjusted using alpha on spheres with\n"
         " r1=%e (on level %d)    r2=%e (on level %d)\n",
         r_ah1, lapse_level1->l, r_ah2, lapse_level2->l);

  /* loop over all levels */
  for (l = g->lmin; l <= g->lmax; l++)
  {
    tL *level = g->level[l];
    double *xp = Ptr(level, "x");
    double *yp = Ptr(level, "y");
    double *zp = Ptr(level, "z");
    double *betax = Ptr(level, "betax");
    double *betay = Ptr(level, "betay");
    double *betaz = Ptr(level, "betaz");
    double *betarotx = Ptr(level, "betarotx");
    double *betaroty = Ptr(level, "betaroty");
    double *betarotz = Ptr(level, "betarotz");
    double dbetx, dbety, dbetz;
    double alpha1, alpha2;

    if(0) printf("adjusting shift on level %d\n", level->l);

    forallpoints(level, i)
    {
      double r1x = xp[i]-x1;
      double r1y = yp[i]-y1;
      double r1z = zp[i]-z1;
      double r1x2= r1x*r1x;
      double r1y2= r1y*r1y;
      double r1z2= r1z*r1z;
      double r1  = sqrt(r1x2 + r1y2 + r1z2);

      double r2x = xp[i]-x2;
      double r2y = yp[i]-y2;
      double r2z = zp[i]-z2;
      double r2x2= r2x*r2x;
      double r2y2= r2y*r2y;
      double r2z2= r2z*r2z;
      double r2  = sqrt(r2x2 + r2y2 + r2z2);

      double dif1, dif2, shape1, shape2;
      double costh1, sinth1, cosph1, sinph1;
      double costh2, sinth2, cosph2, sinph2;
      double s1x,s1y,s1z, s2x,s2y,s2z;   /* coords of points on spheres */
      
      if(r1>0.0)
      {
        costh1 = r1z/r1;
        sinth1 = sin(acos(costh1));
        if(sinth1!=0.0)
        {
          cosph1 = r1x/(r1*sinth1);
          sinph1 = r1y/(r1*sinth1);
        }
        else
        {
          cosph1 = 1.0;
          sinph1 = 0.0;
        }
      }
      else
      { costh1=0.0;  sinth1=1.0;  cosph1=0.0;  sinph1=1.0;}
      
      if(r2>0.0)
      {
        costh2 = r2z/r2;
        sinth2 = sin(acos(costh2));
        if(sinth2!=0.0)
        {
          cosph2 = r2x/(r2*sinth2);
          sinph2 = r2y/(r2*sinth2);
        }
        else
        {
          cosph2 = 1.0;
          sinph2 = 0.0;
        }
      }
      else
      { costh2=0.0;  sinth2=1.0;  cosph2=0.0;  sinph2=1.0;}
      
      /* find alpha on sphere of radius r_ah, around each BH */
      s1x = x1 + r_ah1*sinth1*cosph1;
      s1y = y1 + r_ah1*sinth1*sinph1;
      s1z = z1 + r_ah1*costh1;
      s2x = x2 + r_ah2*sinth2*cosph2;
      s2y = y2 + r_ah2*sinth2*sinph2;
      s2z = z2 + r_ah2*costh2;
      
      if(pr) printf("1: %f %f %f\n", s1x, s1y, s1z);
      alpha1=interpolate_xyz_scalar_withSym(lapse_level1, s1x, s1y, s1z,
                                            Ind("alpha"), pr);

      if(pr) printf("2: %f %f %f\n", s2x, s2y, s2z);
      alpha2=interpolate_xyz_scalar_withSym(lapse_level2, s2x, s2y, s2z,
                                            Ind("alpha"), pr);
      dif1 = bh_lapse - alpha1;
      dif2 = bh_lapse - alpha2;
            
      dif1 = shift_Adjust * dif1;
      dif2 = shift_Adjust * dif2;            

      /* shape of radial shift */
      shape1 = dif1 * radial_shiftshape(r1, m1);
      shape2 = dif2 * radial_shiftshape(r2, m2);

      dbetx = shape1 * r1x/r1 + shape2 * r2x/r2;
      dbety = shape1 * r1y/r1 + shape2 * r2y/r2;
      dbetz = shape1 * r1z/r1 + shape2 * r2z/r2;
      
      betarotx[i] += dbetx;
      betaroty[i] += dbety;
      betarotz[i] += dbetz;
       
      betax[i] += dbetx;
      betay[i] += dbety;
      betaz[i] += dbetz;
      
      if(0) printf("x  %f %f %f\n", xp[i], yp[i], zp[i]);
      if(0) printf("r1 %f %f %f\n", r1x, r1y, r1z);
      if(0) printf("r1=%f\n", r1);
      if(0) printf("s1 %f %f %f\n", s1x, s1y, s1z);
      if(0) printf("costh1=%f th1=%f sinth1=%f\n", 
             costh1, acos(costh1), sin(acos(costh1)));
      if(0) printf("cosph1=%f ph1=%f sinph1=%f\n",
             cosph1, acos(cosph1), sin(acos(cosph1)));
      if(0) printf("r1y/(r1*sinth1)=%f\n", r1y/(r1*sinth1));

      if(0) printf("alpha1=%e  dif1=%e  dbetx=%e\n", alpha1, dif1, dbety);
      if(0) exit(0);
    }
  } /* end for all levels */
  return 0;
}

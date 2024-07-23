/* SphericalHarmonics_precomp.c */
/* mth, 13.6.2012 */
/* Produced with Mathematica */

#include "bam.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))
#define Cal(x,y,z) ((x)?(y):(z))
#define Pi         PI
#define Cos(x)     cos(x)
#define Sin(x)     sin(x)
#define Sqrt(x)    sqrt(x)
#define Csc(x)     (1./sin(x))







void SphericalHarmonicYprecomp(double *rY, double *iY, int l, int m, double phi, double theta)
{

if (l==2) {

  if (m==-2) {
    *rY  = (Sqrt(15/(2.*Pi))*Cos(2*phi)*Power(Sin(theta),2))/4.;
    *iY  = -(Sqrt(15/(2.*Pi))*Sin(2*phi)*Power(Sin(theta),2))/4.;
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(15/(2.*Pi))*Cos(phi)*Cos(theta)*Sin(theta))/2.;
    *iY  = -(Sqrt(15/(2.*Pi))*Cos(theta)*Sin(phi)*Sin(theta))/2.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(5/Pi)*(1 + 3*Cos(2*theta)))/8.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = -(Sqrt(15/(2.*Pi))*Cos(phi)*Cos(theta)*Sin(theta))/2.;
    *iY  = -(Sqrt(15/(2.*Pi))*Cos(theta)*Sin(phi)*Sin(theta))/2.;
    return;
  } else if (m==2) {
    *rY  = (Sqrt(15/(2.*Pi))*Cos(2*phi)*Power(Sin(theta),2))/4.;
    *iY  = (Sqrt(15/(2.*Pi))*Sin(2*phi)*Power(Sin(theta),2))/4.;
    return;
  }

} else if (l==3) {

  if (m==-3) {
    *rY  = (Sqrt(35/Pi)*Cos(3*phi)*Power(Sin(theta),3))/8.;
    *iY  = -(Sqrt(35/Pi)*Sin(3*phi)*Power(Sin(theta),3))/8.;
    return;
  } else if (m==-2) {
    *rY  = (Sqrt(105/(2.*Pi))*Cos(2*phi)*Cos(theta)*Power(Sin(theta),2))/4.;
    *iY  = -(Sqrt(105/(2.*Pi))*Cos(theta)*Sin(2*phi)*Power(Sin(theta),2))/4.;
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(21/Pi)*Cos(phi)*(Sin(theta) + 5*Sin(3*theta)))/32.;
    *iY  = -(Sqrt(21/Pi)*Sin(phi)*(Sin(theta) + 5*Sin(3*theta)))/32.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(7/Pi)*(3*Cos(theta) + 5*Cos(3*theta)))/16.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = -(Sqrt(21/Pi)*Cos(phi)*(Sin(theta) + 5*Sin(3*theta)))/32.;
    *iY  = -(Sqrt(21/Pi)*Sin(phi)*(Sin(theta) + 5*Sin(3*theta)))/32.;
    return;
  } else if (m==2) {
    *rY  = (Sqrt(105/(2.*Pi))*Cos(2*phi)*Cos(theta)*Power(Sin(theta),2))/4.;
    *iY  = (Sqrt(105/(2.*Pi))*Cos(theta)*Sin(2*phi)*Power(Sin(theta),2))/4.;
    return;
  } else if (m==3) {
    *rY  = -(Sqrt(35/Pi)*Cos(3*phi)*Power(Sin(theta),3))/8.;
    *iY  = -(Sqrt(35/Pi)*Sin(3*phi)*Power(Sin(theta),3))/8.;
    return;
  }

} else if (l==4) {

  if (m==-4) {
    *rY  = (3*Sqrt(35/(2.*Pi))*Cos(4*phi)*Power(Sin(theta),4))/16.;
    *iY  = (-3*Sqrt(35/(2.*Pi))*Sin(4*phi)*Power(Sin(theta),4))/16.;
    return;
  } else if (m==-3) {
    *rY  = (3*Sqrt(35/Pi)*Cos(3*phi)*Cos(theta)*Power(Sin(theta),3))/8.;
    *iY  = (-3*Sqrt(35/Pi)*Cos(theta)*Sin(3*phi)*Power(Sin(theta),3))/8.;
    return;
  } else if (m==-2) {
    *rY  = (3*Sqrt(5/(2.*Pi))*Cos(2*phi)*(5 + 7*Cos(2*theta))*Power(Sin(theta),2))/16.;
    *iY  = (-3*Sqrt(5/(2.*Pi))*(5 + 7*Cos(2*theta))*Sin(2*phi)*Power(Sin(theta),2))/16.;
    return;
  } else if (m==-1) {
    *rY  = (3*Sqrt(5/Pi)*Cos(phi)*(2*Sin(2*theta) + 7*Sin(4*theta)))/64.;
    *iY  = (-3*Sqrt(5/Pi)*Sin(phi)*(2*Sin(2*theta) + 7*Sin(4*theta)))/64.;
    return;
  } else if (m==0) {
    *rY  = (3*(9 + 20*Cos(2*theta) + 35*Cos(4*theta)))/(128.*Sqrt(Pi));
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = (-3*Sqrt(5/Pi)*Cos(phi)*(2*Sin(2*theta) + 7*Sin(4*theta)))/64.;
    *iY  = (-3*Sqrt(5/Pi)*Sin(phi)*(2*Sin(2*theta) + 7*Sin(4*theta)))/64.;
    return;
  } else if (m==2) {
    *rY  = (3*Sqrt(5/(2.*Pi))*Cos(2*phi)*(5 + 7*Cos(2*theta))*Power(Sin(theta),2))/16.;
    *iY  = (3*Sqrt(5/(2.*Pi))*(5 + 7*Cos(2*theta))*Sin(2*phi)*Power(Sin(theta),2))/16.;
    return;
  } else if (m==3) {
    *rY  = (-3*Sqrt(35/Pi)*Cos(3*phi)*Cos(theta)*Power(Sin(theta),3))/8.;
    *iY  = (-3*Sqrt(35/Pi)*Cos(theta)*Sin(3*phi)*Power(Sin(theta),3))/8.;
    return;
  } else if (m==4) {
    *rY  = (3*Sqrt(35/(2.*Pi))*Cos(4*phi)*Power(Sin(theta),4))/16.;
    *iY  = (3*Sqrt(35/(2.*Pi))*Sin(4*phi)*Power(Sin(theta),4))/16.;
    return;
  }

} else if (l==5) {

  if (m==-5) {
    *rY  = (3*Sqrt(77/Pi)*Cos(5*phi)*Power(Sin(theta),5))/32.;
    *iY  = (-3*Sqrt(77/Pi)*Sin(5*phi)*Power(Sin(theta),5))/32.;
    return;
  } else if (m==-4) {
    *rY  = (3*Sqrt(385/(2.*Pi))*Cos(4*phi)*Cos(theta)*Power(Sin(theta),4))/16.;
    *iY  = (-3*Sqrt(385/(2.*Pi))*Cos(theta)*Sin(4*phi)*Power(Sin(theta),4))/16.;
    return;
  } else if (m==-3) {
    *rY  = (Sqrt(385/Pi)*Cos(3*phi)*(7 + 9*Cos(2*theta))*Power(Sin(theta),3))/64.;
    *iY  = -(Sqrt(385/Pi)*(7 + 9*Cos(2*theta))*Sin(3*phi)*Power(Sin(theta),3))/64.;
    return;
  } else if (m==-2) {
    *rY  = (Sqrt(1155/(2.*Pi))*Cos(2*phi)*(5*Cos(theta) + 3*Cos(3*theta))*Power(Sin(theta),2))/32.;
    *iY  = -(Sqrt(1155/(2.*Pi))*(5*Cos(theta) + 3*Cos(3*theta))*Sin(2*phi)*Power(Sin(theta),2))/32.;
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(165/(2.*Pi))*Cos(phi)*(2*Sin(theta) + 7*(Sin(3*theta) + 3*Sin(5*theta))))/256.;
    *iY  = -(Sqrt(165/(2.*Pi))*Sin(phi)*(2*Sin(theta) + 7*(Sin(3*theta) + 3*Sin(5*theta))))/256.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(11/Pi)*(30*Cos(theta) + 35*Cos(3*theta) + 63*Cos(5*theta)))/256.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = -(Sqrt(165/(2.*Pi))*Cos(phi)*(2*Sin(theta) + 7*(Sin(3*theta) + 3*Sin(5*theta))))/256.;
    *iY  = -(Sqrt(165/(2.*Pi))*Sin(phi)*(2*Sin(theta) + 7*(Sin(3*theta) + 3*Sin(5*theta))))/256.;
    return;
  } else if (m==2) {
    *rY  = (Sqrt(1155/(2.*Pi))*Cos(2*phi)*(5*Cos(theta) + 3*Cos(3*theta))*Power(Sin(theta),2))/32.;
    *iY  = (Sqrt(1155/(2.*Pi))*(5*Cos(theta) + 3*Cos(3*theta))*Sin(2*phi)*Power(Sin(theta),2))/32.;
    return;
  } else if (m==3) {
    *rY  = -(Sqrt(385/Pi)*Cos(3*phi)*(7 + 9*Cos(2*theta))*Power(Sin(theta),3))/64.;
    *iY  = -(Sqrt(385/Pi)*(7 + 9*Cos(2*theta))*Sin(3*phi)*Power(Sin(theta),3))/64.;
    return;
  } else if (m==4) {
    *rY  = (3*Sqrt(385/(2.*Pi))*Cos(4*phi)*Cos(theta)*Power(Sin(theta),4))/16.;
    *iY  = (3*Sqrt(385/(2.*Pi))*Cos(theta)*Sin(4*phi)*Power(Sin(theta),4))/16.;
    return;
  } else if (m==5) {
    *rY  = (-3*Sqrt(77/Pi)*Cos(5*phi)*Power(Sin(theta),5))/32.;
    *iY  = (-3*Sqrt(77/Pi)*Sin(5*phi)*Power(Sin(theta),5))/32.;
    return;
  }

} else if (l==6) {

  if (m==-6) {
    *rY  = (Sqrt(3003/Pi)*Cos(6*phi)*Power(Sin(theta),6))/64.;
    *iY  = -(Sqrt(3003/Pi)*Sin(6*phi)*Power(Sin(theta),6))/64.;
    return;
  } else if (m==-5) {
    *rY  = (3*Sqrt(1001/Pi)*Cos(5*phi)*Cos(theta)*Power(Sin(theta),5))/32.;
    *iY  = (-3*Sqrt(1001/Pi)*Cos(theta)*Sin(5*phi)*Power(Sin(theta),5))/32.;
    return;
  } else if (m==-4) {
    *rY  = (3*Sqrt(91/(2.*Pi))*Cos(4*phi)*(9 + 11*Cos(2*theta))*Power(Sin(theta),4))/64.;
    *iY  = (-3*Sqrt(91/(2.*Pi))*(9 + 11*Cos(2*theta))*Sin(4*phi)*Power(Sin(theta),4))/64.;
    return;
  } else if (m==-3) {
    *rY  = (Sqrt(1365/Pi)*Cos(3*phi)*(21*Cos(theta) + 11*Cos(3*theta))*Power(Sin(theta),3))/128.;
    *iY  = -(Sqrt(1365/Pi)*(21*Cos(theta) + 11*Cos(3*theta))*Sin(3*phi)*Power(Sin(theta),3))/128.;
    return;
  } else if (m==-2) {
    *rY  = (Sqrt(1365/Pi)*Cos(2*phi)*(35 + 60*Cos(2*theta) + 33*Cos(4*theta))*Power(Sin(theta),2))/512.;
    *iY  = -(Sqrt(1365/Pi)*(35 + 60*Cos(2*theta) + 33*Cos(4*theta))*Sin(2*phi)*Power(Sin(theta),2))/512.;
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(273/(2.*Pi))*Cos(phi)*(5*Sin(2*theta) + 12*Sin(4*theta) + 33*Sin(6*theta)))/512.;
    *iY  = -(Sqrt(273/(2.*Pi))*Sin(phi)*(5*Sin(2*theta) + 12*Sin(4*theta) + 33*Sin(6*theta)))/512.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(13/Pi)*(50 + 105*Cos(2*theta) + 126*Cos(4*theta) + 231*Cos(6*theta)))/1024.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = -(Sqrt(273/(2.*Pi))*Cos(phi)*(5*Sin(2*theta) + 12*Sin(4*theta) + 33*Sin(6*theta)))/512.;
    *iY  = -(Sqrt(273/(2.*Pi))*Sin(phi)*(5*Sin(2*theta) + 12*Sin(4*theta) + 33*Sin(6*theta)))/512.;
    return;
  } else if (m==2) {
    *rY  = (Sqrt(1365/Pi)*Cos(2*phi)*(35 + 60*Cos(2*theta) + 33*Cos(4*theta))*Power(Sin(theta),2))/512.;
    *iY  = (Sqrt(1365/Pi)*(35 + 60*Cos(2*theta) + 33*Cos(4*theta))*Sin(2*phi)*Power(Sin(theta),2))/512.;
    return;
  } else if (m==3) {
    *rY  = -(Sqrt(1365/Pi)*Cos(3*phi)*(21*Cos(theta) + 11*Cos(3*theta))*Power(Sin(theta),3))/128.;
    *iY  = -(Sqrt(1365/Pi)*(21*Cos(theta) + 11*Cos(3*theta))*Sin(3*phi)*Power(Sin(theta),3))/128.;
    return;
  } else if (m==4) {
    *rY  = (3*Sqrt(91/(2.*Pi))*Cos(4*phi)*(9 + 11*Cos(2*theta))*Power(Sin(theta),4))/64.;
    *iY  = (3*Sqrt(91/(2.*Pi))*(9 + 11*Cos(2*theta))*Sin(4*phi)*Power(Sin(theta),4))/64.;
    return;
  } else if (m==5) {
    *rY  = (-3*Sqrt(1001/Pi)*Cos(5*phi)*Cos(theta)*Power(Sin(theta),5))/32.;
    *iY  = (-3*Sqrt(1001/Pi)*Cos(theta)*Sin(5*phi)*Power(Sin(theta),5))/32.;
    return;
  } else if (m==6) {
    *rY  = (Sqrt(3003/Pi)*Cos(6*phi)*Power(Sin(theta),6))/64.;
    *iY  = (Sqrt(3003/Pi)*Sin(6*phi)*Power(Sin(theta),6))/64.;
    return;
  }

} else errorexit("increase lmax in mathematica script");


}




void spinweightedSphericalHarmonicprecomp(double *rY, double *iY, int l, int m, double phi, double theta)
{

if (l==2) {

  if (m==-2) {
    *rY  = (Sqrt(5/Pi)*Cos(2*phi)*Power(Sin(theta/2.),4))/2.;
    *iY  = -(Sqrt(5/Pi)*Cos(phi)*Sin(phi)*Power(Sin(theta/2.),4));
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(5/Pi)*Cos(phi)*Power(Sin(theta/2.),2)*Sin(theta))/2.;
    *iY  = -(Sqrt(5/Pi)*Sin(phi)*Power(Sin(theta/2.),2)*Sin(theta))/2.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(15/(2.*Pi))*Power(Sin(theta),2))/4.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = Sqrt(5/Pi)*Cos(phi)*Power(Cos(theta/2.),3)*Sin(theta/2.);
    *iY  = Sqrt(5/Pi)*Power(Cos(theta/2.),3)*Sin(phi)*Sin(theta/2.);
    return;
  } else if (m==2) {
    *rY  = (Sqrt(5/Pi)*Cos(2*phi)*Power(Cos(theta/2.),4))/2.;
    *iY  = Sqrt(5/Pi)*Cos(phi)*Power(Cos(theta/2.),4)*Sin(phi);
    return;
  }

} else if (l==3) {

  if (m==-3) {
    *rY  = (Sqrt(21/(2.*Pi))*Cos(3*phi)*Power(Sin(theta/2.),4)*Sin(theta))/2.;
    *iY  = -(Sqrt(21/(2.*Pi))*Sin(3*phi)*Power(Sin(theta/2.),4)*Sin(theta))/2.;
    return;
  } else if (m==-2) {
    *rY  = (Sqrt(7/Pi)*Cos(2*phi)*(2 + 3*Cos(theta))*Power(Sin(theta/2.),4))/2.;
    *iY  = -(Sqrt(7/Pi)*Cos(phi)*(2 + 3*Cos(theta))*Sin(phi)*Power(Sin(theta/2.),4));
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(35/(2.*Pi))*Cos(phi)*(5*Cos(theta/2.) + 3*Cos((3*theta)/2.))*Power(Sin(theta/2.),3))/4.;
    *iY  = -(Sqrt(35/(2.*Pi))*(5*Cos(theta/2.) + 3*Cos((3*theta)/2.))*Sin(phi)*Power(Sin(theta/2.),3))/4.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(105/(2.*Pi))*Cos(theta)*Power(Sin(theta),2))/4.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = (Sqrt(35/(2.*Pi))*Cos(phi)*Power(Cos(theta/2.),3)*(-5*Sin(theta/2.) + 3*Sin((3*theta)/2.)))/4.;
    *iY  = (Sqrt(35/(2.*Pi))*Power(Cos(theta/2.),3)*Sin(phi)*(-5*Sin(theta/2.) + 3*Sin((3*theta)/2.)))/4.;
    return;
  } else if (m==2) {
    *rY  = (Sqrt(7/Pi)*Cos(2*phi)*Power(Cos(theta/2.),4)*(-2 + 3*Cos(theta)))/2.;
    *iY  = Sqrt(7/Pi)*Cos(phi)*Power(Cos(theta/2.),4)*(-2 + 3*Cos(theta))*Sin(phi);
    return;
  } else if (m==3) {
    *rY  = -(Sqrt(21/(2.*Pi))*Cos(3*phi)*Power(Cos(theta/2.),5)*Sin(theta/2.));
    *iY  = -(Sqrt(21/(2.*Pi))*Power(Cos(theta/2.),5)*Sin(3*phi)*Sin(theta/2.));
    return;
  }

} else if (l==4) {

  if (m==-4) {
    *rY  = (3*Sqrt(7/Pi)*Cos(4*phi)*Power(Sin(theta/2.),4)*Power(Sin(theta),2))/4.;
    *iY  = (-3*Sqrt(7/Pi)*Sin(4*phi)*Power(Sin(theta/2.),4)*Power(Sin(theta),2))/4.;
    return;
  } else if (m==-3) {
    *rY  = 3*Sqrt(7/(2.*Pi))*Cos(3*phi)*(2*Cos(theta/2.) + Cos((3*theta)/2.))*Power(Sin(theta/2.),5);
    *iY  = -3*Sqrt(7/(2.*Pi))*(2*Cos(theta/2.) + Cos((3*theta)/2.))*Sin(3*phi)*Power(Sin(theta/2.),5);
    return;
  } else if (m==-2) {
    *rY  = (3*Cos(2*phi)*(9 + 14*Cos(theta) + 7*Cos(2*theta))*Power(Sin(theta/2.),4))/(4.*Sqrt(Pi));
    *iY  = (-3*(9 + 14*Cos(theta) + 7*Cos(2*theta))*Sin(2*phi)*Power(Sin(theta/2.),4))/(4.*Sqrt(Pi));
    return;
  } else if (m==-1) {
    *rY  = (3*Cos(phi)*(19*Cos(theta/2.) + 7*(2*Cos((3*theta)/2.) + Cos((5*theta)/2.)))*Power(Sin(theta/2.),3))/(4.*Sqrt(2*Pi));
    *iY  = (-3*(19*Cos(theta/2.) + 7*(2*Cos((3*theta)/2.) + Cos((5*theta)/2.)))*Sin(phi)*Power(Sin(theta/2.),3))/(4.*Sqrt(2*Pi));
    return;
  } else if (m==0) {
    *rY  = (3*Sqrt(5/(2.*Pi))*(5 + 7*Cos(2*theta))*Power(Sin(theta),2))/16.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = (3*Cos(phi)*Power(Cos(theta/2.),3)*(19*Sin(theta/2.) + 7*(-2*Sin((3*theta)/2.) + Sin((5*theta)/2.))))/(4.*Sqrt(2*Pi));
    *iY  = (3*Power(Cos(theta/2.),3)*Sin(phi)*(19*Sin(theta/2.) + 7*(-2*Sin((3*theta)/2.) + Sin((5*theta)/2.))))/(4.*Sqrt(2*Pi));
    return;
  } else if (m==2) {
    *rY  = (3*Cos(2*phi)*Power(Cos(theta/2.),4)*(9 - 14*Cos(theta) + 7*Cos(2*theta)))/(4.*Sqrt(Pi));
    *iY  = (3*Power(Cos(theta/2.),4)*(9 - 14*Cos(theta) + 7*Cos(2*theta))*Sin(2*phi))/(4.*Sqrt(Pi));
    return;
  } else if (m==3) {
    *rY  = 3*Sqrt(7/(2.*Pi))*Cos(3*phi)*Power(Cos(theta/2.),5)*(1 - 2*Cos(theta))*Sin(theta/2.);
    *iY  = 3*Sqrt(7/(2.*Pi))*Power(Cos(theta/2.),5)*(1 - 2*Cos(theta))*Sin(3*phi)*Sin(theta/2.);
    return;
  } else if (m==4) {
    *rY  = 3*Sqrt(7/Pi)*Cos(4*phi)*Power(Cos(theta/2.),6)*Power(Sin(theta/2.),2);
    *iY  = (3*Sqrt(7/Pi)*Power(Csc(theta/2.),4)*Sin(4*phi)*Power(Sin(theta),6))/64.;
    return;
  }

} else if (l==5) {

  if (m==-5) {
    *rY  = Sqrt(330/Pi)*Cos(5*phi)*Power(Cos(theta/2.),3)*Power(Sin(theta/2.),7);
    *iY  = -(Sqrt(330/Pi)*Power(Cos(theta/2.),3)*Sin(5*phi)*Power(Sin(theta/2.),7));
    return;
  } else if (m==-4) {
    *rY  = Sqrt(33/Pi)*Cos(4*phi)*Power(Cos(theta/2.),2)*(2 + 5*Cos(theta))*Power(Sin(theta/2.),6);
    *iY  = -(Sqrt(33/Pi)*Power(Cos(theta/2.),2)*(2 + 5*Cos(theta))*Sin(4*phi)*Power(Sin(theta/2.),6));
    return;
  } else if (m==-3) {
    *rY  = (Sqrt(33/(2.*Pi))*Cos(3*phi)*(58*Cos(theta/2.) + 39*Cos((3*theta)/2.) + 15*Cos((5*theta)/2.))*Power(Sin(theta/2.),5))/8.;
    *iY  = -(Sqrt(33/(2.*Pi))*(58*Cos(theta/2.) + 39*Cos((3*theta)/2.) + 15*Cos((5*theta)/2.))*Sin(3*phi)*Power(Sin(theta/2.),5))/8.;
    return;
  } else if (m==-2) {
    *rY  = (Sqrt(11/Pi)*Cos(2*phi)*(32 + 57*Cos(theta) + 36*Cos(2*theta) + 15*Cos(3*theta))*Power(Sin(theta/2.),4))/8.;
    *iY  = -(Sqrt(11/Pi)*(32 + 57*Cos(theta) + 36*Cos(2*theta) + 15*Cos(3*theta))*Sin(2*phi)*Power(Sin(theta/2.),4))/8.;
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(77/Pi)*Cos(phi)*(61*Cos(theta/2.) + 51*Cos((3*theta)/2.) + 33*Cos((5*theta)/2.) + 15*Cos((7*theta)/2.))*Power(Sin(theta/2.),3))/32.;
    *iY  = -(Sqrt(77/Pi)*(61*Cos(theta/2.) + 51*Cos((3*theta)/2.) + 33*Cos((5*theta)/2.) + 15*Cos((7*theta)/2.))*Sin(phi)*Power(Sin(theta/2.),3))/32.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(1155/(2.*Pi))*(5*Cos(theta) + 3*Cos(3*theta))*Power(Sin(theta),2))/32.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = (Sqrt(77/Pi)*Cos(phi)*Power(Cos(theta/2.),3)*(-61*Sin(theta/2.) + 51*Sin((3*theta)/2.) - 33*Sin((5*theta)/2.) + 15*Sin((7*theta)/2.)))/32.;
    *iY  = (Sqrt(77/Pi)*Power(Cos(theta/2.),3)*Sin(phi)*(-61*Sin(theta/2.) + 51*Sin((3*theta)/2.) - 33*Sin((5*theta)/2.) + 15*Sin((7*theta)/2.)))/32.;
    return;
  } else if (m==2) {
    *rY  = (Sqrt(11/Pi)*Cos(2*phi)*Power(Cos(theta/2.),4)*(-32 + 57*Cos(theta) - 36*Cos(2*theta) + 15*Cos(3*theta)))/8.;
    *iY  = (Sqrt(11/Pi)*Power(Cos(theta/2.),4)*(-32 + 57*Cos(theta) - 36*Cos(2*theta) + 15*Cos(3*theta))*Sin(2*phi))/8.;
    return;
  } else if (m==3) {
    *rY  = -(Sqrt(33/(2.*Pi))*Cos(3*phi)*Power(Cos(theta/2.),5)*(58*Sin(theta/2.) - 39*Sin((3*theta)/2.) + 15*Sin((5*theta)/2.)))/8.;
    *iY  = -(Sqrt(33/(2.*Pi))*Power(Cos(theta/2.),5)*Sin(3*phi)*(58*Sin(theta/2.) - 39*Sin((3*theta)/2.) + 15*Sin((5*theta)/2.)))/8.;
    return;
  } else if (m==4) {
    *rY  = Sqrt(33/Pi)*Cos(4*phi)*Power(Cos(theta/2.),6)*(-2 + 5*Cos(theta))*Power(Sin(theta/2.),2);
    *iY  = Sqrt(33/Pi)*Power(Cos(theta/2.),6)*(-2 + 5*Cos(theta))*Sin(4*phi)*Power(Sin(theta/2.),2);
    return;
  } else if (m==5) {
    *rY  = -(Sqrt(330/Pi)*Cos(5*phi)*Power(Cos(theta/2.),7)*Power(Sin(theta/2.),3));
    *iY  = -(Sqrt(330/Pi)*Power(Cos(theta/2.),7)*Sin(5*phi)*Power(Sin(theta/2.),3));
    return;
  }

} else if (l==6) {

  if (m==-6) {
    *rY  = (3*Sqrt(715/Pi)*Cos(6*phi)*Power(Sin(theta/2.),4)*Power(Sin(theta),4))/32.;
    *iY  = (-3*Sqrt(715/Pi)*Sin(6*phi)*Power(Sin(theta/2.),4)*Power(Sin(theta),4))/32.;
    return;
  } else if (m==-5) {
    *rY  = (Sqrt(2145/Pi)*Cos(5*phi)*Power(Cos(theta/2.),3)*(1 + 3*Cos(theta))*Power(Sin(theta/2.),7))/2.;
    *iY  = -(Sqrt(2145/Pi)*Power(Cos(theta/2.),3)*(1 + 3*Cos(theta))*Sin(5*phi)*Power(Sin(theta/2.),7))/2.;
    return;
  } else if (m==-4) {
    *rY  = (Sqrt(195/(2.*Pi))*Cos(4*phi)*Power(Cos(theta/2.),2)*(35 + 44*Cos(theta) + 33*Cos(2*theta))*Power(Sin(theta/2.),6))/8.;
    *iY  = -(Sqrt(195/(2.*Pi))*Power(Cos(theta/2.),2)*(35 + 44*Cos(theta) + 33*Cos(2*theta))*Sin(4*phi)*Power(Sin(theta/2.),6))/8.;
    return;
  } else if (m==-3) {
    *rY  = (3*Sqrt(13/Pi)*Cos(3*phi)*(381*Cos(theta/2.) + 295*Cos((3*theta)/2.) + 165*Cos((5*theta)/2.) + 55*Cos((7*theta)/2.))*Power(Sin(theta/2.),5))/64.;
    *iY  = (-3*Sqrt(13/Pi)*(381*Cos(theta/2.) + 295*Cos((3*theta)/2.) + 165*Cos((5*theta)/2.) + 55*Cos((7*theta)/2.))*Sin(3*phi)*Power(Sin(theta/2.),5))/64.;
    return;
  } else if (m==-2) {
    *rY  = (Sqrt(13/Pi)*Cos(2*phi)*(1709 + 3096*Cos(theta) + 2340*Cos(2*theta) + 1320*Cos(3*theta) + 495*Cos(4*theta))*Power(Sin(theta/2.),4))/256.;
    *iY  = -(Sqrt(13/Pi)*(1709 + 3096*Cos(theta) + 2340*Cos(2*theta) + 1320*Cos(3*theta) + 495*Cos(4*theta))*Sin(2*phi)*Power(Sin(theta/2.),4))/256.;
    return;
  } else if (m==-1) {
    *rY  = (Sqrt(65/(2.*Pi))*Cos(phi)*(574*Cos(theta/2.) + 504*Cos((3*theta)/2.) + 384*Cos((5*theta)/2.) + 231*Cos((7*theta)/2.) + 99*Cos((9*theta)/2.))*Power(Sin(theta/2.),3))/128.;
    *iY  = -(Sqrt(65/(2.*Pi))*(574*Cos(theta/2.) + 504*Cos((3*theta)/2.) + 384*Cos((5*theta)/2.) + 231*Cos((7*theta)/2.) + 99*Cos((9*theta)/2.))*Sin(phi)*Power(Sin(theta/2.),3))/128.;
    return;
  } else if (m==0) {
    *rY  = (Sqrt(1365/Pi)*(35 + 60*Cos(2*theta) + 33*Cos(4*theta))*Power(Sin(theta),2))/512.;
    *iY  = 0;
    return;
  } else if (m==1) {
    *rY  = (Sqrt(65/(2.*Pi))*Cos(phi)*Power(Cos(theta/2.),3)*(574*Sin(theta/2.) - 504*Sin((3*theta)/2.) + 384*Sin((5*theta)/2.) - 231*Sin((7*theta)/2.) + 99*Sin((9*theta)/2.)))/128.;
    *iY  = (Sqrt(65/(2.*Pi))*Power(Cos(theta/2.),3)*Sin(phi)*(574*Sin(theta/2.) - 504*Sin((3*theta)/2.) + 384*Sin((5*theta)/2.) - 231*Sin((7*theta)/2.) + 99*Sin((9*theta)/2.)))/128.;
    return;
  } else if (m==2) {
    *rY  = (Sqrt(13/Pi)*Cos(2*phi)*Power(Cos(theta/2.),4)*(1709 - 3096*Cos(theta) + 2340*Cos(2*theta) - 1320*Cos(3*theta) + 495*Cos(4*theta)))/256.;
    *iY  = (Sqrt(13/Pi)*Power(Cos(theta/2.),4)*(1709 - 3096*Cos(theta) + 2340*Cos(2*theta) - 1320*Cos(3*theta) + 495*Cos(4*theta))*Sin(2*phi))/256.;
    return;
  } else if (m==3) {
    *rY  = (3*Sqrt(13/Pi)*Cos(3*phi)*Power(Cos(theta/2.),5)*(381*Sin(theta/2.) - 295*Sin((3*theta)/2.) + 165*Sin((5*theta)/2.) - 55*Sin((7*theta)/2.)))/64.;
    *iY  = (3*Sqrt(13/Pi)*Power(Cos(theta/2.),5)*Sin(3*phi)*(381*Sin(theta/2.) - 295*Sin((3*theta)/2.) + 165*Sin((5*theta)/2.) - 55*Sin((7*theta)/2.)))/64.;
    return;
  } else if (m==4) {
    *rY  = (Sqrt(195/(2.*Pi))*Cos(4*phi)*Power(Cos(theta/2.),6)*(35 - 44*Cos(theta) + 33*Cos(2*theta))*Power(Sin(theta/2.),2))/8.;
    *iY  = (Sqrt(195/(2.*Pi))*Power(Cos(theta/2.),6)*(35 - 44*Cos(theta) + 33*Cos(2*theta))*Sin(4*phi)*Power(Sin(theta/2.),2))/8.;
    return;
  } else if (m==5) {
    *rY  = (Sqrt(2145/Pi)*Cos(5*phi)*Power(Cos(theta/2.),7)*(1 - 3*Cos(theta))*Power(Sin(theta/2.),3))/2.;
    *iY  = (Sqrt(2145/Pi)*Power(Cos(theta/2.),7)*(1 - 3*Cos(theta))*Sin(5*phi)*Power(Sin(theta/2.),3))/2.;
    return;
  } else if (m==6) {
    *rY  = (3*Sqrt(715/Pi)*Cos(6*phi)*Power(Cos(theta/2.),8)*Power(Sin(theta/2.),4))/2.;
    *iY  = (3*Sqrt(715/Pi)*Power(Csc(theta/2.),4)*Sin(6*phi)*Power(Sin(theta),8))/512.;
    return;
  }

} else errorexit("increase lmax in mathematica script");


}
/* SphericalHarmonics_precomp.c */

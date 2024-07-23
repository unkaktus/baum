/* utility.c */
/* Bernd Bruegmann 6/02, Wolfgang Tichy 7/2005 */

#include "bam.h"
#include "adm.h"



/* compute determinant */
double detg(double g11, double g12, double g13, 
	    double g22, double g23, double g33)
{
  double gginv11, gginv12, gginv13;

  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  return g11*gginv11 + g12*gginv12 + g13*gginv13;
}




/* compute inverse */
double invg(double g11, double g12, double g13, 
	    double g22, double g23, double g33,
	    double *i11, double *i12, double *i13, 
	    double *i22, double *i23, double *i33)
{
  double detg, gginv11, gginv12, gginv13, gginv22, gginv23, gginv33;

  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  gginv22 = g11*g33 - g13*g13;
  gginv23 = g12*g13 - g11*g23;
  gginv33 = g11*g22 - g12*g12;
  detg = g11*gginv11 + g12*gginv12 + g13*gginv13;
  *i11 = gginv11/detg;
  *i12 = gginv12/detg;
  *i13 = gginv13/detg;
  *i22 = gginv22/detg;
  *i23 = gginv23/detg;
  *i33 = gginv33/detg;
  return detg;
}











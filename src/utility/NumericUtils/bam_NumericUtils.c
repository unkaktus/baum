/* bam_NumericUtils.c */
/* Wolfgang Tichy 2003 */

#include "bam.h"
#include "NumericUtils.h"


void bam_NumericUtils(void) 
{
  printf("Adding NumericUtils\n");

  AddPar("sphere_x0", "0", "origin of coordinate sphere");
  AddPar("sphere_y0", "0", "origin of coordinate sphere");
  AddPar("sphere_z0", "0", "origin of coordinate sphere");
  
  AddPar("integrate_over_sphere", "Simpson", "");
}

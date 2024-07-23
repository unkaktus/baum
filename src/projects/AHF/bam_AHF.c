/* bam_AHF.c */
/* Jose Gonzalez 03/07 */

#include "bam.h"
#include "AHF.h"


void bam_AHF() 
{
  if (!Getv("physics", "AHF")) return;
  printf("Adding AHF\n");

  AddFun(ANALYZE, AHF, "look for apparent horizons");

  AddPar("ahf_verbose"  ,"no" ,"Print information");
  AddPar("ahf_ntheta"   ,"50" ,"Number of circles in theta direction");
  AddPar("ahf_nphi"     ,"50" ,"Number of points in phi direction");
  AddPar("ahf_nhorizons","1"  ,"Number of horizon to search (1 or 2)");
  AddPar("ahf_common"   ,"no" ,"Search for common horizon?");
  AddPar("ahf_common_time","-1.0" ,"Search for common horizon starting at this time");
  AddPar("ahf_flow_iter","2000","Number of flow iterations");
  AddPar("ahf_mass_tol","0.001","mass tolerance");
  AddPar("ahf_hmean_tol","100.0","h mean tolerance");
  AddPar("ahf_time","1.0","how often we search for horizons?");
  AddPar("ahf_interpolation_order","4","order of interpolation onto surface");

  AddPar("ahf_LMAX","10","maximum number of Legendre polynomials");
  AddPar("ahf_output","yes","");
  AddPar("ahf_output_xyt","no","");



  AddVar("ahf_dgd", "ijk+ikj", "first derivative of the metric");
  
  
  AddPar("ahf0_x", "0.0","");
  AddPar("ahf0_y", "0.0","");
  AddPar("ahf0_z", "0.0","");
  AddPar("ahf0_m", "0.0","");
  AddPar("ahf0_r", "-1.0","");
  AddPar("ahf1_x", "0.0","");
  AddPar("ahf1_y", "0.0","");
  AddPar("ahf1_z", "0.0","");
  AddPar("ahf1_m", "0.0","");
  AddPar("ahf1_r", "-1.0","");
  AddPar("ahf2_x", "0.0","");
  AddPar("ahf2_y", "0.0","");
  AddPar("ahf2_z", "0.0","");
  AddPar("ahf2_m", "0.0","");
  AddPar("ahf2_r", "-1.0","");
}

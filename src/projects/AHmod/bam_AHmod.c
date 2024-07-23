/* bam_AHmod.c */
/* Jose Gonzalez 03/07 */
/* Norbert Lages, 02/09 */

#include "bam.h"
#include "AHmod.h"


void bam_AHmod() 
{
  if (!Getv("physics", "AHmod")) return;
  printf("Adding AHmod\n");
  if ((Geti("amr_lmax2") > 0 )&&(Geti("amr_lmax2") != Geti("amr_lmax"))) {
   //necessary if amr_lmax2 is used 
  AddFun(POST_ANALYZE, AHmod, "look for marginally trapped surfaces");
  }
  else {
  AddFun(ANALYZE, AHmod, "look for marginally trapped surfaces");
  }
  AddPar("AHmod_verbose"   ,"no"   ,"Print information");

  AddPar("AHmod_always_search"   ,"no"   ,"Always search for MOTS");
  AddPar("AHmod_start_level", "-1"   ,"Level on which start MOTS search");

  AddPar("AHmod_ntheta"    ,"30"   ,"Number of circles in theta direction");
  AddPar("AHmod_nphi"      ,"60"   ,"Number of points in phi direction");
  AddPar("AHmod_flow_iter" ,"1000" ,"Maximal number of flow iterations");
  AddPar("AHmod_mass_tol"  ,"0.001","mass tolerance, stopping criterion for flow");
  AddPar("AHmod_hmean_tol" ,"100.0","h mean tolerance");
  AddPar("AHmod_time"      ,"1.0"  ,"how often we search for MOTS?");
  AddPar("AHmod_LMAX"      ,"10"   ,"maximum number of Legendre polynomials");
  AddPar("AHmod_uselast"   ,"yes"  ,"use the last surface found as initial guess");
  AddPar("AHmod_initial_guess_expand","1.0","factor to make initial guess larger");
  AddPar("AHmod_UseOptimalLevel",     "no", "Pick only one level, that seems to be optimal, and search within");
  AddPar("AHmod_box_savety_factor",   "1.2","increased initial guess sphere has to fit in box");
  AddPar("AHmod_initial_radius",      "1.0","this is the guess when no puncture data is available");
  
  AddPar("AHmod_nhorizons","0","Such many surfaces (MOTS) will be searched");
  AddPar("AHmod_searchMTS","0","all info about when, which and where to search for MOTS");
  
  AddPar("AHmod_output",     "yes","output information about MOTS");
  AddPar("AHmod_output_xyt", "no", "output location points of MOTS");
  AddPar("AHmod_output_lm",  "no", "output spherical modes of MOTS");
  
  AddVar("AHmod_dgd", "ijk+ikj", "first derivative of the metric");
  AddPar("AHmod_interpolation_order","4","order of interpolation onto surface");
  
  AddPar("AHmod_merger_distance", "0.1", "distance in M at which BHs are considered as merged");

  // more parameters are added automatically within AHmod
  // especially some for each surface 
}

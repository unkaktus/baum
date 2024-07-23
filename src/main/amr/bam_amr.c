/* bam_amr.c */
/* Bernd Bruegmann 12/99 */

#include "bam.h"
#include "amr.h"



void bam_amr() 
{
  printf("Adding amr\n");

  /* variables */
  AddConstantVar("x", "", "x coordinate");
  AddConstantVar("y", "", "y coordinate");
  AddConstantVar("z", "", "z coordinate");
  AddVar("flagregrid",   "", "flag for regridding");
  AddVar("flagrestrict", "", "flag for restriction");
  AddVar("flagprolong",  "", "flag for prolongation");
  
  
  /* parameters */
  AddPar("amr", "onedt", "amr type [onedt,bo]");
  AddPar("amr_bo_dxmax", "0",
	 "if bo and dxmax > 0, do onedt for dx > dxmax/2");
  AddPar("amr_bo_lmin","0","can be used instead of amr_bo_dxmax");

  AddPar("amr_npunctures",        "1",  "amr: number of amr punctures, where we refine");
  AddPar("amr_lmax",              "0",  "amr: maximal number of refinements");
  AddPar("amr_lmax2",       "-1",  
         "amr: maximal number of refinements around puncture0 < amr_lmax");

  AddPar("amr_lfinal",    "0",  "number of levels after amr_ladd_time");
  AddPar("amr_ladd_time", "0",  "time to switch to amr_lfinal levels");
  
  if (Getd("amr_ladd_time") > 0) {
    AddFun(POST_OUTPUT, add_level_while_evolving, "add additional level after amr_ladd_time");
    if(Geti("amr_lfinal")<=Geti("amr_lmax")) errorexit("watch your lmax and lfinal!");
  }

  AddPar("amr_nxyz", "",
	 "list of box sizes, overrides other settings, last one is repeated");

  
  AddPar("grid_half", "", "internal parameter about grid extent");
  AddPar("grid_fournplustwo", "yes",
	 "whether to enforce grid size 4*n+2 automatically");

  /* advection stuff */
  AddPar("advection_lopsided", "1",   // one away from centered
         "for 4-th order stencils: one-sided, lop-sided, centered [0,1,2]");
  AddPar("advection_lopsided6", "3",  // legacy: default is centered
         "for 6-th order stencils: one-sided, 2lop-sided, lop-sided, centered [0,1,2,3]");
  AddPar("advection_lopsided8", "3",  // one away from centered
         "for 8-th order stencils: one-sided, 2lop-sided, lop-sided, centered ...");
  AddPar("advection_lopsided10", "4", // one away from centered
         "for 10-th order stencils: one-sided, 2lop-sided, lop-sided, centered ...");
  /* for arbitrary order, count shifts from centered, which is opposite from the above */
  AddPar("advection_lopsidedN", "1",  // one away from centered
         "for N-th order stencils: shifted by 1, 2, 3, ...");
  
  AddPar("amr_nbuffer", "4", "number of buffer cells");
  AddPar("amr_buffer_addextra", "yes", 
         "add cells rather than using existing cells [no,yes]");
    
  AddPar("amr_move_lcube", "0",
         "only move boxes with l > lcube, use fmr cubes for l <= lcube");

  AddPar("amr_move_nxyz", "0", "size of moving cubes");
  
  
  
  /* special parameter/variables for spherical shells */
  if (Getv("grid","shells")) {
    AddConstantVar("shells_r",     "",   "r coordinate");
    AddConstantVar("shells_R",     "",   "R(r) coordinates ... for fisheye");
    AddConstantVar("shells_phi",   "",   "phi coordinate");
    AddConstantVar("shells_theta", "",   "theta coordinate");
    AddPar("amr_shells_nr",        "-1", "points in the r direction for spherical shells");
    AddPar("amr_shells_nphi",      "-1", "points in the phi/theta direction fspherical shells");
    AddPar("amr_shells_stagger",   "yes","shells are staggered like normal boxes?");
    AddPar("amr_shells_RP",        "yes","do RP?");
    AddPar("amr_shells_optimize",  "yes","mpi optimized version for Nx1x1 parallelization?");
    AddPar("amr_shells_symmetry",  "yes","yes for using the same symmetry as boxes/ no to use full shells");
    AddPar("dissipation_factor_shells", "-1","");
    AddPar("amr_shells_stretch",   "1.",  "par of fisheye coordinates");
    AddPar("amr_shells_r0",        "0.",  "par of fisheye coordinates");
    AddPar("amr_shells_eps",       "1.",  "par of fisheye coordinates");
  }
  AddPar("dissipation_factor_move", "-1","");
  
  /* special development flag */
  AddPar("amr_move_force", "-1", 
	 "should be -1; set to 0, 1, ... to enable forced motion");
  AddPar("amr_move_force_o", "0.7","");
  AddPar("amr_move_oldstable", "no", 
	 "revert to older version, stable but slow");
  if (GetsLax("amr_move_stable"))
    errorexit("Don't use parameter amr_move_stable, it has been removed.\n"
	      "Use amr_move_oldstable, but only if needed. Default is \"no\" "
              "since new version is now stable and faster.");
  
  
  AddPar("storage_verbose", "no", 
         "verbose mode for memory allocation [no,yes,mem]");
  AddPar("amr_verbose", "no", 
         "verbose mode for amr [no,all,move]");
  
  
  
  /* special paramter for different interpolation schemes */
  AddPar("matter_interpolate_vars",   "",     "let this empty!!! this is filled by grhd ");
  AddPar("matter_interpolate_scheme_restriction", "WENO", "lagrange,WENO,WENOZ,linear");
  AddPar("matter_interpolate_scheme_prolongation", "WENO", "lagrange,WENO,WENOZ,linear");
  AddPar("matter_interpolate_order",  "4","   ");
  
  
  
  AddPar("ExitIfNAN",       "no", "exit if NAN is found ");
  AddPar("ExitIfNAN_vars",  "",   "vars to check");
  if (Getv("ExitIfNAN", "yes")) 
    AddFun(POST_EVOLVE, ExitIfNAN, "exit if NAN");
}




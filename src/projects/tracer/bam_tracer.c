/* bam_tracer.c */
/* A Neuweiler 12/22 */

#include "bam.h"
#include "tracer.h"


void bam_tracer() 
{
  if (!Getv("physics", "tracer")) return;
  printf("Adding tracer\n"); 

  AddFun(PRE_EVOLVE, tracer, "set evolution for tracer particles");
  AddFun(OUTPUT, output_tracer, "write output for tracer particles");
  AddFun(POST_MOVEBOX, boxes_tracer, "update box boundaries for communication");

  AddPar("tracer_number", "100000", "number of tracer particles");
  
  AddPar("tracer",    "no", "use tracer [yes/no]");
  AddPar("tracer_init_time",        "0", "time to initialize tracer particles");
  AddPar("tracer_init_distance",   "-1", "proper distance of compact objects to initialize tracer particles");
  AddPar("tracer_init_level",       "0", "level to put tracers");

  AddPar("tracer_init_distribution",   "uniform", "how to distribute particles [uniform/weighted]");
  AddPar("tracer_init_density_threshold",   "0", "threshold density to initialize particles");
  
  AddPar("tracer_interpolation_order", "2", "order of interpolation to Tracer Particels");

  AddPar("tracer_outtime", "-1", "when to output based on time");
  AddPar("tracer_output", "grhd_rho", "variables to output for each Tracer");

}

/* bam_evolve.c */
/* Bernd Bruegmann 6/02 */

#include "bam.h"
#include "evolve.h"




void bam_evolve() 
{
  printf("Adding evolve\n");

  /* functions */
  AddFun(EVOLVE, evolve, "evolve with method of lines like icn, rk, dlf");

  /* variables */

  /* parameters */
  AddPar("evolution_method", "icn", "evolution method [icn,rk,euler,ngs2]");
  AddPar("evolution_method_rk", "rk4", 
	 "Runge-Kutta method [euler,midpoint,rk3,rk3a,rk3b,rk4,rk4a,rk5]");
  AddPar("evolve_ngs_iterations", "5", "number of iterations");

  AddPar("evolve_persist", "yes", "whether additional memory persists");
 // evolve_no_memory = vlalloc(NULL);

  AddPar("evolve_euler_debug", "no",  
	 "obtain rhs in variable after one timestep");
  AddPar("evolve_store_rhs", "no",  "whether to compute rhs for analysis");
  if (Getv("evolve_store_rhs", "yes"))
    AddFun(EVOLVE, evolve_store_rhs, "compute and store rhs");

  AddPar("evolve_compute_change", "", 
	 "list of variables for which change is to be computed");

 
}





void bam_evolve_final() 
{
  printf("Removing evolve\n");
  
  //vlfree(evolve_no_memory);
}



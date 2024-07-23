/* bam_bampi.c */
/* Bernd Bruegmann 5/02 */

#include "bam.h"
#include "bampi.h"



void bam_bampi(void) 
{
  printf("Adding bampi\n");

  if (Geti("amr_lmax") > 0 && bampi_size() > 1) {
    if (0) {
      printf("You are running on more than 1 processor and\n");
      printf("have requested amr_lmax = %d.\n\n", Geti("amr_lmax"));
    }
    if (!Getv("amr", "newfmr"))
      errorexit("Cannot run in parallel, "
		"use \"amr = newfmr\" for parallelized amr");
  // printf("WARNING: this is highly experimental code which will crash!\n\n");
  // what else is new ...
  }

  /* parameters */
  AddPar("bampi_nghosts", "2", "number of ghosts in each direction");

  AddPar("bampi_xsize", "0", "nprocs for x direction, ignored if 0");
  AddPar("bampi_ysize", "0", "nprocs for y direction, ignored if 0");
  AddPar("bampi_zsize", "0", "nprocs for z direction, ignored if 0");

  AddPar("bampi_sizes_constant", "yes",
	 "don't touch: whether to keep xyz-sizes constant "
	 "or to vary them [yes,no]");

  AddPar("bampi_symmpointsextra", "yes",
	 "don't touch: whether to add symmetry ghosts after [no] "
	 "or before [yes] splitting the box");

  AddPar("bampi_lowlatency", "no", 
	 "send many small rather than few large messages");

  AddPar("bampi_verbose", "no", "whether to talk about it");

  AddPar("bampi_timer_on", "no", "turn on timer for various functions");
  AddPar("bampi_timer_barrier", "no", "invoke mpi barrier before timing");
  AddPar("bampi_timer_reset_after_every_iteration", "no", "reset all timers after every iteration");

}

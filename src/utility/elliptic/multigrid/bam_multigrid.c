/* bam_multigrid.c */
/* Bernd Bruegmann 11/02 */

#include "bam.h"
#include "multigrid.h"



void (*multigrid_setbound)(tL *, tVarList *);
	
void bam_multigrid() 
{
  if (!Getv("physics", "multigrid")) return;
  printf("Adding multigrid\n");

  /* functions */
  //AddFun(FLAGREF,     flagsimplegrid, "flags for box/sphere refinement");
  //AddFun(INITIALDATA, wavedata, "compute scalar wave initial data");
  //AddFun(EVOLVE,      waveevolve,   "wave evolution");

  /* variables */
  AddVar("mg_redblack", "", "flag for red-black ordering");
  if (IndLax("flagrestrict") == -1)
    AddVar("flagrestrict", "", "flag for restriction");
  if (IndLax("flagprolong") == -1)
    AddVar("flagprolong",  "", "flag for prolongation");


  /* parameters */
  AddPar("mg_relaxation", "rgbs",   "relaxation scheme");
  AddPar("mg_persist", "no", "whether to keep multigrid stack");
  AddPar("mg_verbose", "no", "whether to be verbose [no,some,yes]");
  AddPar("mg_typeofcycles", "fullV", "type of multigrid cycling");

  AddPar("mg_levels_nmax", "-1", "number of multigrid levels, -1 to max");

  AddPar("mg_size_nmin", "5", "minimal global size of coarsest grid");
  AddPar("mg_size_nmin_local", "4", 
	 "minimal local size [4] (2 works in some cases, but not always)");

  AddPar("mg_relax_npre",  "3", "number of relaxation sweeps before");
  AddPar("mg_relax_npost", "3", "number of relaxation sweeps after");

  AddPar("mg_exact_nmin", "150", "minimum number of exact solve iterations");
  AddPar("mg_exact_nmax", "500", "maximum number of exact solve iterations");

  AddPar("mg_vcycle_nmin", "0",   "minimum number of V cycles");
  AddPar("mg_vcycle_nmax", "100", "maximum number of V cycles");

  AddPar("mg_wcycle_nmin", "0", "minimum number of W cycles");
  AddPar("mg_wcycle_nmax", "0", "maximum number of W cycles");
  AddPar("mg_wcycle_factol", "1.0", "Oversolve factor for W cycles");
  AddPar("mg_setbound", "set_boundary_elliptic",
  	 "function used to set boundaries in multigrid");
}

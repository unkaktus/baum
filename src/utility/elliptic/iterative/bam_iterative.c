/* bam_iterative.c */
/* Bernd Bruegmann 01/00 */

#include "bam.h"
#include "iterative.h"



void bam_iterative() 
{
  printf("Adding elliptic/iterative\n");

  /* global parameter */
  AddPar("ell_verbose", "yes", "talk about it");

  /* check whether there is more to do */
  if (!Getv("physics", "poisson") && 
      !Getv("physics", "Poisson") &&
      !Getv("physics", "punctures")) return;

  /* functions */

  /* variables */
  AddVar("u", "", "scalar field for Poisson equation");
  AddVar("f", "", "right hand side for Poisson equation");
  AddVar("r", "", "residuum for Poisson equation");

  /* tests with poisson equation */
  if (Getv("physics", "poisson") || Getv("physics", "Poisson")) {

    /* functions */
    AddFun(INITIALDATA_SET, poisson, "solve Poisson equation");
    
    /* variables */
    AddVar("au", "", "analytic solution for Poisson equation");
    AddVar("ar", "", "residuum of analytic solution for Poisson equation");
    AddVar("at", "",
	   "truncation error of analytic solution for Poisson equation");

    /* parameters */
    AddPar("poisson_type", "cos*cos*cos", 
	   "type of data [cos*cos*cos,charges]");
    AddPar("poisson_tx", "1", "period in x");
    AddPar("poisson_ty", "1", "period in y");
    AddPar("poisson_tz", "1", "period in z");
    AddPar("poisson_charge0", "0", "charge at outer boundary [#,1/r]");
    AddPar("poisson_charge1", "1", "charge at first inner boundary");
    AddPar("poisson_charge2", "-1", "charge at second inner boundary");

    AddPar("poisson_tol", "1e-10", "l2 tolerance");
    AddPar("poisson_itmax", "100", "maximal number of iterations");
    AddPar("poisson_precon", "none", "preconditioning [none,diag]");
    AddPar("poisson_solver", "bicgstab", 
	   "elliptic solver [bicgstab,cgs,jacobi,sor]");
  }

  /* temporary variables, should be inside if statement to avoid
     cluttering the variable name space 
  */
  AddVar("iterative_p",  "", "temporary variable");
  AddVar("iterative_ph", "", "temporary variable");
  AddVar("iterative_rt", "", "temporary variable");
  AddVar("iterative_s",  "", "temporary variable");
  AddVar("iterative_sh", "", "temporary variable");
  AddVar("iterative_t",  "", "temporary variable");
  AddVar("iterative_v",  "", "temporary variable");
}

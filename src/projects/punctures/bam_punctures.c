/* bam_punctures.c */
/* Bernd Bruegmann 6/02, Wolfgang Tichy 2/2004 */

#include "bam.h"
#include "punctures.h"

int (*punctures_solver)
     (tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c,
      int imax, double tol, double *res,
      void (*lop)(tL *, tVarList *, tVarList *),
      void (*precon)(tL *, tVarList *, tVarList *));

int (*linear_solver)
    (tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c,
    int imax, double tol, double *res,
    void (*lop)(tL *, tVarList *, tVarList *),
    void (*precon)(tL *, tVarList *, tVarList *));


void bam_punctures() 
{
  /* global in some cases: black hole parameters */
  if (!Getv("physics", "punctures")) return;
  printf("Adding punctures\n");
  
  

  if(GetvLax("amr_fmr", "nestedboxes") || GetvLax("amr_fmr", "gradedboxes"))
  {
    AddPar("bh_fmr", "default", 
            "how to do fmr around BHs [default,cover_origin,force_origin]");
    /* search for cover_origin explanation in punctures/FlagRegrid.c */
    AddPar("bh_fmr_ry_o_r_max", "1.0", "max ratio of ry to r");
    AddPar("bh_fmr_ry_o_rx_max", "2.0", "max ratio of ry to rx");
  }

  AddPar("punctures_from_fit", "no",
          "whether to reset puncture pars from sequence fit [no,yes]");
  if(Getv("punctures_from_fit", "yes"))
  {
    AddPar("punctures_from_fit_x", "2.0", 
            "value of separation parameter x=D/(2M) for sequence fit");
    AddFun(PRE_INITIALDATA, ComputePunctureParametersFromFit,
            "set parameters based on fit"); 
  }

  AddPar("punctures_from_PN", "no",
          "whether to reset puncture pars from PN formula [no,yes]");
  if(Getv("punctures_from_PN", "yes"))
  {
    AddPar("punctures_from_PN_D", "7.0", 
            "value of separation parameter D for PN formula");
    AddPar("punctures_from_PN_m1", "0.5", "value for m1");
    AddPar("punctures_from_PN_m2", "0.5", "value for m2");
    AddFun(PRE_INITIALDATA, ComputePunctureParametersFromPN,
            "set parameters based on fit"); 
  }


  /* punctures with a given ADM mass */
  AddPar("punctures_with_given_M_ADM", "no", "whether to reset puncture bare"
          " masses to obtain a specfific ADM mass at each puncture [no,yes]");
  if(Getv("punctures_with_given_M_ADM", "yes"))
  {
    AddPar("puncture_M_ADM_1", "0.5", "ADM mass at puncture 1");
    AddPar("puncture_M_ADM_2", "0.5", "ADM mass at puncture 2");
    AddPar("puncture_M_ADM_tol", "0.0001",
            "tolerance for ADM mass at punctures");
    AddPar("puncture_M_ADM_newton_lnsrch", "no",
            "whether we use newton_lnsrch [yes,no]");
    AddPar("puncture_M_ADM_MAXITS", "200", "max. newton iterations");
  }
  

  /* functions */
  if(Getv("punctures_with_given_M_ADM", "yes"))
    AddFun(INITIALDATA_SET, PunctureData_with_given_M_ADM,
           "compute puncture data with specific ADM mass");

  else /* default */
    AddFun(INITIALDATA_SET, PunctureData,   "compute puncture data");


  /* variables */
  AddVar("punctures_u", "", 
	 "regular correction to Brill-Lindquist conformal factor");
  AddVar("punctures_v", "", "linear correction to u");
  AddVar("punctures_a", "", "coefficient a in Hamiltonian constraint");
  AddVar("punctures_b", "", "coefficient b in Hamiltonian constraint");
  AddVar("punctures_c", "", "coefficient c for maximal slicing");
  AddVar("punctures_cf", "","coefficient cf for maximal slicing");
  AddVar("punctures_d", "", "coefficient d for shift");
  AddVar("punctures_e", "I","coefficient e^i for shift");
  AddVar("punctures_f", "", "right-hand-side in Hamiltonian constraint");
  AddVar("punctures_r", "", "residuum of scalar equation");
  AddVar("punctures_s", "I","residuum of vector equation");
  AddVar("punctures_lu","", "Laplace of function for mass evaluation");

  /* parameters */
  AddPar("punctures_solver", "bicgstab", 
         "solver [bicgstab,jacobi,multgrid,Newton,none]");
  AddPar("punctures_persist", "no", "whether to keep storage on");
  AddPar("punctures_itmax", "10", "maximal number of iterations");
  AddPar("punctures_tolerance", "1e-6", "tolerance in l2-norm of residuum");
  AddPar("punctures_linSolver", "bicgstab",
         "linear solver used [bicgstab, HYPRE]");
  AddPar("punctures_linSolver_itmax", "20", "max num of linSolver iterations");
  AddPar("punctures_linSolver_tolFac","0.1",
         "tol for linSolver is tol * linSolver_tolFac");
  AddPar("punctures_use_HYPRE_tol","yes",
         "whether tolerance from HYPRE_options is used");

  AddPar("punctures_fourth_order","no","Use compact 4th-order stencil");

  AddPar("punctures_Absorbpsi", "no",
         "whether (in the end) psi is absorbed into the 3-metric [no,yes]");

  AddPar("punctures_mass", "no",
	 "compute ADM/Komar mass [yes,no,ZeroCorrection]");

  if (Getv("punctures_solver", "spectral")) {
    AddPar("punctures_ps_order", "8",
	   "order of interpolation from spectral data"); 
    AddPar("punctures_ps_file", "", 
	   "if set, read spectral data from this file");
     AddPar("punctures_ps_interpol_check", "none", 
           "if yes: NR-method for inverse transformation");
    if (strlen(Gets("punctures_ps_file"))) {
      AddFun(PRE_GRID, read_ps_parameters,
	     "read spectral puncture parameters from file");}
  }

  AddPar("punctures_lapse", "none",
	 "lapse [none,maximal,psiBL^(-2),rtoN_atPunc,psipower]");
  AddPar("punctures_alp_psipower", "0", "Pow of psi for initial lapse");
  AddPar("punctures_lapse_at_puncture", "-1", "[(number),MK_equals_MADM]");
  AddPar("punctures_lapse_MK_MADM_diff", "0.001",
	 "allowed deviation for MK_equals_MADM");
  AddPar("punctures_lapse_rPower_atPunc", "4",
         "exponent N for rtoN_atPunc lapse");

  AddPar("punctures_shift", "none", "shift [none,thinsandwich]");
  AddPar("punctures_shift_add_beta_rot", "yes",
	 "add a rotational piece Omega_b x R to shift");
  AddPar("punctures_shift_falloff", "2", "shift falloff [2,3]");
  AddPar("punctures_shift_rhs", "delA",
	 "form of right-hand side of shift equation [delA,Adel]");

  AddPar("punctures_shift_test", "no", "compute diagnostics for shift");
  if (Getv("punctures_shift_test", "yes")) {
    AddVar("punctures_dgdt",  "ij+ji", "time derivative of conformal metric");
    AddVar("punctures_lbeta", "ij+ji", "l-operator of shift");
    AddVar("punctures_lbrhs", "ij+ji", "rhs of lbeta equation");
  } 

  AddPar("punctures_MK_eq_MADM_Sequence", "no",
   	 "create sequence of M_K = Lapse * M_ADM [yes,no]");   
  AddPar("punctures_Sequence_D1", "2", 
  	 "smallest coord sep. for MK_eq_MADM_Sequence");
  AddPar("punctures_Sequence_Dstep", "1.0",
         "coord step size along MK_eq_MADM_Sequence");
  AddPar("punctures_Sequence_PointNumber", "1",
	 "number of points on sequence");
  AddPar("punctures_Sequence_D2", "use_Dstep", 
  	 "largest coord sep. [(some number), use_Dstep]");
  AddPar("punctures_Sequence_newtMAXITS", "200", "max. newton iterations");
  AddPar("punctures_Sequence_newtTOLF",  "1e-4", "newton tolerence");
  AddPar("punctures_solv_included",  "no", 
         "If puncture solver within another IDproject");
  /* Robin boundary */
  AddPar("punctures_boundary", "robin", "boundary condition");
  AddPar("punctures_robin_falloff", "1", "If zero, use dirchlet boundary condition");

  if (Getv("punctures_boundary", "robin")) {
    AddVar("robinindex", "", "index for Robin boundary");
    AddVar("robinflag",  "", "temporary flag for Robin boundary");
  }

}

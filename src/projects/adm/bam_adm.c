/* bam_adm.c */
/* Bernd Bruegmann 6/02 */

#include "bam.h"
#include "adm.h"




void bam_adm() 
{
  if (!Getv("physics", "adm")) return;
  printf("Adding adm\n");  
  
  /* functions */
  AddFun(ANALYZE, computeadmconstraints, "compute ADM constraints");
  
  
  
  /* variables */
  AddVar("adm_g",    "ij+ji", "metric: g_ij");
  AddVar("adm_K",    "ij+ji", "extrinsic curvature: K_ij");
  AddVar("alpha",    "",      "lapse: alpha");
  AddVar("beta",     "I",     "shift: beta^i");
  AddVar("betadot",  "I",     "time derivative of shift: del_t beta^i = B^i");
  
  /* additional matter variables */
  if (Getv("physics","matter")) {
    AddVar("adm_rho",  "",      "ADM rho");
    AddVar("adm_S",    "I",     "ADM S^i");
    AddVar("adm_SS",   "ij+ji", "ADM S_ij");
    AddVar("adm_ST",   "",      "ADM S_i^i");
    AddVar("adm_detg", "",      "determinant of the full 3-metric");
  }

  /* constraints */
  AddVar("ham",      "",      "Hamiltonian constraint");
  AddVar("mom",      "i",     "momentum constraint");
  
  AddPar("adm_normalizedConstraints", "no",
         "whether we compute normalized constraints [no,yes TermByTerm]");
  if (Getv("adm_normalizedConstraints", "yes")) {
    AddVar("normham",  "",      "normalized Hamiltonian constraint");
    AddVar("normmom",  "i",     "normalized momentum constraint");
  }

  /* INTITALDATA conformal factor ... does not exist in evolution */
  AddConstantVar("adm_psi",       "",      "conformal factor");
  AddConstantVar("adm_dpsiopsi",  "i",     "(del_i psi) / psi");
  AddConstantVar("adm_ddpsiopsi", "ij+ji", "(del_i del_j psi) / psi");
  AddFun(PRE_INITIALDATA,    adm_init_conftrans, "adm_init_conftrans");
  AddFun(INITIALDATA_FINISH, adm_undo_conftrans, "adm_undo_conftrans");

  /* 1+log slicing may require the initial value of the trace of K_ij */
  AddVar("adm_K0", "",      "initial value of trace of K_ij");
  AddFun(POST_INITIALDATA, set_K_initial, 
	 "compute and save initial value of the trace of K_ij");
  
  /* check for NANs */
  AddPar("adm_NANcheck", "no", "check if ham contains NAN or INF");
  
  /* simple analytic initial data */
  AddPar("adm_data", "none", "set simple analytic ADM data");
  if (Getv("adm_data", "Minkowski"))
    AddFun(INITIALDATA_SET, adm_Minkowski, "set Minkowski data");
  if (Getv("adm_data", "puncture"))
    AddFun(INITIALDATA_SET, adm_puncture, "set Schwarzschild puncture data");
  if (Getv("adm_data", "gauge_wave")) {
    AddFun(INITIALDATA_SET, adm_gauge_wave, "set Minkowski + gauge wave");
    AddPar("gauge_wave_diameter","5.","");
    AddPar("gauge_wave_amplitude","0.001","");
  }
  
  
}




void bam_adm_final() 
{
  if (!Getv("physics", "adm")) return;
  printf("Removing adm\n");
  
}


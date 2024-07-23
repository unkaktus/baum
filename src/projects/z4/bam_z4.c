/* bam_z4.c */
/* Wolfgang Tichy  4/2004 */

#include "bam.h"
#include "z4.h"




void bam_z4(void) 
{
  if (!Getv("physics", "z4")) return;
  printf("Adding z4\n");

  /* functions */
  AddFun(POST_INITIALDATA, z4_startup, 
	 "initialize z4 system from adm initial data");
  AddFun(POST_EVOLVE, z4_to_adm,
         "fill in ADM variables from z4 variables");

  /* variables */
  AddVar("bssn_g",       "ij+ji", "conformal metric");
  AddVar("bssn_chi",     "",      "conformal factor psi^N");
  AddVar("bssn_A",       "ij+ji", "tracefree extrinsic curvature");
  AddVar("bssn_K",       "",      "trace of extrinsic curvature - 2Theta");
  AddVar("bssn_G",       "I",     "contracted Gamma");
  AddVar("z4_Theta",     "",      "time comp. of Z^{mu}");
  AddVar("z4_Z",         "I",     "spacial comp. of Z^{mu}");

  AddVar("bssn_trA",     "",      "trace of bssn_Aij");
  AddVar("bssn_detg"  ,  "",      "determinant of bssn_gij");
  AddVar("bssn_ginvtmp", "IJ+JI", "inverse metric, temporary variable");

  /* parameters */
  AddPar("z4_chi_psipower",   "-4", "chi = psi^N, N=z4_chi_psipower");
  AddPar("z4_chi_div_floor",  "-1000.", "use max(chi, chi_div_floor) in non-differentiated chi");
  
  AddPar("z4_forceKzero",     "no",  "set K identically to zero");
  AddPar("z4_subtractA",      "yes", "set trace of A identically zero");
  AddPar("z4_normalizedetg",  "yes",  "normalize determinant of gamma to one");
  AddPar("z4_register_algcon","yes",  "algebraic constraint handling after evolution step [no,yes;store]");
  
  if (Getv("z4_register_algcon", "yes"))
    AddFun(EVOLVE, z4_algcon_wrapper, "handle algebraic constraints");

  AddPar("z4_initial_lapse",    "donothing", "initial lapse [donothing,one]");
  AddPar("z4_lapse",            "constant", "lapse equation [constant,1+log,1+log2,harmonic;withshift]");
  AddPar("z4_lapsepsipower",    "0", "power of psi in lapse equation");
  AddPar("z4_lapseharmonicf",   "2",	 "2 for standard 1+log, ignored for harmonic");

  AddPar("z4_initial_shift",    "donothing", "initial shift [donothing,zero]");
  AddPar("z4_shift",            "withShiftadv",  "shift equation [gamma0,gamma2; withGadv; withB; withShiftadv]");
  AddPar("z4_shiftalphapower",  "0", "power of lapse in shift equation");
  AddPar("z4_shiftgammacoeff",  "1.0", "coefficient in shift equation (hardcoded=1 for no B shift)");
  AddPar("z4_shiftdriver",      "0", "coefficient of diffusion term");

  AddPar("z4_kappa1",           "0",  "the kappa1 from David's z4");
  AddPar("z4_kappa2",           "0",  "the kappa2 from David's z4");

  AddPar("z4_shift_use_B",      "no", "use B for beta evolution or not");
  AddPar("z4_bc_use_eta",       "no", "use eta term inside background.m");
  

  AddPar("z4_zeroThetaRHS",     "no", "");
  AddPar("z4_useNPTs",          "no", "");

  AddPar("z4_bc_psi0",          "no"   , "specify incoming radiation (psi0) at boundary ");  
  AddPar("z4_bc_psi0_a",        "0.001", "par a of incoming data");    
  AddPar("z4_bc_psi0_b",        "0.1"  , "par b of incoming data");    
  AddPar("z4_bc_psi0_c",        "200"  , "par c of incoming data");    
}

/* bam_bssn.c */
/* Bernd Bruegmann 6/02 */
/* Wolfgang Tichy  4/2004 */

#include "bam.h"
#include "bssn.h"




void bam_bssn(void) 
{
  if (!Getv("physics", "bssn")) return;
  printf("Adding bssn\n");

  /* functions */
  AddFun(POST_INITIALDATA, bssn_startup, 
	 "initialize bssn system from adm initial data");
  AddFun(POST_EVOLVE, bssn_to_adm,
         "fill in ADM variables from BSSN variables");

  /* variables */
  AddVar("bssn_g",        "ij+ji", "conformal metric");
  AddVar("bssn_chi",      "",      "conformal factor psi^N");
  AddVar("bssn_A",        "ij+ji", "tracefree extrinsic curvature");
  AddVar("bssn_K",        "",      "trace of extrinsic curvature");
  AddVar("bssn_G",        "I",     "contracted Gamma");
  
  AddVar("bssn_trA",      "",      "trace of bssn_Aij");
  AddVar("bssn_detg",     "",      "determinant of bssn_gij");
  AddVar("bssn_ginvtmp",  "IJ+JI", "inverse metric, temporary variable");
  
  AppPar("ExitIfNAN_vars","bssn_gxx bssn_chi bssn_K bssn_Axx bssn_Gx");

  /* parameters */
  AddPar("bssn_chi_psipower", "-4", "chi = psi^N, N=bssn_chi_psipower");
  AddPar("bssn_chi_div_floor","-1000.", "use max(chi, chi_div_floor) in non-differentiated chi");
  
  AddPar("bssn_forceKzero",       "no",  "set K identically to zero");
  AddPar("bssn_subtractA",        "yes", "set trace of A identically zero");
  AddPar("bssn_normalizedetg",    "yes", "normalize determinant of gamma to one");
  AddPar("bssn_register_algcon",  "yes", "register algebraic constraint handling after evolution step [no,yes;store]");
  if (Getv("bssn_register_algcon", "yes"))
    AddFun(EVOLVE, bssn_algcon_wrapper, "handle algebraic constraints");

  AddPar("bssn_initial_lapse", 
         "donothing", "initial lapse [donothing,one,onepluswave]");
  AddPar("bssn_lapse", "constant", 
         "lapse equation [constant,1+log,1+log2,harmonic;withshift]");
  AddPar("bssn_lapseharmonicf", "2",
	 "2 for standard 1+log, ignored for harmonic");
  AddPar("bssn_subtractK0",     "no", 
	 "for harmonic and 1+log, subtract initial traceK or not");
  
  AddPar("bssn_initial_shift", 
         "donothing", "initial shift [donothing,zero]");
  AddPar("bssn_shift", "constant", 
         "shift equation [constant,gamma0,gamma2;withGadv]");
  AddPar("bssn_shiftalphapower", "0", "power of lapse in shift equation");
  AddPar("bssn_shiftgammacoeff", "0.75", "coefficient in shift equation");
  AddPar("bssn_shiftgammaalphapower", "0", "coefficient in shift equation");
  AddPar("bssn_shiftdriver",     "0", "coefficient of diffusion term");

  
  /* eta stuff */  
  AddVar("bssn_eta", "", "loc. dep. shiftdriver");
  AddPar("bssn_use_eta", "no",
         "whether we want location dependent shiftdriver or not");
  if (Getv("bssn_use_eta","yes")) {
    AddPar("bssn_eta_useswitch", "0","0 for psim2, 1 for alpha, 2 for tr K_ij");
    AddPar("bssn_eta_width",
            "1", "test width of gaussian and rational fct. in eta, also used as C in Alcubierre's formula");
    AddPar("bssn_eta_epunc",
            "2.8", "eta_punc parameter in Eq. (3.2) of Alcubierre et al. gr-qc/0411137");
    AddPar("bssn_eta_einf",
            "2.8", "eta_infinity parameter in Eq. (3.2) of Alcubierre et al. gr-qc/0411137");
    AddPar("bssn_eta_rpower",
            "1", "power of r^2 in Eq. etas6 of Mov_Punv_Eta paper");
    AddPar("bssn_eta_width1",
            "1", "test width of rational fct. in eta");
    AddPar("bssn_eta_width2",
            "1", "test width of rational fct. in eta");
    //  used with bssn_eta_useswitch = 9
    AddPar("bssn_eta_rpower0", "3", "isotropic exponent in eta falloff");
    AddPar("bssn_eta_eta0",    "2", "local dimensionless eta around each black hole");
    AddPar("bssn_eta_eta1",    "1", "asymptotic dimensionless eta");
  }

  
  AddPar("bssnRHSto0","no","");
  AddPar("bssnsourceRHSto0","no","");
}




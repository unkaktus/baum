/* bam_BBIDdataReader.c */
/* mth 11/11 */

#include "bam.h"
#include "BBIDdataReader.h"



void bam_BBIDdataReader() 
{
  if (!Getv("physics", "BBIDdataReader")) return;
  printf("Adding BBIDdataReader\n");

  
  AddPar("BBID_modus","NSNS","");
  
  
  // additional variables used for the boost
  AddVar("BBID_psi4",      "",      "psi4 ... only used inside this project");
  AddVar("BBID_dpsi4",     "i",     "derivative of psi4 ... only used inside this project");
  AddVar("BBID_dalpha",    "i",     "derivative of alpha ... only used inside this project");
  AddVar("BBID_Tg",        "ij+ji", "metric");
  AddVar("BBID_TK",        "ij+ji", "extrinsic curvature");
  AddVar("BBID_Talpha",    "",      "lapse");
  AddVar("BBID_Tbeta",     "I",     "shift");
  AddVar("BBID_Trho",      "",      "primitive rho");
  AddVar("BBID_Tepsl",     "",      "primitive epsilon");
  AddVar("BBID_Tv",        "I",     "primitive v");
  AddVar("BBID_Tut",       "", "time component of fluid velocity");
  
  /* set initialdata (punc/TOV) and boost */
  AddFun(INITIALDATA_SET, set_binary_boost, "set binary boosted TOV/BH");
  AddPar("BBID_tov_integrate","h",        "variable to integrate");
  AddPar("BBID_tov_npts",     "100000",   "maximal number of iterations");

  AddPar("BBID_tov_hc",       "0.227869", "central entalpy");
  AddPar("BBID_tov_rhoc",     "0.00128",  "central density");
  AddPar("BBID_tov_R0",       "10.",      "a guess for the radius");
  // 14.01.2013 Added following pars for specifying 2 different TOV
  // - they are employed only in the routine set_TOV()
  // - in all the other cases the (unnumbered) pars above ("old") are employed
  // - the (unnumbered) pars above are kept for backwards compatibility at the moment
  // - this will create confusion, so someone will have to fix it at some point
  AddPar("BBID_tov_hc1",       "0.227869", "central entalpy");
  AddPar("BBID_tov_rhoc1",     "0.00128",  "central density");
  AddPar("BBID_tov_R01",       "10",      "a guess for the radius");
  AddPar("BBID_tov_hc2",       "0.227869", "central entalpy");
  AddPar("BBID_tov_rhoc2",     "0.00128",  "central density");
  AddPar("BBID_tov_R02",       "10",      "a guess for the radius");

  AddPar("BBID_tov_perturb",        "no",       "p,v,vrot");
  AddPar("BBID_tov_perturb_l",      "2",       "");
  AddPar("BBID_tov_perturb_m",      "0",       "");
  AddPar("BBID_tov_perturb_n",      "1",       "");
  AddPar("BBID_tov_perturb_lambda", "0.01",       "");
  AddPar("BBID_tov_perturb_lambda2","0.01",       "");
  


  AddPar("BBID_solve","no", "do we want to solve something");
if (Getv("BBID_solve","yes") && !(checkpoint_checkforfiles("_previous"))) { //don't need to solve if checkpoints are there
    
    AddFun(INITIALDATA_FINISH, BBID_solve, "solve binary boosted TOV/BH");
    
    AddPar("BBID_solve_itmax", 		"1", 		"maximal number of iterations");
    AddPar("BBID_solve_tolerance", 	"1e-6",		"tolerance in l2-norm of residuum");
    AddPar("BBID_Killing", 		"helliptical",	"Assumed approx Killing vector: helliptical, helical");	
    AddPar("BBID_solve_Omega", 		"0", 		"Orbital angular frequency");
    AddPar("BBID_rb", 			"1", 		"Rescaling parameter");
    AddPar("BBID_Econst_0", 		"0.9", 		"Energy constant for first integral in star 0");
    AddPar("BBID_Econst_1", 		"0.9", 		"Energy constant for first integral in star 1");
    AddPar("BBID_solve_soft_ell", 	"0.25", 	"Softening param for elliptic eqs. Gives percentage of new value");
    AddPar("BBID_solve_soft_eta", 	"0.25", 	"Softening param for density. Gives percentage of new value");
    AddPar("BBID_solve_com",		"0.",		"Center of Mass");
    AddPar("BBID_solve_spin", 		"corotational", "Specify corotational or irrotational data");
    AddPar("BBID_solve_ecc",		"0",		"Eccentricity of the orbit");
    AddPar("BBID_solve_fix_at_center",	"eta",		"Fix either eta or M");
    AddPar("BBID_solve_M_central",	"1.625002756",	"Value of M if fix_at_center is M");
    AddPar("BBID_solve_load", 		"yes",		"whether to load the initial data [yes/no]");
    AddPar("BBID_solve_loadDir"	,	"",		"specify directory to load initial data from");
    
    AddVar("BBID_energyconst",       "", "injection energy"); 
    AddVar("BBID_omega_old",       "", "variable for softening storage");  
    AddVar("BBID_apsi_old",       "", "variable for softening storage");  
    AddVar("BBID_beta_old",       "I", "variable for softening storage");  

    
    AddVar("BBID_u_beta",        "I", "temp vector for shift decomp"); //Not true
    AddVar("BBID_u_psi",        "", "temp scalar for shift decomp");	//Not needed
    AddVar("BBID_u_omega",      "", "rescaled psi");
    AddVar("BBID_u_apsi",     "", "apsi");

    AddVar("BBID_v_beta",        "I", "temp vector for shift decomp");
    AddVar("BBID_v_psi",        "", "temp scalar for shift decomp");
    AddVar("BBID_v_omega",      "", "psi");
    AddVar("BBID_v_apsi",     "", "apsi");

    AddVar("BBID_f_beta",        "I", "temp vector for shift decomp");
    AddVar("BBID_f_psi",        "", "temp scalar for shift decomp");
    AddVar("BBID_f_omega",      "", "psi");
    AddVar("BBID_f_apsi",     "", "apsi");

    AddVar("BBID_r_beta",        "I", "temp vector for shift decomp");
    AddVar("BBID_r_psi",        "", "temp scalar for shift decomp");
    AddVar("BBID_r_omega",      "", "psi");
    AddVar("BBID_r_apsi",     "", "apsi");

    AddVar("BBID_rho",		"", "density");
    AddVar("BBID_Si",		"I", "traceless extrinsic curvature");
    AddVar("BBID_S",		"", "traceless extrinsic curvature");
    AddVar("BBID_ut",       "", "traceless extrinsic curvature");
    AddVar("BBID_rhoH",       "", "traceless extrinsic curvature");
    AddVar("BBID_AA",       "", "traceless extrinsic curvature");
    
    AddVar("BBID_eta", 		"", "enthalpy");
    AddVar("BBID_psi",	"", 	"rescaled psi");
    
   AddVar("robinindex", "", "index for robin boundary");
   AddVar("robinflag",  "", "temporary flag for robin boundary");
  
  }
  
  
  /* test of EoS, compute MvsR */
  AddPar("BBID_tov_MvsR","no", "do we want to compute M(R)?");
  AddPar("BBID_tov_MvsR_p", "no", "compute M(R)?");

  if (Getv("BBID_tov_MvsR","yes")) {
    AddFun(PRE_PRE_PRE_INITIALDATA, test_tov_MvsR, "compute M(R)");
    AddPar("BBID_tov_rho0",       "1e-6",      "");
    AddPar("BBID_tov_drho",       "1.1",       "rho_i = rho0*drho^i");
    AddPar("BBID_tov_N",          "120",       "# of i");
    AddPar("BBID_tov_file",       "tov_MvsR.dat",  "file");
  }
  if (Getv("BBID_tov_MvsR_p", "yes")){
    AddFun(PRE_PRE_PRE_INITIALDATA, test_tov_MvsR_p, "compute M(R) with EOS(p)");
    AddPar("BBID_tov_p0",        "1e-6",       "");
    AddPar("BBID_tov_dp",        "1.1",       "p_i = p0*dp^i");
    AddPar("BBID_tov_N",         "120",       "# of i");
    AddPar("BBID_tov_file",       "tov_MvsR.dat",   "file");
  }

    AddPar("BBID_solve_external","no", "If you want to use e.g. punctures to solve the BH?"); 

   if (Getv("BBID_solve_external", "yes")) {      
         
             AddVar("incl_punctures_u",        "", "correction function from punctures"); 
             AddPar("BBID_solve_external_method","multigrid", "Method within puncture"); 

      if (Getv("BBID_solve_external_method", "spectral")) {
          AddPar("punctures_ps_order", "8", "order of interpolation from spectral data"); 
          AddPar("punctures_ps_file", " ", "if set, read spectral data from this file");
          AddPar("punctures_ps_interpol_check", "none","if yes: NR-method for inverse transformation");
      
               AddFun(PRE_GRID, BBID_read_ps_parameters,"read spectral puncture parameters from file");
  }
   

                    }
  
  /* set it here automatically */
  AddPar("mass1", "1.0", "mass");
}


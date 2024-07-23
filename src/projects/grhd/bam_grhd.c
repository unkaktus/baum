/* bam_matter.c */

#include "bam.h"
#include "grhd.h"


void bam_grhd(void){

  if (!Getv("physics", "grhd")) return;
  printf("Adding grhd\n");
  
  // functions
  AddFun(PRE_PRE_INITIALDATA, grhd_startup , "grhd_startup");

  // variables 
  AddVar("grhd_D"   ,"" , "conservative D   * sqrt(gamma)");
  AddVar("grhd_Tau" ,"" , "conservative tau * sqrt(gamma)");
  AddVar("grhd_S"   ,"i", "conservative S_i * sqrt(gamma)");

  AddVar("grhd_rho" ,"" , "primitive rho");
  AddVar("grhd_epsl","" , "primitive epsilon");
  AddVar("grhd_v"   ,"I", "primitive v");
  AddVar("grhd_p"   ,"" , "primitive p");
  AddVar("grhd_v2"  ,"" , "primitive v2");
  
  if (Getv("conservative_amr","yes")){ 
  AddVar("camr_D"   ,"" , "camr correction D   * sqrt(gamma)");
  AddVar("camr_Tau" ,"" , "camr correction tau * sqrt(gamma)");
  AddVar("camr_S"   ,"i", "camr correction S_i * sqrt(gamma)");
  }
  
  //turbulence
  AddPar("grhd_use_turbulence", "no", "yes or no");
  AddPar("grhd_turb_lmix", "0.01", "characteristic turbulence length"); //++ upgrade from par to EOS output
  
  if(Getv("grhd_use_turbulence","yes")){
  AddVar("grhd_turbTau", "ij+ji", "turbulent addition to stress energy tensor");
  }
   
  AppPar("ExitIfNAN_vars","grhd_D grhd_Tau grhd_Sx grhd_rho grhd_epsl grhd_vx grhd_p grhd_v2");
  
  // initialdata 
  AddPar("grhd_use_initialdata_p", "no",  "use initialdata presure or compute by EoS");

  // C2P
  AddPar("grhd_C2P",                 "p_root",  "con2prim routine to be used");
  AddPar("grhd_C2P_NewtonRaphsonTR", "1e-10",   "accuracy of the root finder (con2prim)");
  AddPar("grhd_C2P_NewtonRaphsonNR", "100",     "maximal number of iterations");

  // sources
  AddPar("grhd_source_computation", "prim", "use prim/cons vars to compute source terms");
  AddPar("grhd_cooling", "0", "if we use simple cooling: arxiv1208.5487v2. 0: no 1: yes");
  AddPar("grhd_cooling_time", "1e15", "cooling time parameter");

  AddPar("grhd_vmax",      "0.999", "maximal velocity in fluxes");
  AddPar("grhd_Wlor_max",  "1e8"  , "Wlor corresponding to vmax");
  

  // ATM
  AddPar("grhd_use_atmosphere", "ColdStatic" , "no,ColdStatic,Vacuum");
  AddPar("grhd_atm_factor"    , "100."       , "factor between atm and value where to set atm");
  AddPar("grhd_atm_level"     , "1e-9",        "treshold for cutting matter variables to set atmosphere");

  AddPar("grhd_use_atmosphere_mask"   , "no" , "matter mask");
  AddPar("grhd_use_atmosphere_prerhs" , "yes", " ");
  AddPar("grhd_use_atmosphere_postrhs", "no" , " ");

  /* Vacuum */
  AddPar("grhd_vacuum_prerhs_set", "upre", "for grhd_Vac_set_prerhs "
         "VL in which we set vacuum. Any combination of: upre, ucur");
  AddPar("grhd_vacuum_set", "cons" , "for grhd_Vac_set_u "
         "vars in which we set vacuum. Any combination of: cons, prim");
  AddPar("grhd_vacuum_limit_grhd_D", "no", "set to Vacuum in "
         "grhd_Vac_set_u, if grhd_D is too high inside horizon [no,yes]");
  AddPar("grhd_vacuum_grhd_D_fac", "100", "set to Vacuum in grhd_Vac_set_u "
         "if grhd_D > (grhd_D_fac)*SQRTdetg*Wlor*grhd_rho, "
         "alpha < horizonalpha and grhd_vacuum_limit_grhd_D=yes");
  AddPar("grhd_vacuum_horizonalpha", "0.2", "set to Vacuum in grhd_Vac_set_u "
         "if grhd_D > (grhd_D_fac)*SQRTdetg*Wlor*grhd_rho, "
         "alpha < horizonalpha and grhd_vacuum_limit_grhd_D=yes");

  // excision 
  AddPar("grhd_use_excision"   , "no"  , "excise (do not evolve) matter mask");
  AddPar("grhd_excision_modus" , "atm" , "{cut,lin,atm,vmax,fixr}");
  AddPar("grhd_excision_rmin"  , "0.0" , "excise r<rmin");
  AddPar("grhd_excision_rmax"  , "20"  , "excise r>rmax");
  AddPar("grhd_excision_rfct"  , "0.3" , "excise r<rah*rfct");
  AddPar("grhd_excision_fact"  , "0."  , "multiplicative factor for excision vars");

  // interpolation schemes in R/P
  if(Getv("conservative_amr","yes")) {
     AppPar("matter_interpolate_vars" , "grhd_D grhd_S grhd_Tau"); 
     AppPar("camr_interpolate_vars"   , "camr_D camr_S camr_Tau");}
  else   {
     AppPar("matter_interpolate_vars"  , "grhd_D grhd_S grhd_Tau");
     AppPar("camr_interpolate_vars"  , "");}

  /* how to reconstruct velocity and epsl */
  AddPar("grhd_recvel" , "vx" ,   "vx,bx,rhoWv");
  AddPar("grhd_recepsl", "epsl" , "epsl,rhoepsl");
  AddPar("grhd_HLLC_K", "1.0", "parameter for contact discontinuity detection");
  AddPar("grhd_HLLC_mach_lim", "0.1", "mach number at which start applying correction");

  // checkpointing 
  if (Getv("checkpoint","yes")) {
    AppPar("checkpoint_additional_variables","grhd_p");
  }

  
}



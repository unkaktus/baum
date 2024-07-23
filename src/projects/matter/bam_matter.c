/* bam_matter.c */

#include "bam.h"
#include "matter.h"


void bam_matter(void) 
{
  if (!Getv("physics", "matter")) return;
  printf("Adding matter\n");
  
  if (Getv("physics", "grhd")+Getv("physics", "grhdY")+Getv("physics", "grhdDY")+Getv("physics", "grmhd")+Getv("physics", "grmhdY")+Getv("physics", "grhd2Fluid")!=1) 
         errorexit("add physics: grhd, grhdY, grhdDY, grmhd, grmhdY or grhd2Fluid");

  if (!(Getv("physics", "eos")))  errorexit("add physics: eos");
  

  // functions 
  AddFun(PRE_PRE_INITIALDATA, hrsc_startup, "hrsc_startup");
  AddFun(PRE_PRE_INITIALDATA, matter_startup, "matter_startup");
  AddFun(POST_POST_INITIALDATA, matter_init, "matter_init");
  
  // variables 
  AddVar("matter_mask","" , "bitmask vacuum/matter");

//   // check initial data are set 
//   if (Getv("physics","BBIDdataReader") +
//       Getv("physics","RNSdataReader") +
//       Getv("physics","BNSdataReader") +
//       Getv("physics","LIDdataReader") +
//       Getv("physics","BHNSLIDdataReader") +
//       Getv("physics","Michel") +
//       Getv("physics","shocktube")+ 
//       Getv("physics","srhdtestID") != 1)
//     errorexit("you need exactly one intialdata-project!");
  
  // HRSC 
  AddPar("hrsc_nghosts"         ,"4"     , "number of ghost points for 1d rec ");
  AddPar("hrsc_rec"             ,"LINTVD", "CONST,AVG,LINTVD,PPM,CENO3,WENO5,WENOZ");
  AddPar("hrsc_rec_metric"      ,"LAG4"  , "CONST,LINTVD,AVatter_rhs+G,CENO3,LAG4,WENO5,WENOZ,LAG6");
  AddPar("hrsc_TVD_limiter"     ,"MC2"   , "MM2,MC2");
  AddPar("hrsc_flux"            ,"LLF"   , "LLF,HLL,HO,HO_LLF");
  AddPar("hrsc_flux_switch_rho" ,"5.0"   , "density factor where we employ LLF rho < switch*atm_rho*atm_fac");
  AddPar("hrsc_flux_switch_minrho" ,"0.001",
         "if no atmopshere, switch between HO/LLF at rho=switch*RHOMAX");
  AddPar("hrsc_flux_switch_alpha" ,"0.2"   , "minimum lapse under which we employ LLF instead of HOLLF");
  // HG: adding the safety reconstruction methods
  AddPar("hrsc_rec_safety"      ,"no", "LINTVD");
  AddPar("hrsc_rec_tracer"      ,"LINTVD", "CONST, AVG, LINTVD, PPM, CENO3, WENO5, WENOZ");

  AddPar("hrsc_use_lower"       ,"no", "use lower order for lower densities [yes/no]");
  if (Getv("hrsc_use_lower","yes")) {
   AddPar("hrsc_rec_low"        ,"LINTVD", "falls back on this scheme");
   AddPar("hrsc_rec_lowrho"     ,"1e-10",  "the value below which lower order scheme should be used");
   AddPar("hrsc_rec_lowrho_fatm","-1",  "the value below which lower order scheme should be used as factor of the atmosphere");
  }
  
  //conservative amr
   AddPar("conservative_amr","no", "use conservative amr (yes/no)"); 
  if (Getv("conservative_amr","yes")) {
     AddVar("camr_mask_A","" , "bitmask parent cells (A)");
     AddVar("camr_mask_B","" , "bitmask child cells (B)");
     AddPar("camr_level_min","0", "min level, where c-amr is active");
     AddPar("camr_treshold","1.e5", "do not correct if fluxes are to different");
     AddPar("camr_atm_cut","1e-15", "do not correct rho < camr_atm_cut");
     AddPar("camr_time","-1", "when to activate camr");

     AddFun(POST_MOVEBOX,c_amr_mask, "recompute mask after box moved");

     if(Getd("camr_time")>-0.5) AddFun(POST_ANALYZE, c_amr_time , "activate camr");
  }

  /* Should matter be set to zero on level0, where outer boundary is?
     Usually rhs_matter_null points to matter_rhs_set_null */
  AddPar("rhs_matter_null_ONlevel0", "no",
         "call rhs_matter_null in evolve_rhs [no,yes]");
}






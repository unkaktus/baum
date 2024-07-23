/* bam_eos.c */
/* mth 06/12 */

#include "bam.h"
#include "eos.h"




void bam_eos() 
{
  if (!Getv("physics", "eos")) return;
  printf("Adding eos\n");

  // functions
  AddPar("eos_ID", "no" , "is EOS provided by Initial data project?");
  
  if (Getv("eos_ID","yes")) AddFun(PRE_INITIALDATA, eos_startup, "eos_startup");
  else AddFun(PRE_PRE_PRE_INITIALDATA, eos_startup, "eos_startup");

  //if (Getv("eos","tab3d") && 0)
  //  AddFun(PRE_PRE_PRE_INITIALDATA, eos_tab3d_test, "test eos table");
  
  AddPar("eos",             "none", "EoS");
  AddPar("eos_K"    ,       "100",  "polytropic constant (poly EOS)");
  AddPar("eos_Gamma",       "2"  ,  "polytropic index (ideal gas EOS)");
  AddPar("eos_tab_file",    "none", "EOS tab or pwp file");

  AddPar("eos_tab_file_T0",  "none", "EOS T=0 tab (for 3d tab)");
  
  AddPar("eos_tab_file_beta",  "none", "EOS tab in beta equilibrium"); //HG

  AddPar("eos_tab_file_derivs", "none", "EOS 3d tab with derivs for HO flux"); //HG

  AddPar("eos_Y"    ,       "0.4", "fix constant Ye value for tab3d");

  AddPar("eos_T"    ,       "0.01", "fix constant T value for atm");
 
  AddPar("eos_interp",   "steffen", "cspline/steffen/linear");

  //AddPar("eos_pwp_ID", "none", "pwp used for the ID");

  AddPar("eos_compo", "npe", "npe/npeH/hydro");

  AddPar("grrhd_m1_nls_rates", "no", "use nls rates tab for M1?");

  // EoS checks 
  if (Getv("eos", "none")) 
    errorexit("you have to set an EOS!!!");
  
  if (Getv("eos","tab1d")+Getv("eos","tab1dhot")+
      Getv("eos","pwp")+Getv("eos","pwphot")+Getv("eos","pwphotWT")+
      Getv("eos","tab3d")+Getv("eos","nuclear")+Getv("eos","tab1d_hot")+Getv("eos","tab1d_cold") 
      && Getv("eos_tab_file","none") && !Getv("eos_ID","yes"))
    errorexit(" please pass an EoS table setting eos_tab_file");
  
  
  
}

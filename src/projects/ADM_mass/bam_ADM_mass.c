/* bam_ADM_mass.c */
/* Jose Gonzalez 9/04 */

#include "bam.h"
#include "ADM_mass.h"


void bam_ADM_mass() 
{
  if (!Getv("physics", "ADM_mass")) return;
  printf("Adding ADM_mass\n");

  /* functions: */
  AddFun(ANALYZE, ADM_mass, "ADM_mass");

  /* variables */
  AddVar("ADM_mass_integrand",      "", "Integrand of ADM Mass");
  AddVar("ADM_mass_integrand_chi",  "", "ADM Mass using only chi");
  AddVar("ADM_mass_Pxint",          "", "Integrand of Px");
  AddVar("ADM_mass_Pyint",          "", "Integrand of Py");
  AddVar("ADM_mass_Pzint",          "", "Integrand of Pz");  
  AddVar("ADM_mass_Jxint",          "", "Integrand of Jx");
  AddVar("ADM_mass_Jyint",          "", "Integrand of Jy");
  AddVar("ADM_mass_Jzint",          "", "Integrand of Jz");
  
  /* parameters */
  AddPar("ADM_mass_persist",    "no", "whether additional memory persists");

  AddPar("ADM_mass_ncircles",   "101", "Number of circles");
  AddPar("ADM_mass_npoints",    "80", "Number of points");
  AddPar("ADM_mass_findevery",  "1", "How often to compute the ADM mass");

  AddPar("ADM_mass_lmin",       "0",    "compute ADM mass for levels >= lmin");
  AddPar("ADM_mass_lmax",       "1000", "compute ADM mass for levels <= lmax");
  AddPar("ADM_mass_r",          "20 25 30", "list of radii for ADM mass");


}




/* bam_Invariants.c */
/* Bernd Bruegmann 12/2005 */
/* Jose Gonzalez 01.06 */

#include "bam.h"
#include "Invariants.h"




void bam_Invariants() 
{
  if (!Getv("physics", "Invariants")) return;
  printf("Adding Invariants\n");

  AddFun(ANALYZE, compute_curvature_invariants,
	 "compute curvature invariants");
  AddPar("invariants_order", "4",
	 "order of finite differencing, use to override default");
  AddPar("invariants_EB", "no",
         "whether to compute invariants from E and B, [yes,no]");
  AddPar("invariants_persist", "no", "whether additional memory persists");
  
  AddPar("invariants_tetrad",      "GS", 
         "tetrad to use [GS,kinnersley]");
  AddPar("gauss_codacci_mainardi", "standard",
         "Version of Gauss-Codacci, Mainardi equations [flat,standard]");
  AddPar("invariants_rescale",     "psi4", 
         "rescale to cancel fall-off");
  
  AddPar("invariants_ntheta",      "47", "# points in theta direction on sphere");
  AddPar("invariants_nphi",        "46", "# points in phi direction on sphere");
  
  
  AddPar("invariants_compute_modes",  "yes", "compute modes");
  AddPar("invariants_modes_lmin",     "0",   "compute modes for levels >= lmin");
  AddPar("invariants_modes_lmax",     "100", "compute modes for levels <= lmax");
  AddPar("invariants_modes_r",        "20 25 30 50 100", 
         "genral list of radii psi4 analysis (modes energy r-output)");
  AddPar("invariants_modes_mode_lmin","2",   "compute modes for l >= mode_lmin");
  AddPar("invariants_modes_mode_lmax","4",   "compute modes for l <= mode_lmax");
  AddPar("invariants_modes_output",   "all", "what modes do you want to output (example: 2,0  2,2  2,-2)");
  
  AddPar("invariants_compute_energy", "yes", "compute modes");
  
  AddPar("invariants_shells_output",  "no", "output all shells data");
  
  
  AddPar("invariants_x0", "0", "origin of coordinate sphere");
  AddPar("invariants_y0", "0", "origin of coordinate sphere");
  AddPar("invariants_z0", "0", "origin of coordinate sphere");
  AddPar("invariants_vx0", "0", "speed of coordinate sphere");
  AddPar("invariants_vy0", "0", "speed of coordinate sphere");
  AddPar("invariants_vz0", "0", "speed of coordinate sphere");

  AddPar("invariants_analytic_Psi4","no",
         "compute analytic Psi4 (Teukolsky wave, Kinnersley tetrad)");
  AddPar("invariants_xana", "1.0", 
         "x coordinate of extraction point for analytic");
  AddPar("invariants_yana", "1.0", 
         "y coordinate of extraction point for analytic");
  AddPar("invariants_zana", "1.0", 
         "z coordinate of extraction point for analytic");

  
  /* variables */
  AddVar("rpsi0",    "",      "rpsi0");
  AddVar("ipsi0",    "",      "ipsi0");
  AddVar("rpsi1",    "",      "rpsi1");
  AddVar("ipsi1",    "",      "ipsi1");
  AddVar("rpsi2",    "",      "rpsi2");
  AddVar("ipsi2",    "",      "ipsi2");
  AddVar("rpsi3",    "",      "rpsi3");
  AddVar("ipsi3",    "",      "ipsi3");
  AddVar("rpsi4",    "",      "rpsi4");
  AddVar("ipsi4",    "",      "ipsi4");
  AddVar("rI",       "",      "rI");
  AddVar("iI",       "",      "iI");
  AddVar("rJ",       "",      "rJ");
  AddVar("iJ",       "",      "iJ");
  AddVar("rS",       "",      "rS");
  AddVar("iS",       "",      "iS");
  AddVar("Csqr",     "",      "Csqr");
  
  /* integrands for E, P, J computation (temporary) */
  AddVar("integrand_dEdt","","integrand of dEdt");
  AddVar("integrand_dPxdt","","integrand of dPxdt");
  AddVar("integrand_dPydt","","integrand of dPydt");
  AddVar("integrand_dPzdt","","integrand of dPzdt");
  AddVar("integrand_dJzdt","","integrand of dJzdt");
  AddVar("integrand_One","","integrand for sphere norm"); // XXX quite hacky, not nice
  AddVar("integrand_dA","","integrand of dA, surface area element");

  /* variables for E, P, J computation */
  AddVar("int_rpsi4","","integral in time of rpsi4");
  AddVar("int_ipsi4","","integral in time of ipsi4");
  AddVar("int2_rpsi4","","2nd. integral in time of rpsi4");
  AddVar("int2_ipsi4","","2nd. integral in time of ipsi4");
  
  AddVar("dint_rpsi4dphi","","derivative of the integral in time of rpsi4 with respect to phi");
  AddVar("dint_ipsi4dphi","","derivative of the integral in time of ipsi4 with respect to phi");

  AddVar("dint2_rpsi4dphi","","derivative of the integral in time of rpsi4 with respect to phi");
  AddVar("dint2_ipsi4dphi","","derivative of the integral in time of ipsi4 with respect to phi");

  AddVar("rpsi4_p","","prev rpsi4");
  AddVar("ipsi4_p","","prev ipsi4");

  AddPar("Invariants_output_time","-1","specify time for output explicitly");

  /* variable used for small array storage */
  AddVar("invariants_storage", "", "small array storage");


  /* fix symmetries, some of the above are pseudo scalars */
  VarInvertSymmetry(Ind("ipsi4"));
  /* unfinished ... */

  AddVarModes();
  
}

/* bam_hydroanalysis.c */
/* S Bernuzzi 10/12 */
/* T Dietrich 2014 and 2015*/
/* T Dietrich 2017 */
#include "bam.h"
#include "hydroanalysis.h"


void bam_hydroanalysis() 
{
  if (!Getv("physics", "hydroanalysis")) return;
  printf("Adding hydroanalysis\n"); 

  AddFun(POST_ANALYZE, hydroanalysis, "perform analysis of hydro variables");
  
  /* parameters */
  AddPar("hydroanalysis_lmin", "0", "compute analysis for levels >= lmin");
  AddPar("hydroanalysis_lmax", "20", "compute analysis for levels >= lmin");
  AddPar("hydroanalysis_rATM", "1e-15", "atmosphere (floor) value for r");
  AddPar("hydroanalysis_wMAX", "1e8.", "reset value for Lorentz factor if v>1");

  // get varnames from parfile 
  // => project works with different matter projects (grhd, grmhd, grhdY)
  AddPar("hydroanalysis_D", "grhd_D"   , "tell the code the name of D variable");
  AddPar("hydroanalysis_p", "grhd_p"   , "varname for p in h = 1 + e + p/r");
  AddPar("hydroanalysis_r", "grhd_rho" , "varname for r in h = 1 + e + p/r");
  AddPar("hydroanalysis_e", "grhd_epsl", "varname for r in h = 1 + e + p/r");
  AddPar("hydroanalysis_v", "grhd_vx"  , "varname for 3-velocity");

  /* vars */
  AddVar("hydroa_ut"   , "", "time component of 4-velocity (low index)" );
  AddVar("hydroa_Dh"   , "", "D outside horizon" );
  AddVar("hydroa_Db"   , "", "D bounded" );
  AddVar("hydroa_Du"   , "", "D unbounded" );
  AddVar("hydroa_vorx" , "", "vorticity-x" );
  AddVar("hydroa_vory" , "", "vorticity-y" );
  AddVar("hydroa_vorz" , "", "vorticity-z" );
  AddVar("hydroa_etot" , "", "total energy density of ejecta" );
  AddVar("hydroa_uesc" , "", "total internal energy density of ejecta" );
  AddVar("hydroa_Px" , "", "lin mom hydro-x" );
  AddVar("hydroa_Py" , "", "lin mom hydro-y" );
  AddVar("hydroa_Pz" , "", "lin mom hydro-z" );
  AddVar("hydroa_Pux" , "", "lin mom hydro ejecta-x" );
  AddVar("hydroa_Puy" , "", "lin mom hydro ejecta-y" );
  AddVar("hydroa_Puz" , "", "lin mom hydro ejecta-z" );

  /* tmp vars */
  AddVar("hydroa_ivorx", "", "function to be derived to obtain vorticity" );
  AddVar("hydroa_ivory", "", "function to be derived to obtain vorticity" );
  AddVar("hydroa_ivorz", "", "function to be derived to obtain vorticity" );

  /* mode-projection of D*/
  AddPar("hydroa_mode_projection", "no"   , "no/yes, decide if D_{l,m} projection is wanted");   
  AddPar("hydroanalysis_Mbar_radius", "0", "minimum radius inside which we compute the baryonic mass [0=off]");
  AddPar("hydroanalysis_Mbar_nradius", "3", "number of radii inside we compute the baryonic mass");
  AddPar("hydroanalysis_Mbar_dradius", "1", "spacing of the radius inside which we compute the baryonic mass");
  
  
  /********************************************/
  /* computation of ejected matter through coordinate spheres*/
  
  AddPar("hydroanalysis_ejecta_angular",        "no"    , "no/yes, decide whether to save the direction of the ejecta");    
  AddPar("hydroanalysis_ejecta_angular_radius", "300"   , "at what radius to save the ejecta direction");
  AddPar("hydroanalysis_radiation_angular_radius", "300", "at what radius to save the radiation direction");
  AddPar("hydroanalysis_ejecta_angular_dr",        "100"   , "at what radius to save the ejecta direction");
  AddPar("hydroanalysis_radiation_angular_dr",     "100", "at what radius to save the radiation direction");
  if (Getv("eos", "nuclear"))
    AddPar("hydroanalysis_ejecta_angular_output", "grhd_rho grhd_T grhd_Y hydroa_ut grhd_v2 hydroa_integrand_dDudt hydroa_integrand_dPxdt hydroa_integrand_dPydt hydroa_integrand_dPzdt", "ejecta output on sphere");
  else
    AddPar("hydroanalysis_ejecta_angular_output", "grhd_rho grhd_p grhd_epsl hydroa_ut grhd_v2 hydroa_integrand_dDudt hydroa_integrand_dPxdt hydroa_integrand_dPydt hydroa_integrand_dPzdt", "ejecta output on sphere");  

  AddPar("hydroanalysis_ejecta_spheres", "no"   , "no/yes, decide if ejecta should also be computed on spheres");     
  AddPar("hydroanalysis_ejecta_spheres_radius","100","minimum radius for ejecta computation"); 
  AddPar("hydroanalysis_ejecta_nradius","5","number of radii for which we compute the ejecta");
  AddPar("hydroanalysis_ejecta_dradius","100","spacing of the radius for which we compute the ejecta");   
  
  AddPar("hydroanalysis_sphere_ntheta","47","# points in theta direction on sphere");
  AddPar("hydroanalysis_sphere_nphi","46","# points in phi direction on sphere");
  AddPar("hydroanalysis_sphere_order","4","order of finite differencing, use to override default");  
  
  // integrands for Du,Db,P,Pu 
  AddVar("hydroa_integrand_dDudt","","integrand of dDudt");
  AddVar("hydroa_integrand_dDbdt","","integrand of dDbdt");
  AddVar("hydroa_integrand_dDbndt","","integrand of dDbdt but only for matter moving inwards");
  AddVar("hydroa_integrand_dPuxdt","","integrand of dPuxdt");
  AddVar("hydroa_integrand_dPuydt","","integrand of dPuydt");
  AddVar("hydroa_integrand_dPuzdt","","integrand of dPuzdt");
  AddVar("hydroa_integrand_dPxdt","","integrand of dPxdt");
  AddVar("hydroa_integrand_dPydt","","integrand of dPydt");
  AddVar("hydroa_integrand_dPzdt","","integrand of dPzdt");

  AddPar("hydroanalysis_radiation_spheres", "no"   , "no/yes, decide if radiation should also be computed on spheres");

  if (Getv("physics", "grrhd_m1")) {

//    AddPar("hydroanalysis_radiation_spheres", "no"   , "no/yes, decide if radiation should also be computed on spheres");

    AddVar("hydroa_integrand_dEnuedt","","integrand of dEnuedt");
    AddVar("hydroa_integrand_dEnuadt","","integrand of dEnuedt");
    AddVar("hydroa_integrand_dEnuxdt","","integrand of dEnuedt");

  }

    /* computation of magnetic field quantities (only for grmhd and grmhdY)*/
  if (Getv("physics", "grmhd") || (Getv("physics", "grmhdY"))){
    AddVar("hydroa_Emag","", "electromagnetic energy density");
    AddVar("hydroa_Pmag","", "magnetic pressure");
    AddVar("hydroa_divB","", "divergence of magnetic field");
    AddVar("hydroa_B","", "magnitude of magnetic field");
    AddVar("hydroa_Emagtor","", "electromagnetic energy density for torodial component");
    AddVar("hydroa_Emagpol","", "electromagnetic energy density for poloidal component");
    AddVar("hydroa_Btor","", "magnitude of toroidal magnetic field");
    AddVar("hydroa_Bpol","", "magnitude of poloidal magnetic field");
    AddVar("hydroa_Fpoy","", "Poynting flux on a sphere");
  } 

  // variable used for small array storage
  AddVar("hydroa_storage", "", "small array storage");
   
}

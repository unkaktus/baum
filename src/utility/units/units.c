/* units.c 
   mth 06/10 */

#include "bam.h"
#include "units.h"





void set_units(tU* u)
{
  /***************************************************************************/
  /* Constants in CGS 
   erg = g·cm²/s²
   dyn = g·cm/s²    -> erg/cm = dyn
  */
  /***************************************************************************/
  double clight  = 2.99792458e+10;    // [cm/s]
  double hplanck = 6.62606876e-27;    // [erg·s]
  double Gnewton = 6.67428e-08;       // [cm³/(g·s²)]
  double Msun    = 1.98892e+33;       // [g]
  double hn56Fe  = 1.4919335988e-3;   // []
  double mu      = 1.660538921e-24;   // [g] definition is based on C-12
  // NIST: Atomic Weights and Isotopic Compositions - Version History
  // Version 4.1  July 2015
  // http://physics.nist.gov/Comp
  // mFe = 55.93493633
  // Note: this mass includes the electrons!
  double mn56Fe  = 55.93493633/56 * mu;    // [g] (total mass / number of nucleons) in Fe-56

  
  // set constants in cgs
  u->clight_cgs  = clight;
  u->hplanck_cgs = hplanck; 
  u->Gnewton_cgs = Gnewton;
  u->Msun_cgs    = Msun;
  u->mn56Fe_cgs  = mn56Fe;
  u->hn56Fe_cgs  = hn56Fe;
  u->mu_cgs      = mu;

  // Dimensional factors for conversion from dimensionless "" units
  // (clight = Gnewton = Msun = 1) to cgs units 
  u->Length_cgs  = 0.5*2.953250077*1e5;                             // [cm]
  u->Time_cgs    = u->Length_cgs/clight;                            // [s]
  u->Mass_cgs    = Msun;                                            // [g]
  u->Energy_cgs  = Msun*clight*clight;                              // [erg]
  u->Surface_cgs = u->Length_cgs*u->Length_cgs;                     // [cm²]
  u->Volume_cgs  = u->Length_cgs*u->Length_cgs*u->Length_cgs;       // [cm³]
  u->Press_cgs   = Msun/(u->Time_cgs*u->Time_cgs*u->Length_cgs);    // [g/(s²·cm) = dyn/cm²]
  u->Mdens_cgs   = u->Mass_cgs/u->Volume_cgs;                       // [g/cm³]
  u->Edens_cgs   = u->Energy_cgs/u->Volume_cgs;                     // [g/(s²·cm) = erg/cm³]
  u->Amom_cgs    = u->Energy_cgs*u->Time_cgs;                       // [erg·s]
  
  // Dimension factor for conversion from  length to fermi 
  u->Length_fm = 1e+13*Gnewton*Msun/(clight*clight);                // []
  
  // other units
  u->Length_km    = u->Length_cgs * 1e-5;                           // [km]
  u->Density_kgm3 = u->Mdens_cgs * 1e-3 / pow( 1e-2 ,3);            // [kg/m³]
  u->Energy_MeV   = u->Energy_cgs / 1.602176565e-6;                 // [MeV]
  u->Volume_fm3   = u->Volume_cgs * pow( 1e13 ,3);                  // [fm³]
  
  //u->Mpc_cgs = 3.08568025*1e24; // [cm]    
   
}
























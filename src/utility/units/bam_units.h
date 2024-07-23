/* bam_units.h
   mth 2012 */




typedef struct tUNITS {
  
  double clight_cgs;      // speed of light
  double hplanck_cgs;     // plank constant
  double Gnewton_cgs;     // gravitational constant
  double Msun_cgs;        // mass of the sun
  double mn56Fe_cgs;
  double hn56Fe_cgs;      // [g] (total mass / number of nucleans) in Fe-56
  double mu_cgs;          // atomic mass unit
  
  double Length_cgs;
  double Time_cgs;
  double Mass_cgs;
  double Energy_cgs;
  double Surface_cgs;
  double Volume_cgs;
  double Press_cgs;
  double Mdens_cgs;
  double Edens_cgs;
  double Amom_cgs;
  
  double Length_fm;
    
  double Length_km;
  double Density_kgm3;
  double Energy_MeV;
  double Volume_fm3;
} tU;

typedef tU tUnits;

void set_units(tU* u);

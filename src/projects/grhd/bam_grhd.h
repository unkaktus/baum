/* bam_grhd.h */



typedef struct {
  
  // max vals for speeds in fluxes 
  double HRSC_VMAX, HRSC_WlorMAX;

  // atm parameters
  double ATM_RHOMAX, ATM_PMAX;
  double ATM_RHOATM,ATM_FATM,ATM_PATM,ATM_EPSLATM,ATM_YATM,ATM_TATM;
 
  // turbulence parameters
  double turb_lmix;

  // indices for other vars
  int IND_turbTau;

  double HLLC_K;
  double HLLC_mach_lim;

  // function pointers
  void (*deleptonization)();
  
} tGRHD;

/* infos about all specific function ...*/
extern tGRHD GRHD;














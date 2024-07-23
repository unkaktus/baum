/* grhd.c 
   sbernuz, mth 03/2012 */

#include "bam.h"
#include "grhd.h"

#define PR 0

/* infos about all specific function like grhd ...*/
tGRHD GRHD;



/* settings called within the bam_grhd.c file
   this routine is used mainly to interface the options from parfile 
   and the grhd specific routine with the general setup of the matter* routines 
   and it is stored inside a global struct!
*/
int  grhd_startup(tL *level) 
{
  
  printf("Initializing grhd:\n");
  printf("  provide ghrd varlists and fkt pointers to matter\n");

  // set no vars for grhd
  MATTER.NVq   = MATTER.NVf = 5;
  MATTER.NVw   = 8;  // rho, epsl, v^i, p, v2, detg
  MATTER.NVg   = 11; // alpha, beta^i, g_{ij}
  MATTER.NVo   = 0;

  if (Getv("grhd_use_turbulence","yes")) MATTER.NVo += 6 ;


  // set (local) list with var names 
  // rem order matters !
  MATTER.q_names_list = (char**) malloc( MATTER.NVq * sizeof(char*) );
  MATTER.w_names_list = (char**) malloc( MATTER.NVw * sizeof(char*) );
  MATTER.o_names_list = (char**) malloc( MATTER.NVo * sizeof(char*) );

  MATTER.q_names_list[0] = "grhd_D";
  MATTER.q_names_list[1] = "grhd_Tau";
  MATTER.q_names_list[2] = "grhd_Sx";
  MATTER.q_names_list[3] = "grhd_Sy";
  MATTER.q_names_list[4] = "grhd_Sz";

  MATTER.w_names_list[0] = "grhd_rho";
  MATTER.w_names_list[1] = "grhd_epsl";
  MATTER.w_names_list[2] = "grhd_vx";
  MATTER.w_names_list[3] = "grhd_vy";
  MATTER.w_names_list[4] = "grhd_vz";
  MATTER.w_names_list[5] = "grhd_p";
  MATTER.w_names_list[6] = "grhd_v2";
  MATTER.w_names_list[7] = "adm_detg";

  int v = 0;

  if (Getv("grhd_use_turbulence","yes")) {
     GRHD.IND_turbTau = v;
     MATTER.o_names_list[v  ] = "grhd_turbTauxx";
     MATTER.o_names_list[v+1] = "grhd_turbTauxy";
     MATTER.o_names_list[v+2] = "grhd_turbTauxz";
     MATTER.o_names_list[v+3] = "grhd_turbTauyy";
     MATTER.o_names_list[v+4] = "grhd_turbTauyz";
     MATTER.o_names_list[v+5] = "grhd_turbTauzz";
     v += 6;
   } 

  if (Getv("conservative_amr","yes"))  { 
    MATTER.USECAMR = 1; 
    MATTER.c_names_list = (char**) malloc( MATTER.NVq * sizeof(char*) );
    if(Getd("camr_time")   < 0 ) 
      MATTER.CAMRACTIVE = 1; 
    else 
      MATTER.CAMRACTIVE = 0; 
  } else
    MATTER.USECAMR = 0;

  if (MATTER.USECAMR) {
    MATTER.c_names_list[0] = "camr_D";
    MATTER.c_names_list[1] = "camr_Tau";
    MATTER.c_names_list[2] = "camr_Sx";
    MATTER.c_names_list[3] = "camr_Sy";
    MATTER.c_names_list[4] = "camr_Sz";
  }

  // set ptrs to main routines 
  MATTER.init    = grhd_init;
  MATTER.sources = grhd_sources;
  
  // atmosphere
  if      (Getv("grhd_use_atmosphere","no"))         MATTER.USEATM = 0;
  else if (Getv("grhd_use_atmosphere","ColdStatic")) MATTER.USEATM = GRHD_ATM_COLDSTATIC;
  else if (Getv("grhd_use_atmosphere","Vacuum"))     MATTER.USEATM = GRHD_ATM_VACUUM;
  else errorexit("unknown atmosphere");

  MATTER.USEATM_PRERHS = MATTER.USEATM_POSTRHS = MATTER.USEMASK = 0;
   
  // atmosphere and C2P
  if (MATTER.USEATM) {
    
    if (MATTER.USEATM == GRHD_ATM_COLDSTATIC)
    {
      if (Getv("grhd_use_atmosphere_prerhs","yes"))   MATTER.USEATM_PRERHS  = 1;
      if (Getv("grhd_use_atmosphere_postrhs","yes"))  MATTER.USEATM_POSTRHS = 1;
      if (Getv("grhd_use_atmosphere_mask","yes"))     MATTER.USEMASK        = 1;
      
      if (MATTER.USEATM_PRERHS || MATTER.USEATM_POSTRHS) {
        MATTER.set_atm_prerhs  = grhd_atm_set_prerhs;
        MATTER.set_atm_postrhs = grhd_atm_set_postrhs;
      } else errorexit(" you must enable at least one of: grhd_use_atmosphere_prerhs/grhd_use_atmosphere_postrhs");

      // reset c2p routine     
      if (EOS.COLD) {
        MATTER.c2p = grhd_c2p_rroot_ColdStaticAtm;
     } else { 
      	if  (Getv("grhd_C2P","p_root")) MATTER.c2p = grhd_c2p_proot_ColdStaticAtm_hybrid;	
        else if (Getv("grhd_C2P","p_root_pur")) MATTER.c2p = grhd_c2p_proot_ColdStaticAtm;
        else if (Getv("grhd_C2P","h_root")) MATTER.c2p = grhd_c2p_hroot_ColdStaticAtm;	
      	//else if (Getv("grhd_C2P","r_root")) MATTER.c2p = grhd_c2p_rroot;
        else errorexit(" this c2p method does not exist");
      }
    } 
    else if (MATTER.USEATM == GRHD_ATM_HYDROSTATIC)
    {
      errorexit(" to be implemented, see bam 11.04");
    } 
    else if (MATTER.USEATM == GRHD_ATM_VACUUM)
    {
      MATTER.sources = grhd_sources_Vac;
      if(Getv("grhd_use_turbulence","yes")){errorexit("can't use both \"atmosphere vacuum\" and \"turbulence\"");}
	
      if(Getv("grhd_use_atmosphere_mask","yes"))    MATTER.USEMASK        = 1;
      if(Getv("grhd_use_atmosphere_prerhs","yes"))  MATTER.USEATM_PRERHS  = 1;
      if(Getv("grhd_use_atmosphere_postrhs","yes"))
      {
        MATTER.USEATM_POSTRHS = 1;
        errorexit("not implemented yet");
      }
      
      if(MATTER.USEATM_PRERHS || MATTER.USEATM_POSTRHS)
      {
        MATTER.set_atm_prerhs  = grhd_Vac_set_prerhs;
        MATTER.set_atm_postrhs = NULL; // not implemented
      }
      else
        errorexit("you must enable at least one of: grhd_use_atmosphere_prerhs/grhd_use_atmosphere_postrhs");

      /* reset c2p routine */
      MATTER.c2p = grhd_c2p_proot_Vac; // so far only this one exists
      /* set flux routine */
      MATTER.hrsc_flx1d = grhd_flx1d_fvprb_Vac; // so far only this one exists
    } 
    else errorexit(" this atm method does not exist");
    
  } else {
    
    if      (Getv("grhd_C2P","p_root")) MATTER.c2p = grhd_c2p_proot;
    else if (Getv("grhd_C2P","h_root")) MATTER.c2p = grhd_c2p_hroot;
    //else if (Getv("grhd_C2P","r_root")) matter_c2p = grhd_c2p_rroot;
    else if (Getv("grhd_C2P","v_anal")) MATTER.c2p = grhd_c2p_vanal;
    else errorexit(" C2P does not exist");
    
     // checks
     if ( (Getv("grhd_C2P","v_anal")) && (!Getv("eos","ideal")) )
       errorexit("this C2P works only with Gamma-law EOS (eos=ideal)");


  }
  
  MATTER.set_mask = grhd_atm_set_mask;

  
  //turbulence
  if(Getv("grhd_use_turbulence","yes")){
	  MATTER.sources = grhd_sources_turb;
	  GRHD.turb_lmix = Getd("grhd_turb_lmix");
  }


  // excision
  MATTER.USEXCISION = 0;
  if (Getv("grhd_use_excision","yes")) {
    MATTER.USEXCISION = 1;    
    grhd_excision_init(level);
  } 

  // constant boundary (only for tests)
  if (Getv("grhd_use_atmosphere_postrhs","const_boundary")) {
    MATTER.set_atm_postrhs = grhd_atm_set_postrhs_const;
    MATTER.USEATM_POSTRHS = 1;
  }

  // hrsc stuff
  GRHD.HRSC_VMAX    = Getd("grhd_vmax");
  GRHD.HRSC_WlorMAX = Getd("grhd_Wlor_max");
  
  // flux1d 
  if      (Getv("hrsc_flux","HO"))  MATTER.hrsc_flx1d = grhd_flx1d_ho;
  
  else if (Getv("hrsc_flux","HO_LLF"))  {
    MATTER.hrsc_flx1d = grhd_flx1d_ho_highdens;
    MATTER.flx_LLF_HO_RHO = Getd("hrsc_flux_switch_rho");
    MATTER.flx_LLF_HO_ALPHA = Getd("hrsc_flux_switch_alpha");
  }

  else if (Getv("hrsc_flux","LLF") || Getv("hrsc_flux","HLL")) {
    if (Getv("grhd_recvel", "bx"))   MATTER.hrsc_flx1d = grhd_flx1d_fvprb;
    else                             MATTER.hrsc_flx1d = grhd_flx1d_fvpr;
    if (Getv("hrsc_rec","PPM"))      MATTER.hrsc_flx1d = grhd_flx1d_fvpr_ppm;
  }

  else if (Getv("hrsc_flux","HLLC")) {
    GRHD.HLLC_K = Getd("grhd_HLLC_K");
    GRHD.HLLC_mach_lim = Getd("grhd_HLLC_mach_lim");
    if (Getv("grhd_recvel", "bx"))   MATTER.hrsc_flx1d = grhd_flx1d_fvprb_HLLC;
    else                             MATTER.hrsc_flx1d = grhd_flx1d_fvpr_HLLC;
    if (Getv("hrsc_rec","PPM"))      MATTER.hrsc_flx1d = grhd_flx1d_fvpr_ppm_HLLC;
  }

  else errorexit("this hrsc_flux-function does not exist");

  /* for vaccum we need different fluxes */
  if(Getv("grhd_use_atmosphere","Vacuum"))
  {
    /* set flux routine */
    // we ignore hrsc_flux, because only few options exist for Vacuum
    printf("  Vacuum: evaluating hrsc_flux, ignoring grhd_recvel\n");
    MATTER.hrsc_flx1d = grhd_flx1d_fvprb_Vac;
    if(Getv("grhd_recepsl", "rhoepsl"))
    {
      MATTER.hrsc_flx1d = grhd_flx1d_fvpr_re_rWv_Vac;
      errorexit("  grhd_flx1d_fvpr_re_rWv_Vac gives bad results!!!");
    }
    if (Getv("hrsc_flux","HO_LLF"))
    {
      MATTER.hrsc_flx1d = grhd_flx1d_ho_vac_highdens;
      MATTER.flx_LLF_HO_RHO   = Getd("hrsc_flux_switch_minrho");
      MATTER.flx_LLF_HO_ALPHA = Getd("hrsc_flux_switch_alpha");
    }
  }
 

  if(Getv("grhd_use_turbulence","yes")){ 
  	if (Getv("grhd_recvel", "bx"))  MATTER.hrsc_flx1d = grhd_flx1d_fvprb_turb;
    else                            MATTER.hrsc_flx1d = grhd_flx1d_fvpr_turb;
  }

  return 0;
}


// initialized vars (called in matter_init(), matter.c) 
void grhd_init(tL *level)
{

  double sqrtdetgamma,detgamma;
  double vx,vy,vz,vsqr,vlowx,vlowy,vlowz; // tmp vars for v^i, v^2, W, v_i
  double sqr, W, W2hrho;
  double tmp1,tmp2,tmp3;
  
  // assume these are set by initial data routine
  double *gxx   = Ptr(level, "adm_gxx");
  double *gxy   = Ptr(level, "adm_gxy");
  double *gxz   = Ptr(level, "adm_gxz");
  double *gyy   = Ptr(level, "adm_gyy");
  double *gyz   = Ptr(level, "adm_gyz");
  double *gzz   = Ptr(level, "adm_gzz");
  
  double *rho   = Ptr(level, "grhd_rho");
  double *epsl  = Ptr(level, "grhd_epsl");
  double *p_vx  = Ptr(level, "grhd_vx");
  double *p_vy  = Ptr(level, "grhd_vy");
  double *p_vz  = Ptr(level, "grhd_vz");
  double *p     = Ptr(level, "grhd_p");
  	

  // set the following vars
  double *v2    = Ptr(level, "grhd_v2");
  double *det3g = Ptr(level, "adm_detg");
  
  double *D     = Ptr(level, "grhd_D");
  double *T     = Ptr(level, "grhd_Tau");
  double *Sx    = Ptr(level, "grhd_Sx");
  double *Sy    = Ptr(level, "grhd_Sy");
  double *Sz    = Ptr(level, "grhd_Sz");
  
  double *ADM_rho  = Ptr(level, "adm_rho");
  double *ADM_Sx   = Ptr(level, "adm_Sx");
  double *ADM_Sy   = Ptr(level, "adm_Sy");
  double *ADM_Sz   = Ptr(level, "adm_Sz");
  double *ADM_SSxx = Ptr(level, "adm_SSxx");
  double *ADM_SSxy = Ptr(level, "adm_SSxy");
  double *ADM_SSxz = Ptr(level, "adm_SSxz");
  double *ADM_SSyy = Ptr(level, "adm_SSyy");
  double *ADM_SSyz = Ptr(level, "adm_SSyz");
  double *ADM_SSzz = Ptr(level, "adm_SSzz");
  double *ADM_ST   = Ptr(level, "adm_ST");

  // (if required) set atmosphere levels

  if (MATTER.USEATM) grhd_atm_set_levels(level);
  /* set some things for Vac */
  if (MATTER.USEATM==GRHD_ATM_VACUUM) grhd_Vac_set_levels(level);

  //set index for other variables
  if (!(MATTER.NVo == 0)){
    if (Getv("grhd_use_turbulence","yes") && GRHD.IND_turbTau == 0) MATTER.INDX_VAR_o = Ind("grhd_turbTauxx");
  }

  // compute conservative and ADM vars    
  // (if required) atm is set here for the first time !	
  forallpoints_ijk(level) {
    
    // reset pressure ? 
    if (Getv("grhd_use_initialdata_p","yes")) {
      if (ijk==0) printf("  using pressure from initialdata\n");
    } else {
      EOS.comp("re","","","p","","",
               rho[ijk],epsl[ijk], &(p[ijk]));
    }
    
    // det 3-metric
    detgamma   = detg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk]); 
    det3g[ijk] = detgamma;
    sqrtdetgamma  = sqrt(detgamma);
   
    // (if required) set atm in primitives
    if (MATTER.USEATM && MATTER.USEATM!=GRHD_ATM_VACUUM)
      grhd_atm_set_prim_pt(&(p[ijk]), &(rho[ijk]), &(epsl[ijk]), &(p_vx[ijk]), &(p_vy[ijk]), &(p_vz[ijk]), &(v2[ijk]));
    
    // 3-velocity at ijk, lower 3-velocity v_i , and v2
    vx = p_vx[ijk];
    vy = p_vy[ijk];
    vz = p_vz[ijk];
    
    grhd_compute_v2_pt(vx,vy,vz, 
		       gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk],
		       &vlowx,&vlowy,&vlowz,&vsqr );

    W = 1.0/sqrt(1.0 - vsqr);
    if (vsqr>=1.) errorexit(" v2 > 1 in initial data");
    if (vsqr<0.)  errorexit(" v2 < 0 in initial data");
    
    v2[ijk] = vsqr;
    
    // compute conservative vars
    grhd_compute_q_pt(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk],detgamma, 
		      rho[ijk],epsl[ijk],p[ijk], vlowx,vlowy,vlowz, W,
		      &(D[ijk]), &(T[ijk]), &(Sx[ijk]), &(Sy[ijk]), &(Sz[ijk]));
    
    // compute ADM projections of stress energy tensor
    
    W2hrho = W*W*( rho[ijk] + rho[ijk]*epsl[ijk] + p[ijk] );

    // rho = W^2 h rho - p  
    ADM_rho[ijk]  = W2hrho - p[ijk];
    
    // j^i  = W^2 h rho v^i 
    ADM_Sx[ijk]   = W2hrho *vx;
    ADM_Sy[ijk]   = W2hrho *vy;
    ADM_Sz[ijk]   = W2hrho *vz;
    
    // S_{ij} = h rho W^2 v_i v_j + p g_{ij} 
    ADM_SSxx[ijk] = W2hrho *vlowx*vlowx + gxx[ijk]*p[ijk];
    ADM_SSxy[ijk] = W2hrho *vlowx*vlowy + gxy[ijk]*p[ijk];
    ADM_SSxz[ijk] = W2hrho *vlowx*vlowz + gxz[ijk]*p[ijk];
    ADM_SSyy[ijk] = W2hrho *vlowy*vlowy + gyy[ijk]*p[ijk];
    ADM_SSyz[ijk] = W2hrho *vlowy*vlowz + gyz[ijk]*p[ijk];
    ADM_SSzz[ijk] = W2hrho *vlowz*vlowz + gzz[ijk]*p[ijk];
    
    // S_i^j = trS = W^2 h rho v_i v^i + 3 p 
    ADM_ST[ijk]   = W2hrho * vsqr  + 3.*p[ijk];
    
  } endfor_ijk;

}


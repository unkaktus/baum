/* hrsc.c  */

#include "bam.h"
#include "matter.h"


// hrsc settings 1d flux and rec
int hrsc_startup(tL* level)
{
  
  printf("Initializing hrsc (matter): \n");
  printf("  provide hrsc varlists and fkt pointers to matter\n");

  
  // flux1d 
  if (Getv("hrsc_flux","HLL")) {
    MATTER.hrsc_riemsol1d = hrsc_numflx1d_hll;
    if (Getv("physics", "grmhd")||Getv("physics", "grmhdY")) MATTER.hrsc_riemsol1d = hrsc_numflx1d_hll_grmhd;
  }
  else if (Getv("hrsc_flux","LLF"))
    MATTER.hrsc_riemsol1d = hrsc_numflx1d_llf;
  
      
  // TVD limiter
  if      (Getv("hrsc_TVD_limiter","MM2")) 
    MATTER.hrsc_TVDlim = MM2;
  else if (Getv("hrsc_TVD_limiter","MC2")) 
    MATTER.hrsc_TVDlim = MC2;
  else errorexit(" unknown TVD_limiter");


  // rec for hydro vars  
  if       (Getv("hrsc_rec","CONST")) {
    MATTER.rec1d  = rec1d_godunov;
    MATTER.rec1dl = rec1d_p_godunov;
    MATTER.rec1dr = rec1d_m_godunov; 
  } else if (Getv("hrsc_rec","AVG")) {
    MATTER.rec1d  = rec1d_avg;
    MATTER.rec1dl = rec1d_p_avg;
    MATTER.rec1dr = rec1d_m_avg;  
  } else if (Getv("hrsc_rec","LINTVD")) {
    MATTER.rec1d  = rec1d_lintvd;
    MATTER.rec1dl = rec1d_p_lintvd;
    MATTER.rec1dr = rec1d_m_lintvd;
  } else if (Getv("hrsc_rec","PPM")) {
    MATTER.rec1d  = rec1d_lintvd;
    MATTER.rec1dl = rec1d_p_lintvd;
    MATTER.rec1dr = rec1d_m_lintvd;
  } else if (Getv("hrsc_rec","PPM4")) {
    MATTER.rec1d  = rec1d_ppm4;
    MATTER.rec1dl = rec1d_p_ppm4;
    MATTER.rec1dr = rec1d_m_ppm4;
  } else if (Getv("hrsc_rec","CENO3")) {
    MATTER.rec1d  = rec1d_ceno3;
    MATTER.rec1dl = rec1d_p_ceno3;
    MATTER.rec1dr = rec1d_m_ceno3;
  } else if (Getv("hrsc_rec","MUSCL")) {
    //MATTER.rec1d  = rec1d_muscl;
    //MATTER.rec1dl = rec1d_p_muscl;
    //MATTER.rec1dr = rec1d_m_muscl;
    errorexit("MUSCL reconstruction does not work yet, if you really wanna use it you have to fix it \n");
  }else if (Getv("hrsc_rec","WENO5")) {
    MATTER.rec1d  = rec1d_weno5;
    MATTER.rec1dl = rec1d_p_weno5;
    MATTER.rec1dr = rec1d_m_weno5;
  } else if (Getv("hrsc_rec","WENOZ")) {
    MATTER.rec1d  = rec1d_wenoz;
    MATTER.rec1dl = rec1d_p_wenoz;
    MATTER.rec1dr = rec1d_m_wenoz;
  } else if (Getv("hrsc_rec","MP5")) {
    //errorexit("MP5 routine implemented only for left interface");
    //rec1d  = mp5;
    MATTER.rec1dl = rec1d_p_mp5;
    MATTER.rec1dr = rec1d_m_mp5;
  }  else if (Getv("hrsc_rec","PPM6")) {
    MATTER.rec1d  = rec1d_ppm6;
    MATTER.rec1dl = rec1d_p_ppm6;
    MATTER.rec1dr = rec1d_m_ppm6;
  } else errorexit("unkown reconstruction");  

  // rec for hydro vars if we have to use a lower order
  MATTER.rec1dl_low = rec1d_p_godunov;
  MATTER.rec1dr_low = rec1d_m_godunov;
  if (Getv("hrsc_use_lower","yes")) {
    if        (Getv("hrsc_rec_low","CONST")) {
      MATTER.rec1dl_low = rec1d_p_godunov;
      MATTER.rec1dr_low = rec1d_m_godunov; 
    } else if (Getv("hrsc_rec_low","AVG")) {
      MATTER.rec1dl_low = rec1d_p_avg;
      MATTER.rec1dr_low = rec1d_m_avg;  
    } else if (Getv("hrsc_rec_low","LINTVD")) {
      MATTER.rec1dl_low = rec1d_p_lintvd;
      MATTER.rec1dr_low = rec1d_m_lintvd;
    } else if (Getv("hrsc_rec","PPM")) {
      MATTER.rec1dl_low = rec1d_p_lintvd;
      MATTER.rec1dr_low = rec1d_m_lintvd;
    } else if (Getv("hrsc_rec","PPM4")) {
      MATTER.rec1dl_low = rec1d_p_ppm4;
      MATTER.rec1dr_low = rec1d_m_ppm4;
    } else if (Getv("hrsc_rec_low","CENO3")) {
      MATTER.rec1dl_low = rec1d_p_ceno3;
      MATTER.rec1dr_low = rec1d_m_ceno3;
    } else if (Getv("hrsc_rec","MUSCL")) {
      MATTER.rec1dl_low = rec1d_p_muscl;
      MATTER.rec1dr_low = rec1d_m_muscl;
    } else if (Getv("hrsc_rec_low","WENO5")) {
      MATTER.rec1dl_low = rec1d_p_weno5;
      MATTER.rec1dr_low = rec1d_m_weno5;
    } else if (Getv("hrsc_rec_low","WENOZ")) {
      MATTER.rec1dl_low = rec1d_p_wenoz;
      MATTER.rec1dr_low = rec1d_m_wenoz;
    } else if (Getv("hrsc_rec_low","MP5")) {
      MATTER.rec1dl_low = rec1d_p_mp5;
      MATTER.rec1dr_low = rec1d_m_mp5;
    } else if (Getv("hrsc_rec","PPM6")) {
      MATTER.rec1dl_low = rec1d_p_ppm6;
      MATTER.rec1dr_low = rec1d_m_ppm6;
    } else errorexit("unkown reconstruction for lower order scheme"); 
  }

  // rec for metric vars
  if      (Getv("hrsc_rec_metric","CONST")) {
    MATTER.rec1dm  = rec1d_godunov;
    MATTER.rec1dml = rec1d_p_godunov;
    MATTER.rec1dmr = rec1d_m_godunov; 
  } else if (Getv("hrsc_rec_metric","AVG")) {
    MATTER.rec1dm  = rec1d_avg;
    MATTER.rec1dml = rec1d_p_avg;
    MATTER.rec1dmr = rec1d_m_avg; 
  } else if (Getv("hrsc_rec_metric","LINTVD")) {
    MATTER.rec1dm = rec1d_avg;
    MATTER.rec1dml = rec1d_p_lintvd;
    MATTER.rec1dmr = rec1d_m_lintvd;
  } else if (Getv("hrsc_rec_metric","CENO3")) {
    MATTER.rec1dm  = rec1d_ceno3;
    MATTER.rec1dml = rec1d_p_ceno3;
    MATTER.rec1dmr = rec1d_m_ceno3;
  } else if (Getv("hrsc_rec_metric","LAG4")) {
    MATTER.rec1dm  = rec1d_lag4;
    MATTER.rec1dml = rec1d_p_lag4;
    MATTER.rec1dmr = rec1d_m_lag4;
  } else if (Getv("hrsc_rec_metric","WENO5")) {
    MATTER.rec1dm  = rec1d_weno5;
    MATTER.rec1dml = rec1d_p_weno5;
    MATTER.rec1dmr = rec1d_m_weno5;
  } else if (Getv("hrsc_rec_metric","WENOZ")) {
    MATTER.rec1dm  = rec1d_wenoz;
    MATTER.rec1dml = rec1d_p_wenoz;
    MATTER.rec1dmr = rec1d_m_wenoz;
  } else if (Getv("hrsc_rec_metric","LAG6")) { 
    MATTER.rec1dm  = rec1d_lag6;
    MATTER.rec1dml = rec1d_p_lag6;
    MATTER.rec1dmr = rec1d_m_lag6; 
  } else errorexit(" unkown reconstruction");   


  // rec for tracer (Ye); it is LINTVD by default

  if(Getv("hrsc_rec_tracer","CONST")) {
    MATTER.rec1dsl = rec1d_p_godunov;
    MATTER.rec1dsr = rec1d_m_godunov;
  } else if (Getv("hrsc_rec_tracer","AVG")) {
    MATTER.rec1dsl = rec1d_p_avg;
    MATTER.rec1dsr = rec1d_m_avg;
  } else if (Getv("hrsc_rec_tracer","PPM")) {
    MATTER.rec1dsl = rec1d_p_lintvd;
    MATTER.rec1dsr = rec1d_m_lintvd;
  } else if (Getv("hrsc_rec_tracer","PPM4")) {
    MATTER.rec1dsl = rec1d_p_ppm4;
    MATTER.rec1dsr = rec1d_m_ppm4;
  } else if (Getv("hrsc_rec_tracer","LINTVD")) {
    MATTER.rec1dsl = rec1d_p_lintvd;
    MATTER.rec1dsr = rec1d_m_lintvd;
  } else if (Getv("hrsc_rec_tracer","CENO3")) {
    MATTER.rec1dsl = rec1d_p_ceno3;
    MATTER.rec1dsr = rec1d_m_ceno3;
  } else if (Getv("hrsc_rec_tracer","MUSCL")) {
    //MATTER.rec1dsl = rec1d_p_muscl;
    //MATTER.rec1dsr = rec1d_m_muscl;
    errorexit("MUSCL reconstruction does not work yet, if you really wanna use it you have to fix it \n");
  }else if (Getv("hrsc_rec_tracer","WENO5")) {
    MATTER.rec1dsl = rec1d_p_weno5;
    MATTER.rec1dsr = rec1d_m_weno5;
  } else if (Getv("hrsc_rec_tracer","WENOZ")) {
    MATTER.rec1dsl = rec1d_p_wenoz;
    MATTER.rec1dsr = rec1d_m_wenoz;
  } else if (Getv("hrsc_rec_tracer","MP5")) {
    MATTER.rec1dsl = rec1d_p_mp5;
    MATTER.rec1dsr = rec1d_m_mp5;
  } else if (Getv("hrsc_rec","PPM6")) {
    MATTER.rec1dl_low = rec1d_p_ppm6;
    MATTER.rec1dr_low = rec1d_m_ppm6;
  } else errorexit("unkown reconstruction");

  // ghost point for rec (1d)
  MATTER.HRSC_NGHOST = Geti("hrsc_nghosts");


  // check no ghosts
  if (MATTER.HRSC_NGHOST < 2) errorexit(" you need at least 2 hrsc_nghosts");
  if (Geti("bampi_nghosts") < 4) errorexit(" you need at least 4 bampi_nghosts"); 
  if ( Getv("hrsc_rec","CENO3") || 
       Getv("hrsc_rec","WENOZ") || 
       Getv("hrsc_rec","WENO5") || 
       Getv("hrsc_rec","MP5")   || 
       Getv("hrsc_rec_metric","CENO3") || 
       Getv("hrsc_rec_metric","WENOZ") || 
       Getv("hrsc_rec_metric","WENO5") || 
       Getv("hrsc_rec_metric","LAG6") ) {
    if (MATTER.HRSC_NGHOST < 3) errorexit(" you need at least 3 hrsc_nghosts for high-order rec");
    if (Geti("bampi_nghosts") < 6) errorexit(" you need at least 4 bampi_nghosts for high-order rec ... use 6!");
  }

}










/* eos.c */
/* mth 06/12 */

#include "bam.h"
#include "eos.h"

#define PR 0
#define debug 0

/* infos */
tEOS EOS;







/* settings called within the bam_eos.c file
   this routine is used mainly to interface the options from parfile 
   and the grhd specific routine with the general setup of the matter* routines 
   and it is stored inside a global struct!
*/
int eos_startup(tL *level) 
{
  printf("Initializing eos:\n");
  printf("  provide eos fkt pointers to matter\n");
 
  
  EOS.use1D       = eos_dummy;
  EOS.use1DH      = eos_dummy;
  EOS.use2D       = eos_dummy;
  EOS.use3D       = eos_dummy;
  EOS.use3DT      = eos_dummy;
  EOS.use3DT_D    = eos_dummy;

  // EoS
  if        (Getv("eos","dust"))     {
    EOS.type        = DUST;
    EOS.use1D       = eos_dust;
    EOS.use2D       = eos_dust;
    EOS.COLD        = 1;
  } else if (Getv("eos","poly"))     {
    EOS.type        = POLY;
    EOS.use1D       = eos_poly;
    EOS.use1DH      = eos_poly_H;
    EOS.use2D       = eos_poly;
    EOS.COLD        = 1;
    EOS.Y           = Getd("eos_Y"); //HG: test
  } else if (Getv("eos","pwp"))      {
    EOS.type        = PWP;
    EOS.use1D       = eos_pwp;
    EOS.use1Dp	    = eos_pwp_p;
    EOS.use1DH      = eos_pwp_H;
    EOS.use2D       = eos_pwp;
    EOS.COLD        = 1;
    eos_load_pwp();
  } else if (Getv("eos","tab1d"))    {
    EOS.type        = TAB1D;
    EOS.use1D       = eos_cold_tab1d;
    EOS.use2D       = eos_cold_tab1d;
    EOS.COLD        = 1;
    eos_load_tab1d();
  } else if (Getv("eos","ideal"))    {
    EOS.type        = IDEAL;
    EOS.use1D       = eos_poly;
    EOS.use1DH      = eos_poly_H;
    EOS.use2D       = eos_ideal;
    EOS.COLD        = 0;
    EOS.Y	    = Getd("eos_Y");
    EOS.T 	    = Getd("eos_T");
  } else if (Getv("eos","idealWT"))    {
    EOS.type        = IDEAL;
    EOS.use1D       = eos_poly;   /* <--WT: check if this produces NANs */
    EOS.use1DH      = eos_poly_H; /* <--WT: check if this produces NANs */
    EOS.use2D       = eos_ideal_WT;
    EOS.COLD        = 0; 
  } else if (Getv("eos","pwphot"))   {
    EOS.type        = PWPHOT;
    EOS.use1D       = eos_pwp;
    EOS.use1DH      = eos_pwp_H;
    EOS.use2D       = eos_pwpHot;
    EOS.COLD        = 0;
    EOS.Y           = Getd("eos_Y");
    EOS.T           = Getd("eos_T");
    eos_load_pwp();
    //errorexit("please test me"); // 08/20/2012
  } else if (Getv("eos","pwphotWT")) {
    EOS.type        = PWPHOT;
    EOS.use1D       = eos_pwp;   /* <--WT: check if this produces NANs */
    EOS.use1DH      = eos_pwp_H; /* <--WT: check if this produces NANs */
    EOS.use2D       = eos_pwpHot_WT;
    EOS.COLD        = 0;
    eos_load_pwp(); /* <--WT: I don't expect this one to cause NANs */
  } else if (Getv("eos","tab1dhot")) {
    EOS.type        = TAB1DHOT;
    errorexit("please test me");
  } else if (Getv("eos","tab2d")) {
    EOS.type        = TAB2D;
    errorexit("please test me");
  } else if (Getv("eos","tab3d")) {
    EOS.type        = TAB3D;
    EOS.use1D       = eos_tab3d_shen_T0_Y;
    EOS.use2D       = eos_tab3d_shen_Y;
    EOS.use3D       = eos_tab3d_shen;
    EOS.use3DT      = eos_tab3d_shen_T;
    EOS.COLD        = 0;
    EOS.Y           = Getd("eos_Y");
    eos_load_tab3d();
  } else if (Getv("eos","nuclear")) {
    double dummy;
    EOS.type        = TAB3D;
    if(Getv("eos_interp", "cspline")){
    	EOS.use1D       = peps_1d;
	EOS.betaEoS     = beta_1d;
    }
    else if (Getv("eos_interp","linear")){
	EOS.use1D       = peps_1d_lin;
	EOS.betaEoS	= beta_1d_lin;
    }
    else if (Getv("eos_interp","steffen")){
	EOS.use1D 	= peps_1d_steffen;
	EOS.betaEoS	= beta_1d_steffen;
    }
    EOS.use2D       = eos_tab3d_Y;
    EOS.use3D       = eos_tab3d;
    EOS.use3DT      = eos_tab3d_T;
    if(Getv("hrsc_flux","HO_LLF")) EOS.use3DT_D = eos_tab3d_T_D;
    EOS.COLD        = 0;
    EOS.Y           = Getd("eos_Y");
    EOS.T 	    = Getd("eos_T");
    EOS.micro	    = eos_tab3d_micro; //HG
    EOS.mu	    = eos_tab3d_mu; //FS
    EOS.beta_eq	= eos_tab3d_eq; //FS
    eos_load_complete();
    eos_validity(&dummy, &dummy, &(EOS.rhomin), &(EOS.rhomax), &(EOS.Ymin), &(EOS.Ymax), &dummy,
		 &dummy, &(EOS.Tmin), &(EOS.Tmax), &(EOS.mb), &(EOS.rhomin1D), &(EOS.mb1D));
    EOS.extend      = extend_validity;
    if(Getv("grhd_use_nls", "yes")) {
		EOS.nls_free = nls_free;
		EOS.nls_diff = nls_diff;
		EOS.zeta = nls_zeta;
    }
    if(Getv("physics", "grrhd_m1")) EOS.M1 = eos_m1_rates;
 } else if (Getv("eos","tab1d_cold")) {
    double dummy;
    EOS.type        = TAB1D;
    EOS.COLD	    = 0;
    EOS.Y	    = Getd("eos_Y");
    if(Getv("eos_interp", "cspline")){
    	EOS.use1D       = eos_tab1d;
    	EOS.use2D       = eos_tab1d;
    	EOS.betaEoS     = eos_tab1d_beta;
    	EOS.extend      = eos_tab1d_extend;
	EOS.extrap	= eos_tab1d_extrapolate;
	EOS.use1Dp	= eos_tab1d_p;
    } else if(Getv("eos_interp", "steffen")){
    	EOS.use1D       = eos_tab1d_stf;
    	EOS.use1Dp	= eos_tab1d_stf_p;
    	EOS.use2D       = eos_tab1d_stf;
    	EOS.betaEoS     = eos_tab1d_beta_stf;
    	EOS.extend	= eos_tab1d_extend_stf; 
    } else if(Getv("eos_interp", "linear")){
	EOS.use1D	= eos_tab1d_lin;
	EOS.use1Dp	= eos_tab1d_lin_p;
	EOS.use2D	= eos_tab1d_lin;
	EOS.betaEoS	= eos_tab1d_beta_lin;
	EOS.extend 	= eos_tab1d_extend_lin;
    }
    eos_load_tab1d_cold();
    eos_tab1d_validity(&dummy, &dummy, &(EOS.rhomin), &(EOS.rhomax), 
			       &(EOS.Ymin), &(EOS.Ymax), &(EOS.epsmin), &(EOS.epsmax), &(EOS.mb));
 } else if (Getv("eos","tab1d_hot")){
    double dummy;
    EOS.type        = TAB1DHOT;
    EOS.use1D       = eos_tab1d_hot;
    EOS.use2D       = eos_tab1d_hot;
    EOS.COLD        = 0;
    EOS.Y           = Getd("eos_Y");
    EOS.betaEoS     = eos_tab1d_hot_beta;
    EOS.micro       = eos_tab1d_hot_micro;
    EOS.extend      = eos_tab1d_hot_extend;
    eos_load_tab1d_hot();
    eos_tab1d_hot_validity(&dummy, &dummy, &(EOS.rhomin), &(EOS.rhomax),
			   &(EOS.Ymin), &(EOS.Ymax), &(EOS.epsmin), &(EOS.epsmax), &(EOS.mb));
  }

  else 
    errorexit(" EoS does not exist");
 
  /* set additional parameters */
  EOS.K           = Getd("eos_K");
  EOS.GAMMA       = Getd("eos_Gamma");
  EOS.GAMMAMO     = EOS.GAMMA-1.;

  /* set functions insie the struct */
  EOS.comp        = eos_comp;
  
  return 0;
}



/* EoS Wrapper ... go in with a general interface and depending on the
   letters choose the correct function pointer and give back values 
   which are needed 
=> this is slow, yes, but in compare to the 3d interpolation this should 
   not be the bottleneck
*/
int eos_comp(char*in,char*din,char*ddin, char*out,char*dout,char*ddout, ...)
{
  // all values which are possible
  double p,epsl,rho,Y,T,h,H, cs2, e;
  double dpdrho,dpdepsl, dhdrho;
  double dhdeps;
  double mun, mue, mup, A, Z, nh;
  double dedp, drhodp;
  
  int i;
  int Nin,Ndin,Nddin,Nout,Ndout,Nddout;
  Nin=Ndin=Nddin=Nout=Ndout=Nddout = 0;
  double vars_in[10],vars_din[10],vars_ddin[10];
  double *vars_out[10],*vars_dout[10],*vars_ddout[10];

  /* look for in and outgoing variables */
  va_list args;
  va_start(args, ddout);
  while (in[Nin]!='\0')
    vars_in[Nin++] = va_arg(args,double);
  while (din[Ndin]!='\0')
    vars_din[Ndin++] = va_arg(args,double);
  while (ddin[Nddin]!='\0')
    vars_ddin[Nddin++] = va_arg(args,double);
  while (out[Nout]!='\0')
    vars_out[Nout++] = va_arg(args,double*);
  while (dout[Ndout]!='\0')
    vars_dout[Ndout++] = va_arg(args,double*);
  while (ddout[Nddout]!='\0')
    vars_ddout[Nddout++] = va_arg(args,double*);
  va_end(args);
  
  
  if (debug) {
    printf("------------------------------------\n");
    if (Nin==0 || Nout==0)
      errorexit("eos call is wrong");
    printf(" (%s | %s | %s)   ||   ( %s | %s | %s)\n",in,din,ddin,out,dout,ddout);
    printf(" N = %d %d %d  |  %d %d %d\n", Nin,Ndin,Nddin,Nout,Ndout,Nddout);
    for (i=0; i<Nin; i++)
      printf(" in[%d]    = %e \n",i, vars_in[i]);
    for (i=0; i<Ndin; i++)
      printf(" din[%d]   = %e \n",i, vars_din[i]);
    for (i=0; i<Nddin; i++)
      printf(" ddin[%d]  = %e \n",i, vars_ddin[i]);
  }
  
  int RETURN;
  
  /* *********************************************************************** */
  /* p(rho, ...) */
  if (in[0]=='r' && out[0]=='p') {
    
    /* --------------------------------------------------------------------- */
    /* p(rho) */
    if (Nin==1) {
      
      rho = vars_in[0];
      
      /* output vars */
      if (Nout==3 && out[1]=='e' && out[2]=='c') {
        if (debug) printf("EOS:  rho  --p(rho)-->  p,epsl,cs2");
        epsl = -10.;
        RETURN = EOS.use1D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
        *(vars_out[0]) = p;
        *(vars_out[1]) = epsl;
        *(vars_out[2]) = cs2;
      } else if (Nout==2 && out[1]=='e') {
        if (debug) printf("EOS:  rho  --p(rho)-->  p,epsl");
        epsl = -10.;
        RETURN = EOS.use1D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
        *(vars_out[0]) = p;
        *(vars_out[1]) = epsl;
      } else if (Nout==2 && out[1]=='c') {
        if (debug) printf("EOS:  rho  --p(rho)-->  p,cs2");
        epsl = (epsl<0.)?0.:epsl;
        RETURN = EOS.use1D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
        *(vars_out[0]) = p;
        *(vars_out[1]) = cs2;
      } else if (Nout==1) {
        if (debug) printf("EOS:  rho  --p(rho)-->  p");
        epsl = (epsl<0.)?0.:epsl;
        RETURN = EOS.use1D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
        *(vars_out[0]) = p;
      } else {
        errorexit("EoS");
      }
        
      /* output derivatives ? */
      if (Ndout==2 && dout[0]=='r' && dout[1]=='e') {
        if (debug) printf(",dpdrho,dpdepsl\n");
        *(vars_dout[0]) = dpdrho;
        *(vars_dout[1]) = dpdepsl;
      } else if (Ndout==1 && dout[0]=='r') {
        if (debug) printf(",dpdrho\n");
        *(vars_dout[0]) = dpdrho;
      } else if (Ndout==0) {
        if (debug) printf("\n");
      } else {
        errorexit("EoS");
      }
      
    /* --------------------------------------------------------------------- */
    /* p(rho,T) */

    } else if(Nin == 2 && (in[0] == 'r' && in[1] == 'T') ) {
        rho = vars_in[0];
        T   = vars_in[1];

        RETURN = EOS.use2D(&p,&cs2,&dpdrho,&dpdepsl,&rho,&epsl);

        if(Nout == 1 && out[0]=='p') *(vars_out)[0] = p;

        if(Nout==2 && out[0]=='p' && out[1]=='c'){
                *(vars_out)[0] = p;
                *(vars_out)[1] = cs2;
        }


        if(Nout==2 && out[0]=='p' && out[1]=='e'){
          *(vars_out)[0] = p;
          *(vars_out)[1] = epsl;
	  if(Ndout == 1 && dout[0] == 'r') *(vars_dout)[0] = dpdrho;
        }
   } else if (Nin==2 && in[1]=='e') {
      
      rho  = vars_in[0];
      epsl = vars_in[1];
      
      if (Nout==2 && out[1]=='c' && Ndout==2 && dout[0]=='r' && dout[1]=='e') {
        if (debug) printf("EOS:  rho,epsl  --p(rho,epsl)-->  p,cs2,dpdrho,dpdepsl\n");
        RETURN = EOS.use2D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
        *(vars_out[0])  = p;
        *(vars_out[1])  = cs2;
        *(vars_dout[0]) = dpdrho;
        *(vars_dout[1]) = dpdepsl;
      } else if (Nout==2 && out[1]=='c') {
        if (debug) printf("EOS:  rho,epsl  --p(rho)-->  p,cs2\n");
        RETURN = EOS.use2D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
        *(vars_out[0]) = p;
        *(vars_out[1]) = cs2;
      } else if (Nout==1) {
        if (debug) printf("EOS:  rho,epsl  --p(rho)-->  p\n");
        epsl = (epsl<0.)?0.:epsl;
        RETURN = EOS.use2D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
        *(vars_out[0]) = p;
      } else {
        errorexit("EoS");
      }
      
    /* --------------------------------------------------------------------- */
    /* p(rho,T,Y) */
    } else if (Nin==3 && in[1]=='T' && in[2]=='Y') {
      
      rho  = vars_in[0];
      T    = vars_in[1];
      Y    = vars_in[2];
      
      if (Nout==1){
	if (debug) printf("EOS:  rho,T,Y  --p(rho,T,Y)-->  p\n");
	RETURN = EOS.use3DT(&rho, &epsl, &Y, &p, &T, &cs2, &dpdrho, &dpdepsl);
	*(vars_out[0]) = p;
      }

      else if (Nout==2 && out[1]=='e') {
        if (debug) printf("EOS:  rho,T,Y  --p(rho,T,Y)-->  p,epsl\n");
        RETURN = EOS.use3DT( &rho,&epsl,&Y, &p,&T,&cs2, &dpdrho,&dpdepsl );
        *(vars_out[0]) = p;
        *(vars_out[1]) = epsl;
      }

      else if(Nout==2 && out[1]=='c'){
	if (debug) printf("EOS: rho,T,Y --p(rho,T,Y)--> p, cs2\n");
	RETURN = EOS.use3DT( &rho,&epsl,&Y, &p, &T, &cs2, &dpdrho,&dpdepsl );
	*(vars_out[0]) = p;
	*(vars_out[1]) = cs2;

      } else {
        errorexit("EoS");
      }
    /* --------------------------------------------------------------------- */
    /* p(rho,epsl,Y) */
    } else if (Nin==3 && in[1]=='e' && in[2]=='Y') {
      
      rho  = vars_in[0];
      epsl = vars_in[1];
      Y    = vars_in[2];
      
      if (Nout==3 && out[1]=='c' && out[2]=='T' && Ndout==2 && dout[0]=='r' && dout[1]=='e') {
        if (debug) printf("EOS:  rho,epsl,Y  --p(rho,T,Y)-->  p,cs2,T,dpdrho,dpdepsl\n");
        RETURN = EOS.use3D( &rho,&epsl,&Y, &p,&T,&cs2, &dpdrho,&dpdepsl);
        *(vars_out[0])  = p;
        *(vars_out[1])  = cs2;
        *(vars_out[2])  = T;
        *(vars_dout[0]) = dpdrho;
        *(vars_dout[1]) = dpdepsl;
      } else if (Nout==2 && out[1]=='c' && Ndout==2 && dout[0]=='r' && dout[1]=='e') {
        if (debug) printf("EOS:  rho,epsl,Y  --p(rho,T,Y)-->  p,cs2,dpdrho,dpdepsl\n");
        RETURN = EOS.use3D( &rho,&epsl,&Y, &p,&T,&cs2, &dpdrho,&dpdepsl);
        *(vars_out[0])  = p;
        *(vars_out[1])  = cs2;
        *(vars_dout[0]) = dpdrho;
        *(vars_dout[1]) = dpdepsl;
      } else if (Nout==2 && out[1]=='c' && Ndout==0 ) {
        if (debug) printf("EOS:  rho,epsl,Y  --p(rho,T,Y)-->  p,cs2\n");
        RETURN = EOS.use3D( &rho,&epsl,&Y, &p,&T,&cs2, &dpdrho,&dpdepsl);
        *(vars_out[0])  = p;
        *(vars_out[1])  = cs2;
      } else if (Nout==1 && Ndout==0 ) {
        if (debug) printf("EOS:  rho,epsl,Y  --p(rho,T,Y)-->  p\n");
        RETURN = EOS.use3D( &rho,&epsl,&Y, &p,&T,&cs2, &dpdrho,&dpdepsl);
        *(vars_out[0])  = p;
      } else {
        errorexit("EoS");
      }
      
    } else {
      errorexit("EoS");
    }
    
  /* *********************************************************************** */
  /* epsl(rho, ...) */
  } else if (in[0]=='r' && out[0]=='e') {
    
    rho = vars_in[0];
    
    if (Nout==1) {
      if (debug) printf("EOS:  rho  --epsl(rho)-->  epsl\n");
      epsl = -10.;
      RETURN = EOS.use1D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
      *(vars_out[0]) = epsl;
    } else {
      errorexit("EoS");
    }
    
  /* *********************************************************************** */
  /* h(rho, ...) */
  } else if (in[0]=='r' && out[0]=='h') {
    
    rho = vars_in[0];
    
    if (Nout==2 && out[1]=='e' && Ndout==1 && dout[0]=='r') {
      if (debug) printf("EOS:   rho  --h(rho)-->  h,epsl,dhdrho\n");
      epsl = -10.;
      RETURN = EOS.use1D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
      h       = 1.0 + epsl + p/rho;
      dhdrho  = dpdrho/rho;
      *(vars_out[0])  = h;
      *(vars_out[1])  = epsl;
      *(vars_dout[0]) = dhdrho;
    } else if (Nout==1 && Ndout==1 && dout[0]=='r') {
      if (debug) printf("EOS:   rho  --h(rho)-->  h,dhdrho\n");
      RETURN = EOS.use1D( &p, &cs2, &dpdrho,&dpdepsl, &rho,&epsl);
      h       = 1.0 + epsl + p/rho;
      dhdrho  = dpdrho/rho;
      *(vars_out[0])  = h;
      *(vars_dout[0]) = dhdrho;
    } else {
      printf(" (%s | %s | %s)   ||   ( %s | %s | %s)\n",in,din,ddin,out,dout,ddout);
      errorexit("EoS");
    }
    
  /* *********************************************************************** */
  /* p(H, ...) */
  } else if (in[0]=='H' && out[0]=='p') {
    
    H = vars_in[0];
    
    //printf(" check %c %c %c\n",out[0],out[1],out[2]);

    if (Nout==3 && out[1] == 'r' && out[2] == 'e') {
      if (debug) printf("EOS:   H  --p(H)-->  p,rho,epsl\n");
      H = (H<0.)?0.:H;
      RETURN = EOS.use1DH(&H, &rho,&p,&epsl);
      *(vars_out[0]) = p;
      *(vars_out[1]) = rho;
      *(vars_out[2]) = epsl;
    } else {
      errorexit("EoS");
    }
  /* *********************************************************************** */
  /*HG: Y(r) for beta-equilibrated EoS */
  } else if (in[0]=='r' && out[0]=='Y' && out[1]=='T') {

    rho = vars_in[0];

    if(Nout==2) {
      if (debug) printf("Beta EoS: rho --Y(rho)--> Y");
      RETURN = EOS.betaEoS(&rho, &Y, &T);
      *(vars_out[0]) = Y;
      *(vars_out[1]) = T;

    }
   /*HG: Microphysics variables; the general call is EOS.comp("rTY","","","npeAZh","","") */
  } else if(in[0]=='r' && out[0]=='n') {

      rho = vars_in[0];

      if(Nin==3 && in[1]=='T' && in[2]=='Y' && out[1]=='p' && out[2]=='e'
	 && out[3]=='A' && out[4]=='Z' && out[5]=='h')
      {
         T = vars_in[1];
         Y = vars_in[2];

         RETURN = EOS.micro(&rho, &T, &Y, &mun, &mup, &mue, &A, &Z, &nh);
         *(vars_out[0]) = mun;
         *(vars_out[1]) = mup;
         *(vars_out[2]) = mue;
         *(vars_out[3]) = A;
         *(vars_out[4]) = Z;
         *(vars_out[5]) = nh;
      }

  /*e(p), rho(p), cs2(p) and derivs */
  

  } else if(in[0] =='p') {

      p = vars_in[0];

      RETURN = EOS.use1Dp(&p,&cs2,&drhodp,&dedp,&rho,&e);

      if(Nout==2 && out[0]=='r' && out[1]=='e'){
		*(vars_out)[0] = rho;
		*(vars_out)[1] = e;
      }

      
      if(Ndout==2 && dout[0]=='r' && dout[1]=='e'){
		*(vars_dout)[0] = drhodp;
		*(vars_dout)[1] = dedp;
      }
      else if(Ndout==1 && dout[0]=='r') *(vars_dout)[0] = drhodp;
      else if(Ndout==1 && dout[0]=='e') *(vars_dout)[0] = dedp;


      //} else {
        // errorexit("EoS");
      //}

  /* *********************************************************************** */
  } else if(in[0]=='r' && in[1]=='T' && out[0]=='Y') {

                  RETURN = EOS.betaEoS(&rho, &Y, &T);
                 *(vars_out)[0] = Y;
  } else {
    errorexit("EoS i/o option not recognized");
  }
  
  
  
  
  
  if (debug) {
    
    for (i=0; i<Nout; i++)
      printf(" out[%d]    = %e \n",i, *(vars_out[i]));
    for (i=0; i<Ndout; i++)
      printf(" dout[%d]   = %e \n",i, *(vars_dout[i]));
    for (i=0; i<Nddout; i++)
      printf(" ddout[%d]  = %e \n",i, *(vars_ddout[i]));
    printf("------------------------------------\n");
  }
  
  
  
  
  
  return RETURN;
}




















/* grhd_c2p_v2.c 
 mth, sbernuz 03/12 */

#include "bam.h"
#include "grhd.h"


#define PR 0
#define DEBUG 0
#define PRERROR 0

#define BIG 1e+15

enum {
    COLD,
    ATM,
    OK,
    FAIL,
    NR
};


// c2p using full EoS
int try_hot_part(double conD, double conS2, double conT,
                 double *rho, double *epsl, double *pres, 
                 double init, double errormax,double countmax)
{
    int count;
    
    double temp,tempSQRTarg,tempSQRT, Wlor;
    double pold,peos,pmin;
    double f,df,drhodp,depsldp,kappa,chi,cs2,v2;
    
    double problem,error,lostdigits;
    
    // set minimal p
    pmin = sqrt(conS2)-conD-conT;
    pmin = DMAX(pmin,GRHD.ATM_PATM);
    
    // guess pold from previous step 
    pold = DMAX(init,pmin);
    if (!finite(pold)) {
        if (PR) printf("starting guess for p is not finite");
        return FAIL;
    }
    
    // compute values with first guess
    temp        = conT+pold+conD;
    tempSQRTarg = temp*temp-conS2;
    tempSQRT    = sqrt(tempSQRTarg);
    Wlor        = temp/tempSQRT;
    *rho        = conD/Wlor;
    *epsl       = (tempSQRT - Wlor*pold)/conD - 1.;	
    
     
    // check if cons vars have physical values
    if (tempSQRTarg<=0.0) 
        return COLD;
    
    
    // start Newton Raphson and find pres
    for (count=0,error=2*errormax; count<=countmax; count++) {
        
        // EOS call
        problem = EOS.comp("re","","","pc","re","",
                           *rho,*epsl, &peos,&cs2, &chi,&kappa);
        
        f       = pold - peos; 
        drhodp  = conD*conS2/(tempSQRT *temp*temp);
        depsldp = pold*conS2/(conD*tempSQRTarg*tempSQRT);
        df      = 1. - chi*drhodp - kappa*depsldp;
        
        *pres   = DMAX( pmin, pold-f/df );
        error   = fabs( 1.-pold/(*pres) );
        pold    = *pres;
        
        temp        = conT+pold+conD;
        tempSQRTarg = temp*temp-conS2;
        tempSQRT    = sqrt(tempSQRTarg);
        Wlor        = temp/tempSQRT;
        *rho        = conD/Wlor;
        *epsl       = (tempSQRT - Wlor*pold)/conD - 1.;	
        
        // go out if something
        if (*rho<GRHD.ATM_FATM*GRHD.ATM_RHOATM) return ATM;
        if (problem) return ATM;
        if (df==0.) return NR;
        if (error<errormax) break; 
        
        if ((DEBUG) && (count>2)) {
            printf("    %d   %+2.16e   %+2.16e %+2.16e %2.16e %+2.16e %2.16e\n", count,error,df,f,*pres,*epsl,*rho);
            printf("         => %+2.16e %+2.16e\n", chi*drhodp,kappa*depsldp);
            printf("         => %+2.16e %+2.16e\n",  pold,f/df);
            printf("         => %+2.16e %+2.16e %+2.16e\n", conT,pold,conD);
            printf("         => %+2.16e %+2.16e\n",  temp*temp,conS2);
            printf("         => %+2.16e %+2.16e\n",  tempSQRT,Wlor*pold);
            printf("         => %+2.16e %+2.16e\n",  (tempSQRT - Wlor*pold)/conD,1.);
        }
        
        
    }
    
    
    if (count>=countmax) {
        lostdigits = 0;
        lostdigits = DMAX(lostdigits , fabs(log10(conD)-log10(conT)));
        lostdigits = DMAX(lostdigits , fabs(log10(conD)-log10(*pres)));
        lostdigits = DMAX(lostdigits , fabs(log10(conT)-log10(*pres)));
        //lostdigits = DMAX(lostdigits , fabs(log10(temp*temp)-log10(conS2)));
        lostdigits = DMAX(lostdigits , fabs(log10(tempSQRT)-log10(Wlor*pold)));
        //lostdigits = DMAX(lostdigits , fabs(log10((tempSQRT - Wlor*pold)/conD)-log10(1.)));
        //lostdigits = DMAX(lostdigits , fabs(log10(chi*drhodp)-log10(kappa*depsldp)));
        //printf("  ==> %2.2f    (%2.2f)\n", lostdigits, log10(error)-log10(1e-16));
        
        if ( lostdigits  <  log10(error)-log10(1e-15) || error > 1e-4 ) {
            if (PR) {
                printf("too many Newton-Raphson steps (hot)\n");
                printf("  %d  %e %e    %e %e %e\n",count,error, errormax, f/df, f,df);
                printf("  %+e %+e %+e (values  rho,epsl,p)\n",*rho,*epsl,*pres);
                printf("  %+e %+e %+e (atm)\n",GRHD.ATM_RHOATM,GRHD.ATM_EPSLATM,GRHD.ATM_PATM);
                printf("  %+e %+e %+e (D,S2,T)\n",conD,conS2,conT);
            }
            return NR;
        }
        return NR;
    }
    

    if (!finite(*pres) || !finite(*rho) || !finite(*epsl)) return COLD;
    if (*epsl<0.) return COLD;
    
    // test velocity
    v2 = conS2/((conT+conD+*pres)*(conT+conD+*pres));
    if (v2>GRHD.HRSC_VMAX) 
        return ATM;
    
    return OK;
}


// c2p using cold part soley
int try_cold_part(double conD, double conS2, double *conT,
                  double *rho, double *epsl, double *pres, 
                  double init, double errormax,double countmax)
{
    int count;
    
    double temp,tempSQRTarg,tempSQRT, Wlor;
    double rhomin,rhold;
    double g,dg,h,dhdrho,chi, v2;
    double problem,error;
  
    // check if cons vars have physical values
    if (conD <= 0.) 
        return ATM;
    
    rhomin = GRHD.ATM_RHOATM;
    rhold  = DMAX(init,rhomin);
    
    if (!finite(rhold)) {
        if (PR) printf("starting guess for rho is not finite");
        return FAIL;
    }

    // start Newton Raphson and find rho
    for (count=0,error=2*errormax; count<=countmax; count++) {
        
      // EOS call
      EOS.comp("r","","","pe","r","", rhold, pres,epsl, &chi);	
      h       = 1.0 + *epsl + *pres/rhold;		
      dhdrho  = chi/rhold;				
      
      /*
      // test 04/02/2013
      *epsl   = EOS.K*pow(rhold,EOS.GAMMAMO)/(EOS.GAMMAMO);
      *pres   = EOS.GAMMAMO*(*epsl)*(rhold);
      chi     = (EOS.K*EOS.GAMMA)*pow(rhold,EOS.GAMMAMO);
      h       = 1.0 + *epsl + *pres/rhold;
      dhdrho  = chi/rhold;	
      */
      
      // The following is a BUG: pressure is not set !!!
      /* problem = EOS.comp("r","","","h","r","", rhold, &h, &dhdrho); // do not reset epsl */
      /* problem = EOS.comp("r","","","he","r","", rhold, &h,epsl, &dhdrho); // reset epsl */
      

      if (conD*h==0.)       tempSQRTarg = 1.0 + conS2*BIG;
        else                tempSQRTarg = 1.0 + conS2/((conD*h)*(conD*h));
        Wlor    = sqrt(tempSQRTarg);
        
        g       = Wlor*rhold - conD; 
        dg      = Wlor - (rhold*conS2*dhdrho)/(Wlor*conD*conD*h*h*h);
        
        if ((!finite(dg)) || (dg==0.0)) {
	  if (PR) printf("dg==0\n");
	  return NR;
        }
        
        *rho    = DMAX( rhomin, rhold-g/dg );
        error   = fabs( 1.-rhold/(*rho) );
        rhold   = *rho;
        
        if (error<errormax) break;
    }
    
    if (count>=countmax) {
        if (PR) {
            printf("too many Newton-Raphson steps (cold)\n");
            printf("  %d  %e %e    %e %e %e\n",count,error, errormax, g/dg, g,dg);
            printf("  %+e %+e %+e (values  rho,epsl,p)\n",*rho,*epsl,*pres);
            printf("  %+e %+e %+e (atm)\n",GRHD.ATM_RHOATM,GRHD.ATM_EPSLATM,GRHD.ATM_PATM);
            printf("  %+e %+e %+e (D,S2,T)\n",conD,conS2,*conT);
        }
        return NR;
    }


    // EOS call
    EOS.comp("r","","","pe","r","", rhold, pres,epsl, &chi);	
    h       = 1.0 + *epsl + *pres/rhold;		
    dhdrho  = chi/rhold;				

    /*
    // test 04/02/2013
    *epsl   = EOS.K*pow(rhold,EOS.GAMMAMO)/(EOS.GAMMAMO);
    *pres   = EOS.GAMMAMO*(*epsl)*(rhold);
    chi     = (EOS.K*EOS.GAMMA)*pow(rhold,EOS.GAMMAMO);
    h       = 1.0 + *epsl + *pres/rhold;
    dhdrho  = chi/rhold;
    */

    // The following is a BUG: pressure is not set !!!
    /* problem = EOS.comp("r","","","he","r","", *rho, &h,epsl, &dhdrho); */ 
    
    
    if (conD*h==0.)     tempSQRTarg = 1.0 + conS2*BIG;
    else                tempSQRTarg = 1.0 + conS2/((conD*h)*(conD*h));
    Wlor    = sqrt(tempSQRTarg);
    
    // test for atm
    if (*rho<GRHD.ATM_FATM*GRHD.ATM_RHOATM)
        return ATM;
    
    // update tau
    *conT = Wlor*Wlor*(*rho)*h - *pres - conD;
    
    // test velocity
    v2 = conS2/((*conT+conD+*pres)*(*conT+conD+*pres));
    if (v2>GRHD.HRSC_VMAX) 
        return ATM;
    
    return OK;
}



// Routine for c2p inversion for general EoS p=P(rho,epsl) 
//   combined with cold, one-parameter, EoS p=P(epsl(rho),rho) 
//   Newton-Rapshon algorithm is based on pressure (general EoS) and rho (codl EoS) root   
//   It uses a cold and static atmosphere 
// Method
//   1.) try an inversion with a general EoS based on pressure root
//   2.) if failure assume cold EoS and invert looking for rho root 

void grhd_c2p_proot_ColdStaticAtm_hybrid(tVarList *ucur, tVarList *primVars)
{
    if (DEBUG) printf("# ===> grhd_c2p_proot_ColdStaticAtm_hybrid()\n");
    bampi_openmp_start
    
    tL *level = ucur->level;
  
    double *gxx = level->v[MATTER.INDX_VAR_gxx    ];
    double *gxy = level->v[MATTER.INDX_VAR_gxx + 1];
    double *gxz = level->v[MATTER.INDX_VAR_gxx + 2];
    double *gyy = level->v[MATTER.INDX_VAR_gxx + 3];
    double *gyz = level->v[MATTER.INDX_VAR_gxx + 4];
    double *gzz = level->v[MATTER.INDX_VAR_gxx + 5];

    double *MD  = vldataptr(ucur, MATTER.INDX_VAR_q    );
    double *MT  = vldataptr(ucur, MATTER.INDX_VAR_q + 1);
    double *MSx = vldataptr(ucur, MATTER.INDX_VAR_q + 2);
    double *MSy = vldataptr(ucur, MATTER.INDX_VAR_q + 3);
    double *MSz = vldataptr(ucur, MATTER.INDX_VAR_q + 4);

    double *Mrho  = vldataptr(primVars, 0);
    double *Mepsl = vldataptr(primVars, 1);
    double *Mvx   = vldataptr(primVars, 2);
    double *Mvy   = vldataptr(primVars, 3);
    double *Mvz   = vldataptr(primVars, 4);
    double *Mp    = vldataptr(primVars, 5); 
    double *Mv2   = vldataptr(primVars, 6); 
    double *detg  = vldataptr(primVars, 7);
    
    double conD,conSx,conSy,conSz,conT, conSix,conSiy,conSiz, conS2,SQRTconS2;
    double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg, ooDpT;
    double rho,epsl,pres, tmp1,tmp2,tmp3;

    const double errormax = Getd("grhd_C2P_NewtonRaphsonTR");
    const int countmax    = Getd("grhd_C2P_NewtonRaphsonNR");
    int result;
    

    forallpoints_ijk_openmp(level) {
        
        
      detg[ijk] = invg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
		       &gixx,&gixy,&gixz,&giyy,&giyz,&gizz);
      SQRTdetg = sqrt(fabs(detg[ijk]));
      if (!finite(SQRTdetg)) {
	detg[ijk] = 1.;
	SQRTdetg = 1.;
      }
      ooSQRTdetg = 1./SQRTdetg;

      /* 
      // this is a safer computation (testme):
      grhd_compute_detg_invg_pt(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
      &detg[ijk], 
      &gixx,&gixy,&gixz,&giyy,&giyz,&gizz); 
      SQRTdetg = sqrt(fabs(detg[ijk]));	
        ooSQRTdetg = (SQRTdetg==0.) ? BIG : 1./SQRTdetg;
      */

	conD   = MD[ijk]  *ooSQRTdetg;
        conSx  = MSx[ijk] *ooSQRTdetg;
        conSy  = MSy[ijk] *ooSQRTdetg;
        conSz  = MSz[ijk] *ooSQRTdetg;
        conT   = MT[ijk]  *ooSQRTdetg;
        conSix = gixx*conSx + gixy*conSy + gixz*conSz;
        conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
        conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
        conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
        //SQRTconS2 = sqrt(conS2);
	SQRTconS2 = sqrt(fabs(conS2));
        
        result = OK;
        
	
        // check ingoing variables
	/*
        if (CheckForNANandINF(10,conD,conSx,conSy,conSz,conSix,conSiy,conSiz,conT, SQRTconS2, SQRTdetg)) {
	  printf("conD, conS[xyz], conS_i^i, conT or detg is nan\n");
	  result = FAIL;
        }
	*/
        
	// more detailed checking:	
        // check evolution hydro variables
	if (CheckForNANandINF(5,conD,conSx,conSy,conSz,conT)) {
	  printf("conD, conS[xyz], or conT is nan\n");
	  result = FAIL;
        }
        
        // check ADM matter variables
	if (CheckForNANandINF(13,gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
			      detg[ijk], 
			      gixx,gixy,gixz,giyy,giyz,gizz)) {
	  printf("gxx (1:6), detg (7), or igxx (8:13) is nan\n");
	  result = FAIL; // here igxx is nan, what else to do?
        }
        
        // check evolution hydro variables	
	if (CheckForNANandINF(5,conSix,conSiy,conSiz,conS2,SQRTconS2)) {
	  printf("conSup[xyz], conS2, or SQRT(conS2) is nan, reset to zero, try to go on\n");
	  conSx = conSix = 0.;
	  conSy = conSiy = 0.;
	  conSz = conSiz = 0.;
	  conS2 = SQRTconS2 = 0.;	  
	  result = FAIL;	
	  //result = COLD; // testme
        }

        // try hot part
        if ((result!=FAIL))// && (result!=COLD))
            result = try_hot_part(conD,conS2,conT, &rho,&epsl,&pres, Mp[ijk], errormax,countmax);
        
	// hack used in bam 11	
	/*
	  if (result==NR) {
	  result = ATM;
	  }
	*/
	
        // try cold part
	if (result==COLD) {	  
	//if ((result==COLD) || (result==NR)) {
	  result = try_cold_part(conD,conS2,&conT, &rho,&epsl,&pres, Mrho[ijk], errormax,countmax);
	  MT[ijk] = SQRTdetg*conT;
        }
        
	// try to set ATM if something wrong here
        /* if (result!=OK) 
	   result = ATM; */
	if ((result!=OK) && (result!=ATM)) {
	  result = ATM;
	if(PRERROR)  printf("Problem in cons2prim: (%d)\n",result);	
	  // check detg
	  if (!finite(SQRTdetg)) {
	    detg[ijk] = 1.;
	    SQRTdetg = 1.;
	  }
	  // check if point is in finer level
	  if (ijkinsidefinerlevel(box,ijk)) {
	  if(PRERROR)  printf(" point is inside finer box/ in symmetry\n");                
	  } else {
	  if(PRERROR)  printf(" point is NOT inside finer box NOR in some symmetry area\n");
	  if(PRERROR)  printf(" (D,S2,T) =  %+e %+e %+e \n",conD,conS2,conT);
	    /* errorexit("better stop"); */	    
	  }
	  //printf(" set ATM and go on\n"); 
	}	    

        // set atm
        if (result==ATM) {
            pres     = GRHD.ATM_PATM;
            rho      = GRHD.ATM_RHOATM;
            epsl     = GRHD.ATM_EPSLATM;
            
            conD     = GRHD.ATM_RHOATM;
            conSx    = conSix = 0.;
            conSy    = conSiy = 0.;
            conSz    = conSiz = 0.;
            conS2    = 0.;
            conT     = GRHD.ATM_RHOATM*GRHD.ATM_EPSLATM;
            
            MD[ijk]  = SQRTdetg*conD;
            MSx[ijk] = SQRTdetg*conSx;
            MSy[ijk] = SQRTdetg*conSy;
            MSz[ijk] = SQRTdetg*conSz;
            MT[ijk]  = SQRTdetg*conT;
        }

                
        // set primitives
	ooDpT      = 1.0/(conD+pres+conT);
	//if (!finite(ooDpT)) ooDpT = 1.0/(conD+Mp[ijk]+conT); // try with old val
        Mrho[ijk]  = rho;
        Mepsl[ijk] = epsl;
        Mp[ijk]    = pres;
        Mvx[ijk]   = conSix * ooDpT;
        Mvy[ijk]   = conSiy * ooDpT;
        Mvz[ijk]   = conSiz * ooDpT;
	grhd_compute_v2_pt(Mvx[ijk],Mvy[ijk],Mvz[ijk],
			   gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk],
			   &tmp1,&tmp2,&tmp3, &(Mv2[ijk]) );
	
	// final check for nans
        if (CheckForNANandINF(12,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk],epsl)) 
	  printf("NaN in recovered primitives ! (grhd_c2p_proot_ColdStaticAtm_hybrid())\n");
        
	if ((PR) && result==FAIL) {
	  printf("failure in c2p v2 (grhd_c2p_proot_ColdStaticAtm_hybrid())\n");
	  printf("l=%d  x=%e y=%e x=%e  ijk=%d  detg=%e\n",level->l, Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],ijk,detg[ijk]);
	  if (ijkinsidefinerlevel(box,ijk)==0) {
	    printf("  point is NOT inside finer box NOR in some symmetry area\n"); 
	  } else 
	    printf("  point is inside finer box\n"); 
	}        
        
    } endfor_ijk_openmp;
    
    
    bampi_openmp_stop
}



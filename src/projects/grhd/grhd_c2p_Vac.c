/* grhd_c2p_Vac.c , Wolfgang Tichy 9/2015 */

#include "bam.h"
#include "grhd.h"


#define PR 0
#define DEBUG 0
#define PRERROR 0

#define EXITonNAN 1

enum {
    COLD,
    ATM,
    OK,
    FAIL,
    NR,
    VAC
};


/* structure we pass to rtsafe_itsP root finder */
typedef struct {

  double conD, conT, conS2;  /* cons vars needed */
  double rho, epsl, peos  ;  /* prim vars */
  int ret;                   /* return value */
  
} tCons_rho_epsl;


/* ************************************************************** */
/* WT: Below are some c2p routines that assume true vacuum        */
/* ************************************************************** */


/* Function po-peos. It has to be zero when po is correct pressure.
   used by rtsafe_itsP in c2p_finder_Vac */
void p_minus_peos(double po, double *f, double *df, void *par)
{
  tCons_rho_epsl *Cre = (tCons_rho_epsl *) par;
  int problem;
  double conD = Cre->conD;
  double conT = Cre->conT;
  double conS2 = Cre->conS2;
  double temp,tempSQRTarg,tempSQRT, Wlor;
  double peos, rho, epsl;
  double drhodp,depsldp,kappa,chi,cs2;
  
  temp        = conT+po+conD;      /* temp>0 if po>=pmin */
  tempSQRTarg = temp*temp-conS2;   /* tempSQRTarg > 0 if po>=pmin */
  tempSQRT    = sqrt(tempSQRTarg); /* this should work if po>=pmin */
  Wlor        = temp/tempSQRT;     /* this should work if po>=pmin */
  rho         = conD/Wlor;
  epsl        = (tempSQRT - Wlor*po)/conD - 1.;

  /* if EOS is s.t. p>=0, WEC says epsl >= -1 */
  if(epsl <= -1.0)  epsl = -0.9999999999;

  // EOS call
  problem = EOS.comp("re","","","pc","re","",
                     rho,epsl, &peos,&cs2, &chi,&kappa);
  /* if(problem) errorexit("EOS.comp failed"); */
    
  *f = po - peos;
  /* derivfs for Newton */
  drhodp  = conD*conS2/(tempSQRT *temp*temp);
  depsldp = po*conS2/(conD*tempSQRTarg*tempSQRT);
  *df     = 1. - chi*drhodp - kappa*depsldp;
  if(*df==0.) errorexit("Newton: *df=0");

  /* set some other vars */
  Cre->ret = OK;
  Cre->rho = rho;
  Cre->epsl = epsl;
  Cre->peos = peos;
}


/* c2p using full EoS, using rtsafe_itsP */
int c2p_finder_Vac(double conD, double conS2, double conT,
                   double *rho, double *epsl, double *pres,
                   double *pMin, double *pEOS, 
                   double pinit, double errormax,double countmax)
{
  int ret;
  double temp,tempSQRTarg,tempSQRT;
//  , Wlor;
  double pmin,pmax, pacc;
//  double f,df,drhodp,depsldp,kappa,chi,cs2,v2;
//  double problem,error,lostdigits;
  tCons_rho_epsl Cre[1];

#if(EXITonNAN)  // WT FIXME: remove:
  if (CheckForNANandINF(3,conD,conS2,conT))
  {
    printf("conD, conS2, or conT is nan\n");
    errorexit("found NAN/INF");
  }
#endif
  /* check for vacuum */
  if(conD<=0.0)
  {
    *rho = *epsl = *pres = *pMin = *pEOS = 0.0;
    return VAC;
  }
  /* WEC says: Tau >= (W^2 - 1) p - D and Tau >= -p - D,
     but if  EOS says: epsl >= 0, p >= 0  ==>  Tau >= 0 */
  ///* ensure that conT>=0 */
  /* this seems to lead to NANs in con[D,S,T] after evo step */
  //if(conT<0.0) conT=0.0;


  /* set minimal p, s.t. (conT+pmin+conD)^2 >= conS2, 
     i.e. s.t. Lorentz factor Wlor^(-2) >= 0  */
  if(conS2>=0) pmin = sqrt(conS2)-conT-conD;
  else         errorexit("conS2<0");
  if(pmin<0.0) pmin = 0.0;   /* pmin must be >= 0 */
  pmin += conD*1e-8 + 1e-40; /* add a bit s.t. Wlor^(-2) > 0 */
  *pMin = pmin; /* save min p */

  /* set initial guess for pres */
  *pres = DMAX(pinit, pmin);

  /* test if our *pres works. If not set all to vaccum */
  temp        = conT+*pres+conD;    /* temp>0 if *pres>=pmin */
  tempSQRTarg = temp*temp-conS2;   /* tempSQRTarg > 0 if *pres>=pmin */
  if(tempSQRTarg<=0.0)
  {
    *rho = *epsl = *pres = *pMin = *pEOS = 0.0;
    return VAC;
  }

  /* set up pars for root finder */
  Cre->conD  = conD;
  Cre->conT  = conT;
  Cre->conS2 = conS2;

  // /* get pressure max at epsl=0 to get pressure scale */
  //ret = EOS.comp("re","","","p","","", GRHD.ATM_RHOMAX,0.0, &pmax);
  pmax = 1.009862492202877*GRHD.ATM_PMAX;
  pacc = errormax*pmax; /* set accuracy */

  /* set maximal p */
  pmax = 1e3*DMAX(*pres, pmax);

  /* call root finder */
  ret = rtsafe_itsP(pres, p_minus_peos, pmin,pmax, Cre, countmax,pacc);

  /* set other output vars */
  *rho = Cre->rho;
  *epsl = Cre->epsl;
  *pEOS = Cre->peos;
  *pMin = pmin;
  

#if(EXITonNAN)
  /* FIXME: remove this, we should never get NANs */
  if (!finite(*pres) || !finite(*rho) || !finite(*epsl))
    errorexit("!finite(*pres) || !finite(*rho) || !finite(*epsl)");
#endif

  /* if rtsafe_itsP reported an error */
  if(ret<0)
  {
    return VAC;
  }

  /* too many steps? */
  if(ret>countmax)
  {
    return NR;
  }

  return OK;
}



/* c2p using full EoS, using custom root finder */
int hot_part_Vac(double conD, double conS2, double conT,
                 double *rho, double *epsl, double *pres,
                 double *pMin, double *pEOS, 
                 double pinit, double errormax,double countmax)
{
  int count;
  double temp,tempSQRTarg,tempSQRT, Wlor;
  double pold,peos,pmin;
  double f,df,drhodp,depsldp,kappa,chi,cs2;
  double problem,error,lostdigits;

#if(EXITonNAN)  // WT FIXME: remove:
  if (CheckForNANandINF(3,conD,conS2,conT))
  {
    printf("conD, conS2, or conT is nan\n");
    errorexit("found NAN/INF");
  }
#endif
  /* check for vacuum */
  if(conD<=0.0)
  {
    *rho = *epsl = *pres = *pMin = *pEOS = 0.0;
    return VAC;
  }
  /* WEC says: Tau >= (W^2 - 1) p - D and Tau >= -p - D,
     but if  EOS says: epsl >= 0, p >= 0  ==>  Tau >= 0 */
  ///* ensure that conT>=0 */
  /* this seems to lead to NANs in con[D,S,T] after evo step */
  //if(conT<0.0) conT=0.0;

  /* set minimal p, s.t. (conT+pmin+conD)^2 >= conS2, 
     i.e. s.t. Lorentz factor Wlor^(-2) >= 0  */
  if(conS2>=0) pmin = sqrt(conS2)-conT-conD;
  else         errorexit("conS2<0");
  if(pmin<0.0) pmin = 0.0;   /* pmin must be >= 0 */
  pmin += conD*1e-8 + 1e-40; /* add a bit s.t. Wlor^(-2) > 0 */
  *pMin = pmin; /* save min p */
  /* set pold */
  pold = DMAX(pinit, pmin);

  /* start Newton Raphson and find pres */
  for(count=0, error=2*errormax; count<=countmax; count++)
  {
    temp        = conT+pold+conD;    /* temp>0 if pold>=pmin */
    tempSQRTarg = temp*temp-conS2;   /* tempSQRTarg > 0 if pold>=pmin */
    /* test if our pold works. If not set all to vaccum */
    if( tempSQRTarg<=0.0 || !(isfinite(pold)) )
    {
      *rho = *epsl = *pres = *pMin = *pEOS = 0.0;
      return VAC;
    }
    tempSQRT    = sqrt(tempSQRTarg); /* this should work if pold>=pmin */
    Wlor        = temp/tempSQRT;     /* this should work if pold>=pmin */
    *rho        = conD/Wlor;
    *epsl       = (tempSQRT - Wlor*pold)/conD - 1.;

#if(EXITonNAN) // WT FIXME: remove:
    if (CheckForNANandINF(3,Wlor,*rho,*epsl))
    {
      printf("nan\n");
      errorexit("found NAN/INF");
    }
#endif

    /* if EOS is s.t. p>=0, WEC says epsl >= -1 */
    if(*epsl <= -1.0)  *epsl = -0.9999999999;

    // EOS call
    problem = EOS.comp("re","","","pc","re","",
                       *rho,*epsl, &peos,&cs2, &chi,&kappa);
    if(problem) errorexit("EOS.comp failed");
    
    f       = pold - peos;
    /* derivfs for Newton */
    drhodp  = conD*conS2/(tempSQRT *temp*temp);
    depsldp = pold*conS2/(conD*tempSQRTarg*tempSQRT);
    df      = 1. - chi*drhodp - kappa*depsldp;
    // if(df==0.) errorexit("Newton: df=0");

    /* Newton step */
    *pres   = pold-f/df;
    /* set *pres to pmin if step goes below pmin */
    *pres   = DMAX( pmin, *pres );

    /* diff between this and previous pressure value */
    /* error becomes small if we find solution, or if we get stuck at pmin */
    error   = fabs( 1.-pold/(*pres) );
    pold    = *pres;
    /* break if error is small */
    if(error<errormax) break; 
    
    if((DEBUG) && (count>2))
    {
      printf("    %d   %+2.16e   %+2.16e %+2.16e %2.16e %+2.16e %2.16e\n", count,error,df,f,*pres,*epsl,*rho);
      printf("         => %+2.16e %+2.16e\n", chi*drhodp,kappa*depsldp);
      printf("         => %+2.16e %+2.16e\n",  pold,f/df);
      printf("         => %+2.16e %+2.16e %+2.16e\n", conT,pold,conD);
      printf("         => %+2.16e %+2.16e\n",  temp*temp,conS2);
      printf("         => %+2.16e %+2.16e\n",  tempSQRT,Wlor*pold);
      printf("         => %+2.16e %+2.16e\n",  (tempSQRT - Wlor*pold)/conD,1.);
    }
  } /* end for loop */

  /* p computed from EOS */
  *pEOS = peos;
  
  if (count>=countmax)
  {
    lostdigits = 0;
    lostdigits = DMAX(lostdigits , fabs(log10(conD)-log10(conT)));
    lostdigits = DMAX(lostdigits , fabs(log10(conD)-log10(*pres)));
    lostdigits = DMAX(lostdigits , fabs(log10(conT)-log10(*pres)));
    //lostdigits = DMAX(lostdigits , fabs(log10(temp*temp)-log10(conS2)));
    lostdigits = DMAX(lostdigits , fabs(log10(tempSQRT)-log10(Wlor*pold)));
    //lostdigits = DMAX(lostdigits , fabs(log10((tempSQRT - Wlor*pold)/conD)-log10(1.)));
    //lostdigits = DMAX(lostdigits , fabs(log10(chi*drhodp)-log10(kappa*depsldp)));
    //printf("  ==> %2.2f    (%2.2f)\n", lostdigits, log10(error)-log10(1e-16));
    if ( lostdigits  <  log10(error)-log10(1e-15) || error > 1e-4 )
    {
      if(PR) 
      {
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

#if(EXITonNAN)
  /* FIXME: remove this, we should never get NANs */
  if (!finite(*pres) || !finite(*rho) || !finite(*epsl))
    errorexit("!finite(*pres) || !finite(*rho) || !finite(*epsl)");
#endif
  /* check true error */
  f     = pold - peos;
  error = fabs( f/pold );
  if(error>errormax)
  {
    if(pold == pmin) return VAC;
    else return FAIL;
  }
  //if (*epsl<0.) return VAC;
  
  return OK;
}


/* Routine for c2p inversion for general EoS p=P(rho,epsl) 
   combined with cold, one-parameter, EoS p=P(epsl(rho),rho) 
   Newton-Rapshon algorithm is based on pressure (general EoS) and 
   rho (codl EoS) root   
   It uses a true vacuum!
 Method
   1.) try an inversion with a general EoS based on pressure root
   2.) if failure assume cold EoS and invert looking for rho root */
void grhd_c2p_proot_Vac(tVarList *ucur, tVarList *primVars)
{
  if (DEBUG) printf("# ===> grhd_c2p_proot_Vac()\n");
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
  double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg, DpT,ooDpT;
  double rho,epsl,pres, pmin,peos, tmp1,tmp2,tmp3;

  const double errormax = Getd("grhd_C2P_NewtonRaphsonTR");
  const int countmax    = Getd("grhd_C2P_NewtonRaphsonNR");
  int result;

  forallpoints_ijk_openmp(level)
  {
    detg[ijk] = invg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
		       &gixx,&gixy,&gixz,&giyy,&giyz,&gizz);
    SQRTdetg = sqrt(fabs(detg[ijk]));
    if(detg[ijk]<=0.0)
    {
      printf("grhd_c2p_proot_Vac:  detg[ijk]=%g <= 0\n", detg[ijk]);
      printf(" (x,y,z)=(%g,%g,%g)  ijk=%d  l=%d\n",
             Ptr(level,"x")[ijk], Ptr(level,"y")[ijk], Ptr(level,"z")[ijk],
             ijk, level->l);
      errorexit("detg[ijk]<=0.0");
      printf("  => setting: detg[ijk] = 1\n");
      detg[ijk] = 1.;
      SQRTdetg = 1.;
    }
    ooSQRTdetg = 1./SQRTdetg;

    conD   = MD[ijk]  *ooSQRTdetg;
    conSx  = MSx[ijk] *ooSQRTdetg;
    conSy  = MSy[ijk] *ooSQRTdetg;
    conSz  = MSz[ijk] *ooSQRTdetg;
    conT   = MT[ijk]  *ooSQRTdetg;
    conSix = gixx*conSx + gixy*conSy + gixz*conSz;
    conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
    conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
    conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
    if(conS2<0.0) errorexit("conS2<0.0");
    SQRTconS2 = sqrt(fabs(conS2));
    
    result = OK;

#if(EXITonNAN)
    /* FIXME: disable all CheckForNANandINF later,
       we should never have to deal with this, if we stop putting in NANs
       all over the place! */
    // check evolution hydro variables
    if (CheckForNANandINF(5,conD,conSx,conSy,conSz,conT))
    {
      printf("conD, conS[xyz], or conT is nan\n");
      result = FAIL;
      errorexit("found NAN/INF");
    }
    
    // check ADM matter variables
    if (CheckForNANandINF(13,gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
                          detg[ijk], 
                          gixx,gixy,gixz,giyy,giyz,gizz))
    {
      printf("gxx (1:6), detg (7), or igxx (8:13) is nan\n");
      result = FAIL; // here igxx is nan, what else to do?
      errorexit("found NAN/INF");
    }
    
    // check evolution hydro variables        
    if (CheckForNANandINF(5,conSix,conSiy,conSiz,conS2,SQRTconS2))
    {
      printf("conSup[xyz], conS2, or SQRT(conS2) is nan, reset to zero, try to go on\n");
      conSx = conSix = 0.;
      conSy = conSiy = 0.;
      conSz = conSiz = 0.;
      conS2 = SQRTconS2 = 0.;          
      result = FAIL;        
      //result = COLD; // testme
      errorexit("found NAN/INF");
    }
#endif

    /* find pres (and along with it rhp, epls) with Newton solver */
    if(result!=FAIL)
      result = c2p_finder_Vac(conD,conS2,conT, &rho,&epsl,&pres, 
                              &pmin, &peos,
                              Mp[ijk], errormax,countmax);
    /* result = hot_part_Vac(conD,conS2,conT, &rho,&epsl,&pres, 
                            &pmin, &peos,
                            Mp[ijk], errormax,countmax); */

    /* if we do not get OK, the cons vars cannot be physical */
    if(result!=OK)
    {
      /* now set cons vars to vaccum */
      /* recall: 
         conD   = MD[ijk]  *ooSQRTdetg;
         conSx  = MSx[ijk] *ooSQRTdetg;
         conSy  = MSy[ijk] *ooSQRTdetg;
         conSz  = MSz[ijk] *ooSQRTdetg;
         conT   = MT[ijk]  *ooSQRTdetg;
         conSix = gixx*conSx + gixy*conSy + gixz*conSz;
         conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
         conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
         conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;  */
      MD[ijk] = MT[ijk] = 0.0;
      MSx[ijk] = MSy[ijk] = MSz[ijk] = 0.0;
      conD = conT = 0.0;
      conSx = conSy = conSz = conSix = conSiy = conSiz = conS2 = 0.0;
      /* set prims to vacuum as well */
      rho = epsl = pres = 0.0;
      result = VAC;
    }

    /* set primitives */
    DpT = conD+pres+conT;
    /* if vac set, DpT to dummy val that will ensure vac for all prims below */
    if(result==VAC) DpT = 1.0;
    if(DpT==0.0) DpT = conD+Mp[ijk]+conT; // try with old val
    if(DpT==0.0) errorexit("DpT=0.0");
    ooDpT      = 1.0/DpT;
    Mrho[ijk]  = rho;
    Mepsl[ijk] = epsl;
    Mp[ijk]    = pres;
    Mvx[ijk]   = conSix * ooDpT;
    Mvy[ijk]   = conSiy * ooDpT;
    Mvz[ijk]   = conSiz * ooDpT;
    grhd_compute_v2_pt(Mvx[ijk],Mvy[ijk],Mvz[ijk],
                       gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk],
                       &tmp1,&tmp2,&tmp3, &(Mv2[ijk]) );

#if(EXITonNAN)
    /* FIXME: remove this, we should not have NANs !!! */
    // final check for nans
    if (CheckForNANandINF(12,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk],epsl))
    {
      printf("NaN in recovered primitives ! (grhd_c2p_proot_Vac())\n");
      errorexit("NaN in recovered primitives! (grhd_c2p_proot_Vac())");
    }
#endif
    if(result==FAIL)
    {
      printf("failure in c2p v2 (grhd_c2p_proot_Vac())\n");
      printf("l=%d  x=%e y=%e x=%e  ijk=%d  detg=%e\n",level->l, Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],ijk,detg[ijk]);
      if (ijkinsidefinerlevel(box,ijk)==0)
        printf("  point is NOT inside finer box NOR in some symmetry area\n"); 
      else 
        printf("  point is inside finer box\n");
      errorexit("grhd_c2p_proot_Vac failed!!!");
    }
  } endfor_ijk_openmp;

  bampi_openmp_stop
}

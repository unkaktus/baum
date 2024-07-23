/* grhd_c2p.c 
   sbernuz, mth 03/2012 */




#include "bam.h"
#include "grhd.h"


#define PR 0
#define DEBUG 0
#define ln10 2.302585092994046e+00


// to perform an additional c2p befor output
int grhd_make_c2p(tL *level)
{
  
  if (level->shells)
    return 1;
  
  tVarList *ucur     = get_evolve_vlregister(level);
  tVarList *primVars = give_primVars_vl(level);
  
  // this is needed, I do not know why and if this is correct 
   ucur->level = level;
  
  if (DEBUG) {
    int i;
    for (i=0; i<ucur->n; i++) printf("%d  %s\n",i,VarName(ucur->index[i]));
    //printf("%p %p  %d\n",evolve_rhs_compute_adm, matter_c2p, MATTER_INDX_VAR_q);
    for (i=0; i<primVars->n; i++) printf("%d  %s\n",i,VarName(primVars->index[i]));
    printf("%p %p   %d %d\n",level,ucur->level, level->l,ucur->level->l);
  }
  
  if (PR) printf(" compute ADM   level %d\n",level->l);
  //evolve_rhs_compute_adm(level, ucur);
  
  if (PR) printf(" c2p   level %d\n",level->l);
  MATTER.c2p(ucur, primVars);
  
  vlfree(primVars);
  
  return 1;
}



/* ************************************ 
   macros and functions */


// compute rho from pressure and conservatives 
#define compute_rho(p)							\
  temp = conT+(p)+conD;							\
  tempSQRTarg = temp*temp-conS2;				\
  tempSQRT = sqrt(tempSQRTarg);				\
  Wlor        = (temp)/tempSQRT;				\
  rho         = conD/Wlor;

// compute epsl from pressure and conservatives 
#define  compute_epsl(p)			\
  epsl = (tempSQRT - Wlor*(p))/conD - 1.;  

// compute h(rho) and h(rho)' from rho and conservatives 
// assume cold 1-paramater EoS
// the EoS interface is not called but hardcoded here below for now
#define compute_enthal(rho) \
  EOS.comp("r","","","pec","re","", \
            (rho), &pres,&epsl,&cs2, &chi,&kappa); \
  h       = 1.0 + epsl + pres/(rho);\
  dhdrho  = chi/(rho);\
  
//SB the above line uses thermodynamics at T=0,  
//   general expression is:
//   dhdrho = depsldrho - pres/((rho)*(rho)) + chi/(rho);	

// compute Lorentz factor form rho and conservatives
#define compute_Wlor(rho)                                           \
  if( conD==0. || h==0.)					    \
    tempSQRTarg = 1.0 + conS2*1e+15;				    \
  else								    \
    tempSQRTarg = 1.0 + conS2/((conD*h)*(conD*h));		    \
  Wlor  = sqrt(tempSQRTarg);					    \
  





/* ************************************ 
   c2p routines */





// grhd_c2p_proot
//   - c2p for general EoS 
//   - Root of P with Newthon-Raphson method
void grhd_c2p_proot(tVarList *ucur, tVarList *primVars)
{

  if (PR) printf("# ===> grhd_c2p_proot()\n");


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
  /* double *Mv2   = vldataptr(primVars, 6); */    

  double *detg  = vldataptr(primVars, 7);

  double conD,conSx,conSy,conSz,conT, conSix,conSiy,conSiz, conS2,SQRTconS2;
  double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg;
  double epsl,rho, p,pold,pnew,pmin,peos, kappa,chi, cs2;
  double temp,tempSQRTarg,tempSQRT,Wlor, f,drhodp,depsldp,df;

  int count;
  double error,r;
  
  const double errormax = Getd("grhd_C2P_NewtonRaphsonTR");
  const int countmax    = Getd("grhd_C2P_NewtonRaphsonNR");


  forallpoints_ijk(level) {
        
    detg[ijk] = invg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], &gixx,&gixy,&gixz,&giyy,&giyz,&gizz);
    if ((detg[ijk]<=0.) || (!finite(detg[ijk]))) {
      printf("%e %e %e   %d   %e\n", Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],ijk,detg[ijk]);
      errorexit("problem with determinant");
    }
    SQRTdetg = sqrt(detg[ijk]);
    ooSQRTdetg = 1./SQRTdetg;
    
    conD   = MD[ijk] *ooSQRTdetg;
    conSx  = MSx[ijk]*ooSQRTdetg;
    conSy  = MSy[ijk]*ooSQRTdetg;
    conSz  = MSz[ijk]*ooSQRTdetg;
    conT   = MT[ijk] *ooSQRTdetg;
    
    conSix = gixx*conSx + gixy*conSy + gixz*conSz;
    conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
    conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
    conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
    SQRTconS2 = sqrt(conS2);
    
    if (CheckForNANandINF(6,conD,conSx,conSy,conSz,conT, SQRTconS2)) errorexit("conD, conS[xyz], conS_i^i or conT is nan");
    
    pmin = SQRTconS2-conD-conT; 
    
    // guess pold from previous step 
    pold = Mp[ijk];
    
    pold = DMAX(pold,pmin);
    compute_rho(pold);
    compute_epsl(pold);
    
    if (!finite(pold)) errorexit("starting guess for p is not finite");
    
    // check if cons vars have physical values
    if (tempSQRTarg<0.0)
      errorexit("negative SQRT arg: bad guess or conservative unphysical");
    
    // start Newton Raphson and find p
    for (count=0,error=2*errormax; count<=countmax; count++) {
      
      //problem = GRHD.use_eos(&peos, &cs2, &kappa, &chi, rho,epsl);
      EOS.comp("re","","","pc","re","",
               rho,epsl, &peos,&cs2,&chi,&kappa);
      
      f       = pold - peos; 
      drhodp  = conD*conS2/(tempSQRT *temp*temp);
      depsldp = pold*conS2/(conD*tempSQRTarg*tempSQRT);
      df      = 1. - chi*drhodp - kappa*depsldp;
      if (df==0.0) errorexit("df == zero");
      
      pnew    = DMAX( pmin, pold-f/df );
      error   = fabs( 1.-pold/pnew );
      pold    = pnew;
      
      compute_rho(pold);
      compute_epsl(pold);

      if (error<errormax) break;
    }
    
    if (PR) if (count>=countmax) {
      printf("%d  %e %e    %e %e %e\n",count,error, errormax, f/df, f,df);
      printf("  %e %e %e => r=%e \n", Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],
              sqrt(pow(Ptr(level,"x")[ijk],2)+pow(Ptr(level,"y")[ijk],2)+pow(Ptr(level,"z")[ijk],2)));
      printf("  %+e %+e %+e (values  rho,epsl,p)\n",rho,epsl,pnew);
      printf("  %+e %+e %+e (D,S2,T)\n",conD,conS2,conT);
      printf("  %+e %+e %+e %+e (vx,vy,vz,v)\n",conSix/(conD+pnew+conT),conSix/(conD+pnew+conT),conSix/(conD+pnew+conT),
              sqrt(conSix*conSix+conSiy*conSiy+conSiz*conSiz)/(conD+pnew+conT));
      errorexit("too many Newton-Raphson steps");
    }
    
    // check on internal energy
    if (epsl<0.) {
      printf("  %e %e %e => r=%e \n", Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk], 
              sqrt(pow(Ptr(level,"x")[ijk],2)+pow(Ptr(level,"y")[ijk],2)+pow(Ptr(level,"z")[ijk],2))); 
          errorexit("negative internal energy.");
    }
    
    // set primitive values
    Mrho[ijk]  = rho;
    Mepsl[ijk] = epsl;
    Mp[ijk]    = pnew;
    Mvx[ijk]   = conSix/(conD+pnew+conT);
    Mvy[ijk]   = conSiy/(conD+pnew+conT);
    Mvz[ijk]   = conSiz/(conD+pnew+conT);
    
    CheckForNANandINF(11,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk]);
    
  } endfor_ijk;

}

// grhd_c2p_vanal
//   - c2p for Gamma-law EoS
//   - Analytic recovery based on quartic root
void grhd_c2p_vanal(tVarList *ucur, tVarList *primVars)
{
  if (PR) printf("# ===> grhd_c2p_vanal()\n");


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
  /* double *Mv2   = vldataptr(primVars, 6); */

  double *detg  = vldataptr(primVars, 7);

  double conD,conSx,conSy,conSz,conT, conSix,conSiy,conSiz, conS2,SQRTconS2;
  double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg;

  double S,S2,v,E,SE,D,D2, Gamma,Gammamosq;
  double a0,a1,a2,a3, i1,i2,i3, iR,iS,iT, iB,iC,ix1=0.;
  double div, rho,vx,vy,vz,epsl,p;
  const double gammamo = Getd("eos_Gamma")-1.; // 5/3-1;
  const double oo3 = 1./3.;


  forallpoints_ijk(level) {

    detg[ijk] = invg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], &gixx,&gixy,&gixz,&giyy,&giyz,&gizz);
    if ((detg[ijk]<=0.) || (!finite(detg[ijk]))) {
      printf("%e %e %e   %d   %e\n", Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],ijk,detg[ijk]);
      errorexit("problem with determinant");
    }
    SQRTdetg = sqrt(detg[ijk]);
    ooSQRTdetg = 1./SQRTdetg;

    conD   = MD[ijk] *ooSQRTdetg;
    conSx  = MSx[ijk]*ooSQRTdetg;
    conSy  = MSy[ijk]*ooSQRTdetg;
    conSz  = MSz[ijk]*ooSQRTdetg;
    conT   = MT[ijk] *ooSQRTdetg;

    conSix = gixx*conSx + gixy*conSy + gixz*conSz;
    conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
    conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
    conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
    SQRTconS2 = sqrt(conS2);

    if (CheckForNANandINF(6,conD,conSx,conSy,conSz,conT, SQRTconS2)) errorexit("conD, conS[xyz], conS_i^i or conT is nan");

    D  = conD; // redundant, fixme
    S2 = conS2;
    S  = SQRTconS2;
    v  = 0.0;
    E  = conT+conD;

    if (fabs(S) > 1e-20) {

      SE = S*E;
      D2 = D*D;

      Gamma     = gammamo+1;
      Gammamosq = gammamo*gammamo;

      div = 1.0 / (Gammamosq * (S2 + D2));

      a3 = (-2.0 * Gamma * gammamo * SE) * div;
      a2 = (Gamma*Gamma * E*E + 2.0*gammamo*S2 - Gammamosq*D2)*div;
      a1 = (-2.0 * Gamma * SE) * div;
      a0 = S2 * div;

      i1 = -a2;
      i2 = a3 * a1 - 4.0 * a0;
      i3 = 4.0 * a2 * a0 - a1*a1 - a3*a3 * a0;

      iR = (9.0 * i1 * i2 - 27.0 * i3 - 2.0 * i1*i1 * i1) / 54.0;
      iS = (3.0 * i2 - a2*a2) / 9.0;
      iT = iR*iR + iS*iS * iS;

      if (iT < 0.)
       ix1 = 2.0*pow(sqrt(iR*iR + iT),(oo3))*cos(atan2(sqrt(-iT),iR)*oo3) - i1*oo3;
      else
      ix1 = pow((iR + sqrt(iT)),(oo3)) + pow((iR - sqrt(iT)),(oo3)) - i1*oo3;

      iB = 0.5*(a3 + sqrt(a3*a3 - 4.0*a2 + 4.0*ix1));
      iC = 0.5*(ix1 - sqrt(ix1*ix1 - 4.0*a0));
      v  = (-iB + sqrt(iB*iB - 4.0*iC))*0.5;

    }

    v   = DMAX(v, 0.0);
    v   = DMIN(v, 1.0 - 1.0e-15);
    div = v / S;
    if (fabs(S) < 1e-20)
      div = 0.0;

    rho = sqrt(1.0 - v*v) * D;
    vx  = conSx * div;
    vy  = conSy * div;
    vz  = conSz * div;
    p   = gammamo*((E - conSx*vx - conSy*vy - conSz*vz) - rho);
    epsl =  p/(gammamo*rho);

    // check on internal energy
    if (epsl<0.) {
      printf("  %e %e %e => r=%e \n", Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],
            sqrt(pow(Ptr(level,"x")[ijk],2)+pow(Ptr(level,"y")[ijk],2)+pow(Ptr(level,"z")[ijk],2)));
      errorexit("negative internal energy.");
    }

    // set primitive values
    div = 1./(conD+p+conT);
    Mrho[ijk]  = rho;
    Mepsl[ijk] = epsl;
    Mp[ijk]    = p;
    Mvx[ijk]   = conSix*div;
    Mvy[ijk]   = conSiy*div;
    Mvz[ijk]   = conSiz*div;

    CheckForNANandINF(11,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk]);

  } endfor_ijk;

}



// grhd_c2p_rroot_ColdStaticAtm
//   - c2p for cold eos, cold-static atm treatment
//   - Root of rho with Newthon-Raphson method
void grhd_c2p_rroot_ColdStaticAtm(tVarList *ucur, tVarList *primVars)
{

  if (PR) printf("# ===> grhd_c2p_rroot_ColdStaticAtm()\n");
  
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
  double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg;
  double epsl, cs2, kappa, rhold,rho, pres, chi;
  double tempSQRTarg,Wlor, g,dg, h,dhdrho, vel2, tmp1,tmp2,tmp3;
  
  int count,problem,atm, prob;
  double error;
  char errorstring[10000] = {0};

  const double errormax = Getd("grhd_C2P_NewtonRaphsonTR");
  const int countmax    = Getd("grhd_C2P_NewtonRaphsonNR");  
  const double rhomin   = GRHD.ATM_RHOATM;

  
  forallpoints_ijk_openmp(level) {
    
    prob  = 0;
    atm   = 0;
    
    detg[ijk] = invg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
		     &gixx,&gixy,&gixz,&giyy,&giyz,&gizz);
    if ((detg[ijk]<=0.) || (!finite(detg[ijk]))) {
      aprintf(errorstring,"  problem with determinant\n");
      //detg[ijk] = 1.;
      detg[ijk] =  DETGMIN;
      aprintf(errorstring,"    cured it by detg=1\n");
      prob++;
    }
    SQRTdetg = sqrt(detg[ijk]);
    ooSQRTdetg = 1./SQRTdetg;
    
    conD   = MD[ijk] *ooSQRTdetg;
    conSx  = MSx[ijk]*ooSQRTdetg;
    conSy  = MSy[ijk]*ooSQRTdetg;
    conSz  = MSz[ijk]*ooSQRTdetg;
    conT   = MT[ijk] *ooSQRTdetg;
    conSix = gixx*conSx + gixy*conSy + gixz*conSz;
    conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
    conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
    conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
    SQRTconS2 = sqrt(conS2);
    if (CheckForNANandINF(6,conD,conSx,conSy,conSz,conT, SQRTconS2)) {
      aprintf(errorstring,"  conD, conS[xyz], conS_i^i or conT is nan\n");
      prob++;
    }
    
    // guess rhoold from previous step 
    rhold = Mrho[ijk];
    
    // check if cons vars have physical values
    // try cure => start c2p from atm	 
    if (conD <= 0.)  {
      
      conD     = GRHD.ATM_RHOATM;
      conSx    = conSix = 0.;
      conSy    = conSiy = 0.;
      conSz    = conSiz = 0.;
      conS2    = 0.;	
      conT     = GRHD.ATM_RHOATM*GRHD.ATM_EPSLATM;
      
      //SB 21/05/2010 DO NOT RESET
      /*  
	  MD[ijk]  = SQRTdetg*conD;
	  MSx[ijk] = SQRTdetg*conSx;
	  MSy[ijk] = SQRTdetg*conSy;
	  MSz[ijk] = SQRTdetg*conSz;
	  MT[ijk]  = SQRTdetg*conT;
      */
    }
    
    rhold = DMAX(rhold,rhomin);
    
    if (!finite(rhold)) {
      aprintf(errorstring,"  starting guess for rho is not finite\n");
      prob++;
    }
    
    // start Newton Raphson and find rho
    for (count=0,error=2*errormax; count<=countmax; count++) {
      
      // the EoS call is hardcoded in these macros
      compute_enthal(rhold);
      compute_Wlor(rhold);
      
      g   = Wlor*rhold - conD; 
      dg  = Wlor - (rhold*conS2*dhdrho)/(Wlor*conD*conD*h*h*h);
      if (dg==0.0) {
	aprintf(errorstring,"  dg == zero");
	prob++;
      }
      
      rho     = DMAX( rhomin, rhold-g/dg );
      error   = fabs( 1.-rhold/rho );
      rhold   = rho;
      
      if (error<errormax) break;
    }
    
    compute_enthal(rho);
    compute_Wlor(rho);
    
    if (count>=countmax) {
      aprintf(errorstring,"  too many Newton-Raphson steps\n");
      aprintf(errorstring,"    %d  %e %e    %e %e %e\n",count,error, errormax, g/dg, g,dg);
      aprintf(errorstring,"    %+e %+e %+e (values  rho,epsl,p)\n",rho,epsl,pres);
      aprintf(errorstring,"    %+e %+e %+e (atm)\n",GRHD.ATM_RHOATM,GRHD.ATM_EPSLATM,GRHD.ATM_PATM);
      aprintf(errorstring,"    %+e %+e %+e (D,S2,T)\n",conD,conS2,conT);
      aprintf(errorstring,"    %+e %+e %+e %+e (vx,vy,vz,v)\n",conSix/(conD+pres+conT),conSix/(conD+pres+conT),conSix/(conD+pres+conT),
	      sqrt(conSx*conSix+conSy*conSiy+conSz*conSiz)/(conD+pres+conT));
      prob++;
    }
    
    // check atmosphere for final value
    if (rho<GRHD.ATM_FATM*GRHD.ATM_RHOATM) atm++; 
    
    // re-set Tau
    conT    = Wlor*Wlor*rho*h - pres - conD;
    MT[ijk] = SQRTdetg*conT;
    
    // at this point we found a soultion BUT 
    // if v>=1 it is not physical -> set to ATM
    vel2 = conS2/((conT+conD+pres)*(conT+conD+pres));
    if (vel2>GRHD.HRSC_VMAX) atm++;
    
    // set atm or test computed values
    if (!atm) {
      
      // check for finite values
      if (!finite(rho) || !finite(pres) || !finite(epsl)) {
	aprintf(errorstring,"  non-finite values of primitives\n");
	aprintf(errorstring,"    pres = %e    rho  = %e    epsl = %e     %e\n",pres,rho,epsl,GRHD.ATM_RHOATM);
	prob++;
      }
      
      // check epsl
      if (epsl<0.) {
	aprintf(errorstring,"  negative internal energy   epsl = %e\n",epsl);
	prob++;
      }
      
    } else {
      
      // set atmosphere 
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
   
    // set primitive values
    Mrho[ijk]  = rho;
    Mepsl[ijk] = epsl;
    Mp[ijk]    = pres;
    Mvx[ijk]   = conSix/(conD+pres+conT);
    Mvy[ijk]   = conSiy/(conD+pres+conT);
    Mvz[ijk]   = conSiz/(conD+pres+conT);
    
    grhd_compute_v2_pt(Mvx[ijk],Mvy[ijk],Mvz[ijk], gxx[ijk],
		       gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
		       &tmp1, &tmp2, &tmp3, &( Mv2[ijk]));

    CheckForNANandINF(11,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk]);
    
    if (prob) {
      printf("problem occured at (cons2prim ColdEoS_ColdStaticATM)\n");
      printf("l=%d  x=%e y=%e x=%e  ijk=%d  detg=%e\n",level->l, Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],ijk,detg[ijk]);
      if (ijkinsidefinerlevel(box,ijk)==0) {
	printf("  point is NOT inside finer box NOR in some symmetry area\n"); 
	printf("%s",errorstring);
      } else printf("  point is inside finer box/ in symmetry \n"); 

      /* 
      MD[ijk]=MSx[ijk]=MSy[ijk]=MSz[ijk]=MT[ijk]=Mrho[ijk]=Mepsl[ijk]=Mp[ijk]=Mvx[ijk]=Mvy[ijk]=Mvz[ijk] = 0.;
      */

      /*

      // 27/08/2012 
      // try set atm 
      Mrho[ijk]  = GRHD.ATM_RHOATM;
      Mepsl[ijk] = GRHD.ATM_EPSLATM;
      Mp[ijk]    = GRHD.ATM_PATM;
      Mvx[ijk]   = 0.;
      Mvy[ijk]   = 0.;
      Mvz[ijk]   = 0.;

      MD[ijk]  = SQRTdetg* GRHD.ATM_RHOATM;
      MSx[ijk] = 0.;
      MSy[ijk] = 0.;
      MSz[ijk] = 0.;
      MT[ijk]  = SQRTdetg* GRHD.ATM_RHOATM*GRHD.ATM_EPSLATM;

      */

      // 29/08/2012 
      // do nothing (keep old prim) // not ci 

      sprintf(errorstring,"Problem in c2p.");
      
    }
    
  } endfor_ijk_openmp;
  
  bampi_openmp_stop
}





// grhd_c2p_proot_ColdStaticAtm
//   - c2p for general EoS, cold-static atm treatment   
//   - Root of P with Newthon-Raphson method
void grhd_c2p_proot_ColdStaticAtm(tVarList *ucur, tVarList *primVars)
{

  if (PR) printf("# ===> grhd_c2p_proot_ColdStaticAtm()\n");

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
  double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg;
  double epsl,rho, p,pold,pnew,peos, kappa,chi, cs2;
  double temp,tempSQRTarg,tempSQRT,otemp, Wlor, f,drhodp,depsldp,df;
  double vx,vy,vz, v2, hrho, tmp1,tmp2,tmp3;
  
  int count,problem,atm;
  double error,r;
  
  const double errormax = Getd("grhd_C2P_NewtonRaphsonTR");
  const int countmax    = Getd("grhd_C2P_NewtonRaphsonNR");
  const double pmin = GRHD.ATM_PATM;
  
  
  forallpoints_ijk_openmp(level) {

    // undensitize conservatives
    detg[ijk] = invg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
		     &gixx,&gixy,&gixz,&giyy,&giyz,&gizz);
    if ((detg[ijk]<=0.) || (!finite(detg[ijk]))) {
      printf("%e %e %e   %d   %e\n", 
	     Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],ijk,detg[ijk]);
      errorexit("problem with determinant");
    }
    
    SQRTdetg = sqrt(detg[ijk]);
    ooSQRTdetg = 1./SQRTdetg;
    
    conD   = MD[ijk] *ooSQRTdetg;
    conSx  = MSx[ijk]*ooSQRTdetg;
    conSy  = MSy[ijk]*ooSQRTdetg;
    conSz  = MSz[ijk]*ooSQRTdetg;
    conT   = MT[ijk] *ooSQRTdetg;
    conSix = gixx*conSx + gixy*conSy + gixz*conSz;
    conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
    conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
    conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
    SQRTconS2 = sqrt(conS2);
    
    
    // check if cons vars have physical values
    if (CheckForNANandINF(6,conD,conSx,conSy,conSz,conT, SQRTconS2)) 
      errorexit("conD, conS[xyz], conS_i^i or conT is nan");
    
    if ((conD <= GRHD.ATM_RHOATM) || (conT < GRHD.ATM_EPSLATM) 
	|| (((conT+conD+GRHD.ATM_PATM)*(conT+conD+GRHD.ATM_PATM)*GRHD.HRSC_VMAX)<conS2) ) {
      
      // pt evolved in atm, 		
      // re-set atmosphere and continue 
      MD[ijk]  = SQRTdetg*GRHD.ATM_RHOATM;
      MSx[ijk] = 0.;
      MSy[ijk] = 0.;
      MSz[ijk] = 0.;
      MT[ijk]  = SQRTdetg*GRHD.ATM_RHOATM*GRHD.ATM_EPSLATM;
      
      Mrho[ijk]  = GRHD.ATM_RHOATM;
      Mepsl[ijk] = GRHD.ATM_EPSLATM;
      Mp[ijk]    = GRHD.ATM_PATM;
      Mvx[ijk]   = 0.;
      Mvy[ijk]   = 0.;
      Mvz[ijk]   = 0.;
      Mv2[ijk]   = 0.;
      
      continue;
      
    }
    
    // try inversion
    atm = 0;
    
    // guess pold from previous step
    pold = DMAX(Mp[ijk],pmin);
    
    if (!finite(pold)) errorexit("starting guess for p is not finite");
    
    // Newton-Raphson 
    count = 0;
    error = 1.;
    
    compute_rho(pold);
    compute_epsl(pold);
    problem = EOS.comp("re","","","pc","re","", rho,epsl, &peos,&cs2,&chi,&kappa);
    f       = pold - peos; 
    pnew    = pold;
    
    while ( (!problem) && (error>errormax) && (count<countmax) ) {
      
      count++;
      
      drhodp  = conD*conS2/(tempSQRT *temp*temp);
      depsldp = pold*conS2/(conD*tempSQRTarg*tempSQRT);
      df      = 1. - chi*drhodp - kappa*depsldp;
      
      pold    = pnew;
      pnew    = DMAX(pold-f/df, pmin);
      error   = fabs(pnew-pold)/fabs(pnew);
      
      compute_rho(pnew);
      compute_epsl(pnew);
      
      problem = EOS.comp("re","","","pc","re","", rho,epsl, &peos,&cs2,&chi,&kappa);
      
      f       = pnew - peos; 
      
    }
    
    // velocity 
    otemp = 1./(conD+pnew+conT); 
    vx = conSix *otemp;
    vy = conSiy *otemp;
    vz = conSiz *otemp;
    grhd_compute_v2_pt(vx,vy,vz, gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
		       &tmp1, &tmp2, &tmp3, &v2);

    // look for unphysical vals and cure setting atm 
    if (problem) atm++;
    if (!finite(pnew) || !finite(rho) || !finite(epsl)) atm++;
    if (rho<GRHD.ATM_FATM*GRHD.ATM_RHOATM) atm++;
    if (epsl<=0.) atm++;
    if (v2>=GRHD.HRSC_VMAX) atm++;
    if (count>=countmax) atm++;
    
    if (atm) {
      
      // set atmosphere
      MD[ijk]  = SQRTdetg*GRHD.ATM_RHOATM;
      MSx[ijk] = 0.;
      MSy[ijk] = 0.;
      MSz[ijk] = 0.;
      MT[ijk]  = SQRTdetg*GRHD.ATM_RHOATM*GRHD.ATM_EPSLATM;
      
      Mrho[ijk]  = GRHD.ATM_RHOATM;
      Mepsl[ijk] = GRHD.ATM_EPSLATM;
      Mp[ijk]    = GRHD.ATM_PATM;
      Mvx[ijk]   = 0.;
      Mvy[ijk]   = 0.;
      Mvz[ijk]   = 0.;
      Mv2[ijk]   = 0.;
      
      continue;
      
    }
    
    
    // check no iters
    /*
      if (count>=countmax) {
      if (PR) {
      printf("%d  %e %e    %e %e %.16e\n",count,error, errormax, f/df, f,df);
      printf("  %e %e %e => r=%e \n", Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],
      sqrt(pow(Ptr(level,"x")[ijk],2)+pow(Ptr(level,"y")[ijk],2)+pow(Ptr(level,"z")[ijk],2)));
      printf("  %+e %+e %+e (values  rho,epsl,p)\n",rho,epsl,pnew);
      printf("  %+e , %+e %+e %+e (lev, atm)\n",GRHD.ATM_FATM*GRHD.ATM_RHOATM,GRHD.ATM_RHOATM,GRHD.ATM_EPSLATM,GRHD.ATM_PATM);
      printf("  %+e %+e %+e (D,S2,T)\n",conD,conS2,conT);
      printf("  %+e %+e %+e %+e (vx,vy,vz,v)\n",
      conSix/(conD+pnew+conT),conSix/(conD+pnew+conT),conSix/(conD+pnew+conT),
      sqrt(conSix*conSix+conSiy*conSiy+conSiz*conSiz)/(conD+pnew+conT));
      }
      
      errorexit("Newton-Raphson did not converged");	     
      }
    */
    
    // set primitive values
    Mrho[ijk]  = rho;
    Mepsl[ijk] = epsl;
    Mp[ijk]    = pnew;
    Mvx[ijk]   = vx; 
    Mvy[ijk]   = vy; 
    Mvz[ijk]   = vz; 
    Mv2[ijk]   = v2; 
    
    // check for nans
    CheckForNANandINF(11,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], 
		      Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk]);
    
  } endfor_ijk_openmp;
  
  bampi_openmp_stop
}





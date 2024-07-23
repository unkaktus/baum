/* grhd_c2p_hroot.c 
   sbernuz 05/2012 */

#include "bam.h"
#include "grhd.h"


#define PR 0
#define DEBUG 0
#define ln10 2.302585092994046e+00


/* ************************************ 
   macros and functions */


// f(H) = 0 
double foH(double H, 
	   double conS2, double conE, double conD) {
  double W2, rho, epsl, peos; 
  double f;
  W2   = 0.5*(1.+sqrt(1.+4.*conS2/(H*H)));
  rho  = conD/sqrt(W2);
  epsl = (conE - H*(W2-1.))/rho - 1.;
  
  //problem = GRHD.use_eos(&peos, &tmp1, &tmp2, &tmp3, rho,epsl);	
  EOS.comp("re","","","p","","", rho,epsl, &peos);
  
  f = conE - H*W2 + peos;
  return f;
}


// Brent method (from NR)
double zbrent(double x1, double x2, 
	      double tol, const int itmax, 
	      const double conS2, const double conE, const double conD, 
	      int *error)
{

  const double EPS = 1e-14;

  *error = 0;

  int iter;
  double a=x1;
  double b=x2; 
  double c=x2;
  
  double fa=foH(a, conS2,conE,conD);
  double fb=foH(b, conS2,conE,conD);
  
  double d,e,min1,min2, fc,p,q,r,s,tol1,xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    //errorexit("root must be bracketed in zbrent");
    //printf("root must be bracketed in zbrent (%e,%e) -> (%e,%e)\n",a,b,fa,fb);
    *error=1; 
    //return 0.;
    return x1;
  }
  
  fc=fb;
  for (iter=1;iter<=itmax;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=foH(b, conS2,conE,conD);
  }
  //errorexit("maximum number of iterations exceeded in zbrent");
  *error=2;
  //return 0.;
  return x1;
}


/* ************************************ 
   c2p routines */


// matter_c2p_hroot_notatm
//   - c2p for general EoS
//   - Root of H = rho h with Brent method
void grhd_c2p_hroot(tVarList *ucur, tVarList *primVars)
{

  if (PR) printf("# ===> grhd_c2p_hroot()\n");

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
  double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg, otemp;
  double conE, epsl,rho, peos, tmp1, tmp2, tmp3;
  double vx,vy,vz, v2, Wlor2;
  double H, Hprev;

  const double Hmin  = 1e-18;  
  int c, cmax=12;
  int error, problem,atm;
    
  const double tol   = Getd("grhd_C2P_NewtonRaphsonTR");
  const int countmax = Getd("grhd_C2P_NewtonRaphsonNR");

  
  forallpoints_ijk(level) {
    
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
      
      // check if cons vars have physical values
      if (CheckForNANandINF(6,conD,conSx,conSy,conSz,conT, SQRTconS2)) 
	errorexit("conD, conS[xyz], conS_i^i or conT is nan");
      
      // try inversion
      error = 0;
      atm   = 0;
      
      conSix = gixx*conSx + gixy*conSy + gixz*conSz;
      conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
      conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
      conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
      SQRTconS2 = sqrt(conS2);
      conE   = conD+conT;

      // try to bracket the root with the prev val
      Hprev = Mp[ijk] + Mrho[ijk]*(Mepsl[ijk] + 1.);
      
      // Brent method
      for (c=0; c<cmax;  c++) {	  
	Hprev *= 2.;
	H = zbrent(Hmin,Hprev, tol,countmax, conS2,conE,conD, &error);
	if (!(error==1)) break;
      }      
      
      if (PR) if (error) {

	  printf(" PR=%d Detected problem in c2p, error = %d \n",error,PR);
	  printf("  %e %e %e => r=%e \n", 
		 Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],
		 sqrt(pow(Ptr(level,"x")[ijk],2)+pow(Ptr(level,"y")[ijk],2)+pow(Ptr(level,"z")[ijk],2)));
	  printf("  %+e %+e %+e %+e     (rho,epsl,p,H)\n",rho,epsl,peos,H);
	  printf("  %+e %+e %+e %+e     (D,S2,T,E)\n",conD,conS2,conT,conE);
	  printf("  %+e %+e %+e %+e     (vx,vy,vz,v)\n",
		 conSix/(conD+peos+conT),
		 conSix/(conD+peos+conT),
		 conSix/(conD+peos+conT),
		 sqrt(conSix*conSix+conSiy*conSiy+conSiz*conSiz)/(conD+peos+conT));

      }
      
      // checks
      if (error==1) errorexit("root must be bracketed in zbrent");
      if (error==2) errorexit("maximum number of iterations exceeded in zbrent");
      
      // primitives
      Wlor2   = 0.5*(1.+sqrt(1.+4.*conS2/(H*H)));
      rho     = conD/sqrt(Wlor2);
      epsl    = (conE - H*(Wlor2-1.))/rho - 1.;
      //problem = GRHD.use_eos(&peos, &cs2, &kappa, &chi, rho,epsl);
      EOS.comp("re","","","p","","", rho,epsl, &peos);
      
      
      // velocity 
      otemp = 1./(conD+peos+conT); 
      vx = conSix *otemp;
      vy = conSiy *otemp;
      vz = conSiz *otemp;
      grhd_compute_v2_pt(vx,vy,vz, gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
			 &tmp1, &tmp2, &tmp3, &v2);
      
      // look for unphysical vals and cure setting atm       
      if (!finite(peos) || !finite(rho) || !finite(epsl)) errorexit("unphysical values");
      if (v2>=GRHD.HRSC_VMAX) errorexit("velocity too high");
            
      // set primitive values
      Mrho[ijk]  = rho;
      Mepsl[ijk] = epsl;
      Mp[ijk]    = peos;
      Mvx[ijk]   = vx; 
      Mvy[ijk]   = vy; 
      Mvz[ijk]   = vz; 
      Mv2[ijk]   = v2; 
      
      // check for nans
      CheckForNANandINF(11,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], 
			Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk]);
      
  } endfor_ijk;
  
}





// grhd_c2p_hroot_ColdStaticATM 
//   - c2p for general EoS, cold-static atm treatment   
//   - Root of H = rho h with Brent method
void grhd_c2p_hroot_ColdStaticAtm(tVarList *ucur, tVarList *primVars)
{

  if (PR) printf("# ===> grhd_c2p_hroot_ColdStaticAtm()\n");

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
  double gixx,gixy,gixz,giyy,giyz,gizz, SQRTdetg,ooSQRTdetg, otemp;
  double conE, epsl,rho, peos, tmp1,tmp2,tmp3;
  double vx,vy,vz, v2, Wlor2;
  
  double H, Ha,Hb;
  int c, cmax=12;
  int error, problem,atm;
  
  const double Hatm  = GRHD.ATM_RHOATM*(1. + GRHD.ATM_EPSLATM) + GRHD.ATM_PATM;
  const double ooatmlev = 1./Getd("grhd_atm_level");
  const double tol   = Getd("grhd_C2P_NewtonRaphsonTR");
  const int countmax = Getd("grhd_C2P_NewtonRaphsonNR");
  
  forallpoints_ijk(level) {

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
    
    
    // check if cons vars have physical values
    if (CheckForNANandINF(6,conD,conSx,conSy,conSz,conT, SQRTconS2)) 
      errorexit("conD, conS[xyz], conS_i^i or conT is nan");
    
    
    if (conD < GRHD.ATM_FATM*GRHD.ATM_RHOATM) { 
      
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
    error = 0;
    atm   = 0;
    
    conSix = gixx*conSx + gixy*conSy + gixz*conSz;
    conSiy = gixy*conSx + giyy*conSy + giyz*conSz;
    conSiz = gixz*conSx + giyz*conSy + gizz*conSz;
    conS2  = conSx*conSix + conSy*conSiy + conSz*conSiz;
    SQRTconS2 = sqrt(conS2);
    conE   = conD+conT;
    
    // try to bracket the root with the prev val
    Hb = Mp[ijk] + Mrho[ijk]*(Mepsl[ijk] + 1.);
    Ha = Hatm;
    if (Hb==Hatm) Hb = ooatmlev;
    
    // Brent method
    for (c=0; c<cmax;  c++) {	  
      H = zbrent(Ha,Hb, tol,countmax, conS2,conE,conD, &error);
      Hb *= 2.;
      if (!(error==1)) break;
    }      
    
    // rem here: if (error) => H = Ha !
    
    if (PR) if (error) {
	
	printf(" PR=%d Detected problem in c2p, error = %d \n",error,PR);
	printf("  %e %e %e => r=%e \n", 
	       Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],
	       sqrt(pow(Ptr(level,"x")[ijk],2)+pow(Ptr(level,"y")[ijk],2)+pow(Ptr(level,"z")[ijk],2)));
	printf("  %+e %+e %+e %+e     (rho,epsl,p,H)\n",rho,epsl,peos,H);
	printf("  %+e %+e %+e %+e %+e (atm_lev,rhoatm,p,epsl,H)\n",GRHD.ATM_FATM*GRHD.ATM_RHOATM,GRHD.ATM_RHOATM,GRHD.ATM_EPSLATM,GRHD.ATM_PATM, Hatm);
	printf("  %+e %+e %+e %+e     (D,S2,T,E)\n",conD,conS2,conT,conE);
	printf("  %+e %+e %+e %+e     (vx,vy,vz,v)\n",
	       conSix/(conD+peos+conT),
	       conSix/(conD+peos+conT),
	       conSix/(conD+peos+conT),
	       sqrt(conSix*conSix+conSiy*conSiy+conSiz*conSiz)/(conD+peos+conT));
	
      }
    
    // primitives      
    Wlor2   = 0.5*(1.+sqrt(1.+4.*conS2/(H*H)));
    rho     = conD/sqrt(Wlor2);
    epsl    = (conE - H*(Wlor2-1.))/rho - 1.;
    //problem = GRHD.use_eos(&peos, &cs2, &kappa, &chi, rho,epsl);
    EOS.comp("re","","","p","","", rho,epsl, &peos);
    
    // velocity 
    otemp = 1./(conD+peos+conT); 
    vx = conSix *otemp;
    vy = conSiy *otemp;
    vz = conSiz *otemp;
    grhd_compute_v2_pt(vx,vy,vz, gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
		       &tmp1, &tmp2, &tmp3, &v2);
    
    // look for unphysical vals and cure setting atm       
    if (!finite(peos) || !finite(rho) || !finite(epsl)) atm++;      
    if (rho<GRHD.ATM_FATM*GRHD.ATM_RHOATM) atm++;                                 
    //if (epsl<=0.) atm++;
    if (v2>=GRHD.HRSC_VMAX) atm++;

    // checks
    //if ((!atm) && (error==1)) errorexit("root must be bracketed in zbrent");
    if ((!atm) && (error==2)) errorexit("maximum number of iterations exceeded in zbrent");
    
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
    
    
    // set primitive values
    Mrho[ijk]  = rho;
    Mepsl[ijk] = epsl;
    Mp[ijk]    = peos;
    Mvx[ijk]   = vx; 
    Mvy[ijk]   = vy; 
    Mvz[ijk]   = vz; 
    Mv2[ijk]   = v2; 
    
    
    // check for nans
    CheckForNANandINF(11,MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk], 
		      Mrho[ijk],Mepsl[ijk],Mp[ijk],Mvx[ijk],Mvy[ijk],Mvz[ijk]);
      
  } endfor_ijk;
  
}





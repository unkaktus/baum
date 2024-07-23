/* grhd_Vac.c 
   routines for real vacuum
   Wolfgang Tichy 09/2015 */

#include "bam.h"
#include "grhd.h"


/* compute max of pressure locally */
void grhd_atm_set_PMAX_local(tL *level)
{
  double presmax = 0.0;
  double *pres = Ptr(level, "grhd_p");
  int i;
  int size = bampi_size();
  double send[size], recv[size];

  /* find max pres, but only if pres is not a NULL pointer */
  //printf("level->l=%d pres=%p\n", level->l, pres);
  if(pres) forallpoints_ijk(level)
  {
    if (pres[ijk]>presmax) presmax = pres[ijk];
  } endfor_ijk;

  if (!finite(presmax)) errorexit(" presmax not finite!");

  /* bampi stuff */
  send[0] = presmax;
  bampi_alltoall(send,recv,1);
  for(i=0; i<size; i++)  if(recv[i]>presmax) presmax = recv[i];

  GRHD.ATM_PMAX = DMAX(presmax,GRHD.ATM_PMAX);
}

/* compute PMAX globally */
void grhd_atm_set_PMAX_global(tL *level)
{
  int l;
  tG *g = level->grid;

  GRHD.ATM_PMAX = 0.0;
  for(l = g->lmin; l <= g->lmax; l++)
    grhd_atm_set_PMAX_local(g->level[l]);
  if(!finite(GRHD.ATM_PMAX)) errorexit(" GRHD.ATM_PMAX is not finite!");
}



/* set vaccum mask, 0<= mask <=1 according to MD - MDVaclevel and nn */
void grhd_Vac_set_mask(tVarList *u, tVarList *primVars)
{
  tL *level = u->level;
  double *MD  = vldataptr(u, MATTER.INDX_VAR_q);
  //double *rho  = vldataptr(primVars, 0);
  double *detg = vldataptr(primVars, 7);
  int m;
  int pts = 1; // nn / direction
  double *mask = Ptr(level, "matter_mask");
  int Vac;
  double DVaclevel, MDVaclevel, SQRTdetg;

  DVaclevel = GRHD.ATM_RHOATM * GRHD.ATM_FATM * 0.5;
  // 0         ... matter = evolve anyway
  // 1/(pts+1) ... point is ATM, at least one bounding point is no ATM
  // 2/(pts+1) ... point and NN are ATM
  forallpoints_ijk(level)
  {
    double detgL = DMAX( detg[ijk] , DETGMIN );
    SQRTdetg = sqrt(detgL);
    MDVaclevel = SQRTdetg*DVaclevel;  /* MD[ijk] = SQRTdetg * consD; */

    if(MD[ijk]>MDVaclevel) Vac = 0;
    else 
    {
      Vac = 1;
      for (m=1; m<=pts; m++)
      {
	if ((i-m>=0)    && (MD[ijk-m*di]>MDVaclevel)) break;
	if ((j-m>=0)    && (MD[ijk-m*dj]>MDVaclevel)) break;
	if ((k-m>=0)    && (MD[ijk-m*dk]>MDVaclevel)) break;
	if ((i+m<=imax) && (MD[ijk+m*di]>MDVaclevel)) break;
	if ((j+m<=jmax) && (MD[ijk+m*dj]>MDVaclevel)) break;
	if ((k+m<=kmax) && (MD[ijk+m*dk]>MDVaclevel)) break;
	Vac++;
      }
    }
    mask[ijk] = (double)(Vac)/(double)(pts+1);
  } endfor_ijk;
}


/* set vacuum in evolved vars u or prims, if MD is below some level */
void grhd_Vac_set_u(tVarList *u, tVarList *primVars)
{
  int limit_grhd_D  = Getv("grhd_vacuum_limit_grhd_D","yes");
  double grhd_D_fac = Getd("grhd_vacuum_grhd_D_fac");
  double horizonalpha = Getd("grhd_vacuum_horizonalpha");
  int set_cons = Getv("grhd_vacuum_set","cons");
  int set_prim = Getv("grhd_vacuum_set","prim");
  tL *level = u->level;
  double *MD  = vldataptr(u, MATTER.INDX_VAR_q    );
  double *MT  = vldataptr(u, MATTER.INDX_VAR_q + 1);
  double *MSx = vldataptr(u, MATTER.INDX_VAR_q + 2);
  double *MSy = vldataptr(u, MATTER.INDX_VAR_q + 3);
  double *MSz = vldataptr(u, MATTER.INDX_VAR_q + 4);
  double *alpha = vldataptr(u, MATTER.INDX_VAR_alp);
  
  double *rho  = vldataptr(primVars, 0);
  double *epsl = vldataptr(primVars, 1);
  double *vx   = vldataptr(primVars, 2);
  double *vy   = vldataptr(primVars, 3);
  double *vz   = vldataptr(primVars, 4);
  double *pres = vldataptr(primVars, 5);
  double *v2   = vldataptr(primVars, 6);
  double *detg = vldataptr(primVars, 7);
  double DVaclevel, MDVaclevel, SQRTdetg;

  DVaclevel = GRHD.ATM_RHOATM * GRHD.ATM_FATM;

  /* loop over points */
  forallpoints_ijk(level)
  {
    double Wlor;
    double detgL = DMAX( detg[ijk] , DETGMIN );
    SQRTdetg = sqrt(detgL);
    MDVaclevel = SQRTdetg*DVaclevel;  /* MD[ijk] = SQRTdetg * consD; */

    /* compute Wlor from v^2 */
    if(v2[ijk]>=GRHD.HRSC_VMAX)  Wlor = GRHD.HRSC_WlorMAX;
    else                         Wlor = 1.0/sqrt(1-v2[ijk]);

    /* check which points we set to Vac */
    if( MD[ijk]<MDVaclevel || 
        (limit_grhd_D &&  alpha[ijk] < horizonalpha &&
         MD[ijk] > grhd_D_fac*SQRTdetg*Wlor*rho[ijk]) )
    {
      if(set_cons)
      {
        //detg[ijk] = DMAX( detg[ijk] , DETGMIN );
        //SQRTdetg = sqrt(detg[ijk]);
        MD[ijk]  = 0.0;
        MSx[ijk] = MSy[ijk] = MSz[ijk] = 0.0;
        MT[ijk]  = 0.0;
      }
      if(set_prim)
      {
        rho[ijk] = epsl[ijk] = pres[ijk] = 0.0;
        vx[ijk] = vy[ijk] = vz[ijk] = v2[ijk] = 0.0;
      }
    }
    
    if(CheckForNANandINF(6, MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk],detg[ijk])) {
      MD[ijk] = MSx[ijk] = MSy[ijk] = MSz[ijk] = MT[ijk] =0.;
      detg[ijk] = 1.;
      errorexit("NANs not cured by setting vaccum");
    }
  } endfor_ijk;
}


/* main routine for Vac settings 
   set Vac in previous time step (according to current primitives)
   (if required) set mask in current step */
void grhd_Vac_set_prerhs(tVarList *unew, tVarList *upre, double c,
                         tVarList *ucur, tVarList *primVars)
{
  int set_upre = Getv("grhd_vacuum_prerhs_set","upre");
  int set_ucur = Getv("grhd_vacuum_prerhs_set","ucur");
  tL *level = ucur->level;

  if(set_upre)
    grhd_Vac_set_u(upre, primVars);  /* set vac in prev time step */
  if(set_ucur)
    grhd_Vac_set_u(ucur, primVars);  /* set vac in curr time step */

  /* (if required) set mask in current step */
  if(MATTER.USEMASK)
    grhd_Vac_set_mask(upre, primVars);
}

/* set some atm consts, that may be needed later
   called during initialization (grhd_init(), grhd.c) */
void grhd_Vac_set_levels(tL *level)
{
  /* find the max of rho */
  grhd_atm_set_RHOMAX_global(level);
  grhd_atm_set_PMAX_global(level);

  /* grhd_atm_set_levels has been called before. 
     Should we change some pars to vac afterwards??? */
  if(Getd("grhd_atm_level")<=0.0)
  {
    GRHD.ATM_RHOATM = 0.0;
    GRHD.ATM_EPSLATM = GRHD.ATM_PATM = 0.0;
  }
  printf("grhd_Vac_set_levels:  GRHD.ATM_... vars are:\n");
  printf("RHOMAX=%g PMAX=%g RHOATM=%g FATM=%g EPSLATM=%g PATM=%g\n",
         GRHD.ATM_RHOMAX,GRHD.ATM_PMAX,GRHD.ATM_RHOATM,GRHD.ATM_FATM, 
         GRHD.ATM_EPSLATM,GRHD.ATM_PATM);
}

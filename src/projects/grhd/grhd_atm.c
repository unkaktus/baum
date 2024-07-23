/* grhd_atm.c 
   routines for cold static atmosphere based on rho floor
   sbernuz, mth 02/2012 */

#include "bam.h"
#include "grhd.h"


// compute rhomax locally
void grhd_atm_set_RHOMAX_local(tL *level)
{
  double rhomax = 0.;
  double *rho = Ptr(level, "grhd_rho");
  forallpoints_ijk(level) {
    if (rho[ijk]>rhomax) rhomax = rho[ijk];
  } endfor_ijk;
  
  if (!finite(rhomax)) errorexit(" rhomax not finite");
  
  int i;
  int size = bampi_size();
  double send[size], recv[size];
  
  send[0] = rhomax;
  bampi_alltoall(send,recv,1);
  
  for (i=0; i<size; i++) {
    if (recv[i]>rhomax) rhomax = recv[i];
  }
  /*dtim: We don't want to use always the GRHD.ATM_RHOMAX from the finest level, 
   therefore, compare rhomax with older values*/
  //GRHD.ATM_RHOMAX = rhomax;
  GRHD.ATM_RHOMAX = DMAX(rhomax,GRHD.ATM_RHOMAX);
}


// compute RHOMAX globally
void grhd_atm_set_RHOMAX_global(tL *level)
{
  int l;
  tG *g = level->grid;
  
  GRHD.ATM_RHOMAX = 0.;

  for (l = g->lmin; l <= g->lmax; l++){
    grhd_atm_set_RHOMAX_local(g->level[l]);
  }


  if (!finite(GRHD.ATM_RHOMAX)) errorexit(" GRHD.ATM_RHOMAX is nan");
}

// set atm levels
// called during initialization (grhd_init(), grhd.c) 
void grhd_atm_set_levels(tL *level)
{


  // atm levels
  grhd_atm_set_RHOMAX_global(level);

  GRHD.ATM_RHOATM  = DMAX(1e-20,GRHD.ATM_RHOMAX * Getd("grhd_atm_level"));
  GRHD.ATM_FATM    = Getd("grhd_atm_factor");

  if (0) printf("set ATM:  rhomax=%e rhoatm=%e fatm=%e epslatm=%e patm=%e\n",
		GRHD.ATM_RHOMAX,GRHD.ATM_RHOATM,GRHD.ATM_FATM, 
		GRHD.ATM_EPSLATM,GRHD.ATM_PATM);

  // 27.11.2013 SB
  // the following is wrong.
  // EPSATM = 0 so it cannot be used as input
  // e.g. for ideal it gives PATM=0
  /*
    EOS.comp("re","","", "p","","",
    GRHD.ATM_RHOATM,GRHD.ATM_EPSLATM, &(GRHD.ATM_PATM)); 
  */
  // the idea is to set the ATM according to the cold part of the EOS, ie
  EOS.comp("r","","","pe","","", 
	   GRHD.ATM_RHOATM, &(GRHD.ATM_PATM),&(GRHD.ATM_EPSLATM) );	  

  if (0) printf("set ATM:  rhomax=%e rhoatm=%e fatm=%e epslatm=%e patm=%e\n",
		GRHD.ATM_RHOMAX,GRHD.ATM_RHOATM,GRHD.ATM_FATM, 
		GRHD.ATM_EPSLATM,GRHD.ATM_PATM);
  
  //errorexit(""); 

}



// set primitives to atm at point ijk
// called during initialization (grhd.c) and c2p
void grhd_atm_set_prim_pt(double *p, double *rho, double *epsl, double *vx, double *vy, double *vz, double *v2)
{
  if ( *rho < GRHD.ATM_RHOATM*GRHD.ATM_FATM ) {
    *p    = GRHD.ATM_PATM;
    *rho  = GRHD.ATM_RHOATM;
    *epsl = GRHD.ATM_EPSLATM;
    *vx = *vy = *vz = *v2 = 0.;
  }
}


// set atmosphere maks, 0<= mask <=1 according to rho - atm level and nn
void grhd_atm_set_mask(tL *level)
{
  
  int m;
  int pts = 1; // nn / direction
  double *rho  = Ptr(level, "grhd_rho");
  double *mask = Ptr(level, "matter_mask");
  
  int atm;
  double rhoatmlevel = GRHD.ATM_RHOATM * GRHD.ATM_FATM;
  
  // 0         ... matter = evolve anyway
  // 1/(pts+1) ... point is ATM, at least one bounding point is no ATM
  // 2/(pts+1) ... point and NN are ATM
  forallpoints_ijk(level) {
    if (rho[ijk]>=rhoatmlevel) atm = 0;
    else {
      atm = 1;
      for (m=1; m<=pts; m++) {
	if ((i-m>=0)    && (rho[ijk-m*di]>=rhoatmlevel)) break;
	if ((j-m>=0)    && (rho[ijk-m*dj]>=rhoatmlevel)) break;
	if ((k-m>=0)    && (rho[ijk-m*dk]>=rhoatmlevel)) break;
	if ((i+m<=imax) && (rho[ijk+m*di]>=rhoatmlevel)) break;
	if ((j+m<=jmax) && (rho[ijk+m*dj]>=rhoatmlevel)) break;
	if ((k+m<=kmax) && (rho[ijk+m*dk]>=rhoatmlevel)) break;
	atm++;
      }
    }
    mask[ijk] = (double)(atm)/(double)(pts+1);
  } endfor_ijk;
  
}



// set atmosphere maks, bitmask 0/1 according to rho - atm level
void grhd_atm_set_bitmask(tL *level)
{
  
  int m;
  int pts = 1;
  double *rho  = Ptr(level, "grhd_rho");
  double *mask = Ptr(level, "matter_mask");
  
  int atm;
  double rhoatmlevel = GRHD.ATM_RHOATM * GRHD.ATM_FATM;
  
  // 0 ... matter = evolve anyway
  // 1 ... vacuum = do not evolve 
  forallpoints_ijk(level) {
    if (rho[ijk]>=rhoatmlevel) atm = 0;
    else {
      atm = 0;
      for (m=1; m<=pts; m++) {
	if ((i-m>=0)    && (rho[ijk-m*di]>=rhoatmlevel)) break;
	if ((j-m>=0)    && (rho[ijk-m*dj]>=rhoatmlevel)) break;
	if ((k-m>=0)    && (rho[ijk-m*dk]>=rhoatmlevel)) break;
	if ((i+m<=imax) && (rho[ijk+m*di]>=rhoatmlevel)) break;
	if ((j+m<=jmax) && (rho[ijk+m*dj]>=rhoatmlevel)) break;
	if ((k+m<=kmax) && (rho[ijk+m*dk]>=rhoatmlevel)) break;
	atm = 1;
      }
    }
    mask[ijk] = (double)(atm);
  } endfor_ijk;
  
}


// set rhs = 0 according to Euler timestep
void grhd_atm_set_r0tstep(tVarList *u, tVarList *r)
{
  
  tL *level = u->level;  
  double dt = level->dt;  

  double *MD  = vldataptr(u, MATTER.INDX_VAR_q    );
  double *MT  = vldataptr(u, MATTER.INDX_VAR_q + 1);
  
  double *RD  = vldataptr(r, MATTER.INDX_VAR_q    );
  double *RT  = vldataptr(r, MATTER.INDX_VAR_q + 1);
  double *RSx = vldataptr(r, MATTER.INDX_VAR_q + 2);
  double *RSy = vldataptr(r, MATTER.INDX_VAR_q + 3);
  double *RSz = vldataptr(r, MATTER.INDX_VAR_q + 4);
  
  double *mask = Ptr(level, "matter_mask");
  
  double Deuler, Teuler;


  forallpoints_ijk(level) {
    
    if ( (MATTER.USEMASK) && (mask[ijk]>0.9) ) continue;
    // consider only pts with rhs != 0
    
    // Euler timestep
    Deuler = MD[ijk] + dt*RD[ijk];
    Teuler = MT[ijk] + dt*RT[ijk];
    
    if ( (Deuler < 0.) ||
	 (Teuler < 0.) ) {
      
      RD[ijk] = RT[ijk] = RSx[ijk] = RSy[ijk] = RSz[ijk] = 0.;
      
    }
    
  } endfor_ijk;
  
}


// set atmosphere in evolved vars u
void grhd_atm_set_u(tVarList *u, tVarList *primVars)
{
  
  tL *level = u->level;
  
  double *MD  = vldataptr(u, MATTER.INDX_VAR_q    );
  double *MT  = vldataptr(u, MATTER.INDX_VAR_q + 1);
  double *MSx = vldataptr(u, MATTER.INDX_VAR_q + 2);
  double *MSy = vldataptr(u, MATTER.INDX_VAR_q + 3);
  double *MSz = vldataptr(u, MATTER.INDX_VAR_q + 4);
  
  double *rho  = vldataptr(primVars, 0);
  //double *detg = Ptr(level,"adm_detg");
  double *detg = vldataptr(primVars, 7);

  double SQRTdetg;

  double rhoatmlevel = GRHD.ATM_RHOATM * GRHD.ATM_FATM;
  double Datmlevel   = GRHD.ATM_RHOATM;
  double Tatmlevel   = GRHD.ATM_RHOATM * GRHD.ATM_EPSLATM;

  forallpoints_ijk(level) {
    
    if (rho[ijk]<rhoatmlevel) {
      
      detg[ijk] = DMAX( detg[ijk] , DETGMIN );
      SQRTdetg = sqrt(detg[ijk]);
      if (!finite(SQRTdetg)) SQRTdetg = 1.;
      
      MD[ijk]  = SQRTdetg*Datmlevel;
      MSx[ijk] = 0.;
      MSy[ijk] = 0.;
      MSz[ijk] = 0.;
      MT[ijk]  = SQRTdetg*Tatmlevel;
      
    }
    
    if (CheckForNANandINF(6, MD[ijk],MSx[ijk],MSy[ijk],MSz[ijk],MT[ijk],detg[ijk])) {
      MD[ijk] = MSx[ijk] = MSy[ijk] = MSz[ijk] = MT[ijk] =0.;
      detg[ijk] = 1.;
      errorexit(" nans not cured by setting atmosphere");
    }
    
  } endfor_ijk;
  
}


// main routine for atm settings 
// set atm in previous time step (according to current primitives)
// (if required) set mask in current step
void grhd_atm_set_prerhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars)
{
    
  // set atm in prev time step
  grhd_atm_set_u( upre, primVars);  
  
  // (if required) set mask in current step
  tL *level = ucur->level;
  if (MATTER.USEMASK) grhd_atm_set_mask(level);
  
}


// main routine for atm settings
// set rhs = 0 according to Euler time step
// (if required) set bitmask
void grhd_atm_set_postrhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars)
{

  // set atm in prev time step
  /* grhd_set_atm_u( upre, primVars); */

  // set rhs = 0 if atm in prev time step
  // (if required) set bitmask
  grhd_atm_set_r0tstep( ucur, unew);
  
}


void grhd_atm_set_postrhs_const(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars) {

  bampi_openmp_start

  tL *level = ucur->level;

  const double dx = level->dx;
  const double dy = level->dy;
  const double dz = level->dz;

  double *RD  = vldataptr(unew, MATTER.INDX_VAR_q    );
  double *RT  = vldataptr(unew, MATTER.INDX_VAR_q + 1);
  double *RSx = vldataptr(unew, MATTER.INDX_VAR_q + 2);
  double *RSy = vldataptr(unew, MATTER.INDX_VAR_q + 3);
  double *RSz = vldataptr(unew, MATTER.INDX_VAR_q + 4);

  double x,y,z;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  double boundary;

  forallpoints_ijk_openmp(level) {

    xmin = box->x0 + MATTER.HRSC_NGHOST*dx; 
    ymin = box->y0 + MATTER.HRSC_NGHOST*dy; 
    zmin = box->z0 + MATTER.HRSC_NGHOST*dz;
    xmax = box->x0+(box->m-1)*box->dx - MATTER.HRSC_NGHOST*dx;
    ymax = box->y0+(box->n-1)*box->dy - MATTER.HRSC_NGHOST*dy;
    zmax = box->z0+(box->o-1)*box->dz - MATTER.HRSC_NGHOST*dz;

    // printf("xmin, xmax = %+.6e, %+.6e \n", xmin, xmax);
    // printf("ymin, ymax = %+.6e, %+.6e \n", ymin, ymax);
    // printf("zmin, zmax = %+.6e, %+.6e \n", zmin, zmax);

    x = Ptr(level,"x")[ijk]; 
    y = Ptr(level,"y")[ijk]; 
    z = Ptr(level,"z")[ijk];

    if      ( (x <= xmin) ) boundary=1;
    else if ( (x >= xmax) ) boundary=1;

    else if ( (y <= ymin) ) boundary=1;
    else if ( (y >= ymax) ) boundary=1;

    else if ( (z <= zmin) ) boundary=1;
    else if ( (z >= zmax) ) boundary=1;

    else                    boundary=0;

    if (boundary) RD[ijk] = RT[ijk] = RSx[ijk] = RSy[ijk] = RSz[ijk] = 0.;

  } endfor_ijk_openmp

  bampi_openmp_stop

}
/* grhd_flx1d_fvpr.c 
   sbernuz 0/2012 */

#include "bam.h"
#include "grhd.h"


#define PR 0
#define CureFlxVmax 1

#define HO_TREAT_AVGS 0  // opts : 0, 1, 2
#define HO_TREAT_SPEED 0 // opts : 0, 1, 

#define ZERO 0.
#define HALF 0.5
#define ONE 1.0

// for debug
#define DEBUG 0
#define SetFlxToZero 0


// 2nd order finite volume fluxes based on Riemann solvers (HLL,LLF) 
// and primitive reconstruction
void grhd_flx1d_fvpr_turb(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d)
{

  // reconstructed vars at interfaces
  // rem L (R) i+1/2 => uplus "i" (umins "i+1") 
  
  double rhomins, epslmins, vxmins, vymins, vzmins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus;
  
  double tTauxxmins, tTauxymins, tTauxzmins, tTauyymins, tTauyzmins, tTauzzmins;
  double tTauxxplus, tTauxyplus, tTauxzplus, tTauyyplus, tTauyzplus, tTauzzplus;
  
  double alphaplus, betaplus, detgplus, gupxxplus, gupyyplus, gupzzplus, gupxyplus, gupxzplus;
  double alphamins, betamins, detgmins, gupxxmins, gupyymins, gupzzmins, gupxymins, gupxzmins;
  double gxxplus, gxyplus, gxzplus, gyyplus, gyzplus, gzzplus;
  double gxxmins, gxymins, gxzmins, gyymins, gyzmins, gzzmins;
  
  double pres_L, v2_L, Wlor_L, W2rhoh_L, cs2_L, vlowx_L, vlowy_L, vlowz_L;
  double pres_R, v2_R, Wlor_R, W2rhoh_R, cs2_R, vlowx_R, vlowy_R, vlowz_R;

  double tmp1,tmp2,tmp3,detginv;

  double *qL   = dmalloc(MATTER.NVq);
  double *qR   = dmalloc(MATTER.NVq);  
  double *fL   = dmalloc(MATTER.NVf);  
  double *fR   = dmalloc(MATTER.NVf);  
  double *f_i  = dmalloc(MATTER.NVf);  
  double *lamL = dmalloc(MATTER.NVf);
  double *lamR = dmalloc(MATTER.NVf);  

  int i,k;


  for( i = 0 ; i <= nx ; i++ ) {
    
    // if the point is atm, do nothing EXCEPT
    // set the fluxes to 0 to be sure no random numbers are added
    if ( (MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i]= ZERO;
      continue;
    }
    
    // left state, u_L =  u(i)^+, rec prim     
    alphaplus = MATTER.rec1dml( g1d[0], i );
    betaplus  = MATTER.rec1dml( g1d[1], i );
    gxxplus   = MATTER.rec1dml( g1d[2], i );
    gxyplus   = MATTER.rec1dml( g1d[3], i );
    gxzplus   = MATTER.rec1dml( g1d[4], i );
    gyyplus   = MATTER.rec1dml( g1d[5], i );
    gyzplus   = MATTER.rec1dml( g1d[6], i );
    gzzplus   = MATTER.rec1dml( g1d[7], i );

    rhoplus  = MATTER.rec1dl( w1d[0], i );
    epslplus = MATTER.rec1dl( w1d[1], i );
    vxplus   = MATTER.rec1dl( w1d[2], i );
    vyplus   = MATTER.rec1dl( w1d[3], i );
    vzplus   = MATTER.rec1dl( w1d[4], i );
    
    tTauxxplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau  ], i );
    tTauxyplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+1], i );
    tTauxzplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+2], i );
    tTauyyplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+3], i );
    tTauyzplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+4], i );
    tTauzzplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+5], i );

    // right state, u_R =  u(i+1)^-, rec prim 
    alphamins = MATTER.rec1dmr( g1d[0], i+1 );
    betamins  = MATTER.rec1dmr( g1d[1], i+1 );
    gxxmins   = MATTER.rec1dmr( g1d[2], i+1 );
    gxymins   = MATTER.rec1dmr( g1d[3], i+1 );
    gxzmins   = MATTER.rec1dmr( g1d[4], i+1 );
    gyymins   = MATTER.rec1dmr( g1d[5], i+1 );
    gyzmins   = MATTER.rec1dmr( g1d[6], i+1 );
    gzzmins   = MATTER.rec1dmr( g1d[7], i+1 );

    rhomins  = MATTER.rec1dr( w1d[0], i+1 );
    epslmins = MATTER.rec1dr( w1d[1], i+1 );
    vxmins   = MATTER.rec1dr( w1d[2], i+1 );
    vymins   = MATTER.rec1dr( w1d[3], i+1 );
    vzmins   = MATTER.rec1dr( w1d[4], i+1 );
    
    tTauxxmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau  ], i+1 );
    tTauxymins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+1], i+1 );
    tTauxzmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+2], i+1 );
    tTauyymins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+3], i+1 );
    tTauyzmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+4], i+1 );
    tTauzzmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+5], i+1 );

    // eos
    EOS.comp("re","","","pc","","", rhoplus,epslplus, &pres_L,&cs2_L);
    EOS.comp("re","","","pc","","", rhomins,epslmins, &pres_R,&cs2_R);
    
    // detg
    /*
    grhd_compute_detg_invg_pt(gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus,
    &detgplus, 
    &gupxxplus,&tmp1,&tmp2,&gupyyplus,&tmp3,&gupzzplus);

    grhd_compute_detg_invg_pt(gxxmins,gxymins,gxzmins,gyymins,gyzmins,gzzmins,
    &detgmins,
    &gupxxmins,&tmp1,&tmp2,&gupyymins,&tmp3,&gupzzmins );
    */
    // this saves time (equivalent for rec with symmetric stencils wrt i+1/2)
    detgplus  = detgmins  = detg(gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus);
    detginv   = 1./detgplus;
    gupxxplus = gupxxmins = detginv*(gyyplus*gzzplus - gyzplus*gyzplus);
    gupyyplus = gupyymins = detginv*(gxxplus*gzzplus - gxzplus*gxzplus);
    gupzzplus = gupzzmins = detginv*(gxxplus*gyyplus - gxyplus*gxyplus);
    gupxyplus = gupxzplus = detginv*(gyzplus*gxzplus - gzzplus*gxyplus);
    gupxzplus = gupxzmins = detginv*(gxyplus*gyzplus - gyyplus*gxzplus);

    // v2
    grhd_compute_v2_pt(vxplus,vyplus,vzplus, gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus,
		       &vlowx_L,&vlowy_L,&vlowz_L, &v2_L );

    grhd_compute_v2_pt(vxmins,vymins,vzmins, gxxmins,gxymins,gxzmins,gyymins,gyzmins,gzzmins,
                       &vlowx_R,&vlowy_R,&vlowz_R, &v2_R );

    v2_L = DMIN( v2_L, GRHD.HRSC_VMAX);
    v2_R = DMIN( v2_R, GRHD.HRSC_VMAX);
    
    if (v2_L==GRHD.HRSC_VMAX) Wlor_L = GRHD.HRSC_WlorMAX;
    else Wlor_L = ONE/sqrt( ONE - v2_L );
    
    if (v2_R>=GRHD.HRSC_VMAX) Wlor_R = GRHD.HRSC_WlorMAX;
    else Wlor_R = ONE/sqrt( ONE - v2_R );

#if (CureFlxVmax)
    // important hack used in bam 11                                                                                                                  
    if ( (Wlor_L == GRHD.HRSC_WlorMAX) || (Wlor_R == GRHD.HRSC_WlorMAX) ) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i] = ZERO;
      if (PR) printf(" v>VMAX on both sides, cured flx1d setting them to zero\n");
      continue;
    }
#endif
    
    // compute conservative vars
    grhd_compute_q_pt(gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus, detgplus, 
		      rhoplus,epslplus,pres_L, vlowx_L,vlowy_L,vlowz_L, Wlor_L,
		      &(qL[0]), &(qL[1]), &(qL[2]), &(qL[3]), &(qL[4]));
    
    grhd_compute_q_pt(gxxmins,gxymins,gxzmins,gyymins,gyzmins,gzzmins, detgmins,
                      rhomins,epslmins,pres_R, vlowx_R,vlowy_R,vlowz_R, Wlor_R,
                      &(qR[0]), &(qR[1]), &(qR[2]), &(qR[3]), &(qR[4]));

    // compute phys fluxes
    grhd_compute_fx_pt(alphaplus, betaplus, detgplus, 
		       pres_L, vxplus, 
		       qL[0], qL[1], qL[2], qL[3], qL[4],
		       &(fL[0]), &(fL[1]), &(fL[2]), &(fL[3]), &(fL[4]));
    
    grhd_compute_fx_pt(alphamins, betamins, detgmins,
                       pres_R, vxmins,
                       qR[0], qR[1], qR[2], qR[3], qR[4],
                       &(fR[0]), &(fR[1]), &(fR[2]), &(fR[3]), &(fR[4]));
                       
   grhd_addturb_fx_pt(tTauxxplus, tTauxyplus, tTauxzplus, tTauyyplus, tTauyzplus, tTauzzplus,
			alphaplus, detgplus, gupxxplus, gupxyplus, gupxzplus,
   			&(fL[2]), &(fL[3]), &(fL[4]));
   			
   grhd_addturb_fx_pt(tTauxxmins, tTauxymins, tTauxzmins, tTauyymins, tTauyzmins, tTauzzmins,
			alphamins, detgmins, gupxxmins, gupxymins, gupxzmins,
   			&(fR[2]), &(fR[3]), &(fR[4]));

    // eigenvalues lam^0, lam^-, lam^+
    grhd_eig(alphaplus, betaplus, vxplus,v2_L,gupxxplus, cs2_L,
	     &(lamL[1]),&(lamL[0]),&(lamL[4]));
    lamL[2] = lamL[3] = lamL[1];

    grhd_eig(alphamins, betamins, vxmins,v2_R,gupxxmins, cs2_R,
             &(lamR[1]),&(lamR[0]),&(lamR[4]));
    lamR[2] = lamR[3] = lamR[1];

    // numerical fluxes 

    MATTER.hrsc_riemsol1d(qL,qR, fL,fR, lamL,lamR, f_i, MATTER.NVf);
    for (k=0; k<MATTER.NVf; k++) flux1d[k][i]=f_i[k];
    
    if (EOS.COLD) 
      // set Tau fluxes = 0 
      // (sources are also = 0 and Tau is updated in c2p)
      flux1d[1][i]=ZERO;
      
    
#if (SetFlxToZero)
    // test
    for (k=0; k<MATTER.NVf; k++) flux1d[k][i]=ZERO;
#endif
    
    // info    
    if (PR) {
      
      
      // debug
      if (DEBUG) {
	printf(" ************************** i =%d \n",i);
	printf("  var \t F_i \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %d \t %e \n",k,flux1d[k][i]);
	
	printf(" LEFT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qL[k],fL[k],lamL[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhoplus,epslplus,vxplus,vyplus,vzplus,v2_L,pres_L);
	
	printf(" RIGHT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qR[k],fR[k],lamR[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhomins,epslmins,vxmins,vymins,vzmins,v2_R,pres_R);
	
      }
      
      if ((v2_L==GRHD.HRSC_VMAX) || (v2_R==GRHD.HRSC_VMAX)) {
	printf("               %e %e   (%d %d)\n", Wlor_L,Wlor_R,i,nx);
	printf(" flx1d max speed reached -> corrected L=%e R=%e -> %e\n", v2_L,v2_R,GRHD.HRSC_VMAX);	  
      }
      
      if (CheckForNANandINF(5, flux1d[0][i],flux1d[1][i],flux1d[2][i],flux1d[3][i],flux1d[4][i])) {
	
	printf(" nans inside flx1d\n");
	
	printf("  var \t F_i \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %d \t %e \n",k,flux1d[k][i]);
	
	printf(" LEFT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qL[k],fL[k],lamL[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhoplus,epslplus,vxplus,vyplus,vzplus,v2_L,pres_L);
	
	printf(" RIGHT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qR[k],fR[k],lamR[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhomins,epslmins,vxmins,vymins,vzmins,v2_R,pres_R);
	
	/* errorexit(" (PR=%d) I stop the code",PR); */
	
      }
      
    }
    
  }
    free(qL); free(qR); 
    free(fL); free(fR);
    free(f_i); free(lamL); free(lamR); 
}


// 2nd order finite volume fluxes based on Riemann solvers (HLL,LLF) 
// and primitive reconstruction
// 14.11.2013 in this version the quantity b^i := (Wlor v^i) is reconstructed instead of v^i
void grhd_flx1d_fvprb_turb(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d)
{

  // reconstructed vars at interfaces
  // rem L (R) i+1/2 => uplus "i" (umins "i+1") 
  
  double rhomins, epslmins, vxmins, vymins, vzmins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus;
    
  double tTauxxmins, tTauxymins, tTauxzmins, tTauyymins, tTauyzmins, tTauzzmins;
  double tTauxxplus, tTauxyplus, tTauxzplus, tTauyyplus, tTauyzplus, tTauzzplus;
  
  double alphaplus, betaplus, detgplus, gupxxplus, gupyyplus, gupzzplus, gupxyplus, gupxzplus;
  double alphamins, betamins, detgmins, gupxxmins, gupyymins, gupzzmins, gupxymins, gupxzmins;
  double gxxplus, gxyplus, gxzplus, gyyplus, gyzplus, gzzplus;
  double gxxmins, gxymins, gxzmins, gyymins, gyzmins, gzzmins;
  
  double pres_L, v2_L, Wlor_L, W2rhoh_L, cs2_L, vlowx_L, vlowy_L, vlowz_L;
  double pres_R, v2_R, Wlor_R, W2rhoh_R, cs2_R, vlowx_R, vlowy_R, vlowz_R;

  double tmp1,tmp2,tmp3,detginv;

  // physical conservative and fluxes at interfaces
  /*double qL[MATTER.NVq], qR[MATTER.NVq];
  double fL[MATTER.NVf], fR[MATTER.NVf], f_i[MATTER.NVf]; 
  double lamL[MATTER.NVf], lamR[MATTER.NVf]; */

  double *qL = dmalloc(MATTER.NVq);
  double *qR = dmalloc(MATTER.NVq);  
  double *fL = dmalloc(MATTER.NVf);  
  double *fR = dmalloc(MATTER.NVf);  
  double *f_i = dmalloc(MATTER.NVf);  
  double *lamL = dmalloc(MATTER.NVf);
  double *lamR = dmalloc(MATTER.NVf); 

  int i,k;

  // need to compute the quantities b^i
  const int size = nx + 2*ngx; 
  const int off  = ngx-1;
  double *buf[3];
  buf[0] = dmalloc(size);
  buf[1] = dmalloc(size);
  buf[2] = dmalloc(size);
  double *bx,*by,*bz; 
  double bxplus,byplus,bzplus;
  double bxmins,bymins,bzmins;
  double b2_L,b2_R, ooW_L,ooW_R;

  // b^i at phys pts
  //for( i = 0 ; i <= nx ; i++ ) {
  for( i = -ngx+1 ; i <= (nx+ngx) ; i++ ) {

    // use "plus/L" as tmp vars
    vxplus = w1d[2][i]; // v^i 
    vyplus = w1d[3][i];
    vzplus = w1d[4][i];
    

    gxxplus = g1d[2][i]; // g_{ij}
    gxyplus = g1d[3][i];
    gxzplus = g1d[4][i];
    gyyplus = g1d[5][i];
    gyzplus = g1d[6][i];
    gzzplus = g1d[7][i];

    grhd_compute_v2_pt(vxplus,vyplus,vzplus, gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus,
		       &tmp1,&tmp2,&tmp3, &v2_L );
    v2_L   = DMIN( v2_L, GRHD.HRSC_VMAX);
    Wlor_L = ONE/sqrt( ONE - v2_L );
    
    *(buf[0]+off +i) = Wlor_L * vxplus;
    *(buf[1]+off +i) = Wlor_L * vyplus;
    *(buf[2]+off +i) = Wlor_L * vzplus;

  }
  // extrapolate all the data to the ghost points
  // no, use extrapolated data in loop above
  /*
  for (i=-ngx+1; i<1; i++) { 
    buf[0][i+off] = buf[0][1+off];
    buf[1][i+off] = buf[1][1+off];
    buf[2][i+off] = buf[2][1+off];
  }

  for (i=nx+1; i<=(nx+ngx); i++) { 
    buf[0][i+off] = buf[0][nx+off];
    buf[1][i+off] = buf[1][nx+off];
    buf[2][i+off] = buf[2][nx+off];
  }
  */

  bx = buf[0] + off;
  by = buf[1] + off;
  bz = buf[2] + off;


  //if (PR) printf(" %d %d %d \n",MATTER.NVw,MATTER.NVq,MATTER.NVf);

  for( i = 0 ; i <= nx ; i++ ) {
    
    // if the point is atm, do nothing EXCEPT
    // set the fluxes to 0 to be sure no random numbers are added
    if ( (MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i]= ZERO;
      continue;
    }

    // left state, u_L =  u(i)^+, rec prim     
    alphaplus = MATTER.rec1dml( g1d[0], i );
    betaplus  = MATTER.rec1dml( g1d[1], i );
    gxxplus   = MATTER.rec1dml( g1d[2], i );
    gxyplus   = MATTER.rec1dml( g1d[3], i );
    gxzplus   = MATTER.rec1dml( g1d[4], i );
    gyyplus   = MATTER.rec1dml( g1d[5], i );
    gyzplus   = MATTER.rec1dml( g1d[6], i );
    gzzplus   = MATTER.rec1dml( g1d[7], i );

    rhoplus  = MATTER.rec1dl( w1d[0], i );
    epslplus = MATTER.rec1dl( w1d[1], i );
    bxplus   = MATTER.rec1dl( bx, i ); // note now this is bplus^i
    byplus   = MATTER.rec1dl( by, i );
    bzplus   = MATTER.rec1dl( bz, i );
    
    tTauxxplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau  ], i );
    tTauxyplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+1], i );
    tTauxzplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+2], i );
    tTauyyplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+3], i );
    tTauyzplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+4], i );
    tTauzzplus	= MATTER.rec1dl( o1d[GRHD.IND_turbTau+5], i );
    

    // right state, u_R =  u(i+1)^-, rec prim 
    alphamins = MATTER.rec1dmr( g1d[0], i+1 );
    betamins  = MATTER.rec1dmr( g1d[1], i+1 );
    gxxmins   = MATTER.rec1dmr( g1d[2], i+1 );
    gxymins   = MATTER.rec1dmr( g1d[3], i+1 );
    gxzmins   = MATTER.rec1dmr( g1d[4], i+1 );
    gyymins   = MATTER.rec1dmr( g1d[5], i+1 );
    gyzmins   = MATTER.rec1dmr( g1d[6], i+1 );
    gzzmins   = MATTER.rec1dmr( g1d[7], i+1 );

    rhomins  = MATTER.rec1dr( w1d[0], i+1 );
    epslmins = MATTER.rec1dr( w1d[1], i+1 );
    bxmins   = MATTER.rec1dr( bx, i+1 ); // note now this is bmins^i
    bymins   = MATTER.rec1dr( by, i+1 );
    bzmins   = MATTER.rec1dr( bz, i+1 );
    
    tTauxxmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau  ], i+1 );
    tTauxymins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+1], i+1 );
    tTauxzmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+2], i+1 );
    tTauyymins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+3], i+1 );
    tTauyzmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+4], i+1 );
    tTauzzmins	= MATTER.rec1dr( o1d[GRHD.IND_turbTau+5], i+1 );


    // eos
    EOS.comp("re","","","pc","","", rhoplus,epslplus, &pres_L,&cs2_L);
    EOS.comp("re","","","pc","","", rhomins,epslmins, &pres_R,&cs2_R);
    
    
    // detg
    detgplus  = detgmins  = detg(gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus);
    detginv   = 1./detgplus;
    gupxxplus = gupxxmins = detginv*(gyyplus*gzzplus - gyzplus*gyzplus);
    gupyyplus = gupyymins = detginv*(gxxplus*gzzplus - gxzplus*gxzplus);
    gupzzplus = gupzzmins = detginv*(gxxplus*gyyplus - gxyplus*gxyplus);
    gupxyplus = gupxzplus = detginv*(gyzplus*gxzplus - gzzplus*gxyplus);
    gupxzplus = gupxzmins = detginv*(gxyplus*gyzplus - gyyplus*gxzplus);
    

    // b2 
    // note now this is b^2, and v_i is actually b_i (corrected below)
    grhd_compute_v2_pt(bxplus,byplus,bzplus, gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus, &vlowx_L,&vlowy_L,&vlowz_L, &b2_L );
    grhd_compute_v2_pt(bxmins,bymins,bzmins, gxxmins,gxymins,gxzmins,gyymins,gyzmins,gzzmins, &vlowx_R,&vlowy_R,&vlowz_R, &b2_R );
    
    Wlor_L = sqrt(fabs(b2_L + ONE)); // note Wlor computation now, W^2 = b^2+1 !!!
    Wlor_R = sqrt(fabs(b2_R + ONE));

    ooW_L = (Wlor_L>ONE) ? ONE/Wlor_L : ONE; // optim 1/W , used later
    ooW_R = (Wlor_R>ONE) ? ONE/Wlor_R : ONE;

    vxplus = bxplus * ooW_L; // v^i = b^i / W (maybe not all components neded)
    vyplus = byplus * ooW_L;
    vzplus = bzplus * ooW_L;

    vxmins = bxmins * ooW_R;
    vymins = bymins * ooW_R;
    vzmins = bzmins * ooW_R;

    vlowx_L *= ooW_L; // b_i = W v_i (correct vars from above)
    vlowy_L *= ooW_L;
    vlowz_L *= ooW_L;

    vlowx_R *= ooW_R;
    vlowy_R *= ooW_R;
    vlowz_R *= ooW_R;

    v2_L = b2_L * ooW_L*ooW_L; // v^2 = b^2/W^2
    v2_R = b2_R * ooW_R*ooW_R;

    v2_L = DMIN( v2_L, GRHD.HRSC_VMAX); // needed ? testme !
    v2_R = DMIN( v2_R, GRHD.HRSC_VMAX);
    
#if (0) // (CureFlxVmax)
    // important hack used in bam 11                                                                                                                  
    // do we need it in this routine ? testme !    
    if ( (Wlor_L == GRHD.HRSC_WlorMAX) || (Wlor_R == GRHD.HRSC_WlorMAX) ) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i] = ZERO;      
      if (PR) printf(" v>VMAX on both sides, cured flx1d setting them to zero\n");
      continue;
    }    
#endif

    // compute conservative vars
    grhd_compute_q_pt(gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus, detgplus, 
		      rhoplus,epslplus,pres_L, vlowx_L,vlowy_L,vlowz_L, Wlor_L,
		      &(qL[0]), &(qL[1]), &(qL[2]), &(qL[3]), &(qL[4]));
    
    grhd_compute_q_pt(gxxmins,gxymins,gxzmins,gyymins,gyzmins,gzzmins, detgmins,
                      rhomins,epslmins,pres_R, vlowx_R,vlowy_R,vlowz_R, Wlor_R,
                      &(qR[0]), &(qR[1]), &(qR[2]), &(qR[3]), &(qR[4]));

    // compute phys fluxes
    grhd_compute_fx_pt(alphaplus, betaplus, detgplus, 
		       pres_L, vxplus, 
		       qL[0], qL[1], qL[2], qL[3], qL[4],
		       &(fL[0]), &(fL[1]), &(fL[2]), &(fL[3]), &(fL[4]));
    
    grhd_compute_fx_pt(alphamins, betamins, detgmins,
                       pres_R, vxmins,
                       qR[0], qR[1], qR[2], qR[3], qR[4],
                       &(fR[0]), &(fR[1]), &(fR[2]), &(fR[3]), &(fR[4]));
                       
    grhd_addturb_fx_pt(	tTauxxplus, tTauxyplus, tTauxzplus, tTauyyplus, tTauyzplus, tTauzzplus,
			alphaplus, detgplus, gupxxplus, gupxyplus, gupxzplus,
   			&(fL[2]), &(fL[3]), &(fL[4]));
   			
   grhd_addturb_fx_pt(	tTauxxmins, tTauxymins, tTauxzmins, tTauyymins, tTauyzmins, tTauzzmins,
			alphamins, detgmins, gupxxmins, gupxymins, gupxzmins,
   			&(fR[2]), &(fR[3]), &(fR[4]));
                       

    // eigenvalues lam^0, lam^-, lam^+
    grhd_eig(alphaplus, betaplus, vxplus,v2_L,gupxxplus, cs2_L,
	     &(lamL[1]),&(lamL[0]),&(lamL[4]));
    lamL[2] = lamL[3] = lamL[1];

    grhd_eig(alphamins, betamins, vxmins,v2_R,gupxxmins, cs2_R,
             &(lamR[1]),&(lamR[0]),&(lamR[4]));
    lamR[2] = lamR[3] = lamR[1];

    // numerical fluxes 
    MATTER.hrsc_riemsol1d(qL,qR, fL,fR, lamL,lamR, f_i, MATTER.NVf);
    for (k=0; k<MATTER.NVf; k++) flux1d[k][i]=f_i[k];

    if (EOS.COLD) 
      // set Tau fluxes = 0 
      // (sources are also = 0 and Tau is updated in c2p)
      flux1d[1][i]=ZERO;
         
#if (SetFlxToZero)
    // test
    for (k=0; k<MATTER.NVf; k++) flux1d[k][i]=ZERO;
#endif
    
    // info    
    if (PR) {
      
      
      // debug
      if (DEBUG) {
	printf(" ************************** i =%d \n",i);
	printf("  var \t F_i \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %d \t %e \n",k,flux1d[k][i]);
	
	printf(" LEFT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qL[k],fL[k],lamL[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhoplus,epslplus,vxplus,vyplus,vzplus,v2_L,pres_L);
	
	printf(" RIGHT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qR[k],fR[k],lamR[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhomins,epslmins,vxmins,vymins,vzmins,v2_R,pres_R);
	
      }
      
      if ((v2_L==GRHD.HRSC_VMAX) || (v2_R==GRHD.HRSC_VMAX)) {
	printf("               %e %e   (%d %d)\n", Wlor_L,Wlor_R,i,nx);
	printf(" flx1d max speed reached -> corrected L=%e R=%e -> %e\n", v2_L,v2_R,GRHD.HRSC_VMAX);	  
      }
      
      if (CheckForNANandINF(5, flux1d[0][i],flux1d[1][i],flux1d[2][i],flux1d[3][i],flux1d[4][i])) {
	
	printf(" nans inside flx1d\n");
	
	printf("  var \t F_i \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %d \t %e \n",k,flux1d[k][i]);
	
	printf(" LEFT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qL[k],fL[k],lamL[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhoplus,epslplus,vxplus,vyplus,vzplus,v2_L,pres_L);
	
	printf(" RIGHT\n");
	printf(" q \t f \t lam \n");
	for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qR[k],fR[k],lamR[k]);
	printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e\n", rhomins,epslmins,vxmins,vymins,vzmins,v2_R,pres_R);
	
	/* errorexit(" (PR=%d) I stop the code",PR); */
	
      }
      
    }
    
  }

  // extra memory for b^i (1d)
  free(buf[0]);
  free(buf[1]);
  free(buf[2]);
  free(qL); free(qR); 
  free(fL); free(fR);
  free(f_i); free(lamL); free(lamR); 
}




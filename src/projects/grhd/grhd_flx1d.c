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
void grhd_flx1d_fvpr(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d)
{

  // reconstructed vars at interfaces
  // rem L (R) i+1/2 => uplus "i" (umins "i+1") 
  
  double rhomins, epslmins, vxmins, vymins, vzmins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus;
  
  double alphaplus, betaplus, detgplus, gupxxplus, gupyyplus, gupzzplus;
  double alphamins, betamins, detgmins, gupxxmins, gupyymins, gupzzmins;
  double gxxplus, gxyplus, gxzplus, gyyplus, gyzplus, gzzplus;
  double gxxmins, gxymins, gxzmins, gyymins, gyzmins, gzzmins;
  
  double pres_L, v2_L, Wlor_L, W2rhoh_L, cs2_L, vlowx_L, vlowy_L, vlowz_L;
  double pres_R, v2_R, Wlor_R, W2rhoh_R, cs2_R, vlowx_R, vlowy_R, vlowz_R;

  double tmp1,tmp2,tmp3,detginv;

  // physical conservative and fluxes at interfaces
/*  double qL[MATTER.NVq], qR[MATTER.NVq];
  double fL[MATTER.NVf], fR[MATTER.NVf], f_i[MATTER.NVf]; 
  double lamL[MATTER.NVf], lamR[MATTER.NVf]; 
*/
  double *qL   = dmalloc(MATTER.NVq);
  double *qR   = dmalloc(MATTER.NVq);  
  double *fL   = dmalloc(MATTER.NVf);  
  double *fR   = dmalloc(MATTER.NVf);  
  double *f_i  = dmalloc(MATTER.NVf);  
  double *lamL = dmalloc(MATTER.NVf);
  double *lamR = dmalloc(MATTER.NVf); 

  int i,k;


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
    vxplus   = MATTER.rec1dl( w1d[2], i );
    vyplus   = MATTER.rec1dl( w1d[3], i );
    vzplus   = MATTER.rec1dl( w1d[4], i );

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
void grhd_flx1d_fvprb(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d)
{

  // reconstructed vars at interfaces
  // rem L (R) i+1/2 => uplus "i" (umins "i+1") 
  
  double rhomins, epslmins, vxmins, vymins, vzmins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus;
  
  double alphaplus, betaplus, detgplus, gupxxplus, gupyyplus, gupzzplus;
  double alphamins, betamins, detgmins, gupxxmins, gupyymins, gupzzmins;
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

    // eos
    EOS.comp("re","","","pc","","", rhoplus,epslplus, &pres_L,&cs2_L);
    EOS.comp("re","","","pc","","", rhomins,epslmins, &pres_R,&cs2_R);
    
    // detg
    detgplus  = detgmins  = detg(gxxplus,gxyplus,gxzplus,gyyplus,gyzplus,gzzplus);
    detginv   = 1./detgplus;
    gupxxplus = gupxxmins = detginv*(gyyplus*gzzplus - gyzplus*gyzplus);
    gupyyplus = gupyymins = detginv*(gxxplus*gzzplus - gxzplus*gxzplus);
    gupzzplus = gupzzmins = detginv*(gxxplus*gyyplus - gxyplus*gxyplus);

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


// high-order llf-weno fluxes computed from characteristic vars
void grhd_flx1d_ho_recchange(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d)
{
  
  int i,l; // -ngx+1 ... nx+ngx / 0 ... nx
  int j,k; //  0 ... ne
  
  const int ne     = MATTER.NVf; // no eigenv 
  const int ntot   = nx + 2*ngx;
  const int of     = ngx-1;      // offset
  
  const double disschi = 1.; // dissipation factor [1,1.3]
  const int S = MATTER.HRSC_NGHOST/2; // = 2 // rec stencil 


  // metric
  double *alpha = g1d[0];
  double *beta  = g1d[1]; // beta^(x)
  double *gxx   = g1d[2];
  double *gxy   = g1d[3];
  double *gxz   = g1d[4];
  double *gyy   = g1d[5];
  double *gyz   = g1d[6];
  double *gzz   = g1d[7];

  // conservatives
  double *conD  = q1d[0];
  double *conT  = q1d[1]; 
  double *conSx = q1d[2];
  double *conSy = q1d[3];
  double *conSz = q1d[4];

  // primitives
  double *rho  = w1d[0];
  double *epsl = w1d[1]; 
  double *vx   = w1d[2];
  double *vy   = w1d[3];
  double *vz   = w1d[4];
  double *pres = w1d[5];

  // ptwise metric and primitives at interfaces i+1/2
  double rho_avg, epsl_avg, vx_avg, vy_avg, vz_avg, pres_avg;
  double cs2_avg, h_avg, Wlor_avg; 
  double v2_avg, vlowx_avg, vlowy_avg, vlowz_avg;  
  double kappa_avg, chi_avg;
  double gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg;
  double alpha_avg, betax_avg, detg_avg, gupxx_avg,gupyy_avg,gupzz_avg;

  // ptwise metric and primitives at pt i
  double v2_i, vlowx_i, vlowy_i, vlowz_i, Wlor_i;    
  double alpha_i, beta_i, detg_i, sqrtg_i, detginv, gupxx_i,gupyy_i,gupzz_i;
  double pres_i, kappa_i, chi_i, cs2_i, avb_i;

#if (!(HO_TREAT_AVGS==0))

  // vars to be used for avgs
  double *detga = dmalloc(ntot); double *gupxx = dmalloc(ntot);
  double *vlowx = dmalloc(ntot); double *vlowy = dmalloc(ntot), 
  double *vlowz = dmalloc(ntot); double *v2 = dmalloc(ntot); double *Wlor = dmalloc(ntot);
  double *h = dmalloc(ntot); double *kappa = dmalloc(ntot); double *cs2 = dmalloc(ntot);
  
#endif

  // conservative and fluxes at pts i
  double **Q = ddmalloc(ntot,ne);
  double **F = ddmalloc(ntot,ne);
  // ptwise eigenvalues and speeds

  //double lam[ne], amaxl[ntot][ne], amaxg[ne], a;
  double **amaxl = ddmalloc(ntot,ne);
  double *lam   = dmalloc(ne); 
  double *amaxg = dmalloc(ne);  
  double a;


  // eigenvectors at interfaces
  double *L=dmalloc(ne*ne); double* R = dmalloc(ne*ne); // left/right eigenv, ptwise 5 x 5
  
  // kth-characteristic fields at interfaces
  double *buffer= dmalloc(2*ntot);
  double *fpk, *fmk; 

  // tmp
  double *fsum = dmalloc(ne); double *fterm = dmalloc(ne);
  double tmp1,tmp2,tmp3;
  const double rec_switch = MATTER.flx_LLF_HO_RHO * GRHD.ATM_FATM * GRHD.ATM_RHOATM;

  // init
  fpk = &(buffer[       of]);
  fmk = &(buffer[ntot + of]);
  for( k = 0 ; k < ne ; k++ ) amaxg[k] = ZERO;

  // rem offset all the arrays of dimension = ntot (except fpk fmk)



  // ----------------------------------------------------------
  // step 1: compute phys fluxes and speeds ptwise 
  // ----------------------------------------------------------

  for( i = -of ; i <= nx+ngx; i++ ) {
    
    // inverse metric and detg
    grhd_compute_detg_invg_pt(gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],
			      &detg_i,
			      &gupxx_i,&tmp1,&tmp2,&gupyy_i,&tmp3,&gupzz_i);
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx[i],vy[i],vz[i], 
		       gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],
		       &vlowx_i,&vlowy_i,&vlowz_i, &v2_i );

    // cure here not necessary since eigenvalues are cured anyway
    // and conservatives and fluxes are taken from evolved vars
    v2_i = DMIN(v2_i,GRHD.HRSC_VMAX); 
    if (v2_i==GRHD.HRSC_VMAX) Wlor_i = GRHD.HRSC_WlorMAX;
    else Wlor_i  = ONE/sqrt( ONE - v2_i );

    // EoS
    //GRHD.use_eos(&pres_i,&cs2_i, &kappa_i,&chi_i, rho[i],epsl[i]);
    EOS.comp("re","","","pc","re","", rho[i],epsl[i], &pres_i,&cs2_i,&chi_i,&kappa_i);
    
    
    // eigenvalues
    grhd_eig(alpha[i], beta[i], vx[i],v2_i,gupxx_i, cs2_i,
	     &(lam[0]),&(lam[3]),&(lam[4]));
    lam[1] = lam[2] = lam[0];
    
    // speeds 
    for( k = 0 ; k < ne ; k++ ) {
      amaxl[i+of][k] = fabs(lam[k]);
      amaxg[k] = DMAX(amaxg[k], amaxl[i+of][k]);
    }
        
    // phys fluxes 
    // ordering is changed !
    Q[i+of][0] = conD [i];
    Q[i+of][4] = conT [i];     
    Q[i+of][1] = conSx[i]; 
    Q[i+of][2] = conSy[i]; 
    Q[i+of][3] = conSz[i]; 
    
    grhd_compute_fx_pt(alpha[i],beta[i], detg_i, pres_i, vx[i], 
		       Q[i+of][0], Q[i+of][4], Q[i+of][1], Q[i+of][2], Q[i+of][3],
		       &(F[i+of][0]), &(F[i+of][4]), &(F[i+of][1]), &(F[i+of][2]), &(F[i+of][3]));
    
    // cure ptwise fluxes
    // ...
    
#if (!(HO_TREAT_AVGS==0))
    
    // Store arrays for averages
    // needed ONLY IF opt 1. or 2. below active !
    detga[i+of] = detg_i;
    gupxx[i+of] = gupxx_i;
    vlowx[i+of] = vlowx_i;
    vlowy[i+of] = vlowy_i;
    vlowz[i+of] = vlowz_i;
    v2   [i+of] = v2_i;
    Wlor [i+of] = Wlor_i;
    h    [i+of] = ONE + epsl[i] + pres_i/rho[i];
    kappa[i+of] = kappa_i;
    cs2  [i+of] = cs2_i;

#endif


  } // end loop i
 
 
  // ----------------------------------------------------------
  // step 2: compute numerical fluxes at interfaces 
  // ----------------------------------------------------------

  for( i = 0 ; i <= nx ; i++ ) {
    
    // if pt and NN = atm, set flx = 0, loop
    if ((MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) { 
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i] = ZERO;
      continue;
    }
        
    // compute eigenvectors/eigenvalues at interfaces i+1/2 
    // from Roe avgs
    
    // averages
    rho_avg   = HALF * ( rho  [i] + rho  [i+1] );    
    epsl_avg  = HALF * ( epsl [i] + epsl [i+1] );
    vx_avg    = HALF * ( vx   [i] + vx   [i+1] );
    vy_avg    = HALF * ( vy   [i] + vy   [i+1] );
    vz_avg    = HALF * ( vz   [i] + vz   [i+1] );
    
    alpha_avg = HALF * ( alpha[i] + alpha[i+1] );
    betax_avg = HALF * ( beta [i] + beta [i+1] );
    gxx_avg   = HALF * ( gxx  [i] + gxx  [i+1] );
    gxy_avg   = HALF * ( gxy  [i] + gxy  [i+1] );
    gxz_avg   = HALF * ( gxz  [i] + gxz  [i+1] );
    gyy_avg   = HALF * ( gyy  [i] + gyy  [i+1] );
    gyz_avg   = HALF * ( gyz  [i] + gyz  [i+1] );
    gzz_avg   = HALF * ( gzz  [i] + gzz  [i+1] );
    
    // additional averages are needed   
    // PROBLEM: averages of nonlinear function are different from 
    // nonlinear functions of averages ... what to use ? 
    // Note everything is mixed: metric, primitives, eigenvalues, eigenvectors, 
    // inconsistencies arise !
    
#if (HO_TREAT_AVGS==0)
    
    // opt 0. compute everything else from above avgs 

    // inverse metric and detg
    grhd_compute_detg_invg_pt(gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
			      &detg_avg, 
			      &gupxx_avg,&tmp1,&tmp2,&gupyy_avg,&tmp3,&gupzz_avg);
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx_avg,vy_avg,vz_avg, 
		       gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
		       &vlowx_avg,&vlowy_avg,&vlowz_avg, &v2_avg );

    // EoS
    //GRHD.use_eos(&pres_avg, &cs2_avg,&kappa_avg,&chi_avg, rho_avg,epsl_avg);
    EOS.comp("re","","","pc","re","", rho_avg,epsl_avg, &pres_avg,&cs2_avg,&chi_avg,&kappa_avg);
    
    
    h_avg     = ONE + epsl_avg + pres_avg/rho_avg;
    
    //if ((MATTER_USEMASK) && ((mask[i]>0.4) || (mask[i+1]>0.4))) {
    //  // use bounded vals
    //  pres_avg  = HALF * ( pres [i]    + pres [i+1] );
    //  h_avg     = HALF * ( h    [i+of] + h    [i+1+of] );      
    //  kappa_avg = HALF * ( kappa[i+of] + kappa[i+1+of] );
    //  cs2_avg   = HALF * ( cs2  [i+of] + cs2  [i+1+of] );      
    //}

    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );        
    
#elif (HO_TREAT_AVGS==1)
    
    // opt 1. average only gxxup and detg, others consistently from avgs:
    gupxx_avg = HALF * ( gupxx[i+of] + gupxx[i+1+of] ); 
    detg_avg  = HALF * ( detga[i+of] + detga[i+1+of] ); 
    
    GRHD.use_eos(&pres_avg, &cs2_avg,&kappa_avg,&chi_avg, rho_avg,epsl_avg);
    h_avg     = ONE + epsl_avg + pres_avg/rho_avg;
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx_avg,vy_avg,vz_avg, gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
		       &vlowx_avg,&vlowy_avg,&vlowz_avg, &v2_avg );
    
    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );

#elif (HO_TREAT_AVGS==2)
    
    // opt 2. average all:
    gupxx_avg = HALF * ( gupxx[i+of] + gupxx[i+1+of] ); 
    detg_avg  = HALF * ( detga[i+of] + detga[i+1+of] ); 
    
    pres_avg  = HALF * ( pres [i]    + pres [i+1]    ); 
    h_avg     = HALF * ( h    [i+of] + h    [i+1+of] );      
    kappa_avg = HALF * ( kappa[i+of] + kappa[i+1+of] );
    cs2_avg   = HALF * ( cs2  [i+of] + cs2  [i+1+of] );
    
    v2_avg    = HALF * ( v2   [i+of] + v2   [i+1+of] );
    vlowx_avg = HALF * ( vlowx[i+of] + vlowx[i+1+of] );
    vlowy_avg = HALF * ( vlowy[i+of] + vlowy[i+1+of] );
    vlowz_avg = HALF * ( vlowz[i+of] + vlowz[i+1+of] );
    Wlor_avg  = HALF * ( Wlor [i+of] + Wlor [i+1+of] );

    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );
    
#else
    errorexit("you must choose a way to treat avg in HO flux!");
#endif
          
    // eigenvalues
    grhd_eig( alpha_avg, betax_avg, vx_avg, v2_avg, gupxx_avg, cs2_avg,
	      &lam[0], &lam[3], &lam[4] );
    lam[1] = lam[2] = lam[0];
    
    
    // eigenvectors     
    grhd_eigLR( gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg, 
		alpha_avg, betax_avg,  gupxx_avg, detg_avg, 
		h_avg, Wlor_avg, rho_avg, cs2_avg, kappa_avg, 
		vx_avg, vy_avg, vz_avg, vlowx_avg, vlowy_avg,  vlowz_avg,
		lam[0],  lam[3],  lam[4], 
		L, R );
        
    // compute fluxes at i+1/2      
    for( k = 0 ; k < ne ; k++ ) {          

      // flux contribute from kth char var
      
      // speed
#if (HO_TREAT_SPEED==0)
      // local flux splitting      
      a = ZERO;
      for (l = i-S; l <= i+S; l++) a = DMAX(a,amaxl[l+of][k]); 
#elif (HO_TREAT_SPEED==1)
      // global f-splitting
      a = amaxg[k];        
#else
      a = ONE; 
#endif
      
      //a = disschi * a; // add dissipation

      
      // compute pts needed in the rec (stencil S)
      for (l = i-S; l <= i+S; l++) {	
	fpk[l] = fmk[l] = ZERO;	
	for( j = 0 ; j < ne ; j++ ) {
	  fpk[l] += HALF * L[k*ne+j] * (F[l      +of][j] + a*Q[l      +of][j]);
	  fmk[l] += HALF * L[k*ne+j] * (F[2*i-l+1+of][j] - a*Q[2*i-l+1+of][j]);
	}	
      }    
      
      // reconstruct interfaces values of kth char var
      //fterm[k] = MATTER.rec1dl(fpk,i) + MATTER.rec1dl(fmk,i); 

if (rho_avg < 1000*rec_switch) // better condition would look to rho_avg and its neighbours ...
// LO rec
//fterm[k] = rec1d_p_lintvd(fpk,i) + rec1d_p_lintvd(fmk,i);
 fterm[k] = rec1d_p_wenoz(fpk,i) + rec1d_p_wenoz(fmk,i);
else
// HO rec
fterm[k] = MATTER.rec1dl(fpk,i) + MATTER.rec1dl(fmk,i); // set to a high order MP5

      //if (!finite(fterm[k])) fterm[k] = 0.;

    }
        
    for( j = 0 ; j < ne ; j++ ) {

      // add kth contributes to fluxes 
      fsum[j] = ZERO;   
      for( k = 0 ; k < ne ; k++ ) fsum[j] += fterm[k]*R[k*ne+j];

    }


    // store fluxes at interface i+1/2, 
    // (back to original ordering)
    flux1d[0][i] =  fsum[0];
    flux1d[1][i] =  fsum[4];
    flux1d[2][i] =  fsum[1];
    flux1d[3][i] =  fsum[2];
    flux1d[4][i] =  fsum[3];


    // cure fluxes ?
    if ((!finite(fsum[0])) || 
	(!finite(fsum[1])) || 
	(!finite(fsum[2])) ||
	(!finite(fsum[3])) ||
	(!finite(fsum[4])))
      fsum[0] = fsum[1] = fsum[2] = fsum[3] = fsum[4] = ZERO;


    if (PR) if (CheckForNANandINF(5, flux1d[0][i],flux1d[1][i],flux1d[2][i],flux1d[3][i],flux1d[4][i])) {	

	// redo calculations and display:
	printf("problem inside ho flux\n");

	printf("  var \t F_i \n");
	for (k=0; k<ne; k++) printf(" %d \t %e\n",k,flux1d[k][i]);
	
	printf("----- L\n");
	for( k = 0 ; k < ne ; k++ ) {
	  printf(" k=%d\n",k);
	  for( j = 0 ; j < ne ; j++ ) printf(" %e", L[k*ne + j]);
	  printf("\n"); 
	}
	printf("-----\n");
	printf("----- R\n");
	for( k = 0 ; k < ne ; k++ ) {
	  printf(" k=%d\n",k);
	  for( j = 0 ; j < ne ; j++ ) printf(" %e", R[k*ne + j]);
	  printf("\n"); 
	}
	printf("-----\n");

	printf(" i=%d sound speed avg = %e\n",0,cs2_avg);
	

	for( k = 0 ; k < ne ; k++ ) {
	  
	  a = ZERO;
	  for (l = i-S; l <= i+S; l++) a = DMAX(a,amaxl[l+of][k]); 	  
	  printf(" |==============> char var no k=%d \t speed a=%e \n",k,a);
	  
	  for (l = i-S; l <= i+S; l++) {
	    for( j = 0 ; j < ne ; j++ ) {
	      printf(" i=%d \t j=%d F = %e  F = %e\n",i,j,F[l      +of][j],F[2*i-l+1+of][j]);
              printf(" i=%d \t j=%d F = %e  F = %e\n",i,j,Q[l      +of][j],Q[2*i-l+1+of][j]);
	      fpk[l] += HALF * L[k*ne+j] * (F[l      +of][j] + a*Q[l      +of][j]);
	      fmk[l] += HALF * L[k*ne+j] * (F[2*i-l+1+of][j] - a*Q[2*i-l+1+of][j]);	   
	    }
            printf(" i=%d \t j=%d fp = %e  fm = %e\n",i,j,fpk[l],fmk[l]);
	  }
	  
	  printf(" |==============> char var no k=%d rec: \n",k);
          printf(" i=%d \t j=%d rec fp = %e  rec fm = %e\n",i,j,MATTER.rec1dl(fpk,i),MATTER.rec1dl(fmk,i));
	  
	}
	
	errorexits(" (PR=%d) I stop the code",PR);
	
      }
    
  } // end loop i+1/2
  

#if (!(HO_TREAT_AVGS==0))

  free(detga); free(gupxx); 
  free(vlowx); free(vlowy); free(vlowz);
  free(v2); free(Wlow);
  free(h); free(kappa); free(cs2);

#endif

  free(lam); free(amaxg); 
  free(L); free(R); 
  free(buffer); 
  free(fsum); 
  free(fterm); 
  ddfree(amaxl,ntot); 
  ddfree(Q,ntot); 
  ddfree(F,ntot); 

} 




// high-order llf-weno fluxes computed from characteristic vars
void grhd_flx1d_ho(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d)
{
  
  int i,l; // -ngx+1 ... nx+ngx / 0 ... nx
  int j,k; //  0 ... ne
  
  const int ne     = MATTER.NVf; // no eigenv 
  const int ntot   = nx + 2*ngx;
  const int of     = ngx-1;      // offset
  
  const double disschi = 1.; // dissipation factor [1,1.3]
  const int S = MATTER.HRSC_NGHOST/2; // = 2 // rec stencil 


  // metric
  double *alpha = g1d[0];
  double *beta  = g1d[1]; // beta^(x)
  double *gxx   = g1d[2];
  double *gxy   = g1d[3];
  double *gxz   = g1d[4];
  double *gyy   = g1d[5];
  double *gyz   = g1d[6];
  double *gzz   = g1d[7];

  // conservatives
  double *conD  = q1d[0];
  double *conT  = q1d[1]; 
  double *conSx = q1d[2];
  double *conSy = q1d[3];
  double *conSz = q1d[4];

  // primitives
  double *rho  = w1d[0];
  double *epsl = w1d[1]; 
  double *vx   = w1d[2];
  double *vy   = w1d[3];
  double *vz   = w1d[4];
  double *pres = w1d[5];

  // ptwise metric and primitives at interfaces i+1/2
  double rho_avg, epsl_avg, vx_avg, vy_avg, vz_avg, pres_avg;
  double cs2_avg, h_avg, Wlor_avg; 
  double v2_avg, vlowx_avg, vlowy_avg, vlowz_avg;  
  double kappa_avg, chi_avg;
  double gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg;
  double alpha_avg, betax_avg, detg_avg, gupxx_avg,gupyy_avg,gupzz_avg;

  // ptwise metric and primitives at pt i
  double v2_i, vlowx_i, vlowy_i, vlowz_i, Wlor_i;    
  double alpha_i, beta_i, detg_i, sqrtg_i, detginv, gupxx_i,gupyy_i,gupzz_i;
  double pres_i, kappa_i, chi_i, cs2_i, avb_i;

#if (!(HO_TREAT_AVGS==0))

  // vars to be used for avgs
  double *detga = dmalloc(ntot); double *gupxx = dmalloc(ntot);
  double *vlowx = dmalloc(ntot); double *vlowy = dmalloc(ntot), 
  double *vlowz = dmalloc(ntot); double *v2 = dmalloc(ntot); double *Wlor = dmalloc(ntot);
  double *h = dmalloc(ntot); double *kappa = dmalloc(ntot); double *cs2 = dmalloc(ntot);
  
#endif

  // conservative and fluxes at pts i
  double **Q = ddmalloc(ntot,ne);
  double **F = ddmalloc(ntot,ne);
  // ptwise eigenvalues and speeds

  //double lam[ne], amaxl[ntot][ne], amaxg[ne], a;
  double **amaxl = ddmalloc(ntot,ne);
  double *lam   = dmalloc(ne); 
  double *amaxg = dmalloc(ne);  
  double a;


  // eigenvectors at interfaces
  double *L=dmalloc(ne*ne); double* R = dmalloc(ne*ne); // left/right eigenv, ptwise 5 x 5
  
  // kth-characteristic fields at interfaces
  double *buffer= dmalloc(2*ntot);
  double *fpk, *fmk; 

  // tmp
  double *fsum = dmalloc(ne); double *fterm = dmalloc(ne);
  double tmp1,tmp2,tmp3;

  // init
  fpk = &(buffer[       of]);
  fmk = &(buffer[ntot + of]);
  for( k = 0 ; k < ne ; k++ ) amaxg[k] = ZERO;

  // rem offset all the arrays of dimension = ntot (except fpk fmk)



  // ----------------------------------------------------------
  // step 1: compute phys fluxes and speeds ptwise 
  // ----------------------------------------------------------

  for( i = -of ; i <= nx+ngx; i++ ) {
    
    // inverse metric and detg
    grhd_compute_detg_invg_pt(gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],
			      &detg_i,
			      &gupxx_i,&tmp1,&tmp2,&gupyy_i,&tmp3,&gupzz_i);
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx[i],vy[i],vz[i], 
		       gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],
		       &vlowx_i,&vlowy_i,&vlowz_i, &v2_i );

    // cure here not necessary since eigenvalues are cured anyway
    // and conservatives and fluxes are taken from evolved vars
    v2_i = DMIN(v2_i,GRHD.HRSC_VMAX); 
    if (v2_i==GRHD.HRSC_VMAX) Wlor_i = GRHD.HRSC_WlorMAX;
    else Wlor_i  = ONE/sqrt( ONE - v2_i );

    // EoS
    //GRHD.use_eos(&pres_i,&cs2_i, &kappa_i,&chi_i, rho[i],epsl[i]);
    EOS.comp("re","","","pc","re","", rho[i],epsl[i], &pres_i,&cs2_i,&chi_i,&kappa_i);
    
    
    // eigenvalues
    grhd_eig(alpha[i], beta[i], vx[i],v2_i,gupxx_i, cs2_i,
	     &(lam[0]),&(lam[3]),&(lam[4]));
    lam[1] = lam[2] = lam[0];
    
    // speeds 
    for( k = 0 ; k < ne ; k++ ) {
      amaxl[i+of][k] = fabs(lam[k]);
      amaxg[k] = DMAX(amaxg[k], amaxl[i+of][k]);
    }
        
    // phys fluxes 
    // ordering is changed !
    Q[i+of][0] = conD [i];
    Q[i+of][4] = conT [i];     
    Q[i+of][1] = conSx[i]; 
    Q[i+of][2] = conSy[i]; 
    Q[i+of][3] = conSz[i]; 
    
    grhd_compute_fx_pt(alpha[i],beta[i], detg_i, pres_i, vx[i], 
		       Q[i+of][0], Q[i+of][4], Q[i+of][1], Q[i+of][2], Q[i+of][3],
		       &(F[i+of][0]), &(F[i+of][4]), &(F[i+of][1]), &(F[i+of][2]), &(F[i+of][3]));
    
    // cure ptwise fluxes
    // ...
    
#if (!(HO_TREAT_AVGS==0))
    
    // Store arrays for averages
    // needed ONLY IF opt 1. or 2. below active !
    detga[i+of] = detg_i;
    gupxx[i+of] = gupxx_i;
    vlowx[i+of] = vlowx_i;
    vlowy[i+of] = vlowy_i;
    vlowz[i+of] = vlowz_i;
    v2   [i+of] = v2_i;
    Wlor [i+of] = Wlor_i;
    h    [i+of] = ONE + epsl[i] + pres_i/rho[i];
    kappa[i+of] = kappa_i;
    cs2  [i+of] = cs2_i;

#endif


  } // end loop i
 
 
  // ----------------------------------------------------------
  // step 2: compute numerical fluxes at interfaces 
  // ----------------------------------------------------------

  for( i = 0 ; i <= nx ; i++ ) {
    
    // if pt and NN = atm, set flx = 0, loop
    if ((MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) { 
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i] = ZERO;
      continue;
    }
        
    // compute eigenvectors/eigenvalues at interfaces i+1/2 
    // from Roe avgs
    
    // averages
    rho_avg   = HALF * ( rho  [i] + rho  [i+1] );    
    epsl_avg  = HALF * ( epsl [i] + epsl [i+1] );
    vx_avg    = HALF * ( vx   [i] + vx   [i+1] );
    vy_avg    = HALF * ( vy   [i] + vy   [i+1] );
    vz_avg    = HALF * ( vz   [i] + vz   [i+1] );
    
    alpha_avg = HALF * ( alpha[i] + alpha[i+1] );
    betax_avg = HALF * ( beta [i] + beta [i+1] );
    gxx_avg   = HALF * ( gxx  [i] + gxx  [i+1] );
    gxy_avg   = HALF * ( gxy  [i] + gxy  [i+1] );
    gxz_avg   = HALF * ( gxz  [i] + gxz  [i+1] );
    gyy_avg   = HALF * ( gyy  [i] + gyy  [i+1] );
    gyz_avg   = HALF * ( gyz  [i] + gyz  [i+1] );
    gzz_avg   = HALF * ( gzz  [i] + gzz  [i+1] );
    
    // additional averages are needed   
    // PROBLEM: averages of nonlinear function are different from 
    // nonlinear functions of averages ... what to use ? 
    // Note everything is mixed: metric, primitives, eigenvalues, eigenvectors, 
    // inconsistencies arise !
    
#if (HO_TREAT_AVGS==0)
    
    // opt 0. compute everything else from above avgs 

    // inverse metric and detg
    grhd_compute_detg_invg_pt(gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
			      &detg_avg, 
			      &gupxx_avg,&tmp1,&tmp2,&gupyy_avg,&tmp3,&gupzz_avg);
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx_avg,vy_avg,vz_avg, 
		       gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
		       &vlowx_avg,&vlowy_avg,&vlowz_avg, &v2_avg );

    // EoS
    //GRHD.use_eos(&pres_avg, &cs2_avg,&kappa_avg,&chi_avg, rho_avg,epsl_avg);
    EOS.comp("re","","","pc","re","", rho_avg,epsl_avg, &pres_avg,&cs2_avg,&chi_avg,&kappa_avg);
    
    
    h_avg     = ONE + epsl_avg + pres_avg/rho_avg;
    
    //if ((MATTER_USEMASK) && ((mask[i]>0.4) || (mask[i+1]>0.4))) {
    //  // use bounded vals
    //  pres_avg  = HALF * ( pres [i]    + pres [i+1] );
    //  h_avg     = HALF * ( h    [i+of] + h    [i+1+of] );      
    //  kappa_avg = HALF * ( kappa[i+of] + kappa[i+1+of] );
    //  cs2_avg   = HALF * ( cs2  [i+of] + cs2  [i+1+of] );      
    //}

    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );        
    
#elif (HO_TREAT_AVGS==1)
    
    // opt 1. average only gxxup and detg, others consistently from avgs:
    gupxx_avg = HALF * ( gupxx[i+of] + gupxx[i+1+of] ); 
    detg_avg  = HALF * ( detga[i+of] + detga[i+1+of] ); 
    
    GRHD.use_eos(&pres_avg, &cs2_avg,&kappa_avg,&chi_avg, rho_avg,epsl_avg);
    h_avg     = ONE + epsl_avg + pres_avg/rho_avg;
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx_avg,vy_avg,vz_avg, gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
		       &vlowx_avg,&vlowy_avg,&vlowz_avg, &v2_avg );
    
    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );

#elif (HO_TREAT_AVGS==2)
    
    // opt 2. average all:
    gupxx_avg = HALF * ( gupxx[i+of] + gupxx[i+1+of] ); 
    detg_avg  = HALF * ( detga[i+of] + detga[i+1+of] ); 
    
    pres_avg  = HALF * ( pres [i]    + pres [i+1]    ); 
    h_avg     = HALF * ( h    [i+of] + h    [i+1+of] );      
    kappa_avg = HALF * ( kappa[i+of] + kappa[i+1+of] );
    cs2_avg   = HALF * ( cs2  [i+of] + cs2  [i+1+of] );
    
    v2_avg    = HALF * ( v2   [i+of] + v2   [i+1+of] );
    vlowx_avg = HALF * ( vlowx[i+of] + vlowx[i+1+of] );
    vlowy_avg = HALF * ( vlowy[i+of] + vlowy[i+1+of] );
    vlowz_avg = HALF * ( vlowz[i+of] + vlowz[i+1+of] );
    Wlor_avg  = HALF * ( Wlor [i+of] + Wlor [i+1+of] );

    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );
    
#else
    errorexit("you must choose a way to treat avg in HO flux!");
#endif
          
    // eigenvalues
    grhd_eig( alpha_avg, betax_avg, vx_avg, v2_avg, gupxx_avg, cs2_avg,
	      &lam[0], &lam[3], &lam[4] );
    lam[1] = lam[2] = lam[0];
    
    
    // eigenvectors     
    grhd_eigLR( gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg, 
		alpha_avg, betax_avg,  gupxx_avg, detg_avg, 
		h_avg, Wlor_avg, rho_avg, cs2_avg, kappa_avg, 
		vx_avg, vy_avg, vz_avg, vlowx_avg, vlowy_avg,  vlowz_avg,
		lam[0],  lam[3],  lam[4], 
		L, R );
        
    // compute fluxes at i+1/2      
    for( k = 0 ; k < ne ; k++ ) {          

      // flux contribute from kth char var
      
      // speed
#if (HO_TREAT_SPEED==0)
      // local flux splitting      
      a = ZERO;
      for (l = i-S; l <= i+S; l++) a = DMAX(a,amaxl[l+of][k]); 
#elif (HO_TREAT_SPEED==1)
      // global f-splitting
      a = amaxg[k];        
#else
      a = ONE; 
#endif

      //a = disschi * a; // add dissipation
      
      // compute pts needed in the rec (stencil S)
      for (l = i-S; l <= i+S; l++) {	
	fpk[l] = fmk[l] = ZERO;	
	for( j = 0 ; j < ne ; j++ ) {
	  fpk[l] += HALF * L[k*ne+j] * (F[l      +of][j] + a*Q[l      +of][j]);
	  fmk[l] += HALF * L[k*ne+j] * (F[2*i-l+1+of][j] - a*Q[2*i-l+1+of][j]);
	}	
      }    
      
      // reconstruct interfaces values of kth char var
      fterm[k] = MATTER.rec1dl(fpk,i) + MATTER.rec1dl(fmk,i); 

      //if (!finite(fterm[k])) fterm[k] = 0.;

    }
        
    for( j = 0 ; j < ne ; j++ ) {

      // add kth contributes to fluxes 
      fsum[j] = ZERO;   
      for( k = 0 ; k < ne ; k++ ) fsum[j] += fterm[k]*R[k*ne+j];

    }


    // store fluxes at interface i+1/2, 
    // (back to original ordering)
    flux1d[0][i] =  fsum[0];
    flux1d[1][i] =  fsum[4];
    flux1d[2][i] =  fsum[1];
    flux1d[3][i] =  fsum[2];
    flux1d[4][i] =  fsum[3];


    // cure fluxes ?
    if ((!finite(fsum[0])) || 
	(!finite(fsum[1])) || 
	(!finite(fsum[2])) ||
	(!finite(fsum[3])) ||
	(!finite(fsum[4])))
      fsum[0] = fsum[1] = fsum[2] = fsum[3] = fsum[4] = ZERO;


    if (PR) if (CheckForNANandINF(5, flux1d[0][i],flux1d[1][i],flux1d[2][i],flux1d[3][i],flux1d[4][i])) {	

	// redo calculations and display:
	printf("problem inside ho flux\n");

	printf("  var \t F_i \n");
	for (k=0; k<ne; k++) printf(" %d \t %e\n",k,flux1d[k][i]);
	
	printf("----- L\n");
	for( k = 0 ; k < ne ; k++ ) {
	  printf(" k=%d\n",k);
	  for( j = 0 ; j < ne ; j++ ) printf(" %e", L[k*ne + j]);
	  printf("\n"); 
	}
	printf("-----\n");
	printf("----- R\n");
	for( k = 0 ; k < ne ; k++ ) {
	  printf(" k=%d\n",k);
	  for( j = 0 ; j < ne ; j++ ) printf(" %e", R[k*ne + j]);
	  printf("\n"); 
	}
	printf("-----\n");

	printf(" i=%d sound speed avg = %e\n",0,cs2_avg);
	

	for( k = 0 ; k < ne ; k++ ) {
	  
	  a = ZERO;
	  for (l = i-S; l <= i+S; l++) a = DMAX(a,amaxl[l+of][k]); 	  
	  printf(" |==============> char var no k=%d \t speed a=%e \n",k,a);
	  
	  for (l = i-S; l <= i+S; l++) {
	    for( j = 0 ; j < ne ; j++ ) {
	      printf(" i=%d \t j=%d F = %e  F = %e\n",i,j,F[l      +of][j],F[2*i-l+1+of][j]);
              printf(" i=%d \t j=%d F = %e  F = %e\n",i,j,Q[l      +of][j],Q[2*i-l+1+of][j]);
	      fpk[l] += HALF * L[k*ne+j] * (F[l      +of][j] + a*Q[l      +of][j]);
	      fmk[l] += HALF * L[k*ne+j] * (F[2*i-l+1+of][j] - a*Q[2*i-l+1+of][j]);	   
	    }
            printf(" i=%d \t j=%d fp = %e  fm = %e\n",i,j,fpk[l],fmk[l]);
	  }
	  
	  printf(" |==============> char var no k=%d rec: \n",k);
          printf(" i=%d \t j=%d rec fp = %e  rec fm = %e\n",i,j,MATTER.rec1dl(fpk,i),MATTER.rec1dl(fmk,i));
	  
	}
	
	errorexits(" (PR=%d) I stop the code",PR);
	
      }
    
  } // end loop i+1/2
  

#if (!(HO_TREAT_AVGS==0))

  free(detga); free(gupxx); 
  free(vlowx); free(vlowy); free(vlowz);
  free(v2); free(Wlow);
  free(h); free(kappa); free(cs2);

#endif

  free(lam); free(amaxg); 
  free(L); free(R); 
  free(buffer); 
  free(fsum); 
  free(fterm); 
  ddfree(amaxl,ntot); 
  ddfree(Q,ntot); 
  ddfree(F,ntot); 

}



// high-order llf-weno fluxes computed from characteristic vars
void grhd_flx1d_ho_highdens(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d)
{
  
  int i,l; // -ngx+1 ... nx+ngx / 0 ... nx
  int j,k; //  0 ... ne
  
  const int ne     = MATTER.NVf; // no eigenv 
  const int ntot   = nx + 2*ngx;
  const int of     = ngx-1;      // offset
  
  const double disschi = 1.; // dissipation factor [1,1.3]
  const int S = MATTER.HRSC_NGHOST/2; // = 2 // rec stencil 


  // metric
  double *alpha = g1d[0];
  double *beta  = g1d[1]; // beta^(x)
  double *gxx   = g1d[2];
  double *gxy   = g1d[3];
  double *gxz   = g1d[4];
  double *gyy   = g1d[5];
  double *gyz   = g1d[6];
  double *gzz   = g1d[7];

  // conservatives
  double *conD  = q1d[0];
  double *conT  = q1d[1]; 
  double *conSx = q1d[2];
  double *conSy = q1d[3];
  double *conSz = q1d[4];

  // primitives
  double *rho  = w1d[0];
  double *epsl = w1d[1]; 
  double *vx   = w1d[2];
  double *vy   = w1d[3];
  double *vz   = w1d[4];
  double *pres = w1d[5];

  // ptwise metric and primitives at interfaces i+1/2
  double rho_avg, epsl_avg, vx_avg, vy_avg, vz_avg, pres_avg;
  double cs2_avg, h_avg, Wlor_avg; 
  double v2_avg, vlowx_avg, vlowy_avg, vlowz_avg;  
  double kappa_avg, chi_avg;
  double gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg;
  double alpha_avg, betax_avg, detg_avg, gupxx_avg,gupyy_avg,gupzz_avg;

  // ptwise metric and primitives at pt i
  double v2_i, vlowx_i, vlowy_i, vlowz_i, Wlor_i;    
  double alpha_i, beta_i, detg_i, sqrtg_i, detginv, gupxx_i,gupyy_i,gupzz_i;
  double pres_i, kappa_i, chi_i, cs2_i, avb_i;

#if (!(HO_TREAT_AVGS==0))

  // vars to be used for avgs
  double *detga = dmalloc(ntot); double *gupxx = dmalloc(ntot);
  double *vlowx = dmalloc(ntot); double *vlowy = dmalloc(ntot), 
  double *vlowz = dmalloc(ntot); double *v2 = dmalloc(ntot); double *Wlor = dmalloc(ntot);
  double *h = dmalloc(ntot); double *kappa = dmalloc(ntot); double *cs2 = dmalloc(ntot);
  
#endif

  // conservative and fluxes at pts i
  double **Q = ddmalloc(ntot,ne);
  double **F = ddmalloc(ntot,ne);
  // ptwise eigenvalues and speeds

  //double lam[ne], amaxl[ntot][ne], amaxg[ne], a;
  double **amaxl = ddmalloc(ntot,ne);
  double *lam   = dmalloc(ne); 
  double *amaxg = dmalloc(ne);  
  double a;


  // eigenvectors at interfaces
  double *L=dmalloc(ne*ne); double* R = dmalloc(ne*ne); // left/right eigenv, ptwise 5 x 5
  
  // kth-characteristic fields at interfaces
  double *buffer= dmalloc(2*ntot);
  double *fpk, *fmk; 

  // tmp
  double *fsum = dmalloc(ne); double *fterm = dmalloc(ne);

  double rhomins, epslmins, vxmins, vymins, vzmins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus;
  
  double alphaplus, betaplus, detgplus, gupxxplus, gupyyplus, gupzzplus;
  double alphamins, betamins, detgmins, gupxxmins, gupyymins, gupzzmins;
  double gxxplus, gxyplus, gxzplus, gyyplus, gyzplus, gzzplus;
  double gxxmins, gxymins, gxzmins, gyymins, gyzmins, gzzmins;
  
  double pres_L, v2_L, Wlor_L, W2rhoh_L, cs2_L, vlowx_L, vlowy_L, vlowz_L;
  double pres_R, v2_R, Wlor_R, W2rhoh_R, cs2_R, vlowx_R, vlowy_R, vlowz_R;

  double tmp1,tmp2,tmp3;

  // physical conservative and fluxes at interfaces
/*  double qL[MATTER.NVq], qR[MATTER.NVq];
  double fL[MATTER.NVf], fR[MATTER.NVf], f_i[MATTER.NVf]; 
  double lamL[MATTER.NVf], lamR[MATTER.NVf]; 
*/
  double *qL   = dmalloc(MATTER.NVq);
  double *qR   = dmalloc(MATTER.NVq);  
  double *fL   = dmalloc(MATTER.NVf);  
  double *fR   = dmalloc(MATTER.NVf);  
  double *f_i  = dmalloc(MATTER.NVf);  
  double *lamL = dmalloc(MATTER.NVf);
  double *lamR = dmalloc(MATTER.NVf); 
  const double rec_switch = MATTER.flx_LLF_HO_RHO * GRHD.ATM_FATM * GRHD.ATM_RHOATM;
  const double alpha_switch = MATTER.flx_LLF_HO_ALPHA;

  // init
  fpk = &(buffer[       of]);
  fmk = &(buffer[ntot + of]);
  for( k = 0 ; k < ne ; k++ ) amaxg[k] = ZERO;

  // rem offset all the arrays of dimension = ntot (except fpk fmk)



  // ----------------------------------------------------------
  // step 1: compute phys fluxes and speeds ptwise 
  // ----------------------------------------------------------

  for( i = -of ; i <= nx+ngx; i++ ) {
    
    // inverse metric and detg
    grhd_compute_detg_invg_pt(gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],
			      &detg_i,
			      &gupxx_i,&tmp1,&tmp2,&gupyy_i,&tmp3,&gupzz_i);
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx[i],vy[i],vz[i], 
		       gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],
		       &vlowx_i,&vlowy_i,&vlowz_i, &v2_i );

    // cure here not necessary since eigenvalues are cured anyway
    // and conservatives and fluxes are taken from evolved vars
    v2_i = DMIN(v2_i,GRHD.HRSC_VMAX); 
    if (v2_i==GRHD.HRSC_VMAX) Wlor_i = GRHD.HRSC_WlorMAX;
    else Wlor_i  = ONE/sqrt( ONE - v2_i );

    // EoS
    //GRHD.use_eos(&pres_i,&cs2_i, &kappa_i,&chi_i, rho[i],epsl[i]);
    EOS.comp("re","","","pc","re","", rho[i],epsl[i], &pres_i,&cs2_i,&chi_i,&kappa_i);
    
    
    // eigenvalues
    grhd_eig(alpha[i], beta[i], vx[i],v2_i,gupxx_i, cs2_i,
	     &(lam[0]),&(lam[3]),&(lam[4]));
    lam[1] = lam[2] = lam[0];
    
    // speeds 
    for( k = 0 ; k < ne ; k++ ) {
      amaxl[i+of][k] = fabs(lam[k]);
      amaxg[k] = DMAX(amaxg[k], amaxl[i+of][k]);
    }
        
    // phys fluxes 
    // ordering is changed !
    Q[i+of][0] = conD [i];
    Q[i+of][4] = conT [i];     
    Q[i+of][1] = conSx[i]; 
    Q[i+of][2] = conSy[i]; 
    Q[i+of][3] = conSz[i]; 
    
    grhd_compute_fx_pt(alpha[i],beta[i], detg_i, pres_i, vx[i], 
		       Q[i+of][0], Q[i+of][4], Q[i+of][1], Q[i+of][2], Q[i+of][3],
		       &(F[i+of][0]), &(F[i+of][4]), &(F[i+of][1]), &(F[i+of][2]), &(F[i+of][3]));
    
    // cure ptwise fluxes
    // ...
    
#if (!(HO_TREAT_AVGS==0))
    
    // Store arrays for averages
    // needed ONLY IF opt 1. or 2. below active !
    detga[i+of] = detg_i;
    gupxx[i+of] = gupxx_i;
    vlowx[i+of] = vlowx_i;
    vlowy[i+of] = vlowy_i;
    vlowz[i+of] = vlowz_i;
    v2   [i+of] = v2_i;
    Wlor [i+of] = Wlor_i;
    h    [i+of] = ONE + epsl[i] + pres_i/rho[i];
    kappa[i+of] = kappa_i;
    cs2  [i+of] = cs2_i;

#endif


  } // end loop i
 
 
  // ----------------------------------------------------------
  // step 2: compute numerical fluxes at interfaces 
  // ----------------------------------------------------------

  for( i = 0 ; i <= nx ; i++ ) {


if ((0.5*(rho[i] + rho[i+1])> rec_switch) && (0.5*(alpha[i] + alpha[i+1]) > alpha_switch)) {

    
    // if pt and NN = atm, set flx = 0, loop
    if ((MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) { 
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i] = ZERO;
      continue;
    }
        
    // compute eigenvectors/eigenvalues at interfaces i+1/2 
    // from Roe avgs
    
    // averages
    rho_avg   = HALF * ( rho  [i] + rho  [i+1] );    
    epsl_avg  = HALF * ( epsl [i] + epsl [i+1] );
    vx_avg    = HALF * ( vx   [i] + vx   [i+1] );
    vy_avg    = HALF * ( vy   [i] + vy   [i+1] );
    vz_avg    = HALF * ( vz   [i] + vz   [i+1] );
    
    alpha_avg = HALF * ( alpha[i] + alpha[i+1] );
    betax_avg = HALF * ( beta [i] + beta [i+1] );
    gxx_avg   = HALF * ( gxx  [i] + gxx  [i+1] );
    gxy_avg   = HALF * ( gxy  [i] + gxy  [i+1] );
    gxz_avg   = HALF * ( gxz  [i] + gxz  [i+1] );
    gyy_avg   = HALF * ( gyy  [i] + gyy  [i+1] );
    gyz_avg   = HALF * ( gyz  [i] + gyz  [i+1] );
    gzz_avg   = HALF * ( gzz  [i] + gzz  [i+1] );
    
    // additional averages are needed   
    // PROBLEM: averages of nonlinear function are different from 
    // nonlinear functions of averages ... what to use ? 
    // Note everything is mixed: metric, primitives, eigenvalues, eigenvectors, 
    // inconsistencies arise !
    
#if (HO_TREAT_AVGS==0)
    
    // opt 0. compute everything else from above avgs 

    // inverse metric and detg
    grhd_compute_detg_invg_pt(gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
			      &detg_avg, 
			      &gupxx_avg,&tmp1,&tmp2,&gupyy_avg,&tmp3,&gupzz_avg);
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx_avg,vy_avg,vz_avg, 
		       gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
		       &vlowx_avg,&vlowy_avg,&vlowz_avg, &v2_avg );

    // EoS
    //GRHD.use_eos(&pres_avg, &cs2_avg,&kappa_avg,&chi_avg, rho_avg,epsl_avg);
    EOS.comp("re","","","pc","re","", rho_avg,epsl_avg, &pres_avg,&cs2_avg,&chi_avg,&kappa_avg);
    
    
    h_avg     = ONE + epsl_avg + pres_avg/rho_avg;
    
    //if ((MATTER_USEMASK) && ((mask[i]>0.4) || (mask[i+1]>0.4))) {
    //  // use bounded vals
    //  pres_avg  = HALF * ( pres [i]    + pres [i+1] );
    //  h_avg     = HALF * ( h    [i+of] + h    [i+1+of] );      
    //  kappa_avg = HALF * ( kappa[i+of] + kappa[i+1+of] );
    //  cs2_avg   = HALF * ( cs2  [i+of] + cs2  [i+1+of] );      
    //}

    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );        
    
#elif (HO_TREAT_AVGS==1)
    
    // opt 1. average only gxxup and detg, others consistently from avgs:
    gupxx_avg = HALF * ( gupxx[i+of] + gupxx[i+1+of] ); 
    detg_avg  = HALF * ( detga[i+of] + detga[i+1+of] ); 
    
    GRHD.use_eos(&pres_avg, &cs2_avg,&kappa_avg,&chi_avg, rho_avg,epsl_avg);
    h_avg     = ONE + epsl_avg + pres_avg/rho_avg;
    
    // v_i and v^2 
    grhd_compute_v2_pt(vx_avg,vy_avg,vz_avg, gxx_avg,gxy_avg,gxz_avg,gyy_avg,gyz_avg,gzz_avg,
		       &vlowx_avg,&vlowy_avg,&vlowz_avg, &v2_avg );
    
    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );

#elif (HO_TREAT_AVGS==2)
    
    // opt 2. average all:
    gupxx_avg = HALF * ( gupxx[i+of] + gupxx[i+1+of] ); 
    detg_avg  = HALF * ( detga[i+of] + detga[i+1+of] ); 
    
    pres_avg  = HALF * ( pres [i]    + pres [i+1]    ); 
    h_avg     = HALF * ( h    [i+of] + h    [i+1+of] );      
    kappa_avg = HALF * ( kappa[i+of] + kappa[i+1+of] );
    cs2_avg   = HALF * ( cs2  [i+of] + cs2  [i+1+of] );
    
    v2_avg    = HALF * ( v2   [i+of] + v2   [i+1+of] );
    vlowx_avg = HALF * ( vlowx[i+of] + vlowx[i+1+of] );
    vlowy_avg = HALF * ( vlowy[i+of] + vlowy[i+1+of] );
    vlowz_avg = HALF * ( vlowz[i+of] + vlowz[i+1+of] );
    Wlor_avg  = HALF * ( Wlor [i+of] + Wlor [i+1+of] );

    v2_avg = DMIN(v2_avg,GRHD.HRSC_VMAX); 
    if (v2_avg==GRHD.HRSC_VMAX) Wlor_avg = GRHD.HRSC_WlorMAX;
    else Wlor_avg  = ONE/sqrt( ONE - v2_avg );
    
#else
    errorexit("you must choose a way to treat avg in HO flux!");
#endif
          
    // eigenvalues
    grhd_eig( alpha_avg, betax_avg, vx_avg, v2_avg, gupxx_avg, cs2_avg,
	      &lam[0], &lam[3], &lam[4] );
    lam[1] = lam[2] = lam[0];
    
    
    // eigenvectors     
    grhd_eigLR( gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg, 
		alpha_avg, betax_avg,  gupxx_avg, detg_avg, 
		h_avg, Wlor_avg, rho_avg, cs2_avg, kappa_avg, 
		vx_avg, vy_avg, vz_avg, vlowx_avg, vlowy_avg,  vlowz_avg,
		lam[0],  lam[3],  lam[4], 
		L, R );
        
    // compute fluxes at i+1/2      
    for( k = 0 ; k < ne ; k++ ) {          

      // flux contribute from kth char var
      
      // speed
#if (HO_TREAT_SPEED==0)
      // local flux splitting      
      a = ZERO;
      for (l = i-S; l <= i+S; l++) a = DMAX(a,amaxl[l+of][k]); 
#elif (HO_TREAT_SPEED==1)
      // global f-splitting
      a = amaxg[k];        
#else
      a = ONE; 
#endif
      
      //a = disschi * a; // add dissipation

      
      // compute pts needed in the rec (stencil S)
      for (l = i-S; l <= i+S; l++) {	
	fpk[l] = fmk[l] = ZERO;	
	for( j = 0 ; j < ne ; j++ ) {
	  fpk[l] += HALF * L[k*ne+j] * (F[l      +of][j] + a*Q[l      +of][j]);
	  fmk[l] += HALF * L[k*ne+j] * (F[2*i-l+1+of][j] - a*Q[2*i-l+1+of][j]);
	}	
      }    
      
      // reconstruct interfaces values of kth char var
      fterm[k] = MATTER.rec1dl(fpk,i) + MATTER.rec1dl(fmk,i); 

      //if (!finite(fterm[k])) fterm[k] = 0.;

    }
        
    for( j = 0 ; j < ne ; j++ ) {

      // add kth contributes to fluxes 
      fsum[j] = ZERO;   
      for( k = 0 ; k < ne ; k++ ) fsum[j] += fterm[k]*R[k*ne+j];

    }


    // store fluxes at interface i+1/2, 
    // (back to original ordering)
    flux1d[0][i] =  fsum[0];
    flux1d[1][i] =  fsum[4];
    flux1d[2][i] =  fsum[1];
    flux1d[3][i] =  fsum[2];
    flux1d[4][i] =  fsum[3];


    // cure fluxes ?
    if ((!finite(fsum[0])) || 
	(!finite(fsum[1])) || 
	(!finite(fsum[2])) ||
	(!finite(fsum[3])) ||
	(!finite(fsum[4])))
      fsum[0] = fsum[1] = fsum[2] = fsum[3] = fsum[4] = ZERO;


    if (PR) if (CheckForNANandINF(5, flux1d[0][i],flux1d[1][i],flux1d[2][i],flux1d[3][i],flux1d[4][i])) {	

	// redo calculations and display:
	printf("problem inside ho flux\n");

	printf("  var \t F_i \n");
	for (k=0; k<ne; k++) printf(" %d \t %e\n",k,flux1d[k][i]);
	
	printf("----- L\n");
	for( k = 0 ; k < ne ; k++ ) {
	  printf(" k=%d\n",k);
	  for( j = 0 ; j < ne ; j++ ) printf(" %e", L[k*ne + j]);
	  printf("\n"); 
	}
	printf("-----\n");
	printf("----- R\n");
	for( k = 0 ; k < ne ; k++ ) {
	  printf(" k=%d\n",k);
	  for( j = 0 ; j < ne ; j++ ) printf(" %e", R[k*ne + j]);
	  printf("\n"); 
	}
	printf("-----\n");

	printf(" i=%d sound speed avg = %e\n",0,cs2_avg);
	

	for( k = 0 ; k < ne ; k++ ) {
	  
	  a = ZERO;
	  for (l = i-S; l <= i+S; l++) a = DMAX(a,amaxl[l+of][k]); 	  
	  printf(" |==============> char var no k=%d \t speed a=%e \n",k,a);
	  
	  for (l = i-S; l <= i+S; l++) {
	    for( j = 0 ; j < ne ; j++ ) {
	      printf(" i=%d \t j=%d F = %e  F = %e\n",i,j,F[l      +of][j],F[2*i-l+1+of][j]);
              printf(" i=%d \t j=%d F = %e  F = %e\n",i,j,Q[l      +of][j],Q[2*i-l+1+of][j]);
	      fpk[l] += HALF * L[k*ne+j] * (F[l      +of][j] + a*Q[l      +of][j]);
	      fmk[l] += HALF * L[k*ne+j] * (F[2*i-l+1+of][j] - a*Q[2*i-l+1+of][j]);	   
	    }
            printf(" i=%d \t j=%d fp = %e  fm = %e\n",i,j,fpk[l],fmk[l]);
	  }
	  
	  printf(" |==============> char var no k=%d rec: \n",k);
          printf(" i=%d \t j=%d rec fp = %e  rec fm = %e\n",i,j,MATTER.rec1dl(fpk,i),MATTER.rec1dl(fmk,i));
	  
	}
	
	errorexits(" (PR=%d) I stop the code",PR);
	
      }
    } //switch HO simple LLF    
else { 

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


    // eigenvalues lam^0, lam^-, lam^+
    grhd_eig(alphaplus, betaplus, vxplus,v2_L,gupxxplus, cs2_L,
	     &(lamL[1]),&(lamL[0]),&(lamL[4]));
    lamL[2] = lamL[3] = lamL[1];

    grhd_eig(alphamins, betamins, vxmins,v2_R,gupxxmins, cs2_R,
             &(lamR[1]),&(lamR[0]),&(lamR[4]));
    lamR[2] = lamR[3] = lamR[1];

    // numerical fluxes 
    hrsc_numflx1d_llf(qL,qR, fL,fR, lamL,lamR, f_i, MATTER.NVf);
//    MATTER.hrsc_riemsol1d(qL,qR, fL,fR, lamL,lamR, f_i, MATTER.NVf);

    for (k=0; k<MATTER.NVf; k++) flux1d[k][i]=f_i[k];
    
    if (EOS.COLD) 
      // set Tau fluxes = 0 
      // (sources are also = 0 and Tau is updated in c2p)
      flux1d[1][i]=ZERO;
      
  
#if (SetFlxToZero)
    // test
    for (k=0; k<MATTER.NVf; k++) flux1d[k][i]=ZERO;
#endif

      if ((v2_L==GRHD.HRSC_VMAX) || (v2_R==GRHD.HRSC_VMAX)) {
	printf("               %e %e   (%d %d)\n", Wlor_L,Wlor_R,i,nx);
	printf(" flx1d max speed reached -> corrected L=%e R=%e -> %e\n", v2_L,v2_R,GRHD.HRSC_VMAX);	  
      }



    }
  } // end loop i+1/2
  

#if (!(HO_TREAT_AVGS==0))

  free(detga); free(gupxx); 
  free(vlowx); free(vlowy); free(vlowz);
  free(v2); free(Wlow);
  free(h); free(kappa); free(cs2);

#endif

  free(lam); free(amaxg); 
  free(L); free(R); 
  free(buffer); 
  free(fsum); 
  free(fterm); 
  ddfree(amaxl,ntot); 
  ddfree(Q,ntot); 
  ddfree(F,ntot); 

    free(qL); free(qR); 
    free(fL); free(fR);
    free(f_i); free(lamL); free(lamR); 

} 

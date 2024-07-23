/* grhd_flx1d_fvpr.c 
   sbernuz 0/2012
   fschianchi 02/2023 */

#include "bam.h"
#include "grhd.h"


#define PR 0
#define CureFlxVmax 1

#define HO_TREAT_AVGS 0  // opts : 0, 1, 2
#define HO_TREAT_SPEED 0 // opts : 0, 1, 

#define ZERO 0.
#define HALF 0.5
#define ONE 1.0
#define TINY 1e-30
#define mach_limit GRHD.HLLC_mach_lim
#define K GRHD.HLLC_K

// for debug
#define DEBUG 0
#define SetFlxToZero 0

#define COMPUTE_TETRAD_COMPONENTS   B = 1/sqrt(gupxx); \
                                    C = 1/sqrt(gzz); \
                                    D = 1/sqrt(gzz*(gyy*gzz-gyz*gyz)); \
                                    ett = 1/alpha; \
                                    etx = -betax/alpha; exx = B*gupxx; \
                                    ety = -betay/alpha; exy = B*gupxy; eyy = D*gzz; \
                                    etz = -betaz/alpha; exz = B*gupxz; eyz = -D*gyz; ezz = C; \
                                    dtt = -alpha; dxt = B*betax; dyt = D*(betay_low*gzz-betaz_low*gyz); dzt = C*betaz_low; \
                                                  dxx = B;       dyx = D*(gxy*gzz-gxz*gyz);             dzx = C*gxz; \
                                                                 dyy = D*(gyy*gzz-gyz*gyz);             dzy = C*gyz; \
                                                                                                        dzz = C*gzz; \

void HLLC_central_state(double D, double Tau, double Sx, double Sy, double Sz, 
                        double lambda, double vx, double p, double lambda_c, double Pc, 
                        double *Dc, double *Tauc, double *Sxc, double *Syc, double *Szc) {

  *Dc   =  D  *(lambda - vx)                      /(lambda - lambda_c);
  *Tauc = (Tau*(lambda - vx) + Pc*lambda_c - p*vx)/(lambda - lambda_c);

  *Sxc = (Sx*(lambda - vx) + Pc - p)/(lambda - lambda_c);
  *Syc =  Sy*(lambda - vx)          /(lambda - lambda_c);
  *Szc =  Sz*(lambda - vx)          /(lambda - lambda_c);

}

void HLLC_sided_solver(double vx_L, double vx_R,
                       double pres_L, double pres_R,
                       double *qL, double *qR, 
		                   double *fL, double *fR, 
		                   double *lamL, double *lamR, double v_interface,
		                   double *f, double *q, int nf) {

  double lambda_L,lambda_R,lambda_c,Pc;
  double AL,AR,BL,BR;
  double a,b,c,Delta;

  lambda_L = DMIN(lamL[0],lamR[0]);
  lambda_R = DMAX(lamL[4],lamR[4]);

  if      (lambda_L>=v_interface) for (int k=0; k<MATTER.NVf; k++) {f[k]=fL[k]; q[k]=qL[k];}
  else if (lambda_R<=v_interface) for (int k=0; k<MATTER.NVf; k++) {f[k]=fR[k]; q[k]=qR[k];}

  else { // reconstruct the intermediate state if necessary

    AL = lambda_L*(qL[0]+qL[1]) - qL[2];
    AR = lambda_R*(qR[0]+qR[1]) - qR[2];

    BL = qL[2]*(lambda_L - vx_L) - pres_L;
    BR = qR[2]*(lambda_R - vx_R) - pres_R;

    a = lambda_L*AR - lambda_R*AL;
    b = lambda_R*BL + AL - lambda_L*BR - AR;
    c = BR - BL;

    Delta = b*b - 4*a*c;

    if      (a!=0) lambda_c = 0.5*(-b - sqrt(Delta))/a;
    else if (b!=0) lambda_c = -c/b;
    else           lambda_c = 0.5*(lambda_R+lambda_L); 

    Pc = (lambda_c*(AL+AR)-BR-BL)/(2-lambda_c*(lambda_R+lambda_L));

    if (v_interface<=lambda_c) {
      HLLC_central_state(qL[0], qL[1], qL[2], qL[3], qL[4], 
                          lambda_L, vx_L, pres_L, lambda_c, Pc,
                          &(q[0]), &(q[1]), &(q[2]), &(q[3]), &(q[4]));
      for (int k=0; k<MATTER.NVf; k++) f[k] = fL[k] + lambda_L*(q[k] - qL[k]);        
    }

    else if (v_interface>lambda_c) {
      HLLC_central_state(qR[0], qR[1], qR[2], qR[3], qR[4], 
                          lambda_R, vx_R, pres_R, lambda_c, Pc, 
                          &(q[0]), &(q[1]), &(q[2]), &(q[3]), &(q[4]));
      for (int k=0; k<MATTER.NVf; k++) f[k] = fR[k] + lambda_R*(q[k] - qR[k]);        
    }
        
  }

}

void HLLC_centered_solver(double vx_L, double vx_R, double v_interface, 
                          double pres_L, double pres_R,
                          double cs2_L,  double cs2_R,
                          double *qL, double *qR, 
		                      double *fL, double *fR, 
		                      double *lamL, double *lamR,
		                      double *f, double *q, int nf) {

  double lambda_L,lambda_R,lambda_c,Pc;
  double AL,AR,BL,BR;
  double a,b,c,Delta;

  double *qcL = dmalloc(nf);
  double *qcR = dmalloc(nf);

  lambda_L = DMIN(lamL[0],lamR[0]);
  lambda_R = DMAX(lamL[4],lamR[4]);

  if      (lambda_L>=v_interface) for (int k=0; k<MATTER.NVf; k++) {f[k]=fL[k]; q[k]=qL[k];}
  else if (lambda_R<=v_interface) for (int k=0; k<MATTER.NVf; k++) {f[k]=fR[k]; q[k]=qR[k];}

  else { // reconstruct the intermediate state if necessary

    AL = lambda_L*(qL[0]+qL[1]) - qL[2];
    AR = lambda_R*(qR[0]+qR[1]) - qR[2];

    BL = qL[2]*(lambda_L - vx_L) - pres_L;
    BR = qR[2]*(lambda_R - vx_R) - pres_R;

    a = lambda_L*AR - lambda_R*AL;
    b = lambda_R*BL + AL - lambda_L*BR - AR;
    c = BR - BL;

    Delta = b*b - 4*a*c;

    if      (a!=0) lambda_c = 0.5*(-b - sqrt(Delta))/a;
    else if (b!=0) lambda_c = -c/b;
    else           lambda_c = 0.5*(lambda_R+lambda_L); 

    Pc = (lambda_c*(AL+AR)-BR-BL)/(2-lambda_c*(lambda_R+lambda_L));

    HLLC_central_state(qL[0], qL[1], qL[2], qL[3], qL[4], 
                        lambda_L, vx_L, pres_L, lambda_c, Pc,
                        &(qcL[0]), &(qcL[1]), &(qcL[2]), &(qcL[3]), &(qcL[4]));


    HLLC_central_state(qR[0], qR[1], qR[2], qR[3], qR[4], 
                        lambda_R, vx_R, pres_R, lambda_c, Pc, 
                        &(qcR[0]), &(qcR[1]), &(qcR[2]), &(qcR[3]), &(qcR[4]));
        
    if (v_interface<=lambda_c) {
      for (int k=0; k<nf; k++) {
       f[k] = HALF*(fL[k]+fR[k]) + HALF*(lambda_L*(qcL[k]-qL[k]) + lambda_c*(qcL[k]-qcR[k]) + lambda_R*(qcR[k]-qR[k]));  
       q[k] = qcL[k];      
      }
    }

    else if (v_interface>lambda_c) {
      for (int k=0; k<nf; k++) {
       f[k] = HALF*(fL[k]+fR[k]) + HALF*(lambda_L*(qcL[k]-qL[k]) - lambda_c*(qcL[k]-qcR[k]) + lambda_R*(qcR[k]-qR[k]));  
       q[k] = qcR[k];      
      }    
    }

  }

  free(qcL); free(qcR);

}

void HLLC_LM_solver(double vx_L, double vx_R, double v_interface, 
                    double pres_L, double pres_R,
                    double cs2_L,  double cs2_R,
                    double *qL, double *qR, 
		                double *fL, double *fR, 
		                double *lamL, double *lamR,
		                double *f, double *q, int nf) {

  double lambda_L,lambda_R,lambda_c,Pc;
  double AL,AR,BL,BR;
  double a,b,c,Delta;
  double mach_local,phi;

  double *qcL = dmalloc(nf);
  double *qcR = dmalloc(nf);

  lambda_L = DMIN(lamL[0],lamR[0]);
  lambda_R = DMAX(lamL[4],lamR[4]);

  if      (lambda_L>=v_interface) for (int k=0; k<MATTER.NVf; k++) {f[k]=fL[k]; q[k]=qL[k];}
  else if (lambda_R<=v_interface) for (int k=0; k<MATTER.NVf; k++) {f[k]=fR[k]; q[k]=qR[k];}

  else { // reconstruct the intermediate state if necessary

    AL = lambda_L*(qL[0]+qL[1]) - qL[2];
    AR = lambda_R*(qR[0]+qR[1]) - qR[2];

    BL = qL[2]*(lambda_L - vx_L) - pres_L;
    BR = qR[2]*(lambda_R - vx_R) - pres_R;

    a = lambda_L*AR - lambda_R*AL;
    b = lambda_R*BL + AL - lambda_L*BR - AR;
    c = BR - BL;

    Delta = b*b - 4*a*c;

    if      (a!=0) lambda_c = 0.5*(-b - sqrt(Delta))/a;
    else if (b!=0) lambda_c = -c/b;
    else           lambda_c = 0.5*(lambda_R+lambda_L); 

    Pc = (lambda_c*(AL+AR)-BR-BL)/(2-lambda_c*(lambda_R+lambda_L));

    HLLC_central_state(qL[0], qL[1], qL[2], qL[3], qL[4], 
                        lambda_L, vx_L, pres_L, lambda_c, Pc,
                        &(qcL[0]), &(qcL[1]), &(qcL[2]), &(qcL[3]), &(qcL[4]));


    HLLC_central_state(qR[0], qR[1], qR[2], qR[3], qR[4], 
                        lambda_R, vx_R, pres_R, lambda_c, Pc, 
                        &(qcR[0]), &(qcR[1]), &(qcR[2]), &(qcR[3]), &(qcR[4]));

    mach_local = DMAX(fabs(lambda_c/lambda_R),fabs(lambda_c/lambda_L));  
    phi = sin(DMIN(1,mach_local/mach_limit)*PI/2);

    if (v_interface<=lambda_c) {
      for (int k=0; k<nf; k++) {
       f[k] = HALF*(fL[k]+fR[k]) + HALF*(phi*lambda_L*(qcL[k]-qL[k]) + lambda_c*(qcL[k]-qcR[k]) + phi*lambda_R*(qcR[k]-qR[k]));  
       q[k] = qcL[k];      
      }
    }

    else if (v_interface>lambda_c) {
      for (int k=0; k<nf; k++) {
       f[k] = HALF*(fL[k]+fR[k]) + HALF*(phi*lambda_L*(qcL[k]-qL[k]) - lambda_c*(qcL[k]-qcR[k]) + phi*lambda_R*(qcR[k]-qR[k]));  
       q[k] = qcR[k];      
      }    
    }

  }

  free(qcL); free(qcR);

}

// HLL flux and central state
void HLL_solver(double *uL, double *uR, 
		            double *fL, double *fR, 
		            double *lamL, double *lamR, double v_interface,
		            double *f, double *q, int nf) {

  double lambda_L = DMIN(lamL[0],lamR[0]);
  double lambda_R = DMAX(lamL[4],lamR[4]);
  
  double oda = ONE/(lambda_R - lambda_L + TINY);
  double apm = lambda_R * lambda_L;

  if      (lambda_L>=v_interface) for (int v=0; v<nf; v++) {f[v]=fL[v]; q[v]=uL[v];}
  else if (lambda_R<=v_interface) for (int v=0; v<nf; v++) {f[v]=fR[v]; q[v]=uR[v];}
  else for (int v=0; v<nf; v++) {
    f[v] = oda * ( lambda_R * fL[v] - lambda_L * fR[v] + apm * (uR[v] - uL[v]) );
    q[v] = oda * ( lambda_R * uR[v] - lambda_L * uL[v] +       (fL[v] - fR[v]) );
  }

}


// 2nd order finite volume fluxes based on HLLC Riemann solver in local Minkowski (K. Kiuchi 2022 https://arxiv.org/abs/2205.04487)
// and primitive reconstruction
void grhd_flx1d_fvpr_HLLC(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx, int ngx, int direction, double **flux1d)
{  
  double dummy;

  double rhomins, epslmins, vxmins, vymins, vzmins, Gammamins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus, Gammaplus;
  
  double gxx, gxy, gxz, gyy, gyz, gzz;
  double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
  double betax, betay, betaz, betax_low, betay_low, betaz_low;
  double alpha, detg_i, sqrtdetg;

  double pres_L, v2_L, Wlor_L, cs2_L, vxtet_L, vytet_L, vztet_L;
  double pres_R, v2_R, Wlor_R, cs2_R, vxtet_R, vytet_R, vztet_R;

  double mach_local, phi, contact_indicator;

  double B,C,D;

  // tetrad elements with upper index. first index i tetrad element second one is spacetime index
  double ett, etx, ety, etz;
  double exx, exy, exz;
  double eyy, eyz;
  double ezz;  

  // tetrad elements with lower index. first index i tetrad element second one is spacetime index
  double dtt;
  double dxt, dxx;
  double dyt, dyx, dyy;
  double dzt, dzx, dzy, dzz;

  // Rieman solver vars
  double lambda_L, lambda_R, lambda_c;
  double Pc, v_interface, Delta;
  double AL, AR, BL, BR;
  double a, b, c;

  double *qL     = dmalloc(MATTER.NVq);
  double *qR     = dmalloc(MATTER.NVq);  
  double *qc     = dmalloc(MATTER.NVq);
  double *fL     = dmalloc(MATTER.NVf);  
  double *fR     = dmalloc(MATTER.NVf);  
  double *f_i    = dmalloc(MATTER.NVf);  
  double *f_HLL  = dmalloc(MATTER.NVf);  
  double *q_HLL  = dmalloc(MATTER.NVf);  
  double *f_HLLC = dmalloc(MATTER.NVf);  
  double *q_HLLC = dmalloc(MATTER.NVf);
  double *lamL   = dmalloc(MATTER.NVf);
  double *lamR   = dmalloc(MATTER.NVf); 

  int i,k;

  for( i = 0 ; i <= nx ; i++ ) {
    
    // if the point is atm, do nothing EXCEPT
    // set the fluxes to 0 to be sure no random numbers are added
    if ( (MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i]= ZERO;
      continue;
    }
    
    // rec metric     
    alpha = MATTER.rec1dml( g1d[0], i );
    betax = MATTER.rec1dml( g1d[1], i );
    betay = MATTER.rec1dml( g1d[2], i );
    betaz = MATTER.rec1dml( g1d[3], i );
    gxx   = MATTER.rec1dml( g1d[4], i );
    gxy   = MATTER.rec1dml( g1d[5], i );
    gxz   = MATTER.rec1dml( g1d[6], i );
    gyy   = MATTER.rec1dml( g1d[7], i );
    gyz   = MATTER.rec1dml( g1d[8], i );
    gzz   = MATTER.rec1dml( g1d[9], i );

    // left state, u_L =  u(i)^+, rec prim     
    rhoplus  = MATTER.rec1dl( w1d[0], i );
    epslplus = MATTER.rec1dl( w1d[1], i );
    vxplus   = MATTER.rec1dl( w1d[2], i );
    vyplus   = MATTER.rec1dl( w1d[3], i );
    vzplus   = MATTER.rec1dl( w1d[4], i );

    if (rhoplus<GRHD.ATM_RHOATM)   rhoplus   = MATTER.rec1dsl( w1d[0], i );
    if (epslplus<GRHD.ATM_EPSLATM) epslplus  = MATTER.rec1dsl( w1d[1], i );

    // left state, u_R =  u(i)^-, rec prim     
    rhomins  = MATTER.rec1dr( w1d[0], i+1 );
    epslmins = MATTER.rec1dr( w1d[1], i+1 );
    vxmins   = MATTER.rec1dr( w1d[2], i+1 );
    vymins   = MATTER.rec1dr( w1d[3], i+1 );
    vzmins   = MATTER.rec1dr( w1d[4], i+1 );

    if (rhomins<GRHD.ATM_RHOATM)   rhomins   = MATTER.rec1dsr( w1d[0], i+1 );
    if (epslmins<GRHD.ATM_EPSLATM) epslmins  = MATTER.rec1dsr( w1d[1], i+1 );

    // eos
    EOS.comp("re","","","pc","","", rhoplus,epslplus, &pres_L,&cs2_L);
    EOS.comp("re","","","pc","","", rhomins,epslmins, &pres_R,&cs2_R);
    Gammaplus = EOS_Gamma(rhoplus, pres_L, epslplus, cs2_L);
    Gammamins = EOS_Gamma(rhomins, pres_R, epslmins, cs2_R);
    
    // detg and g^-1
    grhd_compute_detg_invg_pt(gxx,gxy,gxz,gyy,gyz,gzz,
                            &detg_i, 
                            &gupxx,&gupxy,&gupxz,&gupyy,&gupyz,&gupzz);
    sqrtdetg = sqrt(detg_i);

    //beta low
    grhd_compute_v2_pt(betax,betay,betaz, gxx,gxy,gxz,gyy,gyz,gzz,
		                  &betax_low,&betay_low,&betaz_low, &dummy);

    COMPUTE_TETRAD_COMPONENTS
    
    // compute velocity in tetrad frame
    vxtet_L = dxx*vxplus;
    vytet_L = dyx*vxplus + dyy*vyplus;
    vztet_L = dzx*vxplus + dzy*vyplus + dzz*vzplus;

    vxtet_R = dxx*vxmins;
    vytet_R = dyx*vxmins + dyy*vymins;
    vztet_R = dzx*vxmins + dzy*vymins + dzz*vzmins;

    v2_L = pow(vxtet_L,2) + pow(vytet_L,2) + pow(vztet_L,2);
    v2_R = pow(vxtet_R,2) + pow(vytet_R,2) + pow(vztet_R,2);

    v2_L = DMIN(v2_L, GRHD.HRSC_VMAX); // needed ? testme !
    v2_R = DMIN(v2_R, GRHD.HRSC_VMAX);

    Wlor_L = ONE / sqrt(ONE - v2_L);
    Wlor_R = ONE / sqrt(ONE - v2_R);
    
    // compute conservative vars in tetrad frame
    // gamma_ii = alpha = detg = 1 and beta^i = 0 cause we are in the tetrad coordinrates now
    grhd_compute_q_pt(1,0,0,1,0,1,1, 
		      rhoplus,epslplus,pres_L,vxtet_L,vytet_L,vztet_L,Wlor_L,
		      &(qL[0]), &(qL[1]), &(qL[2]), &(qL[3]), &(qL[4]));

    grhd_compute_q_pt(1,0,0,1,0,1,1,
          rhomins,epslmins,pres_R,vxtet_R,vytet_R,vztet_R,Wlor_R,
          &(qR[0]), &(qR[1]), &(qR[2]), &(qR[3]), &(qR[4]));

    // compute L/R fluxes in tetrad frame
    // also here alpha = detg = 1 and beta^i = 0 because of the tetrad frame
    grhd_compute_fx_pt(1,0,1, 
		      pres_L, vxtet_L, 
		      qL[0], qL[1], qL[2], qL[3], qL[4],
		      &(fL[0]), &(fL[1]), &(fL[2]), &(fL[3]), &(fL[4]));
    
    grhd_compute_fx_pt(1,0,1,
          pres_R, vxtet_R,
          qR[0], qR[1], qR[2], qR[3], qR[4],
          &(fR[0]), &(fR[1]), &(fR[2]), &(fR[3]), &(fR[4]));

    // eigenvalues lam^0, lam^-, lam^+ in tetrad frame
    // again alpha = detg = 1 and beta^i = 0 because of the tetrad frame
    grhd_eig(1, 0, vxtet_L,v2_L, 1, cs2_L,
	          &(lamL[1]),&(lamL[0]),&(lamL[4]));
    lamL[2] = lamL[3] = lamL[1];

    grhd_eig(1, 0, vxtet_R,v2_R, 1, cs2_R,
             &(lamR[1]),&(lamR[0]),&(lamR[4]));
    lamR[2] = lamR[3] = lamR[1];

    // now we have to go for the actual Riemann solver
    v_interface = betax/(alpha*sqrt(gupxx));

    HLLC_sided_solver(vxtet_L, vxtet_R,
                      pres_L, pres_R,
                      qL, qR, 
                      fL, fR, 
                      lamL, lamR, v_interface,
                      f_HLLC, q_HLLC, 
                      MATTER.NVf);

    HLL_solver(qL, qR, 
		           fL, fR, 
		           lamL, lamR, v_interface,
		           f_HLL, q_HLL, 
               MATTER.NVf);

    mach_local = DMIN(fabs(vxtet_L)/sqrt(cs2_L),fabs(vxtet_R)/sqrt(cs2_R));  
    phi = sin(DMIN(1,mach_local/mach_limit)*PI/2);
    contact_indicator = K*DMIN(Gammaplus,Gammamins)*fabs(w1d[0][i+1]-w1d[0][i])/DMIN(w1d[0][i+1],w1d[0][i]) 
                                                  - fabs(w1d[5][i+1]-w1d[5][i])/DMIN(w1d[5][i+1],w1d[5][i]);

    if (contact_indicator>=0) phi=1;

    for (k=0; k<MATTER.NVf; k++) {
      f_i[k] = f_HLLC[k] + (1-phi)*(f_HLL[k]-f_HLLC[k]);
      qc[k]  = q_HLLC[k] + (1-phi)*(q_HLL[k]-q_HLLC[k]);
    }

    // finally transform the fluxes back to the original frame and save them into flux1d
    flux1d[0][i] = sqrtdetg*alpha*(etx*qc[0] + exx*f_i[0]);
    flux1d[1][i] = sqrtdetg*alpha*(etx*qc[1] + exx*f_i[1]);

    flux1d[2][i] = sqrtdetg*alpha*(etx*(dxx*qc[2]+dyx*qc[3]+dzx*qc[4]) + exx*(dxx*f_i[2]+dyx*f_i[3]+dzx*f_i[4]));    
    flux1d[3][i] = sqrtdetg*alpha*(etx*(          dyy*qc[3]+dzy*qc[4]) + exx*(           dyy*f_i[3]+dzy*f_i[4]));    
    flux1d[4][i] = sqrtdetg*alpha*(etx*(                    dzz*qc[4]) + exx*(                      dzz*f_i[4]));     



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
      if (DEBUG && direction==1) {

        printf(" ************************** i =%d \n",i);
        printf("  var \t F_i \n");
        for (k=0; k<MATTER.NVf; k++) printf(" %d \t %e \n",k,flux1d[k][i]);
        
        printf(" LEFT\n");
        printf(" q \t f \t lam \n");
        for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qL[k],fL[k],lamL[k]);
        printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e %+.6e \n", rhoplus,epslplus,vxplus,vyplus,vzplus,v2_L,pres_L,lambda_L);
        
        printf(" RIGHT\n");
        printf(" q \t f \t lam \n");
        for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qR[k],fR[k],lamR[k]);
        printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e %+.6e \n", rhomins,epslmins,vxmins,vymins,vzmins,v2_R,pres_R,lambda_R);

        printf(" CENTRAL\n");
        for (k=0; k<MATTER.NVf; k++) printf(" %+.6e \n",qc[k]);
        printf(" w = %+.6e %+.6e %+.6e %+.6e \n", lambda_c,Pc,v_interface,Delta);

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

  free(qL); free(qR); free(qc);  free(q_HLL); free(q_HLLC);
  free(fL); free(fR); free(f_i); free(f_HLL); free(f_HLLC);
  free(lamL); free(lamR); 

}

// 2nd order finite volume fluxes based on HLLC Riemann solver in local Minkowski (K. Kiuchi 2022 https://arxiv.org/abs/2205.04487)
// and primitive reconstruction
// in this version the quantity b^i := (Wlor v^i) is reconstructed instead of v^i
void grhd_flx1d_fvprb_HLLC(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx, int ngx, int direction, double **flux1d)
{ 
  double tmp1, tmp2, tmp3, dummy;

  double rhomins, epslmins, vxmins, vymins, vzmins, Gammamins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus, Gammaplus;
  
  double gxx, gxy, gxz, gyy, gyz, gzz;
  double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
  double betax, betay, betaz, betax_low, betay_low, betaz_low;
  double alpha, detg_i, sqrtdetg;

  double pres_L, v2_L, Wlor_L, cs2_L, vxtet_L, vytet_L, vztet_L;
  double pres_R, v2_R, Wlor_R, cs2_R, vxtet_R, vytet_R, vztet_R;

  double mach_local, phi, contact_indicator;

  double B,C,D;

  // tetrad elements with upper index. first index i tetrad element second one is spacetime index
  double ett, etx, ety, etz;
  double exx, exy, exz;
  double eyy, eyz;
  double ezz;  

  // tetrad elements with lower index. first index i tetrad element second one is spacetime index
  double dtt;
  double dxt, dxx;
  double dyt, dyx, dyy;
  double dzt, dzx, dzy, dzz;

  // Rieman solver vars
  double lambda_L, lambda_R, lambda_c;
  double Pc, v_interface, Delta;
  double AL, AR, BL, BR;
  double a, b, c;

  double *qL     = dmalloc(MATTER.NVq);
  double *qR     = dmalloc(MATTER.NVq);  
  double *qc     = dmalloc(MATTER.NVq);
  double *fL     = dmalloc(MATTER.NVf);  
  double *fR     = dmalloc(MATTER.NVf);  
  double *f_i    = dmalloc(MATTER.NVf);  
  double *f_HLL  = dmalloc(MATTER.NVf);  
  double *q_HLL  = dmalloc(MATTER.NVf);  
  double *f_HLLC = dmalloc(MATTER.NVf);  
  double *q_HLLC = dmalloc(MATTER.NVf);
  double *lamL   = dmalloc(MATTER.NVf);
  double *lamR   = dmalloc(MATTER.NVf); 

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
  for( i = -ngx+1 ; i <= (nx+ngx) ; i++ ) {

    // use "plus/L" as tmp vars
    vxplus = w1d[2][i]; // v^i 
    vyplus = w1d[3][i];
    vzplus = w1d[4][i];

    gxx = g1d[4][i]; // g_{ij}
    gxy = g1d[5][i];
    gxz = g1d[6][i];
    gyy = g1d[7][i];
    gyz = g1d[8][i];
    gzz = g1d[9][i];

    grhd_compute_v2_pt(vxplus,vyplus,vzplus, gxx,gxy,gxz,gyy,gyz,gzz,
		                  &tmp1,&tmp2,&tmp3, &v2_L);
    v2_L   = DMIN( v2_L, GRHD.HRSC_VMAX);
    Wlor_L = ONE/sqrt( ONE - v2_L );
    
    *(buf[0]+off +i) = Wlor_L * vxplus;
    *(buf[1]+off +i) = Wlor_L * vyplus;
    *(buf[2]+off +i) = Wlor_L * vzplus;

  }

  bx = buf[0] + off;
  by = buf[1] + off;
  bz = buf[2] + off;


  for( i = 0 ; i <= nx ; i++ ) {
    
    // if the point is atm, do nothing EXCEPT
    // set the fluxes to 0 to be sure no random numbers are added
    if ( (MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i]= ZERO;
      continue;
    }

    // left state, u_L =  u(i)^+, rec prim     
    // rec metric     
    alpha = MATTER.rec1dml( g1d[0], i );
    betax = MATTER.rec1dml( g1d[1], i );
    betay = MATTER.rec1dml( g1d[2], i );
    betaz = MATTER.rec1dml( g1d[3], i );
    gxx   = MATTER.rec1dml( g1d[4], i );
    gxy   = MATTER.rec1dml( g1d[5], i );
    gxz   = MATTER.rec1dml( g1d[6], i );
    gyy   = MATTER.rec1dml( g1d[7], i );
    gyz   = MATTER.rec1dml( g1d[8], i );
    gzz   = MATTER.rec1dml( g1d[9], i );

    rhoplus  = MATTER.rec1dl( w1d[0], i );
    epslplus = MATTER.rec1dl( w1d[1], i );
    bxplus   = MATTER.rec1dl( bx, i ); // note now this is bplus^i
    byplus   = MATTER.rec1dl( by, i );
    bzplus   = MATTER.rec1dl( bz, i );

    if (rhoplus<GRHD.ATM_RHOATM)   rhoplus   = MATTER.rec1dsl( w1d[0], i );
    if (epslplus<GRHD.ATM_EPSLATM) epslplus  = MATTER.rec1dsl( w1d[1], i );

    rhomins  = MATTER.rec1dr( w1d[0], i+1 );
    epslmins = MATTER.rec1dr( w1d[1], i+1 );
    bxmins   = MATTER.rec1dr( bx, i+1 ); // note now this is bmins^i
    bymins   = MATTER.rec1dr( by, i+1 );
    bzmins   = MATTER.rec1dr( bz, i+1 );

    if (rhomins<GRHD.ATM_RHOATM)   rhomins   = MATTER.rec1dsr( w1d[0], i+1 );
    if (epslmins<GRHD.ATM_EPSLATM) epslmins  = MATTER.rec1dsr( w1d[1], i+1 );

    // eos
    EOS.comp("re","","","pc","","", rhoplus,epslplus, &pres_L,&cs2_L);
    EOS.comp("re","","","pc","","", rhomins,epslmins, &pres_R,&cs2_R);
    Gammaplus = EOS_Gamma(rhoplus, pres_L, epslplus, cs2_L);
    Gammamins = EOS_Gamma(rhomins, pres_R, epslmins, cs2_R);
    
    // detg and g^-1
    grhd_compute_detg_invg_pt(gxx,gxy,gxz,gyy,gyz,gzz,
                            &detg_i, 
                            &gupxx,&gupxy,&gupxz,&gupyy,&gupyz,&gupzz);
    sqrtdetg = sqrt(detg_i);

    //beta low
    grhd_compute_v2_pt(betax,betay,betaz, gxx,gxy,gxz,gyy,gyz,gzz,
		                  &betax_low,&betay_low,&betaz_low, &dummy);

    // b2 
    // note now this is b^2, and v_i is actually b_i (corrected below)
    grhd_compute_v2_pt(bxplus,byplus,bzplus, gxx,gxy,gxz,gyy,gyz,gzz, &tmp1,&tmp2,&tmp3, &b2_L);
    grhd_compute_v2_pt(bxmins,bymins,bzmins, gxx,gxy,gxz,gyy,gyz,gzz, &tmp1,&tmp2,&tmp3, &b2_R);
    
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
    
    COMPUTE_TETRAD_COMPONENTS

    // compute velocity in tetrad frame
    vxtet_L = dxx*vxplus;
    vytet_L = dyx*vxplus + dyy*vyplus;
    vztet_L = dzx*vxplus + dzy*vyplus + dzz*vzplus;

    vxtet_R = dxx*vxmins;
    vytet_R = dyx*vxmins + dyy*vymins;
    vztet_R = dzx*vxmins + dzy*vymins + dzz*vzmins;

    v2_L = pow(vxtet_L,2) + pow(vytet_L,2) + pow(vztet_L,2);
    v2_R = pow(vxtet_R,2) + pow(vytet_R,2) + pow(vztet_R,2);

    v2_L = DMIN(v2_L, GRHD.HRSC_VMAX); // needed ? testme !
    v2_R = DMIN(v2_R, GRHD.HRSC_VMAX);

    Wlor_L = ONE / sqrt(ONE - v2_L);
    Wlor_R = ONE / sqrt(ONE - v2_R);

    #if (0) // (CureFlxVmax)
        // important hack used in bam 11                                                                                                                  
        // do we need it in this routine ? testme !    
        if ( (Wlor_L == GRHD.HRSC_WlorMAX) || (Wlor_R == GRHD.HRSC_WlorMAX) ) {
          for (k=0; k<MATTER.NVf; k++) flux1d[k][i] = ZERO;      
          if (PR) printf(" v>VMAX on both sides, cured flx1d setting them to zero\n");
          continue;
        }    
    #endif

    // compute conservative vars in tetrad frame
    // gamma_ii = alpha = detg = 1 and beta^i = 0 cause we are in the tetrad coordinrates now
    grhd_compute_q_pt(1,0,0,1,0,1,1, 
		      rhoplus,epslplus,pres_L,vxtet_L,vytet_L,vztet_L,Wlor_L,
		      &(qL[0]), &(qL[1]), &(qL[2]), &(qL[3]), &(qL[4]));

    grhd_compute_q_pt(1,0,0,1,0,1,1,
          rhomins,epslmins,pres_R,vxtet_R,vytet_R,vztet_R,Wlor_R,
          &(qR[0]), &(qR[1]), &(qR[2]), &(qR[3]), &(qR[4]));

    // compute L/R fluxes in tetrad frame
    // also here alpha = detg = 1 and beta^i = 0 because of the tetrad frame
    grhd_compute_fx_pt(1,0,1, 
		      pres_L, vxtet_L, 
		      qL[0], qL[1], qL[2], qL[3], qL[4],
		      &(fL[0]), &(fL[1]), &(fL[2]), &(fL[3]), &(fL[4]));
    
    grhd_compute_fx_pt(1,0,1,
          pres_R, vxtet_R,
          qR[0], qR[1], qR[2], qR[3], qR[4],
          &(fR[0]), &(fR[1]), &(fR[2]), &(fR[3]), &(fR[4]));

    // eigenvalues lam^0, lam^-, lam^+ in tetrad frame
    // again alpha = detg = 1 and beta^i = 0 because of the tetrad frame
    grhd_eig(1, 0, vxtet_L,v2_L, 1, cs2_L,
	          &(lamL[1]),&(lamL[0]),&(lamL[4]));
    lamL[2] = lamL[3] = lamL[1];

    grhd_eig(1, 0, vxtet_R,v2_R, 1, cs2_R,
             &(lamR[1]),&(lamR[0]),&(lamR[4]));
    lamR[2] = lamR[3] = lamR[1];

    // now we have to go for the actual Riemann solver
    v_interface = betax/(alpha*sqrt(gupxx));

    HLLC_sided_solver(vxtet_L, vxtet_R,
                      pres_L, pres_R,
                      qL, qR, 
                      fL, fR, 
                      lamL, lamR, v_interface,
                      f_HLLC, q_HLLC, 
                      MATTER.NVf);

    HLL_solver(qL, qR, 
		           fL, fR, 
		           lamL, lamR, v_interface,
		           f_HLL, q_HLL, 
               MATTER.NVf);

    mach_local = DMIN(fabs(vxtet_L)/sqrt(cs2_L),fabs(vxtet_R)/sqrt(cs2_R));  
    phi = sin(DMIN(1,mach_local/mach_limit)*PI/2);
    contact_indicator = K*DMIN(Gammaplus,Gammamins)*fabs(w1d[0][i+1]-w1d[0][i])/DMIN(w1d[0][i+1],w1d[0][i]) 
                                                  - fabs(w1d[5][i+1]-w1d[5][i])/DMIN(w1d[5][i+1],w1d[5][i]);

    if (contact_indicator>=0) phi=1;

    for (k=0; k<MATTER.NVf; k++) {
      f_i[k] = f_HLLC[k] + (1-phi)*(f_HLL[k]-f_HLLC[k]);
      qc[k]  = q_HLLC[k] + (1-phi)*(q_HLL[k]-q_HLLC[k]);
    }

    // finally transform the fluxes back to the original frame and save them into flux1d
    flux1d[0][i] = sqrtdetg*alpha*(etx*qc[0] + exx*f_i[0]);
    flux1d[1][i] = sqrtdetg*alpha*(etx*qc[1] + exx*f_i[1]);

    flux1d[2][i] = sqrtdetg*alpha*(etx*(dxx*qc[2]+dyx*qc[3]+dzx*qc[4]) + exx*(dxx*f_i[2]+dyx*f_i[3]+dzx*f_i[4]));    
    flux1d[3][i] = sqrtdetg*alpha*(etx*(          dyy*qc[3]+dzy*qc[4]) + exx*(           dyy*f_i[3]+dzy*f_i[4]));    
    flux1d[4][i] = sqrtdetg*alpha*(etx*(                    dzz*qc[4]) + exx*(                      dzz*f_i[4]));     

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
  free(qL); free(qR); free(qc);  free(q_HLL); free(q_HLLC);
  free(fL); free(fR); free(f_i); free(f_HLL); free(f_HLLC);
  free(lamL); free(lamR);  

}

// 2nd order finite volume fluxes based on HLL/LLF Riemann solver
// and PPM primitive reconstruction
void grhd_flx1d_fvpr_ppm(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx, int ngx, int direction, double **flux1d)
{

  int Ntot = nx+2*ngx;
  int off  = ngx-1;

  /* contains all informations */
  double **buffer = (double**) malloc (10*sizeof(double*));
  for (int l=0; l<10; l++) 
    buffer[l] = (double*) malloc ((Ntot)*sizeof(double));

  double *rhoL    = buffer[0] + off;
  double *rhoR    = buffer[1] + off;
  double *epslL   = buffer[2] + off;
  double *epslR   = buffer[3] + off;
  double *vxL     = buffer[4] + off;
  double *vxR     = buffer[5] + off;
  double *vyL     = buffer[6] + off;
  double *vyR     = buffer[7] + off;
  double *vzL     = buffer[8] + off;
  double *vzR     = buffer[9] + off;

  double rhomins, epslmins, vxmins, vymins, vzmins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus;

  double alpha, beta, detg_i, detginv, gupxx, gupyy, gupzz;
  double gxx, gxy, gxz, gyy, gyz, gzz;
  
  double pres_L, v2_L, Wlor_L, W2rhoh_L, cs2_L, vlowx_L, vlowy_L, vlowz_L;
  double pres_R, v2_R, Wlor_R, W2rhoh_R, cs2_R, vlowx_R, vlowy_R, vlowz_R;

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

  Rec1D_PPM_grhd(w1d[0], w1d[1], w1d[2], w1d[3], w1d[4], w1d[5], 
		              rhoL, epslL, vxL, vyL, vzL,
		              rhoR, epslR, vxR, vyR, vzR, 
                  nx);


  for( i = 0 ; i <= nx ; i++ ) {
    
    // if the point is atm, do nothing EXCEPT
    // set the fluxes to 0 to be sure no random numbers are added
    if ( (MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i]= ZERO;
      continue;
    }
    
    alpha = MATTER.rec1dml( g1d[0], i );
    beta  = MATTER.rec1dml( g1d[1], i );
    gxx   = MATTER.rec1dml( g1d[2], i );
    gxy   = MATTER.rec1dml( g1d[3], i );
    gxz   = MATTER.rec1dml( g1d[4], i );
    gyy   = MATTER.rec1dml( g1d[5], i );
    gyz   = MATTER.rec1dml( g1d[6], i );
    gzz   = MATTER.rec1dml( g1d[7], i );

    rhomins  = rhoL[i+1];
    epslmins = epslL[i+1];
    vxmins   = vxL[i+1];
    vymins   = vyL[i+1];
    vzmins   = vzL[i+1];

    rhoplus  = rhoR[i];
    epslplus = epslR[i];
    vxplus   = vxR[i];
    vyplus   = vyR[i];
    vzplus   = vzR[i];

    // eos
    EOS.comp("re","","","pc","","", rhoplus,epslplus, &pres_L,&cs2_L);
    EOS.comp("re","","","pc","","", rhomins,epslmins, &pres_R,&cs2_R);
    
    // this saves time (equivalent for rec with symmetric stencils wrt i+1/2)
    detg_i = detg(gxx,gxy,gxz,gyy,gyz,gzz);
    detginv   = 1./detg_i;
    gupxx = detginv*(gyy*gzz - gyz*gyz);
    gupyy = detginv*(gxx*gzz - gxz*gxz);
    gupzz = detginv*(gxx*gyy - gxy*gxy);

    // v2
    grhd_compute_v2_pt(vxplus,vyplus,vzplus, gxx,gxy,gxz,gyy,gyz,gzz,
		                  &vlowx_L,&vlowy_L,&vlowz_L, &v2_L );

    grhd_compute_v2_pt(vxmins,vymins,vzmins, gxx,gxy,gxz,gyy,gyz,gzz,
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
    grhd_compute_q_pt(gxx,gxy,gxz,gyy,gyz,gzz, detg_i, 
		                  rhoplus,epslplus,pres_L, vlowx_L,vlowy_L,vlowz_L, Wlor_L,
		                  &(qL[0]), &(qL[1]), &(qL[2]), &(qL[3]), &(qL[4]));
    
    grhd_compute_q_pt(gxx,gxy,gxz,gyy,gyz,gzz, detg_i,
                      rhomins,epslmins,pres_R, vlowx_R,vlowy_R,vlowz_R, Wlor_R,
                      &(qR[0]), &(qR[1]), &(qR[2]), &(qR[3]), &(qR[4]));

    // compute phys fluxes
    grhd_compute_fx_pt(alpha, beta, detg_i, 
		                   pres_L, vxplus, 
	            	       qL[0], qL[1], qL[2], qL[3], qL[4],
	            	       &(fL[0]), &(fL[1]), &(fL[2]), &(fL[3]), &(fL[4]));
    
    grhd_compute_fx_pt(alpha, beta, detg_i,
                       pres_R, vxmins,
                       qR[0], qR[1], qR[2], qR[3], qR[4],
                       &(fR[0]), &(fR[1]), &(fR[2]), &(fR[3]), &(fR[4]));


    // eigenvalues lam^0, lam^-, lam^+
    grhd_eig(alpha, beta, vxplus,v2_L,gupxx, cs2_L,
	     &(lamL[1]),&(lamL[0]),&(lamL[4]));
    lamL[2] = lamL[3] = lamL[1];

    grhd_eig(alpha, beta, vxmins,v2_R,gupxx, cs2_R,
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
      if (DEBUG && direction==1) {
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

    for (int l=0; l<10; l++) 
      free(buffer[l]);
    free(buffer);

}

// 2nd order finite volume fluxes based on HLLC Riemann solver in local Minkowski (K. Kiuchi 2022 https://arxiv.org/abs/2205.04487)
// and primitive reconstruction with PPM
void grhd_flx1d_fvpr_ppm_HLLC(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx, int ngx, int direction, double **flux1d)
{

  int Ntot = nx+2*ngx;
  int off  = ngx-1;

  /* contains all informations */
  double **buffer = (double**) malloc (10*sizeof(double*));
  for (int l=0; l<10; l++) 
    buffer[l] = (double*) malloc ((Ntot)*sizeof(double));

  double *rhoL    = buffer[0] + off;
  double *rhoR    = buffer[1] + off;
  double *epslL   = buffer[2] + off;
  double *epslR   = buffer[3] + off;
  double *vxL     = buffer[4] + off;
  double *vxR     = buffer[5] + off;
  double *vyL     = buffer[6] + off;
  double *vyR     = buffer[7] + off;
  double *vzL     = buffer[8] + off;
  double *vzR     = buffer[9] + off; 
  
  double dummy;

  double rhomins, epslmins, vxmins, vymins, vzmins, Gammamins;
  double rhoplus, epslplus, vxplus, vyplus, vzplus, Gammaplus;
  
  double gxx, gxy, gxz, gyy, gyz, gzz;
  double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
  double betax, betay, betaz, betax_low, betay_low, betaz_low;
  double alpha, detg_i, sqrtdetg;
  
  double pres_L, v2_L, Wlor_L, cs2_L, vxtet_L, vytet_L, vztet_L;
  double pres_R, v2_R, Wlor_R, cs2_R, vxtet_R, vytet_R, vztet_R;

  double mach_local, phi, contact_indicator;

  double B,C,D;

  // tetrad elements with upper index. first index i tetrad element second one is spacetime index
  double ett, etx, ety, etz;
  double exx, exy, exz;
  double eyy, eyz;
  double ezz;  

  // tetrad elements with lower index. first index i tetrad element second one is spacetime index
  double dtt;
  double dxt, dxx;
  double dyt, dyx, dyy;
  double dzt, dzx, dzy, dzz;

  // Rieman solver vars
  double lambda_L, lambda_R, lambda_c;
  double Pc, v_interface, Delta;
  double AL, AR, BL, BR;
  double a, b, c;

  double *qL     = dmalloc(MATTER.NVq);
  double *qR     = dmalloc(MATTER.NVq);  
  double *qc     = dmalloc(MATTER.NVq);
  double *fL     = dmalloc(MATTER.NVf);  
  double *fR     = dmalloc(MATTER.NVf);  
  double *f_i    = dmalloc(MATTER.NVf); 
  double *f_HLL  = dmalloc(MATTER.NVf);  
  double *q_HLL  = dmalloc(MATTER.NVf);  
  double *f_HLLC = dmalloc(MATTER.NVf);  
  double *q_HLLC = dmalloc(MATTER.NVf);
  double *lamL   = dmalloc(MATTER.NVf);
  double *lamR   = dmalloc(MATTER.NVf); 

  int i,k;

  Rec1D_PPM_grhd(w1d[0], w1d[1], w1d[2], w1d[3], w1d[4], w1d[5], 
		              rhoL, epslL, vxL, vyL, vzL,
		              rhoR, epslR, vxR, vyR, vzR, 
                  nx);

  //if (PR) printf(" %d %d %d \n",MATTER.NVw,MATTER.NVq,MATTER.NVf);

  for( i = 0 ; i <= nx ; i++ ) {
    
    // if the point is atm, do nothing EXCEPT
    // set the fluxes to 0 to be sure no random numbers are added
    if ( (MATTER.USEMASK) && ((mask1d[i]>0.9) || (mask1d[i+1]>0.9))) {
      for (k=0; k<MATTER.NVf; k++) flux1d[k][i]= ZERO;
      continue;
    }
    
    // rec metric     
    alpha = MATTER.rec1dml( g1d[0], i );
    betax = MATTER.rec1dml( g1d[1], i );
    betay = MATTER.rec1dml( g1d[2], i );
    betaz = MATTER.rec1dml( g1d[3], i );
    gxx   = MATTER.rec1dml( g1d[4], i );
    gxy   = MATTER.rec1dml( g1d[5], i );
    gxz   = MATTER.rec1dml( g1d[6], i );
    gyy   = MATTER.rec1dml( g1d[7], i );
    gyz   = MATTER.rec1dml( g1d[8], i );
    gzz   = MATTER.rec1dml( g1d[9], i );

    // left state, u_L =  u(i)^+, rec prim     
    rhomins  = rhoL[i+1];
    epslmins = epslL[i+1];
    vxmins   = vxL[i+1];
    vymins   = vyL[i+1];
    vzmins   = vzL[i+1];

    // left state, u_R =  u(i)^-, rec prim     
    rhoplus  = rhoR[i];
    epslplus = epslR[i];
    vxplus   = vxR[i];
    vyplus   = vyR[i];
    vzplus   = vzR[i];

    // eos
    EOS.comp("re","","","pc","","", rhoplus,epslplus, &pres_L,&cs2_L);
    EOS.comp("re","","","pc","","", rhomins,epslmins, &pres_R,&cs2_R);
    Gammaplus = EOS_Gamma(rhoplus, pres_L, epslplus, cs2_L);
    Gammamins = EOS_Gamma(rhomins, pres_R, epslmins, cs2_R);

    // detg and g^-1
    grhd_compute_detg_invg_pt(gxx,gxy,gxz,gyy,gyz,gzz,
                            &detg_i, 
                            &gupxx,&gupxy,&gupxz,&gupyy,&gupyz,&gupzz);
    sqrtdetg   = sqrt(detg_i);

    //beta low
    grhd_compute_v2_pt(betax,betay,betaz, gxx,gxy,gxz,gyy,gyz,gzz,
		                  &betax_low,&betay_low,&betaz_low, &dummy);

    COMPUTE_TETRAD_COMPONENTS
    
    // compute velocity in tetrad frame
    vxtet_L = dxx*vxplus;
    vytet_L = dyx*vxplus + dyy*vyplus;
    vztet_L = dzx*vxplus + dzy*vyplus + dzz*vzplus;

    vxtet_R = dxx*vxmins;
    vytet_R = dyx*vxmins + dyy*vymins;
    vztet_R = dzx*vxmins + dzy*vymins + dzz*vzmins;

    v2_L = pow(vxtet_L,2) + pow(vytet_L,2) + pow(vztet_L,2);
    v2_R = pow(vxtet_R,2) + pow(vytet_R,2) + pow(vztet_R,2);

    v2_L = DMIN(v2_L, GRHD.HRSC_VMAX); // needed ? testme !
    v2_R = DMIN(v2_R, GRHD.HRSC_VMAX);

    Wlor_L = ONE / sqrt(ONE - v2_L);
    Wlor_R = ONE / sqrt(ONE - v2_R);
    
    // compute conservative vars in tetrad frame
    // gamma_ii = alpha = detg = 1 and beta^i = 0 cause we are in the tetrad coordinrates now
    grhd_compute_q_pt(1,0,0,1,0,1,1, 
		      rhoplus,epslplus,pres_L,vxtet_L,vytet_L,vztet_L,Wlor_L,
		      &(qL[0]), &(qL[1]), &(qL[2]), &(qL[3]), &(qL[4]));

    grhd_compute_q_pt(1,0,0,1,0,1,1,
          rhomins,epslmins,pres_R,vxtet_R,vytet_R,vztet_R,Wlor_R,
          &(qR[0]), &(qR[1]), &(qR[2]), &(qR[3]), &(qR[4]));

    // compute fluxes in tetrad frame
    // also here alpha = detg = 1 and beta^i = 0 because of the tetrad frame
    grhd_compute_fx_pt(1,0,1, 
		      pres_L, vxtet_L, 
		      qL[0], qL[1], qL[2], qL[3], qL[4],
		      &(fL[0]), &(fL[1]), &(fL[2]), &(fL[3]), &(fL[4]));
    
    grhd_compute_fx_pt(1,0,1,
          pres_R, vxtet_R,
          qR[0], qR[1], qR[2], qR[3], qR[4],
          &(fR[0]), &(fR[1]), &(fR[2]), &(fR[3]), &(fR[4]));

    // eigenvalues lam^0, lam^-, lam^+ in tetrad frame
    // again alpha = detg = 1 and beta^i = 0 because of the tetrad frame
    grhd_eig(1, 0, vxtet_L, v2_L, 1, cs2_L,
	          &(lamL[1]),&(lamL[0]),&(lamL[4]));
    lamL[2] = lamL[3] = lamL[1];

    grhd_eig(1, 0, vxtet_R, v2_R, 1, cs2_R,
             &(lamR[1]),&(lamR[0]),&(lamR[4]));
    lamR[2] = lamR[3] = lamR[1];

    // now we have to go for the actual Riemann solver
    v_interface = betax/(alpha*sqrt(gupxx));

    HLLC_sided_solver(vxtet_L, vxtet_R,
                      pres_L, pres_R,
                      qL, qR, 
                      fL, fR, 
                      lamL, lamR, v_interface,
                      f_HLLC, q_HLLC, 
                      MATTER.NVf);

    HLL_solver(qL, qR, 
		           fL, fR, 
		           lamL, lamR, v_interface,
		           f_HLL, q_HLL, 
               MATTER.NVf);

    mach_local = DMIN(fabs(vxtet_L)/sqrt(cs2_L),fabs(vxtet_R)/sqrt(cs2_R));  
    phi = sin(DMIN(1,mach_local/mach_limit)*PI/2);
    contact_indicator = K*DMIN(Gammaplus,Gammamins)*fabs(w1d[0][i+1]-w1d[0][i])/DMIN(w1d[0][i+1],w1d[0][i]) 
                                                  - fabs(w1d[5][i+1]-w1d[5][i])/DMIN(w1d[5][i+1],w1d[5][i]);

    if (contact_indicator>=0) phi=1;

    for (k=0; k<MATTER.NVf; k++) {
      f_i[k] = f_HLLC[k] + (1-phi)*(f_HLL[k]-f_HLLC[k]);
      qc[k]  = q_HLLC[k] + (1-phi)*(q_HLL[k]-q_HLLC[k]);
    }
    // finally transform the fluxes back to the original frame and save them into flux1d
    flux1d[0][i] = sqrtdetg*alpha*(etx*qc[0] + exx*f_i[0]);
    flux1d[1][i] = sqrtdetg*alpha*(etx*qc[1] + exx*f_i[1]);

    flux1d[2][i] = sqrtdetg*alpha*(etx*(dxx*qc[2]+dyx*qc[3]+dzx*qc[4]) + exx*(dxx*f_i[2]+dyx*f_i[3]+dzx*f_i[4]));    
    flux1d[3][i] = sqrtdetg*alpha*(etx*(          dyy*qc[3]+dzy*qc[4]) + exx*(           dyy*f_i[3]+dzy*f_i[4]));    
    flux1d[4][i] = sqrtdetg*alpha*(etx*(                    dzz*qc[4]) + exx*(                      dzz*f_i[4]));  


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
      if (DEBUG && direction==1) {

        printf(" ************************** i =%d \n",i);
        printf("  var \t F_i \n");
        for (k=0; k<MATTER.NVf; k++) printf(" %d \t %e \n",k,flux1d[k][i]);
        
        printf(" LEFT\n");
        printf(" q \t f \t lam \n");
        for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qL[k],fL[k],lamL[k]);
        printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e %+.6e \n", rhoplus,epslplus,vxplus,vyplus,vzplus,v2_L,pres_L,lambda_L);
        
        printf(" RIGHT\n");
        printf(" q \t f \t lam \n");
        for (k=0; k<MATTER.NVf; k++) printf(" %+.6e %+.6e %+.6e\n",qR[k],fR[k],lamR[k]);
        printf(" w = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e  %+.6e %+.6e \n", rhomins,epslmins,vxmins,vymins,vzmins,v2_R,pres_R,lambda_R);

        printf(" CENTRAL\n");
        for (k=0; k<MATTER.NVf; k++) printf(" %+.6e \n",qc[k]);
        printf(" w = %+.6e %+.6e %+.6e %+.6e \n", lambda_c,Pc,v_interface,Delta);

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
  
  free(qL); free(qR); free(qc);  free(q_HLL); free(q_HLLC);
  free(fL); free(fR); free(f_i); free(f_HLL); free(f_HLLC);
  free(lamL); free(lamR);

  for (int l=0; l<10; l++) 
    free(buffer[l]);
  free(buffer); 

}

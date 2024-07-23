/* hrsc_HLLD.c 
   FS 06/20124 */

#include "bam.h"
#include "matter.h"

#define TINY 1e-30
#define ZERO 0.
#define HALF 0.5
#define ONE 1.0
#define DEBUG 0
#define PR 0
#define FLOOR 1e-10

enum {
    COLD,
    ATM,
    OK,
    FAIL,
    NR
};  


// compute conservatives in the aL/aR state given a guess for the pressure
void HLLD_a_state(double Ptot, double lambda, double E,
                  double Bx, double By, double Bz,
                  double RD, double RE, 
                  double RSx, double RSy, double RSz,
                  double RBx, double RBy, double RBz,
                  double *rhoha, double *Da, double *Ea, 
                  double *vxa, double *vya, double *vza,
                  double *Sxa, double *Sya, double *Sza,
                  double *Bxa, double *Bya, double *Bza) {

  double A = RSx - lambda*RE + Ptot*(1-lambda*lambda);
  double G = RBy*RBy + RBz*RBz;
  double C = RSy*RBy + RSz*RBz;
  double Q = -A -G + Bx*Bx*(1-lambda*lambda);
  double X = Bx*(A*lambda*Bx+C) - (A+G)*(lambda*Ptot+RE);

  double vx = ( Bx*(A*Bx+lambda*C)-(A+G)*(Ptot+RSx)) / X;
  double vy = ( Q*RSy + RBy*(C+Bx*(lambda*RSx-RE)) ) / X;
  double vz = ( Q*RSz + RBz*(C+Bx*(lambda*RSx-RE)) ) / X;
  
  *Bxa    = (RBx - Bx*vx)/(lambda-vx);
  *Bya    = (RBy - Bx*vy)/(lambda-vx);
  *Bza    = (RBz - Bx*vz)/(lambda-vx);

  double vB = vx*(*Bxa)+vy*(*Bya)+vz*(*Bza);

  *rhoha  = Ptot + (RE-vx*RSx-vy*RSy-vz*RSz)/(lambda-vx);
  *Da     = RD/(lambda-vx);
  *Ea     = (RE + Ptot*vx - vB*Bx)/(lambda-vx);
  *Sxa    = (*Ea+Ptot)*vx - vB*(*Bxa);
  *Sya    = (*Ea+Ptot)*vy - vB*(*Bya);
  *Sza    = (*Ea+Ptot)*vz - vB*(*Bza);
  *vxa = vx; 
  *vya = vy; 
  *vza = vz;

}

// compute hydro conservatives in the cL/cR state given a guess for the pressure
void HLLD_c_state(double Ptot, double lambda_a,
                  double Bx, double By, double Bz,
                  double Bxc, double Byc, double Bzc,
                  double vxa, double vya, double vza,
                  double vxc, double vyc, double vzc,
                  double Da, double Ea,
                  double Sxa, double Sya, double Sza,
                  double *Dc, double *Ec, 
                  double *Sxc, double *Syc, double *Szc) {

    double vBc = vxc*Bxc + vyc*Byc + vzc*Bzc;
    *Dc = Da * (lambda_a - vxa)/(lambda_a - vxc);
    *Ec = (lambda_a*Ea - Sxa + Ptot*vxc - vBc*Bx) / (lambda_a - vxc);
    *Sxc = (*Ec + Ptot)*vxc - vBc*Bxc;
    *Syc = (*Ec + Ptot)*vyc - vBc*Byc;
    *Szc = (*Ec + Ptot)*vzc - vBc*Bzc;

}

// compute K^i in the cL/cR state given a guess for the pressure
void HLLD_c_K_eta(double Ptot, double lambda, double rhoh, int dir,
                  double Bx,  double RE, 
                  double RSx, double RSy, double RSz,
                  double RBx, double RBy, double RBz,
                  double *eta, double *lambda_a,
                  double *Kx, double *Ky, double *Kz) {

  int sign;
  if (Bx>=0) sign=1; else sign=-1;
  *eta = dir*sign*sqrt(rhoh);
  *Kx = (RSx + Ptot + *eta*RBx)/(lambda*Ptot + RE + *eta*Bx);
  *Ky = (RSy        + *eta*RBy)/(lambda*Ptot + RE + *eta*Bx);
  *Kz = (RSz        + *eta*RBz)/(lambda*Ptot + RE + *eta*Bx);
  *lambda_a = *Kx;

}

// compute B^i in the cL/cR state given a guess for the pressure
void HLLD_c_B(double lambdaL,
              double vxL, double vyL, double vzL,
              double BxL, double ByL, double BzL,
              double lambdaR,
              double vxR, double vyR, double vzR,
              double BxR, double ByR, double BzR,
              double *Bxc, double *Byc, double *Bzc) {

    *Bxc = ( (BxR*(lambdaR-vxR) + BxR*vxR ) - (BxL*(lambdaL-vxL) + BxL*vxL ) ) / (lambdaR - lambdaL);
    *Byc = ( (ByR*(lambdaR-vxR) + BxR*vyR ) - (ByL*(lambdaL-vxL) + BxL*vyL ) ) / (lambdaR - lambdaL);
    *Bzc = ( (BzR*(lambdaR-vxR) + BxR*vzR ) - (BzL*(lambdaL-vxL) + BxL*vzL ) ) / (lambdaR - lambdaL);

}

double HLLD_f(double Ptot, // variable we have to solve for
              
              double lambdaL, double EL,
              double BxL, double ByL, double BzL,
              double RDL, double REL, 
              double RSxL, double RSyL, double RSzL,
              double RBxL, double RByL, double RBzL,
              
              double lambdaR, double ER,
              double BxR, double ByR, double BzR,
              double RDR, double RER, 
              double RSxR, double RSyR, double RSzR,
              double RBxR, double RByR, double RBzR) {

    double rhohaL, DaL, EaL, SxaL, SyaL, SzaL, BxaL, ByaL, BzaL;
    double rhohaR, DaR, EaR, SxaR, SyaR, SzaR, BxaR, ByaR, BzaR;

    double etaL, lambda_aL, KxL, KyL, KzL;
    double etaR, lambda_aR, KxR, KyR, KzR;

    double Bxc, Byc, Bzc;

    double vxaL, vyaL, vzaL;
    double vxaR, vyaR, vzaR;

    HLLD_a_state(Ptot, lambdaL, EL,
                BxL, ByL, BzL,
                RDL, REL, 
                RSxL, RSyL, RSzL,
                RBxL, RByL, RBzL,
                &rhohaL, &DaL, &EaL,
                &vxaL, &vyaL, &vzaL, 
                &SxaL, &SyaL, &SzaL,
                &BxaL, &ByaL, &BzaL);

    HLLD_a_state(Ptot, lambdaR, ER,
                BxR, ByR, BzR,
                RDR, RER, 
                RSxR, RSyR, RSzR,
                RBxR, RByR, RBzR,
                &rhohaR, &DaR, &EaR,
                &vxaR, &vyaR, &vzaR,  
                &SxaR, &SyaR, &SzaR,
                &BxaR, &ByaR, &BzaR);

    HLLD_c_K_eta(Ptot, lambdaL, rhohaL, -1,
                 BxL,  REL, 
                 RSxL, RSyL, RSzL,
                 RBxL, RByL, RBzL,
                 &etaL, &lambda_aL,
                 &KxL, &KyL, &KzL);

    HLLD_c_K_eta(Ptot, lambdaR, rhohaR, 1,
                 BxR,  RER, 
                 RSxR, RSyR, RSzR,
                 RBxR, RByR, RBzR,
                 &etaR, &lambda_aR,
                 &KxR, &KyR, &KzR);

    if (fabs(lambda_aR-lambda_aL)<1e-5*fabs(lambdaR-lambdaL)) return lambda_aR-lambda_aL;
    else {
      HLLD_c_B(lambda_aL,
                vxaL, vyaL, vzaL,
                BxaL, ByaL, BzaL,
                lambda_aR,
                vxaR, vyaR, vzaR,
                BxaR, ByaR, BzaR,
                &Bxc, &Byc, &Bzc);

      // printf("rhohaL, DaL, EaL, vxaL, vyaL, vzaL, SxaL, SyaL, SzaL = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e \n", rhohaL, DaL, EaL, vxaL, vyaL, vzaL, SxaL, SyaL, SzaL);
      // printf("rhohaR, DaR, EaR, vxaR, vyaR, vzaR, SxaR, SyaR, SzaR = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e \n", rhohaR, DaR, EaR, vxaR, vyaR, vzaR, SxaR, SyaR, SzaR);

      // printf("BxaL, ByaL, BzaL = %+.6e %+.6e %+.6e \n", BxaL, ByaL, BzaL);
      // printf("BxaR, ByaR, BzaR = %+.6e %+.6e %+.6e \n", BxaL, ByaL, BzaL);

      // printf("RDL, REL, RSxL, RSyL, RSzL, RBxL, RByL, RBzL = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e \n", RDL, REL, RSxL, RSyL, RSzL, RBxL, RByL, RBzL);
      // printf("RDR, RER, RSxR, RSyR, RSzR, RBxR, RByR, RBzR = %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e %+.6e \n", RDR, RER, RSxR, RSyR, RSzR, RBxR, RByR, RBzR);

      // printf("etaL, lambda_aL, KxL, KyL, KzL = %+.6e %+.6e %+.6e %+.6e %+.6e \n", etaL, lambda_aL, KxL, KyL, KzL);
      // printf("etaR, lambda_aR, KxR, KyR, KzL = %+.6e %+.6e %+.6e %+.6e %+.6e \n", etaR, lambda_aR, KxR, KyR, KzR);

      // printf("BxcR, BycR, BzcR = %+.6e %+.6e %+.6e \n", Bxc, Byc, Bzc);

      double DeltaK = KxR - KxL;
      double K2L = KxL*KxL + KyL*KyL + KzL*KzL;
      double K2R = KxR*KxR + KyR*KyR + KzR*KzR;
      double KBL = KxL*Bxc + KyL*Byc + KzL*Bzc;
      double KBR = KxR*Bxc + KyR*Byc + KzR*Bzc;

      double YR = (1-K2R) / (etaR - KBR + TINY );
      double YL = (1-K2L) / (etaL - KBL + TINY );

      return (DeltaK - Bxc*(YR-YL));
    }

}

// Dekker root finder for central total pressure in the Riemann fan
int HLLD_root_finder(double p0, double *res, 
              
              double lambdaL, double EL,
              double BxL, double ByL, double BzL,
              double RDL, double REL, 
              double RSxL, double RSyL, double RSzL,
              double RBxL, double RByL, double RBzL,
              
              double lambdaR, double ER,
              double BxR, double ByR, double BzR,
              double RDR, double RER, 
              double RSxR, double RSyR, double RSzR,
              double RBxR, double RByR, double RBzR,
              
              double delta_lambda) {

  double fmax, fb, fa, fbplus, fnew;
  double a, b, m, s, bnew;
  double error;
  double Delta = 1e-9;
  double errormax = 1e-7;
  double countmax = 1000;
  double pmax = 2*p0;
  double pmin = 0.1*p0;

  b = p0;

  fb = HLLD_f(b, 
              lambdaL,EL,  BxL,ByL,BzL,  RDL,REL,  RSxL,RSyL,RSzL,  RBxL,RByL,RBzL,
              lambdaR,ER,  BxR,ByR,BzR,  RDR,RER,  RSxR,RSyR,RSzR,  RBxR,RByR,RBzR);

  fmax = HLLD_f(pmax, 
                lambdaL,EL,  BxL,ByL,BzL,  RDL,REL,  RSxL,RSyL,RSzL,  RBxL,RByL,RBzL,
                lambdaR,ER,  BxR,ByR,BzR,  RDR,RER,  RSxR,RSyR,RSzR,  RBxR,RByR,RBzR);

	if (fmax <= 0) {
		if (fb > 0) {a=pmax;}
		else        {a=pmin;}
	}

	else {
		if (fb > 0) {a=pmin;}
		else        {a=pmax;}
	}

  fa = HLLD_f(a, 
              lambdaL,EL,  BxL,ByL,BzL,  RDL,REL,  RSxL,RSyL,RSzL,  RBxL,RByL,RBzL,
              lambdaR,ER,  BxR,ByR,BzR,  RDR,RER,  RSxR,RSyR,RSzR,  RBxR,RByR,RBzR);

	//start searching
  int count;
   for (count=0; count<=countmax; count++) {

    if (DEBUG) printf("a, b, fb, iter = %+.6e %+.6e %+.6e %i \n", a, b, fb, count);

		error = fabs(fb);
	  if (error<=errormax || error!=error) break;

		m = 0.5*(a+b);

    fbplus = HLLD_f(b*(1+Delta), 
                    lambdaL,EL,  BxL,ByL,BzL,  RDL,REL,  RSxL,RSyL,RSzL,  RBxL,RByL,RBzL,
                    lambdaR,ER,  BxR,ByR,BzR,  RDR,RER,  RSxR,RSyR,RSzR,  RBxR,RByR,RBzR);

    if (fbplus!=fb) s = b - (fbplus * b*Delta) / (fbplus - fb);
    else            s = m;

  	// evaluate b_(k+1)
		if      ((m<=s)&(s<=b)) bnew = s;
		else if ((b<=s)&(s<=m)) bnew = s;
		else                    bnew = m;

    fnew = HLLD_f(bnew, 
                  lambdaL,EL,  BxL,ByL,BzL,  RDL,REL,  RSxL,RSyL,RSzL,  RBxL,RByL,RBzL,
                  lambdaR,ER,  BxR,ByR,BzR,  RDR,RER,  RSxR,RSyR,RSzR,  RBxR, RByR, RBzR);

    if (fnew*fb < 0) {a = b; fa = fb;}
  
		b    = bnew;
		fb   = fnew;

  }

  *res = b;

  if (error<10*errormax) return 1;
  //else                return 0;
  else {
    if (PR) printf("HLLD pressure root finder failed, switched to HLL \n");
    if (PR) printf("a, b, fb, iter = %+.6e %+.6e %+.6e %i \n", a, b, fb, count);
    return 0;
  }

}

// this HLL implementation is ment to work in the tetrad frame
double HLL_solver_mag_tet(double *qL, double *qR, 
		                      double *fL, double *fR,
                          double *lamL, double *lamR,
                          double v_interface, double *fHLL, double *qHLL,
                          int nf) {
 
  double amins =  1.;
  double aplus = -1.;

  // max speeds (excluding the eigenvalues for divergence cleaning)
  for(int v=0; v<nf-2; v++) {
    aplus = DMAX( aplus, DMAX(lamR[v], lamL[v]) );
    amins = DMIN( amins, DMIN(lamR[v], lamL[v]) );
  }

  double oda = ONE/(aplus - amins + TINY);
  double apm = aplus * amins;

  if      (amins>=v_interface) for (int v=0; v<nf-1; v++) {fHLL[v]=fL[v]; qHLL[v]=qL[v];}
  else if (aplus<=v_interface) for (int v=0; v<nf-1; v++) {fHLL[v]=fR[v]; qHLL[v]=qR[v];}
  else for (int v=0; v<nf-1; v++) {
    fHLL[v] = oda * ( aplus * fL[v] - amins * fR[v] + apm * (qR[v] - qL[v]) );
    qHLL[v] = oda * ( aplus * qR[v] - amins * qL[v] +       (fL[v] - fR[v]) );
  }

}

double HLLD_initial_guess(double *qL, double *qR, 
		                      double *fL, double *fR,
                          double *lamL, double *lamR, 
                          double *qHLL, double *fHLL, 
                          double Ptot, double mu_in,
                          int nf) {

  if (DEBUG) printf("fD, fTau, fSx, fSy, fSz, fBx, fBy, fBz = %+.6e, %+.6e, %+.6e, %+.6e, %+.6e, %+.6e, %+.6e, %+.6e \n", fHLL[0],fHLL[1],fHLL[2],fHLL[3],fHLL[4],fHLL[5],fHLL[6],fHLL[7]);
  if (DEBUG) printf(" D,  Tau,  Sx,  Sy,  Sz,  Bx,  By,  Bz = %+.6e, %+.6e, %+.6e, %+.6e, %+.6e, %+.6e, %+.6e, %+.6e \n", qHLL[0],qHLL[1],qHLL[2],qHLL[3],qHLL[4],qHLL[5],qHLL[6],qHLL[7]);

  // Now that we have the HLL central states and fluxes,
  // we proceed computing the central pressure, 
  // there are 4 different scenarios based on Bx and vx

  if (pow(qR[5]+qL[5],2)*0.25/Ptot>0.1) {
    // int result;
    // double rho, epsl, press;
    // double mu = mu_in;
    // int converged;
    // /* define rescale variables as in Kastaun 2021 */
    // double bx   = qHLL[5]/sqrt(qHLL[0]);
    // double by   = qHLL[6]/sqrt(qHLL[0]);
    // double bz   = qHLL[7]/sqrt(qHLL[0]);
    // double rx   = qHLL[2]/qHLL[0];
    // double ry   = qHLL[3]/qHLL[0];
    // double rz   = qHLL[4]/qHLL[0];
    // double r2   = (rx*rx+ry*ry+rz*rz);
    // double rb   = (rx*bx+ry*by+rz*bz);
    // double rb2  = rb*rb;
    // double b2   = (bx*bx+by*by+bz*bz);
    // double bro2 = b2*r2 - rb2;
    // double q    = (qHLL[1]-qHLL[0])/qHLL[0];
    // double d    = qHLL[0];

    // result = OK;

    // // check grhd variables for NANs
    // if (CheckForNANandINF(8,qHLL[0],qHLL[1],qHLL[2],qHLL[3],qHLL[4],qHLL[5],qHLL[6],qHLL[7])) {
    //   if (0) printf("conD, conS[xyz], conT, conBup[xyz], conB2, or SQRTconB2 is nan\n");
    //   result = ATM;
    // }

    // if (d<=0)
    //   result = ATM;

    // if ((result==OK)) converged = reprimand(&mu, b2, r2, rb2, bro2, q, d,
    //                                           &rho,&epsl,&press, 1e-6,100);
    // else converged = 0;

    // double x = 1 / (1 + mu * b2);
    // double vx =  mu * x * (rx + (rb*mu) *bx);
    // double vy =  mu * x * (ry + (rb*mu) *by);
    // double vz =  mu * x * (rz + (rb*mu) *bz);

    // double v2 = vx*vx + vy*vy + vz*vz;
    // double W = 1/sqrt(1-v2);

    // double B2 = qHLL[5]*qHLL[5] + qHLL[6]*qHLL[6] + qHLL[7]*qHLL[7];
    // double Bv = vx*qHLL[5] + vy*qHLL[6] + vz*qHLL[7];

    // double mag_press = 0.5*(B2/pow(W,2) + Bv*Bv);

    // //printf("Iinitial guess = %+.6e \n",press + mag_press);
    // return press + mag_press;

    return Ptot;

  } else {
    double b = qHLL[1]-fHLL[2];
    double c = qHLL[2]*fHLL[1] - fHLL[2]*qHLL[1];
    double Delta = b*b - 4*c;
    return (-b + sqrt(Delta))*0.5;
  }
}

void HLLD_solver(double *qL, double *qR, 
		             double *fL, double *fR, 
		             double *lamL, double *lamR, 
                 double v_interface, double *qHLL, double *fHLL, 
                 double Ptot, double mu, double hL, double hR,
		             double *f, double *q, int nf) {

  int consistency = 1;
  double lambdaL =  1.;
  double lambdaR = -1.;
  double vxcL, vxcR;
  double lambdac;

  // max speeds (excluding the eigenvalues for divergence cleaning)
  for(int v=0; v<nf-2; v++) {
    lambdaR = DMAX( lambdaR, DMAX(lamR[v], lamL[v]) );
    lambdaL = DMIN( lambdaL, DMIN(lamR[v], lamL[v]) );
  }

  if      (v_interface<=lambdaL) for (int v=0; v<nf-1; v++) {f[v]=fL[v];q[v]=qL[v];}
  else if (v_interface>=lambdaR) for (int v=0; v<nf-1; v++) {f[v]=fR[v];q[v]=qR[v];}

  else {

    double *qaL = dmalloc(nf);
    double *qaR = dmalloc(nf);
    double rhohaL, rhohaR;
    double vxaL, vyaL, vzaL;
    double vxaR, vyaR, vzaR;
    double lambda_aL, lambda_aR;
    double KxL, KyL, KzL, etaL;
    double KxR, KyR, KzR, etaR;
    double v_avg = (fL[0]/qL[0] + fR[0]/qR[0])*0.5;

    double RDL  = lambdaL*qL[0] - fL[0];   double RDR  = lambdaR*qR[0] - fR[0];
    double REL  = lambdaL*qL[1] - fL[1];   double RER  = lambdaR*qR[1] - fR[1];
    double RSxL = lambdaL*qL[2] - fL[2];   double RSxR = lambdaR*qR[2] - fR[2];
    double RSyL = lambdaL*qL[3] - fL[3];   double RSyR = lambdaR*qR[3] - fR[3];
    double RSzL = lambdaL*qL[4] - fL[4];   double RSzR = lambdaR*qR[4] - fR[4];
    double RBxL = lambdaL*qL[5] - fL[5];   double RBxR = lambdaR*qR[5] - fR[5];
    double RByL = lambdaL*qL[6] - fL[6];   double RByR = lambdaR*qR[6] - fR[6];
    double RBzL = lambdaL*qL[7] - fL[7];   double RBzR = lambdaR*qR[7] - fR[7];

    double p0 = HLLD_initial_guess(qL,qR,fL,fR,lamL,lamR,qHLL,fHLL,Ptot,mu,nf);
    double pc;
    int convergence = HLLD_root_finder(p0, &pc, 
                
                    lambdaL, qL[1],
                    qL[5], qL[6], qL[7],
                    RDL, REL, 
                    RSxL, RSyL, RSzL,
                    RBxL, RByL, RBzL,
                
                    lambdaR, qR[1],
                    qR[5], qR[6], qR[7],
                    RDR, RER, 
                    RSxR, RSyR, RSzR,
                    RBxR, RByR, RBzR,
                    
                    1);

    if (convergence & (hL>=pc) & (hR>=pc)) {

      HLLD_a_state(pc, lambdaL, qL[1],
                  qL[5], qL[6], qL[7],
                  RDL, REL, 
                  RSxL, RSyL, RSzL,
                  RBxL, RByL, RBzL,
                  &rhohaL, &(qaL[0]), &(qaL[1]), 
                  &vxaL, &vyaL, &vzaL,
                  &(qaL[2]), &(qaL[3]), &(qaL[4]),
                  &(qaL[5]), &(qaL[6]), &(qaL[7]));

      HLLD_a_state(pc, lambdaR, qR[1],
                  qR[5], qR[6], qR[7],
                  RDR, RER, 
                  RSxR, RSyR, RSzR,
                  RBxR, RByR, RBzR,
                  &rhohaR, &(qaR[0]), &(qaR[1]), 
                  &vxaR, &vyaR, &vzaR,
                  &(qaR[2]), &(qaR[3]), &(qaR[4]),
                  &(qaR[5]), &(qaR[6]), &(qaR[7]));

      HLLD_c_K_eta(pc, lambdaL, rhohaL, -1,
                  qL[5], REL, 
                  RSxL, RSyL, RSzL,
                  RBxL, RByL, RBzL,
                  &etaL, &lambda_aL,
                  &KxL, &KyL, &KzL);

      HLLD_c_K_eta(pc, lambdaR, rhohaR, 1,
                  qR[5], RER, 
                  RSxR, RSyR, RSzR,
                  RBxR, RByR, RBzR,
                  &etaR, &lambda_aR,
                  &KxR, &KyR, &KzR);

      //if (fabs(lambda_aL-lambdac)<1e-6) {lambda_aL = lambdac; vxcL = lambdac; vxaL = lambdac;}
      //if (fabs(lambda_aR-lambdac)<1e-6) {lambda_aR = lambdac; vxcR = lambdac; vxaR = lambdac;}

      if (fabs(lambda_aR-lambda_aL)<1e-5*fabs(lambdaR-lambdaL)) {lambdac=0.5*(lambda_aR+lambda_aL); lambda_aR=lambdac; lambda_aL=lambdac;}

      if      (v_interface<lambda_aL) for (int v=0; v<nf-1; v++) {f[v]=fL[v]+lambdaL*(qaL[v]-qL[v]);q[v]=qaL[v];}
      else if (v_interface>lambda_aR) for (int v=0; v<nf-1; v++) {f[v]=fR[v]+lambdaR*(qaR[v]-qR[v]);q[v]=qaR[v];}
      else if ((v_interface==lambda_aL) && 
              (v_interface==lambda_aR)) for (int v=0; v<nf-1; v++) {
                                              f[v] = (fR[v] + lambdaR*(qaR[v]-qR[v]) +
                                                      fL[v] + lambdaL*(qaL[v]-qL[v]) ) * 0.5;
                                              q[v] = 0.5*(qaR[v]+qaL[v]);
                                            }

      else {

        double *qc = dmalloc(nf);
        double vyc, vzc;
        double K2L, K2R, KBL, KBR;

        HLLD_c_B(lambdaL,
                  vxaL, vyaL, vzaL,
                  qaL[5], qaL[6], qaL[7],
                  lambdaR,
                  vxaR, vyaR, vzaR,
                  qaR[5], qaR[6], qaR[7],
                  &(qc[5]), &(qc[6]), &(qc[7]));
          
        K2L = KxL*KxL + KyL*KyL + KzL*KzL;
        K2R = KxR*KxR + KyR*KyR + KzR*KzR;

        KBL = KxL*qc[5] + KyL*qc[6] + KzL*qc[7];
        KBR = KxR*qc[5] + KyR*qc[6] + KzR*qc[7];

        vxcL = KxL - qc[5]*(1-K2L)/(etaL-KBL);
        vxcR = KxR - qc[5]*(1-K2R)/(etaR-KBR);

        lambdac = 0.5*(vxcL+vxcR);

        //if (PR) printf("lambda_L, lamda_aL, lambda_c, lambda_aR, lambda_R = %+.6e %+.6e %+.6e %+.6e %+.6e = \n", lambdaL, lambda_aL, lambdac, lambda_aR, lambdaR);

        //if (PR) printf("DL,  EL,  SxL,  SyL,  SzL   = %+.6e %+.6e %+.6e %+.6e %+.6e \n", qL[0], qL[1], qL[2], qL[3], qL[4]);
        //if (PR) printf("DaL, EaL, SxaL, SyaL, SzaL  = %+.6e %+.6e %+.6e %+.6e %+.6e \n", qaL[0], qaL[1], qaL[2], qaL[3], qaL[4]);

        //if (PR) printf("DR,  ER,  SxR,  SyR,  SzR   = %+.6e %+.6e %+.6e %+.6e %+.6e \n", qR[0], qR[1], qR[2], qR[3], qR[4]);
        //if (PR) printf("DaR, EaR, SxaR, SyaR, SzaR  = %+.6e %+.6e %+.6e %+.6e %+.6e \n", qaR[0], qaR[1], qaR[2], qaR[3], qaR[4]);

          if (v_interface<lambdac) {
              vyc = KyL - qc[6]*(1-K2L)/(etaL-KBL);
              vzc = KzL - qc[7]*(1-K2L)/(etaL-KBL);
              HLLD_c_state(pc, lambda_aL,
                          qL[5], qL[6], qL[7],
                          qc[5], qc[6], qc[7],
                          vxaL, vyaL, vzaL,
                          vxcL, vyc, vzc,
                          qaL[0], qaL[1],
                          qaL[2], qaL[3], qaL[4],
                          &(qc[0]), &(qc[1]), 
                          &(qc[2]), &(qc[3]), &(qc[4]));
              for (int v=0; v<nf-1; v++) {
                f[v] = fL[v] + lambdaL*(qaL[v]-qL[v]) + lambda_aL*(qc[v]-qaL[v]);
                q[v] = qc[v];
              }
            }

          else if (v_interface>lambdac) {
              vyc = KyR - qc[6]*(1-K2R)/(etaR-KBR);
              vzc = KzR - qc[7]*(1-K2R)/(etaR-KBR);
              HLLD_c_state(pc, lambda_aR,
                          qR[5], qR[6], qR[7],
                          qc[5], qc[6], qc[7],
                          vxaR, vyaR, vzaR,
                          vxcR, vyc, vzc,
                          qaR[0], qaR[1],
                          qaR[2], qaR[3], qaR[4],
                          &(qc[0]), &(qc[1]), 
                          &(qc[2]), &(qc[3]), &(qc[4]));
              for (int v=0; v<nf-1; v++) {
                f[v] = fR[v] + lambdaR*(qaR[v]-qR[v]) + lambda_aR*(qc[v]-qaR[v]);
                q[v] = qc[v];
              }
            }

          else{ 
              // if v_interface=lambda_c, compute both central states and take the average
              // This can be optimized! But probably not needed
              // I mean we will never end up here, except in some test
              double *qc2 = dmalloc(nf);
              qc2[5]=qc[5]; qc2[6]=qc[6]; qc2[7]=qc[7];
              vyc = KyR - qc[6]*(1-K2R)/(etaR-KBR);
              vzc = KzR - qc[7]*(1-K2R)/(etaR-KBR);
              HLLD_c_state(pc, lambda_aR,
                          qR[5], qR[6], qR[7],
                          qc[5], qc[6], qc[7],
                          vxaR, vyaR, vzaR,
                          vxcR, vyc, vzc,
                          qaR[0], qaR[1],
                          qaR[2], qaR[3], qaR[4],
                          &(qc[0]), &(qc[1]), 
                          &(qc[2]), &(qc[3]), &(qc[4]));
              vyc = KyL - qc[6]*(1-K2L)/(etaR-KBL);
              vzc = KzL - qc[7]*(1-K2L)/(etaR-KBL);
              HLLD_c_state(pc, lambda_aL,
                          qL[5], qL[6], qL[7],
                          qc[5], qc[6], qc[7],
                          vxaL, vyaL, vzaL,
                          vxcL, vyc, vzc,
                          qaL[0], qaL[1],
                          qaL[2], qaL[3], qaL[4],
                          &(qc2[0]), &(qc2[1]), 
                          &(qc2[2]), &(qc2[3]), &(qc2[4]));
              for (int v=0; v<nf-4; v++) {
                f[v] = (fR[v] + lambdaR*(qaR[v]-qR[v]) + lambda_aR*(qc[v] -qaR[v]) +
                        fL[v] + lambdaL*(qaL[v]-qL[v]) + lambda_aL*(qc2[v]-qaL[v])) * 0.5;
                q[v] = 0.5*(qc[v]+qc2[v]);
              }
              free(qc2);
            }
            free(qc);
                  if ((vxaL>=lambdaL) & (vxaR<=lambdaR) & (vxcL>=lambda_aL) & (vxcR<=lambda_aR)) consistency=1;
            else {
              consistency=0;
              if (PR) printf("hL, Pc = %+.6e %+.6e \n",hL,pc);
              if (PR) printf("vax_L, lambda_L, vcx_L, lambda_aL = %+.6e %+.6e %+.6e %+.6e \n",vxaL,lambdaL,vxcL,lambda_aL);
              if (PR) printf("vax_R, lambda_R, vcx_R, lambda_aR = %+.6e %+.6e %+.6e %+.6e \n",vxaR,lambdaR,vxcR,lambda_aR);
            }
      }

      free(qaL); free(qaR);


    } else consistency=0;

    // if consistency check failed, or pressure root finder failed, switch to HLL
    if (consistency==0) {
      if (convergence==0) if (PR) printf("HLLD pressure root finder failed, switched to HLL \n");
      else                if (PR) printf("HLLD found a non-consistent state, switched to HLL \n");
      for (int v=0; v<nf-1; v++) {
        f[v] = fHLL[v];
        q[v] = qHLL[v];
      }
    }

  }

}
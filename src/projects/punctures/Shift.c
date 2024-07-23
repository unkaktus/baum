/* Shift.c */
/* Bernd Bruegmann 2/03 */
/* set shift appropriate for thin sandwich puncture data 
*/

#include "bam.h"
#include "punctures.h"



/* shift for thin sandwich (gtilde_dot = 0)
     (VectorLaplace v)^i = d^i,  v = 0 at boundary
   where
     v^i = beta^i
     d^i = 2 A^ij del_j (alpha psi^-6)
*/


/* set coefficient */
void SetPunctureShiftCoefficient(tL *level)
{
  double *alpha = Ptr(level, "alpha");
  double *psi = Ptr(level, "adm_psi");
  double *K11 = Ptr(level, "adm_Kxx");
  double *K12 = Ptr(level, "adm_Kxy");
  double *K13 = Ptr(level, "adm_Kxz");
  double *K22 = Ptr(level, "adm_Kyy");
  double *K23 = Ptr(level, "adm_Kyz");
  double *K33 = Ptr(level, "adm_Kzz");
  double *u   = Ptr(level, "punctures_u");
  double *v   = Ptr(level, "punctures_v");
  double *e1  = Ptr(level, "punctures_ex");
  double *e2  = Ptr(level, "punctures_ey");
  double *e3  = Ptr(level, "punctures_ez");
  double *s1  = Ptr(level, "punctures_sx");
  double *s2  = Ptr(level, "punctures_sy");
  double *s3  = Ptr(level, "punctures_sz");
  double *d   = PtrEnable(level, "punctures_d");
  double cx = 1/(2*level->dx);
  double cy = 1/(2*level->dy);
  double cz = 1/(2*level->dz);
  double dd1, dd2, dd3, f;
  int i;

  double psiinv, psiinv7;
  double du1, du2, du3;
  double dv1, dv2, dv3;
  double *dpop1 = Ptr(level, "adm_dpsiopsix");
  double *dpop2 = Ptr(level, "adm_dpsiopsiy");
  double *dpop3 = Ptr(level, "adm_dpsiopsiz");

  /* store function that is to be differentiated in temporary variable d */
  forallpoints(level, i) 
    d[i] = alpha[i] * pow(psi[i] + u[i], -6.0);

  /* compute coefficient e^i */
  if (Getv("punctures_shift_rhs", "Adel")) {
    forinner7(level) {
      dd1 = cx * (d[pcc] - d[mcc]);
      dd2 = cy * (d[cpc] - d[cmc]);
      dd3 = cz * (d[ccp] - d[ccm]);

      e1[ccc] = 2 * (K11[ccc]*dd1 + K12[ccc]*dd2 + K13[ccc]*dd3);
      e2[ccc] = 2 * (K12[ccc]*dd1 + K22[ccc]*dd2 + K23[ccc]*dd3);
      e3[ccc] = 2 * (K13[ccc]*dd1 + K23[ccc]*dd2 + K33[ccc]*dd3);
      
      if (0) {
	e1[ccc] = dd1;
	e2[ccc] = dd2;
	e3[ccc] = dd3;
      }
    } endfor;
  }

  /* alternative computation of coefficient e^i:
       pull A inside derivative
       this does not rely on a finite differenced, regular term cancelling
       the divergence in A
     can give factor of 2 difference right next to puncture, but at high
       resolutions the shift doesn't show much of that
  */
  else {
    forallpoints(level, i) {
      s1[i] = 2 * cx * K11[i] * d[i];
      s2[i] = 2 * cy * K12[i] * d[i];
      s3[i] = 2 * cz * K13[i] * d[i];
    }
    forinner7(level) {
      e1[ccc] = s1[pcc]-s1[mcc] + s2[cpc]-s2[cmc] + s3[ccp]-s3[ccm];
    } endfor;
    
    forallpoints(level, i) {
      s1[i] = 2 * cx * K12[i] * d[i];
      s2[i] = 2 * cy * K22[i] * d[i];
      s3[i] = 2 * cz * K23[i] * d[i];
    }
    forinner7(level) {
      e2[ccc] = s1[pcc]-s1[mcc] + s2[cpc]-s2[cmc] + s3[ccp]-s3[ccm];
    } endfor;
    
    forallpoints(level, i) {
      s1[i] = 2 * cx * K13[i] * d[i];
      s2[i] = 2 * cy * K23[i] * d[i];
      s3[i] = 2 * cz * K33[i] * d[i];
    }
    forinner7(level) {
      e3[ccc] = s1[pcc]-s1[mcc] + s2[cpc]-s2[cmc] + s3[ccp]-s3[ccm];
    } endfor;
    
    forallpoints(level, i)
      s1[i] = s2[i] = s3[i] = 0;
  }

  /* synchronize */
  bampi_synchronize(level, Ind("punctures_ex"));

  /* cleanup */
  if (Getv("punctures_persist", "no")) {
    disablevar(level, Ind("punctures_d"));
  }
}




/* set shift */
void SetPunctureShift(tL *level)
{
  double xgrid = level->bbox[1] - level->dx/2.0;
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  double *betax = Ptr(level, "betax");
  double *betay = Ptr(level, "betay");
  double *betaz = Ptr(level, "betaz");
  double *betax_BC0 = Ptr(level, "punctures_ex");
  double *betay_BC0 = Ptr(level, "punctures_ey");
  double *betaz_BC0 = Ptr(level, "punctures_ez");
  double betax0, betay0, betaz0;
  double betax1, betay1, betaz1;
  double betaxyz1_max;
  double Jx, Jy, Jz , I_orb ,   Omega_bx, Omega_by, Omega_bz , Omega_b;
  int i,n;
  FILE *out;
  char filename[10000];

  /* copy the beta computed with BC beta|_\infty=0 into punctures_e=beta_BC0 */
  forallpoints(level, i) 
  {
      betax_BC0[i] = betax[i];
      betay_BC0[i] = betay[i];
      betaz_BC0[i] = betaz[i];
  }

  /* get beta^i at punc. 0 and 1. beta is computed with BC beta|_\infty=0 */
  printf("SetPunctureShift: beta at the punctures, computed with BC beta|_infinity=0:\n");
  betax0=interpolate_xyz_scalar(level, CBL[0][0], CBL[0][1], CBL[0][2],
                                Ind("betax"), Geti("order_RP"),LAGRANGE);
  betax1=interpolate_xyz_scalar(level, CBL[1][0], CBL[1][1], CBL[1][2],
                                Ind("betax"), Geti("order_RP"),LAGRANGE);

  betay0=interpolate_xyz_scalar(level, CBL[0][0], CBL[0][1], CBL[0][2],
                                Ind("betay"), Geti("order_RP"),LAGRANGE);
  betay1=interpolate_xyz_scalar(level, CBL[1][0], CBL[1][1], CBL[1][2],
                                Ind("betay"), Geti("order_RP"),LAGRANGE);

  betaz0=interpolate_xyz_scalar(level, CBL[0][0], CBL[0][1], CBL[0][2],
                                Ind("betaz"), Geti("order_RP"),LAGRANGE);
  betaz1=interpolate_xyz_scalar(level, CBL[1][0], CBL[1][1], CBL[1][2],
                                Ind("betaz"), Geti("order_RP"),LAGRANGE);

  /* angular mom. J */
  Jx=0;  Jy=0;  Jz=0;
  for(n = 0; MBL[n] != 0.0 && n < N_MWBL; n++)
  {
    Jx += CBL[n][1]*PBL[n][2] - CBL[n][2]*PBL[n][1];
    Jy += CBL[n][2]*PBL[n][0] - CBL[n][0]*PBL[n][2];
    Jz += CBL[n][0]*PBL[n][1] - CBL[n][1]*PBL[n][0];
  }

  /* compute Omega_b = J/I_orb,
     such that beta_BC0 + Omega_b x R = 0 at each puncture.
     Here we have to assume that beta_BC0 is of the form -Omega_b x R  */
  /* find largest beta component to compute I_orb */
  I_orb = 0.0;
  betaxyz1_max = 0.0;
  if(betax1 != 0) 
  {
    betaxyz1_max = betax1;
    I_orb = -(Jy * CBL[1][2] - Jz * CBL[1][1])/betax1;
  }
  if( fabs(betay1) > fabs(betaxyz1_max) )
  {
    betaxyz1_max = betay1;
    I_orb = -(Jz * CBL[1][0] - Jx * CBL[1][2])/betay1;
  }
  if( fabs(betaz1) > fabs(betaxyz1_max) )
  {
    betaxyz1_max = betaz1;
    I_orb = -(Jx * CBL[1][1] - Jy * CBL[1][0])/betaz1;
  }
  
  /* now compute Omega_b = J/I_orb */
  if( I_orb != 0.0 )
  {
    Omega_bx = Jx/I_orb;
    Omega_by = Jy/I_orb;  
    Omega_bz = Jz/I_orb;
  }
  else
  {
    Omega_bx = 0.0;
    Omega_by = 0.0;
    Omega_bz = 0.0;
  }
  Omega_b = sqrt(Omega_bx*Omega_bx + Omega_by*Omega_by + Omega_bz*Omega_bz );
  printf("SetPunctureShift: Omega_b = %f\n",Omega_b);
  printf("SetPunctureShift: ");
  printf("Omega_bx = %f  Omega_by = %f  Omega_bz = %f \n"
         ,Omega_bx,Omega_by,Omega_bz);

  /*add beta_rot = Omega x R to the beta computed with BC beta|_\infty=0 */
  if (Getv("punctures_shift_add_beta_rot", "yes"))
  {
    printf("SetPunctureShift: Correcting shift beta by Omega_b cross R:\n");
    forallpoints(level, i) 
    {
       betax[i] += (Omega_by * z[i] - Omega_bz * y[i]);
       betay[i] += (Omega_bz * x[i] - Omega_bx * z[i]);
       betaz[i] += (Omega_bx * y[i] - Omega_by * x[i]);
    }
  }
  
  if (processor0)
  {
    /* write Omega_b in files */
    sprintf(filename, "%s/%s.%d", Gets("outdir"), 
                      "punctures_Omega_b", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, Omega_b);
    fclose(out);
  }

  /* old tests */
  if (0) {
    double Wx, Wy, Wz; 
    double *betax = Ptr(level, "betax");
    double *betay = Ptr(level, "betay");
    double *betaz = Ptr(level, "betaz");
    forallpoints(level, i) {
      BY_Wofxyz(x[i], y[i], z[i], &Wx, &Wy, &Wz);
      betax[i] += Wx;
      betay[i] += Wy;
      betaz[i] += Wz;
    }
  }
  if (0) {
    double *Wx = Ptr(level, "punctures_ex");
    double *Wy = Ptr(level, "punctures_ey");
    double *Wz = Ptr(level, "punctures_ez");
    forallpoints(level, i) {
      BY_Wofxyz(x[i], y[i], z[i], &Wx[i], &Wy[i], &Wz[i]);
      Wx[i] *= 2;
      Wy[i] *= 2;
      Wz[i] *= 2;
    }
  }
}





/* test BY shift */
void TestBYShift(tL *level)
{
  double *betax = Ptr(level, "betax");
  double *betay = Ptr(level, "betay");
  double *betaz = Ptr(level, "betaz");
  double *alpha = Ptr(level, "alpha");
  double *psi = Ptr(level, "adm_psi");
  double *u   = Ptr(level, "punctures_u");
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  double *K11 = Ptr(level, "adm_Kxx");
  double *K12 = Ptr(level, "adm_Kxy");
  double *K13 = Ptr(level, "adm_Kxz");
  double *K22 = Ptr(level, "adm_Kyy");
  double *K23 = Ptr(level, "adm_Kyz");
  double *K33 = Ptr(level, "adm_Kzz");
  double *r11 = Ptr(level, "punctures_lbrhsxx");
  double *r12 = Ptr(level, "punctures_lbrhsxy");
  double *r13 = Ptr(level, "punctures_lbrhsxz");
  double *r22 = Ptr(level, "punctures_lbrhsyy");
  double *r23 = Ptr(level, "punctures_lbrhsyz");
  double *r33 = Ptr(level, "punctures_lbrhszz");
  double p, r, n1, n2, n3, rp;
  int i;

  forallpoints(level, i) {
    if (1) {
      BY_Wofxyz(x[i], y[i], z[i], betax+i, betay+i, betaz+i);
      BY_Kofxyz(x[i], y[i], z[i],
		&K11[i], &K12[i], &K13[i], &K22[i], &K23[i], &K33[i]);
      r = 1;
      r11[i] = r*K11[i];
      r12[i] = r*K12[i];
      r13[i] = r*K13[i];
      r22[i] = r*K22[i];
      r23[i] = r*K23[i];
      r33[i] = r*K33[i];
    }
    else {
      p = 2;
      r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      if (dequal(r, 0.0)) r = EPSZERO;
      rp = pow(r, p);
      n1 = x[i]/r;
      n2 = y[i]/r;
      n3 = z[i]/r;

      betax[i] = rp * n1;
      betay[i] = rp * n2;
      betaz[i] = rp * n3;
      
      p = 2*(p-1)*rp/r;
      r11[i] = p*n1*n1 - p/3;
      r12[i] = p*n1*n2;
      r13[i] = p*n1*n3;
      r22[i] = p*n2*n2 - p/3;
      r23[i] = p*n2*n3;
      r33[i] = p*n3*n3 - p/3;
    }
    //alpha[i] = 0.25;
    //psi[i] = 1;
    //u[i] = 0;
  }
}




/* test shift 
   compute lbeta - 2 alpha psi^-6 A
*/
void TestPunctureShift(tL *level, tVarList *beta, tVarList *vlf, tVarList *vlr)
{
  if (Getv("punctures_shift_test", "yes")) {
    tVarList *dgdt  = VLPtrEnable1(level, "punctures_dgdtxx");
    tVarList *lbeta = VLPtrEnable1(level, "punctures_lbetaxx");
    tVarList *lbrhs = VLPtrEnable1(level, "punctures_lbrhsxx");
    double *alpha = Ptr(level, "alpha");
    double *psi = Ptr(level, "adm_psi");
    double *u   = Ptr(level, "punctures_u");
    double *ap, *lbetap, *lbrhsp, *dgdtp;
    int i, n;
    int testshift = 0;

    if (testshift) TestBYShift(level);

    VectorLieFlat(level, lbeta, beta);
    set_boundary_symmetry(level, lbeta);
    bampi_vlsynchronize(lbeta);

    for (n = 0; n < 6; n++) {
      ap = level->v[Ind("adm_Kxx") + n];
      lbetap = level->v[lbeta->index[n]];
      lbrhsp = level->v[lbrhs->index[n]];
      dgdtp  = level->v[dgdt->index[n]];
      forallpoints(level, i) {
	if (!testshift)
	  lbrhsp[i] = 2 * alpha[i] * pow(psi[i] + u[i], -6.0) * ap[i];
	dgdtp[i] = lbetap[i] - lbrhsp[i];
      }
    }
  }
}





/* apply linear elliptic operator:  lu = VectorLaplace u */
void LPunctureShift(tL *level, tVarList *vllu, tVarList *vlu)
{
  /* apply boundary conditions */
  set_boundary_elliptic(level, vlu);

  /* compute operator */
  VectorLaplaceFlat(level, vllu, vlu);

  /* synchronize */
  bampi_vlsynchronize(vllu);
}




/* linear Gauss-Seidel:  v = u + (f - Lu)/Lii */
/* elliptic operator:  Lu = VectorLaplace u */
void LPunctureShift_GS(tL *level, tVarList *vlv, tVarList *vlu)
{
  tVarList *vlf = VLPtrEnable1(level, "punctures_ex");

  if (0) {
    prdivider(0);
    printf("before\n");
    prdivider(0);
    prvare(level, "punctures_ex");
    prvare(level, "punctures_sx");
  }

  /* compute operator */
  VectorLaplaceFlatGS(level, vlv, vlu, vlf);

  if (0) {
    prdivider(0);
    printf("after\n");
    prdivider(0);
    prdivider(0);
    prvare(level, "betax");
  }

  /* synchronize and set boundary */
  bampi_vlsynchronize(vlv);
  set_boundary_elliptic(level, vlv);
}




/* compute lapse for maximal slicing */
void PunctureShift(tL *level)
{
  tVarList *vlv  = VLPtrEnable1(level, "betax");
  tVarList *vlf  = VLPtrEnable1(level, "punctures_ex");
  tVarList *vlr  = VLPtrEnable1(level, "punctures_sx");
  int itmax = Geti("punctures_itmax");
  double tol = Getd("punctures_tolerance");
  double normres;

  printf("computing shift for puncture data\n");

  /* set coefficient */
  SetPunctureShiftCoefficient(level);

  /* solve */
  /* multigrid */
  if (0 && Getv("punctures_solver", "multigrid")) 
    multigrid(level, vlv, vlf, vlr, 0, itmax, tol, &normres, 
	      LPunctureShift, LPunctureShift_GS);

  /* bicgstab */
  else if (1 || Getv("punctures_solver", "bicgstab"))
    bicgstab(level, vlv, vlf, vlr, 0, itmax, tol, &normres, 
	     LPunctureShift, DPflatlinear);

  /* unknown solver */
  else 
    errorexit("unknown elliptic solver in punctures/Lapse.c");

  /* set shift */
  set_boundary_elliptic(level, vlv);
  SetPunctureShift(level);

  /* test shift */
  TestPunctureShift(level, vlv, vlf, vlr);

  /* clean up */
  if (Getv("punctures_persist", "no")) {
    VLDisableFree(vlf);
    VLDisableFree(vlr);
  }
  printf("computed shift for puncture data\n");
  prdivider(0);
}

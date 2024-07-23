/* PunctureSolve.c */
/* Bernd Bruegmann 4/97, 1/03 */

/* call elliptic solver
   see PunctureData.c for the setup of coefficients 
*/

#include "bam.h"
#include "punctures.h"


extern tPointList *Punc_InnerList;
extern tPointList *Punc_PhyBoundList;



/* non-linear Gauss-Seidel:  v = u + (f - Lu)/Lii */
/* black hole puncture data
   L(u) = Laplace(delta) u + b (1 + a (1+u))^-7 
*/
void Lpuncture_GS(tL *level, tVarList *vlv, tVarList *vlu)
{
  double cx = 1/(level->dx*level->dx);
  double cy = 1/(level->dy*level->dy);
  double cz = 1/(level->dz*level->dz);
  double cc = -2*(cx + cy + cz);
  double *a = Ptr(level, "punctures_a");
  double *b = Ptr(level, "punctures_b");
  double *f = Ptr(level, "punctures_f");
  double *v, *u;
  double lu, lii;

  u = VLPtr(vlu, 0);
  v = VLPtr(vlv, 0);
  
  if (0) bampi_vlsynchronize(vlu);
  if (0) set_boundary_elliptic(level, vlu);

  /* interior */
  forinner7(level) {
    lu = cc * u[ccc] +
	 cx * (u[mcc] + u[pcc]) +
         cy * (u[cmc] + u[cpc]) +
	 cz * (u[ccm] + u[ccp]) +
         b[ccc] * pow(1 + a[ccc]*(1 + u[ccc]), -7.0);

    lii = cc - 7*b[ccc]*a[ccc]*pow(1 + a[ccc]*(1 + u[ccc]), -8.0);

    v[ccc] = u[ccc] + (f[ccc] - lu)/lii; 
  } endfor;

  /* which is actually needed ? */
  bampi_vlsynchronize(vlv);
  set_boundary_elliptic(level, vlv);
}




/* fourth order version */
void Lpuncture_GS4(tL *level, tVarList *vlv, tVarList *vlu)
{
  double cx = 1/(level->dx*level->dx);
  double cc = -4.0*cx;
  double *a = Ptr(level, "punctures_a");
  double *b = Ptr(level, "punctures_b");
  double *f = Ptr(level, "punctures_f");
  double *v, *u;
  double lu, lii;

  u = VLPtr(vlu, 0);
  v = VLPtr(vlv, 0);

  if (0) bampi_vlsynchronize(vlu);
  if (0) set_boundary_elliptic(level, vlu);

  /* interior */
  forinner19(level) {
    lu = cc * u[ccc] +
         cx * (u[mcc] + u[pcc] + u[cmc] + u[cpc] + u[ccm] + u[ccp])/3.0 +
         cx * (u[cpp] + u[cpm] + u[cmp] + u[cmm] + u[pcp] + u[pcm] + u[mcp] + u[mcm] +
               u[ppc] + u[mpc] + u[pmc] + u[mmc] ) / 6.0 +
         b[ccc] * pow(1 + a[ccc]*(1 + u[ccc]), -7.0) / 2.0 +
         ( b[pcc] * pow(1 + a[pcc]*(1 + u[pcc]), -7.0) +
           b[mcc] * pow(1 + a[mcc]*(1 + u[mcc]), -7.0) +
           b[cpc] * pow(1 + a[cpc]*(1 + u[cpc]), -7.0) +
           b[cmc] * pow(1 + a[cmc]*(1 + u[cmc]), -7.0) +
           b[ccp] * pow(1 + a[ccp]*(1 + u[ccp]), -7.0) +
           b[ccm] * pow(1 + a[ccm]*(1 + u[ccm]), -7.0)) / 12.0;

    lii = cc - (7.0/2.0)*b[ccc]*a[ccc]*pow(1 + a[ccc]*(1 + u[ccc]), -8.0);

    v[ccc] = u[ccc] + (f[ccc] - lu)/lii;
  } endfor;

  /* which is actually needed ? */
  bampi_vlsynchronize(vlv);
  set_boundary_elliptic(level, vlv);
}




/* apply non-linear elliptic operator:  lu = L(u) */
/* black hole puncture data
   L(u) = Laplace(delta) u + b (1 + a (1+u))^-7 
*/
void Lpuncture_nonlin(tL *level, tVarList *vllu, tVarList *vlu)
{
  double cx = 1/(level->dx*level->dx);
  double cy = 1/(level->dy*level->dy);
  double cz = 1/(level->dz*level->dz);
  double cc = -2*(cx + cy + cz);
  double *a = Ptr(level, "punctures_a");
  double *b = Ptr(level, "punctures_b");
  double *lu, *u;

  u = VLPtr(vlu, 0);
  lu = VLPtr(vllu, 0);
  
  if (0) bampi_vlsynchronize(vlu);

  /* apply boundary conditions */
  set_boundary_elliptic(level, vlu);

  /* interior */
  forinner7(level) {
    lu[ccc] = cc * u[ccc] +
	      cx * (u[mcc] + u[pcc]) +
              cy * (u[cmc] + u[cpc]) +
	      cz * (u[ccm] + u[ccp]) +
              b[ccc] * pow(1 + a[ccc]*(1 + u[ccc]), -7.0);
  } endfor;

  bampi_vlsynchronize(vllu);

  if (0) set_boundary_elliptic(level, vllu);
}




/* fourth order version */
void Lpuncture_nonlin4(tL *level, tVarList *vllu, tVarList *vlu)
{
  double cx = 1/(level->dx*level->dx);
  double cc = -4.0*cx;
  double *a = Ptr(level, "punctures_a");
  double *b = Ptr(level, "punctures_b");
  double *lu, *u;

  u = VLPtr(vlu, 0);
  lu = VLPtr(vllu, 0);
 
  if (0) bampi_vlsynchronize(vlu);

  /* apply boundary conditions */
  set_boundary_elliptic(level, vlu);

  /* interior */
  forinner19(level) {
    lu[ccc] = cc * u[ccc] +
              cx * (u[mcc] + u[pcc] + u[cmc] + u[cpc] + u[ccm] + u[ccp])/3.0 +
              cx * (u[cpp] + u[cpm] + u[cmp] + u[cmm] + u[pcp] + u[pcm] + u[mcp] + u[mcm] +
                    u[ppc] + u[mpc] + u[pmc] + u[mmc] ) / 6.0 +
                b[ccc] * pow(1 + a[ccc]*(1 + u[ccc]), -7.0) / 2.0 +
              ( b[pcc] * pow(1 + a[pcc]*(1 + u[pcc]), -7.0) +
                b[mcc] * pow(1 + a[mcc]*(1 + u[mcc]), -7.0) +
                b[cpc] * pow(1 + a[cpc]*(1 + u[cpc]), -7.0) +
                b[cmc] * pow(1 + a[cmc]*(1 + u[cmc]), -7.0) +
                b[ccp] * pow(1 + a[ccp]*(1 + u[ccp]), -7.0) +
                b[ccm] * pow(1 + a[ccm]*(1 + u[ccm]), -7.0)) / 12.0;
  } endfor;

  bampi_vlsynchronize(vllu);

  if (0) set_boundary_elliptic(level, vllu);
}




/* apply linearized elliptic operator:  lu = L(u) */
/* black hole puncture data
   L(u) = Laplace(delta) u - 7 a b (1 + a (1+u0))^-8) u 
*/
void Lpuncture_linear(tL *level, tVarList *vllu, tVarList *vlu)
{
  int i;
  double *a = Ptr(level, "punctures_a");
  double *b = Ptr(level, "punctures_b");
  double *u0 = Ptr(level, "punctures_u");
  double *lu, *u;
  double cx = 1/(level->dx*level->dx);
  double cy = 1/(level->dy*level->dy);
  double cz = 1/(level->dz*level->dz);
  double cc = -2*(cx + cy + cz);

  u = VLPtr(vlu, 0);
  lu = VLPtr(vllu, 0);
  
  /* apply boundary conditions */
  set_boundary_elliptic(level, vlu);

  /* interior */
  forinner7(level) {
    lu[ccc] = cc * u[ccc] +
	      cx * (u[mcc] + u[pcc]) +
              cy * (u[cmc] + u[cpc]) +
	      cz * (u[ccm] + u[ccp])
      - 7 * a[ccc] * b[ccc] * pow(1 + a[ccc]*(1 + u0[ccc]), -8.0) * u[ccc];
  } endfor;

  bampi_vlsynchronize(vllu);
}




/* compute linearized elliptic operator:  lu = L(u) , for PointList */
/* black hole puncture data
   L(u) = Laplace(delta) u - 7 a b (1 + a (1+u0))^-8) u 
*/
void Lpuncture_lin_PointList(tPointList *plist,
                             tVarList *vllu, tVarList *vlu)
{
  errorexit("BOX: implement Lpuncture_lin_PointList ... look at older BAM versions");
}




/* compute linearized elliptic operator using:  Lpuncture_lin_PointList
   should do same as:  Lpuncture_linear
*/
void Lpuncture_lin(tL *level, tVarList *vllu, tVarList *vlu)
{
  Lpuncture_lin_PointList(Punc_InnerList, vllu,vlu);
}




/* call solver for puncture form of Hamiltonian constraint */
void PunctureSolve(tL *level)
{
  tG *g = level->grid;
  tVarList *vlu, *vlv, *vlf, *vlr;
  tVarList *vlc = vlalloc(level);
  int itmax = Geti("punctures_itmax");
  double tol = Getd("punctures_tolerance");
  double normres, normresnonlin, u;
  int i, inewton, l;
  int usel;

  double distance;
  double pmass;

  /* storage */
  for (l = g->lmax; l >= level->l; l--) {
    vlu = VLPtrEnable1(g->level[l], "punctures_u");
    vlv = VLPtrEnable1(g->level[l], "punctures_v");
    vlf = VLPtrEnable1(g->level[l], "punctures_f");
    vlr = VLPtrEnable1(g->level[l], "punctures_r");
  }

  /* make variable list for coefficients */
  vlpush(vlc, Ind("punctures_a"));
  vlpush(vlc, Ind("punctures_b"));

  /* no solve */
  if (Getv("punctures_solver", "none"))
    printf("PunctureSolve: no elliptic solve\n");

  /* jacobi */
  else if (Getv("punctures_solver", "jacobi"))
    punctures_solver(level, vlu, vlf, vlr, vlc, itmax, tol, &normres, 
		     Lpuncture_nonlin, DPflatlinear);

  /* multigrid */
  else if (Getv("punctures_solver", "multigrid") || 
	   Getv("punctures_solver", "gaussseidel")) {
    if (Getv("punctures_fourth_order","yes")) 
      punctures_solver(level, vlu, vlf, vlr, vlc, itmax, tol, &normres,
                     Lpuncture_nonlin4, Lpuncture_GS4);
    else
      punctures_solver(level, vlu, vlf, vlr, vlc, itmax, tol, &normres, 
		     Lpuncture_nonlin, Lpuncture_GS);
  }

  /* bicgstab */
  else if (Getv("punctures_solver", "bicgstab")) {
    double *f = Ptr(level, "punctures_f");
    double *u = Ptr(level, "punctures_u");
    double *v = Ptr(level, "punctures_v");

    /* Newton iteration: 
       there now is general purpose Newton solver in src/elliptic, used below
    */
    for (inewton = 0; inewton < itmax; inewton++) {
      Lpuncture_nonlin(level, vlf, vlu);
      forallpoints(level, i) f[i] = - f[i];
      normresnonlin = bampi_allreduce_allnorm2(vlf);
      printf("Newton %d:  norm(res) = %10.3e\n", inewton, normresnonlin);
      if (normresnonlin <= tol) break;

      forallpoints(level, i) v[i] = 0;

      /* solve linear equation */
      bicgstab(level, vlv, vlf, vlr, vlc, 20, tol/10, &normres, 
	       Lpuncture_linear, DPflatlinear);

      set_boundary_elliptic(level, vlv);

      forallpoints(level, i) u[i] += v[i];
    }
  }

  /* Newton solver using specified linear solver */
  else if (Getv("punctures_solver", "Newton"))
  {
    int    linSolver_itmax  = Geti("punctures_linSolver_itmax");
    double linSolver_tolFac = Getd("punctures_linSolver_tolFac");
    int pr=1; /* Newton verbose */

    Newton(level, Lpuncture_nonlin, Lpuncture_lin, vlu, vlf,
           itmax, tol, &normresnonlin, pr,
           linear_solver, DPflatlinear, vlv, vlr,
           vlc,  linSolver_itmax, linSolver_tolFac);
  }

  /* spectral puncture solver, reads from file if available */
  else if (Getv("punctures_solver", "spectral"))
    PuncturesPS(level);

  /* unknown solver */
  else 
    errorexit("unknown elliptic solver in PunctureSolve.c");

  /* info */
  distance = pow(CBL[0][0]-CBL[1][0],2)
    + pow(CBL[0][1]-CBL[1][1],2)
    + pow(CBL[0][2]-CBL[1][2],2);
  distance = sqrt(distance);

  for (i = 0; MBL[i] != 0; i++) {
     
    usel = g->lmaxpunc[i];
    
    u = interpolate_xyz_scalar(g->level[usel], 
                               CBL[i][0], CBL[i][1], CBL[i][2], 
                               Ind("punctures_u"), Geti("order_RP"),LAGRANGE);
    if (i == 0) {
      if (MBL[1])
	pmass = MBL[0]*(1.0 + u + MBL[1]/(2.0*distance));
      else
	pmass = MBL[0]*(1.0 + u);
    }
    if (i == 1)
      pmass = MBL[1]*(1.0 + u + MBL[0]/(2.0*distance));

    printf("ADM mass at puncture%d = %.15g   (u%d = %.15g)\n",
           i+1, pmass, i+1, u);
  }

  /* set boundary for final result */
  if (!Getv("punctures_boundary","none"))
    for (l = g->lmin; l <= g->lmax; l++)
      set_boundary_elliptic(g->level[l], vlu);

  /* fill in overlap (disable for debugging) */
  if (1) restrict_prolong_grid(g, vlu);

  /* clean up */
  vlfree(vlc);
  if (Getv("punctures_persist", "no")) {
    for (l = level->l; l <= g->lmax; l++) {
      vlf->level = g->level[l];
      vlr->level = g->level[l];
      vldisable(vlf);
      vldisable(vlr);
    }
    vlfree(vlf);
    vlfree(vlr);
  }
}





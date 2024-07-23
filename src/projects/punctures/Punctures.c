/* Punctures.c */
/* Bernd Bruegmann 6/02 */
/* based on brbr.c, Bernd Bruegmann, 11/97 */


/* Brill-Lindquist data + Bowen-York momentum and spin 
   = Brandt-Bruegmann data (gr-qc/9703066)

   gphys_ij = phi^4 delta_ij
   Kphys_ij = phi^-2 K_ij,    K_ij = Sum K_ij(P,S)

   1/a = Sum m/(2r)

   psi = 1 + 1/a    -> used for time independent conformal rescaling,
                       identical to Brill-Lindquist

   phi = psi + u    -> u is determined numerically for solution to
                       Hamiltonian constraint

   In terms of u, the Hamiltonian constraint becomes

   H = Laplace(flat) u + b (1 + a (1+u))^-7,   

   where b = 1/8 a^7 K_ij K^ij

   Note that a, b, and u are regular on R^3.

   Note also that we store
     g_ij = psi^-4 gphys_ij
     K_ij =        Kphys_ij
   although for uniformity we should probably store
     K_ij = psi^-4 Kphys_ij
*/


#include "bam.h"
#include "punctures.h"


double MBL[N_MWBL+1], CBL[N_MWBL+1][3], PBL[N_MWBL+1][3], SBL[N_MWBL+1][3];

int fisheye = 0;
double fisheyedx, fisheyedy, fisheyedz;


/* global point lists */
tPointList *Punc_InnerList;
tPointList *Punc_PhyBoundList;
// tPointList *Punc_SymBoundList;




/* set the puncture pars from the fit we produced in gr-qc/0307027 */
int ComputePunctureParametersFromFit(tL *level)
{
  double P, m, Omega, x;
  
  x = Getd("punctures_from_fit_x");
  
  m=0.5 - 0.023870808487626096/pow(x,3.5) + 0.0973014163688868/pow(x,3) -
    0.1393039760269027/pow(x,2.5) + 0.0906936980338675/pow(x,2) -
    0.028008993306827646/pow(x,1.5) - 0.06685488924886536/x;
          
  P=0.10147093802488151/pow(x,3.5) - 0.4070333515125441/pow(x,3) +
    0.5541037717762235/pow(x,2.5) - 0.3176472123574175/pow(x,2) +
    0.2742503126543588/pow(x,1.5) - 0.012334135951504049/x +
    0.1767766952966369/sqrt(x);
              
  Omega=-0.023608084059311016/pow(x,4.5) + 0.09699030947727742/pow(x,4) -
         0.175001025703614/pow(x,3.5) + 0.27549220191298884/pow(x,3) -
         0.34761658037384446/pow(x,2.5) + 0.015800634566454162/pow(x,2) +
         1/(2.*sqrt(2)*pow(x,1.5));
  
  Setd("bhmass1", m);
  Setd("bhmass2", m);

  Setd("bhx1", 0);
  Setd("bhy1", x);
  Setd("bhz1", 0);
  Setd("bhx2", 0);
  Setd("bhy2", -x);
  Setd("bhz2", 0);

  Setd("bhpx1", -P);
  Setd("bhpy1", 0);
  Setd("bhpz1", 0);
  Setd("bhpx2", P);
  Setd("bhpy2", 0);
  Setd("bhpz2", 0);
  
  if (GetsLax("corotation_omega"))
    Setd("corotation_omega", Omega);

  return 0;
}


/* Set the puncture parameters from Kidder's formula */
int ComputePunctureParametersFromPN(tL *level)
{
    double m1,m2,sx1,sy1,sz1,sx2,sy2,sz2,r,y1;
    double a1,a2;
    double sdirnx1,sdirny1,sdirnz1;
    double sdirnx2,sdirny2,sdirnz2;
    double LS1,LS2,Sk,SA,SB,SC,SD,momx,momy,momz,posx,posy,posz,s1s2;
    double m1sqr, m2sqr;
    double M,mu;

    printf("Hello.\n"); fflush(stdout);
 
    sdirnx1 = 0;
    sdirny1 = 0;
    sdirnz1 = 0;
    sdirnx2 = 0;
    sdirny2 = 0;
    sdirnz2 = 0;

/* Be careful to make sure that these are the puncture masses we want */
    m1 = Getd("punctures_from_PN_m1");
    m2 = Getd("punctures_from_PN_m2");

    m1sqr = m1*m1;
    m2sqr = m2*m2;
    M = m1 + m2;
    mu = m1*m2/(m1 + m2);

    sx1 = Getd("bhsx1")/m1sqr;
    sy1 = Getd("bhsy1")/m1sqr;
    sz1 = Getd("bhsz1")/m1sqr;
    sx2 = Getd("bhsx2")/m2sqr;
    sy2 = Getd("bhsy2")/m2sqr;
    sz2 = Getd("bhsz2")/m2sqr;

    a1 = sqrt(sx1*sx1 + sy1*sy1 + sz1*sz1);
    a2 = sqrt(sx2*sx2 + sy2*sy2 + sz2*sz2);

    r = M*Getd("punctures_from_PN_D");

    /* Assume that orbital angular momentum is in the Z direction */
    if (fabs(a1) > 1e-7) LS1 = sz1 / a1; else LS1 = 0;
    if (fabs(a2) > 1e-7) LS2 = sz2 / a2; else LS2 = 0;

    if (fabs(a1) > 1e-7 && fabs(a2) > 1e-7) s1s2 = (sx1*sx2 + sy1*sy2 + sz1*sz2)/(a1*a2); 
    else s1s2 = 0;
    
    if (fabs(a1) > 1e-7) {
	sdirnx1 = sx1 / a1;
	sdirny1 = sy1 / a1;
	sdirnz1 = sz1 / a1;
    }
    else {
        sdirnx1 = 0;
        sdirny1 = 0;
        sdirnz1 = 1.0;
    } 
    if (fabs(a2) > 1e-7) {
	sdirnx2 = sx2 / a2;
	sdirny2 = sy2 / a2;
	sdirnz2 = sz2 / a2;
    }
    else {
        sdirnx2 = 0;
        sdirny2 = 0;
        sdirnz2 = 1.0;
    }

    Sk = mu*sqrt(M*r);
    SA = 1 - 0.25*(a1*LS1*((8*pow(m1,2))/pow(M,2) + (7*mu)/M) + a2*LS2*((8*pow(m2,2))/pow(M,2) + (7*mu)/M))*pow(M/r,1.5) + (2*M)/r +
      (pow(M,2)*(0.5*(5 - (9*mu)/M) - (0.75*a1*a2*mu*(-3*LS1*LS2 + s1s2))/M))/pow(r,2);
    SB = 0.25*pow(M/r,1.5)*(a1*((4*pow(m1,2))/pow(M,2) + mu/M)*sdirnx1 + a2*((4*pow(m2,2))/pow(M,2) + mu/M)*sdirnx2);
    SC = 0.25*pow(M/r,1.5)*(a1*((4*pow(m1,2))/pow(M,2) + mu/M)*sdirny1 + a2*((4*pow(m2,2))/pow(M,2) + mu/M)*sdirny2);
    SD = 0.25*pow(M/r,1.5)*(a1*((4*pow(m1,2))/pow(M,2) + mu/M)*sdirnz1 + a2*((4*pow(m2,2))/pow(M,2) + mu/M)*sdirnz2);

    momx = Sk*(sqrt(SA*SA - SB*SB - SC*SC) + SD) / r;
    momy = 0;
    momz = 0;

    posx = 0;
    posy = 1;
    posz = 0;

    y1 = r*(m2/m1)/(1 + m2/m1);

    Setd("bhmass1",m1*sqrt(1.0 - a1*a1));
    Setd("bhmass2",m2*sqrt(1.0 - a2*a2));

    Setd("bhx1",posx*(r*m2/m1)/(1+m2/m1));
    Setd("bhy1",posy*(r*m2/m1)/(1+m2/m1));
    Setd("bhz1",posz*(r*m2/m1)/(1+m2/m1));
    Setd("bhx2",posx*(y1 - r));
    Setd("bhy2",posy*(y1 - r));
    Setd("bhz2",posz*(y1 - r));

    Setd("bhpx1",-momx);
    Setd("bhpx2",momx);
    Setd("bhpy1",0);
    Setd("bhpy2",0);
    Setd("bhpz1",momz);
    Setd("bhpz2",momz);

    Setd("bhsx1",sx1*m1*m1);
    Setd("bhsy1",sy1*m1*m1);
    Setd("bhsz1",sz1*m1*m1);
    Setd("bhsx2",sx2*m2*m2);
    Setd("bhsy2",sy2*m2*m2);
    Setd("bhsz2",sz2*m2*m2);

    printf("bhp1 = (%e,%e,%e).\n",-momx,momy,momz);
    printf("bhs1 = (%e,%e,%e).\n",sx1*m1*m1,sy1*m1*m1,sz1*m1*m1);
    printf("bhx = (%e,%e,%e).\n",posx*(r*m2/m1)/(1+m2/m1),posy*(r*m2/m1)/(1+m2/m1),posz*(r*m2/m1)/(1+m2/m1));
    printf("bhx = (%e,%e,%e).\n",posx*(y1 - r),posy*(y1 - r),posz*(y1 - r));

    return 0;

}



/* read puncture parameters */
void ReadPunctureParameters(int *ptrKnonzero)
{
  int i, j;
  char *coord[3] = {"x", "y", "z"};
  char t[1024];
  int Knonzero = 0;

  /* get parameters */
  /* old style, change to multiple arguments per line */
  for (i = 0; i < N_MWBL; i++) {
    sprintf(t, "mass%d", i+1);
    MBL[i] = Getd(t);
    if (MBL[i] == 0.0) break;
    for (j = 0; j < 3; j++) {
      sprintf(t, "p%s%d", coord[j], i+1);
      CBL[i][j] = Getd(t);
    }
    for (j = 0; j < 3; j++) {
      sprintf(t, "m%s%d", coord[j], i+1);
      PBL[i][j] = Getd(t);
      if (PBL[i][j] != 0) Knonzero = 1;
    }
    for (j = 0; j < 3; j++) {
      sprintf(t, "s%s%d", coord[j], i+1);
      SBL[i][j] = Getd(t);
      if (SBL[i][j] != 0) Knonzero = 1;
    }
  }

  /* print */
  prdivider(0);
  if (MBL[0] != 0)
    printf("Computing black hole initial data (puncture data):\n");
  for (i = 0; i < N_MWBL; i++) {
    if (MBL[i] == 0.0)
      break;
    printf("m%d  %10.3e\n", i+1, MBL[i]);
    printf("x%d  %10.3e  y%d  %10.3e  z%d  %10.3e\n", 
	   i+1, CBL[i][0], i+1, CBL[i][1], i+1, CBL[i][2]);
    printf("px%d %10.3e  py%d %10.3e  pz%d %10.3e\n",
	   i+1, PBL[i][0], i+1, PBL[i][1], i+1, PBL[i][2]);
    printf("sx%d %10.3e  sy%d %10.3e  sz%d %10.3e\n\n",
	   i+1, SBL[i][0], i+1, SBL[i][1], i+1, SBL[i][2]);
  }
  if (Getv("physics", "punctures") && i == 0) {
    printf("There are 0 black holes: gij = deltaij, Kij = 0, psi = 1\n"); 
  }

  /* never reached, should make this do something again */
  if (i > N_MWBL) {
    printf("\n\npuncture data: there currently is a hard limit");
    printf(" of at most %d punctures\n\n", N_MWBL);
  } 

  /* return whether K is nonzero */
  *ptrKnonzero = Knonzero;
}




/* pick elliptic solver */
void PickEllipticSolver(void)
{
  /* pick elliptic solver */
  if (Getv("punctures_solver", "none"))
    punctures_solver = 0;
  else if (Getv("punctures_solver", "jacobi"))
    punctures_solver = jacobi;
  else if (Getv("punctures_solver", "gaussseidel"))
    punctures_solver = gaussseidel;
  else if (Getv("punctures_solver", "bicgstab"))
    punctures_solver = bicgstab;
  else if (Getv("punctures_solver", "multigrid"))
    punctures_solver = multigrid;
  else if (Getv("punctures_solver", "Newton"))
  {
    if (Getv("punctures_linSolver", "bicgstab"))
      linear_solver=bicgstab;
    /*
    else if (Getv("punctures_linSolver", "HYPRE"))
      linear_solver=Punc_hypreWrapper;
    */
    else
      errorexit("Unknown punctures_linSolver in Punctures.c");
  }
  else if (Getv("punctures_solver", "spectral"))
    punctures_solver = 0;
  else 
    errorexit("Unknown elliptic solver in Punctures.c");
}




/* set coefficients for puncture form of Hamiltonian constraint
   a = (Sum m/(2r))^-1
   b = 1/8 a^7 K_ij K^ij
*/
void SetPunctureCoefficients(tL *level)
{
  double *a = Ptr(level, "punctures_a");
  double *b = Ptr(level, "punctures_b");
  double *psi = Ptr(level, "adm_psi");
  double *K11, *K12, *K13, *K22, *K23, *K33;
  double K2;
  int i;

  i = Ind("adm_Kxx");
  K11 = level->v[i++];
  K12 = level->v[i++];
  K13 = level->v[i++];
  K22 = level->v[i++];
  K23 = level->v[i++];
  K33 = level->v[i++];

#if 0
  if (fisheye) {
    if (0) printf("--> calling set_laplace()\n");
    set_laplace(g);
    setrangeall(g);
    if (0) printf("--> calling set_K2()\n");
    /* from standalone maximal: N = - K_ij K^ij 
       call this utility function since fisheye has a nontrivial g^ij
    */
    set_K2(g, N, M);
  }
#endif

  forallpoints(level, i) {

    a[i] = 1/(psi[i]-1);

    if (fisheye) {
      b[i] = 0.125 * pow(a[i], 7.0) * (-b[i]);
    } else {
      K2 = K11[i]*K11[i] + K22[i]*K22[i] + K33[i]*K33[i] + 
	   2*K12[i]*K12[i] + 2*K13[i]*K13[i] + 2*K23[i]*K23[i];
      b[i] = 0.125 * pow(a[i], 7.0) * K2;
    }
  }
}




/* iterate lapse computation until M_K = M_ADM */
double IterateLapseUntil_MK_equals_MADM(tL *level)
{
  tMasses masses;
  
  double c1, c2, cm;
  double MADMminusMK1, MADMminusMK2, MADMminusMKm;

  c1=-2.0; 
  c2=0.0;
  
  Setd("punctures_lapse_at_puncture",c1);
  PunctureMaximalSlicing(level);
  MADMminusMK1= PunctureMass(level,0, &masses);
  
  Setd("punctures_lapse_at_puncture",c2);
  PunctureMaximalSlicing(level);
  MADMminusMK2= PunctureMass(level,0, &masses);
  
  while(MADMminusMK1*MADMminusMK2<0)
  {
    cm=c1 - MADMminusMK1*(c2-c1)/(MADMminusMK2-MADMminusMK1);
    Setd("punctures_lapse_at_puncture",cm);
    PunctureMaximalSlicing(level);
    MADMminusMKm = PunctureMass(level,0, &masses);
    if( fabs(MADMminusMKm) < Getd("punctures_lapse_MK_MADM_diff") )
      break;
    
    if(MADMminusMKm*MADMminusMK1>0) 
    {
      c1=cm;
      MADMminusMK1 = MADMminusMKm; 
    }
    else
    {
      c2=cm;
      MADMminusMK2 = MADMminusMKm;
    }
  }
  return cm;		
}




/* output punctures_lapse_at_puncture */
void Output_lapse_at_puncture(tL *level)
{
  FILE *out;
  char filename[10000];

  if (processor0)
  {
    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_lapse_at_puncture", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%s\n", Gets("punctures_lapse_at_puncture"));
    fclose(out);
  }
}




/* Compute puncture data.
   Note that the amr driver in main.c calls the initial data routines
   on levels lmin, lmin+1, ..., lmax, which is good for the simplest case
   of analytic data. However, the elliptic solve works on all levels
   simultaneously, so we wait for lmax and do the loop over levels locally.
*/
int PunctureData(tL *level)
{
  tG *g = level->grid;
  int lmin = g->lmin;
  int lmax = g->lmax;
  int Knonzero;
  tMasses masses;
  int sequence = Getv("punctures_MK_eq_MADM_Sequence", "yes");
  int MK_equals_MADM = Getv("punctures_lapse_at_puncture", "MK_equals_MADM");
  int i, l;
  
  enablevar(level,Ind("adm_psi"));
  enablevar(level,Ind("adm_dpsiopsix"));
  enablevar(level,Ind("adm_ddpsiopsixx"));
  
  /* wait until we are called on finest level */
  /* we call the elliptic solver to compute also mixed binaries, here we use the call
     at the coarsest grid*/

  if ((level->l < lmax) && !(Getv("punctures_solv_included", "yes"))) return 0;
  if ((level->l > lmin) &&  (Getv("punctures_solv_included", "yes"))) return 0;


  /* read parameters */
  ReadPunctureParameters(&Knonzero);
  PickEllipticSolver();
  
  /* set information for boundary conditions */
  if (Getv("punctures_boundary", "robin")) {
    double sf = Getd("punctures_shift_falloff");
    double uf = Getd("punctures_robin_falloff");
    VarNameSetBoundaryInfo("punctures_u", 0, uf, 0);
    VarNameSetBoundaryInfo("punctures_v", 0, uf, 0);
    VarNameSetBoundaryInfo("betax", 0, sf, 0);
    VarNameSetBoundaryInfo("betay", 0, sf, 0);
    VarNameSetBoundaryInfo("betaz", 0, sf, 0);
    Appends("boundary", "robin");

    /* only coarsest box is allowed to have physical boundary for now */
    find_robin_normal(g->level[lmin]);
  } else {
    VarNameSetBoundaryInfo("punctures_u", 0, 0, 0);
    VarNameSetBoundaryInfo("punctures_v", 0, 0, 0);
    VarNameSetBoundaryInfo("betax", 0, 0, 0);
    VarNameSetBoundaryInfo("betay", 0, 0, 0);
    VarNameSetBoundaryInfo("betaz", 0, 0, 0);
  }

  /* now do initialization that can be done for each level separately */ 
  for (l = lmin; l <= lmax; l++) {
    level = g->level[l];

    /* temporary: special initialization for point lists and HYPRE */
    PunctureSpecialInit(level);

    /* enable storage */
    enablevar(level, Ind("punctures_u"));
    enablevar(level, Ind("punctures_a"));
    enablevar(level, Ind("punctures_b"));

    /* sequences are done separately in Sequences.c, needs FIX for amr */
    if (sequence) {
      if (lmax > 0) 
	errorexit("PunctureData: sequences not implemented for amr yet");
      Make_MK_eq_MADM_Sequence(level);
    }

    /* no sequence default */
    if (!sequence) {

      /* set Brill-Lindquist data */
      SetMWBL(level);

      /* set Bowen-York extrinsic curvature */
      SetBYK(level);

      /* precompute coefficients for elliptic solver */
      SetPunctureCoefficients(level);
    }
  }

  /* call elliptic solver with coarsest level as argument */
  if (!sequence) {
    level = g->level[lmin];
    PunctureSolve(level);
  }


  /* lapse, needs FIX for AMR */
  if (!sequence) {
    if (Getv("punctures_lapse", "maximal"))
      {
	if (lmax > 0) 
	  errorexit("PunctureData: lapse/shift not implemented for amr yet");
	
	if (MK_equals_MADM)
	  IterateLapseUntil_MK_equals_MADM(level);
	else
	  PunctureMaximalSlicing(level);
	
	if (Getv("punctures_shift", "thinsandwich"))
	  PunctureShift(level);
      }        
    else if (Getv("punctures_lapse", "thinsandwich")) {
      if (lmax > 0) 
	errorexit("PunctureData: lapse/shift not implemented for amr yet");
      if (Getv("punctures_shift", "thinsandwich")) /* unfinished */1;
    }

    /* mass */
    if (!Getv("punctures_mass", "no")) 
      PunctureMass(level, 1, &masses);
  } 

  /* finish each level separately */ 
  for (l = lmin; l <= lmax; l++) {
    level = g->level[l];

    /* given conformal data, set bam's adm data */
    if (Knonzero)
      ConfToPhys(level, "punctures_u");
    else
      ConfToPhys(level, "");


  if (Getv("punctures_solv_included", "yes")) { 

      enablevar(level, Ind("incl_punctures_u"));
      double *inc_u      = Ptr(level, "incl_punctures_u");
      double *ps_u       = Ptr(level, "punctures_u");
       forallpoints(level, i) { 
        inc_u[i] = ps_u[i] ; 
        }
    }

    /* clean up */
    if (Getv("punctures_persist", "no")) {
      disablevar(level, Ind("punctures_u"));
      disablevar(level, Ind("punctures_a"));
      disablevar(level, Ind("punctures_b"));
    }
    PunctureCleanup(level);

    if(Getv("punctures_lapse", "psiBL^(-2)"))  
      Set_alpha_psip(level,-2);
    if(Getv("punctures_lapse", "psipower"))
      Set_alpha_psip(level,Getd("punctures_alp_psipower"));
    if(Getv("punctures_lapse", "rtoN_atPunc"))
      Set_alpha_rtoN_atPunc(level);
    if(Getv("punctures_Absorbpsi", "yes"))
      Absorb_psi(level);
  }

  if (MK_equals_MADM) {
    /* Output punctures_lapse_at_puncture as computed by 
       IterateLapseUntil_MK_equals_MADM(level);
       and then reset it to "MK_equals_MADM"             
       needs FIX for amr
    */
    Output_lapse_at_puncture(g->level[lmin]);
    Sets("punctures_lapse_at_puncture", "MK_equals_MADM");
  }
  
  prdivider(0);
  return 0;
}




/* Cleanup function which is run in POST_INITIALDATA: 
   BB: No! Since PunctureData is called repeatedly on different levels,
       we have to clean up right after we are done (or save the info
       for a full level stack cleanup)
*/
void PunctureCleanup(tL *level)
{
  if (!Getv("punctures_solver", "Newton")) return;

  if (0) prdivider(0);
  if (0) printf("PunctureCleanup: Clean up of Punctures\n");
      
  /* free memory for pointlists */
  FreePointList(Punc_InnerList);
  FreePointList(Punc_PhyBoundList);
  Punc_InnerList = NULL;
  Punc_PhyBoundList = NULL;
  //  FreePointList(Punc_SymBoundList);
            
  /* HYPRE */
  //Punc_hypreCleanup(level);
}




/* special initialization */
void PunctureSpecialInit(tL *level)
{
  int i;

  if (Getv("punctures_solver", "Newton")) {
    FreePointList(Punc_InnerList);
    FreePointList(Punc_PhyBoundList);
    Punc_InnerList = AllocatePointList(level);
    Punc_PhyBoundList = AllocatePointList(level);
    forinner1(level,i)
      AddToPointList(level, Punc_InnerList, i);
    forallpoints(level, i)
      if (boundaryflag(level, i) == PHYBOUND)
	AddToPointList(level, Punc_PhyBoundList, i);
  }    
}


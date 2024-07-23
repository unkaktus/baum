/* z4.c */
/* Wolfgang Tichy  4/2004 */

#include "bam.h"
#include "z4.h"

/* infos about the evolve varlist*/
int METRIC_z4_INDX_VAR;



/* evolve in the interior and 
   for those boundary points set by special evolution
*/
void z4_evolve(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  timer_start(ucur->level, "z4_rhs");

  z4_rhs_movpunc_N(unew, upre, c, ucur);

  timer_stop(ucur->level, "z4_rhs");
}





/* set boundary of Gammas
   - Cactus uses third order extrapolation
   - so far have only linear (second order) extrapolation
   - could also try one sided derivatives
*/
void z4_setbound_G(tL *level, int igi, int iG)
{
  int i;
  
  for (i = 0; i < 3; i++)
    set_boundary_extrapolate(level, iG+i);
}




/* compute trA and detg for z4 variables
   enforce trA = 0 and detg = 1 if requested:
     gtilde_ij = detgtilde^(-1/3) gtilde_ij
     Atilde_ij = Atilde_ij - 1/3 gtilde_ij gtilde^kl Atilde_kl
*/
int z4_algcon_wrapper(tL *level) 
{
  if (0) printf("z4: enforced algebraic constraints on l%d\n", level->l);
  
  tVarList *vl = vlalloc(level);

  vlpush(vl, Ind("bssn_gxx"));
  vlpush(vl, Ind("bssn_Axx"));
  vlpush(vl, Ind("bssn_trA"));
  vlpush(vl, Ind("bssn_detg"));

  /* do we actually want to store the result? */
  if (Getv("z4_register_algcon", "store"))
    enablevarlist(vl);

  /* compute and enforce */
  z4_algcon(vl);

  /* done */
  vlfree(vl);
  return 0;
}




/* derive ADM variables from z4 variables
   g_ij = exp(4 phi) gtilde_ij
   A_ij = exp(4 phi) Atilde_ij
   K_ij = A_ij + 1/3 K g_ij
*/
int z4_to_adm(tL *level) 
{
  tVarList *vl = vlalloc(level);
  vlpush(vl, Ind("bssn_gxx"));
  vlpush(vl, Ind("bssn_chi"));
  vlpush(vl, Ind("bssn_Axx"));
  vlpush(vl, Ind("bssn_K"));
  z4_to_adm_vl(vl,vl,0,vl);
  vlfree(vl);
  
  return 0;
}

void z4_to_adm_vl(tVarList *unew, tVarList *upre, double c, tVarList *ucur) 
{
  if (0) printf("  deriving ADM variables from Z4 variables\n");
  tL* level = ucur->level;
  
  int ig      = Ind("adm_gxx");
  int iK      = Ind("adm_Kxx");
  int igtilde = Indvl("bssn_gxx",ucur);
  int iAtilde = Indvl("bssn_Axx",ucur);
  
  double *trK = level->v[Indvl("bssn_K",ucur)];
  double *chi = level->v[Indvl("bssn_chi",ucur)];
  double chipsipower = Getd("z4_chi_psipower");

  bampi_openmp_start
  forallpoints_ijk_openmp(level) {
    double psi  = pow(chi[ijk],1./chipsipower);
    double psi4 = pow(psi,4.);
    
    for (int n = 0; n < 6; n++) {
      double *g = level->v[ig+n];
      double *K = level->v[iK+n];
      double *gtilde = level->v[igtilde+n];
      double *Atilde = level->v[iAtilde+n];
      
      g[ijk] = psi4 * gtilde[ijk];
      K[ijk] = psi4 * Atilde[ijk] +  1./3. * trK[ijk] * g[ijk] ;
    }
    
  } endfor_ijk_openmp;
  bampi_openmp_stop
}




/* derive z4 variables from ADM variables
   p  = detgbar^(-1/3) 
   p0 = psi^(-4)

   gtilde_ij = p gbar_ij
   Ktilde_ij = p p0 K_ij

   phi = - log(p) / 4
   K   = gtildeinv^ij Ktilde_ij
   Atilde_ij = Ktilde_ij - gtilde_ij K / 3

   G^i = - del_j gtildeinv^ji
*/
void adm_to_z4(tL *level) 
{
  printf("  deriving Z4 variables from ADM variables\n");

  /* need temporary storage for derivative of inverse metric */
  enablevar(level, Ind("bssn_ginvtmpxx"));
  
  tVarList* vl = vlalloc(level);
  vlpush(vl, Ind("adm_gxx"));
  vlpush(vl, Ind("adm_Kxx"));
  vlpush(vl, Ind("bssn_gxx"));
  vlpush(vl, Ind("bssn_chi"));
  vlpush(vl, Ind("bssn_Axx"));
  vlpush(vl, Ind("bssn_K"));
  vlpush(vl, Ind("bssn_Gx"));
  vlpush(vl, Ind("bssn_ginvtmpxx"));
  vlpush(vl, Ind("z4_Theta"));
  
  /* compute */
  z4_init(vl);
  
  /* free varlist */
  vlfree(vl);
  
  /* set boundary for Gammas */
  z4_setbound_G(level, Ind("bssn_ginvtmpxx"), Ind("bssn_Gx"));
  disablevar(level, Ind("bssn_ginvtmpxx"));
  bampi_synchronize(level, Ind("bssn_Gx"));

}


/* initialize z4 after initial data has been computed in POST_INITIALDATA
*/
int z4_startup(tL *level)
{
  int i;
  double vgauge;

  printf("Initializing z4:\n");

  /* get ptr to evolution variables and add bssnvars
  DO NOT free this list, it is stored globally */
  tVarList* z4vars = get_evolve_vlregister(level);
  
  if (level->l==0) {
    
    /* like in MATTER, we better remember where the z4 varaiables are 
       inside the evolution varlist, do it only once */
    METRIC_z4_INDX_VAR = z4vars->n;
    
    /* set gauge speed for lapse and related quantities */
    if (Getv("z4_lapse", "1+log"))
      vgauge = sqrt(Getd("z4_lapseharmonicf"));
    else if (Getv("z4_lapse", "1+log2"))
      vgauge = sqrt(4.0/3.0);
    else
      vgauge = 1;
  
    /* set boundary information for z4 evolution: 
      farlimit, falloff, propagation speed 
    */
    VarNameSetBoundaryInfo("bssn_gxx", 1, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_gxy", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_gxz", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_gyy", 1, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_gyz", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_gzz", 1, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Axx", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Axy", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Axz", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Ayy", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Ayz", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Azz", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_K",   0, 1, vgauge);
    VarNameSetBoundaryInfo("bssn_Gx",  0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Gy",  0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_Gz",  0, 1, 1.0);
    VarNameSetBoundaryInfo("z4_Theta", 0, 1, 1.0);
    VarNameSetBoundaryInfo("bssn_chi", 1, 1, vgauge);
    VarNameSetBoundaryInfo("alpha",    1, 1, vgauge);
    VarNameSetBoundaryInfo("betax",    0, 1, 1.0); // 2./sqrt(3.)); 
    VarNameSetBoundaryInfo("betay",    0, 1, 1.0); // 2./sqrt(3.));
    VarNameSetBoundaryInfo("betaz",    0, 1, 1.0); // 2./sqrt(3.));
    VarNameSetBoundaryInfo("betadotx", 0, 1, 1.0);
    VarNameSetBoundaryInfo("betadoty", 0, 1, 1.0);
    VarNameSetBoundaryInfo("betadotz", 0, 1, 1.0);

  
    /* create a variable list for z4 evolutions 
      note that we include lapse and shift directly
    */
    vlpush(z4vars, Ind("bssn_gxx"));
    vlpush(z4vars, Ind("bssn_chi"));
    vlpush(z4vars, Ind("bssn_Axx"));
    vlpush(z4vars, Ind("bssn_K"));
    vlpush(z4vars, Ind("bssn_Gx"));
    vlpush(z4vars, Ind("z4_Theta"));
    vlpush(z4vars, Ind("alpha"));
    vlpush(z4vars, Ind("betax"));
    vlpush(z4vars, Ind("betadotx"));
    
    /* register evolution routine */
    set_rhs_func_register("geometry",z4_evolve);
    set_rhs_func_register("geometry_adm",z4_to_adm_vl);
    set_rhs_func_register("geometry_source",z4_rhs_add_source);
    set_rhs_func_register("geometry_bound", set_boundary);
    
    /* register special z4 boundary condition */
#if 1
    if ( Getv("grid","shells") )
      // CPBC2 are used on shells if Getv("boundary", "background")
      boundary_register(z4_boundary_shell);
    else {
      // for boxes:
      if ( Getv("boundary", "background") && Getv("boundary", "radcentered") ) {
	// radiative on boxes - centered stencils
	// Nnote there're 2 options: T1 and T7 
	boundary_register(z4_boundary_box_rad);
	printf("  using radiative BCs with centered stencils on box boundary\n");
      } else if ( Getv("boundary", "background") ) { 
	// CPBC2 on boxes - centered stencils
	boundary_register(z4_boundary_box);
	printf("  using CPBC2 BCs with centered stencils on box boundary\n");
      }
    }
#else
    // old stuff
    boundary_register(z4_boundary);
#endif
    
  }

  z4vars->level = level;
  enablevarlist(z4vars);
  if (0) prvarlist(z4vars);
  
  // for debug reasons in movepunc_N
  enablevar(level,Ind("z4_Zx"));


  /* it is not quite clear yet where to set initial gauge */
  if (Getv("z4_initial_lapse", "one")) {
    double *alpha = level->v[Ind("alpha")];

    forallpoints(level, i) 
      alpha[i] = 1;

    printf("  set lapse to one\n");
  } else 
    printf("  left lapse unchanged\n");

  if (Getv("z4_initial_shift", "zero")) {
    double *betax = level->v[Ind("betax")];
    double *betay = level->v[Ind("betay")];
    double *betaz = level->v[Ind("betaz")];
    double *betadotx = level->v[Ind("betadotx")];
    double *betadoty = level->v[Ind("betadoty")];
    double *betadotz = level->v[Ind("betadotz")]; 
    
    forallpoints(level, i) {
      betax[i] = betay[i] = betaz[i] = 0;
      
      betadotx[i] = 0;
      betadoty[i] = 0;
      betadotz[i] = 0;
    }
    printf("  set shift to zero\n");
  } else
    printf("  left shift unchanged\n");

  /* set K identically to zero only if we are doing maximal slicing */
  if (Getv("z4_forceKzero", "yes")) {
    if (!GetvLax("Gauge", "maximal"))
      errorexit("\"z4_forceKzero = yes\" only makes sense in conjunction "
                "with maximal slicing");
  }

  /* translate initial data in ADM variables to z4 variables */
  adm_to_z4(level);
  set_boundary_symmetry(level, z4vars);

  return 0;
}

/* bssn.c */
/* Bernd Bruegmann 6/2002 */
/* Wolfgang Tichy  4/2004 */

#include "bam.h"
#include "bssn.h"

/* infos about the evolve varlist*/
int METRIC_bssn_INDX_VAR;




/* evolve in the interior and 
   for those boundary points set by special evolution
*/ 
void bssn_evolve(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  timer_start(ucur->level, "bssn_rhs");

  bssn_rhs_movpunc_N(unew, upre, c, ucur);
  
  timer_stop(ucur->level, "bssn_rhs");
}

void bssn_rhs_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  set_boundary(unew, upre, c, ucur);

  /* set physical boundary for eta */
  if (!Getv("bssn_use_eta", "no"))
    setbound_eta(ucur->level, Ind("bssn_eta"));
}







/* set boundary of Gammas
   - Cactus uses third order extrapolation
   - so far have only linear (second order) extrapolation
   - could also try one sided derivatives
*/
void bssn_setbound_G(tL *level, int igi, int iG)
{
  int i;
  
  for (i = 0; i < 3; i++)
    set_boundary_extrapolate(level, iG+i);
}

/* set physical boundary of eta */
void setbound_eta(tL *level, int ieta)
{
    if (0) printf("\nset physical boundary for eta\n");
    set_boundary_extrapolate(level, ieta);
}

/* set AMR boundary for eta; needed after ADMtoBSSN */
void synchAMR_eta(tL *level)
{
    tVarList* vl;
    vl = vlalloc(level);
    vlpush(vl, Ind("bssn_eta"));

    if (!parentaligned(level)) return;

    if (0) printf("\nset AMR  boundary for eta\n");
    set_boundary_refinement_cf(level->grid, level->l-1, level->l, vl);

    vlfree(vl);
}






/* compute trA and detg for BSSN variables
   enforce trA = 0 and detg = 1 if requested:
     gtilde_ij = detgtilde^(-1/3) gtilde_ij
     Atilde_ij = Atilde_ij - 1/3 gtilde_ij gtilde^kl Atilde_kl
*/
int bssn_algcon_wrapper(tL *level) 
{
  if (0) printf("BSSN: enforce algebraic constraints on l%d\n", level->l);
  
  tVarList *vl = vlalloc(level);

  vlpush(vl, Ind("bssn_gxx"));
  vlpush(vl, Ind("bssn_Axx"));
  vlpush(vl, Ind("bssn_trA"));
  vlpush(vl, Ind("bssn_detg"));

  /* do we actually want to store the result? */
  if (Getv("bssn_register_algcon", "store"))
    enablevarlist(vl);

  /* compute and enforce */
  bssn_algcon(vl);

  /* done */
  vlfree(vl);
  return 0;
}






/* derive ADM variables from BSSN variables
   g_ij = exp(4 phi) gtilde_ij
   A_ij = exp(4 phi) Atilde_ij
   K_ij = A_ij + 1/3 K g_ij
*/
int bssn_to_adm(tL *level) 
{
    tVarList *vl = vlalloc(level);
    vlpush(vl, Ind("bssn_gxx"));
    vlpush(vl, Ind("bssn_chi"));
    vlpush(vl, Ind("bssn_Axx"));
    vlpush(vl, Ind("bssn_K"));
    bssn_to_adm_vl(vl,vl,0,vl);
    vlfree(vl);
    
    return 1;
}

void bssn_to_adm_vl(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  if (0) printf("  deriving ADM variables from BSSN variables\n");
  tL* level = ucur->level;
  
  int ig      = Ind("adm_gxx");
  int iK      = Ind("adm_Kxx");
  int igtilde = Indvl("bssn_gxx",ucur);
  int iAtilde = Indvl("bssn_Axx",ucur);
  
  double *trK = level->v[Indvl("bssn_K",ucur)];
  double *chi = level->v[Indvl("bssn_chi",ucur)];
  double *g, *K, *gtilde, *Atilde;
  double psi, psi4;
  int n;
  
  double chipsipower = Getd("bssn_chi_psipower");

  forallpoints_ijk(level) {

    psi  = pow(chi[ijk],1./chipsipower);
    psi4 = pow(psi,4.);
    /*
    if (usechi)
      phi[ijk] = log(fabs(chi[ijk])) / chipsipower;
    
    fg = exp(4 * phi[ijk]);
    fK = fg;
    */
    
    
    for (n = 0; n < 6; n++) {
      g = level->v[ig+n];
      K = level->v[iK+n];
      gtilde = level->v[igtilde+n];
      Atilde = level->v[iAtilde+n];
      
      g[ijk] = psi4 * gtilde[ijk];
      K[ijk] = psi4 * Atilde[ijk] +  1./3. * trK[ijk] * g[ijk] ;
    }
    
  } endfor_ijk;
  
  if (0) {
    double dg,dgt;
    forallpoints_ijk(level) {
        double dg  = detg(level->v[ig+0][ijk],level->v[ig+1][ijk],level->v[ig+2][ijk],
                          level->v[ig+3][ijk],level->v[ig+4][ijk],level->v[ig+5][ijk]);
        double dgt = detg(level->v[igtilde+0][ijk],level->v[igtilde+1][ijk],level->v[igtilde+2][ijk],
                          level->v[igtilde+3][ijk],level->v[igtilde+4][ijk],level->v[igtilde+5][ijk]);
        if (dg<1e-10) { 
            printf("Problem with negativ detg in BSSNtoADM\n");
            printf("  %d pts away from boundary\n",boundaryaway(6));
            printf("  %e %e %e    l=%d\n",Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk], level->l);
            printf("  psi^4  = %e,    detg = %e,     detgtilde = %e\n", psi4,dg,dgt);
            printf("  g      = %e %e %e %e %e %e\n", level->v[ig+0][ijk],level->v[ig+1][ijk],level->v[ig+2][ijk],
                    level->v[ig+3][ijk],level->v[ig+4][ijk],level->v[ig+5][ijk]); 
            printf("  gtilde = %e %e %e %e %e %e\n", level->v[igtilde+0][ijk],level->v[igtilde+1][ijk],level->v[igtilde+2][ijk],
                    level->v[igtilde+3][ijk],level->v[igtilde+4][ijk],level->v[igtilde+5][ijk]);
            
            if (ijkinsidefinerlevel(box,ijk)==0) printf("    point is NOT inside finer box NOR in some symmetry area\n"); 
            else printf("    point is inside finer box/ in symmetry \n"); \
        }
    } endfor_ijk;
  }
}






/* derive BSSN variables from ADM variables
   p  = detgbar^(-1/3) 
   p0 = psi^(-4)

   gtilde_ij = p gbar_ij
   Ktilde_ij = p p0 K_ij

   phi = - log(p) / 4
   K   = gtildeinv^ij Ktilde_ij
   Atilde_ij = Ktilde_ij - gtilde_ij K / 3

   G^i = - del_j gtildeinv^ji
*/
void adm_to_bssn(tL *level) 
{
  printf("  deriving BSSN variables from ADM variables\n");

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
  vlpush(vl, Ind("bssn_eta"));
  
  /* compute */
  bssn_init(vl);
  
  /* free varlist */
  vlfree(vl);
  
  /* set boundary for Gammas */
  bssn_setbound_G(level, Ind("bssn_ginvtmpxx"), Ind("bssn_Gx"));
  disablevar(level, Ind("bssn_ginvtmpxx"));
  bampi_synchronize(level, Ind("bssn_Gx"));
  
  /* set boundary for eta */
  if (Getv("bssn_use_eta", "yes")) {
    setbound_eta(level, Ind("bssn_eta"));
    bampi_synchronize(level, Ind("bssn_eta"));
  }
}






/* initialize bssn after initial data has been computed in POST_INITIALDATA
*/
int bssn_startup(tL *level)
{
  int i;
  double vgauge;

  printf("Initializing bssn: \n");

  /* get ptr to evolution variables and add bssnvars
     DO NOT free this list, it is stored globally */
  tVarList* bssnvars = get_evolve_vlregister(level);
  
  if (level->l==0) {
    
    /* like in MATTER, we better remember where the bssn varaiables are 
       inside the evolution varlist, do this only once */
    METRIC_bssn_INDX_VAR = bssnvars->n;
    
    /* set gauge speed for lapse and related quantities */
    if (Getv("bssn_lapse", "1+log"))
      vgauge = sqrt(Getd("bssn_lapseharmonicf"));
    else if (Getv("bssn_lapse", "1+log2"))
      vgauge = sqrt(4.0/3.0);
    else
      vgauge = 1;
  
    /* set boundary information for bssn evolution: 
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
    VarNameSetBoundaryInfo("bssn_chi", 1, 1, vgauge);
    VarNameSetBoundaryInfo("alpha",    1, 1, vgauge);
    VarNameSetBoundaryInfo("betax",    0, 1, 1.0); 
    VarNameSetBoundaryInfo("betay",    0, 1, 1.0);
    VarNameSetBoundaryInfo("betaz",    0, 1, 1.0);
    VarNameSetBoundaryInfo("betadotx", 0, 1, 1.0);
    VarNameSetBoundaryInfo("betadoty", 0, 1, 1.0);
    VarNameSetBoundaryInfo("betadotz", 0, 1, 1.0);
    
    /* create a variable list for bssn evolutions 
       note that we include lapse and shift directly */
    vlpush(bssnvars, Ind("bssn_gxx"));
    vlpush(bssnvars, Ind("bssn_chi"));
    vlpush(bssnvars, Ind("bssn_Axx"));
    vlpush(bssnvars, Ind("bssn_K"));
    vlpush(bssnvars, Ind("bssn_Gx"));
    vlpush(bssnvars, Ind("alpha"));
    vlpush(bssnvars, Ind("betax"));
    vlpush(bssnvars, Ind("betadotx"));
    
    /* register evolution routine */
    set_rhs_func_register("geometry",bssn_evolve);
    set_rhs_func_register("geometry_adm",bssn_to_adm_vl);
    set_rhs_func_register("geometry_source",bssn_rhs_add_source);
    set_rhs_func_register("geometry_bound", bssn_rhs_boundary);
    
    /* register special bssn boundary condition */
    boundary_register(bssn_boundary);
  }
  
  /* enable all evolution vars */
  bssnvars->level = level;
  enablevarlist(bssnvars);
  if (0) prvarlist(bssnvars);
  
  /* enable eta if needed */
  if (Getv("bssn_use_eta", "yes"))
    enablevar(level, Ind("bssn_eta"));
  
  
  
  /* it is not quite clear yet where to set initial gauge */
  if (Getv("bssn_initial_lapse", "one")) {
    double *alpha = level->v[Ind("alpha")];

    forallpoints_ijk(level) { 
      alpha[ijk] = 1;
    } endfor_ijk;

    printf("  set lapse to one\n");

  /* add a gaugewave */
  } else if (Getv("bssn_initial_lapse", "onepluswave")) {
    double *alpha = level->v[Ind("alpha")];
    double *x = level->v[Ind("x")];
    double *y = level->v[Ind("y")];
    double *z = level->v[Ind("z")];

    forallpoints_ijk(level) {
      alpha[ijk] = 1 +
          0.001*exp(-5*(x[ijk]*x[ijk]+y[ijk]*y[ijk]+z[ijk]*z[ijk]));
    } endfor_ijk;
    printf("  set lapse to one + gauge wave\n");
  
  /* precollapsed lapse after we know the initial data 
     might be useful for superposition or puncture_ps-data etc. */
  } else if (Getv("bssn_initial_lapse", "oochi05")) {
    double *alpha = level->v[Ind("alpha")];
    double *bssn_chi = level->v[Ind("adm_gxx")];
    double *x = level->v[Ind("x")];
    double *y = level->v[Ind("y")];
    double *z = level->v[Ind("z")];
 
    forallpoints_ijk(level) { 
      alpha[ijk] =1.0/(sqrt(bssn_chi[ijk]));
    } endfor_ijk;
    printf("  set lapse to psi^-2 \n");    

   /* do nothing */
   }else
    printf("  left lapse unchanged\n");



  /* set initial shift */
  if (Getv("bssn_initial_shift", "zero")) {
    double *betax = level->v[Ind("betax")];
    double *betay = level->v[Ind("betay")];
    double *betaz = level->v[Ind("betaz")];
    double *betadotx = level->v[Ind("betadotx")];
    double *betadoty = level->v[Ind("betadoty")];
    double *betadotz = level->v[Ind("betadotz")];

     forallpoints_ijk(level) {
      betax[ijk] = betay[ijk] = betaz[ijk] = 0;
      betadotx[ijk] = betadoty[ijk] = betadotz[ijk] = 0;
    } endfor_ijk;
    printf("  set shift to zero\n");
  } else
    printf("  left shift unchanged\n");


  /* translate old bam11 par chi_div_floor, so we can use older parfiles */
  if(level->l==0 && GetsLax("chi_div_floor")!=0)
  {
    char *chi_div_floor = Gets("chi_div_floor");
    printf("  WARNING: bssn does not have parameter: chi_div_floor = %s\n",
           chi_div_floor);
    Sets("bssn_chi_div_floor", chi_div_floor);
    printf("    ==> Setting:  bssn_chi_div_floor = %s\n",
           Gets("bssn_chi_div_floor"));
  }

  /* translate initial data in ADM variables to BSSN variables */
  adm_to_bssn(level);
  set_boundary_symmetry(level, bssnvars);

  /* synchronize eta over AMR levels */
  if (Getv("bssn_use_eta", "yes"))
    synchAMR_eta(level);

  return 0;
}





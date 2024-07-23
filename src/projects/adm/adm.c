/* adm.c */
/* Bernd Bruegmann 11/02 */

#include "bam.h"
#include "adm.h"

#define PR 0





/* compute ADM constraints at ANALYZE time */
int computeadmconstraints(tL *level)
{
  int timer = timer_start(0, "admconstraints"); 
  
  tVarList *vl = vlalloc(level);
  vlpush(vl, Ind("ham"));
  vlpush(vl, Ind("momx"));
  if (Getv("adm_normalizedConstraints", "yes")) {
    vlpush(vl, Ind("normham"));
    vlpush(vl, Ind("normmomx"));
  }
  enablevarlist(vl);
  
  /* if time for output */
  if (timeforoutput(level, vl)) {

    if (PR) printf("computing adm constraints on l%d\n", level->l);
    
    // these should be initialized if used
    tVarList *wl = vlalloc(level);
    vlpush(wl, Ind("adm_gxx"));
    vlpush(wl, Ind("adm_Kxx"));
    vlpush(wl, Ind("adm_psi"));
    vlpush(wl, Ind("adm_dpsiopsix"));
    vlpush(wl, Ind("adm_ddpsiopsixx"));
    if (Getv("physics","matter")) {
      vlpush(wl, Ind("adm_rho"));
      vlpush(wl, Ind("adm_Sx"));
    }

    /* compute the adm constraints */
    adm_constraints_N(vl,wl);
    
    vlfree(wl);

    /* set boundaries */
    set_boundary_symmetry(level, vl);
    set_boundary_periodic(level, vl);

    /* synchronize */
    bampi_vlsynchronize(vl);

    /* restrict from l+1 to l and prolong from l to l+1 */
    if (level->l < level->grid->lmax)
      set_boundary_refinement_cf(level->grid, level->l, level->l+1, vl);

    /* Check for NANs and INFs */
    if (Getv("adm_NANcheck", "yes")) {
      if (CheckIfFinite(level, "ham") != 0) {
	printf("Level %d, time %f, iteration %d:\n",
	       level->l, level->time, level->iteration);
	errorexit("CheckIfFinite: NAN/INF error: ham is not finite!");
      }
    }
    
  } 

  vlfree(vl);
  timer_stop(0, "admconstraints"); 
  return 0;
}





/* init conformal transformation by enableling storage */
int adm_init_conftrans(tL *level)
{
  tVarList *vl;
  
  vl = vlalloc(level);
  vlpush(vl, Ind("adm_gxx"));
  vlpush(vl, Ind("adm_psi"));
  vlpush(vl, Ind("adm_dpsiopsix"));
  vlpush(vl, Ind("adm_ddpsiopsixx"));
  vlenable(vl);
  vlfree(vl);

  return 0;
}

/* undo conformal transformation */
int adm_undo_conftrans(tL *level)
{
  double *psi = Ptr(level, "adm_psi");
  double psi4;
  int ig      = Ind("adm_gxx");
  int idpop   = Ind("adm_dpsiopsix");
  int iddpop  = Ind("adm_ddpsiopsixx");
  int i, n;

  printf("adm: undoing conformal transformation\n");

  forallpoints(level, i) {
    psi4 = pow(psi[i], 4.0);
    psi[i] = 1;

    for (n = 0; n < 6; n++) 
      level->v[ig+n][i] *= psi4;

    if (level->v[idpop]) {
      for (n = 0; n < 3; n++) 
        level->v[idpop+n][i] = 0;
      for (n = 0; n < 6; n++) 
        level->v[iddpop+n][i] = 0;
    }
  }
  
  tVarList *vl;
  
  vl = vlalloc(level);
  vlpush(vl, Ind("adm_psi"));
  vlpush(vl, Ind("adm_dpsiopsix"));
  vlpush(vl, Ind("adm_ddpsiopsixx"));
  vldisable(vl);
  vlfree(vl);

  return 0;
}

/* compute and save initial value of the trace of K */
int set_K_initial(tL *level)
{
  double *gxx = Ptr(level, "adm_gxx");
  double *gxy = Ptr(level, "adm_gxy");
  double *gxz = Ptr(level, "adm_gxz");
  double *gyy = Ptr(level, "adm_gyy");
  double *gyz = Ptr(level, "adm_gyz");
  double *gzz = Ptr(level, "adm_gzz");

  double *Kxx = Ptr(level, "adm_Kxx");
  double *Kxy = Ptr(level, "adm_Kxy");
  double *Kxz = Ptr(level, "adm_Kxz");
  double *Kyy = Ptr(level, "adm_Kyy");
  double *Kyz = Ptr(level, "adm_Kyz");
  double *Kzz = Ptr(level, "adm_Kzz");

  double *psi = Ptr(level, "adm_psi");
  double *K_initial = PtrEnable(level, "adm_K0");

  double ixx, ixy, ixz, iyy, iyz, izz;
  int i;

  /* K_initial is zero upon enabling, return if that is what we want */
  if (!Getv("physics", "KerrSchild")) 
    return 0;

  printf("adm: set K_initial\n");
  
  forallpoints(level, i) {
    invg(gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i],
         &ixx, &ixy, &ixz, &iyy, &iyz, &izz);
    
    K_initial[i] = pow(psi[i],-4.0) * 
        (ixx*Kxx[i] + iyy*Kyy[i] + izz*Kzz[i] 
            + 2*(ixy*Kxy[i] + ixz*Kxz[i] + iyz*Kyz[i]));
    
  }

  
  
  
  return 0;
}


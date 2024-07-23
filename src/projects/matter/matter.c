/* matter.c 
   sbernuz */

#include "bam.h"
#include "matter.h"


#define PR 0

#define CONS 0

/* infos about all specific function like grhd ...*/
tMATTER MATTER;







/* ************************************ 
   startup */


tVarList* give_consVars_vl(tL *level)
{
  
  tVarList *consVars = vlalloc(level);
  
  // conservatives (evolved)
  int v;
  for(v=0; v<MATTER.NVq; v++) 
    vlpushone(consVars, Ind(MATTER.q_names_list[v]));
  
  return consVars;
}

tVarList* give_matterADMVars_vl(tL *level)
{
  
  tVarList *matterADMVars = vlalloc(level);
  
  // ADM  
  int v;
  for(v=0; v<MATTER.NVadm; v++) 
    vlpushone(matterADMVars, Ind(MATTER.a_names_list[v]));
  
  return matterADMVars;
}

tVarList* give_camrVars_vl(tL *level)
{ 
  tVarList *camrVars = vlalloc(level);
  
  // corrections for conservatives corresponding to the amr 
  int v;
  for(v=0; v<MATTER.NVq; v++) 
    vlpushone(camrVars, Ind(MATTER.c_names_list[v]));
  
  return camrVars;
}

tVarList* give_primVars_vl(tL *level)
{
  tVarList *primVars = vlalloc(level);
  
  // primitives
  int v;
  for(v=0; v<MATTER.NVw; v++) 
    vlpushone(primVars, Ind(MATTER.w_names_list[v]));
  
  return primVars;
}

tVarList* give_otherVars_vl(tL *level)
{
  tVarList *otherVars = vlalloc(level);

  //other vars
  int v;
  for (v=0; v<MATTER.NVo; v++)
    vlpushone(otherVars, Ind(MATTER.o_names_list[v]));

  return otherVars;
}

/* HG */

tVarList* give_nlsVars_vl(tL *level)
{
  tVarList *nlsVars = vlalloc(level);

  // nls
  int v;
  for(v=0; v<MATTER.NVnls;v++)
    vlpushone(nlsVars, Ind(MATTER.nls_names_list[v]));

  return nlsVars;
}

// NLS optical depth varList
tVarList* give_chiVars_vl(tL *level)
{
  tVarList *chiVars = vlalloc(level);
  int v;
  for(v=0; v<MATTER.NVchi;v++)
	vlpushone(chiVars, Ind(MATTER.chi_names_list[v]));

  return chiVars;

}

int matter_startup(tL *level)
{
  
  printf("Initializing matter:\n");
  printf("  provide matter fkt pointers\n");
  
  // set extra register for the matter part
  // this is ONLY the evolution and ADM part, adding ADM variables to geometry 
  // variables has to be handled by the geometry project!!!
  set_rhs_func_register("matter", matter_rhs);
  set_rhs_func_register("matter_set_null", matter_rhs_set_null);
  set_rhs_func_register("matter_c2p", matter_c2p);

  // other indexes often used
  MATTER.INDX_VAR_gxx  = Ind("adm_gxx");
  MATTER.INDX_VAR_detg = Ind("adm_detg");
  MATTER.INDX_VAR_adm  = Ind("adm_rho");
  MATTER.INDX_VAR_mask = Ind("matter_mask");
  
  // set varlistnames of ADM vars
  MATTER.NVadm   = 11; 
  MATTER.a_names_list = (char**) malloc( MATTER.NVadm * sizeof(char*) );
    
  MATTER.a_names_list[0] = "adm_rho";
  MATTER.a_names_list[1] = "adm_Sx";
  MATTER.a_names_list[2] = "adm_Sy";
  MATTER.a_names_list[3] = "adm_Sz";
  MATTER.a_names_list[4] = "adm_SSxx";
  MATTER.a_names_list[5] = "adm_SSxy";
  MATTER.a_names_list[6] = "adm_SSxz";
  MATTER.a_names_list[7] = "adm_SSyy";
  MATTER.a_names_list[8] = "adm_SSyz";
  MATTER.a_names_list[9] = "adm_SSzz";
  MATTER.a_names_list[10]= "adm_ST";

  MATTER.USENLS = 0;
  if(GetvLax("grhd_use_nls", "yes"))
    MATTER.USENLS = 1;
  
}

int matter_init(tL *level)
{

  printf("Initializing matter:\n");
  printf("  initialize matter varlists and set vars\n");
  int v,l;

  // add these variables to evolved variables
  tVarList *evolve_consVars = give_consVars_vl(level);
  enablevarlist(evolve_consVars);

  // use an extra VarList for ADM vars 
  tVarList *matterADMVars = give_matterADMVars_vl(level);
  enablevarlist(matterADMVars);

  // use an extra VarList for primitives 
  tVarList *primVars = give_primVars_vl(level);
  enablevarlist(primVars);

  //use an extra VarList for other variables
  tVarList *otherVars = give_otherVars_vl(level);
  enablevarlist(otherVars);

  // use an extra VarList for nls
  tVarList *nlsVars;
  tVarList *chiVarsc;
  tVarList *chiVarsf;

  if(MATTER.USENLS) {
    printf(" matter NLS on\n");
    nlsVars = give_nlsVars_vl(level);
    enablevarlist(nlsVars);
    chiVarsf = give_chiVars_vl(level);
    enablevarlist(chiVarsf);
    if(level->l > level->grid->lmin) chiVarsc = give_chiVars_vl((level->grid->level[level->l-1]));
  }

  // (if required) enable extra vars
  // atm mask, ...
  if (MATTER.USEMASK) {
    printf("  matter mask on\n");
    enablevar(level, Ind("matter_mask")); 
  }

   tVarList *camrVars;

  if (MATTER.USECAMR){  

    CAMR.correct_varlist = correct_varlist;
    CAMR.mask = c_amr_mask;

    printf("  add mask for conservative-amr \n");
    enablevar(level, Ind("camr_mask_A"));       
    enablevar(level, Ind("camr_mask_B"));       

    c_amr_mask(level);
    camrVars = give_camrVars_vl(level);
    enablevarlist(camrVars);

   }

  // find out location of alpha beta and consVars inside evolved VarList for faster/safe access
  // do this only once (in 0-level ... better during the firstcall)
  static int firstcall = 1;

  if (firstcall) {
    firstcall = 0;
    
    // take the pointer of evolved variable stack and 
    // add conservative variables to evolved variables ... only once
    tVarList *evolve_vars = get_evolve_vlregister(level); 
    vlpushvl(evolve_vars, evolve_consVars);
    if (MATTER.USECAMR) vlpushvl(evolve_vars, camrVars);
    if (MATTER.USENLS)  vlpushvl(evolve_vars, chiVarsf);

    if (!(MATTER.q_names_list) || !(MATTER.w_names_list))
      errorexit("No q or w varlist defined. Do you use a specifig hydro package?");

    for (MATTER.INDX_VAR_q  =0; MATTER.INDX_VAR_q  <evolve_vars->n; MATTER.INDX_VAR_q++  ) 
      if (Ind(MATTER.q_names_list[0])==evolve_vars->index[MATTER.INDX_VAR_q]) break;

    for (MATTER.INDX_VAR_alp=0; MATTER.INDX_VAR_alp<evolve_vars->n; MATTER.INDX_VAR_alp++) 
      if (Ind("alpha")==evolve_vars->index[MATTER.INDX_VAR_alp]) break;

    for (MATTER.INDX_VAR_bet=0; MATTER.INDX_VAR_bet<evolve_vars->n; MATTER.INDX_VAR_bet++) 
      if (Ind("betax")==evolve_vars->index[MATTER.INDX_VAR_bet]) break;

    if (MATTER.USECAMR){
   for (MATTER.INDX_VAR_camr=0; MATTER.INDX_VAR_camr<evolve_vars->n; MATTER.INDX_VAR_camr++) 
       if (Ind("camr_D")==evolve_vars->index[MATTER.INDX_VAR_camr]) break;}

    if (MATTER.USENLS){
   for(MATTER.INDX_VAR_chi=0; MATTER.INDX_VAR_chi<evolve_vars->n; MATTER.INDX_VAR_chi++)
       if (Ind("nls_chi_nue")==evolve_vars->index[MATTER.INDX_VAR_chi]) break;}

    if (MATTER.INDX_VAR_q==evolve_vars->n)    
       errorexit(" index of conserved quantities not found inside evolved variables \n");
    if (MATTER.INDX_VAR_alp==evolve_vars->n)  
       errorexit(" index of lapse not found inside evolved variables \n");
    if (MATTER.INDX_VAR_bet==evolve_vars->n)  
       errorexit(" index of beta not found inside evolved variables \n");
    if (MATTER.USECAMR) {if (MATTER.INDX_VAR_camr==evolve_vars->n)  
       errorexit(" index of camr not found inside evolved variables \n"); }
    if (MATTER.USENLS) {if (MATTER.INDX_VAR_chi==evolve_vars->n)
        errorexit( "index of chi not found inside evolved variables \n"); }

  }
 
  // copy initial data to matter variables
  // compute conservative and ADM vars    
  // (if required) atm is set here for the first time !	
  printf("  copy adm matter data to evolved matter vars\n");

  MATTER.init(level);
  if(MATTER.USENLS){
	 if(level-> l > level->grid->lmin) MATTER.chi_prolong(chiVarsc, chiVarsf, 1);
	 MATTER.chi_compute_init(level);
	 bampi_vlsynchronize(chiVarsf);
	 set_boundary_symmetry(level, chiVarsf);
	 MATTER.nls_init(level);
	 set_boundary_symmetry(level, nlsVars);
  }
  // free ptrs
  vlfree(matterADMVars);
  vlfree(primVars);
  vlfree(evolve_consVars);
  
  if (MATTER.USENLS) {
	vlfree(nlsVars);
	vlfree(chiVarsf);
	if(level->l > level->grid->lmin) vlfree(chiVarsc);
  }

  if (MATTER.USECAMR) vlfree(camrVars);
   
  return 0;
}









/* ************************************ 
   evolve */
void matter_c2p(tVarList *u_p, tL *level) {

  tVarList *primVars = give_primVars_vl(level);

  // c2p
  if (PR) printf("# ===> c2p previous\n");
  MATTER.c2p(u_p, primVars);

  // free ptrs
  vlfree(primVars);

}


void matter_rhs_prepare(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  tL *level = ucur->level;
  int nnodes = level->nnodes;

  int addlinear = (c != 0.0l);

  // set ptrs to rhs
  bampi_openmp_parallel_for_collapse2
  for(int v=0; v<MATTER.NVq; v++) {
    for (int i=0; i < nnodes; i++) {
      double *n = vldataptr(unew, MATTER.INDX_VAR_q + v);
      double *p = vldataptr(upre, MATTER.INDX_VAR_q + v);
      if (addlinear) 
        n[i] = p[i];
      else
        n[i] = 0.;
    }
  }

  if(MATTER.USENLS){
    bampi_openmp_parallel_for_collapse2
  	for(int v=0;v<MATTER.NVchi;v++) {
      for (int i=0; i < nnodes; i++) {
        double *nc = vldataptr(unew, MATTER.INDX_VAR_chi + v);
        nc[i] = 0.;
      }
  	}
  }
}

void matter_rhs_set_null(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{  
  tL *level = ucur->level;
  int nnodes = level->nnodes;

  // set ptrs to rhs
  bampi_openmp_parallel_for_collapse2
  for(int v=0; v<MATTER.NVq; v++) {
    for (int i=0; i < nnodes; i++) {
      double *r = vldataptr(ucur, MATTER.INDX_VAR_q + v);
      r[i] = 0.;
    }
  }
  

}

void matter_rhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  
  tL *level = ucur->level;
  tVarList *matterADMVars = give_matterADMVars_vl(level);
  tVarList *primVars      = give_primVars_vl(level);
  tVarList *otherVars	    = give_otherVars_vl(level);
  tVarList *camrVars; 
  tVarList *nlsVars;
  tVarList *chiVarsc;

  double *chie, *chia, *chix;

  int ex=0;

  if (MATTER.USECAMR) camrVars     = give_camrVars_vl(level);
  if (MATTER.USENLS) nlsVars       = give_nlsVars_vl(level);

  if (PR) printf("##############################################################\n");
  if (PR) printf("# MATTER:   %p \t %p \t %p \t %d\n",upre,ucur,unew, ucur->level->l);
  timer_start(level, "matter_rhs");

  // (if required) set excision
  if (MATTER.USEXCISION) {
    if (PR) printf("# ===> excision\n");
    MATTER.set_excision(upre,ucur);
  }
  
  // c2p 
  if (PR) printf("# ===> c2p\n");
  timer_start(level, "matter_rhs_c2p");
  MATTER.c2p(ucur, primVars);
  timer_stop(level, "matter_rhs_c2p");  
  
  // (if required) set atm (pre-rhs)
  if ((MATTER.USEATM) && (MATTER.USEATM_PRERHS)) {
    if (PR) printf("# ===> atm (pre-rhs)\n");
    MATTER.set_atm_prerhs(unew, upre, c, ucur, primVars);
  }

  // compute RHS 
  if (PR) printf("# ===> matter rhs\n");
  matter_rhs_prepare(unew, upre, c, ucur);

  if(MATTER.USENLS){
    if(!(EVO.stage)){
      if(PR) printf("# ===> matter NLS chi\n");
      MATTER.chi_compute(unew, upre, c, ucur);
    }
  }

  if (PR) printf("# ===> matter rhs sources\n");
  timer_start(level, "matter_rhs_source");
  MATTER.sources(unew, upre, c, ucur, primVars, otherVars, matterADMVars);
  if(MATTER.USENLS) set_boundary_symmetry(level, nlsVars);
  timer_stop(level, "matter_rhs_source");  

  if (PR) printf("# ===> matter rhs fluxes\n");
  timer_start(level, "matter_rhs_fluxes");
  if (Getv("physics", "grmhd")||Getv("physics", "grmhdY")||Getv("hrsc_flux", "HLLC"))
    matter_fluxes_grmhd(unew, upre, c, ucur, primVars, otherVars);
  else
    matter_fluxes(unew, upre, c, ucur, primVars, otherVars);
  timer_stop(level, "matter_rhs_fluxes");  
  
  // (if required) set atm (post-rhs)
  if ((MATTER.USEATM) && (MATTER.USEATM_POSTRHS)) {
    if (PR) printf("# ===> atm (post-rhs)\n");
    MATTER.set_atm_postrhs(unew, upre, c, ucur, primVars);
  }

  // sync ADM variables, because they are constructed by derivatives
  if (PR) printf("# ===> boundary + sync ADM\n");
  bampi_vlsynchronize(matterADMVars);
  set_boundary_symmetry(level, matterADMVars);
  
  timer_stop(level, "matter_rhs");
  if (PR) printf("# MATTER stop \n");
  if (PR) printf("##############################################################\n");

  // free ptrs
  vlfree(matterADMVars);
  vlfree(primVars);
  vlfree(otherVars);
  if (MATTER.USECAMR) vlfree(camrVars);
  if (MATTER.USENLS) vlfree(nlsVars);

}



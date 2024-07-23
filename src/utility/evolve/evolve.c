/* evolve.c */
/* Bernd Bruegmann 6/02 */
/* based on admc.c for Cactus, Bernd Bruegmann 9/98 */

/* provide the glue for the different evolution schemes */

#include "bam.h"
#include "evolve.h"

#define PR 0

/* variables and their copies for evolution 
   initialized to 0 for clarity since we rely on that
*/
tVarList *u_c = 0, *u_p = 0, *u_q = 0, *u_r = 0;
tVarList *u_pp = 0, *u_p2 = 0, *u_p3 = 0, *u_p4 = 0;
tVarList *u_aux = 0, *u_rhs = 0;

int comp;

/* store function pointer */
void (*rhs_geometry)        	  (tVarList*, tVarList*, double, tVarList*) = 0;
void (*rhs_matter_null)     	  (tVarList*, tVarList*, double, tVarList*) = 0;
void (*rhs_geometry_adm)    	  (tVarList*, tVarList*, double, tVarList*) = 0;
void (*rhs_matter)          	  (tVarList*, tVarList*, double, tVarList*) = 0;
void (*c2p_matter)              (tVarList*, tL*)                          = 0;
void (*rhs_geometry_source) 	  (tVarList*, tVarList*, double, tVarList*) = 0;
void (*rhs_boundary)        	  (tVarList*, tVarList*, double, tVarList*) = 0;
void (*rhs_radiation)	      	  (tVarList*, tVarList*, double, tVarList*) = 0;
void (*rhs_radiation_null)      (tVarList*, tVarList*, double, tVarList*) = 0;
int  (*radiation_check_evolve)  ()                                        = 0;
void (*rhs_radiation_expl)	    (tVarList*, tVarList*, double, tVarList*) = 0;
void (*rhs_radiation_impl)	    (tVarList*, tVarList*, double, tVarList*) = 0;


/* compute right-hand-side */
void evolve_rhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur, int comp) {

  // comp=0: compute rhs of spacetime and matter only
  // comp=1: compute rhs of radiation only
  // comp=2: compute rhs of spacetime, matter and radiation

  timer_start(ucur->level, "evolve_rhs");
  
  if (PR) printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  if (PR) printf("+ EVOLVE:   %p \t %p \t %p \t %d\n",upre,ucur,unew, ucur->level->l);
  
  if (comp==0 || comp==2) {

    /* geometry */
    if (PR) printf("+ ===> rhs_geometry rhs\n");
    rhs_geometry(unew, upre, c, ucur);

    /* matter */
    if (Getv("physics","matter")) {
      
      if(ucur->level->shells ||
        (ucur->level->l==0  && Getv("rhs_matter_null_ONlevel0","yes"))) {
        
        /* do not evolve matter inside spheres or on level0 */
        if (PR) printf("+ ===> rhs_matter_nulls\n");
        rhs_matter_null(unew, upre, c, ucur);
        
      } else {
      
        /* compute adm variables (g_ij,K_ij) */
        if (PR) printf("+ ===> rhs_geometry_adm\n");
        rhs_geometry_adm(unew, upre, c, ucur);
        
        /* compute matter rhs and adm source (rho,S,S^ij) */
        if (PR) printf("+ ===> rhs_matter\n");
        rhs_matter(unew, upre, c, ucur);

        /* compute matter-source (rho,S,S^ij) to geometry-evolution-equations */
        if (PR) printf("+ ===> rhs_geometry_source\n");
        rhs_geometry_source(unew, upre, c, ucur);

      }
    }

    if (Getv("physics","radiation") && comp==0) {
      rhs_radiation_null(unew, unew, 0, unew);
    }

  }

  if (comp==1 || comp==2) {
    if (comp==1) {
      if (PR) printf("+ ===> rhs_nulls\n");
      vlsetconstant(unew, 0); 
    }
    if (PR) printf("+ ===> rhs_radiation\n");
    rhs_radiation(unew, upre, 1, ucur);
  }

  if (comp==0 || comp==2) {
    /* boundaries (physical and symmetrie) ... if not explicitely set use 
      standard boundary function (which is normally called inside funkction) */
    if (PR) printf("+ ===> rhs_boundary\n");
    if (rhs_boundary)
      rhs_boundary(unew, upre, c, ucur);
    else
      set_boundary(unew, upre, c, ucur);
  }

  if (PR) printf("+ EVOLVE stop\n");
  if (PR) printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  
  timer_stop(ucur->level, "evolve_rhs");

}


void evolve_radiation(tVarList *ucur, double dt, tVarList *uk0, tVarList *uk1) {

  timer_start(ucur->level, "evolve_radiation");
  
  if (PR) printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  if (PR) printf("+ EVOLVE RADIATION:   %p \t %d\n",ucur,ucur->level->l);

  if (ucur->level->shells ||
     (ucur->level->l==0 && Getv("rhs_radiation_null_ONlevel0","yes"))) {

      /* do not evolve radiation inside spheres or on level0 */

  } else {
  
    vlsetconstant(uk0, 0); // initialize uk0 to 0
    vlsetconstant(uk1, 0); // initialize uk1 to 0

    rhs_radiation(uk0, u_c, 0, u_c);   // compute radiation rhs and store it in uk0 (does not add backreaction on fluid)
    rhs_matter_null(uk0, uk0, dt, uk0);
    bampi_vlsynchronize(uk0);
    vladd(uk1, 1, u_c, dt, uk0);       // evolve u_c for a dt with uk0 and store it in uk1

    rhs_radiation_expl(uk0, uk1, 0, uk1);  // compute explit part of rhs of uk1
    rhs_matter_null(uk0, uk0, dt, uk0);
    rhs_radiation_impl(uk0, u_c, 1, u_c);  // compute implicit part of rhs of u_c using uk0 as an explicit part (add backreaction on fluid)
    bampi_vlsynchronize(uk0);
    vladdto(u_c, dt, uk0);                 // update u_c for a dt unig uk0 

  }

  if (PR) printf("+ EVOLVE RADIATION stop\n");
  if (PR) printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
  
  timer_stop(ucur->level, "evolve_radiation");

}


/* store variables locally */
tVarList *get_evolve_vlregister(tL *level)
{
  if (!u_c)
    u_c = vlalloc(level);
  u_c->level = level;
  return u_c;
}

/* retrieve variables */
void evolve_vlretrieve(tVarList **vlu_c, tVarList **vlu_p, tVarList **vlu_pp)
{
  if (vlu_c)  *vlu_c  = u_c;
  if (vlu_p)  *vlu_p  = u_p;
  if (vlu_pp) *vlu_pp = u_pp;
}

void evolve_vlretrieve_cpq(tVarList **vlu_c, tVarList **vlu_p, tVarList **vlu_q)
{
  if (vlu_c)  *vlu_c = u_c;
  if (vlu_p)  *vlu_p = u_p;
  if (vlu_q)  *vlu_q = u_q;
}

void evolve_vlretrieve_pn(tVarList **vlu_p2, tVarList **vlu_p3, tVarList **vlu_p4)
{
  if (vlu_p2) *vlu_p2 = u_p2;
  if (vlu_p3) *vlu_p3 = u_p3;
  if (vlu_p4) *vlu_p4 = u_p4;
}


/* store function pointer locally */
void set_rhs_func_register(char *name, void (*f)(tVarList *, tVarList *, double, tVarList *))
{
  if (strcmp(name,"geometry")==0) {
    rhs_geometry = f;
  } else if (strcmp(name,"matter_set_null")==0) {
    rhs_matter_null = f;
  } else if (strcmp(name,"geometry_adm")==0) {
    rhs_geometry_adm = f;
  } else if (strcmp(name,"matter")==0) {
    rhs_matter = f;
  } else if (strcmp(name,"matter_c2p")==0) {
    c2p_matter = f;
  } else if (strcmp(name,"geometry_source")==0) {
    rhs_geometry_source = f;
  } else if (strcmp(name,"geometry_bound")==0) {
    rhs_boundary = f;
  } else if (strcmp(name,"radiation")==0) {
    rhs_radiation = f;
  } else if (strcmp(name,"radiation_impl")==0) {
    rhs_radiation_impl = f;  
  } else if (strcmp(name,"radiation_expl")==0) {
    rhs_radiation_expl = f;
  } else if (strcmp(name,"radiation_set_null")==0) {
    rhs_radiation_null = f;
  } else if (strcmp(name,"radiation_check_evolve")==0) {
    radiation_check_evolve = f;
  }
  else {
    errorexit("this function point does not exist");
  }
}



/* evolve wrapper for function skeleton */
int evolve(tL *level) 
{
  int timer = timer_start(level, "evolve");
  int timeinterpolation = Geti("order_timeinterpolation");
  int i, j;
  int comp;

  /* make sure we have an evolution routine and a list of variables */
  /* fixme: make evolve_ngs registrable */
  if (!rhs_geometry) { 
    return 0;
    // let's assume someone else provides evolution
    // errorexit("evolve: no evolution routine");
  }
  if (!u_c) 
    errorexit("evolve: no list of variables");

  /* store level for this call to evolve */
  u_c->level = level;

  /* add new variables before first time step 
     note that we call AddDuplicate for each level so that level->v is expanded
  */
  if (level->iteration == 0) {
    if (level->grid->lmin == level->l) 
      printf("Adding variables for evolution:\n");
    
    /* cartoon method needs points before first evolution step */
    if ((Getv("grid","1d")) || (Getv("grid","2d")))
      set_boundary_symmetry(level, u_c); 
    
    /* add additional variables */
    if (Getv("evolution_method", "euler")) {
      u_p = AddDuplicate(u_c, "_p");
    }
    if (Getv("evolution_method", "icn")) {
      u_p = AddDuplicate(u_c, "_p");
      u_q = AddDuplicate(u_c, "_q");
    }
    if (Getv("evolution_method", "rk")) {
      u_p = AddDuplicate(u_c, "_p");
    }
    if (Getv("evolution_method", "ngs")) {
      u_p = AddDuplicate(u_c, "_p");
      u_q = AddDuplicate(u_c, "_q");
    }

    /* amr needs storage for interpolation to half step */
    if (Geti("amr_lmax") > 0) {
      if (!u_p) AddDuplicate(u_c, "_p");
      u_pp = AddDuplicate(u_c, "_pp");
      if (timeinterpolation >= 4) u_p2 = AddDuplicate(u_c, "_p2");
      if (timeinterpolation >= 5) u_p3 = AddDuplicate(u_c, "_p3");
      if (timeinterpolation >= 6) u_p4 = AddDuplicate(u_c, "_p4");
    }
  }

  /* store current level in existing variable lists */
  if (u_p) u_p->level = level;
  if (u_q) u_q->level = level;
  if (u_r) u_r->level = level;
  if (u_pp) u_pp->level = level;
  if (u_p2) u_p2->level = level;
  if (u_p3) u_p3->level = level;
  if (u_p4) u_p4->level = level;
  if (u_aux) u_aux->level = level;

  /* turn on memory (if varlist is non null) */
  enablevarlist(u_p);
  enablevarlist(u_q);
  enablevarlist(u_r);
  enablevarlist(u_pp);
  enablevarlist(u_p2);
  enablevarlist(u_p3);
  enablevarlist(u_p4);
  enablevarlist(u_aux);


  /* start forming in-time interpolation to half step, e.g. 3rd order
     u(n+1/2) = (3*u(n+1) + 6*u(n) - u(n-1))/8;
  */
  if ( (level->l < level->grid->lmax)) {
    int o;

    /* reduce order if we don't have yet all the levels available */
    if (level->iteration + 2 > timeinterpolation)
      o = timeinterpolation;
    else
      o = level->iteration + 2;

    if (o == 2) {
      vladd(u_pp, 0.5, u_c, 0.0, u_c);
    }

    if (o == 3) {
      vladd(u_pp, 0.75, u_c, -0.125, u_p);
    }

    if (o == 4) {
      vladd(u_pp, 15.0/16, u_c, -5.0/16, u_p);
      vladd(u_pp,     1.0, u_pp, 1.0/16, u_p2);
    }

    if (o == 5) {
      vladd(u_pp, 140.0/128, u_c, -70.0/128, u_p);
      vladd(u_pp,       1.0, u_pp, 28.0/128, u_p2);
      vladd(u_pp,       1.0, u_pp, -5.0/128, u_p3);
    }

    if (o == 6) {
      vladd(u_pp, 315.0/256, u_c,  -210.0/256, u_p);
      vladd(u_pp,       1.0, u_pp,  126.0/256, u_p2);
      vladd(u_pp,       1.0, u_pp,  -45.0/256, u_p3);
      vladd(u_pp,       1.0, u_pp,    7.0/256, u_p4);
    }
    
    /* shift back time levels */
    if (o >= 6) vlcopy(u_p4, u_p3);
    if (o >= 5) vlcopy(u_p3, u_p2);
    if (o >= 4) vlcopy(u_p2, u_p);
  }


  /* perform one evolution time step
     choose between different method of line integrators
  */
  if (Getv("physics","radiation")) {
    if (Getv("rad_evolution_method","same_as_matter")) comp = 2;
    else                                               comp = 0;
  }
  else comp = 0;

  if (Getv("evolution_method", "icn"))
    evolve_icn(level, comp);
  else if (Getv("evolution_method", "euler"))
    evolve_euler(level,comp);
  else if (Getv("evolution_method", "rk"))
    evolve_rk(level,Gets("evolution_method_rk"),comp);
  else
    errorexit("evolution_method not known, use rk");

  if (comp==0) {
    if (Getv("physics","radiation") && (radiation_check_evolve())) {
      c2p_matter(u_p, level);
      rhs_geometry_adm(u_p, u_p, 0, u_p);
      if      (Getv("rad_evolution_method","imex")) evolve_rad_imex(level);
      else if (Getv("rad_evolution_method","rk"))  evolve_rk(level,Gets("rad_evolution_method_rk"),1);
      else errorexit ("rad_evolution_method not known");
    }
  }

  /* finish time interpolation */
  if ((level->l < level->grid->lmax)) {
    double c[10] = {0, 0, 1.0/2, 3./8, 5./16, 35./128, 63./256};

    if (level->iteration + 2 > timeinterpolation)
      vladd(u_pp, 1.0, u_pp, c[timeinterpolation], u_c);
    else
      vladd(u_pp, 1.0, u_pp, c[level->iteration + 2], u_c);
  }


  /* compute difference between new and old data, store in u_change */
  for (j = 0; j < u_c->n; j++)
  {
    if (Getv("evolve_compute_change", VarName(u_c->index[j])))
    {
      tVarList *var, *change;
      double *old = level->v[u_p->index[j]];
      double *new = level->v[u_c->index[j]];
      double *ch;

      /* make temporary varlist containing only the current variable */
      var = vlalloc(level);
      vlpush(var, u_c->index[j]);
      enablevarlist(var);

      /* duplicate var, to create var_change */
      change = AddDuplicateEnable(var, "_change");
      ch = level->v[change->index[0]];

      /* loop over all points and write into var_change */
      forallpoints(level, i)
      {
        ch[i] = new[i] - old[i];
      }
      /* free temp. var lists */
      vlfree(change);
      vlfree(var);
    }
  }

  /* THIS BLOCK IS OBSOLETE, SINCE CHANGES ARE NOW STORED IN u_change */
  /* compute difference between new and old data, store in u_p */
  if (Geti("amr_lmax") == 0) {
    for (j = 0; j < u_c->n; j++) {
      if (Getv("evolve_compute_change", VarName(u_c->index[j]))) {
	double *old = level->v[u_p->index[j]];
	double *new = level->v[u_c->index[j]];
	forallpoints(level, i) {
	  old[i] = new[i] - old[i];
        }
      }
    }
  }

  /* turn off memory (if varlist is non null) */
  if (Getv("evolve_persist", "no")) {
    if (Geti("amr_lmax") == 0) 
      disablevarlist(u_p);
    disablevarlist(u_q);
    disablevarlist(u_r);
    disablevarlist(u_aux);
  }

  /* turn on memory */
  //enablevarlist(evolve_no_memory);


  /* hack: keep some variables constant */
  if (0) {
    for (j = 0; j < u_c->n; j++) {
      i = u_c->index[j];
      if (!strstr(VarName(i), "alpha") && !strstr(VarName(i), "beta")) {
	double *old = level->v[u_p->index[j]];
	double *new = level->v[u_c->index[j]];
	forallpoints(level, i) {
	  new[i] = old[i];
        }
      }
    }
  }

  timer_stop(level, "evolve");
  return 0;
}




/* store right hand sides, called as analysis function */
int evolve_store_rhs(tL *level, int comp)
{
  if (!Getv("evolve_store_rhs", "yes")) 
    return 0;

  if (!Getv("evolution_method", "euler")) {
      printf("ATTENTION: I hope you know that a non euler\n");
      printf("           time-intergator gives stupid results!\n");
      printf("            At least this implementation.\n");
  }
  
  if (level->iteration == 0) {
    u_rhs = AddDuplicate(u_c, "_rhs");
  }

  u_p->level = u_rhs->level = level;
  enablevarlist(u_rhs);

  evolve_rhs(u_rhs, u_p, 0.0, u_p, comp);
  bampi_vlsynchronize(u_rhs);

  return 0;
}




/* determine globally the size of dissipation */
double get_dissipation_factor(tL *level)
{
  static int order_dissipation;
  static int amr_move_lcube;
  static int z4;
  static double dissipation_factor;
  static double dissipation_factor_level0;
  static double shells_factor = -1.;
  static double move_factor = -1.;
  
  static int firstcall = 1;
  if (firstcall) {
    firstcall = 0;
    
    order_dissipation   = Geti("order_dissipation");
    amr_move_lcube      = Geti("amr_move_lcube");
    z4                  = Getv("physics","Z4d");
    dissipation_factor  = Getd("dissipation_factor");
    dissipation_factor_level0 = Getd("dissipation_factor_level0");
    if (level->shells)
      shells_factor     = Getd("dissipation_factor_shells");
    move_factor         = Getd("dissipation_factor_move");
  }
  
  double dissfactor = dissipation_factor/pow(2.0, (double)order_dissipation);

  if (level->shells && shells_factor != -1.)
    return shells_factor/pow(2.0, (double)order_dissipation);
  
  else if (level->l > amr_move_lcube && move_factor != -1.)
    return move_factor/pow(2.0, (double)order_dissipation);
  
  else if (level->l > amr_move_lcube && !z4)
    return 0.2 * dissipation_factor/pow(2.0, (double)order_dissipation);
 
  else if ( (level->l == 0) &&  (dissipation_factor_level0 != -1.))
    return dissipation_factor_level0/pow(2.0, (double)order_dissipation);   

  else 
    return dissipation_factor/pow(2.0, (double)order_dissipation);
  
}



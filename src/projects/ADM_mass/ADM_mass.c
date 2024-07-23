/* ADM_mass.c */
/* Jose Gonzalez 9/04 */

#include "bam.h"
#include "ADM_mass.h"



/* compute ADM_mass
*/
int ADM_mass(tL *level)
{
  timer_start(0, "adm_mass");
  tVarList *vl;
  tG *grid = level->grid;
  int lmin = level->l;
  int lmax = grid->lmax;
  int l;
  double integral,integral_chi;
  double integralPx,integralPy,integralPz,integralJx,integralJy,integralJz;
  int ncircles = Geti("ADM_mass_ncircles");
  int npoints  = Geti("ADM_mass_npoints");
  int order = Geti("order_centered");
  char name_adm_mass[2001];
  char rname[10];
  
  static int firstcall = 1;
  char outdir[1000];
  sprintf(outdir,"%s/ADM_mass", Gets("outdir"));
  if (processor0 && firstcall) {
    system_mkdir(outdir);
    firstcall = 0;
  }
  
  FILE *fp_adm_mass;
  int i;
  int rank = bampi_rank();

  /* we wait until this level is the coarsest at its time and 
     then deal with the entire time-aligned stack
  */
 
  if (parentaligned(level)) return 0;

  /* variables we may want */
  vl = vlalloc(0);
  vlpush(vl, Ind("ADM_mass_integrand"));
  vlpush(vl, Ind("ADM_mass_integrand_chi"));
  vlpush(vl, Ind("ADM_mass_Pxint"));
  vlpush(vl, Ind("ADM_mass_Pyint"));
  vlpush(vl, Ind("ADM_mass_Pzint"));
  vlpush(vl, Ind("ADM_mass_Jxint"));
  vlpush(vl, Ind("ADM_mass_Jyint"));
  vlpush(vl, Ind("ADM_mass_Jzint"));

  /* for all levels from finest to coarse */
  for (l = lmax; l >= lmin; l--) {
    level = grid->level[l];

    /* make sure storage is enabled */
    vl->level = level;
    enablevarlist(vl);

    /* check whether these variables are wanted, and are wanted now 
       if not, then the coarsest level for output is one finer, l+1
       if lmin > lmax, then no output happens at all
    */

    /* frequency of output is that of scalar output (for backwards compatibility) */
    if (!timeforoutput_di_dt(level, Geti("0doutiter"), Getd("0douttime"))) {
      lmin = l+1;
      break;
    }

    /* ADM mass integrand */
    adm_mass_integrand(level, Ind("x"), Ind("adm_gxx"), Ind("ADM_mass_integrand"));

    /* chi version of ADM mass integrand */
    adm_mass_integrand_chi(level, Ind("x"), Ind("bssn_chi"), Ind("ADM_mass_integrand_chi"));

    /* P and J integrand */
    adm_pj_integrand(level, Ind("x"), Ind("adm_gxx"), Ind("adm_Kxx"), Ind("bssn_K"));

    /* set boundaries */
    set_boundary_symmetry(level, vl);
    
    /* synchronize */
    bampi_vlsynchronize(vl);
    
  }

  /* now we have the ADM mass on all the available time-aligned levels
     - set boundary from parent
     - set interior from children
  */
  for (l = lmax; l > lmin; l--)
    set_boundary_refinement_cf(grid, l-1, l, vl); 

  /* post process */
  for (l = lmax; l >= lmin; l--) {
    level = grid->level[l];
    if (0) printf("ADM mass level %d, time %6.3f\n", level->l, level->time); 

    /* compute and output the ADM mass */
    if (level->l >= Geti("ADM_mass_lmin") &&
	level->l <= Geti("ADM_mass_lmax")) {

      /* do computation only for certain radii */
      char *rstrings = Gets("ADM_mass_r");
      char *rstring;
      double r, rlist[100];
      int n = 0;

      /* make list of radii*/
      while (rstring = NextEntry(rstrings)) {
	r = atof(rstring);
        if (0) printf("rstring %s,  r = %f\n", rstring, r);
	if (0) printbbox(level, level->bbox, 0);
	if (centered_sphere_safely_inside_level(level, r)) {
	  if (0) printf("computing ADM mass for level %d, r %.1f\n", level->l, r);
	  rlist[n++] = r;
	}
      }

      /* we now have all the integrands we need, go ahead and integrate for each radius */
      for (i = 0; i < n; i++) {
	r = rlist[i];
	sprintf(rname, "r%d", i+1);

	/* integrate */
        integral     = integral_over_sphere(level, 0,0,0 ,ncircles, npoints, r, Ind("ADM_mass_integrand"), order);
        integral_chi = integral_over_sphere(level, 0,0,0, ncircles, npoints, r, Ind("ADM_mass_integrand_chi"), order);
        integralPx   = integral_over_sphere(level, 0,0,0, ncircles, npoints, r, Ind("ADM_mass_Pxint"), order);
        integralPy   = integral_over_sphere(level, 0,0,0, ncircles, npoints, r, Ind("ADM_mass_Pyint"), order);
        integralPz   = integral_over_sphere(level, 0,0,0, ncircles, npoints, r, Ind("ADM_mass_Pzint"), order);
        integralJx   = integral_over_sphere(level, 0,0,0, ncircles, npoints, r, Ind("ADM_mass_Jxint"), order);
        integralJy   = integral_over_sphere(level, 0,0,0, ncircles, npoints, r, Ind("ADM_mass_Jyint"), order);
        integralJz   = integral_over_sphere(level, 0,0,0, ncircles, npoints, r, Ind("ADM_mass_Jzint"), order);

        /* mth: I do not know why we do this after integration... this makes no sense at all */
	if (Getv("grid", "box") && !Getv("grid","shells")) {
	  if (Getv("grid", "bitant")){
	    /* z=0 reflection symmetry */
	    integralPz = 0;
	    if (0) printf("bitant\n");
	  } else if (Getv("grid", "quadrant")){
	    /* z=0 reflection and x=y=0 inversion symmetry */
	    integralPx = 0;
	    integralPy = 0;
	    integralPz = 0;
	    integralJx = 0;
	    integralJy = 0;
	    if (0) printf("quadrant\n");
	  } else if (Getv("grid", "qreflect")){
	    /* y=0, and z=0 reflection symmetry */
	    integralPy = 0;
	    integralPz = 0;
	    integralJx = 0;
	    integralJy = 0;
	    integralJz = 0;
	    if (0) printf("qreflect.\n");
	  } else if (Getv("grid", "rotant")) {
	    /* y=0, and z=0 reflection symmetry */
	    integralPx = 0;
	    integralPy = 0;
	    integralJx = 0;
	    integralJy = 0;
	    if (0) printf("rotant.\n");
	  }
	}
	
	/* process 0 does the writing */ 
	if (!rank) {
          // fixme: avoid overflow
	  sprintf(name_adm_mass,   "%s/ADM_mass_%s.l%d",  outdir, rname, level->l);
      
	  fp_adm_mass = fopen(name_adm_mass,   "ab");

	  if (!fp_adm_mass) errorexit("ADM mass: failed opening file");

	  /* write one line header before first output 
	     does that actually make sense?
	  */
	  if (level->iteration == 0) {
	    fprintf(fp_adm_mass, "\" ADM mass:   r = %14.6f \"\n",  r);
	    fprintf(fp_adm_mass, "\"%14s%14s%14s%14s%14s%14s%14s%14s%14s\"\n", 
                    "time    ","mass (full) ","mass (chi)  ",
                    "Px      ","Py      ","Pz      ",
		    "Jx      ","Jy      ","Jz      ");
	  } 
      
	  fprintf(fp_adm_mass, "%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n", 
		  level->time, integral,integral_chi,integralPx,integralPy,integralPz,
		  integralJx,integralJy,integralJz);

	  fclose(fp_adm_mass);
	}

      }
      
    }
    
    /* disable temporary storage */
    if (!Getv("ADM_mass_persist", "yes")) {
      vl->level = level; 
      disablevarlist(vl);
    }
  }

  /* finish */
  vlfree(vl);
  timer_stop(0, "adm_mass");
  return 0;
}












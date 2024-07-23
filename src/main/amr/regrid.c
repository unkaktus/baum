/* regrid.c */
/* Bernd Bruegman 12/99 */

#include "bam.h"
#include "amr.h"




/* given a refinement function, add and remove new grid points */

/* l:   this level will be changed or created
   l-1: the refinement function lives on the coarser level
   l+1: don't remove nodes that have a refinement
*/
void regrid(tG *g, int l, int secondcall) 
{ 
  tL *lc, *lnew, *lold, *levelpoint;
  int nnnew, nnc, nglobal;
  int i, j, k, t, d, di, dj, dk;
  int istart, iend, inew, inewstart, inewstartfine, jstart;
  double *xc, *yc, *zc, x, y, z;
  double *flagregrid;
  int pr = 1;

  /* top level may not be regridded */
  if (l <= g->lmin) return;
  
  /* number of levels is limited */
  if ((l >  Geti("amr_lmax"))   && (!(secondcall))) return;
  if ((l >  Geti("amr_lfinal")) && (secondcall))  return;


  if (pr) prdivider(0);
  if (pr) printf("Regrid l%d\n", l);


  /* find points in lc that will have children in new level 
     return if nothing to do, e.g. fixed mesh refinement
  */

  /* if level l already exists, this is an error for fmr */
  if (l <= g->lmax)
    errorexit("regrid: fmr tries to recreate existing level");

  /* create new level based on nesting function */
  lc = g->level[l-1];
  lnew = make_level_nested_boxes(lc, 0);


  /* set new size of time step */
  lnew->dt = lc->dt / 2;
  lnew->time = lc -> time; 
  lnew->iteration = lc -> iteration; 
  if (Getv("amr", "onedt")) {
    for (i = 0; i < lnew->l; i++) 
      g->level[i]->dt = lnew->dt;
  }
  else if (Getv("amr", "bo")) {
    double dxmax = Getd("amr_bo_dxmax");
    double bolmin = Getd("amr_bo_lmin");

    /* amr_bo_lmin = 0 is a possible value, but this is the 
     * same as absence of amr_bo-Parameters */
    /* use amr_bo_lmin to compute dxmax if specified and not
     * zero */
    if (bolmin){
        dxmax = (g->level[0]->dx)/pow(2,bolmin);
        printf("bolmin= %f, dxmax= %f \n",bolmin,dxmax);
    }

    //if (dxmax > 0 && !dless(lnew->dx, dxmax)) {  // dx >= dxmax
    if (dxmax > 0 && dless(dxmax/2,lnew->dx)) {    // dx > dxmax/2
      for (i = 0; i < lnew->l; i++) 
	g->level[i]->dt = lnew->dt;
      g->lbo = lnew->l;
    }
  }
  if (0) printf("lnew: l%d dt=%f dtpr=%f lbo=%d\n", 
		lnew->l, lnew->dt, g->level[lnew->l-1]->dt, g->lbo);

  if(secondcall) {
    Seti("amr_lmax", Geti("amr_lmax")+1);
  }

  /* info */ 
  if (pr) {
    printf("Created new level:\n");

    forallboxes(lnew) {
      printbbox(lnew, box->bbox, box->ibbox);
    } endforboxes;
   
    if (0) printf("dt = %f, dt/dx = %f\n", lnew->dt, lnew->dt/lnew->dx);
    prdivider(0);
    if (l == Geti("amr_lmax")) {
      prdivider(0);
      printgrid(g);
      prdivider(0);
    }
  }
}

/*create new levels during the evolution and fill them with data from the parent levels
  not tested in all details :) */
int add_level_while_evolving(tL *level) {

  if(level->l < Geti("amr_lmax")) return 0;
  if(!parentaligned(level)) return 0;
  if(level->l == (Geti("amr_lfinal"))) return 0;
  if((level->time) < Getd("amr_ladd_time")) return 0;
  
  int l = level->l; 
  tG* g=level->grid;  
  tVarList *unew, *unewp, *unewpp;
    	
  regrid(g, l+1, 1);  
  realloc_levelvariables(g->level[l+1], g->nvariables);

  evolve_vlretrieve(&unew, &unewp, &unewpp);
  unewpp = unew;  			      
  prolong_varlist(g, l, l+1, unewpp, unew,  -1);			      
  g->level[l+1]->iteration = 2*(level->iteration);
  printgrid(g);

  add_level_while_evolving(g->level[l+1]);  
  
  return 1;
}

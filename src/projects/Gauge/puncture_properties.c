/* puncture_properties.c */
 /* TD 07/13*/
 /* TD 11/15*/

/* this function simply computes the ADM integral on coordinate spheres 
   around the BHs or NSs in a hope to get some quasi-local meassure of the spin. 
   Momenta and Mass are computed as well but even more debatable...*/

 
#include "bam.h"
#include "Gauge.h"

static int firstcall_puncprop = 1;

int puncture_properties(tL *level)
{
  timer_start(0, "puncture_properties");
  tVarList *vl;
  tG *grid = level->grid;
  int lmin = grid->lmin;
  int lmax = grid->lmax;
  int l;
  double integral;
  int ncircles = Geti("puncture_properties_circles");
  int numpunc = Geti("puncture_properties_punc");
  int npoints  = Geti("puncture_properties_npoints");
  int order = Geti("order_centered");
  char name_psp[2001];
  char rname[10],jname[10];
  double  p[2][3];
  char outdir[1000];
  double integralM;
  double integralPx, integralPy, integralPz;
  double integralSx, integralSy, integralSz;
  tB *fbox;

  sprintf(outdir,"%s/puncture_properties", Gets("outdir"));
  if (processor0 && firstcall_puncprop == 1) {
    system_mkdir(outdir);
    firstcall_puncprop = 0;
    }
  
  FILE *fp_psp;

  int i,j,k;
  int rank = bampi_rank();
  double x0,y0,z0;

  /* skip if time or level is wrong */

  if (parentaligned(level)) 
      { timer_stop(0, "puncture_properties"); return 0; }
  if (!timeforoutput_di_dt(level, Geti("0doutiter"), Getd("0douttime"))) 
      { timer_stop(0, "puncture_properties"); return 0;}

  /* define puncture positions */
  for (j = 0; j < numpunc; j++)
    for (i = 0; i < 3; i++) {
        p[j][i] = level->grid->puncpos[j][i];
    } 

  /* enable the storage */
  vl = vlalloc(0);
  vlpush(vl, Ind("puncture_properties_Mint"));
  vlpush(vl, Ind("puncture_properties_Pxint"));
  vlpush(vl, Ind("puncture_properties_Pyint"));
  vlpush(vl, Ind("puncture_properties_Pzint"));
  vlpush(vl, Ind("puncture_properties_Sxint"));
  vlpush(vl, Ind("puncture_properties_Syint"));
  vlpush(vl, Ind("puncture_properties_Szint"));



  for (l = lmax; l >= lmin; l--) {
    level = grid->level[l];
    vl->level = level;
    enablevarlist(vl);
  }

  /* compute puncture properties */
 for (j = 0; j < numpunc; j++) {

        x0=p[j][0];
        y0=p[j][1];
        z0=p[j][2];

  for (l = grid->lmaxpunc[j]; l >= lmin; l--) {
    level = grid->level[l];

    /* note that shells are not included in puncture_properties_integrand*/
    if(level->shells) break; 

      /* computation only for certain radii */
      char *rstrings = Gets("puncture_properties_r");
      char *rstring;
      double r, rlist[100];
      int n = 0;

      fbox = box_containing_xyz(level, x0, y0, z0);

           /* list of radii*/
        while (rstring = NextEntry(rstrings)) {
            r = atof(rstring);
           /* sphere_inside_level is not working if we have more than one box per level, 
              this simple check looks if in one dimension the box is larger than the radius 
              we hope that the symmetry will take care of the other dimensions, but it is not checked*/
           //if (sphere_inside_level(level, x0, y0, z0 , r)) {
            if  ((fabs(fbox->bbox[0]-fbox->bbox[1])>r)||
                 (fabs(fbox->bbox[2]-fbox->bbox[3])>r)||
                 (fabs(fbox->bbox[4]-fbox->bbox[5])>r)) {
               rlist[n++] = r; 
           }
          }

    /*integrand */
    puncture_properties_integrand(level, Ind("x"), Ind("adm_gxx"), Ind("adm_Kxx"), Ind("bssn_K"),x0,y0,z0);

    /* set boundaries */
    set_boundary_symmetry(level, vl);

    /* synchronize */
    bampi_vlsynchronize(vl);

    set_boundary_refinement_cf(grid, l-1, l, vl);

        /* integrate for each radius and puncture */
     for (i = 0; i < n; i++) {  
       r = rlist[i];
       sprintf(rname, "r%d", i+1);
       sprintf(jname, "p%d", j+1); 

       integralM   = integral_over_sphere(level, x0,y0,z0, ncircles, npoints, r, Ind("puncture_properties_Mint"), order);

       integralPx   = integral_over_sphere(level, x0,y0,z0, ncircles, npoints, r, Ind("puncture_properties_Pxint"), order);
       integralPy   = integral_over_sphere(level, x0,y0,z0, ncircles, npoints, r, Ind("puncture_properties_Pyint"), order);
       integralPz   = integral_over_sphere(level, x0,y0,z0, ncircles, npoints, r, Ind("puncture_properties_Pzint"), order);
 
       integralSx   = integral_over_sphere(level, x0,y0,z0, ncircles, npoints, r, Ind("puncture_properties_Sxint"), order);
       integralSy   = integral_over_sphere(level, x0,y0,z0, ncircles, npoints, r, Ind("puncture_properties_Syint"), order);
       integralSz   = integral_over_sphere(level, x0,y0,z0, ncircles, npoints, r, Ind("puncture_properties_Szint"), order);

       /* take care of symmetries*/
	  if (Getv("grid", "bitant")){
	    integralSx = 0;
	    integralSy = 0;
            integralPz = 0;
	  } else if (Getv("grid", "quadrant")){
	    integralPx = 0;
	    integralPy = 0;
	    integralPz = 0;
	    integralSx = 0;
	    integralSy = 0;
	  } else if (Getv("grid", "qreflect")){
	    integralPy = 0;
	    integralPz = 0;
	    integralSx = 0;
	    integralSy = 0;
	    integralSz = 0;
	  } else if (Getv("grid", "rotant")) {
	    integralPx = 0;
	    integralPy = 0;
	    integralSx = 0;
	    integralSy = 0;
	  }
        
	/* process 0 does the writing */ 
         if (!rank) {
	   sprintf(name_psp,   "%s/coordsphere_%s_%s.l%d", outdir, jname, rname, l);
               fp_psp = fopen(name_psp,   "ab");
	   if (!fp_psp) errorexit("puncture_properties: failed opening file");
	    if (level->iteration == 0) {
	     fprintf(fp_psp, " coordinate sphere :   r = %14.6f \"\n",  r);
	     fprintf(fp_psp, " %14s%14s%14s%14s%14s%14s%14s%14s%14s%14s%14s \n", 
                     "time    ","posx    ","posy    ","posz    ","M      ","Px      ","Py      ","Pz      ","Sx      ","Sy      ","Sz      ");
	     } 
      
	   fprintf(fp_psp, "%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n", 
		  level->time, p[j][0],p[j][1],p[j][2],integralM,integralPx,integralPy,integralPz,integralSx,integralSy,integralSz);
           fclose(fp_psp);
         } /* rank 0 */
      } /* close loop over rlist*/ 
    } /* close loop over level */
   } /* close loop over punctures */
 
    /* disable temporary storage */
     for (l = lmax; l >= lmin; l--) {
      vl->level = level; 
      disablevarlist(vl);
     }
  
  /* finish */
  vlfree(vl);
  timer_stop(0, "puncture_properties");
  return 0;
}

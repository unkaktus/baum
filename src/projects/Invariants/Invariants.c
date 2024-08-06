/* Invariants.c */
/* Bernd Bruegmann 12/2005 */
/* Jose Gonzalez 01.06 */
/* Wolfgang Tichy 3/2010 */
/* mth 1/2011 */

#include "bam.h"
#include "Invariants.h"

#define PR 0



/* origin of spherical coordinates (symbol y0 already used in math.h */
double X0, Y0, Z0;


/* this is a part of bam_Invariants.c */
void AddVarModes()
{
    char mode_varR[100], mode_varI[100];
    
    int mode_lmin = Geti("invariants_modes_mode_lmin");
    int mode_lmax = Geti("invariants_modes_mode_lmax");
    if ((mode_lmin < 2) || (mode_lmin < 2) || (mode_lmax < mode_lmin))
        errorexit("check invariants_mode_lmin/max");

    int mode_l;
    int mode_m;
    for(mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
        for(mode_m = -mode_l; mode_m <= mode_l; mode_m++) {
            
            if (!check_mode_for_output(mode_l, mode_m)) continue;
            
            sprintf(mode_varR,"rpsi4model%d%d", mode_l, mode_m);
            sprintf(mode_varI,"ipsi4model%d%d", mode_l, mode_m);

            AddVar(mode_varR, "", mode_varR);
            AddVar(mode_varI, "", mode_varI);

        }
    }
}

/* compute curvature invariants
   registered in ANALYZE, wrapper for Mathematica generated function
*/
int compute_curvature_invariants(tL *level)
{
  int timer = timer_start(0, "invariants");
  tVarList *vl;   // list of rpsi0 etc
  tVarList *wl;   // vl combined with gxx etc
  tG *grid = level->grid;
  int lmin = level->l;
  int lmax = grid->lmax;
  int i, l, n, n0;
  char *rstring,*rstrings;
  double r, rlist[100];
  
  // output directory
  static int firstcall = 1;
  char outdir[1000];
  sprintf(outdir,"%s/Invariants", Gets("outdir"));
  if (processor0 && firstcall==1) {
    system_mkdir(outdir);
    firstcall = 0;
  }
  
  /* since we want to post-process the invariants (compute modes),
     we wait until this level is the coarsest at its time and 
     then deal with the entire time-aligned stack
  */
  if (parentaligned(level)) return 0;

  /* set origin of spherical coordinates */
  set_origin(level->time, &X0, &Y0, &Z0);

  /* variables we may want */
  vl = vlalloc(0);
//   vlpush(vl, Ind("rpsi0"));
//   vlpush(vl, Ind("ipsi0"));
//   vlpush(vl, Ind("rpsi1"));
//   vlpush(vl, Ind("ipsi1"));
//   vlpush(vl, Ind("rpsi2"));
//   vlpush(vl, Ind("ipsi2"));
//   vlpush(vl, Ind("rpsi3"));
//   vlpush(vl, Ind("ipsi3"));
  vlpush(vl, Ind("rpsi4"));
  vlpush(vl, Ind("ipsi4"));
//   vlpush(vl, Ind("rI"));
//   vlpush(vl, Ind("iI"));
//   vlpush(vl, Ind("rJ"));
//   vlpush(vl, Ind("iJ"));
//   vlpush(vl, Ind("rS"));
//   vlpush(vl, Ind("iS"));
  vlpush(vl, Ind("Csqr"));
  
  wl = vlalloc(0);
  vlpush(wl, Ind("adm_gxx"));
  vlpush(wl, Ind("adm_Kxx"));
  vlpush(wl, Ind("alpha"));
  vlpush(wl, Ind("betax"));
  vlpush(wl, Ind("x"));
  vlpush(wl, Ind("y"));
  vlpush(wl, Ind("z"));
  vlpushvl(wl, vl);

  /* enable variable we want to output -> we have to store it */
  for (l = grid->lmax; l >= grid->lmin; l--) {
    level = grid->level[l];
    vl->level = level;
    for (i = 0; i < vl->n; i++)
      enablevar(level, vl->index[i]);
  }

  /* for all levels from finest to coarse */
  for (l = lmax; l >= lmin; l--) {
    level = grid->level[l];
     
    /* make sure storage is enabled */
    vl->level = level;
    wl->level = level;

 
    /* Check if explicitly Invariants output time is enables*/
    if(Getd("Invariants_output_time") > 0){
      if (!timeforoutput_dt(level, Getd("Invariants_output_time"))) {
       lmin = l+1;
       break;
      }
    } else{
      /* check whether these variables are wanted, and are wanted now 
         if not, then the coarsest level for output is one finer, l+1
         if lmin > lmax, then no output happens at all
      */
      if (!timeforoutput(level, vl)) {
        lmin = l+1;
        break;
      }
    }

    /* call the Mathematica generated function to do the computation */
    if (Getv("invariants_EB", "yes")){
      //curvature_invariants_4_EB(wl);
      errorexit("does not exist");
    }

    curvature_invariants_N(wl);

    /* compute the analytic Psi4 for the Teukolsky wave in the 
       Kinnersley tetrad */
    if (Getv("invariants_analytic_Psi4","yes")){
      analytic_Psi4_kinnersley_teukolskywave(level, outdir);
    }

    /* set boundaries */
    /* FIXME: what are the correct symmetries for pseudo scalars? */
    set_boundary_symmetry(level, vl);
    
    /* synchronize */
    bampi_vlsynchronize(vl);
  }


  /* now we have the invariants on all the available time-aligned levels
     - set boundary from parent
     - set interior from children
  */
  for (l = lmax; l > lmin; l--) {
    level = grid->level[l];
    vl->level = level;
    set_boundary_refinement_cf(grid, l-1, l, vl); 
  }

  /* post process */
  for (l = lmax; l >= lmin; l--) {
    level = grid->level[l];
    if (PR) printf("invariants level %d, time %6.3f\n", level->l, level->time); 

    /* rescale with r */
    if (Getv("invariants_rescale", "psi4")) {
      rescale_psi4(level);
    }

    /* go through radii and see what we have to do */
    rstrings = Gets("invariants_modes_r");
    
    n = 0;
    n0= 0;
    while (rstring = NextEntry(rstrings)) {
      r = atof(rstring);
      if (PR) printf("rstring %s,  r = %f\n", rstring, r);
      if (PR) printbbox(level, level->bbox, 0);
      if (sphere_safely_inside_level(level, X0, Y0, Z0, r)) {
        if (PR) printf("computing modes for level %d, r %.1f\n", level->l, r);
        rlist[n0+n++] = r;
      } else {
        if (n==0) n0++;
      }
    }
    
    /* if the stack is not empty */
    if (n) {
      
      if ((Getv("invariants_compute_modes","yes")) &&
           level->l >= Geti("invariants_modes_lmin") &&
           level->l <= Geti("invariants_modes_lmax"))
        compute_modes(level, n0,n0+n, rlist, outdir);
      
      if (Getv("invariants_compute_energy","yes"))
        compute_energy(level, n0,n0+n, rlist, outdir);
      
      output_invariants(level, n0,n0+n, rlist);
      
    }
  }
  
  
  
 
  /* free all related rpsi4 variables --- debug mode */
  for (l = grid->lmax; l >= grid->lmin; l--) {
    level = grid->level[l];
    vl->level = level;
    for (i = 0; i < vl->n; i++)
      if(Getd("Invariants_output_time") > 0){
        if (!timeforoutput_dt(level, Getd("Invariants_output_time"))) 
          disablevar(level, vl->index[i]);
      }else{
        if (!timeforoutput_index(level, vl->index[i]))   
          disablevar(level, vl->index[i]);
      }
  }
  /* finish */
  vlfree(vl);
  vlfree(wl);
  
  timer_stop(0, "invariants");
  return 0;
}












/* determine origin of extraction sphere */
void set_origin(double t, double *x0, double *y0, double *z0)
{
    *x0 = Getd("invariants_x0") + t * Getd("invariants_vx0");
    *y0 = Getd("invariants_y0") + t * Getd("invariants_vy0");
    *z0 = Getd("invariants_z0") + t * Getd("invariants_vz0");
}


/* rescale, should be made more general */
void rescale_psi4(tL *level)
{
    double *xxx = Ptr(level, "x");
    double *yyy = Ptr(level, "y");
    double *zzz = Ptr(level, "z");
    double *rpsi4 = Ptr(level, "rpsi4");
    double *ipsi4 = Ptr(level, "ipsi4");
    double x, y, z, r;
    int i;

    forallpoints(level, i) {
        x = xxx[i] - X0;
        y = yyy[i] - Y0;
        z = zzz[i] - Z0;
        r = sqrt(x*x + y*y + z*z);
        rpsi4[i] *= r;
        ipsi4[i] *= r;
    }

}


/* analytic computation of psi4 in the Kinnersley tetrad */
void analytic_Psi4_kinnersley_teukolskywave(tL *level, char *outdir)
{

  double xana = Getd("invariants_xana");
  double yana = Getd("invariants_yana");
  double zana = Getd("invariants_zana");
  double amp = Getd("Teukolsky_wave_amplitude");
  double lambda = Getd("Teukolsky_wave_width");
  double t=level->time,theta,r,psi4_analytic;
  int rank = bampi_rank();

  char name[201];
  FILE *fp;

  r = sqrt(xana*xana+yana*yana+zana*zana);
  theta = acos(zana/r);

  psi4_analytic = 3.0/8.0*(-1.0+pow(cos(theta),2.0))*amp*(-2640.0*r*r*r*r*r*lambda*
   lambda*lambda*lambda*t*t-288.0*r*r*r*r*lambda*lambda*t*t*t*t*t-64.0*pow(r,11.0)
   -64.0*lambda*lambda*r*r*r*t*t*t*t*t*t+6.0*lambda*lambda*lambda*lambda*lambda*
   lambda*lambda*lambda*r*r*r-1248.0*r*r*r*r*r*r*r*lambda*lambda*lambda*lambda+
   5760.0*r*r*r*r*r*r*r*lambda*lambda*t*t-5440.0*r*r*r*r*r*r*lambda*lambda*t*t*t+
   6.0*lambda*lambda*lambda*lambda*lambda*lambda*lambda*lambda*t*t*t-144.0*lambda*
   lambda*lambda*lambda*lambda*lambda*r*r*t*t*t-24.0*r*r*r*r*lambda*lambda*lambda*
   lambda*lambda*lambda*t+240.0*lambda*lambda*lambda*lambda*r*r*r*t*t*t*t-24.0*
   lambda*lambda*lambda*lambda*lambda*lambda*r*t*t*t*t-2240.0*r*r*r*r*r*r*r*t*t*t*
   t+2240.0*r*r*r*r*r*r*r*r*t*t*t-1344.0*r*r*r*r*r*r*r*r*r*t*t+448.0*pow(r,10.0)*t
   +64.0*r*r*r*r*t*t*t*t*t*t*t-9.0*pow(lambda,10.0)*r-9.0*pow(lambda,10.0)*t-448.0
   *r*r*r*r*r*t*t*t*t*t*t+608.0*r*r*r*r*r*r*r*r*r*lambda*lambda+1344.0*r*r*r*r*r*r
   *t*t*t*t*t+2400.0*r*r*r*r*r*lambda*lambda*t*t*t*t+240.0*r*r*r*r*lambda*lambda*
   lambda*lambda*t*t*t+54.0*lambda*lambda*lambda*lambda*lambda*lambda*lambda*
   lambda*r*r*t-144.0*lambda*lambda*lambda*lambda*lambda*lambda*r*r*r*t*t-2976.0*r
   *r*r*r*r*r*r*r*lambda*lambda*t+54.0*lambda*lambda*lambda*lambda*lambda*lambda*
   lambda*lambda*r*t*t+48.0*lambda*lambda*lambda*lambda*r*r*t*t*t*t*t+3360.0*r*r*r
   *r*r*r*lambda*lambda*lambda*lambda*t+336.0*lambda*lambda*lambda*lambda*lambda*
   lambda*r*r*r*r*r)/pow(lambda,14.0)/(r*r*r*r*r)*exp(-pow(t-r,2.0)/(lambda*lambda
   ))-3.0/8.0*(-1.0+pow(cos(theta),2.0))*amp*(-18.0*lambda*lambda*lambda*lambda*
   lambda*lambda*lambda*lambda*r*r*t-6.0*lambda*lambda*lambda*lambda*lambda*lambda
   *lambda*lambda*r*r*r-6.0*lambda*lambda*lambda*lambda*lambda*lambda*lambda*
   lambda*t*t*t+9.0*pow(lambda,10.0)*r-18.0*lambda*lambda*lambda*lambda*lambda*
   lambda*lambda*lambda*r*t*t+9.0*pow(lambda,10.0)*t)/pow(lambda,14.0)/(r*r*r*r*r)
   *exp(-pow(t+r,2.0)/(lambda*lambda));
    
  if (Getv("invariants_rescale", "psi4")) psi4_analytic*= r;
  
  if (!rank){
	
    sprintf(name, "%s/rpsi4_analytic.l%d", outdir,level->l);
    
    fp = fopen(name, "wb");
    
    fprintf(fp, "%14.6e%14.6e\n", level->time, psi4_analytic);
    
    fclose(fp);
  }

}






int check_mode_for_output(int l, int m) 
{
    if (Getv("invariants_modes_output","all")) return 1;
    
    int a,b;
    char *rstring;
    char *rstrings = Gets("invariants_modes_output");
    while (rstring = NextEntry(rstrings)) {
        sscanf(rstring, "%d,%d",&a,&b);
        if (0) printf("%s  =>  %d %d\n",rstring, a,b);
        if ((l==a) && (m==b)) return 1;
    }
    
    return 0;
}

/* compute modes given psi_4 from l = mode_lmin to l= mode_lmax */
/* P.G. 2009 mth 01/2011 */
void compute_modes(tL *level, int n0,int nr, double *rlist, char *outdir)
{
    int rank = bampi_rank();
    int n;
    int bitant = level->grid->symmetric[2];
    int ntheta = Geti("invariants_ntheta");
    int nphi   = Geti("invariants_nphi");
    
    int order  = Geti("invariants_order");

    int mode_lmin = Geti("invariants_modes_mode_lmin");
    int mode_lmax = Geti("invariants_modes_mode_lmax");
    int lnum = mode_lmax-mode_lmin+1;
    int mnum = mode_lmax*2+1;
    char mode_varR[100], mode_varI[100];
    char mode_var[100];
    int mode_l, mode_m;
    double *rpsi4mode, *ipsi4mode;
    double *rpsi4mode_p,*ipsi4mode_p, *rpsi4mode_m,*ipsi4mode_m;
    
    double integralR[nr][lnum][mnum];
    double integralI[nr][lnum][mnum];

    double *xp = Ptr(level, "x");
    double *yp = Ptr(level, "y");
    double *zp = Ptr(level, "z");
    double *rpsi4 = Ptr(level, "rpsi4");
    double *ipsi4 = Ptr(level, "ipsi4");
    
    
    /* set everything to 0 */
    for (n=n0; n<nr; n++) 
        for (mode_l = 0; mode_l < lnum; mode_l++)
            for (mode_m = 0; mode_m < mnum; mode_m++) 
                integralR[n][mode_l][mode_m] = integralI[n][mode_l][mode_m] = 0;
    
    
    
    /* compute first the integrand */
    for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
        for( mode_m = -mode_l; mode_m <= mode_l; mode_m++) {
    
            if (!check_mode_for_output(mode_l, mode_m)) continue;
            
            sprintf(mode_varR,"rpsi4model%d%d", mode_l, mode_m); 
            rpsi4mode  = PtrEnable(level, mode_varR);

            sprintf(mode_varI,"ipsi4model%d%d", mode_l, mode_m); 
            ipsi4mode  = PtrEnable(level, mode_varI);
            
            // compute 3d integrands 
            // (could be limited to the vicinity of the spheres ...)
            bampi_openmp_start
            forallpoints_ijk_openmp(level) {

                //  X0, Y0 and Z0 !!global variables!! :S 
                double x = xp[ijk] - X0;
                double y = yp[ijk] - Y0;
                double z = zp[ijk] - Z0;
                double r = sqrt(x*x + y*y + z*z);
                
                double costheta = z/r;

                double phi;
                if(x>=0) {
                  phi = atan(y/x);
                } else {
                  phi = PI + atan(y/x);
                }

                // use the new function
                double rYm2, iYm2;
                spinweightedSphericalHarmonic(&rYm2,&iYm2, mode_l,mode_m, phi,costheta);
                
                
                // US, BB 24.3.06
                //   There is a choice of sign here: define the inner product by
                //   (f,g) = int fbar g
                //   and define
                //   psi_mode = (Y, psi)
                rpsi4mode[ijk]  =  rYm2 * rpsi4[ijk] +  iYm2 * ipsi4[ijk];
                ipsi4mode[ijk]  =  rYm2 * ipsi4[ijk] -  iYm2 * rpsi4[ijk];

            } endfor_ijk_openmp;
            bampi_openmp_stop
        }
    }

    /* do something if we have bitant symmetry */
    /* Y^s_{l m}( Pi-th, ph ) = (-1)^{l+s} Y^s_{l -m}(th, ph) */
    if (bitant) {
        double m1p;
        for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
            for( mode_m = 1; mode_m <= mode_l; mode_m++) {
                
                if (!check_mode_for_output(mode_l, mode_m) || 
                    !check_mode_for_output(mode_l, -mode_m)) continue;
                
                sprintf(mode_var,"rpsi4model%d%d", mode_l, mode_m); 
                rpsi4mode_p  = Ptr(level, mode_var);

                sprintf(mode_var,"rpsi4model%d%d", mode_l, -mode_m); 
                rpsi4mode_m  = Ptr(level, mode_var);
                
                sprintf(mode_var,"ipsi4model%d%d", mode_l, mode_m); 
                ipsi4mode_p  = Ptr(level, mode_var);
                
                sprintf(mode_var,"ipsi4model%d%d", mode_l, -mode_m); 
                ipsi4mode_m  = Ptr(level, mode_var);

                m1p = pow(-1., mode_l + (-2) );
                
                /* add something to the integrand */
                forallpoints_ijk(level){
                    rpsi4mode_p[ijk] += m1p * rpsi4mode_m[ijk];
                    ipsi4mode_p[ijk] -= m1p * ipsi4mode_m[ijk];
                } endfor_ijk;
            }
        }
    }
    


    /* now compute the sphere integral */
    /* for symmetry reasons we can use values with positive m 
       => go from l to -l
    */
    for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
        for( mode_m = mode_l; mode_m >= -mode_l; mode_m--) {
            
            if (!check_mode_for_output(mode_l, mode_m)) continue;
            
            sprintf(mode_varR,"rpsi4model%d%d", mode_l, mode_m); 
            rpsi4mode  = Ptr(level, mode_varR);

            sprintf(mode_varI,"ipsi4model%d%d", mode_l, mode_m); 
            ipsi4mode  = Ptr(level, mode_varI);
            
            // we now have all the integrands we need, go ahead and integrate for each radius 
            for (n = n0; n < nr; n++) {
                
                double r = rlist[n];
                
                if ( (mode_m >= 0) || (!bitant) ) {
                    integralR[n][mode_l-mode_lmin][mode_m+mode_l] =  
                            integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind(mode_varR), order);
                    integralI[n][mode_l-mode_lmin][mode_m+mode_l] =  
                            integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind(mode_varI), order);
                } 
                
                if ( (mode_m >  0) && (bitant) ) {
                    integralR[n][mode_l-mode_lmin][mode_m+mode_l]  /= 2;
                    integralI[n][mode_l-mode_lmin][mode_m+mode_l]  /= 2;
                    
                    integralR[n][mode_l-mode_lmin][-mode_m+mode_l]  =  integralR[n][mode_l-mode_lmin][mode_m+mode_l];
                    integralI[n][mode_l-mode_lmin][-mode_m+mode_l]  = -integralI[n][mode_l-mode_lmin][mode_m+mode_l];
                }
                
            }
        }

        
        
    }
    

    /*************************************************************************/
    /* output the integral at each radius in one file */
    if (rank==0) {
        
      /* test if file exist, if not, create and write header */
      if (level->iteration == 0) {
        for (n = n0; n < nr; n++) {
          char filename[1024];
          FILE *fp;

          sprintf(filename,"%s/rpsi4_modes_r%d.l%d",  outdir,n+1, level->l);
          fp = fopen(filename, "ab");
          if (!fp) errorexit("compute_modes: failed opening files");
          fprintf(fp,   "\"Time              ");
          for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
            for( mode_m = -mode_l; mode_m <= mode_l; mode_m++) {
              if (!check_mode_for_output(mode_l, mode_m)) continue;
              fprintf(fp,   "l=%d m=%d    \t",  mode_l, mode_m);
            }
          }
          fprintf(fp,   "\n");
          fclose(fp);
          
          sprintf(filename,"%s/ipsi4_modes_r%d.l%d",  outdir,n+1, level->l);
          fp = fopen(filename, "ab");
          if (!fp) errorexit("compute_modes: failed opening files");
          fprintf(fp,   "\"Time              ");
          for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
            for( mode_m = -mode_l; mode_m <= mode_l; mode_m++) {
              if (!check_mode_for_output(mode_l, mode_m)) continue;
              fprintf(fp,   "l=%d m=%d    \t",  mode_l, mode_m);
            }
          }
          fprintf(fp,   "\n");
          fclose(fp);
        }
      }
      
      /* files should exist, now write all needed modes */
      for (n = n0; n < nr; n++) {
          char rname[10];
          char nameR[1024], nameI[1024];
          
          double r = rlist[n];
          sprintf(rname, "r%d", n+1);
          
          sprintf(nameR,"%s/rpsi4_modes_%s.l%d",  outdir,rname, level->l);
          sprintf(nameI,"%s/ipsi4_modes_%s.l%d",  outdir,rname, level->l);
          
          FILE *fpR   = fopen(nameR, "ab");
          FILE *fpI   = fopen(nameI, "ab");
          if (!fpR || !fpI )
            errorexit("compute_modes: failed opening files");
          
          fprintf(fpR,   "%14.16e\t", level->time);
          fprintf(fpI,   "%14.16e\t", level->time);
          
          for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
              for( mode_m = -mode_l; mode_m <= mode_l; mode_m++) {
                  if (!check_mode_for_output(mode_l, mode_m)) continue;
                  fprintf(fpR,   "%14.16e\t",integralR[n][mode_l-mode_lmin][mode_m+mode_l]);
                  fprintf(fpI,   "%14.16e\t",integralI[n][mode_l-mode_lmin][mode_m+mode_l]);
              }
          }

          fprintf(fpR,   "\n");
          fprintf(fpI,   "\n");
          
          fclose(fpR);
          fclose(fpI);
          
      }
    }
    
    
    /*************************************************************************/
    /* ARNR output format of the integral, one file per mode */
    if (rank==0) {
      for (n = n0; n < nr; n++) {
        char rname[10];
        double r = rlist[n];
        sprintf(rname, "r%d", n+1);
        
        
        for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
          for( mode_m = -mode_l; mode_m <= mode_l; mode_m++) {
            if (!check_mode_for_output(mode_l, mode_m)) continue;
            char filename[1024];
            sprintf(filename,"%s/Rpsi4mode%d",  outdir, mode_l);
            if (mode_m<0) sprintf(filename,"%sm",filename);
            sprintf(filename,"%s%d_%s.l%d", filename,(mode_m>0)?mode_m:-mode_m, rname, level->l);
            
            FILE *fp = fopen(filename, "ab");
            if (!fp) errorexit("compute_modes: failed opening files");
            
            if (level->iteration == 0) {
              fprintf(fp,"\" Rpsi4:   r = %14.6f \"\n",rlist[n]);
              fprintf(fp,"\"%20s%20s%20s\"\n",
                      "time     ","Re r*Psi4  ","Im r*Psi4  "); 
            }
    
            fprintf(fp,   "%20.16e %20.16e %20.16e\n",
                    level->time,
                    integralR[n][mode_l-mode_lmin][mode_m+mode_l],
                    integralI[n][mode_l-mode_lmin][mode_m+mode_l]);
            
            fclose(fp);
          }
        }

      }
    }
    

    /*************************************************************************/
    /* disable temporary storage if wanted */
    if (!Getv("invariants_persist", "yes")) {
        
      for( mode_l = mode_lmin; mode_l <= mode_lmax; mode_l++) {
        for( mode_m = -mode_l; mode_m <= mode_l; mode_m++) {
      
          if (!check_mode_for_output(mode_l, mode_m)) continue;
      
          sprintf(mode_varR,"rpsi4model%d%d", mode_l, mode_m);
          sprintf(mode_varI,"ipsi4model%d%d", mode_l, mode_m);
          
          PtrDisable(level,mode_varR);
          PtrDisable(level,mode_varI);
        }
      }
    }
    
}








/* compute energy given psi_4 */
void compute_energy(tL *level, int n0,int nr, double *rlist, char *outdir)
{
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  
  double *rpsi4 = Ptr(level, "rpsi4");
  double *ipsi4 = Ptr(level, "ipsi4");
  double *int_rpsi4 = PtrEnable(level, "int_rpsi4");
  double *int_ipsi4 = PtrEnable(level, "int_ipsi4");
  double *int2_rpsi4 = PtrEnable(level, "int2_rpsi4");
  double *int2_ipsi4 = PtrEnable(level, "int2_ipsi4");
  double *rpsi4_p = PtrEnable(level, "rpsi4_p");
  double *ipsi4_p = PtrEnable(level, "ipsi4_p");
  double *integrand_dEdt = PtrEnable(level, "integrand_dEdt");
  double *integrand_dPxdt = PtrEnable(level, "integrand_dPxdt");
  double *integrand_dPydt = PtrEnable(level, "integrand_dPydt");
  double *integrand_dPzdt = PtrEnable(level, "integrand_dPzdt");
  double *integrand_dJzdt = PtrEnable(level, "integrand_dJzdt");
  double *integrand_One = PtrEnable(level, "integrand_One");
  double *dint_rpsi4dphi = PtrEnable(level, "dint_rpsi4dphi");
  double *dint_ipsi4dphi = PtrEnable(level, "dint_ipsi4dphi");
  double *dint2_rpsi4dphi = PtrEnable(level, "dint2_rpsi4dphi");
  double *dint2_ipsi4dphi = PtrEnable(level, "dint2_ipsi4dphi");

  
  double t_new,dt;
  double x,y,z, r;
  double dEdt,dPxdt,dPydt,dPzdt,dJzdt;
  char name1[201],name2[201],name3[201],name4[201];
  char name5[201],name6[201],name7[201],name8[201];
  char name9[201],name10[201];
  char nameP[201];
  FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10,*fpP;
  int ntheta = Geti("invariants_ntheta");
  int nphi   = Geti("invariants_nphi");
  int order  = Geti("invariants_order");
  
  int order_centered  = Geti("order_centered");
  
  int rank = bampi_rank();
  double rad,sintheta,costheta,sinphi,cosphi;
  int l = level->l;
  int bitant = level->grid->symmetric[2];
  int rotant = Getv("grid", "rotant") || Getv("grid", "quadrant");
  char rname[10];
  double int_rpsi4_p,int_ipsi4_p;
  int n;
  
  double dx = level->dx;
  double dy = level->dy;
  double oo2dx = 1.0/(2.0*dx);
  double oo2dy = 1.0/(2.0*dy);
  double Rdx,Rdy,Idx,Idy;

  /* there are a number of scalars that have to be stored between calls
     store them all in one standard 3d variable so that checkpointing works
  */
  double *store = PtrEnable(level, "invariants_storage"); 
  int nrmax = 10; /* <<<< SERIOUS PROBLEM */
  int nvars = 12;
  int i = 0;
  double *dEdt_p  = &store[nrmax*(i++)];
  double *dPxdt_p = &store[nrmax*(i++)];
  double *dPydt_p = &store[nrmax*(i++)];
  double *dPzdt_p = &store[nrmax*(i++)];
  double *dJzdt_p = &store[nrmax*(i++)];
  double *Energy  = &store[nrmax*(i++)];
  double *Px      = &store[nrmax*(i++)];
  double *Py      = &store[nrmax*(i++)];
  double *Pz      = &store[nrmax*(i++)];
  double *Jz      = &store[nrmax*(i++)];
  double *t_prev  = &store[nrmax*(i++)];

  double *AverageLapse  = &store[nrmax*(i++)];

  //if (nr >= nrmax) errorexiti("compute_energy: nr >= %d", nrmax);
  if (nvars * nrmax > level->npoints)
    errorexit("compute_energy: problem storing data\n");

  
  
  
  
  /* compute the integrals */
  if (level->iteration == 0) {

    // set all vaues to 0 
    forallpoints_ijk(level) {

      store[ijk] = 0.0;
    
      int_rpsi4[ijk] = 0.0;
      int_ipsi4[ijk] = 0.0;

      int2_rpsi4[ijk] = 0.0;
      int2_ipsi4[ijk] = 0.0;
      
      rpsi4_p[ijk] = rpsi4[ijk];
      ipsi4_p[ijk] = ipsi4[ijk];
      
    } endfor_ijk;

    t_prev[0] = level->time;

  } else {

    t_new = level -> time;
    dt = t_new - t_prev[0];
    
    double rescale = (double)(Getv("invariants_rescale","psi4"));
    
    // set prev psi4, compute both integrals
    forallpoints_ijk(level) {

      x = xp[ijk] - X0;
      y = yp[ijk] - Y0;
      z = zp[ijk] - Z0;
      rad = sqrt(x*x + y*y + z*z);
      sintheta = sqrt(1.0-(z/rad)*(z/rad));
      costheta = z/rad;
      sinphi = y/(rad*sintheta);
      cosphi = x/(rad*sintheta);
      
      int_rpsi4_p = int_rpsi4[ijk];
      int_ipsi4_p = int_ipsi4[ijk];

      // ATTENTION: if we recale rpsi4 (Re Psi_4) is actually rrpsi4 (Re rPsi_4)
      int_rpsi4[ijk] += dt/2.0*(rpsi4_p[ijk] + rpsi4[ijk]) * ( rescale/rad + (1.-rescale) );
      int_ipsi4[ijk] += dt/2.0*(ipsi4_p[ijk] + ipsi4[ijk]) * ( rescale/rad + (1.-rescale) );
      
      int2_rpsi4[ijk] += dt/2.0*(int_rpsi4_p + int_rpsi4[ijk]);
      int2_ipsi4[ijk] += dt/2.0*(int_ipsi4_p + int_ipsi4[ijk]);
      
      rpsi4_p[ijk] = rpsi4[ijk];
      ipsi4_p[ijk] = ipsi4[ijk];
      
    } endfor_ijk;
    
  }
  
  
  
  /* now compute the derivative of one time integral */
#if 0
  /* ------------------------------------------------------ */
  /*  old code, not working with shells, but the rest is ok */
  /* ------------------------------------------------------ */
  if (level->shells) {
    // FIXME: better do nothing, shells are a bit complicated
    forallpoints_ijk(level) {
      dint_rpsi4dphi[ijk] = 0;
      dint_ipsi4dphi[ijk] = 0;
    } endfor_ijk;
  } else {
    forinnerpoints_ijk(level) {
      
      x = xp[ijk] - X0;
      y = yp[ijk] - Y0;
      z = zp[ijk] - Z0;
      rad = sqrt(x*x + y*y + z*z);
      sintheta = sqrt(1.0-(z/rad)*(z/rad));
      sinphi = y/(rad*sintheta);
      cosphi = x/(rad*sintheta);
      
      if (order_centered == 2 || boundary1away) { 
        Rdx = oo2dx*(int_rpsi4[ijk+di] - int_rpsi4[ijk-di]);
        Rdy = oo2dy*(int_rpsi4[ijk+dj] - int_rpsi4[ijk-dj]);
        Idx = oo2dx*(int_ipsi4[ijk+di] - int_ipsi4[ijk-di]);
        Idy = oo2dy*(int_ipsi4[ijk+dj] - int_ipsi4[ijk-dj]);
      } else if (order_centered == 4 || boundaryNaway(2)) { 
        Rdx = 0.16666666666666666667*oo2dx*
            (int_rpsi4[-2*di + ijk] + 8.*(-int_rpsi4[-di + ijk] + int_rpsi4[di + ijk]) - int_rpsi4[2*di + ijk]);
        Rdy = 0.16666666666666666667*oo2dy*
            (int_rpsi4[-2*dj + ijk] + 8.*(-int_rpsi4[-dj + ijk] + int_rpsi4[dj + ijk]) - int_rpsi4[2*dj + ijk]);
        Idx = 0.16666666666666666667*oo2dx*
            (int_ipsi4[-2*di + ijk] + 8.*(-int_ipsi4[-di + ijk] + int_ipsi4[di + ijk]) - int_ipsi4[2*di + ijk]);
        Idy = 0.16666666666666666667*oo2dy*
            (int_ipsi4[-2*dj + ijk] + 8.*(-int_ipsi4[-dj + ijk] + int_ipsi4[dj + ijk]) - int_ipsi4[2*dj + ijk]);
      } else if (order_centered == 6 || boundaryNaway(3)) { 
        Rdx = 0.033333333333333333333*oo2dx*
            (-int_rpsi4[-3*di + ijk] + 45.*(-int_rpsi4[-di + ijk] + int_rpsi4[di + ijk]) + 
              9.*(int_rpsi4[-2*di + ijk] - int_rpsi4[2*di + ijk]) + int_rpsi4[3*di + ijk]);
        Rdy = 0.033333333333333333333*oo2dy*
            (-int_rpsi4[-3*dj + ijk] + 45.*(-int_rpsi4[-dj + ijk] + int_rpsi4[dj + ijk]) + 
            9.*(int_rpsi4[-2*dj + ijk] - int_rpsi4[2*dj + ijk]) + int_rpsi4[3*dj + ijk]);
        Idx = 0.033333333333333333333*oo2dx*
            (-int_ipsi4[-3*di + ijk] + 45.*(-int_ipsi4[-di + ijk] + int_ipsi4[di + ijk]) + 
            9.*(int_ipsi4[-2*di + ijk] - int_ipsi4[2*di + ijk]) + int_ipsi4[3*di + ijk]);
        Idy = 0.033333333333333333333*oo2dy*
            (-int_ipsi4[-3*dj + ijk] + 45.*(-int_ipsi4[-dj + ijk] + int_ipsi4[dj + ijk]) + 
            9.*(int_ipsi4[-2*dj + ijk] - int_ipsi4[2*dj + ijk]) + int_ipsi4[3*dj + ijk]);
      } else errorexit("order not implemented");
      
      // dInt/dphi = dx/dphi * dInt/dx
      dint_rpsi4dphi[ijk] = rad*sintheta*(cosphi*Rdy - sinphi*Rdx);
      dint_ipsi4dphi[ijk] = rad*sintheta*(cosphi*Idy - sinphi*Idx);
      
    } endfor_ijk;
  }
#else 
  /* ------------------------------------------------------ */
  /* new code for shells, to be tested in detail            */
  /* ------------------------------------------------------ */
  tVarList *wl = vlalloc(level);
  vlpush(wl, Ind("int2_rpsi4"));
  vlpush(wl, Ind("int2_ipsi4"));
  vlpush(wl, Ind("dint2_rpsi4dphi"));
  vlpush(wl, Ind("dint2_ipsi4dphi"));
  dphi_int2_psi4_N(wl); 
  vlfree(wl);
#endif
  
  
  /* since we did a derivative we have to sync */
  tVarList *vl = vlalloc(level);
  vlpush(vl, Ind("dint_rpsi4dphi"));
  vlpush(vl, Ind("dint_ipsi4dphi"));
  vlpush(vl, Ind("dint2_rpsi4dphi"));
  vlpush(vl, Ind("dint2_ipsi4dphi"));
  bampi_vlsynchronize(vl);
  set_boundary_symmetry(level, vl);
  vlfree(vl);
  
  
  
  /* now comute all the integrands */
  forallpoints_ijk(level) {

    x = xp[ijk] - X0;
    y = yp[ijk] - Y0;
    z = zp[ijk] - Z0;
    rad = sqrt(x*x + y*y + z*z);
    sintheta = sqrt(1.0-(z/rad)*(z/rad));
    costheta = z/rad;
    sinphi = y/(rad*sintheta);
    cosphi = x/(rad*sintheta);
      
    int_rpsi4_p = int_rpsi4[ijk];
    int_ipsi4_p = int_ipsi4[ijk];

    integrand_dEdt[ijk] = (int_rpsi4[ijk]*int_rpsi4[ijk] + int_ipsi4[ijk]*int_ipsi4[ijk]);
    
    integrand_dPxdt[ijk] = sintheta*cosphi*integrand_dEdt[ijk];
    integrand_dPydt[ijk] = sintheta*sinphi*integrand_dEdt[ijk];
    integrand_dPzdt[ijk] = costheta*integrand_dEdt[ijk];
    
    integrand_dJzdt[ijk] = (dint2_rpsi4dphi[ijk]*int_rpsi4[ijk] + dint2_ipsi4dphi[ijk]*int_ipsi4[ijk]);
    
    integrand_One[ijk] = 1.0;  // For calculating averages

    // symmetry stuff
    if (bitant) 
        integrand_dPzdt[ijk] = 0.0;
      
    if (rotant)
        integrand_dPxdt[ijk] = integrand_dPydt[ijk] = 0.0;
    
  } endfor_ijk;
  

  
  /* integrate over sphere for each radii */
  for (i = n0; i < nr; i++){
     
    r = rlist[i];
    sprintf(rname, "r%d", i+1);
    if (PR) printf("  compute and store inveriants Energy for nr=%d r=%2.2e l=%d\n",i,r,l);

    dEdt  = r*r/(16.*PI)*integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind("integrand_dEdt"), order);
    dPxdt = r*r/(16.*PI)*integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind("integrand_dPxdt"), order);
    dPydt = r*r/(16.*PI)*integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind("integrand_dPydt"), order);
    dPzdt = r*r/(16.*PI)*integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind("integrand_dPzdt"), order);
    dJzdt = r*r/(16.*PI)*integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind("integrand_dJzdt"), order);

    double avg_lapse = integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind("alpha"), order);
    double sphere_norm = integral_over_sphere(level, X0, Y0, Z0, ntheta, nphi, r, Ind("integrand_One"), order);
    AverageLapse[i] = avg_lapse/sphere_norm;
    
    Energy[i] += dt/2.0*(dEdt_p[i]  + dEdt); 
    Px[i]     += dt/2.0*(dPxdt_p[i] + dPxdt); 
    Py[i]     += dt/2.0*(dPydt_p[i] + dPydt); 
    Pz[i]     += dt/2.0*(dPzdt_p[i] + dPzdt); 
    Jz[i]     += dt/2.0*(dJzdt_p[i] + dJzdt); 
    
    dEdt_p[i]  = dEdt;
    dPxdt_p[i] = dPxdt;
    dPydt_p[i] = dPydt;
    dPzdt_p[i] = dPzdt;
    dJzdt_p[i] = dJzdt;
    t_prev[0] = t_new;
 
  
    /* do the output */
    if (!rank){
  
      sprintf(name1, "%s/dEdt_%s.l%d", outdir,rname,l);
      sprintf(name2, "%s/Energy_%s.l%d", outdir,rname,l);
      sprintf(name3, "%s/dPxdt_%s.l%d", outdir,rname,l);
      sprintf(name4, "%s/Px_%s.l%d", outdir,rname,l);
      sprintf(name5, "%s/dPydt_%s.l%d", outdir,rname,l);
      sprintf(name6, "%s/Py_%s.l%d", outdir,rname,l);
      sprintf(name7, "%s/dPzdt_%s.l%d", outdir,rname,l);
      sprintf(name8, "%s/Pz_%s.l%d", outdir,rname,l);
      sprintf(name9, "%s/dJzdt_%s.l%d", outdir,rname,l);
      sprintf(name10, "%s/Jz_%s.l%d", outdir,rname,l);

      char AverageLapseFilename[255];
      sprintf(AverageLapseFilename, "%s/AverageLapse_%s.l%d", outdir, rname, l);


      fp1 = fopen(name1, "ab");
      fp2 = fopen(name2, "ab");
      fp3 = fopen(name3, "ab");
      fp4 = fopen(name4, "ab");
      fp5 = fopen(name5, "ab");
      fp6 = fopen(name6, "ab");
      fp7 = fopen(name7, "ab");
      fp8 = fopen(name8, "ab");
      fp9 = fopen(name9, "ab");
      fp10 = fopen(name10, "ab");

      FILE *fp_AverageLapse = fopen(AverageLapseFilename, "ab");
      
      if (PR) printf("  write: %s\n", name1);
  
      if (!fp1 || !fp2 || !fp3 || !fp4 || !fp5 || !fp6 || !fp7 || !fp8 ||
          !fp9 || !fp10)
          errorexit("compute_energy: failed opening files");
  
      fprintf(fp1, "%14.6e%14.6e\n", level->time, dEdt);
      fprintf(fp2, "%14.6e%14.6e\n", level->time, Energy[i]);
      fprintf(fp3, "%14.6e%14.6e\n", level->time, dPxdt);
      fprintf(fp4, "%14.6e%14.6e\n", level->time, Px[i]);
      fprintf(fp5, "%14.6e%14.6e\n", level->time, dPydt);
      fprintf(fp6, "%14.6e%14.6e\n", level->time, Py[i]);
      fprintf(fp7, "%14.6e%14.6e\n", level->time, dPzdt);
      fprintf(fp8, "%14.6e%14.6e\n", level->time, Pz[i]);
      fprintf(fp9, "%14.6e%14.6e\n", level->time, -dJzdt);
      fprintf(fp10, "%14.6e%14.6e\n", level->time, -Jz[i]);

      fprintf(fp_AverageLapse, "%14.6e%14.6e\n", level->time, AverageLapse[i]);
  
      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      fclose(fp4);
      fclose(fp5);
      fclose(fp6);
      fclose(fp7);
      fclose(fp8);
      fclose(fp9);
      fclose(fp10);

      fclose(fp_AverageLapse);
  
      if (!rotant) {
          double px = Px[i]/4;
          double py = Py[i]/4;
          double pz = Pz[i]/4;
          double pn = sqrt(px*px + py*py);
          double p  = sqrt(px*px + py*py + pz*pz);
  
          sprintf(nameP, "%s/P_%s.l%d", outdir,rname,l);
          fpP = fopen(nameP, "ab");
          fprintf(fpP, "%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
                  level->time, px, py, pz, pn, p);
          fclose(fpP);
      }
  
    }
  }
  
  

  /* disable temporary storage */
  if (!Getv("invariants_persist", "yes")) {
    PtrDisable(level, "integrand_dEdt");
    PtrDisable(level, "integrand_dPxdt");
    PtrDisable(level, "integrand_dPydt");
    PtrDisable(level, "integrand_dPzdt");
    PtrDisable(level, "integrand_dJzdt");
        
    PtrDisable(level, "dint_rpsi4dphi");
    PtrDisable(level, "dint_ipsi4dphi");
    PtrDisable(level, "dint2_rpsi4dphi");
    PtrDisable(level, "dint2_ipsi4dphi");
  }
}







/* psi4 output routine for post analysis */
void output_invariants(tL *level, int n0,int nr, double *rlist)
{
  int *indexarray;
  int numvars;
  int n;
  int ntheta = Geti("invariants_ntheta");
  int nphi   = Geti("invariants_nphi");
  double r;
  char rname[10];
  
  
  /* use makeoutputlist to figure out numvars and indexarray from par 2doutputr */
  makeoutputlist(level, Gets("2doutputr"), Getv("2doutputall", "yes"), &numvars, &indexarray);
  
  
  
  
  /* all radii output ... only for shells */
  if (level->shells && Getv("invariants_shells_output","yes"))
    write_level_shells3d_surface(level, numvars, indexarray);
  
  
  
  
  
  /* radial output */
  if (strlen(Gets("2doutputr")) > 0 ) {

    /* call wrapper function that decides in which format we output */
    for (n = n0; n < nr; n++) {
      r = rlist[n];
      sprintf(rname, "r%d", n+1);
      write_level_sphere(level, numvars, indexarray, ntheta, nphi, r, rname);
    }

  }
  
  
  /* free mem. for indexarray allocated by makeoutputlist */
  free(indexarray);
}





























/* hydroanalysis.c */
/* SBernuzzi 10/12 */
/* dtim      11/13 */
/* dtim      17 */

#include "bam.h"
#include "hydroanalysis.h"


#define PR 0

#define MAXAHF 4


/* perform hydro analysis
   registered in ANALYZE */
int hydroanalysis(tL *level) 
{
  int timer = timer_start(0, "hydroanalysis");

  tVarList *vl;   
  tVarList *wl; // vl combined with grhd_*, gxx etc

  tG *grid = level->grid;
  int lmin = level->l;
  int lmax = grid->lmax;
  int i,l;


  /* skip the shells */
#if (0)
  // TESTME/FIXME
  if (level->shells) {
    timer_stop(0, "hydroanalysis");
    return 0;
  }
#endif
  
  /* analyze only some levels */
#if (0)
  // TESTME/FIXME
  int hydroa_lmin = Geti("hydroanalysis_lmin");
  int hydroa_lmax = Geti("hydroanalysis_lmax");
  if (level->l < hydroa_lmin) return 0;
  if (level->l > hydroa_lmax) return 0;
  if (lmax>hydroa_lmax) lmax=hydroa_lmax;
#endif

  /* wait until this level is the coarsest at its time and 
     then deal with the entire time-aligned stack */
  if (parentaligned(level)) {
    timer_stop(0, "hydroanalysis");
    return 0;
  }
  
   /* check if we want to output at this time!*/  
     if(!timeforoutput_any(level)){
        timer_stop(0, "hydroanalysis");     	
     	  return 0;
      }
      
      
  /* alloc variables */
  vl = vlalloc(0);                  // vl     // wl
  vlpush(vl, Ind("hydroa_Dh"));     // 0      // 20
  vlpush(vl, Ind("hydroa_Db"));     // 1      // 21
  vlpush(vl, Ind("hydroa_Du"));     // 2      // 22

  vlpush(vl, Ind("hydroa_etot"));   // 3      // 23
  vlpush(vl, Ind("hydroa_uesc"));   // 4      // 24

  vlpush(vl, Ind("hydroa_vorx"));   // 5      // 25 
  vlpush(vl, Ind("hydroa_vory"));   // 6      // 26 
  vlpush(vl, Ind("hydroa_vorz"));   // 7      // 27

  vlpush(vl, Ind("hydroa_ivorx"));  // 8      // 28 
  vlpush(vl, Ind("hydroa_ivory"));  // 9      // 29 
  vlpush(vl, Ind("hydroa_ivorz"));  // 10     // 30

  vlpush(vl, Ind("hydroa_Px"));  // 11      // 31 
  vlpush(vl, Ind("hydroa_Py"));  // 12      // 32 
  vlpush(vl, Ind("hydroa_Pz"));  // 13      // 33

  vlpush(vl, Ind("hydroa_Pux"));  // 14      // 34 
  vlpush(vl, Ind("hydroa_Puy"));  // 15      // 35 
  vlpush(vl, Ind("hydroa_Puz"));  // 16      // 36

  vlpush(vl, Ind("hydroa_ut"));   // 17      // 37

  wl = vlalloc(0);
  vlpush(wl, Ind( Gets("hydroanalysis_D") )); // 0
  vlpush(wl, Ind( Gets("hydroanalysis_p") )); // 1
  vlpush(wl, Ind( Gets("hydroanalysis_r") )); // 2
  vlpush(wl, Ind( Gets("hydroanalysis_e") )); // 3 
  vlpush(wl, Ind( Gets("hydroanalysis_v") )); // 4 5 6
  vlpush(wl, Ind("adm_gxx"));                 // 7 
  vlpush(wl, Ind("alpha"));
  vlpush(wl, Ind("betax"));
  vlpush(wl, Ind("x"));
  vlpush(wl, Ind("y"));
  vlpush(wl, Ind("z"));
  vlpushvl(wl, vl);

  // allocate magnetic variables and add to list
  tVarList *vlgrmhd; 
  int grmhd = Getv("physics","grmhd")+ Getv("physics", "grmhdY");
  if(grmhd){
    vlgrmhd = vlalloc(0); 
    vlpush(vlgrmhd,Ind("hydroa_Emag")); // 18
    vlpush(vlgrmhd,Ind("hydroa_Pmag")); // 19
    vlpush(vlgrmhd,Ind("hydroa_B"));    // 20
    vlpush(vlgrmhd,Ind("hydroa_divB")); // 21
    vlpush(vlgrmhd,Ind("hydroa_Emagtor"));  // 22
    vlpush(vlgrmhd,Ind("hydroa_Emagpol"));  // 23
    vlpush(vlgrmhd,Ind("hydroa_Btor"));     // 24
    vlpush(vlgrmhd,Ind("hydroa_Bpol"));     // 25
    vlpush(vlgrmhd,Ind("hydroa_Fpoy"));     // 26

    vlpushvl(vl, vlgrmhd);
  }

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

    if (PR) printf(" ===> hydroanalysis compute vars level %d, time %6.3f\n", level->l, level->time); 

    /* make sure storage is enabled */
    if (PR) printf(" \t hydroanalysis enable vl\n");
    vl->level = level;
    wl->level = level;    

    /* check whether these variables are wanted, and are wanted now 
       if not, then the coarsest level for output is one finer, l+1
       if lmin > lmax, then no output happens at all */
    if (!timeforoutput(level, vl)) {
      lmin = l+1;
      break;
    }
    
    /* compute new hydro variables */
    compute_hydrovars(wl);

    /* compute magnetic variables */        
    if (grmhd)
      compute_magvars(level); 

    /* set boundaries and synchronize */
    set_boundary_symmetry(level, vl);
    bampi_vlsynchronize(vl);
    
  }
  
  
  /* now we have the hydrovars on all the selected time-aligned levels
     - set boundary from parent
     - set interior from children  */
  if (PR) printf(" ===> hydroanalysis set_boundary_refinement_cf\n");
  for (l = lmax; l > lmin; l--) {
    level = grid->level[l];
    vl->level = level;
    set_boundary_refinement_cf(grid, l-1, l, vl); 
  }

  /* post process */
  for (l = grid->lmax; l >= grid->lmin; l--) {
    level = grid->level[l];
    vl->level = level; 

    if (PR) printf(" ===> hydroanalysis post-proc level %d, time %6.3f\n", level->l, level->time); 
    
    /* post-proc hydro var */
    postproc_hydrovars(level);
    
    /* mode_projection of D*/    
    if (Getv("hydroa_mode_projection","yes"))
      modeproject_hydrovars(level);

    /* baryonic mass computation in sphere*/        
    if (Getd("hydroanalysis_Mbar_radius")>1e-6)
      Mbar_spheres_hydrovars(level);    
     
    /* ejecta on spheres */        
    if (Getv("hydroanalysis_ejecta_spheres","yes"))
      compute_ejecta_through_sphere(level);     

    /* ejecta on angles */        
    if (Getv("hydroanalysis_ejecta_angular","yes"))
      compute_ejecta_angular(level);

    /* radiation on spheres */        
    if (Getv("physics","grrhd_m1") & Getv("hydroanalysis_radiation_spheres","yes"))
      compute_radiation_through_sphere(level);

    /* set boundaries and synchronize */    
    set_boundary_symmetry(level, vl);
    bampi_vlsynchronize(vl);      
    
  }
  

  /* disable vars */
  for (l = grid->lmax; l >= grid->lmin; l--) {
    level = grid->level[l];
    vl->level = level;
    for (i = 0; i < vl->n; i++)
      if (!timeforoutput_index(level, vl->index[i]))
        disablevar(level, vl->index[i]);
  }
  

  /* finish */
  vlfree(vl);
  vlfree(wl);
  
  timer_stop(0, "hydroanalysis");
  return 0;

}




/*  macros for D0 drvts */
#define D0x(f) oo2dx*(-f[-di + ijk] + f[di + ijk]);
#define D0y(f) oo2dy*(-f[-dj + ijk] + f[dj + ijk]); 
#define D0z(f) oo2dz*(-f[-dk + ijk] + f[dk + ijk]);


/* post process hydrovars */
void postproc_hydrovars(tL *level)
{

  if (PR) printf(" \t in postproc_hydrovars()\n");

  double ahf_data[MAXAHF][6];

  double *Dh   = Ptr(level, "hydroa_Dh"); // = D at this point 
  double *Db   = Ptr(level, "hydroa_Db");
  double *Du   = Ptr(level, "hydroa_Du");
  
  double *vorx = Ptr(level, "hydroa_vorx");
  double *vory = Ptr(level, "hydroa_vory");
  double *vorz = Ptr(level, "hydroa_vorz");

  double *ivorx = Ptr(level, "hydroa_ivorx");
  double *ivory = Ptr(level, "hydroa_ivory");
  double *ivorz = Ptr(level, "hydroa_ivorz");

  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");

  double *rho  = Ptr(level, Gets("hydroanalysis_r") );
  double rATM  = Getd("hydroanalysis_rATM");
  double *mask = Ptr(level, "matter_mask");

  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
  double xh,yh,zh,rp2;

  int grmhd = Getv("physics", "grmhd")+Getv("physics", "grmhdY");
  double *Emag, *Pmag, *divB, *dBxdx, *dBydy, *dBzdz, *Bnorm, *Emagtor, *Emagpol, *Btor, *Bpol; 
  double *Fpoy; 
  if (grmhd){
    Emag = Ptr(level, "hydroa_Emag");
    Pmag = Ptr(level, "hydroa_Pmag");
    Bnorm = Ptr(level, "hydroa_B"); 
    divB = Ptr(level, "hydroa_divB");
    Emagtor = Ptr(level, "hydroa_Emagtor");
    Emagpol = Ptr(level, "hydroa_Emagpol");
    Btor = Ptr(level, "hydroa_Btor"); 
    Bpol = Ptr(level, "hydroa_Bpol"); 
    Fpoy = Ptr(level, "hydroa_Fpoy");
  }

  /* horizon data */
  char str[100];
  int nAHF = 0; 
  int ihor;

  if (PR) printf(" \t look for AH");

  for (ihor=0; ihor<MAXAHF; ihor++) {    
    
    if (PR) printf(" %d",ihor);    
    
    sprintf(str,"ahf%d_m",ihor);    

    if (ExistPar(str)) {
      
      sprintf(str,"ahf%d_x",ihor);
      ahf_data[ihor][0] = Getd(str);

      sprintf(str,"ahf%d_y",ihor);
      ahf_data[ihor][1] = Getd(str);

      sprintf(str,"ahf%d_z",ihor);
      ahf_data[ihor][2] = Getd(str);

      sprintf(str,"ahf%d_m",ihor);
      ahf_data[ihor][3] = Getd(str);

      sprintf(str,"ahf%d_r",ihor);
      ahf_data[ihor][4] = Getd(str);
      ahf_data[ihor][5] = ahf_data[ihor][4]*ahf_data[ihor][4]; // r^2
      
      // found ? (secure check)
      if ( finite(ahf_data[ihor][4]) && 
	   finite(ahf_data[ihor][3]) &&
	   (ahf_data[ihor][4]>0.) && 
	   (ahf_data[ihor][3]>0.) ) 
	nAHF++;
      
    } 
    
  }
  
  if (PR) printf(" \n");
  
  
  bampi_openmp_start;
  
  forinnerpoints_ijk_openmp(level) {
    
    /* skip atm points */
    if (rho[ijk]<rATM) {

      vorx[ijk] = 0.;
      vory[ijk] = 0.;
      vorz[ijk] = 0.;
      
      continue;
    }
    
    
    /* take derivative of vorticity 
       (D0 could be already oscillatory...) */    
    vorx[ijk] = D0y(ivorz) - D0z(ivory);        
    vory[ijk] = D0z(ivorx) - D0x(ivorz);    
    vorz[ijk] = D0x(ivory) - D0y(ivorx);
    

    /* excise rest-mass inside horizon */
    if (nAHF>0) {
      for (ihor = 0; ihor < nAHF; ihor++){
	
        xh = xp[ijk]-ahf_data[ihor][0];
        yh = yp[ijk]-ahf_data[ihor][1];
        zh = zp[ijk]-ahf_data[ihor][2];
	
        rp2 = xh*xh +yh*yh +zh*zh;
        if (rp2<ahf_data[ihor][5]) {
	  Dh[ijk] = 0.;
	  Du[ijk] = 0.;
	  Db[ijk] = 0.;
    if (grmhd){
      Emag[ijk] = 0.;
      Pmag[ijk] = 0.;
      Bnorm[ijk] = 0.;
      divB[ijk] = 0.;
      Emagtor[ijk] = 0.;
      Emagpol[ijk] = 0.;
      Btor[ijk] = 0.;
      Bpol[ijk] = 0.;
      Fpoy[ijk] = 0.;
    }
	}
	
      }
    }

  } endfor_ijk_openmp; /* loop i, j, k */
  
  bampi_openmp_stop    

}

static int firstcallrhom = 1;    
  
/* mode projection of rho*/
/* NOTE: ORIGIN is used as the center, not the puncture position. 
         What do we actually want? */ 
int modeproject_hydrovars(tL *level)
{

  if (PR) printf(" \t in modeproject_hydrovars()\n");

  /* create output directory */
  char outdir[1000];
  sprintf(outdir,"%s/rho_mode", Gets("outdir"));
  if (processor0 && firstcallrhom==1) {
    system_mkdir(outdir);
    firstcallrhom = 0;  
  }

  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");

  double *D    = Ptr(level, Gets("hydroanalysis_D") );
  double *rho  = Ptr(level, Gets("hydroanalysis_r") );
  double rATM  = Getd("hydroanalysis_rATM");  
          
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  
  double local_real,local_imag, integral_real, integral_imag, realpart, imagpart,integral;
  double rY, iY, oor, phi, costheta;
  int    lspherical, mspherical, sym;
  int    bitant = level->grid->symmetric[2];
  
  
  /*start computation */  
  for (lspherical = 0; lspherical < 5; lspherical++) {
    if(PR) printf("l %d ", lspherical);
    for (mspherical = 0; mspherical <= lspherical; mspherical++){ 
      if(PR) printf("m %d ", mspherical);
      
      integral_real = 0;
      integral_imag = 0;
      local_real = 0;
      local_imag = 0;   	
      if(PR) printf("bitant %d ",bitant);  	
      
      if ((bitant==1) && (((lspherical+mspherical) % 2)==1)){
	integral =0;       /*CHECKME*/      
      } else {   
	forinnerpoints_ijk(level) {  
          if ((rho[ijk]>rATM)&&(zp[ijk]>(-10000.*(1-bitant)))) {	
	    oor        = 1./sqrt(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk]+zp[ijk]*zp[ijk]);     
	    phi        = acos(xp[ijk]/sqrt(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk]));
	    phi        = (yp[ijk]>=0.)? (phi) : (2.*PI-phi);
	    costheta   = zp[ijk]*oor;
	    
	    SphericalHarmonicY( &rY, &iY, lspherical, mspherical, phi, costheta);
	    
	    realpart    =  D[ijk] * rY;
	    imagpart    =  D[ijk] * iY;        
	    local_real += realpart;
	    local_imag += imagpart;		
          }
	} endfor_ijk;
	
	bampi_allreduce_sum_vector(&local_real, &integral_real, 1);
	bampi_allreduce_sum_vector(&local_imag, &integral_imag, 1);  
	integral = (bitant+1)*dx*dy*dz*sqrt(integral_real*integral_real+integral_imag*integral_imag);
        
      }
      
     
      /* process 0 does the writing */ 
      FILE *fp_rhom; 
      char name_rhom[2001];
      int  rank = bampi_rank();
      
      if (!rank) {
	sprintf(name_rhom,   "%s/rho_mode%d_%d.l%d", outdir, lspherical, mspherical, level->l);
        fp_rhom = fopen(name_rhom,"ab");
	if (!fp_rhom) errorexit("hydroanalysis: failed opening file");
	
	if (level->iteration == 0) {
	  fprintf(fp_rhom, "#mode projection of rho for l=%d,m=%d \n", lspherical, mspherical);
	} 
	
	fprintf(fp_rhom, "%14.6e%14.6e\n", level->time, integral);
        fclose(fp_rhom);
      }     
      
    } // lsperical - loop
  } //mspherical - loop
  
  return 1;

}

static int firstcall_Mbar_spheres = 1;    
  
/* compute Mbar inside spheres*/
int Mbar_spheres_hydrovars(tL *level)
{

  if (PR) printf(" \t in Mbar_spheres()\n");

  /* create output directory */
  char outdir[1000];
  sprintf(outdir,"%s/Mbar_spheres", Gets("outdir"));
  if (processor0 && firstcall_Mbar_spheres==1) {
    system_mkdir(outdir);
    firstcall_Mbar_spheres = 0;  
  }

  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");

  double *D    = Ptr(level, Gets("hydroanalysis_D") );
          
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;

  double local_Mbar, integral_Mbar,integral;
  int    bitant = level->grid->symmetric[2];
  double radius = Getd("hydroanalysis_Mbar_radius"); 
  double nradius = Getd("hydroanalysis_Mbar_nradius");
  double dradius = Getd("hydroanalysis_Mbar_dradius"); 
  int numpunc = 2;
  double r;
  int nr,np,i;
  double  p[2][3];

    
  FILE *fp_Mbar; 
  char name_Mbar[2001];
  int  rank = bampi_rank();

  /* define puncture positions */
  for (np = 0; np < numpunc; np++){
    for (i = 0; i < 3; i++)
        p[np][i] = level->grid->puncpos[np][i];

    for (nr=0; nr<=nradius; nr++) {
      r = radius + nr*dradius;
      integral_Mbar = 0;
      local_Mbar = 0;

      forinnerpoints_ijk(level) { 
          if(((xp[ijk]-p[np][0])*(xp[ijk]-p[np][0])+
            (yp[ijk]-p[np][1])*(yp[ijk]-p[np][1])+
            (zp[ijk]-p[np][2])*(zp[ijk]-p[np][2]))<r*r) 	      
	      local_Mbar += D[ijk];
      }endfor_ijk;

      bampi_allreduce_sum_vector(&local_Mbar, &integral_Mbar, 1);
      integral = (bitant+1)*dx*dy*dz*integral_Mbar;
      
      if (!rank) {
	sprintf(name_Mbar,   "%s/Mbar_star%d_r%d.l%d", outdir, np, nr, level->l);
        fp_Mbar = fopen(name_Mbar,"ab");
	if (!fp_Mbar) errorexit("hydroanalysis: failed opening file");
	
	if (level->iteration == 0) {
	  fprintf(fp_Mbar, "#baryonic mass inside radius %e for star %d \n", r, np);
          fprintf(fp_Mbar, "#time x_c y_c z_c Mbar \n");
	} 
	fprintf(fp_Mbar, "%14.6e %14.6e %14.6e %14.6e %14.6e \n", level->time, p[np][0],p[np][1],p[np][2],integral);
        fclose(fp_Mbar);
      }     
      
  }  
 }
  return 1;

}


static int firstcallspheres =1 ;


void compute_ejecta_through_sphere(tL *level)
{

  if (PR) printf(" \t in ejecta on spheres_hydrovars()\n");

  /* create output directory */
  char outdir[1000];
  sprintf(outdir,"%s/ejecta_spheres", Gets("outdir"));
  if (processor0 && firstcallspheres==1) {
    system_mkdir(outdir);
    firstcallspheres = 0;  
  }	
	
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  
  double *Du  = Ptr(level, "hydroa_Du");
  double *Db  = Ptr(level, "hydroa_Db");
  double *Pux = Ptr(level, "hydroa_Pux");
  double *Puy = Ptr(level, "hydroa_Puy");
  double *Puz = Ptr(level, "hydroa_Puz");
  double *Px  = Ptr(level, "hydroa_Px");
  double *Py  = Ptr(level, "hydroa_Py");
  double *Pz  = Ptr(level, "hydroa_Pz");
  
  double *alpha  = Ptr(level, "alpha");
  double *betax  = Ptr(level, "betax");
  double *betay  = Ptr(level, "betay");
  double *betaz  = Ptr(level, "betaz");
  double *vx     = Ptr(level, "grhd_vx");
  double *vy     = Ptr(level, "grhd_vy");
  double *vz     = Ptr(level, "grhd_vz");
    
  double *integrand_dDudt  = PtrEnable(level, "hydroa_integrand_dDudt");
  double *integrand_dDbdt  = PtrEnable(level, "hydroa_integrand_dDbdt");
  double *integrand_dDbndt = PtrEnable(level, "hydroa_integrand_dDbndt");  
  double *integrand_dPuxdt = PtrEnable(level, "hydroa_integrand_dPuxdt");
  double *integrand_dPuydt = PtrEnable(level, "hydroa_integrand_dPuydt");
  double *integrand_dPuzdt = PtrEnable(level, "hydroa_integrand_dPuzdt");
  double *integrand_dPxdt  = PtrEnable(level, "hydroa_integrand_dPxdt");
  double *integrand_dPydt  = PtrEnable(level, "hydroa_integrand_dPydt");
  double *integrand_dPzdt  = PtrEnable(level, "hydroa_integrand_dPzdt");  

  /* there are a number of scalars that have to be stored between calls
     store them all in one standard 3d variable so that checkpointing works
  */
  double *store = PtrEnable(level, "hydroa_storage"); 
  int nvars = 20;
  int i = 0;
  
  // 25 is to be the maximum of allowed spheres
  double *dDudt_p  = &store[25*(i++)];
  double *dDbdt_p  = &store[25*(i++)];
  double *dDbndt_p = &store[25*(i++)];  
  double *dPuxdt_p = &store[25*(i++)];
  double *dPuydt_p = &store[25*(i++)];
  double *dPuzdt_p = &store[25*(i++)];  
  double *dPxdt_p  = &store[25*(i++)];
  double *dPydt_p  = &store[25*(i++)];
  double *dPzdt_p  = &store[25*(i++)];
 
  double *intDu    = &store[25*(i++)];
  double *intDb    = &store[25*(i++)];
  double *intDbn   = &store[25*(i++)];  
  double *intPux   = &store[25*(i++)];
  double *intPuy   = &store[25*(i++)];
  double *intPuz   = &store[25*(i++)];  
  double *intPx    = &store[25*(i++)];
  double *intPy    = &store[25*(i++)];
  double *intPz    = &store[25*(i++)];

  double *t_prev   = &store[25*(i++)];
  
  double t_new = 0; 
  double dt = 0;
  double x,y,z,r;
  double nx,ny,nz,nvbeta,nvbetan,rad;
  double dDudt,dDbdt,dDbndt,dPuxdt,dPuydt,dPuzdt,dPxdt,dPydt,dPzdt;
  double Lpoy;
  char name1[1201],name2[1201],name3[1201],name4[1201];
  char name5[1201],name6[1201],name7[1201];
  char nameP[1201];
  
  FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7;
  
  int ntheta = Geti("hydroanalysis_sphere_ntheta");
  int nphi   = Geti("hydroanalysis_sphere_nphi");
  int order  = Geti("hydroanalysis_sphere_order");
  
  int rank = bampi_rank();
  int l = level->l;
  
  /* compute the integrals */
  if (level->time == 0) {

    // set all values to zero for 0.iteration
    forallpoints_ijk(level) {
    
      store[ijk] = 0.0;    
            
    } endfor_ijk;
  }
    t_new = level->time;
    dt    = t_new - t_prev[0];
    if (PR) printf("times: %e %e %e \n", t_new, dt, t_prev[0]);
    
  /* now compute all the integrands */
  forallpoints_ijk(level) {
    
    
    x = xp[ijk] ;
    y = yp[ijk] ;
    z = zp[ijk] ;
    rad = sqrt(x*x + y*y + z*z);

    nx = (alpha[ijk]*vx[ijk]-betax[ijk])*x/rad;
    ny = (alpha[ijk]*vy[ijk]-betay[ijk])*y/rad;
    nz = (alpha[ijk]*vz[ijk]-betaz[ijk])*z/rad;
    nvbeta  = (nx+ny+nz);
    if (nvbeta > 0) 
      nvbetan = 0;
    else 
      nvbetan = nvbeta;

    integrand_dDudt[ijk]  = Du[ijk]*nvbeta;       
    integrand_dDbdt[ijk]  = Db[ijk]*nvbeta;
    integrand_dDbndt[ijk] = Db[ijk]*nvbetan;  
        
    integrand_dPuxdt[ijk] = Pux[ijk]*nvbeta;
    integrand_dPuydt[ijk] = Puy[ijk]*nvbeta;
    integrand_dPuzdt[ijk] = Puz[ijk]*nvbeta;    
    
    integrand_dPxdt[ijk]  = Px[ijk]*nvbeta;
    integrand_dPydt[ijk]  = Py[ijk]*nvbeta;
    integrand_dPzdt[ijk]  = Pz[ijk]*nvbeta;
    
  } endfor_ijk;
  
  
  double radius  = Getd("hydroanalysis_ejecta_spheres_radius"); 
  double nradius = Getd("hydroanalysis_ejecta_nradius");
  double dradius = Getd("hydroanalysis_ejecta_dradius");   
  int grmhd = Getv("physics", "grmhd")+Getv("physics", "grmhdY");
  
  for (i=0; i<=nradius; i++) {
    
    r      = radius + i*dradius;      
    
    if ( r > level-> bbox[1])    
      break;
      
    dDudt  = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dDudt"),  order);  
    dDbdt  = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dDbdt"),  order); 
    dDbndt = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dDbndt"), order);     
    dPuxdt = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dPuxdt"), order);
    dPuydt = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dPuydt"), order);
    dPuzdt = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dPuzdt"), order);
    dPxdt  = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dPxdt"),  order);
    dPydt  = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dPydt"),  order);
    dPzdt  = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dPzdt"),  order);    
    if (grmhd)
      Lpoy = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_Fpoy"),  order);
    
    //first order integration
    intDu[i]  += dt*dDudt; 
    intDb[i]  += dt*dDbdt; 
    intDbn[i]  += dt*dDbndt; 
    
    intPux[i] += dt*dPuxdt; 
    intPuy[i] += dt*dPuydt; 
    intPuz[i] += dt*dPuzdt; 
    
    intPx[i]  += dt*dPxdt; 
    intPy[i]  += dt*dPydt; 
    intPz[i]  += dt*dPzdt;     
  
    /* do the output */
    if (!rank){

      sprintf(name1, "%s/dDdt_r%d.l%d", outdir,i,l);
      sprintf(name2, "%s/D_r%d.l%d",    outdir,i,l);     
      sprintf(name3, "%s/dPudt_r%d.l%d",outdir,i,l);
      sprintf(name4, "%s/Pu_r%d.l%d",   outdir,i,l);
      sprintf(name5, "%s/P_r%d.l%d",    outdir,i,l);
      sprintf(name6, "%s/dPdt_r%d.l%d", outdir,i,l);
      if (grmhd) sprintf(name7, "%s/Lpoy_r%d.l%d", outdir,i,l);
  
      fp1 = fopen(name1, "ab");
      fp2 = fopen(name2, "ab");
      fp3 = fopen(name3, "ab");
      fp4 = fopen(name4, "ab");
      fp5 = fopen(name5, "ab");
      fp6 = fopen(name6, "ab");
      if (grmhd) fp7 = fopen(name7, "ab");
                                    
      if (PR) printf("  write: %s\n", name1);
  
      if (!fp1 || !fp2 || !fp3 || !fp4 || !fp5 || !fp6 )
          errorexit("compute_ejecta_through_sphere: failed opening files");
      if (grmhd && !fp7 )
          errorexit("compute_ejecta_through_sphere: failed opening file for Lpoy");

  	   if (level->iteration == 0) {
	    fprintf(fp1, "# ejecta computation for r=%e \n", r);
	    fprintf(fp2, "# ejecta computation for r=%e \n", r);
	    fprintf(fp3, "# ejecta computation for r=%e \n", r);
	    fprintf(fp4, "# ejecta computation for r=%e \n", r);	    
	    fprintf(fp5, "# ejecta computation for r=%e \n", r);
	    fprintf(fp6, "# ejecta computation for r=%e \n", r);	 
      if (grmhd)  fprintf(fp7, "# ejecta computation for r=%e \n", r);	
	           
       fprintf(fp1, "#time dDbdt dDudt dDbdt(vr<0)\n");
       fprintf(fp2, "#time Db Du Db(vr<0)\n");
       fprintf(fp3, "#time dPuxdt dPuydt dPuzdt \n");       
       fprintf(fp4, "#time Pux Puy Puz \n");              
       fprintf(fp5, "#time dPxdt dPydt dPzdt \n");       
       fprintf(fp6, "#time Px Py Pz \n");  
       if (grmhd) fprintf(fp7, "#time Lpoy \n");      
	   }
	
      fprintf(fp1,  "%14.6e%14.6e%14.6e%14.6e\n", level->time, dDbdt,dDudt,dDbndt);
      fprintf(fp2,  "%14.6e%14.6e%14.6e%14.6e\n", level->time, intDb[i],intDu[i],intDbn[i]);
      fprintf(fp3,  "%14.6e%14.6e%14.6e%14.6e\n", level->time, dPuxdt,dPuydt,dPuzdt);
      fprintf(fp4,  "%14.6e%14.6e%14.6e%14.6e\n", level->time, intPux[i],intPuy[i],intPuz[i]);   
      fprintf(fp5,  "%14.6e%14.6e%14.6e%14.6e\n", level->time, dPxdt,dPydt,dPzdt);
      fprintf(fp6,  "%14.6e%14.6e%14.6e%14.6e\n", level->time, intPx[i],intPy[i],intPz[i]);    
      if (grmhd)  fprintf(fp7,  "%14.6e%14.6e\n", level->time, Lpoy);  
  
      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      fclose(fp4);
      fclose(fp5);
      fclose(fp6);
      if (grmhd) fclose(fp7);
  
    }
  }
  
        t_prev[0]  = t_new;
  
}


static int firstcallspheres_radiation = 1 ;


void compute_radiation_through_sphere(tL *level)
{

  if (PR) printf(" \t in radiation on spheres_hydrovars()\n");

  /* create output directory */
  char outdir[1000];
  sprintf(outdir,"%s/radiation_spheres", Gets("outdir"));
  if (processor0 && firstcallspheres_radiation==1) {
    system_mkdir(outdir);
    firstcallspheres_radiation = 0;  
  }	
	
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  
  double *nue_E  = Ptr(level, "grrhd_m1_nue_E");
  double *nue_Fx = Ptr(level, "grrhd_m1_nue_Fx");
  double *nue_Fy = Ptr(level, "grrhd_m1_nue_Fy");
  double *nue_Fz = Ptr(level, "grrhd_m1_nue_Fz");

  double *nua_E  = Ptr(level, "grrhd_m1_nua_E");  
  double *nua_Fx = Ptr(level, "grrhd_m1_nua_Fx");
  double *nua_Fy = Ptr(level, "grrhd_m1_nua_Fy");
  double *nua_Fz = Ptr(level, "grrhd_m1_nua_Fz");

  double *nux_E  = Ptr(level, "grrhd_m1_nux_E");
  double *nux_Fx = Ptr(level, "grrhd_m1_nux_Fx");
  double *nux_Fy = Ptr(level, "grrhd_m1_nux_Fy");
  double *nux_Fz = Ptr(level, "grrhd_m1_nux_Fz");

  double *alpha  = Ptr(level, "alpha");
  double *betax  = Ptr(level, "betax");
  double *betay  = Ptr(level, "betay");
  double *betaz  = Ptr(level, "betaz");

  double *gxx    = Ptr(level, "adm_gxx");
  double *gxy    = Ptr(level, "adm_gxy");
  double *gxz    = Ptr(level, "adm_gxz");
  double *gyy    = Ptr(level, "adm_gyy");
  double *gyz    = Ptr(level, "adm_gyz");
  double *gzz    = Ptr(level, "adm_gzz");
  
  double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, detg;

  double *integrand_dEnuedt = PtrEnable(level, "hydroa_integrand_dEnuedt");
  double *integrand_dEnuadt = PtrEnable(level, "hydroa_integrand_dEnuadt");
  double *integrand_dEnuxdt = PtrEnable(level, "hydroa_integrand_dEnuxdt");  

  double *store = PtrEnable(level, "hydroa_storage"); 
  int nvars = 20;
  int i = 0;

  double x,y,z,r,r1,r2,dr;
  double fx,fy,fz,rad;
  double dEnuedt, dEnuadt, dEnuxdt;
  double idEnuedt,idEnuadt,idEnuxdt;
  char name1[201];

  FILE *fp1;
  
  int ntheta = Geti("hydroanalysis_sphere_ntheta");
  int nphi   = Geti("hydroanalysis_sphere_nphi");
  int order  = Geti("hydroanalysis_sphere_order");
  
  int rank = bampi_rank();
  int l = level->l;
  
  /* compute the integrals */
  if (level->time == 0) {

  // set all values to zero for 0.iteration
  forallpoints_ijk(level) {
    
    store[ijk] = 0.0;    
            
    } endfor_ijk;
  }
    
  /* now compute all the integrands */
  forallpoints_ijk(level) {
  
    x = xp[ijk];
    y = yp[ijk];
    z = zp[ijk];
    rad = sqrt(x*x + y*y + z*z);

	  hydroa_compute_detg_invg_pt(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk],
				     &detg,
				     &gupxx,&gupxy,&gupxz,&gupyy,&gupyz,&gupzz);

    fx = (alpha[ijk]*(gupxx*nue_Fx[ijk] + gupxy*nue_Fy[ijk] + gupxz*nue_Fz[ijk]) - nue_E[ijk]*betax[ijk]) * x/rad;
    fy = (alpha[ijk]*(gupxy*nue_Fx[ijk] + gupyy*nue_Fy[ijk] + gupyz*nue_Fz[ijk]) - nue_E[ijk]*betay[ijk]) * y/rad;
    fz = (alpha[ijk]*(gupxz*nue_Fx[ijk] + gupyz*nue_Fy[ijk] + gupzz*nue_Fz[ijk]) - nue_E[ijk]*betaz[ijk]) * z/rad;

    integrand_dEnuedt[ijk] = fx + fy +fz;     

    fx = (alpha[ijk]*(gupxx*nua_Fx[ijk] + gupxy*nua_Fy[ijk] + gupxz*nua_Fz[ijk]) - nua_E[ijk]*betax[ijk]) * x/rad;
    fy = (alpha[ijk]*(gupxy*nua_Fx[ijk] + gupyy*nua_Fy[ijk] + gupyz*nua_Fz[ijk]) - nua_E[ijk]*betay[ijk]) * y/rad;
    fz = (alpha[ijk]*(gupxz*nua_Fx[ijk] + gupyz*nua_Fy[ijk] + gupzz*nua_Fz[ijk]) - nua_E[ijk]*betaz[ijk]) * z/rad;

    integrand_dEnuadt[ijk] = fx + fy +fz;

    fx = (alpha[ijk]*(gupxx*nux_Fx[ijk] + gupxy*nux_Fy[ijk] + gupxz*nux_Fz[ijk]) - nux_E[ijk]*betax[ijk]) * x/rad;
    fy = (alpha[ijk]*(gupxy*nux_Fx[ijk] + gupyy*nux_Fy[ijk] + gupyz*nux_Fz[ijk]) - nux_E[ijk]*betay[ijk]) * y/rad;
    fz = (alpha[ijk]*(gupxz*nux_Fx[ijk] + gupyz*nux_Fy[ijk] + gupzz*nux_Fz[ijk]) - nux_E[ijk]*betaz[ijk]) * z/rad;

    integrand_dEnuxdt[ijk] = fx + fy +fz; 
    
  } endfor_ijk;
  
  
  int indx[9];

  indx[0] = Ind("hydroa_integrand_dEnuedt");
  indx[1] = Ind("hydroa_integrand_dEnuadt");
  indx[2] = Ind("hydroa_integrand_dEnuxdt");

  indx[3] = Ind("grrhd_m1_nue_J");
  indx[4] = Ind("grrhd_m1_nua_J");
  indx[5] = Ind("grrhd_m1_nux_J");

  indx[6] = Ind("grrhd_m1_nue_n");
  indx[7] = Ind("grrhd_m1_nua_n");
  indx[8] = Ind("grrhd_m1_nux_n");

  r1 = Getd("hydroanalysis_radiation_angular_radius");
  dr = Getd("hydroanalysis_radiation_angular_dr");
  r2 = r1 + dr;

  write_level_sphere0(level, 9, indx,     ntheta, nphi, r1,     "r0");
  write_level_sphere0(level, 9, indx,     ntheta, nphi, r2,     "r1");


  double radius  = Getd("hydroanalysis_ejecta_spheres_radius");
  double nradius = Getd("hydroanalysis_ejecta_nradius");
  double dradius = Getd("hydroanalysis_ejecta_dradius");

  for (i=0; i<=nradius; i++) {
    
    r      = radius + i*dradius;      
    
    if ( r > level-> bbox[1])    
      break;
      
    dEnuedt = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dEnuedt"), order);  
    dEnuadt = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dEnuadt"), order); 
    dEnuxdt = r*r*integral_over_sphere(level, 0., 0., 0., ntheta, nphi, r, Ind("hydroa_integrand_dEnuxdt"), order);     
    
    /* do the output */
    if (!rank) {

      sprintf(name1, "%s/dEnudt_r%d.l%d", outdir,i,l);

      fp1 = fopen(name1, "ab");
                                    
      if (PR) printf("  write: %s\n", name1);
  
      if (!fp1)
          errorexit("compute_radiation_through_sphere: failed opening files");

  	  if (level->iteration == 0) {
	      fprintf(fp1, "# ejecta computation for r=%e \n", r);
        fprintf(fp1, "#time dEnuedt dEnuadt dEnuxdt\n");
	    }
	
      fprintf(fp1,  "%14.6e%14.6e%14.6e%14.6e\n", level->time, dEnuedt,dEnuadt,dEnuxdt);
  
      fclose(fp1);
  
    }
  }
  
  
}


static int firstcallspheres_angular = 1 ;

void compute_ejecta_angular(tL *level)
{

  if (PR) printf(" \t in angular hydroanalysis \n");

  int     l = level->l;
  double dx = level->dx;

  double r1     = Getd("hydroanalysis_ejecta_angular_radius");   
  if ( r1 > level-> bbox[1]) return;

  double dr = Getd("hydroanalysis_ejecta_angular_dr");
  double r2 = r1 + dr;

  static char *ou;
  if (firstcallspheres_angular){
    firstcallspheres_angular = 0;
    ou = Gets("hydroanalysis_ejecta_angular_output");
  }

  if (PR) printf(" \t in ejecta on compute ejecta angula\n");

  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  
  double *D   = Ptr(level, "grhd_D");
  double *Px  = Ptr(level, "grhd_Sx");
  double *Py  = Ptr(level, "grhd_Sy");
  double *Pz  = Ptr(level, "grhd_Sz");
  double *u_t = Ptr(level, "hydroa_ut");

  //double *rho = Ptr(level, "grhd_rho");
  //double *T   = Ptr(level, "grhd_T");
  //double *Y   = Ptr(level, "grhd_Y");

  double *vx    = Ptr(level, "grhd_vx");
  double *vy    = Ptr(level, "grhd_vy");
  double *vz    = Ptr(level, "grhd_vz");
  double *v2    = Ptr(level, "grhd_v2");

  double *alpha  = Ptr(level, "alpha");
  double *betax  = Ptr(level, "betax");
  double *betay  = Ptr(level, "betay");
  double *betaz  = Ptr(level, "betaz");

  double *gxx    = Ptr(level, "adm_gxx");
  double *gxy    = Ptr(level, "adm_gxy");
  double *gxz    = Ptr(level, "adm_gxz");
  double *gyy    = Ptr(level, "adm_gyy");
  double *gyz    = Ptr(level, "adm_gyz");
  double *gzz    = Ptr(level, "adm_gzz");

  double *integrand_dDudt = PtrEnable(level, "hydroa_integrand_dDudt");
  double *integrand_dPxdt = PtrEnable(level, "hydroa_integrand_dPxdt");
  double *integrand_dPydt = PtrEnable(level, "hydroa_integrand_dPydt");
  double *integrand_dPzdt = PtrEnable(level, "hydroa_integrand_dPzdt"); 
  double *udownt = PtrEnable(level, "hydroa_ut");

  double uupt, uupx, uupy, uupz;
  double betadownx, betadowny, betadownz, betasq;
  double betaduu;
  double ooalpha, W;

  int np_proc,np_rem;
  double x, y, z;
  double nx, ny, nz, nvbeta, nvbetan, rad;
  double dDudt, dPuxdt, dPuydt, dPuzdt;
  char name[201];
    
  int ncircles  = Geti("hydroanalysis_sphere_ntheta");
  int npoints   = Geti("hydroanalysis_sphere_nphi");
  int order     = Geti("hydroanalysis_sphere_order");
  double dtheta, dphi; 
  int np_total, init_point, end_point;
  double theta_i, phi_i;

  int rank = bampi_rank();
  int micro = Getv("physics", "grhdDY")+Getv("physics", "grmhdY"); 

  /* now compute all the integrands */
  forallpoints_ijk(level) {
    
    x = xp[ijk] ;
    y = yp[ijk] ;
    z = zp[ijk] ;
    rad = sqrt(x*x + y*y + z*z);
    
    nx = (alpha[ijk]*vx[ijk]-betax[ijk])*x/rad;
    ny = (alpha[ijk]*vy[ijk]-betay[ijk])*y/rad;
    nz = (alpha[ijk]*vz[ijk]-betaz[ijk])*z/rad;
    nvbeta  = (nx+ny+nz);

    integrand_dDudt[ijk]  = D[ijk]*nvbeta;       
        
    integrand_dPxdt[ijk] = Px[ijk]*nvbeta;
    integrand_dPydt[ijk] = Py[ijk]*nvbeta;
    integrand_dPzdt[ijk] = Pz[ijk]*nvbeta;    
  
  } endfor_ijk;
  
  int nindex, *index;
  makeoutputlist(level, ou, 0, &nindex, &index);

  /* do the output */
  write_level_sphere0(level, nindex, index, ncircles, npoints, r1, "r0");
  write_level_sphere0(level, nindex, index, ncircles, npoints, r2, "r1");

  // int indx[9];

  // if (Getv("eos", "nuclear")) {
  //   indx[0] = Ind("grhd_rho");
  //   indx[1] = Ind("grhd_T");
  //   indx[2] = Ind("grhd_Y");
  // }
  // else {
  //   indx[0] = Ind("grhd_rho");
  //   indx[1] = Ind("grhd_p");
  //   indx[2] = Ind("grhd_epsl");
  // }
  // indx[3] = Ind("hydroa_integrand_dDudt");
  // indx[4] = Ind("hydroa_integrand_dPxdt");
  // indx[5] = Ind("hydroa_integrand_dPydt");
  // indx[6] = Ind("hydroa_integrand_dPzdt");

  // indx[7] = Ind("hydroa_ut");
  // indx[8] = Ind("grhd_v2");

  // /* do the output */
  // write_level_sphere0(level, 9, indx,     ncircles, npoints, r1,     "r0");
  // write_level_sphere0(level, 9, indx,     ncircles, npoints, r2,     "r1");

}


void hydroa_compute_detg_invg_pt(double g11, double g12, double g13, 
			       double g22, double g23, double g33,
			       double *det,
			       double *i11, double *i12, double *i13, 
			       double *i22, double *i23, double *i33)
{
  
  double detg, oodetg, gginv11, gginv12, gginv13, gginv22, gginv23, gginv33;
  
  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  gginv22 = g11*g33 - g13*g13;
  gginv23 = g12*g13 - g11*g23;
  gginv33 = g11*g22 - g12*g12;
  
  detg    = g11*gginv11 + g12*gginv12 + g13*gginv13;

  oodetg = 1./detg;

  *i11 = gginv11 * oodetg;
  *i12 = gginv12 * oodetg;
  *i13 = gginv13 * oodetg;
  *i22 = gginv22 * oodetg;
  *i23 = gginv23 * oodetg;
  *i33 = gginv33 * oodetg;
  *det = detg;
}

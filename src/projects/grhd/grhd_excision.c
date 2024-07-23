/* matter_excision.c */
/* mth 02/10, sbernuz 08/12 */


#include "bam.h"
#include "grhd.h"

#define PR 0
#define maxAHF 3


double excision_rmin,excision_rmax,excision_rfct,excision_fact;
double ahf_data[maxAHF][5]; //pos1,pos2,pos3,m,r ... simpler storage



// init excision
void grhd_excision_init(tL *level)
{

  if (!Getv("physics","AHmod"))
    errorexit("need AHmod for excision");
  
  // choose type of excision
  if      (Getv("grhd_excision_modus","cut"))  MATTER.set_excision = grhd_set_excision_cut;
  else if (Getv("grhd_excision_modus","lin"))  MATTER.set_excision = grhd_set_excision_lin;
  else if (Getv("grhd_excision_modus","atm"))  MATTER.set_excision = grhd_set_excision_atm;
  else if (Getv("grhd_excision_modus","vmax")) MATTER.set_excision = grhd_set_excision_vmax;
  else if (Getv("grhd_excision_modus","fixr")) MATTER.set_excision = grhd_set_excision_fixr;  
  else errorexit("the chosen excision modus does not exist");     

  // pars
  excision_rmin = Getd("grhd_excision_rmin");
  excision_rmax = Getd("grhd_excision_rmax");
  excision_rfct = Getd("grhd_excision_rfct");
  excision_fact = Getd("grhd_excision_fact");

}


// get values of i-th AHF (if found)
int get_ahf_val(int i)
{
  
  char str[100];
  int found = 0; 

  sprintf(str,"ahf%d_m",i);
  
  if (ExistPar(str)) {
      
    sprintf(str,"ahf%d_x",i);
    ahf_data[i][0] = Getd(str);

    sprintf(str,"ahf%d_y",i);
    ahf_data[i][1] = Getd(str);

    sprintf(str,"ahf%d_z",i);
    ahf_data[i][2] = Getd(str);

    sprintf(str,"ahf%d_m",i);
    ahf_data[i][3] = Getd(str);

    sprintf(str,"ahf%d_r",i);
    ahf_data[i][4] = Getd(str);
    
    if ( finite(ahf_data[i][4]) && 
	 finite(ahf_data[i][3]) &&
	 (ahf_data[i][4]>0.) && 
	 (ahf_data[i][3]>0.) ) 
      found = 1;
    
  } 
    
  return found;

}


/* excision routines */


// re-set old values if r < rmin || r > rmax 
// bad just for testing
void grhd_set_excision_fixr(tVarList *upre, tVarList *ucur)
{
  
  tL *level  = ucur->level;    

  double *PD  = vldataptr(upre, MATTER.INDX_VAR_q    );
  double *PT  = vldataptr(upre, MATTER.INDX_VAR_q + 1);
  double *PSx = vldataptr(upre, MATTER.INDX_VAR_q + 2);
  double *PSy = vldataptr(upre, MATTER.INDX_VAR_q + 3);
  double *PSz = vldataptr(upre, MATTER.INDX_VAR_q + 4);

  double *D  = vldataptr(ucur, MATTER.INDX_VAR_q    );
  double *T  = vldataptr(ucur, MATTER.INDX_VAR_q + 1);
  double *Sx = vldataptr(ucur, MATTER.INDX_VAR_q + 2);
  double *Sy = vldataptr(ucur, MATTER.INDX_VAR_q + 3);
  double *Sz = vldataptr(ucur, MATTER.INDX_VAR_q + 4);

  double *xp    = Ptr(level, "x");
  double *yp    = Ptr(level, "y");
  double *zp    = Ptr(level, "z");
  double r;  
  
  forallpoints_ijk(level) {

    r = sqrt( xp[ijk]*xp[ijk] +yp[ijk]*yp[ijk] +zp[ijk]*zp[ijk] );
    
    if ((r<excision_rmin) || (r>excision_rmax)) {

      D[ijk]  = PD[ijk];
      Sx[ijk] = PSx[ijk];
      Sy[ijk] = PSy[ijk];
      Sz[ijk] = PSz[ijk];
      T[ijk]  = PT[ijk];

    }

    
  } endfor_ijk;
  
}


// set ATM in excised zone
void grhd_set_excision_atm(tVarList *upre, tVarList *ucur)
{
  
  tL *level = ucur->level;    
  double *D  = vldataptr(ucur, MATTER.INDX_VAR_q    );
  double *T  = vldataptr(ucur, MATTER.INDX_VAR_q + 1);
  double *Sx = vldataptr(ucur, MATTER.INDX_VAR_q + 2);
  double *Sy = vldataptr(ucur, MATTER.INDX_VAR_q + 3);
  double *Sz = vldataptr(ucur, MATTER.INDX_VAR_q + 4);

  int n;
  double delta,SQRTdetg;

  double *detg = Ptr(level, "adm_detg");
  double *p    = Ptr(level, "grhd_p");

  double *x    = Ptr(level, "x");
  double *y    = Ptr(level, "y");
  double *z    = Ptr(level, "z");
  
  
  for (n=0; n<maxAHF; n++) {
    
    if (PR) printf("   look for AH%d",n);
    
    if (!(get_ahf_val(n))) {
      if (PR) printf("\n");
      continue;
    }
    
    if (1) printf("  found   => excise region  (%e %e %e   r=%e)\n",ahf_data[n][0],ahf_data[n][1],ahf_data[n][2],ahf_data[n][4]);
    
    forallpoints_ijk(level) {
      
      // look if point is inside => set excision on evolution vars and pressure
      delta = sqrt(pow(x[ijk]-ahf_data[n][0],2)+pow(y[ijk]-ahf_data[n][1],2)+pow(z[ijk]-ahf_data[n][2],2));
      
      if (delta < excision_rfct*ahf_data[n][4]) {
	
	  SQRTdetg = sqrt( DMAX( detg[ijk] , DETGMIN ) );
	  if (!finite(SQRTdetg)) SQRTdetg = 1.;
	  
	  D[ijk]  = SQRTdetg*GRHD.ATM_RHOATM;
	  Sx[ijk] = 0.;
	  Sy[ijk] = 0.;
	  Sz[ijk] = 0.;
	  T[ijk]  = SQRTdetg*GRHD.ATM_RHOATM*GRHD.ATM_EPSLATM;
	  p[ijk]  = GRHD.ATM_PATM;

      } 

    } endfor_ijk;
      
  }
}



// set ATM in excised zone if v2>vmax
void grhd_set_excision_vmax(tVarList *upre, tVarList *ucur)
{
  
  tL *level = ucur->level;    
  double *D  = vldataptr(ucur, MATTER.INDX_VAR_q    );
  double *T  = vldataptr(ucur, MATTER.INDX_VAR_q + 1);
  double *Sx = vldataptr(ucur, MATTER.INDX_VAR_q + 2);
  double *Sy = vldataptr(ucur, MATTER.INDX_VAR_q + 3);
  double *Sz = vldataptr(ucur, MATTER.INDX_VAR_q + 4);

  int n;
  double delta,SQRTdetg;

  double *detg = Ptr(level, "adm_detg");
  double *p    = Ptr(level, "grhd_p");
  double *v2   = Ptr(level, "grhd_v2");

  double *x    = Ptr(level, "x");
  double *y    = Ptr(level, "y");
  double *z    = Ptr(level, "z");
  
  
  for (n=0; n<maxAHF; n++) {
    
    if (PR) printf("   look for AH%d",n);
    
    if (!(get_ahf_val(n))) continue;
    
    if (PR) printf("  found   => excise region  (%e %e %e   r=%e)\n",ahf_data[n][0],ahf_data[n][1],ahf_data[n][2],ahf_data[n][4]);
    
    forallpoints_ijk(level) {
      
      // look if point is inside => set excision on evolution vars and pressure
      delta = sqrt(pow(x[ijk]-ahf_data[n][0],2)+pow(y[ijk]-ahf_data[n][1],2)+pow(z[ijk]-ahf_data[n][2],2));
      
      if (delta < excision_rfct*ahf_data[n][4]) {

	if (v2[ijk] > GRHD.HRSC_VMAX) {
	  
	  SQRTdetg = sqrt( DMAX( detg[ijk] , DETGMIN ) );
	  if (!finite(SQRTdetg)) SQRTdetg = 1.;
	  
	  D[ijk]  = SQRTdetg*GRHD.ATM_RHOATM;
	  Sx[ijk] = 0.;
	  Sy[ijk] = 0.;
	  Sz[ijk] = 0.;
	  T[ijk]  = SQRTdetg*GRHD.ATM_RHOATM*GRHD.ATM_EPSLATM;
	  p[ijk]  = GRHD.ATM_PATM;
	  v2[ijk] = 0.;
	  
	}
	
      } 

    } endfor_ijk;
    
  }
}


// set excision zone to *= excision_fact
void grhd_set_excision_cut(tVarList *upre, tVarList *ucur)
{
  
  tL *level = ucur->level;    
  double *D  = vldataptr(ucur, MATTER.INDX_VAR_q    );
  double *T  = vldataptr(ucur, MATTER.INDX_VAR_q + 1);
  double *Sx = vldataptr(ucur, MATTER.INDX_VAR_q + 2);
  double *Sy = vldataptr(ucur, MATTER.INDX_VAR_q + 3);
  double *Sz = vldataptr(ucur, MATTER.INDX_VAR_q + 4);

  int n;
  double delta;

  double *detg = Ptr(level, "adm_detg");
  double *p    = Ptr(level, "grhd_p");

  double *x    = Ptr(level, "x");
  double *y    = Ptr(level, "y");
  double *z    = Ptr(level, "z");
  
  
  for (n=0; n<maxAHF; n++) {
    
    if (PR) printf("   look for AH%d",n);
    
    if (!(get_ahf_val(n))) continue;
    
    if (PR) printf("  found   => excise region  (%e %e %e   r=%e)\n",ahf_data[n][0],ahf_data[n][1],ahf_data[n][2],ahf_data[n][4]);
    
    forallpoints_ijk(level) {
      
      // look if point is inside => set excision on evolution vars and pressure
      delta = sqrt(pow(x[ijk]-ahf_data[n][0],2)+pow(y[ijk]-ahf_data[n][1],2)+pow(z[ijk]-ahf_data[n][2],2));
      
      if (delta < excision_rfct*ahf_data[n][4]) {
	
	D[ijk]  *= excision_fact;
	Sx[ijk] *= excision_fact;
	Sy[ijk] *= excision_fact;
	Sz[ijk] *= excision_fact;
	T[ijk]  *= excision_fact;
	p[ijk]  *= excision_fact;

      }

    } endfor_ijk;
      
  }
}


// set excision zone to *= r/(excision_rfct rah)
void grhd_set_excision_lin(tVarList *upre, tVarList *ucur)
{
  
  tL *level = ucur->level;    
  double *D  = vldataptr(ucur, MATTER.INDX_VAR_q    );
  double *T  = vldataptr(ucur, MATTER.INDX_VAR_q + 1);
  double *Sx = vldataptr(ucur, MATTER.INDX_VAR_q + 2);
  double *Sy = vldataptr(ucur, MATTER.INDX_VAR_q + 3);
  double *Sz = vldataptr(ucur, MATTER.INDX_VAR_q + 4);

  int n;
  double delta,factor;

  double *detg = Ptr(level, "adm_detg");
  double *p    = Ptr(level, "grhd_p");

  double *x    = Ptr(level, "x");
  double *y    = Ptr(level, "y");
  double *z    = Ptr(level, "z");
  

  for (n=0; n<maxAHF; n++) {
    
    if (PR) printf("   look for AH%d",n);
    
    if (!(get_ahf_val(n))) continue;
    
    if (PR) printf("  found   => excise region  (%e %e %e   r=%e)\n",ahf_data[n][0],ahf_data[n][1],ahf_data[n][2],ahf_data[n][4]);
    
    forallpoints_ijk(level) {
      
      // look if point is inside => set excision on evolution vars and pressure
      delta = sqrt(pow(x[ijk]-ahf_data[n][0],2)+pow(y[ijk]-ahf_data[n][1],2)+pow(z[ijk]-ahf_data[n][2],2));
      factor = delta/(excision_rfct*ahf_data[n][4]);
      
      if (delta < excision_rfct*ahf_data[n][4]) {
	
	D[ijk]  *= factor;
	Sx[ijk] *= factor;
	Sy[ijk] *= factor;
	Sz[ijk] *= factor;
	T[ijk]  *= factor;
	p[ijk]  *= factor;

      } 

    } endfor_ijk;
      
  }
}

































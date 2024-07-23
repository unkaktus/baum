/* lagrange_boxtobox_N.c */
/* Bernd Bruegmann 11/2008 */
/* mth 04/2012 -> only for prolong restrict */

#include "bam.h"
#include "interpolate.h"


#define PR 0
#define NMAX 20

double (*interpolate)(int N, double x, double xmin, double h, double *c, double *u);




/****************************************************************************/
/****************************************************************************/

/* split into two varlists, assume u is from the coarse level */
int split_varlist(tVarList *u, tL *lc, tL *lf,
                  tVarList **uc1, tVarList **uc2, tVarList **uc3, tVarList **uf1, tVarList **uf2, tVarList **uf3)
{
  int i,j;
  int to2,to3,to4;
  char *s;
  
  *uc1 = vlalloc(lc);
  *uc2 = vlalloc(lc);
  *uc3 = vlalloc(lc);
  *uf1 = vlalloc(lf);
  *uf2 = vlalloc(lf);
  *uf3 = vlalloc(lf);

  for (i=0; i<u->n; i++) {

    to2 = 0;
    to3 = 0;
    to4 = 0;

  if (Getv("physics", "matter")) {
    while (s = NextEntry(Gets("matter_interpolate_vars"))) {
      if (strncmp(VarName(u->index[i]),s,strlen(s))==0) 
        to2 = 1;	
    }

    if (Getv("conservative_amr", "yes")){   
      while (s = NextEntry(Gets("camr_interpolate_vars"))) {
        if (strncmp(VarName(u->index[i]),s,strlen(s))==0) 
          to4 = 1;
      }
    }
  }

   if (Getv("physics", "radiation")) {
     while (s = NextEntry(Gets("radiation_interpolate_vars"))) {
        if (strncmp(VarName(u->index[i]),s,strlen(s))==0) 
        to3 = 1;	
     }
   }

   if (to2) {
     vlpushone(*uc2, u->index[i]);
     vlpushone(*uf2, u->index[i]);
   } 
   else if (to3) {
     vlpushone(*uc3, u->index[i]);
     vlpushone(*uf3, u->index[i]);
   } else if (to4) {
      if(PR) printf("%s is a camr variable.\n",VarName(u->index[i]));
   } else {
      vlpushone(*uc1, u->index[i]);
      vlpushone(*uf1, u->index[i]);
   }
    
  }
  
  if (PR) {
    for (i=0; i<(((*uc1)->n > (*uc2)->n)?(*uc1)->n:(*uc2)->n); i++)
      printf("%d    %s  \t\t  %s  \t\t %s \n",i, i<(*uc1)->n?VarName((*uc1)->index[i]):"", i<(*uc2)->n?VarName((*uc2)->index[i]):"",i<(*uc3)->n?VarName((*uc3)->index[i]):"");
  }
  return ((*uc2)->n);
}




/* restrict list of variables from fine to coarse
   nbuffer  = 1:  restrict for Rflag >= 2       (multigrid)
   nbuffer >= 3:  restrict for Rflag >= 3       (evolution)
*/
void restrict_varlist(tG *g, int lfine, int lcoarse,
                      tVarList *uf, tVarList *uc, int nbuffer)
{
  
  if (PR) {
    int i;
    for (i=0; i<uf->n; i++) 
      printf("     restrict: %d  %s   (o=%d)   %d %d\n",i,VarName(uf->index[i]),Geti("order_RP"),lfine,lcoarse);
  }
  
  /* define stuff */
  tL *lf = g->level[lfine];
  tL *lc = g->level[lcoarse];

  /* check */
  if (lfine != lcoarse + 1)
    errorexit("restrict: watch your levels");
  if (nbuffer<3 && nbuffer !=1)
    errorexiti("nbuffer = %d, R needs 1 or >= 3", nbuffer);
  if (!uc || !uc->n || !uf || !uf->n) return;
  

  /* if ghost parent is used, make sure it is initialized */
  enableparent(lf, &lc);

  
  /* divide varlist into parts with different interpolation scheme */
  int order  = Geti("order_RP");
  int scheme = LAGRANGE;
  tVarList *uc1,*uf1,*uc2,*uf2,*uc3,*uf3;
  if (split_varlist(uc, lc,lf, &uf1,&uf2,&uf3,&uc1,&uc2,&uc3)) {
    /* lagrange part */
    restrict_varlist_boxtobox_N(lf, lc, uf1,uc1,nbuffer, order, scheme);
   
    /* matter part */
    order  = Geti("matter_interpolate_order");
    scheme = Getv("matter_interpolate_scheme_restriction","lagrange")*LAGRANGE +
             Getv("matter_interpolate_scheme_restriction","WENO")*WENO +
             Getv("matter_interpolate_scheme_restriction","WENOZ")*WENOZ +
             Getv("matter_interpolate_scheme_restriction","linear")*LIN;

    restrict_varlist_boxtobox_N(lf, lc, uf2,uc2,nbuffer, order, scheme);

    if (Getv("physics", "radiation")) {
    /* radiation part */
    order  = Geti("radiation_interpolate_order");
    scheme = Getv("radiation_interpolate_scheme_restriction","lagrange")*LAGRANGE +
             Getv("radiation_interpolate_scheme_restriction","WENO")*WENO +
             Getv("radiation_interpolate_scheme_restriction","WENOZ")*WENOZ +
             Getv("radiation_interpolate_scheme_restriction","linear")*LIN;

    restrict_varlist_boxtobox_N(lf, lc, uf3,uc3,nbuffer, order, scheme);
    }
    
  } else {
    /* full lagrange part */
    restrict_varlist_boxtobox_N(lf, lc, uf, uc, nbuffer, order, scheme);
  }
  
  
  /* free additional varlists */
  vlfree(uf1);
  vlfree(uc1);
  vlfree(uf2);
  vlfree(uc2);
  vlfree(uf3);
  vlfree(uc3);

  /* synchronize parent */
  if (nbuffer == 3)
    bampi_syncparent_send3(lf, uc);
  else if (nbuffer == 1) 
    bampi_syncparent_send23(lf, uc);


    bampi_vlsynchronize(uf);
    bampi_vlsynchronize(uc);
  
  /* some boundary conditions need to be applied here */
  set_boundary_symmetry(lc, uc);
  set_boundary_symmetry(lf, uf);

}


/* prolong list of variables from coarse to fine
     nbuffer  < 0:  prolong all, Pflag arbitrary (initialization)
     nbuffer  = 1:  prolong for Pflag  = 1       (multigrid)
     nbuffer >= 3:  prolong for Pflag >= 3       (evolution)
   note that P can avoid 3d ghostparent synchronization with recv0123
   if it is preceded by R, so if you call P separately, do recv0123 yourself!
*/
void prolong_varlist(tG *g, int lcoarse, int lfine, 
                     tVarList *uc, tVarList *uf, int nbuffer)
{
  if (PR) {
    int i;
    for (i=0; i<uf->n; i++) 
      printf("     prolong: %d  %s   (o=%d)   %d %d\n",i,VarName(uf->index[i]),Geti("order_RP"),lcoarse,lfine);
  }
  
  tL *lf = g->level[lfine];
  tL *lc = g->level[lcoarse];
  
  /* check */
  if (lfine != lcoarse + 1)	
    errorexit("prolong: watch your levels");
  if (!uc || !uc->n || !uf || !uf->n) return;

  /* if ghost parent is used, make sure it is initialized */
  enableparent(lf, &lc);

  
  /* synchronize parent */
  if (nbuffer >= 3)
    /* assumes that there is a preceding R so that recv0123 is not required */
    bampi_syncparent_recv012(lf, uc);
  else if (nbuffer == 1 || nbuffer == -1 || nbuffer == -3) 
    bampi_syncparent_recv0123(lf, uc);
  else
    errorexiti("nbuffer = %d, P needs >= 3, or = 1, -1, -3", nbuffer);
  if (nbuffer == -3) nbuffer = 3;
  

  /* divide varlist into parts with different interpolation scheme */
  int order = Geti("order_RP");
  int scheme  = LAGRANGE;
  tVarList *uc1,*uf1,*uc2,*uf2,*uc3,*uf3;  
  if (split_varlist(uc, lc,lf, &uf1,&uf2,&uf3,&uc1,&uc2,&uc3)) {
    /* lagrange part */
    prolong_varlist_boxtobox_N(lc, lf, uc1,uf1,nbuffer, order, scheme);
    
    /* matter part */
    order  = Geti("matter_interpolate_order");
    scheme = Getv("matter_interpolate_scheme_prolongation","lagrange")*LAGRANGE +
             Getv("matter_interpolate_scheme_prolongation","WENO")*WENO + 
             Getv("matter_interpolate_scheme_prolongation","WENOZ")*WENOZ + 
             Getv("matter_interpolate_scheme_prolongation","linear")*LIN;
    prolong_varlist_boxtobox_N(lc, lf, uc2,uf2,nbuffer, order, scheme);

    if (Getv("physics", "radiation")) {
    /* radiation part */
    order  = Geti("radiation_interpolate_order");
    scheme = Getv("radiation_interpolate_scheme_prolongation","lagrange")*LAGRANGE +
             Getv("radiation_interpolate_scheme_prolongation","WENO")*WENO + 
             Getv("radiation_interpolate_scheme_prolongation","WENOZ")*WENOZ + 
             Getv("radiation_interpolate_scheme_prolongation","linear")*LIN;
    prolong_varlist_boxtobox_N(lc, lf, uc3,uf3,nbuffer, order, scheme);
    }

  } else {
    /* full lagrange part */
    prolong_varlist_boxtobox_N(lc, lf, uc, uf, nbuffer, order, scheme);
  }
  
  
  /* free additional varlists */
  vlfree(uf1);
  vlfree(uc1);
  vlfree(uf2);
  vlfree(uc2);
  vlfree(uf3);
  vlfree(uc3);

  /* synchronize fine level */
  uf->level = lf;
  bampi_vlsynchronize(uf);

  /* some boundary conditions need to be applied here */
  set_boundary_symmetry(lf, uf);

}



/****************************************************************************/
/****************************************************************************/

/* Here we do N-th order interpolation, 
   see lagrange_boxtobox for 6-th order and for further explanations
*/ 


/* Prolongate box to box given all needed pre-adjusted information

   Example for grid points for 6th order:
       x   x   x   x   x   x   x   x   x   x   x  
                  o o o o o o o o o o  

   The key idea is to do three sweeps, one for each direction, thereby
   reusing 1-d interpolations that were already computed.
*/
void prolong_boxtobox_N_new(int N, int *cbox, int *fbox, int *fbox0,
                        double *cv, double *fv,
                        int cdi, int cdj, int cdk,
                        int fdi, int fdj, int fdk,
                        int scheme)
{
  if (PR) printf("start prolong_boxtobox_N\n");
  
  double *ev,*dv;
  //dv  = dmalloc(8*(cbox[1]-cbox[0])*(cbox[3]-cbox[2])*(cbox[5]-cbox[4]));
  //ev  = dmalloc(8*(cbox[1]-cbox[0])*(cbox[3]-cbox[2])*(cbox[5]-cbox[4]));
  
  
  if (scheme == LIN)
    interpolate = interpolate_linear_lim;
  else if (scheme == WENO)
    interpolate = interpolate_WENO_N;
  else if (scheme == WENOZ)
    interpolate = interpolate_WENOZ_N;
  else if (scheme == LAGRANGE) 
    interpolate = interpolate_lagrange_N;
  else 
    errorexit("interpolate_scheme is not implemented");
  
  
  //bampi_openmp_start
  
  int ndv = 0, nev = 0;
  int a = N;
  double sum;
  double c1[NMAX],c2[NMAX], u[NMAX];
  int d[NMAX+1];
  int i, j, k, p;
  int di, dj, dk;
  int dii, djj, dkk;
  int diii, djjj, dkkk;
  int i0, i1, j0, j1, k0, k1;
  int ii0, ii1, jj0, jj1, kk0, kk1;
  int iii0, iii1, jjj0, jjj1, kkk0, kkk1;
  int ijk, cijk,dijk,eijk;

  
  
  
  /* agrees with FiniteDifferences.nb, checked for N = 2,4,6,8,10,12 */
  coefficients_lagrange_N(N, N/2.0-0.25, 0.0, 1.0, c1);
  coefficients_lagrange_N(N, N/2.0-0.75, 0.0, 1.0, c2);
  
  if (0) {
    for (i = 0; i < N; i++) 
        printf("%.5e %.5e ", c1[i],c2[i]);
    printf("\n");
  }

  /* debug */
  if (0) {
    for (k = cbox[4]; k <= cbox[5]; k++) 
    for (j = cbox[2]; j <= cbox[3]; j++) 
    for (i = cbox[0]; i <= cbox[1]; i++) {
      cijk = i*cdi + j*cdj + k*cdk;
      cv[cijk] = j - 3.5;
      cv[cijk] = k - 3.5;
      cv[cijk] = i - 7.5;
      cv[cijk] = 1;
      if (1) printf("%2d %2d %2d cv[%4d] = %f\n", i, j, k, cijk, cv[cijk]); 
    }
  }

  
  if (PR) printf("0. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", 
      cbox[0],cbox[1],cbox[2],cbox[3],cbox[4],cbox[5], 
      fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);
  
  
  /***************************************************************************/
  /* z direction */
  i0 = 0;
  i1 = cbox[1] - cbox[0];
  j0 = 0;
  j1 = cbox[3] - cbox[2];
  k0 = 0;
  k1 = (cbox[5] - cbox[4] - a) * 2 + 1;

  di = 1;
  dj = (i1-i0+1);
  dk = (i1-i0+1)*(j1-j0+1);
  
  if (PR) printf("1. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", 
      i0,i1,j0,j1,k0,k1, fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);

  ndv = (i1-i0+1)*(j1-j0+1)*(k1-k0+1);
  dv  = dmalloc(ndv);
  /*
  if (bampi_openmp_rank()==0) {
    dv  = dmalloc(ndv);
  }
  bampi_openmp_barrier
  */

  // d1 = cdk; d2 = 2*cdk; d3 = 3*cdk; d4 = 4*cdk; d5 = 5*cdk; d6 = 6*cdk;
  for (p = 0; p <= N; p++) 
    d[p] = p*cdk;

  //bampi_openmp_loop
  for (j = j0; j <= j1; j++) {
  for (i = i0; i <= i1; i++) {

    /* map to the actual coarse box */
    ijk  = (i-i0)*di          + (j-j0)*dj          + (0)*dk;
    cijk = (i-i0+cbox[0])*cdi + (j-j0+cbox[2])*cdj + (cbox[4])*cdk;

    for (k = k0; k <= k1; k+=2) {
      
      if (0) printf("%f %f %f  %f %f %f  %f\n",
	       cv[cijk], cv[cijk+d[1]], cv[cijk+d[2]],  
	       cv[cijk+d[3]], cv[cijk+d[4]], cv[cijk+d[5]], cv[cijk+d[6]]);

      for (p = 0; p < N; p++) u[p] = cv[cijk + d[p]];
      dv[ijk] = interpolate(N,N/2.0-0.25, 0.0, 1.0, c1,u);
      ijk  += dk;
      
      for (p = 0; p < N; p++) u[p] = cv[cijk + d[p+1]];
      dv[ijk] = interpolate(N,N/2.0-0.75, 0.0, 1.0, c2,u);
      ijk  += dk;
      
      cijk += cdk;

      if (0) printf("%2d %2d %2d  %5d %5d  %6d %6d %6d  %f %f\n\n",
		    i, j, k, dk, cdk, ijk-2*dk, ijk-dk, cijk-cdk,
		    dv[ijk-2*dk], dv[ijk-dk]);
    }
  }
  }

  if (0) {
    for (k = k0; k <= k1; k++) 
    for (j = j0; j <= j1; j++) 
    for (i = i0; i <= i1; i++) {
      ijk = i*di + j*dj + k*dk;
      if (1) printf("%2d %2d %2d dv[%4d] = %f\n", i, j, k, ijk, dv[ijk]); 
    }
  }

  
  /***************************************************************************/
  /* y direction */
  ii0 = 0;
  ii1 = cbox[1] - cbox[0];
  jj0 = 0;
  jj1 = (cbox[3] - cbox[2] - a) * 2 + 1;
  kk0 = k0;
  kk1 = k1;
  
  dii = 1;
  djj = (ii1-ii0+1);
  dkk = (ii1-ii0+1)*(jj1-jj0+1);
  
  if (PR) printf("2. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", 
      ii0,ii1,jj0,jj1,kk0,kk1, fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);

  nev = (ii1-ii0+1)*(jj1-jj0+1)*(kk1-kk0+1);
  ev  = dmalloc(nev);
  /*
  if (bampi_openmp_rank()==0) {
    ev  = dmalloc(nev);
  }
  bampi_openmp_barrier
  */

  // d1 = dj; d2 = 2*dj; d3 = 3*dj; d4 = 4*dj; d5 = 5*dj; d6 = 6*dj;
  for (p = 0; p <= N; p++) 
    d[p] = p*dj;
  
  //bampi_openmp_loop
  for (k = kk0; k <= kk1; k++) {
  for (i = ii0; i <= ii1; i++) {

    ijk  = (i-ii0)*dii + (0)*djj + (k-kk0)*dkk;
    dijk = (i-ii0)*di  + (0)*dj  + (k-kk0)*dk;

    for (j = jj0; j <= jj1; j+=2) {
      
      if (0) printf("%f %f %f  %f %f %f  %f\n",
	       dv[dijk], dv[dijk+d[1]], dv[dijk+d[2]],  
	       dv[dijk+d[3]], dv[dijk+d[4]], dv[dijk+d[5]], dv[dijk+d[6]]);

      for (p = 0; p < N; p++) u[p] = dv[dijk + d[p]];
      ev[ijk] = interpolate(N,N/2.0-0.25, 0.0, 1.0, c1,u);
      ijk  += djj;
      
      for (p = 0; p < N; p++) u[p] = dv[dijk + d[p+1]];
      ev[ijk] = interpolate(N,N/2.0-0.75, 0.0, 1.0, c2,u);
      ijk  += djj;
      
      dijk += dj;

      if (0) printf("%2d %2d %2d  %5d %5d  %6d %6d %6d  %f %f\n\n",
		    i, j, k, djj, dj, ijk-2*djj, ijk-djj, cijk-dj,
		    ev[ijk-2*djj], ev[ijk-djj]);
    }
  }
  }

  if (0) {
    for (k = kk0; k <= kk1; k++) 
    for (j = jj0; j <= jj1; j++) 
    for (i = ii0; i <= ii1; i++) {
      ijk = i*dii + j*djj + k*dkk;
      if (1) printf("%2d %2d %2d dv[%4d] = %f\n", i, j, k, ijk, ev[ijk]); 
    }
  }

  
  /***************************************************************************/
  /* x direction */
  iii0 = 0;
  iii1 = (cbox[1] - cbox[0] - a) * 2 + 1;
  jjj0 = jj0;
  jjj1 = jj1;
  kkk0 = kk0;
  kkk1 = kk1;
  
  diii = 1;
  djjj = (iii1-iii0+1);
  dkkk = (iii1-iii0+1)*(jjj1-jjj0+1);
  
  if (PR) printf("3. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", 
      iii0,iii1,jjj0,jjj1,kkk0,kkk1, fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);
  
  // d1 = di; d2 = 2*di; d3 = 3*di; d4 = 4*di; d5 = 5*di; d6 = 6*di;
  for (p = 0; p <= N; p++) 
    d[p] = p*dii;

  //bampi_openmp_loop
  for (k = kkk0; k <= kkk1; k++) {
  for (j = jjj0; j <= jjj1; j++) {

    /* map into actual fine box */
    ijk  = (fbox[0])*diii + (j-jjj0+fbox[2])*djjj + (k-kkk0+fbox[4])*dkkk;
    eijk = (0)*dii        + (j-jjj0)*djj          + (k-kkk0)*dkk;

    for (i = iii0; i <= iii1; i+=2) {
      
      /* inside fbox0 we do not want do interpolate */
      if (i-iii0+fbox[0] >= fbox0[0] && i-iii0+fbox[0] <= fbox0[1] &&
          j-jjj0+fbox[2] >= fbox0[2] && j-jjj0+fbox[2] <= fbox0[3] &&
          k-kkk0+fbox[4] >= fbox0[4] && k-kkk0+fbox[4] <= fbox0[5]) {
	ijk  += 2*diii;
	eijk += dii;
	continue;
      }
      

      if (0) printf("%f %f %f  %f %f %f  %f\n",
	       ev[eijk], ev[eijk+d[1]], ev[eijk+d[2]],  
	       ev[eijk+d[3]], ev[eijk+d[4]], ev[eijk+d[5]], ev[eijk+d[6]]);

      for (p = 0; p < N; p++) u[p] = ev[eijk + d[p]];
      fv[ijk] = interpolate(N,N/2.0-0.25, 0.0, 1.0, c1,u);
      ijk  += diii;
      
      for (p = 0; p < N; p++) u[p] = ev[eijk + d[p+1]];
      fv[ijk] = interpolate(N,N/2.0-0.75, 0.0, 1.0, c2,u);
      ijk  += diii;
      
      eijk += dii;

      if (0) printf("%2d %2d %2d  %5d %5d  %6d %6d %6d  %f %f\n\n",
		    i, j, k, dii, di, ijk-2*dii, ijk-dii, cijk-di,
		    fv[ijk-2*dii], fv[ijk-dii]);
    }
  }
  }
  
  
  if (0) {
    for (k = kkk0; k <= kkk1; k++) 
    for (j = jjj0; j <= jjj1; j++) 
    for (i = iii0; i <= iii1; i++) {
      ijk = i*diii + j*djjj + k*dkkk;
      if (i-iii0+fbox[0] >= fbox0[0] && i-iii0+fbox[0] <= fbox0[1] &&
          j-jjj0+fbox[2] >= fbox0[2] && j-jjj0+fbox[2] <= fbox0[3] &&
          k-kkk0+fbox[4] >= fbox0[4] && k-kkk0+fbox[4] <= fbox0[5]) {
	if (1) printf("skipping %2d %2d %2d fv[%4d] = %f\n",
		      i, j, k, ijk, fv[ijk]); 
      } else {
	if (fv[ijk] != 1)
	  printf("%2d %2d %2d fv[%4d] = %f\n", i, j, k, ijk, fv[ijk]); 
      }
    }
    errorexit("testing");
  }
  
  //bampi_openmp_stop
  
  free(dv);
  free(ev);
  
  if (PR) printf("stop prolong_boxtobox_N\n");
} 

void prolong_boxtobox_N(int N, int *cbox, int *fbox, int *fbox0,
                        double *cv, double *fv,
                        int cdi, int cdj, int cdk,
                        int fdi, int fdj, int fdk,
                        int scheme)
{
  if (PR) printf("start prolong_boxtobox_N\n");
  
  /* global fields */
  double *dv  = dmalloc(8*(cbox[1]-cbox[0])*(cbox[3]-cbox[2])*(cbox[5]-cbox[4])*10);
  double *ev  = dmalloc(8*(cbox[1]-cbox[0])*(cbox[3]-cbox[2])*(cbox[5]-cbox[4])*10);
  
  /* set interpolation scheme */
  if (scheme == LIN)
    interpolate = interpolate_linear_lim;
  else if (scheme == WENO)
    interpolate = interpolate_WENO_N;
  else if (scheme == WENOZ)
    interpolate = interpolate_WENOZ_N;
  else if (scheme == LAGRANGE) 
    interpolate = interpolate_lagrange_N;
  else 
    errorexit("interpolate_scheme is not implemented");
  
  
  bampi_openmp_start
  
  //static int ndv = 0, nev = 0;
  int a = N;
  double sum;
  double c1[NMAX],c2[NMAX], u[NMAX];
  int d[NMAX+1];
  int i, j, k, p;
  int di, dj, dk;
  int ii, jj, kk;
  int dii, djj, dkk;
  int i0, i1, j0, j1, k0, k1;
  int ii0, ii1, jj0, jj1, kk0, kk1;
  int ijk, cijk;

  /* agrees with FiniteDifferences.nb, checked for N = 2,4,6,8,10,12 */
  coefficients_lagrange_N(N, N/2.0-0.25, 0.0, 1.0, c1);
  coefficients_lagrange_N(N, N/2.0-0.75, 0.0, 1.0, c2);
  
  if (0) {
    for (i = 0; i < N; i++) 
        printf("%.5e %.5e ", c1[i],c2[i]);
    printf("\n");
  }

  /* debug */
  if (0) {
    for (k = cbox[4]; k <= cbox[5]; k++) 
    for (j = cbox[2]; j <= cbox[3]; j++) 
    for (i = cbox[0]; i <= cbox[1]; i++) {
      cijk = i*cdi + j*cdj + k*cdk;
      cv[cijk] = j - 3.5;
      cv[cijk] = k - 3.5;
      cv[cijk] = i - 7.5;
      cv[cijk] = 1;
      if (1) printf("%2d %2d %2d cv[%4d] = %f\n", i, j, k, cijk, cv[cijk]); 
    }
  }

  
  /***************************************************************************/
  /* z direction */
  i0 = 0;
  i1 = cbox[1] - cbox[0];
  j0 = 0;
  j1 = cbox[3] - cbox[2];
  k0 = 0;
  k1 = (cbox[5] - cbox[4] - a) * 2 + 1;
  if (PR) printf("1. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", i0,i1,j0,j1,k0,k1, fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);


  i = (i1+1)*(j1+1)*(k1+1);
  //dv = dmalloc(ndv = i);

  di = 1;
  dj = i1+1;
  dk = (i1+1)*(j1+1);

  // d1 = cdk; d2 = 2*cdk; d3 = 3*cdk; d4 = 4*cdk; d5 = 5*cdk; d6 = 6*cdk;
  for (p = 0; p <= N; p++) 
    d[p] = p*cdk;
  
  bampi_openmp_loop
  for (j = j0; j <= j1; j++) {
  for (i = i0; i <= i1; i++) {

    ijk = i*di + j*dj;
    cijk = (i+cbox[0])*cdi + (j+cbox[2])*cdj + cbox[4]*cdk;

    for (k = k0; k <= k1; k+=2) {
      
      if (0) printf("%f %f %f  %f %f %f  %f\n",
	       cv[cijk], cv[cijk+d[1]], cv[cijk+d[2]],  
	       cv[cijk+d[3]], cv[cijk+d[4]], cv[cijk+d[5]], cv[cijk+d[6]]);

      for (p = 0; p < N; p++) u[p] = cv[cijk + d[p]];
      dv[ijk] = interpolate(N,N/2.0-0.25, 0.0, 1.0, c1,u);
      ijk += dk;
      
      for (p = 0; p < N; p++) u[p] = cv[cijk + d[p+1]];
      dv[ijk] = interpolate(N,N/2.0-0.75, 0.0, 1.0, c2,u);
      ijk += dk;
      cijk += cdk;

      if (0) printf("%2d %2d %2d  %5d %5d  %6d %6d %6d  %f %f\n\n",
		    i, j, k, dk, cdk, ijk-2*dk, ijk-dk, cijk-cdk,
		    dv[ijk-2*dk], dv[ijk-dk]);
    }
  }
  }

  if (0) {
    for (k = k0; k <= k1; k++) 
    for (j = j0; j <= j1; j++) 
    for (i = i0; i <= i1; i++) {
      ijk = i*di + j*dj + k*dk;
      if (1) printf("%2d %2d %2d dv[%4d] = %f\n", i, j, k, ijk, dv[ijk]); 
    }
  }

  
  /***************************************************************************/
  /* y direction */
  ii0 = 0;
  ii1 = cbox[1] - cbox[0];
  jj0 = 0;
  jj1 = (cbox[3] - cbox[2] - a) * 2 + 1;
  kk0 = 0;
  kk1 = k1;
  if (PR) printf("2. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", ii0,ii1,jj0,jj1,kk0,kk1, fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);

  i = (ii1+1)*(jj1+1)*(kk1+1);
  //ev = dmalloc(nev = i);
  
  dii = 1;
  djj = ii1+1;
  dkk = (ii1+1)*(jj1+1);

  // d1 = dj; d2 = 2*dj; d3 = 3*dj; d4 = 4*dj; d5 = 5*dj; d6 = 6*dj;
  for (p = 0; p <= N; p++) 
    d[p] = p*dj;
  
  bampi_openmp_loop
  for (k = kk0; k <= kk1; k++) {
  for (i = ii0; i <= ii1; i++) {

    ijk = i*dii + k*dkk;
    cijk = (i+i0)*di + j0*dj + (k+k0)*dk;

    for (j = jj0; j <= jj1; j+=2) {
      
      if (0) printf("%f %f %f  %f %f %f  %f\n",
	       dv[cijk], dv[cijk+d[1]], dv[cijk+d[2]],  
	       dv[cijk+d[3]], dv[cijk+d[4]], dv[cijk+d[5]], dv[cijk+d[6]]);

      for (p = 0; p < N; p++) u[p] = dv[cijk + d[p]];
      ev[ijk] = interpolate(N,N/2.0-0.25, 0.0, 1.0, c1,u);
      ijk += djj;
      
      for (p = 0; p < N; p++) u[p] = dv[cijk + d[p+1]];
      ev[ijk] = interpolate(N,N/2.0-0.75, 0.0, 1.0, c2,u);
      ijk += djj;
      cijk += dj;

      if (0) printf("%2d %2d %2d  %5d %5d  %6d %6d %6d  %f %f\n\n",
		    i, j, k, djj, dj, ijk-2*djj, ijk-djj, cijk-dj,
		    ev[ijk-2*djj], ev[ijk-djj]);
    }
  }
  }

  if (0) {
    for (k = kk0; k <= kk1; k++) 
    for (j = jj0; j <= jj1; j++) 
    for (i = ii0; i <= ii1; i++) {
      ijk = i*dii + j*djj + k*dkk;
      if (1) printf("%2d %2d %2d dv[%4d] = %f\n", i, j, k, ijk, ev[ijk]); 
    }
  }

  
  /***************************************************************************/
  /* x direction */
  ii0 = 0;
  ii1 = (cbox[1] - cbox[0] - a) * 2 + 1;
  i = (ii1+1)*(jj1+1)*(kk1+1);
  if (0) printf("%d %d %d %d %d %d\n", ii0, ii1, jj0, jj1, kk0, kk1);
  
  i0 = ii0;
  j0 = jj0;
  k0 = kk0;

  di = dii; 
  dj = djj; 
  dk = dkk;
  if (PR) printf("3. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", i0,i1,j0,j1,k0,k1, fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);

  /* map into actual fine box */
  ii0 = fbox[0];
  ii1 = fbox[1];
  jj0 = fbox[2];
  jj1 = fbox[3];
  kk0 = fbox[4];
  kk1 = fbox[5];
  dii = fdi;
  djj = fdj;
  dkk = fdk;

  // d1 = di; d2 = 2*di; d3 = 3*di; d4 = 4*di; d5 = 5*di; d6 = 6*di;
  for (p = 0; p <= N; p++) 
    d[p] = p*di;
  
  bampi_openmp_loop
  for (k = kk0; k <= kk1; k++) {
  for (j = jj0; j <= jj1; j++) {

    ijk  = ii0*dii + j*djj + k*dkk;
    cijk = (j+j0-jj0)*dj + (k+k0-kk0)*dk;

    for (i = ii0; i <= ii1; i+=2) {
      
      if (i >= fbox0[0] && i <= fbox0[1] &&
	  j >= fbox0[2] && j <= fbox0[3] &&
	  k >= fbox0[4] && k <= fbox0[5]) {
	ijk += 2*dii;
	cijk += di;
	continue;
      }

      if (0) printf("%f %f %f  %f %f %f  %f\n",
	       ev[cijk], ev[cijk+d[1]], ev[cijk+d[2]],  
	       ev[cijk+d[3]], ev[cijk+d[4]], ev[cijk+d[5]], ev[cijk+d[6]]);

      for (p = 0; p < N; p++) u[p] = ev[cijk + d[p]];
      fv[ijk] = interpolate(N,N/2.0-0.25, 0.0, 1.0, c1,u);
      ijk += dii;
      
      for (p = 0; p < N; p++) u[p] = ev[cijk + d[p+1]];
      fv[ijk] = interpolate(N,N/2.0-0.75, 0.0, 1.0, c2,u);
      ijk += dii;
      cijk += di;

      if (0) printf("%2d %2d %2d  %5d %5d  %6d %6d %6d  %f %f\n\n",
		    i, j, k, dii, di, ijk-2*dii, ijk-dii, cijk-di,
		    fv[ijk-2*dii], fv[ijk-dii]);
    }
  }
  }
  
  
  if (0) {
    for (k = kk0; k <= kk1; k++) 
    for (j = jj0; j <= jj1; j++) 
    for (i = ii0; i <= ii1; i++) {
      ijk = i*dii + j*djj + k*dkk;
      if (i >= fbox0[0] && i <= fbox0[1] &&
	  j >= fbox0[2] && j <= fbox0[3] &&
	  k >= fbox0[4] && k <= fbox0[5]) {
	if (1) printf("skipping %2d %2d %2d fv[%4d] = %f\n",
		      i, j, k, ijk, fv[ijk]); 
      } else {
	if (fv[ijk] != 1)
	  printf("%2d %2d %2d fv[%4d] = %f\n", i, j, k, ijk, fv[ijk]); 
      }
    }
    errorexit("testing");
  }
  bampi_openmp_stop
  
  free(dv);
  free(ev);
  
  if (PR) printf("stop prolong_boxtobox_N\n");
} 



/* prolongate list of variables from coarse to fine
   efficient box to box algorithm to replace cube to point amr algorithm
*/
void prolong_varlist_boxtobox_N(tL *lc, tL *lf, 
                                tVarList *uc, tVarList *uf, int nbuffer,
                                int order, int scheme)
{
  if (PR) printf("start  prolong_varlist_boxtobox_N    %d -> %d\n", lc->l,lf->l);
  
  int nghosts = Geti("bampi_nghosts");
  double *flagprolong = Ptr(lf, "flagprolong");
  double buffer = 1.0*Geti("amr_nbuffer");
  int b = order/2; 
  double fbox[6], fbox0[6], cbox[6];
  int ifbox[6], ifbox0[6], icbox[6];
  int ii, jj, kk;
  int i, n;

  /* for all boxes */
  forallboxes(lf) {

    /* figure out local bounding boxes */



  /*dtim: find region we have to prolong, 
    when levels are not aligned */ 
   double prolong_zone = buffer;
   if (!parentaligned(lf)) 
       { //prolong_zone = 4*(Geti("order_centered")-1); //CHANGEME this is not always correct!
         if ( prolong_zone > buffer)  prolong_zone = buffer;
       }

   // prolong_zone = buffer;


    /* initialize */
    for (i = 0; i < 6; i+=2) ifbox[i] = ifbox0[i] =  INT_MAX;
    for (i = 1; i < 6; i+=2) ifbox[i] = ifbox0[i] =  INT_MIN;

    /* for all points in this box */
    n = 0;
    forallpoints_boxijk(box) {
      
      /* ignore symmetry ghosts */
      if (boundaryflag(lf, ijk) == SYMBOUND) continue;

      /* actually, why bother with ghosts for nbuffer>=0?
         this is what we do for amr:
         if (boundaryflag(lf, ijk) == GHOBOUND && nbuffer < 0) continue;
      */
      if (boundaryflag(lf, ijk) == GHOBOUND) continue;

      /* find bounding box of indices */
   //   if (nbuffer >= 3 && dlessgreater(0,flagprolong[ijk],prolong_zone) ||  // 3 or more buffer nodes
      if ((nbuffer >= 3 && flagprolong[ijk]  > 0 && prolong_zone >= flagprolong[ijk]) ||  // 3 or more buffer nodes
 //       if (nbuffer >= 3 && flagprolong[ijk]  > 0 ||  // 3 or more buffer nodes
	  (nbuffer == 1 && flagprolong[ijk] == 1 ) ||  // 1 buffer node 
	  nbuffer < 0) {                            // everywhere
	/* these need prolongation */
	if (i < ifbox[0]) ifbox[0] = i;
	if (i > ifbox[1]) ifbox[1] = i;
	if (j < ifbox[2]) ifbox[2] = j;
	if (j > ifbox[3]) ifbox[3] = j;
	if (k < ifbox[4]) ifbox[4] = k;
	if (k > ifbox[5]) ifbox[5] = k;
	n++;
      } else {
	/* these do not need prolongation */
	if (i < ifbox0[0]) ifbox0[0] = i;
	if (i > ifbox0[1]) ifbox0[1] = i;
	if (j < ifbox0[2]) ifbox0[2] = j;
	if (j > ifbox0[3]) ifbox0[3] = j;
	if (k < ifbox0[4]) ifbox0[4] = k;
	if (k > ifbox0[5]) ifbox0[5] = k;
      }
    } endfor_ijk;
    
    /* if box not empty */
    if (n) {

      /* derive bounding boxes for coordinates */
      for (i = 0; i < 6; i++) {
	fbox[i]  = box->com->bbox[i-i%2] + box->dd[i/2] * ifbox[i];
	fbox0[i] = box->com->bbox[i-i%2] + box->dd[i/2] * ifbox0[i];
      }
      
      /* now derive bounding boxes on parent level
	 this assumes all the necessary points are available by construction
      */
      box_point_indices(box->pr, box, &ii, &jj, &kk,
			ifbox[0], ifbox[2], ifbox[4]);
      icbox[0] = ii-b+1;
      icbox[2] = jj-b+1;
      icbox[4] = kk-b+1;
      box_point_indices(box->pr, box, &ii, &jj, &kk,
			ifbox[1], ifbox[3], ifbox[5]);
      icbox[1] = ii+b;
      icbox[3] = jj+b;
      icbox[5] = kk+b;
      for (i = 0; i < 6; i++)
	cbox[i] = box->pr->com->bbox[i-i%2] + box->pr->dd[i/2] * icbox[i];
      
      /* info */
      i = box_ainb(cbox, box->pr->com->bbox);
      if (PR) {
	if (!i) {
	  printf("prbox "); 
	  printbbox(lc, box->pr->bbox, box->pr->ibbox);
	}
	printf("cbox  "); printbbox(lc, cbox , icbox ); 
	printf("fbox  "); printbbox(lf, fbox,  ifbox ); 
	printf("fbox0 "); printbbox(lf, fbox0, ifbox0); 
        printf("      interpolation_order=%d\n", order);
      }
      //exit(0);
    }

    /* catch box size errors */
    if (n && !i)
      errorexit("prolongation: boxes are too small, enlarge nxyz");

    /* if there are no points to interpolate into, continue 
       (has to come after bampi calls) */
    if (!n) continue;


    /* now we have the bbox information that we need */

    /* do variables one by one for cache locality, but we should
       group this if it is more efficient to do so */
    for (i = 0; i < uc->n; i++) {
      if (PR) printf("calling p_btb for %d->%d (%s->%s)  %d %d\n", 
		    uc->index[i], uf->index[i],
		    VarName(uc->index[i]), VarName(uf->index[i]),
		    box->pr->noffset, box->noffset);
      
      
      prolong_boxtobox_N(order, icbox, ifbox, ifbox0, 
                        lf->prlocal->v[uc->index[i]] + box->pr->noffset,
                        lf->v[uf->index[i]] + box->noffset,
                        box->pr->di, box->pr->dj, box->pr->dk,
                        box->di, box->dj, box->dk,
                        scheme);
    }

    
  } endforboxes;

  if (PR) printf("stop   prolong_varlist_boxtobox_N\n");
}







/****************************************************************************/
/****************************************************************************/

/* same box to box strategy for restriction
   note that restriction produces results on only roughly 1/8 the number 
   of points that prolongation has to fill in, and some such factor was
   actually visible in profiling; still, after improving prolongation
   the amr restriction was a factor of two slower than btb prolongation,
   which is why we give a fast restriction algorithm, too
*/

/* Restrict box to box given all needed pre-adjusted information
     assumes standard cell-based staggering
     assumes complete child cells
     assumes that all needed fine and coarse points are available
     (i.e. no stencil shifting)
   Example for grid points for 6th order:
       x   x   x   x   x   x   x   x   x   x   x  
                  o o o o o o o o o o  

   The key idea is to do three sweeps, one for each direction, thereby
   reusing 1-d interpolations that were already computed.
*/

/*
    mth: I copied prolong_boxtobox_N and changed some lines
    
    the idea was to change coarse and fine box
    and change resolutions and number of points
*/
void restrict_boxtobox_N(int N, int *fbox, int *cbox,
                        double *fv, double *cv,
                        int fdi, int fdj, int fdk,
                        int cdi, int cdj, int cdk,
                        int scheme)
{
  if (PR) printf("start restrict_boxtobox_N\n");
  
  /* global fields */
  double *dv  = dmalloc((fbox[1]-fbox[0])*(fbox[3]-fbox[2])*(fbox[5]-fbox[4])*10);
  double *ev  = dmalloc((fbox[1]-fbox[0])*(fbox[3]-fbox[2])*(fbox[5]-fbox[4])*10);
  
  /* set interpolation scheme */
  if (scheme == LIN)
    interpolate = interpolate_linear_avg;
  else if (scheme == WENO)
    interpolate = interpolate_WENO_N;
  else if (scheme == WENOZ)
    interpolate = interpolate_WENOZ_N;
  else if (scheme == LAGRANGE) 
    interpolate = interpolate_lagrange_N;
  else 
    errorexit("interpolate_scheme is not implemented");
  
  
  bampi_openmp_start
  
  int a = N;
  double sum;
  double c[NMAX], u[NMAX];
  int d[NMAX+1];
  int i, j, k, p;
  int di, dj, dk;
  int ii, jj, kk;
  int dii, djj, dkk;
  int i0, i1, j0, j1, k0, k1;
  int ii0, ii1, jj0, jj1, kk0, kk1;
  int ijk, fijk;

  /* agrees with FiniteDifferences.nb, checked for N = 2,4,6,8,10,12 */
  coefficients_lagrange_N(N, N/2.0-0.5, 0.0, 1.0, c);
  
  
  /***************************************************************************/
  /* z direction */
  i0 = 0 ;
  i1 = fbox[1] - fbox[0];
  j0 = 0;
  j1 = fbox[3] - fbox[2];
  k0 = 0;
  k1 = (fbox[5] - fbox[4] - a+2 )/2;
  if (PR) printf("1. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", i0,i1,j0,j1,k0,k1, cbox[0],cbox[1],cbox[2],cbox[3],cbox[4],cbox[5]);

  i = (i1+1)*(j1+1)*(k1+1);
  //ev = dmalloc(ndv = i);

  di = 1;
  dj = i1+1;
  dk = (i1+1)*(j1+1);

  for (p = 0; p <= N; p++) 
    d[p] = p*fdk;

  bampi_openmp_loop
  for (j = j0; j <= j1; j++) {
  for (i = i0; i <= i1; i++) {

    ijk  = i*di + j*dj;
    fijk = (i+fbox[0])*fdi + (j+fbox[2])*fdj + fbox[4]*fdk;

    for (k = k0; k <= k1; k+=1) {
      for (p = 0; p < N; p++) u[p] = fv[fijk + d[p]];
      ev[ijk] = interpolate(N,N/2.0-0.5, 0.0, 1.0, c,u);
      
      ijk += dk;
      fijk += 2*fdk;
    }
  }
  }

  
  /***************************************************************************/
  /* y direction */
  ii0 = 0;
  ii1 = fbox[1] - fbox[0];
  jj0 = 0;
  jj1 = (fbox[3] - fbox[2] - a+2 )/2;
  kk0 = 0;
  kk1 = k1;
  if (PR) printf("2. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", ii0,ii1,jj0,jj1,kk0,kk1, cbox[0],cbox[1],cbox[2],cbox[3],cbox[4],cbox[5]);

  i = (ii1+1)*(jj1+1)*(kk1+1);
  //dv = dmalloc(nev = i);

  dii = 1;
  djj = ii1+1;
  dkk = (ii1+1)*(jj1+1);
  
  for (p = 0; p <= N; p++) 
    d[p] = p*dj;
  
  bampi_openmp_loop
  for (k = kk0; k <= kk1; k++) {
  for (i = ii0; i <= ii1; i++) {

    ijk = i*dii + k*dkk;
    fijk = (i+i0)*di + j0*dj + (k+k0)*dk;

    for (j = jj0; j <= jj1; j+=1) {
      for (p = 0; p < N; p++) u[p] = ev[fijk + d[p]];
      dv[ijk] = interpolate(N,N/2.0-0.5, 0.0, 1.0, c,u);
      
      ijk += djj;
      fijk += 2*dj;
    }
  }
  }
  

  /***************************************************************************/
  /* x direction */
  i0 = 0;
  i1 = (fbox[1] - fbox[0] - a+2)/2;
  j0 = 0;
  j1 = jj1;
  k0 = 0;
  k1 = kk1;
  if (PR) printf("3. [%d,%d]x[%d,%d]x[%d,%d]  ->  [%d,%d]x[%d,%d]x[%d,%d]\n", i0,i1,j0,j1,k0,k1, cbox[0],cbox[1],cbox[2],cbox[3],cbox[4],cbox[5]);
  
  di = cdi;
  dj = cdj;
  dk = cdk;
  
  for (p = 0; p <= N; p++) 
    d[p] = p*dii;
  
  bampi_openmp_loop
  for (k = k0; k <= k1; k++) {
  for (j = j0; j <= j1; j++) {

    ijk = cbox[0]*di + (cbox[2]+j)*dj + (cbox[4]+k)*dk;
    fijk = ii0*dii + (j+jj0)*djj + (k+kk0)*dkk;
  
      for (i = i0; i <= i1; i+=1) {
        for (p = 0; p < N; p++) u[p] = dv[fijk + d[p]];
        cv[ijk] = interpolate(N,N/2.0-0.5, 0.0, 1.0, c,u);
        
        ijk += di;
        fijk += 2*dii;

      }
  }
  }
  
  bampi_openmp_stop
  
  
  free(dv);
  free(ev);
  
  if (PR) printf("stop restrict_boxtobox_N\n");
} 




/* restrict list of variables from fine to coarse
   efficient box to box algorithm to replace cube to point amr algorithm
*/
void restrict_varlist_boxtobox_N(tL *lf, tL *lc, 
                                 tVarList *uf, tVarList *uc, int nbuffer,
                                 int order, int scheme)
{
  if (PR) printf("start restrict_varlist_boxtobox_N    %d -> %d\n", lf->l,lc->l);
  
  int nghosts = Geti("bampi_nghosts");
  double *flagrestrict = Ptr(lf->prlocal, "flagrestrict");
  
  int Rflag;
  int b = order/2;
  double fbox[6], cbox[6];
  int ifbox[6], icbox[6];
  int fnbox;
  int ii, jj, kk;
  int i, n;

  
  /* set minimal Rflag based on nbuffer */
  if (nbuffer >= 3) Rflag = 3;
  else if (nbuffer == 1) Rflag = 2;
  else errorexiti("nbuffer = %d, R needs 1 or >= 3", nbuffer);
  
  /* for all boxes
     note that we loop over fine boxes since they know their parent, while
     if we looped over coarse boxes, we would have to figure out which
     fine box is involved
  */
  forallboxes(lf) {
    if (PR) printf("restrict_varlist_boxtobox_N in box  %d\n",nbox);
    /* figure out local bounding boxes */

    /* initialize */
    for (i = 0; i < 6; i+=2) icbox[i] = INT_MAX;
    for (i = 1; i < 6; i+=2) icbox[i] = INT_MIN;

    /* for all points in this box */
    n = 0;
    forallpoints_boxijk(box->pr) {

      /* these need restriction */
      if (flagrestrict[ijk] >= Rflag) {
	if (i < icbox[0]) icbox[0] = i;
	if (i > icbox[1]) icbox[1] = i;
	if (j < icbox[2]) icbox[2] = j;
	if (j > icbox[3]) icbox[3] = j;
	if (k < icbox[4]) icbox[4] = k;
	if (k > icbox[5]) icbox[5] = k;
	n++;
      } 
    } endfor_ijk;
   
    
    /* if box not empty */
    if (n) {

      /* derive bounding boxes for coordinates */
      for (i = 0; i < 6; i++) 
	cbox[i]  = box->pr->com->bbox[i-i%2] + box->pr->dd[i/2] * icbox[i];
      
    
      /* derive bounding boxes on fine level
	 this assumes all the necessary points are available by construction
      */
      box_point_indices(box, box->pr, &ii, &jj, &kk,
			icbox[0], icbox[2], icbox[4]);
      ifbox[0] = ii-b+1;
      ifbox[2] = jj-b+1;
      ifbox[4] = kk-b+1;
      box_point_indices(box, box->pr, &ii, &jj, &kk,
			icbox[1], icbox[3], icbox[5]);
      ifbox[1] = ii+b;
      ifbox[3] = jj+b;
      ifbox[5] = kk+b;
      for (i = 0; i < 6; i++)
        fbox[i] = box->com->bbox[i-i%2] + box->dd[i/2] * ifbox[i];

      /* info */
      i = box_ainb(fbox, box->com->bbox);
      if (PR) {
	printf("cbox  "); printbbox(lc, cbox , icbox ); 
	printf("fbox  "); printbbox(lf, fbox,  ifbox ); 
        printf("      interpolation_order=%d\n", order);
      }
    }

    /* catch box size errors */
    if (n && !i) {
      printf("%d %d\n",n,i);
      errorexit("restriction: boxes are too small, enlarge nxyz");
    }
    /* if there are no points to interpolate into, continue 
       (has to come after bampi calls) */
    if (!n) continue;


    /* now we have the bbox information that we need */

    /* do variables one by one for cache locality, but we should
       group this if it is more efficient to do so */
    for (i = 0; i < uc->n; i++) {
      if (PR) printf("calling p_btb for %d->%d (%s->%s)  %d %d\n", 
		    uc->index[i], uf->index[i],
		    VarName(uc->index[i]), VarName(uf->index[i]),
		    box->pr->noffset, box->noffset);
      
      
      restrict_boxtobox_N(order, ifbox, icbox,
                         lf->v[uf->index[i]] + box->noffset,
                         lf->prlocal->v[uc->index[i]] + box->pr->noffset,
                         box->di, box->dj, box->dk,
                         box->pr->di, box->pr->dj, box->pr->dk,
                         scheme);
    }
    
  } endforboxes;
  
  if (PR) printf("stop  restrict_varlist_boxtobox_N\n");
}





























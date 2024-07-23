/* c_amr.c 
   dtim, sbernuz 10/13*/

#include "bam.h"
#include "matter.h"

#define PR 0
#define DEBUG 0

/* Provide a mask of cells A and B for conservative amr, see 
   grqc 1112.3094. 
   Type A cells-parents, type B cells childs.
   
   boundary , entry in mask
   -x       ,      1        
   +x       ,      2
   -y       ,      3 
   +y       ,      4
   -z       ,      5
   +z       ,      6

   additionally, if are in a region, covered by a finer grid, we set cam_mask_A = 0.7, 
   then we can use it as a mask for the integral computation of e.g. D 
*/

int c_amr_mask(tL* level) {

   /* check if camr is used*/
   if (MATTER.CAMRACTIVE == 0) return 0;
   
   /* For some levels we make nothing!*/
   if (level-> l == 0)                        return 0;
   if (level-> l < Geti("camr_level_min"))    return 0;

   printf("recompute camr mask \n"); 

   /* define stuff */
   int lcur = level->l;
   int lcoarse = lcur-1;
   int i,n, nbuffer = Geti("amr_nbuffer");
   double fbox[6];
   double xminc,xmaxc,yminc,ymaxc,zminc,zmaxc;
   double xmminc,xmmaxc,ymminc,ymmaxc,zmminc,zmmaxc;
   double ldd = level->dx;
   double  dd = 0.5*ldd;
   tG *g = level->grid;
   tL *lc = g->level[lcoarse];

   /* call mask_A before enableparent and after enable parent
   in matter_fluxes we need the mask_Ac, in correct_varlist the mask_A, 
   FIXME: There must be a better way.*/
   
   double *mask_Ac = Ptr(lc,"camr_mask_A"); /* not from enableparent*/ 
   double *x       = Ptr(lc,"x");
   double *y       = Ptr(lc,"y");
   double *z       = Ptr(lc,"z");

   enableparent(level, &lc);
  
   double *mask_A  = Ptr(lc,"camr_mask_A");
   double *mask_B  = Ptr(level,"camr_mask_B");
   double *xf      = Ptr(level,"x");
   double *yf      = Ptr(level,"y");
   double *zf      = Ptr(level,"z");
   double *xc      = Ptr(lc,"x");
   double *yc      = Ptr(lc,"y");
   double *zc      = Ptr(lc,"z");
   double *rflag   = Ptr(lc,"flagrestrict"); 


   /* set A-mask to zero!*/
   forallpoints_ijk(g->level[lcoarse]) {
      mask_Ac[ijk] = 0.7;
   } endfor_ijk; 

   /* go through boxes of child-levels*/
   forallboxes(level) {

      /*save the boundaries*/
      for (n = 0; n < 6; n++) {
         fbox[n]   = box->bbox[n]   + (nbuffer-1)*ldd;
         fbox[n+1] = box->bbox[n+1] - (nbuffer-1)*ldd;
         n+= 1;
      }
   
      /* Where are A-cells*/ 
      xminc = fbox[0]-dd; xmaxc = fbox[1]+dd;
      yminc = fbox[2]-dd; ymaxc = fbox[3]+dd;
      zminc = fbox[4]-dd; zmaxc = fbox[5]+dd;

      /* the boundaries for the mask- computing integralouter*/ 
      xmminc = xminc - nbuffer*ldd; xmmaxc = xmaxc + nbuffer*ldd;
      ymminc = yminc - nbuffer*ldd; ymmaxc = ymaxc + nbuffer*ldd;
      zmminc = zminc - nbuffer*ldd; zmmaxc = zmaxc + nbuffer*ldd;

      if(PR) printf("xminc %e %e %e %e %e %e",xminc, xmaxc, yminc, ymaxc, zminc, zmaxc);

      /* Set type B cells */
      forallpoints_boxijk(box) {
         mask_B[ijk] =  0;
         if (dequal(xf[ijk],fbox[0])) 
            { if (dlessgreater(fbox[2],yf[ijk],fbox[3]) && 
                  dlessgreater(fbox[4],zf[ijk],fbox[5])) mask_B[ijk] = 1; }
         if (dequal(xf[ijk],fbox[1])) 
            { if (dlessgreater(fbox[2],yf[ijk],fbox[3]) && 
                  dlessgreater(fbox[4],zf[ijk],fbox[5])) mask_B[ijk] = 2; }
         if (dequal(yf[ijk],fbox[2]))
            { if (dlessgreater(fbox[0],xf[ijk],fbox[1]) && 
                  dlessgreater(fbox[4],zf[ijk],fbox[5])) mask_B[ijk] = 3; }
         if (dequal(yf[ijk],fbox[3])) 
            { if (dlessgreater(fbox[0],xf[ijk],fbox[1]) && 
                  dlessgreater(fbox[4],zf[ijk],fbox[5])) mask_B[ijk] = 4; }
         if (dequal(zf[ijk],fbox[4])) 
            { if (dlessgreater(fbox[0],xf[ijk],fbox[1]) && 
                  dlessgreater(fbox[2],yf[ijk],fbox[3])) mask_B[ijk] = 5; }
         if (dequal(zf[ijk],fbox[5])) 
            { if (dlessgreater(fbox[0],xf[ijk],fbox[1]) && 
                  dlessgreater(fbox[2],yf[ijk],fbox[3])) mask_B[ijk] = 6; }
      } endfor_ijk;

      /* IMPROVE, for binaries we do this two times, this is not good!!!*/
      /* Actually, is this part needed at all??? */ 
  
      forallpoints_ijk(g->level[lcoarse]) {
         if (dequal(x[ijk],xminc))  
            { if (dlessgreater(yminc,y[ijk],ymaxc) && 
                  dlessgreater(zminc,z[ijk],zmaxc)) mask_Ac[ijk] =  1; }
         if (dequal(x[ijk],xmaxc))           
            { if (dlessgreater(yminc,y[ijk],ymaxc) && 
                  dlessgreater(zminc,z[ijk],zmaxc)) mask_Ac[ijk] =  2; }
         if (dequal(y[ijk],yminc))  
            { if (dlessgreater(xminc,x[ijk],xmaxc) && 
                  dlessgreater(zminc,z[ijk],zmaxc)) mask_Ac[ijk] =  3; }
         if (dequal(y[ijk],ymaxc))  
            { if (dlessgreater(xminc,x[ijk],xmaxc) && 
                  dlessgreater(zminc,z[ijk],zmaxc)) mask_Ac[ijk] =  4; }
         if (dequal(z[ijk],zminc))  
            { if (dlessgreater(xminc,x[ijk],xmaxc) && 
                  dlessgreater(yminc,y[ijk],ymaxc)) mask_Ac[ijk] =  5; }
         if (dequal(z[ijk],zmaxc))  
            { if (dlessgreater(xminc,x[ijk],xmaxc) && 
                  dlessgreater(yminc,y[ijk],ymaxc)) mask_Ac[ijk] =  6; }

         /* this one is only for the output, where we use the integral over different level, 
            hopefully this will help to messure the mass better! */
         if (mask_Ac[ijk] == 0.7 )
            if (!((dlessgreater(xmminc,x[ijk],xmmaxc)) && 
                  (dlessgreater(ymminc,y[ijk],ymmaxc)) && 
                  (dlessgreater(zmminc,z[ijk],zmmaxc))))  mask_Ac[ijk] = 0.;

    /*     if (mask_Ac[ijk] == 0.7 )
            if (!((dlessgreater(xmminc,x[ijk],xmmaxc)) && 
                  (dlessgreater(ymminc,y[ijk],ymmaxc)) && 
                  (dlessgreater(zmminc,z[ijk],zmmaxc))))  mask_Ac[ijk] = 0.25;
    */
         if((PR)&& (mask_Ac[ijk]!=0))  printf(" real %e %e %e %e \n",x[ijk],y[ijk],z[ijk], mask_Ac[ijk]);
            
      } endfor_ijk;


      /* Set type A cells for parent of the actual mpi-job*/
      /* IMPROVE, again half of the computation to much for binaries!*/

      forallpoints_ijk(lc){
         mask_A[ijk] = 0.0;
         if(dequal(rflag[ijk],2.0)){ /* This if ensures that we stay in the correct box. */       
            if (dequal(xc[ijk],xminc))  
               { if (dlessgreater(yminc,yc[ijk],ymaxc) && 
                     dlessgreater(zminc,zc[ijk],zmaxc)) mask_A[ijk] =  1; }
            if (dequal(xc[ijk],xmaxc))          
               { if (dlessgreater(yminc,yc[ijk],ymaxc) && 
                     dlessgreater(zminc,zc[ijk],zmaxc)) mask_A[ijk] =  2; }
            if (dequal(yc[ijk],yminc))  
               { if (dlessgreater(xminc,xc[ijk],xmaxc) && 
                     dlessgreater(zminc,zc[ijk],zmaxc)) mask_A[ijk] =  3; }
            if (dequal(yc[ijk],ymaxc))  
               { if (dlessgreater(xminc,xc[ijk],xmaxc) && 
                     dlessgreater(zminc,zc[ijk],zmaxc)) mask_A[ijk] =  4; }
            if (dequal(zc[ijk],zminc))  
               { if (dlessgreater(xminc,xc[ijk],xmaxc) && 
                     dlessgreater(yminc,yc[ijk],ymaxc)) mask_A[ijk] =  5; }
            if (dequal(zc[ijk],zmaxc))  
               { if (dlessgreater(xminc,xc[ijk],xmaxc) && 
                     dlessgreater(yminc,yc[ijk],ymaxc)) mask_A[ijk] =  6; }           
         }
         if((PR)&&(mask_A[ijk]!=0)) printf("parent %e %e %e %e %e\n",xc[ijk],yc[ijk],zc[ijk], mask_A[ijk],rflag[ijk]);
              
      } endfor_ijk;
      
   } endforboxes; 


   /* make sure everything is set up correct and we do not have any flux stored!*/
   if(level->time > 0) {

      double *f, *c;
      int nvar;

      for (nvar = 0; nvar < MATTER.NVq; nvar++) {

         f     = level->v[Ind("camr_D")+nvar];
         c     = lc->v[Ind("camr_D")+nvar];

         forallpoints_ijk(level) {
               f[ijk]      = 0; 
         } endfor_ijk;


         /*add correction and reset A-cells*/
         forallpoints_ijk(lc) {
               c[ijk]      = 0; 
         } endfor_ijk;
      }
   }
   return 1;
}

/* add the correction stored in camr_XXX to the grhd_XXX variable.
   adn reset cam_XXX-cells */

void correct_varlist(tG *g, int lfine, int lcoarse, tVarList *uf, tVarList *uc , int nbuffer)
{
   
   if(PR) printf("I am in correct_varlist!\n");

   /* define stuff */
   tL *lf = g->level[lfine];
   tL *lc = g->level[lcoarse];


   /* in all these cases we do nothing!*/ 
   if (lf->l > g->lmax) return ; 
   if (lf->l == 0     ) return ;
   if (!parentaligned(g->level[lfine])) return;
   if (dequal(lf->time,0)) return; 

   /* this works, but it is not optimal, CHANGE ME when I am working */ 
   enableparent(lf, &lc);
   bampi_syncparent_recv0123(lf, uc); /*Do we need 0123??? */

   if(PR) printf("I do not skip level \n");
   if(PR) printf("lmax %d lf %d lc %d \n", g->lmax,lf->l, lc->l);

   double *mask_A  = Ptr(lc, "camr_mask_A");
   double *mask_B  = Ptr(lf, "camr_mask_B");
   double *xc      = Ptr(lc, "x");
   double *yc      = Ptr(lc, "y");
   double *zc      = Ptr(lc, "z");
   double *xf      = Ptr(lf, "x");
   double *yf      = Ptr(lf, "y");
   double *zf      = Ptr(lf, "z");
   double fc, fcorrection, camrtresh = Getd("camr_treshold"), camr_cut= Getd("camr_atm_cut");
   double dl = lf->dx;   
   double dlc = 2.*dl;
   int i,j,k,n,nvar,direction,dir; 
   double *f, *c, *cgrhd, *c_rho;
   int fijk, fdi, fdj, fdk;
   tB *fbox;  


   /* go through all camr_variables */	

   for (nvar = 0; nvar < MATTER.NVq; nvar++) {

      if (PR) printf("nvar %d \n",nvar);

      f     = lf->v[Ind(MATTER.c_names_list[nvar])];
      c     = lc->v[Ind(MATTER.c_names_list[nvar])];
      cgrhd = lc->v[Ind(MATTER.q_names_list[nvar])];
      c_rho = lc->v[Ind("grhd_D")];

      forallboxes(lc) {
         forallpoints_boxijk(box) {

            if (mask_A[ijk]) {    
               if(PR) printf("%d %d %d %f \n" , i , j , k, mask_A[ijk]);
               if(PR) printf("%le %le %le %f \n" , xc[ijk] , yc[ijk] , zc[ijk], mask_A[ijk]);

               fbox = box_containing_ijk(lf, box, i, j, k);
               fdi  = fbox->di;
               fdj  = fbox->dj;
               fdk  = fbox->dk;
               fijk = box_child_ijk(fbox, box, i, j, k);


 /*       if (0){
               printf("%le %le %le %le %le\n" , xc[ijk] , yc[ijk] , zc[ijk], mask_A[ijk],c[ijk]);
               printf("%le %le %le %le %le\n" , xf[fijk] , yf[fijk] , zf[fijk], mask_B[fijk], f[fijk]);
               printf("%le %le %le %le %le\n" , xf[fijk+fdi+fdj+fdk] , yf[fijk+fdi+fdj+fdk] , zf[fijk+fdi+fdj+fdk], mask_B[fijk+fdi+fdj+fdk],f[fijk+fdi+fdj+fdk]);
               printf("%le %le %le %le %le\n" , xf[fijk+fdi] , yf[fijk+fdi] , zf[fijk+fdi], mask_B[fijk+fdi],f[fijk+fdi]);
               printf("%le %le %le %le %le\n" , xf[fijk+fdj] , yf[fijk+fdj] , zf[fijk+fdj], mask_B[fijk+fdj],f[fijk+fdj]);
               printf("%le %le %le %le %le\n" , xf[fijk+fdk] , yf[fijk+fdk] , zf[fijk+fdk], mask_B[fijk+fdk],f[fijk+fdk]);
               printf("%le %le %le %le %le\n" , xf[fijk+fdi+fdj] , yf[fijk+fdi+fdj] , zf[fijk+fdi+fdj], mask_B[fijk+fdi+fdj],f[fijk+fdi+fdj]);
               printf("%le %le %le %le %le\n" , xf[fijk+fdi+fdk] , yf[fijk+fdi+fdk] , zf[fijk+fdi+fdk], mask_B[fijk+fdi+fdk],f[fijk+fdi+fdk]);
               printf("%le %le %le %le %le\n \n" , xf[fijk+fdj+fdk] , yf[fijk+fdj+fdk] , zf[fijk+fdj+fdk], mask_B[fijk+fdj+fdk],f[fijk+fdj+fdk]);
               }
*/
               fcorrection  = 0.125*
                             (f[fijk] + f[fijk+fdi+fdj+fdk] +
                              f[fijk+fdi] + f[fijk+fdj] + f[fijk+fdk] +
                              f[fijk+fdi+fdj] + f[fijk+fdi+fdk] + f[fijk+fdj+fdk]);
    
               /* safety checks*/
               if (c[ijk]!=0. &&  fcorrection!=0. &&
                     dlessgreater(1./camrtresh,((fabs(fcorrection)+1.e-15)/(fabs(c[ijk])+1.e-15)),camrtresh))  
                  c[ijk] += fcorrection; 
               else { 
                  c[ijk]  = 0 ; 
                  // printf("Problem with camr: set correction to zero! at (%e,%e,%e)\n",xc[ijk],yc[ijk],zc[ijk]);
               }

               if(c_rho[ijk] < camr_cut) { 
                  if(PR) printf("Do not correct for small densities.\n");
                  c[ijk] = 0; 
               }


               /* for debugging, see if parallelization is correct*/
               if(DEBUG)  {   
            
                  fc = 0;
               
                  if(!dequal(mask_B[fijk],0))              fc+=1;      
                  if(!dequal(mask_B[fijk+fdi+fdj+fdk],0))  fc+=1;          
                  if(!dequal(mask_B[fijk+fdi+fdj],0))      fc+=1;      
                  if(!dequal(mask_B[fijk+fdi+fdk],0))      fc+=1;       
                  if(!dequal(mask_B[fijk+fdk+fdj],0))      fc+=1;      
                  if(!dequal(mask_B[fijk+fdi],0))          fc+=1;      
                  if(!dequal(mask_B[fijk+fdj],0))          fc+=1;      
                  if(!dequal(mask_B[fijk+fdk],0))          fc+=1;                    
               
                  if (!(dequal(fc,4))) {
                    printf("%e   \n", fc);
                    printf("%le %le %le %f  \n", xc[ijk]  , yc[ijk]  , zc[ijk],  mask_A[ijk]);      
                    printf("%le %le %le %le \n", xf[fijk] , yf[fijk] , zf[fijk], mask_B[fijk]);
                    printf("%le %le %le %le \n", xf[fijk+fdi+fdj+fdk] , yf[fijk+fdi+fdj+fdk] , 
                                                  zf[fijk+fdi+fdj+fdk], mask_B[fijk+fdi+fdj+fdk]);
                    printf("%le %le %le %le \n", xf[fijk+fdi] , yf[fijk+fdi] , zf[fijk+fdi], mask_B[fijk+fdi]);
                    printf("%le %le %le %le \n", xf[fijk+fdj] , yf[fijk+fdj] , zf[fijk+fdj], mask_B[fijk+fdj]);
                    printf("%le %le %le %le \n", xf[fijk+fdk] , yf[fijk+fdk] , zf[fijk+fdk], mask_B[fijk+fdk]);
                    printf("%le %le %le %le \n", xf[fijk+fdi+fdj] , yf[fijk+fdi+fdj] , zf[fijk+fdi+fdj], 
                                                  mask_B[fijk+fdi+fdj]);
                    printf("%le %le %le %le \n", xf[fijk+fdi+fdk] , yf[fijk+fdi+fdk] , zf[fijk+fdi+fdk], 
                                                  mask_B[fijk+fdi+fdk]);
                    printf("%le %le %le %le \n \n", xf[fijk+fdj+fdk] , yf[fijk+fdj+fdk] , zf[fijk+fdj+fdk], 
                                                   mask_B[fijk+fdj+fdk]);	
                  }

                 if ((abs(xf[fijk]-xc[ijk])+abs(yf[fijk]-yc[ijk])+abs(zf[fijk]-zc[ijk])) > (3.01*dl)) {
                    printf("UPS! \n");            	
                    printf("%le %le %le %f  \n" , xc[ijk]  , yc[ijk]  , zc[ijk],  mask_A[ijk]);      
                    printf("%le %le %le %le \n" , xc[ijk]  , yc[ijk]  , zc[ijk],  c[ijk]);
                    printf("%le %le %le %le \n" , xf[fijk] , yf[fijk] , zf[fijk], f[fijk]);
                    printf("%e %d %d \n", (abs(xf[fijk]-xc[fijk])+abs(yf[fijk]-yc[ijk])+abs(zf[fijk]-zc[ijk])), fdi, fijk); 
                    errorexit("something in the conservative amr is wrong, try to change amr_bo_dxmax SORRY!");
                 }
     
                 if(PR) printf("%e \n", c[ijk]);
               } /* close if(DEBUG)*/
            }  /* close if (mask_A[ijk]) */
         } endfor_ijk;
      } endforboxes;


      /* reset B-cells*/
      forallpoints_ijk(lf) {
         if (mask_B[ijk])   
            f[ijk]      = 0; 
      } endfor_ijk;


      /*add correction and reset A-cells*/
      forallpoints_ijk(lc) {
         if (mask_A[ijk]) { 
            cgrhd[ijk] += c[ijk]; 
            c[ijk]      = 0; 
         }
      } endfor_ijk;
 
   }/*close loop over variables */

   bampi_syncparent_send0123(lf, uc); /* FIXME: check what is needed, this here is to much and slows things down! */
   set_boundary_symmetry(lc, uc);

   if(PR) printf("I leave correct_varlist.\n"); 

   /* set matter_mask only if it is switched on */
   if(MATTER.USEMASK)
   {
     /* make sure everything is set up correct */
     /* WT: Tim, why do we need this? */
	MATTER.set_mask(lf);
	MATTER.set_mask(lc);
   }
}


int c_amr_time(tL *level) {

   double time      = level->time;
   double camr_time = Getd("camr_time");
   if (!(MATTER.CAMRACTIVE)) printf("camr not switched on \n");
   else  printf("camr ON \n");

   printf("time %e camrtime %e \n",time,camr_time);

   if(dequal(time,camr_time)){
      printf("activate conservative mesh refinement\n");
      c_amr_mask(level); 
      MATTER.CAMRACTIVE = 1;
   }   
}


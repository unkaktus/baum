/* extrapolate.c */
/* Wolfgang Tichy 1/2003 */


#include "bam.h"
#include "boundary.h"



#define centered_lin_extrapol(p,m) (0.5*(varnew[p] + varnew[m]))
#define right_lin_extrapol(p,P) (2.0*varnew[p] - varnew[P])
#define left_lin_extrapol(m,M)  (2.0*varnew[m] - varnew[M])

#define extrapolate(valx,use_x,Mcc,mcc,ccc,pcc,Pcc)    \
if(pcc==ccc)      { valx = left_lin_extrapol(mcc, Mcc);  use_x=1;    } \
else if(mcc==ccc) { valx = right_lin_extrapol(pcc, Pcc); use_x=1;    } \
else              { valx = centered_lin_extrapol(pcc, mcc); use_x=0; }



/* set extrapolate boundary for one variable */
void set_boundary_extrapolate(tL *level, int unew)
{
  double *varnew = level->v[unew];
  double valx, valy, valz, extrval;
  int use_x, use_y, use_z;
  int rep;
  
  /* we go thrice over the boundary, so that edges and corners,
     which are are computed from values on the boundary get computed form
     updated values */
  for(rep=1;rep<=3;rep++) 
    forboundary13(level, PHYBOUND) 
    {
      extrapolate(valx, use_x, Mcc, mcc, ccc, pcc, Pcc);
      extrapolate(valy, use_y, cMc, cmc, ccc, cpc, cPc);
      extrapolate(valz, use_z, ccM, ccm, ccc, ccp, ccP);

      if(use_x+use_y+use_z>0)
        extrval = (use_x*valx + use_y*valy + use_z*valz)/(use_x+use_y+use_z);
      else
        extrval = (valx + valy + valz)/3.0;

      varnew[ccc] = extrval;

    } endfor;
}


/* set extrapolate boundary on flagged boundary for a variable list */
void extrapolate_to_flagged_boundary(tL *level, tVarList *vlu, int Bflag)
{
  double *varnew;
  double valx, valy, valz, extrval;
  int use_x, use_y, use_z;
  int rep, j;
  
  /* we go thrice over the boundary, so that edges and corners,
     which are are computed from values on the boundary get computed form
     updated values */
  for(rep=1;rep<=3;rep++) 
    forboundary13(level, Bflag)
      for(j = 0; j < vlu->n; j++)
      {
        varnew = VLPtr(vlu, j);
    
        extrapolate(valx, use_x, Mcc, mcc, ccc, pcc, Pcc);
        extrapolate(valy, use_y, cMc, cmc, ccc, cpc, cPc);
        extrapolate(valz, use_z, ccM, ccm, ccc, ccp, ccP);

        if(use_x+use_y+use_z>0)
          extrval = (use_x*valx + use_y*valy + use_z*valz)/(use_x+use_y+use_z);
        else
          extrval = (valx + valy + valz)/3.0;

        varnew[ccc] = extrval;

      } endfor;
}


/* set extrapolate boundary for one variable for shells */
void set_boundary_extrapolate_shells(tL *level, tVarList *vlu)
{
  
  int order = Geti("boundary_order_extrapolate");
  int N     = Geti("boundary_N_extrapolate");
  double rb = level->bbox[1] - N*level->dx;
  double*rp = Ptr(level,"shells_R");
  double*rr = Ptr(level,"shells_r");
  if (order<2) return;
  
  bampi_openmp_start
  
  int i,j,k, n,o;
  double v[20],value,r;
  double *var;
  
  bampi_openmp_loop
  for (n=0; n<vlu->n; n++) {
    
    var = VLPtr(vlu, n);
    
    /* here we have luck!!! the radial index goes in that way that 
       we have increasing radius -> we can extrapolate point by point
       in outwards direction */
    /* DO NOT PARALLELIZE WITH OPENMP */
    forallpoints_ijk(level) {
      
      if (dless(rb,rr[ijk])) {
        
        if (i-order<0 || ijk-order*di<0)
          errorexit("no, extrapolate is not working, boxes are too small");
        
        for (o=0; o<order; o++) {
          v[o] = var[ijk-(o+1)*di];
        }
        
        var[ijk] = extrapolate_N(order, v);
      }
      
    } endfor_ijk;
    
  }
  
  bampi_openmp_stop
}







double extrapolate_N(int N, double *v)
{
    double value;
    
    if      (N==2)  value = 2.*v[0] -    v[1];
    else if (N==4)  value = 4.*v[0] - 6.*v[1] +  4.*v[2] -     v[3];
    else if (N==6)  value = 6.*v[0] -15.*v[1] + 20.*v[2] - 15.*v[3] +  6.*v[4] -     v[5];
    else if (N==8)  value = 8.*v[0] -28.*v[1] + 56.*v[2] - 70.*v[3] + 56.*v[4] - 28.*v[5] +  8.*v[6] -    v[7];
    else if (N==10) value =10.*v[0] -45.*v[1] +120.*v[2] -210.*v[3] +252.*v[4] -210.*v[5] +120.*v[6] -45.*v[7] +10.*v[8] -    v[9];
    else errorexit("order is not implemented");
    
    return value;
}

void extrapolate_to_outer_N_pts(tL *level, int ind)
{
    int pr = 0;
    double *var = level->v[ind];
    double v[12];
    int rep, off, o;
    int order = Geti("order_boundary");
    int N = Geti("order_boundary")/2-1;
    
    
    // test if there are enough points
    forallboxes(level) {
        if ((box->m<=2*off+order) || (box->n<=2*off+order) || (box->o<=2*off+order)) 
            errorexit("box to small for boundary condition -> reduce order");
    } endforboxes;
    
    
    // extrapolate layer for layer, beginning with the innermost layer
    for (off=N; off>=0; off--) {
        forallpoints_ijk(level) {
            /* go through all points which are in one layer 
               with off-1 pts away from the boundary
            */
            ifnewPHYBOUND(i,j,k, box, off);
            
            // compute extrapolation
            for (o=0; o<order; o++) 
                v[o] = var[ijk + (o+1)*(di*Di + dj*Dj + dk*Dk)];
                        
            // set extrapolated value
            var[ijk] = extrapolate_N(order, v);
            
        } endfor_ijk;
    }
    
}

void extrapolate_layer(tVarList *ucur)
{
    if (ucur->level->l!=0) return;
    int i;
    
    //first synchronize    
    bampi_vlsynchronize(ucur);
    
    // extrapolate all evolved variables which are outside the NEW boundary
    for (i = 0; i < ucur->n; i++) {
        extrapolate_to_outer_N_pts(ucur->level, ucur->index[i]);
    }
    
    //again synchronize
    bampi_vlsynchronize(ucur);
    
    set_boundary_symmetry(ucur->level, ucur); 
}














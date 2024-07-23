/* radiative.c */
/* Bernd Bruegmann 10/02, Wolfgang Tichy 12/2003 */

/* Radiative boundary condition
   This routine is based on Miguel Alcubierre's "new" radiative boundary
   condition in Cactus.
*/

#include "bam.h"
#include "boundary.h"



#define centeredderiv(p,m) (var[p] - var[m])
#define rightderiv(c,p,P) (4*var[p] - var[P] - 3*var[c])
#define leftderiv(c,m,M)  (3*var[c] + var[M] - 4*var[m])



#define computedel(xp,delx,vx,Mcc,mcc,ccc,pcc,Pcc)    \
if (vx > 0) {\
  if (Mcc != ccc)\
    delx = leftderiv(ccc, mcc, Mcc);\
  else if (mcc != ccc)\
    delx = centeredderiv(pcc, mcc);\
  else\
    delx = rightderiv(ccc, pcc, Pcc);\
} else if (vx < 0) {\
  if (Pcc != ccc) \
    delx = rightderiv(ccc, pcc, Pcc);\
  else if (mcc != ccc)\
    delx = centeredderiv(pcc, mcc);\
  else\
    delx = leftderiv(ccc, mcc, Mcc);\
} else {\
  if ( (pcc != ccc) && (mcc != ccc) ) \
    delx = centeredderiv(pcc, mcc);\
  else if (mcc == ccc)\
    delx = rightderiv(ccc, pcc, Pcc);\
  else\
    delx = leftderiv(ccc, mcc, Mcc);\
}

static double coeff[13][13] = {
  {0.},
  {0.},
  {-1./2.,0.,1./2.},
  {0.},
  {1./12., -2./3., 0., 2./3., -1./12.},
  {0.},
  {6., -1./60., 3./20., -3./4., 0., 3./4., -3./20., 1./60.,},
  {0.},
  {8., 1./280., -4./105., 1./5., -4./5., 0., 4./5., -1./5., 4./105., -1./280.}
};



/* set radiative boundary for one variable */
void set_boundary_radiative(tL *level, 
	                    int unew, int upre, double c, int ucur,
			    double var0, double v) 
{
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  double *var = level->v[ucur];
  double *varnew = level->v[unew];
  double *varpre = level->v[upre];
  double oo2dx = 1/(2*level->dx);
  double oo2dy = 1/(2*level->dy);
  double oo2dz = 1/(2*level->dz);
  double r, vx, vy, vz, x, y, z;
  double delx, dely, delz, rhs;

  forboundary13(level, PHYBOUND) {

    x = xp[ccc];
    y = yp[ccc];
    z = zp[ccc];
    r = sqrt(x*x + y*y + z*z);
    
    vx = v * x/r;
    vy = v * y/r;
    vz = v * z/r;

    rhs = -v*(var[ccc] - var0)/r; 

    computedel(xp, delx, vx, Mcc, mcc, ccc, pcc, Pcc);
    computedel(yp, dely, vy, cMc, cmc, ccc, cpc, cPc);
    computedel(zp, delz, vz, ccM, ccm, ccc, ccp, ccP);

    rhs -= oo2dx*vx*delx + oo2dy*vy*dely + oo2dz*vz*delz;

    if (c != 0.0) 
      varnew[ccc] = varpre[ccc] + c*rhs;
    else
      varnew[ccc] = rhs;

  } endfor;
}


/* WT: set radiative boundary for one variable if shift corresponds 
   to a rigid rotation  */
void set_boundary_rot_rad(tL *level, 
	                  int ivarnew, int ivarpre, double dt, int ivarcur,
	                  double var0, double v,
		  	  int ibetax_cur, int ibetax_far, int rot, 
		  	  int betaAdv) 
{
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];

  double *var;
  double *varcur = level->v[ivarcur];
  double *varnew = level->v[ivarnew];
  double *varpre = level->v[ivarpre];

  double *betax = level->v[ibetax_cur+0];
  double *betay = level->v[ibetax_cur+1];
  double *betaz = level->v[ibetax_cur+2];

  double *betax_far = level->v[ibetax_far+0];
  double *betay_far = level->v[ibetax_far+1];
  double *betaz_far = level->v[ibetax_far+2];
  
  double oo2dx = 1/(2*level->dx);
  double oo2dy = 1/(2*level->dy);
  double oo2dz = 1/(2*level->dz);
  double r, nx, ny, nz, x, y, z;
  double delx, dely, delz, rhs;
  double Wx, Wy, Wz;  /* W^i = v n^i  or  W^i = v n^i - beta^i */
  char *tensorindices;
  int var_is_beta=0;
  int i;

  if(0)
  {
     printf("ivarcur=%d    %s   var0=%g v=%g\n", 
            ivarcur, VarName(ivarcur), var0, v);
     printf("ibetax_cur=%d %s   ibetax_far=%d %s   rot=%d\n",
            ibetax_cur, VarName(ibetax_cur), 
            ibetax_far, VarName(ibetax_far), rot );
  }

  /* are we rotating? */
  if(rot)
  {
    /* get tensor info about variable */
    tensorindices = VarTensorIndices(ivarcur);
    
    /* If we are applying BCs to the shift, subtract beta_far. */
    if(ivarcur==ibetax_cur)
      forallpoints(level, i)
      {
        betax[i] -= betax_far[i];
        var_is_beta = 1;
      }
    if(ivarcur==ibetax_cur+1)
      forallpoints(level, i)
      {
        betay[i] -= betay_far[i];
        var_is_beta = 1;
      }
    if(ivarcur==ibetax_cur+2)
      forallpoints(level, i)
      {
        betaz[i] -= betaz_far[i];
        var_is_beta = 1;
      }
  }
  
  /*
  forallpoints_ijk(level) {
      printf("%d   %d  %e %e\n",ijk,level->boundary[ijk], zp[ijk],level->bbox[5]);
    } endfor_ijk;
  exit(0);
  */
  /* loop over physical boundary */
  forboundary13(level, PHYBOUND)
  {
    x = xp[ccc];
    y = yp[ccc];
    z = zp[ccc];
    r = sqrt(x*x + y*y + z*z);
    
    //printf(" %d %d %d %d %d\n", Mcc, mcc, ccc, pcc, Pcc);
    
    
    Wx = v*x/r;
    Wy = v*y/r;
    Wz = v*z/r;
    
    var=varcur;

    if(!rot || var_is_beta)
    {
      /* standard radiative BC, with advection ideas */
      computedel(xp, delx, Wx, Mcc, mcc, ccc, pcc, Pcc);
      computedel(yp, dely, Wy, cMc, cmc, ccc, cpc, cPc);
      computedel(zp, delz, Wz, ccM, ccm, ccc, ccp, ccP);
      rhs  = -v * ( var[ccc] - var0 )/r;
      rhs += -( Wx*delx*oo2dx + Wy*dely*oo2dy + Wz*delz*oo2dz );
    }
    else  /* are we rotating? , i.e.: if(rot && !var_is_beta) */
    {
      /* which version of stencil do we use??? */
      if(betaAdv)
      {
        Wx -= betax[ccc];
        Wy -= betay[ccc];
        Wz -= betaz[ccc];

        /* radiative BC, with advection along W = v*n - beta */
        computedel(xp, delx, Wx, Mcc, mcc, ccc, pcc, Pcc);
        computedel(yp, dely, Wy, cMc, cmc, ccc, cpc, cPc);
        computedel(zp, delz, Wz, ccM, ccm, ccc, ccp, ccP);
        rhs  = -v * ( var[ccc] - var0 )/r;
        rhs += -( Wx*delx*oo2dx + Wy*dely*oo2dy + Wz*delz*oo2dz );
      }
      else
      {
        /* radiative BC, with advection along W = v*n */
        computedel(xp, delx, Wx, Mcc, mcc, ccc, pcc, Pcc);
        computedel(yp, dely, Wy, cMc, cmc, ccc, cpc, cPc);
        computedel(zp, delz, Wz, ccM, ccm, ccc, ccp, ccP);
        rhs  = -v * ( var[ccc] - var0 )/r;
        rhs += -( Wx*delx*oo2dx + Wy*dely*oo2dy + Wz*delz*oo2dz );

        /* correction to any scalar or tensor if rotational shift,
           using centerd derivs if possible                        */
        computedel(xp, delx, 0, Mcc, mcc, ccc, pcc, Pcc);
        computedel(yp, dely, 0, cMc, cmc, ccc, cpc, cPc);
        computedel(zp, delz, 0, ccM, ccm, ccc, ccp, ccP);
        rhs += betax[ccc]*delx*oo2dx + betay[ccc]*dely*oo2dy + 
               betaz[ccc]*delz*oo2dz;
      }


      /* check if we have to do tensor corrections */
      if(strcmp(tensorindices, "")!=0)
      {
        double dbeta[4][4]; /* dbeta[i][j] = $\beta_{i,j}$ */
       
        /* compute partial derivs of beta  */
        var=betax;
        computedel(xp, delx, 0, Mcc, mcc, ccc, pcc, Pcc);
        computedel(yp, dely, 0, cMc, cmc, ccc, cpc, cPc);
        computedel(zp, delz, 0, ccM, ccm, ccc, ccp, ccP);
        dbeta[1][1] = delx*oo2dx;
        dbeta[1][2] = dely*oo2dy;
        dbeta[1][3] = delz*oo2dz;

        var=betay;
        computedel(xp, delx, 0, Mcc, mcc, ccc, pcc, Pcc);
        computedel(yp, dely, 0, cMc, cmc, ccc, cpc, cPc);
        computedel(zp, delz, 0, ccM, ccm, ccc, ccp, ccP);
        dbeta[2][1] = delx*oo2dx;
        dbeta[2][2] = dely*oo2dy;
        dbeta[2][3] = delz*oo2dz;

        var=betaz;
        computedel(xp, delx, 0, Mcc, mcc, ccc, pcc, Pcc);
        computedel(yp, dely, 0, cMc, cmc, ccc, cpc, cPc);
        computedel(zp, delz, 0, ccM, ccm, ccc, ccp, ccP);
        dbeta[3][1] = delx*oo2dx;
        dbeta[3][2] = dely*oo2dy;
        dbeta[3][3] = delz*oo2dz;
        
        var=varcur;
      
        /* corrections if variable varcur is a tensor 
           $  B^i, B_i, B^{ij}, B_{ij}, ...  $       */
        if(strcmp(tensorindices, "I")==0)
        {
          int k,i;
          int Comp  = VarComponent(ivarcur);
          double B[4];
          
          i=Comp+1;
          
          B[1] = (level->v[ivarcur-Comp])[ccc];
          B[2] = (level->v[ivarcur-Comp+1])[ccc];
          B[3] = (level->v[ivarcur-Comp+2])[ccc];
    
          for(k=1; k<=3; k++)
            rhs -= dbeta[i][k] * B[k];
        }
        else if(strcmp(tensorindices, "i")==0)
        {
          int k,i;
          int Comp  = VarComponent(ivarcur);
          double B[4];
          
          i=Comp+1;
          
          B[1] = (level->v[ivarcur-Comp])[ccc];
          B[2] = (level->v[ivarcur-Comp+1])[ccc];
          B[3] = (level->v[ivarcur-Comp+2])[ccc];
    
          for(k=1; k<=3; k++)
            rhs += dbeta[k][i] * B[k];
        }
        else if(strcmp(tensorindices, "IJ+JI")==0)
        {
          int k,i,j;
          int Comp  = VarComponent(ivarcur);
          double B[4][4];
          
          if(Comp==5)                  { i=3; j=3;}
          else if(Comp==3 || Comp==4)  { i=2; j=Comp-1;}		
          else                         { i=1; j=Comp+1;}		
          
          B[1][1] 	  = (level->v[ivarcur-Comp])[ccc];
          B[1][2] = B[2][1] = (level->v[ivarcur-Comp+1])[ccc];
          B[1][3] = B[3][1] = (level->v[ivarcur-Comp+2])[ccc];
          B[2][2] 	  = (level->v[ivarcur-Comp+3])[ccc];
          B[2][3] = B[3][2] = (level->v[ivarcur-Comp+4])[ccc];
          B[3][3] 	  = (level->v[ivarcur-Comp+5])[ccc];
                  	  
          for(k=1; k<=3; k++)
            rhs -= dbeta[i][k] * B[k][j] + dbeta[j][k] * B[i][k];
        }
        else if(strcmp(tensorindices, "ij+ji")==0)
        {
          int k,i,j;
          int Comp  = VarComponent(ivarcur);
          double B[4][4];
          
          if(Comp==5)                  { i=3; j=3;}
          else if(Comp==3 || Comp==4)  { i=2; j=Comp-1;}		
          else                         { i=1; j=Comp+1;}		
          
          B[1][1] 	  = (level->v[ivarcur-Comp])[ccc];
          B[1][2] = B[2][1] = (level->v[ivarcur-Comp+1])[ccc];
          B[1][3] = B[3][1] = (level->v[ivarcur-Comp+2])[ccc];
          B[2][2] 	  = (level->v[ivarcur-Comp+3])[ccc];
          B[2][3] = B[3][2] = (level->v[ivarcur-Comp+4])[ccc];
          B[3][3] 	  = (level->v[ivarcur-Comp+5])[ccc];

          for(k=1; k<=3; k++)
            rhs += dbeta[k][i] * B[k][j] + dbeta[k][j] * B[i][k];
        }
        else
        {
          printf("set_boundary_rot_rad: don't know how to treat %s "
                 "with tensorindices %s\n", VarName(ivarcur), tensorindices);
          errorexit("set_boundary_rot_rad: can't cope with inidces");
        }
      } /* end: if(strcmp(tensorindices, "")!=0) */
      
      
    }
    
    /* add RHS to var on new time level */    
    if (dt != 0.0) 
      varnew[ccc] = varpre[ccc] + dt*rhs;
    else
      varnew[ccc] = rhs;

  } endfor;

  if(rot)
  {
    /* if we have applied BCs to the shift, add beta_far back in */
    if(ivarcur==ibetax_cur)
      forallpoints(level, i) { betax[i] += betax_far[i]; }

    if(ivarcur==ibetax_cur+1)
      forallpoints(level, i) { betay[i] += betay_far[i]; }

    if(ivarcur==ibetax_cur+2)
      forallpoints(level, i) { betaz[i] += betaz_far[i]; }
  }
}


/* set radiative boundary for one variable speciallized for spherical shells*/
void set_boundary_radiative_shells(tL *level, 
                            int unew, int upre, double c, int ucur,
                            double var0, double v) 
{
  double *rp  = level->v[Ind("shells_R")];
  double *rr  = level->v[Ind("shells_r")];
  double *var = level->v[ucur];
  double *varnew = level->v[unew];
  double *varpre = level->v[upre];
  double oo2dr = 1/(2*level->dx);
  double oo2dp = 1/(2*level->dy);
  double oo2dt = 1/(2*level->dz);
  double delr, rhs;
  
  forboundary13(level, PHYBOUND) {

    rhs = -v*(var[ccc] - var0)/rp[ccc]; 

    computedel(rp[ccc], delr, v, Mcc, mcc, ccc, pcc, Pcc);
    
    // for fisheye 
    delr /= compute_fisheye(rr[ccc],1);
   
    rhs -= oo2dr*v*delr;

    if (c != 0.0) 
      varnew[ccc] = varpre[ccc] + c*rhs;
    else
      varnew[ccc] = rhs;

  } endfor;
}


/* set radiative boundary for one variable speciallized for spherical shells*/
void set_boundary_radiative_shells_centered(tL *level, 
                            int unew, int upre, double c, int ucur,
                            double var0, double v) 
{
  double *var = level->v[ucur];
  double *varnew = level->v[unew];
  double *varpre = level->v[upre];
  double oodr = 1./(level->dx);
 

  //int order = Geti("boundary_order_finitediff");
  int order = Geti("order_centered");
  int N     = Geti("boundary_N_extrapolate");
  double rb = level->bbox[1] - N*level->dx;
  double*rp = Ptr(level,"shells_R");
  double*rr = Ptr(level,"shells_r");
  if (order<2) return;
  

  bampi_openmp_start
  
  double delr, rhs;
  int o;
  
  forallpoints_ijk_openmp(level) {
    
    if (dequal(rb,rr[ijk])) {

      if (i-order/2<0 || i+order/2>imax)
        errorexit("no, extrapolate is not working, boxes are too small");
      
      delr = 0.;
      for (o=0; o<=order; o++) {
        delr += coeff[order][o] * var[ijk + (-order/2+o)*di];
      }
      
      // for fisheye
      delr /= compute_fisheye(rr[ijk],1);
      
      rhs = -v*(var[ijk] - var0)/rp[ijk] - oodr*v*delr;
      
      if (c != 0.0) 
        varnew[ijk] = varpre[ijk] + c*rhs;
      else
        varnew[ijk] = rhs;
    }
    
  } endfor_ijk_openmp;
  bampi_openmp_stop
  
}















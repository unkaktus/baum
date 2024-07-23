/* OutsideBox.c,  3/2003 Wolfgang Tichy
 
   integrate a function func(double x, double y, double z) over a box
   larger than the grid, with the grid excised from the region of integration
*/


#include <stdio.h>
#include <stdlib.h>


/* meine Funktionen */
static double (*func)(double,double,double);
static double func_yxz(double y, double x, double z);
static double func_zxy(double z, double y, double x);

/* domain boundaries */
static double y1_x(double x);
static double y2_x(double x);
static double z1_xy(double x, double y);
static double z2_xy(double x, double y);
static double x1_y(double y);
static double x2_y(double y);
static double z1_yx(double y, double x);
static double z2_yx(double y, double x);
static double x1_z(double z);
static double x2_z(double z);
static double y1_zx(double z, double x);
static double y2_zx(double z, double x);

/* Integrators */ 
extern double integral(double (*f)(double), double a, double b, double s, int max);
extern double rombintegral(double (*f)(double), double a, double b, double s, int max);
extern double integral3D(double (*int_meth)(double (*f_int)(double),
                                            double a, double b,
                                            double eps, int max),
                         double (*func)(double, double, double),
                         double x_limit1, double x_limit2,
                         double (*y_limit1)(double x), double (*y_limit2)(double x),
                         double (*z_limit1)(double x,double y),
                         double (*z_limit2)(double x,double y),
                         double sx, double sy, double sz,
                         int maxx, int maxy, int maxz);
                         
                                                                                                                                                                                                        
/* Variables */
static double xmax, ymax, zmax;


double OutsideBox(double grid_xmax, double grid_ymax, double grid_zmax, 
                         double xout, double (*fn)(double,double,double) )
{
 /* meine Variablen */
 double xp_dom, xm_dom, yp_dom, ym_dom, zp_dom, zm_dom, dom_sum;
 double sx,sy,sz;
 int maxx,maxy,maxz;
 
 double magfac, yout, zout;

 /* set func to the fn which is passed in */
 func=fn;
 xmax=grid_xmax;
 ymax=grid_ymax;
 zmax=grid_zmax;
  
 /* outer boundary of large box */
 if(xmax!=0)
 {
   magfac=xout/xmax;
   yout=ymax*magfac;
   zout=zmax*magfac; 
 }
 else
 {
   yout=xout;
   zout=xout;
 }
        
 
 /* desired rel. errors */ 
 sx=sy=sz=1e-4;
 
 /* log_2 of max. number of points in x,y and z direction */
 maxx=maxy=maxz=20;
 
 
 printf("OutsideBox: \n");
 printf("Box outside grid:  xout=%f  yout=%f  zout=%f\n",xout,yout,zout);
 printf("Grid box:          xmax=%f  ymax=%f  zmax=%f\n",xmax,ymax,zmax);
 printf("fn=%p\n",fn);

 /* integrate over the six domains:
     xp_dom <-> positive x diriection
     zm_dom <-> negative z diriection
 */
 
 /* TRY TRAPEZ */
/*
 printf("Result with integral:\n");
 xp_dom=integral3D(integral, func, xmax,xout, y1_x, y2_x, z1_xy, z2_xy,
                   sx,sy,sz, maxx,maxy,maxz); 
 printf("xp_dom=%f\n",xp_dom);

 xm_dom=integral3D(integral, func, -xout,-xmax, y1_x, y2_x, z1_xy, z2_xy,
                   sx,sy,sz, maxx,maxy,maxz); 
 printf("xm_dom=%f\n",xm_dom);


 yp_dom=integral3D(integral, func_yxz, ymax,yout, x1_y, x2_y, z1_yx, z2_yx,
                   sx,sy,sz, maxx,maxy,maxz); 
 printf("yp_dom=%f\n",yp_dom);

 ym_dom=integral3D(integral, func_yxz, -yout,-ymax, x1_y, x2_y, z1_yx, z2_yx,
                   sx,sy,sz, maxx,maxy,maxz); 
 printf("ym_dom=%f\n",ym_dom);


 zp_dom=integral3D(integral, func_zxy, zmax,zout, x1_z, x2_z, y1_zx, y2_zx,
                   sx,sy,sz, maxx,maxy,maxz); 
 printf("zp_dom=%f\n",zp_dom);

 zm_dom=integral3D(integral, func_zxy, -zout,-zmax, x1_z, x2_z, y1_zx, y2_zx,
                   sx,sy,sz, maxx,maxy,maxz); 
 printf("zm_dom=%f\n",zm_dom);
 
 dom_sum=xp_dom + xm_dom + yp_dom + ym_dom + zp_dom + zm_dom;
 printf("xp_dom + xm_dom + yp_dom + ym_dom + zp_dom + zm_dom = %f\n", dom_sum);
*/                     

 /* TRY ROMBERG: */
 printf("Result with rombintegral:\n");
 xp_dom=integral3D(rombintegral, func, xmax,xout, y1_x, y2_x, z1_xy, z2_xy,
                   sx,sy,sz, maxx,maxy,maxz); 
 xm_dom=integral3D(rombintegral, func, -xout,-xmax, y1_x, y2_x, z1_xy, z2_xy,
                   sx,sy,sz, maxx,maxy,maxz); 

 yp_dom=integral3D(rombintegral, func_yxz, ymax,yout, x1_y, x2_y, z1_yx, z2_yx,
                   sx,sy,sz, maxx,maxy,maxz); 
 ym_dom=integral3D(rombintegral, func_yxz, -yout,-ymax, x1_y, x2_y, z1_yx, z2_yx,
                   sx,sy,sz, maxx,maxy,maxz); 

 zp_dom=integral3D(rombintegral, func_zxy, zmax,zout, x1_z, x2_z, y1_zx, y2_zx,
                   sx,sy,sz, maxx,maxy,maxz); 
 zm_dom=integral3D(rombintegral, func_zxy, -zout,-zmax, x1_z, x2_z, y1_zx, y2_zx,
                   sx,sy,sz, maxx,maxy,maxz); 

 printf("xp_dom=%f  xm_dom=%f\n",xp_dom,xm_dom);
 printf("yp_dom=%f  ym_dom=%f\n",yp_dom,ym_dom);
 printf("zp_dom=%f  zm_dom=%f\n",zp_dom,zm_dom);
 
 dom_sum=xp_dom + xm_dom + yp_dom + ym_dom + zp_dom + zm_dom;
 printf("xp_dom + xm_dom + yp_dom + ym_dom + zp_dom + zm_dom = %f\n", dom_sum);

 return dom_sum;
} 





/* domain decomposition into Pyramidenstuempfe: */

/* func with variables swaped */
static double func_yxz(double y, double x, double z)
{
 return func(x,y,z);
}
static double func_zxy(double z, double x, double y)
{         
 return func(x,y,z);
}       
 
 
/* domain boundaries */ 
static double y1_x(double x)
{
 return -(ymax/xmax)*x;
}

static double y2_x(double x)
{
 return (ymax/xmax)*x;
}


static double z1_xy(double x, double y)
{
 return -(zmax/xmax)*x;
}

static double z2_xy(double x, double y)
{
 return (zmax/xmax)*x;
}

 
static double x1_y(double y)
{
 return -(xmax/ymax)*y;
}

static double x2_y(double y)
{
 return (xmax/ymax)*y;
}


static double z1_yx(double y, double x)
{
 return -(zmax/ymax)*y;
}

static double z2_yx(double y, double x)
{
 return (zmax/ymax)*y;
}


 
static double x1_z(double z)
{
 return -(xmax/zmax)*z;
}

static double x2_z(double z)
{
 return (xmax/zmax)*z;
}


static double y1_zx(double z, double x)
{
 return -(ymax/zmax)*z;
}

static double y2_zx(double z, double x)
{
 return (ymax/zmax)*z;
}


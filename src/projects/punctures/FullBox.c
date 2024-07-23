/* FullBox.c,  3/2003 Wolfgang Tichy
 
   integrate a function func(double x, double y, double z) over a box
   larger than the grid
*/


#include <stdio.h>
#include <stdlib.h>


/* meine Funktionen */
static double (*func)(double,double,double);

/* domain boundaries */
static double y1_x(double x);
static double y2_x(double x);
static double z1_xy(double x, double y);
static double z2_xy(double x, double y);

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
static double ymax, zmax;


double FullBox(double xout, double yout, double zout, 
               double (*fn)(double,double,double) )
{
 /* meine Variablen */
 double I;
 double sx,sy,sz;
 int maxx,maxy,maxz;
 

 /* set func to the fn which is passed in */
 func=fn;
 ymax=yout;
 zmax=zout;
  
 /* desired rel. errors */ 
 sx=sy=sz=1e-4;
 
 /* log_2 of max. number of points in x,y and z direction */
 maxx=maxy=maxz=20;


 printf("FullBox: \n");
 printf("Full box grid:  xout=%f  yout=%f  zout=%f\n",xout,yout,zout);
 printf("fn=%p\n",fn);

 printf("Result with rombintegral:\n");
 I=integral3D(rombintegral, func, -xout,xout, y1_x, y2_x, z1_xy, z2_xy,
                   sx,sy,sz, maxx,maxy,maxz); 
 printf("I = %f\n", I);

 return I;
} 





/* domain boundaries */ 
static double y1_x(double x)
{
 return -ymax;
}

static double y2_x(double x)
{
 return ymax;
}


static double z1_xy(double x, double y)
{
 return -zmax;
}

static double z2_xy(double x, double y)
{
 return zmax;
}


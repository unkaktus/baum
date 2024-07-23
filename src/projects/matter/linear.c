/* linear.c 
   linear interpolation routine for matter RP in AMR */
/* sbernuz 03/13 */


#include "bam.h"
#include "interpolate.h"

#define sqr(x) ((x)*(x))

/* important note

   prototype for the "1d interpolate call" (see interpolate_boxtobox_N.c)

   double  interpolate_xxx(int N, double x, double x0, double h, double *c, double *u)

   assume 
   - equidistant points x_i = x0 + i*h, i = 0 ... N-1
   - N = interpolation order, 
   gives
   P(x) = Sum c_i(x) * u(x_i)
   
   differently from other routines here below the interpolation order is always N=1 
   the variable N however usually contains the value set by "order_RP", 
   and it used e.g. for the largangin interpolation of the metric part.
   to bracket the interpolation point with the two first NN we assume N is even

   e.g. given x = interpolation point, we have always

   N = 2
   0      1   
   o  x   o                

   N = 4
   0      1       2      3    
   o      o    x  o      o       

   N = 6
   0      1       2      3       4      5 
   o      o       o x    o       o      o       

*/

double interpolate_linear(int N, double x, double x0, double h, double *c, double *u)
{
 
  // linear interpolation 
  // => average for conservative restriction 
    
  int i = N/2-1;
  double ooh = 1./h;
  double dx = x - x0 - i*h; // x - x_i
  double ux = u[i] + (u[i+1]  - u[i])*ooh*dx;
  
  return ux;

}


// VanLeer MC limiter 
double MClim( double x, double y ) 
{ 

  double min = DMIN( (2.0*fabs(x)), (2.0*fabs(y)) );
  return( 0.5*( SIGN( 1.0, x ) + SIGN( 1.0, y ) )*DMIN( min, (0.5*fabs(x+y)) ) );
  
}

double interpolate_linear_lim(int N, double x, double x0, double h, double *c, double *u)
{
  
  // linear interpolation with MC limiter
  // 3 pts stencil => shift around nearest pt
  // => prolongation

  int i = N/2-1;
  double ooh = 1./h;

  double dx  = x - x0 - i*h;     // x - x_i
  double dx1 = x - x0 - (i+1)*h; // x - x_i+1

  if (fabs(dx)>fabs(dx1)) {
    i++;
    dx = dx1;
  }  
  
  double slope = MClim( u[i]-u[i-1], u[i+1]-u[i] );
  double ux = u[i] + slope*ooh*dx;
  
  return ux;

}


























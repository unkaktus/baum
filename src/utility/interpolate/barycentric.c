/* barycentric.c */
/* Bernd Bruegmann 2/2007 */

/* Barycentric interpolation based on Berrett (?) and Trefethen 2004
   Has algorithmic advantages over the standard Lagrange formula.

   Intended use: 
     - random access, where there is no additional info about the distribution 
       of interpolation points which would allow additional optimization
     - variable order, although the loops could be written out for specific
       orders for optimization as in lagrange.c
*/

#if 1
#include "bam.h"
#include "interpolate.h"
#else
#include <stdio.h>
#include <math.h>
#endif




/* barycentric interpolation
   x0:   interpolate to this coordinate
   n:    number of points
   x[]:  coordinate of points
   y[]:  values at points
   s:    stride in point list (useful for multidimensional boxes)

   stable even close to x[], say x0 = x[i] + 1e-12
   could cache the omegas
*/
double barycentric_single(double x0, int n, int s, double *x, double *y)
{
  double omega[10];
  double a, b, c, d, o;
  int i, j;

  for (i = 0; i < n; i += s) {

    if (x0 - x[i] == 0) return y[i];

    o = 1;
    for (j = 0; j < n; j += s) {
      if (j != i) {
	o /= (x[i] - x[j]);
      }
    }
    omega[i/s] = o;
  }

  a = b = 0;
  for (i = 0; i < n; i += s) {
    d = x0 - x[i];
    c = omega[i/s]/d;
    b += c;
    a += c * y[i];
  }

  return a/b;
}




/* compute omega[] for barycentric interpolation */
void barycentric_omega(int n, int s, double *x, double *omega)
{
  double o;
  int i, j;

  if (0) printf("%d %d %p %p\n", n, s, x, omega);

  for (i = 0; i < n; i += s) {
    o = 1;
    for (j = 0; j < n; j += s) {
      if (j != i) {
	o /= (x[i] - x[j]);
      }
    }
    omega[i/s] = o;

    if (0) printf("x[%d] = %9.6f omega[%d] = %13.6e\n", i/s, x[i], i/s, o);
  }
}




/* barycentric interpolation with precomputed omega */
double barycentric(double x0, int n, int s, double *x, double *y,
		   double *omega)
{
  double a, b, c, d;
  int i, j;

  if (0) printf("%f %d %d %p %p %p\n", x0, n, s, x, y, omega);

  a = b = 0;
  for (i = 0; i < n; i += s) {
    d = x0 - x[i];
    if (d == 0) return y[i];
    c = omega[i/s]/d;
    b += c;
    a += c * y[i];
  }

  return a/b;
}




/* find index such that xp[i] <= x < xp[i+1]
   uses bisection, which relies on x being ordered
   o is "offset", number of points smaller than x that are required
   returns j = i-(o-1), i.e. if o = 2, then
     xp[j] < xp[j+1] <= x < xp[j+2] < xp[j+3]
   which is useful for interpolation
*/
int find_point_bisection(double x, int n, double *xp, int o)
{
  int i0 = o-1, i1 = n-o;
  int i;

  if (n < 2*o) errorexit("bisection failed");

  if (x <= xp[i0]) return 0;
  if (x >  xp[i1]) return n-2*o;

  while (i0 != i1-1) {
    i = (i0+i1)/2;
    if (x < xp[i]) i1 = i; else i0 = i;
  }

  return i0-o+1;
}




/* three dimensional polynomial interpolation, barycentric */
double interpolate_tri_bar(const int order, double x, double y, double z,
			   int n1, int n2, int n3, 
			   double *x1, double *x2, double *x3, double *yp)
{
  double u, v[order][order], w[order], omega[order];
  int i, j, k, ijk;
  int i1, i2, i3;
  int di = 1, dj = n1, dk = n1*n2;
  int order1 = order > n1 ? n1 : order;
  int order2 = order > n2 ? n2 : order;
  int order3 = order > n3 ? n3 : order;
  
  i1 = find_point_bisection(x, n1, x1, order1/2);
  i2 = find_point_bisection(y, n2, x2, order2/2);
  i3 = find_point_bisection(z, n3, x3, order3/2);
  ijk = i1*di + i2*dj + i3*dk;
  if (0) printf("%d %d %d\n", i1, i2, i3);

  barycentric_omega(order1, 1, &x1[i1], omega);
  for (k = 0; k < order3; k++)
  for (j = 0; j < order2; j++)
    v[k][j] = barycentric(x, order1, 1, &x1[i1], &yp[ijk+j*dj+k*dk], omega);

  if (0) 
  for (k = 0; k < order3; k++)
  for (j = 0; j < order2; j++)
    printf("%2d %2d   %.15f\n", k, j, v[k][j]);


  barycentric_omega(order2, 1, &x2[i2], omega);
  for (k = 0; k < order3; k++)
    w[k] = barycentric(y, order2, 1, &x2[i2], &v[k][0], omega);

  if (0) 
  for (k = 0; k < order3; k++)
    printf("%2d %.15f\n", k, w[k]);


  barycentric_omega(order3, 1, &x3[i3], omega);
  u = barycentric(z, order3, 1, &x3[i3], w, omega);

  return u;
}




#if 0
main()
{
  const int n1 = 20;
  const int n2 = 40;
  const int n3 = 40;
  double x1[n1], x2[n2], x3[n3];
  double y[n1*n2*n3];
  double a, b, c, h, u;
  int i, j, k, n;

  for (k = 0; k < n1; k++) x1[k] = sqrt(1 + 2.123*k/(n1-1.0));//1+k/(n1-1.0);//sqrt(1 + 2.123*k/(n1-1.0));
  for (k = 0; k < n2; k++) x2[k] = sqrt(1 + 2.123*k/(n2-1.0));//2 + k/(n2-1.0);
  for (k = 0; k < n3; k++) x3[k] = sqrt(1 + 2.123*k/(n3-1.0));//3 + k/(n3-1.0);

  if (0) {
    for (k = 0; k < n1; k++) x1[k] = -n1/2 + k;
    for (k = 0; k < n2; k++) x2[k] = -n2/2 + k;
    for (k = 0; k < n3; k++) x3[k] = -n3/2 + k;
  }    

  for (k = 0; k < n3; k++) printf("%d  %f\n", k, x3[k]);

  for (k = 0; k < n3; k++)
  for (j = 0; j < n2; j++)
  for (i = 0; i < n1; i++)
    y[i + j*n1 + k*n1*n2] = sqrt(x1[i]*x1[i] + x3[k]*x3[k]);

  a = 1.25;
  b = 1.25;
  c = 1.5;
  u = interpolate_tri_bar(10, a, b, c, n1, n2, n3, x1, x2, x3, y);

  printf("%.2f %.2f %.2f  u = %.15f\n", a, b, c, u);

  // constant 1.5 in x
  // order 2:  u = 1.499999999999986
  // order 4:  u = 1.500000000058193
  // order 6:  u = 1.500000047475998
  // order 10: u = 1.500000505467133

}
#endif

#if 0

main()
{
  double x[1000];
  double y[1000];
  double h, u;
  int i, j, n;

  if (0) {
    n = 6;
    for (i = 0; i < n; i++) {
      x[i] = i + 1;
      y[i] = sqrt(x[i]); // x[i]*x[i];
      printf("%6.2f %6.2f\n", x[i], y[i]);
    }
    printf("\n");
    
    for (h = 0.5; h < 3; h += 0.1) {
      
      u = barycentric_single(h, n, 1, x, y);
      
      printf("%6.2f %6.2f  %e\n", h, u, u - sqrt(h));
      
    }
  }

  if (1) {
    n = 100;
    for (i = 0; i < n; i++) 
      x[i] = i; //sqrt(i+1);

    h = 99;
    i = 5;
    j = find_point_bisection(h, n, x, i);

    printf("x[%d]=%.2f <= %.2f < x[%d]=%.2f\n",
	   j, x[j], h, j+2*i-1, x[j+2*i-1]);
  }
}
#endif

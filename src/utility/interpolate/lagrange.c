/* lagrange.c */
/* Bernd Bruegmann 3/03 */

#include "bam.h"
#include "interpolate.h"



/* set coefficients for (N-1)-th order, N point Lagrange polynomial
   assumes equidistant points xi = xmin+i*h, P(x) = Sum ci(x) * u(xi)
   for faster evaluation and for x near xi, see barycentric.c
*/
void coefficients_lagrange_N(int n, double x, double xmin, double h, double *c)
{
  int i, j;
  double d;

  for (i = 0; i < n; i++) {
    d = 1;
    for (j = 0; j < n; j++) {
      if (j != i) {
        d *= (x - xmin - j*h) / ((i-j)*h);
      }
    }
    c[i] = d;
  }
} 


double interpolate_lagrange_N(int N, double x, double x0, double h, double *c, double *u)
{
  int i;
  double sum = 0;
    
  for (i=0; i<N; i++)
    sum += c[i]*u[i];
    
  return sum;
}







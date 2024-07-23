/* minimization.c */
/* Bernd Bruegmann 11/03 */

/* test minimization routines from Numerical Recipes 

  cc minimization.c brent.c golden.c nrutil.c linmin.c powell.c mnbrak.c -lm

*/

#include <math.h>
#include <stdio.h>
#include "nrutil.h"

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	     double *xmin);
double golden(double ax, double bx, double cx, double (*f)(double), double tol,
	      double *xmin);
void powell(double p[], double **xi, int n, double ftol, int *iter,
	    double *fret, double (*func)(double []));

int nn;



/* test functions */
double f1(double x)
{
  double f;

  f = pow(x-0.5, 4) + 0;
  f = -sin(x);
  f = -exp(-(x-0.09)*(x-0.09)/0.01);

  printf("%3d:  f(%12.9f) = %12.9f\n", nn++, x, f);

  return f;
}

/* test functions 2d */
double f2(double *x)
{
  double f;
  double x1 = x[1], x2 = x[2];

  f = pow(x1*x1 + x2*x2 + 1, 2);
  f = - exp(-x1*x1 - x2*x2);

  printf("%3d:  f(%19.12e,%19.12e) = %19.12e\n", nn++, x1, x2, f);

  return f;
}




/* main */
int main(int argc, char **argv)
{
  double min, fmin;
  double tol = 1e-6;
  const int n = 2;

  if (1) {
    nn = 0;
    fmin = brent(-2, 0, 3, f1, tol, &min);
    printf("min = %e,  fmin-1 = %e\n", min, fmin-1);

    nn = 0;
    fmin = golden(-2, 0, 3, f1, tol, &min);
    printf("min = %e,  fmin-1 = %e\n", min, fmin-1);

    nn = 0;
    fmin = brent(0.05, 0.10, 0.15, f1, tol, &min);
    printf("min = %e,  fmin-1 = %e\n", min, fmin-1);
  }

  if (0) {
    double *p = dvector(1, n);
    double **xi = dmatrix(1, n, 1, n);
    int i, j;

    for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++) 
      xi[i][j] = (i == j);

    p[1] = p[2] = 1;

    powell(p, xi, n, tol, &i, &fmin, f2);
  }

  return 0;
}

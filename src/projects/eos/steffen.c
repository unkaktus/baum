#include <stdlib.h>
#include <math.h>


//===================================== Steffen interpolation routines (monotonized cubic splines) ========================//
int steffen_init (double *xa, double *ya, int N, double **dya, double **a, double **b, double **c, double **d)
{
	
	int i;
	double k, pi, hi, him1, si, sim1;


	*dya = (double *) malloc (N*sizeof(double));
	*a  = (double *) malloc ((N-1)*sizeof(double));
	*b  = (double *) malloc ((N-1)*sizeof(double));
	*c  = (double *) malloc ((N-1)*sizeof(double));
	*d  = (double *) malloc ((N-1)*sizeof(double));

	double h0 = xa[1] - xa[0];
	double s0 = (ya[1] - ya[0])/h0;
	
	// Left boundary derivative
	(*dya)[0] = s0;
	
	// Compute all necessary coeff. and derivs except the boundaries

	for(i=1;i<N-1;i++){
		hi   = (xa[i+1] - xa[i]  );
		him1 = (xa[i]   - xa[i-1]);

		si   = (ya[i+1] - ya[i]  )/hi;
		sim1 = (ya[i]   - ya[i-1])/him1;

		pi = (sim1*hi + si*him1)/(him1 + hi);

		// Let the checks begin ( Eq. 11 of Steffen )

		if( (sim1*si) <= 0.){
			(*dya)[i] = 0.;
			continue;
		}
		else if( ( fabs(pi) > 2.*fabs(sim1) ) || ( fabs(pi) > 2.*fabs(si) ) ){
			k = (sim1 < 0.) ? -1. : 1.;
			si = (fabs(si) < fabs(sim1)) ? fabs(si) : fabs(sim1);
			(*dya)[i] = 2.*k*si;
			continue;
		}
		else (*dya)[i] = pi;

	}

	// Right boundary

	(*dya)[N-1] = (ya[N-1] - ya[N-2])/(xa[N-1] - xa[N-2]);

	// Compute the interp coefficients for the whole table

	for(i=0;i<N-1;i++){
		hi = ( xa[i+1] - xa[i] );
		si = ( ya[i+1] - ya[i] )/hi;

		(*a)[i] = ( (*dya)[i] + (*dya)[i+1] - 2.*si )/(hi*hi);
		(*b)[i] = ( 3.*si    - 2.*(*dya)[i] -(*dya)[i+1])/hi;
		(*c)[i] = (*dya)[i];
		(*d)[i] = ya[i];
	}

return 0;

}


int steffen_intp(double *xa, double *ya, double *a, double *b, double *c, double *d, int N, double x, double *y, double *dydx)
{

	int klo, khi, k, i=0;
	double h;


	klo = 0;
	khi = N-1;

	// Find position index in tab

	while(khi - klo > 1){
		k = (khi + klo) >> 1;
		if(xa[k] > x) khi = k;
		else klo = k;
	}
	//while ( x > xa[i]) i++;
	//klo = i-1;

	h     = x - xa[klo];
	*y    = d[klo] + h*(c[klo] + h*(b[klo] +h*a[klo]));
	*dydx = c[klo] + h*(2.*b[klo] + 3.*h*a[klo]);

	//printf("xa[%d] = %e, xa[%d] = %e, x = %e\n", klo, xa[klo], khi, xa[khi], x);

	return 0;
}
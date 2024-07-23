/* eos_tab1D_hot.c */
/* hgieg 08/2020 */


#include "bam.h"
#include "eos.h"

#define PR 0
#define TEST 1

struct {

	int N;
	double *Logr, *T, *Y;
	double *eps, *p, *cs2; // Hydrodynamics
	double *mun, *mue, *mup, *Abar, *Zbar, *nh; // Microphysics
	double *d2Ydr2, *d2pdr2, *d2epsdr2, *d2Tdr2, *d2cs2dr2;
	double *d2mundr2, *d2muedr2, *d2mupdr2, *d2Adr2, *d2Zdr2, *d2nhdr2;
	double mb;

} tab1d_hot;


void eos_load_tab1d_hot ()
{
	FILE *fp;

	int key = 1;
	// Read eos into tab1d_hot struct

	eos_read_tab1d_hot(Gets("eos_tab_file"), &(tab1d_hot.N), &(tab1d_hot.T), &(tab1d_hot.Logr), &(tab1d_hot.Y),
			&(tab1d_hot.p), &(tab1d_hot.eps), &(tab1d_hot.cs2), &(tab1d_hot.mun), &(tab1d_hot.mup),
			&(tab1d_hot.mue), &(tab1d_hot.Abar), &(tab1d_hot.Zbar), &(tab1d_hot.nh), &(tab1d_hot.mb), key);

	// Fill in second derivatives with respect to rho for cubic splines interpolation

	spline(tab1d_hot.Logr, tab1d_hot.T, tab1d_hot.N, &(tab1d_hot.d2Tdr2));
	spline(tab1d_hot.Logr, tab1d_hot.Y, tab1d_hot.N, &(tab1d_hot.d2Ydr2));
	spline(tab1d_hot.Logr, tab1d_hot.p, tab1d_hot.N, &(tab1d_hot.d2pdr2));
	spline(tab1d_hot.Logr, tab1d_hot.eps, tab1d_hot.N, &(tab1d_hot.d2epsdr2));
	spline(tab1d_hot.Logr, tab1d_hot.cs2, tab1d_hot.N, &(tab1d_hot.d2cs2dr2));
	spline(tab1d_hot.Logr, tab1d_hot.mun, tab1d_hot.N, &(tab1d_hot.d2mundr2));
	spline(tab1d_hot.Logr, tab1d_hot.mup, tab1d_hot.N, &(tab1d_hot.d2mupdr2));
	spline(tab1d_hot.Logr, tab1d_hot.mue, tab1d_hot.N, &(tab1d_hot.d2muedr2));
	spline(tab1d_hot.Logr, tab1d_hot.Abar, tab1d_hot.N, &(tab1d_hot.d2Adr2));
	spline(tab1d_hot.Logr, tab1d_hot.Zbar, tab1d_hot.N, &(tab1d_hot.d2Zdr2));
	spline(tab1d_hot.Logr, tab1d_hot.nh, tab1d_hot.N, &(tab1d_hot.d2nhdr2));

	if (TEST) {
		
		printf("N = %d\n", tab1d_hot.N);
		tU units;
		set_units(&units);

		int i, N=10000.;
		double rho, p, eps, T, Y, cs2, mun, mup, mue, A, Z, nh, tmp;
		double mn = 939.565413, dr = (tab1d_hot.Logr[tab1d_hot.N -1] - tab1d_hot.Logr[0])/(double) N;

		fp = fopen("EoS_tab1d_hot_test.txt", "w+");
		

		for(i=0;i < N;i++){
			rho = pow(10., (tab1d_hot.Logr[0]+i*dr));
			eos_tab1d_hot(&p, &cs2, &tmp, &tmp, &rho, &eps);
			eos_tab1d_hot_beta(&rho, &Y, &T);
			eos_tab1d_hot_micro(&rho, &Y, &T, &mun, &mup, &mue, &A, &Z, &nh);

			fprintf(fp, "%le %le %le %le %le\n", T, (rho/tab1d_hot.mb*units.Mdens_cgs*1e-39), Y,
				p*units.Energy_MeV/units.Volume_fm3, eps);
		}

		fclose(fp);
	}

}



// EoS validity range


void eos_tab1d_hot_validity (double *rho, double *Y, double *rho_min, double *rho_max, double *Y_min, double *Y_max, 
			     double *eps_min, double *eps_max, double *mb)
{

	int i;
	double Ymin, Ymax, epsmin, epsmax;

	Ymax = Ymin = tab1d_hot.Y[0];
	epsmax = epsmin = tab1d_hot.eps[0];

	for(i=1;i<tab1d_hot.N;i++){
		if(tab1d_hot.Y[i] >= Ymax) Ymax = tab1d_hot.Y[i];
		if(tab1d_hot.Y[i] <= Ymin) Ymin = tab1d_hot.Y[i];
		if(tab1d_hot.eps[i] >= epsmax) epsmax = tab1d_hot.eps[i];
		if(tab1d_hot.eps[i] <= epsmin) epsmin = tab1d_hot.eps[i];
	}

	*Y_min = Ymin; *Y_max = Ymax;
	*eps_min = epsmin; *eps_max = epsmax;
	*rho_min = pow(10.,tab1d_hot.Logr[0]);
	*rho_max = pow(10.,tab1d_hot.Logr[tab1d_hot.N-1]);
	
	*mb = tab1d_hot.mb;

}


// Extend validity

int eos_tab1d_hot_extend (double *rho, double *Y, double *eps)
{
	double rho_tmp, eps_tmp, rhomax, rhomin, dummy, epsl;

	rhomax = pow(10., tab1d_hot.Logr[tab1d_hot.N-1]);
	rhomin = pow(10., tab1d_hot.Logr[0]);

	// Force *rho into the validity range 

	rho_tmp = (*rho <= rhomax) ? *rho : rhomax;
	*rho 	= (rho_tmp >= rhomin) ? rho_tmp : rhomin;

	// Same with eps

	eos_tab1d_hot(&dummy, &dummy, &dummy, &dummy, rho, &epsl);

	*eps = epsl;

	//eps_tmp = (*eps <= epsl) ? *eps : epsl;
	//eps 	= (eps_tmp >= epsl) ? eps_tmp : epsl;

	return 0;

}

int eos_tab1d_hot (double *p, double *cs2, double *dpdrho,  double *dpdeps, double *rho, double *eps)
{

	double dpdLogr, depsdLogr, tmp;
	double Logr = log10(*rho);
	int i=0;

	if(!finite(Logr)){
		if (PR) printf("!finite Logr = %le for rho = %le\n", Logr, *rho);
		return 1;
	}

	splint(tab1d_hot.Logr, tab1d_hot.p, tab1d_hot.d2pdr2, tab1d_hot.N, Logr, p, &dpdLogr);
	splint(tab1d_hot.Logr, tab1d_hot.eps, tab1d_hot.d2epsdr2, tab1d_hot.N, Logr, eps, &depsdLogr);
	splint(tab1d_hot.Logr, tab1d_hot.cs2, tab1d_hot.d2cs2dr2, tab1d_hot.N, Logr, cs2, &tmp);

	/*while(Logr > tab1d_hot.Logr[i]) i++;
	if(i > 0){
		*p   = tab1d_hot.p[i-1] + (tab1d_hot.p[i] - tab1d_hot.p[i-1])/
					(tab1d_hot.Logr[i] - tab1d_hot.Logr[i-1])*(Logr - tab1d_hot.Logr[i-1]);
		*eps = tab1d_hot.eps[i-1] + (tab1d_hot.eps[i] - tab1d_hot.eps[i-1])/(tab1d_hot.Logr[i] - tab1d_hot.Logr[i-1])*(Logr - tab1d_hot.Logr[i-1]);
		*cs2 = tab1d_hot.cs2[i-1] + (tab1d_hot.cs2[i] - tab1d_hot.cs2[i-1])/(tab1d_hot.Logr[i] - tab1d_hot.Logr[i-1])*(Logr - tab1d_hot.Logr[i-1]);
		dpdLogr = (tab1d_hot.p[i] - tab1d_hot.p[i-1])/(tab1d_hot.Logr[i] - tab1d_hot.Logr[i-1]);
	}
	else {
		*p = tab1d_hot.p[i];
		*eps = tab1d_hot.eps[i];
		*cs2 = tab1d_hot.cs2[i];
		 dpdLogr = (tab1d_hot.p[i+1] - tab1d_hot.p[i])/(tab1d_hot.Logr[i+1] - tab1d_hot.Logr[i]);
	}

	*/

	*dpdrho = 1./(*rho*log(10.))*dpdLogr;
	*dpdeps = 0.;

	if (PR) printf("rho = %e  -> eps = %e  p = %e, dpdr = %e  dpdeps = %e -> cs2 = %e\n", *rho, *eps, *p, *dpdrho, *dpdeps, *cs2);
        if (!finite(*rho)|| !finite(*eps) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e, cs2 = %e\n", *rho, *eps, *p, *cs2);
                return 1;
        }
        return 0;


}

int eos_tab1d_hot_micro (double *rho, double *T, double *Y, double *mun, double *mup, double *mue, double *A, double *Z, double *nh)
{

	double Logr = log10(*rho), tmp;

	splint(tab1d_hot.Logr, tab1d_hot.mun, tab1d_hot.d2mundr2, tab1d_hot.N, Logr, mun, &tmp);
	splint(tab1d_hot.Logr, tab1d_hot.mup, tab1d_hot.d2mupdr2, tab1d_hot.N, Logr, mup, &tmp);
	splint(tab1d_hot.Logr, tab1d_hot.mue, tab1d_hot.d2muedr2, tab1d_hot.N, Logr, mue, &tmp);
	splint(tab1d_hot.Logr, tab1d_hot.Abar, tab1d_hot.d2Adr2, tab1d_hot.N, Logr, A, &tmp);
	splint(tab1d_hot.Logr, tab1d_hot.Zbar, tab1d_hot.d2Zdr2, tab1d_hot.N, Logr, Z, &tmp);
	splint(tab1d_hot.Logr, tab1d_hot.nh, tab1d_hot.d2nhdr2, tab1d_hot.N, Logr, nh, &tmp);

	if (PR) printf("rho = %e -> T = %e, Y = %e, mu_n = %e, mu_p = %e, mu_e = %e, A = %e, Z = %e, n_h = %e\n",
			*rho, *T, *Y, *mun, *mup, *mue, *A, *Z, *nh);

	if(!finite(*rho) || !finite(*T) || !finite(*Y) || !finite(*mun) || !finite(*mup) || !finite(*mue)
		|| !finite(*A) || !finite(*Z) || !finite(*nh))
	{
		if (PR) printf("!finite: rho = %e -> T = %e, Y = %e, mu_n = %e, mu_p = %e, mu_e = %e, A = %e, Z = %e, n_h = %e\n",
				*rho, *T, *Y, *mun, *mup, *mue, *A, *Z, *nh);
		return 1;
	}

	return 0; 

}

int eos_tab1d_hot_beta (double *rho, double *Y, double *T)
{

	double dummy;
	double Logr = log10(*rho);
	int i=0;

	splint(tab1d_hot.Logr, tab1d_hot.Y, tab1d_hot.d2Ydr2, tab1d_hot.N, Logr, Y, &dummy);
	splint(tab1d_hot.Logr, tab1d_hot.T, tab1d_hot.d2Tdr2, tab1d_hot.N, Logr, T, &dummy);

	/*while(Logr > tab1d_hot.Logr[i]) i++;
        if(i > 0){
                *Y   = tab1d_hot.Y[i-1] + (tab1d_hot.Y[i] - tab1d_hot.Y[i-1])/(tab1d_hot.Logr[i] - tab1d_hot.Logr[i-1])*(Logr - tab1d_hot.Logr[i-1]);
                *T = tab1d_hot.T[i-1] + (tab1d_hot.T[i] - tab1d_hot.T[i-1])/(tab1d_hot.Logr[i] - tab1d_hot.Logr[i-1])*(Logr - tab1d_hot.Logr[i-1]);
        }
        else {
                *Y = tab1d_hot.Y[i];
                *T = tab1d_hot.T[i];
        }*/


	if(!finite(*Y) || !finite(*T)) return 1;

	return 0;
}

/* eos_tab1D_cold.h */
/* hgieg 03/2020 */

#include "bam.h"
#include "eos.h"


#define PR 0
#define TEST 1

struct {
        int N;
        double *Logr, *T, *Y;
        double *eps, *p;
	double *d2Ydr2, *d2pdr2, *d2epsdr2, *d2Tdr2;
	double *d2rdp2, *d2epsdp2, *d2Tdp2, *d2Ydp2;
	double mb;
} tab1d;

struct {
	double *ap, *bp, *cp, *dp;
	double *aeps, *beps, *ceps, *deps;
	double *aY, *bY, *cY, *dY;
	double *aT, *bT, *cT, *dT;
	double *ar, *br, *cr, *dr;
	double *aep, *bep, *cep, *dep;
	double *dpdLogr, *depsdLogr;
} stf;

/* Read tab1D: 

	columns: T [MeV] | nB [fm^-3] | Ye | p [MeV fm^-3] | internal energy/baryon [MeV]
                        								*/

void eos_read_tab1d(char *fname, int *n, double **x1, double **x2, double **x3, double **v1, double **v2, double *mb)
{
  int N, i, j, P;
  char line[1024];

  tU units;
  set_units(&units);

  printf("\nReading tab1D EoS: \n%s\n", fname);

  // Read for all processors

  for (P=0; P<=bampi_size();P++){
    
    printf(" Reading on proc%d/%d\n", P, bampi_size());

    if(bampi_rank()==P){
	
        FILE *f;
        f = fopen(fname, "r");

        //Count lines
	N = 0;

        while(fgets(line,1024,f)){
                N = N+1;
        }

        rewind(f);
        double **t_temp = malloc(N*sizeof(double*));
        for(i=0;i<N;i++) t_temp[i] = malloc(5*sizeof(double));

        for(i=0;i<N;i++){
                for(j=0;j<5;j++){
                        fscanf(f, "%le", &t_temp[i][j]);
			if (PR) printf("%le ", t_temp[i][j]);
                }
		if (PR) printf("\n");
        }

	// Find the baryonic mass constant mb
	
	*mb = t_temp[0][4];

	for(i=1;i<N;i++){
		if(t_temp[i][4] <= *mb) *mb = t_temp[i][4];
	}

	*mb = *mb/units.Energy_MeV*units.Mass_cgs; //mb in g

        //Setting and filling the vectors

        (*x1) = (double *) malloc(N*sizeof(double));
        (*x2) = (double *) malloc(N*sizeof(double));
 	(*x3) = (double *) malloc(N*sizeof(double));
	(*v1) = (double *) malloc(N*sizeof(double));
	(*v2) = (double *) malloc(N*sizeof(double));

        for(i=0; i<N; i++){
		(*x1)[i] = t_temp[i][0]; 	// temperature in MeV
                (*x2)[i] = log10(t_temp[i][1]*(*mb)*1e+39/units.Mdens_cgs); //log10(rest-mass density)
                (*x3)[i] = t_temp[i][2]; //charge fraction
		(*v1)[i] = t_temp[i][3]/units.Energy_MeV*units.Volume_fm3; //pressure
		(*v2)[i] = (t_temp[i][4]*units.Mass_cgs/units.Energy_MeV)/(*mb) - 1.0; //specific internal energy eps (dimensionless)
        }



	fclose(f);
	for(i=0;i<N;i++) free(t_temp[i]);
  	free(t_temp);

    }

  }

  *n = N;
  printf("N = %d\n", *n);
  printf("mb = %le\n", *mb);
  printf("\nEoS tab1D loaded successfully from file \n%s\n", fname);
}


void eos_load_tab1d_cold ()
{
	double *dummy;
	int key;

	//Read eos tab 1D into the struct tab1d

	eos_read_tab1d(Gets("eos_tab_file"), &(tab1d.N), &(tab1d.T), &(tab1d.Logr), &(tab1d.Y), &(tab1d.p), &(tab1d.eps), &(tab1d.mb));

	//Fill in the second derivatives of Y, p and eps with respect to rho for cubic splines interpolation
	if(Getv("eos_interp", "cspline")){
		spline(tab1d.Logr, tab1d.Y, tab1d.N, &(tab1d.d2Ydr2));
		spline(tab1d.Logr, tab1d.p, tab1d.N, &(tab1d.d2pdr2));
		spline(tab1d.Logr, tab1d.eps, tab1d.N, &(tab1d.d2epsdr2));
		spline(tab1d.Logr, tab1d.T, tab1d.N, &(tab1d.d2Tdr2));

		spline(tab1d.p, tab1d.Logr, tab1d.N, &(tab1d.d2rdp2));
		spline(tab1d.p, tab1d.eps, tab1d.N, &(tab1d.d2epsdp2));
		if(TEST) key = 0;
	}

	else if(Getv("eos_interp", "steffen")){
		steffen_init(tab1d.Logr, tab1d.p, tab1d.N, &(stf.dpdLogr), &(stf.ap), &(stf.bp), &(stf.cp), &(stf.dp));
		steffen_init(tab1d.Logr, tab1d.eps, tab1d.N, &(stf.depsdLogr), &(stf.aeps), &(stf.beps), &(stf.ceps), &(stf.deps));
		steffen_init(tab1d.Logr, tab1d.T, tab1d.N, &dummy, &(stf.aT), &(stf.bT), &(stf.cT), &(stf.dT));
        	steffen_init(tab1d.Logr, tab1d.Y, tab1d.N, &dummy, &(stf.aY), &(stf.bY), &(stf.cY), &(stf.dY));
		
		//Fill constants for interp rho(p), eps(p)

		steffen_init(tab1d.p, tab1d.Logr, tab1d.N, &dummy, &(stf.ar), &(stf.br), &(stf.cr), &(stf.dr));
		steffen_init(tab1d.p, tab1d.eps, tab1d.N, &dummy, &(stf.aep), &(stf.bep), &(stf.cep), &(stf.dep));
		if(TEST) key = 1;
        }


	if(TEST){

		tU units;
		set_units(&units);

		double rho, rhomax, eps, p, dpdr, dpdeps, cs2, drho;
		int i, N = 10000;

		FILE *fp;

		fp = fopen("EoS_test.txt", "w+");
		//rhomax = pow(10.,tab1d.Logr[tab1d.N-1]);
		//drho = (tab1d.Logr[tab1d.N-1] - tab1d.Logr[0])/(N-1);
		drho = (tab1d.Logr[tab1d.N-1] - tab1d.Logr[0])/(double) (N-1);

		for(i=0;i<N;i++){
			rho = tab1d.Logr[0] + i*drho;
			//p = tab1d.p[0]+i*drho;
			rho = pow(10.,rho);
			if(key) eos_tab1d_stf(&p, &cs2, &dpdr,&dpdeps,&rho,&eps);
			//eos_tab1d_stf_p(&p, &cs2, &dpdr, &dpdeps, &rho, &eps);
			//eos_tab1d_lin(&p,&cs2,&dpdr,&dpdeps,&rho,&eps);
			eps = tab1d.mb*(eps+1.);

			rho = rho*units.Mdens_cgs/tab1d.mb*1e-39;
			p = p*units.Energy_MeV/units.Volume_fm3;
			eps = eps*units.Energy_MeV;
			fprintf(fp, "%le %le %le %le %le %le\n", rho, p, eps, cs2, dpdr, dpdeps);
		}

		fclose(fp);
	}

}

// Cubic splines interpolation routines are defined in eos_tab3D_gieg.c


// EoS validity range

void eos_tab1d_validity (double *rho, double *Y, double *rho_min, double *rho_max, double *Y_min, double *Y_max, double *eps_min, double *eps_max,
			 double *mb)
{

	double tmp, Ymin, Ymax, epsmin, epsmax;
	int i=0;
	
	Ymax = Ymin = tab1d.Y[0];
	epsmax = epsmin = tab1d.eps[0];

	for(i=1;i<tab1d.N;i++){
		if(tab1d.Y[i] >= Ymax) Ymax = tab1d.Y[i];
		if(tab1d.Y[i] <= Ymin) Ymin = tab1d.Y[i];
		if(tab1d.eps[i] >= epsmax) epsmax = tab1d.eps[i];
		if(tab1d.eps[i] <= epsmin) epsmin = tab1d.eps[i];
	}

	*Y_min = Ymin;
	*Y_max = Ymax;

	*rho_min = pow(10,tab1d.Logr[0]);
	*rho_max = pow(10,tab1d.Logr[tab1d.N-1]);

	*eps_min = epsmin;
	*eps_max = epsmax;

	*mb = tab1d.mb;
}

// EoS extend validity

int eos_tab1d_extend (double *rho, double *Y, double *eps)
{

	double rho_tmp, eps_tmp, epsl, rhomax, rhomin, dummy;

	rhomax = pow(10,tab1d.Logr[tab1d.N-1]);
	rhomin = pow(10,tab1d.Logr[0]);

	//Force *rho into the validity range

	rho_tmp = (*rho <= rhomax) ? *rho : rhomax;
	*rho    = (rho_tmp >= rhomin) ? rho_tmp: rhomin;

	//Force *eps into the validity range

	eos_tab1d (&dummy, &dummy, &dummy, &dummy, rho, &epsl);

	*eps = epsl;
	return 0;
}


//EoS wrapping functions

int eos_tab1d (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps)
{
	double epsl, depsdr, dpdLogr, depsdLogr;
	double Logr = log10(*rho);
	int i;
	if(!finite(Logr)) {
		if(PR) printf("!finite Logr = %e for rho = %e\n", Logr, *rho);
		return 1;
	}

	splint(tab1d.Logr, tab1d.p, tab1d.d2pdr2, tab1d.N, Logr, p, &dpdLogr);
	splint(tab1d.Logr, tab1d.eps, tab1d.d2epsdr2, tab1d.N, Logr, &epsl, &depsdLogr);
	*eps = epsl;
	*dpdeps = 0.;
	*dpdrho = 1./(*rho*log(10.))*dpdLogr;
	*cs2 = eos_cs2_rep(*rho, epsl, *p, *dpdrho, *dpdeps);
	if(*cs2 < 0.) *cs2 = 0.;
	if (PR) printf("rho = %e  -> eps = %e  p = %e, dpdr = %e  dpdeps = %e -> cs2 = %e\n", *rho, *eps, *p, *dpdrho, *dpdeps, *cs2);
	if (!finite(*rho)|| !finite(*eps) || !finite(*cs2)) {
		if(PR) printf("!finite, rho = %e, eps = %e, p = %e\n", *rho, *eps, *p);
		return 1;
	}
	return 0;
}

int eos_tab1d_beta (double *rho, double *Y, double *T)
{
        double dummy;
        double Logr = log10(*rho);


        splint(tab1d.Logr, tab1d.Y, tab1d.d2Ydr2, tab1d.N, Logr, Y, &dummy);
	splint(tab1d.Logr, tab1d.T, tab1d.d2Tdr2, tab1d.N, Logr, T, &dummy);
        if(!finite(*Y)||!finite(*T)) return 1;

        return 0;
}

int eos_tab1d_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e)
{

	double epsl, Logr, dLogrdp, depsdp, dpdrho, dpdeps;

        if(!finite(*p)) {
                if(PR) printf("!finite p = %e\n", *p);
                return 1;
        }

        splint(tab1d.p, tab1d.Logr, tab1d.d2rdp2, tab1d.N, *p, &Logr, &dLogrdp);
        splint(tab1d.p, tab1d.eps, tab1d.d2epsdp2, tab1d.N, *p, &epsl, &depsdp);

        *rho = pow(10.,Logr);
        *e = *rho*(1.+epsl);
        *drhodp = *rho*log(10.)*dLogrdp;
        *dedp   = (1.+epsl)*(*drhodp) + (*rho)*depsdp;

        dpdeps = 0.;
        dpdrho = 1./(*drhodp);
        *cs2 = eos_cs2_rep(*rho, epsl, *p, dpdrho, dpdeps);
        if (PR) printf("p = %e  -> rho = %e  e = %e, drhodp = %e  dedp = %e -> cs2 = %e\n", *p, *rho, *e, *drhodp, *dedp, *cs2);
        if (!finite(*rho)|| !finite(*e) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e\n", *rho, *e, *p);
                return 1;
        }
        return 0;

}


// tab1D linear extrapolation in p for rho > EOS.rhomax

int eos_tab1d_extrapolate (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps)
{

  	double epsl, depsdr, dpdLogr, depsdLogr, k, pmax, epsmax;
	double pimo, epsimo, gamma, ccc;
	double rhomax = pow(10., tab1d.Logr[tab1d.N - 1]);
	double rhoimo = pow(10., tab1d.Logr[tab1d.N - 2]);
	double Logr = log10(*rho);
        int i;
        if(!finite(Logr)) {
                if(PR) printf("!finite Logr = %e for rho = %e\n", Logr, *rho);
                return 1;
        }

        splint(tab1d.Logr, tab1d.p, tab1d.d2pdr2, tab1d.N, tab1d.Logr[tab1d.N-1], &pmax, &dpdLogr);
	splint(tab1d.Logr, tab1d.p, tab1d.d2pdr2, tab1d.N, tab1d.Logr[tab1d.N-2], &pimo, &dpdLogr);

	gamma = log(pmax/pimo)/log(rhomax/rhoimo);
	k     = pmax/pow(rhomax,gamma);

        splint(tab1d.Logr, tab1d.eps, tab1d.d2epsdr2, tab1d.N, tab1d.Logr[tab1d.N-1], &epsmax, &depsdLogr);
        *dpdrho = gamma*k*pow(*rho, (gamma - 1.));
        *dpdeps = 0.;

        *p = k*pow(*rho, gamma);
	ccc = epsmax - k*pow(rhomax, (gamma - 1.))/(gamma - 1.);
	*eps = k*pow(*rho, (gamma - 1.))/(gamma - 1.) + ccc;
        //k = epsmax - *dpdrho*(1.+log(rhomax)) + pmax/rhomax;
        //eps = k + (*dpdrho * rhomax - pmax)/(*rho) + *dpdrho * log(*rho);
        *cs2 = *dpdrho/(1.+*eps+*p/(*rho));

        if (PR) printf("rho = %e  -> eps = %e  p = %e, dpdr = %e  dpdeps = %e -> cs2 = %e\n", *rho, *eps, *p, *dpdrho, *dpdeps, *cs2);
        if (!finite(*rho)|| !finite(*eps) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e\n", *rho, *eps, *p);
		return 1;
	}

	return 0;

}

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


// Compute usual 1D EOS as function of rho -------> EOS(rho)

int eos_tab1d_stf (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps)
{
	double epsl, depsdr, dpdLogr, depsdLogr;
        double Logr = log10(*rho);
        int i;
        if(!finite(Logr)) {
                if(PR) printf("!finite Logr = %e for rho = %e\n", Logr, *rho);
                return 1;
        }

        steffen_intp(tab1d.Logr, tab1d.p, stf.ap, stf.bp, stf.cp, stf.dp, tab1d.N, Logr, p, &dpdLogr);
        steffen_intp(tab1d.Logr, tab1d.eps, stf.aeps, stf.beps, stf.ceps, stf.deps, tab1d.N, Logr, &epsl, &depsdLogr);
        *eps = epsl;
        *dpdeps = 0.;
        *dpdrho = 1./(*rho*log(10.))*dpdLogr;
        *cs2 = eos_cs2_rep(*rho, epsl, *p, *dpdrho, *dpdeps);
        if (PR) printf("rho = %e  -> eps = %e  p = %e, dpdr = %e  dpdeps = %e -> cs2 = %e\n", *rho, *eps, *p, *dpdrho, *dpdeps, *cs2);
        if (!finite(*rho)|| !finite(*eps) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e\n", *rho, *eps, *p);
                return 1;
        }
        return 0;
}


// Compute 1D EOS as function of p ---> EOS(p)
int eos_tab1d_stf_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e)
{

	double epsl, Logr, dLogrdp, depsdp, dpdrho, dpdeps;

	if(!finite(*p)) {
                if(PR) printf("!finite p = %e\n", *p);
                return 1;
        }

	steffen_intp(tab1d.p, tab1d.Logr, stf.ar, stf.br, stf.cr, stf.dr, tab1d.N, *p, &Logr, &dLogrdp);
	steffen_intp(tab1d.p, tab1d.eps, stf.aep, stf.bep, stf.cep, stf.dep, tab1d.N, *p, &epsl, &depsdp);

	*rho = pow(10.,Logr);
	*e = *rho*(1.+epsl);
	*drhodp = *rho*log(10.)*dLogrdp;
	*dedp   = (1.+epsl)*(*drhodp) + (*rho)*depsdp;
	
	dpdeps = 0.;
	dpdrho = 1./(*drhodp);
	*cs2 = eos_cs2_rep(*rho, epsl, *p, dpdrho, dpdeps);
	if (PR) printf("p = %e  -> rho = %e  e = %e, drhodp = %e  dedp = %e -> cs2 = %e\n", *p, *rho, *e, *drhodp, *dedp, *cs2);
        if (!finite(*rho)|| !finite(*e) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e\n", *rho, *e, *p);
                return 1;
        }
        return 0;

}

int eos_tab1d_beta_stf (double *rho, double *Y, double *T)
{
        double dummy;
        double Logr = log10(*rho);


        steffen_intp(tab1d.Logr, tab1d.Y, stf.aY, stf.bY, stf.cY, stf.dY, tab1d.N, Logr, Y, &dummy);
        steffen_intp(tab1d.Logr, tab1d.T, stf.aT, stf.bT, stf.cT, stf.dT, tab1d.N, Logr, T, &dummy);
        if(!finite(*Y)||!finite(*T)) return 1;

        return 0;
}


int eos_tab1d_extend_stf (double *rho, double *Y, double *eps)
{

        double rho_tmp, eps_tmp, epsl, rhomax, rhomin, dummy;

        rhomax = pow(10,tab1d.Logr[tab1d.N-1]);
        rhomin = pow(10,tab1d.Logr[0]);

        //Force *rho into the validity range

        rho_tmp = (*rho <= rhomax) ? *rho : rhomax;
        *rho    = (rho_tmp >= rhomin) ? rho_tmp: rhomin;

        //Force *eps into the validity range

        eos_tab1d_stf (&dummy, &dummy, &dummy, &dummy, rho, &epsl);

        *eps = epsl;
	return 0;
}


//=========================================== Linear Interpolation Routines =====================================//


int eos_tab1d_lin (double *p, double *cs2, double *dpdrho,  double *dpdeps, double *rho, double *eps)
{

	double dpdLogr, depsdLogr, tmp;
	double Logr = log10(*rho);
	int i=0;

	if(!finite(Logr)){
		if (PR) printf("!finite Logr = %le for rho = %le\n", Logr, *rho);
		return 1;
	}

	//splint(tab1d_hot.Logr, tab1d_hot.p, tab1d_hot.d2pdr2, tab1d_hot.N, Logr, p, &dpdLogr);
	//splint(tab1d_hot.Logr, tab1d_hot.eps, tab1d_hot.d2epsdr2, tab1d_hot.N, Logr, eps, &depsdLogr);
	//splint(tab1d_hot.Logr, tab1d_hot.cs2, tab1d_hot.d2cs2dr2, tab1d_hot.N, Logr, cs2, &tmp);

	while(Logr > tab1d.Logr[i]) i++;
	if(i > 0){
		*p   = tab1d.p[i-1] + (tab1d.p[i] - tab1d.p[i-1])/
					(tab1d.Logr[i] - tab1d.Logr[i-1])*(Logr - tab1d.Logr[i-1]);
		*eps = tab1d.eps[i-1] + (tab1d.eps[i] - tab1d.eps[i-1])/(tab1d.Logr[i] - tab1d.Logr[i-1])*(Logr - tab1d.Logr[i-1]);
		dpdLogr = (tab1d.p[i] - tab1d.p[i-1])/(tab1d.Logr[i] - tab1d.Logr[i-1]);
	}
	else {
		*p = tab1d.p[i];
		*eps = tab1d.eps[i];
		 dpdLogr = (tab1d.p[i+1] - tab1d.p[i])/(tab1d.Logr[i+1] - tab1d.Logr[i]);
	}

	*dpdrho = 1./(*rho*log(10.))*dpdLogr;
	*dpdeps = 0.;
	*cs2 = eos_cs2_rep(*rho, *eps, *p, *dpdrho, *dpdeps);


	if (PR) printf("rho = %e  -> eps = %e  p = %e, dpdr = %e  dpdeps = %e -> cs2 = %e\n", *rho, *eps, *p, *dpdrho, *dpdeps, *cs2);
        if (!finite(*rho)|| !finite(*eps) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e, cs2 = %e\n", *rho, *eps, *p, *cs2);
                return 1;
        }
        return 0;


}


int eos_tab1d_beta_lin (double *rho, double *Y, double *T)
{

	double dummy;
	double Logr = log10(*rho);
	int i=0;

	//splint(tab1d_hot.Logr, tab1d_hot.Y, tab1d_hot.d2Ydr2, tab1d_hot.N, Logr, Y, &dummy);
	//splint(tab1d_hot.Logr, tab1d_hot.T, tab1d_hot.d2Tdr2, tab1d_hot.N, Logr, T, &dummy);

	while(Logr > tab1d.Logr[i]) i++;
        if(i > 0){
                *Y   = tab1d.Y[i-1] + (tab1d.Y[i] - tab1d.Y[i-1])/(tab1d.Logr[i] - tab1d.Logr[i-1])*(Logr - tab1d.Logr[i-1]);
                *T = tab1d.T[i-1] + (tab1d.T[i] - tab1d.T[i-1])/(tab1d.Logr[i] - tab1d.Logr[i-1])*(Logr - tab1d.Logr[i-1]);
        }
        else {
                *Y = tab1d.Y[i];
                *T = tab1d.T[i];
        }


	if(!finite(*Y) || !finite(*T)) return 1;

	return 0;
}

int eos_tab1d_lin_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e)
{

	double epsl, Logr, dLogrdp, depsdp, dpdrho, dpdeps;
	int i;

	if(!finite(*p)) {
                if(PR) printf("!finite p = %e\n", *p);
                return 1;
        }

	while(*p > tab1d.p[i]) i++;
	if(i > 0){
		Logr   = tab1d.Logr[i-1] + (tab1d.Logr[i] - tab1d.Logr[i-1])/
					(tab1d.p[i] - tab1d.p[i-1])*(*p - tab1d.p[i-1]);
		epsl = tab1d.eps[i-1] + (tab1d.eps[i] - tab1d.eps[i-1])/
				(tab1d.p[i] - tab1d.p[i-1])*(*p - tab1d.p[i-1]);
		dLogrdp = (tab1d.Logr[i] - tab1d.Logr[i-1])/(tab1d.p[i] - tab1d.p[i-1]);
		depsdp  = (tab1d.eps[i] - tab1d.eps[i-1])/(tab1d.p[i] - tab1d.p[i-1]);
	}
	else {
		Logr = tab1d.Logr[i];
		epsl = tab1d.eps[i];
		dLogrdp = (tab1d.Logr[i+1] - tab1d.Logr[i])/(tab1d.p[i+1] - tab1d.p[i]);
		depsdp = (tab1d.eps[i+1] - tab1d.eps[i])/(tab1d.p[i+1] - tab1d.p[i]);
	}

	*rho = pow(10.,Logr);
	*e = *rho*(1.+epsl);
	*drhodp = *rho*log(10.)*dLogrdp;
	*dedp   = (1.+epsl)*(*drhodp) + (*rho)*depsdp;
	
	dpdeps = 1./depsdp;
	dpdrho = 1./(*drhodp);
	*cs2 = eos_cs2_rep(*rho, epsl, *p, dpdrho, dpdeps);
	if (PR) printf("p = %e  -> rho = %e  e = %e, drhodp = %e  dedp = %e -> cs2 = %e\n", *p, *rho, *e, *drhodp, *dedp, *cs2);
        if (!finite(*rho)|| !finite(*e) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e\n", *rho, *e, *p);
                return 1;
        }
        return 0;

}

int eos_tab1d_extend_lin (double *rho, double *Y, double *eps)
{

        double rho_tmp, eps_tmp, epsl, rhomax, rhomin, dummy;

        rhomax = pow(10,tab1d.Logr[tab1d.N-1]);
        rhomin = pow(10,tab1d.Logr[0]);

        //Force *rho into the validity range

        rho_tmp = (*rho <= rhomax) ? *rho : rhomax;
        *rho    = (rho_tmp >= rhomin) ? rho_tmp: rhomin;

        //Force *eps into the validity range

        eos_tab1d_lin(&dummy, &dummy, &dummy, &dummy, rho, &epsl);

        *eps = epsl;
	return 0;
}


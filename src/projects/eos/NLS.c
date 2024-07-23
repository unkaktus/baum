/*  Neutrinos Leakage Scheme */ 
/*	hg 04/20	     */


#include "bam.h"
#include "eos.h"

#define PR 0
#define LOC 0
#define DIF 0
#define TEST 0

/******************** Terminology ***********************/
/*						     	*/
/* neutrino species:	nu_e (e): 0	          	*/
/*		      	nu_a (a): 1		     	*/
/*		      	nu_x (x): 2	          	*/
/*						     	*/
/*						     	*/
/* grid indices:	0: central		     	*/
/*			1,2: -x, +x			*/
/*			3,4: -y, +y			*/
/*			5,6: -z, +z			*/
/*							*/
/* direction indices:	0: x				*/
/*			1: y				*/
/*			2: z				*/
/*							*/
/********************************************************/


double const me 	= 0.510998910;				// electron mass in MeV/c^2
double const sigma_0 	= 1.705e-44;				// cm^2
double const alp   	= 1.23;					// dimensionless
double const Cv		= 0.962;				// dimensionless normalized vector const
double const Ca		= 0.5;					// dimensionless normalized axial const
double const fp		= 5.565e-02;				// dimensionless plasmon frequency
double const fsc	= 1.0/137.036;				// fine structure constant

/************************* NLS Wrapper *************************************/

// Receive grid spacings, covariant 3-metric diagonal components (-,0,+) of a given cell,
// the primitives *rho, *T, *Y from a cell and it's six neighbours and outputs the effective
// emission rates, the NLS source terms and the optical depth


void NLS_compute_sources (double dx, double dy, double dz, double *gxx, double *gyy, double *gzz,
			  double *rho, double *T, double *Y, double **Q_eff, double **R_eff, double *sQ,
			  double *sR, double **leak_tau, int keyx, int keyy, int keyz)
{
	tU units;
	set_units(&units);
	double *nn, *np, *nh, *eta_nue, *A, *Z, *eta_e;
	double eta_np, eta_pn, nb;
	double *Q_loc, *R_loc, *Q_diff, *R_diff;
	int i;

	//Allocate EoS vars
	nn = malloc(7*sizeof(double));
	nh = malloc(7*sizeof(double));
	np = malloc(7*sizeof(double));
	A  = malloc(7*sizeof(double));
	Z  = malloc(7*sizeof(double));
	eta_nue = malloc(7*sizeof(double));
	eta_e   = malloc(7*sizeof(double));

	//Allocate emission rates
	Q_loc = malloc(3*sizeof(double));
	R_loc = malloc(3*sizeof(double));
	Q_diff = malloc(3*sizeof(double));
	R_diff = malloc(3*sizeof(double));
	//(*Q_eff) = (double *) malloc(3*sizeof(double));
	//(*R_eff) = (double *) malloc(3*sizeof(double));

	get_EoS_vars (rho, T, Y, &nb, &nn, &np, &nh,
		&eta_np, &eta_pn, &(eta_e), &(eta_nue), &A, &Z);

	compute_local_emission (nb, T[0], eta_np, eta_pn, eta_e[0], eta_nue[0],
			 &(Q_loc), &(R_loc));
	compute_diff_emission (nb, T[0], nn, np, nh, A, Z, eta_nue, eta_e, gxx,
			gyy, gzz, dx, dy, dz, &(Q_diff), &(R_diff), leak_tau, keyx, keyy, keyz);

	if(CheckForNANandINF(12, R_loc[0], R_diff[0], Q_loc[0], Q_diff[0],
				 R_loc[1], R_diff[1], Q_loc[1], Q_diff[1],
				 R_loc[1], R_diff[1], Q_loc[1], Q_diff[1]))
	{
                        printf("R_diff[0] = %le, R_loc[0] = %le, Q_diff[0] = %le, Q_loc[0] = %le\n", R_diff[0], R_loc[0], 
                                                           Q_diff[0], Q_loc[0]);
                        printf("R_diff[1] = %le, R_loc[1] = %le, Q_diff[1] = %le, Q_loc[1] = %le\n", R_diff[1], R_loc[1], 
                                                           Q_diff[1], Q_loc[1]);
                        printf("R_diff[2] = %le, R_loc[2] = %le, Q_diff[2] = %le, Q_loc[2] = %le\n", R_diff[2], R_loc[2], 
                                                           Q_diff[2], Q_loc[2]);
			printf("nb = %le, T = %le, Y = %le, eta_e = %le\n", nb*1e-39, T[0], Y[0], eta_e[0]);
        }

	// Compute effective emission rates
	for(i=0;i<3;i++){
		if (R_diff[i]==0.0) (*R_eff)[i] = 0.0;
		else (*R_eff)[i] = R_loc[i]/(1.+R_loc[i]/(R_diff[i]));
		if (Q_diff[i]==0.0) (*Q_eff)[i] = 0.0;
		else (*Q_eff)[i] = Q_loc[i]/(1.+Q_loc[i]/(Q_diff[i]));
		// Convert units: [Q] = MeV/s -> erg/s, [R] = 1/s
		(*Q_eff)[i] = (*Q_eff)[i]/units.Energy_MeV*units.Energy_cgs;
 	}


 	// Compute source terms of NLS
	
	*sR = ((*R_eff)[1] - (*R_eff)[0])*units.Time_cgs; // in units of G = c = Msun = 1;
	*sQ = ((*Q_eff)[0] + (*Q_eff)[1] + (*Q_eff)[2])/units.Energy_cgs*units.Time_cgs; // in units of G = c = Msun = 1;

        if(CheckForNANandINF(2,*sR,*sQ)){
                        //printf("m = %d, n = %d, o = %d\n", level->box[bi]->m,
                          //      level->box[bi]->n, level->box[bi]->o);
                        printf("gxx[0] = %le, gxx[1] = %le, gxx[2] = %le\n", gxx[0], gxx[1], gxx[2]);
                        printf("gyy[0] = %le, gyy[1] = %le, gyy[2] = %le\n", gyy[0], gyy[1], gyy[2]);
                        printf("gzz[0] = %le, gzz[1] = %le, gzz[2] = %le\n", gzz[0], gzz[1], gzz[2]);
                        for(i=0;i<7;i++){
                                printf("i = %d, rho = %le, T = %le, Y = %le\n", i, rho[i], T[i], Y[i]);
                        }
	}
	// Deallocate arrays
	free(nn);
	free(nh);
	free(np);
	free(A);
	free(Z);
	free(eta_nue);
	free(eta_e);
	free(Q_loc);
	free(R_loc);
	free(Q_diff);
	free(R_diff);
}

/******************** EoS call for parameters ******************************/

void get_EoS_vars (double *rho, double *T, double *Y, double *nb, 
		   double **nn, double **np, double **nh,
                   double *eta_np, double *eta_pn, double **eta_e, 
		   double **eta_nue, double **A, double **Z)
{
	tU units;
	set_units(&units);

	double mu_p, mu_n, mu_e, eta_hat, mu_p0, mu_n0, eta_e0, nb_0;
	int i;

	for(i=0;i<7;i++){
	 *nb = rho[i]*units.Mdens_cgs/EOS.mb;		 			//nb in cm^-3
	 (*np)[i] = Y[i]*(*nb);
	 (*nn)[i] = *nb - (*np)[i];
	 EOS.micro(&(rho[i]), &(T[i]), &(Y[i]), &mu_n, &mu_p, &mu_e, &((*A)[i]), &((*Z)[i]), &((*nh)[i]));
	 
         (*eta_e)[i] = mu_e/T[i];
	 (*eta_nue)[i] = (*eta_e)[i] - mu_n/T[i] + mu_p/T[i]; //electron neutrinos chemical potential come from beta-equilibrium assumption

         if(i==0){
 		mu_p0 = mu_p;
		mu_n0 = mu_n;
		nb_0 = *nb;
	 }

	}
	*nb = nb_0;

       //Extremely important note in order to avoid confusion: here we follow Galeazzi et al. (2013) notation
       //for the nucleons blocking factors bellow. Despite Rosswog & Liebendörfer (2003)
       //interchange the np <-> pn indices, the definitions are still consistent for the application.

	eta_hat = mu_n0/T[0] - mu_p0/T[0] -1.2933/T[0];
	*eta_np = ((*nn)[0]-(*np)[0])/(exp(eta_hat)-1.);
	*eta_np = DMAX(*eta_np, (*np)[0]);
	*eta_pn = ((*np)[0]-(*nn)[0])/(exp(-eta_hat)-1.);
	*eta_pn = DMAX(*eta_pn, (*nn)[0]);

	if(rho[0]<(2e+12/units.Mdens_cgs)){
                *eta_np = (*np)[0];
                *eta_pn = (*nn)[0];
        }

}





/**************************** Local Emission Rates ***************************/

void compute_local_emission (double nb, double T, double eta_np, double eta_pn, double eta_e, 
			     double eta_nue, double **Q_loc, double **R_loc)
{

/*Included processes: 1. Electrons and antielectrons capture in free nucleons
		      2. Electron-positrons pair annihilation
		      3. Transversal plasmon decay */
	
	tU units;
	set_units(&units);
	double h = units.hplanck_cgs*units.Energy_MeV/units.Energy_cgs;  // Planck const in MeV.s
	double c = units.clight_cgs;                     // speed of light in cm/s
	double beta, pair_cons, gamma;
	double block_e, block_a, block_x;
	double eps5_p, eps5_m, eps4_p, eps4_m;
	double R_pair, gamma_cons, R_gamma;
	double eta_nua = -eta_nue;
	double eta_nux = 0.0;
	double F3, F4, F5, fac1, fac2;
	double eps_ratio;
	
	(*Q_loc)[0] = (*Q_loc)[1] = (*Q_loc)[2] = 0.0;
	(*R_loc)[0] = (*R_loc)[1] = (*R_loc)[2] = 0.0;

	// Electrons and anti-electrons capture in free nucleons

	beta = M_PI/(pow(h,3.)*c*c)*(1.+3.*alp*alp)/(me*me)*sigma_0;

	F4 = fermi_integral(4,eta_e);
	F5 = fermi_integral(5,eta_e);
	if (F4==0) fac1 = 5.*(1.+0.0287*exp(0.9257*eta_e))/(1.+0.0147*exp(0.9431*eta_e));
        else       fac1 = F5/F4;

	block_e = 1./(1.+exp(-fac1+eta_nue));

	// Fermi integrals of large negative degeneracies might be a numerical pain in the ass. Compute carefully the limits!

	F4 = fermi_integral(4,-eta_e);
	F5 = fermi_integral(5,-eta_e);
	if(F4 == 0.) fac2 = 5.*(1.+0.0287*exp(-0.9257*eta_e))/(1.+0.0147*exp(-0.9431*eta_e));
	else         fac2 = F5/F4;

	block_a = 1./(1.+exp(-fac2+eta_nua));

	(*Q_loc)[0] = beta*eta_np*pow(T,6.0)*fermi_integral(5,eta_e)*block_e;
	(*Q_loc)[1] = beta*eta_pn*pow(T,6.0)*fermi_integral(5,-eta_e)*block_a;
	
	(*R_loc)[0] = beta*eta_np*pow(T,5.0)*fermi_integral(4,eta_e)*block_e;
	(*R_loc)[1] = beta*eta_pn*pow(T,5.0)*fermi_integral(4,-eta_e)*block_a;

	if(TEST) printf("%le %le %le %le\n", (*Q_loc)[0], (*R_loc)[0], (*Q_loc)[1], (*R_loc)[1]);

	// Electron-antielectron pair annihilation

	//Begin by computing the energy moments of Fermi-Dirac distribution

	eps4_m = 8*M_PI/pow((h*c),3.)*fermi_integral(3,eta_e)*pow(T,4.);
	eps4_p = 8*M_PI/pow((h*c),3.)*fermi_integral(3,-eta_e)*pow(T,4.);
	eps5_m = 8*M_PI/pow((h*c),3.)*fermi_integral(4,eta_e)*pow(T,5.);
	eps5_p = 8*M_PI/pow((h*c),3.)*fermi_integral(4,-eta_e)*pow(T,5.);

	pair_cons = sigma_0*c/(me*me)*(eps4_m)*(eps4_p);


	//Blocking factors
	
	//Again, beware with the ratios of fermi integrals

	F3 = fermi_integral(3,-eta_e);
	F4 = fermi_integral(4,-eta_e);
	if(F3 == 0.) fac2 = 4.*(1.+0.0599*exp(-0.90692*eta_e))/(1.+0.0287*exp(-0.9275*eta_e));
	else         fac2 = F4/F3;

	F3 = fermi_integral(3,eta_e);
	F4 = fermi_integral(4,eta_e);
	if(F3 == 0.) fac1 = 4.*(1.+0.0599*exp(0.90692*eta_e))/(1.+0.0287*exp(0.9275*eta_e));
	else         fac1 = F4/F3;

	block_e = 1./(1.+exp(-(0.5*fac1 + 0.5*fac2 -eta_nue)));

	block_a = 1./(1.+exp(-(0.5*fac1 + 0.5*fac2 -eta_nua)));

	block_x = 1./(1.+exp(-(0.5*fac1 + 0.5*fac2-eta_nux)));

	
	//R and Q of electrons neutrinos and antineutrinos from pair annihilation

	R_pair = pair_cons/36.*block_e*block_a*(pow((Cv-Ca),2.)+pow((Cv+Ca),2.));

	(*R_loc)[0] = (*R_loc)[0] + R_pair;
	(*R_loc)[1] = (*R_loc)[1] + R_pair;

	//if (eps4_p*eps4_m != 0) eps_ratio = (eps5_m*eps4_p+eps4_m*eps5_p)/(eps4_p*eps4_m);
	//else                    eps_ratio = ((1.+0.0287*exp(0.9257*eta_e))/(1.+0.0147*exp(0.9431*eta_e)) + (1.+0.0287*exp(-0.9257*eta_e))/(1.+0.0147*exp(-0.9431*eta_e)));

	eps_ratio = (fac1+fac2)*T;

	if (CheckForNANandINF(1, eps_ratio)) {

		printf("eps4_p = %le, eps4_m = %le, eps5_p = %le, eps5_m = %le \n", eps4_p, eps4_m, eps5_p, eps5_m);
		printf("eps_ratio = %le, \n", eps_ratio);

	};

	(*Q_loc)[0] = (*Q_loc)[0] + R_pair*0.5*eps_ratio; // Galeazzi et al. seems strange, for Eq. (A18) 
	(*Q_loc)[1] = (*Q_loc)[1] + R_pair*0.5*eps_ratio; // expresses the total luminosity per baryon produced
												// in pair annihilation, not dividing it between
												// electron neutrino and antineutrino, like Ruffert
												// et al. 1995 and GR1D.

	//R and Q of heavy leptons neutrinos and antineutrinos

	R_pair = pair_cons/9.*pow(block_x,2.)*(pow((Cv-Ca),2.)+pow((Cv+Ca-2.),2.)); // Another error in Galeazzi? In it's Eq. (A16) the blocking
										   // factor only appears once, but it should be squared,
										   // like in Ruffert et al. 1995 and GR1D.
	
	(*R_loc)[2] = R_pair;
	(*Q_loc)[2] = R_pair*eps_ratio;

	// Transversal Plasmon Decay
	gamma = fp*sqrt(1./3.*(M_PI*M_PI+3.*eta_e*eta_e));
	
	gamma_cons = pow(M_PI,3.)*sigma_0*c/(3.*me*me*fsc*pow((h*c),6.))*pow(T,8.)*pow(gamma,6.)*exp(-gamma)*(1.+gamma);

	block_e = 1./(1.+exp(-(1.+0.5*gamma*gamma/(1.+gamma)-eta_nue)));
	block_a = 1./(1.+exp(-(1.+0.5*gamma*gamma/(1.+gamma)-eta_nua)));
	block_x = 1./(1.+exp(-(1.+0.5*gamma*gamma/(1.+gamma)-eta_nux)));

	R_gamma = gamma_cons*Cv*Cv*block_e*block_a;

	//R and Q of electron neutrinos and antineutrinos

	(*R_loc)[0] = (*R_loc)[0] + R_gamma;
	(*R_loc)[1] = (*R_loc)[1] + R_gamma;

	(*Q_loc)[0] = (*Q_loc)[0] + R_gamma*0.5*T*(2.+gamma*gamma/(1.+gamma));
	(*Q_loc)[1] = (*Q_loc)[1] + R_gamma*0.5*T*(2.+gamma*gamma/(1.+gamma));

	//R and Q of heavy leptons neutrinos

	R_gamma = 4.*gamma_cons*(Cv-1.)*(Cv-1.)*block_x*block_x; // Again in Galeazzi the blocking factor appears with inverse signal AND not squared
								 // like in Rufert et al and GR1D.
	
	(*R_loc)[2] = (*R_loc)[2] + R_gamma;
	(*Q_loc)[2] = (*Q_loc)[2] + R_gamma*0.5*T*(2.+gamma*gamma/(1.+gamma));


}


/********************** Diffusive Emission Rates **********************************************/


void compute_diff_emission (double nb, double T, double *nn, double *np, double *nh, double *A, double *Z, double *eta_nue, double *eta_e, double *gxx,
			    double *gyy, double *gzz, double dx, double dy, double dz, double **Q_diff, double **R_diff, double **leak_tau,
			    int keyx, int keyy, int keyz)
{

	/* Included processes: 1. Neutrinos scattering on free nucleons
			       2. Neutrinos scattering on heavy nuclei
 		               3. Electron-flavor neutrinos absorption on free nucleons */


        tU units;
        set_units(&units);
        double h = units.hplanck_cgs*units.Energy_MeV/units.Energy_cgs;  // Planck const in MeV.s
        double c = units.clight_cgs;                     // speed of light in cm/s
	double zeta[7][3]; // Energy independent opacity: [grid point index][neutrino species]
	double chi[3][3];  // Energy independent optical depth for the i-th direction: [direction index][neutrino species]
	double chi_min[3]; // Minimum energy independent optical depth: [neutrino species]
	double scatt_cons, abs_cons, R_cons;
	double block_nue, block_nua;
	double eta_eq[3];
	int i;
        double F5, F4, F2, fac;


	scatt_cons = 0.25*sigma_0/(me*me);
	abs_cons = (1.+3.*alp*alp)*scatt_cons;

	//Compute zeta for all input grid points
	for(i=0;i<7;i++)
	{
		// Blocking factors
		F4 = fermi_integral(4,eta_nue[i]);
		F5 = fermi_integral(5,eta_nue[i]);
		if (F4 == 0.) block_nue = 1;
		else          block_nue = 1./(1.+exp(-F5/F4+eta_e[i]));

		// This hack is needed to avoid nans in block_nua since for large -eta_nue, F4 -> 0, F5 ->0
		F4 = fermi_integral(4,-eta_nue[i]);
		F5 = fermi_integral(5,-eta_nue[i]);
		if(F4 == 0.) block_nua = 1.; // this is the correct limit
		else         block_nua = 1./(1.+exp(-F5/F4-eta_e[i])); //here eta_nua = -eta_nue 

		zeta[i][0] = (np[i]+nn[i])*scatt_cons + nn[i]*abs_cons*block_nue;
		zeta[i][1] = (np[i]+nn[i])*scatt_cons + np[i]*abs_cons*block_nua;
		zeta[i][2] = (np[i]+nn[i])*scatt_cons;
		if(A[i]!=0.){
			zeta[i][0] = zeta[i][0] + nh[i]*scatt_cons/4.*A[i]*A[i]*(1.-Z[i]/A[i])*(1.-Z[i]/A[i]);
			zeta[i][1] = zeta[i][1] + nh[i]*scatt_cons/4.*A[i]*A[i]*(1.-Z[i]/A[i])*(1.-Z[i]/A[i]);
			zeta[i][2] = zeta[i][2] + nh[i]*scatt_cons/4.*A[i]*A[i]*(1.-Z[i]/A[i])*(1.-Z[i]/A[i]);

		}

	}
	//Integrate zeta along x,y,z directions for each neutrino species

	//Here we have to adopt gxx, gyy, gzz as 3-vectors, with - 0 + -> 1 0 2
	
	dx = dx*units.Length_cgs; dy = dy*units.Length_cgs; dz = dz*units.Length_cgs;

	for(i=0;i<3;i++){
	  if(keyx == 2) chi[0][i] = dx/3.*(zeta[1][i]*sqrt(fabs(gxx[1])) + 4.*zeta[0][i]*sqrt(fabs(gxx[0])) + zeta[2][i]*sqrt(fabs(gxx[2])));
	  else chi[0][i] = dx/3.*(zeta[0][i]*sqrt(fabs(gxx[0])) + 4.*zeta[1][i]*sqrt(fabs(gxx[1])) + zeta[2][i]*sqrt(fabs(gxx[2])));
	  if(keyy == 2)	chi[1][i] = dy/3.*(zeta[3][i]*sqrt(fabs(gyy[1])) + 4.*zeta[0][i]*sqrt(fabs(gyy[0])) + zeta[4][i]*sqrt(fabs(gyy[2])));
	  else chi[1][i] = dy/3.*(zeta[0][i]*sqrt(fabs(gyy[0])) + 4.*zeta[3][i]*sqrt(fabs(gyy[1])) + zeta[4][i]*sqrt(fabs(gyy[2])));
	  if(keyz == 2)	chi[2][i] = dz/3.*(zeta[5][i]*sqrt(fabs(gzz[1])) + 4.*zeta[0][i]*sqrt(fabs(gzz[0])) + zeta[6][i]*sqrt(fabs(gzz[2])));
	  else chi[2][i] = dz/3.*(zeta[0][i]*sqrt(fabs(gzz[0])) + 4.*zeta[5][i]*sqrt(fabs(gzz[1])) + zeta[6][i]*sqrt(fabs(gzz[2])));
		/*if(CheckForNANandINF(3,chi[0][i],chi[1][i],chi[2][i])){
			printf("chi[0][%d] = %le, chi[1][%d] = %le, chi[2][%d] = %le\n", i, chi[0][i], i, chi[1][i], i, chi[2][i]);
			printf("zeta[0][%d] = %le, zeta[1][%d] = %le, zeta[2][%d] = %le\n", i, zeta[0][i], i, zeta[1][i], i, zeta[2][i]);
                        printf("zeta[3][%d] = %le, zeta[4][%d] = %le\n", i, zeta[3][i], i, zeta[4][i]);
                        printf("zeta[5][%d] = %le, zeta[6][%d] = %le\n", i, zeta[5][i], i, zeta[6][i]);
			printf("gxx[0] = %le, gxx[1] = %le, gxx[2] = %le\n", gxx[0], gxx[1], gxx[2]);
                        printf("gyy[0] = %le, gyy[1] = %le, gyy[2] = %le\n", gyy[0], gyy[1], gyy[2]);
                        printf("gzz[0] = %le, gzz[1] = %le, gzz[2] = %le\n", gzz[0], gzz[1], gzz[2]);
		}*/
	}

	//Take the minimum chi along x,y,z directions for each neutrino species

	chi_min[0] = (chi[0][0] < chi[1][0]) ? chi[0][0]:chi[1][0];
	chi_min[0] = (chi_min[0] < chi[2][0]) ? chi_min[0]:chi[2][0];

	chi_min[1] = (chi[0][1] < chi[1][1]) ? chi[0][1]:chi[1][1];
	chi_min[1] = (chi_min[1] < chi[2][1]) ? chi_min[1]:chi[2][1];

	chi_min[2] = (chi[0][2] < chi[1][2]) ? chi[0][2]:chi[1][2];
	chi_min[2] = (chi_min[2] < chi[2][2]) ? chi_min[2]:chi[2][2];

	//Compute emission rates and optical depth (leak_tau)

	R_cons = 1./nb*4.*M_PI*c/(3.*pow((h*c),3.));

	eta_eq[0] = eta_nue[0];
	eta_eq[1] = -eta_nue[0];
	eta_eq[2] = 0.0;

	for(i=0;i<3;i++){
		(*R_diff)[i] = R_cons*zeta[0][i]/(chi_min[i]*chi_min[i])*T*fermi_integral(0,eta_eq[i]);
		(*Q_diff)[i] = R_cons*zeta[0][i]/(chi_min[i]*chi_min[i])*T*T*fermi_integral(1,eta_eq[i]);
		//Again, be careful of the fermi integrals ratio for negative degeneracies
		F2 = fermi_integral(2,eta_eq[i]);
		F4 = fermi_integral(4,eta_eq[i]);
		if(F2 == 0.) fac = 12.*(1.+0.1092*exp(0.8908*eta_eq[i]))/(1.+0.0287*exp(0.9257*eta_eq[i]));
		else fac = F4/F2;
		(*leak_tau)[i] = chi_min[i]*T*T*fac;
		//For heavy lepton neutrinos, multiply by degeneracy factor g = 4
		if(i==2){
			(*R_diff)[i] = 4.*(*R_diff)[i];
			(*Q_diff)[i] = 4.*(*Q_diff)[i];
		}

		//if(CheckForNANandINF(2, (*R_diff)[i], (*Q_diff)[i])){
		//	printf("chi_min[%d] = %le, zeta[0][%d] = %le, T = %le, R_cons = %le, eta = %le, F0(eta) = %le, F1(eta) = %le\n",
		//		i, chi_min[i], i, zeta[0][i], T, R_cons, eta_eq[i], fermi_integral(0,eta_eq[i]), fermi_integral(1,eta_eq[i]));

		//}
	}

}



/********* Fermi Integrals *********/


//Analytical fits taken from Takahashi et al. (1987)

double fermi_integral(int k, double eta)
{

 double F = 0.0;

 if(eta>1e-03){

    switch(k)
	{
		case 0:
		F = eta + log(1.0+exp(-eta));
		break;

		case 1:
		F = (pow(eta,2.0)/2.0+1.6449)/(1.0+exp(-1.6855*eta));
		break;

		case 2:
		F = (pow(eta,3.0)/3.0+3.2899*eta)/(1.0-exp(-1.8246*eta));
		break;

		case 3:
		F = (pow(eta,4.0)/4.0+4.9348*pow(eta,2.0)+11.3644)/(1.0+exp(-1.9039*eta));
		break;

		case 4:
		F = (pow(eta,5.0)/5.0+6.5797*pow(eta,3.0)+45.4576*eta)/(1.0-exp(-1.9484*eta));
		break;

		case 5:
		F = (pow(eta,6.0)/6.0+8.2247*pow(eta,4.0)+113.6439*pow(eta,2.0)+236.5323)/(1.0+exp(-1.9727*eta));
		break;

	}

 } else {

    switch(k)
	{
		case 0:
		F = log(1.0+exp(eta));
        break;

		case 1:
		F = exp(eta)/(1.0+0.2159*exp(0.8857*eta));
		break;

		case 2:
		F = 2.0*exp(eta)/(1.0+0.1092*exp(0.8908*eta));
		break;

		case 3:
		F = 6.0*exp(eta)/(1.0+0.0559*exp(0.9069*eta));
		break;

		case 4:
		F = 24.0*exp(eta)/(1.0+0.0287*exp(0.9257*eta));
		break;

		case 5:
		F = 120.0*exp(eta)/(1.0+0.0147*exp(0.9431*eta));
		break;

	}
 }

  return F;
}

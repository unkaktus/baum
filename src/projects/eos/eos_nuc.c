#include "bam.h"
#include "eos.h"

#define PR 0
#define TEST 0

struct {
	int Nr, NT, NY;
	double *Logr, *T, *Y;
	double *eps, *p, *cs2; //hydrodynamics
	double  *mun, *mup, *mue, *Abar, *Zbar, *nh, *np, *nn; //microphysics
	double mb;
	int USET0;
} tabT0;

struct {
	int Nr, NT, NY;
	double *Logr, *LogT, *Y;
	double *eps, *p, *cs2; //hydrodynamics
	double *mun, *mup, *mue, *Abar, *Zbar, *nh, *np, *nn; //microphysics
	double mb;
	double dLogr, dLogT, dY;
    double *chi, *kappa; // derivs for HO flux
	double *Qe, *Qa, *Qx;  //NLS
	double *Re, *Ra, *Rx;
	double *ze, *za, *zx;
	double *eta_nue;
	double *ka_0_nue, *ka_1_nue, *ks_nue; // opacities for M1
    double *ka_0_nua, *ka_1_nua, *ks_nua;
    double *ka_0_nux, *ka_1_nux, *ks_nux;
	double *Yeq, *Teq;
} tab3d_nuclear;

struct {
	int N;
	double *Logr, *T, *Y;
	double *eps, *p, *cs2;
	double *mun, *mup, *mue, *Abar, *Zbar, *nh;
	double *d2Ydr2, *d2pdr2, *d2epsdr2, *d2Tdr2, *d2cs2dr2;
        double *d2mundr2, *d2muedr2, *d2mupdr2, *d2Adr2, *d2Zdr2, *d2nhdr2;
	double *aY, *bY, *cY, *dY;
	double *aT, *bT, *cT, *dT;
	double *ap, *bp, *cp, *dp;
	double *aeps, *beps, *ceps, *deps;
        double mb;
} tab1d_hot_nuclear;


/************************************************* Eos Loaders *************************************/



void find_mb_T0 (char *fname1, char *fname2, char *fname3, double *mb)
{

  tU units;
  set_units(&units);

  double m1, m2, m3;
  int i,j,countl1, countl2, countl3;
  char line1[1024], line2[1024], line3[1024];

  i=j=countl1=countl2=countl3 = 0;

  printf("  3D EoS (with tab T0): searching mb\n");

  FILE *f1 = fopen(fname1, "r"); //open tab T0
  FILE *f2 = fopen(fname2, "r"); //open tab 3D
  FILE *f3 = fopen(fname3, "r"); //open tab_beta (1D)

  while(fgets(line1,1024,f1)) countl1 = countl1 + 1;

  while(fgets(line2,1024,f2)) countl2 = countl2 + 1;

  while(fgets(line3,1024,f3)) countl3 = countl3 + 1;

  rewind(f1); rewind(f2); rewind(f3);

  // Allocate temporary matrices
  double **t_temp1 = malloc(countl1*sizeof(double*));
  for(i=0;i<countl1;i++) t_temp1[i] = malloc(13*sizeof(double));
  double **t_temp2 = malloc(countl2*sizeof(double*));
  for(i=0;i<countl2;i++) t_temp2[i] = malloc(13*sizeof(double));
  double **t_temp3 = malloc(countl3*sizeof(double*));
  for(i=0;i<countl3;i++) t_temp3[i] = malloc(13*sizeof(double));

  // Read on matrices
  for(i=0;i<countl1;i++){
	for(j=0;j<13;j++) fscanf(f1, "%le", &t_temp1[i][j]);
  }

  for(i=0;i<countl2;i++){
	for(j=0;j<13;j++) fscanf(f2, "%le", &t_temp2[i][j]);
  }

  for(i=0;i<countl3;i++){
	for(j=0;j<13;j++) fscanf(f3, "%le", &t_temp3[i][j]);
  }

  // Find the minimum of mb in each tab

  m1 = t_temp1[0][8]; m2 = t_temp2[0][8]; m3 = t_temp3[0][8];
  for(i=1;i<countl1;i++){
	if(t_temp1[i][8] <= m1) m1 = t_temp1[i][8];
  }
  for(i=1;i<countl2;i++){
        if(t_temp2[i][8] <= m2) m2 = t_temp2[i][8];
  }
  for(i=1;i<countl3;i++){
        if(t_temp3[i][8] <= m3) m3 = t_temp3[i][8];
  }

  m1 = (m1 <= m2) ? m1 : m2;
  m1 = (m1 <= m3) ? m1 : m3;

  *mb = m1*units.Mass_cgs/units.Energy_MeV;

 fclose(f1); fclose(f2); fclose(f3);
 for(i=0;i<countl1;i++) free(t_temp1[i]);
 free(t_temp1);
 for(i=0;i<countl2;i++) free(t_temp2[i]);
 free(t_temp2);
 for(i=0;i<countl3;i++) free(t_temp3[i]);
 free(t_temp3);

}

void find_mb (char *fname1, char *fname2, double *mb, double *mb1d)
{

  tU units;
  set_units(&units);

  double m1, m2;
  int i,j,countl1, countl2, Ncol;
  char line1[1024], line2[1024];

  i=j=countl1=countl2=0;

  printf("  3D EoS: searching mb\n");
  printf(" \n %s \n", fname1);
  printf(" \n %s \n", fname2);

  FILE *f1;
  FILE *f2;

  f1 = fopen(fname1, "r"); //open tab 3D

  f2 = fopen(fname2, "r"); //open tab_beta 1D

  while(fgets(line1,1024,f1)) countl1 = countl1 + 1;

  while(fgets(line2,1024,f2)) countl2 = countl2 + 1;

  if(PR) printf("l1 = %d, l2 = %d\n", countl1, countl2);

  rewind(f1); rewind(f2);

  // Allocate temporary matrices
  double **t_temp1 = malloc(countl1*sizeof(double*));
  if(Getv("eos_compo","npe")){
	for(i=0;i<countl1;i++) t_temp1[i] = malloc(13*sizeof(double));
	Ncol = 13;
  }
  else if(Getv("eos_compo","npeH")){
	for(i=0;i<countl1;i++) t_temp1[i] = malloc(15*sizeof(double));
	Ncol = 15;
  }
  double **t_temp2 = malloc(countl2*sizeof(double*));
  for(i=0;i<countl2;i++) t_temp2[i] = malloc(13*sizeof(double));

  // Read on matrices
  for(i=0;i<countl1;i++){
	for(j=0;j<Ncol;j++){
		fscanf(f1, "%le", &t_temp1[i][j]);
	}
  }

  for(i=0;i<countl2;i++){
	for(j=0;j<13;j++) fscanf(f2, "%le", &t_temp2[i][j]);
  }

  // Find the minimum of mb in each tab

  m1 = t_temp1[0][8]; m2 = t_temp2[0][8];
  for(i=1;i<countl1;i++){
	if(t_temp1[i][8] <= m1) m1 = t_temp1[i][8];
  }
  for(i=1;i<countl2;i++){
        if(t_temp2[i][8] <= m2) m2 = t_temp2[i][8];
  }

  *mb1d = m2*units.Mass_cgs/units.Energy_MeV;

  m1 = (m1 <= m2) ? m1 : m2;

  *mb = m1*units.Mass_cgs/units.Energy_MeV;

  fclose(f1); fclose(f2);
  for(i=0;i<countl1;i++) free(t_temp1[i]);
  free(t_temp1);
  for(i=0;i<countl2;i++) free(t_temp2[i]);
  free(t_temp2);

  
  printf("\n m1 = %.16le, m2 = %.16le\n", m1, m2);
  printf("\n m1 = %.16le, m2 = %.16le\n", m1*units.Mass_cgs/units.Energy_MeV, m2*units.Mass_cgs/units.Energy_MeV);

}

void eos_load_tabT(char *fname, int *n1, int *n2, int *n3,
                        double **x1, double **x2, double **x3,
                        double **v1, double **v2, double **v3, double **v4,
			double **v5, double **v6, double **v7, double **v8,
			double **v9, double **v10, double **v11, double *mb)

{

 tU units;
 set_units(&units);

 int i,j,countl, P, Ncol;
 int N1, N2, N3;
 char line[1024];

 printf("Reading EoS:\n%s\n", fname);
 i=j=countl=0;
 N1=N2=N3=1;


 (*x1) = (double *) malloc(1*sizeof(double));
 (*x2) = (double *) malloc(1*sizeof(double));
 (*x3) = (double *) malloc(1*sizeof(double));
 (*v1) = (double *) malloc(1*sizeof(double));
 (*v2) = (double *) malloc(1*sizeof(double));
 (*v3) = (double *) malloc(1*sizeof(double));
 (*v4) = (double *) malloc(1*sizeof(double));
 (*v5) = (double *) malloc(1*sizeof(double));
 (*v6) = (double *) malloc(1*sizeof(double));
 (*v7) = (double *) malloc(1*sizeof(double));
 (*v8) = (double *) malloc(1*sizeof(double));
 (*v8) = (double *) malloc(1*sizeof(double));
 (*v10) = (double *) malloc(1*sizeof(double));
 (*v11) = (double *) malloc(1*sizeof(double));

 for (P=0; P<=bampi_size();P++){

    printf("Reading on proc%d/%d\n", P, bampi_size());

    if(bampi_rank()==P){
 
	FILE *f=fopen(fname,"r");


	//Read file into the matrix t_temp

	while(fgets(line,1024,f)) countl=countl+1;

 	rewind(f);
	double **t_temp=malloc(countl*sizeof(double*));
	if(Getv("eos_compo","npe")){
	 for(i=0;i<countl;i++) t_temp[i] = malloc(13*sizeof(double));
	 Ncol = 13;
	}
  	else if(Getv("eos_compo","npeH")){
	 for(i=0;i<countl;i++) t_temp[i] = malloc(15*sizeof(double));
	 Ncol = 15;
	}
 
	for(i=0;i<countl;i++){
		for(j=0;j<Ncol;j++){
			fscanf(f, "%le", &t_temp[i][j]);
		}
 	}

 	//Count N3
 	i=1;
 	while(t_temp[i][2]!=t_temp[0][2]){
		N3=N3+1;
		i++;
 	}

 	//Count N2
 	i=1;
 	while(t_temp[i][0]==t_temp[0][0]){
		N2=N2+1;
		i++;
 	}

 	N2=N2/N3;

 	//Bonus: N1
 	N1=countl/(N2*N3);

 	*n1 = N1;
 	*n2 = N2;
 	*n3 = N3;

 	//Filling the vectors

 	*x1 = (double *) realloc (*x1, N1*sizeof(double));
 	*x2 = (double *) realloc (*x2, N2*sizeof(double));
 	*x3 = (double *) realloc (*x3, N3*sizeof(double));
 	*v1 = (double *) realloc (*v1, N1*N2*N3*sizeof(double));
 	*v2 = (double *) realloc (*v2, N1*N2*N3*sizeof(double));
 	*v3 = (double *) realloc (*v3, N1*N2*N3*sizeof(double));
 	*v4 = (double *) realloc (*v4, N1*N2*N3*sizeof(double));
 	*v5 = (double *) realloc (*v5, N1*N2*N3*sizeof(double));
 	*v6 = (double *) realloc (*v6, N1*N2*N3*sizeof(double));
 	*v7 = (double *) realloc (*v7, N1*N2*N3*sizeof(double));
 	*v8 = (double *) realloc (*v8, N1*N2*N3*sizeof(double));
 	*v9 = (double *) realloc (*v9, N1*N2*N3*sizeof(double));
    *v10 = (double *) realloc (*v10, N1*N2*N3*sizeof(double));
    *v11 = (double *) realloc (*v11, N1*N2*N3*sizeof(double));

	double m_n = 939.565413; //neutron mass in MeV, c = 1

 	//Electron fraction (dimensionless)
	for(i=0;i<N3;i++){
		(*x3)[i]=t_temp[i][2];
 	}
 
 	j=0; //log10(rest-mass density) in log10(g*cm^-3)
 	for(i=0;i<N2;i++){
		(*x2)[i]=log10(t_temp[j][1]*(*mb)*1e+39/units.Mdens_cgs);
		j = j+N3;
		//printf("Logr_tab[%d] = %le\n", i, (*x2)[i]);
 	}

 	j=0; //log10(temperature) in log10(MeV)
 	for(i=0;i<N1;i++){
		(*x1)[i]=log10(t_temp[j][0]);
		j = j+N2*N3;
		//printf("LogT_tab[%d] = %le\n", i, (*x1)[i]);
 	}


 	for(i=0;i<N1*N2*N3;i++){
		(*v1)[i] = t_temp[i][3]/units.Energy_MeV*units.Volume_fm3; // Pressure
		(*v2)[i] = (t_temp[i][8]*units.Mass_cgs/units.Energy_MeV)/(*mb) - 1.0; //specific internal energy epsilon (dimensionless)
		(*v3)[i] = t_temp[i][7]; //(speed of sound)^2 in c = 1
		(*v4)[i] = t_temp[i][4]+m_n; //neutrons chemical potential in MeV
		(*v5)[i] = t_temp[i][5]+(*v4)[i]; //protons chemical potential in MeV
		(*v6)[i] = t_temp[i][6]-t_temp[i][5]; //electrons chemical potential in MeV

		if(Getv("eos_compo", "npe")){
			(*v7)[i] = t_temp[i][10]; //avg heavy nuclei mass number Abar (dimensionless)
			(*v8)[i] = t_temp[i][11]; //avg heavy nuclei charge number Zbar (dimensionless)
			(*v9)[i] = t_temp[i][9]*t_temp[i][1]*1e+39; //avg heavy nuclei number density (cm^-3)
			(*v10)[i] = (1.-t_temp[i][2])*t_temp[i][1]*1e+39; // neutrons number density (cm^-3)
			(*v11)[i] = t_temp[i][2]*t_temp[i][1]*1e+39;	  // protons number density (cm^-3)
		}
		else if(Getv("eos_compo", "npeH")){
			(*v7)[i] = t_temp[i][12]; //avg heavy nuclei mass number Abar
			(*v8)[i] = t_temp[i][13]; //avg heavy nuclei charge number Zbar
			(*v9)[i] = DMAX(0.,t_temp[i][11]*t_temp[i][1]*1e+39); // avg heavy nuclei number density (cm^-3)
			(*v10)[i] = DMAX(0.,t_temp[i][9]*t_temp[i][1]*1e+39); // neutrons number density (cm^-3)
			(*v11)[i] = DMAX(0.,t_temp[i][10]*t_temp[i][1]*1e+39); // protons number density (cm^-3)
		}
 	}

 	fclose(f);
	for(i=0;i<N1*N2*N3;i++) free(t_temp[i]);
 	free(t_temp);

    }

 }
	printf("Nr = %d, NT = %d, NY = %d\n", *n2, *n1, *n3);
	printf("mb = %le\n", *mb);
	printf("3D EoS loaded successfully from file \n%s \n", fname);

}

//Load tab for T = 0
void eos_load_tabT0(char *fname, int *n1, int *n2, int *n3,
                        double **x1, double **x2, double **x3,
                        double **v1, double **v2, double **v3, double **v4,
			double **v5, double **v6, double **v7, double **v8,
			double **v9, double *mb)

{

 tU units;
 set_units(&units);

 int i,j,countl, P;
 int N1, N2, N3;
 char line[1024];

 printf("Reading EoS:\n%s\n", fname);
 i=j=countl=0;
 N1=N2=N3=1;


 (*x1) = (double *) malloc(1*sizeof(double));
 (*x2) = (double *) malloc(1*sizeof(double));
 (*x3) = (double *) malloc(1*sizeof(double));
 (*v1) = (double *) malloc(1*sizeof(double));
 (*v2) = (double *) malloc(1*sizeof(double));
 (*v3) = (double *) malloc(1*sizeof(double));
 (*v4) = (double *) malloc(1*sizeof(double));
 (*v5) = (double *) malloc(1*sizeof(double));
 (*v6) = (double *) malloc(1*sizeof(double));
 (*v7) = (double *) malloc(1*sizeof(double));
 (*v8) = (double *) malloc(1*sizeof(double));
 (*v9) = (double *) malloc(1*sizeof(double));

 for (P=0; P<=bampi_size();P++){

    printf("Reading on proc%d/%d\n", P, bampi_size());

    if(bampi_rank()==P){
 
	FILE *f=fopen(fname,"r");


	//Read file into the matrix t_temp

	while(fgets(line,1024,f)) countl=countl+1;

 	rewind(f);
	double **t_temp=malloc(countl*sizeof(double*));
	for(i=0;i<countl;i++) t_temp[i] = malloc(13*sizeof(double));
 
	for(i=0;i<countl;i++){
		for(j=0;j<13;j++){
			fscanf(f, "%le", &t_temp[i][j]);
		}
 	}

 	//Count N3
 	i=1;
 	while(t_temp[i][2]!=t_temp[0][2]){
		N3=N3+1;
		i++;
 	}

 	//Count N2
 	i=1;
 	while(t_temp[i][0]==t_temp[0][0]){
		N2=N2+1;
		i++;
 	}

 	N2=N2/N3;

 	//Bonus: N1
 	N1=countl/(N2*N3);

 	*n1 = N1;
 	*n2 = N2;
 	*n3 = N3;

 	//Filling the vectors

 	*x1 = (double *) realloc (*x1, N1*sizeof(double));
 	*x2 = (double *) realloc (*x2, N2*sizeof(double));
 	*x3 = (double *) realloc (*x3, N3*sizeof(double));
 	*v1 = (double *) realloc (*v1, N1*N2*N3*sizeof(double));
 	*v2 = (double *) realloc (*v2, N1*N2*N3*sizeof(double));
 	*v3 = (double *) realloc (*v3, N1*N2*N3*sizeof(double));
 	*v4 = (double *) realloc (*v4, N1*N2*N3*sizeof(double));
 	*v5 = (double *) realloc (*v5, N1*N2*N3*sizeof(double));
 	*v6 = (double *) realloc (*v6, N1*N2*N3*sizeof(double));
 	*v7 = (double *) realloc (*v7, N1*N2*N3*sizeof(double));
 	*v8 = (double *) realloc (*v8, N1*N2*N3*sizeof(double));
 	*v9 = (double *) realloc (*v9, N1*N2*N3*sizeof(double));


	double m_n = 939.565413; //neutron mass in MeV, c = 1

 	//Electron fraction (dimensionless)
	for(i=0;i<N3;i++){
		(*x3)[i]=t_temp[i][2];
 	}
 
 	j=0; //log10(rest-mass density)
 	for(i=0;i<N2;i++){
		(*x2)[i]=log10(t_temp[j][1]*(*mb)*1e+39/units.Mdens_cgs);
		j = j+N3;
 	}

 	j=0; //log10(temperature) in log10(MeV)
 	for(i=0;i<N1;i++){
		(*x1)[i]=t_temp[j][0];
		j = j+N2*N3;
 	}


 	for(i=0;i<N1*N2*N3;i++){
		(*v1)[i] = t_temp[i][3]/units.Energy_MeV*units.Volume_fm3; // Pressure
		(*v2)[i] = (t_temp[i][8]*units.Mass_cgs/units.Energy_MeV)/(*mb) - 1.0; //specific internal energy epsilon (dimensionless)
		(*v3)[i] = t_temp[i][7]; //(speed of sound)^2 in c = 1
		(*v4)[i] = t_temp[i][4]+m_n; //neutrons chemical potential in MeV
		(*v5)[i] = t_temp[i][5]+(*v4)[i]; //protons chemical potential in MeV
		(*v6)[i] = t_temp[i][6]-t_temp[i][5]; //electrons chemical potential in MeV
		(*v7)[i] = t_temp[i][10]; //avg heavy nuclei mass number Abar (dimensionless)
		(*v8)[i] = t_temp[i][11]; //avg heavy nuclei charge number Zbar (dimensionless)
		(*v9)[i] = t_temp[i][9]*t_temp[i][1]*1e+39; //avg heavy nuclei number density (cm^-3)
 	}

 	fclose(f);
	for(i=0;i<N1*N2*N3;i++) free(t_temp[i]);
 	free(t_temp);

    }

 }
	printf("Nr = %d, NT = %d, NY = %d\n", *n2, *n1, *n3);
	printf("mb = %.16le \n", *mb);
	printf("3D EoS tab T0 loaded successfully from file \n%s \n", fname);

}

// Compose's derivatives are way better. For this one, we require an additional file containing variables 13 & 17,
// which are 'specific heat capacity at constant volume c_V' & 'tension coefficient at constant volume beta_V'

void eos_load_derivs(char *fname, double **v1, double **v2)
{


 tU units;
 set_units(&units);

 int i,j,P;
 int N = tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY;
 double rho;

 //kappa
 (*v1) = (double *) malloc(1*sizeof(double));
 //chi
 (*v2) = (double *) malloc(1*sizeof(double));

 for (P=0; P<=bampi_size();P++){

    printf("Reading on proc%d/%d\n", P, bampi_size());

    if(bampi_rank()==P){

        printf("Reading EoS derivs:\n%s\n", fname);

        FILE *f=fopen(fname,"r");

	*v1 = (double *) realloc (*v1, N*sizeof(double));
        *v2 = (double *) realloc (*v2, N*sizeof(double));

        double **t_temp=malloc(N*sizeof(double*));
        for(i=0;i<N;i++) t_temp[i] = malloc(5*sizeof(double));

        for(i=0;i<N;i++){
                for(j=0;j<5;j++){
                        fscanf(f, "%le", &t_temp[i][j]);
                }
		//kappa 
		(*v1)[i] = (t_temp[i][4]*tab3d_nuclear.mb*1e+39/units.Mdens_cgs)/(t_temp[i][3]); //dimensionless
		rho = t_temp[i][1]*tab3d_nuclear.mb*1e+39/units.Mdens_cgs;
		//chi
		(*v2)[i] = (1. + tab3d_nuclear.eps[i] + tab3d_nuclear.p[i]/rho)*tab3d_nuclear.cs2[i] - tab3d_nuclear.p[i]/(rho*rho)*(*v1)[i]; //dimensionless
        }

        fclose(f);
        for(i=0;i<N;i++) free(t_temp[i]);
        free(t_temp);

    }
 }

 printf("3D EoS tab derivs loaded successfully from file \n%s \n", fname);

}

void eos_load_M1_eq(char *fname, double **v1, double **v2)
{

 tU units;
 set_units(&units);

 int i,j,P;
 int N = tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY;
 double rho;

 (*v1) = (double *) malloc(1*sizeof(double));
 (*v2) = (double *) malloc(1*sizeof(double));

 for (P=0; P<=bampi_size();P++) {

    printf("Reading on proc%d/%d\n", P, bampi_size());

    if(bampi_rank()==P){

        printf("Reading EoS M1 eqilibrium:\n%s\n", fname);

        FILE *f=fopen(fname,"r");

		*v1 = (double *) realloc (*v1, N*sizeof(double));
        *v2 = (double *) realloc (*v2, N*sizeof(double));

        double **t_temp=malloc(N*sizeof(double*));
        for(i=0;i<N;i++) t_temp[i] = malloc(5*sizeof(double));

        for(i=0;i<N;i++){
                for(j=0;j<5;j++){
                        fscanf(f, "%le", &t_temp[i][j]);
                }
		 
			(*v1)[i] = t_temp[i][3];
			(*v2)[i] = t_temp[i][4];
        }

        fclose(f);
        for(i=0;i<N;i++) free(t_temp[i]);
        free(t_temp);

    }
 }

 printf("3D EoS tab M1 equilibrium loaded successfully from file \n%s \n", fname);

}

/************************** Expand T=0 table and set EoS structs *********************************************/

void eos_load_complete()
{

 int i, j, k, ijk_T0, ijk, key=0;

 double en_shift, p, rho, dummy, epsl;
 double *ptr;

 //Load T=0 table

 if(!Getv("eos_tab_file_T0","none")){

		find_mb_T0(Gets("eos_tab_file_T0"),Gets("eos_tab_file"),Gets("eos_tab_file_beta"), &(tab3d_nuclear.mb));

		eos_load_tabT0(Gets("eos_tab_file_T0"), &(tabT0.NT), &(tabT0.Nr), &(tabT0.NY), &(tabT0.T), 
		&(tabT0.Logr), &(tabT0.Y), &(tabT0.p), &(tabT0.eps), &(tabT0.cs2),
		&(tabT0.mun), &(tabT0.mup), &(tabT0.mue), &(tabT0.Abar), &(tabT0.Zbar), &(tabT0.nh), &(tab3d_nuclear.mb));
		tabT0.USET0 = 1;
 } else { 

		tabT0.USET0 = 0;
 		find_mb(Gets("eos_tab_file"),Gets("eos_tab_file_beta"), &(tab3d_nuclear.mb), &(tab1d_hot_nuclear.mb));
 }

 //Load Beta-equilibrated table

 eos_read_tab1d_hot(Gets("eos_tab_file_beta"), &(tab1d_hot_nuclear.N), &(tab1d_hot_nuclear.T), &(tab1d_hot_nuclear.Logr), &(tab1d_hot_nuclear.Y),
                    &(tab1d_hot_nuclear.p), &(tab1d_hot_nuclear.eps), &(tab1d_hot_nuclear.cs2), &(tab1d_hot_nuclear.mun), &(tab1d_hot_nuclear.mup),
                    &(tab1d_hot_nuclear.mue), &(tab1d_hot_nuclear.Abar), &(tab1d_hot_nuclear.Zbar), &(tab1d_hot_nuclear.nh), &(tab3d_nuclear.mb), key);

 if(Getv("eos_interp", "cspline")){
	spline(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.T, tab1d_hot_nuclear.N, &(tab1d_hot_nuclear.d2Tdr2));
 	spline(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.Y, tab1d_hot_nuclear.N, &(tab1d_hot_nuclear.d2Ydr2));
 	spline(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.p, tab1d_hot_nuclear.N, &(tab1d_hot_nuclear.d2pdr2));
 	spline(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.eps, tab1d_hot_nuclear.N, &(tab1d_hot_nuclear.d2epsdr2));
 	spline(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.cs2, tab1d_hot_nuclear.N, &(tab1d_hot_nuclear.d2cs2dr2));
 }
 else if(Getv("eos_interp", "steffen")){
	steffen_init(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.p, tab1d_hot_nuclear.N, &(ptr), &(tab1d_hot_nuclear.ap), &(tab1d_hot_nuclear.bp), &(tab1d_hot_nuclear.cp), &(tab1d_hot_nuclear.dp));
	steffen_init(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.Y, tab1d_hot_nuclear.N, &(ptr), &(tab1d_hot_nuclear.aY), &(tab1d_hot_nuclear.bY), &(tab1d_hot_nuclear.cY), &(tab1d_hot_nuclear.dY));
        steffen_init(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.eps, tab1d_hot_nuclear.N, &(ptr), &(tab1d_hot_nuclear.aeps), &(tab1d_hot_nuclear.beps), &(tab1d_hot_nuclear.ceps), &(tab1d_hot_nuclear.deps));
        steffen_init(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.T, tab1d_hot_nuclear.N, &(ptr), &(tab1d_hot_nuclear.aT), &(tab1d_hot_nuclear.bT), &(tab1d_hot_nuclear.cT), &(tab1d_hot_nuclear.dT));
 }

 //Load T>0 table

 eos_load_tabT(Gets("eos_tab_file"), &(tab3d_nuclear.NT), &(tab3d_nuclear.Nr), &(tab3d_nuclear.NY), &(tab3d_nuclear.LogT), 
		&(tab3d_nuclear.Logr), &(tab3d_nuclear.Y), &(tab3d_nuclear.p), &(tab3d_nuclear.eps), &(tab3d_nuclear.cs2),
		&(tab3d_nuclear.mun), &(tab3d_nuclear.mup), &(tab3d_nuclear.mue), &(tab3d_nuclear.Abar), &(tab3d_nuclear.Zbar), &(tab3d_nuclear.nh), &(tab3d_nuclear.nn), &(tab3d_nuclear.np), &(tab3d_nuclear.mb));

 if(Getv("hrsc_flux", "HO_LLF")){
	if(Getv("eos_tab_file_derivs","none")) errorexit(" Provide a derivs tab for HO_LLF or change hrsc_flux!");
 	else eos_load_derivs(Gets("eos_tab_file_derivs"), &(tab3d_nuclear.kappa), &(tab3d_nuclear.chi));
 }
 // Add T = T1 (lowest tab3d_nuclear temperature) to T = 0 table to make interp routines functional in the interval 0 <= T < T1

 if(!Getv("eos_tab_file_T0","none"))
 {
	tabT0.NT = tabT0.NT + 1;
 	tabT0.T = (double *) realloc (tabT0.T, tabT0.NT*sizeof(double));
 	tabT0.eps = (double *) realloc (tabT0.eps, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.p = (double *) realloc (tabT0.p, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.cs2 = (double *) realloc (tabT0.cs2, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.mun = (double *) realloc (tabT0.mun, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.mup = (double *) realloc (tabT0.mup, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.mue = (double *) realloc (tabT0.mue, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.Abar = (double *) realloc (tabT0.Abar, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.Zbar = (double *) realloc (tabT0.Zbar, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));
 	tabT0.nh = (double *) realloc (tabT0.nh, tabT0.NT*tabT0.Nr*tabT0.NY*sizeof(double));

 	tabT0.T[1] = pow(10,tab3d_nuclear.LogT[0]);

	//Fill in the new tabT0

 	for(i=0;i<tabT0.Nr;i++){
		for(j=0;j<tabT0.NY;j++){
			ijk_T0 = tabT0.NY*tabT0.Nr + i*tabT0.NY + j;
			ijk = i*tab3d_nuclear.NY + j;
			tabT0.eps[ijk_T0] = tab3d_nuclear.eps[ijk];
			tabT0.p[ijk_T0] = tab3d_nuclear.p[ijk];
			tabT0.cs2[ijk_T0] = tab3d_nuclear.cs2[ijk];
                	tabT0.mun[ijk_T0] = tab3d_nuclear.mun[ijk]; 
                	tabT0.mup[ijk_T0] = tab3d_nuclear.mup[ijk]; 
                	tabT0.mue[ijk_T0] = tab3d_nuclear.mue[ijk];
                	tabT0.Abar[ijk_T0] = tab3d_nuclear.Abar[ijk]; 
                	tabT0.Zbar[ijk_T0] = tab3d_nuclear.Zbar[ijk]; 
                	tabT0.nh[ijk_T0] = tab3d_nuclear.nh[ijk]; 
		}

 	}

 }

 tab3d_nuclear.dLogr = (tab3d_nuclear.Logr[tab3d_nuclear.Nr - 1] - tab3d_nuclear.Logr[0])/ (float) (tab3d_nuclear.Nr - 1);
 tab3d_nuclear.dLogT = (tab3d_nuclear.LogT[tab3d_nuclear.NT - 1] - tab3d_nuclear.LogT[0])/ (float) (tab3d_nuclear.NT - 1);
 tab3d_nuclear.dY    = (tab3d_nuclear.Y[tab3d_nuclear.NY - 1]    - tab3d_nuclear.Y[0])/(float) (tab3d_nuclear.NY - 1);

 if(Getv("grhd_use_nls", "yes") || Getv("grrhd_m1_nls_rates","yes")){
	 nls_table(&(tab3d_nuclear.Qe), &(tab3d_nuclear.Qa), &(tab3d_nuclear.Qx), &(tab3d_nuclear.Re), &(tab3d_nuclear.Ra), &(tab3d_nuclear.Rx), 
		   &(tab3d_nuclear.ze), &(tab3d_nuclear.za), &(tab3d_nuclear.zx), &(tab3d_nuclear.eta_nue),
		   &(tab3d_nuclear.ka_0_nue), &(tab3d_nuclear.ka_1_nue), &(tab3d_nuclear.ks_nue),
		   &(tab3d_nuclear.ka_0_nua), &(tab3d_nuclear.ka_1_nua), &(tab3d_nuclear.ks_nua),
		   &(tab3d_nuclear.ka_0_nux), &(tab3d_nuclear.ka_1_nux), &(tab3d_nuclear.ks_nux));
 }

 if(Getv("grrhd_m1_nls_rates", "yes")){
	// Get rid of unnecessary memory
	free(tab3d_nuclear.ze); free(tab3d_nuclear.za); free(tab3d_nuclear.zx);
 }

 char *fname;
 if (Getv("physics", "grrhd_m1")) {
	fname = Gets("grrhd_m1_EqTab_path");
	eos_load_M1_eq(fname, &(tab3d_nuclear.Yeq), &(tab3d_nuclear.Teq));
 }

}

/********************************************** Interpolations *******************************************************************/

//Trilinear interpolation

int intp3d(double x, double y, double z, double *f, int nx, int ny, int nz, 
	   double dx, double dy, double dz, double *f_t, double *x_t, double *y_t, double *z_t,
      	   double *dfdx, double *dfdy, double *dfdz)
{
 double delx, dely, delz, a1, a2, a3, a4, a5, a6, a7, a8;
 double dxi, dyi, dzi, dxyi, dxzi, dyzi, dxyzi;
 int ix, iy, iz, tempx, tempy, tempz, px, py, pz, ijk;
 double fh[8];


 //Determine spacing parameters for *equidistant* table

 dxi = 1./dx;
 dyi = 1./dy;
 dzi = 1./dz;

 dxyi = dxi*dyi;
 dxzi = dxi*dzi;
 dyzi = dyi*dzi;
 dxyzi = dxi*dyi*dzi;

 //Determine position indices

 ix = 1+(int)((x-(x_t)[0]-1e-10)*dxi);
 iy = 1+(int)((y-(y_t)[0]-1e-10)*dyi);
 iz = 1+(int)((z-(z_t)[0]-1e-10)*dzi);

 //tempx = (ix <= nx-1) ? ix:nx-1;
 //tempy = (iy <= ny-1) ? iy:ny-1;
 //tempz = (iz <= nz-1) ? iz:nz-1;

 //ix = (tempx >= 1) ? tempx:1;
 //iy = (tempy >= 1) ? tempy:1;
 //iz = (tempz >= 1) ? tempz:1;

 //Set up arrays for interpolation
 //if(ix < 1 || iy < 1 || iz < 1) printf("ix = %d, iy = %d, iz = %d\n", ix, iy, iz);

 delx = (x_t)[ix]-x;
 dely = (y_t)[iy]-y; 
 delz = (z_t)[iz]-z; 
 
 px = ny*nz;
 py = nz;
 pz = 1;
 
 ijk = ix*px + iy*py + iz*pz;

 fh[0] = (f_t)[ijk];
 fh[1] = (f_t)[ijk - px];
 fh[2] = (f_t)[ijk - py];
 fh[3] = (f_t)[ijk - pz];
 fh[4] = (f_t)[ijk - px - py];
 fh[5] = (f_t)[ijk - px - pz];
 fh[6] = (f_t)[ijk - py - pz];
 fh[7] = (f_t)[ijk - px - py - pz];

 //Set up coefficients of the interpolation polynomial 

 a1 = fh[0];
 a2 = (fh[1] - fh[0])*dxi;
 a3 = (fh[2] - fh[0])*dyi;
 a4 = (fh[3] - fh[0])*dzi; 
 a5 = (fh[4] - fh[1] - fh[2] + fh[0])*dxyi;
 a6 = (fh[5] - fh[1] - fh[3] + fh[0])*dxzi;
 a7 = (fh[6] - fh[2] - fh[3] + fh[0])*dyzi;
 a8 = (fh[7] - fh[0] + fh[1] + fh[2] + fh[3] - fh[4] - fh[5] - fh[6])*dxyzi;

 *f = a1 + a2*delx + a3*dely + a4*delz + a5*delx*dely + a6*delx*delz + a7*dely*delz + a8*delx*dely*delz;

 //Derivatives (first order accurate)

 *dfdx = -a2;
 *dfdy = -a3;
 *dfdz = -a4;


 /*if(CheckForNANandINF(16, fh[0], fh[1], fh[2], fh[3], fh[4], fh[5], fh[6], fh[7], a1, a2, a3, a4, a5, a6, a7, a8)){
	printf("Error in intp3d: fh[0] = %le, fh[1] = %le, fh[2] = %le, fh[3] = %le, fh[4] = %le, fh[5] = %le, fh[6] = %le, fh[7] = %le\n",
		fh[0], fh[1], fh[2], fh[3], fh[4], fh[5], fh[6], fh[7]);
	printf("a1 = %le, a2 = %le, a3 = %le, a4 = %le, a5 = %le, a6 = %le, a7 = %le, a8 = %le\n", a1, a2, a3, a4, a5, a6, a7, a8);
	printf("Input: x = %le, y = %le, z = %le, ix = %d, iy = %d, iz = %d, dx = %le, dy = %le, dz = %le -----> f = %le\n");
 }*/

 return 1;

}
 
//Cubic splines for beta-equilibrated EoS interpolation

void spline(double *x, double *y, int n, double **y2)
{
        //This routine is called just once to produce the *y2 array containing second derivatives for function the evaluation
	//It's called only once!
        int i, k;
        double p, sig;
        double u[n];

	*y2 = (double *) malloc (n*sizeof(double));

        (*y2)[0] = 0.0;
        (*y2)[n-1] = 0.0;
        u[0] = 0.0;

        for(i=1; i<=n-2; i++){
                sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
                p = sig*(*y2)[i-1]+2.0;
                (*y2)[i] = (sig-1.0)/p;
                u[i] = (y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
                u[i] = (6.0*u[i]/(x[i+1]-x[i])-sig*u[i-1])/p;
        }

        for (k = n-2; k>=0; k--){
                (*y2)[k] = (*y2)[k]*(*y2)[k+1]+u[k];
        }

}

//This routine takes the tabulated function from arrays xa and ya, the array of second derivatives y2a
// and computes the function value y and it's derivative

void splint(double *xa, double *ya, double *y2a, int n, double x, double *y, double *dydx)
{
        int klo, khi, k;
        double h, b, a;

        klo = 0;
        khi = n-1;

	if(PR && (x<xa[0] || x>xa[n-1])) printf("Bad x input in splint!\n x = %e, xmin = %e, xmax = %e\n", x, xa[0], xa[n-1]);
        while(khi-klo > 1){
                k = (khi+klo) >> 1;
                if (xa[k] > x) khi=k;
                else klo = k;
        }

        h = xa[khi]-xa[klo];
        if(h==0.0) printf("Bad xa input"); //This should not happen since we're using dinstinc xa's from the table
        a = (xa[khi] - x)/h;
        b = (x-xa[klo])/h;
        *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
        *dydx = (ya[khi]-ya[klo])/h-(3.0*a*a-1.0)/6.0*h*y2a[klo]+(3.0*b*b-1.0)/6.0*h*y2a[khi];
}




/*************************** Validity range for EoS in the (rho, Y, eps) space *********************************************************/

void eos_validity(double *rho, double *Y, double *rho_min, double *rho_max, 
		  double *Y_min, double *Y_max, double *eps_min, double *eps_max, double *Tmin, 
		  double *Tmax, double *mb, double *rhomin1D, double *mb1D)
{
 double dummy, tmp, eps, rho_tmp, Y_tmp;
 double LogTmin, LogTmax, Logr;
 int i;

 *Y_min = tab3d_nuclear.Y[0];
 *Y_max = tab3d_nuclear.Y[tab3d_nuclear.NY-1];

 *rho_min = pow(10., tab3d_nuclear.Logr[0]);
 *rho_max = pow(10., tab3d_nuclear.Logr[tab3d_nuclear.Nr-1]);

 *rhomin1D = pow(10., tab1d_hot_nuclear.Logr[0]);

 *Tmin = pow(10., tab3d_nuclear.LogT[0]);
 *Tmax = pow(10., tab3d_nuclear.LogT[tab3d_nuclear.NT - 1]);

 *mb = tab3d_nuclear.mb;
 *mb1D = tab1d_hot_nuclear.mb;

}

//Extend validity range for con2prim

int extend_validity(double *rho_hat, double *Y, double *eps_hat, double *T)
{
 double Y_tmp, rho_tmp, eps_tmp, eps, epsmax, epsmin;
 double Tmax, LogTmax, Tmin, LogTmin, Logr, dummy;
 int i;

 Y_tmp = (*Y <= tab3d_nuclear.Y[tab3d_nuclear.NY-1]) ? *Y : tab3d_nuclear.Y[tab3d_nuclear.NY-1];
 *Y    = (Y_tmp >= tab3d_nuclear.Y[0]) ? Y_tmp : tab3d_nuclear.Y[0];

 rho_tmp = (*rho_hat <= pow(10.,tab3d_nuclear.Logr[tab3d_nuclear.Nr - 1])) ? *rho_hat : pow(10.,tab3d_nuclear.Logr[tab3d_nuclear.Nr-1]);
 *rho_hat = (rho_tmp >= pow(10.,tab3d_nuclear.Logr[0])) ? rho_tmp: pow(10.,tab3d_nuclear.Logr[0]);
 Logr = log10(*rho_hat);

 //To future implementation: a condition on the use or not of a T=0 tab

 //intp3d(tabT0.T[0], Logr, *Y, eps_min, tabT0.NT, tabT0.Nr, tabT0.NY, tabT0.eps, tabT0.T, tabT0.Logr, tabT0.Y, &dummy, &dummy, &dummy);
 
 //The minimum is determined for the lowest T in tab3d_nuclear. TODO: include for the case of T = 0 a switch case for eps_min
 intp3d(tab3d_nuclear.LogT[0], Logr, *Y, &epsmin, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, 
					tab3d_nuclear.eps, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
 intp3d(tab3d_nuclear.LogT[tab3d_nuclear.NT-1], Logr, *Y, &epsmax, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY,
					tab3d_nuclear.eps, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
 //if(CheckForNANandINF(1,eps)) printf("LogT = %le, Logr = %le, Y = %le -------> eps = %le\n", tab3d_nuclear.LogT[0], Logr, *Y, eps);

 Tmax = pow(10.,tab3d_nuclear.LogT[tab3d_nuclear.NT-1]);
 Tmin = pow(10.,tab3d_nuclear.LogT[0]);

 //Force *eps into the validity range

 eps_tmp = (*eps_hat <= epsmax) ? *eps_hat : epsmax;
 *eps_hat = (eps_tmp >= epsmin) ? eps_tmp : epsmin;

 if(*eps_hat == epsmax) {
	*T = Tmax;
	return 0;
 }
 else if(*eps_hat == epsmin){
	*T = Tmin;
	return 0;
 }

 find_temp(T, rho_hat, Y, eps_hat, &Tmin, &Tmax, &epsmin, &epsmax);
 return 0;
}

/*********************************************** Temperature finder *******************************************************/

double eps_T(double rho, double Y, double T){
	double eps, Logr = log10(rho), LogT = log10(T), dummy;
	intp3d(LogT, Logr, Y, &eps, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
			tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.eps, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
	return eps;
}


int compute_temp(double *T, double rho, double Y, double epsl, double Tmin, double Tmax, double epsmin, double epsmax){

	double Logr = log10(rho);
	double a, b, c, fa, fb, fc, epsa, epsb, epsc, s, m;
	double epsin, dummy, swap, errormax = 1e-10;
	int count, j, countmax = 100;

	epsin = epsl; // input epsl
	
	a = Tmax; b = Tmin;
        epsa = epsmax;
	epsb = epsmin;

	fa = epsin - epsa;
	fb = epsin - epsb;

	// Do you feel lucky, punk? Let's check if the initial guesses are solutions
	if(fabs(fa) <= errormax*epsin || fabs(fb) <= errormax*epsin){
		if(fabs(fa) <= fabs(fb)){
			*T = a;
			return 0;
		}
		else {
			*T = b;
			return 0;
		}
	}


	//Else, ensure fa is the positive one and get down into illinois
	if(fa < 0.){
		swap = fa; fa = fb; fb = swap;
		swap = a; a = b; b = swap;
	}

	count = 0;
	while(count<=countmax){
		if(fa > -fb){
			s = fb/fa;
			m = 1.0 - s;
			c = b/m - s*a/m;
		}
		else {
			s = fa/fb;
			m = 1.0 - s;
			c = a/m - s*b/m;
		}


		// Check the new estimate
		epsc = eps_T(rho, Y, c);
		fc = epsin - epsc;
		if(fabs(fc)<=errormax*epsin){
			a = c;
			b = c;
			*T = c;
			return 0;
		}


		// Update the interval
		
		(fc > 0.0) ? (a=c, fa=fc) : (b=c, fb=fc);
		if(fc > 0.0){
			fb = fb*0.5;
			j = 0;
		}

		if(fc < 0.0){
			fa = fa*0.5;
			j = 1;
		}


		if(j==1) fb = fc;
		if(j==0) fa = fc;
		count += 1;

	}

	if(count >= countmax){
		if(fabs(fb) <= errormax*epsin || fabs(fa) <= errormax*epsin){
			if(fabs(fb)<fabs(fa)) c = b;
			else c = a;
			*T = c;
			return 0;
		}
		else {
			printf("Convergence failure in compute_temp: Tmax = %le, Tmin = %le, eps(Tmax) = %le, eps(Tmin) = %le\n", Tmax, Tmin, epsmax,
					epsmin);
			printf("Last iteration: (a = %le, fa = %le), (b = %le, fb = %le), (c = %le, fc = %le)     errormax = %le\n",
				a, fa, b, fb, c, fc, errormax);
			*T = (fabs(fa) <= fabs(fb)) ? a : b;
		}



	}


	return 1;
}



// Find T given rho, Y and eps. Note that if eps is in the validity range, it almost certainly finds T,
// thus only use after setting eps and rho into the valid region

int find_temp(double *T, double *rho, double *Y, double *epsl, double *Tmin, double *Tmax, double *epsmin, double *epsmax){

 double Logr = log10(*rho);
 double LogT, dLogT, LogTn, eps0, eps1, eps, dummy, LogT1, T0;
 double LogTmin = log10(*Tmin);
 double LogTmax = log10(*Tmax);
 double tol = 1e-10, d1, d2, d3;
 int itmax = 50, i, RETURN;

 eps0 = *epsl;
 eps1 = eps0;
 
 //As an initial guess, adopt the current T
 T0 = *T;
 if(*T < pow(10.,LogTmin)) *T = *Tmin;
 if(*T > pow(10.,LogTmax)) *T = *Tmax;

 LogT = log10(*T);
 LogT1 = LogT;

 intp3d(LogT, Logr, *Y, &eps, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.eps, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &d1, &d2, &d3);

 // First check: is the temperature already correct?

 if(fabs(eps-eps0)<=tol*fabs(eps0)){
	return 0;
 }

 //If not, Newton-Raphson to find the root (eps - eps0). 1D root finding employed for fixed Y and rho

  for(i=0;i<itmax;i++){
	if(d1 == 0.){
		i = itmax;
		break;
	}
	dLogT = -(eps-eps0)/d1;
	LogTn = LogT + dLogT;

	//if(CheckForNANandINF(1,dLogT)) printf("i = %d, bad dLogT = %le, d1 = %le, -(eps-eps0) = %le\n", i, dLogT, d1, -(eps-eps0));

	LogTn = (LogTn<LogTmin) ? LogTmin:LogTn;
	LogTn = (LogTn>LogTmax) ? LogTmax:LogTn;
	
	LogT1 = LogT;
	LogT = LogTn;
	eps1 = eps;
	*T = pow(10.,LogT);

	intp3d(LogT, Logr, *Y, &eps, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
			tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.eps, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &d1, &d2, &d3);
	
	//printf("i = %d, T0 = %le, Logr = %le, LogT = %le, T = %le, Y = %le (LogTmin = %le, LogTmax = %le, epsmin = %le, epsmax = %le)   =>      eps = %le, del = %le, depsdlogT = %le \n",
          //              i, T0, Logr, LogT, *T, *Y, LogTmin, LogTmax, eps, *epsmin, *epsmax, fabs(eps-eps0), d1);

	if(fabs(eps-eps0)<=tol*fabs(eps0)) return 0;

	// Closer than 1e-02 to the root we change to the secant method (crappy tab derivatives)
	if(fabs(eps-eps0)<=1e-03*fabs(eps0)){
		d1 = (eps-eps1)/(LogT - LogT1 + 1e-40); //avoid dividing by zero
		//if(CheckForNANandINF(1,d1)) printf("d1 = %le inside secant! eps = %le, eps1 = %le, LogT = %le, LogT1 = %le\n", d1, eps, eps1, LogT, LogT1);
	}
	
 }

        if(i >= itmax){
			//if(CheckForNANandINF(1,LogT)) printf("Bad LogT! Logr = %le, Y = %le, LogT = %le, T = %le\n", Logr, *Y, LogT, *T);
			if((eps0 - eps) < 0.){
				*Tmax = *T;
				*epsmax = eps;
			}
			else {
				*Tmin = *T;
				*epsmin = eps;
			}

			RETURN = compute_temp(T, *rho, *Y, *epsl, *Tmin, *Tmax, *epsmin, *epsmax);
			if(!RETURN) return 0;
			else {
				printf("T_in = %le, Logr = %le, T = %le, Y = %le, epsin = %le -------> Tmin = %le, Tmax = %le, epsmin = %le, epsmax = %le\n",
					T0, Logr, *T, *Y, *epsl, *Tmin, *Tmax, *epsmin, *epsmax);
			}
        }

return 1;


}


/***************************************** Functions for EoS wrapper *********************************************************/

// Compute cs^2 from rho, Y, eps

void cs2_tab(double rho, double T, double Y, double dpdr_T, double dedr_T,
	     double dpdT_r, double dedT_r, double *cs2, double *dpdr_e, double *dpde_r){
// Note that we're not using the derivatives to compute cs2; instead, we use a tab version of it. Nevertheless this function
// computes useful derivatives for other functions
 double dummy, LogT=log10(T);
 double Logr = log10(rho);

 intp3d(LogT, Logr, Y, cs2, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.cs2, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

 *dpdr_e = dpdr_T - dpdT_r * dedr_T/dedT_r;
 *dpde_r = dpdT_r/dedT_r;

}

// Find p, eps given rho, T, Y

void peps_tab (double *rho, double *T, double *Y, double *p, double *eps,
               double *dpdr_T, double *dedr_T, double *dpdT_r, double *dedT_r) {

	double d3, dpdLogr, dedLogr, dpdLogT, dedLogT;
	double Logr = log10(*rho);
	double LogT = log10(*T);

	intp3d(LogT, Logr, *Y, p, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
			tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.p, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dpdLogT, &dpdLogr, &d3);
    intp3d(LogT, Logr, *Y, eps, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
			tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.eps, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dedLogT, &dedLogr, &d3);
	*dpdr_T = 1./(*rho*log(10.))*dpdLogr;
    *dedr_T = 1./(*rho*log(10.))*dedLogr;
	*dpdT_r = 1./(*T*log(10.))*dpdLogT;
	*dedT_r = 1./(*T*log(10.))*dedLogT;


}


/********** Wrapper functions *************/

//Find T, p, and cs2 given rho, eps and Y

int eos_tab3d (double *rho, double *eps, double *Y, double *p,
               double *T, double *cs2, double *dpdrho, double *dpdeps){

 double dedr_T, dpdr_T, dedT_r, dpdT_r, dummy;

 //WARNING: SHOULD INCLUDE EXTEND_VALIDITY HERE, BUT CURRENTLY WE'RE NOT USING THIS FUNCTION

//Find p given rho, T and Y
 peps_tab(rho, T, Y, p, &dummy, &dpdr_T, &dedr_T, &dpdT_r, &dedT_r);
 if(!finite(*p) || !finite(dpdr_T) || !finite(dedr_T) || !finite(dpdT_r) || !finite(dedT_r)) return 1;

//Find cs2 given rho, T and Y
 cs2_tab(*rho, *T, *Y, dpdr_T, dedr_T, dpdT_r, dedT_r, cs2, dpdrho, dpdeps);
 if(!finite(*cs2)) return 1;

 return 0;
}

//Find p and cs2 given rho, T and Y

int eos_tab3d_T(double *rho, double *eps, double *Y, double *p, double *T, double *cs2, double *dpdrho, double *dpdeps)
{

 double dpdr_T, dpdT_r, dedr_T, dedT_r, dummy;

//Find p given rho, T and Y
 peps_tab(rho, T, Y, p,  eps, &dpdr_T, &dedr_T, &dpdT_r, &dedT_r);
 if(!finite(*p) || !finite(dpdr_T) || !finite(dedr_T) || !finite(dpdT_r) || !finite(dedT_r)) return 1;

//Find cs2 given rho, T and Y
 cs2_tab(*rho, *T, *Y, dpdr_T, dedr_T, dpdT_r, dedT_r, cs2, dpdrho, dpdeps);
 if(!finite(*cs2)) return 1;

 return 0;
}

// Find chi = dpdrho and kappa = dpdeps given rho, T and Y

int eos_tab3d_T_D(double *rho, double *T, double *Y, double *chi, double *kappa)
{

 double dummy;
 double Logr = log10(*rho);
 double LogT = log10(*T);

        intp3d(LogT, Logr, *Y, chi, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.chi, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, *Y, kappa, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.kappa, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

	if(!finite(*chi) || !finite(*kappa)) return 1;

	return 0;
}

//Find p and cs2 given rho, eps and fixed Y

int eos_tab3d_Y(double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps)
{
 
 double dpdT_r, dedT_r, dpdr_T, dedr_T;
 double T, dummy;
 double Y = EOS.Y;
 
 //find_temp(&T, rho, &Y, eps);
 if(!finite(T)) return 1;

//Find p given rho, T and Y
 peps_tab(rho, &T, &Y, p, &dummy, &dpdr_T, &dedr_T, &dpdT_r, &dedT_r);
 if(!finite(*p) || !finite(dpdr_T) || !finite(dedr_T) || !finite(dpdT_r) || !finite(dedT_r)) return 1;

//Find cs2 given rho, T and Y
 cs2_tab(*rho, T, Y, dpdr_T, dedr_T, dpdT_r, dedT_r, cs2, dpdrho, dpdeps);
 if(!finite(*cs2)) return 1;

 return 0;
}

// Find p, eps and cs2 given rho, fixed Y and T = 0; this function is mostly used for atmosphere

int eos_tab3d_T0_Y(double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps)
{
 double dpdr_T, dpdT_r, dedr_T, dedT_r;
 double T = 0.0;
 double Y = EOS.Y;
 
 //Find p, eps and cs2 given rho, fixed Y and T = 0

 peps_tab(rho, &T, &Y, p, eps, &dpdr_T, &dedr_T, &dpdT_r, &dedT_r);
 if(!finite(*p) || !finite(dpdr_T) || !finite(dedr_T) || !finite(dpdT_r) || !finite(dedT_r)) return 1;

//Find cs2 given rho, T and Y
 cs2_tab(*rho, T, Y, dpdr_T, dedr_T, dpdT_r, dedT_r, cs2, dpdrho, dpdeps);
 if(!finite(*cs2)) return 1;

 return 0;

}

int eos_tab3d_micro (double *rho, double *T, double *Y, double *mun, double *mup, double *mue, double *A, double *Z, double *nh)
{
  double Logr = log10(*rho);
  double LogT = log10(*T);
  double dummy;

  intp3d(LogT, Logr, *Y, mun, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.mun, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
  intp3d(LogT, Logr, *Y, mup, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.mup, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
  intp3d(LogT, Logr, *Y, mue, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.mue, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
  intp3d(LogT, Logr, *Y, A, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Abar, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
  intp3d(LogT, Logr, *Y, Z, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Zbar, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
  intp3d(LogT, Logr, *Y, nh, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.nh, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

  if(!finite(*rho) || !finite(*T) || !finite(*Y) || !finite(*mun) || !finite(*mup) || !finite(*mue)
                || !finite(*A) || !finite(*Z) || !finite(*nh))
        {
                if (PR) printf("!finite: rho = %e -> T = %e, Y = %e, mu_n = %e, mu_p = %e, mu_e = %e, A = %e, Z = %e, n_h = %e\n",
                                *rho, *T, *Y, *mun, *mup, *mue, *A, *Z, *nh);
                return 1;
        }

  return 0;
}


/****** wrappers for 1D tab *******************************/


int peps_1d (double *p, double *cs2, double *dpdrho,  double *dpdeps, double *rho, double *eps)
{

        double dpdLogr, depsdLogr, tmp;
        double Logr = log10(*rho);
        int i=0;

	if(!finite(Logr)){
              	if (PR) printf("!finite Logr = %le for rho = %le\n", Logr, *rho);
                return 1;
        }

        splint(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.p, tab1d_hot_nuclear.d2pdr2, tab1d_hot_nuclear.N, Logr, p, &dpdLogr);
        splint(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.eps, tab1d_hot_nuclear.d2epsdr2, tab1d_hot_nuclear.N, Logr, eps, &depsdLogr);
	splint(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.cs2, tab1d_hot_nuclear.d2cs2dr2, tab1d_hot_nuclear.N, Logr, cs2, &tmp);


        *dpdrho = 1./(*rho*log(10.))*dpdLogr;
        *dpdeps = 0.;

        if (PR) printf("rho = %e  -> eps = %e  p = %e, dpdr = %e  dpdeps = %e -> cs2 = %e\n", *rho, *eps, *p, *dpdrho, *dpdeps, *cs2);
        if (!finite(*rho)|| !finite(*eps) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e, cs2 = %e\n", *rho, *eps, *p, *cs2);
                return 1;
        }
        return 0;

}

int beta_1d (double *rho, double *Y, double *T)
{

	double dummy;
        double Logr = log10(*rho);
        int i=0;

        splint(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.Y, tab1d_hot_nuclear.d2Ydr2, tab1d_hot_nuclear.N, Logr, Y, &dummy);
        splint(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.T, tab1d_hot_nuclear.d2Tdr2, tab1d_hot_nuclear.N, Logr, T, &dummy);

	
	if(!finite(*Y) || !finite(*T)) return 1;
	return 0;

}


int peps_1d_lin (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps)
{

	double dpdLogr, depsdLogr, tmp;
	double Logr = log10(*rho);
	int i=0;

	if(!finite(Logr)){
		if (PR) printf("!finite Logr = %le for rho = %le\n", Logr, *rho);
		return 1;
	}

	while(Logr > tab1d_hot_nuclear.Logr[i]) i++;
        if(i > 0){
                *p   = tab1d_hot_nuclear.p[i-1] + (tab1d_hot_nuclear.p[i] - tab1d_hot_nuclear.p[i-1])/
                                        (tab1d_hot_nuclear.Logr[i] - tab1d_hot_nuclear.Logr[i-1])*(Logr - tab1d_hot_nuclear.Logr[i-1]);
                *eps = tab1d_hot_nuclear.eps[i-1] + (tab1d_hot_nuclear.eps[i] - tab1d_hot_nuclear.eps[i-1])/(tab1d_hot_nuclear.Logr[i] - tab1d_hot_nuclear.Logr[i-1])*(Logr - tab1d_hot_nuclear.Logr[i-1]);
                *cs2 = tab1d_hot_nuclear.cs2[i-1] + (tab1d_hot_nuclear.cs2[i] - tab1d_hot_nuclear.cs2[i-1])/(tab1d_hot_nuclear.Logr[i] - tab1d_hot_nuclear.Logr[i-1])*(Logr - tab1d_hot_nuclear.Logr[i-1]);
                dpdLogr = (tab1d_hot_nuclear.p[i] - tab1d_hot_nuclear.p[i-1])/(tab1d_hot_nuclear.Logr[i] - tab1d_hot_nuclear.Logr[i-1]);
        }
        else {
                *p = tab1d_hot_nuclear.p[i];
                *eps = tab1d_hot_nuclear.eps[i];
                *cs2 = tab1d_hot_nuclear.cs2[i];
                 dpdLogr = (tab1d_hot_nuclear.p[i+1] - tab1d_hot_nuclear.p[i])/(tab1d_hot_nuclear.Logr[i+1] - tab1d_hot_nuclear.Logr[i]);
        }

        *dpdrho = 1./(*rho*log(10.))*dpdLogr;
        *dpdeps = 0.;

        if (PR) printf("rho = %e  -> eps = %e  p = %e, dpdr = %e  dpdeps = %e -> cs2 = %e\n", *rho, *eps, *p, *dpdrho, *dpdeps, *cs2);
        if (!finite(*rho)|| !finite(*eps) || !finite(*cs2)) {
                if(PR) printf("!finite, rho = %e, eps = %e, p = %e, cs2 = %e\n", *rho, *eps, *p, *cs2);
                return 1;
        }
        return 0;

}

int beta_1d_lin (double *rho, double *Y, double *T)
{

	double Logr = log10(*rho);
	int i=0;


        while(Logr > tab1d_hot_nuclear.Logr[i]) i++;
        if(i > 0){
                *Y   = tab1d_hot_nuclear.Y[i-1] + (tab1d_hot_nuclear.Y[i] - tab1d_hot_nuclear.Y[i-1])/(tab1d_hot_nuclear.Logr[i] - tab1d_hot_nuclear.Logr[i-1])*(Logr - tab1d_hot_nuclear.Logr[i-1]);
                *T = tab1d_hot_nuclear.T[i-1] + (tab1d_hot_nuclear.T[i] - tab1d_hot_nuclear.T[i-1])/(tab1d_hot_nuclear.Logr[i] - tab1d_hot_nuclear.Logr[i-1])*(Logr - tab1d_hot_nuclear.Logr[i-1]);
        }
        else {
                *Y = tab1d_hot_nuclear.Y[i];
                *T = tab1d_hot_nuclear.T[i];
        }


        if(!finite(*Y) || !finite(*T)) return 1;

        return 0;

}

int peps_1d_steffen(double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps)
{

	double Logr = log10(*rho);
	double dpdLogr, depsdLogr, epsl;

        if(!finite(Logr)) {
                if(PR) printf("!finite Logr = %e for rho = %e\n", Logr, *rho);
                return 1;
        }

	steffen_intp(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.p, tab1d_hot_nuclear.ap, tab1d_hot_nuclear.bp, tab1d_hot_nuclear.cp, tab1d_hot_nuclear.dp, tab1d_hot_nuclear.N, Logr, p, &dpdLogr);
        steffen_intp(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.eps, tab1d_hot_nuclear.aeps, tab1d_hot_nuclear.beps, tab1d_hot_nuclear.ceps, tab1d_hot_nuclear.deps, tab1d_hot_nuclear.N, Logr, &epsl, &depsdLogr);
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

int beta_1d_steffen(double *rho, double *Y, double *T)
{

	double dummy;
	double Logr = log10(*rho);

	if(!finite(Logr)) {
                if(PR) printf("!finite Logr = %e for rho = %e\n", Logr, *rho);
                return 1;
        }


        steffen_intp(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.Y, tab1d_hot_nuclear.aY, tab1d_hot_nuclear.bY, tab1d_hot_nuclear.cY, tab1d_hot_nuclear.dY, tab1d_hot_nuclear.N, Logr, Y, &dummy);
        steffen_intp(tab1d_hot_nuclear.Logr, tab1d_hot_nuclear.T, tab1d_hot_nuclear.aT, tab1d_hot_nuclear.bT, tab1d_hot_nuclear.cT, tab1d_hot_nuclear.dT, tab1d_hot_nuclear.N, Logr, T, &dummy);
        if(!finite(*Y)||!finite(*T)) return 1;

        return 0;

}

/****** Neutrino interactions *******************************/

int eos_tab3d_mu (double Logr, double LogT, double Y, double *mun, double *mup, double *mue)
{

  double dummy;

  intp3d(LogT, Logr, Y, mun, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.mun, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
  intp3d(LogT, Logr, Y, mup, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.mup, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
  intp3d(LogT, Logr, Y, mue, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.mue, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

  if(!finite(Logr) || !finite(LogT) || !finite(Y) || !finite(*mun) || !finite(*mup) || !finite(*mue))
        {
                if (PR) printf("!finite: rho = %e -> T = %e, Y = %e, mu_n = %e, mu_p = %e, mu_e = %e",
                                Logr, LogT, Y, *mun, *mup, *mue);
                return 1;
        }

  return 0;
}

int eos_tab3d_eq (double Logr, double LogT, double Y, double *Yeq, double *Teq)
{
  double dummy;

    intp3d(LogT, Logr, Y, Yeq, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Yeq, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
    intp3d(LogT, Logr, Y, Teq, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
		tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Teq, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

    if(!finite(Logr) || !finite(LogT) || !finite(Y) || !finite(*Yeq) || !finite(*Teq)) {
        if (PR) printf("!finite: rho = %e -> T = %e, Y = %e, Yeq = %e, Teq = %e", Logr, LogT, Y, *Yeq, *Teq);
        return 1;
    }

  return 0;
}

void nls_table(double **Qe, double **Qa, double **Qx, double **Re, double **Ra, double **Rx,
	       double **zetae, double **zetaa, double **zetax, double **eta_nue,
	       double **ka_0_nue, double **ka_1_nue, double **ks_nue, double **ka_0_nua, double **ka_1_nua, double **ks_nua,
	       double **ka_0_nux, double **ka_1_nux, double **ks_nux)
{

    tU units;
    set_units(&units);

    int P, i, j, k, ijk;
    double nb, T, Y, np, nn, nh, eta_e, eta_np, eta_pn, eta_hat, Z, A;
    double *Qloc, *Rloc;
    double F4, F5, F3, F2;
    double scatt_cons, abs_cons, me=0.510998910, sigma_0=1.705e-44, alp=1.23, block_nue, block_nua;
    double Xp, Xn, R_brem, Q_brem, fac1, fac2, fac3, fac4, cp, cn;
    double h = units.hplanck_cgs*units.Energy_MeV/units.Energy_cgs;  // Planck const in MeV.s
    double c = units.clight_cgs;                     // speed of light in cm/s
    scatt_cons = 0.25*sigma_0/(me*me);
    abs_cons = (1.+3.*alp*alp)*scatt_cons;
    cn = (1. + 5.*alp*alp)/6.;
    cp = (0.0064 + 5.*alp*alp)/6.;
    const double g = 4; // degeneracy factor of heavy neutrinos

    printf("\n Constructing NLS free emission and mean free paths table\n");

    for (P=0; P<=bampi_size();P++){

    		printf("Reading on proc%d/%d\n", P, bampi_size());

    	if(bampi_rank()==P){

                Qloc = malloc(3*sizeof(double));
                Rloc = malloc(3*sizeof(double));
				(*Qe) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
				(*Qa) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
				(*Qx) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                (*zetae) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                (*zetaa) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                (*zetax) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                (*Re) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                (*Ra) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                (*Rx) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
				(*eta_nue) = (double*) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));

		if(Getv("grrhd_m1_nls_rates", "yes")) {

                	(*ka_0_nue) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ka_1_nue) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ks_nue) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ka_0_nua) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ka_1_nua) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ks_nua) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ka_0_nux) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ka_1_nux) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
                	(*ks_nux) = (double *) malloc((tab3d_nuclear.Nr*tab3d_nuclear.NT*tab3d_nuclear.NY)*sizeof(double));
		}

		for(i=0;i<tab3d_nuclear.NT;i++){
			T = pow(10.,tab3d_nuclear.LogT[i]);
			for(j=0;j<tab3d_nuclear.Nr;j++){
				nb = pow(10.,tab3d_nuclear.Logr[j])*units.Mdens_cgs/tab3d_nuclear.mb;
				for(k=0;k<tab3d_nuclear.NY;k++){
					ijk = k + tab3d_nuclear.NY*j + tab3d_nuclear.NY*tab3d_nuclear.Nr*i;
					np = tab3d_nuclear.np[ijk]; nn = tab3d_nuclear.nn[ijk];
					nh = tab3d_nuclear.nh[ijk]; A = tab3d_nuclear.Abar[ijk]; Z = tab3d_nuclear.Zbar[ijk];
					eta_hat = (tab3d_nuclear.mun[ijk] - tab3d_nuclear.mup[ijk]-1.2935)/T; //Subtract mass difference here
					eta_np = (nn - np)/(exp(eta_hat)-1.);
					eta_np = DMAX(eta_np, 0.);
					eta_pn = (np - nn)/(exp(-eta_hat)-1.);
					eta_pn = DMAX(eta_pn, 0.);
					if((nb*tab3d_nuclear.mb)<2e+11){
						eta_np = np;
						eta_pn = nn;
					}
					eta_e = tab3d_nuclear.mue[ijk]/T;
					(*eta_nue)[ijk] = eta_e + tab3d_nuclear.mup[ijk]/T - tab3d_nuclear.mun[ijk]/T;
					compute_local_emission(nb, T, eta_np, eta_pn, eta_e, (*eta_nue)[ijk], &(Qloc), &(Rloc));

                                        (*Qe)[ijk] = Qloc[0];
                                        (*Qa)[ijk] = Qloc[1];
                                        (*Qx)[ijk] = Qloc[2];
                                        (*Re)[ijk] = Rloc[0];
                                        (*Ra)[ijk] = Rloc[1];
                                        (*Rx)[ijk] = Rloc[2];
					if(Getv("nls_use_NNBrem","yes")){
					  // Include NN bremsstrahlung; masses took from CompOSE manual
					  Xp = (938.2720813/units.Energy_MeV*np*units.Volume_cgs)/(pow(10., tab3d_nuclear.Logr[j]));
					  Xn = (939.565413/units.Energy_MeV*nn*units.Volume_cgs)/(pow(10., tab3d_nuclear.Logr[j]));
					  fac1 = 0.231*(2.07782e+02/units.Energy_cgs*units.Energy_MeV)*0.5;
				 	  fac2 = (Xn*Xn + Xp*Xp + 28./3.*Xn*Xp)*pow(pow(10., tab3d_nuclear.Logr[j])*units.Mdens_cgs, 2.);
					  R_brem = fac1*fac2*pow(T, 4.5)/nb;
					  Q_brem = R_brem*T/0.231*0.504;

					  (*Qe)[ijk] += Q_brem;
					  (*Qa)[ijk] += Q_brem;
					  (*Qx)[ijk] += 4.*Q_brem;
					  (*Re)[ijk] += R_brem;
					  (*Ra)[ijk] += R_brem;
					  (*Rx)[ijk] += 4.*R_brem;
					}
					// Energy-independent mean free path
                			F4 = fermi_integral(4,(*eta_nue)[ijk]);
                			F5 = fermi_integral(5,(*eta_nue)[ijk]);
                			if (F4 == 0.) block_nue = 1;
                			else          block_nue = 1./(1.+exp(-F5/F4+eta_e));

					// We (and the remaining of the Milky way inhabitants)
					// compute opacities as of Ruffert et al. (1996)
					if(Getv("grrhd_m1_nls_rates", "yes")){
					  F2 = fermi_integral(2, (*eta_nue)[ijk]);
					  F3 = fermi_integral(3, (*eta_nue)[ijk]);
					  if(F2==0.) {
					   fac1 = 12.*(1.+0.1092*exp(0.8908*(*eta_nue)[ijk]))/(1.+0.0287*exp(0.9257*(*eta_nue)[ijk]));
					  } else fac1 = F4/F2;

					  (*ka_0_nue)[ijk] = abs_cons*eta_pn*block_nue*T*T*fac1;

					  if(F3==0.){
					    fac2 = 20.*(1.+0.0559*exp(0.9069*(*eta_nue)[ijk]))/(1.+0.0147*exp(0.9431*(*eta_nue)[ijk]));
					  } else fac2 = F5/F3;

					  (*ka_1_nue)[ijk] = abs_cons*eta_pn*block_nue*T*T*fac2;

					}

                			// This hack is needed to avoid nans in block_nua since for large -eta_nue, F4 -> 0, F5 ->0
                			F4 = fermi_integral(4,-(*eta_nue)[ijk]);
                			F5 = fermi_integral(5,-(*eta_nue)[ijk]);
                			if(F4 == 0.) block_nua = 1.; // this is the correct limit
                			else         block_nua = 1./(1.+exp(-F5/F4-eta_e)); //here eta_nua = -eta_nue 

					if(Getv("grrhd_m1_nls_rates", "yes")){
                                          F2 = fermi_integral(2, -(*eta_nue)[ijk]);
                                          F3 = fermi_integral(3, -(*eta_nue)[ijk]);
                                          if(F2==0.) {
                                           fac3 = 12.*(1.+0.1092*exp(-0.8908*(*eta_nue)[ijk]))/(1.+0.0287*exp(-0.9257*(*eta_nue)[ijk]));
					  } else fac3 = F4/F2;

                                          (*ka_0_nua)[ijk] = abs_cons*eta_np*block_nua*T*T*fac3;

                                          if(F3==0.){
                                            fac4 = 20.*(1.+0.0559*exp(-0.9069*(*eta_nue)[ijk]))/(1.+0.0147*exp(-0.9431*(*eta_nue)[ijk]));
					  } else fac4 = F5/F3;

					  (*ka_1_nua)[ijk] = abs_cons*eta_np*block_nua*T*T*fac4;
					}

					

                	(*zetae)[ijk] = np*scatt_cons*cp + nn*scatt_cons*cn + nn*abs_cons*block_nue;
            		(*zetaa)[ijk] = np*scatt_cons*cp + nn*scatt_cons*cn + np*abs_cons*block_nua;
        			(*zetax)[ijk] = (np*cp + nn*cn)*scatt_cons;
                	if(A!=0.){
                        (*zetae)[ijk] += nh*scatt_cons/4.*A*(1.-Z/A)*(1.-Z/A);
                        (*zetaa)[ijk] += nh*scatt_cons/4.*A*(1.-Z/A)*(1.-Z/A);
            			(*zetax)[ijk] += nh*scatt_cons/4.*A*(1.-Z/A)*(1.-Z/A);
                	}

					// Now the scattering opacities
					// Also in this conditional we set the abs opacities of nux according to
					// Kirchhoff's law
					if(Getv("grrhd_m1_nls_rates", "yes")){
						(*ks_nue)[ijk] = (np*scatt_cons*cp + nn*scatt_cons*cn)*T*T*fac2;
						(*ks_nua)[ijk] = (np*scatt_cons*cp + nn*scatt_cons*cn)*T*T*fac4;
						(*ks_nux)[ijk] = (np*cp + nn*cn)*scatt_cons*T*T*20.8120626786;
						// numerical factor on ks_nux is F5/F3 for eta = 0.;
						if(A!=0.){
						   (*ks_nue)[ijk] += nh*scatt_cons/4.*A*(1.-Z/A)*(1.-Z/A)*T*T*fac2;
						   (*ks_nua)[ijk] += nh*scatt_cons/4.*A*(1.-Z/A)*(1.-Z/A)*T*T*fac4;
						   (*ks_nux)[ijk] += nh*scatt_cons/4.*A*(1.-Z/A)*(1.-Z/A)*T*T*20.8120626786;
						}
						
						// Opacities in cm^-1
						(*ka_0_nux)[ijk] = (*Rx)[ijk]/(g*4.*c*M_PI*pow(T,3.)*fermi_integral(2,0.)/pow(h*c,3.));
						(*ka_1_nux)[ijk] = (*Qx)[ijk]/(g*4.*c*M_PI*pow(T,4.)*fermi_integral(3,0.)/pow(h*c,3.));
						//(*ka_0_nux)[ijk] = 0.;
						//(*ka_1_nux)[ijk] = 0.;

						// Finally, let us set everything to BAM units
						(*Qe)[ijk] *= (units.Volume_cgs*units.Time_cgs)/units.Energy_MeV;
                        (*Qa)[ijk] *= (units.Volume_cgs*units.Time_cgs)/units.Energy_MeV;
                        (*Qx)[ijk] *= (units.Volume_cgs*units.Time_cgs)/units.Energy_MeV;
                        (*Re)[ijk] *= (units.Volume_cgs*units.Time_cgs);
                        (*Ra)[ijk] *= (units.Volume_cgs*units.Time_cgs);
                        (*Rx)[ijk] *= (units.Volume_cgs*units.Time_cgs);
						(*ka_0_nue)[ijk] *= units.Length_cgs;
						(*ka_1_nue)[ijk] *= units.Length_cgs;
						(*ks_nue)[ijk]   *= units.Length_cgs;
                        (*ka_0_nua)[ijk] *= units.Length_cgs;
                        (*ka_1_nua)[ijk] *= units.Length_cgs;
                        (*ks_nua)[ijk]   *= units.Length_cgs;
                        (*ka_0_nux)[ijk] *= units.Length_cgs;
                        (*ka_1_nux)[ijk] *= units.Length_cgs;
                        (*ks_nux)[ijk]   *= units.Length_cgs;
					}

					if(CheckForNANandINF(9, (*Qe)[ijk], (*Qa)[ijk], (*Qx)[ijk],
					(*Re)[ijk], (*Ra)[ijk], (*Rx)[ijk],
					(*zetae)[ijk], (*zetaa)[ijk], (*zetax)[ijk])) 
					printf("Qe = %le, Qa = %le, Qx = %le, Re = %le, Ra = %le, Rx = %le, ze = %le, za = %le, zx = %le\n",
                                                (*Qe)[ijk], (*Qa)[ijk], (*Qx)[ijk], (*Re)[ijk], (*Ra)[ijk], (*Rx)[ijk], (*zetae)[ijk],
                                                (*zetaa)[ijk], (*zetax)[ijk]);

					if(Getv("grrhd_m1_nls_rates","yes")){
					  if(CheckForNANandINF(9, (*ka_0_nue)[ijk], (*ka_1_nue)[ijk], (*ks_nue)[ijk],
                                             (*ka_0_nua)[ijk], (*ka_1_nua)[ijk], (*ks_nua)[ijk],
                                             (*ka_0_nux)[ijk], (*ka_1_nux)[ijk], (*ks_nux)[ijk])){
                                          printf("ka_0_nue = %le, ka_1_nue = %le, ks_nue = %le\n", (*ka_0_nue)[ijk], (*ka_1_nue)[ijk], (*ks_nue)[ijk]);
					  printf("ka_0_nua = %le, ka_1_nua = %le, ks_nua = %le\n", (*ka_0_nua)[ijk], (*ka_1_nua)[ijk], (*ks_nua)[ijk]);
					  printf("ka_0_nux = %le, ka_1_nux = %le, ks_nux = %le\n", (*ka_0_nux)[ijk], (*ka_1_nux)[ijk], (*ks_nux)[ijk]);

					  }
					}
				}

			}

		}

		free(Qloc); free(Rloc);
	}}

}


int nls_free(double rho, double T, double Y, double *Qe, double *Qa, double *Qx,
	     double *Re, double *Ra, double *Rx)
{

	double dummy;
	double LogT = log10(T);
	double Logr = log10(rho);

	//printf("T = %le, rho = %le, Y = %le\n", T, rho, Y);

	intp3d(LogT, Logr, Y, Qe, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Qe, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, Qa, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Qa, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, Qx, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Qx, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, Re, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Re, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, Ra, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Ra, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, Rx, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Rx, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

  	if(!finite(Logr) || !finite(LogT) || !finite(Y) || !finite(*Qe) || !finite(*Qa) || !finite(*Qx))
        {
                if (PR) printf("!finite: rho = %e -> T = %e, Y = %e, Qe = %e, Qa = %e, Qx = %e",
                                Logr, LogT, Y, *Qe, *Qa, *Qx);
                return 1;
        }

  	return 0;


}

int nls_diff (double rho, double T, double Y, double chi_e, double chi_a, double chi_x,
	      double *Qdiff, double *Rdiff, double *leak_tau)
{

	tU units;
	set_units(&units);
        double h = units.hplanck_cgs*units.Energy_MeV/units.Energy_cgs;  // Planck const in MeV.s
        double c = units.clight_cgs;  
	double R_cons, dummy, fac, nb; 
	double zeta[1][3]; //indices: [points]: 0 -> central, 1,2 -> -x, +x; 3,4 -> -y, +y; 5, 6 -> -z,+z; [flavor]: 0, 1, 2 -> e,a,x
	double chi_min[3];
	double eta_eq[3];
	double F2, F4;
	int i;

	intp3d(log10(T), log10(rho), Y, &eta_eq[0], tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.eta_nue, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(log10(T), log10(rho), Y, &(zeta[0][0]), tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ze, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(log10(T), log10(rho), Y, &(zeta[0][1]), tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.za, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(log10(T), log10(rho), Y, &(zeta[0][2]), tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                                tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.zx, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

	//Compute emission rates and optical depth (leak_tau)

	R_cons = 4.*M_PI*c/(6.*pow((h*c),3.));

	eta_eq[1] = -eta_eq[0];
	eta_eq[2] = 0.0;

	chi_min[0] = chi_e;
	chi_min[1] = chi_a;
	chi_min[2] = chi_x;

	for(i=0;i<3;i++){
		Rdiff[i] = R_cons*zeta[0][i]/(chi_min[i]*chi_min[i])*T*fermi_integral(0,eta_eq[i]);
		Qdiff[i] = R_cons*zeta[0][i]/(chi_min[i]*chi_min[i])*T*T*fermi_integral(1,eta_eq[i]);
		//Again, be careful of the fermi integrals ratio for negative degeneracies
                F2 = fermi_integral(2,eta_eq[i]);
                F4 = fermi_integral(4,eta_eq[i]);
                if(F2 == 0.) fac = 12.*(1.+0.1092*exp(0.8908*eta_eq[i]))/(1.+0.0287*exp(0.9257*eta_eq[i]));
                else fac = F4/F2;
                leak_tau[i] = chi_min[i]*T*T*fac;
		//For heavy lepton neutrinos, multiply by degeneracy factor g = 4
		if(i==2){
			Rdiff[i] = 4.*Rdiff[i];
			Qdiff[i] = 4.*Qdiff[i];
		}

		
	}

	return 0;

}

int nls_zeta (double rho, double T, double Y, double *zetae, double *zetaa, double *zetax)
{

	double LogT = log10(T);
	double Logr = log10(rho);
	double dummy;

	intp3d(LogT, Logr, Y, zetae, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ze, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, zetaa, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.za, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, zetax, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY, 
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.zx, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

	return 0;


}

int eos_m1_rates (int ind, double rho, double T, double Y, double *Q, double *R, double *ka_0,
             double *ka_1, double *ks)
{


	double LogT = log10(T);
	double Logr = log10(rho);
	double dummy;

	if(ind==0){
        intp3d(LogT, Logr, Y, Q, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Qe, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, R, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Re, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ka_0, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ka_0_nue, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ka_1, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ka_1_nue, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ks, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ks_nue, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
	}
	else if(ind==1){
        intp3d(LogT, Logr, Y, Q, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Qa, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, R, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Ra, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ka_0, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ka_0_nua, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ka_1, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ka_1_nua, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ks, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ks_nua, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

	}
	else {
        intp3d(LogT, Logr, Y, Q, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Qx, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, R, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.Rx, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ka_0, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ka_0_nux, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ka_1, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ka_1_nux, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);
        intp3d(LogT, Logr, Y, ks, tab3d_nuclear.NT, tab3d_nuclear.Nr, tab3d_nuclear.NY,
                        tab3d_nuclear.dLogT, tab3d_nuclear.dLogr, tab3d_nuclear.dY, tab3d_nuclear.ks_nux, tab3d_nuclear.LogT, tab3d_nuclear.Logr, tab3d_nuclear.Y, &dummy, &dummy, &dummy);

	}
}

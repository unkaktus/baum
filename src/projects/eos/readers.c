#include "bam.h"
#include "eos.h"

#define PR 0

/* Read tab1D hot: */


void eos_read_tab1d_hot(char *fname, int *n, double **x1, double **x2, double **x3,
			double **v1, double **v2, double **v3, double **v4,
			double **v5, double **v6, double **v7, double **v8,
			double **v9, double *mb, int key)
{

	int N, i, j, P;
	char line[1024];

	tU units;
	set_units (&units);
	double m;

	printf("\nReading tab1D hot EoS: \n%s\n", fname);

	// Read for all processors


	for (P=0; P<=bampi_size();P++){

		printf("Reading on proc%d/%d\n", P, bampi_size());

		if(bampi_rank()==P){


			FILE *f;
			f = fopen(fname, "r");

			// Count lines
			N = 0;

			while(fgets(line,1024,f)) N = N+1;

			rewind(f);
			double **t_temp = malloc(N*sizeof(double*));
			for(i=0;i<N;i++) t_temp[i] = malloc(13*sizeof(double));

			for(i=0;i<N;i++){
				for(j=0;j<13;j++){
					fscanf(f, "%le", &t_temp[i][j]);
					if (PR) printf("%le ", t_temp[i][j]);
				}
			if(PR) printf("\n");
			}
		
			// Find baryonic mass constant mb

			m = t_temp[0][8];

			for(i=1;i<N;i++){
				if(t_temp[i][8] <= m) m = t_temp[i][8];
			}

			if(key) *mb = m/units.Energy_MeV*units.Mass_cgs;
			else *mb = (*mb < m) ? *mb : m;

			//Setting and filling the vectors

			(*x1) = (double *) malloc(N*sizeof(double));
			(*x2) = (double *) malloc(N*sizeof(double));
			(*x3) = (double *) malloc(N*sizeof(double));
			(*v1) = (double *) malloc(N*sizeof(double));
			(*v2) = (double *) malloc(N*sizeof(double));
			(*v3) = (double *) malloc(N*sizeof(double));
                        (*v4) = (double *) malloc(N*sizeof(double));
                        (*v5) = (double *) malloc(N*sizeof(double));
                        (*v6) = (double *) malloc(N*sizeof(double));
                        (*v7) = (double *) malloc(N*sizeof(double));
                        (*v8) = (double *) malloc(N*sizeof(double));
                        (*v9) = (double *) malloc(N*sizeof(double));

			double mn = 939.565413;	// Neutron mass in MeV

			for(i=0;i<N;i++){
				(*x1)[i] = t_temp[i][0]; // Temperature in MeV
				(*x2)[i] = log10(t_temp[i][1]*(*mb)*1e+39/units.Mdens_cgs); //log10(rest-mass density)
				(*x3)[i] = t_temp[i][2];	// Charge fraction
				(*v1)[i] = t_temp[i][3]/units.Energy_MeV*units.Volume_fm3; //Pressure
				(*v2)[i] = (t_temp[i][8]*units.Mass_cgs/units.Energy_MeV)/(*mb) - 1.0; // Specific internal energy (dimensionless)
				(*v3)[i] = t_temp[i][7]; // Speed of sound in c = 1
				(*v4)[i] = t_temp[i][4] + mn; // Neutrons chemical potential in MeV
				(*v5)[i] = t_temp[i][5]+(*v4)[i]; // Protons chemical potential in MeV
				(*v6)[i] = t_temp[i][6] - t_temp[i][5]; // Electrons chemical potential in MeV
				(*v7)[i] = t_temp[i][10]; // Avg heavy nuclei mass number Abar (dimensionless)
				(*v8)[i] = t_temp[i][11]; // Avg heavy nuclei charge number Zbar (dimensionless)
				(*v9)[i] = t_temp[i][9]*t_temp[i][1]*1e+39; //Avg heavy nuclei number density (cm^-3)
			}
			
			fclose(f);

			for(i=0;i<N;i++) free(t_temp[i]);
			free(t_temp);

		}


	}

	*n = N;
	printf("N = %d\n", N);
	printf("mb = %.16le \n", *mb);
	printf("\nEoS tab1D_hot loaded successfully from file \n%s\n", fname);

}


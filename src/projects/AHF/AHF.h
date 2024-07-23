/* AHF.h */
/* Jose Gonzalez 03/07 */

static const int AHF_Nmax = 2;
void SetAHF(int N, double t,double x,double y,double z, double r, double mass);



void AHF_1d_1(tVarList *vars, tVarList *ahfvars);
void AHF_1d_2(tVarList *vars, tVarList *ahfvars);

void AHF_output(int ntheta, double dtheta, int nphi, double dphi, 
                double xc, double yc, double zc, 
                double ***rr, char **outdir, int number, double time, double ahf_time);
void AHF_output_xyt(int ntheta, double dtheta, int nphi, double dphi, 
                    double xc, double yc, double zc, 
                    double ***rr, char **outdir, int number, double time, double ahf_time);

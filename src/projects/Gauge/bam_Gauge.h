/* bam_Gauge.h */
/* Bernd Bruegmann, 7/03 */

void add_betarot(tVarList *ucur, double sign);





#define generic_Gauge_VARS \
    double alpha, double betax, double betay, double betaz,\
    double Bx, double By, double Bz,\
    double dadx, double dady, double dadz,\
    double dbxdx, double dbxdy, double dbxdz, \
    double dbydx, double dbydy, double dbydz,\
    double dbzdx, double dbzdy, double dbzdz,\
    double advalpha, double advbetax, double advbetay, double advbetaz,\
    double advBx, double advBy, double advBz,\
    double Gx, double Gy, double Gz, double K,\
    double rGx, double rGy, double rGz,\
    double advGx, double advGy, double advGz,\
    double *ralpha, double *rbetax, double *rbetay, double *rbetaz,\
    double *rBx, double *rBy, double *rBz,\
    int *firstcall

void set_generic_Gauge();
extern void (*generic_Gauge)(generic_Gauge_VARS);


/* moving_punctures */
int Gauge_startup(tL *level);
int track_puncture(tL *level);
int track_moving_puncture_extremum(tL *level, int np,
                                   double *pold, double *pnew);


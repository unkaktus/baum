/* Gauge.h */
/* Bernd Bruegmann, 4/03 , Wolfgang Tichy 12/2003*/


/* from src/iterative/ */
void DPflatlinear(tL *level, tVarList *v, tVarList *u);

/* alpha_center.c */
int write_lapse_between_punctures(tL *level);

/* corotation.c */
int initialize_betarot(tL *level);
void adjust_beta(tG *g, double domega, double drdot);
int monitor_bh_asymmetry(tL *level);
int monitor_bh_lapse(tL *level);

/* commotion.c */
int commotion_instantaneous(tL *level);
int commotion_manual(tL *level);
int commotion_reiterate(tL *level);

/* maximal slicing */
int maximal(tL *level);
void maximal_L(tL *level, tVarList *vlLu, tVarList *vlu,
	       tVarList *vlgi, tVarList *vlGamma, tVarList *vlKK);
void maximal_init(tL *level);

/* moving_punctures.c */
int track_puncture(tL *level);
int compute_moving_punc(tL *level);

/* moving_punctures_integrate.c */
int track_moving_puncture_integrate(tL* level, int np, double *pold, double *pnew);
void compute_moving_puncture_distance(tL *level, double *p1, double *p2, double *length, double *clength);

/* moving_punctures_spin.c */
int moving_puncture_spin(tL *level);

/* scale_shift.c */
int scale_radial_shift(tG *g);
int addto_radial_shift(tG *g);
int adjust_radial_shift(tG *g);

extern int pintoxyplane;
extern int trackmethod;


/* puncture_spin by integrating over coordinate sphere */
int puncture_properties(tL *level);
void puncture_properties_integrand(tL *level, int i_x, int i_g, int i_K,
                                   int i_TrK, double x0,double y0, double z0);


/* old stuff of BB */
struct ps_surface_struct {
  int NumOfPunctures;
  double SearchTimeStart;
  double SearchTimeEnd;
  int *ListPunctures;
};

int pspin_IsMultipleTime(const double time, const double steptime);

int pspin_FileExists(const char * filename);

int pspin_ComputeSlice(const double currenttime, const double timeInverval);

void pspin_GetPunctureLocation(const int number, double *xp, double *yp, double *zp);

double pspin_GetPunctureMass(const int number);

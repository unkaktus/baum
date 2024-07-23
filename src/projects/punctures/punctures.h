/* punctures.h */
/* Bernd Bruegmann 6/02, Wolfgang Tichy 6/2003 */



/* Masses structure,
   originally in projects/punctures/bam_punctures.h [MDH, 11/05] */
typedef struct tmasses
{
  double M_ADM_inf;
  double M_ADM_1;
  double M_ADM_2;
  double M_K_inf;
  double M_K_1;
  double M_K_2;
  double J_ADM;
  double Omega;
  double S1[4];  /* S1[0] is magnitude of S1     */
  double S2[4];  /* S2[i] is x^i-component of S2 */
} tMasses;

extern tMasses GlobalMass;



void PickEllipticSolver(void);
int FlagRegridPunctures(tL *level);

int PunctureData(tL *level);
void PunctureCleanup(tL *level);
void PunctureSpecialInit(tL *level);
void PunctureSolve(tL *level);
void PunctureMaximalSlicing(tL *level);
void PunctureShift(tL *level);

void PuncturesPS(tL *level);

void SetMWBL(tL *level);
void SetBYK(tL *level);
void SetPunctureCoefficients(tL *level);
void ConfToPhys(tL *level, char *use_u);
double PunctureMass(tL *level, int fileoutput, tMasses *masses);

void VectorLaplaceFlat(tL *level, tVarList *vllu, tVarList *vlu);
void VectorLieFlat(tL *level, tVarList *vllu, tVarList *vlu);
void VectorLaplaceFlatGS(tL *level, tVarList *vlv, tVarList *vlu, 
			 tVarList *vlf);

void ThinPunctureLapse(tL *level);
void ThinPunctureShift(tL *level);
void LThinPunctures(tL *level, tVarList *vllv, tVarList *vlv, tVarList *vlu,
		    tVarList *vla, tVarList *vls, 
		    tVarList *vlp, tVarList *vldpop, tVarList *vlddpop);
void LThinPunctures2(tL *level, tVarList *vllv, tVarList *vlv, tVarList *vlu,
		     tVarList *vla, tVarList *vls,
		     tVarList *vlpsi, tVarList *vldpsiopsi,
		     tVarList *vlp, tVarList *vldpop, tVarList *vlddpop);

extern int (*punctures_solver)
     (tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c, 
      int imax, double tol, double *res,
      void (*lop)(tL *, tVarList *, tVarList *), 
      void (*precon)(tL *, tVarList *, tVarList *));

extern int (*linear_solver)
     (tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c, 
      int imax, double tol, double *res,
      void (*lop)(tL *, tVarList *, tVarList *), 
      void (*precon)(tL *, tVarList *, tVarList *));


/* wrapper for HYPRE lin. Solver to be called by Newton */
int Punc_hypreWrapper(tL *level, tVarList *x, tVarList *b, tVarList *r, 
                        tVarList *c, int itmax, double tol, double *res,
                        void (*Atimes)(tL *, tVarList *, tVarList *),
                        void (*precon)(tL *, tVarList *, tVarList *));
int Punc_hypreCleanup(tL *level);


/* from src/iterative/ */
void DPflatlinear(tL *level, tVarList *v, tVarList *u);
void Lflatlinear(tL *level, tVarList *vllu, tVarList *vlu);


#define EPSZERO (dequaleps/1000.0)

#define N_MWBL 2
extern double MBL[N_MWBL+1], CBL[N_MWBL+1][3];
extern double PBL[N_MWBL+1][3], SBL[N_MWBL+1][3];


/* we may want to use transformed coordinates */
#if 0
void coordtrans(double *x, double *y, double *z);
void indextrans_ddsym(double x, double y, double z, 
		      double *txx, double *txy, double *txz, 
		      double *tyy, double *tyz, double *tzz);
#else
#define coordtrans(x, y, z)
#define indextrans_ddsym(x, y, z, txx, txy, txz, tyy, tyz, tyzz)
#endif

extern int fisheye;
extern double fisheyedx, fisheyedy, fisheyedz;


/* Functions to compute Masses outside grid box */
void SetFuncCoeffs(double c_psi, double c_alpha, double c_alphapsi);
double RHS_of_Lapl_punctures_u(double x, double y, double z); 
double RHS_of_Lapl_punctures_v(double x, double y, double z);
double OneOver_a(double x, double y, double z);
double RHS_of_Lapl_alpha(double x, double y, double z);

/* Functions to make Sequences */
void Make_MK_eq_MADM_Sequence(tL *level);
int ComputePunctureParametersFromFit(tL *level);
int ComputePunctureParametersFromPN(tL *level);

/* Functions to prepare moving punctures */
int Absorb_psi(tL *level);
void Move_psi_intoMetric(tL* level, int i_g, 
                         int i_psi, int i_dpsiopsi, int i_ddpsiopsi);
void Set_alpha_psip(tL *level, int psipower);

/* Functions to to acchieve a given ADM mass at each puncture */
int PunctureData_with_given_M_ADM(tL *level);

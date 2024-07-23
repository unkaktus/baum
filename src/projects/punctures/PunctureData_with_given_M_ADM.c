/* PunctureData_with_given_M_ADM.c */
/* Wolfang Tichy 5/2006 */



#include "bam.h"
#include "punctures.h"

#define pow2(x)    ((x)*(x))


/* level to be used in m1_m2_VectorFunc */
tL *level_for_m1_m2_VectorFunc;

/* functions */
void m1_m2_VectorFunc(int n, double *vec, double *fvec);
void m1_VectorFunc(int n, double *vec, double *fvec);
int settoconstant_on_all_levels(tG *grid, int varindex, const double co);


/* Compute puncture data with a given M_ADM_1 and M_ADM_2 */
int PunctureData_with_given_M_ADM(tL *level)
{
  tG *g = level->grid;
  int lmax = g->lmax;
  int u1_eq_u2 = 0;
  int Knonzero;
  int punctures_persist = Getv("punctures_persist", "yes");
  double M_ADM_1 = Getd("puncture_M_ADM_1");
  double M_ADM_2 = Getd("puncture_M_ADM_2");
  double M_ADM_tol = Getd("puncture_M_ADM_tol");
  double diff1, diff2;
  double M1, M2;
  double m1, m2;
  double u1, u2;
  double v1sqr, v2sqr;
  double s1sqr, s2sqr;
  double D;
  double A1, A2, B, temp;

  /* wait until we are called on finest level */
  if (level->l < lmax) return 0;

  /* save level for m1_m2_VectorFunc */
  level_for_m1_m2_VectorFunc=level;
   
  /* take symmetries into account */
  if(Getv("grid", "quadrant")) u1_eq_u2 = 1;
  if(Getv("grid", "octant"))   u1_eq_u2 = 1;
  if(Getv("grid", "rotant"))   u1_eq_u2 = 1;
      
  /* set punctures_persist=yes to keep punctures_u in memory */
  if(!punctures_persist) Sets("punctures_persist", "yes");
  
  /* read parameters */
  // Setd("bhmass1", M_ADM_1);
  // Setd("bhmass2", M_ADM_2);
  ReadPunctureParameters(&Knonzero);
  v1sqr = (PBL[0][0]*PBL[0][0] + PBL[0][1]*PBL[0][1] + PBL[0][2]*PBL[0][2])/
          (M_ADM_1*M_ADM_1);
  v2sqr = (PBL[1][0]*PBL[1][0] + PBL[1][1]*PBL[1][1] + PBL[1][2]*PBL[1][2])/
          (M_ADM_2*M_ADM_2);
  s1sqr = SBL[0][0]*SBL[0][0] + SBL[0][1]*SBL[0][1] + SBL[0][2]*SBL[0][2];
  s2sqr = SBL[1][0]*SBL[1][0] + SBL[1][1]*SBL[1][1] + SBL[1][2]*SBL[1][2];

  /* estimate u1 and u2 */
  u1 = 0.01*v1sqr + 0.25*s1sqr/(M_ADM_1*M_ADM_1*M_ADM_1*M_ADM_1);
  u2 = 0.01*v2sqr + 0.25*s2sqr/(M_ADM_2*M_ADM_2*M_ADM_2*M_ADM_2);

  /* estimate for m1 and m2 */
  m1 = Getd("mass1");
  m2 = Getd("mass2");
    
  /* iterate until tolerance criterion is satisfied */
  if(Getv("puncture_M_ADM_newton_lnsrch", "yes")) /* use newton_lnsrch */
  {
    double m1_m2_vec[3];
    int check;

    /* initial guesses for m1 and m2 */
    m1_m2_vec[1]=m1;
    m1_m2_vec[2]=m2;

    if(u1_eq_u2)
    {
      /* do newton_lnsrch iterations if m1=m2 */
      printf("##############################################################################\n");
      printf("PunctureData_with_given_M_ADM:\n");
      printf("Starting with:\n");
      printf("m1=m2=m1_m2_vec[1]=%.15e\n", m1_m2_vec[1]);
  
      newton_lnsrch(m1_m2_vec, 1, &check, m1_VectorFunc, 
      		    Geti("puncture_M_ADM_MAXITS"),
      		    Getd("puncture_M_ADM_tol") );
  
      printf("##############################################################################\n");
      printf("PunctureData_with_given_M_ADM:\n");
      printf("Finished with:\n");
      printf("m1=m2=m1_m2_vec[1]=%.15e\n", m1_m2_vec[1]);
      printf("newton_lnsrch says: check=%d\n", check);
      printf("##############################################################################\n");
    }
    else
    {
      /* do newton_lnsrch iterations: */
      printf("##############################################################################\n");
      printf("PunctureData_with_given_M_ADM:\n");
      printf("Starting with:\n");
      printf("m1=m1_m2_vec[1]=%.15e  m2=m1_m2_vec[2]=%.15e\n", 
              m1_m2_vec[1], m1_m2_vec[2]);
  
      newton_lnsrch(m1_m2_vec, 2, &check, m1_m2_VectorFunc, 
      		    Geti("puncture_M_ADM_MAXITS"),
      		    Getd("puncture_M_ADM_tol") );
  
      printf("##############################################################################\n");
      printf("PunctureData_with_given_M_ADM:\n");
      printf("Finished with:\n");
      printf("m1=m1_m2_vec[1]=%.15e  m2=m1_m2_vec[2]=%.15e\n", 
              m1_m2_vec[1], m1_m2_vec[2]);
      printf("newton_lnsrch says: check=%d\n", check);
      printf("##############################################################################\n");
    }
  }
  else
  { 
     /* estimate m1 and m2 */
     m1 = M_ADM_1;
     m2 = M_ADM_2;

    /* use iterations of the form x_n+1 = f(x_n) */
    do
    {
      /* compute distance D */    
      D = sqrt( (CBL[0][0]-CBL[1][0])*(CBL[0][0]-CBL[1][0]) + 
                (CBL[0][1]-CBL[1][1])*(CBL[0][1]-CBL[1][1]) +
                (CBL[0][2]-CBL[1][2])*(CBL[0][2]-CBL[1][2])  );
  
      /* set bare masses */
      /*
      A1 = 1.0 + u1;
      A2 = 1.0 + u2;
      B  = 0.5/D;
      m1 = (-(A1*A2) + B*(M1 - M2) +
            sqrt(4.0*A1*A2*B*M2 + pow2((A1*A2 + B*(M1 - M2) )) ))/(2.*A1*B);
      m2 = (-(A1*A2) + B*(M2 - M1) +
            sqrt(4.0*A1*A2*B*M1 + pow2((A1*A2 + B*(M2 - M1) )) ))/(2.*A2*B);
      */
      temp = M_ADM_1/( 1.0 + u1 + m2/(2.0*D) );
      m2   = M_ADM_2/( 1.0 + u2 + m1/(2.0*D) );
      m1   = temp;
      Setd("mass1", m1);
      Setd("mass2", m2);
  
      /* call regular puncture solver */
      settoconstant_on_all_levels(g, Ind("punctures_u"), 0.0);
      settoconstant_on_all_levels(g, Ind("punctures_v"), 0.0);
      settoconstant_on_all_levels(g, Ind("punctures_f"), 0.0);
      settoconstant_on_all_levels(g, Ind("punctures_r"), 0.0);
      PunctureData(level);
  
      u1 = interpolate_xyz_scalar(g->level[g->lmax], 
                                  CBL[0][0], CBL[0][1], CBL[0][2], 
                                  Ind("punctures_u"), Geti("order_RP"),LAGRANGE);
      u2 = interpolate_xyz_scalar(g->level[g->lmax], 
                                  CBL[1][0], CBL[1][1], CBL[1][2],
                                  Ind("punctures_u"), Geti("order_RP"),LAGRANGE);
      if(u1_eq_u2) u2 = u1;
  
      /* compute ADM masses M1 and M2 */    
      M1 = ( 1.0 + u1 + m2/(2.0*D) )*m1;
      M2 = ( 1.0 + u2 + m1/(2.0*D) )*m2;
      printf("ADM masses at punctures:  M1=%.15g,  M2=%.15g\n", M1, M2);
      diff1 = fabs( M1/M_ADM_1 - 1.0 );
      diff2 = fabs( M2/M_ADM_2 - 1.0 );
    } while( diff1 > M_ADM_tol || diff2 > M_ADM_tol );
  }

  printf("ADM masses at punctures are now set to puncture_M_ADM_1 "
         "and puncture_M_ADM_2\n" 
         "(within tolerance puncture_M_ADM_tol=%.15g).\n", M_ADM_tol);

  /* reset punctures_persist to the value specified by the user */
  if(!punctures_persist) Sets("punctures_persist", "no");
    
  prdivider(0);
  return 0;
}


/* function called by newton_lnsrch if m1 != m2 */
void m1_m2_VectorFunc(int n, double *vec, double *fvec)
{
  tG *g = level_for_m1_m2_VectorFunc->grid;
  double m1=vec[1];
  double m2=vec[2];
  double M_ADM_1 = Getd("puncture_M_ADM_1");
  double M_ADM_2 = Getd("puncture_M_ADM_2");
  double M1, M2;
  int u1_eq_u2 = 0;
  double u1, u2;
  double D;

  /* take symmetries into account */
  if(Getv("grid", "quadrant")) u1_eq_u2 = 1;
  if(Getv("grid", "octant"))   u1_eq_u2 = 1;
  if(Getv("grid", "rotant"))   u1_eq_u2 = 1;

  /* set parameters */
  Setd("mass1", m1);
  Setd("mass2", m2);

  /* compute distance D */    
  D = sqrt( (CBL[0][0]-CBL[1][0])*(CBL[0][0]-CBL[1][0]) + 
            (CBL[0][1]-CBL[1][1])*(CBL[0][1]-CBL[1][1]) +
            (CBL[0][2]-CBL[1][2])*(CBL[0][2]-CBL[1][2])  );

  printf("##############################################################################\n");
  printf("m1_m2_VectorFunc:\n");
  printf("m1=vec[1]=%.15e  m2=vec[2]=%.15e\n", vec[1], vec[2]);
  fflush(stdout);
  
  /* call regular puncture solver */
  settoconstant_on_all_levels(g, Ind("punctures_u"), 0.0);
  settoconstant_on_all_levels(g, Ind("punctures_v"), 0.0);
  settoconstant_on_all_levels(g, Ind("punctures_f"), 0.0);
  settoconstant_on_all_levels(g, Ind("punctures_r"), 0.0);
  PunctureData(level_for_m1_m2_VectorFunc);

  u1 = interpolate_xyz_scalar(g->level[g->lmax], 
                              CBL[0][0], CBL[0][1], CBL[0][2], 
                              Ind("punctures_u"), Geti("order_RP"),LAGRANGE);
  u2 = interpolate_xyz_scalar(g->level[g->lmax], 
                              CBL[1][0], CBL[1][1], CBL[1][2], 
                              Ind("punctures_u"), Geti("order_RP"),LAGRANGE);
  if(u1_eq_u2) u2 = u1;

  /* compute ADM masses M1 and M2 */    
  M1 = ( 1.0 + u1 + m2/(2.0*D) )*m1;
  M2 = ( 1.0 + u2 + m1/(2.0*D) )*m2;
  printf("ADM masses at punctures:  M1=%.15g,  M2=%.15g\n", M1, M2);
                   	           
  fvec[1] = M_ADM_1 - M1;
  fvec[2] = M_ADM_2 - M2;

  printf("fvec[1]=M_ADM_1-M1=%.15e  fvec[2]=M_ADM_2-M2=%.15e\n", fvec[1], fvec[2]);
  fflush(stdout);
}


/* function called by newton_lnsrch if m1 = m2 */
void m1_VectorFunc(int n, double *vec, double *fvec)
{
  tG *g = level_for_m1_m2_VectorFunc->grid;
  double m1=vec[1];
  double M_ADM_1 = Getd("puncture_M_ADM_1");
  double M1;
  double u1;
  double D;

  /* assume symmetries, i.e. u1=u2, m1=m2 */
  Setd("mass1", m1);
  Setd("mass2", m1);

  /* compute distance D */    
  D = sqrt( (CBL[0][0]-CBL[1][0])*(CBL[0][0]-CBL[1][0]) + 
            (CBL[0][1]-CBL[1][1])*(CBL[0][1]-CBL[1][1]) +
            (CBL[0][2]-CBL[1][2])*(CBL[0][2]-CBL[1][2])  );

  printf("##############################################################################\n");
  printf("m1_VectorFunc:\n");
  printf("m1=m2=vec[1]=%.15e\n", vec[1]);
  fflush(stdout);
  
  /* call regular puncture solver */
  settoconstant_on_all_levels(g, Ind("punctures_u"), 0.0);
  settoconstant_on_all_levels(g, Ind("punctures_v"), 0.0);
  settoconstant_on_all_levels(g, Ind("punctures_f"), 0.0);
  settoconstant_on_all_levels(g, Ind("punctures_r"), 0.0);
  PunctureData(level_for_m1_m2_VectorFunc);

  u1 = interpolate_xyz_scalar(g->level[g->lmax],
                              CBL[0][0], CBL[0][1], CBL[0][2], 
                              Ind("punctures_u"), Geti("order_RP"),LAGRANGE);

  /* compute ADM mass M1 */    
  M1 = ( 1.0 + u1 + m1/(2.0*D) )*m1;
  printf("ADM masses at punctures:  M1=M2=%.15g\n", M1);
                   	           
  fvec[1] = M_ADM_1 - M1;

  printf("fvec[1]=M_ADM_1-M1=%.15e\n", fvec[1]);
  fflush(stdout);
}


/* set the var with index varindex to zero on all levels if it is allocated */
int settoconstant_on_all_levels(tG *grid, int varindex, const double co)
{
  int lmin = grid->lmin;
  int lmax = grid->lmax;
  int l,i;

  for (l = lmin; l <= lmax; l++)
  {
     tL *level_l = grid->level[l];
     double *var = level_l->v[varindex];
     if(var!=NULL)
       forallpoints(level_l, i) var[i] = co;
  }
  return 0;
}

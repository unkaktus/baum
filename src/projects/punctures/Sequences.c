/* Sequences.c */
/* Wolfgang Tichy 6/2003 */

#include "bam.h"
#include "punctures.h"



/* level to be used in P_c_VectorFunc */
tL *level_for_P_c_VectorFunc;

void P_c_VectorFunc(int n, double *vec, double *fvec);



/* Make a sequence of constant m_1/2 = bhmass1/2 by
   iterating momentum P and lapse_at_puncture c with newton_lnsrch 
   until MK = MADM     at infinity 
   and   MK = -c MADM  at the punctures */
void Make_MK_eq_MADM_Sequence(tL *level)
{
  tMasses masses;
  int N,i,j;
  double D;
  double P_c_vec[3];
  int check;
  double D1=Getd("punctures_Sequence_D1");
  double Dstep;
  double D2;
  int PointNum=Getd("punctures_Sequence_PointNumber");
  double rc1,rc2;
  double mu, M, nu, DoM, PsqrPN2; 
  FILE *out;
  char filename[10000];
  
  if(Getv("punctures_Sequence_D2", "use_Dstep"))
  {
     Dstep=Getd("punctures_Sequence_Dstep");
     D2=D1+Dstep*(PointNum-1.0);
  }
  else
  {
     D2=Getd("punctures_Sequence_D2");    
     Dstep=(D2-D1)/(PointNum-1.0);
  }

  /* save level for P_c_VectorFunc */ 
  level_for_P_c_VectorFunc=level;
  
  /* init files which will hold sequence data */
  if (processor0)
  {
      sprintf(filename, "%s/%s", Gets("outdir"), "punctures_Sequence.txt");
      out=fopen(filename,"a");
      if (out==NULL) errorexits("failed opening %s", filename);
      fprintf(out, 
      "bhr1    \tbhr2    \tbhp     \talpha_punc\tbhmass1 \tbhmass2 \t"
      "M_ADM_inf \tM_ADM_punc\tM_K      \tM_K_punc \tJ/(mu M_ADM_p)\t"
      "M_ADM_p*Omega\t"
      "D            \tM_ADM_1 \tM_ADM_2 \tM_K_1   \tM_K_2     \t"
      "Sx1     \tSy1     \tSz1     \tSx2     \tSy2     \tSz2\n");   
      fclose(out);
      
      sprintf(filename, "%s/%s", Gets("outdir"), "punctures_bareSequence.txt");
      out=fopen(filename,"a");
      if (out==NULL) errorexits("failed opening %s", filename);
      fprintf(out, 
      "bhr1    \tbhr2    \tbhp     \talpha_punc\tbhmass1 \tbhmass2 \t"
      "M_ADM_inf \tM_ADM_punc\tM_K      \tM_K_punc \tJ/(mu M_ADM_p)\t"
      "M_ADM_p*Omega\t"
      "D            \tM_ADM_1 \tM_ADM_2 \tM_K_1   \tM_K_2     \t"
      "Sx1     \tSy1     \tSz1     \tSx2     \tSy2     \tSz2\n");
      fclose(out);
  }

  for(N=1, D=D1;  N<=PointNum;   N++, D=D+Dstep)
  {
    /* masses */
    M = MBL[0] + MBL[1];
    mu = MBL[0] * MBL[1] / M;
    nu = mu/M;
    DoM = D/M;
    
    /* put center of mass into origin */
    rc1=+(MBL[1]/M)*D;
    rc2=-(MBL[0]/M)*D;
    
    /* put punctures on y-axis */ 
    CBL[0][0] = 0;
    CBL[0][1] = rc1;
    CBL[0][2] = 0;
    CBL[1][0] = 0;
    CBL[1][1] = rc2;
    CBL[1][2] = 0;

    /* initial guesses for c and P */
    PsqrPN2 = mu*mu * (74 - 43*nu + 32*DoM + 8*DoM*DoM) / (8 * DoM*DoM*DoM);
    P_c_vec[1]=sqrt(PsqrPN2);     /* use PN momentum */
    P_c_vec[2]=-1.0;              /* lapse at punc */

    /* do newton_lnsrch iterations: */
    printf("##############################################################################\n");
    printf("Sequences:\n");
    printf("D=%e  rc1=%e  rc2=%e\n", D,rc1,rc2);
    printf("P=P_c_vec[1]=%e  lapse_at_puncture=P_c_vec[2]=%e\n", 
            P_c_vec[1], P_c_vec[2]);

    newton_lnsrch(P_c_vec, 2, &check, P_c_VectorFunc, 
    		  Geti("punctures_Sequence_newtMAXITS"),
    		  Getd("punctures_Sequence_newtTOLF") );

    printf("##############################################################################\n");
    printf("Sequences:\n");
    printf("Finished for:  D=%e  rc1=%e  rc2=%e\n", D,rc1,rc2);
    printf("P=P_c_vec[1]=%e  lapse_at_puncture=P_c_vec[2]=%e\n", 
            P_c_vec[1], P_c_vec[2]);
    printf("newton_lnsrch says: check=%d\n", check);
    printf("##############################################################################\n");
        
    /* Save scaled data */
    PunctureMass(level,0, &masses);
    if (processor0)
    {
      double bhmass1, bhmass2, bhr1, bhr2, bhp, c;
      double M_ADM_inf, M_ADM_p, M_K_inf, M_K_p, J, Omega;
      double M_ADM_1, M_ADM_2, M_K_1, M_K_2;
      double Sx1,Sy1,Sz1, Sx2,Sy2,Sz2;

      M_ADM_p=(masses.M_ADM_1 + masses.M_ADM_2);
      bhmass1=MBL[0]/M_ADM_p;
      bhmass2=MBL[1]/M_ADM_p;
      bhr1=rc1/M_ADM_p;
      bhr2=rc2/M_ADM_p;
      bhp=P_c_vec[1]/M_ADM_p;
      c=P_c_vec[2];
      M_ADM_inf=masses.M_ADM_inf/M_ADM_p;
      M_K_inf=masses.M_K_inf/M_ADM_p;
      M_K_p=(masses.M_K_1 + masses.M_K_2)/M_ADM_p;
      J=masses.J_ADM/(masses.M_ADM_1 * masses.M_ADM_2);
      Omega=masses.Omega*M_ADM_p;
      M_ADM_1=masses.M_ADM_1/M_ADM_p;
      M_ADM_2=masses.M_ADM_2/M_ADM_p;
      M_K_1=masses.M_K_1/M_ADM_p;
      M_K_2=masses.M_K_2/M_ADM_p;
      Sx1=masses.S1[1]/(M_ADM_p*M_ADM_p);
      Sy1=masses.S1[2]/(M_ADM_p*M_ADM_p);
      Sz1=masses.S1[3]/(M_ADM_p*M_ADM_p);
      Sx2=masses.S2[1]/(M_ADM_p*M_ADM_p);
      Sy2=masses.S2[2]/(M_ADM_p*M_ADM_p);
      Sz2=masses.S2[3]/(M_ADM_p*M_ADM_p);      

      sprintf(filename, "%s/%s", Gets("outdir"), "punctures_Sequence.txt");
      out=fopen(filename,"a");
      if (out==NULL) errorexits("failed opening %s", filename);
      fprintf(out,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t"
                  "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t"
                  "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e \n",
              bhr1,bhr2, bhp, c, bhmass1,bhmass2, 
              M_ADM_inf, 1.0, M_K_inf, M_K_p, J, Omega, D/M_ADM_p,
              M_ADM_1, M_ADM_2, M_K_1, M_K_2, 
              Sx1,Sy1,Sz1, Sx2,Sy2,Sz2);
      fclose(out);
          
      sprintf(filename, "%s/%s", Gets("outdir"), "punctures_bareSequence.txt");
      out=fopen(filename,"a");
      if (out==NULL) errorexits("failed opening %s", filename);
      fprintf(out,"%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t"
                  "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t"
                  "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e \n",
              rc1,rc2, P_c_vec[1], P_c_vec[2], MBL[0],MBL[1], 
              masses.M_ADM_inf, M_ADM_p,
              masses.M_K_inf, masses.M_K_1 + masses.M_K_2, J, Omega, D,
              masses.M_ADM_1, masses.M_ADM_2, masses.M_K_1, masses.M_K_2,
              masses.S1[1], masses.S1[2], masses.S1[3],
              masses.S2[1], masses.S2[2], masses.S2[3]);
      fclose(out);
    }
  }

  if(Getv("punctures_shift", "thinsandwich"))
    PunctureShift(level);

  /* mass */
  if(!Getv("punctures_mass", "no")) 
    PunctureMass(level,1, &masses);
}


void P_c_VectorFunc(int n, double *vec, double *fvec)
{
  double P0=vec[1];
  double c0=vec[2];
  double MADMminusMK_inf;
  tMasses masses;

  printf("##############################################################################\n");
  printf("P_c_VectorFunc:\n");
  printf("P0=vec[1]=%e  c0=vec[2]=%e\n", vec[1], vec[2]);
  fflush(stdout);
  
  PBL[0][0] = P0;
  PBL[0][1] = 0;
  PBL[0][2] = 0;
  PBL[1][0] = -P0;     
  PBL[1][1] = 0; 
  PBL[1][2] = 0; 
            
  /* set Brill-Lindquist data */
  SetMWBL(level_for_P_c_VectorFunc);
        
  /* set Bowen-York extrinsic curvature */
  SetBYK(level_for_P_c_VectorFunc);

  /* precompute the coefficients */
  SetPunctureCoefficients(level_for_P_c_VectorFunc);

  /* call elliptic solvers */
  PunctureSolve(level_for_P_c_VectorFunc);

  Setd("punctures_lapse_at_puncture",c0);
  PunctureMaximalSlicing(level_for_P_c_VectorFunc);

  PunctureMass(level_for_P_c_VectorFunc,0, &masses);
  printf("##############################################################################\n");
  printf("P_c_VectorFunc:\n");
  printf("masses.M_ADM_inf = %e\nmasses.M_K_inf   = %e\n",
          masses.M_ADM_inf, masses.M_K_inf);
  printf("     masses.M_ADM_1 + masses.M_ADM_2 = %e\n"
         "     masses.M_K_1   + masses.M_K_2   = %e\n",
          masses.M_ADM_1 + masses.M_ADM_2, masses.M_K_1+masses.M_K_2);
  printf("-c0*(masses.M_ADM_1 + masses.M_ADM_2)= %e\n", 
         -c0 * (masses.M_ADM_1 + masses.M_ADM_2) );
                   	           
  fvec[1]=masses.M_ADM_inf - masses.M_K_inf;
  fvec[2]= -c0 * (masses.M_ADM_1 + masses.M_ADM_2) 
  	       - (masses.M_K_1   + masses.M_K_2  );
  /* fvec[2]= -c0 * masses.M_ADM_1 - masses.M_K_1; */
 
  printf("fvec[1]=%e  fvec[2]=%e\n", fvec[1], fvec[2]);
  fflush(stdout);
}

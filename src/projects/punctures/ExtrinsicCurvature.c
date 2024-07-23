/* ExtrinsicCurvature.c */
/* Bernd Bruegmann 1/03 */
/* based on brbr.c, Bernd Bruegmann, 11/97 */


#include "bam.h"
#include "punctures.h"


#define fori     for (i = 1; i <= 3; i++)
#define forsymij fori for (j = i; j <= 3; j++)
#define delta(i,j) (((i)==(j))?1.0:0.0)




/* Bowen-York W^i , such that K^{ij}=LW^{ij} */
void BY_Wofxyz(double x, double y, double z,
               double *W1, double *W2, double *W3)
{
  int n,i;
  double xc, yc, zc, r;
  double *W[4], P[4], S[4], N[4], NP, eNS[4];
  
  W[1] = W1;
  W[2] = W2;
  W[3] = W3;
  fori *W[i]=0.0;
      
  for (n = 0; MBL[n] != 0.0 && n < N_MWBL; n++)
  {
    xc = x-CBL[n][0];
    yc = y-CBL[n][1];
    zc = z-CBL[n][2];
    r  = sqrt(xc*xc+yc*yc+zc*zc);
    if (dequal(r, 0.0)) r = EPSZERO;

    N[1] = xc/r;
    N[2] = yc/r;
    N[3] = zc/r;
    
    fori P[i] = PBL[n][i-1];
    fori S[i] = SBL[n][i-1];
    
    NP = N[1]*P[1] + N[2]*P[2] + N[3]*P[3];

    eNS[1] = N[2]*S[3] - N[3]*S[2];
    eNS[2] = N[3]*S[1] - N[1]*S[3];
    eNS[3] = N[1]*S[2] - N[2]*S[1];
    
    fori *W[i] += -( 7.0*P[i] + N[i]*NP )/(4.0*r) + (eNS[i])/(r*r); 

    //fori *W[i] = pow(r, 2.0) * N[i];
  }
}
                              



/* Bowen-York extrinsic curvature */
void BY_Kofxyz(double x, double y, double z,
	      double *K11, double *K12, double *K13, 
	      double *K22, double *K23, double *K33)
{
  double xold = x, yold = y, zold = z;
  int i, j, n;
  double xc, yc, zc, r;
  double *K[4][4], P[4], S[4], N[4], NP, eNS[4];

  if (fisheye) { coordtrans(&x, &y, &z); }

  K[1][1] = K11;
  K[1][2] = K12;
  K[1][3] = K13;
  K[2][2] = K22;
  K[2][3] = K23;
  K[3][3] = K33;
  forsymij
    *K[i][j] = 0;

  for (n = 0; MBL[n] != 0.0 && n < N_MWBL; n++) {
    xc = x-CBL[n][0];
    yc = y-CBL[n][1];
    zc = z-CBL[n][2];
    r  = sqrt(xc*xc+yc*yc+zc*zc);
    if (dequal(r, 0.0)) r = EPSZERO;

    N[1] = xc/r;
    N[2] = yc/r;
    N[3] = zc/r;
    
    fori P[i] = PBL[n][i-1];
    fori S[i] = SBL[n][i-1];
    
    NP = N[1]*P[1] + N[2]*P[2] + N[3]*P[3];
    
    eNS[1] = N[2]*S[3] - N[3]*S[2];
    eNS[2] = N[3]*S[1] - N[1]*S[3];
    eNS[3] = N[1]*S[2] - N[2]*S[1];
    
    forsymij 
      *K[i][j] +=
      1.5/(r*r) * (N[i]*P[j] + N[j]*P[i] - (delta(i,j) - N[i]*N[j])*NP) -
      3.0/(r*r*r) * (eNS[i]*N[j] + eNS[j]*N[i]);

    // forsymij
    //  *K[i][j] = 2*(2-1)*pow(r,2.0-1.0) * (N[i]*N[j] - delta(i,j)/3);
  }

#if 0
  if (KK);
    *KK = K11[ijk]*K11[ijk] + K22[ijk]*K22[ijk] + K33[ijk]*K33[ijk] + 
          2*K12[ijk]*K12[ijk] + 2*K13[ijk]*K13[ijk] + 2*K23[ijk]*K23[ijk];
#endif

  if (fisheye) {
    indextrans_ddsym(xold, yold, zold, K11, K12, K13, K22, K23, K33);
  }
}




/* Bowen-York KK = K^{ij} K_{ij} from  BY_Kofxyz above*/
double BY_KKofxyz(double x, double y, double z)
{
  double KK;
  double K11, K12, K13, K22, K23, K33;

  BY_Kofxyz(x, y, z, &K11, &K12, &K13, &K22, &K23, &K33);

  KK = K11*K11 + K22*K22 + K33*K33 + 
          2.0*K12*K12 + 2.0*K13*K13 + 2.0*K23*K23;
  /* printf("KK=%e\n",KK); */
  return KK;
}




/* set Bowen-York extrinsic curvature */
void SetBYK(tL *level)
{
  int i;
  double *x = level->v[Ind("x")];
  double *y = level->v[Ind("y")];
  double *z = level->v[Ind("z")];
  double *K11, *K12, *K13, *K22, *K23, *K33;

  i = Ind("adm_Kxx");
  enablevar(level, i);
  K11 = level->v[i++];
  K12 = level->v[i++];
  K13 = level->v[i++];
  K22 = level->v[i++];
  K23 = level->v[i++];
  K33 = level->v[i++];

  forallpoints(level, i) {
    BY_Kofxyz(x[i], y[i], z[i], K11+i, K12+i, K13+i, K22+i, K23+i, K33+i);
  }
}





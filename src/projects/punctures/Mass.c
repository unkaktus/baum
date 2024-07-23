/* Mass.c */
/* Bernd Bruegmann 2/03 + WT 3/2003 */
/* compute different masses for puncture data based on conformal factor
*/

#include "bam.h"
#include "punctures.h"


/* ********************************************************************** */
extern double FullBox(double xmax, double ymax, double zmax,
                      double (*f)(double,double,double) );
extern double OutsideBox(double grid_xmax, double grid_ymax, double grid_zmax,
                         double xout, double (*f)(double,double,double) );


/* Funktion fn(x) to test outsidebox */
double fn(double x, double y, double z)
{
  return 1.0/( (x-5.0)*(x-5.0) * (y-4.0)*(y-4.0) * (z-3.0)*(z-3.0)); 
 /*  return x*(x-1) * y*(y-1) * z*(z-1); */
 /* return exp(1.0*x)*exp(2.0*y)*exp(3.0*z); */
}
/* ********************************************************************** */
                	  

/* compute coefficient of 1/r term in 
     u = b + c/r + O(1/r^2)
   using a volume integral over Laplace(u)
*/
double OneOverRCoefficient(tL *level, char *uname)
{
  tVarList *vllu  = VLPtrEnable1(level, "punctures_lu");  
  double *lu = Ptr(level, "punctures_lu");
  double *u = Ptr(level, uname);
  double cx = 1/(level->dx*level->dx);
  double cy = 1/(level->dy*level->dy);
  double cz = 1/(level->dz*level->dz);
  double cc = -2*(cx + cy + cz);
  double *result, sum, massinbox;
  int nsum;

  if (!u) return 0;

  /* compute Laplace(u) */
  forinner7(level) {
    lu[ccc] = cc * u[ccc] +
	      cx * (u[mcc] + u[pcc]) +
              cy * (u[cmc] + u[cpc]) +
	      cz * (u[ccm] + u[ccp]);
  } endfor;
  
  /* sum up contributions from all points and all processors */
  /* no need to synchronize because the sum does not include ghosts */
  bampi_allreduce_sum(vllu, &result);
  sum = result[0];
  nsum = result[1];
  free(result);
  VLDisableFree(vllu);

  /* integral over inner points with proper weight */
  massinbox = -1/(4*PI) * sum * (level->dx*level->dy*level->dz);

  /* correct for symmetry: assumes these boxes are exactly symmetric, i.e.
     with staggered points and no points on plane of symmetry 
  */
  if (Getv("grid", "bitant"))        massinbox *= 2;
  else if (Getv("grid", "quadrant")) massinbox *= 4;
  else if (Getv("grid", "octant"))   massinbox *= 8;

  if (0) printf("coefficient of 1/r for %11s:  %9.6f\n", uname, massinbox);
  return massinbox;
}




/* compute puncture masses */
double PunctureMass(tL *level, int fileoutput, tMasses *Masses)
{
  double f = -Getd("punctures_lapse_at_puncture");
  double adm, k_alpha, k_alphapsi;
  double c_psi, c_alphapsi;
  double c_u, c_alpha, c_v;
  double M;
  double xgrid = level->bbox[1] - level->dx/2.0;
  double ygrid = level->bbox[3] - level->dy/2.0;
  double zgrid = level->bbox[5] - level->dz/2.0;
  double M_ADM, M_Kalpha, M_Kalphapsi, C_alphapsi;
  double I_Lap_u_outside, I_Lap_v_outside, I_Lap_alpha_outside;
  double M_ADM0, M_ADM1, M_Kalphapsi0, M_Kalphapsi1;
  double u0,u1, v0,v1;
  double x12, y12, z12, r12;
  double Jx,Jy,Jz, Omega, J_ADM;
  FILE *out;
  char filename[10000];
  int n;

  /* compute sum of BL masses */
  M=0; 
  for(n = 0; MBL[n] != 0.0 && n < N_MWBL; n++)
     M += MBL[n];
  
  /* compute coefficient of 1/r term 
     and add in contributions from singular piece */
  if (Getv("punctures_mass", "ZeroCorrection"))
  { 
    printf("setting c_u = c_alpha = c_v = 0 \n");
    c_u     = 0.0;
    c_alpha = 0.0;
    c_v     = 0.0;
  }
  else
  { 
    printf("using punctures_u, alpha, punctures_v to estimate "
	   "c_u, c_alpha, c_v\n");
    c_u        = OneOverRCoefficient(level, "punctures_u");
    c_alpha    = OneOverRCoefficient(level, "alpha");
    c_v        = OneOverRCoefficient(level, "punctures_v");
  }
  c_psi      = c_u + M/2;
  c_alphapsi = c_v - f*M/2;

  adm        = 2 * c_psi;
  k_alpha    = - c_alpha;
  k_alphapsi = adm/2 - c_alphapsi;
  

  printf("Contribution to masses inside grid:\n");
  printf("PunctureMass:  M_ADM       = %9.6f\n", adm);
  printf("PunctureMass:  M_Kalpha    = %9.6f\n", k_alpha);
  printf("PunctureMass:  M_Kalphapsi = %9.6f = M_ADM/2 - C_alphapsi\n", 
	 k_alphapsi);
  printf("PunctureMass:  C_alphapsi  = %9.6f\n", c_alphapsi);

  
  /* run OutsideBox or FullBox to find correction to masses */
  printf("\nComputing corrections to masses:\n");
  SetFuncCoeffs(c_u, c_alpha, c_v);
  if (Getv("punctures_mass", "ZeroCorrection"))
  {
    I_Lap_u_outside = FullBox(10.0 * xgrid*xgrid, 10.0 * ygrid*ygrid, 
                              10.0 * zgrid*zgrid, RHS_of_Lapl_punctures_u);
    M_ADM = adm - I_Lap_u_outside/(2.0*PI);

    I_Lap_v_outside = FullBox(10.0 * xgrid*xgrid, 10.0 * ygrid*ygrid, 
                              10.0 * zgrid*zgrid, RHS_of_Lapl_punctures_v);
    C_alphapsi = c_alphapsi - I_Lap_v_outside/(4.0*PI);
    M_Kalphapsi = M/2.0 - C_alphapsi;

    I_Lap_alpha_outside = FullBox(10.0 * xgrid*xgrid, 10.0 * ygrid*ygrid,
                                  10.0 * zgrid*zgrid, RHS_of_Lapl_alpha);
    M_Kalpha = k_alpha + I_Lap_alpha_outside/(4.0*PI);
  }
  else
  { 
    I_Lap_u_outside = OutsideBox(xgrid,ygrid,zgrid, 10.0 * xgrid*xgrid, 
                                 RHS_of_Lapl_punctures_u);
    M_ADM = adm - I_Lap_u_outside/(2.0*PI);

    I_Lap_v_outside = OutsideBox(xgrid,ygrid,zgrid, 10.0 * xgrid*xgrid,
                                 RHS_of_Lapl_punctures_v);
    C_alphapsi = c_alphapsi - I_Lap_v_outside/(4.0*PI);
    M_Kalphapsi = M_ADM/2.0 - C_alphapsi;

    I_Lap_alpha_outside = OutsideBox(xgrid,ygrid,zgrid, 10.0 * xgrid*xgrid,
                                   RHS_of_Lapl_alpha);
    M_Kalpha = k_alpha + I_Lap_alpha_outside/(4.0*PI);
  }
  
  printf("\nMasses including corrections from OutsideBox or FullBox:\n");
  printf("PunctureMass:  M_ADM        = %9.6f  + %9.6f  = %9.6f\n", 
          adm, M_ADM-adm, M_ADM );
  printf("PunctureMass:  M_Kalpha     = %9.6f  + %9.6f  = %9.6f\n",
          k_alpha, M_Kalpha-k_alpha, M_Kalpha);
  printf("PunctureMass:  M_Kalphapsi  = %9.6f  + %9.6f  = %9.6f\n",
          k_alphapsi, M_Kalphapsi-k_alphapsi, M_Kalphapsi );
  printf("PunctureMass:  C_alphapsi   = %9.6f  + %9.6f  = %9.6f\n",
       	  c_alphapsi, C_alphapsi-c_alphapsi, C_alphapsi );

  /* distance between punctures */
  x12 = CBL[0][0] - CBL[1][0];
  y12 = CBL[0][1] - CBL[1][1];
  z12 = CBL[0][2] - CBL[1][2];
  r12  = sqrt( x12*x12 + y12*y12 + z12*z12);
  /* printf("\nDistance between punctures:  r12 =%9.6f\n", r12);  */

  printf("\npunctures_lapse_at_puncture=%f\n",
            Getd("punctures_lapse_at_puncture"));
		 	  
  /* interpolate u and v to puncture location */
  printf("Interpolating u and v to puncture location...\n");
  u0=interpolate_xyz_scalar(level, CBL[0][0], CBL[0][1], CBL[0][2],
                            Ind("punctures_u"), Geti("order_RP"),LAGRANGE);
  u1=interpolate_xyz_scalar(level, CBL[1][0], CBL[1][1], CBL[1][2],
                            Ind("punctures_u"), Geti("order_RP"),LAGRANGE);
                            
  v0=interpolate_xyz_scalar(level, CBL[0][0], CBL[0][1], CBL[0][2],
                            Ind("punctures_v"), Geti("order_RP"),LAGRANGE);                            
  v1=interpolate_xyz_scalar(level, CBL[1][0], CBL[1][1], CBL[1][2],
                            Ind("punctures_v"), Geti("order_RP"),LAGRANGE);
     
  /* ADM masses in other asympt. flat ends */
  if( r12 != 0.0 )
  {
    M_ADM0 = MBL[0] * (1+u0) + (MBL[0] * MBL[1])/(2.0*r12);
    M_ADM1 = MBL[1] * (1+u1) + (MBL[0] * MBL[1])/(2.0*r12);
  }
  else
  {
    M_ADM0 = MBL[0] * (1+u0);
    M_ADM1 = MBL[1] * (1+u1);
  }

  /* Komar in other asympt. flat ends */
  M_Kalphapsi0 = MBL[0] * ( 1+f + v0+f*u0 )/2.0;
  M_Kalphapsi1 = MBL[1] * ( 1+f + v1+f*u1 )/2.0;

  printf("\nMasses at other asymptotically flat end:\n");
  printf("PunctureMass:  M_ADM0 + M_ADM1             =%9.6f +%9.6f =%9.6f\n",
  	  M_ADM0, M_ADM1, M_ADM0 + M_ADM1 );
  printf("PunctureMass:  M_Kalphapsi0 + M_Kalphapsi1 =%9.6f +%9.6f =%9.6f\n",
  	  M_Kalphapsi0, M_Kalphapsi1, M_Kalphapsi0 + M_Kalphapsi1 );

  /* angular vel. Omega and ang. mom. J_ADM */
  Jx=0;  Jy=0;  Jz=0;
  for(n = 0; MBL[n] != 0.0 && n < N_MWBL; n++)
  {
    Jx += CBL[n][1]*PBL[n][2] - CBL[n][2]*PBL[n][1];
    Jy += CBL[n][2]*PBL[n][0] - CBL[n][0]*PBL[n][2];
    Jz += CBL[n][0]*PBL[n][1] - CBL[n][1]*PBL[n][0];
  }      
  J_ADM = sqrt( Jx*Jx + Jy*Jy +Jz*Jz );
  if( J_ADM != 0.0 )
     Omega = ( M_ADM - (M_Kalphapsi0 + M_Kalphapsi1) )/(2.0*J_ADM);
  else
     Omega = 0.0;

  printf("\nAngular Momentum and velocity at infinity:\n");
  printf("PunctureMass:  J_ADM = %9.6f\n", J_ADM);
  printf("PunctureMass:  Omega = %9.6f\n", Omega);
  
  /* write masses into structure Masses */
  Masses->M_ADM_inf=M_ADM;
  Masses->M_ADM_1=M_ADM0;
  Masses->M_ADM_2=M_ADM1;
  Masses->M_K_inf=M_Kalphapsi;
  Masses->M_K_1=M_Kalphapsi0;  
  Masses->M_K_2=M_Kalphapsi1;
  Masses->Omega=Omega;
  Masses->J_ADM=J_ADM;
  for(n=1; n<=3; n++)
  { 
    Masses->S1[n]=SBL[0][n-1];
    Masses->S2[n]=SBL[1][n-1];
  }
  Masses->S1[0]=sqrt( SBL[0][0]*SBL[0][0] + SBL[0][1]*SBL[0][1]
                     +SBL[0][2]*SBL[0][2] );
  Masses->S2[0]=sqrt( SBL[1][0]*SBL[1][0] + SBL[1][1]*SBL[1][1]
                     +SBL[1][2]*SBL[1][2] );
                       
  /* File output: */
  if (processor0 & fileoutput)
  {
    /* write M_ADM, M_Kalphapsi, M_Kalpha, ... in files */
    sprintf(filename, "%s/%s.%d", Gets("outdir"), "punctures_M_ADM", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, M_ADM);
    fclose(out);
    
    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_M_Kalphapsi", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, M_Kalphapsi);
    fclose(out);
    
    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_M_Kalpha", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, M_Kalpha);
    fclose(out);

    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_M_ADM1", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, M_ADM0);
    fclose(out);

    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_M_ADM2", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, M_ADM1);
    fclose(out);

    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_M_Kalphapsi1", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, M_Kalphapsi0);
    fclose(out);

    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_M_Kalphapsi2", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, M_Kalphapsi1);
    fclose(out);

    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_J_ADM", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, J_ADM);
    fclose(out);

    sprintf(filename, "%s/%s.%d", 
            Gets("outdir"), "punctures_Omega", level->l);
    out=fopen(filename,"a");
    if (out==NULL) errorexits("failed opening %s", filename);
    fprintf(out, "\"Time = %f\"\n", level->time);
    fprintf(out,"%22.15f  %22.15f\n", xgrid, Omega);
    fclose(out);

  }
  
  return (M_ADM - M_Kalphapsi)/M_ADM;
}


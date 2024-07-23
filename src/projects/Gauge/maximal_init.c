/* maximal_init.c */
/* Copyright (C) 1998 Bernd Bruegmann, 12.3.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "Gauge.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define pow4(x)    ((x)*(x)*(x)*(x))
#define Cal(x,y,z) ((x)?(y):(z))




void maximal_init(tL *level)
{

int matter = Getv("physics","matter");

int index_gxx = Ind("gxx");
double *g11 = level->v[index_gxx + 0];
double *g12 = level->v[index_gxx + 1];
double *g13 = level->v[index_gxx + 2];
double *g22 = level->v[index_gxx + 3];
double *g23 = level->v[index_gxx + 4];
double *g33 = level->v[index_gxx + 5];
int index_Kxx = Ind("Kxx");
double *K11 = level->v[index_Kxx + 0];
double *K12 = level->v[index_Kxx + 1];
double *K13 = level->v[index_Kxx + 2];
double *K22 = level->v[index_Kxx + 3];
double *K23 = level->v[index_Kxx + 4];
double *K33 = level->v[index_Kxx + 5];
int index_maximal_gixx = Ind("maximal_gixx");
double *ginv11 = level->v[index_maximal_gixx + 0];
double *ginv12 = level->v[index_maximal_gixx + 1];
double *ginv13 = level->v[index_maximal_gixx + 2];
double *ginv22 = level->v[index_maximal_gixx + 3];
double *ginv23 = level->v[index_maximal_gixx + 4];
double *ginv33 = level->v[index_maximal_gixx + 5];
int index_maximal_Gx = Ind("maximal_Gx");
double *G1 = level->v[index_maximal_Gx + 0];
double *G2 = level->v[index_maximal_Gx + 1];
double *G3 = level->v[index_maximal_Gx + 2];
int index_maximal_KK = Ind("maximal_KK");
double *KK = level->v[index_maximal_KK + 0];
int index_psi = Ind("psi");
double *psi = level->v[index_psi + 0];
int index_dpsiopsix = Ind("dpsiopsix");
double *dpop1 = level->v[index_dpsiopsix + 0];
double *dpop2 = level->v[index_dpsiopsix + 1];
double *dpop3 = level->v[index_dpsiopsix + 2];
int index_rho,index_SSxx;
if (matter) {
index_rho = Ind("rho");
index_SSxx = Ind("SSxx");
} else {
printf(" these are dummy variables, because c has some problems by initialicing variables inside an if statement\n");index_rho = Ind("gxx");
index_rho = Ind("gxx");
}
double *rho = level->v[index_rho + 0];
double *SS11 = level->v[index_SSxx + 0];
double *SS12 = level->v[index_SSxx + 1];
double *SS13 = level->v[index_SSxx + 2];
double *SS22 = level->v[index_SSxx + 3];
double *SS23 = level->v[index_SSxx + 4];
double *SS33 = level->v[index_SSxx + 5];


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double delf1 = 0.;
double delf2 = 0.;
double delf3 = 0.;
double delg111 = 0.;
double delg112 = 0.;
double delg113 = 0.;
double delg122 = 0.;
double delg123 = 0.;
double delg133 = 0.;
double delg211 = 0.;
double delg212 = 0.;
double delg213 = 0.;
double delg222 = 0.;
double delg223 = 0.;
double delg233 = 0.;
double delg311 = 0.;
double delg312 = 0.;
double delg313 = 0.;
double delg322 = 0.;
double delg323 = 0.;
double delg333 = 0.;
double detginvf = 0.;
double f = 0.;
double gamma111 = 0.;
double gamma112 = 0.;
double gamma113 = 0.;
double gamma122 = 0.;
double gamma123 = 0.;
double gamma133 = 0.;
double gamma211 = 0.;
double gamma212 = 0.;
double gamma213 = 0.;
double gamma222 = 0.;
double gamma223 = 0.;
double gamma233 = 0.;
double gamma311 = 0.;
double gamma312 = 0.;
double gamma313 = 0.;
double gamma322 = 0.;
double gamma323 = 0.;
double gamma333 = 0.;
double gammado111 = 0.;
double gammado112 = 0.;
double gammado113 = 0.;
double gammado122 = 0.;
double gammado123 = 0.;
double gammado133 = 0.;
double gammado211 = 0.;
double gammado212 = 0.;
double gammado213 = 0.;
double gammado222 = 0.;
double gammado223 = 0.;
double gammado233 = 0.;
double gammado311 = 0.;
double gammado312 = 0.;
double gammado313 = 0.;
double gammado322 = 0.;
double gammado323 = 0.;
double gammado333 = 0.;



forinner19(level) {


delg111
=
oo2dx*(-g11[-di + ijk] + g11[di + ijk])
;

delg112
=
oo2dx*(-g12[-di + ijk] + g12[di + ijk])
;

delg113
=
oo2dx*(-g13[-di + ijk] + g13[di + ijk])
;

delg122
=
oo2dx*(-g22[-di + ijk] + g22[di + ijk])
;

delg123
=
oo2dx*(-g23[-di + ijk] + g23[di + ijk])
;

delg133
=
oo2dx*(-g33[-di + ijk] + g33[di + ijk])
;

delg211
=
oo2dy*(-g11[-dj + ijk] + g11[dj + ijk])
;

delg212
=
oo2dy*(-g12[-dj + ijk] + g12[dj + ijk])
;

delg213
=
oo2dy*(-g13[-dj + ijk] + g13[dj + ijk])
;

delg222
=
oo2dy*(-g22[-dj + ijk] + g22[dj + ijk])
;

delg223
=
oo2dy*(-g23[-dj + ijk] + g23[dj + ijk])
;

delg233
=
oo2dy*(-g33[-dj + ijk] + g33[dj + ijk])
;

delg311
=
oo2dz*(-g11[-dk + ijk] + g11[dk + ijk])
;

delg312
=
oo2dz*(-g12[-dk + ijk] + g12[dk + ijk])
;

delg313
=
oo2dz*(-g13[-dk + ijk] + g13[dk + ijk])
;

delg322
=
oo2dz*(-g22[-dk + ijk] + g22[dk + ijk])
;

delg323
=
oo2dz*(-g23[-dk + ijk] + g23[dk + ijk])
;

delg333
=
oo2dz*(-g33[-dk + ijk] + g33[dk + ijk])
;

f
=
pow4(psi[ijk])
;

delf1
=
4.*f*dpop1[ijk]
;

delf2
=
4.*f*dpop2[ijk]
;

delf3
=
4.*f*dpop3[ijk]
;

delg111
=
delg111*f + delf1*g11[ijk]
;

delg112
=
delg112*f + delf1*g12[ijk]
;

delg113
=
delg113*f + delf1*g13[ijk]
;

delg122
=
delg122*f + delf1*g22[ijk]
;

delg123
=
delg123*f + delf1*g23[ijk]
;

delg133
=
delg133*f + delf1*g33[ijk]
;

delg211
=
delg211*f + delf2*g11[ijk]
;

delg212
=
delg212*f + delf2*g12[ijk]
;

delg213
=
delg213*f + delf2*g13[ijk]
;

delg222
=
delg222*f + delf2*g22[ijk]
;

delg223
=
delg223*f + delf2*g23[ijk]
;

delg233
=
delg233*f + delf2*g33[ijk]
;

delg311
=
delg311*f + delf3*g11[ijk]
;

delg312
=
delg312*f + delf3*g12[ijk]
;

delg313
=
delg313*f + delf3*g13[ijk]
;

delg322
=
delg322*f + delf3*g22[ijk]
;

delg323
=
delg323*f + delf3*g23[ijk]
;

delg333
=
delg333*f + delf3*g33[ijk]
;

detginvf
=
1/(f*(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
      g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
      g11[ijk]*pow2(g23[ijk])))
;

ginv11[ijk]
=
detginvf*(g22[ijk]*g33[ijk] - pow2(g23[ijk]))
;

ginv12[ijk]
=
detginvf*(g13[ijk]*g23[ijk] - g12[ijk]*g33[ijk])
;

ginv13[ijk]
=
detginvf*(-(g13[ijk]*g22[ijk]) + g12[ijk]*g23[ijk])
;

ginv22[ijk]
=
detginvf*(g11[ijk]*g33[ijk] - pow2(g13[ijk]))
;

ginv23[ijk]
=
detginvf*(g12[ijk]*g13[ijk] - g11[ijk]*g23[ijk])
;

ginv33[ijk]
=
detginvf*(g11[ijk]*g22[ijk] - pow2(g12[ijk]))
;

gammado111
=
0.5*delg111
;

gammado112
=
0.5*delg211
;

gammado113
=
0.5*delg311
;

gammado122
=
-0.5*delg122 + delg212
;

gammado123
=
0.5*(-delg123 + delg213 + delg312)
;

gammado133
=
-0.5*delg133 + delg313
;

gammado211
=
delg112 - 0.5*delg211
;

gammado212
=
0.5*delg122
;

gammado213
=
0.5*(delg123 - delg213 + delg312)
;

gammado222
=
0.5*delg222
;

gammado223
=
0.5*delg322
;

gammado233
=
-0.5*delg233 + delg323
;

gammado311
=
delg113 - 0.5*delg311
;

gammado312
=
0.5*(delg123 + delg213 - delg312)
;

gammado313
=
0.5*delg133
;

gammado322
=
delg223 - 0.5*delg322
;

gammado323
=
0.5*delg233
;

gammado333
=
0.5*delg333
;

gamma111
=
gammado111*ginv11[ijk] + gammado211*ginv12[ijk] + gammado311*ginv13[ijk]
;

gamma112
=
gammado112*ginv11[ijk] + gammado212*ginv12[ijk] + gammado312*ginv13[ijk]
;

gamma113
=
gammado113*ginv11[ijk] + gammado213*ginv12[ijk] + gammado313*ginv13[ijk]
;

gamma122
=
gammado122*ginv11[ijk] + gammado222*ginv12[ijk] + gammado322*ginv13[ijk]
;

gamma123
=
gammado123*ginv11[ijk] + gammado223*ginv12[ijk] + gammado323*ginv13[ijk]
;

gamma133
=
gammado133*ginv11[ijk] + gammado233*ginv12[ijk] + gammado333*ginv13[ijk]
;

gamma211
=
gammado111*ginv12[ijk] + gammado211*ginv22[ijk] + gammado311*ginv23[ijk]
;

gamma212
=
gammado112*ginv12[ijk] + gammado212*ginv22[ijk] + gammado312*ginv23[ijk]
;

gamma213
=
gammado113*ginv12[ijk] + gammado213*ginv22[ijk] + gammado313*ginv23[ijk]
;

gamma222
=
gammado122*ginv12[ijk] + gammado222*ginv22[ijk] + gammado322*ginv23[ijk]
;

gamma223
=
gammado123*ginv12[ijk] + gammado223*ginv22[ijk] + gammado323*ginv23[ijk]
;

gamma233
=
gammado133*ginv12[ijk] + gammado233*ginv22[ijk] + gammado333*ginv23[ijk]
;

gamma311
=
gammado111*ginv13[ijk] + gammado211*ginv23[ijk] + gammado311*ginv33[ijk]
;

gamma312
=
gammado112*ginv13[ijk] + gammado212*ginv23[ijk] + gammado312*ginv33[ijk]
;

gamma313
=
gammado113*ginv13[ijk] + gammado213*ginv23[ijk] + gammado313*ginv33[ijk]
;

gamma322
=
gammado122*ginv13[ijk] + gammado222*ginv23[ijk] + gammado322*ginv33[ijk]
;

gamma323
=
gammado123*ginv13[ijk] + gammado223*ginv23[ijk] + gammado323*ginv33[ijk]
;

gamma333
=
gammado133*ginv13[ijk] + gammado233*ginv23[ijk] + gammado333*ginv33[ijk]
;

G1[ijk]
=
gamma111*ginv11[ijk] + gamma122*ginv22[ijk] + 
  2.*(gamma112*ginv12[ijk] + gamma113*ginv13[ijk] + gamma123*ginv23[ijk]) + 
  gamma133*ginv33[ijk]
;

G2[ijk]
=
gamma211*ginv11[ijk] + gamma222*ginv22[ijk] + 
  2.*(gamma212*ginv12[ijk] + gamma213*ginv13[ijk] + gamma223*ginv23[ijk]) + 
  gamma233*ginv33[ijk]
;

G3[ijk]
=
gamma311*ginv11[ijk] + gamma322*ginv22[ijk] + 
  2.*(gamma312*ginv12[ijk] + gamma313*ginv13[ijk] + gamma323*ginv23[ijk]) + 
  gamma333*ginv33[ijk]
;



/* conditional */
if (matter) {

KK[ijk]
=
8.*PI*(ginv12[ijk]*SS12[ijk] + ginv13[ijk]*SS13[ijk] + 
     ginv23[ijk]*SS23[ijk]) + 4.*
   (ginv23[ijk]*(ginv12[ijk]*K12[ijk] + ginv13[ijk]*K13[ijk])*K23[ijk] + 
     ginv12[ijk]*(ginv13[ijk]*K11[ijk]*K23[ijk] + 
        K13[ijk]*(ginv23[ijk]*K22[ijk] + ginv33[ijk]*K23[ijk])) + 
     K12[ijk]*((ginv12[ijk]*ginv13[ijk] + ginv11[ijk]*ginv23[ijk])*
         K13[ijk] + ginv22[ijk]*
         (ginv12[ijk]*K22[ijk] + ginv13[ijk]*K23[ijk]) + 
        ginv13[ijk]*ginv23[ijk]*K33[ijk]) + 
     PI*(rho[ijk] + ginv33[ijk]*SS33[ijk])) + 
  K33[ijk]*(4.*ginv33[ijk]*(ginv13[ijk]*K13[ijk] + ginv23[ijk]*K23[ijk]) + 
     2.*K22[ijk]*pow2(ginv23[ijk])) + pow2(ginv11[ijk])*pow2(K11[ijk]) + 
  ginv11[ijk]*(4.*(K11[ijk]*(ginv12[ijk]*K12[ijk] + 
           ginv13[ijk]*K13[ijk]) + PI*SS11[ijk]) + 
     2.*(ginv22[ijk]*pow2(K12[ijk]) + ginv33[ijk]*pow2(K13[ijk]))) + 
  pow2(ginv22[ijk])*pow2(K22[ijk]) + 
  ginv22[ijk]*(4.*(ginv23[ijk]*K22[ijk]*K23[ijk] + PI*SS22[ijk]) + 
     2.*ginv33[ijk]*pow2(K23[ijk])) + 
  2.*(pow2(ginv12[ijk])*(K11[ijk]*K22[ijk] + pow2(K12[ijk])) + 
     pow2(ginv13[ijk])*(K11[ijk]*K33[ijk] + pow2(K13[ijk])) + 
     pow2(ginv23[ijk])*pow2(K23[ijk])) + pow2(ginv33[ijk])*pow2(K33[ijk])
;


} else { /* if (!matter) */

KK[ijk]
=
4.*(ginv23[ijk]*(ginv12[ijk]*K12[ijk] + ginv13[ijk]*K13[ijk])*K23[ijk] + 
     ginv12[ijk]*(ginv13[ijk]*K11[ijk]*K23[ijk] + 
        K13[ijk]*(ginv23[ijk]*K22[ijk] + ginv33[ijk]*K23[ijk])) + 
     K12[ijk]*((ginv12[ijk]*ginv13[ijk] + ginv11[ijk]*ginv23[ijk])*
         K13[ijk] + ginv22[ijk]*
         (ginv12[ijk]*K22[ijk] + ginv13[ijk]*K23[ijk]) + 
        ginv13[ijk]*ginv23[ijk]*K33[ijk])) + 
  K33[ijk]*(4.*ginv33[ijk]*(ginv13[ijk]*K13[ijk] + ginv23[ijk]*K23[ijk]) + 
     2.*K22[ijk]*pow2(ginv23[ijk])) + pow2(ginv11[ijk])*pow2(K11[ijk]) + 
  ginv11[ijk]*(4.*K11[ijk]*(ginv12[ijk]*K12[ijk] + ginv13[ijk]*K13[ijk]) + 
     2.*(ginv22[ijk]*pow2(K12[ijk]) + ginv33[ijk]*pow2(K13[ijk]))) + 
  pow2(ginv22[ijk])*pow2(K22[ijk]) + 
  ginv22[ijk]*(4.*ginv23[ijk]*K22[ijk]*K23[ijk] + 
     2.*ginv33[ijk]*pow2(K23[ijk])) + 
  2.*(pow2(ginv12[ijk])*(K11[ijk]*K22[ijk] + pow2(K12[ijk])) + 
     pow2(ginv13[ijk])*(K11[ijk]*K33[ijk] + pow2(K13[ijk])) + 
     pow2(ginv23[ijk])*pow2(K23[ijk])) + pow2(ginv33[ijk])*pow2(K33[ijk])
;

}
/* if (matter) */


ginv11[ijk]
=
f*ginv11[ijk]
;

ginv12[ijk]
=
f*ginv12[ijk]
;

ginv13[ijk]
=
f*ginv13[ijk]
;

ginv22[ijk]
=
f*ginv22[ijk]
;

ginv23[ijk]
=
f*ginv23[ijk]
;

ginv33[ijk]
=
f*ginv33[ijk]
;

G1[ijk]
=
f*G1[ijk]
;

G2[ijk]
=
f*G2[ijk]
;

G3[ijk]
=
f*G3[ijk]
;

KK[ijk]
=
f*KK[ijk]
;



} endfor;  /* loop i, j, k */

}  /* function */

/* maximal_init.c */
/* nvars = 33, nauxs = 59, n* = 363,  n/ = 27,  n+ = 314, n = 704, O = 0 */

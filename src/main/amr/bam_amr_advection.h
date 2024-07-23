/* bam_FiniteDiff_advection.h */
/* Bernd Bruegmann, 10/2002, 2/2006 */

/* Various macros to handle advection stencils
   This originally ended up in amr because the amr loops are supposed
   to hide the index constructions from the user.

   For example, w^x del_x u 
   becomes      w^x_i (3 u_i - 4 u_(i-1) + u_(i-2)) / (2h)
   if w^x is negative (upwinding from the left).

   Note that for us del_t u = + w^x del_x u which advects to the left for
   w^x positive, i.e. the shift vector has a sign opposite to velocity.
*/


/* new box version of advection at various orders */


/* N-th order: 
    advectionlopsidedN = 0 gives stencil oooooox  (one-sided) 
    advectionlopsidedN = 1 gives stencil oooooxo  (2lop-sided)
    advectionlopsidedN = 2 gives stencil ooooxoo  (lop-sided)
    advectionlopsidedN = 3 gives stencil oooxooo  (centered)
*/

/* define the stencils */
#define int_advectionN_stencil					              \
    const double stencil2[3][3] = {					      \
{+1.5, -2.0, +0.5}, {+0.5, 0.0, -0.5}, {-0.5, +2.0, -1.5}};		      \
    const double stencil4[5][5] = {				              \
{+25./12, -48./12, +36./12, -16./12, + 3./12},		  	      \
{+ 3./12, +10./12, -18./12, + 6./12, - 1./12},			      \
{- 1./12, + 8./12,   0./12, - 8./12, + 1./12},			      \
{+ 1./12, - 6./12, +18./12, -10./12, - 3./12},			      \
{- 3./12, +16./12, -36./12, +48./12, -25./12}};			      \
    const double stencil6[7][7] = {					      \
{+147./60, -360./60, +450./60, -400./60, +225./60, - 72./60, + 10./60},   \
{+ 10./60, + 77./60, -150./60, +100./60, - 50./60, + 15./60, -  2./60},   \
{-  2./60, + 24./60, + 35./60, - 80./60, + 30./60, -  8./60, +  1./60},   \
{+  1./60, -  9./60, + 45./60,    0./60, - 45./60, +  9./60, -  1./60},   \
{-  1./60, +  8./60, - 30./60, + 80./60, - 35./60, - 24./60, +  2./60},   \
{+  2./60, - 15./60, + 50./60, -100./60, +150./60, - 77./60, - 10./60},   \
{- 10./60, + 72./60, -225./60, +400./60, -450./60, +360./60, -147./60}};  \
    const double stencil8[9][9] = {					      \
{761./280, -8, 14, -56./3, 35./2, -56./5, 14./3, -8./7, 1./8},            \
{1./8, 223./140, -7./2, 7./2, -35./12, 7./4, -7./10, 1./6, -1./56},       \
{-1./56, 2./7, 19./20, -2, 5./4, -2./3, 1./4, -2./35, 1./168},            \
{1./168, -1./14, 1./2, 9./20, -5./4, 1./2, -1./6, 1./28, -1./280},        \
{-1./280, 4./105, -1./5, 4./5, 0, -4./5, 1./5, -4./105, 1./280},          \
{1./280, -1./28, 1./6, -1./2, 5./4, -9./20, -1./2, 1./14, -1./168},       \
{-1./168, 2./35, -1./4, 2./3, -5./4, 2, -19./20, -2./7, 1./56},           \
{1./56, -1./6, 7./10, -7./4, 35./12, -7./2, 7./2, -223./140, -1./8},      \
{-1./8, 8./7, -14./3, 56./5, -35./2, 56./3, -14, 8, -761./280}};          \
    const double stencil10[11][11] = {					      \
{7381./2520, -10, 45./2, -40, 105./2, -252./5, 35, -120./7, 45./8, -10./9, 1./10},      \
{1./10, 4609./2520, -9./2, 6, -7, 63./10, -21./5, 2, -9./14, 1./8, -1./90},             \
{-1./90, 2./9, 341./280, -8./3, 7./3, -28./15, 7./6, -8./15, 1./6, -2./63, 1./360},     \
{1./360, -1./24, 3./8, 319./420, -7./4, 21./20, -7./12, 1./4, -3./40, 1./72, -1./840},  \
{-1./840, 1./63, -3./28, 4./7, 11./30, -6./5, 1./2, -4./21, 3./56, -1./105, 1./1260},   \
{1./1260, -5./504, 5./84, -5./21, 5./6, 0, -5./6, 5./21, -5./84, 5./504, -1./1260},     \
{-1./1260, 1./105, -3./56, 4./21, -1./2, 6./5, -11./30, -4./7, 3./28, -1./63, 1./840},  \
{1./840, -1./72, 3./40, -1./4, 7./12, -21./20, 7./4, -319./420, -3./8, 1./24, -1./360}, \
{-1./360, 2./63, -1./6, 8./15, -7./6, 28./15, -7./3, 8./3, -341./280, -2./9, 1./90},    \
{1./90, -1./8, 9./14, -2, 21./5, -63./10, 7, -6, 9./2, -4609./2520, -1./10},            \
{-1./10, 10./9, -45./8, 120./7, -35, 252./5, -105./2, 40, -45./2, 10, -7381./2520}}

//due to OpenMP define used variables later
#define int_advectionN_vars	                                              \
    double wsign;							              \
    int offset, rstencil;							      \
    double advu0, advu1, advu2, advu3, advu4, advu5, advu6, advu7, advu8, advu9, advu10; \
    double advv0, advv1, advv2, advv3, advv4, advv5, advv6, advv7, advv8, advv9, advv10; \
    double advw0, advw1, advw2, advw3, advw4, advw5, advw6, advw7, advw8, advw9, advw10; \
    int advi0, advi1, advi2, advi3, advi4, advi5, advi6, advi7, advi8, advi9, advi10; \
    int advj0, advj1, advj2, advj3, advj4, advj5, advj6, advj7, advj8, advj9, advj10; \
    int advk0, advk1, advk2, advk3, advk4, advk5, advk6, advk7, advk8, advk9, advk10





/* figure out the direction of advection 
   again, wx < 0 pointing left means advection to the right using points
   to the left ...

   if we are close to a boundary, we may not be able to use stricly one-sided
   stencils; the strategy is to offset the stencil automatically!
*/
#define set_advection_wsign(i,di,imax,wx)				\
  if (wx < 0) {								\
    wsign = -1;								\
    offset = (i >= rstencil) ? 0 : rstencil - i;			\
  } else if (wx > 0) {							\
    wsign = +1;								\
    offset = (imax - i >= rstencil) ? 0 : rstencil - (imax - i);	\
  } else {								\
    wsign = 0;								\
    offset = 0;								\
  }									\




/* second order: oox */
#define set_advection2_onedir(i,di,imax,wx,				\
			      advi0,advi1,advi2,			\
			      advu0,advu1,advu2)			\
  set_advection_wsign(i,di,imax,wx);					\
  advi0 = ijk + wsign * (0 - offset) * di;				\
  advi1 = ijk + wsign * (1 - offset) * di;				\
  advi2 = ijk + wsign * (2 - offset) * di;				\
  advu0 = - wsign * wx * stencil2[offset][0];				\
  advu1 = - wsign * wx * stencil2[offset][1];				\
  advu2 = - wsign * wx * stencil2[offset][2];				\

#define set_advection2(wx,wy,wz)					\
  rstencil = 2;								\
  set_advection2_onedir(i,di,imax,wx,advi0,advi1,advi2,advu0,advu1,advu2); \
  set_advection2_onedir(j,dj,jmax,wy,advj0,advj1,advj2,advv0,advv1,advv2); \
  set_advection2_onedir(k,dk,kmax,wz,advk0,advk1,advk2,advw0,advw1,advw2); \





/* fourth order: 
   advectionlopsided = 0 gives stencil oooox  (one-sided) 
   advectionlopsided = 1 gives stencil oooxo  (lop-sided)
   advectionlopsided = 2 gives stencil ooxoo  (centered)
*/
#define set_advection4_onedir(i,di,imax,wx,				\
			      advi0,advi1,advi2,advi3,advi4,		\
			      advu0,advu1,advu2,advu3,advu4)		\
  set_advection_wsign(i,di,imax,wx);					\
  offset += advectionlopsided;						\
  advi0 = ijk + wsign * (0 - offset) * di;				\
  advi1 = ijk + wsign * (1 - offset) * di;				\
  advi2 = ijk + wsign * (2 - offset) * di;				\
  advi3 = ijk + wsign * (3 - offset) * di;				\
  advi4 = ijk + wsign * (4 - offset) * di;				\
  advu0 = - wsign * wx * stencil4[offset][0];				\
  advu1 = - wsign * wx * stencil4[offset][1];				\
  advu2 = - wsign * wx * stencil4[offset][2];				\
  advu3 = - wsign * wx * stencil4[offset][3];				\
  advu4 = - wsign * wx * stencil4[offset][4];				\

#define set_advection4(wx,wy,wz)					\
  rstencil = 4 - advectionlopsided;					\
  set_advection4_onedir(i,di,imax,wx,advi0,advi1,advi2,advi3,advi4,	\
			             advu0,advu1,advu2,advu3,advu4);	\
  set_advection4_onedir(j,dj,jmax,wy,advj0,advj1,advj2,advj3,advj4,	\
			             advv0,advv1,advv2,advv3,advv4);	\
  set_advection4_onedir(k,dk,kmax,wz,advk0,advk1,advk2,advk3,advk4,	\
			             advw0,advw1,advw2,advw3,advw4);    \



/* sixth order: 
   advectionlopsided6 = 0 gives stencil oooooox  (one-sided) 
   advectionlopsided6 = 1 gives stencil oooooxo  (2lop-sided)
   advectionlopsided6 = 2 gives stencil ooooxoo  (lop-sided)
   advectionlopsided6 = 3 gives stencil oooxooo  (centered)
*/
#define set_advection6_onedir(i,di,imax,wx,				      \
			      advi0,advi1,advi2,advi3,advi4,advi5,advi6,      \
			      advu0,advu1,advu2,advu3,advu4,advu5,advu6)      \
  set_advection_wsign(i,di,imax,wx);					      \
  offset += advectionlopsided6;						      \
  advi0 = ijk + wsign * (0 - offset) * di;				      \
  advi1 = ijk + wsign * (1 - offset) * di;				      \
  advi2 = ijk + wsign * (2 - offset) * di;				      \
  advi3 = ijk + wsign * (3 - offset) * di;				      \
  advi4 = ijk + wsign * (4 - offset) * di;				      \
  advi5 = ijk + wsign * (5 - offset) * di;				      \
  advi6 = ijk + wsign * (6 - offset) * di;				      \
  advu0 = - wsign * wx * stencil6[offset][0];			      	      \
  advu1 = - wsign * wx * stencil6[offset][1];				      \
  advu2 = - wsign * wx * stencil6[offset][2];				      \
  advu3 = - wsign * wx * stencil6[offset][3];				      \
  advu4 = - wsign * wx * stencil6[offset][4];				      \
  advu5 = - wsign * wx * stencil6[offset][5];				      \
  advu6 = - wsign * wx * stencil6[offset][6];				      \

#define set_advection6(wx,wy,wz)					         \
  rstencil = 6 - advectionlopsided6;					         \
  set_advection6_onedir(i,di,imax,wx,advi0,advi1,advi2,advi3,advi4,advi5,advi6,	 \
			             advu0,advu1,advu2,advu3,advu4,advu5,advu6); \
  set_advection6_onedir(j,dj,jmax,wy,advj0,advj1,advj2,advj3,advj4,advj5,advj6,	 \
			             advv0,advv1,advv2,advv3,advv4,advv5,advv6); \
  set_advection6_onedir(k,dk,kmax,wz,advk0,advk1,advk2,advk3,advk4,advk5,advk6,  \
			             advw0,advw1,advw2,advw3,advw4,advw5,advw6); \







/* eighth order: */
#define set_advection8_onedir(i,di,imax,wx,				      \
			      advi0,advi1,advi2,advi3,advi4,advi5,advi6,advi7,advi8, \
			      advu0,advu1,advu2,advu3,advu4,advu5,advu6,advu7,advu8) \
  set_advection_wsign(i,di,imax,wx);					      \
  offset += advectionlopsided8;						      \
  advi0 = ijk + wsign * (0 - offset) * di;				      \
  advi1 = ijk + wsign * (1 - offset) * di;				      \
  advi2 = ijk + wsign * (2 - offset) * di;				      \
  advi3 = ijk + wsign * (3 - offset) * di;				      \
  advi4 = ijk + wsign * (4 - offset) * di;				      \
  advi5 = ijk + wsign * (5 - offset) * di;				      \
  advi6 = ijk + wsign * (6 - offset) * di;				      \
  advi7 = ijk + wsign * (7 - offset) * di;				      \
  advi8 = ijk + wsign * (8 - offset) * di;				      \
  advu0 = - wsign * wx * stencil8[offset][0];			      	      \
  advu1 = - wsign * wx * stencil8[offset][1];				      \
  advu2 = - wsign * wx * stencil8[offset][2];				      \
  advu3 = - wsign * wx * stencil8[offset][3];				      \
  advu4 = - wsign * wx * stencil8[offset][4];				      \
  advu5 = - wsign * wx * stencil8[offset][5];				      \
  advu6 = - wsign * wx * stencil8[offset][6];				      \
  advu7 = - wsign * wx * stencil8[offset][7];				      \
  advu8 = - wsign * wx * stencil8[offset][8];				      \

#define set_advection8(wx,wy,wz)					         \
  rstencil = 8 - advectionlopsided8;					         \
  set_advection8_onedir(i,di,imax,wx,advi0,advi1,advi2,advi3,advi4,advi5,advi6,advi7,advi8, \
			             advu0,advu1,advu2,advu3,advu4,advu5,advu6,advu7,advu8); \
  set_advection8_onedir(j,dj,jmax,wy,advj0,advj1,advj2,advj3,advj4,advj5,advj6,advj7,advj8, \
			             advv0,advv1,advv2,advv3,advv4,advv5,advv6,advv7,advv8); \
  set_advection8_onedir(k,dk,kmax,wz,advk0,advk1,advk2,advk3,advk4,advk5,advk6,advk7,advk8, \
			             advw0,advw1,advw2,advw3,advw4,advw5,advw6,advw7,advw8); \



/* tenth order: */
#define set_advection10_onedir(i,di,imax,wx,				      \
			       advi0,advi1,advi2,advi3,advi4,advi5,advi6,advi7,advi8,advi9,advi10, \
			       advu0,advu1,advu2,advu3,advu4,advu5,advu6,advu7,advu8,advu9,advu10) \
  set_advection_wsign(i,di,imax,wx);					      \
  offset += advectionlopsided10;						      \
  advi0 = ijk + wsign * (0 - offset) * di;				      \
  advi1 = ijk + wsign * (1 - offset) * di;				      \
  advi2 = ijk + wsign * (2 - offset) * di;				      \
  advi3 = ijk + wsign * (3 - offset) * di;				      \
  advi4 = ijk + wsign * (4 - offset) * di;				      \
  advi5 = ijk + wsign * (5 - offset) * di;				      \
  advi6 = ijk + wsign * (6 - offset) * di;				      \
  advi7 = ijk + wsign * (7 - offset) * di;				      \
  advi8 = ijk + wsign * (8 - offset) * di;				      \
  advi9 = ijk + wsign * (9 - offset) * di;				      \
  advi10 = ijk + wsign * (10 - offset) * di;				      \
  advu0 = - wsign * wx * stencil10[offset][0];			      	      \
  advu1 = - wsign * wx * stencil10[offset][1];				      \
  advu2 = - wsign * wx * stencil10[offset][2];				      \
  advu3 = - wsign * wx * stencil10[offset][3];				      \
  advu4 = - wsign * wx * stencil10[offset][4];				      \
  advu5 = - wsign * wx * stencil10[offset][5];				      \
  advu6 = - wsign * wx * stencil10[offset][6];				      \
  advu7 = - wsign * wx * stencil10[offset][7];				      \
  advu8 = - wsign * wx * stencil10[offset][8];				      \
  advu9 = - wsign * wx * stencil10[offset][9];				      \
  advu10 = - wsign * wx * stencil10[offset][10];			      \


#define set_advection10(wx,wy,wz)					         \
  rstencil = 10 - advectionlopsided10;					         \
  set_advection10_onedir(i,di,imax,wx,advi0,advi1,advi2,advi3,advi4,advi5,advi6,advi7,advi8,advi9,advi10, \
			              advu0,advu1,advu2,advu3,advu4,advu5,advu6,advu7,advu8,advu9,advu10); \
  set_advection10_onedir(j,dj,jmax,wy,advj0,advj1,advj2,advj3,advj4,advj5,advj6,advj7,advj8,advj9,advj10, \
			              advv0,advv1,advv2,advv3,advv4,advv5,advv6,advv7,advv8,advv9,advv10); \
  set_advection10_onedir(k,dk,kmax,wz,advk0,advk1,advk2,advk3,advk4,advk5,advk6,advk7,advk8,advk9,advk10, \
			              advw0,advw1,advw2,advw3,advw4,advw5,advw6,advw7,advw8,advw9,advw10); \
















/* hrsc_rec1d.c 
   sbernuz 02/2012 */

/*   
     reconstruction methods 
     
     harcoded for uniform grids
     
     rec values umins / uplus are relative to the edges of the ith-cell
     the Riemann problem at the ith-interface (point i+1/2) is thus defined by 
     
     uleft  = uplus (i) 
     uright = umins (i+1)
*/


#include "bam.h"
#include "matter.h"

//
//   macros & numbers 
//

#define PRT 0
#define SQ(a) (((a)*(a)))
#define USE_WEIGHTS_OPTIMAL 0

static double cc [32] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 
			  11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 
			  21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 
			  31};
static double oocc [32] = { 1e99, 1., 1./2., 1./3., 1./4., 1./5., 1./6., 1./7., 1./8., 1./9., 1./10.,
			    1./11., 1./12., 1./13., 1./14., 1./15., 1./16., 1./17., 1./18., 1./19., 1./20., 
			    1./21., 1./22., 1./23., 1./24., 1./25., 1./26., 1./27., 1./28., 1./29., 1./30., 
			    1./31};

const double othreeotwo  = 13./12.;
const double EPSL        = 1e-40; //1e-6
const double alpha       = 0.7; // CENO coef
static double optimw[3]  = {1./10., 3./5., 3./10.};// WENO optimal weights

// Lagrangian coeffs for 5th order
static double cl5[5] = {2./60., -13./60., 47./60., 27./60., -3./60.};
//static double cl5[5] = {2.343750000000000e-02, -1.562500000000000e-01, 7.031250000000000e-01, 4.687500000000000e-01, -3.906250000000000e-02};

// Lagrangian coeffs for 6th order
static double cl6[6] = {3./256.,-25./256.,75./128.,75./128.,-25./256.,3./256.};



//
// limiters 
//


// MINMOD 
double MM2( double x, double y ) 
{ 
  double s1 = SIGN( cc[1], x );
  double s2 = SIGN( cc[1], y );
  
  return( oocc[2] * (s1+s2) * DMIN(fabs(x),fabs(y)) );
}


double MM3( double x, double y, double z ) 
{ 
  double s1  = SIGN( cc[1], x );
  double s2  = SIGN( cc[1], y );
  double s3  = SIGN( cc[1], z );
  double s   = fabs(s1+s3);
  double min = DMIN( (fabs(x)), (fabs(y)) );
  
  return( oocc[4]*(s1+s2)*s*DMIN( min, (fabs(z)) ) );
}


double MM4( double x, double y, double z, double w ) 
{ 
  
  double s1 = SIGN( cc[1], x );
  double s2 = SIGN( cc[1], y );
  double s3 = SIGN( cc[1], z );
  double s4 = SIGN( cc[1], w );
  double s  = fabs( (s1+s3)*(s1+s4) );
  double m1 = DMIN( (fabs(x)), (fabs(y)) ); 
  double m2 = DMIN( (fabs(z)), (fabs(w)) ); 
  
  return( oocc[8]*(s1+s2)*s*DMIN( m1, m2 ) );  
}


// VanLeer MC limiter 
double MC2( double x, double y ) 
{ 
  
  double s1  = SIGN( cc[1], x );
  double s2  = SIGN( cc[1], y );
  double min = DMIN( (cc[2]*fabs(x)), (cc[2]*fabs(y)) );
  
  return( oocc[2]*(s1+s2) * DMIN( min, (oocc[2]*fabs(x+y)) ) );  
}

double Median(double a, double b, double c)
{
  return (a + MM2(b-a,c-a));
}


// CENO3 
double ceno3lim( double d[3] )
{
    
    double o3term = 0.0;  
    double absd[3];  
    int kmin;
    
    if ( ((d[0]>=0.) && (d[1]>=0.) && (d[2]>=0.)) || 
	 ((d[0]<0.) && (d[1]<0.) && (d[2]<0.))  ) {
	
	absd[0] = fabs( d[0] );
	absd[1] = fabs( alpha*d[1] );
	absd[2] = fabs( d[2] );
	
	kmin = 0;
	if( absd[1] < absd[kmin] ) kmin = 1;
	if( absd[2] < absd[kmin] ) kmin = 2;
	
	o3term = d[kmin];
	
    } 
    
    return( o3term );
    
}

// MP5
double mp5lim( double u, 
	       double uimt,double uimo,double ui,double uipo,double uipt )
{

  double ump,umin,umax, fUL,fAV,fMD,fLC;
  double d2m,d2c,d2p, dMMm,dMMp; 
  double tmp1,tmp2, mp5;

  /*static double alpha = 2.0;*/ 
  static double alpha = 4.0; 
  static double eps = 1e-20;
  static double fot = 1.3333333333333333333; // 4/3

  mp5 = u;
  ump = ui + MM2( (uipo-ui), alpha*(ui-uimo) ); 
  
  if( ((u-ui)*(u-ump)) <= eps ) return mp5;
  
  d2m = uimt -2.*uimo +ui;
  d2c = uimo -2.*ui   +uipo;
  d2p = ui   -2.*uipo +uipt;
  
  tmp1 = MM2(4.*d2c - d2p, 4.*d2p - d2c);
  tmp2 = MM2(d2c, d2p);
  dMMp = MM2(tmp1,tmp2);   

  tmp1 = MM2(4.*d2m - d2c, 4.*d2c - d2m);
  tmp2 = MM2(d2c, d2m);
  dMMm = MM2(tmp1,tmp2); 
  
  fUL = ui + alpha*(ui-uimo);
  fAV = 0.5*(ui + uipo);
  fMD = fAV - 0.5*dMMp; 
  fLC = 0.5*(3.0*ui - uimo) + fot*dMMm;  
  
  tmp1 = DMIN(ui, uipo); tmp1 = DMIN(tmp1, fMD);
  tmp2 = DMIN(ui, fUL);  tmp2 = DMIN(tmp2, fLC);
  umin = DMAX(tmp1, tmp2); 
  
  tmp1 = DMAX(ui, uipo); tmp1 = DMAX(tmp1, fMD);
  tmp2 = DMAX(ui, fUL);  tmp2 = DMAX(tmp2, fLC);
  umax = DMIN(tmp1, tmp2); 
  
  mp5 = Median(mp5, umin, umax); 
      
  return mp5;  
} 



//
// ptwise rec routines, one side: plus (left interface) 
//


// constant state 
double rec1d_p_godunov(double *u, int i)
{
  return u[i];
}


// average 
double rec1d_p_avg( double *u, int i )
{
  return ( oocc[2]* (u[i] + u[i+1]) );  
}


// lin TVD 
double rec1d_p_lintvd(double *u, int i)
{
  double up;

  double slope = oocc[2] * MATTER.hrsc_TVDlim( ( u [i] - u [i-1] ), ( u [i+1] - u [i] ) ); 
  up = u [i] + slope; 

  return up;
}

// CENO3 
double rec1d_p_ceno3(double *u, int i)
{  
  double up;

  double uipt = u [i+2];
  double uipo = u [i+1];
  double ui   = u [i];
  double uimo = u [i-1];
  double uimt = u [i-2];

  double slope = oocc[2] * MATTER.hrsc_TVDlim( ( ui - uimo ), ( uipo - ui ) );   
  
  double tmpL;
  double tmpd[3];    // these are d^k_i with k = -1,0,1 

  tmpL    = ui + slope;
  tmpd[0] = ( cc[3]*uimt - cc[10]*uimo + cc[15]*ui   )*oocc[8] - tmpL;
  tmpd[1] = ( -     uimo + cc[6] *ui   + cc[3] *uipo )*oocc[8] - tmpL;
  tmpd[2] = ( cc[3]*ui   + cc[6] *uipo -        uipt )*oocc[8] - tmpL;
 
  up  = tmpL + ceno3lim(tmpd); 

  return up;
}

// MUSCL
double rec1d_p_muscl( double *u,int i) {

  double d  = u[i+1] - u[i];
  double dm = u[i]   - u[i-1];

  double rp, rm;

  if (dm!=0) rp = d/dm;
  else       rp = 0;

  if (d!=0) rm = dm/d;
  else      rm = 0;

  return u[i] + oocc[6]*MM2(1,rp)*dm + oocc[3]*MM2(1,rm)*d;

}

// Lagrangian 4th 
double rec1d_p_lag4( double *u, int i)
{
  return ( ( cc[9]*(u[i] + u[i+1]) - u[i-1] - u[i+2] )*oocc[16] );
}

// PPM4
double rec1d_p_ppm4(double *u, int i) 
{
  double delta_u[3];

  delta_u[0] = PPM_MCslope_pt ( u, i-1 );
  delta_u[1] = PPM_MCslope_pt ( u, i   );
  delta_u[2] = PPM_MCslope_pt ( u, i+1 );

  double uplus = oocc[2]*(u[i] + u[i+1]) + (delta_u[1] - delta_u[2])*oocc[6];
  double umins = oocc[2]*(u[i-1] + u[i]) + (delta_u[0] - delta_u[1])*oocc[6];

  PPM_monoton1D_pt ( u, &umins, &uplus, i );

  return uplus;

}

// WENO5 
double rec1d_p_weno5(double *u, int i)
{
  double up;

  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];
  
  double uk[3], a[3], b[3], w[3], dsa;
  int j;
  
#if (USE_WEIGHTS_OPTIMAL)
  for( j = 0 ; j<3; j++) w[j] = optimw[j];
#else
  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uimt-cc[2]*uimo+ui   )) + oocc[4] * SQ((uimt-cc[4]*uimo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uimo-cc[2]*ui  +uipo )) + oocc[4] * SQ((uimo-uipo ) ); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uipo+uipt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uipo+uipt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]/( SQ(( EPSL + b[j])) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
#endif  

  uk[0] = oocc[6]*(   cc[2]*uimt - cc[7]*uimo + cc[11]*ui   );
  uk[1] = oocc[6]*( - cc[1]*uimo + cc[5]*ui   + cc[2] *uipo );
  uk[2] = oocc[6]*(   cc[2]*ui   + cc[5]*uipo -        uipt );
  
  up = w[0]*uk[0] + w[1]*uk[1] + w[2]*uk[2];

  return up;
}

// WENOZ 
double rec1d_p_wenoz(double *u, int i)
{
  double up;

  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];
  
  double uk[3], a[3], b[3], w[3], dsa;
  int j;
  
#if (USE_WEIGHTS_OPTIMAL)
  for( j = 0 ; j<3; j++) w[j] = optimw[j];
#else
  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uimt-cc[2]*uimo+ui   )) + oocc[4] * SQ((uimt-cc[4]*uimo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uimo-cc[2]*ui  +uipo )) + oocc[4] * SQ((uimo-uipo ) ); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uipo+uipt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uipo+uipt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]*( cc[1] + fabs(b[0]-b[2])/( EPSL + b[j]) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
#endif  

  uk[0] = oocc[6]*(   cc[2]*uimt - cc[7]*uimo + cc[11]*ui   );
  uk[1] = oocc[6]*( - cc[1]*uimo + cc[5]*ui   + cc[2] *uipo );
  uk[2] = oocc[6]*(   cc[2]*ui   + cc[5]*uipo -        uipt );
  
  up = w[0]*uk[0] + w[1]*uk[1] + w[2]*uk[2];

  if(PRT &&  (w[0] < 0. || w[1]<0. || w[2]<0.)) printf("WENOZ plus: up = %le, w[0] = %le, w[1] = %le, w[2] = %le, uk[0] = %le, uk[1] = %le, uk[2] = %le\n",
                                    up, w[0], w[1], w[2], uk[0], uk[1], uk[2]);


  return up;
}

// MP5 rec 
double rec1d_p_mp5(double *u, int i)
{
  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];
  
  double ulim  = cl5[0]*uimt + cl5[1]*uimo + cl5[2]*ui + cl5[3]*uipo + cl5[4]*uipt; 
  return mp5lim(ulim, uimt,uimo,ui,uipo,uipt);
}

// PPM6
double rec1d_p_ppm6(double *u, int i) 
{
  double delta_u[5];

  delta_u[0] = PPM_MCslope_pt ( u, i-2 );
  delta_u[1] = PPM_MCslope_pt ( u, i-1 );
  delta_u[2] = PPM_MCslope_pt ( u, i   );
  delta_u[3] = PPM_MCslope_pt ( u, i+1 );
  delta_u[4] = PPM_MCslope_pt ( u, i+2 );

  double uplus = oocc[2]*( u[i] + u[i+1] ) + ( delta_u[2] - delta_u[3] )*oocc[6] - ( 3*(delta_u[3] - delta_u[2]) - (delta_u[4] - delta_u[1]) )*oocc[30];
  double umins = oocc[2]*( u[i-1] + u[i] ) + ( delta_u[1] - delta_u[2] )*oocc[6] - ( 3*(delta_u[2] - delta_u[1]) - (delta_u[3] - delta_u[0]) )*oocc[30];

  PPM_monoton1D_pt ( u, &umins, &uplus, i );

  return uplus;

}

// Lagrangian 6th 
double rec1d_p_lag6( double *u, int i)
{ 
  return( cl6[0]*u[i-2] + cl6[1]*u[i-1] + cl6[2]*u[i] + cl6[3]*u[i+1] + cl6[4]*u[i+2] + cl6[5]*u[i+3] );
}


//
// ptwise rec routines, one side: mins (right interface) 
// do not repeat symmetric ones !
//


// constant state 
double rec1d_m_godunov(double *u, int i)
{
  return rec1d_p_godunov(u, i);
}


// average 
double rec1d_m_avg( double *u, int i )
{
  return rec1d_p_avg( u, i-1 );
}


// lin TVD 
double rec1d_m_lintvd(double *u, int i)
{
  double um;
  
  double slope = oocc[2] * MATTER.hrsc_TVDlim( ( u [i] - u [i-1] ), ( u [i+1] - u [i] ) );
  um = u [i] - slope;

  return um;
}


// CENO3 
double rec1d_m_ceno3(double *u, int i)
{  
  double um;

  double uipt = u [i+2];
  double uipo = u [i+1];
  double ui   = u [i];
  double uimo = u [i-1];
  double uimt = u [i-2];

  double slope = oocc[2] * MATTER.hrsc_TVDlim( ( ui - uimo ), ( uipo - ui ) );   
  
  double tmpL;
  double tmpd[3];    // these are d^k_i with k = -1,0,1 

  tmpL    = ui - slope;
  tmpd[2] = ( cc[3]*uipt - cc[10]*uipo + cc[15]*ui   )*oocc[8] - tmpL;
  tmpd[1] = ( -     uipo + cc[6] *ui   + cc[3] *uimo )*oocc[8] - tmpL;
  tmpd[0] = ( cc[3]*ui   + cc[6] *uimo -        uimt )*oocc[8] - tmpL;   
  
  um  = tmpL + ceno3lim(tmpd); 

  return um;
}

// MUSCL
double rec1d_m_muscl( double *u,int i) {

  double d  = u[i+1] - u[i];
  double dm = u[i]   - u[i-1];

  double rp, rm;

  if (dm!=0) rp = d/dm;
  else       rp = 0;

  if (d!=0) rm = dm/d;
  else      rm = 0;

  return u[i] - oocc[6]*MM2(1,rp)*dm - oocc[3]*MM2(1,rm)*d;

}

// Lagrangian 4th 
double rec1d_m_lag4( double *u, int i)
{
  return rec1d_p_lag4( u, i-1 );
}

// PPM4
double rec1d_m_ppm4(double *u, int i) 
{
  double delta_u[3];

  delta_u[0] = PPM_MCslope_pt ( u, i-1 );
  delta_u[1] = PPM_MCslope_pt ( u, i   );
  delta_u[2] = PPM_MCslope_pt ( u, i+1 );

  double uplus = oocc[2]*(u[i] + u[i+1]) +  (delta_u[1] - delta_u[2])*oocc[6];
  double umins = oocc[2]*(u[i-1] + u[i]) +  (delta_u[0] - delta_u[1])*oocc[6];

  PPM_monoton1D_pt ( u, &umins, &uplus, i );

  return umins;

}

// WENO5 
double rec1d_m_weno5(double *u, int i)
{
  double um;
  
  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];
  
  double uk[3], a[3], b[3], w[3], dsa;
  int j;
  
#if (USE_WEIGHTS_OPTIMAL)
  for( j = 0 ; j<3; j++) w[j] = optimw[j];
#else
  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uipt-cc[2]*uipo+ui   )) + oocc[4] * SQ((uipt-cc[4]*uipo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uipo-cc[2]*ui  +uimo )) + oocc[4] * SQ((uipo-uimo )); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uimo+uimt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uimo+uimt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]/( SQ(( EPSL + b[j])) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
#endif  
  
  uk[0] = oocc[6]*( cc[2]*uipt - cc[7]*uipo + cc[11]*ui   );
  uk[1] = oocc[6]*( -     uipo + cc[5]*ui   + cc[2] *uimo );
  uk[2] = oocc[6]*( cc[2]*ui   + cc[5]*uimo -        uimt );
  
  um = w[0]*uk[0] + w[1]*uk[1] + w[2]*uk[2];

  return um;
}

// WENOZ 
double rec1d_m_wenoz(double *u, int i)
{
  double um;

  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];
  
  double uk[3], a[3], b[3], w[3], dsa;
  int j;
  
#if (USE_WEIGHTS_OPTIMAL)
  for( j = 0 ; j<3; j++) w[j] = optimw[j];
#else
  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uipt-cc[2]*uipo+ui   )) + oocc[4] * SQ((uipt-cc[4]*uipo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uipo-cc[2]*ui  +uimo )) + oocc[4] * SQ((uipo-uimo )); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uimo+uimt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uimo+uimt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]*( cc[1] + fabs(b[0]-b[2])/( EPSL + b[j]) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
#endif  
  
  uk[0] = oocc[6]*( cc[2]*uipt - cc[7]*uipo + cc[11]*ui   );
  uk[1] = oocc[6]*( -     uipo + cc[5]*ui   + cc[2] *uimo );
  uk[2] = oocc[6]*( cc[2]*ui   + cc[5]*uimo -        uimt );
  
  um = w[0]*uk[0] + w[1]*uk[1] + w[2]*uk[2];

  if(PRT && (w[0] <0. || w[1]<0. || w[2]<0.)) printf("WENOZ mins: um = %le, w[0] = %le, w[1] = %le, w[2] = %le, uk[0] = %le, uk[1] = %le, uk[2] = %le\n",
				    um, w[0], w[1], w[2], uk[0], uk[1], uk[2]);
  return um;
}

// MP5 rec 
double rec1d_m_mp5(double *u, int i)
{
  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];

  double ulim  = cl5[0]*uipt + cl5[1]*uipo + cl5[2]*ui + cl5[3]*uimo + cl5[4]*uimt; 
  return mp5lim(ulim, uipt,uipo,ui,uimo,uimt);
}

// PPM6
double rec1d_m_ppm6(double *u, int i) 
{
  double delta_u[5];

  delta_u[0] = PPM_MCslope_pt ( u, i-2 );
  delta_u[1] = PPM_MCslope_pt ( u, i-1 );
  delta_u[2] = PPM_MCslope_pt ( u, i   );
  delta_u[3] = PPM_MCslope_pt ( u, i+1 );
  delta_u[4] = PPM_MCslope_pt ( u, i+2 );

  double uplus = oocc[2]*( u[i] + u[i+1] ) + ( delta_u[2] - delta_u[3] )*oocc[6] - ( 3*(delta_u[3] - delta_u[2]) - (delta_u[4] - delta_u[1]) )*oocc[30];
  double umins = oocc[2]*( u[i-1] + u[i] ) + ( delta_u[1] - delta_u[2] )*oocc[6] - ( 3*(delta_u[2] - delta_u[1]) - (delta_u[3] - delta_u[0]) )*oocc[30];

  PPM_monoton1D_pt ( u, &umins, &uplus, i );

  return umins;

}

// Lagrangian 6th 
double rec1d_m_lag6( double *u, int i)
{
  return rec1d_p_lag6( u, i-1 );
}


//
// ptwise rec routines, plus and mins
//


// constant state 
void rec1d_godunov( double *u, int i, 
	      double *up, double *um )
{
  *up = *um = u[i]; 
}


// average 
void rec1d_avg( double *u, int i, 
		double *up, double *um )
{
  *up = rec1d_p_avg( u, i );
  *um = rec1d_m_avg( u, i );
}


// linear TVD 
void rec1d_lintvd( double *u, int i, 
		   double *up, double *um )
{
  *up = rec1d_p_lintvd( u, i );
  *um = rec1d_m_lintvd( u, i );
}


// CENO3 
void rec1d_ceno3( double *u, int i, 
	    double *up, double *um )
{
  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];

  double tmpL;  
  double tmpd[3];    // these are d^k_i with k = -1,0,1 
  
  double slope = oocc[2] * MATTER.hrsc_TVDlim( ( ui - uimo ), ( uipo - ui ) );   

  tmpL    = ui + slope;
  tmpd[0] = ( cc[3]*uimt - cc[10]*uimo + cc[15]*ui   )*oocc[8] - tmpL;
  tmpd[1] = ( -     uimo + cc[6] *ui   + cc[3] *uipo )*oocc[8] - tmpL;
  tmpd[2] = ( cc[3]*ui   + cc[6] *uipo -        uipt )*oocc[8] - tmpL;
  
  *up  = tmpL + ceno3lim(tmpd); 
  
  tmpL    = ui - slope;
  tmpd[2] = ( cc[3]*uipt - cc[10]*uipo + cc[15]*ui   )*oocc[8] - tmpL;
  tmpd[1] = ( -     uipo + cc[6] *ui   + cc[3] *uimo )*oocc[8] - tmpL;
  tmpd[0] = ( cc[3]*ui   + cc[6] *uimo -        uimt )*oocc[8] - tmpL;   
  
  *um  = tmpL + ceno3lim(tmpd);    
  
}

// MUSCL
void rec1d_muscl( double *u,int i, double *ul,double *ur ) {

  double dp = u[i+2] - u[1+i];
  double d  = u[i+1] - u[i];
  double dm = u[i]   - u[i-1];

  double r1, r2, r3, r4;

  if (dm!=0) r1 = d/dm;
  else       r1 = 0;

  if (d!=0) {r2 = dm/d; r3 = dp/d;}
  else       r2 = r3 = 0;

  if (dp!=0) r4 = d/dp;
  else       r4 = 0;

  *ul = u[i]   + oocc[6]*MM2(1,r1)*dm + oocc[3]*MM2(1,r2)*d;
  *ur = u[i+1] + oocc[6]*MM2(1,r3)*d  + oocc[3]*MM2(1,r4)*dp;

}

// Lagrangian 4th 
void rec1d_lag4( double *u, int i, 
	   double *up, double *um )
{
  *up = rec1d_p_lag4( u, i );
  *um = rec1d_m_lag4( u, i );
}

// PPM4
void rec1d_ppm4( double *u, int i, 
	   double *up, double *um )
{
  *up = rec1d_p_ppm4( u, i );
  *um = rec1d_m_ppm4( u, i );
}

// WENO5 
void rec1d_weno5( double *u, int i, 
	    double *up, double *um )
{
  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];
  
  double a[3], b[3], w[3], dsa;
  int j;

  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uimt-cc[2]*uimo+ui   )) + oocc[4] * SQ((uimt-cc[4]*uimo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uimo-cc[2]*ui  +uipo )) + oocc[4] * SQ((uimo-uipo ) ); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uipo+uipt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uipo+uipt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]/( SQ(( EPSL + b[j])) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
  
  *up = oocc[6]*( w[0]  *( cc[2]*uimt - cc[7]*uimo + cc[11]*ui   )
		     + w[1]*( -     uimo + cc[5]*ui   + cc[2] *uipo )
		     + w[2]*( cc[2]*ui   + cc[5]*uipo -        uipt ) );
  
  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uipt-cc[2]*uipo+ui   )) + oocc[4] * SQ((uipt-cc[4]*uipo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uipo-cc[2]*ui  +uimo )) + oocc[4] * SQ((uipo-uimo )); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uimo+uimt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uimo+uimt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]/( SQ(( EPSL + b[j])) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
  
  *um = oocc[6]*( w[0]  *( cc[2]*uipt - cc[7]*uipo + cc[11]*ui   )
		     + w[1]*( -     uipo + cc[5]*ui   + cc[2] *uimo )
		     + w[2]*( cc[2]*ui   + cc[5]*uimo -        uimt ) );
  
}


// WENOZ 
void rec1d_wenoz( double *u, int i, 
		  double *up, double *um )
{
  double uimt = u [i-2];
  double uimo = u [i-1];
  double ui   = u [i];
  double uipo = u [i+1];
  double uipt = u [i+2];

  double a[3], b[3], w[3], dsa;
  int j;

  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uimt-cc[2]*uimo+ui   )) + oocc[4] * SQ((uimt-cc[4]*uimo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uimo-cc[2]*ui  +uipo )) + oocc[4] * SQ((uimo-uipo ) ); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uipo+uipt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uipo+uipt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]*( cc[1] + fabs(b[0]-b[2])/( EPSL + b[j]) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
  
  *up = oocc[6]*( w[0]  *( cc[2]*uimt - cc[7]*uimo + cc[11]*ui   )
		     + w[1]*( -     uimo + cc[5]*ui   + cc[2] *uipo )
		     + w[2]*( cc[2]*ui   + cc[5]*uipo -        uipt ) );
  
  // smoothness coefs, Jiag & Shu '96
  b[0] = othreeotwo * SQ(( uipt-cc[2]*uipo+ui   )) + oocc[4] * SQ((uipt-cc[4]*uipo+cc[3]*ui)); 
  b[1] = othreeotwo * SQ(( uipo-cc[2]*ui  +uimo )) + oocc[4] * SQ((uipo-uimo )); 
  b[2] = othreeotwo * SQ(( ui  -cc[2]*uimo+uimt )) + oocc[4] * SQ((cc[3]*ui-cc[4]*uimo+uimt)); 
  
  for( j = 0 ; j<3; j++) a[j] = optimw[j]*( cc[1] + fabs(b[0]-b[2])/( EPSL + b[j]) ); 
  dsa = cc[1]/( a[0] + a[1] + a[2] );
  for( j = 0 ; j<3; j++) w[j] = a[j] * dsa;
  
  *um = oocc[6]*( w[0]  *( cc[2]*uipt - cc[7]*uipo + cc[11]*ui   )
		     + w[1]*( -     uipo + cc[5]*ui   + cc[2] *uimo )
		     + w[2]*( cc[2]*ui   + cc[5]*uimo -        uimt ) ); 

}

// MP5
void rec1d_mp5( double *u, int i, 
		double *up, double *um )
{
  *up = rec1d_p_mp5( u, i );
  *um = rec1d_m_mp5( u, i );
}

// PPM6
void rec1d_ppm6( double *u, int i, 
	   double *up, double *um )
{
  *up = rec1d_p_ppm6( u, i );
  *um = rec1d_m_ppm6( u, i );
}

// Lagrangian 6th 
void rec1d_lag6( double *u, int i, 
	   double *up, double *um )
{
  *up = rec1d_p_lag6( u, i );
  *um = rec1d_m_lag6( u, i );
}











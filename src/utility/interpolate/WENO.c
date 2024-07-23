/* interpolate_Weno.c */
/* mth 01/10 */


#include "bam.h"
#include "interpolate.h"

#define sqr(x) ((x)*(x))








double  interpolate_WENO_6(int N, double x, double x0, double h, double *c, double *u)
{
    // J Sci Comput (2008) 35: 219–240
    // Colin B. Macdonald, Steven J. Ruth
    
    int i;
    double p[4],C[4],IS[4],alpha[4],w[6];
    double ooh  = 1./h;
    double ooh2 = ooh/h;
    double ooh3 = ooh2/h;
    
    double dx0 = x - x0;        // x - x_i-2 
    double dx1 = x - x0-1.*h;   // x - x_i-1
    double dx2 = x - x0-2.*h;   // x - x_i
    double dx3 = x - x0-3.*h;   // x - x_i+1
    double dx4 = x - x0-4.*h;   // x - x_i+2
    double dx5 = x - x0-5.*h;   // x - x_i+3
    double epsl = 1e-6;
    
    //i-2 == 0 ...
    p[1] = u[0] + (u[1]-u[0])*ooh*dx0 + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx0*dx1 + (u[3]-3.*u[2]+3.*u[1]-u[0])*ooh3/6.*dx0*dx1*dx2;
    p[2] = u[1] + (u[2]-u[1])*ooh*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx2 + (u[4]-3.*u[3]+3.*u[2]-u[1])*ooh3/6.*dx1*dx2*dx3;
    p[3] = u[2] + (u[3]-u[2])*ooh*dx2 + (u[4]-2.*u[3]+u[2])*ooh2/2.*dx2*dx3 + (u[5]-3.*u[4]+3.*u[3]-u[2])*ooh3/6.*dx2*dx3*dx4;
    
    C[1] = dx4*dx5*ooh2/20.;
    C[2] =-dx5*dx0*ooh2/10.;
    C[3] = dx0*dx1*ooh2/20.;
    
    IS[1] = (814.*sqr(u[3]) + 4326.*sqr(u[2]) + 2976.*sqr(u[1]) + 244.*sqr(u[0]) - 3579.*u[2]*u[3] - 6927.*u[2]*u[1]
            + 1854.*u[2]*u[0] + 2634.*u[3]*u[1] - 683.*u[3]*u[0] - 1659.*u[1]*u[0])/180.;
    IS[2] = (1986.*sqr(u[3]) + 1986.*sqr(u[2]) + 244.*sqr(u[1]) + 244.*sqr(u[4]) + 1074.*u[2]*u[4] - 3777.*u[2]*u[3]
            - 1269.*u[2]*u[1] + 1074.*u[3]*u[1] - 1269.*u[4]*u[3] - 293.*u[4]*u[1])/180.;
    IS[3] = (814.*sqr(u[2]) + 4326.*sqr(u[3]) + 2976.*sqr(u[4]) + 244.*sqr(u[5]) - 683.*u[2]*u[5] + 2634.*u[2]*u[4]
            - 3579.*u[2]*u[3] - 6927.*u[3]*u[4] + 1854.*u[3]*u[5] - 1659.*u[4]*u[5])/180.;

    double sumalpha = 0;
    for (i=1; i<=3; i++) {
        alpha[i] = C[i]/sqr(epsl+IS[i]);
        sumalpha += alpha[i];
    }
    
    for (i=1; i<=3; i++)
        w[i] = alpha[i]/sumalpha;
    
    double sum = 0;
    for (i=1; i<=3; i++)
        sum += w[i]*p[i];
    
    return sum;
}

double  interpolate_WENO_6_withhack(int N, double x, double x0, double h, double *c, double *u)
{
    // J Sci Comput (2008) 35: 219–240
    // Colin B. Macdonald, Steven J. Ruth
    
    int i;
    int hack[4];
    double p[4],C[4],IS[4],alpha[4],w[6], umax,umin;
    double ooh  = 1./h;
    double ooh2 = ooh*ooh;
    double ooh3 = ooh2*ooh;
    
    double dx0 = x - x0;        // x - x_i-2 
    double dx1 = x - x0-1.*h;   // x - x_i-1
    double dx2 = x - x0-2.*h;   // x - x_i
    double dx3 = x - x0-3.*h;   // x - x_i+1
    double dx4 = x - x0-4.*h;   // x - x_i+2
    double dx5 = x - x0-5.*h;   // x - x_i+3
    double epsl = 1e-6;
    
    //i-2 == 0 ...
    p[1] = u[0] + (u[1]-u[0])*ooh*dx0 + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx0*dx1 + (u[3]-3.*u[2]+3.*u[1]-u[0])*ooh3/6.*dx0*dx1*dx2;
    p[2] = u[1] + (u[2]-u[1])*ooh*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx2 + (u[4]-3.*u[3]+3.*u[2]-u[1])*ooh3/6.*dx1*dx2*dx3;
    p[3] = u[2] + (u[3]-u[2])*ooh*dx2 + (u[4]-2.*u[3]+u[2])*ooh2/2.*dx2*dx3 + (u[5]-3.*u[4]+3.*u[3]-u[2])*ooh3/6.*dx2*dx3*dx4;
    
    C[1] = dx4*dx5*ooh2/20.;
    C[2] =-dx5*dx0*ooh2/10.;
    C[3] = dx0*dx1*ooh2/20.;
    
    IS[1] = (814.*sqr(u[3]) + 4326.*sqr(u[2]) + 2976.*sqr(u[1]) + 244.*sqr(u[0]) - 3579.*u[2]*u[3] - 6927.*u[2]*u[1]
            + 1854.*u[2]*u[0] + 2634.*u[3]*u[1] - 683.*u[3]*u[0] - 1659.*u[1]*u[0])/180.;
    IS[2] = (1986.*sqr(u[3]) + 1986.*sqr(u[2]) + 244.*sqr(u[1]) + 244.*sqr(u[4]) + 1074.*u[2]*u[4] - 3777.*u[2]*u[3]
            - 1269.*u[2]*u[1] + 1074.*u[3]*u[1] - 1269.*u[4]*u[3] - 293.*u[4]*u[1])/180.;
    IS[3] = (814.*sqr(u[2]) + 4326.*sqr(u[3]) + 2976.*sqr(u[4]) + 244.*sqr(u[5]) - 683.*u[2]*u[5] + 2634.*u[2]*u[4]
            - 3579.*u[2]*u[3] - 6927.*u[3]*u[4] + 1854.*u[3]*u[5] - 1659.*u[4]*u[5])/180.;

    for (i=1; i<=3; i++) 
        alpha[i] = C[i]/sqr(epsl+IS[i]);
    
    // additional limiter and linear interpolation if necessary
    umax = u[0];
    umin = u[0];
    for (i=1; i<6; i++) {
        umax = DMAX(umax,u[i]);
        umin = DMIN(umin,u[i]);
    }
    
    for (i=1; i<=3; i++) {
        if ( (umax-p[i])*(p[i]-umin) < 0. ) {
            alpha[i] = 0.;
            hack[i] = 1;
        } else {
            hack[i] = 0;
        }
    }
    
    double sumalpha = 0.;
    for (i=1; i<=3; i++)
        sumalpha += alpha[i];
    
    double sum = 0.;
    if (((hack[1]+hack[2]+hack[3] == 3) || (!finite(sumalpha)) || (sumalpha == 0.)) ) {
        sum = u[3] + ooh*dx3*(u[3]-u[2]);
        //if (u[2]!=u[3]) printf("%e %e   %e %e %e\n", x0,x, u[2], sum, u[3]);
        //sum = interpolate_lagrange_N(N, x, x0, h, c,u);
    } else {
        for (i=1; i<=3; i++)
            w[i] = alpha[i]/sumalpha;
        for (i=1; i<=3; i++)
            sum += w[i]*p[i];
    }
    
    if (!(finite(sum))) {
        printf("%e %e   %e %e\n",sum,sumalpha,u[2],u[3]);
        errorexit("WENO interpolation gives nan");
    }
    
    return sum;
}

double  interpolate_WENO_4(int N, double x, double x0, double h, double *c, double *u)
{
    // J Sci Comput (2008) 35: 219–240
    // Colin B. Macdonald, Steven J. Ruth
    
    int i;
    double p[3],C[3],IS[3],alpha[3],w[4];
    double ooh  = 1./h;
    double ooh2 = ooh/h;
    double ooh3 = ooh2/h;
    
    double dx0 = x - x0;        // x - x_i-1 
    double dx1 = x - x0-1.*h;   // x - x_i
    double dx2 = x - x0-2.*h;   // x - x_i+1
    double dx3 = x - x0-3.*h;   // x - x_i+2
    double epsl = 1e-6;
    
    //i-2 == 0 ...
    p[1] = u[1] + (u[2]-u[0])*ooh/2.*dx1             + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx1*dx1;
    p[2] = u[1] + (-u[3]+4.*u[2]-3.*u[1])*ooh/2.*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx1;
    
    C[1] =-dx3*ooh2/3.;
    C[2] = dx0*ooh2/3.;
    
    IS[1] = (25.*sqr(u[2]) + 64.*sqr(u[1]) + 13.*sqr(u[0]) + 
            26.*u[2]*u[0] - 52.*u[1]*u[0] - 76.*u[2]*u[1])/12.;
    IS[2] = (25.*sqr(u[1]) + 64.*sqr(u[2]) + 13.*sqr(u[3]) + 
            26.*u[3]*u[1] - 52.*u[3]*u[2] - 76.*u[2]*u[1])/12.;
    
    double sumalpha = 0;
    for (i=1; i<=2; i++) {
        alpha[i] = C[i]/sqr(epsl+IS[i]);
        sumalpha += alpha[i];
    }
    
    for (i=1; i<=2; i++)
        w[i] = alpha[i]/sumalpha;
    
    double sum = 0;
    for (i=1; i<=2; i++)
        sum += w[i]*p[i];
    
    return sum;
}

double  interpolate_WENO_4_withhack(int N, double x, double x0, double h, double *c, double *u)
{
    // J Sci Comput (2008) 35: 219–240
    // Colin B. Macdonald, Steven J. Ruth
    
    int i;
    int hack[3];
    double p[3],C[3],IS[3],alpha[3],w[4], umax,umin;
    double ooh  = 1./h;
    double ooh2 = ooh/h;
    double ooh3 = ooh2/h;
    
    double dx0 = x - x0;        // x - x_i-1 
    double dx1 = x - x0-1.*h;   // x - x_i
    double dx2 = x - x0-2.*h;   // x - x_i+1
    double dx3 = x - x0-3.*h;   // x - x_i+2
    double epsl = 1e-6;
    
    //i-2 == 0 ...
    p[1] = u[1] + (u[2]-u[0])*ooh/2.*dx1             + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx1*dx1;
    p[2] = u[1] + (-u[3]+4.*u[2]-3.*u[1])*ooh/2.*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx1;
    
    C[1] =-dx3*ooh2/3.;
    C[2] = dx0*ooh2/3.;
    
    IS[1] = (25.*sqr(u[2]) + 64.*sqr(u[1]) + 13.*sqr(u[0]) + 
            26.*u[2]*u[0] - 52.*u[1]*u[0] - 76.*u[2]*u[1])/12.;
    IS[2] = (25.*sqr(u[1]) + 64.*sqr(u[2]) + 13.*sqr(u[3]) + 
            26.*u[3]*u[1] - 52.*u[3]*u[2] - 76.*u[2]*u[1])/12.;
    
    for (i=1; i<=2; i++) 
        alpha[i] = C[i]/sqr(epsl+IS[i]);
    
    // additional limiter and linear interpoaltion if necessary
    umax = u[0];
    umin = u[0];
    for (i=1; i<4; i++) {
        umax = DMAX(umax,u[i]);
        umin = DMIN(umin,u[i]);
    }
    
    for (i=1; i<=2; i++) {
        if ( (umax-p[i])*(p[i]-umin) < 0. ) {
            alpha[i] = 0.;
            hack[i] = 1;
        } else {
            hack[i] = 0;
        }
    }
    
    double sumalpha = 0;
    for (i=1; i<=2; i++)
        sumalpha += alpha[i];
    
    double sum = 0.;
    if (((hack[1]+hack[2] == 2) || (!finite(sumalpha)) || (sumalpha == 0.)) ) {
        sum = u[2] + ooh*dx2*(u[2]-u[1]);
        //if (u[1]!=u[2]) printf("%e %e   %e %e %e\n", x0,x, u[1], sum, u[2]);
        //sum = interpolate_lagrange_N(N, x, x0, h, c,u);
    } else {
        for (i=1; i<=2; i++)
            w[i] = alpha[i]/sumalpha;
        for (i=1; i<=2; i++)
            sum += w[i]*p[i];
    }
    
    
    if (!finite(u[0]) || !finite(u[1]) || !finite(u[2]) || !finite(u[3]))
        printf("%e %e %e %e\n",u[0],u[1],u[2],u[3]);
    
    if (!(finite(sum))) {
        printf("%e %e   %e %e\n",sum,sumalpha,u[2],u[3]);
        errorexit("WENO interpolation gives nan");
    }
    
    return sum;
}

double  interpolate_ENO_3(int N, double x, double x0, double h, double *c, double *u)
{
    // test ... this is only for centered stencils
    // => has to be modified to work with BAM   
    
    double diffleft,diffright, eno1d;
    
    diffleft  = u[0] + u[2] - 2.*u[1];
    diffright = u[1] + u[3] - 2.*u[2];
    
    if (fabs(diffleft) > fabs(diffright))
        eno1d = 0.125 * (-u[0] + 6.*u[1] + 3.*u[2]);
    else
        eno1d = 0.125 * (3.*u[3] + 6.*u[2] - u[3]);
    
    if ( ((eno1d-u[1]) * (u[2]-eno1d) > 0.) || (diffleft*diffright >= 0.) )
        eno1d = 0.5 * (u[1] + u[2]);
    
    return eno1d;
}


/* not tested so far*/
double  interpolate_WENOZ_6(int N, double x, double x0, double h, double *c, double *u)
{
    // J Sci Comput (2008) 227: 3191-3211
    // Borges R., Cormana M., Costa B., Sun Don W. 
    

    int i; 
    double p[4],C[4],IS[4],alpha[4],w[6];
    double ooh  = 1./h;
    double ooh2 = ooh*ooh;
    double ooh3 = ooh2*ooh;
    
    double dx0 = x - x0;        // x - x_i-2 
    double dx1 = x - x0-1.*h;   // x - x_i-1
    double dx2 = x - x0-2.*h;   // x - x_i
    double dx3 = x - x0-3.*h;   // x - x_i+1
    double dx4 = x - x0-4.*h;   // x - x_i+2
    double dx5 = x - x0-5.*h;   // x - x_i+3
    double epsl = 1e-40;
    
    //i-2 == 0 ...
    p[1] = u[0] + (u[1]-u[0])*ooh*dx0 + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx0*dx1 + (u[3]-3.*u[2]+3.*u[1]-u[0])*ooh3/6.*dx0*dx1*dx2;
    p[2] = u[1] + (u[2]-u[1])*ooh*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx2 + (u[4]-3.*u[3]+3.*u[2]-u[1])*ooh3/6.*dx1*dx2*dx3;
    p[3] = u[2] + (u[3]-u[2])*ooh*dx2 + (u[4]-2.*u[3]+u[2])*ooh2/2.*dx2*dx3 + (u[5]-3.*u[4]+3.*u[3]-u[2])*ooh3/6.*dx2*dx3*dx4;
    
    C[1] = dx4*dx5*ooh2/20.;
    C[2] =-dx5*dx0*ooh2/10.;
    C[3] = dx0*dx1*ooh2/20.;
    
    IS[1] = (814.*sqr(u[3]) + 4326.*sqr(u[2]) + 2976.*sqr(u[1]) + 244.*sqr(u[0]) - 3579.*u[2]*u[3] - 6927.*u[2]*u[1]
            + 1854.*u[2]*u[0] + 2634.*u[3]*u[1] - 683.*u[3]*u[0] - 1659.*u[1]*u[0])/180.;
    IS[2] = (1986.*sqr(u[3]) + 1986.*sqr(u[2]) + 244.*sqr(u[1]) + 244.*sqr(u[4]) + 1074.*u[2]*u[4] - 3777.*u[2]*u[3]
            - 1269.*u[2]*u[1] + 1074.*u[3]*u[1] - 1269.*u[4]*u[3] - 293.*u[4]*u[1])/180.;
    IS[3] = (814.*sqr(u[2]) + 4326.*sqr(u[3]) + 2976.*sqr(u[4]) + 244.*sqr(u[5]) - 683.*u[2]*u[5] + 2634.*u[2]*u[4]
            - 3579.*u[2]*u[3] - 6927.*u[3]*u[4] + 1854.*u[3]*u[5] - 1659.*u[4]*u[5])/180.;

    for (i=1; i<=3; i++) 
        alpha[i] = C[i]*(1+fabs(IS[0]-IS[2]))/(epsl+IS[i]);
       
    
    double sumalpha = 0.;
    for (i=1; i<=3; i++)
        sumalpha += alpha[i];
    
    double sum = 0.;

        for (i=1; i<=3; i++)
            w[i] = alpha[i]/sumalpha;
        for (i=1; i<=3; i++)
            sum += w[i]*p[i];
    
    
    if (!(finite(sum))) {
        printf("%e %e   %e %e\n",sum,sumalpha,u[2],u[3]);
        errorexit("WENO interpolation gives nan");
    }
    
    return sum;

}

/* not tested so far*/
double  interpolate_WENOZ_6_withhack(int N, double x, double x0, double h, double *c, double *u)
{
    // J Sci Comput (2008) 227: 3191-3211
    // Borges R., Cormana M., Costa B., Sun Don W. 
    

    int i;
    int hack[4];
    double p[4],C[4],IS[4],alpha[4],w[6], umax,umin;
    double ooh  = 1./h;
    double ooh2 = ooh*ooh;
    double ooh3 = ooh2*ooh;
    
    double dx0 = x - x0;        // x - x_i-2 
    double dx1 = x - x0-1.*h;   // x - x_i-1
    double dx2 = x - x0-2.*h;   // x - x_i
    double dx3 = x - x0-3.*h;   // x - x_i+1
    double dx4 = x - x0-4.*h;   // x - x_i+2
    double dx5 = x - x0-5.*h;   // x - x_i+3
    double epsl = 1e-40;
    
    //i-2 == 0 ...
    p[1] = u[0] + (u[1]-u[0])*ooh*dx0 + (u[2]-2.*u[1]+u[0])*ooh2/2.*dx0*dx1 + (u[3]-3.*u[2]+3.*u[1]-u[0])*ooh3/6.*dx0*dx1*dx2;
    p[2] = u[1] + (u[2]-u[1])*ooh*dx1 + (u[3]-2.*u[2]+u[1])*ooh2/2.*dx1*dx2 + (u[4]-3.*u[3]+3.*u[2]-u[1])*ooh3/6.*dx1*dx2*dx3;
    p[3] = u[2] + (u[3]-u[2])*ooh*dx2 + (u[4]-2.*u[3]+u[2])*ooh2/2.*dx2*dx3 + (u[5]-3.*u[4]+3.*u[3]-u[2])*ooh3/6.*dx2*dx3*dx4;
    
    C[1] = dx4*dx5*ooh2/20.;
    C[2] =-dx5*dx0*ooh2/10.;
    C[3] = dx0*dx1*ooh2/20.;
    
    IS[1] = (814.*sqr(u[3]) + 4326.*sqr(u[2]) + 2976.*sqr(u[1]) + 244.*sqr(u[0]) - 3579.*u[2]*u[3] - 6927.*u[2]*u[1]
            + 1854.*u[2]*u[0] + 2634.*u[3]*u[1] - 683.*u[3]*u[0] - 1659.*u[1]*u[0])/180.;
    IS[2] = (1986.*sqr(u[3]) + 1986.*sqr(u[2]) + 244.*sqr(u[1]) + 244.*sqr(u[4]) + 1074.*u[2]*u[4] - 3777.*u[2]*u[3]
            - 1269.*u[2]*u[1] + 1074.*u[3]*u[1] - 1269.*u[4]*u[3] - 293.*u[4]*u[1])/180.;
    IS[3] = (814.*sqr(u[2]) + 4326.*sqr(u[3]) + 2976.*sqr(u[4]) + 244.*sqr(u[5]) - 683.*u[2]*u[5] + 2634.*u[2]*u[4]
            - 3579.*u[2]*u[3] - 6927.*u[3]*u[4] + 1854.*u[3]*u[5] - 1659.*u[4]*u[5])/180.;

    for (i=1; i<=3; i++) 
        alpha[i] = C[i]*(1+fabs(IS[0]-IS[2]))/(epsl+IS[i]);
    
    // additional limiter and linear interpolation if necessary
    umax = u[0];
    umin = u[0];
    for (i=1; i<6; i++) {
        umax = DMAX(umax,u[i]);
        umin = DMIN(umin,u[i]);
    }
    
    for (i=1; i<=3; i++) {
        if ( (umax-p[i])*(p[i]-umin) < 0. ) {
            alpha[i] = 0.;
            hack[i] = 1;
        } else {
            hack[i] = 0;
        }
    }
    
    double sumalpha = 0.;
    for (i=1; i<=3; i++)
        sumalpha += alpha[i];
    
    double sum = 0.;
    if (((hack[1]+hack[2]+hack[3] == 3) || (!finite(sumalpha)) || (sumalpha == 0.)) ) {
        sum = u[3] + ooh*dx3*(u[3]-u[2]);
        //if (u[2]!=u[3]) printf("%e %e   %e %e %e\n", x0,x, u[2], sum, u[3]);
        //sum = interpolate_lagrange_N(N, x, x0, h, c,u);
    } else {
        for (i=1; i<=3; i++)
            w[i] = alpha[i]/sumalpha;
        for (i=1; i<=3; i++)
            sum += w[i]*p[i];
    }
    
    if (!(finite(sum))) {
        printf("%e %e   %e %e\n",sum,sumalpha,u[2],u[3]);
        errorexit("WENO interpolation gives nan");
    }
    
    return sum;

}



double interpolate_WENO_N(int N, double x, double x0, double h, double *c, double *u)
{
    if (N==4) 
        return interpolate_WENO_4_withhack(N, x, x0, h, c,u);
    else if (N==6) 
        return interpolate_WENO_6_withhack(N, x, x0, h, c,u);
    else
        errorexit("order is not implemented yet");
    
    return 0.;
}


/* not tested so far*/
double interpolate_WENOZ_N(int N, double x, double x0, double h, double *c, double *u)
{
    if (N==6) 
        return interpolate_WENOZ_6_withhack(N, x, x0, h, c,u);
    else
        errorexit("order is not implemented yet");
    
    return 0.;
}




















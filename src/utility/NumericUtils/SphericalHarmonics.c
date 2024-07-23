
#include "bam.h"







inline double Wigner_d_function(int l, int m, int s, double costheta)
{
    int t;
    double dWig = 0;

    int C1 = s > m  ? 0   : m-s;
    int C2 = m < -s ? l+m : l-s;

    double coshtheta = sqrt(0.5*(1.0 + costheta));
    double sinhtheta = sqrt(0.5*(1.0 - costheta));

    int c1 = C1 < C2 ? C1 : C2;
    int c2 = C1 < C2 ? C2 : C1;
    
    if (c1%2==0) {
        
        for( t = c1; t <= c2; t+=2)
            dWig += pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );

        for( t = c1+1; t <= c2; t+=2)
            dWig -= pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );
        
    } else {
        
        for( t = c1; t <= c2; t+=2)
            dWig -= pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );

        for( t = c1+1; t <= c2; t+=2)
            dWig += pow(coshtheta,2*l+m-s-2*t)*pow(sinhtheta,2*t+s-m)/( fact(l+m-t) * fact(l-s-t) * fact(t) * fact(t+s-m) );

    }

    return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * dWig);
}

inline double plegendre(int l, int m, double x)
{
    double fact,pll,pmm,pmmp1,somx2;
    int i,ll;
    
    pmm=1.0; 

    if (m > 0) {
        somx2=sqrt((1.0-x)*(1.0+x));
        fact=1.0;
    
        for (i=1;i<=m;i++) {
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }
  
    if (l == m) {
        return pmm;
    } else { 
        pmmp1=x*(2*m+1)*pmm;

        if (l == (m+1))
            return pmmp1;
        else {
            for (ll=m+2;ll<=l;ll++) {
                pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                pmm=pmmp1;
                pmmp1=pll;
            }
            return pll;
        }
    }
}

// force instantiation of inline function (in case we compile with -g)
double plegendre(int l, int m, double x);



void SphericalHarmonicY(double *rY, double *iY, int l, int m, double phi, double costheta)
{
    if ((l<m) || (l<0) || (m<0)) 
        errorexit("wrong l or m inside SphericalHarmonicY");
    
    if (1) { 
        
        double Plm = sqrt((2.0*l+1)*ffact(l-m)/(4.*PI*ffact(l+m))) * plegendre(l,m, costheta);
        
        *rY = Plm * cos((double)(m)*phi);
        *iY = Plm * sin((double)(m)*phi);
        
    } else {
        
        double theta = acos(costheta);
        
        SphericalHarmonicYprecomp(rY,iY, l,m, phi,theta);
        
    }
    
    
    
}


void spinweightedSphericalHarmonic(double *rY, double *iY, int l, int m, double phi, double costheta)
{
    if ((l<0) || (m<-l) || (m>l)) 
        errorexit("wrong l or m inside spinweightedSphericalHarmonic");
    
    
    if (0) {
        
        double c = sqrt( (2.0*l+1.)/(4.*PI) );
        
        double dWigner = c*Wigner_d_function(l,m,2.,costheta);

        *rY = cos((double)(m)*phi) * dWigner;
        *iY = sin((double)(m)*phi) * dWigner;
        
    } else {
        
        double theta = acos(costheta);
        
        spinweightedSphericalHarmonicprecomp(rY,iY, l,m, phi,theta);
        
    }
    
    
}












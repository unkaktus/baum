#include "bam.h"



double fact(double n);
double ffact(double n);


//fast factorial for small numbers 

double ffact(double n){


    double f[] = {1,         1,          2,     
        6,         24,         120,
        720,       5040,       40320, 
        362880,    3628800,    39916800,   
        479001600, 6227020800, 87178291200};

    if(n<=14)
        return(f[(int)n]);
    else
        return(fact(n));
}


double fact(double n)
{
    if(n < 0)
        errorexit("Computing a negative factorial");
    
    if (n<=14) {
        return(ffact(n));
    } else {
        n = n*fact(n-1);
        return(n);
    }
}



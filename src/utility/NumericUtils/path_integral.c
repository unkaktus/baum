#include "bam.h"

#define PR 0












/* check if both points are inside level
   => interpolate wanted components at the centre
*/
int     check_points(tL* level, double xc, double yc, double zc, int order, int NInd,int Ind[NInd], double v[NInd])
{
    int i,j;
    double sum;
    int size = bampi_size();
    double send[size], recv[size], vtmp[NInd*size];
    double V[NInd];
    
    // look for points
    int check  = check_interpolation_cube_local_withsym(level, xc,yc,zc, order);
    
    // look if check is 0 on all processors -> BAD
    sum = 0.;
    send[0] = (double)(check);
    bampi_alltoall(send,recv,1);
    for (i=0; i<size; i++) sum += recv[i];
    if (sum<0.5) return 0;
    
    // if point found compute central point and interpolate
    if (check) {
        //printf("interpolate  l=%d   %e %e %e\n",level->l,xc,yc,zc);
        interpolate_xyz_local_minimal_withsym(level, xc,yc,zc, NInd,Ind, V, order, LAGRANGE);
    } else {
        for (i=0; i<NInd; i++) V[i] = 0.;
    }
    
    // send all data to everyone
    bampi_alltoall(V,vtmp,NInd);
    for (i=0; i<NInd; i++) {
        v[i] = 0.;
        for (j=0; j<size; j++) {
            v[i] += recv[j]/sum* vtmp[j*NInd+i];
        }
        //printf("     %e  %e\n",v[i],V[i]);
    }
    //printf("       %e %d\n",sum,check);
    
    return 1;
}

double  compute_length(tG* g, int *l, double x1[4], double x2[4], int order, int NInd,int Ind[NInd])
{
    double v[NInd];
    double xc = 0.5*(x2[1]+x1[1]);
    double yc = 0.5*(x2[2]+x1[2]);
    double zc = 0.5*(x2[3]+x1[3]);
    double dx = (x2[1]-x1[1]);
    double dy = (x2[2]-x1[2]);
    double dz = (x2[3]-x1[3]);
    
    *l = (g->lmax<*l+1) ? g->lmax:*l+1;
    while ((*l>=0) && (!check_points(g->level[*l], xc,yc,zc, order, NInd,Ind, v))) (*l)--;
    if (*l<0) {
        
        //errorexit("points are not inside level");
        printf("points are not inside level");
        return 0.;
    }
    
    double length = v[0]*dx*dx + 2.*v[1]*dx*dy + 2.*v[2]*dx*dz +
                    v[3]*dy*dy + 2.*v[4]*dy*dz +    v[5]*dz*dz;
    
    
    if (length<0.) {
        printf("length<0. (%e) ... set to 0.    (level %d)\n",length,*l);
        printf("  g: %e %e %e %e %e %e\n", v[0],v[1],v[2],v[3],v[4],v[5]);
        printf("  detg: %e\n", detg(v[0],v[1],v[2],v[3],v[4],v[5]));
        //errorexit("");
        return 0.;
    } else {
        return sqrt(length);
    }
}

void    build_orthogonal_vector(double x1[4],double x2[4], double o[4],double p[4])
{
    double s[4] = {0, x2[1]-x1[1],x2[2]-x1[2],x2[3]-x1[3]};
    double s2 = s[1]*s[1] + s[2]*s[2] + s[3]*s[3];
    
    int typ = 1;
    if (fabs(s[1])<1e-12) {
        typ++;
        if (fabs(s[2])<1e-12) typ++;
    }
    
    if (typ==1) {
        p[2] = 1.;
        p[3] = 1.;
        p[1] =-(p[2]*s[2] + p[3]*s[3])/s[1];
        
        o[1] = (p[2]*s[3] - p[3]*s[2])/(s2);
        o[2] = (p[2]*s[2]*s[3] + p[3]*(s[1]*s[1] + s[3]*s[3]))/(s[1]*s2);
        o[3] =-(p[3]*s[2]*s[3] + p[2]*(s[1]*s[1] + s[2]*s[2]))/(s[1]*s2);
    } else if (typ==2) {
        p[1] = 1.;
        p[3] = 1.;
        p[2] =-(p[1]*s[1] + p[3]*s[3])/s[2];
        
        o[1] =-(p[1]*s[1]*s[3] + p[3]*(s[2]*s[2] + s[3]*s[3]))/(s[2]*s2);
        o[2] = (p[1]*s[3] - p[3]*s[1])/(s2);
        o[3] = (p[3]*s[1]*s[3] + p[1]*(s[1]*s[1] + s[2]*s[2]))/(s[2]*s2);
    } else if (typ==3) {
        p[1] = 1.;
        p[2] = 1.;
        p[3] =-(p[1]*s[1] + p[2]*s[2])/s[3];
        
        o[1] = (p[1]*s[1]*s[2] + p[2]*(s[2]*s[2] + s[3]*s[3]))/(s[3]*s2);
        o[2] =-(p[2]*s[1]*s[2] + p[1]*(s[1]*s[1] + s[3]*s[3]))/(s[3]*s2);
        o[3] = (p[1]*s[2] - p[2]*s[1])/(s2);
    }
    
    double o2 = sqrt(s2)/sqrt(o[1]*o[1] + o[2]*o[2] + o[3]*o[3]);
    double p2 = sqrt(s2)/sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
    
    o[1] *= o2;
    o[2] *= o2;
    o[3] *= o2;
    
    p[1] *= p2;
    p[2] *= p2;
    p[3] *= p2;
    
    if (0) {
        printf("%e %e %e\n", s[1],s[2],s[3]);
        printf("%e %e %e\n", o[1],o[2],o[3]);
        printf("%e %e %e\n", p[1],p[2],p[3]);
        
        printf("dotproduct:  %e %e %e\n",s[1]*p[1]+s[2]*p[2]+s[3]*p[3], 
            s[1]*o[1]+s[2]*o[2]+s[3]*o[3],o[1]*p[1]+o[2]*p[2]+o[3]*p[3]);
    
    }
}

double  findmin(double vm, double v0, double vp)
{
    /*
    2nd order function to track extremum
    f = ax^2+bx+c   || x = [-1, 0, 1] , f = [vm, v0, vp]
    a = 0.5*(vp+vm)-v0;
    b = 0.5*(vp-vm);
    c = v0;
    f' == 0         =>  x = -b/(2.*a);
    */
    
    double v;
    double a = 0.5*(vp+vm)-v0;
    double b = 0.5*(vp-vm);
    double c = v0;
    
    if (a <= 1e-14) {/*
        if      (vp<v0) v =  1.;
        else if (vm<v0) v = -1.;
        else           */ v =  0.;
    } else {
        v = -0.5*b/a;
        v = (v> 1.)? 1.:v;
        v = (v<-1.)?-1.:v;
    }
    
    return v;
}





/* computing the proper length between two point with several methods
   assuming ADM gxx is computed and active and both points
   are at the grid
*/
double compute_path_integral(tG *g, double x1[4],double x2[4], char *method)
{
    // NO, we have to use order_RP or - again - changing global interpolation order :(
    // THIS FEELS SOOO WRONG
    int N = 500;
    int order = Geti("order_RP");
    double f = 0.005;
    double accuracy = 1e-6;
    
    
    int i,j,k, l,n;
    double length,dlength,length_old;
    double v[6], dx,dy,dz,xc,yc,zc;
    int Indgxx[6];
    for (n=0; n<6; n++) Indgxx[n] = Ind("adm_gxx")+n;
    double *x[N+1],*xnew[N+1];
    
    
    if (PR) {
        printf("compute path integral (%s)\n",method);
        printf("  from   %2.2e %2.2e %2.2e\n",x1[1],x1[2],x1[3]);
        printf("  to     %2.2e %2.2e %2.2e\n",x2[1],x2[2],x2[3]);
    }
    // set initial points (straight line
    for (n=0; n<=N; n++) {
        x[n]    = (double*) malloc (4*sizeof(double));
        xnew[n] = (double*) malloc (4*sizeof(double));
        for (j=1; j<4; j++) 
            x[n][j] = x1[j] + (double)(n)/(double)(N) * (x2[j]-x1[j]);
    }
    
    
    //use a line
    if (strstr(method,"line")) {
        
        length = 0;
        l = g->lmax;
        for (n=1; n<=N; n++) {
            if (PR) printf("  step %d of %d\n",n,N);
            if (PR) printf("    %2.2e %2.2e %2.2e\n",x[n-1][1],x[n-1][2],x[n-1][3]);
            if (PR) printf("    %2.2e %2.2e %2.2e\n",x[n][1],x[n][2],x[n][3]);
            
            dlength = compute_length(g, &l, x[n-1],x[n], order,6,Indgxx);
            length += dlength;
            
            if (PR) printf("    delta = %e\n",dlength);
        }
        if (PR) printf("  l = %e    l_coord = %e\n", length, 
            sqrt((x2[1]-x1[1])*(x2[1]-x1[1]) + (x2[2]-x1[2])*(x2[2]-x1[2]) + (x2[3]-x1[3])*(x2[3]-x1[3])));
        
        
    // try a line and then minimize length ... THIS HAS TO BE TESTED
    } else if (strstr(method,"minimize")) {
        
        double o1[4],o2[4],p[4], dl[5];
        double a,b,c, d1,d2;
        int steps = 2; //test-value
        
        length = 1e10;
        do {
            
            length_old = length;
            length = 0.;
            l = g->lmax;
            for (n=1; n<N; n++) {
                if (PR) printf("  step %d of %d\n",n,N);
                if (PR) printf("    %2.2e %2.2e %2.2e\n",x[n-1][1],x[n-1][2],x[n-1][3]);
                if (PR) printf("    %2.2e %2.2e %2.2e\n",x[n+1][1],x[n+1][2],x[n+1][3]);
                
                // construct orthorgonal vectors o1 and o2
                build_orthogonal_vector(x[n-1],x[n+1], o2,o1);
                
                dl[0] = compute_length(g, &l, x[n-1],x[n], order,6,Indgxx) + compute_length(g, &l, x[n],x[n+1], order,6,Indgxx);
                
                for (k=1; k<4; k++) p[k] = x[n][k] + f*o1[k];
                dl[1] = compute_length(g, &l, x[n-1],p, order,6,Indgxx) + compute_length(g, &l, p,x[n+1], order,6,Indgxx);
                for (k=1; k<4; k++) p[k] = x[n][k] - f*o1[k];
                dl[2] = compute_length(g, &l, x[n-1],p, order,6,Indgxx) + compute_length(g, &l, p,x[n+1], order,6,Indgxx);
                for (k=1; k<4; k++) p[k] = x[n][k] + f*o2[k];
                dl[3] = compute_length(g, &l, x[n-1],p, order,6,Indgxx) + compute_length(g, &l, p,x[n+1], order,6,Indgxx);
                for (k=1; k<4; k++) p[k] = x[n][k] - f*o2[k];
                dl[4] = compute_length(g, &l, x[n-1],p, order,6,Indgxx) + compute_length(g, &l, p,x[n+1], order,6,Indgxx);
                
                // displacement in o1 and o2 direction
                d1 = findmin(dl[2],dl[0],dl[1]);
                d2 = findmin(dl[4],dl[0],dl[3]);
                
                // add displacement
                for (k=1; k<4; k++) 
                    xnew[n][k] = x[n][k] + f*(d1*o1[k] + d2*o2[k]) ;
                
                if (PR) printf("    delta = %e\n",dlength);
                length += (n==1 || n==N-1)?dl[0]*0.75:dl[0]*0.5;
                
            }
            if (1) printf("  l = %2.16e    l_coord = %2.16e\n", length, 
                sqrt((x2[1]-x1[1])*(x2[1]-x1[1]) + (x2[2]-x1[2])*(x2[2]-x1[2]) + (x2[3]-x1[3])*(x2[3]-x1[3])));
            
            // copy all new points to old ones but do not change start and end point
            for (n=1; n<N; n++) 
                for (k=1; k<4; k++) 
                    x[n][k] = xnew[n][k];
            
            steps--;
        } while ((length_old>length*(1.+accuracy)) && (steps));
        
        
    } else {
        errorexit("method does not exist");
    }
    
    
    
    
    
    
    
    
    if (0) {
        printf("=>  l = %2.16e    l_coord = %2.16e\n\n", length, 
               sqrt((x2[1]-x1[1])*(x2[1]-x1[1]) + (x2[2]-x1[2])*(x2[2]-x1[2]) + (x2[3]-x1[3])*(x2[3]-x1[3])));
        for (n=0; n<=N; n++) 
            printf("%e %e %e\n",x[n][1],x[n][2],x[n][3]);
    }
    
    for (n=0; n<=N; n++) {
        free(x[n]);
        free(xnew[n]);
    }
    
    return length;
}







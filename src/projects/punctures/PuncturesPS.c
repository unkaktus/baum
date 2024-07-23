/* PuncturesPS.c */
/* Bernd Bruegmann 2/2007 */

/* interface to spectral puncture data produced with punctures_ps
   
   1a) call the spectral solver 
       currently only on single processor, requires HYPRE
       -> not implemented yet

   1b) read spectral data from file
       multiple processors, independent of HYPRE and punctures_ps

   2) initialize Cartesian AMR boxes by interpolation

   assumes that each process can read the same input file
*/

#include "bam.h"
#include "ctype.h"
#include "punctures.h"



/* the following is copied from punctures_ps/transform.c
   so that we can read from file independently of that module
*/
static int n1, n2, n3, ntotal, pad;
static double ps_b, ps_dx, ps_rxx, ps_rxy, ps_ryx, ps_ryy;

static void set_ABp();
static void xyz_to_ABp(double x, double y, double z,
		       double *A, double *B, double *phi);
static double ps_u_at_xyz(int order, double x, double y, double z, double *v, int varinter);




/* from cartesian to spectral */
void xyz_to_ABp(double x, double y, double z,
		double *A, double *B, double *phi)
{
  const double s2 = sqrt(2.0);
  const double exp = 3.0/2.0;
 
  double r, rr, xx;
  double t, st, u, su, v, sv, w, sw;
  double check1, check2;
  /* rotate onto x-axis if required */
  w = x;
  x = ps_rxx * w + ps_rxy * y;
  y = ps_ryx * w + ps_ryy * y;

  /* center black holes at +b and -b */
  x -= ps_dx;

  /* offset parameter b rescales the coordinates */
  x /= ps_b;
  y /= ps_b;
  z /= ps_b;

  /* helpers */
  r = sqrt(y*y + z*z);
  rr = r*r;
  xx = x*x;


  /* phi as in cylindrical coordinates about x-axis
     acos covers [0,pi], we need [0,2pi)
  */
  if (dless(0.0, r))
    *phi = (dless(z, 0.0)) ? 2*PI - acos(y/r) : acos(y/r);
  else
    *phi = 0;


  /* r > 0 */
  if (dless(0.0, r)) {

      t = (1+rr)*(1+rr) + 2*(-1 + rr)*xx + xx*xx;
      st = sqrt(t);
      u = 1 - xx + rr*(2 + rr + xx + st) + st;
      su = sqrt(u);
      v = 1 + rr*rr - xx + rr*(2 + xx + st) + st;
      sv = sqrt(v);
      w = 1 + rr - s2*su + st;
      sw = sqrt(w);

    /* x != 0, r > 0 */
    
  if (!dequal(xx/1000, 0.0) && !dequal(w/1000,0.0)) {/* treats every x < 10^-5 as equal to zero  */
/* dtim: if (!dequal(xx/(2000.*(rr/(r+1.)+1.)), 0.0) && !dequal(w/(2000.*(rr/(r+1.)+1.))/1000,0.0)) {
        improves accuracy for high resolution, but unnecessary if you use NR-method */
/* MDH: define these auxillary parameters before the conditional, so that we can test them. */
/*
      t = (1+rr)*(1+rr) + 2*(-1 + rr)*xx + xx*xx;
      st = sqrt(t);
      u = 1 - xx + rr*(2 + rr + xx + st) + st;
      su = sqrt(u);
      v = 1 + rr*rr - xx + rr*(2 + xx + st) + st;
      sv = sqrt(v);
      w = 1 + rr - s2*su + st;
      sw = sqrt(w);
*/
      *A = (2*sw*(1 + rr + st - xx) + s2*sv*(-1 - rr + 2*sw + st - xx))
	   /(4.*r*xx);

      *B = -(sw/x);
    }

    /* x == 0, r > 0 */
    else {
	    *A = (sqrt(1 + rr) - 1)/r + ((sqrt(1 + rr) - 1)*xx)/(2*r*pow((1 + rr),exp));

	    *B = -x/(2*sqrt(1 + rr));

	    if (0) printf(" xx %.16e A %.2e B %.2e \n", xx, *A, *B);		
    }
  }

  /* r == 0 */
  else {

    /* x > 1, r == 0 */
    if (dless(1.0, x)) {
      *A = sqrt(x-1)/sqrt(x+1);
      *B = -1;
    }
    
    /* x < 1, r == 0 */
    else if (dless(x, -1.0)) {
      *A = sqrt(-x-1)/sqrt(-x+1);
      *B = +1;
    }

    /* -1 <= x <= 1, r == 0 */
    else {
      *A = 0;

      /* x != 0 */
      if (!dequal(x, 0.0))
	*B = (sqrt(1-xx) - 1)/x;

      /* x == 0 */
      else
	*B = 0;
    }
  }
}

void puncture_ps_ABp_to_xyz(double A, double B, double phi,
		double *x, double *y, double *z)
{

  double AA = A*A;
  double BB = B*B;
  double rho, xt, yt;

  if (dequal(A, 1.0)) A = 1.0 - 1e-12;

  *x  = ps_b * (AA+1.)/(AA-1.) * 2.*B/(BB+1.)+ ps_dx;

  rho = ps_b * 2.*A/(AA-1.) * (BB-1.)/(BB+1.);

  *y = rho * cos(phi);
  *z = rho * sin(phi);

}


/* set 1d arrays for spectral coordinates 
   see Punctures_functions.c for reference
   special: pad phi direction with ghosts for periodicity
*/
void set_ABp(double *coordA, double *coordB, double *coordphi)
{
  int i;
  double Acode;

  for (i = 0; i < n1; i++) {
    Acode = - cos(PIh*(2*i+1)/n1);
    coordA[i] = (Acode+1)/2;
  }

  for (i = 0; i < n2; i++) {
    coordB[i] = - cos(PIh*(2*i+1)/n2);
  }

  for (i = 0; i < n3+2*pad; i++) {
    coordphi[i] = 2*PI*(i-pad)/n3; 
    if (0) printf("coordphi[%2d] = %f  %f\n",
		  i, coordphi[i], coordphi[i]*180/PI); 
  }
}




/* interpolate to point given in Cartesian coordinates
   calls function in utility/interpolation/barycentric.c
*/
double ps_u_at_xyz(int order, double x, double y, double z, double *v, int varinter)
{
  double A, B, phi, u, U;
  double coordA[n1];
  double coordB[n2];
  double coordphi[n3+2*pad];

  set_ABp(coordA,coordB,coordphi);

  xyz_to_ABp(x, y, z, &A, &B, &phi);

if (varinter == 1){

 double rho,rr,xx,AA,AAAA,fA,dfA,cor,res,xyz,xtest,ytest,ztest,bb,dxx,br;
 int i, vb=0;
 double tol=5.e-15;

    
     rr = y*y+z*z;
     rho= sqrt(rr);
     xx = x*x;
     bb = ps_b*ps_b;
    dxx = ps_dx*ps_dx;


   /* dtim: NR-method for computing the inverse transformation, see CoordABTinverse in papers/mpdata/math*/
    for(i=0; i < 15; i++) {

         AA =  A*A;
       AAAA = AA*AA;


       fA = rr + AAAA*AAAA*rr - 4.*AA*(1. + AAAA)*(bb - dxx + 2.*ps_dx*x - xx) - 
            2.*AAAA*(4.*bb + 4.*dxx + rr - 8.*ps_dx*x + 4.*xx);
   
      dfA = 8.*A*(-bb + dxx + AA*AAAA*rr - 2.*ps_dx*x + xx + 3.*AAAA*
            (-bb + dxx - 2.*ps_dx*x + xx) - AA*(4.*bb + 4.*dxx + rr - 8.*ps_dx*x + 4.*xx));
    
     cor = fA/dfA;
     res = sqrt(fA*fA);
      br = sqrt(cor*cor);

       A = A - cor;

    if (vb) printf("%le %d \n",res,i);
        if(br<tol) break;

   }
   
     if (((dless(x, ps_dx)) && (x > 0)) || ((dless(x, ps_dx)) && (x < 0)))
         B = +sqrt((2.*A*ps_b + (-1. + AA)*rho)/(2.*A*ps_b + rho - AA*rho)); 
     else if (((dless(ps_dx, x)) && (x > 0)) || ((dless(ps_dx, x)) && (x < 0)))
         B = -sqrt((2.*A*ps_b + (-1. + AA)*rho)/(2.*A*ps_b + rho - AA*rho)); 
     else if (dequal(x, ps_dx))
         B = 0.0; 
     else
         errorexit("Problem in computing B occured!");
                         

        /*additional check */ 
           puncture_ps_ABp_to_xyz(A, B, phi, &xtest, &ytest, &ztest);
           xyz = sqrt((xtest-x)*(xtest-x)+(ytest-y)*(ytest-y)+(ztest-z)*(ztest-z));
           if (vb)  printf("%le\n",xyz);
           if (xyz > 1.e-9) {printf("Inaccurate transformation on bam grid for:\n");
                              printf("x%le y%le z%le   error: %le  \n",x,y,z,xyz);}

    if (xyz > 1.e-5) errorexit("Bad transformation from spectral to cartesian!");
  

 }


  U = interpolate_tri_bar(order, A, B, phi, n1, n2, n3+2*pad,
			  coordA, coordB, coordphi, v);
  u = 2*(A-1) * U;

  
  return u;
}




/* read black hole parameters from file
   called at PRE_GRID time
   the way the AMR grids are created, we better get these parameters
   before the actual data is read
*/
int read_ps_parameters(tL *level)
{
  char *filename = Gets("punctures_ps_file");
  char s[1000], *t, *u, str[1000];
  FILE *fp;

  /* we need a file name */
  if (strlen(filename) == 0) return 0;

  /* open file */
  fp = fopen(filename, "r");
  if (!fp) 
    errorexits("could not open %s", filename);
  printf("  Reading spectral puncture parameters from %s\n", filename);

  /* read single lines */
  while (fgets(s, 1000, fp)) {

    /* stop at data section */
    t = strstr(s, "data ");
    if (t == s) break;

    /* skip comments */
    if (s[0] == '#') continue;
    
    /* everything else defines a standard parameter
       parse format: ^\S+\s+=\s+\S+
       FIXME: do some error checking
    */
    for (t = s; !isspace(*t); t++);             // find end of first key
    *t = 0;                                     // terminate string
    for (t = t+1; isspace(*t) || *t=='='; t++); // skip white space and = sign
    for (u = t; !isspace(*u); u++);             // find end of value
    *u = 0;                                     // terminate string
    if (0) printf("  %s = %s\n", s, t);
    
    /* set parameter (for old data-> convert parameters) */
    if (strcmp("bhmass1",s)==0) {AddPar("mass1",t,""); Sets("mass1", t);}
    else if (strcmp("bhmass2",s)==0) {AddPar("mass2",t,""); Sets("mass2", t);}
    else if (strcmp("bhx1",s)==0) {AddPar("px1",t,""); Sets("px1", t);}
    else if (strcmp("bhy1",s)==0) {AddPar("py1",t,""); Sets("py1", t);}
    else if (strcmp("bhz1",s)==0) {AddPar("pz1",t,""); Sets("pz1", t);}
    else if (strcmp("bhx2",s)==0) {AddPar("px2",t,""); Sets("px2", t);}
    else if (strcmp("bhy2",s)==0) {AddPar("py2",t,""); Sets("py2", t);}
    else if (strcmp("bhz2",s)==0) {AddPar("pz2",t,""); Sets("pz2", t);}
    else if (strcmp("bhpx1",s)==0) {AddPar("mx1",t,""); Sets("mx1", t);}
    else if (strcmp("bhpy1",s)==0) {AddPar("my1",t,""); Sets("my1", t);}
    else if (strcmp("bhpz1",s)==0) {AddPar("mz1",t,""); Sets("mz1", t);}
    else if (strcmp("bhpx2",s)==0) {AddPar("mx2",t,""); Sets("mx2", t);}
    else if (strcmp("bhpy2",s)==0) {AddPar("my2",t,""); Sets("my2", t);}
    else if (strcmp("bhpz2",s)==0) {AddPar("mz2",t,""); Sets("mz2", t);}
    else if (strcmp("bhsx1",s)==0) {AddPar("sx1",t,""); Sets("sx1", t);}
    else if (strcmp("bhsy1",s)==0) {AddPar("sy1",t,""); Sets("sy1", t);}
    else if (strcmp("bhsz1",s)==0) {AddPar("sz1",t,""); Sets("sz1", t);}
    else if (strcmp("bhsx2",s)==0) {AddPar("sx2",t,""); Sets("sx2", t);}
    else if (strcmp("bhsy2",s)==0) {AddPar("sy2",t,""); Sets("sy2", t);}
    else if (strcmp("bhsz2",s)==0) {AddPar("sz2",t,""); Sets("sz2", t);}
    else
      Sets(s, t);
  }
  
  Seti("nobjects",2);

  /* determine additional transformation of spectral coordinates */
  if (1) {
    double x1 = Getd("px1");
    double y1 = Getd("py1");
    double z1 = Getd("pz1"); 
    double x2 = Getd("px2");
    double y2 = Getd("py2");
    double z2 = Getd("pz2"); 
    double dx = x1 - x2;
    double dy = y1 - y2;

    /* x-axis */
    if (dx != 0 && y1 == 0 && y2 == 0 && z1 == 0 && z2 == 0) {
      ps_b = dx/2;
      ps_dx = (x1+x2)/2;
      ps_rxx = 1;
      ps_rxy = 0;
      ps_ryx = 0;
      ps_ryy = 1;
    } 

    /* y-axis */
    else if (dy != 0 && x1 == 0 && x2 == 0 && z1 == 0 && z2 == 0) {
      ps_b = dy/2;
      ps_dx = (y1+y2)/2;
      ps_rxx =  0;
      ps_rxy = +1;
      ps_ryx = -1;
      ps_ryy =  0;
    } 

    /* else */
    else
      errorexit("puncture location not allowed");
  }

  /* done */
  fclose (fp);
  return 0;
}




/* read spectral data from file
   special: pad phi direction with ghosts for periodic interpolation
            order = 4:    (-2 -1) 0 ...  n-1 (n n+1)
*/
void read_ps_data(char *filename, double **ppu)
{
  char s[1000], *t;
  FILE *fp;
  double *v;
  int nghosts;
  int i;

  /* we need a file name */
  if (strlen(filename) == 0) return;

  /* open file */
  fp = fopen(filename, "r");
  if (!fp) 
    errorexits("could not open %s", filename);
  printf("  reading data from %s\n", filename);

  /* skip to line starting with data, extract size info */
  n1 = n2 = n3 = ntotal = -1;
  while (fgets(s, 1000, fp)) {
    t = strstr(s, "data ");
    if (t != s) continue;
    sscanf(s+5, "%d%d%d", &n1, &n2, &n3);
    ntotal = n1*n2*n3;
    printf("  found data with dimensions %d x %d x %d = %d\n", 
	   n1, n2, n3, ntotal);
    break;
  }
  if (ntotal == -1)
    errorexit("file does not contain the expected data");

  /* get storage if needed */
  nghosts = n1*n2*pad;
  if (!(*ppu))
    *ppu = dmalloc(ntotal+2*nghosts);
  v = *ppu + nghosts;

  /* read data */
  i = 0;
  while (fgets(s, 1000, fp)) {
    if (i < ntotal)
      v[i] = atof(t);
    i++;
  }
  printf("  read %d data lines\n", i);
  if (i < ntotal) errorexit("file contains too few data lines");
  if (i > ntotal) errorexit("file contains too many data lines");

  /* copy data into ghosts */
  for (i = 0; i < nghosts; i++) {
    (*ppu)[i] = v[i + ntotal - nghosts];
    (*ppu)[i + ntotal + nghosts] = v[i];
  }

  if (0)
  for (i = 0; i < ntotal + 2*nghosts; i++)
    printf("yoyo %10d  %.16e\n", i-nghosts, (*ppu)[i]);

  /* done */
  fclose (fp);
}




/* function that wraps the spectral solver
   - will read from file if available 
*/
void PuncturesPS(tL *level) 
{
  tG *g = level->grid;
  int lmin = g->lmin;
  int lmax = g->lmax;
  char *filename = Gets("punctures_ps_file");
  int order = Geti("punctures_ps_order");
  double *pu_ps = 0;
  int i, l;
  int varinter = 0;

  if (Getv("punctures_ps_interpol_check", "yes"))   varinter = 1 ;


  prdivider(0);
  pad = order/2;

  /* if a file is specified, read from file */
  if (strlen(filename)) {
    printf("Reading spectral puncture data:\n");
    read_ps_data(filename, &pu_ps);
  }

  /* if a file is not specified, compute from scratch
     unfinished
  */
  else
    errorexit("PuncturesPS: only works by reading from file for now");


  /* we now have the puncture function u on a spectral grid of size n1*n2*n3
     interpolate onto all the points of all levels on this processor
  */
  printf("Interpolating from spectral to cartesian coordinates\n");
  for (l = lmin; l <= lmax; l++) {
    tL *level = g->level[l];
    double *px = Ptr(level, "x");
    double *py = Ptr(level, "y");
    double *pz = Ptr(level, "z");
    double *pu = Ptr(level, "punctures_u");
    double A, Acode, B, phi;

    bampi_openmp_start
    
    forallpoints_openmp(level, i) {

      pu[i] = ps_u_at_xyz(order, px[i], py[i], pz[i], pu_ps, varinter);

      if (0) 
	printf("ipol %5.2f %5.2f %5.2f  u = %f\n", px[i], py[i], pz[i], pu[i]);
      
      CheckForNANandINF(4,px[i], py[i], pz[i], pu[i]);
    }
    printf("  level %2d,  points %d\n", l, i);
    
    bampi_openmp_stop
  }

  /* make sure that the Robin boundary condition is not applied */
  VarNameSetBoundaryInfo("punctures_u", 0, -10, 0);

  /* done */
  free(pu_ps);
  prdivider(0);
}

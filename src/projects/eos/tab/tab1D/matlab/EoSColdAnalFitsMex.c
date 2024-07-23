/* This is a MATLAB MEX-file */
/* Doc: see EoSColdAnalFits.m */

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "HPFits.c"
#include "ShibFits.c"


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) 
{
  
  /* input */
  char *eosname;
  char *fitname;
  double *rho;    
  
  /* tmp */
  int which_eos, size;
  
  /* output */
  double* pres_p;
  double* DpDrho_p;
  double* DpDepsl_p;
  double* cs2_p;
  double* ene_p;
  double* epsl_p;
  double* DepslDrho_p;
  double* D2epslD2rho_p;
    


  /* check for the proper number of i/o arguments */
  if  (!(nrhs == 3)) mexErrMsgTxt("Three input required.");
  else if (!( nlhs == 5 || nlhs == 8)) mexErrMsgTxt("Five (or eight) output required.");
  


  /* check input arg 0 */

  /* input must be a string */
  if ( mxIsChar(prhs[0])!=1)
    mexErrMsgTxt("Input 1 must be a string.");
  
  /* input must be a row vector */
  if (mxGetM(prhs[0])!=1)
    mexErrMsgTxt("Input 1 must be a row vector.");
  
  /* copy the string data from prhs[0] into a C string input_ buf.    */
  eosname = mxArrayToString(prhs[0]);
  if(eosname == NULL) 
    mexErrMsgTxt("Could not convert input 1 to string.");



  /* check input arg 1 */

  /* input must be a string */
  if ( mxIsChar(prhs[1]) != 1)
    mexErrMsgTxt("Input 2 must be a string.");
  
  /* input must be a row vector */
  if (mxGetM(prhs[1])!=1)
    mexErrMsgTxt("Input 2 must be a row vector.");
  
  /* copy the string data from prhs[0] into a C string input_ buf.    */
  fitname = mxArrayToString(prhs[1]);
  if(fitname == NULL) 
    mexErrMsgTxt("Could not convert input 2 to string.");
    


  /* check input arg 2 */

  /* input must be a row vector */
  if (mxGetM(prhs[2]) != 1)
    mexErrMsgTxt("Input 3 must be a row vector.");

  /* input must be a real double */
  if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
    mexErrMsgTxt("Input 3 must be a noncomplex double vector.");
  
  /* Get the number of elements of rho and epsl */
  size=mxGetNumberOfElements(prhs[2]);

  /* Create a pointer to the input vectors */
  rho = mxGetPr(prhs[2]);



  /* prepare output */

  /* Create the output vectors */
  plhs[0] = mxCreateDoubleMatrix(1,size,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,size,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,size,mxREAL);
  plhs[3] = mxCreateDoubleMatrix(1,size,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,size,mxREAL);
  if ( nlhs == 8 ) {
    plhs[5] = mxCreateDoubleMatrix(1,size,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1,size,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(1,size,mxREAL);
  }
  

  /* Assign a pointer to the output */
  pres_p     = mxGetPr(plhs[0]);
  DpDrho_p   = mxGetPr(plhs[1]);
  DpDepsl_p  = mxGetPr(plhs[2]);
  cs2_p      = mxGetPr(plhs[3]);
  ene_p      = mxGetPr(plhs[4]);
  if ( nlhs == 8 ) {
    /* assign pointers also to additional output */
    epsl_p        = mxGetPr(plhs[5]);
    DepslDrho_p   = mxGetPr(plhs[6]);
    D2epslD2rho_p = mxGetPr(plhs[7]);
  } else {
    /* alloc memory for aux vars  */
    epsl_p        = (double *) malloc(size*sizeof(double));
    DepslDrho_p   = (double *) malloc(size*sizeof(double));
    D2epslD2rho_p = (double *) malloc(size*sizeof(double));
  }
  


  /* set eos index */
  if ( strcmp(eosname,"sly")==0 )      which_eos = 0; /* this eos has both HP and Shib fits */
  else if ( strcmp(eosname,"fps")==0 ) which_eos = 1; /* this eos has both HP and Shib fits */
  else if ( strcmp(eosname,"apr")==0 ) which_eos = 2; /* this eos has only Shib fit */
  else mexErrMsgTxt("unknown EoS name.");

  /* check HP fit for APR is not asked */
  if (strcmp(eosname,"apr")==0 && strcmp(fitname,"shib")!=0 ) 
    mexErrMsgTxt("HP fit for APR EoS not available.");
  


  /* call the C routines */  
  printf(" EoSColsAnalFits   EoS: %s   Fit: %s ...",eosname,fitname);
  if ( strcmp(fitname,"hp")==0 )
    
    HP_EOSFit_Cold(which_eos,size,rho,
		   pres_p,DpDrho_p,DpDepsl_p,cs2_p,ene_p,
		   epsl_p,DepslDrho_p,D2epslD2rho_p); 
  
  else if ( strcmp(fitname,"hps")==0 )  
    
    HPS_EOSFit_Cold(which_eos,size,rho,
		    pres_p,DpDrho_p,DpDepsl_p,cs2_p,ene_p,
		    epsl_p,DepslDrho_p,D2epslD2rho_p); 
  
  else if ( strcmp(fitname,"shib")==0 ) 
    
    Shib_EOSFit_Cold(which_eos,size,rho,
		     pres_p,DpDrho_p,DpDepsl_p,cs2_p,ene_p,
		     epsl_p,DepslDrho_p,D2epslD2rho_p); 
  
  else mexErrMsgTxt("unknown fit name.");
  printf(" done.\n");
  
}




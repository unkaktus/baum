/* bam_AwA.c */

#include "bam.h"
#include "AwA.h"

void bam_AwA(void)
{ 
  if (!Getv("physics", "AwA")) return;
  printf("Adding AwA\n");

  AddPar("AwA", "none", "Choose test: {none, robust_stability, robust_stability_old, linear_wave1, linear_wave2, gauge_wave1}");

  AddPar("AwA_amplitude", "0.1", "");
  AddPar("AwA_sigma", "1.0", "");  
  AddPar("AwA_compute_Dplus", "no", "Compute Dplus norm");

  if (Getv("AwA", "none")) return;

  if (Getv("AwA", "robust_stability")) 
    AddFun(POST_INITIALDATA, AwA_robust_stability, "Set AwA robust stability test");   
 
  if (Getv("AwA", "linear_wave1") || Getv("AwA", "linear_wave2"))
    AddFun(POST_INITIALDATA, AwA_linear_wave, "Set AwA linear wave initial data");

  if (Getv("AwA", "gauge_wave1"))
    AddFun(POST_INITIALDATA, AwA_gauge_wave1, "Set AwA 1D gauge wave initial data");

  if (Getv("AwA", "simple_gauge_wave"))
    AddFun(POST_INITIALDATA, AwA_simple_gauge_wave, "Set AwA Simple 3D gauge wave initial data");
  
  if (Getv("AwA", "robust_stability_old")) 
    AddFun(POST_POST_INITIALDATA, AwA_add_highfrequ, "Set AwA robust stability test, add high-frequency noise to conformal metric");    

  if (Getv("AwA_compute_Dplus", "yes")) {
    /* Note: the implementation of the Dplus norm        
       - assumes periodic BCs and 2 ghosts
       - is valid only for the y=z=0 line
    */
    AddVar("Dplus", "", "");
    AddVar("Dplusx", "", "");
    AddFun(ANALYZE, AwA_compute_Dplus, "Add Dplus norm computation during AwA test");
  }
  
}

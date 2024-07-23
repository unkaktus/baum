/* bam_Gauge.c */
/* Bernd Bruegmann 7/03, Wolfgang Tichy 12/2003 */

#include "bam.h"
#include "Gauge.h"


void (*generic_Gauge)(generic_Gauge_VARS);

void bam_Gauge() 
{
  if (!Getv("physics", "Gauge")) return;
  printf("Adding Gauge\n");
  
  AddPar("Gauge", "donothing", "options handled by Gauge");

  /* maximal slicing */
  if (Getv("Gauge", "maximal")) {
 
    /* eventually, should be registered inside evolve for iterative schemes */
    //AddFun(POST_EVOLVE, maximal, "compute lapse for maximal slicing");

    /* storage */
    AddVar("maximal_f", "", "rhs of maximal slicing equation");
    AddVar("maximal_r", "", "residuum of maximal slicing equation");

    /* precompute coefficients */
    AddVar("maximal_gi", "ij+ji", "inverse metric for maximal slicing");
    AddVar("maximal_G",  "I",     "g^ij Gamma^k_ij for maximal slicing");
    AddVar("maximal_KK", "",      "K_ij K^ij for maximal slicing");

    /* parameters */
    AddPar("maximal_verbose", "no", "whether to talk about it");
    AddPar("maximal_boundary", "robin", "boundary condition [robin]");
    AddPar("maximal_solver", "bicgstab",
	   "which solver to use [bicgstab,multigrid]");
    AddPar("maximal_tolerance", "1e-6", "tolerance for maximal solve");
    AddPar("maximal_itmax", "1000", "maximal number of iterations");
    AddPar("maximal_every", "1", 
	   "compute lapse only every other evolution step");
    AddPar("maximal_use1dstencil", "no", "use 1d stencil");
  }

  /* shift for corotating coordinates */
  if (Getv("Gauge", "corotation")) {
    //Gauge_add_betarot = add_betarot;
    //AddFun(POST_INITIALDATA, initialize_betarot,  "initialize global rotation shift vector");

    AddVar("betarot", "I", "global rotation shift vector");

    AddPar("corotation_adjust_betarot", "yes",
           "whether betarot is also adjusted if  adjust_beta  is "
           "called to adjust beta [no,yes]");
    AddPar("corotation_omega", "0", "angular velocity for global rotation");
    AddPar("corotation_omega_c", "0", "omega attenuation parameter");
    AddPar("corotation_omega_p", "0", "omega attenuation parameter");
    AddPar("corotation_scale", "1.0", "omega = scale*omega");
    AddPar("corotation_puncture_attenuation", "with_psi", 
	   "how shift at puncture is attenuated [with_psi,linear]");
    AddPar("corotation_psipower", "3", "shift /= psi^psipower");

    AddPar("corotation_attenuate", "no",
	   "pick attenuation [no,rational,exponential]");
    AddPar("corotation_attenuate_a", "1", "attenuation parameter a");
    AddPar("corotation_attenuate_b", "1", "attenuation parameter b");
    AddPar("corotation_attenuate_c", "2", "attenuation parameter c");

    AddPar("corotation_rdot", "0", "radial velocity for global motion");
    AddPar("corotation_rdot_c", "0", "radial velocity, attenuation parameter");
    AddPar("corotation_rdot_p", "0", "radial velocity, attenuation parameter");
    AddPar("corotation_exponFarAtt", "no",
           "whether attenuation approaches the farlimit exponentially");

    AddPar("corotation_monitor", "yes", "monitor asymmetry at bh");

    if (Getv("corotation_monitor", "yes"))
    {
      //AddFun(POST_OUTPUT, monitor_bh_asymmetry, "monitor asymmetry at bh");
      //AddFun(POST_OUTPUT, monitor_bh_lapse, "monitor lapse at bh");
    }

    if (Getv("corotation_monitor", "manual")) {
      //AddFun(POST_OUTPUT, commotion_manual, "counter asymmetry at bh by manual adjustment");
      AddPar("commotion_manual_t", "", "times for manual adjustment");
      AddPar("commotion_manual_s", "", "omega scale for manual adjustment");
      AddPar("commotion_manual_r", "", "rdot for manual adjustment");
    }

    if (Getv("corotation_monitor", "adjust")) {
      //AddFun(POST_OUTPUT, commotion_instantaneous, "counter asymmetry at bh by guessing an adjustment");

      AddPar("corotation_adjust_t0", "0.0", "when to start");
      AddPar("corotation_adjust_t1", "-1.0", "when to stop");
      AddPar("corotation_adjust_dt", "1.0", "how often to do it");

      AddPar("corotation_adjust_scale_max", "10000", 
	     "keep scale below this maximum");
      AddPar("corotation_adjust_scale_y_min", "-10000", 
	     "do not adjust scale if y falls below y_min");

      AddPar("corotation_shift_scale", "0", 
             "by how much we scale radial shift, to keep lapse at bh const.");
      AddPar("corotation_bh_lapse", "0.3", "desired lapse at bh");
      AddPar("corotation_shift_additionFactor", "0",
             "how much we add to the radial shift, to keep lapse at bh const.");
      AddPar("corotation_shift_radialAdjustment", "0",
             "by how much we adjust the radial shift, "
             "to keep lapse at r_1/2 = corotation_lapse_radius*m_1/2 const.");
      AddPar("corotation_lapse_radius", "0.5",
             "the radius (in units of m1/2) where we measure the lapse "
             "used to steer the radial shift (with "
             "corotation_shift_radialAdjustment)");

      AddPar("commotion_damping_x",  "0", "damping coefficient");
      AddPar("commotion_damping_y",  "0", "damping coefficient");
      AddPar("commotion_driving_x",  "0", "damping coefficient");
      AddPar("commotion_driving_y",  "0", "damping coefficient");
      AddPar("commotion_offset_x",   "0", "offset in x");
      AddPar("commotion_offset_y",   "0", "offset in y");
      AddPar("commotion_DriftVelocity", "FromPreviousCorrections", 
             "how DriftVelocity is calculated "
             "[FromPreviousCorrections,Instantaneous]");
    }

    if (Getv("corotation_monitor", "reiterate")) {
      //AddFun(POST_OUTPUT, commotion_reiterate, "counter asymmetry at bh by restarting and evolving again");
      AddPar("corotation_reiterate_t0", "0.0", "when to start");
      AddPar("corotation_reiterate_dt", "1.0", "how often to do it");
    }

    AddPar("commotion_x", "0", "measure of asymmetry");
    AddPar("commotion_y", "0", "measure of asymmetry");
  }

  /* moving puncture */
  if (Getv("Gauge", "moving_puncture")) {
    
    AddFun(ANALYZE, track_puncture,  "track puncture");
    AddFun(POST_ANALYZE, compute_moving_punc, "compute distance of punctures");
    
   /* lapse 
    AddFun(ANALYZE, write_lapse_between_punctures,
	   "write lapse between punctures");
    AddPar("alpha_center_merger", "0.3", "rule of thumb for merger");
    AddPar("alpha_center_value",  "1.0", "last value obtained");
    AddPar("alpha_center_time",   "0.0", "time of last measurement");
    */
    
    /* spin 
    AddPar("moving_puncture_spin", "no", "whether to compute puncture spin");
    if (Getv("moving_puncture_spin", "yes")) {
      AddFun(ANALYZE, moving_puncture_spin, "compute spin");
      AddVar("mprpsi2", "", "r psi^2"); 
      AddVar("mpsdh",   "", "sqrt(deth)"); 
      AddVar("mps", "I", "spin integrand"); 
      AddVar("mpt", "I", "spin integrand");
      AddVar("mpu", "I", "spin integrand");
      AddVar("mpv", "I", "spin integrand");
    }
    */
    
    AddPar("compute_moving_puncture_distance","no","use minimize/line inside parfile");
    AddPar("compute_moving_puncture_distance_method","no","AHsimpleNSBH or no");
    
    AddPar("moving_puncture_track",          "yes",          "yes/no");
    AddPar("moving_puncture_track_var",      "alpha" ,       "name of the variable");
    AddPar("moving_puncture_track_mode",     "min",          "max/min");
    AddPar("moving_puncture_track_minmove",  "0.01" ,        "value*dx");
    AddPar("moving_puncture_mir_orgin", "-1" ,   
                          "distance, when the moving puncture tracker results are mirrored");

    AddPar("moving_puncture_finboxfix",       "-1" ,   
                          "distance, when the final box is fixed");

    AddPar("moving_puncture_finboxfixv", "7.925" ,   
                          "if we have one final box, we enlarge the box by some factor");

    /*a=7.925 refers to the position, where we put the two origins of the boxes (a,a,0) and (-a,a,a)
      Thus, we should have the same amount of memory for the new bigger box, as we had it for the inspiral*/

    AddPar("moving_puncture_fixz", "NS", "pin the NS or the BH or both to the x-y-plane");

   
   if (Getv("physics", "matter"))    
    AddPar("moving_puncture_track_method", "ext" , 
           "use track_moving_puncture_extremum for BinNS");
   else 
    AddPar("moving_puncture_track_method", "int" , 
           "use track_moving_puncture_integrate  for BinBH");
 }

  AddPar("moving_puncture_properties", "no" ,
           "use coordinate spheres to determine the properties of object");

    if (Getv("moving_puncture_properties", "yes")){

       printf("Adding puncture_properties \n");
       AddFun(ANALYZE, puncture_properties, 
              "derive Spin/Momentum/Mass of punctures/NSs with BSSN variables");

       AddVar("puncture_properties_Mint","", "Integrand of M");

       AddVar("puncture_properties_Pxint","", "Integrand of Px");
       AddVar("puncture_properties_Pyint","", "Integrand of Py");
       AddVar("puncture_properties_Pzint","", "Integrand of Pz");

       AddVar("puncture_properties_Sxint","", "Integrand of Jx");
       AddVar("puncture_properties_Syint","", "Integrand of Jy");
       AddVar("puncture_properties_Szint","", "Integrand of Jz");

       AddPar("puncture_properties_r", "8 9 10 12 15 20 30", "list of radii for ADM mass");

       AddPar("puncture_properties_punc",   "2", "number of punctures");
       AddPar("puncture_properties_npoints",  " 20", "Number of points");
       AddPar("puncture_properties_circles",   "21", "Number of radial circles");
    }


}

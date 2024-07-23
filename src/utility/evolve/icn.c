/* icn.c */
/* Bernd Bruegmann 6/02 */
/* originally admc_icn.c for Cactus, Bernd Bruegmann 9/98 */

#include "bam.h"
#include "evolve.h"


/* iterative Crank-Nicholson */
void evolve_icn(tL *level, int comp) 
{
  double dt = level->dt;
  int n;

  /* initialize: u_p = u_c */
  vlcopy(u_p, u_c);

  /* iterate CN three times
     some people do not count the first step, which is identical to FT 
     but I like to call it 3 step ICN because the rhs is evaluated 3 times
  */
  for (n = 0; n < 3; n++) {

    /* u_q = (u_c + u_p)/2 */
    vlaverage(u_q, u_c, u_p);

    /* u_c = u_p + k F(u_q) */
    evolve_rhs(u_c, u_p, dt, u_q, comp);

    /* synchronization */
    bampi_vlsynchronize(u_c);
  }
}




/* Euler method */
void evolve_euler(tL *level, int comp) 
{
  double dt = level->dt;

  /* set dt to zero to obtain rhs in variable, evolve for 1 iteration */
  if (Getv("evolve_euler_debug", "yes")) dt = 0;

  /* initialize: u_p = u_c */
  vlcopy(u_p, u_c);

  /* u_c = u_p + k F(u_p) */  
  evolve_rhs(u_c, u_p, dt, u_p, comp);

  /* synchronization */
  bampi_vlsynchronize(u_c);

}

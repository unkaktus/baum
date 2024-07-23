/* ngs.c */
/* Bernd Bruegmann 5/05 */

#include "bam.h"
#include "evolve.h"


/* Newton/non-linear Gauss-Seidel
   two level implicit scheme, e.g. for generalized harmonic gauge (ghg) 
*/
void evolve_ngs(tL *level) 
{
  double dt = level->dt;
  int nmax = Geti("evolve_ngs_iterations");
  int i, n;

  /* initialize: u_q = u_p, u_p = u_c */
  vlcopy(u_q, u_p);
  vlcopy(u_p, u_c);

  /* iterate non-linear Gauss-Seidel */
  for (n = 0; n < nmax; n++) {

    /* update:  u_c = u_c - R(u_c,u_p,u_q)/dRdu_c  */
    evolve_res(u_c, u_p, u_q, u_aux);

    /* set boundary */
    set_boundary(u_c, u_q, 2*dt, u_p);

    /* synchronization */
    bampi_vlsynchronize(u_c);
  }
}




/* normals.c */
/* who? when? */

#include "bam.h"
#include "boundary.h"




/* give unit vector normal to cube boundary 
   'o' is an offset to select the boundary as the 'o'-away from the cube limits 
   return 1 if (i,j,k) is a 'o'-away boundary, 0 otherwise 
   intended for combined use with Naway macros */
int give_cube_normal (int i, int j, int k, 
		      int imax, int jmax, int kmax, int o,
		      double *shat1, double *shat2, double *shat3)
{
  double s1,s2,s3, s,oos=1.;
  int is_o_boundary = 0;
  s1 = s2 = s3 = 0.0;
  
  if(i==o)       /* left x-boundary, i.e. shat = -1,0,0 */
    s1 = -1.0;
  if(i==imax-o)  /* right x-boundary */
    s1 = 1.0;

  if(j==o)       
    s2 = -1.0;
  if(j==jmax-o)  
    s2 = 1.0;

  if(k==o)    
    s3 = -1.0;
  if(k==kmax-o)
    s3 = 1.0;
  
  s = sqrt(s1*s1 + s2*s2 + s3*s3);
  if (s!=0.) {
    oos = 1./s;
    is_o_boundary = 1;
  }

  *shat1 = s1*oos;
  *shat2 = s2*oos;
  *shat3 = s3*oos;

  return is_o_boundary;
}

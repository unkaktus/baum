/* FlagRegrid.c */
/* Bernd Bruegmann 11/03, Wolfgang Tichy 2/2004 */
/* compare src/amr/flag.c */


#include "bam.h"
#include "punctures.h"




/* read puncture positions
   temporarily done this way because FMR regrid is called before INITIALDATA
*/
void ReadPuncturePositions(void)
{
  int i, j;
  char *coord[3] = {"x", "y", "z"};
  char t[1024];

  for (i = 0; i < N_MWBL; i++) {
    sprintf(t, "bhmass%d", i+1);
    MBL[i] = Getd(t);
    if (MBL[i] == 0.0) break;
    for (j = 0; j < 3; j++) {
      sprintf(t, "bh%s%d", coord[j], i+1);
      CBL[i][j] = Getd(t);
    }
  }
}




/* flag points for refinement centered around punctures */
int FlagRegridPunctures(tL *level) 
{
  int i, j, ncube, p;
  double *f = Ptr(level, "flagregrid");
  double x, *px = Ptr(level, "x");
  double y, *py = Ptr(level, "y");
  double z, *pz = Ptr(level, "z");
  double r, r0, rx, ry; 
  double ox, oy, oz;

  /* use the more generic routine FlagRegrid_BoxesAroundCenter? */
  if(Getv("bh_fmr", "BoxesAroundCenter"))
  {
    FlagRegrid_BoxesAroundCenter(level);
    return 0;
  }

  /* make sure we know where the punctures are */
  ReadPuncturePositions();

  /* simple nesting strategy: boxes/spheres of equal size at each level 
     for punctures, they may overlap and form larger boxes
  */
  r0 = level->dx * Geti("nx")/4;

  /* initialize flags to zero and return if nothing else to do */
  forallpoints(level, i) 
    f[i] = 0;
  if (level->l >= Geti("amr_lmax")) 
    return 0;

  /* for all punctures */
  for (p = 0; MBL[p]; p++) {
    ox = CBL[p][0];
    oy = CBL[p][1];
    oz = CBL[p][2];

    /* refinements need buffer to boundary for interpolation box */
    ncube = Geti("order_RP")/2;
    if (dless(level->bbox[3], r0 + oy + (ncube-0.5)*level->dx))
      errorexit("FlagRegridPunctures: "
		"outer boundary is too close to punctures");

    /* nested boxes */
    if (Getv("amr_fmr", "nestedboxes")) {

      /* guard against round-off ?! */
      r = r0 + level->dx/100; 
      /*      + fmod(r0-ox, level->dx/10)/100
              + fmod(r0-oy, level->dy/10)/100
	      + fmod(r0-oz, level->dz/10)/100; */
      if (0) printf("flagregrid: r0 = %f\n", r);

      rx = r; /* width of box (in x-direction) = rx+rx */
      ry = r; /* length of box (in y-direction) = ry+r */
      /* 
        cover_origin explanation:
          if we cover the orgin we get a box like this:
                  ry       r
              -----------+----
              |              |          O = origin
              |              | rx
              O          P   +          P = puncture
              |              |
              |              | rx
              -----------+----
 
          if we do not cover the orgin we get a box like this:
                       r   r
                     ----+----
                     |       |
              O      |   P   | 2r
                     |       |
                     ----+----
      */
      if(Getv("bh_fmr", "cover_origin"))
        if( fabs(oy) > ry - 2.0*level->dx )
        {
          double ry_o_r_max = Getd("bh_fmr_ry_o_r_max");
          double ry_o_rx_max= Getd("bh_fmr_ry_o_rx_max");
          
          /* increase ry to cover origin */
          ry = fabs(oy) + 2.0*level->dx;
          
          /* if box gets too long, don't cover origin and reset ry to r */
          if( ry/r > ry_o_r_max )  ry=r;
          
          /* if ry/rx gets too large increase rx */
          if(ry/rx > ry_o_rx_max)  rx = ry/ry_o_rx_max;
        }
      printf("r = %e   rx = %e   ry = %e\n", r, rx, ry);

      /* set flags */
      forallpoints(level, i)
      {
	x = px[i] - ox;
	y = py[i] - oy;
        z = pz[i] - oz;
	if(x < rx && x > -rx  &&
	   z < r  && z > -r)
	{
	  if(oy > 0)
	  {
 	    if(y < r && y > -ry) f[i] = 1;
 	  }
 	  else
 	  {
	    if(y < ry && y > -r) f[i] = 1;
	  }
	}
      }
    }
  
    /* nested spheres */
    else if (Getv("amr_fmr", "nestedspheres")) {
      forallpoints(level, i) {
	x = px[i] - ox;
	y = py[i] - oy;
        z = pz[i] - oz;
	r = sqrt(x*x + y*y + z*z);
	if (r < r0) 
	  f[i] = 1;
      }
    }

    /* unknown */
    else
      errorexits("FlagRegridPunctures: amr_fmr = %s not known",
		 Gets("amr_fmr"));
  }

  /* adjust flags at boundaries */
  flagregridboundary(level);

  if (0) prvar01(level, "flagregrid");
  return 0;
}

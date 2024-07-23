/* moving_punctures.c */
/* Bernd Bruegmann, 12/2005 */
/* Jose Gonzalez, 6/2006 */
/* Wolfgang Tichy 10/2006 */
/* mth 2012 simplified */
/* dtim 2013 add switch for trackmethods */

#include "bam.h"
#include "Gauge.h"

#define PR 0

int pintoxyplane;
int trackmethod;

/* init moving punctures */
void init_moving_puncture(int n) 
{
  char str1[100],str2[100];
  if (PR) printf("Initializing moving_punctures (Gauge):\n");
  
  sprintf(str1,"moving_puncture%d_x",n);
  sprintf(str2,"px%d",n);
  if (ExistPar(str1)) return;
  AddPar(str1,Gets(str2),"location of moving_puncture in x dir");
  
  sprintf(str1,"moving_puncture%d_y",n);
  sprintf(str2,"py%d",n);
  AddPar(str1,Gets(str2),"location of moving_puncture in y dir");
  
  sprintf(str1,"moving_puncture%d_z",n);
  sprintf(str2,"pz%d",n);
  AddPar(str1,Gets(str2),"location of moving_puncture in z dir");
}
 
/* write result */
void write_tracked_puncture_N(tL *level, char *name, int N, double *result, char *str)
{
  FILE *fp;
  char *outdir = Gets("outdir");
  int n = strlen(outdir) + strlen(name) + 200;
  char *filename;
  int i;

  if (bampi_rank() != 0) return;

  filename = cmalloc(n);
  snprintf(filename, n, "%s/%s.lxyz%d", outdir, name, level->l);
  fp = fopen(filename, "a");
  if (!fp) errorexits("failed opening %s", filename);
  
  if (level->iteration == 0 && strlen(str)>1)
    fprintf(fp, "%s\n", str);
  
  for (i=0; i<N; i++)
    fprintf(fp, "%14.6e", result[i]);
  
  fprintf(fp, "%14.6e\n", level->time);

  fclose(fp);
  free(filename);
}

void write_tracked_puncture_3(tL *level, char *name,
                              double a, double b, double c, double d)
{
  double result[3];
  result[0] = d;
  result[1] = a;
  result[2] = b;
  write_tracked_puncture_N(level, name, 3, result, "");
}

void write_tracked_puncture_6(tL *level, char *name,
                              double a, double b, double c, double d,
                              double e, double f)
{
  double result[6];
  result[0] = a;
  result[1] = b;
  result[2] = c;
  result[3] = d;
  result[4] = e;
  result[5] = f;
  write_tracked_puncture_N(level, name, 6, result, "");
}

/* track moving punctures
   integrate shift at puncture location
   as a check, look for maximum of bssn_phi and minimum of alpha
   there also is the possible to use 0doutputmaxpos = bssn_phi, but 
   that just finds the closest grid point

   to do: handle puncture outside finest level

   dtim: add switch moving_puncture_track_method 
         to decide which puncture tracker we want, it is decided in firstcall to 
         save computation time, despite trackmethod 3, this is set afterwards if they
         are close 
         trackmethod 0 uses the search for the extremum (NS)
         trackmethod 1 uses the integration of the shift (BH)
         trackmethod 2 uses the extremum for object 0 and integration for the rest
         trackmethod 3 (activated if ojbects are close) uses extremum for object 0 and for 
                       the rest we set the coordinates to minus the one of object 0. This is useful 
                       for the HMNS phase and seems to be better than simple trackmethod 0
*/
static int firstcall = 1;
static char str1[1000], str2[1000];

int track_puncture(tL *level)
{ 
  char str[100];
  int i,np;
  double pold[3],pnew[3], moved;  

  if (firstcall) {
    firstcall = 0;
    sprintf(str1,"\"%14s%14s%14s%14s\"",
            "px      ","py      ","pz      ","time     ");
    sprintf(str2,"\"%14s%14s%14s%14s%14s%14s%14s%14s%14s\"",
            "px1     ","py1     ","pz1     ",
            "px2     ","py2     ","pz2     ",
            "d-proper    ","d-coord   ","time     ");

     if (Getv("moving_puncture_track_method","ext"))         trackmethod = 0; 
     else if (Getv("moving_puncture_track_method","int"))    trackmethod = 1;
     else if (Getv("moving_puncture_track_method","extint")) trackmethod = 2;
     else     errorexit("This moving puncture tracker is not implemented!");
     
   //decide if we want to pin the object to the x-y-plane (default is for NS yes)
     if      (Getv("moving_puncture_fixz","NS"))   pintoxyplane = 1;
     else if (Getv("moving_puncture_fixz","BH"))   pintoxyplane = 2;
     else if (Getv("moving_puncture_fixz","BHNS")) pintoxyplane = 3;
     else if (Getv("moving_puncture_fixz","none")) pintoxyplane = 0;
     else errorexit("Do you want to fix z=0 for the compact object?");
  }
  
  /* compute new punctures */
  for (np=1; np<=Geti("nobjects"); np++) {
    
    /* only when we are in the finest level */
     if (level->grid->lmaxpunc[np-1] != level->l) continue;
   


    /* initialize if necessary */
    if (level->iteration == 0) {
      init_moving_puncture(np);
    }
    
    /* get old location */
    for (i=0; i<3; i++) {
      pold[i] = level->grid->puncpos[np-1][i];
    }

    /* find new location by integration or whatever is implemented */
  
    if (trackmethod == 0 ) 
     moved = track_moving_puncture_extremum(level, np-1, pold, pnew);  
    else if (trackmethod == 1 ) 
     moved = track_moving_puncture_integrate(level, np-1, pold, pnew);
    else if (trackmethod == 2 ) { 
    if (np == 1)
     moved = track_moving_puncture_extremum(level, np-1, pold, pnew);
    else
     moved = track_moving_puncture_integrate(level, np-1, pold, pnew);
    } 
    else if (trackmethod == 3 ) { 
      if (np == 1)
      moved = track_moving_puncture_extremum(level, np-1, pold, pnew);
     if (np == 2) {
      for (i=0; i<3; i++) {
      pnew[i] = - (level->grid->puncpos[0][i]);
      moved = 1;  // unsure whether we need this
       }
     }
    }
    else if (trackmethod == 4 ) {
      pnew[0] = (2.*np-3.)*Getd("moving_puncture_finboxfixv");
      pnew[1] = (2.*np-3.)*Getd("moving_puncture_finboxfixv");
      pnew[2] = 0.;
    }

    /* test if new value is inside box, if not, keep the old one */
    if (!xyzinsidelevel(level,pnew[0],pnew[1],pnew[2])  && !Getv("grid", "quadrant"))
      continue;
    
    /* set new location, also needed for checkpointing */
    for (i=0; i<3; i++) {
      sprintf(str,"moving_puncture%d_%c",np,'x'+i);
      Setd(str, pnew[i]);
      level->grid->puncpos[np-1][i] = pnew[i];
    }

    /* write result */
    sprintf(str,"moving_puncture_integrate%d",np);
    write_tracked_puncture_N(level,str, 3, pnew, str1);
    
  }
 return 1;
}

int compute_moving_punc(tL *level){ 
  /* compue distance ... not finished for n punctures*/
  double d, buffer[8];
  int np1,np2,i;
  char str[100];
  if (Getv("physics","matter") && Geti("nobjects")==2 && 
      !Getv("compute_moving_puncture_distance","no")) {
    np1 = 0;
    np2 = 1;
    
    if (level->l == level->grid->lmaxpunc[np1] &&
        level->l == level->grid->lmaxpunc[np1] ) {
      
      for (i=0; i<3; i++) {
        sprintf(str,"moving_puncture%d_%c",np1+1,(i==0)?'x':(i==1)?'y':'z');
        buffer[i] = Getd(str);
      }
      for (i=0; i<3; i++) {
        sprintf(str,"moving_puncture%d_%c",np2+1,(i==0)?'x':(i==1)?'y':'z');
        buffer[3+i] = Getd(str);
      }
      
      compute_moving_puncture_distance(level, &(buffer[0]),&(buffer[3]),&(buffer[6]),&(buffer[7]));
      
      sprintf(str,"moving_puncture_distance");
      write_tracked_puncture_N(level,str, 8, buffer, str2);
    }  
  }
  
  return 1;
}





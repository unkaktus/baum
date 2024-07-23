/* outputVTK.c */
/* Bernd Bruegmann 1/2006 */

#include "bam.h"
#include "output.h"




/* handle binary formats

   VTK and ParaView are by default expecting big endian binary data
   this used to be the default on "workstations", while x86 "PCs"
   typically use little endian

   although the VTK readers can choose, I have not seen this as an option
   for the ParaView binaries, so let's write big endian

   there exists /usr/include/endian.h, but how standard is this?
   also see /usr/include/byteswap.h, which does not do 8 bytes for non-gcc?
   see include/bits/byteswap.h for gcc/x86 optimized code

   On Mac OS endian.h is at another path (?), but the code compiles
   even with "#include <endian.h>" commented out on Mac OS and Linux as well.
   However, for good measure I leave the include intact.
*/
#ifndef __APPLE__
#include <endian.h>
#endif
#if  __BYTE_ORDER == __LITTLE_ENDIAN
#define BYTE_ORDER_LITTLE 1
#else
#define BYTE_ORDER_LITTLE 0
#endif




/* in 3d and 2d it usually comes down to this: write numbers 
   efficiency shouldn't matter much since usually the disk is still slower
   defaults are based on legacy:
   - default is double if float is not specified
   - text/binary has no default (error exit)
*/
void write_raw_vtk(FILE *fp, int n, double *buffer, int nv, int j,
		   int dbl, int flt, int text, int binary)
{
  int swap = BYTE_ORDER_LITTLE;  // if little endian, then swap since
                                 // the vtk default is big endian
  float xfloat;
  double xdouble;
  int i;

  if (!text && !binary)
    errorexit("write_raw(): pick text or binary format");
  if (text && binary)
    errorexit("write_raw(): pick either text or binary format");
  if (dbl && flt)
    errorexit("write_raw(): pick either double or float format");

  if (text) {
    if (dbl) 
      for (i = 0; i < n; i++)
	fprintf(fp, "%22.15e\n", buffer[nv*i+j]);
    else
      for (i = 0; i < n; i++)
	fprintf(fp, "%16.9e\n", (float) buffer[nv*i+j]); // restricts exponent
  }

  if (binary) {
    if (sizeof(char) != 1)
      errorexit("write_raw(): size of char is not 1");

    if (!swap) {
      if (dbl) {
	double xdouble;
	for (i = 0; i < n; i++) {
	  xdouble = buffer[nv*i+j];
	  fwrite(&xdouble, sizeof(double), 1, fp);
	}
      }
      if (flt) {
	float xfloat;
	for (i = 0; i < n; i++) {
	  xfloat = buffer[nv*i+j];
	  fwrite(&xfloat, sizeof(float), 1, fp);
	}
      }
    }

    if (swap) {
      if (dbl && sizeof(double) == 8) {
	double xdouble;
	char c[8], *x;
	for (i = 0; i < n; i++) {
	  xdouble = buffer[nv*i+j];
	  x = (char *) &xdouble;
	  c[0] = x[7];
	  c[1] = x[6];
	  c[2] = x[5];
	  c[3] = x[4];
	  c[4] = x[3];
	  c[5] = x[2];
	  c[6] = x[1];
	  c[7] = x[0];
	  fwrite(c, sizeof(char), 8, fp);
	}
      }
      else if (flt && sizeof(float) == 4) {
	float xfloat;
	char c[8], *x;
	for (i = 0; i < n; i++) {
	  xfloat = buffer[nv*i+j];
	  x = (char *) &xfloat;
	  c[0] = x[3];
	  c[1] = x[2];
	  c[2] = x[1];
	  c[3] = x[0];
	  fwrite(c, sizeof(char), 4, fp);
	}
      }
      else
	errorexit("write_raw(): size of float or double is inappropriate");
    }
  }
}

/* the same for vector data */
void write_raw_vec_vtk(FILE *fp, int n, double *buffer,
                   int dbl, int flt, int text, int binary)
{
  int swap = BYTE_ORDER_LITTLE;  // if little endian, then swap since
                                 // the vtk default is big endian
  float xfloat;
  double xdouble;
  int i,j;

  if (!text && !binary)
    errorexit("write_raw(): pick text or binary format");
  if (text && binary)
    errorexit("write_raw(): pick either text or binary format");
  if (dbl && flt)
    errorexit("write_raw(): pick either double or float format");

  if (text) {
    if (dbl) 
      for (i = 0; i < n; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e\n", buffer[i],buffer[i+n],buffer[i+2*n]);
    else
      for (i = 0; i < n; i++)
        fprintf(fp, "%16.9e %16.9e %16.9e\n", (float)buffer[i],(float)buffer[i+n],(float)buffer[i+2*n]); // restricts exponent
  }

  if (binary) {
    if (sizeof(char) != 1)
      errorexit("write_raw(): size of char is not 1");

    if (!swap) {
      if (dbl) {
        double xdouble;
        for (i = 0; i < n; i++) {
          for (j=0; j<3; j++) {
            xdouble = buffer[i+j*n];
            fwrite(&xdouble, sizeof(double), 1, fp);
          }
        }
      }
      if (flt) {
        float xfloat;
        for (i = 0; i < n; i++) {
          for (j=0; j<3; j++) {
            xfloat = buffer[i+j*n];
            fwrite(&xfloat, sizeof(float), 1, fp);
          }
        }
      }
    }

    if (swap) {
      if (dbl && sizeof(double) == 8) {
        
        double xdouble;
        char c[8], *x;
        for (i = 0; i < n; i++) {
          for (j=0; j<3; j++) {
            xdouble = buffer[i+j*n];
            x = (char *) &xdouble;
            c[0] = x[7];
            c[1] = x[6];
            c[2] = x[5];
            c[3] = x[4];
            c[4] = x[3];
            c[5] = x[2];
            c[6] = x[1];
            c[7] = x[0];
            fwrite(c, sizeof(char), 8, fp);
          }
        }
      }
      else if (flt && sizeof(float) == 4) {
        float xfloat;
        char c[8], *x;
        for (i = 0; i < n; i++) {
          for (j=0; j<3; j++) {
            xfloat = buffer[i+j*n];
            x = (char *) &xfloat;
            c[0] = x[3];
            c[1] = x[2];
            c[2] = x[1];
            c[3] = x[0];
            fwrite(c, sizeof(char), 4, fp);
          }
        }
      }
      else
        errorexit("write_raw(): size of float or double is inappropriate");
    }
  }
}




/* open file for vtk writing */
FILE *fopen_vtk(char *varname, char *dirsuffix, char *suffix, int l, int n)
{
  char filename[1000], str[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp;

  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_vtk", outdir, varname, suffix, l);
  fp = fopen(filename, "r");
  if (!fp)
    mkdir(filename, 0777);  
  else
    fclose(fp);
  
  /* open file */
  snprintf(filename, 1000, "%s/%s.%s%d_vtk/%s.%s%d%s_%04d.vtk", 
           outdir, varname, suffix, l,
	   varname, suffix, l, boxstring, n);
  fp = fopen(filename, "wb");
  if (!fp) 
    errorexits("failed opening %s", filename);
  
  /* return non-null file pointer */
  return fp;
}

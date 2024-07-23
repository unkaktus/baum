/* bamcompare.c */
/* Bernd Bruegmann 10/02 */

/* standalone program to compare 1d output of bam
   - lines containing only white space are ignored
   - lines starting with " are ignored
   - compare other lines assuming there are two columns of numbers
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define NLINE 10000
char buffer0[NLINE];
char buffer1[NLINE];
int m, n;


double abstol = 1e-12;
double reltol = 1e-12;




int readnextdata(FILE *fp, char *s, int *n)
{
  char *p;

  while (1) {

    /* get next line */
    p = fgets(s, NLINE, fp);

    /* exit if error or EOF */
    if (!p) exit(0);

    /* got a new line */
    (*n)++;

    /* read another if this line starts with " */
    if (s[0] == '"') continue;
    
    /* read another if this line contains just white space */
    for (p = s; *p; p++)
      if (!isspace(*p)) break;
    if (*p) break;
  }
}



double relativeerror(double x, double y)
{
  double a = 0.5*(fabs(x) + fabs(y));

  if (a == 0) return 0;
  return (x-y)/a;
}





int main(int argc, char **argv)
{
  char *s, *t;
  double x0, x1, y0, y1;
  double ax, ay, rx, ry;
  FILE *fs, *ft;

  if (argc != 3) {
    printf("Usage:  bamcompare filename filename\n"); 
    exit(0);
  }

  s = argv[1];
  t = argv[2];

  fs = fopen(s, "r");
  if (!fs) {
    printf("Could not open file %s\n", s);
    exit(1);
  }

  ft = fopen(t, "r");
  if (!ft) {
    printf("Could not open file %s\n", t);
    fclose(fs);
    exit(1);
  }

  s = buffer0;
  t = buffer1;

  while (1) {
    readnextdata(fs, s, &m);
    readnextdata(ft, t, &n);

    sscanf(s, "%le%le", &x0, &y0);
    sscanf(t, "%le%le", &x1, &y1);

    if (0) printf("%e %e   %e %e\n", x0, x1, y0, y1);

    if (x0 != x1 || y0 != y1) {
      
      rx = relativeerror(x0, x1);
      ry = relativeerror(y0, y1);

      if (rx > reltol || ry > reltol) {

	ax = fabs(x1-x0);
	ay = fabs(y1-y0);

	if (ax > abstol || ay > abstol) {
	  printf("%4d< %s%4d> %s", m, s, n, t);
	  printf("relative errors %e %e,  absolute errors %e %e\n", 
		 rx, ry, ax, ay);
	}
      }
    }
  }

  fclose(fs);
  fclose(ft);
}

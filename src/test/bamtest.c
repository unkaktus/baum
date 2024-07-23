/* bamtest.c */
/* Bernd Bruegmann 10/02, Wolfgang Tichy 3/2003 */

/* standalone program to test bam
   the test consists in running a given parameter file with bam and to 
   compare with previously obtained output
*/

#define STRLEN 10000

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/wait.h>


int runbam(int n, char *par, int details);
int checkoutput(char *s, char *t, double reltol, double abstol,
		int details, int n);


int main(int argc, char **argv)
{
  int i, l, n = 1;
  int ret1, ret2;
  char *s, *t, *u;
  double reltol=-1.0;
  double abstol=0.0;
  char command[STRLEN];
  int details=0;
  int fails=0, passes=0, tests=0;
  
  if (argc == 1) {
    printf(
    "Usage:  bamtest [-numberofprocs] [-d]"
    "[-rt rel.tolerence] [-at abs.tolerence] name.par [name.par ...]\n"); 
    exit(0);
  }

  /* parse options */
  i = 1;
  for(i=1; (i<argc)&&(argv[i][0] == '-'); i++)
  {
    s = argv[i];

    /* number of processors */
    if (isdigit(s[1])) 
    {
      n = atoi(s+1);
      if (n <= 0) 
      {
	printf(
          "argument %s does not specify a positive number of processors\n", s);
	exit(1);
      }
    }
    else if( (strcmp(s+1,"t")==0)||(strcmp(s+1,"rt")==0) )
    {
      if(i>=argc-1) 
      {
        printf("no rel. tolerence value specified\n");
        exit(1);
      } 
      reltol=atof(argv[i+1]);
      i++;
    } 
    else if( strcmp(s+1,"at")==0 )
    {
      if(i>=argc-1) 
      {
        printf("no abs. tolerence value specified\n");
        exit(1);
      } 
      abstol=atof(argv[i+1]);
      i++;
    } 
    else if( strcmp(s+1,"d")==0 )
    {
      details=1;
    } 
    else 
    {
      printf("unknown argument %s\n", s);
      exit(1);
    }
  }
  if (i >= argc) {
    printf("Error: i >= argc \n");
    exit(1);
  }

  /* remove old log */
  snprintf(command,STRLEN,"rm -f bamtest.log >& /dev/null");
  system(command);
  
  printf("----------------------------------------------"
         "--------------------------------\n");
  printf("bamtest: Running tests on %d processors\n",n);
  printf("----------------------------------------------"
         "--------------------------------\n");
  		     
  /* for all the remaining arguments */
  for (; i < argc; i++) {

    /* strip off extension .par if it exists */
    s = strdup(argv[i]);
    l = strlen(s);
    if (l > 4 && strcmp(s+l-4, ".par") == 0)  s[l-4] = 0;
    else continue;

    /* construct name of test directory */
    u = strdup(s);
    l = strlen(u);
    t = calloc(sizeof(char), l+8);
    sprintf(t, "%s_correct", u);

    /* run bam */
    if(details) printf("==============================================================================\n");
    if(details) printf("\n==>  %s  <==\n\n",s);
    ret1=runbam(n, s, details);
    tests++;
 
    /* copy CVS and stdout files from _correct dir to test dir just created */
    snprintf(command,STRLEN,"cp -r  %s/CVS %s >& /dev/null", t, s);
    system(command);
    snprintf(command,STRLEN,"rm -f  %s/stdout.* >& /dev/null", s);
    system(command);
    snprintf(command,STRLEN,"cp -r  %s/stdout.* %s >& /dev/null", t, s);
    system(command);
    snprintf(command,STRLEN,"rm -f  %s/system.log >& /dev/null", s);
    system(command);
    snprintf(command,STRLEN,"cp -r  %s/system.log %s >& /dev/null", t, s);
    system(command);
    
    if(ret1==0)
    {
      /* compare output with given output */
      ret2=checkoutput(t, s, reltol, abstol, details, n); 
      if( (ret2==1)||(ret2==2) )
      {
        if(details) printf("\n*** ");
        printf("FAILED %s: "
               "There are significant differences!\n", s);
        if(details) printf("\n");
        fails++;
         snprintf(command,STRLEN,"%s%s%s",
          "echo >> bamtest.log;" 
          "echo \"*** ", s, 
          " FAILED: There are significant differences! ***\">> bamtest.log;"
          "echo >> bamtest.log");
         system(command);
      }
      else if( (ret2<0)||(ret2>2) )
      {
        if(details) printf("\n*** ");
        printf("FAILED %s: Cannot open files!\n",s);
        if(details) printf("\n");
        fails++;
        snprintf(command,STRLEN,"%s%s%s",
          "echo >> bamtest.log;" 
          "echo \"*** ", s, 
          " FAILED: Cannot open files! ***\">> bamtest.log;"
          "echo >> bamtest.log");
        system(command);
      }
      else
      {
        if(details) printf("\n*** ");
        printf("Passed %s\n",s);
        if(details) printf("\n");
        passes++;
        snprintf(command,STRLEN,"%s%s%s",
          "echo >> bamtest.log;" 
          "echo \"*** ", s, 
          " Passed ***\">> bamtest.log;"
          "echo >> bamtest.log");
        system(command);
      }
    }
    else
    {
      if(details) printf("\n*** ");
      printf("FAILED %s: Could not run executable!\n",s);
      if(details) printf("\n");
      fails++;
      snprintf(command,STRLEN,"%s%s%s",
          "echo ------------------------------------------------------------"
          "---------------- >> bamtest.log;"
          "echo >> bamtest.log;" 
          "echo \"==>  ", s, "  <==\">> bamtest.log;"
          "echo >> bamtest.log");
      system(command);

      snprintf(command,STRLEN,"%s%s%s",
          "echo >> bamtest.log;" 
          "echo \"*** ", s, 
          " FAILED: Could not run executable! ***\">> bamtest.log;"
          "echo >> bamtest.log");
      system(command);
    }
    /* snprintf(command,STRLEN, "echo \"    on %d Processors\" >> bamtest.log;"
                             "echo >> bamtest.log", n);
       system(command); */ 
    
    free(t);

    snprintf(command,STRLEN,"rm -f  %s/system.log >& /dev/null", s);
    system(command);
    snprintf(command,STRLEN,"rm -f  %s/stdout.* >& /dev/null", s);
    system(command);
    snprintf(command,STRLEN,"rm -rf %s/CVS >& /dev/null", s);
    system(command);
  }
  
  printf("==============================================================================\n");
  printf("\nTest(s) done:   %d\n",tests);
  printf("Test(s) passed: %d\n",passes);
  printf("Test(s) failed: %d\n\n",fails);
  
  if (fails) {
    printf("Do you want to see bamtest.log for detailed test results? ");
    scanf("%s",command);
    if(command[0]=='y') system("less bamtest.log");
  }

  return 0;
}




/* run bam for a given parameter file 
   unclear what we do if bam crashes
*/
int runbam(int n, char *par, int details)
{
  char command[STRLEN];
  char *bam = "../../exe/bam";
  int ret;

  snprintf(command, STRLEN, 
	   "mpirun -np %d %s %s >& %s.log", n, bam, par, par);
  if(details) printf("Executing  %s\n", command);
 
  ret = system(command);

  if(details) printf("Exit value %d\n", ret);
  return ret;
}




/* run diff and floatdiff to compare two output directories
   floatdiff is a separate program which can be run without bam or bamtest
   we could also run bamcompare to compare two output directories
   bamcompare is a separate program which can be run without bam or bamtest
*/
int checkoutput(char *s, char *t, double reltol, double abstol, 
                int details, int n)
{
 char command[STRLEN];
 char tolstr[STRLEN];
 int ret;
 
 if(details) printf("-----------------------------------------------------------------------\n");
 if(details) printf("Comparing  %s  and  %s\n",s,t);
 if(reltol<0.0)
 {
   printf("Please enter relative tolerance: ");
   scanf("%s",tolstr);
   reltol=atof(tolstr);
 }
 
/*  snprintf(command, STRLEN, "diff -r %s %s >& %s_vs_%s.diff",s,t,s,t);
    if(details) printf("\nExecuting  %s\n", command);
    ret = system(command);
    if(details) printf("Exit value %d\n\n\n", ret);      */
   
 if(details) printf("-----------------------------------------------------------------------\n"); 
 snprintf(command,STRLEN,
          "echo ------------------------------------------------------------"
          "---------------- >> bamtest.log;"
          "echo >> bamtest.log;" 
          "echo \"==> %s <==\">> bamtest.log;"
          "echo "
          "\"                                   \""
          "\"                            %d Processors\""
          ">> bamtest.log", t, n);
 system(command);
 snprintf(command,STRLEN,"./floatdiff -rt %e -at %e -d 1 %s %s"
          " >> bamtest.log", reltol, abstol, s, t);
 system(command);
 snprintf(command,STRLEN,"./floatdiff -rt %e -at %e -d 0 %s %s"
          " >& /dev/null",reltol, abstol, s, t);
 /* if(details) printf("Executing  %s\n\n", command); */
 fflush(stdout);
 /* ret = WEXITSTATUS(system(command)); */
 ret = system(command);
 ret = WEXITSTATUS(ret);
 fflush(stdout);
/* if(details) printf("\nExit value %d\n\n", ret); */
/*  snprintf(command, STRLEN, "rm -rf %s_vs_%s.diff >& /dev/null",s,t);
    if(details) printf("Executing  %s\n", command); 			*/
 return ret;
}

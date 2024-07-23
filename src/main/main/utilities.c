/* utilities.c */
/* Bernd Bruegmann, 12/99, Wolfgang Tichy 7/2015 */

#include "bam.h"
#include "main.h"

#define SCPR 1
#define NOSYSTEMCALL 1
#define TESTCALL 0




/* debug */
void yo()  {fprintf(stdout, "Yo!\n");fflush(stdout);}
void yo1() {fprintf(stdout, "Yo1!\n");fflush(stdout);}
void yo2() {fprintf(stdout, "Yo2!\n");fflush(stdout);}
void yo3() {fprintf(stdout, "Yo3!\n");fflush(stdout);}
void yo4() {fprintf(stdout, "Yo4!\n");fflush(stdout);}
void yo5() {fprintf(stdout, "Yo5!\n");fflush(stdout);}
void yo6() {fprintf(stdout, "Yo6!\n");fflush(stdout);}
void yo7() {fprintf(stdout, "Yo7!\n");fflush(stdout);}
void yo8() {fprintf(stdout, "Yo8!\n");fflush(stdout);}
void yo9() {fprintf(stdout, "Yo9!\n");fflush(stdout);}
void yob() {fprintf(stdout, "Yo!\n");fflush(stdout);bampi_barrier();}




int aprintf(char *target, const char *format, ...)
{
  va_list args;
  char temp[1000];
  int result;

  va_start(args, format);
  result = vsprintf(temp, format, args);
  if (result != EOF)
    strcat(target, temp);
  va_end(args);
  return result;
}




/* pretty */
void prdivider(const int n)
{
  printf("------------------------------------------------------------------------------\n");
  fflush(stdout);
}




/* minimum and maximum, also works for integers in most places */
double min2(const double x, const double y)
{
  return (x < y) ? x : y;
}
double min3(const double x, const double y, const double z)
{
  return min2(min2(x, y), z);
}
double max2(const double x, const double y)
{
  return (x > y) ? x : y;
}
double max3(const double x, const double y, const double z)
{
  return max2(max2(x, y), z);
}




/* ugh, but how universal are those built in functions? */
int system2(const char *s1, const char *s2)
{
  return system3(s1, s2, "");
}

int system3(const char *s1, const char *s2, const char *s3)
{
  char command[10000];
  int status = 0;

  // fixme: use snprintf
  if (0) { 
    sprintf(command, "%s %s %s >& /dev/null", s1, s2, s3);
    //if (SCPR) printf("System call for proc %d:  %s\n", bampi_rank(), command);
    status = system(command);
  } else {
    sprintf(command, "%s %s %s", s1, s2, s3);
    //if (SCPR) printf("System call for proc %d:  %s\n", bampi_rank(), command);
    status = system(command);
  }

  return status;
}




/* construct an argv array from a string and return number of args */
int construct_argv(char *str, char ***argv)
{
  char *str1, *token, *saveptr;
  int count;

  *argv = NULL;
  for(count=0, str1=str; ; count++, str1=NULL)
  {
    *argv = (char **) realloc(*argv, sizeof(char *)*(count+1));
    token = strtok_r(str1, " ", &saveptr);
    //printf("token=%p:%s\n", token,token);
    (*argv)[count] = token;
    if(token == NULL) break;
  }
  //printf("saveptr=%p:%s\n", saveptr,saveptr);
  return count;
}

/* run a command, without a shell */
int system_emu(const char *command)
{
  char *com = strdup(command); /* duplicate since construct_argv modifies its args */
  int ret, status;
  pid_t cpid;
  printf("system_emu: running command:\n%s\n", command);

  /* Spawn a child to run the program. */
  cpid = fork();
  if(cpid<0) /* fork failed */
  {
    printf("*** WARNING: fork failed! ***\n");
    status = ret = -911;
  }
  else if(cpid==0) /* child process */
  {
    char **argv;
    construct_argv(com, &argv);
    ret = execv(argv[0], argv);
    printf("*** WARNING: command not found, (execv returned %d) ***\n", ret);
    exit(127); /* exit child, only if execv fails */
  }
  else /* cpid!=0; parent process */
  {
    int wret = waitpid(cpid, &ret, 0); /* wait for child to exit */
    if(wret<0)
    {
      printf("*** WARNING: waitpid failed! ***\n");
      status = ret = -42;
    }  
    else
      status = ret;
    //printf("wret=%d  ret=%d  status=%d\n", wret, ret, status);
  }
  if(status!=0) printf(" -> WARNING: Return value = %d\n", status);
  free(com);
  return status;
}

/* Lock a file from current file position to the end. The lock will be
   released when the file is closed.
   fd is a file descriptor open for writing. */
int lock_curr_til_EOF(FILE *out)
{
  int fd = fileno(out); /* get file dscriptor */
  if(fd==-1) return fd; /* return -1 on error */
  return lockf(fd, F_LOCK, 0);
}
/* Unlock a file from current file position to the end.
   fd is a file descriptor open for writing. */
int unlock_curr_til_EOF(FILE *out)
{
  int fd = fileno(out); /* get file dscriptor */
  if(fd==-1) return fd; /* return -1 on error */
  return lockf(fd, F_ULOCK, 0);
}


/* special functions which are not supported by gcc */
int clear_dir(const char *which_dir)
{
  DIR           *d;
  struct dirent *dir;
  char file[256];
  struct stat s;
  
  d = opendir(which_dir);
  
  if (d) {
    while ((dir = readdir(d)) != NULL) {
      // exclude directories
      if( strcmp( dir->d_name, "." ) == 0 || strcmp( dir->d_name, ".." ) == 0) {
        continue;
      }

      sprintf(file,"%s/%s", which_dir, dir->d_name);
      //printf("*"); //print * for every deleted file

      if (opendir(file)!=NULL)  {
        clear_dir(file);
      } else {
        if (remove(file) == -1) {
          printf("\n%s\n", file);
          perror("Remove failed");
          return 1;
        }
      }
    }

    closedir(d);

    // Deleting directory
    if (rmdir(which_dir) == -1) {
      printf("%s\n", which_dir);
      perror("Remove failed");
      return 1;
    }


  }
  
  return 0;
}




int copy_file(const char *s1, const char *s2)
{
  FILE *from, *to;
  char ch, s3[1000],s4[1000];
  
  // open source file 
  if ((from = fopen(s1, "rb"))==NULL)
    errorexit("Cannot open source file.\n");

  // open destination file 
  int i = strlen(s1)-1;
  while (s1[i] != '/' && i>0) i--;
   
  int j;
  for (j=0; j<strlen(s1)-i-1; j++) 
    s3[j] = s1[j+i+1];
  s3[j] = '\0';
  sprintf(s4,"%s/%s",s2,s3);

  // BB: the above looses the first character from the parameter file name
  //     fixed below
  sprintf(s4, "%s/%s", s2, s1+i);
  if (0) {
    printf("cc=%c s1=%s s2=%s s4=%s\n", s1[i], s1, s2, s4);
  }

  /*
  printf("%d   %s  %s  %s\n",i, s1,s2,s4);
  printf("%s\n",s4);
  system("pwd");
  */
  to = fopen(s4, "wb");
  if (!to) errorexits("Cannot open destination file.",s4);

  // copy the file 
  while (!feof(from)) {
    ch = fgetc(from);
    if(ferror(from)) 
      errorexit("Error reading source file.");
    if(!feof(from)) fputc(ch, to);
    if(ferror(to)) 
      errorexit("Error writing destination file.");
  }

  if (fclose(from)==EOF)
    errorexit("Error closing source file.");

  if (fclose(to)==EOF)
    errorexit("Error closing destination file.");
  
  return 0;
}




int copy_rec(const char *s1, const char *s2)
{
  errorexit("not implemented");
  return 1;
}




/* here are the wrappers for system call or c intern system calls */
int system_isfile(const char *s1)
{
  FILE *fp = fopen(s1,"rw");
  if (fp) {
    fclose(fp);
    return 1;
  } else {
    return 0;
  }
}
int system_isdir(const char *s1)
{
  DIR  *dp = opendir(s1);
  if (dp) {
    closedir(dp);
    return 1;
  } else {
    return 0;
  }
}




int system_mkdir(const char *s1)
{
  int status;
  char command[10000];
  sprintf(command, "mkdir %s ", s1);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m\n", bampi_rank(), command);

#if TESTCALL
  DIR  *dp1 = opendir(s1);
  if ( dp1 ) {
    errorexit("  problem before mkdir");
  } else {
    if (dp1) closedir(dp1);
  }
#endif
  
#if NOSYSTEMCALL
  status = mkdir(s1,(S_IRWXU | S_IRWXG));
#else
  status = system2("mkdir", s1);
#endif
  
#if TESTCALL
  dp1 = opendir(s1);
  if ( !dp1 ) {
    errorexit("  problem after mkdir");
  } else {
    if (dp1) closedir(dp1);
  }
#endif
  return status;
}




int system_cp(const char *s1, const char *s2)
{
  int status;
  char command[10000];
  sprintf(command, "cp %s %s", s1, s2);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m\n", bampi_rank(), command);
  
#if TESTCALL
  FILE *df1 = fopen(s1,"r+");
  DIR  *dp2 = opendir(s2);
  if ( !df1 || !dp2 ) {
    errorexit("  problem before copy");
  } else {
    if (df1) fclose(df1);
    if (dp2) closedir(dp2);
  }
#endif
  
#if NOSYSTEMCALL
  status = copy_file(s1, s2);
#else
  status = system3("cp", s1, s2);
#endif
  
#if TESTCALL
  df1 = fopen(s1,"r+");
  dp2 = opendir(s2);
  if ( !df1 || !dp2 ) {
    errorexit("  problem before copy");
  } else {
    if (df1) fclose(df1);
    if (dp2) closedir(dp2);
  }
#endif
  
  return status;
}




int system_cpr(const char *s1, const char *s2)
{
  int status;
  char command[10000];
  sprintf(command, "cp %s %s", s1, s2);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m\n", bampi_rank(), command);
  
#if TESTCALL
  DIR  *dp1 = opendir(s1);
  DIR  *dp2 = opendir(s2);
  if ( !dp1 || dp2 ) {
    errorexit("  problem before copy");
  } else {
    if (dp1) closedir(dp1);
    if (dp2) closedir(dp2);
  }
#endif
  
#if NOSYSTEMCALL
  status = copy_rec(s1, s2);
#else
  status = system3("cp -r", s1, s2);
#endif
  
#if TESTCALL
  dp2 = opendir(s2);
  if ( !dp2 ) {
    errorexit("  problem before copy");
  } else {
    if (dp2) closedir(dp2);
  }
#endif
  
  return status;
}




int system_move(const char *s1, const char *s2)
{
  int status;
  char command[10000];
  sprintf(command, "mv %s %s", s1, s2);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m\n", bampi_rank(), command);
  
#if TESTCALL
  FILE *fp1 = fopen(s1,"r+");
  DIR  *dp1 = opendir(s1);
  FILE *fp2 = fopen(s2,"r+");
  DIR  *dp2 = opendir(s2);
  if ( (!dp1 && !fp1) || fp2 || dp2 || s2[strlen(s2)-1]=='/' ) {
    printf("%p %p %p %p %c\n", fp1,dp1,fp2,dp2,s2[strlen(s2)-1]);
    errorexit("problem before move");
  } else {
    if (fp1) fclose(fp1);
    if (dp1) closedir(dp1);
    if (dp2) closedir(dp2);
  }
#endif
  
#if NOSYSTEMCALL
  status = rename(s1, s2);
#else
  status = system3("mv", s1, s2);
#endif
  
#if TESTCALL
  fp1 = fopen(s1,"r+");
  dp1 = opendir(s1);
  fp2 = fopen(s2,"r+");
  dp2 = opendir(s2);
  if ( dp1 || fp1 || (!fp2 && !dp2) ) {
    errorexit("  problem after move");
  } else {
    if (dp1) closedir(dp1);
    if (dp2) closedir(dp2);
  }
#endif
  
  return status;
}




int system_rmrf(const char *s1)
{
  int status;
  char command[10000];
  sprintf(command, "rm -rf %s", s1);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m\n", bampi_rank(), command);
  
#if TESTCALL
  DIR  *dp1 = opendir(s1);
  if ( !dp1 ) {
    errorexit("  problem before remove");
  } else {
    if (dp1) closedir(dp1);
  }
#endif
  
#if NOSYSTEMCALL
  status = clear_dir(s1);
#else
  status = system2("rm -rf", s1);
#endif
  
#if TESTCALL
  dp1 = opendir(s1);
  if ( dp1 ) {
    errorexit("  problem after remove");
  } else {
    if (dp1) closedir(dp1);
  }
#endif
  
  return status;
}




int system_chmod(const char *s1, const char *m)
{
  int status;
  char command[10000];
  sprintf(command, "chmod %s %s", s1,m);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m\n", bampi_rank(), command);

#if TESTCALL
  FILE *fp1 = fopen(s1,"r+");
  DIR  *dp1 = opendir(s1);
  if ( !dp1 && !fp1 ) {
    errorexit("  problem before chmod");
  } else {
    if (dp1) closedir(dp1);
  }
#endif
  
  if (strcmp(m,"770")==0) {
    status = chmod(s1, (S_IRWXU | S_IRWXG));
  } else {
    errorexit("implement me");
  }
  
  return status;
}



/* some clusters do not support unix commands; 
int system_move(char *s1, char *s2)
{
#ifdef TESTCALL
  FILE *fp1 = fopen(s1, "r+");
  DIR  *dp1 = opendir(s1);
  FILE *fp2 = fopen(s2, "r+");
  DIR  *dp2 = opendir(s2);
  if ( (!fp1 && !dp1) // file or folder does not exist
    || (fp2)          // dest. file exist -> better stop
  ) {
    printf("  mv %s %s   %p %p %p %p\n",s1,s2,fp1,dp1,fp2,dp2);
    errorexit("  problem before move");
  } else {
    if (fp1) fclose(fp1);
    if (dp1) closedir(dp1);
    if (fp2) fclose(fp2);
    if (dp2) closedir(dp2);
  }
#endif

#ifndef NOSYSTEMCALL
  return system3("mv", s1, s2);
#else
  char command[10000];
  int status = 0;

  sprintf(command, "rename(\"%s\", \"%s\")", s1, s2);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m.\n", bampi_rank(), command);
  status = rename(s1, s2);

  return status;
#endif

#ifdef TESTCALL
  fp1 = fopen(s1, "r+");
  dp1 = opendir(s1);
  fp2 = fopen(s2, "r+");
  dp2 = opendir(s2);
  if ( (fp1 || dp1) || (!fp2 && !dp2)) {
    system("pwd");
    printf("  mv %s %s\n",s1,s2);
    errorexit("  problem after move");
  } else {
    if (fp1) fclose(fp1);
    if (dp1) closedir(dp1);
    if (fp2) fclose(fp2);
    if (dp2) closedir(dp2);
  }
#endif
}

int clear_dir(char *which_dir)
{
  DIR           *d;
  struct dirent *dir;
  char file[256];
  struct stat s;
  
  d = opendir(which_dir);
  
  if (d) {
    while ((dir = readdir(d)) != NULL) {
      // exclude directories
      if( strcmp( dir->d_name, "." ) == 0 || strcmp( dir->d_name, ".." ) == 0) {
        continue;
      }

      sprintf(file,"./%s/%s", which_dir, dir->d_name);
      printf("*"); //print * for every deleted file

      if (opendir(file)!=NULL)  {
        clear_dir(file);
      } else {
        if (remove(file) == -1) {
          printf("\n%s\n", file);
          perror("Remove failed");
          return 1;
        }
      }
    }

    closedir(d);

     // Deleting directory
    if (rmdir(which_dir) == -1) {
      printf("%s\n", which_dir);
      perror("Remove failed");
      return 1;
    }


  }
  
  return 0;
}

int system_rmrf(char *s)
{
#ifndef NOSYSTEMCALL
  return system2("rm -rf", s);
#else
  char command[10000];
  int status = 0;

  sprintf(command, "clear_dir(\"%s\")", s);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m.\n", bampi_rank(), command);
  status = clear_dir(s);

  return status;
#endif

#ifdef TESTCALL
  FILE *fp1 = fopen(s, "r+");
  DIR  *dp1 = opendir(s);
  if ( (!fp1 && !dp1)) {
    printf("  rm -rf %s\n",s);
    errorexit("  problem after rm");
  } else {
    if (fp1) fclose(fp1);
    if (dp1) closedir(dp1);
  }
#endif
}

int copy_file(char *s1, const char *s2)
{
  FILE *from, *to;
  char ch, s3[1000],s4[1000];
  
  // open source file 
  if ((from = fopen(s1, "rb"))==NULL)
    errorexit("Cannot open source file.\n");

  // open destination file 
  int i = strlen(s1)-1;
  while (s1[i] != '/' && i>0) i--;
   
  int j;
  for (j=0; j<strlen(s1)-i-1; j++) 
    s3[j] = s1[j+i+1];
  s3[j] = '\0';
  sprintf(s4,"%s/%s",s2,s3);
  
  printf("%d   %s  %s  %s\n",i, s1,s2,s4);
  printf("%s\n",s4);
  system("pwd");
  
  to = fopen(s4, "wb");
  if (!to) errorexits("Cannot open destination file.",s4);

  // copy the file 
  while (!feof(from)) {
    ch = fgetc(from);
    if(ferror(from)) 
      errorexit("Error reading source file.");
    if(!feof(from)) fputc(ch, to);
    if(ferror(to)) 
      errorexit("Error writing destination file.");
  }

  if (fclose(from)==EOF)
    errorexit("Error closing source file.");

  if (fclose(to)==EOF)
    errorexit("Error closing destination file.");
}

int system_cp(char *s1, char *s2)
{
#ifdef TESTCALL
  FILE *fp1 = fopen(s1, "r+");
  DIR  *dp1 = opendir(s1);
  FILE *fp2 = fopen(s2, "r+");
  DIR  *dp2 = opendir(s2);
  if ( (!fp1) // no file to copy
    || (dp1)  // we do not copy folder
    || (fp2)  // dest file exist
    || (dp2)) {
    system("pwd");
    printf("  cp %s %s\n",s1,s2);
    errorexit("  problem before cp");
  } else {
    if (fp1) fclose(fp1);
    if (dp1) closedir(dp1);
    if (fp2) fclose(fp2);
    if (dp2) closedir(dp2);
  }
#endif

#ifndef NOSYSTEMCALL
  return system3("cp", s1, s2);
#else
  char command[10000];
  int status = 0;

  sprintf(command, "copy(\"%s\",\"%s\")", s1,s2);
  if (SCPR) printf("System call for proc %d:  %s\n", bampi_rank(), command);
  status = copy_file(s1,s2);

  return status;
#endif
}

int system_cpr(char *s1, char *s2)
{
#ifndef NOSYSTEMCALL
  return system3("cp -r", s1, s2);
#else
  errorexit("implement me");
#endif
}

int system_mkdir(char *s)
{
#ifndef NOSYSTEMCALL
  return system2("mkdir", s);
#else
  char command[10000];
  int status = 0;

  status = mkdir(s,(S_IRWXU | S_IRWXG));
  
  sprintf(command, "mkdir(\"%s\")", s);
  if (SCPR) printf("\033[31mSystem call for proc %d:  %s \e[0m.\n", bampi_rank(), command);

  return status;
#endif
}

int system_touch(char *s)
{
#ifndef NOSYSTEMCALL
  return system2("touch", s);
#else
  errorexit("implement me");
#endif
}

int system_chmod(char *s, char *m) 
{
  if (strcmp(m,"770")==0) {
    chmod(s, (S_IRWXU | S_IRWXG));
  } else {
    errorexit("implement me");
  }
}

*/







/* storage wrappers for error checking */
double *dmalloc(const int n)
{
  double *p = (double *) malloc(sizeof(double) * n);
  
  if (!p) errorexiti("out of memory (%d double)", n);
  return p;
}

double **ddmalloc(const int n1, const int n2)
{
  int j; 

  double **p = (double**)malloc(sizeof(double*)*n1);
  for(j=0;j<n1;j++)
    p[j] = (double*)malloc(sizeof(double)*n2);
  
  return p;
}

void ddfree(double **p, const int n){
  int j; 
  for (j=0;j<n;j++) free(p[j]);
  free(p); 
}

int *imalloc(const int n)
{
  int *p = (int *) malloc(sizeof(int) * n);
  
  if (!p) errorexiti("out of memory (%d int)", n);
  return p;
}

char *cmalloc(const int n)
{
  char *p = (char *) malloc(sizeof(char) * n);
  
  if (!p) errorexiti("out of memory (%d char)", n);
  return p;
}

void *pmalloc(const int n)
{
  void *p = malloc(sizeof(void *) * n);
  
  if (!p) errorexiti("out of memory (%d void *)", n);
  return p;
}

void *bcalloc(const int n)
{
  void *p = calloc(n, 1);
  
  if (!p) errorexiti("out of memory (%d byte *)", n);
  return p;
}

double *dcalloc(const int n)
{
  double *p = calloc(n, sizeof(double));
  
  if (!p) errorexiti("out of memory (%d double)", n);
  return p;
}

int *icalloc(const int n)
{
  int *p = calloc(n, sizeof(int));
  
  if (!p) errorexiti("out of memory (%d int)", n);
  return p;
}




/* track all memory operations */
/* we do it ourselves because "man getrusage" says:
       The  above struct was taken from BSD 4.3 Reno.  Not all fields are mean�
       ingful under Linux.  Right now (Linux 2.4)  only  the  fields  ru_utime,
       ru_stime, ru_minflt, ru_majflt, and ru_nswap are maintained.
   there is no standard way to find out about memory usage?!
*/
static size_t totalbytes = 0;
static int npointers = 0;
static int npointersmax = 0;
static int nn = 0;
static void **pointer = 0;
static size_t *bytes = 0;
static int *number = 0;
static int *hide = 0;
static int pr = 1;

#ifdef TRACEMEMORY
#undef malloc
#undef calloc
#undef realloc
#undef free
#undef strdup
#endif

void prmemory_on(void) 
{
  // pr = Getv("trace_memory", "yes"); 
  // printf("Enabling memory tracing: %d\n", pr);
}

void prmemory_off(void) 
{
  pr = 0;
}

void prmemory(char *s)
{
#ifdef TRACEMEMORY
  static size_t previoustotalbytes = 0;

  printf("%20s, memory: %12u, %12d\n", s, 
	 totalbytes, (int) (totalbytes - previoustotalbytes));
  previoustotalbytes = totalbytes;
#endif
}

void prmemory_adv(char *s, tG *g)
{
#ifdef TRACEMEMORY
  static size_t previoustotalbytes = 0;
  int pts = 0;
  int l,b;
  
  if (g) {
    for (l=0; l<=g->lmax; l++) {
      if (g->level[l]->box[0])
      for (b=0; b<g->level[l]->nboxes; b++)
        pts += g->level[l]->box[b]->npoints;
    }
  
    printf("%20s, memory: %12u, %12d %6d %6d %6d %12d (%6d)\n", s, 
           totalbytes, (int) (totalbytes - previoustotalbytes),
           npointers, pts, g->level[0]->nvariables, 
           totalbytes - pts*g->level[0]->nvariables, number[npointers-1]
          );
  } else {
    printf("%20s, memory: %12u, %12d %6d %6d %6d %12d (%6d)\n", s, 
           totalbytes, (int) (totalbytes - previoustotalbytes),
           npointers, 0,0,0, number[npointers-1]
          );
  }
  
  previoustotalbytes = totalbytes;
#endif
}

void prmemory_active(void)
{
#ifdef TRACEMEMORY
  int i;
  printf("active Pointers:\n");
  for (i = npointers-1; i >= 0; i--) {
    printf("%10d)   %12p  %10d   %s\n", i,pointer[i],number[i], hide[i] ? "<-- THIS ptr IS STATIC AND WANTED":"");
  }
  printf("active Pointers:  %d\n", npointers);
#endif
}

void prmemorycond(void)
{
  if (0) prmemory("");
}

void *appendpointer(void *ptr, const size_t size)
{
  if (!ptr) return ptr;
  if (0) printf("append:  %6d  %p %u\n", nn, ptr, (unsigned int) size);

  totalbytes += size; 
  if (npointers == npointersmax) {
    npointersmax = npointers + 1;
    pointer = (void **) realloc(pointer, sizeof(void *) * npointersmax);
    bytes   = (size_t *) realloc(bytes, sizeof(size_t *) * npointersmax);
    number  = (int *) realloc(number, sizeof(int *) * npointersmax);
    hide    = (int *) realloc(hide, sizeof(int *) * npointersmax);
  }
  pointer[npointers] = ptr;
  bytes[npointers]   = size;
  number[npointers]  = nn;
  hide[npointers]    = 0;
  
  /* here you can search for memory leaks ... inserte the number of the final output */
  //if (nn==3026) errorexit("");
  npointers++;
  nn++;
  
  prmemorycond();
  return ptr;
}

void removepointer(const void *ptr)
{
  int i, j;

  if (!ptr) return;

  for (i = npointers-1; i >= 0; i--)
    if (ptr == pointer[i]) break;

  if (0) printf("remove:  %6d  %p\n", number[i],ptr);

  if (i >= 0) {
    totalbytes -= bytes[i];
    for (j = i+1; j < npointers; j++) {
      pointer[j-1] = pointer[j];
      bytes[j-1]   = bytes[j];
      number[j-1]  = number[j];
      hide[j-1]    = hide[j];
    }
    npointers--;
  }
  else 
    printf("Warning: trying to remove pointer that is not in the list.\n");
  
  prmemorycond();
}

void *replacepointer(void *ptrnew, const void *ptr, const size_t size)
{
  int i;

  if (!ptr) return appendpointer(ptrnew, size);
  
  for (i = npointers-1; i >= 0; i--)
    if (ptr == pointer[i]) break;

  if (0) printf("replace:  %6d  %p %u    (->%p)\n",  number[i], ptr, (unsigned int) size,ptrnew);

  if (i >= 0) {
    totalbytes += size - bytes[i];
    bytes[i] = size;
    pointer[i] = ptrnew;
  }
  else 
    printf("Warning: trying to replace pointer that is not in the list.\n");
  
  prmemorycond();
  return ptrnew;
}

void hidepointer(void *ptr)
{
#ifdef TRACEMEMORY
  int i;
  for (i = npointers-1; i >= 0; i--)
    if (ptr == pointer[i]) break;
  hide[i] = 1;
#endif
}


void *bam_malloc(const size_t size)
{
  return appendpointer(malloc(size), size);
}

void *bam_calloc(const size_t nmemb, const size_t size)
{
  return appendpointer(calloc(nmemb, size), nmemb*size);
}

void *bam_realloc(void *ptr, const size_t size)
{
  return replacepointer(realloc(ptr, size), ptr, size);
}

void bam_free(void *ptr)
{
  removepointer(ptr);
  free(ptr);
}

char *bam_strdup(const char *cptr)
{
  return (char *)appendpointer((void *)strdup(cptr), strlen(cptr)+1);
}




/* the one function every program should have */
/* note that bam_main.h defines a macro so that the user does not have
   to specify __FILE__ and __LINE__ for location where the error occured
*/
#undef errorexit
#undef errorexits
#undef errorexiti

#define ABORT 0

#define EXIT 1
#define THE_END if(EXIT) exit(1);  else abort();

void errorexit(const char *file, const int line, const char *s)
{
#ifdef FFLUSH 
    printf("Error(stdout): %s  ", s);
    printf("(%s, line %d)\n", file, line);
#endif
  fprintf(stderr, "Error: %s  ", s);
  fprintf(stderr, "(%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  if (ABORT) bampi_abort();
#ifdef FFLUSH 
    bampi_abort();
#endif
  bampi_finalize(0, 0);
  THE_END
}

void errorexits(const char *file, const int line, const char *s, const char *t)
{
#ifdef FFLUSH 
    printf("Error(stdout): ");
    printf(s, t);
    printf("  (%s, line %d)\n", file, line);
#endif
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, t);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  if (ABORT) bampi_abort();
#ifdef FFLUSH 
    bampi_abort();
#endif
  bampi_finalize(0, 0);
  THE_END
}

void errorexiti(const char *file, const int line, const char *s, const int i)
{
#ifdef FFLUSH 
    printf("Error(stdout): ");
    printf(s, i);
    printf("  (%s, line %d)\n", file, line);
#endif
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, i);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  if (ABORT) bampi_abort();
#ifdef FFLUSH 
    bampi_abort();
#endif
  bampi_finalize(0, 0);
  THE_END
}

/* do not write functions beyond this line because the undef/define 
   method for the errorexit functions means that they should go last */

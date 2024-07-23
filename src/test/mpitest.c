

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define PI  3.1415926535897932

/* call on quadler:
echo "cd /hugedata/mth/tmp; time cvos-launcher mpirun -np 8 /home/mth/bam_11.04/exe/mpitest > xxx.log" | qsub -j oe -l nodes=1:ppn=8 -l walltime=2000:00:00 -q long_eth -N test
*/

/* main */
int main(int argc, char *argv[]) 
{
    int size, rank,sum;
    char file[1000];
    FILE *ofp;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    /* delete old directory and make new one */
    if (rank==0) {
        system("rm -rf mpitest_out");
        system("mkdir mpitest_out");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    
    
    /* create filename */
    sprintf(file,"mpitest_out/out_p%03d",rank);
    printf("write %s\n",file);
    
    
    /* write information into file */
    ofp = fopen(file, "w");
    if (ofp==NULL) {
        printf("problems with output\n");
        exit(0);
    }
    
    fprintf(ofp,"ok!   -> \nI am %d of %d\n",rank,size);
    
    MPI_Reduce(&rank,&sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    fprintf(ofp,"%d\n",sum);
    
    fclose(ofp);
    
    
    
    
    
    /* create filename */
    int i;
    int N1=1000;
    int N2=1000;
    if (argc==3) {
      N1 = atoi(argv[1]);
      N2 = atoi(argv[2]);
    }
    double data[N1];
    sprintf(file,"mpitest_out/writetest_p%03d",rank);
    printf("write %s   (%d x %d times pi in binary)\n",file,N1,N2);
    
    
    /* write information into file */
    ofp = fopen(file, "w");
    if (ofp==NULL) {
      printf("problems with output\n");
      exit(0);
    }
    
    for (i=0; i<N1; i++)
      data[i] = PI;
    
    fprintf(ofp, "started\n");
    for (i=0; i<N2; i++) {
      
      fwrite(data, sizeof(double), N1, ofp);
      fflush(ofp);
      fprintf(ofp, "\n");
    }
    fprintf(ofp, "finished\n");
    fclose(ofp);
    
    
    
    
    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}

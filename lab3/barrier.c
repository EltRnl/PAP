#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"

int main(int argc, char** argv){

    int rank, size;
    int tag = 1;
    double start, end;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    start = MPI_Wtime();
    sleep(rank);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    printf("Process %d, execution time: %f milliseconds\n",rank, (end-start)*1000);
    MPI_Finalize();
    return 0;
}

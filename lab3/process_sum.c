#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char** argv){

    int rank, size, rank_sum;
    int tag = 2;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(tag == 1){
        if(rank == 0) printf("Reduction then broadcast\n");
        MPI_Reduce(&rank, &rank_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&rank_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("I'm process %d, recieved sum: %d\n", rank, rank_sum);
    }
    else if(tag == 2){
        if(rank == 0) printf("Global reduction\n");
        MPI_Allreduce(&rank, &rank_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        printf("I'm process %d, recieved sum: %d\n", rank, rank_sum);
    }
    MPI_Finalize();
    return 0;
}

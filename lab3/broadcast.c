#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"
#define ARRAY_SIZE 10

int main(int argc, char** argv){

    int rank, size;
    int tag = 1;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(time(NULL));
    float array[ARRAY_SIZE];
	if(rank==0){
        for(int i = 0; i < ARRAY_SIZE; i++){
            array[i] = (float)rand()/(float)(RAND_MAX);
        }
	}
    MPI_Bcast(&array, ARRAY_SIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
    printf("I'm process %d, recieved array:\n[ ",rank);
    for(int i = 0; i < ARRAY_SIZE; i++){
        printf("%f ", array[i]);
    }
    printf("]\n");
    MPI_Finalize();
    return 0;
}

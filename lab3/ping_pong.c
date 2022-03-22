#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"
#define ARRAY_SIZE 10

int main(int argc, char** argv){

    int rank, size;
    int tag = 1;
    MPI_Status status;
    double start, end;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(time(NULL));
    float array[ARRAY_SIZE];
	if(rank==0){
        for(int i = 0; i < ARRAY_SIZE; i++){
            array[i] = (float)rand()/(float)(RAND_MAX);
        }
        printf("I'm process %d, sending array\n",rank);
        start = MPI_Wtime();
        MPI_Send(array, ARRAY_SIZE, MPI_FLOAT, 1, tag, MPI_COMM_WORLD);
        MPI_Recv(&array, ARRAY_SIZE, MPI_FLOAT, 1, tag, MPI_COMM_WORLD, &status);
        end = MPI_Wtime();
        printf("I'm process %d, recieved array:\n[ ",rank);
        for(int i = 0; i < ARRAY_SIZE; i++){
            printf("%f ", array[i]);
        }
        printf("]\nRound trip took %f milliseconds\n",(end-start)*1000);
	}
    else if (rank==1){
        MPI_Recv(&array, ARRAY_SIZE, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
        /*printf("I'm process %d, recieved array:\n[ ",rank);
        for(int i = 0; i < ARRAY_SIZE; i++){
            printf("%f ", array[i]);
        }
        printf("]\n");*/
        MPI_Send(array, ARRAY_SIZE, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
    }
    else {
    	printf("Hello, i'm Process %d/%d\n", rank, size);
    }
    MPI_Finalize();
    return 0;
}

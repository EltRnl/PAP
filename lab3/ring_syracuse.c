#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

int trivial_cycle(int n){
    return (n==1 || n==2 || n==4);
}


int main(int argc, char** argv){

    int rank, size;
    int tag = 1;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(time(NULL));
    int token[2];
	if(rank==0){
        token[0] = 0;
        token[1] = rand()%1000;
        printf("Initial value: %d\n",token[1]);
        MPI_Send(token, 2, MPI_INT, 1, tag, MPI_COMM_WORLD);
        while(1){
            MPI_Recv(&token, 2, MPI_INT, size-1, tag, MPI_COMM_WORLD, &status);
            //printf("Process %d recieved %d\n",rank, token[1]);
            if(trivial_cycle(token[1])){
                break;
            }
            token[0]++;
            if(token[1]%2 == 0) token[1] = token[1]/2;
            else token[1] = 3*token[1] + 1;
            MPI_Send(token, 2, MPI_INT, 1, tag, MPI_COMM_WORLD);
        }
        printf("trivial cycle reached in %d iterations\n", token[0]);
        MPI_Send(token, 2, MPI_INT, 1, tag, MPI_COMM_WORLD);
    }
    else {
        while(1){
            MPI_Recv(&token, 2, MPI_INT, rank-1, tag, MPI_COMM_WORLD, &status);
            //printf("Process %d recieved %d\n",rank, token[1]);
            if(trivial_cycle(token[1])){
                break;
            }  
            token[0]++;
            if(token[1]%2 == 0) token[1] = token[1]/2;
            else token[1] = 3*token[1] + 1;
            MPI_Send(token, 2, MPI_INT, (rank+1)%size, tag, MPI_COMM_WORLD);
        }
        MPI_Send(token, 2, MPI_INT, (rank+1)%size, tag, MPI_COMM_WORLD);
        //printf("Process %d reached the end\n",rank);
    }
    MPI_Finalize();
    return 0;
}


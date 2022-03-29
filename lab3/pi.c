#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"

double compute_interval(double n, double x){
    x = (x/n)*(x/n);
    return (1/n)*(4/(1+x));
}

int main(int argc, char** argv){

    int rank, size;
    int tag = 1;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nb_intervals = 0;
    double process_sum = 0;
    double pi;
    if (rank == 0){
        printf("Number of intervals: ");
        fflush(stdout);
        scanf("%d",&nb_intervals);
    }
    MPI_Bcast(&nb_intervals, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < nb_intervals; i++){
        if(rank%size==i%size){
            process_sum += compute_interval((double)nb_intervals,(double)i);
        }
    }
    MPI_Reduce(&process_sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank==0) printf("Pi approximation with %d intervals: %f\n", nb_intervals, pi);
    MPI_Finalize();
    return 0;
}